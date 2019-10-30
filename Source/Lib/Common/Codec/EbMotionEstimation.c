/*
 * Copyright(c) 2019 Intel Corporation
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */
/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

#include <stdio.h>

#include "aom_dsp_rtcd.h"
#include "EbDefinitions.h"

#include "EbPictureControlSet.h"
#include "EbSequenceControlSet.h"
#include "EbMotionEstimation.h"
#include "EbUtility.h"

#include "EbComputeSAD.h"
#include "EbReferenceObject.h"
#include "EbAvcStyleMcp.h"
#include "EbMeSadCalculation.h"

#include "EbIntraPrediction.h"
#include "EbLambdaRateTables.h"
#include "EbPictureOperators.h"
#define OIS_TH_COUNT 4

int32_t OisPointTh[3][MAX_TEMPORAL_LAYERS][OIS_TH_COUNT] = {
    {// Light OIS
     {-20, 50, 150, 200},
     {-20, 50, 150, 200},
     {-20, 50, 100, 150},
     {-20, 50, 200, 300},
     {-20, 50, 200, 300},
     {-20, 50, 200, 300}},
    {// Default OIS
     {-150, 0, 150, 200},
     {-150, 0, 150, 200},
     {-125, 0, 100, 150},
     {-50, 50, 200, 300},
     {-50, 50, 200, 300},
     {-50, 50, 200, 300}},
    {// Heavy OIS
     {-400, -300, -200, 0},
     {-400, -300, -200, 0},
     {-400, -300, -200, 0},
     {-400, -300, -200, 0},
     {-400, -300, -200, 0},
     {-400, -300, -200, 0}}};

#define AVCCODEL
/********************************************
 * Constants
 ********************************************/

#define MAX_INTRA_IN_MD 9
#define REFERENCE_PIC_LIST_0 0
#define REFERENCE_PIC_LIST_1 1

/*******************************************
 * Compute8x4SAD_Default
 *   Unoptimized 8x4 SAD
 *******************************************/
uint32_t compute8x4_sad_kernel_c(
    uint8_t *src,         // input parameter, source samples Ptr
    uint32_t src_stride,  // input parameter, source stride
    uint8_t *ref,         // input parameter, reference samples Ptr
    uint32_t ref_stride)  // input parameter, reference stride
{
    uint32_t rowNumberInBlock8x4;
    uint32_t sadBlock8x4 = 0;

    for (rowNumberInBlock8x4 = 0; rowNumberInBlock8x4 < 4;
         ++rowNumberInBlock8x4) {
        sadBlock8x4 += EB_ABS_DIFF(src[0x00], ref[0x00]);
        sadBlock8x4 += EB_ABS_DIFF(src[0x01], ref[0x01]);
        sadBlock8x4 += EB_ABS_DIFF(src[0x02], ref[0x02]);
        sadBlock8x4 += EB_ABS_DIFF(src[0x03], ref[0x03]);
        sadBlock8x4 += EB_ABS_DIFF(src[0x04], ref[0x04]);
        sadBlock8x4 += EB_ABS_DIFF(src[0x05], ref[0x05]);
        sadBlock8x4 += EB_ABS_DIFF(src[0x06], ref[0x06]);
        sadBlock8x4 += EB_ABS_DIFF(src[0x07], ref[0x07]);
        src += src_stride;
        ref += ref_stride;
    }

    return sadBlock8x4;
}
/*******************************************
 * Compute8x8SAD_Default
 *   Unoptimized 8x8 SAD
 *******************************************/
uint32_t compute8x8_sad_kernel_c(
    uint8_t *src,         // input parameter, source samples Ptr
    uint32_t src_stride,  // input parameter, source stride
    uint8_t *ref,         // input parameter, reference samples Ptr
    uint32_t ref_stride)  // input parameter, reference stride
{
    uint32_t rowNumberInBlock8x8;
    uint32_t sadBlock8x8 = 0;

    for (rowNumberInBlock8x8 = 0; rowNumberInBlock8x8 < 8;
         ++rowNumberInBlock8x8) {
        sadBlock8x8 += EB_ABS_DIFF(src[0x00], ref[0x00]);
        sadBlock8x8 += EB_ABS_DIFF(src[0x01], ref[0x01]);
        sadBlock8x8 += EB_ABS_DIFF(src[0x02], ref[0x02]);
        sadBlock8x8 += EB_ABS_DIFF(src[0x03], ref[0x03]);
        sadBlock8x8 += EB_ABS_DIFF(src[0x04], ref[0x04]);
        sadBlock8x8 += EB_ABS_DIFF(src[0x05], ref[0x05]);
        sadBlock8x8 += EB_ABS_DIFF(src[0x06], ref[0x06]);
        sadBlock8x8 += EB_ABS_DIFF(src[0x07], ref[0x07]);
        src += src_stride;
        ref += ref_stride;
    }

    return sadBlock8x8;
}

/*******************************************
Calcualte SAD for 16x16 and its 8x8 sublcoks
and check if there is improvment, if yes keep
the best SAD+MV
*******************************************/
void ext_sad_calculation_8x8_16x16_c(
    uint8_t *src, uint32_t src_stride, uint8_t *ref, uint32_t ref_stride,
    uint32_t *p_best_sad8x8, uint32_t *p_best_sad16x16, uint32_t *p_best_mv8x8,
    uint32_t *p_best_mv16x16, uint32_t mv, uint32_t *p_sad16x16,
    uint32_t *p_sad8x8, EbBool sub_sad) {
    uint32_t sad16x16;

    if (sub_sad) {
        p_sad8x8[0] = (compute8x4_sad_kernel_c(src + 0 * src_stride + 0,
                                                     2 * src_stride,
                                                     ref + 0 * ref_stride + 0,
                                                     2 * ref_stride))
                      << 1;
        p_sad8x8[1] = (compute8x4_sad_kernel_c(src + 0 * src_stride + 8,
                                                     2 * src_stride,
                                                     ref + 0 * ref_stride + 8,
                                                     2 * ref_stride))
                      << 1;
        p_sad8x8[2] = (compute8x4_sad_kernel_c(src + 8 * src_stride + 0,
                                                     2 * src_stride,
                                                     ref + 8 * ref_stride + 0,
                                                     2 * ref_stride))
                      << 1;
        p_sad8x8[3] = (compute8x4_sad_kernel_c(src + 8 * src_stride + 8,
                                                     2 * src_stride,
                                                     ref + 8 * ref_stride + 8,
                                                     2 * ref_stride))
                      << 1;
    } else {
        p_sad8x8[0] = compute8x8_sad_kernel_c(src + 0 * src_stride + 0,
                                            src_stride,
                                            ref + 0 * ref_stride + 0,
                                            ref_stride);
        p_sad8x8[1] = compute8x8_sad_kernel_c(src + 0 * src_stride + 8,
                                            src_stride,
                                            ref + 0 * ref_stride + 8,
                                            ref_stride);
        p_sad8x8[2] = compute8x8_sad_kernel_c(src + 8 * src_stride + 0,
                                            src_stride,
                                            ref + 8 * ref_stride + 0,
                                            ref_stride);
        p_sad8x8[3] = compute8x8_sad_kernel_c(src + 8 * src_stride + 8,
                                            src_stride,
                                            ref + 8 * ref_stride + 8,
                                            ref_stride);
    }

    if (p_sad8x8[0] < p_best_sad8x8[0]) {
        p_best_sad8x8[0] = (uint32_t)p_sad8x8[0];
        p_best_mv8x8[0] = mv;
    }

    if (p_sad8x8[1] < p_best_sad8x8[1]) {
        p_best_sad8x8[1] = (uint32_t)p_sad8x8[1];
        p_best_mv8x8[1] = mv;
    }

    if (p_sad8x8[2] < p_best_sad8x8[2]) {
        p_best_sad8x8[2] = (uint32_t)p_sad8x8[2];
        p_best_mv8x8[2] = mv;
    }

    if (p_sad8x8[3] < p_best_sad8x8[3]) {
        p_best_sad8x8[3] = (uint32_t)p_sad8x8[3];
        p_best_mv8x8[3] = mv;
    }

    sad16x16 = p_sad8x8[0] + p_sad8x8[1] + p_sad8x8[2] + p_sad8x8[3];
    if (sad16x16 < p_best_sad16x16[0]) {
        p_best_sad16x16[0] = (uint32_t)sad16x16;
        p_best_mv16x16[0] = mv;
    }

    *p_sad16x16 = (uint32_t)sad16x16;
}

/*******************************************
Calcualte SAD for 32x32,64x64 from 16x16
and check if there is improvment, if yes keep
the best SAD+MV
*******************************************/
void ext_sad_calculation_32x32_64x64_c(uint32_t *p_sad16x16,
                                     uint32_t *p_best_sad32x32,
                                     uint32_t *p_best_sad64x64,
                                     uint32_t *p_best_mv32x32,
                                     uint32_t *p_best_mv64x64, uint32_t mv,
                                     uint32_t *p_sad32x32) {
    uint32_t sad32x32_0, sad32x32_1, sad32x32_2, sad32x32_3, sad64x64;

    p_sad32x32[0] = sad32x32_0 =
        p_sad16x16[0] + p_sad16x16[1] + p_sad16x16[2] + p_sad16x16[3];
    if (sad32x32_0 < p_best_sad32x32[0]) {
        p_best_sad32x32[0] = sad32x32_0;
        p_best_mv32x32[0] = mv;
    }

    p_sad32x32[1] = sad32x32_1 =
        p_sad16x16[4] + p_sad16x16[5] + p_sad16x16[6] + p_sad16x16[7];
    if (sad32x32_1 < p_best_sad32x32[1]) {
        p_best_sad32x32[1] = sad32x32_1;
        p_best_mv32x32[1] = mv;
    }

    p_sad32x32[2] = sad32x32_2 =
        p_sad16x16[8] + p_sad16x16[9] + p_sad16x16[10] + p_sad16x16[11];
    if (sad32x32_2 < p_best_sad32x32[2]) {
        p_best_sad32x32[2] = sad32x32_2;
        p_best_mv32x32[2] = mv;
    }

    p_sad32x32[3] = sad32x32_3 =
        p_sad16x16[12] + p_sad16x16[13] + p_sad16x16[14] + p_sad16x16[15];
    if (sad32x32_3 < p_best_sad32x32[3]) {
        p_best_sad32x32[3] = sad32x32_3;
        p_best_mv32x32[3] = mv;
    }
    sad64x64 = sad32x32_0 + sad32x32_1 + sad32x32_2 + sad32x32_3;
    if (sad64x64 < p_best_sad64x64[0]) {
        p_best_sad64x64[0] = sad64x64;
        p_best_mv64x64[0] = mv;
    }
}

/*******************************************
 * GetEightHorizontalSearchPointResults_8x8_16x16_PU
 *******************************************/
void get_eight_horizontal_search_point_results_8x8_16x16_pu_c(
    uint8_t *src, uint32_t src_stride, uint8_t *ref, uint32_t ref_stride,
    uint32_t *p_best_sad8x8, uint32_t *p_best_mv8x8, uint32_t *p_best_sad16x16,
    uint32_t *p_best_mv16x16, uint32_t mv, uint16_t *p_sad16x16,
    EbBool sub_sad) {
    uint32_t xSearchIndex;
    int16_t xMv, yMv;
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

    for (xSearchIndex = 0; xSearchIndex < 8; xSearchIndex++) {
        if (sub_sad) {
            sad8x8[0] = compute8x4_sad_kernel_c(
                            src + 0 * src_stride + 0,
                            2 * src_stride,
                            ref + 0 * ref_stride + 0 + xSearchIndex,
                            2 * ref_stride)
                        << 1;
            sad8x8[1] = compute8x4_sad_kernel_c(
                            src + 0 * src_stride + 8,
                            2 * src_stride,
                            ref + 0 * ref_stride + 8 + xSearchIndex,
                            2 * ref_stride)
                        << 1;
            sad8x8[2] = compute8x4_sad_kernel_c(
                            src + 8 * src_stride + 0,
                            2 * src_stride,
                            ref + 8 * ref_stride + 0 + xSearchIndex,
                            2 * ref_stride)
                        << 1;
            sad8x8[3] = compute8x4_sad_kernel_c(
                            src + 8 * src_stride + 8,
                            2 * src_stride,
                            ref + 8 * ref_stride + 8 + xSearchIndex,
                            2 * ref_stride)
                        << 1;
        } else {
            sad8x8[0] =
                compute8x8_sad_kernel_c(src + 0 * src_stride + 0,
                                      src_stride,
                                      ref + 0 * ref_stride + 0 + xSearchIndex,
                                      ref_stride);
            sad8x8[1] =
                compute8x8_sad_kernel_c(src + 0 * src_stride + 8,
                                      src_stride,
                                      ref + 0 * ref_stride + 8 + xSearchIndex,
                                      ref_stride);
            sad8x8[2] =
                compute8x8_sad_kernel_c(src + 8 * src_stride + 0,
                                      src_stride,
                                      ref + 8 * ref_stride + 0 + xSearchIndex,
                                      ref_stride);
            sad8x8[3] =
                compute8x8_sad_kernel_c(src + 8 * src_stride + 8,
                                      src_stride,
                                      ref + 8 * ref_stride + 8 + xSearchIndex,
                                      ref_stride);
        }

        // 8x8_0
        if (sad8x8[0] < p_best_sad8x8[0]) {
            p_best_sad8x8[0] = sad8x8[0];
            xMv = _MVXT(mv) + (int16_t)xSearchIndex * 4;
            yMv = _MVYT(mv);
            p_best_mv8x8[0] = ((uint16_t)yMv << 16) | ((uint16_t)xMv);
        }

        // 8x8_1
        if (sad8x8[1] < p_best_sad8x8[1]) {
            p_best_sad8x8[1] = sad8x8[1];
            xMv = _MVXT(mv) + (int16_t)xSearchIndex * 4;
            yMv = _MVYT(mv);
            p_best_mv8x8[1] = ((uint16_t)yMv << 16) | ((uint16_t)xMv);
        }

        // 8x8_2
        if (sad8x8[2] < p_best_sad8x8[2]) {
            p_best_sad8x8[2] = sad8x8[2];
            xMv = _MVXT(mv) + (int16_t)xSearchIndex * 4;
            yMv = _MVYT(mv);
            p_best_mv8x8[2] = ((uint16_t)yMv << 16) | ((uint16_t)xMv);
        }

        // 8x8_3
        if (sad8x8[3] < p_best_sad8x8[3]) {
            p_best_sad8x8[3] = sad8x8[3];
            xMv = _MVXT(mv) + (int16_t)xSearchIndex * 4;
            yMv = _MVYT(mv);
            p_best_mv8x8[3] = ((uint16_t)yMv << 16) | ((uint16_t)xMv);
        }

        // 16x16
        sad16x16 = (uint16_t)(sad8x8[0] + sad8x8[1] + sad8x8[2] + sad8x8[3]);
        p_sad16x16[xSearchIndex] =
            sad16x16;  // store the intermediate 16x16 SAD for 32x32.
        if ((uint32_t)(sad16x16) < p_best_sad16x16[0]) {
            p_best_sad16x16[0] = sad16x16;
            xMv = _MVXT(mv) + (int16_t)xSearchIndex * 4;
            yMv = _MVYT(mv);
            p_best_mv16x16[0] = ((uint16_t)yMv << 16) | ((uint16_t)xMv);
        }
    }
}

/*******************************************
Calcualte SAD for 32x32,64x64 from 16x16
and check if there is improvement, if yes keep
the best SAD+MV
*******************************************/
void get_eight_horizontal_search_point_results_32x32_64x64_pu_c(
    uint16_t *p_sad16x16, uint32_t *p_best_sad32x32, uint32_t *p_best_sad64x64,
    uint32_t *p_best_mv32x32, uint32_t *p_best_mv64x64, uint32_t mv) {
    int16_t xMv, yMv;
    uint32_t sad32x32_0, sad32x32_1, sad32x32_2, sad32x32_3, sad64x64;
    uint32_t xSearchIndex;

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

    for (xSearchIndex = 0; xSearchIndex < 8; xSearchIndex++) {
        // 32x32_0
        sad32x32_0 = p_sad16x16[0 * 8 + xSearchIndex] +
                     p_sad16x16[1 * 8 + xSearchIndex] +
                     p_sad16x16[2 * 8 + xSearchIndex] +
                     p_sad16x16[3 * 8 + xSearchIndex];

        if (sad32x32_0 < p_best_sad32x32[0]) {
            p_best_sad32x32[0] = sad32x32_0;
            xMv = _MVXT(mv) + (int16_t)xSearchIndex * 4;
            yMv = _MVYT(mv);
            p_best_mv32x32[0] = ((uint16_t)yMv << 16) | ((uint16_t)xMv);
        }

        // 32x32_1
        sad32x32_1 = p_sad16x16[4 * 8 + xSearchIndex] +
                     p_sad16x16[5 * 8 + xSearchIndex] +
                     p_sad16x16[6 * 8 + xSearchIndex] +
                     p_sad16x16[7 * 8 + xSearchIndex];

        if (sad32x32_1 < p_best_sad32x32[1]) {
            p_best_sad32x32[1] = sad32x32_1;
            xMv = _MVXT(mv) + (int16_t)xSearchIndex * 4;
            yMv = _MVYT(mv);
            p_best_mv32x32[1] = ((uint16_t)yMv << 16) | ((uint16_t)xMv);
        }

        // 32x32_2
        sad32x32_2 = p_sad16x16[8 * 8 + xSearchIndex] +
                     p_sad16x16[9 * 8 + xSearchIndex] +
                     p_sad16x16[10 * 8 + xSearchIndex] +
                     p_sad16x16[11 * 8 + xSearchIndex];

        if (sad32x32_2 < p_best_sad32x32[2]) {
            p_best_sad32x32[2] = sad32x32_2;
            xMv = _MVXT(mv) + (int16_t)xSearchIndex * 4;
            yMv = _MVYT(mv);
            p_best_mv32x32[2] = ((uint16_t)yMv << 16) | ((uint16_t)xMv);
        }

        // 32x32_3
        sad32x32_3 = p_sad16x16[12 * 8 + xSearchIndex] +
                     p_sad16x16[13 * 8 + xSearchIndex] +
                     p_sad16x16[14 * 8 + xSearchIndex] +
                     p_sad16x16[15 * 8 + xSearchIndex];

        if (sad32x32_3 < p_best_sad32x32[3]) {
            p_best_sad32x32[3] = sad32x32_3;
            xMv = _MVXT(mv) + (int16_t)xSearchIndex * 4;
            yMv = _MVYT(mv);
            p_best_mv32x32[3] = ((uint16_t)yMv << 16) | ((uint16_t)xMv);
        }

        // 64x64
        sad64x64 = sad32x32_0 + sad32x32_1 + sad32x32_2 + sad32x32_3;
        if (sad64x64 < p_best_sad64x64[0]) {
            p_best_sad64x64[0] = sad64x64;
            xMv = _MVXT(mv) + (int16_t)xSearchIndex * 4;
            yMv = _MVYT(mv);
            p_best_mv64x64[0] = ((uint16_t)yMv << 16) | ((uint16_t)xMv);
        }
    }
}

/*******************************************
Calcualte SAD for 16x16 and its 8x8 sublcoks
and check if there is improvment, if yes keep
the best SAD+MV
*******************************************/
void sad_calculation_8x8_16x16_c(uint8_t *src, uint32_t src_stride, uint8_t *ref,
                               uint32_t ref_stride, uint32_t *p_best_sad8x8,
                               uint32_t *p_best_sad16x16,
                               uint32_t *p_best_mv8x8, uint32_t *p_best_mv16x16,
                               uint32_t mv, uint32_t *p_sad16x16,
                               EbBool sub_sad) {
    uint64_t sad8x8[4];
    uint64_t sad16x16;

    if (sub_sad) {
        sad8x8[0] = (compute8x4_sad_kernel_c(src + 0 * src_stride + 0,
                                           2 * src_stride,
                                           ref + 0 * ref_stride + 0,
                                           2 * ref_stride))
                    << 1;
        sad8x8[1] = (compute8x4_sad_kernel_c(src + 0 * src_stride + 8,
                                           2 * src_stride,
                                           ref + 0 * ref_stride + 8,
                                           2 * ref_stride))
                    << 1;
        sad8x8[2] = (compute8x4_sad_kernel_c(src + 8 * src_stride + 0,
                                           2 * src_stride,
                                           ref + 8 * ref_stride + 0,
                                           2 * ref_stride))
                    << 1;
        sad8x8[3] = (compute8x4_sad_kernel_c(src + 8 * src_stride + 8,
                                           2 * src_stride,
                                           ref + 8 * ref_stride + 8,
                                           2 * ref_stride))
                    << 1;
    } else {
        sad8x8[0] = compute8x8_sad_kernel_c(src + 0 * src_stride + 0,
                                          src_stride,
                                          ref + 0 * ref_stride + 0,
                                          ref_stride);
        sad8x8[1] = compute8x8_sad_kernel_c(src + 0 * src_stride + 8,
                                          src_stride,
                                          ref + 0 * ref_stride + 8,
                                          ref_stride);
        sad8x8[2] = compute8x8_sad_kernel_c(src + 8 * src_stride + 0,
                                          src_stride,
                                          ref + 8 * ref_stride + 0,
                                          ref_stride);
        sad8x8[3] = compute8x8_sad_kernel_c(src + 8 * src_stride + 8,
                                          src_stride,
                                          ref + 8 * ref_stride + 8,
                                          ref_stride);
    }

    if (sad8x8[0] < p_best_sad8x8[0]) {
        p_best_sad8x8[0] = (uint32_t)sad8x8[0];
        p_best_mv8x8[0] = mv;
    }

    if (sad8x8[1] < p_best_sad8x8[1]) {
        p_best_sad8x8[1] = (uint32_t)sad8x8[1];
        p_best_mv8x8[1] = mv;
    }

    if (sad8x8[2] < p_best_sad8x8[2]) {
        p_best_sad8x8[2] = (uint32_t)sad8x8[2];
        p_best_mv8x8[2] = mv;
    }

    if (sad8x8[3] < p_best_sad8x8[3]) {
        p_best_sad8x8[3] = (uint32_t)sad8x8[3];
        p_best_mv8x8[3] = mv;
    }

    sad16x16 = sad8x8[0] + sad8x8[1] + sad8x8[2] + sad8x8[3];
    if (sad16x16 < p_best_sad16x16[0]) {
        p_best_sad16x16[0] = (uint32_t)sad16x16;
        p_best_mv16x16[0] = mv;
    }

    *p_sad16x16 = (uint32_t)sad16x16;
}

/*******************************************
Calcualte SAD for 32x32,64x64 from 16x16
and check if there is improvment, if yes keep
the best SAD+MV
*******************************************/
void sad_calculation_32x32_64x64_c(uint32_t *p_sad16x16,
                                 uint32_t *p_best_sad32x32,
                                 uint32_t *p_best_sad64x64,
                                 uint32_t *p_best_mv32x32,
                                 uint32_t *p_best_mv64x64, uint32_t mv) {
    uint32_t sad32x32_0, sad32x32_1, sad32x32_2, sad32x32_3, sad64x64;

    sad32x32_0 = p_sad16x16[0] + p_sad16x16[1] + p_sad16x16[2] + p_sad16x16[3];
    if (sad32x32_0 < p_best_sad32x32[0]) {
        p_best_sad32x32[0] = sad32x32_0;
        p_best_mv32x32[0] = mv;
    }

    sad32x32_1 = p_sad16x16[4] + p_sad16x16[5] + p_sad16x16[6] + p_sad16x16[7];
    if (sad32x32_1 < p_best_sad32x32[1]) {
        p_best_sad32x32[1] = sad32x32_1;
        p_best_mv32x32[1] = mv;
    }

    sad32x32_2 =
        p_sad16x16[8] + p_sad16x16[9] + p_sad16x16[10] + p_sad16x16[11];
    if (sad32x32_2 < p_best_sad32x32[2]) {
        p_best_sad32x32[2] = sad32x32_2;
        p_best_mv32x32[2] = mv;
    }

    sad32x32_3 =
        p_sad16x16[12] + p_sad16x16[13] + p_sad16x16[14] + p_sad16x16[15];
    if (sad32x32_3 < p_best_sad32x32[3]) {
        p_best_sad32x32[3] = sad32x32_3;
        p_best_mv32x32[3] = mv;
    }

    sad64x64 = sad32x32_0 + sad32x32_1 + sad32x32_2 + sad32x32_3;
    if (sad64x64 < p_best_sad64x64[0]) {
        p_best_sad64x64[0] = sad64x64;
        p_best_mv64x64[0] = mv;
    }
}

#define BLK_NUM 5
/**********************************************************
Calcualte the best SAD from Rect H, V and H4, V4 partitions

and return the best partition index
***********************************************************/
void nsq_me_analysis(uint32_t *p_sad64x32, uint32_t *p_sad32x16,
                     uint32_t *p_sad16x8, uint32_t *p_sad32x64,
                     uint32_t *p_sad16x32, uint32_t *p_sad8x16,
                     uint32_t *p_sad32x8, uint32_t *p_sad8x32,
                     uint32_t *p_sad64x16, uint32_t *p_sad16x64,
                     uint8_t *p_nsq_64x64, uint8_t *p_nsq_32x32,
                     uint8_t *p_nsq_16x16, uint8_t *p_nsq_8x8) {
    uint32_t sad[BLK_NUM];  // sad_N, sad_H, sad_V, sad_H4, sad_V4, sad_S;
    uint32_t best_nsq_sad;
    uint8_t nsq_index;
    /*64x64*/
    // sad[0] = p_sad64x64;
    sad[1] = p_sad64x32[0] + p_sad64x32[1];
    sad[2] = p_sad32x64[0] + p_sad32x64[1];
    sad[3] = p_sad64x16[0] + p_sad64x16[1] + p_sad64x16[2] + p_sad64x16[3];
    sad[4] = p_sad16x64[0] + p_sad16x64[1] + p_sad16x64[2] + p_sad16x64[3];
    // sad[5] = p_sad32x32[0] + p_sad32x32[1] + p_sad32x32[2] + p_sad32x32[3];
    best_nsq_sad = MAX_SAD_VALUE;
    for (nsq_index = 1; nsq_index < BLK_NUM; nsq_index++) {
        if (sad[nsq_index] < best_nsq_sad) {
            best_nsq_sad = sad[nsq_index];
            *p_nsq_64x64 = nsq_index;
        }
    }
    /*32x32*/
    // 32x32_0
    // sad[0] = p_sad32x32[0];
    sad[1] = p_sad32x16[0] + p_sad32x16[1];
    sad[2] = p_sad16x32[0] + p_sad16x32[1];
    sad[3] = p_sad32x8[0] + p_sad32x8[1] + p_sad32x8[2] + p_sad32x8[3];
    sad[4] = p_sad8x32[0] + p_sad8x32[1] + p_sad8x32[2] + p_sad8x32[3];
    // sad[5] = p_sad16x16[0] + p_sad16x16[1] + p_sad16x16[2] + p_sad16x16[3];
    best_nsq_sad = MAX_SAD_VALUE;
    for (nsq_index = 1; nsq_index < BLK_NUM; nsq_index++) {
        if (sad[nsq_index] < best_nsq_sad) {
            best_nsq_sad = sad[nsq_index];
            p_nsq_32x32[0] = nsq_index;
        }
    }
    // 32x32_1
    // sad[0] = p_sad32x32[1];
    sad[1] = p_sad32x16[2] + p_sad32x16[3];
    sad[2] = p_sad16x32[2] + p_sad16x32[3];
    sad[3] = p_sad32x8[4] + p_sad32x8[5] + p_sad32x8[6] + p_sad32x8[7];
    sad[4] = p_sad8x32[4] + p_sad8x32[5] + p_sad8x32[6] + p_sad8x32[7];
    // sad[5] = p_sad16x16[4] + p_sad16x16[5] + p_sad16x16[6] + p_sad16x16[7];
    best_nsq_sad = MAX_SAD_VALUE;
    for (nsq_index = 1; nsq_index < BLK_NUM; nsq_index++) {
        if (sad[nsq_index] < best_nsq_sad) {
            best_nsq_sad = sad[nsq_index];
            p_nsq_32x32[1] = nsq_index;
        }
    }
    // 32x32_2
    // sad[0] = p_sad32x32[2];
    sad[1] = p_sad32x16[4] + p_sad32x16[5];
    sad[2] = p_sad16x32[4] + p_sad16x32[5];
    sad[3] = p_sad32x8[8] + p_sad32x8[9] + p_sad32x8[10] + p_sad32x8[11];
    sad[4] = p_sad8x32[8] + p_sad8x32[9] + p_sad8x32[10] + p_sad8x32[11];
    // sad[5] = p_sad16x16[8] + p_sad16x16[9] + p_sad16x16[10] + p_sad16x16[11];
    best_nsq_sad = MAX_SAD_VALUE;
    for (nsq_index = 1; nsq_index < BLK_NUM; nsq_index++) {
        if (sad[nsq_index] < best_nsq_sad) {
            best_nsq_sad = sad[nsq_index];
            p_nsq_32x32[2] = nsq_index;
        }
    }
    // 32x32_3
    // sad[0] = p_sad32x32[3];
    sad[1] = p_sad32x16[6] + p_sad32x16[7];
    sad[2] = p_sad16x32[6] + p_sad16x32[7];
    sad[3] = p_sad32x8[12] + p_sad32x8[13] + p_sad32x8[14] + p_sad32x8[15];
    sad[4] = p_sad8x32[12] + p_sad8x32[13] + p_sad8x32[14] + p_sad8x32[15];
    // sad[5] = p_sad16x16[12] + p_sad16x16[13] + p_sad16x16[14] +
    // p_sad16x16[15];
    best_nsq_sad = MAX_SAD_VALUE;
    for (nsq_index = 1; nsq_index < BLK_NUM; nsq_index++) {
        if (sad[nsq_index] < best_nsq_sad) {
            best_nsq_sad = sad[nsq_index];
            p_nsq_32x32[3] = nsq_index;
        }
    }
    /*16x16*/
    // 16x16_0
    // sad[0] = p_sad16x16[0];
    sad[1] = p_sad16x8[0] + p_sad16x8[1];
    sad[2] = p_sad8x16[0] + p_sad8x16[1];
    sad[3] = MAX_SAD_VALUE;
    sad[4] = MAX_SAD_VALUE;
    // sad[5] = p_sad8x8[0] + p_sad8x8[1] + p_sad8x8[2] + p_sad8x8[3];
    best_nsq_sad = MAX_SAD_VALUE;
    for (nsq_index = 1; nsq_index < BLK_NUM; nsq_index++) {
        if (sad[nsq_index] < best_nsq_sad) {
            best_nsq_sad = sad[nsq_index];
            p_nsq_16x16[0] = nsq_index;
        }
    }
    p_nsq_8x8[0] = p_nsq_8x8[1] = p_nsq_8x8[2] = p_nsq_8x8[3] = p_nsq_16x16[0];
    // 16x16_1
    // sad[0] = p_sad16x16[1];
    sad[1] = p_sad16x8[2] + p_sad16x8[3];
    sad[2] = p_sad8x16[2] + p_sad8x16[3];
    sad[3] = MAX_SAD_VALUE;
    sad[4] = MAX_SAD_VALUE;
    // sad[5] = p_sad8x8[4] + p_sad8x8[5] + p_sad8x8[6] + p_sad8x8[7];
    best_nsq_sad = MAX_SAD_VALUE;
    for (nsq_index = 1; nsq_index < BLK_NUM; nsq_index++) {
        if (sad[nsq_index] < best_nsq_sad) {
            best_nsq_sad = sad[nsq_index];
            p_nsq_16x16[1] = nsq_index;
        }
    }
    p_nsq_8x8[4] = p_nsq_8x8[5] = p_nsq_8x8[6] = p_nsq_8x8[7] = p_nsq_16x16[1];
    // 16x16_2
    // sad[0] = p_sad16x16[2];
    sad[1] = p_sad16x8[4] + p_sad16x8[5];
    sad[2] = p_sad8x16[4] + p_sad8x16[5];
    sad[3] = MAX_SAD_VALUE;
    sad[4] = MAX_SAD_VALUE;
    // sad[5] = p_sad8x8[8] + p_sad8x8[9] + p_sad8x8[10] + p_sad8x8[11];
    best_nsq_sad = MAX_SAD_VALUE;
    for (nsq_index = 1; nsq_index < BLK_NUM; nsq_index++) {
        if (sad[nsq_index] < best_nsq_sad) {
            best_nsq_sad = sad[nsq_index];
            p_nsq_16x16[2] = nsq_index;
        }
    }
    p_nsq_8x8[8] = p_nsq_8x8[9] = p_nsq_8x8[10] = p_nsq_8x8[11] =
        p_nsq_16x16[2];
    // 16x16_3
    // sad[0] = p_sad16x16[3];
    sad[1] = p_sad16x8[6] + p_sad16x8[7];
    sad[2] = p_sad8x16[6] + p_sad8x16[7];
    sad[3] = MAX_SAD_VALUE;
    sad[4] = MAX_SAD_VALUE;
    // sad[5] = p_sad8x8[12] + p_sad8x8[13] + p_sad8x8[14] + p_sad8x8[15];
    best_nsq_sad = MAX_SAD_VALUE;
    for (nsq_index = 1; nsq_index < BLK_NUM; nsq_index++) {
        if (sad[nsq_index] < best_nsq_sad) {
            best_nsq_sad = sad[nsq_index];
            p_nsq_16x16[3] = nsq_index;
        }
    }
    p_nsq_8x8[12] = p_nsq_8x8[13] = p_nsq_8x8[14] = p_nsq_8x8[15] =
        p_nsq_16x16[3];
    // 16x16_4
    // sad[0] = p_sad16x16[4];
    sad[1] = p_sad16x8[8] + p_sad16x8[9];
    sad[2] = p_sad8x16[8] + p_sad8x16[9];
    sad[3] = MAX_SAD_VALUE;
    sad[4] = MAX_SAD_VALUE;
    // sad[5] = p_sad8x8[16] + p_sad8x8[17] + p_sad8x8[18] + p_sad8x8[19];
    best_nsq_sad = MAX_SAD_VALUE;
    for (nsq_index = 1; nsq_index < BLK_NUM; nsq_index++) {
        if (sad[nsq_index] < best_nsq_sad) {
            best_nsq_sad = sad[nsq_index];
            p_nsq_16x16[4] = nsq_index;
        }
    }
    p_nsq_8x8[16] = p_nsq_8x8[17] = p_nsq_8x8[18] = p_nsq_8x8[19] =
        p_nsq_16x16[4];
    // 16x16_5
    // sad[0] = p_sad16x16[5];
    sad[1] = p_sad16x8[10] + p_sad16x8[11];
    sad[2] = p_sad8x16[10] + p_sad8x16[11];
    sad[3] = MAX_SAD_VALUE;
    sad[4] = MAX_SAD_VALUE;
    // sad[5] = p_sad8x8[20] + p_sad8x8[21] + p_sad8x8[22] + p_sad8x8[23];
    best_nsq_sad = MAX_SAD_VALUE;
    for (nsq_index = 1; nsq_index < BLK_NUM; nsq_index++) {
        if (sad[nsq_index] < best_nsq_sad) {
            best_nsq_sad = sad[nsq_index];
            p_nsq_16x16[5] = nsq_index;
        }
    }
    p_nsq_8x8[20] = p_nsq_8x8[21] = p_nsq_8x8[22] = p_nsq_8x8[23] =
        p_nsq_16x16[5];
    // 16x16_6
    // sad[0] = p_sad16x16[6];
    sad[1] = p_sad16x8[12] + p_sad16x8[13];
    sad[2] = p_sad8x16[12] + p_sad8x16[13];
    sad[3] = MAX_SAD_VALUE;
    sad[4] = MAX_SAD_VALUE;
    // sad[5] = p_sad8x8[24] + p_sad8x8[25] + p_sad8x8[26] + p_sad8x8[27];
    best_nsq_sad = MAX_SAD_VALUE;
    for (nsq_index = 1; nsq_index < BLK_NUM; nsq_index++) {
        if (sad[nsq_index] < best_nsq_sad) {
            best_nsq_sad = sad[nsq_index];
            p_nsq_16x16[6] = nsq_index;
        }
    }
    p_nsq_8x8[24] = p_nsq_8x8[25] = p_nsq_8x8[26] = p_nsq_8x8[27] =
        p_nsq_16x16[6];
    // 16x16_7
    // sad[0] = p_sad16x16[7];
    sad[1] = p_sad16x8[14] + p_sad16x8[15];
    sad[2] = p_sad8x16[14] + p_sad8x16[15];
    sad[3] = MAX_SAD_VALUE;
    sad[4] = MAX_SAD_VALUE;
    // sad[5] = p_sad8x8[28] + p_sad8x8[29] + p_sad8x8[30] + p_sad8x8[31];
    best_nsq_sad = MAX_SAD_VALUE;
    for (nsq_index = 1; nsq_index < BLK_NUM; nsq_index++) {
        if (sad[nsq_index] < best_nsq_sad) {
            best_nsq_sad = sad[nsq_index];
            p_nsq_16x16[7] = nsq_index;
        }
    }
    p_nsq_8x8[28] = p_nsq_8x8[29] = p_nsq_8x8[30] = p_nsq_8x8[31] =
        p_nsq_16x16[7];
    // 16x16_8
    // sad[0] = p_sad16x16[8];
    sad[1] = p_sad16x8[16] + p_sad16x8[17];
    sad[2] = p_sad8x16[16] + p_sad8x16[17];
    sad[3] = MAX_SAD_VALUE;
    sad[4] = MAX_SAD_VALUE;
    // sad[5] = p_sad8x8[32] + p_sad8x8[33] + p_sad8x8[34] + p_sad8x8[35];
    best_nsq_sad = MAX_SAD_VALUE;
    for (nsq_index = 1; nsq_index < BLK_NUM; nsq_index++) {
        if (sad[nsq_index] < best_nsq_sad) {
            best_nsq_sad = sad[nsq_index];
            p_nsq_16x16[8] = nsq_index;
        }
    }
    p_nsq_8x8[32] = p_nsq_8x8[33] = p_nsq_8x8[34] = p_nsq_8x8[35] =
        p_nsq_16x16[8];
    // 16x16_9
    // sad[0] = p_sad16x16[9];
    sad[1] = p_sad16x8[18] + p_sad16x8[19];
    sad[2] = p_sad8x16[18] + p_sad8x16[19];
    sad[3] = MAX_SAD_VALUE;
    sad[4] = MAX_SAD_VALUE;
    // sad[5] = p_sad8x8[36] + p_sad8x8[37] + p_sad8x8[38] + p_sad8x8[39];
    best_nsq_sad = MAX_SAD_VALUE;
    for (nsq_index = 1; nsq_index < BLK_NUM; nsq_index++) {
        if (sad[nsq_index] < best_nsq_sad) {
            best_nsq_sad = sad[nsq_index];
            p_nsq_16x16[9] = nsq_index;
        }
    }
    p_nsq_8x8[36] = p_nsq_8x8[37] = p_nsq_8x8[38] = p_nsq_8x8[39] =
        p_nsq_16x16[9];
    // 16x16_10
    // sad[0] = p_sad16x16[10];
    sad[1] = p_sad16x8[20] + p_sad16x8[21];
    sad[2] = p_sad8x16[20] + p_sad8x16[21];
    sad[3] = MAX_SAD_VALUE;
    sad[4] = MAX_SAD_VALUE;
    // sad[5] = p_sad8x8[40] + p_sad8x8[41] + p_sad8x8[42] + p_sad8x8[43];
    best_nsq_sad = MAX_SAD_VALUE;
    for (nsq_index = 1; nsq_index < BLK_NUM; nsq_index++) {
        if (sad[nsq_index] < best_nsq_sad) {
            best_nsq_sad = sad[nsq_index];
            p_nsq_16x16[10] = nsq_index;
        }
    }
    p_nsq_8x8[40] = p_nsq_8x8[41] = p_nsq_8x8[42] = p_nsq_8x8[43] =
        p_nsq_16x16[10];
    // 16x16_11
    // sad[0] = p_sad16x16[11];
    sad[1] = p_sad16x8[22] + p_sad16x8[23];
    sad[2] = p_sad8x16[22] + p_sad8x16[23];
    sad[3] = MAX_SAD_VALUE;
    sad[4] = MAX_SAD_VALUE;
    // sad[5] = p_sad8x8[44] + p_sad8x8[45] + p_sad8x8[46] + p_sad8x8[47];
    best_nsq_sad = MAX_SAD_VALUE;
    for (nsq_index = 1; nsq_index < BLK_NUM; nsq_index++) {
        if (sad[nsq_index] < best_nsq_sad) {
            best_nsq_sad = sad[nsq_index];
            p_nsq_16x16[11] = nsq_index;
        }
    }
    p_nsq_8x8[44] = p_nsq_8x8[45] = p_nsq_8x8[46] = p_nsq_8x8[47] =
        p_nsq_16x16[11];
    // 16x16_12
    // sad[0] = p_sad16x16[12];
    sad[1] = p_sad16x8[24] + p_sad16x8[25];
    sad[2] = p_sad8x16[24] + p_sad8x16[25];
    sad[3] = MAX_SAD_VALUE;
    sad[4] = MAX_SAD_VALUE;
    // sad[5] = p_sad8x8[48] + p_sad8x8[49] + p_sad8x8[50] + p_sad8x8[51];
    best_nsq_sad = MAX_SAD_VALUE;
    for (nsq_index = 1; nsq_index < BLK_NUM; nsq_index++) {
        if (sad[nsq_index] < best_nsq_sad) {
            best_nsq_sad = sad[nsq_index];
            p_nsq_16x16[12] = nsq_index;
        }
    }
    p_nsq_8x8[48] = p_nsq_8x8[49] = p_nsq_8x8[50] = p_nsq_8x8[51] =
        p_nsq_16x16[12];
    // 16x16_13
    // sad[0] = p_sad16x16[13];
    sad[1] = p_sad16x8[26] + p_sad16x8[27];
    sad[2] = p_sad8x16[26] + p_sad8x16[27];
    sad[3] = MAX_SAD_VALUE;
    sad[4] = MAX_SAD_VALUE;
    // sad[5] = p_sad8x8[52] + p_sad8x8[53] + p_sad8x8[54] + p_sad8x8[55];
    best_nsq_sad = MAX_SAD_VALUE;
    for (nsq_index = 1; nsq_index < BLK_NUM; nsq_index++) {
        if (sad[nsq_index] < best_nsq_sad) {
            best_nsq_sad = sad[nsq_index];
            p_nsq_16x16[13] = nsq_index;
        }
    }
    p_nsq_8x8[52] = p_nsq_8x8[53] = p_nsq_8x8[54] = p_nsq_8x8[55] =
        p_nsq_16x16[13];
    // 16x16_14
    // sad[0] = p_sad16x16[14];
    sad[1] = p_sad16x8[28] + p_sad16x8[29];
    sad[2] = p_sad8x16[28] + p_sad8x16[29];
    sad[3] = MAX_SAD_VALUE;
    sad[4] = MAX_SAD_VALUE;
    // sad[5] = p_sad8x8[56] + p_sad8x8[57] + p_sad8x8[58] + p_sad8x8[59];
    best_nsq_sad = MAX_SAD_VALUE;
    for (nsq_index = 1; nsq_index < BLK_NUM; nsq_index++) {
        if (sad[nsq_index] < best_nsq_sad) {
            best_nsq_sad = sad[nsq_index];
            p_nsq_16x16[14] = nsq_index;
        }
    }
    p_nsq_8x8[56] = p_nsq_8x8[57] = p_nsq_8x8[58] = p_nsq_8x8[59] =
        p_nsq_16x16[14];
    // 16x16_15
    // sad[0] = p_sad16x16[15];
    sad[1] = p_sad16x8[30] + p_sad16x8[31];
    sad[2] = p_sad8x16[30] + p_sad8x16[31];
    sad[3] = MAX_SAD_VALUE;
    sad[4] = MAX_SAD_VALUE;
    // sad[5] = p_sad8x8[60] + p_sad8x8[61] + p_sad8x8[62] + p_sad8x8[63];
    best_nsq_sad = MAX_SAD_VALUE;
    for (nsq_index = 1; nsq_index < BLK_NUM; nsq_index++) {
        if (sad[nsq_index] < best_nsq_sad) {
            best_nsq_sad = sad[nsq_index];
            p_nsq_16x16[15] = nsq_index;
        }
    }
    p_nsq_8x8[60] = p_nsq_8x8[61] = p_nsq_8x8[62] = p_nsq_8x8[63] =
        p_nsq_16x16[15];
}

/****************************************************
Calcualte SAD for Rect H, V and H4, V4 partitions

and update its Motion info if the result SAD is better
****************************************************/
void ExtSadCalculation(uint32_t *p_sad8x8, uint32_t *p_sad16x16,
                       uint32_t *p_sad32x32, uint32_t *p_best_sad64x32,
                       uint32_t *p_best_mv64x32, uint32_t *p_best_sad32x16,
                       uint32_t *p_best_mv32x16, uint32_t *p_best_sad16x8,
                       uint32_t *p_best_mv16x8, uint32_t *p_best_sad32x64,
                       uint32_t *p_best_mv32x64, uint32_t *p_best_sad16x32,
                       uint32_t *p_best_mv16x32, uint32_t *p_best_sad8x16,
                       uint32_t *p_best_mv8x16, uint32_t *p_best_sad32x8,
                       uint32_t *p_best_mv32x8, uint32_t *p_best_sad8x32,
                       uint32_t *p_best_mv8x32, uint32_t *p_best_sad64x16,
                       uint32_t *p_best_mv64x16, uint32_t *p_best_sad16x64,
                       uint32_t *p_best_mv16x64, uint32_t mv) {
    uint32_t sad;

    uint32_t sad_16x8[32];
    uint32_t sad_8x16[32];
    uint32_t sad_32x16[8];
    uint32_t sad_16x32[8];

    // 64x32
    sad = p_sad32x32[0] + p_sad32x32[1];
    if (sad < p_best_sad64x32[0]) {
        p_best_sad64x32[0] = sad;
        p_best_mv64x32[0] = mv;
    }

    sad = p_sad32x32[2] + p_sad32x32[3];
    if (sad < p_best_sad64x32[1]) {
        p_best_sad64x32[1] = sad;
        p_best_mv64x32[1] = mv;
    }

    // 32x16
    sad_32x16[0] = p_sad16x16[0] + p_sad16x16[1];
    if (sad_32x16[0] < p_best_sad32x16[0]) {
        p_best_sad32x16[0] = sad_32x16[0];
        p_best_mv32x16[0] = mv;
    }

    sad_32x16[1] = p_sad16x16[2] + p_sad16x16[3];
    if (sad_32x16[1] < p_best_sad32x16[1]) {
        p_best_sad32x16[1] = sad_32x16[1];
        p_best_mv32x16[1] = mv;
    }

    sad_32x16[2] = p_sad16x16[4] + p_sad16x16[5];
    if (sad_32x16[2] < p_best_sad32x16[2]) {
        p_best_sad32x16[2] = sad_32x16[2];
        p_best_mv32x16[2] = mv;
    }

    sad_32x16[3] = p_sad16x16[6] + p_sad16x16[7];
    if (sad_32x16[3] < p_best_sad32x16[3]) {
        p_best_sad32x16[3] = sad_32x16[3];
        p_best_mv32x16[3] = mv;
    }

    sad_32x16[4] = p_sad16x16[8] + p_sad16x16[9];
    if (sad_32x16[4] < p_best_sad32x16[4]) {
        p_best_sad32x16[4] = sad_32x16[4];
        p_best_mv32x16[4] = mv;
    }

    sad_32x16[5] = p_sad16x16[10] + p_sad16x16[11];
    if (sad_32x16[5] < p_best_sad32x16[5]) {
        p_best_sad32x16[5] = sad_32x16[5];
        p_best_mv32x16[5] = mv;
    }

    sad_32x16[6] = p_sad16x16[12] + p_sad16x16[13];
    if (sad_32x16[6] < p_best_sad32x16[6]) {
        p_best_sad32x16[6] = sad_32x16[6];
        p_best_mv32x16[6] = mv;
    }

    sad_32x16[7] = p_sad16x16[14] + p_sad16x16[15];
    if (sad_32x16[7] < p_best_sad32x16[7]) {
        p_best_sad32x16[7] = sad_32x16[7];
        p_best_mv32x16[7] = mv;
    }

    // 64x16
    sad = sad_32x16[0] + sad_32x16[2];
    if (sad < p_best_sad64x16[0]) {
        p_best_sad64x16[0] = sad;
        p_best_mv64x16[0] = mv;
    }
    sad = sad_32x16[1] + sad_32x16[3];
    if (sad < p_best_sad64x16[1]) {
        p_best_sad64x16[1] = sad;
        p_best_mv64x16[1] = mv;
    }

    sad = sad_32x16[4] + sad_32x16[6];
    if (sad < p_best_sad64x16[2]) {
        p_best_sad64x16[2] = sad;
        p_best_mv64x16[2] = mv;
    }
    sad = sad_32x16[5] + sad_32x16[7];
    if (sad < p_best_sad64x16[3]) {
        p_best_sad64x16[3] = sad;
        p_best_mv64x16[3] = mv;
    }

    // 16x8
    sad_16x8[0] = p_sad8x8[0] + p_sad8x8[1];
    if (sad_16x8[0] < p_best_sad16x8[0]) {
        p_best_sad16x8[0] = sad_16x8[0];
        p_best_mv16x8[0] = mv;
    }

    sad_16x8[1] = p_sad8x8[2] + p_sad8x8[3];
    if (sad_16x8[1] < p_best_sad16x8[1]) {
        p_best_sad16x8[1] = sad_16x8[1];
        p_best_mv16x8[1] = mv;
    }

    sad_16x8[2] = p_sad8x8[4] + p_sad8x8[5];
    if (sad_16x8[2] < p_best_sad16x8[2]) {
        p_best_sad16x8[2] = sad_16x8[2];
        p_best_mv16x8[2] = mv;
    }

    sad_16x8[3] = p_sad8x8[6] + p_sad8x8[7];
    if (sad_16x8[3] < p_best_sad16x8[3]) {
        p_best_sad16x8[3] = sad_16x8[3];
        p_best_mv16x8[3] = mv;
    }

    sad_16x8[4] = p_sad8x8[8] + p_sad8x8[9];
    if (sad_16x8[4] < p_best_sad16x8[4]) {
        p_best_sad16x8[4] = sad_16x8[4];
        p_best_mv16x8[4] = mv;
    }

    sad_16x8[5] = p_sad8x8[10] + p_sad8x8[11];
    if (sad_16x8[5] < p_best_sad16x8[5]) {
        p_best_sad16x8[5] = sad_16x8[5];
        p_best_mv16x8[5] = mv;
    }

    sad_16x8[6] = p_sad8x8[12] + p_sad8x8[13];
    if (sad_16x8[6] < p_best_sad16x8[6]) {
        p_best_sad16x8[6] = sad_16x8[6];
        p_best_mv16x8[6] = mv;
    }

    sad_16x8[7] = p_sad8x8[14] + p_sad8x8[15];
    if (sad_16x8[7] < p_best_sad16x8[7]) {
        p_best_sad16x8[7] = sad_16x8[7];
        p_best_mv16x8[7] = mv;
    }

    sad_16x8[8] = p_sad8x8[16] + p_sad8x8[17];
    if (sad_16x8[8] < p_best_sad16x8[8]) {
        p_best_sad16x8[8] = sad_16x8[8];
        p_best_mv16x8[8] = mv;
    }

    sad_16x8[9] = p_sad8x8[18] + p_sad8x8[19];
    if (sad_16x8[9] < p_best_sad16x8[9]) {
        p_best_sad16x8[9] = sad_16x8[9];
        p_best_mv16x8[9] = mv;
    }

    sad_16x8[10] = p_sad8x8[20] + p_sad8x8[21];
    if (sad_16x8[10] < p_best_sad16x8[10]) {
        p_best_sad16x8[10] = sad_16x8[10];
        p_best_mv16x8[10] = mv;
    }

    sad_16x8[11] = p_sad8x8[22] + p_sad8x8[23];
    if (sad_16x8[11] < p_best_sad16x8[11]) {
        p_best_sad16x8[11] = sad_16x8[11];
        p_best_mv16x8[11] = mv;
    }

    sad_16x8[12] = p_sad8x8[24] + p_sad8x8[25];
    if (sad_16x8[12] < p_best_sad16x8[12]) {
        p_best_sad16x8[12] = sad_16x8[12];
        p_best_mv16x8[12] = mv;
    }

    sad_16x8[13] = p_sad8x8[26] + p_sad8x8[27];
    if (sad_16x8[13] < p_best_sad16x8[13]) {
        p_best_sad16x8[13] = sad_16x8[13];
        p_best_mv16x8[13] = mv;
    }

    sad_16x8[14] = p_sad8x8[28] + p_sad8x8[29];
    if (sad_16x8[14] < p_best_sad16x8[14]) {
        p_best_sad16x8[14] = sad_16x8[14];
        p_best_mv16x8[14] = mv;
    }

    sad_16x8[15] = p_sad8x8[30] + p_sad8x8[31];
    if (sad_16x8[15] < p_best_sad16x8[15]) {
        p_best_sad16x8[15] = sad_16x8[15];
        p_best_mv16x8[15] = mv;
    }

    sad_16x8[16] = p_sad8x8[32] + p_sad8x8[33];
    if (sad_16x8[16] < p_best_sad16x8[16]) {
        p_best_sad16x8[16] = sad_16x8[16];
        p_best_mv16x8[16] = mv;
    }

    sad_16x8[17] = p_sad8x8[34] + p_sad8x8[35];
    if (sad_16x8[17] < p_best_sad16x8[17]) {
        p_best_sad16x8[17] = sad_16x8[17];
        p_best_mv16x8[17] = mv;
    }

    sad_16x8[18] = p_sad8x8[36] + p_sad8x8[37];
    if (sad_16x8[18] < p_best_sad16x8[18]) {
        p_best_sad16x8[18] = sad_16x8[18];
        p_best_mv16x8[18] = mv;
    }

    sad_16x8[19] = p_sad8x8[38] + p_sad8x8[39];
    if (sad_16x8[19] < p_best_sad16x8[19]) {
        p_best_sad16x8[19] = sad_16x8[19];
        p_best_mv16x8[19] = mv;
    }

    sad_16x8[20] = p_sad8x8[40] + p_sad8x8[41];
    if (sad_16x8[20] < p_best_sad16x8[20]) {
        p_best_sad16x8[20] = sad_16x8[20];
        p_best_mv16x8[20] = mv;
    }

    sad_16x8[21] = p_sad8x8[42] + p_sad8x8[43];
    if (sad_16x8[21] < p_best_sad16x8[21]) {
        p_best_sad16x8[21] = sad_16x8[21];
        p_best_mv16x8[21] = mv;
    }

    sad_16x8[22] = p_sad8x8[44] + p_sad8x8[45];
    if (sad_16x8[22] < p_best_sad16x8[22]) {
        p_best_sad16x8[22] = sad_16x8[22];
        p_best_mv16x8[22] = mv;
    }

    sad_16x8[23] = p_sad8x8[46] + p_sad8x8[47];
    if (sad_16x8[23] < p_best_sad16x8[23]) {
        p_best_sad16x8[23] = sad_16x8[23];
        p_best_mv16x8[23] = mv;
    }

    sad_16x8[24] = p_sad8x8[48] + p_sad8x8[49];
    if (sad_16x8[24] < p_best_sad16x8[24]) {
        p_best_sad16x8[24] = sad_16x8[24];
        p_best_mv16x8[24] = mv;
    }

    sad_16x8[25] = p_sad8x8[50] + p_sad8x8[51];
    if (sad_16x8[25] < p_best_sad16x8[25]) {
        p_best_sad16x8[25] = sad_16x8[25];
        p_best_mv16x8[25] = mv;
    }

    sad_16x8[26] = p_sad8x8[52] + p_sad8x8[53];
    if (sad_16x8[26] < p_best_sad16x8[26]) {
        p_best_sad16x8[26] = sad_16x8[26];
        p_best_mv16x8[26] = mv;
    }

    sad_16x8[27] = p_sad8x8[54] + p_sad8x8[55];
    if (sad_16x8[27] < p_best_sad16x8[27]) {
        p_best_sad16x8[27] = sad_16x8[27];
        p_best_mv16x8[27] = mv;
    }

    sad_16x8[28] = p_sad8x8[56] + p_sad8x8[57];
    if (sad_16x8[28] < p_best_sad16x8[28]) {
        p_best_sad16x8[28] = sad_16x8[28];
        p_best_mv16x8[28] = mv;
    }

    sad_16x8[29] = p_sad8x8[58] + p_sad8x8[59];
    if (sad_16x8[29] < p_best_sad16x8[29]) {
        p_best_sad16x8[29] = sad_16x8[29];
        p_best_mv16x8[29] = mv;
    }

    sad_16x8[30] = p_sad8x8[60] + p_sad8x8[61];
    if (sad_16x8[30] < p_best_sad16x8[30]) {
        p_best_sad16x8[30] = sad_16x8[30];
        p_best_mv16x8[30] = mv;
    }

    sad_16x8[31] = p_sad8x8[62] + p_sad8x8[63];
    if (sad_16x8[31] < p_best_sad16x8[31]) {
        p_best_sad16x8[31] = sad_16x8[31];
        p_best_mv16x8[31] = mv;
    }

    // 32x64
    sad = p_sad32x32[0] + p_sad32x32[2];
    if (sad < p_best_sad32x64[0]) {
        p_best_sad32x64[0] = sad;
        p_best_mv32x64[0] = mv;
    }

    sad = p_sad32x32[1] + p_sad32x32[3];
    if (sad < p_best_sad32x64[1]) {
        p_best_sad32x64[1] = sad;
        p_best_mv32x64[1] = mv;
    }

    // 16x32
    sad_16x32[0] = p_sad16x16[0] + p_sad16x16[2];
    if (sad_16x32[0] < p_best_sad16x32[0]) {
        p_best_sad16x32[0] = sad_16x32[0];
        p_best_mv16x32[0] = mv;
    }

    sad_16x32[1] = p_sad16x16[1] + p_sad16x16[3];
    if (sad_16x32[1] < p_best_sad16x32[1]) {
        p_best_sad16x32[1] = sad_16x32[1];
        p_best_mv16x32[1] = mv;
    }

    sad_16x32[2] = p_sad16x16[4] + p_sad16x16[6];
    if (sad_16x32[2] < p_best_sad16x32[2]) {
        p_best_sad16x32[2] = sad_16x32[2];
        p_best_mv16x32[2] = mv;
    }

    sad_16x32[3] = p_sad16x16[5] + p_sad16x16[7];
    if (sad_16x32[3] < p_best_sad16x32[3]) {
        p_best_sad16x32[3] = sad_16x32[3];
        p_best_mv16x32[3] = mv;
    }

    sad_16x32[4] = p_sad16x16[8] + p_sad16x16[10];
    if (sad_16x32[4] < p_best_sad16x32[4]) {
        p_best_sad16x32[4] = sad_16x32[4];
        p_best_mv16x32[4] = mv;
    }

    sad_16x32[5] = p_sad16x16[9] + p_sad16x16[11];
    if (sad_16x32[5] < p_best_sad16x32[5]) {
        p_best_sad16x32[5] = sad_16x32[5];
        p_best_mv16x32[5] = mv;
    }

    sad_16x32[6] = p_sad16x16[12] + p_sad16x16[14];
    if (sad_16x32[6] < p_best_sad16x32[6]) {
        p_best_sad16x32[6] = sad_16x32[6];
        p_best_mv16x32[6] = mv;
    }

    sad_16x32[7] = p_sad16x16[13] + p_sad16x16[15];
    if (sad_16x32[7] < p_best_sad16x32[7]) {
        p_best_sad16x32[7] = sad_16x32[7];
        p_best_mv16x32[7] = mv;
    }

    sad = sad_16x32[0] + sad_16x32[4];
    if (sad < p_best_sad16x64[0]) {
        p_best_sad16x64[0] = sad;
        p_best_mv16x64[0] = mv;
    }
    sad = sad_16x32[1] + sad_16x32[5];
    if (sad < p_best_sad16x64[1]) {
        p_best_sad16x64[1] = sad;
        p_best_mv16x64[1] = mv;
    }

    sad = sad_16x32[2] + sad_16x32[6];
    if (sad < p_best_sad16x64[2]) {
        p_best_sad16x64[2] = sad;
        p_best_mv16x64[2] = mv;
    }

    sad = sad_16x32[3] + sad_16x32[7];
    if (sad < p_best_sad16x64[3]) {
        p_best_sad16x64[3] = sad;
        p_best_mv16x64[3] = mv;
    }

    // 8x16
    sad_8x16[0] = p_sad8x8[0] + p_sad8x8[2];
    if (sad_8x16[0] < p_best_sad8x16[0]) {
        p_best_sad8x16[0] = sad_8x16[0];
        p_best_mv8x16[0] = mv;
    }

    sad_8x16[1] = p_sad8x8[1] + p_sad8x8[3];
    if (sad_8x16[1] < p_best_sad8x16[1]) {
        p_best_sad8x16[1] = sad_8x16[1];
        p_best_mv8x16[1] = mv;
    }

    sad_8x16[2] = p_sad8x8[4] + p_sad8x8[6];
    if (sad_8x16[2] < p_best_sad8x16[2]) {
        p_best_sad8x16[2] = sad_8x16[2];
        p_best_mv8x16[2] = mv;
    }

    sad_8x16[3] = p_sad8x8[5] + p_sad8x8[7];
    if (sad_8x16[3] < p_best_sad8x16[3]) {
        p_best_sad8x16[3] = sad_8x16[3];
        p_best_mv8x16[3] = mv;
    }

    sad_8x16[4] = p_sad8x8[8] + p_sad8x8[10];
    if (sad_8x16[4] < p_best_sad8x16[4]) {
        p_best_sad8x16[4] = sad_8x16[4];
        p_best_mv8x16[4] = mv;
    }

    sad_8x16[5] = p_sad8x8[9] + p_sad8x8[11];
    if (sad_8x16[5] < p_best_sad8x16[5]) {
        p_best_sad8x16[5] = sad_8x16[5];
        p_best_mv8x16[5] = mv;
    }

    sad_8x16[6] = p_sad8x8[12] + p_sad8x8[14];
    if (sad_8x16[6] < p_best_sad8x16[6]) {
        p_best_sad8x16[6] = sad_8x16[6];
        p_best_mv8x16[6] = mv;
    }

    sad_8x16[7] = p_sad8x8[13] + p_sad8x8[15];
    if (sad_8x16[7] < p_best_sad8x16[7]) {
        p_best_sad8x16[7] = sad_8x16[7];
        p_best_mv8x16[7] = mv;
    }

    sad_8x16[8] = p_sad8x8[16] + p_sad8x8[18];
    if (sad_8x16[8] < p_best_sad8x16[8]) {
        p_best_sad8x16[8] = sad_8x16[8];
        p_best_mv8x16[8] = mv;
    }

    sad_8x16[9] = p_sad8x8[17] + p_sad8x8[19];
    if (sad_8x16[9] < p_best_sad8x16[9]) {
        p_best_sad8x16[9] = sad_8x16[9];
        p_best_mv8x16[9] = mv;
    }

    sad_8x16[10] = p_sad8x8[20] + p_sad8x8[22];
    if (sad_8x16[10] < p_best_sad8x16[10]) {
        p_best_sad8x16[10] = sad_8x16[10];
        p_best_mv8x16[10] = mv;
    }

    sad_8x16[11] = p_sad8x8[21] + p_sad8x8[23];
    if (sad_8x16[11] < p_best_sad8x16[11]) {
        p_best_sad8x16[11] = sad_8x16[11];
        p_best_mv8x16[11] = mv;
    }

    sad_8x16[12] = p_sad8x8[24] + p_sad8x8[26];
    if (sad_8x16[12] < p_best_sad8x16[12]) {
        p_best_sad8x16[12] = sad_8x16[12];
        p_best_mv8x16[12] = mv;
    }

    sad_8x16[13] = p_sad8x8[25] + p_sad8x8[27];
    if (sad_8x16[13] < p_best_sad8x16[13]) {
        p_best_sad8x16[13] = sad_8x16[13];
        p_best_mv8x16[13] = mv;
    }

    sad_8x16[14] = p_sad8x8[28] + p_sad8x8[30];
    if (sad_8x16[14] < p_best_sad8x16[14]) {
        p_best_sad8x16[14] = sad_8x16[14];
        p_best_mv8x16[14] = mv;
    }

    sad_8x16[15] = p_sad8x8[29] + p_sad8x8[31];
    if (sad_8x16[15] < p_best_sad8x16[15]) {
        p_best_sad8x16[15] = sad_8x16[15];
        p_best_mv8x16[15] = mv;
    }

    sad_8x16[16] = p_sad8x8[32] + p_sad8x8[34];
    if (sad_8x16[16] < p_best_sad8x16[16]) {
        p_best_sad8x16[16] = sad_8x16[16];
        p_best_mv8x16[16] = mv;
    }

    sad_8x16[17] = p_sad8x8[33] + p_sad8x8[35];
    if (sad_8x16[17] < p_best_sad8x16[17]) {
        p_best_sad8x16[17] = sad_8x16[17];
        p_best_mv8x16[17] = mv;
    }

    sad_8x16[18] = p_sad8x8[36] + p_sad8x8[38];
    if (sad_8x16[18] < p_best_sad8x16[18]) {
        p_best_sad8x16[18] = sad_8x16[18];
        p_best_mv8x16[18] = mv;
    }

    sad_8x16[19] = p_sad8x8[37] + p_sad8x8[39];
    if (sad_8x16[19] < p_best_sad8x16[19]) {
        p_best_sad8x16[19] = sad_8x16[19];
        p_best_mv8x16[19] = mv;
    }

    sad_8x16[20] = p_sad8x8[40] + p_sad8x8[42];
    if (sad_8x16[20] < p_best_sad8x16[20]) {
        p_best_sad8x16[20] = sad_8x16[20];
        p_best_mv8x16[20] = mv;
    }

    sad_8x16[21] = p_sad8x8[41] + p_sad8x8[43];
    if (sad_8x16[21] < p_best_sad8x16[21]) {
        p_best_sad8x16[21] = sad_8x16[21];
        p_best_mv8x16[21] = mv;
    }

    sad_8x16[22] = p_sad8x8[44] + p_sad8x8[46];
    if (sad_8x16[22] < p_best_sad8x16[22]) {
        p_best_sad8x16[22] = sad_8x16[22];
        p_best_mv8x16[22] = mv;
    }

    sad_8x16[23] = p_sad8x8[45] + p_sad8x8[47];
    if (sad_8x16[23] < p_best_sad8x16[23]) {
        p_best_sad8x16[23] = sad_8x16[23];
        p_best_mv8x16[23] = mv;
    }

    sad_8x16[24] = p_sad8x8[48] + p_sad8x8[50];
    if (sad_8x16[24] < p_best_sad8x16[24]) {
        p_best_sad8x16[24] = sad_8x16[24];
        p_best_mv8x16[24] = mv;
    }

    sad_8x16[25] = p_sad8x8[49] + p_sad8x8[51];
    if (sad_8x16[25] < p_best_sad8x16[25]) {
        p_best_sad8x16[25] = sad_8x16[25];
        p_best_mv8x16[25] = mv;
    }

    sad_8x16[26] = p_sad8x8[52] + p_sad8x8[54];
    if (sad_8x16[26] < p_best_sad8x16[26]) {
        p_best_sad8x16[26] = sad_8x16[26];
        p_best_mv8x16[26] = mv;
    }

    sad_8x16[27] = p_sad8x8[53] + p_sad8x8[55];
    if (sad_8x16[27] < p_best_sad8x16[27]) {
        p_best_sad8x16[27] = sad_8x16[27];
        p_best_mv8x16[27] = mv;
    }

    sad_8x16[28] = p_sad8x8[56] + p_sad8x8[58];
    if (sad_8x16[28] < p_best_sad8x16[28]) {
        p_best_sad8x16[28] = sad_8x16[28];
        p_best_mv8x16[28] = mv;
    }

    sad_8x16[29] = p_sad8x8[57] + p_sad8x8[59];
    if (sad_8x16[29] < p_best_sad8x16[29]) {
        p_best_sad8x16[29] = sad_8x16[29];
        p_best_mv8x16[29] = mv;
    }

    sad_8x16[30] = p_sad8x8[60] + p_sad8x8[62];
    if (sad_8x16[30] < p_best_sad8x16[30]) {
        p_best_sad8x16[30] = sad_8x16[30];
        p_best_mv8x16[30] = mv;
    }

    sad_8x16[31] = p_sad8x8[61] + p_sad8x8[63];
    if (sad_8x16[31] < p_best_sad8x16[31]) {
        p_best_sad8x16[31] = sad_8x16[31];
        p_best_mv8x16[31] = mv;
    }

    // 32x8
    sad = sad_16x8[0] + sad_16x8[2];
    if (sad < p_best_sad32x8[0]) {
        p_best_sad32x8[0] = sad;
        p_best_mv32x8[0] = mv;
    }

    sad = sad_16x8[1] + sad_16x8[3];
    if (sad < p_best_sad32x8[1]) {
        p_best_sad32x8[1] = sad;
        p_best_mv32x8[1] = mv;
    }

    sad = sad_16x8[4] + sad_16x8[6];
    if (sad < p_best_sad32x8[2]) {
        p_best_sad32x8[2] = sad;
        p_best_mv32x8[2] = mv;
    }

    sad = sad_16x8[5] + sad_16x8[7];
    if (sad < p_best_sad32x8[3]) {
        p_best_sad32x8[3] = sad;
        p_best_mv32x8[3] = mv;
    }

    sad = sad_16x8[8] + sad_16x8[10];
    if (sad < p_best_sad32x8[4]) {
        p_best_sad32x8[4] = sad;
        p_best_mv32x8[4] = mv;
    }

    sad = sad_16x8[9] + sad_16x8[11];
    if (sad < p_best_sad32x8[5]) {
        p_best_sad32x8[5] = sad;
        p_best_mv32x8[5] = mv;
    }

    sad = sad_16x8[12] + sad_16x8[14];
    if (sad < p_best_sad32x8[6]) {
        p_best_sad32x8[6] = sad;
        p_best_mv32x8[6] = mv;
    }

    sad = sad_16x8[13] + sad_16x8[15];
    if (sad < p_best_sad32x8[7]) {
        p_best_sad32x8[7] = sad;
        p_best_mv32x8[7] = mv;
    }

    sad = sad_16x8[16] + sad_16x8[18];
    if (sad < p_best_sad32x8[8]) {
        p_best_sad32x8[8] = sad;
        p_best_mv32x8[8] = mv;
    }

    sad = sad_16x8[17] + sad_16x8[19];
    if (sad < p_best_sad32x8[9]) {
        p_best_sad32x8[9] = sad;
        p_best_mv32x8[9] = mv;
    }

    sad = sad_16x8[20] + sad_16x8[22];
    if (sad < p_best_sad32x8[10]) {
        p_best_sad32x8[10] = sad;
        p_best_mv32x8[10] = mv;
    }

    sad = sad_16x8[21] + sad_16x8[23];
    if (sad < p_best_sad32x8[11]) {
        p_best_sad32x8[11] = sad;
        p_best_mv32x8[11] = mv;
    }

    sad = sad_16x8[24] + sad_16x8[26];
    if (sad < p_best_sad32x8[12]) {
        p_best_sad32x8[12] = sad;
        p_best_mv32x8[12] = mv;
    }

    sad = sad_16x8[25] + sad_16x8[27];
    if (sad < p_best_sad32x8[13]) {
        p_best_sad32x8[13] = sad;
        p_best_mv32x8[13] = mv;
    }

    sad = sad_16x8[28] + sad_16x8[30];
    if (sad < p_best_sad32x8[14]) {
        p_best_sad32x8[14] = sad;
        p_best_mv32x8[14] = mv;
    }

    sad = sad_16x8[29] + sad_16x8[31];
    if (sad < p_best_sad32x8[15]) {
        p_best_sad32x8[15] = sad;
        p_best_mv32x8[15] = mv;
    }

    // 8x32
    sad = sad_8x16[0] + sad_8x16[4];
    if (sad < p_best_sad8x32[0]) {
        p_best_sad8x32[0] = sad;
        p_best_mv8x32[0] = mv;
    }

    sad = sad_8x16[1] + sad_8x16[5];
    if (sad < p_best_sad8x32[1]) {
        p_best_sad8x32[1] = sad;
        p_best_mv8x32[1] = mv;
    }

    sad = sad_8x16[2] + sad_8x16[6];
    if (sad < p_best_sad8x32[2]) {
        p_best_sad8x32[2] = sad;
        p_best_mv8x32[2] = mv;
    }

    sad = sad_8x16[3] + sad_8x16[7];
    if (sad < p_best_sad8x32[3]) {
        p_best_sad8x32[3] = sad;
        p_best_mv8x32[3] = mv;
    }

    sad = sad_8x16[8] + sad_8x16[12];
    if (sad < p_best_sad8x32[4]) {
        p_best_sad8x32[4] = sad;
        p_best_mv8x32[4] = mv;
    }

    sad = sad_8x16[9] + sad_8x16[13];
    if (sad < p_best_sad8x32[5]) {
        p_best_sad8x32[5] = sad;
        p_best_mv8x32[5] = mv;
    }

    sad = sad_8x16[10] + sad_8x16[14];
    if (sad < p_best_sad8x32[6]) {
        p_best_sad8x32[6] = sad;
        p_best_mv8x32[6] = mv;
    }

    sad = sad_8x16[11] + sad_8x16[15];
    if (sad < p_best_sad8x32[7]) {
        p_best_sad8x32[7] = sad;
        p_best_mv8x32[7] = mv;
    }

    sad = sad_8x16[16] + sad_8x16[20];
    if (sad < p_best_sad8x32[8]) {
        p_best_sad8x32[8] = sad;
        p_best_mv8x32[8] = mv;
    }

    sad = sad_8x16[17] + sad_8x16[21];
    if (sad < p_best_sad8x32[9]) {
        p_best_sad8x32[9] = sad;
        p_best_mv8x32[9] = mv;
    }

    sad = sad_8x16[18] + sad_8x16[22];
    if (sad < p_best_sad8x32[10]) {
        p_best_sad8x32[10] = sad;
        p_best_mv8x32[10] = mv;
    }

    sad = sad_8x16[19] + sad_8x16[23];
    if (sad < p_best_sad8x32[11]) {
        p_best_sad8x32[11] = sad;
        p_best_mv8x32[11] = mv;
    }

    sad = sad_8x16[24] + sad_8x16[28];
    if (sad < p_best_sad8x32[12]) {
        p_best_sad8x32[12] = sad;
        p_best_mv8x32[12] = mv;
    }

    sad = sad_8x16[25] + sad_8x16[29];
    if (sad < p_best_sad8x32[13]) {
        p_best_sad8x32[13] = sad;
        p_best_mv8x32[13] = mv;
    }

    sad = sad_8x16[26] + sad_8x16[30];
    if (sad < p_best_sad8x32[14]) {
        p_best_sad8x32[14] = sad;
        p_best_mv8x32[14] = mv;
    }

    sad = sad_8x16[27] + sad_8x16[31];
    if (sad < p_best_sad8x32[15]) {
        p_best_sad8x32[15] = sad;
        p_best_mv8x32[15] = mv;
    }
}

/****************************************************
Calcualte SAD for Rect H, V and H4, V4 partitions
and update its Motion info if the result SAD is better
****************************************************/
void ext_eigth_sad_calculation_nsq_c(
    uint32_t p_sad8x8[64][8], uint32_t p_sad16x16[16][8],
    uint32_t p_sad32x32[4][8], uint32_t *p_best_sad64x32,
    uint32_t *p_best_mv64x32, uint32_t *p_best_sad32x16,
    uint32_t *p_best_mv32x16, uint32_t *p_best_sad16x8, uint32_t *p_best_mv16x8,
    uint32_t *p_best_sad32x64, uint32_t *p_best_mv32x64,
    uint32_t *p_best_sad16x32, uint32_t *p_best_mv16x32,
    uint32_t *p_best_sad8x16, uint32_t *p_best_mv8x16, uint32_t *p_best_sad32x8,
    uint32_t *p_best_mv32x8, uint32_t *p_best_sad8x32, uint32_t *p_best_mv8x32,
    uint32_t *p_best_sad64x16, uint32_t *p_best_mv64x16,
    uint32_t *p_best_sad16x64, uint32_t *p_best_mv16x64, uint32_t mv) {
    uint8_t search_index;
    uint32_t sad;
    uint32_t sad_16x8[32];
    uint32_t sad_8x16[32];
    uint32_t sad_32x16[8];
    uint32_t sad_16x32[8];

    int16_t x_mv, y_mv;

    for (search_index = 0; search_index < 8; search_index++) {
        // 64x32
        sad = p_sad32x32[0][search_index] + p_sad32x32[1][search_index];
        if (sad < p_best_sad64x32[0]) {
            p_best_sad64x32[0] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv64x32[0] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = p_sad32x32[2][search_index] + p_sad32x32[3][search_index];
        if (sad < p_best_sad64x32[1]) {
            p_best_sad64x32[1] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv64x32[1] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        // 32x16
        sad_32x16[0] =
            p_sad16x16[0][search_index] + p_sad16x16[1][search_index];
        if (sad_32x16[0] < p_best_sad32x16[0]) {
            p_best_sad32x16[0] = sad_32x16[0];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv32x16[0] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_32x16[1] =
            p_sad16x16[2][search_index] + p_sad16x16[3][search_index];
        if (sad_32x16[1] < p_best_sad32x16[1]) {
            p_best_sad32x16[1] = sad_32x16[1];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv32x16[1] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_32x16[2] =
            p_sad16x16[4][search_index] + p_sad16x16[5][search_index];
        if (sad_32x16[2] < p_best_sad32x16[2]) {
            p_best_sad32x16[2] = sad_32x16[2];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv32x16[2] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_32x16[3] =
            p_sad16x16[6][search_index] + p_sad16x16[7][search_index];
        if (sad_32x16[3] < p_best_sad32x16[3]) {
            p_best_sad32x16[3] = sad_32x16[3];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv32x16[3] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_32x16[4] =
            p_sad16x16[8][search_index] + p_sad16x16[9][search_index];
        if (sad_32x16[4] < p_best_sad32x16[4]) {
            p_best_sad32x16[4] = sad_32x16[4];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv32x16[4] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_32x16[5] =
            p_sad16x16[10][search_index] + p_sad16x16[11][search_index];
        if (sad_32x16[5] < p_best_sad32x16[5]) {
            p_best_sad32x16[5] = sad_32x16[5];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv32x16[5] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_32x16[6] =
            p_sad16x16[12][search_index] + p_sad16x16[13][search_index];
        if (sad_32x16[6] < p_best_sad32x16[6]) {
            p_best_sad32x16[6] = sad_32x16[6];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv32x16[6] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_32x16[7] =
            p_sad16x16[14][search_index] + p_sad16x16[15][search_index];
        if (sad_32x16[7] < p_best_sad32x16[7]) {
            p_best_sad32x16[7] = sad_32x16[7];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv32x16[7] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        // 64x16
        sad = sad_32x16[0] + sad_32x16[2];
        if (sad < p_best_sad64x16[0]) {
            p_best_sad64x16[0] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv64x16[0] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }
        sad = sad_32x16[1] + sad_32x16[3];
        if (sad < p_best_sad64x16[1]) {
            p_best_sad64x16[1] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv64x16[1] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_32x16[4] + sad_32x16[6];
        if (sad < p_best_sad64x16[2]) {
            p_best_sad64x16[2] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv64x16[2] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }
        sad = sad_32x16[5] + sad_32x16[7];
        if (sad < p_best_sad64x16[3]) {
            p_best_sad64x16[3] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv64x16[3] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }
        // 16x8
        sad_16x8[0] = p_sad8x8[0][search_index] + p_sad8x8[1][search_index];
        if (sad_16x8[0] < p_best_sad16x8[0]) {
            p_best_sad16x8[0] = sad_16x8[0];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x8[0] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[1] = p_sad8x8[2][search_index] + p_sad8x8[3][search_index];
        if (sad_16x8[1] < p_best_sad16x8[1]) {
            p_best_sad16x8[1] = sad_16x8[1];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x8[1] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[2] = p_sad8x8[4][search_index] + p_sad8x8[5][search_index];
        if (sad_16x8[2] < p_best_sad16x8[2]) {
            p_best_sad16x8[2] = sad_16x8[2];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x8[2] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[3] = p_sad8x8[6][search_index] + p_sad8x8[7][search_index];
        if (sad_16x8[3] < p_best_sad16x8[3]) {
            p_best_sad16x8[3] = sad_16x8[3];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x8[3] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[4] = p_sad8x8[8][search_index] + p_sad8x8[9][search_index];
        if (sad_16x8[4] < p_best_sad16x8[4]) {
            p_best_sad16x8[4] = sad_16x8[4];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x8[4] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[5] = p_sad8x8[10][search_index] + p_sad8x8[11][search_index];
        if (sad_16x8[5] < p_best_sad16x8[5]) {
            p_best_sad16x8[5] = sad_16x8[5];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x8[5] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[6] = p_sad8x8[12][search_index] + p_sad8x8[13][search_index];
        if (sad_16x8[6] < p_best_sad16x8[6]) {
            p_best_sad16x8[6] = sad_16x8[6];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x8[6] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[7] = p_sad8x8[14][search_index] + p_sad8x8[15][search_index];
        if (sad_16x8[7] < p_best_sad16x8[7]) {
            p_best_sad16x8[7] = sad_16x8[7];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x8[7] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[8] = p_sad8x8[16][search_index] + p_sad8x8[17][search_index];
        if (sad_16x8[8] < p_best_sad16x8[8]) {
            p_best_sad16x8[8] = sad_16x8[8];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x8[8] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[9] = p_sad8x8[18][search_index] + p_sad8x8[19][search_index];
        if (sad_16x8[9] < p_best_sad16x8[9]) {
            p_best_sad16x8[9] = sad_16x8[9];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x8[9] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[10] = p_sad8x8[20][search_index] + p_sad8x8[21][search_index];
        if (sad_16x8[10] < p_best_sad16x8[10]) {
            p_best_sad16x8[10] = sad_16x8[10];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x8[10] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[11] = p_sad8x8[22][search_index] + p_sad8x8[23][search_index];
        if (sad_16x8[11] < p_best_sad16x8[11]) {
            p_best_sad16x8[11] = sad_16x8[11];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x8[11] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[12] = p_sad8x8[24][search_index] + p_sad8x8[25][search_index];
        if (sad_16x8[12] < p_best_sad16x8[12]) {
            p_best_sad16x8[12] = sad_16x8[12];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x8[12] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[13] = p_sad8x8[26][search_index] + p_sad8x8[27][search_index];
        if (sad_16x8[13] < p_best_sad16x8[13]) {
            p_best_sad16x8[13] = sad_16x8[13];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x8[13] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[14] = p_sad8x8[28][search_index] + p_sad8x8[29][search_index];
        if (sad_16x8[14] < p_best_sad16x8[14]) {
            p_best_sad16x8[14] = sad_16x8[14];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x8[14] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[15] = p_sad8x8[30][search_index] + p_sad8x8[31][search_index];
        if (sad_16x8[15] < p_best_sad16x8[15]) {
            p_best_sad16x8[15] = sad_16x8[15];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x8[15] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[16] = p_sad8x8[32][search_index] + p_sad8x8[33][search_index];
        if (sad_16x8[16] < p_best_sad16x8[16]) {
            p_best_sad16x8[16] = sad_16x8[16];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x8[16] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[17] = p_sad8x8[34][search_index] + p_sad8x8[35][search_index];
        if (sad_16x8[17] < p_best_sad16x8[17]) {
            p_best_sad16x8[17] = sad_16x8[17];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x8[17] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[18] = p_sad8x8[36][search_index] + p_sad8x8[37][search_index];
        if (sad_16x8[18] < p_best_sad16x8[18]) {
            p_best_sad16x8[18] = sad_16x8[18];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x8[18] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[19] = p_sad8x8[38][search_index] + p_sad8x8[39][search_index];
        if (sad_16x8[19] < p_best_sad16x8[19]) {
            p_best_sad16x8[19] = sad_16x8[19];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x8[19] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[20] = p_sad8x8[40][search_index] + p_sad8x8[41][search_index];
        if (sad_16x8[20] < p_best_sad16x8[20]) {
            p_best_sad16x8[20] = sad_16x8[20];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x8[20] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[21] = p_sad8x8[42][search_index] + p_sad8x8[43][search_index];
        if (sad_16x8[21] < p_best_sad16x8[21]) {
            p_best_sad16x8[21] = sad_16x8[21];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x8[21] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[22] = p_sad8x8[44][search_index] + p_sad8x8[45][search_index];
        if (sad_16x8[22] < p_best_sad16x8[22]) {
            p_best_sad16x8[22] = sad_16x8[22];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x8[22] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[23] = p_sad8x8[46][search_index] + p_sad8x8[47][search_index];
        if (sad_16x8[23] < p_best_sad16x8[23]) {
            p_best_sad16x8[23] = sad_16x8[23];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x8[23] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[24] = p_sad8x8[48][search_index] + p_sad8x8[49][search_index];
        if (sad_16x8[24] < p_best_sad16x8[24]) {
            p_best_sad16x8[24] = sad_16x8[24];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x8[24] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[25] = p_sad8x8[50][search_index] + p_sad8x8[51][search_index];
        if (sad_16x8[25] < p_best_sad16x8[25]) {
            p_best_sad16x8[25] = sad_16x8[25];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x8[25] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[26] = p_sad8x8[52][search_index] + p_sad8x8[53][search_index];
        if (sad_16x8[26] < p_best_sad16x8[26]) {
            p_best_sad16x8[26] = sad_16x8[26];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x8[26] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[27] = p_sad8x8[54][search_index] + p_sad8x8[55][search_index];
        if (sad_16x8[27] < p_best_sad16x8[27]) {
            p_best_sad16x8[27] = sad_16x8[27];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x8[27] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[28] = p_sad8x8[56][search_index] + p_sad8x8[57][search_index];
        if (sad_16x8[28] < p_best_sad16x8[28]) {
            p_best_sad16x8[28] = sad_16x8[28];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x8[28] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[29] = p_sad8x8[58][search_index] + p_sad8x8[59][search_index];
        if (sad_16x8[29] < p_best_sad16x8[29]) {
            p_best_sad16x8[29] = sad_16x8[29];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x8[29] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[30] = p_sad8x8[60][search_index] + p_sad8x8[61][search_index];
        if (sad_16x8[30] < p_best_sad16x8[30]) {
            p_best_sad16x8[30] = sad_16x8[30];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x8[30] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x8[31] = p_sad8x8[62][search_index] + p_sad8x8[63][search_index];
        if (sad_16x8[31] < p_best_sad16x8[31]) {
            p_best_sad16x8[31] = sad_16x8[31];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x8[31] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        // 32x64
        sad = p_sad32x32[0][search_index] + p_sad32x32[2][search_index];
        if (sad < p_best_sad32x64[0]) {
            p_best_sad32x64[0] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv32x64[0] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = p_sad32x32[1][search_index] + p_sad32x32[3][search_index];
        if (sad < p_best_sad32x64[1]) {
            p_best_sad32x64[1] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv32x64[1] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        // 16x32
        sad_16x32[0] =
            p_sad16x16[0][search_index] + p_sad16x16[2][search_index];
        if (sad_16x32[0] < p_best_sad16x32[0]) {
            p_best_sad16x32[0] = sad_16x32[0];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x32[0] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x32[1] =
            p_sad16x16[1][search_index] + p_sad16x16[3][search_index];
        if (sad_16x32[1] < p_best_sad16x32[1]) {
            p_best_sad16x32[1] = sad_16x32[1];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x32[1] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x32[2] =
            p_sad16x16[4][search_index] + p_sad16x16[6][search_index];
        if (sad_16x32[2] < p_best_sad16x32[2]) {
            p_best_sad16x32[2] = sad_16x32[2];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x32[2] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x32[3] =
            p_sad16x16[5][search_index] + p_sad16x16[7][search_index];
        if (sad_16x32[3] < p_best_sad16x32[3]) {
            p_best_sad16x32[3] = sad_16x32[3];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x32[3] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x32[4] =
            p_sad16x16[8][search_index] + p_sad16x16[10][search_index];
        if (sad_16x32[4] < p_best_sad16x32[4]) {
            p_best_sad16x32[4] = sad_16x32[4];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x32[4] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x32[5] =
            p_sad16x16[9][search_index] + p_sad16x16[11][search_index];
        if (sad_16x32[5] < p_best_sad16x32[5]) {
            p_best_sad16x32[5] = sad_16x32[5];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x32[5] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x32[6] =
            p_sad16x16[12][search_index] + p_sad16x16[14][search_index];
        if (sad_16x32[6] < p_best_sad16x32[6]) {
            p_best_sad16x32[6] = sad_16x32[6];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x32[6] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_16x32[7] =
            p_sad16x16[13][search_index] + p_sad16x16[15][search_index];
        if (sad_16x32[7] < p_best_sad16x32[7]) {
            p_best_sad16x32[7] = sad_16x32[7];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x32[7] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x32[0] + sad_16x32[4];
        if (sad < p_best_sad16x64[0]) {
            p_best_sad16x64[0] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x64[0] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }
        sad = sad_16x32[1] + sad_16x32[5];
        if (sad < p_best_sad16x64[1]) {
            p_best_sad16x64[1] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x64[1] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x32[2] + sad_16x32[6];
        if (sad < p_best_sad16x64[2]) {
            p_best_sad16x64[2] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x64[2] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x32[3] + sad_16x32[7];
        if (sad < p_best_sad16x64[3]) {
            p_best_sad16x64[3] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x64[3] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }
        // 8x16
        sad_8x16[0] = p_sad8x8[0][search_index] + p_sad8x8[2][search_index];
        if (sad_8x16[0] < p_best_sad8x16[0]) {
            p_best_sad8x16[0] = sad_8x16[0];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x16[0] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[1] = p_sad8x8[1][search_index] + p_sad8x8[3][search_index];
        if (sad_8x16[1] < p_best_sad8x16[1]) {
            p_best_sad8x16[1] = sad_8x16[1];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x16[1] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[2] = p_sad8x8[4][search_index] + p_sad8x8[6][search_index];
        if (sad_8x16[2] < p_best_sad8x16[2]) {
            p_best_sad8x16[2] = sad_8x16[2];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x16[2] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[3] = p_sad8x8[5][search_index] + p_sad8x8[7][search_index];
        if (sad_8x16[3] < p_best_sad8x16[3]) {
            p_best_sad8x16[3] = sad_8x16[3];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x16[3] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[4] = p_sad8x8[8][search_index] + p_sad8x8[10][search_index];
        if (sad_8x16[4] < p_best_sad8x16[4]) {
            p_best_sad8x16[4] = sad_8x16[4];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x16[4] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[5] = p_sad8x8[9][search_index] + p_sad8x8[11][search_index];
        if (sad_8x16[5] < p_best_sad8x16[5]) {
            p_best_sad8x16[5] = sad_8x16[5];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x16[5] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[6] = p_sad8x8[12][search_index] + p_sad8x8[14][search_index];
        if (sad_8x16[6] < p_best_sad8x16[6]) {
            p_best_sad8x16[6] = sad_8x16[6];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x16[6] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[7] = p_sad8x8[13][search_index] + p_sad8x8[15][search_index];
        if (sad_8x16[7] < p_best_sad8x16[7]) {
            p_best_sad8x16[7] = sad_8x16[7];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x16[7] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[8] = p_sad8x8[16][search_index] + p_sad8x8[18][search_index];
        if (sad_8x16[8] < p_best_sad8x16[8]) {
            p_best_sad8x16[8] = sad_8x16[8];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x16[8] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[9] = p_sad8x8[17][search_index] + p_sad8x8[19][search_index];
        if (sad_8x16[9] < p_best_sad8x16[9]) {
            p_best_sad8x16[9] = sad_8x16[9];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x16[9] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[10] = p_sad8x8[20][search_index] + p_sad8x8[22][search_index];
        if (sad_8x16[10] < p_best_sad8x16[10]) {
            p_best_sad8x16[10] = sad_8x16[10];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x16[10] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[11] = p_sad8x8[21][search_index] + p_sad8x8[23][search_index];
        if (sad_8x16[11] < p_best_sad8x16[11]) {
            p_best_sad8x16[11] = sad_8x16[11];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x16[11] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[12] = p_sad8x8[24][search_index] + p_sad8x8[26][search_index];
        if (sad_8x16[12] < p_best_sad8x16[12]) {
            p_best_sad8x16[12] = sad_8x16[12];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x16[12] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[13] = p_sad8x8[25][search_index] + p_sad8x8[27][search_index];
        if (sad_8x16[13] < p_best_sad8x16[13]) {
            p_best_sad8x16[13] = sad_8x16[13];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x16[13] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[14] = p_sad8x8[28][search_index] + p_sad8x8[30][search_index];
        if (sad_8x16[14] < p_best_sad8x16[14]) {
            p_best_sad8x16[14] = sad_8x16[14];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x16[14] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[15] = p_sad8x8[29][search_index] + p_sad8x8[31][search_index];
        if (sad_8x16[15] < p_best_sad8x16[15]) {
            p_best_sad8x16[15] = sad_8x16[15];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x16[15] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[16] = p_sad8x8[32][search_index] + p_sad8x8[34][search_index];
        if (sad_8x16[16] < p_best_sad8x16[16]) {
            p_best_sad8x16[16] = sad_8x16[16];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x16[16] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[17] = p_sad8x8[33][search_index] + p_sad8x8[35][search_index];
        if (sad_8x16[17] < p_best_sad8x16[17]) {
            p_best_sad8x16[17] = sad_8x16[17];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x16[17] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[18] = p_sad8x8[36][search_index] + p_sad8x8[38][search_index];
        if (sad_8x16[18] < p_best_sad8x16[18]) {
            p_best_sad8x16[18] = sad_8x16[18];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x16[18] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[19] = p_sad8x8[37][search_index] + p_sad8x8[39][search_index];
        if (sad_8x16[19] < p_best_sad8x16[19]) {
            p_best_sad8x16[19] = sad_8x16[19];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x16[19] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[20] = p_sad8x8[40][search_index] + p_sad8x8[42][search_index];
        if (sad_8x16[20] < p_best_sad8x16[20]) {
            p_best_sad8x16[20] = sad_8x16[20];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x16[20] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[21] = p_sad8x8[41][search_index] + p_sad8x8[43][search_index];
        if (sad_8x16[21] < p_best_sad8x16[21]) {
            p_best_sad8x16[21] = sad_8x16[21];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x16[21] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[22] = p_sad8x8[44][search_index] + p_sad8x8[46][search_index];
        if (sad_8x16[22] < p_best_sad8x16[22]) {
            p_best_sad8x16[22] = sad_8x16[22];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x16[22] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[23] = p_sad8x8[45][search_index] + p_sad8x8[47][search_index];
        if (sad_8x16[23] < p_best_sad8x16[23]) {
            p_best_sad8x16[23] = sad_8x16[23];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x16[23] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[24] = p_sad8x8[48][search_index] + p_sad8x8[50][search_index];
        if (sad_8x16[24] < p_best_sad8x16[24]) {
            p_best_sad8x16[24] = sad_8x16[24];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x16[24] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[25] = p_sad8x8[49][search_index] + p_sad8x8[51][search_index];
        if (sad_8x16[25] < p_best_sad8x16[25]) {
            p_best_sad8x16[25] = sad_8x16[25];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x16[25] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[26] = p_sad8x8[52][search_index] + p_sad8x8[54][search_index];
        if (sad_8x16[26] < p_best_sad8x16[26]) {
            p_best_sad8x16[26] = sad_8x16[26];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x16[26] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[27] = p_sad8x8[53][search_index] + p_sad8x8[55][search_index];
        if (sad_8x16[27] < p_best_sad8x16[27]) {
            p_best_sad8x16[27] = sad_8x16[27];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x16[27] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[28] = p_sad8x8[56][search_index] + p_sad8x8[58][search_index];
        if (sad_8x16[28] < p_best_sad8x16[28]) {
            p_best_sad8x16[28] = sad_8x16[28];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x16[28] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[29] = p_sad8x8[57][search_index] + p_sad8x8[59][search_index];
        if (sad_8x16[29] < p_best_sad8x16[29]) {
            p_best_sad8x16[29] = sad_8x16[29];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x16[29] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[30] = p_sad8x8[60][search_index] + p_sad8x8[62][search_index];
        if (sad_8x16[30] < p_best_sad8x16[30]) {
            p_best_sad8x16[30] = sad_8x16[30];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x16[30] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad_8x16[31] = p_sad8x8[61][search_index] + p_sad8x8[63][search_index];
        if (sad_8x16[31] < p_best_sad8x16[31]) {
            p_best_sad8x16[31] = sad_8x16[31];
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x16[31] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        // 32x8
        sad = sad_16x8[0] + sad_16x8[2];
        if (sad < p_best_sad32x8[0]) {
            p_best_sad32x8[0] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv32x8[0] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x8[1] + sad_16x8[3];
        if (sad < p_best_sad32x8[1]) {
            p_best_sad32x8[1] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv32x8[1] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x8[4] + sad_16x8[6];
        if (sad < p_best_sad32x8[2]) {
            p_best_sad32x8[2] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv32x8[2] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x8[5] + sad_16x8[7];
        if (sad < p_best_sad32x8[3]) {
            p_best_sad32x8[3] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv32x8[3] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x8[8] + sad_16x8[10];
        if (sad < p_best_sad32x8[4]) {
            p_best_sad32x8[4] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv32x8[4] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x8[9] + sad_16x8[11];
        if (sad < p_best_sad32x8[5]) {
            p_best_sad32x8[5] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv32x8[5] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x8[12] + sad_16x8[14];
        if (sad < p_best_sad32x8[6]) {
            p_best_sad32x8[6] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv32x8[6] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x8[13] + sad_16x8[15];
        if (sad < p_best_sad32x8[7]) {
            p_best_sad32x8[7] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv32x8[7] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x8[16] + sad_16x8[18];
        if (sad < p_best_sad32x8[8]) {
            p_best_sad32x8[8] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv32x8[8] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x8[17] + sad_16x8[19];
        if (sad < p_best_sad32x8[9]) {
            p_best_sad32x8[9] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv32x8[9] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x8[20] + sad_16x8[22];
        if (sad < p_best_sad32x8[10]) {
            p_best_sad32x8[10] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv32x8[10] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x8[21] + sad_16x8[23];
        if (sad < p_best_sad32x8[11]) {
            p_best_sad32x8[11] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv32x8[11] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x8[24] + sad_16x8[26];
        if (sad < p_best_sad32x8[12]) {
            p_best_sad32x8[12] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv32x8[12] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x8[25] + sad_16x8[27];
        if (sad < p_best_sad32x8[13]) {
            p_best_sad32x8[13] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv32x8[13] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x8[28] + sad_16x8[30];
        if (sad < p_best_sad32x8[14]) {
            p_best_sad32x8[14] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv32x8[14] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_16x8[29] + sad_16x8[31];
        if (sad < p_best_sad32x8[15]) {
            p_best_sad32x8[15] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv32x8[15] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }
        // 8x32
        sad = sad_8x16[0] + sad_8x16[4];
        if (sad < p_best_sad8x32[0]) {
            p_best_sad8x32[0] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x32[0] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_8x16[1] + sad_8x16[5];
        if (sad < p_best_sad8x32[1]) {
            p_best_sad8x32[1] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x32[1] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_8x16[2] + sad_8x16[6];
        if (sad < p_best_sad8x32[2]) {
            p_best_sad8x32[2] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x32[2] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_8x16[3] + sad_8x16[7];
        if (sad < p_best_sad8x32[3]) {
            p_best_sad8x32[3] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x32[3] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_8x16[8] + sad_8x16[12];
        if (sad < p_best_sad8x32[4]) {
            p_best_sad8x32[4] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x32[4] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_8x16[9] + sad_8x16[13];
        if (sad < p_best_sad8x32[5]) {
            p_best_sad8x32[5] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x32[5] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_8x16[10] + sad_8x16[14];
        if (sad < p_best_sad8x32[6]) {
            p_best_sad8x32[6] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x32[6] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_8x16[11] + sad_8x16[15];
        if (sad < p_best_sad8x32[7]) {
            p_best_sad8x32[7] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x32[7] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_8x16[16] + sad_8x16[20];
        if (sad < p_best_sad8x32[8]) {
            p_best_sad8x32[8] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x32[8] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_8x16[17] + sad_8x16[21];
        if (sad < p_best_sad8x32[9]) {
            p_best_sad8x32[9] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x32[9] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_8x16[18] + sad_8x16[22];
        if (sad < p_best_sad8x32[10]) {
            p_best_sad8x32[10] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x32[10] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_8x16[19] + sad_8x16[23];
        if (sad < p_best_sad8x32[11]) {
            p_best_sad8x32[11] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x32[11] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_8x16[24] + sad_8x16[28];
        if (sad < p_best_sad8x32[12]) {
            p_best_sad8x32[12] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x32[12] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_8x16[25] + sad_8x16[29];
        if (sad < p_best_sad8x32[13]) {
            p_best_sad8x32[13] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x32[13] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_8x16[26] + sad_8x16[30];
        if (sad < p_best_sad8x32[14]) {
            p_best_sad8x32[14] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x32[14] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad = sad_8x16[27] + sad_8x16[31];
        if (sad < p_best_sad8x32[15]) {
            p_best_sad8x32[15] = sad;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x32[15] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }
    }
}

/*******************************************
 * ext_eight_sad_calculation_8x8_16x16
 *******************************************/
static void ext_eight_sad_calculation_8x8_16x16(
    uint8_t *src, uint32_t src_stride, uint8_t *ref, uint32_t ref_stride,
    uint32_t mv, uint32_t start_16x16_pos, uint32_t *p_best_sad8x8,
    uint32_t *p_best_sad16x16, uint32_t *p_best_mv8x8, uint32_t *p_best_mv16x16,
    uint32_t p_eight_sad16x16[16][8], uint32_t p_eight_sad8x8[64][8]) {
    const uint32_t start_8x8_pos = 4 * start_16x16_pos;
    uint32_t sad8x8_0, sad8x8_1, sad8x8_2, sad8x8_3;
    uint32_t sad16x16;
    uint32_t search_index;
    int16_t x_mv, y_mv;
    uint32_t srcStrideSub = (src_stride << 1);
    uint32_t refStrideSub = (ref_stride << 1);

    p_best_sad8x8 += start_8x8_pos;
    p_best_mv8x8 += start_8x8_pos;
    p_best_sad16x16 += start_16x16_pos;
    p_best_mv16x16 += start_16x16_pos;

    for (search_index = 0; search_index < 8; search_index++) {
        p_eight_sad8x8[0 + start_8x8_pos][search_index] = sad8x8_0 =
            (compute8x4_sad_kernel_c(
                src, srcStrideSub, ref + search_index, refStrideSub))
            << 1;
        if (sad8x8_0 < p_best_sad8x8[0]) {
            p_best_sad8x8[0] = (uint32_t)sad8x8_0;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x8[0] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        p_eight_sad8x8[1 + start_8x8_pos][search_index] = sad8x8_1 =
            (compute8x4_sad_kernel_c(
                src + 8, srcStrideSub, ref + 8 + search_index, refStrideSub))
            << 1;
        if (sad8x8_1 < p_best_sad8x8[1]) {
            p_best_sad8x8[1] = (uint32_t)sad8x8_1;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x8[1] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        p_eight_sad8x8[2 + start_8x8_pos][search_index] = sad8x8_2 =
            (compute8x4_sad_kernel_c(
                src + (src_stride << 3),
                srcStrideSub,
                ref + (ref_stride << 3) + search_index,
                refStrideSub))
            << 1;
        if (sad8x8_2 < p_best_sad8x8[2]) {
            p_best_sad8x8[2] = (uint32_t)sad8x8_2;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x8[2] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        p_eight_sad8x8[3 + start_8x8_pos][search_index] = sad8x8_3 =
            (compute8x4_sad_kernel_c(
                src + (src_stride << 3) + 8,
                srcStrideSub,
                ref + (ref_stride << 3) + 8 + search_index,
                refStrideSub))
            << 1;
        if (sad8x8_3 < p_best_sad8x8[3]) {
            p_best_sad8x8[3] = (uint32_t)sad8x8_3;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv8x8[3] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        p_eight_sad16x16[start_16x16_pos][search_index] = sad16x16 =
            sad8x8_0 + sad8x8_1 + sad8x8_2 + sad8x8_3;
        if (sad16x16 < p_best_sad16x16[0]) {
            p_best_sad16x16[0] = (uint32_t)sad16x16;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv16x16[0] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }
    }
}

void ext_all_sad_calculation_8x8_16x16_c(
    uint8_t *src, uint32_t src_stride, uint8_t *ref, uint32_t ref_stride,
    uint32_t mv, uint32_t *p_best_sad8x8, uint32_t *p_best_sad16x16,
    uint32_t *p_best_mv8x8, uint32_t *p_best_mv16x16,
    uint32_t p_eight_sad16x16[16][8], uint32_t p_eight_sad8x8[64][8]) {
    static const char offsets[16] = {
        0, 1, 4, 5, 2, 3, 6, 7, 8, 9, 12, 13, 10, 11, 14, 15};

    //---- 16x16 : 0, 1, 4, 5, 2, 3, 6, 7, 8, 9, 12, 13, 10, 11, 14, 15
    for (int y = 0; y < 4; y++) {
        for (int x = 0; x < 4; x++) {
            const uint32_t blockIndex = 16 * y * src_stride + 16 * x;
            const uint32_t searchPositionIndex = 16 * y * ref_stride + 16 * x;
            ext_eight_sad_calculation_8x8_16x16(src + blockIndex,
                                                src_stride,
                                                ref + searchPositionIndex,
                                                ref_stride,
                                                mv,
                                                offsets[4 * y + x],
                                                p_best_sad8x8,
                                                p_best_sad16x16,
                                                p_best_mv8x8,
                                                p_best_mv16x16,
                                                p_eight_sad16x16,
                                                p_eight_sad8x8);
        }
    }
}

/*******************************************
Calcualte SAD for 32x32,64x64 from 16x16
and check if there is improvment, if yes keep
the best SAD+MV
*******************************************/
void ext_eight_sad_calculation_32x32_64x64_c(
    uint32_t p_sad16x16[16][8], uint32_t *p_best_sad32x32,
    uint32_t *p_best_sad64x64, uint32_t *p_best_mv32x32,
    uint32_t *p_best_mv64x64, uint32_t mv, uint32_t p_sad32x32[4][8]) {
    uint32_t search_index;
    int16_t x_mv, y_mv;
    for (search_index = 0; search_index < 8; search_index++) {
        uint32_t sad32x32_0, sad32x32_1, sad32x32_2, sad32x32_3, sad64x64;

        p_sad32x32[0][search_index] = sad32x32_0 =
            p_sad16x16[0][search_index] + p_sad16x16[1][search_index] +
            p_sad16x16[2][search_index] + p_sad16x16[3][search_index];
        if (sad32x32_0 < p_best_sad32x32[0]) {
            p_best_sad32x32[0] = sad32x32_0;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv32x32[0] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        p_sad32x32[1][search_index] = sad32x32_1 =
            p_sad16x16[4][search_index] + p_sad16x16[5][search_index] +
            p_sad16x16[6][search_index] + p_sad16x16[7][search_index];
        if (sad32x32_1 < p_best_sad32x32[1]) {
            p_best_sad32x32[1] = sad32x32_1;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv32x32[1] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        p_sad32x32[2][search_index] = sad32x32_2 =
            p_sad16x16[8][search_index] + p_sad16x16[9][search_index] +
            p_sad16x16[10][search_index] + p_sad16x16[11][search_index];
        if (sad32x32_2 < p_best_sad32x32[2]) {
            p_best_sad32x32[2] = sad32x32_2;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv32x32[2] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        p_sad32x32[3][search_index] = sad32x32_3 =
            p_sad16x16[12][search_index] + p_sad16x16[13][search_index] +
            p_sad16x16[14][search_index] + p_sad16x16[15][search_index];
        if (sad32x32_3 < p_best_sad32x32[3]) {
            p_best_sad32x32[3] = sad32x32_3;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv32x32[3] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }

        sad64x64 = sad32x32_0 + sad32x32_1 + sad32x32_2 + sad32x32_3;
        if (sad64x64 < p_best_sad64x64[0]) {
            p_best_sad64x64[0] = sad64x64;
            x_mv = _MVXT(mv) + (int16_t)search_index * 4;
            y_mv = _MVYT(mv);
            p_best_mv64x64[0] = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
        }
    }
}

/*******************************************
 * open_loop_me_get_search_point_results_block
 *******************************************/
static void open_loop_me_get_eight_search_point_results_block(
    MeContext
        *context_ptr,    // input parameter, ME context Ptr, used to get SB Ptr
    uint32_t listIndex,  // input parameter, reference list index
    uint32_t ref_pic_index,
    uint32_t searchRegionIndex,  // input parameter, search area origin, used to
                                 // point to reference samples
    int32_t xSearchIndex,  // input parameter, search region position in the
                           // horizontal direction, used to derive xMV
    int32_t ySearchIndex  // input parameter, search region position in the
                           // vertical direction, used to derive yMV
    ) {
    // uint32_t reflumaStride = refPicPtr->stride_y; // NADER
    // uint8_t  *refPtr = refPicPtr->buffer_y; // NADER
    uint32_t reflumaStride =
        context_ptr->interpolated_full_stride[listIndex][ref_pic_index];
    uint8_t *refPtr =
        context_ptr->integer_buffer_ptr[listIndex][ref_pic_index] +
        ((ME_FILTER_TAP >> 1) *
         context_ptr->interpolated_full_stride[listIndex][ref_pic_index]) +
        (ME_FILTER_TAP >> 1) + searchRegionIndex;

    uint32_t currMV1 = (((uint16_t)ySearchIndex) << 18);
    uint16_t currMV2 = (((uint16_t)xSearchIndex << 2));
    uint32_t currMV = currMV1 | currMV2;

    ext_all_sad_calculation_8x8_16x16(
        context_ptr->sb_src_ptr,
        context_ptr->sb_src_stride,
        refPtr,
        reflumaStride,
        currMV,
        context_ptr->p_best_sad8x8,
        context_ptr->p_best_sad16x16,
        context_ptr->p_best_mv8x8,
        context_ptr->p_best_mv16x16,
        context_ptr->p_eight_sad16x16,
        context_ptr->p_eight_sad8x8);

    ext_eight_sad_calculation_32x32_64x64(
        context_ptr->p_eight_sad16x16,
        context_ptr->p_best_sad32x32,
        context_ptr->p_best_sad64x64,
        context_ptr->p_best_mv32x32,
        context_ptr->p_best_mv64x64,
        currMV,
        context_ptr->p_eight_sad32x32);

    ext_eigth_sad_calculation_nsq(
        context_ptr->p_eight_sad8x8,
        context_ptr->p_eight_sad16x16,
        context_ptr->p_eight_sad32x32,
        context_ptr->p_best_sad64x32,
        context_ptr->p_best_mv64x32,
        context_ptr->p_best_sad32x16,
        context_ptr->p_best_mv32x16,
        context_ptr->p_best_sad16x8,
        context_ptr->p_best_mv16x8,
        context_ptr->p_best_sad32x64,
        context_ptr->p_best_mv32x64,
        context_ptr->p_best_sad16x32,
        context_ptr->p_best_mv16x32,
        context_ptr->p_best_sad8x16,
        context_ptr->p_best_mv8x16,
        context_ptr->p_best_sad32x8,
        context_ptr->p_best_mv32x8,
        context_ptr->p_best_sad8x32,
        context_ptr->p_best_mv8x32,
        context_ptr->p_best_sad64x16,
        context_ptr->p_best_mv64x16,
        context_ptr->p_best_sad16x64,
        context_ptr->p_best_mv16x64,
        currMV);
}

/*******************************************
 * nsq_get_analysis_results_block returns the
 * the best partition for each sq_block based
 * on the ME SAD
 *******************************************/
static void nsq_get_analysis_results_block(MeContext *context_ptr) {
    uint32_t *p_best_sad64x32 = context_ptr->p_best_sad64x32;
    uint32_t *p_best_sad32x16 = context_ptr->p_best_sad32x16;
    uint32_t *p_best_sad16x8 = context_ptr->p_best_sad16x8;
    uint32_t *p_best_sad32x64 = context_ptr->p_best_sad32x64;
    uint32_t *p_best_sad16x32 = context_ptr->p_best_sad16x32;
    uint32_t *p_best_sad8x16 = context_ptr->p_best_sad8x16;
    uint32_t *p_best_sad32x8 = context_ptr->p_best_sad32x8;
    uint32_t *p_best_sad8x32 = context_ptr->p_best_sad8x32;
    uint32_t *p_best_sad64x16 = context_ptr->p_best_sad64x16;
    uint32_t *p_best_sad16x64 = context_ptr->p_best_sad16x64;
    uint8_t *p_best_nsq_64x64 = context_ptr->p_best_nsq64x64;
    uint8_t *p_best_nsq_32x32 = context_ptr->p_best_nsq32x32;
    uint8_t *p_best_nsq_16x16 = context_ptr->p_best_nsq16x16;
    uint8_t *p_best_nsq_8x8 = context_ptr->p_best_nsq8x8;

    nsq_me_analysis(p_best_sad64x32,
                    p_best_sad32x16,
                    p_best_sad16x8,
                    p_best_sad32x64,
                    p_best_sad16x32,
                    p_best_sad8x16,
                    p_best_sad32x8,
                    p_best_sad8x32,
                    p_best_sad64x16,
                    p_best_sad16x64,
                    p_best_nsq_64x64,
                    p_best_nsq_32x32,
                    p_best_nsq_16x16,
                    p_best_nsq_8x8);
}

/*******************************************
 * open_loop_me_get_search_point_results_block
 *******************************************/
static void open_loop_me_get_search_point_results_block(
    MeContext
        *context_ptr,    // input parameter, ME context Ptr, used to get SB Ptr
    uint32_t listIndex,  // input parameter, reference list index
    uint32_t ref_pic_index,
    uint32_t searchRegionIndex,  // input parameter, search area origin, used to
                                 // point to reference samples
    int32_t xSearchIndex,  // input parameter, search region position in the
                           // horizontal direction, used to derive xMV
    int32_t ySearchIndex,  // input parameter, search region position in the
                           // vertical direction, used to derive yMV
    EbAsm asm_type) {
    const EbBool sub_sad = (context_ptr->me_search_method == SUB_SAD_SEARCH);
    uint8_t *src_ptr = context_ptr->sb_src_ptr;

    // uint8_t  *refPtr = refPicPtr->buffer_y; // NADER
    uint8_t *refPtr =
        context_ptr->integer_buffer_ptr[listIndex][ref_pic_index] +
        (ME_FILTER_TAP >> 1) +
        ((ME_FILTER_TAP >> 1) *
         context_ptr->interpolated_full_stride[listIndex][ref_pic_index]);
    // uint32_t reflumaStride = refPicPtr->stride_y; // NADER
    uint32_t reflumaStride =
        context_ptr->interpolated_full_stride[listIndex][ref_pic_index];
    uint32_t searchPositionTLIndex = searchRegionIndex;
    uint32_t searchPositionIndex;
    uint32_t blockIndex;
    uint32_t srcNext16x16Offset = (BLOCK_SIZE_64 << 4);
    // uint32_t refNext16x16Offset = (refPicPtr->stride_y << 4); // NADER
    uint32_t refNext16x16Offset = (reflumaStride << 4);
    uint32_t currMV1 = (((uint16_t)ySearchIndex) << 18);
    uint16_t currMV2 = (((uint16_t)xSearchIndex << 2));
    uint32_t currMV = currMV1 | currMV2;
    uint32_t *p_best_sad8x8 = context_ptr->p_best_sad8x8;
    uint32_t *p_best_sad16x16 = context_ptr->p_best_sad16x16;
    uint32_t *p_best_sad32x32 = context_ptr->p_best_sad32x32;
    uint32_t *p_best_sad64x64 = context_ptr->p_best_sad64x64;
    uint32_t *p_best_sad64x32 = context_ptr->p_best_sad64x32;
    uint32_t *p_best_sad32x16 = context_ptr->p_best_sad32x16;
    uint32_t *p_best_sad16x8 = context_ptr->p_best_sad16x8;
    uint32_t *p_best_sad32x64 = context_ptr->p_best_sad32x64;
    uint32_t *p_best_sad16x32 = context_ptr->p_best_sad16x32;
    uint32_t *p_best_sad8x16 = context_ptr->p_best_sad8x16;
    uint32_t *p_best_sad32x8 = context_ptr->p_best_sad32x8;
    uint32_t *p_best_sad8x32 = context_ptr->p_best_sad8x32;
    uint32_t *p_best_sad64x16 = context_ptr->p_best_sad64x16;
    uint32_t *p_best_sad16x64 = context_ptr->p_best_sad16x64;
    uint32_t *p_best_mv8x8 = context_ptr->p_best_mv8x8;
    uint32_t *p_best_mv16x16 = context_ptr->p_best_mv16x16;
    uint32_t *p_best_mv32x32 = context_ptr->p_best_mv32x32;
    uint32_t *p_best_mv64x64 = context_ptr->p_best_mv64x64;
    uint32_t *p_best_mv64x32 = context_ptr->p_best_mv64x32;
    uint32_t *p_best_mv32x16 = context_ptr->p_best_mv32x16;
    uint32_t *p_best_mv16x8 = context_ptr->p_best_mv16x8;
    uint32_t *p_best_mv32x64 = context_ptr->p_best_mv32x64;
    uint32_t *p_best_mv16x32 = context_ptr->p_best_mv16x32;
    uint32_t *p_best_mv8x16 = context_ptr->p_best_mv8x16;
    uint32_t *p_best_mv32x8 = context_ptr->p_best_mv32x8;
    uint32_t *p_best_mv8x32 = context_ptr->p_best_mv8x32;
    uint32_t *p_sad32x32 = context_ptr->p_sad32x32;
    uint32_t *p_sad16x16 = context_ptr->p_sad16x16;
    uint32_t *p_sad8x8 = context_ptr->p_sad8x8;
    uint32_t *p_best_mv64x16 = context_ptr->p_best_mv64x16;
    uint32_t *p_best_mv16x64 = context_ptr->p_best_mv16x64;

    // TODO: blockIndex searchPositionIndex could be removed  + Connect asm_type
    (void)asm_type;

    const uint32_t src_stride = context_ptr->sb_src_stride;
    srcNext16x16Offset = src_stride << 4;

    //---- 16x16 : 0
    blockIndex = 0;
    searchPositionIndex = searchPositionTLIndex;

    ext_sad_calculation_8x8_16x16(
        src_ptr + blockIndex,
        src_stride,
        refPtr + searchPositionIndex,
        reflumaStride,
        &p_best_sad8x8[0],
        &p_best_sad16x16[0],
        &p_best_mv8x8[0],
        &p_best_mv16x16[0],
        currMV,
        &p_sad16x16[0],
        &p_sad8x8[0],
        sub_sad);

    //---- 16x16 : 1
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionTLIndex + 16;
    ext_sad_calculation_8x8_16x16(
        src_ptr + blockIndex,
        src_stride,
        refPtr + searchPositionIndex,
        reflumaStride,
        &p_best_sad8x8[4],
        &p_best_sad16x16[1],
        &p_best_mv8x8[4],
        &p_best_mv16x16[1],
        currMV,
        &p_sad16x16[1],
        &p_sad8x8[4],
        sub_sad);
    //---- 16x16 : 4
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;

    ext_sad_calculation_8x8_16x16(
        src_ptr + blockIndex,
        src_stride,
        refPtr + searchPositionIndex,
        reflumaStride,
        &p_best_sad8x8[16],
        &p_best_sad16x16[4],
        &p_best_mv8x8[16],
        &p_best_mv16x16[4],
        currMV,
        &p_sad16x16[4],
        &p_sad8x8[16],
        sub_sad);

    //---- 16x16 : 5
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    ext_sad_calculation_8x8_16x16(
        src_ptr + blockIndex,
        src_stride,
        refPtr + searchPositionIndex,
        reflumaStride,
        &p_best_sad8x8[20],
        &p_best_sad16x16[5],
        &p_best_mv8x8[20],
        &p_best_mv16x16[5],
        currMV,
        &p_sad16x16[5],
        &p_sad8x8[20],
        sub_sad);

    //---- 16x16 : 2
    blockIndex = srcNext16x16Offset;
    searchPositionIndex = searchPositionTLIndex + refNext16x16Offset;
    ext_sad_calculation_8x8_16x16(
        src_ptr + blockIndex,
        src_stride,
        refPtr + searchPositionIndex,
        reflumaStride,
        &p_best_sad8x8[8],
        &p_best_sad16x16[2],
        &p_best_mv8x8[8],
        &p_best_mv16x16[2],
        currMV,
        &p_sad16x16[2],
        &p_sad8x8[8],
        sub_sad);
    //---- 16x16 : 3
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    ext_sad_calculation_8x8_16x16(
        src_ptr + blockIndex,
        src_stride,
        refPtr + searchPositionIndex,
        reflumaStride,
        &p_best_sad8x8[12],
        &p_best_sad16x16[3],
        &p_best_mv8x8[12],
        &p_best_mv16x16[3],
        currMV,
        &p_sad16x16[3],
        &p_sad8x8[12],
        sub_sad);
    //---- 16x16 : 6
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    ext_sad_calculation_8x8_16x16(
        src_ptr + blockIndex,
        src_stride,
        refPtr + searchPositionIndex,
        reflumaStride,
        &p_best_sad8x8[24],
        &p_best_sad16x16[6],
        &p_best_mv8x8[24],
        &p_best_mv16x16[6],
        currMV,
        &p_sad16x16[6],
        &p_sad8x8[24],
        sub_sad);
    //---- 16x16 : 7
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    ext_sad_calculation_8x8_16x16(
        src_ptr + blockIndex,
        src_stride,
        refPtr + searchPositionIndex,
        reflumaStride,
        &p_best_sad8x8[28],
        &p_best_sad16x16[7],
        &p_best_mv8x8[28],
        &p_best_mv16x16[7],
        currMV,
        &p_sad16x16[7],
        &p_sad8x8[28],
        sub_sad);

    //---- 16x16 : 8
    blockIndex = (srcNext16x16Offset << 1);
    searchPositionIndex = searchPositionTLIndex + (refNext16x16Offset << 1);
    ext_sad_calculation_8x8_16x16(
        src_ptr + blockIndex,
        src_stride,
        refPtr + searchPositionIndex,
        reflumaStride,
        &p_best_sad8x8[32],
        &p_best_sad16x16[8],
        &p_best_mv8x8[32],
        &p_best_mv16x16[8],
        currMV,
        &p_sad16x16[8],
        &p_sad8x8[32],
        sub_sad);
    //---- 16x16 : 9
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    ext_sad_calculation_8x8_16x16(
        src_ptr + blockIndex,
        src_stride,
        refPtr + searchPositionIndex,
        reflumaStride,
        &p_best_sad8x8[36],
        &p_best_sad16x16[9],
        &p_best_mv8x8[36],
        &p_best_mv16x16[9],
        currMV,
        &p_sad16x16[9],
        &p_sad8x8[36],
        sub_sad);
    //---- 16x16 : 12
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    ext_sad_calculation_8x8_16x16(
        src_ptr + blockIndex,
        src_stride,
        refPtr + searchPositionIndex,
        reflumaStride,
        &p_best_sad8x8[48],
        &p_best_sad16x16[12],
        &p_best_mv8x8[48],
        &p_best_mv16x16[12],
        currMV,
        &p_sad16x16[12],
        &p_sad8x8[48],
        sub_sad);
    //---- 16x16 : 13
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    ext_sad_calculation_8x8_16x16(
        src_ptr + blockIndex,
        src_stride,
        refPtr + searchPositionIndex,
        reflumaStride,
        &p_best_sad8x8[52],
        &p_best_sad16x16[13],
        &p_best_mv8x8[52],
        &p_best_mv16x16[13],
        currMV,
        &p_sad16x16[13],
        &p_sad8x8[52],
        sub_sad);

    //---- 16x16 : 10
    blockIndex = (srcNext16x16Offset * 3);
    searchPositionIndex = searchPositionTLIndex + (refNext16x16Offset * 3);
    ext_sad_calculation_8x8_16x16(
        src_ptr + blockIndex,
        src_stride,
        refPtr + searchPositionIndex,
        reflumaStride,
        &p_best_sad8x8[40],
        &p_best_sad16x16[10],
        &p_best_mv8x8[40],
        &p_best_mv16x16[10],
        currMV,
        &p_sad16x16[10],
        &p_sad8x8[40],
        sub_sad);
    //---- 16x16 : 11
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    ext_sad_calculation_8x8_16x16(
        src_ptr + blockIndex,
        src_stride,
        refPtr + searchPositionIndex,
        reflumaStride,
        &p_best_sad8x8[44],
        &p_best_sad16x16[11],
        &p_best_mv8x8[44],
        &p_best_mv16x16[11],
        currMV,
        &p_sad16x16[11],
        &p_sad8x8[44],
        sub_sad);
    //---- 16x16 : 14
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    ext_sad_calculation_8x8_16x16(
        src_ptr + blockIndex,
        src_stride,
        refPtr + searchPositionIndex,
        reflumaStride,
        &p_best_sad8x8[56],
        &p_best_sad16x16[14],
        &p_best_mv8x8[56],
        &p_best_mv16x16[14],
        currMV,
        &p_sad16x16[14],
        &p_sad8x8[56],
        sub_sad);
    //---- 16x16 : 15
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    ext_sad_calculation_8x8_16x16(
        src_ptr + blockIndex,
        src_stride,
        refPtr + searchPositionIndex,
        reflumaStride,
        &p_best_sad8x8[60],
        &p_best_sad16x16[15],
        &p_best_mv8x8[60],
        &p_best_mv16x16[15],
        currMV,
        &p_sad16x16[15],
        &p_sad8x8[60],
        sub_sad);

    ext_sad_calculation_32x32_64x64(p_sad16x16,
                                    p_best_sad32x32,
                                    p_best_sad64x64,
                                    p_best_mv32x32,
                                    p_best_mv64x64,
                                    currMV,
                                    &p_sad32x32[0]);

    ExtSadCalculation(p_sad8x8,
                    p_sad16x16,
                    p_sad32x32,
                    p_best_sad64x32,
                    p_best_mv64x32,
                    p_best_sad32x16,
                    p_best_mv32x16,
                    p_best_sad16x8,
                    p_best_mv16x8,
                    p_best_sad32x64,
                    p_best_mv32x64,
                    p_best_sad16x32,
                    p_best_mv16x32,
                    p_best_sad8x16,
                    p_best_mv8x16,
                    p_best_sad32x8,
                    p_best_mv32x8,
                    p_best_sad8x32,
                    p_best_mv8x32,
                    p_best_sad64x16,
                    p_best_mv64x16,
                    p_best_sad16x64,
                    p_best_mv16x64,
                    currMV);
}

/*******************************************
 * GetSearchPointResults
 *******************************************/
static void GetSearchPointResults(
    MeContext
        *context_ptr,    // input parameter, ME context Ptr, used to get SB Ptr
    uint32_t listIndex,  // input parameter, reference list index
    uint32_t ref_pic_index,
    uint32_t searchRegionIndex,  // input parameter, search area origin, used to
                                 // point to reference samples
    int32_t xSearchIndex,  // input parameter, search region position in the
                           // horizontal direction, used to derive xMV
    int32_t ySearchIndex,  // input parameter, search region position in the
                           // vertical direction, used to derive yMV
    EbAsm asm_type) {
    const EbBool sub_sad = (context_ptr->me_search_method == SUB_SAD_SEARCH);
    uint8_t *src_ptr = context_ptr->sb_src_ptr;

    // uint8_t  *refPtr = refPicPtr->buffer_y; // NADER
    uint8_t *refPtr =
        context_ptr->integer_buffer_ptr[listIndex][ref_pic_index] +
        (ME_FILTER_TAP >> 1) +
        ((ME_FILTER_TAP >> 1) *
         context_ptr->interpolated_full_stride[listIndex][ref_pic_index]);
    // uint32_t reflumaStride = refPicPtr->stride_y; // NADER
    uint32_t reflumaStride =
        context_ptr->interpolated_full_stride[listIndex][ref_pic_index];

    uint32_t searchPositionTLIndex = searchRegionIndex;
    uint32_t searchPositionIndex;
    uint32_t blockIndex;

    uint32_t srcNext16x16Offset = (BLOCK_SIZE_64 << 4);
    // uint32_t refNext16x16Offset = (refPicPtr->stride_y << 4); // NADER
    uint32_t refNext16x16Offset = (reflumaStride << 4);

    uint32_t currMV1 = (((uint16_t)ySearchIndex) << 18);
    uint16_t currMV2 = (((uint16_t)xSearchIndex << 2));
    uint32_t currMV = currMV1 | currMV2;

    uint32_t *p_best_sad8x8 = context_ptr->p_best_sad8x8;
    uint32_t *p_best_sad16x16 = context_ptr->p_best_sad16x16;
    uint32_t *p_best_sad32x32 = context_ptr->p_best_sad32x32;
    uint32_t *p_best_sad64x64 = context_ptr->p_best_sad64x64;

    uint32_t *p_best_mv8x8 = context_ptr->p_best_mv8x8;
    uint32_t *p_best_mv16x16 = context_ptr->p_best_mv16x16;
    uint32_t *p_best_mv32x32 = context_ptr->p_best_mv32x32;
    uint32_t *p_best_mv64x64 = context_ptr->p_best_mv64x64;
    uint32_t *p_sad16x16 = context_ptr->p_sad16x16;

    // TODO: blockIndex searchPositionIndex could be removed  + Connect asm_type
    (void)asm_type;

    const uint32_t src_stride = context_ptr->sb_src_stride;
    srcNext16x16Offset = src_stride << 4;

    //---- 16x16 : 0
    blockIndex = 0;
    searchPositionIndex = searchPositionTLIndex;

    sad_calculation_8x8_16x16(
        src_ptr + blockIndex,
        src_stride,
        refPtr + searchPositionIndex,
        reflumaStride,
        &p_best_sad8x8[0],
        &p_best_sad16x16[0],
        &p_best_mv8x8[0],
        &p_best_mv16x16[0],
        currMV,
        &p_sad16x16[0],
        sub_sad);

    //---- 16x16 : 1
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionTLIndex + 16;
    sad_calculation_8x8_16x16(
        src_ptr + blockIndex,
        src_stride,
        refPtr + searchPositionIndex,
        reflumaStride,
        &p_best_sad8x8[4],
        &p_best_sad16x16[1],
        &p_best_mv8x8[4],
        &p_best_mv16x16[1],
        currMV,
        &p_sad16x16[1],
        sub_sad);
    //---- 16x16 : 4
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;

    sad_calculation_8x8_16x16(
        src_ptr + blockIndex,
        src_stride,
        refPtr + searchPositionIndex,
        reflumaStride,
        &p_best_sad8x8[16],
        &p_best_sad16x16[4],
        &p_best_mv8x8[16],
        &p_best_mv16x16[4],
        currMV,
        &p_sad16x16[4],
        sub_sad);

    //---- 16x16 : 5
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    sad_calculation_8x8_16x16(
        src_ptr + blockIndex,
        src_stride,
        refPtr + searchPositionIndex,
        reflumaStride,
        &p_best_sad8x8[20],
        &p_best_sad16x16[5],
        &p_best_mv8x8[20],
        &p_best_mv16x16[5],
        currMV,
        &p_sad16x16[5],
        sub_sad);

    //---- 16x16 : 2
    blockIndex = srcNext16x16Offset;
    searchPositionIndex = searchPositionTLIndex + refNext16x16Offset;
    sad_calculation_8x8_16x16(
        src_ptr + blockIndex,
        src_stride,
        refPtr + searchPositionIndex,
        reflumaStride,
        &p_best_sad8x8[8],
        &p_best_sad16x16[2],
        &p_best_mv8x8[8],
        &p_best_mv16x16[2],
        currMV,
        &p_sad16x16[2],
        sub_sad);
    //---- 16x16 : 3
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    sad_calculation_8x8_16x16(
        src_ptr + blockIndex,
        src_stride,
        refPtr + searchPositionIndex,
        reflumaStride,
        &p_best_sad8x8[12],
        &p_best_sad16x16[3],
        &p_best_mv8x8[12],
        &p_best_mv16x16[3],
        currMV,
        &p_sad16x16[3],
        sub_sad);
    //---- 16x16 : 6
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    sad_calculation_8x8_16x16(
        src_ptr + blockIndex,
        src_stride,
        refPtr + searchPositionIndex,
        reflumaStride,
        &p_best_sad8x8[24],
        &p_best_sad16x16[6],
        &p_best_mv8x8[24],
        &p_best_mv16x16[6],
        currMV,
        &p_sad16x16[6],
        sub_sad);
    //---- 16x16 : 7
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    sad_calculation_8x8_16x16(
        src_ptr + blockIndex,
        src_stride,
        refPtr + searchPositionIndex,
        reflumaStride,
        &p_best_sad8x8[28],
        &p_best_sad16x16[7],
        &p_best_mv8x8[28],
        &p_best_mv16x16[7],
        currMV,
        &p_sad16x16[7],
        sub_sad);

    //---- 16x16 : 8
    blockIndex = (srcNext16x16Offset << 1);
    searchPositionIndex = searchPositionTLIndex + (refNext16x16Offset << 1);
    sad_calculation_8x8_16x16(
        src_ptr + blockIndex,
        src_stride,
        refPtr + searchPositionIndex,
        reflumaStride,
        &p_best_sad8x8[32],
        &p_best_sad16x16[8],
        &p_best_mv8x8[32],
        &p_best_mv16x16[8],
        currMV,
        &p_sad16x16[8],
        sub_sad);
    //---- 16x16 : 9
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    sad_calculation_8x8_16x16(
        src_ptr + blockIndex,
        src_stride,
        refPtr + searchPositionIndex,
        reflumaStride,
        &p_best_sad8x8[36],
        &p_best_sad16x16[9],
        &p_best_mv8x8[36],
        &p_best_mv16x16[9],
        currMV,
        &p_sad16x16[9],
        sub_sad);
    //---- 16x16 : 12
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    sad_calculation_8x8_16x16(
        src_ptr + blockIndex,
        src_stride,
        refPtr + searchPositionIndex,
        reflumaStride,
        &p_best_sad8x8[48],
        &p_best_sad16x16[12],
        &p_best_mv8x8[48],
        &p_best_mv16x16[12],
        currMV,
        &p_sad16x16[12],
        sub_sad);
    //---- 16x16 : 13
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    sad_calculation_8x8_16x16(
        src_ptr + blockIndex,
        src_stride,
        refPtr + searchPositionIndex,
        reflumaStride,
        &p_best_sad8x8[52],
        &p_best_sad16x16[13],
        &p_best_mv8x8[52],
        &p_best_mv16x16[13],
        currMV,
        &p_sad16x16[13],
        sub_sad);

    //---- 16x16 : 10
    blockIndex = (srcNext16x16Offset * 3);
    searchPositionIndex = searchPositionTLIndex + (refNext16x16Offset * 3);
    sad_calculation_8x8_16x16(
        src_ptr + blockIndex,
        src_stride,
        refPtr + searchPositionIndex,
        reflumaStride,
        &p_best_sad8x8[40],
        &p_best_sad16x16[10],
        &p_best_mv8x8[40],
        &p_best_mv16x16[10],
        currMV,
        &p_sad16x16[10],
        sub_sad);
    //---- 16x16 : 11
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    sad_calculation_8x8_16x16(
        src_ptr + blockIndex,
        src_stride,
        refPtr + searchPositionIndex,
        reflumaStride,
        &p_best_sad8x8[44],
        &p_best_sad16x16[11],
        &p_best_mv8x8[44],
        &p_best_mv16x16[11],
        currMV,
        &p_sad16x16[11],
        sub_sad);
    //---- 16x16 : 14
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    sad_calculation_8x8_16x16(
        src_ptr + blockIndex,
        src_stride,
        refPtr + searchPositionIndex,
        reflumaStride,
        &p_best_sad8x8[56],
        &p_best_sad16x16[14],
        &p_best_mv8x8[56],
        &p_best_mv16x16[14],
        currMV,
        &p_sad16x16[14],
        sub_sad);
    //---- 16x16 : 15
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    sad_calculation_8x8_16x16(
        src_ptr + blockIndex,
        src_stride,
        refPtr + searchPositionIndex,
        reflumaStride,
        &p_best_sad8x8[60],
        &p_best_sad16x16[15],
        &p_best_mv8x8[60],
        &p_best_mv16x16[15],
        currMV,
        &p_sad16x16[15],
        sub_sad);

    sad_calculation_32x32_64x64(p_sad16x16,
                                p_best_sad32x32,
                                p_best_sad64x64,
                                p_best_mv32x32,
                                p_best_mv64x64,
                                currMV);
}

/*******************************************
 * GetEightHorizontalSearchPointResultsAll85CUs
 *******************************************/
static void GetEightHorizontalSearchPointResultsAll85PUs(
    MeContext *context_ptr, uint32_t listIndex,
    uint32_t ref_pic_index,
    uint32_t searchRegionIndex,
    int32_t xSearchIndex,  // input parameter, search region position in the
                           // horizontal direction, used to derive xMV
    int32_t ySearchIndex) {  // input parameter, search region position in the
                             // vertical direction, used to derive yMV
    const EbBool sub_sad = (context_ptr->me_search_method == SUB_SAD_SEARCH);
    uint8_t *src_ptr = context_ptr->sb_src_ptr;
    uint8_t *refPtr =
        context_ptr->integer_buffer_ptr[listIndex][ref_pic_index] +
        (ME_FILTER_TAP >> 1) +
        ((ME_FILTER_TAP >> 1) *
         context_ptr->interpolated_full_stride[listIndex][ref_pic_index]);
    uint32_t reflumaStride =
        context_ptr->interpolated_full_stride[listIndex][ref_pic_index];

    uint32_t searchPositionTLIndex = searchRegionIndex;
    uint32_t searchPositionIndex;
    uint32_t blockIndex;

    uint32_t srcNext16x16Offset = (BLOCK_SIZE_64 << 4);
    uint32_t refNext16x16Offset = (reflumaStride << 4);

    uint32_t currMVy = (((uint16_t)ySearchIndex) << 18);
    uint16_t currMVx = (((uint16_t)xSearchIndex << 2));
    uint32_t currMV = currMVy | currMVx;

    uint32_t *p_best_sad8x8 = context_ptr->p_best_sad8x8;
    uint32_t *p_best_sad16x16 = context_ptr->p_best_sad16x16;
    uint32_t *p_best_sad32x32 = context_ptr->p_best_sad32x32;
    uint32_t *p_best_sad64x64 = context_ptr->p_best_sad64x64;

    uint32_t *p_best_mv8x8 = context_ptr->p_best_mv8x8;
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
    srcNext16x16Offset = src_stride << 4;

    //---- 16x16_0
    blockIndex = 0;
    searchPositionIndex = searchPositionTLIndex;
    get_eight_horizontal_search_point_results_8x8_16x16_pu(
                   src_ptr + blockIndex,
                   context_ptr->sb_src_stride,
                   refPtr + searchPositionIndex,
                   reflumaStride,
                   &p_best_sad8x8[0],
                   &p_best_mv8x8[0],
                   &p_best_sad16x16[0],
                   &p_best_mv16x16[0],
                   currMV,
                   &p_sad16x16[0 * 8],
                   sub_sad);
    //---- 16x16_1
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionTLIndex + 16;
    get_eight_horizontal_search_point_results_8x8_16x16_pu(
                   src_ptr + blockIndex,
                   context_ptr->sb_src_stride,
                   refPtr + searchPositionIndex,
                   reflumaStride,
                   &p_best_sad8x8[4],
                   &p_best_mv8x8[4],
                   &p_best_sad16x16[1],
                   &p_best_mv16x16[1],
                   currMV,
                   &p_sad16x16[1 * 8],
                   sub_sad);
    //---- 16x16_4
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    get_eight_horizontal_search_point_results_8x8_16x16_pu(
                   src_ptr + blockIndex,
                   context_ptr->sb_src_stride,
                   refPtr + searchPositionIndex,
                   reflumaStride,
                   &p_best_sad8x8[16],
                   &p_best_mv8x8[16],
                   &p_best_sad16x16[4],
                   &p_best_mv16x16[4],
                   currMV,
                   &p_sad16x16[4 * 8],
                   sub_sad);
    //---- 16x16_5
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    get_eight_horizontal_search_point_results_8x8_16x16_pu(
                   src_ptr + blockIndex,
                   context_ptr->sb_src_stride,
                   refPtr + searchPositionIndex,
                   reflumaStride,
                   &p_best_sad8x8[20],
                   &p_best_mv8x8[20],
                   &p_best_sad16x16[5],
                   &p_best_mv16x16[5],
                   currMV,
                   &p_sad16x16[5 * 8],
                   sub_sad);

    //---- 16x16_2
    blockIndex = srcNext16x16Offset;
    searchPositionIndex = searchPositionTLIndex + refNext16x16Offset;
    get_eight_horizontal_search_point_results_8x8_16x16_pu(
                   src_ptr + blockIndex,
                   context_ptr->sb_src_stride,
                   refPtr + searchPositionIndex,
                   reflumaStride,
                   &p_best_sad8x8[8],
                   &p_best_mv8x8[8],
                   &p_best_sad16x16[2],
                   &p_best_mv16x16[2],
                   currMV,
                   &p_sad16x16[2 * 8],
                   sub_sad);
    //---- 16x16_3
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    get_eight_horizontal_search_point_results_8x8_16x16_pu(
                   src_ptr + blockIndex,
                   context_ptr->sb_src_stride,
                   refPtr + searchPositionIndex,
                   reflumaStride,
                   &p_best_sad8x8[12],
                   &p_best_mv8x8[12],
                   &p_best_sad16x16[3],
                   &p_best_mv16x16[3],
                   currMV,
                   &p_sad16x16[3 * 8],
                   sub_sad);
    //---- 16x16_6
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    get_eight_horizontal_search_point_results_8x8_16x16_pu(
                   src_ptr + blockIndex,
                   context_ptr->sb_src_stride,
                   refPtr + searchPositionIndex,
                   reflumaStride,
                   &p_best_sad8x8[24],
                   &p_best_mv8x8[24],
                   &p_best_sad16x16[6],
                   &p_best_mv16x16[6],
                   currMV,
                   &p_sad16x16[6 * 8],
                   sub_sad);
    //---- 16x16_7
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    get_eight_horizontal_search_point_results_8x8_16x16_pu(
                   src_ptr + blockIndex,
                   context_ptr->sb_src_stride,
                   refPtr + searchPositionIndex,
                   reflumaStride,
                   &p_best_sad8x8[28],
                   &p_best_mv8x8[28],
                   &p_best_sad16x16[7],
                   &p_best_mv16x16[7],
                   currMV,
                   &p_sad16x16[7 * 8],
                   sub_sad);

    //---- 16x16_8
    blockIndex = (srcNext16x16Offset << 1);
    searchPositionIndex = searchPositionTLIndex + (refNext16x16Offset << 1);
    get_eight_horizontal_search_point_results_8x8_16x16_pu(
                   src_ptr + blockIndex,
                   context_ptr->sb_src_stride,
                   refPtr + searchPositionIndex,
                   reflumaStride,
                   &p_best_sad8x8[32],
                   &p_best_mv8x8[32],
                   &p_best_sad16x16[8],
                   &p_best_mv16x16[8],
                   currMV,
                   &p_sad16x16[8 * 8],
                   sub_sad);
    //---- 16x16_9
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    get_eight_horizontal_search_point_results_8x8_16x16_pu(
                   src_ptr + blockIndex,
                   context_ptr->sb_src_stride,
                   refPtr + searchPositionIndex,
                   reflumaStride,
                   &p_best_sad8x8[36],
                   &p_best_mv8x8[36],
                   &p_best_sad16x16[9],
                   &p_best_mv16x16[9],
                   currMV,
                   &p_sad16x16[9 * 8],
                   sub_sad);
    //---- 16x16_12
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    get_eight_horizontal_search_point_results_8x8_16x16_pu(
                   src_ptr + blockIndex,
                   context_ptr->sb_src_stride,
                   refPtr + searchPositionIndex,
                   reflumaStride,
                   &p_best_sad8x8[48],
                   &p_best_mv8x8[48],
                   &p_best_sad16x16[12],
                   &p_best_mv16x16[12],
                   currMV,
                   &p_sad16x16[12 * 8],
                   sub_sad);
    //---- 16x1_13
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    get_eight_horizontal_search_point_results_8x8_16x16_pu(
                   src_ptr + blockIndex,
                   context_ptr->sb_src_stride,
                   refPtr + searchPositionIndex,
                   reflumaStride,
                   &p_best_sad8x8[52],
                   &p_best_mv8x8[52],
                   &p_best_sad16x16[13],
                   &p_best_mv16x16[13],
                   currMV,
                   &p_sad16x16[13 * 8],
                   sub_sad);

    //---- 16x16_10
    blockIndex = (srcNext16x16Offset * 3);
    searchPositionIndex = searchPositionTLIndex + (refNext16x16Offset * 3);
    get_eight_horizontal_search_point_results_8x8_16x16_pu(
                   src_ptr + blockIndex,
                   context_ptr->sb_src_stride,
                   refPtr + searchPositionIndex,
                   reflumaStride,
                   &p_best_sad8x8[40],
                   &p_best_mv8x8[40],
                   &p_best_sad16x16[10],
                   &p_best_mv16x16[10],
                   currMV,
                   &p_sad16x16[10 * 8],
                   sub_sad);
    //---- 16x16_11
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    get_eight_horizontal_search_point_results_8x8_16x16_pu(
                   src_ptr + blockIndex,
                   context_ptr->sb_src_stride,
                   refPtr + searchPositionIndex,
                   reflumaStride,
                   &p_best_sad8x8[44],
                   &p_best_mv8x8[44],
                   &p_best_sad16x16[11],
                   &p_best_mv16x16[11],
                   currMV,
                   &p_sad16x16[11 * 8],
                   sub_sad);
    //---- 16x16_14
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    get_eight_horizontal_search_point_results_8x8_16x16_pu(
                   src_ptr + blockIndex,
                   context_ptr->sb_src_stride,
                   refPtr + searchPositionIndex,
                   reflumaStride,
                   &p_best_sad8x8[56],
                   &p_best_mv8x8[56],
                   &p_best_sad16x16[14],
                   &p_best_mv16x16[14],
                   currMV,
                   &p_sad16x16[14 * 8],
                   sub_sad);
    //---- 16x16_15
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    get_eight_horizontal_search_point_results_8x8_16x16_pu(
                   src_ptr + blockIndex,
                   context_ptr->sb_src_stride,
                   refPtr + searchPositionIndex,
                   reflumaStride,
                   &p_best_sad8x8[60],
                   &p_best_mv8x8[60],
                   &p_best_sad16x16[15],
                   &p_best_mv16x16[15],
                   currMV,
                   &p_sad16x16[15 * 8],
                   sub_sad);
    // 32x32 and 64x64
    get_eight_horizontal_search_point_results_32x32_64x64_pu(
                   p_sad16x16,
                   p_best_sad32x32,
                   p_best_sad64x64,
                   p_best_mv32x32,
                   p_best_mv64x64,
                   currMV);
}

/*******************************************
 * FullPelSearch_LCU
 *******************************************/
static void FullPelSearch_LCU(MeContext *context_ptr, uint32_t listIndex,
                              uint32_t ref_pic_index,
                              int16_t x_search_area_origin,
                              int16_t y_search_area_origin,
                              uint32_t search_area_width,
                              uint32_t search_area_height, EbAsm asm_type) {
    uint32_t xSearchIndex, ySearchIndex;

    uint32_t searchAreaWidthRest8 = search_area_width & 7;
    uint32_t searchAreaWidthMult8 = search_area_width - searchAreaWidthRest8;

    for (ySearchIndex = 0; ySearchIndex < search_area_height; ySearchIndex++) {
        for (xSearchIndex = 0; xSearchIndex < searchAreaWidthMult8;
             xSearchIndex += 8) {
            // this function will do:  xSearchIndex, +1, +2, ..., +7
            GetEightHorizontalSearchPointResultsAll85PUs(
                context_ptr,
                listIndex,
                ref_pic_index,
                xSearchIndex +
                    ySearchIndex *
                        context_ptr->interpolated_full_stride[listIndex]
                                                             [ref_pic_index],
                (int32_t)xSearchIndex + x_search_area_origin,
                (int32_t)ySearchIndex + y_search_area_origin);
        }

        for (xSearchIndex = searchAreaWidthMult8;
             xSearchIndex < search_area_width;
             xSearchIndex++) {
            GetSearchPointResults(
                context_ptr,
                listIndex,
                ref_pic_index,
                xSearchIndex +
                    ySearchIndex *
                        context_ptr->interpolated_full_stride[listIndex]
                                                             [ref_pic_index],
                (int32_t)xSearchIndex + x_search_area_origin,
                (int32_t)ySearchIndex + y_search_area_origin,
                asm_type);
        }
    }
}
#if OPTIMISED_EX_SUBPEL
/*******************************************
 * PU_HalfPelRefinement
 *   performs Half Pel refinement for one PU
 *******************************************/
static void half_pel_refinement_block(
    MeContext
    *context_ptr,  // input parameter, ME context Ptr, used to get SB Ptr
    uint8_t *ref_buffer, uint32_t ref_stride, uint32_t *p_best_ssd,
    uint32_t src_block_index,  // input parameter, PU origin, used to point to
                               // source samples
    uint8_t *pos_b_buffer,  // input parameter, position "b" interpolated search
                            // area Ptr
    uint8_t *pos_h_buffer,  // input parameter, position "h" interpolated search
                            // area Ptr
    uint8_t *pos_j_buffer,  // input parameter, position "j" interpolated search
                            // area Ptr
    uint32_t pu_width,      // input parameter, PU width
    uint32_t pu_height,     // input parameter, PU height
    int16_t x_search_area_origin,  // input parameter, search area origin in the
                                   // horizontal direction, used to point to
                                   // reference samples
    int16_t y_search_area_origin,  // input parameter, search area origin in the
                                   // vertical direction, used to point to
                                   // reference samples
#if OPTIMISED_EX_SUBPEL
    uint32_t search_area_height,  // input parameter, search area height
    uint32_t search_area_width,  // input parameter, search area width
#endif
    uint32_t *p_best_sad, uint32_t *p_best_mv,
    uint8_t *p_sub_pel_direction, uint32_t *best_pervious_stage_mv,
    uint32_t ineteger_mv) {
    int32_t search_region_index;
    uint64_t distortion_left_position = 0;
    uint64_t distortion_top_position = 0;
    uint64_t distortion_topleft_position = 0;
    uint64_t distortion_topright_position = 0;
    int16_t half_mv_x[8];
    int16_t half_mv_y[8];
    int16_t x_best_mv;
    int16_t y_best_mv;
    int16_t x_mv;
    int16_t y_mv;
    int16_t search_index_x;
    int16_t search_index_y;
    (void)p_sub_pel_direction;
    (void)ineteger_mv;
    // copute distance between best mv and the integer mv candidate
    int16_t offset_x, offset_y;
    for (offset_x = -H_PEL_SEARCH_WIND; offset_x <= H_PEL_SEARCH_WIND; offset_x++) {
        for (offset_y = -H_PEL_SEARCH_WIND; offset_y <= H_PEL_SEARCH_WIND; offset_y++) {
            x_best_mv = _MVXT(*best_pervious_stage_mv);
            y_best_mv = _MVYT(*best_pervious_stage_mv);
            x_mv = x_best_mv + (offset_x * 4);
            y_mv = y_best_mv + (offset_y * 4);
            search_index_x = (x_mv >> 2) - x_search_area_origin;
            search_index_y = (y_mv >> 2) - y_search_area_origin;
            uint32_t integer_mv1 = (((uint16_t)(y_mv >> 2)) << 18);
            uint16_t integer_mv2 = (((uint16_t)(x_mv >> 2) << 2));
            uint32_t integer_mv = integer_mv1 | integer_mv2;
            if (search_index_x < 0 || search_index_x >(int16_t)(search_area_width - 1)) {
                continue;
            }
            if (search_index_y < 0 || search_index_y >(int16_t)(search_area_height - 1)) {
                continue;
            }
            half_mv_x[0] = x_mv - 2;  // L  position
            half_mv_x[1] = x_mv + 2;  // R  position
            half_mv_x[2] = x_mv;      // T  position
            half_mv_x[3] = x_mv;      // B  position
            half_mv_x[4] = x_mv - 2;  // TL position
            half_mv_x[5] = x_mv + 2;  // TR position
            half_mv_x[6] = x_mv + 2;  // BR position
            half_mv_x[7] = x_mv - 2;  // BL position
            half_mv_y[0] = y_mv;      // L  position
            half_mv_y[1] = y_mv;      // R  position
            half_mv_y[2] = y_mv - 2;  // T  position
            half_mv_y[3] = y_mv + 2;  // B  position
            half_mv_y[4] = y_mv - 2;  // TL position
            half_mv_y[5] = y_mv - 2;  // TR position
            half_mv_y[6] = y_mv + 2;  // BR position
            half_mv_y[7] = y_mv + 2;  // BL position
            // Compute SSD for the best full search candidate
            if (context_ptr->fractional_search_method == SSD_SEARCH) {
                uint32_t integer_sse =
                    (uint32_t)spatial_full_distortion_kernel(
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
                    *p_best_mv = integer_mv;
                }
            }
            // L position
            search_region_index =
                search_index_x +
                (int16_t)context_ptr->interpolated_stride * search_index_y;
            if (context_ptr->fractional_search_method == SSD_SEARCH)
                distortion_left_position = spatial_full_distortion_kernel(
                    context_ptr->sb_src_ptr,
                    src_block_index,
                    context_ptr->sb_src_stride,
                    pos_b_buffer,
                    search_region_index,
                    context_ptr->interpolated_stride,
                    pu_width,
                    pu_height);
            else if (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                distortion_left_position = (nxm_sad_kernel(
                    &(context_ptr->sb_src_ptr[src_block_index]),
                    context_ptr->sb_src_stride << 1,
                    &(pos_b_buffer[search_region_index]),
                    context_ptr->interpolated_stride << 1,
                    pu_height >> 1,
                    pu_width)) << 1;
            else
                distortion_left_position = nxm_sad_kernel(
                    &(context_ptr->sb_src_ptr[src_block_index]),
                    context_ptr->sb_src_stride,
                    &(pos_b_buffer[search_region_index]),
                    context_ptr->interpolated_stride,
                    pu_height,
                    pu_width);
            if (context_ptr->fractional_search_method == SSD_SEARCH) {
                if (distortion_left_position < *p_best_ssd) {
                    *p_best_sad = (uint32_t)
                        nxm_sad_kernel(
                            &(context_ptr->sb_src_ptr[src_block_index]),
                            context_ptr->sb_src_stride,
                            &(pos_b_buffer[search_region_index]),
                            context_ptr->interpolated_stride,
                            pu_height,
                            pu_width);
                    *p_best_mv =
                        ((uint16_t)half_mv_y[0] << 16) | ((uint16_t)half_mv_x[0]);
                    *p_best_ssd = (uint32_t)distortion_left_position;
                }
            }
            else {
                if (distortion_left_position < *p_best_sad) {
                    *p_best_sad = (uint32_t)distortion_left_position;
                    *p_best_mv =
                        ((uint16_t)half_mv_y[0] << 16) | ((uint16_t)half_mv_x[0]);
                }
            }
            // T position
            search_region_index =
                search_index_x +
                (int16_t)context_ptr->interpolated_stride * search_index_y;
            if (context_ptr->fractional_search_method == SSD_SEARCH)
                distortion_top_position = spatial_full_distortion_kernel(
                    context_ptr->sb_src_ptr,
                    src_block_index,
                    context_ptr->sb_src_stride,
                    pos_h_buffer,
                    search_region_index,
                    context_ptr->interpolated_stride,
                    pu_width,
                    pu_height);
            else if (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                distortion_top_position = (nxm_sad_kernel(
                    &(context_ptr->sb_src_ptr[src_block_index]),
                    context_ptr->sb_src_stride << 1,
                    &(pos_h_buffer[search_region_index]),
                    context_ptr->interpolated_stride << 1,
                    pu_height >> 1,
                    pu_width)) << 1;
            else
                distortion_top_position = nxm_sad_kernel(
                    &(context_ptr->sb_src_ptr[src_block_index]),
                    context_ptr->sb_src_stride,
                    &(pos_h_buffer[search_region_index]),
                    context_ptr->interpolated_stride,
                    pu_height,
                    pu_width);
            if (context_ptr->fractional_search_method == SSD_SEARCH) {
                if (distortion_top_position < *p_best_ssd) {
                    *p_best_sad = (uint32_t)
                        nxm_sad_kernel(
                            &(context_ptr->sb_src_ptr[src_block_index]),
                            context_ptr->sb_src_stride,
                            &(pos_h_buffer[search_region_index]),
                            context_ptr->interpolated_stride,
                            pu_height,
                            pu_width);
                    *p_best_mv =
                        ((uint16_t)half_mv_y[2] << 16) | ((uint16_t)half_mv_x[2]);
                    *p_best_ssd = (uint32_t)distortion_top_position;
                }
            }
            else {
                if (distortion_top_position < *p_best_sad) {
                    *p_best_sad = (uint32_t)distortion_top_position;
                    *p_best_mv =
                        ((uint16_t)half_mv_y[2] << 16) | ((uint16_t)half_mv_x[2]);
                }
            }
            // TL position
            search_region_index =
                search_index_x +
                (int16_t)context_ptr->interpolated_stride * search_index_y;
            if (context_ptr->fractional_search_method == SSD_SEARCH)
                distortion_topleft_position = spatial_full_distortion_kernel(
                    context_ptr->sb_src_ptr,
                    src_block_index,
                    context_ptr->sb_src_stride,
                    pos_j_buffer,
                    search_region_index,
                    context_ptr->interpolated_stride,
                    pu_width,
                    pu_height);
            else if (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                distortion_topleft_position = (nxm_sad_kernel(
                    &(context_ptr->sb_src_ptr[src_block_index]),
                    context_ptr->sb_src_stride << 1,
                    &(pos_j_buffer[search_region_index]),
                    context_ptr->interpolated_stride << 1,
                    pu_height >> 1,
                    pu_width)) << 1;
            else
                distortion_topleft_position = nxm_sad_kernel(
                    &(context_ptr->sb_src_ptr[src_block_index]),
                    context_ptr->sb_src_stride,
                    &(pos_j_buffer[search_region_index]),
                    context_ptr->interpolated_stride,
                    pu_height,
                    pu_width);
            if (context_ptr->fractional_search_method == SSD_SEARCH) {
                if (distortion_topleft_position < *p_best_ssd) {
                    *p_best_sad = (uint32_t)
                        nxm_sad_kernel(
                            &(context_ptr->sb_src_ptr[src_block_index]),
                            context_ptr->sb_src_stride,
                            &(pos_j_buffer[search_region_index]),
                            context_ptr->interpolated_stride,
                            pu_height,
                            pu_width);
                    *p_best_mv =
                        ((uint16_t)half_mv_y[4] << 16) | ((uint16_t)half_mv_x[4]);
                    *p_best_ssd = (uint32_t)distortion_topleft_position;
                }
            }
            else {
                if (distortion_topleft_position < *p_best_sad) {
                    *p_best_sad = (uint32_t)distortion_topleft_position;
                    *p_best_mv =
                        ((uint16_t)half_mv_y[4] << 16) | ((uint16_t)half_mv_x[4]);
                }
            }
            // TR position
            search_region_index++;
            if (context_ptr->fractional_search_method == SSD_SEARCH)
                distortion_topright_position = spatial_full_distortion_kernel(
                    context_ptr->sb_src_ptr,
                    src_block_index,
                    context_ptr->sb_src_stride,
                    pos_j_buffer,
                    search_region_index,
                    context_ptr->interpolated_stride,
                    pu_width,
                    pu_height);
            else if (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                distortion_topright_position = (nxm_sad_kernel(
                    &(context_ptr->sb_src_ptr[src_block_index]),
                    context_ptr->sb_src_stride << 1,
                    &(pos_j_buffer[search_region_index]),
                    context_ptr->interpolated_stride << 1,
                    pu_height >> 1,
                    pu_width)) << 1;
            else
                distortion_topright_position = nxm_sad_kernel(
                    &(context_ptr->sb_src_ptr[src_block_index]),
                    context_ptr->sb_src_stride,
                    &(pos_j_buffer[search_region_index]),
                    context_ptr->interpolated_stride,
                    pu_height,
                    pu_width);
            if (context_ptr->fractional_search_method == SSD_SEARCH) {
                if (distortion_topright_position < *p_best_ssd) {
                    *p_best_sad = (uint32_t)
                        nxm_sad_kernel(
                            &(context_ptr->sb_src_ptr[src_block_index]),
                            context_ptr->sb_src_stride,
                            &(pos_j_buffer[search_region_index]),
                            context_ptr->interpolated_stride,
                            pu_height,
                            pu_width);
                    *p_best_mv =
                        ((uint16_t)half_mv_y[5] << 16) | ((uint16_t)half_mv_x[5]);
                    *p_best_ssd = (uint32_t)distortion_topright_position;
                }
            }
            else {
                if (distortion_topright_position < *p_best_sad) {
                    *p_best_sad = (uint32_t)distortion_topright_position;
                    *p_best_mv =
                        ((uint16_t)half_mv_y[5] << 16) | ((uint16_t)half_mv_x[5]);
                }
            }
        }
    }
    return;
}
#else
/*******************************************
 * PU_HalfPelRefinement
 *   performs Half Pel refinement for one PU
 *******************************************/
static void half_pel_refinement_block(
    MeContext
        *context_ptr,  // input parameter, ME context Ptr, used to get SB Ptr
    uint8_t *ref_buffer, uint32_t ref_stride, uint32_t *p_best_ssd,
    uint32_t src_block_index,  // input parameter, PU origin, used to point to
                               // source samples
    uint8_t *pos_b_buffer,  // input parameter, position "b" interpolated search
                            // area Ptr
    uint8_t *pos_h_buffer,  // input parameter, position "h" interpolated search
                            // area Ptr
    uint8_t *pos_j_buffer,  // input parameter, position "j" interpolated search
                            // area Ptr
    uint32_t pu_width,      // input parameter, PU width
    uint32_t pu_height,     // input parameter, PU height
    int16_t x_search_area_origin,  // input parameter, search area origin in the
                                   // horizontal direction, used to point to
                                   // reference samples
    int16_t y_search_area_origin,  // input parameter, search area origin in the
                                   // vertical direction, used to point to
                                   // reference samples
    uint32_t *p_best_sad, uint32_t *p_best_mv,
    uint8_t *p_sub_pel_direction, uint32_t *best_pervious_stage_mv,
    uint32_t ineteger_mv) {
    int32_t search_region_index;
    uint64_t distortion_left_position = 0;
    uint64_t distortion_top_position = 0;
    uint64_t distortion_topleft_position = 0;
    uint64_t distortion_topright_position = 0;
    int16_t half_mv_x[8];
    int16_t half_mv_y[8];
    // copute distance between best mv and the integer mv candidate
    int16_t int_x_mv = _MVXT(ineteger_mv);
    int16_t int_y_mv = _MVYT(ineteger_mv);
    int16_t int_search_index_x = (int_x_mv >> 2) - x_search_area_origin;
    int16_t int_search_index_y = (int_y_mv >> 2) - y_search_area_origin;
    int16_t x_best_mv = _MVXT(*best_pervious_stage_mv);
    int16_t y_best_mv = _MVYT(*best_pervious_stage_mv);
    int16_t best_search_index_x = (x_best_mv >> 2) - x_search_area_origin;
    int16_t best_search_index_y = (y_best_mv >> 2) - y_search_area_origin;
    int16_t dis_x = ABS(int_search_index_x - best_search_index_x);
    int16_t dis_y = ABS(int_search_index_y - best_search_index_y);
    // Skip half pel if the integer candidate is not inside the desired window.
    if ((dis_x) > H_PEL_SEARCH_WIND)
        return;
    if ((dis_y) > H_PEL_SEARCH_WIND)
        return;
    int16_t x_mv = _MVXT(ineteger_mv);
    int16_t y_mv = _MVYT(ineteger_mv);
    int16_t search_index_x = (x_mv >> 2) - x_search_area_origin;
    int16_t search_index_y = (y_mv >> 2) - y_search_area_origin;
    (void)p_sub_pel_direction;
    half_mv_x[0] = x_mv - 2;  // L  position
    half_mv_x[1] = x_mv + 2;  // R  position
    half_mv_x[2] = x_mv;      // T  position
    half_mv_x[3] = x_mv;      // B  position
    half_mv_x[4] = x_mv - 2;  // TL position
    half_mv_x[5] = x_mv + 2;  // TR position
    half_mv_x[6] = x_mv + 2;  // BR position
    half_mv_x[7] = x_mv - 2;  // BL position
    half_mv_y[0] = y_mv;      // L  position
    half_mv_y[1] = y_mv;      // R  position
    half_mv_y[2] = y_mv - 2;  // T  position
    half_mv_y[3] = y_mv + 2;  // B  position
    half_mv_y[4] = y_mv - 2;  // TL position
    half_mv_y[5] = y_mv - 2;  // TR position
    half_mv_y[6] = y_mv + 2;  // BR position
    half_mv_y[7] = y_mv + 2;  // BL position
    // Compute SSD for the best full search candidate
    if (context_ptr->fractional_search_method == SSD_SEARCH) {
        uint32_t integer_sse =
            (uint32_t)spatial_full_distortion_kernel(
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
            *p_best_mv = ineteger_mv;
        }
    }
    // L position
    search_region_index =
        search_index_x +
        (int16_t)context_ptr->interpolated_stride * search_index_y;
    distortion_left_position =
        (context_ptr->fractional_search_method == SSD_SEARCH)
            ? spatial_full_distortion_kernel(
                      context_ptr->sb_src_ptr,
                      src_block_index,
                      context_ptr->sb_src_stride,
                      pos_b_buffer,
                      search_region_index,
                      context_ptr->interpolated_stride,
                      pu_width,
                      pu_height)
            : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                  ? (nxm_sad_kernel(
                        &(context_ptr->sb_src_ptr[src_block_index]),
                        context_ptr->sb_src_stride << 1,
                        &(pos_b_buffer[search_region_index]),
                        context_ptr->interpolated_stride << 1,
                        pu_height >> 1,
                        pu_width))
                        << 1
                  : nxm_sad_kernel(
                        &(context_ptr->sb_src_ptr[src_block_index]),
                        context_ptr->sb_src_stride,
                        &(pos_b_buffer[search_region_index]),
                        context_ptr->interpolated_stride,
                        pu_height,
                        pu_width);
    if (context_ptr->fractional_search_method == SSD_SEARCH) {
        if (distortion_left_position < *p_best_ssd) {
            *p_best_sad = (uint32_t)
                nxm_sad_kernel(
                    &(context_ptr->sb_src_ptr[src_block_index]),
                    context_ptr->sb_src_stride,
                    &(pos_b_buffer[search_region_index]),
                    context_ptr->interpolated_stride,
                    pu_height,
                    pu_width);
            *p_best_mv =
                ((uint16_t)half_mv_y[0] << 16) | ((uint16_t)half_mv_x[0]);
            *p_best_ssd = (uint32_t)distortion_left_position;
        }
    } else {
        if (distortion_left_position < *p_best_sad) {
            *p_best_sad = (uint32_t)distortion_left_position;
            *p_best_mv =
                ((uint16_t)half_mv_y[0] << 16) | ((uint16_t)half_mv_x[0]);
        }
    }
#if !HP_REF_OPT
    // R position
    search_region_index++;
    distortion_right_position =
        (context_ptr->fractional_search_method == SSD_SEARCH)
            ? spatial_full_distortion_kernel(
                      &(context_ptr->sb_src_ptr[src_block_index]),
                      context_ptr->sb_src_stride,
                      &(pos_b_buffer[search_region_index]),
                      context_ptr->interpolated_stride,
                      pu_width,
                      pu_height)
            : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                  ? (nxm_sad_kernel(
                        &(context_ptr->sb_src_ptr[src_block_index]),
                        context_ptr->sb_src_stride << 1,
                        &(pos_b_buffer[search_region_index]),
                        context_ptr->interpolated_stride << 1,
                        pu_height >> 1,
                        pu_width))
                        << 1
                  : nxm_sad_kernel(
                        &(context_ptr->sb_src_ptr[src_block_index]),
                        context_ptr->sb_src_stride,
                        &(pos_b_buffer[search_region_index]),
                        context_ptr->interpolated_stride,
                        pu_height,
                        pu_width);
    if (context_ptr->fractional_search_method == SSD_SEARCH) {
        if (distortion_right_position < *p_best_ssd) {
            *p_best_sad = (uint32_t)
                nxm_sad_kernel(
                    &(context_ptr->sb_src_ptr[src_block_index]),
                    context_ptr->sb_src_stride,
                    &(pos_b_buffer[search_region_index]),
                    context_ptr->interpolated_stride,
                    pu_height,
                    pu_width);
            *p_best_mv =
                ((uint16_t)half_mv_y[1] << 16) | ((uint16_t)half_mv_x[1]);
            *p_best_ssd = (uint32_t)distortion_right_position;
        }
    } else {
        if (distortion_right_position < *p_best_sad) {
            *p_best_sad = (uint32_t)distortion_right_position;
            *p_best_mv =
                ((uint16_t)half_mv_y[1] << 16) | ((uint16_t)half_mv_x[1]);
        }
    }
#endif
    // T position
    search_region_index =
        search_index_x +
        (int16_t)context_ptr->interpolated_stride * search_index_y;
    distortion_top_position =
        (context_ptr->fractional_search_method == SSD_SEARCH)
            ? spatial_full_distortion_kernel(
                      context_ptr->sb_src_ptr,
                      src_block_index,
                      context_ptr->sb_src_stride,
                      pos_h_buffer,
                      search_region_index,
                      context_ptr->interpolated_stride,
                      pu_width,
                      pu_height)
            : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                  ? (nxm_sad_kernel(
                        &(context_ptr->sb_src_ptr[src_block_index]),
                        context_ptr->sb_src_stride << 1,
                        &(pos_h_buffer[search_region_index]),
                        context_ptr->interpolated_stride << 1,
                        pu_height >> 1,
                        pu_width))
                        << 1
                  : nxm_sad_kernel(
                        &(context_ptr->sb_src_ptr[src_block_index]),
                        context_ptr->sb_src_stride,
                        &(pos_h_buffer[search_region_index]),
                        context_ptr->interpolated_stride,
                        pu_height,
                        pu_width);
    if (context_ptr->fractional_search_method == SSD_SEARCH) {
        if (distortion_top_position < *p_best_ssd) {
            *p_best_sad = (uint32_t)
                nxm_sad_kernel(
                    &(context_ptr->sb_src_ptr[src_block_index]),
                    context_ptr->sb_src_stride,
                    &(pos_h_buffer[search_region_index]),
                    context_ptr->interpolated_stride,
                    pu_height,
                    pu_width);
            *p_best_mv =
                ((uint16_t)half_mv_y[2] << 16) | ((uint16_t)half_mv_x[2]);
            *p_best_ssd = (uint32_t)distortion_top_position;
        }
    } else {
        if (distortion_top_position < *p_best_sad) {
            *p_best_sad = (uint32_t)distortion_top_position;
            *p_best_mv =
                ((uint16_t)half_mv_y[2] << 16) | ((uint16_t)half_mv_x[2]);
        }
    }
#if !HP_REF_OPT
    // B position
    search_region_index += (int16_t)context_ptr->interpolated_stride;
    distortion_bottom_position =
        (context_ptr->fractional_search_method == SSD_SEARCH)
            ? spatial_full_distortion_kernel(
                      &(context_ptr->sb_src_ptr[src_block_index]),
                      context_ptr->sb_src_stride,
                      &(pos_h_buffer[search_region_index]),
                      context_ptr->interpolated_stride,
                      pu_width,
                      pu_height)
            : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                  ? (nxm_sad_kernel(
                        &(context_ptr->sb_src_ptr[src_block_index]),
                        context_ptr->sb_src_stride << 1,
                        &(pos_h_buffer[search_region_index]),
                        context_ptr->interpolated_stride << 1,
                        pu_height >> 1,
                        pu_width))
                        << 1
                  : nxm_sad_kernel(
                        &(context_ptr->sb_src_ptr[src_block_index]),
                        context_ptr->sb_src_stride,
                        &(pos_h_buffer[search_region_index]),
                        context_ptr->interpolated_stride,
                        pu_height,
                        pu_width);
    if (context_ptr->fractional_search_method == SSD_SEARCH) {
        if (distortion_bottom_position < *p_best_ssd) {
            *p_best_sad = (uint32_t)
                nxm_sad_kernel(
                    &(context_ptr->sb_src_ptr[src_block_index]),
                    context_ptr->sb_src_stride,
                    &(pos_h_buffer[search_region_index]),
                    context_ptr->interpolated_stride,
                    pu_height,
                    pu_width);
            *p_best_mv =
                ((uint16_t)half_mv_y[3] << 16) | ((uint16_t)half_mv_x[3]);
            *p_best_ssd = (uint32_t)distortion_bottom_position;
        }
    } else {
        if (distortion_bottom_position < *p_best_sad) {
            *p_best_sad = (uint32_t)distortionBottomPosition;
            *p_best_mv =
                ((uint16_t)half_mv_y[3] << 16) | ((uint16_t)half_mv_x[3]);
        }
    }
#endif
    // TL position
    search_region_index =
        search_index_x +
        (int16_t)context_ptr->interpolated_stride * search_index_y;
    distortion_topleft_position =
        (context_ptr->fractional_search_method == SSD_SEARCH)
            ? spatial_full_distortion_kernel(
                      context_ptr->sb_src_ptr,
                      src_block_index,
                      context_ptr->sb_src_stride,
                      pos_j_buffer,
                      search_region_index,
                      context_ptr->interpolated_stride,
                      pu_width,
                      pu_height)
            : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                  ? (nxm_sad_kernel(
                        &(context_ptr->sb_src_ptr[src_block_index]),
                        context_ptr->sb_src_stride << 1,
                        &(pos_j_buffer[search_region_index]),
                        context_ptr->interpolated_stride << 1,
                        pu_height >> 1,
                        pu_width))
                        << 1
                  : nxm_sad_kernel(
                        &(context_ptr->sb_src_ptr[src_block_index]),
                        context_ptr->sb_src_stride,
                        &(pos_j_buffer[search_region_index]),
                        context_ptr->interpolated_stride,
                        pu_height,
                        pu_width);
    if (context_ptr->fractional_search_method == SSD_SEARCH) {
        if (distortion_topleft_position < *p_best_ssd) {
            *p_best_sad = (uint32_t)
                nxm_sad_kernel(
                    &(context_ptr->sb_src_ptr[src_block_index]),
                    context_ptr->sb_src_stride,
                    &(pos_j_buffer[search_region_index]),
                    context_ptr->interpolated_stride,
                    pu_height,
                    pu_width);
            *p_best_mv =
                ((uint16_t)half_mv_y[4] << 16) | ((uint16_t)half_mv_x[4]);
            *p_best_ssd = (uint32_t)distortion_topleft_position;
        }
    } else {
        if (distortion_topleft_position < *p_best_sad) {
            *p_best_sad = (uint32_t)distortion_topleft_position;
            *p_best_mv =
                ((uint16_t)half_mv_y[4] << 16) | ((uint16_t)half_mv_x[4]);
        }
    }
    // TR position
    search_region_index++;
    distortion_topright_position =
        (context_ptr->fractional_search_method == SSD_SEARCH)
            ? spatial_full_distortion_kernel(
                      context_ptr->sb_src_ptr,
                      src_block_index,
                      context_ptr->sb_src_stride,
                      pos_j_buffer,
                      search_region_index,
                      context_ptr->interpolated_stride,
                      pu_width,
                      pu_height)
            : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                  ? (nxm_sad_kernel(
                        &(context_ptr->sb_src_ptr[src_block_index]),
                        context_ptr->sb_src_stride << 1,
                        &(pos_j_buffer[search_region_index]),
                        context_ptr->interpolated_stride << 1,
                        pu_height >> 1,
                        pu_width))
                        << 1
                  : nxm_sad_kernel(
                        &(context_ptr->sb_src_ptr[src_block_index]),
                        context_ptr->sb_src_stride,
                        &(pos_j_buffer[search_region_index]),
                        context_ptr->interpolated_stride,
                        pu_height,
                        pu_width);

    if (context_ptr->fractional_search_method == SSD_SEARCH) {
        if (distortion_topright_position < *p_best_ssd) {
            *p_best_sad = (uint32_t)
                nxm_sad_kernel(
                    &(context_ptr->sb_src_ptr[src_block_index]),
                    context_ptr->sb_src_stride,
                    &(pos_j_buffer[search_region_index]),
                    context_ptr->interpolated_stride,
                    pu_height,
                    pu_width);
            *p_best_mv =
                ((uint16_t)half_mv_y[5] << 16) | ((uint16_t)half_mv_x[5]);
            *p_best_ssd = (uint32_t)distortion_topright_position;
        }
    } else {
        if (distortion_topright_position < *p_best_sad) {
            *p_best_sad = (uint32_t)distortion_topright_position;
            *p_best_mv =
                ((uint16_t)half_mv_y[5] << 16) | ((uint16_t)half_mv_x[5]);
        }
    }
#if !HP_REF_OPT
    // BR position
    search_region_index += (int16_t)context_ptr->interpolated_stride;
    distortion_bottomright_position =
        (context_ptr->fractional_search_method == SSD_SEARCH)
            ? spatial_full_distortion_kernel(
                      context_ptr->sb_src_ptr,
                      src_block_index,
                      context_ptr->sb_src_stride,
                      pos_j_buffer,
                      search_region_index,
                      context_ptr->interpolated_stride,
                      pu_width,
                      pu_height)
            : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                  ? (nxm_sad_kernel(
                        &(context_ptr->sb_src_ptr[src_block_index]),
                        context_ptr->sb_src_stride << 1,
                        &(pos_j_buffer[search_region_index]),
                        context_ptr->interpolated_stride << 1,
                        pu_height >> 1,
                        pu_width))
                        << 1
                  : nxm_sad_kernel(
                        &(context_ptr->sb_src_ptr[src_block_index]),
                        context_ptr->sb_src_stride,
                        &(pos_j_buffer[search_region_index]),
                        context_ptr->interpolated_stride,
                        pu_height,
                        pu_width);
    if (context_ptr->fractional_search_method == SSD_SEARCH) {
        if (distortion_bottomright_position < *p_best_ssd) {
            *p_best_sad = (uint32_t)
                nxm_sad_kernel(
                    &(context_ptr->sb_src_ptr[src_block_index]),
                    context_ptr->sb_src_stride,
                    &(pos_j_buffer[search_region_index]),
                    context_ptr->interpolated_stride,
                    pu_height,
                    pu_width);
            *p_best_mv =
                ((uint16_t)half_mv_y[6] << 16) | ((uint16_t)half_mv_x[6]);
            *p_best_ssd = (uint32_t)distortion_bottomright_position;
        }
    } else {
        if (distortion_bottomright_position < *p_best_sad) {
            *p_best_sad = (uint32_t)distortion_bottomright_position;
            *p_best_mv =
                ((uint16_t)half_mv_y[6] << 16) | ((uint16_t)half_mv_x[6]);
        }
    }
    // BL position
    search_region_index--;
    distortion_bottomleft_position =
        (context_ptr->fractional_search_method == SSD_SEARCH)
            ? spatial_full_distortion_kernel(
                      &(context_ptr->sb_src_ptr[src_block_index]),
                      context_ptr->sb_src_stride,
                      &(pos_j_buffer[search_region_index]),
                      context_ptr->interpolated_stride,
                      pu_width,
                      pu_height)
            : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                  ? (nxm_sad_kernel(
                        &(context_ptr->sb_src_ptr[src_block_index]),
                        context_ptr->sb_src_stride << 1,
                        &(pos_j_buffer[search_region_index]),
                        context_ptr->interpolated_stride << 1,
                        pu_height >> 1,
                        pu_width))
                        << 1
                  : (nxm_sad_kernel(
                        &(context_ptr->sb_src_ptr[src_block_index]),
                        context_ptr->sb_src_stride,
                        &(pos_j_buffer[search_region_index]),
                        context_ptr->interpolated_stride,
                        pu_height,
                        pu_width));

    if (context_ptr->fractional_search_method == SSD_SEARCH) {
        if (distortion_bottomleft_position < *p_best_ssd) {
            *p_best_sad = (uint32_t)(
                nxm_sad_kernel(
                    &(context_ptr->sb_src_ptr[src_block_index]),
                    context_ptr->sb_src_stride,
                    &(pos_j_buffer[search_region_index]),
                    context_ptr->interpolated_stride,
                    pu_height,
                    pu_width));
            *p_best_mv =
                ((uint16_t)half_mv_y[7] << 16) | ((uint16_t)half_mv_x[7]);
            *p_best_ssd = (uint32_t)distortion_bottomleft_position;
        }
    } else {
        if (distortion_bottomleft_position < *p_best_sad) {
            *p_best_sad = (uint32_t)distortion_bottomleft_position;
            *p_best_mv =
                ((uint16_t)half_mv_y[7] << 16) | ((uint16_t)half_mv_x[7]);
        }
    }
#endif
    return;
}
#endif
/*******************************************
 * HalfPelSearch_LCU
 *   performs Half Pel refinement for the 85 PUs
 *******************************************/
void half_pel_refinement_sb(
    PictureParentControlSet *picture_control_set_ptr,
    MeContext *context_ptr,  // input/output parameter, ME context Ptr, used to
                             // get/update ME results
    uint8_t *refBuffer, uint32_t ref_stride,
    uint8_t *pos_b_buffer,  // input parameter, position "b" interpolated search
                            // area Ptr
    uint8_t *pos_h_buffer,  // input parameter, position "h" interpolated search
                            // area Ptr
    uint8_t *pos_j_buffer,  // input parameter, position "j" interpolated search
                            // area Ptr
    int16_t x_search_area_origin,  // input parameter, search area origin in the
                                   // horizontal direction, used to point to
                                   // reference samples
    int16_t y_search_area_origin,  // input parameter, search area origin in the
                                   // vertical direction, used to point to
                                   // reference samples
#if OPTIMISED_EX_SUBPEL
    uint32_t search_area_height,  // input parameter, search area height
    uint32_t search_area_width,  // input parameter, search area width
#endif
    uint32_t inetger_mv)
{
    uint32_t idx;
    uint32_t pu_index;
    uint32_t block_index_shift_x;
    uint32_t block_index_shift_y;
    uint32_t src_block_index;
    uint32_t posb_buffer_index;
    uint32_t posh_buffer_index;
    uint32_t posj_buffer_index;
    if (context_ptr->fractional_search64x64)
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
#if OPTIMISED_EX_SUBPEL
                                  search_area_height,
                                  search_area_width,
#endif
                                  context_ptr->p_best_sad64x64,
                                  context_ptr->p_best_mv64x64,
                                  &context_ptr->psub_pel_direction64x64,
                                  context_ptr->p_best_full_pel_mv64x64,
                                  inetger_mv);
    // 32x32 [4 partitions]
    for (pu_index = 0; pu_index < 4; ++pu_index) {
        block_index_shift_x = (pu_index & 0x01) << 5;
        block_index_shift_y = (pu_index >> 1) << 5;
        src_block_index = block_index_shift_x +
                          block_index_shift_y * context_ptr->sb_src_stride;
        posb_buffer_index =
            block_index_shift_x +
            block_index_shift_y * context_ptr->interpolated_stride;
        posh_buffer_index =
            block_index_shift_x +
            block_index_shift_y * context_ptr->interpolated_stride;
        posj_buffer_index =
            block_index_shift_x +
            block_index_shift_y * context_ptr->interpolated_stride;
        half_pel_refinement_block(
            context_ptr,
            &(refBuffer[block_index_shift_y * ref_stride +
                        block_index_shift_x]),
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
#if OPTIMISED_EX_SUBPEL
            search_area_height,
            search_area_width,
#endif
            &context_ptr->p_best_sad32x32[pu_index],
            &context_ptr->p_best_mv32x32[pu_index],
            &context_ptr->psub_pel_direction32x32[pu_index],
            &context_ptr->p_best_full_pel_mv32x32[pu_index],
            inetger_mv);
    }
    // 16x16 [16 partitions]
    for (pu_index = 0; pu_index < 16; ++pu_index) {
        idx = tab16x16[pu_index];
        block_index_shift_x = (pu_index & 0x03) << 4;
        block_index_shift_y = (pu_index >> 2) << 4;
        src_block_index = block_index_shift_x +
                          block_index_shift_y * context_ptr->sb_src_stride;
        posb_buffer_index =
            block_index_shift_x +
            block_index_shift_y * context_ptr->interpolated_stride;
        posh_buffer_index =
            block_index_shift_x +
            block_index_shift_y * context_ptr->interpolated_stride;
        posj_buffer_index =
            block_index_shift_x +
            block_index_shift_y * context_ptr->interpolated_stride;
        half_pel_refinement_block(context_ptr,
                                  &(refBuffer[block_index_shift_y * ref_stride +
                                              block_index_shift_x]),
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
#if OPTIMISED_EX_SUBPEL
                                  search_area_height,
                                  search_area_width,
#endif
                                  &context_ptr->p_best_sad16x16[idx],
                                  &context_ptr->p_best_mv16x16[idx],
                                  &context_ptr->psub_pel_direction16x16[idx],
                                  &context_ptr->p_best_full_pel_mv16x16[idx],
                                  inetger_mv);
    }
    // 8x8   [64 partitions]
    for (pu_index = 0; pu_index < 64; ++pu_index) {
        idx = tab8x8[pu_index];  // TODO bitwise this
        block_index_shift_x = (pu_index & 0x07) << 3;
        block_index_shift_y = (pu_index >> 3) << 3;
        src_block_index = block_index_shift_x +
                          block_index_shift_y * context_ptr->sb_src_stride;
        posb_buffer_index =
            block_index_shift_x +
            block_index_shift_y * context_ptr->interpolated_stride;
        posh_buffer_index =
            block_index_shift_x +
            block_index_shift_y * context_ptr->interpolated_stride;
        posj_buffer_index =
            block_index_shift_x +
            block_index_shift_y * context_ptr->interpolated_stride;
        half_pel_refinement_block(context_ptr,
                                  &(refBuffer[block_index_shift_y * ref_stride +
                                              block_index_shift_x]),
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
#if OPTIMISED_EX_SUBPEL
                                  search_area_height,
                                  search_area_width,
#endif
                                  &context_ptr->p_best_sad8x8[idx],
                                  &context_ptr->p_best_mv8x8[idx],
                                  &context_ptr->psub_pel_direction8x8[idx],
                                  &context_ptr->p_best_full_pel_mv8x8[idx],
                                  inetger_mv);
    }
    if (picture_control_set_ptr->pic_depth_mode <= PIC_ALL_C_DEPTH_MODE) {
        // 64x32
        for (pu_index = 0; pu_index < 2; ++pu_index) {
            block_index_shift_x = 0;
            block_index_shift_y = pu_index << 5;
            src_block_index = block_index_shift_x +
                              block_index_shift_y * context_ptr->sb_src_stride;
            posb_buffer_index =
                block_index_shift_x +
                block_index_shift_y * context_ptr->interpolated_stride;
            posh_buffer_index =
                block_index_shift_x +
                block_index_shift_y * context_ptr->interpolated_stride;
            posj_buffer_index =
                block_index_shift_x +
                block_index_shift_y * context_ptr->interpolated_stride;
            half_pel_refinement_block(
                context_ptr,
                &(refBuffer[block_index_shift_y * ref_stride +
                            block_index_shift_x]),
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
#if OPTIMISED_EX_SUBPEL
                search_area_height,
                search_area_width,
#endif
                &context_ptr->p_best_sad64x32[pu_index],
                &context_ptr->p_best_mv64x32[pu_index],
                &context_ptr->psub_pel_direction64x32[pu_index],
                &context_ptr->p_best_full_pel_mv64x32[pu_index],
                inetger_mv);
        }
        // 32x16
        for (pu_index = 0; pu_index < 8; ++pu_index) {
            idx = tab32x16[pu_index];  // TODO bitwise this
            block_index_shift_x = (pu_index & 0x01) << 5;
            block_index_shift_y = (pu_index >> 1) << 4;
            src_block_index = block_index_shift_x +
                              block_index_shift_y * context_ptr->sb_src_stride;
            posb_buffer_index =
                block_index_shift_x +
                block_index_shift_y * context_ptr->interpolated_stride;
            posh_buffer_index =
                block_index_shift_x +
                block_index_shift_y * context_ptr->interpolated_stride;
            posj_buffer_index =
                block_index_shift_x +
                block_index_shift_y * context_ptr->interpolated_stride;
            half_pel_refinement_block(
                context_ptr,
                &(refBuffer[block_index_shift_y * ref_stride +
                            block_index_shift_x]),
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
#if OPTIMISED_EX_SUBPEL
                search_area_height,
                search_area_width,
#endif
                &context_ptr->p_best_sad32x16[idx],
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
            src_block_index = block_index_shift_x +
                              block_index_shift_y * context_ptr->sb_src_stride;
            posb_buffer_index =
                block_index_shift_x +
                block_index_shift_y * context_ptr->interpolated_stride;
            posh_buffer_index =
                block_index_shift_x +
                block_index_shift_y * context_ptr->interpolated_stride;
            posj_buffer_index =
                block_index_shift_x +
                block_index_shift_y * context_ptr->interpolated_stride;
            half_pel_refinement_block(
                context_ptr,
                &(refBuffer[block_index_shift_y * ref_stride +
                            block_index_shift_x]),
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
#if OPTIMISED_EX_SUBPEL
                search_area_height,
                search_area_width,
#endif
                &context_ptr->p_best_sad16x8[idx],
                &context_ptr->p_best_mv16x8[idx],
                &context_ptr->psub_pel_direction16x8[idx],
                &context_ptr->p_best_full_pel_mv16x8[idx],
                inetger_mv);
        }
        // 32x64
        for (pu_index = 0; pu_index < 2; ++pu_index) {
            block_index_shift_x = pu_index << 5;
            block_index_shift_y = 0;
            src_block_index = block_index_shift_x +
                              block_index_shift_y * context_ptr->sb_src_stride;
            posb_buffer_index =
                block_index_shift_x +
                block_index_shift_y * context_ptr->interpolated_stride;
            posh_buffer_index =
                block_index_shift_x +
                block_index_shift_y * context_ptr->interpolated_stride;
            posj_buffer_index =
                block_index_shift_x +
                block_index_shift_y * context_ptr->interpolated_stride;
            half_pel_refinement_block(
                context_ptr,
                &(refBuffer[block_index_shift_y * ref_stride +
                            block_index_shift_x]),
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
#if OPTIMISED_EX_SUBPEL
                search_area_height,
                search_area_width,
#endif
                &context_ptr->p_best_sad32x64[pu_index],
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
            src_block_index = block_index_shift_x +
                              block_index_shift_y * context_ptr->sb_src_stride;
            posb_buffer_index =
                block_index_shift_x +
                block_index_shift_y * context_ptr->interpolated_stride;
            posh_buffer_index =
                block_index_shift_x +
                block_index_shift_y * context_ptr->interpolated_stride;
            posj_buffer_index =
                block_index_shift_x +
                block_index_shift_y * context_ptr->interpolated_stride;
            half_pel_refinement_block(
                context_ptr,
                &(refBuffer[block_index_shift_y * ref_stride +
                            block_index_shift_x]),
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
#if OPTIMISED_EX_SUBPEL
                search_area_height,
                search_area_width,
#endif
                &context_ptr->p_best_sad16x32[idx],
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
            src_block_index = block_index_shift_x +
                              block_index_shift_y * context_ptr->sb_src_stride;
            posb_buffer_index =
                block_index_shift_x +
                block_index_shift_y * context_ptr->interpolated_stride;
            posh_buffer_index =
                block_index_shift_x +
                block_index_shift_y * context_ptr->interpolated_stride;
            posj_buffer_index =
                block_index_shift_x +
                block_index_shift_y * context_ptr->interpolated_stride;
            half_pel_refinement_block(
                context_ptr,
                &(refBuffer[block_index_shift_y * ref_stride +
                            block_index_shift_x]),
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
#if OPTIMISED_EX_SUBPEL
                search_area_height,
                search_area_width,
#endif
                &context_ptr->p_best_sad8x16[idx],
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
            src_block_index = block_index_shift_x +
                              block_index_shift_y * context_ptr->sb_src_stride;
            posb_buffer_index =
                block_index_shift_x +
                block_index_shift_y * context_ptr->interpolated_stride;
            posh_buffer_index =
                block_index_shift_x +
                block_index_shift_y * context_ptr->interpolated_stride;
            posj_buffer_index =
                block_index_shift_x +
                block_index_shift_y * context_ptr->interpolated_stride;
            half_pel_refinement_block(
                context_ptr,
                &(refBuffer[block_index_shift_y * ref_stride +
                            block_index_shift_x]),
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
#if OPTIMISED_EX_SUBPEL
                search_area_height,
                search_area_width,
#endif
                &context_ptr->p_best_sad32x8[idx],
                &context_ptr->p_best_mv32x8[idx],
                &context_ptr->psub_pel_direction32x8[idx],
                &context_ptr->p_best_full_pel_mv32x8[idx],
                inetger_mv);
        }
        for (pu_index = 0; pu_index < 16; ++pu_index) {
            idx = tab8x32[pu_index];
            block_index_shift_x = (pu_index & 0x07) << 3;
            block_index_shift_y = (pu_index >> 3) << 5;
            src_block_index = block_index_shift_x +
                              block_index_shift_y * context_ptr->sb_src_stride;
            posb_buffer_index =
                block_index_shift_x +
                block_index_shift_y * context_ptr->interpolated_stride;
            posh_buffer_index =
                block_index_shift_x +
                block_index_shift_y * context_ptr->interpolated_stride;
            posj_buffer_index =
                block_index_shift_x +
                block_index_shift_y * context_ptr->interpolated_stride;
            half_pel_refinement_block(
                context_ptr,
                &(refBuffer[block_index_shift_y * ref_stride +
                            block_index_shift_x]),
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
#if OPTIMISED_EX_SUBPEL
                search_area_height,
                search_area_width,
#endif
                &context_ptr->p_best_sad8x32[idx],
                &context_ptr->p_best_mv8x32[idx],
                &context_ptr->psub_pel_direction8x32[idx],
                &context_ptr->p_best_full_pel_mv8x32[idx],
                inetger_mv);
        }
        for (pu_index = 0; pu_index < 4; ++pu_index) {
            idx = pu_index;
            block_index_shift_x = 0;
            block_index_shift_y = pu_index << 4;
            src_block_index = block_index_shift_x +
                              block_index_shift_y * context_ptr->sb_src_stride;
            posb_buffer_index =
                block_index_shift_x +
                block_index_shift_y * context_ptr->interpolated_stride;
            posh_buffer_index =
                block_index_shift_x +
                block_index_shift_y * context_ptr->interpolated_stride;
            posj_buffer_index =
                block_index_shift_x +
                block_index_shift_y * context_ptr->interpolated_stride;
            half_pel_refinement_block(
                context_ptr,
                &(refBuffer[block_index_shift_y * ref_stride +
                            block_index_shift_x]),
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
#if OPTIMISED_EX_SUBPEL
                search_area_height,
                search_area_width,
#endif
                &context_ptr->p_best_sad64x16[idx],
                &context_ptr->p_best_mv64x16[idx],
                &context_ptr->psub_pel_direction64x16[idx],
                &context_ptr->p_best_full_pel_mv64x16[idx],
                inetger_mv);
        }
        for (pu_index = 0; pu_index < 4; ++pu_index) {
            idx = pu_index;
            block_index_shift_x = pu_index << 4;
            block_index_shift_y = 0;
            src_block_index = block_index_shift_x +
                              block_index_shift_y * context_ptr->sb_src_stride;
            posb_buffer_index =
                block_index_shift_x +
                block_index_shift_y * context_ptr->interpolated_stride;
            posh_buffer_index =
                block_index_shift_x +
                block_index_shift_y * context_ptr->interpolated_stride;
            posj_buffer_index =
                block_index_shift_x +
                block_index_shift_y * context_ptr->interpolated_stride;
            half_pel_refinement_block(
                context_ptr,
                &(refBuffer[block_index_shift_y * ref_stride +
                            block_index_shift_x]),
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
#if OPTIMISED_EX_SUBPEL
                search_area_height,
                search_area_width,
#endif
                &context_ptr->p_best_sad16x64[idx],
                &context_ptr->p_best_mv16x64[idx],
                &context_ptr->psub_pel_direction16x64[idx],
                &context_ptr->p_best_full_pel_mv16x64[idx],
                inetger_mv);
        }
    }
    return;
}
/*******************************************
 * open_loop_me_half_pel_search_sblock
 *******************************************/
#if OPTIMISED_EX_SUBPEL
static void open_loop_me_half_pel_search_sblock(
    PictureParentControlSet *picture_control_set_ptr, MeContext *context_ptr,
    uint32_t list_index, uint32_t ref_pic_index, int16_t x_search_area_origin,
    int16_t y_search_area_origin, uint32_t search_area_width,
    uint32_t search_area_height)
{

    half_pel_refinement_sb(
        picture_control_set_ptr,
        context_ptr,
#if M0_HIGH_PRECISION_INTERPOLATION
        context_ptr->integer_buffer_ptr[list_index][ref_pic_index] +
        (ME_FILTER_PAD_DISTANCE >> 1) +
        ((ME_FILTER_PAD_DISTANCE >> 1) *
            context_ptr
            ->interpolated_full_stride[listIndex][ref_pic_index]),
        context_ptr
        ->interpolated_full_stride[list_index][ref_pic_index],
        &(context_ptr->pos_b_buffer[list_index][ref_pic_index]
            [(ME_FILTER_PAD_DISTANCE >> 1) *
            context_ptr->interpolated_stride]),
#else
        context_ptr->integer_buffer_ptr[list_index][ref_pic_index] +
        (ME_FILTER_TAP >> 1) +
        ((ME_FILTER_TAP >> 1) *
            context_ptr
            ->interpolated_full_stride[list_index][ref_pic_index]),
        context_ptr
        ->interpolated_full_stride[list_index][ref_pic_index],
        &(context_ptr->pos_b_buffer[list_index][ref_pic_index]
            [(ME_FILTER_TAP >> 1) *
            context_ptr->interpolated_stride]),
#endif
        &(context_ptr->pos_h_buffer[list_index][ref_pic_index][1]),
        &(context_ptr->pos_j_buffer[list_index][ref_pic_index][0]),
        x_search_area_origin,
        y_search_area_origin,
        search_area_height,
        search_area_width,
        0);
}
#else
static void open_loop_me_half_pel_search_sblock(
    PictureParentControlSet *picture_control_set_ptr, MeContext *context_ptr,
    uint32_t list_index, uint32_t ref_pic_index, int16_t x_search_area_origin,
    int16_t y_search_area_origin, uint32_t search_area_width,
    uint32_t search_area_height)
{
    uint32_t search_index_x, search_index_y;
    for (search_index_y = 0; search_index_y < search_area_height;
         search_index_y++) {
        for (search_index_x = 0; search_index_x < search_area_width;
             search_index_x++) {
            int32_t mvx = (int32_t)search_index_y + x_search_area_origin;
            int32_t mvy = (int32_t)search_index_x + y_search_area_origin;
            uint32_t inetger_mv1 = (((uint16_t)mvy) << 18);
            uint16_t inetger_mv2 = (((uint16_t)mvx << 2));
            uint32_t inetger_mv = inetger_mv1 | inetger_mv2;
            half_pel_refinement_sb(
                picture_control_set_ptr,
                context_ptr,
#if M0_HIGH_PRECISION_INTERPOLATION
                context_ptr->integer_buffer_ptr[list_index][ref_pic_index] +
                    (ME_FILTER_PAD_DISTANCE >> 1) +
                    ((ME_FILTER_PAD_DISTANCE >> 1) *
                     context_ptr
                         ->interpolated_full_stride[listIndex][ref_pic_index]),
                context_ptr
                    ->interpolated_full_stride[list_index][ref_pic_index],
                &(context_ptr->pos_b_buffer[list_index][ref_pic_index]
                                           [(ME_FILTER_PAD_DISTANCE >> 1) *
                                            context_ptr->interpolated_stride]),
#else
                context_ptr->integer_buffer_ptr[list_index][ref_pic_index] +
                    (ME_FILTER_TAP >> 1) +
                    ((ME_FILTER_TAP >> 1) *
                     context_ptr
                         ->interpolated_full_stride[list_index][ref_pic_index]),
                context_ptr
                    ->interpolated_full_stride[list_index][ref_pic_index],
                &(context_ptr->pos_b_buffer[list_index][ref_pic_index]
                                           [(ME_FILTER_TAP >> 1) *
                                            context_ptr->interpolated_stride]),
#endif
                &(context_ptr->pos_h_buffer[list_index][ref_pic_index][1]),
                &(context_ptr->pos_j_buffer[list_index][ref_pic_index][0]),
                x_search_area_origin,
                y_search_area_origin,
                inetger_mv);
        }
    }
}
#endif
static void quarter_pel_refinement_sb(
    MeContext
        *context_ptr,  //[IN/OUT]  ME context Ptr, used to get/update ME results
    uint8_t *pos_full,     //[IN]
    uint32_t full_stride,  //[IN]
    uint8_t *pos_b,        //[IN]
    uint8_t *pos_h,        //[IN]
    uint8_t *pos_j,        //[IN]
    int16_t
        x_search_area_origin,  //[IN] search area origin in the horizontal
                               // direction, used to point to reference samples
    int16_t
        y_search_area_origin,  //[IN] search area origin in the vertical
                               // direction, used to point to reference samples
    uint32_t integer_mv);

/*******************************************
 * open_loop_me_quarter_pel_search_sblock
 *******************************************/
static void open_loop_me_quarter_pel_search_sblock(
    MeContext *context_ptr,
    uint32_t list_index, uint32_t ref_pic_index, int16_t x_search_area_origin,
    int16_t y_search_area_origin, uint32_t search_area_width,
    uint32_t search_area_height)
{
    uint32_t search_index_x, search_index_y;
    for (search_index_y = 0; search_index_y < search_area_height;
         search_index_y++) {
        for (search_index_x = 0; search_index_x < search_area_width;
             search_index_x++) {
            int32_t mvx = (int32_t)search_index_x + x_search_area_origin;
            int32_t mvy = (int32_t)search_index_y + y_search_area_origin;
            uint32_t mv1 = (((uint16_t)mvy) << 18);
            uint16_t mv2 = (((uint16_t)mvx << 2));
            uint32_t mv0 = mv1 | mv2;
            int16_t x_mv = _MVXT(mv0);
            int16_t y_mv = _MVYT(mv0);
            uint32_t inetger_mv = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
            quarter_pel_refinement_sb(
                context_ptr,
#if M0_HIGH_PRECISION_INTERPOLATION
                context_ptr->integer_buffer_ptr[listIndex][ref_pic_index] +
                    (ME_FILTER_PAD_DISTANCE >> 1) +
                    ((ME_FILTER_PAD_DISTANCE >> 1) *
                     context_ptr
                         ->interpolated_full_stride[listIndex][ref_pic_index]),
                context_ptr->interpolated_full_stride[listIndex][ref_pic_index],
                &(context_ptr->pos_b_buffer
                      [listIndex][ref_pic_index]
                      [(ME_FILTER_PAD_DISTANCE >> 1) *
                       context_ptr->interpolated_stride]),  // points to b
                                                            // position of the
                                                            // figure above
#else
                context_ptr->integer_buffer_ptr[list_index][ref_pic_index] +
                    (ME_FILTER_TAP >> 1) +
                    ((ME_FILTER_TAP >> 1) *
                     context_ptr
                         ->interpolated_full_stride[list_index][ref_pic_index]),
                context_ptr
                    ->interpolated_full_stride[list_index][ref_pic_index],
                &(context_ptr->pos_b_buffer
                      [list_index][ref_pic_index]
                      [(ME_FILTER_TAP >> 1) *
                       context_ptr->interpolated_stride]),  // points to b
                                                            // position of the
                                                            // figure above
#endif
                &(context_ptr->pos_h_buffer[list_index][ref_pic_index]
                                           [1]),  // points to h position of the
                                                  // figure above
                &(context_ptr->pos_j_buffer[list_index][ref_pic_index]
                                           [0]),  // points to j position of the
                                                  // figure above
                x_search_area_origin,
                y_search_area_origin,
                inetger_mv);
        }
    }
}
/*******************************************
 * open_loop_me_fullpel_search_sblock
 *******************************************/
static void open_loop_me_fullpel_search_sblock(
    MeContext *context_ptr, uint32_t listIndex,
    uint32_t ref_pic_index,
    int16_t x_search_area_origin, int16_t y_search_area_origin,
    uint32_t search_area_width, uint32_t search_area_height, EbAsm asm_type) {
    uint32_t xSearchIndex, ySearchIndex;
    uint32_t searchAreaWidthRest8 = search_area_width & 7;
    uint32_t searchAreaWidthMult8 = search_area_width - searchAreaWidthRest8;

    for (ySearchIndex = 0; ySearchIndex < search_area_height; ySearchIndex++) {
        for (xSearchIndex = 0; xSearchIndex < searchAreaWidthMult8;
             xSearchIndex += 8) {
            // this function will do:  xSearchIndex, +1, +2, ..., +7
            open_loop_me_get_eight_search_point_results_block(
                context_ptr,
                listIndex,
                ref_pic_index,
                xSearchIndex +
                    ySearchIndex *
                        context_ptr->interpolated_full_stride[listIndex]
                                                             [ref_pic_index],
                (int32_t)xSearchIndex + x_search_area_origin,
                (int32_t)ySearchIndex + y_search_area_origin);
        }

        for (xSearchIndex = searchAreaWidthMult8;
             xSearchIndex < search_area_width;
             xSearchIndex++) {

            open_loop_me_get_search_point_results_block(
                context_ptr,
                listIndex,
                ref_pic_index,
                xSearchIndex +
                    ySearchIndex *
                        context_ptr->interpolated_full_stride[listIndex]
                                                             [ref_pic_index],
                (int32_t)xSearchIndex + x_search_area_origin,
                (int32_t)ySearchIndex + y_search_area_origin,
                asm_type);
        }
    }
}

#ifndef AVCCODEL
/*******************************************
 * HorizontalPelInterpolation
 *   interpolates the search region in the horizontal direction
 *******************************************/
static void HorizontalPelInterpolation(
    uint8_t *src,         // input parameter, input samples Ptr
    uint32_t src_stride,  // input parameter, input stride
    uint32_t width,       // input parameter, input area width
    uint32_t height,      // input parameter, input area height
    const int32_t
        *ifCoeff,  // input parameter, interpolation filter coefficients Ptr
    uint32_t inputBitDepth,  // input parameter, input sample bit depth
    uint32_t dst_stride,     // input parameter, output stride
    uint8_t *dst)            // output parameter, interpolated samples Ptr
{
    uint32_t x, y;
    const int32_t maxSampleValue = (1 << inputBitDepth) - 1;
    const int32_t ifOffset = 1 << (IFShift - 1);
    for (y = 0; y < height; ++y) {
        for (x = 0; x < width; ++x) {
            dst[x] = (uint8_t)CLIP3(
                0,
                (int32_t)maxSampleValue,
                ((((int32_t)src[x] + (int32_t)src[x + 3]) * ifCoeff[0] +
                  ((int32_t)src[x + 1] + (int32_t)src[x + 2]) * ifCoeff[1] +
                  ifOffset) >>
                 IFShift));
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
    uint8_t *src,         // input parameter, input samples ptr
    uint32_t src_stride,  // input parameter, input stride
    uint32_t width,       // input parameter, input area width
    uint32_t height,      // input parameter, input area height
    const int32_t
        ifCoeff[4],  // input parameter, interpolation filter coefficients Ptr
    uint32_t inputBitDepth,  // input parameter, input sample bit depth
    uint32_t dst_stride,     // input parameter, output stride
    uint8_t *dst)            // output parameter, interpolated samples Ptr
{
    uint32_t x, y;

    const int32_t maxSampleValue = (1 << inputBitDepth) - 1;
    const int32_t ifOffset = 1 << (IFShift - 1);

    const uint32_t srcStride2 = src_stride << 1;
    const uint32_t srcStride3 = srcStride2 + src_stride;

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = (uint8_t)CLIP3(
                0,
                maxSampleValue,
                ((((int32_t)src[x] + (int32_t)src[x + srcStride3]) *
                      ifCoeff[0] +
                  ((int32_t)src[x + src_stride] +
                   (int32_t)src[x + srcStride2]) *
                      ifCoeff[1] +
                  ifOffset) >>
                 IFShift));
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
static void AvcStyleInterpolation(
    uint8_t *srcOne,         // input parameter, input samples Ptr
    uint32_t srcOneStride,   // input parameter, input stride
    uint8_t *srcTwo,         // input parameter, input samples Ptr
    uint32_t srcTwoStride,   // input parameter, input stride
    uint32_t width,          // input parameter, input area width
    uint32_t height,         // input parameter, input area height
    uint32_t inputBitDepth,  // input parameter, input sample bit depth
    uint32_t dst_stride,     // input parameter, output stride
    uint8_t *dst)            // output parameter, interpolated samples Ptr
{
    uint32_t x, y;
    int32_t maxSampleValue = POW2(inputBitDepth) - 1;

    for (y = 0; y < height; ++y) {
        for (x = 0; x < width; ++x) {
            dst[x] =
                (uint8_t)CLIP3(0,
                               (int32_t)maxSampleValue,
                               (((int32_t)srcOne[x] + (int32_t)srcTwo[x] + 1) >>
                                IFShiftAvcStyle));
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
void InterpolateSearchRegionAVC(
    MeContext *context_ptr,  // input/output parameter, ME context ptr, used to
                             // get/set interpolated search area Ptr
    uint32_t listIndex,      // Refrence picture list index
    uint32_t ref_pic_index,
    uint8_t *searchRegionBuffer,  // input parameter, search region index, used
                                  // to point to reference samples
    uint32_t lumaStride,          // input parameter, reference Picture stride
    uint32_t search_area_width,   // input parameter, search area width
    uint32_t search_area_height,  // input parameter, search area height
    uint32_t inputBitDepth)       // input parameter, input sample bit depth
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
    uint32_t searchAreaWidthForAsm = ROUND_UP_MUL_8(search_area_width + 2);

#ifdef AVCCODEL

    (void)inputBitDepth;
    // Half pel interpolation of the search region using f1 -> pos_b_buffer
    if (searchAreaWidthForAsm) {
        avc_style_luma_interpolation_filter(
            searchRegionBuffer - (ME_FILTER_TAP >> 1) * lumaStride -
                (ME_FILTER_TAP >> 1) + 1,
            lumaStride,
            context_ptr->pos_b_buffer[listIndex][ref_pic_index],
            context_ptr->interpolated_stride,
            searchAreaWidthForAsm,
            search_area_height + ME_FILTER_TAP,
            context_ptr->avctemp_buffer,
            EB_FALSE,
            2,
            2);
    }

    // Half pel interpolation of the search region using f1 -> pos_h_buffer
    if (searchAreaWidthForAsm) {
        avc_style_luma_interpolation_filter(
            searchRegionBuffer - (ME_FILTER_TAP >> 1) * lumaStride - 1 +
                lumaStride,
            lumaStride,
            context_ptr->pos_h_buffer[listIndex][ref_pic_index],
            context_ptr->interpolated_stride,
            searchAreaWidthForAsm,
            search_area_height + 1,
            context_ptr->avctemp_buffer,
            EB_FALSE,
            2,
            8);
    }

    if (searchAreaWidthForAsm) {
        // Half pel interpolation of the search region using f1 -> pos_j_buffer
        avc_style_luma_interpolation_filter(
            context_ptr->pos_b_buffer[listIndex][ref_pic_index] +
                context_ptr->interpolated_stride,
            context_ptr->interpolated_stride,
            context_ptr->pos_j_buffer[listIndex][ref_pic_index],
            context_ptr->interpolated_stride,
            searchAreaWidthForAsm,
            search_area_height + 1,
            context_ptr->avctemp_buffer,
            EB_FALSE,
            2,
            8);
    }

#else

    // Half pel interpolation of the search region using f1 -> pos_b_buffer
    HorizontalPelInterpolation(searchRegionBuffer -
                                   (ME_FILTER_TAP >> 1) * lumaStride -
                                   (ME_FILTER_TAP >> 1),
                               lumaStride,
                               search_area_width + 1,
                               search_area_height + ME_FILTER_TAP,
                               &(me_if_coeff[F1][0]),
                               inputBitDepth,
                               context_ptr->interpolated_stride,
                               context_ptr->pos_b_buffer);

    // Half pel interpolation of the search region using f1 -> pos_h_buffer
    VerticalPelInterpolation(
        searchRegionBuffer - (ME_FILTER_TAP >> 1) * lumaStride - 1,
        lumaStride,
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
 * InterpolateSearchRegion AVC
 *   interpolates the search area
 *   the whole search area is interpolated 15 times
 *   for each sub position an interpolation is done
 *   15 buffers are required for the storage of the interpolated samples.
 *   F0: {-4, 54, 16, -2}
 *   F1: {-4, 36, 36, -4}
 *   F2: {-2, 16, 54, -4}
 ********************************************/
void interpolate_search_region_AVC_chroma(
    MeContext *context_ptr,  // input/output parameter, ME context ptr, used to
                             // get/set interpolated search area Ptr
    uint8_t *search_region_buffer_cb,  // input parameter, search region buffer
                                       // cb, used to point to reference samples
    uint8_t *search_region_buffer_cr,  // input parameter, search region buffer
                                       // cr, used to point to reference samples
    uint8_t **pos_b_buffer_ch, uint8_t **pos_h_buffer_ch,
    uint8_t **pos_j_buffer_ch, uint32_t interpolated_stride_ch,
    uint32_t interpolated_full_stride_ch,  // input parameter, reference Picture
                                           // stride
    uint32_t search_area_width,            // input parameter, search area width
    uint32_t search_area_height,  // input parameter, search area height
    uint32_t input_bit_depth)     // input parameter, input sample bit depth
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
    uint32_t searchAreaWidthForAsm = ROUND_UP_MUL_8(search_area_width + 2);

    (void)input_bit_depth;
    // Half pel interpolation of the search region using f1 -> pos_b_buffer
    if (searchAreaWidthForAsm) {
        // Cb
        avc_style_luma_interpolation_filter(
            search_region_buffer_cb -
                (ME_FILTER_TAP >> 1) * interpolated_full_stride_ch -
                (ME_FILTER_TAP >> 1) + 1,
            interpolated_full_stride_ch,
            pos_b_buffer_ch[0],
            interpolated_stride_ch,
            searchAreaWidthForAsm,
            search_area_height + ME_FILTER_TAP,
            context_ptr->avctemp_buffer,
            EB_FALSE,
            2,
            2);
        // Cr
        avc_style_luma_interpolation_filter(
            search_region_buffer_cr -
                (ME_FILTER_TAP >> 1) * interpolated_full_stride_ch -
                (ME_FILTER_TAP >> 1) + 1,
            interpolated_full_stride_ch,
            pos_b_buffer_ch[1],
            interpolated_stride_ch,
            searchAreaWidthForAsm,
            search_area_height + ME_FILTER_TAP,
            context_ptr->avctemp_buffer,
            EB_FALSE,
            2,
            2);
    }

    // Half pel interpolation of the search region using f1 -> pos_h_buffer
    if (searchAreaWidthForAsm) {
        // Cb
        avc_style_luma_interpolation_filter(
            search_region_buffer_cb -
                (ME_FILTER_TAP >> 1) * interpolated_full_stride_ch - 1 +
                interpolated_full_stride_ch,
            interpolated_full_stride_ch,
            pos_h_buffer_ch[0],
            interpolated_stride_ch,
            searchAreaWidthForAsm,
            search_area_height + 1,
            context_ptr->avctemp_buffer,
            EB_FALSE,
            2,
            8);
        // Cr
        avc_style_luma_interpolation_filter(
            search_region_buffer_cr -
                (ME_FILTER_TAP >> 1) * interpolated_full_stride_ch - 1 +
                interpolated_full_stride_ch,
            interpolated_full_stride_ch,
            pos_h_buffer_ch[1],
            interpolated_stride_ch,
            searchAreaWidthForAsm,
            search_area_height + 1,
            context_ptr->avctemp_buffer,
            EB_FALSE,
            2,
            8);
    }

    // Half pel interpolation of the search region using f1 -> pos_j_buffer
    if (searchAreaWidthForAsm) {
        // Cb
        avc_style_luma_interpolation_filter(
            pos_b_buffer_ch[0] + interpolated_stride_ch,
            interpolated_stride_ch,
            pos_j_buffer_ch[0],
            interpolated_stride_ch,
            searchAreaWidthForAsm,
            search_area_height + 1,
            context_ptr->avctemp_buffer,
            EB_FALSE,
            2,
            8);
        // Cr
        avc_style_luma_interpolation_filter(
            pos_b_buffer_ch[1] + interpolated_stride_ch,
            interpolated_stride_ch,
            pos_j_buffer_ch[1],
            interpolated_stride_ch,
            searchAreaWidthForAsm,
            search_area_height + 1,
            context_ptr->avctemp_buffer,
            EB_FALSE,
            2,
            8);
    }
}

/*******************************************
 * PU_HalfPelRefinement
 *   performs Half Pel refinement for one PU
 *******************************************/
static void PU_HalfPelRefinement(
    SequenceControlSet
        *sequence_control_set_ptr,  // input parameter, Sequence control set Ptr
    MeContext
        *context_ptr,  // input parameter, ME context Ptr, used to get SB Ptr
    uint8_t *refBuffer, uint32_t ref_stride, uint32_t *pBestSsd,
    uint32_t puLcuBufferIndex,  // input parameter, PU origin, used to point to
                                // source samples
    uint8_t *pos_b_buffer,  // input parameter, position "b" interpolated search
                            // area Ptr
    uint8_t *pos_h_buffer,  // input parameter, position "h" interpolated search
                            // area Ptr
    uint8_t *pos_j_buffer,  // input parameter, position "j" interpolated search
                            // area Ptr
    uint32_t pu_width,      // input parameter, PU width
    uint32_t pu_height,     // input parameter, PU height
    int16_t x_search_area_origin,  // input parameter, search area origin in the
                                   // horizontal direction, used to point to
                                   // reference samples
    int16_t y_search_area_origin,  // input parameter, search area origin in the
                                   // vertical direction, used to point to
                                   // reference samples
    uint32_t *pBestSad, uint32_t *pBestMV,
    uint8_t *psubPelDirection)
{
    EncodeContext *encode_context_ptr =
        sequence_control_set_ptr->encode_context_ptr;

    int32_t searchRegionIndex;
    uint64_t bestHalfSad = 0;
    uint64_t distortionLeftPosition = 0;
    uint64_t distortionRightPosition = 0;
    uint64_t distortionTopPosition = 0;
    uint64_t distortionBottomPosition = 0;
    uint64_t distortionTopLeftPosition = 0;
    uint64_t distortionTopRightPosition = 0;
    uint64_t distortionBottomLeftPosition = 0;
    uint64_t distortionBottomRightPosition = 0;

    int16_t xMvHalf[8];
    int16_t yMvHalf[8];

    int16_t x_mv = _MVXT(*pBestMV);
    int16_t y_mv = _MVYT(*pBestMV);
    int16_t xSearchIndex = (x_mv >> 2) - x_search_area_origin;
    int16_t ySearchIndex = (y_mv >> 2) - y_search_area_origin;

    (void)sequence_control_set_ptr;
    (void)encode_context_ptr;

    // TODO : remove these, and update the MV by just shifts

    xMvHalf[0] = x_mv - 2;  // L  position
    xMvHalf[1] = x_mv + 2;  // R  position
    xMvHalf[2] = x_mv;      // T  position
    xMvHalf[3] = x_mv;      // B  position
    xMvHalf[4] = x_mv - 2;  // TL position
    xMvHalf[5] = x_mv + 2;  // TR position
    xMvHalf[6] = x_mv + 2;  // BR position
    xMvHalf[7] = x_mv - 2;  // BL position

    yMvHalf[0] = y_mv;      // L  position
    yMvHalf[1] = y_mv;      // R  position
    yMvHalf[2] = y_mv - 2;  // T  position
    yMvHalf[3] = y_mv + 2;  // B  position
    yMvHalf[4] = y_mv - 2;  // TL position
    yMvHalf[5] = y_mv - 2;  // TR position
    yMvHalf[6] = y_mv + 2;  // BR position
    yMvHalf[7] = y_mv + 2;  // BL position

    // Compute SSD for the best full search candidate
    if (context_ptr->fractional_search_method == SSD_SEARCH) {
        *pBestSsd = (uint32_t)spatial_full_distortion_kernel(
                context_ptr->sb_src_ptr,
                puLcuBufferIndex,
                context_ptr->sb_src_stride,
                refBuffer,
                ySearchIndex * ref_stride + xSearchIndex,
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
        searchRegionIndex =
            xSearchIndex +
            (int16_t)context_ptr->interpolated_stride * ySearchIndex;
        distortionLeftPosition =
            (context_ptr->fractional_search_method == SSD_SEARCH)
                ? spatial_full_distortion_kernel(
                          context_ptr->sb_src_ptr,
                          puLcuBufferIndex,
                          context_ptr->sb_src_stride,
                          pos_b_buffer,
                          searchRegionIndex,
                          context_ptr->interpolated_stride,
                          pu_width,
                          pu_height)
                : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                      ? (nxm_sad_kernel(
                            &(context_ptr->sb_src_ptr[puLcuBufferIndex]),
                            context_ptr->sb_src_stride << 1,
                            &(pos_b_buffer[searchRegionIndex]),
                            context_ptr->interpolated_stride << 1,
                            pu_height >> 1,
                            pu_width))
                            << 1
                      : nxm_sad_kernel(
                            &(context_ptr->sb_src_ptr[puLcuBufferIndex]),
                            context_ptr->sb_src_stride,
                            &(pos_b_buffer[searchRegionIndex]),
                            context_ptr->interpolated_stride,
                            pu_height,
                            pu_width);
        if (context_ptr->fractional_search_method == SSD_SEARCH) {
            if (distortionLeftPosition < *pBestSsd) {
                *pBestSad = (uint32_t)
                    nxm_sad_kernel(
                        &(context_ptr->sb_src_ptr[puLcuBufferIndex]),
                        context_ptr->sb_src_stride,
                        &(pos_b_buffer[searchRegionIndex]),
                        context_ptr->interpolated_stride,
                        pu_height,
                        pu_width);
                *pBestMV =
                    ((uint16_t)yMvHalf[0] << 16) | ((uint16_t)xMvHalf[0]);
                *pBestSsd = (uint32_t)distortionLeftPosition;
            }
        } else {
            if (distortionLeftPosition < *pBestSad) {
                *pBestSad = (uint32_t)distortionLeftPosition;
                *pBestMV =
                    ((uint16_t)yMvHalf[0] << 16) | ((uint16_t)xMvHalf[0]);
            }
        }
        // R position
        searchRegionIndex++;
        distortionRightPosition =
            (context_ptr->fractional_search_method == SSD_SEARCH)
                ? spatial_full_distortion_kernel(
                          context_ptr->sb_src_ptr,
                          puLcuBufferIndex,
                          context_ptr->sb_src_stride,
                          pos_b_buffer,
                          searchRegionIndex,
                          context_ptr->interpolated_stride,
                          pu_width,
                          pu_height)
                : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                      ? (nxm_sad_kernel(
                            &(context_ptr->sb_src_ptr[puLcuBufferIndex]),
                            context_ptr->sb_src_stride << 1,
                            &(pos_b_buffer[searchRegionIndex]),
                            context_ptr->interpolated_stride << 1,
                            pu_height >> 1,
                            pu_width))
                            << 1
                      : nxm_sad_kernel(
                            &(context_ptr->sb_src_ptr[puLcuBufferIndex]),
                            context_ptr->sb_src_stride,
                            &(pos_b_buffer[searchRegionIndex]),
                            context_ptr->interpolated_stride,
                            pu_height,
                            pu_width);
        if (context_ptr->fractional_search_method == SSD_SEARCH) {
            if (distortionRightPosition < *pBestSsd) {
                *pBestSad = (uint32_t)
                    nxm_sad_kernel(
                        &(context_ptr->sb_src_ptr[puLcuBufferIndex]),
                        context_ptr->sb_src_stride,
                        &(pos_b_buffer[searchRegionIndex]),
                        context_ptr->interpolated_stride,
                        pu_height,
                        pu_width);
                *pBestMV =
                    ((uint16_t)yMvHalf[1] << 16) | ((uint16_t)xMvHalf[1]);
                *pBestSsd = (uint32_t)distortionRightPosition;
            }
        } else {
            if (distortionRightPosition < *pBestSad) {
                *pBestSad = (uint32_t)distortionRightPosition;
                *pBestMV =
                    ((uint16_t)yMvHalf[1] << 16) | ((uint16_t)xMvHalf[1]);
            }
        }
        // T position
        searchRegionIndex =
            xSearchIndex +
            (int16_t)context_ptr->interpolated_stride * ySearchIndex;
        distortionTopPosition =
            (context_ptr->fractional_search_method == SSD_SEARCH)
                ? spatial_full_distortion_kernel(
                          context_ptr->sb_src_ptr,
                          puLcuBufferIndex,
                          context_ptr->sb_src_stride,
                          pos_h_buffer,
                          searchRegionIndex,
                          context_ptr->interpolated_stride,
                          pu_width,
                          pu_height)
                : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                      ? (nxm_sad_kernel(
                            &(context_ptr->sb_src_ptr[puLcuBufferIndex]),
                            context_ptr->sb_src_stride << 1,
                            &(pos_h_buffer[searchRegionIndex]),
                            context_ptr->interpolated_stride << 1,
                            pu_height >> 1,
                            pu_width))
                            << 1
                      : nxm_sad_kernel(
                            &(context_ptr->sb_src_ptr[puLcuBufferIndex]),
                            context_ptr->sb_src_stride,
                            &(pos_h_buffer[searchRegionIndex]),
                            context_ptr->interpolated_stride,
                            pu_height,
                            pu_width);
        if (context_ptr->fractional_search_method == SSD_SEARCH) {
            if (distortionTopPosition < *pBestSsd) {
                *pBestSad = (uint32_t)
                    nxm_sad_kernel(
                        &(context_ptr->sb_src_ptr[puLcuBufferIndex]),
                        context_ptr->sb_src_stride,
                        &(pos_h_buffer[searchRegionIndex]),
                        context_ptr->interpolated_stride,
                        pu_height,
                        pu_width);
                *pBestMV =
                    ((uint16_t)yMvHalf[2] << 16) | ((uint16_t)xMvHalf[2]);
                *pBestSsd = (uint32_t)distortionTopPosition;
            }
        } else {
            if (distortionTopPosition < *pBestSad) {
                *pBestSad = (uint32_t)distortionTopPosition;
                *pBestMV =
                    ((uint16_t)yMvHalf[2] << 16) | ((uint16_t)xMvHalf[2]);
            }
        }

        // B position
        searchRegionIndex += (int16_t)context_ptr->interpolated_stride;
        distortionBottomPosition =
            (context_ptr->fractional_search_method == SSD_SEARCH)
                ? spatial_full_distortion_kernel(
                          context_ptr->sb_src_ptr,
                          puLcuBufferIndex,
                          context_ptr->sb_src_stride,
                          pos_h_buffer,
                          searchRegionIndex,
                          context_ptr->interpolated_stride,
                          pu_width,
                          pu_height)
                : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                      ? (nxm_sad_kernel(
                            &(context_ptr->sb_src_ptr[puLcuBufferIndex]),
                            context_ptr->sb_src_stride << 1,
                            &(pos_h_buffer[searchRegionIndex]),
                            context_ptr->interpolated_stride << 1,
                            pu_height >> 1,
                            pu_width))
                            << 1
                      : nxm_sad_kernel(
                            &(context_ptr->sb_src_ptr[puLcuBufferIndex]),
                            context_ptr->sb_src_stride,
                            &(pos_h_buffer[searchRegionIndex]),
                            context_ptr->interpolated_stride,
                            pu_height,
                            pu_width);
        if (context_ptr->fractional_search_method == SSD_SEARCH) {
            if (distortionBottomPosition < *pBestSsd) {
                *pBestSad = (uint32_t)
                    nxm_sad_kernel(
                        &(context_ptr->sb_src_ptr[puLcuBufferIndex]),
                        context_ptr->sb_src_stride,
                        &(pos_h_buffer[searchRegionIndex]),
                        context_ptr->interpolated_stride,
                        pu_height,
                        pu_width);
                *pBestMV =
                    ((uint16_t)yMvHalf[3] << 16) | ((uint16_t)xMvHalf[3]);
                *pBestSsd = (uint32_t)distortionBottomPosition;
            }
        } else {
            if (distortionBottomPosition < *pBestSad) {
                *pBestSad = (uint32_t)distortionBottomPosition;
                *pBestMV =
                    ((uint16_t)yMvHalf[3] << 16) | ((uint16_t)xMvHalf[3]);
            }
        }

        // TL position
        searchRegionIndex =
            xSearchIndex +
            (int16_t)context_ptr->interpolated_stride * ySearchIndex;
        distortionTopLeftPosition =
            (context_ptr->fractional_search_method == SSD_SEARCH)
                ? spatial_full_distortion_kernel(
                          context_ptr->sb_src_ptr,
                          puLcuBufferIndex,
                          context_ptr->sb_src_stride,
                          pos_j_buffer,
                          searchRegionIndex,
                          context_ptr->interpolated_stride,
                          pu_width,
                          pu_height)
                : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                      ? (nxm_sad_kernel(
                            &(context_ptr->sb_src_ptr[puLcuBufferIndex]),
                            context_ptr->sb_src_stride << 1,
                            &(pos_j_buffer[searchRegionIndex]),
                            context_ptr->interpolated_stride << 1,
                            pu_height >> 1,
                            pu_width))
                            << 1
                      : nxm_sad_kernel(
                            &(context_ptr->sb_src_ptr[puLcuBufferIndex]),
                            context_ptr->sb_src_stride,
                            &(pos_j_buffer[searchRegionIndex]),
                            context_ptr->interpolated_stride,
                            pu_height,
                            pu_width);
        if (context_ptr->fractional_search_method == SSD_SEARCH) {
            if (distortionTopLeftPosition < *pBestSsd) {
                *pBestSad = (uint32_t)
                    nxm_sad_kernel(
                        &(context_ptr->sb_src_ptr[puLcuBufferIndex]),
                        context_ptr->sb_src_stride,
                        &(pos_j_buffer[searchRegionIndex]),
                        context_ptr->interpolated_stride,
                        pu_height,
                        pu_width);
                *pBestMV =
                    ((uint16_t)yMvHalf[4] << 16) | ((uint16_t)xMvHalf[4]);
                *pBestSsd = (uint32_t)distortionTopLeftPosition;
            }
        } else {
            if (distortionTopLeftPosition < *pBestSad) {
                *pBestSad = (uint32_t)distortionTopLeftPosition;
                *pBestMV =
                    ((uint16_t)yMvHalf[4] << 16) | ((uint16_t)xMvHalf[4]);
            }
        }

        // TR position
        searchRegionIndex++;
        distortionTopRightPosition =
            (context_ptr->fractional_search_method == SSD_SEARCH)
                ? spatial_full_distortion_kernel(
                          context_ptr->sb_src_ptr,
                          puLcuBufferIndex,
                          context_ptr->sb_src_stride,
                          pos_j_buffer,
                          searchRegionIndex,
                          context_ptr->interpolated_stride,
                          pu_width,
                          pu_height)
                : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                      ? (nxm_sad_kernel(
                            &(context_ptr->sb_src_ptr[puLcuBufferIndex]),
                            context_ptr->sb_src_stride << 1,
                            &(pos_j_buffer[searchRegionIndex]),
                            context_ptr->interpolated_stride << 1,
                            pu_height >> 1,
                            pu_width))
                            << 1
                      : nxm_sad_kernel(
                            &(context_ptr->sb_src_ptr[puLcuBufferIndex]),
                            context_ptr->sb_src_stride,
                            &(pos_j_buffer[searchRegionIndex]),
                            context_ptr->interpolated_stride,
                            pu_height,
                            pu_width);

        if (context_ptr->fractional_search_method == SSD_SEARCH) {
            if (distortionTopRightPosition < *pBestSsd) {
                *pBestSad = (uint32_t)
                    nxm_sad_kernel(
                        &(context_ptr->sb_src_ptr[puLcuBufferIndex]),
                        context_ptr->sb_src_stride,
                        &(pos_j_buffer[searchRegionIndex]),
                        context_ptr->interpolated_stride,
                        pu_height,
                        pu_width);
                *pBestMV =
                    ((uint16_t)yMvHalf[5] << 16) | ((uint16_t)xMvHalf[5]);
                *pBestSsd = (uint32_t)distortionTopRightPosition;
            }
        } else {
            if (distortionTopRightPosition < *pBestSad) {
                *pBestSad = (uint32_t)distortionTopRightPosition;
                *pBestMV =
                    ((uint16_t)yMvHalf[5] << 16) | ((uint16_t)xMvHalf[5]);
            }
        }

        // BR position
        searchRegionIndex += (int16_t)context_ptr->interpolated_stride;
        distortionBottomRightPosition =
            (context_ptr->fractional_search_method == SSD_SEARCH)
                ? spatial_full_distortion_kernel(
                          context_ptr->sb_src_ptr,
                          puLcuBufferIndex,
                          context_ptr->sb_src_stride,
                          pos_j_buffer,
                          searchRegionIndex,
                          context_ptr->interpolated_stride,
                          pu_width,
                          pu_height)
                : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                      ? (nxm_sad_kernel(
                            &(context_ptr->sb_src_ptr[puLcuBufferIndex]),
                            context_ptr->sb_src_stride << 1,
                            &(pos_j_buffer[searchRegionIndex]),
                            context_ptr->interpolated_stride << 1,
                            pu_height >> 1,
                            pu_width))
                            << 1
                      : nxm_sad_kernel(
                            &(context_ptr->sb_src_ptr[puLcuBufferIndex]),
                            context_ptr->sb_src_stride,
                            &(pos_j_buffer[searchRegionIndex]),
                            context_ptr->interpolated_stride,
                            pu_height,
                            pu_width);
        if (context_ptr->fractional_search_method == SSD_SEARCH) {
            if (distortionBottomRightPosition < *pBestSsd) {
                *pBestSad = (uint32_t)
                    nxm_sad_kernel(
                        &(context_ptr->sb_src_ptr[puLcuBufferIndex]),
                        context_ptr->sb_src_stride,
                        &(pos_j_buffer[searchRegionIndex]),
                        context_ptr->interpolated_stride,
                        pu_height,
                        pu_width);
                *pBestMV =
                    ((uint16_t)yMvHalf[6] << 16) | ((uint16_t)xMvHalf[6]);
                *pBestSsd = (uint32_t)distortionBottomRightPosition;
            }
        } else {
            if (distortionBottomRightPosition < *pBestSad) {
                *pBestSad = (uint32_t)distortionBottomRightPosition;
                *pBestMV =
                    ((uint16_t)yMvHalf[6] << 16) | ((uint16_t)xMvHalf[6]);
            }
        }

        // BL position
        searchRegionIndex--;
        distortionBottomLeftPosition =
            (context_ptr->fractional_search_method == SSD_SEARCH)
                ? spatial_full_distortion_kernel(
                          context_ptr->sb_src_ptr,
                          puLcuBufferIndex,
                          context_ptr->sb_src_stride,
                          pos_j_buffer,
                          searchRegionIndex,
                          context_ptr->interpolated_stride,
                          pu_width,
                          pu_height)
                : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                      ? (nxm_sad_kernel(
                            &(context_ptr->sb_src_ptr[puLcuBufferIndex]),
                            context_ptr->sb_src_stride << 1,
                            &(pos_j_buffer[searchRegionIndex]),
                            context_ptr->interpolated_stride << 1,
                            pu_height >> 1,
                            pu_width))
                            << 1
                      : (nxm_sad_kernel(
                            &(context_ptr->sb_src_ptr[puLcuBufferIndex]),
                            context_ptr->sb_src_stride,
                            &(pos_j_buffer[searchRegionIndex]),
                            context_ptr->interpolated_stride,
                            pu_height,
                            pu_width));

        if (context_ptr->fractional_search_method == SSD_SEARCH) {
            if (distortionBottomLeftPosition < *pBestSsd) {
                *pBestSad = (uint32_t)(
                    nxm_sad_kernel(
                        &(context_ptr->sb_src_ptr[puLcuBufferIndex]),
                        context_ptr->sb_src_stride,
                        &(pos_j_buffer[searchRegionIndex]),
                        context_ptr->interpolated_stride,
                        pu_height,
                        pu_width));
                *pBestMV =
                    ((uint16_t)yMvHalf[7] << 16) | ((uint16_t)xMvHalf[7]);
                *pBestSsd = (uint32_t)distortionBottomLeftPosition;
            }
        } else {
            if (distortionBottomLeftPosition < *pBestSad) {
                *pBestSad = (uint32_t)distortionBottomLeftPosition;
                *pBestMV =
                    ((uint16_t)yMvHalf[7] << 16) | ((uint16_t)xMvHalf[7]);
            }
        }
    }

    bestHalfSad =
        MIN(distortionLeftPosition,
            MIN(distortionRightPosition,
                MIN(distortionTopPosition,
                    MIN(distortionBottomPosition,
                        MIN(distortionTopLeftPosition,
                            MIN(distortionTopRightPosition,
                                MIN(distortionBottomLeftPosition,
                                    distortionBottomRightPosition)))))));

    if (bestHalfSad == distortionLeftPosition)
        *psubPelDirection = LEFT_POSITION;
    else if (bestHalfSad == distortionRightPosition)
        *psubPelDirection = RIGHT_POSITION;
    else if (bestHalfSad == distortionTopPosition)
        *psubPelDirection = TOP_POSITION;
    else if (bestHalfSad == distortionBottomPosition)
        *psubPelDirection = BOTTOM_POSITION;
    else if (bestHalfSad == distortionTopLeftPosition)
        *psubPelDirection = TOP_LEFT_POSITION;
    else if (bestHalfSad == distortionTopRightPosition)
        *psubPelDirection = TOP_RIGHT_POSITION;
    else if (bestHalfSad == distortionBottomLeftPosition)
        *psubPelDirection = BOTTOM_LEFT_POSITION;
    else if (bestHalfSad == distortionBottomRightPosition)
        *psubPelDirection = BOTTOM_RIGHT_POSITION;
    return;
}

/*******************************************
 * HalfPelSearch_LCU
 *   performs Half Pel refinement for the 85 PUs
 *******************************************/
void HalfPelSearch_LCU(
    SequenceControlSet
        *sequence_control_set_ptr,  // input parameter, Sequence control set Ptr
    PictureParentControlSet *picture_control_set_ptr,
    MeContext *context_ptr,  // input/output parameter, ME context Ptr, used to
                             // get/update ME results
    uint8_t *refBuffer, uint32_t ref_stride,
    uint8_t *pos_b_buffer,  // input parameter, position "b" interpolated search
                            // area Ptr
    uint8_t *pos_h_buffer,  // input parameter, position "h" interpolated search
                            // area Ptr
    uint8_t *pos_j_buffer,  // input parameter, position "j" interpolated search
                            // area Ptr
    int16_t x_search_area_origin,  // input parameter, search area origin in the
                                   // horizontal direction, used to point to
                                   // reference samples
    int16_t y_search_area_origin,  // input parameter, search area origin in the
                                   // vertical direction, used to point to
                                   // reference samples
    EbBool disable8x8CuInMeFlag, EbBool enableHalfPel32x32,
    EbBool enableHalfPel16x16, EbBool enableHalfPel8x8)
{
    uint32_t idx;
    uint32_t pu_index;
    uint32_t puShiftXIndex;
    uint32_t puShiftYIndex;
    uint32_t puLcuBufferIndex;
    uint32_t posbBufferIndex;
    uint32_t poshBufferIndex;
    uint32_t posjBufferIndex;

    if (context_ptr->fractional_search64x64)
        PU_HalfPelRefinement(sequence_control_set_ptr,
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
                             context_ptr->p_best_sad64x64,
                             context_ptr->p_best_mv64x64,
                             &context_ptr->psub_pel_direction64x64);

    if (enableHalfPel32x32) {
        // 32x32 [4 partitions]
        for (pu_index = 0; pu_index < 4; ++pu_index) {
            puShiftXIndex = (pu_index & 0x01) << 5;
            puShiftYIndex = (pu_index >> 1) << 5;

            puLcuBufferIndex =
                puShiftXIndex + puShiftYIndex * context_ptr->sb_src_stride;
            posbBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;
            poshBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;
            posjBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;

            PU_HalfPelRefinement(
                sequence_control_set_ptr,
                context_ptr,
                &(refBuffer[puShiftYIndex * ref_stride + puShiftXIndex]),
                ref_stride,
                &context_ptr->p_best_ssd32x32[pu_index],
                puLcuBufferIndex,
                &(pos_b_buffer[posbBufferIndex]),
                &(pos_h_buffer[poshBufferIndex]),
                &(pos_j_buffer[posjBufferIndex]),
                32,
                32,
                x_search_area_origin,
                y_search_area_origin,
                &context_ptr->p_best_sad32x32[pu_index],
                &context_ptr->p_best_mv32x32[pu_index],
                &context_ptr->psub_pel_direction32x32[pu_index]);
        }
    }
    if (enableHalfPel16x16) {
        // 16x16 [16 partitions]
        for (pu_index = 0; pu_index < 16; ++pu_index) {
            idx = tab16x16[pu_index];

            puShiftXIndex = (pu_index & 0x03) << 4;
            puShiftYIndex = (pu_index >> 2) << 4;

            puLcuBufferIndex =
                puShiftXIndex + puShiftYIndex * context_ptr->sb_src_stride;
            posbBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;
            poshBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;
            posjBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;

            PU_HalfPelRefinement(
                sequence_control_set_ptr,
                context_ptr,
                &(refBuffer[puShiftYIndex * ref_stride + puShiftXIndex]),
                ref_stride,
                &context_ptr->p_best_ssd16x16[idx],
                puLcuBufferIndex,
                &(pos_b_buffer[posbBufferIndex]),
                &(pos_h_buffer[poshBufferIndex]),
                &(pos_j_buffer[posjBufferIndex]),
                16,
                16,
                x_search_area_origin,
                y_search_area_origin,
                &context_ptr->p_best_sad16x16[idx],
                &context_ptr->p_best_mv16x16[idx],
                &context_ptr->psub_pel_direction16x16[idx]);
        }
    }
    if (enableHalfPel8x8) {
        // 8x8   [64 partitions]
        if (!disable8x8CuInMeFlag) {
            for (pu_index = 0; pu_index < 64; ++pu_index) {
                idx = tab8x8[pu_index];  // TODO bitwise this

                puShiftXIndex = (pu_index & 0x07) << 3;
                puShiftYIndex = (pu_index >> 3) << 3;

                puLcuBufferIndex =
                    puShiftXIndex + puShiftYIndex * context_ptr->sb_src_stride;

                posbBufferIndex =
                    puShiftXIndex +
                    puShiftYIndex * context_ptr->interpolated_stride;
                poshBufferIndex =
                    puShiftXIndex +
                    puShiftYIndex * context_ptr->interpolated_stride;
                posjBufferIndex =
                    puShiftXIndex +
                    puShiftYIndex * context_ptr->interpolated_stride;

                PU_HalfPelRefinement(
                    sequence_control_set_ptr,
                    context_ptr,
                    &(refBuffer[puShiftYIndex * ref_stride + puShiftXIndex]),
                    ref_stride,
                    &context_ptr->p_best_ssd8x8[idx],
                    puLcuBufferIndex,
                    &(pos_b_buffer[posbBufferIndex]),
                    &(pos_h_buffer[poshBufferIndex]),
                    &(pos_j_buffer[posjBufferIndex]),
                    8,
                    8,
                    x_search_area_origin,
                    y_search_area_origin,
                    &context_ptr->p_best_sad8x8[idx],
                    &context_ptr->p_best_mv8x8[idx],
                    &context_ptr->psub_pel_direction8x8[idx]);
            }
        }
    }
    if (picture_control_set_ptr->pic_depth_mode <= PIC_ALL_C_DEPTH_MODE) {
        // 64x32
        for (pu_index = 0; pu_index < 2; ++pu_index) {
            puShiftXIndex = 0;
            puShiftYIndex = pu_index << 5;

            puLcuBufferIndex =
                puShiftXIndex + puShiftYIndex * context_ptr->sb_src_stride;

            posbBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;
            poshBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;
            posjBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;

            PU_HalfPelRefinement(
                sequence_control_set_ptr,
                context_ptr,
                &(refBuffer[puShiftYIndex * ref_stride + puShiftXIndex]),
                ref_stride,
                &context_ptr->p_best_ssd64x32[pu_index],
                puLcuBufferIndex,
                &(pos_b_buffer[posbBufferIndex]),
                &(pos_h_buffer[poshBufferIndex]),
                &(pos_j_buffer[posjBufferIndex]),
                64,
                32,
                x_search_area_origin,
                y_search_area_origin,
                &context_ptr->p_best_sad64x32[pu_index],
                &context_ptr->p_best_mv64x32[pu_index],
                &context_ptr->psub_pel_direction64x32[pu_index]);
        }

        // 32x16
        for (pu_index = 0; pu_index < 8; ++pu_index) {
            idx = tab32x16[pu_index];  // TODO bitwise this

            puShiftXIndex = (pu_index & 0x01) << 5;
            puShiftYIndex = (pu_index >> 1) << 4;

            puLcuBufferIndex =
                puShiftXIndex + puShiftYIndex * context_ptr->sb_src_stride;

            posbBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;
            poshBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;
            posjBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;

            PU_HalfPelRefinement(
                sequence_control_set_ptr,
                context_ptr,
                &(refBuffer[puShiftYIndex * ref_stride + puShiftXIndex]),
                ref_stride,
                &context_ptr->p_best_ssd32x16[idx],
                puLcuBufferIndex,
                &(pos_b_buffer[posbBufferIndex]),
                &(pos_h_buffer[poshBufferIndex]),
                &(pos_j_buffer[posjBufferIndex]),
                32,
                16,
                x_search_area_origin,
                y_search_area_origin,
                &context_ptr->p_best_sad32x16[idx],
                &context_ptr->p_best_mv32x16[idx],
                &context_ptr->psub_pel_direction32x16[idx]);
        }

        // 16x8
        for (pu_index = 0; pu_index < 32; ++pu_index) {
            idx = tab16x8[pu_index];

            puShiftXIndex = (pu_index & 0x03) << 4;
            puShiftYIndex = (pu_index >> 2) << 3;

            puLcuBufferIndex =
                puShiftXIndex + puShiftYIndex * context_ptr->sb_src_stride;

            posbBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;
            poshBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;
            posjBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;

            PU_HalfPelRefinement(
                sequence_control_set_ptr,
                context_ptr,
                &(refBuffer[puShiftYIndex * ref_stride + puShiftXIndex]),
                ref_stride,
                &context_ptr->p_best_ssd16x8[idx],
                puLcuBufferIndex,
                &(pos_b_buffer[posbBufferIndex]),
                &(pos_h_buffer[poshBufferIndex]),
                &(pos_j_buffer[posjBufferIndex]),
                16,
                8,
                x_search_area_origin,
                y_search_area_origin,
                &context_ptr->p_best_sad16x8[idx],
                &context_ptr->p_best_mv16x8[idx],
                &context_ptr->psub_pel_direction16x8[idx]);
        }

        // 32x64
        for (pu_index = 0; pu_index < 2; ++pu_index) {
            puShiftXIndex = pu_index << 5;
            puShiftYIndex = 0;

            puLcuBufferIndex =
                puShiftXIndex + puShiftYIndex * context_ptr->sb_src_stride;

            posbBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;
            poshBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;
            posjBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;

            PU_HalfPelRefinement(
                sequence_control_set_ptr,
                context_ptr,
                &(refBuffer[puShiftYIndex * ref_stride + puShiftXIndex]),
                ref_stride,
                &context_ptr->p_best_ssd32x64[pu_index],
                puLcuBufferIndex,
                &(pos_b_buffer[posbBufferIndex]),
                &(pos_h_buffer[poshBufferIndex]),
                &(pos_j_buffer[posjBufferIndex]),
                32,
                64,
                x_search_area_origin,
                y_search_area_origin,
                &context_ptr->p_best_sad32x64[pu_index],
                &context_ptr->p_best_mv32x64[pu_index],
                &context_ptr->psub_pel_direction32x64[pu_index]);
        }

        // 16x32
        for (pu_index = 0; pu_index < 8; ++pu_index) {
            idx = tab16x32[pu_index];

            puShiftXIndex = (pu_index & 0x03) << 4;
            puShiftYIndex = (pu_index >> 2) << 5;

            puLcuBufferIndex =
                puShiftXIndex + puShiftYIndex * context_ptr->sb_src_stride;

            posbBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;
            poshBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;
            posjBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;

            PU_HalfPelRefinement(
                sequence_control_set_ptr,
                context_ptr,
                &(refBuffer[puShiftYIndex * ref_stride + puShiftXIndex]),
                ref_stride,
                &context_ptr->p_best_ssd16x32[idx],
                puLcuBufferIndex,
                &(pos_b_buffer[posbBufferIndex]),
                &(pos_h_buffer[poshBufferIndex]),
                &(pos_j_buffer[posjBufferIndex]),
                16,
                32,
                x_search_area_origin,
                y_search_area_origin,
                &context_ptr->p_best_sad16x32[idx],
                &context_ptr->p_best_mv16x32[idx],
                &context_ptr->psub_pel_direction16x32[idx]);
        }

        // 8x16
        for (pu_index = 0; pu_index < 32; ++pu_index) {
            idx = tab8x16[pu_index];

            puShiftXIndex = (pu_index & 0x07) << 3;
            puShiftYIndex = (pu_index >> 3) << 4;

            puLcuBufferIndex =
                puShiftXIndex + puShiftYIndex * context_ptr->sb_src_stride;

            posbBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;
            poshBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;
            posjBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;

            PU_HalfPelRefinement(
                sequence_control_set_ptr,
                context_ptr,
                &(refBuffer[puShiftYIndex * ref_stride + puShiftXIndex]),
                ref_stride,
                &context_ptr->p_best_ssd8x16[idx],
                puLcuBufferIndex,
                &(pos_b_buffer[posbBufferIndex]),
                &(pos_h_buffer[poshBufferIndex]),
                &(pos_j_buffer[posjBufferIndex]),
                8,
                16,
                x_search_area_origin,
                y_search_area_origin,
                &context_ptr->p_best_sad8x16[idx],
                &context_ptr->p_best_mv8x16[idx],
                &context_ptr->psub_pel_direction8x16[idx]);
        }

        // 32x8
        for (pu_index = 0; pu_index < 16; ++pu_index) {
            idx = tab32x8[pu_index];

            puShiftXIndex = (pu_index & 0x01) << 5;
            puShiftYIndex = (pu_index >> 1) << 3;

            puLcuBufferIndex =
                puShiftXIndex + puShiftYIndex * context_ptr->sb_src_stride;

            posbBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;
            poshBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;
            posjBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;

            PU_HalfPelRefinement(
                sequence_control_set_ptr,
                context_ptr,
                &(refBuffer[puShiftYIndex * ref_stride + puShiftXIndex]),
                ref_stride,
                &context_ptr->p_best_ssd32x8[idx],
                puLcuBufferIndex,
                &(pos_b_buffer[posbBufferIndex]),
                &(pos_h_buffer[poshBufferIndex]),
                &(pos_j_buffer[posjBufferIndex]),
                32,
                8,
                x_search_area_origin,
                y_search_area_origin,
                &context_ptr->p_best_sad32x8[idx],
                &context_ptr->p_best_mv32x8[idx],
                &context_ptr->psub_pel_direction32x8[idx]);
        }

        for (pu_index = 0; pu_index < 16; ++pu_index) {
            idx = tab8x32[pu_index];

            puShiftXIndex = (pu_index & 0x07) << 3;
            puShiftYIndex = (pu_index >> 3) << 5;

            puLcuBufferIndex =
                puShiftXIndex + puShiftYIndex * context_ptr->sb_src_stride;
            posbBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;
            poshBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;
            posjBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;

            PU_HalfPelRefinement(
                sequence_control_set_ptr,
                context_ptr,
                &(refBuffer[puShiftYIndex * ref_stride + puShiftXIndex]),
                ref_stride,
                &context_ptr->p_best_ssd8x32[idx],
                puLcuBufferIndex,
                &(pos_b_buffer[posbBufferIndex]),
                &(pos_h_buffer[poshBufferIndex]),
                &(pos_j_buffer[posjBufferIndex]),
                8,
                32,
                x_search_area_origin,
                y_search_area_origin,
                &context_ptr->p_best_sad8x32[idx],
                &context_ptr->p_best_mv8x32[idx],
                &context_ptr->psub_pel_direction8x32[idx]);
        }

        for (pu_index = 0; pu_index < 4; ++pu_index) {
            idx = pu_index;

            puShiftXIndex = 0;
            puShiftYIndex = pu_index << 4;

            puLcuBufferIndex =
                puShiftXIndex + puShiftYIndex * context_ptr->sb_src_stride;
            posbBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;
            poshBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;
            posjBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;

            PU_HalfPelRefinement(
                sequence_control_set_ptr,
                context_ptr,
                &(refBuffer[puShiftYIndex * ref_stride + puShiftXIndex]),
                ref_stride,
                &context_ptr->p_best_ssd64x16[idx],
                puLcuBufferIndex,
                &(pos_b_buffer[posbBufferIndex]),
                &(pos_h_buffer[poshBufferIndex]),
                &(pos_j_buffer[posjBufferIndex]),
                64,
                16,
                x_search_area_origin,
                y_search_area_origin,
                &context_ptr->p_best_sad64x16[idx],
                &context_ptr->p_best_mv64x16[idx],
                &context_ptr->psub_pel_direction64x16[idx]);
        }

        for (pu_index = 0; pu_index < 4; ++pu_index) {
            idx = pu_index;

            puShiftXIndex = pu_index << 4;
            puShiftYIndex = 0;

            puLcuBufferIndex =
                puShiftXIndex + puShiftYIndex * context_ptr->sb_src_stride;
            posbBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;
            poshBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;
            posjBufferIndex = puShiftXIndex +
                              puShiftYIndex * context_ptr->interpolated_stride;

            PU_HalfPelRefinement(
                sequence_control_set_ptr,
                context_ptr,
                &(refBuffer[puShiftYIndex * ref_stride + puShiftXIndex]),
                ref_stride,
                &context_ptr->p_best_ssd16x64[idx],
                puLcuBufferIndex,
                &(pos_b_buffer[posbBufferIndex]),
                &(pos_h_buffer[poshBufferIndex]),
                &(pos_j_buffer[posjBufferIndex]),
                16,
                64,
                x_search_area_origin,
                y_search_area_origin,
                &context_ptr->p_best_sad16x64[idx],
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
uint32_t combined_averaging_ssd_c(uint8_t *src, ptrdiff_t src_stride,
                                  uint8_t *ref1, ptrdiff_t ref1_stride,
                                  uint8_t *ref2, ptrdiff_t ref2_stride,
                                  uint32_t height, uint32_t width) {
    uint32_t x, y;
    uint32_t ssd = 0;
    uint8_t avgpel;
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
 * PU_QuarterPelRefinementOnTheFly
 *   performs Quarter Pel refinement for each PU
 *******************************************/
static void PU_QuarterPelRefinementOnTheFly(
    MeContext *context_ptr,  // [IN] ME context Ptr, used to get SB Ptr
    uint32_t *pBestSsd,
    uint32_t
        puLcuBufferIndex,  // [IN] PU origin, used to point to source samples
    uint8_t **buf1,        // [IN]
    uint32_t *buf1Stride,
    uint8_t **buf2,  // [IN]
    uint32_t *buf2Stride,
    uint32_t pu_width,   // [IN]  PU width
    uint32_t pu_height,  // [IN]  PU height
    int16_t
        x_search_area_origin,  // [IN] search area origin in the horizontal
                               // direction, used to point to reference samples
    int16_t
        y_search_area_origin,  // [IN] search area origin in the vertical
                               // direction, used to point to reference samples
    uint32_t *pBestSad, uint32_t *pBestMV,
    uint8_t sub_pel_direction) {
    int16_t x_mv = _MVXT(*pBestMV);
    int16_t y_mv = _MVYT(*pBestMV);

    int16_t xSearchIndex = ((x_mv + 2) >> 2) - x_search_area_origin;
    int16_t ySearchIndex = ((y_mv + 2) >> 2) - y_search_area_origin;

    uint64_t dist;

    EbBool validTL, validT, validTR, validR, validBR, validB, validBL, validL;

    int16_t xMvQuarter[8];
    int16_t yMvQuarter[8];
    int32_t searchRegionIndex1 = 0;
    int32_t searchRegionIndex2 = 0;
    if (context_ptr->full_quarter_pel_refinement) {
        validTL = EB_TRUE;
        validT = EB_TRUE;
        validTR = EB_TRUE;
        validR = EB_TRUE;
        validBR = EB_TRUE;
        validB = EB_TRUE;
        validBL = EB_TRUE;
        validL = EB_TRUE;
    } else {
        if ((y_mv & 2) + ((x_mv & 2) >> 1)) {
            validTL = (EbBool)(sub_pel_direction == RIGHT_POSITION ||
                               sub_pel_direction == BOTTOM_RIGHT_POSITION ||
                               sub_pel_direction == BOTTOM_POSITION);
            validT = (EbBool)(sub_pel_direction == BOTTOM_RIGHT_POSITION ||
                              sub_pel_direction == BOTTOM_POSITION ||
                              sub_pel_direction == BOTTOM_LEFT_POSITION);
            validTR = (EbBool)(sub_pel_direction == BOTTOM_POSITION ||
                               sub_pel_direction == BOTTOM_LEFT_POSITION ||
                               sub_pel_direction == LEFT_POSITION);
            validR = (EbBool)(sub_pel_direction == BOTTOM_LEFT_POSITION ||
                              sub_pel_direction == LEFT_POSITION ||
                              sub_pel_direction == TOP_LEFT_POSITION);
            validBR = (EbBool)(sub_pel_direction == LEFT_POSITION ||
                               sub_pel_direction == TOP_LEFT_POSITION ||
                               sub_pel_direction == TOP_POSITION);
            validB = (EbBool)(sub_pel_direction == TOP_LEFT_POSITION ||
                              sub_pel_direction == TOP_POSITION ||
                              sub_pel_direction == TOP_RIGHT_POSITION);
            validBL = (EbBool)(sub_pel_direction == TOP_POSITION ||
                               sub_pel_direction == TOP_RIGHT_POSITION ||
                               sub_pel_direction == RIGHT_POSITION);
            validL = (EbBool)(sub_pel_direction == TOP_RIGHT_POSITION ||
                              sub_pel_direction == RIGHT_POSITION ||
                              sub_pel_direction == BOTTOM_RIGHT_POSITION);
        } else {
            validTL = (EbBool)(sub_pel_direction == LEFT_POSITION ||
                               sub_pel_direction == TOP_LEFT_POSITION ||
                               sub_pel_direction == TOP_POSITION);
            validT = (EbBool)(sub_pel_direction == TOP_LEFT_POSITION ||
                              sub_pel_direction == TOP_POSITION ||
                              sub_pel_direction == TOP_RIGHT_POSITION);
            validTR = (EbBool)(sub_pel_direction == TOP_POSITION ||
                               sub_pel_direction == TOP_RIGHT_POSITION ||
                               sub_pel_direction == RIGHT_POSITION);
            validR = (EbBool)(sub_pel_direction == TOP_RIGHT_POSITION ||
                              sub_pel_direction == RIGHT_POSITION ||
                              sub_pel_direction == BOTTOM_RIGHT_POSITION);
            validBR = (EbBool)(sub_pel_direction == RIGHT_POSITION ||
                               sub_pel_direction == BOTTOM_RIGHT_POSITION ||
                               sub_pel_direction == BOTTOM_POSITION);
            validB = (EbBool)(sub_pel_direction == BOTTOM_RIGHT_POSITION ||
                              sub_pel_direction == BOTTOM_POSITION ||
                              sub_pel_direction == BOTTOM_LEFT_POSITION);
            validBL = (EbBool)(sub_pel_direction == BOTTOM_POSITION ||
                               sub_pel_direction == BOTTOM_LEFT_POSITION ||
                               sub_pel_direction == LEFT_POSITION);
            validL = (EbBool)(sub_pel_direction == BOTTOM_LEFT_POSITION ||
                              sub_pel_direction == LEFT_POSITION ||
                              sub_pel_direction == TOP_LEFT_POSITION);
        }
    }
    xMvQuarter[0] = x_mv - 1;  // L  position
    xMvQuarter[1] = x_mv + 1;  // R  position
    xMvQuarter[2] = x_mv;      // T  position
    xMvQuarter[3] = x_mv;      // B  position
    xMvQuarter[4] = x_mv - 1;  // TL position
    xMvQuarter[5] = x_mv + 1;  // TR position
    xMvQuarter[6] = x_mv + 1;  // BR position
    xMvQuarter[7] = x_mv - 1;  // BL position

    yMvQuarter[0] = y_mv;      // L  position
    yMvQuarter[1] = y_mv;      // R  position
    yMvQuarter[2] = y_mv - 1;  // T  position
    yMvQuarter[3] = y_mv + 1;  // B  position
    yMvQuarter[4] = y_mv - 1;  // TL position
    yMvQuarter[5] = y_mv - 1;  // TR position
    yMvQuarter[6] = y_mv + 1;  // BR position
    yMvQuarter[7] = y_mv + 1;  // BL position

    // Use SATD only when QP mod, and RC are OFF
    // QP mod, and RC assume that ME distotion is always SAD.
    // This problem might be solved by computing SAD for the best position after
    // fractional search is done, or by considring the full pel resolution SAD.

    {
        // L position
        if (validL) {
            searchRegionIndex1 = (int32_t)xSearchIndex +
                                 (int32_t)buf1Stride[0] * (int32_t)ySearchIndex;
            searchRegionIndex2 = (int32_t)xSearchIndex +
                                 (int32_t)buf2Stride[0] * (int32_t)ySearchIndex;

            dist =
                (context_ptr->fractional_search_method == SSD_SEARCH)
                    ? combined_averaging_ssd(
                          &(context_ptr->sb_buffer[puLcuBufferIndex]),
                          BLOCK_SIZE_64,
                          buf1[0] + searchRegionIndex1,
                          buf1Stride[0],
                          buf2[0] + searchRegionIndex2,
                          buf2Stride[0],
                          pu_height,
                          pu_width)
                    : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                          ? (nxm_sad_avg_kernel(
                                     &(context_ptr
                                           ->sb_buffer[puLcuBufferIndex]),
                                     BLOCK_SIZE_64 << 1,
                                     buf1[0] + searchRegionIndex1,
                                     buf1Stride[0] << 1,
                                     buf2[0] + searchRegionIndex2,
                                     buf2Stride[0] << 1,
                                     pu_height >> 1,
                                     pu_width))
                                << 1
                          : nxm_sad_avg_kernel(
                                    &(context_ptr->sb_buffer[puLcuBufferIndex]),
                                    BLOCK_SIZE_64,
                                    buf1[0] + searchRegionIndex1,
                                    buf1Stride[0],
                                    buf2[0] + searchRegionIndex2,
                                    buf2Stride[0],
                                    pu_height,
                                    pu_width);
            if (context_ptr->fractional_search_method == SSD_SEARCH) {
                if (dist < *pBestSsd) {
                    *pBestSad =
                        (uint32_t)nxm_sad_avg_kernel(
                                &(context_ptr->sb_buffer[puLcuBufferIndex]),
                                BLOCK_SIZE_64,
                                buf1[0] + searchRegionIndex1,
                                buf1Stride[0],
                                buf2[0] + searchRegionIndex2,
                                buf2Stride[0],
                                pu_height,
                                pu_width);
                    *pBestMV = ((uint16_t)yMvQuarter[0] << 16) |
                               ((uint16_t)xMvQuarter[0]);
                    *pBestSsd = (uint32_t)dist;
                }
            } else {
                if (dist < *pBestSad) {
                    *pBestSad = (uint32_t)dist;
                    *pBestMV = ((uint16_t)yMvQuarter[0] << 16) |
                               ((uint16_t)xMvQuarter[0]);
                }
            }
        }

        // R positions
        if (validR) {
            searchRegionIndex1 = (int32_t)xSearchIndex +
                                 (int32_t)buf1Stride[1] * (int32_t)ySearchIndex;
            searchRegionIndex2 = (int32_t)xSearchIndex +
                                 (int32_t)buf2Stride[1] * (int32_t)ySearchIndex;
            dist =
                (context_ptr->fractional_search_method == SSD_SEARCH)
                    ? combined_averaging_ssd(
                          &(context_ptr->sb_buffer[puLcuBufferIndex]),
                          BLOCK_SIZE_64,
                          buf1[1] + searchRegionIndex1,
                          buf1Stride[1],
                          buf2[1] + searchRegionIndex2,
                          buf2Stride[1],
                          pu_height,
                          pu_width)
                    : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                          ? (nxm_sad_avg_kernel(
                                     &(context_ptr
                                           ->sb_buffer[puLcuBufferIndex]),
                                     BLOCK_SIZE_64 << 1,
                                     buf1[1] + searchRegionIndex1,
                                     buf1Stride[1] << 1,
                                     buf2[1] + searchRegionIndex2,
                                     buf2Stride[1] << 1,
                                     pu_height >> 1,
                                     pu_width))
                                << 1
                          : nxm_sad_avg_kernel(
                                    &(context_ptr->sb_buffer[puLcuBufferIndex]),
                                    BLOCK_SIZE_64,
                                    buf1[1] + searchRegionIndex1,
                                    buf1Stride[1],
                                    buf2[1] + searchRegionIndex2,
                                    buf2Stride[1],
                                    pu_height,
                                    pu_width);
            if (context_ptr->fractional_search_method == SSD_SEARCH) {
                if (dist < *pBestSsd) {
                    *pBestSad =
                        (uint32_t)nxm_sad_avg_kernel(
                                &(context_ptr->sb_buffer[puLcuBufferIndex]),
                                BLOCK_SIZE_64,
                                buf1[1] + searchRegionIndex1,
                                buf1Stride[1],
                                buf2[1] + searchRegionIndex2,
                                buf2Stride[1],
                                pu_height,
                                pu_width);
                    *pBestMV = ((uint16_t)yMvQuarter[1] << 16) |
                               ((uint16_t)xMvQuarter[1]);
                    *pBestSsd = (uint32_t)dist;
                }
            } else {
                if (dist < *pBestSad) {
                    *pBestSad = (uint32_t)dist;
                    *pBestMV = ((uint16_t)yMvQuarter[1] << 16) |
                               ((uint16_t)xMvQuarter[1]);
                }
            }
        }

        // T position
        if (validT) {
            searchRegionIndex1 = (int32_t)xSearchIndex +
                                 (int32_t)buf1Stride[2] * (int32_t)ySearchIndex;
            searchRegionIndex2 = (int32_t)xSearchIndex +
                                 (int32_t)buf2Stride[2] * (int32_t)ySearchIndex;
            dist =
                (context_ptr->fractional_search_method == SSD_SEARCH)
                    ? combined_averaging_ssd(
                          &(context_ptr->sb_buffer[puLcuBufferIndex]),
                          BLOCK_SIZE_64,
                          buf1[2] + searchRegionIndex1,
                          buf1Stride[2],
                          buf2[2] + searchRegionIndex2,
                          buf2Stride[2],
                          pu_height,
                          pu_width)
                    : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                          ? (nxm_sad_avg_kernel(
                                     &(context_ptr
                                           ->sb_buffer[puLcuBufferIndex]),
                                     BLOCK_SIZE_64 << 1,
                                     buf1[2] + searchRegionIndex1,
                                     buf1Stride[2] << 1,
                                     buf2[2] + searchRegionIndex2,
                                     buf2Stride[2] << 1,
                                     pu_height >> 1,
                                     pu_width))
                                << 1
                          : nxm_sad_avg_kernel(
                                    &(context_ptr->sb_buffer[puLcuBufferIndex]),
                                    BLOCK_SIZE_64,
                                    buf1[2] + searchRegionIndex1,
                                    buf1Stride[2],
                                    buf2[2] + searchRegionIndex2,
                                    buf2Stride[2],
                                    pu_height,
                                    pu_width);
            if (context_ptr->fractional_search_method == SSD_SEARCH) {
                if (dist < *pBestSsd) {
                    *pBestSad =
                        (uint32_t)nxm_sad_avg_kernel(
                                &(context_ptr->sb_buffer[puLcuBufferIndex]),
                                BLOCK_SIZE_64,
                                buf1[2] + searchRegionIndex1,
                                buf1Stride[2],
                                buf2[2] + searchRegionIndex2,
                                buf2Stride[2],
                                pu_height,
                                pu_width);
                    *pBestMV = ((uint16_t)yMvQuarter[2] << 16) |
                               ((uint16_t)xMvQuarter[2]);
                    *pBestSsd = (uint32_t)dist;
                }
            } else {
                if (dist < *pBestSad) {
                    *pBestSad = (uint32_t)dist;
                    *pBestMV = ((uint16_t)yMvQuarter[2] << 16) |
                               ((uint16_t)xMvQuarter[2]);
                }
            }
        }

        // B position
        if (validB) {
            searchRegionIndex1 = (int32_t)xSearchIndex +
                                 (int32_t)buf1Stride[3] * (int32_t)ySearchIndex;
            searchRegionIndex2 = (int32_t)xSearchIndex +
                                 (int32_t)buf2Stride[3] * (int32_t)ySearchIndex;
            dist =
                (context_ptr->fractional_search_method == SSD_SEARCH)
                    ? combined_averaging_ssd(
                          &(context_ptr->sb_buffer[puLcuBufferIndex]),
                          BLOCK_SIZE_64,
                          buf1[3] + searchRegionIndex1,
                          buf1Stride[3],
                          buf2[3] + searchRegionIndex2,
                          buf2Stride[3],
                          pu_height,
                          pu_width)
                    : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                          ? (nxm_sad_avg_kernel(
                                     &(context_ptr
                                           ->sb_buffer[puLcuBufferIndex]),
                                     BLOCK_SIZE_64 << 1,
                                     buf1[3] + searchRegionIndex1,
                                     buf1Stride[3] << 1,
                                     buf2[3] + searchRegionIndex2,
                                     buf2Stride[3] << 1,
                                     pu_height >> 1,
                                     pu_width))
                                << 1
                          : nxm_sad_avg_kernel(
                                    &(context_ptr->sb_buffer[puLcuBufferIndex]),
                                    BLOCK_SIZE_64,
                                    buf1[3] + searchRegionIndex1,
                                    buf1Stride[3],
                                    buf2[3] + searchRegionIndex2,
                                    buf2Stride[3],
                                    pu_height,
                                    pu_width);
            if (context_ptr->fractional_search_method == SSD_SEARCH) {
                if (dist < *pBestSsd) {
                    *pBestSad =
                        (uint32_t)nxm_sad_avg_kernel(
                                &(context_ptr->sb_buffer[puLcuBufferIndex]),
                                BLOCK_SIZE_64,
                                buf1[3] + searchRegionIndex1,
                                buf1Stride[3],
                                buf2[3] + searchRegionIndex2,
                                buf2Stride[3],
                                pu_height,
                                pu_width);
                    *pBestMV = ((uint16_t)yMvQuarter[3] << 16) |
                               ((uint16_t)xMvQuarter[3]);
                    *pBestSsd = (uint32_t)dist;
                }
            } else {
                if (dist < *pBestSad) {
                    *pBestSad = (uint32_t)dist;
                    *pBestMV = ((uint16_t)yMvQuarter[3] << 16) |
                               ((uint16_t)xMvQuarter[3]);
                }
            }
        }

        // TL position
        if (validTL) {
            searchRegionIndex1 = (int32_t)xSearchIndex +
                                 (int32_t)buf1Stride[4] * (int32_t)ySearchIndex;
            searchRegionIndex2 = (int32_t)xSearchIndex +
                                 (int32_t)buf2Stride[4] * (int32_t)ySearchIndex;
            dist =
                (context_ptr->fractional_search_method == SSD_SEARCH)
                    ? combined_averaging_ssd(
                          &(context_ptr->sb_buffer[puLcuBufferIndex]),
                          BLOCK_SIZE_64,
                          buf1[4] + searchRegionIndex1,
                          buf1Stride[4],
                          buf2[4] + searchRegionIndex2,
                          buf2Stride[4],
                          pu_height,
                          pu_width)
                    : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                          ? (nxm_sad_avg_kernel(
                                     &(context_ptr
                                           ->sb_buffer[puLcuBufferIndex]),
                                     BLOCK_SIZE_64 << 1,
                                     buf1[4] + searchRegionIndex1,
                                     buf1Stride[4] << 1,
                                     buf2[4] + searchRegionIndex2,
                                     buf2Stride[4] << 1,
                                     pu_height >> 1,
                                     pu_width))
                                << 1
                          : nxm_sad_avg_kernel(
                                    &(context_ptr->sb_buffer[puLcuBufferIndex]),
                                    BLOCK_SIZE_64,
                                    buf1[4] + searchRegionIndex1,
                                    buf1Stride[4],
                                    buf2[4] + searchRegionIndex2,
                                    buf2Stride[4],
                                    pu_height,
                                    pu_width);
            if (context_ptr->fractional_search_method == SSD_SEARCH) {
                if (dist < *pBestSsd) {
                    *pBestSad =
                        (uint32_t)nxm_sad_avg_kernel(
                                &(context_ptr->sb_buffer[puLcuBufferIndex]),
                                BLOCK_SIZE_64,
                                buf1[4] + searchRegionIndex1,
                                buf1Stride[4],
                                buf2[4] + searchRegionIndex2,
                                buf2Stride[4],
                                pu_height,
                                pu_width);
                    *pBestMV = ((uint16_t)yMvQuarter[4] << 16) |
                               ((uint16_t)xMvQuarter[4]);
                    *pBestSsd = (uint32_t)dist;
                }
            } else {
                if (dist < *pBestSad) {
                    *pBestSad = (uint32_t)dist;
                    *pBestMV = ((uint16_t)yMvQuarter[4] << 16) |
                               ((uint16_t)xMvQuarter[4]);
                }
            }
        }

        // TR position
        if (validTR) {
            searchRegionIndex1 = (int32_t)xSearchIndex +
                                 (int32_t)buf1Stride[5] * (int32_t)ySearchIndex;
            searchRegionIndex2 = (int32_t)xSearchIndex +
                                 (int32_t)buf2Stride[5] * (int32_t)ySearchIndex;
            dist =
                (context_ptr->fractional_search_method == SSD_SEARCH)
                    ? combined_averaging_ssd(
                          &(context_ptr->sb_buffer[puLcuBufferIndex]),
                          BLOCK_SIZE_64,
                          buf1[5] + searchRegionIndex1,
                          buf1Stride[5],
                          buf2[5] + searchRegionIndex2,
                          buf2Stride[5],
                          pu_height,
                          pu_width)
                    : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                          ? (nxm_sad_avg_kernel(
                                     &(context_ptr
                                           ->sb_buffer[puLcuBufferIndex]),
                                     BLOCK_SIZE_64 << 1,
                                     buf1[5] + searchRegionIndex1,
                                     buf1Stride[5] << 1,
                                     buf2[5] + searchRegionIndex2,
                                     buf2Stride[5] << 1,
                                     pu_height >> 1,
                                     pu_width))
                                << 1
                          : nxm_sad_avg_kernel(
                                    &(context_ptr->sb_buffer[puLcuBufferIndex]),
                                    BLOCK_SIZE_64,
                                    buf1[5] + searchRegionIndex1,
                                    buf1Stride[5],
                                    buf2[5] + searchRegionIndex2,
                                    buf2Stride[5],
                                    pu_height,
                                    pu_width);
            if (context_ptr->fractional_search_method == SSD_SEARCH) {
                if (dist < *pBestSsd) {
                    *pBestSad =
                        (uint32_t)nxm_sad_avg_kernel(
                                &(context_ptr->sb_buffer[puLcuBufferIndex]),
                                BLOCK_SIZE_64,
                                buf1[5] + searchRegionIndex1,
                                buf1Stride[5],
                                buf2[5] + searchRegionIndex2,
                                buf2Stride[5],
                                pu_height,
                                pu_width);
                    *pBestMV = ((uint16_t)yMvQuarter[5] << 16) |
                               ((uint16_t)xMvQuarter[5]);
                    *pBestSsd = (uint32_t)dist;
                }
            } else {
                if (dist < *pBestSad) {
                    *pBestSad = (uint32_t)dist;
                    *pBestMV = ((uint16_t)yMvQuarter[5] << 16) |
                               ((uint16_t)xMvQuarter[5]);
                }
            }
        }

        // BR position
        if (validBR) {
            searchRegionIndex1 = (int32_t)xSearchIndex +
                                 (int32_t)buf1Stride[6] * (int32_t)ySearchIndex;
            searchRegionIndex2 = (int32_t)xSearchIndex +
                                 (int32_t)buf2Stride[6] * (int32_t)ySearchIndex;
            dist =
                (context_ptr->fractional_search_method == SSD_SEARCH)
                    ? combined_averaging_ssd(
                          &(context_ptr->sb_buffer[puLcuBufferIndex]),
                          BLOCK_SIZE_64,
                          buf1[6] + searchRegionIndex1,
                          buf1Stride[6],
                          buf2[6] + searchRegionIndex2,
                          buf2Stride[6],
                          pu_height,
                          pu_width)
                    : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                          ? (nxm_sad_avg_kernel(
                                     &(context_ptr
                                           ->sb_buffer[puLcuBufferIndex]),
                                     BLOCK_SIZE_64 << 1,
                                     buf1[6] + searchRegionIndex1,
                                     buf1Stride[6] << 1,
                                     buf2[6] + searchRegionIndex2,
                                     buf2Stride[6] << 1,
                                     pu_height >> 1,
                                     pu_width))
                                << 1
                          : nxm_sad_avg_kernel(
                                    &(context_ptr->sb_buffer[puLcuBufferIndex]),
                                    BLOCK_SIZE_64,
                                    buf1[6] + searchRegionIndex1,
                                    buf1Stride[6],
                                    buf2[6] + searchRegionIndex2,
                                    buf2Stride[6],
                                    pu_height,
                                    pu_width);
            if (context_ptr->fractional_search_method == SSD_SEARCH) {
                if (dist < *pBestSsd) {
                    *pBestSad =
                        (uint32_t)nxm_sad_avg_kernel(
                                &(context_ptr->sb_buffer[puLcuBufferIndex]),
                                BLOCK_SIZE_64,
                                buf1[6] + searchRegionIndex1,
                                buf1Stride[6],
                                buf2[6] + searchRegionIndex2,
                                buf2Stride[6],
                                pu_height,
                                pu_width);
                    *pBestMV = ((uint16_t)yMvQuarter[6] << 16) |
                               ((uint16_t)xMvQuarter[6]);
                    *pBestSsd = (uint32_t)dist;
                }
            } else {
                if (dist < *pBestSad) {
                    *pBestSad = (uint32_t)dist;
                    *pBestMV = ((uint16_t)yMvQuarter[6] << 16) |
                               ((uint16_t)xMvQuarter[6]);
                }
            }
        }

        // BL position
        if (validBL) {
            searchRegionIndex1 = (int32_t)xSearchIndex +
                                 (int32_t)buf1Stride[7] * (int32_t)ySearchIndex;
            searchRegionIndex2 = (int32_t)xSearchIndex +
                                 (int32_t)buf2Stride[7] * (int32_t)ySearchIndex;
            dist =
                (context_ptr->fractional_search_method == SSD_SEARCH)
                    ? combined_averaging_ssd(
                          &(context_ptr->sb_buffer[puLcuBufferIndex]),
                          BLOCK_SIZE_64,
                          buf1[7] + searchRegionIndex1,
                          buf1Stride[7],
                          buf2[7] + searchRegionIndex2,
                          buf2Stride[7],
                          pu_height,
                          pu_width)
                    : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                          ? (nxm_sad_avg_kernel(
                                     &(context_ptr
                                           ->sb_buffer[puLcuBufferIndex]),
                                     BLOCK_SIZE_64 << 1,
                                     buf1[7] + searchRegionIndex1,
                                     buf1Stride[7] << 1,
                                     buf2[7] + searchRegionIndex2,
                                     buf2Stride[7] << 1,
                                     pu_height >> 1,
                                     pu_width))
                                << 1
                          : nxm_sad_avg_kernel(
                                    &(context_ptr->sb_buffer[puLcuBufferIndex]),
                                    BLOCK_SIZE_64,
                                    buf1[7] + searchRegionIndex1,
                                    buf1Stride[7],
                                    buf2[7] + searchRegionIndex2,
                                    buf2Stride[7],
                                    pu_height,
                                    pu_width);
            if (context_ptr->fractional_search_method == SSD_SEARCH) {
                if (dist < *pBestSsd) {
                    *pBestSad =
                        (uint32_t)nxm_sad_avg_kernel(
                                &(context_ptr->sb_buffer[puLcuBufferIndex]),
                                BLOCK_SIZE_64,
                                buf1[7] + searchRegionIndex1,
                                buf1Stride[7],
                                buf2[7] + searchRegionIndex2,
                                buf2Stride[7],
                                pu_height,
                                pu_width);
                    *pBestMV = ((uint16_t)yMvQuarter[7] << 16) |
                               ((uint16_t)xMvQuarter[7]);
                    *pBestSsd = (uint32_t)dist;
                }
            } else {
                if (dist < *pBestSad) {
                    *pBestSad = (uint32_t)dist;
                    *pBestMV = ((uint16_t)yMvQuarter[7] << 16) |
                               ((uint16_t)xMvQuarter[7]);
                }
            }
        }
    }

    return;
}

/*******************************************
* SetQuarterPelRefinementInputsOnTheFly
*   determine the 2 half pel buffers to do
averaging for Quarter Pel Refinement
*******************************************/
static void SetQuarterPelRefinementInputsOnTheFly(
    uint8_t *pos_Full,     //[IN] points to A
    uint32_t FullStride,   //[IN]
    uint8_t *pos_b,        //[IN] points to b
    uint8_t *pos_h,        //[IN] points to h
    uint8_t *pos_j,        //[IN] points to j
    uint32_t Stride,       //[IN]
    int16_t x_mv,          //[IN]
    int16_t y_mv,          //[IN]
    uint8_t **buf1,        //[OUT]
    uint32_t *buf1Stride,  //[OUT]
    uint8_t **buf2,        //[OUT]
    uint32_t *buf2Stride   //[OUT]
) {
    uint32_t quarterPelRefinementMethod = (y_mv & 2) + ((x_mv & 2) >> 1);

    // for each one of the 8 postions, we need to determine the 2 half pel
    // buffers to  do averaging

    //     A    a    b    c
    //     d    e    f    g
    //     h    i    j    k
    //     n    p    q    r

    switch (quarterPelRefinementMethod) {
    case EB_QUARTER_IN_FULL:

        /*c=b+A*/ buf1[0] = pos_b;
        buf1Stride[0] = Stride;
        buf2[0] = pos_Full;
        buf2Stride[0] = FullStride;
        /*a=A+b*/ buf1[1] = pos_Full;
        buf1Stride[1] = FullStride;
        buf2[1] = pos_b + 1;
        buf2Stride[1] = Stride;
        /*n=h+A*/ buf1[2] = pos_h;
        buf1Stride[2] = Stride;
        buf2[2] = pos_Full;
        buf2Stride[2] = FullStride;
        /*d=A+h*/ buf1[3] = pos_Full;
        buf1Stride[3] = FullStride;
        buf2[3] = pos_h + Stride;
        buf2Stride[3] = Stride;
        /*r=b+h*/ buf1[4] = pos_b;
        buf1Stride[4] = Stride;
        buf2[4] = pos_h;
        buf2Stride[4] = Stride;
        /*p=h+b*/ buf1[5] = pos_h;
        buf1Stride[5] = Stride;
        buf2[5] = pos_b + 1;
        buf2Stride[5] = Stride;
        /*e=h+b*/ buf1[6] = pos_h + Stride;
        buf1Stride[6] = Stride;
        buf2[6] = pos_b + 1;
        buf2Stride[6] = Stride;
        /*g=b+h*/ buf1[7] = pos_b;
        buf1Stride[7] = Stride;
        buf2[7] = pos_h + Stride;
        buf2Stride[7] = Stride;

        break;

    case EB_QUARTER_IN_HALF_HORIZONTAL:

        /*a=A+b*/ buf1[0] = pos_Full - 1;
        buf1Stride[0] = FullStride;
        buf2[0] = pos_b;
        buf2Stride[0] = Stride;
        /*c=b+A*/ buf1[1] = pos_b;
        buf1Stride[1] = Stride;
        buf2[1] = pos_Full;
        buf2Stride[1] = FullStride;
        /*q=j+b*/ buf1[2] = pos_j;
        buf1Stride[2] = Stride;
        buf2[2] = pos_b;
        buf2Stride[2] = Stride;
        /*f=b+j*/ buf1[3] = pos_b;
        buf1Stride[3] = Stride;
        buf2[3] = pos_j + Stride;
        buf2Stride[3] = Stride;
        /*p=h+b*/ buf1[4] = pos_h - 1;
        buf1Stride[4] = Stride;
        buf2[4] = pos_b;
        buf2Stride[4] = Stride;
        /*r=b+h*/ buf1[5] = pos_b;
        buf1Stride[5] = Stride;
        buf2[5] = pos_h;
        buf2Stride[5] = Stride;
        /*g=b+h*/ buf1[6] = pos_b;
        buf1Stride[6] = Stride;
        buf2[6] = pos_h + Stride;
        buf2Stride[6] = Stride;
        /*e=h+b*/ buf1[7] = pos_h - 1 + Stride;
        buf1Stride[7] = Stride;
        buf2[7] = pos_b;
        buf2Stride[7] = Stride;

        break;

    case EB_QUARTER_IN_HALF_VERTICAL:

        /*k=j+h*/ buf1[0] = pos_j;
        buf1Stride[0] = Stride;
        buf2[0] = pos_h;
        buf2Stride[0] = Stride;
        /*i=h+j*/ buf1[1] = pos_h;
        buf1Stride[1] = Stride;
        buf2[1] = pos_j + 1;
        buf2Stride[1] = Stride;
        /*d=A+h*/ buf1[2] = pos_Full - FullStride;
        buf1Stride[2] = FullStride;
        buf2[2] = pos_h;
        buf2Stride[2] = Stride;
        /*n=h+A*/ buf1[3] = pos_h;
        buf1Stride[3] = Stride;
        buf2[3] = pos_Full;
        buf2Stride[3] = FullStride;
        /*g=b+h*/ buf1[4] = pos_b - Stride;
        buf1Stride[4] = Stride;
        buf2[4] = pos_h;
        buf2Stride[4] = Stride;
        /*e=h+b*/ buf1[5] = pos_h;
        buf1Stride[5] = Stride;
        buf2[5] = pos_b + 1 - Stride;
        buf2Stride[5] = Stride;
        /*p=h+b*/ buf1[6] = pos_h;
        buf1Stride[6] = Stride;
        buf2[6] = pos_b + 1;
        buf2Stride[6] = Stride;
        /*r=b+h*/ buf1[7] = pos_b;
        buf1Stride[7] = Stride;
        buf2[7] = pos_h;
        buf2Stride[7] = Stride;

        break;

    case EB_QUARTER_IN_HALF_DIAGONAL:

        /*i=h+j*/ buf1[0] = pos_h - 1;
        buf1Stride[0] = Stride;
        buf2[0] = pos_j;
        buf2Stride[0] = Stride;
        /*k=j+h*/ buf1[1] = pos_j;
        buf1Stride[1] = Stride;
        buf2[1] = pos_h;
        buf2Stride[1] = Stride;
        /*f=b+j*/ buf1[2] = pos_b - Stride;
        buf1Stride[2] = Stride;
        buf2[2] = pos_j;
        buf2Stride[2] = Stride;
        /*q=j+b*/ buf1[3] = pos_j;
        buf1Stride[3] = Stride;
        buf2[3] = pos_b;
        buf2Stride[3] = Stride;
        /*e=h+b*/ buf1[4] = pos_h - 1;
        buf1Stride[4] = Stride;
        buf2[4] = pos_b - Stride;
        buf2Stride[4] = Stride;
        /*g=b+h*/ buf1[5] = pos_b - Stride;
        buf1Stride[5] = Stride;
        buf2[5] = pos_h;
        buf2Stride[5] = Stride;
        /*r=b+h*/ buf1[6] = pos_b;
        buf1Stride[6] = Stride;
        buf2[6] = pos_h;
        buf2Stride[6] = Stride;
        /*p=h+b*/ buf1[7] = pos_h - 1;
        buf1Stride[7] = Stride;
        buf2[7] = pos_b;
        buf2Stride[7] = Stride;

        break;

    default: break;
    }

    return;
}

/*******************************************
 * QuarterPelSearch_LCU
 *   performs Quarter Pel refinement for the 85 PUs
 *******************************************/
static void QuarterPelSearch_LCU(
    MeContext
        *context_ptr,  //[IN/OUT]  ME context Ptr, used to get/update ME results
    uint8_t *pos_Full,    //[IN]
    uint32_t FullStride,  //[IN]
    uint8_t *pos_b,       //[IN]
    uint8_t *pos_h,       //[IN]
    uint8_t *pos_j,       //[IN]
    int16_t
        x_search_area_origin,  //[IN] search area origin in the horizontal
                               // direction, used to point to reference samples
    int16_t
        y_search_area_origin,  //[IN] search area origin in the vertical
                               // direction, used to point to reference samples
    EbBool disable8x8CuInMeFlag,
    EbBool enable_half_pel32x32, EbBool enable_half_pel16x16,
    EbBool enable_half_pel8x8,
    EbBool enableQuarterPel, EbBool ext_block_flag)
{
    uint32_t pu_index;

    uint32_t puShiftXIndex;
    uint32_t puShiftYIndex;

    uint32_t puLcuBufferIndex;

    // for each one of the 8 positions, we need to determine the 2 buffers to do
    // averaging
    uint8_t *buf1[8];
    uint8_t *buf2[8];

    uint32_t buf1Stride[8];
    uint32_t buf2Stride[8];

    int16_t x_mv, y_mv;
    uint32_t nidx;

    if (context_ptr->fractional_search64x64) {
        x_mv = _MVXT(*context_ptr->p_best_mv64x64);
        y_mv = _MVYT(*context_ptr->p_best_mv64x64);

        SetQuarterPelRefinementInputsOnTheFly(pos_Full,
                                              FullStride,
                                              pos_b,
                                              pos_h,
                                              pos_j,
                                              context_ptr->interpolated_stride,
                                              x_mv,
                                              y_mv,
                                              buf1,
                                              buf1Stride,
                                              buf2,
                                              buf2Stride);

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

        PU_QuarterPelRefinementOnTheFly(context_ptr,
                                        context_ptr->p_best_ssd64x64,
                                        0,
                                        buf1,
                                        buf1Stride,
                                        buf2,
                                        buf2Stride,
                                        64,
                                        64,
                                        x_search_area_origin,
                                        y_search_area_origin,
                                        context_ptr->p_best_sad64x64,
                                        context_ptr->p_best_mv64x64,
                                        context_ptr->psub_pel_direction64x64);
    }
    if (enableQuarterPel && enable_half_pel32x32)
    {
        // 32x32 [4 partitions]
        for (pu_index = 0; pu_index < 4; ++pu_index) {
            x_mv = _MVXT(context_ptr->p_best_mv32x32[pu_index]);
            y_mv = _MVYT(context_ptr->p_best_mv32x32[pu_index]);

            SetQuarterPelRefinementInputsOnTheFly(
                pos_Full,
                FullStride,
                pos_b,
                pos_h,
                pos_j,
                context_ptr->interpolated_stride,
                x_mv,
                y_mv,
                buf1,
                buf1Stride,
                buf2,
                buf2Stride);

            puShiftXIndex = (pu_index & 0x01) << 5;
            puShiftYIndex = (pu_index >> 1) << 5;

            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + puShiftXIndex + puShiftYIndex * buf1Stride[0];
            buf2[0] = buf2[0] + puShiftXIndex + puShiftYIndex * buf2Stride[0];
            buf1[1] = buf1[1] + puShiftXIndex + puShiftYIndex * buf1Stride[1];
            buf2[1] = buf2[1] + puShiftXIndex + puShiftYIndex * buf2Stride[1];
            buf1[2] = buf1[2] + puShiftXIndex + puShiftYIndex * buf1Stride[2];
            buf2[2] = buf2[2] + puShiftXIndex + puShiftYIndex * buf2Stride[2];
            buf1[3] = buf1[3] + puShiftXIndex + puShiftYIndex * buf1Stride[3];
            buf2[3] = buf2[3] + puShiftXIndex + puShiftYIndex * buf2Stride[3];
            buf1[4] = buf1[4] + puShiftXIndex + puShiftYIndex * buf1Stride[4];
            buf2[4] = buf2[4] + puShiftXIndex + puShiftYIndex * buf2Stride[4];
            buf1[5] = buf1[5] + puShiftXIndex + puShiftYIndex * buf1Stride[5];
            buf2[5] = buf2[5] + puShiftXIndex + puShiftYIndex * buf2Stride[5];
            buf1[6] = buf1[6] + puShiftXIndex + puShiftYIndex * buf1Stride[6];
            buf2[6] = buf2[6] + puShiftXIndex + puShiftYIndex * buf2Stride[6];
            buf1[7] = buf1[7] + puShiftXIndex + puShiftYIndex * buf1Stride[7];
            buf2[7] = buf2[7] + puShiftXIndex + puShiftYIndex * buf2Stride[7];

            PU_QuarterPelRefinementOnTheFly(
                context_ptr,
                &context_ptr->p_best_ssd32x32[pu_index],
                puLcuBufferIndex,
                buf1,
                buf1Stride,
                buf2,
                buf2Stride,
                32,
                32,
                x_search_area_origin,
                y_search_area_origin,
                &context_ptr->p_best_sad32x32[pu_index],
                &context_ptr->p_best_mv32x32[pu_index],
                context_ptr->psub_pel_direction32x32[pu_index]);
        }
    }

    if (enableQuarterPel && enable_half_pel16x16)
    {
        // 16x16 [16 partitions]
        for (pu_index = 0; pu_index < 16; ++pu_index) {
            nidx = tab16x16[pu_index];

            x_mv = _MVXT(context_ptr->p_best_mv16x16[nidx]);
            y_mv = _MVYT(context_ptr->p_best_mv16x16[nidx]);

            SetQuarterPelRefinementInputsOnTheFly(
                pos_Full,
                FullStride,
                pos_b,
                pos_h,
                pos_j,
                context_ptr->interpolated_stride,
                x_mv,
                y_mv,
                buf1,
                buf1Stride,
                buf2,
                buf2Stride);

            puShiftXIndex = (pu_index & 0x03) << 4;
            puShiftYIndex = (pu_index >> 2) << 4;

            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + puShiftXIndex + puShiftYIndex * buf1Stride[0];
            buf2[0] = buf2[0] + puShiftXIndex + puShiftYIndex * buf2Stride[0];
            buf1[1] = buf1[1] + puShiftXIndex + puShiftYIndex * buf1Stride[1];
            buf2[1] = buf2[1] + puShiftXIndex + puShiftYIndex * buf2Stride[1];
            buf1[2] = buf1[2] + puShiftXIndex + puShiftYIndex * buf1Stride[2];
            buf2[2] = buf2[2] + puShiftXIndex + puShiftYIndex * buf2Stride[2];
            buf1[3] = buf1[3] + puShiftXIndex + puShiftYIndex * buf1Stride[3];
            buf2[3] = buf2[3] + puShiftXIndex + puShiftYIndex * buf2Stride[3];
            buf1[4] = buf1[4] + puShiftXIndex + puShiftYIndex * buf1Stride[4];
            buf2[4] = buf2[4] + puShiftXIndex + puShiftYIndex * buf2Stride[4];
            buf1[5] = buf1[5] + puShiftXIndex + puShiftYIndex * buf1Stride[5];
            buf2[5] = buf2[5] + puShiftXIndex + puShiftYIndex * buf2Stride[5];
            buf1[6] = buf1[6] + puShiftXIndex + puShiftYIndex * buf1Stride[6];
            buf2[6] = buf2[6] + puShiftXIndex + puShiftYIndex * buf2Stride[6];
            buf1[7] = buf1[7] + puShiftXIndex + puShiftYIndex * buf1Stride[7];
            buf2[7] = buf2[7] + puShiftXIndex + puShiftYIndex * buf2Stride[7];

            PU_QuarterPelRefinementOnTheFly(
                context_ptr,
                &context_ptr->p_best_ssd16x16[nidx],
                puLcuBufferIndex,
                buf1,
                buf1Stride,
                buf2,
                buf2Stride,
                16,
                16,
                x_search_area_origin,
                y_search_area_origin,
                &context_ptr->p_best_sad16x16[nidx],
                &context_ptr->p_best_mv16x16[nidx],
                context_ptr->psub_pel_direction16x16[nidx]);
        }
    }

    if (enableQuarterPel && enable_half_pel8x8)
    {
        // 8x8   [64 partitions]
        if (!disable8x8CuInMeFlag) {
            for (pu_index = 0; pu_index < 64; ++pu_index) {
                nidx = tab8x8[pu_index];

                x_mv = _MVXT(context_ptr->p_best_mv8x8[nidx]);
                y_mv = _MVYT(context_ptr->p_best_mv8x8[nidx]);

                SetQuarterPelRefinementInputsOnTheFly(
                    pos_Full,
                    FullStride,
                    pos_b,
                    pos_h,
                    pos_j,
                    context_ptr->interpolated_stride,
                    x_mv,
                    y_mv,
                    buf1,
                    buf1Stride,
                    buf2,
                    buf2Stride);

                puShiftXIndex = (pu_index & 0x07) << 3;
                puShiftYIndex = (pu_index >> 3) << 3;

                puLcuBufferIndex =
                    puShiftXIndex + puShiftYIndex * BLOCK_SIZE_64;

                buf1[0] =
                    buf1[0] + puShiftXIndex + puShiftYIndex * buf1Stride[0];
                buf2[0] =
                    buf2[0] + puShiftXIndex + puShiftYIndex * buf2Stride[0];
                buf1[1] =
                    buf1[1] + puShiftXIndex + puShiftYIndex * buf1Stride[1];
                buf2[1] =
                    buf2[1] + puShiftXIndex + puShiftYIndex * buf2Stride[1];
                buf1[2] =
                    buf1[2] + puShiftXIndex + puShiftYIndex * buf1Stride[2];
                buf2[2] =
                    buf2[2] + puShiftXIndex + puShiftYIndex * buf2Stride[2];
                buf1[3] =
                    buf1[3] + puShiftXIndex + puShiftYIndex * buf1Stride[3];
                buf2[3] =
                    buf2[3] + puShiftXIndex + puShiftYIndex * buf2Stride[3];
                buf1[4] =
                    buf1[4] + puShiftXIndex + puShiftYIndex * buf1Stride[4];
                buf2[4] =
                    buf2[4] + puShiftXIndex + puShiftYIndex * buf2Stride[4];
                buf1[5] =
                    buf1[5] + puShiftXIndex + puShiftYIndex * buf1Stride[5];
                buf2[5] =
                    buf2[5] + puShiftXIndex + puShiftYIndex * buf2Stride[5];
                buf1[6] =
                    buf1[6] + puShiftXIndex + puShiftYIndex * buf1Stride[6];
                buf2[6] =
                    buf2[6] + puShiftXIndex + puShiftYIndex * buf2Stride[6];
                buf1[7] =
                    buf1[7] + puShiftXIndex + puShiftYIndex * buf1Stride[7];
                buf2[7] =
                    buf2[7] + puShiftXIndex + puShiftYIndex * buf2Stride[7];

                PU_QuarterPelRefinementOnTheFly(
                    context_ptr,
                    &context_ptr->p_best_ssd8x8[nidx],
                    puLcuBufferIndex,
                    buf1,
                    buf1Stride,
                    buf2,
                    buf2Stride,
                    8,
                    8,
                    x_search_area_origin,
                    y_search_area_origin,
                    &context_ptr->p_best_sad8x8[nidx],
                    &context_ptr->p_best_mv8x8[nidx],
                    context_ptr->psub_pel_direction8x8[nidx]);
            }
        }
    }

    if (ext_block_flag) {
        // 64x32
        for (pu_index = 0; pu_index < 2; ++pu_index) {
            puShiftXIndex = 0;
            puShiftYIndex = pu_index << 5;

            x_mv = _MVXT(context_ptr->p_best_mv64x32[pu_index]);
            y_mv = _MVYT(context_ptr->p_best_mv64x32[pu_index]);

            SetQuarterPelRefinementInputsOnTheFly(
                pos_Full,
                FullStride,
                pos_b,
                pos_h,
                pos_j,
                context_ptr->interpolated_stride,
                x_mv,
                y_mv,
                buf1,
                buf1Stride,
                buf2,
                buf2Stride);

            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + puShiftXIndex + puShiftYIndex * buf1Stride[0];
            buf2[0] = buf2[0] + puShiftXIndex + puShiftYIndex * buf2Stride[0];
            buf1[1] = buf1[1] + puShiftXIndex + puShiftYIndex * buf1Stride[1];
            buf2[1] = buf2[1] + puShiftXIndex + puShiftYIndex * buf2Stride[1];
            buf1[2] = buf1[2] + puShiftXIndex + puShiftYIndex * buf1Stride[2];
            buf2[2] = buf2[2] + puShiftXIndex + puShiftYIndex * buf2Stride[2];
            buf1[3] = buf1[3] + puShiftXIndex + puShiftYIndex * buf1Stride[3];
            buf2[3] = buf2[3] + puShiftXIndex + puShiftYIndex * buf2Stride[3];
            buf1[4] = buf1[4] + puShiftXIndex + puShiftYIndex * buf1Stride[4];
            buf2[4] = buf2[4] + puShiftXIndex + puShiftYIndex * buf2Stride[4];
            buf1[5] = buf1[5] + puShiftXIndex + puShiftYIndex * buf1Stride[5];
            buf2[5] = buf2[5] + puShiftXIndex + puShiftYIndex * buf2Stride[5];
            buf1[6] = buf1[6] + puShiftXIndex + puShiftYIndex * buf1Stride[6];
            buf2[6] = buf2[6] + puShiftXIndex + puShiftYIndex * buf2Stride[6];
            buf1[7] = buf1[7] + puShiftXIndex + puShiftYIndex * buf1Stride[7];
            buf2[7] = buf2[7] + puShiftXIndex + puShiftYIndex * buf2Stride[7];

            PU_QuarterPelRefinementOnTheFly(
                context_ptr,
                &context_ptr->p_best_ssd64x32[pu_index],
                puLcuBufferIndex,
                buf1,
                buf1Stride,
                buf2,
                buf2Stride,
                64,
                32,
                x_search_area_origin,
                y_search_area_origin,
                &context_ptr->p_best_sad64x32[pu_index],
                &context_ptr->p_best_mv64x32[pu_index],
                context_ptr->psub_pel_direction64x32[pu_index]);
        }

        // 32x16
        for (pu_index = 0; pu_index < 8; ++pu_index) {
            nidx = tab32x16[pu_index];  // TODO bitwise this

            puShiftXIndex = (pu_index & 0x01) << 5;
            puShiftYIndex = (pu_index >> 1) << 4;

            x_mv = _MVXT(context_ptr->p_best_mv32x16[nidx]);
            y_mv = _MVYT(context_ptr->p_best_mv32x16[nidx]);

            SetQuarterPelRefinementInputsOnTheFly(
                pos_Full,
                FullStride,
                pos_b,
                pos_h,
                pos_j,
                context_ptr->interpolated_stride,
                x_mv,
                y_mv,
                buf1,
                buf1Stride,
                buf2,
                buf2Stride);

            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + puShiftXIndex + puShiftYIndex * buf1Stride[0];
            buf2[0] = buf2[0] + puShiftXIndex + puShiftYIndex * buf2Stride[0];
            buf1[1] = buf1[1] + puShiftXIndex + puShiftYIndex * buf1Stride[1];
            buf2[1] = buf2[1] + puShiftXIndex + puShiftYIndex * buf2Stride[1];
            buf1[2] = buf1[2] + puShiftXIndex + puShiftYIndex * buf1Stride[2];
            buf2[2] = buf2[2] + puShiftXIndex + puShiftYIndex * buf2Stride[2];
            buf1[3] = buf1[3] + puShiftXIndex + puShiftYIndex * buf1Stride[3];
            buf2[3] = buf2[3] + puShiftXIndex + puShiftYIndex * buf2Stride[3];
            buf1[4] = buf1[4] + puShiftXIndex + puShiftYIndex * buf1Stride[4];
            buf2[4] = buf2[4] + puShiftXIndex + puShiftYIndex * buf2Stride[4];
            buf1[5] = buf1[5] + puShiftXIndex + puShiftYIndex * buf1Stride[5];
            buf2[5] = buf2[5] + puShiftXIndex + puShiftYIndex * buf2Stride[5];
            buf1[6] = buf1[6] + puShiftXIndex + puShiftYIndex * buf1Stride[6];
            buf2[6] = buf2[6] + puShiftXIndex + puShiftYIndex * buf2Stride[6];
            buf1[7] = buf1[7] + puShiftXIndex + puShiftYIndex * buf1Stride[7];
            buf2[7] = buf2[7] + puShiftXIndex + puShiftYIndex * buf2Stride[7];

            PU_QuarterPelRefinementOnTheFly(
                context_ptr,
                &context_ptr->p_best_ssd32x16[nidx],
                puLcuBufferIndex,
                buf1,
                buf1Stride,
                buf2,
                buf2Stride,
                32,
                16,
                x_search_area_origin,
                y_search_area_origin,
                &context_ptr->p_best_sad32x16[nidx],
                &context_ptr->p_best_mv32x16[nidx],
                context_ptr->psub_pel_direction32x16[nidx]);
        }

        // 16x8
        for (pu_index = 0; pu_index < 32; ++pu_index) {
            nidx = tab16x8[pu_index];

            puShiftXIndex = (pu_index & 0x03) << 4;
            puShiftYIndex = (pu_index >> 2) << 3;

            x_mv = _MVXT(context_ptr->p_best_mv16x8[nidx]);
            y_mv = _MVYT(context_ptr->p_best_mv16x8[nidx]);

            SetQuarterPelRefinementInputsOnTheFly(
                pos_Full,
                FullStride,
                pos_b,
                pos_h,
                pos_j,
                context_ptr->interpolated_stride,
                x_mv,
                y_mv,
                buf1,
                buf1Stride,
                buf2,
                buf2Stride);

            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + puShiftXIndex + puShiftYIndex * buf1Stride[0];
            buf2[0] = buf2[0] + puShiftXIndex + puShiftYIndex * buf2Stride[0];
            buf1[1] = buf1[1] + puShiftXIndex + puShiftYIndex * buf1Stride[1];
            buf2[1] = buf2[1] + puShiftXIndex + puShiftYIndex * buf2Stride[1];
            buf1[2] = buf1[2] + puShiftXIndex + puShiftYIndex * buf1Stride[2];
            buf2[2] = buf2[2] + puShiftXIndex + puShiftYIndex * buf2Stride[2];
            buf1[3] = buf1[3] + puShiftXIndex + puShiftYIndex * buf1Stride[3];
            buf2[3] = buf2[3] + puShiftXIndex + puShiftYIndex * buf2Stride[3];
            buf1[4] = buf1[4] + puShiftXIndex + puShiftYIndex * buf1Stride[4];
            buf2[4] = buf2[4] + puShiftXIndex + puShiftYIndex * buf2Stride[4];
            buf1[5] = buf1[5] + puShiftXIndex + puShiftYIndex * buf1Stride[5];
            buf2[5] = buf2[5] + puShiftXIndex + puShiftYIndex * buf2Stride[5];
            buf1[6] = buf1[6] + puShiftXIndex + puShiftYIndex * buf1Stride[6];
            buf2[6] = buf2[6] + puShiftXIndex + puShiftYIndex * buf2Stride[6];
            buf1[7] = buf1[7] + puShiftXIndex + puShiftYIndex * buf1Stride[7];
            buf2[7] = buf2[7] + puShiftXIndex + puShiftYIndex * buf2Stride[7];

            PU_QuarterPelRefinementOnTheFly(
                context_ptr,
                &context_ptr->p_best_ssd16x8[nidx],
                puLcuBufferIndex,
                buf1,
                buf1Stride,
                buf2,
                buf2Stride,
                16,
                8,
                x_search_area_origin,
                y_search_area_origin,
                &context_ptr->p_best_sad16x8[nidx],
                &context_ptr->p_best_mv16x8[nidx],
                context_ptr->psub_pel_direction16x8[nidx]);
        }

        // 32x64
        for (pu_index = 0; pu_index < 2; ++pu_index) {
            puShiftXIndex = pu_index << 5;
            puShiftYIndex = 0;

            x_mv = _MVXT(context_ptr->p_best_mv32x64[pu_index]);
            y_mv = _MVYT(context_ptr->p_best_mv32x64[pu_index]);

            SetQuarterPelRefinementInputsOnTheFly(
                pos_Full,
                FullStride,
                pos_b,
                pos_h,
                pos_j,
                context_ptr->interpolated_stride,
                x_mv,
                y_mv,
                buf1,
                buf1Stride,
                buf2,
                buf2Stride);

            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + puShiftXIndex + puShiftYIndex * buf1Stride[0];
            buf2[0] = buf2[0] + puShiftXIndex + puShiftYIndex * buf2Stride[0];
            buf1[1] = buf1[1] + puShiftXIndex + puShiftYIndex * buf1Stride[1];
            buf2[1] = buf2[1] + puShiftXIndex + puShiftYIndex * buf2Stride[1];
            buf1[2] = buf1[2] + puShiftXIndex + puShiftYIndex * buf1Stride[2];
            buf2[2] = buf2[2] + puShiftXIndex + puShiftYIndex * buf2Stride[2];
            buf1[3] = buf1[3] + puShiftXIndex + puShiftYIndex * buf1Stride[3];
            buf2[3] = buf2[3] + puShiftXIndex + puShiftYIndex * buf2Stride[3];
            buf1[4] = buf1[4] + puShiftXIndex + puShiftYIndex * buf1Stride[4];
            buf2[4] = buf2[4] + puShiftXIndex + puShiftYIndex * buf2Stride[4];
            buf1[5] = buf1[5] + puShiftXIndex + puShiftYIndex * buf1Stride[5];
            buf2[5] = buf2[5] + puShiftXIndex + puShiftYIndex * buf2Stride[5];
            buf1[6] = buf1[6] + puShiftXIndex + puShiftYIndex * buf1Stride[6];
            buf2[6] = buf2[6] + puShiftXIndex + puShiftYIndex * buf2Stride[6];
            buf1[7] = buf1[7] + puShiftXIndex + puShiftYIndex * buf1Stride[7];
            buf2[7] = buf2[7] + puShiftXIndex + puShiftYIndex * buf2Stride[7];

            PU_QuarterPelRefinementOnTheFly(
                context_ptr,
                &context_ptr->p_best_ssd32x64[pu_index],
                puLcuBufferIndex,
                buf1,
                buf1Stride,
                buf2,
                buf2Stride,
                32,
                64,
                x_search_area_origin,
                y_search_area_origin,
                &context_ptr->p_best_sad32x64[pu_index],
                &context_ptr->p_best_mv32x64[pu_index],
                context_ptr->psub_pel_direction32x64[pu_index]);
        }

        // 16x32
        for (pu_index = 0; pu_index < 8; ++pu_index) {
            nidx = tab16x32[pu_index];

            puShiftXIndex = (pu_index & 0x03) << 4;
            puShiftYIndex = (pu_index >> 2) << 5;

            x_mv = _MVXT(context_ptr->p_best_mv16x32[nidx]);
            y_mv = _MVYT(context_ptr->p_best_mv16x32[nidx]);

            SetQuarterPelRefinementInputsOnTheFly(
                pos_Full,
                FullStride,
                pos_b,
                pos_h,
                pos_j,
                context_ptr->interpolated_stride,
                x_mv,
                y_mv,
                buf1,
                buf1Stride,
                buf2,
                buf2Stride);

            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + puShiftXIndex + puShiftYIndex * buf1Stride[0];
            buf2[0] = buf2[0] + puShiftXIndex + puShiftYIndex * buf2Stride[0];
            buf1[1] = buf1[1] + puShiftXIndex + puShiftYIndex * buf1Stride[1];
            buf2[1] = buf2[1] + puShiftXIndex + puShiftYIndex * buf2Stride[1];
            buf1[2] = buf1[2] + puShiftXIndex + puShiftYIndex * buf1Stride[2];
            buf2[2] = buf2[2] + puShiftXIndex + puShiftYIndex * buf2Stride[2];
            buf1[3] = buf1[3] + puShiftXIndex + puShiftYIndex * buf1Stride[3];
            buf2[3] = buf2[3] + puShiftXIndex + puShiftYIndex * buf2Stride[3];
            buf1[4] = buf1[4] + puShiftXIndex + puShiftYIndex * buf1Stride[4];
            buf2[4] = buf2[4] + puShiftXIndex + puShiftYIndex * buf2Stride[4];
            buf1[5] = buf1[5] + puShiftXIndex + puShiftYIndex * buf1Stride[5];
            buf2[5] = buf2[5] + puShiftXIndex + puShiftYIndex * buf2Stride[5];
            buf1[6] = buf1[6] + puShiftXIndex + puShiftYIndex * buf1Stride[6];
            buf2[6] = buf2[6] + puShiftXIndex + puShiftYIndex * buf2Stride[6];
            buf1[7] = buf1[7] + puShiftXIndex + puShiftYIndex * buf1Stride[7];
            buf2[7] = buf2[7] + puShiftXIndex + puShiftYIndex * buf2Stride[7];

            PU_QuarterPelRefinementOnTheFly(
                context_ptr,
                &context_ptr->p_best_ssd16x32[nidx],
                puLcuBufferIndex,
                buf1,
                buf1Stride,
                buf2,
                buf2Stride,
                16,
                32,
                x_search_area_origin,
                y_search_area_origin,
                &context_ptr->p_best_sad16x32[nidx],
                &context_ptr->p_best_mv16x32[nidx],
                context_ptr->psub_pel_direction16x32[nidx]);
        }

        // 8x16
        for (pu_index = 0; pu_index < 32; ++pu_index) {
            nidx = tab8x16[pu_index];

            puShiftXIndex = (pu_index & 0x07) << 3;
            puShiftYIndex = (pu_index >> 3) << 4;

            x_mv = _MVXT(context_ptr->p_best_mv8x16[nidx]);
            y_mv = _MVYT(context_ptr->p_best_mv8x16[nidx]);

            SetQuarterPelRefinementInputsOnTheFly(
                pos_Full,
                FullStride,
                pos_b,
                pos_h,
                pos_j,
                context_ptr->interpolated_stride,
                x_mv,
                y_mv,
                buf1,
                buf1Stride,
                buf2,
                buf2Stride);

            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + puShiftXIndex + puShiftYIndex * buf1Stride[0];
            buf2[0] = buf2[0] + puShiftXIndex + puShiftYIndex * buf2Stride[0];
            buf1[1] = buf1[1] + puShiftXIndex + puShiftYIndex * buf1Stride[1];
            buf2[1] = buf2[1] + puShiftXIndex + puShiftYIndex * buf2Stride[1];
            buf1[2] = buf1[2] + puShiftXIndex + puShiftYIndex * buf1Stride[2];
            buf2[2] = buf2[2] + puShiftXIndex + puShiftYIndex * buf2Stride[2];
            buf1[3] = buf1[3] + puShiftXIndex + puShiftYIndex * buf1Stride[3];
            buf2[3] = buf2[3] + puShiftXIndex + puShiftYIndex * buf2Stride[3];
            buf1[4] = buf1[4] + puShiftXIndex + puShiftYIndex * buf1Stride[4];
            buf2[4] = buf2[4] + puShiftXIndex + puShiftYIndex * buf2Stride[4];
            buf1[5] = buf1[5] + puShiftXIndex + puShiftYIndex * buf1Stride[5];
            buf2[5] = buf2[5] + puShiftXIndex + puShiftYIndex * buf2Stride[5];
            buf1[6] = buf1[6] + puShiftXIndex + puShiftYIndex * buf1Stride[6];
            buf2[6] = buf2[6] + puShiftXIndex + puShiftYIndex * buf2Stride[6];
            buf1[7] = buf1[7] + puShiftXIndex + puShiftYIndex * buf1Stride[7];
            buf2[7] = buf2[7] + puShiftXIndex + puShiftYIndex * buf2Stride[7];

            PU_QuarterPelRefinementOnTheFly(
                context_ptr,
                &context_ptr->p_best_ssd8x16[nidx],
                puLcuBufferIndex,
                buf1,
                buf1Stride,
                buf2,
                buf2Stride,
                8,
                16,
                x_search_area_origin,
                y_search_area_origin,
                &context_ptr->p_best_sad8x16[nidx],
                &context_ptr->p_best_mv8x16[nidx],
                context_ptr->psub_pel_direction8x16[nidx]);
        }

        // 32x8
        for (pu_index = 0; pu_index < 16; ++pu_index) {
            nidx = tab32x8[pu_index];

            puShiftXIndex = (pu_index & 0x01) << 5;
            puShiftYIndex = (pu_index >> 1) << 3;

            x_mv = _MVXT(context_ptr->p_best_mv32x8[nidx]);
            y_mv = _MVYT(context_ptr->p_best_mv32x8[nidx]);

            SetQuarterPelRefinementInputsOnTheFly(
                pos_Full,
                FullStride,
                pos_b,
                pos_h,
                pos_j,
                context_ptr->interpolated_stride,
                x_mv,
                y_mv,
                buf1,
                buf1Stride,
                buf2,
                buf2Stride);

            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + puShiftXIndex + puShiftYIndex * buf1Stride[0];
            buf2[0] = buf2[0] + puShiftXIndex + puShiftYIndex * buf2Stride[0];
            buf1[1] = buf1[1] + puShiftXIndex + puShiftYIndex * buf1Stride[1];
            buf2[1] = buf2[1] + puShiftXIndex + puShiftYIndex * buf2Stride[1];
            buf1[2] = buf1[2] + puShiftXIndex + puShiftYIndex * buf1Stride[2];
            buf2[2] = buf2[2] + puShiftXIndex + puShiftYIndex * buf2Stride[2];
            buf1[3] = buf1[3] + puShiftXIndex + puShiftYIndex * buf1Stride[3];
            buf2[3] = buf2[3] + puShiftXIndex + puShiftYIndex * buf2Stride[3];
            buf1[4] = buf1[4] + puShiftXIndex + puShiftYIndex * buf1Stride[4];
            buf2[4] = buf2[4] + puShiftXIndex + puShiftYIndex * buf2Stride[4];
            buf1[5] = buf1[5] + puShiftXIndex + puShiftYIndex * buf1Stride[5];
            buf2[5] = buf2[5] + puShiftXIndex + puShiftYIndex * buf2Stride[5];
            buf1[6] = buf1[6] + puShiftXIndex + puShiftYIndex * buf1Stride[6];
            buf2[6] = buf2[6] + puShiftXIndex + puShiftYIndex * buf2Stride[6];
            buf1[7] = buf1[7] + puShiftXIndex + puShiftYIndex * buf1Stride[7];
            buf2[7] = buf2[7] + puShiftXIndex + puShiftYIndex * buf2Stride[7];

            PU_QuarterPelRefinementOnTheFly(
                context_ptr,
                &context_ptr->p_best_ssd32x8[nidx],
                puLcuBufferIndex,
                buf1,
                buf1Stride,
                buf2,
                buf2Stride,
                32,
                8,
                x_search_area_origin,
                y_search_area_origin,
                &context_ptr->p_best_sad32x8[nidx],
                &context_ptr->p_best_mv32x8[nidx],
                context_ptr->psub_pel_direction32x8[nidx]);
        }

        // 8x32
        for (pu_index = 0; pu_index < 16; ++pu_index) {
            nidx = tab8x32[pu_index];

            puShiftXIndex = (pu_index & 0x07) << 3;
            puShiftYIndex = (pu_index >> 3) << 5;

            x_mv = _MVXT(context_ptr->p_best_mv8x32[nidx]);
            y_mv = _MVYT(context_ptr->p_best_mv8x32[nidx]);

            SetQuarterPelRefinementInputsOnTheFly(
                pos_Full,
                FullStride,
                pos_b,
                pos_h,
                pos_j,
                context_ptr->interpolated_stride,
                x_mv,
                y_mv,
                buf1,
                buf1Stride,
                buf2,
                buf2Stride);

            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + puShiftXIndex + puShiftYIndex * buf1Stride[0];
            buf2[0] = buf2[0] + puShiftXIndex + puShiftYIndex * buf2Stride[0];
            buf1[1] = buf1[1] + puShiftXIndex + puShiftYIndex * buf1Stride[1];
            buf2[1] = buf2[1] + puShiftXIndex + puShiftYIndex * buf2Stride[1];
            buf1[2] = buf1[2] + puShiftXIndex + puShiftYIndex * buf1Stride[2];
            buf2[2] = buf2[2] + puShiftXIndex + puShiftYIndex * buf2Stride[2];
            buf1[3] = buf1[3] + puShiftXIndex + puShiftYIndex * buf1Stride[3];
            buf2[3] = buf2[3] + puShiftXIndex + puShiftYIndex * buf2Stride[3];
            buf1[4] = buf1[4] + puShiftXIndex + puShiftYIndex * buf1Stride[4];
            buf2[4] = buf2[4] + puShiftXIndex + puShiftYIndex * buf2Stride[4];
            buf1[5] = buf1[5] + puShiftXIndex + puShiftYIndex * buf1Stride[5];
            buf2[5] = buf2[5] + puShiftXIndex + puShiftYIndex * buf2Stride[5];
            buf1[6] = buf1[6] + puShiftXIndex + puShiftYIndex * buf1Stride[6];
            buf2[6] = buf2[6] + puShiftXIndex + puShiftYIndex * buf2Stride[6];
            buf1[7] = buf1[7] + puShiftXIndex + puShiftYIndex * buf1Stride[7];
            buf2[7] = buf2[7] + puShiftXIndex + puShiftYIndex * buf2Stride[7];

            PU_QuarterPelRefinementOnTheFly(
                context_ptr,
                &context_ptr->p_best_ssd8x32[nidx],
                puLcuBufferIndex,
                buf1,
                buf1Stride,
                buf2,
                buf2Stride,
                8,
                32,
                x_search_area_origin,
                y_search_area_origin,
                &context_ptr->p_best_sad8x32[nidx],
                &context_ptr->p_best_mv8x32[nidx],
                context_ptr->psub_pel_direction8x32[nidx]);
        }

        // 64x16
        for (pu_index = 0; pu_index < 4; ++pu_index) {
            nidx = pu_index;

            puShiftXIndex = 0;
            puShiftYIndex = pu_index << 4;

            x_mv = _MVXT(context_ptr->p_best_mv64x16[nidx]);
            y_mv = _MVYT(context_ptr->p_best_mv64x16[nidx]);

            SetQuarterPelRefinementInputsOnTheFly(
                pos_Full,
                FullStride,
                pos_b,
                pos_h,
                pos_j,
                context_ptr->interpolated_stride,
                x_mv,
                y_mv,
                buf1,
                buf1Stride,
                buf2,
                buf2Stride);

            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + puShiftXIndex + puShiftYIndex * buf1Stride[0];
            buf2[0] = buf2[0] + puShiftXIndex + puShiftYIndex * buf2Stride[0];
            buf1[1] = buf1[1] + puShiftXIndex + puShiftYIndex * buf1Stride[1];
            buf2[1] = buf2[1] + puShiftXIndex + puShiftYIndex * buf2Stride[1];
            buf1[2] = buf1[2] + puShiftXIndex + puShiftYIndex * buf1Stride[2];
            buf2[2] = buf2[2] + puShiftXIndex + puShiftYIndex * buf2Stride[2];
            buf1[3] = buf1[3] + puShiftXIndex + puShiftYIndex * buf1Stride[3];
            buf2[3] = buf2[3] + puShiftXIndex + puShiftYIndex * buf2Stride[3];
            buf1[4] = buf1[4] + puShiftXIndex + puShiftYIndex * buf1Stride[4];
            buf2[4] = buf2[4] + puShiftXIndex + puShiftYIndex * buf2Stride[4];
            buf1[5] = buf1[5] + puShiftXIndex + puShiftYIndex * buf1Stride[5];
            buf2[5] = buf2[5] + puShiftXIndex + puShiftYIndex * buf2Stride[5];
            buf1[6] = buf1[6] + puShiftXIndex + puShiftYIndex * buf1Stride[6];
            buf2[6] = buf2[6] + puShiftXIndex + puShiftYIndex * buf2Stride[6];
            buf1[7] = buf1[7] + puShiftXIndex + puShiftYIndex * buf1Stride[7];
            buf2[7] = buf2[7] + puShiftXIndex + puShiftYIndex * buf2Stride[7];

            PU_QuarterPelRefinementOnTheFly(
                context_ptr,
                &context_ptr->p_best_ssd64x16[nidx],
                puLcuBufferIndex,
                buf1,
                buf1Stride,
                buf2,
                buf2Stride,
                64,
                16,
                x_search_area_origin,
                y_search_area_origin,
                &context_ptr->p_best_sad64x16[nidx],
                &context_ptr->p_best_mv64x16[nidx],
                context_ptr->psub_pel_direction64x16[nidx]);
        }

        // 16x64
        for (pu_index = 0; pu_index < 4; ++pu_index) {
            nidx = pu_index;

            puShiftXIndex = pu_index << 4;
            puShiftYIndex = 0;

            x_mv = _MVXT(context_ptr->p_best_mv16x64[nidx]);
            y_mv = _MVYT(context_ptr->p_best_mv16x64[nidx]);

            SetQuarterPelRefinementInputsOnTheFly(
                pos_Full,
                FullStride,
                pos_b,
                pos_h,
                pos_j,
                context_ptr->interpolated_stride,
                x_mv,
                y_mv,
                buf1,
                buf1Stride,
                buf2,
                buf2Stride);

            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + puShiftXIndex + puShiftYIndex * buf1Stride[0];
            buf2[0] = buf2[0] + puShiftXIndex + puShiftYIndex * buf2Stride[0];
            buf1[1] = buf1[1] + puShiftXIndex + puShiftYIndex * buf1Stride[1];
            buf2[1] = buf2[1] + puShiftXIndex + puShiftYIndex * buf2Stride[1];
            buf1[2] = buf1[2] + puShiftXIndex + puShiftYIndex * buf1Stride[2];
            buf2[2] = buf2[2] + puShiftXIndex + puShiftYIndex * buf2Stride[2];
            buf1[3] = buf1[3] + puShiftXIndex + puShiftYIndex * buf1Stride[3];
            buf2[3] = buf2[3] + puShiftXIndex + puShiftYIndex * buf2Stride[3];
            buf1[4] = buf1[4] + puShiftXIndex + puShiftYIndex * buf1Stride[4];
            buf2[4] = buf2[4] + puShiftXIndex + puShiftYIndex * buf2Stride[4];
            buf1[5] = buf1[5] + puShiftXIndex + puShiftYIndex * buf1Stride[5];
            buf2[5] = buf2[5] + puShiftXIndex + puShiftYIndex * buf2Stride[5];
            buf1[6] = buf1[6] + puShiftXIndex + puShiftYIndex * buf1Stride[6];
            buf2[6] = buf2[6] + puShiftXIndex + puShiftYIndex * buf2Stride[6];
            buf1[7] = buf1[7] + puShiftXIndex + puShiftYIndex * buf1Stride[7];
            buf2[7] = buf2[7] + puShiftXIndex + puShiftYIndex * buf2Stride[7];

            PU_QuarterPelRefinementOnTheFly(
                context_ptr,
                &context_ptr->p_best_ssd16x64[nidx],
                puLcuBufferIndex,
                buf1,
                buf1Stride,
                buf2,
                buf2Stride,
                16,
                64,
                x_search_area_origin,
                y_search_area_origin,
                &context_ptr->p_best_sad16x64[nidx],
                &context_ptr->p_best_mv16x64[nidx],
                context_ptr->psub_pel_direction16x64[nidx]);
        }
    }

    return;
}
#define QP_REF_OPT 1
/*******************************************
 * quarter_pel_refinemnet_block
 *   performs Quarter Pel refinement for each block
 *******************************************/
static void quarter_pel_refinemnet_block(
    MeContext *context_ptr,  // [IN] ME context Ptr, used to get SB Ptr
    uint32_t *p_best_ssd,
    uint32_t
        src_block_index,  // [IN] PU origin, used to point to source samples
    uint8_t **buf1,       // [IN]
    uint32_t *buf1_stride,
    uint8_t **buf2,  // [IN]
    uint32_t *buf2_stride,
    uint32_t pu_width,   // [IN]  PU width
    uint32_t pu_height,  // [IN]  PU height
    int16_t
        x_search_area_origin,  // [IN] search area origin in the horizontal
                               // direction, used to point to reference samples
    int16_t
        y_search_area_origin,  // [IN] search area origin in the vertical
                               // direction, used to point to reference samples
    uint32_t candidate_mv, uint32_t *p_best_sad,
    uint32_t *p_best_mv, uint16_t is_frac_candidate) {
    int16_t x_mv = _MVXT(candidate_mv);
    int16_t y_mv = _MVYT(candidate_mv);
    int16_t search_Index_x = ((x_mv + 2) >> 2) - x_search_area_origin;
    int16_t search_Index_y = ((y_mv + 2) >> 2) - y_search_area_origin;
    uint64_t dist;
    int16_t quarter_mv_x[8];
    int16_t quarter_mv_y[8];
    int32_t search_region_Index1 = 0;
    int32_t search_region_Index2 = 0;
    quarter_mv_x[0] = x_mv - 1;  // L  position
    quarter_mv_x[1] = x_mv + 1;  // R  position
    quarter_mv_x[2] = x_mv;      // T  position
    quarter_mv_x[3] = x_mv;      // B  position
    quarter_mv_x[4] = x_mv - 1;  // TL position
    quarter_mv_x[5] = x_mv + 1;  // TR position
    quarter_mv_x[6] = x_mv + 1;  // BR position
    quarter_mv_x[7] = x_mv - 1;  // BL position
    quarter_mv_y[0] = y_mv;      // L  position
    quarter_mv_y[1] = y_mv;      // R  position
    quarter_mv_y[2] = y_mv - 1;  // T  position
    quarter_mv_y[3] = y_mv + 1;  // B  position
    quarter_mv_y[4] = y_mv - 1;  // TL position
    quarter_mv_y[5] = y_mv - 1;  // TR position
    quarter_mv_y[6] = y_mv + 1;  // BR position
    quarter_mv_y[7] = y_mv + 1;  // BL position
    // L position
    search_region_Index1 = (int32_t)search_Index_x +
                           (int32_t)buf1_stride[0] * (int32_t)search_Index_y;
    search_region_Index2 = (int32_t)search_Index_x +
                           (int32_t)buf2_stride[0] * (int32_t)search_Index_y;
    dist = (context_ptr->fractional_search_method == SSD_SEARCH)
               ? combined_averaging_ssd(
                     &(context_ptr->sb_buffer[src_block_index]),
                     BLOCK_SIZE_64,
                     buf1[0] + search_region_Index1,
                     buf1_stride[0],
                     buf2[0] + search_region_Index2,
                     buf2_stride[0],
                     pu_height,
                     pu_width)
               : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                     ? (nxm_sad_avg_kernel(
                                &(context_ptr->sb_buffer[src_block_index]),
                                BLOCK_SIZE_64 << 1,
                                buf1[0] + search_region_Index1,
                                buf1_stride[0] << 1,
                                buf2[0] + search_region_Index2,
                                buf2_stride[0] << 1,
                                pu_height >> 1,
                                pu_width))
                           << 1
                     : nxm_sad_avg_kernel(
                               &(context_ptr->sb_buffer[src_block_index]),
                               BLOCK_SIZE_64,
                               buf1[0] + search_region_Index1,
                               buf1_stride[0],
                               buf2[0] + search_region_Index2,
                               buf2_stride[0],
                               pu_height,
                               pu_width);
    if (context_ptr->fractional_search_method == SSD_SEARCH) {
        if (dist < *p_best_ssd) {
            *p_best_sad = (uint32_t)
                nxm_sad_avg_kernel(
                    &(context_ptr->sb_buffer[src_block_index]),
                    BLOCK_SIZE_64,
                    buf1[0] + search_region_Index1,
                    buf1_stride[0],
                    buf2[0] + search_region_Index2,
                    buf2_stride[0],
                    pu_height,
                    pu_width);
            *p_best_mv =
                ((uint16_t)quarter_mv_y[0] << 16) | ((uint16_t)quarter_mv_x[0]);
            *p_best_ssd = (uint32_t)dist;
        }
    } else {
        if (dist < *p_best_sad) {
            *p_best_sad = (uint32_t)dist;
            *p_best_mv =
                ((uint16_t)quarter_mv_y[0] << 16) | ((uint16_t)quarter_mv_x[0]);
        }
    }
    // R positions
    search_region_Index1 = (int32_t)search_Index_x +
                           (int32_t)buf1_stride[1] * (int32_t)search_Index_y;
    search_region_Index2 = (int32_t)search_Index_x +
                           (int32_t)buf2_stride[1] * (int32_t)search_Index_y;
    dist = (context_ptr->fractional_search_method == SSD_SEARCH)
               ? combined_averaging_ssd(
                     &(context_ptr->sb_buffer[src_block_index]),
                     BLOCK_SIZE_64,
                     buf1[1] + search_region_Index1,
                     buf1_stride[1],
                     buf2[1] + search_region_Index2,
                     buf2_stride[1],
                     pu_height,
                     pu_width)
               : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                     ? (nxm_sad_avg_kernel(
                                &(context_ptr->sb_buffer[src_block_index]),
                                BLOCK_SIZE_64 << 1,
                                buf1[1] + search_region_Index1,
                                buf1_stride[1] << 1,
                                buf2[1] + search_region_Index2,
                                buf2_stride[1] << 1,
                                pu_height >> 1,
                                pu_width))
                           << 1
                     : nxm_sad_avg_kernel(
                               &(context_ptr->sb_buffer[src_block_index]),
                               BLOCK_SIZE_64,
                               buf1[1] + search_region_Index1,
                               buf1_stride[1],
                               buf2[1] + search_region_Index2,
                               buf2_stride[1],
                               pu_height,
                               pu_width);
    if (context_ptr->fractional_search_method == SSD_SEARCH) {
        if (dist < *p_best_ssd) {
            *p_best_sad = (uint32_t)
                nxm_sad_avg_kernel(
                    &(context_ptr->sb_buffer[src_block_index]),
                    BLOCK_SIZE_64,
                    buf1[1] + search_region_Index1,
                    buf1_stride[1],
                    buf2[1] + search_region_Index2,
                    buf2_stride[1],
                    pu_height,
                    pu_width);
            *p_best_mv =
                ((uint16_t)quarter_mv_y[1] << 16) | ((uint16_t)quarter_mv_x[1]);
            *p_best_ssd = (uint32_t)dist;
        }
    } else {
        if (dist < *p_best_sad) {
            *p_best_sad = (uint32_t)dist;
            *p_best_mv =
                ((uint16_t)quarter_mv_y[1] << 16) | ((uint16_t)quarter_mv_x[1]);
        }
    }
    // T position
    search_region_Index1 = (int32_t)search_Index_x +
                           (int32_t)buf1_stride[2] * (int32_t)search_Index_y;
    search_region_Index2 = (int32_t)search_Index_x +
                           (int32_t)buf2_stride[2] * (int32_t)search_Index_y;
    dist = (context_ptr->fractional_search_method == SSD_SEARCH)
               ? combined_averaging_ssd(
                     &(context_ptr->sb_buffer[src_block_index]),
                     BLOCK_SIZE_64,
                     buf1[2] + search_region_Index1,
                     buf1_stride[2],
                     buf2[2] + search_region_Index2,
                     buf2_stride[2],
                     pu_height,
                     pu_width)
               : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                     ? (nxm_sad_avg_kernel(
                                &(context_ptr->sb_buffer[src_block_index]),
                                BLOCK_SIZE_64 << 1,
                                buf1[2] + search_region_Index1,
                                buf1_stride[2] << 1,
                                buf2[2] + search_region_Index2,
                                buf2_stride[2] << 1,
                                pu_height >> 1,
                                pu_width))
                           << 1
                     : nxm_sad_avg_kernel(
                               &(context_ptr->sb_buffer[src_block_index]),
                               BLOCK_SIZE_64,
                               buf1[2] + search_region_Index1,
                               buf1_stride[2],
                               buf2[2] + search_region_Index2,
                               buf2_stride[2],
                               pu_height,
                               pu_width);
    if (context_ptr->fractional_search_method == SSD_SEARCH) {
        if (dist < *p_best_ssd) {
            *p_best_sad = (uint32_t)
                nxm_sad_avg_kernel(
                    &(context_ptr->sb_buffer[src_block_index]),
                    BLOCK_SIZE_64,
                    buf1[2] + search_region_Index1,
                    buf1_stride[2],
                    buf2[2] + search_region_Index2,
                    buf2_stride[2],
                    pu_height,
                    pu_width);
            *p_best_mv =
                ((uint16_t)quarter_mv_y[2] << 16) | ((uint16_t)quarter_mv_x[2]);
            *p_best_ssd = (uint32_t)dist;
        }
    } else {
        if (dist < *p_best_sad) {
            *p_best_sad = (uint32_t)dist;
            *p_best_mv =
                ((uint16_t)quarter_mv_y[2] << 16) | ((uint16_t)quarter_mv_x[2]);
        }
    }
    // B position
    search_region_Index1 = (int32_t)search_Index_x +
                           (int32_t)buf1_stride[3] * (int32_t)search_Index_y;
    search_region_Index2 = (int32_t)search_Index_x +
                           (int32_t)buf2_stride[3] * (int32_t)search_Index_y;
    dist = (context_ptr->fractional_search_method == SSD_SEARCH)
               ? combined_averaging_ssd(
                     &(context_ptr->sb_buffer[src_block_index]),
                     BLOCK_SIZE_64,
                     buf1[3] + search_region_Index1,
                     buf1_stride[3],
                     buf2[3] + search_region_Index2,
                     buf2_stride[3],
                     pu_height,
                     pu_width)
               : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                     ? (nxm_sad_avg_kernel(
                                &(context_ptr->sb_buffer[src_block_index]),
                                BLOCK_SIZE_64 << 1,
                                buf1[3] + search_region_Index1,
                                buf1_stride[3] << 1,
                                buf2[3] + search_region_Index2,
                                buf2_stride[3] << 1,
                                pu_height >> 1,
                                pu_width))
                           << 1
                     : nxm_sad_avg_kernel(
                               &(context_ptr->sb_buffer[src_block_index]),
                               BLOCK_SIZE_64,
                               buf1[3] + search_region_Index1,
                               buf1_stride[3],
                               buf2[3] + search_region_Index2,
                               buf2_stride[3],
                               pu_height,
                               pu_width);
    if (context_ptr->fractional_search_method == SSD_SEARCH) {
        if (dist < *p_best_ssd) {
            *p_best_sad = (uint32_t)
                nxm_sad_avg_kernel(
                    &(context_ptr->sb_buffer[src_block_index]),
                    BLOCK_SIZE_64,
                    buf1[3] + search_region_Index1,
                    buf1_stride[3],
                    buf2[3] + search_region_Index2,
                    buf2_stride[3],
                    pu_height,
                    pu_width);
            *p_best_mv =
                ((uint16_t)quarter_mv_y[3] << 16) | ((uint16_t)quarter_mv_x[3]);
            *p_best_ssd = (uint32_t)dist;
        }
    } else {
        if (dist < *p_best_sad) {
            *p_best_sad = (uint32_t)dist;
            *p_best_mv =
                ((uint16_t)quarter_mv_y[3] << 16) | ((uint16_t)quarter_mv_x[3]);
        }
    }
    // TL position
#if QP_REF_OPT
    if (!is_frac_candidate) {
#endif
        search_region_Index1 =
            (int32_t)search_Index_x +
            (int32_t)buf1_stride[4] * (int32_t)search_Index_y;
        search_region_Index2 =
            (int32_t)search_Index_x +
            (int32_t)buf2_stride[4] * (int32_t)search_Index_y;
        dist = (context_ptr->fractional_search_method == SSD_SEARCH)
                   ? combined_averaging_ssd(
                         &(context_ptr->sb_buffer[src_block_index]),
                         BLOCK_SIZE_64,
                         buf1[4] + search_region_Index1,
                         buf1_stride[4],
                         buf2[4] + search_region_Index2,
                         buf2_stride[4],
                         pu_height,
                         pu_width)
                   : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                         ? (nxm_sad_avg_kernel(
                                    &(context_ptr->sb_buffer[src_block_index]),
                                    BLOCK_SIZE_64 << 1,
                                    buf1[4] + search_region_Index1,
                                    buf1_stride[4] << 1,
                                    buf2[4] + search_region_Index2,
                                    buf2_stride[4] << 1,
                                    pu_height >> 1,
                                    pu_width))
                               << 1
                         : nxm_sad_avg_kernel(
                                   &(context_ptr->sb_buffer[src_block_index]),
                                   BLOCK_SIZE_64,
                                   buf1[4] + search_region_Index1,
                                   buf1_stride[4],
                                   buf2[4] + search_region_Index2,
                                   buf2_stride[4],
                                   pu_height,
                                   pu_width);
        if (context_ptr->fractional_search_method == SSD_SEARCH) {
            if (dist < *p_best_ssd) {
                *p_best_sad = (uint32_t)nxm_sad_avg_kernel(
                        &(context_ptr->sb_buffer[src_block_index]),
                        BLOCK_SIZE_64,
                        buf1[4] + search_region_Index1,
                        buf1_stride[4],
                        buf2[4] + search_region_Index2,
                        buf2_stride[4],
                        pu_height,
                        pu_width);
                *p_best_mv = ((uint16_t)quarter_mv_y[4] << 16) |
                             ((uint16_t)quarter_mv_x[4]);
                *p_best_ssd = (uint32_t)dist;
            }
        } else {
            if (dist < *p_best_sad) {
                *p_best_sad = (uint32_t)dist;
                *p_best_mv = ((uint16_t)quarter_mv_y[4] << 16) |
                             ((uint16_t)quarter_mv_x[4]);
            }
        }
        // TR position
        search_region_Index1 =
            (int32_t)search_Index_x +
            (int32_t)buf1_stride[5] * (int32_t)search_Index_y;
        search_region_Index2 =
            (int32_t)search_Index_x +
            (int32_t)buf2_stride[5] * (int32_t)search_Index_y;
        dist = (context_ptr->fractional_search_method == SSD_SEARCH)
                   ? combined_averaging_ssd(
                         &(context_ptr->sb_buffer[src_block_index]),
                         BLOCK_SIZE_64,
                         buf1[5] + search_region_Index1,
                         buf1_stride[5],
                         buf2[5] + search_region_Index2,
                         buf2_stride[5],
                         pu_height,
                         pu_width)
                   : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                         ? (nxm_sad_avg_kernel(
                                    &(context_ptr->sb_buffer[src_block_index]),
                                    BLOCK_SIZE_64 << 1,
                                    buf1[5] + search_region_Index1,
                                    buf1_stride[5] << 1,
                                    buf2[5] + search_region_Index2,
                                    buf2_stride[5] << 1,
                                    pu_height >> 1,
                                    pu_width))
                               << 1
                         : nxm_sad_avg_kernel(
                                   &(context_ptr->sb_buffer[src_block_index]),
                                   BLOCK_SIZE_64,
                                   buf1[5] + search_region_Index1,
                                   buf1_stride[5],
                                   buf2[5] + search_region_Index2,
                                   buf2_stride[5],
                                   pu_height,
                                   pu_width);
        if (context_ptr->fractional_search_method == SSD_SEARCH) {
            if (dist < *p_best_ssd) {
                *p_best_sad = (uint32_t)nxm_sad_avg_kernel(
                        &(context_ptr->sb_buffer[src_block_index]),
                        BLOCK_SIZE_64,
                        buf1[5] + search_region_Index1,
                        buf1_stride[5],
                        buf2[5] + search_region_Index2,
                        buf2_stride[5],
                        pu_height,
                        pu_width);
                *p_best_mv = ((uint16_t)quarter_mv_y[5] << 16) |
                             ((uint16_t)quarter_mv_x[5]);
                *p_best_ssd = (uint32_t)dist;
            }
        } else {
            if (dist < *p_best_sad) {
                *p_best_sad = (uint32_t)dist;
                *p_best_mv = ((uint16_t)quarter_mv_y[5] << 16) |
                             ((uint16_t)quarter_mv_x[5]);
            }
        }
        // BR position
        search_region_Index1 =
            (int32_t)search_Index_x +
            (int32_t)buf1_stride[6] * (int32_t)search_Index_y;
        search_region_Index2 =
            (int32_t)search_Index_x +
            (int32_t)buf2_stride[6] * (int32_t)search_Index_y;
        dist = (context_ptr->fractional_search_method == SSD_SEARCH)
                   ? combined_averaging_ssd(
                         &(context_ptr->sb_buffer[src_block_index]),
                         BLOCK_SIZE_64,
                         buf1[6] + search_region_Index1,
                         buf1_stride[6],
                         buf2[6] + search_region_Index2,
                         buf2_stride[6],
                         pu_height,
                         pu_width)
                   : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                         ? (nxm_sad_avg_kernel(
                                    &(context_ptr->sb_buffer[src_block_index]),
                                    BLOCK_SIZE_64 << 1,
                                    buf1[6] + search_region_Index1,
                                    buf1_stride[6] << 1,
                                    buf2[6] + search_region_Index2,
                                    buf2_stride[6] << 1,
                                    pu_height >> 1,
                                    pu_width))
                               << 1
                         : nxm_sad_avg_kernel(
                                   &(context_ptr->sb_buffer[src_block_index]),
                                   BLOCK_SIZE_64,
                                   buf1[6] + search_region_Index1,
                                   buf1_stride[6],
                                   buf2[6] + search_region_Index2,
                                   buf2_stride[6],
                                   pu_height,
                                   pu_width);
        if (context_ptr->fractional_search_method == SSD_SEARCH) {
            if (dist < *p_best_ssd) {
                *p_best_sad = (uint32_t)nxm_sad_avg_kernel(
                        &(context_ptr->sb_buffer[src_block_index]),
                        BLOCK_SIZE_64,
                        buf1[6] + search_region_Index1,
                        buf1_stride[6],
                        buf2[6] + search_region_Index2,
                        buf2_stride[6],
                        pu_height,
                        pu_width);
                *p_best_mv = ((uint16_t)quarter_mv_y[6] << 16) |
                             ((uint16_t)quarter_mv_x[6]);
                *p_best_ssd = (uint32_t)dist;
            }
        } else {
            if (dist < *p_best_sad) {
                *p_best_sad = (uint32_t)dist;
                *p_best_mv = ((uint16_t)quarter_mv_y[6] << 16) |
                             ((uint16_t)quarter_mv_x[6]);
            }
        }
        // BL position
        search_region_Index1 =
            (int32_t)search_Index_x +
            (int32_t)buf1_stride[7] * (int32_t)search_Index_y;
        search_region_Index2 =
            (int32_t)search_Index_x +
            (int32_t)buf2_stride[7] * (int32_t)search_Index_y;
        dist = (context_ptr->fractional_search_method == SSD_SEARCH)
                   ? combined_averaging_ssd(
                         &(context_ptr->sb_buffer[src_block_index]),
                         BLOCK_SIZE_64,
                         buf1[7] + search_region_Index1,
                         buf1_stride[7],
                         buf2[7] + search_region_Index2,
                         buf2_stride[7],
                         pu_height,
                         pu_width)
                   : (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
                         ? (nxm_sad_avg_kernel(
                                    &(context_ptr->sb_buffer[src_block_index]),
                                    BLOCK_SIZE_64 << 1,
                                    buf1[7] + search_region_Index1,
                                    buf1_stride[7] << 1,
                                    buf2[7] + search_region_Index2,
                                    buf2_stride[7] << 1,
                                    pu_height >> 1,
                                    pu_width))
                               << 1
                         : nxm_sad_avg_kernel(
                                   &(context_ptr->sb_buffer[src_block_index]),
                                   BLOCK_SIZE_64,
                                   buf1[7] + search_region_Index1,
                                   buf1_stride[7],
                                   buf2[7] + search_region_Index2,
                                   buf2_stride[7],
                                   pu_height,
                                   pu_width);
        if (context_ptr->fractional_search_method == SSD_SEARCH) {
            if (dist < *p_best_ssd) {
                *p_best_sad = (uint32_t)nxm_sad_avg_kernel(
                        &(context_ptr->sb_buffer[src_block_index]),
                        BLOCK_SIZE_64,
                        buf1[7] + search_region_Index1,
                        buf1_stride[7],
                        buf2[7] + search_region_Index2,
                        buf2_stride[7],
                        pu_height,
                        pu_width);
                *p_best_mv = ((uint16_t)quarter_mv_y[7] << 16) |
                             ((uint16_t)quarter_mv_x[7]);
                *p_best_ssd = (uint32_t)dist;
            }
        } else {
            if (dist < *p_best_sad) {
                *p_best_sad = (uint32_t)dist;
                *p_best_mv = ((uint16_t)quarter_mv_y[7] << 16) |
                             ((uint16_t)quarter_mv_x[7]);
            }
        }
#if QP_REF_OPT
    }
#endif
    return;
}
/*******************************************
 * quarter_pel_refinement_sb
 *   performs Quarter Pel refinement
 *******************************************/
void quarter_pel_refinement_sb(
    MeContext
        *context_ptr,  //[IN/OUT]  ME context Ptr, used to get/update ME results
    uint8_t *pos_full,     //[IN]
    uint32_t full_stride,  //[IN]
    uint8_t *pos_b,        //[IN]
    uint8_t *pos_h,        //[IN]
    uint8_t *pos_j,        //[IN]
    int16_t
        x_search_area_origin,  //[IN] search area origin in the horizontal
                               // direction, used to point to reference samples
    int16_t
        y_search_area_origin,  //[IN] search area origin in the vertical
                               // direction, used to point to reference samples
    uint32_t integer_mv)
{
    uint32_t pu_index;
    uint32_t block_index_shift_x;
    uint32_t block_index_shift_y;
    uint32_t src_block_index;
    uint8_t *buf1[8];
    uint8_t *buf2[8];
    uint32_t buf1_stride[8];
    uint32_t buf2_stride[8];
    int16_t x_mv, y_mv;
    uint32_t nidx;
    int16_t int_x_mv = _MVXT(integer_mv);
    int16_t int_y_mv = _MVYT(integer_mv);
    int16_t int_xSearchIndex = ((int_x_mv + 2) >> 2) - x_search_area_origin;
    int16_t int_ySearchIndex = ((int_y_mv + 2) >> 2) - y_search_area_origin;
    int16_t x_best_mv;
    int16_t y_best_mv;
    int16_t best_xSearchIndex;
    int16_t best_ySearchIndex;
    int16_t dis_x;
    int16_t dis_y;
    int8_t skip_qp_pel = 0;
    uint32_t testmv;
    int16_t it;
    int16_t num_qp_it = 2;
    if (context_ptr->fractional_search64x64) {
        x_best_mv = _MVXT(*context_ptr->p_best_full_pel_mv64x64);
        y_best_mv = _MVYT(*context_ptr->p_best_full_pel_mv64x64);
        best_xSearchIndex = ((x_best_mv + 2) >> 2) - x_search_area_origin;
        best_ySearchIndex = ((y_best_mv + 2) >> 2) - y_search_area_origin;
        dis_x = ABS(int_xSearchIndex - best_xSearchIndex);
        dis_y = ABS(int_ySearchIndex - best_ySearchIndex);
        skip_qp_pel = 0;
        if ((dis_x) > Q_PEL_SEARCH_WIND)
            skip_qp_pel = 1;
        if ((dis_y) > Q_PEL_SEARCH_WIND)
            skip_qp_pel = 1;
        if (!skip_qp_pel) {
            for (it = 0; it < num_qp_it; it++) {
                x_mv = (int16_t)_MVXT(integer_mv) + (2 * it);
                y_mv = (int16_t)_MVYT(integer_mv) + (2 * it);
                testmv = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
                SetQuarterPelRefinementInputsOnTheFly(
                    pos_full,
                    full_stride,
                    pos_b,
                    pos_h,
                    pos_j,
                    context_ptr->interpolated_stride,
                    x_mv,
                    y_mv,
                    buf1,
                    buf1_stride,
                    buf2,
                    buf2_stride);
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
                quarter_pel_refinemnet_block(context_ptr,
                                             context_ptr->p_best_ssd64x64,
                                             0,
                                             buf1,
                                             buf1_stride,
                                             buf2,
                                             buf2_stride,
                                             64,
                                             64,
                                             x_search_area_origin,
                                             y_search_area_origin,
                                             testmv,
                                             context_ptr->p_best_sad64x64,
                                             context_ptr->p_best_mv64x64,
                                             it);
            }
        }
    }
    // 32x32 [4 partitions]
    for (pu_index = 0; pu_index < 4; ++pu_index) {
        x_best_mv = _MVXT(context_ptr->p_best_full_pel_mv32x32[pu_index]);
        y_best_mv = _MVYT(context_ptr->p_best_full_pel_mv32x32[pu_index]);
        best_xSearchIndex = ((x_best_mv + 2) >> 2) - x_search_area_origin;
        best_ySearchIndex = ((y_best_mv + 2) >> 2) - y_search_area_origin;
        dis_x = ABS(int_xSearchIndex - best_xSearchIndex);
        dis_y = ABS(int_ySearchIndex - best_ySearchIndex);
        skip_qp_pel = 0;
        if ((dis_x) > Q_PEL_SEARCH_WIND)
            skip_qp_pel = 1;
        if ((dis_y) > Q_PEL_SEARCH_WIND)
            skip_qp_pel = 1;
        if (!skip_qp_pel) {
            for (it = 0; it < num_qp_it; it++) {
                x_mv = (int16_t)_MVXT(integer_mv) + (2 * it);
                y_mv = (int16_t)_MVYT(integer_mv) + (2 * it);
                testmv = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
                SetQuarterPelRefinementInputsOnTheFly(
                    pos_full,
                    full_stride,
                    pos_b,
                    pos_h,
                    pos_j,
                    context_ptr->interpolated_stride,
                    x_mv,
                    y_mv,
                    buf1,
                    buf1_stride,
                    buf2,
                    buf2_stride);
                block_index_shift_x = (pu_index & 0x01) << 5;
                block_index_shift_y = (pu_index >> 1) << 5;
                src_block_index =
                    block_index_shift_x + block_index_shift_y * BLOCK_SIZE_64;
                buf1[0] = buf1[0] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[0];
                buf2[0] = buf2[0] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[0];
                buf1[1] = buf1[1] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[1];
                buf2[1] = buf2[1] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[1];
                buf1[2] = buf1[2] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[2];
                buf2[2] = buf2[2] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[2];
                buf1[3] = buf1[3] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[3];
                buf2[3] = buf2[3] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[3];
                buf1[4] = buf1[4] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[4];
                buf2[4] = buf2[4] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[4];
                buf1[5] = buf1[5] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[5];
                buf2[5] = buf2[5] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[5];
                buf1[6] = buf1[6] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[6];
                buf2[6] = buf2[6] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[6];
                buf1[7] = buf1[7] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[7];
                buf2[7] = buf2[7] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[7];
                quarter_pel_refinemnet_block(
                    context_ptr,
                    &context_ptr->p_best_ssd32x32[pu_index],
                    src_block_index,
                    buf1,
                    buf1_stride,
                    buf2,
                    buf2_stride,
                    32,
                    32,
                    x_search_area_origin,
                    y_search_area_origin,
                    testmv,
                    &context_ptr->p_best_sad32x32[pu_index],
                    &context_ptr->p_best_mv32x32[pu_index],
                    it);
            }
        }
    }
    // 16x16 [16 partitions]
    for (pu_index = 0; pu_index < 16; ++pu_index) {
        nidx = tab16x16[pu_index];
        x_best_mv = _MVXT(context_ptr->p_best_full_pel_mv16x16[nidx]);
        y_best_mv = _MVYT(context_ptr->p_best_full_pel_mv16x16[nidx]);
        best_xSearchIndex = ((x_best_mv + 2) >> 2) - x_search_area_origin;
        best_ySearchIndex = ((y_best_mv + 2) >> 2) - y_search_area_origin;
        dis_x = ABS(int_xSearchIndex - best_xSearchIndex);
        dis_y = ABS(int_ySearchIndex - best_ySearchIndex);
        skip_qp_pel = 0;
        if ((dis_x) > Q_PEL_SEARCH_WIND)
            skip_qp_pel = 1;
        if ((dis_y) > Q_PEL_SEARCH_WIND)
            skip_qp_pel = 1;
        if (!skip_qp_pel) {
            for (it = 0; it < num_qp_it; it++) {
                x_mv = (int16_t)_MVXT(integer_mv) + (2 * it);
                y_mv = (int16_t)_MVYT(integer_mv) + (2 * it);
                testmv = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
                SetQuarterPelRefinementInputsOnTheFly(
                    pos_full,
                    full_stride,
                    pos_b,
                    pos_h,
                    pos_j,
                    context_ptr->interpolated_stride,
                    x_mv,
                    y_mv,
                    buf1,
                    buf1_stride,
                    buf2,
                    buf2_stride);
                block_index_shift_x = (pu_index & 0x03) << 4;
                block_index_shift_y = (pu_index >> 2) << 4;
                src_block_index =
                    block_index_shift_x + block_index_shift_y * BLOCK_SIZE_64;
                buf1[0] = buf1[0] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[0];
                buf2[0] = buf2[0] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[0];
                buf1[1] = buf1[1] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[1];
                buf2[1] = buf2[1] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[1];
                buf1[2] = buf1[2] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[2];
                buf2[2] = buf2[2] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[2];
                buf1[3] = buf1[3] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[3];
                buf2[3] = buf2[3] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[3];
                buf1[4] = buf1[4] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[4];
                buf2[4] = buf2[4] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[4];
                buf1[5] = buf1[5] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[5];
                buf2[5] = buf2[5] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[5];
                buf1[6] = buf1[6] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[6];
                buf2[6] = buf2[6] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[6];
                buf1[7] = buf1[7] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[7];
                buf2[7] = buf2[7] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[7];
                quarter_pel_refinemnet_block(
                    context_ptr,
                    &context_ptr->p_best_ssd16x16[nidx],
                    src_block_index,
                    buf1,
                    buf1_stride,
                    buf2,
                    buf2_stride,
                    16,
                    16,
                    x_search_area_origin,
                    y_search_area_origin,
                    testmv,
                    &context_ptr->p_best_sad16x16[nidx],
                    &context_ptr->p_best_mv16x16[nidx],
                    it);
            }
        }
    }
    // 8x8   [64 partitions]
    for (pu_index = 0; pu_index < 64; ++pu_index) {
        nidx = tab8x8[pu_index];
        x_best_mv = _MVXT(context_ptr->p_best_full_pel_mv8x8[nidx]);
        y_best_mv = _MVYT(context_ptr->p_best_full_pel_mv8x8[nidx]);
        best_xSearchIndex = ((x_best_mv + 2) >> 2) - x_search_area_origin;
        best_ySearchIndex = ((y_best_mv + 2) >> 2) - y_search_area_origin;
        dis_x = ABS(int_xSearchIndex - best_xSearchIndex);
        dis_y = ABS(int_ySearchIndex - best_ySearchIndex);
        skip_qp_pel = 0;
        if ((dis_x) > Q_PEL_SEARCH_WIND)
            skip_qp_pel = 1;
        if ((dis_y) > Q_PEL_SEARCH_WIND)
            skip_qp_pel = 1;
        if (!skip_qp_pel) {
            for (it = 0; it < num_qp_it; it++) {
                x_mv = (int16_t)_MVXT(integer_mv) + (2 * it);
                y_mv = (int16_t)_MVYT(integer_mv) + (2 * it);
                testmv = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
                SetQuarterPelRefinementInputsOnTheFly(
                    pos_full,
                    full_stride,
                    pos_b,
                    pos_h,
                    pos_j,
                    context_ptr->interpolated_stride,
                    x_mv,
                    y_mv,
                    buf1,
                    buf1_stride,
                    buf2,
                    buf2_stride);
                block_index_shift_x = (pu_index & 0x07) << 3;
                block_index_shift_y = (pu_index >> 3) << 3;
                src_block_index =
                    block_index_shift_x + block_index_shift_y * BLOCK_SIZE_64;
                buf1[0] = buf1[0] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[0];
                buf2[0] = buf2[0] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[0];
                buf1[1] = buf1[1] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[1];
                buf2[1] = buf2[1] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[1];
                buf1[2] = buf1[2] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[2];
                buf2[2] = buf2[2] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[2];
                buf1[3] = buf1[3] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[3];
                buf2[3] = buf2[3] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[3];
                buf1[4] = buf1[4] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[4];
                buf2[4] = buf2[4] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[4];
                buf1[5] = buf1[5] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[5];
                buf2[5] = buf2[5] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[5];
                buf1[6] = buf1[6] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[6];
                buf2[6] = buf2[6] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[6];
                buf1[7] = buf1[7] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[7];
                buf2[7] = buf2[7] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[7];
                quarter_pel_refinemnet_block(context_ptr,
                                             &context_ptr->p_best_ssd8x8[nidx],
                                             src_block_index,
                                             buf1,
                                             buf1_stride,
                                             buf2,
                                             buf2_stride,
                                             8,
                                             8,
                                             x_search_area_origin,
                                             y_search_area_origin,
                                             testmv,
                                             &context_ptr->p_best_sad8x8[nidx],
                                             &context_ptr->p_best_mv8x8[nidx],
                                             it);
            }
        }
    }
    // 64x32
    for (pu_index = 0; pu_index < 2; ++pu_index) {
        block_index_shift_x = 0;
        block_index_shift_y = pu_index << 5;
        x_best_mv = _MVXT(context_ptr->p_best_full_pel_mv64x32[pu_index]);
        y_best_mv = _MVYT(context_ptr->p_best_full_pel_mv64x32[pu_index]);
        best_xSearchIndex = ((x_best_mv + 2) >> 2) - x_search_area_origin;
        best_ySearchIndex = ((y_best_mv + 2) >> 2) - y_search_area_origin;
        dis_x = ABS(int_xSearchIndex - best_xSearchIndex);
        dis_y = ABS(int_ySearchIndex - best_ySearchIndex);
        skip_qp_pel = 0;
        if ((dis_x) > Q_PEL_SEARCH_WIND)
            skip_qp_pel = 1;
        if ((dis_y) > Q_PEL_SEARCH_WIND)
            skip_qp_pel = 1;
        if (!skip_qp_pel) {
            for (it = 0; it < num_qp_it; it++) {
                x_mv = (int16_t)_MVXT(integer_mv) + (2 * it);
                y_mv = (int16_t)_MVYT(integer_mv) + (2 * it);
                testmv = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
                SetQuarterPelRefinementInputsOnTheFly(
                    pos_full,
                    full_stride,
                    pos_b,
                    pos_h,
                    pos_j,
                    context_ptr->interpolated_stride,
                    x_mv,
                    y_mv,
                    buf1,
                    buf1_stride,
                    buf2,
                    buf2_stride);
                src_block_index =
                    block_index_shift_x + block_index_shift_y * BLOCK_SIZE_64;
                buf1[0] = buf1[0] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[0];
                buf2[0] = buf2[0] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[0];
                buf1[1] = buf1[1] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[1];
                buf2[1] = buf2[1] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[1];
                buf1[2] = buf1[2] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[2];
                buf2[2] = buf2[2] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[2];
                buf1[3] = buf1[3] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[3];
                buf2[3] = buf2[3] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[3];
                buf1[4] = buf1[4] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[4];
                buf2[4] = buf2[4] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[4];
                buf1[5] = buf1[5] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[5];
                buf2[5] = buf2[5] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[5];
                buf1[6] = buf1[6] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[6];
                buf2[6] = buf2[6] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[6];
                buf1[7] = buf1[7] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[7];
                buf2[7] = buf2[7] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[7];
                quarter_pel_refinemnet_block(
                    context_ptr,
                    &context_ptr->p_best_ssd64x32[pu_index],
                    src_block_index,
                    buf1,
                    buf1_stride,
                    buf2,
                    buf2_stride,
                    64,
                    32,
                    x_search_area_origin,
                    y_search_area_origin,
                    testmv,
                    &context_ptr->p_best_sad64x32[pu_index],
                    &context_ptr->p_best_mv64x32[pu_index],
                    it);
            }
        }
    }
    // 32x16
    for (pu_index = 0; pu_index < 8; ++pu_index) {
        nidx = tab32x16[pu_index];
        block_index_shift_x = (pu_index & 0x01) << 5;
        block_index_shift_y = (pu_index >> 1) << 4;
        x_best_mv = _MVXT(context_ptr->p_best_full_pel_mv32x16[nidx]);
        y_best_mv = _MVYT(context_ptr->p_best_full_pel_mv32x16[nidx]);
        best_xSearchIndex = ((x_best_mv + 2) >> 2) - x_search_area_origin;
        best_ySearchIndex = ((y_best_mv + 2) >> 2) - y_search_area_origin;
        dis_x = ABS(int_xSearchIndex - best_xSearchIndex);
        dis_y = ABS(int_ySearchIndex - best_ySearchIndex);
        skip_qp_pel = 0;
        if ((dis_x) > Q_PEL_SEARCH_WIND)
            skip_qp_pel = 1;
        if ((dis_y) > Q_PEL_SEARCH_WIND)
            skip_qp_pel = 1;
        if (!skip_qp_pel) {
            for (it = 0; it < num_qp_it; it++) {
                x_mv = (int16_t)_MVXT(integer_mv) + (2 * it);
                y_mv = (int16_t)_MVYT(integer_mv) + (2 * it);
                testmv = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
                SetQuarterPelRefinementInputsOnTheFly(
                    pos_full,
                    full_stride,
                    pos_b,
                    pos_h,
                    pos_j,
                    context_ptr->interpolated_stride,
                    x_mv,
                    y_mv,
                    buf1,
                    buf1_stride,
                    buf2,
                    buf2_stride);
                src_block_index =
                    block_index_shift_x + block_index_shift_y * BLOCK_SIZE_64;
                buf1[0] = buf1[0] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[0];
                buf2[0] = buf2[0] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[0];
                buf1[1] = buf1[1] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[1];
                buf2[1] = buf2[1] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[1];
                buf1[2] = buf1[2] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[2];
                buf2[2] = buf2[2] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[2];
                buf1[3] = buf1[3] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[3];
                buf2[3] = buf2[3] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[3];
                buf1[4] = buf1[4] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[4];
                buf2[4] = buf2[4] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[4];
                buf1[5] = buf1[5] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[5];
                buf2[5] = buf2[5] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[5];
                buf1[6] = buf1[6] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[6];
                buf2[6] = buf2[6] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[6];
                buf1[7] = buf1[7] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[7];
                buf2[7] = buf2[7] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[7];
                quarter_pel_refinemnet_block(
                    context_ptr,
                    &context_ptr->p_best_ssd32x16[nidx],
                    src_block_index,
                    buf1,
                    buf1_stride,
                    buf2,
                    buf2_stride,
                    32,
                    16,
                    x_search_area_origin,
                    y_search_area_origin,
                    testmv,
                    &context_ptr->p_best_sad32x16[nidx],
                    &context_ptr->p_best_mv32x16[nidx],
                    it);
            }
        }
    }
    // 16x8
    for (pu_index = 0; pu_index < 32; ++pu_index) {
        nidx = tab16x8[pu_index];
        block_index_shift_x = (pu_index & 0x03) << 4;
        block_index_shift_y = (pu_index >> 2) << 3;
        x_best_mv = _MVXT(context_ptr->p_best_full_pel_mv16x8[nidx]);
        y_best_mv = _MVYT(context_ptr->p_best_full_pel_mv16x8[nidx]);
        best_xSearchIndex = ((x_best_mv + 2) >> 2) - x_search_area_origin;
        best_ySearchIndex = ((y_best_mv + 2) >> 2) - y_search_area_origin;
        dis_x = ABS(int_xSearchIndex - best_xSearchIndex);
        dis_y = ABS(int_ySearchIndex - best_ySearchIndex);
        skip_qp_pel = 0;
        if ((dis_x) > Q_PEL_SEARCH_WIND)
            skip_qp_pel = 1;
        if ((dis_y) > Q_PEL_SEARCH_WIND)
            skip_qp_pel = 1;
        if (!skip_qp_pel) {
            for (it = 0; it < num_qp_it; it++) {
                x_mv = (int16_t)_MVXT(integer_mv) + (2 * it);
                y_mv = (int16_t)_MVYT(integer_mv) + (2 * it);
                testmv = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
                SetQuarterPelRefinementInputsOnTheFly(
                    pos_full,
                    full_stride,
                    pos_b,
                    pos_h,
                    pos_j,
                    context_ptr->interpolated_stride,
                    x_mv,
                    y_mv,
                    buf1,
                    buf1_stride,
                    buf2,
                    buf2_stride);
                src_block_index =
                    block_index_shift_x + block_index_shift_y * BLOCK_SIZE_64;
                buf1[0] = buf1[0] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[0];
                buf2[0] = buf2[0] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[0];
                buf1[1] = buf1[1] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[1];
                buf2[1] = buf2[1] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[1];
                buf1[2] = buf1[2] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[2];
                buf2[2] = buf2[2] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[2];
                buf1[3] = buf1[3] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[3];
                buf2[3] = buf2[3] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[3];
                buf1[4] = buf1[4] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[4];
                buf2[4] = buf2[4] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[4];
                buf1[5] = buf1[5] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[5];
                buf2[5] = buf2[5] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[5];
                buf1[6] = buf1[6] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[6];
                buf2[6] = buf2[6] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[6];
                buf1[7] = buf1[7] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[7];
                buf2[7] = buf2[7] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[7];
                quarter_pel_refinemnet_block(context_ptr,
                                             &context_ptr->p_best_ssd16x8[nidx],
                                             src_block_index,
                                             buf1,
                                             buf1_stride,
                                             buf2,
                                             buf2_stride,
                                             16,
                                             8,
                                             x_search_area_origin,
                                             y_search_area_origin,
                                             testmv,
                                             &context_ptr->p_best_sad16x8[nidx],
                                             &context_ptr->p_best_mv16x8[nidx],
                                             it);
            }
        }
    }
    // 32x64
    for (pu_index = 0; pu_index < 2; ++pu_index) {
        block_index_shift_x = pu_index << 5;
        block_index_shift_y = 0;
        x_best_mv = _MVXT(context_ptr->p_best_full_pel_mv32x64[pu_index]);
        y_best_mv = _MVYT(context_ptr->p_best_full_pel_mv32x64[pu_index]);
        best_xSearchIndex = ((x_best_mv + 2) >> 2) - x_search_area_origin;
        best_ySearchIndex = ((y_best_mv + 2) >> 2) - y_search_area_origin;
        dis_x = ABS(int_xSearchIndex - best_xSearchIndex);
        dis_y = ABS(int_ySearchIndex - best_ySearchIndex);
        skip_qp_pel = 0;
        if ((dis_x) > Q_PEL_SEARCH_WIND)
            skip_qp_pel = 1;
        if ((dis_y) > Q_PEL_SEARCH_WIND)
            skip_qp_pel = 1;
        if (!skip_qp_pel) {
            for (it = 0; it < num_qp_it; it++) {
                x_mv = (int16_t)_MVXT(integer_mv) + (2 * it);
                y_mv = (int16_t)_MVYT(integer_mv) + (2 * it);
                testmv = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
                SetQuarterPelRefinementInputsOnTheFly(
                    pos_full,
                    full_stride,
                    pos_b,
                    pos_h,
                    pos_j,
                    context_ptr->interpolated_stride,
                    x_mv,
                    y_mv,
                    buf1,
                    buf1_stride,
                    buf2,
                    buf2_stride);
                src_block_index =
                    block_index_shift_x + block_index_shift_y * BLOCK_SIZE_64;
                buf1[0] = buf1[0] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[0];
                buf2[0] = buf2[0] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[0];
                buf1[1] = buf1[1] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[1];
                buf2[1] = buf2[1] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[1];
                buf1[2] = buf1[2] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[2];
                buf2[2] = buf2[2] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[2];
                buf1[3] = buf1[3] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[3];
                buf2[3] = buf2[3] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[3];
                buf1[4] = buf1[4] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[4];
                buf2[4] = buf2[4] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[4];
                buf1[5] = buf1[5] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[5];
                buf2[5] = buf2[5] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[5];
                buf1[6] = buf1[6] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[6];
                buf2[6] = buf2[6] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[6];
                buf1[7] = buf1[7] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[7];
                buf2[7] = buf2[7] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[7];
                quarter_pel_refinemnet_block(
                    context_ptr,
                    &context_ptr->p_best_ssd32x64[pu_index],
                    src_block_index,
                    buf1,
                    buf1_stride,
                    buf2,
                    buf2_stride,
                    32,
                    64,
                    x_search_area_origin,
                    y_search_area_origin,
                    testmv,
                    &context_ptr->p_best_sad32x64[pu_index],
                    &context_ptr->p_best_mv32x64[pu_index],
                    it);
            }
        }
    }
    // 16x32
    for (pu_index = 0; pu_index < 8; ++pu_index) {
        nidx = tab16x32[pu_index];
        block_index_shift_x = (pu_index & 0x03) << 4;
        block_index_shift_y = (pu_index >> 2) << 5;
        x_best_mv = _MVXT(context_ptr->p_best_full_pel_mv16x32[nidx]);
        y_best_mv = _MVYT(context_ptr->p_best_full_pel_mv16x32[nidx]);
        best_xSearchIndex = ((x_best_mv + 2) >> 2) - x_search_area_origin;
        best_ySearchIndex = ((y_best_mv + 2) >> 2) - y_search_area_origin;
        dis_x = ABS(int_xSearchIndex - best_xSearchIndex);
        dis_y = ABS(int_ySearchIndex - best_ySearchIndex);
        skip_qp_pel = 0;
        if ((dis_x) > Q_PEL_SEARCH_WIND)
            skip_qp_pel = 1;
        if ((dis_y) > Q_PEL_SEARCH_WIND)
            skip_qp_pel = 1;
        if (!skip_qp_pel) {
            for (it = 0; it < num_qp_it; it++) {
                x_mv = (int16_t)_MVXT(integer_mv) + (2 * it);
                y_mv = (int16_t)_MVYT(integer_mv) + (2 * it);
                testmv = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
                SetQuarterPelRefinementInputsOnTheFly(
                    pos_full,
                    full_stride,
                    pos_b,
                    pos_h,
                    pos_j,
                    context_ptr->interpolated_stride,
                    x_mv,
                    y_mv,
                    buf1,
                    buf1_stride,
                    buf2,
                    buf2_stride);
                src_block_index =
                    block_index_shift_x + block_index_shift_y * BLOCK_SIZE_64;
                buf1[0] = buf1[0] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[0];
                buf2[0] = buf2[0] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[0];
                buf1[1] = buf1[1] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[1];
                buf2[1] = buf2[1] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[1];
                buf1[2] = buf1[2] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[2];
                buf2[2] = buf2[2] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[2];
                buf1[3] = buf1[3] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[3];
                buf2[3] = buf2[3] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[3];
                buf1[4] = buf1[4] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[4];
                buf2[4] = buf2[4] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[4];
                buf1[5] = buf1[5] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[5];
                buf2[5] = buf2[5] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[5];
                buf1[6] = buf1[6] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[6];
                buf2[6] = buf2[6] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[6];
                buf1[7] = buf1[7] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[7];
                buf2[7] = buf2[7] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[7];
                quarter_pel_refinemnet_block(
                    context_ptr,
                    &context_ptr->p_best_ssd16x32[nidx],
                    src_block_index,
                    buf1,
                    buf1_stride,
                    buf2,
                    buf2_stride,
                    16,
                    32,
                    x_search_area_origin,
                    y_search_area_origin,
                    testmv,
                    &context_ptr->p_best_sad16x32[nidx],
                    &context_ptr->p_best_mv16x32[nidx],
                    it);
            }
        }
    }
    // 8x16
    for (pu_index = 0; pu_index < 32; ++pu_index) {
        nidx = tab8x16[pu_index];
        block_index_shift_x = (pu_index & 0x07) << 3;
        block_index_shift_y = (pu_index >> 3) << 4;
        x_best_mv = _MVXT(context_ptr->p_best_full_pel_mv8x16[nidx]);
        y_best_mv = _MVYT(context_ptr->p_best_full_pel_mv8x16[nidx]);
        best_xSearchIndex = ((x_best_mv + 2) >> 2) - x_search_area_origin;
        best_ySearchIndex = ((y_best_mv + 2) >> 2) - y_search_area_origin;
        dis_x = ABS(int_xSearchIndex - best_xSearchIndex);
        dis_y = ABS(int_ySearchIndex - best_ySearchIndex);
        skip_qp_pel = 0;
        if ((dis_x) > Q_PEL_SEARCH_WIND)
            skip_qp_pel = 1;
        if ((dis_y) > Q_PEL_SEARCH_WIND)
            skip_qp_pel = 1;
        if (!skip_qp_pel) {
            for (it = 0; it < num_qp_it; it++) {
                x_mv = (int16_t)_MVXT(integer_mv) + (2 * it);
                y_mv = (int16_t)_MVYT(integer_mv) + (2 * it);
                testmv = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
                SetQuarterPelRefinementInputsOnTheFly(
                    pos_full,
                    full_stride,
                    pos_b,
                    pos_h,
                    pos_j,
                    context_ptr->interpolated_stride,
                    x_mv,
                    y_mv,
                    buf1,
                    buf1_stride,
                    buf2,
                    buf2_stride);
                src_block_index =
                    block_index_shift_x + block_index_shift_y * BLOCK_SIZE_64;
                buf1[0] = buf1[0] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[0];
                buf2[0] = buf2[0] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[0];
                buf1[1] = buf1[1] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[1];
                buf2[1] = buf2[1] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[1];
                buf1[2] = buf1[2] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[2];
                buf2[2] = buf2[2] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[2];
                buf1[3] = buf1[3] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[3];
                buf2[3] = buf2[3] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[3];
                buf1[4] = buf1[4] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[4];
                buf2[4] = buf2[4] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[4];
                buf1[5] = buf1[5] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[5];
                buf2[5] = buf2[5] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[5];
                buf1[6] = buf1[6] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[6];
                buf2[6] = buf2[6] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[6];
                buf1[7] = buf1[7] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[7];
                buf2[7] = buf2[7] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[7];
                quarter_pel_refinemnet_block(context_ptr,
                                             &context_ptr->p_best_ssd8x16[nidx],
                                             src_block_index,
                                             buf1,
                                             buf1_stride,
                                             buf2,
                                             buf2_stride,
                                             8,
                                             16,
                                             x_search_area_origin,
                                             y_search_area_origin,
                                             testmv,
                                             &context_ptr->p_best_sad8x16[nidx],
                                             &context_ptr->p_best_mv8x16[nidx],
                                             it);
            }
        }
    }
    // 32x8
    for (pu_index = 0; pu_index < 16; ++pu_index) {
        nidx = tab32x8[pu_index];
        block_index_shift_x = (pu_index & 0x01) << 5;
        block_index_shift_y = (pu_index >> 1) << 3;
        x_best_mv = _MVXT(context_ptr->p_best_full_pel_mv32x8[nidx]);
        y_best_mv = _MVYT(context_ptr->p_best_full_pel_mv32x8[nidx]);
        best_xSearchIndex = ((x_best_mv + 2) >> 2) - x_search_area_origin;
        best_ySearchIndex = ((y_best_mv + 2) >> 2) - y_search_area_origin;
        dis_x = ABS(int_xSearchIndex - best_xSearchIndex);
        dis_y = ABS(int_ySearchIndex - best_ySearchIndex);
        skip_qp_pel = 0;
        if ((dis_x) > Q_PEL_SEARCH_WIND)
            skip_qp_pel = 1;
        if ((dis_y) > Q_PEL_SEARCH_WIND)
            skip_qp_pel = 1;
        if (!skip_qp_pel) {
            for (it = 0; it < num_qp_it; it++) {
                x_mv = (int16_t)_MVXT(integer_mv) + (2 * it);
                y_mv = (int16_t)_MVYT(integer_mv) + (2 * it);
                testmv = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
                SetQuarterPelRefinementInputsOnTheFly(
                    pos_full,
                    full_stride,
                    pos_b,
                    pos_h,
                    pos_j,
                    context_ptr->interpolated_stride,
                    x_mv,
                    y_mv,
                    buf1,
                    buf1_stride,
                    buf2,
                    buf2_stride);
                src_block_index =
                    block_index_shift_x + block_index_shift_y * BLOCK_SIZE_64;
                buf1[0] = buf1[0] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[0];
                buf2[0] = buf2[0] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[0];
                buf1[1] = buf1[1] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[1];
                buf2[1] = buf2[1] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[1];
                buf1[2] = buf1[2] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[2];
                buf2[2] = buf2[2] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[2];
                buf1[3] = buf1[3] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[3];
                buf2[3] = buf2[3] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[3];
                buf1[4] = buf1[4] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[4];
                buf2[4] = buf2[4] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[4];
                buf1[5] = buf1[5] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[5];
                buf2[5] = buf2[5] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[5];
                buf1[6] = buf1[6] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[6];
                buf2[6] = buf2[6] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[6];
                buf1[7] = buf1[7] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[7];
                buf2[7] = buf2[7] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[7];
                quarter_pel_refinemnet_block(context_ptr,
                                             &context_ptr->p_best_ssd32x8[nidx],
                                             src_block_index,
                                             buf1,
                                             buf1_stride,
                                             buf2,
                                             buf2_stride,
                                             32,
                                             8,
                                             x_search_area_origin,
                                             y_search_area_origin,
                                             testmv,
                                             &context_ptr->p_best_sad32x8[nidx],
                                             &context_ptr->p_best_mv32x8[nidx],
                                             it);
            }
        }
    }

    // 8x32
    for (pu_index = 0; pu_index < 16; ++pu_index) {
        nidx = tab8x32[pu_index];
        block_index_shift_x = (pu_index & 0x07) << 3;
        block_index_shift_y = (pu_index >> 3) << 5;
        x_best_mv = _MVXT(context_ptr->p_best_full_pel_mv8x32[nidx]);
        y_best_mv = _MVYT(context_ptr->p_best_full_pel_mv8x32[nidx]);
        best_xSearchIndex = ((x_best_mv + 2) >> 2) - x_search_area_origin;
        best_ySearchIndex = ((y_best_mv + 2) >> 2) - y_search_area_origin;
        dis_x = ABS(int_xSearchIndex - best_xSearchIndex);
        dis_y = ABS(int_ySearchIndex - best_ySearchIndex);
        skip_qp_pel = 0;
        if ((dis_x) > Q_PEL_SEARCH_WIND)
            skip_qp_pel = 1;
        if ((dis_y) > Q_PEL_SEARCH_WIND)
            skip_qp_pel = 1;
        if (!skip_qp_pel) {
            for (it = 0; it < num_qp_it; it++) {
                x_mv = (int16_t)_MVXT(integer_mv) + (2 * it);
                y_mv = (int16_t)_MVYT(integer_mv) + (2 * it);
                testmv = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
                SetQuarterPelRefinementInputsOnTheFly(
                    pos_full,
                    full_stride,
                    pos_b,
                    pos_h,
                    pos_j,
                    context_ptr->interpolated_stride,
                    x_mv,
                    y_mv,
                    buf1,
                    buf1_stride,
                    buf2,
                    buf2_stride);
                src_block_index =
                    block_index_shift_x + block_index_shift_y * BLOCK_SIZE_64;
                buf1[0] = buf1[0] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[0];
                buf2[0] = buf2[0] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[0];
                buf1[1] = buf1[1] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[1];
                buf2[1] = buf2[1] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[1];
                buf1[2] = buf1[2] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[2];
                buf2[2] = buf2[2] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[2];
                buf1[3] = buf1[3] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[3];
                buf2[3] = buf2[3] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[3];
                buf1[4] = buf1[4] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[4];
                buf2[4] = buf2[4] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[4];
                buf1[5] = buf1[5] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[5];
                buf2[5] = buf2[5] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[5];
                buf1[6] = buf1[6] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[6];
                buf2[6] = buf2[6] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[6];
                buf1[7] = buf1[7] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[7];
                buf2[7] = buf2[7] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[7];
                quarter_pel_refinemnet_block(context_ptr,
                                             &context_ptr->p_best_ssd8x32[nidx],
                                             src_block_index,
                                             buf1,
                                             buf1_stride,
                                             buf2,
                                             buf2_stride,
                                             8,
                                             32,
                                             x_search_area_origin,
                                             y_search_area_origin,
                                             testmv,
                                             &context_ptr->p_best_sad8x32[nidx],
                                             &context_ptr->p_best_mv8x32[nidx],
                                             it);
            }
        }
    }

    // 64x16
    for (pu_index = 0; pu_index < 4; ++pu_index) {
        nidx = pu_index;
        block_index_shift_x = 0;
        block_index_shift_y = pu_index << 4;
        x_best_mv = _MVXT(context_ptr->p_best_full_pel_mv64x16[nidx]);
        y_best_mv = _MVYT(context_ptr->p_best_full_pel_mv64x16[nidx]);
        best_xSearchIndex = ((x_best_mv + 2) >> 2) - x_search_area_origin;
        best_ySearchIndex = ((y_best_mv + 2) >> 2) - y_search_area_origin;
        dis_x = ABS(int_xSearchIndex - best_xSearchIndex);
        dis_y = ABS(int_ySearchIndex - best_ySearchIndex);
        skip_qp_pel = 0;
        if ((dis_x) > Q_PEL_SEARCH_WIND)
            skip_qp_pel = 1;
        if ((dis_y) > Q_PEL_SEARCH_WIND)
            skip_qp_pel = 1;
        if (!skip_qp_pel) {
            for (it = 0; it < num_qp_it; it++) {
                x_mv = (int16_t)_MVXT(integer_mv) + (2 * it);
                y_mv = (int16_t)_MVYT(integer_mv) + (2 * it);
                testmv = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
                SetQuarterPelRefinementInputsOnTheFly(
                    pos_full,
                    full_stride,
                    pos_b,
                    pos_h,
                    pos_j,
                    context_ptr->interpolated_stride,
                    x_mv,
                    y_mv,
                    buf1,
                    buf1_stride,
                    buf2,
                    buf2_stride);
                src_block_index =
                    block_index_shift_x + block_index_shift_y * BLOCK_SIZE_64;
                buf1[0] = buf1[0] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[0];
                buf2[0] = buf2[0] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[0];
                buf1[1] = buf1[1] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[1];
                buf2[1] = buf2[1] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[1];
                buf1[2] = buf1[2] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[2];
                buf2[2] = buf2[2] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[2];
                buf1[3] = buf1[3] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[3];
                buf2[3] = buf2[3] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[3];
                buf1[4] = buf1[4] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[4];
                buf2[4] = buf2[4] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[4];
                buf1[5] = buf1[5] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[5];
                buf2[5] = buf2[5] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[5];
                buf1[6] = buf1[6] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[6];
                buf2[6] = buf2[6] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[6];
                buf1[7] = buf1[7] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[7];
                buf2[7] = buf2[7] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[7];
                quarter_pel_refinemnet_block(
                    context_ptr,
                    &context_ptr->p_best_ssd64x16[nidx],
                    src_block_index,
                    buf1,
                    buf1_stride,
                    buf2,
                    buf2_stride,
                    64,
                    16,
                    x_search_area_origin,
                    y_search_area_origin,
                    testmv,
                    &context_ptr->p_best_sad64x16[nidx],
                    &context_ptr->p_best_mv64x16[nidx],
                    it);
            }
        }
    }
    // 16x64
    for (pu_index = 0; pu_index < 4; ++pu_index) {
        nidx = pu_index;
        block_index_shift_x = pu_index << 4;
        block_index_shift_y = 0;
        x_best_mv = _MVXT(context_ptr->p_best_full_pel_mv16x64[nidx]);
        y_best_mv = _MVYT(context_ptr->p_best_full_pel_mv16x64[nidx]);
        best_xSearchIndex = ((x_best_mv + 2) >> 2) - x_search_area_origin;
        best_ySearchIndex = ((y_best_mv + 2) >> 2) - y_search_area_origin;
        dis_x = ABS(int_xSearchIndex - best_xSearchIndex);
        dis_y = ABS(int_ySearchIndex - best_ySearchIndex);
        skip_qp_pel = 0;
        if ((dis_x) > Q_PEL_SEARCH_WIND)
            skip_qp_pel = 1;
        if ((dis_y) > Q_PEL_SEARCH_WIND)
            skip_qp_pel = 1;
        if (!skip_qp_pel) {
            for (it = 0; it < num_qp_it; it++) {
                x_mv = (int16_t)_MVXT(integer_mv) + (2 * it);
                y_mv = (int16_t)_MVYT(integer_mv) + (2 * it);
                testmv = ((uint16_t)y_mv << 16) | ((uint16_t)x_mv);
                SetQuarterPelRefinementInputsOnTheFly(
                    pos_full,
                    full_stride,
                    pos_b,
                    pos_h,
                    pos_j,
                    context_ptr->interpolated_stride,
                    x_mv,
                    y_mv,
                    buf1,
                    buf1_stride,
                    buf2,
                    buf2_stride);
                src_block_index =
                    block_index_shift_x + block_index_shift_y * BLOCK_SIZE_64;
                buf1[0] = buf1[0] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[0];
                buf2[0] = buf2[0] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[0];
                buf1[1] = buf1[1] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[1];
                buf2[1] = buf2[1] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[1];
                buf1[2] = buf1[2] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[2];
                buf2[2] = buf2[2] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[2];
                buf1[3] = buf1[3] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[3];
                buf2[3] = buf2[3] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[3];
                buf1[4] = buf1[4] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[4];
                buf2[4] = buf2[4] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[4];
                buf1[5] = buf1[5] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[5];
                buf2[5] = buf2[5] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[5];
                buf1[6] = buf1[6] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[6];
                buf2[6] = buf2[6] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[6];
                buf1[7] = buf1[7] + block_index_shift_x +
                          block_index_shift_y * buf1_stride[7];
                buf2[7] = buf2[7] + block_index_shift_x +
                          block_index_shift_y * buf2_stride[7];
                quarter_pel_refinemnet_block(
                    context_ptr,
                    &context_ptr->p_best_ssd16x64[nidx],
                    src_block_index,
                    buf1,
                    buf1_stride,
                    buf2,
                    buf2_stride,
                    16,
                    64,
                    x_search_area_origin,
                    y_search_area_origin,
                    testmv,
                    &context_ptr->p_best_sad16x64[nidx],
                    &context_ptr->p_best_mv16x64[nidx],
                    it);
            }
        }
    }
    return;
}
void HmeOneQuadrantLevel0(
    PictureParentControlSet *picture_control_set_ptr,
    MeContext *context_ptr,  // input/output parameter, ME context Ptr, used to
                             // get/update ME results
    int16_t origin_x,        // input parameter, SB position in the horizontal
                             // direction- sixteenth resolution
    int16_t origin_y,        // input parameter, SB position in the vertical
                             // direction- sixteenth resolution
    uint32_t sb_width,   // input parameter, SB pwidth - sixteenth resolution
    uint32_t sb_height,  // input parameter, SB height - sixteenth resolution
    int16_t xHmeSearchCenter,  // input parameter, HME search center in the
                               // horizontal direction
    int16_t yHmeSearchCenter,  // input parameter, HME search center in the
                               // vertical direction
    EbPictureBufferDesc *
        sixteenthRefPicPtr,  // input parameter, sixteenth reference Picture Ptr
    uint64_t *level0BestSad,       // output parameter, Level0 SAD at
                                   // (searchRegionNumberInWidth,
                                   // searchRegionNumberInHeight)
    int16_t *xLevel0SearchCenter,  // output parameter, Level0 xMV at
                                   // (searchRegionNumberInWidth,
                                   // searchRegionNumberInHeight)
    int16_t *yLevel0SearchCenter,  // output parameter, Level0 yMV at
                                   // (searchRegionNumberInWidth,
                                   // searchRegionNumberInHeight)
    uint32_t searchAreaMultiplierX, uint32_t searchAreaMultiplierY,
    EbAsm asm_type) {
    int16_t xTopLeftSearchRegion;
    int16_t yTopLeftSearchRegion;
    uint32_t searchRegionIndex;
    int16_t x_search_area_origin;
    int16_t y_search_area_origin;
    int16_t xSearchRegionDistance;
    int16_t ySearchRegionDistance;

    int16_t padWidth;
    int16_t padHeight;

    (void)picture_control_set_ptr;
    // Round up x_HME_L0 to be a multiple of 16
    int16_t search_area_width =
        (int16_t)((((((context_ptr->hme_level0_total_search_area_width *
                       searchAreaMultiplierX) /
                      100))) +
                   15) &
                  ~0x0F);
    int16_t search_area_height =
        (int16_t)(((context_ptr->hme_level0_total_search_area_height *
                    searchAreaMultiplierY) /
                   100));
    xSearchRegionDistance = xHmeSearchCenter;
    ySearchRegionDistance = yHmeSearchCenter;
    padWidth = (int16_t)(sixteenthRefPicPtr->origin_x) - 1;
    padHeight = (int16_t)(sixteenthRefPicPtr->origin_y) - 1;

    x_search_area_origin =
        -(int16_t)(search_area_width >> 1) + xSearchRegionDistance;
    y_search_area_origin =
        -(int16_t)(search_area_height >> 1) + ySearchRegionDistance;

    // Correct the left edge of the Search Area if it is not on the reference
    // Picture
    x_search_area_origin = ((origin_x + x_search_area_origin) < -padWidth)
                               ? -padWidth - origin_x
                               : x_search_area_origin;

    search_area_width =
        ((origin_x + x_search_area_origin) < -padWidth)
            ? search_area_width -
                  (-padWidth - (origin_x + x_search_area_origin))
            : search_area_width;

    // Correct the right edge of the Search Area if its not on the reference
    // Picture
    x_search_area_origin =
        ((origin_x + x_search_area_origin) >
         (int16_t)sixteenthRefPicPtr->width - 1)
            ? x_search_area_origin - ((origin_x + x_search_area_origin) -
                                      ((int16_t)sixteenthRefPicPtr->width - 1))
            : x_search_area_origin;

    search_area_width =
        ((origin_x + x_search_area_origin + search_area_width) >
         (int16_t)sixteenthRefPicPtr->width)
            ? MAX(1,
                  search_area_width -
                      ((origin_x + x_search_area_origin + search_area_width) -
                       (int16_t)sixteenthRefPicPtr->width))
            : search_area_width;

    // Round down x_HME to be a multiple of 16 as cropping already performed
    search_area_width = (search_area_width < 16) ? search_area_width
                                                 : search_area_width & ~0x0F;
    // Correct the top edge of the Search Area if it is not on the reference
    // Picture
    y_search_area_origin = ((origin_y + y_search_area_origin) < -padHeight)
                               ? -padHeight - origin_y
                               : y_search_area_origin;

    search_area_height =
        ((origin_y + y_search_area_origin) < -padHeight)
            ? search_area_height -
                  (-padHeight - (origin_y + y_search_area_origin))
            : search_area_height;

    // Correct the bottom edge of the Search Area if its not on the reference
    // Picture
    y_search_area_origin =
        ((origin_y + y_search_area_origin) >
         (int16_t)sixteenthRefPicPtr->height - 1)
            ? y_search_area_origin - ((origin_y + y_search_area_origin) -
                                      ((int16_t)sixteenthRefPicPtr->height - 1))
            : y_search_area_origin;

    search_area_height =
        (origin_y + y_search_area_origin + search_area_height >
         (int16_t)sixteenthRefPicPtr->height)
            ? MAX(1,
                  search_area_height -
                      ((origin_y + y_search_area_origin + search_area_height) -
                       (int16_t)sixteenthRefPicPtr->height))
            : search_area_height;

    xTopLeftSearchRegion = ((int16_t)sixteenthRefPicPtr->origin_x + origin_x) +
                           x_search_area_origin;
    yTopLeftSearchRegion = ((int16_t)sixteenthRefPicPtr->origin_y + origin_y) +
                           y_search_area_origin;
    searchRegionIndex = xTopLeftSearchRegion +
                        yTopLeftSearchRegion * sixteenthRefPicPtr->stride_y;

    if (context_ptr->hme_search_type == HME_SPARSE) {
        sad_loop_kernel_sparse(
            &context_ptr->sixteenth_sb_buffer[0],
            context_ptr->sixteenth_sb_buffer_stride,
            &sixteenthRefPicPtr->buffer_y[searchRegionIndex],
            (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                ? sixteenthRefPicPtr->stride_y
                : sixteenthRefPicPtr->stride_y * 2,
            (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                ? sb_height
                : sb_height >> 1,
            sb_width,
            /* results */
            level0BestSad,
            xLevel0SearchCenter,
            yLevel0SearchCenter,
            /* range */
            sixteenthRefPicPtr->stride_y,
            search_area_width,
            search_area_height);
    } else {
        if (((search_area_width & 15) == 0) && (asm_type == ASM_AVX2)) {
            sad_loop_kernel_avx2_hme_l0_intrin(
                &context_ptr->sixteenth_sb_buffer[0],
                context_ptr->sixteenth_sb_buffer_stride,
                &sixteenthRefPicPtr->buffer_y[searchRegionIndex],
                (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                    ? sixteenthRefPicPtr->stride_y
                    : sixteenthRefPicPtr->stride_y * 2,
                (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                    ? sb_height
                    : sb_height >> 1,
                sb_width,
                // results
                level0BestSad,
                xLevel0SearchCenter,
                yLevel0SearchCenter,
                // range
                sixteenthRefPicPtr->stride_y,
                search_area_width,
                search_area_height);
        } else if ((search_area_width & 15) == 0) {
            // Only width equals 16 (LCU equals 64) is updated
            // other width sizes work with the old code as the one
            // in"sad_loop_kernel_sse4_1_intrin"
            sad_loop_kernel_sse4_1_hme_l0_intrin(
                &context_ptr->sixteenth_sb_buffer[0],
                context_ptr->sixteenth_sb_buffer_stride,
                &sixteenthRefPicPtr->buffer_y[searchRegionIndex],
                (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                    ? sixteenthRefPicPtr->stride_y
                    : sixteenthRefPicPtr->stride_y * 2,
                (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                    ? sb_height
                    : sb_height >> 1,
                sb_width,
                /* results */
                level0BestSad,
                xLevel0SearchCenter,
                yLevel0SearchCenter,
                /* range */
                sixteenthRefPicPtr->stride_y,
                search_area_width,
                search_area_height);
        } else {
            // Put the first search location into level0 results
            sad_loop_kernel(
                &context_ptr->sixteenth_sb_buffer[0],
                context_ptr->sixteenth_sb_buffer_stride,
                &sixteenthRefPicPtr->buffer_y[searchRegionIndex],
                (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                    ? sixteenthRefPicPtr->stride_y
                    : sixteenthRefPicPtr->stride_y * 2,
                (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                    ? sb_height
                    : sb_height >> 1,
                sb_width,
                /* results */
                level0BestSad,
                xLevel0SearchCenter,
                yLevel0SearchCenter,
                /* range */
                sixteenthRefPicPtr->stride_y,
                search_area_width,
                search_area_height);
        }
    }

    *level0BestSad =
        (context_ptr->hme_search_method == FULL_SAD_SEARCH)
            ? *level0BestSad
            : *level0BestSad *
                  2;  // Multiply by 2 because considered only ever other line
    *xLevel0SearchCenter += x_search_area_origin;
    *xLevel0SearchCenter *=
        4;  // Multiply by 4 because operating on 1/4 resolution
    *yLevel0SearchCenter += y_search_area_origin;
    *yLevel0SearchCenter *=
        4;  // Multiply by 4 because operating on 1/4 resolution

    return;
}

void HmeLevel0(
    PictureParentControlSet *picture_control_set_ptr,
    MeContext *context_ptr,  // input/output parameter, ME context Ptr, used to
                             // get/update ME results
    int16_t origin_x,        // input parameter, SB position in the horizontal
                             // direction- sixteenth resolution
    int16_t origin_y,        // input parameter, SB position in the vertical
                             // direction- sixteenth resolution
    uint32_t sb_width,   // input parameter, SB pwidth - sixteenth resolution
    uint32_t sb_height,  // input parameter, SB height - sixteenth resolution
    int16_t xHmeSearchCenter,  // input parameter, HME search center in the
                               // horizontal direction
    int16_t yHmeSearchCenter,  // input parameter, HME search center in the
                               // vertical direction
    EbPictureBufferDesc *
        sixteenthRefPicPtr,  // input parameter, sixteenth reference Picture Ptr
    uint32_t searchRegionNumberInWidth,   // input parameter, search region
                                          // number in the horizontal direction
    uint32_t searchRegionNumberInHeight,  // input parameter, search region
                                          // number in the vertical direction
    uint64_t *level0BestSad,              // output parameter, Level0 SAD at
                                          // (searchRegionNumberInWidth,
                                          // searchRegionNumberInHeight)
    int16_t *xLevel0SearchCenter,         // output parameter, Level0 xMV at
                                          // (searchRegionNumberInWidth,
                                          // searchRegionNumberInHeight)
    int16_t *yLevel0SearchCenter,         // output parameter, Level0 yMV at
                                          // (searchRegionNumberInWidth,
                                          // searchRegionNumberInHeight)
    uint32_t searchAreaMultiplierX, uint32_t searchAreaMultiplierY,
    EbAsm asm_type) {
    int16_t xTopLeftSearchRegion;
    int16_t yTopLeftSearchRegion;
    uint32_t searchRegionIndex;
    int16_t x_search_area_origin;
    int16_t y_search_area_origin;
    int16_t xSearchRegionDistance;
    int16_t ySearchRegionDistance;

    int16_t padWidth;
    int16_t padHeight;

    // Adjust SR size based on the searchAreaShift

    (void)picture_control_set_ptr;
    // Round up x_HME_L0 to be a multiple of 16
    int16_t search_area_width =
        (int16_t)((((((context_ptr->hme_level0_search_area_in_width_array
                           [searchRegionNumberInWidth] *
                       searchAreaMultiplierX) /
                      100))) +
                   15) &
                  ~0x0F);
    int16_t search_area_height =
        (int16_t)(((context_ptr->hme_level0_search_area_in_height_array
                        [searchRegionNumberInHeight] *
                    searchAreaMultiplierY) /
                   100));

    xSearchRegionDistance = xHmeSearchCenter;
    ySearchRegionDistance = yHmeSearchCenter;
    padWidth = (int16_t)(sixteenthRefPicPtr->origin_x) - 1;
    padHeight = (int16_t)(sixteenthRefPicPtr->origin_y) - 1;

    while (searchRegionNumberInWidth) {
        searchRegionNumberInWidth--;
        xSearchRegionDistance +=
            (int16_t)(((context_ptr->hme_level0_search_area_in_width_array
                            [searchRegionNumberInWidth] *
                        searchAreaMultiplierX) /
                       100));
    }

    while (searchRegionNumberInHeight) {
        searchRegionNumberInHeight--;
        ySearchRegionDistance +=
            (int16_t)(((context_ptr->hme_level0_search_area_in_height_array
                            [searchRegionNumberInHeight] *
                        searchAreaMultiplierY) /
                       100));
    }
    x_search_area_origin =
        -(int16_t)((((context_ptr->hme_level0_total_search_area_width *
                      searchAreaMultiplierX) /
                     100)) >>
                   1) +
        xSearchRegionDistance;
    y_search_area_origin =
        -(int16_t)((((context_ptr->hme_level0_total_search_area_height *
                      searchAreaMultiplierY) /
                     100)) >>
                   1) +
        ySearchRegionDistance;

    // Correct the left edge of the Search Area if it is not on the reference
    // Picture
    x_search_area_origin = ((origin_x + x_search_area_origin) < -padWidth)
                               ? -padWidth - origin_x
                               : x_search_area_origin;

    search_area_width =
        ((origin_x + x_search_area_origin) < -padWidth)
            ? search_area_width -
                  (-padWidth - (origin_x + x_search_area_origin))
            : search_area_width;

    // Correct the right edge of the Search Area if its not on the reference
    // Picture
    x_search_area_origin =
        ((origin_x + x_search_area_origin) >
         (int16_t)sixteenthRefPicPtr->width - 1)
            ? x_search_area_origin - ((origin_x + x_search_area_origin) -
                                      ((int16_t)sixteenthRefPicPtr->width - 1))
            : x_search_area_origin;

    search_area_width =
        ((origin_x + x_search_area_origin + search_area_width) >
         (int16_t)sixteenthRefPicPtr->width)
            ? MAX(1,
                  search_area_width -
                      ((origin_x + x_search_area_origin + search_area_width) -
                       (int16_t)sixteenthRefPicPtr->width))
            : search_area_width;

    // Round down x_HME to be a multiple of 16 as cropping already performed
    search_area_width = (search_area_width < 16) ? search_area_width
                                                 : search_area_width & ~0x0F;

    // Correct the top edge of the Search Area if it is not on the reference
    // Picture
    y_search_area_origin = ((origin_y + y_search_area_origin) < -padHeight)
                               ? -padHeight - origin_y
                               : y_search_area_origin;

    search_area_height =
        ((origin_y + y_search_area_origin) < -padHeight)
            ? search_area_height -
                  (-padHeight - (origin_y + y_search_area_origin))
            : search_area_height;

    // Correct the bottom edge of the Search Area if its not on the reference
    // Picture
    y_search_area_origin =
        ((origin_y + y_search_area_origin) >
         (int16_t)sixteenthRefPicPtr->height - 1)
            ? y_search_area_origin - ((origin_y + y_search_area_origin) -
                                      ((int16_t)sixteenthRefPicPtr->height - 1))
            : y_search_area_origin;

    search_area_height =
        (origin_y + y_search_area_origin + search_area_height >
         (int16_t)sixteenthRefPicPtr->height)
            ? MAX(1,
                  search_area_height -
                      ((origin_y + y_search_area_origin + search_area_height) -
                       (int16_t)sixteenthRefPicPtr->height))
            : search_area_height;

    xTopLeftSearchRegion = ((int16_t)sixteenthRefPicPtr->origin_x + origin_x) +
                           x_search_area_origin;
    yTopLeftSearchRegion = ((int16_t)sixteenthRefPicPtr->origin_y + origin_y) +
                           y_search_area_origin;
    searchRegionIndex = xTopLeftSearchRegion +
                        yTopLeftSearchRegion * sixteenthRefPicPtr->stride_y;

    if (((sb_width & 7) == 0) || (sb_width == 4)) {
        if (((search_area_width & 15) == 0) && (asm_type == ASM_AVX2)) {
            sad_loop_kernel_avx2_hme_l0_intrin(
                &context_ptr->sixteenth_sb_buffer[0],
                context_ptr->sixteenth_sb_buffer_stride,
                &sixteenthRefPicPtr->buffer_y[searchRegionIndex],
                (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                    ? sixteenthRefPicPtr->stride_y
                    : sixteenthRefPicPtr->stride_y * 2,
                (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                    ? sb_height
                    : sb_height >> 1,
                sb_width,
                // results
                level0BestSad,
                xLevel0SearchCenter,
                yLevel0SearchCenter,
                // range
                sixteenthRefPicPtr->stride_y,
                search_area_width,
                search_area_height);
        } else if ((search_area_width & 15) == 0) {
            // Only width equals 16 (LCU equals 64) is updated
            // other width sizes work with the old code as the one
            // in"sad_loop_kernel_sse4_1_intrin"
            sad_loop_kernel_sse4_1_hme_l0_intrin(
                &context_ptr->sixteenth_sb_buffer[0],
                context_ptr->sixteenth_sb_buffer_stride,
                &sixteenthRefPicPtr->buffer_y[searchRegionIndex],
                (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                    ? sixteenthRefPicPtr->stride_y
                    : sixteenthRefPicPtr->stride_y * 2,
                (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                    ? sb_height
                    : sb_height >> 1,
                sb_width,
                /* results */
                level0BestSad,
                xLevel0SearchCenter,
                yLevel0SearchCenter,
                /* range */
                sixteenthRefPicPtr->stride_y,
                search_area_width,
                search_area_height);
        } else {
            // Put the first search location into level0 results
            sad_loop_kernel(
                &context_ptr->sixteenth_sb_buffer[0],
                context_ptr->sixteenth_sb_buffer_stride,
                &sixteenthRefPicPtr->buffer_y[searchRegionIndex],
                (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                    ? sixteenthRefPicPtr->stride_y
                    : sixteenthRefPicPtr->stride_y * 2,
                (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                    ? sb_height
                    : sb_height >> 1,
                sb_width,
                /* results */
                level0BestSad,
                xLevel0SearchCenter,
                yLevel0SearchCenter,
                /* range */
                sixteenthRefPicPtr->stride_y,
                search_area_width,
                search_area_height);
        }
    } else {
        sad_loop_kernel_c(&context_ptr->sixteenth_sb_buffer[0],
                        context_ptr->sixteenth_sb_buffer_stride,
                        &sixteenthRefPicPtr->buffer_y[searchRegionIndex],
                        (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                            ? sixteenthRefPicPtr->stride_y
                            : sixteenthRefPicPtr->stride_y * 2,
                        (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                            ? sb_height
                            : sb_height >> 1,
                        sb_width,
                        /* results */
                        level0BestSad,
                        xLevel0SearchCenter,
                        yLevel0SearchCenter,
                        /* range */
                        sixteenthRefPicPtr->stride_y,
                        search_area_width,
                        search_area_height);
    }

    *level0BestSad =
        (context_ptr->hme_search_method == FULL_SAD_SEARCH)
            ? *level0BestSad
            : *level0BestSad *
                  2;  // Multiply by 2 because considered only ever other line
    *xLevel0SearchCenter += x_search_area_origin;
    *xLevel0SearchCenter *=
        4;  // Multiply by 4 because operating on 1/4 resolution
    *yLevel0SearchCenter += y_search_area_origin;
    *yLevel0SearchCenter *=
        4;  // Multiply by 4 because operating on 1/4 resolution

    return;
}

void HmeLevel1(
    MeContext *context_ptr,  // input/output parameter, ME context Ptr, used to
                             // get/update ME results
    int16_t origin_x,        // input parameter, SB position in the horizontal
                             // direction - quarter resolution
    int16_t origin_y,  // input parameter, SB position in the vertical direction
                       // - quarter resolution
    uint32_t sb_width,   // input parameter, SB pwidth - quarter resolution
    uint32_t sb_height,  // input parameter, SB height - quarter resolution
    EbPictureBufferDesc
        *quarterRefPicPtr,  // input parameter, quarter reference Picture Ptr
    int16_t hmeLevel1SearchAreaInWidth,   // input parameter, hme level 1 search
                                          // area in width
    int16_t hmeLevel1SearchAreaInHeight,  // input parameter, hme level 1 search
                                          // area in height
    int16_t xLevel0SearchCenter,          // input parameter, best Level0 xMV at
                                          // (searchRegionNumberInWidth,
                                          // searchRegionNumberInHeight)
    int16_t yLevel0SearchCenter,          // input parameter, best Level0 yMV at
                                          // (searchRegionNumberInWidth,
                                          // searchRegionNumberInHeight)
    uint64_t *level1BestSad,              // output parameter, Level1 SAD at
                                          // (searchRegionNumberInWidth,
                                          // searchRegionNumberInHeight)
    int16_t *xLevel1SearchCenter,         // output parameter, Level1 xMV at
                                          // (searchRegionNumberInWidth,
                                          // searchRegionNumberInHeight)
    int16_t *yLevel1SearchCenter         // output parameter, Level1 yMV at
                                          // (searchRegionNumberInWidth,
                                          // searchRegionNumberInHeight)
    ) {
    int16_t xTopLeftSearchRegion;
    int16_t yTopLeftSearchRegion;
    uint32_t searchRegionIndex;
    // Round up x_HME_L0 to be a multiple of 8
    int16_t search_area_width =
        (int16_t)((hmeLevel1SearchAreaInWidth + 7) & ~0x07);
    int16_t search_area_height = hmeLevel1SearchAreaInHeight;

    int16_t x_search_area_origin;
    int16_t y_search_area_origin;

    int16_t padWidth = (int16_t)(quarterRefPicPtr->origin_x) - 1;
    int16_t padHeight = (int16_t)(quarterRefPicPtr->origin_y) - 1;

    x_search_area_origin = -(search_area_width >> 1) + xLevel0SearchCenter;
    y_search_area_origin = -(search_area_height >> 1) + yLevel0SearchCenter;

    // Correct the left edge of the Search Area if it is not on the reference
    // Picture
    x_search_area_origin = ((origin_x + x_search_area_origin) < -padWidth)
                               ? -padWidth - origin_x
                               : x_search_area_origin;

    search_area_width =
        ((origin_x + x_search_area_origin) < -padWidth)
            ? search_area_width -
                  (-padWidth - (origin_x + x_search_area_origin))
            : search_area_width;
    // Correct the right edge of the Search Area if its not on the reference
    // Picture
    x_search_area_origin =
        ((origin_x + x_search_area_origin) >
         (int16_t)quarterRefPicPtr->width - 1)
            ? x_search_area_origin - ((origin_x + x_search_area_origin) -
                                      ((int16_t)quarterRefPicPtr->width - 1))
            : x_search_area_origin;

    search_area_width =
        ((origin_x + x_search_area_origin + search_area_width) >
         (int16_t)quarterRefPicPtr->width)
            ? MAX(1,
                  search_area_width -
                      ((origin_x + x_search_area_origin + search_area_width) -
                       (int16_t)quarterRefPicPtr->width))
            : search_area_width;

    // Constrain x_HME_L1 to be a multiple of 8 (round down as cropping already
    // performed)
    search_area_width =
        (search_area_width < 8) ? search_area_width : search_area_width & ~0x07;
    // Correct the top edge of the Search Area if it is not on the reference
    // Picture
    y_search_area_origin = ((origin_y + y_search_area_origin) < -padHeight)
                               ? -padHeight - origin_y
                               : y_search_area_origin;

    search_area_height =
        ((origin_y + y_search_area_origin) < -padHeight)
            ? search_area_height -
                  (-padHeight - (origin_y + y_search_area_origin))
            : search_area_height;

    // Correct the bottom edge of the Search Area if its not on the reference
    // Picture
    y_search_area_origin =
        ((origin_y + y_search_area_origin) >
         (int16_t)quarterRefPicPtr->height - 1)
            ? y_search_area_origin - ((origin_y + y_search_area_origin) -
                                      ((int16_t)quarterRefPicPtr->height - 1))
            : y_search_area_origin;

    search_area_height =
        (origin_y + y_search_area_origin + search_area_height >
         (int16_t)quarterRefPicPtr->height)
            ? MAX(1,
                  search_area_height -
                      ((origin_y + y_search_area_origin + search_area_height) -
                       (int16_t)quarterRefPicPtr->height))
            : search_area_height;

    // Move to the top left of the search region
    xTopLeftSearchRegion =
        ((int16_t)quarterRefPicPtr->origin_x + origin_x) + x_search_area_origin;
    yTopLeftSearchRegion =
        ((int16_t)quarterRefPicPtr->origin_y + origin_y) + y_search_area_origin;
    searchRegionIndex = xTopLeftSearchRegion +
                        yTopLeftSearchRegion * quarterRefPicPtr->stride_y;

    if (((sb_width & 7) == 0) || (sb_width == 4)) {
        // Put the first search location into level0 results
        sad_loop_kernel(
            &context_ptr->quarter_sb_buffer[0],
            (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                ? context_ptr->quarter_sb_buffer_stride
                : context_ptr->quarter_sb_buffer_stride * 2,
            &quarterRefPicPtr->buffer_y[searchRegionIndex],
            (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                ? quarterRefPicPtr->stride_y
                : quarterRefPicPtr->stride_y * 2,
            (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                ? sb_height
                : sb_height >> 1,
            sb_width,
            /* results */
            level1BestSad,
            xLevel1SearchCenter,
            yLevel1SearchCenter,
            /* range */
            quarterRefPicPtr->stride_y,
            search_area_width,
            search_area_height);
    } else {
        sad_loop_kernel_c(&context_ptr->quarter_sb_buffer[0],
                        (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                            ? context_ptr->quarter_sb_buffer_stride
                            : context_ptr->quarter_sb_buffer_stride * 2,
                        &quarterRefPicPtr->buffer_y[searchRegionIndex],
                        (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                            ? quarterRefPicPtr->stride_y
                            : quarterRefPicPtr->stride_y * 2,
                        (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                            ? sb_height
                            : sb_height >> 1,
                        sb_width,
                        /* results */
                        level1BestSad,
                        xLevel1SearchCenter,
                        yLevel1SearchCenter,
                        /* range */
                        quarterRefPicPtr->stride_y,
                        search_area_width,
                        search_area_height);
    }

    *level1BestSad =
        (context_ptr->hme_search_method == FULL_SAD_SEARCH)
            ? *level1BestSad
            : *level1BestSad *
                  2;  // Multiply by 2 because considered only ever other line
    *xLevel1SearchCenter += x_search_area_origin;
    *xLevel1SearchCenter *=
        2;  // Multiply by 2 because operating on 1/2 resolution
    *yLevel1SearchCenter += y_search_area_origin;
    *yLevel1SearchCenter *=
        2;  // Multiply by 2 because operating on 1/2 resolution

    return;
}

void HmeLevel2(
    PictureParentControlSet
        *picture_control_set_ptr,  // input parameter, Picture control set Ptr
    MeContext *context_ptr,  // input/output parameter, ME context Ptr, used to
                             // get/update ME results
    int16_t
        origin_x,  // input parameter, SB position in the horizontal direction
    int16_t origin_y,  // input parameter, SB position in the vertical direction
    uint32_t sb_width,   // input parameter, SB pwidth - full resolution
    uint32_t sb_height,  // input parameter, SB height - full resolution
    EbPictureBufferDesc *refPicPtr,  // input parameter, reference Picture Ptr
    uint32_t searchRegionNumberInWidth,   // input parameter, search region
                                          // number in the horizontal direction
    uint32_t searchRegionNumberInHeight,  // input parameter, search region
                                          // number in the vertical direction
    int16_t xLevel1SearchCenter,          // input parameter, best Level1 xMV
                                          // at(searchRegionNumberInWidth,
                                          // searchRegionNumberInHeight)
    int16_t yLevel1SearchCenter,          // input parameter, best Level1 yMV
                                          // at(searchRegionNumberInWidth,
                                          // searchRegionNumberInHeight)
    uint64_t *level2BestSad,              // output parameter, Level2 SAD at
                                          // (searchRegionNumberInWidth,
                                          // searchRegionNumberInHeight)
    int16_t *xLevel2SearchCenter,         // output parameter, Level2 xMV at
                                          // (searchRegionNumberInWidth,
                                          // searchRegionNumberInHeight)
    int16_t *yLevel2SearchCenter         // output parameter, Level2 yMV at
                                          // (searchRegionNumberInWidth,
                                          // searchRegionNumberInHeight)
    ) {
    int16_t xTopLeftSearchRegion;
    int16_t yTopLeftSearchRegion;
    uint32_t searchRegionIndex;

    // round the search region width to nearest multiple of 8 if it is less than
    // 8 or non multiple of 8 SAD calculation performance is the same for
    // searchregion width from 1 to 8
    (void)picture_control_set_ptr;
    int16_t hmeLevel2SearchAreaInWidth =
        (int16_t)context_ptr
            ->hme_level2_search_area_in_width_array[searchRegionNumberInWidth];
    // Round up x_HME_L0 to be a multiple of 8
    int16_t search_area_width =
        (int16_t)((hmeLevel2SearchAreaInWidth + 7) & ~0x07);
    int16_t search_area_height =
        (int16_t)context_ptr->hme_level2_search_area_in_height_array
            [searchRegionNumberInHeight];
    int16_t x_search_area_origin;
    int16_t y_search_area_origin;

    int16_t padWidth = (int16_t)BLOCK_SIZE_64 - 1;
    int16_t padHeight = (int16_t)BLOCK_SIZE_64 - 1;

    x_search_area_origin = -(search_area_width >> 1) + xLevel1SearchCenter;
    y_search_area_origin = -(search_area_height >> 1) + yLevel1SearchCenter;

    // Correct the left edge of the Search Area if it is not on the reference
    // Picture
    x_search_area_origin = ((origin_x + x_search_area_origin) < -padWidth)
                               ? -padWidth - origin_x
                               : x_search_area_origin;

    search_area_width =
        ((origin_x + x_search_area_origin) < -padWidth)
            ? search_area_width -
                  (-padWidth - (origin_x + x_search_area_origin))
            : search_area_width;

    // Correct the right edge of the Search Area if its not on the reference
    // Picture
    x_search_area_origin =
        ((origin_x + x_search_area_origin) > (int16_t)refPicPtr->width - 1)
            ? x_search_area_origin - ((origin_x + x_search_area_origin) -
                                      ((int16_t)refPicPtr->width - 1))
            : x_search_area_origin;

    search_area_width =
        ((origin_x + x_search_area_origin + search_area_width) >
         (int16_t)refPicPtr->width)
            ? MAX(1,
                  search_area_width -
                      ((origin_x + x_search_area_origin + search_area_width) -
                       (int16_t)refPicPtr->width))
            : search_area_width;

    // Constrain x_HME_L1 to be a multiple of 8 (round down as cropping already
    // performed)
    search_area_width =
        (search_area_width < 8) ? search_area_width : search_area_width & ~0x07;
    // Correct the top edge of the Search Area if it is not on the reference
    // Picture
    y_search_area_origin = ((origin_y + y_search_area_origin) < -padHeight)
                               ? -padHeight - origin_y
                               : y_search_area_origin;

    search_area_height =
        ((origin_y + y_search_area_origin) < -padHeight)
            ? search_area_height -
                  (-padHeight - (origin_y + y_search_area_origin))
            : search_area_height;

    // Correct the bottom edge of the Search Area if its not on the reference
    // Picture
    y_search_area_origin =
        ((origin_y + y_search_area_origin) > (int16_t)refPicPtr->height - 1)
            ? y_search_area_origin - ((origin_y + y_search_area_origin) -
                                      ((int16_t)refPicPtr->height - 1))
            : y_search_area_origin;

    search_area_height =
        (origin_y + y_search_area_origin + search_area_height >
         (int16_t)refPicPtr->height)
            ? MAX(1,
                  search_area_height -
                      ((origin_y + y_search_area_origin + search_area_height) -
                       (int16_t)refPicPtr->height))
            : search_area_height;

    // Move to the top left of the search region
    xTopLeftSearchRegion =
        ((int16_t)refPicPtr->origin_x + origin_x) + x_search_area_origin;
    yTopLeftSearchRegion =
        ((int16_t)refPicPtr->origin_y + origin_y) + y_search_area_origin;
    searchRegionIndex =
        xTopLeftSearchRegion + yTopLeftSearchRegion * refPicPtr->stride_y;
    if ((((sb_width & 7) == 0) && (sb_width != 40) && (sb_width != 56))) {
        // Put the first search location into level0 results
        sad_loop_kernel(
            context_ptr->sb_src_ptr,
            (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                ? context_ptr->sb_src_stride
                : context_ptr->sb_src_stride * 2,
            &refPicPtr->buffer_y[searchRegionIndex],
            (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                ? refPicPtr->stride_y
                : refPicPtr->stride_y * 2,
            (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                ? sb_height
                : sb_height >> 1,
            sb_width,
            /* results */
            level2BestSad,
            xLevel2SearchCenter,
            yLevel2SearchCenter,
            /* range */
            refPicPtr->stride_y,
            search_area_width,
            search_area_height);
    } else {
        // Put the first search location into level0 results
        sad_loop_kernel_c(context_ptr->sb_src_ptr,
                        (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                            ? context_ptr->sb_src_stride
                            : context_ptr->sb_src_stride * 2,
                        &refPicPtr->buffer_y[searchRegionIndex],
                        (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                            ? refPicPtr->stride_y
                            : refPicPtr->stride_y * 2,
                        (context_ptr->hme_search_method == FULL_SAD_SEARCH)
                            ? sb_height
                            : sb_height >> 1,
                        sb_width,
                        /* results */
                        level2BestSad,
                        xLevel2SearchCenter,
                        yLevel2SearchCenter,
                        /* range */
                        refPicPtr->stride_y,
                        search_area_width,
                        search_area_height);
    }

    *level2BestSad =
        (context_ptr->hme_search_method == FULL_SAD_SEARCH)
            ? *level2BestSad
            : *level2BestSad *
                  2;  // Multiply by 2 because considered only ever other line
    *xLevel2SearchCenter += x_search_area_origin;
    *yLevel2SearchCenter += y_search_area_origin;

    return;
}

static void SelectBuffer(
    uint32_t pu_index,       //[IN]
    uint8_t fracPosition,    //[IN]
    uint32_t pu_width,       //[IN] Refrence picture list index
    uint32_t pu_height,      //[IN] Refrence picture index in the list
    uint8_t *pos_Full,       //[IN]
    uint8_t *pos_b,          //[IN]
    uint8_t *pos_h,          //[IN]
    uint8_t *pos_j,          //[IN]
    uint32_t refHalfStride,  //[IN]
    uint32_t refBufferFullStride,
    uint8_t **dst_ptr,       //[OUT]
    uint32_t *DstPtrStride,  //[OUT]
    EbAsm asm_type) {
    (void)asm_type;
    (void)pu_width;
    (void)pu_height;

    uint32_t puShiftXIndex = pu_search_index_map[pu_index][0];
    uint32_t puShiftYIndex = pu_search_index_map[pu_index][1];
    uint32_t ref_stride = refHalfStride;

    // for each one of the 8 positions, we need to determine the 2 buffers to do
    // averaging
    uint8_t *buf1 = pos_Full;

    switch (fracPosition) {
    case 0:  // integer
        buf1 = pos_Full;
        ref_stride = refBufferFullStride;
        break;
    case 2:  // b
        buf1 = pos_b;
        break;
    case 8:  // h
        buf1 = pos_h;
        break;
    case 10:  // j
        buf1 = pos_j;
        break;
    default: break;
    }

    buf1 = buf1 + puShiftXIndex + puShiftYIndex * ref_stride;

    *dst_ptr = buf1;
    *DstPtrStride = ref_stride;

    return;
}

static void QuarterPelCompensation(
    uint32_t pu_index,       //[IN]
    uint8_t fracPosition,    //[IN]
    uint32_t pu_width,       //[IN] Refrence picture list index
    uint32_t pu_height,      //[IN] Refrence picture index in the list
    uint8_t *pos_Full,       //[IN]
    uint8_t *pos_b,          //[IN]
    uint8_t *pos_h,          //[IN]
    uint8_t *pos_j,          //[IN]
    uint32_t refHalfStride,  //[IN]
    uint32_t refBufferFullStride,
    uint8_t *Dst,        //[IN]
    uint32_t DstStride) { //[IN]

    uint32_t puShiftXIndex = pu_search_index_map[pu_index][0];
    uint32_t puShiftYIndex = pu_search_index_map[pu_index][1];
    uint32_t refStride1 = refHalfStride;
    uint32_t refStride2 = refHalfStride;

    // for each one of the 8 positions, we need to determine the 2 buffers to do
    // averaging
    uint8_t *buf1 = pos_Full;
    uint8_t *buf2 = pos_Full;

    switch (fracPosition) {
    case 1:  // a
        buf1 = pos_Full;
        buf2 = pos_b;
        refStride1 = refBufferFullStride;
        break;

    case 3:  // c
        buf1 = pos_b;
        buf2 = pos_Full + 1;
        refStride2 = refBufferFullStride;
        break;

    case 4:  // d
        buf1 = pos_Full;
        buf2 = pos_h;
        refStride1 = refBufferFullStride;
        break;

    case 5:  // e
        buf1 = pos_b;
        buf2 = pos_h;
        break;

    case 6:  // f
        buf1 = pos_b;
        buf2 = pos_j;
        break;

    case 7:  // g
        buf1 = pos_b;
        buf2 = pos_h + 1;
        break;

    case 9:  // i
        buf1 = pos_h;
        buf2 = pos_j;
        break;

    case 11:  // k
        buf1 = pos_j;
        buf2 = pos_h + 1;
        break;

    case 12:  // L
        buf1 = pos_h;
        buf2 = pos_Full + refBufferFullStride;
        refStride2 = refBufferFullStride;
        break;

    case 13:  // m
        buf1 = pos_h;
        buf2 = pos_b + refHalfStride;
        break;

    case 14:  // n
        buf1 = pos_j;
        buf2 = pos_b + refHalfStride;
        break;
    case 15:  // 0
        buf1 = pos_h + 1;
        buf2 = pos_b + refHalfStride;
        break;
    default: break;
    }

    buf1 = buf1 + puShiftXIndex + puShiftYIndex * refStride1;
    buf2 = buf2 + puShiftXIndex + puShiftYIndex * refStride2;

    picture_average_kernel(buf1,
                            refStride1,
                            buf2,
                            refStride2,
                            Dst,
                            DstStride,
                            pu_width,
                            pu_height);

    return;
}

// TODO: Alt-refs - change previous SelectBuffer and QuarterPelCompensation to
// be applicable for both chroma and luma
static void select_buffer(
    uint32_t pu_index,  //[IN]
    EbBool chroma,
    uint8_t fracPosition,    //[IN]
    uint32_t pu_width,       //[IN] Refrence picture list index
    uint32_t pu_height,      //[IN] Refrence picture index in the list
    uint8_t *pos_Full,       //[IN]
    uint8_t *pos_b,          //[IN]
    uint8_t *pos_h,          //[IN]
    uint8_t *pos_j,          //[IN]
    uint32_t refHalfStride,  //[IN]
    uint32_t refBufferFullStride,
    uint8_t **dst_ptr,       //[OUT]
    uint32_t *DstPtrStride,  //[OUT]
    EbAsm asm_type) {
    (void)asm_type;
    (void)pu_width;
    (void)pu_height;

    uint32_t puShiftXIndex;
    uint32_t puShiftYIndex;

    if (chroma == EB_TRUE) {
        puShiftXIndex = (pu_search_index_map[pu_index][0]) >> 1;
        puShiftYIndex = (pu_search_index_map[pu_index][1]) >> 1;
    } else {
        puShiftXIndex = pu_search_index_map[pu_index][0];
        puShiftYIndex = pu_search_index_map[pu_index][1];
    }

    uint32_t ref_stride = refHalfStride;

    // for each one of the 8 positions, we need to determine the 2 buffers to do
    // averaging
    uint8_t *buf1 = pos_Full;

    switch (fracPosition) {
    case 0:  // integer
        buf1 = pos_Full;
        ref_stride = refBufferFullStride;
        break;
    case 2:  // b
        buf1 = pos_b;
        break;
    case 8:  // h
        buf1 = pos_h;
        break;
    case 10:  // j
        buf1 = pos_j;
        break;
    default: break;
    }

    buf1 = buf1 + puShiftXIndex + puShiftYIndex * ref_stride;

    *dst_ptr = buf1;
    *DstPtrStride = ref_stride;

    return;
}

// TODO: Alt-refs - change previous SelectBuffer and QuarterPelCompensation to
// be applicable for both chroma and luma
static void quarter_pel_compensation(
    uint32_t pu_index,  //[IN]
    EbBool chroma,
    uint8_t fracPosition,    //[IN]
    uint32_t pu_width,       //[IN] Refrence picture list index
    uint32_t pu_height,      //[IN] Refrence picture index in the list
    uint8_t *pos_Full,       //[IN]
    uint8_t *pos_b,          //[IN]
    uint8_t *pos_h,          //[IN]
    uint8_t *pos_j,          //[IN]
    uint32_t refHalfStride,  //[IN]
    uint32_t refBufferFullStride,
    uint8_t *Dst,        //[IN]
    uint32_t DstStride) {  //[IN]
    uint32_t puShiftXIndex;
    uint32_t puShiftYIndex;

    if (chroma == EB_TRUE) {
        puShiftXIndex = (pu_search_index_map[pu_index][0]) >> 1;
        puShiftYIndex = (pu_search_index_map[pu_index][1]) >> 1;
    } else {
        puShiftXIndex = pu_search_index_map[pu_index][0];
        puShiftYIndex = pu_search_index_map[pu_index][1];
    }

    uint32_t refStride1 = refHalfStride;
    uint32_t refStride2 = refHalfStride;

    // for each one of the 8 positions, we need to determine the 2 buffers to do
    // averaging
    uint8_t *buf1 = pos_Full;
    uint8_t *buf2 = pos_Full;

    switch (fracPosition) {
    case 1:  // a
        buf1 = pos_Full;
        buf2 = pos_b;
        refStride1 = refBufferFullStride;
        break;

    case 3:  // c
        buf1 = pos_b;
        buf2 = pos_Full + 1;
        refStride2 = refBufferFullStride;
        break;

    case 4:  // d
        buf1 = pos_Full;
        buf2 = pos_h;
        refStride1 = refBufferFullStride;
        break;

    case 5:  // e
        buf1 = pos_b;
        buf2 = pos_h;
        break;

    case 6:  // f
        buf1 = pos_b;
        buf2 = pos_j;
        break;

    case 7:  // g
        buf1 = pos_b;
        buf2 = pos_h + 1;
        break;

    case 9:  // i
        buf1 = pos_h;
        buf2 = pos_j;
        break;

    case 11:  // k
        buf1 = pos_j;
        buf2 = pos_h + 1;
        break;

    case 12:  // L
        buf1 = pos_h;
        buf2 = pos_Full + refBufferFullStride;
        refStride2 = refBufferFullStride;
        break;

    case 13:  // m
        buf1 = pos_h;
        buf2 = pos_b + refHalfStride;
        break;

    case 14:  // n
        buf1 = pos_j;
        buf2 = pos_b + refHalfStride;
        break;
    case 15:  // 0
        buf1 = pos_h + 1;
        buf2 = pos_b + refHalfStride;
        break;
    default: break;
    }

    buf1 = buf1 + puShiftXIndex + puShiftYIndex * refStride1;
    buf2 = buf2 + puShiftXIndex + puShiftYIndex * refStride2;

    picture_average_kernel(buf1,
                            refStride1,
                            buf2,
                            refStride2,
                            Dst,
                            DstStride,
                            pu_width,
                            pu_height);

    return;
}

/*******************************************************************************
 * Requirement: pu_width      = 8, 16, 24, 32, 48 or 64
 * Requirement: pu_height % 2 = 0
 * Requirement: skip         = 0 or 1
 * Requirement (x86 only): temp_buf % 16 = 0
 * Requirement (x86 only): (dst->buffer_y  + dstLumaIndex  ) % 16 = 0 when
 *pu_width %16 = 0 Requirement (x86 only): (dst->bufferCb + dstChromaIndex) % 16
 *= 0 when pu_width %32 = 0 Requirement (x86 only): (dst->bufferCr +
 *dstChromaIndex) % 16 = 0 when pu_width %32 = 0 Requirement (x86 only):
 *dst->stride_y   % 16 = 0 when pu_width %16 = 0 Requirement (x86 only):
 *dst->chromaStride % 16 = 0 when pu_width %32 = 0
 *******************************************************************************/
void uni_pred_averaging(uint32_t pu_index, EbBool chroma, uint8_t firstFracPos,
                        uint32_t pu_width, uint32_t pu_height,
                        uint8_t *firstRefInteger, uint8_t *firstRefPosB,
                        uint8_t *firstRefPosH, uint8_t *firstRefPosJ,
                        uint32_t refBufferStride,
                        uint32_t refBufferFullList0Stride,
                        uint8_t *firstRefTempDst, uint8_t **comp_blk_ptr,
                        uint32_t *comp_blk_ptr_stride, EbAsm asm_type) {
    // Buffer Selection and quater-pel compensation on the fly
    if (sub_position_type[firstFracPos] != 2) {
        select_buffer(pu_index,
                      chroma,
                      firstFracPos,
                      pu_width,
                      pu_height,
                      firstRefInteger,
                      firstRefPosB,
                      firstRefPosH,
                      firstRefPosJ,
                      refBufferStride,
                      refBufferFullList0Stride,
                      comp_blk_ptr,
                      comp_blk_ptr_stride,
                      asm_type);
    } else {
        quarter_pel_compensation(pu_index,
                                 chroma,
                                 firstFracPos,
                                 pu_width,
                                 pu_height,
                                 firstRefInteger,
                                 firstRefPosB,
                                 firstRefPosH,
                                 firstRefPosJ,
                                 refBufferStride,
                                 refBufferFullList0Stride,
                                 firstRefTempDst,
                                 BLOCK_SIZE_64);

        *comp_blk_ptr = firstRefTempDst;
        *comp_blk_ptr_stride = BLOCK_SIZE_64;
    }
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
 *dst->chromaStride % 16 = 0 when pu_width %32 = 0
 *******************************************************************************/
uint32_t BiPredAverging(
    MeContext *context_ptr, MePredUnit *me_candidate, uint32_t pu_index,
    uint8_t *sourcePic, uint32_t lumaStride, uint8_t firstFracPos,
    uint8_t secondFracPos, uint32_t pu_width, uint32_t pu_height,
    uint8_t *firstRefInteger, uint8_t *firstRefPosB, uint8_t *firstRefPosH,
    uint8_t *firstRefPosJ, uint8_t *secondRefInteger, uint8_t *secondRefPosB,
    uint8_t *secondRefPosH, uint8_t *secondRefPosJ, uint32_t refBufferStride,
    uint32_t refBufferFullList0Stride, uint32_t refBufferFullList1Stride,
    uint8_t *firstRefTempDst, uint8_t *secondRefTempDst, EbAsm asm_type) {
    uint8_t *ptrList0, *ptrList1;
    uint32_t ptrList0Stride, ptrList1Stride;

    // Buffer Selection and quater-pel compensation on the fly
    if (sub_position_type[firstFracPos] != 2) {
        SelectBuffer(pu_index,
                     firstFracPos,
                     pu_width,
                     pu_height,
                     firstRefInteger,
                     firstRefPosB,
                     firstRefPosH,
                     firstRefPosJ,
                     refBufferStride,
                     refBufferFullList0Stride,
                     &ptrList0,
                     &ptrList0Stride,
                     asm_type);
    } else {
        QuarterPelCompensation(pu_index,
                               firstFracPos,
                               pu_width,
                               pu_height,
                               firstRefInteger,
                               firstRefPosB,
                               firstRefPosH,
                               firstRefPosJ,
                               refBufferStride,
                               refBufferFullList0Stride,
                               firstRefTempDst,
                               BLOCK_SIZE_64);

        ptrList0 = firstRefTempDst;
        ptrList0Stride = BLOCK_SIZE_64;
    }

    if (sub_position_type[secondFracPos] != 2) {
        SelectBuffer(pu_index,
                     secondFracPos,
                     pu_width,
                     pu_height,
                     secondRefInteger,
                     secondRefPosB,
                     secondRefPosH,
                     secondRefPosJ,
                     refBufferStride,
                     refBufferFullList1Stride,
                     &ptrList1,
                     &ptrList1Stride,
                     asm_type);
    } else {
        // uni-prediction List1 luma
        // doing the luma interpolation
        QuarterPelCompensation(pu_index,
                               secondFracPos,
                               pu_width,
                               pu_height,
                               secondRefInteger,
                               secondRefPosB,
                               secondRefPosH,
                               secondRefPosJ,
                               refBufferStride,
                               refBufferFullList1Stride,
                               secondRefTempDst,
                               BLOCK_SIZE_64);

        ptrList1 = secondRefTempDst;
        ptrList1Stride = BLOCK_SIZE_64;
    }

    // bi-pred luma
    me_candidate->distortion =
        (context_ptr->fractional_search_method == SUB_SAD_SEARCH)
            ? nxm_sad_avg_kernel(
                  sourcePic,
                  lumaStride << 1,
                  ptrList0,
                  ptrList0Stride << 1,
                  ptrList1,
                  ptrList1Stride << 1,
                  pu_height >> 1,
                  pu_width)
                  << 1
            : nxm_sad_avg_kernel(
                  sourcePic,
                  lumaStride,
                  ptrList0,
                  ptrList0Stride,
                  ptrList1,
                  ptrList1Stride,
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
EbErrorType BiPredictionCompensation(MeContext *context_ptr, uint32_t pu_index,
                                     MePredUnit *me_candidate,
                                     uint32_t firstList,
                                     uint8_t first_list_ref_pic_idx,
                                     uint32_t firstRefMv, uint32_t secondList,
                                     uint8_t second_list_ref_pic_idx,
                                     uint32_t secondRefMv, EbAsm asm_type) {
    EbErrorType return_error = EB_ErrorNone;

    int16_t firstRefPosX;
    int16_t firstRefPosY;
    int16_t firstRefIntegPosx;
    int16_t firstRefIntegPosy;
    uint8_t firstRefFracPosx;
    uint8_t firstRefFracPosy;
    uint8_t firstRefFracPos;
    int32_t xfirstSearchIndex;
    int32_t yfirstSearchIndex;
    int32_t firstSearchRegionIndexPosInteg;
    int32_t firstSearchRegionIndexPosb;
    int32_t firstSearchRegionIndexPosh;
    int32_t firstSearchRegionIndexPosj;

    int16_t secondRefPosX;
    int16_t secondRefPosY;
    int16_t secondRefIntegPosx;
    int16_t secondRefIntegPosy;
    uint8_t secondRefFracPosx;
    uint8_t secondRefFracPosy;
    uint8_t secondRefFracPos;
    int32_t xsecondSearchIndex;
    int32_t ysecondSearchIndex;
    int32_t secondSearchRegionIndexPosInteg;
    int32_t secondSearchRegionIndexPosb;
    int32_t secondSearchRegionIndexPosh;
    int32_t secondSearchRegionIndexPosj;

    uint32_t puShiftXIndex = pu_search_index_map[pu_index][0];
    uint32_t puShiftYIndex = pu_search_index_map[pu_index][1];

    const uint32_t puLcuBufferIndex =
        puShiftXIndex + puShiftYIndex * context_ptr->sb_src_stride;

    me_candidate->prediction_direction = BI_PRED;

    // First refrence
    // Set Candidate information
    firstRefPosX = _MVXT(firstRefMv);
    firstRefPosY = _MVYT(firstRefMv);
    me_candidate->ref_index[0] = (uint8_t)first_list_ref_pic_idx;
    me_candidate->ref0_list = (uint8_t)firstList;

    firstRefIntegPosx = (firstRefPosX >> 2);
    firstRefIntegPosy = (firstRefPosY >> 2);
    firstRefFracPosx = (uint8_t)firstRefPosX & 0x03;
    firstRefFracPosy = (uint8_t)firstRefPosY & 0x03;

    firstRefFracPos = (uint8_t)(firstRefFracPosx + (firstRefFracPosy << 2));
    xfirstSearchIndex =
        (int32_t)firstRefIntegPosx -
        context_ptr->x_search_area_origin[firstList][first_list_ref_pic_idx];
    yfirstSearchIndex =
        (int32_t)firstRefIntegPosy -
        context_ptr->y_search_area_origin[firstList][first_list_ref_pic_idx];
    firstSearchRegionIndexPosInteg =
        (int32_t)(xfirstSearchIndex + (ME_FILTER_TAP >> 1)) +
        (int32_t)context_ptr
                ->interpolated_full_stride[firstList][first_list_ref_pic_idx] *
            (int32_t)(yfirstSearchIndex + (ME_FILTER_TAP >> 1));

    firstSearchRegionIndexPosb =
        (int32_t)(xfirstSearchIndex + (ME_FILTER_TAP >> 1) - 1) +
        (int32_t)context_ptr->interpolated_stride *
            (int32_t)(yfirstSearchIndex + (ME_FILTER_TAP >> 1));
    firstSearchRegionIndexPosh =
        (int32_t)(xfirstSearchIndex + (ME_FILTER_TAP >> 1) - 1) +
        (int32_t)context_ptr->interpolated_stride *
            (int32_t)(yfirstSearchIndex + (ME_FILTER_TAP >> 1) - 1);
    firstSearchRegionIndexPosj =
        (int32_t)(xfirstSearchIndex + (ME_FILTER_TAP >> 1) - 1) +
        (int32_t)context_ptr->interpolated_stride *
            (int32_t)(yfirstSearchIndex + (ME_FILTER_TAP >> 1) - 1);

    // Second refrence

    // Set Candidate information
    secondRefPosX = _MVXT(secondRefMv);
    secondRefPosY = _MVYT(secondRefMv);
    me_candidate->ref_index[1] = (uint8_t)second_list_ref_pic_idx;
    me_candidate->ref1_list = (uint8_t)secondList;
    secondRefIntegPosx = (secondRefPosX >> 2);
    secondRefIntegPosy = (secondRefPosY >> 2);
    secondRefFracPosx = (uint8_t)secondRefPosX & 0x03;
    secondRefFracPosy = (uint8_t)secondRefPosY & 0x03;

    secondRefFracPos = (uint8_t)(secondRefFracPosx + (secondRefFracPosy << 2));
    xsecondSearchIndex =
        secondRefIntegPosx -
        context_ptr->x_search_area_origin[secondList][second_list_ref_pic_idx];
    ysecondSearchIndex =
        secondRefIntegPosy -
        context_ptr->y_search_area_origin[secondList][second_list_ref_pic_idx];
    secondSearchRegionIndexPosInteg =
        (int32_t)(xsecondSearchIndex + (ME_FILTER_TAP >> 1)) +
        (int32_t)
                context_ptr->interpolated_full_stride[secondList]
                                                     [second_list_ref_pic_idx] *
            (int32_t)(ysecondSearchIndex + (ME_FILTER_TAP >> 1));
    secondSearchRegionIndexPosb =
        (int32_t)(xsecondSearchIndex + (ME_FILTER_TAP >> 1) - 1) +
        (int32_t)context_ptr->interpolated_stride *
            (int32_t)(ysecondSearchIndex + (ME_FILTER_TAP >> 1));
    secondSearchRegionIndexPosh =
        (int32_t)(xsecondSearchIndex + (ME_FILTER_TAP >> 1) - 1) +
        (int32_t)context_ptr->interpolated_stride *
            (int32_t)(ysecondSearchIndex + (ME_FILTER_TAP >> 1) - 1);
    secondSearchRegionIndexPosj =
        (int32_t)(xsecondSearchIndex + (ME_FILTER_TAP >> 1) - 1) +
        (int32_t)context_ptr->interpolated_stride *
            (int32_t)(ysecondSearchIndex + (ME_FILTER_TAP >> 1) - 1);

    uint32_t nIndex;

    if (pu_index > 200)
        nIndex = pu_index;
    else if (pu_index > 184)
        nIndex = tab8x32[pu_index - 185] + 185;
    else if (pu_index > 168)
        nIndex = tab32x8[pu_index - 169] + 169;
    else if (pu_index > 136)
        nIndex = tab8x16[pu_index - 137] + 137;
    else if (pu_index > 128)
        nIndex = tab16x32[pu_index - 129] + 129;
    else if (pu_index > 126)
        nIndex = pu_index;
    else if (pu_index > 94)
        nIndex = tab16x8[pu_index - 95] + 95;
    else if (pu_index > 86)
        nIndex = tab32x16[pu_index - 87] + 87;
    else if (pu_index > 84)
        nIndex = pu_index;
    else if (pu_index > 20)
        nIndex = tab8x8[pu_index - 21] + 21;
    else if (pu_index > 4)
        nIndex = tab16x16[pu_index - 5] + 5;
    else
        nIndex = pu_index;
    context_ptr->p_sb_bipred_sad[nIndex] =

        BiPredAverging(
            context_ptr,
            me_candidate,
            pu_index,
            &(context_ptr->sb_src_ptr[puLcuBufferIndex]),
            context_ptr->sb_src_stride,
            firstRefFracPos,
            secondRefFracPos,
            partition_width[pu_index],
            partition_height[pu_index],
            &(context_ptr->integer_buffer_ptr[firstList][first_list_ref_pic_idx]
                                             [firstSearchRegionIndexPosInteg]),
            &(context_ptr->pos_b_buffer[firstList][first_list_ref_pic_idx]
                                       [firstSearchRegionIndexPosb]),
            &(context_ptr->pos_h_buffer[firstList][first_list_ref_pic_idx]
                                       [firstSearchRegionIndexPosh]),
            &(context_ptr->pos_j_buffer[firstList][first_list_ref_pic_idx]
                                       [firstSearchRegionIndexPosj]),
            &(context_ptr
                  ->integer_buffer_ptr[secondList][second_list_ref_pic_idx]
                                      [secondSearchRegionIndexPosInteg]),
            &(context_ptr->pos_b_buffer[secondList][second_list_ref_pic_idx]
                                       [secondSearchRegionIndexPosb]),
            &(context_ptr->pos_h_buffer[secondList][second_list_ref_pic_idx]
                                       [secondSearchRegionIndexPosh]),
            &(context_ptr->pos_j_buffer[secondList][second_list_ref_pic_idx]
                                       [secondSearchRegionIndexPosj]),
            context_ptr->interpolated_stride,
            context_ptr
                ->interpolated_full_stride[firstList][first_list_ref_pic_idx],
            context_ptr
                ->interpolated_full_stride[secondList][second_list_ref_pic_idx],
            &(context_ptr->one_d_intermediate_results_buf0[0]),
            &(context_ptr->one_d_intermediate_results_buf1[0]),
            asm_type);

    return return_error;
}

uint8_t skip_bi_pred(
    PictureParentControlSet *picture_control_set_ptr,
    uint8_t ref_type,
    uint8_t ref_type_table[7]) {

    if (!picture_control_set_ptr->prune_unipred_at_me)
        return 1;

    uint8_t allow_cand = 0;
    uint8_t ref_idx;
    for (ref_idx = 0; ref_idx < PRUNE_REF_ME_TH; ref_idx++) {
        if (ref_type == ref_type_table[ref_idx])
            allow_cand = 1;
    }
    return allow_cand;
}

/*******************************************
 * BiPredictionSearch
 *   performs Bi-Prediction Search (LCU)
 *******************************************/
// This function enables all 16 Bipred candidates when MRP is ON
EbErrorType BiPredictionSearch(
    SequenceControlSet *sequence_control_set_ptr,
    MeContext *context_ptr, uint32_t pu_index, uint8_t candidateIndex,
    uint32_t activeRefPicFirstLisNum, uint32_t activeRefPicSecondLisNum,
    uint8_t *total_me_candidate_index, EbAsm asm_type,
    uint8_t ref_type_table[7],
    PictureParentControlSet *picture_control_set_ptr) {
    EbErrorType return_error = EB_ErrorNone;

    uint32_t firstListRefPictdx;
    uint32_t secondListRefPictdx;

    (void)picture_control_set_ptr;

    uint32_t nIndex;

    if (pu_index > 200)
        nIndex = pu_index;
    else if (pu_index > 184)
        nIndex = tab8x32[pu_index - 185] + 185;
    else if (pu_index > 168)
        nIndex = tab32x8[pu_index - 169] + 169;
    else if (pu_index > 136)
        nIndex = tab8x16[pu_index - 137] + 137;
    else if (pu_index > 128)
        nIndex = tab16x32[pu_index - 129] + 129;
    else if (pu_index > 126)
        nIndex = pu_index;
    else if (pu_index > 94)
        nIndex = tab16x8[pu_index - 95] + 95;
    else if (pu_index > 86)
        nIndex = tab32x16[pu_index - 87] + 87;
    else if (pu_index > 84)
        nIndex = pu_index;
    else if (pu_index > 20)
        nIndex = tab8x8[pu_index - 21] + 21;
    else if (pu_index > 4)
        nIndex = tab16x16[pu_index - 5] + 5;
    else
        nIndex = pu_index;
    // NM: Inter list bipred.
    //(LAST,BWD) , (LAST,ALT)  and (LAST,ALT2)
    //(LAST2,BWD), (LAST2,ALT) and (LAST2,ALT2)
    //(LAST3,BWD), (LAST3,ALT) and (LAST3,ALT2)
    //(GOLD,BWD) , (GOLD,ALT)  and (GOLD,ALT2)
    for (firstListRefPictdx = 0; firstListRefPictdx < activeRefPicFirstLisNum;
         firstListRefPictdx++) {
        for (secondListRefPictdx = 0;
             secondListRefPictdx < activeRefPicSecondLisNum;
             secondListRefPictdx++) {
            {
                     uint8_t to_inject_ref_type_0 = svt_get_ref_frame_type(REF_LIST_0, firstListRefPictdx);
                     uint8_t to_inject_ref_type_1 = svt_get_ref_frame_type(REF_LIST_1, secondListRefPictdx);
                     uint8_t add_bi = skip_bi_pred(
                         picture_control_set_ptr,
                         to_inject_ref_type_0,
                         ref_type_table);
                     add_bi += skip_bi_pred(
                         picture_control_set_ptr,
                         to_inject_ref_type_1,
                         ref_type_table);

                     if (add_bi) {
                BiPredictionCompensation(
                    context_ptr,
                    pu_index,
                    &(context_ptr->me_candidate[candidateIndex].pu[pu_index]),
                    REFERENCE_PIC_LIST_0,
                    firstListRefPictdx,
                    context_ptr->p_sb_best_mv[REFERENCE_PIC_LIST_0]
                                             [firstListRefPictdx][nIndex],
                    REFERENCE_PIC_LIST_1,
                    secondListRefPictdx,
                    context_ptr->p_sb_best_mv[REFERENCE_PIC_LIST_1]
                                             [secondListRefPictdx][nIndex],
                    asm_type);

                candidateIndex++;
                     }
            }
        }
    }

    if (sequence_control_set_ptr->mrp_mode == 0)
    {
        // NM: Within list 0    bipred: (LAST,LAST2)    (LAST,LAST3) (LAST,GOLD)
        for (firstListRefPictdx = 1;
             firstListRefPictdx < activeRefPicFirstLisNum;
             firstListRefPictdx++) {
            uint8_t to_inject_ref_type_0 = svt_get_ref_frame_type(REF_LIST_0, firstListRefPictdx);
            uint8_t add_bi = skip_bi_pred(
                picture_control_set_ptr,
                to_inject_ref_type_0,
                ref_type_table);
            if (add_bi) {
            BiPredictionCompensation(
                context_ptr,
                pu_index,
                &(context_ptr->me_candidate[candidateIndex].pu[pu_index]),
                REFERENCE_PIC_LIST_0,
                0,
                context_ptr->p_sb_best_mv[REFERENCE_PIC_LIST_0][0][nIndex],
                REFERENCE_PIC_LIST_0,
                firstListRefPictdx,
                context_ptr->p_sb_best_mv[REFERENCE_PIC_LIST_0]
                                         [firstListRefPictdx][nIndex],
                asm_type);

            candidateIndex++;
            }
        }
        // NM: Within list 1    bipred: (BWD, ALT)
        for (secondListRefPictdx = 1;
             secondListRefPictdx < MIN(activeRefPicSecondLisNum, 1);
             secondListRefPictdx++) {
            uint8_t to_inject_ref_type_0 = svt_get_ref_frame_type(REF_LIST_0, firstListRefPictdx);
            uint8_t add_bi = skip_bi_pred(
                picture_control_set_ptr,
                to_inject_ref_type_0,
                ref_type_table);
            if (add_bi) {
            BiPredictionCompensation(
                context_ptr,
                pu_index,
                &(context_ptr->me_candidate[candidateIndex].pu[pu_index]),
                REFERENCE_PIC_LIST_1,
                0,
                context_ptr->p_sb_best_mv[REFERENCE_PIC_LIST_1][0][nIndex],
                REFERENCE_PIC_LIST_1,
                secondListRefPictdx,
                context_ptr->p_sb_best_mv[REFERENCE_PIC_LIST_1]
                                         [secondListRefPictdx][nIndex],
                asm_type);

            candidateIndex++;
            }
        }
    }
    *total_me_candidate_index = candidateIndex;

    return return_error;
}

// Nader - to be replaced by loock-up table
/*******************************************
 * get_me_info_index
 *   search the correct index of the motion
 *   info that corresponds to the input
 *   md candidate
 *******************************************/
uint32_t get_me_info_index(uint32_t max_me_block, const BlockGeom *blk_geom,
                           uint32_t geom_offset_x, uint32_t geom_offset_y) {
    // search for motion info
    uint32_t block_index;
    uint32_t me_info_index = 0xFFFFFFF;

    for (block_index = 0; block_index < max_me_block; block_index++) {
        if ((blk_geom->bwidth == partition_width[block_index]) &&
            (blk_geom->bheight == partition_height[block_index]) &&
            ((blk_geom->origin_x - geom_offset_x) ==
             pu_search_index_map[block_index][0]) &&
            ((blk_geom->origin_y - geom_offset_y) ==
             pu_search_index_map[block_index][1])) {
            me_info_index = block_index;
            break;
        }
    }
    return me_info_index;
}

// Nader - to be replaced by loock-up table
/*******************************************
 * get_me_info_index
 *   search the correct index of the motion
 *   info that corresponds to the input
 *   md candidate
 *******************************************/
uint32_t get_in_loop_me_info_index(uint32_t max_me_block, uint8_t is_128_sb,
                                   const BlockGeom *blk_geom) {
    // search for motion info
    uint32_t block_index;
    uint32_t me_info_index = 0xFFFFFFF;
    if (is_128_sb) {
        for (block_index = 0; block_index < max_me_block; block_index++) {
            if (blk_geom->bwidth ==
                    in_loop_me_block_width_128_sb[block_index] &&
                blk_geom->bheight ==
                    in_loop_me_block_height_128_sb[block_index] &&
                blk_geom->origin_x ==
                    in_loop_me_block_index_128_sb[block_index][0] &&
                blk_geom->origin_y ==
                    in_loop_me_block_index_128_sb[block_index][1]) {
                me_info_index = block_index;
                break;
            }
        }
    } else {
        for (block_index = 0; block_index < max_me_block; block_index++) {
            if (blk_geom->bwidth == in_loop_me_block_width[block_index] &&
                blk_geom->bheight == in_loop_me_block_height[block_index] &&
                blk_geom->origin_x == in_loop_me_block_index[block_index][0] &&
                blk_geom->origin_y == in_loop_me_block_index[block_index][1]) {
                me_info_index = block_index;
                break;
            }
        }
    }

    return me_info_index;
}

#define NSET_CAND(mePuResult, num, dist, dir)                      \
    (mePuResult)->distortion_direction[(num)].distortion = (dist); \
    (mePuResult)->distortion_direction[(num)].direction = (dir);

int8_t sort_3_elements(uint32_t a, uint32_t b, uint32_t c) {
    uint8_t sortCode = 0;
    if (a <= b && a <= c) {
        if (b <= c)
            sortCode = a_b_c;
        else
            sortCode = a_c_b;
    } else if (b <= a && b <= c) {
        if (a <= c)
            sortCode = b_a_c;
        else
            sortCode = b_c_a;
    } else if (a <= b)
        sortCode = c_a_b;
    else
        sortCode = c_b_a;
    return sortCode;
}

EbErrorType CheckZeroZeroCenter(EbPictureBufferDesc *refPicPtr,
                                MeContext *context_ptr, uint32_t sb_origin_x,
                                uint32_t sb_origin_y, uint32_t sb_width,
                                uint32_t sb_height, int16_t *x_search_center,
                                int16_t *y_search_center)

{
    EbErrorType return_error = EB_ErrorNone;
    uint32_t searchRegionIndex, zeroMvSad, hmeMvSad, hmeMvdRate;
    uint64_t hmeMvCost, zeroMvCost, searchCenterCost;
    int16_t origin_x = (int16_t)sb_origin_x;
    int16_t origin_y = (int16_t)sb_origin_y;
    uint32_t subsampleSad = 1;
    int16_t pad_width = (int16_t)BLOCK_SIZE_64 - 1;
    int16_t pad_height = (int16_t)BLOCK_SIZE_64 - 1;

    searchRegionIndex =
        (int16_t)refPicPtr->origin_x + origin_x +
        ((int16_t)refPicPtr->origin_y + origin_y) * refPicPtr->stride_y;

    zeroMvSad = nxm_sad_kernel(
        context_ptr->sb_src_ptr,
        context_ptr->sb_src_stride << subsampleSad,
        &(refPicPtr->buffer_y[searchRegionIndex]),
        refPicPtr->stride_y << subsampleSad,
        sb_height >> subsampleSad,
        sb_width);

    zeroMvSad = zeroMvSad << subsampleSad;

    // FIX
    // Correct the left edge of the Search Area if it is not on the reference
    // Picture
    *x_search_center = ((origin_x + *x_search_center) < -pad_width)
                           ? -pad_width - origin_x
                           : *x_search_center;
    // Correct the right edge of the Search Area if its not on the reference
    // Picture
    *x_search_center =
        ((origin_x + *x_search_center) > (int16_t)refPicPtr->width - 1)
            ? *x_search_center - ((origin_x + *x_search_center) -
                                  ((int16_t)refPicPtr->width - 1))
            : *x_search_center;
    // Correct the top edge of the Search Area if it is not on the reference
    // Picture
    *y_search_center = ((origin_y + *y_search_center) < -pad_height)
                           ? -pad_height - origin_y
                           : *y_search_center;
    // Correct the bottom edge of the Search Area if its not on the reference
    // Picture
    *y_search_center =
        ((origin_y + *y_search_center) > (int16_t)refPicPtr->height - 1)
            ? *y_search_center - ((origin_y + *y_search_center) -
                                  ((int16_t)refPicPtr->height - 1))
            : *y_search_center;
    ///

    zeroMvCost = zeroMvSad << COST_PRECISION;
    searchRegionIndex =
        (int16_t)(refPicPtr->origin_x + origin_x) + *x_search_center +
        ((int16_t)(refPicPtr->origin_y + origin_y) + *y_search_center) *
            refPicPtr->stride_y;

    hmeMvSad = nxm_sad_kernel(
        context_ptr->sb_src_ptr,
        context_ptr->sb_src_stride << subsampleSad,
        &(refPicPtr->buffer_y[searchRegionIndex]),
        refPicPtr->stride_y << subsampleSad,
        sb_height >> subsampleSad,
        sb_width);

    hmeMvSad = hmeMvSad << subsampleSad;

    hmeMvdRate = 0;
    // AMIR use AV1 rate estimation functions
    // MeGetMvdFractionBits(
    //    ABS(*x_search_center << 2),
    //    ABS(*y_search_center << 2),
    //    context_ptr->mvd_bits_array,
    //    &hmeMvdRate);

    hmeMvCost = (hmeMvSad << COST_PRECISION) +
                (((context_ptr->lambda * hmeMvdRate) + MD_OFFSET) >> MD_SHIFT);
    searchCenterCost = MIN(zeroMvCost, hmeMvCost);

    *x_search_center = (searchCenterCost == zeroMvCost) ? 0 : *x_search_center;
    *y_search_center = (searchCenterCost == zeroMvCost) ? 0 : *y_search_center;

    return return_error;
}

EbErrorType suPelEnable(MeContext *context_ptr,
                        PictureParentControlSet *picture_control_set_ptr,
                        uint32_t listIndex, uint32_t refPicIndex,
                        EbBool *enableHalfPel32x32, EbBool *enableHalfPel16x16,
                        EbBool *enableHalfPel8x8) {
    EbErrorType return_error = EB_ErrorNone;

    uint32_t mvMag32x32 = 0;
    uint32_t mvMag16x16 = 0;
    uint32_t mvMag8x8 = 0;
    uint32_t avgSad32x32 = 0;
    uint32_t avgSad16x16 = 0;
    uint32_t avgSad8x8 = 0;
    uint32_t avgMvx32x32 = 0;
    uint32_t avgMvy32x32 = 0;
    uint32_t avgMvx16x16 = 0;
    uint32_t avgMvy16x16 = 0;
    uint32_t avgMvx8x8 = 0;
    uint32_t avgMvy8x8 = 0;

    avgMvx32x32 = (_MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                  [ME_TIER_ZERO_PU_32x32_0]) +
                   _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                  [ME_TIER_ZERO_PU_32x32_1]) +
                   _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                  [ME_TIER_ZERO_PU_32x32_2]) +
                   _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                  [ME_TIER_ZERO_PU_32x32_3])) >>
                  2;
    avgMvy32x32 = (_MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                  [ME_TIER_ZERO_PU_32x32_0]) +
                   _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                  [ME_TIER_ZERO_PU_32x32_1]) +
                   _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                  [ME_TIER_ZERO_PU_32x32_2]) +
                   _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                  [ME_TIER_ZERO_PU_32x32_3])) >>
                  2;
    mvMag32x32 = SQR(avgMvx32x32) + SQR(avgMvy32x32);

    avgMvx16x16 =
        (_MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                        [ME_TIER_ZERO_PU_16x16_0]) +
         _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                        [ME_TIER_ZERO_PU_16x16_1]) +
         _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                        [ME_TIER_ZERO_PU_16x16_2]) +
         _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                        [ME_TIER_ZERO_PU_16x16_3]) +
         _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                        [ME_TIER_ZERO_PU_16x16_4]) +
         _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                        [ME_TIER_ZERO_PU_16x16_5]) +
         _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                        [ME_TIER_ZERO_PU_16x16_6]) +
         _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                        [ME_TIER_ZERO_PU_16x16_7]) +
         _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                        [ME_TIER_ZERO_PU_16x16_8]) +
         _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                        [ME_TIER_ZERO_PU_16x16_9]) +
         _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                        [ME_TIER_ZERO_PU_16x16_10]) +
         _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                        [ME_TIER_ZERO_PU_16x16_11]) +
         _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                        [ME_TIER_ZERO_PU_16x16_12]) +
         _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                        [ME_TIER_ZERO_PU_16x16_13]) +
         _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                        [ME_TIER_ZERO_PU_16x16_14]) +
         _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                        [ME_TIER_ZERO_PU_16x16_15])) >>
        4;
    avgMvy16x16 =
        (_MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                        [ME_TIER_ZERO_PU_16x16_0]) +
         _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                        [ME_TIER_ZERO_PU_16x16_1]) +
         _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                        [ME_TIER_ZERO_PU_16x16_2]) +
         _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                        [ME_TIER_ZERO_PU_16x16_3]) +
         _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                        [ME_TIER_ZERO_PU_16x16_4]) +
         _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                        [ME_TIER_ZERO_PU_16x16_5]) +
         _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                        [ME_TIER_ZERO_PU_16x16_6]) +
         _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                        [ME_TIER_ZERO_PU_16x16_7]) +
         _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                        [ME_TIER_ZERO_PU_16x16_8]) +
         _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                        [ME_TIER_ZERO_PU_16x16_9]) +
         _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                        [ME_TIER_ZERO_PU_16x16_10]) +
         _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                        [ME_TIER_ZERO_PU_16x16_11]) +
         _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                        [ME_TIER_ZERO_PU_16x16_12]) +
         _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                        [ME_TIER_ZERO_PU_16x16_13]) +
         _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                        [ME_TIER_ZERO_PU_16x16_14]) +
         _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                        [ME_TIER_ZERO_PU_16x16_15])) >>
        4;
    mvMag16x16 = SQR(avgMvx16x16) + SQR(avgMvy16x16);

    avgMvx8x8 = (_MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_0]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_1]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_2]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_3]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_4]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_5]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_6]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_7]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_8]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_9]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_10]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_11]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_12]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_13]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_14]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_15]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_16]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_17]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_18]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_19]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_20]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_21]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_22]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_23]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_24]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_25]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_26]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_27]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_28]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_29]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_30]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_31]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_32]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_33]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_34]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_35]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_36]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_37]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_38]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_39]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_40]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_41]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_42]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_43]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_44]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_45]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_46]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_47]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_48]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_49]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_50]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_51]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_52]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_53]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_54]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_55]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_56]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_57]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_58]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_59]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_60]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_61]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_62]) +
                 _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_63])) >>
                6;
    avgMvy8x8 = (_MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_0]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_1]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_2]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_3]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_4]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_5]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_6]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_7]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_8]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_9]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_10]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_11]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_12]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_13]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_14]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_15]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_16]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_17]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_18]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_19]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_20]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_21]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_22]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_23]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_24]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_25]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_26]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_27]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_28]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_29]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_30]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_31]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_32]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_33]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_34]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_35]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_36]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_37]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_38]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_39]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_40]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_41]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_42]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_43]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_44]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_45]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_46]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_47]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_48]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_49]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_50]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_51]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_52]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_53]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_54]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_55]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_56]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_57]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_58]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_59]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_60]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_61]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_62]) +
                 _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex]
                                                [ME_TIER_ZERO_PU_8x8_63])) >>
                6;
    mvMag8x8 = SQR(avgMvx8x8) + SQR(avgMvy8x8);

    avgSad32x32 =
        (context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x32_0] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x32_1] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x32_2] +
         context_ptr->p_sb_best_sad[listIndex][refPicIndex]
                                   [ME_TIER_ZERO_PU_32x32_3]) >>
        2;
    avgSad16x16 =
        (context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_0] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_1] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_2] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_3] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_4] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_5] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_6] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_7] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_8] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_9] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_10] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_11] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_12] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_13] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_14] +
         context_ptr->p_sb_best_sad[listIndex][refPicIndex]
                                   [ME_TIER_ZERO_PU_16x16_15]) >>
        4;
    avgSad8x8 =
        (context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_0] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_1] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_2] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_3] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_4] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_5] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_6] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_7] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_8] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_9] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_10] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_11] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_12] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_13] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_14] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_15] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_16] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_17] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_18] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_19] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_20] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_21] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_22] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_23] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_24] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_25] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_26] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_27] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_28] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_29] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_30] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_31] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_32] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_33] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_34] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_35] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_36] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_37] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_38] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_39] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_40] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_41] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_42] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_43] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_44] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_45] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_46] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_47] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_48] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_49] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_50] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_51] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_52] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_53] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_54] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_55] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_56] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_57] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_58] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_59] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_60] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_61] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_62] +
         context_ptr
             ->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_63]) >>
        6;

    if (picture_control_set_ptr->temporal_layer_index == 0) {
        // 32x32
        if ((mvMag32x32 < SQR(48)) && (avgSad32x32 < 32 * 32 * 6))
            *enableHalfPel32x32 = EB_TRUE;  // CLASS_0
        else if ((mvMag32x32 < SQR(48)) && !(avgSad32x32 < 32 * 32 * 6))
            *enableHalfPel32x32 = EB_FALSE;  // CLASS_1
        else if (!(mvMag32x32 < SQR(48)) && (avgSad32x32 < 32 * 32 * 6))
            *enableHalfPel32x32 = EB_TRUE;  // CLASS_2
        else
            *enableHalfPel32x32 = EB_FALSE;  // CLASS_3
        // 16x16
        if ((mvMag16x16 < SQR(48)) && (avgSad16x16 < 16 * 16 * 2))
            *enableHalfPel16x16 = EB_FALSE;  // CLASS_0
        else if ((mvMag16x16 < SQR(48)) && !(avgSad16x16 < 16 * 16 * 2))
            *enableHalfPel16x16 = EB_TRUE;  // CLASS_1
        else if (!(mvMag16x16 < SQR(48)) && (avgSad16x16 < 16 * 16 * 2))
            *enableHalfPel16x16 = EB_FALSE;  // CLASS_2
        else
            *enableHalfPel16x16 = EB_TRUE;  // CLASS_3
        // 8x8
        if ((mvMag8x8 < SQR(48)) && (avgSad8x8 < 8 * 8 * 2))
            *enableHalfPel8x8 = EB_FALSE;  // CLASS_0
        else if ((mvMag8x8 < SQR(48)) && !(avgSad8x8 < 8 * 8 * 2))
            *enableHalfPel8x8 = EB_TRUE;  // CLASS_1
        else if (!(mvMag8x8 < SQR(48)) && (avgSad8x8 < 8 * 8 * 2))
            *enableHalfPel8x8 = EB_FALSE;  // CLASS_2
        else
            *enableHalfPel8x8 = EB_TRUE;  // CLASS_3
    }

    else if (picture_control_set_ptr->temporal_layer_index == 1) {
        // 32x32
        if ((mvMag32x32 < SQR(32)) && (avgSad32x32 < 32 * 32 * 6))
            *enableHalfPel32x32 = EB_TRUE;  // CLASS_0
        else if ((mvMag32x32 < SQR(32)) && !(avgSad32x32 < 32 * 32 * 6))
            *enableHalfPel32x32 = EB_FALSE;  // CLASS_1
        else if (!(mvMag32x32 < SQR(32)) && (avgSad32x32 < 32 * 32 * 6))
            *enableHalfPel32x32 = EB_TRUE;  // CLASS_2
        else
            *enableHalfPel32x32 = EB_TRUE;  // CLASS_3
        // 16x16
        if ((mvMag16x16 < SQR(32)) && (avgSad16x16 < 16 * 16 * 2))
            *enableHalfPel16x16 = EB_FALSE;  // CLASS_0
        else if ((mvMag16x16 < SQR(32)) && !(avgSad16x16 < 16 * 16 * 2))
            *enableHalfPel16x16 = EB_TRUE;  // CLASS_1
        else if (!(mvMag16x16 < SQR(32)) && (avgSad16x16 < 16 * 16 * 2))
            *enableHalfPel16x16 = EB_FALSE;  // CLASS_2
        else
            *enableHalfPel16x16 = EB_TRUE;  // CLASS_3
        // 8x8
        if ((mvMag8x8 < SQR(32)) && (avgSad8x8 < 8 * 8 * 2))
            *enableHalfPel8x8 = EB_FALSE;  // CLASS_0
        else if ((mvMag8x8 < SQR(32)) && !(avgSad8x8 < 8 * 8 * 2))
            *enableHalfPel8x8 = EB_TRUE;  // CLASS_1
        else if (!(mvMag8x8 < SQR(32)) && (avgSad8x8 < 8 * 8 * 2))
            *enableHalfPel8x8 = EB_FALSE;  // CLASS_2
        else
            *enableHalfPel8x8 = EB_TRUE;  // CLASS_3
    } else if (picture_control_set_ptr->temporal_layer_index == 2) {
        // 32x32
        if ((mvMag32x32 < SQR(80)) && (avgSad32x32 < 32 * 32 * 6))
            *enableHalfPel32x32 = EB_TRUE;  // CLASS_0
        else if ((mvMag32x32 < SQR(80)) && !(avgSad32x32 < 32 * 32 * 6))
            *enableHalfPel32x32 = EB_FALSE;  // CLASS_1
        else if (!(mvMag32x32 < SQR(80)) && (avgSad32x32 < 32 * 32 * 6))
            *enableHalfPel32x32 = EB_TRUE;  // CLASS_2
        else
            *enableHalfPel32x32 = EB_FALSE;  // CLASS_3
        // 16x16
        if ((mvMag16x16 < SQR(80)) && (avgSad16x16 < 16 * 16 * 2))
            *enableHalfPel16x16 = EB_FALSE;  // CLASS_0
        else if ((mvMag16x16 < SQR(80)) && !(avgSad16x16 < 16 * 16 * 2))
            *enableHalfPel16x16 = EB_TRUE;  // CLASS_1
        else if (!(mvMag16x16 < SQR(80)) && (avgSad16x16 < 16 * 16 * 2))
            *enableHalfPel16x16 = EB_FALSE;  // CLASS_2
        else
            *enableHalfPel16x16 = EB_TRUE;  // CLASS_3
        // 8x8
        if ((mvMag8x8 < SQR(80)) && (avgSad8x8 < 8 * 8 * 2))
            *enableHalfPel8x8 = EB_FALSE;  // CLASS_0
        else if ((mvMag8x8 < SQR(80)) && !(avgSad8x8 < 8 * 8 * 2))
            *enableHalfPel8x8 = EB_TRUE;  // CLASS_1
        else if (!(mvMag8x8 < SQR(80)) && (avgSad8x8 < 8 * 8 * 2))
            *enableHalfPel8x8 = EB_FALSE;  // CLASS_2
        else
            *enableHalfPel8x8 = EB_TRUE;  // CLASS_3
    } else {
        // 32x32
        if ((mvMag32x32 < SQR(48)) && (avgSad32x32 < 32 * 32 * 6))
            *enableHalfPel32x32 = EB_TRUE;  // CLASS_0
        else if ((mvMag32x32 < SQR(48)) && !(avgSad32x32 < 32 * 32 * 6))
            *enableHalfPel32x32 = EB_TRUE;  // CLASS_1
        else if (!(mvMag32x32 < SQR(48)) && (avgSad32x32 < 32 * 32 * 6))
            *enableHalfPel32x32 = EB_TRUE;  // CLASS_2
        else
            *enableHalfPel32x32 = EB_FALSE;  // CLASS_3
        // 16x16
        if ((mvMag16x16 < SQR(48)) && (avgSad16x16 < 16 * 16 * 2))
            *enableHalfPel16x16 = EB_FALSE;  // CLASS_0
        else if ((mvMag16x16 < SQR(48)) && !(avgSad16x16 < 16 * 16 * 2))
            *enableHalfPel16x16 = EB_TRUE;  // CLASS_1
        else if (!(mvMag16x16 < SQR(48)) && (avgSad16x16 < 16 * 16 * 2))
            *enableHalfPel16x16 = EB_FALSE;  // CLASS_2
        else
            *enableHalfPel16x16 = EB_TRUE;  // CLASS_3
        // 8x8
        if ((mvMag8x8 < SQR(48)) && (avgSad8x8 < 8 * 8 * 2))
            *enableHalfPel8x8 = EB_FALSE;  // CLASS_0
        else if ((mvMag8x8 < SQR(48)) && !(avgSad8x8 < 8 * 8 * 2))
            *enableHalfPel8x8 = EB_TRUE;  // CLASS_1
        else if (!(mvMag8x8 < SQR(48)) && (avgSad8x8 < 8 * 8 * 2))
            *enableHalfPel8x8 = EB_FALSE;  // CLASS_2
        else
            *enableHalfPel8x8 = EB_FALSE;  // EB_TRUE; //CLASS_3
    }

    return return_error;
}

static void hme_mv_center_check(EbPictureBufferDesc *ref_pic_ptr,
                                MeContext *context_ptr, int16_t *xsc,
                                int16_t *ysc, uint32_t list_index,
                                int16_t origin_x, int16_t origin_y,
                                uint32_t sb_width, uint32_t sb_height)
{
    // Search for (-srx/2, 0),  (+srx/2, 0),  (0, -sry/2), (0, +sry/2),
    /*
    |------------C-------------|
    |--------------------------|
    |--------------------------|
    A            0             B
    |--------------------------|
    |--------------------------|
    |------------D-------------|
    */
    uint32_t search_region_index;
    int16_t search_center_x = *xsc;
    int16_t search_center_y = *ysc;
    uint64_t best_cost;
    uint64_t direct_mv_cost = 0xFFFFFFFFFFFFF;
    uint8_t sparce_scale = 1;
    int16_t pad_width = (int16_t)BLOCK_SIZE_64 - 1;
    int16_t pad_height = (int16_t)BLOCK_SIZE_64 - 1;
    // O pos

    search_region_index =
        (int16_t)ref_pic_ptr->origin_x + origin_x +
        ((int16_t)ref_pic_ptr->origin_y + origin_y) * ref_pic_ptr->stride_y;

    uint32_t sub_sampled_sad = 1;
    uint64_t zero_mv_sad = nxm_sad_kernel(
            context_ptr->sb_src_ptr,
            context_ptr->sb_src_stride << sub_sampled_sad,
            &(ref_pic_ptr->buffer_y[search_region_index]),
            ref_pic_ptr->stride_y << sub_sampled_sad,
            sb_height >> sub_sampled_sad,
            sb_width);

    zero_mv_sad = zero_mv_sad << sub_sampled_sad;

    uint64_t zero_mv_cost = zero_mv_sad << COST_PRECISION;

    // A pos
    search_center_x =
        0 - (context_ptr->hme_level0_total_search_area_width * sparce_scale);
    search_center_y = 0;

    // Correct the left edge of the Search Area if it is not on the reference
    // Picture
    search_center_x = ((origin_x + search_center_x) < -pad_width)
                          ? -pad_width - origin_x
                          : search_center_x;
    // Correct the right edge of the Search Area if its not on the reference
    // Picture
    search_center_x =
        ((origin_x + search_center_x) > (int16_t)ref_pic_ptr->width - 1)
            ? search_center_x - ((origin_x + search_center_x) -
                                 ((int16_t)ref_pic_ptr->width - 1))
            : search_center_x;
    // Correct the top edge of the Search Area if it is not on the reference
    // Picture
    search_center_y = ((origin_y + search_center_y) < -pad_height)
                          ? -pad_height - origin_y
                          : search_center_y;

    // Correct the bottom edge of the Search Area if its not on the reference
    // Picture
    search_center_y =
        ((origin_y + search_center_y) > (int16_t)ref_pic_ptr->height - 1)
            ? search_center_y - ((origin_y + search_center_y) -
                                 ((int16_t)ref_pic_ptr->height - 1))
            : search_center_y;

    uint64_t mv_a_sad = nxm_sad_kernel(
        context_ptr->sb_src_ptr,
        context_ptr->sb_src_stride << sub_sampled_sad,
        &(ref_pic_ptr->buffer_y[search_region_index]),
        ref_pic_ptr->stride_y << sub_sampled_sad,
        sb_height >> sub_sampled_sad,
        sb_width);

    mv_a_sad = mv_a_sad << sub_sampled_sad;

    uint64_t mv_a_cost = mv_a_sad << COST_PRECISION;

    // B pos
    search_center_x =
        (context_ptr->hme_level0_total_search_area_width * sparce_scale);
    search_center_y = 0;
    ///////////////// correct
    // Correct the left edge of the Search Area if it is not on the reference
    // Picture
    search_center_x = ((origin_x + search_center_x) < -pad_width)
                          ? -pad_width - origin_x
                          : search_center_x;
    // Correct the right edge of the Search Area if its not on the reference
    // Picture
    search_center_x =
        ((origin_x + search_center_x) > (int16_t)ref_pic_ptr->width - 1)
            ? search_center_x - ((origin_x + search_center_x) -
                                 ((int16_t)ref_pic_ptr->width - 1))
            : search_center_x;
    // Correct the top edge of the Search Area if it is not on the reference
    // Picture
    search_center_y = ((origin_y + search_center_y) < -pad_height)
                          ? -pad_height - origin_y
                          : search_center_y;
    // Correct the bottom edge of the Search Area if its not on the reference
    // Picture
    search_center_y =
        ((origin_y + search_center_y) > (int16_t)ref_pic_ptr->height - 1)
            ? search_center_y - ((origin_y + search_center_y) -
                                 ((int16_t)ref_pic_ptr->height - 1))
            : search_center_y;

    search_region_index =
        (int16_t)(ref_pic_ptr->origin_x + origin_x) + search_center_x +
        ((int16_t)(ref_pic_ptr->origin_y + origin_y) + search_center_y) *
            ref_pic_ptr->stride_y;

    uint64_t mv_b_sad = nxm_sad_kernel(
        context_ptr->sb_src_ptr,
        context_ptr->sb_src_stride << sub_sampled_sad,
        &(ref_pic_ptr->buffer_y[search_region_index]),
        ref_pic_ptr->stride_y << sub_sampled_sad,
        sb_height >> sub_sampled_sad,
        sb_width);

    mv_b_sad = mv_b_sad << sub_sampled_sad;

    uint64_t mv_b_cost = mv_b_sad << COST_PRECISION;
    // C pos
    search_center_x = 0;
    search_center_y =
        0 - (context_ptr->hme_level0_total_search_area_height * sparce_scale);
    ///////////////// correct
    // Correct the left edge of the Search Area if it is not on the reference
    // Picture
    search_center_x = ((origin_x + search_center_x) < -pad_width)
                          ? -pad_width - origin_x
                          : search_center_x;

    // Correct the right edge of the Search Area if its not on the reference
    // Picture
    search_center_x =
        ((origin_x + search_center_x) > (int16_t)ref_pic_ptr->width - 1)
            ? search_center_x - ((origin_x + search_center_x) -
                                 ((int16_t)ref_pic_ptr->width - 1))
            : search_center_x;

    // Correct the top edge of the Search Area if it is not on the reference
    // Picture
    search_center_y = ((origin_y + search_center_y) < -pad_height)
                          ? -pad_height - origin_y
                          : search_center_y;

    // Correct the bottom edge of the Search Area if its not on the reference
    // Picture
    search_center_y =
        ((origin_y + search_center_y) > (int16_t)ref_pic_ptr->height - 1)
            ? search_center_y - ((origin_y + search_center_y) -
                                 ((int16_t)ref_pic_ptr->height - 1))
            : search_center_y;

    search_region_index =
        (int16_t)(ref_pic_ptr->origin_x + origin_x) + search_center_x +
        ((int16_t)(ref_pic_ptr->origin_y + origin_y) + search_center_y) *
            ref_pic_ptr->stride_y;

    uint64_t mv_c_sad = nxm_sad_kernel(
        context_ptr->sb_src_ptr,
        context_ptr->sb_src_stride << sub_sampled_sad,
        &(ref_pic_ptr->buffer_y[search_region_index]),
        ref_pic_ptr->stride_y << sub_sampled_sad,
        sb_height >> sub_sampled_sad,
        sb_width);

    mv_c_sad = mv_c_sad << sub_sampled_sad;

    uint64_t mv_c_cost = mv_c_sad << COST_PRECISION;

    // D pos
    search_center_x = 0;
    search_center_y =
        (context_ptr->hme_level0_total_search_area_height * sparce_scale);
    // Correct the left edge of the Search Area if it is not on the reference
    // Picture
    search_center_x = ((origin_x + search_center_x) < -pad_width)
                          ? -pad_width - origin_x
                          : search_center_x;
    // Correct the right edge of the Search Area if its not on the reference
    // Picture
    search_center_x =
        ((origin_x + search_center_x) > (int16_t)ref_pic_ptr->width - 1)
            ? search_center_x - ((origin_x + search_center_x) -
                                 ((int16_t)ref_pic_ptr->width - 1))
            : search_center_x;
    // Correct the top edge of the Search Area if it is not on the reference
    // Picture
    search_center_y = ((origin_y + search_center_y) < -pad_height)
                          ? -pad_height - origin_y
                          : search_center_y;
    // Correct the bottom edge of the Search Area if its not on the reference
    // Picture
    search_center_y =
        ((origin_y + search_center_y) > (int16_t)ref_pic_ptr->height - 1)
            ? search_center_y - ((origin_y + search_center_y) -
                                 ((int16_t)ref_pic_ptr->height - 1))
            : search_center_y;
    search_region_index =
        (int16_t)(ref_pic_ptr->origin_x + origin_x) + search_center_x +
        ((int16_t)(ref_pic_ptr->origin_y + origin_y) + search_center_y) *
            ref_pic_ptr->stride_y;
    uint64_t mv_d_sad = nxm_sad_kernel(
        context_ptr->sb_src_ptr,
        context_ptr->sb_src_stride << sub_sampled_sad,
        &(ref_pic_ptr->buffer_y[search_region_index]),
        ref_pic_ptr->stride_y << sub_sampled_sad,
        sb_height >> sub_sampled_sad,
        sb_width);

    mv_d_sad = mv_d_sad << sub_sampled_sad;

    uint64_t mv_d_cost = mv_d_sad << COST_PRECISION;

    if (list_index == 1) {
        search_center_x =
            list_index ? 0 - (_MVXT(context_ptr->p_sb_best_mv[0][0][0]) >> 2)
                       : 0;
        search_center_y =
            list_index ? 0 - (_MVYT(context_ptr->p_sb_best_mv[0][0][0]) >> 2)
                       : 0;
        ///////////////// correct
        // Correct the left edge of the Search Area if it is not on the
        // reference Picture
        search_center_x = ((origin_x + search_center_x) < -pad_width)
                              ? -pad_width - origin_x
                              : search_center_x;
        // Correct the right edge of the Search Area if its not on the reference
        // Picture
        search_center_x =
            ((origin_x + search_center_x) > (int16_t)ref_pic_ptr->width - 1)
                ? search_center_x - ((origin_x + search_center_x) -
                                     ((int16_t)ref_pic_ptr->width - 1))
                : search_center_x;
        // Correct the top edge of the Search Area if it is not on the reference
        // Picture
        search_center_y = ((origin_y + search_center_y) < -pad_height)
                              ? -pad_height - origin_y
                              : search_center_y;
        // Correct the bottom edge of the Search Area if its not on the
        // reference Picture
        search_center_y =
            ((origin_y + search_center_y) > (int16_t)ref_pic_ptr->height - 1)
                ? search_center_y - ((origin_y + search_center_y) -
                                     ((int16_t)ref_pic_ptr->height - 1))
                : search_center_y;

        search_region_index =
            (int16_t)(ref_pic_ptr->origin_x + origin_x) + search_center_x +
            ((int16_t)(ref_pic_ptr->origin_y + origin_y) + search_center_y) *
                ref_pic_ptr->stride_y;

        uint64_t direct_mv_sad = nxm_sad_kernel(
                context_ptr->sb_src_ptr,
                context_ptr->sb_src_stride << sub_sampled_sad,
                &(ref_pic_ptr->buffer_y[search_region_index]),
                ref_pic_ptr->stride_y << sub_sampled_sad,
                sb_height >> sub_sampled_sad,
                sb_width);

        direct_mv_sad = direct_mv_sad << sub_sampled_sad;

        direct_mv_cost = (direct_mv_sad << COST_PRECISION);
    }

    best_cost = MIN(
        zero_mv_cost,
        MIN(mv_a_cost,
            MIN(mv_b_cost, MIN(mv_c_cost, MIN(mv_d_cost, direct_mv_cost)))));

    if (best_cost == zero_mv_cost) {
        search_center_x = 0;
        search_center_y = 0;
    } else if (best_cost == mv_a_cost) {
        search_center_x = 0 - (context_ptr->hme_level0_total_search_area_width *
                               sparce_scale);
        search_center_y = 0;
    } else if (best_cost == mv_b_cost) {
        search_center_x =
            (context_ptr->hme_level0_total_search_area_width * sparce_scale);
        search_center_y = 0;
    } else if (best_cost == mv_c_cost) {
        search_center_x = 0;
        search_center_y =
            0 -
            (context_ptr->hme_level0_total_search_area_height * sparce_scale);
    } else if (best_cost == direct_mv_cost) {
        search_center_x =
            list_index ? 0 - (_MVXT(context_ptr->p_sb_best_mv[0][0][0]) >> 2)
                       : 0;
        search_center_y =
            list_index ? 0 - (_MVYT(context_ptr->p_sb_best_mv[0][0][0]) >> 2)
                       : 0;
    } else if (best_cost == mv_d_cost) {
        search_center_x = 0;
        search_center_y =
            (context_ptr->hme_level0_total_search_area_height * sparce_scale);
    }

    else
        SVT_LOG("error no center selected");
    *xsc = search_center_x;
    *ysc = search_center_y;
}

void SwapMeCandidate(MePredUnit *a, MePredUnit *b) {
    MePredUnit tempPtr;
    tempPtr = *a;
    *a = *b;
    *b = tempPtr;
}

/*******************************************
 * motion_estimate_lcu
 *   performs ME (LCU)
 *******************************************/
EbErrorType motion_estimate_lcu(
        PictureParentControlSet   *picture_control_set_ptr,  // input parameter, Picture Control Set Ptr
        uint32_t                   sb_index,              // input parameter, SB Index
        uint32_t                   sb_origin_x,            // input parameter, SB Origin X
        uint32_t                   sb_origin_y,            // input parameter, SB Origin X
        MeContext                 *context_ptr,                        // input parameter, ME Context Ptr, used to store decimated/interpolated LCU/SR
        EbPictureBufferDesc       *input_ptr)              // input parameter, source Picture Ptr

{
    EbErrorType return_error = EB_ErrorNone;

    SequenceControlSet *sequence_control_set_ptr =
        (SequenceControlSet *)picture_control_set_ptr
            ->sequence_control_set_wrapper_ptr->object_ptr;

    int16_t xTopLeftSearchRegion;
    int16_t yTopLeftSearchRegion;
    uint32_t searchRegionIndex;

    int16_t picture_width =
        (int16_t)((SequenceControlSet *)picture_control_set_ptr
                      ->sequence_control_set_wrapper_ptr->object_ptr)
            ->seq_header.max_frame_width;
    int16_t picture_height =
        (int16_t)((SequenceControlSet *)picture_control_set_ptr
                      ->sequence_control_set_wrapper_ptr->object_ptr)
            ->seq_header.max_frame_height;
    uint32_t sb_width = (input_ptr->width - sb_origin_x) < BLOCK_SIZE_64
                            ? input_ptr->width - sb_origin_x
                            : BLOCK_SIZE_64;
    uint32_t sb_height = (input_ptr->height - sb_origin_y) < BLOCK_SIZE_64
                             ? input_ptr->height - sb_origin_y
                             : BLOCK_SIZE_64;

    int16_t padWidth = (int16_t)BLOCK_SIZE_64 - 1;
    int16_t padHeight = (int16_t)BLOCK_SIZE_64 - 1;
    int16_t search_area_width;
    int16_t search_area_height;
    int16_t x_search_area_origin;
    int16_t y_search_area_origin;
    int16_t origin_x = (int16_t)sb_origin_x;
    int16_t origin_y = (int16_t)sb_origin_y;

    // HME
    uint32_t searchRegionNumberInWidth = 0;
    uint32_t searchRegionNumberInHeight = 0;
    int16_t xHmeLevel0SearchCenter[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT]
                                  [EB_HME_SEARCH_AREA_ROW_MAX_COUNT];
    int16_t yHmeLevel0SearchCenter[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT]
                                  [EB_HME_SEARCH_AREA_ROW_MAX_COUNT];
    uint64_t hmeLevel0Sad[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT]
                         [EB_HME_SEARCH_AREA_ROW_MAX_COUNT];
    int16_t xHmeLevel1SearchCenter[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT]
                                  [EB_HME_SEARCH_AREA_ROW_MAX_COUNT];
    int16_t yHmeLevel1SearchCenter[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT]
                                  [EB_HME_SEARCH_AREA_ROW_MAX_COUNT];
    uint64_t hmeLevel1Sad[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT]
                         [EB_HME_SEARCH_AREA_ROW_MAX_COUNT];
    int16_t xHmeLevel2SearchCenter[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT]
                                  [EB_HME_SEARCH_AREA_ROW_MAX_COUNT];
    int16_t yHmeLevel2SearchCenter[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT]
                                  [EB_HME_SEARCH_AREA_ROW_MAX_COUNT];
    uint64_t hmeLevel2Sad[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT]
                         [EB_HME_SEARCH_AREA_ROW_MAX_COUNT];

    // Hierarchical ME Search Center
    int16_t xHmeSearchCenter = 0;
    int16_t yHmeSearchCenter = 0;

    // Final ME Search Center
    int16_t x_search_center = 0;
    int16_t y_search_center = 0;

    // Search Center SADs
    uint64_t hmeMvSad = 0;

    uint32_t pu_index;

    uint32_t max_number_of_pus_per_sb =
        picture_control_set_ptr->max_number_of_pus_per_sb;

    uint32_t numOfListToSearch;
    uint32_t listIndex;
    uint8_t candidateIndex = 0;
    uint8_t total_me_candidate_index = 0;
    EbPaReferenceObject
        *referenceObject;  // input parameter, reference Object Ptr

    uint8_t ref_pic_index;
    uint8_t num_of_ref_pic_to_search;
    uint8_t candidate_index = 0;
    uint32_t next_candidate_index = 0;

    MePredUnit *me_candidate;
    EbPictureBufferDesc *refPicPtr;
    EbPictureBufferDesc *quarterRefPicPtr;
    EbPictureBufferDesc *sixteenthRefPicPtr;

    int16_t tempXHmeSearchCenter = 0;
    int16_t tempYHmeSearchCenter = 0;

    uint32_t numQuadInWidth;
    uint32_t totalMeQuad;
    uint32_t quadIndex;
    uint32_t nextQuadIndex;
    uint64_t tempXHmeSad;

    uint64_t ref0Poc = 0;
    uint64_t ref1Poc = 0;

    EbAsm asm_type = sequence_control_set_ptr->encode_context_ptr->asm_type;
    uint64_t i;

    int16_t hmeLevel1SearchAreaInWidth;
    int16_t hmeLevel1SearchAreaInHeight;
    // Configure HME level 0, level 1 and level 2 from static config parameters
    EbBool enable_hme_level0_flag =
        context_ptr->enable_hme_level0_flag;
    EbBool enable_hme_level1_flag =
        context_ptr->enable_hme_level1_flag;
    EbBool enable_hme_level2_flag =
        context_ptr->enable_hme_level2_flag;

    EbBool enableHalfPel32x32 = EB_FALSE;
    EbBool enableHalfPel16x16 = EB_FALSE;
    EbBool enableHalfPel8x8 = EB_FALSE;
    EbBool enableQuarterPel = EB_FALSE;
    EbBool oneQuadrantHME = EB_FALSE;

    oneQuadrantHME =
        sequence_control_set_ptr->input_resolution < INPUT_SIZE_4K_RANGE
            ? 0
            : oneQuadrantHME;

    numOfListToSearch = (picture_control_set_ptr->slice_type == P_SLICE)
                            ? (uint32_t)REF_LIST_0
                            : (uint32_t)REF_LIST_1;

    EbBool is_nsq_table_used =
        (picture_control_set_ptr->pic_depth_mode <= PIC_ALL_C_DEPTH_MODE &&
         picture_control_set_ptr->nsq_search_level >= NSQ_SEARCH_LEVEL1 &&
         picture_control_set_ptr->nsq_search_level < NSQ_SEARCH_FULL)
            ? EB_TRUE
            : EB_FALSE;

    is_nsq_table_used = picture_control_set_ptr->enc_mode == ENC_M0 ?  EB_FALSE : is_nsq_table_used;

    if (context_ptr->me_alt_ref == EB_TRUE)
        numOfListToSearch = 0;

    // Uni-Prediction motion estimation loop
    // List Loop
    for (listIndex = REF_LIST_0; listIndex <= numOfListToSearch; ++listIndex) {

        if (context_ptr->me_alt_ref == EB_TRUE) {
            num_of_ref_pic_to_search = 1;
        } else {
            num_of_ref_pic_to_search =
                (picture_control_set_ptr->slice_type == P_SLICE)
                    ? picture_control_set_ptr->ref_list0_count
                    : (listIndex == REF_LIST_0)
                          ? picture_control_set_ptr->ref_list0_count
                          : picture_control_set_ptr->ref_list1_count;

            referenceObject = (EbPaReferenceObject *)picture_control_set_ptr
                                  ->ref_pa_pic_ptr_array[0][0]
                                  ->object_ptr;
            ref0Poc = picture_control_set_ptr->ref_pic_poc_array[0][0];
        }

        // Ref Picture Loop
        for (ref_pic_index = 0; ref_pic_index < num_of_ref_pic_to_search;
             ++ref_pic_index)
        {
            if (context_ptr->me_alt_ref == EB_TRUE) {
                referenceObject =
                    (EbPaReferenceObject *)context_ptr->alt_ref_reference_ptr;
            } else {
                if (numOfListToSearch) {
                    referenceObject =
                        (EbPaReferenceObject *)picture_control_set_ptr
                            ->ref_pa_pic_ptr_array[1][0]
                            ->object_ptr;
                    ref1Poc = picture_control_set_ptr->ref_pic_poc_array[1][0];
                }

                referenceObject =
                    (EbPaReferenceObject *)picture_control_set_ptr
                        ->ref_pa_pic_ptr_array[listIndex][ref_pic_index]
                        ->object_ptr;
            }

            refPicPtr = (EbPictureBufferDesc*)referenceObject->input_padded_picture_ptr;
            // Set 1/4 and 1/16 ME reference buffer(s); filtered or decimated
            quarterRefPicPtr = (sequence_control_set_ptr->down_sampling_method_me_search == ME_FILTERED_DOWNSAMPLED) ?
                (EbPictureBufferDesc*)referenceObject->quarter_filtered_picture_ptr :
                (EbPictureBufferDesc*)referenceObject->quarter_decimated_picture_ptr;

            sixteenthRefPicPtr = (sequence_control_set_ptr->down_sampling_method_me_search == ME_FILTERED_DOWNSAMPLED) ?
                (EbPictureBufferDesc*)referenceObject->sixteenth_filtered_picture_ptr:
                (EbPictureBufferDesc*)referenceObject->sixteenth_decimated_picture_ptr;
            if (picture_control_set_ptr->temporal_layer_index > 0 ||
                listIndex == 0) {
                // A - The MV center for Tier0 search could be either (0,0), or
                // HME A - Set HME MV Center
                if (context_ptr->update_hme_search_center_flag)
                    hme_mv_center_check(refPicPtr,
                                        context_ptr,
                                        &x_search_center,
                                        &y_search_center,
                                        listIndex,
                                        origin_x,
                                        origin_y,
                                        sb_width,
                                        sb_height);
                else {
                    x_search_center = 0;
                    y_search_center = 0;
                }
                // B - NO HME in boundaries
                // C - Skip HME

                if (context_ptr->enable_hme_flag &&

                    /*B*/ sb_height ==
                        BLOCK_SIZE_64) {  //(searchCenterSad >
                                          // sequence_control_set_ptr->static_config.skipTier0HmeTh))
                                          //{
                    while (searchRegionNumberInHeight <
                           context_ptr->number_hme_search_region_in_height) {
                        while (searchRegionNumberInWidth <
                               context_ptr->number_hme_search_region_in_width) {
                            xHmeLevel0SearchCenter[searchRegionNumberInWidth]
                                                  [searchRegionNumberInHeight] =
                                                      x_search_center;
                            yHmeLevel0SearchCenter[searchRegionNumberInWidth]
                                                  [searchRegionNumberInHeight] =
                                                      y_search_center;

                            xHmeLevel1SearchCenter[searchRegionNumberInWidth]
                                                  [searchRegionNumberInHeight] =
                                                      x_search_center;
                            yHmeLevel1SearchCenter[searchRegionNumberInWidth]
                                                  [searchRegionNumberInHeight] =
                                                      y_search_center;

                            xHmeLevel2SearchCenter[searchRegionNumberInWidth]
                                                  [searchRegionNumberInHeight] =
                                                      x_search_center;
                            yHmeLevel2SearchCenter[searchRegionNumberInWidth]
                                                  [searchRegionNumberInHeight] =
                                                      y_search_center;

                            searchRegionNumberInWidth++;
                        }
                        searchRegionNumberInWidth = 0;
                        searchRegionNumberInHeight++;
                    }

                    // HME: Level0 search

                    if (enable_hme_level0_flag) {
                        if (oneQuadrantHME && !enable_hme_level1_flag &&
                            !enable_hme_level2_flag) {
                            searchRegionNumberInHeight = 0;
                            searchRegionNumberInWidth = 0;

                            HmeOneQuadrantLevel0(
                                picture_control_set_ptr,
                                context_ptr,
                                origin_x >> 2,
                                origin_y >> 2,
                                sb_width >> 2,
                                sb_height >> 2,
                                x_search_center >> 2,
                                y_search_center >> 2,
                                sixteenthRefPicPtr,
                                &(hmeLevel0Sad[searchRegionNumberInWidth]
                                              [searchRegionNumberInHeight]),
                                &(xHmeLevel0SearchCenter
                                      [searchRegionNumberInWidth]
                                      [searchRegionNumberInHeight]),
                                &(yHmeLevel0SearchCenter
                                      [searchRegionNumberInWidth]
                                      [searchRegionNumberInHeight]),
                                hme_level_0_search_area_multiplier_x
                                    [picture_control_set_ptr
                                         ->hierarchical_levels]
                                    [picture_control_set_ptr
                                         ->temporal_layer_index],
                                hme_level_0_search_area_multiplier_y
                                    [picture_control_set_ptr
                                         ->hierarchical_levels]
                                    [picture_control_set_ptr
                                         ->temporal_layer_index],
                                asm_type);
                        } else {
                            searchRegionNumberInHeight = 0;
                            searchRegionNumberInWidth = 0;
                            {
                                while (
                                    searchRegionNumberInHeight <
                                    context_ptr
                                        ->number_hme_search_region_in_height) {
                                    while (
                                        searchRegionNumberInWidth <
                                        context_ptr
                                            ->number_hme_search_region_in_width) {
                                        HmeLevel0(
                                            picture_control_set_ptr,
                                            context_ptr,
                                            origin_x >> 2,
                                            origin_y >> 2,
                                            sb_width >> 2,
                                            sb_height >> 2,
                                            x_search_center >> 2,
                                            y_search_center >> 2,
                                            sixteenthRefPicPtr,
                                            searchRegionNumberInWidth,
                                            searchRegionNumberInHeight,
                                            &(hmeLevel0Sad
                                                  [searchRegionNumberInWidth]
                                                  [searchRegionNumberInHeight]),
                                            &(xHmeLevel0SearchCenter
                                                  [searchRegionNumberInWidth]
                                                  [searchRegionNumberInHeight]),
                                            &(yHmeLevel0SearchCenter
                                                  [searchRegionNumberInWidth]
                                                  [searchRegionNumberInHeight]),
                                            hme_level_0_search_area_multiplier_x
                                                [picture_control_set_ptr
                                                     ->hierarchical_levels]
                                                [picture_control_set_ptr
                                                     ->temporal_layer_index],
                                            hme_level_0_search_area_multiplier_y
                                                [picture_control_set_ptr
                                                     ->hierarchical_levels]
                                                [picture_control_set_ptr
                                                     ->temporal_layer_index],
                                            asm_type);

                                        searchRegionNumberInWidth++;
                                    }
                                    searchRegionNumberInWidth = 0;
                                    searchRegionNumberInHeight++;
                                }
                            }
                        }
                    }

                    // HME: Level1 search
                    if (enable_hme_level1_flag) {
                        searchRegionNumberInHeight = 0;
                        searchRegionNumberInWidth = 0;

                        {
                            while (searchRegionNumberInHeight <
                                   context_ptr
                                       ->number_hme_search_region_in_height) {
                                while (
                                    searchRegionNumberInWidth <
                                    context_ptr
                                        ->number_hme_search_region_in_width) {
                                    // When HME level 0 has been disabled,
                                    // increase the search area width and height
                                    // for level 1 to (32x12) for Gold only

                                    hmeLevel1SearchAreaInWidth =
                                        (int16_t)context_ptr
                                            ->hme_level1_search_area_in_width_array
                                                [searchRegionNumberInWidth];
                                    hmeLevel1SearchAreaInHeight =
                                        (int16_t)context_ptr
                                            ->hme_level1_search_area_in_height_array
                                                [searchRegionNumberInHeight];

                                    HmeLevel1(
                                        context_ptr,
                                        origin_x >> 1,
                                        origin_y >> 1,
                                        sb_width >> 1,
                                        sb_height >> 1,
                                        quarterRefPicPtr,
                                        hmeLevel1SearchAreaInWidth,
                                        hmeLevel1SearchAreaInHeight,
                                        xHmeLevel0SearchCenter
                                                [searchRegionNumberInWidth]
                                                [searchRegionNumberInHeight] >>
                                            1,
                                        yHmeLevel0SearchCenter
                                                [searchRegionNumberInWidth]
                                                [searchRegionNumberInHeight] >>
                                            1,
                                        &(hmeLevel1Sad
                                              [searchRegionNumberInWidth]
                                              [searchRegionNumberInHeight]),
                                        &(xHmeLevel1SearchCenter
                                              [searchRegionNumberInWidth]
                                              [searchRegionNumberInHeight]),
                                        &(yHmeLevel1SearchCenter
                                              [searchRegionNumberInWidth]
                                              [searchRegionNumberInHeight]));

                                    searchRegionNumberInWidth++;
                                }
                                searchRegionNumberInWidth = 0;
                                searchRegionNumberInHeight++;
                            }
                        }
                    }

                    // HME: Level2 search
                    if (enable_hme_level2_flag) {
                        searchRegionNumberInHeight = 0;
                        searchRegionNumberInWidth = 0;

                        {
                            while (searchRegionNumberInHeight <
                                   context_ptr
                                       ->number_hme_search_region_in_height) {
                                while (
                                    searchRegionNumberInWidth <
                                    context_ptr
                                        ->number_hme_search_region_in_width) {
                                    HmeLevel2(
                                        picture_control_set_ptr,
                                        context_ptr,
                                        origin_x,
                                        origin_y,
                                        sb_width,
                                        sb_height,
                                        refPicPtr,
                                        searchRegionNumberInWidth,
                                        searchRegionNumberInHeight,
                                        xHmeLevel1SearchCenter
                                            [searchRegionNumberInWidth]
                                            [searchRegionNumberInHeight],
                                        yHmeLevel1SearchCenter
                                            [searchRegionNumberInWidth]
                                            [searchRegionNumberInHeight],
                                        &(hmeLevel2Sad
                                              [searchRegionNumberInWidth]
                                              [searchRegionNumberInHeight]),
                                        &(xHmeLevel2SearchCenter
                                              [searchRegionNumberInWidth]
                                              [searchRegionNumberInHeight]),
                                        &(yHmeLevel2SearchCenter
                                              [searchRegionNumberInWidth]
                                              [searchRegionNumberInHeight]));

                                    searchRegionNumberInWidth++;
                                }
                                searchRegionNumberInWidth = 0;
                                searchRegionNumberInHeight++;
                            }
                        }
                    }

                    // Hierarchical ME - Search Center
                    if (enable_hme_level0_flag && !enable_hme_level1_flag &&
                        !enable_hme_level2_flag) {
                        if (oneQuadrantHME) {
                            xHmeSearchCenter = xHmeLevel0SearchCenter[0][0];
                            yHmeSearchCenter = yHmeLevel0SearchCenter[0][0];
                            hmeMvSad = hmeLevel0Sad[0][0];
                        } else {
                            xHmeSearchCenter = xHmeLevel0SearchCenter[0][0];
                            yHmeSearchCenter = yHmeLevel0SearchCenter[0][0];
                            hmeMvSad = hmeLevel0Sad[0][0];

                            searchRegionNumberInWidth = 1;
                            searchRegionNumberInHeight = 0;

                            while (searchRegionNumberInHeight <
                                   context_ptr
                                       ->number_hme_search_region_in_height) {
                                while (
                                    searchRegionNumberInWidth <
                                    context_ptr
                                        ->number_hme_search_region_in_width) {
                                    xHmeSearchCenter =
                                        (hmeLevel0Sad
                                             [searchRegionNumberInWidth]
                                             [searchRegionNumberInHeight] <
                                         hmeMvSad)
                                            ? xHmeLevel0SearchCenter
                                                  [searchRegionNumberInWidth]
                                                  [searchRegionNumberInHeight]
                                            : xHmeSearchCenter;
                                    yHmeSearchCenter =
                                        (hmeLevel0Sad
                                             [searchRegionNumberInWidth]
                                             [searchRegionNumberInHeight] <
                                         hmeMvSad)
                                            ? yHmeLevel0SearchCenter
                                                  [searchRegionNumberInWidth]
                                                  [searchRegionNumberInHeight]
                                            : yHmeSearchCenter;
                                    hmeMvSad =
                                        (hmeLevel0Sad
                                             [searchRegionNumberInWidth]
                                             [searchRegionNumberInHeight] <
                                         hmeMvSad)
                                            ? hmeLevel0Sad
                                                  [searchRegionNumberInWidth]
                                                  [searchRegionNumberInHeight]
                                            : hmeMvSad;
                                    searchRegionNumberInWidth++;
                                }
                                searchRegionNumberInWidth = 0;
                                searchRegionNumberInHeight++;
                            }
                        }
                    }

                    if (enable_hme_level1_flag && !enable_hme_level2_flag) {
                        xHmeSearchCenter = xHmeLevel1SearchCenter[0][0];
                        yHmeSearchCenter = yHmeLevel1SearchCenter[0][0];
                        hmeMvSad = hmeLevel1Sad[0][0];

                        searchRegionNumberInWidth = 1;
                        searchRegionNumberInHeight = 0;

                        while (
                            searchRegionNumberInHeight <
                            context_ptr->number_hme_search_region_in_height) {
                            while (searchRegionNumberInWidth <
                                   context_ptr
                                       ->number_hme_search_region_in_width) {
                                xHmeSearchCenter =
                                    (hmeLevel1Sad[searchRegionNumberInWidth]
                                                 [searchRegionNumberInHeight] <
                                     hmeMvSad)
                                        ? xHmeLevel1SearchCenter
                                              [searchRegionNumberInWidth]
                                              [searchRegionNumberInHeight]
                                        : xHmeSearchCenter;
                                yHmeSearchCenter =
                                    (hmeLevel1Sad[searchRegionNumberInWidth]
                                                 [searchRegionNumberInHeight] <
                                     hmeMvSad)
                                        ? yHmeLevel1SearchCenter
                                              [searchRegionNumberInWidth]
                                              [searchRegionNumberInHeight]
                                        : yHmeSearchCenter;
                                hmeMvSad =
                                    (hmeLevel1Sad[searchRegionNumberInWidth]
                                                 [searchRegionNumberInHeight] <
                                     hmeMvSad)
                                        ? hmeLevel1Sad
                                              [searchRegionNumberInWidth]
                                              [searchRegionNumberInHeight]
                                        : hmeMvSad;
                                searchRegionNumberInWidth++;
                            }
                            searchRegionNumberInWidth = 0;
                            searchRegionNumberInHeight++;
                        }
                    }

                    if (enable_hme_level2_flag) {
                        xHmeSearchCenter = xHmeLevel2SearchCenter[0][0];
                        yHmeSearchCenter = yHmeLevel2SearchCenter[0][0];
                        hmeMvSad = hmeLevel2Sad[0][0];

                        searchRegionNumberInWidth = 1;
                        searchRegionNumberInHeight = 0;

                        while (
                            searchRegionNumberInHeight <
                            context_ptr->number_hme_search_region_in_height) {
                            while (searchRegionNumberInWidth <
                                   context_ptr
                                       ->number_hme_search_region_in_width) {
                                xHmeSearchCenter =
                                    (hmeLevel2Sad[searchRegionNumberInWidth]
                                                 [searchRegionNumberInHeight] <
                                     hmeMvSad)
                                        ? xHmeLevel2SearchCenter
                                              [searchRegionNumberInWidth]
                                              [searchRegionNumberInHeight]
                                        : xHmeSearchCenter;
                                yHmeSearchCenter =
                                    (hmeLevel2Sad[searchRegionNumberInWidth]
                                                 [searchRegionNumberInHeight] <
                                     hmeMvSad)
                                        ? yHmeLevel2SearchCenter
                                              [searchRegionNumberInWidth]
                                              [searchRegionNumberInHeight]
                                        : yHmeSearchCenter;
                                hmeMvSad =
                                    (hmeLevel2Sad[searchRegionNumberInWidth]
                                                 [searchRegionNumberInHeight] <
                                     hmeMvSad)
                                        ? hmeLevel2Sad
                                              [searchRegionNumberInWidth]
                                              [searchRegionNumberInHeight]
                                        : hmeMvSad;
                                searchRegionNumberInWidth++;
                            }
                            searchRegionNumberInWidth = 0;
                            searchRegionNumberInHeight++;
                        }

                        numQuadInWidth =
                            context_ptr->number_hme_search_region_in_width;
                        totalMeQuad =
                            context_ptr->number_hme_search_region_in_height *
                            context_ptr->number_hme_search_region_in_width;

                        if ((ref0Poc == ref1Poc) && (listIndex == 1) &&
                            (totalMeQuad > 1)) {
                            for (quadIndex = 0; quadIndex < totalMeQuad - 1;
                                 ++quadIndex) {
                                for (nextQuadIndex = quadIndex + 1;
                                     nextQuadIndex < totalMeQuad;
                                     ++nextQuadIndex) {
                                    if (hmeLevel2Sad[quadIndex / numQuadInWidth]
                                                    [quadIndex %
                                                     numQuadInWidth] >
                                        hmeLevel2Sad[nextQuadIndex /
                                                     numQuadInWidth]
                                                    [nextQuadIndex %
                                                     numQuadInWidth]) {
                                        tempXHmeSearchCenter =
                                            xHmeLevel2SearchCenter
                                                [quadIndex / numQuadInWidth]
                                                [quadIndex % numQuadInWidth];
                                        tempYHmeSearchCenter =
                                            yHmeLevel2SearchCenter
                                                [quadIndex / numQuadInWidth]
                                                [quadIndex % numQuadInWidth];
                                        tempXHmeSad =
                                            hmeLevel2Sad[quadIndex /
                                                         numQuadInWidth]
                                                        [quadIndex %
                                                         numQuadInWidth];

                                        xHmeLevel2SearchCenter
                                            [quadIndex / numQuadInWidth]
                                            [quadIndex % numQuadInWidth] =
                                                xHmeLevel2SearchCenter
                                                    [nextQuadIndex /
                                                     numQuadInWidth]
                                                    [nextQuadIndex %
                                                     numQuadInWidth];
                                        yHmeLevel2SearchCenter
                                            [quadIndex / numQuadInWidth]
                                            [quadIndex % numQuadInWidth] =
                                                yHmeLevel2SearchCenter
                                                    [nextQuadIndex /
                                                     numQuadInWidth]
                                                    [nextQuadIndex %
                                                     numQuadInWidth];
                                        hmeLevel2Sad
                                            [quadIndex / numQuadInWidth]
                                            [quadIndex % numQuadInWidth] =
                                                hmeLevel2Sad[nextQuadIndex /
                                                             numQuadInWidth]
                                                            [nextQuadIndex %
                                                             numQuadInWidth];

                                        xHmeLevel2SearchCenter
                                            [nextQuadIndex / numQuadInWidth]
                                            [nextQuadIndex % numQuadInWidth] =
                                                tempXHmeSearchCenter;
                                        yHmeLevel2SearchCenter
                                            [nextQuadIndex / numQuadInWidth]
                                            [nextQuadIndex % numQuadInWidth] =
                                                tempYHmeSearchCenter;
                                        hmeLevel2Sad[nextQuadIndex /
                                                     numQuadInWidth]
                                                    [nextQuadIndex %
                                                     numQuadInWidth] =
                                                        tempXHmeSad;
                                    }
                                }
                            }

                            xHmeSearchCenter = xHmeLevel2SearchCenter[0][1];
                            yHmeSearchCenter = yHmeLevel2SearchCenter[0][1];
                        }
                    }

                    x_search_center = xHmeSearchCenter;
                    y_search_center = yHmeSearchCenter;
                }
            }

            else {
                x_search_center = 0;
                y_search_center = 0;
            }
            // Constrain x_ME to be a multiple of 8 (round up)
            search_area_width = (context_ptr->search_area_width + 7) & ~0x07;
            search_area_height = context_ptr->search_area_height;
            if ((x_search_center != 0 || y_search_center != 0) &&
                (picture_control_set_ptr->is_used_as_reference_flag ==
                 EB_TRUE)) {
                CheckZeroZeroCenter(refPicPtr,
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

            if(sequence_control_set_ptr->static_config.unrestricted_motion_vector == 0)
            {
                int tile_start_x = sequence_control_set_ptr->sb_params_array[sb_index].tile_start_x;
                int tile_end_x   = sequence_control_set_ptr->sb_params_array[sb_index].tile_end_x;

                // Correct the left edge of the Search Area if it is not on the
                // reference Picture
                x_search_area_origin =
                    ((origin_x + x_search_area_origin) < tile_start_x)
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
                        ? x_search_area_origin - ((origin_x + x_search_area_origin) - (tile_end_x - 1))
                        : x_search_area_origin;

                search_area_width =
                    ((origin_x + x_search_area_origin + search_area_width) > tile_end_x)
                        ? MAX(1, search_area_width - ((origin_x + x_search_area_origin + search_area_width) - tile_end_x))
                        : search_area_width;

                // Constrain x_ME to be a multiple of 8 (round down as cropping
                // already performed)
                search_area_width = (search_area_width < 8)
                                        ? search_area_width
                                        : search_area_width & ~0x07;
            } else {
            // Correct the left edge of the Search Area if it is not on the
            // reference Picture
            x_search_area_origin =
                ((origin_x + x_search_area_origin) < -padWidth)
                    ? -padWidth - origin_x
                    : x_search_area_origin;

            search_area_width =
                ((origin_x + x_search_area_origin) < -padWidth)
                    ? search_area_width -
                          (-padWidth - (origin_x + x_search_area_origin))
                    : search_area_width;

            // Correct the right edge of the Search Area if its not on the
            // reference Picture
            x_search_area_origin =
                ((origin_x + x_search_area_origin) > picture_width - 1)
                    ? x_search_area_origin -
                          ((origin_x + x_search_area_origin) -
                           (picture_width - 1))
                    : x_search_area_origin;

            search_area_width =
                ((origin_x + x_search_area_origin + search_area_width) >
                 picture_width)
                    ? MAX(1,
                          search_area_width -
                              ((origin_x + x_search_area_origin +
                                search_area_width) -
                               picture_width))
                    : search_area_width;

            // Constrain x_ME to be a multiple of 8 (round down as cropping
            // already performed)
            search_area_width = (search_area_width < 8)
                                    ? search_area_width
                                    : search_area_width & ~0x07;
            }

            if(sequence_control_set_ptr->static_config.unrestricted_motion_vector == 0)
            {
                int tile_start_y = sequence_control_set_ptr->sb_params_array[sb_index].tile_start_y;
                int tile_end_y   = sequence_control_set_ptr->sb_params_array[sb_index].tile_end_y;

                // Correct the top edge of the Search Area if it is not on the
                // reference Picture
                y_search_area_origin =
                    ((origin_y + y_search_area_origin) < tile_start_y)
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
                        ? y_search_area_origin - ((origin_y + y_search_area_origin) - (tile_end_y - 1))
                        : y_search_area_origin;

                search_area_height =
                    (origin_y + y_search_area_origin + search_area_height > tile_end_y)
                        ? MAX(1, search_area_height - ((origin_y + y_search_area_origin + search_area_height) - tile_end_y))
                        : search_area_height;
            } else {
            // Correct the top edge of the Search Area if it is not on the
            // reference Picture
            y_search_area_origin =
                ((origin_y + y_search_area_origin) < -padHeight)
                    ? -padHeight - origin_y
                    : y_search_area_origin;

            search_area_height =
                ((origin_y + y_search_area_origin) < -padHeight)
                    ? search_area_height -
                          (-padHeight - (origin_y + y_search_area_origin))
                    : search_area_height;

            // Correct the bottom edge of the Search Area if its not on the
            // reference Picture
            y_search_area_origin =
                ((origin_y + y_search_area_origin) > picture_height - 1)
                    ? y_search_area_origin -
                          ((origin_y + y_search_area_origin) -
                           (picture_height - 1))
                    : y_search_area_origin;

            search_area_height =
                (origin_y + y_search_area_origin + search_area_height >
                 picture_height)
                    ? MAX(1,
                          search_area_height -
                              ((origin_y + y_search_area_origin +
                                search_area_height) -
                               picture_height))
                    : search_area_height;
            }
            context_ptr->x_search_area_origin[listIndex][ref_pic_index] =
                x_search_area_origin;
            context_ptr->y_search_area_origin[listIndex][ref_pic_index] =
                y_search_area_origin;

            context_ptr->adj_search_area_width = search_area_width;
            context_ptr->adj_search_area_height = search_area_height;

            xTopLeftSearchRegion =
                (int16_t)(refPicPtr->origin_x + sb_origin_x) -
                (ME_FILTER_TAP >> 1) + x_search_area_origin;
            yTopLeftSearchRegion =
                (int16_t)(refPicPtr->origin_y + sb_origin_y) -
                (ME_FILTER_TAP >> 1) + y_search_area_origin;
            searchRegionIndex = (xTopLeftSearchRegion) +
                                (yTopLeftSearchRegion)*refPicPtr->stride_y;
            context_ptr->integer_buffer_ptr[listIndex][ref_pic_index] =
                &(refPicPtr->buffer_y[searchRegionIndex]);
            context_ptr->interpolated_full_stride[listIndex][ref_pic_index] =
                refPicPtr->stride_y;

            // Move to the top left of the search region
            xTopLeftSearchRegion =
                (int16_t)(refPicPtr->origin_x + sb_origin_x) +
                x_search_area_origin;
            yTopLeftSearchRegion =
                (int16_t)(refPicPtr->origin_y + sb_origin_y) +
                y_search_area_origin;
            searchRegionIndex = xTopLeftSearchRegion +
                                yTopLeftSearchRegion * refPicPtr->stride_y;

            {
                {
                    if (picture_control_set_ptr->pic_depth_mode <=
                        PIC_ALL_C_DEPTH_MODE) {
                        initialize_buffer_32bits(
                            context_ptr
                                ->p_sb_best_sad[listIndex][ref_pic_index],
                            52,
                            1,
                            MAX_SAD_VALUE);

                        context_ptr->p_best_sad64x64 = &(
                            context_ptr->p_sb_best_sad[listIndex][ref_pic_index]
                                                      [ME_TIER_ZERO_PU_64x64]);
                        context_ptr->p_best_sad32x32 =
                            &(context_ptr
                                  ->p_sb_best_sad[listIndex][ref_pic_index]
                                                 [ME_TIER_ZERO_PU_32x32_0]);
                        context_ptr->p_best_sad16x16 =
                            &(context_ptr
                                  ->p_sb_best_sad[listIndex][ref_pic_index]
                                                 [ME_TIER_ZERO_PU_16x16_0]);
                        context_ptr->p_best_sad8x8 = &(
                            context_ptr->p_sb_best_sad[listIndex][ref_pic_index]
                                                      [ME_TIER_ZERO_PU_8x8_0]);
                        context_ptr->p_best_sad64x32 =
                            &(context_ptr
                                  ->p_sb_best_sad[listIndex][ref_pic_index]
                                                 [ME_TIER_ZERO_PU_64x32_0]);
                        context_ptr->p_best_sad32x16 =
                            &(context_ptr
                                  ->p_sb_best_sad[listIndex][ref_pic_index]
                                                 [ME_TIER_ZERO_PU_32x16_0]);
                        context_ptr->p_best_sad16x8 = &(
                            context_ptr->p_sb_best_sad[listIndex][ref_pic_index]
                                                      [ME_TIER_ZERO_PU_16x8_0]);
                        context_ptr->p_best_sad32x64 =
                            &(context_ptr
                                  ->p_sb_best_sad[listIndex][ref_pic_index]
                                                 [ME_TIER_ZERO_PU_32x64_0]);
                        context_ptr->p_best_sad16x32 =
                            &(context_ptr
                                  ->p_sb_best_sad[listIndex][ref_pic_index]
                                                 [ME_TIER_ZERO_PU_16x32_0]);
                        context_ptr->p_best_sad8x16 = &(
                            context_ptr->p_sb_best_sad[listIndex][ref_pic_index]
                                                      [ME_TIER_ZERO_PU_8x16_0]);
                        context_ptr->p_best_sad32x8 = &(
                            context_ptr->p_sb_best_sad[listIndex][ref_pic_index]
                                                      [ME_TIER_ZERO_PU_32x8_0]);
                        context_ptr->p_best_sad8x32 = &(
                            context_ptr->p_sb_best_sad[listIndex][ref_pic_index]
                                                      [ME_TIER_ZERO_PU_8x32_0]);
                        context_ptr->p_best_sad64x16 =
                            &(context_ptr
                                  ->p_sb_best_sad[listIndex][ref_pic_index]
                                                 [ME_TIER_ZERO_PU_64x16_0]);
                        context_ptr->p_best_sad16x64 =
                            &(context_ptr
                                  ->p_sb_best_sad[listIndex][ref_pic_index]
                                                 [ME_TIER_ZERO_PU_16x64_0]);

                        context_ptr->p_best_mv64x64 = &(
                            context_ptr->p_sb_best_mv[listIndex][ref_pic_index]
                                                     [ME_TIER_ZERO_PU_64x64]);
                        context_ptr->p_best_mv32x32 = &(
                            context_ptr->p_sb_best_mv[listIndex][ref_pic_index]
                                                     [ME_TIER_ZERO_PU_32x32_0]);
                        context_ptr->p_best_mv16x16 = &(
                            context_ptr->p_sb_best_mv[listIndex][ref_pic_index]
                                                     [ME_TIER_ZERO_PU_16x16_0]);
                        context_ptr->p_best_mv8x8 = &(
                            context_ptr->p_sb_best_mv[listIndex][ref_pic_index]
                                                     [ME_TIER_ZERO_PU_8x8_0]);
                        context_ptr->p_best_mv64x32 = &(
                            context_ptr->p_sb_best_mv[listIndex][ref_pic_index]
                                                     [ME_TIER_ZERO_PU_64x32_0]);
                        context_ptr->p_best_mv32x16 = &(
                            context_ptr->p_sb_best_mv[listIndex][ref_pic_index]
                                                     [ME_TIER_ZERO_PU_32x16_0]);
                        context_ptr->p_best_mv16x8 = &(
                            context_ptr->p_sb_best_mv[listIndex][ref_pic_index]
                                                     [ME_TIER_ZERO_PU_16x8_0]);
                        context_ptr->p_best_mv32x64 = &(
                            context_ptr->p_sb_best_mv[listIndex][ref_pic_index]
                                                     [ME_TIER_ZERO_PU_32x64_0]);
                        context_ptr->p_best_mv16x32 = &(
                            context_ptr->p_sb_best_mv[listIndex][ref_pic_index]
                                                     [ME_TIER_ZERO_PU_16x32_0]);
                        context_ptr->p_best_mv8x16 = &(
                            context_ptr->p_sb_best_mv[listIndex][ref_pic_index]
                                                     [ME_TIER_ZERO_PU_8x16_0]);
                        context_ptr->p_best_mv32x8 = &(
                            context_ptr->p_sb_best_mv[listIndex][ref_pic_index]
                                                     [ME_TIER_ZERO_PU_32x8_0]);
                        context_ptr->p_best_mv8x32 = &(
                            context_ptr->p_sb_best_mv[listIndex][ref_pic_index]
                                                     [ME_TIER_ZERO_PU_8x32_0]);
                        context_ptr->p_best_mv64x16 = &(
                            context_ptr->p_sb_best_mv[listIndex][ref_pic_index]
                                                     [ME_TIER_ZERO_PU_64x16_0]);
                        context_ptr->p_best_mv16x64 = &(
                            context_ptr->p_sb_best_mv[listIndex][ref_pic_index]
                                                     [ME_TIER_ZERO_PU_16x64_0]);

                        context_ptr->p_best_ssd64x64 = &(
                            context_ptr->p_sb_best_ssd[listIndex][ref_pic_index]
                                                      [ME_TIER_ZERO_PU_64x64]);
                        context_ptr->p_best_ssd32x32 =
                            &(context_ptr
                                  ->p_sb_best_ssd[listIndex][ref_pic_index]
                                                 [ME_TIER_ZERO_PU_32x32_0]);
                        context_ptr->p_best_ssd16x16 =
                            &(context_ptr
                                  ->p_sb_best_ssd[listIndex][ref_pic_index]
                                                 [ME_TIER_ZERO_PU_16x16_0]);
                        context_ptr->p_best_ssd8x8 = &(
                            context_ptr->p_sb_best_ssd[listIndex][ref_pic_index]
                                                      [ME_TIER_ZERO_PU_8x8_0]);
                        context_ptr->p_best_ssd64x32 =
                            &(context_ptr
                                  ->p_sb_best_ssd[listIndex][ref_pic_index]
                                                 [ME_TIER_ZERO_PU_64x32_0]);
                        context_ptr->p_best_ssd32x16 =
                            &(context_ptr
                                  ->p_sb_best_ssd[listIndex][ref_pic_index]
                                                 [ME_TIER_ZERO_PU_32x16_0]);
                        context_ptr->p_best_ssd16x8 = &(
                            context_ptr->p_sb_best_ssd[listIndex][ref_pic_index]
                                                      [ME_TIER_ZERO_PU_16x8_0]);
                        context_ptr->p_best_ssd32x64 =
                            &(context_ptr
                                  ->p_sb_best_ssd[listIndex][ref_pic_index]
                                                 [ME_TIER_ZERO_PU_32x64_0]);
                        context_ptr->p_best_ssd16x32 =
                            &(context_ptr
                                  ->p_sb_best_ssd[listIndex][ref_pic_index]
                                                 [ME_TIER_ZERO_PU_16x32_0]);
                        context_ptr->p_best_ssd8x16 = &(
                            context_ptr->p_sb_best_ssd[listIndex][ref_pic_index]
                                                      [ME_TIER_ZERO_PU_8x16_0]);
                        context_ptr->p_best_ssd32x8 = &(
                            context_ptr->p_sb_best_ssd[listIndex][ref_pic_index]
                                                      [ME_TIER_ZERO_PU_32x8_0]);
                        context_ptr->p_best_ssd8x32 = &(
                            context_ptr->p_sb_best_ssd[listIndex][ref_pic_index]
                                                      [ME_TIER_ZERO_PU_8x32_0]);
                        context_ptr->p_best_ssd64x16 =
                            &(context_ptr
                                  ->p_sb_best_ssd[listIndex][ref_pic_index]
                                                 [ME_TIER_ZERO_PU_64x16_0]);
                        context_ptr->p_best_ssd16x64 =
                            &(context_ptr
                                  ->p_sb_best_ssd[listIndex][ref_pic_index]
                                                 [ME_TIER_ZERO_PU_16x64_0]);

                        open_loop_me_fullpel_search_sblock(context_ptr,
                                                           listIndex,
                                                           ref_pic_index,
                                                           x_search_area_origin,
                                                           y_search_area_origin,
                                                           search_area_width,
                                                           search_area_height,
                                                           asm_type);
                        context_ptr->full_quarter_pel_refinement = 0;

                        if (context_ptr->half_pel_mode ==
                            EX_HP_MODE) {
                            // Move to the top left of the search region
                            xTopLeftSearchRegion =
                                (int16_t)(refPicPtr->origin_x + sb_origin_x) +
                                x_search_area_origin;
                            yTopLeftSearchRegion =
                                (int16_t)(refPicPtr->origin_y + sb_origin_y) +
                                y_search_area_origin;
                            searchRegionIndex =
                                xTopLeftSearchRegion +
                                yTopLeftSearchRegion * refPicPtr->stride_y;
                            // Interpolate the search region for Half-Pel
                            // Refinements H - AVC Style
                            InterpolateSearchRegionAVC(
                                context_ptr,
                                listIndex,
                                ref_pic_index,
                                context_ptr->integer_buffer_ptr[listIndex]
                                                               [ref_pic_index] +
                                    (ME_FILTER_TAP >> 1) +
                                    ((ME_FILTER_TAP >> 1) *
                                     context_ptr->interpolated_full_stride
                                         [listIndex][ref_pic_index]),
                                context_ptr
                                    ->interpolated_full_stride[listIndex]
                                                              [ref_pic_index],
                                (uint32_t)search_area_width +
                                    (BLOCK_SIZE_64 - 1),
                                (uint32_t)search_area_height +
                                    (BLOCK_SIZE_64 - 1),
                                8);

                            initialize_buffer_32bits(
                                context_ptr
                                    ->p_sb_best_ssd[listIndex][ref_pic_index],
                                52,
                                1,
                                MAX_SSE_VALUE);
                            memcpy(context_ptr
                                       ->p_sb_best_full_pel_mv[listIndex]
                                                              [ref_pic_index],
                                   context_ptr
                                       ->p_sb_best_mv[listIndex][ref_pic_index],
                                   MAX_ME_PU_COUNT * sizeof(uint32_t));
                            context_ptr->full_quarter_pel_refinement = 1;
                            context_ptr->p_best_full_pel_mv64x64 =
                                &(context_ptr->p_sb_best_full_pel_mv
                                      [listIndex][ref_pic_index]
                                      [ME_TIER_ZERO_PU_64x64]);
                            context_ptr->p_best_full_pel_mv32x32 =
                                &(context_ptr->p_sb_best_full_pel_mv
                                      [listIndex][ref_pic_index]
                                      [ME_TIER_ZERO_PU_32x32_0]);
                            context_ptr->p_best_full_pel_mv16x16 =
                                &(context_ptr->p_sb_best_full_pel_mv
                                      [listIndex][ref_pic_index]
                                      [ME_TIER_ZERO_PU_16x16_0]);
                            context_ptr->p_best_full_pel_mv8x8 =
                                &(context_ptr->p_sb_best_full_pel_mv
                                      [listIndex][ref_pic_index]
                                      [ME_TIER_ZERO_PU_8x8_0]);
                            context_ptr->p_best_full_pel_mv64x32 =
                                &(context_ptr->p_sb_best_full_pel_mv
                                      [listIndex][ref_pic_index]
                                      [ME_TIER_ZERO_PU_64x32_0]);
                            context_ptr->p_best_full_pel_mv32x16 =
                                &(context_ptr->p_sb_best_full_pel_mv
                                      [listIndex][ref_pic_index]
                                      [ME_TIER_ZERO_PU_32x16_0]);
                            context_ptr->p_best_full_pel_mv16x8 =
                                &(context_ptr->p_sb_best_full_pel_mv
                                      [listIndex][ref_pic_index]
                                      [ME_TIER_ZERO_PU_16x8_0]);
                            context_ptr->p_best_full_pel_mv32x64 =
                                &(context_ptr->p_sb_best_full_pel_mv
                                      [listIndex][ref_pic_index]
                                      [ME_TIER_ZERO_PU_32x64_0]);
                            context_ptr->p_best_full_pel_mv16x32 =
                                &(context_ptr->p_sb_best_full_pel_mv
                                      [listIndex][ref_pic_index]
                                      [ME_TIER_ZERO_PU_16x32_0]);
                            context_ptr->p_best_full_pel_mv8x16 =
                                &(context_ptr->p_sb_best_full_pel_mv
                                      [listIndex][ref_pic_index]
                                      [ME_TIER_ZERO_PU_8x16_0]);
                            context_ptr->p_best_full_pel_mv32x8 =
                                &(context_ptr->p_sb_best_full_pel_mv
                                      [listIndex][ref_pic_index]
                                      [ME_TIER_ZERO_PU_32x8_0]);
                            context_ptr->p_best_full_pel_mv8x32 =
                                &(context_ptr->p_sb_best_full_pel_mv
                                      [listIndex][ref_pic_index]
                                      [ME_TIER_ZERO_PU_8x32_0]);
                            context_ptr->p_best_full_pel_mv64x16 =
                                &(context_ptr->p_sb_best_full_pel_mv
                                      [listIndex][ref_pic_index]
                                      [ME_TIER_ZERO_PU_64x16_0]);
                            context_ptr->p_best_full_pel_mv16x64 =
                                &(context_ptr->p_sb_best_full_pel_mv
                                      [listIndex][ref_pic_index]
                                      [ME_TIER_ZERO_PU_16x64_0]);
                            // half-Pel search
                            open_loop_me_half_pel_search_sblock(
                                picture_control_set_ptr,
                                context_ptr,
                                listIndex,
                                ref_pic_index,
                                x_search_area_origin,
                                y_search_area_origin,
                                search_area_width,
                                search_area_height);
                        }

                        if (context_ptr->quarter_pel_mode ==
                            EX_QP_MODE) {
                            // Quarter-Pel search
                            memcpy(context_ptr
                                       ->p_sb_best_full_pel_mv[listIndex]
                                                              [ref_pic_index],
                                   context_ptr
                                       ->p_sb_best_mv[listIndex][ref_pic_index],
                                   MAX_ME_PU_COUNT * sizeof(uint32_t));
                            open_loop_me_quarter_pel_search_sblock(
                                context_ptr,
                                listIndex,
                                ref_pic_index,
                                x_search_area_origin,
                                y_search_area_origin,
                                search_area_width,
                                search_area_height);
                        }
                    } else {
                        initialize_buffer_32bits(
                            context_ptr
                                ->p_sb_best_sad[listIndex][ref_pic_index],
                            21,
                            1,
                            MAX_SAD_VALUE);

                        context_ptr->full_quarter_pel_refinement = 0;

                        context_ptr->p_best_sad64x64 = &(
                            context_ptr->p_sb_best_sad[listIndex][ref_pic_index]
                                                      [ME_TIER_ZERO_PU_64x64]);
                        context_ptr->p_best_sad32x32 =
                            &(context_ptr
                                  ->p_sb_best_sad[listIndex][ref_pic_index]
                                                 [ME_TIER_ZERO_PU_32x32_0]);
                        context_ptr->p_best_sad16x16 =
                            &(context_ptr
                                  ->p_sb_best_sad[listIndex][ref_pic_index]
                                                 [ME_TIER_ZERO_PU_16x16_0]);
                        context_ptr->p_best_sad8x8 = &(
                            context_ptr->p_sb_best_sad[listIndex][ref_pic_index]
                                                      [ME_TIER_ZERO_PU_8x8_0]);

                        context_ptr->p_best_mv64x64 = &(
                            context_ptr->p_sb_best_mv[listIndex][ref_pic_index]
                                                     [ME_TIER_ZERO_PU_64x64]);
                        context_ptr->p_best_mv32x32 = &(
                            context_ptr->p_sb_best_mv[listIndex][ref_pic_index]
                                                     [ME_TIER_ZERO_PU_32x32_0]);
                        context_ptr->p_best_mv16x16 = &(
                            context_ptr->p_sb_best_mv[listIndex][ref_pic_index]
                                                     [ME_TIER_ZERO_PU_16x16_0]);
                        context_ptr->p_best_mv8x8 = &(
                            context_ptr->p_sb_best_mv[listIndex][ref_pic_index]
                                                     [ME_TIER_ZERO_PU_8x8_0]);

                        context_ptr->p_best_ssd64x64 = &(
                            context_ptr->p_sb_best_ssd[listIndex][ref_pic_index]
                                                      [ME_TIER_ZERO_PU_64x64]);
                        context_ptr->p_best_ssd32x32 =
                            &(context_ptr
                                  ->p_sb_best_ssd[listIndex][ref_pic_index]
                                                 [ME_TIER_ZERO_PU_32x32_0]);
                        context_ptr->p_best_ssd16x16 =
                            &(context_ptr
                                  ->p_sb_best_ssd[listIndex][ref_pic_index]
                                                 [ME_TIER_ZERO_PU_16x16_0]);
                        context_ptr->p_best_ssd8x8 = &(
                            context_ptr->p_sb_best_ssd[listIndex][ref_pic_index]
                                                      [ME_TIER_ZERO_PU_8x8_0]);
                        FullPelSearch_LCU(context_ptr,
                                          listIndex,
                                          ref_pic_index,
                                          x_search_area_origin,
                                          y_search_area_origin,
                                          search_area_width,
                                          search_area_height,
                                          asm_type);
                    }
                }

                if (context_ptr->fractional_search_model == 0) {
                    enableHalfPel32x32 = EB_TRUE;
                    enableHalfPel16x16 = EB_TRUE;
                    enableHalfPel8x8 = EB_TRUE;
                    enableQuarterPel = EB_TRUE;
                } else if (context_ptr->fractional_search_model == 1) {
                    suPelEnable(context_ptr,
                                picture_control_set_ptr,
                                listIndex,
                                0,
                                &enableHalfPel32x32,
                                &enableHalfPel16x16,
                                &enableHalfPel8x8);
                    enableQuarterPel = EB_TRUE;
                } else {
                    enableHalfPel32x32 = EB_FALSE;
                    enableHalfPel16x16 = EB_FALSE;
                    enableHalfPel8x8 = EB_FALSE;
                    enableQuarterPel = EB_FALSE;
                }
                if (enableHalfPel32x32 || enableHalfPel16x16 ||
                    enableHalfPel8x8 || enableQuarterPel) {
                    // if((picture_control_set_ptr->is_used_as_reference_flag ==
                    // EB_TRUE)) {
                    // Move to the top left of the search region
                    xTopLeftSearchRegion =
                        (int16_t)(refPicPtr->origin_x + sb_origin_x) +
                        x_search_area_origin;
                    yTopLeftSearchRegion =
                        (int16_t)(refPicPtr->origin_y + sb_origin_y) +
                        y_search_area_origin;
                    searchRegionIndex =
                        xTopLeftSearchRegion +
                        yTopLeftSearchRegion * refPicPtr->stride_y;

                    // Interpolate the search region for Half-Pel Refinements
                    // H - AVC Style

                    if (context_ptr->half_pel_mode ==
                        REFINMENT_HP_MODE) {
                        InterpolateSearchRegionAVC(
                            context_ptr,
                            listIndex,
                            ref_pic_index,
                            context_ptr->integer_buffer_ptr[listIndex]
                                                           [ref_pic_index] +
                                (ME_FILTER_TAP >> 1) +
                                ((ME_FILTER_TAP >> 1) *
                                 context_ptr
                                     ->interpolated_full_stride[listIndex]
                                                               [ref_pic_index]),
                            context_ptr
                                ->interpolated_full_stride[listIndex]
                                                          [ref_pic_index],
                            (uint32_t)search_area_width + (BLOCK_SIZE_64 - 1),
                            (uint32_t)search_area_height + (BLOCK_SIZE_64 - 1),
                            8);

                        // Half-Pel Refinement [8 search positions]
                        HalfPelSearch_LCU(
                            sequence_control_set_ptr,
                            picture_control_set_ptr,
                            context_ptr,
#if M0_HIGH_PRECISION_INTERPOLATION
                            context_ptr->integer_buffer_ptr[listIndex]
                                                           [ref_pic_index] +
                                (ME_FILTER_PAD_DISTANCE >> 1) +
                                ((ME_FILTER_PAD_DISTANCE >> 1) *
                                 context_ptr
                                     ->interpolated_full_stride[listIndex]
                                                               [ref_pic_index]),
                            context_ptr
                                ->interpolated_full_stride[listIndex]
                                                          [ref_pic_index],
                            &(context_ptr->pos_b_buffer
                                  [listIndex][ref_pic_index]
                                  [(ME_FILTER_PAD_DISTANCE >> 1) *
                                   context_ptr->interpolated_stride]),
#else
                            context_ptr->integer_buffer_ptr[listIndex]
                                                           [ref_pic_index] +
                                (ME_FILTER_TAP >> 1) +
                                ((ME_FILTER_TAP >> 1) *
                                 context_ptr
                                     ->interpolated_full_stride[listIndex]
                                                               [ref_pic_index]),
                            context_ptr
                                ->interpolated_full_stride[listIndex]
                                                          [ref_pic_index],
                            &(context_ptr->pos_b_buffer
                                  [listIndex][ref_pic_index]
                                  [(ME_FILTER_TAP >> 1) *
                                   context_ptr->interpolated_stride]),
#endif
                            &(context_ptr
                                  ->pos_h_buffer[listIndex][ref_pic_index][1]),
                            &(context_ptr
                                  ->pos_j_buffer[listIndex][ref_pic_index][0]),
                            x_search_area_origin,
                            y_search_area_origin,
                            picture_control_set_ptr->cu8x8_mode ==
                                CU_8x8_MODE_1,
                            enableHalfPel32x32,
                            enableHalfPel16x16,
                            enableHalfPel8x8);
                    }

                    if (context_ptr->quarter_pel_mode ==
                        REFINMENT_QP_MODE) {
                        // Quarter-Pel Refinement [8 search positions]
                        QuarterPelSearch_LCU(
                            context_ptr,
#if M0_HIGH_PRECISION_INTERPOLATION
                            context_ptr->integer_buffer_ptr[listIndex]
                                                           [ref_pic_index] +
                                (ME_FILTER_PAD_DISTANCE >> 1) +
                                ((ME_FILTER_PAD_DISTANCE >> 1) *
                                 context_ptr
                                     ->interpolated_full_stride[listIndex]
                                                               [ref_pic_index]),
                            context_ptr
                                ->interpolated_full_stride[listIndex]
                                                          [ref_pic_index],
                            &(context_ptr->pos_b_buffer
                                  [listIndex][ref_pic_index]
                                  [(ME_FILTER_PAD_DISTANCE >> 1) *
                                   context_ptr
                                       ->interpolated_stride]),  // points to b
                                                                 // position of
                                                                 // the figure
                                                                 // above
#else
                            context_ptr->integer_buffer_ptr[listIndex]
                                                           [ref_pic_index] +
                                (ME_FILTER_TAP >> 1) +
                                ((ME_FILTER_TAP >> 1) *
                                 context_ptr
                                     ->interpolated_full_stride[listIndex]
                                                               [ref_pic_index]),
                            context_ptr
                                ->interpolated_full_stride[listIndex]
                                                          [ref_pic_index],
                            &(context_ptr->pos_b_buffer
                                  [listIndex][ref_pic_index]
                                  [(ME_FILTER_TAP >> 1) *
                                   context_ptr
                                       ->interpolated_stride]),  // points to b
                                                                 // position of
                                                                 // the figure
                                                                 // above
#endif
                            &(context_ptr
                                  ->pos_h_buffer[listIndex][ref_pic_index]
                                                [1]),  // points to h position
                                                       // of the figure above
                            &(context_ptr
                                  ->pos_j_buffer[listIndex][ref_pic_index]
                                                [0]),  // points to j position
                                                       // of the figure above
                            x_search_area_origin,
                            y_search_area_origin,
                            picture_control_set_ptr->cu8x8_mode ==
                                CU_8x8_MODE_1,
                            enableHalfPel32x32,
                            enableHalfPel16x16,
                            enableHalfPel8x8,
                            enableQuarterPel,
                            picture_control_set_ptr->pic_depth_mode <=
                                PIC_ALL_C_DEPTH_MODE);
                    }
                }
                if (is_nsq_table_used && ref_pic_index == 0) {
                    context_ptr->p_best_nsq64x64 =
                        &(context_ptr->p_sb_best_nsq[listIndex][0]
                                                    [ME_TIER_ZERO_PU_64x64]);
                    context_ptr->p_best_nsq32x32 =
                        &(context_ptr->p_sb_best_nsq[listIndex][0]
                                                    [ME_TIER_ZERO_PU_32x32_0]);
                    context_ptr->p_best_nsq16x16 =
                        &(context_ptr->p_sb_best_nsq[listIndex][0]
                                                    [ME_TIER_ZERO_PU_16x16_0]);
                    context_ptr->p_best_nsq8x8 =
                        &(context_ptr->p_sb_best_nsq[listIndex][0]
                                                    [ME_TIER_ZERO_PU_8x8_0]);
                    nsq_get_analysis_results_block(context_ptr);
                }
        }
    }
}

if (context_ptr->me_alt_ref == EB_FALSE) {

    // Bi-Prediction motion estimation loop
    for (pu_index = 0; pu_index < max_number_of_pus_per_sb; ++pu_index) {
        candidateIndex = 0;

        uint32_t nIdx;

        if (pu_index > 200)
            nIdx = pu_index;
        else if (pu_index > 184)
            nIdx = tab8x32[pu_index - 185] + 185;
        else if (pu_index > 168)
            nIdx = tab32x8[pu_index - 169] + 169;
        else if (pu_index > 136)
            nIdx = tab8x16[pu_index - 137] + 137;
        else if (pu_index > 128)
            nIdx = tab16x32[pu_index - 129] + 129;
        else if (pu_index > 126)
            nIdx = pu_index;
        else if (pu_index > 94)
            nIdx = tab16x8[pu_index - 95] + 95;
        else if (pu_index > 86)
            nIdx = tab32x16[pu_index - 87] + 87;
        else if (pu_index > 84)
            nIdx = pu_index;
        else if (pu_index > 20)
            nIdx = tab8x8[pu_index - 21] + 21;
        else if (pu_index > 4)
            nIdx = tab16x16[pu_index - 5] + 5;
        else
            nIdx = pu_index;
        for (listIndex = REF_LIST_0; listIndex <= numOfListToSearch;
             ++listIndex) {
            num_of_ref_pic_to_search =
                (picture_control_set_ptr->slice_type == P_SLICE)
                    ? picture_control_set_ptr->ref_list0_count
                    : (listIndex == REF_LIST_0)
                          ? picture_control_set_ptr->ref_list0_count
                          : picture_control_set_ptr->ref_list1_count;

            // Ref Picture Loop
            for (ref_pic_index = 0; ref_pic_index < num_of_ref_pic_to_search;
                 ++ref_pic_index) {
                me_candidate =
                    &(context_ptr->me_candidate[candidateIndex].pu[pu_index]);
                me_candidate->prediction_direction = listIndex;
                me_candidate->ref_index[listIndex] = ref_pic_index;
                me_candidate->ref0_list = me_candidate->prediction_direction == 0 ? listIndex : 24;
                me_candidate->ref1_list = me_candidate->prediction_direction == 1 ? listIndex : 24;
                me_candidate->distortion =
                    context_ptr->p_sb_best_sad[listIndex][ref_pic_index][nIdx];
                candidateIndex++;
            }
        }

        total_me_candidate_index = candidateIndex;
        uint8_t ref_type_table[7];
        if (picture_control_set_ptr->prune_unipred_at_me) {
            // Sorting of the ME candidates
            for (candidate_index = 0;
                candidate_index < total_me_candidate_index - 1;
                ++candidate_index) {
                for (next_candidate_index = candidate_index + 1;
                    next_candidate_index < total_me_candidate_index;
                    ++next_candidate_index) {
                    if (context_ptr->me_candidate[candidate_index]
                        .pu[pu_index]
                        .distortion >
                        context_ptr->me_candidate[next_candidate_index]
                        .pu[pu_index]
                        .distortion) {
                        SwapMeCandidate(
                            &(context_ptr->me_candidate[candidate_index]
                                .pu[pu_index]),
                            &(context_ptr->me_candidate[next_candidate_index]
                                .pu[pu_index]));
                    }
                }
            }
            for (candidate_index = 0;
                candidate_index < total_me_candidate_index;
                ++candidate_index) {

                me_candidate =
                    &(context_ptr->me_candidate[candidate_index].pu[pu_index]);

                if (me_candidate->prediction_direction == 0)
                    ref_type_table[candidate_index] = svt_get_ref_frame_type(me_candidate->ref0_list, me_candidate->ref_index[0]);
                else
                    ref_type_table[candidate_index] = svt_get_ref_frame_type(me_candidate->ref1_list, me_candidate->ref_index[1]);

            }
        }
        if (numOfListToSearch) {
            if (picture_control_set_ptr->cu8x8_mode == CU_8x8_MODE_0 ||
                pu_index < 21 ||
                (picture_control_set_ptr->pic_depth_mode <=
                 PIC_ALL_C_DEPTH_MODE)) {
                BiPredictionSearch(
                    sequence_control_set_ptr,
                    context_ptr,
                    pu_index,
                    candidateIndex,
                    picture_control_set_ptr->ref_list0_count,
                    picture_control_set_ptr->ref_list1_count,
                    &total_me_candidate_index,
                    asm_type,
                    ref_type_table,
                    picture_control_set_ptr);
            }
        }

        // Sorting of the ME candidates
        for (candidate_index = 0;
             candidate_index < total_me_candidate_index - 1;
             ++candidate_index) {
            for (next_candidate_index = candidate_index + 1;
                 next_candidate_index < total_me_candidate_index;
                 ++next_candidate_index) {
                if (context_ptr->me_candidate[candidate_index]
                        .pu[pu_index]
                        .distortion >
                    context_ptr->me_candidate[next_candidate_index]
                        .pu[pu_index]
                        .distortion) {
                    SwapMeCandidate(
                        &(context_ptr->me_candidate[candidate_index]
                              .pu[pu_index]),
                        &(context_ptr->me_candidate[next_candidate_index]
                              .pu[pu_index]));
                }
            }
        }

        MeLcuResults *mePuResult =
            picture_control_set_ptr->me_results[sb_index];
        mePuResult->total_me_candidate_index[pu_index] =
            total_me_candidate_index;

        uint8_t l0_nsq =
            is_nsq_table_used ? context_ptr->p_sb_best_nsq[0][0][nIdx] : 0;
        uint8_t l1_nsq =
            is_nsq_table_used ? context_ptr->p_sb_best_nsq[1][0][nIdx] : 0;
        mePuResult->me_nsq_0[pu_index] = l0_nsq;
        mePuResult->me_nsq_1[pu_index] = l1_nsq;

        mePuResult->total_me_candidate_index[pu_index] =
            MIN(total_me_candidate_index, ME_RES_CAND_MRP_MODE_0);
        // Assining the ME candidates to the me Results buffer
        for (candidateIndex = 0; candidateIndex < total_me_candidate_index;
             ++candidateIndex) {
            me_candidate =
                &(context_ptr->me_candidate[candidateIndex].pu[pu_index]);
            picture_control_set_ptr->me_results[sb_index]
                ->me_candidate[pu_index][candidateIndex]
                .distortion = me_candidate->distortion;
            picture_control_set_ptr->me_results[sb_index]
                ->me_candidate[pu_index][candidateIndex]
                .direction = me_candidate->prediction_direction;
            picture_control_set_ptr->me_results[sb_index]
                ->me_candidate[pu_index][candidateIndex]
                .ref_idx_l0 = me_candidate->ref_index[0];
            picture_control_set_ptr->me_results[sb_index]
                ->me_candidate[pu_index][candidateIndex]
                .ref_idx_l1 = me_candidate->ref_index[1];
            picture_control_set_ptr->me_results[sb_index]
                ->me_candidate[pu_index][candidateIndex]
                .ref0_list = me_candidate->ref0_list;
            picture_control_set_ptr->me_results[sb_index]
                ->me_candidate[pu_index][candidateIndex]
                .ref1_list = me_candidate->ref1_list;
        }

        for (listIndex = REF_LIST_0; listIndex <= numOfListToSearch;
             ++listIndex) {
            num_of_ref_pic_to_search =
                (picture_control_set_ptr->slice_type == P_SLICE)
                    ? picture_control_set_ptr->ref_list0_count
                    : (listIndex == REF_LIST_0)
                          ? picture_control_set_ptr->ref_list0_count
                          : picture_control_set_ptr->ref_list1_count;

            // Ref Picture Loop
            for (ref_pic_index = 0; ref_pic_index < num_of_ref_pic_to_search;
                 ++ref_pic_index) {
                picture_control_set_ptr->me_results[sb_index]
                    ->me_mv_array[pu_index]
                                 [((listIndex &&
                                    sequence_control_set_ptr->mrp_mode == 0)
                                       ? 4
                                       : listIndex ? 2 : 0) +
                                  ref_pic_index]
                    .x_mv = _MVXT(
                    context_ptr->p_sb_best_mv[listIndex][ref_pic_index][nIdx]);
                picture_control_set_ptr->me_results[sb_index]
                    ->me_mv_array[pu_index]
                                 [((listIndex &&
                                    sequence_control_set_ptr->mrp_mode == 0)
                                       ? 4
                                       : listIndex ? 2 : 0) +
                                  ref_pic_index]
                    .y_mv = _MVYT(
                    context_ptr->p_sb_best_mv[listIndex][ref_pic_index][nIdx]);
            }
        }
    }
    {
        // Compute the sum of the distortion of all 16 16x16 (best) blocks
        // in the LCU
        picture_control_set_ptr->rc_me_distortion[sb_index] = 0;
        for (i = 0; i < 16; i++)
            picture_control_set_ptr->rc_me_distortion[sb_index] +=
                picture_control_set_ptr->me_results[sb_index]
                    ->me_candidate[5 + i][0]
                    .distortion;
    }

}

return return_error;
}

/*******************************************
 * SixteenthDecimatedSearch
 *  performs a 1/16 decimated search
 *******************************************/
uint64_t SixteenthDecimatedSearch(MeContext *context_ptr, int16_t origin_x,
                                  int16_t origin_y, uint32_t sb_width,
                                  uint32_t sb_height,
                                  EbPictureBufferDesc *sixteenthRefPicPtr,
                                  int16_t search_area_width,
                                  int16_t search_area_height, EbAsm asm_type) {
    int16_t xTopLeftSearchRegion;
    int16_t yTopLeftSearchRegion;
    uint32_t searchRegionIndex;
    int16_t x_search_area_origin;
    int16_t y_search_area_origin;

    int16_t padWidth = (int16_t)(sixteenthRefPicPtr->origin_x) - 1;
    int16_t padHeight = (int16_t)(sixteenthRefPicPtr->origin_y) - 1;

    uint64_t best_sad;
    int16_t x_search_center;
    int16_t y_search_center;

    x_search_area_origin = -(search_area_width >> 1);
    y_search_area_origin = -(search_area_height >> 1);

    // Correct the left edge of the Search Area if it is not on the reference
    // Picture
    x_search_area_origin = ((origin_x + x_search_area_origin) < -padWidth)
                               ? -padWidth - origin_x
                               : x_search_area_origin;

    search_area_width =
        ((origin_x + x_search_area_origin) < -padWidth)
            ? search_area_width -
                  (-padWidth - (origin_x + x_search_area_origin))
            : search_area_width;

    // Correct the right edge of the Search Area if its not on the reference
    // Picture
    x_search_area_origin =
        ((origin_x + x_search_area_origin) >
         (int16_t)sixteenthRefPicPtr->width - 1)
            ? x_search_area_origin - ((origin_x + x_search_area_origin) -
                                      ((int16_t)sixteenthRefPicPtr->width - 1))
            : x_search_area_origin;

    search_area_width =
        ((origin_x + x_search_area_origin + search_area_width) >
         (int16_t)sixteenthRefPicPtr->width)
            ? MAX(1,
                  search_area_width -
                      ((origin_x + x_search_area_origin + search_area_width) -
                       (int16_t)sixteenthRefPicPtr->width))
            : search_area_width;

    // Correct the top edge of the Search Area if it is not on the reference
    // Picture
    y_search_area_origin = ((origin_y + y_search_area_origin) < -padHeight)
                               ? -padHeight - origin_y
                               : y_search_area_origin;

    search_area_height =
        ((origin_y + y_search_area_origin) < -padHeight)
            ? search_area_height -
                  (-padHeight - (origin_y + y_search_area_origin))
            : search_area_height;

    // Correct the bottom edge of the Search Area if its not on the reference
    // Picture
    y_search_area_origin =
        ((origin_y + y_search_area_origin) >
         (int16_t)sixteenthRefPicPtr->height - 1)
            ? y_search_area_origin - ((origin_y + y_search_area_origin) -
                                      ((int16_t)sixteenthRefPicPtr->height - 1))
            : y_search_area_origin;

    search_area_height =
        (origin_y + y_search_area_origin + search_area_height >
         (int16_t)sixteenthRefPicPtr->height)
            ? MAX(1,
                  search_area_height -
                      ((origin_y + y_search_area_origin + search_area_height) -
                       (int16_t)sixteenthRefPicPtr->height))
            : search_area_height;

    xTopLeftSearchRegion = ((int16_t)sixteenthRefPicPtr->origin_x + origin_x) +
                           x_search_area_origin;
    yTopLeftSearchRegion = ((int16_t)sixteenthRefPicPtr->origin_y + origin_y) +
                           y_search_area_origin;
    searchRegionIndex = xTopLeftSearchRegion +
                        yTopLeftSearchRegion * sixteenthRefPicPtr->stride_y;

    if (((search_area_width & 15) == 0) && (asm_type == ASM_AVX2)) {
        sad_loop_kernel_avx2_hme_l0_intrin(
            &context_ptr->sixteenth_sb_buffer[0],
            context_ptr->sixteenth_sb_buffer_stride * 2,
            &sixteenthRefPicPtr->buffer_y[searchRegionIndex],
            sixteenthRefPicPtr->stride_y * 2,
            sb_height >> 1,
            sb_width,
            // results
            &best_sad,
            &x_search_center,
            &y_search_center,
            // range
            sixteenthRefPicPtr->stride_y,
            search_area_width,
            search_area_height);
    } else if ((search_area_width & 15) == 0) {
        // Only width equals 16 (LCU equals 64) is updated
        // other width sizes work with the old code as the one
        // in"sad_loop_kernel_sse4_1_intrin"
        sad_loop_kernel_sse4_1_hme_l0_intrin(
            &context_ptr->sixteenth_sb_buffer[0],
            context_ptr->sixteenth_sb_buffer_stride * 2,
            &sixteenthRefPicPtr->buffer_y[searchRegionIndex],
            sixteenthRefPicPtr->stride_y * 2,
            sb_height >> 1,
            sb_width,
            /* results */
            &best_sad,
            &x_search_center,
            &y_search_center,
            /* range */
            sixteenthRefPicPtr->stride_y,
            search_area_width,
            search_area_height);
    } else {
        // Put the first search location into level0 results
        sad_loop_kernel(
            &context_ptr->sixteenth_sb_buffer[0],
            context_ptr->sixteenth_sb_buffer_stride * 2,
            &sixteenthRefPicPtr->buffer_y[searchRegionIndex],
            sixteenthRefPicPtr->stride_y * 2,
            sb_height >> 1,
            sb_width,
            /* results */
            &best_sad,
            &x_search_center,
            &y_search_center,
            /* range */
            sixteenthRefPicPtr->stride_y,
            search_area_width,
            search_area_height);
    }

    return (best_sad);
}

/*******************************************
 * IsComplexLcu
 *   returns true is the SB has a high spatial & temporal complexity
 *******************************************/
EbBool IsComplexLcu(PictureParentControlSet *previousParentPcs,
                    PictureParentControlSet *currentParentPcs,
                    PictureParentControlSet *plusOneParentPcs,
                    uint32_t pictureWidthInLcus, uint32_t lcuAdrr,
                    uint32_t sb_origin_x, uint32_t sb_origin_y,
                    uint32_t sb_width, uint32_t sb_height,
                    uint32_t lcuCollocatedSad) {
    uint32_t availableLcusCount = 0;
    uint32_t highVarianceLcusCount = 0;

    // Check the variance of the current LCU
    if ((currentParentPcs->variance[lcuAdrr][ME_TIER_ZERO_PU_64x64]) >
        IS_COMPLEX_LCU_VARIANCE_TH) {
        availableLcusCount++;
        highVarianceLcusCount++;
    }

    // Check the variance of left SB if available
    if (sb_origin_x != 0) {
        availableLcusCount++;
        if ((currentParentPcs->variance[lcuAdrr - 1][ME_TIER_ZERO_PU_64x64]) >
            IS_COMPLEX_LCU_VARIANCE_TH)
            highVarianceLcusCount++;
    }

    // Check the variance of right SB if available
    if ((sb_origin_x + BLOCK_SIZE_64) <
        currentParentPcs->enhanced_picture_ptr->width) {
        availableLcusCount++;
        if ((currentParentPcs->variance[lcuAdrr + 1][ME_TIER_ZERO_PU_64x64]) >
            IS_COMPLEX_LCU_VARIANCE_TH)
            highVarianceLcusCount++;
    }

    // Check the variance of top SB if available
    if (sb_origin_y != 0) {
        availableLcusCount++;
        if ((currentParentPcs->variance[lcuAdrr - pictureWidthInLcus]
                                       [ME_TIER_ZERO_PU_64x64]) >
            IS_COMPLEX_LCU_VARIANCE_TH)
            highVarianceLcusCount++;
    }

    // Check the variance of bottom LCU
    if ((sb_origin_y + BLOCK_SIZE_64) <
        currentParentPcs->enhanced_picture_ptr->height) {
        availableLcusCount++;
        if ((currentParentPcs->variance[lcuAdrr + pictureWidthInLcus]
                                       [ME_TIER_ZERO_PU_64x64]) >
            IS_COMPLEX_LCU_VARIANCE_TH)
            highVarianceLcusCount++;
    }

    // Check the variance of top-left LCU
    if ((sb_origin_x >= BLOCK_SIZE_64) && (sb_origin_y >= BLOCK_SIZE_64)) {
        availableLcusCount++;
        if ((currentParentPcs->variance[lcuAdrr - pictureWidthInLcus - 1]
                                       [ME_TIER_ZERO_PU_64x64]) >
            IS_COMPLEX_LCU_VARIANCE_TH)
            highVarianceLcusCount++;
    }

    // Check the variance of top-right LCU
    if ((sb_origin_x <
         currentParentPcs->enhanced_picture_ptr->width - BLOCK_SIZE_64) &&
        (sb_origin_y >= BLOCK_SIZE_64)) {
        availableLcusCount++;
        if ((currentParentPcs->variance[lcuAdrr - pictureWidthInLcus + 1]
                                       [ME_TIER_ZERO_PU_64x64]) >
            IS_COMPLEX_LCU_VARIANCE_TH)
            highVarianceLcusCount++;
    }

    // Check the variance of bottom-left LCU
    if ((sb_origin_x >= BLOCK_SIZE_64) &&
        (sb_origin_y <
         currentParentPcs->enhanced_picture_ptr->height - BLOCK_SIZE_64)) {
        availableLcusCount++;
        if ((currentParentPcs->variance[lcuAdrr + pictureWidthInLcus - 1]
                                       [ME_TIER_ZERO_PU_64x64]) >
            IS_COMPLEX_LCU_VARIANCE_TH)
            highVarianceLcusCount++;
    }

    // Check the variance of bottom-right LCU
    if ((sb_origin_x <
         currentParentPcs->enhanced_picture_ptr->width - BLOCK_SIZE_64) &&
        (sb_origin_y <
         currentParentPcs->enhanced_picture_ptr->height - BLOCK_SIZE_64)) {
        availableLcusCount++;
        if ((currentParentPcs->variance[lcuAdrr + pictureWidthInLcus + 1]
                                       [ME_TIER_ZERO_PU_64x64]) >
            IS_COMPLEX_LCU_VARIANCE_TH)
            highVarianceLcusCount++;
    }

    EbBool varianceFluctuateFlag = EB_FALSE;

    if ((previousParentPcs->variance[lcuAdrr][ME_TIER_ZERO_PU_64x64]) >
            IS_COMPLEX_LCU_FLAT_VARIANCE_TH &&
        (currentParentPcs->variance[lcuAdrr][ME_TIER_ZERO_PU_64x64]) >
            IS_COMPLEX_LCU_FLAT_VARIANCE_TH &&
        (plusOneParentPcs->variance[lcuAdrr][ME_TIER_ZERO_PU_64x64]) >
            IS_COMPLEX_LCU_FLAT_VARIANCE_TH) {
        varianceFluctuateFlag = (EbBool)(
            (((ABS((int32_t)currentParentPcs
                       ->variance[lcuAdrr][ME_TIER_ZERO_PU_64x64] -
                   (int32_t)previousParentPcs
                       ->variance[lcuAdrr][ME_TIER_ZERO_PU_64x64]) *
               100) /
              (int32_t)previousParentPcs
                  ->variance[lcuAdrr][ME_TIER_ZERO_PU_64x64]) >=
             IS_COMPLEX_LCU_VARIANCE_DEVIATION_TH) &&
            (((ABS((int32_t)currentParentPcs
                       ->variance[lcuAdrr][ME_TIER_ZERO_PU_64x64] -
                   (int32_t)plusOneParentPcs
                       ->variance[lcuAdrr][ME_TIER_ZERO_PU_64x64]) *
               100) /
              (int32_t)
                  plusOneParentPcs->variance[lcuAdrr][ME_TIER_ZERO_PU_64x64]) >=
             IS_COMPLEX_LCU_VARIANCE_DEVIATION_TH));
    }

    if (lcuCollocatedSad >=
            ((sb_width * sb_height) * IS_COMPLEX_LCU_ZZ_SAD_FACTOR_TH) &&
        highVarianceLcusCount >= (availableLcusCount >> 1) &&
        varianceFluctuateFlag) {
        return EB_TRUE;
    }

    return EB_FALSE;
}

EbErrorType open_loop_intra_search_sb(
    PictureParentControlSet *picture_control_set_ptr, uint32_t sb_index,
    MotionEstimationContext_t *context_ptr, EbPictureBufferDesc *input_ptr)
{
    EbErrorType return_error = EB_ErrorNone;
    SequenceControlSet *sequence_control_set_ptr =
        (SequenceControlSet *)picture_control_set_ptr
            ->sequence_control_set_wrapper_ptr->object_ptr;

    uint32_t cu_origin_x;
    uint32_t cu_origin_y;
    uint32_t pa_blk_index = 0;
#if !PAETH_HBD
    uint8_t is_16_bit =
        (sequence_control_set_ptr->static_config.encoder_bit_depth > EB_8BIT);
#endif
    SbParams *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];
    OisSbResults *ois_sb_results_ptr =
        picture_control_set_ptr->ois_sb_results[sb_index];
    uint8_t *above_row;
    uint8_t *left_col;

    DECLARE_ALIGNED(16, uint8_t, left_data[MAX_TX_SIZE * 2 + 32]);
    DECLARE_ALIGNED(16, uint8_t, above_data[MAX_TX_SIZE * 2 + 32]);
    while (pa_blk_index < CU_MAX_COUNT) {
        const CodedUnitStats *blk_stats_ptr;
        blk_stats_ptr = get_coded_unit_stats(pa_blk_index);
        uint8_t bsize = blk_stats_ptr->size;
        TxSize tx_size =
            bsize == 8
                ? TX_8X8
                : bsize == 16 ? TX_16X16 : bsize == 32 ? TX_32X32 : TX_64X64;
        if (sb_params->raster_scan_cu_validity
                [md_scan_to_raster_scan[pa_blk_index]]) {
            OisCandidate *ois_blk_ptr =
                ois_sb_results_ptr->ois_candidate_array[pa_blk_index];
            cu_origin_x = sb_params->origin_x + blk_stats_ptr->origin_x;
            cu_origin_y = sb_params->origin_y + blk_stats_ptr->origin_y;
            above_row = above_data + 16;
            left_col = left_data + 16;

            // Fill Neighbor Arrays
            update_neighbor_samples_array_open_loop(above_row - 1,
                                                    left_col - 1,
                                                    input_ptr,
                                                    input_ptr->stride_y,
                                                    cu_origin_x,
                                                    cu_origin_y,
                                                    bsize,
                                                    bsize);
            uint8_t ois_intra_mode;
            uint8_t ois_intra_count = 0;
            uint8_t best_intra_ois_index = 0;
            uint32_t best_intra_ois_distortion = 64 * 64 * 255;
            uint8_t intra_mode_start = DC_PRED;
#if PAETH_HBD
            uint8_t intra_mode_end = PAETH_PRED;
#else
            uint8_t intra_mode_end = is_16_bit ? SMOOTH_H_PRED : PAETH_PRED;
#endif
            uint8_t angle_delta_counter = 0;
            uint8_t angle_delta_shift = 1;
            EbBool use_angle_delta = (bsize >= 8);
            uint8_t angle_delta_candidate_count = use_angle_delta ? 7 : 1;
            uint8_t disable_angular_prediction = 0;
            if (picture_control_set_ptr->intra_pred_mode == 5) {
                intra_mode_end =
                    (picture_control_set_ptr->is_used_as_reference_flag == 0)
                        ? DC_PRED
                        : intra_mode_end;
                disable_angular_prediction =
                    picture_control_set_ptr->temporal_layer_index > 0
                        ? 1
                        : (bsize > 16) ? 1 : 0;
                angle_delta_candidate_count =
                    disable_angular_prediction ? 1 : use_angle_delta ? 5 : 1;
                angle_delta_shift = 1;
            }
            else if (picture_control_set_ptr->intra_pred_mode == 6) {
                intra_mode_end =
                    (picture_control_set_ptr->is_used_as_reference_flag == 0)
                        ? DC_PRED
                        : intra_mode_end;
                disable_angular_prediction =
                    picture_control_set_ptr->temporal_layer_index > 0
                        ? 1
                        : (bsize > 16) ? 1 : 0;
                angle_delta_candidate_count = 1;
                angle_delta_shift = 1;
            } else {
                if (picture_control_set_ptr->slice_type == I_SLICE) {
#if PAETH_HBD
                    intra_mode_end = /*is_16_bit ? SMOOTH_H_PRED :*/ PAETH_PRED;
#else
                    intra_mode_end = is_16_bit ? SMOOTH_H_PRED : PAETH_PRED;
#endif
                    angle_delta_candidate_count = use_angle_delta ? 5 : 1;
                    disable_angular_prediction = 0;
                    angle_delta_shift = 1;
                } else if (picture_control_set_ptr->temporal_layer_index == 0) {
#if PAETH_HBD
                    intra_mode_end = /*is_16_bit ? SMOOTH_H_PRED :*/ PAETH_PRED;
#else
                    intra_mode_end = is_16_bit ? SMOOTH_H_PRED : PAETH_PRED;
#endif
                    angle_delta_candidate_count =
                        (bsize > 16) ? 1 : use_angle_delta ? 2 : 1;
                    disable_angular_prediction = 0;
                    angle_delta_shift = 3;
                } else {
                    intra_mode_end = DC_PRED;
                    disable_angular_prediction = 1;
                    angle_delta_candidate_count = 1;
                    angle_delta_shift = 1;
                }
            }
            for (ois_intra_mode = intra_mode_start;
                 ois_intra_mode <= intra_mode_end;
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
                                     : angle_delta_counter -
                                           (angle_delta_candidate_count >> 1));
                            int32_t p_angle =
                                mode_to_angle_map[(
                                    PredictionMode)ois_intra_mode] +
                                angle_delta * ANGLE_STEP;
                            // PRED
                            intra_prediction_open_loop(p_angle,
                                                       ois_intra_mode,
                                                       cu_origin_x,
                                                       cu_origin_y,
                                                       tx_size,
                                                       above_row,
                                                       left_col,
                                                       context_ptr);
                            // Distortion
                            ois_blk_ptr[ois_intra_count].distortion =
                                (uint32_t)nxm_sad_kernel(  // Always SAD without weighting
                                        &(input_ptr->buffer_y
                                              [(input_ptr->origin_y +
                                                cu_origin_y) *
                                                   input_ptr->stride_y +
                                               (input_ptr->origin_x +
                                                cu_origin_x)]),
                                        input_ptr->stride_y,
                                        &(context_ptr->me_context_ptr
                                              ->sb_buffer[0]),
                                        BLOCK_SIZE_64,
                                        bsize,
                                        bsize);
                            // kepp track of best SAD
                            if (ois_blk_ptr[ois_intra_count].distortion <
                                best_intra_ois_distortion) {
                                best_intra_ois_index = ois_intra_count;
                                best_intra_ois_distortion =
                                    ois_blk_ptr[ois_intra_count].distortion;
                            }
                            ois_blk_ptr[ois_intra_count].intra_mode =
                                ois_intra_mode;
                            ois_blk_ptr[ois_intra_count].valid_distortion =
                                EB_TRUE;
                            ois_blk_ptr[ois_intra_count++].angle_delta =
                                angle_delta;
                        }
                    }
                } else {
                    // PRED
                    intra_prediction_open_loop(0,
                                               ois_intra_mode,
                                               cu_origin_x,
                                               cu_origin_y,
                                               tx_size,
                                               above_row,
                                               left_col,
                                               context_ptr);
                    // Distortion
                    ois_blk_ptr[ois_intra_count]
                        .distortion = (uint32_t)nxm_sad_kernel(  // Always SAD without weighting
                            &(input_ptr->buffer_y
                                  [(input_ptr->origin_y + cu_origin_y) *
                                       input_ptr->stride_y +
                                   (input_ptr->origin_x + cu_origin_x)]),
                            input_ptr->stride_y,
                            &(context_ptr->me_context_ptr->sb_buffer[0]),
                            BLOCK_SIZE_64,
                            bsize,
                            bsize);
                    // kepp track of best SAD
                    if (ois_blk_ptr[ois_intra_count].distortion <
                        best_intra_ois_distortion) {
                        best_intra_ois_index = ois_intra_count;
                        best_intra_ois_distortion =
                            ois_blk_ptr[ois_intra_count].distortion;
                    }
                    ois_blk_ptr[ois_intra_count].intra_mode = ois_intra_mode;
                    ois_blk_ptr[ois_intra_count].valid_distortion = EB_TRUE;
                    ois_blk_ptr[ois_intra_count++].angle_delta = 0;
                }
            }
            ois_sb_results_ptr->best_distortion_index[pa_blk_index] =
                best_intra_ois_index;
            ois_sb_results_ptr->total_ois_intra_candidate[pa_blk_index] =
                ois_intra_count;
        }
        pa_blk_index++;
    }
    return return_error;
}
