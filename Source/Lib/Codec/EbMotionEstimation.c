/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdio.h>

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
#include <math.h>
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
#include "EbPictureOperators.h"
#endif
#define OIS_TH_COUNT    4

int32_t OisPointTh[3][MAX_TEMPORAL_LAYERS][OIS_TH_COUNT] = {
    {
        // Light OIS
        { -20, 50, 150, 200 },
        { -20, 50, 150, 200 },
        { -20, 50, 100, 150 },
        { -20, 50, 200, 300 },
        { -20, 50, 200, 300 },
        { -20, 50, 200, 300 }
    },
    {
        // Default OIS
        { -150, 0, 150, 200 },
        { -150, 0, 150, 200 },
        { -125, 0, 100, 150 },
        { -50, 50, 200, 300 },
        { -50, 50, 200, 300 },
        { -50, 50, 200, 300 }
    }
    ,
    {
        // Heavy OIS
        { -400, -300, -200, 0 },
        { -400, -300, -200, 0 },
        { -400, -300, -200, 0 },
        { -400, -300, -200, 0 },
        { -400, -300, -200, 0 },
        { -400, -300, -200, 0 }
    }
};





#define AVCCODEL
/********************************************
* Constants
********************************************/

static const uint8_t numberOfOisModePoints[5/*IOS poitnt*/] = {
    1, 3, 5, 7, 9
};

#define MAX_INTRA_IN_MD   9
//CHKN: Note that the max number set in this table should be referenced by the upper macro!!!
static uint8_t intraSearchInMd[5][4] = {
    /*depth0*/    /*depth1*/    /*depth2*/    /*depth3*/
    { 0, 1, 1, 1 },  // Very fast
    { 0, 3, 3, 3 },  // fast
    { 0, 5, 5, 5 },  // Meduim
    { 0, 7, 7, 7 },  // Complex
    { 0, MAX_INTRA_IN_MD, MAX_INTRA_IN_MD, MAX_INTRA_IN_MD }  // Very Complex
};

// Intra Open Loop
uint32_t iSliceModesArray[11] = { EB_INTRA_PLANAR, EB_INTRA_DC, EB_INTRA_HORIZONTAL, EB_INTRA_VERTICAL, EB_INTRA_MODE_2, EB_INTRA_MODE_18, EB_INTRA_MODE_34, EB_INTRA_MODE_6, EB_INTRA_MODE_14, EB_INTRA_MODE_22, EB_INTRA_MODE_30 };
uint32_t stage1ModesArray[9] = { EB_INTRA_HORIZONTAL, EB_INTRA_VERTICAL, EB_INTRA_MODE_2, EB_INTRA_MODE_18, EB_INTRA_MODE_34, EB_INTRA_MODE_6, EB_INTRA_MODE_14, EB_INTRA_MODE_22, EB_INTRA_MODE_30 };

#define REFERENCE_PIC_LIST_0  0
#define REFERENCE_PIC_LIST_1  1


/*******************************************
* Compute8x4SAD_Default
*   Unoptimized 8x4 SAD
*******************************************/
uint32_t Compute8x4SAD_Kernel(
    uint8_t  *src,                            // input parameter, source samples Ptr
    uint32_t  src_stride,                      // input parameter, source stride
    uint8_t  *ref,                            // input parameter, reference samples Ptr
    uint32_t  refStride)                      // input parameter, reference stride
{
    uint32_t rowNumberInBlock8x4;
    uint32_t sadBlock8x4 = 0;

    for (rowNumberInBlock8x4 = 0; rowNumberInBlock8x4 < 4; ++rowNumberInBlock8x4) {
        sadBlock8x4 += EB_ABS_DIFF(src[0x00], ref[0x00]);
        sadBlock8x4 += EB_ABS_DIFF(src[0x01], ref[0x01]);
        sadBlock8x4 += EB_ABS_DIFF(src[0x02], ref[0x02]);
        sadBlock8x4 += EB_ABS_DIFF(src[0x03], ref[0x03]);
        sadBlock8x4 += EB_ABS_DIFF(src[0x04], ref[0x04]);
        sadBlock8x4 += EB_ABS_DIFF(src[0x05], ref[0x05]);
        sadBlock8x4 += EB_ABS_DIFF(src[0x06], ref[0x06]);
        sadBlock8x4 += EB_ABS_DIFF(src[0x07], ref[0x07]);
        src += src_stride;
        ref += refStride;
    }

    return sadBlock8x4;
}
static EB_COMPUTE8X4SAD_TYPE FUNC_TABLE compute8x4SAD_funcPtrArray[ASM_TYPE_TOTAL] =// [C_DEFAULT/ASM]
{
    // C_DEFAULT
    Compute8x4SAD_Kernel,
    // SSE2
    Compute8x4SAD_Kernel,
};
/***************************************
* Function Tables
***************************************/
static EB_EXTSADCALCULATION8X8AND16X16_TYPE ExtSadCalculation_8x8_16x16_funcPtrArray[ASM_TYPE_TOTAL] = {
    // NON_AVX2
    ExtSadCalculation_8x8_16x16,
    // AVX2
    ExtSadCalculation_8x8_16x16_SSE4_INTRIN
};
static EB_EXTSADCALCULATION32X32AND64X64_TYPE ExtSadCalculation_32x32_64x64_funcPtrArray[ASM_TYPE_TOTAL] = {
    // NON_AVX2
    ExtSadCalculation_32x32_64x64,
    // AVX2
    ExtSadCalculation_32x32_64x64_SSE4_INTRIN
};
static EB_SADCALCULATION8X8AND16X16_TYPE SadCalculation_8x8_16x16_funcPtrArray[ASM_TYPE_TOTAL] = {
    // NON_AVX2
    SadCalculation_8x8_16x16_SSE2_INTRIN,
    // AVX2
    SadCalculation_8x8_16x16_SSE2_INTRIN,
};
static EB_SADCALCULATION32X32AND64X64_TYPE SadCalculation_32x32_64x64_funcPtrArray[ASM_TYPE_TOTAL] = {
    // NON_AVX2
    SadCalculation_32x32_64x64_SSE2_INTRIN,
    // AVX2
    SadCalculation_32x32_64x64_SSE2_INTRIN,
};

/*******************************************
Calcualte SAD for 16x16 and its 8x8 sublcoks
and check if there is improvment, if yes keep
the best SAD+MV
*******************************************/
void ExtSadCalculation_8x8_16x16(
    uint8_t   *src,
    uint32_t   src_stride,
    uint8_t   *ref,
    uint32_t   refStride,
    uint32_t  *p_best_sad8x8,
    uint32_t  *p_best_sad16x16,
    uint32_t  *p_best_mv8x8,
    uint32_t  *p_best_mv16x16,
    uint32_t   mv,
    uint32_t  *p_sad16x16,
    uint32_t  *p_sad8x8)
{
    uint32_t sad8x8_0, sad8x8_1, sad8x8_2, sad8x8_3;
    uint32_t sad16x16;

    uint32_t   srcStrideSub = (src_stride << 1); //TODO get these from outside
    uint32_t   refStrideSub = (refStride << 1);


    p_sad8x8[0] = sad8x8_0 = (compute8x4SAD_funcPtrArray[0](src, srcStrideSub, ref, refStrideSub)) << 1;
    if (sad8x8_0 < p_best_sad8x8[0]) {
        p_best_sad8x8[0] = (uint32_t)sad8x8_0;
        p_best_mv8x8[0] = mv;
    }

    p_sad8x8[1] = sad8x8_1 = (compute8x4SAD_funcPtrArray[0](src + 8, srcStrideSub, ref + 8, refStrideSub)) << 1;
    if (sad8x8_1 < p_best_sad8x8[1]) {
        p_best_sad8x8[1] = (uint32_t)sad8x8_1;
        p_best_mv8x8[1] = mv;
    }

    p_sad8x8[2] = sad8x8_2 = (compute8x4SAD_funcPtrArray[0](src + (src_stride << 3), srcStrideSub, ref + (refStride << 3), refStrideSub)) << 1;
    if (sad8x8_2 < p_best_sad8x8[2]) {
        p_best_sad8x8[2] = (uint32_t)sad8x8_2;
        p_best_mv8x8[2] = mv;
    }

    p_sad8x8[3] = sad8x8_3 = (compute8x4SAD_funcPtrArray[0](src + (src_stride << 3) + 8, srcStrideSub, ref + (refStride << 3) + 8, refStrideSub)) << 1;
    if (sad8x8_3 < p_best_sad8x8[3]) {
        p_best_sad8x8[3] = (uint32_t)sad8x8_3;
        p_best_mv8x8[3] = mv;
    }

    sad16x16 = sad8x8_0 + sad8x8_1 + sad8x8_2 + sad8x8_3;
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
void ExtSadCalculation_32x32_64x64(
    uint32_t  *p_sad16x16,
    uint32_t  *p_best_sad32x32,
    uint32_t  *p_best_sad64x64,
    uint32_t  *p_best_mv32x32,
    uint32_t  *p_best_mv64x64,
    uint32_t   mv,
    uint32_t  *p_sad32x32)
{

    uint32_t sad32x32_0, sad32x32_1, sad32x32_2, sad32x32_3, sad64x64;

    p_sad32x32[0] = sad32x32_0 = p_sad16x16[0] + p_sad16x16[1] + p_sad16x16[2] + p_sad16x16[3];
    if (sad32x32_0 < p_best_sad32x32[0]) {
        p_best_sad32x32[0] = sad32x32_0;
        p_best_mv32x32[0] = mv;
    }

    p_sad32x32[1] = sad32x32_1 = p_sad16x16[4] + p_sad16x16[5] + p_sad16x16[6] + p_sad16x16[7];
    if (sad32x32_1 < p_best_sad32x32[1]) {
        p_best_sad32x32[1] = sad32x32_1;
        p_best_mv32x32[1] = mv;
    }

    p_sad32x32[2] = sad32x32_2 = p_sad16x16[8] + p_sad16x16[9] + p_sad16x16[10] + p_sad16x16[11];
    if (sad32x32_2 < p_best_sad32x32[2]) {
        p_best_sad32x32[2] = sad32x32_2;
        p_best_mv32x32[2] = mv;
    }

    p_sad32x32[3] = sad32x32_3 = p_sad16x16[12] + p_sad16x16[13] + p_sad16x16[14] + p_sad16x16[15];
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

/****************************************************
Calcualte SAD for Rect H, V and H4, V4 partitions

and update its Motion info if the result SAD is better
****************************************************/
void ExtSadCalculation(
    uint32_t  *p_sad8x8,
    uint32_t  *p_sad16x16,
    uint32_t  *p_sad32x32,
    uint32_t  *p_best_sad64x32,
    uint32_t  *p_best_mv64x32,
    uint32_t  *p_best_sad32x16,
    uint32_t  *p_best_mv32x16,
    uint32_t  *p_best_sad16x8,
    uint32_t  *p_best_mv16x8,
    uint32_t  *p_best_sad32x64,
    uint32_t  *p_best_mv32x64,
    uint32_t  *p_best_sad16x32,
    uint32_t  *p_best_mv16x32,
    uint32_t  *p_best_sad8x16,
    uint32_t  *p_best_mv8x16,
    uint32_t  *p_best_sad32x8,
    uint32_t  *p_best_mv32x8,
    uint32_t  *p_best_sad8x32,
    uint32_t  *p_best_mv8x32,
    uint32_t  *p_best_sad64x16,
    uint32_t  *p_best_mv64x16,
    uint32_t  *p_best_sad16x64,
    uint32_t  *p_best_mv16x64,
    uint32_t   mv)
{
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
    if (sad < p_best_sad32x16[5]) {
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
static EB_EXTSADCALCULATION_TYPE ExtSadCalculation_funcPtrArray[ASM_TYPE_TOTAL] = {
    // Should be written in Assembly
    // C_DEFAULT
    ExtSadCalculation,
    // Assembly
    ExtSadCalculation
};

/*******************************************
* open_loop_me_get_search_point_results_block
*******************************************/
static void open_loop_me_get_search_point_results_block(
    MeContext_t             *context_ptr,                    // input parameter, ME context Ptr, used to get SB Ptr
    uint32_t                   listIndex,                     // input parameter, reference list index
    uint32_t                   searchRegionIndex,             // input parameter, search area origin, used to point to reference samples
    int32_t                   xSearchIndex,                  // input parameter, search region position in the horizontal direction, used to derive xMV
    int32_t                   ySearchIndex,                  // input parameter, search region position in the vertical direction, used to derive yMV
    EbAsm                   asm_type)
{
    uint8_t  *srcPtr = context_ptr->sb_src_ptr;

    // uint8_t  *refPtr = refPicPtr->bufferY; // NADER
    uint8_t  *refPtr = context_ptr->integer_buffer_ptr[listIndex][0] + (ME_FILTER_TAP >> 1) + ((ME_FILTER_TAP >> 1) * context_ptr->interpolated_full_stride[listIndex][0]);

    // uint32_t reflumaStride = refPicPtr->strideY; // NADER
    uint32_t reflumaStride = context_ptr->interpolated_full_stride[listIndex][0];
    uint32_t searchPositionTLIndex = searchRegionIndex;
    uint32_t searchPositionIndex;
    uint32_t blockIndex;
    uint32_t srcNext16x16Offset = (BLOCK_SIZE_64 << 4);
    //uint32_t refNext16x16Offset = (refPicPtr->strideY << 4); // NADER
    uint32_t refNext16x16Offset = (reflumaStride << 4);
    uint32_t   currMV1 = (((uint16_t)ySearchIndex) << 18);
    uint16_t   currMV2 = (((uint16_t)xSearchIndex << 2));
    uint32_t   currMV = currMV1 | currMV2;
    uint32_t  *p_best_sad8x8 = context_ptr->p_best_sad8x8;
    uint32_t  *p_best_sad16x16 = context_ptr->p_best_sad16x16;
    uint32_t  *p_best_sad32x32 = context_ptr->p_best_sad32x32;
    uint32_t  *p_best_sad64x64 = context_ptr->p_best_sad64x64;
    uint32_t  *p_best_sad64x32 = context_ptr->p_best_sad64x32;
    uint32_t  *p_best_sad32x16 = context_ptr->p_best_sad32x16;
    uint32_t  *p_best_sad16x8 = context_ptr->p_best_sad16x8;
    uint32_t  *p_best_sad32x64 = context_ptr->p_best_sad32x64;
    uint32_t  *p_best_sad16x32 = context_ptr->p_best_sad16x32;
    uint32_t  *p_best_sad8x16 = context_ptr->p_best_sad8x16;
    uint32_t  *p_best_sad32x8 = context_ptr->p_best_sad32x8;
    uint32_t  *p_best_sad8x32 = context_ptr->p_best_sad8x32;
    uint32_t  *p_best_sad64x16 = context_ptr->p_best_sad64x16;
    uint32_t  *p_best_sad16x64 = context_ptr->p_best_sad16x64;
    uint32_t  *p_best_mv8x8 = context_ptr->p_best_mv8x8;
    uint32_t  *p_best_mv16x16 = context_ptr->p_best_mv16x16;
    uint32_t  *p_best_mv32x32 = context_ptr->p_best_mv32x32;
    uint32_t  *p_best_mv64x64 = context_ptr->p_best_mv64x64;
    uint32_t  *p_best_mv64x32 = context_ptr->p_best_mv64x32;
    uint32_t  *p_best_mv32x16 = context_ptr->p_best_mv32x16;
    uint32_t  *p_best_mv16x8 = context_ptr->p_best_mv16x8;
    uint32_t  *p_best_mv32x64 = context_ptr->p_best_mv32x64;
    uint32_t  *p_best_mv16x32 = context_ptr->p_best_mv16x32;
    uint32_t  *p_best_mv8x16 = context_ptr->p_best_mv8x16;
    uint32_t  *p_best_mv32x8 = context_ptr->p_best_mv32x8;
    uint32_t  *p_best_mv8x32 = context_ptr->p_best_mv8x32;
    uint32_t  *p_sad32x32 = context_ptr->p_sad32x32;
    uint32_t  *p_sad16x16 = context_ptr->p_sad16x16;
    uint32_t  *p_sad8x8 = context_ptr->p_sad8x8;
    uint32_t  *p_best_mv64x16 = context_ptr->p_best_mv64x16;
    uint32_t  *p_best_mv16x64 = context_ptr->p_best_mv16x64;

    //TODO: blockIndex searchPositionIndex could be removed  + Connect asm_type
    (void)asm_type;

    const uint32_t  src_stride = context_ptr->sb_src_stride;
    srcNext16x16Offset = src_stride << 4;

    //---- 16x16 : 0
    blockIndex = 0;
    searchPositionIndex = searchPositionTLIndex;

    ExtSadCalculation_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[0], &p_best_sad16x16[0], &p_best_mv8x8[0], &p_best_mv16x16[0], currMV, &p_sad16x16[0], &p_sad8x8[0]);

    //---- 16x16 : 1
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionTLIndex + 16;
    ExtSadCalculation_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[4], &p_best_sad16x16[1], &p_best_mv8x8[4], &p_best_mv16x16[1], currMV, &p_sad16x16[1], &p_sad8x8[4]);
    //---- 16x16 : 4
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;

    ExtSadCalculation_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[16], &p_best_sad16x16[4], &p_best_mv8x8[16], &p_best_mv16x16[4], currMV, &p_sad16x16[4], &p_sad8x8[16]);


    //---- 16x16 : 5
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    ExtSadCalculation_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[20], &p_best_sad16x16[5], &p_best_mv8x8[20], &p_best_mv16x16[5], currMV, &p_sad16x16[5], &p_sad8x8[20]);


    //---- 16x16 : 2
    blockIndex = srcNext16x16Offset;
    searchPositionIndex = searchPositionTLIndex + refNext16x16Offset;
    ExtSadCalculation_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[8], &p_best_sad16x16[2], &p_best_mv8x8[8], &p_best_mv16x16[2], currMV, &p_sad16x16[2], &p_sad8x8[8]);
    //---- 16x16 : 3
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    ExtSadCalculation_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[12], &p_best_sad16x16[3], &p_best_mv8x8[12], &p_best_mv16x16[3], currMV, &p_sad16x16[3], &p_sad8x8[12]);
    //---- 16x16 : 6
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    ExtSadCalculation_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[24], &p_best_sad16x16[6], &p_best_mv8x8[24], &p_best_mv16x16[6], currMV, &p_sad16x16[6], &p_sad8x8[24]);
    //---- 16x16 : 7
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    ExtSadCalculation_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[28], &p_best_sad16x16[7], &p_best_mv8x8[28], &p_best_mv16x16[7], currMV, &p_sad16x16[7], &p_sad8x8[28]);


    //---- 16x16 : 8
    blockIndex = (srcNext16x16Offset << 1);
    searchPositionIndex = searchPositionTLIndex + (refNext16x16Offset << 1);
    ExtSadCalculation_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[32], &p_best_sad16x16[8], &p_best_mv8x8[32], &p_best_mv16x16[8], currMV, &p_sad16x16[8], &p_sad8x8[32]);
    //---- 16x16 : 9
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    ExtSadCalculation_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[36], &p_best_sad16x16[9], &p_best_mv8x8[36], &p_best_mv16x16[9], currMV, &p_sad16x16[9], &p_sad8x8[36]);
    //---- 16x16 : 12
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    ExtSadCalculation_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[48], &p_best_sad16x16[12], &p_best_mv8x8[48], &p_best_mv16x16[12], currMV, &p_sad16x16[12], &p_sad8x8[48]);
    //---- 16x16 : 13
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    ExtSadCalculation_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[52], &p_best_sad16x16[13], &p_best_mv8x8[52], &p_best_mv16x16[13], currMV, &p_sad16x16[13], &p_sad8x8[52]);


    //---- 16x16 : 10
    blockIndex = (srcNext16x16Offset * 3);
    searchPositionIndex = searchPositionTLIndex + (refNext16x16Offset * 3);
    ExtSadCalculation_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[40], &p_best_sad16x16[10], &p_best_mv8x8[40], &p_best_mv16x16[10], currMV, &p_sad16x16[10], &p_sad8x8[40]);
    //---- 16x16 : 11
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    ExtSadCalculation_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[44], &p_best_sad16x16[11], &p_best_mv8x8[44], &p_best_mv16x16[11], currMV, &p_sad16x16[11], &p_sad8x8[44]);
    //---- 16x16 : 14
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    ExtSadCalculation_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[56], &p_best_sad16x16[14], &p_best_mv8x8[56], &p_best_mv16x16[14], currMV, &p_sad16x16[14], &p_sad8x8[56]);
    //---- 16x16 : 15
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    ExtSadCalculation_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[60], &p_best_sad16x16[15], &p_best_mv8x8[60], &p_best_mv16x16[15], currMV, &p_sad16x16[15], &p_sad8x8[60]);

    ExtSadCalculation_32x32_64x64_funcPtrArray[asm_type](p_sad16x16, p_best_sad32x32, p_best_sad64x64, p_best_mv32x32, p_best_mv64x64, currMV, &p_sad32x32[0]);

    ExtSadCalculation_funcPtrArray[asm_type](
        p_sad8x8,
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
    MeContext_t             *context_ptr,                    // input parameter, ME context Ptr, used to get SB Ptr
    uint32_t                   listIndex,                     // input parameter, reference list index
    uint32_t                   searchRegionIndex,             // input parameter, search area origin, used to point to reference samples
    int32_t                   xSearchIndex,                  // input parameter, search region position in the horizontal direction, used to derive xMV
    int32_t                   ySearchIndex,                  // input parameter, search region position in the vertical direction, used to derive yMV
    EbAsm                   asm_type)
{
    uint8_t  *srcPtr = context_ptr->sb_src_ptr;

    // uint8_t  *refPtr = refPicPtr->bufferY; // NADER
    uint8_t  *refPtr = context_ptr->integer_buffer_ptr[listIndex][0] + (ME_FILTER_TAP >> 1) + ((ME_FILTER_TAP >> 1) * context_ptr->interpolated_full_stride[listIndex][0]);

    // uint32_t reflumaStride = refPicPtr->strideY; // NADER
    uint32_t reflumaStride = context_ptr->interpolated_full_stride[listIndex][0];

    uint32_t searchPositionTLIndex = searchRegionIndex;
    uint32_t searchPositionIndex;
    uint32_t blockIndex;

    uint32_t srcNext16x16Offset = (BLOCK_SIZE_64 << 4);
    //uint32_t refNext16x16Offset = (refPicPtr->strideY << 4); // NADER
    uint32_t refNext16x16Offset = (reflumaStride << 4);

    uint32_t currMV1 = (((uint16_t)ySearchIndex) << 18);
    uint16_t currMV2 = (((uint16_t)xSearchIndex << 2));
    uint32_t currMV = currMV1 | currMV2;


    uint32_t  *p_best_sad8x8 = context_ptr->p_best_sad8x8;
    uint32_t  *p_best_sad16x16 = context_ptr->p_best_sad16x16;
    uint32_t  *p_best_sad32x32 = context_ptr->p_best_sad32x32;
    uint32_t  *p_best_sad64x64 = context_ptr->p_best_sad64x64;

    uint32_t  *p_best_mv8x8 = context_ptr->p_best_mv8x8;
    uint32_t  *p_best_mv16x16 = context_ptr->p_best_mv16x16;
    uint32_t  *p_best_mv32x32 = context_ptr->p_best_mv32x32;
    uint32_t  *p_best_mv64x64 = context_ptr->p_best_mv64x64;
    uint32_t  *p_sad16x16 = context_ptr->p_sad16x16;


    //TODO: blockIndex searchPositionIndex could be removed  + Connect asm_type
    (void)asm_type;

    const uint32_t  src_stride = context_ptr->sb_src_stride;
    srcNext16x16Offset = src_stride << 4;

    //---- 16x16 : 0
    blockIndex = 0;
    searchPositionIndex = searchPositionTLIndex;

    SadCalculation_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[0], &p_best_sad16x16[0], &p_best_mv8x8[0], &p_best_mv16x16[0], currMV, &p_sad16x16[0]);

    //---- 16x16 : 1
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionTLIndex + 16;
    SadCalculation_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[4], &p_best_sad16x16[1], &p_best_mv8x8[4], &p_best_mv16x16[1], currMV, &p_sad16x16[1]);
    //---- 16x16 : 4
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;

    SadCalculation_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[16], &p_best_sad16x16[4], &p_best_mv8x8[16], &p_best_mv16x16[4], currMV, &p_sad16x16[4]);


    //---- 16x16 : 5
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    SadCalculation_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[20], &p_best_sad16x16[5], &p_best_mv8x8[20], &p_best_mv16x16[5], currMV, &p_sad16x16[5]);


    //---- 16x16 : 2
    blockIndex = srcNext16x16Offset;
    searchPositionIndex = searchPositionTLIndex + refNext16x16Offset;
    SadCalculation_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[8], &p_best_sad16x16[2], &p_best_mv8x8[8], &p_best_mv16x16[2], currMV, &p_sad16x16[2]);
    //---- 16x16 : 3
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    SadCalculation_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[12], &p_best_sad16x16[3], &p_best_mv8x8[12], &p_best_mv16x16[3], currMV, &p_sad16x16[3]);
    //---- 16x16 : 6
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    SadCalculation_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[24], &p_best_sad16x16[6], &p_best_mv8x8[24], &p_best_mv16x16[6], currMV, &p_sad16x16[6]);
    //---- 16x16 : 7
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    SadCalculation_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[28], &p_best_sad16x16[7], &p_best_mv8x8[28], &p_best_mv16x16[7], currMV, &p_sad16x16[7]);


    //---- 16x16 : 8
    blockIndex = (srcNext16x16Offset << 1);
    searchPositionIndex = searchPositionTLIndex + (refNext16x16Offset << 1);
    SadCalculation_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[32], &p_best_sad16x16[8], &p_best_mv8x8[32], &p_best_mv16x16[8], currMV, &p_sad16x16[8]);
    //---- 16x16 : 9
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    SadCalculation_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[36], &p_best_sad16x16[9], &p_best_mv8x8[36], &p_best_mv16x16[9], currMV, &p_sad16x16[9]);
    //---- 16x16 : 12
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    SadCalculation_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[48], &p_best_sad16x16[12], &p_best_mv8x8[48], &p_best_mv16x16[12], currMV, &p_sad16x16[12]);
    //---- 16x16 : 13
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    SadCalculation_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[52], &p_best_sad16x16[13], &p_best_mv8x8[52], &p_best_mv16x16[13], currMV, &p_sad16x16[13]);


    //---- 16x16 : 10
    blockIndex = (srcNext16x16Offset * 3);
    searchPositionIndex = searchPositionTLIndex + (refNext16x16Offset * 3);
    SadCalculation_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[40], &p_best_sad16x16[10], &p_best_mv8x8[40], &p_best_mv16x16[10], currMV, &p_sad16x16[10]);
    //---- 16x16 : 11
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    SadCalculation_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[44], &p_best_sad16x16[11], &p_best_mv8x8[44], &p_best_mv16x16[11], currMV, &p_sad16x16[11]);
    //---- 16x16 : 14
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    SadCalculation_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[56], &p_best_sad16x16[14], &p_best_mv8x8[56], &p_best_mv16x16[14], currMV, &p_sad16x16[14]);
    //---- 16x16 : 15
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    SadCalculation_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[60], &p_best_sad16x16[15], &p_best_mv8x8[60], &p_best_mv16x16[15], currMV, &p_sad16x16[15]);



    SadCalculation_32x32_64x64_funcPtrArray[asm_type](p_sad16x16, p_best_sad32x32, p_best_sad64x64, p_best_mv32x32, p_best_mv64x64, currMV);

}

/*******************************************
* GetEightHorizontalSearchPointResultsAll85CUs
*******************************************/
static void GetEightHorizontalSearchPointResultsAll85PUs(
    MeContext_t             *context_ptr,
    uint32_t                   listIndex,
    uint32_t                   searchRegionIndex,
    int32_t                   xSearchIndex,                  // input parameter, search region position in the horizontal direction, used to derive xMV
    int32_t                   ySearchIndex,                  // input parameter, search region position in the vertical direction, used to derive yMV
    EbAsm                   asm_type
)
{
    uint8_t  *srcPtr = context_ptr->sb_src_ptr;
    uint8_t  *refPtr = context_ptr->integer_buffer_ptr[listIndex][0] + (ME_FILTER_TAP >> 1) + ((ME_FILTER_TAP >> 1) * context_ptr->interpolated_full_stride[listIndex][0]);
    uint32_t reflumaStride = context_ptr->interpolated_full_stride[listIndex][0];

    uint32_t searchPositionTLIndex = searchRegionIndex;
    uint32_t searchPositionIndex;
    uint32_t blockIndex;

    uint32_t srcNext16x16Offset = (BLOCK_SIZE_64 << 4);
    uint32_t refNext16x16Offset = (reflumaStride << 4);

    uint32_t currMVy = (((uint16_t)ySearchIndex) << 18);
    uint16_t currMVx = (((uint16_t)xSearchIndex << 2));
    uint32_t currMV = currMVy | currMVx;

    uint32_t  *p_best_sad8x8 = context_ptr->p_best_sad8x8;
    uint32_t  *p_best_sad16x16 = context_ptr->p_best_sad16x16;
    uint32_t  *p_best_sad32x32 = context_ptr->p_best_sad32x32;
    uint32_t  *p_best_sad64x64 = context_ptr->p_best_sad64x64;

    uint32_t  *p_best_mv8x8 = context_ptr->p_best_mv8x8;
    uint32_t  *p_best_mv16x16 = context_ptr->p_best_mv16x16;
    uint32_t  *p_best_mv32x32 = context_ptr->p_best_mv32x32;
    uint32_t  *p_best_mv64x64 = context_ptr->p_best_mv64x64;

    uint16_t  *p_sad16x16 = context_ptr->p_eight_pos_sad16x16;

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

    const uint32_t  src_stride = context_ptr->sb_src_stride;
    srcNext16x16Offset = src_stride << 4;

    //---- 16x16_0
    blockIndex = 0;
    searchPositionIndex = searchPositionTLIndex;
    GetEightHorizontalSearchPointResults_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, context_ptr->sb_src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[0], &p_best_mv8x8[0], &p_best_sad16x16[0], &p_best_mv16x16[0], currMV, &p_sad16x16[0 * 8]);
    //---- 16x16_1
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionTLIndex + 16;
    GetEightHorizontalSearchPointResults_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, context_ptr->sb_src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[4], &p_best_mv8x8[4], &p_best_sad16x16[1], &p_best_mv16x16[1], currMV, &p_sad16x16[1 * 8]);
    //---- 16x16_4
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    GetEightHorizontalSearchPointResults_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, context_ptr->sb_src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[16], &p_best_mv8x8[16], &p_best_sad16x16[4], &p_best_mv16x16[4], currMV, &p_sad16x16[4 * 8]);
    //---- 16x16_5
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    GetEightHorizontalSearchPointResults_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, context_ptr->sb_src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[20], &p_best_mv8x8[20], &p_best_sad16x16[5], &p_best_mv16x16[5], currMV, &p_sad16x16[5 * 8]);



    //---- 16x16_2
    blockIndex = srcNext16x16Offset;
    searchPositionIndex = searchPositionTLIndex + refNext16x16Offset;
    GetEightHorizontalSearchPointResults_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, context_ptr->sb_src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[8], &p_best_mv8x8[8], &p_best_sad16x16[2], &p_best_mv16x16[2], currMV, &p_sad16x16[2 * 8]);
    //---- 16x16_3
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    GetEightHorizontalSearchPointResults_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, context_ptr->sb_src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[12], &p_best_mv8x8[12], &p_best_sad16x16[3], &p_best_mv16x16[3], currMV, &p_sad16x16[3 * 8]);
    //---- 16x16_6
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    GetEightHorizontalSearchPointResults_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, context_ptr->sb_src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[24], &p_best_mv8x8[24], &p_best_sad16x16[6], &p_best_mv16x16[6], currMV, &p_sad16x16[6 * 8]);
    //---- 16x16_7
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    GetEightHorizontalSearchPointResults_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, context_ptr->sb_src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[28], &p_best_mv8x8[28], &p_best_sad16x16[7], &p_best_mv16x16[7], currMV, &p_sad16x16[7 * 8]);


    //---- 16x16_8
    blockIndex = (srcNext16x16Offset << 1);
    searchPositionIndex = searchPositionTLIndex + (refNext16x16Offset << 1);
    GetEightHorizontalSearchPointResults_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, context_ptr->sb_src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[32], &p_best_mv8x8[32], &p_best_sad16x16[8], &p_best_mv16x16[8], currMV, &p_sad16x16[8 * 8]);
    //---- 16x16_9
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    GetEightHorizontalSearchPointResults_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, context_ptr->sb_src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[36], &p_best_mv8x8[36], &p_best_sad16x16[9], &p_best_mv16x16[9], currMV, &p_sad16x16[9 * 8]);
    //---- 16x16_12
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    GetEightHorizontalSearchPointResults_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, context_ptr->sb_src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[48], &p_best_mv8x8[48], &p_best_sad16x16[12], &p_best_mv16x16[12], currMV, &p_sad16x16[12 * 8]);
    //---- 16x1_13
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    GetEightHorizontalSearchPointResults_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, context_ptr->sb_src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[52], &p_best_mv8x8[52], &p_best_sad16x16[13], &p_best_mv16x16[13], currMV, &p_sad16x16[13 * 8]);



    //---- 16x16_10
    blockIndex = (srcNext16x16Offset * 3);
    searchPositionIndex = searchPositionTLIndex + (refNext16x16Offset * 3);
    GetEightHorizontalSearchPointResults_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, context_ptr->sb_src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[40], &p_best_mv8x8[40], &p_best_sad16x16[10], &p_best_mv16x16[10], currMV, &p_sad16x16[10 * 8]);
    //---- 16x16_11
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    GetEightHorizontalSearchPointResults_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, context_ptr->sb_src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[44], &p_best_mv8x8[44], &p_best_sad16x16[11], &p_best_mv16x16[11], currMV, &p_sad16x16[11 * 8]);
    //---- 16x16_14
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    GetEightHorizontalSearchPointResults_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, context_ptr->sb_src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[56], &p_best_mv8x8[56], &p_best_sad16x16[14], &p_best_mv16x16[14], currMV, &p_sad16x16[14 * 8]);
    //---- 16x16_15
    blockIndex = blockIndex + 16;
    searchPositionIndex = searchPositionIndex + 16;
    GetEightHorizontalSearchPointResults_8x8_16x16_funcPtrArray[asm_type](srcPtr + blockIndex, context_ptr->sb_src_stride, refPtr + searchPositionIndex, reflumaStride, &p_best_sad8x8[60], &p_best_mv8x8[60], &p_best_sad16x16[15], &p_best_mv16x16[15], currMV, &p_sad16x16[15 * 8]);




    //32x32 and 64x64
    GetEightHorizontalSearchPointResults_32x32_64x64_funcPtrArray[asm_type](p_sad16x16, p_best_sad32x32, p_best_sad64x64, p_best_mv32x32, p_best_mv64x64, currMV);

}

/*******************************************
* FullPelSearch_LCU
*******************************************/
static void FullPelSearch_LCU(
    MeContext_t             *context_ptr,
    uint32_t                   listIndex,
    int16_t                   x_search_area_origin,
    int16_t                     y_search_area_origin,
    uint32_t                   search_area_width,
    uint32_t                   search_area_height,
    EbAsm                   asm_type
)
{

    uint32_t  xSearchIndex, ySearchIndex;

    uint32_t  searchAreaWidthRest8 = search_area_width & 7;
    uint32_t  searchAreaWidthMult8 = search_area_width - searchAreaWidthRest8;

    for (ySearchIndex = 0; ySearchIndex < search_area_height; ySearchIndex++) {

        for (xSearchIndex = 0; xSearchIndex < searchAreaWidthMult8; xSearchIndex += 8) {

            //this function will do:  xSearchIndex, +1, +2, ..., +7
            GetEightHorizontalSearchPointResultsAll85PUs(
                context_ptr,
                listIndex,
                xSearchIndex + ySearchIndex * context_ptr->interpolated_full_stride[listIndex][0],
                (int32_t)xSearchIndex + x_search_area_origin,
                (int32_t)ySearchIndex + y_search_area_origin,
                asm_type
            );
        }

        for (xSearchIndex = searchAreaWidthMult8; xSearchIndex < search_area_width; xSearchIndex++) {

            GetSearchPointResults(
                context_ptr,
                listIndex,
                xSearchIndex + ySearchIndex * context_ptr->interpolated_full_stride[listIndex][0],
                (int32_t)xSearchIndex + x_search_area_origin,
                (int32_t)ySearchIndex + y_search_area_origin,
                asm_type);



        }

    }

}

/*******************************************
* open_loop_me_fullpel_search_sblock
*******************************************/
static void open_loop_me_fullpel_search_sblock(
    MeContext_t             *context_ptr,
    uint32_t                   listIndex,
    int16_t                   x_search_area_origin,
    int16_t                     y_search_area_origin,
    uint32_t                   search_area_width,
    uint32_t                   search_area_height,
    EbAsm                   asm_type)
{

    uint32_t xSearchIndex, ySearchIndex;
    for (ySearchIndex = 0; ySearchIndex < search_area_height; ySearchIndex++) {

        //for (xSearchIndex = 0; xSearchIndex < searchAreaWidthMult8; xSearchIndex += 8){

        //    //this function will do:  xSearchIndex, +1, +2, ..., +7
        //    GetEightHorizontalSearchPointResultsAll85PUs(
        //        context_ptr,
        //        listIndex,
        //        xSearchIndex + ySearchIndex * context_ptr->interpolated_full_stride[listIndex][0],
        //        xSearchIndex + x_search_area_origin,
        //        ySearchIndex + y_search_area_origin,
        //        asm_type
        //        );
        //}

        for (xSearchIndex = 0/*searchAreaWidthMult8*/; xSearchIndex < search_area_width; xSearchIndex++) {

            open_loop_me_get_search_point_results_block(
                context_ptr,
                listIndex,
                xSearchIndex + ySearchIndex * context_ptr->interpolated_full_stride[listIndex][0],
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
    uint8_t   *src,           // input parameter, input samples Ptr
    uint32_t   src_stride,     // input parameter, input stride
    uint32_t   width,         // input parameter, input area width
    uint32_t   height,        // input parameter, input area height
    const int32_t  *ifCoeff, // input parameter, interpolation filter coefficients Ptr
    uint32_t   inputBitDepth, // input parameter, input sample bit depth
    uint32_t   dst_stride,     // input parameter, output stride
    uint8_t   *dst)           // output parameter, interpolated samples Ptr
{
    uint32_t x, y;
    const int32_t maxSampleValue = (1 << inputBitDepth) - 1;
    const int32_t ifOffset = 1 << (IFShift - 1);
    for (y = 0; y < height; ++y) {
        for (x = 0; x < width; ++x) {
            dst[x] = (uint8_t)CLIP3(0, (int32_t)maxSampleValue, ((
                ((int32_t)src[x] + (int32_t)src[x + 3])  * ifCoeff[0] +
                ((int32_t)src[x + 1] + (int32_t)src[x + 2])  * ifCoeff[1] + ifOffset) >> IFShift));
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
    uint8_t   *src,           // input parameter, input samples ptr
    uint32_t   src_stride,     // input parameter, input stride
    uint32_t   width,         // input parameter, input area width
    uint32_t   height,        // input parameter, input area height
    const int32_t   ifCoeff[4],    // input parameter, interpolation filter coefficients Ptr
    uint32_t   inputBitDepth, // input parameter, input sample bit depth
    uint32_t   dst_stride,     // input parameter, output stride
    uint8_t   *dst)           // output parameter, interpolated samples Ptr
{
    uint32_t x, y;

    const int32_t maxSampleValue = (1 << inputBitDepth) - 1;
    const int32_t ifOffset = 1 << (IFShift - 1);

    const uint32_t srcStride2 = src_stride << 1;
    const uint32_t srcStride3 = srcStride2 + src_stride;

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            dst[x] = (uint8_t)CLIP3(0, maxSampleValue, ((
                ((int32_t)src[x] + (int32_t)src[x + srcStride3]) * ifCoeff[0] +
                ((int32_t)src[x + src_stride] + (int32_t)src[x + srcStride2]) * ifCoeff[1] + ifOffset) >> IFShift));
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
    uint8_t   *srcOne,        // input parameter, input samples Ptr
    uint32_t   srcOneStride,  // input parameter, input stride
    uint8_t   *srcTwo,        // input parameter, input samples Ptr
    uint32_t   srcTwoStride,  // input parameter, input stride
    uint32_t   width,         // input parameter, input area width
    uint32_t   height,        // input parameter, input area height
    uint32_t   inputBitDepth, // input parameter, input sample bit depth
    uint32_t   dst_stride,     // input parameter, output stride
    uint8_t   *dst)           // output parameter, interpolated samples Ptr
{
    uint32_t x, y;
    int32_t maxSampleValue = POW2(inputBitDepth) - 1;

    for (y = 0; y < height; ++y) {
        for (x = 0; x < width; ++x) {
            dst[x] = (uint8_t)CLIP3(0, (int32_t)maxSampleValue, ((
                (int32_t)srcOne[x] +
                (int32_t)srcTwo[x] + 1) >> IFShiftAvcStyle));
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
    MeContext_t             *context_ptr,           // input/output parameter, ME context ptr, used to get/set interpolated search area Ptr
    uint32_t                   listIndex,            // Refrence picture list index
    uint8_t                   *searchRegionBuffer,   // input parameter, search region index, used to point to reference samples
    uint32_t                   lumaStride,           // input parameter, reference Picture stride
    uint32_t                   search_area_width,      // input parameter, search area width
    uint32_t                   search_area_height,     // input parameter, search area height
    uint32_t                   inputBitDepth,           // input parameter, input sample bit depth
    EbAsm                     asm_type)
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
    // Also the search area must be oversized by 2 to account for edge conditions
    uint32_t searchAreaWidthForAsm = ROUND_UP_MUL_8(search_area_width + 2);

#ifdef AVCCODEL

    (void)inputBitDepth;
    // Half pel interpolation of the search region using f1 -> pos_b_buffer
    if (searchAreaWidthForAsm) {

        avc_style_uni_pred_luma_if_function_ptr_array[asm_type][2](
            searchRegionBuffer - (ME_FILTER_TAP >> 1) * lumaStride - (ME_FILTER_TAP >> 1) + 1,
            lumaStride,
            context_ptr->pos_b_buffer[listIndex][0],
            context_ptr->interpolated_stride,
            searchAreaWidthForAsm,
            search_area_height + ME_FILTER_TAP,
            context_ptr->avctemp_buffer,
            EB_FALSE,
            2);
    }

    // Half pel interpolation of the search region using f1 -> pos_h_buffer
    if (searchAreaWidthForAsm) {
        avc_style_uni_pred_luma_if_function_ptr_array[asm_type][8](
            searchRegionBuffer - (ME_FILTER_TAP >> 1) * lumaStride - 1 + lumaStride,
            lumaStride,
            context_ptr->pos_h_buffer[listIndex][0],
            context_ptr->interpolated_stride,
            searchAreaWidthForAsm,
            search_area_height + 1,
            context_ptr->avctemp_buffer,
            EB_FALSE,
            2);
    }

    if (searchAreaWidthForAsm) {
        // Half pel interpolation of the search region using f1 -> pos_j_buffer
        avc_style_uni_pred_luma_if_function_ptr_array[asm_type][8](
            context_ptr->pos_b_buffer[listIndex][0] + context_ptr->interpolated_stride,
            context_ptr->interpolated_stride,
            context_ptr->pos_j_buffer[listIndex][0],
            context_ptr->interpolated_stride,
            searchAreaWidthForAsm,
            search_area_height + 1,
            context_ptr->avctemp_buffer,
            EB_FALSE,
            2);
    }

#else

    // Half pel interpolation of the search region using f1 -> pos_b_buffer
    HorizontalPelInterpolation(
        searchRegionBuffer - (ME_FILTER_TAP >> 1) * lumaStride - (ME_FILTER_TAP >> 1),
        lumaStride,
        search_area_width + 1,
        search_area_height + ME_FILTER_TAP,
        &(MeIFCoeff[F1][0]),
        inputBitDepth,
        context_ptr->interpolated_stride,
        context_ptr->pos_b_buffer);


    // Half pel interpolation of the search region using f1 -> pos_h_buffer
    VerticalPelInterpolation(
        searchRegionBuffer - (ME_FILTER_TAP >> 1) * lumaStride - 1,
        lumaStride,
        search_area_width + 2,
        search_area_height + 1,
        &(MeIFCoeff[F1][0]),
        inputBitDepth,
        context_ptr->interpolated_stride,
        context_ptr->pos_h_buffer);


    // Half pel interpolation of the search region using f1 -> pos_j_buffer
    VerticalPelInterpolation(
        context_ptr->pos_b_buffer,
        context_ptr->interpolated_stride,
        search_area_width + 1,
        search_area_height + 1,
        &(MeIFCoeff[F1][0]),
        inputBitDepth,
        context_ptr->interpolated_stride,
        context_ptr->pos_j_buffer);

#endif


    return;
}


/*******************************************
* PU_HalfPelRefinement
*   performs Half Pel refinement for one PU
*******************************************/
static void PU_HalfPelRefinement(
    SequenceControlSet_t    *sequence_control_set_ptr,             // input parameter, Sequence control set Ptr
    MeContext_t             *context_ptr,                        // input parameter, ME context Ptr, used to get SB Ptr
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
    uint8_t                   *refBuffer,
    uint32_t                   refStride,
    uint32_t                  *pBestSsd,
#endif
    uint32_t                   puLcuBufferIndex,                  // input parameter, PU origin, used to point to source samples
    uint8_t                   *pos_b_buffer,                        // input parameter, position "b" interpolated search area Ptr
    uint8_t                   *pos_h_buffer,                        // input parameter, position "h" interpolated search area Ptr
    uint8_t                   *pos_j_buffer,                        // input parameter, position "j" interpolated search area Ptr
    uint32_t                   pu_width,                           // input parameter, PU width
    uint32_t                   pu_height,                          // input parameter, PU height
    int16_t                   x_search_area_origin,                 // input parameter, search area origin in the horizontal direction, used to point to reference samples
    int16_t                   y_search_area_origin,                 // input parameter, search area origin in the vertical direction, used to point to reference samples
    EbAsm                   asm_type,
    uint32_t                  *pBestSad,
    uint32_t                  *pBestMV,
    uint8_t                   *psubPelDirection
)
{

    EncodeContext_t         *encode_context_ptr = sequence_control_set_ptr->encode_context_ptr;

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

    //TODO : remove these, and update the MV by just shifts

    xMvHalf[0] = x_mv - 2; // L  position
    xMvHalf[1] = x_mv + 2; // R  position
    xMvHalf[2] = x_mv;     // T  position
    xMvHalf[3] = x_mv;     // B  position
    xMvHalf[4] = x_mv - 2; // TL position
    xMvHalf[5] = x_mv + 2; // TR position
    xMvHalf[6] = x_mv + 2; // BR position
    xMvHalf[7] = x_mv - 2; // BL position

    yMvHalf[0] = y_mv;     // L  position
    yMvHalf[1] = y_mv;     // R  position
    yMvHalf[2] = y_mv - 2; // T  position
    yMvHalf[3] = y_mv + 2; // B  position
    yMvHalf[4] = y_mv - 2; // TL position
    yMvHalf[5] = y_mv - 2; // TR position
    yMvHalf[6] = y_mv + 2; // BR position
    yMvHalf[7] = y_mv + 2; // BL position

#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
    // Compute SSD for the best full search candidate
    if (context_ptr->fractionalSearchMethod == SSD_SEARCH) {
        *pBestSsd = (uint32_t)SpatialFullDistortionKernel_funcPtrArray[asm_type][Log2f(pu_width) - 2](
            &(context_ptr->sb_src_ptr[puLcuBufferIndex]),
            context_ptr->sb_src_stride,
            &(refBuffer[ySearchIndex * refStride + xSearchIndex]),
            refStride,
            pu_width,
            pu_height);
    }
#endif
    // Use SATD only when QP mod, and RC are OFF
    // QP mod, and RC assume that ME distotion is always SAD.
    // This problem might be solved by computing SAD for the best position after fractional search is done, or by considring the full pel resolution SAD.
    {
        // L position
        searchRegionIndex = xSearchIndex + (int16_t)context_ptr->interpolated_stride * ySearchIndex;
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
        distortionLeftPosition = (context_ptr->fractionalSearchMethod == SSD_SEARCH) ?
            SpatialFullDistortionKernel_funcPtrArray[asm_type][Log2f(pu_width) - 2](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride, &(pos_b_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_width, pu_height) :
            (context_ptr->fractionalSearchMethod == SUB_SAD_SEARCH) ?
            (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride << 1, &(pos_b_buffer[searchRegionIndex]), context_ptr->interpolated_stride << 1, pu_height >> 1, pu_width)) << 1 :
            NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride, &(pos_b_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_height, pu_width);
#elif M0_SAD_HALF_QUARTER_PEL_BIPRED_SEARCH
        distortionLeftPosition = (context_ptr->useSubSadFracBipredSearch) ?
            (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride << 1, &(pos_b_buffer[searchRegionIndex]), context_ptr->interpolated_stride << 1, pu_height >> 1, pu_width)) << 1 :
            NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride, &(pos_b_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_height, pu_width);
#else
        distortionLeftPosition = (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride << 1, &(pos_b_buffer[searchRegionIndex]), context_ptr->interpolated_stride << 1, pu_height >> 1, pu_width)) << 1;
#endif
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
        if (context_ptr->fractionalSearchMethod == SSD_SEARCH) {
            if (distortionLeftPosition < *pBestSsd) {
                *pBestSad = (uint32_t)NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride, &(pos_b_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_height, pu_width);
                *pBestMV = ((uint16_t)yMvHalf[0] << 16) | ((uint16_t)xMvHalf[0]);
                *pBestSsd = (uint32_t)distortionLeftPosition;
            }
        }
        else {
            if (distortionLeftPosition < *pBestSad) {
                *pBestSad = (uint32_t)distortionLeftPosition;
                *pBestMV = ((uint16_t)yMvHalf[0] << 16) | ((uint16_t)xMvHalf[0]);
            }
        }
#else
        if (distortionLeftPosition < *pBestSad) {
            *pBestSad = (uint32_t)distortionLeftPosition;
            *pBestMV = ((uint16_t)yMvHalf[0] << 16) | ((uint16_t)xMvHalf[0]);
        }
#endif
        // R position
        searchRegionIndex++;
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
        distortionRightPosition = (context_ptr->fractionalSearchMethod == SSD_SEARCH) ?
            SpatialFullDistortionKernel_funcPtrArray[asm_type][Log2f(pu_width) - 2](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride, &(pos_b_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_width, pu_height) :
            (context_ptr->fractionalSearchMethod == SUB_SAD_SEARCH) ?
            (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride << 1, &(pos_b_buffer[searchRegionIndex]), context_ptr->interpolated_stride << 1, pu_height >> 1, pu_width)) << 1 :
            NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride, &(pos_b_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_height, pu_width);
#elif M0_SAD_HALF_QUARTER_PEL_BIPRED_SEARCH
        distortionRightPosition = (context_ptr->useSubSadFracBipredSearch) ?
            (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride << 1, &(pos_b_buffer[searchRegionIndex]), context_ptr->interpolated_stride << 1, pu_height >> 1, pu_width)) << 1 :
            NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride, &(pos_b_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_height, pu_width);
#else
        distortionRightPosition = (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride << 1, &(pos_b_buffer[searchRegionIndex]), context_ptr->interpolated_stride << 1, pu_height >> 1, pu_width)) << 1;
#endif
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
        if (context_ptr->fractionalSearchMethod == SSD_SEARCH) {
            if (distortionRightPosition < *pBestSsd) {
                *pBestSad = (uint32_t)NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride, &(pos_b_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_height, pu_width);
                *pBestMV = ((uint16_t)yMvHalf[1] << 16) | ((uint16_t)xMvHalf[1]);
                *pBestSsd = (uint32_t)distortionRightPosition;
            }
        }
        else {
            if (distortionRightPosition < *pBestSad) {
                *pBestSad = (uint32_t)distortionRightPosition;
                *pBestMV = ((uint16_t)yMvHalf[1] << 16) | ((uint16_t)xMvHalf[1]);
            }
        }
#else
        if (distortionRightPosition < *pBestSad) {
            *pBestSad = (uint32_t)distortionRightPosition;
            *pBestMV = ((uint16_t)yMvHalf[1] << 16) | ((uint16_t)xMvHalf[1]);
        }
#endif
        // T position
        searchRegionIndex = xSearchIndex + (int16_t)context_ptr->interpolated_stride * ySearchIndex;
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
        distortionTopPosition = (context_ptr->fractionalSearchMethod == SSD_SEARCH) ?
            SpatialFullDistortionKernel_funcPtrArray[asm_type][Log2f(pu_width) - 2](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride, &(pos_h_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_width, pu_height) :
            (context_ptr->fractionalSearchMethod == SUB_SAD_SEARCH) ?
            (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride << 1, &(pos_h_buffer[searchRegionIndex]), context_ptr->interpolated_stride << 1, pu_height >> 1, pu_width)) << 1 :
            NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride, &(pos_h_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_height, pu_width);
#elif M0_SAD_HALF_QUARTER_PEL_BIPRED_SEARCH
        distortionTopPosition = (context_ptr->useSubSadFracBipredSearch) ?
            (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride << 1, &(pos_h_buffer[searchRegionIndex]), context_ptr->interpolated_stride << 1, pu_height >> 1, pu_width)) << 1 :
            NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride, &(pos_h_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_height, pu_width);
#else
        distortionTopPosition = (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride << 1, &(pos_h_buffer[searchRegionIndex]), context_ptr->interpolated_stride << 1, pu_height >> 1, pu_width)) << 1;
#endif
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
        if (context_ptr->fractionalSearchMethod == SSD_SEARCH) {
            if (distortionTopPosition < *pBestSsd) {
                *pBestSad = (uint32_t)NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride, &(pos_h_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_height, pu_width);
                *pBestMV = ((uint16_t)yMvHalf[2] << 16) | ((uint16_t)xMvHalf[2]);
                *pBestSsd = (uint32_t)distortionTopPosition;
            }
        }
        else {
            if (distortionTopPosition < *pBestSad) {
                *pBestSad = (uint32_t)distortionTopPosition;
                *pBestMV = ((uint16_t)yMvHalf[2] << 16) | ((uint16_t)xMvHalf[2]);
            }
        }
#else
        if (distortionTopPosition < *pBestSad) {
            *pBestSad = (uint32_t)distortionTopPosition;
            *pBestMV = ((uint16_t)yMvHalf[2] << 16) | ((uint16_t)xMvHalf[2]);
        }
#endif

        // B position
        searchRegionIndex += (int16_t)context_ptr->interpolated_stride;
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
        distortionBottomPosition = (context_ptr->fractionalSearchMethod == SSD_SEARCH) ?
            SpatialFullDistortionKernel_funcPtrArray[asm_type][Log2f(pu_width) - 2](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride, &(pos_h_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_width, pu_height) :
            (context_ptr->fractionalSearchMethod == SUB_SAD_SEARCH) ?
            (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride << 1, &(pos_h_buffer[searchRegionIndex]), context_ptr->interpolated_stride << 1, pu_height >> 1, pu_width)) << 1 :
            NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride, &(pos_h_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_height, pu_width);
#elif M0_SAD_HALF_QUARTER_PEL_BIPRED_SEARCH
        distortionBottomPosition = (context_ptr->useSubSadFracBipredSearch) ?
            (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride << 1, &(pos_h_buffer[searchRegionIndex]), context_ptr->interpolated_stride << 1, pu_height >> 1, pu_width)) << 1 :
            NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride, &(pos_h_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_height, pu_width);
#else
        distortionBottomPosition = (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride << 1, &(pos_h_buffer[searchRegionIndex]), context_ptr->interpolated_stride << 1, pu_height >> 1, pu_width)) << 1;
#endif
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
        if (context_ptr->fractionalSearchMethod == SSD_SEARCH) {
            if (distortionBottomPosition < *pBestSsd) {
                *pBestSad = (uint32_t)NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride, &(pos_h_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_height, pu_width);
                *pBestMV = ((uint16_t)yMvHalf[3] << 16) | ((uint16_t)xMvHalf[3]);
                *pBestSsd = (uint32_t)distortionBottomPosition;
            }
        }
        else {
            if (distortionBottomPosition < *pBestSad) {
                *pBestSad = (uint32_t)distortionBottomPosition;
                *pBestMV = ((uint16_t)yMvHalf[3] << 16) | ((uint16_t)xMvHalf[3]);
            }
        }
#else
        if (distortionBottomPosition < *pBestSad) {
            *pBestSad = (uint32_t)distortionBottomPosition;
            *pBestMV = ((uint16_t)yMvHalf[3] << 16) | ((uint16_t)xMvHalf[3]);
        }
#endif

        //TL position
        searchRegionIndex = xSearchIndex + (int16_t)context_ptr->interpolated_stride * ySearchIndex;
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
        distortionTopLeftPosition = (context_ptr->fractionalSearchMethod == SSD_SEARCH) ?
            SpatialFullDistortionKernel_funcPtrArray[asm_type][Log2f(pu_width) - 2](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_width, pu_height) :
            (context_ptr->fractionalSearchMethod == SUB_SAD_SEARCH) ?
            (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride << 1, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride << 1, pu_height >> 1, pu_width)) << 1 :
            NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_height, pu_width);
#elif M0_SAD_HALF_QUARTER_PEL_BIPRED_SEARCH
        distortionTopLeftPosition = (context_ptr->useSubSadFracBipredSearch) ?
            (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride << 1, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride << 1, pu_height >> 1, pu_width)) << 1 :
            NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_height, pu_width);
#else
        distortionTopLeftPosition = (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride << 1, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride << 1, pu_height >> 1, pu_width)) << 1;
#endif
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
        if (context_ptr->fractionalSearchMethod == SSD_SEARCH) {
            if (distortionTopLeftPosition < *pBestSsd) {
                *pBestSad = (uint32_t)NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_height, pu_width);
                *pBestMV = ((uint16_t)yMvHalf[4] << 16) | ((uint16_t)xMvHalf[4]);
                *pBestSsd = (uint32_t)distortionTopLeftPosition;
            }
        }
        else {
            if (distortionTopLeftPosition < *pBestSad) {
                *pBestSad = (uint32_t)distortionTopLeftPosition;
                *pBestMV = ((uint16_t)yMvHalf[4] << 16) | ((uint16_t)xMvHalf[4]);
            }
        }
#else
        if (distortionTopLeftPosition < *pBestSad) {
            *pBestSad = (uint32_t)distortionTopLeftPosition;
            *pBestMV = ((uint16_t)yMvHalf[4] << 16) | ((uint16_t)xMvHalf[4]);
        }
#endif

        //TR position
        searchRegionIndex++;
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
        distortionTopRightPosition = (context_ptr->fractionalSearchMethod == SSD_SEARCH) ?
            SpatialFullDistortionKernel_funcPtrArray[asm_type][Log2f(pu_width) - 2](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_width, pu_height) :
            (context_ptr->fractionalSearchMethod == SUB_SAD_SEARCH) ?
            (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride << 1, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride << 1, pu_height >> 1, pu_width)) << 1 :
            NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_height, pu_width);
#elif M0_SAD_HALF_QUARTER_PEL_BIPRED_SEARCH
        distortionTopRightPosition = (context_ptr->useSubSadFracBipredSearch) ?
            (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride << 1, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride << 1, pu_height >> 1, pu_width)) << 1 :
            NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_height, pu_width);
#else
        distortionTopRightPosition = (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride << 1, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride << 1, pu_height >> 1, pu_width)) << 1;
#endif
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
        if (context_ptr->fractionalSearchMethod == SSD_SEARCH) {
            if (distortionTopRightPosition < *pBestSsd) {
                *pBestSad = (uint32_t)NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_height, pu_width);
                *pBestMV = ((uint16_t)yMvHalf[5] << 16) | ((uint16_t)xMvHalf[5]);
                *pBestSsd = (uint32_t)distortionTopRightPosition;
            }
        }
        else {
            if (distortionTopRightPosition < *pBestSad) {
                *pBestSad = (uint32_t)distortionTopRightPosition;
                *pBestMV = ((uint16_t)yMvHalf[5] << 16) | ((uint16_t)xMvHalf[5]);
            }
        }
#else
        if (distortionTopRightPosition < *pBestSad) {
            *pBestSad = (uint32_t)distortionTopRightPosition;
            *pBestMV = ((uint16_t)yMvHalf[5] << 16) | ((uint16_t)xMvHalf[5]);
        }
#endif

        //BR position
        searchRegionIndex += (int16_t)context_ptr->interpolated_stride;
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
        distortionBottomRightPosition = (context_ptr->fractionalSearchMethod == SSD_SEARCH) ?
            SpatialFullDistortionKernel_funcPtrArray[asm_type][Log2f(pu_width) - 2](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_width, pu_height) :
            (context_ptr->fractionalSearchMethod == SUB_SAD_SEARCH) ?
            (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride << 1, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride << 1, pu_height >> 1, pu_width)) << 1 :
            NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_height, pu_width);
#elif M0_SAD_HALF_QUARTER_PEL_BIPRED_SEARCH
        distortionBottomRightPosition = (context_ptr->useSubSadFracBipredSearch) ?
            (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride << 1, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride << 1, pu_height >> 1, pu_width)) << 1 :
            NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_height, pu_width);
#else
        distortionBottomRightPosition = (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride << 1, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride << 1, pu_height >> 1, pu_width)) << 1;
#endif
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
        if (context_ptr->fractionalSearchMethod == SSD_SEARCH) {
            if (distortionBottomRightPosition < *pBestSsd) {
                *pBestSad = (uint32_t)NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_height, pu_width);
                *pBestMV = ((uint16_t)yMvHalf[6] << 16) | ((uint16_t)xMvHalf[6]);
                *pBestSsd = (uint32_t)distortionBottomRightPosition;
            }
        }
        else {
            if (distortionBottomRightPosition < *pBestSad) {
                *pBestSad = (uint32_t)distortionBottomRightPosition;
                *pBestMV = ((uint16_t)yMvHalf[6] << 16) | ((uint16_t)xMvHalf[6]);
            }
        }
#else
        if (distortionBottomRightPosition < *pBestSad) {
            *pBestSad = (uint32_t)distortionBottomRightPosition;
            *pBestMV = ((uint16_t)yMvHalf[6] << 16) | ((uint16_t)xMvHalf[6]);
        }
#endif

        //BL position
        searchRegionIndex--;
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
        distortionBottomLeftPosition = (context_ptr->fractionalSearchMethod == SSD_SEARCH) ?
            SpatialFullDistortionKernel_funcPtrArray[asm_type][Log2f(pu_width) - 2](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_width, pu_height) :
            (context_ptr->fractionalSearchMethod == SUB_SAD_SEARCH) ?
            (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride << 1, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride << 1, pu_height >> 1, pu_width)) << 1 :
            (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_height, pu_width));
#elif M0_SAD_HALF_QUARTER_PEL_BIPRED_SEARCH
        distortionBottomLeftPosition = (context_ptr->useSubSadFracBipredSearch) ?
            (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride << 1, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride << 1, pu_height >> 1, pu_width)) << 1 :
            (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_height, pu_width));
#else
        distortionBottomLeftPosition = (NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride << 1, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride << 1, pu_height >> 1, pu_width)) << 1;
#endif
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
        if (context_ptr->fractionalSearchMethod == SSD_SEARCH) {
            if (distortionBottomLeftPosition < *pBestSsd) {
                *pBestSad = (uint32_t)(NxMSadKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_src_ptr[puLcuBufferIndex]), context_ptr->sb_src_stride, &(pos_j_buffer[searchRegionIndex]), context_ptr->interpolated_stride, pu_height, pu_width));
                *pBestMV = ((uint16_t)yMvHalf[7] << 16) | ((uint16_t)xMvHalf[7]);
                *pBestSsd = (uint32_t)distortionBottomLeftPosition;
            }
        }
        else {
            if (distortionBottomLeftPosition < *pBestSad) {
                *pBestSad = (uint32_t)distortionBottomLeftPosition;
                *pBestMV = ((uint16_t)yMvHalf[7] << 16) | ((uint16_t)xMvHalf[7]);
            }
        }
#else
        if (distortionBottomLeftPosition < *pBestSad) {
            *pBestSad = (uint32_t)distortionBottomLeftPosition;
            *pBestMV = ((uint16_t)yMvHalf[7] << 16) | ((uint16_t)xMvHalf[7]);
        }
#endif
    }

    bestHalfSad = MIN(distortionLeftPosition, MIN(distortionRightPosition, MIN(distortionTopPosition, MIN(distortionBottomPosition, MIN(distortionTopLeftPosition, MIN(distortionTopRightPosition, MIN(distortionBottomLeftPosition, distortionBottomRightPosition)))))));


    if (bestHalfSad == distortionLeftPosition) {
        *psubPelDirection = LEFT_POSITION;
    }
    else if (bestHalfSad == distortionRightPosition) {
        *psubPelDirection = RIGHT_POSITION;
    }
    else if (bestHalfSad == distortionTopPosition) {
        *psubPelDirection = TOP_POSITION;
    }
    else if (bestHalfSad == distortionBottomPosition) {
        *psubPelDirection = BOTTOM_POSITION;
    }
    else if (bestHalfSad == distortionTopLeftPosition) {
        *psubPelDirection = TOP_LEFT_POSITION;
    }
    else if (bestHalfSad == distortionTopRightPosition) {
        *psubPelDirection = TOP_RIGHT_POSITION;
    }
    else if (bestHalfSad == distortionBottomLeftPosition) {
        *psubPelDirection = BOTTOM_LEFT_POSITION;
    }
    else if (bestHalfSad == distortionBottomRightPosition) {
        *psubPelDirection = BOTTOM_RIGHT_POSITION;
    }
    return;
}

/*******************************************
* HalfPelSearch_LCU
*   performs Half Pel refinement for the 85 PUs
*******************************************/
void HalfPelSearch_LCU(
    SequenceControlSet_t    *sequence_control_set_ptr,             // input parameter, Sequence control set Ptr
#if DISABLE_NSQ_FOR_NON_REF || DISABLE_NSQ
    PictureParentControlSet_t *picture_control_set_ptr,
#endif
    MeContext_t             *context_ptr,                        // input/output parameter, ME context Ptr, used to get/update ME results
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
    uint8_t                   *refBuffer,
    uint32_t                   refStride,
#endif
    uint8_t                   *pos_b_buffer,                        // input parameter, position "b" interpolated search area Ptr
    uint8_t                   *pos_h_buffer,                        // input parameter, position "h" interpolated search area Ptr
    uint8_t                   *pos_j_buffer,                        // input parameter, position "j" interpolated search area Ptr
    int16_t                   x_search_area_origin,                 // input parameter, search area origin in the horizontal direction, used to point to reference samples
    int16_t                   y_search_area_origin,                 // input parameter, search area origin in the vertical direction, used to point to reference samples
    EbAsm                   asm_type,
    EbBool                     disable8x8CuInMeFlag,
    EbBool                    enableHalfPel32x32,
    EbBool                    enableHalfPel16x16,
    EbBool                    enableHalfPel8x8) {

    uint32_t idx;
    uint32_t pu_index;
    uint32_t puShiftXIndex;
    uint32_t puShiftYIndex;
    uint32_t puLcuBufferIndex;
    uint32_t posbBufferIndex;
    uint32_t poshBufferIndex;
    uint32_t posjBufferIndex;

#if M0_64x64_32x32_HALF_QUARTER_PEL
    if (context_ptr->fractional_search64x64)
        PU_HalfPelRefinement(
            sequence_control_set_ptr,
            context_ptr,
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
            &(refBuffer[0]),
            refStride,
            context_ptr->p_best_ssd64x64,
#endif
            0,
            &(pos_b_buffer[0]),
            &(pos_h_buffer[0]),
            &(pos_j_buffer[0]),
            64,
            64,
            x_search_area_origin,
            y_search_area_origin,
            asm_type,
            context_ptr->p_best_sad64x64,
            context_ptr->p_best_mv64x64,
            &context_ptr->psub_pel_direction64x64);
#else
    //no PU 64x64, Half Pel Refinement
#endif
    if (enableHalfPel32x32)
    {
        // 32x32 [4 partitions]
        for (pu_index = 0; pu_index < 4; ++pu_index) {

            puShiftXIndex = (pu_index & 0x01) << 5;
            puShiftYIndex = (pu_index >> 1) << 5;

            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->sb_src_stride;
            posbBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;
            poshBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;
            posjBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;


            PU_HalfPelRefinement(
                sequence_control_set_ptr,
                context_ptr,
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
                &(refBuffer[puShiftYIndex * refStride + puShiftXIndex]),
                refStride,
                &context_ptr->p_best_ssd32x32[pu_index],
#endif
                puLcuBufferIndex,
                &(pos_b_buffer[posbBufferIndex]),
                &(pos_h_buffer[poshBufferIndex]),
                &(pos_j_buffer[posjBufferIndex]),
                32,
                32,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad32x32[pu_index],
                &context_ptr->p_best_mv32x32[pu_index],
                &context_ptr->psub_pel_direction32x32[pu_index]);
        }
    }
    if (enableHalfPel16x16)
    {
        // 16x16 [16 partitions]
        for (pu_index = 0; pu_index < 16; ++pu_index) {

            idx = tab16x16[pu_index];

            puShiftXIndex = (pu_index & 0x03) << 4;
            puShiftYIndex = (pu_index >> 2) << 4;

            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->sb_src_stride;
            posbBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;
            poshBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;
            posjBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;

            PU_HalfPelRefinement(
                sequence_control_set_ptr,
                context_ptr,
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
                &(refBuffer[puShiftYIndex * refStride + puShiftXIndex]),
                refStride,
                &context_ptr->p_best_ssd16x16[idx],
#endif
                puLcuBufferIndex,
                &(pos_b_buffer[posbBufferIndex]),
                &(pos_h_buffer[poshBufferIndex]),
                &(pos_j_buffer[posjBufferIndex]),
                16,
                16,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad16x16[idx],
                &context_ptr->p_best_mv16x16[idx],
                &context_ptr->psub_pel_direction16x16[idx]);

        }
    }
    if (enableHalfPel8x8)
    {
        // 8x8   [64 partitions]
        if (!disable8x8CuInMeFlag) {
            for (pu_index = 0; pu_index < 64; ++pu_index) {

                idx = tab8x8[pu_index];  //TODO bitwise this

                puShiftXIndex = (pu_index & 0x07) << 3;
                puShiftYIndex = (pu_index >> 3) << 3;

                puLcuBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->sb_src_stride;

                posbBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;
                poshBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;
                posjBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;

                PU_HalfPelRefinement(
                    sequence_control_set_ptr,
                    context_ptr,
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
                    &(refBuffer[puShiftYIndex * refStride + puShiftXIndex]),
                    refStride,
                    &context_ptr->p_best_ssd8x8[idx],
#endif
                    puLcuBufferIndex,
                    &(pos_b_buffer[posbBufferIndex]),
                    &(pos_h_buffer[poshBufferIndex]),
                    &(pos_j_buffer[posjBufferIndex]),
                    8,
                    8,
                    x_search_area_origin,
                    y_search_area_origin,
                    asm_type,
                    &context_ptr->p_best_sad8x8[idx],
                    &context_ptr->p_best_mv8x8[idx],
                    &context_ptr->psub_pel_direction8x8[idx]);

            }
        }
    }
#if DISABLE_NSQ_FOR_NON_REF || DISABLE_NSQ
#if ENCODER_MODE_CLEANUP
    if (picture_control_set_ptr->pic_depth_mode <= PIC_ALL_C_DEPTH_MODE) {
#else
    if (picture_control_set_ptr->non_square_block_flag) {
#endif
#else
    if (sequence_control_set_ptr->static_config.ext_block_flag) {
#endif

        // 64x32
        for (pu_index = 0; pu_index < 2; ++pu_index) {


            puShiftXIndex = 0;
            puShiftYIndex = pu_index << 5;

            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->sb_src_stride;

            posbBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;
            poshBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;
            posjBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;

            PU_HalfPelRefinement(
                sequence_control_set_ptr,
                context_ptr,
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
                &(refBuffer[puShiftYIndex * refStride + puShiftXIndex]),
                refStride,
                &context_ptr->p_best_ssd64x32[pu_index],
#endif
                puLcuBufferIndex,
                &(pos_b_buffer[posbBufferIndex]),
                &(pos_h_buffer[poshBufferIndex]),
                &(pos_j_buffer[posjBufferIndex]),
                64,
                32,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad64x32[pu_index],
                &context_ptr->p_best_mv64x32[pu_index],
                &context_ptr->psub_pel_direction64x32[pu_index]);

        }

        // 32x16
        for (pu_index = 0; pu_index < 8; ++pu_index) {

            idx = tab32x16[pu_index];  //TODO bitwise this

            puShiftXIndex = (pu_index & 0x01) << 5;
            puShiftYIndex = (pu_index >> 1) << 4;


            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->sb_src_stride;

            posbBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;
            poshBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;
            posjBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;

            PU_HalfPelRefinement(
                sequence_control_set_ptr,
                context_ptr,
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
                &(refBuffer[puShiftYIndex * refStride + puShiftXIndex]),
                refStride,
                &context_ptr->p_best_ssd32x16[idx],
#endif
                puLcuBufferIndex,
                &(pos_b_buffer[posbBufferIndex]),
                &(pos_h_buffer[poshBufferIndex]),
                &(pos_j_buffer[posjBufferIndex]),
                32,
                16,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad32x16[idx],
                &context_ptr->p_best_mv32x16[idx],
                &context_ptr->psub_pel_direction32x16[idx]);

        }

        // 16x8
        for (pu_index = 0; pu_index < 32; ++pu_index) {


            idx = tab16x8[pu_index];

            puShiftXIndex = (pu_index & 0x03) << 4;
            puShiftYIndex = (pu_index >> 2) << 3;

            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->sb_src_stride;

            posbBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;
            poshBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;
            posjBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;

            PU_HalfPelRefinement(
                sequence_control_set_ptr,
                context_ptr,
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
                &(refBuffer[puShiftYIndex * refStride + puShiftXIndex]),
                refStride,
                &context_ptr->p_best_ssd16x8[idx],
#endif
                puLcuBufferIndex,
                &(pos_b_buffer[posbBufferIndex]),
                &(pos_h_buffer[poshBufferIndex]),
                &(pos_j_buffer[posjBufferIndex]),
                16,
                8,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad16x8[idx],
                &context_ptr->p_best_mv16x8[idx],
                &context_ptr->psub_pel_direction16x8[idx]);

        }

        // 32x64
        for (pu_index = 0; pu_index < 2; ++pu_index) {

            puShiftXIndex = pu_index << 5;
            puShiftYIndex = 0;

            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->sb_src_stride;

            posbBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;
            poshBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;
            posjBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;

            PU_HalfPelRefinement(
                sequence_control_set_ptr,
                context_ptr,
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
                &(refBuffer[puShiftYIndex * refStride + puShiftXIndex]),
                refStride,
                &context_ptr->p_best_ssd32x64[pu_index],
#endif
                puLcuBufferIndex,
                &(pos_b_buffer[posbBufferIndex]),
                &(pos_h_buffer[poshBufferIndex]),
                &(pos_j_buffer[posjBufferIndex]),
                32,
                64,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad32x64[pu_index],
                &context_ptr->p_best_mv32x64[pu_index],
                &context_ptr->psub_pel_direction32x64[pu_index]);

        }

        // 16x32
        for (pu_index = 0; pu_index < 8; ++pu_index) {

            idx = tab16x32[pu_index];

            puShiftXIndex = (pu_index & 0x03) << 4;
            puShiftYIndex = (pu_index >> 2) << 5;

            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->sb_src_stride;

            posbBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;
            poshBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;
            posjBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;

            PU_HalfPelRefinement(
                sequence_control_set_ptr,
                context_ptr,
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
                &(refBuffer[puShiftYIndex * refStride + puShiftXIndex]),
                refStride,
                &context_ptr->p_best_ssd16x32[idx],
#endif
                puLcuBufferIndex,
                &(pos_b_buffer[posbBufferIndex]),
                &(pos_h_buffer[poshBufferIndex]),
                &(pos_j_buffer[posjBufferIndex]),
                16,
                32,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad16x32[idx],
                &context_ptr->p_best_mv16x32[idx],
                &context_ptr->psub_pel_direction16x32[idx]);
        }

        // 8x16
        for (pu_index = 0; pu_index < 32; ++pu_index) {

            idx = tab8x16[pu_index];

            puShiftXIndex = (pu_index & 0x07) << 3;
            puShiftYIndex = (pu_index >> 3) << 4;

            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->sb_src_stride;

            posbBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;
            poshBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;
            posjBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;

            PU_HalfPelRefinement(
                sequence_control_set_ptr,
                context_ptr,
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
                &(refBuffer[puShiftYIndex * refStride + puShiftXIndex]),
                refStride,
                &context_ptr->p_best_ssd8x16[idx],
#endif
                puLcuBufferIndex,
                &(pos_b_buffer[posbBufferIndex]),
                &(pos_h_buffer[poshBufferIndex]),
                &(pos_j_buffer[posjBufferIndex]),
                8,
                16,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad8x16[idx],
                &context_ptr->p_best_mv8x16[idx],
                &context_ptr->psub_pel_direction8x16[idx]);
        }

        // 32x8
        for (pu_index = 0; pu_index < 16; ++pu_index) {

            idx = tab32x8[pu_index];

            puShiftXIndex = (pu_index & 0x01) << 5;
            puShiftYIndex = (pu_index >> 1) << 3;

            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->sb_src_stride;

            posbBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;
            poshBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;
            posjBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;

            PU_HalfPelRefinement(
                sequence_control_set_ptr,
                context_ptr,
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
                &(refBuffer[puShiftYIndex * refStride + puShiftXIndex]),
                refStride,
                &context_ptr->p_best_ssd32x8[idx],
#endif
                puLcuBufferIndex,
                &(pos_b_buffer[posbBufferIndex]),
                &(pos_h_buffer[poshBufferIndex]),
                &(pos_j_buffer[posjBufferIndex]),
                32,
                8,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad32x8[idx],
                &context_ptr->p_best_mv32x8[idx],
                &context_ptr->psub_pel_direction32x8[idx]);
        }

        for (pu_index = 0; pu_index < 16; ++pu_index) {

            idx = tab8x32[pu_index];

            puShiftXIndex = (pu_index & 0x07) << 3;
            puShiftYIndex = (pu_index >> 3) << 5;

            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->sb_src_stride;
            posbBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;
            poshBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;
            posjBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;

            PU_HalfPelRefinement(
                sequence_control_set_ptr,
                context_ptr,
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
                &(refBuffer[puShiftYIndex * refStride + puShiftXIndex]),
                refStride,
                &context_ptr->p_best_ssd8x32[idx],
#endif
                puLcuBufferIndex,
                &(pos_b_buffer[posbBufferIndex]),
                &(pos_h_buffer[poshBufferIndex]),
                &(pos_j_buffer[posjBufferIndex]),
                8,
                32,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad8x32[idx],
                &context_ptr->p_best_mv8x32[idx],
                &context_ptr->psub_pel_direction8x32[idx]);
        }

        for (pu_index = 0; pu_index < 4; ++pu_index) {

            idx = pu_index;

            puShiftXIndex = 0;
            puShiftYIndex = pu_index << 4;

            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->sb_src_stride;
            posbBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;
            poshBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;
            posjBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;

            PU_HalfPelRefinement(
                sequence_control_set_ptr,
                context_ptr,
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
                &(refBuffer[puShiftYIndex * refStride + puShiftXIndex]),
                refStride,
                &context_ptr->p_best_ssd64x16[idx],
#endif
                puLcuBufferIndex,
                &(pos_b_buffer[posbBufferIndex]),
                &(pos_h_buffer[poshBufferIndex]),
                &(pos_j_buffer[posjBufferIndex]),
                64,
                16,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad64x16[idx],
                &context_ptr->p_best_mv64x16[idx],
                &context_ptr->psub_pel_direction64x16[idx]);
        }

        for (pu_index = 0; pu_index < 4; ++pu_index) {

            idx = pu_index;

            puShiftXIndex = pu_index << 4;
            puShiftYIndex = 0;

            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->sb_src_stride;
            posbBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;
            poshBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;
            posjBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->interpolated_stride;

            PU_HalfPelRefinement(
                sequence_control_set_ptr,
                context_ptr,
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
                &(refBuffer[puShiftYIndex * refStride + puShiftXIndex]),
                refStride,
                &context_ptr->p_best_ssd16x64[idx],
#endif
                puLcuBufferIndex,
                &(pos_b_buffer[posbBufferIndex]),
                &(pos_h_buffer[poshBufferIndex]),
                &(pos_j_buffer[posjBufferIndex]),
                16,
                64,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad16x64[idx],
                &context_ptr->p_best_mv16x64[idx],
                &context_ptr->psub_pel_direction16x64[idx]);
        }

    }

    return;
}
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
/*******************************************
* combined_averaging_ssd
*
*******************************************/
uint32_t combined_averaging_ssd_c(
    uint8_t   *src,
    ptrdiff_t  src_stride,
    uint8_t   *ref1,
    ptrdiff_t  ref1_stride,
    uint8_t   *ref2,
    ptrdiff_t  ref2_stride,
    uint32_t   height,
    uint32_t   width)
{
    uint32_t x, y;
    uint32_t ssd = 0;
    uint8_t avgpel;
    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            avgpel = (ref1[x] + ref2[x] + 1) >> 1;
            ssd += SQR((int64_t)(src[x]) - (avgpel));
        }
        src += src_stride;
        ref1 += ref1_stride;
        ref2 += ref2_stride;
    }
    return ssd;
}
#endif
#if M0_ME_QUARTER_PEL_SEARCH
/*******************************************
* PU_QuarterPelRefinementOnTheFly
*   performs Quarter Pel refinement for each PU
*******************************************/
static void PU_QuarterPelRefinementOnTheFly(
    MeContext_t           *context_ptr,                      // [IN] ME context Ptr, used to get SB Ptr
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
    uint32_t                *pBestSsd,
#endif
    uint32_t                 puLcuBufferIndex,                // [IN] PU origin, used to point to source samples
    uint8_t                **buf1,                            // [IN]
    uint32_t                *buf1Stride,
    uint8_t                **buf2,                            // [IN]
    uint32_t                *buf2Stride,
    uint32_t                 pu_width,                         // [IN]  PU width
    uint32_t                 pu_height,                        // [IN]  PU height
    int16_t                 x_search_area_origin,               // [IN] search area origin in the horizontal direction, used to point to reference samples
    int16_t                 y_search_area_origin,               // [IN] search area origin in the vertical direction, used to point to reference samples
    EbAsm                 asm_type,
    uint32_t                *pBestSad,
    uint32_t                *pBestMV,
    uint8_t                  sub_pel_direction)
{

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

    if ((y_mv & 2) + ((x_mv & 2) >> 1)) {

        validTL = (EbBool)(sub_pel_direction == RIGHT_POSITION || sub_pel_direction == BOTTOM_RIGHT_POSITION || sub_pel_direction == BOTTOM_POSITION);
        validT = (EbBool)(sub_pel_direction == BOTTOM_RIGHT_POSITION || sub_pel_direction == BOTTOM_POSITION || sub_pel_direction == BOTTOM_LEFT_POSITION);
        validTR = (EbBool)(sub_pel_direction == BOTTOM_POSITION || sub_pel_direction == BOTTOM_LEFT_POSITION || sub_pel_direction == LEFT_POSITION);
        validR = (EbBool)(sub_pel_direction == BOTTOM_LEFT_POSITION || sub_pel_direction == LEFT_POSITION || sub_pel_direction == TOP_LEFT_POSITION);
        validBR = (EbBool)(sub_pel_direction == LEFT_POSITION || sub_pel_direction == TOP_LEFT_POSITION || sub_pel_direction == TOP_POSITION);
        validB = (EbBool)(sub_pel_direction == TOP_LEFT_POSITION || sub_pel_direction == TOP_POSITION || sub_pel_direction == TOP_RIGHT_POSITION);
        validBL = (EbBool)(sub_pel_direction == TOP_POSITION || sub_pel_direction == TOP_RIGHT_POSITION || sub_pel_direction == RIGHT_POSITION);
        validL = (EbBool)(sub_pel_direction == TOP_RIGHT_POSITION || sub_pel_direction == RIGHT_POSITION || sub_pel_direction == BOTTOM_RIGHT_POSITION);

    }
    else {

        validTL = (EbBool)(sub_pel_direction == LEFT_POSITION || sub_pel_direction == TOP_LEFT_POSITION || sub_pel_direction == TOP_POSITION);
        validT = (EbBool)(sub_pel_direction == TOP_LEFT_POSITION || sub_pel_direction == TOP_POSITION || sub_pel_direction == TOP_RIGHT_POSITION);
        validTR = (EbBool)(sub_pel_direction == TOP_POSITION || sub_pel_direction == TOP_RIGHT_POSITION || sub_pel_direction == RIGHT_POSITION);
        validR = (EbBool)(sub_pel_direction == TOP_RIGHT_POSITION || sub_pel_direction == RIGHT_POSITION || sub_pel_direction == BOTTOM_RIGHT_POSITION);
        validBR = (EbBool)(sub_pel_direction == RIGHT_POSITION || sub_pel_direction == BOTTOM_RIGHT_POSITION || sub_pel_direction == BOTTOM_POSITION);
        validB = (EbBool)(sub_pel_direction == BOTTOM_RIGHT_POSITION || sub_pel_direction == BOTTOM_POSITION || sub_pel_direction == BOTTOM_LEFT_POSITION);
        validBL = (EbBool)(sub_pel_direction == BOTTOM_POSITION || sub_pel_direction == BOTTOM_LEFT_POSITION || sub_pel_direction == LEFT_POSITION);
        validL = (EbBool)(sub_pel_direction == BOTTOM_LEFT_POSITION || sub_pel_direction == LEFT_POSITION || sub_pel_direction == TOP_LEFT_POSITION);
    }

    xMvQuarter[0] = x_mv - 1; // L  position
    xMvQuarter[1] = x_mv + 1; // R  position
    xMvQuarter[2] = x_mv;     // T  position
    xMvQuarter[3] = x_mv;     // B  position
    xMvQuarter[4] = x_mv - 1; // TL position
    xMvQuarter[5] = x_mv + 1; // TR position
    xMvQuarter[6] = x_mv + 1; // BR position
    xMvQuarter[7] = x_mv - 1; // BL position

    yMvQuarter[0] = y_mv;     // L  position
    yMvQuarter[1] = y_mv;     // R  position
    yMvQuarter[2] = y_mv - 1; // T  position
    yMvQuarter[3] = y_mv + 1; // B  position
    yMvQuarter[4] = y_mv - 1; // TL position
    yMvQuarter[5] = y_mv - 1; // TR position
    yMvQuarter[6] = y_mv + 1; // BR position
    yMvQuarter[7] = y_mv + 1; // BL position

    // Use SATD only when QP mod, and RC are OFF
    // QP mod, and RC assume that ME distotion is always SAD.
    // This problem might be solved by computing SAD for the best position after fractional search is done, or by considring the full pel resolution SAD.

    {
        // L position
        if (validL) {

            searchRegionIndex1 = (int32_t)xSearchIndex + (int32_t)buf1Stride[0] * (int32_t)ySearchIndex;
            searchRegionIndex2 = (int32_t)xSearchIndex + (int32_t)buf2Stride[0] * (int32_t)ySearchIndex;

#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
            dist = (context_ptr->fractionalSearchMethod == SSD_SEARCH) ?
                combined_averaging_ssd_func_ptr_array[asm_type](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64, buf1[0] + searchRegionIndex1, buf1Stride[0], buf2[0] + searchRegionIndex2, buf2Stride[0], pu_height, pu_width) :
                (context_ptr->fractionalSearchMethod == SUB_SAD_SEARCH) ?
                (NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64 << 1, buf1[0] + searchRegionIndex1, buf1Stride[0] << 1, buf2[0] + searchRegionIndex2, buf2Stride[0] << 1, pu_height >> 1, pu_width)) << 1 :
                NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64, buf1[0] + searchRegionIndex1, buf1Stride[0], buf2[0] + searchRegionIndex2, buf2Stride[0], pu_height, pu_width);
#elif M0_SAD_HALF_QUARTER_PEL_BIPRED_SEARCH
            dist = (context_ptr->useSubSadFracBipredSearch) ?
                (NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64 << 1, buf1[0] + searchRegionIndex1, buf1Stride[0] << 1, buf2[0] + searchRegionIndex2, buf2Stride[0] << 1, pu_height >> 1, pu_width)) << 1 :
                NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64, buf1[0] + searchRegionIndex1, buf1Stride[0], buf2[0] + searchRegionIndex2, buf2Stride[0], pu_height, pu_width);
#else
            dist = NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64 << 1, buf1[0] + searchRegionIndex1, buf1Stride[0] << 1, buf2[0] + searchRegionIndex2, buf2Stride[0] << 1, pu_height >> 1, pu_width);
            dist = dist << 1;
#endif
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
            if (context_ptr->fractionalSearchMethod == SSD_SEARCH) {
                if (dist < *pBestSsd) {
                    *pBestSad = (uint32_t)NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64, buf1[0] + searchRegionIndex1, buf1Stride[0], buf2[0] + searchRegionIndex2, buf2Stride[0], pu_height, pu_width);
                    *pBestMV = ((uint16_t)yMvQuarter[0] << 16) | ((uint16_t)xMvQuarter[0]);
                    *pBestSsd = (uint32_t)dist;
                }
            }
            else {
                if (dist < *pBestSad) {
                    *pBestSad = (uint32_t)dist;
                    *pBestMV = ((uint16_t)yMvQuarter[0] << 16) | ((uint16_t)xMvQuarter[0]);
                }
            }
#else
            if (dist < *pBestSad) {
                *pBestSad = (uint32_t)dist;
                *pBestMV = ((uint16_t)yMvQuarter[0] << 16) | ((uint16_t)xMvQuarter[0]);
            }
#endif
        }

        // R positions
        if (validR) {

            searchRegionIndex1 = (int32_t)xSearchIndex + (int32_t)buf1Stride[1] * (int32_t)ySearchIndex;
            searchRegionIndex2 = (int32_t)xSearchIndex + (int32_t)buf2Stride[1] * (int32_t)ySearchIndex;
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
            dist = (context_ptr->fractionalSearchMethod == SSD_SEARCH) ?
                combined_averaging_ssd_func_ptr_array[asm_type](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64, buf1[1] + searchRegionIndex1, buf1Stride[1], buf2[1] + searchRegionIndex2, buf2Stride[1], pu_height, pu_width) :
                (context_ptr->fractionalSearchMethod == SUB_SAD_SEARCH) ?
                (NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64 << 1, buf1[1] + searchRegionIndex1, buf1Stride[1] << 1, buf2[1] + searchRegionIndex2, buf2Stride[1] << 1, pu_height >> 1, pu_width)) << 1 :
                NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64, buf1[1] + searchRegionIndex1, buf1Stride[1], buf2[1] + searchRegionIndex2, buf2Stride[1], pu_height, pu_width);
#elif M0_SAD_HALF_QUARTER_PEL_BIPRED_SEARCH
            dist = (context_ptr->useSubSadFracBipredSearch) ?
                (NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64 << 1, buf1[1] + searchRegionIndex1, buf1Stride[1] << 1, buf2[1] + searchRegionIndex2, buf2Stride[1] << 1, pu_height >> 1, pu_width)) << 1 :
                NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64, buf1[1] + searchRegionIndex1, buf1Stride[1], buf2[1] + searchRegionIndex2, buf2Stride[1], pu_height, pu_width);
#else
            dist = NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64 << 1, buf1[1] + searchRegionIndex1, buf1Stride[1] << 1, buf2[1] + searchRegionIndex2, buf2Stride[1] << 1, pu_height >> 1, pu_width);
            dist = dist << 1;
#endif
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
            if (context_ptr->fractionalSearchMethod == SSD_SEARCH) {
                if (dist < *pBestSsd) {
                    *pBestSad = (uint32_t)NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64, buf1[1] + searchRegionIndex1, buf1Stride[1], buf2[1] + searchRegionIndex2, buf2Stride[1], pu_height, pu_width);
                    *pBestMV = ((uint16_t)yMvQuarter[1] << 16) | ((uint16_t)xMvQuarter[1]);
                    *pBestSsd = (uint32_t)dist;
                }
            }
            else {
                if (dist < *pBestSad) {
                    *pBestSad = (uint32_t)dist;
                    *pBestMV = ((uint16_t)yMvQuarter[1] << 16) | ((uint16_t)xMvQuarter[1]);
                }
            }
#else
            if (dist < *pBestSad) {
                *pBestSad = (uint32_t)dist;
                *pBestMV = ((uint16_t)yMvQuarter[1] << 16) | ((uint16_t)xMvQuarter[1]);
            }
#endif
        }

        // T position
        if (validT) {

            searchRegionIndex1 = (int32_t)xSearchIndex + (int32_t)buf1Stride[2] * (int32_t)ySearchIndex;
            searchRegionIndex2 = (int32_t)xSearchIndex + (int32_t)buf2Stride[2] * (int32_t)ySearchIndex;
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
            dist = (context_ptr->fractionalSearchMethod == SSD_SEARCH) ?
                combined_averaging_ssd_func_ptr_array[asm_type](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64, buf1[2] + searchRegionIndex1, buf1Stride[2], buf2[2] + searchRegionIndex2, buf2Stride[2], pu_height, pu_width) :
                (context_ptr->fractionalSearchMethod == SUB_SAD_SEARCH) ?
                (NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64 << 1, buf1[2] + searchRegionIndex1, buf1Stride[2] << 1, buf2[2] + searchRegionIndex2, buf2Stride[2] << 1, pu_height >> 1, pu_width)) << 1 :
                NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64, buf1[2] + searchRegionIndex1, buf1Stride[2], buf2[2] + searchRegionIndex2, buf2Stride[2], pu_height, pu_width);
#elif M0_SAD_HALF_QUARTER_PEL_BIPRED_SEARCH
            dist = (context_ptr->useSubSadFracBipredSearch) ?
                (NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64 << 1, buf1[2] + searchRegionIndex1, buf1Stride[2] << 1, buf2[2] + searchRegionIndex2, buf2Stride[2] << 1, pu_height >> 1, pu_width)) << 1 :
                NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64, buf1[2] + searchRegionIndex1, buf1Stride[2], buf2[2] + searchRegionIndex2, buf2Stride[2], pu_height, pu_width);
#else
            dist = NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64 << 1, buf1[2] + searchRegionIndex1, buf1Stride[2] << 1, buf2[2] + searchRegionIndex2, buf2Stride[2] << 1, pu_height >> 1, pu_width);
            dist = dist << 1;
#endif
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
            if (context_ptr->fractionalSearchMethod == SSD_SEARCH) {
                if (dist < *pBestSsd) {
                    *pBestSad = (uint32_t)NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64, buf1[2] + searchRegionIndex1, buf1Stride[2], buf2[2] + searchRegionIndex2, buf2Stride[2], pu_height, pu_width);
                    *pBestMV = ((uint16_t)yMvQuarter[2] << 16) | ((uint16_t)xMvQuarter[2]);
                    *pBestSsd = (uint32_t)dist;
                }
            }
            else {
                if (dist < *pBestSad) {
                    *pBestSad = (uint32_t)dist;
                    *pBestMV = ((uint16_t)yMvQuarter[2] << 16) | ((uint16_t)xMvQuarter[2]);
                }
            }
#else
            if (dist < *pBestSad) {
                *pBestSad = (uint32_t)dist;
                *pBestMV = ((uint16_t)yMvQuarter[2] << 16) | ((uint16_t)xMvQuarter[2]);
            }
#endif
        }

        // B position
        if (validB) {

            searchRegionIndex1 = (int32_t)xSearchIndex + (int32_t)buf1Stride[3] * (int32_t)ySearchIndex;
            searchRegionIndex2 = (int32_t)xSearchIndex + (int32_t)buf2Stride[3] * (int32_t)ySearchIndex;
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
            dist = (context_ptr->fractionalSearchMethod == SSD_SEARCH) ?
                combined_averaging_ssd_func_ptr_array[asm_type](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64, buf1[3] + searchRegionIndex1, buf1Stride[3], buf2[3] + searchRegionIndex2, buf2Stride[3], pu_height, pu_width) :
                (context_ptr->fractionalSearchMethod == SUB_SAD_SEARCH) ?
                (NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64 << 1, buf1[3] + searchRegionIndex1, buf1Stride[3] << 1, buf2[3] + searchRegionIndex2, buf2Stride[3] << 1, pu_height >> 1, pu_width)) << 1 :
                NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64, buf1[3] + searchRegionIndex1, buf1Stride[3], buf2[3] + searchRegionIndex2, buf2Stride[3], pu_height, pu_width);
#elif M0_SAD_HALF_QUARTER_PEL_BIPRED_SEARCH
            dist = (context_ptr->useSubSadFracBipredSearch) ?
                (NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64 << 1, buf1[3] + searchRegionIndex1, buf1Stride[3] << 1, buf2[3] + searchRegionIndex2, buf2Stride[3] << 1, pu_height >> 1, pu_width)) << 1 :
                NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64, buf1[3] + searchRegionIndex1, buf1Stride[3], buf2[3] + searchRegionIndex2, buf2Stride[3], pu_height, pu_width);
#else
            dist = NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64 << 1, buf1[3] + searchRegionIndex1, buf1Stride[3] << 1, buf2[3] + searchRegionIndex2, buf2Stride[3] << 1, pu_height >> 1, pu_width);
            dist = dist << 1;
#endif
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
            if (context_ptr->fractionalSearchMethod == SSD_SEARCH) {
                if (dist < *pBestSsd) {
                    *pBestSad = (uint32_t)NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64, buf1[3] + searchRegionIndex1, buf1Stride[3], buf2[3] + searchRegionIndex2, buf2Stride[3], pu_height, pu_width);
                    *pBestMV = ((uint16_t)yMvQuarter[3] << 16) | ((uint16_t)xMvQuarter[3]);
                    *pBestSsd = (uint32_t)dist;
                }
            }
            else {
                if (dist < *pBestSad) {
                    *pBestSad = (uint32_t)dist;
                    *pBestMV = ((uint16_t)yMvQuarter[3] << 16) | ((uint16_t)xMvQuarter[3]);
                }
            }
#else
            if (dist < *pBestSad) {
                *pBestSad = (uint32_t)dist;
                *pBestMV = ((uint16_t)yMvQuarter[3] << 16) | ((uint16_t)xMvQuarter[3]);
            }
#endif
        }

        //TL position
        if (validTL) {

            searchRegionIndex1 = (int32_t)xSearchIndex + (int32_t)buf1Stride[4] * (int32_t)ySearchIndex;
            searchRegionIndex2 = (int32_t)xSearchIndex + (int32_t)buf2Stride[4] * (int32_t)ySearchIndex;
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
            dist = (context_ptr->fractionalSearchMethod == SSD_SEARCH) ?
                combined_averaging_ssd_func_ptr_array[asm_type](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64, buf1[4] + searchRegionIndex1, buf1Stride[4], buf2[4] + searchRegionIndex2, buf2Stride[4], pu_height, pu_width) :
                (context_ptr->fractionalSearchMethod == SUB_SAD_SEARCH) ?
                (NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64 << 1, buf1[4] + searchRegionIndex1, buf1Stride[4] << 1, buf2[4] + searchRegionIndex2, buf2Stride[4] << 1, pu_height >> 1, pu_width)) << 1 :
                NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64, buf1[4] + searchRegionIndex1, buf1Stride[4], buf2[4] + searchRegionIndex2, buf2Stride[4], pu_height, pu_width);
#elif M0_SAD_HALF_QUARTER_PEL_BIPRED_SEARCH
            dist = (context_ptr->useSubSadFracBipredSearch) ?
                (NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64 << 1, buf1[4] + searchRegionIndex1, buf1Stride[4] << 1, buf2[4] + searchRegionIndex2, buf2Stride[4] << 1, pu_height >> 1, pu_width)) << 1 :
                NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64, buf1[4] + searchRegionIndex1, buf1Stride[4], buf2[4] + searchRegionIndex2, buf2Stride[4], pu_height, pu_width);
#else
            dist = NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64 << 1, buf1[4] + searchRegionIndex1, buf1Stride[4] << 1, buf2[4] + searchRegionIndex2, buf2Stride[4] << 1, pu_height >> 1, pu_width);
            dist = dist << 1;
#endif
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
            if (context_ptr->fractionalSearchMethod == SSD_SEARCH) {
                if (dist < *pBestSsd) {
                    *pBestSad = (uint32_t)NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64, buf1[4] + searchRegionIndex1, buf1Stride[4], buf2[4] + searchRegionIndex2, buf2Stride[4], pu_height, pu_width);
                    *pBestMV = ((uint16_t)yMvQuarter[4] << 16) | ((uint16_t)xMvQuarter[4]);
                    *pBestSsd = (uint32_t)dist;
                }
            }
            else {
                if (dist < *pBestSad) {
                    *pBestSad = (uint32_t)dist;
                    *pBestMV = ((uint16_t)yMvQuarter[4] << 16) | ((uint16_t)xMvQuarter[4]);
                }

            }
#else
            if (dist < *pBestSad) {
                *pBestSad = (uint32_t)dist;
                *pBestMV = ((uint16_t)yMvQuarter[4] << 16) | ((uint16_t)xMvQuarter[4]);
            }
#endif
        }

        //TR position
        if (validTR) {

            searchRegionIndex1 = (int32_t)xSearchIndex + (int32_t)buf1Stride[5] * (int32_t)ySearchIndex;
            searchRegionIndex2 = (int32_t)xSearchIndex + (int32_t)buf2Stride[5] * (int32_t)ySearchIndex;
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
            dist = (context_ptr->fractionalSearchMethod == SSD_SEARCH) ?
                combined_averaging_ssd_func_ptr_array[asm_type](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64, buf1[5] + searchRegionIndex1, buf1Stride[5], buf2[5] + searchRegionIndex2, buf2Stride[5], pu_height, pu_width) : \
                (context_ptr->fractionalSearchMethod == SUB_SAD_SEARCH) ?
                (NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64 << 1, buf1[5] + searchRegionIndex1, buf1Stride[5] << 1, buf2[5] + searchRegionIndex2, buf2Stride[5] << 1, pu_height >> 1, pu_width)) << 1 :
                NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64, buf1[5] + searchRegionIndex1, buf1Stride[5], buf2[5] + searchRegionIndex2, buf2Stride[5], pu_height, pu_width);
#elif M0_SAD_HALF_QUARTER_PEL_BIPRED_SEARCH
            dist = (context_ptr->useSubSadFracBipredSearch) ?
                (NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64 << 1, buf1[5] + searchRegionIndex1, buf1Stride[5] << 1, buf2[5] + searchRegionIndex2, buf2Stride[5] << 1, pu_height >> 1, pu_width)) << 1 :
                NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64, buf1[5] + searchRegionIndex1, buf1Stride[5], buf2[5] + searchRegionIndex2, buf2Stride[5], pu_height, pu_width);
#else
            dist = NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64 << 1, buf1[5] + searchRegionIndex1, buf1Stride[5] << 1, buf2[5] + searchRegionIndex2, buf2Stride[5] << 1, pu_height >> 1, pu_width);
            dist = dist << 1;
#endif
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
            if (context_ptr->fractionalSearchMethod == SSD_SEARCH) {
                if (dist < *pBestSsd) {
                    *pBestSad = (uint32_t)NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64, buf1[5] + searchRegionIndex1, buf1Stride[5], buf2[5] + searchRegionIndex2, buf2Stride[5], pu_height, pu_width);
                    *pBestMV = ((uint16_t)yMvQuarter[5] << 16) | ((uint16_t)xMvQuarter[5]);
                    *pBestSsd = (uint32_t)dist;
                }
            }
            else {
                if (dist < *pBestSad) {
                    *pBestSad = (uint32_t)dist;
                    *pBestMV = ((uint16_t)yMvQuarter[5] << 16) | ((uint16_t)xMvQuarter[5]);
                }
            }
#else
            if (dist < *pBestSad) {
                *pBestSad = (uint32_t)dist;
                *pBestMV = ((uint16_t)yMvQuarter[5] << 16) | ((uint16_t)xMvQuarter[5]);
            }
#endif
        }

        //BR position
        if (validBR) {

            searchRegionIndex1 = (int32_t)xSearchIndex + (int32_t)buf1Stride[6] * (int32_t)ySearchIndex;
            searchRegionIndex2 = (int32_t)xSearchIndex + (int32_t)buf2Stride[6] * (int32_t)ySearchIndex;
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
            dist = (context_ptr->fractionalSearchMethod == SSD_SEARCH) ?
                combined_averaging_ssd_func_ptr_array[asm_type](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64, buf1[6] + searchRegionIndex1, buf1Stride[6], buf2[6] + searchRegionIndex2, buf2Stride[6], pu_height, pu_width) :
                (context_ptr->fractionalSearchMethod == SUB_SAD_SEARCH) ?
                (NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64 << 1, buf1[6] + searchRegionIndex1, buf1Stride[6] << 1, buf2[6] + searchRegionIndex2, buf2Stride[6] << 1, pu_height >> 1, pu_width)) << 1 :
                NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64, buf1[6] + searchRegionIndex1, buf1Stride[6], buf2[6] + searchRegionIndex2, buf2Stride[6], pu_height, pu_width);
#elif M0_SAD_HALF_QUARTER_PEL_BIPRED_SEARCH
            dist = (context_ptr->useSubSadFracBipredSearch) ?
                (NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64 << 1, buf1[6] + searchRegionIndex1, buf1Stride[6] << 1, buf2[6] + searchRegionIndex2, buf2Stride[6] << 1, pu_height >> 1, pu_width)) << 1 :
                NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64, buf1[6] + searchRegionIndex1, buf1Stride[6], buf2[6] + searchRegionIndex2, buf2Stride[6], pu_height, pu_width);
#else
            dist = NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64 << 1, buf1[6] + searchRegionIndex1, buf1Stride[6] << 1, buf2[6] + searchRegionIndex2, buf2Stride[6] << 1, pu_height >> 1, pu_width);
            dist = dist << 1;
#endif
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
            if (context_ptr->fractionalSearchMethod == SSD_SEARCH) {
                if (dist < *pBestSsd) {
                    *pBestSad = (uint32_t)NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64, buf1[6] + searchRegionIndex1, buf1Stride[6], buf2[6] + searchRegionIndex2, buf2Stride[6], pu_height, pu_width);
                    *pBestMV = ((uint16_t)yMvQuarter[6] << 16) | ((uint16_t)xMvQuarter[6]);
                    *pBestSsd = (uint32_t)dist;
                }
            }
            else {
                if (dist < *pBestSad) {
                    *pBestSad = (uint32_t)dist;
                    *pBestMV = ((uint16_t)yMvQuarter[6] << 16) | ((uint16_t)xMvQuarter[6]);
                }
            }
#else
            if (dist < *pBestSad) {
                *pBestSad = (uint32_t)dist;
                *pBestMV = ((uint16_t)yMvQuarter[6] << 16) | ((uint16_t)xMvQuarter[6]);
            }
#endif
        }

        //BL position
        if (validBL) {

            searchRegionIndex1 = (int32_t)xSearchIndex + (int32_t)buf1Stride[7] * (int32_t)ySearchIndex;
            searchRegionIndex2 = (int32_t)xSearchIndex + (int32_t)buf2Stride[7] * (int32_t)ySearchIndex;
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
            dist = (context_ptr->fractionalSearchMethod == SSD_SEARCH) ?
                combined_averaging_ssd_func_ptr_array[asm_type](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64, buf1[7] + searchRegionIndex1, buf1Stride[7], buf2[7] + searchRegionIndex2, buf2Stride[7], pu_height, pu_width) :
                (context_ptr->fractionalSearchMethod == SUB_SAD_SEARCH) ?
                (NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64 << 1, buf1[7] + searchRegionIndex1, buf1Stride[7] << 1, buf2[7] + searchRegionIndex2, buf2Stride[7] << 1, pu_height >> 1, pu_width)) << 1 :
                NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64, buf1[7] + searchRegionIndex1, buf1Stride[7], buf2[7] + searchRegionIndex2, buf2Stride[7], pu_height, pu_width);
#elif M0_SAD_HALF_QUARTER_PEL_BIPRED_SEARCH
            dist = (context_ptr->useSubSadFracBipredSearch) ?
                (NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64 << 1, buf1[7] + searchRegionIndex1, buf1Stride[7] << 1, buf2[7] + searchRegionIndex2, buf2Stride[7] << 1, pu_height >> 1, pu_width)) << 1 :
                NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64, buf1[7] + searchRegionIndex1, buf1Stride[7], buf2[7] + searchRegionIndex2, buf2Stride[7], pu_height, pu_width);
#else
            dist = NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64 << 1, buf1[7] + searchRegionIndex1, buf1Stride[7] << 1, buf2[7] + searchRegionIndex2, buf2Stride[7] << 1, pu_height >> 1, pu_width);
            dist = dist << 1;
#endif
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
            if (context_ptr->fractionalSearchMethod == SSD_SEARCH) {
                if (dist < *pBestSsd) {
                    *pBestSad = (uint32_t)NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](&(context_ptr->sb_buffer[puLcuBufferIndex]), BLOCK_SIZE_64, buf1[7] + searchRegionIndex1, buf1Stride[7], buf2[7] + searchRegionIndex2, buf2Stride[7], pu_height, pu_width);
                    *pBestMV = ((uint16_t)yMvQuarter[7] << 16) | ((uint16_t)xMvQuarter[7]);
                    *pBestSsd = (uint32_t)dist;
                }
            }
            else {
                if (dist < *pBestSad) {
                    *pBestSad = (uint32_t)dist;
                    *pBestMV = ((uint16_t)yMvQuarter[7] << 16) | ((uint16_t)xMvQuarter[7]);
                }
            }
#else
            if (dist < *pBestSad) {
                *pBestSad = (uint32_t)dist;
                *pBestMV = ((uint16_t)yMvQuarter[7] << 16) | ((uint16_t)xMvQuarter[7]);
            }
#endif
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
    uint8_t  * pos_Full,   //[IN] points to A
    uint32_t   FullStride, //[IN]
    uint8_t  * pos_b,     //[IN] points to b
    uint8_t  * pos_h,     //[IN] points to h
    uint8_t  * pos_j,     //[IN] points to j
    uint32_t   Stride,    //[IN]
    int16_t   x_mv,        //[IN]
    int16_t   y_mv,        //[IN]
    uint8_t  **buf1,       //[OUT]
    uint32_t * buf1Stride, //[OUT]
    uint8_t  **buf2,       //[OUT]
    uint32_t*  buf2Stride  //[OUT]
)
{

    uint32_t  quarterPelRefinementMethod = (y_mv & 2) + ((x_mv & 2) >> 1);

    //for each one of the 8 postions, we need to determine the 2 half pel buffers to  do averaging

    //     A    a    b    c
    //     d    e    f    g
    //     h    i    j    k
    //     n    p    q    r

    switch (quarterPelRefinementMethod) {

    case EB_QUARTER_IN_FULL:

        /*c=b+A*/ buf1[0] = pos_b;                     buf1Stride[0] = Stride;        buf2[0] = pos_Full;             buf2Stride[0] = FullStride;
        /*a=A+b*/ buf1[1] = pos_Full;                  buf1Stride[1] = FullStride;    buf2[1] = pos_b + 1;             buf2Stride[1] = Stride;
        /*n=h+A*/ buf1[2] = pos_h;                      buf1Stride[2] = Stride;        buf2[2] = pos_Full;              buf2Stride[2] = FullStride;
        /*d=A+h*/ buf1[3] = pos_Full;                   buf1Stride[3] = FullStride;    buf2[3] = pos_h + Stride;        buf2Stride[3] = Stride;
        /*r=b+h*/ buf1[4] = pos_b;                      buf1Stride[4] = Stride;        buf2[4] = pos_h;                 buf2Stride[4] = Stride;
        /*p=h+b*/ buf1[5] = pos_h;                      buf1Stride[5] = Stride;        buf2[5] = pos_b + 1;             buf2Stride[5] = Stride;
        /*e=h+b*/ buf1[6] = pos_h + Stride;             buf1Stride[6] = Stride;        buf2[6] = pos_b + 1;             buf2Stride[6] = Stride;
        /*g=b+h*/ buf1[7] = pos_b;                      buf1Stride[7] = Stride;        buf2[7] = pos_h + Stride;        buf2Stride[7] = Stride;

        break;

    case EB_QUARTER_IN_HALF_HORIZONTAL:

        /*a=A+b*/ buf1[0] = pos_Full - 1;               buf1Stride[0] = FullStride;    buf2[0] = pos_b;                buf2Stride[0] = Stride;
        /*c=b+A*/ buf1[1] = pos_b;                     buf1Stride[1] = Stride;        buf2[1] = pos_Full;             buf2Stride[1] = FullStride;
        /*q=j+b*/ buf1[2] = pos_j;                     buf1Stride[2] = Stride;        buf2[2] = pos_b;                buf2Stride[2] = Stride;
        /*f=b+j*/ buf1[3] = pos_b;                     buf1Stride[3] = Stride;        buf2[3] = pos_j + Stride;        buf2Stride[3] = Stride;
        /*p=h+b*/ buf1[4] = pos_h - 1;                  buf1Stride[4] = Stride;        buf2[4] = pos_b;                buf2Stride[4] = Stride;
        /*r=b+h*/ buf1[5] = pos_b;                     buf1Stride[5] = Stride;        buf2[5] = pos_h;                buf2Stride[5] = Stride;
        /*g=b+h*/ buf1[6] = pos_b;                     buf1Stride[6] = Stride;        buf2[6] = pos_h + Stride;        buf2Stride[6] = Stride;
        /*e=h+b*/ buf1[7] = pos_h - 1 + Stride;         buf1Stride[7] = Stride;        buf2[7] = pos_b;                buf2Stride[7] = Stride;

        break;

    case EB_QUARTER_IN_HALF_VERTICAL:

        /*k=j+h*/buf1[0] = pos_j;                      buf1Stride[0] = Stride;        buf2[0] = pos_h;                 buf2Stride[0] = Stride;
        /*i=h+j*/buf1[1] = pos_h;                      buf1Stride[1] = Stride;        buf2[1] = pos_j + 1;              buf2Stride[1] = Stride;
        /*d=A+h*/buf1[2] = pos_Full - FullStride;      buf1Stride[2] = FullStride;    buf2[2] = pos_h;                  buf2Stride[2] = Stride;
        /*n=h+A*/buf1[3] = pos_h;                       buf1Stride[3] = Stride;        buf2[3] = pos_Full;               buf2Stride[3] = FullStride;
        /*g=b+h*/buf1[4] = pos_b - Stride;              buf1Stride[4] = Stride;        buf2[4] = pos_h;                  buf2Stride[4] = Stride;
        /*e=h+b*/buf1[5] = pos_h;                      buf1Stride[5] = Stride;        buf2[5] = pos_b + 1 - Stride;     buf2Stride[5] = Stride;
        /*p=h+b*/buf1[6] = pos_h;                      buf1Stride[6] = Stride;        buf2[6] = pos_b + 1;              buf2Stride[6] = Stride;
        /*r=b+h*/buf1[7] = pos_b;                      buf1Stride[7] = Stride;        buf2[7] = pos_h;                 buf2Stride[7] = Stride;

        break;

    case EB_QUARTER_IN_HALF_DIAGONAL:

        /*i=h+j*/buf1[0] = pos_h - 1;                   buf1Stride[0] = Stride;        buf2[0] = pos_j;                  buf2Stride[0] = Stride;
        /*k=j+h*/buf1[1] = pos_j;                       buf1Stride[1] = Stride;        buf2[1] = pos_h;                  buf2Stride[1] = Stride;
        /*f=b+j*/buf1[2] = pos_b - Stride;              buf1Stride[2] = Stride;        buf2[2] = pos_j;                  buf2Stride[2] = Stride;
        /*q=j+b*/buf1[3] = pos_j;                       buf1Stride[3] = Stride;        buf2[3] = pos_b;                  buf2Stride[3] = Stride;
        /*e=h+b*/buf1[4] = pos_h - 1;                   buf1Stride[4] = Stride;        buf2[4] = pos_b - Stride;         buf2Stride[4] = Stride;
        /*g=b+h*/buf1[5] = pos_b - Stride;              buf1Stride[5] = Stride;        buf2[5] = pos_h;                  buf2Stride[5] = Stride;
        /*r=b+h*/buf1[6] = pos_b;                      buf1Stride[6] = Stride;        buf2[6] = pos_h;                  buf2Stride[6] = Stride;
        /*p=h+b*/buf1[7] = pos_h - 1;                   buf1Stride[7] = Stride;        buf2[7] = pos_b;                  buf2Stride[7] = Stride;

        break;

    default:
        break;

    }

    return;
}

/*******************************************
* QuarterPelSearch_LCU
*   performs Quarter Pel refinement for the 85 PUs
*******************************************/
static void QuarterPelSearch_LCU(
    MeContext_t                    *context_ptr,                     //[IN/OUT]  ME context Ptr, used to get/update ME results
    uint8_t                        *pos_Full,                       //[IN]
    uint32_t                        FullStride,                      //[IN]
    uint8_t                        *pos_b,                          //[IN]
    uint8_t                        *pos_h,                          //[IN]
    uint8_t                        *pos_j,                          //[IN]
    int16_t                        x_search_area_origin,               //[IN] search area origin in the horizontal direction, used to point to reference samples
    int16_t                        y_search_area_origin,               //[IN] search area origin in the vertical direction, used to point to reference samples
    EbAsm                        asm_type,
    EbBool                        disable8x8CuInMeFlag,
    EbBool                        enableQuarterPel,
    EbBool                   ext_block_flag)
{
    uint32_t  pu_index;

    uint32_t  puShiftXIndex;
    uint32_t  puShiftYIndex;

    uint32_t  puLcuBufferIndex;

    //for each one of the 8 positions, we need to determine the 2 buffers to  do averaging
    uint8_t  *buf1[8];
    uint8_t  *buf2[8];

    uint32_t  buf1Stride[8];
    uint32_t  buf2Stride[8];

    int16_t  x_mv, y_mv;
    uint32_t  nidx;

#if M0_64x64_32x32_HALF_QUARTER_PEL
    if (context_ptr->fractional_search64x64) {
        x_mv = _MVXT(*context_ptr->p_best_mv64x64);
        y_mv = _MVYT(*context_ptr->p_best_mv64x64);

        SetQuarterPelRefinementInputsOnTheFly(
            pos_Full,
            FullStride,
            pos_b,
            pos_h,
            pos_j,
            context_ptr->interpolated_stride,
            x_mv,
            y_mv,
            buf1, buf1Stride,
            buf2, buf2Stride);

        buf1[0] = buf1[0];              buf2[0] = buf2[0];
        buf1[1] = buf1[1];              buf2[1] = buf2[1];
        buf1[2] = buf1[2];              buf2[2] = buf2[2];
        buf1[3] = buf1[3];              buf2[3] = buf2[3];
        buf1[4] = buf1[4];              buf2[4] = buf2[4];
        buf1[5] = buf1[5];              buf2[5] = buf2[5];
        buf1[6] = buf1[6];              buf2[6] = buf2[6];
        buf1[7] = buf1[7];              buf2[7] = buf2[7];


        PU_QuarterPelRefinementOnTheFly(
            context_ptr,
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
            context_ptr->p_best_ssd64x64,
#endif
            0,
            buf1, buf1Stride,
            buf2, buf2Stride,
            32, 32,
            x_search_area_origin,
            y_search_area_origin,
            asm_type,
            context_ptr->p_best_sad64x64,
            context_ptr->p_best_mv64x64,
            context_ptr->psub_pel_direction64x64);
    }
#else
    //no PU 64x64, Half Pel Refinement
#endif
    if (enableQuarterPel)
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
                buf1, buf1Stride,
                buf2, buf2Stride);


            puShiftXIndex = (pu_index & 0x01) << 5;
            puShiftYIndex = (pu_index >> 1) << 5;

            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + puShiftXIndex + puShiftYIndex * buf1Stride[0];              buf2[0] = buf2[0] + puShiftXIndex + puShiftYIndex * buf2Stride[0];
            buf1[1] = buf1[1] + puShiftXIndex + puShiftYIndex * buf1Stride[1];              buf2[1] = buf2[1] + puShiftXIndex + puShiftYIndex * buf2Stride[1];
            buf1[2] = buf1[2] + puShiftXIndex + puShiftYIndex * buf1Stride[2];              buf2[2] = buf2[2] + puShiftXIndex + puShiftYIndex * buf2Stride[2];
            buf1[3] = buf1[3] + puShiftXIndex + puShiftYIndex * buf1Stride[3];              buf2[3] = buf2[3] + puShiftXIndex + puShiftYIndex * buf2Stride[3];
            buf1[4] = buf1[4] + puShiftXIndex + puShiftYIndex * buf1Stride[4];              buf2[4] = buf2[4] + puShiftXIndex + puShiftYIndex * buf2Stride[4];
            buf1[5] = buf1[5] + puShiftXIndex + puShiftYIndex * buf1Stride[5];              buf2[5] = buf2[5] + puShiftXIndex + puShiftYIndex * buf2Stride[5];
            buf1[6] = buf1[6] + puShiftXIndex + puShiftYIndex * buf1Stride[6];              buf2[6] = buf2[6] + puShiftXIndex + puShiftYIndex * buf2Stride[6];
            buf1[7] = buf1[7] + puShiftXIndex + puShiftYIndex * buf1Stride[7];              buf2[7] = buf2[7] + puShiftXIndex + puShiftYIndex * buf2Stride[7];


            PU_QuarterPelRefinementOnTheFly(
                context_ptr,
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
                &context_ptr->p_best_ssd32x32[pu_index],
#endif
                puLcuBufferIndex,
                buf1, buf1Stride,
                buf2, buf2Stride,
                32, 32,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad32x32[pu_index],
                &context_ptr->p_best_mv32x32[pu_index],
                context_ptr->psub_pel_direction32x32[pu_index]);

        }
    }

    if (enableQuarterPel)
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
                buf1, buf1Stride,
                buf2, buf2Stride);


            puShiftXIndex = (pu_index & 0x03) << 4;
            puShiftYIndex = (pu_index >> 2) << 4;

            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + puShiftXIndex + puShiftYIndex * buf1Stride[0];              buf2[0] = buf2[0] + puShiftXIndex + puShiftYIndex * buf2Stride[0];
            buf1[1] = buf1[1] + puShiftXIndex + puShiftYIndex * buf1Stride[1];              buf2[1] = buf2[1] + puShiftXIndex + puShiftYIndex * buf2Stride[1];
            buf1[2] = buf1[2] + puShiftXIndex + puShiftYIndex * buf1Stride[2];              buf2[2] = buf2[2] + puShiftXIndex + puShiftYIndex * buf2Stride[2];
            buf1[3] = buf1[3] + puShiftXIndex + puShiftYIndex * buf1Stride[3];              buf2[3] = buf2[3] + puShiftXIndex + puShiftYIndex * buf2Stride[3];
            buf1[4] = buf1[4] + puShiftXIndex + puShiftYIndex * buf1Stride[4];              buf2[4] = buf2[4] + puShiftXIndex + puShiftYIndex * buf2Stride[4];
            buf1[5] = buf1[5] + puShiftXIndex + puShiftYIndex * buf1Stride[5];              buf2[5] = buf2[5] + puShiftXIndex + puShiftYIndex * buf2Stride[5];
            buf1[6] = buf1[6] + puShiftXIndex + puShiftYIndex * buf1Stride[6];              buf2[6] = buf2[6] + puShiftXIndex + puShiftYIndex * buf2Stride[6];
            buf1[7] = buf1[7] + puShiftXIndex + puShiftYIndex * buf1Stride[7];              buf2[7] = buf2[7] + puShiftXIndex + puShiftYIndex * buf2Stride[7];

            PU_QuarterPelRefinementOnTheFly(
                context_ptr,
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
                &context_ptr->p_best_ssd16x16[nidx],
#endif
                puLcuBufferIndex,
                buf1, buf1Stride,
                buf2, buf2Stride,
                16, 16,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad16x16[nidx],
                &context_ptr->p_best_mv16x16[nidx],
                context_ptr->psub_pel_direction16x16[nidx]);
        }
    }

    if (enableQuarterPel)
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
                    buf1, buf1Stride,
                    buf2, buf2Stride);


                puShiftXIndex = (pu_index & 0x07) << 3;
                puShiftYIndex = (pu_index >> 3) << 3;

                puLcuBufferIndex = puShiftXIndex + puShiftYIndex * BLOCK_SIZE_64;

                buf1[0] = buf1[0] + puShiftXIndex + puShiftYIndex * buf1Stride[0];              buf2[0] = buf2[0] + puShiftXIndex + puShiftYIndex * buf2Stride[0];
                buf1[1] = buf1[1] + puShiftXIndex + puShiftYIndex * buf1Stride[1];              buf2[1] = buf2[1] + puShiftXIndex + puShiftYIndex * buf2Stride[1];
                buf1[2] = buf1[2] + puShiftXIndex + puShiftYIndex * buf1Stride[2];              buf2[2] = buf2[2] + puShiftXIndex + puShiftYIndex * buf2Stride[2];
                buf1[3] = buf1[3] + puShiftXIndex + puShiftYIndex * buf1Stride[3];              buf2[3] = buf2[3] + puShiftXIndex + puShiftYIndex * buf2Stride[3];
                buf1[4] = buf1[4] + puShiftXIndex + puShiftYIndex * buf1Stride[4];              buf2[4] = buf2[4] + puShiftXIndex + puShiftYIndex * buf2Stride[4];
                buf1[5] = buf1[5] + puShiftXIndex + puShiftYIndex * buf1Stride[5];              buf2[5] = buf2[5] + puShiftXIndex + puShiftYIndex * buf2Stride[5];
                buf1[6] = buf1[6] + puShiftXIndex + puShiftYIndex * buf1Stride[6];              buf2[6] = buf2[6] + puShiftXIndex + puShiftYIndex * buf2Stride[6];
                buf1[7] = buf1[7] + puShiftXIndex + puShiftYIndex * buf1Stride[7];              buf2[7] = buf2[7] + puShiftXIndex + puShiftYIndex * buf2Stride[7];

                PU_QuarterPelRefinementOnTheFly(
                    context_ptr,
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
                    &context_ptr->p_best_ssd8x8[nidx],
#endif
                    puLcuBufferIndex,
                    buf1, buf1Stride,
                    buf2, buf2Stride,
                    8, 8,
                    x_search_area_origin,
                    y_search_area_origin,
                    asm_type,
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
                buf1, buf1Stride,
                buf2, buf2Stride);

            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + puShiftXIndex + puShiftYIndex * buf1Stride[0];              buf2[0] = buf2[0] + puShiftXIndex + puShiftYIndex * buf2Stride[0];
            buf1[1] = buf1[1] + puShiftXIndex + puShiftYIndex * buf1Stride[1];              buf2[1] = buf2[1] + puShiftXIndex + puShiftYIndex * buf2Stride[1];
            buf1[2] = buf1[2] + puShiftXIndex + puShiftYIndex * buf1Stride[2];              buf2[2] = buf2[2] + puShiftXIndex + puShiftYIndex * buf2Stride[2];
            buf1[3] = buf1[3] + puShiftXIndex + puShiftYIndex * buf1Stride[3];              buf2[3] = buf2[3] + puShiftXIndex + puShiftYIndex * buf2Stride[3];
            buf1[4] = buf1[4] + puShiftXIndex + puShiftYIndex * buf1Stride[4];              buf2[4] = buf2[4] + puShiftXIndex + puShiftYIndex * buf2Stride[4];
            buf1[5] = buf1[5] + puShiftXIndex + puShiftYIndex * buf1Stride[5];              buf2[5] = buf2[5] + puShiftXIndex + puShiftYIndex * buf2Stride[5];
            buf1[6] = buf1[6] + puShiftXIndex + puShiftYIndex * buf1Stride[6];              buf2[6] = buf2[6] + puShiftXIndex + puShiftYIndex * buf2Stride[6];
            buf1[7] = buf1[7] + puShiftXIndex + puShiftYIndex * buf1Stride[7];              buf2[7] = buf2[7] + puShiftXIndex + puShiftYIndex * buf2Stride[7];

            PU_QuarterPelRefinementOnTheFly(
                context_ptr,
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
                &context_ptr->p_best_ssd64x32[pu_index],
#endif
                puLcuBufferIndex,
                buf1, buf1Stride,
                buf2, buf2Stride,
                64,
                32,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad64x32[pu_index],
                &context_ptr->p_best_mv64x32[pu_index],
                context_ptr->psub_pel_direction64x32[pu_index]);

        }

        // 32x16
        for (pu_index = 0; pu_index < 8; ++pu_index) {

            nidx = tab32x16[pu_index];  //TODO bitwise this

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
                buf1, buf1Stride,
                buf2, buf2Stride);

            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + puShiftXIndex + puShiftYIndex * buf1Stride[0];              buf2[0] = buf2[0] + puShiftXIndex + puShiftYIndex * buf2Stride[0];
            buf1[1] = buf1[1] + puShiftXIndex + puShiftYIndex * buf1Stride[1];              buf2[1] = buf2[1] + puShiftXIndex + puShiftYIndex * buf2Stride[1];
            buf1[2] = buf1[2] + puShiftXIndex + puShiftYIndex * buf1Stride[2];              buf2[2] = buf2[2] + puShiftXIndex + puShiftYIndex * buf2Stride[2];
            buf1[3] = buf1[3] + puShiftXIndex + puShiftYIndex * buf1Stride[3];              buf2[3] = buf2[3] + puShiftXIndex + puShiftYIndex * buf2Stride[3];
            buf1[4] = buf1[4] + puShiftXIndex + puShiftYIndex * buf1Stride[4];              buf2[4] = buf2[4] + puShiftXIndex + puShiftYIndex * buf2Stride[4];
            buf1[5] = buf1[5] + puShiftXIndex + puShiftYIndex * buf1Stride[5];              buf2[5] = buf2[5] + puShiftXIndex + puShiftYIndex * buf2Stride[5];
            buf1[6] = buf1[6] + puShiftXIndex + puShiftYIndex * buf1Stride[6];              buf2[6] = buf2[6] + puShiftXIndex + puShiftYIndex * buf2Stride[6];
            buf1[7] = buf1[7] + puShiftXIndex + puShiftYIndex * buf1Stride[7];              buf2[7] = buf2[7] + puShiftXIndex + puShiftYIndex * buf2Stride[7];

            PU_QuarterPelRefinementOnTheFly(
                context_ptr,
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
                &context_ptr->p_best_ssd32x16[nidx],
#endif
                puLcuBufferIndex,
                buf1, buf1Stride,
                buf2, buf2Stride,
                32,
                16,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
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
                buf1, buf1Stride,
                buf2, buf2Stride);

            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + puShiftXIndex + puShiftYIndex * buf1Stride[0];              buf2[0] = buf2[0] + puShiftXIndex + puShiftYIndex * buf2Stride[0];
            buf1[1] = buf1[1] + puShiftXIndex + puShiftYIndex * buf1Stride[1];              buf2[1] = buf2[1] + puShiftXIndex + puShiftYIndex * buf2Stride[1];
            buf1[2] = buf1[2] + puShiftXIndex + puShiftYIndex * buf1Stride[2];              buf2[2] = buf2[2] + puShiftXIndex + puShiftYIndex * buf2Stride[2];
            buf1[3] = buf1[3] + puShiftXIndex + puShiftYIndex * buf1Stride[3];              buf2[3] = buf2[3] + puShiftXIndex + puShiftYIndex * buf2Stride[3];
            buf1[4] = buf1[4] + puShiftXIndex + puShiftYIndex * buf1Stride[4];              buf2[4] = buf2[4] + puShiftXIndex + puShiftYIndex * buf2Stride[4];
            buf1[5] = buf1[5] + puShiftXIndex + puShiftYIndex * buf1Stride[5];              buf2[5] = buf2[5] + puShiftXIndex + puShiftYIndex * buf2Stride[5];
            buf1[6] = buf1[6] + puShiftXIndex + puShiftYIndex * buf1Stride[6];              buf2[6] = buf2[6] + puShiftXIndex + puShiftYIndex * buf2Stride[6];
            buf1[7] = buf1[7] + puShiftXIndex + puShiftYIndex * buf1Stride[7];              buf2[7] = buf2[7] + puShiftXIndex + puShiftYIndex * buf2Stride[7];

            PU_QuarterPelRefinementOnTheFly(
                context_ptr,
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
                &context_ptr->p_best_ssd16x8[nidx],
#endif
                puLcuBufferIndex,
                buf1, buf1Stride,
                buf2, buf2Stride,
                16,
                8,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
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
                buf1, buf1Stride,
                buf2, buf2Stride);

            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + puShiftXIndex + puShiftYIndex * buf1Stride[0];              buf2[0] = buf2[0] + puShiftXIndex + puShiftYIndex * buf2Stride[0];
            buf1[1] = buf1[1] + puShiftXIndex + puShiftYIndex * buf1Stride[1];              buf2[1] = buf2[1] + puShiftXIndex + puShiftYIndex * buf2Stride[1];
            buf1[2] = buf1[2] + puShiftXIndex + puShiftYIndex * buf1Stride[2];              buf2[2] = buf2[2] + puShiftXIndex + puShiftYIndex * buf2Stride[2];
            buf1[3] = buf1[3] + puShiftXIndex + puShiftYIndex * buf1Stride[3];              buf2[3] = buf2[3] + puShiftXIndex + puShiftYIndex * buf2Stride[3];
            buf1[4] = buf1[4] + puShiftXIndex + puShiftYIndex * buf1Stride[4];              buf2[4] = buf2[4] + puShiftXIndex + puShiftYIndex * buf2Stride[4];
            buf1[5] = buf1[5] + puShiftXIndex + puShiftYIndex * buf1Stride[5];              buf2[5] = buf2[5] + puShiftXIndex + puShiftYIndex * buf2Stride[5];
            buf1[6] = buf1[6] + puShiftXIndex + puShiftYIndex * buf1Stride[6];              buf2[6] = buf2[6] + puShiftXIndex + puShiftYIndex * buf2Stride[6];
            buf1[7] = buf1[7] + puShiftXIndex + puShiftYIndex * buf1Stride[7];              buf2[7] = buf2[7] + puShiftXIndex + puShiftYIndex * buf2Stride[7];

            PU_QuarterPelRefinementOnTheFly(
                context_ptr,
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
                &context_ptr->p_best_ssd32x64[pu_index],
#endif
                puLcuBufferIndex,
                buf1, buf1Stride,
                buf2, buf2Stride,
                32,
                64,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
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
                buf1, buf1Stride,
                buf2, buf2Stride);

            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + puShiftXIndex + puShiftYIndex * buf1Stride[0];              buf2[0] = buf2[0] + puShiftXIndex + puShiftYIndex * buf2Stride[0];
            buf1[1] = buf1[1] + puShiftXIndex + puShiftYIndex * buf1Stride[1];              buf2[1] = buf2[1] + puShiftXIndex + puShiftYIndex * buf2Stride[1];
            buf1[2] = buf1[2] + puShiftXIndex + puShiftYIndex * buf1Stride[2];              buf2[2] = buf2[2] + puShiftXIndex + puShiftYIndex * buf2Stride[2];
            buf1[3] = buf1[3] + puShiftXIndex + puShiftYIndex * buf1Stride[3];              buf2[3] = buf2[3] + puShiftXIndex + puShiftYIndex * buf2Stride[3];
            buf1[4] = buf1[4] + puShiftXIndex + puShiftYIndex * buf1Stride[4];              buf2[4] = buf2[4] + puShiftXIndex + puShiftYIndex * buf2Stride[4];
            buf1[5] = buf1[5] + puShiftXIndex + puShiftYIndex * buf1Stride[5];              buf2[5] = buf2[5] + puShiftXIndex + puShiftYIndex * buf2Stride[5];
            buf1[6] = buf1[6] + puShiftXIndex + puShiftYIndex * buf1Stride[6];              buf2[6] = buf2[6] + puShiftXIndex + puShiftYIndex * buf2Stride[6];
            buf1[7] = buf1[7] + puShiftXIndex + puShiftYIndex * buf1Stride[7];              buf2[7] = buf2[7] + puShiftXIndex + puShiftYIndex * buf2Stride[7];

            PU_QuarterPelRefinementOnTheFly(
                context_ptr,
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
                &context_ptr->p_best_ssd16x32[nidx],
#endif
                puLcuBufferIndex,
                buf1, buf1Stride,
                buf2, buf2Stride,
                16,
                32,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
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
                buf1, buf1Stride,
                buf2, buf2Stride);

            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + puShiftXIndex + puShiftYIndex * buf1Stride[0];              buf2[0] = buf2[0] + puShiftXIndex + puShiftYIndex * buf2Stride[0];
            buf1[1] = buf1[1] + puShiftXIndex + puShiftYIndex * buf1Stride[1];              buf2[1] = buf2[1] + puShiftXIndex + puShiftYIndex * buf2Stride[1];
            buf1[2] = buf1[2] + puShiftXIndex + puShiftYIndex * buf1Stride[2];              buf2[2] = buf2[2] + puShiftXIndex + puShiftYIndex * buf2Stride[2];
            buf1[3] = buf1[3] + puShiftXIndex + puShiftYIndex * buf1Stride[3];              buf2[3] = buf2[3] + puShiftXIndex + puShiftYIndex * buf2Stride[3];
            buf1[4] = buf1[4] + puShiftXIndex + puShiftYIndex * buf1Stride[4];              buf2[4] = buf2[4] + puShiftXIndex + puShiftYIndex * buf2Stride[4];
            buf1[5] = buf1[5] + puShiftXIndex + puShiftYIndex * buf1Stride[5];              buf2[5] = buf2[5] + puShiftXIndex + puShiftYIndex * buf2Stride[5];
            buf1[6] = buf1[6] + puShiftXIndex + puShiftYIndex * buf1Stride[6];              buf2[6] = buf2[6] + puShiftXIndex + puShiftYIndex * buf2Stride[6];
            buf1[7] = buf1[7] + puShiftXIndex + puShiftYIndex * buf1Stride[7];              buf2[7] = buf2[7] + puShiftXIndex + puShiftYIndex * buf2Stride[7];

            PU_QuarterPelRefinementOnTheFly(
                context_ptr,
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
                &context_ptr->p_best_ssd8x16[nidx],
#endif
                puLcuBufferIndex,
                buf1, buf1Stride,
                buf2, buf2Stride,
                8,
                16,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
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
                buf1, buf1Stride,
                buf2, buf2Stride);

            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + puShiftXIndex + puShiftYIndex * buf1Stride[0];              buf2[0] = buf2[0] + puShiftXIndex + puShiftYIndex * buf2Stride[0];
            buf1[1] = buf1[1] + puShiftXIndex + puShiftYIndex * buf1Stride[1];              buf2[1] = buf2[1] + puShiftXIndex + puShiftYIndex * buf2Stride[1];
            buf1[2] = buf1[2] + puShiftXIndex + puShiftYIndex * buf1Stride[2];              buf2[2] = buf2[2] + puShiftXIndex + puShiftYIndex * buf2Stride[2];
            buf1[3] = buf1[3] + puShiftXIndex + puShiftYIndex * buf1Stride[3];              buf2[3] = buf2[3] + puShiftXIndex + puShiftYIndex * buf2Stride[3];
            buf1[4] = buf1[4] + puShiftXIndex + puShiftYIndex * buf1Stride[4];              buf2[4] = buf2[4] + puShiftXIndex + puShiftYIndex * buf2Stride[4];
            buf1[5] = buf1[5] + puShiftXIndex + puShiftYIndex * buf1Stride[5];              buf2[5] = buf2[5] + puShiftXIndex + puShiftYIndex * buf2Stride[5];
            buf1[6] = buf1[6] + puShiftXIndex + puShiftYIndex * buf1Stride[6];              buf2[6] = buf2[6] + puShiftXIndex + puShiftYIndex * buf2Stride[6];
            buf1[7] = buf1[7] + puShiftXIndex + puShiftYIndex * buf1Stride[7];              buf2[7] = buf2[7] + puShiftXIndex + puShiftYIndex * buf2Stride[7];

            PU_QuarterPelRefinementOnTheFly(
                context_ptr,
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
                &context_ptr->p_best_ssd32x8[nidx],
#endif
                puLcuBufferIndex,
                buf1, buf1Stride,
                buf2, buf2Stride,
                32,
                8,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
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
                buf1, buf1Stride,
                buf2, buf2Stride);

            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + puShiftXIndex + puShiftYIndex * buf1Stride[0];              buf2[0] = buf2[0] + puShiftXIndex + puShiftYIndex * buf2Stride[0];
            buf1[1] = buf1[1] + puShiftXIndex + puShiftYIndex * buf1Stride[1];              buf2[1] = buf2[1] + puShiftXIndex + puShiftYIndex * buf2Stride[1];
            buf1[2] = buf1[2] + puShiftXIndex + puShiftYIndex * buf1Stride[2];              buf2[2] = buf2[2] + puShiftXIndex + puShiftYIndex * buf2Stride[2];
            buf1[3] = buf1[3] + puShiftXIndex + puShiftYIndex * buf1Stride[3];              buf2[3] = buf2[3] + puShiftXIndex + puShiftYIndex * buf2Stride[3];
            buf1[4] = buf1[4] + puShiftXIndex + puShiftYIndex * buf1Stride[4];              buf2[4] = buf2[4] + puShiftXIndex + puShiftYIndex * buf2Stride[4];
            buf1[5] = buf1[5] + puShiftXIndex + puShiftYIndex * buf1Stride[5];              buf2[5] = buf2[5] + puShiftXIndex + puShiftYIndex * buf2Stride[5];
            buf1[6] = buf1[6] + puShiftXIndex + puShiftYIndex * buf1Stride[6];              buf2[6] = buf2[6] + puShiftXIndex + puShiftYIndex * buf2Stride[6];
            buf1[7] = buf1[7] + puShiftXIndex + puShiftYIndex * buf1Stride[7];              buf2[7] = buf2[7] + puShiftXIndex + puShiftYIndex * buf2Stride[7];

            PU_QuarterPelRefinementOnTheFly(
                context_ptr,
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
                &context_ptr->p_best_ssd8x32[nidx],
#endif
                puLcuBufferIndex,
                buf1, buf1Stride,
                buf2, buf2Stride,
                8,
                32,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
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
                buf1, buf1Stride,
                buf2, buf2Stride);

            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + puShiftXIndex + puShiftYIndex * buf1Stride[0];              buf2[0] = buf2[0] + puShiftXIndex + puShiftYIndex * buf2Stride[0];
            buf1[1] = buf1[1] + puShiftXIndex + puShiftYIndex * buf1Stride[1];              buf2[1] = buf2[1] + puShiftXIndex + puShiftYIndex * buf2Stride[1];
            buf1[2] = buf1[2] + puShiftXIndex + puShiftYIndex * buf1Stride[2];              buf2[2] = buf2[2] + puShiftXIndex + puShiftYIndex * buf2Stride[2];
            buf1[3] = buf1[3] + puShiftXIndex + puShiftYIndex * buf1Stride[3];              buf2[3] = buf2[3] + puShiftXIndex + puShiftYIndex * buf2Stride[3];
            buf1[4] = buf1[4] + puShiftXIndex + puShiftYIndex * buf1Stride[4];              buf2[4] = buf2[4] + puShiftXIndex + puShiftYIndex * buf2Stride[4];
            buf1[5] = buf1[5] + puShiftXIndex + puShiftYIndex * buf1Stride[5];              buf2[5] = buf2[5] + puShiftXIndex + puShiftYIndex * buf2Stride[5];
            buf1[6] = buf1[6] + puShiftXIndex + puShiftYIndex * buf1Stride[6];              buf2[6] = buf2[6] + puShiftXIndex + puShiftYIndex * buf2Stride[6];
            buf1[7] = buf1[7] + puShiftXIndex + puShiftYIndex * buf1Stride[7];              buf2[7] = buf2[7] + puShiftXIndex + puShiftYIndex * buf2Stride[7];

            PU_QuarterPelRefinementOnTheFly(
                context_ptr,
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
                &context_ptr->p_best_ssd64x16[nidx],
#endif
                puLcuBufferIndex,
                buf1, buf1Stride,
                buf2, buf2Stride,
                64,
                16,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
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
                buf1, buf1Stride,
                buf2, buf2Stride);

            puLcuBufferIndex = puShiftXIndex + puShiftYIndex * BLOCK_SIZE_64;

            buf1[0] = buf1[0] + puShiftXIndex + puShiftYIndex * buf1Stride[0];              buf2[0] = buf2[0] + puShiftXIndex + puShiftYIndex * buf2Stride[0];
            buf1[1] = buf1[1] + puShiftXIndex + puShiftYIndex * buf1Stride[1];              buf2[1] = buf2[1] + puShiftXIndex + puShiftYIndex * buf2Stride[1];
            buf1[2] = buf1[2] + puShiftXIndex + puShiftYIndex * buf1Stride[2];              buf2[2] = buf2[2] + puShiftXIndex + puShiftYIndex * buf2Stride[2];
            buf1[3] = buf1[3] + puShiftXIndex + puShiftYIndex * buf1Stride[3];              buf2[3] = buf2[3] + puShiftXIndex + puShiftYIndex * buf2Stride[3];
            buf1[4] = buf1[4] + puShiftXIndex + puShiftYIndex * buf1Stride[4];              buf2[4] = buf2[4] + puShiftXIndex + puShiftYIndex * buf2Stride[4];
            buf1[5] = buf1[5] + puShiftXIndex + puShiftYIndex * buf1Stride[5];              buf2[5] = buf2[5] + puShiftXIndex + puShiftYIndex * buf2Stride[5];
            buf1[6] = buf1[6] + puShiftXIndex + puShiftYIndex * buf1Stride[6];              buf2[6] = buf2[6] + puShiftXIndex + puShiftYIndex * buf2Stride[6];
            buf1[7] = buf1[7] + puShiftXIndex + puShiftYIndex * buf1Stride[7];              buf2[7] = buf2[7] + puShiftXIndex + puShiftYIndex * buf2Stride[7];

            PU_QuarterPelRefinementOnTheFly(
                context_ptr,
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
                &context_ptr->p_best_ssd16x64[nidx],
#endif
                puLcuBufferIndex,
                buf1, buf1Stride,
                buf2, buf2Stride,
                16,
                64,
                x_search_area_origin,
                y_search_area_origin,
                asm_type,
                &context_ptr->p_best_sad16x64[nidx],
                &context_ptr->p_best_mv16x64[nidx],
                context_ptr->psub_pel_direction16x64[nidx]);
        }
    }

    return;
}
#endif
void HmeOneQuadrantLevel0(
    PictureParentControlSet_t   *picture_control_set_ptr,
    MeContext_t             *context_ptr,                        // input/output parameter, ME context Ptr, used to get/update ME results
    int16_t                   origin_x,                           // input parameter, SB position in the horizontal direction- sixteenth resolution
    int16_t                   origin_y,                           // input parameter, SB position in the vertical direction- sixteenth resolution
    uint32_t                   sb_width,                          // input parameter, SB pwidth - sixteenth resolution
    uint32_t                   sb_height,                         // input parameter, SB height - sixteenth resolution
    int16_t                   xHmeSearchCenter,                  // input parameter, HME search center in the horizontal direction
    int16_t                   yHmeSearchCenter,                  // input parameter, HME search center in the vertical direction
    EbPictureBufferDesc_t   *sixteenthRefPicPtr,                // input parameter, sixteenth reference Picture Ptr
    uint64_t                  *level0BestSad,                     // output parameter, Level0 SAD at (searchRegionNumberInWidth, searchRegionNumberInHeight)
    int16_t                  *xLevel0SearchCenter,               // output parameter, Level0 xMV at (searchRegionNumberInWidth, searchRegionNumberInHeight)
    int16_t                  *yLevel0SearchCenter,               // output parameter, Level0 yMV at (searchRegionNumberInWidth, searchRegionNumberInHeight)
    uint32_t                   searchAreaMultiplierX,
    uint32_t                   searchAreaMultiplierY,
    EbAsm                   asm_type)
{

    int16_t xTopLeftSearchRegion;
    int16_t yTopLeftSearchRegion;
    uint32_t searchRegionIndex;
    int16_t x_search_area_origin;
    int16_t y_search_area_origin;
    int16_t xSearchRegionDistance;
    int16_t ySearchRegionDistance;

    int16_t padWidth;
    int16_t padHeight;


#if ME_HME_OQ
    (void)picture_control_set_ptr;
    int16_t search_area_width = (int16_t)(((context_ptr->hme_level0_total_search_area_width  * searchAreaMultiplierX) / 100));
    int16_t search_area_height = (int16_t)(((context_ptr->hme_level0_total_search_area_height * searchAreaMultiplierY) / 100));
#else
    int16_t search_area_width = (int16_t)(((picture_control_set_ptr->hme_level0_total_search_area_width  * searchAreaMultiplierX) / 100));
    int16_t search_area_height = (int16_t)(((picture_control_set_ptr->hme_level0_total_search_area_height * searchAreaMultiplierY) / 100));

#endif
    if (context_ptr->hme_search_type == HME_SPARSE)
        search_area_width = ((search_area_width + 4) >> 3) << 3;  //round down/up the width to the nearest multiple of 8.

    xSearchRegionDistance = xHmeSearchCenter;
    ySearchRegionDistance = yHmeSearchCenter;
    padWidth = (int16_t)(sixteenthRefPicPtr->origin_x) - 1;
    padHeight = (int16_t)(sixteenthRefPicPtr->origin_y) - 1;


    x_search_area_origin = -(int16_t)(search_area_width >> 1) + xSearchRegionDistance;
    y_search_area_origin = -(int16_t)(search_area_height >> 1) + ySearchRegionDistance;

    // Correct the left edge of the Search Area if it is not on the reference Picture
    x_search_area_origin = ((origin_x + x_search_area_origin) < -padWidth) ?
        -padWidth - origin_x :
        x_search_area_origin;

    search_area_width = ((origin_x + x_search_area_origin) < -padWidth) ?
        search_area_width - (-padWidth - (origin_x + x_search_area_origin)) :
        search_area_width;

    // Correct the right edge of the Search Area if its not on the reference Picture
    x_search_area_origin = ((origin_x + x_search_area_origin) > (int16_t)sixteenthRefPicPtr->width - 1) ?
        x_search_area_origin - ((origin_x + x_search_area_origin) - ((int16_t)sixteenthRefPicPtr->width - 1)) :
        x_search_area_origin;

    search_area_width = ((origin_x + x_search_area_origin + search_area_width) > (int16_t)sixteenthRefPicPtr->width) ?
        MAX(1, search_area_width - ((origin_x + x_search_area_origin + search_area_width) - (int16_t)sixteenthRefPicPtr->width)) :
        search_area_width;

    // Correct the top edge of the Search Area if it is not on the reference Picture
    y_search_area_origin = ((origin_y + y_search_area_origin) < -padHeight) ?
        -padHeight - origin_y :
        y_search_area_origin;

    search_area_height = ((origin_y + y_search_area_origin) < -padHeight) ?
        search_area_height - (-padHeight - (origin_y + y_search_area_origin)) :
        search_area_height;

    // Correct the bottom edge of the Search Area if its not on the reference Picture
    y_search_area_origin = ((origin_y + y_search_area_origin) > (int16_t)sixteenthRefPicPtr->height - 1) ?
        y_search_area_origin - ((origin_y + y_search_area_origin) - ((int16_t)sixteenthRefPicPtr->height - 1)) :
        y_search_area_origin;

    search_area_height = (origin_y + y_search_area_origin + search_area_height > (int16_t)sixteenthRefPicPtr->height) ?
        MAX(1, search_area_height - ((origin_y + y_search_area_origin + search_area_height) - (int16_t)sixteenthRefPicPtr->height)) :
        search_area_height;

    xTopLeftSearchRegion = ((int16_t)sixteenthRefPicPtr->origin_x + origin_x) + x_search_area_origin;
    yTopLeftSearchRegion = ((int16_t)sixteenthRefPicPtr->origin_y + origin_y) + y_search_area_origin;
    searchRegionIndex = xTopLeftSearchRegion + yTopLeftSearchRegion * sixteenthRefPicPtr->strideY;

    if (context_ptr->hme_search_type == HME_SPARSE)
    {
        //ensure that search area is multiple of 8.
        search_area_width = ((search_area_width >> 3) << 3);

        NxMSadLoopKernelSparse_funcPtrArray[asm_type](
            &context_ptr->sixteenth_sb_buffer[0],
            context_ptr->sixteenth_sb_buffer_stride,
            &sixteenthRefPicPtr->bufferY[searchRegionIndex],
            sixteenthRefPicPtr->strideY * 2,
            sb_height >> 1, sb_width,
            /* results */
            level0BestSad,
            xLevel0SearchCenter,
            yLevel0SearchCenter,
            /* range */
            sixteenthRefPicPtr->strideY,
            search_area_width,
            search_area_height
            );


    }
    else {

        if ((search_area_width & 15) != 0)
        {
            search_area_width = (int16_t)(floor((double)((search_area_width >> 4) << 4)));
        }

        if (((search_area_width & 15) == 0) && (asm_type == ASM_AVX2))
        {
            SadLoopKernel_AVX2_HmeL0_INTRIN(
                &context_ptr->sixteenth_sb_buffer[0],
                context_ptr->sixteenth_sb_buffer_stride,
                &sixteenthRefPicPtr->bufferY[searchRegionIndex],
                sixteenthRefPicPtr->strideY * 2,
                sb_height >> 1, sb_width,
                // results
                level0BestSad,
                xLevel0SearchCenter,
                yLevel0SearchCenter,
                // range
                sixteenthRefPicPtr->strideY,
                search_area_width,
                search_area_height
            );
        }
        else if ((search_area_width & 15) == 0)
        {
            // Only width equals 16 (LCU equals 64) is updated
            // other width sizes work with the old code as the one in"SadLoopKernel_SSE4_1_INTRIN"
            SadLoopKernel_SSE4_1_HmeL0_INTRIN(
                &context_ptr->sixteenth_sb_buffer[0],
                context_ptr->sixteenth_sb_buffer_stride,
                &sixteenthRefPicPtr->bufferY[searchRegionIndex],
                sixteenthRefPicPtr->strideY * 2,
                sb_height >> 1, sb_width,
                /* results */
                level0BestSad,
                xLevel0SearchCenter,
                yLevel0SearchCenter,
                /* range */
                sixteenthRefPicPtr->strideY,
                search_area_width,
                search_area_height
            );
        }
        else
        {
            // Put the first search location into level0 results
            NxMSadLoopKernel_funcPtrArray[asm_type](
                &context_ptr->sixteenth_sb_buffer[0],
                context_ptr->sixteenth_sb_buffer_stride,
                &sixteenthRefPicPtr->bufferY[searchRegionIndex],
                sixteenthRefPicPtr->strideY * 2,
                sb_height >> 1, sb_width,
                /* results */
                level0BestSad,
                xLevel0SearchCenter,
                yLevel0SearchCenter,
                /* range */
                sixteenthRefPicPtr->strideY,
                search_area_width,
                search_area_height
                );
        }
    }

    *level0BestSad *= 2; // Multiply by 2 because considered only ever other line
    *xLevel0SearchCenter += x_search_area_origin;
    *xLevel0SearchCenter *= 4; // Multiply by 4 because operating on 1/4 resolution
    *yLevel0SearchCenter += y_search_area_origin;
    *yLevel0SearchCenter *= 4; // Multiply by 4 because operating on 1/4 resolution

    return;
}


void HmeLevel0(
    PictureParentControlSet_t   *picture_control_set_ptr,
    MeContext_t             *context_ptr,                        // input/output parameter, ME context Ptr, used to get/update ME results
    int16_t                   origin_x,                           // input parameter, SB position in the horizontal direction- sixteenth resolution
    int16_t                   origin_y,                           // input parameter, SB position in the vertical direction- sixteenth resolution
    uint32_t                   sb_width,                          // input parameter, SB pwidth - sixteenth resolution
    uint32_t                   sb_height,                         // input parameter, SB height - sixteenth resolution
    int16_t                   xHmeSearchCenter,                  // input parameter, HME search center in the horizontal direction
    int16_t                   yHmeSearchCenter,                  // input parameter, HME search center in the vertical direction
    EbPictureBufferDesc_t   *sixteenthRefPicPtr,                // input parameter, sixteenth reference Picture Ptr
    uint32_t                   searchRegionNumberInWidth,         // input parameter, search region number in the horizontal direction
    uint32_t                   searchRegionNumberInHeight,        // input parameter, search region number in the vertical direction
    uint64_t                  *level0BestSad,                     // output parameter, Level0 SAD at (searchRegionNumberInWidth, searchRegionNumberInHeight)
    int16_t                  *xLevel0SearchCenter,               // output parameter, Level0 xMV at (searchRegionNumberInWidth, searchRegionNumberInHeight)
    int16_t                  *yLevel0SearchCenter,               // output parameter, Level0 yMV at (searchRegionNumberInWidth, searchRegionNumberInHeight)
    uint32_t                   searchAreaMultiplierX,
    uint32_t                   searchAreaMultiplierY,
    EbAsm                   asm_type)
{

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

#if ME_HME_OQ
    (void)picture_control_set_ptr;
    int16_t search_area_width = (int16_t)(((context_ptr->hme_level0_search_area_in_width_array[searchRegionNumberInWidth] * searchAreaMultiplierX) / 100));
    int16_t search_area_height = (int16_t)(((context_ptr->hme_level0_search_area_in_height_array[searchRegionNumberInHeight] * searchAreaMultiplierY) / 100));
#else
    int16_t search_area_width = (int16_t)(((picture_control_set_ptr->hme_level0_search_area_in_width_array[searchRegionNumberInWidth] * searchAreaMultiplierX) / 100));
    int16_t search_area_height = (int16_t)(((picture_control_set_ptr->hme_level0_search_area_in_height_array[searchRegionNumberInHeight] * searchAreaMultiplierY) / 100));
#endif
    xSearchRegionDistance = xHmeSearchCenter;
    ySearchRegionDistance = yHmeSearchCenter;
    padWidth = (int16_t)(sixteenthRefPicPtr->origin_x) - 1;
    padHeight = (int16_t)(sixteenthRefPicPtr->origin_y) - 1;

    while (searchRegionNumberInWidth) {
        searchRegionNumberInWidth--;
#if ME_HME_OQ
        xSearchRegionDistance += (int16_t)(((context_ptr->hme_level0_search_area_in_width_array[searchRegionNumberInWidth] * searchAreaMultiplierX) / 100));
#else
        xSearchRegionDistance += (int16_t)(((picture_control_set_ptr->hme_level0_search_area_in_width_array[searchRegionNumberInWidth] * searchAreaMultiplierX) / 100));
#endif
    }

    while (searchRegionNumberInHeight) {
        searchRegionNumberInHeight--;
#if ME_HME_OQ
        ySearchRegionDistance += (int16_t)(((context_ptr->hme_level0_search_area_in_height_array[searchRegionNumberInHeight] * searchAreaMultiplierY) / 100));
#else
        ySearchRegionDistance += (int16_t)(((picture_control_set_ptr->hme_level0_search_area_in_height_array[searchRegionNumberInHeight] * searchAreaMultiplierY) / 100));
#endif
    }
#if ME_HME_OQ
    x_search_area_origin = -(int16_t)((((context_ptr->hme_level0_total_search_area_width * searchAreaMultiplierX) / 100)) >> 1) + xSearchRegionDistance;
    y_search_area_origin = -(int16_t)((((context_ptr->hme_level0_total_search_area_height * searchAreaMultiplierY) / 100)) >> 1) + ySearchRegionDistance;
#else
    x_search_area_origin = -(int16_t)((((picture_control_set_ptr->hme_level0_total_search_area_width * searchAreaMultiplierX) / 100)) >> 1) + xSearchRegionDistance;
    y_search_area_origin = -(int16_t)((((picture_control_set_ptr->hme_level0_total_search_area_height * searchAreaMultiplierY) / 100)) >> 1) + ySearchRegionDistance;

#endif
    // Correct the left edge of the Search Area if it is not on the reference Picture
    x_search_area_origin = ((origin_x + x_search_area_origin) < -padWidth) ?
        -padWidth - origin_x :
        x_search_area_origin;

    search_area_width = ((origin_x + x_search_area_origin) < -padWidth) ?
        search_area_width - (-padWidth - (origin_x + x_search_area_origin)) :
        search_area_width;

    // Correct the right edge of the Search Area if its not on the reference Picture
    x_search_area_origin = ((origin_x + x_search_area_origin) > (int16_t)sixteenthRefPicPtr->width - 1) ?
        x_search_area_origin - ((origin_x + x_search_area_origin) - ((int16_t)sixteenthRefPicPtr->width - 1)) :
        x_search_area_origin;

    search_area_width = ((origin_x + x_search_area_origin + search_area_width) > (int16_t)sixteenthRefPicPtr->width) ?
        MAX(1, search_area_width - ((origin_x + x_search_area_origin + search_area_width) - (int16_t)sixteenthRefPicPtr->width)) :
        search_area_width;

    // Correct the top edge of the Search Area if it is not on the reference Picture
    y_search_area_origin = ((origin_y + y_search_area_origin) < -padHeight) ?
        -padHeight - origin_y :
        y_search_area_origin;

    search_area_height = ((origin_y + y_search_area_origin) < -padHeight) ?
        search_area_height - (-padHeight - (origin_y + y_search_area_origin)) :
        search_area_height;

    // Correct the bottom edge of the Search Area if its not on the reference Picture
    y_search_area_origin = ((origin_y + y_search_area_origin) > (int16_t)sixteenthRefPicPtr->height - 1) ?
        y_search_area_origin - ((origin_y + y_search_area_origin) - ((int16_t)sixteenthRefPicPtr->height - 1)) :
        y_search_area_origin;

    search_area_height = (origin_y + y_search_area_origin + search_area_height > (int16_t)sixteenthRefPicPtr->height) ?
        MAX(1, search_area_height - ((origin_y + y_search_area_origin + search_area_height) - (int16_t)sixteenthRefPicPtr->height)) :
        search_area_height;

    xTopLeftSearchRegion = ((int16_t)sixteenthRefPicPtr->origin_x + origin_x) + x_search_area_origin;
    yTopLeftSearchRegion = ((int16_t)sixteenthRefPicPtr->origin_y + origin_y) + y_search_area_origin;
    searchRegionIndex = xTopLeftSearchRegion + yTopLeftSearchRegion * sixteenthRefPicPtr->strideY;

    if (((sb_width & 7) == 0) || (sb_width == 4))
    {
        if (((search_area_width & 15) == 0) && (asm_type == ASM_AVX2))
        {
            SadLoopKernel_AVX2_HmeL0_INTRIN(
                &context_ptr->sixteenth_sb_buffer[0],
                context_ptr->sixteenth_sb_buffer_stride,
                &sixteenthRefPicPtr->bufferY[searchRegionIndex],
                sixteenthRefPicPtr->strideY * 2,
                sb_height >> 1, sb_width,
                // results
                level0BestSad,
                xLevel0SearchCenter,
                yLevel0SearchCenter,
                // range
                sixteenthRefPicPtr->strideY,
                search_area_width,
                search_area_height
            );
        }
        else if ((search_area_width & 15) == 0)
        {
            // Only width equals 16 (LCU equals 64) is updated
            // other width sizes work with the old code as the one in"SadLoopKernel_SSE4_1_INTRIN"
            SadLoopKernel_SSE4_1_HmeL0_INTRIN(
                &context_ptr->sixteenth_sb_buffer[0],
                context_ptr->sixteenth_sb_buffer_stride,
                &sixteenthRefPicPtr->bufferY[searchRegionIndex],
                sixteenthRefPicPtr->strideY * 2,
                sb_height >> 1, sb_width,
                /* results */
                level0BestSad,
                xLevel0SearchCenter,
                yLevel0SearchCenter,
                /* range */
                sixteenthRefPicPtr->strideY,
                search_area_width,
                search_area_height
            );
        }
        else
        {
            // Put the first search location into level0 results
            NxMSadLoopKernel_funcPtrArray[asm_type](
                &context_ptr->sixteenth_sb_buffer[0],
                context_ptr->sixteenth_sb_buffer_stride,
                &sixteenthRefPicPtr->bufferY[searchRegionIndex],
                sixteenthRefPicPtr->strideY * 2,
                sb_height >> 1, sb_width,
                /* results */
                level0BestSad,
                xLevel0SearchCenter,
                yLevel0SearchCenter,
                /* range */
                sixteenthRefPicPtr->strideY,
                search_area_width,
                search_area_height
                );
        }
    }
    else
    {
        SadLoopKernel(
            &context_ptr->sixteenth_sb_buffer[0],
            context_ptr->sixteenth_sb_buffer_stride,
            &sixteenthRefPicPtr->bufferY[searchRegionIndex],
            sixteenthRefPicPtr->strideY * 2,
            sb_height >> 1, sb_width,
            /* results */
            level0BestSad,
            xLevel0SearchCenter,
            yLevel0SearchCenter,
            /* range */
            sixteenthRefPicPtr->strideY,
            search_area_width,
            search_area_height
        );
    }

    *level0BestSad *= 2; // Multiply by 2 because considered only ever other line
    *xLevel0SearchCenter += x_search_area_origin;
    *xLevel0SearchCenter *= 4; // Multiply by 4 because operating on 1/4 resolution
    *yLevel0SearchCenter += y_search_area_origin;
    *yLevel0SearchCenter *= 4; // Multiply by 4 because operating on 1/4 resolution

    return;
}

void HmeLevel1(
    MeContext_t             *context_ptr,                        // input/output parameter, ME context Ptr, used to get/update ME results
    int16_t                   origin_x,                           // input parameter, SB position in the horizontal direction - quarter resolution
    int16_t                   origin_y,                           // input parameter, SB position in the vertical direction - quarter resolution
    uint32_t                   sb_width,                          // input parameter, SB pwidth - quarter resolution
    uint32_t                   sb_height,                         // input parameter, SB height - quarter resolution
    EbPictureBufferDesc_t   *quarterRefPicPtr,                  // input parameter, quarter reference Picture Ptr
    int16_t                   hmeLevel1SearchAreaInWidth,         // input parameter, hme level 1 search area in width
    int16_t                   hmeLevel1SearchAreaInHeight,        // input parameter, hme level 1 search area in height
    int16_t                   xLevel0SearchCenter,               // input parameter, best Level0 xMV at (searchRegionNumberInWidth, searchRegionNumberInHeight)
    int16_t                   yLevel0SearchCenter,               // input parameter, best Level0 yMV at (searchRegionNumberInWidth, searchRegionNumberInHeight)
    uint64_t                  *level1BestSad,                     // output parameter, Level1 SAD at (searchRegionNumberInWidth, searchRegionNumberInHeight)
    int16_t                  *xLevel1SearchCenter,               // output parameter, Level1 xMV at (searchRegionNumberInWidth, searchRegionNumberInHeight)
    int16_t                  *yLevel1SearchCenter,               // output parameter, Level1 yMV at (searchRegionNumberInWidth, searchRegionNumberInHeight)
    EbAsm                   asm_type)
{


    int16_t xTopLeftSearchRegion;
    int16_t yTopLeftSearchRegion;
    uint32_t searchRegionIndex;
    // round the search region width to nearest multiple of 8 if it is less than 8 or non multiple of 8
    // SAD calculation performance is the same for searchregion width from 1 to 8
    int16_t search_area_width = (hmeLevel1SearchAreaInWidth < 8) ? 8 : (hmeLevel1SearchAreaInWidth & 7) ? hmeLevel1SearchAreaInWidth + (hmeLevel1SearchAreaInWidth - ((hmeLevel1SearchAreaInWidth >> 3) << 3)) : hmeLevel1SearchAreaInWidth;
    int16_t search_area_height = hmeLevel1SearchAreaInHeight;

    int16_t x_search_area_origin;
    int16_t y_search_area_origin;

    int16_t padWidth = (int16_t)(quarterRefPicPtr->origin_x) - 1;
    int16_t padHeight = (int16_t)(quarterRefPicPtr->origin_y) - 1;

    x_search_area_origin = -(search_area_width >> 1) + xLevel0SearchCenter;
    y_search_area_origin = -(search_area_height >> 1) + yLevel0SearchCenter;

    // Correct the left edge of the Search Area if it is not on the reference Picture
    x_search_area_origin = ((origin_x + x_search_area_origin) < -padWidth) ?
        -padWidth - origin_x :
        x_search_area_origin;

    search_area_width = ((origin_x + x_search_area_origin) < -padWidth) ?
        search_area_width - (-padWidth - (origin_x + x_search_area_origin)) :
        search_area_width;
    // Correct the right edge of the Search Area if its not on the reference Picture
    x_search_area_origin = ((origin_x + x_search_area_origin) > (int16_t)quarterRefPicPtr->width - 1) ?
        x_search_area_origin - ((origin_x + x_search_area_origin) - ((int16_t)quarterRefPicPtr->width - 1)) :
        x_search_area_origin;

    search_area_width = ((origin_x + x_search_area_origin + search_area_width) > (int16_t)quarterRefPicPtr->width) ?
        MAX(1, search_area_width - ((origin_x + x_search_area_origin + search_area_width) - (int16_t)quarterRefPicPtr->width)) :
        search_area_width;

    // Correct the top edge of the Search Area if it is not on the reference Picture
    y_search_area_origin = ((origin_y + y_search_area_origin) < -padHeight) ?
        -padHeight - origin_y :
        y_search_area_origin;

    search_area_height = ((origin_y + y_search_area_origin) < -padHeight) ?
        search_area_height - (-padHeight - (origin_y + y_search_area_origin)) :
        search_area_height;

    // Correct the bottom edge of the Search Area if its not on the reference Picture
    y_search_area_origin = ((origin_y + y_search_area_origin) > (int16_t)quarterRefPicPtr->height - 1) ?
        y_search_area_origin - ((origin_y + y_search_area_origin) - ((int16_t)quarterRefPicPtr->height - 1)) :
        y_search_area_origin;

    search_area_height = (origin_y + y_search_area_origin + search_area_height > (int16_t)quarterRefPicPtr->height) ?
        MAX(1, search_area_height - ((origin_y + y_search_area_origin + search_area_height) - (int16_t)quarterRefPicPtr->height)) :
        search_area_height;

    // Move to the top left of the search region
    xTopLeftSearchRegion = ((int16_t)quarterRefPicPtr->origin_x + origin_x) + x_search_area_origin;
    yTopLeftSearchRegion = ((int16_t)quarterRefPicPtr->origin_y + origin_y) + y_search_area_origin;
    searchRegionIndex = xTopLeftSearchRegion + yTopLeftSearchRegion * quarterRefPicPtr->strideY;

    if (((sb_width & 7) == 0) || (sb_width == 4))
    {
        // Put the first search location into level0 results
        NxMSadLoopKernel_funcPtrArray[asm_type](
            &context_ptr->quarter_sb_buffer[0],
            context_ptr->quarter_sb_buffer_stride * 2,
            &quarterRefPicPtr->bufferY[searchRegionIndex],
            quarterRefPicPtr->strideY * 2,
            sb_height >> 1, sb_width,
            /* results */
            level1BestSad,
            xLevel1SearchCenter,
            yLevel1SearchCenter,
            /* range */
            quarterRefPicPtr->strideY,
            search_area_width,
            search_area_height
            );
    }
    else
    {
        SadLoopKernel(
            &context_ptr->quarter_sb_buffer[0],
            context_ptr->quarter_sb_buffer_stride * 2,
            &quarterRefPicPtr->bufferY[searchRegionIndex],
            quarterRefPicPtr->strideY * 2,
            sb_height >> 1, sb_width,
            /* results */
            level1BestSad,
            xLevel1SearchCenter,
            yLevel1SearchCenter,
            /* range */
            quarterRefPicPtr->strideY,
            search_area_width,
            search_area_height
        );
    }

    *level1BestSad *= 2; // Multiply by 2 because considered only ever other line
    *xLevel1SearchCenter += x_search_area_origin;
    *xLevel1SearchCenter *= 2; // Multiply by 2 because operating on 1/2 resolution
    *yLevel1SearchCenter += y_search_area_origin;
    *yLevel1SearchCenter *= 2; // Multiply by 2 because operating on 1/2 resolution

    return;
}

void HmeLevel2(
    PictureParentControlSet_t   *picture_control_set_ptr,            // input parameter, Picture control set Ptr
    MeContext_t             *context_ptr,                        // input/output parameter, ME context Ptr, used to get/update ME results
    int16_t                   origin_x,                           // input parameter, SB position in the horizontal direction
    int16_t                   origin_y,                           // input parameter, SB position in the vertical direction
    uint32_t                   sb_width,                          // input parameter, SB pwidth - full resolution
    uint32_t                   sb_height,                         // input parameter, SB height - full resolution
    EbPictureBufferDesc_t   *refPicPtr,                         // input parameter, reference Picture Ptr
    uint32_t                   searchRegionNumberInWidth,         // input parameter, search region number in the horizontal direction
    uint32_t                   searchRegionNumberInHeight,        // input parameter, search region number in the vertical direction
    int16_t                   xLevel1SearchCenter,               // input parameter, best Level1 xMV at(searchRegionNumberInWidth, searchRegionNumberInHeight)
    int16_t                   yLevel1SearchCenter,               // input parameter, best Level1 yMV at(searchRegionNumberInWidth, searchRegionNumberInHeight)
    uint64_t                  *level2BestSad,                     // output parameter, Level2 SAD at (searchRegionNumberInWidth, searchRegionNumberInHeight)
    int16_t                  *xLevel2SearchCenter,               // output parameter, Level2 xMV at (searchRegionNumberInWidth, searchRegionNumberInHeight)
    int16_t                  *yLevel2SearchCenter,               // output parameter, Level2 yMV at (searchRegionNumberInWidth, searchRegionNumberInHeight)
    EbAsm                   asm_type)
{


    int16_t xTopLeftSearchRegion;
    int16_t yTopLeftSearchRegion;
    uint32_t searchRegionIndex;

    // round the search region width to nearest multiple of 8 if it is less than 8 or non multiple of 8
    // SAD calculation performance is the same for searchregion width from 1 to 8
#if ME_HME_OQ
    (void)picture_control_set_ptr;
    int16_t hmeLevel2SearchAreaInWidth = (int16_t)context_ptr->hme_level2_search_area_in_width_array[searchRegionNumberInWidth];
#else
    int16_t hmeLevel2SearchAreaInWidth = (int16_t)picture_control_set_ptr->hme_level2_search_area_in_width_array[searchRegionNumberInWidth];
#endif
    int16_t search_area_width = (hmeLevel2SearchAreaInWidth < 8) ? 8 : (hmeLevel2SearchAreaInWidth & 7) ? hmeLevel2SearchAreaInWidth + (hmeLevel2SearchAreaInWidth - ((hmeLevel2SearchAreaInWidth >> 3) << 3)) : hmeLevel2SearchAreaInWidth;
#if ME_HME_OQ
    int16_t search_area_height = (int16_t)context_ptr->hme_level2_search_area_in_height_array[searchRegionNumberInHeight];
#else
    int16_t search_area_height = (int16_t)picture_control_set_ptr->hme_level2_search_area_in_height_array[searchRegionNumberInHeight];
#endif
    int16_t x_search_area_origin;
    int16_t y_search_area_origin;

    int16_t padWidth = (int16_t)BLOCK_SIZE_64 - 1;
    int16_t padHeight = (int16_t)BLOCK_SIZE_64 - 1;

    x_search_area_origin = -(search_area_width >> 1) + xLevel1SearchCenter;
    y_search_area_origin = -(search_area_height >> 1) + yLevel1SearchCenter;

    // Correct the left edge of the Search Area if it is not on the reference Picture
    x_search_area_origin = ((origin_x + x_search_area_origin) < -padWidth) ?
        -padWidth - origin_x :
        x_search_area_origin;

    search_area_width = ((origin_x + x_search_area_origin) < -padWidth) ?
        search_area_width - (-padWidth - (origin_x + x_search_area_origin)) :
        search_area_width;

    // Correct the right edge of the Search Area if its not on the reference Picture
    x_search_area_origin = ((origin_x + x_search_area_origin) > (int16_t)refPicPtr->width - 1) ?
        x_search_area_origin - ((origin_x + x_search_area_origin) - ((int16_t)refPicPtr->width - 1)) :
        x_search_area_origin;

    search_area_width = ((origin_x + x_search_area_origin + search_area_width) > (int16_t)refPicPtr->width) ?
        MAX(1, search_area_width - ((origin_x + x_search_area_origin + search_area_width) - (int16_t)refPicPtr->width)) :
        search_area_width;

    // Correct the top edge of the Search Area if it is not on the reference Picture
    y_search_area_origin = ((origin_y + y_search_area_origin) < -padHeight) ?
        -padHeight - origin_y :
        y_search_area_origin;

    search_area_height = ((origin_y + y_search_area_origin) < -padHeight) ?
        search_area_height - (-padHeight - (origin_y + y_search_area_origin)) :
        search_area_height;

    // Correct the bottom edge of the Search Area if its not on the reference Picture
    y_search_area_origin = ((origin_y + y_search_area_origin) > (int16_t)refPicPtr->height - 1) ?
        y_search_area_origin - ((origin_y + y_search_area_origin) - ((int16_t)refPicPtr->height - 1)) :
        y_search_area_origin;

    search_area_height = (origin_y + y_search_area_origin + search_area_height > (int16_t)refPicPtr->height) ?
        MAX(1, search_area_height - ((origin_y + y_search_area_origin + search_area_height) - (int16_t)refPicPtr->height)) :
        search_area_height;

    // Move to the top left of the search region
    xTopLeftSearchRegion = ((int16_t)refPicPtr->origin_x + origin_x) + x_search_area_origin;
    yTopLeftSearchRegion = ((int16_t)refPicPtr->origin_y + origin_y) + y_search_area_origin;
    searchRegionIndex = xTopLeftSearchRegion + yTopLeftSearchRegion * refPicPtr->strideY;
    if ((((sb_width & 7) == 0) && (sb_width != 40) && (sb_width != 56)))
    {
        // Put the first search location into level0 results
        NxMSadLoopKernel_funcPtrArray[asm_type](
            context_ptr->sb_src_ptr,
            context_ptr->sb_src_stride * 2,
            &refPicPtr->bufferY[searchRegionIndex],
            refPicPtr->strideY * 2,
            sb_height >> 1, sb_width,
            /* results */
            level2BestSad,
            xLevel2SearchCenter,
            yLevel2SearchCenter,
            /* range */
            refPicPtr->strideY,
            search_area_width,
            search_area_height
            );
    }
    else
    {
        // Put the first search location into level0 results
        SadLoopKernel(
            context_ptr->sb_src_ptr,
            context_ptr->sb_src_stride * 2,
            &refPicPtr->bufferY[searchRegionIndex],
            refPicPtr->strideY * 2,
            sb_height >> 1, sb_width,
            /* results */
            level2BestSad,
            xLevel2SearchCenter,
            yLevel2SearchCenter,
            /* range */
            refPicPtr->strideY,
            search_area_width,
            search_area_height
        );

    }

    *level2BestSad *= 2; // Multiply by 2 because considered only every other line
    *xLevel2SearchCenter += x_search_area_origin;
    *yLevel2SearchCenter += y_search_area_origin;

    return;
}



static void SelectBuffer(
    uint32_t                 pu_index,                         //[IN]
    uint8_t                  fracPosition,                    //[IN]
    uint32_t                 pu_width,                         //[IN] Refrence picture list index
    uint32_t                 pu_height,                        //[IN] Refrence picture index in the list
    uint8_t                 *pos_Full,                        //[IN]
    uint8_t                 *pos_b,                           //[IN]
    uint8_t                 *pos_h,                           //[IN]
    uint8_t                 *pos_j,                           //[IN]
    uint32_t                 refHalfStride,                       //[IN]
    uint32_t                 refBufferFullStride,
    uint8_t                 **DstPtr,                             //[OUT]
    uint32_t                 *DstPtrStride,                       //[OUT]
    EbAsm                 asm_type)
{

    (void)asm_type;
    (void)pu_width;
    (void)pu_height;

    uint32_t puShiftXIndex = puSearchIndexMap[pu_index][0];
    uint32_t puShiftYIndex = puSearchIndexMap[pu_index][1];
    uint32_t refStride = refHalfStride;

    //for each one of the 8 positions, we need to determine the 2 buffers to  do averaging
    uint8_t  *buf1 = pos_Full;

    switch (fracPosition)
    {
    case 0: // integer
        buf1 = pos_Full;
        refStride = refBufferFullStride;
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
    default:
        break;
    }

    buf1 = buf1 + puShiftXIndex + puShiftYIndex * refStride;

    *DstPtr = buf1;
    *DstPtrStride = refStride;


    return;
}


static void QuarterPelCompensation(
    uint32_t                 pu_index,                         //[IN]
    uint8_t                  fracPosition,                    //[IN]
    uint32_t                 pu_width,                         //[IN] Refrence picture list index
    uint32_t                 pu_height,                        //[IN] Refrence picture index in the list
    uint8_t                 *pos_Full,                        //[IN]
    uint8_t                 *pos_b,                           //[IN]
    uint8_t                 *pos_h,                           //[IN]
    uint8_t                 *pos_j,                           //[IN]
    uint32_t                 refHalfStride,                       //[IN]
    uint32_t                 refBufferFullStride,
    uint8_t                 *Dst,                             //[IN]
    uint32_t                 DstStride,                       //[IN]
    EbAsm                 asm_type)
{


    uint32_t puShiftXIndex = puSearchIndexMap[pu_index][0];
    uint32_t puShiftYIndex = puSearchIndexMap[pu_index][1];
    uint32_t refStride1 = refHalfStride;
    uint32_t refStride2 = refHalfStride;

    //for each one of the 8 positions, we need to determine the 2 buffers to  do averaging
    uint8_t  *buf1 = pos_Full;
    uint8_t  *buf2 = pos_Full;

    switch (fracPosition)
    {
    case 1: // a
        buf1 = pos_Full;
        buf2 = pos_b;
        refStride1 = refBufferFullStride;
        break;

    case 3: // c
        buf1 = pos_b;
        buf2 = pos_Full + 1;
        refStride2 = refBufferFullStride;
        break;

    case 4: // d
        buf1 = pos_Full;
        buf2 = pos_h;
        refStride1 = refBufferFullStride;
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
        buf1 = pos_h;
        buf2 = pos_Full + refBufferFullStride;
        refStride2 = refBufferFullStride;
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
    default:
        break;
    }



    buf1 = buf1 + puShiftXIndex + puShiftYIndex * refStride1;
    buf2 = buf2 + puShiftXIndex + puShiftYIndex * refStride2;

    picture_average_array[asm_type](buf1, refStride1, buf2, refStride2, Dst, DstStride, pu_width, pu_height);

    return;
}

/*******************************************************************************
* Requirement: pu_width      = 8, 16, 24, 32, 48 or 64
* Requirement: pu_height % 2 = 0
* Requirement: skip         = 0 or 1
* Requirement (x86 only): tempBuf % 16 = 0
* Requirement (x86 only): (dst->bufferY  + dstLumaIndex  ) % 16 = 0 when pu_width %16 = 0
* Requirement (x86 only): (dst->bufferCb + dstChromaIndex) % 16 = 0 when pu_width %32 = 0
* Requirement (x86 only): (dst->bufferCr + dstChromaIndex) % 16 = 0 when pu_width %32 = 0
* Requirement (x86 only): dst->strideY   % 16 = 0 when pu_width %16 = 0
* Requirement (x86 only): dst->chromaStride % 16 = 0 when pu_width %32 = 0
*******************************************************************************/
uint32_t BiPredAverging(
#if M0_SAD_HALF_QUARTER_PEL_BIPRED_SEARCH
    MeContext_t           *context_ptr,
#else
    // MeContext_t           *context_ptr,
#endif
    MePredUnit_t          *me_candidate,
    uint32_t                 pu_index,
    uint8_t                 *sourcePic,
    uint32_t                 lumaStride,
    uint8_t                  firstFracPos,
    uint8_t                  secondFracPos,
    uint32_t                 pu_width,
    uint32_t                 pu_height,
    uint8_t                  *firstRefInteger,
    uint8_t                  *firstRefPosB,
    uint8_t                  *firstRefPosH,
    uint8_t                  *firstRefPosJ,
    uint8_t                  *secondRefInteger,
    uint8_t                  *secondRefPosB,
    uint8_t                  *secondRefPosH,
    uint8_t                  *secondRefPosJ,
    uint32_t                 refBufferStride,
    uint32_t                 refBufferFullList0Stride,
    uint32_t                 refBufferFullList1Stride,
    uint8_t                  *firstRefTempDst,
    uint8_t                  *secondRefTempDst,
    EbAsm                 asm_type)
{

    uint8_t                  *ptrList0, *ptrList1;
    uint32_t                  ptrList0Stride, ptrList1Stride;

    // Buffer Selection and quater-pel compensation on the fly
    if (subPositionType[firstFracPos] != 2) {

        SelectBuffer(
            pu_index,
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

    }
    else {

        QuarterPelCompensation(
            pu_index,
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
            BLOCK_SIZE_64,
            asm_type);

        ptrList0 = firstRefTempDst;
        ptrList0Stride = BLOCK_SIZE_64;
    }

    if (subPositionType[secondFracPos] != 2) {

        SelectBuffer(
            pu_index,
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

    }
    else {
        //uni-prediction List1 luma
        //doing the luma interpolation
        QuarterPelCompensation(
            pu_index,
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
            BLOCK_SIZE_64,
            asm_type);

        ptrList1 = secondRefTempDst;
        ptrList1Stride = BLOCK_SIZE_64;

    }

    // bi-pred luma
#if M0_SAD_HALF_QUARTER_PEL_BIPRED_SEARCH
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
    me_candidate->distortion = (context_ptr->fractionalSearchMethod == SUB_SAD_SEARCH) ?
#else
    me_candidate->distortion = (context_ptr->useSubSadFracBipredSearch) ?
#endif
        NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](sourcePic,
            lumaStride << 1,
            ptrList0,
            ptrList0Stride << 1,
            ptrList1,
            ptrList1Stride << 1,
            pu_height >> 1,
            pu_width) << 1 :
    NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](sourcePic,
        lumaStride,
        ptrList0,
        ptrList0Stride,
        ptrList1,
        ptrList1Stride,
        pu_height,
        pu_width);
#else
    me_candidate->distortion = NxMSadAveragingKernel_funcPtrArray[asm_type][pu_width >> 3](sourcePic,
        lumaStride << 1,
        ptrList0,
        ptrList0Stride << 1,
        ptrList1,
        ptrList1Stride << 1,
        pu_height >> 1,
        pu_width) << 1;
#endif
    return me_candidate->distortion;
}


/*******************************************
* BiPredictionComponsation
*   performs componsation fro List 0 and
*   List1 Candidates and then compute the
*   average
*******************************************/
EbErrorType  BiPredictionCompensation(
    MeContext_t             *context_ptr,
    uint32_t                  pu_index,
    MePredUnit_t            *me_candidate,
    uint32_t                  firstList,
    uint32_t                  firstRefMv,
    uint32_t                  secondList,
    uint32_t                  secondRefMv,
    EbAsm                  asm_type)
{
    EbErrorType                    return_error = EB_ErrorNone;

    int16_t                           firstRefPosX;
    int16_t                           firstRefPosY;
    int16_t                           firstRefIntegPosx;
    int16_t                           firstRefIntegPosy;
    uint8_t                            firstRefFracPosx;
    uint8_t                            firstRefFracPosy;
    uint8_t                            firstRefFracPos;
    int32_t                           xfirstSearchIndex;
    int32_t                           yfirstSearchIndex;
    int32_t                           firstSearchRegionIndexPosInteg;
    int32_t                           firstSearchRegionIndexPosb;
    int32_t                           firstSearchRegionIndexPosh;
    int32_t                           firstSearchRegionIndexPosj;

    int16_t                           secondRefPosX;
    int16_t                           secondRefPosY;
    int16_t                           secondRefIntegPosx;
    int16_t                           secondRefIntegPosy;
    uint8_t                            secondRefFracPosx;
    uint8_t                            secondRefFracPosy;
    uint8_t                            secondRefFracPos;
    int32_t                           xsecondSearchIndex;
    int32_t                           ysecondSearchIndex;
    int32_t                           secondSearchRegionIndexPosInteg;
    int32_t                           secondSearchRegionIndexPosb;
    int32_t                           secondSearchRegionIndexPosh;
    int32_t                           secondSearchRegionIndexPosj;


    uint32_t puShiftXIndex = puSearchIndexMap[pu_index][0];
    uint32_t puShiftYIndex = puSearchIndexMap[pu_index][1];

    const uint32_t puLcuBufferIndex = puShiftXIndex + puShiftYIndex * context_ptr->sb_src_stride;

    me_candidate->prediction_direction = BI_PRED;


    // First refrence
    // Set Candidate information

    firstRefPosX = _MVXT(firstRefMv),
        firstRefPosY = _MVYT(firstRefMv),
        me_candidate->mv[0] = firstRefMv;

    firstRefIntegPosx = (firstRefPosX >> 2);
    firstRefIntegPosy = (firstRefPosY >> 2);
    firstRefFracPosx = (uint8_t)firstRefPosX & 0x03;
    firstRefFracPosy = (uint8_t)firstRefPosY & 0x03;

    firstRefFracPos = (uint8_t)(firstRefFracPosx + (firstRefFracPosy << 2));

    xfirstSearchIndex = (int32_t)firstRefIntegPosx - context_ptr->x_search_area_origin[firstList][0];
    yfirstSearchIndex = (int32_t)firstRefIntegPosy - context_ptr->y_search_area_origin[firstList][0];
    firstSearchRegionIndexPosInteg = (int32_t)(xfirstSearchIndex + (ME_FILTER_TAP >> 1)) + (int32_t)context_ptr->interpolated_full_stride[firstList][0] * (int32_t)(yfirstSearchIndex + (ME_FILTER_TAP >> 1));
    firstSearchRegionIndexPosb = (int32_t)(xfirstSearchIndex + (ME_FILTER_TAP >> 1) - 1) + (int32_t)context_ptr->interpolated_stride * (int32_t)(yfirstSearchIndex + (ME_FILTER_TAP >> 1));
    firstSearchRegionIndexPosh = (int32_t)(xfirstSearchIndex + (ME_FILTER_TAP >> 1) - 1) + (int32_t)context_ptr->interpolated_stride * (int32_t)(yfirstSearchIndex + (ME_FILTER_TAP >> 1) - 1);
    firstSearchRegionIndexPosj = (int32_t)(xfirstSearchIndex + (ME_FILTER_TAP >> 1) - 1) + (int32_t)context_ptr->interpolated_stride * (int32_t)(yfirstSearchIndex + (ME_FILTER_TAP >> 1) - 1);

    // Second refrence

    // Set Candidate information
    secondRefPosX = _MVXT(secondRefMv),
        secondRefPosY = _MVYT(secondRefMv),
        me_candidate->mv[1] = secondRefMv;

    secondRefIntegPosx = (secondRefPosX >> 2);
    secondRefIntegPosy = (secondRefPosY >> 2);
    secondRefFracPosx = (uint8_t)secondRefPosX & 0x03;
    secondRefFracPosy = (uint8_t)secondRefPosY & 0x03;

    secondRefFracPos = (uint8_t)(secondRefFracPosx + (secondRefFracPosy << 2));

    xsecondSearchIndex = secondRefIntegPosx - context_ptr->x_search_area_origin[secondList][0];
    ysecondSearchIndex = secondRefIntegPosy - context_ptr->y_search_area_origin[secondList][0];
    secondSearchRegionIndexPosInteg = (int32_t)(xsecondSearchIndex + (ME_FILTER_TAP >> 1)) + (int32_t)context_ptr->interpolated_full_stride[secondList][0] * (int32_t)(ysecondSearchIndex + (ME_FILTER_TAP >> 1));
    secondSearchRegionIndexPosb = (int32_t)(xsecondSearchIndex + (ME_FILTER_TAP >> 1) - 1) + (int32_t)context_ptr->interpolated_stride * (int32_t)(ysecondSearchIndex + (ME_FILTER_TAP >> 1));
    secondSearchRegionIndexPosh = (int32_t)(xsecondSearchIndex + (ME_FILTER_TAP >> 1) - 1) + (int32_t)context_ptr->interpolated_stride * (int32_t)(ysecondSearchIndex + (ME_FILTER_TAP >> 1) - 1);
    secondSearchRegionIndexPosj = (int32_t)(xsecondSearchIndex + (ME_FILTER_TAP >> 1) - 1) + (int32_t)context_ptr->interpolated_stride * (int32_t)(ysecondSearchIndex + (ME_FILTER_TAP >> 1) - 1);

    uint32_t nIndex;

    if (pu_index > 200) {
        nIndex = pu_index;
    }
    else if (pu_index > 184) {

        nIndex = tab8x32[pu_index - 185] + 185;
    }
    else if (pu_index > 168) {
        nIndex = tab32x8[pu_index - 169] + 169;
    }
    else if (pu_index > 136) {
        nIndex = tab8x16[pu_index - 137] + 137;
    }
    else if (pu_index > 128) {
        nIndex = tab16x32[pu_index - 129] + 129;
    }
    else if (pu_index > 126) {
        nIndex = pu_index;
    }
    else if (pu_index > 94) {
        nIndex = tab16x8[pu_index - 95] + 95;
    }
    else if (pu_index > 86) {
        nIndex = tab32x16[pu_index - 87] + 87;
    }
    else if (pu_index > 84) {
        nIndex = pu_index;
    }
    else if (pu_index > 20) {
        nIndex = tab8x8[pu_index - 21] + 21;
    }
    else if (pu_index > 4) {
        nIndex = tab16x16[pu_index - 5] + 5;
    }
    else {
        nIndex = pu_index;
    }


    context_ptr->p_sb_bipred_sad[nIndex] =

        BiPredAverging(
#if M0_SAD_HALF_QUARTER_PEL_BIPRED_SEARCH
            context_ptr,
#else
            // context_ptr,
#endif
            me_candidate,
            pu_index,
            &(context_ptr->sb_src_ptr[puLcuBufferIndex]),
            context_ptr->sb_src_stride,
            firstRefFracPos,
            secondRefFracPos,
            partitionWidth[pu_index],
            partitionHeight[pu_index],
            &(context_ptr->integer_buffer_ptr[firstList][0][firstSearchRegionIndexPosInteg]),
            &(context_ptr->pos_b_buffer[firstList][0][firstSearchRegionIndexPosb]),
            &(context_ptr->pos_h_buffer[firstList][0][firstSearchRegionIndexPosh]),
            &(context_ptr->pos_j_buffer[firstList][0][firstSearchRegionIndexPosj]),
            &(context_ptr->integer_buffer_ptr[secondList][0][secondSearchRegionIndexPosInteg]),
            &(context_ptr->pos_b_buffer[secondList][0][secondSearchRegionIndexPosb]),
            &(context_ptr->pos_h_buffer[secondList][0][secondSearchRegionIndexPosh]),
            &(context_ptr->pos_j_buffer[secondList][0][secondSearchRegionIndexPosj]),
            context_ptr->interpolated_stride,
            context_ptr->interpolated_full_stride[firstList][0],
            context_ptr->interpolated_full_stride[secondList][0],
            &(context_ptr->one_d_intermediate_results_buf0[0]),
            &(context_ptr->one_d_intermediate_results_buf1[0]),
            asm_type);

    return return_error;
}

/*******************************************
* BiPredictionSearch
*   performs Bi-Prediction Search (LCU)
*******************************************/
// This function enables all 16 Bipred candidates when MRP is ON
EbErrorType  BiPredictionSearch(
    MeContext_t                    *context_ptr,
    uint32_t                        pu_index,
    uint8_t                        candidateIndex,
    uint32_t                        activeRefPicFirstLisNum,
    uint32_t                        activeRefPicSecondLisNum,
    uint8_t                        *totalMeCandidateIndex,
    EbAsm                        asm_type,
    PictureParentControlSet_t   *picture_control_set_ptr)
{
    EbErrorType                 return_error = EB_ErrorNone;

    uint32_t firstListRefPictdx;
    uint32_t secondListRefPictdx;

    (void)picture_control_set_ptr;

    uint32_t nIndex;

    if (pu_index > 200) {
        nIndex = pu_index;
    }
    else if (pu_index > 184) {

        nIndex = tab8x32[pu_index - 185] + 185;
    }
    else if (pu_index > 168) {
        nIndex = tab32x8[pu_index - 169] + 169;
    }
    else if (pu_index > 136) {
        nIndex = tab8x16[pu_index - 137] + 137;
    }
    else if (pu_index > 128) {
        nIndex = tab16x32[pu_index - 129] + 129;
    }
    else if (pu_index > 126) {
        nIndex = pu_index;
    }
    else if (pu_index > 94) {
        nIndex = tab16x8[pu_index - 95] + 95;
    }
    else if (pu_index > 86) {
        nIndex = tab32x16[pu_index - 87] + 87;
    }
    else if (pu_index > 84) {
        nIndex = pu_index;
    }
    else if (pu_index > 20) {
        nIndex = tab8x8[pu_index - 21] + 21;
    }
    else if (pu_index > 4) {
        nIndex = tab16x16[pu_index - 5] + 5;
    }
    else {
        nIndex = pu_index;
    }

    for (firstListRefPictdx = 0; firstListRefPictdx < activeRefPicFirstLisNum; firstListRefPictdx++) {
        for (secondListRefPictdx = 0; secondListRefPictdx < activeRefPicSecondLisNum; secondListRefPictdx++) {

            {

                BiPredictionCompensation(
                    context_ptr,
                    pu_index,
                    &(context_ptr->me_candidate[candidateIndex].pu[pu_index]),
                    REFERENCE_PIC_LIST_0,
                    context_ptr->p_sb_best_mv[REFERENCE_PIC_LIST_0][firstListRefPictdx][nIndex],
                    REFERENCE_PIC_LIST_1,
                    context_ptr->p_sb_best_mv[REFERENCE_PIC_LIST_1][secondListRefPictdx][nIndex],
                    asm_type);

                candidateIndex++;

            }
        }
    }

    *totalMeCandidateIndex = candidateIndex;

    return return_error;
}

// Nader - to be replaced by loock-up table
/*******************************************
* get_me_info_index
*   search the correct index of the motion
*   info that corresponds to the input
*   md candidate
*******************************************/
uint32_t get_me_info_index(
    uint32_t                  max_me_block,
    const BlockGeom        *blk_geom,
    uint32_t               geom_offset_x,
    uint32_t                 geom_offset_y

)
{
    // search for motion info
    uint32_t                  block_index;
    uint32_t                  me_info_index = 0xFFFFFFF;


    for (block_index = 0; block_index < max_me_block; block_index++) {

        if ((blk_geom->bwidth == partitionWidth[block_index]) &&
            (blk_geom->bheight == partitionHeight[block_index]) &&
            ((blk_geom->origin_x - geom_offset_x) == puSearchIndexMap[block_index][0]) &&
            ((blk_geom->origin_y - geom_offset_y) == puSearchIndexMap[block_index][1])) {


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
uint32_t get_in_loop_me_info_index(
    uint32_t                  max_me_block,
    uint8_t                   is_128_sb,
    const BlockGeom        *blk_geom
)
{
    // search for motion info
    uint32_t                  block_index;
    uint32_t                  me_info_index = 0xFFFFFFF;
    if (is_128_sb) {

        for (block_index = 0; block_index < max_me_block; block_index++) {

            if (blk_geom->bwidth == in_loop_me_block_width_128_sb[block_index] && blk_geom->bheight == in_loop_me_block_height_128_sb[block_index] &&
                blk_geom->origin_x == in_loop_me_block_index_128_sb[block_index][0] && blk_geom->origin_y == in_loop_me_block_index_128_sb[block_index][1]) {

                me_info_index = block_index;
                break;

            }

        }

    }
    else {
        for (block_index = 0; block_index < max_me_block; block_index++) {

            if (blk_geom->bwidth == in_loop_me_block_width[block_index] && blk_geom->bheight == in_loop_me_block_height[block_index] &&
                blk_geom->origin_x == in_loop_me_block_index[block_index][0] && blk_geom->origin_y == in_loop_me_block_index[block_index][1]) {

                me_info_index = block_index;
                break;

            }

        }
    }

    return me_info_index;
}

#define NSET_CAND(mePuResult, num, dist, dir) \
     (mePuResult)->distortionDirection[(num)].distortion = (dist); \
     (mePuResult)->distortionDirection[(num)].direction = (dir)  ;


int8_t Sort3Elements(uint32_t a, uint32_t b, uint32_t c) {

    uint8_t sortCode = 0;
    if (a <= b && a <= c) {
        if (b <= c) {
            sortCode = a_b_c;
        }
        else {
            sortCode = a_c_b;
        }
    }
    else if (b <= a && b <= c) {
        if (a <= c) {
            sortCode = b_a_c;
        }
        else {
            sortCode = b_c_a;
        }

    }
    else if (a <= b) {
        sortCode = c_a_b;
    }
    else {
        sortCode = c_b_a;
    }


    return sortCode;
}


EbErrorType CheckZeroZeroCenter(
    EbPictureBufferDesc_t        *refPicPtr,
    MeContext_t                  *context_ptr,
    uint32_t                       sb_origin_x,
    uint32_t                       sb_origin_y,
    uint32_t                       sb_width,
    uint32_t                       sb_height,
    int16_t                       *xSearchCenter,
    int16_t                       *ySearchCenter,
    EbAsm                       asm_type)

{
    EbErrorType return_error = EB_ErrorNone;
    uint32_t       searchRegionIndex, zeroMvSad, hmeMvSad, hmeMvdRate;
    uint64_t       hmeMvCost, zeroMvCost, searchCenterCost;
    int16_t        origin_x = (int16_t)sb_origin_x;
    int16_t        origin_y = (int16_t)sb_origin_y;
    uint32_t       subsampleSad = 1;
#if HME_ENHANCED_CENTER_SEARCH
    int16_t        pad_width = (int16_t)BLOCK_SIZE_64 - 1;
    int16_t        pad_height = (int16_t)BLOCK_SIZE_64 - 1;
#endif
    searchRegionIndex = (int16_t)refPicPtr->origin_x + origin_x +
        ((int16_t)refPicPtr->origin_y + origin_y) * refPicPtr->strideY;

    zeroMvSad = NxMSadKernel_funcPtrArray[asm_type][sb_width >> 3](
        context_ptr->sb_src_ptr,
        context_ptr->sb_src_stride << subsampleSad,
        &(refPicPtr->bufferY[searchRegionIndex]),
        refPicPtr->strideY << subsampleSad,
        sb_height >> subsampleSad,
        sb_width);

    zeroMvSad = zeroMvSad << subsampleSad;

#if HME_ENHANCED_CENTER_SEARCH
    // FIX
    // Correct the left edge of the Search Area if it is not on the reference Picture
    *xSearchCenter = ((origin_x + *xSearchCenter) < -pad_width) ?
        -pad_width - origin_x :
        *xSearchCenter;
    // Correct the right edge of the Search Area if its not on the reference Picture
    *xSearchCenter = ((origin_x + *xSearchCenter) > (int16_t)refPicPtr->width - 1) ?
        *xSearchCenter - ((origin_x + *xSearchCenter) - ((int16_t)refPicPtr->width - 1)) :
        *xSearchCenter;
    // Correct the top edge of the Search Area if it is not on the reference Picture
    *ySearchCenter = ((origin_y + *ySearchCenter) < -pad_height) ?
        -pad_height - origin_y :
        *ySearchCenter;
    // Correct the bottom edge of the Search Area if its not on the reference Picture
    *ySearchCenter = ((origin_y + *ySearchCenter) > (int16_t)refPicPtr->height - 1) ?
        *ySearchCenter - ((origin_y + *ySearchCenter) - ((int16_t)refPicPtr->height - 1)) :
        *ySearchCenter;
    ///
#endif

    zeroMvCost = zeroMvSad << COST_PRECISION;
    searchRegionIndex = (int16_t)(refPicPtr->origin_x + origin_x) + *xSearchCenter +
        ((int16_t)(refPicPtr->origin_y + origin_y) + *ySearchCenter) * refPicPtr->strideY;

    hmeMvSad = NxMSadKernel_funcPtrArray[asm_type][sb_width >> 3](
        context_ptr->sb_src_ptr,
        context_ptr->sb_src_stride << subsampleSad,
        &(refPicPtr->bufferY[searchRegionIndex]),
        refPicPtr->strideY << subsampleSad,
        sb_height >> subsampleSad,
        sb_width);

    hmeMvSad = hmeMvSad << subsampleSad;


    hmeMvdRate = 0;
    // AMIR use AV1 rate estimation functions
    //MeGetMvdFractionBits(
    //    ABS(*xSearchCenter << 2),
    //    ABS(*ySearchCenter << 2),
    //    context_ptr->mvd_bits_array,
    //    &hmeMvdRate);

    hmeMvCost = (hmeMvSad << COST_PRECISION) + (((context_ptr->lambda * hmeMvdRate) + MD_OFFSET) >> MD_SHIFT);
    searchCenterCost = MIN(zeroMvCost, hmeMvCost);

    *xSearchCenter = (searchCenterCost == zeroMvCost) ? 0 : *xSearchCenter;
    *ySearchCenter = (searchCenterCost == zeroMvCost) ? 0 : *ySearchCenter;

    return return_error;
}

EbErrorType     suPelEnable(
    MeContext_t                 *context_ptr,
    PictureParentControlSet_t   *picture_control_set_ptr,
    uint32_t listIndex,
    uint32_t refPicIndex,
    EbBool *enableHalfPel32x32,
    EbBool *enableHalfPel16x16,
    EbBool *enableHalfPel8x8)
{
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


    avgMvx32x32 = (_MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x32_0]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x32_1]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x32_2]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x32_3])) >> 2;
    avgMvy32x32 = (_MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x32_0]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x32_1]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x32_2]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x32_3])) >> 2;
    mvMag32x32 = SQR(avgMvx32x32) + SQR(avgMvy32x32);

    avgMvx16x16 = (_MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_0]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_1]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_2]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_3])
        + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_4]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_5]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_6]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_7])
        + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_8]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_9]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_10]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_11])
        + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_12]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_13]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_14]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_15])
        ) >> 4;
    avgMvy16x16 = (_MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_0]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_1]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_2]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_3])
        + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_4]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_5]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_6]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_7])
        + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_8]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_9]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_10]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_11])
        + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_12]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_13]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_14]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_15])
        ) >> 4;
    mvMag16x16 = SQR(avgMvx16x16) + SQR(avgMvy16x16);

    avgMvx8x8 = (_MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_0]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_1]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_2]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_3])
        + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_4]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_5]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_6]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_7])
        + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_8]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_9]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_10]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_11])
        + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_12]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_13]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_14]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_15])
        + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_16]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_17]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_18]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_19])
        + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_20]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_21]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_22]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_23])
        + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_24]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_25]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_26]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_27])
        + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_28]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_29]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_30]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_31])
        + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_32]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_33]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_34]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_35])
        + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_36]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_37]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_38]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_39])
        + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_40]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_41]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_42]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_43])
        + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_44]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_45]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_46]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_47])
        + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_48]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_49]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_50]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_51])
        + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_52]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_53]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_54]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_55])
        + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_56]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_57]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_58]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_59])
        + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_60]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_61]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_62]) + _MVXT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_63])
        ) >> 6;
    avgMvy8x8 = (_MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_0]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_1]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_2]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_3])
        + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_4]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_5]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_6]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_7])
        + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_8]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_9]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_10]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_11])
        + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_12]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_13]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_14]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_15])
        + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_16]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_17]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_18]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_19])
        + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_20]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_21]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_22]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_23])
        + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_24]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_25]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_26]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_27])
        + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_28]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_29]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_30]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_31])
        + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_32]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_33]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_34]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_35])
        + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_36]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_37]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_38]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_39])
        + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_40]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_41]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_42]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_43])
        + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_44]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_45]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_46]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_47])
        + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_48]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_49]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_50]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_51])
        + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_52]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_53]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_54]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_55])
        + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_56]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_57]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_58]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_59])
        + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_60]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_61]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_62]) + _MVYT(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_63])
        ) >> 6;
    mvMag8x8 = SQR(avgMvx8x8) + SQR(avgMvy8x8);

    avgSad32x32 = (context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x32_0] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x32_1] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x32_2] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x32_3]) >> 2;
    avgSad16x16 = (context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_0] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_1] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_2] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_3]
        + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_4] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_5] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_6] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_7]
        + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_8] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_9] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_10] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_11]
        + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_12] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_13] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_14] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_15]
        ) >> 4;
    avgSad8x8 = (context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_0] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_1] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_2] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_3]
        + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_4] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_5] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_6] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_7]
        + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_8] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_9] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_10] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_11]
        + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_12] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_13] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_14] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_15]
        + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_16] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_17] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_18] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_19]
        + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_20] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_21] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_22] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_23]
        + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_24] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_25] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_26] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_27]
        + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_28] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_29] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_30] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_31]
        + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_32] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_33] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_34] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_35]
        + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_36] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_37] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_38] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_39]
        + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_40] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_41] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_42] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_43]
        + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_44] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_45] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_46] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_47]
        + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_48] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_49] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_50] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_51]
        + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_52] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_53] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_54] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_55]
        + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_56] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_57] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_58] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_59]
        + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_60] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_61] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_62] + context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_63]
        ) >> 6;


    if (picture_control_set_ptr->temporal_layer_index == 0)
    {
        //32x32
        if ((mvMag32x32 < SQR(48)) && (avgSad32x32 < 32 * 32 * 6))
        {
            *enableHalfPel32x32 = EB_TRUE; //CLASS_0
        }
        else if ((mvMag32x32 < SQR(48)) && !(avgSad32x32 < 32 * 32 * 6))
        {
            *enableHalfPel32x32 = EB_FALSE; //CLASS_1
        }
        else if (!(mvMag32x32 < SQR(48)) && (avgSad32x32 < 32 * 32 * 6))
        {
            *enableHalfPel32x32 = EB_TRUE; //CLASS_2
        }
        else
        {
            *enableHalfPel32x32 = EB_FALSE; //CLASS_3
        }
        //16x16
        if ((mvMag16x16 < SQR(48)) && (avgSad16x16 < 16 * 16 * 2))
        {
            *enableHalfPel16x16 = EB_FALSE; //CLASS_0
        }
        else if ((mvMag16x16 < SQR(48)) && !(avgSad16x16 < 16 * 16 * 2))
        {
            *enableHalfPel16x16 = EB_TRUE; //CLASS_1
        }
        else if (!(mvMag16x16 < SQR(48)) && (avgSad16x16 < 16 * 16 * 2))
        {
            *enableHalfPel16x16 = EB_FALSE; //CLASS_2
        }
        else
        {
            *enableHalfPel16x16 = EB_TRUE; //CLASS_3
        }
        //8x8
        if ((mvMag8x8 < SQR(48)) && (avgSad8x8 < 8 * 8 * 2))
        {
            *enableHalfPel8x8 = EB_FALSE; //CLASS_0
        }
        else if ((mvMag8x8 < SQR(48)) && !(avgSad8x8 < 8 * 8 * 2))
        {
            *enableHalfPel8x8 = EB_TRUE; //CLASS_1
        }
        else if (!(mvMag8x8 < SQR(48)) && (avgSad8x8 < 8 * 8 * 2))
        {
            *enableHalfPel8x8 = EB_FALSE; //CLASS_2
        }
        else
        {
            *enableHalfPel8x8 = EB_TRUE; //CLASS_3
        }

    }

    else if (picture_control_set_ptr->temporal_layer_index == 1)
    {
        //32x32
        if ((mvMag32x32 < SQR(32)) && (avgSad32x32 < 32 * 32 * 6))
        {
            *enableHalfPel32x32 = EB_TRUE; //CLASS_0
        }
        else if ((mvMag32x32 < SQR(32)) && !(avgSad32x32 < 32 * 32 * 6))
        {
            *enableHalfPel32x32 = EB_FALSE; //CLASS_1
        }
        else if (!(mvMag32x32 < SQR(32)) && (avgSad32x32 < 32 * 32 * 6))
        {
            *enableHalfPel32x32 = EB_TRUE; //CLASS_2
        }
        else
        {
            *enableHalfPel32x32 = EB_TRUE; //CLASS_3
        }
        //16x16
        if ((mvMag16x16 < SQR(32)) && (avgSad16x16 < 16 * 16 * 2))
        {
            *enableHalfPel16x16 = EB_FALSE; //CLASS_0
        }
        else if ((mvMag16x16 < SQR(32)) && !(avgSad16x16 < 16 * 16 * 2))
        {
            *enableHalfPel16x16 = EB_TRUE; //CLASS_1
        }
        else if (!(mvMag16x16 < SQR(32)) && (avgSad16x16 < 16 * 16 * 2))
        {
            *enableHalfPel16x16 = EB_FALSE; //CLASS_2
        }
        else
        {
            *enableHalfPel16x16 = EB_TRUE; //CLASS_3
        }
        //8x8
        if ((mvMag8x8 < SQR(32)) && (avgSad8x8 < 8 * 8 * 2))
        {
            *enableHalfPel8x8 = EB_FALSE; //CLASS_0
        }
        else if ((mvMag8x8 < SQR(32)) && !(avgSad8x8 < 8 * 8 * 2))
        {
            *enableHalfPel8x8 = EB_TRUE; //CLASS_1
        }
        else if (!(mvMag8x8 < SQR(32)) && (avgSad8x8 < 8 * 8 * 2))
        {
            *enableHalfPel8x8 = EB_FALSE; //CLASS_2
        }
        else
        {
            *enableHalfPel8x8 = EB_TRUE; //CLASS_3
        }

    }
    else if (picture_control_set_ptr->temporal_layer_index == 2)
    {
        //32x32
        if ((mvMag32x32 < SQR(80)) && (avgSad32x32 < 32 * 32 * 6))
        {
            *enableHalfPel32x32 = EB_TRUE; //CLASS_0
        }
        else if ((mvMag32x32 < SQR(80)) && !(avgSad32x32 < 32 * 32 * 6))
        {
            *enableHalfPel32x32 = EB_FALSE; //CLASS_1

        }
        else if (!(mvMag32x32 < SQR(80)) && (avgSad32x32 < 32 * 32 * 6))
        {
            *enableHalfPel32x32 = EB_TRUE; //CLASS_2
        }
        else
        {
            *enableHalfPel32x32 = EB_FALSE; //CLASS_3
        }
        //16x16
        if ((mvMag16x16 < SQR(80)) && (avgSad16x16 < 16 * 16 * 2))
        {
            *enableHalfPel16x16 = EB_FALSE; //CLASS_0
        }
        else if ((mvMag16x16 < SQR(80)) && !(avgSad16x16 < 16 * 16 * 2))
        {
            *enableHalfPel16x16 = EB_TRUE; //CLASS_1
        }
        else if (!(mvMag16x16 < SQR(80)) && (avgSad16x16 < 16 * 16 * 2))
        {
            *enableHalfPel16x16 = EB_FALSE; //CLASS_2
        }
        else
        {
            *enableHalfPel16x16 = EB_TRUE; //CLASS_3
        }
        //8x8
        if ((mvMag8x8 < SQR(80)) && (avgSad8x8 < 8 * 8 * 2))
        {
            *enableHalfPel8x8 = EB_FALSE; //CLASS_0
        }
        else if ((mvMag8x8 < SQR(80)) && !(avgSad8x8 < 8 * 8 * 2))
        {
            *enableHalfPel8x8 = EB_TRUE; //CLASS_1
        }
        else if (!(mvMag8x8 < SQR(80)) && (avgSad8x8 < 8 * 8 * 2))
        {
            *enableHalfPel8x8 = EB_FALSE; //CLASS_2
        }
        else
        {
            *enableHalfPel8x8 = EB_TRUE; //CLASS_3
        }

    }
    else
    {
        //32x32
        if ((mvMag32x32 < SQR(48)) && (avgSad32x32 < 32 * 32 * 6))
        {
            *enableHalfPel32x32 = EB_TRUE; //CLASS_0
        }
        else if ((mvMag32x32 < SQR(48)) && !(avgSad32x32 < 32 * 32 * 6))
        {
            *enableHalfPel32x32 = EB_TRUE; //CLASS_1
        }
        else if (!(mvMag32x32 < SQR(48)) && (avgSad32x32 < 32 * 32 * 6))
        {
            *enableHalfPel32x32 = EB_TRUE; //CLASS_2
        }
        else
        {
            *enableHalfPel32x32 = EB_FALSE; //CLASS_3
        }
        //16x16
        if ((mvMag16x16 < SQR(48)) && (avgSad16x16 < 16 * 16 * 2))
        {
            *enableHalfPel16x16 = EB_FALSE; //CLASS_0
        }
        else if ((mvMag16x16 < SQR(48)) && !(avgSad16x16 < 16 * 16 * 2))
        {
            *enableHalfPel16x16 = EB_TRUE; //CLASS_1
        }
        else if (!(mvMag16x16 < SQR(48)) && (avgSad16x16 < 16 * 16 * 2))
        {
            *enableHalfPel16x16 = EB_FALSE; //CLASS_2
        }
        else
        {
            *enableHalfPel16x16 = EB_TRUE; //CLASS_3
        }
        //8x8
        if ((mvMag8x8 < SQR(48)) && (avgSad8x8 < 8 * 8 * 2))
        {
            *enableHalfPel8x8 = EB_FALSE; //CLASS_0
        }
        else if ((mvMag8x8 < SQR(48)) && !(avgSad8x8 < 8 * 8 * 2))
        {
            *enableHalfPel8x8 = EB_TRUE; //CLASS_1
        }
        else if (!(mvMag8x8 < SQR(48)) && (avgSad8x8 < 8 * 8 * 2))
        {
            *enableHalfPel8x8 = EB_FALSE; //CLASS_2
        }
        else
        {
            *enableHalfPel8x8 = EB_FALSE;// EB_TRUE; //CLASS_3
        }

    }

    return return_error;
}
#if HME_ENHANCED_CENTER_SEARCH
static void hme_mv_center_check(
    EbPictureBufferDesc_t       *ref_pic_ptr,
    MeContext_t                    *context_ptr,
    int16_t                     *xsc,
    int16_t                     *ysc,
    uint32_t                     list_index,
    int16_t                      origin_x,
    int16_t                      origin_y,
    uint32_t                     sb_width,
    uint32_t                     sb_height,
    EbAsm                        asm_type)
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
    uint8_t  sparce_scale = 1;
    int16_t pad_width = (int16_t)BLOCK_SIZE_64 - 1;
    int16_t pad_height = (int16_t)BLOCK_SIZE_64 - 1;
    // O pos

    search_region_index = (int16_t)ref_pic_ptr->origin_x + origin_x +
        ((int16_t)ref_pic_ptr->origin_y + origin_y) * ref_pic_ptr->strideY;

    uint32_t sub_sampled_sad = 1;
    uint64_t zero_mv_sad = NxMSadKernel_funcPtrArray[asm_type][sb_width >> 3](
        context_ptr->sb_src_ptr,
        context_ptr->sb_src_stride << sub_sampled_sad,
        &(ref_pic_ptr->bufferY[search_region_index]),
        ref_pic_ptr->strideY << sub_sampled_sad,
        sb_height >> sub_sampled_sad,
        sb_width);

    zero_mv_sad = zero_mv_sad << sub_sampled_sad;

    uint64_t zero_mv_cost = zero_mv_sad << COST_PRECISION;

    //A pos
    search_center_x = 0 - (context_ptr->hme_level0_total_search_area_width*sparce_scale);
    search_center_y = 0;

    // Correct the left edge of the Search Area if it is not on the reference Picture
    search_center_x = ((origin_x + search_center_x) < -pad_width) ?
        -pad_width - origin_x :
        search_center_x;
    // Correct the right edge of the Search Area if its not on the reference Picture
    search_center_x = ((origin_x + search_center_x) > (int16_t)ref_pic_ptr->width - 1) ?
        search_center_x - ((origin_x + search_center_x) - ((int16_t)ref_pic_ptr->width - 1)) :
        search_center_x;
    // Correct the top edge of the Search Area if it is not on the reference Picture
    search_center_y = ((origin_y + search_center_y) < -pad_height) ?
        -pad_height - origin_y :
        search_center_y;

    // Correct the bottom edge of the Search Area if its not on the reference Picture
    search_center_y = ((origin_y + search_center_y) > (int16_t)ref_pic_ptr->height - 1) ?
        search_center_y - ((origin_y + search_center_y) - ((int16_t)ref_pic_ptr->height - 1)) :
        search_center_y;


    uint64_t mv_a_sad = NxMSadKernel_funcPtrArray[asm_type][sb_width >> 3](
        context_ptr->sb_src_ptr,
        context_ptr->sb_src_stride << sub_sampled_sad,
        &(ref_pic_ptr->bufferY[search_region_index]),
        ref_pic_ptr->strideY << sub_sampled_sad,
        sb_height >> sub_sampled_sad,
        sb_width);

    mv_a_sad = mv_a_sad << sub_sampled_sad;

    uint64_t mv_a_cost = mv_a_sad << COST_PRECISION;

    //B pos
    search_center_x = (context_ptr->hme_level0_total_search_area_width*sparce_scale);
    search_center_y = 0;
    ///////////////// correct
    // Correct the left edge of the Search Area if it is not on the reference Picture
    search_center_x = ((origin_x + search_center_x) < -pad_width) ?
        -pad_width - origin_x :
        search_center_x;
    // Correct the right edge of the Search Area if its not on the reference Picture
    search_center_x = ((origin_x + search_center_x) > (int16_t)ref_pic_ptr->width - 1) ?
        search_center_x - ((origin_x + search_center_x) - ((int16_t)ref_pic_ptr->width - 1)) :
        search_center_x;
    // Correct the top edge of the Search Area if it is not on the reference Picture
    search_center_y = ((origin_y + search_center_y) < -pad_height) ?
        -pad_height - origin_y :
        search_center_y;
    // Correct the bottom edge of the Search Area if its not on the reference Picture
    search_center_y = ((origin_y + search_center_y) > (int16_t)ref_pic_ptr->height - 1) ?
        search_center_y - ((origin_y + search_center_y) - ((int16_t)ref_pic_ptr->height - 1)) :
        search_center_y;


    search_region_index = (int16_t)(ref_pic_ptr->origin_x + origin_x) + search_center_x +
        ((int16_t)(ref_pic_ptr->origin_y + origin_y) + search_center_y) * ref_pic_ptr->strideY;

    uint64_t mv_b_sad = NxMSadKernel_funcPtrArray[asm_type][sb_width >> 3](
        context_ptr->sb_src_ptr,
        context_ptr->sb_src_stride << sub_sampled_sad,
        &(ref_pic_ptr->bufferY[search_region_index]),
        ref_pic_ptr->strideY << sub_sampled_sad,
        sb_height >> sub_sampled_sad,
        sb_width);

    mv_b_sad = mv_b_sad << sub_sampled_sad;


    uint64_t mv_b_cost = mv_b_sad << COST_PRECISION;
    //C pos
    search_center_x = 0;
    search_center_y = 0 - (context_ptr->hme_level0_total_search_area_height * sparce_scale);
    ///////////////// correct
    // Correct the left edge of the Search Area if it is not on the reference Picture
    search_center_x = ((origin_x + search_center_x) < -pad_width) ?
        -pad_width - origin_x :
        search_center_x;

    // Correct the right edge of the Search Area if its not on the reference Picture
    search_center_x = ((origin_x + search_center_x) > (int16_t)ref_pic_ptr->width - 1) ?
        search_center_x - ((origin_x + search_center_x) - ((int16_t)ref_pic_ptr->width - 1)) :
        search_center_x;

    // Correct the top edge of the Search Area if it is not on the reference Picture
    search_center_y = ((origin_y + search_center_y) < -pad_height) ?
        -pad_height - origin_y :
        search_center_y;

    // Correct the bottom edge of the Search Area if its not on the reference Picture
    search_center_y = ((origin_y + search_center_y) > (int16_t)ref_pic_ptr->height - 1) ?
        search_center_y - ((origin_y + search_center_y) - ((int16_t)ref_pic_ptr->height - 1)) :
        search_center_y;

    search_region_index = (int16_t)(ref_pic_ptr->origin_x + origin_x) + search_center_x +
        ((int16_t)(ref_pic_ptr->origin_y + origin_y) + search_center_y) * ref_pic_ptr->strideY;

    uint64_t mv_c_sad = NxMSadKernel_funcPtrArray[asm_type][sb_width >> 3](
        context_ptr->sb_src_ptr,
        context_ptr->sb_src_stride << sub_sampled_sad,
        &(ref_pic_ptr->bufferY[search_region_index]),
        ref_pic_ptr->strideY << sub_sampled_sad,
        sb_height >> sub_sampled_sad,
        sb_width);

    mv_c_sad = mv_c_sad << sub_sampled_sad;

    uint64_t mv_c_cost = mv_c_sad << COST_PRECISION;

    //D pos
    search_center_x = 0;
    search_center_y = (context_ptr->hme_level0_total_search_area_height * sparce_scale);
    // Correct the left edge of the Search Area if it is not on the reference Picture
    search_center_x = ((origin_x + search_center_x) < -pad_width) ?
        -pad_width - origin_x :
        search_center_x;
    // Correct the right edge of the Search Area if its not on the reference Picture
    search_center_x = ((origin_x + search_center_x) > (int16_t)ref_pic_ptr->width - 1) ?
        search_center_x - ((origin_x + search_center_x) - ((int16_t)ref_pic_ptr->width - 1)) :
        search_center_x;
    // Correct the top edge of the Search Area if it is not on the reference Picture
    search_center_y = ((origin_y + search_center_y) < -pad_height) ?
        -pad_height - origin_y :
        search_center_y;
    // Correct the bottom edge of the Search Area if its not on the reference Picture
    search_center_y = ((origin_y + search_center_y) > (int16_t)ref_pic_ptr->height - 1) ?
        search_center_y - ((origin_y + search_center_y) - ((int16_t)ref_pic_ptr->height - 1)) :
        search_center_y;
    search_region_index = (int16_t)(ref_pic_ptr->origin_x + origin_x) + search_center_x +
        ((int16_t)(ref_pic_ptr->origin_y + origin_y) + search_center_y) * ref_pic_ptr->strideY;
    uint64_t mv_d_sad = NxMSadKernel_funcPtrArray[asm_type][sb_width >> 3](
        context_ptr->sb_src_ptr,
        context_ptr->sb_src_stride << sub_sampled_sad,
        &(ref_pic_ptr->bufferY[search_region_index]),
        ref_pic_ptr->strideY << sub_sampled_sad,
        sb_height >> sub_sampled_sad,
        sb_width);

    mv_d_sad = mv_d_sad << sub_sampled_sad;

    uint64_t mv_d_cost = mv_d_sad << COST_PRECISION;

    if (list_index == 1) {

        search_center_x = list_index ? 0 - (_MVXT(context_ptr->p_sb_best_mv[0][0][0]) >> 2) : 0;
        search_center_y = list_index ? 0 - (_MVYT(context_ptr->p_sb_best_mv[0][0][0]) >> 2) : 0;
        ///////////////// correct
        // Correct the left edge of the Search Area if it is not on the reference Picture
        search_center_x = ((origin_x + search_center_x) < -pad_width) ?
            -pad_width - origin_x :
            search_center_x;
        // Correct the right edge of the Search Area if its not on the reference Picture
        search_center_x = ((origin_x + search_center_x) > (int16_t)ref_pic_ptr->width - 1) ?
            search_center_x - ((origin_x + search_center_x) - ((int16_t)ref_pic_ptr->width - 1)) :
            search_center_x;
        // Correct the top edge of the Search Area if it is not on the reference Picture
        search_center_y = ((origin_y + search_center_y) < -pad_height) ?
            -pad_height - origin_y :
            search_center_y;
        // Correct the bottom edge of the Search Area if its not on the reference Picture
        search_center_y = ((origin_y + search_center_y) > (int16_t)ref_pic_ptr->height - 1) ?
            search_center_y - ((origin_y + search_center_y) - ((int16_t)ref_pic_ptr->height - 1)) :
            search_center_y;

        search_region_index = (int16_t)(ref_pic_ptr->origin_x + origin_x) + search_center_x +
            ((int16_t)(ref_pic_ptr->origin_y + origin_y) + search_center_y) * ref_pic_ptr->strideY;

        uint64_t direct_mv_sad = NxMSadKernel_funcPtrArray[asm_type][sb_width >> 3](
            context_ptr->sb_src_ptr,
            context_ptr->sb_src_stride << sub_sampled_sad,
            &(ref_pic_ptr->bufferY[search_region_index]),
            ref_pic_ptr->strideY << sub_sampled_sad,
            sb_height >> sub_sampled_sad,
            sb_width);

        direct_mv_sad = direct_mv_sad << sub_sampled_sad;

        direct_mv_cost = (direct_mv_sad << COST_PRECISION);
    }

    best_cost = MIN(zero_mv_cost, MIN(mv_a_cost, MIN(mv_b_cost, MIN(mv_c_cost, MIN(mv_d_cost, direct_mv_cost)))));

    if (best_cost == zero_mv_cost) {
        search_center_x = 0;
        search_center_y = 0;
    }
    else if (best_cost == mv_a_cost) {
        search_center_x = 0 - (context_ptr->hme_level0_total_search_area_width*sparce_scale);
        search_center_y = 0;
    }
    else if (best_cost == mv_b_cost) {
        search_center_x = (context_ptr->hme_level0_total_search_area_width*sparce_scale);
        search_center_y = 0;
    }
    else if (best_cost == mv_c_cost) {
        search_center_x = 0;
        search_center_y = 0 - (context_ptr->hme_level0_total_search_area_height * sparce_scale);
    }
    else if (best_cost == direct_mv_cost) {
        search_center_x = list_index ? 0 - (_MVXT(context_ptr->p_sb_best_mv[0][0][0]) >> 2) : 0;
        search_center_y = list_index ? 0 - (_MVYT(context_ptr->p_sb_best_mv[0][0][0]) >> 2) : 0;
    }
    else if (best_cost == mv_d_cost) {
        search_center_x = 0;
        search_center_y = (context_ptr->hme_level0_total_search_area_height * sparce_scale);
    }

    else {
        SVT_LOG("error no center selected");
    }

    *xsc = search_center_x;
    *ysc = search_center_y;
}
#endif

/*******************************************
* MotionEstimateLcu
*   performs ME (LCU)
*******************************************/
EbErrorType MotionEstimateLcu(
    PictureParentControlSet_t   *picture_control_set_ptr,  // input parameter, Picture Control Set Ptr
    uint32_t                       sb_index,              // input parameter, SB Index
    uint32_t                       sb_origin_x,            // input parameter, SB Origin X
    uint32_t                       sb_origin_y,            // input parameter, SB Origin X
    MeContext_t                    *context_ptr,                        // input parameter, ME Context Ptr, used to store decimated/interpolated LCU/SR
    EbPictureBufferDesc_t       *input_ptr)              // input parameter, source Picture Ptr
{
    EbErrorType return_error = EB_ErrorNone;

    SequenceControlSet_t    *sequence_control_set_ptr = (SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->objectPtr;

    int16_t                  xTopLeftSearchRegion;
    int16_t                  yTopLeftSearchRegion;
    uint32_t                  searchRegionIndex;
    int16_t                  picture_width = (int16_t)((SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->objectPtr)->luma_width;
    int16_t                  picture_height = (int16_t)((SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->objectPtr)->luma_height;
    uint32_t                  sb_width = (input_ptr->width - sb_origin_x) < BLOCK_SIZE_64 ? input_ptr->width - sb_origin_x : BLOCK_SIZE_64;
    uint32_t                  sb_height = (input_ptr->height - sb_origin_y) < BLOCK_SIZE_64 ? input_ptr->height - sb_origin_y : BLOCK_SIZE_64;
    int16_t                  padWidth = (int16_t)BLOCK_SIZE_64 - 1;
    int16_t                  padHeight = (int16_t)BLOCK_SIZE_64 - 1;
    int16_t                  search_area_width;
    int16_t                  search_area_height;
    int16_t                  x_search_area_origin;
    int16_t                  y_search_area_origin;
    int16_t                  origin_x = (int16_t)sb_origin_x;
    int16_t                  origin_y = (int16_t)sb_origin_y;

    // HME
    uint32_t                  searchRegionNumberInWidth = 0;
    uint32_t                  searchRegionNumberInHeight = 0;
    int16_t                  xHmeLevel0SearchCenter[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT][EB_HME_SEARCH_AREA_ROW_MAX_COUNT];
    int16_t                  yHmeLevel0SearchCenter[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT][EB_HME_SEARCH_AREA_ROW_MAX_COUNT];
    uint64_t                  hmeLevel0Sad[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT][EB_HME_SEARCH_AREA_ROW_MAX_COUNT];
    int16_t                  xHmeLevel1SearchCenter[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT][EB_HME_SEARCH_AREA_ROW_MAX_COUNT];
    int16_t                  yHmeLevel1SearchCenter[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT][EB_HME_SEARCH_AREA_ROW_MAX_COUNT];
    uint64_t                  hmeLevel1Sad[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT][EB_HME_SEARCH_AREA_ROW_MAX_COUNT];
    int16_t                  xHmeLevel2SearchCenter[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT][EB_HME_SEARCH_AREA_ROW_MAX_COUNT];
    int16_t                  yHmeLevel2SearchCenter[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT][EB_HME_SEARCH_AREA_ROW_MAX_COUNT];
    uint64_t                  hmeLevel2Sad[EB_HME_SEARCH_AREA_COLUMN_MAX_COUNT][EB_HME_SEARCH_AREA_ROW_MAX_COUNT];

    // Hierarchical ME Search Center
    int16_t                  xHmeSearchCenter = 0;
    int16_t                  yHmeSearchCenter = 0;

    // Final ME Search Center
    int16_t                  xSearchCenter = 0;
    int16_t                  ySearchCenter = 0;

    // Search Center SADs
    uint64_t                  hmeMvSad = 0;

    uint32_t                  pu_index;

    uint32_t                  max_number_of_pus_per_sb = picture_control_set_ptr->max_number_of_pus_per_sb;

    uint32_t                  numOfListToSearch;
    uint32_t                  listIndex;
    uint8_t                   candidateIndex = 0;
    uint8_t                   totalMeCandidateIndex = 0;
    EbPaReferenceObject_t  *referenceObject;  // input parameter, reference Object Ptr

    EbPictureBufferDesc_t  *refPicPtr;
    EbPictureBufferDesc_t  *quarterRefPicPtr;
    EbPictureBufferDesc_t  *sixteenthRefPicPtr;

    int16_t                  tempXHmeSearchCenter = 0;
    int16_t                  tempYHmeSearchCenter = 0;

    uint32_t                  numQuadInWidth;
    uint32_t                  totalMeQuad;
    uint32_t                  quadIndex;
    uint32_t                  nextQuadIndex;
    uint64_t                  tempXHmeSad;

    uint64_t                  ref0Poc = 0;
    uint64_t                  ref1Poc = 0;

    EbAsm                  asm_type = sequence_control_set_ptr->encode_context_ptr->asm_type;
    uint64_t                  i;

    int16_t                  hmeLevel1SearchAreaInWidth;
    int16_t                  hmeLevel1SearchAreaInHeight;

    uint32_t                  adjustSearchAreaDirection = 0;

    // Configure HME level 0, level 1 and level 2 from static config parameters
    EbBool                 enable_hme_level0_flag = picture_control_set_ptr->enable_hme_level0_flag;
    EbBool                 enable_hme_level1_flag = picture_control_set_ptr->enable_hme_level1_flag;
    EbBool                 enable_hme_level2_flag = picture_control_set_ptr->enable_hme_level2_flag;
    EbBool                    enableHalfPel32x32 = EB_FALSE;
    EbBool                    enableHalfPel16x16 = EB_FALSE;
    EbBool                    enableHalfPel8x8 = EB_FALSE;
    EbBool                    enableQuarterPel = EB_FALSE;
#if ENCODER_MODE_CLEANUP
    EbBool                 oneQuadrantHME =  EB_FALSE;
#else
    EbBool                 oneQuadrantHME = (picture_control_set_ptr->enc_mode >= ENC_M3) ? EB_TRUE : EB_FALSE;
#endif

#if M0_SAD_HALF_QUARTER_PEL_BIPRED_SEARCH || M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
    context_ptr->fractionalSearchMethod = (picture_control_set_ptr->enc_mode >= ENC_M3) ? FULL_SAD_SEARCH : SSD_SEARCH;
#else
    context_ptr->fractionalSearchMethod = SUB_SAD_SEARCH;
#endif
#endif
#if M0_64x64_32x32_HALF_QUARTER_PEL
#if ENCODER_MODE_CLEANUP
    context_ptr->fractional_search64x64 = EB_TRUE;
#else
    context_ptr->fractional_search64x64 = (picture_control_set_ptr->enc_mode <= ENC_M1 /*&& sequence_control_set_ptr->static_config.tune != TUNE_VQ*/) ?
        EB_TRUE :
        EB_FALSE;
#endif
#endif
    oneQuadrantHME = sequence_control_set_ptr->input_resolution < INPUT_SIZE_4K_RANGE ? 0 : oneQuadrantHME;
#if M0_ME_SEARCH_BASE
#if ENCODER_MODE_CLEANUP
    numOfListToSearch = (picture_control_set_ptr->slice_type == P_SLICE ) ? (uint32_t)REF_LIST_0 : (uint32_t)REF_LIST_1;

#else
    numOfListToSearch = (picture_control_set_ptr->slice_type == P_SLICE || (picture_control_set_ptr->temporal_layer_index == 0 && picture_control_set_ptr->enc_mode > ENC_M1)) ? (uint32_t)REF_LIST_0 : (uint32_t)REF_LIST_1;
#endif
#else
    numOfListToSearch = (picture_control_set_ptr->slice_type == P_SLICE) || (picture_control_set_ptr->temporal_layer_index == 0) ? (uint32_t)REF_LIST_0 : (uint32_t)REF_LIST_1;
#endif
    referenceObject = (EbPaReferenceObject_t*)picture_control_set_ptr->ref_pa_pic_ptr_array[0]->objectPtr;
    ref0Poc = picture_control_set_ptr->ref_pic_poc_array[0];

    if (numOfListToSearch) {

        referenceObject = (EbPaReferenceObject_t*)picture_control_set_ptr->ref_pa_pic_ptr_array[1]->objectPtr;
        ref1Poc = picture_control_set_ptr->ref_pic_poc_array[1];
    }

    // Uni-Prediction motion estimation loop
    // List Loop
    for (listIndex = REF_LIST_0; listIndex <= numOfListToSearch; ++listIndex) {

        // Ref Picture Loop
        {

            referenceObject = (EbPaReferenceObject_t*)picture_control_set_ptr->ref_pa_pic_ptr_array[listIndex]->objectPtr;
            refPicPtr = (EbPictureBufferDesc_t*)referenceObject->inputPaddedPicturePtr;
            quarterRefPicPtr = (EbPictureBufferDesc_t*)referenceObject->quarterDecimatedPicturePtr;
            sixteenthRefPicPtr = (EbPictureBufferDesc_t*)referenceObject->sixteenthDecimatedPicturePtr;

            if (picture_control_set_ptr->temporal_layer_index > 0 || listIndex == 0) {
                // A - The MV center for Tier0 search could be either (0,0), or HME
                // A - Set HME MV Center
#if HME_ENHANCED_CENTER_SEARCH 
                hme_mv_center_check(
                    refPicPtr,
                    context_ptr,
                    &xSearchCenter,
                    &ySearchCenter,
                    listIndex,
                    origin_x,
                    origin_y,
                    sb_width,
                    sb_height,
                    asm_type);

#else
                xSearchCenter = 0;
                ySearchCenter = 0;
#endif
                // B - NO HME in boundaries
                // C - Skip HME

                if (picture_control_set_ptr->enable_hme_flag && /*B*/sb_height == BLOCK_SIZE_64) {//(searchCenterSad > sequence_control_set_ptr->static_config.skipTier0HmeTh)) {
#if ME_HME_OQ
                    while (searchRegionNumberInHeight < context_ptr->number_hme_search_region_in_height) {
                        while (searchRegionNumberInWidth < context_ptr->number_hme_search_region_in_width) {
#else
                    while (searchRegionNumberInHeight < picture_control_set_ptr->number_hme_search_region_in_height) {
                        while (searchRegionNumberInWidth < picture_control_set_ptr->number_hme_search_region_in_width) {
#endif

                            xHmeLevel0SearchCenter[searchRegionNumberInWidth][searchRegionNumberInHeight] = xSearchCenter;
                            yHmeLevel0SearchCenter[searchRegionNumberInWidth][searchRegionNumberInHeight] = ySearchCenter;

                            xHmeLevel1SearchCenter[searchRegionNumberInWidth][searchRegionNumberInHeight] = xSearchCenter;
                            yHmeLevel1SearchCenter[searchRegionNumberInWidth][searchRegionNumberInHeight] = ySearchCenter;

                            xHmeLevel2SearchCenter[searchRegionNumberInWidth][searchRegionNumberInHeight] = xSearchCenter;
                            yHmeLevel2SearchCenter[searchRegionNumberInWidth][searchRegionNumberInHeight] = ySearchCenter;

                            searchRegionNumberInWidth++;
                        }
                        searchRegionNumberInWidth = 0;
                        searchRegionNumberInHeight++;
                    }

                    // HME: Level0 search


                    if (enable_hme_level0_flag) {

                        if (oneQuadrantHME && !enable_hme_level1_flag && !enable_hme_level2_flag) {

                            searchRegionNumberInHeight = 0;
                            searchRegionNumberInWidth = 0;

                            HmeOneQuadrantLevel0(
                                picture_control_set_ptr,
                                context_ptr,
                                origin_x >> 2,
                                origin_y >> 2,
                                sb_width >> 2,
                                sb_height >> 2,
                                xSearchCenter >> 2,
                                ySearchCenter >> 2,
                                sixteenthRefPicPtr,
                                &(hmeLevel0Sad[searchRegionNumberInWidth][searchRegionNumberInHeight]),
                                &(xHmeLevel0SearchCenter[searchRegionNumberInWidth][searchRegionNumberInHeight]),
                                &(yHmeLevel0SearchCenter[searchRegionNumberInWidth][searchRegionNumberInHeight]),
                                HME_LEVEL_0_SEARCH_AREA_MULTIPLIER_X[picture_control_set_ptr->hierarchical_levels][picture_control_set_ptr->temporal_layer_index],
                                HME_LEVEL_0_SEARCH_AREA_MULTIPLIER_Y[picture_control_set_ptr->hierarchical_levels][picture_control_set_ptr->temporal_layer_index],
                                asm_type);



                        }
                        else
                        {

                            searchRegionNumberInHeight = 0;
                            searchRegionNumberInWidth = 0;
                            {
#if ME_HME_OQ
                                while (searchRegionNumberInHeight < context_ptr->number_hme_search_region_in_height) {
                                    while (searchRegionNumberInWidth < context_ptr->number_hme_search_region_in_width) {
#else
                                while (searchRegionNumberInHeight < picture_control_set_ptr->number_hme_search_region_in_height) {
                                    while (searchRegionNumberInWidth < picture_control_set_ptr->number_hme_search_region_in_width) {

#endif
                                        HmeLevel0(
                                            picture_control_set_ptr,
                                            context_ptr,
                                            origin_x >> 2,
                                            origin_y >> 2,
                                            sb_width >> 2,
                                            sb_height >> 2,
                                            xSearchCenter >> 2,
                                            ySearchCenter >> 2,
                                            sixteenthRefPicPtr,
                                            searchRegionNumberInWidth,
                                            searchRegionNumberInHeight,
                                            &(hmeLevel0Sad[searchRegionNumberInWidth][searchRegionNumberInHeight]),
                                            &(xHmeLevel0SearchCenter[searchRegionNumberInWidth][searchRegionNumberInHeight]),
                                            &(yHmeLevel0SearchCenter[searchRegionNumberInWidth][searchRegionNumberInHeight]),
                                            HME_LEVEL_0_SEARCH_AREA_MULTIPLIER_X[picture_control_set_ptr->hierarchical_levels][picture_control_set_ptr->temporal_layer_index],
                                            HME_LEVEL_0_SEARCH_AREA_MULTIPLIER_Y[picture_control_set_ptr->hierarchical_levels][picture_control_set_ptr->temporal_layer_index],
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
#if ME_HME_OQ
                            while (searchRegionNumberInHeight < context_ptr->number_hme_search_region_in_height) {
                                while (searchRegionNumberInWidth < context_ptr->number_hme_search_region_in_width) {

                                    // When HME level 0 has been disabled, increase the search area width and height for level 1 to (32x12) for Gold only
                                    hmeLevel1SearchAreaInWidth = (int16_t)context_ptr->hme_level1_search_area_in_width_array[searchRegionNumberInWidth];
                                    hmeLevel1SearchAreaInHeight = (int16_t)context_ptr->hme_level1_search_area_in_height_array[searchRegionNumberInHeight];
#else

                            while (searchRegionNumberInHeight < picture_control_set_ptr->number_hme_search_region_in_height) {
                                while (searchRegionNumberInWidth < picture_control_set_ptr->number_hme_search_region_in_width) {

                                    // When HME level 0 has been disabled, increase the search area width and height for level 1 to (32x12) for Gold only
                                    hmeLevel1SearchAreaInWidth = (int16_t)picture_control_set_ptr->hme_level1_search_area_in_width_array[searchRegionNumberInWidth];
                                    hmeLevel1SearchAreaInHeight = (int16_t)picture_control_set_ptr->hme_level1_search_area_in_height_array[searchRegionNumberInHeight];

#endif
                                    HmeLevel1(
                                        context_ptr,
                                        origin_x >> 1,
                                        origin_y >> 1,
                                        sb_width >> 1,
                                        sb_height >> 1,
                                        quarterRefPicPtr,
                                        hmeLevel1SearchAreaInWidth,
                                        hmeLevel1SearchAreaInHeight,
                                        xHmeLevel0SearchCenter[searchRegionNumberInWidth][searchRegionNumberInHeight] >> 1,
                                        yHmeLevel0SearchCenter[searchRegionNumberInWidth][searchRegionNumberInHeight] >> 1,
                                        &(hmeLevel1Sad[searchRegionNumberInWidth][searchRegionNumberInHeight]),
                                        &(xHmeLevel1SearchCenter[searchRegionNumberInWidth][searchRegionNumberInHeight]),
                                        &(yHmeLevel1SearchCenter[searchRegionNumberInWidth][searchRegionNumberInHeight]),
                                        asm_type);


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
#if ME_HME_OQ
                            while (searchRegionNumberInHeight < context_ptr->number_hme_search_region_in_height) {
                                while (searchRegionNumberInWidth < context_ptr->number_hme_search_region_in_width) {
#else
                            while (searchRegionNumberInHeight < picture_control_set_ptr->number_hme_search_region_in_height) {
                                while (searchRegionNumberInWidth < picture_control_set_ptr->number_hme_search_region_in_width) {

#endif

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
                                        xHmeLevel1SearchCenter[searchRegionNumberInWidth][searchRegionNumberInHeight],
                                        yHmeLevel1SearchCenter[searchRegionNumberInWidth][searchRegionNumberInHeight],
                                        &(hmeLevel2Sad[searchRegionNumberInWidth][searchRegionNumberInHeight]),
                                        &(xHmeLevel2SearchCenter[searchRegionNumberInWidth][searchRegionNumberInHeight]),
                                        &(yHmeLevel2SearchCenter[searchRegionNumberInWidth][searchRegionNumberInHeight]),
                                        asm_type);


                                    searchRegionNumberInWidth++;
                                }
                                searchRegionNumberInWidth = 0;
                                searchRegionNumberInHeight++;
                            }
                                }
                            }

                    // Hierarchical ME - Search Center
                    if (enable_hme_level0_flag && !enable_hme_level1_flag && !enable_hme_level2_flag) {

                        if (oneQuadrantHME)
                        {
                            xHmeSearchCenter = xHmeLevel0SearchCenter[0][0];
                            yHmeSearchCenter = yHmeLevel0SearchCenter[0][0];
                            hmeMvSad = hmeLevel0Sad[0][0];
                        }
                        else {

                            xHmeSearchCenter = xHmeLevel0SearchCenter[0][0];
                            yHmeSearchCenter = yHmeLevel0SearchCenter[0][0];
                            hmeMvSad = hmeLevel0Sad[0][0];

                            searchRegionNumberInWidth = 1;
                            searchRegionNumberInHeight = 0;
#if ME_HME_OQ
                            while (searchRegionNumberInHeight < context_ptr->number_hme_search_region_in_height) {
                                while (searchRegionNumberInWidth < context_ptr->number_hme_search_region_in_width) {
#else
                            while (searchRegionNumberInHeight < picture_control_set_ptr->number_hme_search_region_in_height) {
                                while (searchRegionNumberInWidth < picture_control_set_ptr->number_hme_search_region_in_width) {

#endif

                                    xHmeSearchCenter = (hmeLevel0Sad[searchRegionNumberInWidth][searchRegionNumberInHeight] < hmeMvSad) ? xHmeLevel0SearchCenter[searchRegionNumberInWidth][searchRegionNumberInHeight] : xHmeSearchCenter;
                                    yHmeSearchCenter = (hmeLevel0Sad[searchRegionNumberInWidth][searchRegionNumberInHeight] < hmeMvSad) ? yHmeLevel0SearchCenter[searchRegionNumberInWidth][searchRegionNumberInHeight] : yHmeSearchCenter;
                                    hmeMvSad = (hmeLevel0Sad[searchRegionNumberInWidth][searchRegionNumberInHeight] < hmeMvSad) ? hmeLevel0Sad[searchRegionNumberInWidth][searchRegionNumberInHeight] : hmeMvSad;
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
#if ME_HME_OQ
                        while (searchRegionNumberInHeight < context_ptr->number_hme_search_region_in_height) {
                            while (searchRegionNumberInWidth < context_ptr->number_hme_search_region_in_width) {
#else

                        while (searchRegionNumberInHeight < picture_control_set_ptr->number_hme_search_region_in_height) {
                            while (searchRegionNumberInWidth < picture_control_set_ptr->number_hme_search_region_in_width) {
#endif

                                xHmeSearchCenter = (hmeLevel1Sad[searchRegionNumberInWidth][searchRegionNumberInHeight] < hmeMvSad) ? xHmeLevel1SearchCenter[searchRegionNumberInWidth][searchRegionNumberInHeight] : xHmeSearchCenter;
                                yHmeSearchCenter = (hmeLevel1Sad[searchRegionNumberInWidth][searchRegionNumberInHeight] < hmeMvSad) ? yHmeLevel1SearchCenter[searchRegionNumberInWidth][searchRegionNumberInHeight] : yHmeSearchCenter;
                                hmeMvSad = (hmeLevel1Sad[searchRegionNumberInWidth][searchRegionNumberInHeight] < hmeMvSad) ? hmeLevel1Sad[searchRegionNumberInWidth][searchRegionNumberInHeight] : hmeMvSad;
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
#if ME_HME_OQ
                        while (searchRegionNumberInHeight < context_ptr->number_hme_search_region_in_height) {
                            while (searchRegionNumberInWidth < context_ptr->number_hme_search_region_in_width) {
#else

                        while (searchRegionNumberInHeight < picture_control_set_ptr->number_hme_search_region_in_height) {
                            while (searchRegionNumberInWidth < picture_control_set_ptr->number_hme_search_region_in_width) {

#endif
                                xHmeSearchCenter = (hmeLevel2Sad[searchRegionNumberInWidth][searchRegionNumberInHeight] < hmeMvSad) ? xHmeLevel2SearchCenter[searchRegionNumberInWidth][searchRegionNumberInHeight] : xHmeSearchCenter;
                                yHmeSearchCenter = (hmeLevel2Sad[searchRegionNumberInWidth][searchRegionNumberInHeight] < hmeMvSad) ? yHmeLevel2SearchCenter[searchRegionNumberInWidth][searchRegionNumberInHeight] : yHmeSearchCenter;
                                hmeMvSad = (hmeLevel2Sad[searchRegionNumberInWidth][searchRegionNumberInHeight] < hmeMvSad) ? hmeLevel2Sad[searchRegionNumberInWidth][searchRegionNumberInHeight] : hmeMvSad;
                                searchRegionNumberInWidth++;
                            }
                            searchRegionNumberInWidth = 0;
                            searchRegionNumberInHeight++;
                        }

#if ME_HME_OQ
                        numQuadInWidth = context_ptr->number_hme_search_region_in_width;
                        totalMeQuad = context_ptr->number_hme_search_region_in_height * context_ptr->number_hme_search_region_in_width;
#else
                        numQuadInWidth = picture_control_set_ptr->number_hme_search_region_in_width;
                        totalMeQuad = picture_control_set_ptr->number_hme_search_region_in_height * picture_control_set_ptr->number_hme_search_region_in_width;

#endif
                        if ((ref0Poc == ref1Poc) && (listIndex == 1) && (totalMeQuad > 1)) {

                            for (quadIndex = 0; quadIndex < totalMeQuad - 1; ++quadIndex) {
                                for (nextQuadIndex = quadIndex + 1; nextQuadIndex < totalMeQuad; ++nextQuadIndex) {


                                    if (hmeLevel2Sad[quadIndex / numQuadInWidth][quadIndex%numQuadInWidth] > hmeLevel2Sad[nextQuadIndex / numQuadInWidth][nextQuadIndex%numQuadInWidth]) {

                                        tempXHmeSearchCenter = xHmeLevel2SearchCenter[quadIndex / numQuadInWidth][quadIndex%numQuadInWidth];
                                        tempYHmeSearchCenter = yHmeLevel2SearchCenter[quadIndex / numQuadInWidth][quadIndex%numQuadInWidth];
                                        tempXHmeSad = hmeLevel2Sad[quadIndex / numQuadInWidth][quadIndex%numQuadInWidth];

                                        xHmeLevel2SearchCenter[quadIndex / numQuadInWidth][quadIndex%numQuadInWidth] = xHmeLevel2SearchCenter[nextQuadIndex / numQuadInWidth][nextQuadIndex%numQuadInWidth];
                                        yHmeLevel2SearchCenter[quadIndex / numQuadInWidth][quadIndex%numQuadInWidth] = yHmeLevel2SearchCenter[nextQuadIndex / numQuadInWidth][nextQuadIndex%numQuadInWidth];
                                        hmeLevel2Sad[quadIndex / numQuadInWidth][quadIndex%numQuadInWidth] = hmeLevel2Sad[nextQuadIndex / numQuadInWidth][nextQuadIndex%numQuadInWidth];

                                        xHmeLevel2SearchCenter[nextQuadIndex / numQuadInWidth][nextQuadIndex%numQuadInWidth] = tempXHmeSearchCenter;
                                        yHmeLevel2SearchCenter[nextQuadIndex / numQuadInWidth][nextQuadIndex%numQuadInWidth] = tempYHmeSearchCenter;
                                        hmeLevel2Sad[nextQuadIndex / numQuadInWidth][nextQuadIndex%numQuadInWidth] = tempXHmeSad;
                                    }
                                }
                            }

                            xHmeSearchCenter = xHmeLevel2SearchCenter[0][1];
                            yHmeSearchCenter = yHmeLevel2SearchCenter[0][1];
                        }
                            }

                    xSearchCenter = xHmeSearchCenter;
                    ySearchCenter = yHmeSearchCenter;
                        }
                    }

            else {
                xSearchCenter = 0;
                ySearchCenter = 0;
            }
#if ME_HME_OQ
            search_area_width = (int16_t)MIN(context_ptr->search_area_width, 127);
            search_area_height = (int16_t)MIN(context_ptr->search_area_height, 127);
#else
            search_area_width = (int16_t)MIN(picture_control_set_ptr->search_area_width, 127);
            search_area_height = (int16_t)MIN(picture_control_set_ptr->search_area_height, 127);
#endif
#if ENCODER_MODE_CLEANUP
            if ((xSearchCenter != 0 || ySearchCenter != 0) && (picture_control_set_ptr->is_used_as_reference_flag == EB_TRUE)) {

#else
            if ((xSearchCenter != 0 || ySearchCenter != 0) && (picture_control_set_ptr->is_used_as_reference_flag == EB_TRUE || picture_control_set_ptr->enc_mode == ENC_M0)) {
#endif
                CheckZeroZeroCenter(
                    refPicPtr,
                    context_ptr,
                    sb_origin_x,
                    sb_origin_y,
                    sb_width,
                    sb_height,
                    &xSearchCenter,
                    &ySearchCenter,
                    asm_type);
            }
            x_search_area_origin = xSearchCenter - (search_area_width >> 1);
            y_search_area_origin = ySearchCenter - (search_area_height >> 1);

            if (listIndex == 1 && sb_width == BLOCK_SIZE_64 && sb_height == BLOCK_SIZE_64)
            {
                {
                    if (adjustSearchAreaDirection == 1) {
                        search_area_width = (int16_t)8;
                        search_area_height = (int16_t)5;
                    }
                    else if (adjustSearchAreaDirection == 2) {
                        search_area_width = (int16_t)16;
                        search_area_height = (int16_t)9;
                    }

                    x_search_area_origin = xSearchCenter - (search_area_width >> 1);
                    y_search_area_origin = ySearchCenter - (search_area_height >> 1);

                }
            }


            // Correct the left edge of the Search Area if it is not on the reference Picture
            x_search_area_origin = ((origin_x + x_search_area_origin) < -padWidth) ?
                -padWidth - origin_x :
                x_search_area_origin;

            search_area_width = ((origin_x + x_search_area_origin) < -padWidth) ?
                search_area_width - (-padWidth - (origin_x + x_search_area_origin)) :
                search_area_width;

            // Correct the right edge of the Search Area if its not on the reference Picture
            x_search_area_origin = ((origin_x + x_search_area_origin) > picture_width - 1) ?
                x_search_area_origin - ((origin_x + x_search_area_origin) - (picture_width - 1)) :
                x_search_area_origin;

            search_area_width = ((origin_x + x_search_area_origin + search_area_width) > picture_width) ?
                MAX(1, search_area_width - ((origin_x + x_search_area_origin + search_area_width) - picture_width)) :
                search_area_width;

            // Correct the top edge of the Search Area if it is not on the reference Picture
            y_search_area_origin = ((origin_y + y_search_area_origin) < -padHeight) ?
                -padHeight - origin_y :
                y_search_area_origin;

            search_area_height = ((origin_y + y_search_area_origin) < -padHeight) ?
                search_area_height - (-padHeight - (origin_y + y_search_area_origin)) :
                search_area_height;

            // Correct the bottom edge of the Search Area if its not on the reference Picture
            y_search_area_origin = ((origin_y + y_search_area_origin) > picture_height - 1) ?
                y_search_area_origin - ((origin_y + y_search_area_origin) - (picture_height - 1)) :
                y_search_area_origin;

            search_area_height = (origin_y + y_search_area_origin + search_area_height > picture_height) ?
                MAX(1, search_area_height - ((origin_y + y_search_area_origin + search_area_height) - picture_height)) :
                search_area_height;

            context_ptr->x_search_area_origin[listIndex][0] = x_search_area_origin;
            context_ptr->y_search_area_origin[listIndex][0] = y_search_area_origin;

            xTopLeftSearchRegion = (int16_t)(refPicPtr->origin_x + sb_origin_x) - (ME_FILTER_TAP >> 1) + x_search_area_origin;
            yTopLeftSearchRegion = (int16_t)(refPicPtr->origin_y + sb_origin_y) - (ME_FILTER_TAP >> 1) + y_search_area_origin;
            searchRegionIndex = (xTopLeftSearchRegion)+(yTopLeftSearchRegion)* refPicPtr->strideY;

            context_ptr->integer_buffer_ptr[listIndex][0] = &(refPicPtr->bufferY[searchRegionIndex]);
            context_ptr->interpolated_full_stride[listIndex][0] = refPicPtr->strideY;

            // Move to the top left of the search region
            xTopLeftSearchRegion = (int16_t)(refPicPtr->origin_x + sb_origin_x) + x_search_area_origin;
            yTopLeftSearchRegion = (int16_t)(refPicPtr->origin_y + sb_origin_y) + y_search_area_origin;
            searchRegionIndex = xTopLeftSearchRegion + yTopLeftSearchRegion * refPicPtr->strideY;

            {
                {

#if DISABLE_NSQ_FOR_NON_REF || DISABLE_NSQ
#if ENCODER_MODE_CLEANUP
                    if (picture_control_set_ptr->pic_depth_mode <= PIC_ALL_C_DEPTH_MODE) {
#else
                    if (picture_control_set_ptr->non_square_block_flag) {
#endif
#else
                    if (sequence_control_set_ptr->static_config.ext_block_flag) {
#endif
                        uint8_t refPicIndex = 0;

                        InitializeBuffer_32bits_funcPtrArray[asm_type](context_ptr->p_sb_best_sad[listIndex][refPicIndex], 52, 1, MAX_SAD_VALUE);

                        context_ptr->p_best_sad64x64 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_64x64]);
                        context_ptr->p_best_sad32x32 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x32_0]);
                        context_ptr->p_best_sad16x16 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_0]);
                        context_ptr->p_best_sad8x8 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_0]);
                        context_ptr->p_best_sad64x32 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_64x32_0]);
                        context_ptr->p_best_sad32x16 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x16_0]);
                        context_ptr->p_best_sad16x8 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x8_0]);
                        context_ptr->p_best_sad32x64 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x64_0]);
                        context_ptr->p_best_sad16x32 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x32_0]);
                        context_ptr->p_best_sad8x16 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x16_0]);
                        context_ptr->p_best_sad32x8 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x8_0]);
                        context_ptr->p_best_sad8x32 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x32_0]);
                        context_ptr->p_best_sad64x16 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_64x16_0]);
                        context_ptr->p_best_sad16x64 = &(context_ptr->p_sb_best_sad[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x64_0]);

                        context_ptr->p_best_mv64x64 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_64x64]);
                        context_ptr->p_best_mv32x32 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x32_0]);
                        context_ptr->p_best_mv16x16 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_0]);
                        context_ptr->p_best_mv8x8 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_0]);
                        context_ptr->p_best_mv64x32 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_64x32_0]);
                        context_ptr->p_best_mv32x16 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x16_0]);
                        context_ptr->p_best_mv16x8 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x8_0]);
                        context_ptr->p_best_mv32x64 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x64_0]);
                        context_ptr->p_best_mv16x32 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x32_0]);
                        context_ptr->p_best_mv8x16 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x16_0]);
                        context_ptr->p_best_mv32x8 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x8_0]);
                        context_ptr->p_best_mv8x32 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x32_0]);
                        context_ptr->p_best_mv64x16 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_64x16_0]);
                        context_ptr->p_best_mv16x64 = &(context_ptr->p_sb_best_mv[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x64_0]);

#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH

                        context_ptr->p_best_ssd64x64 = &(context_ptr->p_sb_best_ssd[listIndex][refPicIndex][ME_TIER_ZERO_PU_64x64]);
                        context_ptr->p_best_ssd32x32 = &(context_ptr->p_sb_best_ssd[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x32_0]);
                        context_ptr->p_best_ssd16x16 = &(context_ptr->p_sb_best_ssd[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x16_0]);
                        context_ptr->p_best_ssd8x8 = &(context_ptr->p_sb_best_ssd[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x8_0]);
                        context_ptr->p_best_ssd64x32 = &(context_ptr->p_sb_best_ssd[listIndex][refPicIndex][ME_TIER_ZERO_PU_64x32_0]);
                        context_ptr->p_best_ssd32x16 = &(context_ptr->p_sb_best_ssd[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x16_0]);
                        context_ptr->p_best_ssd16x8 = &(context_ptr->p_sb_best_ssd[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x8_0]);
                        context_ptr->p_best_ssd32x64 = &(context_ptr->p_sb_best_ssd[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x64_0]);
                        context_ptr->p_best_ssd16x32 = &(context_ptr->p_sb_best_ssd[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x32_0]);
                        context_ptr->p_best_ssd8x16 = &(context_ptr->p_sb_best_ssd[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x16_0]);
                        context_ptr->p_best_ssd32x8 = &(context_ptr->p_sb_best_ssd[listIndex][refPicIndex][ME_TIER_ZERO_PU_32x8_0]);
                        context_ptr->p_best_ssd8x32 = &(context_ptr->p_sb_best_ssd[listIndex][refPicIndex][ME_TIER_ZERO_PU_8x32_0]);
                        context_ptr->p_best_ssd64x16 = &(context_ptr->p_sb_best_ssd[listIndex][refPicIndex][ME_TIER_ZERO_PU_64x16_0]);
                        context_ptr->p_best_ssd16x64 = &(context_ptr->p_sb_best_ssd[listIndex][refPicIndex][ME_TIER_ZERO_PU_16x64_0]);
#endif


                        open_loop_me_fullpel_search_sblock(
                            context_ptr,
                            listIndex,
                            x_search_area_origin,
                            y_search_area_origin,
                            search_area_width,
                            search_area_height,
                            asm_type);



                    }
                    else {


                        InitializeBuffer_32bits_funcPtrArray[asm_type](context_ptr->p_sb_best_sad[listIndex][0], 21, 1, MAX_SAD_VALUE);
                        context_ptr->p_best_sad64x64 = &(context_ptr->p_sb_best_sad[listIndex][0][ME_TIER_ZERO_PU_64x64]);
                        context_ptr->p_best_sad32x32 = &(context_ptr->p_sb_best_sad[listIndex][0][ME_TIER_ZERO_PU_32x32_0]);
                        context_ptr->p_best_sad16x16 = &(context_ptr->p_sb_best_sad[listIndex][0][ME_TIER_ZERO_PU_16x16_0]);
                        context_ptr->p_best_sad8x8 = &(context_ptr->p_sb_best_sad[listIndex][0][ME_TIER_ZERO_PU_8x8_0]);

                        context_ptr->p_best_mv64x64 = &(context_ptr->p_sb_best_mv[listIndex][0][ME_TIER_ZERO_PU_64x64]);
                        context_ptr->p_best_mv32x32 = &(context_ptr->p_sb_best_mv[listIndex][0][ME_TIER_ZERO_PU_32x32_0]);
                        context_ptr->p_best_mv16x16 = &(context_ptr->p_sb_best_mv[listIndex][0][ME_TIER_ZERO_PU_16x16_0]);
                        context_ptr->p_best_mv8x8 = &(context_ptr->p_sb_best_mv[listIndex][0][ME_TIER_ZERO_PU_8x8_0]);

#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
                        context_ptr->p_best_ssd64x64 = &(context_ptr->p_sb_best_ssd[listIndex][0][ME_TIER_ZERO_PU_64x64]);
                        context_ptr->p_best_ssd32x32 = &(context_ptr->p_sb_best_ssd[listIndex][0][ME_TIER_ZERO_PU_32x32_0]);
                        context_ptr->p_best_ssd16x16 = &(context_ptr->p_sb_best_ssd[listIndex][0][ME_TIER_ZERO_PU_16x16_0]);
                        context_ptr->p_best_ssd8x8 = &(context_ptr->p_sb_best_ssd[listIndex][0][ME_TIER_ZERO_PU_8x8_0]);
#endif
                        FullPelSearch_LCU(
                            context_ptr,
                            listIndex,
                            x_search_area_origin,
                            y_search_area_origin,
                            search_area_width,
                            search_area_height,
                            asm_type
                        );

                    }

                }

                enableHalfPel32x32 = EB_TRUE;
                enableHalfPel16x16 = EB_TRUE;
                enableHalfPel8x8 = EB_TRUE;
#if M0_ME_QUARTER_PEL_SEARCH
                enableQuarterPel = EB_TRUE;
#endif
                if (picture_control_set_ptr->use_subpel_flag == 1) {
#if ENCODER_MODE_CLEANUP
                    if (0) {
#else
                    if (picture_control_set_ptr->enc_mode > ENC_M1) {
#endif
                        suPelEnable(
                            context_ptr,
                            picture_control_set_ptr,
                            listIndex,
                            0,
                            &enableHalfPel32x32,
                            &enableHalfPel16x16,
                            &enableHalfPel8x8);
#if M0_ME_QUARTER_PEL_SEARCH
                        enableQuarterPel = EB_FALSE;
#endif
                    }
#if M0_ME_QUARTER_PEL_SEARCH
                    else {
                        enableQuarterPel = EB_TRUE; // AMIR enable in M1
                    }
#else
                    enableQuarterPel = EB_FALSE;
#endif

                    if (enableHalfPel32x32 || enableHalfPel16x16 || enableHalfPel8x8 || enableQuarterPel) {
                        //if((picture_control_set_ptr->is_used_as_reference_flag == EB_TRUE)) {

                        // Move to the top left of the search region
                        xTopLeftSearchRegion = (int16_t)(refPicPtr->origin_x + sb_origin_x) + x_search_area_origin;
                        yTopLeftSearchRegion = (int16_t)(refPicPtr->origin_y + sb_origin_y) + y_search_area_origin;
                        searchRegionIndex = xTopLeftSearchRegion + yTopLeftSearchRegion * refPicPtr->strideY;

                        // Interpolate the search region for Half-Pel Refinements
                        // H - AVC Style

                        InterpolateSearchRegionAVC(
                            context_ptr,
                            listIndex,
                            context_ptr->integer_buffer_ptr[listIndex][0] + (ME_FILTER_TAP >> 1) + ((ME_FILTER_TAP >> 1) * context_ptr->interpolated_full_stride[listIndex][0]),
                            context_ptr->interpolated_full_stride[listIndex][0],
                            (uint32_t)search_area_width + (BLOCK_SIZE_64 - 1),
                            (uint32_t)search_area_height + (BLOCK_SIZE_64 - 1),
                            8,
                            asm_type);


                        // Half-Pel Refinement [8 search positions]
                        HalfPelSearch_LCU(
                            sequence_control_set_ptr,
#if DISABLE_NSQ_FOR_NON_REF || DISABLE_NSQ
                            picture_control_set_ptr,
#endif
                            context_ptr,
#if M0_HIGH_PRECISION_INTERPOLATION
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
                            context_ptr->integer_buffer_ptr[listIndex][0] + (ME_FILTER_PAD_DISTANCE >> 1) + ((ME_FILTER_PAD_DISTANCE >> 1) * context_ptr->interpolated_full_stride[listIndex][0]),
                            context_ptr->interpolated_full_stride[listIndex][0],
#endif
                            &(context_ptr->pos_b_buffer[listIndex][0][(ME_FILTER_PAD_DISTANCE >> 1) * context_ptr->interpolated_stride]),
#else
#if M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
                            context_ptr->integer_buffer_ptr[listIndex][0] + (ME_FILTER_TAP >> 1) + ((ME_FILTER_TAP >> 1) * context_ptr->interpolated_full_stride[listIndex][0]),
                            context_ptr->interpolated_full_stride[listIndex][0],
#endif
                            &(context_ptr->pos_b_buffer[listIndex][0][(ME_FILTER_TAP >> 1) * context_ptr->interpolated_stride]),
#endif
                            &(context_ptr->pos_h_buffer[listIndex][0][1]),
                            &(context_ptr->pos_j_buffer[listIndex][0][0]),
                            x_search_area_origin,
                            y_search_area_origin,
                            asm_type,
                            picture_control_set_ptr->cu8x8_mode == CU_8x8_MODE_1,
                            enableHalfPel32x32,
                            enableHalfPel16x16,
                            enableHalfPel8x8);

#if M0_ME_QUARTER_PEL_SEARCH
                        // Quarter-Pel Refinement [8 search positions]
                        QuarterPelSearch_LCU(
                            context_ptr,
#if M0_HIGH_PRECISION_INTERPOLATION
                            context_ptr->integer_buffer_ptr[listIndex][0] + (ME_FILTER_PAD_DISTANCE >> 1) + ((ME_FILTER_PAD_DISTANCE >> 1) * context_ptr->interpolated_full_stride[listIndex][0]),
                            context_ptr->interpolated_full_stride[listIndex][0],
                            &(context_ptr->pos_b_buffer[listIndex][0][(ME_FILTER_PAD_DISTANCE >> 1) * context_ptr->interpolated_stride]),  //points to b position of the figure above
#else
                            context_ptr->integer_buffer_ptr[listIndex][0] + (ME_FILTER_TAP >> 1) + ((ME_FILTER_TAP >> 1) * context_ptr->interpolated_full_stride[listIndex][0]),
                            context_ptr->interpolated_full_stride[listIndex][0],
                            &(context_ptr->pos_b_buffer[listIndex][0][(ME_FILTER_TAP >> 1) * context_ptr->interpolated_stride]),  //points to b position of the figure above
#endif
                            &(context_ptr->pos_h_buffer[listIndex][0][1]),                                                      //points to h position of the figure above
                            &(context_ptr->pos_j_buffer[listIndex][0][0]),                                                      //points to j position of the figure above
                            x_search_area_origin,
                            y_search_area_origin,
                            asm_type,
                            picture_control_set_ptr->cu8x8_mode == CU_8x8_MODE_1,
                            enableQuarterPel,
#if DISABLE_NSQ_FOR_NON_REF || DISABLE_NSQ
#if ENCODER_MODE_CLEANUP
                            picture_control_set_ptr->pic_depth_mode <= PIC_ALL_C_DEPTH_MODE);
#else
                            picture_control_set_ptr->non_square_block_flag);
#endif
#else
                            sequence_control_set_ptr->static_config.ext_block_flag);
#endif
#endif

                    }
                }

            }
                        }
                    }

    // Bi-Prediction motion estimation loop
    for (pu_index = 0; pu_index < max_number_of_pus_per_sb; ++pu_index) {

        candidateIndex = 0;

        uint32_t nIdx;

        if (pu_index > 200) {
            nIdx = pu_index;
        }
        else if (pu_index > 184) {
            nIdx = tab8x32[pu_index - 185] + 185;
        }
        else if (pu_index > 168) {
            nIdx = tab32x8[pu_index - 169] + 169;
        }
        else if (pu_index > 136) {
            nIdx = tab8x16[pu_index - 137] + 137;
        }
        else if (pu_index > 128) {
            nIdx = tab16x32[pu_index - 129] + 129;
        }
        else if (pu_index > 126) {
            nIdx = pu_index;
        }
        else if (pu_index > 94) {
            nIdx = tab16x8[pu_index - 95] + 95;
        }
        else if (pu_index > 86) {
            nIdx = tab32x16[pu_index - 87] + 87;
        }
        else if (pu_index > 84) {
            nIdx = pu_index;
        }
        else if (pu_index > 20) {
            nIdx = tab8x8[pu_index - 21] + 21;
        }
        else if (pu_index > 4) {
            nIdx = tab16x16[pu_index - 5] + 5;
        }
        else {
            nIdx = pu_index;
        }


        for (listIndex = REF_LIST_0; listIndex <= numOfListToSearch; ++listIndex) {
            candidateIndex++;
        }


        totalMeCandidateIndex = candidateIndex;

        if (numOfListToSearch) {
#if DISABLE_NSQ_FOR_NON_REF || DISABLE_NSQ
#if ENCODER_MODE_CLEANUP
            if (picture_control_set_ptr->cu8x8_mode == CU_8x8_MODE_0 || pu_index < 21 || (picture_control_set_ptr->pic_depth_mode <= PIC_ALL_C_DEPTH_MODE)) {
#else
            if (picture_control_set_ptr->cu8x8_mode == CU_8x8_MODE_0 || pu_index < 21 || picture_control_set_ptr->non_square_block_flag) {
#endif
#else
            if (picture_control_set_ptr->cu8x8_mode == CU_8x8_MODE_0 || pu_index < 21 || sequence_control_set_ptr->static_config.ext_block_flag) {
#endif
                BiPredictionSearch(
                    context_ptr,
                    pu_index,
                    candidateIndex,
                    1,
                    1,
                    &totalMeCandidateIndex,
                    asm_type,
                    picture_control_set_ptr);
            }
        }

        MeCuResults_t * mePuResult = &picture_control_set_ptr->me_results[sb_index][pu_index];
        mePuResult->totalMeCandidateIndex = totalMeCandidateIndex;

        if (totalMeCandidateIndex == 3) {

            uint32_t L0Sad = context_ptr->p_sb_best_sad[0][0][nIdx];
            uint32_t L1Sad = context_ptr->p_sb_best_sad[1][0][nIdx];
            uint32_t biSad = context_ptr->p_sb_bipred_sad[nIdx];



            mePuResult->xMvL0 = _MVXT(context_ptr->p_sb_best_mv[0][0][nIdx]);
            mePuResult->yMvL0 = _MVYT(context_ptr->p_sb_best_mv[0][0][nIdx]);
            mePuResult->xMvL1 = _MVXT(context_ptr->p_sb_best_mv[1][0][nIdx]);
            mePuResult->yMvL1 = _MVYT(context_ptr->p_sb_best_mv[1][0][nIdx]);

            uint32_t order = Sort3Elements(L0Sad, L1Sad, biSad);

            switch (order) {
                // a = l0Sad, b= l1Sad, c= biSad
            case a_b_c:


                NSET_CAND(mePuResult, 0, L0Sad, UNI_PRED_LIST_0)
                    NSET_CAND(mePuResult, 1, L1Sad, UNI_PRED_LIST_1)
                    NSET_CAND(mePuResult, 2, biSad, BI_PRED)
                    break;

            case a_c_b:

                NSET_CAND(mePuResult, 0, L0Sad, UNI_PRED_LIST_0)
                    NSET_CAND(mePuResult, 1, biSad, BI_PRED)
                    NSET_CAND(mePuResult, 2, L1Sad, UNI_PRED_LIST_1)
                    break;

            case b_a_c:

                NSET_CAND(mePuResult, 0, L1Sad, UNI_PRED_LIST_1)
                    NSET_CAND(mePuResult, 1, L0Sad, UNI_PRED_LIST_0)
                    NSET_CAND(mePuResult, 2, biSad, BI_PRED)
                    break;

            case b_c_a:

                NSET_CAND(mePuResult, 0, L1Sad, UNI_PRED_LIST_1)
                    NSET_CAND(mePuResult, 1, biSad, BI_PRED)
                    NSET_CAND(mePuResult, 2, L0Sad, UNI_PRED_LIST_0)
                    break;

            case c_a_b:

                NSET_CAND(mePuResult, 0, biSad, BI_PRED)
                    NSET_CAND(mePuResult, 1, L0Sad, UNI_PRED_LIST_0)
                    NSET_CAND(mePuResult, 2, L1Sad, UNI_PRED_LIST_1)
                    break;

            case c_b_a:

                NSET_CAND(mePuResult, 0, biSad, BI_PRED)
                    NSET_CAND(mePuResult, 1, L1Sad, UNI_PRED_LIST_1)
                    NSET_CAND(mePuResult, 2, L0Sad, UNI_PRED_LIST_0)
                    break;

            default:
                printf("Err in sorting");
                break;
            }

        }
        else if (totalMeCandidateIndex == 2) {

            uint32_t L0Sad = context_ptr->p_sb_best_sad[0][0][nIdx];
            uint32_t L1Sad = context_ptr->p_sb_best_sad[1][0][nIdx];

            mePuResult->xMvL0 = _MVXT(context_ptr->p_sb_best_mv[0][0][nIdx]);
            mePuResult->yMvL0 = _MVYT(context_ptr->p_sb_best_mv[0][0][nIdx]);
            mePuResult->xMvL1 = _MVXT(context_ptr->p_sb_best_mv[1][0][nIdx]);
            mePuResult->yMvL1 = _MVYT(context_ptr->p_sb_best_mv[1][0][nIdx]);

            if (L0Sad <= L1Sad) {
                NSET_CAND(mePuResult, 0, L0Sad, UNI_PRED_LIST_0)
                    NSET_CAND(mePuResult, 1, L1Sad, UNI_PRED_LIST_1)
            }
            else {
                NSET_CAND(mePuResult, 0, L1Sad, UNI_PRED_LIST_1)
                    NSET_CAND(mePuResult, 1, L0Sad, UNI_PRED_LIST_0)
            }

        }
        else {
            uint32_t L0Sad = context_ptr->p_sb_best_sad[0][0][nIdx];
            mePuResult->xMvL0 = _MVXT(context_ptr->p_sb_best_mv[0][0][nIdx]);
            mePuResult->yMvL0 = _MVYT(context_ptr->p_sb_best_mv[0][0][nIdx]);
            mePuResult->xMvL1 = _MVXT(context_ptr->p_sb_best_mv[1][0][nIdx]);
            mePuResult->yMvL1 = _MVYT(context_ptr->p_sb_best_mv[1][0][nIdx]);
            NSET_CAND(mePuResult, 0, L0Sad, UNI_PRED_LIST_0)
        }


    }



    if (sequence_control_set_ptr->static_config.rate_control_mode) {

        // Compute the sum of the distortion of all 16 16x16 (best) blocks in the LCU
        picture_control_set_ptr->rc_me_distortion[sb_index] = 0;
        for (i = 0; i < 16; i++) {
            picture_control_set_ptr->rc_me_distortion[sb_index] += picture_control_set_ptr->me_results[sb_index][5 + i].distortionDirection[0].distortion;
        }

    }


    return return_error;
                        }


/*******************************************
* SixteenthDecimatedSearch
*  performs a 1/16 decimated search
*******************************************/
uint64_t SixteenthDecimatedSearch(
    MeContext_t             *context_ptr,
    int16_t                   origin_x,
    int16_t                   origin_y,
    uint32_t                   sb_width,
    uint32_t                   sb_height,
    EbPictureBufferDesc_t   *sixteenthRefPicPtr,
    int16_t                     search_area_width,
    int16_t                     search_area_height,
    EbAsm                   asm_type)
{

    int16_t xTopLeftSearchRegion;
    int16_t yTopLeftSearchRegion;
    uint32_t searchRegionIndex;
    int16_t x_search_area_origin;
    int16_t y_search_area_origin;

    int16_t padWidth = (int16_t)(sixteenthRefPicPtr->origin_x) - 1;
    int16_t padHeight = (int16_t)(sixteenthRefPicPtr->origin_y) - 1;

    uint64_t bestSad;
    int16_t xSearchCenter;
    int16_t ySearchCenter;

    x_search_area_origin = -(search_area_width >> 1);
    y_search_area_origin = -(search_area_height >> 1);

    // Correct the left edge of the Search Area if it is not on the reference Picture
    x_search_area_origin = ((origin_x + x_search_area_origin) < -padWidth) ?
        -padWidth - origin_x :
        x_search_area_origin;

    search_area_width = ((origin_x + x_search_area_origin) < -padWidth) ?
        search_area_width - (-padWidth - (origin_x + x_search_area_origin)) :
        search_area_width;

    // Correct the right edge of the Search Area if its not on the reference Picture
    x_search_area_origin = ((origin_x + x_search_area_origin) > (int16_t)sixteenthRefPicPtr->width - 1) ?
        x_search_area_origin - ((origin_x + x_search_area_origin) - ((int16_t)sixteenthRefPicPtr->width - 1)) :
        x_search_area_origin;

    search_area_width = ((origin_x + x_search_area_origin + search_area_width) > (int16_t)sixteenthRefPicPtr->width) ?
        MAX(1, search_area_width - ((origin_x + x_search_area_origin + search_area_width) - (int16_t)sixteenthRefPicPtr->width)) :
        search_area_width;

    // Correct the top edge of the Search Area if it is not on the reference Picture
    y_search_area_origin = ((origin_y + y_search_area_origin) < -padHeight) ?
        -padHeight - origin_y :
        y_search_area_origin;

    search_area_height = ((origin_y + y_search_area_origin) < -padHeight) ?
        search_area_height - (-padHeight - (origin_y + y_search_area_origin)) :
        search_area_height;

    // Correct the bottom edge of the Search Area if its not on the reference Picture
    y_search_area_origin = ((origin_y + y_search_area_origin) > (int16_t)sixteenthRefPicPtr->height - 1) ?
        y_search_area_origin - ((origin_y + y_search_area_origin) - ((int16_t)sixteenthRefPicPtr->height - 1)) :
        y_search_area_origin;

    search_area_height = (origin_y + y_search_area_origin + search_area_height > (int16_t)sixteenthRefPicPtr->height) ?
        MAX(1, search_area_height - ((origin_y + y_search_area_origin + search_area_height) - (int16_t)sixteenthRefPicPtr->height)) :
        search_area_height;

    xTopLeftSearchRegion = ((int16_t)sixteenthRefPicPtr->origin_x + origin_x) + x_search_area_origin;
    yTopLeftSearchRegion = ((int16_t)sixteenthRefPicPtr->origin_y + origin_y) + y_search_area_origin;
    searchRegionIndex = xTopLeftSearchRegion + yTopLeftSearchRegion * sixteenthRefPicPtr->strideY;

    if (((search_area_width & 15) == 0) && (asm_type == ASM_AVX2))
    {
        SadLoopKernel_AVX2_HmeL0_INTRIN(
            &context_ptr->sixteenth_sb_buffer[0],
            context_ptr->sixteenth_sb_buffer_stride * 2,
            &sixteenthRefPicPtr->bufferY[searchRegionIndex],
            sixteenthRefPicPtr->strideY * 2,
            sb_height >> 1, sb_width,
            // results
            &bestSad,
            &xSearchCenter,
            &ySearchCenter,
            // range
            sixteenthRefPicPtr->strideY,
            search_area_width,
            search_area_height
        );
    }
    else if ((search_area_width & 15) == 0)
    {
        // Only width equals 16 (LCU equals 64) is updated
        // other width sizes work with the old code as the one in"SadLoopKernel_SSE4_1_INTRIN"
        SadLoopKernel_SSE4_1_HmeL0_INTRIN(
            &context_ptr->sixteenth_sb_buffer[0],
            context_ptr->sixteenth_sb_buffer_stride * 2,
            &sixteenthRefPicPtr->bufferY[searchRegionIndex],
            sixteenthRefPicPtr->strideY * 2,
            sb_height >> 1, sb_width,
            /* results */
            &bestSad,
            &xSearchCenter,
            &ySearchCenter,
            /* range */
            sixteenthRefPicPtr->strideY,
            search_area_width,
            search_area_height
        );
    }
    else
    {
        // Put the first search location into level0 results
        NxMSadLoopKernel_funcPtrArray[asm_type](
            &context_ptr->sixteenth_sb_buffer[0],
            context_ptr->sixteenth_sb_buffer_stride * 2,
            &sixteenthRefPicPtr->bufferY[searchRegionIndex],
            sixteenthRefPicPtr->strideY * 2,
            sb_height >> 1, sb_width,
            /* results */
            &bestSad,
            &xSearchCenter,
            &ySearchCenter,
            /* range */
            sixteenthRefPicPtr->strideY,
            search_area_width,
            search_area_height
            );
    }

    return(bestSad);
}

/*******************************************
* IsComplexLcu
*   returns true is the SB has a high spatial & temporal complexity
*******************************************/
EbBool IsComplexLcu(
    PictureParentControlSet_t    *previousParentPcs,
    PictureParentControlSet_t    *currentParentPcs,
    PictureParentControlSet_t    *plusOneParentPcs,
    uint32_t pictureWidthInLcus,
    uint32_t lcuAdrr,
    uint32_t sb_origin_x,
    uint32_t sb_origin_y,
    uint32_t sb_width,
    uint32_t sb_height,
    uint32_t lcuCollocatedSad)
{

    uint32_t availableLcusCount = 0;
    uint32_t highVarianceLcusCount = 0;

    // Check the variance of the current LCU
    if ((currentParentPcs->variance[lcuAdrr][ME_TIER_ZERO_PU_64x64]) > IS_COMPLEX_LCU_VARIANCE_TH) {
        availableLcusCount++;
        highVarianceLcusCount++;
    }

    // Check the variance of left SB if available
    if (sb_origin_x != 0) {
        availableLcusCount++;
        if ((currentParentPcs->variance[lcuAdrr - 1][ME_TIER_ZERO_PU_64x64]) > IS_COMPLEX_LCU_VARIANCE_TH) {
            highVarianceLcusCount++;
        }
    }

    // Check the variance of right SB if available
    if ((sb_origin_x + BLOCK_SIZE_64) < currentParentPcs->enhanced_picture_ptr->width) {
        availableLcusCount++;
        if ((currentParentPcs->variance[lcuAdrr + 1][ME_TIER_ZERO_PU_64x64]) > IS_COMPLEX_LCU_VARIANCE_TH) {
            highVarianceLcusCount++;
        }

    }

    // Check the variance of top SB if available
    if (sb_origin_y != 0) {
        availableLcusCount++;
        if ((currentParentPcs->variance[lcuAdrr - pictureWidthInLcus][ME_TIER_ZERO_PU_64x64]) > IS_COMPLEX_LCU_VARIANCE_TH) {
            highVarianceLcusCount++;
        }
    }

    // Check the variance of bottom LCU
    if ((sb_origin_y + BLOCK_SIZE_64) < currentParentPcs->enhanced_picture_ptr->height) {
        availableLcusCount++;
        if ((currentParentPcs->variance[lcuAdrr + pictureWidthInLcus][ME_TIER_ZERO_PU_64x64]) > IS_COMPLEX_LCU_VARIANCE_TH) {
            highVarianceLcusCount++;
        }

    }

    // Check the variance of top-left LCU
    if ((sb_origin_x >= BLOCK_SIZE_64) && (sb_origin_y >= BLOCK_SIZE_64)) {
        availableLcusCount++;
        if ((currentParentPcs->variance[lcuAdrr - pictureWidthInLcus - 1][ME_TIER_ZERO_PU_64x64]) > IS_COMPLEX_LCU_VARIANCE_TH) {
            highVarianceLcusCount++;
        }

    }

    // Check the variance of top-right LCU
    if ((sb_origin_x < currentParentPcs->enhanced_picture_ptr->width - BLOCK_SIZE_64) && (sb_origin_y >= BLOCK_SIZE_64)) {
        availableLcusCount++;
        if ((currentParentPcs->variance[lcuAdrr - pictureWidthInLcus + 1][ME_TIER_ZERO_PU_64x64]) > IS_COMPLEX_LCU_VARIANCE_TH) {
            highVarianceLcusCount++;
        }
    }

    // Check the variance of bottom-left LCU
    if ((sb_origin_x >= BLOCK_SIZE_64) && (sb_origin_y < currentParentPcs->enhanced_picture_ptr->height - BLOCK_SIZE_64)) {
        availableLcusCount++;
        if ((currentParentPcs->variance[lcuAdrr + pictureWidthInLcus - 1][ME_TIER_ZERO_PU_64x64]) > IS_COMPLEX_LCU_VARIANCE_TH) {
            highVarianceLcusCount++;
        }

    }

    // Check the variance of bottom-right LCU
    if ((sb_origin_x < currentParentPcs->enhanced_picture_ptr->width - BLOCK_SIZE_64) && (sb_origin_y < currentParentPcs->enhanced_picture_ptr->height - BLOCK_SIZE_64)) {
        availableLcusCount++;
        if ((currentParentPcs->variance[lcuAdrr + pictureWidthInLcus + 1][ME_TIER_ZERO_PU_64x64]) > IS_COMPLEX_LCU_VARIANCE_TH) {
            highVarianceLcusCount++;
        }
    }

    EbBool varianceFluctuateFlag = EB_FALSE;

    if ((previousParentPcs->variance[lcuAdrr][ME_TIER_ZERO_PU_64x64]) > IS_COMPLEX_LCU_FLAT_VARIANCE_TH &&
        (currentParentPcs->variance[lcuAdrr][ME_TIER_ZERO_PU_64x64]) > IS_COMPLEX_LCU_FLAT_VARIANCE_TH &&
        (plusOneParentPcs->variance[lcuAdrr][ME_TIER_ZERO_PU_64x64]) > IS_COMPLEX_LCU_FLAT_VARIANCE_TH) {

        varianceFluctuateFlag = (EbBool)
            ((((ABS((int32_t)currentParentPcs->variance[lcuAdrr][ME_TIER_ZERO_PU_64x64] - (int32_t)previousParentPcs->variance[lcuAdrr][ME_TIER_ZERO_PU_64x64]) * 100) / (int32_t)previousParentPcs->variance[lcuAdrr][ME_TIER_ZERO_PU_64x64]) >= IS_COMPLEX_LCU_VARIANCE_DEVIATION_TH) &&
            (((ABS((int32_t)currentParentPcs->variance[lcuAdrr][ME_TIER_ZERO_PU_64x64] - (int32_t)plusOneParentPcs->variance[lcuAdrr][ME_TIER_ZERO_PU_64x64]) * 100) / (int32_t)plusOneParentPcs->variance[lcuAdrr][ME_TIER_ZERO_PU_64x64]) >= IS_COMPLEX_LCU_VARIANCE_DEVIATION_TH));

    }

    if (lcuCollocatedSad >= ((sb_width * sb_height) * IS_COMPLEX_LCU_ZZ_SAD_FACTOR_TH) &&
        highVarianceLcusCount >= (availableLcusCount >> 1) &&
        varianceFluctuateFlag) {

        return EB_TRUE;

    }

    return EB_FALSE;

}

/** IntraOpenLoopSearchTheseModesOutputBest

*/
void IntraOpenLoopSearchTheseModesOutputBest(
    uint32_t                               cu_size,
    MotionEstimationContext_t           *context_ptr,
    EbAsm                               asm_type,
    uint8_t                               *src,
    uint32_t                               src_stride,
    uint8_t                                NumOfModesToTest,
    uint32_t                  *stage1ModesArray,
    uint32_t                              *sadArray,
    uint32_t                  *bestMode
)
{
    uint32_t   candidateIndex;
    uint32_t   mode;
    uint32_t   bestSAD = 32 * 32 * 255;

    for (candidateIndex = 0; candidateIndex < NumOfModesToTest; candidateIndex++) {

        mode = stage1ModesArray[candidateIndex];

        // Intra Prediction
        IntraPredictionOpenLoop(
            cu_size,
            context_ptr,
            (uint32_t)mode,
            asm_type);

        //Distortion
        sadArray[candidateIndex] = (uint32_t)NxMSadKernel_funcPtrArray[asm_type][cu_size >> 3]( // Always SAD without weighting
            src,
            src_stride,
            &(context_ptr->me_context_ptr->sb_buffer[0]),
            BLOCK_SIZE_64,
            cu_size,
            cu_size);

        //kepp track of best SAD
        if (sadArray[candidateIndex] < bestSAD) {
            *bestMode = (uint32_t)mode;
            bestSAD = sadArray[candidateIndex];
        }

    }

}

void InjectIntraCandidatesBasedOnBestModeIslice(
    OisCandidate_t              *OisCuPtr,
    uint32_t                        *stage1SadArray,
    uint32_t                       bestMode,
    uint8_t                       *count)
{

    OisCuPtr[(*count)].valid_distortion = EB_TRUE;
    OisCuPtr[(*count)].distortion = stage1SadArray[0];
    OisCuPtr[(*count)++].intra_mode = EB_INTRA_PLANAR;
    OisCuPtr[(*count)++].intra_mode = EB_INTRA_DC;

    switch (bestMode) {

    case EB_INTRA_PLANAR:
    case EB_INTRA_DC:

        break;
    case EB_INTRA_MODE_2:

        OisCuPtr[(*count)++].intra_mode = EB_INTRA_MODE_2;
        OisCuPtr[(*count)++].intra_mode = EB_INTRA_MODE_4;
        OisCuPtr[(*count)++].intra_mode = EB_INTRA_MODE_6;

        break;

    case EB_INTRA_HORIZONTAL:

        OisCuPtr[(*count)++].intra_mode = EB_INTRA_HORIZONTAL;
        OisCuPtr[(*count)++].intra_mode = EB_INTRA_MODE_6;
        OisCuPtr[(*count)++].intra_mode = EB_INTRA_MODE_14;

        break;

    case EB_INTRA_MODE_18:

        OisCuPtr[(*count)++].intra_mode = EB_INTRA_MODE_18;
        OisCuPtr[(*count)++].intra_mode = EB_INTRA_MODE_14;
        OisCuPtr[(*count)++].intra_mode = EB_INTRA_MODE_22;

        break;
    case EB_INTRA_VERTICAL:

        OisCuPtr[(*count)++].intra_mode = EB_INTRA_VERTICAL;
        OisCuPtr[(*count)++].intra_mode = EB_INTRA_MODE_22;
        OisCuPtr[(*count)++].intra_mode = EB_INTRA_MODE_30;

        break;

    case EB_INTRA_MODE_34:
    default:

        OisCuPtr[(*count)++].intra_mode = EB_INTRA_MODE_34;
        OisCuPtr[(*count)++].intra_mode = EB_INTRA_MODE_32;
        OisCuPtr[(*count)++].intra_mode = EB_INTRA_MODE_30;

        break;
    }
}

void InjectIntraCandidatesBasedOnBestMode(
    PictureParentControlSet_t   *picture_control_set_ptr,
    OisCandidate_t              *OisCuPtr,
    uint32_t                        *stage1SadArray,
    uint8_t                        temporal_layer_index,
    uint32_t                       bestMode)
{
    uint32_t count = 0;
    switch (bestMode) {

    case EB_INTRA_MODE_2:
        OisCuPtr[count].distortion = stage1SadArray[2];
#if ENCODER_MODE_CLEANUP
        if (1) 
#else
        if (picture_control_set_ptr->enc_mode <= ENC_M1)
#endif
            OisCuPtr[count].valid_distortion = (temporal_layer_index > 1) ? EB_TRUE : EB_FALSE;
        else
            OisCuPtr[count].valid_distortion = EB_TRUE;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_2;
        OisCuPtr[count++].intra_mode = EB_INTRA_DC;
        OisCuPtr[count++].intra_mode = EB_INTRA_PLANAR;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_3;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_4;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_5;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_7;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_8;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_9;

        break;

    case EB_INTRA_HORIZONTAL:

        OisCuPtr[count].distortion = stage1SadArray[0];

#if ENCODER_MODE_CLEANUP
        if (1)
#else
        if (picture_control_set_ptr->enc_mode <= ENC_M1)
#endif
            OisCuPtr[count].valid_distortion = (temporal_layer_index > 1) ? EB_TRUE : EB_FALSE;
        else
            OisCuPtr[count].valid_distortion = EB_TRUE;
        OisCuPtr[count++].intra_mode = EB_INTRA_HORIZONTAL;
        OisCuPtr[count++].intra_mode = EB_INTRA_DC;
        OisCuPtr[count++].intra_mode = EB_INTRA_PLANAR;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_9;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_11;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_8;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_12;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_7;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_13;


        break;

    case EB_INTRA_MODE_18:

        OisCuPtr[count].distortion = stage1SadArray[3];
#if ENCODER_MODE_CLEANUP
        if (1)
#else
        if (picture_control_set_ptr->enc_mode <= ENC_M1)
#endif
            OisCuPtr[count].valid_distortion = (temporal_layer_index > 1) ? EB_TRUE : EB_FALSE;
        else
            OisCuPtr[count].valid_distortion = EB_TRUE;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_18;
        OisCuPtr[count++].intra_mode = EB_INTRA_DC;
        OisCuPtr[count++].intra_mode = EB_INTRA_PLANAR;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_17;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_19;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_16;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_20;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_15;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_21;

        break;
    case EB_INTRA_VERTICAL:

        OisCuPtr[count].distortion = stage1SadArray[1];

        if (picture_control_set_ptr->enc_mode <= ENC_M1)

            OisCuPtr[count].valid_distortion = (temporal_layer_index > 1) ? EB_TRUE : EB_FALSE;
        else
            OisCuPtr[count].valid_distortion = EB_TRUE;
        OisCuPtr[count++].intra_mode = EB_INTRA_VERTICAL;
        OisCuPtr[count++].intra_mode = EB_INTRA_DC;
        OisCuPtr[count++].intra_mode = EB_INTRA_PLANAR;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_25;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_27;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_24;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_28;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_23;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_29;



        break;

    case EB_INTRA_MODE_34:

        OisCuPtr[count].distortion = stage1SadArray[4];
#if ENCODER_MODE_CLEANUP
        if (1)
#else
        if (picture_control_set_ptr->enc_mode <= ENC_M1)
#endif
            OisCuPtr[count].valid_distortion = (temporal_layer_index > 1) ? EB_TRUE : EB_FALSE;
        else
            OisCuPtr[count].valid_distortion = EB_TRUE;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_34;
        OisCuPtr[count++].intra_mode = EB_INTRA_DC;
        OisCuPtr[count++].intra_mode = EB_INTRA_PLANAR;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_33;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_32;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_29;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_31;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_27;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_28;

        break;

    case EB_INTRA_MODE_6:

        OisCuPtr[count].distortion = stage1SadArray[5];
#if ENCODER_MODE_CLEANUP
        if (1)
#else
        if (picture_control_set_ptr->enc_mode <= ENC_M1)
#endif
            OisCuPtr[count].valid_distortion = (temporal_layer_index > 1) ? EB_TRUE : EB_FALSE;
        else
            OisCuPtr[count].valid_distortion = EB_TRUE;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_6;
        OisCuPtr[count++].intra_mode = EB_INTRA_DC;
        OisCuPtr[count++].intra_mode = EB_INTRA_PLANAR;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_7;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_5;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_4;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_8;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_3;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_9;

        break;

    case EB_INTRA_MODE_14:

        OisCuPtr[count].distortion = stage1SadArray[6];
#if ENCODER_MODE_CLEANUP
        if (1)
#else
        if (picture_control_set_ptr->enc_mode <= ENC_M1)
#endif
            OisCuPtr[count].valid_distortion = (temporal_layer_index > 1) ? EB_TRUE : EB_FALSE;
        else
            OisCuPtr[count].valid_distortion = EB_TRUE;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_14;
        OisCuPtr[count++].intra_mode = EB_INTRA_DC;
        OisCuPtr[count++].intra_mode = EB_INTRA_PLANAR;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_13;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_15;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_12;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_16;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_11;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_17;

        break;

    case EB_INTRA_MODE_22:

        OisCuPtr[count].distortion = stage1SadArray[7];
#if ENCODER_MODE_CLEANUP
        if (1)
#else
        if (picture_control_set_ptr->enc_mode <= ENC_M1)
#endif
            OisCuPtr[count].valid_distortion = (temporal_layer_index > 1) ? EB_TRUE : EB_FALSE;
        else
            OisCuPtr[count].valid_distortion = EB_TRUE;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_22;
        OisCuPtr[count++].intra_mode = EB_INTRA_DC;
        OisCuPtr[count++].intra_mode = EB_INTRA_PLANAR;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_21;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_23;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_20;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_24;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_19;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_25;

        break;

    case EB_INTRA_MODE_30:
    default:

        OisCuPtr[count].distortion = stage1SadArray[8];
#if ENCODER_MODE_CLEANUP
        if (1)
#else
        if (picture_control_set_ptr->enc_mode <= ENC_M1)
#endif
            OisCuPtr[count].valid_distortion = (temporal_layer_index > 1) ? EB_TRUE : EB_FALSE;
        else
            OisCuPtr[count].valid_distortion = EB_TRUE;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_30;
        OisCuPtr[count++].intra_mode = EB_INTRA_DC;
        OisCuPtr[count++].intra_mode = EB_INTRA_PLANAR;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_29;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_31;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_28;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_32;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_27;
        OisCuPtr[count++].intra_mode = EB_INTRA_MODE_33;

        break;


    }
}

int32_t GetInterIntraSadDistance(
    MotionEstimationContext_t   *context_ptr,
    EbPictureBufferDesc_t       *input_ptr,
    uint32_t                       cu_size,
    uint32_t                      *stage1SadArray,
    uint32_t                       meSad,
    uint32_t                       cu_origin_x,
    uint32_t                       cu_origin_y,
    EbAsm                       asm_type
)

{
    int32_t   sadDiff = 0;
    uint8_t   *src = &(input_ptr->bufferY[(input_ptr->origin_y + cu_origin_y) * input_ptr->strideY + (input_ptr->origin_x + cu_origin_x)]);
    // Compute Prediction & SAD for Intra Planer

    // Intra Prediction
    IntraPredictionOpenLoop(
        cu_size,
        context_ptr,
        (uint32_t)1,
        asm_type);

    //Distortion
    stage1SadArray[0] = (uint32_t)NxMSadKernel_funcPtrArray[asm_type][cu_size >> 3]( // Always SAD without weighting
        src,
        input_ptr->strideY,
        &(context_ptr->me_context_ptr->sb_buffer[0]),
        BLOCK_SIZE_64,
        cu_size,
        cu_size);

    // Perform Inter-Intra Comparision
    sadDiff = (int32_t)(meSad - stage1SadArray[0]) * 100;

    return (stage1SadArray[0] ? (sadDiff / (int32_t)stage1SadArray[0]) : 0);

}

void InitValidDistortion(
    OisCandidate_t * oisCuPtr
)


{
    uint8_t intraCandidateIndex;
    for (intraCandidateIndex = 0; intraCandidateIndex < MAX_INTRA_IN_MD; intraCandidateIndex++) {
        oisCuPtr[intraCandidateIndex].valid_distortion = EB_FALSE;
    }
    return;
}

EB_IOS_POINT GetOisPoint(
    uint8_t           oisThSet,
    uint32_t          meSad,
    uint8_t           temporal_layer_index,
    int32_t          interIntraSadDistance,
    uint32_t         *stage1SadArray)

{
    EB_IOS_POINT            oisPoint = OIS_VERY_COMPLEX_MODE;
    // Intra points switcher
    if (stage1SadArray[0] == 0 || meSad == 0 || interIntraSadDistance <= OisPointTh[oisThSet][temporal_layer_index][0]) {
        oisPoint = OIS_VERY_FAST_MODE;
    }
    else
    {
        if (interIntraSadDistance <= OisPointTh[oisThSet][temporal_layer_index][1]) {
            oisPoint = OIS_FAST_MODE;
        }
        else if (interIntraSadDistance <= OisPointTh[oisThSet][temporal_layer_index][2]) {
            oisPoint = OIS_MEDUIM_MODE;
        }
        else if (interIntraSadDistance <= OisPointTh[oisThSet][temporal_layer_index][3]) {
            oisPoint = OIS_COMPLEX_MODE;
        }
    }

    return oisPoint;
}

EbErrorType SortOisCandidateOpenLoop(
    OisCandidate_t                 *oisCandidate)   // input OIS candidate array
{
    EbErrorType                return_error = EB_ErrorNone;
    uint32_t   index1;
    uint32_t   index2;
    uint32_t   intraCandidateMode;
    uint64_t   intraSadDistortion;

    for (index1 = 0; index1 < MAX_OPEN_LOOP_INTRA_CANDIDATES; ++index1)
    {
        for (index2 = index1; index2 < MAX_OPEN_LOOP_INTRA_CANDIDATES; ++index2)
        {
            if (oisCandidate[index1].distortion > oisCandidate[index2].distortion)
            {
                intraCandidateMode = oisCandidate[index1].intra_mode;
                oisCandidate[index1].intra_mode = oisCandidate[index2].intra_mode;
                oisCandidate[index2].intra_mode = intraCandidateMode;

                intraSadDistortion = oisCandidate[index1].distortion;
                oisCandidate[index1].distortion = oisCandidate[index2].distortion;
                oisCandidate[index2].distortion = (uint32_t)intraSadDistortion;
            }
        }
    }
    return return_error;
}

EbErrorType SortIntraModesOpenLoop(
    PictureParentControlSet_t   *picture_control_set_ptr,          // input parameter, pointer to the current lcu
    uint32_t                       sb_index,                      // input parameter, lcu Index
    uint32_t                       cu_index,                       // input parameter, cu index
    uint32_t                       sadDistortion,                 // input parameter, SAD
    uint32_t           openLoopIntraCandidateIndex)   // input parameter, intra mode
{
    EbErrorType                return_error = EB_ErrorNone;

    uint32_t  worstIntraCandidateIndex;
    uint32_t  bestIntraCandidateIndex;
    uint64_t    worstSadDistortion; // could be uint32_t


    OisCu32Cu16Results_t            *oisCu32Cu16ResultsPtr = picture_control_set_ptr->ois_cu32_cu16_results[sb_index];
    OisCu8Results_t                   *oisCu8ResultsPtr = picture_control_set_ptr->ois_cu8_results[sb_index];

    OisCandidate_t * OisCuPtr = cu_index < RASTER_SCAN_CU_INDEX_8x8_0 ?
        oisCu32Cu16ResultsPtr->sorted_ois_candidate[cu_index] : oisCu8ResultsPtr->sorted_ois_candidate[cu_index - RASTER_SCAN_CU_INDEX_8x8_0];


    // Set the best MAX_OPEN_LOOP_INTRA_CANDIDATES modes
    if (openLoopIntraCandidateIndex < MAX_OPEN_LOOP_INTRA_CANDIDATES) {

        OisCuPtr[openLoopIntraCandidateIndex].distortion = sadDistortion;
        OisCuPtr[openLoopIntraCandidateIndex].intra_mode = openLoopIntraCandidateIndex;

        // Get a copy of the OIS SAD and mode - This array will be sorted
        // This array is not used in RC and DeltaQP
        OisCuPtr[openLoopIntraCandidateIndex].distortion = sadDistortion;
        OisCuPtr[openLoopIntraCandidateIndex].intra_mode = openLoopIntraCandidateIndex;

    }
    else {
        uint32_t intraIndex;
        // Determine max SAD distortion
        worstSadDistortion = OisCuPtr[0].distortion;
        worstIntraCandidateIndex = EB_INTRA_PLANAR;

        for (intraIndex = 1; intraIndex < MAX_OPEN_LOOP_INTRA_CANDIDATES; ++intraIndex) {
            bestIntraCandidateIndex = (uint32_t)intraIndex;
            if (OisCuPtr[bestIntraCandidateIndex].distortion > worstSadDistortion) {

                worstSadDistortion = OisCuPtr[bestIntraCandidateIndex].distortion;
                worstIntraCandidateIndex = bestIntraCandidateIndex;
            }
        }

        // Update the best MAX_OPEN_LOOP_INTRA_CANDIDATES modes
        if (sadDistortion < worstSadDistortion) {

            OisCuPtr[worstIntraCandidateIndex].distortion = sadDistortion;
            OisCuPtr[worstIntraCandidateIndex].intra_mode = openLoopIntraCandidateIndex;
        }
    }

    return return_error;
}

uint32_t UpdateNeighborDcIntraPred(
    MotionEstimationContext_t       *context_ptr,
    EbPictureBufferDesc_t           *input_ptr,
    uint32_t                           cu_origin_x,
    uint32_t                           cu_origin_y,
    uint32_t                           cu_size,
    EbAsm                             asm_type)
{
    uint32_t distortion;
    //    printf("cu_size=%i  x=%i  y=%i  rasterScanCuIndex=%i   mdScanCuIndex=%i \n", cu_size, RASTER_SCAN_CU_X[rasterScanCuIndex], RASTER_SCAN_CU_Y[rasterScanCuIndex],rasterScanCuIndex, mdScanCuIndex );
    // Fill Neighbor Arrays
    UpdateNeighborSamplesArrayOpenLoop(
        context_ptr->intra_ref_ptr,
        input_ptr,
        input_ptr->strideY,
        cu_origin_x,
        cu_origin_y,
        cu_size);

    // Intra Prediction
    IntraPredictionOpenLoop(
        cu_size,
        context_ptr,
        (uint32_t)INTRA_DC_MODE,
        asm_type);

    distortion = (uint32_t)NxMSadKernel_funcPtrArray[asm_type][cu_size >> 3]( // Always SAD without weighting
        &(input_ptr->bufferY[(input_ptr->origin_y + cu_origin_y) * input_ptr->strideY + (input_ptr->origin_x + cu_origin_x)]),
        input_ptr->strideY,
        &(context_ptr->me_context_ptr->sb_buffer[0]),
        BLOCK_SIZE_64,
        cu_size,
        cu_size);
    return(distortion);
}

EbErrorType OpenLoopIntraDC(
    PictureParentControlSet_t   *picture_control_set_ptr,
    uint32_t                       sb_index,
    MotionEstimationContext_t   *context_ptr,
    EbPictureBufferDesc_t       *input_ptr,
    EbAsm                       asm_type,
    uint32_t                         cu_origin_x,
    uint32_t                         cu_origin_y,
    uint32_t                         rasterScanCuIndex)
{
    EbErrorType return_error = EB_ErrorNone;

    OisCu32Cu16Results_t            *oisCu32Cu16ResultsPtr = picture_control_set_ptr->ois_cu32_cu16_results[sb_index];
    OisCu8Results_t                   *oisCu8ResultsPtr = picture_control_set_ptr->ois_cu8_results[sb_index];
    OisCandidate_t *OisCuPtr = rasterScanCuIndex < RASTER_SCAN_CU_INDEX_8x8_0 ?
        oisCu32Cu16ResultsPtr->sorted_ois_candidate[rasterScanCuIndex] : oisCu8ResultsPtr->sorted_ois_candidate[rasterScanCuIndex - RASTER_SCAN_CU_INDEX_8x8_0];
    uint32_t cu_size = RASTER_SCAN_CU_SIZE[rasterScanCuIndex];

    if ((cu_size == 32) || (cu_size == 16) || (cu_size == 8))
    {
        if (asm_type == ASM_AVX2)
        {
            OisCuPtr[0].distortion = (uint32_t)UpdateNeighborDcIntraPred_AVX2_INTRIN(
                context_ptr->intra_ref_ptr->y_intra_reference_array_reverse,
                input_ptr->height,
                input_ptr->strideY,
                input_ptr->bufferY,
                input_ptr->origin_y,
                input_ptr->origin_x,
                cu_origin_x,
                cu_origin_y,
                cu_size,
                asm_type);
        }
        else
        {
            OisCuPtr[0].distortion = (uint32_t)UpdateNeighborDcIntraPred(
                context_ptr,
                input_ptr,
                cu_origin_x,
                cu_origin_y,
                cu_size,
                asm_type);
        }

        OisCuPtr[0].intra_mode = INTRA_DC_MODE;
    }
    else {

        OisCuPtr[0].distortion = UpdateNeighborDcIntraPred(
            context_ptr,
            input_ptr,
            cu_origin_x,
            cu_origin_y,
            cu_size,
            asm_type);

        OisCuPtr[0].intra_mode = INTRA_DC_MODE;
    }

    OisCuPtr[0].valid_distortion = 1;

    if (rasterScanCuIndex < RASTER_SCAN_CU_INDEX_8x8_0)
        oisCu32Cu16ResultsPtr->total_intra_luma_mode[rasterScanCuIndex] = 1;
    else
        oisCu8ResultsPtr->total_intra_luma_mode[rasterScanCuIndex - RASTER_SCAN_CU_INDEX_8x8_0] = 1;


    return return_error;
}

uint8_t GetNumOfIntraModesFromOisPoint(
    PictureParentControlSet_t   *picture_control_set_ptr,
    uint32_t                       meSad,
    uint32_t                       oisDcSad
)
{

    int32_t sadDiff = (int32_t)(meSad - oisDcSad) * 100;
    int32_t interIntraSadDistance = oisDcSad ? (sadDiff / (int32_t)oisDcSad) : 0;

    uint8_t oisPoint = GetOisPoint(
        0,
        meSad,
        picture_control_set_ptr->temporal_layer_index,
        interIntraSadDistance,
        &oisDcSad);


    return  numberOfOisModePoints[oisPoint];

}




EbErrorType OpenLoopIntraSearchLcu(
    PictureParentControlSet_t   *picture_control_set_ptr,
    uint32_t                       sb_index,
    MotionEstimationContext_t   *context_ptr,
    EbPictureBufferDesc_t       *input_ptr,
    EbAsm                       asm_type)
{
    EbErrorType return_error = EB_ErrorNone;
    SequenceControlSet_t    *sequence_control_set_ptr = (SequenceControlSet_t*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->objectPtr;
    uint32_t                   rasterScanCuIndex;
    uint32_t                   meSad = 0xFFFFFFFF;
    uint8_t                    stage1NumOfModes = 0;
    int32_t                   interIntraSadDistance = 0;
    uint32_t                     cu_origin_x;
    uint32_t                   cu_origin_y;
    uint32_t                   cu_size;
    uint32_t                   cu_depth;
    uint32_t                     stage1SadArray[11] = { 0 };
    uint32_t                     openLoopIntraCandidateIndex;
    uint32_t                     sadDistortion;
    uint32_t                     intraCandidateIndex;
    uint32_t                   bestMode = EB_INTRA_PLANAR;

    SbParams_t             *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];


    OisCu32Cu16Results_t            *oisCu32Cu16ResultsPtr = picture_control_set_ptr->ois_cu32_cu16_results[sb_index];
    OisCu8Results_t                   *oisCu8ResultsPtr = picture_control_set_ptr->ois_cu8_results[sb_index];


    if (picture_control_set_ptr->slice_type == I_SLICE) {

        for (rasterScanCuIndex = RASTER_SCAN_CU_INDEX_32x32_0; rasterScanCuIndex <= RASTER_SCAN_CU_INDEX_8x8_63; rasterScanCuIndex++) {

            cu_size = RASTER_SCAN_CU_SIZE[rasterScanCuIndex];

            OisCandidate_t *OisCuPtr = rasterScanCuIndex < RASTER_SCAN_CU_INDEX_8x8_0 ?
                oisCu32Cu16ResultsPtr->sorted_ois_candidate[rasterScanCuIndex] : oisCu8ResultsPtr->sorted_ois_candidate[rasterScanCuIndex - RASTER_SCAN_CU_INDEX_8x8_0];

            // Init Valid Distortion to EB_FALSE
            InitValidDistortion(
                OisCuPtr);

            if (sb_params->raster_scan_cu_validity[rasterScanCuIndex]) {

                cu_origin_x = sb_params->origin_x + RASTER_SCAN_CU_X[rasterScanCuIndex];
                cu_origin_y = sb_params->origin_y + RASTER_SCAN_CU_Y[rasterScanCuIndex];
                // Fill Neighbor Arrays
                UpdateNeighborSamplesArrayOpenLoop(
                    context_ptr->intra_ref_ptr,
                    input_ptr,
                    input_ptr->strideY,
                    cu_origin_x,
                    cu_origin_y,
                    cu_size);

                if (cu_size == 32) {

                    // Intra Prediction
                    IntraPredictionOpenLoop(
                        cu_size,
                        context_ptr,
                        (uint32_t)EB_INTRA_PLANAR,
                        asm_type);

                    //Distortion
                    OisCuPtr[0].distortion = (uint32_t)NxMSadKernel_funcPtrArray[asm_type][cu_size >> 3]( // Always SAD without weighting
                        &(input_ptr->bufferY[(input_ptr->origin_y + cu_origin_y) * input_ptr->strideY + (input_ptr->origin_x + cu_origin_x)]),
                        input_ptr->strideY,
                        &(context_ptr->me_context_ptr->sb_buffer[0]),
                        BLOCK_SIZE_64,
                        cu_size,
                        cu_size);


                    OisCuPtr[0].intra_mode = EB_INTRA_PLANAR;
                    OisCuPtr[0].valid_distortion = EB_TRUE;

                }
                else {
                    uint8_t count = 0;
                    IntraOpenLoopSearchTheseModesOutputBest(cu_size,
                        context_ptr,
                        asm_type,
                        &(input_ptr->bufferY[(input_ptr->origin_y + cu_origin_y) * input_ptr->strideY + (input_ptr->origin_x + cu_origin_x)]),
                        input_ptr->strideY,
                        7, // PL ,DC ,  2 , H , 18, V , 34
                        iSliceModesArray,
                        stage1SadArray,
                        &bestMode);

                    InjectIntraCandidatesBasedOnBestModeIslice(
                        OisCuPtr,
                        stage1SadArray,
                        bestMode,
                        &count);


                    if (rasterScanCuIndex < RASTER_SCAN_CU_INDEX_8x8_0)
                        oisCu32Cu16ResultsPtr->total_intra_luma_mode[rasterScanCuIndex] = count;
                    else
                        oisCu8ResultsPtr->total_intra_luma_mode[rasterScanCuIndex - RASTER_SCAN_CU_INDEX_8x8_0] = count;


                }
            }
        }
    }
    else {
        uint8_t oisThSet;

        if (sequence_control_set_ptr->input_resolution == INPUT_SIZE_4K_RANGE) {
#if ENCODER_MODE_CLEANUP
            oisThSet = (  (picture_control_set_ptr->temporal_layer_index == 0 || picture_control_set_ptr->is_used_as_reference_flag == EB_TRUE)) ? 2 : 1;

#else
            oisThSet = (picture_control_set_ptr->enc_mode >= ENC_M2) ?
                0 : // Light
                1; // Default
            oisThSet = ((picture_control_set_ptr->enc_mode <= ENC_M1) && (picture_control_set_ptr->temporal_layer_index == 0 || picture_control_set_ptr->is_used_as_reference_flag == EB_TRUE)) ? 2 : oisThSet;
#endif
        }
        else {
#if ENCODER_MODE_CLEANUP
            oisThSet = 2;

#else
            oisThSet = (picture_control_set_ptr->enc_mode <= ENC_M2) ?
                2 : //Heavy
                (picture_control_set_ptr->enc_mode == ENC_M3) ? 1 :// Default
                0; // Light
#endif
        }



#if ENCODER_MODE_CLEANUP
        EbBool  use16x16Stat = EB_FALSE;
#else
        EbBool  use16x16Stat = (sequence_control_set_ptr->input_resolution == INPUT_SIZE_4K_RANGE
            && picture_control_set_ptr->enc_mode >= ENC_M3);
#endif
        uint32_t   maxCuIndex = use16x16Stat ? RASTER_SCAN_CU_INDEX_16x16_15 : RASTER_SCAN_CU_INDEX_8x8_63;

        for (rasterScanCuIndex = RASTER_SCAN_CU_INDEX_32x32_0; rasterScanCuIndex <= maxCuIndex; rasterScanCuIndex++) {

            EB_IOS_POINT            oisPoint = OIS_VERY_COMPLEX_MODE;

            cu_size = RASTER_SCAN_CU_SIZE[rasterScanCuIndex];
            cu_depth = RASTER_SCAN_CU_DEPTH[rasterScanCuIndex];

            OisCandidate_t * OisCuPtr = rasterScanCuIndex < RASTER_SCAN_CU_INDEX_8x8_0 ?
                oisCu32Cu16ResultsPtr->sorted_ois_candidate[rasterScanCuIndex] : oisCu8ResultsPtr->sorted_ois_candidate[rasterScanCuIndex - RASTER_SCAN_CU_INDEX_8x8_0];


            // Init Valid Distortion to EB_FALSE
            InitValidDistortion(
                OisCuPtr);
            if (sb_params->raster_scan_cu_validity[rasterScanCuIndex] && !((picture_control_set_ptr->cu8x8_mode == CU_8x8_MODE_1) && cu_size == 8)) {
                cu_origin_x = sb_params->origin_x + RASTER_SCAN_CU_X[rasterScanCuIndex];
                cu_origin_y = sb_params->origin_y + RASTER_SCAN_CU_Y[rasterScanCuIndex];
                if (picture_control_set_ptr->limit_ois_to_dc_mode_flag == EB_FALSE) {
                    //    printf("cu_size=%i  x=%i  y=%i  rasterScanCuIndex=%i   mdScanCuIndex=%i \n", cu_size, RASTER_SCAN_CU_X[rasterScanCuIndex], RASTER_SCAN_CU_Y[rasterScanCuIndex],rasterScanCuIndex, mdScanCuIndex );
                    // Fill Neighbor Arrays
                    UpdateNeighborSamplesArrayOpenLoop(
                        context_ptr->intra_ref_ptr,
                        input_ptr,
                        input_ptr->strideY,
                        cu_origin_x,
                        cu_origin_y,
                        cu_size);
                }
#if ENCODER_MODE_CLEANUP
                if ((picture_control_set_ptr->temporal_layer_index == 0) && (sequence_control_set_ptr->input_resolution < INPUT_SIZE_4K_RANGE)) {

#else
                if ((picture_control_set_ptr->temporal_layer_index == 0) && (picture_control_set_ptr->enc_mode <= ENC_M1 && sequence_control_set_ptr->input_resolution < INPUT_SIZE_4K_RANGE)) {
#endif
                    for (intraCandidateIndex = 0; intraCandidateIndex < MAX_OPEN_LOOP_INTRA_CANDIDATES; intraCandidateIndex++) {
                        OisCuPtr[intraCandidateIndex].valid_distortion = EB_FALSE;
                    }

                    uint32_t oisIndex;
                    for (oisIndex = 0; oisIndex < MAX_INTRA_MODES; ++oisIndex) {

                        openLoopIntraCandidateIndex = (uint32_t)oisIndex;
                        // Intra Prediction
                        IntraPredictionOpenLoop(
                            cu_size,
                            context_ptr,
                            openLoopIntraCandidateIndex,
                            asm_type);

                        //Distortion
                        sadDistortion = (uint32_t)NxMSadKernel_funcPtrArray[asm_type][cu_size >> 3](
                            &(input_ptr->bufferY[(input_ptr->origin_y + cu_origin_y) * input_ptr->strideY + (input_ptr->origin_x + cu_origin_x)]),
                            input_ptr->strideY,
                            &(context_ptr->me_context_ptr->sb_buffer[0]),
                            BLOCK_SIZE_64,
                            cu_size,
                            cu_size);

                        // BEST MAX_OPEN_LOOP_INTRA_CANDIDATES
                        SortIntraModesOpenLoop(
                            picture_control_set_ptr,
                            sb_index,
                            rasterScanCuIndex,
                            sadDistortion,
                            openLoopIntraCandidateIndex);


                    }



                    // The sorted array is not used in RC and DeltaQP

                    SortOisCandidateOpenLoop(OisCuPtr);

                    if (rasterScanCuIndex < RASTER_SCAN_CU_INDEX_8x8_0)
                        oisCu32Cu16ResultsPtr->total_intra_luma_mode[rasterScanCuIndex] = MAX_OPEN_LOOP_INTRA_CANDIDATES;
                    else
                        oisCu8ResultsPtr->total_intra_luma_mode[rasterScanCuIndex - RASTER_SCAN_CU_INDEX_8x8_0] = MAX_OPEN_LOOP_INTRA_CANDIDATES;


                }
                else
                {

                    if (picture_control_set_ptr->limit_ois_to_dc_mode_flag == EB_TRUE) {

                        OpenLoopIntraDC(
                            picture_control_set_ptr,
                            sb_index,
                            context_ptr,
                            input_ptr,
                            asm_type,
                            cu_origin_x,
                            cu_origin_y,
                            rasterScanCuIndex);
                    }
                    else {
                        // Set ME distortion

                        meSad = picture_control_set_ptr->me_results[sb_index][rasterScanCuIndex].distortionDirection[0].distortion;

                        interIntraSadDistance = GetInterIntraSadDistance(
                            context_ptr,
                            input_ptr,
                            cu_size,
                            stage1SadArray,
                            meSad,
                            cu_origin_x,
                            cu_origin_y,
                            asm_type);

                        oisPoint = GetOisPoint(
                            oisThSet,
                            meSad,
                            picture_control_set_ptr->temporal_layer_index,
                            interIntraSadDistance,
                            stage1SadArray);



                        if (oisPoint == OIS_VERY_FAST_MODE) {

                            OisCuPtr[0].intra_mode = INTRA_DC_MODE;
                            OisCuPtr[0].distortion = stage1SadArray[0];


                            if (rasterScanCuIndex < RASTER_SCAN_CU_INDEX_8x8_0)
                                oisCu32Cu16ResultsPtr->total_intra_luma_mode[rasterScanCuIndex] = 1;
                            else
                                oisCu8ResultsPtr->total_intra_luma_mode[rasterScanCuIndex - RASTER_SCAN_CU_INDEX_8x8_0] = 1;


                        }
                        else {
                            stage1NumOfModes = numberOfOisModePoints[oisPoint];

                            IntraOpenLoopSearchTheseModesOutputBest(cu_size,
                                context_ptr,
                                asm_type,
                                &(input_ptr->bufferY[(input_ptr->origin_y + cu_origin_y) * input_ptr->strideY + (input_ptr->origin_x + cu_origin_x)]),
                                input_ptr->strideY,
                                stage1NumOfModes,
                                stage1ModesArray,
                                stage1SadArray,
                                &bestMode);


                            InjectIntraCandidatesBasedOnBestMode(
                                picture_control_set_ptr,
                                OisCuPtr,
                                stage1SadArray,
                                picture_control_set_ptr->temporal_layer_index,
                                bestMode);



                            if (rasterScanCuIndex < RASTER_SCAN_CU_INDEX_8x8_0)
                                oisCu32Cu16ResultsPtr->total_intra_luma_mode[rasterScanCuIndex] = intraSearchInMd[oisPoint][cu_depth];
                            else
                                oisCu8ResultsPtr->total_intra_luma_mode[rasterScanCuIndex - RASTER_SCAN_CU_INDEX_8x8_0] = intraSearchInMd[oisPoint][cu_depth];




                        }
                    }
                }

            }
        }
    }

    return return_error;
}





