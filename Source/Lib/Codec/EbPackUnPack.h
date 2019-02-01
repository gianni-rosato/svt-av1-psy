/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPackUnPack_h
#define EbPackUnPack_h
#ifdef __cplusplus
extern "C" {
#endif

#include "EbPackUnPack_C.h"
#include "EbPackUnPack_SSE2.h"
#include "EbPictureOperators.h"

    typedef void(*EB_ENC_Pack2D_TYPE)(
        uint8_t     *in8BitBuffer,
        uint32_t     in8Stride,
        uint8_t     *innBitBuffer,
        uint16_t    *out16BitBuffer,
        uint32_t     innStride,
        uint32_t     outStride,
        uint32_t     width,
        uint32_t     height);

    EB_ENC_Pack2D_TYPE Pack2D_funcPtrArray_16Bit_SRC[2][ASM_TYPE_TOTAL] =
    {
        {
            // NON_AVX2
            EB_ENC_msbPack2D,
            // AVX2
            EB_ENC_msbPack2D,
        },
        {
            // NON_AVX2
            EB_ENC_msbPack2D_SSE2_INTRIN,
            // AVX2
            EB_ENC_msbPack2D_AVX2_INTRIN_AL,//EB_ENC_msbPack2D_AVX2
        }
    };


    EB_ENC_Pack2D_TYPE CompressedPack_funcPtrArray[ASM_TYPE_TOTAL] =
    {
        // NON_AVX2
        CompressedPackmsb,
        // AVX2
        CompressedPackmsb_AVX2_INTRIN,

    };

    typedef void(*COMPPack_TYPE)(
        const uint8_t     *innBitBuffer,
        uint32_t     innStride,
        uint8_t     *inCompnBitBuffer,
        uint32_t     outStride,
        uint8_t    *localCache,
        uint32_t     width,
        uint32_t     height);

    COMPPack_TYPE  Convert_Unpack_CPack_funcPtrArray[ASM_TYPE_TOTAL] =
    {
        // NON_AVX2
        CPack_C,
        // AVX2
        CPack_AVX2_INTRIN,

    };

    typedef void(*EB_ENC_UnPack2D_TYPE)(
        uint16_t      *in16BitBuffer,
        uint32_t       inStride,
        uint8_t       *out8BitBuffer,
        uint8_t       *outnBitBuffer,
        uint32_t       out8Stride,
        uint32_t       outnStride,
        uint32_t       width,
        uint32_t       height);

    EB_ENC_UnPack2D_TYPE UnPack2D_funcPtrArray_16Bit[2][ASM_TYPE_TOTAL] =
    {
        {
            // NON_AVX2
            EB_ENC_msbUnPack2D,
            // AVX2
            EB_ENC_msbUnPack2D,
        },
        {
            // NON_AVX2
            EB_ENC_msbUnPack2D_SSE2_INTRIN,
            // AVX2
            EB_ENC_msbUnPack2D_SSE2_INTRIN,
        }
    };

    typedef void(*EB_ENC_UnpackAvg_TYPE)(
        uint16_t *ref16L0,
        uint32_t  refL0Stride,
        uint16_t *ref16L1,
        uint32_t  refL1Stride,
        uint8_t  *dstPtr,
        uint32_t  dst_stride,
        uint32_t  width,
        uint32_t  height);
    EB_ENC_UnpackAvg_TYPE UnPackAvg_funcPtrArray[ASM_TYPE_TOTAL] =
    {
        // NON_AVX2
        UnpackAvg,
        // AVX2
        UnpackAvg_AVX2_INTRIN,//UnpackAvg_SSE2_INTRIN,

    };
    typedef void(*EB_ENC_UnpackAvgSub_TYPE)(
        uint16_t *ref16L0,
        uint32_t  refL0Stride,
        uint16_t *ref16L1,
        uint32_t  refL1Stride,
        uint8_t  *dstPtr,
        uint32_t  dst_stride,
        EbBool      subPred,
        uint32_t  width,
        uint32_t  height);
    EB_ENC_UnpackAvgSub_TYPE UnPackAvgSafeSub_funcPtrArray[ASM_TYPE_TOTAL] =
    {
        // NON_AVX2
        UnpackAvgSafeSub,
        // AVX2  SafeSub
        UnpackAvgSafeSub_AVX2_INTRIN,//UnpackAvg_SSE2_INTRIN,

    };

    typedef void(*EB_ENC_UnPack8BitData_TYPE)(
        uint16_t      *in16BitBuffer,
        uint32_t       inStride,
        uint8_t       *out8BitBuffer,
        uint32_t       out8Stride,
        uint32_t       width,
        uint32_t       height);
    EB_ENC_UnPack8BitData_TYPE UnPack8BIT_funcPtrArray_16Bit[2][ASM_TYPE_TOTAL] =
    {
        {
           UnPack8BitData,
           UnPack8BitData,
        },
        {
            // NON_AVX2
            EB_ENC_UnPack8BitData_SSE2_INTRIN,
            // AVX2
            EB_ENC_UnPack8BitData_SSE2_INTRIN,
        }
    };
    typedef void(*EB_ENC_UnPack8BitDataSUB_TYPE)(
        uint16_t      *in16BitBuffer,
        uint32_t       inStride,
        uint8_t       *out8BitBuffer,
        uint32_t       out8Stride,
        uint32_t       width,
        uint32_t       height,
        EbBool      subPred
        );
    EB_ENC_UnPack8BitDataSUB_TYPE UnPack8BITSafeSub_funcPtrArray_16Bit[ASM_TYPE_TOTAL] =
    {
        // NON_AVX2
        EB_ENC_UnPack8BitDataSafeSub_SSE2_INTRIN,
        // AVX2
        EB_ENC_UnPack8BitDataSafeSub_SSE2_INTRIN,

    };

#ifdef __cplusplus
}
#endif
#endif // EbPackUnPack_h