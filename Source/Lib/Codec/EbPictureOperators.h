/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPictureOperators_h
#define EbPictureOperators_h

#include "EbPictureOperators_C.h"
#include "EbPictureOperators_SSE2.h"
#include "EbPictureOperators_SSE4_1.h"
#include "EbPictureOperators_AVX2.h"
#include "EbHmCode.h"
#include "EbDefinitions.h"
#include "EbPictureBufferDesc.h"
#include "EbSequenceControlSet.h"
#ifdef __cplusplus
extern "C" {
#endif

    extern void PictureAddition(
        uint8_t  *predPtr,
        uint32_t  predStride,
        int16_t *residual_ptr,
        uint32_t  residualStride,
        uint8_t  *reconPtr,
        uint32_t  reconStride,
        uint32_t  width,
        uint32_t  height,
        EbAsm  asm_type);

    extern EbErrorType PictureCopy8Bit(
        EbPictureBufferDesc_t   *src,
        uint32_t                   srcLumaOriginIndex,
        uint32_t                   srcChromaOriginIndex,
        EbPictureBufferDesc_t   *dst,
        uint32_t                   dstLumaOriginIndex,
        uint32_t                   dstChromaOriginIndex,
        uint32_t                   areaWidth,
        uint32_t                   areaHeight,
        uint32_t                   chromaAreaWidth,
        uint32_t                   chromaAreaHeight,
        uint32_t                   component_mask,
        EbAsm                   asm_type);
    extern EbErrorType PictureFullDistortion32Bits(
        EbPictureBufferDesc_t   *coeff,
        uint32_t                   coeffLumaOriginIndex,
        uint32_t                   coeffChromaOriginIndex,
        EbPictureBufferDesc_t   *reconCoeff,
        uint32_t                   reconCoeffLumaOriginIndex,
        uint32_t                   reconCoeffChromaOriginIndex,
        uint32_t                   bwidth,
        uint32_t                   bheight,
        uint32_t                   bwidth_uv,
        uint32_t                   bheight_uv,
        uint64_t                   y_distortion[DIST_CALC_TOTAL],
        uint64_t                   cb_distortion[DIST_CALC_TOTAL],
        uint64_t                   cr_distortion[DIST_CALC_TOTAL],
        uint32_t                   y_count_non_zero_coeffs,
        uint32_t                   cbCountNonZeroCoeffs,
        uint32_t                   crCountNonZeroCoeffs,
        COMPONENT_TYPE            componentType,
        EbAsm                   asm_type);

    extern uint64_t ComputeNxMSatdSadLCU(
        uint8_t  *src,        // input parameter, source samples Ptr
        uint32_t  src_stride,  // input parameter, source stride
        uint32_t  width,      // input parameter, block width (N)
        uint32_t  height,     // input parameter, block height (M)
        EbAsm  asm_type);

    //Residual Data

    void CompressedPackLcu(
        uint8_t     *in8BitBuffer,
        uint32_t     in8Stride,
        uint8_t     *innBitBuffer,
        uint32_t     innStride,
        uint16_t    *out16BitBuffer,
        uint32_t     outStride,
        uint32_t     width,
        uint32_t     height,
        EbAsm     asm_type
    );
    void Conv2bToCPackLcu(
        const uint8_t     *innBitBuffer,
        uint32_t     innStride,
        uint8_t     *inCompnBitBuffer,
        uint32_t     outStride,
        uint8_t    *localCache,
        uint32_t     width,
        uint32_t     height,
        EbAsm     asm_type);


    void Pack2D_SRC(
        uint8_t     *in8BitBuffer,
        uint32_t     in8Stride,
        uint8_t     *innBitBuffer,
        uint32_t     innStride,
        uint16_t    *out16BitBuffer,
        uint32_t     outStride,
        uint32_t     width,
        uint32_t     height,
        EbAsm      asm_type);

    void UnPack2D(
        uint16_t      *in16BitBuffer,
        uint32_t       inStride,
        uint8_t       *out8BitBuffer,
        uint32_t       out8Stride,
        uint8_t       *outnBitBuffer,
        uint32_t       outnStride,
        uint32_t       width,
        uint32_t       height,
        EbAsm          asm_type);

    void extract_8bit_data(
        uint16_t      *in16BitBuffer,
        uint32_t       inStride,
        uint8_t       *out8BitBuffer,
        uint32_t       out8Stride,
        uint32_t       width,
        uint32_t       height,
        EbAsm       asm_type
    );
    void unpack_l0l1_avg(
        uint16_t *ref16L0,
        uint32_t  refL0Stride,
        uint16_t *ref16L1,
        uint32_t  refL1Stride,
        uint8_t  *dstPtr,
        uint32_t  dst_stride,
        uint32_t  width,
        uint32_t  height,
        EbAsm  asm_type);
    void extract8_bitdata_safe_sub(
        uint16_t      *in16BitBuffer,
        uint32_t       inStride,
        uint8_t       *out8BitBuffer,
        uint32_t       out8Stride,
        uint32_t       width,
        uint32_t       height,
        EbBool      subPred,
        EbAsm       asm_type
    );
    void unpack_l0l1_avg_safe_sub(
        uint16_t *ref16L0,
        uint32_t  refL0Stride,
        uint16_t *ref16L1,
        uint32_t  refL1Stride,
        uint8_t  *dstPtr,
        uint32_t  dst_stride,
        uint32_t  width,
        uint32_t  height,
        EbBool      subPred,
        EbAsm  asm_type);

    void memcpy16bit(
        uint16_t                     * outPtr,
        uint16_t                     * inPtr,
        uint64_t                       numOfElements);

    void memcpy32bit(
        uint32_t                     * outPtr,
        uint32_t                     * inPtr,
        uint64_t                       numOfElements);

    void memset16bit(
        uint16_t                     * inPtr,
        uint16_t                       value,
        uint64_t                       numOfElements);


    void memset32bit(
        uint32_t                     * inPtr,
        uint32_t                       value,
        uint64_t                       numOfElements);

    static void PictureAdditionVoidFunc() {}
    static void PicResdVoidFunc() {}
    static void PicZeroOutCoefVoidFunc() {}

    int32_t  sumResidual(int16_t * inPtr,
        uint32_t   size,
        uint32_t   strideIn);


    typedef int32_t(*EB_SUM_RES)(
        int16_t * inPtr,
        uint32_t   size,
        uint32_t   strideIn);

    static EB_SUM_RES FUNC_TABLE SumResidual_funcPtrArray[ASM_TYPE_TOTAL] = {
        // NON_AVX2
        sumResidual,
        // AVX2
        sumResidual8bit_AVX2_INTRIN,
    };

    void memset16bitBlock(
        int16_t * inPtr,
        uint32_t   strideIn,
        uint32_t   size,
        int16_t   value
    );

    typedef void(*EB_MEMSET16bitBLK)(
        int16_t * inPtr,
        uint32_t   strideIn,
        uint32_t   size,
        int16_t   value);

    static EB_MEMSET16bitBLK FUNC_TABLE memset16bitBlock_funcPtrArray[ASM_TYPE_TOTAL] = {
        // NON_AVX2
        memset16bitBlock,
        // AVX2
        memset16bitBlock_AVX2_INTRIN,
    };

    void FullDistortionKernelCbfZero32Bits(
        int32_t  *coeff,
        uint32_t   coeffStride,
        int32_t  *reconCoeff,
        uint32_t   reconCoeffStride,
        uint64_t   distortionResult[DIST_CALC_TOTAL],
        uint32_t   areaWidth,
        uint32_t   areaHeight);

    void FullDistortionKernel32Bits(
        int32_t  *coeff,
        uint32_t   coeffStride,
        int32_t  *reconCoeff,
        uint32_t   reconCoeffStride,
        uint64_t   distortionResult[DIST_CALC_TOTAL],
        uint32_t   areaWidth,
        uint32_t   areaHeight);

    typedef void(*EB_FUllDISTORTIONKERNELCBFZERO32BITS)(
        int32_t  *coeff,
        uint32_t   coeffStride,
        int32_t  *reconCoeff,
        uint32_t   reconCoeffStride,
        uint64_t   distortionResult[DIST_CALC_TOTAL],
        uint32_t   areaWidth,
        uint32_t   areaHeight);

    typedef void(*EB_FUllDISTORTIONKERNEL32BITS)(
        int32_t  *coeff,
        uint32_t   coeffStride,
        int32_t  *reconCoeff,
        uint32_t   reconCoeffStride,
        uint64_t   distortionResult[DIST_CALC_TOTAL],
        uint32_t   areaWidth,
        uint32_t   areaHeight);

    static EB_FUllDISTORTIONKERNELCBFZERO32BITS FUNC_TABLE FullDistortionKernelCbfZero32Bits_funcPtrArray[ASM_TYPE_TOTAL] = {
        // NON_AVX2
        FullDistortionKernelCbfZero32Bits,
        // AVX2
        FullDistortionKernelCbfZero32Bits_AVX2,
    };

    static EB_FUllDISTORTIONKERNEL32BITS FUNC_TABLE FullDistortionKernel32Bits_funcPtrArray[ASM_TYPE_TOTAL] = {
        // NON_AVX2
        FullDistortionKernel32Bits,
        // AVX2
        FullDistortionKernel32Bits_AVX2,
    };

    /***************************************
    * Function Types
    ***************************************/
    typedef void(*EB_ADDDKERNEL_TYPE)(
        uint8_t  *predPtr,
        uint32_t  predStride,
        int16_t *residual_ptr,
        uint32_t  residualStride,
        uint8_t  *reconPtr,
        uint32_t  reconStride,
        uint32_t  width,
        uint32_t  height);

    typedef void(*EB_ADDDKERNEL_AV1_TYPE)(
        uint8_t  *predPtr,
        uint32_t  predStride,
        int32_t *residual_ptr,
        uint32_t  residualStride,
        uint8_t  *reconPtr,
        uint32_t  reconStride,
        uint32_t  width,
        uint32_t  height,
        int32_t     bd);

    typedef void(*EB_ADDDKERNEL_TYPE_16BIT)(
        uint16_t  *predPtr,
        uint32_t  predStride,
        int16_t *residual_ptr,
        uint32_t  residualStride,
        uint16_t  *reconPtr,
        uint32_t  reconStride,
        uint32_t  width,
        uint32_t  height);

    typedef void(*EB_PICCOPY_TYPE)(
        EbByte                  src,
        uint32_t                   src_stride,
        EbByte                  dst,
        uint32_t                   dst_stride,
        uint32_t                   areaWidth,
        uint32_t                   areaHeight);



    typedef void(*EB_RESDKERNEL_TYPE)(
        uint8_t   *input,
        uint32_t   inputStride,
        uint8_t   *pred,
        uint32_t   predStride,
        int16_t  *residual,
        uint32_t   residualStride,
        uint32_t   areaWidth,
        uint32_t   areaHeight);

    typedef void(*EB_RESDKERNEL_TYPE_16BIT)(
        uint16_t   *input,
        uint32_t   inputStride,
        uint16_t   *pred,
        uint32_t   predStride,
        int16_t  *residual,
        uint32_t   residualStride,
        uint32_t   areaWidth,
        uint32_t   areaHeight);

    typedef void(*EB_ZEROCOEFF_TYPE)(
        int16_t*                  coeffbuffer,
        uint32_t                   coeffStride,
        uint32_t                   coeffOriginIndex,
        uint32_t                   areaWidth,
        uint32_t                   areaHeight);

    typedef void(*EB_FULLDIST_TYPE)(
        int16_t  *coeff,
        uint32_t   coeffStride,
        int16_t  *reconCoeff,
        uint32_t   reconCoeffStride,
        uint64_t   distortionResult[DIST_CALC_TOTAL],
        uint32_t   areaWidth,
        uint32_t   areaHeight);

    typedef uint64_t(*EB_SATD_TYPE)(
        int16_t *diff);

    typedef uint64_t(*EB_SATD_U8_TYPE)(
        uint8_t *diff,
        uint64_t *dcValue,
        uint32_t  src_stride);

    /***************************************
    * Function Tables
    ***************************************/
    static EB_ADDDKERNEL_TYPE FUNC_TABLE AdditionKernel_funcPtrArray[ASM_TYPE_TOTAL][9] = {
        // NON_AVX2
        {
            /*0 4x4   */    PictureAdditionKernel4x4_SSE_INTRIN,
            /*1 8x8   */    PictureAdditionKernel8x8_SSE2_INTRIN,
            /*2 16x16 */    PictureAdditionKernel16x16_SSE2_INTRIN,
            /*3       */    (EB_ADDDKERNEL_TYPE)PictureAdditionVoidFunc,
            /*4 32x32 */    PictureAdditionKernel32x32_SSE2_INTRIN,
            /*5       */    (EB_ADDDKERNEL_TYPE)PictureAdditionVoidFunc,
            /*6       */    (EB_ADDDKERNEL_TYPE)PictureAdditionVoidFunc,
            /*7       */    (EB_ADDDKERNEL_TYPE)PictureAdditionVoidFunc,
            /*8 64x64 */    PictureAdditionKernel64x64_SSE2_INTRIN,
        },
        // AVX2
        {
            /*0 4x4   */    PictureAdditionKernel4x4_SSE_INTRIN,
            /*1 8x8   */    PictureAdditionKernel8x8_SSE2_INTRIN,
            /*2 16x16 */    PictureAdditionKernel16x16_SSE2_INTRIN,
            /*3       */    (EB_ADDDKERNEL_TYPE)PictureAdditionVoidFunc,
            /*4 32x32 */    PictureAdditionKernel32x32_SSE2_INTRIN,
            /*5       */    (EB_ADDDKERNEL_TYPE)PictureAdditionVoidFunc,
            /*6       */    (EB_ADDDKERNEL_TYPE)PictureAdditionVoidFunc,
            /*7       */    (EB_ADDDKERNEL_TYPE)PictureAdditionVoidFunc,
            /*8 64x64 */    PictureAdditionKernel64x64_SSE2_INTRIN,
        },
    };


    static EB_ADDDKERNEL_AV1_TYPE FUNC_TABLE AdditionKernel_Av1_funcPtrArray[ASM_TYPE_TOTAL][9] = {
        // NON_AVX2
        {
            /*0 4x4   */    PictureAdditionKernel,
            /*1 8x8   */    PictureAdditionKernel,
            /*2 16x16 */    PictureAdditionKernel,
            /*3       */    (EB_ADDDKERNEL_AV1_TYPE)PictureAdditionVoidFunc,
            /*4 32x32 */    PictureAdditionKernel,
            /*5       */    (EB_ADDDKERNEL_AV1_TYPE)PictureAdditionVoidFunc,
            /*6       */    (EB_ADDDKERNEL_AV1_TYPE)PictureAdditionVoidFunc,
            /*7       */    (EB_ADDDKERNEL_AV1_TYPE)PictureAdditionVoidFunc,
            /*8 64x64 */    PictureAdditionKernel,
        },
        // AVX2
        {
            /*0 4x4   */    PictureAdditionKernel4x4_AV1_SSE2_INTRIN,
            /*1 8x8   */    PictureAdditionKernel8x8_AV1_SSE2_INTRIN,
            /*2 16x16 */    PictureAdditionKernel16x16_AV1_SSE2_INTRIN,
            /*3       */    (EB_ADDDKERNEL_AV1_TYPE)PictureAdditionVoidFunc,
            /*4 32x32 */    PictureAdditionKernel32x32_AV1_SSE2_INTRIN,
            /*5       */    (EB_ADDDKERNEL_AV1_TYPE)PictureAdditionVoidFunc,
            /*6       */    (EB_ADDDKERNEL_AV1_TYPE)PictureAdditionVoidFunc,
            /*7       */    (EB_ADDDKERNEL_AV1_TYPE)PictureAdditionVoidFunc,
            /*8 64x64 */    PictureAdditionKernel64x64_AV1_SSE2_INTRIN,
        },
    };

    static EB_ADDDKERNEL_TYPE_16BIT FUNC_TABLE AdditionKernel_funcPtrArray16bit[ASM_TYPE_TOTAL] = {
        // NON_AVX2
        PictureAdditionKernel16bit_SSE2_INTRIN,
        // AVX2
        PictureAdditionKernel16bit_SSE2_INTRIN,
    };

    typedef void(*EB_RESDKERNELSUBSAMPLED_TYPE)(
        uint8_t   *input,
        uint32_t   inputStride,
        uint8_t   *pred,
        uint32_t   predStride,
        int16_t  *residual,
        uint32_t   residualStride,
        uint32_t   areaWidth,
        uint32_t   areaHeight,
        uint8_t    lastLine
        );
    static EB_RESDKERNELSUBSAMPLED_TYPE FUNC_TABLE ResidualKernelSubSampled_funcPtrArray[ASM_TYPE_TOTAL][9] = {
        // NON_AVX2
        {
            /*0 4x4  */     ResidualKernelSubSampled4x4_SSE_INTRIN,
            /*1 8x8  */     ResidualKernelSubSampled8x8_SSE2_INTRIN,
            /*2 16x16 */    ResidualKernelSubSampled16x16_SSE2_INTRIN,
            /*3  */         (EB_RESDKERNELSUBSAMPLED_TYPE)PicResdVoidFunc,
            /*4 32x32 */    ResidualKernelSubSampled32x32_SSE2_INTRIN,
            /*5      */     (EB_RESDKERNELSUBSAMPLED_TYPE)PicResdVoidFunc,
            /*6  */         (EB_RESDKERNELSUBSAMPLED_TYPE)PicResdVoidFunc,
            /*7      */     (EB_RESDKERNELSUBSAMPLED_TYPE)PicResdVoidFunc,
            /*8 64x64 */    ResidualKernelSubSampled64x64_SSE2_INTRIN,
        },
        // AVX2
        {
            /*0 4x4  */     ResidualKernelSubSampled4x4_SSE_INTRIN,
            /*1 8x8  */     ResidualKernelSubSampled8x8_SSE2_INTRIN,
            /*2 16x16 */    ResidualKernelSubSampled16x16_SSE2_INTRIN,
            /*3  */         (EB_RESDKERNELSUBSAMPLED_TYPE)PicResdVoidFunc,
            /*4 32x32 */    ResidualKernelSubSampled32x32_SSE2_INTRIN,
            /*5      */     (EB_RESDKERNELSUBSAMPLED_TYPE)PicResdVoidFunc,
            /*6  */         (EB_RESDKERNELSUBSAMPLED_TYPE)PicResdVoidFunc,
            /*7      */     (EB_RESDKERNELSUBSAMPLED_TYPE)PicResdVoidFunc,
            /*8 64x64 */    ResidualKernelSubSampled64x64_SSE2_INTRIN,
        },
    };

    void ResidualKernel16bit(
        uint16_t   *input,
        uint32_t   inputStride,
        uint16_t   *pred,
        uint32_t   predStride,
        int16_t  *residual,
        uint32_t   residualStride,
        uint32_t   areaWidth,
        uint32_t   areaHeight);

    void ResidualKernel_c(
        uint8_t   *input,
        uint32_t   inputStride,
        uint8_t   *pred,
        uint32_t   predStride,
        int16_t  *residual,
        uint32_t   residualStride,
        uint32_t   areaWidth,
        uint32_t   areaHeight);

    static EB_RESDKERNEL_TYPE_16BIT FUNC_TABLE ResidualKernel_funcPtrArray16Bit[ASM_TYPE_TOTAL] = {
        // NON_AVX2
        ResidualKernel16bit_SSE2_INTRIN,
        // AVX2
        ResidualKernel16bit_SSE2_INTRIN,
    };

    static EB_ZEROCOEFF_TYPE FUNC_TABLE PicZeroOutCoef_funcPtrArray[ASM_TYPE_TOTAL][5] = {
        // NON_AVX2
        {
            /*0 4x4   */     ZeroOutCoeff4x4_SSE,
            /*1 8x8   */     ZeroOutCoeff8x8_SSE2,
            /*2 16x16 */     ZeroOutCoeff16x16_SSE2,
            /*3       */     (EB_ZEROCOEFF_TYPE)PicZeroOutCoefVoidFunc,
            /*4 32x32 */     ZeroOutCoeff32x32_SSE2
        },
        // AVX2
        {
            /*0 4x4   */     ZeroOutCoeff4x4_SSE,
            /*1 8x8   */     ZeroOutCoeff8x8_SSE2,
            /*2 16x16 */     ZeroOutCoeff16x16_SSE2,
            /*3       */     (EB_ZEROCOEFF_TYPE)PicZeroOutCoefVoidFunc,
            /*4 32x32 */     ZeroOutCoeff32x32_SSE2
        },
    };

    
    static EB_SATD_TYPE FUNC_TABLE Compute8x8Satd_funcPtrArray[ASM_TYPE_TOTAL] = {
        // NON_AVX2
        Compute8x8Satd_SSE4,
        // ASM_AVX2
        Compute8x8Satd_SSE4
    };

    static EB_SATD_U8_TYPE FUNC_TABLE Compute8x8Satd_U8_funcPtrArray[ASM_TYPE_TOTAL] = {
        // NON_AVX2
        Compute8x8Satd_U8_SSE4,
        // ASM_AVX2
        Compute8x8Satd_U8_SSE4
    };

#if  M0_SPATIAL_SSE || SPATIAL_SSE_I_B_SLICES || M0_SSD_HALF_QUARTER_PEL_BIPRED_SEARCH
    typedef uint64_t(*EB_SPATIALFULLDIST_TYPE)(
        uint8_t   *input,
        uint32_t   inputStride,
        uint8_t   *recon,
        uint32_t   reconStride,
        uint32_t   areaWidth,
        uint32_t   areaHeight);

    static EB_SPATIALFULLDIST_TYPE FUNC_TABLE SpatialFullDistortionKernel_funcPtrArray[ASM_TYPE_TOTAL][5] = {
        // NON_AVX2
        {
            // 4x4
            SpatialFullDistortionKernel4x4_SSSE3_INTRIN,
            // 8x8
            SpatialFullDistortionKernel8x8_SSSE3_INTRIN,
            // 16x16
            SpatialFullDistortionKernel16MxN_SSSE3_INTRIN,
            // 32x32
            SpatialFullDistortionKernel16MxN_SSSE3_INTRIN,
            // 64x64
            SpatialFullDistortionKernel16MxN_SSSE3_INTRIN
        },
        // ASM_AVX2
        {
            // 4x4
            SpatialFullDistortionKernel4x4_SSSE3_INTRIN,
            // 8x8
            SpatialFullDistortionKernel8x8_SSSE3_INTRIN,
            // 16x16
            SpatialFullDistortionKernel16MxN_SSSE3_INTRIN,
            // 32x32
            SpatialFullDistortionKernel16MxN_SSSE3_INTRIN,
            // 64x64
            SpatialFullDistortionKernel16MxN_SSSE3_INTRIN
        },
    };
#endif

    void PictureAdditionKernel16Bit(
        uint16_t  *predPtr,
        uint32_t  predStride,
        int32_t *residual_ptr,
        uint32_t  residualStride,
        uint16_t  *reconPtr,
        uint32_t  reconStride,
        uint32_t  width,
        uint32_t  height,
        int32_t     bd);

#ifdef __cplusplus
}
#endif
#endif // EbPictureOperators_h