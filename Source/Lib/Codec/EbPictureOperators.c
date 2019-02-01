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

/*********************************
 * Includes
 *********************************/

#include "EbPictureOperators.h"
#include "EbPackUnPack.h"
#include <immintrin.h>

#define VARIANCE_PRECISION      16
#define MEAN_PRECISION      (VARIANCE_PRECISION >> 1)

void *aom_memset16(void *dest, int32_t val, size_t length);

/************************************************
* faster copy for <= 64B blocks
************************************************/
#if _WIN32
static void eb_memcpy_small(void* dst_ptr, void const* src_ptr, size_t size)
#else
void eb_memcpy_small(void* dst_ptr, void const* src_ptr, size_t size)
#endif
{
    const char* src = (const char*)src_ptr;
    char*       dst = (char*)dst_ptr;
    size_t      i = 0;

    while ((i + 16) <= size)
    {
        _mm_storeu_ps((float*)(dst + i), _mm_loadu_ps((const float*)(src + i)));
        i += 16;
    }

    if ((i + 8) <= size)
    {
        _mm_store_sd((double*)(dst + i), _mm_load_sd((const double*)(src + i)));
        i += 8;
    }

    for (; i < size; ++i)
        dst[i] = src[i];
}


/************************************************
* mem copy that allows NT loads/stores to avoid cache pollution
************************************************/

#if _WIN32
FORCE_INLINE  void eb_memcpySSE_intel(void* dst_ptr, void const* src_ptr, size_t size)
#else
void eb_memcpySSE_intel(void* dst_ptr, void const* src_ptr, size_t size)
#endif
{
    const char* src = (const char*)src_ptr;
    char*       dst = (char*)dst_ptr;
    size_t      i = 0;
    size_t align_cnt = MIN((64 - ((size_t)dst & 63)), size);



    // align dest to $line
    if (align_cnt != 64)
    {
        eb_memcpy_small(dst, src, align_cnt);
        dst += align_cnt;
        src += align_cnt;
        size -= align_cnt;
    }

    // Copy $line at a time
    // dst aligned to $line
    size_t cline_cnt = (size & ~(size_t)63);
    for (i = 0; i < cline_cnt; i += 64)
    {

        __m128 c0 = _mm_loadu_ps((const float*)(src + i));
        __m128 c1 = _mm_loadu_ps((const float*)(src + i + sizeof(c0)));
        __m128 c2 = _mm_loadu_ps((const float*)(src + i + sizeof(c0) * 2));
        __m128 c3 = _mm_loadu_ps((const float*)(src + i + sizeof(c0) * 3));

        _mm_storeu_ps((float*)(dst + i), c0);
        _mm_storeu_ps((float*)(dst + i + sizeof(c0)), c1);
        _mm_storeu_ps((float*)(dst + i + sizeof(c0) * 2), c2);
        _mm_storeu_ps((float*)(dst + i + sizeof(c0) * 3), c3);

    }

    // Finish tail
    if (i < size)
        eb_memcpy_small(dst + i, src + i, size - i);
}

void  DRDmemcpy(void  *dstPtr, void  *srcPtr, uint32_t  cnt)
{
    if (cnt > 64)
    {
        eb_memcpySSE_intel(dstPtr, srcPtr, cnt);
    }
    else {
        eb_memcpy_small(dstPtr, srcPtr, cnt);
    }
}


/*********************************
 * x86 implememtation of Picture Addition
 *********************************/
void PictureAddition(
    uint8_t  *predPtr,
    uint32_t  predStride,
    int16_t *residual_ptr,
    uint32_t  residualStride,
    uint8_t  *reconPtr,
    uint32_t  reconStride,
    uint32_t  width,
    uint32_t  height,
    EbAsm  asm_type)
{

    AdditionKernel_funcPtrArray[asm_type][width >> 3](
        predPtr,
        predStride,
        residual_ptr,
        residualStride,
        reconPtr,
        reconStride,
        width,
        height
        );

    return;
}

void pic_copy_kernel(
    EbByte                  src,
    uint32_t                   src_stride,
    EbByte                  dst,
    uint32_t                   dst_stride,
    uint32_t                   areaWidth,
    uint32_t                   areaHeight)
{
    uint32_t   j;

    for (j = 0; j < areaHeight; j++)
        memcpy(dst + j * dst_stride, src + j * src_stride, areaWidth);

}
/*********************************
 * Picture Copy 8bit Elements
 *********************************/
EbErrorType PictureCopy8Bit(
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
    EbAsm                   asm_type)
{
    UNUSED(asm_type);
    EbErrorType return_error = EB_ErrorNone;

    // Execute the Kernels
    if (component_mask & PICTURE_BUFFER_DESC_Y_FLAG) {

        pic_copy_kernel(
            &(src->bufferY[srcLumaOriginIndex]),
            src->strideY,
            &(dst->bufferY[dstLumaOriginIndex]),
            dst->strideY,
            areaWidth,
            areaHeight);
    }

    if (component_mask & PICTURE_BUFFER_DESC_Cb_FLAG) {

        pic_copy_kernel(
            &(src->bufferCb[srcChromaOriginIndex]),
            src->strideCb,
            &(dst->bufferCb[dstChromaOriginIndex]),
            dst->strideCb,
            chromaAreaWidth,
            chromaAreaHeight);
    }

    if (component_mask & PICTURE_BUFFER_DESC_Cr_FLAG) {

        pic_copy_kernel(
            &(src->bufferCr[srcChromaOriginIndex]),
            src->strideCr,
            &(dst->bufferCr[dstChromaOriginIndex]),
            dst->strideCr,
            chromaAreaWidth,
            chromaAreaHeight);
    }

    return return_error;
}

/*******************************************
* Residual Kernel 16bit
Computes the residual data
*******************************************/
void ResidualKernel16bit(
    uint16_t   *input,
    uint32_t   inputStride,
    uint16_t   *pred,
    uint32_t   predStride,
    int16_t  *residual,
    uint32_t   residualStride,
    uint32_t   areaWidth,
    uint32_t   areaHeight)
{
    uint32_t  columnIndex;
    uint32_t  rowIndex = 0;

    while (rowIndex < areaHeight) {
        columnIndex = 0;
        while (columnIndex < areaWidth) {
            residual[columnIndex] = ((int16_t)input[columnIndex]) - ((int16_t)pred[columnIndex]);
            ++columnIndex;
        }

        input += inputStride;
        pred += predStride;
        residual += residualStride;
        ++rowIndex;
    }

    return;
}
/*******************************************
* Residual Kernel
Computes the residual data
*******************************************/
void ResidualKernel_c(
    uint8_t   *input,
    uint32_t   inputStride,
    uint8_t   *pred,
    uint32_t   predStride,
    int16_t  *residual,
    uint32_t   residualStride,
    uint32_t   areaWidth,
    uint32_t   areaHeight)
{
    uint32_t  columnIndex;
    uint32_t  rowIndex = 0;

    while (rowIndex < areaHeight) {
        columnIndex = 0;
        while (columnIndex < areaWidth) {
            residual[columnIndex] = ((int16_t)input[columnIndex]) - ((int16_t)pred[columnIndex]);
            ++columnIndex;
        }

        input += inputStride;
        pred += predStride;
        residual += residualStride;
        ++rowIndex;
    }

    return;
}

uint64_t ComputeNxMSatd8x8Units_U8(
    uint8_t  *src,      //int16_t *diff,       // input parameter, diff samples Ptr
    uint32_t  src_stride, //uint32_t  diffStride, // input parameter, source stride
    uint32_t  width,      // input parameter, block width (N)
    uint32_t  height,     // input parameter, block height (M)
    EbAsm  asm_type,
    uint64_t *dcValue)
{
    uint64_t satd = 0;
    uint32_t blockIndexInWidth;
    uint32_t blockIndexInHeight;
    EB_SATD_U8_TYPE Compute8x8SatdFunction = Compute8x8Satd_U8_funcPtrArray[asm_type];

    for (blockIndexInHeight = 0; blockIndexInHeight < height >> 3; ++blockIndexInHeight) {
        for (blockIndexInWidth = 0; blockIndexInWidth < width >> 3; ++blockIndexInWidth) {
            satd += Compute8x8SatdFunction(&(src[(blockIndexInWidth << 3) + (blockIndexInHeight << 3) * src_stride]), dcValue, src_stride);
        }
    }

    return satd;
}


uint64_t ComputeNxMSatd4x4Units_U8(
    uint8_t  *src,       //int16_t *diff,       // input parameter, diff samples Ptr
    uint32_t  src_stride, //uint32_t  diffStride, // input parameter, source stride
    uint32_t  width,      // input parameter, block width (N)
    uint32_t  height,     // input parameter, block height (M)
    uint64_t *dcValue)
{

    uint64_t satd = 0;
    uint32_t blockIndexInWidth;
    uint32_t blockIndexInHeight;

    for (blockIndexInHeight = 0; blockIndexInHeight < height >> 2; ++blockIndexInHeight) {
        for (blockIndexInWidth = 0; blockIndexInWidth < width >> 2; ++blockIndexInWidth) {
            satd += Compute4x4Satd_U8(&(src[(blockIndexInWidth << 2) + (blockIndexInHeight << 2) * src_stride]), dcValue, src_stride);

        }
    }

    return satd;
}
/*******************************************
 *   returns NxM Sum of Absolute Transformed Differences using Compute4x4Satd
 *******************************************/
uint64_t ComputeNxMSatdSadLCU(
    uint8_t  *src,        // input parameter, source samples Ptr
    uint32_t  src_stride,  // input parameter, source stride
    uint32_t  width,      // input parameter, block width (N)
    uint32_t  height,     // input parameter, block height (M)
    EbAsm  asm_type)
{
    uint64_t satd = 0;
    uint64_t  dcValue = 0;
    uint64_t  acValue = 0;

    if (width >= 8 && height >= 8) {
        satd = ComputeNxMSatd8x8Units_U8(
            src,
            src_stride,
            width,
            height,
            asm_type,
            &dcValue);
    }
    else {
        satd =
            ComputeNxMSatd4x4Units_U8(
                src,
                src_stride,
                width,
                height,
                &dcValue);

    }

    acValue = satd - (dcValue >> 2);

    return acValue;
}

/*******************************************
* Picture Full Distortion
*  Used in the Full Mode Decision Loop for the only case of a MVP-SKIP candidate
*******************************************/

void FullDistortionKernel32Bits(
    int32_t  *coeff,
    uint32_t   coeffStride,
    int32_t  *reconCoeff,
    uint32_t   reconCoeffStride,
    uint64_t   distortionResult[DIST_CALC_TOTAL],
    uint32_t   areaWidth,
    uint32_t   areaHeight)
{
    uint32_t  columnIndex;
    uint32_t  rowIndex = 0;
    uint64_t  residualDistortion = 0;
    uint64_t  predictionDistortion = 0;

    while (rowIndex < areaHeight) {

        columnIndex = 0;
        while (columnIndex < areaWidth) {
            residualDistortion += (int64_t)SQR((int64_t)(coeff[columnIndex]) - (reconCoeff[columnIndex]));
            predictionDistortion += (int64_t)SQR((int64_t)(coeff[columnIndex]));
            ++columnIndex;
        }

        coeff += coeffStride;
        reconCoeff += reconCoeffStride;
        ++rowIndex;
    }

    distortionResult[DIST_CALC_RESIDUAL] = residualDistortion;
    distortionResult[DIST_CALC_PREDICTION] = predictionDistortion;
}

/*******************************************
* Picture Distortion Full Kernel CbfZero
*******************************************/
void FullDistortionKernelCbfZero32Bits(
    int32_t  *coeff,
    uint32_t   coeffStride,
    int32_t  *reconCoeff,
    uint32_t   reconCoeffStride,
    uint64_t   distortionResult[DIST_CALC_TOTAL],
    uint32_t   areaWidth,
    uint32_t   areaHeight)
{
    uint32_t  columnIndex;
    uint32_t  rowIndex = 0;
    uint64_t  predictionDistortion = 0;
    (void)reconCoeff;
    (void)reconCoeffStride;

    while (rowIndex < areaHeight) {

        columnIndex = 0;
        while (columnIndex < areaWidth) {
            predictionDistortion += (int64_t)SQR((int64_t)(coeff[columnIndex]));
            ++columnIndex;
        }

        coeff += coeffStride;
        ++rowIndex;
    }

    distortionResult[DIST_CALC_RESIDUAL] = predictionDistortion;
    distortionResult[DIST_CALC_PREDICTION] = predictionDistortion;
}

EbErrorType PictureFullDistortion32Bits(
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
    EbAsm                   asm_type)
{
    EbErrorType return_error = EB_ErrorNone;

    //TODO due to a change in full kernel distortion , ASM has to be updated to not accumulate the input distortion by the output

    if (componentType == COMPONENT_LUMA || componentType == COMPONENT_ALL) {

        y_distortion[0] = 0;
        y_distortion[1] = 0;


        bwidth = bwidth < 64 ? bwidth : 32;
        bheight = bheight < 64 ? bheight : 32;

        if (y_count_non_zero_coeffs) {
            FullDistortionKernel32Bits_funcPtrArray[asm_type](
                &(((int32_t*)coeff->bufferY)[coeffLumaOriginIndex]),
                bwidth,
                &(((int32_t*)reconCoeff->bufferY)[reconCoeffLumaOriginIndex]),
                bwidth,
                y_distortion,
                bwidth,
                bheight);
        }
        else {
            FullDistortionKernelCbfZero32Bits_funcPtrArray[asm_type](
                &(((int32_t*)coeff->bufferY)[coeffLumaOriginIndex]),
                bwidth,
                &(((int32_t*)reconCoeff->bufferY)[reconCoeffLumaOriginIndex]),
                bwidth,
                y_distortion,
                bwidth,
                bheight);
        }
    }

    if (componentType == COMPONENT_CHROMA_CB || componentType == COMPONENT_CHROMA || componentType == COMPONENT_ALL) {
        cb_distortion[0] = 0;
        cb_distortion[1] = 0;

        // CB
        if (cbCountNonZeroCoeffs) {
            FullDistortionKernel32Bits_funcPtrArray[asm_type](
                &(((int32_t*)coeff->bufferCb)[coeffChromaOriginIndex]),
                bwidth_uv,
                &(((int32_t*)reconCoeff->bufferCb)[reconCoeffChromaOriginIndex]),
                bwidth_uv,
                cb_distortion,
                bwidth_uv,
                bheight_uv);
        }
        else {
            FullDistortionKernelCbfZero32Bits_funcPtrArray[asm_type](
                &(((int32_t*)coeff->bufferCb)[coeffChromaOriginIndex]),
                bwidth_uv,
                &(((int32_t*)reconCoeff->bufferCb)[reconCoeffChromaOriginIndex]),
                bwidth_uv,
                cb_distortion,
                bwidth_uv,
                bheight_uv);
        }
    }
    if (componentType == COMPONENT_CHROMA_CR || componentType == COMPONENT_CHROMA || componentType == COMPONENT_ALL) {
        cr_distortion[0] = 0;
        cr_distortion[1] = 0;
        // CR
        if (crCountNonZeroCoeffs) {
            FullDistortionKernel32Bits_funcPtrArray[asm_type](
                &(((int32_t*)coeff->bufferCr)[coeffChromaOriginIndex]),
                bwidth_uv,
                &(((int32_t*)reconCoeff->bufferCr)[reconCoeffChromaOriginIndex]),
                bwidth_uv,
                cr_distortion,
                bwidth_uv,
                bheight_uv);
        }
        else {
            FullDistortionKernelCbfZero32Bits_funcPtrArray[asm_type](
                &(((int32_t*)coeff->bufferCr)[coeffChromaOriginIndex]),
                bwidth_uv,
                &(((int32_t*)reconCoeff->bufferCr)[reconCoeffChromaOriginIndex]),
                bwidth_uv,
                cr_distortion,
                bwidth_uv,
                bheight_uv);
        }
    }

    return return_error;
}

void extract_8bit_data(
    uint16_t      *in16BitBuffer,
    uint32_t       inStride,
    uint8_t       *out8BitBuffer,
    uint32_t       out8Stride,
    uint32_t       width,
    uint32_t       height,
    EbAsm       asm_type
)
{

    UnPack8BIT_funcPtrArray_16Bit[((width & 3) == 0) && ((height & 1) == 0)][asm_type](
        in16BitBuffer,
        inStride,
        out8BitBuffer,
        out8Stride,
        width,
        height);
}
void unpack_l0l1_avg(
    uint16_t *ref16L0,
    uint32_t  refL0Stride,
    uint16_t *ref16L1,
    uint32_t  refL1Stride,
    uint8_t  *dstPtr,
    uint32_t  dst_stride,
    uint32_t  width,
    uint32_t  height,
    EbAsm  asm_type)
{

    UnPackAvg_funcPtrArray[asm_type](
        ref16L0,
        refL0Stride,
        ref16L1,
        refL1Stride,
        dstPtr,
        dst_stride,
        width,
        height);


}
void extract8_bitdata_safe_sub(
    uint16_t      *in16BitBuffer,
    uint32_t       inStride,
    uint8_t       *out8BitBuffer,
    uint32_t       out8Stride,
    uint32_t       width,
    uint32_t       height,
    EbBool      subPred,
    EbAsm       asm_type
)
{

    UnPack8BITSafeSub_funcPtrArray_16Bit[asm_type](
        in16BitBuffer,
        inStride,
        out8BitBuffer,
        out8Stride,
        width,
        height,
        subPred
        );
}
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
    EbAsm  asm_type)
{
    //fix C

    UnPackAvgSafeSub_funcPtrArray[asm_type](
        ref16L0,
        refL0Stride,
        ref16L1,
        refL1Stride,
        dstPtr,
        dst_stride,
        subPred,
        width,
        height);


}
void UnPack2D(
    uint16_t      *in16BitBuffer,
    uint32_t       inStride,
    uint8_t       *out8BitBuffer,
    uint32_t       out8Stride,
    uint8_t       *outnBitBuffer,
    uint32_t       outnStride,
    uint32_t       width,
    uint32_t       height,
    EbAsm       asm_type
)
{

    UnPack2D_funcPtrArray_16Bit[((width & 3) == 0) && ((height & 1) == 0)][asm_type](
        in16BitBuffer,
        inStride,
        out8BitBuffer,
        outnBitBuffer,
        out8Stride,
        outnStride,
        width,
        height);
}

void Pack2D_SRC(
    uint8_t     *in8BitBuffer,
    uint32_t     in8Stride,
    uint8_t     *innBitBuffer,
    uint32_t     innStride,
    uint16_t    *out16BitBuffer,
    uint32_t     outStride,
    uint32_t     width,
    uint32_t     height,
    EbAsm     asm_type
)
{

    Pack2D_funcPtrArray_16Bit_SRC[((width & 3) == 0) && ((height & 1) == 0)][asm_type](
        in8BitBuffer,
        in8Stride,
        innBitBuffer,
        out16BitBuffer,
        innStride,
        outStride,
        width,
        height);
}

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
)
{


    CompressedPack_funcPtrArray[(width == 64 || width == 32) ? asm_type : ASM_NON_AVX2](
        in8BitBuffer,
        in8Stride,
        innBitBuffer,
        out16BitBuffer,
        innStride,
        outStride,
        width,
        height);

}

void Conv2bToCPackLcu(
    const uint8_t     *innBitBuffer,
    uint32_t     innStride,
    uint8_t     *inCompnBitBuffer,
    uint32_t     outStride,
    uint8_t    *localCache,
    uint32_t     width,
    uint32_t     height,
    EbAsm     asm_type)
{

    Convert_Unpack_CPack_funcPtrArray[(width == 64 || width == 32) ? asm_type : ASM_NON_AVX2](
        innBitBuffer,
        innStride,
        inCompnBitBuffer,
        outStride,
        localCache,
        width,
        height);

}


/*******************************************
 * memset16bit
 *******************************************/
void memset16bit(
    uint16_t                     * inPtr,
    uint16_t                       value,
    uint64_t                       numOfElements)
{
    uint64_t i;

    for (i = 0; i < numOfElements; i++) {
        inPtr[i] = value;
    }
}
void memset32bit(
    uint32_t                     * inPtr,
    uint32_t                       value,
    uint64_t                       numOfElements)
{
    uint64_t i;

    for (i = 0; i < numOfElements; i++) {
        inPtr[i] = value;
    }
}


/*******************************************
 * memcpy16bit
 *******************************************/
void memcpy16bit(
    uint16_t                     * outPtr,
    uint16_t                     * inPtr,
    uint64_t                       numOfElements)
{
    uint64_t i;

    for (i = 0; i < numOfElements; i++) {
        outPtr[i] = inPtr[i];
    }
}

/*******************************************
 * memcpy32bit
 *******************************************/
void memcpy32bit(
    uint32_t                     * outPtr,
    uint32_t                     * inPtr,
    uint64_t                       numOfElements)
{
    uint64_t i;

    for (i = 0; i < numOfElements; i++) {
        outPtr[i] = inPtr[i];
    }
}

int32_t  sumResidual(int16_t * inPtr,
    uint32_t   size,
    uint32_t   strideIn)
{

    int32_t sumBlock = 0;
    uint32_t i, j;

    for (j = 0; j < size; j++)
        for (i = 0; i < size; i++)
            sumBlock += inPtr[j*strideIn + i];

    return sumBlock;

}

void memset16bitBlock(
    int16_t * inPtr,
    uint32_t   strideIn,
    uint32_t   size,
    int16_t   value)
{

    uint32_t i;
    for (i = 0; i < size; i++)
        memset16bit((uint16_t*)inPtr + i * strideIn, value, size);

}

static void extend_plane(uint8_t *const src, int32_t src_stride, int32_t width,
    int32_t height, int32_t extend_top, int32_t extend_left,
    int32_t extend_bottom, int32_t extend_right) {
    int32_t i;
    const int32_t linesize = extend_left + extend_right + width;

    /* copy the left and right most columns out */
    uint8_t *src_ptr1 = src;
    uint8_t *src_ptr2 = src + width - 1;
    uint8_t *dst_ptr1 = src - extend_left;
    uint8_t *dst_ptr2 = src + width;

    for (i = 0; i < height; ++i) {
        memset(dst_ptr1, src_ptr1[0], extend_left);
        memset(dst_ptr2, src_ptr2[0], extend_right);
        src_ptr1 += src_stride;
        src_ptr2 += src_stride;
        dst_ptr1 += src_stride;
        dst_ptr2 += src_stride;
    }

    /* Now copy the top and bottom lines into each line of the respective
    * borders
    */
    src_ptr1 = src - extend_left;
    src_ptr2 = src + src_stride * (height - 1) - extend_left;
    dst_ptr1 = src + src_stride * -extend_top - extend_left;
    dst_ptr2 = src + src_stride * height - extend_left;

    for (i = 0; i < extend_top; ++i) {
        memcpy(dst_ptr1, src_ptr1, linesize);
        dst_ptr1 += src_stride;
    }

    for (i = 0; i < extend_bottom; ++i) {
        memcpy(dst_ptr2, src_ptr2, linesize);
        dst_ptr2 += src_stride;
    }
}

static void extend_plane_high(uint8_t *const src8, int32_t src_stride, int32_t width,
    int32_t height, int32_t extend_top, int32_t extend_left,
    int32_t extend_bottom, int32_t extend_right) {
    int32_t i;
    const int32_t linesize = extend_left + extend_right + width;
    uint16_t *src = CONVERT_TO_SHORTPTR(src8);

    /* copy the left and right most columns out */
    uint16_t *src_ptr1 = src;
    uint16_t *src_ptr2 = src + width - 1;
    uint16_t *dst_ptr1 = src - extend_left;
    uint16_t *dst_ptr2 = src + width;

    for (i = 0; i < height; ++i) {
        aom_memset16(dst_ptr1, src_ptr1[0], extend_left);
        aom_memset16(dst_ptr2, src_ptr2[0], extend_right);
        src_ptr1 += src_stride;
        src_ptr2 += src_stride;
        dst_ptr1 += src_stride;
        dst_ptr2 += src_stride;
    }

    /* Now copy the top and bottom lines into each line of the respective
    * borders
    */
    src_ptr1 = src - extend_left;
    src_ptr2 = src + src_stride * (height - 1) - extend_left;
    dst_ptr1 = src + src_stride * -extend_top - extend_left;
    dst_ptr2 = src + src_stride * height - extend_left;

    for (i = 0; i < extend_top; ++i) {
        memcpy(dst_ptr1, src_ptr1, linesize * sizeof(uint16_t));
        dst_ptr1 += src_stride;
    }

    for (i = 0; i < extend_bottom; ++i) {
        memcpy(dst_ptr2, src_ptr2, linesize * sizeof(uint16_t));
        dst_ptr2 += src_stride;
    }
}

void aom_yv12_extend_frame_borders_c(Yv12BufferConfig *ybf,
    const int32_t num_planes) {
    assert(ybf->border % 2 == 0);
    assert(ybf->y_height - ybf->y_crop_height < 16);
    assert(ybf->y_width - ybf->y_crop_width < 16);
    assert(ybf->y_height - ybf->y_crop_height >= 0);
    assert(ybf->y_width - ybf->y_crop_width >= 0);

    if (ybf->flags & YV12_FLAG_HIGHBITDEPTH) {
        for (int32_t plane = 0; plane < num_planes; ++plane) {
            const int32_t is_uv = plane > 0;
            const int32_t plane_border = ybf->border >> is_uv;
            extend_plane_high(
                ybf->buffers[plane], ybf->strides[is_uv], ybf->crop_widths[is_uv],
                ybf->crop_heights[is_uv], plane_border, plane_border,
                plane_border + ybf->heights[is_uv] - ybf->crop_heights[is_uv],
                plane_border + ybf->widths[is_uv] - ybf->crop_widths[is_uv]);
        }
        return;
    }
    for (int32_t plane = 0; plane < num_planes; ++plane) {
        const int32_t is_uv = plane > 0;
        const int32_t plane_border = ybf->border >> is_uv;
        extend_plane(ybf->buffers[plane], ybf->strides[is_uv],
            ybf->crop_widths[is_uv], ybf->crop_heights[is_uv],
            plane_border, plane_border,
            plane_border + ybf->heights[is_uv] - ybf->crop_heights[is_uv],
            plane_border + ybf->widths[is_uv] - ybf->crop_widths[is_uv]);
    }
}



static void memcpy_short_addr(uint8_t *dst8, const uint8_t *src8, int32_t num) {
    uint16_t *dst = CONVERT_TO_SHORTPTR(dst8);
    uint16_t *src = CONVERT_TO_SHORTPTR(src8);
    memcpy(dst, src, num * sizeof(uint16_t));
}

// Copies the source image into the destination image and updates the
// destination's UMV borders.
// Note: The frames are assumed to be identical in size.
void aom_yv12_copy_frame_c(const Yv12BufferConfig *src_bc,
    Yv12BufferConfig *dst_bc, const int32_t num_planes) {

    assert((src_bc->flags & YV12_FLAG_HIGHBITDEPTH) ==
        (dst_bc->flags & YV12_FLAG_HIGHBITDEPTH));

    if (src_bc->flags & YV12_FLAG_HIGHBITDEPTH) {
        for (int32_t plane = 0; plane < num_planes; ++plane) {
            const uint8_t *plane_src = src_bc->buffers[plane];
            uint8_t *plane_dst = dst_bc->buffers[plane];
            const int32_t is_uv = plane > 0;

            for (int32_t row = 0; row < src_bc->heights[is_uv]; ++row) {
                memcpy_short_addr(plane_dst, plane_src, src_bc->widths[is_uv]);
                plane_src += src_bc->strides[is_uv];
                plane_dst += dst_bc->strides[is_uv];
            }
        }
        aom_yv12_extend_frame_borders_c(dst_bc, num_planes);
        return;
    }
    for (int32_t plane = 0; plane < num_planes; ++plane) {
        const uint8_t *plane_src = src_bc->buffers[plane];
        uint8_t *plane_dst = dst_bc->buffers[plane];
        const int32_t is_uv = plane > 0;

        for (int32_t row = 0; row < src_bc->heights[is_uv]; ++row) {
            memcpy(plane_dst, plane_src, src_bc->widths[is_uv]);
            plane_src += src_bc->strides[is_uv];
            plane_dst += dst_bc->strides[is_uv];
        }
    }
    aom_yv12_extend_frame_borders_c(dst_bc, num_planes);
}


void aom_yv12_copy_y_c(const Yv12BufferConfig *src_ybc,
    Yv12BufferConfig *dst_ybc) {
    int32_t row;
    const uint8_t *src = src_ybc->y_buffer;
    uint8_t *dst = dst_ybc->y_buffer;

    if (src_ybc->flags & YV12_FLAG_HIGHBITDEPTH) {
        const uint16_t *src16 = CONVERT_TO_SHORTPTR(src);
        uint16_t *dst16 = CONVERT_TO_SHORTPTR(dst);
        for (row = 0; row < src_ybc->y_height; ++row) {
            memcpy(dst16, src16, src_ybc->y_width * sizeof(uint16_t));
            src16 += src_ybc->y_stride;
            dst16 += dst_ybc->y_stride;
        }
        return;
    }

    for (row = 0; row < src_ybc->y_height; ++row) {
        memcpy(dst, src, src_ybc->y_width);
        src += src_ybc->y_stride;
        dst += dst_ybc->y_stride;
    }
}

void aom_yv12_copy_u_c(const Yv12BufferConfig *src_bc,
    Yv12BufferConfig *dst_bc) {
    int32_t row;
    const uint8_t *src = src_bc->u_buffer;
    uint8_t *dst = dst_bc->u_buffer;

    if (src_bc->flags & YV12_FLAG_HIGHBITDEPTH) {
        const uint16_t *src16 = CONVERT_TO_SHORTPTR(src);
        uint16_t *dst16 = CONVERT_TO_SHORTPTR(dst);
        for (row = 0; row < src_bc->uv_height; ++row) {
            memcpy(dst16, src16, src_bc->uv_width * sizeof(uint16_t));
            src16 += src_bc->uv_stride;
            dst16 += dst_bc->uv_stride;
        }
        return;
    }

    for (row = 0; row < src_bc->uv_height; ++row) {
        memcpy(dst, src, src_bc->uv_width);
        src += src_bc->uv_stride;
        dst += dst_bc->uv_stride;
    }
}

void aom_yv12_copy_v_c(const Yv12BufferConfig *src_bc,
    Yv12BufferConfig *dst_bc) {
    int32_t row;
    const uint8_t *src = src_bc->v_buffer;
    uint8_t *dst = dst_bc->v_buffer;

    if (src_bc->flags & YV12_FLAG_HIGHBITDEPTH) {
        const uint16_t *src16 = CONVERT_TO_SHORTPTR(src);
        uint16_t *dst16 = CONVERT_TO_SHORTPTR(dst);
        for (row = 0; row < src_bc->uv_height; ++row) {
            memcpy(dst16, src16, src_bc->uv_width * sizeof(uint16_t));
            src16 += src_bc->uv_stride;
            dst16 += dst_bc->uv_stride;
        }
        return;
    }

    for (row = 0; row < src_bc->uv_height; ++row) {
        memcpy(dst, src, src_bc->uv_width);
        src += src_bc->uv_stride;
        dst += dst_bc->uv_stride;
    }
}

