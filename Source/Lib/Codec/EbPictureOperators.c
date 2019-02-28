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

void  DRDmemcpy(void  *dst_ptr, void  *src_ptr, uint32_t  cnt)
{
    if (cnt > 64)
    {
        eb_memcpySSE_intel(dst_ptr, src_ptr, cnt);
    }
    else {
        eb_memcpy_small(dst_ptr, src_ptr, cnt);
    }
}


/*********************************
 * x86 implememtation of Picture Addition
 *********************************/
void picture_addition(
    uint8_t  *pred_ptr,
    uint32_t  pred_stride,
    int16_t *residual_ptr,
    uint32_t  residual_stride,
    uint8_t  *recon_ptr,
    uint32_t  recon_stride,
    uint32_t  width,
    uint32_t  height,
    EbAsm  asm_type)
{

    addition_kernel_func_ptr_array[asm_type][width >> 3](
        pred_ptr,
        pred_stride,
        residual_ptr,
        residual_stride,
        recon_ptr,
        recon_stride,
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
    uint32_t                   area_width,
    uint32_t                   area_height)
{
    uint32_t   j;

    for (j = 0; j < area_height; j++)
        memcpy(dst + j * dst_stride, src + j * src_stride, area_width);

}
/*********************************
 * Picture Copy 8bit Elements
 *********************************/
EbErrorType picture_copy8_bit(
    EbPictureBufferDesc_t   *src,
    uint32_t                   src_luma_origin_index,
    uint32_t                   src_chroma_origin_index,
    EbPictureBufferDesc_t   *dst,
    uint32_t                   dst_luma_origin_index,
    uint32_t                   dst_chroma_origin_index,
    uint32_t                   area_width,
    uint32_t                   area_height,
    uint32_t                   chroma_area_width,
    uint32_t                   chroma_area_height,
    uint32_t                   component_mask,
    EbAsm                   asm_type)
{
    UNUSED(asm_type);
    EbErrorType return_error = EB_ErrorNone;

    // Execute the Kernels
    if (component_mask & PICTURE_BUFFER_DESC_Y_FLAG) {

        pic_copy_kernel(
            &(src->buffer_y[src_luma_origin_index]),
            src->stride_y,
            &(dst->buffer_y[dst_luma_origin_index]),
            dst->stride_y,
            area_width,
            area_height);
    }

    if (component_mask & PICTURE_BUFFER_DESC_Cb_FLAG) {

        pic_copy_kernel(
            &(src->bufferCb[src_chroma_origin_index]),
            src->strideCb,
            &(dst->bufferCb[dst_chroma_origin_index]),
            dst->strideCb,
            chroma_area_width,
            chroma_area_height);
    }

    if (component_mask & PICTURE_BUFFER_DESC_Cr_FLAG) {

        pic_copy_kernel(
            &(src->bufferCr[src_chroma_origin_index]),
            src->strideCr,
            &(dst->bufferCr[dst_chroma_origin_index]),
            dst->strideCr,
            chroma_area_width,
            chroma_area_height);
    }

    return return_error;
}

/*******************************************
* Residual Kernel 16bit
Computes the residual data
*******************************************/
void residual_kernel16bit(
    uint16_t   *input,
    uint32_t   input_stride,
    uint16_t   *pred,
    uint32_t   pred_stride,
    int16_t  *residual,
    uint32_t   residual_stride,
    uint32_t   area_width,
    uint32_t   area_height)
{
    uint32_t  columnIndex;
    uint32_t  row_index = 0;

    while (row_index < area_height) {
        columnIndex = 0;
        while (columnIndex < area_width) {
            residual[columnIndex] = ((int16_t)input[columnIndex]) - ((int16_t)pred[columnIndex]);
            ++columnIndex;
        }

        input += input_stride;
        pred += pred_stride;
        residual += residual_stride;
        ++row_index;
    }

    return;
}
/*******************************************
* Residual Kernel
Computes the residual data
*******************************************/
void residual_kernel_c(
    uint8_t   *input,
    uint32_t   input_stride,
    uint8_t   *pred,
    uint32_t   pred_stride,
    int16_t  *residual,
    uint32_t   residual_stride,
    uint32_t   area_width,
    uint32_t   area_height)
{
    uint32_t  columnIndex;
    uint32_t  row_index = 0;

    while (row_index < area_height) {
        columnIndex = 0;
        while (columnIndex < area_width) {
            residual[columnIndex] = ((int16_t)input[columnIndex]) - ((int16_t)pred[columnIndex]);
            ++columnIndex;
        }

        input += input_stride;
        pred += pred_stride;
        residual += residual_stride;
        ++row_index;
    }

    return;
}

uint64_t ComputeNxMSatd8x8Units_U8(
    uint8_t  *src,      //int16_t *diff,       // input parameter, diff samples Ptr
    uint32_t  src_stride, //uint32_t  diffStride, // input parameter, source stride
    uint32_t  width,      // input parameter, block width (N)
    uint32_t  height,     // input parameter, block height (M)
    EbAsm  asm_type,
    uint64_t *dc_value)
{
    uint64_t satd = 0;
    uint32_t blockIndexInWidth;
    uint32_t blockIndexInHeight;
    EB_SATD_U8_TYPE Compute8x8SatdFunction = compute8x8_satd_u8_func_ptr_array[asm_type];

    for (blockIndexInHeight = 0; blockIndexInHeight < height >> 3; ++blockIndexInHeight) {
        for (blockIndexInWidth = 0; blockIndexInWidth < width >> 3; ++blockIndexInWidth) {
            satd += Compute8x8SatdFunction(&(src[(blockIndexInWidth << 3) + (blockIndexInHeight << 3) * src_stride]), dc_value, src_stride);
        }
    }

    return satd;
}


uint64_t ComputeNxMSatd4x4Units_U8(
    uint8_t  *src,       //int16_t *diff,       // input parameter, diff samples Ptr
    uint32_t  src_stride, //uint32_t  diffStride, // input parameter, source stride
    uint32_t  width,      // input parameter, block width (N)
    uint32_t  height,     // input parameter, block height (M)
    uint64_t *dc_value)
{

    uint64_t satd = 0;
    uint32_t blockIndexInWidth;
    uint32_t blockIndexInHeight;

    for (blockIndexInHeight = 0; blockIndexInHeight < height >> 2; ++blockIndexInHeight) {
        for (blockIndexInWidth = 0; blockIndexInWidth < width >> 2; ++blockIndexInWidth) {
            satd += Compute4x4Satd_U8(&(src[(blockIndexInWidth << 2) + (blockIndexInHeight << 2) * src_stride]), dc_value, src_stride);

        }
    }

    return satd;
}
/*******************************************
 *   returns NxM Sum of Absolute Transformed Differences using Compute4x4Satd
 *******************************************/
uint64_t compute_nx_m_satd_sad_lcu(
    uint8_t  *src,        // input parameter, source samples Ptr
    uint32_t  src_stride,  // input parameter, source stride
    uint32_t  width,      // input parameter, block width (N)
    uint32_t  height,     // input parameter, block height (M)
    EbAsm  asm_type)
{
    uint64_t satd = 0;
    uint64_t  dc_value = 0;
    uint64_t  acValue = 0;

    if (width >= 8 && height >= 8) {
        satd = ComputeNxMSatd8x8Units_U8(
            src,
            src_stride,
            width,
            height,
            asm_type,
            &dc_value);
    }
    else {
        satd =
            ComputeNxMSatd4x4Units_U8(
                src,
                src_stride,
                width,
                height,
                &dc_value);

    }

    acValue = satd - (dc_value >> 2);

    return acValue;
}

/*******************************************
* Picture Full Distortion
*  Used in the Full Mode Decision Loop for the only case of a MVP-SKIP candidate
*******************************************/

void full_distortion_kernel32_bits(
    int32_t  *coeff,
    uint32_t   coeff_stride,
    int32_t  *recon_coeff,
    uint32_t   recon_coeff_stride,
    uint64_t   distortion_result[DIST_CALC_TOTAL],
    uint32_t   area_width,
    uint32_t   area_height)
{
    uint32_t  columnIndex;
    uint32_t  row_index = 0;
    uint64_t  residualDistortion = 0;
    uint64_t  predictionDistortion = 0;

    while (row_index < area_height) {

        columnIndex = 0;
        while (columnIndex < area_width) {
            residualDistortion += (int64_t)SQR((int64_t)(coeff[columnIndex]) - (recon_coeff[columnIndex]));
            predictionDistortion += (int64_t)SQR((int64_t)(coeff[columnIndex]));
            ++columnIndex;
        }

        coeff += coeff_stride;
        recon_coeff += recon_coeff_stride;
        ++row_index;
    }

    distortion_result[DIST_CALC_RESIDUAL] = residualDistortion;
    distortion_result[DIST_CALC_PREDICTION] = predictionDistortion;
}

/*******************************************
* Picture Distortion Full Kernel CbfZero
*******************************************/
void full_distortion_kernel_cbf_zero32_bits(
    int32_t  *coeff,
    uint32_t   coeff_stride,
    int32_t  *recon_coeff,
    uint32_t   recon_coeff_stride,
    uint64_t   distortion_result[DIST_CALC_TOTAL],
    uint32_t   area_width,
    uint32_t   area_height)
{
    uint32_t  columnIndex;
    uint32_t  row_index = 0;
    uint64_t  predictionDistortion = 0;
    (void)recon_coeff;
    (void)recon_coeff_stride;

    while (row_index < area_height) {

        columnIndex = 0;
        while (columnIndex < area_width) {
            predictionDistortion += (int64_t)SQR((int64_t)(coeff[columnIndex]));
            ++columnIndex;
        }

        coeff += coeff_stride;
        ++row_index;
    }

    distortion_result[DIST_CALC_RESIDUAL] = predictionDistortion;
    distortion_result[DIST_CALC_PREDICTION] = predictionDistortion;
}

EbErrorType picture_full_distortion32_bits(
    EbPictureBufferDesc_t   *coeff,
    uint32_t                   coeff_luma_origin_index,
    uint32_t                   coeff_chroma_origin_index,
    EbPictureBufferDesc_t   *recon_coeff,
    uint32_t                   recon_coeff_luma_origin_index,
    uint32_t                   recon_coeff_chroma_origin_index,
    uint32_t                   bwidth,
    uint32_t                   bheight,
    uint32_t                   bwidth_uv,
    uint32_t                   bheight_uv,
    uint64_t                   y_distortion[DIST_CALC_TOTAL],
    uint64_t                   cb_distortion[DIST_CALC_TOTAL],
    uint64_t                   cr_distortion[DIST_CALC_TOTAL],
    uint32_t                   y_count_non_zero_coeffs,
    uint32_t                   cb_count_non_zero_coeffs,
    uint32_t                   cr_count_non_zero_coeffs,
    COMPONENT_TYPE            component_type,
    EbAsm                   asm_type)
{
    EbErrorType return_error = EB_ErrorNone;

    //TODO due to a change in full kernel distortion , ASM has to be updated to not accumulate the input distortion by the output

    if (component_type == COMPONENT_LUMA || component_type == COMPONENT_ALL) {

        y_distortion[0] = 0;
        y_distortion[1] = 0;


        bwidth = bwidth < 64 ? bwidth : 32;
        bheight = bheight < 64 ? bheight : 32;

        if (y_count_non_zero_coeffs) {
            full_distortion_kernel32_bits_func_ptr_array[asm_type](
                &(((int32_t*)coeff->buffer_y)[coeff_luma_origin_index]),
                bwidth,
                &(((int32_t*)recon_coeff->buffer_y)[recon_coeff_luma_origin_index]),
                bwidth,
                y_distortion,
                bwidth,
                bheight);
        }
        else {
            full_distortion_kernel_cbf_zero32_bits_func_ptr_array[asm_type](
                &(((int32_t*)coeff->buffer_y)[coeff_luma_origin_index]),
                bwidth,
                &(((int32_t*)recon_coeff->buffer_y)[recon_coeff_luma_origin_index]),
                bwidth,
                y_distortion,
                bwidth,
                bheight);
        }
    }

    if (component_type == COMPONENT_CHROMA_CB || component_type == COMPONENT_CHROMA || component_type == COMPONENT_ALL) {
        cb_distortion[0] = 0;
        cb_distortion[1] = 0;

        // CB
        if (cb_count_non_zero_coeffs) {
            full_distortion_kernel32_bits_func_ptr_array[asm_type](
                &(((int32_t*)coeff->bufferCb)[coeff_chroma_origin_index]),
                bwidth_uv,
                &(((int32_t*)recon_coeff->bufferCb)[recon_coeff_chroma_origin_index]),
                bwidth_uv,
                cb_distortion,
                bwidth_uv,
                bheight_uv);
        }
        else {
            full_distortion_kernel_cbf_zero32_bits_func_ptr_array[asm_type](
                &(((int32_t*)coeff->bufferCb)[coeff_chroma_origin_index]),
                bwidth_uv,
                &(((int32_t*)recon_coeff->bufferCb)[recon_coeff_chroma_origin_index]),
                bwidth_uv,
                cb_distortion,
                bwidth_uv,
                bheight_uv);
        }
    }
    if (component_type == COMPONENT_CHROMA_CR || component_type == COMPONENT_CHROMA || component_type == COMPONENT_ALL) {
        cr_distortion[0] = 0;
        cr_distortion[1] = 0;
        // CR
        if (cr_count_non_zero_coeffs) {
            full_distortion_kernel32_bits_func_ptr_array[asm_type](
                &(((int32_t*)coeff->bufferCr)[coeff_chroma_origin_index]),
                bwidth_uv,
                &(((int32_t*)recon_coeff->bufferCr)[recon_coeff_chroma_origin_index]),
                bwidth_uv,
                cr_distortion,
                bwidth_uv,
                bheight_uv);
        }
        else {
            full_distortion_kernel_cbf_zero32_bits_func_ptr_array[asm_type](
                &(((int32_t*)coeff->bufferCr)[coeff_chroma_origin_index]),
                bwidth_uv,
                &(((int32_t*)recon_coeff->bufferCr)[recon_coeff_chroma_origin_index]),
                bwidth_uv,
                cr_distortion,
                bwidth_uv,
                bheight_uv);
        }
    }

    return return_error;
}

void extract_8bit_data(
    uint16_t      *in16_bit_buffer,
    uint32_t       in_stride,
    uint8_t       *out8_bit_buffer,
    uint32_t       out8_stride,
    uint32_t       width,
    uint32_t       height,
    EbAsm       asm_type
)
{

    unpack8_bit_func_ptr_array_16_bit[((width & 3) == 0) && ((height & 1) == 0)][asm_type](
        in16_bit_buffer,
        in_stride,
        out8_bit_buffer,
        out8_stride,
        width,
        height);
}
void unpack_l0l1_avg(
    uint16_t *ref16_l0,
    uint32_t  ref_l0_stride,
    uint16_t *ref16_l1,
    uint32_t  ref_l1_stride,
    uint8_t  *dst_ptr,
    uint32_t  dst_stride,
    uint32_t  width,
    uint32_t  height,
    EbAsm  asm_type)
{

    un_pack_avg_func_ptr_array[asm_type](
        ref16_l0,
        ref_l0_stride,
        ref16_l1,
        ref_l1_stride,
        dst_ptr,
        dst_stride,
        width,
        height);


}
void extract8_bitdata_safe_sub(
    uint16_t      *in16_bit_buffer,
    uint32_t       in_stride,
    uint8_t       *out8_bit_buffer,
    uint32_t       out8_stride,
    uint32_t       width,
    uint32_t       height,
    EbBool      sub_pred,
    EbAsm       asm_type
)
{
    /* sub_pred not implemented */
    (void)sub_pred;

    unpack8_bit_safe_sub_func_ptr_array_16_bit[asm_type](
        in16_bit_buffer,
        in_stride,
        out8_bit_buffer,
        out8_stride,
        width,
        height
        );
}
void unpack_l0l1_avg_safe_sub(
    uint16_t *ref16_l0,
    uint32_t  ref_l0_stride,
    uint16_t *ref16_l1,
    uint32_t  ref_l1_stride,
    uint8_t  *dst_ptr,
    uint32_t  dst_stride,
    uint32_t  width,
    uint32_t  height,
    EbBool      sub_pred,
    EbAsm  asm_type)
{
    //fix C

    unpack_avg_safe_sub_func_ptr_array[asm_type](
        ref16_l0,
        ref_l0_stride,
        ref16_l1,
        ref_l1_stride,
        dst_ptr,
        dst_stride,
        sub_pred,
        width,
        height);


}
void un_pack2d(
    uint16_t      *in16_bit_buffer,
    uint32_t       in_stride,
    uint8_t       *out8_bit_buffer,
    uint32_t       out8_stride,
    uint8_t       *outn_bit_buffer,
    uint32_t       outn_stride,
    uint32_t       width,
    uint32_t       height,
    EbAsm       asm_type
)
{

    un_pack2d_func_ptr_array_16_bit[((width & 3) == 0) && ((height & 1) == 0)][asm_type](
        in16_bit_buffer,
        in_stride,
        out8_bit_buffer,
        outn_bit_buffer,
        out8_stride,
        outn_stride,
        width,
        height);
}

void pack2d_src(
    uint8_t     *in8_bit_buffer,
    uint32_t     in8_stride,
    uint8_t     *inn_bit_buffer,
    uint32_t     inn_stride,
    uint16_t    *out16_bit_buffer,
    uint32_t     out_stride,
    uint32_t     width,
    uint32_t     height,
    EbAsm     asm_type
)
{

    pack2d_func_ptr_array_16_bit_src[((width & 3) == 0) && ((height & 1) == 0)][asm_type](
        in8_bit_buffer,
        in8_stride,
        inn_bit_buffer,
        out16_bit_buffer,
        inn_stride,
        out_stride,
        width,
        height);
}

void compressed_pack_lcu(
    uint8_t     *in8_bit_buffer,
    uint32_t     in8_stride,
    uint8_t     *inn_bit_buffer,
    uint32_t     inn_stride,
    uint16_t    *out16_bit_buffer,
    uint32_t     out_stride,
    uint32_t     width,
    uint32_t     height,
    EbAsm     asm_type
)
{


    compressed_pack_func_ptr_array[(width == 64 || width == 32) ? asm_type : ASM_NON_AVX2](
        in8_bit_buffer,
        in8_stride,
        inn_bit_buffer,
        out16_bit_buffer,
        inn_stride,
        out_stride,
        width,
        height);

}

void conv2b_to_c_pack_lcu(
    const uint8_t     *inn_bit_buffer,
    uint32_t     inn_stride,
    uint8_t     *in_compn_bit_buffer,
    uint32_t     out_stride,
    uint8_t    *local_cache,
    uint32_t     width,
    uint32_t     height,
    EbAsm     asm_type)
{

    convert_unpack_c_pack_func_ptr_array[(width == 64 || width == 32) ? asm_type : ASM_NON_AVX2](
        inn_bit_buffer,
        inn_stride,
        in_compn_bit_buffer,
        out_stride,
        local_cache,
        width,
        height);

}


/*******************************************
 * memset16bit
 *******************************************/
void memset16bit(
    uint16_t                     * in_ptr,
    uint16_t                       value,
    uint64_t                       num_of_elements)
{
    uint64_t i;

    for (i = 0; i < num_of_elements; i++) {
        in_ptr[i] = value;
    }
}
void memset32bit(
    uint32_t                     * in_ptr,
    uint32_t                       value,
    uint64_t                       num_of_elements)
{
    uint64_t i;

    for (i = 0; i < num_of_elements; i++) {
        in_ptr[i] = value;
    }
}


/*******************************************
 * memcpy16bit
 *******************************************/
void memcpy16bit(
    uint16_t                     * out_ptr,
    uint16_t                     * in_ptr,
    uint64_t                       num_of_elements)
{
    uint64_t i;

    for (i = 0; i < num_of_elements; i++) {
        out_ptr[i] = in_ptr[i];
    }
}

/*******************************************
 * memcpy32bit
 *******************************************/
void memcpy32bit(
    uint32_t                     * out_ptr,
    uint32_t                     * in_ptr,
    uint64_t                       num_of_elements)
{
    uint64_t i;

    for (i = 0; i < num_of_elements; i++) {
        out_ptr[i] = in_ptr[i];
    }
}

int32_t  sum_residual(int16_t * in_ptr,
    uint32_t   size,
    uint32_t   stride_in)
{

    int32_t sumBlock = 0;
    uint32_t i, j;

    for (j = 0; j < size; j++)
        for (i = 0; i < size; i++)
            sumBlock += in_ptr[j*stride_in + i];

    return sumBlock;

}

void memset16bit_block(
    int16_t * in_ptr,
    uint32_t   stride_in,
    uint32_t   size,
    int16_t   value)
{

    uint32_t i;
    for (i = 0; i < size; i++)
        memset16bit((uint16_t*)in_ptr + i * stride_in, value, size);

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

