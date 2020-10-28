/*
* Copyright(c) 2019 Intel Corporation
* Copyright (c) 2016, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

/*********************************
 * Includes
 *********************************/

#include "EbPictureOperators.h"
#include "EbPackUnPack.h"
#include "common_dsp_rtcd.h"
#include "EbUtility.h"

#define VARIANCE_PRECISION 16
#define MEAN_PRECISION (VARIANCE_PRECISION >> 1)

void *svt_aom_memset16(void *dest, int32_t val, size_t length);

/*********************************
 * Picture Copy
 *********************************/

void pic_copy_kernel_8bit(EbByte src, uint32_t src_stride, EbByte dst, uint32_t dst_stride,
                          uint32_t area_width, uint32_t area_height) {
    for (uint32_t j = 0; j < area_height; j++)
        svt_memcpy(dst + j * dst_stride, src + j * src_stride, area_width);
}
void pic_copy_kernel_16bit(uint16_t *src, uint32_t src_stride, uint16_t *dst, uint32_t dst_stride,
                           uint32_t width, uint32_t height) {
    for (uint32_t j = 0; j < height; j++)
        svt_memcpy(dst + j * dst_stride, src + j * src_stride, sizeof(uint16_t) * width);
}

EbErrorType picture_copy(EbPictureBufferDesc *src, uint32_t src_luma_origin_index,
                         uint32_t src_chroma_origin_index, EbPictureBufferDesc *dst,
                         uint32_t dst_luma_origin_index, uint32_t dst_chroma_origin_index,
                         uint32_t area_width, uint32_t area_height, uint32_t chroma_area_width,
                         uint32_t chroma_area_height, uint32_t component_mask, EbBool hbd) {
    EbErrorType return_error = EB_ErrorNone;

    if (hbd) {
        if (component_mask & PICTURE_BUFFER_DESC_Y_FLAG)
            pic_copy_kernel_16bit(((uint16_t *)src->buffer_y) + src_luma_origin_index,
                                  src->stride_y,
                                  ((uint16_t *)dst->buffer_y) + dst_luma_origin_index,
                                  dst->stride_y,
                                  area_width,
                                  area_height);

        if (component_mask & PICTURE_BUFFER_DESC_Cb_FLAG)
            pic_copy_kernel_16bit(((uint16_t *)src->buffer_cb) + src_chroma_origin_index,
                                  src->stride_cb,
                                  ((uint16_t *)dst->buffer_cb) + dst_chroma_origin_index,
                                  dst->stride_cb,
                                  chroma_area_width,
                                  chroma_area_height);

        if (component_mask & PICTURE_BUFFER_DESC_Cr_FLAG)
            pic_copy_kernel_16bit(((uint16_t *)src->buffer_cr) + src_chroma_origin_index,
                                  src->stride_cr,
                                  ((uint16_t *)dst->buffer_cr) + dst_chroma_origin_index,
                                  dst->stride_cr,
                                  chroma_area_width,
                                  chroma_area_height);
    } else {
        if (component_mask & PICTURE_BUFFER_DESC_Y_FLAG)
            pic_copy_kernel_8bit(&(src->buffer_y[src_luma_origin_index]),
                                 src->stride_y,
                                 &(dst->buffer_y[dst_luma_origin_index]),
                                 dst->stride_y,
                                 area_width,
                                 area_height);

        if (component_mask & PICTURE_BUFFER_DESC_Cb_FLAG)
            pic_copy_kernel_8bit(&(src->buffer_cb[src_chroma_origin_index]),
                                 src->stride_cb,
                                 &(dst->buffer_cb[dst_chroma_origin_index]),
                                 dst->stride_cb,
                                 chroma_area_width,
                                 chroma_area_height);

        if (component_mask & PICTURE_BUFFER_DESC_Cr_FLAG)
            pic_copy_kernel_8bit(&(src->buffer_cr[src_chroma_origin_index]),
                                 src->stride_cr,
                                 &(dst->buffer_cr[dst_chroma_origin_index]),
                                 dst->stride_cr,
                                 chroma_area_width,
                                 chroma_area_height);
    }

    return return_error;
}

/*******************************************
* Residual Kernel 16bit
Computes the residual data
*******************************************/
void svt_residual_kernel16bit_c(uint16_t *input, uint32_t input_stride, uint16_t *pred,
                                uint32_t pred_stride, int16_t *residual, uint32_t residual_stride,
                                uint32_t area_width, uint32_t area_height) {
    uint32_t row_index = 0;

    while (row_index < area_height) {
        uint32_t column_index = 0;
        while (column_index < area_width) {
            residual[column_index] = ((int16_t)input[column_index]) - ((int16_t)pred[column_index]);
            ++column_index;
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
void svt_residual_kernel8bit_c(uint8_t *input, uint32_t input_stride, uint8_t *pred,
                               uint32_t pred_stride, int16_t *residual, uint32_t residual_stride,
                               uint32_t area_width, uint32_t area_height) {
    uint32_t row_index = 0;

    while (row_index < area_height) {
        uint32_t column_index = 0;
        while (column_index < area_width) {
            residual[column_index] = ((int16_t)input[column_index]) - ((int16_t)pred[column_index]);
            ++column_index;
        }

        input += input_stride;
        pred += pred_stride;
        residual += residual_stride;
        ++row_index;
    }

    return;
}

/*******************************************
* Picture Full Distortion
*  Used in the Full Mode Decision Loop for the only case of a MVP-SKIP candidate
*******************************************/

void svt_full_distortion_kernel32_bits_c(int32_t *coeff, uint32_t coeff_stride,
                                         int32_t *recon_coeff, uint32_t recon_coeff_stride,
                                         uint64_t distortion_result[DIST_CALC_TOTAL],
                                         uint32_t area_width, uint32_t area_height) {
    uint32_t row_index             = 0;
    uint64_t residual_distortion   = 0;
    uint64_t prediction_distortion = 0;

    while (row_index < area_height) {
        uint32_t column_index = 0;
        while (column_index < area_width) {
            residual_distortion +=
                (int64_t)SQR((int64_t)(coeff[column_index]) - (recon_coeff[column_index]));
            prediction_distortion += (int64_t)SQR((int64_t)(coeff[column_index]));
            ++column_index;
        }

        coeff += coeff_stride;
        recon_coeff += recon_coeff_stride;
        ++row_index;
    }

    distortion_result[DIST_CALC_RESIDUAL]   = residual_distortion;
    distortion_result[DIST_CALC_PREDICTION] = prediction_distortion;
}

uint64_t svt_full_distortion_kernel16_bits_c(uint8_t *input, uint32_t input_offset,
                                             uint32_t input_stride, uint8_t *pred,
                                             int32_t pred_offset, uint32_t pred_stride,
                                             uint32_t area_width, uint32_t area_height) {
    uint32_t row_index      = 0;
    uint64_t sse_distortion = 0;

    uint16_t *input_16bit = (uint16_t *)input;
    uint16_t *pred_16bit  = (uint16_t *)pred;
    input_16bit += input_offset;
    pred_16bit += pred_offset;

    while (row_index < area_height) {
        uint32_t column_index = 0;
        while (column_index < area_width) {
            sse_distortion +=
                (int64_t)SQR((int64_t)(input_16bit[column_index]) - (pred_16bit[column_index]));
            ++column_index;
        }
        input_16bit += input_stride;
        pred_16bit += pred_stride;
        ++row_index;
    }

    return sse_distortion;
}

/*******************************************
* Picture Distortion Full Kernel CbfZero
*******************************************/
void svt_full_distortion_kernel_cbf_zero32_bits_c(int32_t *coeff, uint32_t coeff_stride,
                                                  uint64_t distortion_result[DIST_CALC_TOTAL],
                                                  uint32_t area_width, uint32_t area_height) {
    uint32_t row_index             = 0;
    uint64_t prediction_distortion = 0;

    while (row_index < area_height) {
        uint32_t column_index = 0;
        while (column_index < area_width) {
            prediction_distortion += (int64_t)SQR((int64_t)(coeff[column_index]));
            ++column_index;
        }

        coeff += coeff_stride;
        ++row_index;
    }

    distortion_result[DIST_CALC_RESIDUAL]   = prediction_distortion;
    distortion_result[DIST_CALC_PREDICTION] = prediction_distortion;
}

EbErrorType picture_full_distortion32_bits(
    EbPictureBufferDesc *coeff, uint32_t coeff_luma_origin_index,
    uint32_t coeff_chroma_origin_index, EbPictureBufferDesc *recon_coeff,
    uint32_t recon_coeff_luma_origin_index, uint32_t recon_coeff_chroma_origin_index,
    uint32_t bwidth, uint32_t bheight, uint32_t bwidth_uv, uint32_t bheight_uv,
    uint64_t y_distortion[DIST_CALC_TOTAL], uint64_t cb_distortion[DIST_CALC_TOTAL],
    uint64_t cr_distortion[DIST_CALC_TOTAL], uint32_t y_count_non_zero_coeffs,
    uint32_t cb_count_non_zero_coeffs, uint32_t cr_count_non_zero_coeffs,
    COMPONENT_TYPE component_type) {
    EbErrorType return_error = EB_ErrorNone;

    //TODO due to a change in full kernel distortion , ASM has to be updated to not accumulate the input distortion by the output

    if (component_type == COMPONENT_LUMA || component_type == COMPONENT_ALL) {
        y_distortion[0] = 0;
        y_distortion[1] = 0;

        bwidth  = bwidth < 64 ? bwidth : 32;
        bheight = bheight < 64 ? bheight : 32;

        if (y_count_non_zero_coeffs) {
            svt_full_distortion_kernel32_bits(
                &(((int32_t *)coeff->buffer_y)[coeff_luma_origin_index]),
                bwidth,
                &(((int32_t *)recon_coeff->buffer_y)[recon_coeff_luma_origin_index]),
                bwidth,
                y_distortion,
                bwidth,
                bheight);
        } else {
            svt_full_distortion_kernel_cbf_zero32_bits(
                &(((int32_t *)coeff->buffer_y)[coeff_luma_origin_index]),
                bwidth,
                y_distortion,
                bwidth,
                bheight);
        }
    }

    if (component_type == COMPONENT_CHROMA_CB || component_type == COMPONENT_CHROMA ||
        component_type == COMPONENT_ALL) {
        cb_distortion[0] = 0;
        cb_distortion[1] = 0;

        // CB
        if (cb_count_non_zero_coeffs) {
            svt_full_distortion_kernel32_bits(
                &(((int32_t *)coeff->buffer_cb)[coeff_chroma_origin_index]),
                bwidth_uv,
                &(((int32_t *)recon_coeff->buffer_cb)[recon_coeff_chroma_origin_index]),
                bwidth_uv,
                cb_distortion,
                bwidth_uv,
                bheight_uv);
        } else {
            svt_full_distortion_kernel_cbf_zero32_bits(
                &(((int32_t *)coeff->buffer_cb)[coeff_chroma_origin_index]),
                bwidth_uv,
                cb_distortion,
                bwidth_uv,
                bheight_uv);
        }
    }
    if (component_type == COMPONENT_CHROMA_CR || component_type == COMPONENT_CHROMA ||
        component_type == COMPONENT_ALL) {
        cr_distortion[0] = 0;
        cr_distortion[1] = 0;
        // CR
        if (cr_count_non_zero_coeffs) {
            svt_full_distortion_kernel32_bits(
                &(((int32_t *)coeff->buffer_cr)[coeff_chroma_origin_index]),
                bwidth_uv,
                &(((int32_t *)recon_coeff->buffer_cr)[recon_coeff_chroma_origin_index]),
                bwidth_uv,
                cr_distortion,
                bwidth_uv,
                bheight_uv);
        } else {
            svt_full_distortion_kernel_cbf_zero32_bits(
                &(((int32_t *)coeff->buffer_cr)[coeff_chroma_origin_index]),
                bwidth_uv,
                cr_distortion,
                bwidth_uv,
                bheight_uv);
        }
    }

    return return_error;
}
void un_pack2d(uint16_t *in16_bit_buffer, uint32_t in_stride, uint8_t *out8_bit_buffer,
               uint32_t out8_stride, uint8_t *outn_bit_buffer, uint32_t outn_stride, uint32_t width,
               uint32_t height) {
    if (((width & 3) == 0) && ((height & 1) == 0)) {
        svt_un_pack2d_16_bit_src_mul4(in16_bit_buffer,
                                      in_stride,
                                      out8_bit_buffer,
                                      outn_bit_buffer,
                                      out8_stride,
                                      outn_stride,
                                      width,
                                      height);
    } else {
        svt_enc_msb_un_pack2_d(in16_bit_buffer,
                               in_stride,
                               out8_bit_buffer,
                               outn_bit_buffer,
                               out8_stride,
                               outn_stride,
                               width,
                               height);
    }
}

void pack2d_src(uint8_t *in8_bit_buffer, uint32_t in8_stride, uint8_t *inn_bit_buffer,
                uint32_t inn_stride, uint16_t *out16_bit_buffer, uint32_t out_stride,
                uint32_t width, uint32_t height) {
    if (((width & 3) == 0) && ((height & 1) == 0)) {
        svt_pack2d_16_bit_src_mul4(in8_bit_buffer,
                                   in8_stride,
                                   inn_bit_buffer,
                                   out16_bit_buffer,
                                   inn_stride,
                                   out_stride,
                                   width,
                                   height);
    } else {
        svt_enc_msb_pack2_d(in8_bit_buffer,
                            in8_stride,
                            inn_bit_buffer,
                            out16_bit_buffer,
                            inn_stride,
                            out_stride,
                            width,
                            height);
    }
}

void compressed_pack_sb(uint8_t *in8_bit_buffer, uint32_t in8_stride, uint8_t *inn_bit_buffer,
                        uint32_t inn_stride, uint16_t *out16_bit_buffer, uint32_t out_stride,
                        uint32_t width, uint32_t height) {
    if (width == 64 || width == 32) {
        svt_compressed_packmsb(in8_bit_buffer,
                           in8_stride,
                           inn_bit_buffer,
                           out16_bit_buffer,
                           inn_stride,
                           out_stride,
                           width,
                           height);
    } else {
        svt_compressed_packmsb_c(in8_bit_buffer,
                             in8_stride,
                             inn_bit_buffer,
                             out16_bit_buffer,
                             inn_stride,
                             out_stride,
                             width,
                             height);
    }
}
// Copies the source image into the destination image and updates the
// destination's UMV borders.
// Note: The frames are assumed to be identical in size.
void svt_aom_yv12_copy_y_c(const Yv12BufferConfig *src_ybc, Yv12BufferConfig *dst_ybc) {
    int32_t        row;
    const uint8_t *src = src_ybc->y_buffer;
    uint8_t *      dst = dst_ybc->y_buffer;

    if (src_ybc->flags & YV12_FLAG_HIGHBITDEPTH) {
        const uint16_t *src16 = CONVERT_TO_SHORTPTR(src);
        uint16_t *      dst16 = CONVERT_TO_SHORTPTR(dst);
        for (row = 0; row < src_ybc->y_height; ++row) {
            svt_memcpy(dst16, src16, src_ybc->y_width * sizeof(uint16_t));
            src16 += src_ybc->y_stride;
            dst16 += dst_ybc->y_stride;
        }
        return;
    }

    for (row = 0; row < src_ybc->y_height; ++row) {
        svt_memcpy(dst, src, src_ybc->y_width);
        src += src_ybc->y_stride;
        dst += dst_ybc->y_stride;
    }
}

void svt_aom_yv12_copy_u_c(const Yv12BufferConfig *src_bc, Yv12BufferConfig *dst_bc) {
    int32_t        row;
    const uint8_t *src = src_bc->u_buffer;
    uint8_t *      dst = dst_bc->u_buffer;

    if (src_bc->flags & YV12_FLAG_HIGHBITDEPTH) {
        const uint16_t *src16 = CONVERT_TO_SHORTPTR(src);
        uint16_t *      dst16 = CONVERT_TO_SHORTPTR(dst);
        for (row = 0; row < src_bc->uv_height; ++row) {
            svt_memcpy(dst16, src16, src_bc->uv_width * sizeof(uint16_t));
            src16 += src_bc->uv_stride;
            dst16 += dst_bc->uv_stride;
        }
        return;
    }

    for (row = 0; row < src_bc->uv_height; ++row) {
        svt_memcpy(dst, src, src_bc->uv_width);
        src += src_bc->uv_stride;
        dst += dst_bc->uv_stride;
    }
}

void svt_aom_yv12_copy_v_c(const Yv12BufferConfig *src_bc, Yv12BufferConfig *dst_bc) {
    int32_t        row;
    const uint8_t *src = src_bc->v_buffer;
    uint8_t *      dst = dst_bc->v_buffer;

    if (src_bc->flags & YV12_FLAG_HIGHBITDEPTH) {
        const uint16_t *src16 = CONVERT_TO_SHORTPTR(src);
        uint16_t *      dst16 = CONVERT_TO_SHORTPTR(dst);
        for (row = 0; row < src_bc->uv_height; ++row) {
            svt_memcpy(dst16, src16, src_bc->uv_width * sizeof(uint16_t));
            src16 += src_bc->uv_stride;
            dst16 += dst_bc->uv_stride;
        }
        return;
    }

    for (row = 0; row < src_bc->uv_height; ++row) {
        svt_memcpy(dst, src, src_bc->uv_width);
        src += src_bc->uv_stride;
        dst += dst_bc->uv_stride;
    }
}
