/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EBAVCSTYLEMCP_H
#define EBAVCSTYLEMCP_H
#include "EbDefinitions.h"
#include "EbPictureBufferDesc.h"

#ifdef __cplusplus
extern "C" {
#endif
#if !REMOVE_UNUSED_CODE
void avc_style_copy_c(EbByte refPic, uint32_t srcStride, EbByte dst, uint32_t dstStride,
                      uint32_t puWidth, uint32_t puHeight, EbByte tempBuf, uint32_t fracPos);
void avc_style_luma_interpolation_filter_horizontal_c(EbByte ref_pic, uint32_t src_stride,
                                                      EbByte dst, uint32_t dst_stride,
                                                      uint32_t pu_width, uint32_t pu_height,
                                                      EbByte temp_buf, uint32_t frac_pos);
void avc_style_luma_interpolation_filter_vertical_c(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                    uint32_t dst_stride, uint32_t pu_width,
                                                    uint32_t pu_height, EbByte temp_buf,
                                                    uint32_t frac_pos);
void avc_style_luma_interpolation_filter_pose_c(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                uint32_t dst_stride, uint32_t pu_width,
                                                uint32_t pu_height, EbByte temp_buf,
                                                uint32_t frac_pos);
void avc_style_luma_interpolation_filter_posf_c(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                uint32_t dst_stride, uint32_t pu_width,
                                                uint32_t pu_height, EbByte temp_buf,
                                                uint32_t frac_pos);
void avc_style_luma_interpolation_filter_posg_c(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                uint32_t dst_stride, uint32_t pu_width,
                                                uint32_t pu_height, EbByte temp_buf,
                                                uint32_t frac_pos);
void avc_style_luma_interpolation_filter_posi_c(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                uint32_t dst_stride, uint32_t pu_width,
                                                uint32_t pu_height, EbByte temp_buf,
                                                uint32_t frac_pos);
void avc_style_luma_interpolation_filter_posj_c(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                uint32_t dst_stride, uint32_t pu_width,
                                                uint32_t pu_height, EbByte temp_buf,
                                                uint32_t frac_pos);
void avc_style_luma_interpolation_filter_posk_c(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                uint32_t dst_stride, uint32_t pu_width,
                                                uint32_t pu_height, EbByte temp_buf,
                                                uint32_t frac_pos);
void avc_style_luma_interpolation_filter_posp_c(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                uint32_t dst_stride, uint32_t pu_width,
                                                uint32_t pu_height, EbByte temp_buf,
                                                uint32_t frac_pos);
void avc_style_luma_interpolation_filter_posq_c(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                uint32_t dst_stride, uint32_t pu_width,
                                                uint32_t pu_height, EbByte temp_buf,
                                                uint32_t frac_pos);
void avc_style_luma_interpolation_filter_posr_c(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                uint32_t dst_stride, uint32_t pu_width,
                                                uint32_t pu_height, EbByte temp_buf,
                                                uint32_t frac_pos);

void avc_style_luma_interpolation_filter_helper_c(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                  uint32_t dst_stride, uint32_t pu_width,
                                                  uint32_t pu_height, EbByte temp_buf,
                                                  uint32_t frac_pos, uint8_t fractional_position);
#endif
#ifdef __cplusplus
}
#endif
#endif //EBAVCSTYLEMCP_H
