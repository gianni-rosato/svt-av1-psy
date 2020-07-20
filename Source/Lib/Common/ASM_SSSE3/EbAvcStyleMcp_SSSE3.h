/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EBAVCSTYLEMCP_SSSE3_H
#define EBAVCSTYLEMCP_SSSE3_H

#include "EbDefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif
#if !REMOVE_UNUSED_CODE
void avc_style_luma_interpolation_filter_pose_ssse3(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                    uint32_t dst_stride, uint32_t pu_width,
                                                    uint32_t pu_height, EbByte temp_buf,
                                                    uint32_t frac_pos);

void avc_style_luma_interpolation_filter_posf_ssse3(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                    uint32_t dst_stride, uint32_t pu_width,
                                                    uint32_t pu_height, EbByte temp_buf,
                                                    uint32_t frac_pos);

void avc_style_luma_interpolation_filter_posg_ssse3(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                    uint32_t dst_stride, uint32_t pu_width,
                                                    uint32_t pu_height, EbByte temp_buf,
                                                    uint32_t frac_pos);

void avc_style_luma_interpolation_filter_posi_ssse3(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                    uint32_t dst_stride, uint32_t pu_width,
                                                    uint32_t pu_height, EbByte temp_buf,
                                                    uint32_t frac_pos);

void avc_style_luma_interpolation_filter_posj_ssse3(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                    uint32_t dst_stride, uint32_t pu_width,
                                                    uint32_t pu_height, EbByte temp_buf,
                                                    uint32_t frac_pos);

void avc_style_luma_interpolation_filter_posk_ssse3(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                    uint32_t dst_stride, uint32_t pu_width,
                                                    uint32_t pu_height, EbByte temp_buf,
                                                    uint32_t frac_pos);

void avc_style_luma_interpolation_filter_posp_ssse3(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                    uint32_t dst_stride, uint32_t pu_width,
                                                    uint32_t pu_height, EbByte temp_buf,
                                                    uint32_t frac_pos);

void avc_style_luma_interpolation_filter_posq_ssse3(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                    uint32_t dst_stride, uint32_t pu_width,
                                                    uint32_t pu_height, EbByte temp_buf,
                                                    uint32_t frac_pos);

void avc_style_luma_interpolation_filter_posr_ssse3(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                                                    uint32_t dst_stride, uint32_t pu_width,
                                                    uint32_t pu_height, EbByte temp_buf,
                                                    uint32_t frac_pos);

void avc_style_luma_interpolation_filter_horizontal_ssse3_intrin(
    EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width,
    uint32_t pu_height, EbByte temp_buf, uint32_t frac_pos);

void avc_style_luma_interpolation_filter_vertical_ssse3_intrin(EbByte ref_pic, uint32_t src_stride,
                                                               EbByte dst, uint32_t dst_stride,
                                                               uint32_t pu_width,
                                                               uint32_t pu_height, EbByte temp_buf,
                                                               uint32_t frac_pos);
#endif
#ifdef __cplusplus
}
#endif
#endif //EBAVCSTYLEMCP_H
