/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EBMCP_SSSE3_H
#define EBMCP_SSSE3_H

#include "EbDefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif

    // SSSE3 functions
    void chroma_interpolation_copy_ssse3(
        EbByte   ref_pic, 
        uint32_t src_stride, 
        EbByte   dst, 
        uint32_t dst_stride, 
        uint32_t pu_width, 
        uint32_t pu_height, 
        int16_t *first_pass_if_dst, 
        uint32_t frac_pos_x, 
        uint32_t frac_pos_y);

    void chroma_interpolation_filter_one_d_horizontal_ssse3(
        EbByte   ref_pic, 
        uint32_t src_stride, 
        EbByte   dst, 
        uint32_t dst_stride, 
        uint32_t pu_width, 
        uint32_t pu_height, 
        int16_t *first_pass_if_dst, 
        uint32_t frac_pos_x, 
        uint32_t frac_pos_y);

    void chroma_interpolation_filter_one_d_out_raw_horizontal_ssse3(
        EbByte   ref_pic, 
        uint32_t src_stride, 
        int16_t *dst, 
        uint32_t pu_width, 
        uint32_t pu_height, 
        int16_t *first_pass_if_dst, 
        uint32_t frac_pos_x, 
        uint32_t frac_pos_y);

    void chroma_interpolation_filter_one_d_vertical_ssse3(
        EbByte   ref_pic, 
        uint32_t src_stride, 
        EbByte   dst, 
        uint32_t dst_stride, 
        uint32_t pu_width,
        uint32_t pu_height, 
        int16_t *first_pass_if_dst,
        uint32_t frac_pos_x, 
        uint32_t frac_pos_y);

    void chroma_interpolation_filter_one_d_out_raw_vertical_ssse3(
        EbByte   ref_pic, 
        uint32_t src_stride, 
        int16_t *dst, 
        uint32_t pu_width, 
        uint32_t pu_height, 
        int16_t *first_pass_if_dst, 
        uint32_t frac_pos_x, 
        uint32_t frac_pos_y);

    void chroma_interpolation_filter_two_d_ssse3(
        EbByte   ref_pic, 
        uint32_t src_stride, 
        EbByte   dst, 
        uint32_t dst_stride,
        uint32_t pu_width, 
        uint32_t pu_height, 
        int16_t *first_pass_if_dst, 
        uint32_t frac_pos_x, 
        uint32_t frac_pos_y);

    void chroma_interpolation_filter_two_d_out_raw_ssse3(
        EbByte   ref_pic, 
        uint32_t src_stride, 
        int16_t *dst,
        uint32_t pu_width, 
        uint32_t pu_height, 
        int16_t *first_pass_if_dst,
        uint32_t frac_pos_x, 
        uint32_t frac_pos_y);

    void chroma_interpolation_copy_out_raw_ssse3(
        EbByte   ref_pic, 
        uint32_t src_stride,
        int16_t *dst, 
        uint32_t pu_width, 
        uint32_t pu_height, 
        int16_t *first_pass_if_dst,
        uint32_t frac_pos_x, 
        uint32_t frac_pos_y);

    void luma_interpolation_copy_ssse3(
        EbByte   ref_pic, 
        uint32_t src_stride, 
        EbByte   dst, 
        uint32_t dst_stride, 
        uint32_t pu_width, 
        uint32_t pu_height, 
        int16_t *first_pass_if_dst);

    void luma_interpolation_filter_posa_ssse3(
        EbByte   ref_pic, 
        uint32_t src_stride, 
        EbByte   dst, 
        uint32_t dst_stride, 
        uint32_t pu_width, 
        uint32_t pu_height, 
        int16_t *first_pass_if_dst);

    void luma_interpolation_filter_posb_ssse3(
        EbByte   ref_pic,
        uint32_t src_stride,
        EbByte   dst,
        uint32_t dst_stride,
        uint32_t pu_width,
        uint32_t pu_height,
        int16_t *first_pass_if_dst);

    void luma_interpolation_filter_posc_ssse3(
        EbByte   ref_pic,
        uint32_t src_stride,
        EbByte   dst,
        uint32_t dst_stride,
        uint32_t pu_width,
        uint32_t pu_height,
        int16_t *first_pass_if_dst);

    void luma_interpolation_filter_posd_ssse3(
        EbByte   ref_pic,
        uint32_t src_stride,
        EbByte   dst,
        uint32_t dst_stride,
        uint32_t pu_width,
        uint32_t pu_height,
        int16_t *first_pass_if_dst);

    void luma_interpolation_filter_pose_ssse3(
        EbByte   ref_pic,
        uint32_t src_stride,
        EbByte   dst,
        uint32_t dst_stride,
        uint32_t pu_width,
        uint32_t pu_height,
        int16_t *first_pass_if_dst);

    void luma_interpolation_filter_posf_ssse3(
        EbByte   ref_pic,
        uint32_t src_stride,
        EbByte   dst,
        uint32_t dst_stride,
        uint32_t pu_width,
        uint32_t pu_height,
        int16_t *first_pass_if_dst);

    void luma_interpolation_filter_posg_ssse3(
        EbByte   ref_pic,
        uint32_t src_stride,
        EbByte   dst,
        uint32_t dst_stride,
        uint32_t pu_width,
        uint32_t pu_height,
        int16_t *first_pass_if_dst);

    void luma_interpolation_filter_posh_ssse3(
        EbByte   ref_pic,
        uint32_t src_stride,
        EbByte   dst,
        uint32_t dst_stride,
        uint32_t pu_width,
        uint32_t pu_height,
        int16_t *first_pass_if_dst);

    void luma_interpolation_filter_posi_ssse3(
        EbByte   ref_pic,
        uint32_t src_stride,
        EbByte   dst,
        uint32_t dst_stride,
        uint32_t pu_width,
        uint32_t pu_height,
        int16_t *first_pass_if_dst);

    void luma_interpolation_filter_posj_ssse3(
        EbByte   ref_pic,
        uint32_t src_stride,
        EbByte   dst,
        uint32_t dst_stride,
        uint32_t pu_width,
        uint32_t pu_height,
        int16_t *first_pass_if_dst);

    void luma_interpolation_filter_posk_ssse3(
        EbByte   ref_pic,
        uint32_t src_stride,
        EbByte   dst,
        uint32_t dst_stride,
        uint32_t pu_width,
        uint32_t pu_height,
        int16_t *first_pass_if_dst);

    void luma_interpolation_filter_posn_ssse3(
        EbByte   ref_pic,
        uint32_t src_stride,
        EbByte   dst,
        uint32_t dst_stride,
        uint32_t pu_width,
        uint32_t pu_height,
        int16_t *first_pass_if_dst);

    void luma_interpolation_filter_posp_ssse3(
        EbByte   ref_pic,
        uint32_t src_stride,
        EbByte   dst,
        uint32_t dst_stride,
        uint32_t pu_width,
        uint32_t pu_height,
        int16_t *first_pass_if_dst);

    void luma_interpolation_filter_posq_ssse3(
        EbByte   ref_pic,
        uint32_t src_stride,
        EbByte   dst,
        uint32_t dst_stride,
        uint32_t pu_width,
        uint32_t pu_height,
        int16_t *first_pass_if_dst);

    void luma_interpolation_filter_posr_ssse3(
        EbByte   ref_pic,
        uint32_t src_stride,
        EbByte   dst,
        uint32_t dst_stride,
        uint32_t pu_width,
        uint32_t pu_height,
        int16_t *first_pass_if_dst);

    void luma_interpolation_copy_out_raw_ssse3(
        EbByte   ref_pic, 
        uint32_t src_stride, 
        int16_t *dst, 
        uint32_t pu_width, 
        uint32_t pu_height, 
        int16_t *first_pass_if_dst);

    void luma_interpolation_filter_posa_out_raw_ssse3(
        EbByte   ref_pic, 
        uint32_t src_stride, 
        int16_t *dst, 
        uint32_t pu_width, 
        uint32_t pu_height, 
        int16_t *first_pass_if_dst);

    void luma_interpolation_filter_posb_out_raw_ssse3(
        EbByte   ref_pic,
        uint32_t src_stride,
        int16_t *dst,
        uint32_t pu_width,
        uint32_t pu_height,
        int16_t *first_pass_if_dst);

    void luma_interpolation_filter_posc_out_raw_ssse3(
        EbByte   ref_pic,
        uint32_t src_stride,
        int16_t *dst,
        uint32_t pu_width,
        uint32_t pu_height,
        int16_t *first_pass_if_dst);

    void luma_interpolation_filter_posd_out_raw_ssse3(
        EbByte   ref_pic,
        uint32_t src_stride,
        int16_t *dst,
        uint32_t pu_width,
        uint32_t pu_height,
        int16_t *first_pass_if_dst);

    void luma_interpolation_filter_pose_out_raw_ssse3(
        EbByte   ref_pic,
        uint32_t src_stride,
        int16_t *dst,
        uint32_t pu_width,
        uint32_t pu_height,
        int16_t *first_pass_if_dst);

    void luma_interpolation_filter_posf_out_raw_ssse3(
        EbByte   ref_pic,
        uint32_t src_stride,
        int16_t *dst,
        uint32_t pu_width,
        uint32_t pu_height,
        int16_t *first_pass_if_dst);

    void luma_interpolation_filter_posg_out_raw_ssse3(
        EbByte   ref_pic,
        uint32_t src_stride,
        int16_t *dst,
        uint32_t pu_width,
        uint32_t pu_height,
        int16_t *first_pass_if_dst);

    void luma_interpolation_filter_posh_out_raw_ssse3(
        EbByte   ref_pic,
        uint32_t src_stride,
        int16_t *dst,
        uint32_t pu_width,
        uint32_t pu_height,
        int16_t *first_pass_if_dst);

    void luma_interpolation_filter_posi_out_raw_ssse3(
        EbByte   ref_pic,
        uint32_t src_stride,
        int16_t *dst,
        uint32_t pu_width,
        uint32_t pu_height,
        int16_t *first_pass_if_dst);

    void luma_interpolation_filter_posj_out_raw_ssse3(
        EbByte   ref_pic,
        uint32_t src_stride,
        int16_t *dst,
        uint32_t pu_width,
        uint32_t pu_height,
        int16_t *first_pass_if_dst);

    void luma_interpolation_filter_posk_out_raw_ssse3(
        EbByte   ref_pic,
        uint32_t src_stride,
        int16_t *dst,
        uint32_t pu_width,
        uint32_t pu_height,
        int16_t *first_pass_if_dst);

    void luma_interpolation_filter_posn_out_raw_ssse3(
        EbByte   ref_pic,
        uint32_t src_stride,
        int16_t *dst,
        uint32_t pu_width,
        uint32_t pu_height,
        int16_t *first_pass_if_dst);

    void luma_interpolation_filter_posp_out_raw_ssse3(
        EbByte   ref_pic,
        uint32_t src_stride,
        int16_t *dst,
        uint32_t pu_width,
        uint32_t pu_height,
        int16_t *first_pass_if_dst);

    void luma_interpolation_filter_posq_out_raw_ssse3(
        EbByte   ref_pic,
        uint32_t src_stride,
        int16_t *dst,
        uint32_t pu_width,
        uint32_t pu_height,
        int16_t *first_pass_if_dst);

    void luma_interpolation_filter_posr_out_raw_ssse3(
        EbByte   ref_pic,
        uint32_t src_stride,
        int16_t *dst,
        uint32_t pu_width,
        uint32_t pu_height,
        int16_t *first_pass_if_dst);

#ifdef __cplusplus
}
#endif
#endif //EBMCP_SSSE3_H