/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EBAVCSTYLEMCP_H
#define EBAVCSTYLEMCP_H
#include "EbPictureOperators.h"
#include "EbMcp.h"
#include "EbDefinitions.h"
#include "EbPictureBufferDesc.h"
#include "EbPictureControlSet.h"

#ifdef __cplusplus
extern "C" {
#endif

    void estimate_bi_pred_interpolation_avc_luma(
        EbPictureBufferDesc *ref_pic_list0,
        EbPictureBufferDesc *ref_pic_list1,
        uint32_t                 ref_list0_pos_x,
        uint32_t                 ref_list0_pos_y,
        uint32_t                 ref_list1_pos_x,
        uint32_t                 ref_list1_pos_y,
        uint32_t                 pu_width,
        uint32_t                 pu_height,
        EbPictureBufferDesc *bi_dst,
        uint32_t                 dst_luma_index,
        uint32_t                 dst_chroma_index,
        uint32_t                 component_mask,
        EbByte                ref_list0_temp_dst,
        EbByte                ref_list1_temp_dst,
        EbByte                first_pass_if_temp_dst,
        EbBool                sub_sample_pred_flag);

    void estimate_uni_pred_interpolation_avc_luma(
        EbPictureBufferDesc *ref_pic,
        uint32_t                 pos_x,
        uint32_t                 pos_y,
        uint32_t                 pu_width,
        uint32_t                 pu_height,
        EbPictureBufferDesc *dst,
        uint32_t                 dst_luma_index,
        uint32_t                 dst_chroma_index,
        uint32_t                 component_mask,
        EbByte                temp_buf,
        EbBool                sub_sample_pred_flag);

    void estimate_bi_pred_interpolation_unpacked_avc_style(
        EbPictureBufferDesc *ref_pic_list0,
        EbPictureBufferDesc *ref_pic_list1,
        uint32_t                 ref_list0_pos_x,
        uint32_t                 ref_list0_pos_y,
        uint32_t                 ref_list1_pos_x,
        uint32_t                 ref_list1_pos_y,
        uint32_t                 pu_width,
        uint32_t                 pu_height,
        EbPictureBufferDesc *bi_dst,
        uint32_t                 dst_luma_index,
        uint32_t                 dst_chroma_index,
        uint32_t                 component_mask,
        EbByte                ref_list0_temp_dst,
        EbByte                ref_list1_temp_dst,
        EbByte                first_pass_if_temp_dst,
        EbBool                sub_sample_pred_flag);

    void estimate_uni_pred_interpolation_unpacked_avc_style(
        EbPictureBufferDesc *ref_pic,
        uint32_t                 pos_x,
        uint32_t                 pos_y,
        uint32_t                 pu_width,
        uint32_t                 pu_height,
        EbPictureBufferDesc *dst,
        uint32_t                 dst_luma_index,
        uint32_t                 dst_chroma_index,
        uint32_t                 component_mask,
        EbByte                temp_buf,
        EbBool                sub_sample_pred_flag);

    void estimate_uni_pred_interpolation_avc_chroma_ref10_bit(
        EbPictureBufferDesc *ref_frame_pic_list0,
        uint32_t                 pos_x,
        uint32_t                 pos_y,
        uint32_t                 pu_width,
        uint32_t                 pu_height,
        EbPictureBufferDesc *dst,
        uint32_t                 dst_luma_index,
        uint32_t                 dst_chroma_index,          //input parameter, please refer to the detailed explanation above.
        uint32_t                 component_mask,
        EbByte                temp_buf,
        EbBool                sub_pred);

    void estimate_bi_pred_interpolation_avc_chroma_ref10_bit(
        EbPictureBufferDesc *ref_frame_pic_list0,
        EbPictureBufferDesc *ref_frame_pic_list1,
        uint32_t                 ref_list0_pos_x,
        uint32_t                 ref_list0_pos_y,
        uint32_t                 ref_list1_pos_x,
        uint32_t                 ref_list1_pos_y,
        uint32_t                 pu_width,
        uint32_t                 pu_height,
        EbPictureBufferDesc *bi_dst,
        uint32_t                 dst_luma_index,
        uint32_t                 dst_chroma_index,
        uint32_t                 component_mask,
        EbByte                ref_list0_temp_dst,
        EbByte                ref_list1_temp_dst,
        EbByte                first_pass_if_temp_dst,
        EbBool                sub_pred);

    void estimate_uni_pred_interpolation_avc_lumaRef10Bit(
        EbPictureBufferDesc *ref_frame_pic_list0,
        uint32_t                 pos_x,
        uint32_t                 pos_y,
        uint32_t                 pu_width,
        uint32_t                 pu_height,
        EbPictureBufferDesc *dst,
        uint32_t                 dst_luma_index,
        uint32_t                 dst_chroma_index,          //input parameter, please refer to the detailed explanation above.
        uint32_t                 component_mask,
        EbByte                temp_buf,
        EbBool                sub_pred,
        EbBool                sub_pred_chroma);

    void estimate_bi_pred_interpolation_avc_luma_ref10_bit(
        EbPictureBufferDesc *ref_frame_pic_list0,
        EbPictureBufferDesc *ref_frame_pic_list1,
        uint32_t                 ref_list0_pos_x,
        uint32_t                 ref_list0_pos_y,
        uint32_t                 ref_list1_pos_x,
        uint32_t                 ref_list1_pos_y,
        uint32_t                 pu_width,
        uint32_t                 pu_height,
        EbPictureBufferDesc *bi_dst,
        uint32_t                 dst_luma_index,
        uint32_t                 dst_chroma_index,
        uint32_t                 component_mask,
        EbByte                ref_list0_temp_dst,
        EbByte                ref_list1_temp_dst,
        EbByte                first_pass_if_temp_dst,
        EbBool                sub_pred,
        EbBool                sub_pred_chroma);

    void uni_pred_i_free_ref8_bit(
        EbPictureBufferDesc *ref_pic,
        uint32_t                 pos_x,
        uint32_t                 pos_y,
        uint32_t                 pu_width,
        uint32_t                 pu_height,
        EbPictureBufferDesc *dst,
        uint32_t                 dst_luma_index,
        uint32_t                 dst_chroma_index,          //input parameter, please refer to the detailed explanation above.
        uint32_t                 component_mask,
        EbByte                temp_buf,
        EbBool                sub_sample_pred_flag,
        EbBool                subSamplePredFlagChroma);

    void bi_pred_i_free_ref8_bit(
        EbPictureBufferDesc *ref_pic_list0,
        EbPictureBufferDesc *ref_pic_list1,
        uint32_t                 ref_list0_pos_x,
        uint32_t                 ref_list0_pos_y,
        uint32_t                 ref_list1_pos_x,
        uint32_t                 ref_list1_pos_y,
        uint32_t                 pu_width,
        uint32_t                 pu_height,
        EbPictureBufferDesc *bi_dst,
        uint32_t                 dst_luma_index,
        uint32_t                 dst_chroma_index,
        uint32_t                 component_mask,
        EbByte                ref_list0_temp_dst,
        EbByte                ref_list1_temp_dst,
        EbByte                first_pass_if_temp_dst,
        EbBool                sub_sample_pred_flag,
        EbBool                subSamplePredFlagChroma);

    void avc_style_copy_c(EbByte refPic, uint32_t srcStride, EbByte dst,
                    uint32_t dstStride, uint32_t puWidth, uint32_t puHeight,
                    EbByte tempBuf, uint32_t fracPos);
    void avc_style_luma_interpolation_filter_horizontal_c(
        EbByte ref_pic, uint32_t src_stride,
        EbByte dst, uint32_t dst_stride,
        uint32_t pu_width, uint32_t pu_height,
        EbByte temp_buf, uint32_t frac_pos);
    void avc_style_luma_interpolation_filter_vertical_c(
        EbByte ref_pic, uint32_t src_stride,
        EbByte dst, uint32_t dst_stride,
        uint32_t pu_width, uint32_t pu_height,
        EbByte temp_buf, uint32_t frac_pos);
    void avc_style_luma_interpolation_filter_pose_c(EbByte ref_pic, uint32_t src_stride,
        EbByte dst, uint32_t dst_stride,
        uint32_t pu_width, uint32_t pu_height,
        EbByte temp_buf, uint32_t frac_pos);
    void avc_style_luma_interpolation_filter_posf_c(EbByte ref_pic, uint32_t src_stride,
        EbByte dst, uint32_t dst_stride,
        uint32_t pu_width, uint32_t pu_height,
        EbByte temp_buf, uint32_t frac_pos);
    void avc_style_luma_interpolation_filter_posg_c(EbByte ref_pic, uint32_t src_stride,
        EbByte dst, uint32_t dst_stride,
        uint32_t pu_width, uint32_t pu_height,
        EbByte temp_buf, uint32_t frac_pos);
    void avc_style_luma_interpolation_filter_posi_c(EbByte ref_pic, uint32_t src_stride,
        EbByte dst, uint32_t dst_stride,
        uint32_t pu_width, uint32_t pu_height,
        EbByte temp_buf, uint32_t frac_pos);
    void avc_style_luma_interpolation_filter_posj_c(EbByte ref_pic, uint32_t src_stride,
        EbByte dst, uint32_t dst_stride,
        uint32_t pu_width, uint32_t pu_height,
        EbByte temp_buf, uint32_t frac_pos);
    void avc_style_luma_interpolation_filter_posk_c(EbByte ref_pic, uint32_t src_stride,
        EbByte dst, uint32_t dst_stride,
        uint32_t pu_width, uint32_t pu_height,
        EbByte temp_buf, uint32_t frac_pos);
    void avc_style_luma_interpolation_filter_posp_c(EbByte ref_pic, uint32_t src_stride,
        EbByte dst, uint32_t dst_stride,
        uint32_t pu_width, uint32_t pu_height,
        EbByte temp_buf, uint32_t frac_pos);
    void avc_style_luma_interpolation_filter_posq_c(EbByte ref_pic, uint32_t src_stride,
        EbByte dst, uint32_t dst_stride,
        uint32_t pu_width, uint32_t pu_height,
        EbByte temp_buf, uint32_t frac_pos);
    void avc_style_luma_interpolation_filter_posr_c(EbByte ref_pic, uint32_t src_stride,
        EbByte dst, uint32_t dst_stride,
        uint32_t pu_width, uint32_t pu_height,
        EbByte temp_buf, uint32_t frac_pos);

    void avc_style_luma_interpolation_filter_helper_c(
        EbByte ref_pic,
        uint32_t src_stride,
        EbByte dst,
        uint32_t dst_stride,
        uint32_t pu_width,
        uint32_t pu_height,
        EbByte temp_buf,
        EbBool skip,
        uint32_t frac_pos,
        uint8_t fractional_position);

#ifdef __cplusplus
}
#endif
#endif //EBAVCSTYLEMCP_H
