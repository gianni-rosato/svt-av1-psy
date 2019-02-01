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
    void ChromaInterpolationCopy_SSSE3(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst, uint32_t frac_pos_x, uint32_t frac_pos_y);
    void ChromaInterpolationFilterOneDHorizontal_SSSE3(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst, uint32_t frac_pos_x, uint32_t frac_pos_y);
    void ChromaInterpolationFilterOneDOutRawHorizontal_SSSE3(EbByte ref_pic, uint32_t src_stride, int16_t *dst, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst, uint32_t frac_pos_x, uint32_t frac_pos_y);
    void ChromaInterpolationFilterOneDVertical_SSSE3(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst, uint32_t frac_pos_x, uint32_t frac_pos_y);
    void ChromaInterpolationFilterOneDOutRawVertical_SSSE3(EbByte ref_pic, uint32_t src_stride, int16_t *dst, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst, uint32_t frac_pos_x, uint32_t frac_pos_y);
    void ChromaInterpolationFilterTwoD_SSSE3(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst, uint32_t frac_pos_x, uint32_t frac_pos_y);
    void ChromaInterpolationFilterTwoDOutRaw_SSSE3(EbByte ref_pic, uint32_t src_stride, int16_t *dst, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst, uint32_t frac_pos_x, uint32_t frac_pos_y);
    void ChromaInterpolationCopyOutRaw_SSSE3(EbByte ref_pic, uint32_t src_stride, int16_t *dst, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst, uint32_t frac_pos_x, uint32_t frac_pos_y);
    void LumaInterpolationCopy_SSSE3(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst);
    void LumaInterpolationFilterPosa_SSSE3(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst);
    void LumaInterpolationFilterPosb_SSSE3(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst);
    void LumaInterpolationFilterPosc_SSSE3(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst);
    void LumaInterpolationFilterPosd_SSSE3(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst);
    void LumaInterpolationFilterPose_SSSE3(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst);
    void LumaInterpolationFilterPosf_SSSE3(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst);
    void LumaInterpolationFilterPosg_SSSE3(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst);
    void LumaInterpolationFilterPosh_SSSE3(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst);
    void LumaInterpolationFilterPosi_SSSE3(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst);
    void LumaInterpolationFilterPosj_SSSE3(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst);
    void LumaInterpolationFilterPosk_SSSE3(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst);
    void LumaInterpolationFilterPosn_SSSE3(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst);
    void LumaInterpolationFilterPosp_SSSE3(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst);
    void LumaInterpolationFilterPosq_SSSE3(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst);
    void LumaInterpolationFilterPosr_SSSE3(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst);

    void LumaInterpolationCopyOutRaw_SSSE3(EbByte ref_pic, uint32_t src_stride, int16_t *dst, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst);
    void LumaInterpolationFilterPosaOutRaw_SSSE3(EbByte ref_pic, uint32_t src_stride, int16_t* dst, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst);
    void LumaInterpolationFilterPosbOutRaw_SSSE3(EbByte ref_pic, uint32_t src_stride, int16_t* dst, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst);
    void LumaInterpolationFilterPoscOutRaw_SSSE3(EbByte ref_pic, uint32_t src_stride, int16_t* dst, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst);
    void LumaInterpolationFilterPosdOutRaw_SSSE3(EbByte ref_pic, uint32_t src_stride, int16_t* dst, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst);
    void LumaInterpolationFilterPoseOutRaw_SSSE3(EbByte ref_pic, uint32_t src_stride, int16_t* dst, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst);
    void LumaInterpolationFilterPosfOutRaw_SSSE3(EbByte ref_pic, uint32_t src_stride, int16_t* dst, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst);
    void LumaInterpolationFilterPosgOutRaw_SSSE3(EbByte ref_pic, uint32_t src_stride, int16_t* dst, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst);
    void LumaInterpolationFilterPoshOutRaw_SSSE3(EbByte ref_pic, uint32_t src_stride, int16_t* dst, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst);
    void LumaInterpolationFilterPosiOutRaw_SSSE3(EbByte ref_pic, uint32_t src_stride, int16_t* dst, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst);
    void LumaInterpolationFilterPosjOutRaw_SSSE3(EbByte ref_pic, uint32_t src_stride, int16_t* dst, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst);
    void LumaInterpolationFilterPoskOutRaw_SSSE3(EbByte ref_pic, uint32_t src_stride, int16_t* dst, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst);
    void LumaInterpolationFilterPosnOutRaw_SSSE3(EbByte ref_pic, uint32_t src_stride, int16_t* dst, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst);
    void LumaInterpolationFilterPospOutRaw_SSSE3(EbByte ref_pic, uint32_t src_stride, int16_t* dst, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst);
    void LumaInterpolationFilterPosqOutRaw_SSSE3(EbByte ref_pic, uint32_t src_stride, int16_t* dst, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst);
    void LumaInterpolationFilterPosrOutRaw_SSSE3(EbByte ref_pic, uint32_t src_stride, int16_t* dst, uint32_t pu_width, uint32_t pu_height, int16_t *first_pass_if_dst);

#ifdef __cplusplus
}
#endif
#endif //EBMCP_SSSE3_H