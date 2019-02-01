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

    void AvcStyleLumaInterpolationFilterPose_SSSE3(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, EbByte tempBuf, EbBool skip, uint32_t fracPos);
    void AvcStyleLumaInterpolationFilterPosf_SSSE3(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, EbByte tempBuf, EbBool skip, uint32_t fracPos);
    void AvcStyleLumaInterpolationFilterPosg_SSSE3(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, EbByte tempBuf, EbBool skip, uint32_t fracPos);
    void AvcStyleLumaInterpolationFilterPosi_SSSE3(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, EbByte tempBuf, EbBool skip, uint32_t fracPos);
    void AvcStyleLumaInterpolationFilterPosj_SSSE3(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, EbByte tempBuf, EbBool skip, uint32_t fracPos);
    void AvcStyleLumaInterpolationFilterPosk_SSSE3(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, EbByte tempBuf, EbBool skip, uint32_t fracPos);
    void AvcStyleLumaInterpolationFilterPosp_SSSE3(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, EbByte tempBuf, EbBool skip, uint32_t fracPos);
    void AvcStyleLumaInterpolationFilterPosq_SSSE3(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, EbByte tempBuf, EbBool skip, uint32_t fracPos);
    void AvcStyleLumaInterpolationFilterPosr_SSSE3(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, EbByte tempBuf, EbBool skip, uint32_t fracPos);
    void AvcStyleLumaInterpolationFilterHorizontal_SSSE3_INTRIN(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, EbByte tempBuf, EbBool skip, uint32_t fracPos);
    void AvcStyleLumaInterpolationFilterVertical_SSSE3_INTRIN(EbByte ref_pic, uint32_t src_stride, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, EbByte tempBuf, EbBool skip, uint32_t fracPos);
#ifdef __cplusplus
}
#endif
#endif //EBAVCSTYLEMCP_H
