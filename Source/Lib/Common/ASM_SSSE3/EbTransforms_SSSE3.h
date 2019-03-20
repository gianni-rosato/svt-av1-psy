/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/
#ifndef EbTransforms_SSSE3_h
#define EbTransforms_SSSE3_h

#include "EbDefinitions.h"
#ifdef __cplusplus
extern "C" {
#endif

    void quantize_inv_quantize_nx_n_sse3(
        int16_t        *coeff,
        const uint32_t  coeff_stride,
        int16_t        *quant_coeff,
        int16_t        *recon_coeff,
        const uint32_t  q_func,
        const uint32_t  q_offset,
        const int32_t   shifted_q_bits,
        const int32_t   shifted_f_func,
        const int32_t   iq_offset,
        const int32_t   shift_num,
        const uint32_t  area_size,
        uint32_t       *nonzerocoeff);

    void quantize_inv_quantize8x8_sse3(
        int16_t        *coeff,
        const uint32_t  coeff_stride,
        int16_t        *quant_coeff,
        int16_t        *recon_coeff,
        const uint32_t  q_func,
        const uint32_t  q_offset,
        const int32_t   shifted_q_bits,
        const int32_t   shifted_f_func,
        const int32_t   iq_offset,
        const int32_t   shift_num,
        const uint32_t  area_size,
        uint32_t       *nonzerocoeff);

    void quantize_inv_quantize4x4_sse3(
        int16_t        *coeff,
        const uint32_t  coeff_stride,
        int16_t        *quant_coeff,
        int16_t        *recon_coeff,
        const uint32_t  q_func,
        const uint32_t  q_offset,
        const int32_t   shifted_q_bits,
        const int32_t   shifted_f_func,
        const int32_t   iq_offset,
        const int32_t   shift_num,
        const uint32_t  area_size,
        uint32_t       *nonzerocoeff);


    void low_precision_transform16x16_ssse3(
        int16_t *src, 
        uint32_t src_stride, 
        int16_t *dst, 
        uint32_t dst_stride, 
        int16_t *intermediate, 
        uint32_t addshift);

    void low_precision_transform32x32_ssse3(
        int16_t *src, 
        uint32_t src_stride, 
        int16_t *dst, 
        uint32_t dst_stride, 
        int16_t *intermediate, 
        uint32_t addshift);

    void p_finv_transform32x32_ssse3(
        int16_t *src, 
        uint32_t src_stride, 
        int16_t *dst, 
        uint32_t dst_stride, 
        int16_t *intermediate, 
        uint32_t addshift);

    void p_finv_transform16x16_ssse3(
        int16_t *src, 
        uint32_t src_stride, 
        int16_t *dst, 
        uint32_t dst_stride, 
        int16_t *intermediate, 
        uint32_t addshift);

#ifdef __cplusplus
}
#endif
#endif // EbTransforms_SSSE3_h

