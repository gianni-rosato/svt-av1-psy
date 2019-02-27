/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/
#ifndef EbTransforms_AVX2_h
#define EbTransforms_AVX2_h

#include "EbDefinitions.h"
#ifdef __cplusplus
extern "C" {
#endif

    void quantize_inv_quantize8x8_avx2_intrin(
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

    void quantize_inv_quantize_nxn_avx2_intrin(
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

    void low_precision_transform16x16_avx2_intrin(
        int16_t *src, 
        uint32_t src_stride, 
        int16_t *dst, 
        uint32_t dst_stride, 
        int16_t *intermediate, 
        uint32_t addshift);

    void low_precision_transform32x32_avx2_intrin(
        int16_t *src, 
        uint32_t src_stride, 
        int16_t *dst, 
        uint32_t dst_stride, 
        int16_t *intermediate, 
        uint32_t addshift);

    void pfreq_transform32x32_avx2_intrin(
        int16_t       *src,
        const uint32_t src_stride,
        int16_t       *dst,
        const uint32_t dst_stride,
        int16_t       *intermediate,
        uint32_t       addshift);

    void pfreq_n4_transform32x32_avx2_intrin(
        int16_t       *src,
        const uint32_t src_stride,
        int16_t       *dst,
        const uint32_t dst_stride,
        int16_t       *intermediate,
        uint32_t       addshift);

    void mat_mult4x4_avx2_intrin(
        int16_t*        coeff,
        const uint32_t  coeff_stride,
        const uint16_t *masking_matrix,
        const uint32_t  masking_matrix_stride,  //Matrix size
        const uint32_t  compute_size,           //Computation area size
        const int32_t   offset,                 //(PMP_MAX >> 1)
        const int32_t   shift_num,              //PMP_PRECISION
        uint32_t*       nonzerocoeff);

    void mat_mult4x4_out_buff_avx2_intrin(
        int16_t*        coeff,
        const uint32_t  coeff_stride,
        int16_t*        coeff_out,
        const uint32_t  coeff_out_stride,
        const uint16_t *masking_matrix,
        const uint32_t  masking_matrix_stride,
        const uint32_t  compute_size,
        const int32_t   offset,
        const int32_t   shift_num,
        uint32_t*       nonzerocoeff);

    void mat_mult8x8_avx2_intrin(
        int16_t*        coeff,
        const uint32_t  coeff_stride,
        const uint16_t *masking_matrix,
        const uint32_t  masking_matrix_stride,  //Matrix size
        const uint32_t  compute_size,           //Computation area size
        const int32_t   offset,                 //(PMP_MAX >> 1)
        const int32_t   shift_num,              //PMP_PRECISION
        uint32_t*       nonzerocoeff);

    void mat_mult_nxn_avx2_intrin(
        int16_t*        coeff,
        const uint32_t  coeff_stride,
        const uint16_t *masking_matrix,
        const uint32_t  masking_matrix_stride,  //Matrix size
        const uint32_t  compute_size,           //Computation area size
        const int32_t   offset,                 //(PMP_MAX >> 1)
        const int32_t   shift_num,              //PMP_PRECISION
        uint32_t*       nonzerocoeff);

#ifdef __cplusplus
}
#endif
#endif // EbTransforms_AVX2_h


