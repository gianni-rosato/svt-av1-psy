/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPictureOperators_AVX2
#define EbPictureOperators_AVX2

#include "EbDefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif

    extern void eb_enc_msb_pack2d_avx2_intrin_al(
        uint8_t  *in8_bit_buffer,
        uint32_t  in8_stride,
        uint8_t  *inn_bit_buffer,
        uint16_t *out16_bit_buffer,
        uint32_t  inn_stride,
        uint32_t  out_stride,
        uint32_t  width,
        uint32_t  height);
                 
    extern void compressed_packmsb_avx2_intrin(
        uint8_t  *in8_bit_buffer,
        uint32_t  in8_stride,
        uint8_t  *inn_bit_buffer,
        uint16_t *out16_bit_buffer,
        uint32_t  inn_stride,
        uint32_t  out_stride,
        uint32_t  width,
        uint32_t  height);

    void c_pack_avx2_intrin(
        const uint8_t *inn_bit_buffer,
        uint32_t       inn_stride,
        uint8_t       *in_compn_bit_buffer,
        uint32_t       out_stride,
        uint8_t       *local_cache,
        uint32_t       width,
        uint32_t       height);

    void unpack_avg_avx2_intrin(
        uint16_t *ref16_l0,
        uint32_t  ref_l0_stride,
        uint16_t *ref16_l1,
        uint32_t  ref_l1_stride,
        uint8_t  *dst_ptr,
        uint32_t  dst_stride,
        uint32_t  width,
        uint32_t  height);

    int32_t sum_residual8bit_avx2_intrin(
        int16_t  *in_ptr,
        uint32_t  size,
        uint32_t  stride_in);

    void memset16bit_block_avx2_intrin(
        int16_t *in_ptr,
        uint32_t stride_in,
        uint32_t size,
        int16_t  value);

    void unpack_avg_safe_sub_avx2_intrin(
        uint16_t *ref16_l0,
        uint32_t  ref_l0_stride,
        uint16_t *ref16_l1,
        uint32_t  ref_l1_stride,
        uint8_t  *dst_ptr,
        uint32_t  dst_stride,
        EbBool    sub_pred,
        uint32_t  width,
        uint32_t  height);

    void picture_addition_kernel4x4_av1_sse2_intrin(
        uint8_t  *pred_ptr,
        uint32_t  pred_stride,
        int32_t  *residual_ptr,
        uint32_t  residual_stride,
        uint8_t  *recon_ptr,
        uint32_t  recon_stride,
        uint32_t  width,
        uint32_t  height,
        int32_t   bd);

    void picture_addition_kernel8x8_av1_sse2_intrin(
        uint8_t  *pred_ptr,
        uint32_t  pred_stride,
        int32_t  *residual_ptr,
        uint32_t  residual_stride,
        uint8_t  *recon_ptr,
        uint32_t  recon_stride,
        uint32_t  width,
        uint32_t  height,
        int32_t   bd);

    void picture_addition_kernel16x16_av1_sse2_intrin(
        uint8_t  *pred_ptr,
        uint32_t  pred_stride,
        int32_t  *residual_ptr,
        uint32_t  residual_stride,
        uint8_t  *recon_ptr,
        uint32_t  recon_stride,
        uint32_t  width,
        uint32_t  height,
        int32_t   bd);

    void picture_addition_kernel32x32_av1_sse2_intrin(
        uint8_t  *pred_ptr,
        uint32_t  pred_stride,
        int32_t  *residual_ptr,
        uint32_t  residual_stride,
        uint8_t  *recon_ptr,
        uint32_t  recon_stride,
        uint32_t  width,
        uint32_t  height,
        int32_t   bd);

    void picture_addition_kernel64x64_av1_sse2_intrin(
        uint8_t  *pred_ptr,
        uint32_t  pred_stride,
        int32_t  *residual_ptr,
        uint32_t  residual_stride,
        uint8_t  *recon_ptr,
        uint32_t  recon_stride,
        uint32_t  width,
        uint32_t  height,
        int32_t   bd);

    void full_distortion_kernel_cbf_zero32_bits_avx2(
        int32_t  *coeff,
        uint32_t  coeff_stride,
        int32_t  *recon_coeff,
        uint32_t  recon_coeff_stride,
        uint64_t  distortion_result[DIST_CALC_TOTAL],
        uint32_t  area_width,
        uint32_t  area_height);

    void full_distortion_kernel32_bits_avx2(
        int32_t  *coeff,
        uint32_t  coeff_stride,
        int32_t  *recon_coeff,
        uint32_t  recon_coeff_stride,
        uint64_t  distortion_result[DIST_CALC_TOTAL],
        uint32_t  area_width,
        uint32_t  area_height);


#ifdef __cplusplus
}
#endif
#endif // EbPictureOperators_AVX2