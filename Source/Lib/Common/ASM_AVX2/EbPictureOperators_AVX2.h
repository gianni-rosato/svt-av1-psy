/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPictureOperators_AVX2
#define EbPictureOperators_AVX2

#include <immintrin.h>
#include "EbDefinitions.h"
#include "EbPictureOperators_SSE2.h"

#ifdef __cplusplus
extern "C" {
#endif

    static INLINE void Distortion_AVX2_INTRIN(const __m256i input,
        const __m256i recon, __m256i *const sum) {
        const __m256i in = _mm256_unpacklo_epi8(input, _mm256_setzero_si256());
        const __m256i re = _mm256_unpacklo_epi8(recon, _mm256_setzero_si256());
        const __m256i diff = _mm256_sub_epi16(in, re);
        const __m256i dist = _mm256_madd_epi16(diff, diff);
        *sum = _mm256_add_epi32(*sum, dist);
    }

    static INLINE void SpatialFullDistortionKernel16_AVX2_INTRIN(
        const uint8_t *const input, const uint8_t *const recon,
        __m256i *const sum)
    {
        const __m128i in8 = _mm_loadu_si128((__m128i *)input);
        const __m128i re8 = _mm_loadu_si128((__m128i *)recon);
        const __m256i in16 = _mm256_cvtepu8_epi16(in8);
        const __m256i re16 = _mm256_cvtepu8_epi16(re8);
        const __m256i diff = _mm256_sub_epi16(in16, re16);
        const __m256i dist = _mm256_madd_epi16(diff, diff);
        *sum = _mm256_add_epi32(*sum, dist);
    }

    static INLINE void SpatialFullDistortionKernel32Leftover_AVX2_INTRIN(
        const uint8_t *const input, const uint8_t *const recon, __m256i *const sum0,
        __m256i *const sum1)
    {
        const __m256i in = _mm256_loadu_si256((__m256i *)input);
        const __m256i re = _mm256_loadu_si256((__m256i *)recon);
        const __m256i max = _mm256_max_epu8(in, re);
        const __m256i min = _mm256_min_epu8(in, re);
        const __m256i diff = _mm256_sub_epi8(max, min);
        const __m256i diff_L = _mm256_unpacklo_epi8(diff, _mm256_setzero_si256());
        const __m256i diff_H = _mm256_unpackhi_epi8(diff, _mm256_setzero_si256());
        const __m256i dist_L = _mm256_madd_epi16(diff_L, diff_L);
        const __m256i dist_H = _mm256_madd_epi16(diff_H, diff_H);
        *sum0 = _mm256_add_epi32(*sum0, dist_L);
        *sum1 = _mm256_add_epi32(*sum1, dist_H);
    }

    static INLINE void SpatialFullDistortionKernel32_AVX2_INTRIN(
        const uint8_t *const input, const uint8_t *const recon, __m256i *const sum)
    {
        const __m256i in = _mm256_loadu_si256((__m256i *)input);
        const __m256i re = _mm256_loadu_si256((__m256i *)recon);
        const __m256i max = _mm256_max_epu8(in, re);
        const __m256i min = _mm256_min_epu8(in, re);
        const __m256i diff = _mm256_sub_epi8(max, min);
        const __m256i diff_L = _mm256_unpacklo_epi8(diff, _mm256_setzero_si256());
        const __m256i diff_H = _mm256_unpackhi_epi8(diff, _mm256_setzero_si256());
        const __m256i dist_L = _mm256_madd_epi16(diff_L, diff_L);
        const __m256i dist_H = _mm256_madd_epi16(diff_H, diff_H);
        const __m256i dist = _mm256_add_epi32(dist_L, dist_H);
        *sum = _mm256_add_epi32(*sum, dist);
    }

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

    static INLINE int32_t Hadd32_AVX2_INTRIN(const __m256i src) {
        const __m128i src_L = _mm256_extracti128_si256(src, 0);
        const __m128i src_H = _mm256_extracti128_si256(src, 1);
        const __m128i sum = _mm_add_epi32(src_L, src_H);

        return Hadd32_SSE2_INTRIN(sum);
    }

    uint64_t spatial_full_distortion_kernel_avx2(
        uint8_t *input,
        uint32_t input_offset,
        uint32_t input_stride,
        uint8_t *recon,
        uint32_t recon_offset,
        uint32_t recon_stride,
        uint32_t area_width,
        uint32_t area_height);

    uint64_t spatial_full_distortion_kernel4x_n_avx2_intrin(
        uint8_t   *input,
        uint32_t   input_offset,
        uint32_t   input_stride,
        uint8_t   *recon,
        uint32_t   recon_offset,
        uint32_t   recon_stride,
        uint32_t   area_width,
        uint32_t   area_height);

    uint64_t spatial_full_distortion_kernel8x_n_avx2_intrin(
        uint8_t   *input,
        uint32_t   input_offset,
        uint32_t   input_stride,
        uint8_t   *recon,
        uint32_t   recon_offset,
        uint32_t   recon_stride,
        uint32_t   area_width,
        uint32_t   area_height);

    uint64_t spatial_full_distortion_kernel16x_n_avx2_intrin(
        uint8_t   *input,
        uint32_t   input_offset,
        uint32_t   input_stride,
        uint8_t   *recon,
        uint32_t   recon_offset,
        uint32_t   recon_stride,
        uint32_t   area_width,
        uint32_t   area_height);

    uint64_t spatial_full_distortion_kernel32x_n_avx2_intrin(
        uint8_t   *input,
        uint32_t   input_offset,
        uint32_t   input_stride,
        uint8_t   *recon,
        uint32_t   recon_offset,
        uint32_t   recon_stride,
        uint32_t   area_width,
        uint32_t   area_height);

    uint64_t spatial_full_distortion_kernel64x_n_avx2_intrin(
        uint8_t   *input,
        uint32_t   input_offset,
        uint32_t   input_stride,
        uint8_t   *recon,
        uint32_t   recon_offset,
        uint32_t   recon_stride,
        uint32_t   area_width,
        uint32_t   area_height);

    uint64_t spatial_full_distortion_kernel128x_n_avx2_intrin(
        uint8_t   *input,
        uint32_t   input_offset,
        uint32_t   input_stride,
        uint8_t   *recon,
        uint32_t   recon_offset,
        uint32_t   recon_stride,
        uint32_t   area_width,
        uint32_t   area_height);

#ifdef __cplusplus
}
#endif
#endif // EbPictureOperators_AVX2
