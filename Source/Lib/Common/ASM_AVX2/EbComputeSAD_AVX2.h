/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/
#ifndef EbComputeSAD_AVX2_h
#define EbComputeSAD_AVX2_h

#include "EbDefinitions.h"
#ifdef __cplusplus
extern "C" {
#endif

    void ext_sad_calculation_8x8_16x16_avx2_intrin(
        uint8_t  *src,
        uint32_t  src_stride,
        uint8_t  *ref,
        uint32_t  ref_stride,
        uint32_t *p_best_sad8x8,
        uint32_t *p_best_sad16x16,
        uint32_t *p_best_mv8x8,
        uint32_t *p_best_mv16x16,
        uint32_t  mv,
        uint32_t *p_sad16x16,
        uint32_t *p_sad8x8,
        EbBool    sub_sad);

    uint32_t eb_compute4x_m_sad_avx2_intrin(
        const uint8_t  *src,                      // input parameter, source samples Ptr
        uint32_t  src_stride,                     // input parameter, source stride
        const uint8_t  *ref,                      // input parameter, reference samples Ptr
        uint32_t  ref_stride,                     // input parameter, reference stride
        uint32_t  height,                         // input parameter, block height (M)
        uint32_t  width);                         // input parameter, block width (N)

    uint32_t eb_compute8x_m_sad_avx2_intrin(
        const uint8_t  *src,                      // input parameter, source samples Ptr
        uint32_t  src_stride,                     // input parameter, source stride
        const uint8_t  *ref,                      // input parameter, reference samples Ptr
        uint32_t  ref_stride,                     // input parameter, reference stride
        uint32_t  height,                         // input parameter, block height (M)
        uint32_t  width);                         // input parameter, block width (N)

    uint32_t eb_compute16x_m_sad_avx2_intrin(
        const uint8_t  *src,                      // input parameter, source samples Ptr
        uint32_t  src_stride,                     // input parameter, source stride
        const uint8_t  *ref,                      // input parameter, reference samples Ptr
        uint32_t  ref_stride,                     // input parameter, reference stride
        uint32_t  height,                         // input parameter, block height (M)
        uint32_t  width);                         // input parameter, block width (N)

    uint32_t eb_compute24x_m_sad_avx2_intrin(
        const uint8_t  *src,                      // input parameter, source samples Ptr
        uint32_t  src_stride,                     // input parameter, source stride
        const uint8_t  *ref,                      // input parameter, reference samples Ptr
        uint32_t  ref_stride,                     // input parameter, reference stride
        uint32_t  height,                         // input parameter, block height (M)
        uint32_t  width);                         // input parameter, block width (N)

    uint32_t eb_compute32x_m_sad_avx2_intrin(
        const uint8_t  *src,                      // input parameter, source samples Ptr
        uint32_t  src_stride,                     // input parameter, source stride
        const uint8_t  *ref,                      // input parameter, reference samples Ptr
        uint32_t  ref_stride,                     // input parameter, reference stride
        uint32_t  height,                         // input parameter, block height (M)
        uint32_t  width);                         // input parameter, block width (N)

    uint32_t eb_compute48x_m_sad_avx2_intrin(
        const uint8_t  *src,                      // input parameter, source samples Ptr
        uint32_t  src_stride,                     // input parameter, source stride
        const uint8_t  *ref,                      // input parameter, reference samples Ptr
        uint32_t  ref_stride,                     // input parameter, reference stride
        uint32_t  height,                         // input parameter, block height (M)
        uint32_t  width);                         // input parameter, block width (N)

    uint32_t eb_compute64x_m_sad_avx2_intrin(
        const uint8_t  *src,                      // input parameter, source samples Ptr
        uint32_t  src_stride,                     // input parameter, source stride
        const uint8_t  *ref,                      // input parameter, reference samples Ptr
        uint32_t  ref_stride,                     // input parameter, reference stride
        uint32_t  height,                         // input parameter, block height (M)
        uint32_t  width);                         // input parameter, block width (N)

    void sad_loop_kernel_avx2_hme_l0_intrin(
        uint8_t  *src,                            // input parameter, source samples Ptr
        uint32_t  src_stride,                     // input parameter, source stride
        uint8_t  *ref,                            // input parameter, reference samples Ptr
        uint32_t  ref_stride,                     // input parameter, reference stride
        uint32_t  block_height,                   // input parameter, block height (M)
        uint32_t  block_width,                    // input parameter, block width (N)
        uint64_t *best_sad,
        int16_t  *x_search_center,
        int16_t  *y_search_center,
        uint32_t  src_stride_raw,                 // input parameter, source stride (no line skipping)
        int16_t   search_area_width,
        int16_t   search_area_height);

    void get_eight_horizontal_search_point_results_8x8_16x16_pu_avx2_intrin(
        uint8_t   *src,
        uint32_t   src_stride,
        uint8_t   *ref,
        uint32_t   ref_stride,
        uint32_t  *p_best_sad8x8,
        uint32_t  *p_best_mv8x8,
        uint32_t  *p_best_sad16x16,
        uint32_t  *p_best_mv16x16,
        uint32_t   mv,
        uint16_t  *p_sad16x16,
        EbBool     sub_sad);

    void get_eight_horizontal_search_point_results_32x32_64x64_pu_avx2_intrin(
        uint16_t  *p_sad16x16,
        uint32_t  *p_best_sad32x32,
        uint32_t  *p_best_sad64x64,
        uint32_t  *p_best_mv32x32,
        uint32_t  *p_best_mv64x64,
        uint32_t   mv);

    void sad_loop_kernel_sparse_avx2_intrin(
        uint8_t  *src,                            // input parameter, source samples Ptr
        uint32_t  src_stride,                     // input parameter, source stride
        uint8_t  *ref,                            // input parameter, reference samples Ptr
        uint32_t  ref_stride,                     // input parameter, reference stride
        uint32_t  block_height,                   // input parameter, block height (M)
        uint32_t  block_width,                    // input parameter, block width (N)
        uint64_t *best_sad,
        int16_t  *x_search_center,
        int16_t  *y_search_center,
        uint32_t  src_stride_raw,                 // input parameter, source stride (no line skipping)
        int16_t   search_area_width,
        int16_t   search_area_height);

    void ext_all_sad_calculation_8x8_16x16_avx2(
        uint8_t   *src,
        uint32_t   src_stride,
        uint8_t   *ref,
        uint32_t   ref_stride,
        uint32_t   mv,
        uint32_t  *p_best_sad8x8,
        uint32_t  *p_best_sad16x16,
        uint32_t  *p_best_mv8x8,
        uint32_t  *p_best_mv16x16,
        uint32_t   p_eight_sad16x16[16][8],
        uint32_t   p_eight_sad8x8[64][8]);

    void ext_eigth_sad_calculation_nsq_avx2(
        uint32_t   p_sad8x8[64][8],
        uint32_t   p_sad16x16[16][8],
        uint32_t   p_sad32x32[4][8],
        uint32_t  *p_best_sad64x32,
        uint32_t  *p_best_mv64x32,
        uint32_t  *p_best_sad32x16,
        uint32_t  *p_best_mv32x16,
        uint32_t  *p_best_sad16x8,
        uint32_t  *p_best_mv16x8,
        uint32_t  *p_best_sad32x64,
        uint32_t  *p_best_mv32x64,
        uint32_t  *p_best_sad16x32,
        uint32_t  *p_best_mv16x32,
        uint32_t  *p_best_sad8x16,
        uint32_t  *p_best_mv8x16,
        uint32_t  *p_best_sad32x8,
        uint32_t  *p_best_mv32x8,
        uint32_t  *p_best_sad8x32,
        uint32_t  *p_best_mv8x32,
        uint32_t  *p_best_sad64x16,
        uint32_t  *p_best_mv64x16,
        uint32_t  *p_best_sad16x64,
        uint32_t  *p_best_mv16x64,
        uint32_t   mv);

    void ext_eight_sad_calculation_32x32_64x64_avx2(
        uint32_t  p_sad16x16[16][8],
        uint32_t *p_best_sad32x32,
        uint32_t *p_best_sad64x64,
        uint32_t *p_best_mv32x32,
        uint32_t *p_best_mv64x64,
        uint32_t  mv,
        uint32_t  p_sad32x32[4][8]);

    uint32_t nxm_sad_kernel_sub_sampled_helper_avx2(
        const uint8_t  *src,
        uint32_t  src_stride,
        const uint8_t  *ref,
        uint32_t  ref_stride,
        uint32_t  height,
        uint32_t  width);

    uint32_t nxm_sad_kernel_helper_avx2(
        const uint8_t  *src,
        uint32_t  src_stride,
        const uint8_t  *ref,
        uint32_t  ref_stride,
        uint32_t  height,
        uint32_t  width);

#ifdef __cplusplus
}
#endif
#endif // EbComputeSAD_AVX2_h
