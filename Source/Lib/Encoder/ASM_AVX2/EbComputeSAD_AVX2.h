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


uint32_t eb_compute8x_m_sad_avx2_intrin(
    const uint8_t *src, // input parameter, source samples Ptr
    uint32_t       src_stride, // input parameter, source stride
    const uint8_t *ref, // input parameter, reference samples Ptr
    uint32_t       ref_stride, // input parameter, reference stride
    uint32_t       height, // input parameter, block height (M)
    uint32_t       width); // input parameter, block width (N)

uint32_t eb_compute16x_m_sad_avx2_intrin(
    const uint8_t *src, // input parameter, source samples Ptr
    uint32_t       src_stride, // input parameter, source stride
    const uint8_t *ref, // input parameter, reference samples Ptr
    uint32_t       ref_stride, // input parameter, reference stride
    uint32_t       height, // input parameter, block height (M)
    uint32_t       width); // input parameter, block width (N)

uint32_t eb_compute24x_m_sad_avx2_intrin(
    const uint8_t *src, // input parameter, source samples Ptr
    uint32_t       src_stride, // input parameter, source stride
    const uint8_t *ref, // input parameter, reference samples Ptr
    uint32_t       ref_stride, // input parameter, reference stride
    uint32_t       height, // input parameter, block height (M)
    uint32_t       width); // input parameter, block width (N)

uint32_t eb_compute32x_m_sad_avx2_intrin(
    const uint8_t *src, // input parameter, source samples Ptr
    uint32_t       src_stride, // input parameter, source stride
    const uint8_t *ref, // input parameter, reference samples Ptr
    uint32_t       ref_stride, // input parameter, reference stride
    uint32_t       height, // input parameter, block height (M)
    uint32_t       width); // input parameter, block width (N)

uint32_t eb_compute48x_m_sad_avx2_intrin(
    const uint8_t *src, // input parameter, source samples Ptr
    uint32_t       src_stride, // input parameter, source stride
    const uint8_t *ref, // input parameter, reference samples Ptr
    uint32_t       ref_stride, // input parameter, reference stride
    uint32_t       height, // input parameter, block height (M)
    uint32_t       width); // input parameter, block width (N)

uint32_t eb_compute64x_m_sad_avx2_intrin(
    const uint8_t *src, // input parameter, source samples Ptr
    uint32_t       src_stride, // input parameter, source stride
    const uint8_t *ref, // input parameter, reference samples Ptr
    uint32_t       ref_stride, // input parameter, reference stride
    uint32_t       height, // input parameter, block height (M)
    uint32_t       width); // input parameter, block width (N)

uint32_t eb_compute128x_m_sad_avx2_intrin(
    const uint8_t *src, // input parameter, source samples Ptr
    uint32_t       src_stride, // input parameter, source stride
    const uint8_t *ref, // input parameter, reference samples Ptr
    uint32_t       ref_stride, // input parameter, reference stride
    uint32_t       height, // input parameter, block height (M)
    uint32_t       width); // input parameter, block width (N)


#if RESTRUCTURE_SAD
void pme_sad_loop_kernel_avx2(uint8_t * src, // input parameter, source samples Ptr
    uint32_t  src_stride, // input parameter, source stride
    uint8_t * ref, // input parameter, reference samples Ptr
    uint32_t  ref_stride, // input parameter, reference stride
    uint32_t  block_height, // input parameter, block height (M)
    uint32_t  block_width, // input parameter, block width (N)
    uint32_t *best_sad, int16_t *best_mvx, int16_t *best_mvy,
    int16_t search_position_start_x, int16_t search_position_start_y,
    int16_t search_area_width, int16_t search_area_height,
    int16_t search_step, int16_t mvx, int16_t mvy);
#endif

#ifdef __cplusplus
}
#endif
#endif // EbComputeSAD_AVX2_h
