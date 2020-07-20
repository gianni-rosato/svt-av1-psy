/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/
#ifndef EbComputeSAD_C_h
#define EbComputeSAD_C_h

#include "EbDefinitions.h"
#ifdef __cplusplus
extern "C" {
#endif

uint32_t fast_loop_nxm_sad_kernel(const uint8_t *src, // input parameter, source samples Ptr
                                  uint32_t       src_stride, // input parameter, source stride
                                  const uint8_t *ref, // input parameter, reference samples Ptr
                                  uint32_t       ref_stride, // input parameter, reference stride
                                  uint32_t       height, // input parameter, block height (M)
                                  uint32_t       width); // input parameter, block width (N)
#if !REMOVE_ME_SUBPEL_CODE

uint32_t combined_averaging_sad(uint8_t *src, uint32_t src_stride, uint8_t *ref1,
                                uint32_t ref1_stride, uint8_t *ref2, uint32_t ref2_stride,
                                uint32_t height, uint32_t width);
#endif
#if !REMOVE_UNUSED_CODE
void sad_loop_kernel_sparse_c(
    uint8_t * src, // input parameter, source samples Ptr
    uint32_t  srcStride, // input parameter, source stride
    uint8_t * ref, // input parameter, reference samples Ptr
    uint32_t  refStride, // input parameter, reference stride
    uint32_t  block_height, // input parameter, block height (M)
    uint32_t  block_width, // input parameter, block width (N)
    uint64_t *best_sad, int16_t *x_search_center, int16_t *y_search_center,
    uint32_t src_stride_raw, // input parameter, source stride (no line skipping)
    int16_t search_area_width, int16_t search_area_height);
#endif
void sad_loop_kernel_c(uint8_t * src, // input parameter, source samples Ptr
                       uint32_t  src_stride, // input parameter, source stride
                       uint8_t * ref, // input parameter, reference samples Ptr
                       uint32_t  ref_stride, // input parameter, reference stride
                       uint32_t  block_height, // input parameter, block height (M)
                       uint32_t  block_width, // input parameter, block width (N)
                       uint64_t *best_sad, int16_t *x_search_center, int16_t *y_search_center,
                       uint32_t src_stride_raw, // input parameter, source stride (no line skipping)
                       int16_t search_area_width, int16_t search_area_height);

#if RESTRUCTURE_SAD
void pme_sad_loop_kernel_c(uint8_t * src, // input parameter, source samples Ptr
                           uint32_t  src_stride, // input parameter, source stride
                           uint8_t * ref, // input parameter, reference samples Ptr
                           uint32_t  ref_stride, // input parameter, reference stride
                           uint32_t  block_height, // input parameter, block height (M)
                           uint32_t  block_width, // input parameter, block width (N)
                           uint32_t *best_sad, int16_t *best_mvx, int16_t *best_mvy,
                           int16_t search_position_start_x, int16_t search_position_start_y,
                           int16_t search_area_width, int16_t search_area_height, int16_t search_step, int16_t mvx,
                           int16_t mvy);
#endif

uint32_t nxm_sad_kernel_helper_c(const uint8_t *src, uint32_t src_stride, const uint8_t *ref,
                                 uint32_t ref_stride, uint32_t height, uint32_t width);

#if !REMOVE_ME_SUBPEL_CODE
uint32_t nxm_sad_avg_kernel_helper_c(uint8_t *src, uint32_t src_stride, uint8_t *ref1,
                                     uint32_t ref1_stride, uint8_t *ref2, uint32_t ref2_stride,
                                     uint32_t height, uint32_t width);
#endif
#ifdef __cplusplus
}
#endif
#endif // EbComputeSAD_C_h
