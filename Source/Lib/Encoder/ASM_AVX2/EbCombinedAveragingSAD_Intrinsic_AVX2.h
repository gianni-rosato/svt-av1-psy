/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/
#ifndef EbCombinedAveragingSAD_AVX2_h
#define EbCombinedAveragingSAD_AVX2_h

#include "EbDefinitions.h"
#ifdef __cplusplus
extern "C" {
#endif

    uint32_t combined_averaging_8xm_sad_avx2_intrin(
        uint8_t  *src,
        uint32_t  src_stride,
        uint8_t  *ref1,
        uint32_t  ref1_stride,
        uint8_t  *ref2,
        uint32_t  ref2_stride,
        uint32_t  height,
        uint32_t  width);

    uint32_t combined_averaging_16xm_sad_avx2_intrin(
        uint8_t  *src,
        uint32_t  src_stride,
        uint8_t  *ref1,
        uint32_t  ref1_stride,
        uint8_t  *ref2,
        uint32_t  ref2_stride,
        uint32_t  height,
        uint32_t  width);

    uint32_t combined_averaging_24xm_sad_avx2_intrin(
        uint8_t  *src,
        uint32_t  src_stride,
        uint8_t  *ref1,
        uint32_t  ref1_stride,
        uint8_t  *ref2,
        uint32_t  ref2_stride,
        uint32_t  height,
        uint32_t  width);

    uint32_t combined_averaging_32xm_sad_avx2_intrin(
        uint8_t  *src,
        uint32_t  src_stride,
        uint8_t  *ref1,
        uint32_t  ref1_stride,
        uint8_t  *ref2,
        uint32_t  ref2_stride,
        uint32_t  height,
        uint32_t  width);

    uint32_t combined_averaging_48xm_sad_avx2_intrin(
        uint8_t  *src,
        uint32_t  src_stride,
        uint8_t  *ref1,
        uint32_t  ref1_stride,
        uint8_t  *ref2,
        uint32_t  ref2_stride,
        uint32_t  height,
        uint32_t  width);

    uint32_t combined_averaging_64xm_sad_avx2_intrin(
        uint8_t  *src,
        uint32_t  src_stride,
        uint8_t  *ref1,
        uint32_t  ref1_stride,
        uint8_t  *ref2,
        uint32_t  ref2_stride,
        uint32_t  height,
        uint32_t  width);

    uint64_t compute_mean8x8_avx2_intrin(
        uint8_t  *input_samples,      // input parameter, input samples Ptr
        uint32_t  input_stride,       // input parameter, input stride
        uint32_t  input_area_width,   // input parameter, input area width
        uint32_t  input_area_height);

    void compute_interm_var_four8x8_avx2_intrin(
        uint8_t  *input_samples,
        uint16_t  input_stride,
        uint64_t *mean_of8x8_blocks,      // mean of four  8x8
        uint64_t *mean_of_squared8x8_blocks);

    uint32_t nxm_sad_avg_kernel_helper_avx2(
        uint8_t  *src,
        uint32_t  src_stride,
        uint8_t  *ref1,
        uint32_t  ref1_stride,
        uint8_t  *ref2,
        uint32_t  ref2_stride,
        uint32_t  height,
        uint32_t  width);

#ifdef __cplusplus
}
#endif
#endif
