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

uint32_t combined_averaging_8xm_sad_avx2_intrin(uint8_t *src, uint32_t src_stride, uint8_t *ref1,
                                                uint32_t ref1_stride, uint8_t *ref2,
                                                uint32_t ref2_stride, uint32_t height,
                                                uint32_t width);

uint32_t combined_averaging_16xm_sad_avx2_intrin(uint8_t *src, uint32_t src_stride, uint8_t *ref1,
                                                 uint32_t ref1_stride, uint8_t *ref2,
                                                 uint32_t ref2_stride, uint32_t height,
                                                 uint32_t width);

uint32_t combined_averaging_24xm_sad_avx2_intrin(uint8_t *src, uint32_t src_stride, uint8_t *ref1,
                                                 uint32_t ref1_stride, uint8_t *ref2,
                                                 uint32_t ref2_stride, uint32_t height,
                                                 uint32_t width);

uint32_t combined_averaging_32xm_sad_avx2_intrin(uint8_t *src, uint32_t src_stride, uint8_t *ref1,
                                                 uint32_t ref1_stride, uint8_t *ref2,
                                                 uint32_t ref2_stride, uint32_t height,
                                                 uint32_t width);

uint32_t combined_averaging_48xm_sad_avx2_intrin(uint8_t *src, uint32_t src_stride, uint8_t *ref1,
                                                 uint32_t ref1_stride, uint8_t *ref2,
                                                 uint32_t ref2_stride, uint32_t height,
                                                 uint32_t width);

uint32_t combined_averaging_64xm_sad_avx2_intrin(uint8_t *src, uint32_t src_stride, uint8_t *ref1,
                                                 uint32_t ref1_stride, uint8_t *ref2,
                                                 uint32_t ref2_stride, uint32_t height,
                                                 uint32_t width);

#ifdef __cplusplus
}
#endif
#endif
