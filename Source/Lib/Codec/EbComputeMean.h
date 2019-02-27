/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbComputeMean_h
#define EbComputeMean_h

#include "EbComputeMean_SSE2.h"
#include "EbCombinedAveragingSAD_Intrinsic_AVX2.h"
#include "EbDefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif


    typedef uint64_t(*EB_COMPUTE_MEAN_FUNC)(
        uint8_t *input_samples,
        uint32_t input_stride,
        uint32_t input_area_width,
        uint32_t input_area_height);

    static const EB_COMPUTE_MEAN_FUNC ComputeMeanFunc[2][ASM_TYPE_TOTAL] = {
        {
            // NON_AVX2
            compute_mean8x8_sse2_intrin,
            // AVX2
            compute_mean8x8_avx2_intrin
        },
        {
            // NON_AVX2
            compute_mean_of_squared_values8x8_sse2_intrin,
            // AVX2
            compute_mean_of_squared_values8x8_sse2_intrin
        }
    };

#ifdef __cplusplus
}
#endif

#endif