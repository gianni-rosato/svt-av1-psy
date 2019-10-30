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

    uint64_t compute_sub_mean_c(
        uint8_t* input_samples,     /**< input parameter, input samples Ptr */
        uint32_t input_stride,      /**< input parameter, input stride */
        uint32_t input_area_width,  /**< input parameter, input area width */
        uint32_t input_area_height);

    uint64_t compute_sub_mean_squared_values_c(
        uint8_t* input_samples,     /**< input parameter, input samples Ptr */
        uint32_t input_stride,      /**< input parameter, input stride */
        uint32_t input_area_width,  /**< input parameter, input area width */
        uint32_t input_area_height); /**< input parameter, input area height */

    void compute_interm_var_four8x8_c(
        uint8_t *  input_samples,
        uint16_t   input_stride,
        uint64_t * mean_of8x8_blocks,      // mean of four  8x8
        uint64_t * mean_of_squared8x8_blocks);  // meanSquared

#ifdef __cplusplus
}
#endif

#endif
