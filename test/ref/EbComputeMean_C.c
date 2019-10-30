/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

#include "EbDefinitions.h"

#define VARIANCE_PRECISION 16

typedef uint64_t (*EbComputeMeanFunc)(uint8_t* input_samples,
                                      uint32_t input_stride,
                                      uint32_t input_area_width,
                                      uint32_t input_area_height);

/*******************************************
 * compute_mean
 *   returns the mean of a block
 *******************************************/
uint64_t compute_mean(
    uint8_t* input_samples,     /**< input parameter, input samples Ptr */
    uint32_t input_stride,      /**< input parameter, input stride */
    uint32_t input_area_width,  /**< input parameter, input area width */
    uint32_t input_area_height) /**< input parameter, input area height */
{
    uint32_t hi, vi;
    uint64_t block_mean = 0;

    for (vi = 0; vi < input_area_height; vi++) {
        for (hi = 0; hi < input_area_width; hi++) {
            block_mean += input_samples[hi];
        }
        input_samples += input_stride;
    }

    block_mean = (block_mean << (VARIANCE_PRECISION >> 1)) /
                 (input_area_width * input_area_height);

    return block_mean;
}

/*******************************************
 * compute_mean_squared_values
 *   returns the Mean of Squared Values
 *******************************************/
uint64_t compute_mean_squared_values(
    uint8_t* input_samples,     /**< input parameter, input samples Ptr */
    uint32_t input_stride,      /**< input parameter, input stride */
    uint32_t input_area_width,  /**< input parameter, input area width */
    uint32_t input_area_height) /**< input parameter, input area height */
{
    uint32_t hi, vi;
    uint64_t block_mean = 0;

    for (vi = 0; vi < input_area_height; vi++) {
        for (hi = 0; hi < input_area_width; hi++) {
            block_mean += input_samples[hi] * input_samples[hi];
        }
        input_samples += input_stride;
    }

    block_mean = (block_mean << VARIANCE_PRECISION) /
                 (input_area_width * input_area_height);

    return block_mean;
}


