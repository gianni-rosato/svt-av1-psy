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
 * ComputeMean
 *   returns the mean of a block
 *******************************************/
uint64_t ComputeMean(
    uint8_t* inputSamples,     // input parameter, input samples Ptr
    uint32_t inputStride,      // input parameter, input stride
    uint32_t inputAreaWidth,   // input parameter, input area width
    uint32_t inputAreaHeight)  // input parameter, input area height
{
    uint32_t horizontalIndex;
    uint32_t verticalIndex;
    uint64_t blockMean = 0;

    for (verticalIndex = 0; verticalIndex < inputAreaHeight; verticalIndex++) {
        for (horizontalIndex = 0; horizontalIndex < inputAreaWidth;
             horizontalIndex++) {
            blockMean += inputSamples[horizontalIndex];
        }
        inputSamples += inputStride;
    }

    blockMean = (blockMean << (VARIANCE_PRECISION >> 1)) /
                (inputAreaWidth * inputAreaHeight);

    return blockMean;
}

/*******************************************
 * ComputeMeanOfSquaredValues
 *   returns the Mean of Squared Values
 *******************************************/
uint64_t ComputeMeanOfSquaredValues(
    uint8_t* inputSamples,     // input parameter, input samples Ptr
    uint32_t inputStride,      // input parameter, input stride
    uint32_t inputAreaWidth,   // input parameter, input area width
    uint32_t inputAreaHeight)  // input parameter, input area height
{
    uint32_t horizontalIndex;
    uint32_t verticalIndex;
    uint64_t blockMean = 0;

    for (verticalIndex = 0; verticalIndex < inputAreaHeight; verticalIndex++) {
        for (horizontalIndex = 0; horizontalIndex < inputAreaWidth;
             horizontalIndex++) {
            blockMean +=
                inputSamples[horizontalIndex] * inputSamples[horizontalIndex];
        }
        inputSamples += inputStride;
    }

    blockMean =
        (blockMean << VARIANCE_PRECISION) / (inputAreaWidth * inputAreaHeight);

    return blockMean;
}

uint64_t ComputeSubMeanOfSquaredValues(
    uint8_t* inputSamples,     // input parameter, input samples Ptr
    uint32_t inputStride,      // input parameter, input stride
    uint32_t inputAreaWidth,   // input parameter, input area width
    uint32_t inputAreaHeight)  // input parameter, input area height
{
    uint32_t horizontalIndex;
    uint32_t verticalIndex = 0;
    uint64_t blockMean = 0;
    uint16_t skip = 0;

    for (verticalIndex = 0; skip < inputAreaHeight;
         skip = verticalIndex + verticalIndex) {
        for (horizontalIndex = 0; horizontalIndex < inputAreaWidth;
             horizontalIndex++) {
            blockMean +=
                inputSamples[horizontalIndex] * inputSamples[horizontalIndex];
        }
        inputSamples += 2 * inputStride;
        verticalIndex++;
    }

    blockMean =
        blockMean
        << 11;  // VARIANCE_PRECISION) / (inputAreaWidth * inputAreaHeight);

    return blockMean;
}

uint64_t ComputeSubMean8x8(
    uint8_t* inputSamples,     // input parameter, input samples Ptr
    uint32_t inputStride,      // input parameter, input stride
    uint32_t inputAreaWidth,   // input parameter, input area width
    uint32_t inputAreaHeight)  // input parameter, input area height
{
    uint32_t horizontalIndex;
    uint32_t verticalIndex = 0;
    uint64_t blockMean = 0;
    uint16_t skip = 0;

    for (verticalIndex = 0; skip < inputAreaHeight;
         skip = verticalIndex + verticalIndex) {
        for (horizontalIndex = 0; horizontalIndex < inputAreaWidth;
             horizontalIndex++) {
            blockMean += inputSamples[horizontalIndex];
        }
        inputSamples += 2 * inputStride;
        verticalIndex++;
    }

    blockMean = blockMean << 3;  // (VARIANCE_PRECISION >> 1)) /
                                 // (inputAreaWidth * inputAreaHeight/2)

    return blockMean;
}
