/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include <string.h>

#include "aom_dsp_rtcd.h"
#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#include "EbPictureControlSet.h"
#include "EbSequenceControlSet.h"
#include "EbPictureBufferDesc.h"

#include "EbResourceCoordinationResults.h"
#include "EbPictureAnalysisProcess.h"
#include "EbPictureAnalysisResults.h"
#include "EbMcp.h"
#include "EbMotionEstimation.h"
#include "EbReferenceObject.h"

#include "EbComputeMean.h"
#include "EbMeSadCalculation.h"
#include "EbComputeMean_SSE2.h"
#include "EbCombinedAveragingSAD_Intrinsic_AVX2.h"

#define VARIANCE_PRECISION        16
#define  LCU_LOW_VAR_TH                5
#define  PIC_LOW_VAR_PERCENTAGE_TH    60
#define    FLAT_MAX_VAR            50
#define FLAT_MAX_VAR_DECIM        (50-00)
#define    NOISE_MIN_LEVEL            70000//120000
#define NOISE_MIN_LEVEL_DECIM   (70000+000000)//(120000+000000)
#define    NOISE_MIN_LEVEL_M6_M7       120000
#define NOISE_MIN_LEVEL_DECIM_M6_M7    (120000+000000)
#define DENOISER_QP_TH            29
#define DENOISER_BITRATE_TH        14000000
#define SAMPLE_THRESHOLD_PRECENT_BORDER_LINE      15
#define SAMPLE_THRESHOLD_PRECENT_TWO_BORDER_LINES 10

static void picture_analysis_context_dctor(EbPtr p)
{
    PictureAnalysisContext *obj = (PictureAnalysisContext*)p;
    EB_DELETE(obj->noise_picture_ptr);
    EB_DELETE(obj->denoised_picture_ptr);
}
/************************************************
* Picture Analysis Context Constructor
************************************************/
EbErrorType picture_analysis_context_ctor(
    PictureAnalysisContext *context_ptr,
    EbPictureBufferDescInitData * input_picture_buffer_desc_init_data,
    EbBool                         denoise_flag,
    EbFifo *resource_coordination_results_input_fifo_ptr,
    EbFifo *picture_analysis_results_output_fifo_ptr)
{
    context_ptr->resource_coordination_results_input_fifo_ptr = resource_coordination_results_input_fifo_ptr;
    context_ptr->picture_analysis_results_output_fifo_ptr = picture_analysis_results_output_fifo_ptr;

    context_ptr->dctor = picture_analysis_context_dctor;

    if (denoise_flag == EB_TRUE) {
        //denoised
        // If 420/422, re-use luma for chroma
        // If 444, re-use luma for Cr
        if (input_picture_buffer_desc_init_data->color_format != EB_YUV444) {
            input_picture_buffer_desc_init_data->buffer_enable_mask = PICTURE_BUFFER_DESC_Y_FLAG;
        } else
            input_picture_buffer_desc_init_data->buffer_enable_mask = PICTURE_BUFFER_DESC_Y_FLAG | PICTURE_BUFFER_DESC_Cb_FLAG;
        EB_NEW(
            context_ptr->denoised_picture_ptr,
            eb_picture_buffer_desc_ctor,
            (EbPtr)input_picture_buffer_desc_init_data);

        if (input_picture_buffer_desc_init_data->color_format != EB_YUV444) {
            context_ptr->denoised_picture_ptr->buffer_cb = context_ptr->denoised_picture_ptr->buffer_y;
            context_ptr->denoised_picture_ptr->buffer_cr = context_ptr->denoised_picture_ptr->buffer_y + context_ptr->denoised_picture_ptr->chroma_size;
        } else
            context_ptr->denoised_picture_ptr->buffer_cr = context_ptr->denoised_picture_ptr->buffer_y;
        // noise
        input_picture_buffer_desc_init_data->max_height = BLOCK_SIZE_64;
        input_picture_buffer_desc_init_data->buffer_enable_mask = PICTURE_BUFFER_DESC_Y_FLAG;

        EB_NEW(
            context_ptr->noise_picture_ptr,
            eb_picture_buffer_desc_ctor,
            (EbPtr)input_picture_buffer_desc_init_data);
    }
    return EB_ErrorNone;
}
void DownSampleChroma(EbPictureBufferDesc* input_picture_ptr, EbPictureBufferDesc* outputPicturePtr)
{
    uint32_t input_color_format = input_picture_ptr->color_format;
    const uint16_t input_subsampling_x = (input_color_format == EB_YUV444 ? 1 : 2) - 1;
    const uint16_t input_subsampling_y = (input_color_format >= EB_YUV422 ? 1 : 2) - 1;

    uint32_t output_color_format = outputPicturePtr->color_format;
    const uint16_t output_subsampling_x = (output_color_format == EB_YUV444 ? 1 : 2) - 1;
    const uint16_t output_subsampling_y = (output_color_format >= EB_YUV422 ? 1 : 2) - 1;

    uint32_t stride_in, strideOut;
    uint32_t inputOriginIndex, outputOriginIndex;

    uint8_t *ptrIn;
    uint8_t *ptrOut;

    uint32_t ii, jj;

    //Cb
    {
        stride_in = input_picture_ptr->stride_cb;
        inputOriginIndex = (input_picture_ptr->origin_x >> input_subsampling_x) +
            (input_picture_ptr->origin_y >> input_subsampling_y)  * input_picture_ptr->stride_cb;
        ptrIn = &(input_picture_ptr->buffer_cb[inputOriginIndex]);

        strideOut = outputPicturePtr->stride_cb;
        outputOriginIndex = (outputPicturePtr->origin_x >> output_subsampling_x) +
            (outputPicturePtr->origin_y >> output_subsampling_y)  * outputPicturePtr->stride_cb;
        ptrOut = &(outputPicturePtr->buffer_cb[outputOriginIndex]);

        for (jj = 0; jj < (uint32_t)(outputPicturePtr->height >> output_subsampling_y); jj++) {
            for (ii = 0; ii < (uint32_t)(outputPicturePtr->width >> output_subsampling_x); ii++) {
                ptrOut[ii + jj * strideOut] =
                    ptrIn[(ii << (1 - input_subsampling_x)) +
                    (jj << (1 - input_subsampling_y)) * stride_in];
            }
        }
    }

    //Cr
    {
        stride_in = input_picture_ptr->stride_cr;
        inputOriginIndex = (input_picture_ptr->origin_x >> input_subsampling_x) + (input_picture_ptr->origin_y >> input_subsampling_y)  * input_picture_ptr->stride_cr;
        ptrIn = &(input_picture_ptr->buffer_cr[inputOriginIndex]);

        strideOut = outputPicturePtr->stride_cr;
        outputOriginIndex = (outputPicturePtr->origin_x >> output_subsampling_x) + (outputPicturePtr->origin_y >> output_subsampling_y)  * outputPicturePtr->stride_cr;
        ptrOut = &(outputPicturePtr->buffer_cr[outputOriginIndex]);

        for (jj = 0; jj < (uint32_t)(outputPicturePtr->height >> output_subsampling_y); jj++) {
            for (ii = 0; ii < (uint32_t)(outputPicturePtr->width >> output_subsampling_x); ii++) {
                ptrOut[ii + jj * strideOut] =
                    ptrIn[(ii << (1 - input_subsampling_x)) +
                    (jj << (1 - input_subsampling_y)) * stride_in];
            }
        }
    }
}

/************************************************
 * Picture Analysis Context Destructor
 ************************************************/
  /********************************************
    * decimation_2d
    *      decimates the input
    ********************************************/
void decimation_2d(
    uint8_t *  input_samples,      // input parameter, input samples Ptr
    uint32_t   input_stride,       // input parameter, input stride
    uint32_t   input_area_width,   // input parameter, input area width
    uint32_t   input_area_height,  // input parameter, input area height
    uint8_t *  decim_samples,      // output parameter, decimated samples Ptr
    uint32_t   decim_stride,       // input parameter, output stride
    uint32_t   decim_step)         // input parameter, decimation amount in pixels
{
    uint32_t horizontal_index;
    uint32_t vertical_index;
    uint32_t input_stripe_stride = input_stride * decim_step;

    for (vertical_index = 0; vertical_index < input_area_height; vertical_index += decim_step) {
        for (horizontal_index = 0; horizontal_index < input_area_width; horizontal_index += decim_step)
            decim_samples[(horizontal_index >> (decim_step >> 1))] = input_samples[horizontal_index];

        input_samples += input_stripe_stride;
        decim_samples += decim_stride;
    }

    return;
}

/********************************************
 * downsample_2d
 *      downsamples the input
 * Alternative implementation to decimation_2d that performs filtering (2x2, 0-phase)
 ********************************************/
void downsample_2d(
    uint8_t *  input_samples,      // input parameter, input samples Ptr
    uint32_t   input_stride,       // input parameter, input stride
    uint32_t   input_area_width,    // input parameter, input area width
    uint32_t   input_area_height,   // input parameter, input area height
    uint8_t *  decim_samples,      // output parameter, decimated samples Ptr
    uint32_t   decim_stride,       // input parameter, output stride
    uint32_t   decim_step)        // input parameter, decimation amount in pixels
{

    uint32_t horizontal_index;
    uint32_t vertical_index;
    uint32_t input_stripe_stride = input_stride * decim_step;
    uint32_t decim_horizontal_index;
    const uint32_t half_decim_step = decim_step >> 1;

    for (input_samples += half_decim_step * input_stride, vertical_index = half_decim_step; vertical_index < input_area_height; vertical_index += decim_step) {
        uint8_t *prev_input_line = input_samples - input_stride;
        for (horizontal_index = half_decim_step, decim_horizontal_index = 0; horizontal_index < input_area_width; horizontal_index += decim_step, decim_horizontal_index++) {
            uint32_t sum = (uint32_t)prev_input_line[horizontal_index - 1] + (uint32_t)prev_input_line[horizontal_index] + (uint32_t)input_samples[horizontal_index - 1] + (uint32_t)input_samples[horizontal_index];
            decim_samples[decim_horizontal_index] = (sum + 2) >> 2;

        }
        input_samples += input_stripe_stride;
        decim_samples += decim_stride;
    }

    return;
}

/********************************************
* CalculateHistogram
*      creates n-bins histogram for the input
********************************************/
void CalculateHistogram(
    uint8_t *  input_samples,      // input parameter, input samples Ptr
    uint32_t   input_area_width,    // input parameter, input area width
    uint32_t   input_area_height,   // input parameter, input area height
    uint32_t   stride,            // input parameter, input stride
    uint8_t    decim_step,         // input parameter, area height
    uint32_t  *histogram,            // output parameter, output histogram
    uint64_t  *sum)
{
    uint32_t horizontal_index;
    uint32_t vertical_index;
    *sum = 0;

    for (vertical_index = 0; vertical_index < input_area_height; vertical_index += decim_step) {
        for (horizontal_index = 0; horizontal_index < input_area_width; horizontal_index += decim_step) {
            ++(histogram[input_samples[horizontal_index]]);
            *sum += input_samples[horizontal_index];
        }
        input_samples += (stride << (decim_step >> 1));
    }

    return;
}

uint64_t ComputeVariance32x32(
    EbPictureBufferDesc       *input_padded_picture_ptr,         // input parameter, Input Padded Picture
    uint32_t                       inputLumaOriginIndex,          // input parameter, SB index, used to point to source/reference samples
    uint64_t                        *variance8x8) {
    uint32_t blockIndex;

    uint64_t mean_of8x8_blocks[16];
    uint64_t meanOf8x8SquaredValuesBlocks[16];

    uint64_t meanOf16x16Blocks[4];
    uint64_t meanOf16x16SquaredValuesBlocks[4];

    uint64_t meanOf32x32Blocks;
    uint64_t meanOf32x32SquaredValuesBlocks;
    /////////////////////////////////////////////
    // (0,0)
    blockIndex = inputLumaOriginIndex;

    mean_of8x8_blocks[0] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
    meanOf8x8SquaredValuesBlocks[0] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

    // (0,1)
    blockIndex = blockIndex + 8;
    mean_of8x8_blocks[1] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
    meanOf8x8SquaredValuesBlocks[1] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

    // (0,2)
    blockIndex = blockIndex + 8;
    mean_of8x8_blocks[2] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
    meanOf8x8SquaredValuesBlocks[2] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

    // (0,3)
    blockIndex = blockIndex + 8;
    mean_of8x8_blocks[3] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
    meanOf8x8SquaredValuesBlocks[3] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

    // (1,0)
    blockIndex = inputLumaOriginIndex + (input_padded_picture_ptr->stride_y << 3);
    mean_of8x8_blocks[4] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
    meanOf8x8SquaredValuesBlocks[4] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

    // (1,1)
    blockIndex = blockIndex + 8;
    mean_of8x8_blocks[5] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
    meanOf8x8SquaredValuesBlocks[5] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

    // (1,2)
    blockIndex = blockIndex + 8;
    mean_of8x8_blocks[6] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
    meanOf8x8SquaredValuesBlocks[6] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

    // (1,3)
    blockIndex = blockIndex + 8;
    mean_of8x8_blocks[7] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
    meanOf8x8SquaredValuesBlocks[7] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

    // (2,0)
    blockIndex = inputLumaOriginIndex + (input_padded_picture_ptr->stride_y << 4);
    mean_of8x8_blocks[8] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
    meanOf8x8SquaredValuesBlocks[8] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

    // (2,1)
    blockIndex = blockIndex + 8;
    mean_of8x8_blocks[9] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
    meanOf8x8SquaredValuesBlocks[9] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

    // (2,2)
    blockIndex = blockIndex + 8;
    mean_of8x8_blocks[10] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
    meanOf8x8SquaredValuesBlocks[10] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

    // (2,3)
    blockIndex = blockIndex + 8;
    mean_of8x8_blocks[11] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
    meanOf8x8SquaredValuesBlocks[11] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

    // (3,0)
    blockIndex = inputLumaOriginIndex + (input_padded_picture_ptr->stride_y << 3) + (input_padded_picture_ptr->stride_y << 4);
    mean_of8x8_blocks[12] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
    meanOf8x8SquaredValuesBlocks[12] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

    // (3,1)
    blockIndex = blockIndex + 8;
    mean_of8x8_blocks[13] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
    meanOf8x8SquaredValuesBlocks[13] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

    // (3,2)
    blockIndex = blockIndex + 8;
    mean_of8x8_blocks[14] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
    meanOf8x8SquaredValuesBlocks[14] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

    // (3,3)
    blockIndex = blockIndex + 8;
    mean_of8x8_blocks[15] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
    meanOf8x8SquaredValuesBlocks[15] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

    /////////////////////////////////////////////

    variance8x8[0] = meanOf8x8SquaredValuesBlocks[0] - (mean_of8x8_blocks[0] * mean_of8x8_blocks[0]);
    variance8x8[1] = meanOf8x8SquaredValuesBlocks[1] - (mean_of8x8_blocks[1] * mean_of8x8_blocks[1]);
    variance8x8[2] = meanOf8x8SquaredValuesBlocks[2] - (mean_of8x8_blocks[2] * mean_of8x8_blocks[2]);
    variance8x8[3] = meanOf8x8SquaredValuesBlocks[3] - (mean_of8x8_blocks[3] * mean_of8x8_blocks[3]);
    variance8x8[4] = meanOf8x8SquaredValuesBlocks[4] - (mean_of8x8_blocks[4] * mean_of8x8_blocks[4]);
    variance8x8[5] = meanOf8x8SquaredValuesBlocks[5] - (mean_of8x8_blocks[5] * mean_of8x8_blocks[5]);
    variance8x8[6] = meanOf8x8SquaredValuesBlocks[6] - (mean_of8x8_blocks[6] * mean_of8x8_blocks[6]);
    variance8x8[7] = meanOf8x8SquaredValuesBlocks[7] - (mean_of8x8_blocks[7] * mean_of8x8_blocks[7]);
    variance8x8[8] = meanOf8x8SquaredValuesBlocks[8] - (mean_of8x8_blocks[8] * mean_of8x8_blocks[8]);
    variance8x8[9] = meanOf8x8SquaredValuesBlocks[9] - (mean_of8x8_blocks[9] * mean_of8x8_blocks[9]);
    variance8x8[10] = meanOf8x8SquaredValuesBlocks[10] - (mean_of8x8_blocks[10] * mean_of8x8_blocks[10]);
    variance8x8[11] = meanOf8x8SquaredValuesBlocks[11] - (mean_of8x8_blocks[11] * mean_of8x8_blocks[11]);
    variance8x8[12] = meanOf8x8SquaredValuesBlocks[12] - (mean_of8x8_blocks[12] * mean_of8x8_blocks[12]);
    variance8x8[13] = meanOf8x8SquaredValuesBlocks[13] - (mean_of8x8_blocks[13] * mean_of8x8_blocks[13]);
    variance8x8[14] = meanOf8x8SquaredValuesBlocks[14] - (mean_of8x8_blocks[14] * mean_of8x8_blocks[14]);
    variance8x8[15] = meanOf8x8SquaredValuesBlocks[15] - (mean_of8x8_blocks[15] * mean_of8x8_blocks[15]);

    // 16x16
    meanOf16x16Blocks[0] = (mean_of8x8_blocks[0] + mean_of8x8_blocks[1] + mean_of8x8_blocks[8] + mean_of8x8_blocks[9]) >> 2;
    meanOf16x16Blocks[1] = (mean_of8x8_blocks[2] + mean_of8x8_blocks[3] + mean_of8x8_blocks[10] + mean_of8x8_blocks[11]) >> 2;
    meanOf16x16Blocks[2] = (mean_of8x8_blocks[4] + mean_of8x8_blocks[5] + mean_of8x8_blocks[12] + mean_of8x8_blocks[13]) >> 2;
    meanOf16x16Blocks[3] = (mean_of8x8_blocks[6] + mean_of8x8_blocks[7] + mean_of8x8_blocks[14] + mean_of8x8_blocks[15]) >> 2;

    meanOf16x16SquaredValuesBlocks[0] = (meanOf8x8SquaredValuesBlocks[0] + meanOf8x8SquaredValuesBlocks[1] + meanOf8x8SquaredValuesBlocks[8] + meanOf8x8SquaredValuesBlocks[9]) >> 2;
    meanOf16x16SquaredValuesBlocks[1] = (meanOf8x8SquaredValuesBlocks[2] + meanOf8x8SquaredValuesBlocks[3] + meanOf8x8SquaredValuesBlocks[10] + meanOf8x8SquaredValuesBlocks[11]) >> 2;
    meanOf16x16SquaredValuesBlocks[2] = (meanOf8x8SquaredValuesBlocks[4] + meanOf8x8SquaredValuesBlocks[5] + meanOf8x8SquaredValuesBlocks[12] + meanOf8x8SquaredValuesBlocks[13]) >> 2;
    meanOf16x16SquaredValuesBlocks[3] = (meanOf8x8SquaredValuesBlocks[6] + meanOf8x8SquaredValuesBlocks[7] + meanOf8x8SquaredValuesBlocks[14] + meanOf8x8SquaredValuesBlocks[15]) >> 2;

    // 32x32
    meanOf32x32Blocks = (meanOf16x16Blocks[0] + meanOf16x16Blocks[1] + meanOf16x16Blocks[2] + meanOf16x16Blocks[3]) >> 2;

    meanOf32x32SquaredValuesBlocks = (meanOf16x16SquaredValuesBlocks[0] + meanOf16x16SquaredValuesBlocks[1] + meanOf16x16SquaredValuesBlocks[2] + meanOf16x16SquaredValuesBlocks[3]) >> 2;

    return (meanOf32x32SquaredValuesBlocks - (meanOf32x32Blocks * meanOf32x32Blocks));
}

uint64_t ComputeVariance16x16(
    EbPictureBufferDesc       *input_padded_picture_ptr,         // input parameter, Input Padded Picture
    uint32_t                       inputLumaOriginIndex,          // input parameter, SB index, used to point to source/reference samples
    uint64_t                        *variance8x8)
{
    uint32_t blockIndex;

    uint64_t mean_of8x8_blocks[4];
    uint64_t meanOf8x8SquaredValuesBlocks[4];

    uint64_t meanOf16x16Blocks;
    uint64_t meanOf16x16SquaredValuesBlocks;

    // (0,0)
    blockIndex = inputLumaOriginIndex;

    mean_of8x8_blocks[0] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
    meanOf8x8SquaredValuesBlocks[0] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

    // (0,1)
    blockIndex = blockIndex + 8;
    mean_of8x8_blocks[1] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
    meanOf8x8SquaredValuesBlocks[1] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

    // (1,0)
    blockIndex = inputLumaOriginIndex + (input_padded_picture_ptr->stride_y << 3);
    mean_of8x8_blocks[2] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
    meanOf8x8SquaredValuesBlocks[2] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

    // (1,1)
    blockIndex = blockIndex + 8;
    mean_of8x8_blocks[3] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
    meanOf8x8SquaredValuesBlocks[3] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

    variance8x8[0] = meanOf8x8SquaredValuesBlocks[0] - (mean_of8x8_blocks[0] * mean_of8x8_blocks[0]);
    variance8x8[1] = meanOf8x8SquaredValuesBlocks[1] - (mean_of8x8_blocks[1] * mean_of8x8_blocks[1]);
    variance8x8[2] = meanOf8x8SquaredValuesBlocks[2] - (mean_of8x8_blocks[2] * mean_of8x8_blocks[2]);
    variance8x8[3] = meanOf8x8SquaredValuesBlocks[3] - (mean_of8x8_blocks[3] * mean_of8x8_blocks[3]);

    // 16x16
    meanOf16x16Blocks = (mean_of8x8_blocks[0] + mean_of8x8_blocks[1] + mean_of8x8_blocks[2] + mean_of8x8_blocks[3]) >> 2;
    meanOf16x16SquaredValuesBlocks = (meanOf8x8SquaredValuesBlocks[0] + meanOf8x8SquaredValuesBlocks[1] + meanOf8x8SquaredValuesBlocks[2] + meanOf8x8SquaredValuesBlocks[3]) >> 2;

    return (meanOf16x16SquaredValuesBlocks - (meanOf16x16Blocks * meanOf16x16Blocks));
}


/*******************************************
 * compute_mean
 *   returns the mean of a block
 *******************************************/
uint64_t compute_mean_c(
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
 * compute_mean_squared_values_c
 *   returns the Mean of Squared Values
 *******************************************/
uint64_t compute_mean_squared_values_c(
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

uint64_t compute_sub_mean_c(
    uint8_t* input_samples,     /**< input parameter, input samples Ptr */
    uint32_t input_stride,      /**< input parameter, input stride */
    uint32_t input_area_width,  /**< input parameter, input area width */
    uint32_t input_area_height) /**< input parameter, input area height */
{
    uint32_t hi, vi;
    uint64_t block_mean = 0;
    uint16_t skip = 0;

    for (vi = 0; skip < input_area_height; skip = vi + vi) {
        for (hi = 0; hi < input_area_width; hi++) {
            block_mean += input_samples[hi];
        }
        input_samples += 2 * input_stride;
        vi++;
    }

    block_mean = block_mean << 3;  // (VARIANCE_PRECISION >> 1)) /
                                   // (input_area_width * input_area_height/2)

    return block_mean;
}

uint64_t compute_sub_mean_squared_values_c(
    uint8_t* input_samples,     /**< input parameter, input samples Ptr */
    uint32_t input_stride,      /**< input parameter, input stride */
    uint32_t input_area_width,  /**< input parameter, input area width */
    uint32_t input_area_height) /**< input parameter, input area height */
{
    uint32_t hi, vi;
    uint64_t block_mean = 0;
    uint16_t skip = 0;

    for (vi = 0; skip < input_area_height; skip = vi + vi) {
        for (hi = 0; hi < input_area_width; hi++) {
            block_mean += input_samples[hi] * input_samples[hi];
        }
        input_samples += 2 * input_stride;
        vi++;
    }

    block_mean =
        block_mean
        << 11;  // VARIANCE_PRECISION) / (input_area_width * input_area_height);

    return block_mean;
}

void compute_interm_var_four8x8_c(
    uint8_t *  input_samples,
    uint16_t   input_stride,
    uint64_t * mean_of8x8_blocks,      // mean of four  8x8
    uint64_t * mean_of_squared8x8_blocks)  // meanSquared
{
    uint32_t blockIndex = 0;
    // (0,1)
    mean_of8x8_blocks[0] = compute_sub_mean_c(
        input_samples + blockIndex, input_stride, 8, 8);
    mean_of_squared8x8_blocks[0] = compute_sub_mean_squared_values_c(input_samples + blockIndex, input_stride, 8, 8);

    // (0,2)
    blockIndex = blockIndex + 8;
    mean_of8x8_blocks[1] = compute_sub_mean_c(input_samples + blockIndex, input_stride, 8, 8);
    mean_of_squared8x8_blocks[1] = compute_sub_mean_squared_values_c(input_samples + blockIndex, input_stride, 8, 8);

    // (0,3)
    blockIndex = blockIndex + 8;
    mean_of8x8_blocks[2] = compute_sub_mean_c(input_samples + blockIndex, input_stride, 8, 8);
    mean_of_squared8x8_blocks[2] = compute_sub_mean_squared_values_c(input_samples + blockIndex, input_stride, 8, 8);

    // (0,4)
    blockIndex = blockIndex + 8;
    mean_of8x8_blocks[3] = compute_sub_mean_c(input_samples + blockIndex, input_stride, 8, 8);
    mean_of_squared8x8_blocks[3] = compute_sub_mean_squared_values_c(input_samples + blockIndex, input_stride, 8, 8);
}

/*******************************************
ComputeVariance64x64
this function is exactly same as
PictureAnalysisComputeVarianceLcu excpet it
does not store data for every block,
just returns the 64x64 data point
*******************************************/
uint64_t ComputeVariance64x64(
    SequenceControlSet        *sequence_control_set_ptr,
    EbPictureBufferDesc       *input_padded_picture_ptr,         // input parameter, Input Padded Picture
    uint32_t                       inputLumaOriginIndex,          // input parameter, SB index, used to point to source/reference samples
    uint64_t                        *variance32x32)
{
    uint32_t blockIndex;

    uint64_t mean_of8x8_blocks[64];
    uint64_t meanOf8x8SquaredValuesBlocks[64];

    uint64_t meanOf16x16Blocks[16];
    uint64_t meanOf16x16SquaredValuesBlocks[16];

    uint64_t meanOf32x32Blocks[4];
    uint64_t meanOf32x32SquaredValuesBlocks[4];

    uint64_t meanOf64x64Blocks;
    uint64_t meanOf64x64SquaredValuesBlocks;

    // (0,0)
    blockIndex = inputLumaOriginIndex;
    const uint16_t stride_y = input_padded_picture_ptr->stride_y;
    if (sequence_control_set_ptr->block_mean_calc_prec == BLOCK_MEAN_PREC_FULL) {
        mean_of8x8_blocks[0] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[0] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (0,1)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[1] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[1] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (0,2)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[2] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[2] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (0,3)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[3] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[3] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (0,4)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[4] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[4] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (0,5)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[5] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[5] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (0,6)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[6] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[6] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (0,7)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[7] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[7] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (1,0)
        blockIndex = inputLumaOriginIndex + (input_padded_picture_ptr->stride_y << 3);
        mean_of8x8_blocks[8] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[8] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (1,1)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[9] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[9] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (1,2)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[10] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[10] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (1,3)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[11] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[11] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (1,4)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[12] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[12] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (1,5)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[13] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[13] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (1,6)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[14] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[14] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (1,7)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[15] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[15] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (2,0)
        blockIndex = inputLumaOriginIndex + (input_padded_picture_ptr->stride_y << 4);
        mean_of8x8_blocks[16] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[16] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (2,1)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[17] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[17] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (2,2)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[18] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[18] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (2,3)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[19] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[19] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        /// (2,4)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[20] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[20] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (2,5)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[21] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[21] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (2,6)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[22] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[22] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (2,7)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[23] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[23] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (3,0)
        blockIndex = inputLumaOriginIndex + (input_padded_picture_ptr->stride_y << 3) + (input_padded_picture_ptr->stride_y << 4);
        mean_of8x8_blocks[24] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[24] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (3,1)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[25] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[25] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (3,2)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[26] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[26] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (3,3)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[27] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[27] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (3,4)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[28] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[28] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (3,5)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[29] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[29] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (3,6)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[30] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[30] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (3,7)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[31] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[31] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (4,0)
        blockIndex = inputLumaOriginIndex + (input_padded_picture_ptr->stride_y << 5);
        mean_of8x8_blocks[32] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[32] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (4,1)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[33] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[33] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (4,2)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[34] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[34] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (4,3)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[35] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[35] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (4,4)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[36] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[36] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (4,5)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[37] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[37] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (4,6)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[38] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[38] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (4,7)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[39] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[39] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (5,0)
        blockIndex = inputLumaOriginIndex + (input_padded_picture_ptr->stride_y << 3) + (input_padded_picture_ptr->stride_y << 5);
        mean_of8x8_blocks[40] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[40] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (5,1)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[41] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[41] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (5,2)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[42] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[42] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (5,3)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[43] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[43] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (5,4)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[44] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[44] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (5,5)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[45] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[45] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (5,6)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[46] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[46] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (5,7)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[47] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[47] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (6,0)
        blockIndex = inputLumaOriginIndex + (input_padded_picture_ptr->stride_y << 4) + (input_padded_picture_ptr->stride_y << 5);
        mean_of8x8_blocks[48] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[48] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (6,1)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[49] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[49] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (6,2)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[50] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[50] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (6,3)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[51] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[51] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (6,4)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[52] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[52] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (6,5)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[53] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[53] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (6,6)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[54] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[54] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (6,7)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[55] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[55] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (7,0)
        blockIndex = inputLumaOriginIndex + (input_padded_picture_ptr->stride_y << 3) + (input_padded_picture_ptr->stride_y << 4) + (input_padded_picture_ptr->stride_y << 5);
        mean_of8x8_blocks[56] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[56] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (7,1)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[57] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[57] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (7,2)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[58] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[58] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (7,3)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[59] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[59] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (7,4)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[60] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[60] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (7,5)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[61] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[61] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (7,6)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[62] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[62] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (7,7)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[63] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[63] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
    }

    else {

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), stride_y, &mean_of8x8_blocks[0], &meanOf8x8SquaredValuesBlocks[0]);

        // (0,1)
        blockIndex = blockIndex + 32;

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), stride_y, &mean_of8x8_blocks[4], &meanOf8x8SquaredValuesBlocks[4]);
        // (0,5)
        blockIndex = blockIndex + 24;

        // (1,0)
        blockIndex = inputLumaOriginIndex + (stride_y << 3);

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), stride_y, &mean_of8x8_blocks[8], &meanOf8x8SquaredValuesBlocks[8]);

        // (1,1)
        blockIndex = blockIndex + 32;

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), stride_y, &mean_of8x8_blocks[12], &meanOf8x8SquaredValuesBlocks[12]);

        // (1,5)
        blockIndex = blockIndex + 24;

        // (2,0)
        blockIndex = inputLumaOriginIndex + (stride_y << 4);

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), stride_y, &mean_of8x8_blocks[16], &meanOf8x8SquaredValuesBlocks[16]);

        // (2,1)
        blockIndex = blockIndex + 32;

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), stride_y, &mean_of8x8_blocks[20], &meanOf8x8SquaredValuesBlocks[20]);

        // (2,5)
        blockIndex = blockIndex + 24;

        // (3,0)
        blockIndex = inputLumaOriginIndex + (stride_y << 3) + (stride_y << 4);

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), stride_y, &mean_of8x8_blocks[24], &meanOf8x8SquaredValuesBlocks[24]);

        // (3,1)
        blockIndex = blockIndex + 32;

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), stride_y, &mean_of8x8_blocks[28], &meanOf8x8SquaredValuesBlocks[28]);

        // (3,5)
        blockIndex = blockIndex + 24;

        // (4,0)
        blockIndex = inputLumaOriginIndex + (stride_y << 5);

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), stride_y, &mean_of8x8_blocks[32], &meanOf8x8SquaredValuesBlocks[32]);

        // (4,1)
        blockIndex = blockIndex + 32;

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), stride_y, &mean_of8x8_blocks[36], &meanOf8x8SquaredValuesBlocks[36]);

        // (4,5)
        blockIndex = blockIndex + 24;

        // (5,0)
        blockIndex = inputLumaOriginIndex + (stride_y << 3) + (stride_y << 5);

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), stride_y, &mean_of8x8_blocks[40], &meanOf8x8SquaredValuesBlocks[40]);

        // (5,1)
        blockIndex = blockIndex + 32;

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), stride_y, &mean_of8x8_blocks[44], &meanOf8x8SquaredValuesBlocks[44]);

        // (5,5)
        blockIndex = blockIndex + 24;

        // (6,0)
        blockIndex = inputLumaOriginIndex + (stride_y << 4) + (stride_y << 5);

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), stride_y, &mean_of8x8_blocks[48], &meanOf8x8SquaredValuesBlocks[48]);

        // (6,1)
        blockIndex = blockIndex + 32;

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), stride_y, &mean_of8x8_blocks[52], &meanOf8x8SquaredValuesBlocks[52]);

        // (6,5)
        blockIndex = blockIndex + 24;

        // (7,0)
        blockIndex = inputLumaOriginIndex + (stride_y << 3) + (stride_y << 4) + (stride_y << 5);

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), stride_y, &mean_of8x8_blocks[56], &meanOf8x8SquaredValuesBlocks[56]);

        // (7,1)
        blockIndex = blockIndex + 32;

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), stride_y, &mean_of8x8_blocks[60], &meanOf8x8SquaredValuesBlocks[60]);

    }

    // 16x16
    meanOf16x16Blocks[0] = (mean_of8x8_blocks[0] + mean_of8x8_blocks[1] + mean_of8x8_blocks[8] + mean_of8x8_blocks[9]) >> 2;
    meanOf16x16Blocks[1] = (mean_of8x8_blocks[2] + mean_of8x8_blocks[3] + mean_of8x8_blocks[10] + mean_of8x8_blocks[11]) >> 2;
    meanOf16x16Blocks[2] = (mean_of8x8_blocks[4] + mean_of8x8_blocks[5] + mean_of8x8_blocks[12] + mean_of8x8_blocks[13]) >> 2;
    meanOf16x16Blocks[3] = (mean_of8x8_blocks[6] + mean_of8x8_blocks[7] + mean_of8x8_blocks[14] + mean_of8x8_blocks[15]) >> 2;

    meanOf16x16Blocks[4] = (mean_of8x8_blocks[16] + mean_of8x8_blocks[17] + mean_of8x8_blocks[24] + mean_of8x8_blocks[25]) >> 2;
    meanOf16x16Blocks[5] = (mean_of8x8_blocks[18] + mean_of8x8_blocks[19] + mean_of8x8_blocks[26] + mean_of8x8_blocks[27]) >> 2;
    meanOf16x16Blocks[6] = (mean_of8x8_blocks[20] + mean_of8x8_blocks[21] + mean_of8x8_blocks[28] + mean_of8x8_blocks[29]) >> 2;
    meanOf16x16Blocks[7] = (mean_of8x8_blocks[22] + mean_of8x8_blocks[23] + mean_of8x8_blocks[30] + mean_of8x8_blocks[31]) >> 2;

    meanOf16x16Blocks[8] = (mean_of8x8_blocks[32] + mean_of8x8_blocks[33] + mean_of8x8_blocks[40] + mean_of8x8_blocks[41]) >> 2;
    meanOf16x16Blocks[9] = (mean_of8x8_blocks[34] + mean_of8x8_blocks[35] + mean_of8x8_blocks[42] + mean_of8x8_blocks[43]) >> 2;
    meanOf16x16Blocks[10] = (mean_of8x8_blocks[36] + mean_of8x8_blocks[37] + mean_of8x8_blocks[44] + mean_of8x8_blocks[45]) >> 2;
    meanOf16x16Blocks[11] = (mean_of8x8_blocks[38] + mean_of8x8_blocks[39] + mean_of8x8_blocks[46] + mean_of8x8_blocks[47]) >> 2;

    meanOf16x16Blocks[12] = (mean_of8x8_blocks[48] + mean_of8x8_blocks[49] + mean_of8x8_blocks[56] + mean_of8x8_blocks[57]) >> 2;
    meanOf16x16Blocks[13] = (mean_of8x8_blocks[50] + mean_of8x8_blocks[51] + mean_of8x8_blocks[58] + mean_of8x8_blocks[59]) >> 2;
    meanOf16x16Blocks[14] = (mean_of8x8_blocks[52] + mean_of8x8_blocks[53] + mean_of8x8_blocks[60] + mean_of8x8_blocks[61]) >> 2;
    meanOf16x16Blocks[15] = (mean_of8x8_blocks[54] + mean_of8x8_blocks[55] + mean_of8x8_blocks[62] + mean_of8x8_blocks[63]) >> 2;

    meanOf16x16SquaredValuesBlocks[0] = (meanOf8x8SquaredValuesBlocks[0] + meanOf8x8SquaredValuesBlocks[1] + meanOf8x8SquaredValuesBlocks[8] + meanOf8x8SquaredValuesBlocks[9]) >> 2;
    meanOf16x16SquaredValuesBlocks[1] = (meanOf8x8SquaredValuesBlocks[2] + meanOf8x8SquaredValuesBlocks[3] + meanOf8x8SquaredValuesBlocks[10] + meanOf8x8SquaredValuesBlocks[11]) >> 2;
    meanOf16x16SquaredValuesBlocks[2] = (meanOf8x8SquaredValuesBlocks[4] + meanOf8x8SquaredValuesBlocks[5] + meanOf8x8SquaredValuesBlocks[12] + meanOf8x8SquaredValuesBlocks[13]) >> 2;
    meanOf16x16SquaredValuesBlocks[3] = (meanOf8x8SquaredValuesBlocks[6] + meanOf8x8SquaredValuesBlocks[7] + meanOf8x8SquaredValuesBlocks[14] + meanOf8x8SquaredValuesBlocks[15]) >> 2;

    meanOf16x16SquaredValuesBlocks[4] = (meanOf8x8SquaredValuesBlocks[16] + meanOf8x8SquaredValuesBlocks[17] + meanOf8x8SquaredValuesBlocks[24] + meanOf8x8SquaredValuesBlocks[25]) >> 2;
    meanOf16x16SquaredValuesBlocks[5] = (meanOf8x8SquaredValuesBlocks[18] + meanOf8x8SquaredValuesBlocks[19] + meanOf8x8SquaredValuesBlocks[26] + meanOf8x8SquaredValuesBlocks[27]) >> 2;
    meanOf16x16SquaredValuesBlocks[6] = (meanOf8x8SquaredValuesBlocks[20] + meanOf8x8SquaredValuesBlocks[21] + meanOf8x8SquaredValuesBlocks[28] + meanOf8x8SquaredValuesBlocks[29]) >> 2;
    meanOf16x16SquaredValuesBlocks[7] = (meanOf8x8SquaredValuesBlocks[22] + meanOf8x8SquaredValuesBlocks[23] + meanOf8x8SquaredValuesBlocks[30] + meanOf8x8SquaredValuesBlocks[31]) >> 2;

    meanOf16x16SquaredValuesBlocks[8] = (meanOf8x8SquaredValuesBlocks[32] + meanOf8x8SquaredValuesBlocks[33] + meanOf8x8SquaredValuesBlocks[40] + meanOf8x8SquaredValuesBlocks[41]) >> 2;
    meanOf16x16SquaredValuesBlocks[9] = (meanOf8x8SquaredValuesBlocks[34] + meanOf8x8SquaredValuesBlocks[35] + meanOf8x8SquaredValuesBlocks[42] + meanOf8x8SquaredValuesBlocks[43]) >> 2;
    meanOf16x16SquaredValuesBlocks[10] = (meanOf8x8SquaredValuesBlocks[36] + meanOf8x8SquaredValuesBlocks[37] + meanOf8x8SquaredValuesBlocks[44] + meanOf8x8SquaredValuesBlocks[45]) >> 2;
    meanOf16x16SquaredValuesBlocks[11] = (meanOf8x8SquaredValuesBlocks[38] + meanOf8x8SquaredValuesBlocks[39] + meanOf8x8SquaredValuesBlocks[46] + meanOf8x8SquaredValuesBlocks[47]) >> 2;

    meanOf16x16SquaredValuesBlocks[12] = (meanOf8x8SquaredValuesBlocks[48] + meanOf8x8SquaredValuesBlocks[49] + meanOf8x8SquaredValuesBlocks[56] + meanOf8x8SquaredValuesBlocks[57]) >> 2;
    meanOf16x16SquaredValuesBlocks[13] = (meanOf8x8SquaredValuesBlocks[50] + meanOf8x8SquaredValuesBlocks[51] + meanOf8x8SquaredValuesBlocks[58] + meanOf8x8SquaredValuesBlocks[59]) >> 2;
    meanOf16x16SquaredValuesBlocks[14] = (meanOf8x8SquaredValuesBlocks[52] + meanOf8x8SquaredValuesBlocks[53] + meanOf8x8SquaredValuesBlocks[60] + meanOf8x8SquaredValuesBlocks[61]) >> 2;
    meanOf16x16SquaredValuesBlocks[15] = (meanOf8x8SquaredValuesBlocks[54] + meanOf8x8SquaredValuesBlocks[55] + meanOf8x8SquaredValuesBlocks[62] + meanOf8x8SquaredValuesBlocks[63]) >> 2;

    // 32x32
    meanOf32x32Blocks[0] = (meanOf16x16Blocks[0] + meanOf16x16Blocks[1] + meanOf16x16Blocks[4] + meanOf16x16Blocks[5]) >> 2;
    meanOf32x32Blocks[1] = (meanOf16x16Blocks[2] + meanOf16x16Blocks[3] + meanOf16x16Blocks[6] + meanOf16x16Blocks[7]) >> 2;
    meanOf32x32Blocks[2] = (meanOf16x16Blocks[8] + meanOf16x16Blocks[9] + meanOf16x16Blocks[12] + meanOf16x16Blocks[13]) >> 2;
    meanOf32x32Blocks[3] = (meanOf16x16Blocks[10] + meanOf16x16Blocks[11] + meanOf16x16Blocks[14] + meanOf16x16Blocks[15]) >> 2;

    meanOf32x32SquaredValuesBlocks[0] = (meanOf16x16SquaredValuesBlocks[0] + meanOf16x16SquaredValuesBlocks[1] + meanOf16x16SquaredValuesBlocks[4] + meanOf16x16SquaredValuesBlocks[5]) >> 2;
    meanOf32x32SquaredValuesBlocks[1] = (meanOf16x16SquaredValuesBlocks[2] + meanOf16x16SquaredValuesBlocks[3] + meanOf16x16SquaredValuesBlocks[6] + meanOf16x16SquaredValuesBlocks[7]) >> 2;
    meanOf32x32SquaredValuesBlocks[2] = (meanOf16x16SquaredValuesBlocks[8] + meanOf16x16SquaredValuesBlocks[9] + meanOf16x16SquaredValuesBlocks[12] + meanOf16x16SquaredValuesBlocks[13]) >> 2;
    meanOf32x32SquaredValuesBlocks[3] = (meanOf16x16SquaredValuesBlocks[10] + meanOf16x16SquaredValuesBlocks[11] + meanOf16x16SquaredValuesBlocks[14] + meanOf16x16SquaredValuesBlocks[15]) >> 2;

    variance32x32[0] = meanOf32x32SquaredValuesBlocks[0] - (meanOf32x32Blocks[0] * meanOf32x32Blocks[0]);
    variance32x32[1] = meanOf32x32SquaredValuesBlocks[1] - (meanOf32x32Blocks[1] * meanOf32x32Blocks[1]);
    variance32x32[2] = meanOf32x32SquaredValuesBlocks[2] - (meanOf32x32Blocks[2] * meanOf32x32Blocks[2]);
    variance32x32[3] = meanOf32x32SquaredValuesBlocks[3] - (meanOf32x32Blocks[3] * meanOf32x32Blocks[3]);

    // 64x64
    meanOf64x64Blocks = (meanOf32x32Blocks[0] + meanOf32x32Blocks[1] + meanOf32x32Blocks[2] + meanOf32x32Blocks[3]) >> 2;
    meanOf64x64SquaredValuesBlocks = (meanOf32x32SquaredValuesBlocks[0] + meanOf32x32SquaredValuesBlocks[1] + meanOf32x32SquaredValuesBlocks[2] + meanOf32x32SquaredValuesBlocks[3]) >> 2;

    return (meanOf64x64SquaredValuesBlocks - (meanOf64x64Blocks * meanOf64x64Blocks));
}

uint8_t  getFilteredTypes(uint8_t  *ptr,
    uint32_t  stride,
    uint8_t   filterType)
{
    uint8_t *p = ptr - 1 - stride;

    uint32_t a = 0;

    if (filterType == 0) {
        //Luma
        a = (p[1] +
            p[0 + stride] + 4 * p[1 + stride] + p[2 + stride] +
            p[1 + 2 * stride]) / 8;
    }
    else if (filterType == 1) {
        a = (2 * p[1] +
            2 * p[0 + stride] + 4 * p[1 + stride] + 2 * p[2 + stride] +
            2 * p[1 + 2 * stride]);

        a = (((uint32_t)((a * 2730) >> 14) + 1) >> 1) & 0xFFFF;

        //fixed point version of a=a/12 to mimic x86 instruction _mm256_mulhrs_epi16;
        //a= (a*2730)>>15;
    }
    else if (filterType == 2) {
        a = (4 * p[1] +
            4 * p[0 + stride] + 4 * p[1 + stride] + 4 * p[2 + stride] +
            4 * p[1 + 2 * stride]) / 20;
    }
    else if (filterType == 3) {
        a = (1 * p[0] + 1 * p[1] + 1 * p[2] +
            1 * p[0 + stride] + 4 * p[1 + stride] + 1 * p[2 + stride] +
            1 * p[0 + 2 * stride] + 1 * p[1 + 2 * stride] + 1 * p[2 + 2 * stride]) / 12;
    }
    else if (filterType == 4) {
        //gaussian matrix(Chroma)
        a = (1 * p[0] + 2 * p[1] + 1 * p[2] +
            2 * p[0 + stride] + 4 * p[1 + stride] + 2 * p[2 + stride] +
            1 * p[0 + 2 * stride] + 2 * p[1 + 2 * stride] + 1 * p[2 + 2 * stride]) / 16;
    }
    else if (filterType == 5) {
        a = (2 * p[0] + 2 * p[1] + 2 * p[2] +
            2 * p[0 + stride] + 4 * p[1 + stride] + 2 * p[2 + stride] +
            2 * p[0 + 2 * stride] + 2 * p[1 + 2 * stride] + 2 * p[2 + 2 * stride]) / 20;
    }
    else if (filterType == 6) {
        a = (4 * p[0] + 4 * p[1] + 4 * p[2] +
            4 * p[0 + stride] + 4 * p[1 + stride] + 4 * p[2 + stride] +
            4 * p[0 + 2 * stride] + 4 * p[1 + 2 * stride] + 4 * p[2 + 2 * stride]) / 36;
    }

    return  (uint8_t)CLIP3EQ(0, 255, a);
}

/*******************************************
* noise_extract_luma_strong
*  strong filter Luma.
*******************************************/
void noise_extract_luma_strong_c(
    EbPictureBufferDesc       *input_picture_ptr,
    EbPictureBufferDesc       *denoised_picture_ptr,
    uint32_t                       sb_origin_y
    , uint32_t                       sb_origin_x
)
{
    uint32_t  ii, jj;
    uint32_t  picHeight, sb_height;
    uint32_t  picWidth;
    uint32_t  inputOriginIndex;
    uint32_t  inputOriginIndexPad;

    uint8_t *ptrIn;
    uint32_t stride_in;
    uint8_t *ptr_denoised;

    uint32_t strideOut;
    uint32_t idx = (sb_origin_x + BLOCK_SIZE_64 > input_picture_ptr->width) ? sb_origin_x : 0;

    //Luma
    {
        picHeight = input_picture_ptr->height;
        picWidth = input_picture_ptr->width;
        sb_height = MIN(BLOCK_SIZE_64, picHeight - sb_origin_y);

        stride_in = input_picture_ptr->stride_y;
        inputOriginIndex = input_picture_ptr->origin_x + (input_picture_ptr->origin_y + sb_origin_y)* input_picture_ptr->stride_y;
        ptrIn = &(input_picture_ptr->buffer_y[inputOriginIndex]);

        inputOriginIndexPad = denoised_picture_ptr->origin_x + (denoised_picture_ptr->origin_y + sb_origin_y) * denoised_picture_ptr->stride_y;
        strideOut = denoised_picture_ptr->stride_y;
        ptr_denoised = &(denoised_picture_ptr->buffer_y[inputOriginIndexPad]);

        for (jj = 0; jj < sb_height; jj++) {
            for (ii = idx; ii < picWidth; ii++) {
                if ((jj > 0 || sb_origin_y > 0) && (jj < sb_height - 1 || sb_origin_y + sb_height < picHeight) && ii > 0 && ii < picWidth - 1)
                    ptr_denoised[ii + jj * strideOut] = getFilteredTypes(&ptrIn[ii + jj * stride_in], stride_in, 4);
                else
                    ptr_denoised[ii + jj * strideOut] = ptrIn[ii + jj * stride_in];
            }
        }
    }
}
/*******************************************
* noise_extract_chroma_strong
*  strong filter chroma.
*******************************************/
void noise_extract_chroma_strong_c(
    EbPictureBufferDesc       *input_picture_ptr,
    EbPictureBufferDesc       *denoised_picture_ptr,
    uint32_t                       sb_origin_y
    , uint32_t                       sb_origin_x
)
{
    uint32_t  ii, jj;
    uint32_t  picHeight, sb_height;
    uint32_t  picWidth;
    uint32_t  inputOriginIndex;
    uint32_t  inputOriginIndexPad;

    uint8_t *ptrIn;
    uint32_t stride_in;
    uint8_t *ptr_denoised;

    uint32_t strideOut;
    uint32_t idx = (sb_origin_x + BLOCK_SIZE_64 > input_picture_ptr->width) ? sb_origin_x : 0;

    uint32_t color_format = input_picture_ptr->color_format;
    const uint16_t subsampling_x = (color_format == EB_YUV444 ? 1 : 2) - 1;
    const uint16_t subsampling_y = (color_format >= EB_YUV422 ? 1 : 2) - 1;

    //Cb
    {
        picHeight = input_picture_ptr->height >> subsampling_y;
        picWidth = input_picture_ptr->width >> subsampling_x;
        sb_height = MIN(BLOCK_SIZE_64 >> subsampling_y, picHeight - sb_origin_y);

        stride_in = input_picture_ptr->stride_cb;
        inputOriginIndex = (input_picture_ptr->origin_x >> subsampling_x) + ((input_picture_ptr->origin_y >> subsampling_y) + sb_origin_y) * input_picture_ptr->stride_cb;
        ptrIn = &(input_picture_ptr->buffer_cb[inputOriginIndex]);

        inputOriginIndexPad = (denoised_picture_ptr->origin_x >> subsampling_x) + ((denoised_picture_ptr->origin_y >> subsampling_y) + sb_origin_y) * denoised_picture_ptr->stride_cb;
        strideOut = denoised_picture_ptr->stride_cb;
        ptr_denoised = &(denoised_picture_ptr->buffer_cb[inputOriginIndexPad]);

        for (jj = 0; jj < sb_height; jj++) {
            for (ii = idx; ii < picWidth; ii++) {
                if ((jj > 0 || sb_origin_y > 0) && (jj < sb_height - 1 || (sb_origin_y + sb_height) < picHeight) && ii > 0 && ii < picWidth - 1)
                    ptr_denoised[ii + jj * strideOut] = getFilteredTypes(&ptrIn[ii + jj * stride_in], stride_in, 6);
                else
                    ptr_denoised[ii + jj * strideOut] = ptrIn[ii + jj * stride_in];
            }
        }
    }

    //Cr
    {
        picHeight = input_picture_ptr->height >> subsampling_y;
        picWidth = input_picture_ptr->width >> subsampling_x;
        sb_height = MIN(BLOCK_SIZE_64 >> subsampling_y, picHeight - sb_origin_y);

        stride_in = input_picture_ptr->stride_cr;
        inputOriginIndex = (input_picture_ptr->origin_x >> subsampling_x) + ((input_picture_ptr->origin_y >> subsampling_y) + sb_origin_y) * input_picture_ptr->stride_cr;

        ptrIn = &(input_picture_ptr->buffer_cr[inputOriginIndex]);

        inputOriginIndexPad = (denoised_picture_ptr->origin_x >> subsampling_x) + ((denoised_picture_ptr->origin_y >> subsampling_y) + sb_origin_y) * denoised_picture_ptr->stride_cr;
        strideOut = denoised_picture_ptr->stride_cr;
        ptr_denoised = &(denoised_picture_ptr->buffer_cr[inputOriginIndexPad]);

        for (jj = 0; jj < sb_height; jj++) {
            for (ii = idx; ii < picWidth; ii++) {
                if ((jj > 0 || sb_origin_y > 0) && (jj < sb_height - 1 || (sb_origin_y + sb_height) < picHeight) && ii > 0 && ii < picWidth - 1)
                    ptr_denoised[ii + jj * strideOut] = getFilteredTypes(&ptrIn[ii + jj * stride_in], stride_in, 6);
                else
                    ptr_denoised[ii + jj * strideOut] = ptrIn[ii + jj * stride_in];
            }
        }
    }
}

/*******************************************
* noise_extract_chroma_weak
*  weak filter chroma.
*******************************************/
void noise_extract_chroma_weak_c(
    EbPictureBufferDesc       *input_picture_ptr,
    EbPictureBufferDesc       *denoised_picture_ptr,
    uint32_t                       sb_origin_y
    , uint32_t                       sb_origin_x
)
{
    uint32_t  ii, jj;
    uint32_t  picHeight, sb_height;
    uint32_t  picWidth;
    uint32_t  inputOriginIndex;
    uint32_t  inputOriginIndexPad;

    uint8_t *ptrIn;
    uint32_t stride_in;
    uint8_t *ptr_denoised;

    uint32_t strideOut;

    uint32_t idx = (sb_origin_x + BLOCK_SIZE_64 > input_picture_ptr->width) ? sb_origin_x : 0;

    uint32_t color_format = input_picture_ptr->color_format;
    const uint16_t subsampling_x = (color_format == EB_YUV444 ? 1 : 2) - 1;
    const uint16_t subsampling_y = (color_format >= EB_YUV422 ? 1 : 2) - 1;

    //Cb
    {
        picHeight = input_picture_ptr->height >> subsampling_y;
        picWidth = input_picture_ptr->width >> subsampling_x;
        sb_height = MIN(BLOCK_SIZE_64 >> subsampling_y, picHeight - sb_origin_y);

        stride_in = input_picture_ptr->stride_cb;
        inputOriginIndex = (input_picture_ptr->origin_x >> subsampling_x) + ((input_picture_ptr->origin_y >> subsampling_y) + sb_origin_y) * input_picture_ptr->stride_cb;

        ptrIn = &(input_picture_ptr->buffer_cb[inputOriginIndex]);

        inputOriginIndexPad = (denoised_picture_ptr->origin_x >> subsampling_x) + ((denoised_picture_ptr->origin_y >> subsampling_y) + sb_origin_y) * denoised_picture_ptr->stride_cb;

        strideOut = denoised_picture_ptr->stride_cb;
        ptr_denoised = &(denoised_picture_ptr->buffer_cb[inputOriginIndexPad]);

        for (jj = 0; jj < sb_height; jj++) {
            for (ii = idx; ii < picWidth; ii++) {
                if ((jj > 0 || sb_origin_y > 0) && (jj < sb_height - 1 || (sb_origin_y + sb_height) < picHeight) && ii > 0 && ii < picWidth - 1)
                    ptr_denoised[ii + jj * strideOut] = getFilteredTypes(&ptrIn[ii + jj * stride_in], stride_in, 4);
                else
                    ptr_denoised[ii + jj * strideOut] = ptrIn[ii + jj * stride_in];
            }
        }
    }

    //Cr
    {
        picHeight = input_picture_ptr->height >> subsampling_y;
        picWidth = input_picture_ptr->width >> subsampling_x;
        sb_height = MIN(BLOCK_SIZE_64 >> subsampling_y, picHeight - sb_origin_y);

        stride_in = input_picture_ptr->stride_cr;
        inputOriginIndex = (input_picture_ptr->origin_x >> subsampling_x) + ((input_picture_ptr->origin_y >> subsampling_y) + sb_origin_y) * input_picture_ptr->stride_cr;
        ptrIn = &(input_picture_ptr->buffer_cr[inputOriginIndex]);

        inputOriginIndexPad = (denoised_picture_ptr->origin_x >> subsampling_x) + ((denoised_picture_ptr->origin_y >> subsampling_y) + sb_origin_y) * denoised_picture_ptr->stride_cr;
        strideOut = denoised_picture_ptr->stride_cr;
        ptr_denoised = &(denoised_picture_ptr->buffer_cr[inputOriginIndexPad]);

        for (jj = 0; jj < sb_height; jj++) {
            for (ii = idx; ii < picWidth; ii++) {
                if ((jj > 0 || sb_origin_y > 0) && (jj < sb_height - 1 || (sb_origin_y + sb_height) < picHeight) && ii > 0 && ii < picWidth - 1)
                    ptr_denoised[ii + jj * strideOut] = getFilteredTypes(&ptrIn[ii + jj * stride_in], stride_in, 4);
                else
                    ptr_denoised[ii + jj * strideOut] = ptrIn[ii + jj * stride_in];
            }
        }
    }
}

/*******************************************
* noise_extract_luma_weak
*  weak filter Luma and store noise.
*******************************************/
void noise_extract_luma_weak_c(
    EbPictureBufferDesc       *input_picture_ptr,
    EbPictureBufferDesc       *denoised_picture_ptr,
    EbPictureBufferDesc       *noise_picture_ptr,
    uint32_t                       sb_origin_y
    , uint32_t                         sb_origin_x
)
{
    uint32_t  ii, jj;
    uint32_t  picHeight, sb_height;
    uint32_t  picWidth;
    uint32_t  inputOriginIndex;
    uint32_t  inputOriginIndexPad;
    uint32_t  noiseOriginIndex;

    uint8_t *ptrIn;
    uint32_t stride_in;
    uint8_t *ptr_denoised;

    uint8_t *ptr_noise;
    uint32_t strideOut;

    uint32_t idx = (sb_origin_x + BLOCK_SIZE_64 > input_picture_ptr->width) ? sb_origin_x : 0;

    //Luma
    {
        picHeight = input_picture_ptr->height;
        picWidth = input_picture_ptr->width;
        sb_height = MIN(BLOCK_SIZE_64, picHeight - sb_origin_y);

        stride_in = input_picture_ptr->stride_y;
        inputOriginIndex = input_picture_ptr->origin_x + (input_picture_ptr->origin_y + sb_origin_y) * input_picture_ptr->stride_y;
        ptrIn = &(input_picture_ptr->buffer_y[inputOriginIndex]);

        inputOriginIndexPad = denoised_picture_ptr->origin_x + (denoised_picture_ptr->origin_y + sb_origin_y) * denoised_picture_ptr->stride_y;
        strideOut = denoised_picture_ptr->stride_y;
        ptr_denoised = &(denoised_picture_ptr->buffer_y[inputOriginIndexPad]);

        noiseOriginIndex = noise_picture_ptr->origin_x + noise_picture_ptr->origin_y * noise_picture_ptr->stride_y;
        ptr_noise = &(noise_picture_ptr->buffer_y[noiseOriginIndex]);

        for (jj = 0; jj < sb_height; jj++) {
            for (ii = idx; ii < picWidth; ii++) {
                if ((jj > 0 || sb_origin_y > 0) && (jj < sb_height - 1 || sb_origin_y + sb_height < picHeight) && ii > 0 && ii < picWidth - 1) {
                    ptr_denoised[ii + jj * strideOut] = getFilteredTypes(&ptrIn[ii + jj * stride_in], stride_in, 0);
                    ptr_noise[ii + jj * strideOut] = CLIP3EQ(0, 255, ptrIn[ii + jj * stride_in] - ptr_denoised[ii + jj * strideOut]);
                }
                else {
                    ptr_denoised[ii + jj * strideOut] = ptrIn[ii + jj * stride_in];
                    ptr_noise[ii + jj * strideOut] = 0;
                }
            }
        }
    }
}

void noise_extract_luma_weak_lcu_c(
    EbPictureBufferDesc       *input_picture_ptr,
    EbPictureBufferDesc       *denoised_picture_ptr,
    EbPictureBufferDesc       *noise_picture_ptr,
    uint32_t                       sb_origin_y
    , uint32_t                         sb_origin_x
)
{
    uint32_t  ii, jj;
    uint32_t  picHeight, sb_height;
    uint32_t  picWidth, sb_width;
    uint32_t  inputOriginIndex;
    uint32_t  inputOriginIndexPad;
    uint32_t  noiseOriginIndex;

    uint8_t *ptrIn;
    uint32_t stride_in;
    uint8_t *ptr_denoised;

    uint8_t *ptr_noise;
    uint32_t strideOut;

    uint32_t idx = (sb_origin_x + BLOCK_SIZE_64 > input_picture_ptr->width) ? sb_origin_x : 0;

    //Luma
    {
        picHeight = input_picture_ptr->height;
        picWidth = input_picture_ptr->width;
        sb_height = MIN(BLOCK_SIZE_64, picHeight - sb_origin_y);
        sb_width = MIN(BLOCK_SIZE_64, picWidth - sb_origin_x);

        stride_in = input_picture_ptr->stride_y;
        inputOriginIndex = input_picture_ptr->origin_x + sb_origin_x + (input_picture_ptr->origin_y + sb_origin_y) * input_picture_ptr->stride_y;
        ptrIn = &(input_picture_ptr->buffer_y[inputOriginIndex]);

        inputOriginIndexPad = denoised_picture_ptr->origin_x + sb_origin_x + (denoised_picture_ptr->origin_y + sb_origin_y) * denoised_picture_ptr->stride_y;
        strideOut = denoised_picture_ptr->stride_y;
        ptr_denoised = &(denoised_picture_ptr->buffer_y[inputOriginIndexPad]);

        noiseOriginIndex = noise_picture_ptr->origin_x + sb_origin_x + noise_picture_ptr->origin_y * noise_picture_ptr->stride_y;
        ptr_noise = &(noise_picture_ptr->buffer_y[noiseOriginIndex]);

        for (jj = 0; jj < sb_height; jj++) {
            for (ii = idx; ii < sb_width; ii++) {
                if ((jj > 0 || sb_origin_y > 0) && (jj < sb_height - 1 || sb_origin_y + sb_height < picHeight) && (ii > 0 || sb_origin_x > 0) && (ii + sb_origin_x) < picWidth - 1/* & ii < sb_width - 1*/) {
                    ptr_denoised[ii + jj * strideOut] = getFilteredTypes(&ptrIn[ii + jj * stride_in], stride_in, 0);
                    ptr_noise[ii + jj * strideOut] = CLIP3EQ(0, 255, ptrIn[ii + jj * stride_in] - ptr_denoised[ii + jj * strideOut]);
                }
                else {
                    ptr_denoised[ii + jj * strideOut] = ptrIn[ii + jj * stride_in];
                    ptr_noise[ii + jj * strideOut] = 0;
                }
            }
        }
    }
}

EbErrorType ZeroOutChromaBlockMean(
    PictureParentControlSet   *picture_control_set_ptr,          // input parameter, Picture Control Set Ptr
    uint32_t                       lcuCodingOrder                // input parameter, SB address
)
{
    EbErrorType return_error = EB_ErrorNone;
    // 16x16 mean
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_0] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_1] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_2] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_3] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_4] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_5] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_6] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_7] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_8] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_9] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_10] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_11] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_12] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_13] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_14] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_15] = 0;

    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_0] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_1] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_2] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_3] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_4] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_5] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_6] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_7] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_8] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_9] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_10] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_11] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_12] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_13] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_14] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_15] = 0;

    // 32x32 mean
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_32x32_0] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_32x32_1] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_32x32_2] = 0;
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_32x32_3] = 0;

    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_32x32_0] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_32x32_1] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_32x32_2] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_32x32_3] = 0;

    // 64x64 mean
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_64x64] = 0;
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_64x64] = 0;

    return return_error;
}
/*******************************************
* ComputeChromaBlockMean
*   computes the chroma block mean for 64x64, 32x32 and 16x16 CUs inside the tree block
*******************************************/
EbErrorType ComputeChromaBlockMean(
    SequenceControlSet        *sequence_control_set_ptr,
    PictureParentControlSet   *picture_control_set_ptr,          // input parameter, Picture Control Set Ptr
    EbPictureBufferDesc       *input_padded_picture_ptr,         // input parameter, Input Padded Picture
    uint32_t                       lcuCodingOrder,                // input parameter, SB address
    uint32_t                       inputCbOriginIndex,            // input parameter, SB index, used to point to source/reference samples
    uint32_t                       inputCrOriginIndex)            // input parameter, SB index, used to point to source/reference samples
{
    EbErrorType return_error = EB_ErrorNone;

    uint32_t cbBlockIndex, crBlockIndex;

    uint64_t cbMeanOf16x16Blocks[16];
    uint64_t crMeanOf16x16Blocks[16];

    uint64_t cbMeanOf32x32Blocks[4];
    uint64_t crMeanOf32x32Blocks[4];

    uint64_t cbMeanOf64x64Blocks;
    uint64_t crMeanOf64x64Blocks;

    // (0,0) 16x16 block
    cbBlockIndex = inputCbOriginIndex;
    crBlockIndex = inputCrOriginIndex;
    if (sequence_control_set_ptr->block_mean_calc_prec == BLOCK_MEAN_PREC_FULL) {
        cbMeanOf16x16Blocks[0] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_cb[cbBlockIndex]), input_padded_picture_ptr->stride_cb, 8, 8);
        crMeanOf16x16Blocks[0] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_cr[crBlockIndex]), input_padded_picture_ptr->stride_cr, 8, 8);

        // (0,1)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[1] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_cb[cbBlockIndex]), input_padded_picture_ptr->stride_cb, 8, 8);
        crMeanOf16x16Blocks[1] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_cr[crBlockIndex]), input_padded_picture_ptr->stride_cr, 8, 8);

        // (0,2)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[2] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_cb[cbBlockIndex]), input_padded_picture_ptr->stride_cb, 8, 8);
        crMeanOf16x16Blocks[2] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_cr[crBlockIndex]), input_padded_picture_ptr->stride_cr, 8, 8);

        // (0,3)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[3] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_cb[cbBlockIndex]), input_padded_picture_ptr->stride_cb, 8, 8);
        crMeanOf16x16Blocks[3] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_cr[crBlockIndex]), input_padded_picture_ptr->stride_cr, 8, 8);

        // (1,0)
        cbBlockIndex = inputCbOriginIndex + (input_padded_picture_ptr->stride_cb << 3);
        crBlockIndex = inputCrOriginIndex + (input_padded_picture_ptr->stride_cr << 3);
        cbMeanOf16x16Blocks[4] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_cb[cbBlockIndex]), input_padded_picture_ptr->stride_cb, 8, 8);
        crMeanOf16x16Blocks[4] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_cr[crBlockIndex]), input_padded_picture_ptr->stride_cr, 8, 8);

        // (1,1)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[5] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_cb[cbBlockIndex]), input_padded_picture_ptr->stride_cb, 8, 8);
        crMeanOf16x16Blocks[5] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_cr[crBlockIndex]), input_padded_picture_ptr->stride_cr, 8, 8);

        // (1,2)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[6] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_cb[cbBlockIndex]), input_padded_picture_ptr->stride_cb, 8, 8);
        crMeanOf16x16Blocks[6] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_cr[crBlockIndex]), input_padded_picture_ptr->stride_cr, 8, 8);

        // (1,3)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[7] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_cb[cbBlockIndex]), input_padded_picture_ptr->stride_cb, 8, 8);
        crMeanOf16x16Blocks[7] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_cr[crBlockIndex]), input_padded_picture_ptr->stride_cr, 8, 8);

        // (2,0)
        cbBlockIndex = inputCbOriginIndex + (input_padded_picture_ptr->stride_cb << 4);
        crBlockIndex = inputCrOriginIndex + (input_padded_picture_ptr->stride_cr << 4);
        cbMeanOf16x16Blocks[8] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_cb[cbBlockIndex]), input_padded_picture_ptr->stride_cb, 8, 8);
        crMeanOf16x16Blocks[8] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_cr[crBlockIndex]), input_padded_picture_ptr->stride_cr, 8, 8);

        // (2,1)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[9] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_cb[cbBlockIndex]), input_padded_picture_ptr->stride_cb, 8, 8);
        crMeanOf16x16Blocks[9] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_cr[crBlockIndex]), input_padded_picture_ptr->stride_cr, 8, 8);

        // (2,2)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[10] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_cb[cbBlockIndex]), input_padded_picture_ptr->stride_cb, 8, 8);
        crMeanOf16x16Blocks[10] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_cr[crBlockIndex]), input_padded_picture_ptr->stride_cr, 8, 8);

        // (2,3)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[11] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_cb[cbBlockIndex]), input_padded_picture_ptr->stride_cb, 8, 8);
        crMeanOf16x16Blocks[11] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_cr[crBlockIndex]), input_padded_picture_ptr->stride_cr, 8, 8);

        // (3,0)
        cbBlockIndex = inputCbOriginIndex + (input_padded_picture_ptr->stride_cb * 24);
        crBlockIndex = inputCrOriginIndex + (input_padded_picture_ptr->stride_cr * 24);
        cbMeanOf16x16Blocks[12] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_cb[cbBlockIndex]), input_padded_picture_ptr->stride_cb, 8, 8);
        crMeanOf16x16Blocks[12] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_cr[crBlockIndex]), input_padded_picture_ptr->stride_cr, 8, 8);

        // (3,1)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[13] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_cb[cbBlockIndex]), input_padded_picture_ptr->stride_cb, 8, 8);
        crMeanOf16x16Blocks[13] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_cr[crBlockIndex]), input_padded_picture_ptr->stride_cr, 8, 8);

        // (3,2)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[14] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_cb[cbBlockIndex]), input_padded_picture_ptr->stride_cb, 8, 8);
        crMeanOf16x16Blocks[14] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_cr[crBlockIndex]), input_padded_picture_ptr->stride_cr, 8, 8);

        // (3,3)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[15] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_cb[cbBlockIndex]), input_padded_picture_ptr->stride_cb, 8, 8);
        crMeanOf16x16Blocks[15] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_cr[crBlockIndex]), input_padded_picture_ptr->stride_cr, 8, 8);
    }
    else {
        const uint16_t stride_cb = input_padded_picture_ptr->stride_cb;
        const uint16_t stride_cr = input_padded_picture_ptr->stride_cr;

        cbMeanOf16x16Blocks[0] = compute_sub_mean8x8_sse2_intrin(&(input_padded_picture_ptr->buffer_cb[cbBlockIndex]), stride_cb);
        crMeanOf16x16Blocks[0] = compute_sub_mean8x8_sse2_intrin(&(input_padded_picture_ptr->buffer_cr[crBlockIndex]), stride_cr);

        // (0,1)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[1] = compute_sub_mean8x8_sse2_intrin(&(input_padded_picture_ptr->buffer_cb[cbBlockIndex]), stride_cb);
        crMeanOf16x16Blocks[1] = compute_sub_mean8x8_sse2_intrin(&(input_padded_picture_ptr->buffer_cr[crBlockIndex]), stride_cr);

        // (0,2)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[2] = compute_sub_mean8x8_sse2_intrin(&(input_padded_picture_ptr->buffer_cb[cbBlockIndex]), stride_cb);
        crMeanOf16x16Blocks[2] = compute_sub_mean8x8_sse2_intrin(&(input_padded_picture_ptr->buffer_cr[crBlockIndex]), stride_cr);

        // (0,3)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[3] = compute_sub_mean8x8_sse2_intrin(&(input_padded_picture_ptr->buffer_cb[cbBlockIndex]), stride_cb);
        crMeanOf16x16Blocks[3] = compute_sub_mean8x8_sse2_intrin(&(input_padded_picture_ptr->buffer_cr[crBlockIndex]), stride_cr);

        // (1,0)
        cbBlockIndex = inputCbOriginIndex + (stride_cb << 3);
        crBlockIndex = inputCrOriginIndex + (stride_cr << 3);
        cbMeanOf16x16Blocks[4] = compute_sub_mean8x8_sse2_intrin(&(input_padded_picture_ptr->buffer_cb[cbBlockIndex]), stride_cb);
        crMeanOf16x16Blocks[4] = compute_sub_mean8x8_sse2_intrin(&(input_padded_picture_ptr->buffer_cr[crBlockIndex]), stride_cr);

        // (1,1)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[5] = compute_sub_mean8x8_sse2_intrin(&(input_padded_picture_ptr->buffer_cb[cbBlockIndex]), stride_cb);
        crMeanOf16x16Blocks[5] = compute_sub_mean8x8_sse2_intrin(&(input_padded_picture_ptr->buffer_cr[crBlockIndex]), stride_cr);

        // (1,2)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[6] = compute_sub_mean8x8_sse2_intrin(&(input_padded_picture_ptr->buffer_cb[cbBlockIndex]), stride_cb);
        crMeanOf16x16Blocks[6] = compute_sub_mean8x8_sse2_intrin(&(input_padded_picture_ptr->buffer_cr[crBlockIndex]), stride_cr);

        // (1,3)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[7] = compute_sub_mean8x8_sse2_intrin(&(input_padded_picture_ptr->buffer_cb[cbBlockIndex]), stride_cb);
        crMeanOf16x16Blocks[7] = compute_sub_mean8x8_sse2_intrin(&(input_padded_picture_ptr->buffer_cr[crBlockIndex]), stride_cr);

        // (2,0)
        cbBlockIndex = inputCbOriginIndex + (stride_cb << 4);
        crBlockIndex = inputCrOriginIndex + (stride_cr << 4);
        cbMeanOf16x16Blocks[8] = compute_sub_mean8x8_sse2_intrin(&(input_padded_picture_ptr->buffer_cb[cbBlockIndex]), stride_cb);
        crMeanOf16x16Blocks[8] = compute_sub_mean8x8_sse2_intrin(&(input_padded_picture_ptr->buffer_cr[crBlockIndex]), stride_cr);

        // (2,1)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[9] = compute_sub_mean8x8_sse2_intrin(&(input_padded_picture_ptr->buffer_cb[cbBlockIndex]), stride_cb);
        crMeanOf16x16Blocks[9] = compute_sub_mean8x8_sse2_intrin(&(input_padded_picture_ptr->buffer_cr[crBlockIndex]), stride_cr);

        // (2,2)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[10] = compute_sub_mean8x8_sse2_intrin(&(input_padded_picture_ptr->buffer_cb[cbBlockIndex]), stride_cb);
        crMeanOf16x16Blocks[10] = compute_sub_mean8x8_sse2_intrin(&(input_padded_picture_ptr->buffer_cr[crBlockIndex]), stride_cr);

        // (2,3)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[11] = compute_sub_mean8x8_sse2_intrin(&(input_padded_picture_ptr->buffer_cb[cbBlockIndex]), stride_cb);
        crMeanOf16x16Blocks[11] = compute_sub_mean8x8_sse2_intrin(&(input_padded_picture_ptr->buffer_cr[crBlockIndex]), stride_cr);

        // (3,0)
        cbBlockIndex = inputCbOriginIndex + (stride_cb * 24);
        crBlockIndex = inputCrOriginIndex + (stride_cr * 24);
        cbMeanOf16x16Blocks[12] = compute_sub_mean8x8_sse2_intrin(&(input_padded_picture_ptr->buffer_cb[cbBlockIndex]), stride_cb);
        crMeanOf16x16Blocks[12] = compute_sub_mean8x8_sse2_intrin(&(input_padded_picture_ptr->buffer_cr[crBlockIndex]), stride_cr);

        // (3,1)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[13] = compute_sub_mean8x8_sse2_intrin(&(input_padded_picture_ptr->buffer_cb[cbBlockIndex]), stride_cb);
        crMeanOf16x16Blocks[13] = compute_sub_mean8x8_sse2_intrin(&(input_padded_picture_ptr->buffer_cr[crBlockIndex]), stride_cr);

        // (3,2)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[14] = compute_sub_mean8x8_sse2_intrin(&(input_padded_picture_ptr->buffer_cb[cbBlockIndex]), stride_cb);
        crMeanOf16x16Blocks[14] = compute_sub_mean8x8_sse2_intrin(&(input_padded_picture_ptr->buffer_cr[crBlockIndex]), stride_cr);

        // (3,3)
        cbBlockIndex = cbBlockIndex + 8;
        crBlockIndex = crBlockIndex + 8;
        cbMeanOf16x16Blocks[15] = compute_sub_mean8x8_sse2_intrin(&(input_padded_picture_ptr->buffer_cb[cbBlockIndex]), stride_cb);
        crMeanOf16x16Blocks[15] = compute_sub_mean8x8_sse2_intrin(&(input_padded_picture_ptr->buffer_cr[crBlockIndex]), stride_cr);
    }

    // 32x32
    cbMeanOf32x32Blocks[0] = (cbMeanOf16x16Blocks[0] + cbMeanOf16x16Blocks[1] + cbMeanOf16x16Blocks[4] + cbMeanOf16x16Blocks[5]) >> 2;
    crMeanOf32x32Blocks[0] = (crMeanOf16x16Blocks[0] + crMeanOf16x16Blocks[1] + crMeanOf16x16Blocks[4] + crMeanOf16x16Blocks[5]) >> 2;

    cbMeanOf32x32Blocks[1] = (cbMeanOf16x16Blocks[2] + cbMeanOf16x16Blocks[3] + cbMeanOf16x16Blocks[6] + cbMeanOf16x16Blocks[7]) >> 2;
    crMeanOf32x32Blocks[1] = (crMeanOf16x16Blocks[2] + crMeanOf16x16Blocks[3] + crMeanOf16x16Blocks[6] + crMeanOf16x16Blocks[7]) >> 2;

    cbMeanOf32x32Blocks[2] = (cbMeanOf16x16Blocks[8] + cbMeanOf16x16Blocks[9] + cbMeanOf16x16Blocks[12] + cbMeanOf16x16Blocks[13]) >> 2;
    crMeanOf32x32Blocks[2] = (crMeanOf16x16Blocks[8] + crMeanOf16x16Blocks[9] + crMeanOf16x16Blocks[12] + crMeanOf16x16Blocks[13]) >> 2;

    cbMeanOf32x32Blocks[3] = (cbMeanOf16x16Blocks[10] + cbMeanOf16x16Blocks[11] + cbMeanOf16x16Blocks[14] + cbMeanOf16x16Blocks[15]) >> 2;
    crMeanOf32x32Blocks[3] = (crMeanOf16x16Blocks[10] + crMeanOf16x16Blocks[11] + crMeanOf16x16Blocks[14] + crMeanOf16x16Blocks[15]) >> 2;

    // 64x64
    cbMeanOf64x64Blocks = (cbMeanOf32x32Blocks[0] + cbMeanOf32x32Blocks[1] + cbMeanOf32x32Blocks[3] + cbMeanOf32x32Blocks[3]) >> 2;
    crMeanOf64x64Blocks = (crMeanOf32x32Blocks[0] + crMeanOf32x32Blocks[1] + crMeanOf32x32Blocks[3] + crMeanOf32x32Blocks[3]) >> 2;
    // 16x16 mean
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_0] = (uint8_t)(cbMeanOf16x16Blocks[0] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_1] = (uint8_t)(cbMeanOf16x16Blocks[1] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_2] = (uint8_t)(cbMeanOf16x16Blocks[2] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_3] = (uint8_t)(cbMeanOf16x16Blocks[3] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_4] = (uint8_t)(cbMeanOf16x16Blocks[4] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_5] = (uint8_t)(cbMeanOf16x16Blocks[5] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_6] = (uint8_t)(cbMeanOf16x16Blocks[6] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_7] = (uint8_t)(cbMeanOf16x16Blocks[7] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_8] = (uint8_t)(cbMeanOf16x16Blocks[8] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_9] = (uint8_t)(cbMeanOf16x16Blocks[9] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_10] = (uint8_t)(cbMeanOf16x16Blocks[10] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_11] = (uint8_t)(cbMeanOf16x16Blocks[11] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_12] = (uint8_t)(cbMeanOf16x16Blocks[12] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_13] = (uint8_t)(cbMeanOf16x16Blocks[13] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_14] = (uint8_t)(cbMeanOf16x16Blocks[14] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_15] = (uint8_t)(cbMeanOf16x16Blocks[15] >> MEAN_PRECISION);

    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_0] = (uint8_t)(crMeanOf16x16Blocks[0] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_1] = (uint8_t)(crMeanOf16x16Blocks[1] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_2] = (uint8_t)(crMeanOf16x16Blocks[2] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_3] = (uint8_t)(crMeanOf16x16Blocks[3] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_4] = (uint8_t)(crMeanOf16x16Blocks[4] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_5] = (uint8_t)(crMeanOf16x16Blocks[5] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_6] = (uint8_t)(crMeanOf16x16Blocks[6] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_7] = (uint8_t)(crMeanOf16x16Blocks[7] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_8] = (uint8_t)(crMeanOf16x16Blocks[8] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_9] = (uint8_t)(crMeanOf16x16Blocks[9] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_10] = (uint8_t)(crMeanOf16x16Blocks[10] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_11] = (uint8_t)(crMeanOf16x16Blocks[11] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_12] = (uint8_t)(crMeanOf16x16Blocks[12] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_13] = (uint8_t)(crMeanOf16x16Blocks[13] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_14] = (uint8_t)(crMeanOf16x16Blocks[14] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_16x16_15] = (uint8_t)(crMeanOf16x16Blocks[15] >> MEAN_PRECISION);

    // 32x32 mean
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_32x32_0] = (uint8_t)(cbMeanOf32x32Blocks[0] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_32x32_1] = (uint8_t)(cbMeanOf32x32Blocks[1] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_32x32_2] = (uint8_t)(cbMeanOf32x32Blocks[2] >> MEAN_PRECISION);
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_32x32_3] = (uint8_t)(cbMeanOf32x32Blocks[3] >> MEAN_PRECISION);

    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_32x32_0] = (uint8_t)(crMeanOf32x32Blocks[0] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_32x32_1] = (uint8_t)(crMeanOf32x32Blocks[1] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_32x32_2] = (uint8_t)(crMeanOf32x32Blocks[2] >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_32x32_3] = (uint8_t)(crMeanOf32x32Blocks[3] >> MEAN_PRECISION);

    // 64x64 mean
    picture_control_set_ptr->cbMean[lcuCodingOrder][ME_TIER_ZERO_PU_64x64] = (uint8_t)(cbMeanOf64x64Blocks >> MEAN_PRECISION);
    picture_control_set_ptr->crMean[lcuCodingOrder][ME_TIER_ZERO_PU_64x64] = (uint8_t)(crMeanOf64x64Blocks >> MEAN_PRECISION);

    return return_error;
}

/*******************************************
* ComputeBlockMeanComputeVariance
*   computes the variance and the block mean of all CUs inside the tree block
*******************************************/
EbErrorType ComputeBlockMeanComputeVariance(
    SequenceControlSet        *sequence_control_set_ptr,
    PictureParentControlSet   *picture_control_set_ptr,          // input parameter, Picture Control Set Ptr
    EbPictureBufferDesc       *input_padded_picture_ptr,         // input parameter, Input Padded Picture
    uint32_t                       sb_index,                // input parameter, SB address
    uint32_t                       inputLumaOriginIndex)          // input parameter, SB index, used to point to source/reference samples
{
    EbErrorType return_error = EB_ErrorNone;

    uint32_t blockIndex;

    uint64_t mean_of8x8_blocks[64];
    uint64_t meanOf8x8SquaredValuesBlocks[64];

    uint64_t meanOf16x16Blocks[16];
    uint64_t meanOf16x16SquaredValuesBlocks[16];

    uint64_t meanOf32x32Blocks[4];
    uint64_t meanOf32x32SquaredValuesBlocks[4];

    uint64_t meanOf64x64Blocks;
    uint64_t meanOf64x64SquaredValuesBlocks;

    // (0,0)
    blockIndex = inputLumaOriginIndex;
    if (sequence_control_set_ptr->block_mean_calc_prec == BLOCK_MEAN_PREC_FULL) {
        mean_of8x8_blocks[0] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[0] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (0,1)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[1] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[1] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (0,2)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[2] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[2] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (0,3)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[3] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[3] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (0,4)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[4] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[4] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (0,5)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[5] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[5] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (0,6)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[6] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[6] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (0,7)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[7] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[7] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (1,0)
        blockIndex = inputLumaOriginIndex + (input_padded_picture_ptr->stride_y << 3);
        mean_of8x8_blocks[8] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[8] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (1,1)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[9] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[9] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (1,2)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[10] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[10] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (1,3)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[11] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[11] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (1,4)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[12] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[12] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (1,5)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[13] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[13] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (1,6)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[14] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[14] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (1,7)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[15] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[15] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (2,0)
        blockIndex = inputLumaOriginIndex + (input_padded_picture_ptr->stride_y << 4);
        mean_of8x8_blocks[16] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[16] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (2,1)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[17] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[17] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (2,2)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[18] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[18] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (2,3)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[19] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[19] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        /// (2,4)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[20] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[20] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (2,5)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[21] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[21] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (2,6)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[22] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[22] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (2,7)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[23] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[23] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (3,0)
        blockIndex = inputLumaOriginIndex + (input_padded_picture_ptr->stride_y << 3) + (input_padded_picture_ptr->stride_y << 4);
        mean_of8x8_blocks[24] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[24] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (3,1)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[25] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[25] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (3,2)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[26] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[26] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (3,3)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[27] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[27] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (3,4)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[28] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[28] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (3,5)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[29] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[29] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (3,6)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[30] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[30] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (3,7)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[31] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[31] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (4,0)
        blockIndex = inputLumaOriginIndex + (input_padded_picture_ptr->stride_y << 5);
        mean_of8x8_blocks[32] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[32] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (4,1)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[33] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[33] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (4,2)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[34] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[34] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (4,3)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[35] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[35] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (4,4)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[36] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[36] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (4,5)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[37] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[37] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (4,6)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[38] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[38] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (4,7)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[39] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[39] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (5,0)
        blockIndex = inputLumaOriginIndex + (input_padded_picture_ptr->stride_y << 3) + (input_padded_picture_ptr->stride_y << 5);
        mean_of8x8_blocks[40] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[40] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (5,1)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[41] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[41] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (5,2)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[42] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[42] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (5,3)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[43] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[43] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (5,4)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[44] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[44] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (5,5)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[45] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[45] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (5,6)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[46] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[46] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (5,7)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[47] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[47] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (6,0)
        blockIndex = inputLumaOriginIndex + (input_padded_picture_ptr->stride_y << 4) + (input_padded_picture_ptr->stride_y << 5);
        mean_of8x8_blocks[48] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[48] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (6,1)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[49] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[49] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (6,2)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[50] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[50] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (6,3)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[51] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[51] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (6,4)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[52] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[52] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (6,5)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[53] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[53] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (6,6)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[54] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[54] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (6,7)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[55] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[55] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (7,0)
        blockIndex = inputLumaOriginIndex + (input_padded_picture_ptr->stride_y << 3) + (input_padded_picture_ptr->stride_y << 4) + (input_padded_picture_ptr->stride_y << 5);
        mean_of8x8_blocks[56] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[56] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (7,1)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[57] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[57] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (7,2)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[58] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[58] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (7,3)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[59] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[59] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (7,4)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[60] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[60] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (7,5)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[61] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[61] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (7,6)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[62] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[62] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);

        // (7,7)
        blockIndex = blockIndex + 8;
        mean_of8x8_blocks[63] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
        meanOf8x8SquaredValuesBlocks[63] = compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), input_padded_picture_ptr->stride_y, 8, 8);
    }
    else {
        const uint16_t stride_y = input_padded_picture_ptr->stride_y;

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), stride_y, &mean_of8x8_blocks[0], &meanOf8x8SquaredValuesBlocks[0]);

        // (0,1)
        blockIndex = blockIndex + 32;

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), stride_y, &mean_of8x8_blocks[4], &meanOf8x8SquaredValuesBlocks[4]);

        // (0,5)
        blockIndex = blockIndex + 24;

        // (1,0)
        blockIndex = inputLumaOriginIndex + (stride_y << 3);

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), stride_y, &mean_of8x8_blocks[8], &meanOf8x8SquaredValuesBlocks[8]);

        // (1,1)
        blockIndex = blockIndex + 32;

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), stride_y, &mean_of8x8_blocks[12], &meanOf8x8SquaredValuesBlocks[12]);

        // (1,5)
        blockIndex = blockIndex + 24;

        // (2,0)
        blockIndex = inputLumaOriginIndex + (stride_y << 4);

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), stride_y, &mean_of8x8_blocks[16], &meanOf8x8SquaredValuesBlocks[16]);

        // (2,1)
        blockIndex = blockIndex + 32;

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), stride_y, &mean_of8x8_blocks[20], &meanOf8x8SquaredValuesBlocks[20]);

        // (2,5)
        blockIndex = blockIndex + 24;

        // (3,0)
        blockIndex = inputLumaOriginIndex + (stride_y << 3) + (stride_y << 4);

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), stride_y, &mean_of8x8_blocks[24], &meanOf8x8SquaredValuesBlocks[24]);

        // (3,1)
        blockIndex = blockIndex + 32;

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), stride_y, &mean_of8x8_blocks[28], &meanOf8x8SquaredValuesBlocks[28]);

        // (3,5)
        blockIndex = blockIndex + 24;

        // (4,0)
        blockIndex = inputLumaOriginIndex + (stride_y << 5);

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), stride_y, &mean_of8x8_blocks[32], &meanOf8x8SquaredValuesBlocks[32]);

        // (4,1)
        blockIndex = blockIndex + 32;

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), stride_y, &mean_of8x8_blocks[36], &meanOf8x8SquaredValuesBlocks[36]);

        // (4,5)
        blockIndex = blockIndex + 24;

        // (5,0)
        blockIndex = inputLumaOriginIndex + (stride_y << 3) + (stride_y << 5);

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), stride_y, &mean_of8x8_blocks[40], &meanOf8x8SquaredValuesBlocks[40]);

        // (5,1)
        blockIndex = blockIndex + 32;

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), stride_y, &mean_of8x8_blocks[44], &meanOf8x8SquaredValuesBlocks[44]);

        // (5,5)
        blockIndex = blockIndex + 24;

        // (6,0)
        blockIndex = inputLumaOriginIndex + (stride_y << 4) + (stride_y << 5);

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), stride_y, &mean_of8x8_blocks[48], &meanOf8x8SquaredValuesBlocks[48]);

        // (6,1)
        blockIndex = blockIndex + 32;

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), stride_y, &mean_of8x8_blocks[52], &meanOf8x8SquaredValuesBlocks[52]);

        // (6,5)
        blockIndex = blockIndex + 24;

        // (7,0)
        blockIndex = inputLumaOriginIndex + (stride_y << 3) + (stride_y << 4) + (stride_y << 5);

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), stride_y, &mean_of8x8_blocks[56], &meanOf8x8SquaredValuesBlocks[56]);

        // (7,1)
        blockIndex = blockIndex + 32;

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[blockIndex]), stride_y, &mean_of8x8_blocks[60], &meanOf8x8SquaredValuesBlocks[60]);

    }

    // 16x16
    meanOf16x16Blocks[0] = (mean_of8x8_blocks[0] + mean_of8x8_blocks[1] + mean_of8x8_blocks[8] + mean_of8x8_blocks[9]) >> 2;
    meanOf16x16Blocks[1] = (mean_of8x8_blocks[2] + mean_of8x8_blocks[3] + mean_of8x8_blocks[10] + mean_of8x8_blocks[11]) >> 2;
    meanOf16x16Blocks[2] = (mean_of8x8_blocks[4] + mean_of8x8_blocks[5] + mean_of8x8_blocks[12] + mean_of8x8_blocks[13]) >> 2;
    meanOf16x16Blocks[3] = (mean_of8x8_blocks[6] + mean_of8x8_blocks[7] + mean_of8x8_blocks[14] + mean_of8x8_blocks[15]) >> 2;

    meanOf16x16Blocks[4] = (mean_of8x8_blocks[16] + mean_of8x8_blocks[17] + mean_of8x8_blocks[24] + mean_of8x8_blocks[25]) >> 2;
    meanOf16x16Blocks[5] = (mean_of8x8_blocks[18] + mean_of8x8_blocks[19] + mean_of8x8_blocks[26] + mean_of8x8_blocks[27]) >> 2;
    meanOf16x16Blocks[6] = (mean_of8x8_blocks[20] + mean_of8x8_blocks[21] + mean_of8x8_blocks[28] + mean_of8x8_blocks[29]) >> 2;
    meanOf16x16Blocks[7] = (mean_of8x8_blocks[22] + mean_of8x8_blocks[23] + mean_of8x8_blocks[30] + mean_of8x8_blocks[31]) >> 2;

    meanOf16x16Blocks[8] = (mean_of8x8_blocks[32] + mean_of8x8_blocks[33] + mean_of8x8_blocks[40] + mean_of8x8_blocks[41]) >> 2;
    meanOf16x16Blocks[9] = (mean_of8x8_blocks[34] + mean_of8x8_blocks[35] + mean_of8x8_blocks[42] + mean_of8x8_blocks[43]) >> 2;
    meanOf16x16Blocks[10] = (mean_of8x8_blocks[36] + mean_of8x8_blocks[37] + mean_of8x8_blocks[44] + mean_of8x8_blocks[45]) >> 2;
    meanOf16x16Blocks[11] = (mean_of8x8_blocks[38] + mean_of8x8_blocks[39] + mean_of8x8_blocks[46] + mean_of8x8_blocks[47]) >> 2;

    meanOf16x16Blocks[12] = (mean_of8x8_blocks[48] + mean_of8x8_blocks[49] + mean_of8x8_blocks[56] + mean_of8x8_blocks[57]) >> 2;
    meanOf16x16Blocks[13] = (mean_of8x8_blocks[50] + mean_of8x8_blocks[51] + mean_of8x8_blocks[58] + mean_of8x8_blocks[59]) >> 2;
    meanOf16x16Blocks[14] = (mean_of8x8_blocks[52] + mean_of8x8_blocks[53] + mean_of8x8_blocks[60] + mean_of8x8_blocks[61]) >> 2;
    meanOf16x16Blocks[15] = (mean_of8x8_blocks[54] + mean_of8x8_blocks[55] + mean_of8x8_blocks[62] + mean_of8x8_blocks[63]) >> 2;

    meanOf16x16SquaredValuesBlocks[0] = (meanOf8x8SquaredValuesBlocks[0] + meanOf8x8SquaredValuesBlocks[1] + meanOf8x8SquaredValuesBlocks[8] + meanOf8x8SquaredValuesBlocks[9]) >> 2;
    meanOf16x16SquaredValuesBlocks[1] = (meanOf8x8SquaredValuesBlocks[2] + meanOf8x8SquaredValuesBlocks[3] + meanOf8x8SquaredValuesBlocks[10] + meanOf8x8SquaredValuesBlocks[11]) >> 2;
    meanOf16x16SquaredValuesBlocks[2] = (meanOf8x8SquaredValuesBlocks[4] + meanOf8x8SquaredValuesBlocks[5] + meanOf8x8SquaredValuesBlocks[12] + meanOf8x8SquaredValuesBlocks[13]) >> 2;
    meanOf16x16SquaredValuesBlocks[3] = (meanOf8x8SquaredValuesBlocks[6] + meanOf8x8SquaredValuesBlocks[7] + meanOf8x8SquaredValuesBlocks[14] + meanOf8x8SquaredValuesBlocks[15]) >> 2;

    meanOf16x16SquaredValuesBlocks[4] = (meanOf8x8SquaredValuesBlocks[16] + meanOf8x8SquaredValuesBlocks[17] + meanOf8x8SquaredValuesBlocks[24] + meanOf8x8SquaredValuesBlocks[25]) >> 2;
    meanOf16x16SquaredValuesBlocks[5] = (meanOf8x8SquaredValuesBlocks[18] + meanOf8x8SquaredValuesBlocks[19] + meanOf8x8SquaredValuesBlocks[26] + meanOf8x8SquaredValuesBlocks[27]) >> 2;
    meanOf16x16SquaredValuesBlocks[6] = (meanOf8x8SquaredValuesBlocks[20] + meanOf8x8SquaredValuesBlocks[21] + meanOf8x8SquaredValuesBlocks[28] + meanOf8x8SquaredValuesBlocks[29]) >> 2;
    meanOf16x16SquaredValuesBlocks[7] = (meanOf8x8SquaredValuesBlocks[22] + meanOf8x8SquaredValuesBlocks[23] + meanOf8x8SquaredValuesBlocks[30] + meanOf8x8SquaredValuesBlocks[31]) >> 2;

    meanOf16x16SquaredValuesBlocks[8] = (meanOf8x8SquaredValuesBlocks[32] + meanOf8x8SquaredValuesBlocks[33] + meanOf8x8SquaredValuesBlocks[40] + meanOf8x8SquaredValuesBlocks[41]) >> 2;
    meanOf16x16SquaredValuesBlocks[9] = (meanOf8x8SquaredValuesBlocks[34] + meanOf8x8SquaredValuesBlocks[35] + meanOf8x8SquaredValuesBlocks[42] + meanOf8x8SquaredValuesBlocks[43]) >> 2;
    meanOf16x16SquaredValuesBlocks[10] = (meanOf8x8SquaredValuesBlocks[36] + meanOf8x8SquaredValuesBlocks[37] + meanOf8x8SquaredValuesBlocks[44] + meanOf8x8SquaredValuesBlocks[45]) >> 2;
    meanOf16x16SquaredValuesBlocks[11] = (meanOf8x8SquaredValuesBlocks[38] + meanOf8x8SquaredValuesBlocks[39] + meanOf8x8SquaredValuesBlocks[46] + meanOf8x8SquaredValuesBlocks[47]) >> 2;

    meanOf16x16SquaredValuesBlocks[12] = (meanOf8x8SquaredValuesBlocks[48] + meanOf8x8SquaredValuesBlocks[49] + meanOf8x8SquaredValuesBlocks[56] + meanOf8x8SquaredValuesBlocks[57]) >> 2;
    meanOf16x16SquaredValuesBlocks[13] = (meanOf8x8SquaredValuesBlocks[50] + meanOf8x8SquaredValuesBlocks[51] + meanOf8x8SquaredValuesBlocks[58] + meanOf8x8SquaredValuesBlocks[59]) >> 2;
    meanOf16x16SquaredValuesBlocks[14] = (meanOf8x8SquaredValuesBlocks[52] + meanOf8x8SquaredValuesBlocks[53] + meanOf8x8SquaredValuesBlocks[60] + meanOf8x8SquaredValuesBlocks[61]) >> 2;
    meanOf16x16SquaredValuesBlocks[15] = (meanOf8x8SquaredValuesBlocks[54] + meanOf8x8SquaredValuesBlocks[55] + meanOf8x8SquaredValuesBlocks[62] + meanOf8x8SquaredValuesBlocks[63]) >> 2;

    // 32x32
    meanOf32x32Blocks[0] = (meanOf16x16Blocks[0] + meanOf16x16Blocks[1] + meanOf16x16Blocks[4] + meanOf16x16Blocks[5]) >> 2;
    meanOf32x32Blocks[1] = (meanOf16x16Blocks[2] + meanOf16x16Blocks[3] + meanOf16x16Blocks[6] + meanOf16x16Blocks[7]) >> 2;
    meanOf32x32Blocks[2] = (meanOf16x16Blocks[8] + meanOf16x16Blocks[9] + meanOf16x16Blocks[12] + meanOf16x16Blocks[13]) >> 2;
    meanOf32x32Blocks[3] = (meanOf16x16Blocks[10] + meanOf16x16Blocks[11] + meanOf16x16Blocks[14] + meanOf16x16Blocks[15]) >> 2;

    meanOf32x32SquaredValuesBlocks[0] = (meanOf16x16SquaredValuesBlocks[0] + meanOf16x16SquaredValuesBlocks[1] + meanOf16x16SquaredValuesBlocks[4] + meanOf16x16SquaredValuesBlocks[5]) >> 2;
    meanOf32x32SquaredValuesBlocks[1] = (meanOf16x16SquaredValuesBlocks[2] + meanOf16x16SquaredValuesBlocks[3] + meanOf16x16SquaredValuesBlocks[6] + meanOf16x16SquaredValuesBlocks[7]) >> 2;
    meanOf32x32SquaredValuesBlocks[2] = (meanOf16x16SquaredValuesBlocks[8] + meanOf16x16SquaredValuesBlocks[9] + meanOf16x16SquaredValuesBlocks[12] + meanOf16x16SquaredValuesBlocks[13]) >> 2;
    meanOf32x32SquaredValuesBlocks[3] = (meanOf16x16SquaredValuesBlocks[10] + meanOf16x16SquaredValuesBlocks[11] + meanOf16x16SquaredValuesBlocks[14] + meanOf16x16SquaredValuesBlocks[15]) >> 2;

    // 64x64
    meanOf64x64Blocks = (meanOf32x32Blocks[0] + meanOf32x32Blocks[1] + meanOf32x32Blocks[2] + meanOf32x32Blocks[3]) >> 2;
    meanOf64x64SquaredValuesBlocks = (meanOf32x32SquaredValuesBlocks[0] + meanOf32x32SquaredValuesBlocks[1] + meanOf32x32SquaredValuesBlocks[2] + meanOf32x32SquaredValuesBlocks[3]) >> 2;

    // 8x8 means
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_0] = (uint8_t)(mean_of8x8_blocks[0] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_1] = (uint8_t)(mean_of8x8_blocks[1] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_2] = (uint8_t)(mean_of8x8_blocks[2] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_3] = (uint8_t)(mean_of8x8_blocks[3] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_4] = (uint8_t)(mean_of8x8_blocks[4] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_5] = (uint8_t)(mean_of8x8_blocks[5] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_6] = (uint8_t)(mean_of8x8_blocks[6] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_7] = (uint8_t)(mean_of8x8_blocks[7] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_8] = (uint8_t)(mean_of8x8_blocks[8] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_9] = (uint8_t)(mean_of8x8_blocks[9] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_10] = (uint8_t)(mean_of8x8_blocks[10] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_11] = (uint8_t)(mean_of8x8_blocks[11] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_12] = (uint8_t)(mean_of8x8_blocks[12] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_13] = (uint8_t)(mean_of8x8_blocks[13] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_14] = (uint8_t)(mean_of8x8_blocks[14] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_15] = (uint8_t)(mean_of8x8_blocks[15] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_16] = (uint8_t)(mean_of8x8_blocks[16] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_17] = (uint8_t)(mean_of8x8_blocks[17] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_18] = (uint8_t)(mean_of8x8_blocks[18] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_19] = (uint8_t)(mean_of8x8_blocks[19] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_20] = (uint8_t)(mean_of8x8_blocks[20] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_21] = (uint8_t)(mean_of8x8_blocks[21] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_22] = (uint8_t)(mean_of8x8_blocks[22] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_23] = (uint8_t)(mean_of8x8_blocks[23] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_24] = (uint8_t)(mean_of8x8_blocks[24] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_25] = (uint8_t)(mean_of8x8_blocks[25] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_26] = (uint8_t)(mean_of8x8_blocks[26] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_27] = (uint8_t)(mean_of8x8_blocks[27] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_28] = (uint8_t)(mean_of8x8_blocks[28] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_29] = (uint8_t)(mean_of8x8_blocks[29] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_30] = (uint8_t)(mean_of8x8_blocks[30] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_31] = (uint8_t)(mean_of8x8_blocks[31] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_32] = (uint8_t)(mean_of8x8_blocks[32] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_33] = (uint8_t)(mean_of8x8_blocks[33] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_34] = (uint8_t)(mean_of8x8_blocks[34] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_35] = (uint8_t)(mean_of8x8_blocks[35] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_36] = (uint8_t)(mean_of8x8_blocks[36] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_37] = (uint8_t)(mean_of8x8_blocks[37] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_38] = (uint8_t)(mean_of8x8_blocks[38] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_39] = (uint8_t)(mean_of8x8_blocks[39] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_40] = (uint8_t)(mean_of8x8_blocks[40] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_41] = (uint8_t)(mean_of8x8_blocks[41] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_42] = (uint8_t)(mean_of8x8_blocks[42] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_43] = (uint8_t)(mean_of8x8_blocks[43] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_44] = (uint8_t)(mean_of8x8_blocks[44] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_45] = (uint8_t)(mean_of8x8_blocks[45] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_46] = (uint8_t)(mean_of8x8_blocks[46] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_47] = (uint8_t)(mean_of8x8_blocks[47] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_48] = (uint8_t)(mean_of8x8_blocks[48] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_49] = (uint8_t)(mean_of8x8_blocks[49] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_50] = (uint8_t)(mean_of8x8_blocks[50] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_51] = (uint8_t)(mean_of8x8_blocks[51] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_52] = (uint8_t)(mean_of8x8_blocks[52] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_53] = (uint8_t)(mean_of8x8_blocks[53] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_54] = (uint8_t)(mean_of8x8_blocks[54] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_55] = (uint8_t)(mean_of8x8_blocks[55] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_56] = (uint8_t)(mean_of8x8_blocks[56] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_57] = (uint8_t)(mean_of8x8_blocks[57] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_58] = (uint8_t)(mean_of8x8_blocks[58] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_59] = (uint8_t)(mean_of8x8_blocks[59] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_60] = (uint8_t)(mean_of8x8_blocks[60] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_61] = (uint8_t)(mean_of8x8_blocks[61] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_62] = (uint8_t)(mean_of8x8_blocks[62] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_63] = (uint8_t)(mean_of8x8_blocks[63] >> MEAN_PRECISION);

    // 16x16 mean
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_16x16_0] = (uint8_t)(meanOf16x16Blocks[0] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_16x16_1] = (uint8_t)(meanOf16x16Blocks[1] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_16x16_2] = (uint8_t)(meanOf16x16Blocks[2] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_16x16_3] = (uint8_t)(meanOf16x16Blocks[3] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_16x16_4] = (uint8_t)(meanOf16x16Blocks[4] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_16x16_5] = (uint8_t)(meanOf16x16Blocks[5] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_16x16_6] = (uint8_t)(meanOf16x16Blocks[6] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_16x16_7] = (uint8_t)(meanOf16x16Blocks[7] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_16x16_8] = (uint8_t)(meanOf16x16Blocks[8] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_16x16_9] = (uint8_t)(meanOf16x16Blocks[9] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_16x16_10] = (uint8_t)(meanOf16x16Blocks[10] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_16x16_11] = (uint8_t)(meanOf16x16Blocks[11] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_16x16_12] = (uint8_t)(meanOf16x16Blocks[12] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_16x16_13] = (uint8_t)(meanOf16x16Blocks[13] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_16x16_14] = (uint8_t)(meanOf16x16Blocks[14] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_16x16_15] = (uint8_t)(meanOf16x16Blocks[15] >> MEAN_PRECISION);

    // 32x32 mean
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_32x32_0] = (uint8_t)(meanOf32x32Blocks[0] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_32x32_1] = (uint8_t)(meanOf32x32Blocks[1] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_32x32_2] = (uint8_t)(meanOf32x32Blocks[2] >> MEAN_PRECISION);
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_32x32_3] = (uint8_t)(meanOf32x32Blocks[3] >> MEAN_PRECISION);

    // 64x64 mean
    picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_64x64] = (uint8_t)(meanOf64x64Blocks >> MEAN_PRECISION);

    // 8x8 variances
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_0] = (uint16_t)((meanOf8x8SquaredValuesBlocks[0] - (mean_of8x8_blocks[0] * mean_of8x8_blocks[0])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_1] = (uint16_t)((meanOf8x8SquaredValuesBlocks[1] - (mean_of8x8_blocks[1] * mean_of8x8_blocks[1])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_2] = (uint16_t)((meanOf8x8SquaredValuesBlocks[2] - (mean_of8x8_blocks[2] * mean_of8x8_blocks[2])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_3] = (uint16_t)((meanOf8x8SquaredValuesBlocks[3] - (mean_of8x8_blocks[3] * mean_of8x8_blocks[3])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_4] = (uint16_t)((meanOf8x8SquaredValuesBlocks[4] - (mean_of8x8_blocks[4] * mean_of8x8_blocks[4])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_5] = (uint16_t)((meanOf8x8SquaredValuesBlocks[5] - (mean_of8x8_blocks[5] * mean_of8x8_blocks[5])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_6] = (uint16_t)((meanOf8x8SquaredValuesBlocks[6] - (mean_of8x8_blocks[6] * mean_of8x8_blocks[6])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_7] = (uint16_t)((meanOf8x8SquaredValuesBlocks[7] - (mean_of8x8_blocks[7] * mean_of8x8_blocks[7])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_8] = (uint16_t)((meanOf8x8SquaredValuesBlocks[8] - (mean_of8x8_blocks[8] * mean_of8x8_blocks[8])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_9] = (uint16_t)((meanOf8x8SquaredValuesBlocks[9] - (mean_of8x8_blocks[9] * mean_of8x8_blocks[9])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_10] = (uint16_t)((meanOf8x8SquaredValuesBlocks[10] - (mean_of8x8_blocks[10] * mean_of8x8_blocks[10])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_11] = (uint16_t)((meanOf8x8SquaredValuesBlocks[11] - (mean_of8x8_blocks[11] * mean_of8x8_blocks[11])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_12] = (uint16_t)((meanOf8x8SquaredValuesBlocks[12] - (mean_of8x8_blocks[12] * mean_of8x8_blocks[12])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_13] = (uint16_t)((meanOf8x8SquaredValuesBlocks[13] - (mean_of8x8_blocks[13] * mean_of8x8_blocks[13])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_14] = (uint16_t)((meanOf8x8SquaredValuesBlocks[14] - (mean_of8x8_blocks[14] * mean_of8x8_blocks[14])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_15] = (uint16_t)((meanOf8x8SquaredValuesBlocks[15] - (mean_of8x8_blocks[15] * mean_of8x8_blocks[15])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_16] = (uint16_t)((meanOf8x8SquaredValuesBlocks[16] - (mean_of8x8_blocks[16] * mean_of8x8_blocks[16])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_17] = (uint16_t)((meanOf8x8SquaredValuesBlocks[17] - (mean_of8x8_blocks[17] * mean_of8x8_blocks[17])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_18] = (uint16_t)((meanOf8x8SquaredValuesBlocks[18] - (mean_of8x8_blocks[18] * mean_of8x8_blocks[18])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_19] = (uint16_t)((meanOf8x8SquaredValuesBlocks[19] - (mean_of8x8_blocks[19] * mean_of8x8_blocks[19])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_20] = (uint16_t)((meanOf8x8SquaredValuesBlocks[20] - (mean_of8x8_blocks[20] * mean_of8x8_blocks[20])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_21] = (uint16_t)((meanOf8x8SquaredValuesBlocks[21] - (mean_of8x8_blocks[21] * mean_of8x8_blocks[21])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_22] = (uint16_t)((meanOf8x8SquaredValuesBlocks[22] - (mean_of8x8_blocks[22] * mean_of8x8_blocks[22])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_23] = (uint16_t)((meanOf8x8SquaredValuesBlocks[23] - (mean_of8x8_blocks[23] * mean_of8x8_blocks[23])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_24] = (uint16_t)((meanOf8x8SquaredValuesBlocks[24] - (mean_of8x8_blocks[24] * mean_of8x8_blocks[24])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_25] = (uint16_t)((meanOf8x8SquaredValuesBlocks[25] - (mean_of8x8_blocks[25] * mean_of8x8_blocks[25])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_26] = (uint16_t)((meanOf8x8SquaredValuesBlocks[26] - (mean_of8x8_blocks[26] * mean_of8x8_blocks[26])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_27] = (uint16_t)((meanOf8x8SquaredValuesBlocks[27] - (mean_of8x8_blocks[27] * mean_of8x8_blocks[27])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_28] = (uint16_t)((meanOf8x8SquaredValuesBlocks[28] - (mean_of8x8_blocks[28] * mean_of8x8_blocks[28])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_29] = (uint16_t)((meanOf8x8SquaredValuesBlocks[29] - (mean_of8x8_blocks[29] * mean_of8x8_blocks[29])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_30] = (uint16_t)((meanOf8x8SquaredValuesBlocks[30] - (mean_of8x8_blocks[30] * mean_of8x8_blocks[30])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_31] = (uint16_t)((meanOf8x8SquaredValuesBlocks[31] - (mean_of8x8_blocks[31] * mean_of8x8_blocks[31])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_32] = (uint16_t)((meanOf8x8SquaredValuesBlocks[32] - (mean_of8x8_blocks[32] * mean_of8x8_blocks[32])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_33] = (uint16_t)((meanOf8x8SquaredValuesBlocks[33] - (mean_of8x8_blocks[33] * mean_of8x8_blocks[33])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_34] = (uint16_t)((meanOf8x8SquaredValuesBlocks[34] - (mean_of8x8_blocks[34] * mean_of8x8_blocks[34])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_35] = (uint16_t)((meanOf8x8SquaredValuesBlocks[35] - (mean_of8x8_blocks[35] * mean_of8x8_blocks[35])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_36] = (uint16_t)((meanOf8x8SquaredValuesBlocks[36] - (mean_of8x8_blocks[36] * mean_of8x8_blocks[36])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_37] = (uint16_t)((meanOf8x8SquaredValuesBlocks[37] - (mean_of8x8_blocks[37] * mean_of8x8_blocks[37])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_38] = (uint16_t)((meanOf8x8SquaredValuesBlocks[38] - (mean_of8x8_blocks[38] * mean_of8x8_blocks[38])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_39] = (uint16_t)((meanOf8x8SquaredValuesBlocks[39] - (mean_of8x8_blocks[39] * mean_of8x8_blocks[39])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_40] = (uint16_t)((meanOf8x8SquaredValuesBlocks[40] - (mean_of8x8_blocks[40] * mean_of8x8_blocks[40])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_41] = (uint16_t)((meanOf8x8SquaredValuesBlocks[41] - (mean_of8x8_blocks[41] * mean_of8x8_blocks[41])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_42] = (uint16_t)((meanOf8x8SquaredValuesBlocks[42] - (mean_of8x8_blocks[42] * mean_of8x8_blocks[42])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_43] = (uint16_t)((meanOf8x8SquaredValuesBlocks[43] - (mean_of8x8_blocks[43] * mean_of8x8_blocks[43])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_44] = (uint16_t)((meanOf8x8SquaredValuesBlocks[44] - (mean_of8x8_blocks[44] * mean_of8x8_blocks[44])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_45] = (uint16_t)((meanOf8x8SquaredValuesBlocks[45] - (mean_of8x8_blocks[45] * mean_of8x8_blocks[45])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_46] = (uint16_t)((meanOf8x8SquaredValuesBlocks[46] - (mean_of8x8_blocks[46] * mean_of8x8_blocks[46])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_47] = (uint16_t)((meanOf8x8SquaredValuesBlocks[47] - (mean_of8x8_blocks[47] * mean_of8x8_blocks[47])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_48] = (uint16_t)((meanOf8x8SquaredValuesBlocks[48] - (mean_of8x8_blocks[48] * mean_of8x8_blocks[48])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_49] = (uint16_t)((meanOf8x8SquaredValuesBlocks[49] - (mean_of8x8_blocks[49] * mean_of8x8_blocks[49])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_50] = (uint16_t)((meanOf8x8SquaredValuesBlocks[50] - (mean_of8x8_blocks[50] * mean_of8x8_blocks[50])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_51] = (uint16_t)((meanOf8x8SquaredValuesBlocks[51] - (mean_of8x8_blocks[51] * mean_of8x8_blocks[51])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_52] = (uint16_t)((meanOf8x8SquaredValuesBlocks[52] - (mean_of8x8_blocks[52] * mean_of8x8_blocks[52])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_53] = (uint16_t)((meanOf8x8SquaredValuesBlocks[53] - (mean_of8x8_blocks[53] * mean_of8x8_blocks[53])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_54] = (uint16_t)((meanOf8x8SquaredValuesBlocks[54] - (mean_of8x8_blocks[54] * mean_of8x8_blocks[54])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_55] = (uint16_t)((meanOf8x8SquaredValuesBlocks[55] - (mean_of8x8_blocks[55] * mean_of8x8_blocks[55])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_56] = (uint16_t)((meanOf8x8SquaredValuesBlocks[56] - (mean_of8x8_blocks[56] * mean_of8x8_blocks[56])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_57] = (uint16_t)((meanOf8x8SquaredValuesBlocks[57] - (mean_of8x8_blocks[57] * mean_of8x8_blocks[57])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_58] = (uint16_t)((meanOf8x8SquaredValuesBlocks[58] - (mean_of8x8_blocks[58] * mean_of8x8_blocks[58])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_59] = (uint16_t)((meanOf8x8SquaredValuesBlocks[59] - (mean_of8x8_blocks[59] * mean_of8x8_blocks[59])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_60] = (uint16_t)((meanOf8x8SquaredValuesBlocks[60] - (mean_of8x8_blocks[60] * mean_of8x8_blocks[60])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_61] = (uint16_t)((meanOf8x8SquaredValuesBlocks[61] - (mean_of8x8_blocks[61] * mean_of8x8_blocks[61])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_62] = (uint16_t)((meanOf8x8SquaredValuesBlocks[62] - (mean_of8x8_blocks[62] * mean_of8x8_blocks[62])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_63] = (uint16_t)((meanOf8x8SquaredValuesBlocks[63] - (mean_of8x8_blocks[63] * mean_of8x8_blocks[63])) >> VARIANCE_PRECISION);

    // 16x16 variances
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_0] = (uint16_t)((meanOf16x16SquaredValuesBlocks[0] - (meanOf16x16Blocks[0] * meanOf16x16Blocks[0])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_1] = (uint16_t)((meanOf16x16SquaredValuesBlocks[1] - (meanOf16x16Blocks[1] * meanOf16x16Blocks[1])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_2] = (uint16_t)((meanOf16x16SquaredValuesBlocks[2] - (meanOf16x16Blocks[2] * meanOf16x16Blocks[2])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_3] = (uint16_t)((meanOf16x16SquaredValuesBlocks[3] - (meanOf16x16Blocks[3] * meanOf16x16Blocks[3])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_4] = (uint16_t)((meanOf16x16SquaredValuesBlocks[4] - (meanOf16x16Blocks[4] * meanOf16x16Blocks[4])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_5] = (uint16_t)((meanOf16x16SquaredValuesBlocks[5] - (meanOf16x16Blocks[5] * meanOf16x16Blocks[5])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_6] = (uint16_t)((meanOf16x16SquaredValuesBlocks[6] - (meanOf16x16Blocks[6] * meanOf16x16Blocks[6])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_7] = (uint16_t)((meanOf16x16SquaredValuesBlocks[7] - (meanOf16x16Blocks[7] * meanOf16x16Blocks[7])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_8] = (uint16_t)((meanOf16x16SquaredValuesBlocks[8] - (meanOf16x16Blocks[8] * meanOf16x16Blocks[8])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_9] = (uint16_t)((meanOf16x16SquaredValuesBlocks[9] - (meanOf16x16Blocks[9] * meanOf16x16Blocks[9])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_10] = (uint16_t)((meanOf16x16SquaredValuesBlocks[10] - (meanOf16x16Blocks[10] * meanOf16x16Blocks[10])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_11] = (uint16_t)((meanOf16x16SquaredValuesBlocks[11] - (meanOf16x16Blocks[11] * meanOf16x16Blocks[11])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_12] = (uint16_t)((meanOf16x16SquaredValuesBlocks[12] - (meanOf16x16Blocks[12] * meanOf16x16Blocks[12])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_13] = (uint16_t)((meanOf16x16SquaredValuesBlocks[13] - (meanOf16x16Blocks[13] * meanOf16x16Blocks[13])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_14] = (uint16_t)((meanOf16x16SquaredValuesBlocks[14] - (meanOf16x16Blocks[14] * meanOf16x16Blocks[14])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_15] = (uint16_t)((meanOf16x16SquaredValuesBlocks[15] - (meanOf16x16Blocks[15] * meanOf16x16Blocks[15])) >> VARIANCE_PRECISION);

    // 32x32 variances
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_32x32_0] = (uint16_t)((meanOf32x32SquaredValuesBlocks[0] - (meanOf32x32Blocks[0] * meanOf32x32Blocks[0])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_32x32_1] = (uint16_t)((meanOf32x32SquaredValuesBlocks[1] - (meanOf32x32Blocks[1] * meanOf32x32Blocks[1])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_32x32_2] = (uint16_t)((meanOf32x32SquaredValuesBlocks[2] - (meanOf32x32Blocks[2] * meanOf32x32Blocks[2])) >> VARIANCE_PRECISION);
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_32x32_3] = (uint16_t)((meanOf32x32SquaredValuesBlocks[3] - (meanOf32x32Blocks[3] * meanOf32x32Blocks[3])) >> VARIANCE_PRECISION);

    // 64x64 variance
    picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_64x64] = (uint16_t)((meanOf64x64SquaredValuesBlocks - (meanOf64x64Blocks * meanOf64x64Blocks)) >> VARIANCE_PRECISION);

    return return_error;
}

EbErrorType DenoiseInputPicture(
    PictureAnalysisContext    *context_ptr,
    SequenceControlSet        *sequence_control_set_ptr,
    PictureParentControlSet   *picture_control_set_ptr,
    uint32_t                       sb_total_count,
    EbPictureBufferDesc       *input_picture_ptr,
    EbPictureBufferDesc       *denoised_picture_ptr,
    uint32_t                         picture_width_in_sb)
{
    EbErrorType return_error = EB_ErrorNone;

    uint32_t         lcuCodingOrder;
    uint32_t       sb_origin_x;
    uint32_t       sb_origin_y;
    uint16_t       verticalIdx;

    uint32_t color_format = input_picture_ptr->color_format;
    const uint16_t subsampling_x = (color_format == EB_YUV444 ? 1 : 2) - 1;
    const uint16_t subsampling_y = (color_format >= EB_YUV422 ? 1 : 2) - 1;

    //use denoised input if the source is extremly noisy
    if (picture_control_set_ptr->pic_noise_class >= PIC_NOISE_CLASS_4) {
        uint32_t inLumaOffSet = input_picture_ptr->origin_x + input_picture_ptr->origin_y      * input_picture_ptr->stride_y;
        uint32_t inChromaOffSet = (input_picture_ptr->origin_x >> subsampling_x) + (input_picture_ptr->origin_y >> subsampling_y) * input_picture_ptr->stride_cb;
        uint32_t denLumaOffSet = denoised_picture_ptr->origin_x + denoised_picture_ptr->origin_y   * denoised_picture_ptr->stride_y;
        uint32_t denChromaOffSet = (denoised_picture_ptr->origin_x >> subsampling_x) + (denoised_picture_ptr->origin_y >> subsampling_y) * denoised_picture_ptr->stride_cb;

        //filter Luma
        for (lcuCodingOrder = 0; lcuCodingOrder < sb_total_count; ++lcuCodingOrder) {
            sb_origin_x = (lcuCodingOrder % picture_width_in_sb) * sequence_control_set_ptr->sb_sz;
            sb_origin_y = (lcuCodingOrder / picture_width_in_sb) * sequence_control_set_ptr->sb_sz;

            if (sb_origin_x == 0)
                noise_extract_luma_strong(
                    input_picture_ptr,
                    denoised_picture_ptr,
                    sb_origin_y,
                    sb_origin_x);

            if (sb_origin_x + BLOCK_SIZE_64 > input_picture_ptr->width)
            {
                noise_extract_luma_strong_c(
                    input_picture_ptr,
                    denoised_picture_ptr,
                    sb_origin_y,
                    sb_origin_x);
            }
        }

        //copy Luma
        for (verticalIdx = 0; verticalIdx < input_picture_ptr->height; ++verticalIdx) {
            EB_MEMCPY(input_picture_ptr->buffer_y + inLumaOffSet + verticalIdx * input_picture_ptr->stride_y,
                denoised_picture_ptr->buffer_y + denLumaOffSet + verticalIdx * denoised_picture_ptr->stride_y,
                sizeof(uint8_t) * input_picture_ptr->width);
        }

        //copy chroma
        for (lcuCodingOrder = 0; lcuCodingOrder < sb_total_count; ++lcuCodingOrder) {
            sb_origin_x = (lcuCodingOrder % picture_width_in_sb) * sequence_control_set_ptr->sb_sz;
            sb_origin_y = (lcuCodingOrder / picture_width_in_sb) * sequence_control_set_ptr->sb_sz;

            if (sb_origin_x == 0)
                noise_extract_chroma_strong(
                    input_picture_ptr,
                    denoised_picture_ptr,
                    sb_origin_y >> subsampling_y,
                    sb_origin_x >> subsampling_x);

            if (sb_origin_x + BLOCK_SIZE_64 > input_picture_ptr->width)
            {
                noise_extract_chroma_strong_c(
                    input_picture_ptr,
                    denoised_picture_ptr,
                    sb_origin_y >> subsampling_y,
                    sb_origin_x >> subsampling_x);
            }
        }

        //copy chroma
        for (verticalIdx = 0; verticalIdx < input_picture_ptr->height >> subsampling_y; ++verticalIdx) {
            EB_MEMCPY(input_picture_ptr->buffer_cb + inChromaOffSet + verticalIdx * input_picture_ptr->stride_cb,
                denoised_picture_ptr->buffer_cb + denChromaOffSet + verticalIdx * denoised_picture_ptr->stride_cb,
                sizeof(uint8_t) * input_picture_ptr->width >> subsampling_x);

            EB_MEMCPY(input_picture_ptr->buffer_cr + inChromaOffSet + verticalIdx * input_picture_ptr->stride_cr,
                denoised_picture_ptr->buffer_cr + denChromaOffSet + verticalIdx * denoised_picture_ptr->stride_cr,
                sizeof(uint8_t) * input_picture_ptr->width >> subsampling_x);
        }
    }
    else if (picture_control_set_ptr->pic_noise_class >= PIC_NOISE_CLASS_3_1) {
        uint32_t inLumaOffSet = input_picture_ptr->origin_x + input_picture_ptr->origin_y      * input_picture_ptr->stride_y;
        uint32_t inChromaOffSet = (input_picture_ptr->origin_x >> subsampling_x) + (input_picture_ptr->origin_y >> subsampling_y) * input_picture_ptr->stride_cb;
        uint32_t denLumaOffSet = denoised_picture_ptr->origin_x + denoised_picture_ptr->origin_y   * denoised_picture_ptr->stride_y;
        uint32_t denChromaOffSet = (denoised_picture_ptr->origin_x >> subsampling_x) + (denoised_picture_ptr->origin_y >> subsampling_y) * denoised_picture_ptr->stride_cb;

        for (verticalIdx = 0; verticalIdx < input_picture_ptr->height; ++verticalIdx) {
            EB_MEMCPY(input_picture_ptr->buffer_y + inLumaOffSet + verticalIdx * input_picture_ptr->stride_y,
                denoised_picture_ptr->buffer_y + denLumaOffSet + verticalIdx * denoised_picture_ptr->stride_y,
                sizeof(uint8_t) * input_picture_ptr->width);
        }

        //copy chroma
        for (lcuCodingOrder = 0; lcuCodingOrder < sb_total_count; ++lcuCodingOrder) {
            sb_origin_x = (lcuCodingOrder % picture_width_in_sb) * sequence_control_set_ptr->sb_sz;
            sb_origin_y = (lcuCodingOrder / picture_width_in_sb) * sequence_control_set_ptr->sb_sz;

            if (sb_origin_x == 0)
                noise_extract_chroma_weak(
                    input_picture_ptr,
                    denoised_picture_ptr,
                    sb_origin_y >> subsampling_y,
                    sb_origin_x >> subsampling_x);

            if (sb_origin_x + BLOCK_SIZE_64 > input_picture_ptr->width)
            {
                noise_extract_chroma_weak_c(
                    input_picture_ptr,
                    denoised_picture_ptr,
                    sb_origin_y >> subsampling_y,
                    sb_origin_x >> subsampling_x);
            }
        }

        for (verticalIdx = 0; verticalIdx < input_picture_ptr->height >> subsampling_y; ++verticalIdx) {
            EB_MEMCPY(input_picture_ptr->buffer_cb + inChromaOffSet + verticalIdx * input_picture_ptr->stride_cb,
                denoised_picture_ptr->buffer_cb + denChromaOffSet + verticalIdx * denoised_picture_ptr->stride_cb,
                sizeof(uint8_t) * input_picture_ptr->width >> subsampling_x);

            EB_MEMCPY(input_picture_ptr->buffer_cr + inChromaOffSet + verticalIdx * input_picture_ptr->stride_cr,
                denoised_picture_ptr->buffer_cr + denChromaOffSet + verticalIdx * denoised_picture_ptr->stride_cr,
                sizeof(uint8_t) * input_picture_ptr->width >> subsampling_x);
        }
    }
    else if (context_ptr->pic_noise_variance_float >= 1.0) {
        //Luma : use filtered only for flatNoise LCUs
        for (lcuCodingOrder = 0; lcuCodingOrder < sb_total_count; ++lcuCodingOrder) {
            sb_origin_x = (lcuCodingOrder % picture_width_in_sb) * sequence_control_set_ptr->sb_sz;
            sb_origin_y = (lcuCodingOrder / picture_width_in_sb) * sequence_control_set_ptr->sb_sz;
            uint32_t  sb_height = MIN(BLOCK_SIZE_64, input_picture_ptr->height - sb_origin_y);
            uint32_t  sb_width = MIN(BLOCK_SIZE_64, input_picture_ptr->width - sb_origin_x);

            uint32_t inLumaOffSet = input_picture_ptr->origin_x + sb_origin_x + (input_picture_ptr->origin_y + sb_origin_y) * input_picture_ptr->stride_y;
            uint32_t denLumaOffSet = denoised_picture_ptr->origin_x + sb_origin_x + (denoised_picture_ptr->origin_y + sb_origin_y) * denoised_picture_ptr->stride_y;

            if (picture_control_set_ptr->sb_flat_noise_array[lcuCodingOrder] == 1) {
                for (verticalIdx = 0; verticalIdx < sb_height; ++verticalIdx) {
                    EB_MEMCPY(input_picture_ptr->buffer_y + inLumaOffSet + verticalIdx * input_picture_ptr->stride_y,
                        denoised_picture_ptr->buffer_y + denLumaOffSet + verticalIdx * denoised_picture_ptr->stride_y,
                        sizeof(uint8_t) * sb_width);
                }
            }
        }
    }

    return return_error;
}

EbErrorType DetectInputPictureNoise(
    PictureAnalysisContext    *context_ptr,
    SequenceControlSet        *sequence_control_set_ptr,
    PictureParentControlSet   *picture_control_set_ptr,
    uint32_t                       sb_total_count,
    EbPictureBufferDesc       *input_picture_ptr,
    EbPictureBufferDesc       *noise_picture_ptr,
    EbPictureBufferDesc       *denoised_picture_ptr,
    uint32_t                     picture_width_in_sb)
{
    EbErrorType                 return_error = EB_ErrorNone;
    uint32_t                    lcuCodingOrder;

    uint64_t                    picNoiseVariance;

    uint32_t                    totLcuCount, noiseTh;

    uint32_t                    sb_origin_x;
    uint32_t                    sb_origin_y;
    uint32_t                    inputLumaOriginIndex;

    picNoiseVariance = 0;
    totLcuCount = 0;

    //Variance calc for noise picture
    for (lcuCodingOrder = 0; lcuCodingOrder < sb_total_count; ++lcuCodingOrder) {
        sb_origin_x = (lcuCodingOrder % picture_width_in_sb) * sequence_control_set_ptr->sb_sz;
        sb_origin_y = (lcuCodingOrder / picture_width_in_sb) * sequence_control_set_ptr->sb_sz;
        inputLumaOriginIndex = (noise_picture_ptr->origin_y + sb_origin_y) * noise_picture_ptr->stride_y +
            noise_picture_ptr->origin_x + sb_origin_x;

        uint32_t  noiseOriginIndex = noise_picture_ptr->origin_x + sb_origin_x + noise_picture_ptr->origin_y * noise_picture_ptr->stride_y;

        if (sb_origin_x == 0)
            noise_extract_luma_weak(
                input_picture_ptr,
                denoised_picture_ptr,
                noise_picture_ptr,
                sb_origin_y,
                sb_origin_x);

        if (sb_origin_x + BLOCK_SIZE_64 > input_picture_ptr->width)
        {
            noise_extract_luma_weak_c(
                input_picture_ptr,
                denoised_picture_ptr,
                noise_picture_ptr,
                sb_origin_y,
                sb_origin_x);
        }

        //do it only for complete 64x64 blocks
        if (sb_origin_x + 64 <= input_picture_ptr->width && sb_origin_y + 64 <= input_picture_ptr->height)
        {
            uint64_t noiseBlkVar32x32[4], denoiseBlkVar32x32[4];

            uint64_t noiseBlkVar = ComputeVariance64x64(
                sequence_control_set_ptr,
                noise_picture_ptr,
                noiseOriginIndex,
                noiseBlkVar32x32);

            uint64_t noiseBlkVarTh;
            uint64_t denBlkVarTh = FLAT_MAX_VAR;
            noiseBlkVarTh = NOISE_MIN_LEVEL_M6_M7;

            picNoiseVariance += (noiseBlkVar >> 16);

            uint64_t denBlkVar = ComputeVariance64x64(
                sequence_control_set_ptr,
                denoised_picture_ptr,
                inputLumaOriginIndex,
                denoiseBlkVar32x32) >> 16;

            if (denBlkVar < denBlkVarTh && noiseBlkVar > noiseBlkVarTh)
                picture_control_set_ptr->sb_flat_noise_array[lcuCodingOrder] = 1;
            totLcuCount++;
        }
    }

    if (totLcuCount > 0) {
        context_ptr->pic_noise_variance_float = (double)picNoiseVariance / (double)totLcuCount;
        picNoiseVariance = picNoiseVariance / totLcuCount;
    }

    //the variance of a 64x64 noise area tends to be bigger for small resolutions.
    if (sequence_control_set_ptr->seq_header.max_frame_height <= 720)
        noiseTh = 25;
    else
        noiseTh = 0;

    if (picNoiseVariance >= 80 + noiseTh)
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_10;
    else if (picNoiseVariance >= 70 + noiseTh)
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_9;
    else if (picNoiseVariance >= 60 + noiseTh)
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_8;
    else if (picNoiseVariance >= 50 + noiseTh)
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_7;
    else if (picNoiseVariance >= 40 + noiseTh)
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_6;
    else if (picNoiseVariance >= 30 + noiseTh)
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_5;
    else if (picNoiseVariance >= 20 + noiseTh)
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_4;
    else if (picNoiseVariance >= 17 + noiseTh)
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_3_1;
    else if (picNoiseVariance >= 10 + noiseTh)
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_3;
    else if (picNoiseVariance >= 5 + noiseTh)
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_2;
    else
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_1;

    if (picture_control_set_ptr->pic_noise_class >= PIC_NOISE_CLASS_4)
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_3_1;

    return return_error;
}

static int32_t apply_denoise_2d(SequenceControlSet        *scs_ptr,
    PictureParentControlSet   *pcs_ptr,
    EbPictureBufferDesc *inputPicturePointer)
{
    if (eb_aom_denoise_and_model_run(pcs_ptr->denoise_and_model, inputPicturePointer,
        &pcs_ptr->frm_hdr.film_grain_params,
        scs_ptr->static_config.encoder_bit_depth > EB_8BIT)){
    }
    return 0;
}

EbErrorType denoise_estimate_film_grain(
    SequenceControlSet        *sequence_control_set_ptr,
    PictureParentControlSet   *picture_control_set_ptr)
{
    EbErrorType return_error = EB_ErrorNone;

    FrameHeader *frm_hdr = &picture_control_set_ptr->frm_hdr;

    EbPictureBufferDesc    *input_picture_ptr = picture_control_set_ptr->enhanced_picture_ptr;
    frm_hdr->film_grain_params.apply_grain = 0;

    if (sequence_control_set_ptr->film_grain_denoise_strength) {
        if (apply_denoise_2d(sequence_control_set_ptr, picture_control_set_ptr, input_picture_ptr) < 0)
            return 1;
    }

    sequence_control_set_ptr->seq_header.film_grain_params_present |= frm_hdr->film_grain_params.apply_grain;

    return return_error;  //todo: add proper error handling
}

EbErrorType FullSampleDenoise(
    PictureAnalysisContext    *context_ptr,
    SequenceControlSet        *sequence_control_set_ptr,
    PictureParentControlSet   *picture_control_set_ptr,
    uint32_t                     sb_total_count,
    EbBool                       denoise_flag,
    uint32_t                     picture_width_in_sb)
{
    EbErrorType return_error = EB_ErrorNone;

    uint32_t                     lcuCodingOrder;
    EbPictureBufferDesc    *input_picture_ptr = picture_control_set_ptr->enhanced_picture_ptr;
    EbPictureBufferDesc    *denoised_picture_ptr = context_ptr->denoised_picture_ptr;
    EbPictureBufferDesc    *noise_picture_ptr = context_ptr->noise_picture_ptr;

    //Reset the flat noise flag array to False for both RealTime/HighComplexity Modes
    for (lcuCodingOrder = 0; lcuCodingOrder < sb_total_count; ++lcuCodingOrder)
        picture_control_set_ptr->sb_flat_noise_array[lcuCodingOrder] = 0;
    picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_INV; //this init is for both REAL-TIME and BEST-QUALITY

    DetectInputPictureNoise(
        context_ptr,
        sequence_control_set_ptr,
        picture_control_set_ptr,
        sb_total_count,
        input_picture_ptr,
        noise_picture_ptr,
        denoised_picture_ptr,
        picture_width_in_sb);

    if (denoise_flag == EB_TRUE)
    {
        DenoiseInputPicture(
            context_ptr,
            sequence_control_set_ptr,
            picture_control_set_ptr,
            sb_total_count,
            input_picture_ptr,
            denoised_picture_ptr,
            picture_width_in_sb);
    }

    return return_error;
}

EbErrorType SubSampleFilterNoise(
    SequenceControlSet        *sequence_control_set_ptr,
    PictureParentControlSet   *picture_control_set_ptr,
    uint32_t                       sb_total_count,
    EbPictureBufferDesc       *input_picture_ptr,
    EbPictureBufferDesc       *noise_picture_ptr,
    EbPictureBufferDesc       *denoised_picture_ptr,
    uint32_t                         picture_width_in_sb)
{
    EbErrorType return_error = EB_ErrorNone;

    uint32_t         lcuCodingOrder;
    uint32_t       sb_origin_x;
    uint32_t       sb_origin_y;
    uint16_t        verticalIdx;

    uint32_t color_format = input_picture_ptr->color_format;
    const uint16_t subsampling_x = (color_format == EB_YUV444 ? 1 : 2) - 1;
    const uint16_t subsampling_y = (color_format >= EB_YUV422 ? 1 : 2) - 1;

    if (picture_control_set_ptr->pic_noise_class == PIC_NOISE_CLASS_3_1) {
        uint32_t inLumaOffSet = input_picture_ptr->origin_x + input_picture_ptr->origin_y      * input_picture_ptr->stride_y;
        uint32_t inChromaOffSet = (input_picture_ptr->origin_x >> subsampling_x) + (input_picture_ptr->origin_y >> subsampling_y) * input_picture_ptr->stride_cb;
        uint32_t denLumaOffSet = denoised_picture_ptr->origin_x + denoised_picture_ptr->origin_y   * denoised_picture_ptr->stride_y;
        uint32_t denChromaOffSet = (denoised_picture_ptr->origin_x >> subsampling_x) + (denoised_picture_ptr->origin_y >> subsampling_y) * denoised_picture_ptr->stride_cb;

        //filter Luma
        for (lcuCodingOrder = 0; lcuCodingOrder < sb_total_count; ++lcuCodingOrder) {
            sb_origin_x = (lcuCodingOrder % picture_width_in_sb) * sequence_control_set_ptr->sb_sz;
            sb_origin_y = (lcuCodingOrder / picture_width_in_sb) * sequence_control_set_ptr->sb_sz;

            if (sb_origin_x == 0)
                noise_extract_luma_weak(
                    input_picture_ptr,
                    denoised_picture_ptr,
                    noise_picture_ptr,
                    sb_origin_y,
                    sb_origin_x);

            if (sb_origin_x + BLOCK_SIZE_64 > input_picture_ptr->width)
            {
                noise_extract_luma_weak_c(
                    input_picture_ptr,
                    denoised_picture_ptr,
                    noise_picture_ptr,
                    sb_origin_y,
                    sb_origin_x);
            }
        }

        //copy luma
        for (verticalIdx = 0; verticalIdx < input_picture_ptr->height; ++verticalIdx) {
            EB_MEMCPY(input_picture_ptr->buffer_y + inLumaOffSet + verticalIdx * input_picture_ptr->stride_y,
                denoised_picture_ptr->buffer_y + denLumaOffSet + verticalIdx * denoised_picture_ptr->stride_y,
                sizeof(uint8_t) * input_picture_ptr->width);
        }

        //filter chroma
        for (lcuCodingOrder = 0; lcuCodingOrder < sb_total_count; ++lcuCodingOrder) {
            sb_origin_x = (lcuCodingOrder % picture_width_in_sb) * sequence_control_set_ptr->sb_sz;
            sb_origin_y = (lcuCodingOrder / picture_width_in_sb) * sequence_control_set_ptr->sb_sz;

            if (sb_origin_x == 0)
                noise_extract_chroma_weak(
                    input_picture_ptr,
                    denoised_picture_ptr,
                    sb_origin_y >> subsampling_y,
                    sb_origin_x >> subsampling_x);

            if (sb_origin_x + BLOCK_SIZE_64 > input_picture_ptr->width)
            {
                noise_extract_chroma_weak_c(
                    input_picture_ptr,
                    denoised_picture_ptr,
                    sb_origin_y >> subsampling_y,
                    sb_origin_x >> subsampling_x);
            }
        }

        //copy chroma
        for (verticalIdx = 0; verticalIdx < input_picture_ptr->height >> subsampling_y; ++verticalIdx) {
            EB_MEMCPY(input_picture_ptr->buffer_cb + inChromaOffSet + verticalIdx * input_picture_ptr->stride_cb,
                denoised_picture_ptr->buffer_cb + denChromaOffSet + verticalIdx * denoised_picture_ptr->stride_cb,
                sizeof(uint8_t) * input_picture_ptr->width >> subsampling_x);

            EB_MEMCPY(input_picture_ptr->buffer_cr + inChromaOffSet + verticalIdx * input_picture_ptr->stride_cr,
                denoised_picture_ptr->buffer_cr + denChromaOffSet + verticalIdx * denoised_picture_ptr->stride_cr,
                sizeof(uint8_t) * input_picture_ptr->width >> subsampling_x);
        }
    }
    else if (picture_control_set_ptr->pic_noise_class == PIC_NOISE_CLASS_2) {
        uint32_t newTotFN = 0;

        //for each SB ,re check the FN information for only the FNdecim ones
        for (lcuCodingOrder = 0; lcuCodingOrder < sb_total_count; ++lcuCodingOrder) {
            sb_origin_x = (lcuCodingOrder % picture_width_in_sb) * sequence_control_set_ptr->sb_sz;
            sb_origin_y = (lcuCodingOrder / picture_width_in_sb) * sequence_control_set_ptr->sb_sz;
            uint32_t  inputLumaOriginIndex = noise_picture_ptr->origin_x + sb_origin_x + (noise_picture_ptr->origin_y + sb_origin_y) * noise_picture_ptr->stride_y;
            uint32_t  noiseOriginIndex = noise_picture_ptr->origin_x + sb_origin_x + (noise_picture_ptr->origin_y * noise_picture_ptr->stride_y);

            if (sb_origin_x + 64 <= input_picture_ptr->width && sb_origin_y + 64 <= input_picture_ptr->height && picture_control_set_ptr->sb_flat_noise_array[lcuCodingOrder] == 1)
            {
                noise_extract_luma_weak_lcu(
                    input_picture_ptr,
                    denoised_picture_ptr,
                    noise_picture_ptr,
                    sb_origin_y,
                    sb_origin_x);

                if (sb_origin_x + BLOCK_SIZE_64 > input_picture_ptr->width)
                {
                    noise_extract_luma_weak_lcu_c(
                        input_picture_ptr,
                        denoised_picture_ptr,
                        noise_picture_ptr,
                        sb_origin_y,
                        sb_origin_x);
                }

                uint64_t noiseBlkVar32x32[4], denoiseBlkVar32x32[4];
                uint64_t noiseBlkVar = ComputeVariance64x64(
                    sequence_control_set_ptr,
                    noise_picture_ptr,
                    noiseOriginIndex,
                    noiseBlkVar32x32);
                uint64_t denBlkVar = ComputeVariance64x64(
                    sequence_control_set_ptr,
                    denoised_picture_ptr,
                    inputLumaOriginIndex,
                    denoiseBlkVar32x32) >> 16;

                uint64_t noiseBlkVarTh;
                uint64_t denBlkVarTh = FLAT_MAX_VAR;
                noiseBlkVarTh = NOISE_MIN_LEVEL_M6_M7;

                if (denBlkVar<denBlkVarTh && noiseBlkVar> noiseBlkVarTh) {
                    picture_control_set_ptr->sb_flat_noise_array[lcuCodingOrder] = 1;
                    //printf("POC %i (%i,%i) denBlkVar: %i  noiseBlkVar :%i\n", picture_control_set_ptr->picture_number,sb_origin_x,sb_origin_y, denBlkVar, noiseBlkVar);
                    newTotFN++;
                }
                else
                    picture_control_set_ptr->sb_flat_noise_array[lcuCodingOrder] = 0;
            }
        }

        for (lcuCodingOrder = 0; lcuCodingOrder < sb_total_count; ++lcuCodingOrder) {
            sb_origin_x = (lcuCodingOrder % picture_width_in_sb) * sequence_control_set_ptr->sb_sz;
            sb_origin_y = (lcuCodingOrder / picture_width_in_sb) * sequence_control_set_ptr->sb_sz;

            if (sb_origin_x + 64 <= input_picture_ptr->width && sb_origin_y + 64 <= input_picture_ptr->height)
            {
                //use the denoised for FN LCUs
                if (picture_control_set_ptr->sb_flat_noise_array[lcuCodingOrder] == 1) {
                    uint32_t  sb_height = MIN(BLOCK_SIZE_64, input_picture_ptr->height - sb_origin_y);
                    uint32_t  sb_width = MIN(BLOCK_SIZE_64, input_picture_ptr->width - sb_origin_x);

                    uint32_t inLumaOffSet = input_picture_ptr->origin_x + sb_origin_x + (input_picture_ptr->origin_y + sb_origin_y) * input_picture_ptr->stride_y;
                    uint32_t denLumaOffSet = denoised_picture_ptr->origin_x + sb_origin_x + (denoised_picture_ptr->origin_y + sb_origin_y) * denoised_picture_ptr->stride_y;

                    for (verticalIdx = 0; verticalIdx < sb_height; ++verticalIdx) {
                        EB_MEMCPY(input_picture_ptr->buffer_y + inLumaOffSet + verticalIdx * input_picture_ptr->stride_y,
                            denoised_picture_ptr->buffer_y + denLumaOffSet + verticalIdx * denoised_picture_ptr->stride_y,
                            sizeof(uint8_t) * sb_width);
                    }
                }
            }
        }
    }
    return return_error;
}

EbErrorType QuarterSampleDetectNoise(
    PictureAnalysisContext    *context_ptr,
    PictureParentControlSet   *picture_control_set_ptr,
    EbPictureBufferDesc       *quarter_decimated_picture_ptr,
    EbPictureBufferDesc       *noise_picture_ptr,
    EbPictureBufferDesc       *denoised_picture_ptr,
    uint32_t                         picture_width_in_sb)
{
    EbErrorType return_error = EB_ErrorNone;

    uint64_t                   picNoiseVariance;

    uint32_t                     totLcuCount, noiseTh;

    uint32_t                     blockIndex;

    picNoiseVariance = 0;
    totLcuCount = 0;

    uint16_t vert64x64Index;
    uint16_t horz64x64Index;
    uint32_t block64x64X;
    uint32_t block64x64Y;
    uint32_t vert32x32Index;
    uint32_t horz32x32Index;
    uint32_t block32x32X;
    uint32_t block32x32Y;
    uint32_t noiseOriginIndex;
    uint32_t lcuCodingOrder;

    // Loop over 64x64 blocks on the downsampled domain (each block would contain 16 LCUs on the full sampled domain)
    for (vert64x64Index = 0; vert64x64Index < (quarter_decimated_picture_ptr->height / 64); vert64x64Index++) {
        for (horz64x64Index = 0; horz64x64Index < (quarter_decimated_picture_ptr->width / 64); horz64x64Index++) {
            block64x64X = horz64x64Index * 64;
            block64x64Y = vert64x64Index * 64;

            if (block64x64X == 0)
                noise_extract_luma_weak(
                    quarter_decimated_picture_ptr,
                    denoised_picture_ptr,
                    noise_picture_ptr,
                    block64x64Y,
                    block64x64X);

            if (block64x64Y + BLOCK_SIZE_64 > quarter_decimated_picture_ptr->width)
            {
                noise_extract_luma_weak_c(
                    quarter_decimated_picture_ptr,
                    denoised_picture_ptr,
                    noise_picture_ptr,
                    block64x64Y,
                    block64x64X);
            }

            // Loop over 32x32 blocks (i.e, 64x64 blocks in full resolution)
            for (vert32x32Index = 0; vert32x32Index < 2; vert32x32Index++) {
                for (horz32x32Index = 0; horz32x32Index < 2; horz32x32Index++) {
                    block32x32X = block64x64X + horz32x32Index * 32;
                    block32x32Y = block64x64Y + vert32x32Index * 32;

                    //do it only for complete 32x32 blocks (i.e, complete 64x64 blocks in full resolution)
                    if ((block32x32X + 32 <= quarter_decimated_picture_ptr->width) && (block32x32Y + 32 <= quarter_decimated_picture_ptr->height))
                    {
                        lcuCodingOrder = ((vert64x64Index * 2) + vert32x32Index) * picture_width_in_sb + ((horz64x64Index * 2) + horz32x32Index);

                        uint64_t noiseBlkVar8x8[16], denoiseBlkVar8x8[16];

                        noiseOriginIndex = noise_picture_ptr->origin_x + block32x32X + noise_picture_ptr->origin_y * noise_picture_ptr->stride_y;

                        uint64_t noiseBlkVar = ComputeVariance32x32(
                            noise_picture_ptr,
                            noiseOriginIndex,
                            noiseBlkVar8x8);

                        picNoiseVariance += (noiseBlkVar >> 16);

                        blockIndex = (noise_picture_ptr->origin_y + block32x32Y) * noise_picture_ptr->stride_y + noise_picture_ptr->origin_x + block32x32X;

                        uint64_t denBlkVar = ComputeVariance32x32(
                            denoised_picture_ptr,
                            blockIndex,
                            denoiseBlkVar8x8) >> 16;

                        uint64_t denBlkVarDecTh;
                        denBlkVarDecTh = NOISE_MIN_LEVEL_DECIM_M6_M7;
                        if (denBlkVar < FLAT_MAX_VAR_DECIM && noiseBlkVar> denBlkVarDecTh)
                            picture_control_set_ptr->sb_flat_noise_array[lcuCodingOrder] = 1;
                        totLcuCount++;
                    }
                }
            }
        }
    }

    if (totLcuCount > 0) {
        context_ptr->pic_noise_variance_float = (double)picNoiseVariance / (double)totLcuCount;
        picNoiseVariance = picNoiseVariance / totLcuCount;
    }

    //the variance of a 64x64 noise area tends to be bigger for small resolutions.
    //if (sequence_control_set_ptr->seq_header.max_frame_height <= 720)
    //    noiseTh = 25;
    //else if (sequence_control_set_ptr->seq_header.max_frame_height <= 1080)
    //    noiseTh = 10;
    //else
    noiseTh = 0;

    //look for extreme noise or big enough flat noisy area to be denoised.
    if (picNoiseVariance > 60)
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_3_1; //Noise+Edge information is too big, so may be this is all noise (action: frame based denoising)
    else if (picNoiseVariance >= 10 + noiseTh)
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_3;   //Noise+Edge information is big enough, so there is no big enough flat noisy area (action : no denoising)
    else if (picNoiseVariance >= 5 + noiseTh)
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_2;   //Noise+Edge information is relatively small, so there might be a big enough flat noisy area(action : denoising only for FN blocks)
    else
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_1;   //Noise+Edge information is very small, so no noise nor edge area (action : no denoising)

    return return_error;
}

EbErrorType SubSampleDetectNoise(
    PictureAnalysisContext    *context_ptr,
    SequenceControlSet        *sequence_control_set_ptr,
    PictureParentControlSet   *picture_control_set_ptr,
    EbPictureBufferDesc       *sixteenth_decimated_picture_ptr,
    EbPictureBufferDesc       *noise_picture_ptr,
    EbPictureBufferDesc       *denoised_picture_ptr,
    uint32_t                         picture_width_in_sb)
{
    EbErrorType return_error = EB_ErrorNone;

    uint64_t                   picNoiseVariance;

    uint32_t                     totLcuCount, noiseTh;

    uint32_t                     blockIndex;

    picNoiseVariance = 0;
    totLcuCount = 0;

    uint16_t vert64x64Index;
    uint16_t horz64x64Index;
    uint32_t block64x64X;
    uint32_t block64x64Y;
    uint32_t vert16x16Index;
    uint32_t horz16x16Index;
    uint32_t block16x16X;
    uint32_t block16x16Y;
    uint32_t noiseOriginIndex;
    uint32_t lcuCodingOrder;

    // Loop over 64x64 blocks on the downsampled domain (each block would contain 16 LCUs on the full sampled domain)
    for (vert64x64Index = 0; vert64x64Index < (sixteenth_decimated_picture_ptr->height / 64); vert64x64Index++) {
        for (horz64x64Index = 0; horz64x64Index < (sixteenth_decimated_picture_ptr->width / 64); horz64x64Index++) {
            block64x64X = horz64x64Index * 64;
            block64x64Y = vert64x64Index * 64;

            if (block64x64X == 0)
                noise_extract_luma_weak(
                    sixteenth_decimated_picture_ptr,
                    denoised_picture_ptr,
                    noise_picture_ptr,
                    block64x64Y,
                    block64x64X);

            if (block64x64Y + BLOCK_SIZE_64 > sixteenth_decimated_picture_ptr->width)
            {
                noise_extract_luma_weak_c(
                    sixteenth_decimated_picture_ptr,
                    denoised_picture_ptr,
                    noise_picture_ptr,
                    block64x64Y,
                    block64x64X);
            }

            // Loop over 16x16 blocks (i.e, 64x64 blocks in full resolution)
            for (vert16x16Index = 0; vert16x16Index < 4; vert16x16Index++) {
                for (horz16x16Index = 0; horz16x16Index < 4; horz16x16Index++) {
                    block16x16X = block64x64X + horz16x16Index * 16;
                    block16x16Y = block64x64Y + vert16x16Index * 16;

                    //do it only for complete 16x16 blocks (i.e, complete 64x64 blocks in full resolution)
                    if (block16x16X + 16 <= sixteenth_decimated_picture_ptr->width && block16x16Y + 16 <= sixteenth_decimated_picture_ptr->height)
                    {
                        lcuCodingOrder = ((vert64x64Index * 4) + vert16x16Index) * picture_width_in_sb + ((horz64x64Index * 4) + horz16x16Index);

                        uint64_t noiseBlkVar8x8[4], denoiseBlkVar8x8[4];

                        noiseOriginIndex = noise_picture_ptr->origin_x + block16x16X + noise_picture_ptr->origin_y * noise_picture_ptr->stride_y;

                        uint64_t noiseBlkVar = ComputeVariance16x16(
                            noise_picture_ptr,
                            noiseOriginIndex,
                            noiseBlkVar8x8);

                        picNoiseVariance += (noiseBlkVar >> 16);

                        blockIndex = (noise_picture_ptr->origin_y + block16x16Y) * noise_picture_ptr->stride_y + noise_picture_ptr->origin_x + block16x16X;

                        uint64_t denBlkVar = ComputeVariance16x16(
                            denoised_picture_ptr,
                            blockIndex,
                            denoiseBlkVar8x8) >> 16;

                        uint64_t  noiseBlkVarDecTh;
                        uint64_t denBlkVarDecTh = FLAT_MAX_VAR_DECIM;

                        noiseBlkVarDecTh = NOISE_MIN_LEVEL_DECIM_M6_M7;
                        if (denBlkVar < denBlkVarDecTh && noiseBlkVar> noiseBlkVarDecTh)
                            picture_control_set_ptr->sb_flat_noise_array[lcuCodingOrder] = 1;
                        totLcuCount++;
                    }
                }
            }
        }
    }

    if (totLcuCount > 0) {
        context_ptr->pic_noise_variance_float = (double)picNoiseVariance / (double)totLcuCount;
        picNoiseVariance = picNoiseVariance / totLcuCount;
    }

    //the variance of a 64x64 noise area tends to be bigger for small resolutions.
    if (sequence_control_set_ptr->seq_header.max_frame_height <= 720)
        noiseTh = 25;
    else if (sequence_control_set_ptr->seq_header.max_frame_height <= 1080)
        noiseTh = 10;
    else
        noiseTh = 0;

    //look for extreme noise or big enough flat noisy area to be denoised.
    if (picNoiseVariance >= 55 + noiseTh)
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_3_1; //Noise+Edge information is too big, so may be this is all noise (action: frame based denoising)
    else if (picNoiseVariance >= 10 + noiseTh)
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_3;   //Noise+Edge information is big enough, so there is no big enough flat noisy area (action : no denoising)
    else if (picNoiseVariance >= 5 + noiseTh)
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_2;   //Noise+Edge information is relatively small, so there might be a big enough flat noisy area(action : denoising only for FN blocks)
    else
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_1;   //Noise+Edge information is very small, so no noise nor edge area (action : no denoising)

    return return_error;
}

EbErrorType QuarterSampleDenoise(
    PictureAnalysisContext    *context_ptr,
    SequenceControlSet        *sequence_control_set_ptr,
    PictureParentControlSet   *picture_control_set_ptr,
    EbPictureBufferDesc        *quarter_decimated_picture_ptr,
    uint32_t                       sb_total_count,
    EbBool                      denoise_flag,
    uint32_t                         picture_width_in_sb)
{
    EbErrorType return_error = EB_ErrorNone;

    uint32_t                     lcuCodingOrder;
    EbPictureBufferDesc    *input_picture_ptr = picture_control_set_ptr->enhanced_picture_ptr;
    EbPictureBufferDesc    *denoised_picture_ptr = context_ptr->denoised_picture_ptr;
    EbPictureBufferDesc    *noise_picture_ptr = context_ptr->noise_picture_ptr;

    //Reset the flat noise flag array to False for both RealTime/HighComplexity Modes
    for (lcuCodingOrder = 0; lcuCodingOrder < sb_total_count; ++lcuCodingOrder)
        picture_control_set_ptr->sb_flat_noise_array[lcuCodingOrder] = 0;
    picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_INV; //this init is for both REAL-TIME and BEST-QUALITY

    decimation_2d(
        &input_picture_ptr->buffer_y[input_picture_ptr->origin_x + input_picture_ptr->origin_y * input_picture_ptr->stride_y],
        input_picture_ptr->stride_y,
        input_picture_ptr->width,
        input_picture_ptr->height,
        &quarter_decimated_picture_ptr->buffer_y[quarter_decimated_picture_ptr->origin_x + (quarter_decimated_picture_ptr->origin_y * quarter_decimated_picture_ptr->stride_y)],
        quarter_decimated_picture_ptr->stride_y,
        2);

    QuarterSampleDetectNoise(
        context_ptr,
        picture_control_set_ptr,
        quarter_decimated_picture_ptr,
        noise_picture_ptr,
        denoised_picture_ptr,
        picture_width_in_sb);

    if (denoise_flag == EB_TRUE) {
        // Turn OFF the de-noiser for Class 2 at QP=29 and lower (for Fixed_QP) and at the target rate of 14Mbps and higher (for RC=ON)
        if ((picture_control_set_ptr->pic_noise_class == PIC_NOISE_CLASS_3_1) ||
            ((picture_control_set_ptr->pic_noise_class == PIC_NOISE_CLASS_2) && ((sequence_control_set_ptr->static_config.rate_control_mode == 0 && sequence_control_set_ptr->qp > DENOISER_QP_TH) || (sequence_control_set_ptr->static_config.rate_control_mode != 0 && sequence_control_set_ptr->static_config.target_bit_rate < DENOISER_BITRATE_TH)))) {
            SubSampleFilterNoise(
                sequence_control_set_ptr,
                picture_control_set_ptr,
                sb_total_count,
                input_picture_ptr,
                noise_picture_ptr,
                denoised_picture_ptr,
                picture_width_in_sb);
        }
    }

    return return_error;
}

EbErrorType SubSampleDenoise(
    PictureAnalysisContext    *context_ptr,
    SequenceControlSet        *sequence_control_set_ptr,
    PictureParentControlSet   *picture_control_set_ptr,
    EbPictureBufferDesc        *sixteenth_decimated_picture_ptr,
    uint32_t                       sb_total_count,
    EbBool                      denoise_flag,
    uint32_t                         picture_width_in_sb)
{
    EbErrorType return_error = EB_ErrorNone;

    uint32_t                     lcuCodingOrder;
    EbPictureBufferDesc    *input_picture_ptr = picture_control_set_ptr->enhanced_picture_ptr;
    EbPictureBufferDesc    *denoised_picture_ptr = context_ptr->denoised_picture_ptr;
    EbPictureBufferDesc    *noise_picture_ptr = context_ptr->noise_picture_ptr;

    //Reset the flat noise flag array to False for both RealTime/HighComplexity Modes
    for (lcuCodingOrder = 0; lcuCodingOrder < sb_total_count; ++lcuCodingOrder)
        picture_control_set_ptr->sb_flat_noise_array[lcuCodingOrder] = 0;
    picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_INV; //this init is for both REAL-TIME and BEST-QUALITY

    decimation_2d(
        &input_picture_ptr->buffer_y[input_picture_ptr->origin_x + input_picture_ptr->origin_y * input_picture_ptr->stride_y],
        input_picture_ptr->stride_y,
        input_picture_ptr->width,
        input_picture_ptr->height,
        &sixteenth_decimated_picture_ptr->buffer_y[sixteenth_decimated_picture_ptr->origin_x + (sixteenth_decimated_picture_ptr->origin_y * sixteenth_decimated_picture_ptr->stride_y)],
        sixteenth_decimated_picture_ptr->stride_y,
        4);

    SubSampleDetectNoise(
        context_ptr,
        sequence_control_set_ptr,
        picture_control_set_ptr,
        sixteenth_decimated_picture_ptr,
        noise_picture_ptr,
        denoised_picture_ptr,
        picture_width_in_sb);

    if (denoise_flag == EB_TRUE) {
        // Turn OFF the de-noiser for Class 2 at QP=29 and lower (for Fixed_QP) and at the target rate of 14Mbps and higher (for RC=ON)
        if ((picture_control_set_ptr->pic_noise_class == PIC_NOISE_CLASS_3_1) ||
            ((picture_control_set_ptr->pic_noise_class == PIC_NOISE_CLASS_2) && ((sequence_control_set_ptr->static_config.rate_control_mode == 0 && sequence_control_set_ptr->qp > DENOISER_QP_TH) || (sequence_control_set_ptr->static_config.rate_control_mode != 0 && sequence_control_set_ptr->static_config.target_bit_rate < DENOISER_BITRATE_TH)))) {
            SubSampleFilterNoise(
                sequence_control_set_ptr,
                picture_control_set_ptr,
                sb_total_count,
                input_picture_ptr,
                noise_picture_ptr,
                denoised_picture_ptr,
                picture_width_in_sb);
        }
    }

    return return_error;
}

/************************************************
 * Set Picture Parameters based on input configuration
 ** Setting Number of regions per resolution
 ** Setting width and height for subpicture and when picture scan type is 1
 ************************************************/
void SetPictureParametersForStatisticsGathering(
    SequenceControlSet            *sequence_control_set_ptr
)
{
    sequence_control_set_ptr->picture_analysis_number_of_regions_per_width = HIGHER_THAN_CLASS_1_REGION_SPLIT_PER_WIDTH;
    sequence_control_set_ptr->picture_analysis_number_of_regions_per_height = HIGHER_THAN_CLASS_1_REGION_SPLIT_PER_HEIGHT;

    return;
}
/************************************************
 * Picture Pre Processing Operations *
 *** A function that groups all of the Pre proceesing
 * operations performed on the input picture
 *** Operations included at this point:
 ***** Borders preprocessing
 ***** Denoising
 ************************************************/
void PicturePreProcessingOperations(
    PictureParentControlSet       *picture_control_set_ptr,
    SequenceControlSet            *sequence_control_set_ptr,
    uint32_t                       sb_total_count)
{
    if (sequence_control_set_ptr->film_grain_denoise_strength) {
        denoise_estimate_film_grain(
            sequence_control_set_ptr,
            picture_control_set_ptr);
    }
    else {
        //Reset the flat noise flag array to False for both RealTime/HighComplexity Modes
        for (uint32_t lcuCodingOrder = 0; lcuCodingOrder < sb_total_count; ++lcuCodingOrder)
            picture_control_set_ptr->sb_flat_noise_array[lcuCodingOrder] = 0;
        picture_control_set_ptr->pic_noise_class = PIC_NOISE_CLASS_INV; //this init is for both REAL-TIME and BEST-QUALITY
    }
    return;
}

/**************************************************************
* Generate picture histogram bins for YUV pixel intensity *
* Calculation is done on a region based (Set previously, resolution dependent)
**************************************************************/
void SubSampleLumaGeneratePixelIntensityHistogramBins(
    SequenceControlSet            *sequence_control_set_ptr,
    PictureParentControlSet       *picture_control_set_ptr,
    EbPictureBufferDesc           *input_picture_ptr,
    uint64_t                          *sumAverageIntensityTotalRegionsLuma) {
    uint32_t                          regionWidth;
    uint32_t                          regionHeight;
    uint32_t                          regionWidthOffset;
    uint32_t                          regionHeightOffset;
    uint32_t                          regionInPictureWidthIndex;
    uint32_t                          regionInPictureHeightIndex;
    uint32_t                            histogramBin;
    uint64_t                          sum;

    regionWidth = input_picture_ptr->width / sequence_control_set_ptr->picture_analysis_number_of_regions_per_width;
    regionHeight = input_picture_ptr->height / sequence_control_set_ptr->picture_analysis_number_of_regions_per_height;

    // Loop over regions inside the picture
    for (regionInPictureWidthIndex = 0; regionInPictureWidthIndex < sequence_control_set_ptr->picture_analysis_number_of_regions_per_width; regionInPictureWidthIndex++) {  // loop over horizontal regions
        for (regionInPictureHeightIndex = 0; regionInPictureHeightIndex < sequence_control_set_ptr->picture_analysis_number_of_regions_per_height; regionInPictureHeightIndex++) { // loop over vertical regions

            // Initialize bins to 1
            initialize_buffer_32bits(picture_control_set_ptr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][0], 64, 0, 1);

            regionWidthOffset = (regionInPictureWidthIndex == sequence_control_set_ptr->picture_analysis_number_of_regions_per_width - 1) ?
                input_picture_ptr->width - (sequence_control_set_ptr->picture_analysis_number_of_regions_per_width * regionWidth) :
                0;

            regionHeightOffset = (regionInPictureHeightIndex == sequence_control_set_ptr->picture_analysis_number_of_regions_per_height - 1) ?
                input_picture_ptr->height - (sequence_control_set_ptr->picture_analysis_number_of_regions_per_height * regionHeight) :
                0;

            // Y Histogram
            CalculateHistogram(
                &input_picture_ptr->buffer_y[(input_picture_ptr->origin_x + regionInPictureWidthIndex * regionWidth) + ((input_picture_ptr->origin_y + regionInPictureHeightIndex * regionHeight) * input_picture_ptr->stride_y)],
                regionWidth + regionWidthOffset,
                regionHeight + regionHeightOffset,
                input_picture_ptr->stride_y,
                1,
                picture_control_set_ptr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][0],
                &sum);

            picture_control_set_ptr->average_intensity_per_region[regionInPictureWidthIndex][regionInPictureHeightIndex][0] = (uint8_t)((sum + (((regionWidth + regionWidthOffset)*(regionHeight + regionHeightOffset)) >> 1)) / ((regionWidth + regionWidthOffset)*(regionHeight + regionHeightOffset)));
            (*sumAverageIntensityTotalRegionsLuma) += (sum << 4);
            for (histogramBin = 0; histogramBin < HISTOGRAM_NUMBER_OF_BINS; histogramBin++) { // Loop over the histogram bins
                picture_control_set_ptr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][0][histogramBin] =
                    picture_control_set_ptr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][0][histogramBin] << 4;
            }
        }
    }

    return;
}

void SubSampleChromaGeneratePixelIntensityHistogramBins(
    SequenceControlSet            *sequence_control_set_ptr,
    PictureParentControlSet       *picture_control_set_ptr,
    EbPictureBufferDesc           *input_picture_ptr,
    uint64_t                          *sumAverageIntensityTotalRegionsCb,
    uint64_t                          *sumAverageIntensityTotalRegionsCr) {
    uint64_t                          sum;
    uint32_t                          regionWidth;
    uint32_t                          regionHeight;
    uint32_t                          regionWidthOffset;
    uint32_t                          regionHeightOffset;
    uint32_t                          regionInPictureWidthIndex;
    uint32_t                          regionInPictureHeightIndex;

    uint16_t                          histogramBin;
    uint8_t                           decim_step = 4;

    regionWidth = input_picture_ptr->width / sequence_control_set_ptr->picture_analysis_number_of_regions_per_width;
    regionHeight = input_picture_ptr->height / sequence_control_set_ptr->picture_analysis_number_of_regions_per_height;

    // Loop over regions inside the picture
    for (regionInPictureWidthIndex = 0; regionInPictureWidthIndex < sequence_control_set_ptr->picture_analysis_number_of_regions_per_width; regionInPictureWidthIndex++) {  // loop over horizontal regions
        for (regionInPictureHeightIndex = 0; regionInPictureHeightIndex < sequence_control_set_ptr->picture_analysis_number_of_regions_per_height; regionInPictureHeightIndex++) { // loop over vertical regions

            // Initialize bins to 1
            initialize_buffer_32bits(picture_control_set_ptr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][1], 64, 0, 1);
            initialize_buffer_32bits(picture_control_set_ptr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][2], 64, 0, 1);

            regionWidthOffset = (regionInPictureWidthIndex == sequence_control_set_ptr->picture_analysis_number_of_regions_per_width - 1) ?
                input_picture_ptr->width - (sequence_control_set_ptr->picture_analysis_number_of_regions_per_width * regionWidth) :
                0;

            regionHeightOffset = (regionInPictureHeightIndex == sequence_control_set_ptr->picture_analysis_number_of_regions_per_height - 1) ?
                input_picture_ptr->height - (sequence_control_set_ptr->picture_analysis_number_of_regions_per_height * regionHeight) :
                0;

            // U Histogram
            CalculateHistogram(
                &input_picture_ptr->buffer_cb[((input_picture_ptr->origin_x + regionInPictureWidthIndex * regionWidth) >> 1) + (((input_picture_ptr->origin_y + regionInPictureHeightIndex * regionHeight) >> 1) * input_picture_ptr->stride_cb)],
                (regionWidth + regionWidthOffset) >> 1,
                (regionHeight + regionHeightOffset) >> 1,
                input_picture_ptr->stride_cb,
                decim_step,
                picture_control_set_ptr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][1],
                &sum);

            sum = (sum << decim_step);
            *sumAverageIntensityTotalRegionsCb += sum;
            picture_control_set_ptr->average_intensity_per_region[regionInPictureWidthIndex][regionInPictureHeightIndex][1] = (uint8_t)((sum + (((regionWidth + regionWidthOffset)*(regionHeight + regionHeightOffset)) >> 3)) / (((regionWidth + regionWidthOffset)*(regionHeight + regionHeightOffset)) >> 2));

            for (histogramBin = 0; histogramBin < HISTOGRAM_NUMBER_OF_BINS; histogramBin++) { // Loop over the histogram bins
                picture_control_set_ptr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][1][histogramBin] =
                    picture_control_set_ptr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][1][histogramBin] << decim_step;
            }

            // V Histogram
            CalculateHistogram(
                &input_picture_ptr->buffer_cr[((input_picture_ptr->origin_x + regionInPictureWidthIndex * regionWidth) >> 1) + (((input_picture_ptr->origin_y + regionInPictureHeightIndex * regionHeight) >> 1) * input_picture_ptr->stride_cr)],
                (regionWidth + regionWidthOffset) >> 1,
                (regionHeight + regionHeightOffset) >> 1,
                input_picture_ptr->stride_cr,
                decim_step,
                picture_control_set_ptr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][2],
                &sum);

            sum = (sum << decim_step);
            *sumAverageIntensityTotalRegionsCr += sum;
            picture_control_set_ptr->average_intensity_per_region[regionInPictureWidthIndex][regionInPictureHeightIndex][2] = (uint8_t)((sum + (((regionWidth + regionWidthOffset)*(regionHeight + regionHeightOffset)) >> 3)) / (((regionWidth + regionWidthOffset)*(regionHeight + regionHeightOffset)) >> 2));

            for (histogramBin = 0; histogramBin < HISTOGRAM_NUMBER_OF_BINS; histogramBin++) { // Loop over the histogram bins
                picture_control_set_ptr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][2][histogramBin] =
                    picture_control_set_ptr->picture_histogram[regionInPictureWidthIndex][regionInPictureHeightIndex][2][histogramBin] << decim_step;
            }
        }
    }
    return;
}

void EdgeDetectionMeanLumaChroma16x16(
    SequenceControlSet        *sequence_control_set_ptr,
    PictureParentControlSet   *picture_control_set_ptr,
    uint32_t                       totalLcuCount)
{
    uint32_t               sb_index;

    uint32_t maxGrad = 1;

    // The values are calculated for every 4th frame
    if ((picture_control_set_ptr->picture_number & 3) == 0) {
        for (sb_index = 0; sb_index < totalLcuCount; sb_index++) {
            SbStat *sb_stat_ptr = &picture_control_set_ptr->sb_stat_array[sb_index];

            EB_MEMSET(sb_stat_ptr, 0, sizeof(SbStat));
            SbParams     *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];
            if (sb_params->potential_logo_sb &&sb_params->is_complete_sb)

            {
                uint8_t *y_mean_ptr = picture_control_set_ptr->y_mean[sb_index];
                uint8_t *cr_mean_ptr = picture_control_set_ptr->crMean[sb_index];
                uint8_t *cb_mean_ptr = picture_control_set_ptr->cbMean[sb_index];

                uint8_t rasterScanCuIndex;

                for (rasterScanCuIndex = RASTER_SCAN_CU_INDEX_16x16_0; rasterScanCuIndex <= RASTER_SCAN_CU_INDEX_16x16_15; rasterScanCuIndex++) {
                    uint8_t cu_index = rasterScanCuIndex - 5;
                    uint8_t x = cu_index & 3;
                    uint8_t y = (cu_index >> 2);
                    int32_t gradx = 0;
                    int32_t grady = 0;
                    int32_t nbcompx = 0;
                    int32_t nbcompy = 0;
                    if (x != 0)
                    {
                        gradx += ABS((int32_t)(y_mean_ptr[rasterScanCuIndex]) - (int32_t)(y_mean_ptr[rasterScanCuIndex - 1]));
                        gradx += ABS((int32_t)(cr_mean_ptr[rasterScanCuIndex]) - (int32_t)(cr_mean_ptr[rasterScanCuIndex - 1]));
                        gradx += ABS((int32_t)(cb_mean_ptr[rasterScanCuIndex]) - (int32_t)(cb_mean_ptr[rasterScanCuIndex - 1]));
                        nbcompx++;
                    }
                    if (x != 3)
                    {
                        gradx += ABS((int32_t)(y_mean_ptr[rasterScanCuIndex + 1]) - (int32_t)(y_mean_ptr[rasterScanCuIndex]));
                        gradx += ABS((int32_t)(cr_mean_ptr[rasterScanCuIndex + 1]) - (int32_t)(cr_mean_ptr[rasterScanCuIndex]));
                        gradx += ABS((int32_t)(cb_mean_ptr[rasterScanCuIndex + 1]) - (int32_t)(cb_mean_ptr[rasterScanCuIndex]));
                        nbcompx++;
                    }
                    gradx = gradx / nbcompx;

                    if (y != 0)
                    {
                        grady += ABS((int32_t)(y_mean_ptr[rasterScanCuIndex]) - (int32_t)(y_mean_ptr[rasterScanCuIndex - 4]));
                        grady += ABS((int32_t)(cr_mean_ptr[rasterScanCuIndex]) - (int32_t)(cr_mean_ptr[rasterScanCuIndex - 4]));
                        grady += ABS((int32_t)(cb_mean_ptr[rasterScanCuIndex]) - (int32_t)(cb_mean_ptr[rasterScanCuIndex - 4]));
                        nbcompy++;
                    }
                    if (y != 3)
                    {
                        grady += ABS((int32_t)(y_mean_ptr[rasterScanCuIndex + 4]) - (int32_t)(y_mean_ptr[rasterScanCuIndex]));
                        grady += ABS((int32_t)(cr_mean_ptr[rasterScanCuIndex + 4]) - (int32_t)(cr_mean_ptr[rasterScanCuIndex]));
                        grady += ABS((int32_t)(cb_mean_ptr[rasterScanCuIndex + 4]) - (int32_t)(cb_mean_ptr[rasterScanCuIndex]));

                        nbcompy++;
                    }

                    grady = grady / nbcompy;
                    sb_stat_ptr->cu_stat_array[rasterScanCuIndex].grad = (uint32_t)ABS(gradx) + ABS(grady);
                    if (sb_stat_ptr->cu_stat_array[rasterScanCuIndex].grad > maxGrad)
                        maxGrad = sb_stat_ptr->cu_stat_array[rasterScanCuIndex].grad;
                }
            }
        }

        for (sb_index = 0; sb_index < totalLcuCount; sb_index++) {
            SbParams *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];
            if (sb_params->potential_logo_sb &&sb_params->is_complete_sb) {
                SbStat *sb_stat_ptr = &picture_control_set_ptr->sb_stat_array[sb_index];

                uint32_t rasterScanCuIndex;
                for (rasterScanCuIndex = RASTER_SCAN_CU_INDEX_16x16_0; rasterScanCuIndex <= RASTER_SCAN_CU_INDEX_16x16_15; rasterScanCuIndex++)
                    sb_stat_ptr->cu_stat_array[rasterScanCuIndex].edge_cu = (uint16_t)MIN(((sb_stat_ptr->cu_stat_array[rasterScanCuIndex].grad * (255 * 3)) / maxGrad), 255) < 30 ? 0 : 1;
            }
        }
    }
    else {
        for (sb_index = 0; sb_index < totalLcuCount; sb_index++) {
            SbStat *sb_stat_ptr = &picture_control_set_ptr->sb_stat_array[sb_index];

            EB_MEMSET(sb_stat_ptr, 0, sizeof(SbStat));
        }
    }
}

/******************************************************
* Edge map derivation
******************************************************/
void EdgeDetection(
    SequenceControlSet            *sequence_control_set_ptr,
    PictureParentControlSet       *picture_control_set_ptr)
{
    uint16_t  *variancePtr;
    uint32_t sb_total_count = picture_control_set_ptr->sb_total_count;
    uint64_t thrsldLevel0 = (picture_control_set_ptr->pic_avg_variance * 70) / 100;
    uint8_t  *meanPtr;
    uint32_t picture_width_in_sb = (sequence_control_set_ptr->seq_header.max_frame_width + sequence_control_set_ptr->sb_sz - 1) / sequence_control_set_ptr->sb_sz;
    uint32_t picture_height_in_sb = (sequence_control_set_ptr->seq_header.max_frame_height + sequence_control_set_ptr->sb_sz - 1) / sequence_control_set_ptr->sb_sz;
    uint32_t neighbourLcuIndex = 0;
    uint64_t similarityCount = 0;
    uint64_t similarityCount0 = 0;
    uint64_t similarityCount1 = 0;
    uint64_t similarityCount2 = 0;
    uint64_t similarityCount3 = 0;
    uint32_t sb_x = 0;
    uint32_t sb_y = 0;
    uint32_t sb_index;
    EbBool highVarianceLucFlag;

    uint32_t rasterScanCuIndex = 0;
    uint32_t numberOfEdgeLcu = 0;
    EbBool highIntensityLcuFlag;

    uint64_t neighbourLcuMean;
    int32_t i, j;

    uint8_t highIntensityTh = 180;
    uint8_t lowIntensityTh = 120;
    uint8_t highIntensityTh1 = 200;
    uint8_t veryLowIntensityTh = 20;

    for (sb_index = 0; sb_index < sb_total_count; ++sb_index) {
        sb_x = sb_index % picture_width_in_sb;
        sb_y = sb_index / picture_width_in_sb;

        EdgeLcuResults *edge_results_ptr = picture_control_set_ptr->edge_results_ptr;
        picture_control_set_ptr->edge_results_ptr[sb_index].edge_block_num = 0;
        picture_control_set_ptr->edge_results_ptr[sb_index].isolated_high_intensity_sb = 0;
        picture_control_set_ptr->sharp_edge_sb_flag[sb_index] = 0;

        if (sb_x > 0 && sb_x < (uint32_t)(picture_width_in_sb - 1) && sb_y >  0 && sb_y < (uint32_t)(picture_height_in_sb - 1)) {
            variancePtr = picture_control_set_ptr->variance[sb_index];
            meanPtr = picture_control_set_ptr->y_mean[sb_index];

            similarityCount = 0;

            highVarianceLucFlag =
                (variancePtr[RASTER_SCAN_CU_INDEX_64x64] > thrsldLevel0) ? EB_TRUE : EB_FALSE;
            edge_results_ptr[sb_index].edge_block_num = highVarianceLucFlag;
            if (variancePtr[0] > highIntensityTh1) {
                uint8_t sharpEdge = 0;
                for (rasterScanCuIndex = RASTER_SCAN_CU_INDEX_16x16_0; rasterScanCuIndex <= RASTER_SCAN_CU_INDEX_16x16_15; rasterScanCuIndex++)
                    sharpEdge = (variancePtr[rasterScanCuIndex] < veryLowIntensityTh) ? sharpEdge + 1 : sharpEdge;
                if (sharpEdge > 4)
                    picture_control_set_ptr->sharp_edge_sb_flag[sb_index] = 1;
            }

            if (sb_x > 3 && sb_x < (uint32_t)(picture_width_in_sb - 4) && sb_y >  3 && sb_y < (uint32_t)(picture_height_in_sb - 4)) {
                highIntensityLcuFlag =
                    (meanPtr[RASTER_SCAN_CU_INDEX_64x64] > highIntensityTh) ? EB_TRUE : EB_FALSE;

                if (highIntensityLcuFlag) {
                    neighbourLcuIndex = sb_index - 1;
                    neighbourLcuMean = picture_control_set_ptr->y_mean[neighbourLcuIndex][RASTER_SCAN_CU_INDEX_64x64];

                    similarityCount0 = (neighbourLcuMean < lowIntensityTh) ? 1 : 0;

                    neighbourLcuIndex = sb_index + 1;

                    neighbourLcuMean = picture_control_set_ptr->y_mean[neighbourLcuIndex][RASTER_SCAN_CU_INDEX_64x64];
                    similarityCount1 = (neighbourLcuMean < lowIntensityTh) ? 1 : 0;

                    neighbourLcuIndex = sb_index - picture_width_in_sb;
                    neighbourLcuMean = picture_control_set_ptr->y_mean[neighbourLcuIndex][RASTER_SCAN_CU_INDEX_64x64];
                    similarityCount2 = (neighbourLcuMean < lowIntensityTh) ? 1 : 0;

                    neighbourLcuIndex = sb_index + picture_width_in_sb;
                    neighbourLcuMean = picture_control_set_ptr->y_mean[neighbourLcuIndex][RASTER_SCAN_CU_INDEX_64x64];
                    similarityCount3 = (neighbourLcuMean < lowIntensityTh) ? 1 : 0;

                    similarityCount = similarityCount0 + similarityCount1 + similarityCount2 + similarityCount3;

                    if (similarityCount > 0) {
                        for (i = -4; i < 5; i++) {
                            for (j = -4; j < 5; j++) {
                                neighbourLcuIndex = sb_index + (i * picture_width_in_sb) + j;
                                picture_control_set_ptr->edge_results_ptr[neighbourLcuIndex].isolated_high_intensity_sb = 1;
                            }
                        }
                    }
                }
            }

            if (highVarianceLucFlag)
                numberOfEdgeLcu += edge_results_ptr[sb_index].edge_block_num;
        }
    }
    return;
}

/******************************************************
* Calculate the variance of variance to determine Homogeneous regions. Note: Variance calculation should be on.
******************************************************/
void DetermineHomogeneousRegionInPicture(
    SequenceControlSet            *sequence_control_set_ptr,
    PictureParentControlSet       *picture_control_set_ptr)
{
    uint16_t  *variancePtr;
    uint32_t sb_index;
    uint64_t nullVarCnt = 0;
    uint64_t veryLowVarCnt = 0;
    uint64_t varLcuCnt = 0;
    uint32_t sb_total_count = picture_control_set_ptr->sb_total_count;

    for (sb_index = 0; sb_index < sb_total_count; ++sb_index) {
        SbParams sb_params = sequence_control_set_ptr->sb_params_array[sb_index];
        variancePtr = picture_control_set_ptr->variance[sb_index];

        if (sb_params.is_complete_sb) {
            nullVarCnt += (variancePtr[ME_TIER_ZERO_PU_64x64] == 0) ? 1 : 0;

            varLcuCnt++;

            veryLowVarCnt += ((variancePtr[ME_TIER_ZERO_PU_64x64]) < LCU_LOW_VAR_TH) ? 1 : 0;
        }
    }
    picture_control_set_ptr->very_low_var_pic_flag = EB_FALSE;
    if ((varLcuCnt > 0) && (((veryLowVarCnt * 100) / varLcuCnt) > PIC_LOW_VAR_PERCENTAGE_TH))
        picture_control_set_ptr->very_low_var_pic_flag = EB_TRUE;
    picture_control_set_ptr->logo_pic_flag = EB_FALSE;
    if ((varLcuCnt > 0) && (((veryLowVarCnt * 100) / varLcuCnt) > 80))
        picture_control_set_ptr->logo_pic_flag = EB_TRUE;
    return;
}
/************************************************
 * ComputePictureSpatialStatistics
 ** Compute Block Variance
 ** Compute Picture Variance
 ** Compute Block Mean for all blocks in the picture
 ************************************************/
void ComputePictureSpatialStatistics(
    SequenceControlSet            *sequence_control_set_ptr,
    PictureParentControlSet       *picture_control_set_ptr,
    EbPictureBufferDesc           *input_picture_ptr,
    EbPictureBufferDesc           *input_padded_picture_ptr,
    uint32_t                           sb_total_count)
{
    uint32_t sb_index;
    uint32_t sb_origin_x;        // to avoid using child PCS
    uint32_t sb_origin_y;
    uint32_t inputLumaOriginIndex;
    uint32_t inputCbOriginIndex;
    uint32_t inputCrOriginIndex;
    uint64_t picTotVariance;

    // Variance
    picTotVariance = 0;

    for (sb_index = 0; sb_index < picture_control_set_ptr->sb_total_count; ++sb_index) {
        SbParams   *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];

        sb_origin_x = sb_params->origin_x;
        sb_origin_y = sb_params->origin_y;
        inputLumaOriginIndex = (input_padded_picture_ptr->origin_y + sb_origin_y) * input_padded_picture_ptr->stride_y +
            input_padded_picture_ptr->origin_x + sb_origin_x;

        inputCbOriginIndex = ((input_picture_ptr->origin_y + sb_origin_y) >> 1) * input_picture_ptr->stride_cb + ((input_picture_ptr->origin_x + sb_origin_x) >> 1);
        inputCrOriginIndex = ((input_picture_ptr->origin_y + sb_origin_y) >> 1) * input_picture_ptr->stride_cr + ((input_picture_ptr->origin_x + sb_origin_x) >> 1);

        ComputeBlockMeanComputeVariance(
            sequence_control_set_ptr,
            picture_control_set_ptr,
            input_padded_picture_ptr,
            sb_index,
            inputLumaOriginIndex);

        if (sb_params->is_complete_sb) {
            ComputeChromaBlockMean(
                sequence_control_set_ptr,
                picture_control_set_ptr,
                input_picture_ptr,
                sb_index,
                inputCbOriginIndex,
                inputCrOriginIndex);
        }
        else {
            ZeroOutChromaBlockMean(
                picture_control_set_ptr,
                sb_index);
        }

        picTotVariance += (picture_control_set_ptr->variance[sb_index][RASTER_SCAN_CU_INDEX_64x64]);
    }

    picture_control_set_ptr->pic_avg_variance = (uint16_t)(picTotVariance / sb_total_count);

    // Calculate the variance of variance to determine Homogeneous regions. Note: Variance calculation should be on.
    DetermineHomogeneousRegionInPicture(
        sequence_control_set_ptr,
        picture_control_set_ptr);

    EdgeDetectionMeanLumaChroma16x16(
        sequence_control_set_ptr,
        picture_control_set_ptr,
        sequence_control_set_ptr->sb_total_count);

    EdgeDetection(
        sequence_control_set_ptr,
        picture_control_set_ptr);

    return;
}

void CalculateInputAverageIntensity(
    SequenceControlSet            *sequence_control_set_ptr,
    PictureParentControlSet       *picture_control_set_ptr,
    EbPictureBufferDesc           *input_picture_ptr,
    uint64_t                           sumAverageIntensityTotalRegionsLuma,
    uint64_t                           sumAverageIntensityTotalRegionsCb,
    uint64_t                           sumAverageIntensityTotalRegionsCr)
{
    if (sequence_control_set_ptr->scd_mode == SCD_MODE_0) {
        uint16_t blockIndexInWidth;
        uint16_t blockIndexInHeight;
        uint64_t mean = 0;

        const uint16_t stride_y = input_picture_ptr->stride_y;
        // Loop over 8x8 blocks and calculates the mean value
        if (sequence_control_set_ptr->block_mean_calc_prec == BLOCK_MEAN_PREC_FULL) {
            for (blockIndexInHeight = 0; blockIndexInHeight < input_picture_ptr->height >> 3; ++blockIndexInHeight) {
                for (blockIndexInWidth = 0; blockIndexInWidth < input_picture_ptr->width >> 3; ++blockIndexInWidth)
                    mean += compute_mean_8x8(&(input_picture_ptr->buffer_y[(blockIndexInWidth << 3) + (blockIndexInHeight << 3) * input_picture_ptr->stride_y]), input_picture_ptr->stride_y, 8, 8);
            }
        }
        else {
            for (blockIndexInHeight = 0; blockIndexInHeight < input_picture_ptr->height >> 3; ++blockIndexInHeight) {
                for (blockIndexInWidth = 0; blockIndexInWidth < input_picture_ptr->width >> 3; ++blockIndexInWidth)
                    mean += compute_sub_mean8x8_sse2_intrin(&(input_picture_ptr->buffer_y[(blockIndexInWidth << 3) + (blockIndexInHeight << 3) * stride_y]), stride_y);
            }
        }
        mean = ((mean + ((input_picture_ptr->height* input_picture_ptr->width) >> 7)) / ((input_picture_ptr->height* input_picture_ptr->width) >> 6));
        mean = (mean + (1 << (MEAN_PRECISION - 1))) >> MEAN_PRECISION;
        picture_control_set_ptr->average_intensity[0] = (uint8_t)mean;
    }

    else {
        picture_control_set_ptr->average_intensity[0] = (uint8_t)((sumAverageIntensityTotalRegionsLuma + ((input_picture_ptr->width*input_picture_ptr->height) >> 1)) / (input_picture_ptr->width*input_picture_ptr->height));
        picture_control_set_ptr->average_intensity[1] = (uint8_t)((sumAverageIntensityTotalRegionsCb + ((input_picture_ptr->width*input_picture_ptr->height) >> 3)) / ((input_picture_ptr->width*input_picture_ptr->height) >> 2));
        picture_control_set_ptr->average_intensity[2] = (uint8_t)((sumAverageIntensityTotalRegionsCr + ((input_picture_ptr->width*input_picture_ptr->height) >> 3)) / ((input_picture_ptr->width*input_picture_ptr->height) >> 2));
    }

    return;
}

/************************************************
 * Gathering statistics per picture
 ** Calculating the pixel intensity histogram bins per picture needed for SCD
 ** Computing Picture Variance
 ************************************************/
void GatheringPictureStatistics(
    SequenceControlSet            *sequence_control_set_ptr,
    PictureParentControlSet       *picture_control_set_ptr,
    EbPictureBufferDesc           *input_picture_ptr,
    EbPictureBufferDesc           *input_padded_picture_ptr,
    EbPictureBufferDesc            *sixteenth_decimated_picture_ptr,
    uint32_t                           sb_total_count)
{
    uint64_t                          sumAverageIntensityTotalRegionsLuma = 0;
    uint64_t                          sumAverageIntensityTotalRegionsCb = 0;
    uint64_t                          sumAverageIntensityTotalRegionsCr = 0;

    // Histogram bins
        // Use 1/16 Luma for Histogram generation
        // 1/16 input ready
    SubSampleLumaGeneratePixelIntensityHistogramBins(
        sequence_control_set_ptr,
        picture_control_set_ptr,
        sixteenth_decimated_picture_ptr,
        &sumAverageIntensityTotalRegionsLuma);

    // Use 1/4 Chroma for Histogram generation
    // 1/4 input not ready => perform operation on the fly
    SubSampleChromaGeneratePixelIntensityHistogramBins(
        sequence_control_set_ptr,
        picture_control_set_ptr,
        input_picture_ptr,
        &sumAverageIntensityTotalRegionsCb,
        &sumAverageIntensityTotalRegionsCr);
    //
    // Calculate the LUMA average intensity
    CalculateInputAverageIntensity(
        sequence_control_set_ptr,
        picture_control_set_ptr,
        input_picture_ptr,
        sumAverageIntensityTotalRegionsLuma,
        sumAverageIntensityTotalRegionsCb,
        sumAverageIntensityTotalRegionsCr);

    ComputePictureSpatialStatistics(
        sequence_control_set_ptr,
        picture_control_set_ptr,
        input_picture_ptr,
        input_padded_picture_ptr,
        sb_total_count);

    return;
}
/************************************************
 * Pad Picture at the right and bottom sides
 ** To match a multiple of min CU size in width and height
 ************************************************/
void PadPictureToMultipleOfMinCuSizeDimensions(
    SequenceControlSet            *sequence_control_set_ptr,
    EbPictureBufferDesc           *input_picture_ptr)
{
    EbBool                          is16BitInput = (EbBool)(sequence_control_set_ptr->static_config.encoder_bit_depth > EB_8BIT);

    uint32_t color_format = input_picture_ptr->color_format;
    const uint16_t subsampling_x = (color_format == EB_YUV444 ? 1 : 2) - 1;
    const uint16_t subsampling_y = (color_format >= EB_YUV422 ? 1 : 2) - 1;

    // Input Picture Padding
    pad_input_picture(
        &input_picture_ptr->buffer_y[input_picture_ptr->origin_x + (input_picture_ptr->origin_y * input_picture_ptr->stride_y)],
        input_picture_ptr->stride_y,
        (input_picture_ptr->width - sequence_control_set_ptr->pad_right),
        (input_picture_ptr->height - sequence_control_set_ptr->pad_bottom),
        sequence_control_set_ptr->pad_right,
        sequence_control_set_ptr->pad_bottom);

    pad_input_picture(
        &input_picture_ptr->buffer_cb[(input_picture_ptr->origin_x >> subsampling_x) + ((input_picture_ptr->origin_y >> subsampling_y) * input_picture_ptr->stride_cb)],
        input_picture_ptr->stride_cb,
        (input_picture_ptr->width - sequence_control_set_ptr->pad_right) >> subsampling_x,
        (input_picture_ptr->height - sequence_control_set_ptr->pad_bottom) >> subsampling_y,
        sequence_control_set_ptr->pad_right >> subsampling_x,
        sequence_control_set_ptr->pad_bottom >> subsampling_y);

    pad_input_picture(
        &input_picture_ptr->buffer_cr[(input_picture_ptr->origin_x >> subsampling_x) + ((input_picture_ptr->origin_y >> subsampling_y) * input_picture_ptr->stride_cb)],
        input_picture_ptr->stride_cr,
        (input_picture_ptr->width - sequence_control_set_ptr->pad_right) >> subsampling_x,
        (input_picture_ptr->height - sequence_control_set_ptr->pad_bottom) >> subsampling_y,
        sequence_control_set_ptr->pad_right >> subsampling_x,
        sequence_control_set_ptr->pad_bottom >> subsampling_y);

    if (is16BitInput)
    {
        pad_input_picture(
            &input_picture_ptr->buffer_bit_inc_y[input_picture_ptr->origin_x + (input_picture_ptr->origin_y * input_picture_ptr->stride_bit_inc_y)],
            input_picture_ptr->stride_bit_inc_y,
            (input_picture_ptr->width - sequence_control_set_ptr->pad_right),
            (input_picture_ptr->height - sequence_control_set_ptr->pad_bottom),
            sequence_control_set_ptr->pad_right,
            sequence_control_set_ptr->pad_bottom);

        pad_input_picture(
            &input_picture_ptr->buffer_bit_inc_cb[(input_picture_ptr->origin_x >> subsampling_x) + ((input_picture_ptr->origin_y >> subsampling_y) * input_picture_ptr->stride_bit_inc_cb)],
            input_picture_ptr->stride_bit_inc_cb,
            (input_picture_ptr->width - sequence_control_set_ptr->pad_right) >> subsampling_x,
            (input_picture_ptr->height - sequence_control_set_ptr->pad_bottom) >> subsampling_y,
            sequence_control_set_ptr->pad_right >> subsampling_x,
            sequence_control_set_ptr->pad_bottom >> subsampling_y);

        pad_input_picture(
            &input_picture_ptr->buffer_bit_inc_cr[(input_picture_ptr->origin_x >> subsampling_x) + ((input_picture_ptr->origin_y >> subsampling_y) * input_picture_ptr->stride_bit_inc_cb)],
            input_picture_ptr->stride_bit_inc_cr,
            (input_picture_ptr->width - sequence_control_set_ptr->pad_right) >> subsampling_x,
            (input_picture_ptr->height - sequence_control_set_ptr->pad_bottom) >> subsampling_y,
            sequence_control_set_ptr->pad_right >> subsampling_x,
            sequence_control_set_ptr->pad_bottom >> subsampling_y);
    }

    return;
}

/************************************************
 * Pad Picture at the right and bottom sides
 ** To complete border SB smaller than SB size
 ************************************************/
void PadPictureToMultipleOfLcuDimensions(
    EbPictureBufferDesc           *input_padded_picture_ptr
)
{
    // Generate Padding
    generate_padding(
        &input_padded_picture_ptr->buffer_y[0],
        input_padded_picture_ptr->stride_y,
        input_padded_picture_ptr->width,
        input_padded_picture_ptr->height,
        input_padded_picture_ptr->origin_x,
        input_padded_picture_ptr->origin_y);

    return;
}

/************************************************
* 1/4 & 1/16 input picture decimation
************************************************/
void DownsampleDecimationInputPicture(
    PictureParentControlSet       *picture_control_set_ptr,
    EbPictureBufferDesc           *input_padded_picture_ptr,
    EbPictureBufferDesc           *quarter_decimated_picture_ptr,
    EbPictureBufferDesc           *sixteenth_decimated_picture_ptr) {
    // Decimate input picture for HME L0 and L1
    if (picture_control_set_ptr->enable_hme_flag || picture_control_set_ptr->tf_enable_hme_flag) {
        if (picture_control_set_ptr->enable_hme_level1_flag || picture_control_set_ptr->tf_enable_hme_level1_flag) {
            decimation_2d(
                &input_padded_picture_ptr->buffer_y[input_padded_picture_ptr->origin_x + input_padded_picture_ptr->origin_y * input_padded_picture_ptr->stride_y],
                input_padded_picture_ptr->stride_y,
                input_padded_picture_ptr->width,
                input_padded_picture_ptr->height,
                &quarter_decimated_picture_ptr->buffer_y[quarter_decimated_picture_ptr->origin_x + quarter_decimated_picture_ptr->origin_x*quarter_decimated_picture_ptr->stride_y],
                quarter_decimated_picture_ptr->stride_y,
                2);
            generate_padding(
                &quarter_decimated_picture_ptr->buffer_y[0],
                quarter_decimated_picture_ptr->stride_y,
                quarter_decimated_picture_ptr->width,
                quarter_decimated_picture_ptr->height,
                quarter_decimated_picture_ptr->origin_x,
                quarter_decimated_picture_ptr->origin_y);
        }
    }

    // Always perform 1/16th decimation as
    // Sixteenth Input Picture Decimation
    decimation_2d(
        &input_padded_picture_ptr->buffer_y[input_padded_picture_ptr->origin_x + input_padded_picture_ptr->origin_y * input_padded_picture_ptr->stride_y],
        input_padded_picture_ptr->stride_y,
        input_padded_picture_ptr->width,
        input_padded_picture_ptr->height,
        &sixteenth_decimated_picture_ptr->buffer_y[sixteenth_decimated_picture_ptr->origin_x + sixteenth_decimated_picture_ptr->origin_x*sixteenth_decimated_picture_ptr->stride_y],
        sixteenth_decimated_picture_ptr->stride_y,
        4);

    generate_padding(
        &sixteenth_decimated_picture_ptr->buffer_y[0],
        sixteenth_decimated_picture_ptr->stride_y,
        sixteenth_decimated_picture_ptr->width,
        sixteenth_decimated_picture_ptr->height,
        sixteenth_decimated_picture_ptr->origin_x,
        sixteenth_decimated_picture_ptr->origin_y);

}
#if PAL_SUP
int av1_count_colors_highbd(uint16_t *src, int stride, int rows, int cols,
    int bit_depth, int *val_count) {
    assert(bit_depth <= 12);
    const int max_pix_val = 1 << bit_depth;
   // const uint16_t *src = CONVERT_TO_SHORTPTR(src8);
    memset(val_count, 0, max_pix_val * sizeof(val_count[0]));
    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            const int this_val = src[r * stride + c];
            assert(this_val < max_pix_val);
            if (this_val >= max_pix_val) return 0;
            ++val_count[this_val];
        }
    }
    int n = 0;
    for (int i = 0; i < max_pix_val; ++i) {
        if (val_count[i]) ++n;
    }
    return n;
}
#endif
int eb_av1_count_colors(const uint8_t *src, int stride, int rows, int cols,
    int *val_count) {
    const int max_pix_val = 1 << 8;
    memset(val_count, 0, max_pix_val * sizeof(val_count[0]));
    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            const int this_val = src[r * stride + c];
            assert(this_val < max_pix_val);
            ++val_count[this_val];
        }
    }
    int n = 0;
    for (int i = 0; i < max_pix_val; ++i)
        if (val_count[i]) ++n;
    return n;
}
extern aom_variance_fn_ptr_t mefn_ptr[BlockSizeS_ALL];

// This is used as a reference when computing the source variance for the
//  purposes of activity masking.
// Eventually this should be replaced by custom no-reference routines,
//  which will be faster.
const uint8_t eb_AV1_VAR_OFFS[MAX_SB_SIZE] = {
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128
};

unsigned int eb_av1_get_sby_perpixel_variance(const aom_variance_fn_ptr_t *fn_ptr, //const AV1_COMP *cpi,
                                           const uint8_t *src,int stride,//const struct buf_2d *ref,
                                           BlockSize bs) {
  unsigned int sse;
  const unsigned int var =
      //cpi->fn_ptr[bs].vf(ref->buf, ref->stride, eb_AV1_VAR_OFFS, 0, &sse);
     fn_ptr->vf(src,  stride, eb_AV1_VAR_OFFS, 0, &sse);
  return ROUND_POWER_OF_TWO(var, num_pels_log2_lookup[bs]);
}

// Estimate if the source frame is screen content, based on the portion of
// blocks that have no more than 4 (experimentally selected) luma colors.
static void is_screen_content(
    PictureParentControlSet     *picture_control_set_ptr,
    const uint8_t               *src,
    int                          use_hbd,
    int                          stride,
    int                         width,
    int                         height) {
    assert(src != NULL);
    const int blk_w = 16;
    const int blk_h = 16;
    // These threshold values are selected experimentally.
    const int color_thresh = 4;
    const unsigned int var_thresh = 0;
    // Counts of blocks with no more than color_thresh colors.
    int counts_1 = 0;
    // Counts of blocks with no more than color_thresh colors and variance larger
    // than var_thresh.
    int counts_2 = 0;

    for (int r = 0; r + blk_h <= height; r += blk_h) {
        for (int c = 0; c + blk_w <= width; c += blk_w) {
            int count_buf[1 << 12];  // Maximum (1 << 12) color levels.
            const int n_colors =
                use_hbd ? 0 /*av1_count_colors_highbd(src + r * stride + c, stride, blk_w,
                    blk_h, bd, count_buf)*/
                : eb_av1_count_colors(src + r * stride + c, stride, blk_w, blk_h,
                    count_buf);
            if (n_colors > 1 && n_colors <= color_thresh) {
                ++counts_1;
                //struct buf_2d buf;
                //buf.stride = stride;
                //buf.buf = (uint8_t *)src;
                const aom_variance_fn_ptr_t *fn_ptr = &mefn_ptr[BLOCK_16X16];

                const unsigned int var = eb_av1_get_sby_perpixel_variance(fn_ptr, src + r * stride + c,stride, BLOCK_16X16);
                               /* use_hbd
                ? av1_high_get_sby_perpixel_variance(cpi, &buf, BLOCK_16X16, bd)
                : */
                if (var > var_thresh) ++counts_2;
            }
        }
    }

    picture_control_set_ptr->sc_content_detected =
        (counts_1 * blk_h * blk_w * 10 > width * height) &&
        ( counts_2 * blk_h * blk_w * 15 > width * height) ;
}


/************************************************
 * 1/4 & 1/16 input picture downsampling (filtering)
 ************************************************/
void DownsampleFilteringInputPicture(
    PictureParentControlSet       *picture_control_set_ptr,
    EbPictureBufferDesc           *input_padded_picture_ptr,
    EbPictureBufferDesc           *quarter_picture_ptr,
    EbPictureBufferDesc           *sixteenth_picture_ptr) {

    // Downsample input picture for HME L0 and L1
    if (picture_control_set_ptr->enable_hme_flag || picture_control_set_ptr->tf_enable_hme_flag) {
        if (picture_control_set_ptr->enable_hme_level1_flag || picture_control_set_ptr->tf_enable_hme_level1_flag) {

            downsample_2d(
                &input_padded_picture_ptr->buffer_y[input_padded_picture_ptr->origin_x + input_padded_picture_ptr->origin_y * input_padded_picture_ptr->stride_y],
                input_padded_picture_ptr->stride_y,
                input_padded_picture_ptr->width,
                input_padded_picture_ptr->height,
                &quarter_picture_ptr->buffer_y[quarter_picture_ptr->origin_x + quarter_picture_ptr->origin_x * quarter_picture_ptr->stride_y],
                quarter_picture_ptr->stride_y,
                2);
            generate_padding(
                &quarter_picture_ptr->buffer_y[0],
                quarter_picture_ptr->stride_y,
                quarter_picture_ptr->width,
                quarter_picture_ptr->height,
                quarter_picture_ptr->origin_x,
                quarter_picture_ptr->origin_y);

        }

        if (picture_control_set_ptr->enable_hme_level0_flag || picture_control_set_ptr->tf_enable_hme_level0_flag) {
            // Sixteenth Input Picture Downsampling
            if (picture_control_set_ptr->enable_hme_level1_flag || picture_control_set_ptr->tf_enable_hme_level1_flag)
                downsample_2d(
                    &quarter_picture_ptr->buffer_y[quarter_picture_ptr->origin_x + quarter_picture_ptr->origin_y * quarter_picture_ptr->stride_y],
                    quarter_picture_ptr->stride_y,
                    quarter_picture_ptr->width,
                    quarter_picture_ptr->height,
                    &sixteenth_picture_ptr->buffer_y[sixteenth_picture_ptr->origin_x + sixteenth_picture_ptr->origin_x*sixteenth_picture_ptr->stride_y],
                    sixteenth_picture_ptr->stride_y,
                    2);
            else
                downsample_2d(
                    &input_padded_picture_ptr->buffer_y[input_padded_picture_ptr->origin_x + input_padded_picture_ptr->origin_y * input_padded_picture_ptr->stride_y],
                    input_padded_picture_ptr->stride_y,
                    input_padded_picture_ptr->width,
                    input_padded_picture_ptr->height,
                    &sixteenth_picture_ptr->buffer_y[sixteenth_picture_ptr->origin_x + sixteenth_picture_ptr->origin_x*sixteenth_picture_ptr->stride_y],
                    sixteenth_picture_ptr->stride_y,
                    4);

            generate_padding(
                &sixteenth_picture_ptr->buffer_y[0],
                sixteenth_picture_ptr->stride_y,
                sixteenth_picture_ptr->width,
                sixteenth_picture_ptr->height,
                sixteenth_picture_ptr->origin_x,
                sixteenth_picture_ptr->origin_y);

        }
    }
}

/************************************************
 * Picture Analysis Kernel
 * The Picture Analysis Process pads & decimates the input pictures.
 * The Picture Analysis also includes creating an n-bin Histogram,
 * gathering picture 1st and 2nd moment statistics for each 8x8 block,
 * which are used to compute variance.
 * The Picture Analysis process is multithreaded, so pictures can be
 * processed out of order as long as all inputs are available.
 ************************************************/
void* picture_analysis_kernel(void *input_ptr)
{
    PictureAnalysisContext        *context_ptr = (PictureAnalysisContext*)input_ptr;
    PictureParentControlSet       *picture_control_set_ptr;
    SequenceControlSet            *sequence_control_set_ptr;

    EbObjectWrapper               *inputResultsWrapperPtr;
    ResourceCoordinationResults   *inputResultsPtr;
    EbObjectWrapper               *outputResultsWrapperPtr;
    PictureAnalysisResults        *outputResultsPtr;
    EbPaReferenceObject           *paReferenceObject;

    EbPictureBufferDesc           *input_padded_picture_ptr;
    EbPictureBufferDesc           *input_picture_ptr;

    // Variance
    uint32_t                        picture_width_in_sb;
    uint32_t                        pictureHeighInLcu;
    uint32_t                        sb_total_count;

    for (;;) {
        // Get Input Full Object
        eb_get_full_object(
            context_ptr->resource_coordination_results_input_fifo_ptr,
            &inputResultsWrapperPtr);

        inputResultsPtr = (ResourceCoordinationResults*)inputResultsWrapperPtr->object_ptr;
        picture_control_set_ptr = (PictureParentControlSet*)inputResultsPtr->picture_control_set_wrapper_ptr->object_ptr;

        // There is no need to do processing for overlay picture. Overlay and AltRef share the same results.
        if (!picture_control_set_ptr->is_overlay)
        {
            sequence_control_set_ptr = (SequenceControlSet*)picture_control_set_ptr->sequence_control_set_wrapper_ptr->object_ptr;
            input_picture_ptr = picture_control_set_ptr->enhanced_picture_ptr;

            paReferenceObject = (EbPaReferenceObject*)picture_control_set_ptr->pa_reference_picture_wrapper_ptr->object_ptr;
            input_padded_picture_ptr = (EbPictureBufferDesc*)paReferenceObject->input_padded_picture_ptr;
            // Variance
            picture_width_in_sb = (sequence_control_set_ptr->seq_header.max_frame_width + sequence_control_set_ptr->sb_sz - 1) / sequence_control_set_ptr->sb_sz;
            pictureHeighInLcu = (sequence_control_set_ptr->seq_header.max_frame_height + sequence_control_set_ptr->sb_sz - 1) / sequence_control_set_ptr->sb_sz;
            sb_total_count = picture_width_in_sb * pictureHeighInLcu;

            // Set picture parameters to account for subpicture, picture scantype, and set regions by resolutions
            SetPictureParametersForStatisticsGathering(
                sequence_control_set_ptr);

            // Pad pictures to multiple min cu size
            PadPictureToMultipleOfMinCuSizeDimensions(
                sequence_control_set_ptr,
                input_picture_ptr);

            // Pre processing operations performed on the input picture
            PicturePreProcessingOperations(
                picture_control_set_ptr,
                sequence_control_set_ptr,
                sb_total_count);
            if (input_picture_ptr->color_format >= EB_YUV422) {
                // Jing: Do the conversion of 422/444=>420 here since it's multi-threaded kernel
                //       Reuse the Y, only add cb/cr in the newly created buffer desc
                //       NOTE: since denoise may change the src, so this part is after PicturePreProcessingOperations()
                picture_control_set_ptr->chroma_downsampled_picture_ptr->buffer_y = input_picture_ptr->buffer_y;
                DownSampleChroma(input_picture_ptr, picture_control_set_ptr->chroma_downsampled_picture_ptr);
            }
            else
                picture_control_set_ptr->chroma_downsampled_picture_ptr = input_picture_ptr;
            // Pad input picture to complete border LCUs
            PadPictureToMultipleOfLcuDimensions(
                input_padded_picture_ptr);
            // 1/4 & 1/16 input picture decimation
            DownsampleDecimationInputPicture(
                picture_control_set_ptr,
                input_padded_picture_ptr,
                (EbPictureBufferDesc*)paReferenceObject->quarter_decimated_picture_ptr,
                (EbPictureBufferDesc*)paReferenceObject->sixteenth_decimated_picture_ptr);

            // 1/4 & 1/16 input picture downsampling through filtering
            if (sequence_control_set_ptr->down_sampling_method_me_search == ME_FILTERED_DOWNSAMPLED) {
                DownsampleFilteringInputPicture(
                    picture_control_set_ptr,
                    input_padded_picture_ptr,
                    (EbPictureBufferDesc*)paReferenceObject->quarter_filtered_picture_ptr,
                    (EbPictureBufferDesc*)paReferenceObject->sixteenth_filtered_picture_ptr);
            }
           // Gathering statistics of input picture, including Variance Calculation, Histogram Bins
            GatheringPictureStatistics(
                sequence_control_set_ptr,
                picture_control_set_ptr,
                picture_control_set_ptr->chroma_downsampled_picture_ptr, //420 input_picture_ptr
                input_padded_picture_ptr,
                (EbPictureBufferDesc*)paReferenceObject->sixteenth_decimated_picture_ptr, // Hsan: always use decimated until studying the trade offs
                sb_total_count);

            if (sequence_control_set_ptr->static_config.screen_content_mode == 2){ // auto detect
                is_screen_content(
                    picture_control_set_ptr,
                    input_picture_ptr->buffer_y + input_picture_ptr->origin_x + input_picture_ptr->origin_y*input_picture_ptr->stride_y,
                    0,
                    input_picture_ptr->stride_y,
                    sequence_control_set_ptr->seq_header.max_frame_width, sequence_control_set_ptr->seq_header.max_frame_height);
            }
            else // off / on
                picture_control_set_ptr->sc_content_detected = sequence_control_set_ptr->static_config.screen_content_mode;

            // Hold the 64x64 variance and mean in the reference frame
            uint32_t sb_index;
            for (sb_index = 0; sb_index < picture_control_set_ptr->sb_total_count; ++sb_index) {
                paReferenceObject->variance[sb_index] = picture_control_set_ptr->variance[sb_index][ME_TIER_ZERO_PU_64x64];
                paReferenceObject->y_mean[sb_index] = picture_control_set_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_64x64];
            }
        }
        // Get Empty Results Object
        eb_get_empty_object(
            context_ptr->picture_analysis_results_output_fifo_ptr,
            &outputResultsWrapperPtr);

        outputResultsPtr = (PictureAnalysisResults*)outputResultsWrapperPtr->object_ptr;
        outputResultsPtr->picture_control_set_wrapper_ptr = inputResultsPtr->picture_control_set_wrapper_ptr;

        // Release the Input Results
        eb_release_object(inputResultsWrapperPtr);

        // Post the Full Results Object
        eb_post_full_object(outputResultsWrapperPtr);
    }
    return EB_NULL;
}
