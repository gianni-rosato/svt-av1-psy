/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <stdlib.h>
#include <string.h>

#include "aom_dsp_rtcd.h"
#include "EbDefinitions.h"
#include "EbEncHandle.h"
#include "EbSystemResourceManager.h"
#include "EbPictureControlSet.h"
#include "EbSequenceControlSet.h"
#include "EbPictureBufferDesc.h"

#include "EbResourceCoordinationResults.h"
#include "EbPictureAnalysisProcess.h"
#include "EbPictureAnalysisResults.h"
#include "EbMcp.h"
#include "EbReferenceObject.h"
//#include "EbComputeMean_SSE2.h"
#include "EbUtility.h"
#include "EbMotionEstimationContext.h"
#include "EbPictureOperators.h"
#include "EbResize.h"

#define VARIANCE_PRECISION 16
#define SB_LOW_VAR_TH 5
#define PIC_LOW_VAR_PERCENTAGE_TH 60
#define FLAT_MAX_VAR 50
#define FLAT_MAX_VAR_DECIM (50 - 00)
#define NOISE_MIN_LEVEL 70000 //120000
#define NOISE_MIN_LEVEL_DECIM (70000 + 000000) //(120000+000000)
#define NOISE_MIN_LEVEL_M6_M7 120000
#define NOISE_MIN_LEVEL_DECIM_M6_M7 (120000 + 000000)
#define DENOISER_QP_TH 29
#define DENOISER_BITRATE_TH 14000000
#define SAMPLE_THRESHOLD_PRECENT_BORDER_LINE 15
#define SAMPLE_THRESHOLD_PRECENT_TWO_BORDER_LINES 10

/**************************************
 * Context
 **************************************/
typedef struct PictureAnalysisContext {
    EB_ALIGN(64) uint8_t local_cache[64];
    EbFifo *             resource_coordination_results_input_fifo_ptr;
    EbFifo *             picture_analysis_results_output_fifo_ptr;
    EbPictureBufferDesc *denoised_picture_ptr;
    EbPictureBufferDesc *noise_picture_ptr;
    double               pic_noise_variance_float;
} PictureAnalysisContext;

static void picture_analysis_context_dctor(EbPtr p) {
    EbThreadContext *       thread_context_ptr = (EbThreadContext *)p;
    PictureAnalysisContext *obj                = (PictureAnalysisContext *)thread_context_ptr->priv;
    EB_DELETE(obj->noise_picture_ptr);
    EB_DELETE(obj->denoised_picture_ptr);
    EB_FREE_ARRAY(obj);
}
/************************************************
* Picture Analysis Context Constructor
************************************************/
EbErrorType picture_analysis_context_ctor(EbThreadContext *  thread_context_ptr,
                                          const EbEncHandle *enc_handle_ptr, int index) {
    EbBool denoise_flag = EB_TRUE;

    PictureAnalysisContext *context_ptr;
    EB_CALLOC_ARRAY(context_ptr, 1);
    thread_context_ptr->priv  = context_ptr;
    thread_context_ptr->dctor = picture_analysis_context_dctor;

    context_ptr->resource_coordination_results_input_fifo_ptr =
        eb_system_resource_get_consumer_fifo(
            enc_handle_ptr->resource_coordination_results_resource_ptr, index);
    context_ptr->picture_analysis_results_output_fifo_ptr = eb_system_resource_get_producer_fifo(
        enc_handle_ptr->picture_analysis_results_resource_ptr, index);

    if (denoise_flag == EB_TRUE) {
        EbPictureBufferDescInitData desc;
        const SequenceControlSet *  scs_ptr = enc_handle_ptr->scs_instance_array[0]->scs_ptr;

        memset(&desc, 0, sizeof(desc));
        desc.color_format       = scs_ptr->static_config.encoder_color_format;
        desc.max_width          = scs_ptr->max_input_luma_width;
        desc.max_height         = scs_ptr->max_input_luma_height;
        desc.bit_depth          = EB_8BIT;
        desc.buffer_enable_mask = PICTURE_BUFFER_DESC_Y_FLAG;
        //denoised
        // If 420/422, re-use luma for chroma
        // If 444, re-use luma for Cr
        if (desc.color_format != EB_YUV444) {
            desc.buffer_enable_mask = PICTURE_BUFFER_DESC_Y_FLAG;
        } else
            desc.buffer_enable_mask = PICTURE_BUFFER_DESC_Y_FLAG | PICTURE_BUFFER_DESC_Cb_FLAG;
        EB_NEW(context_ptr->denoised_picture_ptr, eb_picture_buffer_desc_ctor, (EbPtr)&desc);

        if (desc.color_format != EB_YUV444) {
            context_ptr->denoised_picture_ptr->buffer_cb =
                context_ptr->denoised_picture_ptr->buffer_y;
            context_ptr->denoised_picture_ptr->buffer_cr =
                context_ptr->denoised_picture_ptr->buffer_y +
                context_ptr->denoised_picture_ptr->chroma_size;
        } else
            context_ptr->denoised_picture_ptr->buffer_cr =
                context_ptr->denoised_picture_ptr->buffer_y;
        // noise
        desc.max_height         = BLOCK_SIZE_64;
        desc.buffer_enable_mask = PICTURE_BUFFER_DESC_Y_FLAG;

        EB_NEW(context_ptr->noise_picture_ptr, eb_picture_buffer_desc_ctor, (EbPtr)&desc);
    }
    return EB_ErrorNone;
}
void down_sample_chroma(EbPictureBufferDesc *input_picture_ptr,
                        EbPictureBufferDesc *outputPicturePtr) {
    uint32_t       input_color_format  = input_picture_ptr->color_format;
    const uint16_t input_subsampling_x = (input_color_format == EB_YUV444 ? 1 : 2) - 1;
    const uint16_t input_subsampling_y = (input_color_format >= EB_YUV422 ? 1 : 2) - 1;

    uint32_t       output_color_format  = outputPicturePtr->color_format;
    const uint16_t output_subsampling_x = (output_color_format == EB_YUV444 ? 1 : 2) - 1;
    const uint16_t output_subsampling_y = (output_color_format >= EB_YUV422 ? 1 : 2) - 1;

    uint32_t stride_in, stride_out;
    uint32_t input_origin_index, output_origin_index;

    uint8_t *ptr_in;
    uint8_t *ptr_out;

    uint32_t ii, jj;

    //Cb
    {
        stride_in = input_picture_ptr->stride_cb;
        input_origin_index =
            (input_picture_ptr->origin_x >> input_subsampling_x) +
            (input_picture_ptr->origin_y >> input_subsampling_y) * input_picture_ptr->stride_cb;
        ptr_in = &(input_picture_ptr->buffer_cb[input_origin_index]);

        stride_out = outputPicturePtr->stride_cb;
        output_origin_index =
            (outputPicturePtr->origin_x >> output_subsampling_x) +
            (outputPicturePtr->origin_y >> output_subsampling_y) * outputPicturePtr->stride_cb;
        ptr_out = &(outputPicturePtr->buffer_cb[output_origin_index]);

        for (jj = 0; jj < (uint32_t)(outputPicturePtr->height >> output_subsampling_y); jj++) {
            for (ii = 0; ii < (uint32_t)(outputPicturePtr->width >> output_subsampling_x); ii++) {
                ptr_out[ii + jj * stride_out] =
                    ptr_in[(ii << (1 - input_subsampling_x)) +
                           (jj << (1 - input_subsampling_y)) * stride_in];
            }
        }
    }

    //Cr
    {
        stride_in = input_picture_ptr->stride_cr;
        input_origin_index =
            (input_picture_ptr->origin_x >> input_subsampling_x) +
            (input_picture_ptr->origin_y >> input_subsampling_y) * input_picture_ptr->stride_cr;
        ptr_in = &(input_picture_ptr->buffer_cr[input_origin_index]);

        stride_out = outputPicturePtr->stride_cr;
        output_origin_index =
            (outputPicturePtr->origin_x >> output_subsampling_x) +
            (outputPicturePtr->origin_y >> output_subsampling_y) * outputPicturePtr->stride_cr;
        ptr_out = &(outputPicturePtr->buffer_cr[output_origin_index]);

        for (jj = 0; jj < (uint32_t)(outputPicturePtr->height >> output_subsampling_y); jj++) {
            for (ii = 0; ii < (uint32_t)(outputPicturePtr->width >> output_subsampling_x); ii++) {
                ptr_out[ii + jj * stride_out] =
                    ptr_in[(ii << (1 - input_subsampling_x)) +
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
void decimation_2d(uint8_t *input_samples, // input parameter, input samples Ptr
                   uint32_t input_stride, // input parameter, input stride
                   uint32_t input_area_width, // input parameter, input area width
                   uint32_t input_area_height, // input parameter, input area height
                   uint8_t *decim_samples, // output parameter, decimated samples Ptr
                   uint32_t decim_stride, // input parameter, output stride
                   uint32_t decim_step) // input parameter, decimation amount in pixels
{
    uint32_t horizontal_index;
    uint32_t vertical_index;
    uint32_t input_stripe_stride = input_stride * decim_step;

    for (vertical_index = 0; vertical_index < input_area_height; vertical_index += decim_step) {
        for (horizontal_index = 0; horizontal_index < input_area_width;
             horizontal_index += decim_step)
            decim_samples[(horizontal_index >> (decim_step >> 1))] =
                input_samples[horizontal_index];

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
void downsample_2d(uint8_t *input_samples, // input parameter, input samples Ptr
                   uint32_t input_stride, // input parameter, input stride
                   uint32_t input_area_width, // input parameter, input area width
                   uint32_t input_area_height, // input parameter, input area height
                   uint8_t *decim_samples, // output parameter, decimated samples Ptr
                   uint32_t decim_stride, // input parameter, output stride
                   uint32_t decim_step) // input parameter, decimation amount in pixels
{
    uint32_t       horizontal_index;
    uint32_t       vertical_index;
    uint32_t       input_stripe_stride = input_stride * decim_step;
    uint32_t       decim_horizontal_index;
    const uint32_t half_decim_step = decim_step >> 1;

    for (input_samples += half_decim_step * input_stride, vertical_index = half_decim_step;
         vertical_index < input_area_height;
         vertical_index += decim_step) {
        uint8_t *prev_input_line = input_samples - input_stride;
        for (horizontal_index = half_decim_step, decim_horizontal_index = 0;
             horizontal_index < input_area_width;
             horizontal_index += decim_step, decim_horizontal_index++) {
            uint32_t sum = (uint32_t)prev_input_line[horizontal_index - 1] +
                           (uint32_t)prev_input_line[horizontal_index] +
                           (uint32_t)input_samples[horizontal_index - 1] +
                           (uint32_t)input_samples[horizontal_index];
            decim_samples[decim_horizontal_index] = (sum + 2) >> 2;
        }
        input_samples += input_stripe_stride;
        decim_samples += decim_stride;
    }

    return;
}

/********************************************
* calculate_histogram
*      creates n-bins histogram for the input
********************************************/
void calculate_histogram(uint8_t * input_samples, // input parameter, input samples Ptr
                         uint32_t  input_area_width, // input parameter, input area width
                         uint32_t  input_area_height, // input parameter, input area height
                         uint32_t  stride, // input parameter, input stride
                         uint8_t   decim_step, // input parameter, area height
                         uint32_t *histogram, // output parameter, output histogram
                         uint64_t *sum) {
    uint32_t horizontal_index;
    uint32_t vertical_index;
    *sum = 0;

    for (vertical_index = 0; vertical_index < input_area_height; vertical_index += decim_step) {
        for (horizontal_index = 0; horizontal_index < input_area_width;
             horizontal_index += decim_step) {
            ++(histogram[input_samples[horizontal_index]]);
            *sum += input_samples[horizontal_index];
        }
        input_samples += (stride << (decim_step >> 1));
    }

    return;
}
/*******************************************
 * compute_mean
 *   returns the mean of a block
 *******************************************/
uint64_t compute_mean_c(uint8_t *input_samples, /**< input parameter, input samples Ptr */
                        uint32_t input_stride, /**< input parameter, input stride */
                        uint32_t input_area_width, /**< input parameter, input area width */
                        uint32_t input_area_height) /**< input parameter, input area height */
{
    uint32_t hi, vi;
    uint64_t block_mean = 0;
    assert(input_area_width > 0);
    assert(input_area_height > 0);

    for (vi = 0; vi < input_area_height; vi++) {
        for (hi = 0; hi < input_area_width; hi++) { block_mean += input_samples[hi]; }
        input_samples += input_stride;
    }

    block_mean = (block_mean << (VARIANCE_PRECISION >> 1)) / (input_area_width * input_area_height);

    return block_mean;
}

/*******************************************
 * compute_mean_squared_values_c
 *   returns the Mean of Squared Values
 *******************************************/
uint64_t compute_mean_squared_values_c(
    uint8_t *input_samples, /**< input parameter, input samples Ptr */
    uint32_t input_stride, /**< input parameter, input stride */
    uint32_t input_area_width, /**< input parameter, input area width */
    uint32_t input_area_height) /**< input parameter, input area height */
{
    uint32_t hi, vi;
    uint64_t block_mean = 0;
    assert(input_area_width > 0);
    assert(input_area_height > 0);

    for (vi = 0; vi < input_area_height; vi++) {
        for (hi = 0; hi < input_area_width; hi++) {
            block_mean += input_samples[hi] * input_samples[hi];
        }
        input_samples += input_stride;
    }

    block_mean = (block_mean << VARIANCE_PRECISION) / (input_area_width * input_area_height);

    return block_mean;
}

uint64_t compute_sub_mean_8x8_c(uint8_t *input_samples, /**< input parameter, input samples Ptr */
                            uint16_t input_stride)/**< input parameter, input stride */
{
    uint32_t hi, vi;
    uint64_t block_mean = 0;
    uint16_t skip       = 0;

    for (vi = 0; skip < 8; skip = vi + vi) {
        for (hi = 0; hi < 8; hi++) { block_mean += input_samples[hi]; }
        input_samples += 2 * input_stride;
        vi++;
    }

    block_mean = block_mean << 3; // (VARIANCE_PRECISION >> 1)) /
        // (input_area_width * input_area_height/2)

    return block_mean;
}

uint64_t compute_sub_mean_squared_values_c(
    uint8_t *input_samples, /**< input parameter, input samples Ptr */
    uint32_t input_stride, /**< input parameter, input stride */
    uint32_t input_area_width, /**< input parameter, input area width */
    uint32_t input_area_height) /**< input parameter, input area height */
{
    uint32_t hi, vi;
    uint64_t block_mean = 0;
    uint16_t skip       = 0;

    for (vi = 0; skip < input_area_height; skip = vi + vi) {
        for (hi = 0; hi < input_area_width; hi++) {
            block_mean += input_samples[hi] * input_samples[hi];
        }
        input_samples += 2 * input_stride;
        vi++;
    }

    block_mean = block_mean << 11; // VARIANCE_PRECISION) / (input_area_width * input_area_height);

    return block_mean;
}

void compute_interm_var_four8x8_c(uint8_t *input_samples, uint16_t input_stride,
                                  uint64_t *mean_of8x8_blocks, // mean of four  8x8
                                  uint64_t *mean_of_squared8x8_blocks) // meanSquared
{
    uint32_t block_index = 0;
    // (0,1)
    mean_of8x8_blocks[0] = compute_sub_mean_8x8_c(input_samples + block_index, input_stride);
    mean_of_squared8x8_blocks[0] =
        compute_sub_mean_squared_values_c(input_samples + block_index, input_stride, 8, 8);

    // (0,2)
    block_index          = block_index + 8;
    mean_of8x8_blocks[1] = compute_sub_mean_8x8_c(input_samples + block_index, input_stride);
    mean_of_squared8x8_blocks[1] =
        compute_sub_mean_squared_values_c(input_samples + block_index, input_stride, 8, 8);

    // (0,3)
    block_index          = block_index + 8;
    mean_of8x8_blocks[2] = compute_sub_mean_8x8_c(input_samples + block_index, input_stride);
    mean_of_squared8x8_blocks[2] =
        compute_sub_mean_squared_values_c(input_samples + block_index, input_stride, 8, 8);

    // (0,4)
    block_index          = block_index + 8;
    mean_of8x8_blocks[3] = compute_sub_mean_8x8_c(input_samples + block_index, input_stride);
    mean_of_squared8x8_blocks[3] =
        compute_sub_mean_squared_values_c(input_samples + block_index, input_stride, 8, 8);
}

uint8_t get_filtered_types(uint8_t *ptr, uint32_t stride, uint8_t filter_type) {
    uint8_t *p = ptr - 1 - stride;

    uint32_t a = 0;

    if (filter_type == 0) {
        //Luma
        a = (p[1] + p[0 + stride] + 4 * p[1 + stride] + p[2 + stride] + p[1 + 2 * stride]) / 8;
    } else if (filter_type == 1) {
        a = (2 * p[1] + 2 * p[0 + stride] + 4 * p[1 + stride] + 2 * p[2 + stride] +
             2 * p[1 + 2 * stride]);

        a = (((uint32_t)((a * 2730) >> 14) + 1) >> 1) & 0xFFFF;

        //fixed point version of a=a/12 to mimic x86 instruction _mm256_mulhrs_epi16;
        //a= (a*2730)>>15;
    } else if (filter_type == 2) {
        a = (4 * p[1] + 4 * p[0 + stride] + 4 * p[1 + stride] + 4 * p[2 + stride] +
             4 * p[1 + 2 * stride]) /
            20;
    } else if (filter_type == 3) {
        a = (1 * p[0] + 1 * p[1] + 1 * p[2] + 1 * p[0 + stride] + 4 * p[1 + stride] +
             1 * p[2 + stride] + 1 * p[0 + 2 * stride] + 1 * p[1 + 2 * stride] +
             1 * p[2 + 2 * stride]) /
            12;
    } else if (filter_type == 4) {
        //gaussian matrix(Chroma)
        a = (1 * p[0] + 2 * p[1] + 1 * p[2] + 2 * p[0 + stride] + 4 * p[1 + stride] +
             2 * p[2 + stride] + 1 * p[0 + 2 * stride] + 2 * p[1 + 2 * stride] +
             1 * p[2 + 2 * stride]) /
            16;
    } else if (filter_type == 5) {
        a = (2 * p[0] + 2 * p[1] + 2 * p[2] + 2 * p[0 + stride] + 4 * p[1 + stride] +
             2 * p[2 + stride] + 2 * p[0 + 2 * stride] + 2 * p[1 + 2 * stride] +
             2 * p[2 + 2 * stride]) /
            20;
    } else if (filter_type == 6) {
        a = (4 * p[0] + 4 * p[1] + 4 * p[2] + 4 * p[0 + stride] + 4 * p[1 + stride] +
             4 * p[2 + stride] + 4 * p[0 + 2 * stride] + 4 * p[1 + 2 * stride] +
             4 * p[2 + 2 * stride]) /
            36;
    }

    return (uint8_t)CLIP3EQ(0, 255, a);
}

EbErrorType zero_out_chroma_block_mean(
    PictureParentControlSet *pcs_ptr, // input parameter, Picture Control Set Ptr
    uint32_t                 sb_coding_order // input parameter, SB address
) {
    EbErrorType return_error = EB_ErrorNone;
    // 16x16 mean
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_0]  = 0;
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_1]  = 0;
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_2]  = 0;
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_3]  = 0;
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_4]  = 0;
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_5]  = 0;
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_6]  = 0;
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_7]  = 0;
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_8]  = 0;
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_9]  = 0;
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_10] = 0;
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_11] = 0;
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_12] = 0;
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_13] = 0;
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_14] = 0;
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_15] = 0;

    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_0]  = 0;
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_1]  = 0;
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_2]  = 0;
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_3]  = 0;
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_4]  = 0;
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_5]  = 0;
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_6]  = 0;
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_7]  = 0;
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_8]  = 0;
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_9]  = 0;
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_10] = 0;
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_11] = 0;
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_12] = 0;
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_13] = 0;
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_14] = 0;
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_15] = 0;

    // 32x32 mean
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_32x32_0] = 0;
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_32x32_1] = 0;
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_32x32_2] = 0;
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_32x32_3] = 0;

    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_32x32_0] = 0;
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_32x32_1] = 0;
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_32x32_2] = 0;
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_32x32_3] = 0;

    // 64x64 mean
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_64x64] = 0;
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_64x64] = 0;

    return return_error;
}
/*******************************************
* compute_chroma_block_mean
*   computes the chroma block mean for 64x64, 32x32 and 16x16 CUs inside the tree block
*******************************************/
EbErrorType compute_chroma_block_mean(
    SequenceControlSet *     scs_ptr,
    PictureParentControlSet *pcs_ptr, // input parameter, Picture Control Set Ptr
    EbPictureBufferDesc *    input_padded_picture_ptr, // input parameter, Input Padded Picture
    uint32_t                 sb_coding_order, // input parameter, SB address
    uint32_t
        input_cb_origin_index, // input parameter, SB index, used to point to source/reference samples
    uint32_t
        input_cr_origin_index) // input parameter, SB index, used to point to source/reference samples
{
    EbErrorType return_error = EB_ErrorNone;

    uint32_t cb_block_index, cr_block_index;

    uint64_t cb_mean_of_16x16_blocks[16];
    uint64_t cr_mean_of_16x16_blocks[16];

    uint64_t cb_mean_of_32x32_blocks[4];
    uint64_t cr_mean_of_32x32_blocks[4];

    uint64_t cb_mean_of_64x64_blocks;
    uint64_t cr_mean_of_64x64_blocks;

    // (0,0) 16x16 block
    cb_block_index = input_cb_origin_index;
    cr_block_index = input_cr_origin_index;
    if (scs_ptr->block_mean_calc_prec == BLOCK_MEAN_PREC_FULL) {
        cb_mean_of_16x16_blocks[0] =
            compute_mean_8x8(&(input_padded_picture_ptr->buffer_cb[cb_block_index]),
                             input_padded_picture_ptr->stride_cb,
                             8,
                             8);
        cr_mean_of_16x16_blocks[0] =
            compute_mean_8x8(&(input_padded_picture_ptr->buffer_cr[cr_block_index]),
                             input_padded_picture_ptr->stride_cr,
                             8,
                             8);

        // (0,1)
        cb_block_index = cb_block_index + 8;
        cr_block_index = cr_block_index + 8;
        cb_mean_of_16x16_blocks[1] =
            compute_mean_8x8(&(input_padded_picture_ptr->buffer_cb[cb_block_index]),
                             input_padded_picture_ptr->stride_cb,
                             8,
                             8);
        cr_mean_of_16x16_blocks[1] =
            compute_mean_8x8(&(input_padded_picture_ptr->buffer_cr[cr_block_index]),
                             input_padded_picture_ptr->stride_cr,
                             8,
                             8);

        // (0,2)
        cb_block_index = cb_block_index + 8;
        cr_block_index = cr_block_index + 8;
        cb_mean_of_16x16_blocks[2] =
            compute_mean_8x8(&(input_padded_picture_ptr->buffer_cb[cb_block_index]),
                             input_padded_picture_ptr->stride_cb,
                             8,
                             8);
        cr_mean_of_16x16_blocks[2] =
            compute_mean_8x8(&(input_padded_picture_ptr->buffer_cr[cr_block_index]),
                             input_padded_picture_ptr->stride_cr,
                             8,
                             8);

        // (0,3)
        cb_block_index = cb_block_index + 8;
        cr_block_index = cr_block_index + 8;
        cb_mean_of_16x16_blocks[3] =
            compute_mean_8x8(&(input_padded_picture_ptr->buffer_cb[cb_block_index]),
                             input_padded_picture_ptr->stride_cb,
                             8,
                             8);
        cr_mean_of_16x16_blocks[3] =
            compute_mean_8x8(&(input_padded_picture_ptr->buffer_cr[cr_block_index]),
                             input_padded_picture_ptr->stride_cr,
                             8,
                             8);

        // (1,0)
        cb_block_index = input_cb_origin_index + (input_padded_picture_ptr->stride_cb << 3);
        cr_block_index = input_cr_origin_index + (input_padded_picture_ptr->stride_cr << 3);
        cb_mean_of_16x16_blocks[4] =
            compute_mean_8x8(&(input_padded_picture_ptr->buffer_cb[cb_block_index]),
                             input_padded_picture_ptr->stride_cb,
                             8,
                             8);
        cr_mean_of_16x16_blocks[4] =
            compute_mean_8x8(&(input_padded_picture_ptr->buffer_cr[cr_block_index]),
                             input_padded_picture_ptr->stride_cr,
                             8,
                             8);

        // (1,1)
        cb_block_index = cb_block_index + 8;
        cr_block_index = cr_block_index + 8;
        cb_mean_of_16x16_blocks[5] =
            compute_mean_8x8(&(input_padded_picture_ptr->buffer_cb[cb_block_index]),
                             input_padded_picture_ptr->stride_cb,
                             8,
                             8);
        cr_mean_of_16x16_blocks[5] =
            compute_mean_8x8(&(input_padded_picture_ptr->buffer_cr[cr_block_index]),
                             input_padded_picture_ptr->stride_cr,
                             8,
                             8);

        // (1,2)
        cb_block_index = cb_block_index + 8;
        cr_block_index = cr_block_index + 8;
        cb_mean_of_16x16_blocks[6] =
            compute_mean_8x8(&(input_padded_picture_ptr->buffer_cb[cb_block_index]),
                             input_padded_picture_ptr->stride_cb,
                             8,
                             8);
        cr_mean_of_16x16_blocks[6] =
            compute_mean_8x8(&(input_padded_picture_ptr->buffer_cr[cr_block_index]),
                             input_padded_picture_ptr->stride_cr,
                             8,
                             8);

        // (1,3)
        cb_block_index = cb_block_index + 8;
        cr_block_index = cr_block_index + 8;
        cb_mean_of_16x16_blocks[7] =
            compute_mean_8x8(&(input_padded_picture_ptr->buffer_cb[cb_block_index]),
                             input_padded_picture_ptr->stride_cb,
                             8,
                             8);
        cr_mean_of_16x16_blocks[7] =
            compute_mean_8x8(&(input_padded_picture_ptr->buffer_cr[cr_block_index]),
                             input_padded_picture_ptr->stride_cr,
                             8,
                             8);

        // (2,0)
        cb_block_index = input_cb_origin_index + (input_padded_picture_ptr->stride_cb << 4);
        cr_block_index = input_cr_origin_index + (input_padded_picture_ptr->stride_cr << 4);
        cb_mean_of_16x16_blocks[8] =
            compute_mean_8x8(&(input_padded_picture_ptr->buffer_cb[cb_block_index]),
                             input_padded_picture_ptr->stride_cb,
                             8,
                             8);
        cr_mean_of_16x16_blocks[8] =
            compute_mean_8x8(&(input_padded_picture_ptr->buffer_cr[cr_block_index]),
                             input_padded_picture_ptr->stride_cr,
                             8,
                             8);

        // (2,1)
        cb_block_index = cb_block_index + 8;
        cr_block_index = cr_block_index + 8;
        cb_mean_of_16x16_blocks[9] =
            compute_mean_8x8(&(input_padded_picture_ptr->buffer_cb[cb_block_index]),
                             input_padded_picture_ptr->stride_cb,
                             8,
                             8);
        cr_mean_of_16x16_blocks[9] =
            compute_mean_8x8(&(input_padded_picture_ptr->buffer_cr[cr_block_index]),
                             input_padded_picture_ptr->stride_cr,
                             8,
                             8);

        // (2,2)
        cb_block_index = cb_block_index + 8;
        cr_block_index = cr_block_index + 8;
        cb_mean_of_16x16_blocks[10] =
            compute_mean_8x8(&(input_padded_picture_ptr->buffer_cb[cb_block_index]),
                             input_padded_picture_ptr->stride_cb,
                             8,
                             8);
        cr_mean_of_16x16_blocks[10] =
            compute_mean_8x8(&(input_padded_picture_ptr->buffer_cr[cr_block_index]),
                             input_padded_picture_ptr->stride_cr,
                             8,
                             8);

        // (2,3)
        cb_block_index = cb_block_index + 8;
        cr_block_index = cr_block_index + 8;
        cb_mean_of_16x16_blocks[11] =
            compute_mean_8x8(&(input_padded_picture_ptr->buffer_cb[cb_block_index]),
                             input_padded_picture_ptr->stride_cb,
                             8,
                             8);
        cr_mean_of_16x16_blocks[11] =
            compute_mean_8x8(&(input_padded_picture_ptr->buffer_cr[cr_block_index]),
                             input_padded_picture_ptr->stride_cr,
                             8,
                             8);

        // (3,0)
        cb_block_index = input_cb_origin_index + (input_padded_picture_ptr->stride_cb * 24);
        cr_block_index = input_cr_origin_index + (input_padded_picture_ptr->stride_cr * 24);
        cb_mean_of_16x16_blocks[12] =
            compute_mean_8x8(&(input_padded_picture_ptr->buffer_cb[cb_block_index]),
                             input_padded_picture_ptr->stride_cb,
                             8,
                             8);
        cr_mean_of_16x16_blocks[12] =
            compute_mean_8x8(&(input_padded_picture_ptr->buffer_cr[cr_block_index]),
                             input_padded_picture_ptr->stride_cr,
                             8,
                             8);

        // (3,1)
        cb_block_index = cb_block_index + 8;
        cr_block_index = cr_block_index + 8;
        cb_mean_of_16x16_blocks[13] =
            compute_mean_8x8(&(input_padded_picture_ptr->buffer_cb[cb_block_index]),
                             input_padded_picture_ptr->stride_cb,
                             8,
                             8);
        cr_mean_of_16x16_blocks[13] =
            compute_mean_8x8(&(input_padded_picture_ptr->buffer_cr[cr_block_index]),
                             input_padded_picture_ptr->stride_cr,
                             8,
                             8);

        // (3,2)
        cb_block_index = cb_block_index + 8;
        cr_block_index = cr_block_index + 8;
        cb_mean_of_16x16_blocks[14] =
            compute_mean_8x8(&(input_padded_picture_ptr->buffer_cb[cb_block_index]),
                             input_padded_picture_ptr->stride_cb,
                             8,
                             8);
        cr_mean_of_16x16_blocks[14] =
            compute_mean_8x8(&(input_padded_picture_ptr->buffer_cr[cr_block_index]),
                             input_padded_picture_ptr->stride_cr,
                             8,
                             8);

        // (3,3)
        cb_block_index = cb_block_index + 8;
        cr_block_index = cr_block_index + 8;
        cb_mean_of_16x16_blocks[15] =
            compute_mean_8x8(&(input_padded_picture_ptr->buffer_cb[cb_block_index]),
                             input_padded_picture_ptr->stride_cb,
                             8,
                             8);
        cr_mean_of_16x16_blocks[15] =
            compute_mean_8x8(&(input_padded_picture_ptr->buffer_cr[cr_block_index]),
                             input_padded_picture_ptr->stride_cr,
                             8,
                             8);
    } else {
        const uint16_t stride_cb = input_padded_picture_ptr->stride_cb;
        const uint16_t stride_cr = input_padded_picture_ptr->stride_cr;

        cb_mean_of_16x16_blocks[0] = compute_sub_mean_8x8(
            &(input_padded_picture_ptr->buffer_cb[cb_block_index]), stride_cb);
        cr_mean_of_16x16_blocks[0] = compute_sub_mean_8x8(
            &(input_padded_picture_ptr->buffer_cr[cr_block_index]), stride_cr);

        // (0,1)
        cb_block_index             = cb_block_index + 8;
        cr_block_index             = cr_block_index + 8;
        cb_mean_of_16x16_blocks[1] = compute_sub_mean_8x8(
            &(input_padded_picture_ptr->buffer_cb[cb_block_index]), stride_cb);
        cr_mean_of_16x16_blocks[1] = compute_sub_mean_8x8(
            &(input_padded_picture_ptr->buffer_cr[cr_block_index]), stride_cr);

        // (0,2)
        cb_block_index             = cb_block_index + 8;
        cr_block_index             = cr_block_index + 8;
        cb_mean_of_16x16_blocks[2] = compute_sub_mean_8x8(
            &(input_padded_picture_ptr->buffer_cb[cb_block_index]), stride_cb);
        cr_mean_of_16x16_blocks[2] = compute_sub_mean_8x8(
            &(input_padded_picture_ptr->buffer_cr[cr_block_index]), stride_cr);

        // (0,3)
        cb_block_index             = cb_block_index + 8;
        cr_block_index             = cr_block_index + 8;
        cb_mean_of_16x16_blocks[3] = compute_sub_mean_8x8(
            &(input_padded_picture_ptr->buffer_cb[cb_block_index]), stride_cb);
        cr_mean_of_16x16_blocks[3] = compute_sub_mean_8x8(
            &(input_padded_picture_ptr->buffer_cr[cr_block_index]), stride_cr);

        // (1,0)
        cb_block_index             = input_cb_origin_index + (stride_cb << 3);
        cr_block_index             = input_cr_origin_index + (stride_cr << 3);
        cb_mean_of_16x16_blocks[4] = compute_sub_mean_8x8(
            &(input_padded_picture_ptr->buffer_cb[cb_block_index]), stride_cb);
        cr_mean_of_16x16_blocks[4] = compute_sub_mean_8x8(
            &(input_padded_picture_ptr->buffer_cr[cr_block_index]), stride_cr);

        // (1,1)
        cb_block_index             = cb_block_index + 8;
        cr_block_index             = cr_block_index + 8;
        cb_mean_of_16x16_blocks[5] = compute_sub_mean_8x8(
            &(input_padded_picture_ptr->buffer_cb[cb_block_index]), stride_cb);
        cr_mean_of_16x16_blocks[5] = compute_sub_mean_8x8(
            &(input_padded_picture_ptr->buffer_cr[cr_block_index]), stride_cr);

        // (1,2)
        cb_block_index             = cb_block_index + 8;
        cr_block_index             = cr_block_index + 8;
        cb_mean_of_16x16_blocks[6] = compute_sub_mean_8x8(
            &(input_padded_picture_ptr->buffer_cb[cb_block_index]), stride_cb);
        cr_mean_of_16x16_blocks[6] = compute_sub_mean_8x8(
            &(input_padded_picture_ptr->buffer_cr[cr_block_index]), stride_cr);

        // (1,3)
        cb_block_index             = cb_block_index + 8;
        cr_block_index             = cr_block_index + 8;
        cb_mean_of_16x16_blocks[7] = compute_sub_mean_8x8(
            &(input_padded_picture_ptr->buffer_cb[cb_block_index]), stride_cb);
        cr_mean_of_16x16_blocks[7] = compute_sub_mean_8x8(
            &(input_padded_picture_ptr->buffer_cr[cr_block_index]), stride_cr);

        // (2,0)
        cb_block_index             = input_cb_origin_index + (stride_cb << 4);
        cr_block_index             = input_cr_origin_index + (stride_cr << 4);
        cb_mean_of_16x16_blocks[8] = compute_sub_mean_8x8(
            &(input_padded_picture_ptr->buffer_cb[cb_block_index]), stride_cb);
        cr_mean_of_16x16_blocks[8] = compute_sub_mean_8x8(
            &(input_padded_picture_ptr->buffer_cr[cr_block_index]), stride_cr);

        // (2,1)
        cb_block_index             = cb_block_index + 8;
        cr_block_index             = cr_block_index + 8;
        cb_mean_of_16x16_blocks[9] = compute_sub_mean_8x8(
            &(input_padded_picture_ptr->buffer_cb[cb_block_index]), stride_cb);
        cr_mean_of_16x16_blocks[9] = compute_sub_mean_8x8(
            &(input_padded_picture_ptr->buffer_cr[cr_block_index]), stride_cr);

        // (2,2)
        cb_block_index              = cb_block_index + 8;
        cr_block_index              = cr_block_index + 8;
        cb_mean_of_16x16_blocks[10] = compute_sub_mean_8x8(
            &(input_padded_picture_ptr->buffer_cb[cb_block_index]), stride_cb);
        cr_mean_of_16x16_blocks[10] = compute_sub_mean_8x8(
            &(input_padded_picture_ptr->buffer_cr[cr_block_index]), stride_cr);

        // (2,3)
        cb_block_index              = cb_block_index + 8;
        cr_block_index              = cr_block_index + 8;
        cb_mean_of_16x16_blocks[11] = compute_sub_mean_8x8(
            &(input_padded_picture_ptr->buffer_cb[cb_block_index]), stride_cb);
        cr_mean_of_16x16_blocks[11] = compute_sub_mean_8x8(
            &(input_padded_picture_ptr->buffer_cr[cr_block_index]), stride_cr);

        // (3,0)
        cb_block_index              = input_cb_origin_index + (stride_cb * 24);
        cr_block_index              = input_cr_origin_index + (stride_cr * 24);
        cb_mean_of_16x16_blocks[12] = compute_sub_mean_8x8(
            &(input_padded_picture_ptr->buffer_cb[cb_block_index]), stride_cb);
        cr_mean_of_16x16_blocks[12] = compute_sub_mean_8x8(
            &(input_padded_picture_ptr->buffer_cr[cr_block_index]), stride_cr);

        // (3,1)
        cb_block_index              = cb_block_index + 8;
        cr_block_index              = cr_block_index + 8;
        cb_mean_of_16x16_blocks[13] = compute_sub_mean_8x8(
            &(input_padded_picture_ptr->buffer_cb[cb_block_index]), stride_cb);
        cr_mean_of_16x16_blocks[13] = compute_sub_mean_8x8(
            &(input_padded_picture_ptr->buffer_cr[cr_block_index]), stride_cr);

        // (3,2)
        cb_block_index              = cb_block_index + 8;
        cr_block_index              = cr_block_index + 8;
        cb_mean_of_16x16_blocks[14] = compute_sub_mean_8x8(
            &(input_padded_picture_ptr->buffer_cb[cb_block_index]), stride_cb);
        cr_mean_of_16x16_blocks[14] = compute_sub_mean_8x8(
            &(input_padded_picture_ptr->buffer_cr[cr_block_index]), stride_cr);

        // (3,3)
        cb_block_index              = cb_block_index + 8;
        cr_block_index              = cr_block_index + 8;
        cb_mean_of_16x16_blocks[15] = compute_sub_mean_8x8(
            &(input_padded_picture_ptr->buffer_cb[cb_block_index]), stride_cb);
        cr_mean_of_16x16_blocks[15] = compute_sub_mean_8x8(
            &(input_padded_picture_ptr->buffer_cr[cr_block_index]), stride_cr);
    }

    // 32x32
    cb_mean_of_32x32_blocks[0] = (cb_mean_of_16x16_blocks[0] + cb_mean_of_16x16_blocks[1] +
                                  cb_mean_of_16x16_blocks[4] + cb_mean_of_16x16_blocks[5]) >>
                                 2;
    cr_mean_of_32x32_blocks[0] = (cr_mean_of_16x16_blocks[0] + cr_mean_of_16x16_blocks[1] +
                                  cr_mean_of_16x16_blocks[4] + cr_mean_of_16x16_blocks[5]) >>
                                 2;

    cb_mean_of_32x32_blocks[1] = (cb_mean_of_16x16_blocks[2] + cb_mean_of_16x16_blocks[3] +
                                  cb_mean_of_16x16_blocks[6] + cb_mean_of_16x16_blocks[7]) >>
                                 2;
    cr_mean_of_32x32_blocks[1] = (cr_mean_of_16x16_blocks[2] + cr_mean_of_16x16_blocks[3] +
                                  cr_mean_of_16x16_blocks[6] + cr_mean_of_16x16_blocks[7]) >>
                                 2;

    cb_mean_of_32x32_blocks[2] = (cb_mean_of_16x16_blocks[8] + cb_mean_of_16x16_blocks[9] +
                                  cb_mean_of_16x16_blocks[12] + cb_mean_of_16x16_blocks[13]) >>
                                 2;
    cr_mean_of_32x32_blocks[2] = (cr_mean_of_16x16_blocks[8] + cr_mean_of_16x16_blocks[9] +
                                  cr_mean_of_16x16_blocks[12] + cr_mean_of_16x16_blocks[13]) >>
                                 2;

    cb_mean_of_32x32_blocks[3] = (cb_mean_of_16x16_blocks[10] + cb_mean_of_16x16_blocks[11] +
                                  cb_mean_of_16x16_blocks[14] + cb_mean_of_16x16_blocks[15]) >>
                                 2;
    cr_mean_of_32x32_blocks[3] = (cr_mean_of_16x16_blocks[10] + cr_mean_of_16x16_blocks[11] +
                                  cr_mean_of_16x16_blocks[14] + cr_mean_of_16x16_blocks[15]) >>
                                 2;

    // 64x64
    cb_mean_of_64x64_blocks = (cb_mean_of_32x32_blocks[0] + cb_mean_of_32x32_blocks[1] +
                               cb_mean_of_32x32_blocks[3] + cb_mean_of_32x32_blocks[3]) >>
                              2;
    cr_mean_of_64x64_blocks = (cr_mean_of_32x32_blocks[0] + cr_mean_of_32x32_blocks[1] +
                               cr_mean_of_32x32_blocks[3] + cr_mean_of_32x32_blocks[3]) >>
                              2;
    // 16x16 mean
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_0] =
        (uint8_t)(cb_mean_of_16x16_blocks[0] >> MEAN_PRECISION);
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_1] =
        (uint8_t)(cb_mean_of_16x16_blocks[1] >> MEAN_PRECISION);
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_2] =
        (uint8_t)(cb_mean_of_16x16_blocks[2] >> MEAN_PRECISION);
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_3] =
        (uint8_t)(cb_mean_of_16x16_blocks[3] >> MEAN_PRECISION);
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_4] =
        (uint8_t)(cb_mean_of_16x16_blocks[4] >> MEAN_PRECISION);
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_5] =
        (uint8_t)(cb_mean_of_16x16_blocks[5] >> MEAN_PRECISION);
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_6] =
        (uint8_t)(cb_mean_of_16x16_blocks[6] >> MEAN_PRECISION);
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_7] =
        (uint8_t)(cb_mean_of_16x16_blocks[7] >> MEAN_PRECISION);
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_8] =
        (uint8_t)(cb_mean_of_16x16_blocks[8] >> MEAN_PRECISION);
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_9] =
        (uint8_t)(cb_mean_of_16x16_blocks[9] >> MEAN_PRECISION);
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_10] =
        (uint8_t)(cb_mean_of_16x16_blocks[10] >> MEAN_PRECISION);
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_11] =
        (uint8_t)(cb_mean_of_16x16_blocks[11] >> MEAN_PRECISION);
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_12] =
        (uint8_t)(cb_mean_of_16x16_blocks[12] >> MEAN_PRECISION);
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_13] =
        (uint8_t)(cb_mean_of_16x16_blocks[13] >> MEAN_PRECISION);
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_14] =
        (uint8_t)(cb_mean_of_16x16_blocks[14] >> MEAN_PRECISION);
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_15] =
        (uint8_t)(cb_mean_of_16x16_blocks[15] >> MEAN_PRECISION);

    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_0] =
        (uint8_t)(cr_mean_of_16x16_blocks[0] >> MEAN_PRECISION);
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_1] =
        (uint8_t)(cr_mean_of_16x16_blocks[1] >> MEAN_PRECISION);
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_2] =
        (uint8_t)(cr_mean_of_16x16_blocks[2] >> MEAN_PRECISION);
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_3] =
        (uint8_t)(cr_mean_of_16x16_blocks[3] >> MEAN_PRECISION);
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_4] =
        (uint8_t)(cr_mean_of_16x16_blocks[4] >> MEAN_PRECISION);
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_5] =
        (uint8_t)(cr_mean_of_16x16_blocks[5] >> MEAN_PRECISION);
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_6] =
        (uint8_t)(cr_mean_of_16x16_blocks[6] >> MEAN_PRECISION);
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_7] =
        (uint8_t)(cr_mean_of_16x16_blocks[7] >> MEAN_PRECISION);
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_8] =
        (uint8_t)(cr_mean_of_16x16_blocks[8] >> MEAN_PRECISION);
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_9] =
        (uint8_t)(cr_mean_of_16x16_blocks[9] >> MEAN_PRECISION);
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_10] =
        (uint8_t)(cr_mean_of_16x16_blocks[10] >> MEAN_PRECISION);
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_11] =
        (uint8_t)(cr_mean_of_16x16_blocks[11] >> MEAN_PRECISION);
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_12] =
        (uint8_t)(cr_mean_of_16x16_blocks[12] >> MEAN_PRECISION);
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_13] =
        (uint8_t)(cr_mean_of_16x16_blocks[13] >> MEAN_PRECISION);
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_14] =
        (uint8_t)(cr_mean_of_16x16_blocks[14] >> MEAN_PRECISION);
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_16x16_15] =
        (uint8_t)(cr_mean_of_16x16_blocks[15] >> MEAN_PRECISION);

    // 32x32 mean
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_32x32_0] =
        (uint8_t)(cb_mean_of_32x32_blocks[0] >> MEAN_PRECISION);
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_32x32_1] =
        (uint8_t)(cb_mean_of_32x32_blocks[1] >> MEAN_PRECISION);
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_32x32_2] =
        (uint8_t)(cb_mean_of_32x32_blocks[2] >> MEAN_PRECISION);
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_32x32_3] =
        (uint8_t)(cb_mean_of_32x32_blocks[3] >> MEAN_PRECISION);

    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_32x32_0] =
        (uint8_t)(cr_mean_of_32x32_blocks[0] >> MEAN_PRECISION);
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_32x32_1] =
        (uint8_t)(cr_mean_of_32x32_blocks[1] >> MEAN_PRECISION);
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_32x32_2] =
        (uint8_t)(cr_mean_of_32x32_blocks[2] >> MEAN_PRECISION);
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_32x32_3] =
        (uint8_t)(cr_mean_of_32x32_blocks[3] >> MEAN_PRECISION);

    // 64x64 mean
    pcs_ptr->cb_mean[sb_coding_order][ME_TIER_ZERO_PU_64x64] =
        (uint8_t)(cb_mean_of_64x64_blocks >> MEAN_PRECISION);
    pcs_ptr->cr_mean[sb_coding_order][ME_TIER_ZERO_PU_64x64] =
        (uint8_t)(cr_mean_of_64x64_blocks >> MEAN_PRECISION);

    return return_error;
}

/*******************************************
* compute_block_mean_compute_variance
*   computes the variance and the block mean of all CUs inside the tree block
*******************************************/
EbErrorType compute_block_mean_compute_variance(
    SequenceControlSet *     scs_ptr,
    PictureParentControlSet *pcs_ptr, // input parameter, Picture Control Set Ptr
    EbPictureBufferDesc *    input_padded_picture_ptr, // input parameter, Input Padded Picture
    uint32_t                 sb_index, // input parameter, SB address
    uint32_t
        input_luma_origin_index) // input parameter, SB index, used to point to source/reference samples
{
    EbErrorType return_error = EB_ErrorNone;

    uint32_t block_index;

    uint64_t mean_of8x8_blocks[64];
    uint64_t mean_of_8x8_squared_values_blocks[64];

    uint64_t mean_of_16x16_blocks[16];
    uint64_t mean_of16x16_squared_values_blocks[16];

    uint64_t mean_of_32x32_blocks[4];
    uint64_t mean_of32x32_squared_values_blocks[4];

    uint64_t mean_of_64x64_blocks;
    uint64_t mean_of64x64_squared_values_blocks;

    // (0,0)
    block_index = input_luma_origin_index;
    if (scs_ptr->block_mean_calc_prec == BLOCK_MEAN_PREC_FULL) {
        mean_of8x8_blocks[0] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                input_padded_picture_ptr->stride_y,
                                                8,
                                                8);
        mean_of_8x8_squared_values_blocks[0] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (0,1)
        block_index          = block_index + 8;
        mean_of8x8_blocks[1] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                input_padded_picture_ptr->stride_y,
                                                8,
                                                8);
        mean_of_8x8_squared_values_blocks[1] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (0,2)
        block_index          = block_index + 8;
        mean_of8x8_blocks[2] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                input_padded_picture_ptr->stride_y,
                                                8,
                                                8);
        mean_of_8x8_squared_values_blocks[2] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (0,3)
        block_index          = block_index + 8;
        mean_of8x8_blocks[3] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                input_padded_picture_ptr->stride_y,
                                                8,
                                                8);
        mean_of_8x8_squared_values_blocks[3] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (0,4)
        block_index          = block_index + 8;
        mean_of8x8_blocks[4] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                input_padded_picture_ptr->stride_y,
                                                8,
                                                8);
        mean_of_8x8_squared_values_blocks[4] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (0,5)
        block_index          = block_index + 8;
        mean_of8x8_blocks[5] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                input_padded_picture_ptr->stride_y,
                                                8,
                                                8);
        mean_of_8x8_squared_values_blocks[5] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (0,6)
        block_index          = block_index + 8;
        mean_of8x8_blocks[6] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                input_padded_picture_ptr->stride_y,
                                                8,
                                                8);
        mean_of_8x8_squared_values_blocks[6] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (0,7)
        block_index          = block_index + 8;
        mean_of8x8_blocks[7] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                input_padded_picture_ptr->stride_y,
                                                8,
                                                8);
        mean_of_8x8_squared_values_blocks[7] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (1,0)
        block_index          = input_luma_origin_index + (input_padded_picture_ptr->stride_y << 3);
        mean_of8x8_blocks[8] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                input_padded_picture_ptr->stride_y,
                                                8,
                                                8);
        mean_of_8x8_squared_values_blocks[8] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (1,1)
        block_index          = block_index + 8;
        mean_of8x8_blocks[9] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                input_padded_picture_ptr->stride_y,
                                                8,
                                                8);
        mean_of_8x8_squared_values_blocks[9] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (1,2)
        block_index           = block_index + 8;
        mean_of8x8_blocks[10] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[10] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (1,3)
        block_index           = block_index + 8;
        mean_of8x8_blocks[11] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[11] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (1,4)
        block_index           = block_index + 8;
        mean_of8x8_blocks[12] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[12] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (1,5)
        block_index           = block_index + 8;
        mean_of8x8_blocks[13] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[13] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (1,6)
        block_index           = block_index + 8;
        mean_of8x8_blocks[14] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[14] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (1,7)
        block_index           = block_index + 8;
        mean_of8x8_blocks[15] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[15] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (2,0)
        block_index           = input_luma_origin_index + (input_padded_picture_ptr->stride_y << 4);
        mean_of8x8_blocks[16] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[16] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (2,1)
        block_index           = block_index + 8;
        mean_of8x8_blocks[17] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[17] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (2,2)
        block_index           = block_index + 8;
        mean_of8x8_blocks[18] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[18] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (2,3)
        block_index           = block_index + 8;
        mean_of8x8_blocks[19] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[19] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        /// (2,4)
        block_index           = block_index + 8;
        mean_of8x8_blocks[20] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[20] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (2,5)
        block_index           = block_index + 8;
        mean_of8x8_blocks[21] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[21] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (2,6)
        block_index           = block_index + 8;
        mean_of8x8_blocks[22] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[22] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (2,7)
        block_index           = block_index + 8;
        mean_of8x8_blocks[23] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[23] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (3,0)
        block_index = input_luma_origin_index + (input_padded_picture_ptr->stride_y << 3) +
                      (input_padded_picture_ptr->stride_y << 4);
        mean_of8x8_blocks[24] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[24] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (3,1)
        block_index           = block_index + 8;
        mean_of8x8_blocks[25] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[25] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (3,2)
        block_index           = block_index + 8;
        mean_of8x8_blocks[26] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[26] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (3,3)
        block_index           = block_index + 8;
        mean_of8x8_blocks[27] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[27] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (3,4)
        block_index           = block_index + 8;
        mean_of8x8_blocks[28] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[28] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (3,5)
        block_index           = block_index + 8;
        mean_of8x8_blocks[29] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[29] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (3,6)
        block_index           = block_index + 8;
        mean_of8x8_blocks[30] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[30] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (3,7)
        block_index           = block_index + 8;
        mean_of8x8_blocks[31] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[31] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (4,0)
        block_index           = input_luma_origin_index + (input_padded_picture_ptr->stride_y << 5);
        mean_of8x8_blocks[32] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[32] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (4,1)
        block_index           = block_index + 8;
        mean_of8x8_blocks[33] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[33] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (4,2)
        block_index           = block_index + 8;
        mean_of8x8_blocks[34] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[34] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (4,3)
        block_index           = block_index + 8;
        mean_of8x8_blocks[35] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[35] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (4,4)
        block_index           = block_index + 8;
        mean_of8x8_blocks[36] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[36] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (4,5)
        block_index           = block_index + 8;
        mean_of8x8_blocks[37] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[37] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (4,6)
        block_index           = block_index + 8;
        mean_of8x8_blocks[38] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[38] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (4,7)
        block_index           = block_index + 8;
        mean_of8x8_blocks[39] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[39] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (5,0)
        block_index = input_luma_origin_index + (input_padded_picture_ptr->stride_y << 3) +
                      (input_padded_picture_ptr->stride_y << 5);
        mean_of8x8_blocks[40] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[40] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (5,1)
        block_index           = block_index + 8;
        mean_of8x8_blocks[41] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[41] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (5,2)
        block_index           = block_index + 8;
        mean_of8x8_blocks[42] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[42] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (5,3)
        block_index           = block_index + 8;
        mean_of8x8_blocks[43] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[43] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (5,4)
        block_index           = block_index + 8;
        mean_of8x8_blocks[44] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[44] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (5,5)
        block_index           = block_index + 8;
        mean_of8x8_blocks[45] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[45] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (5,6)
        block_index           = block_index + 8;
        mean_of8x8_blocks[46] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[46] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (5,7)
        block_index           = block_index + 8;
        mean_of8x8_blocks[47] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[47] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (6,0)
        block_index = input_luma_origin_index + (input_padded_picture_ptr->stride_y << 4) +
                      (input_padded_picture_ptr->stride_y << 5);
        mean_of8x8_blocks[48] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[48] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (6,1)
        block_index           = block_index + 8;
        mean_of8x8_blocks[49] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[49] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (6,2)
        block_index           = block_index + 8;
        mean_of8x8_blocks[50] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[50] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (6,3)
        block_index           = block_index + 8;
        mean_of8x8_blocks[51] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[51] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (6,4)
        block_index           = block_index + 8;
        mean_of8x8_blocks[52] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[52] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (6,5)
        block_index           = block_index + 8;
        mean_of8x8_blocks[53] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[53] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (6,6)
        block_index           = block_index + 8;
        mean_of8x8_blocks[54] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[54] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (6,7)
        block_index           = block_index + 8;
        mean_of8x8_blocks[55] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[55] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (7,0)
        block_index = input_luma_origin_index + (input_padded_picture_ptr->stride_y << 3) +
                      (input_padded_picture_ptr->stride_y << 4) +
                      (input_padded_picture_ptr->stride_y << 5);
        mean_of8x8_blocks[56] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[56] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (7,1)
        block_index           = block_index + 8;
        mean_of8x8_blocks[57] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[57] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (7,2)
        block_index           = block_index + 8;
        mean_of8x8_blocks[58] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[58] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (7,3)
        block_index           = block_index + 8;
        mean_of8x8_blocks[59] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[59] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (7,4)
        block_index           = block_index + 8;
        mean_of8x8_blocks[60] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[60] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (7,5)
        block_index           = block_index + 8;
        mean_of8x8_blocks[61] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[61] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (7,6)
        block_index           = block_index + 8;
        mean_of8x8_blocks[62] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[62] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);

        // (7,7)
        block_index           = block_index + 8;
        mean_of8x8_blocks[63] = compute_mean_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                                 input_padded_picture_ptr->stride_y,
                                                 8,
                                                 8);
        mean_of_8x8_squared_values_blocks[63] =
            compute_mean_square_values_8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                           input_padded_picture_ptr->stride_y,
                                           8,
                                           8);
    } else {
        const uint16_t stride_y = input_padded_picture_ptr->stride_y;

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                   stride_y,
                                   &mean_of8x8_blocks[0],
                                   &mean_of_8x8_squared_values_blocks[0]);

        // (0,1)
        block_index = block_index + 32;

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                   stride_y,
                                   &mean_of8x8_blocks[4],
                                   &mean_of_8x8_squared_values_blocks[4]);

        // (0,5)
        block_index = block_index + 24;

        // (1,0)
        block_index = input_luma_origin_index + (stride_y << 3);

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                   stride_y,
                                   &mean_of8x8_blocks[8],
                                   &mean_of_8x8_squared_values_blocks[8]);

        // (1,1)
        block_index = block_index + 32;

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                   stride_y,
                                   &mean_of8x8_blocks[12],
                                   &mean_of_8x8_squared_values_blocks[12]);

        // (1,5)
        block_index = block_index + 24;

        // (2,0)
        block_index = input_luma_origin_index + (stride_y << 4);

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                   stride_y,
                                   &mean_of8x8_blocks[16],
                                   &mean_of_8x8_squared_values_blocks[16]);

        // (2,1)
        block_index = block_index + 32;

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                   stride_y,
                                   &mean_of8x8_blocks[20],
                                   &mean_of_8x8_squared_values_blocks[20]);

        // (2,5)
        block_index = block_index + 24;

        // (3,0)
        block_index = input_luma_origin_index + (stride_y << 3) + (stride_y << 4);

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                   stride_y,
                                   &mean_of8x8_blocks[24],
                                   &mean_of_8x8_squared_values_blocks[24]);

        // (3,1)
        block_index = block_index + 32;

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                   stride_y,
                                   &mean_of8x8_blocks[28],
                                   &mean_of_8x8_squared_values_blocks[28]);

        // (3,5)
        block_index = block_index + 24;

        // (4,0)
        block_index = input_luma_origin_index + (stride_y << 5);

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                   stride_y,
                                   &mean_of8x8_blocks[32],
                                   &mean_of_8x8_squared_values_blocks[32]);

        // (4,1)
        block_index = block_index + 32;

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                   stride_y,
                                   &mean_of8x8_blocks[36],
                                   &mean_of_8x8_squared_values_blocks[36]);

        // (4,5)
        block_index = block_index + 24;

        // (5,0)
        block_index = input_luma_origin_index + (stride_y << 3) + (stride_y << 5);

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                   stride_y,
                                   &mean_of8x8_blocks[40],
                                   &mean_of_8x8_squared_values_blocks[40]);

        // (5,1)
        block_index = block_index + 32;

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                   stride_y,
                                   &mean_of8x8_blocks[44],
                                   &mean_of_8x8_squared_values_blocks[44]);

        // (5,5)
        block_index = block_index + 24;

        // (6,0)
        block_index = input_luma_origin_index + (stride_y << 4) + (stride_y << 5);

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                   stride_y,
                                   &mean_of8x8_blocks[48],
                                   &mean_of_8x8_squared_values_blocks[48]);

        // (6,1)
        block_index = block_index + 32;

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                   stride_y,
                                   &mean_of8x8_blocks[52],
                                   &mean_of_8x8_squared_values_blocks[52]);

        // (6,5)
        block_index = block_index + 24;

        // (7,0)
        block_index = input_luma_origin_index + (stride_y << 3) + (stride_y << 4) + (stride_y << 5);

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                   stride_y,
                                   &mean_of8x8_blocks[56],
                                   &mean_of_8x8_squared_values_blocks[56]);

        // (7,1)
        block_index = block_index + 32;

        compute_interm_var_four8x8(&(input_padded_picture_ptr->buffer_y[block_index]),
                                   stride_y,
                                   &mean_of8x8_blocks[60],
                                   &mean_of_8x8_squared_values_blocks[60]);
    }

    // 16x16
    mean_of_16x16_blocks[0] = (mean_of8x8_blocks[0] + mean_of8x8_blocks[1] + mean_of8x8_blocks[8] +
                               mean_of8x8_blocks[9]) >>
                              2;
    mean_of_16x16_blocks[1] = (mean_of8x8_blocks[2] + mean_of8x8_blocks[3] + mean_of8x8_blocks[10] +
                               mean_of8x8_blocks[11]) >>
                              2;
    mean_of_16x16_blocks[2] = (mean_of8x8_blocks[4] + mean_of8x8_blocks[5] + mean_of8x8_blocks[12] +
                               mean_of8x8_blocks[13]) >>
                              2;
    mean_of_16x16_blocks[3] = (mean_of8x8_blocks[6] + mean_of8x8_blocks[7] + mean_of8x8_blocks[14] +
                               mean_of8x8_blocks[15]) >>
                              2;

    mean_of_16x16_blocks[4] = (mean_of8x8_blocks[16] + mean_of8x8_blocks[17] +
                               mean_of8x8_blocks[24] + mean_of8x8_blocks[25]) >>
                              2;
    mean_of_16x16_blocks[5] = (mean_of8x8_blocks[18] + mean_of8x8_blocks[19] +
                               mean_of8x8_blocks[26] + mean_of8x8_blocks[27]) >>
                              2;
    mean_of_16x16_blocks[6] = (mean_of8x8_blocks[20] + mean_of8x8_blocks[21] +
                               mean_of8x8_blocks[28] + mean_of8x8_blocks[29]) >>
                              2;
    mean_of_16x16_blocks[7] = (mean_of8x8_blocks[22] + mean_of8x8_blocks[23] +
                               mean_of8x8_blocks[30] + mean_of8x8_blocks[31]) >>
                              2;

    mean_of_16x16_blocks[8] = (mean_of8x8_blocks[32] + mean_of8x8_blocks[33] +
                               mean_of8x8_blocks[40] + mean_of8x8_blocks[41]) >>
                              2;
    mean_of_16x16_blocks[9] = (mean_of8x8_blocks[34] + mean_of8x8_blocks[35] +
                               mean_of8x8_blocks[42] + mean_of8x8_blocks[43]) >>
                              2;
    mean_of_16x16_blocks[10] = (mean_of8x8_blocks[36] + mean_of8x8_blocks[37] +
                                mean_of8x8_blocks[44] + mean_of8x8_blocks[45]) >>
                               2;
    mean_of_16x16_blocks[11] = (mean_of8x8_blocks[38] + mean_of8x8_blocks[39] +
                                mean_of8x8_blocks[46] + mean_of8x8_blocks[47]) >>
                               2;

    mean_of_16x16_blocks[12] = (mean_of8x8_blocks[48] + mean_of8x8_blocks[49] +
                                mean_of8x8_blocks[56] + mean_of8x8_blocks[57]) >>
                               2;
    mean_of_16x16_blocks[13] = (mean_of8x8_blocks[50] + mean_of8x8_blocks[51] +
                                mean_of8x8_blocks[58] + mean_of8x8_blocks[59]) >>
                               2;
    mean_of_16x16_blocks[14] = (mean_of8x8_blocks[52] + mean_of8x8_blocks[53] +
                                mean_of8x8_blocks[60] + mean_of8x8_blocks[61]) >>
                               2;
    mean_of_16x16_blocks[15] = (mean_of8x8_blocks[54] + mean_of8x8_blocks[55] +
                                mean_of8x8_blocks[62] + mean_of8x8_blocks[63]) >>
                               2;

    mean_of16x16_squared_values_blocks[0] =
        (mean_of_8x8_squared_values_blocks[0] + mean_of_8x8_squared_values_blocks[1] +
         mean_of_8x8_squared_values_blocks[8] + mean_of_8x8_squared_values_blocks[9]) >>
        2;
    mean_of16x16_squared_values_blocks[1] =
        (mean_of_8x8_squared_values_blocks[2] + mean_of_8x8_squared_values_blocks[3] +
         mean_of_8x8_squared_values_blocks[10] + mean_of_8x8_squared_values_blocks[11]) >>
        2;
    mean_of16x16_squared_values_blocks[2] =
        (mean_of_8x8_squared_values_blocks[4] + mean_of_8x8_squared_values_blocks[5] +
         mean_of_8x8_squared_values_blocks[12] + mean_of_8x8_squared_values_blocks[13]) >>
        2;
    mean_of16x16_squared_values_blocks[3] =
        (mean_of_8x8_squared_values_blocks[6] + mean_of_8x8_squared_values_blocks[7] +
         mean_of_8x8_squared_values_blocks[14] + mean_of_8x8_squared_values_blocks[15]) >>
        2;

    mean_of16x16_squared_values_blocks[4] =
        (mean_of_8x8_squared_values_blocks[16] + mean_of_8x8_squared_values_blocks[17] +
         mean_of_8x8_squared_values_blocks[24] + mean_of_8x8_squared_values_blocks[25]) >>
        2;
    mean_of16x16_squared_values_blocks[5] =
        (mean_of_8x8_squared_values_blocks[18] + mean_of_8x8_squared_values_blocks[19] +
         mean_of_8x8_squared_values_blocks[26] + mean_of_8x8_squared_values_blocks[27]) >>
        2;
    mean_of16x16_squared_values_blocks[6] =
        (mean_of_8x8_squared_values_blocks[20] + mean_of_8x8_squared_values_blocks[21] +
         mean_of_8x8_squared_values_blocks[28] + mean_of_8x8_squared_values_blocks[29]) >>
        2;
    mean_of16x16_squared_values_blocks[7] =
        (mean_of_8x8_squared_values_blocks[22] + mean_of_8x8_squared_values_blocks[23] +
         mean_of_8x8_squared_values_blocks[30] + mean_of_8x8_squared_values_blocks[31]) >>
        2;

    mean_of16x16_squared_values_blocks[8] =
        (mean_of_8x8_squared_values_blocks[32] + mean_of_8x8_squared_values_blocks[33] +
         mean_of_8x8_squared_values_blocks[40] + mean_of_8x8_squared_values_blocks[41]) >>
        2;
    mean_of16x16_squared_values_blocks[9] =
        (mean_of_8x8_squared_values_blocks[34] + mean_of_8x8_squared_values_blocks[35] +
         mean_of_8x8_squared_values_blocks[42] + mean_of_8x8_squared_values_blocks[43]) >>
        2;
    mean_of16x16_squared_values_blocks[10] =
        (mean_of_8x8_squared_values_blocks[36] + mean_of_8x8_squared_values_blocks[37] +
         mean_of_8x8_squared_values_blocks[44] + mean_of_8x8_squared_values_blocks[45]) >>
        2;
    mean_of16x16_squared_values_blocks[11] =
        (mean_of_8x8_squared_values_blocks[38] + mean_of_8x8_squared_values_blocks[39] +
         mean_of_8x8_squared_values_blocks[46] + mean_of_8x8_squared_values_blocks[47]) >>
        2;

    mean_of16x16_squared_values_blocks[12] =
        (mean_of_8x8_squared_values_blocks[48] + mean_of_8x8_squared_values_blocks[49] +
         mean_of_8x8_squared_values_blocks[56] + mean_of_8x8_squared_values_blocks[57]) >>
        2;
    mean_of16x16_squared_values_blocks[13] =
        (mean_of_8x8_squared_values_blocks[50] + mean_of_8x8_squared_values_blocks[51] +
         mean_of_8x8_squared_values_blocks[58] + mean_of_8x8_squared_values_blocks[59]) >>
        2;
    mean_of16x16_squared_values_blocks[14] =
        (mean_of_8x8_squared_values_blocks[52] + mean_of_8x8_squared_values_blocks[53] +
         mean_of_8x8_squared_values_blocks[60] + mean_of_8x8_squared_values_blocks[61]) >>
        2;
    mean_of16x16_squared_values_blocks[15] =
        (mean_of_8x8_squared_values_blocks[54] + mean_of_8x8_squared_values_blocks[55] +
         mean_of_8x8_squared_values_blocks[62] + mean_of_8x8_squared_values_blocks[63]) >>
        2;

    // 32x32
    mean_of_32x32_blocks[0] = (mean_of_16x16_blocks[0] + mean_of_16x16_blocks[1] +
                               mean_of_16x16_blocks[4] + mean_of_16x16_blocks[5]) >>
                              2;
    mean_of_32x32_blocks[1] = (mean_of_16x16_blocks[2] + mean_of_16x16_blocks[3] +
                               mean_of_16x16_blocks[6] + mean_of_16x16_blocks[7]) >>
                              2;
    mean_of_32x32_blocks[2] = (mean_of_16x16_blocks[8] + mean_of_16x16_blocks[9] +
                               mean_of_16x16_blocks[12] + mean_of_16x16_blocks[13]) >>
                              2;
    mean_of_32x32_blocks[3] = (mean_of_16x16_blocks[10] + mean_of_16x16_blocks[11] +
                               mean_of_16x16_blocks[14] + mean_of_16x16_blocks[15]) >>
                              2;

    mean_of32x32_squared_values_blocks[0] =
        (mean_of16x16_squared_values_blocks[0] + mean_of16x16_squared_values_blocks[1] +
         mean_of16x16_squared_values_blocks[4] + mean_of16x16_squared_values_blocks[5]) >>
        2;
    mean_of32x32_squared_values_blocks[1] =
        (mean_of16x16_squared_values_blocks[2] + mean_of16x16_squared_values_blocks[3] +
         mean_of16x16_squared_values_blocks[6] + mean_of16x16_squared_values_blocks[7]) >>
        2;
    mean_of32x32_squared_values_blocks[2] =
        (mean_of16x16_squared_values_blocks[8] + mean_of16x16_squared_values_blocks[9] +
         mean_of16x16_squared_values_blocks[12] + mean_of16x16_squared_values_blocks[13]) >>
        2;
    mean_of32x32_squared_values_blocks[3] =
        (mean_of16x16_squared_values_blocks[10] + mean_of16x16_squared_values_blocks[11] +
         mean_of16x16_squared_values_blocks[14] + mean_of16x16_squared_values_blocks[15]) >>
        2;

    // 64x64
    mean_of_64x64_blocks = (mean_of_32x32_blocks[0] + mean_of_32x32_blocks[1] +
                            mean_of_32x32_blocks[2] + mean_of_32x32_blocks[3]) >>
                           2;
    mean_of64x64_squared_values_blocks =
        (mean_of32x32_squared_values_blocks[0] + mean_of32x32_squared_values_blocks[1] +
         mean_of32x32_squared_values_blocks[2] + mean_of32x32_squared_values_blocks[3]) >>
        2;

    // 8x8 means
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_0] =
        (uint8_t)(mean_of8x8_blocks[0] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_1] =
        (uint8_t)(mean_of8x8_blocks[1] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_2] =
        (uint8_t)(mean_of8x8_blocks[2] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_3] =
        (uint8_t)(mean_of8x8_blocks[3] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_4] =
        (uint8_t)(mean_of8x8_blocks[4] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_5] =
        (uint8_t)(mean_of8x8_blocks[5] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_6] =
        (uint8_t)(mean_of8x8_blocks[6] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_7] =
        (uint8_t)(mean_of8x8_blocks[7] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_8] =
        (uint8_t)(mean_of8x8_blocks[8] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_9] =
        (uint8_t)(mean_of8x8_blocks[9] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_10] =
        (uint8_t)(mean_of8x8_blocks[10] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_11] =
        (uint8_t)(mean_of8x8_blocks[11] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_12] =
        (uint8_t)(mean_of8x8_blocks[12] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_13] =
        (uint8_t)(mean_of8x8_blocks[13] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_14] =
        (uint8_t)(mean_of8x8_blocks[14] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_15] =
        (uint8_t)(mean_of8x8_blocks[15] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_16] =
        (uint8_t)(mean_of8x8_blocks[16] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_17] =
        (uint8_t)(mean_of8x8_blocks[17] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_18] =
        (uint8_t)(mean_of8x8_blocks[18] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_19] =
        (uint8_t)(mean_of8x8_blocks[19] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_20] =
        (uint8_t)(mean_of8x8_blocks[20] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_21] =
        (uint8_t)(mean_of8x8_blocks[21] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_22] =
        (uint8_t)(mean_of8x8_blocks[22] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_23] =
        (uint8_t)(mean_of8x8_blocks[23] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_24] =
        (uint8_t)(mean_of8x8_blocks[24] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_25] =
        (uint8_t)(mean_of8x8_blocks[25] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_26] =
        (uint8_t)(mean_of8x8_blocks[26] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_27] =
        (uint8_t)(mean_of8x8_blocks[27] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_28] =
        (uint8_t)(mean_of8x8_blocks[28] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_29] =
        (uint8_t)(mean_of8x8_blocks[29] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_30] =
        (uint8_t)(mean_of8x8_blocks[30] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_31] =
        (uint8_t)(mean_of8x8_blocks[31] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_32] =
        (uint8_t)(mean_of8x8_blocks[32] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_33] =
        (uint8_t)(mean_of8x8_blocks[33] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_34] =
        (uint8_t)(mean_of8x8_blocks[34] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_35] =
        (uint8_t)(mean_of8x8_blocks[35] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_36] =
        (uint8_t)(mean_of8x8_blocks[36] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_37] =
        (uint8_t)(mean_of8x8_blocks[37] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_38] =
        (uint8_t)(mean_of8x8_blocks[38] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_39] =
        (uint8_t)(mean_of8x8_blocks[39] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_40] =
        (uint8_t)(mean_of8x8_blocks[40] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_41] =
        (uint8_t)(mean_of8x8_blocks[41] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_42] =
        (uint8_t)(mean_of8x8_blocks[42] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_43] =
        (uint8_t)(mean_of8x8_blocks[43] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_44] =
        (uint8_t)(mean_of8x8_blocks[44] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_45] =
        (uint8_t)(mean_of8x8_blocks[45] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_46] =
        (uint8_t)(mean_of8x8_blocks[46] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_47] =
        (uint8_t)(mean_of8x8_blocks[47] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_48] =
        (uint8_t)(mean_of8x8_blocks[48] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_49] =
        (uint8_t)(mean_of8x8_blocks[49] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_50] =
        (uint8_t)(mean_of8x8_blocks[50] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_51] =
        (uint8_t)(mean_of8x8_blocks[51] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_52] =
        (uint8_t)(mean_of8x8_blocks[52] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_53] =
        (uint8_t)(mean_of8x8_blocks[53] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_54] =
        (uint8_t)(mean_of8x8_blocks[54] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_55] =
        (uint8_t)(mean_of8x8_blocks[55] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_56] =
        (uint8_t)(mean_of8x8_blocks[56] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_57] =
        (uint8_t)(mean_of8x8_blocks[57] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_58] =
        (uint8_t)(mean_of8x8_blocks[58] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_59] =
        (uint8_t)(mean_of8x8_blocks[59] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_60] =
        (uint8_t)(mean_of8x8_blocks[60] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_61] =
        (uint8_t)(mean_of8x8_blocks[61] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_62] =
        (uint8_t)(mean_of8x8_blocks[62] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_8x8_63] =
        (uint8_t)(mean_of8x8_blocks[63] >> MEAN_PRECISION);

    // 16x16 mean
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_16x16_0] =
        (uint8_t)(mean_of_16x16_blocks[0] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_16x16_1] =
        (uint8_t)(mean_of_16x16_blocks[1] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_16x16_2] =
        (uint8_t)(mean_of_16x16_blocks[2] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_16x16_3] =
        (uint8_t)(mean_of_16x16_blocks[3] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_16x16_4] =
        (uint8_t)(mean_of_16x16_blocks[4] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_16x16_5] =
        (uint8_t)(mean_of_16x16_blocks[5] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_16x16_6] =
        (uint8_t)(mean_of_16x16_blocks[6] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_16x16_7] =
        (uint8_t)(mean_of_16x16_blocks[7] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_16x16_8] =
        (uint8_t)(mean_of_16x16_blocks[8] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_16x16_9] =
        (uint8_t)(mean_of_16x16_blocks[9] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_16x16_10] =
        (uint8_t)(mean_of_16x16_blocks[10] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_16x16_11] =
        (uint8_t)(mean_of_16x16_blocks[11] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_16x16_12] =
        (uint8_t)(mean_of_16x16_blocks[12] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_16x16_13] =
        (uint8_t)(mean_of_16x16_blocks[13] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_16x16_14] =
        (uint8_t)(mean_of_16x16_blocks[14] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_16x16_15] =
        (uint8_t)(mean_of_16x16_blocks[15] >> MEAN_PRECISION);

    // 32x32 mean
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_32x32_0] =
        (uint8_t)(mean_of_32x32_blocks[0] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_32x32_1] =
        (uint8_t)(mean_of_32x32_blocks[1] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_32x32_2] =
        (uint8_t)(mean_of_32x32_blocks[2] >> MEAN_PRECISION);
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_32x32_3] =
        (uint8_t)(mean_of_32x32_blocks[3] >> MEAN_PRECISION);

    // 64x64 mean
    pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_64x64] =
        (uint8_t)(mean_of_64x64_blocks >> MEAN_PRECISION);

    // 8x8 variances
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_0] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[0] - (mean_of8x8_blocks[0] * mean_of8x8_blocks[0])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_1] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[1] - (mean_of8x8_blocks[1] * mean_of8x8_blocks[1])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_2] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[2] - (mean_of8x8_blocks[2] * mean_of8x8_blocks[2])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_3] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[3] - (mean_of8x8_blocks[3] * mean_of8x8_blocks[3])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_4] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[4] - (mean_of8x8_blocks[4] * mean_of8x8_blocks[4])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_5] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[5] - (mean_of8x8_blocks[5] * mean_of8x8_blocks[5])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_6] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[6] - (mean_of8x8_blocks[6] * mean_of8x8_blocks[6])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_7] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[7] - (mean_of8x8_blocks[7] * mean_of8x8_blocks[7])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_8] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[8] - (mean_of8x8_blocks[8] * mean_of8x8_blocks[8])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_9] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[9] - (mean_of8x8_blocks[9] * mean_of8x8_blocks[9])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_10] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[10] - (mean_of8x8_blocks[10] * mean_of8x8_blocks[10])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_11] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[11] - (mean_of8x8_blocks[11] * mean_of8x8_blocks[11])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_12] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[12] - (mean_of8x8_blocks[12] * mean_of8x8_blocks[12])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_13] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[13] - (mean_of8x8_blocks[13] * mean_of8x8_blocks[13])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_14] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[14] - (mean_of8x8_blocks[14] * mean_of8x8_blocks[14])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_15] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[15] - (mean_of8x8_blocks[15] * mean_of8x8_blocks[15])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_16] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[16] - (mean_of8x8_blocks[16] * mean_of8x8_blocks[16])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_17] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[17] - (mean_of8x8_blocks[17] * mean_of8x8_blocks[17])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_18] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[18] - (mean_of8x8_blocks[18] * mean_of8x8_blocks[18])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_19] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[19] - (mean_of8x8_blocks[19] * mean_of8x8_blocks[19])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_20] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[20] - (mean_of8x8_blocks[20] * mean_of8x8_blocks[20])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_21] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[21] - (mean_of8x8_blocks[21] * mean_of8x8_blocks[21])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_22] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[22] - (mean_of8x8_blocks[22] * mean_of8x8_blocks[22])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_23] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[23] - (mean_of8x8_blocks[23] * mean_of8x8_blocks[23])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_24] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[24] - (mean_of8x8_blocks[24] * mean_of8x8_blocks[24])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_25] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[25] - (mean_of8x8_blocks[25] * mean_of8x8_blocks[25])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_26] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[26] - (mean_of8x8_blocks[26] * mean_of8x8_blocks[26])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_27] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[27] - (mean_of8x8_blocks[27] * mean_of8x8_blocks[27])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_28] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[28] - (mean_of8x8_blocks[28] * mean_of8x8_blocks[28])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_29] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[29] - (mean_of8x8_blocks[29] * mean_of8x8_blocks[29])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_30] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[30] - (mean_of8x8_blocks[30] * mean_of8x8_blocks[30])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_31] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[31] - (mean_of8x8_blocks[31] * mean_of8x8_blocks[31])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_32] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[32] - (mean_of8x8_blocks[32] * mean_of8x8_blocks[32])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_33] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[33] - (mean_of8x8_blocks[33] * mean_of8x8_blocks[33])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_34] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[34] - (mean_of8x8_blocks[34] * mean_of8x8_blocks[34])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_35] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[35] - (mean_of8x8_blocks[35] * mean_of8x8_blocks[35])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_36] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[36] - (mean_of8x8_blocks[36] * mean_of8x8_blocks[36])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_37] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[37] - (mean_of8x8_blocks[37] * mean_of8x8_blocks[37])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_38] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[38] - (mean_of8x8_blocks[38] * mean_of8x8_blocks[38])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_39] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[39] - (mean_of8x8_blocks[39] * mean_of8x8_blocks[39])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_40] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[40] - (mean_of8x8_blocks[40] * mean_of8x8_blocks[40])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_41] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[41] - (mean_of8x8_blocks[41] * mean_of8x8_blocks[41])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_42] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[42] - (mean_of8x8_blocks[42] * mean_of8x8_blocks[42])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_43] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[43] - (mean_of8x8_blocks[43] * mean_of8x8_blocks[43])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_44] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[44] - (mean_of8x8_blocks[44] * mean_of8x8_blocks[44])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_45] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[45] - (mean_of8x8_blocks[45] * mean_of8x8_blocks[45])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_46] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[46] - (mean_of8x8_blocks[46] * mean_of8x8_blocks[46])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_47] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[47] - (mean_of8x8_blocks[47] * mean_of8x8_blocks[47])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_48] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[48] - (mean_of8x8_blocks[48] * mean_of8x8_blocks[48])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_49] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[49] - (mean_of8x8_blocks[49] * mean_of8x8_blocks[49])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_50] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[50] - (mean_of8x8_blocks[50] * mean_of8x8_blocks[50])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_51] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[51] - (mean_of8x8_blocks[51] * mean_of8x8_blocks[51])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_52] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[52] - (mean_of8x8_blocks[52] * mean_of8x8_blocks[52])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_53] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[53] - (mean_of8x8_blocks[53] * mean_of8x8_blocks[53])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_54] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[54] - (mean_of8x8_blocks[54] * mean_of8x8_blocks[54])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_55] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[55] - (mean_of8x8_blocks[55] * mean_of8x8_blocks[55])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_56] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[56] - (mean_of8x8_blocks[56] * mean_of8x8_blocks[56])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_57] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[57] - (mean_of8x8_blocks[57] * mean_of8x8_blocks[57])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_58] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[58] - (mean_of8x8_blocks[58] * mean_of8x8_blocks[58])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_59] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[59] - (mean_of8x8_blocks[59] * mean_of8x8_blocks[59])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_60] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[60] - (mean_of8x8_blocks[60] * mean_of8x8_blocks[60])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_61] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[61] - (mean_of8x8_blocks[61] * mean_of8x8_blocks[61])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_62] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[62] - (mean_of8x8_blocks[62] * mean_of8x8_blocks[62])) >>
        VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_8x8_63] = (uint16_t)(
        (mean_of_8x8_squared_values_blocks[63] - (mean_of8x8_blocks[63] * mean_of8x8_blocks[63])) >>
        VARIANCE_PRECISION);

    // 16x16 variances
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_0] =
        (uint16_t)((mean_of16x16_squared_values_blocks[0] -
                    (mean_of_16x16_blocks[0] * mean_of_16x16_blocks[0])) >>
                   VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_1] =
        (uint16_t)((mean_of16x16_squared_values_blocks[1] -
                    (mean_of_16x16_blocks[1] * mean_of_16x16_blocks[1])) >>
                   VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_2] =
        (uint16_t)((mean_of16x16_squared_values_blocks[2] -
                    (mean_of_16x16_blocks[2] * mean_of_16x16_blocks[2])) >>
                   VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_3] =
        (uint16_t)((mean_of16x16_squared_values_blocks[3] -
                    (mean_of_16x16_blocks[3] * mean_of_16x16_blocks[3])) >>
                   VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_4] =
        (uint16_t)((mean_of16x16_squared_values_blocks[4] -
                    (mean_of_16x16_blocks[4] * mean_of_16x16_blocks[4])) >>
                   VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_5] =
        (uint16_t)((mean_of16x16_squared_values_blocks[5] -
                    (mean_of_16x16_blocks[5] * mean_of_16x16_blocks[5])) >>
                   VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_6] =
        (uint16_t)((mean_of16x16_squared_values_blocks[6] -
                    (mean_of_16x16_blocks[6] * mean_of_16x16_blocks[6])) >>
                   VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_7] =
        (uint16_t)((mean_of16x16_squared_values_blocks[7] -
                    (mean_of_16x16_blocks[7] * mean_of_16x16_blocks[7])) >>
                   VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_8] =
        (uint16_t)((mean_of16x16_squared_values_blocks[8] -
                    (mean_of_16x16_blocks[8] * mean_of_16x16_blocks[8])) >>
                   VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_9] =
        (uint16_t)((mean_of16x16_squared_values_blocks[9] -
                    (mean_of_16x16_blocks[9] * mean_of_16x16_blocks[9])) >>
                   VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_10] =
        (uint16_t)((mean_of16x16_squared_values_blocks[10] -
                    (mean_of_16x16_blocks[10] * mean_of_16x16_blocks[10])) >>
                   VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_11] =
        (uint16_t)((mean_of16x16_squared_values_blocks[11] -
                    (mean_of_16x16_blocks[11] * mean_of_16x16_blocks[11])) >>
                   VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_12] =
        (uint16_t)((mean_of16x16_squared_values_blocks[12] -
                    (mean_of_16x16_blocks[12] * mean_of_16x16_blocks[12])) >>
                   VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_13] =
        (uint16_t)((mean_of16x16_squared_values_blocks[13] -
                    (mean_of_16x16_blocks[13] * mean_of_16x16_blocks[13])) >>
                   VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_14] =
        (uint16_t)((mean_of16x16_squared_values_blocks[14] -
                    (mean_of_16x16_blocks[14] * mean_of_16x16_blocks[14])) >>
                   VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_16x16_15] =
        (uint16_t)((mean_of16x16_squared_values_blocks[15] -
                    (mean_of_16x16_blocks[15] * mean_of_16x16_blocks[15])) >>
                   VARIANCE_PRECISION);

    // 32x32 variances
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_32x32_0] =
        (uint16_t)((mean_of32x32_squared_values_blocks[0] -
                    (mean_of_32x32_blocks[0] * mean_of_32x32_blocks[0])) >>
                   VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_32x32_1] =
        (uint16_t)((mean_of32x32_squared_values_blocks[1] -
                    (mean_of_32x32_blocks[1] * mean_of_32x32_blocks[1])) >>
                   VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_32x32_2] =
        (uint16_t)((mean_of32x32_squared_values_blocks[2] -
                    (mean_of_32x32_blocks[2] * mean_of_32x32_blocks[2])) >>
                   VARIANCE_PRECISION);
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_32x32_3] =
        (uint16_t)((mean_of32x32_squared_values_blocks[3] -
                    (mean_of_32x32_blocks[3] * mean_of_32x32_blocks[3])) >>
                   VARIANCE_PRECISION);

    // 64x64 variance
    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_64x64] = (uint16_t)(
        (mean_of64x64_squared_values_blocks - (mean_of_64x64_blocks * mean_of_64x64_blocks)) >>
        VARIANCE_PRECISION);

    return return_error;
}

static int32_t apply_denoise_2d(SequenceControlSet *scs_ptr, PictureParentControlSet *pcs_ptr,
                                EbPictureBufferDesc *inputPicturePointer) {
    if (eb_aom_denoise_and_model_run(pcs_ptr->denoise_and_model,
                                     inputPicturePointer,
                                     &pcs_ptr->frm_hdr.film_grain_params,
                                     scs_ptr->static_config.encoder_bit_depth > EB_8BIT)) {}
    return 0;
}

EbErrorType denoise_estimate_film_grain(SequenceControlSet *     scs_ptr,
                                        PictureParentControlSet *pcs_ptr) {
    EbErrorType return_error = EB_ErrorNone;

    FrameHeader *frm_hdr = &pcs_ptr->frm_hdr;

    EbPictureBufferDesc *input_picture_ptr = pcs_ptr->enhanced_picture_ptr;
    frm_hdr->film_grain_params.apply_grain = 0;

    if (scs_ptr->film_grain_denoise_strength) {
        if (apply_denoise_2d(scs_ptr, pcs_ptr, input_picture_ptr) < 0) return 1;
    }

    scs_ptr->seq_header.film_grain_params_present |= frm_hdr->film_grain_params.apply_grain;

    return return_error; //todo: add proper error handling
}

/************************************************
 * Set Picture Parameters based on input configuration
 ** Setting Number of regions per resolution
 ** Setting width and height for subpicture and when picture scan type is 1
 ************************************************/
void set_picture_parameters_for_statistics_gathering(SequenceControlSet *scs_ptr) {
    scs_ptr->picture_analysis_number_of_regions_per_width =
        HIGHER_THAN_CLASS_1_REGION_SPLIT_PER_WIDTH;
    scs_ptr->picture_analysis_number_of_regions_per_height =
        HIGHER_THAN_CLASS_1_REGION_SPLIT_PER_HEIGHT;

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
void picture_pre_processing_operations(PictureParentControlSet *pcs_ptr,
                                       SequenceControlSet *scs_ptr, uint32_t sb_total_count) {
    if (scs_ptr->film_grain_denoise_strength) {
        denoise_estimate_film_grain(scs_ptr, pcs_ptr);
    } else {
        //Reset the flat noise flag array to False for both RealTime/HighComplexity Modes
        for (uint32_t sb_coding_order = 0; sb_coding_order < sb_total_count; ++sb_coding_order)
            pcs_ptr->sb_flat_noise_array[sb_coding_order] = 0;
        pcs_ptr->pic_noise_class =
            PIC_NOISE_CLASS_INV; //this init is for both REAL-TIME and BEST-QUALITY
    }
    return;
}

/**************************************************************
* Generate picture histogram bins for YUV pixel intensity *
* Calculation is done on a region based (Set previously, resolution dependent)
**************************************************************/
void sub_sample_luma_generate_pixel_intensity_histogram_bins(
    SequenceControlSet *scs_ptr, PictureParentControlSet *pcs_ptr,
    EbPictureBufferDesc *input_picture_ptr, uint64_t *sum_avg_intensity_ttl_regions_luma) {
    uint32_t region_width;
    uint32_t region_height;
    uint32_t region_width_offset;
    uint32_t region_height_offset;
    uint32_t region_in_picture_width_index;
    uint32_t region_in_picture_height_index;
    uint32_t histogram_bin;
    uint64_t sum;

    region_width = input_picture_ptr->width / scs_ptr->picture_analysis_number_of_regions_per_width;
    region_height =
        input_picture_ptr->height / scs_ptr->picture_analysis_number_of_regions_per_height;

    // Loop over regions inside the picture
    for (region_in_picture_width_index = 0;
         region_in_picture_width_index < scs_ptr->picture_analysis_number_of_regions_per_width;
         region_in_picture_width_index++) { // loop over horizontal regions
        for (region_in_picture_height_index = 0;
             region_in_picture_height_index <
             scs_ptr->picture_analysis_number_of_regions_per_height;
             region_in_picture_height_index++) { // loop over vertical regions

            // Initialize bins to 1
            initialize_buffer_32bits(pcs_ptr->picture_histogram[region_in_picture_width_index]
                                                               [region_in_picture_height_index][0],
                                     64,
                                     0,
                                     1);

            region_width_offset =
                (region_in_picture_width_index ==
                 scs_ptr->picture_analysis_number_of_regions_per_width - 1)
                    ? input_picture_ptr->width -
                          (scs_ptr->picture_analysis_number_of_regions_per_width * region_width)
                    : 0;

            region_height_offset =
                (region_in_picture_height_index ==
                 scs_ptr->picture_analysis_number_of_regions_per_height - 1)
                    ? input_picture_ptr->height -
                          (scs_ptr->picture_analysis_number_of_regions_per_height * region_height)
                    : 0;

            // Y Histogram
            calculate_histogram(
                &input_picture_ptr->buffer_y[(input_picture_ptr->origin_x +
                                              region_in_picture_width_index * region_width) +
                                             ((input_picture_ptr->origin_y +
                                               region_in_picture_height_index * region_height) *
                                              input_picture_ptr->stride_y)],
                region_width + region_width_offset,
                region_height + region_height_offset,
                input_picture_ptr->stride_y,
                1,
                pcs_ptr->picture_histogram[region_in_picture_width_index]
                                          [region_in_picture_height_index][0],
                &sum);

            pcs_ptr->average_intensity_per_region[region_in_picture_width_index]
                                                 [region_in_picture_height_index][0] = (uint8_t)(
                (sum +
                 (((region_width + region_width_offset) * (region_height + region_height_offset)) >>
                  1)) /
                ((region_width + region_width_offset) * (region_height + region_height_offset)));
            (*sum_avg_intensity_ttl_regions_luma) += (sum << 4);
            for (histogram_bin = 0; histogram_bin < HISTOGRAM_NUMBER_OF_BINS;
                 histogram_bin++) { // Loop over the histogram bins
                pcs_ptr->picture_histogram[region_in_picture_width_index]
                                          [region_in_picture_height_index][0][histogram_bin] =
                    pcs_ptr->picture_histogram[region_in_picture_width_index]
                                              [region_in_picture_height_index][0][histogram_bin]
                    << 4;
            }
        }
    }

    return;
}

void sub_sample_chroma_generate_pixel_intensity_histogram_bins(
    SequenceControlSet *scs_ptr, PictureParentControlSet *pcs_ptr,
    EbPictureBufferDesc *input_picture_ptr, uint64_t *sum_avg_intensity_ttl_regions_cb,
    uint64_t *sum_avg_intensity_ttl_regions_cr) {
    uint64_t sum;
    uint32_t region_width;
    uint32_t region_height;
    uint32_t region_width_offset;
    uint32_t region_height_offset;
    uint32_t region_in_picture_width_index;
    uint32_t region_in_picture_height_index;

    uint16_t histogram_bin;
    uint8_t  decim_step = 4;

    region_width = input_picture_ptr->width / scs_ptr->picture_analysis_number_of_regions_per_width;
    region_height =
        input_picture_ptr->height / scs_ptr->picture_analysis_number_of_regions_per_height;

    // Loop over regions inside the picture
    for (region_in_picture_width_index = 0;
         region_in_picture_width_index < scs_ptr->picture_analysis_number_of_regions_per_width;
         region_in_picture_width_index++) { // loop over horizontal regions
        for (region_in_picture_height_index = 0;
             region_in_picture_height_index <
             scs_ptr->picture_analysis_number_of_regions_per_height;
             region_in_picture_height_index++) { // loop over vertical regions

            // Initialize bins to 1
            initialize_buffer_32bits(pcs_ptr->picture_histogram[region_in_picture_width_index]
                                                               [region_in_picture_height_index][1],
                                     64,
                                     0,
                                     1);
            initialize_buffer_32bits(pcs_ptr->picture_histogram[region_in_picture_width_index]
                                                               [region_in_picture_height_index][2],
                                     64,
                                     0,
                                     1);

            region_width_offset =
                (region_in_picture_width_index ==
                 scs_ptr->picture_analysis_number_of_regions_per_width - 1)
                    ? input_picture_ptr->width -
                          (scs_ptr->picture_analysis_number_of_regions_per_width * region_width)
                    : 0;

            region_height_offset =
                (region_in_picture_height_index ==
                 scs_ptr->picture_analysis_number_of_regions_per_height - 1)
                    ? input_picture_ptr->height -
                          (scs_ptr->picture_analysis_number_of_regions_per_height * region_height)
                    : 0;

            // U Histogram
            calculate_histogram(
                &input_picture_ptr->buffer_cb[((input_picture_ptr->origin_x +
                                                region_in_picture_width_index * region_width) >>
                                               1) +
                                              (((input_picture_ptr->origin_y +
                                                 region_in_picture_height_index * region_height) >>
                                                1) *
                                               input_picture_ptr->stride_cb)],
                (region_width + region_width_offset) >> 1,
                (region_height + region_height_offset) >> 1,
                input_picture_ptr->stride_cb,
                decim_step,
                pcs_ptr->picture_histogram[region_in_picture_width_index]
                                          [region_in_picture_height_index][1],
                &sum);

            sum = (sum << decim_step);
            *sum_avg_intensity_ttl_regions_cb += sum;
            pcs_ptr->average_intensity_per_region[region_in_picture_width_index]
                                                 [region_in_picture_height_index][1] = (uint8_t)(
                (sum +
                 (((region_width + region_width_offset) * (region_height + region_height_offset)) >>
                  3)) /
                (((region_width + region_width_offset) * (region_height + region_height_offset)) >>
                 2));

            for (histogram_bin = 0; histogram_bin < HISTOGRAM_NUMBER_OF_BINS;
                 histogram_bin++) { // Loop over the histogram bins
                pcs_ptr->picture_histogram[region_in_picture_width_index]
                                          [region_in_picture_height_index][1][histogram_bin] =
                    pcs_ptr->picture_histogram[region_in_picture_width_index]
                                              [region_in_picture_height_index][1][histogram_bin]
                    << decim_step;
            }

            // V Histogram
            calculate_histogram(
                &input_picture_ptr->buffer_cr[((input_picture_ptr->origin_x +
                                                region_in_picture_width_index * region_width) >>
                                               1) +
                                              (((input_picture_ptr->origin_y +
                                                 region_in_picture_height_index * region_height) >>
                                                1) *
                                               input_picture_ptr->stride_cr)],
                (region_width + region_width_offset) >> 1,
                (region_height + region_height_offset) >> 1,
                input_picture_ptr->stride_cr,
                decim_step,
                pcs_ptr->picture_histogram[region_in_picture_width_index]
                                          [region_in_picture_height_index][2],
                &sum);

            sum = (sum << decim_step);
            *sum_avg_intensity_ttl_regions_cr += sum;
            pcs_ptr->average_intensity_per_region[region_in_picture_width_index]
                                                 [region_in_picture_height_index][2] = (uint8_t)(
                (sum +
                 (((region_width + region_width_offset) * (region_height + region_height_offset)) >>
                  3)) /
                (((region_width + region_width_offset) * (region_height + region_height_offset)) >>
                 2));

            for (histogram_bin = 0; histogram_bin < HISTOGRAM_NUMBER_OF_BINS;
                 histogram_bin++) { // Loop over the histogram bins
                pcs_ptr->picture_histogram[region_in_picture_width_index]
                                          [region_in_picture_height_index][2][histogram_bin] =
                    pcs_ptr->picture_histogram[region_in_picture_width_index]
                                              [region_in_picture_height_index][2][histogram_bin]
                    << decim_step;
            }
        }
    }
    return;
}
/************************************************
 * compute_picture_spatial_statistics
 ** Compute Block Variance
 ** Compute Picture Variance
 ** Compute Block Mean for all blocks in the picture
 ************************************************/
void compute_picture_spatial_statistics(SequenceControlSet *     scs_ptr,
                                        PictureParentControlSet *pcs_ptr,
                                        EbPictureBufferDesc *    input_picture_ptr,
                                        EbPictureBufferDesc *    input_padded_picture_ptr,
                                        uint32_t                 sb_total_count) {
    uint32_t sb_index;
    uint32_t sb_origin_x; // to avoid using child PCS
    uint32_t sb_origin_y;
    uint32_t input_luma_origin_index;
    uint32_t input_cb_origin_index;
    uint32_t input_cr_origin_index;
    uint64_t pic_tot_variance;

    // Variance
    pic_tot_variance = 0;

    for (sb_index = 0; sb_index < pcs_ptr->sb_total_count; ++sb_index) {
        SbParams *sb_params = &pcs_ptr->sb_params_array[sb_index];

        sb_origin_x             = sb_params->origin_x;
        sb_origin_y             = sb_params->origin_y;
        input_luma_origin_index = (input_padded_picture_ptr->origin_y + sb_origin_y) *
                                      input_padded_picture_ptr->stride_y +
                                  input_padded_picture_ptr->origin_x + sb_origin_x;

        input_cb_origin_index =
            ((input_picture_ptr->origin_y + sb_origin_y) >> 1) * input_picture_ptr->stride_cb +
            ((input_picture_ptr->origin_x + sb_origin_x) >> 1);
        input_cr_origin_index =
            ((input_picture_ptr->origin_y + sb_origin_y) >> 1) * input_picture_ptr->stride_cr +
            ((input_picture_ptr->origin_x + sb_origin_x) >> 1);

        compute_block_mean_compute_variance(
            scs_ptr, pcs_ptr, input_padded_picture_ptr, sb_index, input_luma_origin_index);

        if (sb_params->is_complete_sb) {
            compute_chroma_block_mean(scs_ptr,
                                      pcs_ptr,
                                      input_picture_ptr,
                                      sb_index,
                                      input_cb_origin_index,
                                      input_cr_origin_index);
        } else {
            zero_out_chroma_block_mean(pcs_ptr, sb_index);
        }

        pic_tot_variance += (pcs_ptr->variance[sb_index][RASTER_SCAN_CU_INDEX_64x64]);
    }

    pcs_ptr->pic_avg_variance = (uint16_t)(pic_tot_variance / sb_total_count);

    return;
}

void calculate_input_average_intensity(SequenceControlSet *     scs_ptr,
                                       PictureParentControlSet *pcs_ptr,
                                       EbPictureBufferDesc *    input_picture_ptr,
                                       uint64_t                 sum_avg_intensity_ttl_regions_luma,
                                       uint64_t                 sum_avg_intensity_ttl_regions_cb,
                                       uint64_t                 sum_avg_intensity_ttl_regions_cr) {
    if (scs_ptr->scd_mode == SCD_MODE_0) {
        uint16_t block_index_in_width;
        uint16_t block_index_in_height;
        uint64_t mean = 0;

        const uint16_t stride_y = input_picture_ptr->stride_y;
        // Loop over 8x8 blocks and calculates the mean value
        if (scs_ptr->block_mean_calc_prec == BLOCK_MEAN_PREC_FULL) {
            for (block_index_in_height = 0; input_picture_ptr->height >> 3 > block_index_in_height;
                 ++block_index_in_height) {
                for (block_index_in_width = 0; input_picture_ptr->width >> 3 > block_index_in_width;
                     ++block_index_in_width)
                    mean += compute_mean_8x8(
                        &(input_picture_ptr->buffer_y[(block_index_in_width << 3) +
                                                      (block_index_in_height << 3) *
                                                          input_picture_ptr->stride_y]),
                        input_picture_ptr->stride_y,
                        8,
                        8);
            }
        } else {
            for (block_index_in_height = 0; input_picture_ptr->height >> 3 > block_index_in_height;
                 ++block_index_in_height) {
                for (block_index_in_width = 0; input_picture_ptr->width >> 3 > block_index_in_width;
                     ++block_index_in_width)
                    mean += compute_sub_mean_8x8(
                        &(input_picture_ptr->buffer_y[(block_index_in_width << 3) +
                                                      (block_index_in_height << 3) * stride_y]),
                        stride_y);
            }
        }
        mean = ((mean + ((input_picture_ptr->height * input_picture_ptr->width) >> 7)) /
                ((input_picture_ptr->height * input_picture_ptr->width) >> 6));
        mean = (mean + (1 << (MEAN_PRECISION - 1))) >> MEAN_PRECISION;
        pcs_ptr->average_intensity[0] = (uint8_t)mean;
    }

    else {
        pcs_ptr->average_intensity[0] =
            (uint8_t)((sum_avg_intensity_ttl_regions_luma +
                       ((input_picture_ptr->width * input_picture_ptr->height) >> 1)) /
                      (input_picture_ptr->width * input_picture_ptr->height));
        pcs_ptr->average_intensity[1] =
            (uint8_t)((sum_avg_intensity_ttl_regions_cb +
                       ((input_picture_ptr->width * input_picture_ptr->height) >> 3)) /
                      ((input_picture_ptr->width * input_picture_ptr->height) >> 2));
        pcs_ptr->average_intensity[2] =
            (uint8_t)((sum_avg_intensity_ttl_regions_cr +
                       ((input_picture_ptr->width * input_picture_ptr->height) >> 3)) /
                      ((input_picture_ptr->width * input_picture_ptr->height) >> 2));
    }

    return;
}

/************************************************
 * Gathering statistics per picture
 ** Calculating the pixel intensity histogram bins per picture needed for SCD
 ** Computing Picture Variance
 ************************************************/
void gathering_picture_statistics(SequenceControlSet *scs_ptr, PictureParentControlSet *pcs_ptr,
                                  EbPictureBufferDesc *input_picture_ptr,
                                  EbPictureBufferDesc *input_padded_picture_ptr,
                                  EbPictureBufferDesc *sixteenth_decimated_picture_ptr,
                                  uint32_t             sb_total_count) {
    uint64_t sum_avg_intensity_ttl_regions_luma = 0;
    uint64_t sum_avg_intensity_ttl_regions_cb   = 0;
    uint64_t sum_avg_intensity_ttl_regions_cr   = 0;

    // Histogram bins
    // Use 1/16 Luma for Histogram generation
    // 1/16 input ready
    sub_sample_luma_generate_pixel_intensity_histogram_bins(
        scs_ptr, pcs_ptr, sixteenth_decimated_picture_ptr, &sum_avg_intensity_ttl_regions_luma);

    // Use 1/4 Chroma for Histogram generation
    // 1/4 input not ready => perform operation on the fly
    sub_sample_chroma_generate_pixel_intensity_histogram_bins(scs_ptr,
                                                              pcs_ptr,
                                                              input_picture_ptr,
                                                              &sum_avg_intensity_ttl_regions_cb,
                                                              &sum_avg_intensity_ttl_regions_cr);
    //
    // Calculate the LUMA average intensity
    calculate_input_average_intensity(scs_ptr,
                                      pcs_ptr,
                                      input_picture_ptr,
                                      sum_avg_intensity_ttl_regions_luma,
                                      sum_avg_intensity_ttl_regions_cb,
                                      sum_avg_intensity_ttl_regions_cr);

    compute_picture_spatial_statistics(
        scs_ptr, pcs_ptr, input_picture_ptr, input_padded_picture_ptr, sb_total_count);

    return;
}

/************************************************
 * Pad Picture at the right and bottom sides
 ** To match a multiple of min CU size in width and height
 ************************************************/
void pad_picture_to_multiple_of_min_blk_size_dimensions(SequenceControlSet * scs_ptr,
                                                        EbPictureBufferDesc *input_picture_ptr) {
    EbBool is16_bit_input = (EbBool)(scs_ptr->static_config.encoder_bit_depth > EB_8BIT);

    uint32_t       color_format  = input_picture_ptr->color_format;
    const uint16_t subsampling_x = (color_format == EB_YUV444 ? 1 : 2) - 1;
    const uint16_t subsampling_y = (color_format >= EB_YUV422 ? 1 : 2) - 1;

    // Input Picture Padding
    pad_input_picture(
        &input_picture_ptr->buffer_y[input_picture_ptr->origin_x +
                                     (input_picture_ptr->origin_y * input_picture_ptr->stride_y)],
        input_picture_ptr->stride_y,
        (input_picture_ptr->width - scs_ptr->pad_right),
        (input_picture_ptr->height - scs_ptr->pad_bottom),
        scs_ptr->pad_right,
        scs_ptr->pad_bottom);

    pad_input_picture(
        &input_picture_ptr->buffer_cb[(input_picture_ptr->origin_x >> subsampling_x) +
                                      ((input_picture_ptr->origin_y >> subsampling_y) *
                                       input_picture_ptr->stride_cb)],
        input_picture_ptr->stride_cb,
        (input_picture_ptr->width - scs_ptr->pad_right) >> subsampling_x,
        (input_picture_ptr->height - scs_ptr->pad_bottom) >> subsampling_y,
        scs_ptr->pad_right >> subsampling_x,
        scs_ptr->pad_bottom >> subsampling_y);

    pad_input_picture(
        &input_picture_ptr->buffer_cr[(input_picture_ptr->origin_x >> subsampling_x) +
                                      ((input_picture_ptr->origin_y >> subsampling_y) *
                                       input_picture_ptr->stride_cb)],
        input_picture_ptr->stride_cr,
        (input_picture_ptr->width - scs_ptr->pad_right) >> subsampling_x,
        (input_picture_ptr->height - scs_ptr->pad_bottom) >> subsampling_y,
        scs_ptr->pad_right >> subsampling_x,
        scs_ptr->pad_bottom >> subsampling_y);

    if (is16_bit_input) {
        pad_input_picture(
            &input_picture_ptr->buffer_bit_inc_y[input_picture_ptr->origin_x +
                                                 (input_picture_ptr->origin_y *
                                                  input_picture_ptr->stride_bit_inc_y)],
            input_picture_ptr->stride_bit_inc_y,
            (input_picture_ptr->width - scs_ptr->pad_right),
            (input_picture_ptr->height - scs_ptr->pad_bottom),
            scs_ptr->pad_right,
            scs_ptr->pad_bottom);

        pad_input_picture(
            &input_picture_ptr->buffer_bit_inc_cb[(input_picture_ptr->origin_x >> subsampling_x) +
                                                  ((input_picture_ptr->origin_y >> subsampling_y) *
                                                   input_picture_ptr->stride_bit_inc_cb)],
            input_picture_ptr->stride_bit_inc_cb,
            (input_picture_ptr->width - scs_ptr->pad_right) >> subsampling_x,
            (input_picture_ptr->height - scs_ptr->pad_bottom) >> subsampling_y,
            scs_ptr->pad_right >> subsampling_x,
            scs_ptr->pad_bottom >> subsampling_y);

        pad_input_picture(
            &input_picture_ptr->buffer_bit_inc_cr[(input_picture_ptr->origin_x >> subsampling_x) +
                                                  ((input_picture_ptr->origin_y >> subsampling_y) *
                                                   input_picture_ptr->stride_bit_inc_cb)],
            input_picture_ptr->stride_bit_inc_cr,
            (input_picture_ptr->width - scs_ptr->pad_right) >> subsampling_x,
            (input_picture_ptr->height - scs_ptr->pad_bottom) >> subsampling_y,
            scs_ptr->pad_right >> subsampling_x,
            scs_ptr->pad_bottom >> subsampling_y);
    }

    return;
}

/************************************************
 * Pad Picture at the right and bottom sides
 ** To complete border SB smaller than SB size
 ************************************************/
void pad_picture_to_multiple_of_sb_dimensions(EbPictureBufferDesc *input_padded_picture_ptr) {
    // Generate Padding
    generate_padding(&input_padded_picture_ptr->buffer_y[0],
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
void downsample_decimation_input_picture(PictureParentControlSet *pcs_ptr,
                                         EbPictureBufferDesc *    input_padded_picture_ptr,
                                         EbPictureBufferDesc *    quarter_decimated_picture_ptr,
                                         EbPictureBufferDesc *    sixteenth_decimated_picture_ptr) {
    // Decimate input picture for HME L0 and L1
    if (pcs_ptr->enable_hme_flag || pcs_ptr->tf_enable_hme_flag) {
        if (pcs_ptr->enable_hme_level1_flag || pcs_ptr->tf_enable_hme_level1_flag) {
            decimation_2d(
                &input_padded_picture_ptr->buffer_y[input_padded_picture_ptr->origin_x +
                                                    input_padded_picture_ptr->origin_y *
                                                        input_padded_picture_ptr->stride_y],
                input_padded_picture_ptr->stride_y,
                input_padded_picture_ptr->width,
                input_padded_picture_ptr->height,
                &quarter_decimated_picture_ptr
                     ->buffer_y[quarter_decimated_picture_ptr->origin_x +
                                quarter_decimated_picture_ptr->origin_x *
                                    quarter_decimated_picture_ptr->stride_y],
                quarter_decimated_picture_ptr->stride_y,
                2);
            generate_padding(&quarter_decimated_picture_ptr->buffer_y[0],
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
        &input_padded_picture_ptr
             ->buffer_y[input_padded_picture_ptr->origin_x +
                        input_padded_picture_ptr->origin_y * input_padded_picture_ptr->stride_y],
        input_padded_picture_ptr->stride_y,
        input_padded_picture_ptr->width,
        input_padded_picture_ptr->height,
        &sixteenth_decimated_picture_ptr->buffer_y[sixteenth_decimated_picture_ptr->origin_x +
                                                   sixteenth_decimated_picture_ptr->origin_x *
                                                       sixteenth_decimated_picture_ptr->stride_y],
        sixteenth_decimated_picture_ptr->stride_y,
        4);

    generate_padding(&sixteenth_decimated_picture_ptr->buffer_y[0],
                     sixteenth_decimated_picture_ptr->stride_y,
                     sixteenth_decimated_picture_ptr->width,
                     sixteenth_decimated_picture_ptr->height,
                     sixteenth_decimated_picture_ptr->origin_x,
                     sixteenth_decimated_picture_ptr->origin_y);
}

int av1_count_colors_highbd(uint16_t *src, int stride, int rows, int cols, int bit_depth,
                            int *val_count) {
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

int eb_av1_count_colors(const uint8_t *src, int stride, int rows, int cols, int *val_count) {
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
extern AomVarianceFnPtr mefn_ptr[BlockSizeS_ALL];

// This is used as a reference when computing the source variance for the
//  purposes of activity masking.
// Eventually this should be replaced by custom no-reference routines,
//  which will be faster.
const uint8_t eb_av1_var_offs[MAX_SB_SIZE] = {
    128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
    128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128};

static const uint16_t eb_AV1_HIGH_VAR_OFFS_10[MAX_SB_SIZE] = {
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4,
  128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4, 128 * 4
};

unsigned int eb_av1_get_sby_perpixel_variance(const AomVarianceFnPtr *fn_ptr, //const AV1_COMP *cpi,
                                              const uint8_t *         src,
                                              int       stride, //const struct Buf2D *ref,
                                              BlockSize bs) {
    unsigned int       sse;
    const unsigned int var =
        //cpi->fn_ptr[bs].vf(ref->buf, ref->stride, eb_av1_var_offs, 0, &sse);
        fn_ptr->vf(src, stride, eb_av1_var_offs, 0, &sse);
    return ROUND_POWER_OF_TWO(var, num_pels_log2_lookup[bs]);
}

unsigned int eb_av1_high_get_sby_perpixel_variance(const AomVarianceFnPtr *fn_ptr,
                                                   const uint16_t *src, int stride,
                                                   BlockSize bs) {
  unsigned int sse;
  const unsigned int var =
     fn_ptr->vf_hbd_10(CONVERT_TO_BYTEPTR(src), stride,
                       CONVERT_TO_BYTEPTR(eb_AV1_HIGH_VAR_OFFS_10),
                       0, &sse);
  return ROUND_POWER_OF_TWO(var, num_pels_log2_lookup[bs]);
}


// Check if the number of color of a block is superior to 1 and inferior
// to a given threshold.
EbBool is_valid_palette_nb_colors(const uint8_t *src, int stride,
                                  int rows, int cols,
                                  int nb_colors_threshold)
{
    EbBool has_color[1 << 8];  // Maximum (1 << 8) color levels.
    memset(has_color, 0, (1 << 8) * sizeof(*has_color));
    int nb_colors = 0;

    for (int r = 0; r < rows; ++r)
    {
        for (int c = 0; c < cols; ++c)
        {
            const int this_val = src[r * stride + c];
            if (has_color[this_val] == 0)
            {
                has_color[this_val] = 1;
                nb_colors++;
                if (nb_colors > nb_colors_threshold)
                    return EB_FALSE;
            }
        }
    }
    if (nb_colors <= 1)
        return EB_FALSE;

    return EB_TRUE;
}

EbBool is_valid_palette_nb_colors_highbd(const uint16_t *src, int stride,
                                         int rows, int cols,
                                         int nb_colors_threshold,
                                         int bit_depth)
{
    EbBool has_color[1 << 12]; // Maximum (1 << 12) color levels.
    memset(has_color, 0, ((size_t)1 << bit_depth) * sizeof(*has_color));
    int nb_colors = 0;

    for (int r = 0; r < rows; ++r)
    {
        for (int c = 0; c < cols; ++c)
        {
            const int this_val = src[r * stride + c];
            if (has_color[this_val] == 0)
            {
                has_color[this_val] = 1;
                nb_colors++;
                if (nb_colors > nb_colors_threshold)
                    return EB_FALSE;
            }
        }
    }
    if (nb_colors <= 1)
        return EB_FALSE;

    return EB_TRUE;
}


// Estimate if the source frame is screen content, based on the portion of
// blocks that have no more than 4 (experimentally selected) luma colors.
static void is_screen_content(PictureParentControlSet *pcs_ptr, int bit_depth) {
    const int blk_w = 16;
    const int blk_h = 16;
    // These threshold values are selected experimentally.
    const int          color_thresh = 4;
    const int          var_thresh = 0;
    // Counts of blocks with no more than color_thresh colors.
    int counts_1 = 0;
    // Counts of blocks with no more than color_thresh colors and variance larger
    // than var_thresh.
    int counts_2 = 0;

    const AomVarianceFnPtr *fn_ptr = &mefn_ptr[BLOCK_16X16];

    EbBool is16bit = bit_depth > EB_8BIT;
    EbPictureBufferDesc *input_picture_ptr =
        pcs_ptr->enhanced_picture_ptr;

    for (int r = 0; r + blk_h <= input_picture_ptr->height; r += blk_h)
    {
        for (int c = 0; c + blk_w <= input_picture_ptr->width; c += blk_w)
        {
            if (is16bit)
            {
                uint16_t src[BLOCK_SIZE_64 * BLOCK_SIZE_64];

                // In 16 bit mode, we must pack the split input in a temporary buffer.
                const uint32_t input_luma_offset =
                    (r + input_picture_ptr->origin_y) * input_picture_ptr->stride_y +
                    c + input_picture_ptr->origin_x;
                const uint32_t input_bit_inc_luma_offset =
                    (r + input_picture_ptr->origin_y) * input_picture_ptr->stride_bit_inc_y +
                    c + input_picture_ptr->origin_x;

                pack2d_src(
                    input_picture_ptr->buffer_y + input_luma_offset,
                    input_picture_ptr->stride_y,
                    input_picture_ptr->buffer_bit_inc_y + input_bit_inc_luma_offset,
                    input_picture_ptr->stride_bit_inc_y,
                    src,
                    blk_w,
                    blk_w,
                    blk_h);

                if (is_valid_palette_nb_colors_highbd((uint16_t *)src, blk_w,
                                                      blk_w, blk_h, color_thresh,
                                                      bit_depth))
                {
                    ++counts_1;
                    int var = eb_av1_high_get_sby_perpixel_variance(fn_ptr, src, blk_w,
                                                                    BLOCK_16X16);
                    if (var > var_thresh)
                        ++counts_2;
                }
            }
            else
            {
                uint8_t *src = input_picture_ptr->buffer_y +
                    (input_picture_ptr->origin_y + r) * input_picture_ptr->stride_y +
                    input_picture_ptr->origin_x + c;

                if (is_valid_palette_nb_colors(src, input_picture_ptr->stride_y,
                                               blk_w, blk_h, color_thresh))
                {
                    ++counts_1;
                    int var = eb_av1_get_sby_perpixel_variance(fn_ptr, src, input_picture_ptr->stride_y,
                                                               BLOCK_16X16);
                    if (var > var_thresh)
                        ++counts_2;
                }
            }
        }
    }

    pcs_ptr->sc_content_detected =
        (counts_1 * blk_h * blk_w * 10 > input_picture_ptr->width * input_picture_ptr->height) &&
        (counts_2 * blk_h * blk_w * 15 > input_picture_ptr->width * input_picture_ptr->height);
}

/************************************************
 * 1/4 & 1/16 input picture downsampling (filtering)
 ************************************************/
void downsample_filtering_input_picture(PictureParentControlSet *pcs_ptr,
                                        EbPictureBufferDesc *    input_padded_picture_ptr,
                                        EbPictureBufferDesc *    quarter_picture_ptr,
                                        EbPictureBufferDesc *    sixteenth_picture_ptr) {
    // Downsample input picture for HME L0 and L1
    if (pcs_ptr->enable_hme_flag || pcs_ptr->tf_enable_hme_flag) {
        if (pcs_ptr->enable_hme_level1_flag || pcs_ptr->tf_enable_hme_level1_flag) {
            downsample_2d(
                &input_padded_picture_ptr->buffer_y[input_padded_picture_ptr->origin_x +
                                                    input_padded_picture_ptr->origin_y *
                                                        input_padded_picture_ptr->stride_y],
                input_padded_picture_ptr->stride_y,
                input_padded_picture_ptr->width,
                input_padded_picture_ptr->height,
                &quarter_picture_ptr
                     ->buffer_y[quarter_picture_ptr->origin_x +
                                quarter_picture_ptr->origin_x * quarter_picture_ptr->stride_y],
                quarter_picture_ptr->stride_y,
                2);
            generate_padding(&quarter_picture_ptr->buffer_y[0],
                             quarter_picture_ptr->stride_y,
                             quarter_picture_ptr->width,
                             quarter_picture_ptr->height,
                             quarter_picture_ptr->origin_x,
                             quarter_picture_ptr->origin_y);
        }

        if (pcs_ptr->enable_hme_level0_flag || pcs_ptr->tf_enable_hme_level0_flag) {
            // Sixteenth Input Picture Downsampling
            if (pcs_ptr->enable_hme_level1_flag || pcs_ptr->tf_enable_hme_level1_flag)
                downsample_2d(
                    &quarter_picture_ptr
                         ->buffer_y[quarter_picture_ptr->origin_x +
                                    quarter_picture_ptr->origin_y * quarter_picture_ptr->stride_y],
                    quarter_picture_ptr->stride_y,
                    quarter_picture_ptr->width,
                    quarter_picture_ptr->height,
                    &sixteenth_picture_ptr->buffer_y[sixteenth_picture_ptr->origin_x +
                                                     sixteenth_picture_ptr->origin_x *
                                                         sixteenth_picture_ptr->stride_y],
                    sixteenth_picture_ptr->stride_y,
                    2);
            else
                downsample_2d(
                    &input_padded_picture_ptr->buffer_y[input_padded_picture_ptr->origin_x +
                                                        input_padded_picture_ptr->origin_y *
                                                            input_padded_picture_ptr->stride_y],
                    input_padded_picture_ptr->stride_y,
                    input_padded_picture_ptr->width,
                    input_padded_picture_ptr->height,
                    &sixteenth_picture_ptr->buffer_y[sixteenth_picture_ptr->origin_x +
                                                     sixteenth_picture_ptr->origin_x *
                                                         sixteenth_picture_ptr->stride_y],
                    sixteenth_picture_ptr->stride_y,
                    4);

            generate_padding(&sixteenth_picture_ptr->buffer_y[0],
                             sixteenth_picture_ptr->stride_y,
                             sixteenth_picture_ptr->width,
                             sixteenth_picture_ptr->height,
                             sixteenth_picture_ptr->origin_x,
                             sixteenth_picture_ptr->origin_y);
        }
    }
}

/* Picture Analysis Kernel */

/*********************************************************************************
*
* @brief
*  The Picture Analysis processes perform the first stage of encoder pre-processing analysis
*  as well as any intra-picture image conversion procedures, such as resampling, color space conversion, or tone mapping.
*
* @par Description:
*  The Picture Analysis processes can be multithreaded and as such can process multiple input pictures at a time.
*  The Picture Analysis also includes creating an n-bin Histogram, gathering 1st and 2nd moment statistics for each 8x8 block
*  in the picture, which are used in variance calculations. Since the Picture Analysis process is multithreaded,
*  the pictures can be processed out of order as long as all image-modifying functions are completed before any
*  statistics-gathering functions begin.
*
* @param[in] Pictures
*  The Picture Analysis Kernel performs pre-processing analysis as well as any intra-picture image conversion,
*  color space conversion or tone mapping on the pictures that it was given.
*
* @param[out] statistics
*  n-bin histogram is created to gather 1st and 2nd moment statistics for each 8x8 block which is then used to compute statistics
*
********************************************************************************/

void *picture_analysis_kernel(void *input_ptr) {
    EbThreadContext *        thread_context_ptr = (EbThreadContext *)input_ptr;
    PictureAnalysisContext * context_ptr = (PictureAnalysisContext *)thread_context_ptr->priv;
    PictureParentControlSet *pcs_ptr;
    SequenceControlSet *     scs_ptr;

    EbObjectWrapper *            in_results_wrapper_ptr;
    ResourceCoordinationResults *in_results_ptr;
    EbObjectWrapper *            out_results_wrapper_ptr;
    PictureAnalysisResults *     out_results_ptr;
    EbPaReferenceObject *        pa_ref_obj_;

    EbPictureBufferDesc *input_padded_picture_ptr;
    EbPictureBufferDesc *input_picture_ptr;

    // Variance
    uint32_t pic_width_in_sb;
    uint32_t pic_height_in_sb;
    uint32_t sb_total_count;

    for (;;) {
        // Get Input Full Object
        EB_GET_FULL_OBJECT(context_ptr->resource_coordination_results_input_fifo_ptr,
                           &in_results_wrapper_ptr);

        in_results_ptr = (ResourceCoordinationResults *)in_results_wrapper_ptr->object_ptr;
        pcs_ptr        = (PictureParentControlSet *)in_results_ptr->pcs_wrapper_ptr->object_ptr;

        // Mariana : save enhanced picture ptr, move this from here
        pcs_ptr->enhanced_unscaled_picture_ptr = pcs_ptr->enhanced_picture_ptr;

        // There is no need to do processing for overlay picture. Overlay and AltRef share the same results.
        if (!pcs_ptr->is_overlay) {
            scs_ptr           = (SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr;
            input_picture_ptr = pcs_ptr->enhanced_picture_ptr;

            pa_ref_obj_ =
                (EbPaReferenceObject *)pcs_ptr->pa_reference_picture_wrapper_ptr->object_ptr;
            input_padded_picture_ptr = (EbPictureBufferDesc *)pa_ref_obj_->input_padded_picture_ptr;
            // Variance
            pic_width_in_sb = (pcs_ptr->aligned_width + scs_ptr->sb_sz - 1) / scs_ptr->sb_sz;
            pic_height_in_sb = (pcs_ptr->aligned_height + scs_ptr->sb_sz - 1) / scs_ptr->sb_sz;
            sb_total_count = pic_width_in_sb * pic_height_in_sb;
            generate_padding(input_picture_ptr->buffer_y,
                             input_picture_ptr->stride_y,
                             input_picture_ptr->width,
                             input_picture_ptr->height,
                             input_picture_ptr->origin_x,
                             input_picture_ptr->origin_y);
            {
                uint8_t *pa =
                    input_padded_picture_ptr->buffer_y + input_padded_picture_ptr->origin_x +
                    input_padded_picture_ptr->origin_y * input_padded_picture_ptr->stride_y;
                uint8_t *in = input_picture_ptr->buffer_y + input_picture_ptr->origin_x +
                              input_picture_ptr->origin_y * input_picture_ptr->stride_y;
                for (uint32_t row = 0; row < input_picture_ptr->height; row++)
                    eb_memcpy(pa + row * input_padded_picture_ptr->stride_y,
                              in + row * input_picture_ptr->stride_y,
                              sizeof(uint8_t) * input_picture_ptr->width);
            }

            // Set picture parameters to account for subpicture, picture scantype, and set regions by resolutions
            set_picture_parameters_for_statistics_gathering(scs_ptr);

            // Pad pictures to multiple min cu size
            pad_picture_to_multiple_of_min_blk_size_dimensions(scs_ptr, input_picture_ptr);

            // Pre processing operations performed on the input picture
            picture_pre_processing_operations(pcs_ptr, scs_ptr, sb_total_count);
            if (input_picture_ptr->color_format >= EB_YUV422) {
                // Jing: Do the conversion of 422/444=>420 here since it's multi-threaded kernel
                //       Reuse the Y, only add cb/cr in the newly created buffer desc
                //       NOTE: since denoise may change the src, so this part is after picture_pre_processing_operations()
                pcs_ptr->chroma_downsampled_picture_ptr->buffer_y = input_picture_ptr->buffer_y;
                down_sample_chroma(input_picture_ptr, pcs_ptr->chroma_downsampled_picture_ptr);
            } else
                pcs_ptr->chroma_downsampled_picture_ptr = input_picture_ptr;
            // Pad input picture to complete border SBs
            pad_picture_to_multiple_of_sb_dimensions(input_padded_picture_ptr);
            // 1/4 & 1/16 input picture decimation
            downsample_decimation_input_picture(
                pcs_ptr,
                input_padded_picture_ptr,
                (EbPictureBufferDesc *)pa_ref_obj_->quarter_decimated_picture_ptr,
                (EbPictureBufferDesc *)pa_ref_obj_->sixteenth_decimated_picture_ptr);

            // 1/4 & 1/16 input picture downsampling through filtering
            if (scs_ptr->down_sampling_method_me_search == ME_FILTERED_DOWNSAMPLED) {
                downsample_filtering_input_picture(
                    pcs_ptr,
                    input_padded_picture_ptr,
                    (EbPictureBufferDesc *)pa_ref_obj_->quarter_filtered_picture_ptr,
                    (EbPictureBufferDesc *)pa_ref_obj_->sixteenth_filtered_picture_ptr);
            }

            // Gathering statistics of input picture, including Variance Calculation, Histogram Bins
            gathering_picture_statistics(
                scs_ptr,
                pcs_ptr,
                pcs_ptr->chroma_downsampled_picture_ptr, //420 input_picture_ptr
                input_padded_picture_ptr,
                (EbPictureBufferDesc *)pa_ref_obj_
                    ->sixteenth_decimated_picture_ptr, // Hsan: always use decimated until studying the trade offs
                sb_total_count);

            if (scs_ptr->static_config.screen_content_mode == 2) { // auto detect
                is_screen_content(pcs_ptr,
                                  scs_ptr->static_config.encoder_bit_depth);
            } else // off / on
                pcs_ptr->sc_content_detected = scs_ptr->static_config.screen_content_mode;

            // Hold the 64x64 variance and mean in the reference frame
            uint32_t sb_index;
            for (sb_index = 0; sb_index < pcs_ptr->sb_total_count; ++sb_index) {
                pa_ref_obj_->variance[sb_index] =
                    pcs_ptr->variance[sb_index][ME_TIER_ZERO_PU_64x64];
                pa_ref_obj_->y_mean[sb_index] = pcs_ptr->y_mean[sb_index][ME_TIER_ZERO_PU_64x64];
            }
        }
        // Get Empty Results Object
        eb_get_empty_object(context_ptr->picture_analysis_results_output_fifo_ptr,
                            &out_results_wrapper_ptr);

        out_results_ptr = (PictureAnalysisResults *)out_results_wrapper_ptr->object_ptr;
        out_results_ptr->pcs_wrapper_ptr = in_results_ptr->pcs_wrapper_ptr;

        // Release the Input Results
        eb_release_object(in_results_wrapper_ptr);

        // Post the Full Results Object
        eb_post_full_object(out_results_wrapper_ptr);
    }
    return EB_NULL;
}
