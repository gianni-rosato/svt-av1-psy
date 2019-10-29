/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPictureAnalysis_h
#define EbPictureAnalysis_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#include "EbNoiseExtractAVX2.h"
#include "EbObject.h"
#include "EbPictureControlSet.h"

/**************************************
 * Context
 **************************************/
typedef struct PictureAnalysisContext
{
    EbDctor                     dctor;
    EB_ALIGN(64) uint8_t            local_cache[64];
    EbFifo                     *resource_coordination_results_input_fifo_ptr;
    EbFifo                     *picture_analysis_results_output_fifo_ptr;
    EbPictureBufferDesc        *denoised_picture_ptr;
    EbPictureBufferDesc        *noise_picture_ptr;
    double                          pic_noise_variance_float;
} PictureAnalysisContext;

/***************************************
 * Extern Function Declaration
 ***************************************/
extern EbErrorType picture_analysis_context_ctor(
    PictureAnalysisContext     *context_ptr,
    EbPictureBufferDescInitData *input_picture_buffer_desc_init_data,
    EbBool                         denoise_flag,
    EbFifo                      *resource_coordination_results_input_fifo_ptr,
    EbFifo                      *picture_analysis_results_output_fifo_ptr);

extern void* picture_analysis_kernel(void *input_ptr);

void noise_extract_luma_weak_c(
    EbPictureBufferDesc *input_picture_ptr,
    EbPictureBufferDesc *denoised_picture_ptr,
    EbPictureBufferDesc *noise_picture_ptr,
    uint32_t               sb_origin_y,
    uint32_t               sb_origin_x);

void DownsampleFilteringInputPicture(
    PictureParentControlSet       *picture_control_set_ptr,
    EbPictureBufferDesc           *input_padded_picture_ptr,
    EbPictureBufferDesc           *quarter_picture_ptr,
    EbPictureBufferDesc *sixteenth_picture_ptr);

void noise_extract_luma_weak_lcu_c(
    EbPictureBufferDesc *input_picture_ptr,
    EbPictureBufferDesc *denoised_picture_ptr,
    EbPictureBufferDesc *noise_picture_ptr,
    uint32_t               sb_origin_y,
    uint32_t               sb_origin_x);

void noise_extract_luma_strong_c(
    EbPictureBufferDesc *input_picture_ptr,
    EbPictureBufferDesc *denoised_picture_ptr,
    uint32_t               sb_origin_y,
    uint32_t               sb_origin_x);

void noise_extract_chroma_strong_c(
    EbPictureBufferDesc *input_picture_ptr,
    EbPictureBufferDesc *denoised_picture_ptr,
    uint32_t               sb_origin_y,
    uint32_t               sb_origin_x);

void noise_extract_chroma_weak_c(
    EbPictureBufferDesc *input_picture_ptr,
    EbPictureBufferDesc *denoised_picture_ptr,
    uint32_t               sb_origin_y,
    uint32_t               sb_origin_x);

#endif // EbPictureAnalysis_h
