/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPictureAnalysis_h
#define EbPictureAnalysis_h

#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#include "EbNoiseExtractAVX2.h"

/**************************************
 * Context
 **************************************/
typedef struct PictureAnalysisContext
{
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
    EbPictureBufferDescInitData *input_picture_buffer_desc_init_data,
    EbBool                         denoise_flag,
    PictureAnalysisContext     **context_dbl_ptr,
    EbFifo                      *resource_coordination_results_input_fifo_ptr,
    EbFifo                      *picture_analysis_results_output_fifo_ptr);

extern void* picture_analysis_kernel(void *input_ptr);

void noise_extract_luma_weak(
    EbPictureBufferDesc *input_picture_ptr,
    EbPictureBufferDesc *denoised_picture_ptr,
    EbPictureBufferDesc *noise_picture_ptr,
    uint32_t               sb_origin_y,
    uint32_t               sb_origin_x);

typedef void(*EbWeakLumaFilterType)(
    EbPictureBufferDesc *input_picture_ptr,
    EbPictureBufferDesc *denoised_picture_ptr,
    EbPictureBufferDesc *noise_picture_ptr,
    uint32_t               sb_origin_y,
    uint32_t               sb_origin_x);

static EbWeakLumaFilterType FUNC_TABLE weak_luma_filter_func_ptr_array[ASM_TYPE_TOTAL] =
{
    // NON_AVX2
    noise_extract_luma_weak,
    // AVX2
    noise_extract_luma_weak_avx2_intrin,
};

void noise_extract_luma_weak_lcu(
    EbPictureBufferDesc *input_picture_ptr,
    EbPictureBufferDesc *denoised_picture_ptr,
    EbPictureBufferDesc *noise_picture_ptr,
    uint32_t               sb_origin_y,
    uint32_t               sb_origin_x);

static EbWeakLumaFilterType FUNC_TABLE weak_luma_filter_lcu_func_ptr_array[ASM_TYPE_TOTAL] =
{
    // NON_AVX2
    noise_extract_luma_weak_lcu,
    // AVX2
    noise_extract_luma_weak_lcu_avx2_intrin,
};

void noise_extract_luma_strong(
    EbPictureBufferDesc *input_picture_ptr,
    EbPictureBufferDesc *denoised_picture_ptr,
    uint32_t               sb_origin_y,
    uint32_t               sb_origin_x);

typedef void(*EbStrongLumaFilterType)(
    EbPictureBufferDesc *input_picture_ptr,
    EbPictureBufferDesc *denoised_picture_ptr,
    uint32_t               sb_origin_y,
    uint32_t               sb_origin_x);

static EbStrongLumaFilterType FUNC_TABLE strong_luma_filter_func_ptr_array[ASM_TYPE_TOTAL] =
{
    // NON_AVX2
    noise_extract_luma_strong,
    // AVX2
    noise_extract_luma_strong_avx2_intrin,
};
void noise_extract_chroma_strong(
    EbPictureBufferDesc *input_picture_ptr,
    EbPictureBufferDesc *denoised_picture_ptr,
    uint32_t               sb_origin_y,
    uint32_t               sb_origin_x);

typedef void(*EbStrongChromaFilterType)(
    EbPictureBufferDesc *input_picture_ptr,
    EbPictureBufferDesc *denoised_picture_ptr,
    uint32_t               sb_origin_y,
    uint32_t               sb_origin_x);

static EbStrongChromaFilterType FUNC_TABLE strong_chroma_filter_func_ptr_array[ASM_TYPE_TOTAL] =
{
    // NON_AVX2
    noise_extract_chroma_strong,
    // AVX2
    noise_extract_chroma_strong_avx2_intrin,
};

void noise_extract_chroma_weak(
    EbPictureBufferDesc *input_picture_ptr,
    EbPictureBufferDesc *denoised_picture_ptr,
    uint32_t               sb_origin_y,
    uint32_t               sb_origin_x);

typedef void(*EbWeakChromaFilterType)(
    EbPictureBufferDesc *input_picture_ptr,
    EbPictureBufferDesc *denoised_picture_ptr,
    uint32_t               sb_origin_y,
    uint32_t               sb_origin_x);

static EbWeakChromaFilterType FUNC_TABLE weak_chroma_filter_func_ptr_array[ASM_TYPE_TOTAL] =
{
    // NON_AVX2
    noise_extract_chroma_weak,
    // AVX2
    noise_extract_chroma_weak_avx2_intrin,
};

#endif // EbPictureAnalysis_h
