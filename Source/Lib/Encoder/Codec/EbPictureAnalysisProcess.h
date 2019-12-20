// clang-format off
/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPictureAnalysis_h
#define EbPictureAnalysis_h

#include "EbDefinitions.h"
#include "EbNoiseExtractAVX2.h"
#include "EbPictureControlSet.h"

/***************************************
 * Extern Function Declaration
 ***************************************/
EbErrorType picture_analysis_context_ctor(
    EbThreadContext   *thread_context_ptr,
    const EbEncHandle *enc_handle_ptr,
    int               index);

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
// clang-format on
