/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbNoiseExtractAVX2_h
#define EbNoiseExtractAVX2_h

#include "immintrin.h"
#include "EbDefinitions.h"
#include "EbPictureBufferDesc.h"
#ifdef __cplusplus
extern "C" {
#endif

/*******************************************
    * noise_extract_luma_weak
    *  weak filter Luma and store noise.
    *******************************************/
void noise_extract_luma_weak_avx2_intrin(EbPictureBufferDesc *input_picture_ptr,
                                         EbPictureBufferDesc *denoised_picture_ptr,
                                         EbPictureBufferDesc *noise_picture_ptr,
                                         uint32_t sb_origin_y, uint32_t sb_origin_x);

void noise_extract_luma_weak_sb_avx2_intrin(EbPictureBufferDesc *input_picture_ptr,
                                            EbPictureBufferDesc *denoised_picture_ptr,
                                            EbPictureBufferDesc *noise_picture_ptr,
                                            uint32_t sb_origin_y, uint32_t sb_origin_x);

void noise_extract_chroma_strong_avx2_intrin(EbPictureBufferDesc *input_picture_ptr,
                                             EbPictureBufferDesc *denoised_picture_ptr,
                                             uint32_t sb_origin_y, uint32_t sb_origin_x);

void noise_extract_chroma_weak_avx2_intrin(EbPictureBufferDesc *input_picture_ptr,
                                           EbPictureBufferDesc *denoised_picture_ptr,
                                           uint32_t sb_origin_y, uint32_t sb_origin_x);

void noise_extract_luma_strong_avx2_intrin(EbPictureBufferDesc *input_picture_ptr,
                                           EbPictureBufferDesc *denoised_picture_ptr,
                                           uint32_t sb_origin_y, uint32_t sb_origin_x);

#ifdef __cplusplus
}
#endif
#endif // EbNoiseExtractAVX2_h
