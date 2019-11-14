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
    void noise_extract_luma_weak_avx2_intrin(
        EbPictureBufferDesc *input_picture_ptr,
        EbPictureBufferDesc *denoised_picture_ptr,
        EbPictureBufferDesc *noise_picture_ptr,
        uint32_t               sb_origin_y,
        uint32_t               sb_origin_x);

    void noise_extract_luma_weak_lcu_avx2_intrin(
        EbPictureBufferDesc *input_picture_ptr,
        EbPictureBufferDesc *denoised_picture_ptr,
        EbPictureBufferDesc *noise_picture_ptr,
        uint32_t               sb_origin_y,
        uint32_t               sb_origin_x);

    void noise_extract_chroma_strong_avx2_intrin(
        EbPictureBufferDesc *input_picture_ptr,
        EbPictureBufferDesc *denoised_picture_ptr,
        uint32_t               sb_origin_y,
        uint32_t               sb_origin_x);

    void noise_extract_chroma_weak_avx2_intrin(
        EbPictureBufferDesc *input_picture_ptr,
        EbPictureBufferDesc *denoised_picture_ptr,
        uint32_t               sb_origin_y,
        uint32_t               sb_origin_x);

    void noise_extract_luma_strong_avx2_intrin(
        EbPictureBufferDesc *input_picture_ptr,
        EbPictureBufferDesc *denoised_picture_ptr,
        uint32_t               sb_origin_y,
        uint32_t               sb_origin_x);

    void chroma_strong_avx2_intrin(
        __m256i   top,
        __m256i   curr,
        __m256i   bottom,
        __m256i   curr_prev,
        __m256i   curr_next,
        __m256i   top_prev,
        __m256i   top_next,
        __m256i   bottom_prev,
        __m256i   bottom_next,
        uint8_t  *ptr_denoised);

    void luma_weak_filter_avx2_intrin(
        __m256i  top,
        __m256i  curr,
        __m256i  bottom,
        __m256i  curr_prev,
        __m256i  curr_next,
        uint8_t *ptr_denoised,
        uint8_t *ptr_noise);

    void chroma_weak_luma_strong_filter_avx2_intrin(
        __m256i  top,
        __m256i  curr,
        __m256i  bottom,
        __m256i  curr_prev,
        __m256i  curr_next,
        __m256i  top_prev,
        __m256i  top_next,
        __m256i  bottom_prev,
        __m256i  bottom_next,
        uint8_t *ptr_denoised);

    void luma_weak_filter_128_avx2_intrin(
        __m128i                        top,
        __m128i                        curr,
        __m128i                        bottom,
        __m128i                        curr_prev,
        __m128i                        curr_next,
        uint8_t* ptr_denoised,
        uint8_t* ptr_noise);

    void chroma_strong_128_avx2_intrin(
        __m128i                        top,
        __m128i                        curr,
        __m128i                        bottom,
        __m128i                        curr_prev,
        __m128i                        curr_next,
        __m128i                        top_prev,
        __m128i                        top_next,
        __m128i                        bottom_prev,
        __m128i                        bottom_next,
        uint8_t* ptr_denoised);

    void chroma_weak_luma_strong_filter_128_avx2_intrin(
        __m128i                        top,
        __m128i                        curr,
        __m128i                        bottom,
        __m128i                        curr_prev,
        __m128i                        curr_next,
        __m128i                        top_prev,
        __m128i                        top_next,
        __m128i                        bottom_prev,
        __m128i                        bottom_next,
        uint8_t* ptr_denoised);

#ifdef __cplusplus
}
#endif
#endif // EbNoiseExtractAVX2_h
