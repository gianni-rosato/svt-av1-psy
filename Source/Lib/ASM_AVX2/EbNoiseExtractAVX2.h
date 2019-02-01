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
    * noiseExtractLumaWeak
    *  weak filter Luma and store noise.
    *******************************************/
    void noiseExtractLumaWeak_AVX2_INTRIN(
        EbPictureBufferDesc_t       *inputPicturePtr,
        EbPictureBufferDesc_t       *denoisedPicturePtr,
        EbPictureBufferDesc_t       *noisePicturePtr,
        uint32_t                       sb_origin_y,
        uint32_t                       sb_origin_x
    );

    void noiseExtractLumaWeakLcu_AVX2_INTRIN(
        EbPictureBufferDesc_t       *inputPicturePtr,
        EbPictureBufferDesc_t       *denoisedPicturePtr,
        EbPictureBufferDesc_t       *noisePicturePtr,
        uint32_t                       sb_origin_y,
        uint32_t                       sb_origin_x
    );

    void noiseExtractChromaStrong_AVX2_INTRIN(
        EbPictureBufferDesc_t       *inputPicturePtr,
        EbPictureBufferDesc_t       *denoisedPicturePtr,
        uint32_t                       sb_origin_y,
        uint32_t                         sb_origin_x);

    void noiseExtractChromaWeak_AVX2_INTRIN(
        EbPictureBufferDesc_t       *inputPicturePtr,
        EbPictureBufferDesc_t       *denoisedPicturePtr,
        uint32_t                       sb_origin_y,
        uint32_t                         sb_origin_x);

    void noiseExtractLumaStrong_AVX2_INTRIN(
        EbPictureBufferDesc_t       *inputPicturePtr,
        EbPictureBufferDesc_t       *denoisedPicturePtr,
        uint32_t                       sb_origin_y,
        uint32_t                       sb_origin_x);

    void ChromaStrong_AVX2_INTRIN(
        __m256i                        top,
        __m256i                        curr,
        __m256i                        bottom,
        __m256i                        currPrev,
        __m256i                        currNext,
        __m256i                        topPrev,
        __m256i                        topNext,
        __m256i                        bottomPrev,
        __m256i                        bottomNext,
        uint8_t                       *ptrDenoised);

    void lumaWeakFilter_AVX2_INTRIN(
        __m256i                        top,
        __m256i                        curr,
        __m256i                        bottom,
        __m256i                        currPrev,
        __m256i                        currNext,
        uint8_t                       *ptrDenoised,
        uint8_t                        *ptrNoise);


    void chromaWeakLumaStrongFilter_AVX2_INTRIN(
        __m256i                        top,
        __m256i                        curr,
        __m256i                        bottom,
        __m256i                        currPrev,
        __m256i                        currNext,
        __m256i                        topPrev,
        __m256i                        topNext,
        __m256i                        bottomPrev,
        __m256i                        bottomNext,
        uint8_t                       *ptrDenoised);

#ifdef __cplusplus
}
#endif
#endif // EbNoiseExtractAVX2_h

