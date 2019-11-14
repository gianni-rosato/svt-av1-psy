/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbNoiseExtractAVX2.h"
#include "EbDefinitions.h"
#include "immintrin.h"
#include "EbUtility.h"

EB_EXTERN EB_ALIGN(16) const uint8_t filterType[] = {
    1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4
};

EB_EXTERN EB_ALIGN(16) const uint8_t WeakChromafilter[2][32] = {
        { 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4 },
        { 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2 },
};

inline void luma_weak_filter_avx2_intrin(
    __m256i                        top,
    __m256i                        curr,
    __m256i                        bottom,
    __m256i                        curr_prev,
    __m256i                        curr_next,
    uint8_t                       *ptr_denoised,
    uint8_t                        *ptr_noise
)
{
    __m256i  topFirstHalf, bottomFirstHalf,
        filterFirstHalf, filterSecondHalf,
        currNextFirstHalf, currNextSecondHalf,
        weights, currLeftMidFirstHalfWeight,
        currLeftMidFirstHalflo, currLeftMidFirstHalfhi, currPrevPermutation, currPermutation, currNextPermutation,
        topPermutation, bottomPermutation;

    currPrevPermutation = _mm256_permute4x64_epi64(curr_prev, 216);
    currPermutation = _mm256_permute4x64_epi64(curr, 216);
    currLeftMidFirstHalflo = _mm256_unpacklo_epi8(currPrevPermutation, currPermutation);
    weights = _mm256_loadu_si256((__m256i*)filterType);
    currLeftMidFirstHalfWeight = _mm256_maddubs_epi16(currLeftMidFirstHalflo, weights);
    currNextPermutation = _mm256_permute4x64_epi64(curr_next, 88);
    currNextFirstHalf = _mm256_unpacklo_epi8(currNextPermutation, _mm256_setzero_si256());
    currLeftMidFirstHalflo = _mm256_add_epi16(currNextFirstHalf, currLeftMidFirstHalfWeight);

    currLeftMidFirstHalfhi = _mm256_unpackhi_epi8(currPrevPermutation, currPermutation);
    currLeftMidFirstHalfWeight = _mm256_maddubs_epi16(currLeftMidFirstHalfhi, weights);
    currNextPermutation = _mm256_permute4x64_epi64(curr_next, 216);
    currNextSecondHalf = _mm256_unpackhi_epi8(currNextPermutation, _mm256_setzero_si256());
    currLeftMidFirstHalfhi = _mm256_add_epi16(currNextSecondHalf, currLeftMidFirstHalfWeight);

    topPermutation = _mm256_permute4x64_epi64(top, 216);
    topFirstHalf = _mm256_unpacklo_epi8(topPermutation, _mm256_setzero_si256());
    bottomPermutation = _mm256_permute4x64_epi64(bottom, 216);
    bottomFirstHalf = _mm256_unpacklo_epi8(bottomPermutation, _mm256_setzero_si256());
    filterFirstHalf = _mm256_adds_epi16(_mm256_adds_epi16(bottomFirstHalf, topFirstHalf), currLeftMidFirstHalflo);
    filterFirstHalf = _mm256_srli_epi16(filterFirstHalf, 3);

    topFirstHalf = _mm256_unpackhi_epi8(topPermutation, _mm256_setzero_si256());
    bottomFirstHalf = _mm256_unpackhi_epi8(bottomPermutation, _mm256_setzero_si256());
    filterSecondHalf = _mm256_adds_epi16(_mm256_adds_epi16(bottomFirstHalf, topFirstHalf), currLeftMidFirstHalfhi);
    filterSecondHalf = _mm256_srli_epi16(filterSecondHalf, 3);

    filterFirstHalf = _mm256_permute4x64_epi64(_mm256_packus_epi16(filterFirstHalf, filterSecondHalf), 216);
    _mm256_storeu_si256((__m256i *)(ptr_denoised), filterFirstHalf);

    _mm256_storeu_si256((__m256i *)(ptr_noise), _mm256_subs_epu8(curr, filterFirstHalf));
}
inline void luma_weak_filter_128_avx2_intrin(
    __m128i                        top,
    __m128i                        curr,
    __m128i                        bottom,
    __m128i                        curr_prev,
    __m128i                        curr_next,
    uint8_t                       *ptr_denoised,
    uint8_t                       *ptr_noise
)
{
    __m128i topFirstHalf, bottomFirstHalf,
        filterFirstHalf, filterSecondHalf,
        currNextFirstHalf, currNextSecondHalf,
        weights, currLeftMidFirstHalfWeight,
        currLeftMidFirstHalflo, currLeftMidFirstHalfhi;

    currLeftMidFirstHalflo = _mm_unpacklo_epi8(curr_prev, curr);
    weights = _mm_loadu_si128((__m128i*)filterType);
    currLeftMidFirstHalfWeight = _mm_maddubs_epi16(currLeftMidFirstHalflo, weights);
    currNextFirstHalf = _mm_unpacklo_epi8(curr_next, _mm_setzero_si128());
    currLeftMidFirstHalflo = _mm_add_epi16(currNextFirstHalf, currLeftMidFirstHalfWeight);

    currLeftMidFirstHalfhi = _mm_unpackhi_epi8(curr_prev, curr);
    currLeftMidFirstHalfWeight = _mm_maddubs_epi16(currLeftMidFirstHalfhi, weights);
    currNextSecondHalf = _mm_unpackhi_epi8(curr_next, _mm_setzero_si128());
    currLeftMidFirstHalfhi = _mm_add_epi16(currNextSecondHalf, currLeftMidFirstHalfWeight);

    topFirstHalf = _mm_unpacklo_epi8(top, _mm_setzero_si128());
    bottomFirstHalf = _mm_unpacklo_epi8(bottom, _mm_setzero_si128());
    filterFirstHalf = _mm_adds_epi16(_mm_adds_epi16(bottomFirstHalf, topFirstHalf), currLeftMidFirstHalflo);
    filterFirstHalf = _mm_srli_epi16(filterFirstHalf, 3);

    topFirstHalf = _mm_unpackhi_epi8(top, _mm_setzero_si128());
    bottomFirstHalf = _mm_unpackhi_epi8(bottom, _mm_setzero_si128());
    filterSecondHalf = _mm_adds_epi16(_mm_adds_epi16(bottomFirstHalf, topFirstHalf), currLeftMidFirstHalfhi);
    filterSecondHalf = _mm_srli_epi16(filterSecondHalf, 3);

    filterFirstHalf = _mm_packus_epi16(filterFirstHalf, filterSecondHalf);
    _mm_storel_epi64((__m128i*)(ptr_denoised), filterFirstHalf);

    _mm_storel_epi64((__m128i*)(ptr_noise), _mm_subs_epu8(curr, filterFirstHalf));
}
inline void chroma_weak_luma_strong_filter_avx2_intrin(
    __m256i                        top,
    __m256i                        curr,
    __m256i                        bottom,
    __m256i                        curr_prev,
    __m256i                        curr_next,
    __m256i                        top_prev,
    __m256i                        top_next,
    __m256i                        bottom_prev,
    __m256i                        bottom_next,
    uint8_t                       *ptr_denoised
)
{
    __m256i filterFirstHalf, filterSecondHalf,
        currNextFirstHalf, currNextSecondHalf,
        weights, currLeftMidFirstHalfWeight,
        currLeftMidFirstHalflo, currLeftMidFirstHalfhi, currPrevPermutation, currPermutation, currNextPermutation,
        topPermutation, bottomPermutation,
        topPrevPermutation, topLeftMidFirstHalflo, topLeftMidFirstHalfWeight, topNextFirstHalf,
        topNextPermutation, topLeftMidFirstHalfhi, topNextSecondHalf,
        bottomPrevPermutation, bottomLeftMidFirstHalflo, bottomLeftMidFirstHalfWeight, bottomNextPermutation,
        bottomNextFirstHalf, bottomLeftMidFirstHalfhi, bottomNextSecondHalf;

    //  Curr
    currPrevPermutation = _mm256_permute4x64_epi64(curr_prev, 216);
    currPermutation = _mm256_permute4x64_epi64(curr, 216);
    currLeftMidFirstHalflo = _mm256_unpacklo_epi8(currPrevPermutation, currPermutation);
    weights = _mm256_loadu_si256((__m256i*)WeakChromafilter[0]);
    currLeftMidFirstHalfWeight = _mm256_maddubs_epi16(currLeftMidFirstHalflo, weights);
    currNextPermutation = _mm256_permute4x64_epi64(curr_next, 88);
    currNextFirstHalf = _mm256_unpacklo_epi8(currNextPermutation, _mm256_setzero_si256());
    currNextFirstHalf = _mm256_slli_epi16(currNextFirstHalf, 1);
    currLeftMidFirstHalflo = _mm256_add_epi16(currNextFirstHalf, currLeftMidFirstHalfWeight);

    currLeftMidFirstHalfhi = _mm256_unpackhi_epi8(currPrevPermutation, currPermutation);
    currLeftMidFirstHalfWeight = _mm256_maddubs_epi16(currLeftMidFirstHalfhi, weights);
    currNextPermutation = _mm256_permute4x64_epi64(curr_next, 216);
    currNextSecondHalf = _mm256_unpackhi_epi8(currNextPermutation, _mm256_setzero_si256());
    currNextSecondHalf = _mm256_slli_epi16(currNextSecondHalf, 1);
    currLeftMidFirstHalfhi = _mm256_add_epi16(currNextSecondHalf, currLeftMidFirstHalfWeight);

    // Top
    topPrevPermutation = _mm256_permute4x64_epi64(top_prev, 216);
    topPermutation = _mm256_permute4x64_epi64(top, 216);
    topLeftMidFirstHalflo = _mm256_unpacklo_epi8(topPrevPermutation, topPermutation);
    weights = _mm256_loadu_si256((__m256i*)WeakChromafilter[1]);
    topLeftMidFirstHalfWeight = _mm256_maddubs_epi16(topLeftMidFirstHalflo, weights);
    topNextPermutation = _mm256_permute4x64_epi64(top_next, 88);
    topNextFirstHalf = _mm256_unpacklo_epi8(topNextPermutation, _mm256_setzero_si256());
    topLeftMidFirstHalflo = _mm256_add_epi16(topNextFirstHalf, topLeftMidFirstHalfWeight);

    topLeftMidFirstHalfhi = _mm256_unpackhi_epi8(topPrevPermutation, topPermutation);
    topLeftMidFirstHalfWeight = _mm256_maddubs_epi16(topLeftMidFirstHalfhi, weights);
    topNextPermutation = _mm256_permute4x64_epi64(top_next, 216);
    topNextSecondHalf = _mm256_unpackhi_epi8(topNextPermutation, _mm256_setzero_si256());
    topLeftMidFirstHalfhi = _mm256_add_epi16(topNextSecondHalf, topLeftMidFirstHalfWeight);

    // Bottom
    bottomPrevPermutation = _mm256_permute4x64_epi64(bottom_prev, 216);
    bottomPermutation = _mm256_permute4x64_epi64(bottom, 216);
    bottomLeftMidFirstHalflo = _mm256_unpacklo_epi8(bottomPrevPermutation, bottomPermutation);
    weights = _mm256_loadu_si256((__m256i*)WeakChromafilter[1]);
    bottomLeftMidFirstHalfWeight = _mm256_maddubs_epi16(bottomLeftMidFirstHalflo, weights);
    bottomNextPermutation = _mm256_permute4x64_epi64(bottom_next, 88);
    bottomNextFirstHalf = _mm256_unpacklo_epi8(bottomNextPermutation, _mm256_setzero_si256());
    bottomLeftMidFirstHalflo = _mm256_add_epi16(bottomNextFirstHalf, bottomLeftMidFirstHalfWeight);

    bottomLeftMidFirstHalfhi = _mm256_unpackhi_epi8(bottomPrevPermutation, bottomPermutation);
    bottomLeftMidFirstHalfWeight = _mm256_maddubs_epi16(bottomLeftMidFirstHalfhi, weights);
    bottomNextPermutation = _mm256_permute4x64_epi64(bottom_next, 216);
    bottomNextSecondHalf = _mm256_unpackhi_epi8(bottomNextPermutation, _mm256_setzero_si256());
    bottomLeftMidFirstHalfhi = _mm256_add_epi16(bottomNextSecondHalf, bottomLeftMidFirstHalfWeight);

    filterFirstHalf = _mm256_adds_epi16(_mm256_adds_epi16(bottomLeftMidFirstHalflo, topLeftMidFirstHalflo), currLeftMidFirstHalflo);
    filterFirstHalf = _mm256_srli_epi16(filterFirstHalf, 4);
    filterSecondHalf = _mm256_adds_epi16(_mm256_adds_epi16(bottomLeftMidFirstHalfhi, topLeftMidFirstHalfhi), currLeftMidFirstHalfhi);
    filterSecondHalf = _mm256_srli_epi16(filterSecondHalf, 4);

    filterFirstHalf = _mm256_permute4x64_epi64(_mm256_packus_epi16(filterFirstHalf, filterSecondHalf), 216);
    _mm256_storeu_si256((__m256i *)(ptr_denoised), filterFirstHalf);
}

inline void chroma_weak_luma_strong_filter_128_avx2_intrin(
    __m128i                        top,
    __m128i                        curr,
    __m128i                        bottom,
    __m128i                        curr_prev,
    __m128i                        curr_next,
    __m128i                        top_prev,
    __m128i                        top_next,
    __m128i                        bottom_prev,
    __m128i                        bottom_next,
    uint8_t                       *ptr_denoised
)
{
    __m128i filterFirstHalf, filterSecondHalf,
        currNextFirstHalf, currNextSecondHalf,
        weights, currLeftMidFirstHalfWeight,
        currLeftMidFirstHalflo, currLeftMidFirstHalfhi,
        topLeftMidFirstHalflo, topLeftMidFirstHalfWeight, topNextFirstHalf,
        topLeftMidFirstHalfhi, topNextSecondHalf,
        bottomLeftMidFirstHalflo, bottomLeftMidFirstHalfWeight,
        bottomNextFirstHalf, bottomLeftMidFirstHalfhi, bottomNextSecondHalf;

    //  Curr
    currLeftMidFirstHalflo = _mm_unpacklo_epi8(curr_prev, curr);
    weights = _mm_loadu_si128((__m128i*)WeakChromafilter[0]);
    currLeftMidFirstHalfWeight = _mm_maddubs_epi16(currLeftMidFirstHalflo, weights);
    currNextFirstHalf = _mm_unpacklo_epi8(curr_next, _mm_setzero_si128());
    currNextFirstHalf = _mm_slli_epi16(currNextFirstHalf, 1);
    currLeftMidFirstHalflo = _mm_add_epi16(currNextFirstHalf, currLeftMidFirstHalfWeight);

    currLeftMidFirstHalfhi = _mm_unpackhi_epi8(curr_prev, curr);
    currLeftMidFirstHalfWeight = _mm_maddubs_epi16(currLeftMidFirstHalfhi, weights);
    currNextSecondHalf = _mm_unpackhi_epi8(curr_next, _mm_setzero_si128());
    currNextSecondHalf = _mm_slli_epi16(currNextSecondHalf, 1);
    currLeftMidFirstHalfhi = _mm_add_epi16(currNextSecondHalf, currLeftMidFirstHalfWeight);

    // Top
    topLeftMidFirstHalflo = _mm_unpacklo_epi8(top_prev, top);
    weights = _mm_loadu_si128((__m128i*)WeakChromafilter[1]);
    topLeftMidFirstHalfWeight = _mm_maddubs_epi16(topLeftMidFirstHalflo, weights);
    topNextFirstHalf = _mm_unpacklo_epi8(top_next, _mm_setzero_si128());
    topLeftMidFirstHalflo = _mm_add_epi16(topNextFirstHalf, topLeftMidFirstHalfWeight);

    topLeftMidFirstHalfhi = _mm_unpackhi_epi8(top_prev, top);
    topLeftMidFirstHalfWeight = _mm_maddubs_epi16(topLeftMidFirstHalfhi, weights);
    topNextSecondHalf = _mm_unpackhi_epi8(top_next, _mm_setzero_si128());
    topLeftMidFirstHalfhi = _mm_add_epi16(topNextSecondHalf, topLeftMidFirstHalfWeight);

    // Bottom
    bottomLeftMidFirstHalflo = _mm_unpacklo_epi8(bottom_prev, bottom);
    weights = _mm_loadu_si128((__m128i*)WeakChromafilter[1]);
    bottomLeftMidFirstHalfWeight = _mm_maddubs_epi16(bottomLeftMidFirstHalflo, weights);
    bottomNextFirstHalf = _mm_unpacklo_epi8(bottom_next, _mm_setzero_si128());
    bottomLeftMidFirstHalflo = _mm_add_epi16(bottomNextFirstHalf, bottomLeftMidFirstHalfWeight);

    bottomLeftMidFirstHalfhi = _mm_unpackhi_epi8(bottom_prev, bottom);
    bottomLeftMidFirstHalfWeight = _mm_maddubs_epi16(bottomLeftMidFirstHalfhi, weights);
    bottomNextSecondHalf = _mm_unpackhi_epi8(bottom_next, _mm_setzero_si128());
    bottomLeftMidFirstHalfhi = _mm_add_epi16(bottomNextSecondHalf, bottomLeftMidFirstHalfWeight);

    filterFirstHalf = _mm_adds_epi16(_mm_adds_epi16(bottomLeftMidFirstHalflo, topLeftMidFirstHalflo), currLeftMidFirstHalflo);
    filterFirstHalf = _mm_srli_epi16(filterFirstHalf, 4);
    filterSecondHalf = _mm_adds_epi16(_mm_adds_epi16(bottomLeftMidFirstHalfhi, topLeftMidFirstHalfhi), currLeftMidFirstHalfhi);
    filterSecondHalf = _mm_srli_epi16(filterSecondHalf, 4);

    filterFirstHalf = _mm_packus_epi16(filterFirstHalf, filterSecondHalf);
    _mm_storel_epi64((__m128i*)(ptr_denoised), filterFirstHalf);
}

inline void chroma_strong_avx2_intrin(
    __m256i                        top,
    __m256i                        curr,
    __m256i                        bottom,
    __m256i                        curr_prev,
    __m256i                        curr_next,
    __m256i                        top_prev,
    __m256i                        top_next,
    __m256i                        bottom_prev,
    __m256i                        bottom_next,
    uint8_t                       *ptr_denoised
)
{
    __m256i   currLeftMidFirstHalflo, currLeftMidFirstHalfhi, currPrevPermutation, currPermutation, currNextPermutation,
        topPermutation, topPrevPermutation, topLeftMidFirstHalflo, topNextPermutation, topLeftMidFirstHalfhi,
        bottomPermutation, bottomPrevPermutation, bottomLeftMidFirstHalflo, bottomNextPermutation, bottomLeftMidFirstHalfhi;

    currPrevPermutation = _mm256_permute4x64_epi64(curr_prev, 216);
    currPermutation = _mm256_permute4x64_epi64(curr, 216);
    currNextPermutation = _mm256_permute4x64_epi64(curr_next, 216);

    currLeftMidFirstHalflo = _mm256_add_epi16(_mm256_unpacklo_epi8(currPermutation, _mm256_setzero_si256()),
        _mm256_unpacklo_epi8(currPrevPermutation, _mm256_setzero_si256()));
    currLeftMidFirstHalflo = _mm256_add_epi16(_mm256_unpacklo_epi8(currNextPermutation, _mm256_setzero_si256()), currLeftMidFirstHalflo);

    currLeftMidFirstHalfhi = _mm256_add_epi16(_mm256_unpackhi_epi8(currPermutation, _mm256_setzero_si256()),
        _mm256_unpackhi_epi8(currPrevPermutation, _mm256_setzero_si256()));
    currLeftMidFirstHalfhi = _mm256_add_epi16(_mm256_unpackhi_epi8(currNextPermutation, _mm256_setzero_si256()), currLeftMidFirstHalfhi);

    topPrevPermutation = _mm256_permute4x64_epi64(top_prev, 216);
    topPermutation = _mm256_permute4x64_epi64(top, 216);
    topNextPermutation = _mm256_permute4x64_epi64(top_next, 216);

    topLeftMidFirstHalflo = _mm256_add_epi16(_mm256_unpacklo_epi8(topPermutation, _mm256_setzero_si256()),
        _mm256_unpacklo_epi8(topPrevPermutation, _mm256_setzero_si256()));
    topLeftMidFirstHalflo = _mm256_add_epi16(_mm256_unpacklo_epi8(topNextPermutation, _mm256_setzero_si256()), topLeftMidFirstHalflo);

    topLeftMidFirstHalfhi = _mm256_add_epi16(_mm256_unpackhi_epi8(topPermutation, _mm256_setzero_si256()),
        _mm256_unpackhi_epi8(topPrevPermutation, _mm256_setzero_si256()));
    topLeftMidFirstHalfhi = _mm256_add_epi16(_mm256_unpackhi_epi8(topNextPermutation, _mm256_setzero_si256()), topLeftMidFirstHalfhi);

    bottomPrevPermutation = _mm256_permute4x64_epi64(bottom_prev, 216);
    bottomPermutation = _mm256_permute4x64_epi64(bottom, 216);
    bottomNextPermutation = _mm256_permute4x64_epi64(bottom_next, 216);

    bottomLeftMidFirstHalflo = _mm256_add_epi16(_mm256_unpacklo_epi8(bottomPermutation, _mm256_setzero_si256()),
        _mm256_unpacklo_epi8(bottomPrevPermutation, _mm256_setzero_si256()));
    bottomLeftMidFirstHalflo = _mm256_add_epi16(_mm256_unpacklo_epi8(bottomNextPermutation, _mm256_setzero_si256()), bottomLeftMidFirstHalflo);

    bottomLeftMidFirstHalfhi = _mm256_add_epi16(_mm256_unpackhi_epi8(bottomPermutation, _mm256_setzero_si256()),
        _mm256_unpackhi_epi8(bottomPrevPermutation, _mm256_setzero_si256()));
    bottomLeftMidFirstHalfhi = _mm256_add_epi16(_mm256_unpackhi_epi8(bottomNextPermutation, _mm256_setzero_si256()), bottomLeftMidFirstHalfhi);

    currLeftMidFirstHalflo = _mm256_add_epi16(_mm256_add_epi16(currLeftMidFirstHalflo, topLeftMidFirstHalflo), bottomLeftMidFirstHalflo);
    currLeftMidFirstHalfhi = _mm256_add_epi16(_mm256_add_epi16(currLeftMidFirstHalfhi, topLeftMidFirstHalfhi), bottomLeftMidFirstHalfhi);

    topLeftMidFirstHalflo = _mm256_unpacklo_epi16(currLeftMidFirstHalflo, _mm256_setzero_si256());
    topLeftMidFirstHalflo = _mm256_mullo_epi32(topLeftMidFirstHalflo, _mm256_set1_epi32(7282));
    topLeftMidFirstHalflo = _mm256_srli_epi32(topLeftMidFirstHalflo, 16);
    bottomLeftMidFirstHalflo = _mm256_unpackhi_epi16(currLeftMidFirstHalflo, _mm256_setzero_si256());
    bottomLeftMidFirstHalflo = _mm256_mullo_epi32(bottomLeftMidFirstHalflo, _mm256_set1_epi32(7282));
    bottomLeftMidFirstHalflo = _mm256_srli_epi32(bottomLeftMidFirstHalflo, 16);
    currLeftMidFirstHalflo = _mm256_packus_epi32(topLeftMidFirstHalflo, bottomLeftMidFirstHalflo);

    currLeftMidFirstHalflo = _mm256_insertf128_si256(_mm256_setzero_si256(), _mm_packus_epi16(_mm256_extracti128_si256(currLeftMidFirstHalflo, 0), _mm256_extracti128_si256(currLeftMidFirstHalflo, 1)), 0);

    topLeftMidFirstHalfhi = _mm256_unpacklo_epi16(currLeftMidFirstHalfhi, _mm256_setzero_si256());
    topLeftMidFirstHalfhi = _mm256_mullo_epi32(topLeftMidFirstHalfhi, _mm256_set1_epi32(7282));
    topLeftMidFirstHalfhi = _mm256_srli_epi32(topLeftMidFirstHalfhi, 16);

    bottomLeftMidFirstHalfhi = _mm256_unpackhi_epi16(currLeftMidFirstHalfhi, _mm256_setzero_si256());
    bottomLeftMidFirstHalfhi = _mm256_mullo_epi32(bottomLeftMidFirstHalfhi, _mm256_set1_epi32(7282));
    bottomLeftMidFirstHalfhi = _mm256_srli_epi32(bottomLeftMidFirstHalfhi, 16);
    currLeftMidFirstHalfhi = _mm256_packus_epi32(topLeftMidFirstHalfhi, bottomLeftMidFirstHalfhi);

    currLeftMidFirstHalflo = _mm256_insertf128_si256(currLeftMidFirstHalflo, _mm_packus_epi16(_mm256_extracti128_si256(currLeftMidFirstHalfhi, 0), _mm256_extracti128_si256(currLeftMidFirstHalfhi, 1)), 1);
    _mm256_storeu_si256((__m256i *)(ptr_denoised), currLeftMidFirstHalflo);
}

inline void chroma_strong_128_avx2_intrin(
    __m128i                        top,
    __m128i                        curr,
    __m128i                        bottom,
    __m128i                        curr_prev,
    __m128i                        curr_next,
    __m128i                        top_prev,
    __m128i                        top_next,
    __m128i                        bottom_prev,
    __m128i                        bottom_next,
    uint8_t                       *ptr_denoised
)
{
    __m128i   currLeftMidFirstHalflo, currLeftMidFirstHalfhi, topLeftMidFirstHalflo, topLeftMidFirstHalfhi,
              bottomLeftMidFirstHalflo, bottomLeftMidFirstHalfhi;

    currLeftMidFirstHalflo = _mm_add_epi16(_mm_unpacklo_epi8(curr, _mm_setzero_si128()),
        _mm_unpacklo_epi8(curr_prev, _mm_setzero_si128()));
    currLeftMidFirstHalflo = _mm_add_epi16(_mm_unpacklo_epi8(curr_next, _mm_setzero_si128()), currLeftMidFirstHalflo);

    currLeftMidFirstHalfhi = _mm_add_epi16(_mm_unpackhi_epi8(curr, _mm_setzero_si128()),
        _mm_unpackhi_epi8(curr_prev, _mm_setzero_si128()));
    currLeftMidFirstHalfhi = _mm_add_epi16(_mm_unpackhi_epi8(curr_next, _mm_setzero_si128()), currLeftMidFirstHalfhi);

    topLeftMidFirstHalflo = _mm_add_epi16(_mm_unpacklo_epi8(top, _mm_setzero_si128()),
        _mm_unpacklo_epi8(top_prev, _mm_setzero_si128()));
    topLeftMidFirstHalflo = _mm_add_epi16(_mm_unpacklo_epi8(top_next, _mm_setzero_si128()), topLeftMidFirstHalflo);

    topLeftMidFirstHalfhi = _mm_add_epi16(_mm_unpackhi_epi8(top, _mm_setzero_si128()),
        _mm_unpackhi_epi8(top_prev, _mm_setzero_si128()));
    topLeftMidFirstHalfhi = _mm_add_epi16(_mm_unpackhi_epi8(top_next, _mm_setzero_si128()), topLeftMidFirstHalfhi);

    bottomLeftMidFirstHalflo = _mm_add_epi16(_mm_unpacklo_epi8(bottom, _mm_setzero_si128()),
        _mm_unpacklo_epi8(bottom_prev, _mm_setzero_si128()));
    bottomLeftMidFirstHalflo = _mm_add_epi16(_mm_unpacklo_epi8(bottom_next, _mm_setzero_si128()), bottomLeftMidFirstHalflo);

    bottomLeftMidFirstHalfhi = _mm_add_epi16(_mm_unpackhi_epi8(bottom, _mm_setzero_si128()),
        _mm_unpackhi_epi8(bottom_prev, _mm_setzero_si128()));
    bottomLeftMidFirstHalfhi = _mm_add_epi16(_mm_unpackhi_epi8(bottom_next, _mm_setzero_si128()), bottomLeftMidFirstHalfhi);

    currLeftMidFirstHalflo = _mm_add_epi16(_mm_add_epi16(currLeftMidFirstHalflo, topLeftMidFirstHalflo), bottomLeftMidFirstHalflo);
    currLeftMidFirstHalfhi = _mm_add_epi16(_mm_add_epi16(currLeftMidFirstHalfhi, topLeftMidFirstHalfhi), bottomLeftMidFirstHalfhi);

    topLeftMidFirstHalflo = _mm_unpacklo_epi16(currLeftMidFirstHalflo, _mm_setzero_si128());
    topLeftMidFirstHalflo = _mm_mullo_epi32(topLeftMidFirstHalflo, _mm_set1_epi32(7282));
    topLeftMidFirstHalflo = _mm_srli_epi32(topLeftMidFirstHalflo, 16);
    bottomLeftMidFirstHalflo = _mm_unpackhi_epi16(currLeftMidFirstHalflo, _mm_setzero_si128());
    bottomLeftMidFirstHalflo = _mm_mullo_epi32(bottomLeftMidFirstHalflo, _mm_set1_epi32(7282));
    bottomLeftMidFirstHalflo = _mm_srli_epi32(bottomLeftMidFirstHalflo, 16);
    currLeftMidFirstHalflo = _mm_packus_epi32(topLeftMidFirstHalflo, bottomLeftMidFirstHalflo);

    currLeftMidFirstHalflo = _mm_packus_epi16(currLeftMidFirstHalflo, currLeftMidFirstHalflo);

    _mm_storel_epi64((__m128i*)(ptr_denoised), currLeftMidFirstHalflo);
}
/*******************************************
* noise_extract_luma_weak
*  weak filter Luma and store noise.
*******************************************/
void noise_extract_luma_weak_avx2_intrin(
    EbPictureBufferDesc       *input_picture_ptr,
    EbPictureBufferDesc       *denoised_picture_ptr,
    EbPictureBufferDesc       *noise_picture_ptr,
    uint32_t                       sb_origin_y,
    uint32_t                         sb_origin_x
)
{
    uint32_t  ii, jj, kk;
    uint32_t  picHeight, sb_height;
    uint32_t  picWidth;
    uint32_t  inputOriginIndex;
    uint32_t  inputOriginIndexPad;
    uint32_t  noiseOriginIndex;

    uint8_t *ptrIn;
    uint32_t stride_in;
    uint8_t *ptr_denoised, *ptrDenoisedInterm;

    uint8_t *ptr_noise, *ptrNoiseInterm;
    uint32_t strideOut;

    __m256i top, curr, bottom, curr_prev, curr_next,
        secondtop, secondcurr, secondbottom, secondcurrPrev, secondcurrNext;
    __m128i top_128, curr_128, bottom_128, curr_prev_128, curr_next_128,
        secondtop_128, secondcurr_128, secondbottom_128, secondcurrPrev_128, secondcurrNext_128;
    uint32_t idx = (sb_origin_x + BLOCK_SIZE_64 > input_picture_ptr->width) ? sb_origin_x : 0;

    //Luma
    {
        picHeight = input_picture_ptr->height;
        picWidth = input_picture_ptr->width;
        sb_height = MIN(BLOCK_SIZE_64, picHeight - sb_origin_y);
        sb_height = ((sb_origin_y + BLOCK_SIZE_64 >= picHeight) || (sb_origin_y == 0)) ? sb_height - 1 : sb_height;
        stride_in = input_picture_ptr->stride_y;
        inputOriginIndex = input_picture_ptr->origin_x + (input_picture_ptr->origin_y + sb_origin_y) * input_picture_ptr->stride_y;
        ptrIn = &(input_picture_ptr->buffer_y[inputOriginIndex]);

        inputOriginIndexPad = denoised_picture_ptr->origin_x + (denoised_picture_ptr->origin_y + sb_origin_y) * denoised_picture_ptr->stride_y;
        strideOut = denoised_picture_ptr->stride_y;
        ptr_denoised = &(denoised_picture_ptr->buffer_y[inputOriginIndexPad]);
        ptrDenoisedInterm = ptr_denoised;

        noiseOriginIndex = noise_picture_ptr->origin_x + noise_picture_ptr->origin_y * noise_picture_ptr->stride_y;
        ptr_noise = &(noise_picture_ptr->buffer_y[noiseOriginIndex]);
        ptrNoiseInterm = ptr_noise;

        ////Luma
        //a = (p[1] +
        //    p[0 + stride] + 4 * p[1 + stride] + p[2 + stride] +
        //    p[1 + 2 * stride]) / 8;

        top = curr = secondtop = secondcurr = _mm256_setzero_si256();

        for (kk = idx; kk + BLOCK_SIZE_64 <= picWidth; kk += BLOCK_SIZE_64)
        {
            for (jj = 0; jj < sb_height; jj++)
            {
                if (sb_origin_y == 0)
                {
                    if (jj == 0)
                    {
                        top = _mm256_loadu_si256((__m256i*)(ptrIn + kk + jj * stride_in));
                        secondtop = _mm256_loadu_si256((__m256i*)(ptrIn + kk + 32 + jj * stride_in));
                        curr = _mm256_loadu_si256((__m256i*)(ptrIn + kk + (1 + jj)*stride_in));
                        secondcurr = _mm256_loadu_si256((__m256i*)((ptrIn + kk + 32) + (1 + jj)*stride_in));
                        _mm256_storeu_si256((__m256i *)(ptr_denoised + kk), top);
                        _mm256_storeu_si256((__m256i *)(ptr_denoised + kk + 32), secondtop);
                        _mm256_storeu_si256((__m256i *)(ptr_noise + kk), _mm256_setzero_si256());
                        _mm256_storeu_si256((__m256i *)(ptr_noise + kk + 32), _mm256_setzero_si256());
                    }
                    curr_prev = _mm256_loadu_si256((__m256i*)(ptrIn - 1 + kk + ((1 + jj)*stride_in)));
                    curr_next = _mm256_loadu_si256((__m256i*)(ptrIn + 1 + kk + ((1 + jj)*stride_in)));
                    secondcurrPrev = _mm256_loadu_si256((__m256i*)((ptrIn + kk + 32) - 1 + ((1 + jj)*stride_in)));
                    secondcurrNext = _mm256_loadu_si256((__m256i*)((ptrIn + kk + 32) + 1 + ((1 + jj)*stride_in)));
                    bottom = _mm256_loadu_si256((__m256i*)((ptrIn + kk) + (2 + jj)* stride_in));
                    secondbottom = _mm256_loadu_si256((__m256i*)((ptrIn + kk + 32) + (2 + jj)* stride_in));
                    ptrDenoisedInterm = ptr_denoised + kk + ((1 + jj)*strideOut);
                    ptrNoiseInterm = ptr_noise + kk + ((1 + jj)*strideOut);
                }
                else
                {
                    if (jj == 0)
                    {
                        top = _mm256_loadu_si256((__m256i*)(ptrIn + kk + jj * stride_in - stride_in));
                        secondtop = _mm256_loadu_si256((__m256i*)(ptrIn + kk + 32 + jj * stride_in - stride_in));
                        curr = _mm256_loadu_si256((__m256i*)(ptrIn + kk + (1 + jj)*stride_in - stride_in));
                        secondcurr = _mm256_loadu_si256((__m256i*)((ptrIn + kk + 32) + (1 + jj)*stride_in - stride_in));
                    }
                    curr_prev = _mm256_loadu_si256((__m256i*)(ptrIn - 1 + kk + ((1 + jj)*stride_in - stride_in)));
                    curr_next = _mm256_loadu_si256((__m256i*)(ptrIn + 1 + kk + ((1 + jj)*stride_in - stride_in)));
                    secondcurrPrev = _mm256_loadu_si256((__m256i*)((ptrIn + kk + 32) - 1 + ((1 + jj)*stride_in - stride_in)));
                    secondcurrNext = _mm256_loadu_si256((__m256i*)((ptrIn + kk + 32) + 1 + ((1 + jj)*stride_in - stride_in)));
                    bottom = _mm256_loadu_si256((__m256i*)((ptrIn + kk) + (2 + jj)* stride_in - stride_in));
                    secondbottom = _mm256_loadu_si256((__m256i*)((ptrIn + kk + 32) + (2 + jj)* stride_in - stride_in));
                    ptrDenoisedInterm = ptr_denoised + kk + ((1 + jj)*strideOut - strideOut);
                    ptrNoiseInterm = ptr_noise + kk + jj * strideOut;
                }

                luma_weak_filter_avx2_intrin(
                    top,
                    curr,
                    bottom,
                    curr_prev,
                    curr_next,
                    ptrDenoisedInterm,
                    ptrNoiseInterm);

                luma_weak_filter_avx2_intrin(
                    secondtop,
                    secondcurr,
                    secondbottom,
                    secondcurrPrev,
                    secondcurrNext,
                    ptrDenoisedInterm + 32,
                    ptrNoiseInterm + 32);

                top = curr;
                curr = bottom;
                secondtop = secondcurr;
                secondcurr = secondbottom;
            }
        }

        for (; kk + 16 <= picWidth; kk += 16)
        {
            for (jj = 0; jj < sb_height; jj++)
            {
                if (sb_origin_y == 0)
                {
                    if (jj == 0)
                    {
                        top_128 = _mm_loadl_epi64((__m128i*)(ptrIn + kk + jj * stride_in));
                        secondtop_128 = _mm_loadl_epi64((__m128i*)(ptrIn + kk + 8 + jj * stride_in));
                        curr_128 = _mm_loadl_epi64((__m128i*)(ptrIn + kk + (1 + jj) * stride_in));
                        secondcurr_128 = _mm_loadl_epi64((__m128i*)((ptrIn + kk + 8) + (1 + jj) * stride_in));
                        _mm_storel_epi64((__m128i*)(ptr_denoised + kk), top_128);
                        _mm_storel_epi64((__m128i*)(ptr_denoised + kk + 8), secondtop_128);
                        _mm_storel_epi64((__m128i*)(ptr_noise + kk), _mm_setzero_si128());
                        _mm_storel_epi64((__m128i*)(ptr_noise + kk + 8), _mm_setzero_si128());
                    }
                    curr_prev_128 = _mm_loadl_epi64((__m128i*)(ptrIn - 1 + kk + ((1 + jj) * stride_in)));
                    curr_next_128 = _mm_loadl_epi64((__m128i*)(ptrIn + 1 + kk + ((1 + jj) * stride_in)));
                    secondcurrPrev_128 = _mm_loadl_epi64((__m128i*)((ptrIn + kk + 8) - 1 + ((1 + jj) * stride_in)));
                    secondcurrNext_128 = _mm_loadl_epi64((__m128i*)((ptrIn + kk + 8) + 1 + ((1 + jj) * stride_in)));
                    bottom_128 = _mm_loadl_epi64((__m128i*)((ptrIn + kk) + (2 + jj) * stride_in));
                    secondbottom_128 = _mm_loadl_epi64((__m128i*)((ptrIn + kk + 8) + (2 + jj) * stride_in));
                    ptrDenoisedInterm = ptr_denoised + kk + ((1 + jj) * strideOut);
                    ptrNoiseInterm = ptr_noise + kk + ((1 + jj) * strideOut);
                }
                else
                {
                    if (jj == 0)
                    {
                        top_128 = _mm_loadl_epi64((__m128i*)(ptrIn + kk + jj * stride_in - stride_in));
                        secondtop_128 = _mm_loadl_epi64((__m128i*)(ptrIn + kk + 8 + jj * stride_in - stride_in));
                        curr_128 = _mm_loadl_epi64((__m128i*)(ptrIn + kk + (1 + jj) * stride_in - stride_in));
                        secondcurr_128 = _mm_loadl_epi64((__m128i*)((ptrIn + kk + 8) + (1 + jj) * stride_in - stride_in));
                    }
                    curr_prev_128 = _mm_loadl_epi64((__m128i*)(ptrIn - 1 + kk + ((1 + jj) * stride_in - stride_in)));
                    curr_next_128 = _mm_loadl_epi64((__m128i*)(ptrIn + 1 + kk + ((1 + jj) * stride_in - stride_in)));
                    secondcurrPrev_128 = _mm_loadl_epi64((__m128i*)((ptrIn + kk + 8) - 1 + ((1 + jj) * stride_in - stride_in)));
                    secondcurrNext_128 = _mm_loadl_epi64((__m128i*)((ptrIn + kk + 8) + 1 + ((1 + jj) * stride_in - stride_in)));
                    bottom_128 = _mm_loadl_epi64((__m128i*)((ptrIn + kk) + (2 + jj) * stride_in - stride_in));
                    secondbottom_128 = _mm_loadl_epi64((__m128i*)((ptrIn + kk + 8) + (2 + jj) * stride_in - stride_in));
                    ptrDenoisedInterm = ptr_denoised + kk + ((1 + jj) * strideOut - strideOut);
                    ptrNoiseInterm = ptr_noise + kk + jj * strideOut;
                }

                luma_weak_filter_128_avx2_intrin(
                    top_128,
                    curr_128,
                    bottom_128,
                    curr_prev_128,
                    curr_next_128,
                    ptrDenoisedInterm,
                    ptrNoiseInterm);

                luma_weak_filter_128_avx2_intrin(
                    secondtop_128,
                    secondcurr_128,
                    secondbottom_128,
                    secondcurrPrev_128,
                    secondcurrNext_128,
                    ptrDenoisedInterm + 8,
                    ptrNoiseInterm + 8);

                top_128 = curr_128;
                curr_128 = bottom_128;
                secondtop_128 = secondcurr_128;
                secondcurr_128 = secondbottom_128;
            }
        }

        sb_height = MIN(BLOCK_SIZE_64, picHeight - sb_origin_y);

        for (jj = 0; jj < sb_height; jj++) {
            for (ii = 0; ii < picWidth; ii++) {
                if (!((jj < sb_height - 1 || sb_origin_y + sb_height < picHeight) && ii > 0 && ii < picWidth - 1)) {
                    ptr_denoised[ii + jj * strideOut] = ptrIn[ii + jj * stride_in];
                    ptr_noise[ii + jj * strideOut] = 0;
                }
            }
        }
    }
}

void noise_extract_luma_weak_lcu_avx2_intrin(
    EbPictureBufferDesc       *input_picture_ptr,
    EbPictureBufferDesc       *denoised_picture_ptr,
    EbPictureBufferDesc       *noise_picture_ptr,
    uint32_t                       sb_origin_y,
    uint32_t                         sb_origin_x
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
    uint8_t *ptr_denoised, *ptrDenoisedInterm;

    uint8_t *ptr_noise, *ptrNoiseInterm;
    uint32_t strideOut;

    __m256i top, curr, bottom, curr_prev, curr_next,
        secondtop, secondcurr, secondbottom, secondcurrPrev, secondcurrNext;
    (void)sb_origin_x;

    //Luma
    {
        picHeight = input_picture_ptr->height;
        picWidth = input_picture_ptr->width;
        sb_height = MIN(BLOCK_SIZE_64, picHeight - sb_origin_y);
        sb_width = MIN(BLOCK_SIZE_64, picWidth - sb_origin_x);
        sb_height = ((sb_origin_y + BLOCK_SIZE_64 >= picHeight) || (sb_origin_y == 0)) ? sb_height - 1 : sb_height;
        stride_in = input_picture_ptr->stride_y;
        inputOriginIndex = input_picture_ptr->origin_x + sb_origin_x + (input_picture_ptr->origin_y + sb_origin_y) * input_picture_ptr->stride_y;
        ptrIn = &(input_picture_ptr->buffer_y[inputOriginIndex]);

        inputOriginIndexPad = denoised_picture_ptr->origin_x + sb_origin_x + (denoised_picture_ptr->origin_y + sb_origin_y) * denoised_picture_ptr->stride_y;
        strideOut = denoised_picture_ptr->stride_y;
        ptr_denoised = &(denoised_picture_ptr->buffer_y[inputOriginIndexPad]);
        ptrDenoisedInterm = ptr_denoised;

        noiseOriginIndex = noise_picture_ptr->origin_x + sb_origin_x + noise_picture_ptr->origin_y * noise_picture_ptr->stride_y;
        ptr_noise = &(noise_picture_ptr->buffer_y[noiseOriginIndex]);
        ptrNoiseInterm = ptr_noise;

        ////Luma
        //a = (p[1] +
        //    p[0 + stride] + 4 * p[1 + stride] + p[2 + stride] +
        //    p[1 + 2 * stride]) / 8;

        top = curr = secondtop = secondcurr = _mm256_setzero_si256();

        //for (kk = 0; kk + BLOCK_SIZE_64 <= picWidth; kk += BLOCK_SIZE_64)
        {
            for (jj = 0; jj < sb_height; jj++)
            {
                if (sb_origin_y == 0)
                {
                    if (jj == 0)
                    {
                        top = _mm256_loadu_si256((__m256i*)(ptrIn + jj * stride_in));
                        secondtop = _mm256_loadu_si256((__m256i*)(ptrIn + 32 + jj * stride_in));
                        curr = _mm256_loadu_si256((__m256i*)(ptrIn + (1 + jj)*stride_in));
                        secondcurr = _mm256_loadu_si256((__m256i*)((ptrIn + 32) + (1 + jj)*stride_in));
                        _mm256_storeu_si256((__m256i *)(ptr_denoised), top);
                        _mm256_storeu_si256((__m256i *)(ptr_denoised + 32), secondtop);
                        _mm256_storeu_si256((__m256i *)(ptr_noise), _mm256_setzero_si256());
                        _mm256_storeu_si256((__m256i *)(ptr_noise + 32), _mm256_setzero_si256());
                    }
                    curr_prev = _mm256_loadu_si256((__m256i*)(ptrIn - 1 + ((1 + jj)*stride_in)));
                    curr_next = _mm256_loadu_si256((__m256i*)(ptrIn + 1 + ((1 + jj)*stride_in)));
                    secondcurrPrev = _mm256_loadu_si256((__m256i*)((ptrIn + 32) - 1 + ((1 + jj)*stride_in)));
                    secondcurrNext = _mm256_loadu_si256((__m256i*)((ptrIn + 32) + 1 + ((1 + jj)*stride_in)));
                    bottom = _mm256_loadu_si256((__m256i*)((ptrIn)+(2 + jj)* stride_in));
                    secondbottom = _mm256_loadu_si256((__m256i*)((ptrIn + 32) + (2 + jj)* stride_in));
                    ptrDenoisedInterm = ptr_denoised + ((1 + jj)*strideOut);
                    ptrNoiseInterm = ptr_noise + ((1 + jj)*strideOut);
                }
                else
                {
                    if (jj == 0)
                    {
                        top = _mm256_loadu_si256((__m256i*)(ptrIn + jj * stride_in - stride_in));
                        secondtop = _mm256_loadu_si256((__m256i*)(ptrIn + 32 + jj * stride_in - stride_in));
                        curr = _mm256_loadu_si256((__m256i*)(ptrIn + (1 + jj)*stride_in - stride_in));
                        secondcurr = _mm256_loadu_si256((__m256i*)((ptrIn + 32) + (1 + jj)*stride_in - stride_in));
                    }
                    curr_prev = _mm256_loadu_si256((__m256i*)(ptrIn - 1 + ((1 + jj)*stride_in - stride_in)));
                    curr_next = _mm256_loadu_si256((__m256i*)(ptrIn + 1 + ((1 + jj)*stride_in - stride_in)));
                    secondcurrPrev = _mm256_loadu_si256((__m256i*)((ptrIn + 32) - 1 + ((1 + jj)*stride_in - stride_in)));
                    secondcurrNext = _mm256_loadu_si256((__m256i*)((ptrIn + 32) + 1 + ((1 + jj)*stride_in - stride_in)));
                    bottom = _mm256_loadu_si256((__m256i*)((ptrIn)+(2 + jj)* stride_in - stride_in));
                    secondbottom = _mm256_loadu_si256((__m256i*)((ptrIn + 32) + (2 + jj)* stride_in - stride_in));
                    ptrDenoisedInterm = ptr_denoised + ((1 + jj)*strideOut - strideOut);
                    ptrNoiseInterm = ptr_noise + jj * strideOut;
                }

                luma_weak_filter_avx2_intrin(
                    top,
                    curr,
                    bottom,
                    curr_prev,
                    curr_next,
                    ptrDenoisedInterm,
                    ptrNoiseInterm);

                luma_weak_filter_avx2_intrin(
                    secondtop,
                    secondcurr,
                    secondbottom,
                    secondcurrPrev,
                    secondcurrNext,
                    ptrDenoisedInterm + 32,
                    ptrNoiseInterm + 32);

                top = curr;
                curr = bottom;
                secondtop = secondcurr;
                secondcurr = secondbottom;
            }
        }

        sb_height = MIN(BLOCK_SIZE_64, picHeight - sb_origin_y);

        for (jj = 0; jj < sb_height; jj++) {
            for (ii = 0; ii < sb_width; ii++) {
                if (!((jj > 0 || sb_origin_y > 0) && (jj < sb_height - 1 || sb_origin_y + sb_height < picHeight) && (ii > 0 || sb_origin_x > 0) && (ii + sb_origin_x) < picWidth - 1)) {
                    ptr_denoised[ii + jj * strideOut] = ptrIn[ii + jj * stride_in];
                    ptr_noise[ii + jj * strideOut] = 0;
                }
            }
        }
    }
}
/*******************************************
* noise_extract_luma_strong
*  strong filter Luma.
*******************************************/
void noise_extract_luma_strong_avx2_intrin(
    EbPictureBufferDesc       *input_picture_ptr,
    EbPictureBufferDesc       *denoised_picture_ptr,
    uint32_t                       sb_origin_y,
    uint32_t                       sb_origin_x
)
{
    uint32_t  ii, jj, kk;
    uint32_t  picHeight, sb_height;
    uint32_t  picWidth;
    uint32_t  inputOriginIndex;
    uint32_t  inputOriginIndexPad;

    uint8_t *ptrIn;
    uint32_t stride_in;
    uint8_t *ptr_denoised, *ptrDenoisedInterm;

    uint32_t strideOut;
    __m256i top, curr, bottom, curr_prev, curr_next, top_prev, top_next, bottom_prev, bottom_next,
        secondtop, secondcurr, secondcurrPrev, secondcurrNext, secondbottom, secondtopPrev, secondtopNext, secondbottomPrev, secondbottomNext;

    __m128i top_128, curr_128, bottom_128, curr_prev_128, curr_next_128, top_prev_128, top_next_128, bottom_prev_128, bottom_next_128,
        secondtop_128, secondcurr_128, secondcurrPrev_128, secondcurrNext_128, secondbottom_128, secondtopPrev_128, secondtopNext_128,
        secondbottomPrev_128, secondbottomNext_128;

    uint32_t idx = (sb_origin_x + BLOCK_SIZE_64 > input_picture_ptr->width) ? sb_origin_x : 0;
    //Luma
    {
        picHeight = input_picture_ptr->height;
        picWidth = input_picture_ptr->width;
        sb_height = MIN(BLOCK_SIZE_64, picHeight - sb_origin_y);

        sb_height = ((sb_origin_y + BLOCK_SIZE_64 >= picHeight) || (sb_origin_y == 0)) ? sb_height - 1 : sb_height;
        stride_in = input_picture_ptr->stride_y;
        inputOriginIndex = input_picture_ptr->origin_x + (input_picture_ptr->origin_y + sb_origin_y)* input_picture_ptr->stride_y;
        ptrIn = &(input_picture_ptr->buffer_y[inputOriginIndex]);

        inputOriginIndexPad = denoised_picture_ptr->origin_x + (denoised_picture_ptr->origin_y + sb_origin_y) * denoised_picture_ptr->stride_y;
        strideOut = denoised_picture_ptr->stride_y;
        ptr_denoised = &(denoised_picture_ptr->buffer_y[inputOriginIndexPad]);
        ptrDenoisedInterm = ptr_denoised;

        top = curr = secondtop = secondcurr = top_next = top_prev = curr_next = curr_prev = secondcurrPrev = secondcurrNext = secondtopPrev = secondtopNext = _mm256_setzero_si256();
        for (kk = idx; kk + BLOCK_SIZE_64 <= picWidth; kk += BLOCK_SIZE_64)
        {
            for (jj = 0; jj < sb_height; jj++)
            {
                if (sb_origin_y == 0)
                {
                    if (jj == 0)
                    {
                        top = _mm256_loadu_si256((__m256i*)(ptrIn + kk + jj * stride_in));
                        secondtop = _mm256_loadu_si256((__m256i*)(ptrIn + kk + 32 + jj * stride_in));

                        curr = _mm256_loadu_si256((__m256i*)(ptrIn + kk + (1 + jj)*stride_in));
                        secondcurr = _mm256_loadu_si256((__m256i*)((ptrIn + kk + 32) + (1 + jj)*stride_in));

                        top_prev = _mm256_loadu_si256((__m256i*)(ptrIn - 1 + kk + ((jj)*stride_in)));
                        secondtopPrev = _mm256_loadu_si256((__m256i*)(ptrIn - 1 + kk + 32 + ((jj)*stride_in)));

                        top_next = _mm256_loadu_si256((__m256i*)(ptrIn + 1 + kk + ((jj)*stride_in)));
                        secondtopNext = _mm256_loadu_si256((__m256i*)(ptrIn + 1 + kk + 32 + ((jj)*stride_in)));

                        curr_prev = _mm256_loadu_si256((__m256i*)(ptrIn - 1 + kk + ((1 + jj)*stride_in)));
                        secondcurrPrev = _mm256_loadu_si256((__m256i*)(ptrIn - 1 + kk + 32 + ((1 + jj)*stride_in)));

                        curr_next = _mm256_loadu_si256((__m256i*)(ptrIn + 1 + kk + ((1 + jj)*stride_in)));
                        secondcurrNext = _mm256_loadu_si256((__m256i*)(ptrIn + 1 + kk + 32 + ((1 + jj)*stride_in)));

                        _mm256_storeu_si256((__m256i *)(ptr_denoised + kk), top);
                        _mm256_storeu_si256((__m256i *)(ptr_denoised + kk + 32), secondtop);
                    }
                    bottom_prev = _mm256_loadu_si256((__m256i*)(ptrIn - 1 + kk + ((2 + jj)*stride_in)));
                    secondbottomPrev = _mm256_loadu_si256((__m256i*)(ptrIn - 1 + kk + 32 + ((2 + jj)*stride_in)));

                    bottom_next = _mm256_loadu_si256((__m256i*)(ptrIn + 1 + kk + ((2 + jj)*stride_in)));
                    secondbottomNext = _mm256_loadu_si256((__m256i*)(ptrIn + 1 + kk + 32 + ((2 + jj)*stride_in)));

                    bottom = _mm256_loadu_si256((__m256i*)((ptrIn + kk) + (2 + jj)* stride_in));
                    secondbottom = _mm256_loadu_si256((__m256i*)((ptrIn + kk + 32) + (2 + jj)* stride_in));
                    ptrDenoisedInterm = ptr_denoised + kk + ((1 + jj)*strideOut);
                }
                else
                {
                    if (jj == 0)
                    {
                        top = _mm256_loadu_si256((__m256i*)(ptrIn + kk + jj * stride_in - stride_in));
                        secondtop = _mm256_loadu_si256((__m256i*)(ptrIn + kk + 32 + jj * stride_in - stride_in));
                        curr = _mm256_loadu_si256((__m256i*)(ptrIn + kk + (1 + jj)*stride_in - stride_in));
                        secondcurr = _mm256_loadu_si256((__m256i*)((ptrIn + kk + 32) + (1 + jj)*stride_in - stride_in));
                        top_prev = _mm256_loadu_si256((__m256i*)(ptrIn - 1 + kk + ((jj)*stride_in) - stride_in));
                        secondtopPrev = _mm256_loadu_si256((__m256i*)(ptrIn - 1 + kk + 32 + ((jj)*stride_in) - stride_in));

                        top_next = _mm256_loadu_si256((__m256i*)(ptrIn + 1 + kk + ((jj)*stride_in) - stride_in));
                        secondtopNext = _mm256_loadu_si256((__m256i*)(ptrIn + 1 + kk + 32 + ((jj)*stride_in) - stride_in));

                        curr_prev = _mm256_loadu_si256((__m256i*)(ptrIn - 1 + kk + ((1 + jj)*stride_in - stride_in)));
                        secondcurrPrev = _mm256_loadu_si256((__m256i*)(ptrIn - 1 + kk + 32 + ((1 + jj)*stride_in - stride_in)));

                        curr_next = _mm256_loadu_si256((__m256i*)(ptrIn + 1 + kk + ((1 + jj)*stride_in - stride_in)));
                        secondcurrNext = _mm256_loadu_si256((__m256i*)(ptrIn + 1 + kk + 32 + ((1 + jj)*stride_in - stride_in)));
                    }
                    bottom_prev = _mm256_loadu_si256((__m256i*)(ptrIn - 1 + kk + ((2 + jj)*stride_in) - stride_in));
                    secondbottomPrev = _mm256_loadu_si256((__m256i*)(ptrIn - 1 + kk + 32 + ((2 + jj)*stride_in - stride_in)));

                    bottom_next = _mm256_loadu_si256((__m256i*)(ptrIn + 1 + kk + ((2 + jj)*stride_in) - stride_in));
                    secondbottomNext = _mm256_loadu_si256((__m256i*)(ptrIn + 1 + kk + 32 + ((2 + jj)*stride_in - stride_in)));

                    bottom = _mm256_loadu_si256((__m256i*)((ptrIn + kk) + (2 + jj)* stride_in - stride_in));
                    secondbottom = _mm256_loadu_si256((__m256i*)((ptrIn + kk + 32) + (2 + jj)* stride_in - stride_in));

                    ptrDenoisedInterm = ptr_denoised + kk + ((1 + jj)*strideOut - strideOut);
                }

                chroma_weak_luma_strong_filter_avx2_intrin(
                    top,
                    curr,
                    bottom,
                    curr_prev,
                    curr_next,
                    top_prev,
                    top_next,
                    bottom_prev,
                    bottom_next,
                    ptrDenoisedInterm);

                chroma_weak_luma_strong_filter_avx2_intrin(
                    secondtop,
                    secondcurr,
                    secondbottom,
                    secondcurrPrev,
                    secondcurrNext,
                    secondtopPrev,
                    secondtopNext,
                    secondbottomPrev,
                    secondbottomNext,
                    ptrDenoisedInterm + 32);

                top = curr;
                curr = bottom;
                top_prev = curr_prev;
                top_next = curr_next;
                curr_prev = bottom_prev;
                curr_next = bottom_next;
                secondtop = secondcurr;
                secondcurr = secondbottom;
                secondtopPrev = secondcurrPrev;
                secondtopNext = secondcurrNext;
                secondcurrPrev = secondbottomPrev;
                secondcurrNext = secondbottomNext;
            }
        }

        for (; kk + 16 <= picWidth; kk += 16)
        {
            for (jj = 0; jj < sb_height; jj++)
            {
                if (sb_origin_y == 0)
                {
                    if (jj == 0)
                    {
                        top_128 = _mm_loadl_epi64((__m128i*)(ptrIn + kk + jj * stride_in));
                        secondtop_128 = _mm_loadl_epi64((__m128i*)(ptrIn + kk + 8 + jj * stride_in));

                        curr_128 = _mm_loadl_epi64((__m128i*)(ptrIn + kk + (1 + jj) * stride_in));
                        secondcurr_128 = _mm_loadl_epi64((__m128i*)((ptrIn + kk + 8) + (1 + jj) * stride_in));

                        top_prev_128 = _mm_loadl_epi64((__m128i*)(ptrIn - 1 + kk + ((jj)*stride_in)));
                        secondtopPrev_128 = _mm_loadl_epi64((__m128i*)(ptrIn - 1 + kk + 8 + ((jj)*stride_in)));

                        top_next_128 = _mm_loadl_epi64((__m128i*)(ptrIn + 1 + kk + ((jj)*stride_in)));
                        secondtopNext_128 = _mm_loadl_epi64((__m128i*)(ptrIn + 1 + kk + 8 + ((jj)*stride_in)));

                        curr_prev_128 = _mm_loadl_epi64((__m128i*)(ptrIn - 1 + kk + ((1 + jj) * stride_in)));
                        secondcurrPrev_128 = _mm_loadl_epi64((__m128i*)(ptrIn - 1 + kk + 8 + ((1 + jj) * stride_in)));

                        curr_next_128 = _mm_loadl_epi64((__m128i*)(ptrIn + 1 + kk + ((1 + jj) * stride_in)));
                        secondcurrNext_128 = _mm_loadl_epi64((__m128i*)(ptrIn + 1 + kk + 8 + ((1 + jj) * stride_in)));

                        _mm_storel_epi64((__m128i*)(ptr_denoised + kk), top_128);
                        _mm_storel_epi64((__m128i*)(ptr_denoised + kk + 8), secondtop_128);
                    }
                    bottom_prev_128 = _mm_loadl_epi64((__m128i*)(ptrIn - 1 + kk + ((2 + jj) * stride_in)));
                    secondbottomPrev_128 = _mm_loadl_epi64((__m128i*)(ptrIn - 1 + kk + 8 + ((2 + jj) * stride_in)));

                    bottom_next_128 = _mm_loadl_epi64((__m128i*)(ptrIn + 1 + kk + ((2 + jj) * stride_in)));
                    secondbottomNext_128 = _mm_loadl_epi64((__m128i*)(ptrIn + 1 + kk + 8 + ((2 + jj) * stride_in)));

                    bottom_128 = _mm_loadl_epi64((__m128i*)((ptrIn + kk) + (2 + jj) * stride_in));
                    secondbottom_128 = _mm_loadl_epi64((__m128i*)((ptrIn + kk + 8) + (2 + jj) * stride_in));
                    ptrDenoisedInterm = ptr_denoised + kk + ((1 + jj) * strideOut);
                }
                else
                {
                    if (jj == 0)
                    {
                        top_128 = _mm_loadl_epi64((__m128i*)(ptrIn + kk + jj * stride_in - stride_in));
                        secondtop_128 = _mm_loadl_epi64((__m128i*)(ptrIn + kk + 8 + jj * stride_in - stride_in));
                        curr_128 = _mm_loadl_epi64((__m128i*)(ptrIn + kk + (1 + jj) * stride_in - stride_in));
                        secondcurr_128 = _mm_loadl_epi64((__m128i*)((ptrIn + kk + 8) + (1 + jj) * stride_in - stride_in));
                        top_prev_128 = _mm_loadl_epi64((__m128i*)(ptrIn - 1 + kk + ((jj)*stride_in) - stride_in));
                        secondtopPrev_128 = _mm_loadl_epi64((__m128i*)(ptrIn - 1 + kk + 8 + ((jj)*stride_in) - stride_in));

                        top_next_128 = _mm_loadl_epi64((__m128i*)(ptrIn + 1 + kk + ((jj)*stride_in) - stride_in));
                        secondtopNext_128 = _mm_loadl_epi64((__m128i*)(ptrIn + 1 + kk + 8 + ((jj)*stride_in) - stride_in));

                       curr_prev_128 = _mm_loadl_epi64((__m128i*)(ptrIn - 1 + kk + ((1 + jj) * stride_in - stride_in)));
                       secondcurrPrev_128 = _mm_loadl_epi64((__m128i*)(ptrIn - 1 + kk + 8 + ((1 + jj) * stride_in - stride_in)));

                       curr_next_128 = _mm_loadl_epi64((__m128i*)(ptrIn + 1 + kk + ((1 + jj) * stride_in - stride_in)));
                       secondcurrNext_128 = _mm_loadl_epi64((__m128i*)(ptrIn + 1 + kk + 8 + ((1 + jj) * stride_in - stride_in)));
                    }
                    bottom_prev_128 = _mm_loadl_epi64((__m128i*)(ptrIn - 1 + kk + ((2 + jj) * stride_in) - stride_in));
                    secondbottomPrev_128 = _mm_loadl_epi64((__m128i*)(ptrIn - 1 + kk + 8 + ((2 + jj) * stride_in - stride_in)));

                    bottom_next_128 = _mm_loadl_epi64((__m128i*)(ptrIn + 1 + kk + ((2 + jj) * stride_in) - stride_in));
                    secondbottomNext_128 = _mm_loadl_epi64((__m128i*)(ptrIn + 1 + kk + 8 + ((2 + jj) * stride_in - stride_in)));

                    bottom_128 = _mm_loadl_epi64((__m128i*)((ptrIn + kk) + (2 + jj) * stride_in - stride_in));
                    secondbottom_128 = _mm_loadl_epi64((__m128i*)((ptrIn + kk + 8) + (2 + jj) * stride_in - stride_in));

                    ptrDenoisedInterm = ptr_denoised + kk + ((1 + jj) * strideOut - strideOut);
                }

                chroma_weak_luma_strong_filter_128_avx2_intrin(
                    top_128,
                    curr_128,
                    bottom_128,
                    curr_prev_128,
                    curr_next_128,
                    top_prev_128,
                    top_next_128,
                    bottom_prev_128,
                    bottom_next_128,
                    ptrDenoisedInterm);

                 chroma_weak_luma_strong_filter_128_avx2_intrin(
                    secondtop_128,
                    secondcurr_128,
                    secondbottom_128,
                    secondcurrPrev_128,
                    secondcurrNext_128,
                    secondtopPrev_128,
                    secondtopNext_128,
                    secondbottomPrev_128,
                    secondbottomNext_128,
                    ptrDenoisedInterm + 8);

                top_128 = curr_128;
                curr_128 = bottom_128;
                top_prev_128 = curr_prev_128;
                top_next_128 = curr_next_128;
                curr_prev_128 = bottom_prev_128;
                curr_next_128 = bottom_next_128;
                secondtop_128 = secondcurr_128;
                secondcurr_128 = secondbottom_128;
                secondtopPrev_128 = secondcurrPrev_128;
                secondtopNext_128 = secondcurrNext_128;
                secondcurrPrev_128 = secondbottomPrev_128;
                secondcurrNext_128 = secondbottomNext_128;
            }
        }

        sb_height = MIN(BLOCK_SIZE_64, picHeight - sb_origin_y);

        for (jj = 0; jj < sb_height; jj++) {
            for (ii = 0; ii < picWidth; ii++) {
                if (!((jj < sb_height - 1 || sb_origin_y + sb_height < picHeight) && ii > 0 && ii < picWidth - 1))
                    ptr_denoised[ii + jj * strideOut] = ptrIn[ii + jj * stride_in];
            }
        }
    }
}

/*******************************************
* noise_extract_chroma_strong
*  strong filter chroma.
*******************************************/
void noise_extract_chroma_strong_avx2_intrin(
    EbPictureBufferDesc       *input_picture_ptr,
    EbPictureBufferDesc       *denoised_picture_ptr,
    uint32_t                       sb_origin_y,
    uint32_t                       sb_origin_x
)
{
    uint32_t  ii, jj, kk;
    uint32_t  picHeight, sb_height;
    uint32_t  picWidth;
    uint32_t  inputOriginIndex;
    uint32_t  inputOriginIndexPad;

    uint8_t *ptrIn, *ptrInCr;
    uint32_t stride_in, strideInCr;
    uint8_t *ptr_denoised, *ptrDenoisedInterm, *ptrDenoisedCr, *ptrDenoisedIntermCr;

    uint32_t strideOut, strideOutCr;
    __m256i top, curr, bottom, curr_prev, curr_next, top_prev, top_next, bottom_prev, bottom_next,
        topCr, currCr, bottomCr, currPrevCr, currNextCr, topPrevCr, topNextCr, bottomPrevCr, bottomNextCr;

    __m128i top_128, curr_128, bottom_128, curr_prev_128, curr_next_128, top_prev_128, top_next_128,
        bottom_prev_128, bottom_next_128, topCr_128, currCr_128, bottomCr_128, currPrevCr_128,
        currNextCr_128, topPrevCr_128, topNextCr_128, bottomPrevCr_128, bottomNextCr_128;
    uint32_t idx = (sb_origin_x + BLOCK_SIZE_64 > input_picture_ptr->width) ? sb_origin_x : 0;
    {
        picHeight = input_picture_ptr->height / 2;
        picWidth = input_picture_ptr->width / 2;
        sb_height = MIN(BLOCK_SIZE_64 / 2, picHeight - sb_origin_y);

        sb_height = ((sb_origin_y + BLOCK_SIZE_64 / 2 >= picHeight) || (sb_origin_y == 0)) ? sb_height - 1 : sb_height;

        stride_in = input_picture_ptr->stride_cb;
        inputOriginIndex = input_picture_ptr->origin_x / 2 + (input_picture_ptr->origin_y / 2 + sb_origin_y)  * input_picture_ptr->stride_cb;
        ptrIn = &(input_picture_ptr->buffer_cb[inputOriginIndex]);

        inputOriginIndexPad = denoised_picture_ptr->origin_x / 2 + (denoised_picture_ptr->origin_y / 2 + sb_origin_y)  * denoised_picture_ptr->stride_cb;
        strideOut = denoised_picture_ptr->stride_cb;
        ptr_denoised = &(denoised_picture_ptr->buffer_cb[inputOriginIndexPad]);
        ptrDenoisedInterm = ptr_denoised;

        strideInCr = input_picture_ptr->stride_cr;
        inputOriginIndex = input_picture_ptr->origin_x / 2 + (input_picture_ptr->origin_y / 2 + sb_origin_y)  * input_picture_ptr->stride_cr;
        ptrInCr = &(input_picture_ptr->buffer_cr[inputOriginIndex]);

        inputOriginIndexPad = denoised_picture_ptr->origin_x / 2 + (denoised_picture_ptr->origin_y / 2 + sb_origin_y)  * denoised_picture_ptr->stride_cr;
        strideOutCr = denoised_picture_ptr->stride_cr;
        ptrDenoisedCr = &(denoised_picture_ptr->buffer_cr[inputOriginIndexPad]);
        ptrDenoisedIntermCr = ptrDenoisedCr;
        ////Chroma
        //a = (4 * p[0] + 4 * p[1] + 4 * p[2] +
        //    4 * p[0 + stride] + 4 * p[1 + stride] + 4 * p[2 + stride] +
        //    4 * p[0 + 2 * stride] + 4 * p[1 + 2 * stride] + 4 * p[2 + 2 * stride]) / 36;

        top = curr = top_next = top_prev = curr_next = curr_prev = topCr = currCr = topNextCr = topPrevCr = currNextCr = currPrevCr = _mm256_setzero_si256();

        for (kk = idx; kk + BLOCK_SIZE_64 / 2 <= picWidth; kk += BLOCK_SIZE_64 / 2)
        {
            for (jj = 0; jj < sb_height; jj++)
            {
                if (sb_origin_y == 0)
                {
                    if (jj == 0)
                    {
                        top = _mm256_loadu_si256((__m256i*)(ptrIn + kk + jj * stride_in));
                        curr = _mm256_loadu_si256((__m256i*)(ptrIn + kk + (1 + jj)*stride_in));
                        top_prev = _mm256_loadu_si256((__m256i*)(ptrIn - 1 + kk + ((jj)*stride_in)));
                        top_next = _mm256_loadu_si256((__m256i*)(ptrIn + 1 + kk + ((jj)*stride_in)));
                        curr_prev = _mm256_loadu_si256((__m256i*)(ptrIn - 1 + kk + ((1 + jj)*stride_in)));
                        curr_next = _mm256_loadu_si256((__m256i*)(ptrIn + 1 + kk + ((1 + jj)*stride_in)));
                        topCr = _mm256_loadu_si256((__m256i*)(ptrInCr + kk + jj * strideInCr));
                        currCr = _mm256_loadu_si256((__m256i*)(ptrInCr + kk + (1 + jj)*strideInCr));
                        topPrevCr = _mm256_loadu_si256((__m256i*)(ptrInCr - 1 + kk + ((jj)*strideInCr)));
                        topNextCr = _mm256_loadu_si256((__m256i*)(ptrInCr + 1 + kk + ((jj)*strideInCr)));
                        currPrevCr = _mm256_loadu_si256((__m256i*)(ptrInCr - 1 + kk + ((1 + jj)*strideInCr)));
                        currNextCr = _mm256_loadu_si256((__m256i*)(ptrInCr + 1 + kk + ((1 + jj)*strideInCr)));
                        _mm256_storeu_si256((__m256i *)(ptr_denoised + kk), top);
                        _mm256_storeu_si256((__m256i *)(ptrDenoisedCr + kk), topCr);
                    }
                    bottom_prev = _mm256_loadu_si256((__m256i*)(ptrIn - 1 + kk + ((2 + jj)*stride_in)));
                    bottom_next = _mm256_loadu_si256((__m256i*)(ptrIn + 1 + kk + ((2 + jj)*stride_in)));
                    bottom = _mm256_loadu_si256((__m256i*)((ptrIn + kk) + (2 + jj)* stride_in));
                    bottomPrevCr = _mm256_loadu_si256((__m256i*)(ptrInCr - 1 + kk + ((2 + jj)*strideInCr)));
                    bottomNextCr = _mm256_loadu_si256((__m256i*)(ptrInCr + 1 + kk + ((2 + jj)*strideInCr)));
                    bottomCr = _mm256_loadu_si256((__m256i*)((ptrInCr + kk) + (2 + jj)* strideInCr));
                    ptrDenoisedInterm = ptr_denoised + kk + ((1 + jj)*strideOut);
                    ptrDenoisedIntermCr = ptrDenoisedCr + kk + ((1 + jj)*strideOutCr);
                }
                else
                {
                    if (jj == 0)
                    {
                        top = _mm256_loadu_si256((__m256i*)(ptrIn + kk + jj * stride_in - stride_in));
                        curr = _mm256_loadu_si256((__m256i*)(ptrIn + kk + (1 + jj)*stride_in - stride_in));
                        top_prev = _mm256_loadu_si256((__m256i*)(ptrIn - 1 + kk + ((jj)*stride_in) - stride_in));
                        top_next = _mm256_loadu_si256((__m256i*)(ptrIn + 1 + kk + ((jj)*stride_in) - stride_in));
                        curr_prev = _mm256_loadu_si256((__m256i*)(ptrIn - 1 + kk + ((1 + jj)*stride_in - stride_in)));
                        curr_next = _mm256_loadu_si256((__m256i*)(ptrIn + 1 + kk + ((1 + jj)*stride_in - stride_in)));
                        topCr = _mm256_loadu_si256((__m256i*)(ptrInCr + kk + jj * strideInCr - strideInCr));
                        currCr = _mm256_loadu_si256((__m256i*)(ptrInCr + kk + (1 + jj)*strideInCr - strideInCr));
                        topPrevCr = _mm256_loadu_si256((__m256i*)(ptrInCr - 1 + kk + ((jj)*strideInCr) - strideInCr));
                        topNextCr = _mm256_loadu_si256((__m256i*)(ptrInCr + 1 + kk + ((jj)*strideInCr) - strideInCr));
                        currPrevCr = _mm256_loadu_si256((__m256i*)(ptrInCr - 1 + kk + ((1 + jj)*strideInCr - strideInCr)));
                        currNextCr = _mm256_loadu_si256((__m256i*)(ptrInCr + 1 + kk + ((1 + jj)*strideInCr - strideInCr)));
                    }
                    bottom_prev = _mm256_loadu_si256((__m256i*)(ptrIn - 1 + kk + ((2 + jj)*stride_in) - stride_in));
                    bottom_next = _mm256_loadu_si256((__m256i*)(ptrIn + 1 + kk + ((2 + jj)*stride_in) - stride_in));
                    bottom = _mm256_loadu_si256((__m256i*)((ptrIn + kk) + (2 + jj)* stride_in - stride_in));
                    ptrDenoisedInterm = ptr_denoised + kk + ((1 + jj)*strideOut - strideOut);
                    bottomPrevCr = _mm256_loadu_si256((__m256i*)(ptrInCr - 1 + kk + ((2 + jj)*strideInCr) - strideInCr));
                    bottomNextCr = _mm256_loadu_si256((__m256i*)(ptrInCr + 1 + kk + ((2 + jj)*strideInCr) - strideInCr));
                    bottomCr = _mm256_loadu_si256((__m256i*)((ptrInCr + kk) + (2 + jj)* strideInCr - strideInCr));
                    ptrDenoisedIntermCr = ptrDenoisedCr + kk + ((1 + jj)*strideOutCr - strideOutCr);
                }

                chroma_strong_avx2_intrin(
                    top,
                    curr,
                    bottom,
                    curr_prev,
                    curr_next,
                    top_prev,
                    top_next,
                    bottom_prev,
                    bottom_next,
                    ptrDenoisedInterm);

                chroma_strong_avx2_intrin(
                    topCr,
                    currCr,
                    bottomCr,
                    currPrevCr,
                    currNextCr,
                    topPrevCr,
                    topNextCr,
                    bottomPrevCr,
                    bottomNextCr,
                    ptrDenoisedIntermCr);

                top = curr;
                curr = bottom;
                top_prev = curr_prev;
                top_next = curr_next;
                curr_prev = bottom_prev;
                curr_next = bottom_next;
                topCr = currCr;
                currCr = bottomCr;
                topPrevCr = currPrevCr;
                topNextCr = currNextCr;
                currPrevCr = bottomPrevCr;
                currNextCr = bottomNextCr;
            }
        }

        for (; kk + 8 <= picWidth; kk += 8)
        {
            for (jj = 0; jj < sb_height; jj++)
            {
                if (sb_origin_y == 0)
                {
                    if (jj == 0)
                    {
                        top_128 = _mm_loadl_epi64((__m128i*)(ptrIn + kk + jj * stride_in));
                        curr_128 = _mm_loadl_epi64((__m128i*)(ptrIn + kk + (1 + jj) * stride_in));
                        top_prev_128 = _mm_loadl_epi64((__m128i*)(ptrIn - 1 + kk + ((jj)*stride_in)));
                        top_next_128 = _mm_loadl_epi64((__m128i*)(ptrIn + 1 + kk + ((jj)*stride_in)));
                        curr_prev_128 = _mm_loadl_epi64((__m128i*)(ptrIn - 1 + kk + ((1 + jj) * stride_in)));
                        curr_next_128 = _mm_loadl_epi64((__m128i*)(ptrIn + 1 + kk + ((1 + jj) * stride_in)));
                        topCr_128 = _mm_loadl_epi64((__m128i*)(ptrInCr + kk + jj * strideInCr));
                        currCr_128 = _mm_loadl_epi64((__m128i*)(ptrInCr + kk + (1 + jj) * strideInCr));
                        topPrevCr_128 = _mm_loadl_epi64((__m128i*)(ptrInCr - 1 + kk + ((jj)*strideInCr)));
                        topNextCr_128 = _mm_loadl_epi64((__m128i*)(ptrInCr + 1 + kk + ((jj)*strideInCr)));
                        currPrevCr_128 = _mm_loadl_epi64((__m128i*)(ptrInCr - 1 + kk + ((1 + jj) * strideInCr)));
                        currNextCr_128 = _mm_loadl_epi64((__m128i*)(ptrInCr + 1 + kk + ((1 + jj) * strideInCr)));
                        _mm_storel_epi64((__m128i*)(ptr_denoised + kk), top_128);
                        _mm_storel_epi64((__m128i*)(ptrDenoisedCr + kk), topCr_128);
                    }
                    bottom_prev_128 = _mm_loadl_epi64((__m128i*)(ptrIn - 1 + kk + ((2 + jj) * stride_in)));
                    bottom_next_128 = _mm_loadl_epi64((__m128i*)(ptrIn + 1 + kk + ((2 + jj) * stride_in)));
                    bottom_128 = _mm_loadl_epi64((__m128i*)((ptrIn + kk) + (2 + jj) * stride_in));
                    bottomPrevCr_128 = _mm_loadl_epi64((__m128i*)(ptrInCr - 1 + kk + ((2 + jj) * strideInCr)));
                    bottomNextCr_128 = _mm_loadl_epi64((__m128i*)(ptrInCr + 1 + kk + ((2 + jj) * strideInCr)));
                    bottomCr_128 = _mm_loadl_epi64((__m128i*)((ptrInCr + kk) + (2 + jj) * strideInCr));
                    ptrDenoisedInterm = ptr_denoised + kk + ((1 + jj) * strideOut);
                    ptrDenoisedIntermCr = ptrDenoisedCr + kk + ((1 + jj) * strideOutCr);
                }
                else
                {
                    if (jj == 0)
                    {
                        top_128 = _mm_loadl_epi64((__m128i*)(ptrIn + kk + jj * stride_in - stride_in));
                        curr_128 = _mm_loadl_epi64((__m128i*)(ptrIn + kk + (1 + jj) * stride_in - stride_in));
                        top_prev_128 = _mm_loadl_epi64((__m128i*)(ptrIn - 1 + kk + ((jj)*stride_in) - stride_in));
                        top_next_128 = _mm_loadl_epi64((__m128i*)(ptrIn + 1 + kk + ((jj)*stride_in) - stride_in));
                        curr_prev_128 = _mm_loadl_epi64((__m128i*)(ptrIn - 1 + kk + ((1 + jj) * stride_in - stride_in)));
                        curr_next_128 = _mm_loadl_epi64((__m128i*)(ptrIn + 1 + kk + ((1 + jj) * stride_in - stride_in)));
                        topCr_128 = _mm_loadl_epi64((__m128i*)(ptrInCr + kk + jj * strideInCr - strideInCr));
                        currCr_128 = _mm_loadl_epi64((__m128i*)(ptrInCr + kk + (1 + jj) * strideInCr - strideInCr));
                        topPrevCr_128 = _mm_loadl_epi64((__m128i*)(ptrInCr - 1 + kk + ((jj)*strideInCr) - strideInCr));
                        topNextCr_128 = _mm_loadl_epi64((__m128i*)(ptrInCr + 1 + kk + ((jj)*strideInCr) - strideInCr));
                        currPrevCr_128 = _mm_loadl_epi64((__m128i*)(ptrInCr - 1 + kk + ((1 + jj) * strideInCr - strideInCr)));
                        currNextCr_128 = _mm_loadl_epi64((__m128i*)(ptrInCr + 1 + kk + ((1 + jj) * strideInCr - strideInCr)));
                    }
                    bottom_prev_128 = _mm_loadl_epi64((__m128i*)(ptrIn - 1 + kk + ((2 + jj) * stride_in) - stride_in));
                    bottom_next_128 = _mm_loadl_epi64((__m128i*)(ptrIn + 1 + kk + ((2 + jj) * stride_in) - stride_in));
                    bottom_128 = _mm_loadl_epi64((__m128i*)((ptrIn + kk) + (2 + jj) * stride_in - stride_in));
                    ptrDenoisedInterm = ptr_denoised + kk + ((1 + jj) * strideOut - strideOut);
                    bottomPrevCr_128 = _mm_loadl_epi64((__m128i*)(ptrInCr - 1 + kk + ((2 + jj) * strideInCr) - strideInCr));
                    bottomNextCr_128 = _mm_loadl_epi64((__m128i*)(ptrInCr + 1 + kk + ((2 + jj) * strideInCr) - strideInCr));
                    bottomCr_128 = _mm_loadl_epi64((__m128i*)((ptrInCr + kk) + (2 + jj) * strideInCr - strideInCr));
                    ptrDenoisedIntermCr = ptrDenoisedCr + kk + ((1 + jj) * strideOutCr - strideOutCr);
                }

                chroma_strong_128_avx2_intrin(
                    top_128,
                    curr_128,
                    bottom_128,
                    curr_prev_128,
                    curr_next_128,
                    top_prev_128,
                    top_next_128,
                    bottom_prev_128,
                    bottom_next_128,
                    ptrDenoisedInterm);

                chroma_strong_128_avx2_intrin(
                    topCr_128,
                    currCr_128,
                    bottomCr_128,
                    currPrevCr_128,
                    currNextCr_128,
                    topPrevCr_128,
                    topNextCr_128,
                    bottomPrevCr_128,
                    bottomNextCr_128,
                    ptrDenoisedIntermCr);

                top_128 = curr_128;
                curr_128 = bottom_128;
                top_prev_128 = curr_prev_128;
                top_next_128 = curr_next_128;
                curr_prev_128 = bottom_prev_128;
                curr_next_128 = bottom_next_128;
                topCr_128 = currCr_128;
                currCr_128 = bottomCr_128;
                topPrevCr_128 = currPrevCr_128;
                topNextCr_128 = currNextCr_128;
                currPrevCr_128 = bottomPrevCr_128;
                currNextCr_128 = bottomNextCr_128;
            }
        }

        sb_height = MIN(BLOCK_SIZE_64 / 2, picHeight - sb_origin_y);

        for (jj = 0; jj < sb_height; jj++) {
            for (ii = idx; ii < picWidth; ii++) {
                if (!((jj < sb_height - 1 || (sb_origin_y + sb_height) < picHeight) && ii > 0 && ii < picWidth - 1)) {
                    ptr_denoised[ii + jj * strideOut] = ptrIn[ii + jj * stride_in];
                    ptrDenoisedCr[ii + jj * strideOut] = ptrInCr[ii + jj * stride_in];
                }
            }
        }
    }
}

/*******************************************
* noise_extract_chroma_weak
*  weak filter chroma.
*******************************************/
void noise_extract_chroma_weak_avx2_intrin(
    EbPictureBufferDesc       *input_picture_ptr,
    EbPictureBufferDesc       *denoised_picture_ptr,
    uint32_t                       sb_origin_y,
    uint32_t                       sb_origin_x
)
{
    uint32_t  ii, jj, kk;
    uint32_t  picHeight, sb_height;
    uint32_t  picWidth;
    uint32_t  inputOriginIndex;
    uint32_t  inputOriginIndexPad;

    uint8_t *ptrIn, *ptrInCr;
    uint32_t stride_in, strideInCr;
    uint8_t *ptr_denoised, *ptrDenoisedInterm, *ptrDenoisedCr, *ptrDenoisedIntermCr;

    uint32_t strideOut, strideOutCr;

    __m256i top, curr, bottom, curr_prev, curr_next, top_prev, top_next, bottom_prev, bottom_next,
        topCr, currCr, bottomCr, currPrevCr, currNextCr, topPrevCr, topNextCr, bottomPrevCr, bottomNextCr;

    __m128i top_128, curr_128, bottom_128, curr_prev_128, curr_next_128, top_prev_128, top_next_128,
        bottom_prev_128, bottom_next_128, topCr_128, currCr_128, bottomCr_128, currPrevCr_128,
        currNextCr_128, topPrevCr_128, topNextCr_128, bottomPrevCr_128, bottomNextCr_128;
    uint32_t idx = (sb_origin_x + BLOCK_SIZE_64 > input_picture_ptr->width) ? sb_origin_x : 0;
    ////gaussian matrix(Chroma)
    //a = (1 * p[0] + 2 * p[1] + 1 * p[2] +
    //    2 * p[0 + stride] + 4 * p[1 + stride] + 2 * p[2 + stride] +
    //    1 * p[0 + 2 * stride] + 2 * p[1 + 2 * stride] + 1 * p[2 + 2 * stride]) / 16;

    {
        picHeight = input_picture_ptr->height / 2;
        picWidth = input_picture_ptr->width / 2;

        sb_height = MIN(BLOCK_SIZE_64 / 2, picHeight - sb_origin_y);

        sb_height = ((sb_origin_y + BLOCK_SIZE_64 / 2 >= picHeight) || (sb_origin_y == 0)) ? sb_height - 1 : sb_height;
        stride_in = input_picture_ptr->stride_cb;
        inputOriginIndex = input_picture_ptr->origin_x / 2 + (input_picture_ptr->origin_y / 2 + sb_origin_y)* input_picture_ptr->stride_cb;
        ptrIn = &(input_picture_ptr->buffer_cb[inputOriginIndex]);

        inputOriginIndexPad = denoised_picture_ptr->origin_x / 2 + (denoised_picture_ptr->origin_y / 2 + sb_origin_y)* denoised_picture_ptr->stride_cb;
        strideOut = denoised_picture_ptr->stride_cb;
        ptr_denoised = &(denoised_picture_ptr->buffer_cb[inputOriginIndexPad]);
        ptrDenoisedInterm = ptr_denoised;

        strideInCr = input_picture_ptr->stride_cr;
        inputOriginIndex = input_picture_ptr->origin_x / 2 + (input_picture_ptr->origin_y / 2 + sb_origin_y)  * input_picture_ptr->stride_cr;
        ptrInCr = &(input_picture_ptr->buffer_cr[inputOriginIndex]);

        inputOriginIndexPad = denoised_picture_ptr->origin_x / 2 + (denoised_picture_ptr->origin_y / 2 + sb_origin_y)  * denoised_picture_ptr->stride_cr;
        strideOutCr = denoised_picture_ptr->stride_cr;
        ptrDenoisedCr = &(denoised_picture_ptr->buffer_cr[inputOriginIndexPad]);
        ptrDenoisedIntermCr = ptrDenoisedCr;

        top = curr = top_next = top_prev = curr_next = curr_prev = topCr = currCr = topNextCr = topPrevCr = currNextCr = currPrevCr = _mm256_setzero_si256();
        for (kk = idx; kk + BLOCK_SIZE_64 / 2 <= picWidth; kk += BLOCK_SIZE_64 / 2)
        {
            for (jj = 0; jj < sb_height; jj++)
            {
                if (sb_origin_y == 0)
                {
                    if (jj == 0)
                    {
                        top = _mm256_loadu_si256((__m256i*)(ptrIn + kk + jj * stride_in));
                        curr = _mm256_loadu_si256((__m256i*)(ptrIn + kk + (1 + jj)*stride_in));
                        top_prev = _mm256_loadu_si256((__m256i*)(ptrIn - 1 + kk + ((jj)*stride_in)));
                        top_next = _mm256_loadu_si256((__m256i*)(ptrIn + 1 + kk + ((jj)*stride_in)));
                        curr_prev = _mm256_loadu_si256((__m256i*)(ptrIn - 1 + kk + ((1 + jj)*stride_in)));
                        curr_next = _mm256_loadu_si256((__m256i*)(ptrIn + 1 + kk + ((1 + jj)*stride_in)));
                        _mm256_storeu_si256((__m256i *)(ptr_denoised + kk), top);
                        topCr = _mm256_loadu_si256((__m256i*)(ptrInCr + kk + jj * strideInCr));
                        currCr = _mm256_loadu_si256((__m256i*)(ptrInCr + kk + (1 + jj)*strideInCr));
                        topPrevCr = _mm256_loadu_si256((__m256i*)(ptrInCr - 1 + kk + ((jj)*strideInCr)));
                        topNextCr = _mm256_loadu_si256((__m256i*)(ptrInCr + 1 + kk + ((jj)*strideInCr)));
                        currPrevCr = _mm256_loadu_si256((__m256i*)(ptrInCr - 1 + kk + ((1 + jj)*strideInCr)));
                        currNextCr = _mm256_loadu_si256((__m256i*)(ptrInCr + 1 + kk + ((1 + jj)*strideInCr)));
                        _mm256_storeu_si256((__m256i *)(ptrDenoisedCr + kk), topCr);
                    }
                    bottom_prev = _mm256_loadu_si256((__m256i*)(ptrIn - 1 + kk + ((2 + jj)*stride_in)));
                    bottom_next = _mm256_loadu_si256((__m256i*)(ptrIn + 1 + kk + ((2 + jj)*stride_in)));
                    bottom = _mm256_loadu_si256((__m256i*)((ptrIn + kk) + (2 + jj)* stride_in));
                    ptrDenoisedInterm = ptr_denoised + kk + ((1 + jj)*strideOut);
                    bottomPrevCr = _mm256_loadu_si256((__m256i*)(ptrInCr - 1 + kk + ((2 + jj)*strideInCr)));
                    bottomNextCr = _mm256_loadu_si256((__m256i*)(ptrInCr + 1 + kk + ((2 + jj)*strideInCr)));
                    bottomCr = _mm256_loadu_si256((__m256i*)((ptrInCr + kk) + (2 + jj)* strideInCr));
                    ptrDenoisedIntermCr = ptrDenoisedCr + kk + ((1 + jj)*strideOutCr);
                }
                else
                {
                    if (jj == 0)
                    {
                        top = _mm256_loadu_si256((__m256i*)(ptrIn + kk + jj * stride_in - stride_in));
                        curr = _mm256_loadu_si256((__m256i*)(ptrIn + kk + (1 + jj)*stride_in - stride_in));
                        top_prev = _mm256_loadu_si256((__m256i*)(ptrIn - 1 + kk + ((jj)*stride_in) - stride_in));
                        top_next = _mm256_loadu_si256((__m256i*)(ptrIn + 1 + kk + ((jj)*stride_in) - stride_in));
                        curr_prev = _mm256_loadu_si256((__m256i*)(ptrIn - 1 + kk + ((1 + jj)*stride_in - stride_in)));
                        curr_next = _mm256_loadu_si256((__m256i*)(ptrIn + 1 + kk + ((1 + jj)*stride_in - stride_in)));
                        topCr = _mm256_loadu_si256((__m256i*)(ptrInCr + kk + jj * strideInCr - strideInCr));
                        currCr = _mm256_loadu_si256((__m256i*)(ptrInCr + kk + (1 + jj)*strideInCr - strideInCr));
                        topPrevCr = _mm256_loadu_si256((__m256i*)(ptrInCr - 1 + kk + ((jj)*strideInCr) - strideInCr));
                        topNextCr = _mm256_loadu_si256((__m256i*)(ptrInCr + 1 + kk + ((jj)*strideInCr) - strideInCr));
                        currPrevCr = _mm256_loadu_si256((__m256i*)(ptrInCr - 1 + kk + ((1 + jj)*strideInCr - strideInCr)));
                        currNextCr = _mm256_loadu_si256((__m256i*)(ptrInCr + 1 + kk + ((1 + jj)*strideInCr - strideInCr)));
                    }
                    bottom_prev = _mm256_loadu_si256((__m256i*)(ptrIn - 1 + kk + ((2 + jj)*stride_in) - stride_in));
                    bottom_next = _mm256_loadu_si256((__m256i*)(ptrIn + 1 + kk + ((2 + jj)*stride_in) - stride_in));
                    bottom = _mm256_loadu_si256((__m256i*)((ptrIn + kk) + (2 + jj)* stride_in - stride_in));
                    ptrDenoisedInterm = ptr_denoised + kk + ((1 + jj)*strideOut - strideOut);
                    bottomPrevCr = _mm256_loadu_si256((__m256i*)(ptrInCr - 1 + kk + ((2 + jj)*strideInCr) - strideInCr));
                    bottomNextCr = _mm256_loadu_si256((__m256i*)(ptrInCr + 1 + kk + ((2 + jj)*strideInCr) - strideInCr));
                    bottomCr = _mm256_loadu_si256((__m256i*)((ptrInCr + kk) + (2 + jj)* strideInCr - strideInCr));
                    ptrDenoisedIntermCr = ptrDenoisedCr + kk + ((1 + jj)*strideOutCr - strideOutCr);
                }

                chroma_weak_luma_strong_filter_avx2_intrin(
                    top,
                    curr,
                    bottom,
                    curr_prev,
                    curr_next,
                    top_prev,
                    top_next,
                    bottom_prev,
                    bottom_next,
                    ptrDenoisedInterm);

                chroma_weak_luma_strong_filter_avx2_intrin(
                    topCr,
                    currCr,
                    bottomCr,
                    currPrevCr,
                    currNextCr,
                    topPrevCr,
                    topNextCr,
                    bottomPrevCr,
                    bottomNextCr,
                    ptrDenoisedIntermCr);

                top = curr;
                curr = bottom;
                top_prev = curr_prev;
                top_next = curr_next;
                curr_prev = bottom_prev;
                curr_next = bottom_next;
                topCr = currCr;
                currCr = bottomCr;
                topPrevCr = currPrevCr;
                topNextCr = currNextCr;
                currPrevCr = bottomPrevCr;
                currNextCr = bottomNextCr;
            }
        }

        for (; kk + 8 <= picWidth; kk += 8)
        {
            for (jj = 0; jj < sb_height; jj++)
            {
                if (sb_origin_y == 0)
                {
                    if (jj == 0)
                    {
                        top_128 = _mm_loadl_epi64((__m128i*)(ptrIn + kk + jj * stride_in));
                        curr_128 = _mm_loadl_epi64((__m128i*)(ptrIn + kk + (1 + jj) * stride_in));
                        top_prev_128 = _mm_loadl_epi64((__m128i*)(ptrIn - 1 + kk + ((jj)*stride_in)));
                        top_next_128 = _mm_loadl_epi64((__m128i*)(ptrIn + 1 + kk + ((jj)*stride_in)));
                        curr_prev_128 = _mm_loadl_epi64((__m128i*)(ptrIn - 1 + kk + ((1 + jj) * stride_in)));
                        curr_next_128 = _mm_loadl_epi64((__m128i*)(ptrIn + 1 + kk + ((1 + jj) * stride_in)));
                        _mm_storel_epi64((__m128i*)(ptr_denoised + kk), top_128);
                        topCr_128 = _mm_loadl_epi64((__m128i*)(ptrInCr + kk + jj * strideInCr));
                        currCr_128 = _mm_loadl_epi64((__m128i*)(ptrInCr + kk + (1 + jj) * strideInCr));
                        topPrevCr_128 = _mm_loadl_epi64((__m128i*)(ptrInCr - 1 + kk + ((jj)*strideInCr)));
                        topNextCr_128 = _mm_loadl_epi64((__m128i*)(ptrInCr + 1 + kk + ((jj)*strideInCr)));
                        currPrevCr_128 = _mm_loadl_epi64((__m128i*)(ptrInCr - 1 + kk + ((1 + jj) * strideInCr)));
                        currNextCr_128 = _mm_loadl_epi64((__m128i*)(ptrInCr + 1 + kk + ((1 + jj) * strideInCr)));
                        _mm_storel_epi64((__m128i*)(ptrDenoisedCr + kk), topCr_128);
                    }
                    bottom_prev_128 = _mm_loadl_epi64((__m128i*)(ptrIn - 1 + kk + ((2 + jj) * stride_in)));
                    bottom_next_128 = _mm_loadl_epi64((__m128i*)(ptrIn + 1 + kk + ((2 + jj) * stride_in)));
                    bottom_128 = _mm_loadl_epi64((__m128i*)((ptrIn + kk) + (2 + jj) * stride_in));
                    ptrDenoisedInterm = ptr_denoised + kk + ((1 + jj) * strideOut);
                    bottomPrevCr_128 = _mm_loadl_epi64((__m128i*)(ptrInCr - 1 + kk + ((2 + jj) * strideInCr)));
                    bottomNextCr_128 = _mm_loadl_epi64((__m128i*)(ptrInCr + 1 + kk + ((2 + jj) * strideInCr)));
                    bottomCr_128 = _mm_loadl_epi64((__m128i*)((ptrInCr + kk) + (2 + jj) * strideInCr));
                    ptrDenoisedIntermCr = ptrDenoisedCr + kk + ((1 + jj) * strideOutCr);
                }
                else
                {
                    if (jj == 0)
                    {
                        top_128 = _mm_loadl_epi64((__m128i*)(ptrIn + kk + jj * stride_in - stride_in));
                        curr_128 = _mm_loadl_epi64((__m128i*)(ptrIn + kk + (1 + jj) * stride_in - stride_in));
                        top_prev_128 = _mm_loadl_epi64((__m128i*)(ptrIn - 1 + kk + ((jj)*stride_in) - stride_in));
                        top_next_128 = _mm_loadl_epi64((__m128i*)(ptrIn + 1 + kk + ((jj)*stride_in) - stride_in));
                        curr_prev_128 = _mm_loadl_epi64((__m128i*)(ptrIn - 1 + kk + ((1 + jj) * stride_in - stride_in)));
                        curr_next_128 = _mm_loadl_epi64((__m128i*)(ptrIn + 1 + kk + ((1 + jj) * stride_in - stride_in)));
                        topCr_128 = _mm_loadl_epi64((__m128i*)(ptrInCr + kk + jj * strideInCr - strideInCr));
                        currCr_128 = _mm_loadl_epi64((__m128i*)(ptrInCr + kk + (1 + jj) * strideInCr - strideInCr));
                        topPrevCr_128 = _mm_loadl_epi64((__m128i*)(ptrInCr - 1 + kk + ((jj)*strideInCr) - strideInCr));
                        topNextCr_128 = _mm_loadl_epi64((__m128i*)(ptrInCr + 1 + kk + ((jj)*strideInCr) - strideInCr));
                        currPrevCr_128 = _mm_loadl_epi64((__m128i*)(ptrInCr - 1 + kk + ((1 + jj) * strideInCr - strideInCr)));
                        currNextCr_128 = _mm_loadl_epi64((__m128i*)(ptrInCr + 1 + kk + ((1 + jj) * strideInCr - strideInCr)));
                    }
                    bottom_prev_128 = _mm_loadl_epi64((__m128i*)(ptrIn - 1 + kk + ((2 + jj) * stride_in) - stride_in));
                    bottom_next_128 = _mm_loadl_epi64((__m128i*)(ptrIn + 1 + kk + ((2 + jj) * stride_in) - stride_in));
                    bottom_128 = _mm_loadl_epi64((__m128i*)((ptrIn + kk) + (2 + jj) * stride_in - stride_in));
                    ptrDenoisedInterm = ptr_denoised + kk + ((1 + jj) * strideOut - strideOut);
                    bottomPrevCr_128 = _mm_loadl_epi64((__m128i*)(ptrInCr - 1 + kk + ((2 + jj) * strideInCr) - strideInCr));
                    bottomNextCr_128 = _mm_loadl_epi64((__m128i*)(ptrInCr + 1 + kk + ((2 + jj) * strideInCr) - strideInCr));
                    bottomCr_128 = _mm_loadl_epi64((__m128i*)((ptrInCr + kk) + (2 + jj) * strideInCr - strideInCr));
                    ptrDenoisedIntermCr = ptrDenoisedCr + kk + ((1 + jj) * strideOutCr - strideOutCr);
                }

                chroma_weak_luma_strong_filter_128_avx2_intrin(
                    top_128,
                    curr_128,
                    bottom_128,
                    curr_prev_128,
                    curr_next_128,
                    top_prev_128,
                    top_next_128,
                    bottom_prev_128,
                    bottom_next_128,
                    ptrDenoisedInterm);

                chroma_weak_luma_strong_filter_128_avx2_intrin(
                    topCr_128,
                    currCr_128,
                    bottomCr_128,
                    currPrevCr_128,
                    currNextCr_128,
                    topPrevCr_128,
                    topNextCr_128,
                    bottomPrevCr_128,
                    bottomNextCr_128,
                    ptrDenoisedIntermCr);

                top_128 = curr_128;
                curr_128 = bottom_128;
                top_prev_128 = curr_prev_128;
                top_next_128 = curr_next_128;
                curr_prev_128 = bottom_prev_128;
                curr_next_128 = bottom_next_128;
                topCr_128 = currCr_128;
                currCr_128 = bottomCr_128;
                topPrevCr_128 = currPrevCr_128;
                topNextCr_128 = currNextCr_128;
                currPrevCr_128 = bottomPrevCr_128;
                currNextCr_128 = bottomNextCr_128;
            }
        }

        sb_height = MIN(BLOCK_SIZE_64 / 2, picHeight - sb_origin_y);
        for (jj = 0; jj < sb_height; jj++) {
            for (ii = idx; ii < picWidth; ii++) {
                if (!((jj < sb_height - 1 || (sb_origin_y + sb_height) < picHeight) && ii > 0 && ii < picWidth - 1)) {
                    ptr_denoised[ii + jj * strideOut] = ptrIn[ii + jj * stride_in];
                    ptrDenoisedCr[ii + jj * strideOut] = ptrInCr[ii + jj * strideInCr];
                }
            }
        }
    }
}
