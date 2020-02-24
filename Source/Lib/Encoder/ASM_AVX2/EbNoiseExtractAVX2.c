/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbNoiseExtractAVX2.h"
#include "EbDefinitions.h"
#include "immintrin.h"
#include "EbUtility.h"

EB_EXTERN EB_ALIGN(16) const uint8_t filter_type[] = {
    1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4};

EB_EXTERN EB_ALIGN(16) const uint8_t weak_chroma_filter[2][32] = {
    {2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4,
     2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4},
    {1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2,
     1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2},
};

static inline void luma_weak_filter_avx2_intrin(__m256i top, __m256i curr, __m256i bottom,
                                         __m256i curr_prev, __m256i curr_next,
                                         uint8_t *ptr_denoised, uint8_t *ptr_noise) {
    __m256i top_first_half, bottom_first_half, filter_first_half, filter_second_half,
        curr_next_first_half, curr_next_second_half, weights, curr_left_mid_first_half_weight,
        curr_left_mid_first_halflo, curr_left_mid_first_halfhi, curr_prev_permutation,
        curr_permutation, curr_next_permutation, top_permutation, bottom_permutation;

    curr_prev_permutation           = _mm256_permute4x64_epi64(curr_prev, 216);
    curr_permutation                = _mm256_permute4x64_epi64(curr, 216);
    curr_left_mid_first_halflo      = _mm256_unpacklo_epi8(curr_prev_permutation, curr_permutation);
    weights                         = _mm256_loadu_si256((__m256i *)filter_type);
    curr_left_mid_first_half_weight = _mm256_maddubs_epi16(curr_left_mid_first_halflo, weights);
    curr_next_permutation           = _mm256_permute4x64_epi64(curr_next, 88);
    curr_next_first_half = _mm256_unpacklo_epi8(curr_next_permutation, _mm256_setzero_si256());
    curr_left_mid_first_halflo =
        _mm256_add_epi16(curr_next_first_half, curr_left_mid_first_half_weight);

    curr_left_mid_first_halfhi      = _mm256_unpackhi_epi8(curr_prev_permutation, curr_permutation);
    curr_left_mid_first_half_weight = _mm256_maddubs_epi16(curr_left_mid_first_halfhi, weights);
    curr_next_permutation           = _mm256_permute4x64_epi64(curr_next, 216);
    curr_next_second_half = _mm256_unpackhi_epi8(curr_next_permutation, _mm256_setzero_si256());
    curr_left_mid_first_halfhi =
        _mm256_add_epi16(curr_next_second_half, curr_left_mid_first_half_weight);

    top_permutation    = _mm256_permute4x64_epi64(top, 216);
    top_first_half     = _mm256_unpacklo_epi8(top_permutation, _mm256_setzero_si256());
    bottom_permutation = _mm256_permute4x64_epi64(bottom, 216);
    bottom_first_half  = _mm256_unpacklo_epi8(bottom_permutation, _mm256_setzero_si256());
    filter_first_half  = _mm256_adds_epi16(_mm256_adds_epi16(bottom_first_half, top_first_half),
                                          curr_left_mid_first_halflo);
    filter_first_half  = _mm256_srli_epi16(filter_first_half, 3);

    top_first_half     = _mm256_unpackhi_epi8(top_permutation, _mm256_setzero_si256());
    bottom_first_half  = _mm256_unpackhi_epi8(bottom_permutation, _mm256_setzero_si256());
    filter_second_half = _mm256_adds_epi16(_mm256_adds_epi16(bottom_first_half, top_first_half),
                                           curr_left_mid_first_halfhi);
    filter_second_half = _mm256_srli_epi16(filter_second_half, 3);

    filter_first_half =
        _mm256_permute4x64_epi64(_mm256_packus_epi16(filter_first_half, filter_second_half), 216);
    _mm256_storeu_si256((__m256i *)(ptr_denoised), filter_first_half);

    _mm256_storeu_si256((__m256i *)(ptr_noise), _mm256_subs_epu8(curr, filter_first_half));
}
static inline void luma_weak_filter_128_avx2_intrin(__m128i top, __m128i curr, __m128i bottom,
                                             __m128i curr_prev, __m128i curr_next,
                                             uint8_t *ptr_denoised, uint8_t *ptr_noise) {

    __m128i top_first_half, bottom_first_half, filter_first_half, filter_second_half,
        curr_next_first_half, curr_next_second_half, weights, curr_left_mid_first_half_weight,
        curr_left_mid_first_halflo, curr_left_mid_first_halfhi;

    curr_left_mid_first_halflo      = _mm_unpacklo_epi8(curr_prev, curr);
    weights                         = _mm_loadu_si128((__m128i *)filter_type);
    curr_left_mid_first_half_weight = _mm_maddubs_epi16(curr_left_mid_first_halflo, weights);
    curr_next_first_half            = _mm_unpacklo_epi8(curr_next, _mm_setzero_si128());
    curr_left_mid_first_halflo =
        _mm_add_epi16(curr_next_first_half, curr_left_mid_first_half_weight);

    curr_left_mid_first_halfhi      = _mm_unpackhi_epi8(curr_prev, curr);
    curr_left_mid_first_half_weight = _mm_maddubs_epi16(curr_left_mid_first_halfhi, weights);
    curr_next_second_half           = _mm_unpackhi_epi8(curr_next, _mm_setzero_si128());
    curr_left_mid_first_halfhi =
        _mm_add_epi16(curr_next_second_half, curr_left_mid_first_half_weight);

    top_first_half    = _mm_unpacklo_epi8(top, _mm_setzero_si128());
    bottom_first_half = _mm_unpacklo_epi8(bottom, _mm_setzero_si128());
    filter_first_half = _mm_adds_epi16(_mm_adds_epi16(bottom_first_half, top_first_half),
                                       curr_left_mid_first_halflo);
    filter_first_half = _mm_srli_epi16(filter_first_half, 3);

    top_first_half     = _mm_unpackhi_epi8(top, _mm_setzero_si128());
    bottom_first_half  = _mm_unpackhi_epi8(bottom, _mm_setzero_si128());
    filter_second_half = _mm_adds_epi16(_mm_adds_epi16(bottom_first_half, top_first_half),
                                        curr_left_mid_first_halfhi);
    filter_second_half = _mm_srli_epi16(filter_second_half, 3);

    filter_first_half = _mm_packus_epi16(filter_first_half, filter_second_half);
    _mm_storel_epi64((__m128i *)(ptr_denoised), filter_first_half);

    _mm_storel_epi64((__m128i *)(ptr_noise), _mm_subs_epu8(curr, filter_first_half));
}
static inline void chroma_weak_luma_strong_filter_avx2_intrin(__m256i top, __m256i curr, __m256i bottom,
                                                       __m256i curr_prev, __m256i curr_next,
                                                       __m256i top_prev, __m256i top_next,
                                                       __m256i bottom_prev, __m256i bottom_next,
                                                       uint8_t *ptr_denoised) {
    __m256i filter_first_half, filter_second_half, curr_next_first_half, curr_next_second_half,
        weights, curr_left_mid_first_half_weight, curr_left_mid_first_halflo,
        curr_left_mid_first_halfhi, curr_prev_permutation, curr_permutation, curr_next_permutation,
        top_permutation, bottom_permutation, top_prev_permutation, top_left_mid_first_halflo,
        top_left_mid_first_half_weight, top_next_first_half, top_next_permutation,
        top_left_mid_first_halfhi, top_next_second_half, bottom_prev_permutation,
        bottom_left_mid_first_halflo, bottom_left_mid_first_half_weight, bottom_next_permutation,
        bottom_next_first_half, bottom_left_mid_first_halfhi, bottom_next_second_half;

    //  Curr
    curr_prev_permutation           = _mm256_permute4x64_epi64(curr_prev, 216);
    curr_permutation                = _mm256_permute4x64_epi64(curr, 216);
    curr_left_mid_first_halflo      = _mm256_unpacklo_epi8(curr_prev_permutation, curr_permutation);
    weights                         = _mm256_loadu_si256((__m256i *)weak_chroma_filter[0]);
    curr_left_mid_first_half_weight = _mm256_maddubs_epi16(curr_left_mid_first_halflo, weights);
    curr_next_permutation           = _mm256_permute4x64_epi64(curr_next, 88);
    curr_next_first_half = _mm256_unpacklo_epi8(curr_next_permutation, _mm256_setzero_si256());
    curr_next_first_half = _mm256_slli_epi16(curr_next_first_half, 1);
    curr_left_mid_first_halflo =
        _mm256_add_epi16(curr_next_first_half, curr_left_mid_first_half_weight);

    curr_left_mid_first_halfhi      = _mm256_unpackhi_epi8(curr_prev_permutation, curr_permutation);
    curr_left_mid_first_half_weight = _mm256_maddubs_epi16(curr_left_mid_first_halfhi, weights);
    curr_next_permutation           = _mm256_permute4x64_epi64(curr_next, 216);
    curr_next_second_half = _mm256_unpackhi_epi8(curr_next_permutation, _mm256_setzero_si256());
    curr_next_second_half = _mm256_slli_epi16(curr_next_second_half, 1);
    curr_left_mid_first_halfhi =
        _mm256_add_epi16(curr_next_second_half, curr_left_mid_first_half_weight);

    // Top
    top_prev_permutation           = _mm256_permute4x64_epi64(top_prev, 216);
    top_permutation                = _mm256_permute4x64_epi64(top, 216);
    top_left_mid_first_halflo      = _mm256_unpacklo_epi8(top_prev_permutation, top_permutation);
    weights                        = _mm256_loadu_si256((__m256i *)weak_chroma_filter[1]);
    top_left_mid_first_half_weight = _mm256_maddubs_epi16(top_left_mid_first_halflo, weights);
    top_next_permutation           = _mm256_permute4x64_epi64(top_next, 88);
    top_next_first_half = _mm256_unpacklo_epi8(top_next_permutation, _mm256_setzero_si256());
    top_left_mid_first_halflo =
        _mm256_add_epi16(top_next_first_half, top_left_mid_first_half_weight);

    top_left_mid_first_halfhi      = _mm256_unpackhi_epi8(top_prev_permutation, top_permutation);
    top_left_mid_first_half_weight = _mm256_maddubs_epi16(top_left_mid_first_halfhi, weights);
    top_next_permutation           = _mm256_permute4x64_epi64(top_next, 216);
    top_next_second_half = _mm256_unpackhi_epi8(top_next_permutation, _mm256_setzero_si256());
    top_left_mid_first_halfhi =
        _mm256_add_epi16(top_next_second_half, top_left_mid_first_half_weight);

    // Bottom
    bottom_prev_permutation = _mm256_permute4x64_epi64(bottom_prev, 216);
    bottom_permutation      = _mm256_permute4x64_epi64(bottom, 216);
    bottom_left_mid_first_halflo =
        _mm256_unpacklo_epi8(bottom_prev_permutation, bottom_permutation);
    weights                           = _mm256_loadu_si256((__m256i *)weak_chroma_filter[1]);
    bottom_left_mid_first_half_weight = _mm256_maddubs_epi16(bottom_left_mid_first_halflo, weights);
    bottom_next_permutation           = _mm256_permute4x64_epi64(bottom_next, 88);
    bottom_next_first_half = _mm256_unpacklo_epi8(bottom_next_permutation, _mm256_setzero_si256());
    bottom_left_mid_first_halflo =
        _mm256_add_epi16(bottom_next_first_half, bottom_left_mid_first_half_weight);

    bottom_left_mid_first_halfhi =
        _mm256_unpackhi_epi8(bottom_prev_permutation, bottom_permutation);
    bottom_left_mid_first_half_weight = _mm256_maddubs_epi16(bottom_left_mid_first_halfhi, weights);
    bottom_next_permutation           = _mm256_permute4x64_epi64(bottom_next, 216);
    bottom_next_second_half = _mm256_unpackhi_epi8(bottom_next_permutation, _mm256_setzero_si256());
    bottom_left_mid_first_halfhi =
        _mm256_add_epi16(bottom_next_second_half, bottom_left_mid_first_half_weight);

    filter_first_half = _mm256_adds_epi16(
        _mm256_adds_epi16(bottom_left_mid_first_halflo, top_left_mid_first_halflo),
        curr_left_mid_first_halflo);
    filter_first_half  = _mm256_srli_epi16(filter_first_half, 4);
    filter_second_half = _mm256_adds_epi16(
        _mm256_adds_epi16(bottom_left_mid_first_halfhi, top_left_mid_first_halfhi),
        curr_left_mid_first_halfhi);
    filter_second_half = _mm256_srli_epi16(filter_second_half, 4);

    filter_first_half =
        _mm256_permute4x64_epi64(_mm256_packus_epi16(filter_first_half, filter_second_half), 216);
    _mm256_storeu_si256((__m256i *)(ptr_denoised), filter_first_half);
}

static inline void chroma_weak_luma_strong_filter_128_avx2_intrin(__m128i top, __m128i curr,
                                                           __m128i bottom, __m128i curr_prev,
                                                           __m128i curr_next, __m128i top_prev,
                                                           __m128i top_next, __m128i bottom_prev,
                                                           __m128i  bottom_next,
                                                           uint8_t *ptr_denoised) {
    __m128i filter_first_half, filter_second_half, curr_next_first_half, curr_next_second_half,
        weights, curr_left_mid_first_half_weight, curr_left_mid_first_halflo,
        curr_left_mid_first_halfhi, top_left_mid_first_halflo, top_left_mid_first_half_weight,
        top_next_first_half, top_left_mid_first_halfhi, top_next_second_half,
        bottom_left_mid_first_halflo, bottom_left_mid_first_half_weight, bottom_next_first_half,
        bottom_left_mid_first_halfhi, bottom_next_second_half;

    //  Curr
    curr_left_mid_first_halflo      = _mm_unpacklo_epi8(curr_prev, curr);
    weights                         = _mm_loadu_si128((__m128i *)weak_chroma_filter[0]);
    curr_left_mid_first_half_weight = _mm_maddubs_epi16(curr_left_mid_first_halflo, weights);
    curr_next_first_half            = _mm_unpacklo_epi8(curr_next, _mm_setzero_si128());
    curr_next_first_half            = _mm_slli_epi16(curr_next_first_half, 1);
    curr_left_mid_first_halflo =
        _mm_add_epi16(curr_next_first_half, curr_left_mid_first_half_weight);

    curr_left_mid_first_halfhi      = _mm_unpackhi_epi8(curr_prev, curr);
    curr_left_mid_first_half_weight = _mm_maddubs_epi16(curr_left_mid_first_halfhi, weights);
    curr_next_second_half           = _mm_unpackhi_epi8(curr_next, _mm_setzero_si128());
    curr_next_second_half           = _mm_slli_epi16(curr_next_second_half, 1);
    curr_left_mid_first_halfhi =
        _mm_add_epi16(curr_next_second_half, curr_left_mid_first_half_weight);

    // Top
    top_left_mid_first_halflo      = _mm_unpacklo_epi8(top_prev, top);
    weights                        = _mm_loadu_si128((__m128i *)weak_chroma_filter[1]);
    top_left_mid_first_half_weight = _mm_maddubs_epi16(top_left_mid_first_halflo, weights);
    top_next_first_half            = _mm_unpacklo_epi8(top_next, _mm_setzero_si128());
    top_left_mid_first_halflo = _mm_add_epi16(top_next_first_half, top_left_mid_first_half_weight);

    top_left_mid_first_halfhi      = _mm_unpackhi_epi8(top_prev, top);
    top_left_mid_first_half_weight = _mm_maddubs_epi16(top_left_mid_first_halfhi, weights);
    top_next_second_half           = _mm_unpackhi_epi8(top_next, _mm_setzero_si128());
    top_left_mid_first_halfhi = _mm_add_epi16(top_next_second_half, top_left_mid_first_half_weight);

    // Bottom
    bottom_left_mid_first_halflo      = _mm_unpacklo_epi8(bottom_prev, bottom);
    weights                           = _mm_loadu_si128((__m128i *)weak_chroma_filter[1]);
    bottom_left_mid_first_half_weight = _mm_maddubs_epi16(bottom_left_mid_first_halflo, weights);
    bottom_next_first_half            = _mm_unpacklo_epi8(bottom_next, _mm_setzero_si128());
    bottom_left_mid_first_halflo =
        _mm_add_epi16(bottom_next_first_half, bottom_left_mid_first_half_weight);

    bottom_left_mid_first_halfhi      = _mm_unpackhi_epi8(bottom_prev, bottom);
    bottom_left_mid_first_half_weight = _mm_maddubs_epi16(bottom_left_mid_first_halfhi, weights);
    bottom_next_second_half           = _mm_unpackhi_epi8(bottom_next, _mm_setzero_si128());
    bottom_left_mid_first_halfhi =
        _mm_add_epi16(bottom_next_second_half, bottom_left_mid_first_half_weight);

    filter_first_half =
        _mm_adds_epi16(_mm_adds_epi16(bottom_left_mid_first_halflo, top_left_mid_first_halflo),
                       curr_left_mid_first_halflo);
    filter_first_half = _mm_srli_epi16(filter_first_half, 4);
    filter_second_half =
        _mm_adds_epi16(_mm_adds_epi16(bottom_left_mid_first_halfhi, top_left_mid_first_halfhi),
                       curr_left_mid_first_halfhi);
    filter_second_half = _mm_srli_epi16(filter_second_half, 4);

    filter_first_half = _mm_packus_epi16(filter_first_half, filter_second_half);
    _mm_storel_epi64((__m128i *)(ptr_denoised), filter_first_half);
}

static inline void chroma_strong_avx2_intrin(__m256i top, __m256i curr, __m256i bottom, __m256i curr_prev,
                                      __m256i curr_next, __m256i top_prev, __m256i top_next,
                                      __m256i bottom_prev, __m256i bottom_next,
                                      uint8_t *ptr_denoised) {
    __m256i curr_left_mid_first_halflo, curr_left_mid_first_halfhi, curr_prev_permutation,
        curr_permutation, curr_next_permutation, top_permutation, top_prev_permutation,
        top_left_mid_first_halflo, top_next_permutation, top_left_mid_first_halfhi,
        bottom_permutation, bottom_prev_permutation, bottom_left_mid_first_halflo,
        bottom_next_permutation, bottom_left_mid_first_halfhi;

    curr_prev_permutation = _mm256_permute4x64_epi64(curr_prev, 216);
    curr_permutation      = _mm256_permute4x64_epi64(curr, 216);
    curr_next_permutation = _mm256_permute4x64_epi64(curr_next, 216);

    curr_left_mid_first_halflo =
        _mm256_add_epi16(_mm256_unpacklo_epi8(curr_permutation, _mm256_setzero_si256()),
                         _mm256_unpacklo_epi8(curr_prev_permutation, _mm256_setzero_si256()));
    curr_left_mid_first_halflo =
        _mm256_add_epi16(_mm256_unpacklo_epi8(curr_next_permutation, _mm256_setzero_si256()),
                         curr_left_mid_first_halflo);

    curr_left_mid_first_halfhi =
        _mm256_add_epi16(_mm256_unpackhi_epi8(curr_permutation, _mm256_setzero_si256()),
                         _mm256_unpackhi_epi8(curr_prev_permutation, _mm256_setzero_si256()));
    curr_left_mid_first_halfhi =
        _mm256_add_epi16(_mm256_unpackhi_epi8(curr_next_permutation, _mm256_setzero_si256()),
                         curr_left_mid_first_halfhi);

    top_prev_permutation = _mm256_permute4x64_epi64(top_prev, 216);
    top_permutation      = _mm256_permute4x64_epi64(top, 216);
    top_next_permutation = _mm256_permute4x64_epi64(top_next, 216);

    top_left_mid_first_halflo =
        _mm256_add_epi16(_mm256_unpacklo_epi8(top_permutation, _mm256_setzero_si256()),
                         _mm256_unpacklo_epi8(top_prev_permutation, _mm256_setzero_si256()));
    top_left_mid_first_halflo =
        _mm256_add_epi16(_mm256_unpacklo_epi8(top_next_permutation, _mm256_setzero_si256()),
                         top_left_mid_first_halflo);

    top_left_mid_first_halfhi =
        _mm256_add_epi16(_mm256_unpackhi_epi8(top_permutation, _mm256_setzero_si256()),
                         _mm256_unpackhi_epi8(top_prev_permutation, _mm256_setzero_si256()));
    top_left_mid_first_halfhi =
        _mm256_add_epi16(_mm256_unpackhi_epi8(top_next_permutation, _mm256_setzero_si256()),
                         top_left_mid_first_halfhi);

    bottom_prev_permutation = _mm256_permute4x64_epi64(bottom_prev, 216);
    bottom_permutation      = _mm256_permute4x64_epi64(bottom, 216);
    bottom_next_permutation = _mm256_permute4x64_epi64(bottom_next, 216);

    bottom_left_mid_first_halflo =
        _mm256_add_epi16(_mm256_unpacklo_epi8(bottom_permutation, _mm256_setzero_si256()),
                         _mm256_unpacklo_epi8(bottom_prev_permutation, _mm256_setzero_si256()));
    bottom_left_mid_first_halflo =
        _mm256_add_epi16(_mm256_unpacklo_epi8(bottom_next_permutation, _mm256_setzero_si256()),
                         bottom_left_mid_first_halflo);

    bottom_left_mid_first_halfhi =
        _mm256_add_epi16(_mm256_unpackhi_epi8(bottom_permutation, _mm256_setzero_si256()),
                         _mm256_unpackhi_epi8(bottom_prev_permutation, _mm256_setzero_si256()));
    bottom_left_mid_first_halfhi =
        _mm256_add_epi16(_mm256_unpackhi_epi8(bottom_next_permutation, _mm256_setzero_si256()),
                         bottom_left_mid_first_halfhi);

    curr_left_mid_first_halflo =
        _mm256_add_epi16(_mm256_add_epi16(curr_left_mid_first_halflo, top_left_mid_first_halflo),
                         bottom_left_mid_first_halflo);
    curr_left_mid_first_halfhi =
        _mm256_add_epi16(_mm256_add_epi16(curr_left_mid_first_halfhi, top_left_mid_first_halfhi),
                         bottom_left_mid_first_halfhi);

    top_left_mid_first_halflo =
        _mm256_unpacklo_epi16(curr_left_mid_first_halflo, _mm256_setzero_si256());
    top_left_mid_first_halflo =
        _mm256_mullo_epi32(top_left_mid_first_halflo, _mm256_set1_epi32(7282));
    top_left_mid_first_halflo = _mm256_srli_epi32(top_left_mid_first_halflo, 16);
    bottom_left_mid_first_halflo =
        _mm256_unpackhi_epi16(curr_left_mid_first_halflo, _mm256_setzero_si256());
    bottom_left_mid_first_halflo =
        _mm256_mullo_epi32(bottom_left_mid_first_halflo, _mm256_set1_epi32(7282));
    bottom_left_mid_first_halflo = _mm256_srli_epi32(bottom_left_mid_first_halflo, 16);
    curr_left_mid_first_halflo =
        _mm256_packus_epi32(top_left_mid_first_halflo, bottom_left_mid_first_halflo);

    curr_left_mid_first_halflo = _mm256_insertf128_si256(
        _mm256_setzero_si256(),
        _mm_packus_epi16(_mm256_extracti128_si256(curr_left_mid_first_halflo, 0),
                         _mm256_extracti128_si256(curr_left_mid_first_halflo, 1)),
        0);

    top_left_mid_first_halfhi =
        _mm256_unpacklo_epi16(curr_left_mid_first_halfhi, _mm256_setzero_si256());
    top_left_mid_first_halfhi =
        _mm256_mullo_epi32(top_left_mid_first_halfhi, _mm256_set1_epi32(7282));
    top_left_mid_first_halfhi = _mm256_srli_epi32(top_left_mid_first_halfhi, 16);

    bottom_left_mid_first_halfhi =
        _mm256_unpackhi_epi16(curr_left_mid_first_halfhi, _mm256_setzero_si256());
    bottom_left_mid_first_halfhi =
        _mm256_mullo_epi32(bottom_left_mid_first_halfhi, _mm256_set1_epi32(7282));
    bottom_left_mid_first_halfhi = _mm256_srli_epi32(bottom_left_mid_first_halfhi, 16);
    curr_left_mid_first_halfhi =
        _mm256_packus_epi32(top_left_mid_first_halfhi, bottom_left_mid_first_halfhi);

    curr_left_mid_first_halflo = _mm256_insertf128_si256(
        curr_left_mid_first_halflo,
        _mm_packus_epi16(_mm256_extracti128_si256(curr_left_mid_first_halfhi, 0),
                         _mm256_extracti128_si256(curr_left_mid_first_halfhi, 1)),
        1);
    _mm256_storeu_si256((__m256i *)(ptr_denoised), curr_left_mid_first_halflo);
}

static inline void chroma_strong_128_avx2_intrin(__m128i top, __m128i curr, __m128i bottom,
                                          __m128i curr_prev, __m128i curr_next, __m128i top_prev,
                                          __m128i top_next, __m128i bottom_prev,
                                          __m128i bottom_next, uint8_t *ptr_denoised) {
    __m128i curr_left_mid_first_halflo, curr_left_mid_first_halfhi, top_left_mid_first_halflo,
        top_left_mid_first_halfhi, bottom_left_mid_first_halflo, bottom_left_mid_first_halfhi;

    curr_left_mid_first_halflo = _mm_add_epi16(_mm_unpacklo_epi8(curr, _mm_setzero_si128()),
                                               _mm_unpacklo_epi8(curr_prev, _mm_setzero_si128()));
    curr_left_mid_first_halflo = _mm_add_epi16(_mm_unpacklo_epi8(curr_next, _mm_setzero_si128()),
                                               curr_left_mid_first_halflo);

    curr_left_mid_first_halfhi = _mm_add_epi16(_mm_unpackhi_epi8(curr, _mm_setzero_si128()),
                                               _mm_unpackhi_epi8(curr_prev, _mm_setzero_si128()));
    curr_left_mid_first_halfhi = _mm_add_epi16(_mm_unpackhi_epi8(curr_next, _mm_setzero_si128()),
                                               curr_left_mid_first_halfhi);

    top_left_mid_first_halflo = _mm_add_epi16(_mm_unpacklo_epi8(top, _mm_setzero_si128()),
                                              _mm_unpacklo_epi8(top_prev, _mm_setzero_si128()));
    top_left_mid_first_halflo =
        _mm_add_epi16(_mm_unpacklo_epi8(top_next, _mm_setzero_si128()), top_left_mid_first_halflo);

    top_left_mid_first_halfhi = _mm_add_epi16(_mm_unpackhi_epi8(top, _mm_setzero_si128()),
                                              _mm_unpackhi_epi8(top_prev, _mm_setzero_si128()));
    top_left_mid_first_halfhi =
        _mm_add_epi16(_mm_unpackhi_epi8(top_next, _mm_setzero_si128()), top_left_mid_first_halfhi);

    bottom_left_mid_first_halflo =
        _mm_add_epi16(_mm_unpacklo_epi8(bottom, _mm_setzero_si128()),
                      _mm_unpacklo_epi8(bottom_prev, _mm_setzero_si128()));
    bottom_left_mid_first_halflo = _mm_add_epi16(
        _mm_unpacklo_epi8(bottom_next, _mm_setzero_si128()), bottom_left_mid_first_halflo);

    bottom_left_mid_first_halfhi =
        _mm_add_epi16(_mm_unpackhi_epi8(bottom, _mm_setzero_si128()),
                      _mm_unpackhi_epi8(bottom_prev, _mm_setzero_si128()));
    bottom_left_mid_first_halfhi = _mm_add_epi16(
        _mm_unpackhi_epi8(bottom_next, _mm_setzero_si128()), bottom_left_mid_first_halfhi);

    curr_left_mid_first_halflo =
        _mm_add_epi16(_mm_add_epi16(curr_left_mid_first_halflo, top_left_mid_first_halflo),
                      bottom_left_mid_first_halflo);
    curr_left_mid_first_halfhi =
        _mm_add_epi16(_mm_add_epi16(curr_left_mid_first_halfhi, top_left_mid_first_halfhi),
                      bottom_left_mid_first_halfhi);

    top_left_mid_first_halflo = _mm_unpacklo_epi16(curr_left_mid_first_halflo, _mm_setzero_si128());
    top_left_mid_first_halflo = _mm_mullo_epi32(top_left_mid_first_halflo, _mm_set1_epi32(7282));
    top_left_mid_first_halflo = _mm_srli_epi32(top_left_mid_first_halflo, 16);
    bottom_left_mid_first_halflo =
        _mm_unpackhi_epi16(curr_left_mid_first_halflo, _mm_setzero_si128());
    bottom_left_mid_first_halflo =
        _mm_mullo_epi32(bottom_left_mid_first_halflo, _mm_set1_epi32(7282));
    bottom_left_mid_first_halflo = _mm_srli_epi32(bottom_left_mid_first_halflo, 16);
    curr_left_mid_first_halflo =
        _mm_packus_epi32(top_left_mid_first_halflo, bottom_left_mid_first_halflo);

    curr_left_mid_first_halflo =
        _mm_packus_epi16(curr_left_mid_first_halflo, curr_left_mid_first_halflo);

    _mm_storel_epi64((__m128i *)(ptr_denoised), curr_left_mid_first_halflo);
}
/*******************************************
* noise_extract_luma_weak
*  weak filter Luma and store noise.
*******************************************/
void noise_extract_luma_weak_avx2_intrin(EbPictureBufferDesc *input_picture_ptr,
                                         EbPictureBufferDesc *denoised_picture_ptr,
                                         EbPictureBufferDesc *noise_picture_ptr,
                                         uint32_t sb_origin_y, uint32_t sb_origin_x) {
    uint32_t ii, jj, kk;
    uint32_t pic_height, sb_height;
    uint32_t pic_width;
    uint32_t input_origin_index;
    uint32_t input_origin_index_pad;
    uint32_t noise_origin_index;

    uint8_t *ptr_in;
    uint32_t stride_in;
    uint8_t *ptr_denoised, *ptr_denoised_interm;

    uint8_t *ptr_noise, *ptr_noise_interm;
    uint32_t stride_out;

    __m256i top, curr, bottom, curr_prev, curr_next, second_top, second_curr, second_bottom,
        second_curr_prev, second_curr_next;
    __m128i top_128, curr_128, bottom_128, curr_prev_128, curr_next_128, second_top_128,
        second_curr_128, second_bottom_128, second_curr_prev_128, second_curr_next_128;
    uint32_t idx = (sb_origin_x + BLOCK_SIZE_64 > input_picture_ptr->width) ? sb_origin_x : 0;

    //Luma
    {
        pic_height = input_picture_ptr->height;
        pic_width  = input_picture_ptr->width;
        sb_height  = MIN(BLOCK_SIZE_64, pic_height - sb_origin_y);
        sb_height  = ((sb_origin_y + BLOCK_SIZE_64 >= pic_height) || (sb_origin_y == 0))
                        ? sb_height - 1
                        : sb_height;
        stride_in = input_picture_ptr->stride_y;
        input_origin_index =
            input_picture_ptr->origin_x +
            (input_picture_ptr->origin_y + sb_origin_y) * input_picture_ptr->stride_y;
        ptr_in = &(input_picture_ptr->buffer_y[input_origin_index]);

        input_origin_index_pad =
            denoised_picture_ptr->origin_x +
            (denoised_picture_ptr->origin_y + sb_origin_y) * denoised_picture_ptr->stride_y;
        stride_out          = denoised_picture_ptr->stride_y;
        ptr_denoised        = &(denoised_picture_ptr->buffer_y[input_origin_index_pad]);
        ptr_denoised_interm = ptr_denoised;

        noise_origin_index =
            noise_picture_ptr->origin_x + noise_picture_ptr->origin_y * noise_picture_ptr->stride_y;
        ptr_noise        = &(noise_picture_ptr->buffer_y[noise_origin_index]);
        ptr_noise_interm = ptr_noise;

        ////Luma
        //a = (p[1] +
        //    p[0 + stride] + 4 * p[1 + stride] + p[2 + stride] +
        //    p[1 + 2 * stride]) / 8;

        top = curr = second_top = second_curr = _mm256_setzero_si256();

        for (kk = idx; kk + BLOCK_SIZE_64 <= pic_width; kk += BLOCK_SIZE_64) {
            for (jj = 0; jj < sb_height; jj++) {
                if (sb_origin_y == 0) {
                    if (jj == 0) {
                        top = _mm256_loadu_si256((__m256i *)(ptr_in + kk + jj * stride_in));
                        second_top =
                            _mm256_loadu_si256((__m256i *)(ptr_in + kk + 32 + jj * stride_in));
                        curr = _mm256_loadu_si256((__m256i *)(ptr_in + kk + (1 + jj) * stride_in));
                        second_curr = _mm256_loadu_si256(
                            (__m256i *)((ptr_in + kk + 32) + (1 + jj) * stride_in));
                        _mm256_storeu_si256((__m256i *)(ptr_denoised + kk), top);
                        _mm256_storeu_si256((__m256i *)(ptr_denoised + kk + 32), second_top);
                        _mm256_storeu_si256((__m256i *)(ptr_noise + kk), _mm256_setzero_si256());
                        _mm256_storeu_si256((__m256i *)(ptr_noise + kk + 32),
                                            _mm256_setzero_si256());
                    }
                    curr_prev =
                        _mm256_loadu_si256((__m256i *)(ptr_in - 1 + kk + ((1 + jj) * stride_in)));
                    curr_next =
                        _mm256_loadu_si256((__m256i *)(ptr_in + 1 + kk + ((1 + jj) * stride_in)));
                    second_curr_prev = _mm256_loadu_si256(
                        (__m256i *)((ptr_in + kk + 32) - 1 + ((1 + jj) * stride_in)));
                    second_curr_next = _mm256_loadu_si256(
                        (__m256i *)((ptr_in + kk + 32) + 1 + ((1 + jj) * stride_in)));
                    bottom = _mm256_loadu_si256((__m256i *)((ptr_in + kk) + (2 + jj) * stride_in));
                    second_bottom =
                        _mm256_loadu_si256((__m256i *)((ptr_in + kk + 32) + (2 + jj) * stride_in));
                    ptr_denoised_interm = ptr_denoised + kk + ((1 + jj) * stride_out);
                    ptr_noise_interm    = ptr_noise + kk + ((1 + jj) * stride_out);
                } else {
                    if (jj == 0) {
                        top = _mm256_loadu_si256(
                            (__m256i *)(ptr_in + kk + jj * stride_in - stride_in));
                        second_top = _mm256_loadu_si256(
                            (__m256i *)(ptr_in + kk + 32 + jj * stride_in - stride_in));
                        curr = _mm256_loadu_si256(
                            (__m256i *)(ptr_in + kk + (1 + jj) * stride_in - stride_in));
                        second_curr = _mm256_loadu_si256(
                            (__m256i *)((ptr_in + kk + 32) + (1 + jj) * stride_in - stride_in));
                    }
                    curr_prev = _mm256_loadu_si256(
                        (__m256i *)(ptr_in - 1 + kk + ((1 + jj) * stride_in - stride_in)));
                    curr_next = _mm256_loadu_si256(
                        (__m256i *)(ptr_in + 1 + kk + ((1 + jj) * stride_in - stride_in)));
                    second_curr_prev = _mm256_loadu_si256(
                        (__m256i *)((ptr_in + kk + 32) - 1 + ((1 + jj) * stride_in - stride_in)));
                    second_curr_next = _mm256_loadu_si256(
                        (__m256i *)((ptr_in + kk + 32) + 1 + ((1 + jj) * stride_in - stride_in)));
                    bottom = _mm256_loadu_si256(
                        (__m256i *)((ptr_in + kk) + (2 + jj) * stride_in - stride_in));
                    second_bottom = _mm256_loadu_si256(
                        (__m256i *)((ptr_in + kk + 32) + (2 + jj) * stride_in - stride_in));
                    ptr_denoised_interm = ptr_denoised + kk + ((1 + jj) * stride_out - stride_out);
                    ptr_noise_interm    = ptr_noise + kk + jj * stride_out;
                }

                luma_weak_filter_avx2_intrin(
                    top, curr, bottom, curr_prev, curr_next, ptr_denoised_interm, ptr_noise_interm);

                luma_weak_filter_avx2_intrin(second_top,
                                             second_curr,
                                             second_bottom,
                                             second_curr_prev,
                                             second_curr_next,
                                             ptr_denoised_interm + 32,
                                             ptr_noise_interm + 32);

                top         = curr;
                curr        = bottom;
                second_top  = second_curr;
                second_curr = second_bottom;
            }
        }

        for (; kk + 16 <= pic_width; kk += 16) {
            for (jj = 0; jj < sb_height; jj++) {
                if (sb_origin_y == 0) {
                    if (jj == 0) {
                        top_128 = _mm_loadl_epi64((__m128i *)(ptr_in + kk + jj * stride_in));
                        second_top_128 =
                            _mm_loadl_epi64((__m128i *)(ptr_in + kk + 8 + jj * stride_in));
                        curr_128 = _mm_loadl_epi64((__m128i *)(ptr_in + kk + (1 + jj) * stride_in));
                        second_curr_128 =
                            _mm_loadl_epi64((__m128i *)((ptr_in + kk + 8) + (1 + jj) * stride_in));
                        _mm_storel_epi64((__m128i *)(ptr_denoised + kk), top_128);
                        _mm_storel_epi64((__m128i *)(ptr_denoised + kk + 8), second_top_128);
                        _mm_storel_epi64((__m128i *)(ptr_noise + kk), _mm_setzero_si128());
                        _mm_storel_epi64((__m128i *)(ptr_noise + kk + 8), _mm_setzero_si128());
                    }
                    curr_prev_128 =
                        _mm_loadl_epi64((__m128i *)(ptr_in - 1 + kk + ((1 + jj) * stride_in)));
                    curr_next_128 =
                        _mm_loadl_epi64((__m128i *)(ptr_in + 1 + kk + ((1 + jj) * stride_in)));
                    second_curr_prev_128 = _mm_loadl_epi64(
                        (__m128i *)((ptr_in + kk + 8) - 1 + ((1 + jj) * stride_in)));
                    second_curr_next_128 = _mm_loadl_epi64(
                        (__m128i *)((ptr_in + kk + 8) + 1 + ((1 + jj) * stride_in)));
                    bottom_128 = _mm_loadl_epi64((__m128i *)((ptr_in + kk) + (2 + jj) * stride_in));
                    second_bottom_128 =
                        _mm_loadl_epi64((__m128i *)((ptr_in + kk + 8) + (2 + jj) * stride_in));
                    ptr_denoised_interm = ptr_denoised + kk + ((1 + jj) * stride_out);
                    ptr_noise_interm    = ptr_noise + kk + ((1 + jj) * stride_out);
                } else {
                    if (jj == 0) {
                        top_128 =
                            _mm_loadl_epi64((__m128i *)(ptr_in + kk + jj * stride_in - stride_in));
                        second_top_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in + kk + 8 + jj * stride_in - stride_in));
                        curr_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in + kk + (1 + jj) * stride_in - stride_in));
                        second_curr_128 = _mm_loadl_epi64(
                            (__m128i *)((ptr_in + kk + 8) + (1 + jj) * stride_in - stride_in));
                    }
                    curr_prev_128 = _mm_loadl_epi64(
                        (__m128i *)(ptr_in - 1 + kk + ((1 + jj) * stride_in - stride_in)));
                    curr_next_128 = _mm_loadl_epi64(
                        (__m128i *)(ptr_in + 1 + kk + ((1 + jj) * stride_in - stride_in)));
                    second_curr_prev_128 = _mm_loadl_epi64(
                        (__m128i *)((ptr_in + kk + 8) - 1 + ((1 + jj) * stride_in - stride_in)));
                    second_curr_next_128 = _mm_loadl_epi64(
                        (__m128i *)((ptr_in + kk + 8) + 1 + ((1 + jj) * stride_in - stride_in)));
                    bottom_128 = _mm_loadl_epi64(
                        (__m128i *)((ptr_in + kk) + (2 + jj) * stride_in - stride_in));
                    second_bottom_128 = _mm_loadl_epi64(
                        (__m128i *)((ptr_in + kk + 8) + (2 + jj) * stride_in - stride_in));
                    ptr_denoised_interm = ptr_denoised + kk + ((1 + jj) * stride_out - stride_out);
                    ptr_noise_interm    = ptr_noise + kk + jj * stride_out;
                }

                luma_weak_filter_128_avx2_intrin(top_128,
                                                 curr_128,
                                                 bottom_128,
                                                 curr_prev_128,
                                                 curr_next_128,
                                                 ptr_denoised_interm,
                                                 ptr_noise_interm);

                luma_weak_filter_128_avx2_intrin(second_top_128,
                                                 second_curr_128,
                                                 second_bottom_128,
                                                 second_curr_prev_128,
                                                 second_curr_next_128,
                                                 ptr_denoised_interm + 8,
                                                 ptr_noise_interm + 8);

                top_128         = curr_128;
                curr_128        = bottom_128;
                second_top_128  = second_curr_128;
                second_curr_128 = second_bottom_128;
            }
        }

        sb_height = MIN(BLOCK_SIZE_64, pic_height - sb_origin_y);

        for (jj = 0; jj < sb_height; jj++) {
            for (ii = 0; ii < pic_width; ii++) {
                if (!((jj < sb_height - 1 || sb_origin_y + sb_height < pic_height) && ii > 0 &&
                      ii < pic_width - 1)) {
                    ptr_denoised[ii + jj * stride_out] = ptr_in[ii + jj * stride_in];
                    ptr_noise[ii + jj * stride_out]    = 0;
                }
            }
        }
    }
}

void noise_extract_luma_weak_sb_avx2_intrin(EbPictureBufferDesc *input_picture_ptr,
                                            EbPictureBufferDesc *denoised_picture_ptr,
                                            EbPictureBufferDesc *noise_picture_ptr,
                                            uint32_t sb_origin_y, uint32_t sb_origin_x) {
    uint32_t ii, jj;
    uint32_t pic_height, sb_height;
    uint32_t pic_width, sb_width;
    uint32_t input_origin_index;
    uint32_t input_origin_index_pad;
    uint32_t noise_origin_index;

    uint8_t *ptr_in;
    uint32_t stride_in;
    uint8_t *ptr_denoised, *ptr_denoised_interm;

    uint8_t *ptr_noise, *ptr_noise_interm;
    uint32_t stride_out;

    __m256i top, curr, bottom, curr_prev, curr_next, second_top, second_curr, second_bottom,
        second_curr_prev, second_curr_next;
    (void)sb_origin_x;

    //Luma
    {
        pic_height = input_picture_ptr->height;
        pic_width  = input_picture_ptr->width;
        sb_height  = MIN(BLOCK_SIZE_64, pic_height - sb_origin_y);
        sb_width   = MIN(BLOCK_SIZE_64, pic_width - sb_origin_x);
        sb_height  = ((sb_origin_y + BLOCK_SIZE_64 >= pic_height) || (sb_origin_y == 0))
                        ? sb_height - 1
                        : sb_height;
        stride_in = input_picture_ptr->stride_y;
        input_origin_index =
            input_picture_ptr->origin_x + sb_origin_x +
            (input_picture_ptr->origin_y + sb_origin_y) * input_picture_ptr->stride_y;
        ptr_in = &(input_picture_ptr->buffer_y[input_origin_index]);

        input_origin_index_pad =
            denoised_picture_ptr->origin_x + sb_origin_x +
            (denoised_picture_ptr->origin_y + sb_origin_y) * denoised_picture_ptr->stride_y;
        stride_out          = denoised_picture_ptr->stride_y;
        ptr_denoised        = &(denoised_picture_ptr->buffer_y[input_origin_index_pad]);
        ptr_denoised_interm = ptr_denoised;

        noise_origin_index = noise_picture_ptr->origin_x + sb_origin_x +
                             noise_picture_ptr->origin_y * noise_picture_ptr->stride_y;
        ptr_noise        = &(noise_picture_ptr->buffer_y[noise_origin_index]);
        ptr_noise_interm = ptr_noise;

        ////Luma
        //a = (p[1] +
        //    p[0 + stride] + 4 * p[1 + stride] + p[2 + stride] +
        //    p[1 + 2 * stride]) / 8;

        top = curr = second_top = second_curr = _mm256_setzero_si256();

        //for (kk = 0; kk + BLOCK_SIZE_64 <= pic_width; kk += BLOCK_SIZE_64)
        {
            for (jj = 0; jj < sb_height; jj++) {
                if (sb_origin_y == 0) {
                    if (jj == 0) {
                        top        = _mm256_loadu_si256((__m256i *)(ptr_in + jj * stride_in));
                        second_top = _mm256_loadu_si256((__m256i *)(ptr_in + 32 + jj * stride_in));
                        curr       = _mm256_loadu_si256((__m256i *)(ptr_in + (1 + jj) * stride_in));
                        second_curr =
                            _mm256_loadu_si256((__m256i *)((ptr_in + 32) + (1 + jj) * stride_in));
                        _mm256_storeu_si256((__m256i *)(ptr_denoised), top);
                        _mm256_storeu_si256((__m256i *)(ptr_denoised + 32), second_top);
                        _mm256_storeu_si256((__m256i *)(ptr_noise), _mm256_setzero_si256());
                        _mm256_storeu_si256((__m256i *)(ptr_noise + 32), _mm256_setzero_si256());
                    }
                    curr_prev =
                        _mm256_loadu_si256((__m256i *)(ptr_in - 1 + ((1 + jj) * stride_in)));
                    curr_next =
                        _mm256_loadu_si256((__m256i *)(ptr_in + 1 + ((1 + jj) * stride_in)));
                    second_curr_prev =
                        _mm256_loadu_si256((__m256i *)((ptr_in + 32) - 1 + ((1 + jj) * stride_in)));
                    second_curr_next =
                        _mm256_loadu_si256((__m256i *)((ptr_in + 32) + 1 + ((1 + jj) * stride_in)));
                    bottom = _mm256_loadu_si256((__m256i *)((ptr_in) + (2 + jj) * stride_in));
                    second_bottom =
                        _mm256_loadu_si256((__m256i *)((ptr_in + 32) + (2 + jj) * stride_in));
                    ptr_denoised_interm = ptr_denoised + ((1 + jj) * stride_out);
                    ptr_noise_interm    = ptr_noise + ((1 + jj) * stride_out);
                } else {
                    if (jj == 0) {
                        top = _mm256_loadu_si256((__m256i *)(ptr_in + jj * stride_in - stride_in));
                        second_top = _mm256_loadu_si256(
                            (__m256i *)(ptr_in + 32 + jj * stride_in - stride_in));
                        curr = _mm256_loadu_si256(
                            (__m256i *)(ptr_in + (1 + jj) * stride_in - stride_in));
                        second_curr = _mm256_loadu_si256(
                            (__m256i *)((ptr_in + 32) + (1 + jj) * stride_in - stride_in));
                    }
                    curr_prev = _mm256_loadu_si256(
                        (__m256i *)(ptr_in - 1 + ((1 + jj) * stride_in - stride_in)));
                    curr_next = _mm256_loadu_si256(
                        (__m256i *)(ptr_in + 1 + ((1 + jj) * stride_in - stride_in)));
                    second_curr_prev = _mm256_loadu_si256(
                        (__m256i *)((ptr_in + 32) - 1 + ((1 + jj) * stride_in - stride_in)));
                    second_curr_next = _mm256_loadu_si256(
                        (__m256i *)((ptr_in + 32) + 1 + ((1 + jj) * stride_in - stride_in)));
                    bottom = _mm256_loadu_si256(
                        (__m256i *)((ptr_in) + (2 + jj) * stride_in - stride_in));
                    second_bottom = _mm256_loadu_si256(
                        (__m256i *)((ptr_in + 32) + (2 + jj) * stride_in - stride_in));
                    ptr_denoised_interm = ptr_denoised + ((1 + jj) * stride_out - stride_out);
                    ptr_noise_interm    = ptr_noise + jj * stride_out;
                }

                luma_weak_filter_avx2_intrin(
                    top, curr, bottom, curr_prev, curr_next, ptr_denoised_interm, ptr_noise_interm);

                luma_weak_filter_avx2_intrin(second_top,
                                             second_curr,
                                             second_bottom,
                                             second_curr_prev,
                                             second_curr_next,
                                             ptr_denoised_interm + 32,
                                             ptr_noise_interm + 32);

                top         = curr;
                curr        = bottom;
                second_top  = second_curr;
                second_curr = second_bottom;
            }
        }

        sb_height = MIN(BLOCK_SIZE_64, pic_height - sb_origin_y);

        for (jj = 0; jj < sb_height; jj++) {
            for (ii = 0; ii < sb_width; ii++) {
                if (!((jj > 0 || sb_origin_y > 0) &&
                      (jj < sb_height - 1 || sb_origin_y + sb_height < pic_height) &&
                      (ii > 0 || sb_origin_x > 0) && (ii + sb_origin_x) < pic_width - 1)) {
                    ptr_denoised[ii + jj * stride_out] = ptr_in[ii + jj * stride_in];
                    ptr_noise[ii + jj * stride_out]    = 0;
                }
            }
        }
    }
}
/*******************************************
* noise_extract_luma_strong
*  strong filter Luma.
*******************************************/
void noise_extract_luma_strong_avx2_intrin(EbPictureBufferDesc *input_picture_ptr,
                                           EbPictureBufferDesc *denoised_picture_ptr,
                                           uint32_t sb_origin_y, uint32_t sb_origin_x) {
    uint32_t ii, jj, kk;
    uint32_t pic_height, sb_height;
    uint32_t pic_width;
    uint32_t input_origin_index;
    uint32_t input_origin_index_pad;

    uint8_t *ptr_in;
    uint32_t stride_in;
    uint8_t *ptr_denoised, *ptr_denoised_interm;

    uint32_t stride_out;
    __m256i  top, curr, bottom, curr_prev, curr_next, top_prev, top_next, bottom_prev, bottom_next,
        second_top, second_curr, second_curr_prev, second_curr_next, second_bottom, second_top_prev,
        second_top_next, second_bottom_prev, second_bottom_next;

    __m128i top_128, curr_128, bottom_128, curr_prev_128, curr_next_128, top_prev_128, top_next_128,
        bottom_prev_128, bottom_next_128, second_top_128, second_curr_128, second_curr_prev_128,
        second_curr_next_128, second_bottom_128, second_top_prev_128, second_top_next_128,
        second_bottom_prev_128, second_bottom_next_128;

    uint32_t idx = (sb_origin_x + BLOCK_SIZE_64 > input_picture_ptr->width) ? sb_origin_x : 0;
    //Luma
    {
        pic_height = input_picture_ptr->height;
        pic_width  = input_picture_ptr->width;
        sb_height  = MIN(BLOCK_SIZE_64, pic_height - sb_origin_y);

        sb_height = ((sb_origin_y + BLOCK_SIZE_64 >= pic_height) || (sb_origin_y == 0))
                        ? sb_height - 1
                        : sb_height;
        stride_in = input_picture_ptr->stride_y;
        input_origin_index =
            input_picture_ptr->origin_x +
            (input_picture_ptr->origin_y + sb_origin_y) * input_picture_ptr->stride_y;
        ptr_in = &(input_picture_ptr->buffer_y[input_origin_index]);

        input_origin_index_pad =
            denoised_picture_ptr->origin_x +
            (denoised_picture_ptr->origin_y + sb_origin_y) * denoised_picture_ptr->stride_y;
        stride_out          = denoised_picture_ptr->stride_y;
        ptr_denoised        = &(denoised_picture_ptr->buffer_y[input_origin_index_pad]);
        ptr_denoised_interm = ptr_denoised;

        top = curr = second_top = second_curr = top_next = top_prev = curr_next = curr_prev =
            second_curr_prev = second_curr_next = second_top_prev = second_top_next =
                _mm256_setzero_si256();
        for (kk = idx; kk + BLOCK_SIZE_64 <= pic_width; kk += BLOCK_SIZE_64) {
            for (jj = 0; jj < sb_height; jj++) {
                if (sb_origin_y == 0) {
                    if (jj == 0) {
                        top = _mm256_loadu_si256((__m256i *)(ptr_in + kk + jj * stride_in));
                        second_top =
                            _mm256_loadu_si256((__m256i *)(ptr_in + kk + 32 + jj * stride_in));

                        curr = _mm256_loadu_si256((__m256i *)(ptr_in + kk + (1 + jj) * stride_in));
                        second_curr = _mm256_loadu_si256(
                            (__m256i *)((ptr_in + kk + 32) + (1 + jj) * stride_in));

                        top_prev =
                            _mm256_loadu_si256((__m256i *)(ptr_in - 1 + kk + ((jj)*stride_in)));
                        second_top_prev = _mm256_loadu_si256(
                            (__m256i *)(ptr_in - 1 + kk + 32 + ((jj)*stride_in)));

                        top_next =
                            _mm256_loadu_si256((__m256i *)(ptr_in + 1 + kk + ((jj)*stride_in)));
                        second_top_next = _mm256_loadu_si256(
                            (__m256i *)(ptr_in + 1 + kk + 32 + ((jj)*stride_in)));

                        curr_prev = _mm256_loadu_si256(
                            (__m256i *)(ptr_in - 1 + kk + ((1 + jj) * stride_in)));
                        second_curr_prev = _mm256_loadu_si256(
                            (__m256i *)(ptr_in - 1 + kk + 32 + ((1 + jj) * stride_in)));

                        curr_next = _mm256_loadu_si256(
                            (__m256i *)(ptr_in + 1 + kk + ((1 + jj) * stride_in)));
                        second_curr_next = _mm256_loadu_si256(
                            (__m256i *)(ptr_in + 1 + kk + 32 + ((1 + jj) * stride_in)));

                        _mm256_storeu_si256((__m256i *)(ptr_denoised + kk), top);
                        _mm256_storeu_si256((__m256i *)(ptr_denoised + kk + 32), second_top);
                    }
                    bottom_prev =
                        _mm256_loadu_si256((__m256i *)(ptr_in - 1 + kk + ((2 + jj) * stride_in)));
                    second_bottom_prev = _mm256_loadu_si256(
                        (__m256i *)(ptr_in - 1 + kk + 32 + ((2 + jj) * stride_in)));

                    bottom_next =
                        _mm256_loadu_si256((__m256i *)(ptr_in + 1 + kk + ((2 + jj) * stride_in)));
                    second_bottom_next = _mm256_loadu_si256(
                        (__m256i *)(ptr_in + 1 + kk + 32 + ((2 + jj) * stride_in)));

                    bottom = _mm256_loadu_si256((__m256i *)((ptr_in + kk) + (2 + jj) * stride_in));
                    second_bottom =
                        _mm256_loadu_si256((__m256i *)((ptr_in + kk + 32) + (2 + jj) * stride_in));
                    ptr_denoised_interm = ptr_denoised + kk + ((1 + jj) * stride_out);
                } else {
                    if (jj == 0) {
                        top = _mm256_loadu_si256(
                            (__m256i *)(ptr_in + kk + jj * stride_in - stride_in));
                        second_top = _mm256_loadu_si256(
                            (__m256i *)(ptr_in + kk + 32 + jj * stride_in - stride_in));
                        curr = _mm256_loadu_si256(
                            (__m256i *)(ptr_in + kk + (1 + jj) * stride_in - stride_in));
                        second_curr = _mm256_loadu_si256(
                            (__m256i *)((ptr_in + kk + 32) + (1 + jj) * stride_in - stride_in));
                        top_prev = _mm256_loadu_si256(
                            (__m256i *)(ptr_in - 1 + kk + ((jj)*stride_in) - stride_in));
                        second_top_prev = _mm256_loadu_si256(
                            (__m256i *)(ptr_in - 1 + kk + 32 + ((jj)*stride_in) - stride_in));

                        top_next = _mm256_loadu_si256(
                            (__m256i *)(ptr_in + 1 + kk + ((jj)*stride_in) - stride_in));
                        second_top_next = _mm256_loadu_si256(
                            (__m256i *)(ptr_in + 1 + kk + 32 + ((jj)*stride_in) - stride_in));

                        curr_prev = _mm256_loadu_si256(
                            (__m256i *)(ptr_in - 1 + kk + ((1 + jj) * stride_in - stride_in)));
                        second_curr_prev = _mm256_loadu_si256(
                            (__m256i *)(ptr_in - 1 + kk + 32 + ((1 + jj) * stride_in - stride_in)));

                        curr_next = _mm256_loadu_si256(
                            (__m256i *)(ptr_in + 1 + kk + ((1 + jj) * stride_in - stride_in)));
                        second_curr_next = _mm256_loadu_si256(
                            (__m256i *)(ptr_in + 1 + kk + 32 + ((1 + jj) * stride_in - stride_in)));
                    }
                    bottom_prev = _mm256_loadu_si256(
                        (__m256i *)(ptr_in - 1 + kk + ((2 + jj) * stride_in) - stride_in));
                    second_bottom_prev = _mm256_loadu_si256(
                        (__m256i *)(ptr_in - 1 + kk + 32 + ((2 + jj) * stride_in - stride_in)));

                    bottom_next = _mm256_loadu_si256(
                        (__m256i *)(ptr_in + 1 + kk + ((2 + jj) * stride_in) - stride_in));
                    second_bottom_next = _mm256_loadu_si256(
                        (__m256i *)(ptr_in + 1 + kk + 32 + ((2 + jj) * stride_in - stride_in)));

                    bottom = _mm256_loadu_si256(
                        (__m256i *)((ptr_in + kk) + (2 + jj) * stride_in - stride_in));
                    second_bottom = _mm256_loadu_si256(
                        (__m256i *)((ptr_in + kk + 32) + (2 + jj) * stride_in - stride_in));

                    ptr_denoised_interm = ptr_denoised + kk + ((1 + jj) * stride_out - stride_out);
                }

                chroma_weak_luma_strong_filter_avx2_intrin(top,
                                                           curr,
                                                           bottom,
                                                           curr_prev,
                                                           curr_next,
                                                           top_prev,
                                                           top_next,
                                                           bottom_prev,
                                                           bottom_next,
                                                           ptr_denoised_interm);

                chroma_weak_luma_strong_filter_avx2_intrin(second_top,
                                                           second_curr,
                                                           second_bottom,
                                                           second_curr_prev,
                                                           second_curr_next,
                                                           second_top_prev,
                                                           second_top_next,
                                                           second_bottom_prev,
                                                           second_bottom_next,
                                                           ptr_denoised_interm + 32);

                top              = curr;
                curr             = bottom;
                top_prev         = curr_prev;
                top_next         = curr_next;
                curr_prev        = bottom_prev;
                curr_next        = bottom_next;
                second_top       = second_curr;
                second_curr      = second_bottom;
                second_top_prev    = second_curr_prev;
                second_top_next    = second_curr_next;
                second_curr_prev = second_bottom_prev;
                second_curr_next = second_bottom_next;
            }
        }

        for (; kk + 16 <= pic_width; kk += 16) {
            for (jj = 0; jj < sb_height; jj++) {
                if (sb_origin_y == 0) {
                    if (jj == 0) {
                        top_128 = _mm_loadl_epi64((__m128i *)(ptr_in + kk + jj * stride_in));
                        second_top_128 =
                            _mm_loadl_epi64((__m128i *)(ptr_in + kk + 8 + jj * stride_in));

                        curr_128 = _mm_loadl_epi64((__m128i *)(ptr_in + kk + (1 + jj) * stride_in));
                        second_curr_128 =
                            _mm_loadl_epi64((__m128i *)((ptr_in + kk + 8) + (1 + jj) * stride_in));

                        top_prev_128 =
                            _mm_loadl_epi64((__m128i *)(ptr_in - 1 + kk + ((jj)*stride_in)));
                        second_top_prev_128 =
                            _mm_loadl_epi64((__m128i *)(ptr_in - 1 + kk + 8 + ((jj)*stride_in)));

                        top_next_128 =
                            _mm_loadl_epi64((__m128i *)(ptr_in + 1 + kk + ((jj)*stride_in)));
                        second_top_next_128 =
                            _mm_loadl_epi64((__m128i *)(ptr_in + 1 + kk + 8 + ((jj)*stride_in)));

                        curr_prev_128 =
                            _mm_loadl_epi64((__m128i *)(ptr_in - 1 + kk + ((1 + jj) * stride_in)));
                        second_curr_prev_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in - 1 + kk + 8 + ((1 + jj) * stride_in)));

                        curr_next_128 =
                            _mm_loadl_epi64((__m128i *)(ptr_in + 1 + kk + ((1 + jj) * stride_in)));
                        second_curr_next_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in + 1 + kk + 8 + ((1 + jj) * stride_in)));

                        _mm_storel_epi64((__m128i *)(ptr_denoised + kk), top_128);
                        _mm_storel_epi64((__m128i *)(ptr_denoised + kk + 8), second_top_128);
                    }
                    bottom_prev_128 =
                        _mm_loadl_epi64((__m128i *)(ptr_in - 1 + kk + ((2 + jj) * stride_in)));
                    second_bottom_prev_128 =
                        _mm_loadl_epi64((__m128i *)(ptr_in - 1 + kk + 8 + ((2 + jj) * stride_in)));

                    bottom_next_128 =
                        _mm_loadl_epi64((__m128i *)(ptr_in + 1 + kk + ((2 + jj) * stride_in)));
                    second_bottom_next_128 =
                        _mm_loadl_epi64((__m128i *)(ptr_in + 1 + kk + 8 + ((2 + jj) * stride_in)));

                    bottom_128 = _mm_loadl_epi64((__m128i *)((ptr_in + kk) + (2 + jj) * stride_in));
                    second_bottom_128 =
                        _mm_loadl_epi64((__m128i *)((ptr_in + kk + 8) + (2 + jj) * stride_in));
                    ptr_denoised_interm = ptr_denoised + kk + ((1 + jj) * stride_out);
                } else {
                    if (jj == 0) {
                        top_128 =
                            _mm_loadl_epi64((__m128i *)(ptr_in + kk + jj * stride_in - stride_in));
                        second_top_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in + kk + 8 + jj * stride_in - stride_in));
                        curr_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in + kk + (1 + jj) * stride_in - stride_in));
                        second_curr_128 = _mm_loadl_epi64(
                            (__m128i *)((ptr_in + kk + 8) + (1 + jj) * stride_in - stride_in));
                        top_prev_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in - 1 + kk + ((jj)*stride_in) - stride_in));
                        second_top_prev_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in - 1 + kk + 8 + ((jj)*stride_in) - stride_in));

                        top_next_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in + 1 + kk + ((jj)*stride_in) - stride_in));
                        second_top_next_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in + 1 + kk + 8 + ((jj)*stride_in) - stride_in));

                        curr_prev_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in - 1 + kk + ((1 + jj) * stride_in - stride_in)));
                        second_curr_prev_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in - 1 + kk + 8 + ((1 + jj) * stride_in - stride_in)));

                        curr_next_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in + 1 + kk + ((1 + jj) * stride_in - stride_in)));
                        second_curr_next_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in + 1 + kk + 8 + ((1 + jj) * stride_in - stride_in)));
                    }
                    bottom_prev_128 = _mm_loadl_epi64(
                        (__m128i *)(ptr_in - 1 + kk + ((2 + jj) * stride_in) - stride_in));
                    second_bottom_prev_128 = _mm_loadl_epi64(
                        (__m128i *)(ptr_in - 1 + kk + 8 + ((2 + jj) * stride_in - stride_in)));

                    bottom_next_128 = _mm_loadl_epi64(
                        (__m128i *)(ptr_in + 1 + kk + ((2 + jj) * stride_in) - stride_in));
                    second_bottom_next_128 = _mm_loadl_epi64(
                        (__m128i *)(ptr_in + 1 + kk + 8 + ((2 + jj) * stride_in - stride_in)));

                    bottom_128 = _mm_loadl_epi64(
                        (__m128i *)((ptr_in + kk) + (2 + jj) * stride_in - stride_in));
                    second_bottom_128 = _mm_loadl_epi64(
                        (__m128i *)((ptr_in + kk + 8) + (2 + jj) * stride_in - stride_in));

                    ptr_denoised_interm = ptr_denoised + kk + ((1 + jj) * stride_out - stride_out);
                }

                chroma_weak_luma_strong_filter_128_avx2_intrin(top_128,
                                                               curr_128,
                                                               bottom_128,
                                                               curr_prev_128,
                                                               curr_next_128,
                                                               top_prev_128,
                                                               top_next_128,
                                                               bottom_prev_128,
                                                               bottom_next_128,
                                                               ptr_denoised_interm);

                chroma_weak_luma_strong_filter_128_avx2_intrin(second_top_128,
                                                               second_curr_128,
                                                               second_bottom_128,
                                                               second_curr_prev_128,
                                                               second_curr_next_128,
                                                               second_top_prev_128,
                                                               second_top_next_128,
                                                               second_bottom_prev_128,
                                                               second_bottom_next_128,
                                                               ptr_denoised_interm + 8);

                top_128              = curr_128;
                curr_128             = bottom_128;
                top_prev_128         = curr_prev_128;
                top_next_128         = curr_next_128;
                curr_prev_128        = bottom_prev_128;
                curr_next_128        = bottom_next_128;
                second_top_128       = second_curr_128;
                second_curr_128      = second_bottom_128;
                second_top_prev_128  = second_curr_prev_128;
                second_top_next_128  = second_curr_next_128;
                second_curr_prev_128 = second_bottom_prev_128;
                second_curr_next_128 = second_bottom_next_128;
            }
        }

        sb_height = MIN(BLOCK_SIZE_64, pic_height - sb_origin_y);

        for (jj = 0; jj < sb_height; jj++) {
            for (ii = 0; ii < pic_width; ii++) {
                if (!((jj < sb_height - 1 || sb_origin_y + sb_height < pic_height) && ii > 0 &&
                      ii < pic_width - 1))
                    ptr_denoised[ii + jj * stride_out] = ptr_in[ii + jj * stride_in];
            }
        }
    }
}

/*******************************************
* noise_extract_chroma_strong
*  strong filter chroma.
*******************************************/
void noise_extract_chroma_strong_avx2_intrin(EbPictureBufferDesc *input_picture_ptr,
                                             EbPictureBufferDesc *denoised_picture_ptr,
                                             uint32_t sb_origin_y, uint32_t sb_origin_x) {
    uint32_t ii, jj, kk;
    uint32_t pic_height, sb_height;
    uint32_t pic_width;
    uint32_t input_origin_index;
    uint32_t input_origin_index_pad;

    uint8_t *ptr_in, *ptr_in_cr;
    uint32_t stride_in, stride_in_cr;
    uint8_t *ptr_denoised, *ptr_denoised_interm, *ptr_denoised_cr, *ptr_denoised_interm_cr;

    uint32_t stride_out, stride_out_cr;
    __m256i  top, curr, bottom, curr_prev, curr_next, top_prev, top_next, bottom_prev, bottom_next,
        top_cr, curr_cr, bottom_cr, curr_prev_cr, curr_next_cr, top_prev_cr, top_next_cr,
        bottom_prev_cr, bottom_next_cr;

    __m128i top_128, curr_128, bottom_128, curr_prev_128, curr_next_128, top_prev_128, top_next_128,
        bottom_prev_128, bottom_next_128, top_cr_128, curr_cr_128, bottom_cr_128, curr_prev_cr_128,
        curr_next_cr_128, top_prev_cr_128, top_next_cr_128, bottom_prev_cr_128, bottom_next_cr_128;
    uint32_t idx = (sb_origin_x + BLOCK_SIZE_64 > input_picture_ptr->width) ? sb_origin_x : 0;
    {
        pic_height = input_picture_ptr->height / 2;
        pic_width  = input_picture_ptr->width / 2;
        sb_height  = MIN(BLOCK_SIZE_64 / 2, pic_height - sb_origin_y);

        sb_height = ((sb_origin_y + BLOCK_SIZE_64 / 2 >= pic_height) || (sb_origin_y == 0))
                        ? sb_height - 1
                        : sb_height;

        stride_in = input_picture_ptr->stride_cb;
        input_origin_index =
            input_picture_ptr->origin_x / 2 +
            (input_picture_ptr->origin_y / 2 + sb_origin_y) * input_picture_ptr->stride_cb;
        ptr_in = &(input_picture_ptr->buffer_cb[input_origin_index]);

        input_origin_index_pad =
            denoised_picture_ptr->origin_x / 2 +
            (denoised_picture_ptr->origin_y / 2 + sb_origin_y) * denoised_picture_ptr->stride_cb;
        stride_out          = denoised_picture_ptr->stride_cb;
        ptr_denoised        = &(denoised_picture_ptr->buffer_cb[input_origin_index_pad]);
        ptr_denoised_interm = ptr_denoised;

        stride_in_cr = input_picture_ptr->stride_cr;
        input_origin_index =
            input_picture_ptr->origin_x / 2 +
            (input_picture_ptr->origin_y / 2 + sb_origin_y) * input_picture_ptr->stride_cr;
        ptr_in_cr = &(input_picture_ptr->buffer_cr[input_origin_index]);

        input_origin_index_pad =
            denoised_picture_ptr->origin_x / 2 +
            (denoised_picture_ptr->origin_y / 2 + sb_origin_y) * denoised_picture_ptr->stride_cr;
        stride_out_cr          = denoised_picture_ptr->stride_cr;
        ptr_denoised_cr        = &(denoised_picture_ptr->buffer_cr[input_origin_index_pad]);
        ptr_denoised_interm_cr = ptr_denoised_cr;
        ////Chroma
        //a = (4 * p[0] + 4 * p[1] + 4 * p[2] +
        //    4 * p[0 + stride] + 4 * p[1 + stride] + 4 * p[2 + stride] +
        //    4 * p[0 + 2 * stride] + 4 * p[1 + 2 * stride] + 4 * p[2 + 2 * stride]) / 36;

        top = curr = top_next = top_prev = curr_next = curr_prev = top_cr = curr_cr = top_next_cr =
            top_prev_cr = curr_next_cr = curr_prev_cr = _mm256_setzero_si256();

        for (kk = idx; kk + BLOCK_SIZE_64 / 2 <= pic_width; kk += BLOCK_SIZE_64 / 2) {
            for (jj = 0; jj < sb_height; jj++) {
                if (sb_origin_y == 0) {
                    if (jj == 0) {
                        top  = _mm256_loadu_si256((__m256i *)(ptr_in + kk + jj * stride_in));
                        curr = _mm256_loadu_si256((__m256i *)(ptr_in + kk + (1 + jj) * stride_in));
                        top_prev =
                            _mm256_loadu_si256((__m256i *)(ptr_in - 1 + kk + ((jj)*stride_in)));
                        top_next =
                            _mm256_loadu_si256((__m256i *)(ptr_in + 1 + kk + ((jj)*stride_in)));
                        curr_prev = _mm256_loadu_si256(
                            (__m256i *)(ptr_in - 1 + kk + ((1 + jj) * stride_in)));
                        curr_next = _mm256_loadu_si256(
                            (__m256i *)(ptr_in + 1 + kk + ((1 + jj) * stride_in)));
                        top_cr =
                            _mm256_loadu_si256((__m256i *)(ptr_in_cr + kk + jj * stride_in_cr));
                        curr_cr = _mm256_loadu_si256(
                            (__m256i *)(ptr_in_cr + kk + (1 + jj) * stride_in_cr));
                        top_prev_cr = _mm256_loadu_si256(
                            (__m256i *)(ptr_in_cr - 1 + kk + ((jj)*stride_in_cr)));
                        top_next_cr = _mm256_loadu_si256(
                            (__m256i *)(ptr_in_cr + 1 + kk + ((jj)*stride_in_cr)));
                        curr_prev_cr = _mm256_loadu_si256(
                            (__m256i *)(ptr_in_cr - 1 + kk + ((1 + jj) * stride_in_cr)));
                        curr_next_cr = _mm256_loadu_si256(
                            (__m256i *)(ptr_in_cr + 1 + kk + ((1 + jj) * stride_in_cr)));
                        _mm256_storeu_si256((__m256i *)(ptr_denoised + kk), top);
                        _mm256_storeu_si256((__m256i *)(ptr_denoised_cr + kk), top_cr);
                    }
                    bottom_prev =
                        _mm256_loadu_si256((__m256i *)(ptr_in - 1 + kk + ((2 + jj) * stride_in)));
                    bottom_next =
                        _mm256_loadu_si256((__m256i *)(ptr_in + 1 + kk + ((2 + jj) * stride_in)));
                    bottom = _mm256_loadu_si256((__m256i *)((ptr_in + kk) + (2 + jj) * stride_in));
                    bottom_prev_cr = _mm256_loadu_si256(
                        (__m256i *)(ptr_in_cr - 1 + kk + ((2 + jj) * stride_in_cr)));
                    bottom_next_cr = _mm256_loadu_si256(
                        (__m256i *)(ptr_in_cr + 1 + kk + ((2 + jj) * stride_in_cr)));
                    bottom_cr =
                        _mm256_loadu_si256((__m256i *)((ptr_in_cr + kk) + (2 + jj) * stride_in_cr));
                    ptr_denoised_interm    = ptr_denoised + kk + ((1 + jj) * stride_out);
                    ptr_denoised_interm_cr = ptr_denoised_cr + kk + ((1 + jj) * stride_out_cr);
                } else {
                    if (jj == 0) {
                        top = _mm256_loadu_si256(
                            (__m256i *)(ptr_in + kk + jj * stride_in - stride_in));
                        curr = _mm256_loadu_si256(
                            (__m256i *)(ptr_in + kk + (1 + jj) * stride_in - stride_in));
                        top_prev = _mm256_loadu_si256(
                            (__m256i *)(ptr_in - 1 + kk + ((jj)*stride_in) - stride_in));
                        top_next = _mm256_loadu_si256(
                            (__m256i *)(ptr_in + 1 + kk + ((jj)*stride_in) - stride_in));
                        curr_prev = _mm256_loadu_si256(
                            (__m256i *)(ptr_in - 1 + kk + ((1 + jj) * stride_in - stride_in)));
                        curr_next = _mm256_loadu_si256(
                            (__m256i *)(ptr_in + 1 + kk + ((1 + jj) * stride_in - stride_in)));
                        top_cr = _mm256_loadu_si256(
                            (__m256i *)(ptr_in_cr + kk + jj * stride_in_cr - stride_in_cr));
                        curr_cr = _mm256_loadu_si256(
                            (__m256i *)(ptr_in_cr + kk + (1 + jj) * stride_in_cr - stride_in_cr));
                        top_prev_cr = _mm256_loadu_si256(
                            (__m256i *)(ptr_in_cr - 1 + kk + ((jj)*stride_in_cr) - stride_in_cr));
                        top_next_cr = _mm256_loadu_si256(
                            (__m256i *)(ptr_in_cr + 1 + kk + ((jj)*stride_in_cr) - stride_in_cr));
                        curr_prev_cr = _mm256_loadu_si256(
                            (__m256i *)(ptr_in_cr - 1 + kk +
                                        ((1 + jj) * stride_in_cr - stride_in_cr)));
                        curr_next_cr = _mm256_loadu_si256(
                            (__m256i *)(ptr_in_cr + 1 + kk +
                                        ((1 + jj) * stride_in_cr - stride_in_cr)));
                    }
                    bottom_prev = _mm256_loadu_si256(
                        (__m256i *)(ptr_in - 1 + kk + ((2 + jj) * stride_in) - stride_in));
                    bottom_next = _mm256_loadu_si256(
                        (__m256i *)(ptr_in + 1 + kk + ((2 + jj) * stride_in) - stride_in));
                    bottom = _mm256_loadu_si256(
                        (__m256i *)((ptr_in + kk) + (2 + jj) * stride_in - stride_in));
                    ptr_denoised_interm = ptr_denoised + kk + ((1 + jj) * stride_out - stride_out);
                    bottom_prev_cr      = _mm256_loadu_si256(
                        (__m256i *)(ptr_in_cr - 1 + kk + ((2 + jj) * stride_in_cr) - stride_in_cr));
                    bottom_next_cr = _mm256_loadu_si256(
                        (__m256i *)(ptr_in_cr + 1 + kk + ((2 + jj) * stride_in_cr) - stride_in_cr));
                    bottom_cr = _mm256_loadu_si256(
                        (__m256i *)((ptr_in_cr + kk) + (2 + jj) * stride_in_cr - stride_in_cr));
                    ptr_denoised_interm_cr =
                        ptr_denoised_cr + kk + ((1 + jj) * stride_out_cr - stride_out_cr);
                }

                chroma_strong_avx2_intrin(top,
                                          curr,
                                          bottom,
                                          curr_prev,
                                          curr_next,
                                          top_prev,
                                          top_next,
                                          bottom_prev,
                                          bottom_next,
                                          ptr_denoised_interm);

                chroma_strong_avx2_intrin(top_cr,
                                          curr_cr,
                                          bottom_cr,
                                          curr_prev_cr,
                                          curr_next_cr,
                                          top_prev_cr,
                                          top_next_cr,
                                          bottom_prev_cr,
                                          bottom_next_cr,
                                          ptr_denoised_interm_cr);

                top          = curr;
                curr         = bottom;
                top_prev     = curr_prev;
                top_next     = curr_next;
                curr_prev    = bottom_prev;
                curr_next    = bottom_next;
                top_cr       = curr_cr;
                curr_cr      = bottom_cr;
                top_prev_cr  = curr_prev_cr;
                top_next_cr  = curr_next_cr;
                curr_prev_cr = bottom_prev_cr;
                curr_next_cr = bottom_next_cr;
            }
        }

        for (; kk + 8 <= pic_width; kk += 8) {
            for (jj = 0; jj < sb_height; jj++) {
                if (sb_origin_y == 0) {
                    if (jj == 0) {
                        top_128  = _mm_loadl_epi64((__m128i *)(ptr_in + kk + jj * stride_in));
                        curr_128 = _mm_loadl_epi64((__m128i *)(ptr_in + kk + (1 + jj) * stride_in));
                        top_prev_128 =
                            _mm_loadl_epi64((__m128i *)(ptr_in - 1 + kk + ((jj)*stride_in)));
                        top_next_128 =
                            _mm_loadl_epi64((__m128i *)(ptr_in + 1 + kk + ((jj)*stride_in)));
                        curr_prev_128 =
                            _mm_loadl_epi64((__m128i *)(ptr_in - 1 + kk + ((1 + jj) * stride_in)));
                        curr_next_128 =
                            _mm_loadl_epi64((__m128i *)(ptr_in + 1 + kk + ((1 + jj) * stride_in)));
                        top_cr_128 =
                            _mm_loadl_epi64((__m128i *)(ptr_in_cr + kk + jj * stride_in_cr));
                        curr_cr_128 =
                            _mm_loadl_epi64((__m128i *)(ptr_in_cr + kk + (1 + jj) * stride_in_cr));
                        top_prev_cr_128 =
                            _mm_loadl_epi64((__m128i *)(ptr_in_cr - 1 + kk + ((jj)*stride_in_cr)));
                        top_next_cr_128 =
                            _mm_loadl_epi64((__m128i *)(ptr_in_cr + 1 + kk + ((jj)*stride_in_cr)));
                        curr_prev_cr_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in_cr - 1 + kk + ((1 + jj) * stride_in_cr)));
                        curr_next_cr_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in_cr + 1 + kk + ((1 + jj) * stride_in_cr)));
                        _mm_storel_epi64((__m128i *)(ptr_denoised + kk), top_128);
                        _mm_storel_epi64((__m128i *)(ptr_denoised_cr + kk), top_cr_128);
                    }
                    bottom_prev_128 =
                        _mm_loadl_epi64((__m128i *)(ptr_in - 1 + kk + ((2 + jj) * stride_in)));
                    bottom_next_128 =
                        _mm_loadl_epi64((__m128i *)(ptr_in + 1 + kk + ((2 + jj) * stride_in)));
                    bottom_128 = _mm_loadl_epi64((__m128i *)((ptr_in + kk) + (2 + jj) * stride_in));
                    bottom_prev_cr_128 = _mm_loadl_epi64(
                        (__m128i *)(ptr_in_cr - 1 + kk + ((2 + jj) * stride_in_cr)));
                    bottom_next_cr_128 = _mm_loadl_epi64(
                        (__m128i *)(ptr_in_cr + 1 + kk + ((2 + jj) * stride_in_cr)));
                    bottom_cr_128 =
                        _mm_loadl_epi64((__m128i *)((ptr_in_cr + kk) + (2 + jj) * stride_in_cr));
                    ptr_denoised_interm    = ptr_denoised + kk + ((1 + jj) * stride_out);
                    ptr_denoised_interm_cr = ptr_denoised_cr + kk + ((1 + jj) * stride_out_cr);
                } else {
                    if (jj == 0) {
                        top_128 =
                            _mm_loadl_epi64((__m128i *)(ptr_in + kk + jj * stride_in - stride_in));
                        curr_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in + kk + (1 + jj) * stride_in - stride_in));
                        top_prev_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in - 1 + kk + ((jj)*stride_in) - stride_in));
                        top_next_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in + 1 + kk + ((jj)*stride_in) - stride_in));
                        curr_prev_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in - 1 + kk + ((1 + jj) * stride_in - stride_in)));
                        curr_next_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in + 1 + kk + ((1 + jj) * stride_in - stride_in)));
                        top_cr_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in_cr + kk + jj * stride_in_cr - stride_in_cr));
                        curr_cr_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in_cr + kk + (1 + jj) * stride_in_cr - stride_in_cr));
                        top_prev_cr_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in_cr - 1 + kk + ((jj)*stride_in_cr) - stride_in_cr));
                        top_next_cr_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in_cr + 1 + kk + ((jj)*stride_in_cr) - stride_in_cr));
                        curr_prev_cr_128 =
                            _mm_loadl_epi64((__m128i *)(ptr_in_cr - 1 + kk +
                                                        ((1 + jj) * stride_in_cr - stride_in_cr)));
                        curr_next_cr_128 =
                            _mm_loadl_epi64((__m128i *)(ptr_in_cr + 1 + kk +
                                                        ((1 + jj) * stride_in_cr - stride_in_cr)));
                    }
                    bottom_prev_128 = _mm_loadl_epi64(
                        (__m128i *)(ptr_in - 1 + kk + ((2 + jj) * stride_in) - stride_in));
                    bottom_next_128 = _mm_loadl_epi64(
                        (__m128i *)(ptr_in + 1 + kk + ((2 + jj) * stride_in) - stride_in));
                    bottom_128 = _mm_loadl_epi64(
                        (__m128i *)((ptr_in + kk) + (2 + jj) * stride_in - stride_in));
                    ptr_denoised_interm = ptr_denoised + kk + ((1 + jj) * stride_out - stride_out);
                    bottom_prev_cr_128  = _mm_loadl_epi64(
                        (__m128i *)(ptr_in_cr - 1 + kk + ((2 + jj) * stride_in_cr) - stride_in_cr));
                    bottom_next_cr_128 = _mm_loadl_epi64(
                        (__m128i *)(ptr_in_cr + 1 + kk + ((2 + jj) * stride_in_cr) - stride_in_cr));
                    bottom_cr_128 = _mm_loadl_epi64(
                        (__m128i *)((ptr_in_cr + kk) + (2 + jj) * stride_in_cr - stride_in_cr));
                    ptr_denoised_interm_cr =
                        ptr_denoised_cr + kk + ((1 + jj) * stride_out_cr - stride_out_cr);
                }

                chroma_strong_128_avx2_intrin(top_128,
                                              curr_128,
                                              bottom_128,
                                              curr_prev_128,
                                              curr_next_128,
                                              top_prev_128,
                                              top_next_128,
                                              bottom_prev_128,
                                              bottom_next_128,
                                              ptr_denoised_interm);

                chroma_strong_128_avx2_intrin(top_cr_128,
                                              curr_cr_128,
                                              bottom_cr_128,
                                              curr_prev_cr_128,
                                              curr_next_cr_128,
                                              top_prev_cr_128,
                                              top_next_cr_128,
                                              bottom_prev_cr_128,
                                              bottom_next_cr_128,
                                              ptr_denoised_interm_cr);

                top_128          = curr_128;
                curr_128         = bottom_128;
                top_prev_128     = curr_prev_128;
                top_next_128     = curr_next_128;
                curr_prev_128    = bottom_prev_128;
                curr_next_128    = bottom_next_128;
                top_cr_128       = curr_cr_128;
                curr_cr_128      = bottom_cr_128;
                top_prev_cr_128  = curr_prev_cr_128;
                top_next_cr_128  = curr_next_cr_128;
                curr_prev_cr_128 = bottom_prev_cr_128;
                curr_next_cr_128 = bottom_next_cr_128;
            }
        }

        sb_height = MIN(BLOCK_SIZE_64 / 2, pic_height - sb_origin_y);

        for (jj = 0; jj < sb_height; jj++) {
            for (ii = idx; ii < pic_width; ii++) {
                if (!((jj < sb_height - 1 || (sb_origin_y + sb_height) < pic_height) && ii > 0 &&
                      ii < pic_width - 1)) {
                    ptr_denoised[ii + jj * stride_out]    = ptr_in[ii + jj * stride_in];
                    ptr_denoised_cr[ii + jj * stride_out] = ptr_in_cr[ii + jj * stride_in];
                }
            }
        }
    }
}

/*******************************************
* noise_extract_chroma_weak
*  weak filter chroma.
*******************************************/
void noise_extract_chroma_weak_avx2_intrin(EbPictureBufferDesc *input_picture_ptr,
                                           EbPictureBufferDesc *denoised_picture_ptr,
                                           uint32_t sb_origin_y, uint32_t sb_origin_x) {
    uint32_t ii, jj, kk;
    uint32_t pic_height, sb_height;
    uint32_t pic_width;
    uint32_t input_origin_index;
    uint32_t input_origin_index_pad;

    uint8_t *ptr_in, *ptr_in_cr;
    uint32_t stride_in, stride_in_cr;
    uint8_t *ptr_denoised, *ptr_denoised_interm, *ptr_denoised_cr, *ptr_denoised_interm_cr;

    uint32_t stride_out, stride_out_cr;

    __m256i top, curr, bottom, curr_prev, curr_next, top_prev, top_next, bottom_prev, bottom_next,
        top_cr, curr_cr, bottom_cr, curr_prev_cr, curr_next_cr, top_prev_cr, top_next_cr,
        bottom_prev_cr, bottom_next_cr;

    __m128i top_128, curr_128, bottom_128, curr_prev_128, curr_next_128, top_prev_128, top_next_128,
        bottom_prev_128, bottom_next_128, top_cr_128, curr_cr_128, bottom_cr_128, curr_prev_cr_128,
        curr_next_cr_128, top_prev_cr_128, top_next_cr_128, bottom_prev_cr_128, bottom_next_cr_128;
    uint32_t idx = (sb_origin_x + BLOCK_SIZE_64 > input_picture_ptr->width) ? sb_origin_x : 0;
    ////gaussian matrix(Chroma)
    //a = (1 * p[0] + 2 * p[1] + 1 * p[2] +
    //    2 * p[0 + stride] + 4 * p[1 + stride] + 2 * p[2 + stride] +
    //    1 * p[0 + 2 * stride] + 2 * p[1 + 2 * stride] + 1 * p[2 + 2 * stride]) / 16;

    {
        pic_height = input_picture_ptr->height / 2;
        pic_width  = input_picture_ptr->width / 2;

        sb_height = MIN(BLOCK_SIZE_64 / 2, pic_height - sb_origin_y);

        sb_height = ((sb_origin_y + BLOCK_SIZE_64 / 2 >= pic_height) || (sb_origin_y == 0))
                        ? sb_height - 1
                        : sb_height;
        stride_in = input_picture_ptr->stride_cb;
        input_origin_index =
            input_picture_ptr->origin_x / 2 +
            (input_picture_ptr->origin_y / 2 + sb_origin_y) * input_picture_ptr->stride_cb;
        ptr_in = &(input_picture_ptr->buffer_cb[input_origin_index]);

        input_origin_index_pad =
            denoised_picture_ptr->origin_x / 2 +
            (denoised_picture_ptr->origin_y / 2 + sb_origin_y) * denoised_picture_ptr->stride_cb;
        stride_out          = denoised_picture_ptr->stride_cb;
        ptr_denoised        = &(denoised_picture_ptr->buffer_cb[input_origin_index_pad]);
        ptr_denoised_interm = ptr_denoised;

        stride_in_cr = input_picture_ptr->stride_cr;
        input_origin_index =
            input_picture_ptr->origin_x / 2 +
            (input_picture_ptr->origin_y / 2 + sb_origin_y) * input_picture_ptr->stride_cr;
        ptr_in_cr = &(input_picture_ptr->buffer_cr[input_origin_index]);

        input_origin_index_pad =
            denoised_picture_ptr->origin_x / 2 +
            (denoised_picture_ptr->origin_y / 2 + sb_origin_y) * denoised_picture_ptr->stride_cr;
        stride_out_cr          = denoised_picture_ptr->stride_cr;
        ptr_denoised_cr        = &(denoised_picture_ptr->buffer_cr[input_origin_index_pad]);
        ptr_denoised_interm_cr = ptr_denoised_cr;

        top = curr = top_next = top_prev = curr_next = curr_prev = top_cr = curr_cr = top_next_cr =
            top_prev_cr = curr_next_cr = curr_prev_cr = _mm256_setzero_si256();
        for (kk = idx; kk + BLOCK_SIZE_64 / 2 <= pic_width; kk += BLOCK_SIZE_64 / 2) {
            for (jj = 0; jj < sb_height; jj++) {
                if (sb_origin_y == 0) {
                    if (jj == 0) {
                        top  = _mm256_loadu_si256((__m256i *)(ptr_in + kk + jj * stride_in));
                        curr = _mm256_loadu_si256((__m256i *)(ptr_in + kk + (1 + jj) * stride_in));
                        top_prev =
                            _mm256_loadu_si256((__m256i *)(ptr_in - 1 + kk + ((jj)*stride_in)));
                        top_next =
                            _mm256_loadu_si256((__m256i *)(ptr_in + 1 + kk + ((jj)*stride_in)));
                        curr_prev = _mm256_loadu_si256(
                            (__m256i *)(ptr_in - 1 + kk + ((1 + jj) * stride_in)));
                        curr_next = _mm256_loadu_si256(
                            (__m256i *)(ptr_in + 1 + kk + ((1 + jj) * stride_in)));
                        _mm256_storeu_si256((__m256i *)(ptr_denoised + kk), top);
                        top_cr =
                            _mm256_loadu_si256((__m256i *)(ptr_in_cr + kk + jj * stride_in_cr));
                        curr_cr = _mm256_loadu_si256(
                            (__m256i *)(ptr_in_cr + kk + (1 + jj) * stride_in_cr));
                        top_prev_cr = _mm256_loadu_si256(
                            (__m256i *)(ptr_in_cr - 1 + kk + ((jj)*stride_in_cr)));
                        top_next_cr = _mm256_loadu_si256(
                            (__m256i *)(ptr_in_cr + 1 + kk + ((jj)*stride_in_cr)));
                        curr_prev_cr = _mm256_loadu_si256(
                            (__m256i *)(ptr_in_cr - 1 + kk + ((1 + jj) * stride_in_cr)));
                        curr_next_cr = _mm256_loadu_si256(
                            (__m256i *)(ptr_in_cr + 1 + kk + ((1 + jj) * stride_in_cr)));
                        _mm256_storeu_si256((__m256i *)(ptr_denoised_cr + kk), top_cr);
                    }
                    bottom_prev =
                        _mm256_loadu_si256((__m256i *)(ptr_in - 1 + kk + ((2 + jj) * stride_in)));
                    bottom_next =
                        _mm256_loadu_si256((__m256i *)(ptr_in + 1 + kk + ((2 + jj) * stride_in)));
                    bottom = _mm256_loadu_si256((__m256i *)((ptr_in + kk) + (2 + jj) * stride_in));
                    ptr_denoised_interm = ptr_denoised + kk + ((1 + jj) * stride_out);
                    bottom_prev_cr      = _mm256_loadu_si256(
                        (__m256i *)(ptr_in_cr - 1 + kk + ((2 + jj) * stride_in_cr)));
                    bottom_next_cr = _mm256_loadu_si256(
                        (__m256i *)(ptr_in_cr + 1 + kk + ((2 + jj) * stride_in_cr)));
                    bottom_cr =
                        _mm256_loadu_si256((__m256i *)((ptr_in_cr + kk) + (2 + jj) * stride_in_cr));
                    ptr_denoised_interm_cr = ptr_denoised_cr + kk + ((1 + jj) * stride_out_cr);
                } else {
                    if (jj == 0) {
                        top = _mm256_loadu_si256(
                            (__m256i *)(ptr_in + kk + jj * stride_in - stride_in));
                        curr = _mm256_loadu_si256(
                            (__m256i *)(ptr_in + kk + (1 + jj) * stride_in - stride_in));
                        top_prev = _mm256_loadu_si256(
                            (__m256i *)(ptr_in - 1 + kk + ((jj)*stride_in) - stride_in));
                        top_next = _mm256_loadu_si256(
                            (__m256i *)(ptr_in + 1 + kk + ((jj)*stride_in) - stride_in));
                        curr_prev = _mm256_loadu_si256(
                            (__m256i *)(ptr_in - 1 + kk + ((1 + jj) * stride_in - stride_in)));
                        curr_next = _mm256_loadu_si256(
                            (__m256i *)(ptr_in + 1 + kk + ((1 + jj) * stride_in - stride_in)));
                        top_cr = _mm256_loadu_si256(
                            (__m256i *)(ptr_in_cr + kk + jj * stride_in_cr - stride_in_cr));
                        curr_cr = _mm256_loadu_si256(
                            (__m256i *)(ptr_in_cr + kk + (1 + jj) * stride_in_cr - stride_in_cr));
                        top_prev_cr = _mm256_loadu_si256(
                            (__m256i *)(ptr_in_cr - 1 + kk + ((jj)*stride_in_cr) - stride_in_cr));
                        top_next_cr = _mm256_loadu_si256(
                            (__m256i *)(ptr_in_cr + 1 + kk + ((jj)*stride_in_cr) - stride_in_cr));
                        curr_prev_cr = _mm256_loadu_si256(
                            (__m256i *)(ptr_in_cr - 1 + kk +
                                        ((1 + jj) * stride_in_cr - stride_in_cr)));
                        curr_next_cr = _mm256_loadu_si256(
                            (__m256i *)(ptr_in_cr + 1 + kk +
                                        ((1 + jj) * stride_in_cr - stride_in_cr)));
                    }
                    bottom_prev = _mm256_loadu_si256(
                        (__m256i *)(ptr_in - 1 + kk + ((2 + jj) * stride_in) - stride_in));
                    bottom_next = _mm256_loadu_si256(
                        (__m256i *)(ptr_in + 1 + kk + ((2 + jj) * stride_in) - stride_in));
                    bottom = _mm256_loadu_si256(
                        (__m256i *)((ptr_in + kk) + (2 + jj) * stride_in - stride_in));
                    ptr_denoised_interm = ptr_denoised + kk + ((1 + jj) * stride_out - stride_out);
                    bottom_prev_cr      = _mm256_loadu_si256(
                        (__m256i *)(ptr_in_cr - 1 + kk + ((2 + jj) * stride_in_cr) - stride_in_cr));
                    bottom_next_cr = _mm256_loadu_si256(
                        (__m256i *)(ptr_in_cr + 1 + kk + ((2 + jj) * stride_in_cr) - stride_in_cr));
                    bottom_cr = _mm256_loadu_si256(
                        (__m256i *)((ptr_in_cr + kk) + (2 + jj) * stride_in_cr - stride_in_cr));
                    ptr_denoised_interm_cr =
                        ptr_denoised_cr + kk + ((1 + jj) * stride_out_cr - stride_out_cr);
                }

                chroma_weak_luma_strong_filter_avx2_intrin(top,
                                                           curr,
                                                           bottom,
                                                           curr_prev,
                                                           curr_next,
                                                           top_prev,
                                                           top_next,
                                                           bottom_prev,
                                                           bottom_next,
                                                           ptr_denoised_interm);

                chroma_weak_luma_strong_filter_avx2_intrin(top_cr,
                                                           curr_cr,
                                                           bottom_cr,
                                                           curr_prev_cr,
                                                           curr_next_cr,
                                                           top_prev_cr,
                                                           top_next_cr,
                                                           bottom_prev_cr,
                                                           bottom_next_cr,
                                                           ptr_denoised_interm_cr);

                top          = curr;
                curr         = bottom;
                top_prev     = curr_prev;
                top_next     = curr_next;
                curr_prev    = bottom_prev;
                curr_next    = bottom_next;
                top_cr       = curr_cr;
                curr_cr      = bottom_cr;
                top_prev_cr  = curr_prev_cr;
                top_next_cr  = curr_next_cr;
                curr_prev_cr = bottom_prev_cr;
                curr_next_cr = bottom_next_cr;
            }
        }

        for (; kk + 8 <= pic_width; kk += 8) {
            for (jj = 0; jj < sb_height; jj++) {
                if (sb_origin_y == 0) {
                    if (jj == 0) {
                        top_128  = _mm_loadl_epi64((__m128i *)(ptr_in + kk + jj * stride_in));
                        curr_128 = _mm_loadl_epi64((__m128i *)(ptr_in + kk + (1 + jj) * stride_in));
                        top_prev_128 =
                            _mm_loadl_epi64((__m128i *)(ptr_in - 1 + kk + ((jj)*stride_in)));
                        top_next_128 =
                            _mm_loadl_epi64((__m128i *)(ptr_in + 1 + kk + ((jj)*stride_in)));
                        curr_prev_128 =
                            _mm_loadl_epi64((__m128i *)(ptr_in - 1 + kk + ((1 + jj) * stride_in)));
                        curr_next_128 =
                            _mm_loadl_epi64((__m128i *)(ptr_in + 1 + kk + ((1 + jj) * stride_in)));
                        _mm_storel_epi64((__m128i *)(ptr_denoised + kk), top_128);
                        top_cr_128 =
                            _mm_loadl_epi64((__m128i *)(ptr_in_cr + kk + jj * stride_in_cr));
                        curr_cr_128 =
                            _mm_loadl_epi64((__m128i *)(ptr_in_cr + kk + (1 + jj) * stride_in_cr));
                        top_prev_cr_128 =
                            _mm_loadl_epi64((__m128i *)(ptr_in_cr - 1 + kk + ((jj)*stride_in_cr)));
                        top_next_cr_128 =
                            _mm_loadl_epi64((__m128i *)(ptr_in_cr + 1 + kk + ((jj)*stride_in_cr)));
                        curr_prev_cr_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in_cr - 1 + kk + ((1 + jj) * stride_in_cr)));
                        curr_next_cr_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in_cr + 1 + kk + ((1 + jj) * stride_in_cr)));
                        _mm_storel_epi64((__m128i *)(ptr_denoised_cr + kk), top_cr_128);
                    }
                    bottom_prev_128 =
                        _mm_loadl_epi64((__m128i *)(ptr_in - 1 + kk + ((2 + jj) * stride_in)));
                    bottom_next_128 =
                        _mm_loadl_epi64((__m128i *)(ptr_in + 1 + kk + ((2 + jj) * stride_in)));
                    bottom_128 = _mm_loadl_epi64((__m128i *)((ptr_in + kk) + (2 + jj) * stride_in));
                    ptr_denoised_interm = ptr_denoised + kk + ((1 + jj) * stride_out);
                    bottom_prev_cr_128  = _mm_loadl_epi64(
                        (__m128i *)(ptr_in_cr - 1 + kk + ((2 + jj) * stride_in_cr)));
                    bottom_next_cr_128 = _mm_loadl_epi64(
                        (__m128i *)(ptr_in_cr + 1 + kk + ((2 + jj) * stride_in_cr)));
                    bottom_cr_128 =
                        _mm_loadl_epi64((__m128i *)((ptr_in_cr + kk) + (2 + jj) * stride_in_cr));
                    ptr_denoised_interm_cr = ptr_denoised_cr + kk + ((1 + jj) * stride_out_cr);
                } else {
                    if (jj == 0) {
                        top_128 =
                            _mm_loadl_epi64((__m128i *)(ptr_in + kk + jj * stride_in - stride_in));
                        curr_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in + kk + (1 + jj) * stride_in - stride_in));
                        top_prev_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in - 1 + kk + ((jj)*stride_in) - stride_in));
                        top_next_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in + 1 + kk + ((jj)*stride_in) - stride_in));
                        curr_prev_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in - 1 + kk + ((1 + jj) * stride_in - stride_in)));
                        curr_next_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in + 1 + kk + ((1 + jj) * stride_in - stride_in)));
                        top_cr_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in_cr + kk + jj * stride_in_cr - stride_in_cr));
                        curr_cr_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in_cr + kk + (1 + jj) * stride_in_cr - stride_in_cr));
                        top_prev_cr_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in_cr - 1 + kk + ((jj)*stride_in_cr) - stride_in_cr));
                        top_next_cr_128 = _mm_loadl_epi64(
                            (__m128i *)(ptr_in_cr + 1 + kk + ((jj)*stride_in_cr) - stride_in_cr));
                        curr_prev_cr_128 =
                            _mm_loadl_epi64((__m128i *)(ptr_in_cr - 1 + kk +
                                                        ((1 + jj) * stride_in_cr - stride_in_cr)));
                        curr_next_cr_128 =
                            _mm_loadl_epi64((__m128i *)(ptr_in_cr + 1 + kk +
                                                        ((1 + jj) * stride_in_cr - stride_in_cr)));
                    }
                    bottom_prev_128 = _mm_loadl_epi64(
                        (__m128i *)(ptr_in - 1 + kk + ((2 + jj) * stride_in) - stride_in));
                    bottom_next_128 = _mm_loadl_epi64(
                        (__m128i *)(ptr_in + 1 + kk + ((2 + jj) * stride_in) - stride_in));
                    bottom_128 = _mm_loadl_epi64(
                        (__m128i *)((ptr_in + kk) + (2 + jj) * stride_in - stride_in));
                    ptr_denoised_interm = ptr_denoised + kk + ((1 + jj) * stride_out - stride_out);
                    bottom_prev_cr_128  = _mm_loadl_epi64(
                        (__m128i *)(ptr_in_cr - 1 + kk + ((2 + jj) * stride_in_cr) - stride_in_cr));
                    bottom_next_cr_128 = _mm_loadl_epi64(
                        (__m128i *)(ptr_in_cr + 1 + kk + ((2 + jj) * stride_in_cr) - stride_in_cr));
                    bottom_cr_128 = _mm_loadl_epi64(
                        (__m128i *)((ptr_in_cr + kk) + (2 + jj) * stride_in_cr - stride_in_cr));
                    ptr_denoised_interm_cr =
                        ptr_denoised_cr + kk + ((1 + jj) * stride_out_cr - stride_out_cr);
                }

                chroma_weak_luma_strong_filter_128_avx2_intrin(top_128,
                                                               curr_128,
                                                               bottom_128,
                                                               curr_prev_128,
                                                               curr_next_128,
                                                               top_prev_128,
                                                               top_next_128,
                                                               bottom_prev_128,
                                                               bottom_next_128,
                                                               ptr_denoised_interm);

                chroma_weak_luma_strong_filter_128_avx2_intrin(top_cr_128,
                                                               curr_cr_128,
                                                               bottom_cr_128,
                                                               curr_prev_cr_128,
                                                               curr_next_cr_128,
                                                               top_prev_cr_128,
                                                               top_next_cr_128,
                                                               bottom_prev_cr_128,
                                                               bottom_next_cr_128,
                                                               ptr_denoised_interm_cr);

                top_128          = curr_128;
                curr_128         = bottom_128;
                top_prev_128     = curr_prev_128;
                top_next_128     = curr_next_128;
                curr_prev_128    = bottom_prev_128;
                curr_next_128    = bottom_next_128;
                top_cr_128       = curr_cr_128;
                curr_cr_128      = bottom_cr_128;
                top_prev_cr_128  = curr_prev_cr_128;
                top_next_cr_128  = curr_next_cr_128;
                curr_prev_cr_128 = bottom_prev_cr_128;
                curr_next_cr_128 = bottom_next_cr_128;
            }
        }

        sb_height = MIN(BLOCK_SIZE_64 / 2, pic_height - sb_origin_y);
        for (jj = 0; jj < sb_height; jj++) {
            for (ii = idx; ii < pic_width; ii++) {
                if (!((jj < sb_height - 1 || (sb_origin_y + sb_height) < pic_height) && ii > 0 &&
                      ii < pic_width - 1)) {
                    ptr_denoised[ii + jj * stride_out]    = ptr_in[ii + jj * stride_in];
                    ptr_denoised_cr[ii + jj * stride_out] = ptr_in_cr[ii + jj * stride_in_cr];
                }
            }
        }
    }
}
