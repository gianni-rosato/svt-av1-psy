/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#define MACRO_VERTICAL_LUMA_8(A, B, C)\
    _mm_storel_epi64((__m128i*)prediction_ptr, _mm_or_si128(_mm_and_si128(A, B), C)); \
    A = _mm_srli_si128(A, 1); \
    _mm_storel_epi64((__m128i*)(prediction_ptr + pStride), _mm_or_si128(_mm_and_si128(A, B), C)); \
    A = _mm_srli_si128(A, 1); \
    _mm_storel_epi64((__m128i*)(prediction_ptr + 2*pStride), _mm_or_si128(_mm_and_si128(A, B), C)); \
    A = _mm_srli_si128(A, 1); \
    _mm_storel_epi64((__m128i*)(prediction_ptr + 3*pStride), _mm_or_si128(_mm_and_si128(A, B), C)); \
    A = _mm_srli_si128(A, 1);

#define MACRO_VERTICAL_LUMA_16(A, B, C)\
    _mm_storeu_si128((__m128i*)prediction_ptr, _mm_or_si128(_mm_and_si128(A, B), C)); \
    A = _mm_srli_si128(A, 1); \
    _mm_storeu_si128((__m128i*)(prediction_ptr + pStride), _mm_or_si128(_mm_and_si128(A, B), C)); \
    A = _mm_srli_si128(A, 1); \
    _mm_storeu_si128((__m128i*)(prediction_ptr + 2*pStride), _mm_or_si128(_mm_and_si128(A, B), C)); \
    A = _mm_srli_si128(A, 1); \
    _mm_storeu_si128((__m128i*)(prediction_ptr + 3*pStride), _mm_or_si128(_mm_and_si128(A, B), C)); \
    A = _mm_srli_si128(A, 1);
