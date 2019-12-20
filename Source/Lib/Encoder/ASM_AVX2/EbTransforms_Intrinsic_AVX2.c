/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbTransforms_AVX2.h"

#include <emmintrin.h>
#include <immintrin.h>

// Coefficients for forward 16/32-point transform
#ifdef __GNUC__
__attribute__((aligned(16)))
#endif
EB_ALIGN(32) const int16_t coeff_tbl_AVX2[48 * 16] =
{
    64, 64, 89, 75, 83, 36, 75, -18, 64, 64, 89, 75, 83, 36, 75, -18, 64, -64, 50, -89, 36, -83, 18, -50, 64, -64, 50, -89, 36, -83, 18, -50,
    64, 64, 50, 18, -36, -83, -89, -50, 64, 64, 50, 18, -36, -83, -89, -50, -64, 64, 18, 75, 83, -36, 75, -89, -64, 64, 18, 75, 83, -36, 75, -89,
    64, 64, -18, -50, -83, -36, 50, 89, 64, 64, -18, -50, -83, -36, 50, 89, 64, -64, -75, -18, -36, 83, 89, -75, 64, -64, -75, -18, -36, 83, 89, -75,
    64, 64, -75, -89, 36, 83, 18, -75, 64, 64, -75, -89, 36, 83, 18, -75, -64, 64, 89, -50, -83, 36, 50, -18, -64, 64, 89, -50, -83, 36, 50, -18,
    90, 87, 87, 57, 80, 9, 70, -43, 90, 87, 87, 57, 80, 9, 70, -43, 57, -80, 43, -90, 25, -70, 9, -25, 57, -80, 43, -90, 25, -70, 9, -25,
    80, 70, 9, -43, -70, -87, -87, 9, 80, 70, 9, -43, -70, -87, -87, 9, -25, 90, 57, 25, 90, -80, 43, -57, -25, 90, 57, 25, 90, -80, 43, -57,
    57, 43, -80, -90, -25, 57, 90, 25, 57, 43, -80, -90, -25, 57, 90, 25, -9, -87, -87, 70, 43, 9, 70, -80, -9, -87, -87, 70, 43, 9, 70, -80,
    25, 9, -70, -25, 90, 43, -80, -57, 25, 9, -70, -25, 90, 43, -80, -57, 43, 70, 9, -80, -57, 87, 87, -90, 43, 70, 9, -80, -57, 87, 87, -90,
    90, 90, 90, 82, 88, 67, 85, 46, 90, 90, 90, 82, 88, 67, 85, 46, 82, 22, 78, -4, 73, -31, 67, -54, 82, 22, 78, -4, 73, -31, 67, -54,
    61, -73, 54, -85, 46, -90, 38, -88, 61, -73, 54, -85, 46, -90, 38, -88, 31, -78, 22, -61, 13, -38, 4, -13, 31, -78, 22, -61, 13, -38, 4, -13,
    88, 85, 67, 46, 31, -13, -13, -67, 88, 85, 67, 46, 31, -13, -13, -67, -54, -90, -82, -73, -90, -22, -78, 38, -54, -90, -82, -73, -90, -22, -78, 38,
    -46, 82, -4, 88, 38, 54, 73, -4, -46, 82, -4, 88, 38, 54, 73, -4, 90, -61, 85, -90, 61, -78, 22, -31, 90, -61, 85, -90, 61, -78, 22, -31,
    82, 78, 22, -4, -54, -82, -90, -73, 82, 78, 22, -4, -54, -82, -90, -73, -61, 13, 13, 85, 78, 67, 85, -22, -61, 13, 13, 85, 78, 67, 85, -22,
    31, -88, -46, -61, -90, 31, -67, 90, 31, -88, -46, -61, -90, 31, -67, 90, 4, 54, 73, -38, 88, -90, 38, -46, 4, 54, 73, -38, 88, -90, 38, -46,
    73, 67, -31, -54, -90, -78, -22, 38, 73, 67, -31, -54, -90, -78, -22, 38, 78, 85, 67, -22, -38, -90, -90, 4, 78, 85, 67, -22, -38, -90, -90, 4,
    -13, 90, 82, 13, 61, -88, -46, -31, -13, 90, 82, 13, 61, -88, -46, -31, -88, 82, -4, 46, 85, -73, 54, -61, -88, 82, -4, 46, 85, -73, 54, -61,
    61, 54, -73, -85, -46, -4, 82, 88, 61, 54, -73, -85, -46, -4, 82, 88, 31, -46, -88, -61, -13, 82, 90, 13, 31, -46, -88, -61, -13, 82, 90, 13,
    -4, -90, -90, 38, 22, 67, 85, -78, -4, -90, -90, 38, 22, 67, 85, -78, -38, -22, -78, 90, 54, -31, 67, -73, -38, -22, -78, 90, 54, -31, 67, -73,
    46, 38, -90, -88, 38, 73, 54, -4, 46, 38, -90, -88, 38, 73, 54, -4, -90, -67, 31, 90, 61, -46, -88, -31, -90, -67, 31, 90, 61, -46, -88, -31,
    22, 85, 67, -78, -85, 13, 13, 61, 22, 85, 67, -78, -85, 13, 13, 61, 73, -90, -82, 54, 4, 22, 78, -82, 73, -90, -82, 54, 4, 22, 78, -82,
    31, 22, -78, -61, 90, 85, -61, -90, 31, 22, -78, -61, 90, 85, -61, -90, 4, 73, 54, -38, -88, -4, 82, 46, 4, 73, 54, -38, -88, -4, 82, 46,
    -38, -78, -22, 90, 73, -82, -90, 54, -38, -78, -22, 90, 73, -82, -90, 54, 67, -13, -13, -31, -46, 67, 85, -88, 67, -13, -13, -31, -46, 67, 85, -88,
    13, 4, -38, -13, 61, 22, -78, -31, 13, 4, -38, -13, 61, 22, -78, -31, 88, 38, -90, -46, 85, 54, -73, -61, 88, 38, -90, -46, 85, 54, -73, -61,
    54, 67, -31, -73, 4, 78, 22, -82, 54, 67, -31, -73, 4, 78, 22, -82, -46, 85, 67, -88, -82, 90, 90, -90, -46, 85, 67, -88, -82, 90, 90, -90
};

/*******************************************************************************
* Requirement: area_size = 4, 8, or area_size % 16 = 0
*******************************************************************************/

// transpose 16x16 block of data
void transpose16_AVX2_INTRIN(int16_t *src, uint32_t src_stride, int16_t *dst, uint32_t dst_stride)
{
    uint32_t i;
    for (i = 0; i < 2; i++)
    {
        __m256i a0, a1, a2, a3, a4, a5, a6, a7;
        __m256i b0, b1, b2, b3, b4, b5, b6, b7;

        a0 = _mm256_loadu_si256((const __m256i *)(src + (8 * i + 0)*src_stride));
        a1 = _mm256_loadu_si256((const __m256i *)(src + (8 * i + 1)*src_stride));
        a2 = _mm256_loadu_si256((const __m256i *)(src + (8 * i + 2)*src_stride));
        a3 = _mm256_loadu_si256((const __m256i *)(src + (8 * i + 3)*src_stride));
        a4 = _mm256_loadu_si256((const __m256i *)(src + (8 * i + 4)*src_stride));
        a5 = _mm256_loadu_si256((const __m256i *)(src + (8 * i + 5)*src_stride));
        a6 = _mm256_loadu_si256((const __m256i *)(src + (8 * i + 6)*src_stride));
        a7 = _mm256_loadu_si256((const __m256i *)(src + (8 * i + 7)*src_stride));

        b0 = _mm256_unpacklo_epi16(a0, a4);
        b1 = _mm256_unpacklo_epi16(a1, a5);
        b2 = _mm256_unpacklo_epi16(a2, a6);
        b3 = _mm256_unpacklo_epi16(a3, a7);
        b4 = _mm256_unpackhi_epi16(a0, a4);
        b5 = _mm256_unpackhi_epi16(a1, a5);
        b6 = _mm256_unpackhi_epi16(a2, a6);
        b7 = _mm256_unpackhi_epi16(a3, a7);

        a0 = _mm256_unpacklo_epi16(b0, b2);
        a1 = _mm256_unpacklo_epi16(b1, b3);
        a2 = _mm256_unpackhi_epi16(b0, b2);
        a3 = _mm256_unpackhi_epi16(b1, b3);
        a4 = _mm256_unpacklo_epi16(b4, b6);
        a5 = _mm256_unpacklo_epi16(b5, b7);
        a6 = _mm256_unpackhi_epi16(b4, b6);
        a7 = _mm256_unpackhi_epi16(b5, b7);

        b0 = _mm256_unpacklo_epi16(a0, a1);
        b1 = _mm256_unpackhi_epi16(a0, a1);
        b2 = _mm256_unpacklo_epi16(a2, a3);
        b3 = _mm256_unpackhi_epi16(a2, a3);
        b4 = _mm256_unpacklo_epi16(a4, a5);
        b5 = _mm256_unpackhi_epi16(a4, a5);
        b6 = _mm256_unpacklo_epi16(a6, a7);
        b7 = _mm256_unpackhi_epi16(a6, a7);

        _mm_storeu_si128((__m128i *)(dst + 0 * dst_stride + 8 * i), _mm256_extracti128_si256(b0, 0));
        _mm_storeu_si128((__m128i *)(dst + 1 * dst_stride + 8 * i), _mm256_extracti128_si256(b1, 0));
        _mm_storeu_si128((__m128i *)(dst + 2 * dst_stride + 8 * i), _mm256_extracti128_si256(b2, 0));
        _mm_storeu_si128((__m128i *)(dst + 3 * dst_stride + 8 * i), _mm256_extracti128_si256(b3, 0));
        _mm_storeu_si128((__m128i *)(dst + 4 * dst_stride + 8 * i), _mm256_extracti128_si256(b4, 0));
        _mm_storeu_si128((__m128i *)(dst + 5 * dst_stride + 8 * i), _mm256_extracti128_si256(b5, 0));
        _mm_storeu_si128((__m128i *)(dst + 6 * dst_stride + 8 * i), _mm256_extracti128_si256(b6, 0));
        _mm_storeu_si128((__m128i *)(dst + 7 * dst_stride + 8 * i), _mm256_extracti128_si256(b7, 0));
        _mm_storeu_si128((__m128i *)(dst + 8 * dst_stride + 8 * i), _mm256_extracti128_si256(b0, 1));
        _mm_storeu_si128((__m128i *)(dst + 9 * dst_stride + 8 * i), _mm256_extracti128_si256(b1, 1));
        _mm_storeu_si128((__m128i *)(dst + 10 * dst_stride + 8 * i), _mm256_extracti128_si256(b2, 1));
        _mm_storeu_si128((__m128i *)(dst + 11 * dst_stride + 8 * i), _mm256_extracti128_si256(b3, 1));
        _mm_storeu_si128((__m128i *)(dst + 12 * dst_stride + 8 * i), _mm256_extracti128_si256(b4, 1));
        _mm_storeu_si128((__m128i *)(dst + 13 * dst_stride + 8 * i), _mm256_extracti128_si256(b5, 1));
        _mm_storeu_si128((__m128i *)(dst + 14 * dst_stride + 8 * i), _mm256_extracti128_si256(b6, 1));
        _mm_storeu_si128((__m128i *)(dst + 15 * dst_stride + 8 * i), _mm256_extracti128_si256(b7, 1));
    }
}

// transpose 32x32 block of data
void transpose32_AVX2_INTRIN(int16_t *src, uint32_t src_stride, int16_t *dst, uint32_t dst_stride)
{
    uint32_t i, j;
    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 2; j++)
        {
            __m256i a0, a1, a2, a3, a4, a5, a6, a7;
            __m256i b0, b1, b2, b3, b4, b5, b6, b7;

            a0 = _mm256_loadu_si256((const __m256i *)(src + (8 * i + 0)*src_stride + 16 * j));
            a1 = _mm256_loadu_si256((const __m256i *)(src + (8 * i + 1)*src_stride + 16 * j));
            a2 = _mm256_loadu_si256((const __m256i *)(src + (8 * i + 2)*src_stride + 16 * j));
            a3 = _mm256_loadu_si256((const __m256i *)(src + (8 * i + 3)*src_stride + 16 * j));
            a4 = _mm256_loadu_si256((const __m256i *)(src + (8 * i + 4)*src_stride + 16 * j));
            a5 = _mm256_loadu_si256((const __m256i *)(src + (8 * i + 5)*src_stride + 16 * j));
            a6 = _mm256_loadu_si256((const __m256i *)(src + (8 * i + 6)*src_stride + 16 * j));
            a7 = _mm256_loadu_si256((const __m256i *)(src + (8 * i + 7)*src_stride + 16 * j));

            b0 = _mm256_unpacklo_epi16(a0, a4);
            b1 = _mm256_unpacklo_epi16(a1, a5);
            b2 = _mm256_unpacklo_epi16(a2, a6);
            b3 = _mm256_unpacklo_epi16(a3, a7);
            b4 = _mm256_unpackhi_epi16(a0, a4);
            b5 = _mm256_unpackhi_epi16(a1, a5);
            b6 = _mm256_unpackhi_epi16(a2, a6);
            b7 = _mm256_unpackhi_epi16(a3, a7);

            a0 = _mm256_unpacklo_epi16(b0, b2);
            a1 = _mm256_unpacklo_epi16(b1, b3);
            a2 = _mm256_unpackhi_epi16(b0, b2);
            a3 = _mm256_unpackhi_epi16(b1, b3);
            a4 = _mm256_unpacklo_epi16(b4, b6);
            a5 = _mm256_unpacklo_epi16(b5, b7);
            a6 = _mm256_unpackhi_epi16(b4, b6);
            a7 = _mm256_unpackhi_epi16(b5, b7);

            b0 = _mm256_unpacklo_epi16(a0, a1);
            b1 = _mm256_unpackhi_epi16(a0, a1);
            b2 = _mm256_unpacklo_epi16(a2, a3);
            b3 = _mm256_unpackhi_epi16(a2, a3);
            b4 = _mm256_unpacklo_epi16(a4, a5);
            b5 = _mm256_unpackhi_epi16(a4, a5);
            b6 = _mm256_unpacklo_epi16(a6, a7);
            b7 = _mm256_unpackhi_epi16(a6, a7);

            _mm_storeu_si128((__m128i *)(dst + (16 * j + 0)*dst_stride + 8 * i), _mm256_extracti128_si256(b0, 0));
            _mm_storeu_si128((__m128i *)(dst + (16 * j + 1)*dst_stride + 8 * i), _mm256_extracti128_si256(b1, 0));
            _mm_storeu_si128((__m128i *)(dst + (16 * j + 2)*dst_stride + 8 * i), _mm256_extracti128_si256(b2, 0));
            _mm_storeu_si128((__m128i *)(dst + (16 * j + 3)*dst_stride + 8 * i), _mm256_extracti128_si256(b3, 0));
            _mm_storeu_si128((__m128i *)(dst + (16 * j + 4)*dst_stride + 8 * i), _mm256_extracti128_si256(b4, 0));
            _mm_storeu_si128((__m128i *)(dst + (16 * j + 5)*dst_stride + 8 * i), _mm256_extracti128_si256(b5, 0));
            _mm_storeu_si128((__m128i *)(dst + (16 * j + 6)*dst_stride + 8 * i), _mm256_extracti128_si256(b6, 0));
            _mm_storeu_si128((__m128i *)(dst + (16 * j + 7)*dst_stride + 8 * i), _mm256_extracti128_si256(b7, 0));
            _mm_storeu_si128((__m128i *)(dst + (16 * j + 8)*dst_stride + 8 * i), _mm256_extracti128_si256(b0, 1));
            _mm_storeu_si128((__m128i *)(dst + (16 * j + 9)*dst_stride + 8 * i), _mm256_extracti128_si256(b1, 1));
            _mm_storeu_si128((__m128i *)(dst + (16 * j + 10)*dst_stride + 8 * i), _mm256_extracti128_si256(b2, 1));
            _mm_storeu_si128((__m128i *)(dst + (16 * j + 11)*dst_stride + 8 * i), _mm256_extracti128_si256(b3, 1));
            _mm_storeu_si128((__m128i *)(dst + (16 * j + 12)*dst_stride + 8 * i), _mm256_extracti128_si256(b4, 1));
            _mm_storeu_si128((__m128i *)(dst + (16 * j + 13)*dst_stride + 8 * i), _mm256_extracti128_si256(b5, 1));
            _mm_storeu_si128((__m128i *)(dst + (16 * j + 14)*dst_stride + 8 * i), _mm256_extracti128_si256(b6, 1));
            _mm_storeu_si128((__m128i *)(dst + (16 * j + 15)*dst_stride + 8 * i), _mm256_extracti128_si256(b7, 1));
        }
    }
}

// 32-point forward transform (32 rows)
void transform32_AVX2_INTRIN(int16_t *src, uint32_t src_stride, int16_t *dst, uint32_t dst_stride, uint32_t shift)
{
    uint32_t i;
    __m128i s0;
    __m256i o0;
    const __m256i *coeff32 = (const __m256i *)coeff_tbl_AVX2;

    shift &= 0x0000FFFF; // Redundant code to fix Visual Studio 2012 AVX2 compiler error
    s0 = _mm_cvtsi32_si128(shift);
    o0 = _mm256_set1_epi32(1 << (shift - 1));

    for (i = 0; i < 16; i++)
    {
        __m256i x0, x1, x2, x3;
        __m256i y0, y1, y2, y3;
        __m256i a0, a1, a2, a3, a4, a5, a6, a7;
        __m256i b0, b1, b2, b3, b4, b5, b6, b7;

        x0 = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_loadu_si128((const __m128i *)(src + 0x00))), _mm_loadu_si128((const __m128i *)(src + src_stride + 0x00)), 0x1);
        x1 = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_loadu_si128((const __m128i *)(src + 0x08))), _mm_loadu_si128((const __m128i *)(src + src_stride + 0x08)), 0x1);
        x2 = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_loadu_si128((const __m128i *)(src + 0x10))), _mm_loadu_si128((const __m128i *)(src + src_stride + 0x10)), 0x1);
        x3 = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_loadu_si128((const __m128i *)(src + 0x18))), _mm_loadu_si128((const __m128i *)(src + src_stride + 0x18)), 0x1);

        // 32-point butterfly
        x2 = _mm256_shuffle_epi8(x2, _mm256_setr_epi8(14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1, 14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1));
        x3 = _mm256_shuffle_epi8(x3, _mm256_setr_epi8(14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1, 14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1));

        y0 = _mm256_add_epi16(x0, x3);
        y1 = _mm256_add_epi16(x1, x2);

        y2 = _mm256_sub_epi16(x0, x3);
        y3 = _mm256_sub_epi16(x1, x2);

        // 16-point butterfly
        y1 = _mm256_shuffle_epi8(y1, _mm256_setr_epi8(14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1, 14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1));

        x0 = _mm256_add_epi16(y0, y1);
        x1 = _mm256_sub_epi16(y0, y1);

        x2 = y2;
        x3 = y3;

        a0 = _mm256_madd_epi16(_mm256_shuffle_epi32(x0, 0x00), coeff32[0]);
        a0 = _mm256_add_epi32(a0, _mm256_madd_epi16(_mm256_shuffle_epi32(x0, 0x55), coeff32[2]));
        a0 = _mm256_add_epi32(a0, _mm256_madd_epi16(_mm256_shuffle_epi32(x0, 0xaa), coeff32[4]));
        a0 = _mm256_add_epi32(a0, _mm256_madd_epi16(_mm256_shuffle_epi32(x0, 0xff), coeff32[6]));

        a1 = _mm256_madd_epi16(_mm256_shuffle_epi32(x0, 0x00), coeff32[1]);
        a1 = _mm256_add_epi32(a1, _mm256_madd_epi16(_mm256_shuffle_epi32(x0, 0x55), coeff32[3]));
        a1 = _mm256_add_epi32(a1, _mm256_madd_epi16(_mm256_shuffle_epi32(x0, 0xaa), coeff32[5]));
        a1 = _mm256_add_epi32(a1, _mm256_madd_epi16(_mm256_shuffle_epi32(x0, 0xff), coeff32[7]));

        a2 = _mm256_madd_epi16(_mm256_shuffle_epi32(x1, 0x00), coeff32[8]);
        a2 = _mm256_add_epi32(a2, _mm256_madd_epi16(_mm256_shuffle_epi32(x1, 0x55), coeff32[10]));
        a2 = _mm256_add_epi32(a2, _mm256_madd_epi16(_mm256_shuffle_epi32(x1, 0xaa), coeff32[12]));
        a2 = _mm256_add_epi32(a2, _mm256_madd_epi16(_mm256_shuffle_epi32(x1, 0xff), coeff32[14]));

        a3 = _mm256_madd_epi16(_mm256_shuffle_epi32(x1, 0x00), coeff32[9]);
        a3 = _mm256_add_epi32(a3, _mm256_madd_epi16(_mm256_shuffle_epi32(x1, 0x55), coeff32[11]));
        a3 = _mm256_add_epi32(a3, _mm256_madd_epi16(_mm256_shuffle_epi32(x1, 0xaa), coeff32[13]));
        a3 = _mm256_add_epi32(a3, _mm256_madd_epi16(_mm256_shuffle_epi32(x1, 0xff), coeff32[15]));

        a4 = _mm256_madd_epi16(_mm256_shuffle_epi32(x2, 0x00), coeff32[16]);
        a4 = _mm256_add_epi32(a4, _mm256_madd_epi16(_mm256_shuffle_epi32(x2, 0x55), coeff32[20]));
        a4 = _mm256_add_epi32(a4, _mm256_madd_epi16(_mm256_shuffle_epi32(x2, 0xaa), coeff32[24]));
        a4 = _mm256_add_epi32(a4, _mm256_madd_epi16(_mm256_shuffle_epi32(x2, 0xff), coeff32[28]));
        a4 = _mm256_add_epi32(a4, _mm256_madd_epi16(_mm256_shuffle_epi32(x3, 0x00), coeff32[32]));
        a4 = _mm256_add_epi32(a4, _mm256_madd_epi16(_mm256_shuffle_epi32(x3, 0x55), coeff32[36]));
        a4 = _mm256_add_epi32(a4, _mm256_madd_epi16(_mm256_shuffle_epi32(x3, 0xaa), coeff32[40]));
        a4 = _mm256_add_epi32(a4, _mm256_madd_epi16(_mm256_shuffle_epi32(x3, 0xff), coeff32[44]));

        a5 = _mm256_madd_epi16(_mm256_shuffle_epi32(x2, 0x00), coeff32[17]);
        a5 = _mm256_add_epi32(a5, _mm256_madd_epi16(_mm256_shuffle_epi32(x2, 0x55), coeff32[21]));
        a5 = _mm256_add_epi32(a5, _mm256_madd_epi16(_mm256_shuffle_epi32(x2, 0xaa), coeff32[25]));
        a5 = _mm256_add_epi32(a5, _mm256_madd_epi16(_mm256_shuffle_epi32(x2, 0xff), coeff32[29]));
        a5 = _mm256_add_epi32(a5, _mm256_madd_epi16(_mm256_shuffle_epi32(x3, 0x00), coeff32[33]));
        a5 = _mm256_add_epi32(a5, _mm256_madd_epi16(_mm256_shuffle_epi32(x3, 0x55), coeff32[37]));
        a5 = _mm256_add_epi32(a5, _mm256_madd_epi16(_mm256_shuffle_epi32(x3, 0xaa), coeff32[41]));
        a5 = _mm256_add_epi32(a5, _mm256_madd_epi16(_mm256_shuffle_epi32(x3, 0xff), coeff32[45]));

        a6 = _mm256_madd_epi16(_mm256_shuffle_epi32(x2, 0x00), coeff32[18]);
        a6 = _mm256_add_epi32(a6, _mm256_madd_epi16(_mm256_shuffle_epi32(x2, 0x55), coeff32[22]));
        a6 = _mm256_add_epi32(a6, _mm256_madd_epi16(_mm256_shuffle_epi32(x2, 0xaa), coeff32[26]));
        a6 = _mm256_add_epi32(a6, _mm256_madd_epi16(_mm256_shuffle_epi32(x2, 0xff), coeff32[30]));
        a6 = _mm256_add_epi32(a6, _mm256_madd_epi16(_mm256_shuffle_epi32(x3, 0x00), coeff32[34]));
        a6 = _mm256_add_epi32(a6, _mm256_madd_epi16(_mm256_shuffle_epi32(x3, 0x55), coeff32[38]));
        a6 = _mm256_add_epi32(a6, _mm256_madd_epi16(_mm256_shuffle_epi32(x3, 0xaa), coeff32[42]));
        a6 = _mm256_add_epi32(a6, _mm256_madd_epi16(_mm256_shuffle_epi32(x3, 0xff), coeff32[46]));

        a7 = _mm256_madd_epi16(_mm256_shuffle_epi32(x2, 0x00), coeff32[19]);
        a7 = _mm256_add_epi32(a7, _mm256_madd_epi16(_mm256_shuffle_epi32(x2, 0x55), coeff32[23]));
        a7 = _mm256_add_epi32(a7, _mm256_madd_epi16(_mm256_shuffle_epi32(x2, 0xaa), coeff32[27]));
        a7 = _mm256_add_epi32(a7, _mm256_madd_epi16(_mm256_shuffle_epi32(x2, 0xff), coeff32[31]));
        a7 = _mm256_add_epi32(a7, _mm256_madd_epi16(_mm256_shuffle_epi32(x3, 0x00), coeff32[35]));
        a7 = _mm256_add_epi32(a7, _mm256_madd_epi16(_mm256_shuffle_epi32(x3, 0x55), coeff32[39]));
        a7 = _mm256_add_epi32(a7, _mm256_madd_epi16(_mm256_shuffle_epi32(x3, 0xaa), coeff32[43]));
        a7 = _mm256_add_epi32(a7, _mm256_madd_epi16(_mm256_shuffle_epi32(x3, 0xff), coeff32[47]));

        b0 = _mm256_sra_epi32(_mm256_add_epi32(a0, o0), s0);
        b1 = _mm256_sra_epi32(_mm256_add_epi32(a1, o0), s0);
        b2 = _mm256_sra_epi32(_mm256_add_epi32(a2, o0), s0);
        b3 = _mm256_sra_epi32(_mm256_add_epi32(a3, o0), s0);
        b4 = _mm256_sra_epi32(_mm256_add_epi32(a4, o0), s0);
        b5 = _mm256_sra_epi32(_mm256_add_epi32(a5, o0), s0);
        b6 = _mm256_sra_epi32(_mm256_add_epi32(a6, o0), s0);
        b7 = _mm256_sra_epi32(_mm256_add_epi32(a7, o0), s0);

        x0 = _mm256_packs_epi32(b0, b1);
        x1 = _mm256_packs_epi32(b2, b3);
        x2 = _mm256_packs_epi32(b4, b5);
        x3 = _mm256_packs_epi32(b6, b7);

        y0 = _mm256_unpacklo_epi16(x0, x1);
        y1 = _mm256_unpackhi_epi16(x0, x1);
        y2 = x2;
        y3 = x3;
        x0 = _mm256_unpacklo_epi16(y0, y2);
        x1 = _mm256_unpackhi_epi16(y0, y2);
        x2 = _mm256_unpacklo_epi16(y1, y3);
        x3 = _mm256_unpackhi_epi16(y1, y3);

        y0 = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm256_extracti128_si256(x0, 0)), _mm256_extracti128_si256(x1, 0), 0x1);
        y1 = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm256_extracti128_si256(x2, 0)), _mm256_extracti128_si256(x3, 0), 0x1);
        y2 = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm256_extracti128_si256(x0, 1)), _mm256_extracti128_si256(x1, 1), 0x1);
        y3 = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm256_extracti128_si256(x2, 1)), _mm256_extracti128_si256(x3, 1), 0x1);
        _mm256_storeu_si256((__m256i *)(dst + 0x00), y0);
        _mm256_storeu_si256((__m256i *)(dst + 0x10), y1);
        _mm256_storeu_si256((__m256i *)(dst + dst_stride + 0x00), y2);
        _mm256_storeu_si256((__m256i *)(dst + dst_stride + 0x10), y3);

        src += 2 * src_stride;
        dst += 2 * dst_stride;
    }
}

void mat_mult4x4_out_buff_avx2_intrin(
    int16_t*              coeff,
    const uint32_t         coeff_stride,
    int16_t*              coeff_out,
    const uint32_t         coeff_out_stride,
    const uint16_t        *maskingMatrix,
    const uint32_t         maskingMatrixStride,
    const uint32_t         compute_size,
    const int32_t         offset,
    const int32_t         shift_num,
    uint32_t*              nonzerocoeff)

{
    __m256i z = _mm256_setzero_si256();
    __m128i a, b;
    __m256i coeffTemp, a0, a1, b0, b1, ymm_computed, MaskingMatrix, offsetREG, coeffTempORG;
    (void)compute_size;

    coeffTemp = a0 = a1 = b0 = b1 = ymm_computed = MaskingMatrix = offsetREG = _mm256_setzero_si256();

    // prepare Shift REG
    __m128i PMP_PRECISION_REG = _mm_set_epi16(0, 0, 0, 0, 0, 0, 0, (int16_t)shift_num); //_mm_set1_epi16((uint16_t)shift_num);//_mm_set1_epi32(shift_num);

    //prepare the offset
    offsetREG = _mm256_set1_epi32(offset);

    //load maskingMatrix_new
    MaskingMatrix = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_unpacklo_epi64(_mm_loadl_epi64((__m128i*)maskingMatrix), _mm_loadl_epi64((__m128i*)(maskingMatrix + maskingMatrixStride)))), _mm_unpacklo_epi64(_mm_loadl_epi64((__m128i*)(maskingMatrix + 2 * maskingMatrixStride)), _mm_loadl_epi64((__m128i*)(maskingMatrix + 3 * maskingMatrixStride))), 0x1);

    //load coefftemp
    a = _mm_unpacklo_epi64(_mm_loadl_epi64((__m128i*)coeff), _mm_loadl_epi64((__m128i*)(coeff + coeff_stride))); // 1st and 2nd row of the 4x4 block
    b = _mm_unpacklo_epi64(_mm_loadl_epi64((__m128i*)(coeff + 2 * coeff_stride)), _mm_loadl_epi64((__m128i*)(coeff + 3 * coeff_stride))); // 3rd and 4th row of the 4x4 block
    coeffTemp = _mm256_insertf128_si256(_mm256_castsi128_si256(a), b, 0x1); // the 4x4 block is now loaded

    coeffTempORG = coeffTemp;
    //Absolute val
    coeffTemp = _mm256_abs_epi16(coeffTemp);

    a0 = _mm256_mullo_epi16(coeffTemp, MaskingMatrix);
    a1 = _mm256_mulhi_epi16(coeffTemp, MaskingMatrix);

    b0 = _mm256_unpacklo_epi16(a0, a1);
    b1 = _mm256_unpackhi_epi16(a0, a1);

    b0 = _mm256_add_epi32(b0, offsetREG);
    b1 = _mm256_add_epi32(b1, offsetREG);

    //Shift right by PMP_PRECISION_REG
    b0 = _mm256_sra_epi32(b0, PMP_PRECISION_REG);
    b1 = _mm256_sra_epi32(b1, PMP_PRECISION_REG);

    //coefftemp in c
    ymm_computed = _mm256_packs_epi32(b0, b1);//Convert packed 32-bit integers from a and b to packed 16-bit integers using signed saturation, and store the results in dst.
    z = _mm256_sub_epi16(z, _mm256_cmpgt_epi16(ymm_computed, _mm256_setzero_si256())); //coeffTemp = (coeff[coeffLocation] < 0)? -coeffTemp : coeffTemp;

    ymm_computed = _mm256_sign_epi16(ymm_computed, coeffTempORG);// coeffTemp);

    a = _mm256_extracti128_si256(ymm_computed, 0);
    b = _mm256_extracti128_si256(ymm_computed, 1);
    _mm_storel_epi64((__m128i *) coeff_out, a);
    _mm_storel_epi64((__m128i *)(coeff_out + 1 * coeff_out_stride), _mm_srli_si128(a, 8));
    _mm_storel_epi64((__m128i *)(coeff_out + 2 * coeff_out_stride), b);
    _mm_storel_epi64((__m128i *)(coeff_out + 3 * coeff_out_stride), _mm_srli_si128(b, 8));

    z = _mm256_sad_epu8(z, _mm256_srli_si256(z, 8));
    *nonzerocoeff = _mm_cvtsi128_si32(_mm_add_epi32(_mm256_extracti128_si256(z, 0), _mm256_extracti128_si256(z, 1)));
}

/*****************************************************************************************************************************************************************/
void mat_mult4x4_avx2_intrin(
    int16_t*              coeff,
    const uint32_t         coeff_stride,
    const uint16_t        *maskingMatrix,
    const uint32_t         maskingMatrixStride,  //Matrix size
    const uint32_t         compute_size,  //Computation area size
    const int32_t         offset,     //(PMP_MAX >> 1)
    const int32_t         shift_num, //PMP_PRECISION
    uint32_t*              nonzerocoeff)

{
    __m256i z = _mm256_setzero_si256();
    __m128i a, b;
    __m256i coeffTemp, a0, a1, b0, b1, ymm_computed, MaskingMatrix, offsetREG, coeffTempORG;
    (void)compute_size;

    coeffTemp = a0 = a1 = b0 = b1 = ymm_computed = MaskingMatrix = offsetREG = _mm256_setzero_si256();

    // prepare Shift REG
    __m128i PMP_PRECISION_REG = _mm_set_epi16(0, 0, 0, 0, 0, 0, 0, (int16_t)shift_num); //_mm_set1_epi16((uint16_t)shift_num);//_mm_set1_epi32(shift_num);

    //prepare the offset
    offsetREG = _mm256_set1_epi32(offset);

    //load maskingMatrix_new
    MaskingMatrix = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_unpacklo_epi64(_mm_loadl_epi64((__m128i*)maskingMatrix), _mm_loadl_epi64((__m128i*)(maskingMatrix + maskingMatrixStride)))), _mm_unpacklo_epi64(_mm_loadl_epi64((__m128i*)(maskingMatrix + 2 * maskingMatrixStride)), _mm_loadl_epi64((__m128i*)(maskingMatrix + 3 * maskingMatrixStride))), 0x1);

    //load coefftemp
    a = _mm_unpacklo_epi64(_mm_loadl_epi64((__m128i*)coeff), _mm_loadl_epi64((__m128i*)(coeff + coeff_stride))); // 1st and 2nd row of the 4x4 block
    b = _mm_unpacklo_epi64(_mm_loadl_epi64((__m128i*)(coeff + 2 * coeff_stride)), _mm_loadl_epi64((__m128i*)(coeff + 3 * coeff_stride))); // 3rd and 4th row of the 4x4 block
    coeffTemp = _mm256_insertf128_si256(_mm256_castsi128_si256(a), b, 0x1); // the 4x4 block is now loaded

    coeffTempORG = coeffTemp;
    //Absolute val
    coeffTemp = _mm256_abs_epi16(coeffTemp);

    a0 = _mm256_mullo_epi16(coeffTemp, MaskingMatrix);
    a1 = _mm256_mulhi_epi16(coeffTemp, MaskingMatrix);

    b0 = _mm256_unpacklo_epi16(a0, a1);
    b1 = _mm256_unpackhi_epi16(a0, a1);

    b0 = _mm256_add_epi32(b0, offsetREG);
    b1 = _mm256_add_epi32(b1, offsetREG);

    //Shift right by PMP_PRECISION_REG
    b0 = _mm256_sra_epi32(b0, PMP_PRECISION_REG);
    b1 = _mm256_sra_epi32(b1, PMP_PRECISION_REG);

    //coefftemp in c
    ymm_computed = _mm256_packs_epi32(b0, b1);//Convert packed 32-bit integers from a and b to packed 16-bit integers using signed saturation, and store the results in dst.
    z = _mm256_sub_epi16(z, _mm256_cmpgt_epi16(ymm_computed, _mm256_setzero_si256())); //coeffTemp = (coeff[coeffLocation] < 0)? -coeffTemp : coeffTemp;

    ymm_computed = _mm256_sign_epi16(ymm_computed, coeffTempORG);// coeffTemp);

    a = _mm256_extracti128_si256(ymm_computed, 0);
    b = _mm256_extracti128_si256(ymm_computed, 1);
    _mm_storel_epi64((__m128i *)coeff, a);
    _mm_storel_epi64((__m128i *)(coeff + coeff_stride), _mm_srli_si128(a, 8));
    _mm_storel_epi64((__m128i *)(coeff + 2 * coeff_stride), b);
    _mm_storel_epi64((__m128i *)(coeff + 3 * coeff_stride), _mm_srli_si128(b, 8));

    z = _mm256_sad_epu8(z, _mm256_srli_si256(z, 8));
    *nonzerocoeff = _mm_cvtsi128_si32(_mm_add_epi32(_mm256_extracti128_si256(z, 0), _mm256_extracti128_si256(z, 1)));
}
/*******************************************mat_mult8x8_avx2_intrin**************************************************/
void mat_mult8x8_avx2_intrin(
    int16_t*              coeff,
    const uint32_t         coeff_stride,
    const uint16_t        *maskingMatrix,
    const uint32_t         maskingMatrixStride,  //Matrix size
    const uint32_t         compute_size,  //Computation area size
    const int32_t         offset,     //(PMP_MAX >> 1)
    const int32_t         shift_num, //PMP_PRECISION
    uint32_t*              nonzerocoeff)
{
    unsigned row;
    __m256i z = _mm256_setzero_si256();
    //__m128i a, b;
    __m256i coeffTemp, a0, a1, b0, b1, ymm_computed, MaskingMatrix, coeffTempORG;

    coeffTemp = a0 = a1 = b0 = b1 = ymm_computed = MaskingMatrix = _mm256_setzero_si256();

    // prepare Shift REG
    __m128i PMP_PRECISION_REG = _mm_set_epi16(0, 0, 0, 0, 0, 0, 0, (int16_t)shift_num);//_mm_set1_epi32(shift_num);

    //prepare the offset
    __m256i offsetREG = _mm256_set1_epi32(offset);
    row = 0;
    do {
        //load maskingMatrix_new
        MaskingMatrix = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_loadu_si128((__m128i*)(maskingMatrix + maskingMatrixStride * row))), _mm_loadu_si128((__m128i*)(maskingMatrix + maskingMatrixStride * (row + 1))), 0x1);

        //load coefftemp
        coeffTemp = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_loadu_si128((__m128i*)(coeff + coeff_stride * row))),
            _mm_loadu_si128((__m128i*)(coeff + coeff_stride * (row + 1))), 0x1);

        coeffTempORG = coeffTemp;
        //Absolute val
        coeffTemp = _mm256_abs_epi16(coeffTemp);

        //Multiply
        a0 = _mm256_mullo_epi16(coeffTemp, MaskingMatrix);
        a1 = _mm256_mulhi_epi16(coeffTemp, MaskingMatrix);

        b0 = _mm256_unpacklo_epi16(a0, a1);
        b1 = _mm256_unpackhi_epi16(a0, a1);

        //Add
        b0 = _mm256_add_epi32(b0, offsetREG);
        b1 = _mm256_add_epi32(b1, offsetREG);

        //Shift right by PMP_PRECISION_REG
        b0 = _mm256_sra_epi32(b0, PMP_PRECISION_REG);
        b1 = _mm256_sra_epi32(b1, PMP_PRECISION_REG);

        //coefftemp in c
        ymm_computed = _mm256_packs_epi32(b0, b1);//Convert packed 32-bit integers from b0 and b1 to packed 16-bit integers using signed saturation, and store the results in dst.
        z = _mm256_sub_epi16(z, _mm256_cmpgt_epi16(ymm_computed, _mm256_setzero_si256())); //coeffTemp = (coeff[coeffLocation] < 0)? -coeffTemp : coeffTemp;

        ymm_computed = _mm256_sign_epi16(ymm_computed, coeffTempORG);// coeffTemp);

        _mm_storeu_si128((__m128i *)(coeff + coeff_stride * row), _mm256_extracti128_si256(ymm_computed, 0));
        _mm_storeu_si128((__m128i *)(coeff + coeff_stride * (row + 1)), _mm256_extracti128_si256(ymm_computed, 1));

        row += 2;
    } while (row < compute_size);

    z = _mm256_sad_epu8(z, _mm256_srli_si256(z, 7));
    *nonzerocoeff = _mm_cvtsi128_si32(_mm_add_epi32(_mm256_extracti128_si256(z, 0), _mm256_extracti128_si256(z, 1)));
}
/***************************************mat_mult_nxn_avx2_intrin****************************************************/
void mat_mult_nxn_avx2_intrin(
    int16_t*              coeff,
    const uint32_t         coeff_stride,
    const uint16_t        *maskingMatrix,
    const uint32_t         maskingMatrixStride,  //Matrix size
    const uint32_t         compute_size,  //Computation area size
    const int32_t         offset,     //(PMP_MAX >> 1)
    const int32_t         shift_num, //PMP_PRECISION
    uint32_t*              nonzerocoeff)
{
    unsigned row, col;
    __m256i z = _mm256_setzero_si256();
    //__m128i a, b;
    __m256i coeffTemp, a0, a1, b0, b1, ymm_computed, MaskingMatrix, coeffTempORG;

    coeffTemp = a0 = a1 = b0 = b1 = ymm_computed = MaskingMatrix = _mm256_setzero_si256();

    // prepare Shift REG
    __m128i PMP_PRECISION_REG = _mm_set_epi16(0, 0, 0, 0, 0, 0, 0, (int16_t)shift_num);//_mm_set1_epi32(shift_num);

    //prepare the offset
    __m256i offsetREG = _mm256_set1_epi32(offset);

    row = 0;
    do {
        col = 0;
        do {
            //load coefftemp
            coeffTemp = _mm256_loadu_si256((__m256i *)(coeff + coeff_stride * row + col));

            //load maskingMatrix_new
            MaskingMatrix = _mm256_loadu_si256((__m256i *) (maskingMatrix + maskingMatrixStride * row + col));

            coeffTempORG = coeffTemp;

            //Absolute val
            coeffTemp = _mm256_abs_epi16(coeffTemp);

            //Multiply
            a0 = _mm256_mullo_epi16(coeffTemp, MaskingMatrix);
            a1 = _mm256_mulhi_epi16(coeffTemp, MaskingMatrix);

            b0 = _mm256_unpacklo_epi16(a0, a1);
            b1 = _mm256_unpackhi_epi16(a0, a1);

            //Add
            b0 = _mm256_add_epi32(b0, offsetREG);
            b1 = _mm256_add_epi32(b1, offsetREG);

            //Shift right by PMP_PRECISION_REG
            b0 = _mm256_sra_epi32(b0, PMP_PRECISION_REG);
            b1 = _mm256_sra_epi32(b1, PMP_PRECISION_REG);

            //coefftemp in c
            ymm_computed = _mm256_packs_epi32(b0, b1);//Convert packed 32-bit integers from b0 and b1 to packed 16-bit integers using signed saturation, and store the results in dst.
            z = _mm256_sub_epi16(z, _mm256_cmpgt_epi16(ymm_computed, _mm256_setzero_si256())); //coeffTemp = (coeff[coeffLocation] < 0)? -coeffTemp : coeffTemp;

            ymm_computed = _mm256_sign_epi16(ymm_computed, coeffTempORG);// coeffTemp);

            _mm256_storeu_si256((__m256i *)(coeff + coeff_stride * row + col), ymm_computed);

            col += 16;
        } while (col < compute_size);

        row++;
    } while (row < compute_size);

    z = _mm256_sad_epu8(z, _mm256_srli_si256(z, 7));
    *nonzerocoeff = _mm_cvtsi128_si32(_mm_add_epi32(_mm256_extracti128_si256(z, 0), _mm256_extracti128_si256(z, 1)));
}

static INLINE void EnergyComputation_kernel_avx2(const int32_t *const in,
    __m256i *const sum256)
{
    const __m256i zero = _mm256_setzero_si256();
    const __m256i input = _mm256_load_si256((__m256i *)in);
    const __m256i in_lo = _mm256_unpacklo_epi32(input, zero);
    const __m256i in_hi = _mm256_unpackhi_epi32(input, zero);
    const __m256i energy_lo = _mm256_mul_epi32(in_lo, in_lo);
    const __m256i energy_hi = _mm256_mul_epi32(in_hi, in_hi);
    *sum256 = _mm256_add_epi64(*sum256, energy_lo);
    *sum256 = _mm256_add_epi64(*sum256, energy_hi);
}

static INLINE uint64_t hadd64_avx2(const __m256i sum256) {
    const __m128i sum256_lo = _mm256_castsi256_si128(sum256);
    const __m128i sum256_hi = _mm256_extracti128_si256(sum256, 1);
    const __m128i sum128 = _mm_add_epi64(sum256_lo, sum256_hi);
    const __m128i sum128_hi = _mm_srli_si128(sum128, 8);
    const __m128i sum = _mm_add_epi64(sum128, sum128_hi);

    return _mm_extract_epi64(sum, 0);
}

static INLINE uint64_t EnergyComputation_avx2(const int32_t *const in,
    const uint32_t size)
{
    const __m256i zero = _mm256_setzero_si256();
    uint32_t i = 0;
    __m256i sum = zero;

    do {
        EnergyComputation_kernel_avx2(in + i, &sum);
        i += 8;
    } while (i < size);

    return hadd64_avx2(sum);
}

static INLINE uint64_t EnergyComputation64_avx2(const int32_t *in,
    const uint32_t height)
{
    const __m256i zero = _mm256_setzero_si256();
    uint32_t i = height;
    __m256i sum = zero;

    do {
        EnergyComputation_kernel_avx2(in + 0 * 8, &sum);
        EnergyComputation_kernel_avx2(in + 1 * 8, &sum);
        EnergyComputation_kernel_avx2(in + 2 * 8, &sum);
        EnergyComputation_kernel_avx2(in + 3 * 8, &sum);
        in += 64;
    } while (--i);

    return hadd64_avx2(sum);
}

static INLINE void clean_256_bytes_avx2(int32_t *buf, const uint32_t height) {
    const __m256i zero = _mm256_setzero_si256();
    uint32_t h = height;

    do {
        _mm256_store_si256((__m256i *)(buf + 0 * 8), zero);
        _mm256_store_si256((__m256i *)(buf + 1 * 8), zero);
        _mm256_store_si256((__m256i *)(buf + 2 * 8), zero);
        _mm256_store_si256((__m256i *)(buf + 3 * 8), zero);
        buf += 64;
    } while (--h);
}

static INLINE void copy_32_bytes_avx2(const int32_t *src, int32_t *dst) {
    const __m256i val = _mm256_load_si256((__m256i *)(src + 0 * 8));
    _mm256_store_si256((__m256i *)(dst + 0 * 8), val);
}

static INLINE void copy_256x_bytes_avx2(const int32_t *src, int32_t *dst,
    const uint32_t height) {
    uint32_t h = height;

    do {
        copy_32_bytes_avx2(src + 0 * 8, dst + 0 * 8);
        copy_32_bytes_avx2(src + 1 * 8, dst + 1 * 8);
        copy_32_bytes_avx2(src + 2 * 8, dst + 2 * 8);
        copy_32_bytes_avx2(src + 3 * 8, dst + 3 * 8);
        src += 64;
        dst += 32;
    } while (--h);
}

uint64_t HandleTransform16x64_avx2(int32_t *output) {
    //bottom 16x32 area.
    const uint64_t three_quad_energy =
        EnergyComputation_avx2(output + 16 * 32, 16 * 32);

    // Zero out the bottom 16x32 area.
    memset(output + 16 * 32, 0, 16 * 32 * sizeof(*output));

    return three_quad_energy;
}

uint64_t HandleTransform32x64_avx2(int32_t *output) {
    //bottom 32x32 area.
    const uint64_t three_quad_energy =
        EnergyComputation_avx2(output + 32 * 32, 32 * 32);

    // Zero out the bottom 32x32 area.
    memset(output + 32 * 32, 0, 32 * 32 * sizeof(*output));

    return three_quad_energy;
}

uint64_t HandleTransform64x16_avx2(int32_t *output) {
    // top - right 32x16 area.
    const uint64_t three_quad_energy =
        EnergyComputation64_avx2(output + 32, 16);

    // Zero out right 32x16 area.
    clean_256_bytes_avx2(output + 32, 16);

    // Re-pack non-zero coeffs in the first 32x16 indices.
    copy_256x_bytes_avx2(output + 64, output + 32, 15);

    return three_quad_energy;
}

uint64_t HandleTransform64x32_avx2(int32_t *output) {
    // top - right 32x32 area.
    const uint64_t three_quad_energy =
        EnergyComputation64_avx2(output + 32, 32);

    // Zero out right 32x32 area.
    clean_256_bytes_avx2(output + 32, 32);

    // Re-pack non-zero coeffs in the first 32x32 indices.
    copy_256x_bytes_avx2(output + 64, output + 32, 31);

    return three_quad_energy;
}

uint64_t HandleTransform64x64_avx2(int32_t *output) {
    uint64_t three_quad_energy;

    // top - right 32x32 area.
    three_quad_energy = EnergyComputation64_avx2(output + 32, 32);
    //bottom 64x32 area.
    three_quad_energy += EnergyComputation_avx2(output + 32 * 64, 64 * 32);

    // Zero out top-right 32x32 area.
    clean_256_bytes_avx2(output + 32, 32);

    // Zero out the bottom 64x32 area.
    memset(output + 32 * 64, 0, 32 * 64 * sizeof(*output));

    // Re-pack non-zero coeffs in the first 32x32 indices.
    copy_256x_bytes_avx2(output + 64, output + 32, 31);

    return three_quad_energy;
}
