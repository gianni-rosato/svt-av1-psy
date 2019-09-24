/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbDefinitions.h"

#include "emmintrin.h"

#ifndef PREFETCH
#define PREFETCH 0 // prefetching: enables prefetching of data before interpolation
#endif

#include "tmmintrin.h"

#ifdef __GNUC__
#ifndef __cplusplus
__attribute__((visibility("hidden")))
#endif
#endif
const int16_t lumaFilterCoeff[4][8] =
{
  { 0, 0,  0, 64,  0,  0, 0,  0},
  {-1, 4,-10, 58, 17, -5, 1,  0},
  {-1, 4,-11, 40, 40,-11, 4, -1},
  { 0, 1, -5, 17, 58,-10, 4, -1}
};

#ifdef __GNUC__
#ifndef __cplusplus
__attribute__((visibility("hidden")))
#endif
#endif
const int16_t lumaFilterCoeff7[4][8] =
{
  { 0, 0,  0, 64,  0,  0, 0,  0},
  {-1, 4,-10, 58, 17, -5, 1,  0},
  {-1, 4,-11, 40, 40,-11, 4, -1},
  { 1, -5, 17, 58,-10, 4, -1, 0}
};

#ifdef __GNUC__
#ifndef __cplusplus
__attribute__((visibility("hidden")))
#endif
#endif
const int16_t chromaFilterCoeff[8][4] =
{
  { 0, 64,  0,  0},
  {-2, 58, 10, -2},
  {-4, 54, 16, -2},
  {-6, 46, 28, -4},
  {-4, 36, 36, -4},
  {-4, 28, 46, -6},
  {-2, 16, 54, -4},
  {-2, 10, 58, -2},
};

static void PrefetchBlock(uint8_t *src, uint32_t src_stride, uint32_t blkWidth, uint32_t blkHeight)
{
#if PREFETCH
    uint32_t row_count = blkHeight;

    do {
        uint8_t *addr0 = src;
        uint8_t *addr1 = addr0 + blkWidth - 1;
        src += src_stride;

        _mm_prefetch((char*)addr0, _MM_HINT_T0);
        _mm_prefetch((char*)addr1, _MM_HINT_T0);
    } while (--row_count != 0);
#else
    (void)src;
    (void)src_stride;
    (void)blkWidth;
    (void)blkHeight;
#endif
}

void LumaInterpolationFilterTwoDInRaw7_SSSE3(int16_t *first_pass_if_dst, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, uint32_t frac_pos_y)
{
    int32_t row_count, col_count;
    __m128i c0, c1, c2;
    __m128i a0, a1, a2, a3, a4, a5, a6;
    __m128i sum0, sum1;
    __m128i b0l, b0h, b1l, b1h, b2l, b2h;

    EbByte qtr;

    c0 = _mm_loadu_si128((__m128i *)lumaFilterCoeff7[frac_pos_y]);
    c2 = _mm_shuffle_epi32(c0, 0xaa);
    c1 = _mm_shuffle_epi32(c0, 0x55);
    c0 = _mm_shuffle_epi32(c0, 0x00);

    if (pu_width & 4)
    {
        row_count = pu_height;

        qtr = dst;

        do {
            a0 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 0 * 4));
            a1 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 1 * 4));
            a2 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 2 * 4));
            a3 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 3 * 4));
            a4 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 4 * 4));
            a5 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 5 * 4));
            a6 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 6 * 4));
            a0 = _mm_sub_epi16(a0, a6);

            sum0 = _mm_set1_epi32(257 << 11);
            sum1 = _mm_set1_epi32(257 << 11);

            b0l = _mm_unpacklo_epi16(a0, a1);
            b0h = _mm_unpackhi_epi16(a0, a1);
            b1l = _mm_unpacklo_epi16(a2, a3);
            b1h = _mm_unpackhi_epi16(a2, a3);
            b2l = _mm_unpacklo_epi16(a4, a5);
            b2h = _mm_unpackhi_epi16(a4, a5);

            sum0 = _mm_add_epi32(sum0, _mm_madd_epi16(b0l, c0));
            sum1 = _mm_add_epi32(sum1, _mm_madd_epi16(b0h, c0));
            sum0 = _mm_add_epi32(sum0, _mm_madd_epi16(b1l, c1));
            sum1 = _mm_add_epi32(sum1, _mm_madd_epi16(b1h, c1));
            sum0 = _mm_add_epi32(sum0, _mm_madd_epi16(b2l, c2));
            sum1 = _mm_add_epi32(sum1, _mm_madd_epi16(b2h, c2));

            sum0 = _mm_srai_epi32(sum0, 12);
            sum1 = _mm_srai_epi32(sum1, 12);
            sum0 = _mm_packs_epi32(sum0, sum1);
            sum0 = _mm_packus_epi16(sum0, sum0);

            *(uint32_t *)qtr = _mm_cvtsi128_si32(sum0); qtr += dst_stride;
            *(uint32_t *)qtr = _mm_cvtsi128_si32(_mm_srli_si128(sum0, 4)); qtr += dst_stride;

            first_pass_if_dst += 8;
            row_count -= 2;
        } while (row_count > 0);

        pu_width -= 4;
        if (pu_width == 0)
            return;
        first_pass_if_dst += (frac_pos_y == 2) ? 32 : 24;
        dst += 4;
    }

    col_count = pu_width;
    do {
        EbByte qtr = dst;

        row_count = pu_height;
        do {
            a0 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 0 * 8));
            a1 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 1 * 8));
            a2 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 2 * 8));
            a3 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 3 * 8));
            a4 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 4 * 8));
            a5 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 5 * 8));
            a6 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 6 * 8));
            a0 = _mm_sub_epi16(a0, a6);

            sum0 = _mm_set1_epi32(257 << 11);
            sum1 = _mm_set1_epi32(257 << 11);
            sum0 = _mm_add_epi32(sum0, _mm_madd_epi16(_mm_unpacklo_epi16(a0, a1), c0));
            sum1 = _mm_add_epi32(sum1, _mm_madd_epi16(_mm_unpackhi_epi16(a0, a1), c0));
            sum0 = _mm_add_epi32(sum0, _mm_madd_epi16(_mm_unpacklo_epi16(a2, a3), c1));
            sum1 = _mm_add_epi32(sum1, _mm_madd_epi16(_mm_unpackhi_epi16(a2, a3), c1));
            sum0 = _mm_add_epi32(sum0, _mm_madd_epi16(_mm_unpacklo_epi16(a4, a5), c2));
            sum1 = _mm_add_epi32(sum1, _mm_madd_epi16(_mm_unpackhi_epi16(a4, a5), c2));

            sum0 = _mm_srai_epi32(sum0, 12);
            sum1 = _mm_srai_epi32(sum1, 12);
            sum0 = _mm_packs_epi32(sum0, sum1);
            sum0 = _mm_packus_epi16(sum0, sum0);

            _mm_storel_epi64((__m128i *)qtr, sum0); qtr += dst_stride;

            first_pass_if_dst += 8;
            row_count--;
        } while (row_count > 0);

        first_pass_if_dst += (frac_pos_y == 2) ? 56 : 48;
        dst += 8;
        col_count -= 8;
    } while (col_count > 0);
}

void LumaInterpolationFilterTwoDInRawOutRaw7_SSSE3(int16_t *first_pass_if_dst, int16_t *dst, uint32_t pu_width, uint32_t pu_height, uint32_t frac_pos_y)
{
    int32_t row_count, col_count;

    __m128i a0, a1, a2, a3, a4, a5, a6;
    __m128i c0, c1, c2;
    c0 = _mm_loadu_si128((__m128i *)lumaFilterCoeff7[frac_pos_y]);
    c2 = _mm_shuffle_epi32(c0, 0xaa);
    c1 = _mm_shuffle_epi32(c0, 0x55);
    c0 = _mm_shuffle_epi32(c0, 0x00);

    if (pu_width & 4)
    {
        row_count = pu_height;

        do {
            __m128i sum0, sum1;
            a0 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 0 * 4));
            a1 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 1 * 4));
            a2 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 2 * 4));
            a3 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 3 * 4));
            a4 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 4 * 4));
            a5 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 5 * 4));
            a6 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 6 * 4));
            a0 = _mm_sub_epi16(a0, a6);

            sum0 = _mm_madd_epi16(_mm_unpacklo_epi16(a0, a1), c0);
            sum1 = _mm_madd_epi16(_mm_unpackhi_epi16(a0, a1), c0);
            sum0 = _mm_add_epi32(sum0, _mm_madd_epi16(_mm_unpacklo_epi16(a2, a3), c1));
            sum1 = _mm_add_epi32(sum1, _mm_madd_epi16(_mm_unpackhi_epi16(a2, a3), c1));
            sum0 = _mm_add_epi32(sum0, _mm_madd_epi16(_mm_unpacklo_epi16(a4, a5), c2));
            sum1 = _mm_add_epi32(sum1, _mm_madd_epi16(_mm_unpackhi_epi16(a4, a5), c2));

            sum0 = _mm_srai_epi32(sum0, 6);
            sum1 = _mm_srai_epi32(sum1, 6);
            sum0 = _mm_packs_epi32(sum0, sum1);

            _mm_storeu_si128((__m128i *)dst, sum0);
            dst += 8;

            first_pass_if_dst += 8;
            row_count -= 2;
        } while (row_count > 0);

        pu_width -= 4;
        if (pu_width == 0)
            return;
        first_pass_if_dst += (frac_pos_y == 2) ? 32 : 24;
    }

    col_count = pu_width;
    do {
        row_count = pu_height;
        do {
            __m128i b0l, b0h, b1l, b1h, b2l, b2h;
            __m128i sum0, sum1;

            a0 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 0 * 8));
            a1 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 1 * 8));
            a2 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 2 * 8));
            a3 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 3 * 8));
            a4 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 4 * 8));
            a5 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 5 * 8));
            a6 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 6 * 8));
            a0 = _mm_sub_epi16(a0, a6);

            b0l = _mm_unpacklo_epi16(a0, a1);
            b0h = _mm_unpackhi_epi16(a0, a1);
            b1l = _mm_unpacklo_epi16(a2, a3);
            b1h = _mm_unpackhi_epi16(a2, a3);
            b2l = _mm_unpacklo_epi16(a4, a5);
            b2h = _mm_unpackhi_epi16(a4, a5);

            sum0 = _mm_madd_epi16(b0l, c0);
            sum1 = _mm_madd_epi16(b0h, c0);
            sum0 = _mm_add_epi32(sum0, _mm_madd_epi16(b1l, c1));
            sum1 = _mm_add_epi32(sum1, _mm_madd_epi16(b1h, c1));
            sum0 = _mm_add_epi32(sum0, _mm_madd_epi16(b2l, c2));
            sum1 = _mm_add_epi32(sum1, _mm_madd_epi16(b2h, c2));

            sum0 = _mm_srai_epi32(sum0, 6);
            sum1 = _mm_srai_epi32(sum1, 6);
            sum0 = _mm_packs_epi32(sum0, sum1);

            _mm_storeu_si128((__m128i *)dst, sum0);
            dst += 8;

            first_pass_if_dst += 8;
            row_count--;
        } while (row_count > 0);

        first_pass_if_dst += (frac_pos_y == 2) ? 56 : 48;
        col_count -= 8;
    } while (col_count > 0);
}

void LumaInterpolationFilterTwoDInRawM_SSSE3(int16_t *first_pass_if_dst, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height)
{
    int32_t row_count, col_count;

    __m128i c0, c1;
    __m128i a0, a1, a2, a3, a4, a5, a6, a7;
    __m128i sum0, sum1;

    EbByte qtr;

    c0 = _mm_loadu_si128((__m128i *)lumaFilterCoeff7[2]);
    c1 = _mm_shuffle_epi32(c0, 0x55);
    c0 = _mm_shuffle_epi32(c0, 0x00);

    if (pu_width & 4){
        row_count = pu_height;
        qtr = dst;

        do {
            a0 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 0 * 4));
            a1 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 1 * 4));
            a2 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 2 * 4));
            a3 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 3 * 4));
            a4 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 4 * 4));
            a5 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 5 * 4));
            a6 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 6 * 4));
            a7 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 7 * 4));

            sum0 = _mm_set1_epi32(257 << 11);
            sum1 = _mm_set1_epi32(257 << 11);

            a0 = _mm_add_epi16(a0, a7);
            a1 = _mm_add_epi16(a1, a6);
            a2 = _mm_add_epi16(a2, a5);
            a3 = _mm_add_epi16(a3, a4);
            sum0 = _mm_add_epi32(sum0, _mm_madd_epi16(_mm_unpacklo_epi16(a0, a1), c0));
            sum1 = _mm_add_epi32(sum1, _mm_madd_epi16(_mm_unpackhi_epi16(a0, a1), c0));
            sum0 = _mm_add_epi32(sum0, _mm_madd_epi16(_mm_unpacklo_epi16(a2, a3), c1));
            sum1 = _mm_add_epi32(sum1, _mm_madd_epi16(_mm_unpackhi_epi16(a2, a3), c1));

            sum0 = _mm_srai_epi32(sum0, 12);
            sum1 = _mm_srai_epi32(sum1, 12);
            sum0 = _mm_packs_epi32(sum0, sum1);
            sum0 = _mm_packus_epi16(sum0, sum0);

            *(uint32_t *)qtr = _mm_cvtsi128_si32(sum0); qtr += dst_stride;
            *(uint32_t *)qtr = _mm_cvtsi128_si32(_mm_srli_si128(sum0, 4)); qtr += dst_stride;
            first_pass_if_dst += 8;
            row_count -= 2;
        } while (row_count > 0);

        pu_width -= 4;
        if (pu_width == 0)
            return;
        first_pass_if_dst += 32;
        dst += 4;
    }

    col_count = pu_width;
    do {
        qtr = dst;

        row_count = pu_height;
        do {
            a0 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 0 * 8));
            a1 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 1 * 8));
            a2 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 2 * 8));
            a3 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 3 * 8));
            a4 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 4 * 8));
            a5 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 5 * 8));
            a6 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 6 * 8));
            a7 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 7 * 8));

            sum0 = _mm_set1_epi32(257 << 11);
            sum1 = _mm_set1_epi32(257 << 11);
            a0 = _mm_add_epi16(a0, a7);
            a1 = _mm_add_epi16(a1, a6);
            a2 = _mm_add_epi16(a2, a5);
            a3 = _mm_add_epi16(a3, a4);
            sum0 = _mm_add_epi32(sum0, _mm_madd_epi16(_mm_unpacklo_epi16(a0, a1), c0));
            sum1 = _mm_add_epi32(sum1, _mm_madd_epi16(_mm_unpackhi_epi16(a0, a1), c0));
            sum0 = _mm_add_epi32(sum0, _mm_madd_epi16(_mm_unpacklo_epi16(a2, a3), c1));
            sum1 = _mm_add_epi32(sum1, _mm_madd_epi16(_mm_unpackhi_epi16(a2, a3), c1));

            sum0 = _mm_srai_epi32(sum0, 12);
            sum1 = _mm_srai_epi32(sum1, 12);
            sum0 = _mm_packs_epi32(sum0, sum1);
            sum0 = _mm_packus_epi16(sum0, sum0);

            _mm_storel_epi64((__m128i *)qtr, sum0); qtr += dst_stride;
            first_pass_if_dst += 8;
        } while (--row_count > 0);

        first_pass_if_dst += 56;
        dst += 8;
        col_count -= 8;
    } while (col_count > 0);
}

void LumaInterpolationFilterTwoDInRawOutRawM_SSSE3(int16_t *first_pass_if_dst, int16_t *dst, uint32_t pu_width, uint32_t pu_height){
    int32_t row_count, col_count;

    __m128i a0, a1, a2, a3, a4, a5, a6, a7;
    __m128i c0, c1;
    c0 = _mm_loadu_si128((__m128i *)lumaFilterCoeff7[2]);
    c1 = _mm_shuffle_epi32(c0, 0x55);
    c0 = _mm_shuffle_epi32(c0, 0x00);

    if (pu_width & 4) {
        row_count = pu_height;

        do {
            __m128i sum0, sum1;
            a0 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 0 * 4));
            a1 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 1 * 4));
            a2 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 2 * 4));
            a3 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 3 * 4));
            a4 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 4 * 4));
            a5 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 5 * 4));
            a6 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 6 * 4));
            a7 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 7 * 4));

            a0 = _mm_add_epi16(a0, a7);
            a1 = _mm_add_epi16(a1, a6);
            a2 = _mm_add_epi16(a2, a5);
            a3 = _mm_add_epi16(a3, a4);
            sum0 = _mm_madd_epi16(_mm_unpacklo_epi16(a0, a1), c0);
            sum1 = _mm_madd_epi16(_mm_unpackhi_epi16(a0, a1), c0);
            sum0 = _mm_add_epi32(sum0, _mm_madd_epi16(_mm_unpacklo_epi16(a2, a3), c1));
            sum1 = _mm_add_epi32(sum1, _mm_madd_epi16(_mm_unpackhi_epi16(a2, a3), c1));

            sum0 = _mm_srai_epi32(sum0, 6);
            sum1 = _mm_srai_epi32(sum1, 6);
            sum0 = _mm_packs_epi32(sum0, sum1);

            _mm_storeu_si128((__m128i *)dst, sum0);
            dst += 8;
            first_pass_if_dst += 8;
            row_count -= 2;
        } while (row_count > 0);

        pu_width -= 4;
        if (pu_width == 0)
            return;
        first_pass_if_dst += 32;
    }

    col_count = pu_width;
    do {
        row_count = pu_height;
        do {
            __m128i sum0, sum1;
            a0 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 0 * 8));
            a1 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 1 * 8));
            a2 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 2 * 8));
            a3 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 3 * 8));
            a4 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 4 * 8));
            a5 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 5 * 8));
            a6 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 6 * 8));
            a7 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 7 * 8));

            a0 = _mm_add_epi16(a0, a7);
            a1 = _mm_add_epi16(a1, a6);
            a2 = _mm_add_epi16(a2, a5);
            a3 = _mm_add_epi16(a3, a4);
            sum0 = _mm_madd_epi16(_mm_unpacklo_epi16(a0, a1), c0);
            sum1 = _mm_madd_epi16(_mm_unpackhi_epi16(a0, a1), c0);
            sum0 = _mm_add_epi32(sum0, _mm_madd_epi16(_mm_unpacklo_epi16(a2, a3), c1));
            sum1 = _mm_add_epi32(sum1, _mm_madd_epi16(_mm_unpackhi_epi16(a2, a3), c1));

            sum0 = _mm_srai_epi32(sum0, 6);
            sum1 = _mm_srai_epi32(sum1, 6);
            sum0 = _mm_packs_epi32(sum0, sum1);

            _mm_storeu_si128((__m128i *)dst, sum0);
            dst += 8;
            first_pass_if_dst += 8;
        } while (--row_count > 0);

        first_pass_if_dst += 56;
        col_count -= 8;
    } while (col_count > 0);
}

void PictureCopyKernelOutRaw_SSSE3(
    EbByte                  ref_pic,
    uint32_t                   src_stride,
    int16_t                  *dst,
    uint32_t                   pu_width,
    uint32_t                   pu_height,
    int16_t                   offset)
{
    uint32_t row_count, col_count;
    __m128i o;

    PrefetchBlock(ref_pic, src_stride, pu_width, pu_height);

    /*__m128i*/ o = _mm_set1_epi16(offset);

    if (pu_width & 2) {
        __m128i a0;
        EbByte ptr = ref_pic;
        row_count = pu_height;
        /*__m128i*/ a0 = _mm_setzero_si128();
        do {
            a0 = _mm_insert_epi16(a0, *(uint16_t *)ptr, 0); ptr += src_stride;
            a0 = _mm_insert_epi16(a0, *(uint16_t *)ptr, 1); ptr += src_stride;
            a0 = _mm_insert_epi16(a0, *(uint16_t *)ptr, 2); ptr += src_stride;
            a0 = _mm_insert_epi16(a0, *(uint16_t *)ptr, 3); ptr += src_stride;
            a0 = _mm_unpacklo_epi8(a0, _mm_setzero_si128());
            a0 = _mm_slli_epi16(a0, 6);
            a0 = _mm_sub_epi16(a0, o);
            _mm_storeu_si128((__m128i *)dst, a0);

            dst += 8;
            row_count -= 4;
        } while (row_count != 0);

        pu_width -= 2;
        if (pu_width == 0)
            return;
        ref_pic += 2;
    }

    if (pu_width & 4) {
        EbByte ptr = ref_pic;
        row_count = pu_height;
        do {
            __m128i a0, a1;
            a0 = _mm_cvtsi32_si128(*(uint32_t *)ptr); ptr += src_stride;
            a1 = _mm_cvtsi32_si128(*(uint32_t *)ptr); ptr += src_stride;
            a0 = _mm_unpacklo_epi32(a0, a1);
            a0 = _mm_unpacklo_epi8(a0, _mm_setzero_si128());
            a0 = _mm_slli_epi16(a0, 6);
            a0 = _mm_sub_epi16(a0, o);
            _mm_storeu_si128((__m128i *)dst, a0);

            dst += 8;
            row_count -= 2;
        } while (row_count != 0);

        pu_width -= 4;
        if (pu_width == 0)
            return;
        ref_pic += 4;
    }

    col_count = pu_width;
    do {
        __m128i a0;
        EbByte ptr = ref_pic;
        row_count = pu_height;
        do {
            /*__m128i*/ a0 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
            a0 = _mm_unpacklo_epi8(a0, _mm_setzero_si128());
            a0 = _mm_slli_epi16(a0, 6);
            a0 = _mm_sub_epi16(a0, o);
            _mm_storeu_si128((__m128i *)dst, a0);
            dst += 8;
        } while (--row_count != 0);

        col_count -= 8;
        ref_pic += 8;
    } while (col_count != 0);
}
