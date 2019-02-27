/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/


#include "EbMcp_SSSE3.h"
#include "EbDefinitions.h"

#include "emmintrin.h"

#ifndef PREFETCH
#define PREFETCH 0 // prefetching: enables prefetching of data before interpolation
#endif

#include "tmmintrin.h"
#include "smmintrin.h"
// Note: _mm_extract_epi32 & _mm_extract_epi64 are SSE4 functions


#if defined(__linux__) || defined(__APPLE__)
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

#if defined(__linux__) || defined(__APPLE__)
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

#if defined(__linux__) || defined(__APPLE__)
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

    do
    {
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

void PictureCopyKernel_SSSE3(
    EbByte                  src,
    uint32_t                   src_stride,
    EbByte                  dst,
    uint32_t                   dst_stride,
    uint32_t                   area_width,
    uint32_t                   area_height,
    uint32_t                   bytes_per_sample)
{

    uint32_t row_count, colCount;
    (void)bytes_per_sample;


    PrefetchBlock(src, src_stride, area_width, area_height);



    if (area_width & 2)
    {
        EbByte ptr = src;
        EbByte qtr = dst;
        row_count = area_height;
        do
        {
            // Note: do only one at a time to minimize number of registers used
            uint16_t a0;
            a0 = *(uint16_t *)ptr; ptr += src_stride;
            *(uint16_t *)qtr = a0; qtr += dst_stride;
            row_count--;
        } while (row_count != 0);
        area_width -= 2;
        if (area_width == 0)
        {
            return;
        }
        src += 2;
        dst += 2;
    }

    if (area_width & 4)
    {
        EbByte ptr = src;
        EbByte qtr = dst;
        row_count = area_height;
        do
        {
            __m128i a0, a1;
            a0 = _mm_cvtsi32_si128(*(uint32_t *)ptr); ptr += src_stride;
            a1 = _mm_cvtsi32_si128(*(uint32_t *)ptr); ptr += src_stride;
            *(uint32_t *)qtr = _mm_cvtsi128_si32(a0); qtr += dst_stride;
            *(uint32_t *)qtr = _mm_cvtsi128_si32(a1); qtr += dst_stride;
            row_count -= 2;
        } while (row_count != 0);
        area_width -= 4;
        if (area_width == 0)
        {
            return;
        }
        src += 4;
        dst += 4;
    }

    if (area_width & 8)
    {
        EbByte ptr = src;
        EbByte qtr = dst;
        row_count = area_height;
        do
        {
            __m128i a0, a1;
            a0 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
            a1 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
            _mm_storel_epi64((__m128i *)qtr, a0); qtr += dst_stride;
            _mm_storel_epi64((__m128i *)qtr, a1); qtr += dst_stride;
            row_count -= 2;
        } while (row_count != 0);
        area_width -= 8;
        if (area_width == 0)
        {
            return;
        }
        src += 8;
        dst += 8;
    }

    colCount = area_width;
    do
    {
        EbByte ptr = src;
        EbByte qtr = dst;
        row_count = area_height;
        do
        {
            __m128i a0, a1;
            a0 = _mm_loadu_si128((__m128i *)ptr); ptr += src_stride;
            a1 = _mm_loadu_si128((__m128i *)ptr); ptr += src_stride;
            _mm_storeu_si128((__m128i *)qtr, a0); qtr += dst_stride;
            _mm_storeu_si128((__m128i *)qtr, a1); qtr += dst_stride;
            row_count -= 2;
        } while (row_count != 0);
        colCount -= 16;
        src += 16;
        dst += 16;
    } while (colCount != 0);
}

void chroma_interpolation_copy_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    EbByte               dst,
    uint32_t                dst_stride,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst,
    uint32_t                frac_pos_x,
    uint32_t                frac_pos_y)
{
    (void)first_pass_if_dst;
    (void)frac_pos_x;
    (void)frac_pos_y;
    PictureCopyKernel_SSSE3(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, 1);
}

void luma_interpolation_copy_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    EbByte               dst,
    uint32_t                dst_stride,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst)
{
    (void)first_pass_if_dst;
    PictureCopyKernel_SSSE3(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, 1);
}

void LumaInterpolationFilterTwoDInRaw7_SSSE3(int16_t *first_pass_if_dst, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, uint32_t frac_pos_y)
{
    int32_t row_count, colCount;
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

        do
        {
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

            *(uint32_t *)qtr = _mm_extract_epi32(sum0, 0); qtr += dst_stride;
            *(uint32_t *)qtr = _mm_extract_epi32(sum0, 1); qtr += dst_stride;

            first_pass_if_dst += 8;
            row_count -= 2;
        } while (row_count > 0);

        pu_width -= 4;
        if (pu_width == 0)
        {
            return;
        }

        first_pass_if_dst += (frac_pos_y == 2) ? 32 : 24;
        dst += 4;
    }

    colCount = pu_width;
    do
    {
        EbByte qtr = dst;

        row_count = pu_height;
        do
        {
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
        colCount -= 8;
    } while (colCount > 0);
}

void LumaInterpolationFilterTwoDInRawOutRaw7_SSSE3(int16_t *first_pass_if_dst, int16_t *dst, uint32_t pu_width, uint32_t pu_height, uint32_t frac_pos_y)
{
    int32_t row_count, colCount;

    __m128i a0, a1, a2, a3, a4, a5, a6;
    __m128i c0, c1, c2;
    c0 = _mm_loadu_si128((__m128i *)lumaFilterCoeff7[frac_pos_y]);
    c2 = _mm_shuffle_epi32(c0, 0xaa);
    c1 = _mm_shuffle_epi32(c0, 0x55);
    c0 = _mm_shuffle_epi32(c0, 0x00);

    if (pu_width & 4)
    {
        row_count = pu_height;

        do
        {
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
        {
            return;
        }

        first_pass_if_dst += (frac_pos_y == 2) ? 32 : 24;
    }

    colCount = pu_width;
    do
    {
        row_count = pu_height;
        do
        {
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
        colCount -= 8;
    } while (colCount > 0);
}

void LumaInterpolationFilterTwoDInRawM_SSSE3(int16_t *first_pass_if_dst, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height)
{
    int32_t row_count, colCount;

    __m128i c0, c1;
    __m128i a0, a1, a2, a3, a4, a5, a6, a7;
    __m128i sum0, sum1;

    EbByte qtr;

    c0 = _mm_loadu_si128((__m128i *)lumaFilterCoeff7[2]);
    c1 = _mm_shuffle_epi32(c0, 0x55);
    c0 = _mm_shuffle_epi32(c0, 0x00);



    if (pu_width & 4)
    {
        row_count = pu_height;
        qtr = dst;

        do
        {
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

            *(uint32_t *)qtr = _mm_extract_epi32(sum0, 0); qtr += dst_stride;
            *(uint32_t *)qtr = _mm_extract_epi32(sum0, 1); qtr += dst_stride;
            first_pass_if_dst += 8;
            row_count -= 2;
        } while (row_count > 0);

        pu_width -= 4;
        if (pu_width == 0)
        {
            return;
        }

        first_pass_if_dst += 32;
        dst += 4;
    }

    colCount = pu_width;
    do
    {
        qtr = dst;

        row_count = pu_height;
        do
        {
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
        colCount -= 8;
    } while (colCount > 0);
}

void LumaInterpolationFilterTwoDInRawOutRawM_SSSE3(int16_t *first_pass_if_dst, int16_t *dst, uint32_t pu_width, uint32_t pu_height)
{
    int32_t row_count, colCount;

    __m128i a0, a1, a2, a3, a4, a5, a6, a7;
    __m128i c0, c1;
    c0 = _mm_loadu_si128((__m128i *)lumaFilterCoeff7[2]);
    c1 = _mm_shuffle_epi32(c0, 0x55);
    c0 = _mm_shuffle_epi32(c0, 0x00);

    if (pu_width & 4)
    {
        row_count = pu_height;

        do
        {
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
        {
            return;
        }

        first_pass_if_dst += 32;
    }

    colCount = pu_width;
    do
    {
        row_count = pu_height;
        do
        {
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
        colCount -= 8;
    } while (colCount > 0);
}

void PictureCopyKernelOutRaw_SSSE3(
    EbByte                  ref_pic,
    uint32_t                   src_stride,
    int16_t                  *dst,
    uint32_t                   pu_width,
    uint32_t                   pu_height,
    int16_t                   offset)
{

    uint32_t row_count, colCount;
    __m128i o;

    PrefetchBlock(ref_pic, src_stride, pu_width, pu_height);



    /*__m128i*/ o = _mm_set1_epi16(offset);

    if (pu_width & 2)
    {
        __m128i a0;
        EbByte ptr = ref_pic;
        row_count = pu_height;
        /*__m128i*/ a0 = _mm_setzero_si128();
        do
        {
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
        {
            return;
        }

        ref_pic += 2;
    }

    if (pu_width & 4)
    {
        EbByte ptr = ref_pic;
        row_count = pu_height;
        do
        {
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
        {
            return;
        }

        ref_pic += 4;
    }

    colCount = pu_width;
    do
    {
        __m128i a0;
        EbByte ptr = ref_pic;
        row_count = pu_height;
        do
        {
            /*__m128i*/ a0 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
            a0 = _mm_unpacklo_epi8(a0, _mm_setzero_si128());
            a0 = _mm_slli_epi16(a0, 6);
            a0 = _mm_sub_epi16(a0, o);
            _mm_storeu_si128((__m128i *)dst, a0);
            dst += 8;
        } while (--row_count != 0);

        colCount -= 8;
        ref_pic += 8;
    } while (colCount != 0);
}

void chroma_interpolation_copy_out_raw_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    int16_t               *dst,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst,
    uint32_t                frac_pos_x,
    uint32_t                frac_pos_y)
{
    (void)first_pass_if_dst;
    (void)frac_pos_x;
    (void)frac_pos_y;
    PictureCopyKernelOutRaw_SSSE3(ref_pic, src_stride, dst, pu_width, pu_height, 0);
}
void ChromaInterpolationFilterTwoDInRaw_SSSE3(int16_t *first_pass_if_dst, EbByte dst, uint32_t dst_stride, uint32_t pu_width, uint32_t pu_height, uint32_t frac_pos_y)
{
    uint32_t row_count, colCount;
    __m128i c0, c1; // coeffs
    __m128i offset;
    EbByte qtr;
    __m128i a0, a1, a2, a3, b0, b1, b2, b3;
    __m128i sum0, sum1;

    c0 = _mm_loadl_epi64((__m128i *)chromaFilterCoeff[frac_pos_y]);
    c1 = _mm_shuffle_epi32(c0, 0x55);
    c0 = _mm_shuffle_epi32(c0, 0);
    offset = _mm_set1_epi32(1 << 11);

    if (pu_width & 2)
    {

        row_count = pu_height;
        qtr = dst;
        do
        {
            a0 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 0));
            a1 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 2));
            a2 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 4));
            a3 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 6));

            b0 = _mm_unpacklo_epi16(a0, a1);
            b1 = _mm_unpacklo_epi16(a2, a3);
            b2 = _mm_unpackhi_epi16(a0, a1);
            b3 = _mm_unpackhi_epi16(a2, a3);

            sum0 = _mm_madd_epi16(b0, c0);
            sum0 = _mm_add_epi32(sum0, _mm_madd_epi16(b1, c1));
            sum1 = _mm_madd_epi16(b2, c0);
            sum1 = _mm_add_epi32(sum1, _mm_madd_epi16(b3, c1));
            sum0 = _mm_add_epi32(sum0, offset);
            sum1 = _mm_add_epi32(sum1, offset);
            sum0 = _mm_srai_epi32(sum0, 12);
            sum1 = _mm_srai_epi32(sum1, 12);
            sum0 = _mm_packs_epi32(sum0, sum1);
            sum0 = _mm_packus_epi16(sum0, sum0);

            *(uint16_t *)qtr = (uint16_t)(_mm_extract_epi16(sum0, 0)); qtr += dst_stride;
            *(uint16_t *)qtr = (uint16_t)(_mm_extract_epi16(sum0, 1)); qtr += dst_stride;
            *(uint16_t *)qtr = (uint16_t)(_mm_extract_epi16(sum0, 2)); qtr += dst_stride;
            *(uint16_t *)qtr = (uint16_t)(_mm_extract_epi16(sum0, 3)); qtr += dst_stride;

            first_pass_if_dst += 8;
            row_count -= 4;
        } while (row_count != 0);

        pu_width -= 2;
        if (pu_width == 0)
        {
            return;
        }
        first_pass_if_dst += 8;
        dst += 2;
    }

    if (pu_width & 4)
    {
        row_count = pu_height;
        qtr = dst;
        do
        {
            a0 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 0));
            a1 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 4));
            a2 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 8));
            a3 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 12));
            b0 = _mm_unpacklo_epi16(a0, a1);
            b1 = _mm_unpacklo_epi16(a2, a3);
            b2 = _mm_unpackhi_epi16(a0, a1);
            b3 = _mm_unpackhi_epi16(a2, a3);

            sum0 = _mm_madd_epi16(b0, c0);
            sum0 = _mm_add_epi32(sum0, _mm_madd_epi16(b1, c1));
            sum1 = _mm_madd_epi16(b2, c0);
            sum1 = _mm_add_epi32(sum1, _mm_madd_epi16(b3, c1));
            sum0 = _mm_add_epi32(sum0, offset);
            sum1 = _mm_add_epi32(sum1, offset);
            sum0 = _mm_srai_epi32(sum0, 12);
            sum1 = _mm_srai_epi32(sum1, 12);
            sum0 = _mm_packs_epi32(sum0, sum1);
            sum0 = _mm_packus_epi16(sum0, sum0);

            *(uint32_t *)qtr = _mm_extract_epi32(sum0, 0); qtr += dst_stride;
            *(uint32_t *)qtr = _mm_extract_epi32(sum0, 1); qtr += dst_stride;

            first_pass_if_dst += 8;
            row_count -= 2;
        } while (row_count != 0);

        pu_width -= 4;
        if (pu_width == 0)
        {
            return;
        }

        first_pass_if_dst += 16;
        dst += 4;
    }

    colCount = pu_width;
    do
    {
        qtr = dst;
        row_count = pu_height;
        a0 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 0 * 8));
        a1 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 1 * 8));
        a2 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 2 * 8));

        do
        {
            a3 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 3 * 8));

            b0 = _mm_unpacklo_epi16(a0, a1);
            b1 = _mm_unpacklo_epi16(a2, a3);
            b2 = _mm_unpackhi_epi16(a0, a1);
            b3 = _mm_unpackhi_epi16(a2, a3);

            sum0 = _mm_madd_epi16(b0, c0);
            sum0 = _mm_add_epi32(sum0, _mm_madd_epi16(b1, c1));
            sum1 = _mm_madd_epi16(b2, c0);
            sum1 = _mm_add_epi32(sum1, _mm_madd_epi16(b3, c1));
            sum0 = _mm_add_epi32(sum0, offset);
            sum1 = _mm_add_epi32(sum1, offset);
            sum0 = _mm_srai_epi32(sum0, 12);
            sum1 = _mm_srai_epi32(sum1, 12);
            sum0 = _mm_packs_epi32(sum0, sum1);
            sum0 = _mm_packus_epi16(sum0, sum0);
            _mm_storel_epi64((__m128i *)qtr, sum0); qtr += dst_stride;

            a0 = a1;
            a1 = a2;
            a2 = a3;
            first_pass_if_dst += 8;
        } while (--row_count != 0);

        first_pass_if_dst += 24;
        colCount -= 8;
        dst += 8;
    } while (colCount != 0);
}


void chroma_interpolation_filter_one_d_horizontal_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    EbByte               dst,
    uint32_t                dst_stride,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst,
    uint32_t                frac_pos_x,
    uint32_t                frac_pos_y)
{
    int32_t row_count, colCount;
    __m128i c0, c2; // coeffs
    __m128i a0, a1, a2, a3;
    __m128i b0, b1;
    __m128i sum;

    EbByte ptr, qtr;

    (void)first_pass_if_dst;
    (void)frac_pos_y;


    ref_pic--;

    PrefetchBlock(ref_pic, src_stride, pu_width + 8, pu_height);

    c0 = _mm_loadl_epi64((__m128i *)chromaFilterCoeff[frac_pos_x]);
    c0 = _mm_packs_epi16(c0, c0);
    c0 = _mm_unpacklo_epi16(c0, c0);
    c2 = _mm_shuffle_epi32(c0, 0x55);
    c0 = _mm_shuffle_epi32(c0, 0x00);

    if (pu_width & 2)
    {
        row_count = pu_height;
        ptr = ref_pic;
        qtr = dst;

        do
        {
            // Need 5 samples, load 8
            a0 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
            a1 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
            a2 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
            a3 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;

            a0 = _mm_unpacklo_epi8(a0, a1); // 10 samples
            a2 = _mm_unpacklo_epi8(a2, a3); // 10 samples

            b0 = _mm_unpacklo_epi16(a0, a2);

            sum = _mm_set1_epi16(32);
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(_mm_shuffle_epi8(b0, _mm_setr_epi8(0, 4, 4, 8, 1, 5, 5, 9, 2, 6, 6, 10, 3, 7, 7, 11)), c0));
            b0 = _mm_unpacklo_epi16(_mm_srli_si128(a0, 4), _mm_srli_si128(a2, 4));
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(_mm_shuffle_epi8(b0, _mm_setr_epi8(0, 4, 4, 8, 1, 5, 5, 9, 2, 6, 6, 10, 3, 7, 7, 11)), c2));
            sum = _mm_srai_epi16(sum, 6);
            sum = _mm_packus_epi16(sum, sum);

            *(uint16_t *)qtr = (uint16_t)(_mm_extract_epi16(sum, 0)); qtr += dst_stride;
            *(uint16_t *)qtr = (uint16_t)(_mm_extract_epi16(sum, 1)); qtr += dst_stride;
            *(uint16_t *)qtr = (uint16_t)(_mm_extract_epi16(sum, 2)); qtr += dst_stride;
            *(uint16_t *)qtr = (uint16_t)(_mm_extract_epi16(sum, 3)); qtr += dst_stride;

            row_count -= 4;
        } while (row_count > 0);

        pu_width -= 2;
        if (pu_width == 0)
        {
            return;
        }

        ref_pic += 2;
        dst += 2;
    }


    if (pu_width & 4)
    {
        row_count = pu_height;
        ptr = ref_pic;
        qtr = dst;
        do
        {
            // Need 7 samples, load 8
            a0 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
            a1 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;

            sum = _mm_set1_epi16(32);
            a0 = _mm_unpacklo_epi64(a0, a1);
            b0 = _mm_shuffle_epi8(a0, _mm_setr_epi8(0, 1, 1, 2, 2, 3, 3, 4, 8, 9, 9, 10, 10, 11, 11, 12));
            b1 = _mm_shuffle_epi8(a0, _mm_setr_epi8(2, 3, 3, 4, 4, 5, 5, 6, 10, 11, 11, 12, 12, 13, 13, 14));
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(b0, c0));
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(b1, c2));

            sum = _mm_srai_epi16(sum, 6);
            sum = _mm_packus_epi16(sum, sum);

            *(uint32_t *)qtr = _mm_extract_epi32(sum, 0); qtr += dst_stride;
            *(uint32_t *)qtr = _mm_extract_epi32(sum, 1); qtr += dst_stride;

            row_count -= 2;
        } while (row_count > 0);

        pu_width -= 4;
        if (pu_width == 0)
        {
            return;
        }

        ref_pic += 4;
        dst += 4;
    }

    colCount = pu_width;
    do
    {
        row_count = pu_height;
        ptr = ref_pic;
        qtr = dst;
        do
        {
            // Need 11 samples, load 16
            a0 = _mm_loadu_si128((__m128i *)ptr); ptr += src_stride;

            sum = _mm_set1_epi16(32);
            a2 = _mm_shuffle_epi8(a0, _mm_setr_epi8(2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10));
            a0 = _mm_shuffle_epi8(a0, _mm_setr_epi8(0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8));
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(a0, c0));
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(a2, c2));
            sum = _mm_srai_epi16(sum, 6);
            sum = _mm_packus_epi16(sum, sum);

            _mm_storel_epi64((__m128i *)qtr, sum); qtr += dst_stride;
        } while (--row_count > 0);

        ref_pic += 8;
        dst += 8;
        colCount -= 8;
    } while (colCount > 0);
}


void chroma_interpolation_filter_one_d_out_raw_horizontal_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    int16_t               *dst,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst,
    uint32_t                frac_pos_x,
    uint32_t                frac_pos_y)
{
    int32_t row_count, colCount;
    __m128i c0, c2; // coeffs
    __m128i a0, a1, a2, a3;
    __m128i b0, b1;
    __m128i sum;

    EbByte ptr;
    (void)first_pass_if_dst;
    (void)frac_pos_y;

    ref_pic--;
    PrefetchBlock(ref_pic, src_stride, pu_width + 8, pu_height);

    c0 = _mm_loadl_epi64((__m128i *)chromaFilterCoeff[frac_pos_x]);
    c0 = _mm_packs_epi16(c0, c0);
    c0 = _mm_unpacklo_epi16(c0, c0);
    c2 = _mm_shuffle_epi32(c0, 0x55);
    c0 = _mm_shuffle_epi32(c0, 0x00);

    if (pu_width & 2)
    {
        row_count = pu_height;
        ptr = ref_pic;
        do
        {
            // Need 5 samples, load 8
            a0 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
            a1 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
            a2 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
            a3 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;

            a0 = _mm_unpacklo_epi8(a0, a1); // 10 samples
            a2 = _mm_unpacklo_epi8(a2, a3); // 10 samples
            b0 = _mm_unpacklo_epi16(a0, a2);

            sum = _mm_maddubs_epi16(_mm_shuffle_epi8(b0, _mm_setr_epi8(0, 4, 4, 8, 1, 5, 5, 9, 2, 6, 6, 10, 3, 7, 7, 11)), c0);
            b0 = _mm_unpacklo_epi16(_mm_srli_si128(a0, 4), _mm_srli_si128(a2, 4));
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(_mm_shuffle_epi8(b0, _mm_setr_epi8(0, 4, 4, 8, 1, 5, 5, 9, 2, 6, 6, 10, 3, 7, 7, 11)), c2));

            _mm_storeu_si128((__m128i *)dst, sum);
            dst += 8;

            row_count -= 4;
        } while (row_count > 0);

        pu_width -= 2;
        if (pu_width == 0)
        {
            return;
        }

        ref_pic += 2;
    }

    if (pu_width & 4)
    {
        row_count = pu_height;
        ptr = ref_pic;
        do
        {

            // Need 7 samples, load 8
            a0 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
            a1 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;

            a0 = _mm_unpacklo_epi64(a0, a1);
            b0 = _mm_shuffle_epi8(a0, _mm_setr_epi8(0, 1, 1, 2, 2, 3, 3, 4, 8, 9, 9, 10, 10, 11, 11, 12));
            b1 = _mm_shuffle_epi8(a0, _mm_setr_epi8(2, 3, 3, 4, 4, 5, 5, 6, 10, 11, 11, 12, 12, 13, 13, 14));
            sum = _mm_maddubs_epi16(b0, c0);
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(b1, c2));

            _mm_storeu_si128((__m128i *)dst, sum);
            dst += 8;

            row_count -= 2;
        } while (row_count > 0);

        pu_width -= 4;
        if (pu_width == 0)
        {
            return;
        }

        ref_pic += 4;
    }


    colCount = pu_width;
    do
    {
        ptr = ref_pic;
        row_count = pu_height;
        do
        {

            // Need 11 samples, load 16
            a0 = _mm_loadu_si128((__m128i *)ptr); ptr += src_stride;
            a2 = _mm_shuffle_epi8(a0, _mm_setr_epi8(2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10));
            a0 = _mm_shuffle_epi8(a0, _mm_setr_epi8(0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8));
            sum = _mm_maddubs_epi16(a0, c0);
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(a2, c2));
            _mm_storeu_si128((__m128i *)dst, sum);
            dst += 8;
        } while (--row_count > 0);

        colCount -= 8;
        ref_pic += 8;
    } while (colCount > 0);
}



void chroma_interpolation_filter_one_d_vertical_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    EbByte               dst,
    uint32_t                dst_stride,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst,
    uint32_t                frac_pos_x,
    uint32_t                frac_pos_y)
{
    int32_t row_count, colCount;
    __m128i c0, c2; // coeffs
    __m128i a0, a1, a2, a3;
    __m128i b0, b1, b2, b3;
    __m128i sum;

    EbByte ptr, qtr;

    (void)first_pass_if_dst;
    (void)frac_pos_x;

    ref_pic -= src_stride;
    PrefetchBlock(ref_pic, src_stride, pu_width, pu_height + 3);


    c0 = _mm_loadl_epi64((__m128i *)chromaFilterCoeff[frac_pos_y]);
    c0 = _mm_packs_epi16(c0, c0);
    c0 = _mm_unpacklo_epi16(c0, c0);
    c2 = _mm_shuffle_epi32(c0, 0x55);
    c0 = _mm_shuffle_epi32(c0, 0x00);

    if (pu_width & 2)
    {
        row_count = pu_height;
        ptr = ref_pic;
        qtr = dst;

        a0 = _mm_setzero_si128();
        a0 = _mm_insert_epi16(a0, *(uint16_t *)ptr, 0); ptr += src_stride;
        a0 = _mm_insert_epi16(a0, *(uint16_t *)ptr, 1); ptr += src_stride;
        a0 = _mm_insert_epi16(a0, *(uint16_t *)ptr, 2); ptr += src_stride;

        do
        {
            a0 = _mm_insert_epi16(a0, *(uint16_t *)ptr, 3); ptr += src_stride;
            a0 = _mm_insert_epi16(a0, *(uint16_t *)ptr, 4); ptr += src_stride;
            a0 = _mm_insert_epi16(a0, *(uint16_t *)ptr, 5); ptr += src_stride;
            a0 = _mm_insert_epi16(a0, *(uint16_t *)ptr, 6); ptr += src_stride;

            sum = _mm_set1_epi16(32);
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(_mm_shuffle_epi8(a0, _mm_setr_epi8(0, 2, 1, 3, 2, 4, 3, 5, 4, 6, 5, 7, 6, 8, 7, 9)), c0));
            a0 = _mm_srli_si128(a0, 4);
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(_mm_shuffle_epi8(a0, _mm_setr_epi8(0, 2, 1, 3, 2, 4, 3, 5, 4, 6, 5, 7, 6, 8, 7, 9)), c2));
            a0 = _mm_srli_si128(a0, 4);
            sum = _mm_srai_epi16(sum, 6);
            sum = _mm_packus_epi16(sum, sum);

            *(uint16_t *)qtr = (uint16_t)(_mm_extract_epi16(sum, 0)); qtr += dst_stride;
            *(uint16_t *)qtr = (uint16_t)(_mm_extract_epi16(sum, 1)); qtr += dst_stride;
            *(uint16_t *)qtr = (uint16_t)(_mm_extract_epi16(sum, 2)); qtr += dst_stride;
            *(uint16_t *)qtr = (uint16_t)(_mm_extract_epi16(sum, 3)); qtr += dst_stride;

            row_count -= 4;
        } while (row_count > 0);

        pu_width -= 2;
        if (pu_width == 0)
        {
            return;
        }

        ref_pic += 2;
        dst += 2;
    }

    if (pu_width & 4)
    {
        row_count = pu_height;
        ptr = ref_pic;
        qtr = dst;
        b0 = _mm_cvtsi32_si128(*(uint32_t *)ptr); ptr += src_stride;
        b1 = _mm_cvtsi32_si128(*(uint32_t *)ptr); ptr += src_stride;
        b2 = _mm_cvtsi32_si128(*(uint32_t *)ptr); ptr += src_stride;
        a0 = _mm_unpacklo_epi32(b2, b1);
        a0 = _mm_unpacklo_epi64(a0, b0);

        do
        {
            sum = _mm_set1_epi16(32);
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(_mm_shuffle_epi8(a0, _mm_setr_epi8(8, 4, 9, 5, 10, 6, 11, 7, 4, 0, 5, 1, 6, 2, 7, 3)), c0));
            b3 = _mm_cvtsi32_si128(*(uint32_t *)ptr); ptr += src_stride;
            b0 = _mm_cvtsi32_si128(*(uint32_t *)ptr); ptr += src_stride;
            b3 = _mm_unpacklo_epi32(b0, b3);
            a0 = _mm_unpacklo_epi64(b3, a0);
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(_mm_shuffle_epi8(a0, _mm_setr_epi8(8, 4, 9, 5, 10, 6, 11, 7, 4, 0, 5, 1, 6, 2, 7, 3)), c2));
            sum = _mm_srai_epi16(sum, 6);
            sum = _mm_packus_epi16(sum, sum);

            *(uint32_t *)qtr = _mm_extract_epi32(sum, 0); qtr += dst_stride;
            *(uint32_t *)qtr = _mm_extract_epi32(sum, 1); qtr += dst_stride;

            row_count -= 2;
        } while (row_count > 0);

        pu_width -= 4;
        if (pu_width == 0)
        {
            return;
        }

        ref_pic += 4;
        dst += 4;
    }

    colCount = pu_width;
    do
    {
        row_count = pu_height;
        ptr = ref_pic;
        qtr = dst;
        a0 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
        a1 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
        a2 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;

        do
        {

            sum = _mm_set1_epi16(32);
            a3 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(_mm_unpacklo_epi8(a0, a1), c0));
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(_mm_unpacklo_epi8(a2, a3), c2));

            a0 = a1;
            a1 = a2;
            a2 = a3;

            sum = _mm_srai_epi16(sum, 6);
            sum = _mm_packus_epi16(sum, sum);

            _mm_storel_epi64((__m128i *)qtr, sum); qtr += dst_stride;

        } while (--row_count > 0);

        ref_pic += 8;
        dst += 8;
        colCount -= 8;
    } while (colCount > 0);
}


void chroma_interpolation_filter_one_d_out_raw_vertical_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    int16_t               *dst,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst,
    uint32_t                frac_pos_x,
    uint32_t                frac_pos_y)
{
    int32_t row_count, colCount;
    __m128i c0, c2; // coeffs
    __m128i a0, a1, a2, a3;
    __m128i b0, b1, b2, b3;
    __m128i sum;

    EbByte ptr;

    (void)first_pass_if_dst;
    (void)frac_pos_x;

    ref_pic -= src_stride;
    PrefetchBlock(ref_pic, src_stride, pu_width, pu_height + 3);

    c0 = _mm_loadl_epi64((__m128i *)chromaFilterCoeff[frac_pos_y]);
    c0 = _mm_packs_epi16(c0, c0);
    c0 = _mm_unpacklo_epi16(c0, c0);
    c2 = _mm_shuffle_epi32(c0, 0x55);
    c0 = _mm_shuffle_epi32(c0, 0x00);

    if (pu_width & 2)
    {
        row_count = pu_height;
        ptr = ref_pic;
        a0 = _mm_setzero_si128();
        a0 = _mm_insert_epi16(a0, *(uint16_t *)ptr, 0); ptr += src_stride;
        a0 = _mm_insert_epi16(a0, *(uint16_t *)ptr, 1); ptr += src_stride;
        a0 = _mm_insert_epi16(a0, *(uint16_t *)ptr, 2); ptr += src_stride;

        do
        {
            a0 = _mm_insert_epi16(a0, *(uint16_t *)ptr, 3); ptr += src_stride;
            a0 = _mm_insert_epi16(a0, *(uint16_t *)ptr, 4); ptr += src_stride;
            a0 = _mm_insert_epi16(a0, *(uint16_t *)ptr, 5); ptr += src_stride;
            a0 = _mm_insert_epi16(a0, *(uint16_t *)ptr, 6); ptr += src_stride;

            sum = _mm_maddubs_epi16(_mm_shuffle_epi8(a0, _mm_setr_epi8(0, 2, 1, 3, 2, 4, 3, 5, 4, 6, 5, 7, 6, 8, 7, 9)), c0);
            a0 = _mm_srli_si128(a0, 4);
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(_mm_shuffle_epi8(a0, _mm_setr_epi8(0, 2, 1, 3, 2, 4, 3, 5, 4, 6, 5, 7, 6, 8, 7, 9)), c2));
            a0 = _mm_srli_si128(a0, 4);

            _mm_storeu_si128((__m128i *)dst, sum);
            dst += 8;

            row_count -= 4;
        } while (row_count > 0);

        pu_width -= 2;
        if (pu_width == 0)
        {
            return;
        }

        ref_pic += 2;
    }

    if (pu_width & 4)
    {
        row_count = pu_height;
        ptr = ref_pic;
        b0 = _mm_cvtsi32_si128(*(uint32_t *)ptr); ptr += src_stride;
        b1 = _mm_cvtsi32_si128(*(uint32_t *)ptr); ptr += src_stride;
        b2 = _mm_cvtsi32_si128(*(uint32_t *)ptr); ptr += src_stride;
        a0 = _mm_unpacklo_epi32(b2, b1);
        a0 = _mm_unpacklo_epi64(a0, b0);

        do
        {
            sum = _mm_maddubs_epi16(_mm_shuffle_epi8(a0, _mm_setr_epi8(8, 4, 9, 5, 10, 6, 11, 7, 4, 0, 5, 1, 6, 2, 7, 3)), c0);
            b3 = _mm_cvtsi32_si128(*(uint32_t *)ptr); ptr += src_stride;
            b0 = _mm_cvtsi32_si128(*(uint32_t *)ptr); ptr += src_stride;
            b3 = _mm_unpacklo_epi32(b0, b3);
            a0 = _mm_unpacklo_epi64(b3, a0);
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(_mm_shuffle_epi8(a0, _mm_setr_epi8(8, 4, 9, 5, 10, 6, 11, 7, 4, 0, 5, 1, 6, 2, 7, 3)), c2));

            _mm_storeu_si128((__m128i *)dst, sum);
            dst += 8;

            row_count -= 2;
        } while (row_count > 0);

        pu_width -= 4;
        if (pu_width == 0)
        {
            return;
        }

        ref_pic += 4;
    }

    colCount = pu_width;
    do
    {
        row_count = pu_height;
        ptr = ref_pic;
        a0 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
        a1 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
        a2 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;

        do
        {

            a3 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
            sum = _mm_maddubs_epi16(_mm_unpacklo_epi8(a0, a1), c0);
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(_mm_unpacklo_epi8(a2, a3), c2));
            a0 = a1;
            a1 = a2;
            a2 = a3;

            _mm_storeu_si128((__m128i *)dst, sum);
            dst += 8;
        } while (--row_count > 0);

        ref_pic += 8;
        colCount -= 8;
    } while (colCount > 0);
}



void chroma_interpolation_filter_two_d_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    EbByte               dst,
    uint32_t                dst_stride,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst,
    uint32_t                frac_pos_x,
    uint32_t                frac_pos_y)
{
    chroma_interpolation_filter_one_d_out_raw_horizontal_ssse3(ref_pic - src_stride, src_stride, first_pass_if_dst, pu_width, pu_height + 3, NULL, frac_pos_x, 0);
    ChromaInterpolationFilterTwoDInRaw_SSSE3(first_pass_if_dst, dst, dst_stride, pu_width, pu_height, frac_pos_y);
}

void ChromaInterpolationFilterTwoDInRawOutRaw_SSSE3(int16_t *first_pass_if_dst, int16_t *dst, uint32_t pu_width, uint32_t pu_height, uint32_t frac_pos_y)
{
    int32_t row_count, colCount;
    __m128i c0, c1; // coeffs

    c0 = _mm_loadl_epi64((__m128i *)chromaFilterCoeff[frac_pos_y]);
    c1 = _mm_shuffle_epi32(c0, 0x55);
    c0 = _mm_shuffle_epi32(c0, 0);

    if (pu_width & 2)
    {
        row_count = pu_height;
        do
        {
            __m128i a0 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 0));
            __m128i a1 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 2));
            __m128i a2 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 4));
            __m128i a3 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 6));
            __m128i b0, b1, b2, b3;
            __m128i sum0, sum1;

            b0 = _mm_unpacklo_epi16(a0, a1);
            b1 = _mm_unpacklo_epi16(a2, a3);
            b2 = _mm_unpackhi_epi16(a0, a1);
            b3 = _mm_unpackhi_epi16(a2, a3);
            sum0 = _mm_madd_epi16(b0, c0);
            sum0 = _mm_add_epi32(sum0, _mm_madd_epi16(b1, c1));
            sum1 = _mm_madd_epi16(b2, c0);
            sum1 = _mm_add_epi32(sum1, _mm_madd_epi16(b3, c1));
            sum0 = _mm_srai_epi32(sum0, 6);
            sum1 = _mm_srai_epi32(sum1, 6);
            sum0 = _mm_packs_epi32(sum0, sum1);
            _mm_storeu_si128((__m128i *)dst, sum0);

            first_pass_if_dst += 8;
            dst += 8;
            row_count -= 4;
        } while (row_count > 0);

        pu_width -= 2;
        if (pu_width == 0)
        {
            return;
        }
        first_pass_if_dst += 8;
    }

    if (pu_width & 4)
    {
        row_count = pu_height;
        do
        {
            __m128i a0 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 0));
            __m128i a1 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 4));
            __m128i a2 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 8));
            __m128i a3 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 12));
            __m128i b0, b1, b2, b3;
            __m128i sum0, sum1;

            b0 = _mm_unpacklo_epi16(a0, a1);
            b1 = _mm_unpacklo_epi16(a2, a3);
            b2 = _mm_unpackhi_epi16(a0, a1);
            b3 = _mm_unpackhi_epi16(a2, a3);
            sum0 = _mm_madd_epi16(b0, c0);
            sum0 = _mm_add_epi32(sum0, _mm_madd_epi16(b1, c1));
            sum1 = _mm_madd_epi16(b2, c0);
            sum1 = _mm_add_epi32(sum1, _mm_madd_epi16(b3, c1));
            sum0 = _mm_srai_epi32(sum0, 6);
            sum1 = _mm_srai_epi32(sum1, 6);
            sum0 = _mm_packs_epi32(sum0, sum1);
            _mm_storeu_si128((__m128i *)dst, sum0);

            first_pass_if_dst += 8;
            dst += 8;
            row_count -= 2;
        } while (row_count > 0);

        pu_width -= 4;
        if (pu_width == 0)
        {
            return;
        }

        first_pass_if_dst += 16;
    }

    colCount = pu_width;
    do
    {
        __m128i a0 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 0 * 8));
        __m128i a1 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 1 * 8));
        __m128i a2 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 2 * 8));
        row_count = pu_height;

        do
        {
            __m128i a3 = _mm_loadu_si128((__m128i *)(first_pass_if_dst + 3 * 8));
            __m128i b0, b1, b2, b3;
            __m128i sum0, sum1;
            b0 = _mm_unpacklo_epi16(a0, a1);
            b1 = _mm_unpacklo_epi16(a2, a3);
            b2 = _mm_unpackhi_epi16(a0, a1);
            b3 = _mm_unpackhi_epi16(a2, a3);
            sum0 = _mm_madd_epi16(b0, c0);
            sum0 = _mm_add_epi32(sum0, _mm_madd_epi16(b1, c1));
            sum1 = _mm_madd_epi16(b2, c0);
            sum1 = _mm_add_epi32(sum1, _mm_madd_epi16(b3, c1));
            sum0 = _mm_srai_epi32(sum0, 6);
            sum1 = _mm_srai_epi32(sum1, 6);
            sum0 = _mm_packs_epi32(sum0, sum1);
            _mm_storeu_si128((__m128i *)dst, sum0);

            a0 = a1;
            a1 = a2;
            a2 = a3;
            first_pass_if_dst += 8;
            dst += 8;
        } while (--row_count > 0);

        first_pass_if_dst += 24;
        colCount -= 8;
    } while (colCount > 0);
}


void chroma_interpolation_filter_two_d_out_raw_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    int16_t               *dst,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst,
    uint32_t                frac_pos_x,
    uint32_t                frac_pos_y)
{
    chroma_interpolation_filter_one_d_out_raw_horizontal_ssse3(ref_pic - src_stride, src_stride, first_pass_if_dst, pu_width, pu_height + 3, NULL, frac_pos_x, 0);
    ChromaInterpolationFilterTwoDInRawOutRaw_SSSE3(first_pass_if_dst, dst, pu_width, pu_height, frac_pos_y);
}



void LumaInterpolationFilterOneDHorizontal_SSSE3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    EbByte               dst,
    uint32_t                dst_stride,
    uint32_t                pu_width,
    uint32_t                pu_height,
    uint32_t                frac_pos_x)
{
    int32_t row_count, colCount;
    __m128i c0, c1, c2, c3; // coeffs
    __m128i a0, a1;
    __m128i b0;
    __m128i sum;


    ref_pic -= 3;

    PrefetchBlock(ref_pic, src_stride, (pu_width == 4) ? 16 : pu_width + 8, (pu_width == 4) ? ((pu_height + 1)&~1) : pu_height);

    c0 = _mm_loadu_si128((__m128i *)lumaFilterCoeff[frac_pos_x]);
    c0 = _mm_packs_epi16(c0, c0);
    c0 = _mm_unpacklo_epi16(c0, c0);
    c3 = _mm_shuffle_epi32(c0, 0xff);
    c2 = _mm_shuffle_epi32(c0, 0xaa);
    c1 = _mm_shuffle_epi32(c0, 0x55);
    c0 = _mm_shuffle_epi32(c0, 0x00);

    if (pu_width & 4)
    {
        EbByte ptr = ref_pic;
        EbByte qtr = dst;
        row_count = pu_height;
        do
        {
            a0 = _mm_loadu_si128((__m128i *)ptr); ptr += src_stride;
            a1 = _mm_loadu_si128((__m128i *)ptr); ptr += src_stride;


            sum = _mm_set1_epi16(32);

            b0 = _mm_unpacklo_epi64(a0, a1);
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(_mm_shuffle_epi8(b0, _mm_setr_epi8(0, 1, 1, 2, 2, 3, 3, 4, 8, 9, 9, 10, 10, 11, 11, 12)), c0));
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(_mm_shuffle_epi8(b0, _mm_setr_epi8(2, 3, 3, 4, 4, 5, 5, 6, 10, 11, 11, 12, 12, 13, 13, 14)), c1));
            b0 = _mm_unpacklo_epi64(_mm_srli_si128(a0, 4), _mm_srli_si128(a1, 4));
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(_mm_shuffle_epi8(b0, _mm_setr_epi8(0, 1, 1, 2, 2, 3, 3, 4, 8, 9, 9, 10, 10, 11, 11, 12)), c2));
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(_mm_shuffle_epi8(b0, _mm_setr_epi8(2, 3, 3, 4, 4, 5, 5, 6, 10, 11, 11, 12, 12, 13, 13, 14)), c3));
            sum = _mm_srai_epi16(sum, 6);
            sum = _mm_packus_epi16(sum, sum);

            *(uint32_t *)qtr = _mm_extract_epi32(sum, 0); qtr += dst_stride;
            *(uint32_t *)qtr = _mm_extract_epi32(sum, 1); qtr += dst_stride;

            row_count -= 2;
        } while (row_count > 0);

        pu_width -= 4;
        if (pu_width == 0)
        {
            return;
        }

        ref_pic += 4;
        dst += 4;
    }

    colCount = pu_width;
    do
    {
        EbByte ptr = ref_pic;
        EbByte qtr = dst;

        row_count = pu_height;
        do
        {
            a0 = _mm_loadu_si128((__m128i *)ptr); ptr += src_stride;

            sum = _mm_set1_epi16(32);
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(_mm_shuffle_epi8(a0, _mm_setr_epi8(0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8)), c0));
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(_mm_shuffle_epi8(a0, _mm_setr_epi8(2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10)), c1));
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(_mm_shuffle_epi8(a0, _mm_setr_epi8(4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12)), c2));
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(_mm_shuffle_epi8(a0, _mm_setr_epi8(6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14)), c3));
            sum = _mm_srai_epi16(sum, 6);
            sum = _mm_packus_epi16(sum, sum);

            _mm_storel_epi64((__m128i *)qtr, sum); qtr += dst_stride;
        } while (--row_count > 0);

        ref_pic += 8;
        dst += 8;

        colCount -= 8;
    } while (colCount > 0);
}

void LumaInterpolationFilterOneDOutRawHorizontal_SSSE3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    int16_t               *dst,
    uint32_t                pu_width,
    uint32_t                pu_height,
    uint32_t                frac_pos_x)
{
    int32_t row_count, colCount;
    __m128i c0, c1, c2, c3; // coeffs
    __m128i a0, a1;
    __m128i b0;
    __m128i sum;
    EbByte ptr;

    ref_pic -= 3;

    PrefetchBlock(ref_pic, src_stride, (pu_width == 4) ? 16 : pu_width + 8, (pu_width == 4) ? ((pu_height + 1)&~1) : pu_height);


    c0 = _mm_loadu_si128((__m128i *)lumaFilterCoeff[frac_pos_x]);
    c0 = _mm_packs_epi16(c0, c0);
    c0 = _mm_unpacklo_epi16(c0, c0);
    c3 = _mm_shuffle_epi32(c0, 0xff);
    c2 = _mm_shuffle_epi32(c0, 0xaa);
    c1 = _mm_shuffle_epi32(c0, 0x55);
    c0 = _mm_shuffle_epi32(c0, 0x00);

    if (pu_width & 4)
    {
        ptr = ref_pic;
        row_count = pu_height;
        do
        {
            a0 = _mm_loadu_si128((__m128i *)ptr); ptr += src_stride;
            a1 = _mm_loadu_si128((__m128i *)ptr); ptr += src_stride;


            b0 = _mm_unpacklo_epi64(a0, a1);
            sum = _mm_maddubs_epi16(_mm_shuffle_epi8(b0, _mm_setr_epi8(0, 1, 1, 2, 2, 3, 3, 4, 8, 9, 9, 10, 10, 11, 11, 12)), c0);
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(_mm_shuffle_epi8(b0, _mm_setr_epi8(2, 3, 3, 4, 4, 5, 5, 6, 10, 11, 11, 12, 12, 13, 13, 14)), c1));
            b0 = _mm_unpacklo_epi64(_mm_srli_si128(a0, 4), _mm_srli_si128(a1, 4));
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(_mm_shuffle_epi8(b0, _mm_setr_epi8(0, 1, 1, 2, 2, 3, 3, 4, 8, 9, 9, 10, 10, 11, 11, 12)), c2));
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(_mm_shuffle_epi8(b0, _mm_setr_epi8(2, 3, 3, 4, 4, 5, 5, 6, 10, 11, 11, 12, 12, 13, 13, 14)), c3));

            sum = _mm_sub_epi16(sum, _mm_set1_epi16(128 * 64));

            _mm_storeu_si128((__m128i *)dst, sum);
            dst += 8;

            row_count -= 2;
        } while (row_count > 0);

        pu_width -= 4;
        if (pu_width == 0)
        {
            return;
        }

        ref_pic += 4;
    }

    colCount = pu_width;
    do
    {
        ptr = ref_pic;
        row_count = pu_height;
        do
        {
            a0 = _mm_loadu_si128((__m128i *)ptr); ptr += src_stride;

            sum = _mm_maddubs_epi16(_mm_shuffle_epi8(a0, _mm_setr_epi8(0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8)), c0);
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(_mm_shuffle_epi8(a0, _mm_setr_epi8(2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10)), c1));
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(_mm_shuffle_epi8(a0, _mm_setr_epi8(4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12)), c2));
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(_mm_shuffle_epi8(a0, _mm_setr_epi8(6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14)), c3));

            sum = _mm_sub_epi16(sum, _mm_set1_epi16(128 * 64));

            _mm_storeu_si128((__m128i *)dst, sum);
            dst += 8;
        } while (--row_count > 0);

        ref_pic += 8;
        colCount -= 8;
    } while (colCount > 0);
}

void LumaInterpolationFilterOneDVertical_SSSE3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    EbByte               dst,
    uint32_t                dst_stride,
    uint32_t                pu_width,
    uint32_t                pu_height,
    uint32_t                frac_pos_x)
{
    int32_t row_count, colCount;
    __m128i c0, c1, c2, c3; // coeffs
    __m128i a0, a1, a2, a3, a4, a5, a6, a7, a8;
    __m128i b0, b1, b2, b3, b4, b5, b6, b7;
    __m128i sum;

    ref_pic -= 3 * src_stride;

    PrefetchBlock(ref_pic, src_stride, pu_width, pu_height + 7);

    c0 = _mm_loadu_si128((__m128i *)lumaFilterCoeff[frac_pos_x]);
    c0 = _mm_packs_epi16(c0, c0); // Convert 16-bit coefficients to 8 bits
    c0 = _mm_unpacklo_epi16(c0, c0);
    c3 = _mm_shuffle_epi32(c0, 0xff);
    c2 = _mm_shuffle_epi32(c0, 0xaa);
    c1 = _mm_shuffle_epi32(c0, 0x55);
    c0 = _mm_shuffle_epi32(c0, 0x00);

    if (pu_width & 4)
    {
        EbByte ptr = ref_pic;
        EbByte qtr = dst;

        row_count = pu_height;
        a0 = _mm_cvtsi32_si128(*(uint32_t *)ptr); ptr += src_stride;
        a1 = _mm_cvtsi32_si128(*(uint32_t *)ptr); ptr += src_stride;
        a2 = _mm_cvtsi32_si128(*(uint32_t *)ptr); ptr += src_stride;
        a3 = _mm_cvtsi32_si128(*(uint32_t *)ptr); ptr += src_stride;
        a4 = _mm_cvtsi32_si128(*(uint32_t *)ptr); ptr += src_stride;
        a5 = _mm_cvtsi32_si128(*(uint32_t *)ptr); ptr += src_stride;
        a6 = _mm_cvtsi32_si128(*(uint32_t *)ptr); ptr += src_stride;

        b0 = _mm_unpacklo_epi32(_mm_setzero_si128(), a0);
        b1 = _mm_unpacklo_epi32(a1, a2);
        b2 = _mm_unpacklo_epi32(a3, a4);
        b3 = _mm_unpacklo_epi32(a5, a6);
        b0 = _mm_unpacklo_epi64(b0, b1);
        b1 = _mm_unpacklo_epi64(b2, b3);

        do
        {
            sum = _mm_set1_epi16(32);

            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(_mm_shuffle_epi8(b0, _mm_setr_epi8(4, 8, 5, 9, 6, 10, 7, 11, 8, 12, 9, 13, 10, 14, 11, 15)), c0));
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(_mm_shuffle_epi8(b1, _mm_setr_epi8(4, 8, 5, 9, 6, 10, 7, 11, 8, 12, 9, 13, 10, 14, 11, 15)), c2));

            a7 = _mm_cvtsi32_si128(*(uint32_t *)ptr); ptr += src_stride;
            a8 = _mm_cvtsi32_si128(*(uint32_t *)ptr); ptr += src_stride;
            b2 = _mm_unpacklo_epi32(a7, a8);
            b0 = _mm_alignr_epi8(b1, b0, 8);
            b1 = _mm_alignr_epi8(b2, b1, 8);
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(_mm_shuffle_epi8(b0, _mm_setr_epi8(4, 8, 5, 9, 6, 10, 7, 11, 8, 12, 9, 13, 10, 14, 11, 15)), c1));
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(_mm_shuffle_epi8(b1, _mm_setr_epi8(4, 8, 5, 9, 6, 10, 7, 11, 8, 12, 9, 13, 10, 14, 11, 15)), c3));
            sum = _mm_srai_epi16(sum, 6);
            sum = _mm_packus_epi16(sum, sum);

            *(uint32_t *)qtr = _mm_extract_epi32(sum, 0); qtr += dst_stride;
            *(uint32_t *)qtr = _mm_extract_epi32(sum, 1); qtr += dst_stride;

            row_count -= 2;
        } while (row_count > 0);

        pu_width -= 4;
        if (pu_width == 0)
        {
            return;
        }

        ref_pic += 4;
        dst += 4;
    }

    colCount = pu_width;
    do
    {
        EbByte ptr = ref_pic;
        EbByte qtr = dst;

        row_count = pu_height;
        a0 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
        a1 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
        a2 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
        a3 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
        a4 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
        a5 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
        a6 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
        b0 = _mm_unpacklo_epi8(a0, a1);
        b1 = _mm_unpacklo_epi8(a1, a2);
        b2 = _mm_unpacklo_epi8(a2, a3);
        b3 = _mm_unpacklo_epi8(a3, a4);
        b4 = _mm_unpacklo_epi8(a4, a5);
        b5 = _mm_unpacklo_epi8(a5, a6);
        do
        {
            sum = _mm_set1_epi16(32);
            a7 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
            a8 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
            b6 = _mm_unpacklo_epi8(a6, a7);
            b7 = _mm_unpacklo_epi8(a7, a8);

            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(b0, c0));
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(b2, c1));
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(b4, c2));
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(b6, c3));

            sum = _mm_srai_epi16(sum, 6);
            sum = _mm_packus_epi16(sum, sum);

            _mm_storel_epi64((__m128i *)qtr, sum); qtr += dst_stride;

            sum = _mm_set1_epi16(32);
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(b1, c0));
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(b3, c1));
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(b5, c2));
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(b7, c3));
            sum = _mm_srai_epi16(sum, 6);
            sum = _mm_packus_epi16(sum, sum);

            _mm_storel_epi64((__m128i *)qtr, sum); qtr += dst_stride;

            b0 = b2;
            b1 = b3;
            b2 = b4;
            b3 = b5;
            b4 = b6;
            b5 = b7;
            a6 = a8;
            row_count -= 2;
        } while (row_count > 0);

        ref_pic += 8;
        dst += 8;

        colCount -= 8;
    } while (colCount > 0);
}

void LumaInterpolationFilterOneDOutRawVertical_SSSE3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    int16_t*               dst,
    uint32_t                pu_width,
    uint32_t                pu_height,
    uint32_t                frac_pos_x)
{
    int32_t row_count, colCount;
    __m128i c0, c1, c2, c3; // coeffs
    __m128i a0, a1, a2, a3, a4, a5, a6, a7, a8;
    __m128i b0, b1, b2, b3, b4, b5, b6, b7;
    __m128i sum;
    EbByte ptr;

    ref_pic -= 3 * src_stride;

    PrefetchBlock(ref_pic, src_stride, pu_width, pu_height + 7);

    c0 = _mm_loadu_si128((__m128i *)lumaFilterCoeff[frac_pos_x]);
    c0 = _mm_packs_epi16(c0, c0); // Convert 16-bit coefficients to 8 bits
    c0 = _mm_unpacklo_epi16(c0, c0);
    c3 = _mm_shuffle_epi32(c0, 0xff);
    c2 = _mm_shuffle_epi32(c0, 0xaa);
    c1 = _mm_shuffle_epi32(c0, 0x55);
    c0 = _mm_shuffle_epi32(c0, 0x00);

    if (pu_width & 4)
    {
        row_count = pu_height;
        ptr = ref_pic;
        a0 = _mm_cvtsi32_si128(*(uint32_t *)ptr); ptr += src_stride;
        a1 = _mm_cvtsi32_si128(*(uint32_t *)ptr); ptr += src_stride;
        a2 = _mm_cvtsi32_si128(*(uint32_t *)ptr); ptr += src_stride;
        a3 = _mm_cvtsi32_si128(*(uint32_t *)ptr); ptr += src_stride;
        a4 = _mm_cvtsi32_si128(*(uint32_t *)ptr); ptr += src_stride;
        a5 = _mm_cvtsi32_si128(*(uint32_t *)ptr); ptr += src_stride;
        a6 = _mm_cvtsi32_si128(*(uint32_t *)ptr); ptr += src_stride;

        b0 = _mm_unpacklo_epi32(_mm_setzero_si128(), a0);
        b1 = _mm_unpacklo_epi32(a1, a2);
        b2 = _mm_unpacklo_epi32(a3, a4);
        b3 = _mm_unpacklo_epi32(a5, a6);
        b0 = _mm_unpacklo_epi64(b0, b1);
        b1 = _mm_unpacklo_epi64(b2, b3);

        do
        {

            sum = _mm_maddubs_epi16(_mm_shuffle_epi8(b0, _mm_setr_epi8(4, 8, 5, 9, 6, 10, 7, 11, 8, 12, 9, 13, 10, 14, 11, 15)), c0);
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(_mm_shuffle_epi8(b1, _mm_setr_epi8(4, 8, 5, 9, 6, 10, 7, 11, 8, 12, 9, 13, 10, 14, 11, 15)), c2));

            a7 = _mm_cvtsi32_si128(*(uint32_t *)ptr); ptr += src_stride;
            a8 = _mm_cvtsi32_si128(*(uint32_t *)ptr); ptr += src_stride;
            b2 = _mm_unpacklo_epi32(a7, a8);
            b0 = _mm_alignr_epi8(b1, b0, 8);
            b1 = _mm_alignr_epi8(b2, b1, 8);

            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(_mm_shuffle_epi8(b0, _mm_setr_epi8(4, 8, 5, 9, 6, 10, 7, 11, 8, 12, 9, 13, 10, 14, 11, 15)), c1));
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(_mm_shuffle_epi8(b1, _mm_setr_epi8(4, 8, 5, 9, 6, 10, 7, 11, 8, 12, 9, 13, 10, 14, 11, 15)), c3));
            sum = _mm_sub_epi16(sum, _mm_set1_epi16(128 * 64));

            _mm_storeu_si128((__m128i *)dst, sum);
            dst += 8;

            row_count -= 2;
        } while (row_count > 0);

        pu_width -= 4;
        if (pu_width == 0)
        {
            return;
        }

        ref_pic += 4;
    }

    colCount = pu_width;
    do
    {
        ptr = ref_pic;
        row_count = pu_height;
        a0 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
        a1 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
        a2 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
        a3 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
        a4 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
        a5 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
        a6 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
        b0 = _mm_unpacklo_epi8(a0, a1);
        b1 = _mm_unpacklo_epi8(a1, a2);
        b2 = _mm_unpacklo_epi8(a2, a3);
        b3 = _mm_unpacklo_epi8(a3, a4);
        b4 = _mm_unpacklo_epi8(a4, a5);
        b5 = _mm_unpacklo_epi8(a5, a6);
        do
        {
            a7 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
            a8 = _mm_loadl_epi64((__m128i *)ptr); ptr += src_stride;
            b6 = _mm_unpacklo_epi8(a6, a7);
            b7 = _mm_unpacklo_epi8(a7, a8);

            sum = _mm_maddubs_epi16(b0, c0);
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(b2, c1));
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(b4, c2));
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(b6, c3));

            sum = _mm_sub_epi16(sum, _mm_set1_epi16(128 * 64));
            _mm_storeu_si128((__m128i *)dst, sum);
            dst += 8;

            sum = _mm_maddubs_epi16(b1, c0);
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(b3, c1));
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(b5, c2));
            sum = _mm_add_epi16(sum, _mm_maddubs_epi16(b7, c3));

            sum = _mm_sub_epi16(sum, _mm_set1_epi16(128 * 64));

            _mm_storeu_si128((__m128i *)dst, sum);
            dst += 8;

            b0 = b2;
            b1 = b3;
            b2 = b4;
            b3 = b5;
            b4 = b6;
            b5 = b7;
            a6 = a8;
            row_count -= 2;
        } while (row_count > 0);

        ref_pic += 8;
        colCount -= 8;
    } while (colCount > 0);
}


void luma_interpolation_filter_posa_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    EbByte               dst,
    uint32_t                dst_stride,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst)
{
    (void)first_pass_if_dst;

    LumaInterpolationFilterOneDHorizontal_SSSE3(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, 1);
}

void luma_interpolation_filter_posb_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    EbByte               dst,
    uint32_t                dst_stride,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst)
{
    (void)first_pass_if_dst;

    LumaInterpolationFilterOneDHorizontal_SSSE3(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, 2);
}

void luma_interpolation_filter_posc_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    EbByte               dst,
    uint32_t                dst_stride,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst)
{
    (void)first_pass_if_dst;

    LumaInterpolationFilterOneDHorizontal_SSSE3(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, 3);
}

void luma_interpolation_filter_posd_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    EbByte               dst,
    uint32_t                dst_stride,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst)
{
    (void)first_pass_if_dst;

    LumaInterpolationFilterOneDVertical_SSSE3(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, 1);
}

void luma_interpolation_filter_posh_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    EbByte               dst,
    uint32_t                dst_stride,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst)
{
    (void)first_pass_if_dst;

    LumaInterpolationFilterOneDVertical_SSSE3(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, 2);
}

void luma_interpolation_filter_posn_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    EbByte               dst,
    uint32_t                dst_stride,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst)
{
    (void)first_pass_if_dst;

    LumaInterpolationFilterOneDVertical_SSSE3(ref_pic, src_stride, dst, dst_stride, pu_width, pu_height, 3);
}

void luma_interpolation_filter_pose_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    EbByte               dst,
    uint32_t                dst_stride,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst)
{
    LumaInterpolationFilterOneDOutRawHorizontal_SSSE3(ref_pic - 3 * src_stride, src_stride, first_pass_if_dst, pu_width, pu_height + 6, 1);
    LumaInterpolationFilterTwoDInRaw7_SSSE3(first_pass_if_dst, dst, dst_stride, pu_width, pu_height, 1);
}


void luma_interpolation_filter_posf_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    EbByte               dst,
    uint32_t                dst_stride,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst)
{
    uint32_t puHeight1 = pu_height + 6;
    EbByte refPic1 = ref_pic - 3 * src_stride;
    LumaInterpolationFilterOneDOutRawHorizontal_SSSE3(refPic1, src_stride, first_pass_if_dst, pu_width, puHeight1, 2);
    LumaInterpolationFilterTwoDInRaw7_SSSE3(first_pass_if_dst, dst, dst_stride, pu_width, pu_height, 1);
}

void luma_interpolation_filter_posg_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    EbByte               dst,
    uint32_t                dst_stride,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst)
{
    LumaInterpolationFilterOneDOutRawHorizontal_SSSE3(ref_pic - 3 * src_stride, src_stride, first_pass_if_dst, pu_width, pu_height + 6, 3);
    LumaInterpolationFilterTwoDInRaw7_SSSE3(first_pass_if_dst, dst, dst_stride, pu_width, pu_height, 1);
}

void luma_interpolation_filter_posi_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    EbByte               dst,
    uint32_t                dst_stride,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst)
{
    LumaInterpolationFilterOneDOutRawHorizontal_SSSE3(ref_pic - 3 * src_stride, src_stride, first_pass_if_dst, pu_width, pu_height + 7, 1);
    LumaInterpolationFilterTwoDInRawM_SSSE3(first_pass_if_dst, dst, dst_stride, pu_width, pu_height);
}



void luma_interpolation_filter_posj_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    EbByte               dst,
    uint32_t                dst_stride,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst)
{
    LumaInterpolationFilterOneDOutRawHorizontal_SSSE3(ref_pic - 3 * src_stride, src_stride, first_pass_if_dst, pu_width, pu_height + 7, 2);
    LumaInterpolationFilterTwoDInRawM_SSSE3(first_pass_if_dst, dst, dst_stride, pu_width, pu_height);
}

void luma_interpolation_filter_posk_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    EbByte               dst,
    uint32_t                dst_stride,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst)
{
    LumaInterpolationFilterOneDOutRawHorizontal_SSSE3(ref_pic - 3 * src_stride, src_stride, first_pass_if_dst, pu_width, pu_height + 7, 3);
    LumaInterpolationFilterTwoDInRawM_SSSE3(first_pass_if_dst, dst, dst_stride, pu_width, pu_height);
}

void luma_interpolation_filter_posp_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    EbByte               dst,
    uint32_t                dst_stride,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst)
{
    LumaInterpolationFilterOneDOutRawHorizontal_SSSE3(ref_pic - 2 * src_stride, src_stride, first_pass_if_dst, pu_width, pu_height + 6, 1);
    LumaInterpolationFilterTwoDInRaw7_SSSE3(first_pass_if_dst, dst, dst_stride, pu_width, pu_height, 3);
}

void luma_interpolation_filter_posq_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    EbByte               dst,
    uint32_t                dst_stride,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst)
{
    uint32_t puHeight1 = pu_height + 6;
    EbByte refPic1 = ref_pic - 2 * src_stride;
    LumaInterpolationFilterOneDOutRawHorizontal_SSSE3(refPic1, src_stride, first_pass_if_dst, pu_width, puHeight1, 2);
    LumaInterpolationFilterTwoDInRaw7_SSSE3(first_pass_if_dst, dst, dst_stride, pu_width, pu_height, 3);
}

void luma_interpolation_filter_posr_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    EbByte               dst,
    uint32_t                dst_stride,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst)
{
    LumaInterpolationFilterOneDOutRawHorizontal_SSSE3(ref_pic - 2 * src_stride, src_stride, first_pass_if_dst, pu_width, pu_height + 6, 3);
    LumaInterpolationFilterTwoDInRaw7_SSSE3(first_pass_if_dst, dst, dst_stride, pu_width, pu_height, 3);
}


void luma_interpolation_copy_out_raw_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    int16_t               *dst,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst)
{
    (void)first_pass_if_dst;
    PictureCopyKernelOutRaw_SSSE3(ref_pic, src_stride, dst, pu_width, pu_height, 128 * 64);
}



void luma_interpolation_filter_posa_out_raw_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    int16_t*               dst,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst)
{
    (void)first_pass_if_dst;
    LumaInterpolationFilterOneDOutRawHorizontal_SSSE3(ref_pic, src_stride, dst, pu_width, pu_height, 1);
}

void luma_interpolation_filter_posb_out_raw_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    int16_t*               dst,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst)
{
    (void)first_pass_if_dst;
    //LumaInterpolationFilterOneDOutRawHorizontalOut_SSSE3(ref_pic, src_stride, dst, pu_width, pu_height, 2);
    LumaInterpolationFilterOneDOutRawHorizontal_SSSE3(ref_pic, src_stride, dst, pu_width, pu_height, 2);
}

void luma_interpolation_filter_posc_out_raw_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    int16_t*               dst,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst)
{
    (void)first_pass_if_dst;
    //LumaInterpolationFilterOneDOutRawHorizontalOut_SSSE3(ref_pic, src_stride, dst, pu_width, pu_height, 3);
    LumaInterpolationFilterOneDOutRawHorizontal_SSSE3(ref_pic, src_stride, dst, pu_width, pu_height, 3);
}

void luma_interpolation_filter_posd_out_raw_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    int16_t*               dst,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst)
{
    (void)first_pass_if_dst;

    LumaInterpolationFilterOneDOutRawVertical_SSSE3(ref_pic, src_stride, dst, pu_width, pu_height, 1);
}

void luma_interpolation_filter_posh_out_raw_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    int16_t*               dst,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst)
{
    (void)first_pass_if_dst;

    LumaInterpolationFilterOneDOutRawVertical_SSSE3(ref_pic, src_stride, dst, pu_width, pu_height, 2);
}



void luma_interpolation_filter_posn_out_raw_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    int16_t*               dst,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst)
{
    (void)first_pass_if_dst;

    LumaInterpolationFilterOneDOutRawVertical_SSSE3(ref_pic, src_stride, dst, pu_width, pu_height, 3);
}

void luma_interpolation_filter_pose_out_raw_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    int16_t*               dst,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst)
{
    LumaInterpolationFilterOneDOutRawHorizontal_SSSE3(ref_pic - 3 * src_stride, src_stride, first_pass_if_dst, pu_width, pu_height + 6, 1);
    LumaInterpolationFilterTwoDInRawOutRaw7_SSSE3(first_pass_if_dst, dst, pu_width, pu_height, 1);
}

void luma_interpolation_filter_posf_out_raw_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    int16_t*               dst,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst)
{
    uint32_t puHeight1 = pu_height + 6;
    EbByte refPic1 = ref_pic - 3 * src_stride;
    LumaInterpolationFilterOneDOutRawHorizontal_SSSE3(refPic1, src_stride, first_pass_if_dst, pu_width, puHeight1, 2);
    LumaInterpolationFilterTwoDInRawOutRaw7_SSSE3(first_pass_if_dst, dst, pu_width, pu_height, 1);
}

void luma_interpolation_filter_posg_out_raw_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    int16_t*               dst,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst)
{
    LumaInterpolationFilterOneDOutRawHorizontal_SSSE3(ref_pic - 3 * src_stride, src_stride, first_pass_if_dst, pu_width, pu_height + 6, 3);
    LumaInterpolationFilterTwoDInRawOutRaw7_SSSE3(first_pass_if_dst, dst, pu_width, pu_height, 1);
}

void luma_interpolation_filter_posi_out_raw_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    int16_t*               dst,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst)
{
    LumaInterpolationFilterOneDOutRawHorizontal_SSSE3(ref_pic - 3 * src_stride, src_stride, first_pass_if_dst, pu_width, pu_height + 7, 1);
    LumaInterpolationFilterTwoDInRawOutRawM_SSSE3(first_pass_if_dst, dst, pu_width, pu_height);
}

void luma_interpolation_filter_posj_out_raw_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    int16_t*               dst,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst)
{
    LumaInterpolationFilterOneDOutRawHorizontal_SSSE3(ref_pic - 3 * src_stride, src_stride, first_pass_if_dst, pu_width, pu_height + 7, 2);
    LumaInterpolationFilterTwoDInRawOutRawM_SSSE3(first_pass_if_dst, dst, pu_width, pu_height);
}

void luma_interpolation_filter_posk_out_raw_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    int16_t*               dst,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst)
{
    LumaInterpolationFilterOneDOutRawHorizontal_SSSE3(ref_pic - 3 * src_stride, src_stride, first_pass_if_dst, pu_width, pu_height + 7, 3);
    LumaInterpolationFilterTwoDInRawOutRawM_SSSE3(first_pass_if_dst, dst, pu_width, pu_height);
}

void luma_interpolation_filter_posp_out_raw_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    int16_t*               dst,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst)
{
    LumaInterpolationFilterOneDOutRawHorizontal_SSSE3(ref_pic - 2 * src_stride, src_stride, first_pass_if_dst, pu_width, pu_height + 6, 1);
    LumaInterpolationFilterTwoDInRawOutRaw7_SSSE3(first_pass_if_dst, dst, pu_width, pu_height, 3);
}

void luma_interpolation_filter_posq_out_raw_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    int16_t*               dst,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst)
{
    uint32_t puHeight1 = pu_height + 6;
    EbByte refPic1 = ref_pic - 2 * src_stride;
    LumaInterpolationFilterOneDOutRawHorizontal_SSSE3(refPic1, src_stride, first_pass_if_dst, pu_width, puHeight1, 2);
    LumaInterpolationFilterTwoDInRawOutRaw7_SSSE3(first_pass_if_dst, dst, pu_width, pu_height, 3);
}

void luma_interpolation_filter_posr_out_raw_ssse3(
    EbByte               ref_pic,
    uint32_t                src_stride,
    int16_t*               dst,
    uint32_t                pu_width,
    uint32_t                pu_height,
    int16_t               *first_pass_if_dst)
{
    LumaInterpolationFilterOneDOutRawHorizontal_SSSE3(ref_pic - 2 * src_stride, src_stride, first_pass_if_dst, pu_width, pu_height + 6, 3);
    LumaInterpolationFilterTwoDInRawOutRaw7_SSSE3(first_pass_if_dst, dst, pu_width, pu_height, 3);
}


