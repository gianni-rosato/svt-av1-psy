/*
 * Copyright (c) 2018, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include "EbDefinitions.h"
#include <immintrin.h>
#include <assert.h>

#include "aom_dsp_rtcd.h"
#include "convolve.h"
#include "synonyms.h"
#include "synonyms_avx2.h"


DECLARE_ALIGNED(32, static const uint8_t, filt_center_global_avx2[32]) = {
  3, 255, 4, 255, 5, 255, 6, 255, 7, 255, 8, 255, 9, 255, 10, 255,
  3, 255, 4, 255, 5, 255, 6, 255, 7, 255, 8, 255, 9, 255, 10, 255
};

DECLARE_ALIGNED(32, static const uint8_t, filt1_global_avx2[32]) = {
  0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8,
  0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8
};

DECLARE_ALIGNED(32, static const uint8_t, filt2_global_avx2[32]) = {
  2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10,
  2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10
};

DECLARE_ALIGNED(32, static const uint8_t, filt3_global_avx2[32]) = {
  4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12,
  4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12
};

DECLARE_ALIGNED(32, static const uint8_t, filt4_global_avx2[32]) = {
  6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14,
  6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14
};

static INLINE __m256i convolve_lowbd(const __m256i *const s,
    const __m256i *const coeffs) {
    const __m256i res_01 = _mm256_maddubs_epi16(s[0], coeffs[0]);
    const __m256i res_23 = _mm256_maddubs_epi16(s[1], coeffs[1]);
    const __m256i res_45 = _mm256_maddubs_epi16(s[2], coeffs[2]);
    const __m256i res_67 = _mm256_maddubs_epi16(s[3], coeffs[3]);

    // order: 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
    const __m256i res = _mm256_add_epi16(_mm256_add_epi16(res_01, res_45),
        _mm256_add_epi16(res_23, res_67));

    return res;
}

static INLINE __m256i convolve_lowbd_x(const __m256i data,
    const __m256i *const coeffs,
    const __m256i *const filt) {
    __m256i s[4];

    s[0] = _mm256_shuffle_epi8(data, filt[0]);
    s[1] = _mm256_shuffle_epi8(data, filt[1]);
    s[2] = _mm256_shuffle_epi8(data, filt[2]);
    s[3] = _mm256_shuffle_epi8(data, filt[3]);

    return convolve_lowbd(s, coeffs);
}

static INLINE __m256i convolve(const __m256i *const s,
    const __m256i *const coeffs) {
    const __m256i res_0 = _mm256_madd_epi16(s[0], coeffs[0]);
    const __m256i res_1 = _mm256_madd_epi16(s[1], coeffs[1]);
    const __m256i res_2 = _mm256_madd_epi16(s[2], coeffs[2]);
    const __m256i res_3 = _mm256_madd_epi16(s[3], coeffs[3]);

    const __m256i res = _mm256_add_epi32(_mm256_add_epi32(res_0, res_1),
        _mm256_add_epi32(res_2, res_3));

    return res;
}

 // 128-bit xmmwords are written as [ ... ] with the MSB on the left.
 // 256-bit ymmwords are written as two xmmwords, [ ... ][ ... ] with the MSB
 // on the left.
 // A row of, say, 8-bit pixels with values p0, p1, p2, ..., p30, p31 will be
 // loaded and stored as [ p31 ... p17 p16 ][ p15 ... p1 p0 ].

// Exploiting the range of wiener filter coefficients,
// horizontal filtering can be done in 16 bit intermediate precision.
// The details are as follows :
// Consider the horizontal wiener filter coefficients of the following form :
//      [C0, C1, C2, 2^(FILTER_BITS) -2 * (C0 + C1 + C2), C2, C1, C0]
// Subtracting  2^(FILTER_BITS) from the centre tap we get the following  :
//      [C0, C1, C2,     -2 * (C0 + C1 + C2),             C2, C1, C0]
// The sum of the product "C0 * p0 + C1 * p1 + C2 * p2 -2 * (C0 + C1 + C2) * p3
// + C2 * p4 + C1 * p5 + C0 * p6" would be in the range of signed 16 bit
// precision. Finally, after rounding the above result by round_0, we multiply
// the centre pixel by 2^(FILTER_BITS - round_0) and add it to get the
// horizontal filter output.

void av1_wiener_convolve_add_src_avx2(const uint8_t *src, ptrdiff_t src_stride,
    uint8_t *dst, ptrdiff_t dst_stride,
    const int16_t *filter_x, int32_t x_step_q4,
    const int16_t *filter_y, int32_t y_step_q4,
    int32_t w, int32_t h,
    const ConvolveParams *conv_params) {
    assert(x_step_q4 == 16 && y_step_q4 == 16);
    assert(!(w & 7));
    (void)x_step_q4;
    (void)y_step_q4;

    DECLARE_ALIGNED(32, int16_t, im_block[(MAX_SB_SIZE + SUBPEL_TAPS) * 8]);
    int32_t im_h = h + SUBPEL_TAPS - 1;
    int32_t im_stride = 8;
    int32_t i, j;
    const int32_t center_tap = (SUBPEL_TAPS - 1) / 2;
    const uint8_t *const src_ptr = src - center_tap * src_stride - center_tap;

    __m256i filt[4], coeffs_h[4], coeffs_v[4], filt_center;

    assert(conv_params->round_0 > 0);

    filt[0] = _mm256_load_si256((__m256i const *)filt1_global_avx2);
    filt[1] = _mm256_load_si256((__m256i const *)filt2_global_avx2);
    filt[2] = _mm256_load_si256((__m256i const *)filt3_global_avx2);
    filt[3] = _mm256_load_si256((__m256i const *)filt4_global_avx2);

    filt_center = _mm256_load_si256((__m256i const *)filt_center_global_avx2);

    const __m128i coeffs_x = _mm_loadu_si128((__m128i *)filter_x);
    const __m256i filter_coeffs_x = _mm256_broadcastsi128_si256(coeffs_x);

    // coeffs 0 1 0 1 0 1 0 1
    coeffs_h[0] =
        _mm256_shuffle_epi8(filter_coeffs_x, _mm256_set1_epi16(0x0200u));
    // coeffs 2 3 2 3 2 3 2 3
    coeffs_h[1] =
        _mm256_shuffle_epi8(filter_coeffs_x, _mm256_set1_epi16(0x0604u));
    // coeffs 4 5 4 5 4 5 4 5
    coeffs_h[2] =
        _mm256_shuffle_epi8(filter_coeffs_x, _mm256_set1_epi16(0x0a08u));
    // coeffs 6 7 6 7 6 7 6 7
    coeffs_h[3] =
        _mm256_shuffle_epi8(filter_coeffs_x, _mm256_set1_epi16(0x0e0cu));

    const __m256i round_const_h =
        _mm256_set1_epi16((1 << (conv_params->round_0 - 1)));
    const __m128i round_shift_h = _mm_cvtsi32_si128(conv_params->round_0);

    // Add an offset to account for the "add_src" part of the convolve function.
    const __m128i zero_128 = _mm_setzero_si128();
    const __m128i offset_0 = _mm_insert_epi16(zero_128, 1 << FILTER_BITS, 3);
    const __m128i coeffs_y = _mm_add_epi16(xx_loadu_128(filter_y), offset_0);

    const __m256i filter_coeffs_y = _mm256_broadcastsi128_si256(coeffs_y);

    // coeffs 0 1 0 1 0 1 0 1
    coeffs_v[0] = _mm256_shuffle_epi32(filter_coeffs_y, 0x00);
    // coeffs 2 3 2 3 2 3 2 3
    coeffs_v[1] = _mm256_shuffle_epi32(filter_coeffs_y, 0x55);
    // coeffs 4 5 4 5 4 5 4 5
    coeffs_v[2] = _mm256_shuffle_epi32(filter_coeffs_y, 0xaa);
    // coeffs 6 7 6 7 6 7 6 7
    coeffs_v[3] = _mm256_shuffle_epi32(filter_coeffs_y, 0xff);

    const __m256i round_const_v =
        _mm256_set1_epi32((1 << (conv_params->round_1 - 1)));
    const __m128i round_shift_v = _mm_cvtsi32_si128(conv_params->round_1);

    for (j = 0; j < w; j += 8) {
        for (i = 0; i < im_h; i += 2) {
            __m256i data = _mm256_castsi128_si256(
                _mm_loadu_si128((__m128i *)&src_ptr[(i * src_stride) + j]));

            // Load the next line
            if (i + 1 < im_h)
                data = _mm256_inserti128_si256(
                    data,
                    _mm_loadu_si128(
                    (__m128i *)&src_ptr[(i * src_stride) + j + src_stride]),
                    1);

            __m256i res = convolve_lowbd_x(data, coeffs_h, filt);

            res =
                _mm256_sra_epi16(_mm256_add_epi16(res, round_const_h), round_shift_h);

            __m256i data_0 = _mm256_shuffle_epi8(data, filt_center);

            // multiply the center pixel by 2^(FILTER_BITS - round_0) and add it to
            // the result
            data_0 = _mm256_slli_epi16(data_0, FILTER_BITS - conv_params->round_0);
            res = _mm256_add_epi16(res, data_0);

            _mm256_store_si256((__m256i *)&im_block[i * im_stride], res);
        }

        /* Vertical filter */
        {
            __m256i src_0 = _mm256_loadu_si256((__m256i *)(im_block + 0 * im_stride));
            __m256i src_1 = _mm256_loadu_si256((__m256i *)(im_block + 1 * im_stride));
            __m256i src_2 = _mm256_loadu_si256((__m256i *)(im_block + 2 * im_stride));
            __m256i src_3 = _mm256_loadu_si256((__m256i *)(im_block + 3 * im_stride));
            __m256i src_4 = _mm256_loadu_si256((__m256i *)(im_block + 4 * im_stride));
            __m256i src_5 = _mm256_loadu_si256((__m256i *)(im_block + 5 * im_stride));

            __m256i s[8];
            s[0] = _mm256_unpacklo_epi16(src_0, src_1);
            s[1] = _mm256_unpacklo_epi16(src_2, src_3);
            s[2] = _mm256_unpacklo_epi16(src_4, src_5);

            s[4] = _mm256_unpackhi_epi16(src_0, src_1);
            s[5] = _mm256_unpackhi_epi16(src_2, src_3);
            s[6] = _mm256_unpackhi_epi16(src_4, src_5);

            for (i = 0; i < h - 1; i += 2) {
                const int16_t *data = &im_block[i * im_stride];

                const __m256i s6 =
                    _mm256_loadu_si256((__m256i *)(data + 6 * im_stride));
                const __m256i s7 =
                    _mm256_loadu_si256((__m256i *)(data + 7 * im_stride));

                s[3] = _mm256_unpacklo_epi16(s6, s7);
                s[7] = _mm256_unpackhi_epi16(s6, s7);

                __m256i res_a = convolve(s, coeffs_v);
                __m256i res_b = convolve(s + 4, coeffs_v);

                const __m256i res_a_round = _mm256_sra_epi32(
                    _mm256_add_epi32(res_a, round_const_v), round_shift_v);
                const __m256i res_b_round = _mm256_sra_epi32(
                    _mm256_add_epi32(res_b, round_const_v), round_shift_v);

                /* rounding code */
                // 16 bit conversion
                const __m256i res_16bit = _mm256_packs_epi32(res_a_round, res_b_round);
                // 8 bit conversion and saturation to uint8
                const __m256i res_8b = _mm256_packus_epi16(res_16bit, res_16bit);

                const __m128i res_0 = _mm256_castsi256_si128(res_8b);
                const __m128i res_1 = _mm256_extracti128_si256(res_8b, 1);

                // Store values into the destination buffer
                __m128i *const p_0 = (__m128i *)&dst[i * dst_stride + j];
                __m128i *const p_1 = (__m128i *)&dst[i * dst_stride + j + dst_stride];

                _mm_storel_epi64(p_0, res_0);
                _mm_storel_epi64(p_1, res_1);

                s[0] = s[1];
                s[1] = s[2];
                s[2] = s[3];

                s[4] = s[5];
                s[5] = s[6];
                s[6] = s[7];
            }
            if (h - i) {
                s[0] = _mm256_permute2x128_si256(s[0], s[4], 0x20);
                s[1] = _mm256_permute2x128_si256(s[1], s[5], 0x20);
                s[2] = _mm256_permute2x128_si256(s[2], s[6], 0x20);

                const int16_t *data = &im_block[i * im_stride];
                const __m128i s6_ = _mm_loadu_si128((__m128i *)(data + 6 * im_stride));
                const __m128i s7_ = _mm_loadu_si128((__m128i *)(data + 7 * im_stride));

                __m128i s3 = _mm_unpacklo_epi16(s6_, s7_);
                __m128i s7 = _mm_unpackhi_epi16(s6_, s7_);

                s[3] = _mm256_inserti128_si256(_mm256_castsi128_si256(s3), s7, 1);
                __m256i convolveres = convolve(s, coeffs_v);

                const __m256i res_round = _mm256_sra_epi32(
                    _mm256_add_epi32(convolveres, round_const_v), round_shift_v);

                /* rounding code */
                // 16 bit conversion
                __m128i reslo = _mm256_castsi256_si128(res_round);
                __m128i reshi = _mm256_extracti128_si256(res_round, 1);
                const __m128i res_16bit = _mm_packus_epi32(reslo, reshi);

                // 8 bit conversion and saturation to uint8
                const __m128i res_8b = _mm_packus_epi16(res_16bit, res_16bit);
                __m128i *const p_0 = (__m128i *)&dst[i * dst_stride + j];
                _mm_storel_epi64(p_0, res_8b);
            }
        }
    }
}

// 128-bit xmmwords are written as [ ... ] with the MSB on the left.
// 256-bit ymmwords are written as two xmmwords, [ ... ][ ... ] with the MSB
// on the left.
// A row of, say, 16-bit pixels with values p0, p1, p2, ..., p14, p15 will be
// loaded and stored as [ p15 ... p9 p8 ][ p7 ... p1 p0 ].
void av1_highbd_wiener_convolve_add_src_avx2(
    const uint8_t *src8, ptrdiff_t src_stride, uint8_t *dst8,
    ptrdiff_t dst_stride, const int16_t *filter_x, int32_t x_step_q4,
    const int16_t *filter_y, int32_t y_step_q4, int32_t w, int32_t h,
    const ConvolveParams *conv_params, int32_t bd) {
    assert(x_step_q4 == 16 && y_step_q4 == 16);
    assert(!(w & 7));
    assert(bd + FILTER_BITS - conv_params->round_0 + 2 <= 16);
    (void)x_step_q4;
    (void)y_step_q4;

    const uint16_t *const src = CONVERT_TO_SHORTPTR(src8);
    uint16_t *const dst = CONVERT_TO_SHORTPTR(dst8);

    DECLARE_ALIGNED(32, uint16_t,
    temp[(MAX_SB_SIZE + SUBPEL_TAPS - 1) * MAX_SB_SIZE]);
    int32_t intermediate_height = h + SUBPEL_TAPS - 1;
    const int32_t center_tap = ((SUBPEL_TAPS - 1) / 2);
    const uint16_t *const src_ptr = src - center_tap * src_stride - center_tap;

    const __m128i zero_128 = _mm_setzero_si128();
    const __m256i zero_256 = _mm256_setzero_si256();

    // Add an offset to account for the "add_src" part of the convolve function.
    const __m128i offset = _mm_insert_epi16(zero_128, 1 << FILTER_BITS, 3);

    const __m256i clamp_low = zero_256;

    /* Horizontal filter */
    {
        const __m256i clamp_high_ep =
            _mm256_set1_epi16(WIENER_CLAMP_LIMIT(conv_params->round_0, bd) - 1);

        // coeffs [ f7 f6 f5 f4 f3 f2 f1 f0 ]
        const __m128i coeffs_x = _mm_add_epi16(xx_loadu_128(filter_x), offset);

        // coeffs [ f3 f2 f3 f2 f1 f0 f1 f0 ]
        const __m128i coeffs_0123 = _mm_unpacklo_epi32(coeffs_x, coeffs_x);
        // coeffs [ f7 f6 f7 f6 f5 f4 f5 f4 ]
        const __m128i coeffs_4567 = _mm_unpackhi_epi32(coeffs_x, coeffs_x);

        // coeffs [ f1 f0 f1 f0 f1 f0 f1 f0 ]
        const __m128i coeffs_01_128 = _mm_unpacklo_epi64(coeffs_0123, coeffs_0123);
        // coeffs [ f3 f2 f3 f2 f3 f2 f3 f2 ]
        const __m128i coeffs_23_128 = _mm_unpackhi_epi64(coeffs_0123, coeffs_0123);
        // coeffs [ f5 f4 f5 f4 f5 f4 f5 f4 ]
        const __m128i coeffs_45_128 = _mm_unpacklo_epi64(coeffs_4567, coeffs_4567);
        // coeffs [ f7 f6 f7 f6 f7 f6 f7 f6 ]
        const __m128i coeffs_67_128 = _mm_unpackhi_epi64(coeffs_4567, coeffs_4567);

        // coeffs [ f1 f0 f1 f0 f1 f0 f1 f0 ][ f1 f0 f1 f0 f1 f0 f1 f0 ]
        const __m256i coeffs_01 = yy_set_m128i(coeffs_01_128, coeffs_01_128);
        // coeffs [ f3 f2 f3 f2 f3 f2 f3 f2 ][ f3 f2 f3 f2 f3 f2 f3 f2 ]
        const __m256i coeffs_23 = yy_set_m128i(coeffs_23_128, coeffs_23_128);
        // coeffs [ f5 f4 f5 f4 f5 f4 f5 f4 ][ f5 f4 f5 f4 f5 f4 f5 f4 ]
        const __m256i coeffs_45 = yy_set_m128i(coeffs_45_128, coeffs_45_128);
        // coeffs [ f7 f6 f7 f6 f7 f6 f7 f6 ][ f7 f6 f7 f6 f7 f6 f7 f6 ]
        const __m256i coeffs_67 = yy_set_m128i(coeffs_67_128, coeffs_67_128);

        const __m256i round_const = _mm256_set1_epi32(
            (1 << (conv_params->round_0 - 1)) + (1 << (bd + FILTER_BITS - 1)));

        for (int32_t i = 0; i < intermediate_height; ++i) {
            for (int32_t j = 0; j < w; j += 16) {
                const uint16_t *src_ij = src_ptr + i * src_stride + j;

                // Load 16-bit src data
                const __m256i src_0 = yy_loadu_256(src_ij + 0);
                const __m256i src_1 = yy_loadu_256(src_ij + 1);
                const __m256i src_2 = yy_loadu_256(src_ij + 2);
                const __m256i src_3 = yy_loadu_256(src_ij + 3);
                const __m256i src_4 = yy_loadu_256(src_ij + 4);
                const __m256i src_5 = yy_loadu_256(src_ij + 5);
                const __m256i src_6 = yy_loadu_256(src_ij + 6);
                const __m256i src_7 = yy_loadu_256(src_ij + 7);

                // Multiply src data by filter coeffs and sum pairs
                const __m256i res_0 = _mm256_madd_epi16(src_0, coeffs_01);
                const __m256i res_1 = _mm256_madd_epi16(src_1, coeffs_01);
                const __m256i res_2 = _mm256_madd_epi16(src_2, coeffs_23);
                const __m256i res_3 = _mm256_madd_epi16(src_3, coeffs_23);
                const __m256i res_4 = _mm256_madd_epi16(src_4, coeffs_45);
                const __m256i res_5 = _mm256_madd_epi16(src_5, coeffs_45);
                const __m256i res_6 = _mm256_madd_epi16(src_6, coeffs_67);
                const __m256i res_7 = _mm256_madd_epi16(src_7, coeffs_67);

                // Calculate scalar product for even- and odd-indices separately,
                // increasing to 32-bit precision
                const __m256i res_even_sum = _mm256_add_epi32(
                    _mm256_add_epi32(res_0, res_4), _mm256_add_epi32(res_2, res_6));
                const __m256i res_even = _mm256_srai_epi32(
                    _mm256_add_epi32(res_even_sum, round_const), conv_params->round_0);

                const __m256i res_odd_sum = _mm256_add_epi32(
                    _mm256_add_epi32(res_1, res_5), _mm256_add_epi32(res_3, res_7));
                const __m256i res_odd = _mm256_srai_epi32(
                    _mm256_add_epi32(res_odd_sum, round_const), conv_params->round_0);

                // Reduce to 16-bit precision and pack even- and odd-index results
                // back into one register. The _mm256_packs_epi32 intrinsic returns
                // a register with the pixels ordered as follows:
                // [ 15 13 11 9 14 12 10 8 ] [ 7 5 3 1 6 4 2 0 ]
                const __m256i res = _mm256_packs_epi32(res_even, res_odd);
                const __m256i res_clamped =
                    _mm256_min_epi16(_mm256_max_epi16(res, clamp_low), clamp_high_ep);

                // Store in a temporary array
                yy_storeu_256(temp + i * MAX_SB_SIZE + j, res_clamped);
            }
        }
    }

    /* Vertical filter */
    {
        const __m256i clamp_high = _mm256_set1_epi16((1 << bd) - 1);

        // coeffs [ f7 f6 f5 f4 f3 f2 f1 f0 ]
        const __m128i coeffs_y = _mm_add_epi16(xx_loadu_128(filter_y), offset);

        // coeffs [ f3 f2 f3 f2 f1 f0 f1 f0 ]
        const __m128i coeffs_0123 = _mm_unpacklo_epi32(coeffs_y, coeffs_y);
        // coeffs [ f7 f6 f7 f6 f5 f4 f5 f4 ]
        const __m128i coeffs_4567 = _mm_unpackhi_epi32(coeffs_y, coeffs_y);

        // coeffs [ f1 f0 f1 f0 f1 f0 f1 f0 ]
        const __m128i coeffs_01_128 = _mm_unpacklo_epi64(coeffs_0123, coeffs_0123);
        // coeffs [ f3 f2 f3 f2 f3 f2 f3 f2 ]
        const __m128i coeffs_23_128 = _mm_unpackhi_epi64(coeffs_0123, coeffs_0123);
        // coeffs [ f5 f4 f5 f4 f5 f4 f5 f4 ]
        const __m128i coeffs_45_128 = _mm_unpacklo_epi64(coeffs_4567, coeffs_4567);
        // coeffs [ f7 f6 f7 f6 f7 f6 f7 f6 ]
        const __m128i coeffs_67_128 = _mm_unpackhi_epi64(coeffs_4567, coeffs_4567);

        // coeffs [ f1 f0 f1 f0 f1 f0 f1 f0 ][ f1 f0 f1 f0 f1 f0 f1 f0 ]
        const __m256i coeffs_01 = yy_set_m128i(coeffs_01_128, coeffs_01_128);
        // coeffs [ f3 f2 f3 f2 f3 f2 f3 f2 ][ f3 f2 f3 f2 f3 f2 f3 f2 ]
        const __m256i coeffs_23 = yy_set_m128i(coeffs_23_128, coeffs_23_128);
        // coeffs [ f5 f4 f5 f4 f5 f4 f5 f4 ][ f5 f4 f5 f4 f5 f4 f5 f4 ]
        const __m256i coeffs_45 = yy_set_m128i(coeffs_45_128, coeffs_45_128);
        // coeffs [ f7 f6 f7 f6 f7 f6 f7 f6 ][ f7 f6 f7 f6 f7 f6 f7 f6 ]
        const __m256i coeffs_67 = yy_set_m128i(coeffs_67_128, coeffs_67_128);

        const __m256i round_const =
            _mm256_set1_epi32((1 << (conv_params->round_1 - 1)) -
            (1 << (bd + conv_params->round_1 - 1)));

        for (int32_t i = 0; i < h; ++i) {
            for (int32_t j = 0; j < w; j += 16) {
                const uint16_t *temp_ij = temp + i * MAX_SB_SIZE + j;

                // Load 16-bit data from the output of the horizontal filter in
                // which the pixels are ordered as follows:
                // [ 15 13 11 9 14 12 10 8 ] [ 7 5 3 1 6 4 2 0 ]
                const __m256i data_0 = yy_loadu_256(temp_ij + 0 * MAX_SB_SIZE);
                const __m256i data_1 = yy_loadu_256(temp_ij + 1 * MAX_SB_SIZE);
                const __m256i data_2 = yy_loadu_256(temp_ij + 2 * MAX_SB_SIZE);
                const __m256i data_3 = yy_loadu_256(temp_ij + 3 * MAX_SB_SIZE);
                const __m256i data_4 = yy_loadu_256(temp_ij + 4 * MAX_SB_SIZE);
                const __m256i data_5 = yy_loadu_256(temp_ij + 5 * MAX_SB_SIZE);
                const __m256i data_6 = yy_loadu_256(temp_ij + 6 * MAX_SB_SIZE);
                const __m256i data_7 = yy_loadu_256(temp_ij + 7 * MAX_SB_SIZE);

                // Filter the even-indices, increasing to 32-bit precision
                const __m256i src_0 = _mm256_unpacklo_epi16(data_0, data_1);
                const __m256i src_2 = _mm256_unpacklo_epi16(data_2, data_3);
                const __m256i src_4 = _mm256_unpacklo_epi16(data_4, data_5);
                const __m256i src_6 = _mm256_unpacklo_epi16(data_6, data_7);

                const __m256i res_0 = _mm256_madd_epi16(src_0, coeffs_01);
                const __m256i res_2 = _mm256_madd_epi16(src_2, coeffs_23);
                const __m256i res_4 = _mm256_madd_epi16(src_4, coeffs_45);
                const __m256i res_6 = _mm256_madd_epi16(src_6, coeffs_67);

                const __m256i res_even = _mm256_add_epi32(
                    _mm256_add_epi32(res_0, res_2), _mm256_add_epi32(res_4, res_6));

                // Filter the odd-indices, increasing to 32-bit precision
                const __m256i src_1 = _mm256_unpackhi_epi16(data_0, data_1);
                const __m256i src_3 = _mm256_unpackhi_epi16(data_2, data_3);
                const __m256i src_5 = _mm256_unpackhi_epi16(data_4, data_5);
                const __m256i src_7 = _mm256_unpackhi_epi16(data_6, data_7);

                const __m256i res_1 = _mm256_madd_epi16(src_1, coeffs_01);
                const __m256i res_3 = _mm256_madd_epi16(src_3, coeffs_23);
                const __m256i res_5 = _mm256_madd_epi16(src_5, coeffs_45);
                const __m256i res_7 = _mm256_madd_epi16(src_7, coeffs_67);

                const __m256i res_odd = _mm256_add_epi32(
                    _mm256_add_epi32(res_1, res_3), _mm256_add_epi32(res_5, res_7));

                // Pixels are currently in the following order:
                // res_even order: [ 14 12 10 8 ] [ 6 4 2 0 ]
                // res_odd order:  [ 15 13 11 9 ] [ 7 5 3 1 ]
                //
                // Rearrange the pixels into the following order:
                // res_lo order: [ 11 10  9  8 ] [ 3 2 1 0 ]
                // res_hi order: [ 15 14 13 12 ] [ 7 6 5 4 ]
                const __m256i res_lo = _mm256_unpacklo_epi32(res_even, res_odd);
                const __m256i res_hi = _mm256_unpackhi_epi32(res_even, res_odd);

                const __m256i res_lo_round = _mm256_srai_epi32(
                    _mm256_add_epi32(res_lo, round_const), conv_params->round_1);
                const __m256i res_hi_round = _mm256_srai_epi32(
                    _mm256_add_epi32(res_hi, round_const), conv_params->round_1);

                // Reduce to 16-bit precision and pack into the correct order:
                // [ 15 14 13 12 11 10 9 8 ][ 7 6 5 4 3 2 1 0 ]
                const __m256i res_16bit =
                    _mm256_packs_epi32(res_lo_round, res_hi_round);
                const __m256i res_16bit_clamped = _mm256_min_epi16(
                    _mm256_max_epi16(res_16bit, clamp_low), clamp_high);

                // Store in the dst array
                yy_storeu_256(dst + i * dst_stride + j, res_16bit_clamped);
            }
        }
    }
}



