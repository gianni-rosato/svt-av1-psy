/*
 * Copyright (c) 2017, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <immintrin.h>
#include "EbDefinitions.h"
#include "EbMemory_SSE4_1.h"
#include "aom_dsp_rtcd.h"
#include "convolve.h"
#include "convolve_avx2.h"
#include "synonyms.h"

void eb_av1_convolve_2d_sr_avx2(const uint8_t *src, int32_t src_stride,
    uint8_t *dst, int32_t dst_stride, int32_t w,
    int32_t h, InterpFilterParams *filter_params_x,
    InterpFilterParams *filter_params_y,
    const int32_t subpel_x_q4,
    const int32_t subpel_y_q4,
    ConvolveParams *conv_params) {
    const int32_t tap_y = get_convolve_tap(filter_params_y->filter_ptr);
    int32_t y;
    const uint8_t *src_ptr;
    int16_t *im;
    // Note: im_block is 8-pixel interlaced for width 32 and up, to avoid data
    //       permutation.
    DECLARE_ALIGNED(
    32, int16_t, im_block[(MAX_SB_SIZE + MAX_FILTER_TAP) * MAX_SB_SIZE]);

    assert(conv_params->round_0 == 3);
    assert(conv_params->round_1 == 11);

    (void)conv_params;

    // horizontal filter
    __m128i coeffs_128[4];
    __m256i coeffs_256[4];

    src_ptr = src + ((8 - tap_y) / 2 - 3) * src_stride;
    im = im_block;
    // Have to calculate 1 more row for small widths, since 2 lines are
    // calculated in each loop for them.
    y = h + tap_y - (w >= 32);

    if (is_convolve_2tap(filter_params_x->filter_ptr)) {
        // horz_filt as 2 tap
        if (w <= 8) {
            prepare_half_coeffs_2tap_ssse3(
                filter_params_x, subpel_x_q4, coeffs_128);

            if (w == 2) {
                const __m128i c = _mm_setr_epi8(
                    0, 1, 1, 2, 4, 5, 5, 6, 0, 1, 1, 2, 4, 5, 5, 6);

                do {
                    __m128i s_128;

                    s_128 = load_u8_4x2_sse4_1(src_ptr, src_stride);
                    s_128 = _mm_shuffle_epi8(s_128, c);
                    const __m128i res = convolve_2tap_ssse3(&s_128, coeffs_128);
                    const __m128i d = convolve_2dx_round_sse2(res);
                    _mm_storel_epi64((__m128i *)im, d);

                    src_ptr += 2 * src_stride;
                    im += 2 * 2;
                    y -= 2;
                } while (y);
            }
            else if (w == 4) {
                const __m128i c = _mm_setr_epi8(
                    0, 1, 1, 2, 2, 3, 3, 4, 8, 9, 9, 10, 10, 11, 11, 12);

                do {
                    __m128i s_128;

                    s_128 = load_u8_8x2_sse2(src_ptr, src_stride);
                    s_128 = _mm_shuffle_epi8(s_128, c);
                    const __m128i res = convolve_2tap_ssse3(&s_128, coeffs_128);
                    const __m128i d = convolve_2dx_round_sse2(res);
                    _mm_store_si128((__m128i *)im, d);

                    src_ptr += 2 * src_stride;
                    im += 2 * 4;
                    y -= 2;
                } while (y);
            }
            else {
                assert(w == 8);

                do {
                    __m128i s_128[2], res[2];

                    const __m128i s00 = _mm_loadu_si128((__m128i *)src_ptr);
                    const __m128i s10 =
                        _mm_loadu_si128((__m128i *)(src_ptr + src_stride));
                    const __m128i s01 = _mm_srli_si128(s00, 1);
                    const __m128i s11 = _mm_srli_si128(s10, 1);
                    s_128[0] = _mm_unpacklo_epi8(s00, s01);
                    s_128[1] = _mm_unpacklo_epi8(s10, s11);

                    res[0] = convolve_2tap_ssse3(&s_128[0], coeffs_128);
                    res[1] = convolve_2tap_ssse3(&s_128[1], coeffs_128);
                    res[0] = convolve_2dx_round_sse2(res[0]);
                    res[1] = convolve_2dx_round_sse2(res[1]);
                    _mm_store_si128((__m128i *)im, res[0]);
                    _mm_store_si128((__m128i *)(im + 8), res[1]);

                    src_ptr += 2 * src_stride;
                    im += 2 * 8;
                    y -= 2;
                } while (y);
            }
        }
        else {
            prepare_half_coeffs_2tap_avx2(
                filter_params_x, subpel_x_q4, coeffs_256);

            if (w == 16) {
                do {
                    __m128i s_128[2][2];
                    __m256i s_256[2];

                    s_128[0][0] = _mm_loadu_si128((__m128i *)src_ptr);
                    s_128[0][1] = _mm_loadu_si128((__m128i *)(src_ptr + 1));
                    s_128[1][0] =
                        _mm_loadu_si128((__m128i *)(src_ptr + src_stride));
                    s_128[1][1] =
                        _mm_loadu_si128((__m128i *)(src_ptr + src_stride + 1));
                    s_256[0] = _mm256_setr_m128i(s_128[0][0], s_128[1][0]);
                    s_256[1] = _mm256_setr_m128i(s_128[0][1], s_128[1][1]);
                    const __m256i s0 = _mm256_unpacklo_epi8(s_256[0], s_256[1]);
                    const __m256i s1 = _mm256_unpackhi_epi8(s_256[0], s_256[1]);
                    const __m256i res0 = convolve_2tap_avx2(&s0, coeffs_256);
                    const __m256i res1 = convolve_2tap_avx2(&s1, coeffs_256);
                    const __m256i r0 = convolve_2dx_round_avx2(res0);
                    const __m256i r1 = convolve_2dx_round_avx2(res1);
                    const __m256i d0 = _mm256_inserti128_si256(
                        r0, _mm256_extracti128_si256(r1, 0), 1);
                    const __m256i d1 = _mm256_inserti128_si256(
                        r1, _mm256_extracti128_si256(r0, 1), 0);
                    _mm256_store_si256((__m256i *)im, d0);
                    _mm256_store_si256((__m256i *)(im + 16), d1);

                    src_ptr += 2 * src_stride;
                    im += 2 * 16;
                    y -= 2;
                } while (y);
            }
            else if (w == 32) {
                do {
                    convolve_2dx_2tap_32_avx2(src_ptr, coeffs_256, im);
                    src_ptr += src_stride;
                    im += 32;
                } while (--y);
            }
            else if (w == 64) {
                do {
                    convolve_2dx_2tap_32_avx2(
                        src_ptr + 0 * 32, coeffs_256, im + 0 * 32);
                    convolve_2dx_2tap_32_avx2(
                        src_ptr + 1 * 32, coeffs_256, im + 1 * 32);
                    src_ptr += src_stride;
                    im += 64;
                } while (--y);
            }
            else {
                assert(w == 128);

                do {
                    convolve_2dx_2tap_32_avx2(
                        src_ptr + 0 * 32, coeffs_256, im + 0 * 32);
                    convolve_2dx_2tap_32_avx2(
                        src_ptr + 1 * 32, coeffs_256, im + 1 * 32);
                    convolve_2dx_2tap_32_avx2(
                        src_ptr + 2 * 32, coeffs_256, im + 2 * 32);
                    convolve_2dx_2tap_32_avx2(
                        src_ptr + 3 * 32, coeffs_256, im + 3 * 32);
                    src_ptr += src_stride;
                    im += 128;
                } while (--y);
            }
        }
    }
    else if (is_convolve_4tap(filter_params_x->filter_ptr)) {
        // horz_filt as 4 tap
        src_ptr -= 1;
        const __m128i c0 =
            _mm_setr_epi8(0, 1, 1, 2, 2, 3, 3, 4, 8, 9, 9, 10, 10, 11, 11, 12);
        const __m128i c1 = _mm_setr_epi8(
            2, 3, 3, 4, 4, 5, 5, 6, 10, 11, 11, 12, 12, 13, 13, 14);
        __m128i t, s_128[2];

        prepare_half_coeffs_4tap_ssse3(
            filter_params_x, subpel_x_q4, coeffs_128);

        if (w == 2) {
            do {
                t = load_u8_8x2_sse2(src_ptr, src_stride);
                s_128[0] = _mm_shuffle_epi8(t, c0);
                s_128[1] = _mm_shuffle_epi8(t, c1);
                const __m128i res = convolve_4tap_ssse3(s_128, coeffs_128);
                const __m128i d = convolve_2dx_round_sse2(res);
                *(int32_t *)im = _mm_cvtsi128_si32(d);
                *(int32_t *)(im + 2) = _mm_extract_epi32(d, 2);

                src_ptr += 2 * src_stride;
                im += 2 * 2;
                y -= 2;
            } while (y);
        }
        else {
            assert(w == 4);

            do {
                t = load_u8_8x2_sse2(src_ptr, src_stride);
                s_128[0] = _mm_shuffle_epi8(t, c0);
                s_128[1] = _mm_shuffle_epi8(t, c1);
                const __m128i res = convolve_4tap_ssse3(s_128, coeffs_128);
                const __m128i d = convolve_2dx_round_sse2(res);
                _mm_store_si128((__m128i *)im, d);

                src_ptr += 2 * src_stride;
                im += 2 * 4;
                y -= 2;
            } while (y);
        }
    }
    else {
        __m256i filt_256[4];

        filt_256[0] = _mm256_load_si256((__m256i const *)filt1_global_avx2);
        filt_256[1] = _mm256_load_si256((__m256i const *)filt2_global_avx2);
        filt_256[2] = _mm256_load_si256((__m256i const *)filt3_global_avx2);

        if (is_convolve_6tap(filter_params_x->filter_ptr)) {
            // horz_filt as 6 tap
            src_ptr -= 2;

            prepare_half_coeffs_6tap_avx2(
                filter_params_x, subpel_x_q4, coeffs_256);

            if (w == 8) {
                do {
                    const __m128i s0 = _mm_loadu_si128((__m128i *)src_ptr);
                    const __m128i s1 =
                        _mm_loadu_si128((__m128i *)(src_ptr + src_stride));
                    const __m256i s_256 = _mm256_setr_m128i(s0, s1);
                    const __m256i res =
                        convolve_x_6tap_avx2(s_256, coeffs_256, filt_256);
                    const __m256i d = convolve_2dx_round_avx2(res);
                    _mm256_store_si256((__m256i *)im, d);
                    src_ptr += 2 * src_stride;
                    im += 2 * 8;
                    y -= 2;
                } while (y);
            }
            else if (w == 16) {
                do {
                    convolve_2dx_6tap_16x2_avx2(
                        src_ptr, src_stride, coeffs_256, filt_256, im);
                    src_ptr += 2 * src_stride;
                    im += 2 * 16;
                    y -= 2;
                } while (y);
            }
            else if (w == 32) {
                do {
                    convolve_2dx_6tap_32_avx2(
                        src_ptr, 16, coeffs_256, filt_256, im);
                    src_ptr += src_stride;
                    im += 32;
                } while (--y);
            }
            else if (w == 64) {
                do {
                    convolve_2dx_6tap_32_avx2(
                        src_ptr, 16, coeffs_256, filt_256, im);
                    convolve_2dx_6tap_32_avx2(
                        src_ptr + 32, 16, coeffs_256, filt_256, im + 32);
                    src_ptr += src_stride;
                    im += 64;
                } while (--y);
            }
            else {
                assert(w == 128);

                do {
                    convolve_2dx_6tap_32_avx2(
                        src_ptr, 16, coeffs_256, filt_256, im);
                    convolve_2dx_6tap_32_avx2(
                        src_ptr + 32, 16, coeffs_256, filt_256, im + 32);
                    convolve_2dx_6tap_32_avx2(
                        src_ptr + 64, 16, coeffs_256, filt_256, im + 64);
                    convolve_2dx_6tap_32_avx2(
                        src_ptr + 96, 16, coeffs_256, filt_256, im + 96);
                    src_ptr += src_stride;
                    im += 128;
                } while (--y);
            }
        }
        else {
            // horz_filt as 8 tap
            src_ptr -= 3;

            filt_256[3] = _mm256_load_si256((__m256i const *)filt4_global_avx2);

            prepare_half_coeffs_8tap_avx2(
                filter_params_x, subpel_x_q4, coeffs_256);

            if (w == 8) {
                do {
                    const __m128i s0 = _mm_loadu_si128((__m128i *)src_ptr);
                    const __m128i s1 =
                        _mm_loadu_si128((__m128i *)(src_ptr + src_stride));
                    const __m256i s_256 = _mm256_setr_m128i(s0, s1);
                    const __m256i res =
                        convolve_x_8tap_avx2(s_256, coeffs_256, filt_256);
                    const __m256i d = convolve_2dx_round_avx2(res);
                    _mm256_store_si256((__m256i *)im, d);
                    src_ptr += 2 * src_stride;
                    im += 2 * 8;
                    y -= 2;
                } while (y);
            }
            else if (w == 16) {
                do {
                    convolve_2dx_8tap_16x2_avx2(
                        src_ptr, src_stride, coeffs_256, filt_256, im);
                    src_ptr += 2 * src_stride;
                    im += 2 * 16;
                    y -= 2;
                } while (y);
            }
            else if (w == 32) {
                do {
                    convolve_2dx_8tap_32_avx2(
                        src_ptr, 16, coeffs_256, filt_256, im);
                    src_ptr += src_stride;
                    im += 32;
                } while (--y);
            }
            else if (w == 64) {
                do {
                    convolve_2dx_8tap_32_avx2(
                        src_ptr, 16, coeffs_256, filt_256, im);
                    convolve_2dx_8tap_32_avx2(
                        src_ptr + 32, 16, coeffs_256, filt_256, im + 32);
                    src_ptr += src_stride;
                    im += 64;
                } while (--y);
            }
            else {
                assert(w == 128);

                do {
                    convolve_2dx_8tap_32_avx2(
                        src_ptr, 16, coeffs_256, filt_256, im);
                    convolve_2dx_8tap_32_avx2(
                        src_ptr + 32, 16, coeffs_256, filt_256, im + 32);
                    convolve_2dx_8tap_32_avx2(
                        src_ptr + 64, 16, coeffs_256, filt_256, im + 64);
                    convolve_2dx_8tap_32_avx2(
                        src_ptr + 96, 16, coeffs_256, filt_256, im + 96);
                    src_ptr += src_stride;
                    im += 128;
                } while (--y);
            }
        }
    }

    // vertical filter
    im = im_block;

    if (tap_y == 2) {
        // vert_filt as 2 tap
        y = h;

        if (subpel_y_q4 != 8) {
            if (w <= 4) {
                prepare_coeffs_2tap_sse2(
                    filter_params_y, subpel_y_q4, coeffs_128);

                if (w == 2) {
                    __m128i s32[2], s_128[2];

                    s32[0] = _mm_cvtsi32_si128(*(int32_t *)im);

                    do {
                        s32[1] = _mm_cvtsi32_si128(*(int32_t *)(im + 2));
                        s_128[0] = _mm_unpacklo_epi32(s32[0], s32[1]);
                        s32[0] = _mm_cvtsi32_si128(*(int32_t *)(im + 2 * 2));
                        s_128[1] = _mm_unpacklo_epi32(s32[1], s32[0]);
                        const __m128i ss =
                            _mm_unpacklo_epi16(s_128[0], s_128[1]);
                        const __m128i res =
                            convolve16_2tap_sse2(&ss, coeffs_128);
                        const __m128i r = convolve_2dy_round_sse2(res);
                        const __m128i rr = _mm_packs_epi32(r, r);
                        convolve_store_2x2_sse2(rr, dst, dst_stride);
                        im += 2 * 2;
                        dst += 2 * dst_stride;
                        y -= 2;
                    } while (y);
                }
                else {
                    __m128i s_64[2], s_128[2];

                    assert(w == 4);

                    s_64[0] = _mm_loadl_epi64((__m128i *)im);

                    do {
                        s_64[1] = _mm_loadl_epi64((__m128i *)(im + 4));
                        s_128[0] = _mm_unpacklo_epi64(s_64[0], s_64[1]);
                        s_64[0] = _mm_loadl_epi64((__m128i *)(im + 2 * 4));
                        s_128[1] = _mm_unpacklo_epi64(s_64[1], s_64[0]);
                        const __m128i ss0 =
                            _mm_unpacklo_epi16(s_128[0], s_128[1]);
                        const __m128i ss1 =
                            _mm_unpackhi_epi16(s_128[0], s_128[1]);
                        const __m128i res0 =
                            convolve16_2tap_sse2(&ss0, coeffs_128);
                        const __m128i res1 =
                            convolve16_2tap_sse2(&ss1, coeffs_128);
                        const __m128i r0 = convolve_2dy_round_sse2(res0);
                        const __m128i r1 = convolve_2dy_round_sse2(res1);
                        const __m128i r = _mm_packs_epi32(r0, r1);
                        convolve_store_4x2_sse2(r, dst, dst_stride);
                        im += 2 * 4;
                        dst += 2 * dst_stride;
                        y -= 2;
                    } while (y);
                }
            }
            else {
                prepare_coeffs_2tap_avx2(
                    filter_params_y, subpel_y_q4, coeffs_256);

                if (w == 8) {
                    __m128i s_128[2];
                    __m256i s_256[2];

                    s_128[0] = _mm_load_si128((__m128i *)im);

                    do {
                        s_128[1] = _mm_load_si128((__m128i *)(im + 8));
                        s_256[0] = _mm256_setr_m128i(s_128[0], s_128[1]);
                        s_128[0] = _mm_load_si128((__m128i *)(im + 2 * 8));
                        s_256[1] = _mm256_setr_m128i(s_128[1], s_128[0]);
                        const __m256i ss0 =
                            _mm256_unpacklo_epi16(s_256[0], s_256[1]);
                        const __m256i ss1 =
                            _mm256_unpackhi_epi16(s_256[0], s_256[1]);
                        const __m256i res0 =
                            convolve16_2tap_avx2(&ss0, coeffs_256);
                        const __m256i res1 =
                            convolve16_2tap_avx2(&ss1, coeffs_256);
                        const __m256i r0 = convolve_2dy_round_avx2(res0);
                        const __m256i r1 = convolve_2dy_round_avx2(res1);
                        const __m256i r = _mm256_packs_epi32(r0, r1);
                        convolve_store_8x2_avx2(r, dst, dst_stride);
                        im += 2 * 8;
                        dst += 2 * dst_stride;
                        y -= 2;
                    } while (y);
                }
                else if (w == 16) {
                    __m256i s_256[2];

                    s_256[0] = _mm256_load_si256((__m256i *)im);

                    do {
                        s_256[1] = _mm256_load_si256((__m256i *)(im + 16));
                        const __m256i r0 = convolve_2dy_2tap_16_avx2(
                            s_256[0], s_256[1], coeffs_256);
                        s_256[0] = _mm256_load_si256((__m256i *)(im + 2 * 16));
                        const __m256i r1 = convolve_2dy_2tap_16_avx2(
                            s_256[1], s_256[0], coeffs_256);
                        convolve_2dy_store_16x2_avx2(r0, r1, dst, dst_stride);
                        im += 2 * 16;
                        dst += 2 * dst_stride;
                        y -= 2;
                    } while (y);
                }
                else if (w == 32) {
                    __m256i s_256[2][2];

                    s_256[0][0] = _mm256_load_si256((__m256i *)(im + 0 * 16));
                    s_256[0][1] = _mm256_load_si256((__m256i *)(im + 1 * 16));

                    do {
                        convolve_2dy_2tap_32_avx2(
                            im + 32, s_256[0], s_256[1], coeffs_256, dst);
                        convolve_2dy_2tap_32_avx2(im + 2 * 32,
                            s_256[1],
                            s_256[0],
                            coeffs_256,
                            dst + dst_stride);
                        im += 2 * 32;
                        dst += 2 * dst_stride;
                        y -= 2;
                    } while (y);
                }
                else if (w == 64) {
                    __m256i s_256[2][4];

                    s_256[0][0] = _mm256_load_si256((__m256i *)(im + 0 * 16));
                    s_256[0][1] = _mm256_load_si256((__m256i *)(im + 1 * 16));
                    s_256[0][2] = _mm256_load_si256((__m256i *)(im + 2 * 16));
                    s_256[0][3] = _mm256_load_si256((__m256i *)(im + 3 * 16));

                    do {
                        convolve_2dy_2tap_32_avx2(im + 64,
                            s_256[0] + 0,
                            s_256[1] + 0,
                            coeffs_256,
                            dst);
                        convolve_2dy_2tap_32_avx2(im + 96,
                            s_256[0] + 2,
                            s_256[1] + 2,
                            coeffs_256,
                            dst + 32);
                        im += 2 * 64;
                        convolve_2dy_2tap_32_avx2(im,
                            s_256[1] + 0,
                            s_256[0] + 0,
                            coeffs_256,
                            dst + dst_stride);
                        convolve_2dy_2tap_32_avx2(im + 32,
                            s_256[1] + 2,
                            s_256[0] + 2,
                            coeffs_256,
                            dst + dst_stride + 32);
                        dst += 2 * dst_stride;
                        y -= 2;
                    } while (y);
                }
                else {
                    __m256i s_256[2][8];

                    assert(w == 128);

                    load_16bit_8rows_avx2(im, 16, s_256[0]);

                    do {
                        convolve_2dy_2tap_32_avx2(im + 128,
                            s_256[0] + 0,
                            s_256[1] + 0,
                            coeffs_256,
                            dst);
                        convolve_2dy_2tap_32_avx2(im + 160,
                            s_256[0] + 2,
                            s_256[1] + 2,
                            coeffs_256,
                            dst + 1 * 32);
                        convolve_2dy_2tap_32_avx2(im + 192,
                            s_256[0] + 4,
                            s_256[1] + 4,
                            coeffs_256,
                            dst + 2 * 32);
                        convolve_2dy_2tap_32_avx2(im + 224,
                            s_256[0] + 6,
                            s_256[1] + 6,
                            coeffs_256,
                            dst + 3 * 32);
                        im += 2 * 128;
                        convolve_2dy_2tap_32_avx2(im,
                            s_256[1] + 0,
                            s_256[0] + 0,
                            coeffs_256,
                            dst + dst_stride);
                        convolve_2dy_2tap_32_avx2(im + 32,
                            s_256[1] + 2,
                            s_256[0] + 2,
                            coeffs_256,
                            dst + dst_stride + 1 * 32);
                        convolve_2dy_2tap_32_avx2(im + 64,
                            s_256[1] + 4,
                            s_256[0] + 4,
                            coeffs_256,
                            dst + dst_stride + 2 * 32);
                        convolve_2dy_2tap_32_avx2(im + 96,
                            s_256[1] + 6,
                            s_256[0] + 6,
                            coeffs_256,
                            dst + dst_stride + 3 * 32);
                        dst += 2 * dst_stride;
                        y -= 2;
                    } while (y);
                }
            }
        }
        else {
            if (w == 2) {
                __m128i s32[2], s_128[2];

                s32[0] = _mm_cvtsi32_si128(*(int32_t *)im);

                do {
                    s32[1] = _mm_cvtsi32_si128(*(int32_t *)(im + 2));
                    s_128[0] = _mm_unpacklo_epi32(s32[0], s32[1]);
                    s32[0] = _mm_cvtsi32_si128(*(int32_t *)(im + 2 * 2));
                    s_128[1] = _mm_unpacklo_epi32(s32[1], s32[0]);
                    const __m128i r =
                        convolve_2dy_avg_round_sse2(s_128[0], s_128[1]);
                    convolve_store_2x2_sse2(r, dst, dst_stride);
                    im += 2 * 2;
                    dst += 2 * dst_stride;
                    y -= 2;
                } while (y);
            }
            else if (w == 4) {
                __m128i s_64[2], s_128[2];

                s_64[0] = _mm_loadl_epi64((__m128i *)im);

                do {
                    s_64[1] = _mm_loadl_epi64((__m128i *)(im + 4));
                    s_128[0] = _mm_unpacklo_epi64(s_64[0], s_64[1]);
                    s_64[0] = _mm_loadl_epi64((__m128i *)(im + 2 * 4));
                    s_128[1] = _mm_unpacklo_epi64(s_64[1], s_64[0]);
                    const __m128i r =
                        convolve_2dy_avg_round_sse2(s_128[0], s_128[1]);
                    convolve_store_4x2_sse2(r, dst, dst_stride);
                    im += 2 * 4;
                    dst += 2 * dst_stride;
                    y -= 2;
                } while (y);
            }
            else if (w == 8) {
                __m128i s_128[2];
                __m256i s_256[2];

                s_128[0] = _mm_load_si128((__m128i *)im);

                do {
                    s_128[1] = _mm_load_si128((__m128i *)(im + 8));
                    s_256[0] = _mm256_setr_m128i(s_128[0], s_128[1]);
                    s_128[0] = _mm_load_si128((__m128i *)(im + 2 * 8));
                    s_256[1] = _mm256_setr_m128i(s_128[1], s_128[0]);
                    const __m256i r =
                        convolve_2dy_avg_round_avx2(s_256[0], s_256[1]);
                    convolve_store_8x2_avx2(r, dst, dst_stride);
                    im += 2 * 8;
                    dst += 2 * dst_stride;
                    y -= 2;
                } while (y);
            }
            else if (w == 16) {
                __m256i s_256[2];

                s_256[0] = _mm256_load_si256((__m256i *)im);

                do {
                    s_256[1] = _mm256_load_si256((__m256i *)(im + 16));
                    const __m256i r0 =
                        convolve_2dy_avg_round_avx2(s_256[0], s_256[1]);
                    s_256[0] = _mm256_load_si256((__m256i *)(im + 2 * 16));
                    const __m256i r1 =
                        convolve_2dy_avg_round_avx2(s_256[1], s_256[0]);
                    convolve_2dy_store_16x2_avx2(r0, r1, dst, dst_stride);
                    im += 2 * 16;
                    dst += 2 * dst_stride;
                    y -= 2;
                } while (y);
            }
            else if (w == 32) {
                __m256i s_256[2][2];

                s_256[0][0] = _mm256_load_si256((__m256i *)(im + 0 * 16));
                s_256[0][1] = _mm256_load_si256((__m256i *)(im + 1 * 16));

                do {
                    convolve_2dy_2tap_32_avg_avx2(
                        im + 32, s_256[0], s_256[1], dst);
                    convolve_2dy_2tap_32_avg_avx2(
                        im + 2 * 32, s_256[1], s_256[0], dst + dst_stride);
                    im += 2 * 32;
                    dst += 2 * dst_stride;
                    y -= 2;
                } while (y);
            }
            else if (w == 64) {
                __m256i s_256[2][4];

                s_256[0][0] = _mm256_load_si256((__m256i *)(im + 0 * 16));
                s_256[0][1] = _mm256_load_si256((__m256i *)(im + 1 * 16));
                s_256[0][2] = _mm256_load_si256((__m256i *)(im + 2 * 16));
                s_256[0][3] = _mm256_load_si256((__m256i *)(im + 3 * 16));

                do {
                    convolve_2dy_2tap_32_avg_avx2(
                        im + 64, s_256[0] + 0, s_256[1] + 0, dst);
                    convolve_2dy_2tap_32_avg_avx2(
                        im + 96, s_256[0] + 2, s_256[1] + 2, dst + 32);
                    im += 2 * 64;
                    convolve_2dy_2tap_32_avg_avx2(
                        im, s_256[1] + 0, s_256[0] + 0, dst + dst_stride);
                    convolve_2dy_2tap_32_avg_avx2(im + 32,
                        s_256[1] + 2,
                        s_256[0] + 2,
                        dst + dst_stride + 32);
                    dst += 2 * dst_stride;
                    y -= 2;
                } while (y);
            }
            else {
                __m256i s_256[2][8];

                assert(w == 128);

                load_16bit_8rows_avx2(im, 16, s_256[0]);

                do {
                    convolve_2dy_2tap_32_avg_avx2(
                        im + 128, s_256[0] + 0, s_256[1] + 0, dst);
                    convolve_2dy_2tap_32_avg_avx2(
                        im + 160, s_256[0] + 2, s_256[1] + 2, dst + 1 * 32);
                    convolve_2dy_2tap_32_avg_avx2(
                        im + 192, s_256[0] + 4, s_256[1] + 4, dst + 2 * 32);
                    convolve_2dy_2tap_32_avg_avx2(
                        im + 224, s_256[0] + 6, s_256[1] + 6, dst + 3 * 32);
                    im += 2 * 128;
                    convolve_2dy_2tap_32_avg_avx2(
                        im, s_256[1] + 0, s_256[0] + 0, dst + dst_stride);
                    convolve_2dy_2tap_32_avg_avx2(im + 32,
                        s_256[1] + 2,
                        s_256[0] + 2,
                        dst + dst_stride + 1 * 32);
                    convolve_2dy_2tap_32_avg_avx2(im + 64,
                        s_256[1] + 4,
                        s_256[0] + 4,
                        dst + dst_stride + 2 * 32);
                    convolve_2dy_2tap_32_avg_avx2(im + 96,
                        s_256[1] + 6,
                        s_256[0] + 6,
                        dst + dst_stride + 3 * 32);
                    dst += 2 * dst_stride;
                    y -= 2;
                } while (y);
            }
        }
    }
    else if (tap_y == 4) {
        // vert_filt as 4 tap
        y = h;

        if (w == 2) {
            __m128i s_32[4], ss_128[2];

            prepare_coeffs_4tap_sse2(filter_params_y, subpel_y_q4, coeffs_128);

            s_32[0] = _mm_cvtsi32_si128(*(int32_t *)(im + 0 * 2));
            s_32[1] = _mm_cvtsi32_si128(*(int32_t *)(im + 1 * 2));
            s_32[2] = _mm_cvtsi32_si128(*(int32_t *)(im + 2 * 2));

            const __m128i src01 = _mm_unpacklo_epi32(s_32[0], s_32[1]);
            const __m128i src12 = _mm_unpacklo_epi32(s_32[1], s_32[2]);

            ss_128[0] = _mm_unpacklo_epi16(src01, src12);

            do {
                s_32[3] = _mm_cvtsi32_si128(*(int32_t *)(im + 3 * 2));
                const __m128i src45 = _mm_unpacklo_epi32(s_32[2], s_32[3]);
                s_32[2] = _mm_cvtsi32_si128(*(int32_t *)(im + 4 * 2));
                const __m128i src56 = _mm_unpacklo_epi32(s_32[3], s_32[2]);
                ss_128[1] = _mm_unpacklo_epi16(src45, src56);

                const __m128i res = convolve16_4tap_sse2(ss_128, coeffs_128);
                const __m128i r = convolve_2dy_round_sse2(res);
                const __m128i rr = _mm_packs_epi32(r, r);
                convolve_store_2x2_sse2(rr, dst, dst_stride);

                ss_128[0] = ss_128[1];
                im += 2 * 2;
                dst += 2 * dst_stride;
                y -= 2;
            } while (y);
        }
        else {
            prepare_coeffs_4tap_avx2(filter_params_y, subpel_y_q4, coeffs_256);

            if (w == 4) {
                __m128i s_64[4];
                __m256i s_256[4], ss_256[2];

                s_64[0] = _mm_loadl_epi64((__m128i *)(im + 0 * 4));
                s_64[1] = _mm_loadl_epi64((__m128i *)(im + 1 * 4));
                s_64[2] = _mm_loadl_epi64((__m128i *)(im + 2 * 4));

                // Load lines a and b. Line a to lower 128, line b to upper 128
                s_256[0] = _mm256_setr_m128i(s_64[0], s_64[1]);
                s_256[1] = _mm256_setr_m128i(s_64[1], s_64[2]);

                ss_256[0] = _mm256_unpacklo_epi16(s_256[0], s_256[1]);

                do {
                    s_64[3] = _mm_loadl_epi64((__m128i *)(im + 3 * 4));
                    s_256[2] = _mm256_setr_m128i(s_64[2], s_64[3]);
                    s_64[2] = _mm_loadl_epi64((__m128i *)(im + 4 * 4));
                    s_256[3] = _mm256_setr_m128i(s_64[3], s_64[2]);
                    ss_256[1] = _mm256_unpacklo_epi16(s_256[2], s_256[3]);

                    const __m256i res =
                        convolve16_4tap_avx2(ss_256, coeffs_256);
                    const __m256i r = convolve_2dy_round_avx2(res);
                    const __m256i rr = _mm256_packs_epi32(r, r);
                    convolve_store_4x2_avx2(rr, dst, dst_stride);

                    ss_256[0] = ss_256[1];
                    im += 2 * 4;
                    dst += 2 * dst_stride;
                    y -= 2;
                } while (y);
            }
            else if (w == 8) {
                __m256i s_256[4];

                s_256[0] = _mm256_loadu_si256((__m256i *)(im + 0 * 8));
                s_256[1] = _mm256_loadu_si256((__m256i *)(im + 1 * 8));

                if (subpel_y_q4 != 8) {
                    __m256i ss_256[4];

                    ss_256[0] = _mm256_unpacklo_epi16(s_256[0], s_256[1]);
                    ss_256[2] = _mm256_unpackhi_epi16(s_256[0], s_256[1]);

                    do {
                        s_256[2] = _mm256_loadu_si256((__m256i *)(im + 2 * 8));
                        s_256[3] = _mm256_loadu_si256((__m256i *)(im + 3 * 8));
                        ss_256[1] = _mm256_unpacklo_epi16(s_256[2], s_256[3]);
                        ss_256[3] = _mm256_unpackhi_epi16(s_256[2], s_256[3]);

                        const __m256i r =
                            convolve_2dy_4tap_16_avx2(ss_256, coeffs_256);
                        convolve_store_8x2_avx2(r, dst, dst_stride);

                        ss_256[0] = ss_256[1];
                        ss_256[2] = ss_256[3];
                        im += 2 * 8;
                        dst += 2 * dst_stride;
                        y -= 2;
                    } while (y);
                }
                else {
                    __m256i a_256[2];

                    do {
                        s_256[2] = _mm256_loadu_si256((__m256i *)(im + 2 * 8));
                        s_256[3] = _mm256_loadu_si256((__m256i *)(im + 3 * 8));

                        a_256[0] = _mm256_add_epi16(s_256[0], s_256[3]);
                        a_256[1] = _mm256_add_epi16(s_256[1], s_256[2]);
                        const __m256i r = convolve_2dy_2tap_16_avx2(
                            a_256[0], a_256[1], coeffs_256);
                        convolve_store_8x2_avx2(r, dst, dst_stride);

                        s_256[0] = s_256[2];
                        s_256[1] = s_256[3];
                        im += 2 * 8;
                        dst += 2 * dst_stride;
                        y -= 2;
                    } while (y);
                }
            }
            else {
                __m256i s_256[5];

                assert(w == 16);

                s_256[0] = _mm256_loadu_si256((__m256i *)(im + 0 * 16));
                s_256[1] = _mm256_loadu_si256((__m256i *)(im + 1 * 16));
                s_256[2] = _mm256_loadu_si256((__m256i *)(im + 2 * 16));

                if (subpel_y_q4 != 8) {
                    __m256i ss_256[4], tt_256[4];

                    ss_256[0] = _mm256_unpacklo_epi16(s_256[0], s_256[1]);
                    ss_256[2] = _mm256_unpackhi_epi16(s_256[0], s_256[1]);

                    tt_256[0] = _mm256_unpacklo_epi16(s_256[1], s_256[2]);
                    tt_256[2] = _mm256_unpackhi_epi16(s_256[1], s_256[2]);

                    do {
                        s_256[3] = _mm256_loadu_si256((__m256i *)(im + 3 * 16));
                        ss_256[1] = _mm256_unpacklo_epi16(s_256[2], s_256[3]);
                        ss_256[3] = _mm256_unpackhi_epi16(s_256[2], s_256[3]);
                        s_256[2] = _mm256_loadu_si256((__m256i *)(im + 4 * 16));
                        tt_256[1] = _mm256_unpacklo_epi16(s_256[3], s_256[2]);
                        tt_256[3] = _mm256_unpackhi_epi16(s_256[3], s_256[2]);

                        const __m256i r0 =
                            convolve_2dy_4tap_16_avx2(ss_256, coeffs_256);
                        const __m256i r1 =
                            convolve_2dy_4tap_16_avx2(tt_256, coeffs_256);
                        convolve_2dy_store_16x2_avx2(r0, r1, dst, dst_stride);

                        ss_256[0] = ss_256[1];
                        ss_256[2] = ss_256[3];

                        tt_256[0] = tt_256[1];
                        tt_256[2] = tt_256[3];
                        im += 2 * 16;
                        dst += 2 * dst_stride;
                        y -= 2;
                    } while (y);
                }
                else {
                    __m256i a_256[2];

                    do {
                        s_256[3] = _mm256_loadu_si256((__m256i *)(im + 3 * 16));
                        s_256[4] = _mm256_loadu_si256((__m256i *)(im + 4 * 16));

                        a_256[0] = _mm256_add_epi16(s_256[0], s_256[3]);
                        a_256[1] = _mm256_add_epi16(s_256[1], s_256[2]);
                        const __m256i r0 = convolve_2dy_2tap_16_avx2(
                            a_256[0], a_256[1], coeffs_256);

                        a_256[0] = _mm256_add_epi16(s_256[1], s_256[4]);
                        a_256[1] = _mm256_add_epi16(s_256[2], s_256[3]);
                        const __m256i r1 = convolve_2dy_2tap_16_avx2(
                            a_256[0], a_256[1], coeffs_256);
                        convolve_2dy_store_16x2_avx2(r0, r1, dst, dst_stride);

                        s_256[0] = s_256[2];
                        s_256[1] = s_256[3];
                        s_256[2] = s_256[4];
                        im += 2 * 16;
                        dst += 2 * dst_stride;
                        y -= 2;
                    } while (y);
                }
            }
        }
    }
    else if (tap_y == 6) {
        // vert_filt as 6 tap
        if (w == 2) {
            __m128i s_32[6], ss_128[3];

            prepare_coeffs_6tap_ssse3(filter_params_y, subpel_y_q4, coeffs_128);

            s_32[0] = _mm_cvtsi32_si128(*(int32_t *)(im + 0 * 2));
            s_32[1] = _mm_cvtsi32_si128(*(int32_t *)(im + 1 * 2));
            s_32[2] = _mm_cvtsi32_si128(*(int32_t *)(im + 2 * 2));
            s_32[3] = _mm_cvtsi32_si128(*(int32_t *)(im + 3 * 2));
            s_32[4] = _mm_cvtsi32_si128(*(int32_t *)(im + 4 * 2));

            const __m128i src01 = _mm_unpacklo_epi32(s_32[0], s_32[1]);
            const __m128i src12 = _mm_unpacklo_epi32(s_32[1], s_32[2]);
            const __m128i src23 = _mm_unpacklo_epi32(s_32[2], s_32[3]);
            const __m128i src34 = _mm_unpacklo_epi32(s_32[3], s_32[4]);

            ss_128[0] = _mm_unpacklo_epi16(src01, src12);
            ss_128[1] = _mm_unpacklo_epi16(src23, src34);

            y = h;
            do {
                s_32[5] = _mm_cvtsi32_si128(*(int32_t *)(im + 5 * 2));
                const __m128i src45 = _mm_unpacklo_epi32(s_32[4], s_32[5]);
                s_32[4] = _mm_cvtsi32_si128(*(int32_t *)(im + 6 * 2));
                const __m128i src56 = _mm_unpacklo_epi32(s_32[5], s_32[4]);
                ss_128[2] = _mm_unpacklo_epi16(src45, src56);

                const __m128i res = convolve16_6tap_sse2(ss_128, coeffs_128);
                const __m128i r = convolve_2dy_round_sse2(res);
                const __m128i rr = _mm_packs_epi32(r, r);
                convolve_store_2x2_sse2(rr, dst, dst_stride);

                ss_128[0] = ss_128[1];
                ss_128[1] = ss_128[2];
                im += 2 * 2;
                dst += 2 * dst_stride;
                y -= 2;
            } while (y);
        }
        else {
            prepare_coeffs_6tap_avx2(filter_params_y, subpel_y_q4, coeffs_256);

            if (w == 4) {
                __m128i s_64[6];
                __m256i s_256[6], ss_256[3];

                s_64[0] = _mm_loadl_epi64((__m128i *)(im + 0 * 4));
                s_64[1] = _mm_loadl_epi64((__m128i *)(im + 1 * 4));
                s_64[2] = _mm_loadl_epi64((__m128i *)(im + 2 * 4));
                s_64[3] = _mm_loadl_epi64((__m128i *)(im + 3 * 4));
                s_64[4] = _mm_loadl_epi64((__m128i *)(im + 4 * 4));

                // Load lines a and b. Line a to lower 128, line b to upper 128
                s_256[0] = _mm256_setr_m128i(s_64[0], s_64[1]);
                s_256[1] = _mm256_setr_m128i(s_64[1], s_64[2]);
                s_256[2] = _mm256_setr_m128i(s_64[2], s_64[3]);
                s_256[3] = _mm256_setr_m128i(s_64[3], s_64[4]);

                ss_256[0] = _mm256_unpacklo_epi16(s_256[0], s_256[1]);
                ss_256[1] = _mm256_unpacklo_epi16(s_256[2], s_256[3]);

                y = h;
                do {
                    s_64[5] = _mm_loadl_epi64((__m128i *)(im + 5 * 4));
                    s_256[4] = _mm256_setr_m128i(s_64[4], s_64[5]);
                    s_64[4] = _mm_loadl_epi64((__m128i *)(im + 6 * 4));
                    s_256[5] = _mm256_setr_m128i(s_64[5], s_64[4]);
                    ss_256[2] = _mm256_unpacklo_epi16(s_256[4], s_256[5]);

                    const __m256i res =
                        convolve16_6tap_avx2(ss_256, coeffs_256);
                    const __m256i r = convolve_2dy_round_avx2(res);
                    const __m256i rr = _mm256_packs_epi32(r, r);
                    convolve_store_4x2_avx2(rr, dst, dst_stride);

                    ss_256[0] = ss_256[1];
                    ss_256[1] = ss_256[2];
                    im += 2 * 4;
                    dst += 2 * dst_stride;
                    y -= 2;
                } while (y);
            }
            else if (w == 8) {
                __m256i s_256[6];

                s_256[0] = _mm256_loadu_si256((__m256i *)(im + 0 * 8));
                s_256[1] = _mm256_loadu_si256((__m256i *)(im + 1 * 8));
                s_256[2] = _mm256_loadu_si256((__m256i *)(im + 2 * 8));
                s_256[3] = _mm256_loadu_si256((__m256i *)(im + 3 * 8));
                y = h;

                if (subpel_y_q4 != 8) {
                    __m256i ss_256[6];

                    ss_256[0] = _mm256_unpacklo_epi16(s_256[0], s_256[1]);
                    ss_256[1] = _mm256_unpacklo_epi16(s_256[2], s_256[3]);

                    ss_256[3] = _mm256_unpackhi_epi16(s_256[0], s_256[1]);
                    ss_256[4] = _mm256_unpackhi_epi16(s_256[2], s_256[3]);

                    do {
                        s_256[4] = _mm256_loadu_si256((__m256i *)(im + 4 * 8));
                        s_256[5] = _mm256_loadu_si256((__m256i *)(im + 5 * 8));
                        ss_256[2] = _mm256_unpacklo_epi16(s_256[4], s_256[5]);
                        ss_256[5] = _mm256_unpackhi_epi16(s_256[4], s_256[5]);

                        const __m256i r =
                            convolve_2dy_6tap_16_avx2(ss_256, coeffs_256);
                        convolve_store_8x2_avx2(r, dst, dst_stride);

                        ss_256[0] = ss_256[1];
                        ss_256[1] = ss_256[2];
                        ss_256[3] = ss_256[4];
                        ss_256[4] = ss_256[5];
                        im += 2 * 8;
                        dst += 2 * dst_stride;
                        y -= 2;
                    } while (y);
                }
                else {
                    __m256i a_256[4], ss_256[4];

                    do {
                        s_256[4] = _mm256_loadu_si256((__m256i *)(im + 4 * 8));
                        s_256[5] = _mm256_loadu_si256((__m256i *)(im + 5 * 8));

                        a_256[0] = _mm256_add_epi16(s_256[0], s_256[5]);
                        a_256[1] = _mm256_add_epi16(s_256[1], s_256[4]);
                        ss_256[0] = _mm256_unpacklo_epi16(a_256[0], a_256[1]);
                        ss_256[1] = _mm256_unpacklo_epi16(s_256[2], s_256[3]);
                        ss_256[2] = _mm256_unpackhi_epi16(a_256[0], a_256[1]);
                        ss_256[3] = _mm256_unpackhi_epi16(s_256[2], s_256[3]);

                        const __m256i r =
                            convolve_2dy_4tap_16_avx2(ss_256, coeffs_256);
                        convolve_store_8x2_avx2(r, dst, dst_stride);

                        s_256[0] = s_256[2];
                        s_256[1] = s_256[3];
                        s_256[2] = s_256[4];
                        s_256[3] = s_256[5];
                        im += 2 * 8;
                        dst += 2 * dst_stride;
                        y -= 2;
                    } while (y);
                }
            }
            else if (w == 16) {
                __m256i s_256[6];

                s_256[0] = _mm256_loadu_si256((__m256i *)(im + 0 * 16));
                s_256[1] = _mm256_loadu_si256((__m256i *)(im + 1 * 16));
                s_256[2] = _mm256_loadu_si256((__m256i *)(im + 2 * 16));
                s_256[3] = _mm256_loadu_si256((__m256i *)(im + 3 * 16));
                s_256[4] = _mm256_loadu_si256((__m256i *)(im + 4 * 16));
                y = h;

                if (subpel_y_q4 != 8) {
                    __m256i ss_256[6], tt_256[6];

                    ss_256[0] = _mm256_unpacklo_epi16(s_256[0], s_256[1]);
                    ss_256[1] = _mm256_unpacklo_epi16(s_256[2], s_256[3]);
                    ss_256[3] = _mm256_unpackhi_epi16(s_256[0], s_256[1]);
                    ss_256[4] = _mm256_unpackhi_epi16(s_256[2], s_256[3]);

                    tt_256[0] = _mm256_unpacklo_epi16(s_256[1], s_256[2]);
                    tt_256[1] = _mm256_unpacklo_epi16(s_256[3], s_256[4]);
                    tt_256[3] = _mm256_unpackhi_epi16(s_256[1], s_256[2]);
                    tt_256[4] = _mm256_unpackhi_epi16(s_256[3], s_256[4]);

                    do {
                        s_256[5] = _mm256_loadu_si256((__m256i *)(im + 5 * 16));
                        ss_256[2] = _mm256_unpacklo_epi16(s_256[4], s_256[5]);
                        ss_256[5] = _mm256_unpackhi_epi16(s_256[4], s_256[5]);
                        s_256[4] = _mm256_loadu_si256((__m256i *)(im + 6 * 16));
                        tt_256[2] = _mm256_unpacklo_epi16(s_256[5], s_256[4]);
                        tt_256[5] = _mm256_unpackhi_epi16(s_256[5], s_256[4]);

                        const __m256i r0 =
                            convolve_2dy_6tap_16_avx2(ss_256, coeffs_256);
                        const __m256i r1 =
                            convolve_2dy_6tap_16_avx2(tt_256, coeffs_256);
                        convolve_2dy_store_16x2_avx2(r0, r1, dst, dst_stride);

                        ss_256[0] = ss_256[1];
                        ss_256[1] = ss_256[2];
                        ss_256[3] = ss_256[4];
                        ss_256[4] = ss_256[5];

                        tt_256[0] = tt_256[1];
                        tt_256[1] = tt_256[2];
                        tt_256[3] = tt_256[4];
                        tt_256[4] = tt_256[5];
                        im += 2 * 16;
                        dst += 2 * dst_stride;
                        y -= 2;
                    } while (y);
                }
                else {
                    __m256i a_256[2], ss_256[4];

                    do {
                        s_256[5] = _mm256_loadu_si256((__m256i *)(im + 5 * 16));
                        a_256[0] = _mm256_add_epi16(s_256[0], s_256[5]);
                        a_256[1] = _mm256_add_epi16(s_256[1], s_256[4]);
                        ss_256[0] = _mm256_unpacklo_epi16(a_256[0], a_256[1]);
                        ss_256[1] = _mm256_unpacklo_epi16(s_256[2], s_256[3]);
                        ss_256[2] = _mm256_unpackhi_epi16(a_256[0], a_256[1]);
                        ss_256[3] = _mm256_unpackhi_epi16(s_256[2], s_256[3]);
                        const __m256i r0 =
                            convolve_2dy_4tap_16_avx2(ss_256, coeffs_256);

                        a_256[1] = _mm256_add_epi16(s_256[2], s_256[5]);
                        s_256[0] = s_256[2];
                        s_256[2] = s_256[4];
                        s_256[4] = _mm256_loadu_si256((__m256i *)(im + 6 * 16));
                        a_256[0] = _mm256_add_epi16(s_256[1], s_256[4]);
                        s_256[1] = s_256[3];
                        s_256[3] = s_256[5];
                        ss_256[0] = _mm256_unpacklo_epi16(a_256[0], a_256[1]);
                        ss_256[1] = _mm256_unpacklo_epi16(s_256[1], s_256[2]);
                        ss_256[2] = _mm256_unpackhi_epi16(a_256[0], a_256[1]);
                        ss_256[3] = _mm256_unpackhi_epi16(s_256[1], s_256[2]);
                        const __m256i r1 =
                            convolve_2dy_4tap_16_avx2(ss_256, coeffs_256);
                        convolve_2dy_store_16x2_avx2(r0, r1, dst, dst_stride);

                        im += 2 * 16;
                        dst += 2 * dst_stride;
                        y -= 2;
                    } while (y);
                }
            }
            else {
                int32_t x = 0;

                assert(!(w % 32));

                if (subpel_y_q4 != 8) {
                    __m256i s_256[2][6], ss_256[2][6], tt_256[2][6];

                    do {
                        const int16_t *s = im + x;
                        uint8_t *d = dst + x;

                        s_256[0][0] =
                            _mm256_loadu_si256((__m256i *)(s + 0 * w));
                        s_256[0][1] =
                            _mm256_loadu_si256((__m256i *)(s + 1 * w));
                        s_256[0][2] =
                            _mm256_loadu_si256((__m256i *)(s + 2 * w));
                        s_256[0][3] =
                            _mm256_loadu_si256((__m256i *)(s + 3 * w));
                        s_256[0][4] =
                            _mm256_loadu_si256((__m256i *)(s + 4 * w));

                        ss_256[0][0] =
                            _mm256_unpacklo_epi16(s_256[0][0], s_256[0][1]);
                        ss_256[0][1] =
                            _mm256_unpacklo_epi16(s_256[0][2], s_256[0][3]);
                        ss_256[0][3] =
                            _mm256_unpackhi_epi16(s_256[0][0], s_256[0][1]);
                        ss_256[0][4] =
                            _mm256_unpackhi_epi16(s_256[0][2], s_256[0][3]);

                        tt_256[0][0] =
                            _mm256_unpacklo_epi16(s_256[0][1], s_256[0][2]);
                        tt_256[0][1] =
                            _mm256_unpacklo_epi16(s_256[0][3], s_256[0][4]);
                        tt_256[0][3] =
                            _mm256_unpackhi_epi16(s_256[0][1], s_256[0][2]);
                        tt_256[0][4] =
                            _mm256_unpackhi_epi16(s_256[0][3], s_256[0][4]);

                        s_256[1][0] =
                            _mm256_loadu_si256((__m256i *)(s + 0 * w + 16));
                        s_256[1][1] =
                            _mm256_loadu_si256((__m256i *)(s + 1 * w + 16));
                        s_256[1][2] =
                            _mm256_loadu_si256((__m256i *)(s + 2 * w + 16));
                        s_256[1][3] =
                            _mm256_loadu_si256((__m256i *)(s + 3 * w + 16));
                        s_256[1][4] =
                            _mm256_loadu_si256((__m256i *)(s + 4 * w + 16));

                        ss_256[1][0] =
                            _mm256_unpacklo_epi16(s_256[1][0], s_256[1][1]);
                        ss_256[1][1] =
                            _mm256_unpacklo_epi16(s_256[1][2], s_256[1][3]);
                        ss_256[1][3] =
                            _mm256_unpackhi_epi16(s_256[1][0], s_256[1][1]);
                        ss_256[1][4] =
                            _mm256_unpackhi_epi16(s_256[1][2], s_256[1][3]);

                        tt_256[1][0] =
                            _mm256_unpacklo_epi16(s_256[1][1], s_256[1][2]);
                        tt_256[1][1] =
                            _mm256_unpacklo_epi16(s_256[1][3], s_256[1][4]);
                        tt_256[1][3] =
                            _mm256_unpackhi_epi16(s_256[1][1], s_256[1][2]);
                        tt_256[1][4] =
                            _mm256_unpackhi_epi16(s_256[1][3], s_256[1][4]);

                        y = h;
                        do {
                            s_256[0][5] =
                                _mm256_loadu_si256((__m256i *)(s + 5 * w));
                            ss_256[0][2] =
                                _mm256_unpacklo_epi16(s_256[0][4], s_256[0][5]);
                            ss_256[0][5] =
                                _mm256_unpackhi_epi16(s_256[0][4], s_256[0][5]);
                            s_256[0][4] =
                                _mm256_loadu_si256((__m256i *)(s + 6 * w));
                            tt_256[0][2] =
                                _mm256_unpacklo_epi16(s_256[0][5], s_256[0][4]);
                            tt_256[0][5] =
                                _mm256_unpackhi_epi16(s_256[0][5], s_256[0][4]);

                            const __m256i r0 = convolve_2dy_6tap_16_avx2(
                                ss_256[0], coeffs_256);
                            const __m256i r1 = convolve_2dy_6tap_16_avx2(
                                tt_256[0], coeffs_256);

                            s_256[1][5] =
                                _mm256_loadu_si256((__m256i *)(s + 5 * w + 16));
                            ss_256[1][2] =
                                _mm256_unpacklo_epi16(s_256[1][4], s_256[1][5]);
                            ss_256[1][5] =
                                _mm256_unpackhi_epi16(s_256[1][4], s_256[1][5]);
                            s_256[1][4] =
                                _mm256_loadu_si256((__m256i *)(s + 6 * w + 16));
                            tt_256[1][2] =
                                _mm256_unpacklo_epi16(s_256[1][5], s_256[1][4]);
                            tt_256[1][5] =
                                _mm256_unpackhi_epi16(s_256[1][5], s_256[1][4]);

                            const __m256i x0 = convolve_2dy_6tap_16_avx2(
                                ss_256[1], coeffs_256);
                            const __m256i x1 = convolve_2dy_6tap_16_avx2(
                                tt_256[1], coeffs_256);

                            convolve_2dy_store_32_avx2(r0, x0, d);
                            convolve_2dy_store_32_avx2(r1, x1, d + dst_stride);

                            ss_256[0][0] = ss_256[0][1];
                            ss_256[0][1] = ss_256[0][2];
                            ss_256[0][3] = ss_256[0][4];
                            ss_256[0][4] = ss_256[0][5];

                            tt_256[0][0] = tt_256[0][1];
                            tt_256[0][1] = tt_256[0][2];
                            tt_256[0][3] = tt_256[0][4];
                            tt_256[0][4] = tt_256[0][5];

                            ss_256[1][0] = ss_256[1][1];
                            ss_256[1][1] = ss_256[1][2];
                            ss_256[1][3] = ss_256[1][4];
                            ss_256[1][4] = ss_256[1][5];

                            tt_256[1][0] = tt_256[1][1];
                            tt_256[1][1] = tt_256[1][2];
                            tt_256[1][3] = tt_256[1][4];
                            tt_256[1][4] = tt_256[1][5];
                            s += 2 * w;
                            d += 2 * dst_stride;
                            y -= 2;
                        } while (y);

                        x += 32;
                    } while (x < w);
                }
                else {
                    __m256i s_256[2][6], a_256[2][2], ss_256[2][4];

                    do {
                        const int16_t *s = im + x;
                        uint8_t *d = dst + x;

                        s_256[0][0] =
                            _mm256_loadu_si256((__m256i *)(s + 0 * w));
                        s_256[0][1] =
                            _mm256_loadu_si256((__m256i *)(s + 1 * w));
                        s_256[0][2] =
                            _mm256_loadu_si256((__m256i *)(s + 2 * w));
                        s_256[0][3] =
                            _mm256_loadu_si256((__m256i *)(s + 3 * w));
                        s_256[0][4] =
                            _mm256_loadu_si256((__m256i *)(s + 4 * w));

                        s_256[1][0] =
                            _mm256_loadu_si256((__m256i *)(s + 0 * w + 16));
                        s_256[1][1] =
                            _mm256_loadu_si256((__m256i *)(s + 1 * w + 16));
                        s_256[1][2] =
                            _mm256_loadu_si256((__m256i *)(s + 2 * w + 16));
                        s_256[1][3] =
                            _mm256_loadu_si256((__m256i *)(s + 3 * w + 16));
                        s_256[1][4] =
                            _mm256_loadu_si256((__m256i *)(s + 4 * w + 16));

                        y = h;
                        do {
                            s_256[0][5] =
                                _mm256_loadu_si256((__m256i *)(s + 5 * w));
                            a_256[0][0] =
                                _mm256_add_epi16(s_256[0][0], s_256[0][5]);
                            a_256[0][1] =
                                _mm256_add_epi16(s_256[0][1], s_256[0][4]);
                            ss_256[0][0] =
                                _mm256_unpacklo_epi16(a_256[0][0], a_256[0][1]);
                            ss_256[0][1] =
                                _mm256_unpacklo_epi16(s_256[0][2], s_256[0][3]);
                            ss_256[0][2] =
                                _mm256_unpackhi_epi16(a_256[0][0], a_256[0][1]);
                            ss_256[0][3] =
                                _mm256_unpackhi_epi16(s_256[0][2], s_256[0][3]);
                            const __m256i r0 = convolve_2dy_4tap_16_avx2(
                                ss_256[0], coeffs_256);

                            a_256[0][1] =
                                _mm256_add_epi16(s_256[0][2], s_256[0][5]);
                            s_256[0][0] = s_256[0][2];
                            s_256[0][2] = s_256[0][4];
                            s_256[0][4] =
                                _mm256_loadu_si256((__m256i *)(s + 6 * w));
                            a_256[0][0] =
                                _mm256_add_epi16(s_256[0][1], s_256[0][4]);
                            s_256[0][1] = s_256[0][3];
                            s_256[0][3] = s_256[0][5];
                            ss_256[0][0] =
                                _mm256_unpacklo_epi16(a_256[0][0], a_256[0][1]);
                            ss_256[0][1] =
                                _mm256_unpacklo_epi16(s_256[0][1], s_256[0][2]);
                            ss_256[0][2] =
                                _mm256_unpackhi_epi16(a_256[0][0], a_256[0][1]);
                            ss_256[0][3] =
                                _mm256_unpackhi_epi16(s_256[0][1], s_256[0][2]);
                            const __m256i r1 = convolve_2dy_4tap_16_avx2(
                                ss_256[0], coeffs_256);

                            s_256[1][5] =
                                _mm256_loadu_si256((__m256i *)(s + 5 * w + 16));
                            a_256[1][0] =
                                _mm256_add_epi16(s_256[1][0], s_256[1][5]);
                            a_256[1][1] =
                                _mm256_add_epi16(s_256[1][1], s_256[1][4]);
                            ss_256[1][0] =
                                _mm256_unpacklo_epi16(a_256[1][0], a_256[1][1]);
                            ss_256[1][1] =
                                _mm256_unpacklo_epi16(s_256[1][2], s_256[1][3]);
                            ss_256[1][2] =
                                _mm256_unpackhi_epi16(a_256[1][0], a_256[1][1]);
                            ss_256[1][3] =
                                _mm256_unpackhi_epi16(s_256[1][2], s_256[1][3]);
                            const __m256i x0 = convolve_2dy_4tap_16_avx2(
                                ss_256[1], coeffs_256);

                            a_256[1][1] =
                                _mm256_add_epi16(s_256[1][2], s_256[1][5]);
                            s_256[1][0] = s_256[1][2];
                            s_256[1][2] = s_256[1][4];
                            s_256[1][4] =
                                _mm256_loadu_si256((__m256i *)(s + 6 * w + 16));
                            a_256[1][0] =
                                _mm256_add_epi16(s_256[1][1], s_256[1][4]);
                            s_256[1][1] = s_256[1][3];
                            s_256[1][3] = s_256[1][5];
                            ss_256[1][0] =
                                _mm256_unpacklo_epi16(a_256[1][0], a_256[1][1]);
                            ss_256[1][1] =
                                _mm256_unpacklo_epi16(s_256[1][1], s_256[1][2]);
                            ss_256[1][2] =
                                _mm256_unpackhi_epi16(a_256[1][0], a_256[1][1]);
                            ss_256[1][3] =
                                _mm256_unpackhi_epi16(s_256[1][1], s_256[1][2]);
                            const __m256i x1 = convolve_2dy_4tap_16_avx2(
                                ss_256[1], coeffs_256);

                            convolve_2dy_store_32_avx2(r0, x0, d);
                            convolve_2dy_store_32_avx2(r1, x1, d + dst_stride);

                            s += 2 * w;
                            d += 2 * dst_stride;
                            y -= 2;
                        } while (y);

                        x += 32;
                    } while (x < w);
                }
            }
        }
    }
    else {
        // vert_filt as 8 tap
        if (w == 2) {
            __m128i s_32[8], ss_128[4];

            prepare_coeffs_8tap_sse2(filter_params_y, subpel_y_q4, coeffs_128);

            s_32[0] = _mm_cvtsi32_si128(*(int32_t *)(im + 0 * 2));
            s_32[1] = _mm_cvtsi32_si128(*(int32_t *)(im + 1 * 2));
            s_32[2] = _mm_cvtsi32_si128(*(int32_t *)(im + 2 * 2));
            s_32[3] = _mm_cvtsi32_si128(*(int32_t *)(im + 3 * 2));
            s_32[4] = _mm_cvtsi32_si128(*(int32_t *)(im + 4 * 2));
            s_32[5] = _mm_cvtsi32_si128(*(int32_t *)(im + 5 * 2));
            s_32[6] = _mm_cvtsi32_si128(*(int32_t *)(im + 6 * 2));

            const __m128i src01 = _mm_unpacklo_epi32(s_32[0], s_32[1]);
            const __m128i src12 = _mm_unpacklo_epi32(s_32[1], s_32[2]);
            const __m128i src23 = _mm_unpacklo_epi32(s_32[2], s_32[3]);
            const __m128i src34 = _mm_unpacklo_epi32(s_32[3], s_32[4]);
            const __m128i src45 = _mm_unpacklo_epi32(s_32[4], s_32[5]);
            const __m128i src56 = _mm_unpacklo_epi32(s_32[5], s_32[6]);

            ss_128[0] = _mm_unpacklo_epi16(src01, src12);
            ss_128[1] = _mm_unpacklo_epi16(src23, src34);
            ss_128[2] = _mm_unpacklo_epi16(src45, src56);

            y = h;
            do {
                s_32[7] = _mm_cvtsi32_si128(*(int32_t *)(im + 7 * 2));
                const __m128i src67 = _mm_unpacklo_epi32(s_32[6], s_32[7]);
                s_32[6] = _mm_cvtsi32_si128(*(int32_t *)(im + 8 * 2));
                const __m128i src78 = _mm_unpacklo_epi32(s_32[7], s_32[6]);
                ss_128[3] = _mm_unpacklo_epi16(src67, src78);

                const __m128i res = convolve16_8tap_sse2(ss_128, coeffs_128);
                const __m128i r = convolve_2dy_round_sse2(res);
                const __m128i rr = _mm_packs_epi32(r, r);
                convolve_store_2x2_sse2(rr, dst, dst_stride);

                ss_128[0] = ss_128[1];
                ss_128[1] = ss_128[2];
                ss_128[2] = ss_128[3];
                im += 2 * 2;
                dst += 2 * dst_stride;
                y -= 2;
            } while (y);
        }
        else {
            prepare_coeffs_8tap_avx2(filter_params_y, subpel_y_q4, coeffs_256);

            if (w == 4) {
                __m128i s_64[8];
                __m256i s_256[8], ss_256[4];

                s_64[0] = _mm_loadl_epi64((__m128i *)(im + 0 * 4));
                s_64[1] = _mm_loadl_epi64((__m128i *)(im + 1 * 4));
                s_64[2] = _mm_loadl_epi64((__m128i *)(im + 2 * 4));
                s_64[3] = _mm_loadl_epi64((__m128i *)(im + 3 * 4));
                s_64[4] = _mm_loadl_epi64((__m128i *)(im + 4 * 4));
                s_64[5] = _mm_loadl_epi64((__m128i *)(im + 5 * 4));
                s_64[6] = _mm_loadl_epi64((__m128i *)(im + 6 * 4));

                // Load lines a and b. Line a to lower 128, line b to upper 128
                s_256[0] = _mm256_setr_m128i(s_64[0], s_64[1]);
                s_256[1] = _mm256_setr_m128i(s_64[1], s_64[2]);
                s_256[2] = _mm256_setr_m128i(s_64[2], s_64[3]);
                s_256[3] = _mm256_setr_m128i(s_64[3], s_64[4]);
                s_256[4] = _mm256_setr_m128i(s_64[4], s_64[5]);
                s_256[5] = _mm256_setr_m128i(s_64[5], s_64[6]);

                ss_256[0] = _mm256_unpacklo_epi16(s_256[0], s_256[1]);
                ss_256[1] = _mm256_unpacklo_epi16(s_256[2], s_256[3]);
                ss_256[2] = _mm256_unpacklo_epi16(s_256[4], s_256[5]);

                y = h;
                do {
                    s_64[7] = _mm_loadl_epi64((__m128i *)(im + 7 * 4));
                    s_256[6] = _mm256_setr_m128i(s_64[6], s_64[7]);
                    s_64[6] = _mm_loadl_epi64((__m128i *)(im + 8 * 4));
                    s_256[7] = _mm256_setr_m128i(s_64[7], s_64[6]);
                    ss_256[3] = _mm256_unpacklo_epi16(s_256[6], s_256[7]);

                    const __m256i res =
                        convolve16_8tap_avx2(ss_256, coeffs_256);
                    const __m256i r = convolve_2dy_round_avx2(res);
                    const __m256i rr = _mm256_packs_epi32(r, r);
                    convolve_store_4x2_avx2(rr, dst, dst_stride);

                    ss_256[0] = ss_256[1];
                    ss_256[1] = ss_256[2];
                    ss_256[2] = ss_256[3];
                    im += 2 * 4;
                    dst += 2 * dst_stride;
                    y -= 2;
                } while (y);
            }
            else if (w == 8) {
                __m256i s_256[8];

                s_256[0] = _mm256_loadu_si256((__m256i *)(im + 0 * 8));
                s_256[1] = _mm256_loadu_si256((__m256i *)(im + 1 * 8));
                s_256[2] = _mm256_loadu_si256((__m256i *)(im + 2 * 8));
                s_256[3] = _mm256_loadu_si256((__m256i *)(im + 3 * 8));
                s_256[4] = _mm256_loadu_si256((__m256i *)(im + 4 * 8));
                s_256[5] = _mm256_loadu_si256((__m256i *)(im + 5 * 8));
                y = h;

                if (subpel_y_q4 != 8) {
                    __m256i ss_256[8];

                    convolve_8tap_unapck_avx2(s_256, ss_256);

                    do {
                        s_256[6] = _mm256_loadu_si256((__m256i *)(im + 6 * 8));
                        s_256[7] = _mm256_loadu_si256((__m256i *)(im + 7 * 8));
                        ss_256[3] = _mm256_unpacklo_epi16(s_256[6], s_256[7]);
                        ss_256[7] = _mm256_unpackhi_epi16(s_256[6], s_256[7]);

                        const __m256i r =
                            convolve_2dy_8tap_16_avx2(ss_256, coeffs_256);
                        convolve_store_8x2_avx2(r, dst, dst_stride);

                        ss_256[0] = ss_256[1];
                        ss_256[1] = ss_256[2];
                        ss_256[2] = ss_256[3];
                        ss_256[4] = ss_256[5];
                        ss_256[5] = ss_256[6];
                        ss_256[6] = ss_256[7];
                        im += 2 * 8;
                        dst += 2 * dst_stride;
                        y -= 2;
                    } while (y);
                }
                else {
                    __m256i a_256[4], ss_256[4];

                    do {
                        s_256[6] = _mm256_loadu_si256((__m256i *)(im + 6 * 8));
                        s_256[7] = _mm256_loadu_si256((__m256i *)(im + 7 * 8));

                        a_256[0] = _mm256_add_epi16(s_256[0], s_256[7]);
                        a_256[1] = _mm256_add_epi16(s_256[1], s_256[6]);
                        a_256[2] = _mm256_add_epi16(s_256[2], s_256[5]);
                        a_256[3] = _mm256_add_epi16(s_256[3], s_256[4]);
                        ss_256[0] = _mm256_unpacklo_epi16(a_256[0], a_256[1]);
                        ss_256[1] = _mm256_unpacklo_epi16(a_256[2], a_256[3]);
                        ss_256[2] = _mm256_unpackhi_epi16(a_256[0], a_256[1]);
                        ss_256[3] = _mm256_unpackhi_epi16(a_256[2], a_256[3]);

                        const __m256i r =
                            convolve_2dy_4tap_16_avx2(ss_256, coeffs_256);
                        convolve_store_8x2_avx2(r, dst, dst_stride);

                        s_256[0] = s_256[2];
                        s_256[1] = s_256[3];
                        s_256[2] = s_256[4];
                        s_256[3] = s_256[5];
                        s_256[4] = s_256[6];
                        s_256[5] = s_256[7];
                        im += 2 * 8;
                        dst += 2 * dst_stride;
                        y -= 2;
                    } while (y);
                }
            }
            else if (w == 16) {
                __m256i s_256[8];

                load_16bit_7rows_avx2(im, 16, s_256);
                y = h;

                if (subpel_y_q4 != 8) {
                    __m256i ss_256[8], tt_256[8];

                    convolve_8tap_unapck_avx2(s_256, ss_256);
                    convolve_8tap_unapck_avx2(s_256 + 1, tt_256);

                    do {
                        __m256i r[2];

                        convolve_2dy_8tap_16x2_avx2(
                            im, 16, coeffs_256, s_256, ss_256, tt_256, r);
                        convolve_2dy_store_16x2_avx2(
                            r[0], r[1], dst, dst_stride);

                        im += 2 * 16;
                        dst += 2 * dst_stride;
                        y -= 2;
                    } while (y);
                }
                else {
                    do {
                        __m256i r[2];

                        convolve_2dy_8tap_16x2_half_pel_avx2(
                            im, 16, coeffs_256, s_256, r);
                        convolve_2dy_store_16x2_avx2(
                            r[0], r[1], dst, dst_stride);

                        im += 2 * 16;
                        dst += 2 * dst_stride;
                        y -= 2;
                    } while (y);
                }
            }
            else {
                int32_t x = 0;

                assert(!(w % 32));

                if (subpel_y_q4 != 8) {
                    __m256i s_256[2][8], ss_256[2][8], tt_256[2][8];

                    do {
                        const int16_t *s = im + x;
                        uint8_t *d = dst + x;

                        load_16bit_7rows_avx2(s, w, s_256[0]);
                        convolve_8tap_unapck_avx2(s_256[0], ss_256[0]);
                        convolve_8tap_unapck_avx2(s_256[0] + 1, tt_256[0]);

                        load_16bit_7rows_avx2(s + 16, w, s_256[1]);
                        convolve_8tap_unapck_avx2(s_256[1], ss_256[1]);
                        convolve_8tap_unapck_avx2(s_256[1] + 1, tt_256[1]);

                        y = h;
                        do {
                            __m256i r[2][2];

                            convolve_2dy_8tap_16x2_avx2(s,
                                w,
                                coeffs_256,
                                s_256[0],
                                ss_256[0],
                                tt_256[0],
                                r[0]);
                            convolve_2dy_8tap_16x2_avx2(s + 16,
                                w,
                                coeffs_256,
                                s_256[1],
                                ss_256[1],
                                tt_256[1],
                                r[1]);
                            convolve_2dy_store_32_avx2(r[0][0], r[1][0], d);
                            convolve_2dy_store_32_avx2(
                                r[0][1], r[1][1], d + dst_stride);

                            s += 2 * w;
                            d += 2 * dst_stride;
                            y -= 2;
                        } while (y);

                        x += 32;
                    } while (x < w);
                }
                else {
                    __m256i s_256[2][8];

                    do {
                        const int16_t *s = im + x;
                        uint8_t *d = dst + x;

                        load_16bit_7rows_avx2(s, w, s_256[0]);
                        load_16bit_7rows_avx2(s + 16, w, s_256[1]);

                        y = h;
                        do {
                            __m256i r[2][2];

                            convolve_2dy_8tap_16x2_half_pel_avx2(
                                s, w, coeffs_256, s_256[0], r[0]);
                            convolve_2dy_8tap_16x2_half_pel_avx2(
                                s + 16, w, coeffs_256, s_256[1], r[1]);
                            convolve_2dy_store_32_avx2(r[0][0], r[1][0], d);
                            convolve_2dy_store_32_avx2(
                                r[0][1], r[1][1], d + dst_stride);

                            s += 2 * w;
                            d += 2 * dst_stride;
                            y -= 2;
                        } while (y);

                        x += 32;
                    } while (x < w);
                }
            }
        }
    }
}

SIMD_INLINE void copy_128(const uint8_t *src, uint8_t *dst) {
    __m256i s[4];
    s[0] = _mm256_loadu_si256((__m256i *)(src + 0 * 32));
    s[1] = _mm256_loadu_si256((__m256i *)(src + 1 * 32));
    s[2] = _mm256_loadu_si256((__m256i *)(src + 2 * 32));
    s[3] = _mm256_loadu_si256((__m256i *)(src + 3 * 32));
    _mm256_storeu_si256((__m256i *)(dst + 0 * 32), s[0]);
    _mm256_storeu_si256((__m256i *)(dst + 1 * 32), s[1]);
    _mm256_storeu_si256((__m256i *)(dst + 2 * 32), s[2]);
    _mm256_storeu_si256((__m256i *)(dst + 3 * 32), s[3]);
}

void eb_av1_convolve_2d_copy_sr_avx2(const uint8_t *src, int32_t src_stride,
    uint8_t *dst, int32_t dst_stride, int32_t w, int32_t h,
    InterpFilterParams *filter_params_x,
    InterpFilterParams *filter_params_y,
    const int32_t subpel_x_q4, const int32_t subpel_y_q4,
    ConvolveParams *conv_params) {
    (void)filter_params_x;
    (void)filter_params_y;
    (void)subpel_x_q4;
    (void)subpel_y_q4;
    (void)conv_params;

    if (w == 2) {
        do {
            memcpy(dst, src, 2 * sizeof(*src));
            src += src_stride;
            dst += dst_stride;
            memcpy(dst, src, 2 * sizeof(*src));
            src += src_stride;
            dst += dst_stride;
            h -= 2;
        } while (h);
    }
    else if (w == 4) {
        do {
            memcpy(dst, src, 4 * sizeof(*src));
            src += src_stride;
            dst += dst_stride;
            memcpy(dst, src, 4 * sizeof(*src));
            src += src_stride;
            dst += dst_stride;
            h -= 2;
        } while (h);
    }
    else if (w == 8) {
        do {
            __m128i s[2];
            s[0] = _mm_loadl_epi64((__m128i *)src);
            src += src_stride;
            s[1] = _mm_loadl_epi64((__m128i *)src);
            src += src_stride;
            _mm_storel_epi64((__m128i *)dst, s[0]);
            dst += dst_stride;
            _mm_storel_epi64((__m128i *)dst, s[1]);
            dst += dst_stride;
            h -= 2;
        } while (h);
    }
    else if (w == 16) {
        do {
            __m128i s[2];
            s[0] = _mm_loadu_si128((__m128i *)src);
            src += src_stride;
            s[1] = _mm_loadu_si128((__m128i *)src);
            src += src_stride;
            _mm_storeu_si128((__m128i *)dst, s[0]);
            dst += dst_stride;
            _mm_storeu_si128((__m128i *)dst, s[1]);
            dst += dst_stride;
            h -= 2;
        } while (h);
    }
    else if (w == 32) {
        do {
            __m256i s[2];
            s[0] = _mm256_loadu_si256((__m256i *)src);
            src += src_stride;
            s[1] = _mm256_loadu_si256((__m256i *)src);
            src += src_stride;
            _mm256_storeu_si256((__m256i *)dst, s[0]);
            dst += dst_stride;
            _mm256_storeu_si256((__m256i *)dst, s[1]);
            dst += dst_stride;
            h -= 2;
        } while (h);
    }
    else if (w == 64) {
        do {
            __m256i s[4];
            s[0] = _mm256_loadu_si256((__m256i *)(src + 0 * 32));
            s[1] = _mm256_loadu_si256((__m256i *)(src + 1 * 32));
            src += src_stride;
            s[2] = _mm256_loadu_si256((__m256i *)(src + 0 * 32));
            s[3] = _mm256_loadu_si256((__m256i *)(src + 1 * 32));
            src += src_stride;
            _mm256_storeu_si256((__m256i *)(dst + 0 * 32), s[0]);
            _mm256_storeu_si256((__m256i *)(dst + 1 * 32), s[1]);
            dst += dst_stride;
            _mm256_storeu_si256((__m256i *)(dst + 0 * 32), s[2]);
            _mm256_storeu_si256((__m256i *)(dst + 1 * 32), s[3]);
            dst += dst_stride;
            h -= 2;
        } while (h);
    }
    else {
        do {
            copy_128(src, dst);
            src += src_stride;
            dst += dst_stride;
            copy_128(src, dst);
            src += src_stride;
            dst += dst_stride;
            h -= 2;
        } while (h);
    }
}
