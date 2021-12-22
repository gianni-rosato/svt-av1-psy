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

#include <emmintrin.h>
#include "EbDefinitions.h"
#include "common_dsp_rtcd.h"
#include "filter.h"

static INLINE void svt_prepare_coeffs(const InterpFilterParams *const filter_params,
                                      const int subpel_q4, __m128i *const coeffs /* [4] */) {
    const int16_t *const y_filter = av1_get_interp_filter_subpel_kernel(*filter_params,
                                                                        subpel_q4 & SUBPEL_MASK);
    const __m128i        coeffs_y = _mm_loadu_si128((__m128i *)y_filter);
    // coeffs 0 1 0 1 2 3 2 3
    const __m128i tmp_0 = _mm_unpacklo_epi32(coeffs_y, coeffs_y);
    // coeffs 4 5 4 5 6 7 6 7
    const __m128i tmp_1 = _mm_unpackhi_epi32(coeffs_y, coeffs_y);

    coeffs[0] = _mm_unpacklo_epi64(tmp_0, tmp_0); // coeffs 0 1 0 1 0 1 0 1
    coeffs[1] = _mm_unpackhi_epi64(tmp_0, tmp_0); // coeffs 2 3 2 3 2 3 2 3
    coeffs[2] = _mm_unpacklo_epi64(tmp_1, tmp_1); // coeffs 4 5 4 5 4 5 4 5
    coeffs[3] = _mm_unpackhi_epi64(tmp_1, tmp_1); // coeffs 6 7 6 7 6 7 6 7
}

static INLINE __m128i convolve(const __m128i *const s, const __m128i *const coeffs) {
    const __m128i d0 = _mm_madd_epi16(s[0], coeffs[0]);
    const __m128i d1 = _mm_madd_epi16(s[1], coeffs[1]);
    const __m128i d2 = _mm_madd_epi16(s[2], coeffs[2]);
    const __m128i d3 = _mm_madd_epi16(s[3], coeffs[3]);
    const __m128i d  = _mm_add_epi32(_mm_add_epi32(d0, d1), _mm_add_epi32(d2, d3));
    return d;
}

static INLINE __m128i convolve_lo_x(const __m128i *const s, const __m128i *const coeffs) {
    __m128i ss[4];
    ss[0] = _mm_unpacklo_epi8(s[0], _mm_setzero_si128());
    ss[1] = _mm_unpacklo_epi8(s[1], _mm_setzero_si128());
    ss[2] = _mm_unpacklo_epi8(s[2], _mm_setzero_si128());
    ss[3] = _mm_unpacklo_epi8(s[3], _mm_setzero_si128());
    return convolve(ss, coeffs);
}

static INLINE __m128i convolve_lo_y(const __m128i *const s, const __m128i *const coeffs) {
    __m128i ss[4];
    ss[0] = _mm_unpacklo_epi8(s[0], _mm_setzero_si128());
    ss[1] = _mm_unpacklo_epi8(s[2], _mm_setzero_si128());
    ss[2] = _mm_unpacklo_epi8(s[4], _mm_setzero_si128());
    ss[3] = _mm_unpacklo_epi8(s[6], _mm_setzero_si128());
    return convolve(ss, coeffs);
}

static INLINE __m128i convolve_hi_y(const __m128i *const s, const __m128i *const coeffs) {
    __m128i ss[4];
    ss[0] = _mm_unpackhi_epi8(s[0], _mm_setzero_si128());
    ss[1] = _mm_unpackhi_epi8(s[2], _mm_setzero_si128());
    ss[2] = _mm_unpackhi_epi8(s[4], _mm_setzero_si128());
    ss[3] = _mm_unpackhi_epi8(s[6], _mm_setzero_si128());
    return convolve(ss, coeffs);
}

static INLINE __m128i convolve_12tap(const __m128i *s, const __m128i *coeffs) {
    const __m128i d0     = _mm_madd_epi16(s[0], coeffs[0]);
    const __m128i d1     = _mm_madd_epi16(s[1], coeffs[1]);
    const __m128i d2     = _mm_madd_epi16(s[2], coeffs[2]);
    const __m128i d3     = _mm_madd_epi16(s[3], coeffs[3]);
    const __m128i d4     = _mm_madd_epi16(s[4], coeffs[4]);
    const __m128i d5     = _mm_madd_epi16(s[5], coeffs[5]);
    const __m128i d_0123 = _mm_add_epi32(_mm_add_epi32(d0, d1), _mm_add_epi32(d2, d3));
    const __m128i d      = _mm_add_epi32(_mm_add_epi32(d4, d5), d_0123);
    return d;
}

static INLINE __m128i convolve_lo_x_12tap(const __m128i *s, const __m128i *coeffs,
                                          const __m128i zero) {
    __m128i ss[6];
    ss[0] = _mm_unpacklo_epi8(s[0], zero); //  0  1  1  2  2  3  3  4
    ss[1] = _mm_unpacklo_epi8(s[1], zero); //  2  3  3  4  4  5  5  6
    ss[2] = _mm_unpacklo_epi8(s[2], zero); //  4  5  5  6  6  7  7  8
    ss[3] = _mm_unpacklo_epi8(s[3], zero); //  6  7  7  8  8  9  9 10
    ss[4] = _mm_unpackhi_epi8(s[2], zero); //  8  9  9 10 10 11 11 12
    ss[5] = _mm_unpackhi_epi8(s[3], zero); // 10 11 11 12 12 13 13 14
    return convolve_12tap(ss, coeffs);
}

static INLINE __m128i convolve_lo_y_12tap(const __m128i *s, const __m128i *coeffs) {
    __m128i       ss[6];
    const __m128i zero = _mm_setzero_si128();
    ss[0]              = _mm_unpacklo_epi8(s[0], zero);
    ss[1]              = _mm_unpacklo_epi8(s[2], zero);
    ss[2]              = _mm_unpacklo_epi8(s[4], zero);
    ss[3]              = _mm_unpacklo_epi8(s[6], zero);
    ss[4]              = _mm_unpacklo_epi8(s[8], zero);
    ss[5]              = _mm_unpacklo_epi8(s[10], zero);
    return convolve_12tap(ss, coeffs);
}

static INLINE __m128i convolve_hi_y_12tap(const __m128i *s, const __m128i *coeffs) {
    __m128i       ss[6];
    const __m128i zero = _mm_setzero_si128();
    ss[0]              = _mm_unpackhi_epi8(s[0], zero);
    ss[1]              = _mm_unpackhi_epi8(s[2], zero);
    ss[2]              = _mm_unpackhi_epi8(s[4], zero);
    ss[3]              = _mm_unpackhi_epi8(s[6], zero);
    ss[4]              = _mm_unpackhi_epi8(s[8], zero);
    ss[5]              = _mm_unpackhi_epi8(s[10], zero);
    return convolve_12tap(ss, coeffs);
}

static INLINE void prepare_coeffs_12tap(const InterpFilterParams *filter_params, int subpel_q4,
                                        __m128i *coeffs /* [6] */) {
    const int16_t *const y_filter = av1_get_interp_filter_subpel_kernel(*filter_params,
                                                                        subpel_q4 & SUBPEL_MASK);

    __m128i coeffs_y = _mm_loadu_si128((__m128i *)y_filter);

    coeffs[0] = _mm_shuffle_epi32(coeffs_y, 0); // coeffs 0 1 0 1 0 1 0 1
    coeffs[1] = _mm_shuffle_epi32(coeffs_y, 85); // coeffs 2 3 2 3 2 3 2 3
    coeffs[2] = _mm_shuffle_epi32(coeffs_y, 170); // coeffs 4 5 4 5 4 5 4 5
    coeffs[3] = _mm_shuffle_epi32(coeffs_y, 255); // coeffs 6 7 6 7 6 7 6 7

    coeffs_y = _mm_loadl_epi64((__m128i *)(y_filter + 8));

    coeffs[4] = _mm_shuffle_epi32(coeffs_y, 0); // coeffs 8 9 8 9 8 9 8 9
    coeffs[5] = _mm_shuffle_epi32(coeffs_y, 85); // coeffs 10 11 10 11 10 11 10 11
}

void svt_av1_convolve_y_sr_12tap_sse2(const uint8_t *src, int src_stride, uint8_t *dst,
                                      int dst_stride, int w, int h,
                                      const InterpFilterParams *filter_params_y, int subpel_y_qn) {
    const int      fo_vert     = filter_params_y->taps / 2 - 1;
    const uint8_t *src_ptr     = src - fo_vert * src_stride;
    const __m128i  round_const = _mm_set1_epi32((1 << FILTER_BITS) >> 1);
    const __m128i  round_shift = _mm_cvtsi32_si128(FILTER_BITS);
    __m128i        coeffs[6];

    prepare_coeffs_12tap(filter_params_y, subpel_y_qn, coeffs);

    int j = 0;
    do {
        __m128i        s[12], src10, res_lo, res_hi;
        __m128i        res_lo_round, res_hi_round, res16, res;
        const uint8_t *data = &src_ptr[j];

        src10 = _mm_loadl_epi64((__m128i *)(data + 10 * src_stride));
        s[0]  = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(data + 0 * src_stride)),
                                 _mm_loadl_epi64((__m128i *)(data + 1 * src_stride)));
        s[1]  = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(data + 1 * src_stride)),
                                 _mm_loadl_epi64((__m128i *)(data + 2 * src_stride)));
        s[2]  = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(data + 2 * src_stride)),
                                 _mm_loadl_epi64((__m128i *)(data + 3 * src_stride)));
        s[3]  = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(data + 3 * src_stride)),
                                 _mm_loadl_epi64((__m128i *)(data + 4 * src_stride)));
        s[4]  = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(data + 4 * src_stride)),
                                 _mm_loadl_epi64((__m128i *)(data + 5 * src_stride)));
        s[5]  = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(data + 5 * src_stride)),
                                 _mm_loadl_epi64((__m128i *)(data + 6 * src_stride)));
        s[6]  = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(data + 6 * src_stride)),
                                 _mm_loadl_epi64((__m128i *)(data + 7 * src_stride)));
        s[7]  = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(data + 7 * src_stride)),
                                 _mm_loadl_epi64((__m128i *)(data + 8 * src_stride)));
        s[8]  = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(data + 8 * src_stride)),
                                 _mm_loadl_epi64((__m128i *)(data + 9 * src_stride)));
        s[9]  = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(data + 9 * src_stride)), src10);

        int i = 0;
        do {
            data  = &src_ptr[i * src_stride + j];
            s[10] = _mm_unpacklo_epi8(src10, _mm_loadl_epi64((__m128i *)(data + 11 * src_stride)));
            src10 = _mm_loadl_epi64((__m128i *)(data + 12 * src_stride));
            s[11] = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(data + 11 * src_stride)), src10);

            res_lo = convolve_lo_y_12tap(s, coeffs); // Filter low index pixels
            res_hi = convolve_hi_y_12tap(s, coeffs); // Filter high index pixels

            res_lo_round = _mm_sra_epi32(_mm_add_epi32(res_lo, round_const), round_shift);
            res_hi_round = _mm_sra_epi32(_mm_add_epi32(res_hi, round_const), round_shift);

            res16 = _mm_packs_epi32(res_lo_round, res_hi_round);
            res   = _mm_packus_epi16(res16, res16);

            _mm_storel_epi64((__m128i *)(dst + i * dst_stride + j), res);
            i++;

            res_lo = convolve_lo_y_12tap(s + 1, coeffs); // Filter low index pixels
            res_hi = convolve_hi_y_12tap(s + 1, coeffs); // Filter high index pixels

            res_lo_round = _mm_sra_epi32(_mm_add_epi32(res_lo, round_const), round_shift);
            res_hi_round = _mm_sra_epi32(_mm_add_epi32(res_hi, round_const), round_shift);

            res16 = _mm_packs_epi32(res_lo_round, res_hi_round);
            res   = _mm_packus_epi16(res16, res16);

            _mm_storel_epi64((__m128i *)(dst + i * dst_stride + j), res);
            i++;

            s[0] = s[2];
            s[1] = s[3];
            s[2] = s[4];
            s[3] = s[5];
            s[4] = s[6];
            s[5] = s[7];
            s[6] = s[8];
            s[7] = s[9];
            s[8] = s[10];
            s[9] = s[11];
        } while (i < h);
        j += 8;
    } while (j < w);
}

void svt_av1_convolve_y_sr_sse2(const uint8_t *src, int32_t src_stride, uint8_t *dst,
                                int32_t dst_stride, int32_t w, int32_t h,
                                InterpFilterParams *filter_params_x,
                                InterpFilterParams *filter_params_y, const int32_t subpel_x_q4,
                                const int32_t subpel_y_q4, ConvolveParams *conv_params) {
    if (filter_params_y->taps > 8) {
        if (w < 8) {
            svt_av1_convolve_y_sr_c(src,
                                    src_stride,
                                    dst,
                                    dst_stride,
                                    w,
                                    h,
                                    filter_params_x,
                                    filter_params_y,
                                    subpel_x_q4,
                                    subpel_y_q4,
                                    conv_params);
        } else {
            svt_av1_convolve_y_sr_12tap_sse2(
                src, src_stride, dst, dst_stride, w, h, filter_params_y, subpel_y_q4);
        }
    } else {
        const int      fo_vert     = filter_params_y->taps / 2 - 1;
        const uint8_t *src_ptr     = src - fo_vert * src_stride;
        const __m128i  round_const = _mm_set1_epi32((1 << FILTER_BITS) >> 1);
        const __m128i  round_shift = _mm_cvtsi32_si128(FILTER_BITS);
        __m128i        coeffs[4];

        svt_prepare_coeffs(filter_params_y, subpel_y_q4, coeffs);

        if (w <= 4) {
            __m128i s[8], src6, res, res_round, res16;
            src6 = _mm_cvtsi32_si128(*(uint32_t *)(src_ptr + 6 * src_stride));
            s[0] = _mm_unpacklo_epi8(_mm_cvtsi32_si128(*(uint32_t *)(src_ptr + 0 * src_stride)),
                                     _mm_cvtsi32_si128(*(uint32_t *)(src_ptr + 1 * src_stride)));
            s[1] = _mm_unpacklo_epi8(_mm_cvtsi32_si128(*(uint32_t *)(src_ptr + 1 * src_stride)),
                                     _mm_cvtsi32_si128(*(uint32_t *)(src_ptr + 2 * src_stride)));
            s[2] = _mm_unpacklo_epi8(_mm_cvtsi32_si128(*(uint32_t *)(src_ptr + 2 * src_stride)),
                                     _mm_cvtsi32_si128(*(uint32_t *)(src_ptr + 3 * src_stride)));
            s[3] = _mm_unpacklo_epi8(_mm_cvtsi32_si128(*(uint32_t *)(src_ptr + 3 * src_stride)),
                                     _mm_cvtsi32_si128(*(uint32_t *)(src_ptr + 4 * src_stride)));
            s[4] = _mm_unpacklo_epi8(_mm_cvtsi32_si128(*(uint32_t *)(src_ptr + 4 * src_stride)),
                                     _mm_cvtsi32_si128(*(uint32_t *)(src_ptr + 5 * src_stride)));
            s[5] = _mm_unpacklo_epi8(_mm_cvtsi32_si128(*(uint32_t *)(src_ptr + 5 * src_stride)),
                                     src6);

            do {
                s[6] = _mm_unpacklo_epi8(
                    src6, _mm_cvtsi32_si128(*(uint32_t *)(src_ptr + 7 * src_stride)));
                src6 = _mm_cvtsi32_si128(*(uint32_t *)(src_ptr + 8 * src_stride));
                s[7] = _mm_unpacklo_epi8(_mm_cvtsi32_si128(*(uint32_t *)(src_ptr + 7 * src_stride)),
                                         src6);

                res              = convolve_lo_y(s + 0, coeffs);
                res_round        = _mm_sra_epi32(_mm_add_epi32(res, round_const), round_shift);
                res16            = _mm_packs_epi32(res_round, res_round);
                uint32_t res_int = _mm_cvtsi128_si32(_mm_packus_epi16(res16, res16));

                if (w == 2)
                    *(uint16_t *)dst = (uint16_t)res_int;
                else
                    *(uint32_t *)dst = res_int;

                src_ptr += src_stride;
                dst += dst_stride;

                res       = convolve_lo_y(s + 1, coeffs);
                res_round = _mm_sra_epi32(_mm_add_epi32(res, round_const), round_shift);
                res16     = _mm_packs_epi32(res_round, res_round);
                res_int   = _mm_cvtsi128_si32(_mm_packus_epi16(res16, res16));

                if (w == 2)
                    *(uint16_t *)dst = (uint16_t)res_int;
                else
                    *(uint32_t *)dst = res_int;

                src_ptr += src_stride;
                dst += dst_stride;

                s[0] = s[2];
                s[1] = s[3];
                s[2] = s[4];
                s[3] = s[5];
                s[4] = s[6];
                s[5] = s[7];
                h -= 2;
            } while (h);
        } else {
            assert(!(w % 8));
            int j = 0;
            do {
                __m128i        s[8], src6, res_lo, res_hi;
                __m128i        res_lo_round, res_hi_round, res16, res;
                const uint8_t *data = &src_ptr[j];

                src6 = _mm_loadl_epi64((__m128i *)(data + 6 * src_stride));
                s[0] = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(data + 0 * src_stride)),
                                         _mm_loadl_epi64((__m128i *)(data + 1 * src_stride)));
                s[1] = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(data + 1 * src_stride)),
                                         _mm_loadl_epi64((__m128i *)(data + 2 * src_stride)));
                s[2] = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(data + 2 * src_stride)),
                                         _mm_loadl_epi64((__m128i *)(data + 3 * src_stride)));
                s[3] = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(data + 3 * src_stride)),
                                         _mm_loadl_epi64((__m128i *)(data + 4 * src_stride)));
                s[4] = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(data + 4 * src_stride)),
                                         _mm_loadl_epi64((__m128i *)(data + 5 * src_stride)));
                s[5] = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(data + 5 * src_stride)), src6);

                int i = 0;
                do {
                    data = &src_ptr[i * src_stride + j];
                    s[6] = _mm_unpacklo_epi8(src6,
                                             _mm_loadl_epi64((__m128i *)(data + 7 * src_stride)));
                    src6 = _mm_loadl_epi64((__m128i *)(data + 8 * src_stride));
                    s[7] = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(data + 7 * src_stride)),
                                             src6);

                    res_lo = convolve_lo_y(s, coeffs); // Filter low index pixels
                    res_hi = convolve_hi_y(s, coeffs); // Filter high index pixels

                    res_lo_round = _mm_sra_epi32(_mm_add_epi32(res_lo, round_const), round_shift);
                    res_hi_round = _mm_sra_epi32(_mm_add_epi32(res_hi, round_const), round_shift);

                    res16 = _mm_packs_epi32(res_lo_round, res_hi_round);
                    res   = _mm_packus_epi16(res16, res16);

                    _mm_storel_epi64((__m128i *)(dst + i * dst_stride + j), res);
                    i++;

                    res_lo = convolve_lo_y(s + 1, coeffs); // Filter low index pixels
                    res_hi = convolve_hi_y(s + 1, coeffs); // Filter high index pixels

                    res_lo_round = _mm_sra_epi32(_mm_add_epi32(res_lo, round_const), round_shift);
                    res_hi_round = _mm_sra_epi32(_mm_add_epi32(res_hi, round_const), round_shift);

                    res16 = _mm_packs_epi32(res_lo_round, res_hi_round);
                    res   = _mm_packus_epi16(res16, res16);

                    _mm_storel_epi64((__m128i *)(dst + i * dst_stride + j), res);
                    i++;

                    s[0] = s[2];
                    s[1] = s[3];
                    s[2] = s[4];
                    s[3] = s[5];
                    s[4] = s[6];
                    s[5] = s[7];
                } while (i < h);
                j += 8;
            } while (j < w);
        }
    }
}

void svt_av1_convolve_x_sr_12tap_sse2(const uint8_t *src, int src_stride, uint8_t *dst,
                                      int dst_stride, int w, int h,
                                      const InterpFilterParams *filter_params_x, int subpel_x_qn,
                                      ConvolveParams *conv_params) {
    const int      fo_horiz      = filter_params_x->taps / 2 - 1;
    const uint8_t *src_ptr       = src - fo_horiz;
    const int      bits          = FILTER_BITS - conv_params->round_0;
    const __m128i  round_0_const = _mm_set1_epi32((1 << conv_params->round_0) >> 1);
    const __m128i  round_const   = _mm_set1_epi32((1 << bits) >> 1);
    const __m128i  round_0_shift = _mm_cvtsi32_si128(conv_params->round_0);
    const __m128i  round_shift   = _mm_cvtsi32_si128(bits);
    const __m128i  zero          = _mm_setzero_si128();
    __m128i        coeffs[6];

    assert(bits >= 0);
    assert((FILTER_BITS - conv_params->round_1) >= 0 ||
           ((conv_params->round_0 + conv_params->round_1) == 2 * FILTER_BITS));

    prepare_coeffs_12tap(filter_params_x, subpel_x_qn, coeffs);

    int i = 0;
    do {
        int j = 0;
        do {
            const __m128i data = _mm_loadu_si128((__m128i *)&src_ptr[i * src_stride + j]);
            __m128i       s[4];

            s[0] = _mm_unpacklo_epi16(data, _mm_srli_si128(data, 1));
            s[1] = _mm_unpacklo_epi16(_mm_srli_si128(data, 2), _mm_srli_si128(data, 3));
            s[2] = _mm_unpacklo_epi16(_mm_srli_si128(data, 4), _mm_srli_si128(data, 5));
            s[3] = _mm_unpacklo_epi16(_mm_srli_si128(data, 6), _mm_srli_si128(data, 7));

            const __m128i res32 = convolve_lo_x_12tap(s, coeffs, zero);

            __m128i res32_round = _mm_sra_epi32(_mm_add_epi32(res32, round_0_const), round_0_shift);
            res32_round = _mm_sra_epi32(_mm_add_epi32(res32_round, round_const), round_shift);

            const __m128i res16 = _mm_packs_epi32(res32_round, zero);
            const __m128i res   = _mm_packus_epi16(res16, zero);

            const int val = _mm_cvtsi128_si32(res);
            memcpy((dst + i * dst_stride + j), &val, sizeof(val));
            j += 4;
        } while (j < w);
    } while (++i < h);
}

void svt_av1_convolve_x_sr_sse2(const uint8_t *src, int32_t src_stride, uint8_t *dst,
                                int32_t dst_stride, int32_t w, int32_t h,
                                InterpFilterParams *filter_params_x,
                                InterpFilterParams *filter_params_y, const int32_t subpel_x_q4,
                                const int32_t subpel_y_q4, ConvolveParams *conv_params) {
    if (filter_params_x->taps > 8) {
        if (w < 4) {
            svt_av1_convolve_x_sr_c(src,
                                    src_stride,
                                    dst,
                                    dst_stride,
                                    w,
                                    h,
                                    filter_params_x,
                                    filter_params_y,
                                    subpel_x_q4,
                                    subpel_y_q4,
                                    conv_params);
        } else {
            svt_av1_convolve_x_sr_12tap_sse2(
                src, src_stride, dst, dst_stride, w, h, filter_params_x, subpel_x_q4, conv_params);
        }
    } else {
        const int      fo_horiz      = filter_params_x->taps / 2 - 1;
        const uint8_t *src_ptr       = src - fo_horiz;
        const int      bits          = FILTER_BITS - conv_params->round_0;
        const __m128i  round_0_const = _mm_set1_epi32((1 << conv_params->round_0) >> 1);
        const __m128i  round_const   = _mm_set1_epi32((1 << bits) >> 1);
        const __m128i  round_0_shift = _mm_cvtsi32_si128(conv_params->round_0);
        const __m128i  round_shift   = _mm_cvtsi32_si128(bits);
        __m128i        coeffs[4];

        assert(bits >= 0);
        assert((FILTER_BITS - conv_params->round_1) >= 0 ||
               ((conv_params->round_0 + conv_params->round_1) == 2 * FILTER_BITS));

        svt_prepare_coeffs(filter_params_x, subpel_x_q4, coeffs);

        if (w <= 4) {
            do {
                const __m128i data = _mm_loadu_si128((__m128i *)src_ptr);
                __m128i       s[4];

                s[0] = _mm_unpacklo_epi8(data, _mm_srli_si128(data, 1));
                s[1] = _mm_unpacklo_epi8(_mm_srli_si128(data, 2), _mm_srli_si128(data, 3));
                s[2] = _mm_unpacklo_epi8(_mm_srli_si128(data, 4), _mm_srli_si128(data, 5));
                s[3] = _mm_unpacklo_epi8(_mm_srli_si128(data, 6), _mm_srli_si128(data, 7));
                const __m128i res_lo       = convolve_lo_x(s, coeffs);
                __m128i       res_lo_round = _mm_sra_epi32(_mm_add_epi32(res_lo, round_0_const),
                                                     round_0_shift);
                res_lo_round = _mm_sra_epi32(_mm_add_epi32(res_lo_round, round_const), round_shift);

                const __m128i res16 = _mm_packs_epi32(res_lo_round, res_lo_round);
                const __m128i res   = _mm_packus_epi16(res16, res16);

                uint32_t r = _mm_cvtsi128_si32(res);
                if (w == 2)
                    *(uint16_t *)dst = (uint16_t)r;
                else
                    *(uint32_t *)dst = r;

                src_ptr += src_stride;
                dst += dst_stride;
            } while (--h);
        } else {
            assert(!(w % 8));
            int i = 0;
            do {
                int j = 0;
                do {
                    const __m128i data = _mm_loadu_si128((__m128i *)&src_ptr[i * src_stride + j]);
                    __m128i       s[4];

                    // Filter even-index pixels
                    s[0]                   = data;
                    s[1]                   = _mm_srli_si128(data, 2);
                    s[2]                   = _mm_srli_si128(data, 4);
                    s[3]                   = _mm_srli_si128(data, 6);
                    const __m128i res_even = convolve_lo_x(s, coeffs);

                    // Filter odd-index pixels
                    s[0]                  = _mm_srli_si128(data, 1);
                    s[1]                  = _mm_srli_si128(data, 3);
                    s[2]                  = _mm_srli_si128(data, 5);
                    s[3]                  = _mm_srli_si128(data, 7);
                    const __m128i res_odd = convolve_lo_x(s, coeffs);

                    // Rearrange pixels back into the order 0 ... 7
                    const __m128i res_lo       = _mm_unpacklo_epi32(res_even, res_odd);
                    const __m128i res_hi       = _mm_unpackhi_epi32(res_even, res_odd);
                    __m128i       res_lo_round = _mm_sra_epi32(_mm_add_epi32(res_lo, round_0_const),
                                                         round_0_shift);
                    res_lo_round         = _mm_sra_epi32(_mm_add_epi32(res_lo_round, round_const),
                                                 round_shift);
                    __m128i res_hi_round = _mm_sra_epi32(_mm_add_epi32(res_hi, round_0_const),
                                                         round_0_shift);
                    res_hi_round         = _mm_sra_epi32(_mm_add_epi32(res_hi_round, round_const),
                                                 round_shift);

                    const __m128i res16 = _mm_packs_epi32(res_lo_round, res_hi_round);
                    const __m128i res   = _mm_packus_epi16(res16, res16);

                    _mm_storel_epi64((__m128i *)(dst + i * dst_stride + j), res);
                    j += 8;
                } while (j < w);
            } while (++i < h);
        }
    }
}
