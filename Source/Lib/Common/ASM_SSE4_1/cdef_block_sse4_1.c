/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */
#include "EbDefinitions.h"
#include <smmintrin.h>
#include "EbCdef.h"
#include "EbBitstreamUnit.h"

typedef struct {
    __m128i val[2];
} v256;

SIMD_INLINE __m128i v128_shr_u8(__m128i a, unsigned int c) {
    return _mm_and_si128(_mm_set1_epi8((char)(0xff >> c)), _mm_srl_epi16(a, _mm_cvtsi32_si128(c)));
}

SIMD_INLINE v256 v256_from_v128(__m128i hi, __m128i lo) {
    v256 t;
    t.val[1] = hi;
    t.val[0] = lo;
    return t;
}

SIMD_INLINE v256 v256_dup_16(uint16_t x) {
    __m128i t = _mm_set1_epi16(x);
    return v256_from_v128(t, t);
}

SIMD_INLINE v256 v256_zero(void) {
    return v256_from_v128(_mm_setzero_si128(), _mm_setzero_si128());
}

SIMD_INLINE v256 v256_max_s16(v256 a, v256 b) {
    return v256_from_v128(_mm_max_epi16(a.val[1], b.val[1]), _mm_max_epi16(a.val[0], b.val[0]));
}

SIMD_INLINE v256 v256_andn(v256 a, v256 b) {
    return v256_from_v128(_mm_andnot_si128(b.val[1], a.val[1]), _mm_andnot_si128(b.val[0], a.val[0]));
}

SIMD_INLINE v256 v256_cmpeq_16(v256 a, v256 b) {
    return v256_from_v128(_mm_cmpeq_epi16(a.val[1], b.val[1]), _mm_cmpeq_epi16(a.val[0], b.val[0]));
}

SIMD_INLINE v256 v256_min_s16(v256 a, v256 b) {
    return v256_from_v128(_mm_min_epi16(a.val[1], b.val[1]), _mm_min_epi16(a.val[0], b.val[0]));
}

SIMD_INLINE v256 v256_sub_16(v256 a, v256 b) {
    return v256_from_v128(_mm_sub_epi16(a.val[1], b.val[1]), _mm_sub_epi16(a.val[0], b.val[0]));
}

SIMD_INLINE v256 v256_add_16(v256 a, v256 b) {
    return v256_from_v128(_mm_add_epi16(a.val[1], b.val[1]), _mm_add_epi16(a.val[0], b.val[0]));
}

SIMD_INLINE v256 v256_madd_us8(v256 a, v256 b) {
    return v256_from_v128(_mm_maddubs_epi16(a.val[1], b.val[1]), _mm_maddubs_epi16(a.val[0], b.val[0]));
}

SIMD_INLINE v256 v256_dup_8(uint8_t x) {
    __m128i t = _mm_set1_epi8(x);
    return v256_from_v128(t, t);
}

SIMD_INLINE v256 v256_cmplt_s16(v256 a, v256 b) {
    return v256_from_v128(_mm_cmplt_epi16(a.val[1], b.val[1]), _mm_cmplt_epi16(a.val[0], b.val[0]));
}

SIMD_INLINE v256 v256_pack_s16_u8(v256 a, v256 b) {
    return v256_from_v128(_mm_packus_epi16(a.val[0], a.val[1]),
                          _mm_packus_epi16(b.val[0], b.val[1]));
}

SIMD_INLINE v256 v256_from_v64(__m128i a, __m128i b, __m128i c, __m128i d) {
    return v256_from_v128(_mm_unpacklo_epi64(b, a), _mm_unpacklo_epi64(d, c));
}

SIMD_INLINE v256 v256_mullo_s16(v256 a, v256 b) {
    return v256_from_v128(_mm_mullo_epi16(a.val[1], b.val[1]), _mm_mullo_epi16(a.val[0], b.val[0]));
}

SIMD_INLINE v256 v256_abs_s16(v256 a) {
    return v256_from_v128(_mm_abs_epi16(a.val[1]), _mm_abs_epi16(a.val[0]));
}

SIMD_INLINE v256 v256_ssub_u16(v256 a, v256 b) {
    return v256_from_v128(_mm_subs_epu16(a.val[1], b.val[1]), _mm_subs_epu16(a.val[0], b.val[0]));
}

SIMD_INLINE v256 v256_shr_u16(v256 a, const unsigned int c) {
    return v256_from_v128(_mm_srli_epi16(a.val[1], c), _mm_srli_epi16(a.val[0], c));
}

SIMD_INLINE v256 v256_xor(v256 a, v256 b) {
    return v256_from_v128(_mm_xor_si128(a.val[1], b.val[1]), _mm_xor_si128(a.val[0], b.val[0]));
}

#define v128_shr_n_s16(a, c) _mm_srai_epi16(a, c)
#define v256_shr_n_s16(a, n) \
    v256_from_v128(v128_shr_n_s16(a.val[1], n), v128_shr_n_s16(a.val[0], n))

// sign(a - b) * min(abs(a - b), max(0, strength - (abs(a - b) >> adjdamp)))
SIMD_INLINE __m128i constrain(v256 a, v256 b, unsigned int strength, unsigned int adjdamp) {
    const v256 diff16 = v256_sub_16(a, b);
    __m128i    diff   = _mm_packs_epi16(diff16.val[0], diff16.val[1]);
    const __m128i sign   = _mm_cmplt_epi8(diff, _mm_setzero_si128());
    diff              = _mm_abs_epi8(diff);
    return _mm_xor_si128(
        _mm_add_epi8(
            sign,
            _mm_min_epu8(diff, _mm_subs_epu8(_mm_set1_epi8(strength), v128_shr_u8(diff, adjdamp)))),
        sign);
}

SIMD_INLINE v256 constrain16(v256 a, v256 b, unsigned int threshold, unsigned int adjdamp) {
    v256       diff = v256_sub_16(a, b);
    const v256 sign = v256_shr_n_s16(diff, 15);
    diff            = v256_abs_s16(diff);
    const v256 s    = v256_ssub_u16(v256_dup_16(threshold), v256_shr_u16(diff, adjdamp));
    return v256_xor(v256_add_16(sign, v256_min_s16(diff, s)), sign);
}

void svt_av1_cdef_filter_block_8xn_8_sse4_1(uint8_t *dst, int dstride, const uint16_t *in,
                                        int pri_strength, int sec_strength, int dir,
                                            int pri_damping, int sec_damping, int coeff_shift,
                                            uint8_t height, uint8_t subsampling_factor) {
    int  i;
    __m128i p0, p1, p2, p3;
    v256 sum, row, res, tap;
    v256 max, min, large = v256_dup_16(CDEF_VERY_LARGE);
    int  po1  = eb_cdef_directions[dir][0];
    int  po2  = eb_cdef_directions[dir][1];
    int  s1o1 = eb_cdef_directions[(dir + 2) & 7][0];
    int  s1o2 = eb_cdef_directions[(dir + 2) & 7][1];
    int  s2o1 = eb_cdef_directions[(dir + 6) & 7][0];
    int  s2o2 = eb_cdef_directions[(dir + 6) & 7][1];

    const int *pri_taps = eb_cdef_pri_taps[(pri_strength >> coeff_shift) & 1];
    const int *sec_taps = eb_cdef_sec_taps[0];

    if (pri_strength)
        pri_damping = AOMMAX(0, pri_damping - get_msb(pri_strength));
    if (sec_strength)
        sec_damping = AOMMAX(0, sec_damping - get_msb(sec_strength));
    for (i = 0; i < height; i += (2 * subsampling_factor)) {
        sum = v256_zero();
        row = v256_from_v128(_mm_loadu_si128((__m128i *)(in + i * CDEF_BSTRIDE)),
                             _mm_loadu_si128((__m128i *)(in + (i + subsampling_factor) * CDEF_BSTRIDE)));

        max = min = row;
        // Primary near taps
        tap = v256_from_v128(_mm_loadu_si128((__m128i *)(in+i * CDEF_BSTRIDE + po1)),
                             _mm_loadu_si128((__m128i *)(in+(i + subsampling_factor) * CDEF_BSTRIDE + po1)));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p0  = constrain(tap, row, pri_strength, pri_damping);
        tap = v256_from_v128(_mm_loadu_si128((__m128i *)(in+i * CDEF_BSTRIDE - po1)),
                             _mm_loadu_si128((__m128i *)(in+(i + subsampling_factor) * CDEF_BSTRIDE - po1)));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p1  = constrain(tap, row, pri_strength, pri_damping);

        // sum += pri_taps[0] * (p0 + p1)
        sum = v256_add_16(
            sum,
            v256_madd_us8(v256_dup_8(pri_taps[0]),
                          v256_from_v128(_mm_unpackhi_epi8(p1, p0), _mm_unpacklo_epi8(p1, p0))));

        // Primary far taps
        tap = v256_from_v128(_mm_loadu_si128((__m128i *)(in+i * CDEF_BSTRIDE + po2)),
                             _mm_loadu_si128((__m128i *)(in+(i + subsampling_factor) * CDEF_BSTRIDE + po2)));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p0  = constrain(tap, row, pri_strength, pri_damping);
        tap = v256_from_v128(_mm_loadu_si128((__m128i *)(in+i * CDEF_BSTRIDE - po2)),
                             _mm_loadu_si128((__m128i *)(in+(i + subsampling_factor) * CDEF_BSTRIDE - po2)));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p1  = constrain(tap, row, pri_strength, pri_damping);

        // sum += pri_taps[1] * (p0 + p1)
        sum = v256_add_16(
            sum,
            v256_madd_us8(v256_dup_8(pri_taps[1]),
                          v256_from_v128(_mm_unpackhi_epi8(p1, p0), _mm_unpacklo_epi8(p1, p0))));

        // Secondary near taps
        tap = v256_from_v128(_mm_loadu_si128((__m128i *)(in+i * CDEF_BSTRIDE + s1o1)),
                             _mm_loadu_si128((__m128i *)(in+(i + subsampling_factor) * CDEF_BSTRIDE + s1o1)));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p0  = constrain(tap, row, sec_strength, sec_damping);
        tap = v256_from_v128(_mm_loadu_si128((__m128i *)(in+i * CDEF_BSTRIDE - s1o1)),
                             _mm_loadu_si128((__m128i *)(in+(i + subsampling_factor) * CDEF_BSTRIDE - s1o1)));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p1  = constrain(tap, row, sec_strength, sec_damping);
        tap = v256_from_v128(_mm_loadu_si128((__m128i *)(in+i * CDEF_BSTRIDE + s2o1)),
                             _mm_loadu_si128((__m128i *)(in+(i + subsampling_factor) * CDEF_BSTRIDE + s2o1)));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p2  = constrain(tap, row, sec_strength, sec_damping);
        tap = v256_from_v128(_mm_loadu_si128((__m128i *)(in+i * CDEF_BSTRIDE - s2o1)),
                             _mm_loadu_si128((__m128i *)(in+(i + subsampling_factor) * CDEF_BSTRIDE - s2o1)));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p3  = constrain(tap, row, sec_strength, sec_damping);

        // sum += sec_taps[0] * (p0 + p1 + p2 + p3)
        p0  = _mm_add_epi8(p0, p1);
        p2  = _mm_add_epi8(p2, p3);
        sum = v256_add_16(
            sum,
            v256_madd_us8(v256_dup_8(sec_taps[0]),
                          v256_from_v128(_mm_unpackhi_epi8(p2, p0), _mm_unpacklo_epi8(p2, p0))));

        // Secondary far taps
        tap = v256_from_v128(_mm_loadu_si128((__m128i *)(in+i * CDEF_BSTRIDE + s1o2)),
                             _mm_loadu_si128((__m128i *)(in+(i + subsampling_factor) * CDEF_BSTRIDE + s1o2)));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p0  = constrain(tap, row, sec_strength, sec_damping);
        tap = v256_from_v128(_mm_loadu_si128((__m128i *)(in+i * CDEF_BSTRIDE - s1o2)),
                             _mm_loadu_si128((__m128i *)(in+(i + subsampling_factor) * CDEF_BSTRIDE - s1o2)));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p1  = constrain(tap, row, sec_strength, sec_damping);
        tap = v256_from_v128(_mm_loadu_si128((__m128i *)(in+i * CDEF_BSTRIDE + s2o2)),
                             _mm_loadu_si128((__m128i *)(in+(i + subsampling_factor) * CDEF_BSTRIDE + s2o2)));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p2  = constrain(tap, row, sec_strength, sec_damping);
        tap = v256_from_v128(_mm_loadu_si128((__m128i *)(in+i * CDEF_BSTRIDE - s2o2)),
                             _mm_loadu_si128((__m128i *)(in+(i + subsampling_factor) * CDEF_BSTRIDE - s2o2)));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p3  = constrain(tap, row, sec_strength, sec_damping);

        // sum += sec_taps[1] * (p0 + p1 + p2 + p3)
        p0  = _mm_add_epi8(p0, p1);
        p2  = _mm_add_epi8(p2, p3);
        sum = v256_add_16(
            sum,
            v256_madd_us8(v256_dup_8(sec_taps[1]),
                          v256_from_v128(_mm_unpackhi_epi8(p2, p0), _mm_unpacklo_epi8(p2, p0))));

        // res = row + ((sum - (sum < 0) + 8) >> 4)
        sum = v256_add_16(sum, v256_cmplt_s16(sum, v256_zero()));
        res = v256_add_16(sum, v256_dup_16(8));
        res = v256_shr_n_s16(res, 4);
        res = v256_add_16(row, res);
        res = v256_min_s16(v256_max_s16(res, min), max);
        res = v256_pack_s16_u8(res, res);

        p0 = res.val[0];
        _mm_storel_epi64((__m128i *)(dst + i * dstride), _mm_srli_si128(p0, 8));
        _mm_storel_epi64((__m128i *)(dst + (i + subsampling_factor) * dstride), p0);
    }
}

void svt_av1_cdef_filter_block_4xn_8_sse4_1(uint8_t *dst, int dstride, const uint16_t *in,
                                        int pri_strength, int sec_strength, int dir,
                                            int pri_damping, int sec_damping, int coeff_shift,
                                            uint8_t height, uint8_t subsampling_factor) {
    __m128i p0, p1, p2, p3;
    v256 sum, row, tap, res;
    v256 max, min, large = v256_dup_16(CDEF_VERY_LARGE);
    int  po1  = eb_cdef_directions[dir][0];
    int  po2  = eb_cdef_directions[dir][1];
    int  s1o1 = eb_cdef_directions[(dir + 2) & 7][0];
    int  s1o2 = eb_cdef_directions[(dir + 2) & 7][1];
    int  s2o1 = eb_cdef_directions[(dir + 6) & 7][0];
    int  s2o2 = eb_cdef_directions[(dir + 6) & 7][1];

    const int *pri_taps = eb_cdef_pri_taps[(pri_strength >> coeff_shift) & 1];
    const int *sec_taps = eb_cdef_sec_taps[(pri_strength >> coeff_shift) & 1];

    if (pri_strength)
        pri_damping = AOMMAX(0, pri_damping - get_msb(pri_strength));
    if (sec_strength)
        sec_damping = AOMMAX(0, sec_damping - get_msb(sec_strength));

    for (uint32_t i = 0; i < height; i += (4 * subsampling_factor)) {
        sum = v256_zero();
        row = v256_from_v64(_mm_loadl_epi64((__m128i *)(in +i * CDEF_BSTRIDE)),
                            _mm_loadl_epi64((__m128i *)(in +(i + (1 * subsampling_factor)) * CDEF_BSTRIDE)),
                            _mm_loadl_epi64((__m128i *)(in +(i + (2 * subsampling_factor)) * CDEF_BSTRIDE)),
                            _mm_loadl_epi64((__m128i *)(in +(i + (3 * subsampling_factor)) * CDEF_BSTRIDE)));
        max = min = row;

        if (pri_strength) {
            // Primary near taps
            tap = v256_from_v64(_mm_loadl_epi64((__m128i *)(in +i * CDEF_BSTRIDE + po1)),
                                _mm_loadl_epi64((__m128i *)(in +(i + (1 * subsampling_factor)) * CDEF_BSTRIDE + po1)),
                                _mm_loadl_epi64((__m128i *)(in +(i + (2 * subsampling_factor)) * CDEF_BSTRIDE + po1)),
                                _mm_loadl_epi64((__m128i *)(in +(i + (3 * subsampling_factor)) * CDEF_BSTRIDE + po1)));
            max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
            min = v256_min_s16(min, tap);
            p0  = constrain(tap, row, pri_strength, pri_damping);
            tap = v256_from_v64(_mm_loadl_epi64((__m128i *)(in +i * CDEF_BSTRIDE - po1)),
                                _mm_loadl_epi64((__m128i *)(in +(i + (1 * subsampling_factor)) * CDEF_BSTRIDE - po1)),
                                _mm_loadl_epi64((__m128i *)(in +(i + (2 * subsampling_factor)) * CDEF_BSTRIDE - po1)),
                                _mm_loadl_epi64((__m128i *)(in +(i + (3 * subsampling_factor)) * CDEF_BSTRIDE - po1)));
            max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
            min = v256_min_s16(min, tap);
            p1  = constrain(tap, row, pri_strength, pri_damping);

            // sum += pri_taps[0] * (p0 + p1)
            sum = v256_add_16(
                sum,
                v256_madd_us8(v256_dup_8(pri_taps[0]),
                              v256_from_v128(_mm_unpackhi_epi8(p1, p0), _mm_unpacklo_epi8(p1, p0))));

            // Primary far taps
            tap = v256_from_v64(_mm_loadl_epi64((__m128i *)(in +i * CDEF_BSTRIDE + po2)),
                                _mm_loadl_epi64((__m128i *)(in +(i + (1 * subsampling_factor)) * CDEF_BSTRIDE + po2)),
                                _mm_loadl_epi64((__m128i *)(in +(i + (2 * subsampling_factor)) * CDEF_BSTRIDE + po2)),
                                _mm_loadl_epi64((__m128i *)(in +(i + (3 * subsampling_factor)) * CDEF_BSTRIDE + po2)));
            max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
            min = v256_min_s16(min, tap);
            p0  = constrain(tap, row, pri_strength, pri_damping);
            tap = v256_from_v64(_mm_loadl_epi64((__m128i *)(in +i * CDEF_BSTRIDE - po2)),
                                _mm_loadl_epi64((__m128i *)(in +(i + (1 * subsampling_factor)) * CDEF_BSTRIDE - po2)),
                                _mm_loadl_epi64((__m128i *)(in +(i + (2 * subsampling_factor)) * CDEF_BSTRIDE - po2)),
                                _mm_loadl_epi64((__m128i *)(in +(i + (3 * subsampling_factor)) * CDEF_BSTRIDE - po2)));
            max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
            min = v256_min_s16(min, tap);
            p1  = constrain(tap, row, pri_strength, pri_damping);

            // sum += pri_taps[1] * (p0 + p1)
            sum = v256_add_16(
                sum,
                v256_madd_us8(v256_dup_8(pri_taps[1]),
                              v256_from_v128(_mm_unpackhi_epi8(p1, p0), _mm_unpacklo_epi8(p1, p0))));
        }

        if (sec_strength) {
            // Secondary near taps
            tap = v256_from_v64(_mm_loadl_epi64((__m128i *)(in + i * CDEF_BSTRIDE + s1o1)),
                                _mm_loadl_epi64((__m128i *)(in + (i + (1 * subsampling_factor)) * CDEF_BSTRIDE + s1o1)),
                                _mm_loadl_epi64((__m128i *)(in + (i + (2 * subsampling_factor)) * CDEF_BSTRIDE + s1o1)),
                                _mm_loadl_epi64((__m128i *)(in + (i + (3 * subsampling_factor)) * CDEF_BSTRIDE + s1o1)));
            max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
            min = v256_min_s16(min, tap);
            p0  = constrain(tap, row, sec_strength, sec_damping);
            tap = v256_from_v64(_mm_loadl_epi64((__m128i *)(in + i * CDEF_BSTRIDE - s1o1)),
                                _mm_loadl_epi64((__m128i *)(in + (i + (1 * subsampling_factor)) * CDEF_BSTRIDE - s1o1)),
                                _mm_loadl_epi64((__m128i *)(in + (i + (2 * subsampling_factor)) * CDEF_BSTRIDE - s1o1)),
                                _mm_loadl_epi64((__m128i *)(in + (i + (3 * subsampling_factor)) * CDEF_BSTRIDE - s1o1)));
            max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
            min = v256_min_s16(min, tap);
            p1  = constrain(tap, row, sec_strength, sec_damping);
            tap = v256_from_v64(_mm_loadl_epi64((__m128i *)(in + i * CDEF_BSTRIDE + s2o1)),
                                _mm_loadl_epi64((__m128i *)(in + (i + (1 * subsampling_factor)) * CDEF_BSTRIDE + s2o1)),
                                _mm_loadl_epi64((__m128i *)(in + (i + (2 * subsampling_factor)) * CDEF_BSTRIDE + s2o1)),
                                _mm_loadl_epi64((__m128i *)(in + (i + (3 * subsampling_factor)) * CDEF_BSTRIDE + s2o1)));
            max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
            min = v256_min_s16(min, tap);
            p2  = constrain(tap, row, sec_strength, sec_damping);
            tap = v256_from_v64(_mm_loadl_epi64((__m128i *)(in + i * CDEF_BSTRIDE - s2o1)),
                                _mm_loadl_epi64((__m128i *)(in + (i + (1 * subsampling_factor)) * CDEF_BSTRIDE - s2o1)),
                                _mm_loadl_epi64((__m128i *)(in + (i + (2 * subsampling_factor)) * CDEF_BSTRIDE - s2o1)),
                                _mm_loadl_epi64((__m128i *)(in + (i + (3 * subsampling_factor)) * CDEF_BSTRIDE - s2o1)));
            max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
            min = v256_min_s16(min, tap);
            p3  = constrain(tap, row, sec_strength, sec_damping);

            // sum += sec_taps[0] * (p0 + p1 + p2 + p3)
            p0  = _mm_add_epi8(p0, p1);
            p2  = _mm_add_epi8(p2, p3);
            sum = v256_add_16(
                sum,
                v256_madd_us8(v256_dup_8(sec_taps[0]),
                              v256_from_v128(_mm_unpackhi_epi8(p2, p0), _mm_unpacklo_epi8(p2, p0))));

            // Secondary far taps
            tap = v256_from_v64(_mm_loadl_epi64((__m128i *)(in +i * CDEF_BSTRIDE + s1o2)),
                                _mm_loadl_epi64((__m128i *)(in +(i + (1 * subsampling_factor)) * CDEF_BSTRIDE + s1o2)),
                                _mm_loadl_epi64((__m128i *)(in +(i + (2 * subsampling_factor)) * CDEF_BSTRIDE + s1o2)),
                                _mm_loadl_epi64((__m128i *)(in +(i + (3 * subsampling_factor)) * CDEF_BSTRIDE + s1o2)));
            max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
            min = v256_min_s16(min, tap);
            p0  = constrain(tap, row, sec_strength, sec_damping);
            tap = v256_from_v64(_mm_loadl_epi64((__m128i *)(in +i * CDEF_BSTRIDE - s1o2)),
                                _mm_loadl_epi64((__m128i *)(in +(i + (1 * subsampling_factor)) * CDEF_BSTRIDE - s1o2)),
                                _mm_loadl_epi64((__m128i *)(in +(i + (2 * subsampling_factor)) * CDEF_BSTRIDE - s1o2)),
                                _mm_loadl_epi64((__m128i *)(in +(i + (3 * subsampling_factor)) * CDEF_BSTRIDE - s1o2)));
            max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
            min = v256_min_s16(min, tap);
            p1  = constrain(tap, row, sec_strength, sec_damping);
            tap = v256_from_v64(_mm_loadl_epi64((__m128i *)(in +i * CDEF_BSTRIDE + s2o2)),
                                _mm_loadl_epi64((__m128i *)(in +(i + (1 * subsampling_factor)) * CDEF_BSTRIDE + s2o2)),
                                _mm_loadl_epi64((__m128i *)(in +(i + (2 * subsampling_factor)) * CDEF_BSTRIDE + s2o2)),
                                _mm_loadl_epi64((__m128i *)(in +(i + (3 * subsampling_factor)) * CDEF_BSTRIDE + s2o2)));
            max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
            min = v256_min_s16(min, tap);
            p2  = constrain(tap, row, sec_strength, sec_damping);
            tap = v256_from_v64(_mm_loadl_epi64((__m128i *)(in +i * CDEF_BSTRIDE - s2o2)),
                                _mm_loadl_epi64((__m128i *)(in +(i + (1 * subsampling_factor)) * CDEF_BSTRIDE - s2o2)),
                                _mm_loadl_epi64((__m128i *)(in +(i + (2 * subsampling_factor)) * CDEF_BSTRIDE - s2o2)),
                                _mm_loadl_epi64((__m128i *)(in +(i + (3 * subsampling_factor)) * CDEF_BSTRIDE - s2o2)));
            max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
            min = v256_min_s16(min, tap);
            p3  = constrain(tap, row, sec_strength, sec_damping);

            // sum += sec_taps[1] * (p0 + p1 + p2 + p3)
            p0 = _mm_add_epi8(p0, p1);
            p2 = _mm_add_epi8(p2, p3);

            sum = v256_add_16(
                sum,
                v256_madd_us8(v256_dup_8(sec_taps[1]),
                              v256_from_v128(_mm_unpackhi_epi8(p2, p0), _mm_unpacklo_epi8(p2, p0))));
        }

        // res = row + ((sum - (sum < 0) + 8) >> 4)
        sum = v256_add_16(sum, v256_cmplt_s16(sum, v256_zero()));
        res = v256_add_16(sum, v256_dup_16(8));
        res = v256_shr_n_s16(res, 4);
        res = v256_add_16(row, res);
        res = v256_min_s16(v256_max_s16(res, min), max);
        res = v256_pack_s16_u8(res, res);

        p0 = res.val[0];
        *((uint32_t *)(dst + i * dstride)) = _mm_cvtsi128_si32(_mm_srli_si128(p0, 12));
        *((uint32_t *)(dst + (i + (1 * subsampling_factor)) * dstride)) = _mm_cvtsi128_si32(_mm_srli_si128(p0, 8));
        *((uint32_t *)(dst + (i + (2 * subsampling_factor)) * dstride)) = _mm_cvtsi128_si32(_mm_srli_si128(p0, 4));
        *((uint32_t *)(dst + (i + (3 * subsampling_factor)) * dstride)) = _mm_cvtsi128_si32(p0);
    }
}

void svt_av1_cdef_filter_block_8xn_16_sse4_1(uint16_t *dst, int dstride, const uint16_t *in,
                                         int pri_strength, int sec_strength, int dir,
                                             int pri_damping, int sec_damping, int coeff_shift,
                                             uint8_t height, uint8_t subsampling_factor) {
    int  i;
    v256 sum, p0, p1, p2, p3, row, res;
    v256 max, min, large = v256_dup_16(CDEF_VERY_LARGE);
    int  po1  = eb_cdef_directions[dir][0];
    int  po2  = eb_cdef_directions[dir][1];
    int  s1o1 = eb_cdef_directions[(dir + 2) & 7][0];
    int  s1o2 = eb_cdef_directions[(dir + 2) & 7][1];
    int  s2o1 = eb_cdef_directions[(dir + 6) & 7][0];
    int  s2o2 = eb_cdef_directions[(dir + 6) & 7][1];

    const int *pri_taps = eb_cdef_pri_taps[(pri_strength >> coeff_shift) & 1];
    const int *sec_taps = eb_cdef_sec_taps[(pri_strength >> coeff_shift) & 1];

    if (pri_strength)
        pri_damping = AOMMAX(0, pri_damping - get_msb(pri_strength));
    if (sec_strength)
        sec_damping = AOMMAX(0, sec_damping - get_msb(sec_strength));

    for (i = 0; i < height; i += (2 * subsampling_factor)) {
        sum = v256_zero();
        row = v256_from_v128(_mm_loadu_si128((__m128i *)(in + i * CDEF_BSTRIDE)),
                             _mm_loadu_si128((__m128i *)(in + (i + subsampling_factor) * CDEF_BSTRIDE)));

        min = max = row;
        // Primary near taps
        p0  = v256_from_v128(_mm_loadu_si128((__m128i *)(in+i * CDEF_BSTRIDE + po1)),
                             _mm_loadu_si128((__m128i *)(in+(i + subsampling_factor) * CDEF_BSTRIDE + po1)));
        p1  = v256_from_v128(_mm_loadu_si128((__m128i *)(in+i * CDEF_BSTRIDE - po1)),
                             _mm_loadu_si128((__m128i *)(in+(i + subsampling_factor) * CDEF_BSTRIDE - po1)));
        max = v256_max_s16(v256_max_s16(max, v256_andn(p0, v256_cmpeq_16(p0, large))),
                           v256_andn(p1, v256_cmpeq_16(p1, large)));
        min = v256_min_s16(v256_min_s16(min, p0), p1);
        p0  = constrain16(p0, row, pri_strength, pri_damping);
        p1  = constrain16(p1, row, pri_strength, pri_damping);

        // sum += pri_taps[0] * (p0 + p1)
        sum = v256_add_16(sum, v256_mullo_s16(v256_dup_16(pri_taps[0]), v256_add_16(p0, p1)));

        // Primary far taps
        p0  = v256_from_v128(_mm_loadu_si128((__m128i *)(in+i * CDEF_BSTRIDE + po2)),
                             _mm_loadu_si128((__m128i *)(in+(i + subsampling_factor) * CDEF_BSTRIDE + po2)));
        p1  = v256_from_v128(_mm_loadu_si128((__m128i *)(in+i * CDEF_BSTRIDE - po2)),
                             _mm_loadu_si128((__m128i *)(in+(i + subsampling_factor) * CDEF_BSTRIDE - po2)));
        max = v256_max_s16(v256_max_s16(max, v256_andn(p0, v256_cmpeq_16(p0, large))),
                           v256_andn(p1, v256_cmpeq_16(p1, large)));
        min = v256_min_s16(v256_min_s16(min, p0), p1);
        p0  = constrain16(p0, row, pri_strength, pri_damping);
        p1  = constrain16(p1, row, pri_strength, pri_damping);

        // sum += pri_taps[1] * (p0 + p1)
        sum = v256_add_16(sum, v256_mullo_s16(v256_dup_16(pri_taps[1]), v256_add_16(p0, p1)));

        // Secondary near taps
        p0  = v256_from_v128(_mm_loadu_si128((__m128i *)(in+i * CDEF_BSTRIDE + s1o1)),
                             _mm_loadu_si128((__m128i *)(in+(i + subsampling_factor) * CDEF_BSTRIDE + s1o1)));
        p1  = v256_from_v128(_mm_loadu_si128((__m128i *)(in+i * CDEF_BSTRIDE - s1o1)),
                             _mm_loadu_si128((__m128i *)(in+(i + subsampling_factor) * CDEF_BSTRIDE - s1o1)));
        p2  = v256_from_v128(_mm_loadu_si128((__m128i *)(in+i * CDEF_BSTRIDE + s2o1)),
                             _mm_loadu_si128((__m128i *)(in+(i + subsampling_factor) * CDEF_BSTRIDE + s2o1)));
        p3  = v256_from_v128(_mm_loadu_si128((__m128i *)(in+i * CDEF_BSTRIDE - s2o1)),
                             _mm_loadu_si128((__m128i *)(in+(i + subsampling_factor) * CDEF_BSTRIDE - s2o1)));
        max = v256_max_s16(v256_max_s16(max, v256_andn(p0, v256_cmpeq_16(p0, large))),
                           v256_andn(p1, v256_cmpeq_16(p1, large)));
        max = v256_max_s16(v256_max_s16(max, v256_andn(p2, v256_cmpeq_16(p2, large))),
                           v256_andn(p3, v256_cmpeq_16(p3, large)));
        min = v256_min_s16(v256_min_s16(v256_min_s16(v256_min_s16(min, p0), p1), p2), p3);
        p0  = constrain16(p0, row, sec_strength, sec_damping);
        p1  = constrain16(p1, row, sec_strength, sec_damping);
        p2  = constrain16(p2, row, sec_strength, sec_damping);
        p3  = constrain16(p3, row, sec_strength, sec_damping);

        // sum += sec_taps[0] * (p0 + p1 + p2 + p3)
        sum = v256_add_16(sum,
                          v256_mullo_s16(v256_dup_16(sec_taps[0]),
                                         v256_add_16(v256_add_16(p0, p1), v256_add_16(p2, p3))));

        // Secondary far taps
        p0  = v256_from_v128(_mm_loadu_si128((__m128i *)(in+i * CDEF_BSTRIDE + s1o2)),
                             _mm_loadu_si128((__m128i *)(in+(i + subsampling_factor) * CDEF_BSTRIDE + s1o2)));
        p1  = v256_from_v128(_mm_loadu_si128((__m128i *)(in+i * CDEF_BSTRIDE - s1o2)),
                             _mm_loadu_si128((__m128i *)(in+(i + subsampling_factor) * CDEF_BSTRIDE - s1o2)));
        p2  = v256_from_v128(_mm_loadu_si128((__m128i *)(in+i * CDEF_BSTRIDE + s2o2)),
                             _mm_loadu_si128((__m128i *)(in+(i + subsampling_factor) * CDEF_BSTRIDE + s2o2)));
        p3  = v256_from_v128(_mm_loadu_si128((__m128i *)(in+i * CDEF_BSTRIDE - s2o2)),
                             _mm_loadu_si128((__m128i *)(in+(i + subsampling_factor) * CDEF_BSTRIDE - s2o2)));
        max = v256_max_s16(v256_max_s16(max, v256_andn(p0, v256_cmpeq_16(p0, large))),
                           v256_andn(p1, v256_cmpeq_16(p1, large)));
        max = v256_max_s16(v256_max_s16(max, v256_andn(p2, v256_cmpeq_16(p2, large))),
                           v256_andn(p3, v256_cmpeq_16(p3, large)));
        min = v256_min_s16(v256_min_s16(v256_min_s16(v256_min_s16(min, p0), p1), p2), p3);
        p0  = constrain16(p0, row, sec_strength, sec_damping);
        p1  = constrain16(p1, row, sec_strength, sec_damping);
        p2  = constrain16(p2, row, sec_strength, sec_damping);
        p3  = constrain16(p3, row, sec_strength, sec_damping);

        // sum += sec_taps[1] * (p0 + p1 + p2 + p3)
        sum = v256_add_16(sum,
                          v256_mullo_s16(v256_dup_16(sec_taps[1]),
                                         v256_add_16(v256_add_16(p0, p1), v256_add_16(p2, p3))));

        // res = row + ((sum - (sum < 0) + 8) >> 4)
        sum = v256_add_16(sum, v256_cmplt_s16(sum, v256_zero()));
        res = v256_add_16(sum, v256_dup_16(8));
        res = v256_shr_n_s16(res, 4);
        res = v256_add_16(row, res);
        res = v256_min_s16(v256_max_s16(res, min), max);

        _mm_storeu_si128((__m128i *)(dst + i * dstride), res.val[1]);
        _mm_storeu_si128((__m128i *)(dst + (i + subsampling_factor) * dstride), res.val[0]);

    }
}

void svt_av1_cdef_filter_block_4xn_16_sse4_1(uint16_t *dst, int dstride, const uint16_t *in,
                                         int pri_strength, int sec_strength, int dir,
                                         int pri_damping, int sec_damping, int coeff_shift, uint8_t height, uint8_t subsampling_factor) {
    int  i;
    v256 p0, p1, p2, p3, sum, row, res;
    v256 max, min, large = v256_dup_16(CDEF_VERY_LARGE);
    int  po1  = eb_cdef_directions[dir][0];
    int  po2  = eb_cdef_directions[dir][1];
    int  s1o1 = eb_cdef_directions[(dir + 2) & 7][0];
    int  s1o2 = eb_cdef_directions[(dir + 2) & 7][1];
    int  s2o1 = eb_cdef_directions[(dir + 6) & 7][0];
    int  s2o2 = eb_cdef_directions[(dir + 6) & 7][1];

    const int *pri_taps = eb_cdef_pri_taps[(pri_strength >> coeff_shift) & 1];
    const int *sec_taps = eb_cdef_sec_taps[(pri_strength >> coeff_shift) & 1];

    if (pri_strength)
        pri_damping = AOMMAX(0, pri_damping - get_msb(pri_strength));
    if (sec_strength)
        sec_damping = AOMMAX(0, sec_damping - get_msb(sec_strength));
    for (i = 0; i < height; i += (4 * subsampling_factor)) {
        sum = v256_zero();
        row = v256_from_v64(_mm_loadl_epi64((__m128i *)(in + i * CDEF_BSTRIDE)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 1 * subsampling_factor) * CDEF_BSTRIDE)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 2 * subsampling_factor) * CDEF_BSTRIDE)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 3 * subsampling_factor) * CDEF_BSTRIDE)));
        min = max = row;

        // Primary near taps
        p0  = v256_from_v64(_mm_loadl_epi64((__m128i *)(in + i * CDEF_BSTRIDE + po1)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 1 * subsampling_factor) * CDEF_BSTRIDE + po1)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 2 * subsampling_factor) * CDEF_BSTRIDE + po1)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 3 * subsampling_factor) * CDEF_BSTRIDE + po1)));
        p1  = v256_from_v64(_mm_loadl_epi64((__m128i *)(in + i * CDEF_BSTRIDE - po1)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 1 * subsampling_factor) * CDEF_BSTRIDE - po1)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 2 * subsampling_factor) * CDEF_BSTRIDE - po1)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 3 * subsampling_factor) * CDEF_BSTRIDE - po1)));
        max = v256_max_s16(v256_max_s16(max, v256_andn(p0, v256_cmpeq_16(p0, large))),
                           v256_andn(p1, v256_cmpeq_16(p1, large)));
        min = v256_min_s16(v256_min_s16(min, p0), p1);
        p0  = constrain16(p0, row, pri_strength, pri_damping);
        p1  = constrain16(p1, row, pri_strength, pri_damping);

        // sum += pri_taps[0] * (p0 + p1)
        sum = v256_add_16(sum, v256_mullo_s16(v256_dup_16(pri_taps[0]), v256_add_16(p0, p1)));

        // Primary far taps
        p0  = v256_from_v64(_mm_loadl_epi64((__m128i *)(in + i * CDEF_BSTRIDE + po2)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 1 * subsampling_factor) * CDEF_BSTRIDE + po2)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 2 * subsampling_factor) * CDEF_BSTRIDE + po2)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 3 * subsampling_factor) * CDEF_BSTRIDE + po2)));
        p1  = v256_from_v64(_mm_loadl_epi64((__m128i *)(in + i * CDEF_BSTRIDE - po2)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 1 * subsampling_factor) * CDEF_BSTRIDE - po2)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 2 * subsampling_factor) * CDEF_BSTRIDE - po2)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 3 * subsampling_factor) * CDEF_BSTRIDE - po2)));
        max = v256_max_s16(v256_max_s16(max, v256_andn(p0, v256_cmpeq_16(p0, large))),
                           v256_andn(p1, v256_cmpeq_16(p1, large)));
        min = v256_min_s16(v256_min_s16(min, p0), p1);
        p0  = constrain16(p0, row, pri_strength, pri_damping);
        p1  = constrain16(p1, row, pri_strength, pri_damping);

        // sum += pri_taps[1] * (p0 + p1)
        sum = v256_add_16(sum, v256_mullo_s16(v256_dup_16(pri_taps[1]), v256_add_16(p0, p1)));

        // Secondary near taps
        p0  = v256_from_v64(_mm_loadl_epi64((__m128i *)(in + i * CDEF_BSTRIDE + s1o1)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 1 * subsampling_factor) * CDEF_BSTRIDE + s1o1)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 2 * subsampling_factor) * CDEF_BSTRIDE + s1o1)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 3 * subsampling_factor) * CDEF_BSTRIDE + s1o1)));
        p1  = v256_from_v64(_mm_loadl_epi64((__m128i *)(in + i * CDEF_BSTRIDE - s1o1)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 1 * subsampling_factor) * CDEF_BSTRIDE - s1o1)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 2 * subsampling_factor) * CDEF_BSTRIDE - s1o1)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 3 * subsampling_factor) * CDEF_BSTRIDE - s1o1)));
        p2  = v256_from_v64(_mm_loadl_epi64((__m128i *)(in + i * CDEF_BSTRIDE + s2o1)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 1 * subsampling_factor) * CDEF_BSTRIDE + s2o1)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 2 * subsampling_factor) * CDEF_BSTRIDE + s2o1)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 3 * subsampling_factor) * CDEF_BSTRIDE + s2o1)));
        p3  = v256_from_v64(_mm_loadl_epi64((__m128i *)(in + i * CDEF_BSTRIDE - s2o1)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 1 * subsampling_factor) * CDEF_BSTRIDE - s2o1)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 2 * subsampling_factor) * CDEF_BSTRIDE - s2o1)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 3 * subsampling_factor) * CDEF_BSTRIDE - s2o1)));
        max = v256_max_s16(v256_max_s16(max, v256_andn(p0, v256_cmpeq_16(p0, large))),
                           v256_andn(p1, v256_cmpeq_16(p1, large)));
        max = v256_max_s16(v256_max_s16(max, v256_andn(p2, v256_cmpeq_16(p2, large))),
                           v256_andn(p3, v256_cmpeq_16(p3, large)));
        min = v256_min_s16(v256_min_s16(v256_min_s16(v256_min_s16(min, p0), p1), p2), p3);
        p0  = constrain16(p0, row, sec_strength, sec_damping);
        p1  = constrain16(p1, row, sec_strength, sec_damping);
        p2  = constrain16(p2, row, sec_strength, sec_damping);
        p3  = constrain16(p3, row, sec_strength, sec_damping);

        // sum += sec_taps[0] * (p0 + p1 + p2 + p3)
        sum = v256_add_16(sum,
                          v256_mullo_s16(v256_dup_16(sec_taps[0]),
                                         v256_add_16(v256_add_16(p0, p1), v256_add_16(p2, p3))));

        // Secondary far taps
        p0  = v256_from_v64(_mm_loadl_epi64((__m128i *)(in + i * CDEF_BSTRIDE + s1o2)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 1 * subsampling_factor) * CDEF_BSTRIDE + s1o2)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 2 * subsampling_factor) * CDEF_BSTRIDE + s1o2)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 3 * subsampling_factor) * CDEF_BSTRIDE + s1o2)));
        p1  = v256_from_v64(_mm_loadl_epi64((__m128i *)(in + i * CDEF_BSTRIDE - s1o2)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 1 * subsampling_factor) * CDEF_BSTRIDE - s1o2)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 2 * subsampling_factor) * CDEF_BSTRIDE - s1o2)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 3 * subsampling_factor) * CDEF_BSTRIDE - s1o2)));
        p2  = v256_from_v64(_mm_loadl_epi64((__m128i *)(in + i * CDEF_BSTRIDE + s2o2)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 1 * subsampling_factor) * CDEF_BSTRIDE + s2o2)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 2 * subsampling_factor) * CDEF_BSTRIDE + s2o2)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 3 * subsampling_factor) * CDEF_BSTRIDE + s2o2)));
        p3  = v256_from_v64(_mm_loadl_epi64((__m128i *)(in + i * CDEF_BSTRIDE - s2o2)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 1 * subsampling_factor) * CDEF_BSTRIDE - s2o2)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 2 * subsampling_factor) * CDEF_BSTRIDE - s2o2)),
                            _mm_loadl_epi64((__m128i *)(in + (i + 3 * subsampling_factor) * CDEF_BSTRIDE - s2o2)));
        max = v256_max_s16(v256_max_s16(max, v256_andn(p0, v256_cmpeq_16(p0, large))),
                           v256_andn(p1, v256_cmpeq_16(p1, large)));
        max = v256_max_s16(v256_max_s16(max, v256_andn(p2, v256_cmpeq_16(p2, large))),
                           v256_andn(p3, v256_cmpeq_16(p3, large)));
        min = v256_min_s16(v256_min_s16(v256_min_s16(v256_min_s16(min, p0), p1), p2), p3);
        p0  = constrain16(p0, row, sec_strength, sec_damping);
        p1  = constrain16(p1, row, sec_strength, sec_damping);
        p2  = constrain16(p2, row, sec_strength, sec_damping);
        p3  = constrain16(p3, row, sec_strength, sec_damping);

        // sum += sec_taps[1] * (p0 + p1 + p2 + p3)
        sum = v256_add_16(sum,
                          v256_mullo_s16(v256_dup_16(sec_taps[1]),
                                         v256_add_16(v256_add_16(p0, p1), v256_add_16(p2, p3))));

        // res = row + ((sum - (sum < 0) + 8) >> 4)
        sum = v256_add_16(sum, v256_cmplt_s16(sum, v256_zero()));
        res = v256_add_16(sum, v256_dup_16(8));
        res = v256_shr_n_s16(res, 4);
        res = v256_add_16(row, res);
        res = v256_min_s16(v256_max_s16(res, min), max);

        _mm_storel_epi64((__m128i *)(dst + i * dstride), _mm_srli_si128(res.val[1], 8));
        _mm_storel_epi64((__m128i *)(dst + (i + 1 * subsampling_factor) * dstride), res.val[1]);
        _mm_storel_epi64((__m128i *)(dst + (i + 2 * subsampling_factor) * dstride), _mm_srli_si128(res.val[0], 8));
        _mm_storel_epi64((__m128i *)(dst + (i + 3 * subsampling_factor) * dstride), res.val[0]);
    }
}

void svt_av1_cdef_filter_block_sse4_1(uint8_t *dst8, uint16_t *dst16, int dstride,
                                  const uint16_t *in, int pri_strength,
                                  int sec_strength, int dir, int pri_damping, int sec_damping, int bsize,
                                      int coeff_shift, uint8_t subsampling_factor) {
  if (dst8) {
    if (bsize == BLOCK_8X8) {
      svt_av1_cdef_filter_block_8xn_8_sse4_1
      (dst8, dstride, in, pri_strength, sec_strength, dir, pri_damping,
       sec_damping, coeff_shift, 8, subsampling_factor);
    } else if (bsize == BLOCK_4X8) {
      svt_av1_cdef_filter_block_4xn_8_sse4_1
      (dst8, dstride, in, pri_strength, sec_strength, dir, pri_damping,
       sec_damping, coeff_shift, 8, subsampling_factor);
    } else if (bsize == BLOCK_8X4) {
      svt_av1_cdef_filter_block_8xn_8_sse4_1
      (dst8, dstride, in, pri_strength, sec_strength, dir, pri_damping,
       sec_damping, coeff_shift, 4, subsampling_factor);
    } else {
      svt_av1_cdef_filter_block_4xn_8_sse4_1
      (dst8, dstride, in, pri_strength, sec_strength, dir, pri_damping,
       sec_damping, coeff_shift, 4, 1);
    }
  } else {
    if (bsize == BLOCK_8X8) {
      svt_av1_cdef_filter_block_8xn_16_sse4_1
      (dst16, dstride, in, pri_strength, sec_strength, dir, pri_damping,
       sec_damping, coeff_shift, 8, subsampling_factor);
    } else if (bsize == BLOCK_4X8) {
      svt_av1_cdef_filter_block_4xn_16_sse4_1
      (dst16, dstride, in, pri_strength, sec_strength, dir, pri_damping,
       sec_damping, coeff_shift, 8, subsampling_factor);
    } else if (bsize == BLOCK_8X4) {
      svt_av1_cdef_filter_block_8xn_16_sse4_1
      (dst16, dstride, in, pri_strength, sec_strength, dir, pri_damping,
       sec_damping, coeff_shift, 4, subsampling_factor);
    } else {
      assert(bsize == BLOCK_4X4);
      svt_av1_cdef_filter_block_4xn_16_sse4_1
      (dst16, dstride, in, pri_strength, sec_strength, dir, pri_damping,
       sec_damping, coeff_shift, 4, 1);
    }
  }
}
