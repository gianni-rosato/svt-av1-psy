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
#include <immintrin.h>

#include "aom_dsp_rtcd.h"
#include "EbCdef.h"
#include "EbBitstreamUnit.h"

#include "v128_intrinsics_x86.h"
#include "v256_intrinsics_x86.h"

#include "aom_dsp_rtcd.h"

#define SIMD_FUNC(name) name##_avx2

#if defined(__SSE4_1__)
#undef CDEF_AVX_OPT
#define CDEF_AVX_OPT 0
#endif

/* partial A is a 16-bit vector of the form:
[x8 x7 x6 x5 x4 x3 x2 x1] and partial B has the form:
[0  y1 y2 y3 y4 y5 y6 y7].
This function computes (x1^2+y1^2)*C1 + (x2^2+y2^2)*C2 + ...
(x7^2+y2^7)*C7 + (x8^2+0^2)*C8 where the C1..C8 constants are in const1
and const2. */
#if CDEF_AVX_OPT
static INLINE v256 fold_mul_and_sum(v256 partial, v256 const_var) {
    partial = _mm256_shuffle_epi8(partial, _mm256_set_epi32(
        0x0f0e0100, 0x03020504, 0x07060908, 0x0b0a0d0c,
        0x0f0e0d0c, 0x0b0a0908, 0x07060504, 0x03020100));
    partial = _mm256_permutevar8x32_epi32(partial,
        _mm256_set_epi32(7, 3, 6, 2, 5, 1, 4, 0));
    partial = _mm256_shuffle_epi8(partial, _mm256_set_epi32(
        0x0f0e0b0a, 0x0d0c0908, 0x07060302, 0x05040100,
        0x0f0e0b0a, 0x0d0c0908, 0x07060302, 0x05040100));
    partial = _mm256_madd_epi16(partial, partial);
    partial = _mm256_mullo_epi32(partial, const_var);
    return partial;
}

#else
static INLINE v128 fold_mul_and_sum(v128 partiala, v128 partialb, v128 const1,
    v128 const2) {
    v128 tmp;
    /* Reverse partial B. */
    partialb = v128_shuffle_8(
        partialb, v128_from_32(0x0f0e0100, 0x03020504, 0x07060908, 0x0b0a0d0c));
    /* Interleave the x and y values of identical indices and pair x8 with 0. */
    tmp = partiala;
    partiala = v128_ziplo_16(partialb, partiala);
    partialb = v128_ziphi_16(partialb, tmp);
    /* Square and add the corresponding x and y values. */
    partiala = v128_madd_s16(partiala, partiala);
    partialb = v128_madd_s16(partialb, partialb);
    /* Multiply by constant. */
    partiala = v128_mullo_s32(partiala, const1);
    partialb = v128_mullo_s32(partialb, const2);
    /* Sum all results. */
    partiala = v128_add_32(partiala, partialb);
    return partiala;
}
#endif
static INLINE v128 hsum4(v128 x0, v128 x1, v128 x2, v128 x3) {
    v128 t0, t1, t2, t3;
    t0 = v128_ziplo_32(x1, x0);
    t1 = v128_ziplo_32(x3, x2);
    t2 = v128_ziphi_32(x1, x0);
    t3 = v128_ziphi_32(x3, x2);
    x0 = v128_ziplo_64(t1, t0);
    x1 = v128_ziphi_64(t1, t0);
    x2 = v128_ziplo_64(t3, t2);
    x3 = v128_ziphi_64(t3, t2);
    return v128_add_32(v128_add_32(x0, x1), v128_add_32(x2, x3));
}

/* Computes cost for directions 0, 5, 6 and 7. We can call this function again
to compute the remaining directions. */
#if CDEF_AVX_OPT
static INLINE void compute_directions(v128 lines[8], int32_t tmp_cost1[4]) {
    v128 partial6;
    v128 tmp;

    v256 partial4;
    v256 partial5;
    v256 partial7;
    v256 tmp_avx2;
#else
static INLINE v128 compute_directions(v128 lines[8], int32_t tmp_cost1[4]) {
    v128 partial4a, partial4b, partial5a, partial5b, partial7a, partial7b;
    v128 partial6;
    v128 tmp;
#endif
    /* Partial sums for lines 0 and 1. */
#if CDEF_AVX_OPT
    partial4 = _mm256_insertf128_si256(_mm256_castsi128_si256(
        v128_shl_n_byte(lines[0], 14)), v128_shr_n_byte(lines[0], 2), 0x1);
    tmp_avx2 = _mm256_insertf128_si256(_mm256_castsi128_si256(
        v128_shl_n_byte(lines[1], 12)), v128_shr_n_byte(lines[1], 4), 0x1);
    partial4 = _mm256_add_epi16(partial4, tmp_avx2);
#else
    partial4a = v128_shl_n_byte(lines[0], 14);
    partial4b = v128_shr_n_byte(lines[0], 2);
    partial4a = v128_add_16(partial4a, v128_shl_n_byte(lines[1], 12));
    partial4b = v128_add_16(partial4b, v128_shr_n_byte(lines[1], 4));
#endif
    tmp = v128_add_16(lines[0], lines[1]);
#if CDEF_AVX_OPT
    partial5 = _mm256_insertf128_si256(_mm256_castsi128_si256(
        v128_shl_n_byte(tmp, 10)), v128_shr_n_byte(tmp, 6), 0x1);
    partial7 = _mm256_insertf128_si256(_mm256_castsi128_si256(
        v128_shl_n_byte(tmp, 4)), v128_shr_n_byte(tmp, 12), 0x1);
#else
    partial5a = v128_shl_n_byte(tmp, 10);
    partial5b = v128_shr_n_byte(tmp, 6);
    partial7a = v128_shl_n_byte(tmp, 4);
    partial7b = v128_shr_n_byte(tmp, 12);
#endif
    partial6 = tmp;

    /* Partial sums for lines 2 and 3. */
#if CDEF_AVX_OPT
    tmp_avx2 = _mm256_insertf128_si256(_mm256_castsi128_si256(
        v128_shl_n_byte(lines[2], 10)), v128_shr_n_byte(lines[2], 6), 0x1);
    partial4 = _mm256_add_epi16(partial4, tmp_avx2);
    tmp_avx2 = _mm256_insertf128_si256(_mm256_castsi128_si256(
        v128_shl_n_byte(lines[3], 8)), v128_shr_n_byte(lines[3], 8), 0x1);
    partial4 = _mm256_add_epi16(partial4, tmp_avx2);
#else
    partial4a = v128_add_16(partial4a, v128_shl_n_byte(lines[2], 10));
    partial4b = v128_add_16(partial4b, v128_shr_n_byte(lines[2], 6));
    partial4a = v128_add_16(partial4a, v128_shl_n_byte(lines[3], 8));
    partial4b = v128_add_16(partial4b, v128_shr_n_byte(lines[3], 8));
#endif
    tmp = v128_add_16(lines[2], lines[3]);
#if CDEF_AVX_OPT
    tmp_avx2 = _mm256_insertf128_si256(_mm256_castsi128_si256(
        v128_shl_n_byte(tmp, 8)), v128_shr_n_byte(tmp, 8), 0x1);
    partial5 = _mm256_add_epi16(partial5, tmp_avx2);
    tmp_avx2 = _mm256_insertf128_si256(_mm256_castsi128_si256(
        v128_shl_n_byte(tmp, 6)), v128_shr_n_byte(tmp, 10), 0x1);
    partial7 = _mm256_add_epi16(partial7, tmp_avx2);
#else
    partial5a = v128_add_16(partial5a, v128_shl_n_byte(tmp, 8));
    partial5b = v128_add_16(partial5b, v128_shr_n_byte(tmp, 8));
    partial7a = v128_add_16(partial7a, v128_shl_n_byte(tmp, 6));
    partial7b = v128_add_16(partial7b, v128_shr_n_byte(tmp, 10));
#endif
    partial6 = v128_add_16(partial6, tmp);

    /* Partial sums for lines 4 and 5. */
#if CDEF_AVX_OPT
    tmp_avx2 = _mm256_insertf128_si256(_mm256_castsi128_si256(
        v128_shl_n_byte(lines[4], 6)), v128_shr_n_byte(lines[4], 10), 0x1);
    partial4 = _mm256_add_epi16(partial4, tmp_avx2);
    tmp_avx2 = _mm256_insertf128_si256(_mm256_castsi128_si256(
        v128_shl_n_byte(lines[5], 4)), v128_shr_n_byte(lines[5], 12), 0x1);
    partial4 = _mm256_add_epi16(partial4, tmp_avx2);
#else
    partial4a = v128_add_16(partial4a, v128_shl_n_byte(lines[4], 6));
    partial4b = v128_add_16(partial4b, v128_shr_n_byte(lines[4], 10));
    partial4a = v128_add_16(partial4a, v128_shl_n_byte(lines[5], 4));
    partial4b = v128_add_16(partial4b, v128_shr_n_byte(lines[5], 12));
#endif
    tmp = v128_add_16(lines[4], lines[5]);
#if CDEF_AVX_OPT
    tmp_avx2 = _mm256_insertf128_si256(_mm256_castsi128_si256(
        v128_shl_n_byte(tmp, 6)), v128_shr_n_byte(tmp, 10), 0x1);
    partial5 = _mm256_add_epi16(partial5, tmp_avx2);
    tmp_avx2 = _mm256_insertf128_si256(_mm256_castsi128_si256(
        v128_shl_n_byte(tmp, 8)), v128_shr_n_byte(tmp, 8), 0x1);
    partial7 = _mm256_add_epi16(partial7, tmp_avx2);

#else
    partial5a = v128_add_16(partial5a, v128_shl_n_byte(tmp, 6));
    partial5b = v128_add_16(partial5b, v128_shr_n_byte(tmp, 10));
    partial7a = v128_add_16(partial7a, v128_shl_n_byte(tmp, 8));
    partial7b = v128_add_16(partial7b, v128_shr_n_byte(tmp, 8));
#endif
    partial6 = v128_add_16(partial6, tmp);

    /* Partial sums for lines 6 and 7. */
#if CDEF_AVX_OPT
    tmp_avx2 = _mm256_insertf128_si256(_mm256_castsi128_si256(
        v128_shl_n_byte(lines[6], 2)), v128_shr_n_byte(lines[6], 14), 0x1);
    partial4 = _mm256_add_epi16(partial4, tmp_avx2);
    tmp_avx2 = _mm256_insertf128_si256(_mm256_setzero_si256(), lines[7], 0x0);
    partial4 = _mm256_add_epi16(partial4, tmp_avx2);
#else
    partial4a = v128_add_16(partial4a, v128_shl_n_byte(lines[6], 2));
    partial4b = v128_add_16(partial4b, v128_shr_n_byte(lines[6], 14));
    partial4a = v128_add_16(partial4a, lines[7]);
#endif
    tmp = v128_add_16(lines[6], lines[7]);
#if  CDEF_AVX_OPT
    tmp_avx2 = _mm256_insertf128_si256(_mm256_castsi128_si256(
        v128_shl_n_byte(tmp, 4)), v128_shr_n_byte(tmp, 12), 0x1);
    partial5 = _mm256_add_epi16(partial5, tmp_avx2);
    tmp_avx2 = _mm256_insertf128_si256(_mm256_castsi128_si256(
        v128_shl_n_byte(tmp, 10)), v128_shr_n_byte(tmp, 6), 0x1);
    partial7 = _mm256_add_epi16(partial7, tmp_avx2);
#else
    partial5a = v128_add_16(partial5a, v128_shl_n_byte(tmp, 4));
    partial5b = v128_add_16(partial5b, v128_shr_n_byte(tmp, 12));
    partial7a = v128_add_16(partial7a, v128_shl_n_byte(tmp, 10));
    partial7b = v128_add_16(partial7b, v128_shr_n_byte(tmp, 6));
#endif
    partial6 = v128_add_16(partial6, tmp);

    /* Compute costs in terms of partial sums. */
#if CDEF_AVX_OPT
    partial4 = fold_mul_and_sum(partial4, _mm256_set_epi32(
        105, 120, 140, 168, 210, 280, 420, 840));
    partial7 = fold_mul_and_sum(partial7, _mm256_set_epi32(
        105, 105, 105, 140, 210, 420, 0, 0));
    partial5 = fold_mul_and_sum(partial5, _mm256_set_epi32(
        105, 105, 105, 140, 210, 420, 0, 0));
#else
    partial4a =
        fold_mul_and_sum(partial4a, partial4b, v128_from_32(210, 280, 420, 840),
            v128_from_32(105, 120, 140, 168));
    partial7a =
        fold_mul_and_sum(partial7a, partial7b, v128_from_32(210, 420, 0, 0),
            v128_from_32(105, 105, 105, 140));
    partial5a =
        fold_mul_and_sum(partial5a, partial5b, v128_from_32(210, 420, 0, 0),
            v128_from_32(105, 105, 105, 140));
#endif
    partial6 = v128_madd_s16(partial6, partial6);
    partial6 = v128_mullo_s32(partial6, v128_dup_32(105));
#if CDEF_AVX_OPT
    v128 a, b, c;
    a = _mm_add_epi32(_mm256_castsi256_si128(partial4),
        _mm256_extracti128_si256(partial4, 1));
    b = _mm_add_epi32(_mm256_castsi256_si128(partial5),
        _mm256_extracti128_si256(partial5, 1));
    c = _mm_add_epi32(_mm256_castsi256_si128(partial7),
        _mm256_extracti128_si256(partial7, 1));

    v128_store_unaligned(tmp_cost1, hsum4(a, b, partial6, c));
#else
    partial4a = hsum4(partial4a, partial5a, partial6, partial7a);
    v128_store_unaligned(tmp_cost1, partial4a);
    return partial4a;
#endif
}

/* transpose and reverse the order of the lines -- equivalent to a 90-degree
counter-clockwise rotation of the pixels. */
static INLINE void array_reverse_transpose_8x8(v128 *in, v128 *res) {
    const v128 tr0_0 = v128_ziplo_16(in[1], in[0]);
    const v128 tr0_1 = v128_ziplo_16(in[3], in[2]);
    const v128 tr0_2 = v128_ziphi_16(in[1], in[0]);
    const v128 tr0_3 = v128_ziphi_16(in[3], in[2]);
    const v128 tr0_4 = v128_ziplo_16(in[5], in[4]);
    const v128 tr0_5 = v128_ziplo_16(in[7], in[6]);
    const v128 tr0_6 = v128_ziphi_16(in[5], in[4]);
    const v128 tr0_7 = v128_ziphi_16(in[7], in[6]);

    const v128 tr1_0 = v128_ziplo_32(tr0_1, tr0_0);
    const v128 tr1_1 = v128_ziplo_32(tr0_5, tr0_4);
    const v128 tr1_2 = v128_ziphi_32(tr0_1, tr0_0);
    const v128 tr1_3 = v128_ziphi_32(tr0_5, tr0_4);
    const v128 tr1_4 = v128_ziplo_32(tr0_3, tr0_2);
    const v128 tr1_5 = v128_ziplo_32(tr0_7, tr0_6);
    const v128 tr1_6 = v128_ziphi_32(tr0_3, tr0_2);
    const v128 tr1_7 = v128_ziphi_32(tr0_7, tr0_6);

    res[7] = v128_ziplo_64(tr1_1, tr1_0);
    res[6] = v128_ziphi_64(tr1_1, tr1_0);
    res[5] = v128_ziplo_64(tr1_3, tr1_2);
    res[4] = v128_ziphi_64(tr1_3, tr1_2);
    res[3] = v128_ziplo_64(tr1_5, tr1_4);
    res[2] = v128_ziphi_64(tr1_5, tr1_4);
    res[1] = v128_ziplo_64(tr1_7, tr1_6);
    res[0] = v128_ziphi_64(tr1_7, tr1_6);
}

int32_t SIMD_FUNC(cdef_find_dir)(const uint16_t *img, int32_t stride, int32_t *var,
    int32_t coeff_shift) {
    int32_t i;
    int32_t cost[8];
    int32_t best_cost = 0;
    int32_t best_dir = 0;
    v128 lines[8];
#if CDEF_AVX_OPT
    v128 const_128 = v128_dup_16(128);
#endif
    for (i = 0; i < 8; i++) {
        lines[i] = v128_load_unaligned(&img[i * stride]);
        lines[i] =
#if CDEF_AVX_OPT
            v128_sub_16(v128_shr_s16(lines[i], coeff_shift), const_128);
#else
            v128_sub_16(v128_shr_s16(lines[i], coeff_shift), v128_dup_16(128));
#endif
    }

#if defined(__SSE4_1__)
    /* Compute "mostly vertical" directions. */
    __m128i dir47 = compute_directions(lines, cost + 4);

    array_reverse_transpose_8x8(lines, lines);

    /* Compute "mostly horizontal" directions. */
    __m128i dir03 = compute_directions(lines, cost);

    __m128i max = _mm_max_epi32(dir03, dir47);
    max = _mm_max_epi32(max, _mm_shuffle_epi32(max, _MM_SHUFFLE(1, 0, 3, 2)));
    max = _mm_max_epi32(max, _mm_shuffle_epi32(max, _MM_SHUFFLE(2, 3, 0, 1)));
    best_cost = _mm_cvtsi128_si32(max);
    __m128i t =
        _mm_packs_epi32(_mm_cmpeq_epi32(max, dir03), _mm_cmpeq_epi32(max, dir47));
    best_dir = _mm_movemask_epi8(_mm_packs_epi16(t, t));
    best_dir = get_msb(best_dir ^ (best_dir - 1));  // Count trailing zeros
#else

    /* Compute "mostly vertical" directions. */
    compute_directions(lines, cost + 4);

    array_reverse_transpose_8x8(lines, lines);

    /* Compute "mostly horizontal" directions. */
    compute_directions(lines, cost);

    for (i = 0; i < 8; i++) {
        if (cost[i] > best_cost) {
            best_cost = cost[i];
            best_dir = i;
        }
    }
#endif

    /* Difference between the optimal variance and the variance along the
    orthogonal direction. Again, the sum(x^2) terms cancel out. */
    *var = best_cost - cost[(best_dir + 4) & 7];
    /* We'd normally divide by 840, but dividing by 1024 is close enough
    for what we're going to do with this. */
    *var >>= 10;
    return best_dir;
}

// sign(a-b) * min(abs(a-b), max(0, threshold - (abs(a-b) >> adjdamp)))
#if CDEF_AVX_OPT
SIMD_INLINE v256 constrain16(v256 a, v256 b, v256 threshold,
#else
SIMD_INLINE v256 constrain16(v256 a, v256 b, uint32_t threshold,
#endif
    uint32_t adjdamp) {
    v256 diff = v256_sub_16(a, b);
    const v256 sign = v256_shr_n_s16(diff, 15);
    diff = v256_abs_s16(diff);
    const v256 s =
#if CDEF_AVX_OPT
        v256_ssub_u16(threshold, v256_shr_u16(diff, adjdamp));
#else
        v256_ssub_u16(v256_dup_16(threshold), v256_shr_u16(diff, adjdamp));
#endif
    return v256_xor(v256_add_16(sign, v256_min_s16(diff, s)), sign);
}

// sign(a - b) * min(abs(a - b), max(0, strength - (abs(a - b) >> adjdamp)))
SIMD_INLINE v128 constrain(v256 a, v256 b, uint32_t strength,
    uint32_t adjdamp) {
    const v256 diff16 = v256_sub_16(a, b);
    v128 diff = v128_pack_s16_s8(v256_high_v128(diff16), v256_low_v128(diff16));
    const v128 sign = v128_cmplt_s8(diff, v128_zero());
    diff = v128_abs_s8(diff);
    return v128_xor(
        v128_add_8(sign,
            v128_min_u8(diff, v128_ssub_u8(v128_dup_8(strength),
                v128_shr_u8(diff, adjdamp)))),
        sign);
}

void SIMD_FUNC(cdef_filter_block_4x4_8)(uint8_t *dst, int32_t dstride,
    const uint16_t *in, int32_t pri_strength,
    int32_t sec_strength, int32_t dir,
    int32_t pri_damping, int32_t sec_damping,
    /* AOM_UNUSED*/ int32_t max_unused,
    int32_t coeff_shift) {
    (void)max_unused;
    v128 p0, p1, p2, p3;
    v256 sum, row, tap, res;
    v256 max, min, large = v256_dup_16(CDEF_VERY_LARGE);
    int32_t po1 = cdef_directions[dir][0];
    int32_t po2 = cdef_directions[dir][1];
    int32_t s1o1 = cdef_directions[(dir + 2) & 7][0];
    int32_t s1o2 = cdef_directions[(dir + 2) & 7][1];
    int32_t s2o1 = cdef_directions[(dir + 6) & 7][0];
    int32_t s2o2 = cdef_directions[(dir + 6) & 7][1];

    const int32_t *pri_taps = cdef_pri_taps[(pri_strength >> coeff_shift) & 1];
    const int32_t *sec_taps = cdef_sec_taps[(pri_strength >> coeff_shift) & 1];

    if (pri_strength)
        pri_damping = AOMMAX(0, pri_damping - get_msb(pri_strength));
    if (sec_strength)
        sec_damping = AOMMAX(0, sec_damping - get_msb(sec_strength));

    sum = v256_zero();
#if CDEF_AVX_OPT
    row = _mm256_set_epi64x(*(uint64_t*)(in),
        *(uint64_t*)(in + CDEF_BSTRIDE),
        *(uint64_t*)(in + 2 * CDEF_BSTRIDE),
        *(uint64_t*)(in + 3 * CDEF_BSTRIDE));
#else
    row = v256_from_v64(v64_load_aligned(&in[0 * CDEF_BSTRIDE]),
        v64_load_aligned(&in[1 * CDEF_BSTRIDE]),
        v64_load_aligned(&in[2 * CDEF_BSTRIDE]),
        v64_load_aligned(&in[3 * CDEF_BSTRIDE]));
#endif
    max = min = row;

    if (pri_strength) {
#if CDEF_AVX_OPT
        tap = _mm256_set_epi64x(*(uint64_t*)(in + po1),
            *(uint64_t*)(in + CDEF_BSTRIDE + po1),
            *(uint64_t*)(in + 2 * CDEF_BSTRIDE + po1),
            *(uint64_t*)(in + 3 * CDEF_BSTRIDE + po1));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p0 = constrain(tap, row, pri_strength, pri_damping);
        tap = _mm256_set_epi64x(*(uint64_t*)(in - po1),
            *(uint64_t*)(in + CDEF_BSTRIDE - po1),
            *(uint64_t*)(in + 2 * CDEF_BSTRIDE - po1),
            *(uint64_t*)(in + 3 * CDEF_BSTRIDE - po1));
#else
        // Primary near taps
        tap = v256_from_v64(v64_load_unaligned(&in[0 * CDEF_BSTRIDE + po1]),
            v64_load_unaligned(&in[1 * CDEF_BSTRIDE + po1]),
            v64_load_unaligned(&in[2 * CDEF_BSTRIDE + po1]),
            v64_load_unaligned(&in[3 * CDEF_BSTRIDE + po1]));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p0 = constrain(tap, row, pri_strength, pri_damping);
        tap = v256_from_v64(v64_load_unaligned(&in[0 * CDEF_BSTRIDE - po1]),
            v64_load_unaligned(&in[1 * CDEF_BSTRIDE - po1]),
            v64_load_unaligned(&in[2 * CDEF_BSTRIDE - po1]),
            v64_load_unaligned(&in[3 * CDEF_BSTRIDE - po1]));

#endif
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p1 = constrain(tap, row, pri_strength, pri_damping);

        // sum += pri_taps[0] * (p0 + p1)
        sum = v256_add_16(sum, v256_madd_us8(v256_dup_8(pri_taps[0]),
            v256_from_v128(v128_ziphi_8(p0, p1),
                v128_ziplo_8(p0, p1))));
#if CDEF_AVX_OPT
        // Primary far taps
        tap = _mm256_set_epi64x(*(uint64_t*)(in + po2),
            *(uint64_t*)(in + CDEF_BSTRIDE + po2),
            *(uint64_t*)(in + 2 * CDEF_BSTRIDE + po2),
            *(uint64_t*)(in + 3 * CDEF_BSTRIDE + po2));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p0 = constrain(tap, row, pri_strength, pri_damping);
            tap = _mm256_set_epi64x(*(uint64_t*)(in - po2),
                *(uint64_t*)(in + CDEF_BSTRIDE - po2),
                *(uint64_t*)(in + 2 * CDEF_BSTRIDE - po2),
                *(uint64_t*)(in + 3 * CDEF_BSTRIDE - po2));
#else
        // Primary far taps
        tap = v256_from_v64(v64_load_unaligned(&in[0 * CDEF_BSTRIDE + po2]),
            v64_load_unaligned(&in[1 * CDEF_BSTRIDE + po2]),
            v64_load_unaligned(&in[2 * CDEF_BSTRIDE + po2]),
            v64_load_unaligned(&in[3 * CDEF_BSTRIDE + po2]));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p0 = constrain(tap, row, pri_strength, pri_damping);
        tap = v256_from_v64(v64_load_unaligned(&in[0 * CDEF_BSTRIDE - po2]),
            v64_load_unaligned(&in[1 * CDEF_BSTRIDE - po2]),
            v64_load_unaligned(&in[2 * CDEF_BSTRIDE - po2]),
            v64_load_unaligned(&in[3 * CDEF_BSTRIDE - po2]));
#endif
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p1 = constrain(tap, row, pri_strength, pri_damping);

        // sum += pri_taps[1] * (p0 + p1)
        sum = v256_add_16(sum, v256_madd_us8(v256_dup_8(pri_taps[1]),
            v256_from_v128(v128_ziphi_8(p0, p1),
                v128_ziplo_8(p0, p1))));
    }

    if (sec_strength) {
        // Secondary near taps

#if CDEF_AVX_OPT
        tap = _mm256_set_epi64x(*(uint64_t*)(in + s1o1),
            *(uint64_t*)(in + CDEF_BSTRIDE + s1o1),
            *(uint64_t*)(in + 2 * CDEF_BSTRIDE + s1o1),
            *(uint64_t*)(in + 3 * CDEF_BSTRIDE + s1o1));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p0 = constrain(tap, row, sec_strength, sec_damping);
        tap = _mm256_set_epi64x(*(uint64_t*)(in - s1o1),
            *(uint64_t*)(in + CDEF_BSTRIDE - s1o1),
            *(uint64_t*)(in + 2 * CDEF_BSTRIDE - s1o1),
            *(uint64_t*)(in + 3 * CDEF_BSTRIDE - s1o1));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p1 = constrain(tap, row, sec_strength, sec_damping);
        tap = _mm256_set_epi64x(*(uint64_t*)(in + s2o1),
            *(uint64_t*)(in + CDEF_BSTRIDE + s2o1),
            *(uint64_t*)(in + 2 * CDEF_BSTRIDE + s2o1),
            *(uint64_t*)(in + 3 * CDEF_BSTRIDE + s2o1));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p2 = constrain(tap, row, sec_strength, sec_damping);
        tap = _mm256_set_epi64x(*(uint64_t*)(in - s2o1),
            *(uint64_t*)(in + CDEF_BSTRIDE - s2o1),
            *(uint64_t*)(in + 2 * CDEF_BSTRIDE - s2o1),
            *(uint64_t*)(in + 3 * CDEF_BSTRIDE - s2o1));
#else
        tap = v256_from_v64(v64_load_unaligned(&in[0 * CDEF_BSTRIDE + s1o1]),
            v64_load_unaligned(&in[1 * CDEF_BSTRIDE + s1o1]),
            v64_load_unaligned(&in[2 * CDEF_BSTRIDE + s1o1]),
            v64_load_unaligned(&in[3 * CDEF_BSTRIDE + s1o1]));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p0 = constrain(tap, row, sec_strength, sec_damping);
        tap = v256_from_v64(v64_load_unaligned(&in[0 * CDEF_BSTRIDE - s1o1]),
            v64_load_unaligned(&in[1 * CDEF_BSTRIDE - s1o1]),
            v64_load_unaligned(&in[2 * CDEF_BSTRIDE - s1o1]),
            v64_load_unaligned(&in[3 * CDEF_BSTRIDE - s1o1]));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p1 = constrain(tap, row, sec_strength, sec_damping);
        tap = v256_from_v64(v64_load_unaligned(&in[0 * CDEF_BSTRIDE + s2o1]),
            v64_load_unaligned(&in[1 * CDEF_BSTRIDE + s2o1]),
            v64_load_unaligned(&in[2 * CDEF_BSTRIDE + s2o1]),
            v64_load_unaligned(&in[3 * CDEF_BSTRIDE + s2o1]));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p2 = constrain(tap, row, sec_strength, sec_damping);
        tap = v256_from_v64(v64_load_unaligned(&in[0 * CDEF_BSTRIDE - s2o1]),
            v64_load_unaligned(&in[1 * CDEF_BSTRIDE - s2o1]),
            v64_load_unaligned(&in[2 * CDEF_BSTRIDE - s2o1]),
            v64_load_unaligned(&in[3 * CDEF_BSTRIDE - s2o1]));
#endif
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p3 = constrain(tap, row, sec_strength, sec_damping);

        // sum += sec_taps[0] * (p0 + p1 + p2 + p3)
        p0 = v128_add_8(p0, p1);
        p2 = v128_add_8(p2, p3);
        sum = v256_add_16(sum, v256_madd_us8(v256_dup_8(sec_taps[0]),
            v256_from_v128(v128_ziphi_8(p0, p2),
                v128_ziplo_8(p0, p2))));

        // Secondary far taps
#if CDEF_AVX_OPT
        tap = _mm256_set_epi64x(*(uint64_t*)(in + s1o2),
            *(uint64_t*)(in + CDEF_BSTRIDE + s1o2),
            *(uint64_t*)(in + 2 * CDEF_BSTRIDE + s1o2),
            *(uint64_t*)(in + 3 * CDEF_BSTRIDE + s1o2));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p0 = constrain(tap, row, sec_strength, sec_damping);
        tap = _mm256_set_epi64x(*(uint64_t*)(in - s1o2),
            *(uint64_t*)(in + CDEF_BSTRIDE - s1o2),
            *(uint64_t*)(in + 2 * CDEF_BSTRIDE - s1o2),
            *(uint64_t*)(in + 3 * CDEF_BSTRIDE - s1o2));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p1 = constrain(tap, row, sec_strength, sec_damping);
        tap = _mm256_set_epi64x(*(uint64_t*)(in + s2o2),
            *(uint64_t*)(in + CDEF_BSTRIDE + s2o2),
            *(uint64_t*)(in + 2 * CDEF_BSTRIDE + s2o2),
            *(uint64_t*)(in + 3 * CDEF_BSTRIDE + s2o2));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p2 = constrain(tap, row, sec_strength, sec_damping);
        tap = _mm256_set_epi64x(*(uint64_t*)(in - s2o2),
            *(uint64_t*)(in + CDEF_BSTRIDE - s2o2),
            *(uint64_t*)(in + 2 * CDEF_BSTRIDE - s2o2),
            *(uint64_t*)(in + 3 * CDEF_BSTRIDE - s2o2));

#else
        tap = v256_from_v64(v64_load_unaligned(&in[0 * CDEF_BSTRIDE + s1o2]),
            v64_load_unaligned(&in[1 * CDEF_BSTRIDE + s1o2]),
            v64_load_unaligned(&in[2 * CDEF_BSTRIDE + s1o2]),
            v64_load_unaligned(&in[3 * CDEF_BSTRIDE + s1o2]));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p0 = constrain(tap, row, sec_strength, sec_damping);
        tap = v256_from_v64(v64_load_unaligned(&in[0 * CDEF_BSTRIDE - s1o2]),
            v64_load_unaligned(&in[1 * CDEF_BSTRIDE - s1o2]),
            v64_load_unaligned(&in[2 * CDEF_BSTRIDE - s1o2]),
            v64_load_unaligned(&in[3 * CDEF_BSTRIDE - s1o2]));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p1 = constrain(tap, row, sec_strength, sec_damping);
        tap = v256_from_v64(v64_load_unaligned(&in[0 * CDEF_BSTRIDE + s2o2]),
            v64_load_unaligned(&in[1 * CDEF_BSTRIDE + s2o2]),
            v64_load_unaligned(&in[2 * CDEF_BSTRIDE + s2o2]),
            v64_load_unaligned(&in[3 * CDEF_BSTRIDE + s2o2]));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p2 = constrain(tap, row, sec_strength, sec_damping);
        tap = v256_from_v64(v64_load_unaligned(&in[0 * CDEF_BSTRIDE - s2o2]),
            v64_load_unaligned(&in[1 * CDEF_BSTRIDE - s2o2]),
            v64_load_unaligned(&in[2 * CDEF_BSTRIDE - s2o2]),
            v64_load_unaligned(&in[3 * CDEF_BSTRIDE - s2o2]));
#endif
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p3 = constrain(tap, row, sec_strength, sec_damping);

        // sum += sec_taps[1] * (p0 + p1 + p2 + p3)
        p0 = v128_add_8(p0, p1);
        p2 = v128_add_8(p2, p3);

        sum = v256_add_16(sum, v256_madd_us8(v256_dup_8(sec_taps[1]),
            v256_from_v128(v128_ziphi_8(p0, p2),
                v128_ziplo_8(p0, p2))));
    }

    // res = row + ((sum - (sum < 0) + 8) >> 4)
    sum = v256_add_16(sum, v256_cmplt_s16(sum, v256_zero()));
    res = v256_add_16(sum, v256_dup_16(8));
    res = v256_shr_n_s16(res, 4);
    res = v256_add_16(row, res);
    res = v256_min_s16(v256_max_s16(res, min), max);
    res = v256_pack_s16_u8(res, res);

    p0 = v256_low_v128(res);
    u32_store_aligned(&dst[0 * dstride], v64_high_u32(v128_high_v64(p0)));
    u32_store_aligned(&dst[1 * dstride], v64_low_u32(v128_high_v64(p0)));
    u32_store_aligned(&dst[2 * dstride], v64_high_u32(v128_low_v64(p0)));
    u32_store_aligned(&dst[3 * dstride], v64_low_u32(v128_low_v64(p0)));
}

void SIMD_FUNC(cdef_filter_block_8x8_8)(uint8_t *dst, int32_t dstride,
    const uint16_t *in, int32_t pri_strength,
    int32_t sec_strength, int32_t dir,
    int32_t pri_damping, int32_t sec_damping,
    /*AOM_UNUSED*/ int32_t max_unused,
    int32_t coeff_shift) {
    (void)max_unused;
    int32_t i;
    v128 p0, p1, p2, p3;
    v256 sum, row, res, tap;
    v256 max, min, large = v256_dup_16(CDEF_VERY_LARGE);
    int32_t po1 = cdef_directions[dir][0];
    int32_t po2 = cdef_directions[dir][1];
    int32_t s1o1 = cdef_directions[(dir + 2) & 7][0];
    int32_t s1o2 = cdef_directions[(dir + 2) & 7][1];
    int32_t s2o1 = cdef_directions[(dir + 6) & 7][0];
    int32_t s2o2 = cdef_directions[(dir + 6) & 7][1];

    const int32_t *pri_taps = cdef_pri_taps[(pri_strength >> coeff_shift) & 1];
    const int32_t *sec_taps = cdef_sec_taps[(pri_strength >> coeff_shift) & 1];
#if CDEF_AVX_OPT
    v256 pri_taps_0 = v256_dup_8(pri_taps[0]);
    v256 pri_taps_1 = v256_dup_8(pri_taps[1]);
    v256 sec_taps_0 = v256_dup_8(sec_taps[0]);
    v256 sec_taps_1 = v256_dup_8(sec_taps[1]);
    v256 duplicate_8 = v256_dup_16(8);

#endif
    if (pri_strength)
        pri_damping = AOMMAX(0, pri_damping - get_msb(pri_strength));
    if (sec_strength)
        sec_damping = AOMMAX(0, sec_damping - get_msb(sec_strength));
    for (i = 0; i < 8; i += 2) {
        sum = v256_zero();
#if CDEF_AVX_OPT
        row = _mm256_insertf128_si256(_mm256_castsi128_si256(
            _mm_loadu_si128((__m128i *)(in + (i + 1) * CDEF_BSTRIDE))),
            _mm_loadu_si128((__m128i *)(in + i * CDEF_BSTRIDE)), 0x1);
#else
        row = v256_from_v128(v128_load_aligned(&in[i * CDEF_BSTRIDE]),
            v128_load_aligned(&in[(i + 1) * CDEF_BSTRIDE]));
#endif

        max = min = row;
        // Primary near taps
#if CDEF_AVX_OPT
        tap = _mm256_insertf128_si256(_mm256_castsi128_si256(
            _mm_loadu_si128((__m128i *)(in + (i + 1) * CDEF_BSTRIDE + po1))),
            _mm_loadu_si128((__m128i *)(in + i * CDEF_BSTRIDE + po1)), 0x1);
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p0 = constrain(tap, row, pri_strength, pri_damping);
        tap = _mm256_insertf128_si256(_mm256_castsi128_si256(
            _mm_loadu_si128((__m128i *)(in + (i + 1) * CDEF_BSTRIDE - po1))),
            _mm_loadu_si128((__m128i *)(in + i * CDEF_BSTRIDE - po1)), 0x1);
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p1 = constrain(tap, row, pri_strength, pri_damping);

        // sum += pri_taps[0] * (p0 + p1)
        sum = v256_add_16(sum, v256_madd_us8(pri_taps_0,
#else
        tap =
            v256_from_v128(v128_load_unaligned(&in[i * CDEF_BSTRIDE + po1]),
                v128_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE + po1]));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p0 = constrain(tap, row, pri_strength, pri_damping);
        tap =
            v256_from_v128(v128_load_unaligned(&in[i * CDEF_BSTRIDE - po1]),
                v128_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE - po1]));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p1 = constrain(tap, row, pri_strength, pri_damping);

        // sum += pri_taps[0] * (p0 + p1)
        sum = v256_add_16(sum, v256_madd_us8(v256_dup_8(pri_taps[0]),
#endif
            v256_from_v128(v128_ziphi_8(p0, p1),
                v128_ziplo_8(p0, p1))));

        // Primary far taps
#if CDEF_AVX_OPT
        tap = _mm256_insertf128_si256(_mm256_castsi128_si256(
            _mm_loadu_si128((__m128i *)(in + (i + 1) * CDEF_BSTRIDE + po2))),
            _mm_loadu_si128((__m128i *)(in + i * CDEF_BSTRIDE + po2)), 0x1);
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p0 = constrain(tap, row, pri_strength, pri_damping);
        tap = _mm256_insertf128_si256(_mm256_castsi128_si256(
            _mm_loadu_si128((__m128i *)(in + (i + 1) * CDEF_BSTRIDE - po2))),
            _mm_loadu_si128((__m128i *)(in + i * CDEF_BSTRIDE - po2)), 0x1);
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p1 = constrain(tap, row, pri_strength, pri_damping);
        // sum += pri_taps[1] * (p0 + p1)
        sum = v256_add_16(sum, v256_madd_us8(pri_taps_1,
#else
        tap =
            v256_from_v128(v128_load_unaligned(&in[i * CDEF_BSTRIDE + po2]),
                v128_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE + po2]));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p0 = constrain(tap, row, pri_strength, pri_damping);
        tap =
            v256_from_v128(v128_load_unaligned(&in[i * CDEF_BSTRIDE - po2]),
                v128_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE - po2]));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p1 = constrain(tap, row, pri_strength, pri_damping);

        // sum += pri_taps[1] * (p0 + p1)
        sum = v256_add_16(sum, v256_madd_us8(v256_dup_8(pri_taps[1]),
#endif
            v256_from_v128(v128_ziphi_8(p0, p1),
                v128_ziplo_8(p0, p1))));

        // Secondary near taps
#if CDEF_AVX_OPT
        tap = _mm256_insertf128_si256(_mm256_castsi128_si256(
            _mm_loadu_si128((__m128i *)(in + (i + 1) * CDEF_BSTRIDE + s1o1))),
            _mm_loadu_si128((__m128i *)(in + i * CDEF_BSTRIDE + s1o1)), 0x1);
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p0 = constrain(tap, row, sec_strength, sec_damping);
        tap = _mm256_insertf128_si256(_mm256_castsi128_si256(
            _mm_loadu_si128((__m128i *)(in + (i + 1) * CDEF_BSTRIDE - s1o1))),
            _mm_loadu_si128((__m128i *)(in + i * CDEF_BSTRIDE - s1o1)), 0x1);
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p1 = constrain(tap, row, sec_strength, sec_damping);
        tap = _mm256_insertf128_si256(_mm256_castsi128_si256(
            _mm_loadu_si128((__m128i *)(in + (i + 1) * CDEF_BSTRIDE + s2o1))),
            _mm_loadu_si128((__m128i *)(in + i * CDEF_BSTRIDE + s2o1)), 0x1);
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p2 = constrain(tap, row, sec_strength, sec_damping);
        tap = _mm256_insertf128_si256(_mm256_castsi128_si256(
            _mm_loadu_si128((__m128i *)(in + (i + 1) * CDEF_BSTRIDE - s2o1))),
            _mm_loadu_si128((__m128i *)(in + i * CDEF_BSTRIDE - s2o1)), 0x1);
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p3 = constrain(tap, row, sec_strength, sec_damping);

        // sum += sec_taps[0] * (p0 + p1 + p2 + p3)
        p0 = v128_add_8(p0, p1);
        p2 = v128_add_8(p2, p3);
        sum = v256_add_16(sum, v256_madd_us8(sec_taps_0,
#else
        tap =
            v256_from_v128(v128_load_unaligned(&in[i * CDEF_BSTRIDE + s1o1]),
                v128_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE + s1o1]));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p0 = constrain(tap, row, sec_strength, sec_damping);
        tap =
            v256_from_v128(v128_load_unaligned(&in[i * CDEF_BSTRIDE - s1o1]),
                v128_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE - s1o1]));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p1 = constrain(tap, row, sec_strength, sec_damping);
        tap =
            v256_from_v128(v128_load_unaligned(&in[i * CDEF_BSTRIDE + s2o1]),
                v128_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE + s2o1]));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p2 = constrain(tap, row, sec_strength, sec_damping);
        tap =
            v256_from_v128(v128_load_unaligned(&in[i * CDEF_BSTRIDE - s2o1]),
                v128_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE - s2o1]));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p3 = constrain(tap, row, sec_strength, sec_damping);

        // sum += sec_taps[0] * (p0 + p1 + p2 + p3)
        p0 = v128_add_8(p0, p1);
        p2 = v128_add_8(p2, p3);
        sum = v256_add_16(sum, v256_madd_us8(v256_dup_8(sec_taps[0]),
#endif
            v256_from_v128(v128_ziphi_8(p0, p2),
                v128_ziplo_8(p0, p2))));

        // Secondary far taps
#if CDEF_AVX_OPT
        tap = _mm256_insertf128_si256(_mm256_castsi128_si256(
            _mm_loadu_si128((__m128i *)(in + (i + 1) * CDEF_BSTRIDE + s1o2))),
            _mm_loadu_si128((__m128i *)(in + i * CDEF_BSTRIDE + s1o2)), 0x1);
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p0 = constrain(tap, row, sec_strength, sec_damping);
        tap = _mm256_insertf128_si256(_mm256_castsi128_si256(
            _mm_loadu_si128((__m128i *)(in + (i + 1) * CDEF_BSTRIDE - s1o2))),
            _mm_loadu_si128((__m128i *)(in + i * CDEF_BSTRIDE - s1o2)), 0x1);
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p1 = constrain(tap, row, sec_strength, sec_damping);
        tap = _mm256_insertf128_si256(_mm256_castsi128_si256(
            _mm_loadu_si128((__m128i *)(in + (i + 1) * CDEF_BSTRIDE + s2o2))),
            _mm_loadu_si128((__m128i *)(in + i * CDEF_BSTRIDE + s2o2)), 0x1);
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p2 = constrain(tap, row, sec_strength, sec_damping);
        tap = _mm256_insertf128_si256(_mm256_castsi128_si256(
            _mm_loadu_si128((__m128i *)(in + (i + 1) * CDEF_BSTRIDE - s2o2))),
            _mm_loadu_si128((__m128i *)(in + i * CDEF_BSTRIDE - s2o2)), 0x1);
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p3 = constrain(tap, row, sec_strength, sec_damping);

        // sum += sec_taps[1] * (p0 + p1 + p2 + p3)
        p0 = v128_add_8(p0, p1);
        p2 = v128_add_8(p2, p3);
        sum = v256_add_16(sum, v256_madd_us8(sec_taps_1,
            v256_from_v128(v128_ziphi_8(p0, p2),
                v128_ziplo_8(p0, p2))));

        // res = row + ((sum - (sum < 0) + 8) >> 4)
        sum = v256_add_16(sum, v256_cmplt_s16(sum, v256_zero()));
        res = v256_add_16(sum, duplicate_8);
        res = v256_shr_n_s16(res, 4);
        res = v256_add_16(row, res);
        res = v256_min_s16(v256_max_s16(res, min), max);
        res = v256_pack_s16_u8(res, res);
        *(uint64_t*)(dst + i * dstride) = _mm256_extract_epi64(res, 1);
        *(uint64_t*)(dst + (i + 1) * dstride) = _mm256_extract_epi64(res, 0);
#else
        tap =
            v256_from_v128(v128_load_unaligned(&in[i * CDEF_BSTRIDE + s1o2]),
                v128_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE + s1o2]));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p0 = constrain(tap, row, sec_strength, sec_damping);
        tap =
            v256_from_v128(v128_load_unaligned(&in[i * CDEF_BSTRIDE - s1o2]),
                v128_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE - s1o2]));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p1 = constrain(tap, row, sec_strength, sec_damping);
        tap =
            v256_from_v128(v128_load_unaligned(&in[i * CDEF_BSTRIDE + s2o2]),
                v128_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE + s2o2]));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p2 = constrain(tap, row, sec_strength, sec_damping);
        tap =
            v256_from_v128(v128_load_unaligned(&in[i * CDEF_BSTRIDE - s2o2]),
                v128_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE - s2o2]));
        max = v256_max_s16(max, v256_andn(tap, v256_cmpeq_16(tap, large)));
        min = v256_min_s16(min, tap);
        p3 = constrain(tap, row, sec_strength, sec_damping);

        // sum += sec_taps[1] * (p0 + p1 + p2 + p3)
        p0 = v128_add_8(p0, p1);
        p2 = v128_add_8(p2, p3);
        sum = v256_add_16(sum, v256_madd_us8(v256_dup_8(sec_taps[1]),
            v256_from_v128(v128_ziphi_8(p0, p2),
                v128_ziplo_8(p0, p2))));

        // res = row + ((sum - (sum < 0) + 8) >> 4)
        sum = v256_add_16(sum, v256_cmplt_s16(sum, v256_zero()));
        res = v256_add_16(sum, v256_dup_16(8));

        res = v256_shr_n_s16(res, 4);
        res = v256_add_16(row, res);
        res = v256_min_s16(v256_max_s16(res, min), max);
        res = v256_pack_s16_u8(res, res);

        p0 = v256_low_v128(res);
        v64_store_aligned(&dst[i * dstride], v128_high_v64(p0));
        v64_store_aligned(&dst[(i + 1) * dstride], v128_low_v64(p0));
#endif
    }
}

void SIMD_FUNC(cdef_filter_block_4x4_16)(uint16_t *dst, int32_t dstride,
    const uint16_t *in, int32_t pri_strength,
    int32_t sec_strength, int32_t dir,
    int32_t pri_damping, int32_t sec_damping,
    /*AOM_UNUSED*/ int32_t max_unused,
    int32_t coeff_shift) {
    (void)max_unused;
#if ! CDEF_AVX_OPT
    int32_t i;
#endif
    v256 p0, p1, p2, p3, sum, row, res;
    v256 max, min, large = v256_dup_16(CDEF_VERY_LARGE);
    int32_t po1 = cdef_directions[dir][0];
    int32_t po2 = cdef_directions[dir][1];
    int32_t s1o1 = cdef_directions[(dir + 2) & 7][0];
    int32_t s1o2 = cdef_directions[(dir + 2) & 7][1];
    int32_t s2o1 = cdef_directions[(dir + 6) & 7][0];
    int32_t s2o2 = cdef_directions[(dir + 6) & 7][1];

    const int32_t *pri_taps = cdef_pri_taps[(pri_strength >> coeff_shift) & 1];
    const int32_t *sec_taps = cdef_sec_taps[(pri_strength >> coeff_shift) & 1];
#if  CDEF_AVX_OPT
    v256 pri_strength_256 = v256_dup_16(pri_strength);
    v256 sec_strength_256 = v256_dup_16(sec_strength);

#endif
    if (pri_strength)
        pri_damping = AOMMAX(0, pri_damping - get_msb(pri_strength));
    if (sec_strength)
        sec_damping = AOMMAX(0, sec_damping - get_msb(sec_strength));
#if  CDEF_AVX_OPT
    sum = v256_zero();
    row = _mm256_set_epi64x(*(uint64_t*)(in),
        *(uint64_t*)(in + 1 * CDEF_BSTRIDE),
        *(uint64_t*)(in + 2 * CDEF_BSTRIDE),
        *(uint64_t*)(in + 3 * CDEF_BSTRIDE));
    min = max = row;

    // Primary near taps
    p0 = _mm256_set_epi64x(*(uint64_t*)(in + po1),
        *(uint64_t*)(in + 1 * CDEF_BSTRIDE + po1),
        *(uint64_t*)(in + 2 * CDEF_BSTRIDE + po1),
        *(uint64_t*)(in + 3 * CDEF_BSTRIDE + po1));
    p1 = _mm256_set_epi64x(*(uint64_t*)(in - po1),
        *(uint64_t*)(in + 1 * CDEF_BSTRIDE - po1),
        *(uint64_t*)(in + 2 * CDEF_BSTRIDE - po1),
        *(uint64_t*)(in + 3 * CDEF_BSTRIDE - po1));

    max =
        v256_max_s16(v256_max_s16(max, v256_andn(p0, v256_cmpeq_16(p0, large))),
            v256_andn(p1, v256_cmpeq_16(p1, large)));
    min = v256_min_s16(v256_min_s16(min, p0), p1);
    p0 = constrain16(p0, row, pri_strength_256, pri_damping);
    p1 = constrain16(p1, row, pri_strength_256, pri_damping);

    // sum += pri_taps[0] * (p0 + p1)
    sum = v256_add_16(
        sum, v256_mullo_s16(v256_dup_16(pri_taps[0]), v256_add_16(p0, p1)));

    // Primary far taps
    p0 = _mm256_set_epi64x(*(uint64_t*)(in + po2),
        *(uint64_t*)(in + 1 * CDEF_BSTRIDE + po2),
        *(uint64_t*)(in + 2 * CDEF_BSTRIDE + po2),
        *(uint64_t*)(in + 3 * CDEF_BSTRIDE + po2));
    p1 = _mm256_set_epi64x(*(uint64_t*)(in - po2),
        *(uint64_t*)(in + 1 * CDEF_BSTRIDE - po2),
        *(uint64_t*)(in + 2 * CDEF_BSTRIDE - po2),
        *(uint64_t*)(in + 3 * CDEF_BSTRIDE - po2));
    max =
        v256_max_s16(v256_max_s16(max, v256_andn(p0, v256_cmpeq_16(p0, large))),
            v256_andn(p1, v256_cmpeq_16(p1, large)));
    min = v256_min_s16(v256_min_s16(min, p0), p1);
    p0 = constrain16(p0, row, pri_strength_256, pri_damping);
    p1 = constrain16(p1, row, pri_strength_256, pri_damping);

    // sum += pri_taps[1] * (p0 + p1)
    sum = v256_add_16(
        sum, v256_mullo_s16(v256_dup_16(pri_taps[1]), v256_add_16(p0, p1)));

    // Secondary near taps
    p0 = _mm256_set_epi64x(*(uint64_t*)(in + s1o1),
        *(uint64_t*)(in + 1 * CDEF_BSTRIDE + s1o1),
        *(uint64_t*)(in + 2 * CDEF_BSTRIDE + s1o1),
        *(uint64_t*)(in + 3 * CDEF_BSTRIDE + s1o1));
    p1 = _mm256_set_epi64x(*(uint64_t*)(in - s1o1),
        *(uint64_t*)(in + 1 * CDEF_BSTRIDE - s1o1),
        *(uint64_t*)(in + 2 * CDEF_BSTRIDE - s1o1),
        *(uint64_t*)(in + 3 * CDEF_BSTRIDE - s1o1));
    p2 = _mm256_set_epi64x(*(uint64_t*)(in + s2o1),
        *(uint64_t*)(in + 1 * CDEF_BSTRIDE + s2o1),
        *(uint64_t*)(in + 2 * CDEF_BSTRIDE + s2o1),
        *(uint64_t*)(in + 3 * CDEF_BSTRIDE + s2o1));
    p3 = _mm256_set_epi64x(*(uint64_t*)(in - s2o1),
        *(uint64_t*)(in + 1 * CDEF_BSTRIDE - s2o1),
        *(uint64_t*)(in + 2 * CDEF_BSTRIDE - s2o1),
        *(uint64_t*)(in + 3 * CDEF_BSTRIDE - s2o1));
    max =
        v256_max_s16(v256_max_s16(max, v256_andn(p0, v256_cmpeq_16(p0, large))),
            v256_andn(p1, v256_cmpeq_16(p1, large)));
    max =
        v256_max_s16(v256_max_s16(max, v256_andn(p2, v256_cmpeq_16(p2, large))),
            v256_andn(p3, v256_cmpeq_16(p3, large)));
    min = v256_min_s16(
        v256_min_s16(v256_min_s16(v256_min_s16(min, p0), p1), p2), p3);
    p0 = constrain16(p0, row, sec_strength_256, sec_damping);
    p1 = constrain16(p1, row, sec_strength_256, sec_damping);
    p2 = constrain16(p2, row, sec_strength_256, sec_damping);
    p3 = constrain16(p3, row, sec_strength_256, sec_damping);

    // sum += sec_taps[0] * (p0 + p1 + p2 + p3)
    sum = v256_add_16(sum, v256_mullo_s16(v256_dup_16(sec_taps[0]),
        v256_add_16(v256_add_16(p0, p1),
            v256_add_16(p2, p3))));

    // Secondary far taps
    p0 = _mm256_set_epi64x(*(uint64_t*)(in + s1o2),
        *(uint64_t*)(in + 1 * CDEF_BSTRIDE + s1o2),
        *(uint64_t*)(in + 2 * CDEF_BSTRIDE + s1o2),
        *(uint64_t*)(in + 3 * CDEF_BSTRIDE + s1o2));
    p1 = _mm256_set_epi64x(*(uint64_t*)(in - s1o2),
        *(uint64_t*)(in + 1 * CDEF_BSTRIDE - s1o2),
        *(uint64_t*)(in + 2 * CDEF_BSTRIDE - s1o2),
        *(uint64_t*)(in + 3 * CDEF_BSTRIDE - s1o2));
    p2 = _mm256_set_epi64x(*(uint64_t*)(in + s2o2),
        *(uint64_t*)(in + 1 * CDEF_BSTRIDE + s2o2),
        *(uint64_t*)(in + 2 * CDEF_BSTRIDE + s2o2),
        *(uint64_t*)(in + 3 * CDEF_BSTRIDE + s2o2));
    p3 = _mm256_set_epi64x(*(uint64_t*)(in - s2o2),
        *(uint64_t*)(in + 1 * CDEF_BSTRIDE - s2o2),
        *(uint64_t*)(in + 2 * CDEF_BSTRIDE - s2o2),
        *(uint64_t*)(in + 3 * CDEF_BSTRIDE - s2o2));
    max =
        v256_max_s16(v256_max_s16(max, v256_andn(p0, v256_cmpeq_16(p0, large))),
            v256_andn(p1, v256_cmpeq_16(p1, large)));
    max =
        v256_max_s16(v256_max_s16(max, v256_andn(p2, v256_cmpeq_16(p2, large))),
            v256_andn(p3, v256_cmpeq_16(p3, large)));
    min = v256_min_s16(
        v256_min_s16(v256_min_s16(v256_min_s16(min, p0), p1), p2), p3);
    p0 = constrain16(p0, row, sec_strength_256, sec_damping);
    p1 = constrain16(p1, row, sec_strength_256, sec_damping);
    p2 = constrain16(p2, row, sec_strength_256, sec_damping);
    p3 = constrain16(p3, row, sec_strength_256, sec_damping);

    // sum += sec_taps[1] * (p0 + p1 + p2 + p3)
    sum = v256_add_16(sum, v256_mullo_s16(v256_dup_16(sec_taps[1]),
        v256_add_16(v256_add_16(p0, p1),
            v256_add_16(p2, p3))));

#else
    for (i = 0; i < 4; i += 4) {
        sum = v256_zero();
        row = v256_from_v64(v64_load_aligned(&in[i * CDEF_BSTRIDE]),
            v64_load_aligned(&in[(i + 1) * CDEF_BSTRIDE]),
            v64_load_aligned(&in[(i + 2) * CDEF_BSTRIDE]),
            v64_load_aligned(&in[(i + 3) * CDEF_BSTRIDE]));
        min = max = row;

        // Primary near taps
        p0 = v256_from_v64(v64_load_unaligned(&in[i * CDEF_BSTRIDE + po1]),
            v64_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE + po1]),
            v64_load_unaligned(&in[(i + 2) * CDEF_BSTRIDE + po1]),
            v64_load_unaligned(&in[(i + 3) * CDEF_BSTRIDE + po1]));
        p1 = v256_from_v64(v64_load_unaligned(&in[i * CDEF_BSTRIDE - po1]),
            v64_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE - po1]),
            v64_load_unaligned(&in[(i + 2) * CDEF_BSTRIDE - po1]),
            v64_load_unaligned(&in[(i + 3) * CDEF_BSTRIDE - po1]));
        max =
            v256_max_s16(v256_max_s16(max, v256_andn(p0, v256_cmpeq_16(p0, large))),
                v256_andn(p1, v256_cmpeq_16(p1, large)));
        min = v256_min_s16(v256_min_s16(min, p0), p1);
        p0 = constrain16(p0, row, pri_strength, pri_damping);
        p1 = constrain16(p1, row, pri_strength, pri_damping);

        // sum += pri_taps[0] * (p0 + p1)
        sum = v256_add_16(
            sum, v256_mullo_s16(v256_dup_16(pri_taps[0]), v256_add_16(p0, p1)));

        // Primary far taps
        p0 = v256_from_v64(v64_load_unaligned(&in[i * CDEF_BSTRIDE + po2]),
            v64_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE + po2]),
            v64_load_unaligned(&in[(i + 2) * CDEF_BSTRIDE + po2]),
            v64_load_unaligned(&in[(i + 3) * CDEF_BSTRIDE + po2]));
        p1 = v256_from_v64(v64_load_unaligned(&in[i * CDEF_BSTRIDE - po2]),
            v64_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE - po2]),
            v64_load_unaligned(&in[(i + 2) * CDEF_BSTRIDE - po2]),
            v64_load_unaligned(&in[(i + 3) * CDEF_BSTRIDE - po2]));
        max =
            v256_max_s16(v256_max_s16(max, v256_andn(p0, v256_cmpeq_16(p0, large))),
                v256_andn(p1, v256_cmpeq_16(p1, large)));
        min = v256_min_s16(v256_min_s16(min, p0), p1);
        p0 = constrain16(p0, row, pri_strength, pri_damping);
        p1 = constrain16(p1, row, pri_strength, pri_damping);

        // sum += pri_taps[1] * (p0 + p1)
        sum = v256_add_16(
            sum, v256_mullo_s16(v256_dup_16(pri_taps[1]), v256_add_16(p0, p1)));

        // Secondary near taps
        p0 = v256_from_v64(v64_load_unaligned(&in[i * CDEF_BSTRIDE + s1o1]),
            v64_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE + s1o1]),
            v64_load_unaligned(&in[(i + 2) * CDEF_BSTRIDE + s1o1]),
            v64_load_unaligned(&in[(i + 3) * CDEF_BSTRIDE + s1o1]));
        p1 = v256_from_v64(v64_load_unaligned(&in[i * CDEF_BSTRIDE - s1o1]),
            v64_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE - s1o1]),
            v64_load_unaligned(&in[(i + 2) * CDEF_BSTRIDE - s1o1]),
            v64_load_unaligned(&in[(i + 3) * CDEF_BSTRIDE - s1o1]));
        p2 = v256_from_v64(v64_load_unaligned(&in[i * CDEF_BSTRIDE + s2o1]),
            v64_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE + s2o1]),
            v64_load_unaligned(&in[(i + 2) * CDEF_BSTRIDE + s2o1]),
            v64_load_unaligned(&in[(i + 3) * CDEF_BSTRIDE + s2o1]));
        p3 = v256_from_v64(v64_load_unaligned(&in[i * CDEF_BSTRIDE - s2o1]),
            v64_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE - s2o1]),
            v64_load_unaligned(&in[(i + 2) * CDEF_BSTRIDE - s2o1]),
            v64_load_unaligned(&in[(i + 3) * CDEF_BSTRIDE - s2o1]));
        max =
            v256_max_s16(v256_max_s16(max, v256_andn(p0, v256_cmpeq_16(p0, large))),
                v256_andn(p1, v256_cmpeq_16(p1, large)));
        max =
            v256_max_s16(v256_max_s16(max, v256_andn(p2, v256_cmpeq_16(p2, large))),
                v256_andn(p3, v256_cmpeq_16(p3, large)));
        min = v256_min_s16(
            v256_min_s16(v256_min_s16(v256_min_s16(min, p0), p1), p2), p3);
        p0 = constrain16(p0, row, sec_strength, sec_damping);
        p1 = constrain16(p1, row, sec_strength, sec_damping);
        p2 = constrain16(p2, row, sec_strength, sec_damping);
        p3 = constrain16(p3, row, sec_strength, sec_damping);

        // sum += sec_taps[0] * (p0 + p1 + p2 + p3)
        sum = v256_add_16(sum, v256_mullo_s16(v256_dup_16(sec_taps[0]),
            v256_add_16(v256_add_16(p0, p1),
                v256_add_16(p2, p3))));
        // Secondary far taps
        p0 = v256_from_v64(v64_load_unaligned(&in[i * CDEF_BSTRIDE + s1o2]),
            v64_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE + s1o2]),
            v64_load_unaligned(&in[(i + 2) * CDEF_BSTRIDE + s1o2]),
            v64_load_unaligned(&in[(i + 3) * CDEF_BSTRIDE + s1o2]));
        p1 = v256_from_v64(v64_load_unaligned(&in[i * CDEF_BSTRIDE - s1o2]),
            v64_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE - s1o2]),
            v64_load_unaligned(&in[(i + 2) * CDEF_BSTRIDE - s1o2]),
            v64_load_unaligned(&in[(i + 3) * CDEF_BSTRIDE - s1o2]));
        p2 = v256_from_v64(v64_load_unaligned(&in[i * CDEF_BSTRIDE + s2o2]),
            v64_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE + s2o2]),
            v64_load_unaligned(&in[(i + 2) * CDEF_BSTRIDE + s2o2]),
            v64_load_unaligned(&in[(i + 3) * CDEF_BSTRIDE + s2o2]));
        p3 = v256_from_v64(v64_load_unaligned(&in[i * CDEF_BSTRIDE - s2o2]),
            v64_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE - s2o2]),
            v64_load_unaligned(&in[(i + 2) * CDEF_BSTRIDE - s2o2]),
            v64_load_unaligned(&in[(i + 3) * CDEF_BSTRIDE - s2o2]));
        max =
            v256_max_s16(v256_max_s16(max, v256_andn(p0, v256_cmpeq_16(p0, large))),
                v256_andn(p1, v256_cmpeq_16(p1, large)));
        max =
            v256_max_s16(v256_max_s16(max, v256_andn(p2, v256_cmpeq_16(p2, large))),
                v256_andn(p3, v256_cmpeq_16(p3, large)));
        min = v256_min_s16(
            v256_min_s16(v256_min_s16(v256_min_s16(min, p0), p1), p2), p3);
        p0 = constrain16(p0, row, sec_strength, sec_damping);
        p1 = constrain16(p1, row, sec_strength, sec_damping);
        p2 = constrain16(p2, row, sec_strength, sec_damping);
        p3 = constrain16(p3, row, sec_strength, sec_damping);

        // sum += sec_taps[1] * (p0 + p1 + p2 + p3)
        sum = v256_add_16(sum, v256_mullo_s16(v256_dup_16(sec_taps[1]),
            v256_add_16(v256_add_16(p0, p1),
                v256_add_16(p2, p3))));
#endif

#if  CDEF_AVX_OPT
        // res = row + ((sum - (sum < 0) + 8) >> 4)
        sum = v256_add_16(sum, v256_cmplt_s16(sum, v256_zero()));
        res = v256_add_16(sum, v256_dup_16(8));
        res = v256_shr_n_s16(res, 4);
        res = v256_add_16(row, res);
        res = v256_min_s16(v256_max_s16(res, min), max);
        *(uint64_t*)(dst) = _mm256_extract_epi64(res, 3);
        *(uint64_t*)(dst + 1 * dstride) = _mm256_extract_epi64(res, 2);
        *(uint64_t*)(dst + 2 * dstride) = _mm256_extract_epi64(res, 1);
        *(uint64_t*)(dst + 3 * dstride) = _mm256_extract_epi64(res, 0);

#else
        // res = row + ((sum - (sum < 0) + 8) >> 4)
        sum = v256_add_16(sum, v256_cmplt_s16(sum, v256_zero()));
        res = v256_add_16(sum, v256_dup_16(8));
        res = v256_shr_n_s16(res, 4);
        res = v256_add_16(row, res);
        res = v256_min_s16(v256_max_s16(res, min), max);

        v64_store_aligned(&dst[i * dstride], v128_high_v64(v256_high_v128(res)));
        v64_store_aligned(&dst[(i + 1) * dstride],
            v128_low_v64(v256_high_v128(res)));
        v64_store_aligned(&dst[(i + 2) * dstride],
            v128_high_v64(v256_low_v128(res)));
        v64_store_aligned(&dst[(i + 3) * dstride],
            v128_low_v64(v256_low_v128(res)));
    }
#endif
}

void SIMD_FUNC(cdef_filter_block_8x8_16)(uint16_t *dst, int32_t dstride,
    const uint16_t *in, int32_t pri_strength,
    int32_t sec_strength, int32_t dir,
    int32_t pri_damping, int32_t sec_damping,
    /*AOM_UNUSED*/ int32_t max_unused,
    int32_t coeff_shift) {
    (void)max_unused;
    int32_t i;
    v256 sum, p0, p1, p2, p3, row, res;
    v256 max, min, large = v256_dup_16(CDEF_VERY_LARGE);
    int32_t po1 = cdef_directions[dir][0];
    int32_t po2 = cdef_directions[dir][1];
    int32_t s1o1 = cdef_directions[(dir + 2) & 7][0];
    int32_t s1o2 = cdef_directions[(dir + 2) & 7][1];
    int32_t s2o1 = cdef_directions[(dir + 6) & 7][0];
    int32_t s2o2 = cdef_directions[(dir + 6) & 7][1];
    //SSE CHKN
    const int32_t *pri_taps = cdef_pri_taps[(pri_strength >> coeff_shift) & 1];
    const int32_t *sec_taps = cdef_sec_taps[(pri_strength >> coeff_shift) & 1];
#if CDEF_AVX_OPT
    v256 pri_taps_0 = v256_dup_16(pri_taps[0]);
    v256 pri_taps_1 = v256_dup_16(pri_taps[1]);
    v256 sec_taps_0 = v256_dup_16(sec_taps[0]);
    v256 sec_taps_1 = v256_dup_16(sec_taps[1]);
    v256 duplicate_8 = v256_dup_16(8);
    v256 pri_strength_256 = v256_dup_16(pri_strength);
    v256 sec_strength_256 = v256_dup_16(sec_strength);
#endif
    if (pri_strength)
        pri_damping = AOMMAX(0, pri_damping - get_msb(pri_strength));
    if (sec_strength)
        sec_damping = AOMMAX(0, sec_damping - get_msb(sec_strength));

    for (i = 0; i < 8; i += 2) {
        sum = v256_zero();
#if CDEF_AVX_OPT
        row = _mm256_insertf128_si256(_mm256_castsi128_si256(
            _mm_loadu_si128((__m128i *)(in + (i + 1) * CDEF_BSTRIDE))),
            _mm_loadu_si128((__m128i *)(in + i * CDEF_BSTRIDE)), 0x1);

#else
        row = v256_from_v128(v128_load_aligned(&in[i * CDEF_BSTRIDE]),
            v128_load_aligned(&in[(i + 1) * CDEF_BSTRIDE]));
#endif

        min = max = row;
        // Primary near taps
#if CDEF_AVX_OPT
        p0 = _mm256_insertf128_si256(_mm256_castsi128_si256(
            _mm_loadu_si128((__m128i *)(in + (i + 1) * CDEF_BSTRIDE + po1))),
            _mm_loadu_si128((__m128i *)(in + i * CDEF_BSTRIDE + po1)), 0x1);
        p1 = _mm256_insertf128_si256(_mm256_castsi128_si256(
            _mm_loadu_si128((__m128i *)(in + (i + 1) * CDEF_BSTRIDE - po1))),
            _mm_loadu_si128((__m128i *)(in + i * CDEF_BSTRIDE - po1)), 0x1);
#else
        p0 = v256_from_v128(v128_load_unaligned(&in[i * CDEF_BSTRIDE + po1]),
            v128_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE + po1]));
        p1 = v256_from_v128(v128_load_unaligned(&in[i * CDEF_BSTRIDE - po1]),
            v128_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE - po1]));
#endif
        max =
            v256_max_s16(v256_max_s16(max, v256_andn(p0, v256_cmpeq_16(p0, large))),
                v256_andn(p1, v256_cmpeq_16(p1, large)));
        min = v256_min_s16(v256_min_s16(min, p0), p1);
#if CDEF_AVX_OPT
        p0 = constrain16(p0, row, pri_strength_256, pri_damping);
        p1 = constrain16(p1, row, pri_strength_256, pri_damping);
#else
        p0 = constrain16(p0, row, pri_strength, pri_damping);
        p1 = constrain16(p1, row, pri_strength, pri_damping);
#endif

        // sum += pri_taps[0] * (p0 + p1)
        sum = v256_add_16(
#if CDEF_AVX_OPT
            sum, v256_mullo_s16(pri_taps_0, v256_add_16(p0, p1)));
#else
            sum, v256_mullo_s16(v256_dup_16(pri_taps[0]), v256_add_16(p0, p1)));
#endif

        // Primary far taps
#if CDEF_AVX_OPT
        p0 = _mm256_insertf128_si256(_mm256_castsi128_si256(
            _mm_loadu_si128((__m128i *)(in + (i + 1) * CDEF_BSTRIDE + po2))),
            _mm_loadu_si128((__m128i *)(in + i * CDEF_BSTRIDE + po2)), 0x1);
        p1 = _mm256_insertf128_si256(_mm256_castsi128_si256(
            _mm_loadu_si128((__m128i *)(in + (i + 1) * CDEF_BSTRIDE - po2))),
            _mm_loadu_si128((__m128i *)(in + i * CDEF_BSTRIDE - po2)), 0x1);
#else
        p0 = v256_from_v128(v128_load_unaligned(&in[i * CDEF_BSTRIDE + po2]),
            v128_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE + po2]));
        p1 = v256_from_v128(v128_load_unaligned(&in[i * CDEF_BSTRIDE - po2]),
            v128_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE - po2]));
#endif
        max =
            v256_max_s16(v256_max_s16(max, v256_andn(p0, v256_cmpeq_16(p0, large))),
                v256_andn(p1, v256_cmpeq_16(p1, large)));
        min = v256_min_s16(v256_min_s16(min, p0), p1);
#if CDEF_AVX_OPT
        p0 = constrain16(p0, row, pri_strength_256, pri_damping);
        p1 = constrain16(p1, row, pri_strength_256, pri_damping);
#else
        p0 = constrain16(p0, row, pri_strength, pri_damping);
        p1 = constrain16(p1, row, pri_strength, pri_damping);
#endif

        // sum += pri_taps[1] * (p0 + p1)
        sum = v256_add_16(
#if CDEF_AVX_OPT
            sum, v256_mullo_s16(pri_taps_1, v256_add_16(p0, p1)));
#else
            sum, v256_mullo_s16(v256_dup_16(pri_taps[1]), v256_add_16(p0, p1)));
#endif

        // Secondary near taps
#if CDEF_AVX_OPT
        p0 = _mm256_insertf128_si256(_mm256_castsi128_si256(
            _mm_loadu_si128((__m128i *)(in + (i + 1) * CDEF_BSTRIDE + s1o1))),
            _mm_loadu_si128((__m128i *)(in + i * CDEF_BSTRIDE + s1o1)), 0x1);
        p1 = _mm256_insertf128_si256(_mm256_castsi128_si256(
            _mm_loadu_si128((__m128i *)(in + (i + 1) * CDEF_BSTRIDE - s1o1))),
            _mm_loadu_si128((__m128i *)(in + i * CDEF_BSTRIDE - s1o1)), 0x1);
        p2 = _mm256_insertf128_si256(_mm256_castsi128_si256(
            _mm_loadu_si128((__m128i *)(in + (i + 1) * CDEF_BSTRIDE + s2o1))),
            _mm_loadu_si128((__m128i *)(in + i * CDEF_BSTRIDE + s2o1)), 0x1);
        p3 = _mm256_insertf128_si256(_mm256_castsi128_si256(
            _mm_loadu_si128((__m128i *)(in + (i + 1) * CDEF_BSTRIDE - s2o1))),
            _mm_loadu_si128((__m128i *)(in + i * CDEF_BSTRIDE - s2o1)), 0x1);
#else
        p0 =
            v256_from_v128(v128_load_unaligned(&in[i * CDEF_BSTRIDE + s1o1]),
                v128_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE + s1o1]));
        p1 =
            v256_from_v128(v128_load_unaligned(&in[i * CDEF_BSTRIDE - s1o1]),
                v128_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE - s1o1]));
        p2 =
            v256_from_v128(v128_load_unaligned(&in[i * CDEF_BSTRIDE + s2o1]),
                v128_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE + s2o1]));
        p3 =
            v256_from_v128(v128_load_unaligned(&in[i * CDEF_BSTRIDE - s2o1]),
                v128_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE - s2o1]));
#endif
        max =
            v256_max_s16(v256_max_s16(max, v256_andn(p0, v256_cmpeq_16(p0, large))),
                v256_andn(p1, v256_cmpeq_16(p1, large)));
        max =
            v256_max_s16(v256_max_s16(max, v256_andn(p2, v256_cmpeq_16(p2, large))),
                v256_andn(p3, v256_cmpeq_16(p3, large)));
        min = v256_min_s16(
            v256_min_s16(v256_min_s16(v256_min_s16(min, p0), p1), p2), p3);
#if CDEF_AVX_OPT
        p0 = constrain16(p0, row, sec_strength_256, sec_damping);
        p1 = constrain16(p1, row, sec_strength_256, sec_damping);
        p2 = constrain16(p2, row, sec_strength_256, sec_damping);
        p3 = constrain16(p3, row, sec_strength_256, sec_damping);

#else
        p0 = constrain16(p0, row, sec_strength, sec_damping);
        p1 = constrain16(p1, row, sec_strength, sec_damping);
        p2 = constrain16(p2, row, sec_strength, sec_damping);
        p3 = constrain16(p3, row, sec_strength, sec_damping);
#endif
        // sum += sec_taps[0] * (p0 + p1 + p2 + p3)
#if CDEF_AVX_OPT
        sum = v256_add_16(sum, v256_mullo_s16(sec_taps_0,
#else
        sum = v256_add_16(sum, v256_mullo_s16(v256_dup_16(sec_taps[0]),
#endif
            v256_add_16(v256_add_16(p0, p1),
                v256_add_16(p2, p3))));

        // Secondary far taps
#if CDEF_AVX_OPT
        p0 = _mm256_insertf128_si256(_mm256_castsi128_si256(
            _mm_loadu_si128((__m128i *)(in + (i + 1) * CDEF_BSTRIDE + s1o2))),
            _mm_loadu_si128((__m128i *)(in + i * CDEF_BSTRIDE + s1o2)), 0x1);
        p1 = _mm256_insertf128_si256(_mm256_castsi128_si256(
            _mm_loadu_si128((__m128i *)(in + (i + 1) * CDEF_BSTRIDE - s1o2))),
            _mm_loadu_si128((__m128i *)(in + i * CDEF_BSTRIDE - s1o2)), 0x1);
        p2 = _mm256_insertf128_si256(_mm256_castsi128_si256(
            _mm_loadu_si128((__m128i *)(in + (i + 1)* CDEF_BSTRIDE + s2o2))),
            _mm_loadu_si128((__m128i *)(in + i * CDEF_BSTRIDE + s2o2)), 0x1);
        p3 = _mm256_insertf128_si256(_mm256_castsi128_si256(
            _mm_loadu_si128((__m128i *)(in + (i + 1) * CDEF_BSTRIDE - s2o2))),
            _mm_loadu_si128((__m128i *)(in + i * CDEF_BSTRIDE - s2o2)), 0x1);
#else
        p0 =
            v256_from_v128(v128_load_unaligned(&in[i * CDEF_BSTRIDE + s1o2]),
                v128_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE + s1o2]));
        p1 =
            v256_from_v128(v128_load_unaligned(&in[i * CDEF_BSTRIDE - s1o2]),
                v128_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE - s1o2]));
        p2 =
            v256_from_v128(v128_load_unaligned(&in[i * CDEF_BSTRIDE + s2o2]),
                v128_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE + s2o2]));
        p3 =
            v256_from_v128(v128_load_unaligned(&in[i * CDEF_BSTRIDE - s2o2]),
                v128_load_unaligned(&in[(i + 1) * CDEF_BSTRIDE - s2o2]));
#endif
        max =
            v256_max_s16(v256_max_s16(max, v256_andn(p0, v256_cmpeq_16(p0, large))),
                v256_andn(p1, v256_cmpeq_16(p1, large)));
        max =
            v256_max_s16(v256_max_s16(max, v256_andn(p2, v256_cmpeq_16(p2, large))),
                v256_andn(p3, v256_cmpeq_16(p3, large)));
        min = v256_min_s16(
            v256_min_s16(v256_min_s16(v256_min_s16(min, p0), p1), p2), p3);
#if CDEF_AVX_OPT
        p0 = constrain16(p0, row, sec_strength_256, sec_damping);
        p1 = constrain16(p1, row, sec_strength_256, sec_damping);
        p2 = constrain16(p2, row, sec_strength_256, sec_damping);
        p3 = constrain16(p3, row, sec_strength_256, sec_damping);
#else
        p0 = constrain16(p0, row, sec_strength, sec_damping);
        p1 = constrain16(p1, row, sec_strength, sec_damping);
        p2 = constrain16(p2, row, sec_strength, sec_damping);
        p3 = constrain16(p3, row, sec_strength, sec_damping);
#endif
        // sum += sec_taps[1] * (p0 + p1 + p2 + p3)
#if CDEF_AVX_OPT
        sum = v256_add_16(sum, v256_mullo_s16(sec_taps_1,
#else
        sum = v256_add_16(sum, v256_mullo_s16(v256_dup_16(sec_taps[1]),
#endif
            v256_add_16(v256_add_16(p0, p1),
#if CDEF_AVX_OPT
                v256_add_16(p2, p3))));
#else
                v256_add_16(p2, p3))));
#endif

        // res = row + ((sum - (sum < 0) + 8) >> 4)
        sum = v256_add_16(sum, v256_cmplt_s16(sum, v256_zero()));
#if CDEF_AVX_OPT
        res = v256_add_16(sum, duplicate_8);
#else
        res = v256_add_16(sum, v256_dup_16(8));
#endif
        res = v256_shr_n_s16(res, 4);
        res = v256_add_16(row, res);
        res = v256_min_s16(v256_max_s16(res, min), max);
        v128_store_unaligned(&dst[i * dstride], v256_high_v128(res));
#if CDEF_AVX_OPT
        v128_store_unaligned(&dst[(i + 1) * dstride], _mm256_castsi256_si128(res));
#else
        v128_store_unaligned(&dst[(i + 1) * dstride], v256_low_v128(res));
#endif
    }
}

void SIMD_FUNC(cdef_filter_block)(uint8_t *dst8, uint16_t *dst16, int32_t dstride,
    const uint16_t *in, int32_t pri_strength,
    int32_t sec_strength, int32_t dir, int32_t pri_damping,
    int32_t sec_damping, int32_t bsize, int32_t max,
    int32_t coeff_shift) {
    if (dst8) {
        if (bsize == BLOCK_8X8) {
            SIMD_FUNC(cdef_filter_block_8x8_8)
                (dst8, dstride, in, pri_strength, sec_strength, dir, pri_damping,
                    sec_damping, max, coeff_shift);
        }
        else if (bsize == BLOCK_4X8) {
            SIMD_FUNC(cdef_filter_block_4x4_8)
                (dst8, dstride, in, pri_strength, sec_strength, dir, pri_damping,
                    sec_damping, max, coeff_shift);
            SIMD_FUNC(cdef_filter_block_4x4_8)
                (dst8 + 4 * dstride, dstride, in + 4 * CDEF_BSTRIDE, pri_strength,
                    sec_strength, dir, pri_damping, sec_damping, max, coeff_shift);
        }
        else if (bsize == BLOCK_8X4) {
            SIMD_FUNC(cdef_filter_block_4x4_8)
                (dst8, dstride, in, pri_strength, sec_strength, dir, pri_damping,
                    sec_damping, max, coeff_shift);
            SIMD_FUNC(cdef_filter_block_4x4_8)
                (dst8 + 4, dstride, in + 4, pri_strength, sec_strength, dir, pri_damping,
                    sec_damping, max, coeff_shift);
        }
        else {
            SIMD_FUNC(cdef_filter_block_4x4_8)
                (dst8, dstride, in, pri_strength, sec_strength, dir, pri_damping,
                    sec_damping, max, coeff_shift);
        }
    }
    else {
        if (bsize == BLOCK_8X8) {
            SIMD_FUNC(cdef_filter_block_8x8_16)
                (dst16, dstride, in, pri_strength, sec_strength, dir, pri_damping,
                    sec_damping, max, coeff_shift);
        }
        else if (bsize == BLOCK_4X8) {
            SIMD_FUNC(cdef_filter_block_4x4_16)
                (dst16, dstride, in, pri_strength, sec_strength, dir, pri_damping,
                    sec_damping, max, coeff_shift);
            SIMD_FUNC(cdef_filter_block_4x4_16)
                (dst16 + 4 * dstride, dstride, in + 4 * CDEF_BSTRIDE, pri_strength,
                    sec_strength, dir, pri_damping, sec_damping, max, coeff_shift);
        }
        else if (bsize == BLOCK_8X4) {
            SIMD_FUNC(cdef_filter_block_4x4_16)
                (dst16, dstride, in, pri_strength, sec_strength, dir, pri_damping,
                    sec_damping, max, coeff_shift);
            SIMD_FUNC(cdef_filter_block_4x4_16)
                (dst16 + 4, dstride, in + 4, pri_strength, sec_strength, dir, pri_damping,
                    sec_damping, max, coeff_shift);
        }
        else {
            assert(bsize == BLOCK_4X4);
            SIMD_FUNC(cdef_filter_block_4x4_16)
                (dst16, dstride, in, pri_strength, sec_strength, dir, pri_damping,
                    sec_damping, max, coeff_shift);
        }
    }
}

void SIMD_FUNC(copy_rect8_8bit_to_16bit)(uint16_t *dst, int32_t dstride,
    const uint8_t *src, int32_t sstride, int32_t v,
    int32_t h) {
    int32_t i, j;
    for (i = 0; i < v; i++) {
        for (j = 0; j < (h & ~0x7); j += 8) {
            v64 row = v64_load_unaligned(&src[i * sstride + j]);
            v128_store_unaligned(&dst[i * dstride + j], v128_unpack_u8_s16(row));
        }
        for (; j < h; j++)
            dst[i * dstride + j] = src[i * sstride + j];
    }
}

void SIMD_FUNC(copy_rect8_16bit_to_16bit)(uint16_t *dst, int32_t dstride,
    const uint16_t *src, int32_t sstride,
    int32_t v, int32_t h) {
    int32_t i, j;
    for (i = 0; i < v; i++) {
        for (j = 0; j < (h & ~0x7); j += 8) {
            v128 row = v128_load_unaligned(&src[i * sstride + j]);
            v128_store_unaligned(&dst[i * dstride + j], row);
        }
        for (; j < h; j++)
            dst[i * dstride + j] = src[i * sstride + j];
    }
}
