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

#include "EbDefinitions.h"
#include <immintrin.h>

#include "aom_dsp_rtcd.h"

// Note: TranHigh is the datatype used for intermediate transform stages.
typedef int64_t TranHigh;

static INLINE void init_one_qp(const __m128i *p, __m256i *qp) {
    const __m128i sign = _mm_srai_epi16(*p, 15);
    const __m128i dc   = _mm_unpacklo_epi16(*p, sign);
    const __m128i ac   = _mm_unpackhi_epi16(*p, sign);
    *qp                = _mm256_insertf128_si256(_mm256_castsi128_si256(dc), ac, 1);
}

static INLINE void update_qp(__m256i *qp) {
    int32_t i;
    for (i = 0; i < 5; ++i) qp[i] = _mm256_permute2x128_si256(qp[i], qp[i], 0x11);
}

static INLINE void init_qp(const int16_t *zbin_ptr, const int16_t *round_ptr,
                           const int16_t *quant_ptr, const int16_t *dequant_ptr,
                           const int16_t *quant_shift_ptr, __m256i *qp) {
    const __m128i zbin        = _mm_loadu_si128((const __m128i *)zbin_ptr);
    const __m128i round       = _mm_loadu_si128((const __m128i *)round_ptr);
    const __m128i quant       = _mm_loadu_si128((const __m128i *)quant_ptr);
    const __m128i dequant     = _mm_loadu_si128((const __m128i *)dequant_ptr);
    const __m128i quant_shift = _mm_loadu_si128((const __m128i *)quant_shift_ptr);
    init_one_qp(&zbin, &qp[0]);
    init_one_qp(&round, &qp[1]);
    init_one_qp(&quant, &qp[2]);
    init_one_qp(&dequant, &qp[3]);
    init_one_qp(&quant_shift, &qp[4]);
}

// Note:
// *x is vector multiplied by *y which is 16 int32_t parallel multiplication
// and right shift 16.  The output, 16 int32_t is save in *p.
static INLINE void mm256_mul_shift_epi32(const __m256i *x, const __m256i *y, __m256i *p) {
    __m256i       prod_lo = _mm256_mul_epi32(*x, *y);
    __m256i       prod_hi = _mm256_srli_epi64(*x, 32);
    const __m256i mult_hi = _mm256_srli_epi64(*y, 32);
    prod_hi               = _mm256_mul_epi32(prod_hi, mult_hi);

    prod_lo            = _mm256_srli_epi64(prod_lo, 16);
    const __m256i mask = _mm256_set_epi32(0, -1, 0, -1, 0, -1, 0, -1);
    prod_lo            = _mm256_and_si256(prod_lo, mask);
    prod_hi            = _mm256_srli_epi64(prod_hi, 16);

    prod_hi = _mm256_slli_epi64(prod_hi, 32);
    *p      = _mm256_or_si256(prod_lo, prod_hi);
}

static INLINE void clamp_epi32(__m256i *x, __m256i min, __m256i max) {
    *x = _mm256_min_epi32(*x, max);
    *x = _mm256_max_epi32(*x, min);
}

static INLINE void quantize(const __m256i *qp, __m256i *c, const int16_t *iscan_ptr,
                            TranLow *qcoeff, TranLow *dqcoeff, __m256i *eob, __m256i min,
                            __m256i max) {
    const __m256i abs    = _mm256_abs_epi32(*c);
    const __m256i flag1  = _mm256_cmpgt_epi32(abs, qp[0]);
    __m256i       flag2  = _mm256_cmpeq_epi32(abs, qp[0]);
    flag2                = _mm256_or_si256(flag1, flag2);
    const int32_t nzflag = _mm256_movemask_epi8(flag2);

    if (LIKELY(nzflag)) {
        __m256i q = _mm256_add_epi32(abs, qp[1]);
        clamp_epi32(&q, min, max);
        __m256i tmp;
        mm256_mul_shift_epi32(&q, &qp[2], &tmp);
        q = _mm256_add_epi32(tmp, q);

        mm256_mul_shift_epi32(&q, &qp[4], &q);
        __m256i dq = _mm256_mullo_epi32(q, qp[3]);

        q  = _mm256_sign_epi32(q, *c);
        dq = _mm256_sign_epi32(dq, *c);
        q  = _mm256_and_si256(q, flag2);
        dq = _mm256_and_si256(dq, flag2);

        _mm256_storeu_si256((__m256i *)qcoeff, q);
        _mm256_storeu_si256((__m256i *)dqcoeff, dq);

        const __m128i isc   = _mm_loadu_si128((const __m128i *)iscan_ptr);
        const __m128i zr    = _mm_setzero_si128();
        const __m128i lo    = _mm_unpacklo_epi16(isc, zr);
        const __m128i hi    = _mm_unpackhi_epi16(isc, zr);
        const __m256i iscan = _mm256_insertf128_si256(_mm256_castsi128_si256(lo), hi, 1);

        const __m256i zero    = _mm256_setzero_si256();
        const __m256i zc      = _mm256_cmpeq_epi32(dq, zero);
        const __m256i nz      = _mm256_cmpeq_epi32(zc, zero);
        __m256i       cur_eob = _mm256_sub_epi32(iscan, nz);
        cur_eob               = _mm256_and_si256(cur_eob, nz);
        *eob                  = _mm256_max_epi32(cur_eob, *eob);
    } else {
        const __m256i zero = _mm256_setzero_si256();
        _mm256_storeu_si256((__m256i *)qcoeff, zero);
        _mm256_storeu_si256((__m256i *)dqcoeff, zero);
    }
}

void eb_aom_quantize_b_avx2(const TranLow *coeff_ptr, intptr_t n_coeffs, int32_t skip_block,
                            const int16_t *zbin_ptr, const int16_t *round_ptr,
                            const int16_t *quant_ptr, const int16_t *quant_shift_ptr,
                            TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr, const int16_t *dequant_ptr,
                            uint16_t *eob_ptr, const int16_t *scan, const int16_t *iscan) {
    (void)scan;
    const uint32_t step = 8;

    if (LIKELY(!skip_block)) {
        __m256i qp[5], coeff;
        init_qp(zbin_ptr, round_ptr, quant_ptr, dequant_ptr, quant_shift_ptr, qp);
        coeff = _mm256_loadu_si256((const __m256i *)coeff_ptr);

        __m256i eob = _mm256_setzero_si256();
        __m256i min = _mm256_set1_epi32(INT16_MIN);
        __m256i max = _mm256_set1_epi32(INT16_MAX);
        quantize(qp, &coeff, iscan, qcoeff_ptr, dqcoeff_ptr, &eob, min, max);

        coeff_ptr += step;
        qcoeff_ptr += step;
        dqcoeff_ptr += step;
        iscan += step;
        n_coeffs -= step;

        update_qp(qp);

        while (n_coeffs > 0) {
            coeff = _mm256_loadu_si256((const __m256i *)coeff_ptr);
            quantize(qp, &coeff, iscan, qcoeff_ptr, dqcoeff_ptr, &eob, min, max);

            coeff_ptr += step;
            qcoeff_ptr += step;
            dqcoeff_ptr += step;
            iscan += step;
            n_coeffs -= step;
        }
        {
            __m256i eob_s;
            eob_s = _mm256_shuffle_epi32(eob, 0xe);
            eob   = _mm256_max_epi16(eob, eob_s);
            eob_s = _mm256_shufflelo_epi16(eob, 0xe);
            eob   = _mm256_max_epi16(eob, eob_s);
            eob_s = _mm256_shufflelo_epi16(eob, 1);
            eob   = _mm256_max_epi16(eob, eob_s);
            const __m128i final_eob =
                _mm_max_epi16(_mm256_castsi256_si128(eob), _mm256_extractf128_si256(eob, 1));
            *eob_ptr = _mm_extract_epi16(final_eob, 0);
        }
    } else {
        do {
            const __m256i zero = _mm256_setzero_si256();
            _mm256_storeu_si256((__m256i *)qcoeff_ptr, zero);
            _mm256_storeu_si256((__m256i *)dqcoeff_ptr, zero);
            qcoeff_ptr += step;
            dqcoeff_ptr += step;
            n_coeffs -= step;
        } while (n_coeffs > 0);
        *eob_ptr = 0;
    }
}

static INLINE void quantize_highbd(const __m256i *qp, __m256i *c, const int16_t *iscan_ptr,
                                   TranLow *qcoeff, TranLow *dqcoeff, __m256i *eob) {
    const __m256i abs    = _mm256_abs_epi32(*c);
    const __m256i flag1  = _mm256_cmpgt_epi32(abs, qp[0]);
    __m256i       flag2  = _mm256_cmpeq_epi32(abs, qp[0]);
    flag2                = _mm256_or_si256(flag1, flag2);
    const int32_t nzflag = _mm256_movemask_epi8(flag2);

    if (LIKELY(nzflag)) {
        __m256i q = _mm256_add_epi32(abs, qp[1]);
        __m256i tmp;
        mm256_mul_shift_epi32(&q, &qp[2], &tmp);
        q = _mm256_add_epi32(tmp, q);

        mm256_mul_shift_epi32(&q, &qp[4], &q);
        __m256i dq = _mm256_mullo_epi32(q, qp[3]);

        q  = _mm256_sign_epi32(q, *c);
        dq = _mm256_sign_epi32(dq, *c);
        q  = _mm256_and_si256(q, flag2);
        dq = _mm256_and_si256(dq, flag2);

        _mm256_storeu_si256((__m256i *)qcoeff, q);
        _mm256_storeu_si256((__m256i *)dqcoeff, dq);

        const __m128i isc   = _mm_loadu_si128((const __m128i *)iscan_ptr);
        const __m128i zr    = _mm_setzero_si128();
        const __m128i lo    = _mm_unpacklo_epi16(isc, zr);
        const __m128i hi    = _mm_unpackhi_epi16(isc, zr);
        const __m256i iscan = _mm256_insertf128_si256(_mm256_castsi128_si256(lo), hi, 1);

        const __m256i zero    = _mm256_setzero_si256();
        const __m256i zc      = _mm256_cmpeq_epi32(dq, zero);
        const __m256i nz      = _mm256_cmpeq_epi32(zc, zero);
        __m256i       cur_eob = _mm256_sub_epi32(iscan, nz);
        cur_eob               = _mm256_and_si256(cur_eob, nz);
        *eob                  = _mm256_max_epi32(cur_eob, *eob);
    } else {
        const __m256i zero = _mm256_setzero_si256();
        _mm256_storeu_si256((__m256i *)qcoeff, zero);
        _mm256_storeu_si256((__m256i *)dqcoeff, zero);
    }
}

void eb_aom_highbd_quantize_b_avx2(const TranLow *coeff_ptr, intptr_t n_coeffs, int32_t skip_block,
                                   const int16_t *zbin_ptr, const int16_t *round_ptr,
                                   const int16_t *quant_ptr, const int16_t *quant_shift_ptr,
                                   TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr,
                                   const int16_t *dequant_ptr, uint16_t *eob_ptr,
                                   const int16_t *scan, const int16_t *iscan) {
    (void)scan;
    const uint32_t step = 8;

    if (LIKELY(!skip_block)) {
        __m256i qp[5], coeff;
        init_qp(zbin_ptr, round_ptr, quant_ptr, dequant_ptr, quant_shift_ptr, qp);
        coeff = _mm256_loadu_si256((const __m256i *)coeff_ptr);

        __m256i eob = _mm256_setzero_si256();
        quantize_highbd(qp, &coeff, iscan, qcoeff_ptr, dqcoeff_ptr, &eob);

        coeff_ptr += step;
        qcoeff_ptr += step;
        dqcoeff_ptr += step;
        iscan += step;
        n_coeffs -= step;

        update_qp(qp);

        while (n_coeffs > 0) {
            coeff = _mm256_loadu_si256((const __m256i *)coeff_ptr);
            quantize_highbd(qp, &coeff, iscan, qcoeff_ptr, dqcoeff_ptr, &eob);

            coeff_ptr += step;
            qcoeff_ptr += step;
            dqcoeff_ptr += step;
            iscan += step;
            n_coeffs -= step;
        }
        {
            __m256i eob_s;
            eob_s = _mm256_shuffle_epi32(eob, 0xe);
            eob   = _mm256_max_epi16(eob, eob_s);
            eob_s = _mm256_shufflelo_epi16(eob, 0xe);
            eob   = _mm256_max_epi16(eob, eob_s);
            eob_s = _mm256_shufflelo_epi16(eob, 1);
            eob   = _mm256_max_epi16(eob, eob_s);
            const __m128i final_eob =
                _mm_max_epi16(_mm256_castsi256_si128(eob), _mm256_extractf128_si256(eob, 1));
            *eob_ptr = _mm_extract_epi16(final_eob, 0);
        }
    } else {
        do {
            const __m256i zero = _mm256_setzero_si256();
            _mm256_storeu_si256((__m256i *)qcoeff_ptr, zero);
            _mm256_storeu_si256((__m256i *)dqcoeff_ptr, zero);
            qcoeff_ptr += step;
            dqcoeff_ptr += step;
            n_coeffs -= step;
        } while (n_coeffs > 0);
        *eob_ptr = 0;
    }
}

static INLINE void init_qp_64x64(const int16_t *zbin_ptr, const int16_t *round_ptr,
                                 const int16_t *quant_ptr, const int16_t *dequant_ptr,
                                 const int16_t *quant_shift_ptr, __m256i *qp) {
    __m128i       zbin        = _mm_loadu_si128((const __m128i *)zbin_ptr);
    __m128i       round       = _mm_loadu_si128((const __m128i *)round_ptr);
    const __m128i quant       = _mm_loadu_si128((const __m128i *)quant_ptr);
    const __m128i dequant     = _mm_loadu_si128((const __m128i *)dequant_ptr);
    const __m128i quant_shift = _mm_loadu_si128((const __m128i *)quant_shift_ptr);
    const __m128i add         = _mm_set_epi16(2, 2, 2, 2, 2, 2, 2, 2);
    zbin                      = _mm_add_epi16(zbin, add);
    round                     = _mm_add_epi16(round, add);
    zbin                      = _mm_srli_epi16(zbin, 2);
    round                     = _mm_srli_epi16(round, 2);
    init_one_qp(&zbin, &qp[0]);
    init_one_qp(&round, &qp[1]);
    init_one_qp(&quant, &qp[2]);
    init_one_qp(&dequant, &qp[3]);
    init_one_qp(&quant_shift, &qp[4]);
}

static INLINE void mm256_mul_shift_epi32_64x64(const __m256i *x, const __m256i *y, __m256i *p) {
    __m256i       prod_lo = _mm256_mul_epi32(*x, *y);
    __m256i       prod_hi = _mm256_srli_epi64(*x, 32);
    const __m256i mult_hi = _mm256_srli_epi64(*y, 32);
    prod_hi               = _mm256_mul_epi32(prod_hi, mult_hi);

    prod_lo            = _mm256_srli_epi64(prod_lo, 14);
    const __m256i mask = _mm256_set_epi32(0, -1, 0, -1, 0, -1, 0, -1);
    prod_lo            = _mm256_and_si256(prod_lo, mask);
    prod_hi            = _mm256_srli_epi64(prod_hi, 14);

    prod_hi = _mm256_slli_epi64(prod_hi, 32);
    *p      = _mm256_or_si256(prod_lo, prod_hi);
}

static INLINE void quantize_64x64(const __m256i *qp, __m256i *c, const int16_t *iscan_ptr,
                                  TranLow *qcoeff, TranLow *dqcoeff, __m256i *eob, __m256i min,
                                  __m256i max) {
    const __m256i abs    = _mm256_abs_epi32(*c);
    const __m256i flag1  = _mm256_cmpgt_epi32(abs, qp[0]);
    __m256i       flag2  = _mm256_cmpeq_epi32(abs, qp[0]);
    flag2                = _mm256_or_si256(flag1, flag2);
    const int32_t nzflag = _mm256_movemask_epi8(flag2);

    if (LIKELY(nzflag)) {
        __m256i q = _mm256_add_epi32(abs, qp[1]);
        clamp_epi32(&q, min, max);
        __m256i tmp;
        mm256_mul_shift_epi32(&q, &qp[2], &tmp);
        q = _mm256_add_epi32(tmp, q);

        mm256_mul_shift_epi32_64x64(&q, &qp[4], &q);
        __m256i dq = _mm256_mullo_epi32(q, qp[3]);
        dq         = _mm256_srli_epi32(dq, 2);

        q  = _mm256_sign_epi32(q, *c);
        dq = _mm256_sign_epi32(dq, *c);
        q  = _mm256_and_si256(q, flag2);
        dq = _mm256_and_si256(dq, flag2);

        _mm256_storeu_si256((__m256i *)qcoeff, q);
        _mm256_storeu_si256((__m256i *)dqcoeff, dq);

        const __m128i isc   = _mm_loadu_si128((const __m128i *)iscan_ptr);
        const __m128i zr    = _mm_setzero_si128();
        const __m128i lo    = _mm_unpacklo_epi16(isc, zr);
        const __m128i hi    = _mm_unpackhi_epi16(isc, zr);
        const __m256i iscan = _mm256_insertf128_si256(_mm256_castsi128_si256(lo), hi, 1);

        const __m256i zero    = _mm256_setzero_si256();
        const __m256i zc      = _mm256_cmpeq_epi32(dq, zero);
        const __m256i nz      = _mm256_cmpeq_epi32(zc, zero);
        __m256i       cur_eob = _mm256_sub_epi32(iscan, nz);
        cur_eob               = _mm256_and_si256(cur_eob, nz);
        *eob                  = _mm256_max_epi32(cur_eob, *eob);
    } else {
        const __m256i zero = _mm256_setzero_si256();
        _mm256_storeu_si256((__m256i *)qcoeff, zero);
        _mm256_storeu_si256((__m256i *)dqcoeff, zero);
    }
}

void eb_aom_quantize_b_64x64_avx2(const TranLow *coeff_ptr, intptr_t n_coeffs, int32_t skip_block,
                                  const int16_t *zbin_ptr, const int16_t *round_ptr,
                                  const int16_t *quant_ptr, const int16_t *quant_shift_ptr,
                                  TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr,
                                  const int16_t *dequant_ptr, uint16_t *eob_ptr,
                                  const int16_t *scan, const int16_t *iscan) {
    (void)scan;
    const uint32_t step = 8;

    if (LIKELY(!skip_block)) {
        __m256i qp[5], coeff;
        init_qp_64x64(zbin_ptr, round_ptr, quant_ptr, dequant_ptr, quant_shift_ptr, qp);
        coeff = _mm256_loadu_si256((const __m256i *)coeff_ptr);

        __m256i eob = _mm256_setzero_si256();
        __m256i min = _mm256_set1_epi32(INT16_MIN);
        __m256i max = _mm256_set1_epi32(INT16_MAX);
        quantize_64x64(qp, &coeff, iscan, qcoeff_ptr, dqcoeff_ptr, &eob, min, max);

        coeff_ptr += step;
        qcoeff_ptr += step;
        dqcoeff_ptr += step;
        iscan += step;
        n_coeffs -= step;

        update_qp(qp);

        while (n_coeffs > 0) {
            coeff = _mm256_loadu_si256((const __m256i *)coeff_ptr);
            quantize_64x64(qp, &coeff, iscan, qcoeff_ptr, dqcoeff_ptr, &eob, min, max);

            coeff_ptr += step;
            qcoeff_ptr += step;
            dqcoeff_ptr += step;
            iscan += step;
            n_coeffs -= step;
        }
        {
            __m256i eob_s;
            eob_s = _mm256_shuffle_epi32(eob, 0xe);
            eob   = _mm256_max_epi16(eob, eob_s);
            eob_s = _mm256_shufflelo_epi16(eob, 0xe);
            eob   = _mm256_max_epi16(eob, eob_s);
            eob_s = _mm256_shufflelo_epi16(eob, 1);
            eob   = _mm256_max_epi16(eob, eob_s);
            const __m128i final_eob =
                _mm_max_epi16(_mm256_castsi256_si128(eob), _mm256_extractf128_si256(eob, 1));
            *eob_ptr = _mm_extract_epi16(final_eob, 0);
        }
    } else {
        do {
            const __m256i zero = _mm256_setzero_si256();
            _mm256_storeu_si256((__m256i *)qcoeff_ptr, zero);
            _mm256_storeu_si256((__m256i *)dqcoeff_ptr, zero);
            qcoeff_ptr += step;
            dqcoeff_ptr += step;
            n_coeffs -= step;
        } while (n_coeffs > 0);
        *eob_ptr = 0;
    }
}
static INLINE void quantize_64x64_highbd(const __m256i *qp, __m256i *c, const int16_t *iscan_ptr,
                                         TranLow *qcoeff, TranLow *dqcoeff, __m256i *eob) {
    const __m256i abs    = _mm256_abs_epi32(*c);
    const __m256i flag1  = _mm256_cmpgt_epi32(abs, qp[0]);
    __m256i       flag2  = _mm256_cmpeq_epi32(abs, qp[0]);
    flag2                = _mm256_or_si256(flag1, flag2);
    const int32_t nzflag = _mm256_movemask_epi8(flag2);

    if (LIKELY(nzflag)) {
        __m256i q = _mm256_add_epi32(abs, qp[1]);
        __m256i tmp;
        mm256_mul_shift_epi32(&q, &qp[2], &tmp);
        q = _mm256_add_epi32(tmp, q);

        mm256_mul_shift_epi32_64x64(&q, &qp[4], &q);
        __m256i dq = _mm256_mullo_epi32(q, qp[3]);
        dq         = _mm256_srli_epi32(dq, 2);

        q  = _mm256_sign_epi32(q, *c);
        dq = _mm256_sign_epi32(dq, *c);
        q  = _mm256_and_si256(q, flag2);
        dq = _mm256_and_si256(dq, flag2);

        _mm256_storeu_si256((__m256i *)qcoeff, q);
        _mm256_storeu_si256((__m256i *)dqcoeff, dq);

        const __m128i isc   = _mm_loadu_si128((const __m128i *)iscan_ptr);
        const __m128i zr    = _mm_setzero_si128();
        const __m128i lo    = _mm_unpacklo_epi16(isc, zr);
        const __m128i hi    = _mm_unpackhi_epi16(isc, zr);
        const __m256i iscan = _mm256_insertf128_si256(_mm256_castsi128_si256(lo), hi, 1);

        const __m256i zero    = _mm256_setzero_si256();
        const __m256i zc      = _mm256_cmpeq_epi32(dq, zero);
        const __m256i nz      = _mm256_cmpeq_epi32(zc, zero);
        __m256i       cur_eob = _mm256_sub_epi32(iscan, nz);
        cur_eob               = _mm256_and_si256(cur_eob, nz);
        *eob                  = _mm256_max_epi32(cur_eob, *eob);
    } else {
        const __m256i zero = _mm256_setzero_si256();
        _mm256_storeu_si256((__m256i *)qcoeff, zero);
        _mm256_storeu_si256((__m256i *)dqcoeff, zero);
    }
}

void eb_aom_highbd_quantize_b_64x64_avx2(const TranLow *coeff_ptr, intptr_t n_coeffs,
                                         int32_t skip_block, const int16_t *zbin_ptr,
                                         const int16_t *round_ptr, const int16_t *quant_ptr,
                                         const int16_t *quant_shift_ptr, TranLow *qcoeff_ptr,
                                         TranLow *dqcoeff_ptr, const int16_t *dequant_ptr,
                                         uint16_t *eob_ptr, const int16_t *scan,
                                         const int16_t *iscan) {
    (void)scan;
    const uint32_t step = 8;

    if (LIKELY(!skip_block)) {
        __m256i qp[5], coeff;
        init_qp_64x64(zbin_ptr, round_ptr, quant_ptr, dequant_ptr, quant_shift_ptr, qp);
        coeff = _mm256_loadu_si256((const __m256i *)coeff_ptr);

        __m256i eob = _mm256_setzero_si256();
        quantize_64x64_highbd(qp, &coeff, iscan, qcoeff_ptr, dqcoeff_ptr, &eob);

        coeff_ptr += step;
        qcoeff_ptr += step;
        dqcoeff_ptr += step;
        iscan += step;
        n_coeffs -= step;

        update_qp(qp);

        while (n_coeffs > 0) {
            coeff = _mm256_loadu_si256((const __m256i *)coeff_ptr);
            quantize_64x64_highbd(qp, &coeff, iscan, qcoeff_ptr, dqcoeff_ptr, &eob);

            coeff_ptr += step;
            qcoeff_ptr += step;
            dqcoeff_ptr += step;
            iscan += step;
            n_coeffs -= step;
        }
        {
            __m256i eob_s;
            eob_s = _mm256_shuffle_epi32(eob, 0xe);
            eob   = _mm256_max_epi16(eob, eob_s);
            eob_s = _mm256_shufflelo_epi16(eob, 0xe);
            eob   = _mm256_max_epi16(eob, eob_s);
            eob_s = _mm256_shufflelo_epi16(eob, 1);
            eob   = _mm256_max_epi16(eob, eob_s);
            const __m128i final_eob =
                _mm_max_epi16(_mm256_castsi256_si128(eob), _mm256_extractf128_si256(eob, 1));
            *eob_ptr = _mm_extract_epi16(final_eob, 0);
        }
    } else {
        do {
            const __m256i zero = _mm256_setzero_si256();
            _mm256_storeu_si256((__m256i *)qcoeff_ptr, zero);
            _mm256_storeu_si256((__m256i *)dqcoeff_ptr, zero);
            qcoeff_ptr += step;
            dqcoeff_ptr += step;
            n_coeffs -= step;
        } while (n_coeffs > 0);
        *eob_ptr = 0;
    }
}
static INLINE void init_qp_32x32(const int16_t *zbin_ptr, const int16_t *round_ptr,
                                 const int16_t *quant_ptr, const int16_t *dequant_ptr,
                                 const int16_t *quant_shift_ptr, __m256i *qp) {
    __m128i       zbin        = _mm_loadu_si128((const __m128i *)zbin_ptr);
    __m128i       round       = _mm_loadu_si128((const __m128i *)round_ptr);
    const __m128i quant       = _mm_loadu_si128((const __m128i *)quant_ptr);
    const __m128i dequant     = _mm_loadu_si128((const __m128i *)dequant_ptr);
    const __m128i quant_shift = _mm_loadu_si128((const __m128i *)quant_shift_ptr);
    const __m128i add         = _mm_set_epi16(1, 1, 1, 1, 1, 1, 1, 1);
    zbin                      = _mm_add_epi16(zbin, add);
    round                     = _mm_add_epi16(round, add);
    zbin                      = _mm_srli_epi16(zbin, 1);
    round                     = _mm_srli_epi16(round, 1);
    init_one_qp(&zbin, &qp[0]);
    init_one_qp(&round, &qp[1]);
    init_one_qp(&quant, &qp[2]);
    init_one_qp(&dequant, &qp[3]);
    init_one_qp(&quant_shift, &qp[4]);
}

static INLINE void mm256_mul_shift_epi32_32x32(const __m256i *x, const __m256i *y, __m256i *p) {
    __m256i       prod_lo = _mm256_mul_epi32(*x, *y);
    __m256i       prod_hi = _mm256_srli_epi64(*x, 32);
    const __m256i mult_hi = _mm256_srli_epi64(*y, 32);
    prod_hi               = _mm256_mul_epi32(prod_hi, mult_hi);

    prod_lo            = _mm256_srli_epi64(prod_lo, 15);
    const __m256i mask = _mm256_set_epi32(0, -1, 0, -1, 0, -1, 0, -1);
    prod_lo            = _mm256_and_si256(prod_lo, mask);
    prod_hi            = _mm256_srli_epi64(prod_hi, 15);

    prod_hi = _mm256_slli_epi64(prod_hi, 32);
    *p      = _mm256_or_si256(prod_lo, prod_hi);
}

static INLINE void quantize_32x32(const __m256i *qp, __m256i *c, const int16_t *iscan_ptr,
                                  TranLow *qcoeff, TranLow *dqcoeff, __m256i *eob, __m256i min,
                                  __m256i max) {
    const __m256i abs    = _mm256_abs_epi32(*c);
    const __m256i flag1  = _mm256_cmpgt_epi32(abs, qp[0]);
    __m256i       flag2  = _mm256_cmpeq_epi32(abs, qp[0]);
    flag2                = _mm256_or_si256(flag1, flag2);
    const int32_t nzflag = _mm256_movemask_epi8(flag2);

    if (LIKELY(nzflag)) {
        __m256i q = _mm256_add_epi32(abs, qp[1]);
        clamp_epi32(&q, min, max);
        __m256i tmp;
        mm256_mul_shift_epi32(&q, &qp[2], &tmp);
        q = _mm256_add_epi32(tmp, q);

        mm256_mul_shift_epi32_32x32(&q, &qp[4], &q);
        __m256i dq = _mm256_mullo_epi32(q, qp[3]);
        dq         = _mm256_srli_epi32(dq, 1);

        q  = _mm256_sign_epi32(q, *c);
        dq = _mm256_sign_epi32(dq, *c);
        q  = _mm256_and_si256(q, flag2);
        dq = _mm256_and_si256(dq, flag2);

        _mm256_storeu_si256((__m256i *)qcoeff, q);
        _mm256_storeu_si256((__m256i *)dqcoeff, dq);

        const __m128i isc   = _mm_loadu_si128((const __m128i *)iscan_ptr);
        const __m128i zr    = _mm_setzero_si128();
        const __m128i lo    = _mm_unpacklo_epi16(isc, zr);
        const __m128i hi    = _mm_unpackhi_epi16(isc, zr);
        const __m256i iscan = _mm256_insertf128_si256(_mm256_castsi128_si256(lo), hi, 1);

        const __m256i zero    = _mm256_setzero_si256();
        const __m256i zc      = _mm256_cmpeq_epi32(dq, zero);
        const __m256i nz      = _mm256_cmpeq_epi32(zc, zero);
        __m256i       cur_eob = _mm256_sub_epi32(iscan, nz);
        cur_eob               = _mm256_and_si256(cur_eob, nz);
        *eob                  = _mm256_max_epi32(cur_eob, *eob);
    } else {
        const __m256i zero = _mm256_setzero_si256();
        _mm256_storeu_si256((__m256i *)qcoeff, zero);
        _mm256_storeu_si256((__m256i *)dqcoeff, zero);
    }
}

void eb_aom_quantize_b_32x32_avx2(const TranLow *coeff_ptr, intptr_t n_coeffs, int skip_block,
                                  const int16_t *zbin_ptr, const int16_t *round_ptr,
                                  const int16_t *quant_ptr, const int16_t *quant_shift_ptr,
                                  TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr,
                                  const int16_t *dequant_ptr, uint16_t *eob_ptr,
                                  const int16_t *scan, const int16_t *iscan) {
    (void)scan;
    const unsigned int step = 8;

    if (LIKELY(!skip_block)) {
        __m256i qp[5], coeff;
        init_qp_32x32(zbin_ptr, round_ptr, quant_ptr, dequant_ptr, quant_shift_ptr, qp);
        coeff = _mm256_loadu_si256((const __m256i *)coeff_ptr);

        __m256i eob = _mm256_setzero_si256();
        __m256i min = _mm256_set1_epi32(INT16_MIN);
        __m256i max = _mm256_set1_epi32(INT16_MAX);
        quantize_32x32(qp, &coeff, iscan, qcoeff_ptr, dqcoeff_ptr, &eob, min, max);

        coeff_ptr += step;
        qcoeff_ptr += step;
        dqcoeff_ptr += step;
        iscan += step;
        n_coeffs -= step;

        update_qp(qp);

        while (n_coeffs > 0) {
            coeff = _mm256_loadu_si256((const __m256i *)coeff_ptr);
            quantize_32x32(qp, &coeff, iscan, qcoeff_ptr, dqcoeff_ptr, &eob, min, max);

            coeff_ptr += step;
            qcoeff_ptr += step;
            dqcoeff_ptr += step;
            iscan += step;
            n_coeffs -= step;
        }
        {
            __m256i eob_s;
            eob_s = _mm256_shuffle_epi32(eob, 0xe);
            eob   = _mm256_max_epi16(eob, eob_s);
            eob_s = _mm256_shufflelo_epi16(eob, 0xe);
            eob   = _mm256_max_epi16(eob, eob_s);
            eob_s = _mm256_shufflelo_epi16(eob, 1);
            eob   = _mm256_max_epi16(eob, eob_s);
            const __m128i final_eob =
                _mm_max_epi16(_mm256_castsi256_si128(eob), _mm256_extractf128_si256(eob, 1));
            *eob_ptr = _mm_extract_epi16(final_eob, 0);
        }
    } else {
        do {
            const __m256i zero = _mm256_setzero_si256();
            _mm256_storeu_si256((__m256i *)qcoeff_ptr, zero);
            _mm256_storeu_si256((__m256i *)dqcoeff_ptr, zero);
            qcoeff_ptr += step;
            dqcoeff_ptr += step;
            n_coeffs -= step;
        } while (n_coeffs > 0);
        *eob_ptr = 0;
    }
}

static INLINE void quantize_32x32_highbd(const __m256i *qp, __m256i *c, const int16_t *iscan_ptr,
                                         TranLow *qcoeff, TranLow *dqcoeff, __m256i *eob) {
    const __m256i abs    = _mm256_abs_epi32(*c);
    const __m256i flag1  = _mm256_cmpgt_epi32(abs, qp[0]);
    __m256i       flag2  = _mm256_cmpeq_epi32(abs, qp[0]);
    flag2                = _mm256_or_si256(flag1, flag2);
    const int32_t nzflag = _mm256_movemask_epi8(flag2);

    if (LIKELY(nzflag)) {
        __m256i q = _mm256_add_epi32(abs, qp[1]);
        __m256i tmp;
        mm256_mul_shift_epi32(&q, &qp[2], &tmp);
        q = _mm256_add_epi32(tmp, q);

        mm256_mul_shift_epi32_32x32(&q, &qp[4], &q);
        __m256i dq = _mm256_mullo_epi32(q, qp[3]);
        dq         = _mm256_srli_epi32(dq, 1);

        q  = _mm256_sign_epi32(q, *c);
        dq = _mm256_sign_epi32(dq, *c);
        q  = _mm256_and_si256(q, flag2);
        dq = _mm256_and_si256(dq, flag2);

        _mm256_storeu_si256((__m256i *)qcoeff, q);
        _mm256_storeu_si256((__m256i *)dqcoeff, dq);

        const __m128i isc   = _mm_loadu_si128((const __m128i *)iscan_ptr);
        const __m128i zr    = _mm_setzero_si128();
        const __m128i lo    = _mm_unpacklo_epi16(isc, zr);
        const __m128i hi    = _mm_unpackhi_epi16(isc, zr);
        const __m256i iscan = _mm256_insertf128_si256(_mm256_castsi128_si256(lo), hi, 1);

        const __m256i zero    = _mm256_setzero_si256();
        const __m256i zc      = _mm256_cmpeq_epi32(dq, zero);
        const __m256i nz      = _mm256_cmpeq_epi32(zc, zero);
        __m256i       cur_eob = _mm256_sub_epi32(iscan, nz);
        cur_eob               = _mm256_and_si256(cur_eob, nz);
        *eob                  = _mm256_max_epi32(cur_eob, *eob);
    } else {
        const __m256i zero = _mm256_setzero_si256();
        _mm256_storeu_si256((__m256i *)qcoeff, zero);
        _mm256_storeu_si256((__m256i *)dqcoeff, zero);
    }
}

void eb_aom_highbd_quantize_b_32x32_avx2(const TranLow *coeff_ptr, intptr_t n_coeffs,
                                         int skip_block, const int16_t *zbin_ptr,
                                         const int16_t *round_ptr, const int16_t *quant_ptr,
                                         const int16_t *quant_shift_ptr, TranLow *qcoeff_ptr,
                                         TranLow *dqcoeff_ptr, const int16_t *dequant_ptr,
                                         uint16_t *eob_ptr, const int16_t *scan,
                                         const int16_t *iscan) {
    (void)scan;
    const unsigned int step = 8;

    if (LIKELY(!skip_block)) {
        __m256i qp[5], coeff;
        init_qp_32x32(zbin_ptr, round_ptr, quant_ptr, dequant_ptr, quant_shift_ptr, qp);
        coeff = _mm256_loadu_si256((const __m256i *)coeff_ptr);

        __m256i eob = _mm256_setzero_si256();
        quantize_32x32_highbd(qp, &coeff, iscan, qcoeff_ptr, dqcoeff_ptr, &eob);

        coeff_ptr += step;
        qcoeff_ptr += step;
        dqcoeff_ptr += step;
        iscan += step;
        n_coeffs -= step;

        update_qp(qp);

        while (n_coeffs > 0) {
            coeff = _mm256_loadu_si256((const __m256i *)coeff_ptr);
            quantize_32x32_highbd(qp, &coeff, iscan, qcoeff_ptr, dqcoeff_ptr, &eob);

            coeff_ptr += step;
            qcoeff_ptr += step;
            dqcoeff_ptr += step;
            iscan += step;
            n_coeffs -= step;
        }
        {
            __m256i eob_s;
            eob_s = _mm256_shuffle_epi32(eob, 0xe);
            eob   = _mm256_max_epi16(eob, eob_s);
            eob_s = _mm256_shufflelo_epi16(eob, 0xe);
            eob   = _mm256_max_epi16(eob, eob_s);
            eob_s = _mm256_shufflelo_epi16(eob, 1);
            eob   = _mm256_max_epi16(eob, eob_s);
            const __m128i final_eob =
                _mm_max_epi16(_mm256_castsi256_si128(eob), _mm256_extractf128_si256(eob, 1));
            *eob_ptr = _mm_extract_epi16(final_eob, 0);
        }
    } else {
        do {
            const __m256i zero = _mm256_setzero_si256();
            _mm256_storeu_si256((__m256i *)qcoeff_ptr, zero);
            _mm256_storeu_si256((__m256i *)dqcoeff_ptr, zero);
            qcoeff_ptr += step;
            dqcoeff_ptr += step;
            n_coeffs -= step;
        } while (n_coeffs > 0);
        *eob_ptr = 0;
    }
}

static INLINE void init_one_qp_fp(const __m128i *p, __m256i *qp) {
    const __m128i zero = _mm_setzero_si128();
    const __m128i dc   = _mm_unpacklo_epi16(*p, zero);
    const __m128i ac   = _mm_unpackhi_epi16(*p, zero);
    *qp                = _mm256_insertf128_si256(_mm256_castsi128_si256(dc), ac, 1);
}

static INLINE void init_qp_fp(const int16_t *round_ptr, const int16_t *quant_ptr,
                              const int16_t *dequant_ptr, int log_scale, __m256i *qp) {
    __m128i round = _mm_loadu_si128((const __m128i *)round_ptr);
    if (log_scale) {
        const __m128i round_scale = _mm_set1_epi16(1 << (15 - log_scale));
        round                     = _mm_mulhrs_epi16(round, round_scale);
    }
    const __m128i quant   = _mm_loadu_si128((const __m128i *)quant_ptr);
    const __m128i dequant = _mm_loadu_si128((const __m128i *)dequant_ptr);

    init_one_qp_fp(&round, &qp[0]);
    init_one_qp_fp(&quant, &qp[1]);
    init_one_qp_fp(&dequant, &qp[2]);
}

static INLINE void quantize_highbd_fp(const __m256i *qp, __m256i *c, const int16_t *iscan_ptr,
                                      int log_scale, TranLow *qcoeff, TranLow *dqcoeff,
                                      __m256i *eob) {
    const __m256i abs_coeff = _mm256_abs_epi32(*c);
    __m256i       q         = _mm256_add_epi32(abs_coeff, qp[0]);

    __m256i       q_lo  = _mm256_mul_epi32(q, qp[1]);
    __m256i       q_hi  = _mm256_srli_epi64(q, 32);
    const __m256i qp_hi = _mm256_srli_epi64(qp[1], 32);
    q_hi                = _mm256_mul_epi32(q_hi, qp_hi);
    q_lo                = _mm256_srli_epi64(q_lo, 16 - log_scale);
    q_hi                = _mm256_srli_epi64(q_hi, 16 - log_scale);
    q_hi                = _mm256_slli_epi64(q_hi, 32);
    q                   = _mm256_or_si256(q_lo, q_hi);
    const __m256i abs_s = _mm256_slli_epi32(abs_coeff, 1 + log_scale);
    const __m256i mask  = _mm256_cmpgt_epi32(qp[2], abs_s);
    q                   = _mm256_andnot_si256(mask, q);

    __m256i dq = _mm256_mullo_epi32(q, qp[2]);
    dq         = _mm256_srai_epi32(dq, log_scale);
    q          = _mm256_sign_epi32(q, *c);
    dq         = _mm256_sign_epi32(dq, *c);

    _mm256_storeu_si256((__m256i *)qcoeff, q);
    _mm256_storeu_si256((__m256i *)dqcoeff, dq);

    const __m128i isc   = _mm_loadu_si128((const __m128i *)iscan_ptr);
    const __m128i zr    = _mm_setzero_si128();
    const __m128i lo    = _mm_unpacklo_epi16(isc, zr);
    const __m128i hi    = _mm_unpackhi_epi16(isc, zr);
    const __m256i iscan = _mm256_insertf128_si256(_mm256_castsi128_si256(lo), hi, 1);

    const __m256i zero    = _mm256_setzero_si256();
    const __m256i zc      = _mm256_cmpeq_epi32(dq, zero);
    const __m256i nz      = _mm256_cmpeq_epi32(zc, zero);
    __m256i       cur_eob = _mm256_sub_epi32(iscan, nz);
    cur_eob               = _mm256_and_si256(cur_eob, nz);
    *eob                  = _mm256_max_epi32(cur_eob, *eob);
}

static INLINE void update_qp_fp(__m256i *qp) {
    qp[0] = _mm256_permute2x128_si256(qp[0], qp[0], 0x11);
    qp[1] = _mm256_permute2x128_si256(qp[1], qp[1], 0x11);
    qp[2] = _mm256_permute2x128_si256(qp[2], qp[2], 0x11);
}

void eb_av1_highbd_quantize_fp_avx2(const TranLow *coeff_ptr, intptr_t n_coeffs,
                                    const int16_t *zbin_ptr, const int16_t *round_ptr,
                                    const int16_t *quant_ptr, const int16_t *quant_shift_ptr,
                                    TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr,
                                    const int16_t *dequant_ptr, uint16_t *eob_ptr,
                                    const int16_t *scan, const int16_t *iscan, int16_t log_scale) {
    (void)scan;
    (void)zbin_ptr;
    (void)quant_shift_ptr;
    const unsigned int step = 8;
    __m256i            qp[3], coeff;

    init_qp_fp(round_ptr, quant_ptr, dequant_ptr, log_scale, qp);
    coeff = _mm256_loadu_si256((const __m256i *)coeff_ptr);

    __m256i eob = _mm256_setzero_si256();
    quantize_highbd_fp(qp, &coeff, iscan, log_scale, qcoeff_ptr, dqcoeff_ptr, &eob);

    coeff_ptr += step;
    qcoeff_ptr += step;
    dqcoeff_ptr += step;
    iscan += step;
    n_coeffs -= step;

    update_qp_fp(qp);
    while (n_coeffs > 0) {
        coeff = _mm256_loadu_si256((const __m256i *)coeff_ptr);
        quantize_highbd_fp(qp, &coeff, iscan, log_scale, qcoeff_ptr, dqcoeff_ptr, &eob);

        coeff_ptr += step;
        qcoeff_ptr += step;
        dqcoeff_ptr += step;
        iscan += step;
        n_coeffs -= step;
    }

    {
        __m256i eob_s;
        eob_s = _mm256_shuffle_epi32(eob, 0xe);
        eob   = _mm256_max_epi16(eob, eob_s);
        eob_s = _mm256_shufflelo_epi16(eob, 0xe);
        eob   = _mm256_max_epi16(eob, eob_s);
        eob_s = _mm256_shufflelo_epi16(eob, 1);
        eob   = _mm256_max_epi16(eob, eob_s);
        const __m128i final_eob =
            _mm_max_epi16(_mm256_castsi256_si128(eob), _mm256_extractf128_si256(eob, 1));
        *eob_ptr = _mm_extract_epi16(final_eob, 0);
    }
}
