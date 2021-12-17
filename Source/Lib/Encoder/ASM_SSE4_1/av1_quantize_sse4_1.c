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
#include "aom_dsp_rtcd.h"
#include <emmintrin.h>
#include <smmintrin.h>
#include <xmmintrin.h>

static INLINE void read_coeff(const TranLow *coeff, intptr_t offset,
                              __m128i *c0, __m128i *c1) {
    const TranLow *addr = coeff + offset;
    if (sizeof(TranLow) == 4) {
    const __m128i x0 = _mm_loadu_si128((const __m128i *)addr);
    const __m128i x1 = _mm_loadu_si128((const __m128i *)addr + 1);
    const __m128i x2 = _mm_loadu_si128((const __m128i *)addr + 2);
    const __m128i x3 = _mm_loadu_si128((const __m128i *)addr + 3);
    *c0 = _mm_packs_epi32(x0, x1);
    *c1 = _mm_packs_epi32(x2, x3);
  } else {
    *c0 = _mm_loadu_si128((const __m128i *)addr);
    *c1 = _mm_loadu_si128((const __m128i *)addr + 1);
  }
}

static INLINE void write_qcoeff(const __m128i *qc0, const __m128i *qc1, TranLow *qcoeff,
                                intptr_t offset) {
    TranLow *addr = qcoeff + offset;
    if (sizeof(TranLow) == 4) {
    const __m128i zero = _mm_setzero_si128();
    __m128i sign_bits = _mm_cmplt_epi16(*qc0, zero);
    __m128i y0 = _mm_unpacklo_epi16(*qc0, sign_bits);
    __m128i y1 = _mm_unpackhi_epi16(*qc0, sign_bits);
    _mm_storeu_si128((__m128i *)addr, y0);
    _mm_storeu_si128((__m128i *)addr + 1, y1);

    sign_bits = _mm_cmplt_epi16(*qc1, zero);
    y0 = _mm_unpacklo_epi16(*qc1, sign_bits);
    y1 = _mm_unpackhi_epi16(*qc1, sign_bits);
    _mm_storeu_si128((__m128i *)addr + 2, y0);
    _mm_storeu_si128((__m128i *)addr + 3, y1);
  } else {
    _mm_storeu_si128((__m128i *)addr, *qc0);
    _mm_storeu_si128((__m128i *)addr + 1, *qc1);
  }
}

static INLINE void write_zero(TranLow *qcoeff, intptr_t offset) {
  const __m128i zero = _mm_setzero_si128();
  TranLow *     addr = qcoeff + offset;
  if (sizeof(TranLow) == 4) {
    _mm_storeu_si128((__m128i *)addr, zero);
    _mm_storeu_si128((__m128i *)addr + 1, zero);
    _mm_storeu_si128((__m128i *)addr + 2, zero);
    _mm_storeu_si128((__m128i *)addr + 3, zero);
  } else {
    _mm_storeu_si128((__m128i *)addr, zero);
    _mm_storeu_si128((__m128i *)addr + 1, zero);
  }
}

static INLINE void quantize(const int16_t *iscan_ptr, const TranLow *coeff_ptr, intptr_t n_coeffs,
                            TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr,
                            const __m128i *round0, const __m128i *round1,
                            const __m128i *quant0, const __m128i *quant1,
                            const __m128i *dequant0, const __m128i *dequant1,
                            const __m128i *thr0, const __m128i *thr1,
                            __m128i *eob) {
  __m128i coeff0, coeff1;
  // Do DC and first 15 AC
  read_coeff(coeff_ptr, n_coeffs, &coeff0, &coeff1);

  // Poor man's sign extract
  const __m128i coeff0_sign = _mm_srai_epi16(coeff0, 15);
  const __m128i coeff1_sign = _mm_srai_epi16(coeff1, 15);
  __m128i qcoeff0 = _mm_xor_si128(coeff0, coeff0_sign);
  __m128i qcoeff1 = _mm_xor_si128(coeff1, coeff1_sign);
  qcoeff0 = _mm_sub_epi16(qcoeff0, coeff0_sign);
  qcoeff1 = _mm_sub_epi16(qcoeff1, coeff1_sign);
  const __m128i mask0 = _mm_or_si128(_mm_cmpgt_epi16(qcoeff0, *thr0),
                                     _mm_cmpeq_epi16(qcoeff0, *thr0));
  const __m128i mask1 = _mm_or_si128(_mm_cmpgt_epi16(qcoeff1, *thr1),
                                     _mm_cmpeq_epi16(qcoeff1, *thr1));
  const int nzflag = _mm_movemask_epi8(mask0) | _mm_movemask_epi8(mask1);

  if (nzflag) {
    qcoeff0 = _mm_adds_epi16(qcoeff0, *round0);
    qcoeff1 = _mm_adds_epi16(qcoeff1, *round1);
    const __m128i qtmp0 = _mm_mulhi_epi16(qcoeff0, *quant0);
    const __m128i qtmp1 = _mm_mulhi_epi16(qcoeff1, *quant1);

    // Reinsert signs
    qcoeff0 = _mm_xor_si128(qtmp0, coeff0_sign);
    qcoeff1 = _mm_xor_si128(qtmp1, coeff1_sign);
    qcoeff0 = _mm_sub_epi16(qcoeff0, coeff0_sign);
    qcoeff1 = _mm_sub_epi16(qcoeff1, coeff1_sign);

    write_qcoeff(&qcoeff0, &qcoeff1, qcoeff_ptr, n_coeffs);

    coeff0 = _mm_mullo_epi16(qcoeff0, *dequant0);
    coeff1 = _mm_mullo_epi16(qcoeff1, *dequant1);

    write_qcoeff(&coeff0, &coeff1, dqcoeff_ptr, n_coeffs);

    const __m128i zero = _mm_setzero_si128();
    // Scan for eob
    const __m128i zero_coeff0 = _mm_cmpeq_epi16(coeff0, zero);
    const __m128i zero_coeff1 = _mm_cmpeq_epi16(coeff1, zero);
    const __m128i nzero_coeff0 = _mm_cmpeq_epi16(zero_coeff0, zero);
    const __m128i nzero_coeff1 = _mm_cmpeq_epi16(zero_coeff1, zero);
    const __m128i iscan0 =
        _mm_loadu_si128((const __m128i *)(iscan_ptr + n_coeffs));
    const __m128i iscan1 =
        _mm_loadu_si128((const __m128i *)(iscan_ptr + n_coeffs) + 1);
    // Add one to convert from indices to counts
    const __m128i iscan0_nz = _mm_sub_epi16(iscan0, nzero_coeff0);
    const __m128i iscan1_nz = _mm_sub_epi16(iscan1, nzero_coeff1);
    const __m128i eob0 = _mm_and_si128(iscan0_nz, nzero_coeff0);
    const __m128i eob1 = _mm_and_si128(iscan1_nz, nzero_coeff1);
    const __m128i eob2 = _mm_max_epi16(eob0, eob1);
    *eob = _mm_max_epi16(*eob, eob2);
  } else {
    write_zero(qcoeff_ptr, n_coeffs);
    write_zero(dqcoeff_ptr, n_coeffs);
  }
}

void svt_av1_quantize_fp_sse4_1(const TranLow *coeff_ptr, intptr_t n_coeffs, const int16_t *zbin_ptr,
                              const int16_t *round_ptr, const int16_t *quant_ptr,
                              const int16_t *quant_shift_ptr, TranLow *qcoeff_ptr,
                              TranLow *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr,
                              const int16_t *scan_ptr, const int16_t *iscan_ptr) {
  (void)scan_ptr;
  (void)zbin_ptr;
  (void)quant_shift_ptr;

  coeff_ptr += n_coeffs;
  iscan_ptr += n_coeffs;
  qcoeff_ptr += n_coeffs;
  dqcoeff_ptr += n_coeffs;
  n_coeffs = -n_coeffs;

  const __m128i round0 = _mm_loadu_si128((const __m128i *)round_ptr);
  const __m128i round1 = _mm_unpackhi_epi64(round0, round0);
  const __m128i quant0 = _mm_loadu_si128((const __m128i *)quant_ptr);
  const __m128i quant1 = _mm_unpackhi_epi64(quant0, quant0);
  const __m128i dequant0 = _mm_loadu_si128((const __m128i *)dequant_ptr);
  const __m128i dequant1 = _mm_unpackhi_epi64(dequant0, dequant0);
  const __m128i thr0 = _mm_srai_epi16(dequant0, 1);
  const __m128i thr1 = _mm_srai_epi16(dequant1, 1);
  __m128i eob = _mm_setzero_si128();

  quantize(iscan_ptr, coeff_ptr, n_coeffs, qcoeff_ptr, dqcoeff_ptr, &round0,
           &round1, &quant0, &quant1, &dequant0, &dequant1, &thr0, &thr1, &eob);

  n_coeffs += 8 * 2;

  // AC only loop
  while (n_coeffs < 0) {
    quantize(iscan_ptr, coeff_ptr, n_coeffs, qcoeff_ptr, dqcoeff_ptr, &round1,
             &round1, &quant1, &quant1, &dequant1, &dequant1, &thr1, &thr1,
             &eob);
    n_coeffs += 8 * 2;
  }

  // Accumulate EOB
  {
    __m128i eob_shuffled;
    eob_shuffled = _mm_shuffle_epi32(eob, 0xe);
    eob = _mm_max_epi16(eob, eob_shuffled);
    eob_shuffled = _mm_shufflelo_epi16(eob, 0xe);
    eob = _mm_max_epi16(eob, eob_shuffled);
    eob_shuffled = _mm_shufflelo_epi16(eob, 0x1);
    eob = _mm_max_epi16(eob, eob_shuffled);
    *eob_ptr = _mm_extract_epi16(eob, 1);
  }
}


static INLINE void init_qp_add_shift(const int16_t *zbin_ptr, const int16_t *round_ptr,
                                     const int16_t *quant_ptr, const int16_t *dequant_ptr,
                                     const int16_t *quant_shift_ptr, __m128i *qp,
                                     const int add_shift) {
    __m128i       zbin        = _mm_loadl_epi64((const __m128i *)zbin_ptr);
    __m128i       round       = _mm_loadl_epi64((const __m128i *)round_ptr);
    const __m128i quant       = _mm_loadl_epi64((const __m128i *)quant_ptr);
    const __m128i dequant     = _mm_loadl_epi64((const __m128i *)dequant_ptr);
    const __m128i quant_shift = _mm_loadl_epi64((const __m128i *)quant_shift_ptr);
    if (add_shift) {
        const __m128i add = _mm_set1_epi16((int16_t)add_shift);
        zbin              = _mm_add_epi16(zbin, add);
        round             = _mm_add_epi16(round, add);
        zbin              = _mm_srli_epi16(zbin, add_shift);
        round             = _mm_srli_epi16(round, add_shift);
    }
    qp[0] = _mm_cvtepi16_epi32(zbin);
    qp[1] = _mm_cvtepi16_epi32(round);
    qp[2] = _mm_cvtepi16_epi32(quant);
    qp[3] = _mm_cvtepi16_epi32(dequant);
    qp[4] = _mm_cvtepi16_epi32(quant_shift);
}

static INLINE void mm_mul_shift_epi32(const __m128i *x, const __m128i *y, __m128i *p,
                                         int shift) {
    __m128i prod_lo       = _mm_mul_epi32(*x, *y);
    prod_lo               = _mm_srli_epi64(prod_lo, shift);
    __m128i       prod_hi = _mm_srli_epi64(*x, 32);
    const __m128i mult_hi = _mm_srli_epi64(*y, 32);
    prod_hi               = _mm_mul_epi32(prod_hi, mult_hi);

    prod_hi = _mm_srli_epi64(prod_hi, shift);

    prod_hi = _mm_slli_epi64(prod_hi, 32);
    *p      = _mm_blend_epi16(prod_lo, prod_hi, 0xCC); //interleave prod_lo, prod_hi
}

static INLINE void clamp_epi32(__m128i *x, __m128i min, __m128i max) {
    *x = _mm_min_epi32(*x, max);
    *x = _mm_max_epi32(*x, min);
}

static INLINE void quantize_128(const __m128i *qp, __m128i c, const int16_t *iscan_ptr, TranLow *qcoeff,
                            TranLow *dqcoeff, __m128i *eob, __m128i min, __m128i max,
                            int shift_dq) {
    const __m128i zero   = _mm_setzero_si128();
    const __m128i abs    = _mm_abs_epi32(c);
    const __m128i flag1  = _mm_cmpgt_epi32(qp[0], abs);
    const int32_t nzflag = _mm_movemask_epi8(flag1);

    if (LIKELY(~nzflag)) {
        __m128i q = _mm_add_epi32(abs, qp[1]);
        clamp_epi32(&q, min, max);
        __m128i tmp;
        mm_mul_shift_epi32(&q, &qp[2], &tmp, 16);
        q = _mm_add_epi32(tmp, q);

        mm_mul_shift_epi32(&q, &qp[4], &q, 16 - shift_dq);
        __m128i dq = _mm_mullo_epi32(q, qp[3]);
        dq         = _mm_srli_epi32(dq, shift_dq);

        q  = _mm_sign_epi32(q, c);
        dq = _mm_sign_epi32(dq, c);
        q  = _mm_andnot_si128(flag1, q);
        dq = _mm_andnot_si128(flag1, dq);

        _mm_storeu_si128((__m128i *)qcoeff, q);
        _mm_storeu_si128((__m128i *)dqcoeff, dq);

        const __m128i isc   = _mm_loadl_epi64((const __m128i *)iscan_ptr);
        const __m128i iscan = _mm_cvtepi16_epi32(isc);

        const __m128i zc      = _mm_cmpeq_epi32(dq, zero);
        const __m128i nz      = _mm_cmpeq_epi32(zc, zero);
        __m128i       cur_eob = _mm_sub_epi32(iscan, nz);
        cur_eob               = _mm_and_si128(cur_eob, nz);
        *eob                  = _mm_max_epi32(cur_eob, *eob);
    } else {
        _mm_storeu_si128((__m128i *)qcoeff, zero);
        _mm_storeu_si128((__m128i *)dqcoeff, zero);
    }
}

void svt_aom_quantize_b_sse4_1(const TranLow *coeff_ptr, intptr_t n_coeffs, const int16_t *zbin_ptr,
                             const int16_t *round_ptr, const int16_t *quant_ptr,
                             const int16_t *quant_shift_ptr, TranLow *qcoeff_ptr,
                             TranLow *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr,
                             const int16_t *scan, const int16_t *iscan, const QmVal *qm_ptr,
                             const QmVal *iqm_ptr, const int32_t log_scale) {
    (void)qm_ptr;
    (void)iqm_ptr;
    (void)scan;
    const uint32_t step = 4;

    __m128i qp[5], coeff;
    init_qp_add_shift(zbin_ptr, round_ptr, quant_ptr, dequant_ptr, quant_shift_ptr, qp, log_scale);
    coeff = _mm_loadu_si128((const __m128i *)coeff_ptr);

    __m128i eob = _mm_setzero_si128();
    __m128i min = _mm_set1_epi32(INT16_MIN);
    __m128i max = _mm_set1_epi32(INT16_MAX);
    quantize_128(qp, coeff, iscan, qcoeff_ptr, dqcoeff_ptr, &eob, min, max, log_scale);
    //Equivalent to update_qp(qp) from svt_aom_quantize_b_avx2 is calling init_qp_add_shift() with pointers + 4
    init_qp_add_shift(zbin_ptr + 4, round_ptr + 4, quant_ptr + 4, dequant_ptr + 4, quant_shift_ptr + 4, qp, log_scale);

    while (n_coeffs > step) {
        coeff_ptr += step;
        qcoeff_ptr += step;
        dqcoeff_ptr += step;
        iscan += step;
        n_coeffs -= step;

        coeff = _mm_loadu_si128((const __m128i *)coeff_ptr);
        quantize_128(qp, coeff, iscan, qcoeff_ptr, dqcoeff_ptr, &eob, min, max, log_scale);
    }
    {
        __m128i eob_s;
        eob_s                   = _mm_shuffle_epi32(eob, 0xe);
        eob                     = _mm_max_epi16(eob, eob_s);
        eob_s                   = _mm_shufflelo_epi16(eob, 0xe);
        eob                     = _mm_max_epi16(eob, eob_s);
        eob_s                   = _mm_shufflelo_epi16(eob, 1);
        const __m128i final_eob = _mm_max_epi16(eob, eob_s);
        *eob_ptr                = _mm_extract_epi16(final_eob, 0);
    }
}
