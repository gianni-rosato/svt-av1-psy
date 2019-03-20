/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbDefinitions.h"
#include <immintrin.h>
#include <math.h>
#include "synonyms.h"

#define REDUCED_PRI_STRENGTHS 8
#define REDUCED_TOTAL_STRENGTHS (REDUCED_PRI_STRENGTHS * CDEF_SEC_STRENGTHS)
#define TOTAL_STRENGTHS (CDEF_PRI_STRENGTHS * CDEF_SEC_STRENGTHS)


/* Search for the best luma+chroma strength to add as an option, knowing we
already selected nb_strengths options. */
#if FAST_CDEF
uint64_t search_one_dual_avx2(int *lev0, int *lev1, int nb_strengths,
	uint64_t(**mse)[TOTAL_STRENGTHS], int sb_count,
	int fast, int start_gi, int end_gi) {
#else
uint64_t search_one_dual_avx2(int *lev0, int *lev1, int nb_strengths,
                              uint64_t(**mse)[TOTAL_STRENGTHS], int sb_count,
                              int fast) {
#endif

  DECLARE_ALIGNED(32, uint64_t, tot_mse[TOTAL_STRENGTHS][TOTAL_STRENGTHS]);
  int i, j;
  uint64_t best_tot_mse = (uint64_t)1 << 62;
  int best_id0 = 0;
  int best_id1 = 0;
#if FAST_CDEF
  (void)fast;
  const int total_strengths = end_gi;
#else
  const int total_strengths = fast ? REDUCED_TOTAL_STRENGTHS : TOTAL_STRENGTHS;
#endif
  __m256i best_mse_;
  __m256i curr;
  __m256i v_tot;
  __m256i v_mse;
  __m256i mask;
  __m256i tmp;

  memset(tot_mse, 0, sizeof(tot_mse));

  for (i = 0; i < sb_count; i++) {
    int gi;
    uint64_t best_mse = (uint64_t)1 << 62;
    /* Find best mse among already selected options. */
    for (gi = 0; gi < nb_strengths; gi++) {
      uint64_t curr = mse[0][i][lev0[gi]];
      curr += mse[1][i][lev1[gi]];
      if (curr < best_mse) {
        best_mse = curr;
      }
    }
    best_mse_ = _mm256_set1_epi64x(best_mse);
    /* Find best mse when adding each possible new option. */
    assert(~total_strengths % 4);
#if FAST_CDEF
	for (int j = start_gi; j < total_strengths; ++j) { // process by 4x4
		tmp = _mm256_set1_epi64x(mse[0][i][j]);
		for (int k = 0; k < total_strengths; k += 4) {
#else
    for (int j = 0; j < total_strengths; ++j) { // process by 4x4
      tmp = _mm256_set1_epi64x(mse[0][i][j]);
      for (int k = 0; k < total_strengths; k += 4) {
#endif
        v_mse = _mm256_loadu_si256((const __m256i*)&mse[1][i][k]);
        v_tot = _mm256_loadu_si256((const __m256i*)&tot_mse[j][k]);
        curr = _mm256_add_epi64(tmp, v_mse);
        mask = _mm256_cmpgt_epi64(best_mse_, curr);
        v_tot = _mm256_add_epi64(v_tot, _mm256_or_si256(
          _mm256_andnot_si256(mask, best_mse_),
          _mm256_and_si256(mask, curr)));
        _mm256_storeu_si256((__m256i*)&tot_mse[j][k], v_tot);
      }
    }
  }
#if FAST_CDEF
  for (j = start_gi; j < total_strengths; j++) {
	  int k;
	  for (k = start_gi; k < total_strengths; k++) {
#else
  for (j = 0; j < total_strengths; j++) {
    int k;
    for (k = 0; k < total_strengths; k++) {
#endif
      if (tot_mse[j][k] < best_tot_mse) {
        best_tot_mse = tot_mse[j][k];
        best_id0 = j;
        best_id1 = k;
      }
    }
  }
  lev0[nb_strengths] = best_id0;
  lev1[nb_strengths] = best_id1;

  return best_tot_mse;
}

uint64_t mse_4x4_16bit_avx2(uint16_t *dst, int dstride, uint16_t *src,
                            int sstride) {
  __m256i src_t, dst_t, low_1, low_2, high_1, high_2;
  uint64_t ret;

  src_t = _mm256_set_epi64x(*(uint64_t*)(src + 0 * sstride),
      *(uint64_t*)(src + 2 * sstride),
      *(uint64_t*)(src + 1 * sstride),
      *(uint64_t*)(src + 3 * sstride));
  dst_t = _mm256_set_epi64x(*(uint64_t*)(dst + 0 * dstride),
      *(uint64_t*)(dst + 2 * dstride),
      *(uint64_t*)(dst + 1 * dstride),
      *(uint64_t*)(dst + 3 * dstride));

  dst_t = _mm256_sub_epi16(dst_t, src_t);

  src_t = _mm256_mullo_epi16(dst_t, dst_t);

  src_t = _mm256_hadd_epi16(src_t, src_t);//indexes 0...8
  src_t = _mm256_hadd_epi16(src_t, src_t);//indexes 0,1,4,5
  src_t = _mm256_hadd_epi16(src_t, src_t);//indexes 0,4
  ret = _mm256_extract_epi16(src_t, 0) + _mm256_extract_epi16(src_t, 8);

  return ret;
}

static INLINE uint32_t sum32(__m256i v) {
  return _mm256_extract_epi32(v, 0) +
      _mm256_extract_epi32(v, 1) +
      _mm256_extract_epi32(v, 2) +
      _mm256_extract_epi32(v, 3) +
      _mm256_extract_epi32(v, 4) +
      _mm256_extract_epi32(v, 5) +
      _mm256_extract_epi32(v, 6) +
      _mm256_extract_epi32(v, 7);
}

uint64_t dist_8x8_16bit_avx2(uint16_t *dst, int dstride, uint16_t *src,
                             int sstride, int coeff_shift) {
  __m256i m_sum_s = _mm256_cvtepu16_epi32(_mm_load_si128((const __m128i*)src));
  __m256i m_sum_d = _mm256_cvtepu16_epi32(_mm_load_si128((const __m128i*)dst));
  __m256i m_sum_s2 = _mm256_mullo_epi32(m_sum_s, m_sum_s);
  __m256i m_sum_sd = _mm256_mullo_epi32(m_sum_s, m_sum_d);
  __m256i m_sum_d2 = _mm256_mullo_epi32(m_sum_d, m_sum_d);
  for (unsigned r = 1; r < 8; r++) {
    __m256i s = _mm256_cvtepu16_epi32(
        _mm_load_si128((const __m128i*)(src + r * sstride)));
    __m256i d = _mm256_cvtepu16_epi32(
        _mm_load_si128((const __m128i*)(dst + r * dstride)));
    m_sum_s = _mm256_add_epi32(m_sum_s, s);
    m_sum_d = _mm256_add_epi32(m_sum_d, d);
    m_sum_s2 = _mm256_add_epi32(m_sum_s2, _mm256_mullo_epi32(s, s));
    m_sum_sd = _mm256_add_epi32(m_sum_sd, _mm256_mullo_epi32(s, d));
    m_sum_d2 = _mm256_add_epi32(m_sum_d2, _mm256_mullo_epi32(d, d));
  }
  /* Compute the variance -- the calculation cannot go negative. */
  uint64_t sum_s = sum32(m_sum_s);
  uint64_t sum_d = sum32(m_sum_d);
  uint64_t sum_s2 = sum32(m_sum_s2);
  uint64_t sum_d2 = sum32(m_sum_d2);
  uint64_t sum_sd = sum32(m_sum_sd);
  /* Compute the variance -- the calculation cannot go negative. */
  uint64_t svar = sum_s2 - ((sum_s * sum_s + 32) >> 6);
  uint64_t dvar = sum_d2 - ((sum_d * sum_d + 32) >> 6);
  uint64_t ret = (uint64_t)floor(
    .5 + (sum_d2 + sum_s2 - 2 * sum_sd) * .5 *
    (svar + dvar + (400 << 2 * coeff_shift)) /
    (sqrt((20000 << 4 * coeff_shift) + svar * (double)dvar)));

  return ret;
}
