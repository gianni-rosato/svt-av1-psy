/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbDefinitions.h"
#include <immintrin.h>
#include <math.h>

#define REDUCED_PRI_STRENGTHS 8
#define REDUCED_TOTAL_STRENGTHS (REDUCED_PRI_STRENGTHS * CDEF_SEC_STRENGTHS)
#define TOTAL_STRENGTHS (CDEF_PRI_STRENGTHS * CDEF_SEC_STRENGTHS)

#ifndef _mm256_set_m128i
#define _mm256_set_m128i(/* __m128i */ hi, /* __m128i */ lo) \
    _mm256_insertf128_si256(_mm256_castsi128_si256(lo), (hi), 0x1)
#endif

#ifndef _mm256_setr_m128i
#define _mm256_setr_m128i(/* __m128i */ lo, /* __m128i */ hi) _mm256_set_m128i((hi), (lo))
#endif

/* Search for the best luma+chroma strength to add as an option, knowing we
already selected nb_strengths options. */
uint64_t search_one_dual_avx2(int *lev0, int *lev1, int nb_strengths,
                              uint64_t (**mse)[TOTAL_STRENGTHS], int sb_count, int fast,
                              int start_gi, int end_gi) {
    DECLARE_ALIGNED(32, uint64_t, tot_mse[TOTAL_STRENGTHS][TOTAL_STRENGTHS]);
    int      i, j;
    uint64_t best_tot_mse = (uint64_t)1 << 62;
    int      best_id0     = 0;
    int      best_id1     = 0;
    (void)fast;
    const int total_strengths = end_gi;
    __m256i   best_mse_;
    __m256i   curr;
    __m256i   v_tot;
    __m256i   v_mse;
    __m256i   mask;
    __m256i   tmp;

    memset(tot_mse, 0, sizeof(tot_mse));

    for (i = 0; i < sb_count; i++) {
        int      gi;
        uint64_t best_mse = (uint64_t)1 << 62;
        /* Find best mse among already selected options. */
        for (gi = 0; gi < nb_strengths; gi++) {
            uint64_t curr = mse[0][i][lev0[gi]];
            curr += mse[1][i][lev1[gi]];
            if (curr < best_mse) best_mse = curr;
        }
        best_mse_ = _mm256_set1_epi64x(best_mse);
        /* Find best mse when adding each possible new option. */
        //assert(~total_strengths % 4);
        for (int j = start_gi; j < total_strengths; ++j) { // process by 4x4
            tmp = _mm256_set1_epi64x(mse[0][i][j]);
            for (int k = 0; k < total_strengths; k += 4) {
                v_mse = _mm256_loadu_si256((const __m256i *)&mse[1][i][k]);
                v_tot = _mm256_loadu_si256((const __m256i *)&tot_mse[j][k]);
                curr  = _mm256_add_epi64(tmp, v_mse);
                mask  = _mm256_cmpgt_epi64(best_mse_, curr);
                v_tot = _mm256_add_epi64(v_tot,
                                         _mm256_or_si256(_mm256_andnot_si256(mask, best_mse_),
                                                         _mm256_and_si256(mask, curr)));
                _mm256_storeu_si256((__m256i *)&tot_mse[j][k], v_tot);
            }
        }
    }
    for (j = start_gi; j < total_strengths; j++) {
        int k;
        for (k = start_gi; k < total_strengths; k++) {
            if (tot_mse[j][k] < best_tot_mse) {
                best_tot_mse = tot_mse[j][k];
                best_id0     = j;
                best_id1     = k;
            }
        }
    }
    lev0[nb_strengths] = best_id0;
    lev1[nb_strengths] = best_id1;

    return best_tot_mse;
}

static INLINE void mse_4x4_16bit_avx2(const uint16_t **src, const uint16_t *dst,
                                      const int32_t dstride, __m256i *sum) {
    const __m256i s    = _mm256_loadu_si256((const __m256i *)*src);
    const __m256i d    = _mm256_setr_epi64x(*(uint64_t *)(dst + 0 * dstride),
                                         *(uint64_t *)(dst + 1 * dstride),
                                         *(uint64_t *)(dst + 2 * dstride),
                                         *(uint64_t *)(dst + 3 * dstride));
    const __m256i diff = _mm256_sub_epi16(d, s);
    const __m256i mse  = _mm256_madd_epi16(diff, diff);
    *sum               = _mm256_add_epi32(*sum, mse);
    *src += 16;
}

static INLINE void mse_4x4_8bit_avx2(const uint8_t **src, const uint8_t *dst, const int32_t dstride,
                                     __m256i *sum) {
    const __m128i s = _mm_loadu_si128((const __m128i *)*src);
    const __m128i d = _mm_setr_epi32(*(uint32_t *)(dst + 0 * dstride),
                                     *(uint32_t *)(dst + 1 * dstride),
                                     *(uint32_t *)(dst + 2 * dstride),
                                     *(uint32_t *)(dst + 3 * dstride));

    const __m256i s_16 = _mm256_cvtepu8_epi16(s);
    const __m256i d_16 = _mm256_cvtepu8_epi16(d);

    const __m256i diff = _mm256_sub_epi16(d_16, s_16);
    const __m256i mse  = _mm256_madd_epi16(diff, diff);
    *sum               = _mm256_add_epi32(*sum, mse);
    *src += 16;
}

static INLINE void mse_8x2_16bit_avx2(const uint16_t **src, const uint16_t *dst,
                                      const int32_t dstride, __m256i *sum) {
    const __m256i s    = _mm256_loadu_si256((const __m256i *)*src);
    const __m128i d0   = _mm_loadu_si128((const __m128i *)(dst + 0 * dstride));
    const __m128i d1   = _mm_loadu_si128((const __m128i *)(dst + 1 * dstride));
    const __m256i d    = _mm256_setr_m128i(d0, d1);
    const __m256i diff = _mm256_sub_epi16(d, s);
    const __m256i mse  = _mm256_madd_epi16(diff, diff);
    *sum               = _mm256_add_epi32(*sum, mse);
    *src += 16;
}

static INLINE void mse_8x2_8bit_avx2(const uint8_t **src, const uint8_t *dst, const int32_t dstride,
                                     __m256i *sum) {
    const __m128i s = _mm_loadu_si128((const __m128i *)*src);
    const __m128i d =
        _mm_set_epi64x(*(uint64_t *)(dst + 1 * dstride), *(uint64_t *)(dst + 0 * dstride));

    const __m256i s_16 = _mm256_cvtepu8_epi16(s);
    const __m256i d_16 = _mm256_cvtepu8_epi16(d);

    const __m256i diff = _mm256_sub_epi16(d_16, s_16);
    const __m256i mse  = _mm256_madd_epi16(diff, diff);
    *sum               = _mm256_add_epi32(*sum, mse);
    *src += 16;
}

static INLINE void mse_8x4_16bit_avx2(const uint16_t **src, const uint16_t *dst,
                                      const int32_t dstride, __m256i *sum) {
    mse_8x2_16bit_avx2(src, dst + 0 * dstride, dstride, sum);
    mse_8x2_16bit_avx2(src, dst + 2 * dstride, dstride, sum);
}

static INLINE void mse_8x4_8bit_avx2(const uint8_t **src, const uint8_t *dst, const int32_t dstride,
                                     __m256i *sum) {
    mse_8x2_8bit_avx2(src, dst + 0 * dstride, dstride, sum);
    mse_8x2_8bit_avx2(src, dst + 2 * dstride, dstride, sum);
}

static INLINE uint32_t sum32(const __m256i src) {
    const __m128i src_l = _mm256_extracti128_si256(src, 0);
    const __m128i src_h = _mm256_extracti128_si256(src, 1);
    const __m128i s     = _mm_add_epi32(src_l, src_h);
    __m128i       dst;

    dst = _mm_hadd_epi32(s, s);
    dst = _mm_hadd_epi32(dst, dst);

    return (uint32_t)_mm_cvtsi128_si32(dst);
}

static INLINE uint64_t dist_8x8_16bit_avx2(const uint16_t **src, const uint16_t *dst,
                                           const int32_t dstride, const int32_t coeff_shift) {
    __m256i ss = _mm256_setzero_si256();
    __m256i dd = _mm256_setzero_si256();
    __m256i s2 = _mm256_setzero_si256();
    __m256i sd = _mm256_setzero_si256();
    __m256i d2 = _mm256_setzero_si256();
    __m256i ssdd;
    __m128i sum;

    for (int32_t r = 0; r < 4; r++) {
        const __m256i s  = _mm256_loadu_si256((const __m256i *)*src);
        const __m128i d0 = _mm_loadu_si128((const __m128i *)(dst + 2 * r * dstride + 0 * dstride));
        const __m128i d1 = _mm_loadu_si128((const __m128i *)(dst + 2 * r * dstride + 1 * dstride));
        const __m256i d  = _mm256_setr_m128i(d0, d1);
        ss               = _mm256_add_epi16(ss, s);
        dd               = _mm256_add_epi16(dd, d);
        s2               = _mm256_add_epi32(s2, _mm256_madd_epi16(s, s));
        sd               = _mm256_add_epi32(sd, _mm256_madd_epi16(s, d));
        d2               = _mm256_add_epi32(d2, _mm256_madd_epi16(d, d));
        *src += 16;
    }

    ssdd                 = _mm256_hadd_epi16(ss, dd);
    ssdd                 = _mm256_hadd_epi16(ssdd, ssdd);
    ssdd                 = _mm256_unpacklo_epi16(ssdd, _mm256_setzero_si256());
    const __m128i ssdd_l = _mm256_extracti128_si256(ssdd, 0);
    const __m128i ssdd_h = _mm256_extracti128_si256(ssdd, 1);
    sum                  = _mm_add_epi32(ssdd_l, ssdd_h);
    sum                  = _mm_hadd_epi32(sum, sum);

    /* Compute the variance -- the calculation cannot go negative. */
    uint64_t sum_s  = _mm_cvtsi128_si32(sum);
    uint64_t sum_d  = _mm_extract_epi32(sum, 1);
    uint64_t sum_s2 = sum32(s2);
    uint64_t sum_d2 = sum32(d2);
    uint64_t sum_sd = sum32(sd);

    /* Compute the variance -- the calculation cannot go negative. */
    uint64_t svar = sum_s2 - ((sum_s * sum_s + 32) >> 6);
    uint64_t dvar = sum_d2 - ((sum_d * sum_d + 32) >> 6);
    return (uint64_t)floor(.5 + (sum_d2 + sum_s2 - 2 * sum_sd) * .5 *
                                    (svar + dvar + (400 << 2 * coeff_shift)) /
                                    (sqrt((20000 << 4 * coeff_shift) + svar * (double)dvar)));
}

static INLINE uint64_t dist_8x8_8bit_avx2(const uint8_t **src, const uint8_t *dst,
                                          const int32_t dstride, const int32_t coeff_shift) {
    __m256i ss = _mm256_setzero_si256();
    __m256i dd = _mm256_setzero_si256();
    __m256i s2 = _mm256_setzero_si256();
    __m256i sd = _mm256_setzero_si256();
    __m256i d2 = _mm256_setzero_si256();
    __m256i ssdd;
    __m128i sum;

    for (int32_t r = 0; r < 4; r++) {
        const __m128i s = _mm_loadu_si128((const __m128i *)*src);
        const __m128i d = _mm_set_epi64x(*(uint64_t *)(dst + 2 * r * dstride + 1 * dstride),
                                         *(uint64_t *)(dst + 2 * r * dstride + 0 * dstride));

        const __m256i s_16 = _mm256_cvtepu8_epi16(s);
        const __m256i d_16 = _mm256_cvtepu8_epi16(d);

        ss = _mm256_add_epi16(ss, s_16);
        dd = _mm256_add_epi16(dd, d_16);
        s2 = _mm256_add_epi32(s2, _mm256_madd_epi16(s_16, s_16));
        sd = _mm256_add_epi32(sd, _mm256_madd_epi16(s_16, d_16));
        d2 = _mm256_add_epi32(d2, _mm256_madd_epi16(d_16, d_16));
        *src += 16;
    }

    ssdd                 = _mm256_hadd_epi16(ss, dd);
    ssdd                 = _mm256_hadd_epi16(ssdd, ssdd);
    ssdd                 = _mm256_unpacklo_epi16(ssdd, _mm256_setzero_si256());
    const __m128i ssdd_l = _mm256_extracti128_si256(ssdd, 0);
    const __m128i ssdd_h = _mm256_extracti128_si256(ssdd, 1);
    sum                  = _mm_add_epi32(ssdd_l, ssdd_h);
    sum                  = _mm_hadd_epi32(sum, sum);

    /* Compute the variance -- the calculation cannot go negative. */
    uint64_t sum_s  = _mm_cvtsi128_si32(sum);
    uint64_t sum_d  = _mm_extract_epi32(sum, 1);
    uint64_t sum_s2 = sum32(s2);
    uint64_t sum_d2 = sum32(d2);
    uint64_t sum_sd = sum32(sd);

    /* Compute the variance -- the calculation cannot go negative. */
    uint64_t svar = sum_s2 - ((sum_s * sum_s + 32) >> 6);
    uint64_t dvar = sum_d2 - ((sum_d * sum_d + 32) >> 6);
    return (uint64_t)floor(.5 + (sum_d2 + sum_s2 - 2 * sum_sd) * .5 *
                                    (svar + dvar + (400 << 2 * coeff_shift)) /
                                    (sqrt((20000 << 4 * coeff_shift) + svar * (double)dvar)));
}

static INLINE void sum_32_to_64(const __m256i src, __m256i *dst) {
    const __m256i src_l = _mm256_unpacklo_epi32(src, _mm256_setzero_si256());
    const __m256i src_h = _mm256_unpackhi_epi32(src, _mm256_setzero_si256());
    *dst                = _mm256_add_epi64(*dst, src_l);
    *dst                = _mm256_add_epi64(*dst, src_h);
}

static INLINE uint64_t sum64(const __m256i src) {
    const __m128i src_l = _mm256_extracti128_si256(src, 0);
    const __m128i src_h = _mm256_extracti128_si256(src, 1);
    const __m128i s     = _mm_add_epi64(src_l, src_h);
    const __m128i dst   = _mm_add_epi64(s, _mm_srli_si128(s, 8));

    return (uint64_t)_mm_cvtsi128_si64(dst);
}

/* Compute MSE only on the blocks we filtered. */
uint64_t compute_cdef_dist_avx2(const uint16_t *dst, int32_t dstride, const uint16_t *src,
                                const CdefList *dlist, int32_t cdef_count, BlockSize bsize,
                                int32_t coeff_shift, int32_t pli) {
    uint64_t sum;
    int32_t  bi, bx, by;

    if ((bsize == BLOCK_8X8) && (pli == 0)) {
        sum = 0;
        for (bi = 0; bi < cdef_count; bi++) {
            by = dlist[bi].by;
            bx = dlist[bi].bx;
            sum += dist_8x8_16bit_avx2(&src, dst + 8 * by * dstride + 8 * bx, dstride, coeff_shift);
        }
    } else {
        __m256i mse64 = _mm256_setzero_si256();

        if (bsize == BLOCK_8X8) {
            for (bi = 0; bi < cdef_count; bi++) {
                __m256i mse32 = _mm256_setzero_si256();
                by            = dlist[bi].by;
                bx            = dlist[bi].bx;
                mse_8x4_16bit_avx2(&src, dst + (8 * by + 0) * dstride + 8 * bx, dstride, &mse32);
                mse_8x4_16bit_avx2(&src, dst + (8 * by + 4) * dstride + 8 * bx, dstride, &mse32);
                sum_32_to_64(mse32, &mse64);
            }
        } else if (bsize == BLOCK_4X8) {
            for (bi = 0; bi < cdef_count; bi++) {
                __m256i mse32 = _mm256_setzero_si256();
                by            = dlist[bi].by;
                bx            = dlist[bi].bx;
                mse_4x4_16bit_avx2(&src, dst + (8 * by + 0) * dstride + 4 * bx, dstride, &mse32);
                mse_4x4_16bit_avx2(&src, dst + (8 * by + 4) * dstride + 4 * bx, dstride, &mse32);
                sum_32_to_64(mse32, &mse64);
            }
        } else if (bsize == BLOCK_8X4) {
            for (bi = 0; bi < cdef_count; bi++) {
                __m256i mse32 = _mm256_setzero_si256();
                by            = dlist[bi].by;
                bx            = dlist[bi].bx;
                mse_8x4_16bit_avx2(&src, dst + 4 * by * dstride + 8 * bx, dstride, &mse32);
                sum_32_to_64(mse32, &mse64);
            }
        } else {
            assert(bsize == BLOCK_4X4);
            for (bi = 0; bi < cdef_count; bi++) {
                __m256i mse32 = _mm256_setzero_si256();
                by            = dlist[bi].by;
                bx            = dlist[bi].bx;
                mse_4x4_16bit_avx2(&src, dst + 4 * by * dstride + 4 * bx, dstride, &mse32);
                sum_32_to_64(mse32, &mse64);
            }
        }

        sum = sum64(mse64);
    }

    return sum >> 2 * coeff_shift;
}

uint64_t compute_cdef_dist_8bit_avx2(const uint8_t *dst8, int32_t dstride, const uint8_t *src8,
                                     const CdefList *dlist, int32_t cdef_count, BlockSize bsize,
                                     int32_t coeff_shift, int32_t pli) {
    uint64_t sum;
    int32_t  bi, bx, by;

    if ((bsize == BLOCK_8X8) && (pli == 0)) {
        sum = 0;
        for (bi = 0; bi < cdef_count; bi++) {
            by = dlist[bi].by;
            bx = dlist[bi].bx;
            sum +=
                dist_8x8_8bit_avx2(&src8, dst8 + 8 * by * dstride + 8 * bx, dstride, coeff_shift);
        }
    } else {
        __m256i mse64 = _mm256_setzero_si256();

        if (bsize == BLOCK_8X8) {
            for (bi = 0; bi < cdef_count; bi++) {
                __m256i mse32 = _mm256_setzero_si256();
                by            = dlist[bi].by;
                bx            = dlist[bi].bx;
                mse_8x4_8bit_avx2(&src8, dst8 + (8 * by + 0) * dstride + 8 * bx, dstride, &mse32);
                mse_8x4_8bit_avx2(&src8, dst8 + (8 * by + 4) * dstride + 8 * bx, dstride, &mse32);
                sum_32_to_64(mse32, &mse64);
            }
        } else if (bsize == BLOCK_4X8) {
            for (bi = 0; bi < cdef_count; bi++) {
                __m256i mse32 = _mm256_setzero_si256();
                by            = dlist[bi].by;
                bx            = dlist[bi].bx;
                mse_4x4_8bit_avx2(&src8, dst8 + (8 * by + 0) * dstride + 4 * bx, dstride, &mse32);
                mse_4x4_8bit_avx2(&src8, dst8 + (8 * by + 4) * dstride + 4 * bx, dstride, &mse32);
                sum_32_to_64(mse32, &mse64);
            }
        } else if (bsize == BLOCK_8X4) {
            for (bi = 0; bi < cdef_count; bi++) {
                __m256i mse32 = _mm256_setzero_si256();
                by            = dlist[bi].by;
                bx            = dlist[bi].bx;
                mse_8x4_8bit_avx2(&src8, dst8 + 4 * by * dstride + 8 * bx, dstride, &mse32);
                sum_32_to_64(mse32, &mse64);
            }
        } else {
            assert(bsize == BLOCK_4X4);
            for (bi = 0; bi < cdef_count; bi++) {
                __m256i mse32 = _mm256_setzero_si256();
                by            = dlist[bi].by;
                bx            = dlist[bi].bx;
                mse_4x4_8bit_avx2(&src8, dst8 + 4 * by * dstride + 4 * bx, dstride, &mse32);
                sum_32_to_64(mse32, &mse64);
            }
        }

        sum = sum64(mse64);
    }
    return sum >> 2 * coeff_shift;
}
