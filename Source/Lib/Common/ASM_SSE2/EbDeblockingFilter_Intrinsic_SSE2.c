/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/ 

#include "emmintrin.h"
#include "EbDefinitions.h"

#include "EbDeblockingFilter_SSE2.h"


/***********************************************************************************/
// synonyms.h
/**
* Various reusable shorthands for x86 SIMD intrinsics.
*
* Intrinsics prefixed with xx_ operate on or return 128bit XMM registers.
* Intrinsics prefixed with yy_ operate on or return 256bit YMM registers.
*/
// Loads and stores to do away with the tedium of casting the address
// to the right type.

static INLINE void xx_storel_32(void *const a, const __m128i v) {
    *(uint32_t *)a = _mm_cvtsi128_si32(v);
}

/***********************************************************************************/
// lpf_common_sse2.h


static INLINE void highbd_transpose6x6_sse2(__m128i *x0, __m128i *x1,
    __m128i *x2, __m128i *x3,
    __m128i *x4, __m128i *x5,
    __m128i *d0, __m128i *d1,
    __m128i *d2, __m128i *d3,
    __m128i *d4, __m128i *d5) {
    __m128i w0, w1, w2, w3, w4, w5, ww0;

    // 00 01 02 03 04 05 xx xx
    // 10 11 12 13 14 15 xx xx
    // 20 21 22 23 24 25 xx xx
    // 30 31 32 33 34 35 xx xx
    // 40 41 42 43 44 45 xx xx
    // 50 51 52 53 54 55 xx xx

    w0 = _mm_unpacklo_epi16(*x0, *x1);  // 00 10 01 11 02 12 03 13
    w1 = _mm_unpacklo_epi16(*x2, *x3);  // 20 30 21 31 22 32 23 33
    w2 = _mm_unpacklo_epi16(*x4, *x5);  // 40 50 41 51 42 52 43 53

    ww0 = _mm_unpacklo_epi32(w0, w1);   // 00 10 20 30 01 11 21 31
    *d0 = _mm_unpacklo_epi64(ww0, w2);  // 00 10 20 30 40 50 41 51
    *d1 = _mm_unpackhi_epi64(ww0,
        _mm_srli_si128(w2, 4));  // 01 11 21 31 41 51 xx xx

    ww0 = _mm_unpackhi_epi32(w0, w1);  // 02 12 22 32 03 13 23 33
    *d2 = _mm_unpacklo_epi64(ww0,
        _mm_srli_si128(w2, 8));  // 02 12 22 32 42 52 xx xx

    w3 = _mm_unpackhi_epi16(*x0, *x1);  // 04 14 05 15 xx xx xx xx
    w4 = _mm_unpackhi_epi16(*x2, *x3);  // 24 34 25 35 xx xx xx xx
    w5 = _mm_unpackhi_epi16(*x4, *x5);  // 44 54 45 55 xx xx xx xx

    *d3 = _mm_unpackhi_epi64(ww0, _mm_srli_si128(w2, 4));  // 03 13 23 33 43 53

    ww0 = _mm_unpacklo_epi32(w3, w4);   //  04 14 24 34 05 15 25 35
    *d4 = _mm_unpacklo_epi64(ww0, w5);  //  04 14 24 34 44 54 45 55
    *d5 = _mm_unpackhi_epi64(ww0,
        _mm_slli_si128(w5, 4));  // 05 15 25 35 45 55 xx xx
}

static INLINE void highbd_transpose4x8_8x4_sse2(__m128i *x0, __m128i *x1,
    __m128i *x2, __m128i *x3,
    __m128i *d0, __m128i *d1,
    __m128i *d2, __m128i *d3,
    __m128i *d4, __m128i *d5,
    __m128i *d6, __m128i *d7) {
    __m128i w0, w1, ww0, ww1, ww2, ww3;
    // input
    // x0 00 01 02 03 04 05 06 07
    // x1 10 11 12 13 14 15 16 17
    // x2 20 21 22 23 24 25 26 27
    // x3 30 31 32 33 34 35 36 37
    // output
    // 00 10 20 30 xx xx xx xx
    // 01 11 21 31 xx xx xx xx
    // 02 12 22 32 xx xx xx xx
    // 03 13 23 33 xx xx xx xx
    // 04 14 24 34 xx xx xx xx
    // 05 15 25 35 xx xx xx xx
    // 06 16 26 36 xx xx xx xx
    // 07 17 27 37 xx xx xx xx

    __m128i zero = _mm_setzero_si128();

    w0 = _mm_unpacklo_epi16(*x0, *x1);  // 00 10 01 11 02 12 03 13
    w1 = _mm_unpacklo_epi16(*x2, *x3);  // 20 30 21 31 22 32 23 33

    ww0 = _mm_unpacklo_epi32(w0, w1);  // 00 10 20 30 01 11 21 31
    ww1 = _mm_unpackhi_epi32(w0, w1);  // 02 12 22 32 03 13 23 33

    w0 = _mm_unpackhi_epi16(*x0, *x1);  // 04 14 05 15 06 16 07 17
    w1 = _mm_unpackhi_epi16(*x2, *x3);  // 24 34 25 35 26 36 27 37

    ww2 = _mm_unpacklo_epi32(w0, w1);  //  04 14 24 34 05 15 25 35
    ww3 = _mm_unpackhi_epi32(w0, w1);  //  06 16 26 36 07 17 27 37

    *d0 = _mm_unpacklo_epi64(ww0, zero);  // 00 10 20 30 xx xx xx xx
    *d1 = _mm_unpackhi_epi64(ww0, zero);  // 01 11 21 31 xx xx xx xx
    *d2 = _mm_unpacklo_epi64(ww1, zero);  // 02 12 22 32 xx xx xx xx
    *d3 = _mm_unpackhi_epi64(ww1, zero);  // 03 13 23 33 xx xx xx xx
    *d4 = _mm_unpacklo_epi64(ww2, zero);  // 04 14 24 34 xx xx xx xx
    *d5 = _mm_unpackhi_epi64(ww2, zero);  // 05 15 25 35 xx xx xx xx
    *d6 = _mm_unpacklo_epi64(ww3, zero);  // 06 16 26 36 xx xx xx xx
    *d7 = _mm_unpackhi_epi64(ww3, zero);  // 07 17 27 37 xx xx xx xx
}

static INLINE void highbd_transpose8x8_low_sse2(__m128i *x0, __m128i *x1,
    __m128i *x2, __m128i *x3,
    __m128i *x4, __m128i *x5,
    __m128i *x6, __m128i *x7,
    __m128i *d0, __m128i *d1,
    __m128i *d2, __m128i *d3) {
    __m128i w0, w1, w2, w3, ww0, ww1;
    // x0 00 01 02 03 04 05 06 07
    // x1 10 11 12 13 14 15 16 17
    // x2 20 21 22 23 24 25 26 27
    // x3 30 31 32 33 34 35 36 37
    // x4 40 41 42 43 44 45 46 47
    // x5 50 51 52 53 54 55 56 57
    // x6 60 61 62 63 64 65 66 67
    // x7 70 71 72 73 74 75 76 77

    w0 = _mm_unpacklo_epi16(*x0, *x1);  // 00 10 01 11 02 12 03 13
    w1 = _mm_unpacklo_epi16(*x2, *x3);  // 20 30 21 31 22 32 23 33
    w2 = _mm_unpacklo_epi16(*x4, *x5);  // 40 50 41 51 42 52 43 53
    w3 = _mm_unpacklo_epi16(*x6, *x7);  // 60 70 61 71 62 72 63 73

    ww0 = _mm_unpacklo_epi32(w0, w1);  // 00 10 20 30 01 11 21 31
    ww1 = _mm_unpacklo_epi32(w2, w3);  // 40 50 60 70 41 51 61 71

    *d0 = _mm_unpacklo_epi64(ww0, ww1);  // 00 10 20 30 40 50 60 70
    *d1 = _mm_unpackhi_epi64(ww0, ww1);  // 01 11 21 31 41 51 61 71

    ww0 = _mm_unpackhi_epi32(w0, w1);  // 02 12 22 32 03 13 23 33
    ww1 = _mm_unpackhi_epi32(w2, w3);  // 42 52 62 72 43 53 63 73

    *d2 = _mm_unpacklo_epi64(ww0, ww1);  // 02 12 22 32 42 52 62 72
    *d3 = _mm_unpackhi_epi64(ww0, ww1);  // 03 13 23 33 43 53 63 73
}

static INLINE void highbd_transpose8x8_high_sse2(__m128i *x0, __m128i *x1,
    __m128i *x2, __m128i *x3,
    __m128i *x4, __m128i *x5,
    __m128i *x6, __m128i *x7,
    __m128i *d4, __m128i *d5,
    __m128i *d6, __m128i *d7) {
    __m128i w0, w1, w2, w3, ww0, ww1;
    // x0 00 01 02 03 04 05 06 07
    // x1 10 11 12 13 14 15 16 17
    // x2 20 21 22 23 24 25 26 27
    // x3 30 31 32 33 34 35 36 37
    // x4 40 41 42 43 44 45 46 47
    // x5 50 51 52 53 54 55 56 57
    // x6 60 61 62 63 64 65 66 67
    // x7 70 71 72 73 74 75 76 77
    w0 = _mm_unpackhi_epi16(*x0, *x1);  // 04 14 05 15 06 16 07 17
    w1 = _mm_unpackhi_epi16(*x2, *x3);  // 24 34 25 35 26 36 27 37
    w2 = _mm_unpackhi_epi16(*x4, *x5);  // 44 54 45 55 46 56 47 57
    w3 = _mm_unpackhi_epi16(*x6, *x7);  // 64 74 65 75 66 76 67 77

    ww0 = _mm_unpacklo_epi32(w0, w1);  // 04 14 24 34 05 15 25 35
    ww1 = _mm_unpacklo_epi32(w2, w3);  // 44 54 64 74 45 55 65 75

    *d4 = _mm_unpacklo_epi64(ww0, ww1);  // 04 14 24 34 44 54 64 74
    *d5 = _mm_unpackhi_epi64(ww0, ww1);  // 05 15 25 35 45 55 65 75

    ww0 = _mm_unpackhi_epi32(w0, w1);  // 06 16 26 36 07 17 27 37
    ww1 = _mm_unpackhi_epi32(w2, w3);  // 46 56 66 76 47 57 67 77

    *d6 = _mm_unpacklo_epi64(ww0, ww1);  // 06 16 26 36 46 56 66 76
    *d7 = _mm_unpackhi_epi64(ww0, ww1);  // 07 17 27 37 47 57 67 77
}

// here in and out pointers (x and d) should be different! we don't store their
// values inside
static INLINE void highbd_transpose8x8_sse2(
    __m128i *x0, __m128i *x1, __m128i *x2, __m128i *x3, __m128i *x4,
    __m128i *x5, __m128i *x6, __m128i *x7, __m128i *d0, __m128i *d1,
    __m128i *d2, __m128i *d3, __m128i *d4, __m128i *d5, __m128i *d6,
    __m128i *d7) {
    highbd_transpose8x8_low_sse2(x0, x1, x2, x3, x4, x5, x6, x7, d0, d1, d2, d3);
    highbd_transpose8x8_high_sse2(x0, x1, x2, x3, x4, x5, x6, x7, d4, d5, d6, d7);
}

/***********************************************************************************/
// loopfilter_sse2.c


static INLINE __m128i abs_diff(__m128i a, __m128i b) {
    return _mm_or_si128(_mm_subs_epu8(a, b), _mm_subs_epu8(b, a));
}

// filter_mask and hev_mask
#define FILTER_HEV_MASK4                                                      \
  do {                                                                        \
    /* (abs(q1 - q0), abs(p1 - p0) */                                         \
    __m128i flat = abs_diff(q1p1, q0p0);                                      \
    /* abs(p1 - q1), abs(p0 - q0) */                                          \
    const __m128i abs_p1q1p0q0 = abs_diff(p1p0, q1q0);                        \
    __m128i abs_p0q0, abs_p1q1;                                               \
                                                                              \
    /* const uint8_t hev = hev_mask(thresh, *op1, *op0, *oq0, *oq1); */       \
    hev =                                                                     \
        _mm_unpacklo_epi8(_mm_max_epu8(flat, _mm_srli_si128(flat, 8)), zero); \
    hev = _mm_cmpgt_epi16(hev, thresh);                                       \
    hev = _mm_packs_epi16(hev, hev);                                          \
                                                                              \
    /* const int8_t mask = filter_mask2(*limit, *blimit, */                   \
    /*                                  p1, p0, q0, q1); */                   \
    abs_p0q0 =                                                                \
        _mm_adds_epu8(abs_p1q1p0q0, abs_p1q1p0q0); /* abs(p0 - q0) * 2 */     \
    abs_p1q1 =                                                                \
        _mm_unpackhi_epi8(abs_p1q1p0q0, abs_p1q1p0q0); /* abs(p1 - q1) */     \
    abs_p1q1 = _mm_srli_epi16(abs_p1q1, 9);                                   \
    abs_p1q1 = _mm_packs_epi16(abs_p1q1, abs_p1q1); /* abs(p1 - q1) / 2 */    \
    /* abs(p0 - q0) * 2 + abs(p1 - q1) / 2 */                                 \
    mask = _mm_adds_epu8(abs_p0q0, abs_p1q1);                                 \
    flat = _mm_max_epu8(flat, _mm_srli_si128(flat, 8));                       \
    mask = _mm_unpacklo_epi64(mask, flat);                                    \
    mask = _mm_subs_epu8(mask, limit);                                        \
    mask = _mm_cmpeq_epi8(mask, zero);                                        \
    mask = _mm_and_si128(mask, _mm_srli_si128(mask, 8));                      \
    } while (0)

AOM_FORCE_INLINE void filter4_sse2(__m128i *p1p0, __m128i *q1q0, __m128i *hev,
    __m128i *mask, __m128i *qs1qs0,
    __m128i *ps1ps0) {
    const __m128i t3t4 =
        _mm_set_epi8(3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4);
    const __m128i t80 = _mm_set1_epi8((char)0x80);
    __m128i filter, filter2filter1, work;
    __m128i ps1ps0_work, qs1qs0_work;
    const __m128i ff = _mm_cmpeq_epi8(t80, t80);

    ps1ps0_work = _mm_xor_si128(*p1p0, t80); /* ^ 0x80 */
    qs1qs0_work = _mm_xor_si128(*q1q0, t80);

    /* int8_t filter = signed_char_clamp(ps1 - qs1) & hev; */
    work = _mm_subs_epi8(ps1ps0_work, qs1qs0_work);
    filter = _mm_and_si128(_mm_srli_si128(work, 8), *hev);
    /* filter = signed_char_clamp(filter + 3 * (qs0 - ps0)) & mask; */
    filter = _mm_subs_epi8(filter, work);
    filter = _mm_subs_epi8(filter, work);
    filter = _mm_subs_epi8(filter, work);  /* + 3 * (qs0 - ps0) */
    filter = _mm_and_si128(filter, *mask); /* & mask */
    filter = _mm_unpacklo_epi64(filter, filter);

    /* filter1 = signed_char_clamp(filter + 4) >> 3; */
    /* filter2 = signed_char_clamp(filter + 3) >> 3; */
    filter2filter1 = _mm_adds_epi8(filter, t3t4); /* signed_char_clamp */
    filter = _mm_unpackhi_epi8(filter2filter1, filter2filter1);
    filter2filter1 = _mm_unpacklo_epi8(filter2filter1, filter2filter1);
    filter2filter1 = _mm_srai_epi16(filter2filter1, 11); /* >> 3 */
    filter = _mm_srai_epi16(filter, 11);                 /* >> 3 */
    filter2filter1 = _mm_packs_epi16(filter2filter1, filter);

    /* filter = ROUND_POWER_OF_TWO(filter1, 1) & ~hev; */
    filter = _mm_subs_epi8(filter2filter1, ff); /* + 1 */
    filter = _mm_unpacklo_epi8(filter, filter);
    filter = _mm_srai_epi16(filter, 9); /* round */
    filter = _mm_packs_epi16(filter, filter);
    filter = _mm_andnot_si128(*hev, filter);

    *hev = _mm_unpackhi_epi64(filter2filter1, filter);
    filter2filter1 = _mm_unpacklo_epi64(filter2filter1, filter);

    /* signed_char_clamp(qs1 - filter), signed_char_clamp(qs0 - filter1) */
    qs1qs0_work = _mm_subs_epi8(qs1qs0_work, filter2filter1);
    /* signed_char_clamp(ps1 + filter), signed_char_clamp(ps0 + filter2) */
    ps1ps0_work = _mm_adds_epi8(ps1ps0_work, *hev);
    *qs1qs0 = _mm_xor_si128(qs1qs0_work, t80); /* ^ 0x80 */
    *ps1ps0 = _mm_xor_si128(ps1ps0_work, t80); /* ^ 0x80 */
}

void aom_lpf_horizontal_4_sse2(uint8_t *s, int32_t p /* pitch */,
    const uint8_t *_blimit, const uint8_t *_limit,
    const uint8_t *_thresh) {
    const __m128i zero = _mm_setzero_si128();
    const __m128i limit =
        _mm_unpacklo_epi64(_mm_loadl_epi64((const __m128i *)_blimit),
            _mm_loadl_epi64((const __m128i *)_limit));
    const __m128i thresh =
        _mm_unpacklo_epi8(_mm_loadl_epi64((const __m128i *)_thresh), zero);
    __m128i q1p1, q0p0, p1p0, q1q0, ps1ps0, qs1qs0;
    __m128i mask, hev;
    q1p1 = _mm_unpacklo_epi64(_mm_cvtsi32_si128(*(int32_t *)(s - 2 * p)),
        _mm_cvtsi32_si128(*(int32_t *)(s + 1 * p)));
    q0p0 = _mm_unpacklo_epi64(_mm_cvtsi32_si128(*(int32_t *)(s - 1 * p)),
        _mm_cvtsi32_si128(*(int32_t *)(s + 0 * p)));
    p1p0 = _mm_unpacklo_epi64(q0p0, q1p1);
    q1q0 = _mm_unpackhi_epi64(q0p0, q1p1);
    FILTER_HEV_MASK4;
    filter4_sse2(&p1p0, &q1q0, &hev, &mask, &qs1qs0, &ps1ps0);

    xx_storel_32(s - 1 * p, ps1ps0);
    xx_storel_32(s - 2 * p, _mm_srli_si128(ps1ps0, 8));
    xx_storel_32(s + 0 * p, qs1qs0);
    xx_storel_32(s + 1 * p, _mm_srli_si128(qs1qs0, 8));
}

void aom_lpf_vertical_4_sse2(uint8_t *s, int32_t p /* pitch */,
    const uint8_t *_blimit, const uint8_t *_limit,
    const uint8_t *_thresh) {
    const __m128i zero = _mm_setzero_si128();
    const __m128i limit =
        _mm_unpacklo_epi64(_mm_loadl_epi64((const __m128i *)_blimit),
            _mm_loadl_epi64((const __m128i *)_limit));
    const __m128i thresh =
        _mm_unpacklo_epi8(_mm_loadl_epi64((const __m128i *)_thresh), zero);

    __m128i x0, x1, x2, x3;
    __m128i q1p1, q0p0, p1p0, q1q0, ps1ps0, qs1qs0;
    __m128i mask, hev;

    // 00 10 01 11 02 12 03 13 04 14 05 15 06 16 07 17
    q1q0 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(s + 0 * p - 4)),
        _mm_loadl_epi64((__m128i *)(s + 1 * p - 4)));

    // 20 30 21 31 22 32 23 33 24 34 25 35 26 36 27 37
    x1 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(s + 2 * p - 4)),
        _mm_loadl_epi64((__m128i *)(s + 3 * p - 4)));

    x2 = _mm_setzero_si128();
    x3 = _mm_setzero_si128();
    // Transpose 8x8
    // 00 10 20 30 01 11 21 31  02 12 22 32 03 13 23 33
    p1p0 = _mm_unpacklo_epi16(q1q0, x1);
    // 40 50 60 70 41 51 61 71  42 52 62 72 43 53 63 73
    x0 = _mm_unpacklo_epi16(x2, x3);
    // 02 12 22 32 42 52 62 72  03 13 23 33 43 53 63 73
    p1p0 = _mm_unpackhi_epi32(p1p0, x0);
    p1p0 = _mm_unpackhi_epi64(p1p0, _mm_slli_si128(p1p0, 8));  // swap lo and high

    // 04 14 24 34 05 15 25 35  06 16 26 36 07 17 27 37
    q1q0 = _mm_unpackhi_epi16(q1q0, x1);
    // 44 54 64 74 45 55 65 75  46 56 66 76 47 57 67 77
    x2 = _mm_unpackhi_epi16(x2, x3);
    // 04 14 24 34 44 54 64 74  05 15 25 35 45 55 65 75
    q1q0 = _mm_unpacklo_epi32(q1q0, x2);

    q0p0 = _mm_unpacklo_epi64(p1p0, q1q0);
    q1p1 = _mm_unpackhi_epi64(p1p0, q1q0);
    p1p0 = _mm_unpacklo_epi64(q0p0, q1p1);
    FILTER_HEV_MASK4;
    filter4_sse2(&p1p0, &q1q0, &hev, &mask, &qs1qs0, &ps1ps0);

    // Transpose 8x4 to 4x8
    // qs1qs0: 20 21 22 23 24 25 26 27  30 31 32 33 34 34 36 37
    // ps1ps0: 10 11 12 13 14 15 16 17  00 01 02 03 04 05 06 07
    // 00 01 02 03 04 05 06 07  10 11 12 13 14 15 16 17
    ps1ps0 = _mm_unpackhi_epi64(ps1ps0, _mm_slli_si128(ps1ps0, 8));
    // 10 30 11 31 12 32 13 33  14 34 15 35 16 36 17 37
    x0 = _mm_unpackhi_epi8(ps1ps0, qs1qs0);
    // 00 20 01 21 02 22 03 23  04 24 05 25 06 26 07 27
    ps1ps0 = _mm_unpacklo_epi8(ps1ps0, qs1qs0);
    // 00 10 20 30 01 11 21 31  02 12 22 32 03 13 23 33
    ps1ps0 = _mm_unpacklo_epi8(ps1ps0, x0);

    xx_storel_32(s + 0 * p - 2, ps1ps0);
    xx_storel_32(s + 1 * p - 2, _mm_srli_si128(ps1ps0, 4));
    xx_storel_32(s + 2 * p - 2, _mm_srli_si128(ps1ps0, 8));
    xx_storel_32(s + 3 * p - 2, _mm_srli_si128(ps1ps0, 12));
}

static INLINE void store_buffer_horz_8(__m128i x, int32_t p, int32_t num, uint8_t *s) {
    xx_storel_32(s - (num + 1) * p, x);
    xx_storel_32(s + num * p, _mm_srli_si128(x, 8));
}

static AOM_FORCE_INLINE void lpf_internal_14_sse2(
    __m128i *q6p6, __m128i *q5p5, __m128i *q4p4, __m128i *q3p3, __m128i *q2p2,
    __m128i *q1p1, __m128i *q0p0, const uint8_t *_blimit,
    const uint8_t *_limit, const uint8_t *_thresh) {
    const __m128i zero = _mm_setzero_si128();
    const __m128i one = _mm_set1_epi8(1);
    const __m128i blimit = _mm_load_si128((const __m128i *)_blimit);
    const __m128i limit = _mm_load_si128((const __m128i *)_limit);
    const __m128i thresh = _mm_load_si128((const __m128i *)_thresh);
    __m128i mask, hev, flat, flat2;
    __m128i p0q0, p1q1;

    __m128i abs_p1p0;

    p0q0 = _mm_shuffle_epi32(*q0p0, _MM_SHUFFLE(1, 0, 3, 2));
    p1q1 = _mm_shuffle_epi32(*q1p1, _MM_SHUFFLE(1, 0, 3, 2));

    {
        __m128i abs_p1q1, abs_p0q0, abs_q1q0, fe, ff, work;
        abs_p1p0 = abs_diff(*q1p1, *q0p0);
        abs_q1q0 = _mm_srli_si128(abs_p1p0, 8);
        fe = _mm_set1_epi8((char)0xfe);
        ff = _mm_cmpeq_epi8(abs_p1p0, abs_p1p0);
        abs_p0q0 = abs_diff(*q0p0, p0q0);
        abs_p1q1 = abs_diff(*q1p1, p1q1);
        flat = _mm_max_epu8(abs_p1p0, abs_q1q0);
        hev = _mm_subs_epu8(flat, thresh);
        hev = _mm_xor_si128(_mm_cmpeq_epi8(hev, zero), ff);

        abs_p0q0 = _mm_adds_epu8(abs_p0q0, abs_p0q0);
        abs_p1q1 = _mm_srli_epi16(_mm_and_si128(abs_p1q1, fe), 1);
        mask = _mm_subs_epu8(_mm_adds_epu8(abs_p0q0, abs_p1q1), blimit);
        mask = _mm_xor_si128(_mm_cmpeq_epi8(mask, zero), ff);
        // mask |= (abs(p0 - q0) * 2 + abs(p1 - q1) / 2  > blimit) * -1;
        mask = _mm_max_epu8(abs_p1p0, mask);
        // mask |= (abs(p1 - p0) > limit) * -1;
        // mask |= (abs(q1 - q0) > limit) * -1;

        work = _mm_max_epu8(abs_diff(*q2p2, *q1p1), abs_diff(*q3p3, *q2p2));
        mask = _mm_max_epu8(work, mask);
        mask = _mm_max_epu8(mask, _mm_srli_si128(mask, 8));
        mask = _mm_subs_epu8(mask, limit);
        mask = _mm_cmpeq_epi8(mask, zero);
    }

    // lp filter
    {
        const __m128i t4 = _mm_set1_epi8(4);
        const __m128i t3 = _mm_set1_epi8(3);
        const __m128i t80 = _mm_set1_epi8((char)0x80);
        const __m128i t1 = _mm_set1_epi16(0x1);
        __m128i qs1ps1 = _mm_xor_si128(*q1p1, t80);
        __m128i qs0ps0 = _mm_xor_si128(*q0p0, t80);
        __m128i qs0 = _mm_xor_si128(p0q0, t80);
        __m128i qs1 = _mm_xor_si128(p1q1, t80);
        __m128i filt;
        __m128i work_a;
        __m128i filter1, filter2;
        __m128i flat2_q5p5, flat2_q4p4, flat2_q3p3, flat2_q2p2;
        __m128i flat2_q1p1, flat2_q0p0, flat_q2p2, flat_q1p1, flat_q0p0;

        filt = _mm_and_si128(_mm_subs_epi8(qs1ps1, qs1), hev);
        work_a = _mm_subs_epi8(qs0, qs0ps0);
        filt = _mm_adds_epi8(filt, work_a);
        filt = _mm_adds_epi8(filt, work_a);
        filt = _mm_adds_epi8(filt, work_a);
        // (aom_filter + 3 * (qs0 - ps0)) & mask
        filt = _mm_and_si128(filt, mask);

        filter1 = _mm_adds_epi8(filt, t4);
        filter2 = _mm_adds_epi8(filt, t3);

        filter1 = _mm_unpacklo_epi8(zero, filter1);
        filter1 = _mm_srai_epi16(filter1, 0xB);
        filter2 = _mm_unpacklo_epi8(zero, filter2);
        filter2 = _mm_srai_epi16(filter2, 0xB);

        // Filter1 >> 3
        filt = _mm_packs_epi16(filter2, _mm_subs_epi16(zero, filter1));
        qs0ps0 = _mm_xor_si128(_mm_adds_epi8(qs0ps0, filt), t80);

        // filt >> 1
        filt = _mm_adds_epi16(filter1, t1);
        filt = _mm_srai_epi16(filt, 1);
        filt = _mm_andnot_si128(_mm_srai_epi16(_mm_unpacklo_epi8(zero, hev), 0x8),
            filt);
        filt = _mm_packs_epi16(filt, _mm_subs_epi16(zero, filt));
        qs1ps1 = _mm_xor_si128(_mm_adds_epi8(qs1ps1, filt), t80);
        // loopfilter done

        {
            __m128i work;
            flat = _mm_max_epu8(abs_diff(*q2p2, *q0p0), abs_diff(*q3p3, *q0p0));
            flat = _mm_max_epu8(abs_p1p0, flat);
            flat = _mm_max_epu8(flat, _mm_srli_si128(flat, 8));
            flat = _mm_subs_epu8(flat, one);
            flat = _mm_cmpeq_epi8(flat, zero);
            flat = _mm_and_si128(flat, mask);

            flat2 = _mm_max_epu8(abs_diff(*q4p4, *q0p0), abs_diff(*q5p5, *q0p0));
            work = abs_diff(*q6p6, *q0p0);
            flat2 = _mm_max_epu8(work, flat2);
            flat2 = _mm_max_epu8(flat2, _mm_srli_si128(flat2, 8));
            flat2 = _mm_subs_epu8(flat2, one);
            flat2 = _mm_cmpeq_epi8(flat2, zero);
            flat2 = _mm_and_si128(flat2, flat);  // flat2 & flat & mask
        }
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // flat and wide flat calculations
        {
            const __m128i eight = _mm_set1_epi16(8);
            const __m128i four = _mm_set1_epi16(4);
            __m128i p6_16, p5_16, p4_16, p3_16, p2_16, p1_16, p0_16;
            __m128i q6_16, q5_16, q4_16, q3_16, q2_16, q1_16, q0_16;
            __m128i pixelFilter_p, pixelFilter_q;
            __m128i pixetFilter_p2p1p0, pixetFilter_q2q1q0;
            __m128i sum_p6, sum_q6;
            __m128i sum_p3, sum_q3, res_p, res_q;

            p6_16 = _mm_unpacklo_epi8(*q6p6, zero);
            p5_16 = _mm_unpacklo_epi8(*q5p5, zero);
            p4_16 = _mm_unpacklo_epi8(*q4p4, zero);
            p3_16 = _mm_unpacklo_epi8(*q3p3, zero);
            p2_16 = _mm_unpacklo_epi8(*q2p2, zero);
            p1_16 = _mm_unpacklo_epi8(*q1p1, zero);
            p0_16 = _mm_unpacklo_epi8(*q0p0, zero);
            q0_16 = _mm_unpackhi_epi8(*q0p0, zero);
            q1_16 = _mm_unpackhi_epi8(*q1p1, zero);
            q2_16 = _mm_unpackhi_epi8(*q2p2, zero);
            q3_16 = _mm_unpackhi_epi8(*q3p3, zero);
            q4_16 = _mm_unpackhi_epi8(*q4p4, zero);
            q5_16 = _mm_unpackhi_epi8(*q5p5, zero);
            q6_16 = _mm_unpackhi_epi8(*q6p6, zero);
            pixelFilter_p = _mm_add_epi16(p5_16, _mm_add_epi16(p4_16, p3_16));
            pixelFilter_q = _mm_add_epi16(q5_16, _mm_add_epi16(q4_16, q3_16));

            pixetFilter_p2p1p0 = _mm_add_epi16(p0_16, _mm_add_epi16(p2_16, p1_16));
            pixelFilter_p = _mm_add_epi16(pixelFilter_p, pixetFilter_p2p1p0);

            pixetFilter_q2q1q0 = _mm_add_epi16(q0_16, _mm_add_epi16(q2_16, q1_16));
            pixelFilter_q = _mm_add_epi16(pixelFilter_q, pixetFilter_q2q1q0);
            pixelFilter_p =
                _mm_add_epi16(eight, _mm_add_epi16(pixelFilter_p, pixelFilter_q));
            pixetFilter_p2p1p0 = _mm_add_epi16(
                four, _mm_add_epi16(pixetFilter_p2p1p0, pixetFilter_q2q1q0));
            res_p = _mm_srli_epi16(
                _mm_add_epi16(pixelFilter_p,
                    _mm_add_epi16(_mm_add_epi16(p6_16, p0_16),
                        _mm_add_epi16(p1_16, q0_16))),
                4);
            res_q = _mm_srli_epi16(
                _mm_add_epi16(pixelFilter_p,
                    _mm_add_epi16(_mm_add_epi16(q6_16, q0_16),
                        _mm_add_epi16(p0_16, q1_16))),
                4);
            flat2_q0p0 = _mm_packus_epi16(res_p, res_q);

            res_p = _mm_srli_epi16(
                _mm_add_epi16(pixetFilter_p2p1p0, _mm_add_epi16(p3_16, p0_16)), 3);
            res_q = _mm_srli_epi16(
                _mm_add_epi16(pixetFilter_p2p1p0, _mm_add_epi16(q3_16, q0_16)), 3);

            flat_q0p0 = _mm_packus_epi16(res_p, res_q);

            sum_p6 = _mm_add_epi16(p6_16, p6_16);
            sum_q6 = _mm_add_epi16(q6_16, q6_16);
            sum_p3 = _mm_add_epi16(p3_16, p3_16);
            sum_q3 = _mm_add_epi16(q3_16, q3_16);

            pixelFilter_q = _mm_sub_epi16(pixelFilter_p, p5_16);
            pixelFilter_p = _mm_sub_epi16(pixelFilter_p, q5_16);

            res_p = _mm_srli_epi16(
                _mm_add_epi16(
                    pixelFilter_p,
                    _mm_add_epi16(sum_p6,
                        _mm_add_epi16(p1_16, _mm_add_epi16(p2_16, p0_16)))),
                4);
            res_q = _mm_srli_epi16(
                _mm_add_epi16(
                    pixelFilter_q,
                    _mm_add_epi16(sum_q6,
                        _mm_add_epi16(q1_16, _mm_add_epi16(q0_16, q2_16)))),
                4);
            flat2_q1p1 = _mm_packus_epi16(res_p, res_q);

            pixetFilter_q2q1q0 = _mm_sub_epi16(pixetFilter_p2p1p0, p2_16);
            pixetFilter_p2p1p0 = _mm_sub_epi16(pixetFilter_p2p1p0, q2_16);
            res_p = _mm_srli_epi16(
                _mm_add_epi16(pixetFilter_p2p1p0, _mm_add_epi16(sum_p3, p1_16)), 3);
            res_q = _mm_srli_epi16(
                _mm_add_epi16(pixetFilter_q2q1q0, _mm_add_epi16(sum_q3, q1_16)), 3);
            flat_q1p1 = _mm_packus_epi16(res_p, res_q);

            sum_p6 = _mm_add_epi16(sum_p6, p6_16);
            sum_q6 = _mm_add_epi16(sum_q6, q6_16);
            sum_p3 = _mm_add_epi16(sum_p3, p3_16);
            sum_q3 = _mm_add_epi16(sum_q3, q3_16);

            pixelFilter_p = _mm_sub_epi16(pixelFilter_p, q4_16);
            pixelFilter_q = _mm_sub_epi16(pixelFilter_q, p4_16);

            res_p = _mm_srli_epi16(
                _mm_add_epi16(
                    pixelFilter_p,
                    _mm_add_epi16(sum_p6,
                        _mm_add_epi16(p2_16, _mm_add_epi16(p3_16, p1_16)))),
                4);
            res_q = _mm_srli_epi16(
                _mm_add_epi16(
                    pixelFilter_q,
                    _mm_add_epi16(sum_q6,
                        _mm_add_epi16(q2_16, _mm_add_epi16(q1_16, q3_16)))),
                4);
            flat2_q2p2 = _mm_packus_epi16(res_p, res_q);

            pixetFilter_p2p1p0 = _mm_sub_epi16(pixetFilter_p2p1p0, q1_16);
            pixetFilter_q2q1q0 = _mm_sub_epi16(pixetFilter_q2q1q0, p1_16);

            res_p = _mm_srli_epi16(
                _mm_add_epi16(pixetFilter_p2p1p0, _mm_add_epi16(sum_p3, p2_16)), 3);
            res_q = _mm_srli_epi16(
                _mm_add_epi16(pixetFilter_q2q1q0, _mm_add_epi16(sum_q3, q2_16)), 3);
            flat_q2p2 = _mm_packus_epi16(res_p, res_q);

            sum_p6 = _mm_add_epi16(sum_p6, p6_16);
            sum_q6 = _mm_add_epi16(sum_q6, q6_16);

            pixelFilter_p = _mm_sub_epi16(pixelFilter_p, q3_16);
            pixelFilter_q = _mm_sub_epi16(pixelFilter_q, p3_16);

            res_p = _mm_srli_epi16(
                _mm_add_epi16(
                    pixelFilter_p,
                    _mm_add_epi16(sum_p6,
                        _mm_add_epi16(p3_16, _mm_add_epi16(p4_16, p2_16)))),
                4);
            res_q = _mm_srli_epi16(
                _mm_add_epi16(
                    pixelFilter_q,
                    _mm_add_epi16(sum_q6,
                        _mm_add_epi16(q3_16, _mm_add_epi16(q2_16, q4_16)))),
                4);
            flat2_q3p3 = _mm_packus_epi16(res_p, res_q);

            sum_p6 = _mm_add_epi16(sum_p6, p6_16);
            sum_q6 = _mm_add_epi16(sum_q6, q6_16);

            pixelFilter_p = _mm_sub_epi16(pixelFilter_p, q2_16);
            pixelFilter_q = _mm_sub_epi16(pixelFilter_q, p2_16);

            res_p = _mm_srli_epi16(
                _mm_add_epi16(
                    pixelFilter_p,
                    _mm_add_epi16(sum_p6,
                        _mm_add_epi16(p4_16, _mm_add_epi16(p5_16, p3_16)))),
                4);
            res_q = _mm_srli_epi16(
                _mm_add_epi16(
                    pixelFilter_q,
                    _mm_add_epi16(sum_q6,
                        _mm_add_epi16(q4_16, _mm_add_epi16(q3_16, q5_16)))),
                4);
            flat2_q4p4 = _mm_packus_epi16(res_p, res_q);

            sum_p6 = _mm_add_epi16(sum_p6, p6_16);
            sum_q6 = _mm_add_epi16(sum_q6, q6_16);
            pixelFilter_p = _mm_sub_epi16(pixelFilter_p, q1_16);
            pixelFilter_q = _mm_sub_epi16(pixelFilter_q, p1_16);

            res_p = _mm_srli_epi16(
                _mm_add_epi16(
                    pixelFilter_p,
                    _mm_add_epi16(sum_p6,
                        _mm_add_epi16(p5_16, _mm_add_epi16(p6_16, p4_16)))),
                4);
            res_q = _mm_srli_epi16(
                _mm_add_epi16(
                    pixelFilter_q,
                    _mm_add_epi16(sum_q6,
                        _mm_add_epi16(q5_16, _mm_add_epi16(q6_16, q4_16)))),
                4);
            flat2_q5p5 = _mm_packus_epi16(res_p, res_q);
        }
        // wide flat
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        flat = _mm_shuffle_epi32(flat, 68);
        flat2 = _mm_shuffle_epi32(flat2, 68);

        *q2p2 = _mm_andnot_si128(flat, *q2p2);
        flat_q2p2 = _mm_and_si128(flat, flat_q2p2);
        *q2p2 = _mm_or_si128(*q2p2, flat_q2p2);

        qs1ps1 = _mm_andnot_si128(flat, qs1ps1);
        flat_q1p1 = _mm_and_si128(flat, flat_q1p1);
        *q1p1 = _mm_or_si128(qs1ps1, flat_q1p1);

        qs0ps0 = _mm_andnot_si128(flat, qs0ps0);
        flat_q0p0 = _mm_and_si128(flat, flat_q0p0);
        *q0p0 = _mm_or_si128(qs0ps0, flat_q0p0);

        *q5p5 = _mm_andnot_si128(flat2, *q5p5);
        flat2_q5p5 = _mm_and_si128(flat2, flat2_q5p5);
        *q5p5 = _mm_or_si128(*q5p5, flat2_q5p5);

        *q4p4 = _mm_andnot_si128(flat2, *q4p4);
        flat2_q4p4 = _mm_and_si128(flat2, flat2_q4p4);
        *q4p4 = _mm_or_si128(*q4p4, flat2_q4p4);

        *q3p3 = _mm_andnot_si128(flat2, *q3p3);
        flat2_q3p3 = _mm_and_si128(flat2, flat2_q3p3);
        *q3p3 = _mm_or_si128(*q3p3, flat2_q3p3);

        *q2p2 = _mm_andnot_si128(flat2, *q2p2);
        flat2_q2p2 = _mm_and_si128(flat2, flat2_q2p2);
        *q2p2 = _mm_or_si128(*q2p2, flat2_q2p2);

        *q1p1 = _mm_andnot_si128(flat2, *q1p1);
        flat2_q1p1 = _mm_and_si128(flat2, flat2_q1p1);
        *q1p1 = _mm_or_si128(*q1p1, flat2_q1p1);

        *q0p0 = _mm_andnot_si128(flat2, *q0p0);
        flat2_q0p0 = _mm_and_si128(flat2, flat2_q0p0);
        *q0p0 = _mm_or_si128(*q0p0, flat2_q0p0);
    }
}

void aom_lpf_horizontal_14_sse2(uint8_t *s, int32_t p,
    const uint8_t *_blimit,
    const uint8_t *_limit,
    const uint8_t *_thresh) {
    __m128i q6p6, q5p5, q4p4, q3p3, q2p2, q1p1, q0p0;

    q4p4 = _mm_unpacklo_epi64(_mm_cvtsi32_si128(*(int32_t *)(s - 5 * p)),
        _mm_cvtsi32_si128(*(int32_t *)(s + 4 * p)));
    q3p3 = _mm_unpacklo_epi64(_mm_cvtsi32_si128(*(int32_t *)(s - 4 * p)),
        _mm_cvtsi32_si128(*(int32_t *)(s + 3 * p)));
    q2p2 = _mm_unpacklo_epi64(_mm_cvtsi32_si128(*(int32_t *)(s - 3 * p)),
        _mm_cvtsi32_si128(*(int32_t *)(s + 2 * p)));
    q1p1 = _mm_unpacklo_epi64(_mm_cvtsi32_si128(*(int32_t *)(s - 2 * p)),
        _mm_cvtsi32_si128(*(int32_t *)(s + 1 * p)));

    q0p0 = _mm_unpacklo_epi64(_mm_cvtsi32_si128(*(int32_t *)(s - 1 * p)),
        _mm_cvtsi32_si128(*(int32_t *)(s - 0 * p)));

    q5p5 = _mm_unpacklo_epi64(_mm_cvtsi32_si128(*(int32_t *)(s - 6 * p)),
        _mm_cvtsi32_si128(*(int32_t *)(s + 5 * p)));

    q6p6 = _mm_unpacklo_epi64(_mm_cvtsi32_si128(*(int32_t *)(s - 7 * p)),
        _mm_cvtsi32_si128(*(int32_t *)(s + 6 * p)));

    lpf_internal_14_sse2(&q6p6, &q5p5, &q4p4, &q3p3, &q2p2, &q1p1, &q0p0, _blimit,
        _limit, _thresh);

    store_buffer_horz_8(q0p0, p, 0, s);
    store_buffer_horz_8(q1p1, p, 1, s);
    store_buffer_horz_8(q2p2, p, 2, s);
    store_buffer_horz_8(q3p3, p, 3, s);
    store_buffer_horz_8(q4p4, p, 4, s);
    store_buffer_horz_8(q5p5, p, 5, s);
}

static AOM_FORCE_INLINE void lpf_internal_6_sse2(
    __m128i *p2, __m128i *q2, __m128i *p1, __m128i *q1, __m128i *p0,
    __m128i *q0, __m128i *q1q0, __m128i *p1p0, const uint8_t *_blimit,
    const uint8_t *_limit, const uint8_t *_thresh) {
    const __m128i zero = _mm_setzero_si128();
    const __m128i blimit = _mm_load_si128((const __m128i *)_blimit);
    const __m128i limit = _mm_load_si128((const __m128i *)_limit);
    const __m128i thresh = _mm_load_si128((const __m128i *)_thresh);
    __m128i mask, hev, flat;
    __m128i q2p2, q1p1, q0p0, p1q1, p0q0, flat_p1p0, flat_q0q1;
    __m128i p2_16, q2_16, p1_16, q1_16, p0_16, q0_16;
    __m128i ps1ps0, qs1qs0;

    q2p2 = _mm_unpacklo_epi64(*p2, *q2);
    q1p1 = _mm_unpacklo_epi64(*p1, *q1);
    q0p0 = _mm_unpacklo_epi64(*p0, *q0);

    p1q1 = _mm_shuffle_epi32(q1p1, _MM_SHUFFLE(1, 0, 3, 2));
    p0q0 = _mm_shuffle_epi32(q0p0, _MM_SHUFFLE(1, 0, 3, 2));

    const __m128i one = _mm_set1_epi8(1);
    const __m128i fe = _mm_set1_epi8((char)0xfe);
    const __m128i ff = _mm_cmpeq_epi8(fe, fe);
    {
        // filter_mask and hev_mask
        __m128i abs_p1q1, abs_p0q0, abs_q1q0, abs_p1p0, work;
        abs_p1p0 = abs_diff(q1p1, q0p0);
        abs_q1q0 = _mm_srli_si128(abs_p1p0, 8);

        abs_p0q0 = abs_diff(q0p0, p0q0);
        abs_p1q1 = abs_diff(q1p1, p1q1);

        // considering sse doesn't have unsigned elements comparison the idea is
        // to find at least one case when X > limit, it means the corresponding
        // mask bit is set.
        // to achieve that we find global max value of all inputs of abs(x-y) or
        // (abs(p0 - q0) * 2 + abs(p1 - q1) / 2 If it is > limit the mask is set
        // otherwise - not

        flat = _mm_max_epu8(abs_p1p0, abs_q1q0);
        hev = _mm_subs_epu8(flat, thresh);
        hev = _mm_xor_si128(_mm_cmpeq_epi8(hev, zero), ff);
        // replicate for the further "merged variables" usage
        hev = _mm_unpacklo_epi64(hev, hev);

        abs_p0q0 = _mm_adds_epu8(abs_p0q0, abs_p0q0);
        abs_p1q1 = _mm_srli_epi16(_mm_and_si128(abs_p1q1, fe), 1);
        mask = _mm_subs_epu8(_mm_adds_epu8(abs_p0q0, abs_p1q1), blimit);
        mask = _mm_xor_si128(_mm_cmpeq_epi8(mask, zero), ff);
        // mask |= (abs(p0 - q0) * 2 + abs(p1 - q1) / 2  > blimit) * -1;
        mask = _mm_max_epu8(abs_p1p0, mask);
        // mask |= (abs(p1 - p0) > limit) * -1;
        // mask |= (abs(q1 - q0) > limit) * -1;

        work = abs_diff(q2p2, q1p1);
        mask = _mm_max_epu8(work, mask);
        mask = _mm_max_epu8(mask, _mm_srli_si128(mask, 8));
        mask = _mm_subs_epu8(mask, limit);
        mask = _mm_cmpeq_epi8(mask, zero);
        // replicate for the further "merged variables" usage
        mask = _mm_unpacklo_epi64(mask, mask);

        // flat_mask
        flat = _mm_max_epu8(abs_diff(q2p2, q0p0), abs_p1p0);
        flat = _mm_max_epu8(flat, _mm_srli_si128(flat, 8));
        flat = _mm_subs_epu8(flat, one);
        flat = _mm_cmpeq_epi8(flat, zero);
        flat = _mm_and_si128(flat, mask);
        // replicate for the further "merged variables" usage
        flat = _mm_unpacklo_epi64(flat, flat);
    }

    // 5 tap filter
    {
        const __m128i four = _mm_set1_epi16(4);

        __m128i workp_a, workp_b, workp_shft0, workp_shft1;
        p2_16 = _mm_unpacklo_epi8(*p2, zero);
        p1_16 = _mm_unpacklo_epi8(*p1, zero);
        p0_16 = _mm_unpacklo_epi8(*p0, zero);
        q0_16 = _mm_unpacklo_epi8(*q0, zero);
        q1_16 = _mm_unpacklo_epi8(*q1, zero);
        q2_16 = _mm_unpacklo_epi8(*q2, zero);

        // op1
        workp_a = _mm_add_epi16(_mm_add_epi16(p0_16, p0_16),
            _mm_add_epi16(p1_16, p1_16));  // p0 *2 + p1 * 2
        workp_a = _mm_add_epi16(_mm_add_epi16(workp_a, four),
            p2_16);  // p2 + p0 * 2 + p1 * 2 + 4

        workp_b = _mm_add_epi16(_mm_add_epi16(p2_16, p2_16), q0_16);
        workp_shft0 = _mm_srli_epi16(_mm_add_epi16(workp_a, workp_b),
            3);  // p2 * 3 + p1 * 2 + p0 * 2 + q0 + 4

        // op0
        workp_b = _mm_add_epi16(_mm_add_epi16(q0_16, q0_16), q1_16);  // q0 * 2 + q1
        workp_a = _mm_add_epi16(workp_a,
            workp_b);  // p2 + p0 * 2 + p1 * 2 + q0 * 2 + q1 + 4
        workp_shft1 = _mm_srli_epi16(workp_a, 3);

        flat_p1p0 = _mm_packus_epi16(workp_shft1, workp_shft0);

        // oq0
        workp_a = _mm_sub_epi16(_mm_sub_epi16(workp_a, p2_16),
            p1_16);  // p0 * 2 + p1  + q0 * 2 + q1 + 4
        workp_b = _mm_add_epi16(q1_16, q2_16);
        workp_a = _mm_add_epi16(
            workp_a, workp_b);  // p0 * 2 + p1  + q0 * 2 + q1 * 2 + q2 + 4
        workp_shft0 = _mm_srli_epi16(workp_a, 3);

        // oq1
        workp_a = _mm_sub_epi16(_mm_sub_epi16(workp_a, p1_16),
            p0_16);  // p0   + q0 * 2 + q1 * 2 + q2 + 4
        workp_b = _mm_add_epi16(q2_16, q2_16);
        workp_shft1 = _mm_srli_epi16(_mm_add_epi16(workp_a, workp_b),
            3);  // p0  + q0 * 2 + q1 * 2 + q2 * 3 + 4

        flat_q0q1 = _mm_packus_epi16(workp_shft0, workp_shft1);
    }

    // lp filter
    *p1p0 = _mm_unpacklo_epi64(q0p0, q1p1);
    *q1q0 = _mm_unpackhi_epi64(q0p0, q1p1);

    filter4_sse2(p1p0, q1q0, &hev, &mask, &qs1qs0, &ps1ps0);

    qs1qs0 = _mm_andnot_si128(flat, qs1qs0);
    *q1q0 = _mm_and_si128(flat, flat_q0q1);
    *q1q0 = _mm_or_si128(qs1qs0, *q1q0);

    ps1ps0 = _mm_andnot_si128(flat, ps1ps0);
    *p1p0 = _mm_and_si128(flat, flat_p1p0);
    *p1p0 = _mm_or_si128(ps1ps0, *p1p0);
}

void aom_lpf_horizontal_6_sse2(uint8_t *s, int32_t p,
    const uint8_t *_blimit,
    const uint8_t *_limit,
    const uint8_t *_thresh) {
    __m128i p2, p1, p0, q0, q1, q2;
    __m128i p1p0, q1q0;

    p2 = _mm_cvtsi32_si128(*(int32_t *)(s - 3 * p));
    p1 = _mm_cvtsi32_si128(*(int32_t *)(s - 2 * p));
    p0 = _mm_cvtsi32_si128(*(int32_t *)(s - 1 * p));
    q0 = _mm_cvtsi32_si128(*(int32_t *)(s - 0 * p));
    q1 = _mm_cvtsi32_si128(*(int32_t *)(s + 1 * p));
    q2 = _mm_cvtsi32_si128(*(int32_t *)(s + 2 * p));

    lpf_internal_6_sse2(&p2, &q2, &p1, &q1, &p0, &q0, &q1q0, &p1p0, _blimit,
        _limit, _thresh);

    xx_storel_32(s - 1 * p, p1p0);
    xx_storel_32(s - 2 * p, _mm_srli_si128(p1p0, 8));
    xx_storel_32(s + 0 * p, q1q0);
    xx_storel_32(s + 1 * p, _mm_srli_si128(q1q0, 8));
}

static AOM_FORCE_INLINE void lpf_internal_8_sse2(
    __m128i *p3_8, __m128i *q3_8, __m128i *p2_8, __m128i *q2_8, __m128i *p1_8,
    __m128i *q1_8, __m128i *p0_8, __m128i *q0_8, __m128i *q1q0_out,
    __m128i *p1p0_out, __m128i *p2_out, __m128i *q2_out,
    const uint8_t *_blimit, const uint8_t *_limit,
    const uint8_t *_thresh) {
    const __m128i zero = _mm_setzero_si128();
    const __m128i blimit = _mm_load_si128((const __m128i *)_blimit);
    const __m128i limit = _mm_load_si128((const __m128i *)_limit);
    const __m128i thresh = _mm_load_si128((const __m128i *)_thresh);
    __m128i mask, hev, flat;
    __m128i p2, q2, p1, p0, q0, q1, p3, q3, q3p3, flat_p1p0, flat_q0q1;
    __m128i q2p2, q1p1, q0p0, p1q1, p0q0;
    __m128i q1q0, p1p0, ps1ps0, qs1qs0;
    __m128i work_a, op2, oq2;

    q3p3 = _mm_unpacklo_epi64(*p3_8, *q3_8);
    q2p2 = _mm_unpacklo_epi64(*p2_8, *q2_8);
    q1p1 = _mm_unpacklo_epi64(*p1_8, *q1_8);
    q0p0 = _mm_unpacklo_epi64(*p0_8, *q0_8);

    p1q1 = _mm_shuffle_epi32(q1p1, _MM_SHUFFLE(1, 0, 3, 2));
    p0q0 = _mm_shuffle_epi32(q0p0, _MM_SHUFFLE(1, 0, 3, 2));

    {
        // filter_mask and hev_mask

        // considering sse doesn't have unsigned elements comparison the idea is to
        // find at least one case when X > limit, it means the corresponding  mask
        // bit is set.
        // to achieve that we find global max value of all inputs of abs(x-y) or
        // (abs(p0 - q0) * 2 + abs(p1 - q1) / 2 If it is > limit the mask is set
        // otherwise - not

        const __m128i one = _mm_set1_epi8(1);
        const __m128i fe = _mm_set1_epi8((char)0xfe);
        const __m128i ff = _mm_cmpeq_epi8(fe, fe);
        __m128i abs_p1q1, abs_p0q0, abs_q1q0, abs_p1p0, work;

        abs_p1p0 = abs_diff(q1p1, q0p0);
        abs_q1q0 = _mm_srli_si128(abs_p1p0, 8);

        abs_p0q0 = abs_diff(q0p0, p0q0);
        abs_p1q1 = abs_diff(q1p1, p1q1);
        flat = _mm_max_epu8(abs_p1p0, abs_q1q0);
        hev = _mm_subs_epu8(flat, thresh);
        hev = _mm_xor_si128(_mm_cmpeq_epi8(hev, zero), ff);
        // replicate for the further "merged variables" usage
        hev = _mm_unpacklo_epi64(hev, hev);

        abs_p0q0 = _mm_adds_epu8(abs_p0q0, abs_p0q0);
        abs_p1q1 = _mm_srli_epi16(_mm_and_si128(abs_p1q1, fe), 1);
        mask = _mm_subs_epu8(_mm_adds_epu8(abs_p0q0, abs_p1q1), blimit);
        mask = _mm_xor_si128(_mm_cmpeq_epi8(mask, zero), ff);
        // mask |= (abs(p0 - q0) * 2 + abs(p1 - q1) / 2  > blimit) * -1;
        mask = _mm_max_epu8(abs_p1p0, mask);
        // mask |= (abs(p1 - p0) > limit) * -1;
        // mask |= (abs(q1 - q0) > limit) * -1;

        work = _mm_max_epu8(abs_diff(q2p2, q1p1), abs_diff(q3p3, q2p2));

        mask = _mm_max_epu8(work, mask);
        mask = _mm_max_epu8(mask, _mm_srli_si128(mask, 8));
        mask = _mm_subs_epu8(mask, limit);
        mask = _mm_cmpeq_epi8(mask, zero);
        // replicate for the further "merged variables" usage
        mask = _mm_unpacklo_epi64(mask, mask);

        // flat_mask4

        flat = _mm_max_epu8(abs_diff(q2p2, q0p0), abs_diff(q3p3, q0p0));
        flat = _mm_max_epu8(abs_p1p0, flat);

        flat = _mm_max_epu8(flat, _mm_srli_si128(flat, 8));
        flat = _mm_subs_epu8(flat, one);
        flat = _mm_cmpeq_epi8(flat, zero);
        flat = _mm_and_si128(flat, mask);
        // replicate for the further "merged variables" usage
        flat = _mm_unpacklo_epi64(flat, flat);
    }

    // filter8
    {
        const __m128i four = _mm_set1_epi16(4);

        __m128i workp_a, workp_b, workp_shft0, workp_shft1;
        p2 = _mm_unpacklo_epi8(*p2_8, zero);
        p1 = _mm_unpacklo_epi8(*p1_8, zero);
        p0 = _mm_unpacklo_epi8(*p0_8, zero);
        q0 = _mm_unpacklo_epi8(*q0_8, zero);
        q1 = _mm_unpacklo_epi8(*q1_8, zero);
        q2 = _mm_unpacklo_epi8(*q2_8, zero);
        p3 = _mm_unpacklo_epi8(*p3_8, zero);
        q3 = _mm_unpacklo_epi8(*q3_8, zero);

        // op2
        workp_a = _mm_add_epi16(_mm_add_epi16(p3, p3), _mm_add_epi16(p2, p1));
        workp_a = _mm_add_epi16(_mm_add_epi16(workp_a, four), p0);
        workp_b = _mm_add_epi16(_mm_add_epi16(q0, p2), p3);
        workp_shft0 = _mm_srli_epi16(_mm_add_epi16(workp_a, workp_b), 3);
        op2 = _mm_packus_epi16(workp_shft0, workp_shft0);

        // op1
        workp_b = _mm_add_epi16(_mm_add_epi16(q0, q1), p1);
        workp_shft0 = _mm_srli_epi16(_mm_add_epi16(workp_a, workp_b), 3);

        // op0
        workp_a = _mm_add_epi16(_mm_sub_epi16(workp_a, p3), q2);
        workp_b = _mm_add_epi16(_mm_sub_epi16(workp_b, p1), p0);
        workp_shft1 = _mm_srli_epi16(_mm_add_epi16(workp_a, workp_b), 3);

        flat_p1p0 = _mm_packus_epi16(workp_shft1, workp_shft0);

        // oq0
        workp_a = _mm_add_epi16(_mm_sub_epi16(workp_a, p3), q3);
        workp_b = _mm_add_epi16(_mm_sub_epi16(workp_b, p0), q0);
        workp_shft0 = _mm_srli_epi16(_mm_add_epi16(workp_a, workp_b), 3);

        // oq1
        workp_a = _mm_add_epi16(_mm_sub_epi16(workp_a, p2), q3);
        workp_b = _mm_add_epi16(_mm_sub_epi16(workp_b, q0), q1);
        workp_shft1 = _mm_srli_epi16(_mm_add_epi16(workp_a, workp_b), 3);

        flat_q0q1 = _mm_packus_epi16(workp_shft0, workp_shft1);

        // oq2
        workp_a = _mm_add_epi16(_mm_sub_epi16(workp_a, p1), q3);
        workp_b = _mm_add_epi16(_mm_sub_epi16(workp_b, q1), q2);
        workp_shft1 = _mm_srli_epi16(_mm_add_epi16(workp_a, workp_b), 3);
        oq2 = _mm_packus_epi16(workp_shft1, workp_shft1);
    }

    // lp filter - the same for 8 and 6 versions
    p1p0 = _mm_unpacklo_epi64(q0p0, q1p1);
    q1q0 = _mm_unpackhi_epi64(q0p0, q1p1);

    filter4_sse2(&p1p0, &q1q0, &hev, &mask, &qs1qs0, &ps1ps0);

    qs1qs0 = _mm_andnot_si128(flat, qs1qs0);
    q1q0 = _mm_and_si128(flat, flat_q0q1);
    *q1q0_out = _mm_or_si128(qs1qs0, q1q0);

    ps1ps0 = _mm_andnot_si128(flat, ps1ps0);
    p1p0 = _mm_and_si128(flat, flat_p1p0);
    *p1p0_out = _mm_or_si128(ps1ps0, p1p0);

    work_a = _mm_andnot_si128(flat, *q2_8);
    q2 = _mm_and_si128(flat, oq2);
    *q2_out = _mm_or_si128(work_a, q2);

    work_a = _mm_andnot_si128(flat, *p2_8);
    p2 = _mm_and_si128(flat, op2);
    *p2_out = _mm_or_si128(work_a, p2);
}

void aom_lpf_horizontal_8_sse2(uint8_t *s, int32_t p,
    const uint8_t *_blimit,
    const uint8_t *_limit,
    const uint8_t *_thresh) {
    __m128i p2_8, p1_8, p0_8, q0_8, q1_8, q2_8, p3_8, q3_8;
    __m128i q1q0, p1p0, p2, q2;
    p3_8 = _mm_cvtsi32_si128(*(int32_t *)(s - 4 * p));
    p2_8 = _mm_cvtsi32_si128(*(int32_t *)(s - 3 * p));
    p1_8 = _mm_cvtsi32_si128(*(int32_t *)(s - 2 * p));
    p0_8 = _mm_cvtsi32_si128(*(int32_t *)(s - 1 * p));
    q0_8 = _mm_cvtsi32_si128(*(int32_t *)(s - 0 * p));
    q1_8 = _mm_cvtsi32_si128(*(int32_t *)(s + 1 * p));
    q2_8 = _mm_cvtsi32_si128(*(int32_t *)(s + 2 * p));
    q3_8 = _mm_cvtsi32_si128(*(int32_t *)(s + 3 * p));

    lpf_internal_8_sse2(&p3_8, &q3_8, &p2_8, &q2_8, &p1_8, &q1_8, &p0_8, &q0_8,
        &q1q0, &p1p0, &p2, &q2, _blimit, _limit, _thresh);

    xx_storel_32(s - 1 * p, p1p0);
    xx_storel_32(s - 2 * p, _mm_srli_si128(p1p0, 8));
    xx_storel_32(s + 0 * p, q1q0);
    xx_storel_32(s + 1 * p, _mm_srli_si128(q1q0, 8));
    xx_storel_32(s - 3 * p, p2);
    xx_storel_32(s + 2 * p, q2);
}
void aom_lpf_horizontal_14_dual_sse2(uint8_t *s, int32_t p,
    const uint8_t *_blimit,
    const uint8_t *_limit,
    const uint8_t *_thresh) {
    aom_lpf_horizontal_14_sse2(s, p, _blimit, _limit, _thresh);
    aom_lpf_horizontal_14_sse2(s + 4, p, _blimit, _limit, _thresh);
}

void aom_lpf_horizontal_8_dual_sse2(uint8_t *s, int32_t p, const uint8_t *_blimit0,
    const uint8_t *_limit0,
    const uint8_t *_thresh0,
    const uint8_t *_blimit1,
    const uint8_t *_limit1,
    const uint8_t *_thresh1) {
    DECLARE_ALIGNED(16, uint8_t, flat_op2[16]);
    DECLARE_ALIGNED(16, uint8_t, flat_op1[16]);
    DECLARE_ALIGNED(16, uint8_t, flat_op0[16]);
    DECLARE_ALIGNED(16, uint8_t, flat_oq2[16]);
    DECLARE_ALIGNED(16, uint8_t, flat_oq1[16]);
    DECLARE_ALIGNED(16, uint8_t, flat_oq0[16]);
    const __m128i zero = _mm_setzero_si128();
    const __m128i blimit =
        _mm_unpacklo_epi64(_mm_load_si128((const __m128i *)_blimit0),
            _mm_load_si128((const __m128i *)_blimit1));
    const __m128i limit =
        _mm_unpacklo_epi64(_mm_load_si128((const __m128i *)_limit0),
            _mm_load_si128((const __m128i *)_limit1));
    const __m128i thresh =
        _mm_unpacklo_epi64(_mm_load_si128((const __m128i *)_thresh0),
            _mm_load_si128((const __m128i *)_thresh1));

    __m128i mask, hev, flat;
    __m128i p3, p2, p1, p0, q0, q1, q2, q3;

    p3 = _mm_loadu_si128((__m128i *)(s - 4 * p));
    p2 = _mm_loadu_si128((__m128i *)(s - 3 * p));
    p1 = _mm_loadu_si128((__m128i *)(s - 2 * p));
    p0 = _mm_loadu_si128((__m128i *)(s - 1 * p));
    q0 = _mm_loadu_si128((__m128i *)(s - 0 * p));
    q1 = _mm_loadu_si128((__m128i *)(s + 1 * p));
    q2 = _mm_loadu_si128((__m128i *)(s + 2 * p));
    q3 = _mm_loadu_si128((__m128i *)(s + 3 * p));
    {
        const __m128i abs_p1p0 =
            _mm_or_si128(_mm_subs_epu8(p1, p0), _mm_subs_epu8(p0, p1));
        const __m128i abs_q1q0 =
            _mm_or_si128(_mm_subs_epu8(q1, q0), _mm_subs_epu8(q0, q1));
        const __m128i one = _mm_set1_epi8(1);
        const __m128i fe = _mm_set1_epi8((char)0xfe);
        const __m128i ff = _mm_cmpeq_epi8(abs_p1p0, abs_p1p0);
        __m128i abs_p0q0 =
            _mm_or_si128(_mm_subs_epu8(p0, q0), _mm_subs_epu8(q0, p0));
        __m128i abs_p1q1 =
            _mm_or_si128(_mm_subs_epu8(p1, q1), _mm_subs_epu8(q1, p1));
        __m128i work;

        // filter_mask and hev_mask
        flat = _mm_max_epu8(abs_p1p0, abs_q1q0);
        hev = _mm_subs_epu8(flat, thresh);
        hev = _mm_xor_si128(_mm_cmpeq_epi8(hev, zero), ff);

        abs_p0q0 = _mm_adds_epu8(abs_p0q0, abs_p0q0);
        abs_p1q1 = _mm_srli_epi16(_mm_and_si128(abs_p1q1, fe), 1);
        mask = _mm_subs_epu8(_mm_adds_epu8(abs_p0q0, abs_p1q1), blimit);
        mask = _mm_xor_si128(_mm_cmpeq_epi8(mask, zero), ff);
        // mask |= (abs(p0 - q0) * 2 + abs(p1 - q1) / 2  > blimit) * -1;
        mask = _mm_max_epu8(flat, mask);
        // mask |= (abs(p1 - p0) > limit) * -1;
        // mask |= (abs(q1 - q0) > limit) * -1;
        work = _mm_max_epu8(
            _mm_or_si128(_mm_subs_epu8(p2, p1), _mm_subs_epu8(p1, p2)),
            _mm_or_si128(_mm_subs_epu8(p3, p2), _mm_subs_epu8(p2, p3)));
        mask = _mm_max_epu8(work, mask);
        work = _mm_max_epu8(
            _mm_or_si128(_mm_subs_epu8(q2, q1), _mm_subs_epu8(q1, q2)),
            _mm_or_si128(_mm_subs_epu8(q3, q2), _mm_subs_epu8(q2, q3)));
        mask = _mm_max_epu8(work, mask);
        mask = _mm_subs_epu8(mask, limit);
        mask = _mm_cmpeq_epi8(mask, zero);

        // flat_mask4
        work = _mm_max_epu8(
            _mm_or_si128(_mm_subs_epu8(p2, p0), _mm_subs_epu8(p0, p2)),
            _mm_or_si128(_mm_subs_epu8(q2, q0), _mm_subs_epu8(q0, q2)));
        flat = _mm_max_epu8(work, flat);
        work = _mm_max_epu8(
            _mm_or_si128(_mm_subs_epu8(p3, p0), _mm_subs_epu8(p0, p3)),
            _mm_or_si128(_mm_subs_epu8(q3, q0), _mm_subs_epu8(q0, q3)));
        flat = _mm_max_epu8(work, flat);
        flat = _mm_subs_epu8(flat, one);
        flat = _mm_cmpeq_epi8(flat, zero);
        flat = _mm_and_si128(flat, mask);
    }
    {
        const __m128i four = _mm_set1_epi16(4);
        uint8_t *src = s;
        int32_t i = 0;

        do {
            __m128i workp_a, workp_b, workp_shft;
            p3 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(src - 4 * p)), zero);
            p2 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(src - 3 * p)), zero);
            p1 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(src - 2 * p)), zero);
            p0 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(src - 1 * p)), zero);
            q0 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(src - 0 * p)), zero);
            q1 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(src + 1 * p)), zero);
            q2 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(src + 2 * p)), zero);
            q3 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(src + 3 * p)), zero);

            workp_a = _mm_add_epi16(_mm_add_epi16(p3, p3), _mm_add_epi16(p2, p1));
            workp_a = _mm_add_epi16(_mm_add_epi16(workp_a, four), p0);
            workp_b = _mm_add_epi16(_mm_add_epi16(q0, p2), p3);
            workp_shft = _mm_srli_epi16(_mm_add_epi16(workp_a, workp_b), 3);
            _mm_storel_epi64((__m128i *)&flat_op2[i * 8],
                _mm_packus_epi16(workp_shft, workp_shft));

            workp_b = _mm_add_epi16(_mm_add_epi16(q0, q1), p1);
            workp_shft = _mm_srli_epi16(_mm_add_epi16(workp_a, workp_b), 3);
            _mm_storel_epi64((__m128i *)&flat_op1[i * 8],
                _mm_packus_epi16(workp_shft, workp_shft));

            workp_a = _mm_add_epi16(_mm_sub_epi16(workp_a, p3), q2);
            workp_b = _mm_add_epi16(_mm_sub_epi16(workp_b, p1), p0);
            workp_shft = _mm_srli_epi16(_mm_add_epi16(workp_a, workp_b), 3);
            _mm_storel_epi64((__m128i *)&flat_op0[i * 8],
                _mm_packus_epi16(workp_shft, workp_shft));

            workp_a = _mm_add_epi16(_mm_sub_epi16(workp_a, p3), q3);
            workp_b = _mm_add_epi16(_mm_sub_epi16(workp_b, p0), q0);
            workp_shft = _mm_srli_epi16(_mm_add_epi16(workp_a, workp_b), 3);
            _mm_storel_epi64((__m128i *)&flat_oq0[i * 8],
                _mm_packus_epi16(workp_shft, workp_shft));

            workp_a = _mm_add_epi16(_mm_sub_epi16(workp_a, p2), q3);
            workp_b = _mm_add_epi16(_mm_sub_epi16(workp_b, q0), q1);
            workp_shft = _mm_srli_epi16(_mm_add_epi16(workp_a, workp_b), 3);
            _mm_storel_epi64((__m128i *)&flat_oq1[i * 8],
                _mm_packus_epi16(workp_shft, workp_shft));

            workp_a = _mm_add_epi16(_mm_sub_epi16(workp_a, p1), q3);
            workp_b = _mm_add_epi16(_mm_sub_epi16(workp_b, q1), q2);
            workp_shft = _mm_srli_epi16(_mm_add_epi16(workp_a, workp_b), 3);
            _mm_storel_epi64((__m128i *)&flat_oq2[i * 8],
                _mm_packus_epi16(workp_shft, workp_shft));

            src += 8;
        } while (++i < 2);
    }
    // lp filter
    {
        const __m128i t4 = _mm_set1_epi8(4);
        const __m128i t3 = _mm_set1_epi8(3);
        const __m128i t80 = _mm_set1_epi8((char)0x80);
        const __m128i te0 = _mm_set1_epi8((char)0xe0);
        const __m128i t1f = _mm_set1_epi8(0x1f);
        const __m128i t1 = _mm_set1_epi8(0x1);
        const __m128i t7f = _mm_set1_epi8(0x7f);

        const __m128i ps1 =
            _mm_xor_si128(_mm_loadu_si128((__m128i *)(s - 2 * p)), t80);
        const __m128i ps0 =
            _mm_xor_si128(_mm_loadu_si128((__m128i *)(s - 1 * p)), t80);
        const __m128i qs0 =
            _mm_xor_si128(_mm_loadu_si128((__m128i *)(s + 0 * p)), t80);
        const __m128i qs1 =
            _mm_xor_si128(_mm_loadu_si128((__m128i *)(s + 1 * p)), t80);
        __m128i filt;
        __m128i work_a;
        __m128i filter1, filter2;

        filt = _mm_and_si128(_mm_subs_epi8(ps1, qs1), hev);
        work_a = _mm_subs_epi8(qs0, ps0);
        filt = _mm_adds_epi8(filt, work_a);
        filt = _mm_adds_epi8(filt, work_a);
        filt = _mm_adds_epi8(filt, work_a);
        // (aom_filter + 3 * (qs0 - ps0)) & mask
        filt = _mm_and_si128(filt, mask);

        filter1 = _mm_adds_epi8(filt, t4);
        filter2 = _mm_adds_epi8(filt, t3);

        // Filter1 >> 3
        work_a = _mm_cmpgt_epi8(zero, filter1);
        filter1 = _mm_srli_epi16(filter1, 3);
        work_a = _mm_and_si128(work_a, te0);
        filter1 = _mm_and_si128(filter1, t1f);
        filter1 = _mm_or_si128(filter1, work_a);

        // Filter2 >> 3
        work_a = _mm_cmpgt_epi8(zero, filter2);
        filter2 = _mm_srli_epi16(filter2, 3);
        work_a = _mm_and_si128(work_a, te0);
        filter2 = _mm_and_si128(filter2, t1f);
        filter2 = _mm_or_si128(filter2, work_a);

        // filt >> 1
        filt = _mm_adds_epi8(filter1, t1);
        work_a = _mm_cmpgt_epi8(zero, filt);
        filt = _mm_srli_epi16(filt, 1);
        work_a = _mm_and_si128(work_a, t80);
        filt = _mm_and_si128(filt, t7f);
        filt = _mm_or_si128(filt, work_a);

        filt = _mm_andnot_si128(hev, filt);

        work_a = _mm_xor_si128(_mm_subs_epi8(qs0, filter1), t80);
        q0 = _mm_load_si128((__m128i *)flat_oq0);
        work_a = _mm_andnot_si128(flat, work_a);
        q0 = _mm_and_si128(flat, q0);
        q0 = _mm_or_si128(work_a, q0);

        work_a = _mm_xor_si128(_mm_subs_epi8(qs1, filt), t80);
        q1 = _mm_load_si128((__m128i *)flat_oq1);
        work_a = _mm_andnot_si128(flat, work_a);
        q1 = _mm_and_si128(flat, q1);
        q1 = _mm_or_si128(work_a, q1);

        work_a = _mm_loadu_si128((__m128i *)(s + 2 * p));
        q2 = _mm_load_si128((__m128i *)flat_oq2);
        work_a = _mm_andnot_si128(flat, work_a);
        q2 = _mm_and_si128(flat, q2);
        q2 = _mm_or_si128(work_a, q2);

        work_a = _mm_xor_si128(_mm_adds_epi8(ps0, filter2), t80);
        p0 = _mm_load_si128((__m128i *)flat_op0);
        work_a = _mm_andnot_si128(flat, work_a);
        p0 = _mm_and_si128(flat, p0);
        p0 = _mm_or_si128(work_a, p0);

        work_a = _mm_xor_si128(_mm_adds_epi8(ps1, filt), t80);
        p1 = _mm_load_si128((__m128i *)flat_op1);
        work_a = _mm_andnot_si128(flat, work_a);
        p1 = _mm_and_si128(flat, p1);
        p1 = _mm_or_si128(work_a, p1);

        work_a = _mm_loadu_si128((__m128i *)(s - 3 * p));
        p2 = _mm_load_si128((__m128i *)flat_op2);
        work_a = _mm_andnot_si128(flat, work_a);
        p2 = _mm_and_si128(flat, p2);
        p2 = _mm_or_si128(work_a, p2);

        _mm_storeu_si128((__m128i *)(s - 3 * p), p2);
        _mm_storeu_si128((__m128i *)(s - 2 * p), p1);
        _mm_storeu_si128((__m128i *)(s - 1 * p), p0);
        _mm_storeu_si128((__m128i *)(s + 0 * p), q0);
        _mm_storeu_si128((__m128i *)(s + 1 * p), q1);
        _mm_storeu_si128((__m128i *)(s + 2 * p), q2);
    }
}

void aom_lpf_horizontal_4_dual_sse2(uint8_t *s, int32_t p,
    const uint8_t *_blimit0,
    const uint8_t *_limit0,
    const uint8_t *_thresh0,
    const uint8_t *_blimit1,
    const uint8_t *_limit1,
    const uint8_t *_thresh1) {
    const __m128i blimit =
        _mm_unpacklo_epi64(_mm_load_si128((const __m128i *)_blimit0),
            _mm_load_si128((const __m128i *)_blimit1));
    const __m128i limit =
        _mm_unpacklo_epi64(_mm_load_si128((const __m128i *)_limit0),
            _mm_load_si128((const __m128i *)_limit1));
    const __m128i thresh =
        _mm_unpacklo_epi64(_mm_load_si128((const __m128i *)_thresh0),
            _mm_load_si128((const __m128i *)_thresh1));
    const __m128i zero = _mm_setzero_si128();
    __m128i p1, p0, q0, q1;
    __m128i mask, hev, flat;
    p1 = _mm_loadu_si128((__m128i *)(s - 2 * p));
    p0 = _mm_loadu_si128((__m128i *)(s - 1 * p));
    q0 = _mm_loadu_si128((__m128i *)(s - 0 * p));
    q1 = _mm_loadu_si128((__m128i *)(s + 1 * p));
    // filter_mask and hev_mask
    {
        const __m128i abs_p1p0 =
            _mm_or_si128(_mm_subs_epu8(p1, p0), _mm_subs_epu8(p0, p1));
        const __m128i abs_q1q0 =
            _mm_or_si128(_mm_subs_epu8(q1, q0), _mm_subs_epu8(q0, q1));
        const __m128i fe = _mm_set1_epi8((char)0xfe);
        const __m128i ff = _mm_cmpeq_epi8(abs_p1p0, abs_p1p0);
        __m128i abs_p0q0 =
            _mm_or_si128(_mm_subs_epu8(p0, q0), _mm_subs_epu8(q0, p0));
        __m128i abs_p1q1 =
            _mm_or_si128(_mm_subs_epu8(p1, q1), _mm_subs_epu8(q1, p1));
        flat = _mm_max_epu8(abs_p1p0, abs_q1q0);
        hev = _mm_subs_epu8(flat, thresh);
        hev = _mm_xor_si128(_mm_cmpeq_epi8(hev, zero), ff);

        abs_p0q0 = _mm_adds_epu8(abs_p0q0, abs_p0q0);
        abs_p1q1 = _mm_srli_epi16(_mm_and_si128(abs_p1q1, fe), 1);
        mask = _mm_subs_epu8(_mm_adds_epu8(abs_p0q0, abs_p1q1), blimit);
        mask = _mm_xor_si128(_mm_cmpeq_epi8(mask, zero), ff);
        // mask |= (abs(p0 - q0) * 2 + abs(p1 - q1) / 2  > blimit) * -1;
        mask = _mm_max_epu8(flat, mask);
        mask = _mm_subs_epu8(mask, limit);
        mask = _mm_cmpeq_epi8(mask, zero);
    }

    // filter4
    {
        const __m128i t4 = _mm_set1_epi8(4);
        const __m128i t3 = _mm_set1_epi8(3);
        const __m128i t80 = _mm_set1_epi8((char)0x80);
        const __m128i te0 = _mm_set1_epi8((char)0xe0);
        const __m128i t1f = _mm_set1_epi8(0x1f);
        const __m128i t1 = _mm_set1_epi8(0x1);
        const __m128i t7f = _mm_set1_epi8(0x7f);

        const __m128i ps1 =
            _mm_xor_si128(_mm_loadu_si128((__m128i *)(s - 2 * p)), t80);
        const __m128i ps0 =
            _mm_xor_si128(_mm_loadu_si128((__m128i *)(s - 1 * p)), t80);
        const __m128i qs0 =
            _mm_xor_si128(_mm_loadu_si128((__m128i *)(s + 0 * p)), t80);
        const __m128i qs1 =
            _mm_xor_si128(_mm_loadu_si128((__m128i *)(s + 1 * p)), t80);
        __m128i filt;
        __m128i work_a;
        __m128i filter1, filter2;

        filt = _mm_and_si128(_mm_subs_epi8(ps1, qs1), hev);
        work_a = _mm_subs_epi8(qs0, ps0);
        filt = _mm_adds_epi8(filt, work_a);
        filt = _mm_adds_epi8(filt, work_a);
        filt = _mm_adds_epi8(filt, work_a);
        // (aom_filter + 3 * (qs0 - ps0)) & mask
        filt = _mm_and_si128(filt, mask);

        filter1 = _mm_adds_epi8(filt, t4);
        filter2 = _mm_adds_epi8(filt, t3);

        // Filter1 >> 3
        work_a = _mm_cmpgt_epi8(zero, filter1);
        filter1 = _mm_srli_epi16(filter1, 3);
        work_a = _mm_and_si128(work_a, te0);
        filter1 = _mm_and_si128(filter1, t1f);
        filter1 = _mm_or_si128(filter1, work_a);

        // Filter2 >> 3
        work_a = _mm_cmpgt_epi8(zero, filter2);
        filter2 = _mm_srli_epi16(filter2, 3);
        work_a = _mm_and_si128(work_a, te0);
        filter2 = _mm_and_si128(filter2, t1f);
        filter2 = _mm_or_si128(filter2, work_a);

        // filt >> 1
        filt = _mm_adds_epi8(filter1, t1);
        work_a = _mm_cmpgt_epi8(zero, filt);
        filt = _mm_srli_epi16(filt, 1);
        work_a = _mm_and_si128(work_a, t80);
        filt = _mm_and_si128(filt, t7f);
        filt = _mm_or_si128(filt, work_a);

        filt = _mm_andnot_si128(hev, filt);

        q0 = _mm_xor_si128(_mm_subs_epi8(qs0, filter1), t80);
        q1 = _mm_xor_si128(_mm_subs_epi8(qs1, filt), t80);
        p0 = _mm_xor_si128(_mm_adds_epi8(ps0, filter2), t80);
        p1 = _mm_xor_si128(_mm_adds_epi8(ps1, filt), t80);

        _mm_storeu_si128((__m128i *)(s - 2 * p), p1);
        _mm_storeu_si128((__m128i *)(s - 1 * p), p0);
        _mm_storeu_si128((__m128i *)(s + 0 * p), q0);
        _mm_storeu_si128((__m128i *)(s + 1 * p), q1);
    }
}

static INLINE void transpose8x16(uint8_t *in0, uint8_t *in1,
    int32_t in_p, uint8_t *out, int32_t out_p) {
    __m128i x0, x1, x2, x3, x4, x5, x6, x7;
    __m128i x8, x9, x10, x11, x12, x13, x14, x15;

    // 2-way interleave w/hoisting of unpacks
    x0 = _mm_loadl_epi64((__m128i *)in0);           // 1
    x1 = _mm_loadl_epi64((__m128i *)(in0 + in_p));  // 3
    x0 = _mm_unpacklo_epi8(x0, x1);                 // 1

    x2 = _mm_loadl_epi64((__m128i *)(in0 + 2 * in_p));  // 5
    x3 = _mm_loadl_epi64((__m128i *)(in0 + 3 * in_p));  // 7
    x1 = _mm_unpacklo_epi8(x2, x3);                     // 2

    x4 = _mm_loadl_epi64((__m128i *)(in0 + 4 * in_p));  // 9
    x5 = _mm_loadl_epi64((__m128i *)(in0 + 5 * in_p));  // 11
    x2 = _mm_unpacklo_epi8(x4, x5);                     // 3

    x6 = _mm_loadl_epi64((__m128i *)(in0 + 6 * in_p));  // 13
    x7 = _mm_loadl_epi64((__m128i *)(in0 + 7 * in_p));  // 15
    x3 = _mm_unpacklo_epi8(x6, x7);                     // 4
    x4 = _mm_unpacklo_epi16(x0, x1);                    // 9

    x8 = _mm_loadl_epi64((__m128i *)in1);           // 2
    x9 = _mm_loadl_epi64((__m128i *)(in1 + in_p));  // 4
    x8 = _mm_unpacklo_epi8(x8, x9);                 // 5
    x5 = _mm_unpacklo_epi16(x2, x3);                // 10

    x10 = _mm_loadl_epi64((__m128i *)(in1 + 2 * in_p));  // 6
    x11 = _mm_loadl_epi64((__m128i *)(in1 + 3 * in_p));  // 8
    x9 = _mm_unpacklo_epi8(x10, x11);                    // 6

    x12 = _mm_loadl_epi64((__m128i *)(in1 + 4 * in_p));  // 10
    x13 = _mm_loadl_epi64((__m128i *)(in1 + 5 * in_p));  // 12
    x10 = _mm_unpacklo_epi8(x12, x13);                   // 7
    x12 = _mm_unpacklo_epi16(x8, x9);                    // 11

    x14 = _mm_loadl_epi64((__m128i *)(in1 + 6 * in_p));  // 14
    x15 = _mm_loadl_epi64((__m128i *)(in1 + 7 * in_p));  // 16
    x11 = _mm_unpacklo_epi8(x14, x15);                   // 8
    x13 = _mm_unpacklo_epi16(x10, x11);                  // 12

    x6 = _mm_unpacklo_epi32(x4, x5);     // 13
    x7 = _mm_unpackhi_epi32(x4, x5);     // 14
    x14 = _mm_unpacklo_epi32(x12, x13);  // 15
    x15 = _mm_unpackhi_epi32(x12, x13);  // 16

    // Store first 4-line result
    _mm_storeu_si128((__m128i *)out, _mm_unpacklo_epi64(x6, x14));
    _mm_storeu_si128((__m128i *)(out + out_p), _mm_unpackhi_epi64(x6, x14));
    _mm_storeu_si128((__m128i *)(out + 2 * out_p), _mm_unpacklo_epi64(x7, x15));
    _mm_storeu_si128((__m128i *)(out + 3 * out_p), _mm_unpackhi_epi64(x7, x15));

    x4 = _mm_unpackhi_epi16(x0, x1);
    x5 = _mm_unpackhi_epi16(x2, x3);
    x12 = _mm_unpackhi_epi16(x8, x9);
    x13 = _mm_unpackhi_epi16(x10, x11);

    x6 = _mm_unpacklo_epi32(x4, x5);
    x7 = _mm_unpackhi_epi32(x4, x5);
    x14 = _mm_unpacklo_epi32(x12, x13);
    x15 = _mm_unpackhi_epi32(x12, x13);

    // Store second 4-line result
    _mm_storeu_si128((__m128i *)(out + 4 * out_p), _mm_unpacklo_epi64(x6, x14));
    _mm_storeu_si128((__m128i *)(out + 5 * out_p), _mm_unpackhi_epi64(x6, x14));
    _mm_storeu_si128((__m128i *)(out + 6 * out_p), _mm_unpacklo_epi64(x7, x15));
    _mm_storeu_si128((__m128i *)(out + 7 * out_p), _mm_unpackhi_epi64(x7, x15));
}

#define movq(p) _mm_loadl_epi64((const __m128i *)(p))
#define punpcklbw(r0, r1) _mm_unpacklo_epi8(r0, r1)
#define punpcklwd(r0, r1) _mm_unpacklo_epi16(r0, r1)
#define punpckhwd(r0, r1) _mm_unpackhi_epi16(r0, r1)
#define pshufd(r, imm) _mm_shuffle_epi32(r, imm)
enum { ROTATE_DWORD_RIGHT = 0x39 };
static INLINE void transpose16x4(uint8_t *pDst, const ptrdiff_t dst_stride,
    const uint8_t *pSrc,
    const ptrdiff_t src_stride) {
    for (uint32_t idx = 0; idx < 2; idx += 1) {
        __m128i r0, r1, r2, r3;
        // load data
        r0 = movq(pSrc);
        r1 = movq(pSrc + src_stride);
        r2 = movq(pSrc + src_stride * 2);
        r3 = movq(pSrc + src_stride * 3);
        // transpose
        r0 = punpcklbw(r0, r1);
        r2 = punpcklbw(r2, r3);
        r1 = punpckhwd(r0, r2);
        r0 = punpcklwd(r0, r2);
        // store data
        xx_storel_32(pDst, r0);
        r0 = pshufd(r0, ROTATE_DWORD_RIGHT);
        xx_storel_32(pDst + dst_stride, r0);
        r0 = pshufd(r0, ROTATE_DWORD_RIGHT);
        xx_storel_32(pDst + dst_stride * 2, r0);
        r0 = pshufd(r0, ROTATE_DWORD_RIGHT);
        xx_storel_32(pDst + dst_stride * 3, r0);
        xx_storel_32(pDst + dst_stride * 4, r1);
        r1 = pshufd(r1, ROTATE_DWORD_RIGHT);
        xx_storel_32(pDst + dst_stride * 5, r1);
        r1 = pshufd(r1, ROTATE_DWORD_RIGHT);
        xx_storel_32(pDst + dst_stride * 6, r1);
        r1 = pshufd(r1, ROTATE_DWORD_RIGHT);
        xx_storel_32(pDst + dst_stride * 7, r1);
        // advance the pointers
        pDst += dst_stride * 8;
        pSrc += 8;
    }
}

static INLINE void transpose6x6_sse2(__m128i *x0, __m128i *x1, __m128i *x2,
    __m128i *x3, __m128i *x4, __m128i *x5,
    __m128i *d0d1, __m128i *d2d3,
    __m128i *d4d5) {
    __m128i w0, w1, w2, w4, w5;
    // x0  00 01 02 03 04 05 xx xx
    // x1  10 11 12 13 14 15 xx xx

    w0 = _mm_unpacklo_epi8(*x0, *x1);
    // 00 10 01 11 02 12 03 13  04 14 05 15 xx xx  xx xx

    // x2 20 21 22 23 24 25
    // x3 30 31 32 33 34 35

    w1 = _mm_unpacklo_epi8(
        *x2, *x3);  // 20 30 21 31 22 32 23 33  24 34 25 35 xx xx  xx xx

    // x4 40 41 42 43 44 45
    // x5 50 51 52 53 54 55

    w2 = _mm_unpacklo_epi8(*x4, *x5);  // 40 50 41 51 42 52 43 53 44 54 45 55

    w4 = _mm_unpacklo_epi16(
        w0, w1);  // 00 10 20 30 01 11 21 31 02 12 22 32 03 13 23 33
    w5 = _mm_unpacklo_epi16(
        w2, w0);  // 40 50 xx xx 41 51 xx xx 42 52 xx xx 43 53 xx xx

    *d0d1 = _mm_unpacklo_epi32(
        w4, w5);  // 00 10 20 30 40 50 xx xx 01 11 21 31 41 51 xx xx

    *d2d3 = _mm_unpackhi_epi32(
        w4, w5);  // 02 12 22 32 42 52 xx xx 03 13 23 33 43 53 xx xx

    w4 = _mm_unpackhi_epi16(
        w0, w1);  // 04 14 24 34 05 15 25 35 xx xx xx xx xx xx xx xx
    w5 = _mm_unpackhi_epi16(
        w2, *x3);  // 44 54 xx xx 45 55 xx xx xx xx xx xx xx xx xx xx
    *d4d5 = _mm_unpacklo_epi32(
        w4, w5);  // 04 14 24 34 44 54 xx xx 05 15 25 35 45 55 xx xx
}

static INLINE void transpose8x8_sse2(__m128i *x0, __m128i *x1, __m128i *x2,
    __m128i *x3, __m128i *x4, __m128i *x5,
    __m128i *x6, __m128i *x7, __m128i *d0d1,
    __m128i *d2d3, __m128i *d4d5,
    __m128i *d6d7) {
    __m128i w0, w1, w2, w3, w4, w5, w6, w7;
    // x0 00 01 02 03 04 05 06 07
    // x1 10 11 12 13 14 15 16 17
    w0 = _mm_unpacklo_epi8(
        *x0, *x1);  // 00 10 01 11 02 12 03 13 04 14 05 15 06 16 07 17

    // x2 20 21 22 23 24 25 26 27
    // x3 30 31 32 33 34 35 36 37
    w1 = _mm_unpacklo_epi8(
        *x2, *x3);  // 20 30 21 31 22 32 23 33 24 34 25 35 26 36 27 37

    // x4 40 41 42 43 44 45 46 47
    // x5  50 51 52 53 54 55 56 57
    w2 = _mm_unpacklo_epi8(
        *x4, *x5);  // 40 50 41 51 42 52 43 53 44 54 45 55 46 56 47 57

    // x6  60 61 62 63 64 65 66 67
    // x7 70 71 72 73 74 75 76 77
    w3 = _mm_unpacklo_epi8(
        *x6, *x7);  // 60 70 61 71 62 72 63 73 64 74 65 75 66 76 67 77

    w4 = _mm_unpacklo_epi16(
        w0, w1);  // 00 10 20 30 01 11 21 31 02 12 22 32 03 13 23 33
    w5 = _mm_unpacklo_epi16(
        w2, w3);  // 40 50 60 70 41 51 61 71 42 52 62 72 43 53 63 73

    *d0d1 = _mm_unpacklo_epi32(
        w4, w5);  // 00 10 20 30 40 50 60 70 01 11 21 31 41 51 61 71
    *d2d3 = _mm_unpackhi_epi32(
        w4, w5);  // 02 12 22 32 42 52 62 72 03 13 23 33 43 53 63 73

    w6 = _mm_unpackhi_epi16(
        w0, w1);  // 04 14 24 34 05 15 25 35 06 16 26 36 07 17 27 37
    w7 = _mm_unpackhi_epi16(
        w2, w3);  // 44 54 64 74 45 55 65 75 46 56 66 76 47 57 67 77

    *d4d5 = _mm_unpacklo_epi32(
        w6, w7);  // 04 14 24 34 44 54 64 74 05 15 25 35 45 55 65 75
    *d6d7 = _mm_unpackhi_epi32(
        w6, w7);  // 06 16 26 36 46 56 66 76 07 17 27 37 47 57 67 77
}

void aom_lpf_vertical_4_dual_sse2(uint8_t *s, int32_t p, const uint8_t *blimit0,
    const uint8_t *limit0, const uint8_t *thresh0,
    const uint8_t *blimit1, const uint8_t *limit1,
    const uint8_t *thresh1) {
    DECLARE_ALIGNED(16, uint8_t, t_dst[16 * 8]);
    // Transpose 8x16
    transpose8x16(s - 4, s - 4 + p * 8, p, t_dst, 16);

    // Loop filtering
    aom_lpf_horizontal_4_dual_sse2(t_dst + 4 * 16, 16, blimit0, limit0, thresh0,
        blimit1, limit1, thresh1);
    transpose16x4(s - 2, p, t_dst + 16 * 2, 16);
}

void aom_lpf_vertical_6_sse2(uint8_t *s, int32_t p,
    const uint8_t *blimit,
    const uint8_t *limit,
    const uint8_t *thresh) {
    __m128i d0d1, d2d3, d4d5;
    __m128i d1, d3, d5;

    __m128i p2, p1, p0, q0, q1, q2;
    __m128i p1p0, q1q0;
    DECLARE_ALIGNED(16, uint8_t, temp_dst[16]);

    p2 = _mm_loadl_epi64((__m128i *)((s - 3) + 0 * p));
    p1 = _mm_loadl_epi64((__m128i *)((s - 3) + 1 * p));
    p0 = _mm_loadl_epi64((__m128i *)((s - 3) + 2 * p));
    q0 = _mm_loadl_epi64((__m128i *)((s - 3) + 3 * p));
    q1 = _mm_setzero_si128();
    q2 = _mm_setzero_si128();

    transpose6x6_sse2(&p2, &p1, &p0, &q0, &q1, &q2, &d0d1, &d2d3, &d4d5);

    d1 = _mm_srli_si128(d0d1, 8);
    d3 = _mm_srli_si128(d2d3, 8);
    d5 = _mm_srli_si128(d4d5, 8);

    // Loop filtering
    lpf_internal_6_sse2(&d0d1, &d5, &d1, &d4d5, &d2d3, &d3, &q1q0, &p1p0, blimit,
        limit, thresh);

    p0 = _mm_srli_si128(p1p0, 8);
    q0 = _mm_srli_si128(q1q0, 8);

    transpose6x6_sse2(&d0d1, &p0, &p1p0, &q1q0, &q0, &d5, &d0d1, &d2d3, &d4d5);

    _mm_store_si128((__m128i *)(temp_dst), d0d1);
    memcpy((s - 3) + 0 * p, temp_dst, 6);
    memcpy((s - 3) + 1 * p, temp_dst + 8, 6);
    _mm_store_si128((__m128i *)(temp_dst), d2d3);
    memcpy((s - 3) + 2 * p, temp_dst, 6);
    memcpy((s - 3) + 3 * p, temp_dst + 8, 6);
}

void aom_lpf_vertical_8_sse2(uint8_t *s, int32_t p,
    const uint8_t *blimit,
    const uint8_t *limit,
    const uint8_t *thresh) {
    __m128i d0d1, d2d3, d4d5, d6d7;
    __m128i d1, d3, d5, d7;
    __m128i p2_8, p1_8, p0_8, q0_8, q1_8, q2_8, p3_8, q3_8;
    __m128i q1q0, p1p0, p0, q0, p2, q2;

    p3_8 = _mm_loadl_epi64((__m128i *)((s - 4) + 0 * p));
    p2_8 = _mm_loadl_epi64((__m128i *)((s - 4) + 1 * p));
    p1_8 = _mm_loadl_epi64((__m128i *)((s - 4) + 2 * p));
    p0_8 = _mm_loadl_epi64((__m128i *)((s - 4) + 3 * p));
    q0_8 = _mm_setzero_si128();
    q1_8 = _mm_setzero_si128();
    q2_8 = _mm_setzero_si128();
    q3_8 = _mm_setzero_si128();

    transpose8x8_sse2(&p3_8, &p2_8, &p1_8, &p0_8, &q0_8, &q1_8, &q2_8, &q3_8,
        &d0d1, &d2d3, &d4d5, &d6d7);

    d1 = _mm_srli_si128(d0d1, 8);
    d3 = _mm_srli_si128(d2d3, 8);
    d5 = _mm_srli_si128(d4d5, 8);
    d7 = _mm_srli_si128(d6d7, 8);

    // Loop filtering
    lpf_internal_8_sse2(&d0d1, &d7, &d1, &d6d7, &d2d3, &d5, &d3, &d4d5, &q1q0,
        &p1p0, &p2, &q2, blimit, limit, thresh);

    p0 = _mm_srli_si128(p1p0, 8);
    q0 = _mm_srli_si128(q1q0, 8);

    transpose8x8_sse2(&d0d1, &p2, &p0, &p1p0, &q1q0, &q0, &q2, &d7, &d0d1, &d2d3,
        &d4d5, &d6d7);

    _mm_storel_epi64((__m128i *)(s - 4 + 0 * p), d0d1);
    _mm_storel_epi64((__m128i *)(s - 4 + 1 * p), _mm_srli_si128(d0d1, 8));
    _mm_storel_epi64((__m128i *)(s - 4 + 2 * p), d2d3);
    _mm_storel_epi64((__m128i *)(s - 4 + 3 * p), _mm_srli_si128(d2d3, 8));
}

void aom_lpf_vertical_14_sse2(uint8_t *s, int32_t p,
    const uint8_t *blimit,
    const uint8_t *limit,
    const uint8_t *thresh) {
    __m128i q6p6, q5p5, q4p4, q3p3, q2p2, q1p1, q0p0;
    __m128i p6, p5, p4, p3, p2, p1, p0, q0;
    __m128i p6_2, p5_2, p4_2, p3_2, p2_2, p1_2, q0_2, p0_2;
    __m128i d1, d3, d5, d7, d1_2, d3_2, d5_2, d7_2;
    __m128i d0d1, d2d3, d4d5, d6d7;
    __m128i d0d1_2, d2d3_2, d4d5_2, d6d7_2;

    p6 = _mm_loadl_epi64((__m128i *)((s - 8) + 0 * p));
    p5 = _mm_loadl_epi64((__m128i *)((s - 8) + 1 * p));
    p4 = _mm_loadl_epi64((__m128i *)((s - 8) + 2 * p));
    p3 = _mm_loadl_epi64((__m128i *)((s - 8) + 3 * p));
    p2 = _mm_setzero_si128();
    p1 = _mm_setzero_si128();
    p0 = _mm_setzero_si128();
    q0 = _mm_setzero_si128();

    transpose8x8_sse2(&p6, &p5, &p4, &p3, &p2, &p1, &p0, &q0, &d0d1, &d2d3, &d4d5,
        &d6d7);

    p6_2 = _mm_loadl_epi64((__m128i *)(s + 0 * p));
    p5_2 = _mm_loadl_epi64((__m128i *)(s + 1 * p));
    p4_2 = _mm_loadl_epi64((__m128i *)(s + 2 * p));
    p3_2 = _mm_loadl_epi64((__m128i *)(s + 3 * p));
    p2_2 = _mm_setzero_si128();
    p1_2 = _mm_setzero_si128();
    p0_2 = _mm_setzero_si128();
    q0_2 = _mm_setzero_si128();

    transpose8x8_sse2(&p6_2, &p5_2, &p4_2, &p3_2, &p2_2, &p1_2, &p0_2, &q0_2,
        &d0d1_2, &d2d3_2, &d4d5_2, &d6d7_2);

    d1 = _mm_srli_si128(d0d1, 8);
    d3 = _mm_srli_si128(d2d3, 8);
    d5 = _mm_srli_si128(d4d5, 8);
    d7 = _mm_srli_si128(d6d7, 8);

    d1_2 = _mm_srli_si128(d0d1_2, 8);
    d3_2 = _mm_srli_si128(d2d3_2, 8);
    d5_2 = _mm_srli_si128(d4d5_2, 8);
    d7_2 = _mm_srli_si128(d6d7_2, 8);

    q6p6 = _mm_unpacklo_epi64(d1, d6d7_2);
    q5p5 = _mm_unpacklo_epi64(d2d3, d5_2);
    q4p4 = _mm_unpacklo_epi64(d3, d4d5_2);
    q3p3 = _mm_unpacklo_epi64(d4d5, d3_2);
    q2p2 = _mm_unpacklo_epi64(d5, d2d3_2);
    q1p1 = _mm_unpacklo_epi64(d6d7, d1_2);
    q0p0 = _mm_unpacklo_epi64(d7, d0d1_2);

    lpf_internal_14_sse2(&q6p6, &q5p5, &q4p4, &q3p3, &q2p2, &q1p1, &q0p0, blimit,
        limit, thresh);

    transpose8x8_sse2(&d0d1, &d1, &q5p5, &q4p4, &q3p3, &q2p2, &q1p1, &q0p0, &d0d1,
        &d2d3, &d4d5, &d6d7);

    p0 = _mm_srli_si128(q0p0, 8);
    p1 = _mm_srli_si128(q1p1, 8);
    p2 = _mm_srli_si128(q2p2, 8);
    p3 = _mm_srli_si128(q3p3, 8);
    p4 = _mm_srli_si128(q4p4, 8);
    p5 = _mm_srli_si128(q5p5, 8);
    p6 = _mm_srli_si128(q6p6, 8);

    transpose8x8_sse2(&p0, &p1, &p2, &p3, &p4, &p5, &p6, &d7_2, &d0d1_2, &d2d3_2,
        &d4d5_2, &d6d7_2);

    _mm_storel_epi64((__m128i *)(s - 8 + 0 * p), d0d1);
    _mm_storel_epi64((__m128i *)(s - 8 + 1 * p), _mm_srli_si128(d0d1, 8));
    _mm_storel_epi64((__m128i *)(s - 8 + 2 * p), d2d3);
    _mm_storel_epi64((__m128i *)(s - 8 + 3 * p), _mm_srli_si128(d2d3, 8));

    _mm_storel_epi64((__m128i *)(s + 0 * p), d0d1_2);
    _mm_storel_epi64((__m128i *)(s + 1 * p), _mm_srli_si128(d0d1_2, 8));
    _mm_storel_epi64((__m128i *)(s + 2 * p), d2d3_2);
    _mm_storel_epi64((__m128i *)(s + 3 * p), _mm_srli_si128(d2d3_2, 8));
}

void aom_lpf_vertical_14_dual_sse2(uint8_t *s, int32_t p,
    const uint8_t *blimit, const uint8_t *limit,
    const uint8_t *thresh) {
    DECLARE_ALIGNED(16, uint8_t, t_dst[256]);
    // Transpose 16x16
    transpose8x16(s - 8, s - 8 + 8 * p, p, t_dst, 16);
    transpose8x16(s, s + 8 * p, p, t_dst + 8 * 16, 16);
    // Loop filtering
    aom_lpf_horizontal_14_dual_sse2(t_dst + 8 * 16, 16, blimit, limit, thresh);
    // Transpose back
    transpose8x16(t_dst, t_dst + 8 * 16, 16, s - 8, p);
    transpose8x16(t_dst + 8, t_dst + 8 + 8 * 16, 16, s - 8 + 8 * p, p);
}
/*********************************************************************************************/
//highbd_loopfilter_sse2.c
static AOM_FORCE_INLINE void pixel_clamp(const __m128i *min, const __m128i *max,
    __m128i *pixel) {
    *pixel = _mm_min_epi16(*pixel, *max);
    *pixel = _mm_max_epi16(*pixel, *min);
}

static AOM_FORCE_INLINE __m128i abs_diff16(__m128i a, __m128i b) {
    return _mm_or_si128(_mm_subs_epu16(a, b), _mm_subs_epu16(b, a));
}

static INLINE void get_limit(const uint8_t *bl, const uint8_t *l,
    const uint8_t *t, int32_t bd, __m128i *blt,
    __m128i *lt, __m128i *thr, __m128i *t80_out) {
    const int32_t shift = bd - 8;
    const __m128i zero = _mm_setzero_si128();

    __m128i x = _mm_unpacklo_epi8(_mm_load_si128((const __m128i *)bl), zero);
    *blt = _mm_slli_epi16(x, shift);

    x = _mm_unpacklo_epi8(_mm_load_si128((const __m128i *)l), zero);
    *lt = _mm_slli_epi16(x, shift);

    x = _mm_unpacklo_epi8(_mm_load_si128((const __m128i *)t), zero);
    *thr = _mm_slli_epi16(x, shift);

    *t80_out = _mm_set1_epi16(1 << (bd - 1));
}

static INLINE void get_limit_dual(
    const uint8_t *_blimit0, const uint8_t *_limit0, const uint8_t *_thresh0,
    const uint8_t *_blimit1, const uint8_t *_limit1, const uint8_t *_thresh1,
    int32_t bd, __m128i *blt_out, __m128i *lt_out, __m128i *thr_out,
    __m128i *t80_out) {
    const int32_t shift = bd - 8;
    const __m128i zero = _mm_setzero_si128();

    __m128i x0 =
        _mm_unpacklo_epi8(_mm_load_si128((const __m128i *)_blimit0), zero);
    __m128i x1 =
        _mm_unpacklo_epi8(_mm_load_si128((const __m128i *)_blimit1), zero);
    x0 = _mm_unpacklo_epi64(x0, x1);
    *blt_out = _mm_slli_epi16(x0, shift);

    x0 = _mm_unpacklo_epi8(_mm_load_si128((const __m128i *)_limit0), zero);
    x1 = _mm_unpacklo_epi8(_mm_load_si128((const __m128i *)_limit1), zero);
    x0 = _mm_unpacklo_epi64(x0, x1);
    *lt_out = _mm_slli_epi16(x0, shift);

    x0 = _mm_unpacklo_epi8(_mm_load_si128((const __m128i *)_thresh0), zero);
    x1 = _mm_unpacklo_epi8(_mm_load_si128((const __m128i *)_thresh1), zero);
    x0 = _mm_unpacklo_epi64(x0, x1);
    *thr_out = _mm_slli_epi16(x0, shift);

    *t80_out = _mm_set1_epi16(1 << (bd - 1));
}

static INLINE void load_highbd_pixel(const uint16_t *s, int32_t size, int32_t pitch,
    __m128i *p, __m128i *q) {
    int32_t i;
    for (i = 0; i < size; i++) {
        p[i] = _mm_loadu_si128((__m128i *)(s - (i + 1) * pitch));
        q[i] = _mm_loadu_si128((__m128i *)(s + i * pitch));
    }
}

static INLINE void highbd_hev_mask(const __m128i *p0q0, const __m128i *p1q1,
    const __m128i *t, __m128i *abs_p1p0,
    __m128i *hev) {
    *abs_p1p0 = abs_diff16(*p1q1, *p0q0);
    __m128i abs_q1q0 = _mm_srli_si128(*abs_p1p0, 8);
    __m128i h = _mm_max_epi16(*abs_p1p0, abs_q1q0);
    h = _mm_subs_epu16(h, *t);

    const __m128i ffff = _mm_set1_epi16((short)0xFFFF);
    const __m128i zero = _mm_setzero_si128();
    *hev = _mm_xor_si128(_mm_cmpeq_epi16(h, zero), ffff);
    // replicate for the further "merged variables" usage
    *hev = _mm_unpacklo_epi64(*hev, *hev);
}

static INLINE void highbd_filter_mask(const __m128i *p, const __m128i *q,
    const __m128i *l, const __m128i *bl,
    __m128i *mask) {
    __m128i abs_p0q0 = abs_diff16(p[0], q[0]);
    __m128i abs_p1q1 = abs_diff16(p[1], q[1]);
    abs_p0q0 = _mm_adds_epu16(abs_p0q0, abs_p0q0);
    abs_p1q1 = _mm_srli_epi16(abs_p1q1, 1);

    const __m128i zero = _mm_setzero_si128();
    const __m128i one = _mm_set1_epi16(1);
    const __m128i ffff = _mm_set1_epi16((short)0xFFFF);
    __m128i max = _mm_subs_epu16(_mm_adds_epu16(abs_p0q0, abs_p1q1), *bl);
    max = _mm_xor_si128(_mm_cmpeq_epi16(max, zero), ffff);
    max = _mm_and_si128(max, _mm_adds_epu16(*l, one));

    int32_t i;
    for (i = 1; i < 4; ++i) {
        max = _mm_max_epi16(max, abs_diff16(p[i], p[i - 1]));
        max = _mm_max_epi16(max, abs_diff16(q[i], q[i - 1]));
    }
    max = _mm_subs_epu16(max, *l);
    *mask = _mm_cmpeq_epi16(max, zero);  // return ~mask
}

static INLINE void flat_mask_internal(const __m128i *th, const __m128i *p,
    const __m128i *q, int32_t bd, int32_t start,
    int32_t end, __m128i *flat) {
    int32_t i;
    __m128i max =
        _mm_max_epi16(abs_diff16(q[start], q[0]), abs_diff16(p[start], p[0]));

    for (i = start + 1; i < end; ++i) {
        max = _mm_max_epi16(max, abs_diff16(p[i], p[0]));
        max = _mm_max_epi16(max, abs_diff16(q[i], q[0]));
    }

    __m128i ft;
    if (bd == 8)
        ft = _mm_subs_epu16(max, *th);
    else if (bd == 10)
        ft = _mm_subs_epu16(max, _mm_slli_epi16(*th, 2));
    else  // bd == 12
        ft = _mm_subs_epu16(max, _mm_slli_epi16(*th, 4));

    const __m128i zero = _mm_setzero_si128();
    *flat = _mm_cmpeq_epi16(ft, zero);
}

// Note:
//  Access p[3-1], p[0], and q[3-1], q[0]
static INLINE void highbd_flat_mask4(const __m128i *th, const __m128i *p,
    const __m128i *q, __m128i *flat, int32_t bd) {
    // check the distance 1,2,3 against 0
    flat_mask_internal(th, p, q, bd, 1, 4, flat);
}

// Note:
//  access p[6-4], p[0], and q[6-4], q[0]
static INLINE void highbd_flat_mask4_13(const __m128i *th, const __m128i *p,
    const __m128i *q, __m128i *flat,
    int32_t bd) {
    flat_mask_internal(th, p, q, bd, 4, 7, flat);
}

static AOM_FORCE_INLINE void highbd_filter4_sse2(__m128i *p1p0, __m128i *q1q0,
    __m128i *hev, __m128i *mask,
    __m128i *qs1qs0,
    __m128i *ps1ps0, __m128i *t80,
    int32_t bd) {
    const __m128i zero = _mm_setzero_si128();
    const __m128i one = _mm_set1_epi16(1);
    const __m128i pmax =
        _mm_subs_epi16(_mm_subs_epi16(_mm_slli_epi16(one, bd), one), *t80);
    const __m128i pmin = _mm_subs_epi16(zero, *t80);

    const __m128i t3t4 = _mm_set_epi16(3, 3, 3, 3, 4, 4, 4, 4);
    __m128i ps1ps0_work, qs1qs0_work, work;
    __m128i filt, filter2filter1, filter2filt, filter1filt;

    ps1ps0_work = _mm_subs_epi16(*p1p0, *t80);
    qs1qs0_work = _mm_subs_epi16(*q1q0, *t80);

    work = _mm_subs_epi16(ps1ps0_work, qs1qs0_work);
    pixel_clamp(&pmin, &pmax, &work);
    filt = _mm_and_si128(_mm_srli_si128(work, 8), *hev);

    filt = _mm_subs_epi16(filt, work);
    filt = _mm_subs_epi16(filt, work);
    filt = _mm_subs_epi16(filt, work);
    // (aom_filter + 3 * (qs0 - ps0)) & mask
    pixel_clamp(&pmin, &pmax, &filt);
    filt = _mm_and_si128(filt, *mask);
    filt = _mm_unpacklo_epi64(filt, filt);

    filter2filter1 = _mm_adds_epi16(filt, t3t4); /* signed_short_clamp */
    pixel_clamp(&pmin, &pmax, &filter2filter1);
    filter2filter1 = _mm_srai_epi16(filter2filter1, 3); /* >> 3 */

    filt = _mm_unpacklo_epi64(filter2filter1, filter2filter1);

    // filt >> 1
    filt = _mm_adds_epi16(filt, one);
    filt = _mm_srai_epi16(filt, 1);
    filt = _mm_andnot_si128(*hev, filt);

    filter2filt = _mm_unpackhi_epi64(filter2filter1, filt);
    filter1filt = _mm_unpacklo_epi64(filter2filter1, filt);

    qs1qs0_work = _mm_subs_epi16(qs1qs0_work, filter1filt);
    ps1ps0_work = _mm_adds_epi16(ps1ps0_work, filter2filt);

    pixel_clamp(&pmin, &pmax, &qs1qs0_work);
    pixel_clamp(&pmin, &pmax, &ps1ps0_work);

    *qs1qs0 = _mm_adds_epi16(qs1qs0_work, *t80);
    *ps1ps0 = _mm_adds_epi16(ps1ps0_work, *t80);
}

static INLINE void highbd_filter4_dual_sse2(__m128i *p, __m128i *q, __m128i *ps,
    __m128i *qs, const __m128i *mask,
    const __m128i *th, int32_t bd,

    __m128i *t80) {
    __m128i ps0 = _mm_subs_epi16(p[0], *t80);
    __m128i ps1 = _mm_subs_epi16(p[1], *t80);
    __m128i qs0 = _mm_subs_epi16(q[0], *t80);
    __m128i qs1 = _mm_subs_epi16(q[1], *t80);
    const __m128i one = _mm_set1_epi16(1);
    const __m128i pmax =
        _mm_subs_epi16(_mm_subs_epi16(_mm_slli_epi16(one, bd), one), *t80);

    const __m128i zero = _mm_setzero_si128();
    const __m128i pmin = _mm_subs_epi16(zero, *t80);
    __m128i filter = _mm_subs_epi16(ps1, qs1);
    pixel_clamp(&pmin, &pmax, &filter);

    // hev_filter
    __m128i hev;
    const __m128i abs_p1p0 = abs_diff16(p[1], p[0]);
    const __m128i abs_q1q0 = abs_diff16(q[1], q[0]);
    __m128i h = _mm_max_epi16(abs_p1p0, abs_q1q0);
    h = _mm_subs_epu16(h, *th);
    const __m128i ffff = _mm_cmpeq_epi16(h, h);
    hev = _mm_xor_si128(_mm_cmpeq_epi16(h, zero), ffff);

    filter = _mm_and_si128(filter, hev);

    const __m128i x = _mm_subs_epi16(qs0, ps0);
    filter = _mm_adds_epi16(filter, x);
    filter = _mm_adds_epi16(filter, x);
    filter = _mm_adds_epi16(filter, x);
    pixel_clamp(&pmin, &pmax, &filter);
    filter = _mm_and_si128(filter, *mask);
    const __m128i t3 = _mm_set1_epi16(3);
    const __m128i t4 = _mm_set1_epi16(4);
    __m128i filter1 = _mm_adds_epi16(filter, t4);
    __m128i filter2 = _mm_adds_epi16(filter, t3);
    pixel_clamp(&pmin, &pmax, &filter1);
    pixel_clamp(&pmin, &pmax, &filter2);
    filter1 = _mm_srai_epi16(filter1, 3);
    filter2 = _mm_srai_epi16(filter2, 3);
    qs0 = _mm_subs_epi16(qs0, filter1);
    pixel_clamp(&pmin, &pmax, &qs0);
    ps0 = _mm_adds_epi16(ps0, filter2);
    pixel_clamp(&pmin, &pmax, &ps0);
    qs[0] = _mm_adds_epi16(qs0, *t80);
    ps[0] = _mm_adds_epi16(ps0, *t80);
    filter = _mm_adds_epi16(filter1, one);
    filter = _mm_srai_epi16(filter, 1);
    filter = _mm_andnot_si128(hev, filter);
    qs1 = _mm_subs_epi16(qs1, filter);
    pixel_clamp(&pmin, &pmax, &qs1);
    ps1 = _mm_adds_epi16(ps1, filter);
    pixel_clamp(&pmin, &pmax, &ps1);
    qs[1] = _mm_adds_epi16(qs1, *t80);
    ps[1] = _mm_adds_epi16(ps1, *t80);
}

static AOM_FORCE_INLINE void highbd_lpf_internal_14_sse2(
    __m128i *p, __m128i *q, __m128i *pq, const uint8_t *blt,
    const uint8_t *lt, const uint8_t *thr, int32_t bd) {
    int32_t i;
    __m128i blimit, limit, thresh;
    __m128i t80;
    get_limit(blt, lt, thr, bd, &blimit, &limit, &thresh, &t80);

    __m128i mask;
    highbd_filter_mask(p, q, &limit, &blimit, &mask);

    __m128i flat, flat2;
    const __m128i one = _mm_set1_epi16(1);
    highbd_flat_mask4(&one, p, q, &flat, bd);
    highbd_flat_mask4_13(&one, p, q, &flat2, bd);

    flat = _mm_and_si128(flat, mask);
    flat2 = _mm_and_si128(flat2, flat);

    // replicate for the further "merged variables" usage
    flat = _mm_unpacklo_epi64(flat, flat);
    flat2 = _mm_unpacklo_epi64(flat2, flat2);

    __m128i ps0ps1, qs0qs1, p1p0, q1q0;

    // filters - hev and filter4
    __m128i hevhev;
    __m128i abs_p1p0;
    for (i = 0; i < 6; i++) {
        pq[i] = _mm_unpacklo_epi64(p[i], q[i]);
    }

    highbd_hev_mask(&pq[0], &pq[1], &thresh, &abs_p1p0, &hevhev);

    p1p0 = _mm_unpacklo_epi64(p[0], p[1]);
    q1q0 = _mm_unpacklo_epi64(q[0], q[1]);
    highbd_filter4_sse2(&p1p0, &q1q0, &hevhev, &mask, &qs0qs1, &ps0ps1, &t80, bd);

    // flat and wide flat calculations
    __m128i flat_p[3], flat_q[3], flat_pq[3];
    __m128i flat2_p[6], flat2_q[6];
    __m128i flat2_pq[6];

    {
        const __m128i eight = _mm_set1_epi16(8);
        const __m128i four = _mm_set1_epi16(4);
        __m128i sum_p = _mm_add_epi16(p[5], _mm_add_epi16(p[4], p[3]));
        __m128i sum_q = _mm_add_epi16(q[5], _mm_add_epi16(q[4], q[3]));

        __m128i sum_lp = _mm_add_epi16(p[0], _mm_add_epi16(p[2], p[1]));
        sum_p = _mm_add_epi16(sum_p, sum_lp);

        __m128i sum_lq = _mm_add_epi16(q[0], _mm_add_epi16(q[2], q[1]));
        sum_q = _mm_add_epi16(sum_q, sum_lq);

        sum_p = _mm_add_epi16(eight, _mm_add_epi16(sum_p, sum_q));
        sum_lp = _mm_add_epi16(four, _mm_add_epi16(sum_lp, sum_lq));

        flat2_p[0] = _mm_add_epi16(sum_p, _mm_add_epi16(_mm_add_epi16(p[6], p[0]),
            _mm_add_epi16(p[1], q[0])));
        flat2_q[0] = _mm_add_epi16(sum_p, _mm_add_epi16(_mm_add_epi16(q[6], q[0]),
            _mm_add_epi16(p[0], q[1])));

        flat_p[0] = _mm_add_epi16(sum_lp, _mm_add_epi16(p[3], p[0]));
        flat_q[0] = _mm_add_epi16(sum_lp, _mm_add_epi16(q[3], q[0]));
        __m128i sum_p6 = _mm_add_epi16(p[6], p[6]);
        __m128i sum_q6 = _mm_add_epi16(q[6], q[6]);
        __m128i sum_p3 = _mm_add_epi16(p[3], p[3]);
        __m128i sum_q3 = _mm_add_epi16(q[3], q[3]);

        sum_q = _mm_sub_epi16(sum_p, p[5]);
        sum_p = _mm_sub_epi16(sum_p, q[5]);

        flat2_p[1] = _mm_add_epi16(
            sum_p,
            _mm_add_epi16(sum_p6, _mm_add_epi16(p[1], _mm_add_epi16(p[2], p[0]))));
        flat2_q[1] = _mm_add_epi16(
            sum_q,
            _mm_add_epi16(sum_q6, _mm_add_epi16(q[1], _mm_add_epi16(q[0], q[2]))));

        sum_lq = _mm_sub_epi16(sum_lp, p[2]);
        sum_lp = _mm_sub_epi16(sum_lp, q[2]);

        flat_p[1] = _mm_add_epi16(sum_lp, _mm_add_epi16(sum_p3, p[1]));
        flat_q[1] = _mm_add_epi16(sum_lq, _mm_add_epi16(sum_q3, q[1]));

        flat_pq[0] = _mm_srli_epi16(_mm_unpacklo_epi64(flat_p[0], flat_q[0]), 3);
        flat_pq[1] = _mm_srli_epi16(_mm_unpacklo_epi64(flat_p[1], flat_q[1]), 3);

        flat2_pq[0] = _mm_srli_epi16(_mm_unpacklo_epi64(flat2_p[0], flat2_q[0]), 4);
        flat2_pq[1] = _mm_srli_epi16(_mm_unpacklo_epi64(flat2_p[1], flat2_q[1]), 4);

        sum_p6 = _mm_add_epi16(sum_p6, p[6]);
        sum_q6 = _mm_add_epi16(sum_q6, q[6]);
        sum_p3 = _mm_add_epi16(sum_p3, p[3]);
        sum_q3 = _mm_add_epi16(sum_q3, q[3]);
        sum_p = _mm_sub_epi16(sum_p, q[4]);
        sum_q = _mm_sub_epi16(sum_q, p[4]);
        flat2_p[2] = _mm_add_epi16(
            sum_p,
            _mm_add_epi16(sum_p6, _mm_add_epi16(p[2], _mm_add_epi16(p[3], p[1]))));
        flat2_q[2] = _mm_add_epi16(
            sum_q,
            _mm_add_epi16(sum_q6, _mm_add_epi16(q[2], _mm_add_epi16(q[1], q[3]))));
        flat2_pq[2] = _mm_srli_epi16(_mm_unpacklo_epi64(flat2_p[2], flat2_q[2]), 4);

        sum_lp = _mm_sub_epi16(sum_lp, q[1]);
        sum_lq = _mm_sub_epi16(sum_lq, p[1]);
        flat_p[2] = _mm_add_epi16(sum_lp, _mm_add_epi16(sum_p3, p[2]));
        flat_q[2] = _mm_add_epi16(sum_lq, _mm_add_epi16(sum_q3, q[2]));
        flat_pq[2] = _mm_srli_epi16(_mm_unpacklo_epi64(flat_p[2], flat_q[2]), 3);

        sum_p6 = _mm_add_epi16(sum_p6, p[6]);
        sum_q6 = _mm_add_epi16(sum_q6, q[6]);
        sum_p = _mm_sub_epi16(sum_p, q[3]);
        sum_q = _mm_sub_epi16(sum_q, p[3]);
        flat2_p[3] = _mm_add_epi16(
            sum_p,
            _mm_add_epi16(sum_p6, _mm_add_epi16(p[3], _mm_add_epi16(p[4], p[2]))));
        flat2_q[3] = _mm_add_epi16(
            sum_q,
            _mm_add_epi16(sum_q6, _mm_add_epi16(q[3], _mm_add_epi16(q[2], q[4]))));
        flat2_pq[3] = _mm_srli_epi16(_mm_unpacklo_epi64(flat2_p[3], flat2_q[3]), 4);

        sum_p6 = _mm_add_epi16(sum_p6, p[6]);
        sum_q6 = _mm_add_epi16(sum_q6, q[6]);
        sum_p = _mm_sub_epi16(sum_p, q[2]);
        sum_q = _mm_sub_epi16(sum_q, p[2]);
        flat2_p[4] = _mm_add_epi16(
            sum_p,
            _mm_add_epi16(sum_p6, _mm_add_epi16(p[4], _mm_add_epi16(p[5], p[3]))));
        flat2_q[4] = _mm_add_epi16(
            sum_q,
            _mm_add_epi16(sum_q6, _mm_add_epi16(q[4], _mm_add_epi16(q[3], q[5]))));
        flat2_pq[4] = _mm_srli_epi16(_mm_unpacklo_epi64(flat2_p[4], flat2_q[4]), 4);

        sum_p6 = _mm_add_epi16(sum_p6, p[6]);
        sum_q6 = _mm_add_epi16(sum_q6, q[6]);
        sum_p = _mm_sub_epi16(sum_p, q[1]);
        sum_q = _mm_sub_epi16(sum_q, p[1]);
        flat2_p[5] = _mm_add_epi16(
            sum_p,
            _mm_add_epi16(sum_p6, _mm_add_epi16(p[5], _mm_add_epi16(p[6], p[4]))));
        flat2_q[5] = _mm_add_epi16(
            sum_q,
            _mm_add_epi16(sum_q6, _mm_add_epi16(q[5], _mm_add_epi16(q[4], q[6]))));

        flat2_pq[5] = _mm_srli_epi16(_mm_unpacklo_epi64(flat2_p[5], flat2_q[5]), 4);
    }

    // highbd_filter8
    pq[0] = _mm_unpacklo_epi64(ps0ps1, qs0qs1);
    pq[1] = _mm_unpackhi_epi64(ps0ps1, qs0qs1);

    for (i = 0; i < 3; i++) {
        pq[i] = _mm_andnot_si128(flat, pq[i]);
        flat_pq[i] = _mm_and_si128(flat, flat_pq[i]);
        pq[i] = _mm_or_si128(pq[i], flat_pq[i]);
    }

    // highbd_filter16
    for (i = 5; i >= 0; i--) {
        //  p[i] remains unchanged if !(flat2 && flat && mask)
        pq[i] = _mm_andnot_si128(flat2, pq[i]);
        flat2_pq[i] = _mm_and_si128(flat2, flat2_pq[i]);
        //  get values for when (flat2 && flat && mask)
        pq[i] = _mm_or_si128(pq[i], flat2_pq[i]);  // full list of pq values
    }
}

void aom_highbd_lpf_horizontal_14_sse2(uint16_t *s, int32_t pitch,
    const uint8_t *blt, const uint8_t *lt,
    const uint8_t *thr, int32_t bd) {
    __m128i p[7], q[7], pq[7];
    int32_t i;

    for (i = 0; i < 7; i++) {
        p[i] = _mm_loadl_epi64((__m128i *)(s - (i + 1) * pitch));
        q[i] = _mm_loadl_epi64((__m128i *)(s + i * pitch));
    }

    highbd_lpf_internal_14_sse2(p, q, pq, blt, lt, thr, bd);

    for (i = 0; i < 6; i++) {
        _mm_storel_epi64((__m128i *)(s - (i + 1) * pitch), pq[i]);
        _mm_storel_epi64((__m128i *)(s + i * pitch), _mm_srli_si128(pq[i], 8));
    }
}

static AOM_FORCE_INLINE void highbd_lpf_internal_14_dual_sse2(
    __m128i *p, __m128i *q, const uint8_t *blt, const uint8_t *lt,
    const uint8_t *thr, int32_t bd) {
    __m128i blimit, limit, thresh, t80;
    get_limit(blt, lt, thr, bd, &blimit, &limit, &thresh, &t80);
    __m128i mask;
    highbd_filter_mask(p, q, &limit, &blimit, &mask);
    __m128i flat, flat2;
    const __m128i one = _mm_set1_epi16(1);
    highbd_flat_mask4(&one, p, q, &flat, bd);
    highbd_flat_mask4_13(&one, p, q, &flat2, bd);
    flat = _mm_and_si128(flat, mask);
    flat2 = _mm_and_si128(flat2, flat);
    __m128i ps[2], qs[2];
    highbd_filter4_dual_sse2(p, q, ps, qs, &mask, &thresh, bd, &t80);
    // flat and wide flat calculations
    __m128i flat_p[3], flat_q[3];
    __m128i flat2_p[6], flat2_q[6];
    {
        const __m128i eight = _mm_set1_epi16(8);
        const __m128i four = _mm_set1_epi16(4);
        __m128i sum_p = _mm_add_epi16(p[5], _mm_add_epi16(p[4], p[3]));
        __m128i sum_q = _mm_add_epi16(q[5], _mm_add_epi16(q[4], q[3]));
        __m128i sum_lp = _mm_add_epi16(p[0], _mm_add_epi16(p[2], p[1]));
        sum_p = _mm_add_epi16(sum_p, sum_lp);
        __m128i sum_lq = _mm_add_epi16(q[0], _mm_add_epi16(q[2], q[1]));
        sum_q = _mm_add_epi16(sum_q, sum_lq);
        sum_p = _mm_add_epi16(eight, _mm_add_epi16(sum_p, sum_q));
        sum_lp = _mm_add_epi16(four, _mm_add_epi16(sum_lp, sum_lq));
        flat2_p[0] = _mm_srli_epi16(
            _mm_add_epi16(sum_p, _mm_add_epi16(_mm_add_epi16(p[6], p[0]),
                _mm_add_epi16(p[1], q[0]))),
            4);
        flat2_q[0] = _mm_srli_epi16(
            _mm_add_epi16(sum_p, _mm_add_epi16(_mm_add_epi16(q[6], q[0]),
                _mm_add_epi16(p[0], q[1]))),
            4);
        flat_p[0] =
            _mm_srli_epi16(_mm_add_epi16(sum_lp, _mm_add_epi16(p[3], p[0])), 3);
        flat_q[0] =
            _mm_srli_epi16(_mm_add_epi16(sum_lp, _mm_add_epi16(q[3], q[0])), 3);
        __m128i sum_p6 = _mm_add_epi16(p[6], p[6]);
        __m128i sum_q6 = _mm_add_epi16(q[6], q[6]);
        __m128i sum_p3 = _mm_add_epi16(p[3], p[3]);
        __m128i sum_q3 = _mm_add_epi16(q[3], q[3]);
        sum_q = _mm_sub_epi16(sum_p, p[5]);
        sum_p = _mm_sub_epi16(sum_p, q[5]);
        flat2_p[1] = _mm_srli_epi16(
            _mm_add_epi16(
                sum_p, _mm_add_epi16(
                    sum_p6, _mm_add_epi16(p[1], _mm_add_epi16(p[2], p[0])))),
            4);
        flat2_q[1] = _mm_srli_epi16(
            _mm_add_epi16(
                sum_q, _mm_add_epi16(
                    sum_q6, _mm_add_epi16(q[1], _mm_add_epi16(q[0], q[2])))),
            4);
        sum_lq = _mm_sub_epi16(sum_lp, p[2]);
        sum_lp = _mm_sub_epi16(sum_lp, q[2]);
        flat_p[1] =
            _mm_srli_epi16(_mm_add_epi16(sum_lp, _mm_add_epi16(sum_p3, p[1])), 3);
        flat_q[1] =
            _mm_srli_epi16(_mm_add_epi16(sum_lq, _mm_add_epi16(sum_q3, q[1])), 3);
        sum_p6 = _mm_add_epi16(sum_p6, p[6]);
        sum_q6 = _mm_add_epi16(sum_q6, q[6]);
        sum_p3 = _mm_add_epi16(sum_p3, p[3]);
        sum_q3 = _mm_add_epi16(sum_q3, q[3]);
        sum_p = _mm_sub_epi16(sum_p, q[4]);
        sum_q = _mm_sub_epi16(sum_q, p[4]);
        flat2_p[2] = _mm_srli_epi16(
            _mm_add_epi16(
                sum_p, _mm_add_epi16(
                    sum_p6, _mm_add_epi16(p[2], _mm_add_epi16(p[3], p[1])))),
            4);
        flat2_q[2] = _mm_srli_epi16(
            _mm_add_epi16(
                sum_q, _mm_add_epi16(
                    sum_q6, _mm_add_epi16(q[2], _mm_add_epi16(q[1], q[3])))),
            4);
        sum_lp = _mm_sub_epi16(sum_lp, q[1]);
        sum_lq = _mm_sub_epi16(sum_lq, p[1]);
        flat_p[2] =
            _mm_srli_epi16(_mm_add_epi16(sum_lp, _mm_add_epi16(sum_p3, p[2])), 3);
        flat_q[2] =
            _mm_srli_epi16(_mm_add_epi16(sum_lq, _mm_add_epi16(sum_q3, q[2])), 3);
        sum_p6 = _mm_add_epi16(sum_p6, p[6]);
        sum_q6 = _mm_add_epi16(sum_q6, q[6]);
        sum_p = _mm_sub_epi16(sum_p, q[3]);
        sum_q = _mm_sub_epi16(sum_q, p[3]);
        flat2_p[3] = _mm_srli_epi16(
            _mm_add_epi16(
                sum_p, _mm_add_epi16(
                    sum_p6, _mm_add_epi16(p[3], _mm_add_epi16(p[4], p[2])))),
            4);
        flat2_q[3] = _mm_srli_epi16(
            _mm_add_epi16(
                sum_q, _mm_add_epi16(
                    sum_q6, _mm_add_epi16(q[3], _mm_add_epi16(q[2], q[4])))),
            4);
        sum_p6 = _mm_add_epi16(sum_p6, p[6]);
        sum_q6 = _mm_add_epi16(sum_q6, q[6]);
        sum_p = _mm_sub_epi16(sum_p, q[2]);
        sum_q = _mm_sub_epi16(sum_q, p[2]);
        flat2_p[4] = _mm_srli_epi16(
            _mm_add_epi16(
                sum_p, _mm_add_epi16(
                    sum_p6, _mm_add_epi16(p[4], _mm_add_epi16(p[5], p[3])))),
            4);
        flat2_q[4] = _mm_srli_epi16(
            _mm_add_epi16(
                sum_q, _mm_add_epi16(
                    sum_q6, _mm_add_epi16(q[4], _mm_add_epi16(q[3], q[5])))),
            4);
        sum_p6 = _mm_add_epi16(sum_p6, p[6]);
        sum_q6 = _mm_add_epi16(sum_q6, q[6]);
        sum_p = _mm_sub_epi16(sum_p, q[1]);
        sum_q = _mm_sub_epi16(sum_q, p[1]);
        flat2_p[5] = _mm_srli_epi16(
            _mm_add_epi16(
                sum_p, _mm_add_epi16(
                    sum_p6, _mm_add_epi16(p[5], _mm_add_epi16(p[6], p[4])))),
            4);
        flat2_q[5] = _mm_srli_epi16(
            _mm_add_epi16(
                sum_q, _mm_add_epi16(
                    sum_q6, _mm_add_epi16(q[5], _mm_add_epi16(q[4], q[6])))),
            4);
    }
    // highbd_filter8
    p[2] = _mm_andnot_si128(flat, p[2]);
    //  p2 remains unchanged if !(flat && mask)
    flat_p[2] = _mm_and_si128(flat, flat_p[2]);
    //  when (flat && mask)
    p[2] = _mm_or_si128(p[2], flat_p[2]);  // full list of p2 values
    q[2] = _mm_andnot_si128(flat, q[2]);
    flat_q[2] = _mm_and_si128(flat, flat_q[2]);
    q[2] = _mm_or_si128(q[2], flat_q[2]);  // full list of q2 values
    int32_t i;
    for (i = 1; i >= 0; i--) {
        ps[i] = _mm_andnot_si128(flat, ps[i]);
        flat_p[i] = _mm_and_si128(flat, flat_p[i]);
        p[i] = _mm_or_si128(ps[i], flat_p[i]);
        qs[i] = _mm_andnot_si128(flat, qs[i]);
        flat_q[i] = _mm_and_si128(flat, flat_q[i]);
        q[i] = _mm_or_si128(qs[i], flat_q[i]);
    }
    // highbd_filter16
    for (i = 5; i >= 0; i--) {
        //  p[i] remains unchanged if !(flat2 && flat && mask)
        p[i] = _mm_andnot_si128(flat2, p[i]);
        flat2_p[i] = _mm_and_si128(flat2, flat2_p[i]);
        //  get values for when (flat2 && flat && mask)
        p[i] = _mm_or_si128(p[i], flat2_p[i]);  // full list of p values
        q[i] = _mm_andnot_si128(flat2, q[i]);
        flat2_q[i] = _mm_and_si128(flat2, flat2_q[i]);
        q[i] = _mm_or_si128(q[i], flat2_q[i]);
    }
}

void aom_highbd_lpf_horizontal_14_dual_sse2(uint16_t *s, int32_t pitch,
    const uint8_t *_blimit,
    const uint8_t *_limit,
    const uint8_t *_thresh, int32_t bd) {
    __m128i p[7], q[7];
    int32_t i;
    load_highbd_pixel(s, 7, pitch, p, q);

    highbd_lpf_internal_14_dual_sse2(p, q, _blimit, _limit, _thresh, bd);
    for (i = 0; i < 6; i++) {
        _mm_store_si128((__m128i *)(s - (i + 1) * pitch), p[i]);
        _mm_store_si128((__m128i *)(s + i * pitch), q[i]);
    }
}

static AOM_FORCE_INLINE void highbd_lpf_internal_6_sse2(
    __m128i *p2, __m128i *p1, __m128i *p0, __m128i *q0, __m128i *q1,
    __m128i *q2, __m128i *p1p0_out, __m128i *q1q0_out, const uint8_t *_blimit,
    const uint8_t *_limit, const uint8_t *_thresh, int32_t bd) {
    const __m128i zero = _mm_setzero_si128();
    __m128i blimit, limit, thresh;
    __m128i mask, hev, flat;
    __m128i q2p2, q1p1, q0p0, p1q1, p0q0;
    __m128i p1p0, q1q0, ps1ps0, qs1qs0;
    __m128i flat_p1p0, flat_q0q1;

    q2p2 = _mm_unpacklo_epi64(*p2, *q2);
    q1p1 = _mm_unpacklo_epi64(*p1, *q1);
    q0p0 = _mm_unpacklo_epi64(*p0, *q0);

    p1q1 = _mm_shuffle_epi32(q1p1, _MM_SHUFFLE(1, 0, 3, 2));
    p0q0 = _mm_shuffle_epi32(q0p0, _MM_SHUFFLE(1, 0, 3, 2));

    __m128i abs_p1q1, abs_p0q0, abs_p1p0, work;

    const __m128i four = _mm_set1_epi16(4);
    __m128i t80;
    const __m128i one = _mm_set1_epi16(0x1);
    const __m128i ffff = _mm_cmpeq_epi16(one, one);

    get_limit(_blimit, _limit, _thresh, bd, &blimit, &limit, &thresh, &t80);

    // filter_mask and hev_mask
    highbd_hev_mask(&p0q0, &p1q1, &thresh, &abs_p1p0, &hev);

    abs_p0q0 = abs_diff16(q0p0, p0q0);
    abs_p1q1 = abs_diff16(q1p1, p1q1);

    abs_p0q0 = _mm_adds_epu16(abs_p0q0, abs_p0q0);
    abs_p1q1 = _mm_srli_epi16(abs_p1q1, 1);
    mask = _mm_subs_epu16(_mm_adds_epu16(abs_p0q0, abs_p1q1), blimit);
    mask = _mm_xor_si128(_mm_cmpeq_epi16(mask, zero), ffff);
    // mask |= (abs(*p0 - *q0) * 2 + abs(*p1 - *q1) / 2  > blimit) * -1;
    // So taking maximums continues to work:
    mask = _mm_and_si128(mask, _mm_adds_epu16(limit, one));
    mask = _mm_max_epi16(abs_p1p0, mask);
    // mask |= (abs(*p1 - *p0) > limit) * -1;
    // mask |= (abs(*q1 - *q0) > limit) * -1;

    work = abs_diff16(q2p2, q1p1);

    mask = _mm_max_epi16(work, mask);
    mask = _mm_max_epi16(mask, _mm_srli_si128(mask, 8));
    mask = _mm_subs_epu16(mask, limit);
    mask = _mm_cmpeq_epi16(mask, zero);

    // flat_mask
    flat = _mm_max_epi16(abs_diff16(q2p2, q0p0), abs_p1p0);
    flat = _mm_max_epi16(flat, _mm_srli_si128(flat, 8));

    if (bd == 8)
        flat = _mm_subs_epu16(flat, one);
    else if (bd == 10)
        flat = _mm_subs_epu16(flat, _mm_slli_epi16(one, 2));
    else  // bd == 12
        flat = _mm_subs_epu16(flat, _mm_slli_epi16(one, 4));

    flat = _mm_cmpeq_epi16(flat, zero);
    flat = _mm_and_si128(
        flat, mask);  // flat & mask
    // replicate for the further "merged variables" usage
    flat = _mm_unpacklo_epi64(flat, flat);

    {
        __m128i workp_a, workp_b, workp_shft0, workp_shft1;

        // op1
        workp_a = _mm_add_epi16(_mm_add_epi16(*p0, *p0),
            _mm_add_epi16(*p1, *p1));  // *p0 *2 + *p1 * 2
        workp_a = _mm_add_epi16(_mm_add_epi16(workp_a, four),
            *p2);  // *p2 + *p0 * 2 + *p1 * 2 + 4

        workp_b = _mm_add_epi16(_mm_add_epi16(*p2, *p2), *q0);
        workp_shft0 = _mm_add_epi16(
            workp_a, workp_b);  // *p2 * 3 + *p1 * 2 + *p0 * 2 + *q0 + 4

        // op0
        workp_b = _mm_add_epi16(_mm_add_epi16(*q0, *q0), *q1);  // *q0 * 2 + *q1
        workp_a =
            _mm_add_epi16(workp_a,
                workp_b);  // *p2 + *p0 * 2 + *p1 * 2 + *q0 * 2 + *q1 + 4

        flat_p1p0 = _mm_srli_epi16(_mm_unpacklo_epi64(workp_a, workp_shft0), 3);

        // oq0
        workp_a = _mm_sub_epi16(_mm_sub_epi16(workp_a, *p2),
            *p1);  // *p0 * 2 + *p1  + *q0 * 2 + *q1 + 4
        workp_b = _mm_add_epi16(*q1, *q2);
        workp_shft0 = _mm_add_epi16(
            workp_a, workp_b);  // *p0 * 2 + *p1  + *q0 * 2 + *q1 * 2 + *q2 + 4

        // oq1
        workp_a = _mm_sub_epi16(_mm_sub_epi16(workp_shft0, *p1),
            *p0);  // *p0   + *q0 * 2 + *q1 * 2 + *q2 + 4
        workp_b = _mm_add_epi16(*q2, *q2);
        workp_shft1 = _mm_add_epi16(
            workp_a, workp_b);  // *p0  + *q0 * 2 + *q1 * 2 + *q2 * 3 + 4

        flat_q0q1 = _mm_srli_epi16(_mm_unpacklo_epi64(workp_shft0, workp_shft1), 3);
    }
    // lp filter
    {
        p1p0 = _mm_unpacklo_epi64(q0p0, q1p1);
        q1q0 = _mm_unpackhi_epi64(q0p0, q1p1);

        highbd_filter4_sse2(&p1p0, &q1q0, &hev, &mask, &qs1qs0, &ps1ps0, &t80, bd);
    }

    qs1qs0 = _mm_andnot_si128(flat, qs1qs0);
    q1q0 = _mm_and_si128(flat, flat_q0q1);
    *q1q0_out = _mm_or_si128(qs1qs0, q1q0);

    ps1ps0 = _mm_andnot_si128(flat, ps1ps0);
    p1p0 = _mm_and_si128(flat, flat_p1p0);
    *p1p0_out = _mm_or_si128(ps1ps0, p1p0);
}

void aom_highbd_lpf_horizontal_6_sse2(uint16_t *s, int32_t p,
    const uint8_t *_blimit,
    const uint8_t *_limit,
    const uint8_t *_thresh, int32_t bd) {
    __m128i p2, p1, p0, q0, q1, q2, p1p0_out, q1q0_out;

    p2 = _mm_loadl_epi64((__m128i *)(s - 3 * p));
    p1 = _mm_loadl_epi64((__m128i *)(s - 2 * p));
    p0 = _mm_loadl_epi64((__m128i *)(s - 1 * p));
    q0 = _mm_loadl_epi64((__m128i *)(s + 0 * p));
    q1 = _mm_loadl_epi64((__m128i *)(s + 1 * p));
    q2 = _mm_loadl_epi64((__m128i *)(s + 2 * p));

    highbd_lpf_internal_6_sse2(&p2, &p1, &p0, &q0, &q1, &q2, &p1p0_out, &q1q0_out,
        _blimit, _limit, _thresh, bd);

    _mm_storel_epi64((__m128i *)(s - 2 * p), _mm_srli_si128(p1p0_out, 8));
    _mm_storel_epi64((__m128i *)(s - 1 * p), p1p0_out);
    _mm_storel_epi64((__m128i *)(s + 0 * p), q1q0_out);
    _mm_storel_epi64((__m128i *)(s + 1 * p), _mm_srli_si128(q1q0_out, 8));
}

static AOM_FORCE_INLINE void highbd_lpf_internal_8_sse2(
    __m128i *p3, __m128i *q3, __m128i *p2, __m128i *q2, __m128i *p1,
    __m128i *q1, __m128i *p0, __m128i *q0, __m128i *q1q0_out, __m128i *p1p0_out,
    const uint8_t *_blimit, const uint8_t *_limit,
    const uint8_t *_thresh, int32_t bd) {
    const __m128i zero = _mm_setzero_si128();
    __m128i blimit, limit, thresh;
    __m128i mask, hev, flat;
    __m128i q2p2, q1p1, q0p0, p1q1, p0q0, q3p3;
    __m128i p1p0, q1q0, ps1ps0, qs1qs0;
    __m128i work_a, op2, oq2, flat_p1p0, flat_q0q1;

    q3p3 = _mm_unpacklo_epi64(*p3, *q3);
    q2p2 = _mm_unpacklo_epi64(*p2, *q2);
    q1p1 = _mm_unpacklo_epi64(*p1, *q1);
    q0p0 = _mm_unpacklo_epi64(*p0, *q0);

    p1q1 = _mm_shuffle_epi32(q1p1, _MM_SHUFFLE(1, 0, 3, 2));
    p0q0 = _mm_shuffle_epi32(q0p0, _MM_SHUFFLE(1, 0, 3, 2));

    __m128i abs_p1q1, abs_p0q0, abs_p1p0, work;

    const __m128i four = _mm_set1_epi16(4);
    __m128i t80;
    const __m128i one = _mm_set1_epi16(0x1);
    const __m128i ffff = _mm_cmpeq_epi16(one, one);

    get_limit(_blimit, _limit, _thresh, bd, &blimit, &limit, &thresh, &t80);

    // filter_mask and hev_mask
    highbd_hev_mask(&p0q0, &p1q1, &thresh, &abs_p1p0, &hev);

    abs_p0q0 = abs_diff16(q0p0, p0q0);
    abs_p1q1 = abs_diff16(q1p1, p1q1);

    abs_p0q0 = _mm_adds_epu16(abs_p0q0, abs_p0q0);
    abs_p1q1 = _mm_srli_epi16(abs_p1q1, 1);
    mask = _mm_subs_epu16(_mm_adds_epu16(abs_p0q0, abs_p1q1), blimit);
    mask = _mm_xor_si128(_mm_cmpeq_epi16(mask, zero), ffff);
    // mask |= (abs(*p0 - q0) * 2 + abs(*p1 - q1) / 2  > blimit) * -1;
    // So taking maximums continues to work:
    mask = _mm_and_si128(mask, _mm_adds_epu16(limit, one));
    mask = _mm_max_epi16(abs_p1p0, mask);
    // mask |= (abs(*p1 - *p0) > limit) * -1;
    // mask |= (abs(q1 - q0) > limit) * -1;

    work = _mm_max_epi16(abs_diff16(q2p2, q1p1), abs_diff16(q3p3, q2p2));
    mask = _mm_max_epi16(work, mask);
    mask = _mm_max_epi16(mask, _mm_srli_si128(mask, 8));
    mask = _mm_subs_epu16(mask, limit);
    mask = _mm_cmpeq_epi16(mask, zero);

    // flat_mask4
    flat = _mm_max_epi16(abs_diff16(q2p2, q0p0), abs_diff16(q3p3, q0p0));
    flat = _mm_max_epi16(abs_p1p0, flat);
    flat = _mm_max_epi16(flat, _mm_srli_si128(flat, 8));

    if (bd == 8)
        flat = _mm_subs_epu16(flat, one);
    else if (bd == 10)
        flat = _mm_subs_epu16(flat, _mm_slli_epi16(one, 2));
    else  // bd == 12
        flat = _mm_subs_epu16(flat, _mm_slli_epi16(one, 4));

    flat = _mm_cmpeq_epi16(flat, zero);
    flat = _mm_and_si128(
        flat, mask);  // flat & mask
    // replicate for the further "merged variables" usage
    flat = _mm_unpacklo_epi64(flat, flat);

    {
        __m128i workp_a, workp_b, workp_shft0, workp_shft1;
        // Added before shift for rounding part of ROUND_POWER_OF_TWO

        // o*p2
        workp_a = _mm_add_epi16(_mm_add_epi16(*p3, *p3), _mm_add_epi16(*p2, *p1));
        workp_a = _mm_add_epi16(_mm_add_epi16(workp_a, four), *p0);
        workp_b = _mm_add_epi16(_mm_add_epi16(*q0, *p2), *p3);
        op2 = _mm_srli_epi16(_mm_add_epi16(workp_a, workp_b), 3);

        // o*p1
        workp_b = _mm_add_epi16(_mm_add_epi16(*q0, *q1), *p1);
        workp_shft0 = _mm_add_epi16(workp_a, workp_b);

        // o*p0
        workp_a = _mm_add_epi16(_mm_sub_epi16(workp_a, *p3), *q2);
        workp_b = _mm_add_epi16(_mm_sub_epi16(workp_b, *p1), *p0);
        workp_shft1 = _mm_add_epi16(workp_a, workp_b);

        flat_p1p0 = _mm_srli_epi16(_mm_unpacklo_epi64(workp_shft1, workp_shft0), 3);

        // oq0
        workp_a = _mm_add_epi16(_mm_sub_epi16(workp_a, *p3), *q3);
        workp_b = _mm_add_epi16(_mm_sub_epi16(workp_b, *p0), *q0);
        workp_shft0 = _mm_add_epi16(workp_a, workp_b);

        // oq1
        workp_a = _mm_add_epi16(_mm_sub_epi16(workp_a, *p2), *q3);
        workp_b = _mm_add_epi16(_mm_sub_epi16(workp_b, *q0), *q1);
        workp_shft1 = _mm_add_epi16(workp_a, workp_b);

        flat_q0q1 = _mm_srli_epi16(_mm_unpacklo_epi64(workp_shft0, workp_shft1), 3);

        // oq2
        workp_a = _mm_add_epi16(_mm_sub_epi16(workp_a, *p1), *q3);
        workp_b = _mm_add_epi16(_mm_sub_epi16(workp_b, *q1), *q2);
        oq2 = _mm_srli_epi16(_mm_add_epi16(workp_a, workp_b), 3);
    }

    // lp filter
    {
        p1p0 = _mm_unpacklo_epi64(q0p0, q1p1);
        q1q0 = _mm_unpackhi_epi64(q0p0, q1p1);

        highbd_filter4_sse2(&p1p0, &q1q0, &hev, &mask, &qs1qs0, &ps1ps0, &t80, bd);
    }

    qs1qs0 = _mm_andnot_si128(flat, qs1qs0);
    q1q0 = _mm_and_si128(flat, flat_q0q1);
    *q1q0_out = _mm_or_si128(qs1qs0, q1q0);

    ps1ps0 = _mm_andnot_si128(flat, ps1ps0);
    p1p0 = _mm_and_si128(flat, flat_p1p0);
    *p1p0_out = _mm_or_si128(ps1ps0, p1p0);

    work_a = _mm_andnot_si128(flat, *q2);
    *q2 = _mm_and_si128(flat, oq2);
    *q2 = _mm_or_si128(work_a, *q2);

    work_a = _mm_andnot_si128(flat, *p2);
    *p2 = _mm_and_si128(flat, op2);
    *p2 = _mm_or_si128(work_a, *p2);
}

static AOM_FORCE_INLINE void highbd_lpf_internal_8_dual_sse2(
    __m128i *p3, __m128i *q3, __m128i *p2, __m128i *q2, __m128i *p1,
    __m128i *q1, __m128i *p0, __m128i *q0, const uint8_t *_blimit0,
    const uint8_t *_limit0, const uint8_t *_thresh0,
    const uint8_t *_blimit1, const uint8_t *_limit1,
    const uint8_t *_thresh1, int32_t bd) {
    __m128i blimit0, limit0, thresh0;
    __m128i t80;
    __m128i mask, flat;
    __m128i work_a, op2, oq2, op1, op0, oq0, oq1;
    __m128i abs_p1q1, abs_p0q0, work0, work1, work2;

    const __m128i zero = _mm_setzero_si128();
    const __m128i four = _mm_set1_epi16(4);
    const __m128i one = _mm_set1_epi16(0x1);
    const __m128i ffff = _mm_cmpeq_epi16(one, one);

    get_limit_dual(_blimit0, _limit0, _thresh0, _blimit1, _limit1, _thresh1, bd,
        &blimit0, &limit0, &thresh0, &t80);

    abs_p0q0 = abs_diff16(*p0, *q0);
    abs_p1q1 = abs_diff16(*p1, *q1);

    abs_p0q0 = _mm_adds_epu16(abs_p0q0, abs_p0q0);
    abs_p1q1 = _mm_srli_epi16(abs_p1q1, 1);
    mask = _mm_subs_epu16(_mm_adds_epu16(abs_p0q0, abs_p1q1), blimit0);
    mask = _mm_xor_si128(_mm_cmpeq_epi16(mask, zero), ffff);
    // mask |= (abs(*p0 - q0) * 2 + abs(*p1 - q1) / 2  > blimit) * -1;

    // So taking maximums continues to work:
    mask = _mm_and_si128(mask, _mm_adds_epu16(limit0, one));

    work0 = _mm_max_epi16(abs_diff16(*p3, *p2), abs_diff16(*p2, *p1));
    work1 =
        _mm_max_epi16(abs_diff16(*p1, *p0), abs_diff16(*q1, *q0));  // tbu 4 flat
    work0 = _mm_max_epi16(work0, work1);
    work2 = _mm_max_epi16(abs_diff16(*q2, *q1), abs_diff16(*q2, *q3));
    work2 = _mm_max_epi16(work2, work0);
    mask = _mm_max_epi16(work2, mask);

    mask = _mm_subs_epu16(mask, limit0);
    mask = _mm_cmpeq_epi16(mask, zero);

    flat = _mm_max_epi16(abs_diff16(*p2, *p0), abs_diff16(*q2, *q0));
    flat = _mm_max_epi16(work1, flat);
    work0 = _mm_max_epi16(abs_diff16(*p3, *p0), abs_diff16(*q3, *q0));
    flat = _mm_max_epi16(work0, flat);

    if (bd == 8)
        flat = _mm_subs_epu16(flat, one);
    else if (bd == 10)
        flat = _mm_subs_epu16(flat, _mm_slli_epi16(one, 2));
    else  // bd == 12
        flat = _mm_subs_epu16(flat, _mm_slli_epi16(one, 4));

    flat = _mm_cmpeq_epi16(flat, zero);
    flat = _mm_and_si128(flat, mask);  // flat & mask

    {
        __m128i workp_a, workp_b;
        // Added before shift for rounding part of ROUND_POWER_OF_TWO

        // o*p2
        workp_a = _mm_add_epi16(_mm_add_epi16(*p3, *p3), _mm_add_epi16(*p2, *p1));
        workp_a = _mm_add_epi16(_mm_add_epi16(workp_a, four), *p0);
        workp_b = _mm_add_epi16(_mm_add_epi16(*q0, *p2), *p3);
        op2 = _mm_srli_epi16(_mm_add_epi16(workp_a, workp_b), 3);

        // o*p1
        workp_b = _mm_add_epi16(_mm_add_epi16(*q0, *q1), *p1);
        op1 = _mm_srli_epi16(_mm_add_epi16(workp_a, workp_b), 3);

        // o*p0
        workp_a = _mm_add_epi16(_mm_sub_epi16(workp_a, *p3), *q2);
        workp_b = _mm_add_epi16(_mm_sub_epi16(workp_b, *p1), *p0);
        op0 = _mm_srli_epi16(_mm_add_epi16(workp_a, workp_b), 3);

        // oq0
        workp_a = _mm_add_epi16(_mm_sub_epi16(workp_a, *p3), *q3);
        workp_b = _mm_add_epi16(_mm_sub_epi16(workp_b, *p0), *q0);
        oq0 = _mm_srli_epi16(_mm_add_epi16(workp_a, workp_b), 3);

        // oq1
        workp_a = _mm_add_epi16(_mm_sub_epi16(workp_a, *p2), *q3);
        workp_b = _mm_add_epi16(_mm_sub_epi16(workp_b, *q0), *q1);
        oq1 = _mm_srli_epi16(_mm_add_epi16(workp_a, workp_b), 3);

        // oq2
        workp_a = _mm_add_epi16(_mm_sub_epi16(workp_a, *p1), *q3);
        workp_b = _mm_add_epi16(_mm_sub_epi16(workp_b, *q1), *q2);
        oq2 = _mm_srli_epi16(_mm_add_epi16(workp_a, workp_b), 3);
    }

    // lp filter
    __m128i ps[2], qs[2], p[2], q[2];
    {
        p[0] = *p0;
        p[1] = *p1;
        q[0] = *q0;
        q[1] = *q1;
        // filter_mask and hev_mask
        highbd_filter4_dual_sse2(p, q, ps, qs, &mask, &thresh0, bd, &t80);
    }

    qs[0] = _mm_andnot_si128(flat, qs[0]);
    oq0 = _mm_and_si128(flat, oq0);
    *q0 = _mm_or_si128(qs[0], oq0);

    qs[1] = _mm_andnot_si128(flat, qs[1]);
    oq1 = _mm_and_si128(flat, oq1);
    *q1 = _mm_or_si128(qs[1], oq1);

    ps[0] = _mm_andnot_si128(flat, ps[0]);
    op0 = _mm_and_si128(flat, op0);
    *p0 = _mm_or_si128(ps[0], op0);

    ps[1] = _mm_andnot_si128(flat, ps[1]);
    op1 = _mm_and_si128(flat, op1);
    *p1 = _mm_or_si128(ps[1], op1);

    work_a = _mm_andnot_si128(flat, *q2);
    *q2 = _mm_and_si128(flat, oq2);
    *q2 = _mm_or_si128(work_a, *q2);

    work_a = _mm_andnot_si128(flat, *p2);
    *p2 = _mm_and_si128(flat, op2);
    *p2 = _mm_or_si128(work_a, *p2);
}

void aom_highbd_lpf_horizontal_8_sse2(uint16_t *s, int32_t p,
    const uint8_t *_blimit,
    const uint8_t *_limit,
    const uint8_t *_thresh, int32_t bd) {
    __m128i p2, p1, p0, q0, q1, q2, p3, q3;
    __m128i q1q0, p1p0;

    p3 = _mm_loadl_epi64((__m128i *)(s - 4 * p));
    q3 = _mm_loadl_epi64((__m128i *)(s + 3 * p));
    p2 = _mm_loadl_epi64((__m128i *)(s - 3 * p));
    q2 = _mm_loadl_epi64((__m128i *)(s + 2 * p));
    p1 = _mm_loadl_epi64((__m128i *)(s - 2 * p));
    q1 = _mm_loadl_epi64((__m128i *)(s + 1 * p));
    p0 = _mm_loadl_epi64((__m128i *)(s - 1 * p));
    q0 = _mm_loadl_epi64((__m128i *)(s + 0 * p));

    highbd_lpf_internal_8_sse2(&p3, &q3, &p2, &q2, &p1, &q1, &p0, &q0, &q1q0,
        &p1p0, _blimit, _limit, _thresh, bd);

    _mm_storel_epi64((__m128i *)(s - 3 * p), p2);
    _mm_storel_epi64((__m128i *)(s - 2 * p), _mm_srli_si128(p1p0, 8));
    _mm_storel_epi64((__m128i *)(s - 1 * p), p1p0);
    _mm_storel_epi64((__m128i *)(s + 0 * p), q1q0);
    _mm_storel_epi64((__m128i *)(s + 1 * p), _mm_srli_si128(q1q0, 8));
    _mm_storel_epi64((__m128i *)(s + 2 * p), q2);
}

void aom_highbd_lpf_horizontal_8_dual_sse2(
    uint16_t *s, int32_t p, const uint8_t *_blimit0, const uint8_t *_limit0,
    const uint8_t *_thresh0, const uint8_t *_blimit1, const uint8_t *_limit1,
    const uint8_t *_thresh1, int32_t bd) {
    __m128i p2, p1, p0, q0, q1, q2, p3, q3;

    p3 = _mm_loadu_si128((__m128i *)(s - 4 * p));
    q3 = _mm_loadu_si128((__m128i *)(s + 3 * p));
    p2 = _mm_loadu_si128((__m128i *)(s - 3 * p));
    q2 = _mm_loadu_si128((__m128i *)(s + 2 * p));
    p1 = _mm_loadu_si128((__m128i *)(s - 2 * p));
    q1 = _mm_loadu_si128((__m128i *)(s + 1 * p));
    p0 = _mm_loadu_si128((__m128i *)(s - 1 * p));
    q0 = _mm_loadu_si128((__m128i *)(s + 0 * p));

    highbd_lpf_internal_8_dual_sse2(&p3, &q3, &p2, &q2, &p1, &q1, &p0, &q0,
        _blimit0, _limit0, _thresh0, _blimit1,
        _limit1, _thresh1, bd);

    _mm_storeu_si128((__m128i *)(s - 3 * p), p2);
    _mm_storeu_si128((__m128i *)(s - 2 * p), p1);
    _mm_storeu_si128((__m128i *)(s - 1 * p), p0);
    _mm_storeu_si128((__m128i *)(s + 0 * p), q0);
    _mm_storeu_si128((__m128i *)(s + 1 * p), q1);
    _mm_storeu_si128((__m128i *)(s + 2 * p), q2);
}

static AOM_FORCE_INLINE void highbd_lpf_internal_4_sse2(
    __m128i *p1, __m128i *p0, __m128i *q0, __m128i *q1, __m128i *q1q0_out,
    __m128i *p1p0_out, const uint8_t *_blimit, const uint8_t *_limit,
    const uint8_t *_thresh, int32_t bd) {
    __m128i blimit, limit, thresh;
    __m128i mask, hev, flat;
    __m128i p1p0, q1q0;

    const __m128i zero = _mm_setzero_si128();

    __m128i abs_p0q0, abs_p1q1, abs_p1p0, abs_q1q0;

    const __m128i ffff = _mm_cmpeq_epi16(zero, zero);
    const __m128i one = _mm_set1_epi16(1);

    __m128i t80;
    get_limit(_blimit, _limit, _thresh, bd, &blimit, &limit, &thresh, &t80);

    p1p0 = _mm_unpacklo_epi64(*p0, *p1);
    q1q0 = _mm_unpacklo_epi64(*q0, *q1);

    abs_p1p0 = abs_diff16(*p1, *p0);
    abs_q1q0 = abs_diff16(*q1, *q0);

    abs_p0q0 = abs_diff16(p1p0, q1q0);
    abs_p1q1 = _mm_srli_si128(abs_p0q0, 8);

    // filter_mask and hev_mask
    flat = _mm_max_epi16(abs_p1p0, abs_q1q0);
    hev = _mm_subs_epu16(flat, thresh);
    hev = _mm_xor_si128(_mm_cmpeq_epi16(hev, zero), ffff);

    abs_p0q0 = _mm_adds_epu16(abs_p0q0, abs_p0q0);
    abs_p1q1 = _mm_srli_epi16(abs_p1q1, 1);
    mask = _mm_subs_epu16(_mm_adds_epu16(abs_p0q0, abs_p1q1), blimit);
    mask = _mm_xor_si128(_mm_cmpeq_epi16(mask, zero), ffff);
    // mask |= (abs(*p0 - *q0) * 2 + abs(*p1 - *q1) / 2  > blimit) * -1;
    // So taking maximums continues to work:
    mask = _mm_and_si128(mask, _mm_adds_epu16(limit, one));
    mask = _mm_max_epi16(flat, mask);

    mask = _mm_subs_epu16(mask, limit);
    mask = _mm_cmpeq_epi16(mask, zero);

    mask = _mm_unpacklo_epi64(mask, mask);
    hev = _mm_unpacklo_epi64(hev, hev);

    highbd_filter4_sse2(&p1p0, &q1q0, &hev, &mask, q1q0_out, p1p0_out, &t80, bd);
}

static AOM_FORCE_INLINE void highbd_lpf_internal_4_dual_sse2(
    __m128i *p1, __m128i *p0, __m128i *q0, __m128i *q1, __m128i *ps,
    __m128i *qs, const uint8_t *_blimit0, const uint8_t *_limit0,
    const uint8_t *_thresh0, const uint8_t *_blimit1, const uint8_t *_limit1,
    const uint8_t *_thresh1, int32_t bd) {
    __m128i blimit0, limit0, thresh0;
    __m128i mask, flat;
    __m128i p[2], q[2];

    const __m128i zero = _mm_setzero_si128();
    __m128i abs_p0q0 = abs_diff16(*q0, *p0);
    __m128i abs_p1q1 = abs_diff16(*q1, *p1);

    __m128i abs_p1p0 = abs_diff16(*p1, *p0);
    __m128i abs_q1q0 = abs_diff16(*q1, *q0);

    const __m128i ffff = _mm_cmpeq_epi16(abs_p1p0, abs_p1p0);
    const __m128i one = _mm_set1_epi16(1);

    __m128i t80;

    get_limit_dual(_blimit0, _limit0, _thresh0, _blimit1, _limit1, _thresh1, bd,
        &blimit0, &limit0, &thresh0, &t80);

    // filter_mask and hev_mask
    flat = _mm_max_epi16(abs_p1p0, abs_q1q0);

    abs_p0q0 = _mm_adds_epu16(abs_p0q0, abs_p0q0);
    abs_p1q1 = _mm_srli_epi16(abs_p1q1, 1);

    mask = _mm_subs_epu16(_mm_adds_epu16(abs_p0q0, abs_p1q1), blimit0);
    mask = _mm_xor_si128(_mm_cmpeq_epi16(mask, zero), ffff);
    // mask |= (abs(*p0 - *q0) * 2 + abs(*p1 - *q1) / 2  > blimit) * -1;
    // So taking maximums continues to work:
    mask = _mm_and_si128(mask, _mm_adds_epu16(limit0, one));
    mask = _mm_max_epi16(flat, mask);

    mask = _mm_subs_epu16(mask, limit0);
    mask = _mm_cmpeq_epi16(mask, zero);

    p[0] = *p0;
    p[1] = *p1;
    q[0] = *q0;
    q[1] = *q1;

    highbd_filter4_dual_sse2(p, q, ps, qs, &mask, &thresh0, bd, &t80);
}

void aom_highbd_lpf_horizontal_4_sse2(uint16_t *s, int32_t p,
    const uint8_t *_blimit,
    const uint8_t *_limit,
    const uint8_t *_thresh, int32_t bd) {
    __m128i p1p0, q1q0;
    __m128i p1 = _mm_loadl_epi64((__m128i *)(s - 2 * p));
    __m128i p0 = _mm_loadl_epi64((__m128i *)(s - 1 * p));
    __m128i q0 = _mm_loadl_epi64((__m128i *)(s - 0 * p));
    __m128i q1 = _mm_loadl_epi64((__m128i *)(s + 1 * p));

    highbd_lpf_internal_4_sse2(&p1, &p0, &q0, &q1, &q1q0, &p1p0, _blimit, _limit,
        _thresh, bd);

    _mm_storel_epi64((__m128i *)(s - 2 * p), _mm_srli_si128(p1p0, 8));
    _mm_storel_epi64((__m128i *)(s - 1 * p), p1p0);
    _mm_storel_epi64((__m128i *)(s + 0 * p), q1q0);
    _mm_storel_epi64((__m128i *)(s + 1 * p), _mm_srli_si128(q1q0, 8));
}

void aom_highbd_lpf_horizontal_4_dual_sse2(
    uint16_t *s, int32_t p, const uint8_t *_blimit0, const uint8_t *_limit0,
    const uint8_t *_thresh0, const uint8_t *_blimit1, const uint8_t *_limit1,
    const uint8_t *_thresh1, int32_t bd) {
    __m128i p1 = _mm_loadu_si128((__m128i *)(s - 2 * p));
    __m128i p0 = _mm_loadu_si128((__m128i *)(s - 1 * p));
    __m128i q0 = _mm_loadu_si128((__m128i *)(s - 0 * p));
    __m128i q1 = _mm_loadu_si128((__m128i *)(s + 1 * p));
    __m128i ps[2], qs[2];
    highbd_lpf_internal_4_dual_sse2(&p1, &p0, &q0, &q1, ps, qs, _blimit0, _limit0,
        _thresh0, _blimit1, _limit1, _thresh1, bd);

    _mm_storeu_si128((__m128i *)(s - 2 * p), ps[1]);
    _mm_storeu_si128((__m128i *)(s - 1 * p), ps[0]);
    _mm_storeu_si128((__m128i *)(s + 0 * p), qs[0]);
    _mm_storeu_si128((__m128i *)(s + 1 * p), qs[1]);
}

void aom_highbd_lpf_vertical_4_sse2(uint16_t *s, int32_t p, const uint8_t *blimit,
    const uint8_t *limit, const uint8_t *thresh,
    int32_t bd) {
    __m128i x0, x1, x2, x3, d0, d1, d2, d3, d4, d5, d6, d7;
    __m128i p1p0, q1q0;
    __m128i p1, q1;

    x0 = _mm_loadu_si128((__m128i *)(s - 4 + 0 * p));
    x1 = _mm_loadu_si128((__m128i *)(s - 4 + 1 * p));
    x2 = _mm_loadu_si128((__m128i *)(s - 4 + 2 * p));
    x3 = _mm_loadu_si128((__m128i *)(s - 4 + 3 * p));

    highbd_transpose4x8_8x4_sse2(&x0, &x1, &x2, &x3, &d0, &d1, &d2, &d3, &d4, &d5,
        &d6, &d7);

    highbd_lpf_internal_4_sse2(&d2, &d3, &d4, &d5, &q1q0, &p1p0, blimit, limit,
        thresh, bd);

    // transpose from 8x4 to 4x8
    p1 = _mm_srli_si128(p1p0, 8);
    q1 = _mm_srli_si128(q1q0, 8);
    highbd_transpose8x8_low_sse2(&d0, &d1, &p1, &p1p0, &q1q0, &q1, &d6, &d7, &d0,
        &d1, &d2, &d3);

    _mm_storeu_si128((__m128i *)(s - 4 + 0 * p), d0);
    _mm_storeu_si128((__m128i *)(s - 4 + 1 * p), d1);
    _mm_storeu_si128((__m128i *)(s - 4 + 2 * p), d2);
    _mm_storeu_si128((__m128i *)(s - 4 + 3 * p), d3);
}

void aom_highbd_lpf_vertical_4_dual_sse2(
    uint16_t *s, int32_t p, const uint8_t *blimit0, const uint8_t *limit0,
    const uint8_t *thresh0, const uint8_t *blimit1, const uint8_t *limit1,
    const uint8_t *thresh1, int32_t bd) {
    __m128i x0, x1, x2, x3, x4, x5, x6, x7;
    __m128i d0, d1, d2, d3, d4, d5, d6, d7;
    __m128i ps[2], qs[2];

    x0 = _mm_loadu_si128((__m128i *)(s - 4 + 0 * p));
    x1 = _mm_loadu_si128((__m128i *)(s - 4 + 1 * p));
    x2 = _mm_loadu_si128((__m128i *)(s - 4 + 2 * p));
    x3 = _mm_loadu_si128((__m128i *)(s - 4 + 3 * p));
    x4 = _mm_loadu_si128((__m128i *)(s - 4 + 4 * p));
    x5 = _mm_loadu_si128((__m128i *)(s - 4 + 5 * p));
    x6 = _mm_loadu_si128((__m128i *)(s - 4 + 6 * p));
    x7 = _mm_loadu_si128((__m128i *)(s - 4 + 7 * p));

    highbd_transpose8x8_sse2(&x0, &x1, &x2, &x3, &x4, &x5, &x6, &x7, &d0, &d1,
        &d2, &d3, &d4, &d5, &d6, &d7);

    highbd_lpf_internal_4_dual_sse2(&d2, &d3, &d4, &d5, ps, qs, blimit0, limit0,
        thresh0, blimit1, limit1, thresh1, bd);

    highbd_transpose8x8_sse2(&d0, &d1, &ps[1], &ps[0], &qs[0], &qs[1], &d6, &d7,
        &x0, &x1, &x2, &x3, &x4, &x5, &x6, &x7);

    _mm_storeu_si128((__m128i *)(s - 4 + 0 * p), x0);
    _mm_storeu_si128((__m128i *)(s - 4 + 1 * p), x1);
    _mm_storeu_si128((__m128i *)(s - 4 + 2 * p), x2);
    _mm_storeu_si128((__m128i *)(s - 4 + 3 * p), x3);
    _mm_storeu_si128((__m128i *)(s - 4 + 4 * p), x4);
    _mm_storeu_si128((__m128i *)(s - 4 + 5 * p), x5);
    _mm_storeu_si128((__m128i *)(s - 4 + 6 * p), x6);
    _mm_storeu_si128((__m128i *)(s - 4 + 7 * p), x7);
}

void aom_highbd_lpf_vertical_6_sse2(uint16_t *s, int32_t p, const uint8_t *blimit,
    const uint8_t *limit, const uint8_t *thresh,
    int32_t bd) {
    __m128i d0, d1, d2, d3, d4, d5;
    __m128i p2, p1, p0, q0, q1, q2;
    __m128i p1p0, q1q0;

    p2 = _mm_loadu_si128((__m128i *)((s - 3) + 0 * p));
    p1 = _mm_loadu_si128((__m128i *)((s - 3) + 1 * p));
    p0 = _mm_loadu_si128((__m128i *)((s - 3) + 2 * p));
    q0 = _mm_loadu_si128((__m128i *)((s - 3) + 3 * p));
    q1 = _mm_setzero_si128();
    q2 = _mm_setzero_si128();

    highbd_transpose6x6_sse2(&p2, &p1, &p0, &q0, &q1, &q2, &d0, &d1, &d2, &d3,
        &d4, &d5);

    highbd_lpf_internal_6_sse2(&d0, &d1, &d2, &d3, &d4, &d5, &p1p0, &q1q0, blimit,
        limit, thresh, bd);

    p0 = _mm_srli_si128(p1p0, 8);
    q0 = _mm_srli_si128(q1q0, 8);

    highbd_transpose6x6_sse2(&d0, &p0, &p1p0, &q1q0, &q0, &d5, &d0, &d1, &d2, &d3,
        &d4, &d5);

    _mm_storel_epi64((__m128i *)((s - 3) + 0 * p), d0);
    *(int32_t *)((s - 3) + 0 * p + 4) = _mm_cvtsi128_si32(_mm_srli_si128(d0, 8));
    _mm_storel_epi64((__m128i *)((s - 3) + 1 * p), d1);
    *(int32_t *)((s - 3) + 1 * p + 4) = _mm_cvtsi128_si32(_mm_srli_si128(d0, 12));
    _mm_storel_epi64((__m128i *)((s - 3) + 2 * p), d2);
    *(int32_t *)((s - 3) + 2 * p + 4) = _mm_cvtsi128_si32(_mm_srli_si128(d2, 8));
    _mm_storel_epi64((__m128i *)((s - 3) + 3 * p), d3);
    *(int32_t *)((s - 3) + 3 * p + 4) = _mm_cvtsi128_si32(_mm_srli_si128(d3, 8));
}

void aom_highbd_lpf_vertical_8_sse2(uint16_t *s, int32_t p, const uint8_t *blimit,
    const uint8_t *limit, const uint8_t *thresh,
    int32_t bd) {
    __m128i d0, d1, d2, d3, d4, d5, d6, d7;
    __m128i p2, p1, p0, p3, q0;
    __m128i q1q0, p1p0;

    p3 = _mm_loadu_si128((__m128i *)((s - 4) + 0 * p));
    p2 = _mm_loadu_si128((__m128i *)((s - 4) + 1 * p));
    p1 = _mm_loadu_si128((__m128i *)((s - 4) + 2 * p));
    p0 = _mm_loadu_si128((__m128i *)((s - 4) + 3 * p));

    highbd_transpose4x8_8x4_sse2(&p3, &p2, &p1, &p0, &d0, &d1, &d2, &d3, &d4, &d5,
        &d6, &d7);

    // Loop filtering
    highbd_lpf_internal_8_sse2(&d0, &d7, &d1, &d6, &d2, &d5, &d3, &d4, &q1q0,
        &p1p0, blimit, limit, thresh, bd);

    p0 = _mm_srli_si128(p1p0, 8);
    q0 = _mm_srli_si128(q1q0, 8);

    highbd_transpose8x8_low_sse2(&d0, &d1, &p0, &p1p0, &q1q0, &q0, &d6, &d7, &d0,
        &d1, &d2, &d3);

    _mm_storeu_si128((__m128i *)(s - 4 + 0 * p), d0);
    _mm_storeu_si128((__m128i *)(s - 4 + 1 * p), d1);
    _mm_storeu_si128((__m128i *)(s - 4 + 2 * p), d2);
    _mm_storeu_si128((__m128i *)(s - 4 + 3 * p), d3);
}

void aom_highbd_lpf_vertical_8_dual_sse2(
    uint16_t *s, int32_t p, const uint8_t *blimit0, const uint8_t *limit0,
    const uint8_t *thresh0, const uint8_t *blimit1, const uint8_t *limit1,
    const uint8_t *thresh1, int32_t bd) {
    __m128i x0, x1, x2, x3, x4, x5, x6, x7;
    __m128i d0, d1, d2, d3, d4, d5, d6, d7;

    x0 = _mm_loadu_si128((__m128i *)(s - 4 + 0 * p));
    x1 = _mm_loadu_si128((__m128i *)(s - 4 + 1 * p));
    x2 = _mm_loadu_si128((__m128i *)(s - 4 + 2 * p));
    x3 = _mm_loadu_si128((__m128i *)(s - 4 + 3 * p));
    x4 = _mm_loadu_si128((__m128i *)(s - 4 + 4 * p));
    x5 = _mm_loadu_si128((__m128i *)(s - 4 + 5 * p));
    x6 = _mm_loadu_si128((__m128i *)(s - 4 + 6 * p));
    x7 = _mm_loadu_si128((__m128i *)(s - 4 + 7 * p));

    highbd_transpose8x8_sse2(&x0, &x1, &x2, &x3, &x4, &x5, &x6, &x7, &d0, &d1,
        &d2, &d3, &d4, &d5, &d6, &d7);

    highbd_lpf_internal_8_dual_sse2(&d0, &d7, &d1, &d6, &d2, &d5, &d3, &d4,
        blimit0, limit0, thresh0, blimit1, limit1,
        thresh1, bd);

    highbd_transpose8x8_sse2(&d0, &d1, &d2, &d3, &d4, &d5, &d6, &d7, &x0, &x1,
        &x2, &x3, &x4, &x5, &x6, &x7);

    _mm_storeu_si128((__m128i *)(s - 4 + 0 * p), x0);
    _mm_storeu_si128((__m128i *)(s - 4 + 1 * p), x1);
    _mm_storeu_si128((__m128i *)(s - 4 + 2 * p), x2);
    _mm_storeu_si128((__m128i *)(s - 4 + 3 * p), x3);
    _mm_storeu_si128((__m128i *)(s - 4 + 4 * p), x4);
    _mm_storeu_si128((__m128i *)(s - 4 + 5 * p), x5);
    _mm_storeu_si128((__m128i *)(s - 4 + 6 * p), x6);
    _mm_storeu_si128((__m128i *)(s - 4 + 7 * p), x7);
}

void aom_highbd_lpf_vertical_14_sse2(uint16_t *s, int32_t pitch,
    const uint8_t *blimit,
    const uint8_t *limit,
    const uint8_t *thresh, int32_t bd) {
    __m128i q[7], p[7], pq[7];
    __m128i p6, p5, p4, p3;
    __m128i p6_2, p5_2, p4_2, p3_2;
    __m128i d0, d1, d2, d3;
    __m128i d0_2, d1_2, d2_2, d3_2, d7_2;

    p6 = _mm_loadu_si128((__m128i *)((s - 8) + 0 * pitch));
    p5 = _mm_loadu_si128((__m128i *)((s - 8) + 1 * pitch));
    p4 = _mm_loadu_si128((__m128i *)((s - 8) + 2 * pitch));
    p3 = _mm_loadu_si128((__m128i *)((s - 8) + 3 * pitch));

    highbd_transpose4x8_8x4_sse2(&p6, &p5, &p4, &p3, &d0, &p[6], &p[5], &p[4],
        &p[3], &p[2], &p[1], &p[0]);

    p6_2 = _mm_loadu_si128((__m128i *)(s + 0 * pitch));
    p5_2 = _mm_loadu_si128((__m128i *)(s + 1 * pitch));
    p4_2 = _mm_loadu_si128((__m128i *)(s + 2 * pitch));
    p3_2 = _mm_loadu_si128((__m128i *)(s + 3 * pitch));

    highbd_transpose4x8_8x4_sse2(&p6_2, &p5_2, &p4_2, &p3_2, &q[0], &q[1], &q[2],
        &q[3], &q[4], &q[5], &q[6], &d7_2);

    highbd_lpf_internal_14_sse2(p, q, pq, blimit, limit, thresh, bd);

    highbd_transpose8x8_low_sse2(&d0, &p[6], &pq[5], &pq[4], &pq[3], &pq[2],
        &pq[1], &pq[0], &d0, &d1, &d2, &d3);

    q[0] = _mm_srli_si128(pq[0], 8);
    q[1] = _mm_srli_si128(pq[1], 8);
    q[2] = _mm_srli_si128(pq[2], 8);
    q[3] = _mm_srli_si128(pq[3], 8);
    q[4] = _mm_srli_si128(pq[4], 8);
    q[5] = _mm_srli_si128(pq[5], 8);

    highbd_transpose8x8_low_sse2(&q[0], &q[1], &q[2], &q[3], &q[4], &q[5], &q[6],
        &d7_2, &d0_2, &d1_2, &d2_2, &d3_2);

    _mm_storeu_si128((__m128i *)(s - 8 + 0 * pitch), d0);
    _mm_storeu_si128((__m128i *)(s - 8 + 1 * pitch), d1);
    _mm_storeu_si128((__m128i *)(s - 8 + 2 * pitch), d2);
    _mm_storeu_si128((__m128i *)(s - 8 + 3 * pitch), d3);

    _mm_storeu_si128((__m128i *)(s + 0 * pitch), d0_2);
    _mm_storeu_si128((__m128i *)(s + 1 * pitch), d1_2);
    _mm_storeu_si128((__m128i *)(s + 2 * pitch), d2_2);
    _mm_storeu_si128((__m128i *)(s + 3 * pitch), d3_2);
}

void aom_highbd_lpf_vertical_14_dual_sse2(uint16_t *s, int32_t pitch,
    const uint8_t *blimit,
    const uint8_t *limit,
    const uint8_t *thresh, int32_t bd) {
    __m128i q[7], p[7];
    __m128i p6, p5, p4, p3, p2, p1, p0, q0;
    __m128i p6_2, p5_2, p4_2, p3_2, p2_2, p1_2, q0_2, p0_2;
    __m128i d0, d7;
    __m128i d0_out, d1_out, d2_out, d3_out, d4_out, d5_out, d6_out, d7_out;

    p6 = _mm_loadu_si128((__m128i *)((s - 8) + 0 * pitch));
    p5 = _mm_loadu_si128((__m128i *)((s - 8) + 1 * pitch));
    p4 = _mm_loadu_si128((__m128i *)((s - 8) + 2 * pitch));
    p3 = _mm_loadu_si128((__m128i *)((s - 8) + 3 * pitch));
    p2 = _mm_loadu_si128((__m128i *)((s - 8) + 4 * pitch));
    p1 = _mm_loadu_si128((__m128i *)((s - 8) + 5 * pitch));
    p0 = _mm_loadu_si128((__m128i *)((s - 8) + 6 * pitch));
    q0 = _mm_loadu_si128((__m128i *)((s - 8) + 7 * pitch));

    highbd_transpose8x8_sse2(&p6, &p5, &p4, &p3, &p2, &p1, &p0, &q0, &d0, &p[6],
        &p[5], &p[4], &p[3], &p[2], &p[1], &p[0]);

    p6_2 = _mm_loadu_si128((__m128i *)(s + 0 * pitch));
    p5_2 = _mm_loadu_si128((__m128i *)(s + 1 * pitch));
    p4_2 = _mm_loadu_si128((__m128i *)(s + 2 * pitch));
    p3_2 = _mm_loadu_si128((__m128i *)(s + 3 * pitch));
    p2_2 = _mm_loadu_si128((__m128i *)(s + 4 * pitch));
    p1_2 = _mm_loadu_si128((__m128i *)(s + 5 * pitch));
    p0_2 = _mm_loadu_si128((__m128i *)(s + 6 * pitch));
    q0_2 = _mm_loadu_si128((__m128i *)(s + 7 * pitch));

    highbd_transpose8x8_sse2(&p6_2, &p5_2, &p4_2, &p3_2, &p2_2, &p1_2, &p0_2,
        &q0_2, &q[0], &q[1], &q[2], &q[3], &q[4], &q[5],
        &q[6], &d7);

    highbd_lpf_internal_14_dual_sse2(p, q, blimit, limit, thresh, bd);

    highbd_transpose8x8_sse2(&d0, &p[6], &p[5], &p[4], &p[3], &p[2], &p[1], &p[0],
        &d0_out, &d1_out, &d2_out, &d3_out, &d4_out, &d5_out,
        &d6_out, &d7_out);

    _mm_storeu_si128((__m128i *)(s - 8 + 0 * pitch), d0_out);
    _mm_storeu_si128((__m128i *)(s - 8 + 1 * pitch), d1_out);
    _mm_storeu_si128((__m128i *)(s - 8 + 2 * pitch), d2_out);
    _mm_storeu_si128((__m128i *)(s - 8 + 3 * pitch), d3_out);
    _mm_storeu_si128((__m128i *)(s - 8 + 4 * pitch), d4_out);
    _mm_storeu_si128((__m128i *)(s - 8 + 5 * pitch), d5_out);
    _mm_storeu_si128((__m128i *)(s - 8 + 6 * pitch), d6_out);
    _mm_storeu_si128((__m128i *)(s - 8 + 7 * pitch), d7_out);

    highbd_transpose8x8_sse2(&q[0], &q[1], &q[2], &q[3], &q[4], &q[5], &q[6], &d7,
        &d0_out, &d1_out, &d2_out, &d3_out, &d4_out, &d5_out,
        &d6_out, &d7_out);

    _mm_storeu_si128((__m128i *)(s + 0 * pitch), d0_out);
    _mm_storeu_si128((__m128i *)(s + 1 * pitch), d1_out);
    _mm_storeu_si128((__m128i *)(s + 2 * pitch), d2_out);
    _mm_storeu_si128((__m128i *)(s + 3 * pitch), d3_out);
    _mm_storeu_si128((__m128i *)(s + 4 * pitch), d4_out);
    _mm_storeu_si128((__m128i *)(s + 5 * pitch), d5_out);
    _mm_storeu_si128((__m128i *)(s + 6 * pitch), d6_out);
    _mm_storeu_si128((__m128i *)(s + 7 * pitch), d7_out);
}





