/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include <string.h>

#include "EbDefinitions.h"
#include "immintrin.h"
#include "EbIntrinMacros_SSE2.h"
#include "EbIntraPrediction_AVX2.h"
#include "lpf_common_sse2.h"
#if INTRA_ASM
#include "aom_dsp_rtcd.h"
#endif
#ifndef _mm256_setr_m128i
#define _mm256_setr_m128i(/* __m128i */ hi, /* __m128i */ lo) \
    _mm256_insertf128_si256(_mm256_castsi128_si256(lo), (hi), 0x1)
#endif

#define MACRO_VERTICAL_LUMA_4(A, B, C) \
    *(uint32_t*)prediction_ptr = _mm_cvtsi128_si32(_mm_or_si128(_mm_and_si128(A, B), C)); \
    A = _mm_srli_si128(A, 1); \
    *(uint32_t*)(prediction_ptr + pStride) = _mm_cvtsi128_si32(_mm_or_si128(_mm_and_si128(A, B), C)); \
    A = _mm_srli_si128(A, 1);

#define _mm256_set_m128i(/* __m128i */ hi, /* __m128i */ lo) \
    _mm256_insertf128_si256(_mm256_castsi128_si256(lo), (hi), 0x1)

static INLINE void highbd_transpose16x4_8x8_sse2(__m128i *x, __m128i *d) {
    __m128i r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15;

    r0 = _mm_unpacklo_epi16(x[0], x[1]);
    r1 = _mm_unpacklo_epi16(x[2], x[3]);
    r2 = _mm_unpacklo_epi16(x[4], x[5]);
    r3 = _mm_unpacklo_epi16(x[6], x[7]);

    r4 = _mm_unpacklo_epi16(x[8], x[9]);
    r5 = _mm_unpacklo_epi16(x[10], x[11]);
    r6 = _mm_unpacklo_epi16(x[12], x[13]);
    r7 = _mm_unpacklo_epi16(x[14], x[15]);

    r8 = _mm_unpacklo_epi32(r0, r1);
    r9 = _mm_unpackhi_epi32(r0, r1);
    r10 = _mm_unpacklo_epi32(r2, r3);
    r11 = _mm_unpackhi_epi32(r2, r3);

    r12 = _mm_unpacklo_epi32(r4, r5);
    r13 = _mm_unpackhi_epi32(r4, r5);
    r14 = _mm_unpacklo_epi32(r6, r7);
    r15 = _mm_unpackhi_epi32(r6, r7);

    r0 = _mm_unpacklo_epi64(r8, r9);
    r1 = _mm_unpackhi_epi64(r8, r9);
    r2 = _mm_unpacklo_epi64(r10, r11);
    r3 = _mm_unpackhi_epi64(r10, r11);

    r4 = _mm_unpacklo_epi64(r12, r13);
    r5 = _mm_unpackhi_epi64(r12, r13);
    r6 = _mm_unpacklo_epi64(r14, r15);
    r7 = _mm_unpackhi_epi64(r14, r15);

    d[0] = _mm_unpacklo_epi64(r0, r2);
    d[1] = _mm_unpacklo_epi64(r4, r6);
    d[2] = _mm_unpacklo_epi64(r1, r3);
    d[3] = _mm_unpacklo_epi64(r5, r7);

    d[4] = _mm_unpackhi_epi64(r0, r2);
    d[5] = _mm_unpackhi_epi64(r4, r6);
    d[6] = _mm_unpackhi_epi64(r1, r3);
    d[7] = _mm_unpackhi_epi64(r5, r7);
}
static INLINE void highbd_transpose4x16_avx2(__m256i *x, __m256i *d) {
    __m256i w0, w1, w2, w3, ww0, ww1;

    w0 = _mm256_unpacklo_epi16(x[0], x[1]);  // 00 10 01 11 02 12 03 13
    w1 = _mm256_unpacklo_epi16(x[2], x[3]);  // 20 30 21 31 22 32 23 33
    w2 = _mm256_unpackhi_epi16(x[0], x[1]);  // 40 50 41 51 42 52 43 53
    w3 = _mm256_unpackhi_epi16(x[2], x[3]);  // 60 70 61 71 62 72 63 73

    ww0 = _mm256_unpacklo_epi32(w0, w1);  // 00 10 20 30 01 11 21 31
    ww1 = _mm256_unpacklo_epi32(w2, w3);  // 40 50 60 70 41 51 61 71

    d[0] = _mm256_unpacklo_epi64(ww0, ww1);  // 00 10 20 30 40 50 60 70
    d[1] = _mm256_unpackhi_epi64(ww0, ww1);  // 01 11 21 31 41 51 61 71

    ww0 = _mm256_unpackhi_epi32(w0, w1);  // 02 12 22 32 03 13 23 33
    ww1 = _mm256_unpackhi_epi32(w2, w3);  // 42 52 62 72 43 53 63 73

    d[2] = _mm256_unpacklo_epi64(ww0, ww1);  // 02 12 22 32 42 52 62 72
    d[3] = _mm256_unpackhi_epi64(ww0, ww1);  // 03 13 23 33 43 53 63 73
}

static INLINE void highbd_transpose8x16_16x8_avx2(__m256i *x, __m256i *d) {
    __m256i w0, w1, w2, w3, ww0, ww1;

    w0 = _mm256_unpacklo_epi16(x[0], x[1]);  // 00 10 01 11 02 12 03 13
    w1 = _mm256_unpacklo_epi16(x[2], x[3]);  // 20 30 21 31 22 32 23 33
    w2 = _mm256_unpacklo_epi16(x[4], x[5]);  // 40 50 41 51 42 52 43 53
    w3 = _mm256_unpacklo_epi16(x[6], x[7]);  // 60 70 61 71 62 72 63 73

    ww0 = _mm256_unpacklo_epi32(w0, w1);  // 00 10 20 30 01 11 21 31
    ww1 = _mm256_unpacklo_epi32(w2, w3);  // 40 50 60 70 41 51 61 71

    d[0] = _mm256_unpacklo_epi64(ww0, ww1);  // 00 10 20 30 40 50 60 70
    d[1] = _mm256_unpackhi_epi64(ww0, ww1);  // 01 11 21 31 41 51 61 71

    ww0 = _mm256_unpackhi_epi32(w0, w1);  // 02 12 22 32 03 13 23 33
    ww1 = _mm256_unpackhi_epi32(w2, w3);  // 42 52 62 72 43 53 63 73

    d[2] = _mm256_unpacklo_epi64(ww0, ww1);  // 02 12 22 32 42 52 62 72
    d[3] = _mm256_unpackhi_epi64(ww0, ww1);  // 03 13 23 33 43 53 63 73

    w0 = _mm256_unpackhi_epi16(x[0], x[1]);  // 04 14 05 15 06 16 07 17
    w1 = _mm256_unpackhi_epi16(x[2], x[3]);  // 24 34 25 35 26 36 27 37
    w2 = _mm256_unpackhi_epi16(x[4], x[5]);  // 44 54 45 55 46 56 47 57
    w3 = _mm256_unpackhi_epi16(x[6], x[7]);  // 64 74 65 75 66 76 67 77

    ww0 = _mm256_unpacklo_epi32(w0, w1);  // 04 14 24 34 05 15 25 35
    ww1 = _mm256_unpacklo_epi32(w2, w3);  // 44 54 64 74 45 55 65 75

    d[4] = _mm256_unpacklo_epi64(ww0, ww1);  // 04 14 24 34 44 54 64 74
    d[5] = _mm256_unpackhi_epi64(ww0, ww1);  // 05 15 25 35 45 55 65 75

    ww0 = _mm256_unpackhi_epi32(w0, w1);  // 06 16 26 36 07 17 27 37
    ww1 = _mm256_unpackhi_epi32(w2, w3);  // 46 56 66 76 47 57 67 77

    d[6] = _mm256_unpacklo_epi64(ww0, ww1);  // 06 16 26 36 46 56 66 76
    d[7] = _mm256_unpackhi_epi64(ww0, ww1);  // 07 17 27 37 47 57 67 77
}

static INLINE void highbd_transpose16x16_avx2(__m256i *x, __m256i *d) {
    __m256i w0, w1, w2, w3, ww0, ww1;
    __m256i dd[16];
    w0 = _mm256_unpacklo_epi16(x[0], x[1]);
    w1 = _mm256_unpacklo_epi16(x[2], x[3]);
    w2 = _mm256_unpacklo_epi16(x[4], x[5]);
    w3 = _mm256_unpacklo_epi16(x[6], x[7]);

    ww0 = _mm256_unpacklo_epi32(w0, w1);  //
    ww1 = _mm256_unpacklo_epi32(w2, w3);  //

    dd[0] = _mm256_unpacklo_epi64(ww0, ww1);
    dd[1] = _mm256_unpackhi_epi64(ww0, ww1);

    ww0 = _mm256_unpackhi_epi32(w0, w1);  //
    ww1 = _mm256_unpackhi_epi32(w2, w3);  //

    dd[2] = _mm256_unpacklo_epi64(ww0, ww1);
    dd[3] = _mm256_unpackhi_epi64(ww0, ww1);

    w0 = _mm256_unpackhi_epi16(x[0], x[1]);
    w1 = _mm256_unpackhi_epi16(x[2], x[3]);
    w2 = _mm256_unpackhi_epi16(x[4], x[5]);
    w3 = _mm256_unpackhi_epi16(x[6], x[7]);

    ww0 = _mm256_unpacklo_epi32(w0, w1);  //
    ww1 = _mm256_unpacklo_epi32(w2, w3);  //

    dd[4] = _mm256_unpacklo_epi64(ww0, ww1);
    dd[5] = _mm256_unpackhi_epi64(ww0, ww1);

    ww0 = _mm256_unpackhi_epi32(w0, w1);  //
    ww1 = _mm256_unpackhi_epi32(w2, w3);  //

    dd[6] = _mm256_unpacklo_epi64(ww0, ww1);
    dd[7] = _mm256_unpackhi_epi64(ww0, ww1);

    w0 = _mm256_unpacklo_epi16(x[8], x[9]);
    w1 = _mm256_unpacklo_epi16(x[10], x[11]);
    w2 = _mm256_unpacklo_epi16(x[12], x[13]);
    w3 = _mm256_unpacklo_epi16(x[14], x[15]);

    ww0 = _mm256_unpacklo_epi32(w0, w1);
    ww1 = _mm256_unpacklo_epi32(w2, w3);

    dd[8] = _mm256_unpacklo_epi64(ww0, ww1);
    dd[9] = _mm256_unpackhi_epi64(ww0, ww1);

    ww0 = _mm256_unpackhi_epi32(w0, w1);
    ww1 = _mm256_unpackhi_epi32(w2, w3);

    dd[10] = _mm256_unpacklo_epi64(ww0, ww1);
    dd[11] = _mm256_unpackhi_epi64(ww0, ww1);

    w0 = _mm256_unpackhi_epi16(x[8], x[9]);
    w1 = _mm256_unpackhi_epi16(x[10], x[11]);
    w2 = _mm256_unpackhi_epi16(x[12], x[13]);
    w3 = _mm256_unpackhi_epi16(x[14], x[15]);

    ww0 = _mm256_unpacklo_epi32(w0, w1);
    ww1 = _mm256_unpacklo_epi32(w2, w3);

    dd[12] = _mm256_unpacklo_epi64(ww0, ww1);
    dd[13] = _mm256_unpackhi_epi64(ww0, ww1);

    ww0 = _mm256_unpackhi_epi32(w0, w1);
    ww1 = _mm256_unpackhi_epi32(w2, w3);

    dd[14] = _mm256_unpacklo_epi64(ww0, ww1);
    dd[15] = _mm256_unpackhi_epi64(ww0, ww1);

    for (int32_t i = 0; i < 8; i++) {
        d[i] = _mm256_insertf128_si256(dd[i], _mm256_castsi256_si128(dd[i + 8]), 1);
        d[i + 8] = _mm256_insertf128_si256(dd[i + 8],
            _mm256_extracti128_si256(dd[i], 1), 0);
    }
}
static void _mm_storeh_epi64(__m128i * p, __m128i x)
{
    _mm_storeh_pd((double *)p, _mm_castsi128_pd(x));
}

EB_EXTERN void intra_mode_planar_avx2_intrin(
    const uint32_t   size,                       //input parameter, denotes the size of the current PU
    uint8_t         *ref_samples,                 //input parameter, pointer to the reference samples
    uint8_t         *prediction_ptr,              //output parameter, pointer to the prediction
    const uint32_t   prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
    const EbBool  skip)                       //skip half rows
{
    uint32_t topOffset = (size << 1) + 1;
    uint8_t  topRightPel, bottomLeftPel;
    uint32_t x, y;
    __m128i left, delta0, delta1, a0, a1, a2, a3, a4;
    __m256i lleft, ddelta0, ddelta1, aa0, aa1, aa2, aa3, aa4, cc0, cc1;

    // --------- Reference Samples Structure ---------
    // (ref_samples are similar as vertical mode)
    // ref_samples[0]        = Left[0]
    // ref_samples[1]        = Left[1]
    // ...
    // ref_samples[size]     = Left[size]
    // ... (arbitrary value)
    // ref_samples[2*size+1] = Top[0]
    // ref_samples[2*size+2] = Top[1]
    // ...
    // ref_samples[3*size+1] = Top[size]
    // -----------------------------------------------

    // Get above and left reference samples
    topRightPel = ref_samples[topOffset + size];
    bottomLeftPel = ref_samples[size];

    // Generate prediction
    if (size == 4) {
        a0 = _mm_set1_epi8(topRightPel);
        a1 = _mm_set1_epi16(bottomLeftPel);
        left = _mm_cvtsi32_si128(*(uint32_t *)ref_samples); // leftPel
        a2 = _mm_cvtsi32_si128(*(uint32_t *)(ref_samples + topOffset)); // topPel
        delta0 = _mm_unpacklo_epi8(a2, _mm_setzero_si128());
        delta0 = _mm_sub_epi16(a1, delta0);
        a2 = _mm_unpacklo_epi8(a2, a0);
        a2 = _mm_maddubs_epi16(a2, _mm_setr_epi8(4, 1, 4, 2, 4, 3, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0)); // (x + 1)*topRightPel + size * ref_samples[topOffset + x]
        a2 = _mm_add_epi16(a2, _mm_set1_epi16(4)); // (x + 1)*topRightPel + size * ref_samples[topOffset + x] + size
        if (skip) {
            a2 = _mm_add_epi16(a2, delta0);
            delta0 = _mm_slli_epi16(delta0, 1);
            a4 = _mm_shuffle_epi8(left, _mm_setzero_si128());
            a0 = _mm_add_epi16(a2, _mm_maddubs_epi16(a4, _mm_setr_epi8(3, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)));
            a0 = _mm_srai_epi16(a0, 3);
            a0 = _mm_packus_epi16(a0, a0);
            *(uint32_t *)prediction_ptr = _mm_cvtsi128_si32(a0);
            a2 = _mm_add_epi16(a2, delta0);
            left = _mm_srli_si128(left, 2);
            a4 = _mm_shuffle_epi8(left, _mm_setzero_si128());
            a0 = _mm_add_epi16(a2, _mm_maddubs_epi16(a4, _mm_setr_epi8(3, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)));
            a0 = _mm_srai_epi16(a0, 3);
            a0 = _mm_packus_epi16(a0, a0);
            *(uint32_t *)(prediction_ptr + 2 * prediction_buffer_stride) = _mm_cvtsi128_si32(a0);
        }
        else {
            for (y = 0; y < 4; y++) {
                a2 = _mm_add_epi16(a2, delta0);
                a4 = _mm_shuffle_epi8(left, _mm_setzero_si128());
                a0 = _mm_add_epi16(a2, _mm_maddubs_epi16(a4, _mm_setr_epi8(3, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)));
                a0 = _mm_srai_epi16(a0, 3);
                a0 = _mm_packus_epi16(a0, a0);
                *(uint32_t *)prediction_ptr = _mm_cvtsi128_si32(a0);
                prediction_ptr += prediction_buffer_stride;
                left = _mm_srli_si128(left, 1);
            }
        }
    }
    else if (size == 8) {
        a0 = _mm_set1_epi8(topRightPel);
        a1 = _mm_set1_epi16(bottomLeftPel);
        left = _mm_loadl_epi64((__m128i *)ref_samples); // leftPel
        a2 = _mm_loadl_epi64((__m128i *)(ref_samples + topOffset)); // topPel
        delta0 = _mm_unpacklo_epi8(a2, _mm_setzero_si128());
        delta0 = _mm_sub_epi16(a1, delta0);
        a2 = _mm_unpacklo_epi8(a2, a0);
        a2 = _mm_maddubs_epi16(a2, _mm_setr_epi8(8, 1, 8, 2, 8, 3, 8, 4, 8, 5, 8, 6, 8, 7, 8, 8)); // (x + 1)*topRightPel + size * ref_samples[topOffset + x]
        a2 = _mm_add_epi16(a2, _mm_set1_epi16(8)); // (x + 1)*topRightPel + size * ref_samples[topOffset + x] + size
        if (skip) {
            a2 = _mm_add_epi16(a2, delta0);
            delta0 = _mm_slli_epi16(delta0, 1);
            for (y = 0; y < 8; y += 2) {
                a4 = _mm_shuffle_epi8(left, _mm_setzero_si128());
                a0 = _mm_add_epi16(a2, _mm_maddubs_epi16(a4, _mm_setr_epi8(7, 0, 6, 0, 5, 0, 4, 0, 3, 0, 2, 0, 1, 0, 0, 0)));
                a0 = _mm_srai_epi16(a0, 4);
                a0 = _mm_packus_epi16(a0, a0);
                _mm_storel_epi64((__m128i *)prediction_ptr, a0);
                prediction_ptr += 2 * prediction_buffer_stride;
                a2 = _mm_add_epi16(a2, delta0);
                left = _mm_srli_si128(left, 2);
            }
        }
        else {
            for (y = 0; y < 8; y++) {
                a2 = _mm_add_epi16(a2, delta0);
                a4 = _mm_shuffle_epi8(left, _mm_setzero_si128());
                a0 = _mm_add_epi16(a2, _mm_maddubs_epi16(a4, _mm_setr_epi8(7, 0, 6, 0, 5, 0, 4, 0, 3, 0, 2, 0, 1, 0, 0, 0)));
                a0 = _mm_srai_epi16(a0, 4);
                a0 = _mm_packus_epi16(a0, a0);
                _mm_storel_epi64((__m128i *)prediction_ptr, a0);
                prediction_ptr += prediction_buffer_stride;
                left = _mm_srli_si128(left, 1);
            }
        }
    }
    else if (size == 16) {
        a0 = _mm_set1_epi8(topRightPel);
        a1 = _mm_set1_epi16(bottomLeftPel);
        left = _mm_loadu_si128((__m128i *)ref_samples); // leftPel
        a2 = _mm_loadu_si128((__m128i *)(ref_samples + topOffset)); // topPel
        delta0 = _mm_unpacklo_epi8(a2, _mm_setzero_si128());
        delta1 = _mm_unpackhi_epi8(a2, _mm_setzero_si128());
        delta0 = _mm_sub_epi16(a1, delta0);
        delta1 = _mm_sub_epi16(a1, delta1);
        a3 = _mm_unpackhi_epi8(a2, a0);
        a2 = _mm_unpacklo_epi8(a2, a0);
        a2 = _mm_maddubs_epi16(a2, _mm_setr_epi8(16, 1, 16, 2, 16, 3, 16, 4, 16, 5, 16, 6, 16, 7, 16, 8)); // (x + 1)*topRightPel + size * ref_samples[topOffset + x]
        a3 = _mm_maddubs_epi16(a3, _mm_setr_epi8(16, 9, 16, 10, 16, 11, 16, 12, 16, 13, 16, 14, 16, 15, 16, 16)); // (x + 1)*topRightPel + size * ref_samples[topOffset + x]
        a2 = _mm_add_epi16(a2, _mm_set1_epi16(16)); // (x + 1)*topRightPel + size * ref_samples[topOffset + x] + size
        a3 = _mm_add_epi16(a3, _mm_set1_epi16(16));
        if (skip) {
            a2 = _mm_add_epi16(a2, delta0);
            a3 = _mm_add_epi16(a3, delta1);
            delta0 = _mm_slli_epi16(delta0, 1);
            delta1 = _mm_slli_epi16(delta1, 1);
            for (y = 0; y < 16; y += 2) {
                a4 = _mm_shuffle_epi8(left, _mm_setzero_si128());
                a0 = _mm_add_epi16(a2, _mm_maddubs_epi16(a4, _mm_setr_epi8(15, 0, 14, 0, 13, 0, 12, 0, 11, 0, 10, 0, 9, 0, 8, 0)));
                a1 = _mm_add_epi16(a3, _mm_maddubs_epi16(a4, _mm_setr_epi8(7, 0, 6, 0, 5, 0, 4, 0, 3, 0, 2, 0, 1, 0, 0, 0)));
                a0 = _mm_srai_epi16(a0, 5);
                a1 = _mm_srai_epi16(a1, 5);
                a0 = _mm_packus_epi16(a0, a1);
                _mm_storeu_si128((__m128i *)prediction_ptr, a0);
                prediction_ptr += 2 * prediction_buffer_stride;
                a2 = _mm_add_epi16(a2, delta0);
                a3 = _mm_add_epi16(a3, delta1);
                left = _mm_srli_si128(left, 2);
            }
        }
        else {
            for (y = 0; y < 16; y++) {
                a2 = _mm_add_epi16(a2, delta0);
                a3 = _mm_add_epi16(a3, delta1);
                a4 = _mm_shuffle_epi8(left, _mm_setzero_si128());
                a0 = _mm_add_epi16(a2, _mm_maddubs_epi16(a4, _mm_setr_epi8(15, 0, 14, 0, 13, 0, 12, 0, 11, 0, 10, 0, 9, 0, 8, 0)));
                a1 = _mm_add_epi16(a3, _mm_maddubs_epi16(a4, _mm_setr_epi8(7, 0, 6, 0, 5, 0, 4, 0, 3, 0, 2, 0, 1, 0, 0, 0)));
                a0 = _mm_srai_epi16(a0, 5);
                a1 = _mm_srai_epi16(a1, 5);
                a0 = _mm_packus_epi16(a0, a1);
                _mm_storeu_si128((__m128i *)prediction_ptr, a0);
                prediction_ptr += prediction_buffer_stride;
                left = _mm_srli_si128(left, 1);
            }
        }
    }
    else if (size == 32) {
        aa0 = _mm256_set1_epi8(topRightPel);
        aa1 = _mm256_set1_epi16(bottomLeftPel);
        //lleft = _mm256_loadu_si256((__m256i *)ref_samples); // leftPel
        //lleft = _mm256_permute4x64_epi64(lleft, 0x44);
        aa2 = _mm256_loadu_si256((__m256i *)(ref_samples + topOffset)); // topPel
        ddelta0 = _mm256_unpacklo_epi8(aa2, _mm256_setzero_si256());
        ddelta1 = _mm256_unpackhi_epi8(aa2, _mm256_setzero_si256());
        ddelta0 = _mm256_sub_epi16(aa1, ddelta0);
        ddelta1 = _mm256_sub_epi16(aa1, ddelta1);
        aa3 = _mm256_unpackhi_epi8(aa2, aa0);
        aa2 = _mm256_unpacklo_epi8(aa2, aa0);
        aa2 = _mm256_maddubs_epi16(aa2, _mm256_setr_epi8(32, 1, 32, 2, 32, 3, 32, 4, 32, 5, 32, 6, 32, 7, 32, 8, 32, 17, 32, 18, 32, 19, 32, 20, 32, 21, 32, 22, 32, 23, 32, 24)); // (x + 1)*topRightPel + size * ref_samples[topOffset + x]
        aa3 = _mm256_maddubs_epi16(aa3, _mm256_setr_epi8(32, 9, 32, 10, 32, 11, 32, 12, 32, 13, 32, 14, 32, 15, 32, 16, 32, 25, 32, 26, 32, 27, 32, 28, 32, 29, 32, 30, 32, 31, 32, 32)); // (x + 1)*topRightPel + size * ref_samples[topOffset + x]
        aa2 = _mm256_add_epi16(aa2, _mm256_set1_epi16(32)); // (x + 1)*topRightPel + size * ref_samples[topOffset + x] + size
        aa3 = _mm256_add_epi16(aa3, _mm256_set1_epi16(32));
        cc0 = _mm256_setr_epi8(31, 0, 30, 0, 29, 0, 28, 0, 27, 0, 26, 0, 25, 0, 24, 0, 15, 0, 14, 0, 13, 0, 12, 0, 11, 0, 10, 0, 9, 0, 8, 0);
        cc1 = _mm256_setr_epi8(23, 0, 22, 0, 21, 0, 20, 0, 19, 0, 18, 0, 17, 0, 16, 0, 7, 0, 6, 0, 5, 0, 4, 0, 3, 0, 2, 0, 1, 0, 0, 0);
        if (skip) {
            lleft = _mm256_loadu_si256((__m256i *)ref_samples); // leftPel
            lleft = _mm256_shuffle_epi8(lleft, _mm256_setr_epi8(0, 2, 4, 6, 8, 10, 12, 14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 4, 6, 8, 10, 12, 14, 0, 0, 0, 0, 0, 0, 0, 0));
            lleft = _mm256_permute4x64_epi64(lleft, 0x88);
            aa2 = _mm256_add_epi16(aa2, ddelta0);
            aa3 = _mm256_add_epi16(aa3, ddelta1);
            ddelta0 = _mm256_slli_epi16(ddelta0, 1);
            ddelta1 = _mm256_slli_epi16(ddelta1, 1);
            for (y = 0; y < 16; y++) {
                aa4 = _mm256_shuffle_epi8(lleft, _mm256_setzero_si256());
                aa0 = _mm256_maddubs_epi16(aa4, cc0);
                aa1 = _mm256_maddubs_epi16(aa4, cc1);
                aa0 = _mm256_add_epi16(aa0, aa2);
                aa1 = _mm256_add_epi16(aa1, aa3);
                aa0 = _mm256_srai_epi16(aa0, 6);
                aa1 = _mm256_srai_epi16(aa1, 6);
                aa0 = _mm256_packus_epi16(aa0, aa1);
                _mm256_storeu_si256((__m256i *)prediction_ptr, aa0);
                prediction_ptr += 2 * prediction_buffer_stride;
                aa2 = _mm256_add_epi16(aa2, ddelta0);
                aa3 = _mm256_add_epi16(aa3, ddelta1);
                lleft = _mm256_srli_si256(lleft, 1);
            }
        }
        else {
            for (x = 0; x < 2; x++) {
                lleft = _mm256_loadu_si256((__m256i *)ref_samples); // leftPel
                lleft = _mm256_permute4x64_epi64(lleft, 0x44);
                ref_samples += 16;
                for (y = 0; y < 16; y++) {
                    aa2 = _mm256_add_epi16(aa2, ddelta0);
                    aa3 = _mm256_add_epi16(aa3, ddelta1);
                    aa4 = _mm256_shuffle_epi8(lleft, _mm256_setzero_si256());
                    aa0 = _mm256_maddubs_epi16(aa4, cc0);
                    aa1 = _mm256_maddubs_epi16(aa4, cc1);
                    aa0 = _mm256_add_epi16(aa0, aa2);
                    aa1 = _mm256_add_epi16(aa1, aa3);
                    aa0 = _mm256_srai_epi16(aa0, 6);
                    aa1 = _mm256_srai_epi16(aa1, 6);
                    aa0 = _mm256_packus_epi16(aa0, aa1);
                    _mm256_storeu_si256((__m256i *)prediction_ptr, aa0);
                    prediction_ptr += prediction_buffer_stride;
                    lleft = _mm256_srli_si256(lleft, 1);
                }
            }
        }
    }
}

// -----------------------------------------------------------------------------
// SMOOTH_PRED

// pixels[0]: above and below_pred interleave vector
// pixels[1]: left vector
// pixels[2]: right_pred vector
static INLINE void load_pixel_w4(const uint8_t *above, const uint8_t *left,
    int32_t height, __m128i *pixels) {
    __m128i d = _mm_loadl_epi64((const __m128i *)above);
    pixels[2] = _mm_set1_epi16((uint16_t)above[3]);
    pixels[1] = _mm_loadl_epi64((const __m128i *)left);

    const __m128i bp = _mm_set1_epi16((uint16_t)left[height - 1]);
    const __m128i zero = _mm_setzero_si128();
    d = _mm_unpacklo_epi8(d, zero);
    pixels[0] = _mm_unpacklo_epi16(d, bp);
}
// weights[0]: weights_h vector
// weights[1]: scale - weights_h vecotr
// weights[2]: weights_w and scale - weights_w interleave vector
static INLINE void load_weight_w4(const uint8_t *weight_array, int32_t height,
    __m128i *weights) {
    __m128i t = _mm_loadu_si128((const __m128i *)&weight_array[4]);
    const __m128i zero = _mm_setzero_si128();

    weights[0] = _mm_unpacklo_epi8(t, zero);
    const __m128i d = _mm_set1_epi16((uint16_t)(1 << sm_weight_log2_scale));
    weights[1] = _mm_sub_epi16(d, weights[0]);
    weights[2] = _mm_unpacklo_epi16(weights[0], weights[1]);

    if (height == 8) {
        t = _mm_srli_si128(t, 4);
        weights[0] = _mm_unpacklo_epi8(t, zero);
        weights[1] = _mm_sub_epi16(d, weights[0]);
    }
}
static INLINE void smooth_pred_4xh(const __m128i *pixel, const __m128i *weight,
    int32_t h, uint8_t *dst, ptrdiff_t stride) {
    const __m128i round = _mm_set1_epi32((1 << sm_weight_log2_scale));
    const __m128i one = _mm_set1_epi16(1);
    const __m128i inc = _mm_set1_epi16(0x202);
    const __m128i gat = _mm_set1_epi32(0xc080400);
    __m128i rep = _mm_set1_epi16((short)0x8000);
    __m128i d = _mm_set1_epi16(0x100);

    int32_t i;
    for (i = 0; i < h; ++i) {
        const __m128i wg_wg = _mm_shuffle_epi8(weight[0], d);
        const __m128i sc_sc = _mm_shuffle_epi8(weight[1], d);
        const __m128i wh_sc = _mm_unpacklo_epi16(wg_wg, sc_sc);
        __m128i s = _mm_madd_epi16(pixel[0], wh_sc);

        __m128i b = _mm_shuffle_epi8(pixel[1], rep);
        b = _mm_unpacklo_epi16(b, pixel[2]);
        __m128i sum = _mm_madd_epi16(b, weight[2]);

        sum = _mm_add_epi32(s, sum);
        sum = _mm_add_epi32(sum, round);
        sum = _mm_srai_epi32(sum, 1 + sm_weight_log2_scale);

        sum = _mm_shuffle_epi8(sum, gat);
        *(uint32_t *)dst = _mm_cvtsi128_si32(sum);
        dst += stride;

        rep = _mm_add_epi16(rep, one);
        d = _mm_add_epi16(d, inc);
    }
}
void intra_mode_planar_av1_avx2_intrin(
    const uint32_t   size,                       //input parameter, denotes the size of the current PU
    uint8_t         *ref_samples,                 //input parameter, pointer to the reference samples
    uint8_t         *dst,              //output parameter, pointer to the prediction
    const uint32_t   prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
    const EbBool  skip)
{

    uint32_t rowStride = skip ? 2 : 1;
    uint32_t leftOffset = 0;
    uint32_t topOffset = (size << 1) + 1;

    if (size == 4)
    {
        __m128i pixels[3];
        load_pixel_w4(&ref_samples[topOffset], &ref_samples[leftOffset], 4, pixels);

        __m128i weights[3];
        load_weight_w4(sm_weight_arrays, 4, weights);

        smooth_pred_4xh(pixels, weights, 4, dst, rowStride * prediction_buffer_stride);

    }
    else
    {
        const uint8_t *const sm_weights_w = sm_weight_arrays + size;
        const uint8_t *const sm_weights_h = sm_weight_arrays + size;
        const __m128i zero = _mm_setzero_si128();
        const __m128i scale_value =
            _mm_set1_epi16((uint16_t)(1 << sm_weight_log2_scale));
        const __m128i bottom_left = _mm_cvtsi32_si128((uint32_t)ref_samples[leftOffset + size - 1]);
        const __m128i dup16 =
            _mm_set_epi32(0x01000100, 0x01000100, 0x01000100, 0x01000100);
        const __m128i top_right =
            _mm_shuffle_epi8(_mm_cvtsi32_si128((uint32_t)ref_samples[topOffset + size - 1]), dup16);
        const __m128i gat = _mm_set_epi32(0, 0, 0xe0c0a08, 0x6040200);
        const __m128i round = _mm_set1_epi32((uint16_t)(1 << sm_weight_log2_scale));

        for (uint32_t y = 0; y < size; ++y) {
            const __m128i weights_y = _mm_cvtsi32_si128((uint32_t)sm_weights_h[y]);
            const __m128i left_y = _mm_cvtsi32_si128((uint32_t)ref_samples[leftOffset + y]);
            const __m128i scale_m_weights_y = _mm_sub_epi16(scale_value, weights_y);
            __m128i pred_scaled_bl = _mm_mullo_epi16(scale_m_weights_y, bottom_left);
            const __m128i wl_y =
                _mm_shuffle_epi32(_mm_unpacklo_epi16(weights_y, left_y), 0);
            pred_scaled_bl = _mm_add_epi32(pred_scaled_bl, round);
            pred_scaled_bl = _mm_shuffle_epi32(pred_scaled_bl, 0);

            for (uint32_t x = 0; x < size; x += 8) {
                const __m128i top_x = _mm_loadu_si128((const __m128i *)((ref_samples + topOffset + x)));
                const __m128i weights_x =
                    _mm_loadl_epi64((const __m128i *)(sm_weights_w + x));
                const __m128i tw_x = _mm_unpacklo_epi8(top_x, weights_x);
                const __m128i tw_x_lo = _mm_unpacklo_epi8(tw_x, zero);
                const __m128i tw_x_hi = _mm_unpackhi_epi8(tw_x, zero);

                __m128i pred_lo = _mm_madd_epi16(tw_x_lo, wl_y);
                __m128i pred_hi = _mm_madd_epi16(tw_x_hi, wl_y);

                const __m128i scale_m_weights_x =
                    _mm_sub_epi16(scale_value, _mm_unpacklo_epi8(weights_x, zero));
                const __m128i swxtr = _mm_mullo_epi16(scale_m_weights_x, top_right);
                const __m128i swxtr_lo = _mm_unpacklo_epi16(swxtr, zero);
                const __m128i swxtr_hi = _mm_unpackhi_epi16(swxtr, zero);

                pred_lo = _mm_add_epi32(pred_lo, pred_scaled_bl);
                pred_hi = _mm_add_epi32(pred_hi, pred_scaled_bl);

                pred_lo = _mm_add_epi32(pred_lo, swxtr_lo);
                pred_hi = _mm_add_epi32(pred_hi, swxtr_hi);

                pred_lo = _mm_srai_epi32(pred_lo, (1 + sm_weight_log2_scale));
                pred_hi = _mm_srai_epi32(pred_hi, (1 + sm_weight_log2_scale));

                __m128i pred = _mm_packus_epi16(pred_lo, pred_hi);
                pred = _mm_shuffle_epi8(pred, gat);
                _mm_storel_epi64((__m128i *)(dst + x), pred);
            }
            dst += rowStride * prediction_buffer_stride;
        }
    }
}


EB_EXTERN void intra_mode_angular_vertical_kernel_avx2_intrin(
    uint32_t         size,                       //input parameter, denotes the size of the current PU
    uint8_t         *ref_samp_main,                //input parameter, pointer to the reference samples
    uint8_t         *prediction_ptr,              //output parameter, pointer to the prediction
    uint32_t   prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
    const EbBool  skip,
    int32_t   intra_pred_angle)
{
    uint32_t row_index;
    uint32_t height = size;
    int32_t deltaSum = intra_pred_angle;
    int32_t deltaInt;
    uint32_t deltaFract;
    __m128i top0, top1, top2, sum0, sum1, a0, a1;
    __m256i ttop0, ttop1, ttop2, ssum0, ssum1, aa0;

    // --------- Reference Samples Structure ---------
    // ref_samp_main[-size+1] to ref_samp_main[-1] must be prepared (from bottom to top) for mode 19 to 25 (not required for mode 27 to 33)
    // ref_samp_main[0]      = TopLeft[0]
    // ref_samp_main[1]      = Top[0]
    // ref_samp_main[2]      = Top[1]
    // ...
    // ref_samp_main[size]   = Top[size-1]
    // ref_samp_main[size+1] = Top[size]     for mode 27 to 33 (not required for mode 19 to 25)
    // ...
    // ref_samp_main[2*size] = Top[2*size-1] for mode 27 to 33 (not required for mode 19 to 25)
    // -----------------------------------------------

    // Compute the prediction
    ref_samp_main += 1; // top0 sample
    if (skip) {
        height >>= 1;
        prediction_buffer_stride <<= 1;
        intra_pred_angle <<= 1;
    }

    if (size == 4) {
        for (row_index = 0; row_index < height; row_index += 2) {
            deltaInt = deltaSum >> 5;
            deltaFract = deltaSum & 31;
            a0 = _mm_unpacklo_epi8(_mm_cvtsi32_si128(32 - deltaFract), _mm_cvtsi32_si128(deltaFract));
            top0 = _mm_loadl_epi64((__m128i *)(ref_samp_main + deltaInt));
            deltaSum += intra_pred_angle;
            deltaInt = deltaSum >> 5;
            deltaFract = deltaSum & 31;
            a1 = _mm_unpacklo_epi8(_mm_cvtsi32_si128(32 - deltaFract), _mm_cvtsi32_si128(deltaFract));
            a0 = _mm_unpacklo_epi16(a0, a1);
            a0 = _mm_shuffle_epi8(a0, _mm_setr_epi8(0, 1, 0, 1, 0, 1, 0, 1, 2, 3, 2, 3, 2, 3, 2, 3));
            top0 = _mm_castps_si128(_mm_loadh_pi(_mm_castsi128_ps(top0), (__m64 *)(ref_samp_main + deltaInt)));
            top0 = _mm_shuffle_epi8(top0, _mm_setr_epi8(0, 1, 1, 2, 2, 3, 3, 4, 8, 9, 9, 10, 10, 11, 11, 12));
            sum0 = _mm_add_epi16(_mm_set1_epi16(16), _mm_maddubs_epi16(top0, a0));
            sum0 = _mm_srai_epi16(sum0, 5);
            sum0 = _mm_packus_epi16(sum0, sum0);
            *(uint32_t *)prediction_ptr = _mm_cvtsi128_si32(sum0);
            *(uint32_t *)(prediction_ptr + prediction_buffer_stride) = _mm_cvtsi128_si32(_mm_srli_si128(sum0, 4));
            prediction_ptr += 2 * prediction_buffer_stride;
            deltaSum += intra_pred_angle;
        }
    }
    else if (size == 8) {
        for (row_index = 0; row_index < height; row_index++) {
            deltaInt = deltaSum >> 5;
            deltaFract = deltaSum & 31;
            a0 = _mm_unpacklo_epi8(_mm_cvtsi32_si128(32 - deltaFract), _mm_cvtsi32_si128(deltaFract));
            a0 = _mm_shuffle_epi8(a0, _mm_setr_epi8(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1));
            top0 = _mm_loadu_si128((__m128i *)(ref_samp_main + deltaInt));
            top0 = _mm_shuffle_epi8(top0, _mm_setr_epi8(0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8));
            sum0 = _mm_add_epi16(_mm_set1_epi16(16), _mm_maddubs_epi16(top0, a0));
            sum0 = _mm_srai_epi16(sum0, 5);
            sum0 = _mm_packus_epi16(sum0, sum0);
            _mm_storel_epi64((__m128i *)prediction_ptr, sum0);
            prediction_ptr += prediction_buffer_stride;
            deltaSum += intra_pred_angle;
        }
    }
    else if (size == 16) {
        for (row_index = 0; row_index < height; row_index++) {
            deltaInt = deltaSum >> 5;
            deltaFract = deltaSum & 31;
            a0 = _mm_unpacklo_epi8(_mm_cvtsi32_si128(32 - deltaFract), _mm_cvtsi32_si128(deltaFract));
            a0 = _mm_shuffle_epi8(a0, _mm_setr_epi8(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1));
            top0 = _mm_loadu_si128((__m128i *)(ref_samp_main + deltaInt));
            top1 = _mm_loadu_si128((__m128i *)(ref_samp_main + deltaInt + 1));
            top2 = _mm_unpacklo_epi8(top0, top1);
            top0 = _mm_unpackhi_epi8(top0, top1);
            sum0 = _mm_add_epi16(_mm_set1_epi16(16), _mm_maddubs_epi16(top2, a0));
            sum1 = _mm_add_epi16(_mm_set1_epi16(16), _mm_maddubs_epi16(top0, a0));
            sum0 = _mm_srai_epi16(sum0, 5);
            sum1 = _mm_srai_epi16(sum1, 5);
            sum0 = _mm_packus_epi16(sum0, sum1);
            _mm_storeu_si128((__m128i *)prediction_ptr, sum0);
            prediction_ptr += prediction_buffer_stride;
            deltaSum += intra_pred_angle;
        }
    }
    else { // size == 32
        for (row_index = 0; row_index < height; row_index++) {
            deltaInt = deltaSum >> 5;
            deltaFract = deltaSum & 31;
            a0 = _mm_unpacklo_epi8(_mm_cvtsi32_si128(32 - deltaFract), _mm_cvtsi32_si128(deltaFract));
            a0 = _mm_shuffle_epi8(a0, _mm_setr_epi8(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1));
            aa0 = _mm256_setr_m128i(a0, a0);
            ttop0 = _mm256_loadu_si256((__m256i *)(ref_samp_main + deltaInt));
            ttop1 = _mm256_loadu_si256((__m256i *)(ref_samp_main + deltaInt + 1));
            ttop2 = _mm256_unpacklo_epi8(ttop0, ttop1);
            ttop0 = _mm256_unpackhi_epi8(ttop0, ttop1);
            ssum0 = _mm256_add_epi16(_mm256_set1_epi16(16), _mm256_maddubs_epi16(ttop2, aa0));
            ssum1 = _mm256_add_epi16(_mm256_set1_epi16(16), _mm256_maddubs_epi16(ttop0, aa0));
            ssum0 = _mm256_srai_epi16(ssum0, 5);
            ssum1 = _mm256_srai_epi16(ssum1, 5);
            ssum0 = _mm256_packus_epi16(ssum0, ssum1);
            _mm256_storeu_si256((__m256i *)prediction_ptr, ssum0);
            prediction_ptr += prediction_buffer_stride;
            deltaSum += intra_pred_angle;
        }
    }
}

EB_EXTERN void intra_mode_angular_horizontal_kernel_avx2_intrin(
    uint32_t         size,                       //input parameter, denotes the size of the current PU
    uint8_t         *ref_samp_main,                //input parameter, pointer to the reference samples
    uint8_t         *prediction_ptr,              //output parameter, pointer to the prediction
    uint32_t         prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
    const EbBool  skip,
    int32_t         intra_pred_angle)
{
    uint32_t row_index, colIndex;
    //uint32_t rowStride = skip ? 2 : 1;
    int32_t deltaSum = 0;
    int32_t deltaInt;
    uint32_t deltaFract;
    __m128i top0, top1, top2, sum0, sum1, a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11;
    __m256i ttop0, ttop1, ttop2, ssum0, ssum1, aa0;
    //__m256i aa1, aa2, aa3, aa4, aa5, aa6, aa7, aa8, aa9, aa10, aa11;
    uint8_t temp_buf[32 * 32];
    uint8_t *p = temp_buf;

    // --------- Reference Samples Structure ---------
    // ref_samp_main[-size+1] to ref_samp_main[-1] must be prepared (from right to left) for mode 11 to 17 (not required for mode 3 to 9)
    // ref_samp_main[0]      = TopLeft[0]
    // ref_samp_main[1]      = Left[0]
    // ref_samp_main[2]      = Left[1]
    // ...
    // ref_samp_main[size]   = Left[size-1]
    // ref_samp_main[size+1] = Left[size]     for mode 3 to 9 (not required for mode 11 to 17)
    // ...
    // ref_samp_main[2*size] = Left[2*size-1] for mode 3 to 9 (not required for mode 11 to 17)
    // -----------------------------------------------

    // Compute the prediction
    ref_samp_main += 1; // left sample

    if (skip) {
        prediction_buffer_stride <<= 1;
        if (size == 4) {
            deltaSum += intra_pred_angle;
            deltaInt = deltaSum >> 5;
            deltaFract = deltaSum & 31;
            a0 = _mm_unpacklo_epi8(_mm_cvtsi32_si128(32 - deltaFract), _mm_cvtsi32_si128(deltaFract));
            top0 = _mm_cvtsi32_si128(*(uint32_t *)(ref_samp_main + deltaInt));

            deltaSum += intra_pred_angle;
            deltaInt = deltaSum >> 5;
            deltaFract = deltaSum & 31;
            a1 = _mm_unpacklo_epi8(_mm_cvtsi32_si128(32 - deltaFract), _mm_cvtsi32_si128(deltaFract));
            a0 = _mm_unpacklo_epi16(a0, a1);
            top0 = _mm_unpacklo_epi32(top0, _mm_cvtsi32_si128(*(uint32_t *)(ref_samp_main + deltaInt)));

            deltaSum += intra_pred_angle;
            deltaInt = deltaSum >> 5;
            deltaFract = deltaSum & 31;
            a2 = _mm_unpacklo_epi8(_mm_cvtsi32_si128(32 - deltaFract), _mm_cvtsi32_si128(deltaFract));
            a4 = _mm_cvtsi32_si128(*(uint32_t *)(ref_samp_main + deltaInt));

            deltaSum += intra_pred_angle;
            deltaInt = deltaSum >> 5;
            deltaFract = deltaSum & 31;
            a3 = _mm_unpacklo_epi8(_mm_cvtsi32_si128(32 - deltaFract), _mm_cvtsi32_si128(deltaFract));
            a2 = _mm_unpacklo_epi16(a2, a3);
            a0 = _mm_unpacklo_epi32(a0, a2);
            a0 = _mm_unpacklo_epi16(a0, a0);
            a4 = _mm_unpacklo_epi32(a4, _mm_cvtsi32_si128(*(uint32_t *)(ref_samp_main + deltaInt)));

            top0 = _mm_unpacklo_epi64(top0, a4);
            sum0 = _mm_add_epi16(_mm_set1_epi16(16), _mm_maddubs_epi16(top0, a0));
            sum0 = _mm_srai_epi16(sum0, 5);
            sum0 = _mm_packus_epi16(sum0, sum0);
            sum0 = _mm_shuffle_epi8(sum0, _mm_setr_epi8(0, 2, 4, 6, 1, 3, 5, 7, 2, 6, 10, 14, 3, 7, 11, 15));
            *(uint32_t *)prediction_ptr = _mm_cvtsi128_si32(sum0); sum0 = _mm_srli_si128(sum0, 4); prediction_ptr += prediction_buffer_stride;
            *(uint32_t *)prediction_ptr = _mm_cvtsi128_si32(sum0);
        }
        else if (size == 8) {
            for (colIndex = 0; colIndex < size; colIndex++) {
                deltaSum += intra_pred_angle;
                deltaInt = deltaSum >> 5;
                deltaFract = deltaSum & 31;
                a0 = _mm_unpacklo_epi8(_mm_cvtsi32_si128(32 - deltaFract), _mm_cvtsi32_si128(deltaFract));
                a0 = _mm_shuffle_epi8(a0, _mm_setr_epi8(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1));
                top0 = _mm_loadl_epi64((__m128i *)(ref_samp_main + deltaInt));
                sum0 = _mm_add_epi16(_mm_set1_epi16(16), _mm_maddubs_epi16(top0, a0));
                sum0 = _mm_srai_epi16(sum0, 5);
                sum0 = _mm_packus_epi16(sum0, sum0);
                *(uint32_t *)p = _mm_cvtsi128_si32(sum0);
                p += 4;
            }
            a0 = _mm_loadu_si128((__m128i *)temp_buf);
            a0 = _mm_shuffle_epi8(a0, _mm_setr_epi8(0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15));
            a1 = _mm_loadu_si128((__m128i *)(temp_buf + 16));
            a1 = _mm_shuffle_epi8(a1, _mm_setr_epi8(0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15));
            a2 = _mm_unpackhi_epi32(a0, a1);
            a0 = _mm_unpacklo_epi32(a0, a1);
            a1 = _mm_unpackhi_epi64(a0, a2);
            a0 = _mm_unpacklo_epi64(a0, a2);
            _mm_storel_epi64((__m128i *)prediction_ptr, a0); prediction_ptr += prediction_buffer_stride;
            _mm_storel_epi64((__m128i *)prediction_ptr, a1); prediction_ptr += prediction_buffer_stride;
            _mm_storeh_epi64((__m128i *)prediction_ptr, a0); prediction_ptr += prediction_buffer_stride;
            _mm_storeh_epi64((__m128i *)prediction_ptr, a1);
        }
        else if (size == 16) {
            for (colIndex = 0; colIndex < size; colIndex++) {
                deltaSum += intra_pred_angle;
                deltaInt = deltaSum >> 5;
                deltaFract = deltaSum & 31;
                a0 = _mm_unpacklo_epi8(_mm_cvtsi32_si128(32 - deltaFract), _mm_cvtsi32_si128(deltaFract));
                a0 = _mm_shuffle_epi8(a0, _mm_setr_epi8(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1));
                top0 = _mm_loadu_si128((__m128i *)(ref_samp_main + deltaInt));
                top0 = _mm_xor_si128(top0, _mm_setzero_si128()); // Redundant code. Visual Studio Express 2013 Release build messes up the order of operands in _mm_maddubs_epi16(), and hard to work around.
                sum0 = _mm_add_epi16(_mm_set1_epi16(16), _mm_maddubs_epi16(top0, a0));
                sum0 = _mm_srai_epi16(sum0, 5);
                sum0 = _mm_packus_epi16(sum0, sum0);
                _mm_storel_epi64((__m128i *)p, sum0);
                p += 8;
            }
            p = temp_buf;
            a0 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x00)), _mm_loadl_epi64((__m128i *)(p + 0x08)));
            a1 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x10)), _mm_loadl_epi64((__m128i *)(p + 0x18)));
            a2 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x20)), _mm_loadl_epi64((__m128i *)(p + 0x28)));
            a3 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x30)), _mm_loadl_epi64((__m128i *)(p + 0x38)));
            a4 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x40)), _mm_loadl_epi64((__m128i *)(p + 0x48)));
            a5 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x50)), _mm_loadl_epi64((__m128i *)(p + 0x58)));
            a6 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x60)), _mm_loadl_epi64((__m128i *)(p + 0x68)));
            a7 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x70)), _mm_loadl_epi64((__m128i *)(p + 0x78)));

            a8 = _mm_unpackhi_epi16(a0, a1);
            a0 = _mm_unpacklo_epi16(a0, a1);
            a9 = _mm_unpackhi_epi16(a2, a3);
            a2 = _mm_unpacklo_epi16(a2, a3);
            a10 = _mm_unpackhi_epi16(a4, a5);
            a4 = _mm_unpacklo_epi16(a4, a5);
            a11 = _mm_unpackhi_epi16(a6, a7);
            a6 = _mm_unpacklo_epi16(a6, a7);

            a1 = _mm_unpackhi_epi32(a0, a2);
            a0 = _mm_unpacklo_epi32(a0, a2);
            a3 = _mm_unpackhi_epi32(a4, a6);
            a4 = _mm_unpacklo_epi32(a4, a6);
            a5 = _mm_unpackhi_epi32(a8, a9);
            a8 = _mm_unpacklo_epi32(a8, a9);
            a7 = _mm_unpackhi_epi32(a10, a11);
            a10 = _mm_unpacklo_epi32(a10, a11);

            a2 = _mm_unpackhi_epi64(a0, a4);
            a0 = _mm_unpacklo_epi64(a0, a4);
            a6 = _mm_unpackhi_epi64(a8, a10);
            a8 = _mm_unpacklo_epi64(a8, a10);
            a9 = _mm_unpackhi_epi64(a1, a3);
            a1 = _mm_unpacklo_epi64(a1, a3);
            a11 = _mm_unpackhi_epi64(a5, a7);
            a5 = _mm_unpacklo_epi64(a5, a7);

            _mm_storeu_si128((__m128i *)prediction_ptr, a0);  prediction_ptr += prediction_buffer_stride;
            _mm_storeu_si128((__m128i *)prediction_ptr, a2);  prediction_ptr += prediction_buffer_stride;
            _mm_storeu_si128((__m128i *)prediction_ptr, a1);  prediction_ptr += prediction_buffer_stride;
            _mm_storeu_si128((__m128i *)prediction_ptr, a9);  prediction_ptr += prediction_buffer_stride;
            _mm_storeu_si128((__m128i *)prediction_ptr, a8);  prediction_ptr += prediction_buffer_stride;
            _mm_storeu_si128((__m128i *)prediction_ptr, a6);  prediction_ptr += prediction_buffer_stride;
            _mm_storeu_si128((__m128i *)prediction_ptr, a5);  prediction_ptr += prediction_buffer_stride;
            _mm_storeu_si128((__m128i *)prediction_ptr, a11); prediction_ptr += prediction_buffer_stride;
        }
        else { // size == 32
            for (colIndex = 0; colIndex < size; colIndex++) {
                deltaSum += intra_pred_angle;
                deltaInt = deltaSum >> 5;
                deltaFract = deltaSum & 31;
                a0 = _mm_unpacklo_epi8(_mm_cvtsi32_si128(32 - deltaFract), _mm_cvtsi32_si128(deltaFract));
                a0 = _mm_shuffle_epi8(a0, _mm_setr_epi8(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1));
                top0 = _mm_loadu_si128((__m128i *)(ref_samp_main + deltaInt));
                top1 = _mm_loadu_si128((__m128i *)(ref_samp_main + deltaInt + 16));
                top0 = _mm_xor_si128(top0, _mm_setzero_si128()); // Redundant code. Visual Studio Express 2013 Release build messes up the order of operands in _mm_maddubs_epi16(), and hard to work around.
                top1 = _mm_xor_si128(top1, _mm_setzero_si128()); // Redundant code. Visual Studio Express 2013 Release build messes up the order of operands in _mm_maddubs_epi16(), and hard to work around.
                sum0 = _mm_add_epi16(_mm_set1_epi16(16), _mm_maddubs_epi16(top0, a0));
                sum1 = _mm_add_epi16(_mm_set1_epi16(16), _mm_maddubs_epi16(top1, a0));
                sum0 = _mm_srai_epi16(sum0, 5);
                sum1 = _mm_srai_epi16(sum1, 5);
                sum0 = _mm_packus_epi16(sum0, sum1);
                _mm_storeu_si128((__m128i *)p, sum0);
                p += 16;
            }
            p = temp_buf;
            for (colIndex = 0; colIndex < 2; colIndex++) {
                for (row_index = 0; row_index < 2; row_index++) {
                    a0 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x00)), _mm_loadl_epi64((__m128i *)(p + 0x10)));
                    a1 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x20)), _mm_loadl_epi64((__m128i *)(p + 0x30)));
                    a2 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x40)), _mm_loadl_epi64((__m128i *)(p + 0x50)));
                    a3 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x60)), _mm_loadl_epi64((__m128i *)(p + 0x70)));
                    a4 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x80)), _mm_loadl_epi64((__m128i *)(p + 0x90)));
                    a5 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0xA0)), _mm_loadl_epi64((__m128i *)(p + 0xB0)));
                    a6 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0xC0)), _mm_loadl_epi64((__m128i *)(p + 0xD0)));
                    a7 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0xE0)), _mm_loadl_epi64((__m128i *)(p + 0xF0)));

                    a8 = _mm_unpackhi_epi16(a0, a1);
                    a0 = _mm_unpacklo_epi16(a0, a1);
                    a9 = _mm_unpackhi_epi16(a2, a3);
                    a2 = _mm_unpacklo_epi16(a2, a3);
                    a10 = _mm_unpackhi_epi16(a4, a5);
                    a4 = _mm_unpacklo_epi16(a4, a5);
                    a11 = _mm_unpackhi_epi16(a6, a7);
                    a6 = _mm_unpacklo_epi16(a6, a7);

                    a1 = _mm_unpackhi_epi32(a0, a2);
                    a0 = _mm_unpacklo_epi32(a0, a2);
                    a3 = _mm_unpackhi_epi32(a4, a6);
                    a4 = _mm_unpacklo_epi32(a4, a6);
                    a5 = _mm_unpackhi_epi32(a8, a9);
                    a8 = _mm_unpacklo_epi32(a8, a9);
                    a7 = _mm_unpackhi_epi32(a10, a11);
                    a10 = _mm_unpacklo_epi32(a10, a11);

                    a2 = _mm_unpackhi_epi64(a0, a4);
                    a0 = _mm_unpacklo_epi64(a0, a4);
                    a6 = _mm_unpackhi_epi64(a8, a10);
                    a8 = _mm_unpacklo_epi64(a8, a10);
                    a9 = _mm_unpackhi_epi64(a1, a3);
                    a1 = _mm_unpacklo_epi64(a1, a3);
                    a11 = _mm_unpackhi_epi64(a5, a7);
                    a5 = _mm_unpacklo_epi64(a5, a7);

                    _mm_storeu_si128((__m128i *)prediction_ptr, a0);  prediction_ptr += prediction_buffer_stride;
                    _mm_storeu_si128((__m128i *)prediction_ptr, a2);  prediction_ptr += prediction_buffer_stride;
                    _mm_storeu_si128((__m128i *)prediction_ptr, a1);  prediction_ptr += prediction_buffer_stride;
                    _mm_storeu_si128((__m128i *)prediction_ptr, a9);  prediction_ptr += prediction_buffer_stride;
                    _mm_storeu_si128((__m128i *)prediction_ptr, a8);  prediction_ptr += prediction_buffer_stride;
                    _mm_storeu_si128((__m128i *)prediction_ptr, a6);  prediction_ptr += prediction_buffer_stride;
                    _mm_storeu_si128((__m128i *)prediction_ptr, a5);  prediction_ptr += prediction_buffer_stride;
                    _mm_storeu_si128((__m128i *)prediction_ptr, a11); prediction_ptr += prediction_buffer_stride;
                    p += 8;
                }
                p = temp_buf + 0x100;
                prediction_ptr -= 16 * prediction_buffer_stride - 16;
            }
        }
    }
    else {
        if (size == 4) {
            for (colIndex = 0; colIndex < size; colIndex += 2) {
                deltaSum += intra_pred_angle;
                deltaInt = deltaSum >> 5;
                deltaFract = deltaSum & 31;
                a0 = _mm_unpacklo_epi8(_mm_cvtsi32_si128(32 - deltaFract), _mm_cvtsi32_si128(deltaFract));
                top0 = _mm_loadl_epi64((__m128i *)(ref_samp_main + deltaInt));
                deltaSum += intra_pred_angle;
                deltaInt = deltaSum >> 5;
                deltaFract = deltaSum & 31;
                a1 = _mm_unpacklo_epi8(_mm_cvtsi32_si128(32 - deltaFract), _mm_cvtsi32_si128(deltaFract));
                a0 = _mm_unpacklo_epi16(a0, a1);
                a0 = _mm_shuffle_epi8(a0, _mm_setr_epi8(0, 1, 0, 1, 0, 1, 0, 1, 2, 3, 2, 3, 2, 3, 2, 3));
                top0 = _mm_castps_si128(_mm_loadh_pi(_mm_castsi128_ps(top0), (__m64 *)(ref_samp_main + deltaInt)));
                top0 = _mm_shuffle_epi8(top0, _mm_setr_epi8(0, 1, 1, 2, 2, 3, 3, 4, 8, 9, 9, 10, 10, 11, 11, 12));
                sum0 = _mm_add_epi16(_mm_set1_epi16(16), _mm_maddubs_epi16(top0, a0));
                sum0 = _mm_srai_epi16(sum0, 5);
                sum0 = _mm_packus_epi16(sum0, sum0);
                *(uint32_t *)p = _mm_cvtsi128_si32(sum0);
                *(uint32_t *)(p + 4) = _mm_cvtsi128_si32(_mm_srli_si128(sum0, 4));
                p += 8;
            }
            a0 = _mm_loadu_si128((__m128i *)temp_buf);
            a0 = _mm_shuffle_epi8(a0, _mm_setr_epi8(0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15));
            *(uint32_t *)prediction_ptr = _mm_cvtsi128_si32(a0); a0 = _mm_srli_si128(a0, 4); prediction_ptr += prediction_buffer_stride;
            *(uint32_t *)prediction_ptr = _mm_cvtsi128_si32(a0); a0 = _mm_srli_si128(a0, 4); prediction_ptr += prediction_buffer_stride;
            *(uint32_t *)prediction_ptr = _mm_cvtsi128_si32(a0); a0 = _mm_srli_si128(a0, 4); prediction_ptr += prediction_buffer_stride;
            *(uint32_t *)prediction_ptr = _mm_cvtsi128_si32(a0);
        }
        else if (size == 8) {
            for (colIndex = 0; colIndex < size; colIndex++) {
                deltaSum += intra_pred_angle;
                deltaInt = deltaSum >> 5;
                deltaFract = deltaSum & 31;
                a0 = _mm_unpacklo_epi8(_mm_cvtsi32_si128(32 - deltaFract), _mm_cvtsi32_si128(deltaFract));
                a0 = _mm_shuffle_epi8(a0, _mm_setr_epi8(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1));
                top0 = _mm_loadu_si128((__m128i *)(ref_samp_main + deltaInt));
                top0 = _mm_shuffle_epi8(top0, _mm_setr_epi8(0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8));
                sum0 = _mm_add_epi16(_mm_set1_epi16(16), _mm_maddubs_epi16(top0, a0));
                sum0 = _mm_srai_epi16(sum0, 5);
                sum0 = _mm_packus_epi16(sum0, sum0);
                _mm_storel_epi64((__m128i *)p, sum0);
                p += 8;
            }
            a0 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(temp_buf + 0x00)), _mm_loadl_epi64((__m128i *)(temp_buf + 0x08)));
            a1 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(temp_buf + 0x10)), _mm_loadl_epi64((__m128i *)(temp_buf + 0x18)));
            a2 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(temp_buf + 0x20)), _mm_loadl_epi64((__m128i *)(temp_buf + 0x28)));
            a3 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(temp_buf + 0x30)), _mm_loadl_epi64((__m128i *)(temp_buf + 0x38)));
            a4 = _mm_unpackhi_epi16(a0, a1);
            a0 = _mm_unpacklo_epi16(a0, a1);
            a5 = _mm_unpackhi_epi16(a2, a3);
            a2 = _mm_unpacklo_epi16(a2, a3);
            a1 = _mm_unpackhi_epi32(a0, a2);
            a0 = _mm_unpacklo_epi32(a0, a2);
            a3 = _mm_unpackhi_epi32(a4, a5);
            a4 = _mm_unpacklo_epi32(a4, a5);
            _mm_storel_epi64((__m128i *)prediction_ptr, a0); prediction_ptr += prediction_buffer_stride;
            _mm_storeh_epi64((__m128i *)prediction_ptr, a0); prediction_ptr += prediction_buffer_stride;
            _mm_storel_epi64((__m128i *)prediction_ptr, a1); prediction_ptr += prediction_buffer_stride;
            _mm_storeh_epi64((__m128i *)prediction_ptr, a1); prediction_ptr += prediction_buffer_stride;
            _mm_storel_epi64((__m128i *)prediction_ptr, a4); prediction_ptr += prediction_buffer_stride;
            _mm_storeh_epi64((__m128i *)prediction_ptr, a4); prediction_ptr += prediction_buffer_stride;
            _mm_storel_epi64((__m128i *)prediction_ptr, a3); prediction_ptr += prediction_buffer_stride;
            _mm_storeh_epi64((__m128i *)prediction_ptr, a3);
        }
        else if (size == 16) {
            for (colIndex = 0; colIndex < size; colIndex++) {
                deltaSum += intra_pred_angle;
                deltaInt = deltaSum >> 5;
                deltaFract = deltaSum & 31;
                a0 = _mm_unpacklo_epi8(_mm_cvtsi32_si128(32 - deltaFract), _mm_cvtsi32_si128(deltaFract));
                a0 = _mm_shuffle_epi8(a0, _mm_setr_epi8(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1));
                top0 = _mm_loadu_si128((__m128i *)(ref_samp_main + deltaInt));
                top1 = _mm_loadu_si128((__m128i *)(ref_samp_main + deltaInt + 1));
                top2 = _mm_unpacklo_epi8(top0, top1);
                top0 = _mm_unpackhi_epi8(top0, top1);
                sum0 = _mm_add_epi16(_mm_set1_epi16(16), _mm_maddubs_epi16(top2, a0));
                sum1 = _mm_add_epi16(_mm_set1_epi16(16), _mm_maddubs_epi16(top0, a0));
                sum0 = _mm_srai_epi16(sum0, 5);
                sum1 = _mm_srai_epi16(sum1, 5);
                sum0 = _mm_packus_epi16(sum0, sum1);
                _mm_storeu_si128((__m128i *)p, sum0);
                p += 16;
            }
            p = temp_buf;
            for (colIndex = 0; colIndex < 2; colIndex++) {
                a0 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x00)), _mm_loadl_epi64((__m128i *)(p + 0x10)));
                a1 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x20)), _mm_loadl_epi64((__m128i *)(p + 0x30)));
                a2 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x40)), _mm_loadl_epi64((__m128i *)(p + 0x50)));
                a3 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x60)), _mm_loadl_epi64((__m128i *)(p + 0x70)));
                a4 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x80)), _mm_loadl_epi64((__m128i *)(p + 0x90)));
                a5 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0xA0)), _mm_loadl_epi64((__m128i *)(p + 0xB0)));
                a6 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0xC0)), _mm_loadl_epi64((__m128i *)(p + 0xD0)));
                a7 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0xE0)), _mm_loadl_epi64((__m128i *)(p + 0xF0)));

                a8 = _mm_unpackhi_epi16(a0, a1);
                a0 = _mm_unpacklo_epi16(a0, a1);
                a9 = _mm_unpackhi_epi16(a2, a3);
                a2 = _mm_unpacklo_epi16(a2, a3);
                a10 = _mm_unpackhi_epi16(a4, a5);
                a4 = _mm_unpacklo_epi16(a4, a5);
                a11 = _mm_unpackhi_epi16(a6, a7);
                a6 = _mm_unpacklo_epi16(a6, a7);

                a1 = _mm_unpackhi_epi32(a0, a2);
                a0 = _mm_unpacklo_epi32(a0, a2);
                a3 = _mm_unpackhi_epi32(a4, a6);
                a4 = _mm_unpacklo_epi32(a4, a6);
                a5 = _mm_unpackhi_epi32(a8, a9);
                a8 = _mm_unpacklo_epi32(a8, a9);
                a7 = _mm_unpackhi_epi32(a10, a11);
                a10 = _mm_unpacklo_epi32(a10, a11);

                a2 = _mm_unpackhi_epi64(a0, a4);
                a0 = _mm_unpacklo_epi64(a0, a4);
                a6 = _mm_unpackhi_epi64(a8, a10);
                a8 = _mm_unpacklo_epi64(a8, a10);
                a9 = _mm_unpackhi_epi64(a1, a3);
                a1 = _mm_unpacklo_epi64(a1, a3);
                a11 = _mm_unpackhi_epi64(a5, a7);
                a5 = _mm_unpacklo_epi64(a5, a7);

                _mm_storeu_si128((__m128i *)prediction_ptr, a0);  prediction_ptr += prediction_buffer_stride;
                _mm_storeu_si128((__m128i *)prediction_ptr, a2);  prediction_ptr += prediction_buffer_stride;
                _mm_storeu_si128((__m128i *)prediction_ptr, a1);  prediction_ptr += prediction_buffer_stride;
                _mm_storeu_si128((__m128i *)prediction_ptr, a9);  prediction_ptr += prediction_buffer_stride;
                _mm_storeu_si128((__m128i *)prediction_ptr, a8);  prediction_ptr += prediction_buffer_stride;
                _mm_storeu_si128((__m128i *)prediction_ptr, a6);  prediction_ptr += prediction_buffer_stride;
                _mm_storeu_si128((__m128i *)prediction_ptr, a5);  prediction_ptr += prediction_buffer_stride;
                _mm_storeu_si128((__m128i *)prediction_ptr, a11); prediction_ptr += prediction_buffer_stride;
                p += 8;
            }
        }
        else { // size == 32
            for (colIndex = 0; colIndex < size; colIndex++) {
                deltaSum += intra_pred_angle;
                deltaInt = deltaSum >> 5;
                deltaFract = deltaSum & 31;
                a0 = _mm_unpacklo_epi8(_mm_cvtsi32_si128(32 - deltaFract), _mm_cvtsi32_si128(deltaFract));
                a0 = _mm_shuffle_epi8(a0, _mm_setr_epi8(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1));
                aa0 = _mm256_setr_m128i(a0, a0);
                ttop0 = _mm256_loadu_si256((__m256i *)(ref_samp_main + deltaInt));
                ttop1 = _mm256_loadu_si256((__m256i *)(ref_samp_main + deltaInt + 1));
                ttop2 = _mm256_unpacklo_epi8(ttop0, ttop1);
                ttop0 = _mm256_unpackhi_epi8(ttop0, ttop1);
                ssum0 = _mm256_add_epi16(_mm256_set1_epi16(16), _mm256_maddubs_epi16(ttop2, aa0));
                ssum1 = _mm256_add_epi16(_mm256_set1_epi16(16), _mm256_maddubs_epi16(ttop0, aa0));
                ssum0 = _mm256_srai_epi16(ssum0, 5);
                ssum1 = _mm256_srai_epi16(ssum1, 5);
                ssum0 = _mm256_packus_epi16(ssum0, ssum1);
                _mm256_storeu_si256((__m256i *)p, ssum0);
                p += 32;
            }
            p = temp_buf;
            for (colIndex = 0; colIndex < 2; colIndex++) {
                for (row_index = 0; row_index < 4; row_index++) {
                    a0 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x000)), _mm_loadl_epi64((__m128i *)(p + 0x020)));
                    a1 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x040)), _mm_loadl_epi64((__m128i *)(p + 0x060)));
                    a2 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x080)), _mm_loadl_epi64((__m128i *)(p + 0x0A0)));
                    a3 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x0C0)), _mm_loadl_epi64((__m128i *)(p + 0x0E0)));
                    a4 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x100)), _mm_loadl_epi64((__m128i *)(p + 0x120)));
                    a5 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x140)), _mm_loadl_epi64((__m128i *)(p + 0x160)));
                    a6 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x180)), _mm_loadl_epi64((__m128i *)(p + 0x1A0)));
                    a7 = _mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(p + 0x1C0)), _mm_loadl_epi64((__m128i *)(p + 0x1E0)));

                    a8 = _mm_unpackhi_epi16(a0, a1);
                    a0 = _mm_unpacklo_epi16(a0, a1);
                    a9 = _mm_unpackhi_epi16(a2, a3);
                    a2 = _mm_unpacklo_epi16(a2, a3);
                    a10 = _mm_unpackhi_epi16(a4, a5);
                    a4 = _mm_unpacklo_epi16(a4, a5);
                    a11 = _mm_unpackhi_epi16(a6, a7);
                    a6 = _mm_unpacklo_epi16(a6, a7);

                    a1 = _mm_unpackhi_epi32(a0, a2);
                    a0 = _mm_unpacklo_epi32(a0, a2);
                    a3 = _mm_unpackhi_epi32(a4, a6);
                    a4 = _mm_unpacklo_epi32(a4, a6);
                    a5 = _mm_unpackhi_epi32(a8, a9);
                    a8 = _mm_unpacklo_epi32(a8, a9);
                    a7 = _mm_unpackhi_epi32(a10, a11);
                    a10 = _mm_unpacklo_epi32(a10, a11);

                    a2 = _mm_unpackhi_epi64(a0, a4);
                    a0 = _mm_unpacklo_epi64(a0, a4);
                    a6 = _mm_unpackhi_epi64(a8, a10);
                    a8 = _mm_unpacklo_epi64(a8, a10);
                    a9 = _mm_unpackhi_epi64(a1, a3);
                    a1 = _mm_unpacklo_epi64(a1, a3);
                    a11 = _mm_unpackhi_epi64(a5, a7);
                    a5 = _mm_unpacklo_epi64(a5, a7);

                    _mm_storeu_si128((__m128i *)prediction_ptr, a0);  prediction_ptr += prediction_buffer_stride;
                    _mm_storeu_si128((__m128i *)prediction_ptr, a2);  prediction_ptr += prediction_buffer_stride;
                    _mm_storeu_si128((__m128i *)prediction_ptr, a1);  prediction_ptr += prediction_buffer_stride;
                    _mm_storeu_si128((__m128i *)prediction_ptr, a9);  prediction_ptr += prediction_buffer_stride;
                    _mm_storeu_si128((__m128i *)prediction_ptr, a8);  prediction_ptr += prediction_buffer_stride;
                    _mm_storeu_si128((__m128i *)prediction_ptr, a6);  prediction_ptr += prediction_buffer_stride;
                    _mm_storeu_si128((__m128i *)prediction_ptr, a5);  prediction_ptr += prediction_buffer_stride;
                    _mm_storeu_si128((__m128i *)prediction_ptr, a11); prediction_ptr += prediction_buffer_stride;
                    p += 8;
                }
                p = temp_buf + 0x200;
                prediction_ptr -= 32 * prediction_buffer_stride - 16;
            }
        }
    }
}

/*  ************************************************************  Vertical IntraPred***********************
                                                                   NOTE: this function has been updated only when size == 32******************************************/
void intra_mode_vertical_luma_avx2_intrin(
    const uint32_t      size,                   //input parameter, denotes the size of the current PU
    uint8_t            *ref_samples,             //input parameter, pointer to the reference samples
    uint8_t            *prediction_ptr,          //output parameter, pointer to the prediction
    const uint32_t      prediction_buffer_stride, //input parameter, denotes the stride for the prediction ptr
    const EbBool     skip                    //skip one row 
    )
{
    uint32_t topLeftOffset = size << 1;
    uint32_t leftOffset = 0;
    uint32_t topOffset = (size << 1) + 1;
    __m128i xmm0;
    uint32_t pStride = prediction_buffer_stride;

    if (size != 32) {
        __m128i xmm_mask1, xmm_mask2, xmm_topLeft, xmm_topLeft_lo, xmm_topLeft_hi, xmm_mask_skip, xmm_top, xmm_left, xmm_left_lo, xmm_left_hi;

        xmm0 = _mm_setzero_si128();
        xmm_mask1 = _mm_slli_si128(_mm_set1_epi8((signed char)0xFF), 1);
        xmm_mask2 = _mm_srli_si128(xmm_mask1, 15);

        xmm_topLeft = _mm_set_epi16((signed short)0xffff, (signed short)0xffff, (signed short)0xffff, (signed short)0xffff, (signed short)0xffff, (signed short)0xffff, (signed short)0xffff, *(uint16_t*)(ref_samples + topLeftOffset));
        xmm_topLeft = _mm_unpacklo_epi8(xmm_topLeft, xmm0);
        xmm_topLeft = _mm_unpacklo_epi16(xmm_topLeft, xmm_topLeft);
        xmm_topLeft = _mm_unpacklo_epi32(xmm_topLeft, xmm_topLeft);
        xmm_topLeft_hi = _mm_unpackhi_epi64(xmm_topLeft, xmm_topLeft);
        xmm_topLeft_lo = _mm_unpacklo_epi64(xmm_topLeft, xmm_topLeft);

        if (!skip) {

            if (size == 8) {

                xmm_left = _mm_packus_epi16(_mm_add_epi16(_mm_srai_epi16(_mm_sub_epi16(_mm_unpacklo_epi8(_mm_loadl_epi64((__m128i *)(ref_samples + leftOffset)), xmm0), xmm_topLeft_lo), 1), xmm_topLeft_hi), xmm0);
                xmm_top = _mm_and_si128(_mm_loadl_epi64((__m128i *)(ref_samples + topOffset)), xmm_mask1);
                MACRO_VERTICAL_LUMA_8(xmm_left, xmm_mask2, xmm_top)
                    prediction_ptr += (pStride << 2);
                MACRO_VERTICAL_LUMA_8(xmm_left, xmm_mask2, xmm_top)
            }
            else if (size == 16) {
                xmm_left = _mm_loadu_si128((__m128i *)(ref_samples + leftOffset));
                xmm_left_lo = _mm_add_epi16(_mm_srai_epi16(_mm_sub_epi16(_mm_unpacklo_epi8(xmm_left, xmm0), xmm_topLeft_lo), 1), xmm_topLeft_hi);
                xmm_left_hi = _mm_add_epi16(_mm_srai_epi16(_mm_sub_epi16(_mm_unpackhi_epi8(xmm_left, xmm0), xmm_topLeft_lo), 1), xmm_topLeft_hi);
                xmm_left = _mm_packus_epi16(xmm_left_lo, xmm_left_hi);
                xmm_top = _mm_and_si128(_mm_loadu_si128((__m128i *)(ref_samples + topOffset)), xmm_mask1);
                MACRO_VERTICAL_LUMA_16(xmm_left, xmm_mask2, xmm_top)
                    prediction_ptr += (pStride << 2);
                MACRO_VERTICAL_LUMA_16(xmm_left, xmm_mask2, xmm_top)
                    prediction_ptr += (pStride << 2);
                MACRO_VERTICAL_LUMA_16(xmm_left, xmm_mask2, xmm_top)
                    prediction_ptr += (pStride << 2);
                MACRO_VERTICAL_LUMA_16(xmm_left, xmm_mask2, xmm_top)
            }
            else {
                xmm_left = _mm_packus_epi16(_mm_add_epi16(_mm_srai_epi16(_mm_sub_epi16(_mm_unpacklo_epi8(_mm_cvtsi32_si128(*(uint32_t*)(ref_samples + leftOffset)), xmm0), xmm_topLeft_lo), 1), xmm_topLeft_hi), xmm0);
                xmm_top = _mm_and_si128(_mm_loadl_epi64((__m128i *)(ref_samples + topOffset)), xmm_mask1);
                MACRO_VERTICAL_LUMA_4(xmm_left, xmm_mask2, xmm_top)
                    prediction_ptr += (pStride << 1);
                MACRO_VERTICAL_LUMA_4(xmm_left, xmm_mask2, xmm_top)
            }
        }
        else {
            pStride <<= 1;
            xmm_mask_skip = _mm_set1_epi16(0x00FF);

            if (size == 8) {
                xmm_left = _mm_packus_epi16(_mm_add_epi16(_mm_srai_epi16(_mm_sub_epi16(_mm_and_si128(_mm_loadl_epi64((__m128i *)(ref_samples + leftOffset)), xmm_mask_skip), xmm_topLeft_lo), 1), xmm_topLeft_hi), xmm0);
                xmm_top = _mm_and_si128(_mm_loadl_epi64((__m128i *)(ref_samples + topOffset)), xmm_mask1);
                MACRO_VERTICAL_LUMA_8(xmm_left, xmm_mask2, xmm_top);
            }
            else if (size == 16) {
                xmm_left = _mm_packus_epi16(_mm_add_epi16(_mm_srai_epi16(_mm_sub_epi16(_mm_and_si128(_mm_loadu_si128((__m128i *)(ref_samples + leftOffset)), xmm_mask_skip), xmm_topLeft_lo), 1), xmm_topLeft_hi), xmm0);
                xmm_top = _mm_and_si128(_mm_loadu_si128((__m128i *)(ref_samples + topOffset)), xmm_mask1);
                MACRO_VERTICAL_LUMA_16(xmm_left, xmm_mask2, xmm_top)
                    prediction_ptr += (pStride << 2);
                MACRO_VERTICAL_LUMA_16(xmm_left, xmm_mask2, xmm_top)
            }
            else {
                xmm_left = _mm_packus_epi16(_mm_add_epi16(_mm_srai_epi16(_mm_sub_epi16(_mm_and_si128(_mm_cvtsi32_si128(*(uint32_t*)(ref_samples + leftOffset)), xmm_mask_skip), xmm_topLeft_lo), 1), xmm_topLeft_hi), xmm0);
                xmm_top = _mm_and_si128(_mm_cvtsi32_si128(*(uint32_t*)(ref_samples + topOffset)), xmm_mask1);
                MACRO_VERTICAL_LUMA_4(xmm_left, xmm_mask2, xmm_top)
            }
        }
    }
    else {
        __m256i xmm0;
        uint64_t size_to_write;
        uint32_t count;

        // Each storeu calls stores 32 bytes. Hence each iteration stores 8 * 32 bytes.
        // Depending on skip, we need 4 or 2 iterations to store 32x32 bytes.
        size_to_write = 4 >> (skip ? 1 : 0);
        pStride = pStride << (skip ? 1 : 0);

        xmm0 = _mm256_loadu_si256((__m256i *)(ref_samples + topOffset));

        for (count = 0; count < size_to_write; count++) {
            _mm256_storeu_si256((__m256i *)(prediction_ptr), xmm0);
            _mm256_storeu_si256((__m256i *)(prediction_ptr + pStride), xmm0);
            _mm256_storeu_si256((__m256i *)(prediction_ptr + 2 * pStride), xmm0);
            _mm256_storeu_si256((__m256i *)(prediction_ptr + 3 * pStride), xmm0);

            prediction_ptr += (pStride << 2);
            _mm256_storeu_si256((__m256i *)(prediction_ptr), xmm0);
            _mm256_storeu_si256((__m256i *)(prediction_ptr + pStride), xmm0);
            _mm256_storeu_si256((__m256i *)(prediction_ptr + 2 * pStride), xmm0);
            _mm256_storeu_si256((__m256i *)(prediction_ptr + 3 * pStride), xmm0);

            prediction_ptr += (pStride << 2);
        }
    }

    return;
}

void intra_mode_dc_luma_avx2_intrin(
    const uint32_t      size,                       //input parameter, denotes the size of the current PU
    uint8_t            *ref_samples,                 //input parameter, pointer to the reference samples
    uint8_t            *prediction_ptr,              //output parameter, pointer to the prediction
    const uint32_t      prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
    const EbBool     skip)                       //skip one row 
{
    __m128i xmm0 = _mm_setzero_si128();
    uint32_t pStride = prediction_buffer_stride;
    uint32_t topOffset = (size << 1) + 1;
    uint32_t leftOffset = 0;

    if (size != 32) {

        __m128i xmm_mask1 = _mm_slli_si128(_mm_set1_epi8((signed char)0xFF), 1);
        __m128i xmm_mask2 = _mm_srli_si128(xmm_mask1, 15);
        __m128i xmm_C2 = _mm_set1_epi16(0x0002);

        if (!skip) {

            if (size == 16) {
                __m128i xmm_predictionDcValue, xmm_top, xmm_left, xmm_sum, xmm_prediction_ptr_0;
                __m128i xmm_top_lo, xmm_top_hi, xmm_left_lo, xmm_left_hi, xmm_predictionDcValue_16, xmm_predictionDcValue_16_x2, xmm_predictionDcValue_16_x3;
                xmm_top = _mm_loadu_si128((__m128i *)(ref_samples + topOffset));
                xmm_left = _mm_loadu_si128((__m128i *)(ref_samples + leftOffset));
                xmm_top_lo = _mm_unpacklo_epi8(xmm_top, xmm0);
                xmm_top_hi = _mm_unpackhi_epi8(xmm_top, xmm0);
                xmm_left_lo = _mm_unpacklo_epi8(xmm_left, xmm0);
                xmm_left_hi = _mm_unpackhi_epi8(xmm_left, xmm0);

                xmm_sum = _mm_add_epi32(_mm_sad_epu8(xmm_top, xmm0), _mm_sad_epu8(xmm_left, xmm0));

                xmm_predictionDcValue = _mm_srli_epi32(_mm_add_epi32(_mm_add_epi32(_mm_srli_si128(xmm_sum, 8), xmm_sum), _mm_cvtsi32_si128(16)), 5);
                xmm_predictionDcValue = _mm_unpacklo_epi8(xmm_predictionDcValue, xmm_predictionDcValue);
                xmm_predictionDcValue = _mm_unpacklo_epi16(xmm_predictionDcValue, xmm_predictionDcValue);
                xmm_predictionDcValue = _mm_unpacklo_epi32(xmm_predictionDcValue, xmm_predictionDcValue);
                xmm_predictionDcValue = _mm_unpacklo_epi64(xmm_predictionDcValue, xmm_predictionDcValue);

                xmm_predictionDcValue_16 = _mm_srli_epi16(xmm_predictionDcValue, 8);
                xmm_predictionDcValue = _mm_and_si128(xmm_predictionDcValue, xmm_mask1);
                xmm_predictionDcValue_16_x2 = _mm_add_epi16(xmm_predictionDcValue_16, xmm_predictionDcValue_16);
                xmm_predictionDcValue_16_x3 = _mm_add_epi16(xmm_predictionDcValue_16_x2, xmm_predictionDcValue_16);

                xmm_top = _mm_packus_epi16(_mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(xmm_top_lo, xmm_predictionDcValue_16_x3), xmm_C2), 2),
                    _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(xmm_top_hi, xmm_predictionDcValue_16_x3), xmm_C2), 2));

                xmm_left = _mm_srli_si128(_mm_packus_epi16(_mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(xmm_left_lo, xmm_predictionDcValue_16_x3), xmm_C2), 2),
                    _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(xmm_left_hi, xmm_predictionDcValue_16_x3), xmm_C2), 2)), 1);

                xmm_prediction_ptr_0 = _mm_or_si128(_mm_and_si128(_mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(_mm_add_epi16(xmm_top_lo, xmm_left_lo), xmm_predictionDcValue_16_x2), xmm_C2), 2), xmm_mask2), _mm_and_si128(xmm_top, xmm_mask1));
                _mm_storeu_si128((__m128i *)(prediction_ptr), xmm_prediction_ptr_0);

                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue));
                xmm_left = _mm_srli_si128(xmm_left, 1);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue));
                xmm_left = _mm_srli_si128(xmm_left, 1);
                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue));
                xmm_left = _mm_srli_si128(xmm_left, 1);

                prediction_ptr += (pStride << 2);
                MACRO_VERTICAL_LUMA_16(xmm_left, xmm_mask2, xmm_predictionDcValue)
                    prediction_ptr += (pStride << 2);
                MACRO_VERTICAL_LUMA_16(xmm_left, xmm_mask2, xmm_predictionDcValue)
                    prediction_ptr += (pStride << 2);
                MACRO_VERTICAL_LUMA_16(xmm_left, xmm_mask2, xmm_predictionDcValue)
            }
            else if (size == 8) {

                __m128i xmm_left, xmm_top, xmm_top_lo, xmm_left_lo, xmm_predictionDcValue, xmm_predictionDcValue_16;
                __m128i xmm_predictionDcValue_16_x2, xmm_predictionDcValue_16_x3, xmm_prediction_ptr_0;

                xmm_top = _mm_loadl_epi64((__m128i *)(ref_samples + topOffset));
                xmm_left = _mm_loadl_epi64((__m128i *)(ref_samples + leftOffset));

                xmm_top_lo = _mm_unpacklo_epi8(xmm_top, xmm0);
                xmm_left_lo = _mm_unpacklo_epi8(xmm_left, xmm0);

                xmm_predictionDcValue = _mm_srli_epi32(_mm_add_epi32(_mm_add_epi32(_mm_sad_epu8(xmm_top, xmm0), _mm_sad_epu8(xmm_left, xmm0)), _mm_cvtsi32_si128(8)), 4);
                xmm_predictionDcValue = _mm_unpacklo_epi8(xmm_predictionDcValue, xmm_predictionDcValue);
                xmm_predictionDcValue = _mm_unpacklo_epi16(xmm_predictionDcValue, xmm_predictionDcValue);
                xmm_predictionDcValue = _mm_unpacklo_epi32(xmm_predictionDcValue, xmm_predictionDcValue);
                xmm_predictionDcValue = _mm_unpacklo_epi64(xmm_predictionDcValue, xmm_predictionDcValue);

                xmm_predictionDcValue_16 = _mm_srli_epi16(xmm_predictionDcValue, 8);
                xmm_predictionDcValue = _mm_and_si128(xmm_predictionDcValue, xmm_mask1);

                xmm_predictionDcValue_16_x2 = _mm_add_epi16(xmm_predictionDcValue_16, xmm_predictionDcValue_16);
                xmm_predictionDcValue_16_x3 = _mm_add_epi16(xmm_predictionDcValue_16_x2, xmm_predictionDcValue_16);

                xmm_top = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(xmm_top_lo, xmm_predictionDcValue_16_x3), xmm_C2), 2);
                xmm_left = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(xmm_left_lo, xmm_predictionDcValue_16_x3), xmm_C2), 2);

                xmm_left = _mm_srli_si128(_mm_packus_epi16(xmm_left, xmm_left), 1);
                xmm_top = _mm_and_si128(_mm_packus_epi16(xmm_top, xmm_top), xmm_mask1);

                xmm_prediction_ptr_0 = _mm_or_si128(_mm_and_si128(_mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(_mm_add_epi16(xmm_top_lo, xmm_left_lo), xmm_predictionDcValue_16_x2), xmm_C2), 2), xmm_mask2), xmm_top);
                _mm_storel_epi64((__m128i *)(prediction_ptr), xmm_prediction_ptr_0);

                _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue));
                xmm_left = _mm_srli_si128(xmm_left, 1);
                _mm_storel_epi64((__m128i *)(prediction_ptr + 2 * pStride), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue));
                xmm_left = _mm_srli_si128(xmm_left, 1);
                _mm_storel_epi64((__m128i *)(prediction_ptr + 3 * pStride), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue));
                xmm_left = _mm_srli_si128(xmm_left, 1);

                prediction_ptr += (pStride << 2);
                MACRO_VERTICAL_LUMA_8(xmm_left, xmm_mask2, xmm_predictionDcValue)

            }
            else {
                __m128i xmm_left, xmm_top, xmm_top_lo, xmm_left_lo, xmm_predictionDcValue, xmm_predictionDcValue_16;
                __m128i xmm_predictionDcValue_16_x2, xmm_predictionDcValue_16_x3, xmm_prediction_ptr_0;

                xmm_top = _mm_cvtsi32_si128(*(uint32_t*)(ref_samples + topOffset));
                xmm_left = _mm_cvtsi32_si128(*(uint32_t*)(ref_samples + leftOffset));

                xmm_top_lo = _mm_unpacklo_epi8(xmm_top, xmm0);
                xmm_left_lo = _mm_unpacklo_epi8(xmm_left, xmm0);
                xmm_predictionDcValue = _mm_srli_epi32(_mm_add_epi32(_mm_add_epi32(_mm_sad_epu8(xmm_top, xmm0), _mm_sad_epu8(xmm_left, xmm0)), _mm_cvtsi32_si128(4)), 3);
                xmm_predictionDcValue = _mm_unpacklo_epi8(xmm_predictionDcValue, xmm_predictionDcValue);
                xmm_predictionDcValue = _mm_unpacklo_epi16(xmm_predictionDcValue, xmm_predictionDcValue);
                xmm_predictionDcValue = _mm_unpacklo_epi32(xmm_predictionDcValue, xmm_predictionDcValue);

                xmm_predictionDcValue_16 = _mm_srli_epi16(xmm_predictionDcValue, 8);
                xmm_predictionDcValue = _mm_and_si128(xmm_predictionDcValue, xmm_mask1);

                xmm_predictionDcValue_16_x2 = _mm_add_epi16(xmm_predictionDcValue_16, xmm_predictionDcValue_16);
                xmm_predictionDcValue_16_x3 = _mm_add_epi16(xmm_predictionDcValue_16_x2, xmm_predictionDcValue_16);

                xmm_top = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(xmm_top_lo, xmm_predictionDcValue_16_x3), xmm_C2), 2);
                xmm_left = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(xmm_left_lo, xmm_predictionDcValue_16_x3), xmm_C2), 2);
                xmm_top = _mm_and_si128(_mm_packus_epi16(xmm_top, xmm_top), xmm_mask1);
                xmm_left = _mm_srli_si128(_mm_packus_epi16(xmm_left, xmm_left), 1);

                xmm_prediction_ptr_0 = _mm_or_si128(_mm_and_si128(_mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(_mm_add_epi16(xmm_top_lo, xmm_left_lo), xmm_predictionDcValue_16_x2), xmm_C2), 2), xmm_mask2), xmm_top);
                *(uint32_t*)prediction_ptr = _mm_cvtsi128_si32(xmm_prediction_ptr_0);

                *(uint32_t*)(prediction_ptr + pStride) = _mm_cvtsi128_si32(_mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue));
                xmm_left = _mm_srli_si128(xmm_left, 1);

                *(uint32_t*)(prediction_ptr + 2 * pStride) = _mm_cvtsi128_si32(_mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue));
                *(uint32_t*)(prediction_ptr + 3 * pStride) = _mm_cvtsi128_si32(_mm_or_si128(_mm_and_si128(_mm_srli_si128(xmm_left, 1), xmm_mask2), xmm_predictionDcValue));
            }
        }
        else {
            pStride <<= 1;

            __m128i xmm_skip_mask = _mm_set1_epi16(0x00FF);
            __m128i xmm_left, xmm_sum, xmm_top, xmm_top_lo, xmm_top_hi, xmm_left_skipped, xmm_predictionDcValue, xmm_predictionDcValue_16;
            __m128i xmm_predictionDcValue_16_x2, xmm_predictionDcValue_16_x3, xmm_prediction_ptr_0;

            if (size == 16) {
                xmm_top = _mm_loadu_si128((__m128i *)(ref_samples + topOffset));
                xmm_left = _mm_loadu_si128((__m128i *)(ref_samples + leftOffset));

                xmm_top_lo = _mm_unpacklo_epi8(xmm_top, xmm0);
                xmm_top_hi = _mm_unpackhi_epi8(xmm_top, xmm0);

                xmm_left_skipped = _mm_and_si128(xmm_skip_mask, xmm_left);

                xmm_sum = _mm_add_epi32(_mm_sad_epu8(xmm_top, xmm0), _mm_sad_epu8(xmm_left, xmm0));

                xmm_predictionDcValue = _mm_srli_epi32(_mm_add_epi32(_mm_add_epi32(_mm_srli_si128(xmm_sum, 8), xmm_sum), _mm_cvtsi32_si128(16)), 5);
                xmm_predictionDcValue = _mm_unpacklo_epi8(xmm_predictionDcValue, xmm_predictionDcValue);
                xmm_predictionDcValue = _mm_unpacklo_epi16(xmm_predictionDcValue, xmm_predictionDcValue);
                xmm_predictionDcValue = _mm_unpacklo_epi32(xmm_predictionDcValue, xmm_predictionDcValue);
                xmm_predictionDcValue = _mm_unpacklo_epi64(xmm_predictionDcValue, xmm_predictionDcValue);

                xmm_predictionDcValue_16 = _mm_srli_epi16(xmm_predictionDcValue, 8);
                xmm_predictionDcValue = _mm_and_si128(xmm_predictionDcValue, xmm_mask1);

                xmm_predictionDcValue_16_x2 = _mm_add_epi16(xmm_predictionDcValue_16, xmm_predictionDcValue_16);
                xmm_predictionDcValue_16_x3 = _mm_add_epi16(xmm_predictionDcValue_16_x2, xmm_predictionDcValue_16);

                xmm_left = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(xmm_left_skipped, xmm_predictionDcValue_16_x3), xmm_C2), 2);
                xmm_left = _mm_srli_si128(_mm_packus_epi16(xmm_left, xmm_left), 1);

                xmm_top = _mm_and_si128(_mm_packus_epi16(_mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(xmm_top_lo, xmm_predictionDcValue_16_x3), xmm_C2), 2),
                    _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(xmm_top_hi, xmm_predictionDcValue_16_x3), xmm_C2), 2)), xmm_mask1);

                xmm_prediction_ptr_0 = _mm_or_si128(_mm_and_si128(_mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(_mm_add_epi16(xmm_top_lo, xmm_left_skipped), xmm_predictionDcValue_16_x2), xmm_C2), 2), xmm_mask2), xmm_top);
                _mm_storeu_si128((__m128i *)prediction_ptr, xmm_prediction_ptr_0);

                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue));
                xmm_left = _mm_srli_si128(xmm_left, 1);

                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue));
                xmm_left = _mm_srli_si128(xmm_left, 1);

                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue));
                xmm_left = _mm_srli_si128(xmm_left, 1);
                prediction_ptr += (pStride << 2);
                MACRO_VERTICAL_LUMA_16(xmm_left, xmm_mask2, xmm_predictionDcValue)
            }
            else {

                xmm_top = _mm_loadl_epi64((__m128i *)(ref_samples + topOffset));
                xmm_left = _mm_loadl_epi64((__m128i *)(ref_samples + leftOffset));

                xmm_top_lo = _mm_unpacklo_epi8(xmm_top, xmm0);
                xmm_left_skipped = _mm_and_si128(xmm_skip_mask, xmm_left);

                xmm_predictionDcValue = _mm_srli_epi32(_mm_add_epi32(_mm_add_epi32(_mm_sad_epu8(xmm_top, xmm0), _mm_sad_epu8(xmm_left, xmm0)), _mm_cvtsi32_si128(8)), 4);
                xmm_predictionDcValue = _mm_unpacklo_epi8(xmm_predictionDcValue, xmm_predictionDcValue);
                xmm_predictionDcValue = _mm_unpacklo_epi16(xmm_predictionDcValue, xmm_predictionDcValue);
                xmm_predictionDcValue = _mm_unpacklo_epi32(xmm_predictionDcValue, xmm_predictionDcValue);
                xmm_predictionDcValue = _mm_unpacklo_epi64(xmm_predictionDcValue, xmm_predictionDcValue);

                xmm_predictionDcValue_16 = _mm_srli_epi16(xmm_predictionDcValue, 8);
                xmm_predictionDcValue = _mm_and_si128(xmm_predictionDcValue, xmm_mask1);

                xmm_predictionDcValue_16_x2 = _mm_add_epi16(xmm_predictionDcValue_16, xmm_predictionDcValue_16);
                xmm_predictionDcValue_16_x3 = _mm_add_epi16(xmm_predictionDcValue_16_x2, xmm_predictionDcValue_16);

                xmm_top = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(xmm_top_lo, xmm_predictionDcValue_16_x3), xmm_C2), 2);
                xmm_left = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(xmm_left_skipped, xmm_predictionDcValue_16_x3), xmm_C2), 2);
                xmm_prediction_ptr_0 = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(_mm_add_epi16(xmm_top_lo, xmm_left_skipped), xmm_predictionDcValue_16_x2), xmm_C2), 2);

                xmm_left = _mm_srli_si128(_mm_packus_epi16(xmm_left, xmm_left), 1);
                xmm_prediction_ptr_0 = _mm_or_si128(_mm_and_si128(xmm_prediction_ptr_0, xmm_mask2), _mm_and_si128(_mm_packus_epi16(xmm_top, xmm_top), xmm_mask1));
                _mm_storel_epi64((__m128i *)prediction_ptr, xmm_prediction_ptr_0);

                _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue));
                xmm_left = _mm_srli_si128(xmm_left, 1);

                _mm_storel_epi64((__m128i *)(prediction_ptr + 2 * pStride), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue));

                _mm_storel_epi64((__m128i *)(prediction_ptr + 3 * pStride), _mm_or_si128(_mm_and_si128(_mm_srli_si128(xmm_left, 1), xmm_mask2), xmm_predictionDcValue));

            }
        }
    }



    /*************************************************************************************************************************************************************************************************************/
    else {
        __m256i xmm_sum, xmm_sadleft, xmm_sadtop, xmm_toptmp, xmm_lefttmp, xmm_set, xmm_sum128_2, xmm_sum256, xmm_predictionDcValue;
        __m256i xmm1 = _mm256_setzero_si256();
        __m128i xmm_sumhi, xmm_sumlo, xmm_sum1, xmm_sum128, xmm_sumhitmp, xmm_sumlotmp, xmm_movelotmp, xmm_movehitmp;

        xmm_sumhi = xmm_sumlo = xmm_sum128 = xmm_sumhitmp = xmm_sumlotmp = _mm_setzero_si128();
        xmm_sum = xmm_sadleft = xmm_sadtop = xmm_toptmp = xmm_lefttmp = _mm256_setzero_si256();

        xmm_toptmp = _mm256_sad_epu8(_mm256_set_m128i(_mm_loadu_si128((__m128i *)(ref_samples + topOffset + 16)), _mm_loadu_si128((__m128i *)(ref_samples + topOffset))), xmm1);
        xmm_lefttmp = _mm256_sad_epu8(_mm256_set_m128i(_mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 16)), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset))), xmm1);

        xmm_sum = _mm256_add_epi32(xmm_toptmp, xmm_lefttmp);
        xmm_sum = _mm256_hadd_epi32(xmm_sum, xmm_sum);
        xmm_sumlo = _mm256_extracti128_si256(xmm_sum, 0);
        xmm_sumhi = _mm256_extracti128_si256(xmm_sum, 1);

        xmm_movelotmp = _mm_move_epi64(xmm_sumlo);
        xmm_movehitmp = _mm_move_epi64(xmm_sumhi);

        xmm_sum1 = _mm_add_epi32(xmm_movelotmp, xmm_movehitmp);

        xmm_sum1 = _mm_hadd_epi32(xmm_sum1, xmm_sum1);

        xmm_sum256 = _mm256_castsi128_si256(xmm_sum1);


        //#endif      
        xmm_set = _mm256_castsi128_si256(_mm_set1_epi32(32));

        xmm_sum128_2 = _mm256_add_epi32(xmm_sum256, xmm_set); // add offset
        xmm_predictionDcValue = _mm256_srli_epi32(xmm_sum128_2, 6); //_mm256_srli_epi32


        __m128i dc128 = _mm256_castsi256_si128(xmm_predictionDcValue);

        uint8_t dc = _mm_cvtsi128_si32(dc128);
        xmm_predictionDcValue = _mm256_set1_epi8(dc);//_mm_broadcastb_epi8


        uint32_t count;

        for (count = 0; count < 2; ++count) {

            _mm256_storeu_si256((__m256i *) prediction_ptr, xmm_predictionDcValue);
            _mm256_storeu_si256((__m256i *)(prediction_ptr + 1 * pStride), xmm_predictionDcValue);
            _mm256_storeu_si256((__m256i *)(prediction_ptr + 2 * pStride), xmm_predictionDcValue);
            _mm256_storeu_si256((__m256i *)(prediction_ptr + 3 * pStride), xmm_predictionDcValue);

            prediction_ptr += (pStride << 2);

            _mm256_storeu_si256((__m256i *) prediction_ptr, xmm_predictionDcValue);
            _mm256_storeu_si256((__m256i *)(prediction_ptr + 1 * pStride), xmm_predictionDcValue);
            _mm256_storeu_si256((__m256i *)(prediction_ptr + 2 * pStride), xmm_predictionDcValue);
            _mm256_storeu_si256((__m256i *)(prediction_ptr + 3 * pStride), xmm_predictionDcValue);

            prediction_ptr += (pStride << 2);

            _mm256_storeu_si256((__m256i *) prediction_ptr, xmm_predictionDcValue);
            _mm256_storeu_si256((__m256i *)(prediction_ptr + 1 * pStride), xmm_predictionDcValue);
            _mm256_storeu_si256((__m256i *)(prediction_ptr + 2 * pStride), xmm_predictionDcValue);
            _mm256_storeu_si256((__m256i *)(prediction_ptr + 3 * pStride), xmm_predictionDcValue);

            prediction_ptr += (pStride << 2);

            _mm256_storeu_si256((__m256i *) prediction_ptr, xmm_predictionDcValue);
            _mm256_storeu_si256((__m256i *)(prediction_ptr + 1 * pStride), xmm_predictionDcValue);
            _mm256_storeu_si256((__m256i *)(prediction_ptr + 2 * pStride), xmm_predictionDcValue);
            _mm256_storeu_si256((__m256i *)(prediction_ptr + 3 * pStride), xmm_predictionDcValue);

            prediction_ptr += (pStride << 2);


        }


    }


}

// TODO(luoyi) The following two functions are shared with intrapred_sse2.c.
// Use a header file, intrapred_common_x86.h
static INLINE __m128i dc_sum_16_sse2(const uint8_t *ref) {
    __m128i x = _mm_load_si128((__m128i const *)ref);
    const __m128i zero = _mm_setzero_si128();
    x = _mm_sad_epu8(x, zero);
    const __m128i high = _mm_unpackhi_epi64(x, x);
    return _mm_add_epi16(x, high);
}

static INLINE __m128i dc_sum_32_sse2(const uint8_t *ref) {
    __m128i x0 = _mm_load_si128((__m128i const *)ref);
    __m128i x1 = _mm_load_si128((__m128i const *)(ref + 16));
    const __m128i zero = _mm_setzero_si128();
    x0 = _mm_sad_epu8(x0, zero);
    x1 = _mm_sad_epu8(x1, zero);
    x0 = _mm_add_epi16(x0, x1);
    const __m128i high = _mm_unpackhi_epi64(x0, x0);
    return _mm_add_epi16(x0, high);
}

static INLINE __m256i dc_sum_32(const uint8_t *ref) {
    const __m256i x = _mm256_loadu_si256((const __m256i *)ref);
    const __m256i zero = _mm256_setzero_si256();
    __m256i y = _mm256_sad_epu8(x, zero);
    __m256i u = _mm256_permute2x128_si256(y, y, 1);
    y = _mm256_add_epi64(u, y);
    u = _mm256_unpackhi_epi64(y, y);
    return _mm256_add_epi16(y, u);
}
static INLINE void row_store_32xh(const __m256i *r, int32_t height, uint8_t *dst,
    ptrdiff_t stride) {
    for (int32_t i = 0; i < height; ++i) {
        _mm256_storeu_si256((__m256i *)dst, *r);
        dst += stride;
    }
}
void intra_mode_dc_32x32_av1_avx2_intrin(
    EbBool        isLeftAvailble,
    EbBool        isAboveAvailble,
    const uint32_t   size,                       //input parameter, denotes the size of the current PU
    uint8_t         *ref_samples,                 //input parameter, pointer to the reference samples
    uint8_t         *dst,              //output parameter, pointer to the prediction
    const uint32_t   prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
    const EbBool  skip)                       //skip half rows
{

    uint32_t leftOffset = 0;
    uint32_t topOffset = (size << 1) + 1;
    uint32_t rowStride = skip ? 2 : 1;

    if (isLeftAvailble && !isAboveAvailble) {
        __m256i sum = dc_sum_32(&ref_samples[leftOffset]);
        const __m256i sixteen = _mm256_set1_epi16(16);
        sum = _mm256_add_epi16(sum, sixteen);
        sum = _mm256_srai_epi16(sum, 5);
        const __m256i zero = _mm256_setzero_si256();
        __m256i row = _mm256_shuffle_epi8(sum, zero);
        row_store_32xh(&row, 32, dst, rowStride * prediction_buffer_stride);
    }
    else if (isAboveAvailble && !isLeftAvailble) {
        __m256i sum = dc_sum_32(&ref_samples[topOffset]);
        const __m256i sixteen = _mm256_set1_epi16(16);
        sum = _mm256_add_epi16(sum, sixteen);
        sum = _mm256_srai_epi16(sum, 5);
        const __m256i zero = _mm256_setzero_si256();
        __m256i row = _mm256_shuffle_epi8(sum, zero);
        row_store_32xh(&row, 32, dst, rowStride * prediction_buffer_stride);
    }
    else {
        const __m256i sum_above = dc_sum_32(&ref_samples[topOffset]);
        __m256i sum_left = dc_sum_32(&ref_samples[leftOffset]);
        sum_left = _mm256_add_epi16(sum_left, sum_above);
        const __m256i thirtytwo = _mm256_set1_epi16(32);
        sum_left = _mm256_add_epi16(sum_left, thirtytwo);
        sum_left = _mm256_srai_epi16(sum_left, 6);
        const __m256i zero = _mm256_setzero_si256();
        __m256i row = _mm256_shuffle_epi8(sum_left, zero);
        row_store_32xh(&row, 32, dst, rowStride * prediction_buffer_stride);
    }
}

static INLINE void row_store_64xh(const __m256i *r, int32_t height, uint8_t *dst,
    ptrdiff_t stride) {
    for (int32_t i = 0; i < height; ++i) {
        _mm256_storeu_si256((__m256i *)dst, *r);
        _mm256_storeu_si256((__m256i *)(dst + 32), *r);
        dst += stride;
    }
}
static INLINE __m256i dc_sum_64(const uint8_t *ref) {
    const __m256i x0 = _mm256_loadu_si256((const __m256i *)ref);
    const __m256i x1 = _mm256_loadu_si256((const __m256i *)(ref + 32));
    const __m256i zero = _mm256_setzero_si256();
    __m256i y0 = _mm256_sad_epu8(x0, zero);
    __m256i y1 = _mm256_sad_epu8(x1, zero);
    y0 = _mm256_add_epi64(y0, y1);
    __m256i u0 = _mm256_permute2x128_si256(y0, y0, 1);
    y0 = _mm256_add_epi64(u0, y0);
    u0 = _mm256_unpackhi_epi64(y0, y0);
    return _mm256_add_epi16(y0, u0);
}
void aom_dc_predictor_64x64_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above, const uint8_t *left) {
    const __m256i sum_above = dc_sum_64(above);
    __m256i sum_left = dc_sum_64(left);
    sum_left = _mm256_add_epi16(sum_left, sum_above);
    uint32_t sum = _mm_cvtsi128_si32(_mm256_castsi256_si128(sum_left));
    sum += 64;
    sum /= 128;
    const __m256i row = _mm256_set1_epi8((uint8_t)sum);
    row_store_64xh(&row, 64, dst, stride);
}

void aom_dc_left_predictor_64x64_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above,
    const uint8_t *left) {
    __m256i sum = dc_sum_64(left);
    (void)above;

    const __m256i thirtytwo = _mm256_set1_epi16(32);
    sum = _mm256_add_epi16(sum, thirtytwo);
    sum = _mm256_srai_epi16(sum, 6);
    const __m256i zero = _mm256_setzero_si256();
    __m256i row = _mm256_shuffle_epi8(sum, zero);
    row_store_64xh(&row, 64, dst, stride);
}
void aom_dc_top_predictor_64x64_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above,
    const uint8_t *left) {
    __m256i sum = dc_sum_64(above);
    (void)left;

    const __m256i thirtytwo = _mm256_set1_epi16(32);
    sum = _mm256_add_epi16(sum, thirtytwo);
    sum = _mm256_srai_epi16(sum, 6);
    const __m256i zero = _mm256_setzero_si256();
    __m256i row = _mm256_shuffle_epi8(sum, zero);
    row_store_64xh(&row, 64, dst, stride);
}
#if INTRA_ASM
void aom_dc_top_predictor_32x32_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above,
    const uint8_t *left) {
    __m256i sum = dc_sum_32(above);
    (void)left;

    const __m256i sixteen = _mm256_set1_epi16(16);
    sum = _mm256_add_epi16(sum, sixteen);
    sum = _mm256_srai_epi16(sum, 5);
    const __m256i zero = _mm256_setzero_si256();
    __m256i row = _mm256_shuffle_epi8(sum, zero);
    row_store_32xh(&row, 32, dst, stride);
}
void aom_dc_left_predictor_32x32_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above,
    const uint8_t *left) {
    __m256i sum = dc_sum_32(left);
    (void)above;

    const __m256i sixteen = _mm256_set1_epi16(16);
    sum = _mm256_add_epi16(sum, sixteen);
    sum = _mm256_srai_epi16(sum, 5);
    const __m256i zero = _mm256_setzero_si256();
    __m256i row = _mm256_shuffle_epi8(sum, zero);
    row_store_32xh(&row, 32, dst, stride);
}
void aom_dc_128_predictor_64x64_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above,
    const uint8_t *left) {
    (void)above;
    (void)left;
    const __m256i row = _mm256_set1_epi8((uint8_t)0x80);
    row_store_64xh(&row, 64, dst, stride);
}
void aom_dc_128_predictor_32x32_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above,
    const uint8_t *left) {
    (void)above;
    (void)left;
    const __m256i row = _mm256_set1_epi8((uint8_t)0x80);
    row_store_32xh(&row, 32, dst, stride);
}

void aom_dc_predictor_32x16_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above, const uint8_t *left) {
    const __m128i top_sum = dc_sum_32_sse2(above);
    __m128i left_sum = dc_sum_16_sse2(left);
    left_sum = _mm_add_epi16(top_sum, left_sum);
    uint32_t sum = _mm_cvtsi128_si32(left_sum);
    sum += 24;
    sum /= 48;
    const __m256i row = _mm256_set1_epi8((uint8_t)sum);
    row_store_32xh(&row, 16, dst, stride);
}

void aom_dc_predictor_32x64_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above, const uint8_t *left) {
    const __m256i sum_above = dc_sum_32(above);
    __m256i sum_left = dc_sum_64(left);
    sum_left = _mm256_add_epi16(sum_left, sum_above);
    uint32_t sum = _mm_cvtsi128_si32(_mm256_castsi256_si128(sum_left));
    sum += 48;
    sum /= 96;
    const __m256i row = _mm256_set1_epi8((uint8_t)sum);
    row_store_32xh(&row, 64, dst, stride);
}

void aom_dc_predictor_64x32_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above, const uint8_t *left) {
    const __m256i sum_above = dc_sum_64(above);
    __m256i sum_left = dc_sum_32(left);
    sum_left = _mm256_add_epi16(sum_left, sum_above);
    uint32_t sum = _mm_cvtsi128_si32(_mm256_castsi256_si128(sum_left));
    sum += 48;
    sum /= 96;
    const __m256i row = _mm256_set1_epi8((uint8_t)sum);
    row_store_64xh(&row, 32, dst, stride);
}

void aom_dc_predictor_64x16_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above, const uint8_t *left) {
    const __m256i sum_above = dc_sum_64(above);
    __m256i sum_left = _mm256_castsi128_si256(dc_sum_16_sse2(left));
    sum_left = _mm256_add_epi16(sum_left, sum_above);
    uint32_t sum = _mm_cvtsi128_si32(_mm256_castsi256_si128(sum_left));
    sum += 40;
    sum /= 80;
    const __m256i row = _mm256_set1_epi8((uint8_t)sum);
    row_store_64xh(&row, 16, dst, stride);
}


void aom_dc_left_predictor_32x16_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above,
    const uint8_t *left) {
    __m128i sum = dc_sum_16_sse2(left);
    (void)above;

    const __m128i eight = _mm_set1_epi16(8);
    sum = _mm_add_epi16(sum, eight);
    sum = _mm_srai_epi16(sum, 4);
    const __m128i zero = _mm_setzero_si128();
    const __m128i r = _mm_shuffle_epi8(sum, zero);
    const __m256i row = _mm256_inserti128_si256(_mm256_castsi128_si256(r), r, 1);
    row_store_32xh(&row, 16, dst, stride);
}

void aom_dc_left_predictor_32x64_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above,
    const uint8_t *left) {
    __m256i sum = dc_sum_64(left);
    (void)above;

    const __m256i thirtytwo = _mm256_set1_epi16(32);
    sum = _mm256_add_epi16(sum, thirtytwo);
    sum = _mm256_srai_epi16(sum, 6);
    const __m256i zero = _mm256_setzero_si256();
    __m256i row = _mm256_shuffle_epi8(sum, zero);
    row_store_32xh(&row, 64, dst, stride);
}

void aom_dc_left_predictor_64x32_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above,
    const uint8_t *left) {
    __m256i sum = dc_sum_32(left);
    (void)above;

    const __m256i sixteen = _mm256_set1_epi16(16);
    sum = _mm256_add_epi16(sum, sixteen);
    sum = _mm256_srai_epi16(sum, 5);
    const __m256i zero = _mm256_setzero_si256();
    __m256i row = _mm256_shuffle_epi8(sum, zero);
    row_store_64xh(&row, 32, dst, stride);
}

void aom_dc_left_predictor_64x16_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above,
    const uint8_t *left) {
    __m128i sum = dc_sum_16_sse2(left);
    (void)above;

    const __m128i eight = _mm_set1_epi16(8);
    sum = _mm_add_epi16(sum, eight);
    sum = _mm_srai_epi16(sum, 4);
    const __m128i zero = _mm_setzero_si128();
    const __m128i r = _mm_shuffle_epi8(sum, zero);
    const __m256i row = _mm256_inserti128_si256(_mm256_castsi128_si256(r), r, 1);
    row_store_64xh(&row, 16, dst, stride);
}


void aom_dc_top_predictor_32x16_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above,
    const uint8_t *left) {
    __m256i sum = dc_sum_32(above);
    (void)left;

    const __m256i sixteen = _mm256_set1_epi16(16);
    sum = _mm256_add_epi16(sum, sixteen);
    sum = _mm256_srai_epi16(sum, 5);
    const __m256i zero = _mm256_setzero_si256();
    __m256i row = _mm256_shuffle_epi8(sum, zero);
    row_store_32xh(&row, 16, dst, stride);
}

void aom_dc_top_predictor_32x64_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above,
    const uint8_t *left) {
    __m256i sum = dc_sum_32(above);
    (void)left;

    const __m256i sixteen = _mm256_set1_epi16(16);
    sum = _mm256_add_epi16(sum, sixteen);
    sum = _mm256_srai_epi16(sum, 5);
    const __m256i zero = _mm256_setzero_si256();
    __m256i row = _mm256_shuffle_epi8(sum, zero);
    row_store_32xh(&row, 64, dst, stride);
}

void aom_dc_top_predictor_64x32_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above,
    const uint8_t *left) {
    __m256i sum = dc_sum_64(above);
    (void)left;

    const __m256i thirtytwo = _mm256_set1_epi16(32);
    sum = _mm256_add_epi16(sum, thirtytwo);
    sum = _mm256_srai_epi16(sum, 6);
    const __m256i zero = _mm256_setzero_si256();
    __m256i row = _mm256_shuffle_epi8(sum, zero);
    row_store_64xh(&row, 32, dst, stride);
}

void aom_dc_top_predictor_64x16_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above,
    const uint8_t *left) {
    __m256i sum = dc_sum_64(above);
    (void)left;

    const __m256i thirtytwo = _mm256_set1_epi16(32);
    sum = _mm256_add_epi16(sum, thirtytwo);
    sum = _mm256_srai_epi16(sum, 6);
    const __m256i zero = _mm256_setzero_si256();
    __m256i row = _mm256_shuffle_epi8(sum, zero);
    row_store_64xh(&row, 16, dst, stride);
}



void aom_dc_128_predictor_32x16_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above,
    const uint8_t *left) {
    (void)above;
    (void)left;
    const __m256i row = _mm256_set1_epi8((uint8_t)0x80);
    row_store_32xh(&row, 16, dst, stride);
}
void aom_dc_128_predictor_32x64_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above,
    const uint8_t *left) {
    (void)above;
    (void)left;
    const __m256i row = _mm256_set1_epi8((uint8_t)0x80);
    row_store_32xh(&row, 64, dst, stride);
}
void aom_dc_128_predictor_64x16_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above,
    const uint8_t *left) {
    (void)above;
    (void)left;
    const __m256i row = _mm256_set1_epi8((uint8_t)0x80);
    row_store_64xh(&row, 16, dst, stride);
}
void aom_dc_128_predictor_64x32_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above,
    const uint8_t *left) {
    (void)above;
    (void)left;
    const __m256i row = _mm256_set1_epi8((uint8_t)0x80);
    row_store_64xh(&row, 32, dst, stride);
}


// There are 32 rows togeter. This function does line:
// 0,1,2,3, and 16,17,18,19. The next call would do
// 4,5,6,7, and 20,21,22,23. So 4 times of calling
// would finish 32 rows.
static INLINE void h_predictor_32x8line(const __m256i *row, uint8_t *dst,
    ptrdiff_t stride) {
    __m256i t[4];
    __m256i m = _mm256_setzero_si256();
    const __m256i inc = _mm256_set1_epi8(4);
    int32_t i;

    for (i = 0; i < 4; i++) {
        t[i] = _mm256_shuffle_epi8(*row, m);
        __m256i r0 = _mm256_permute2x128_si256(t[i], t[i], 0);
        __m256i r1 = _mm256_permute2x128_si256(t[i], t[i], 0x11);
        _mm256_storeu_si256((__m256i *)dst, r0);
        _mm256_storeu_si256((__m256i *)(dst + (stride << 4)), r1);
        dst += stride;
        m = _mm256_add_epi8(m, inc);
    }
}

void aom_h_predictor_32x32_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above, const uint8_t *left) {
    (void)above;
    const __m256i left_col = _mm256_loadu_si256((__m256i const *)left);

    __m256i u = _mm256_unpacklo_epi8(left_col, left_col);

    __m256i v = _mm256_unpacklo_epi8(u, u);
    h_predictor_32x8line(&v, dst, stride);
    dst += stride << 2;

    v = _mm256_unpackhi_epi8(u, u);
    h_predictor_32x8line(&v, dst, stride);
    dst += stride << 2;

    u = _mm256_unpackhi_epi8(left_col, left_col);

    v = _mm256_unpacklo_epi8(u, u);
    h_predictor_32x8line(&v, dst, stride);
    dst += stride << 2;

    v = _mm256_unpackhi_epi8(u, u);
    h_predictor_32x8line(&v, dst, stride);
}
static INLINE void row_store_32x2xh(const __m256i *r0, const __m256i *r1,
    int32_t height, uint8_t *dst,
    ptrdiff_t stride) {
    for (int32_t i = 0; i < height; ++i) {
        _mm256_storeu_si256((__m256i *)dst, *r0);
        _mm256_storeu_si256((__m256i *)(dst + 32), *r1);
        dst += stride;
    }
}
void aom_v_predictor_64x64_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above, const uint8_t *left) {
    const __m256i row0 = _mm256_loadu_si256((const __m256i *)above);
    const __m256i row1 = _mm256_loadu_si256((const __m256i *)(above + 32));
    (void)left;
    row_store_32x2xh(&row0, &row1, 64, dst, stride);
}
void aom_v_predictor_32x32_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above, const uint8_t *left) {
    const __m256i row = _mm256_loadu_si256((const __m256i *)above);
    (void)left;
    row_store_32xh(&row, 32, dst, stride);
}

void aom_v_predictor_32x16_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above, const uint8_t *left) {
    const __m256i row = _mm256_loadu_si256((const __m256i *)above);
    (void)left;
    row_store_32xh(&row, 16, dst, stride);
}
void aom_v_predictor_32x64_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above, const uint8_t *left) {
    const __m256i row = _mm256_loadu_si256((const __m256i *)above);
    (void)left;
    row_store_32xh(&row, 64, dst, stride);
}
void aom_v_predictor_64x16_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above, const uint8_t *left) {
    const __m256i row0 = _mm256_loadu_si256((const __m256i *)above);
    const __m256i row1 = _mm256_loadu_si256((const __m256i *)(above + 32));
    (void)left;
    row_store_32x2xh(&row0, &row1, 16, dst, stride);
}
void aom_v_predictor_64x32_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above, const uint8_t *left) {
    const __m256i row0 = _mm256_loadu_si256((const __m256i *)above);
    const __m256i row1 = _mm256_loadu_si256((const __m256i *)(above + 32));
    (void)left;
    row_store_32x2xh(&row0, &row1, 32, dst, stride);
}

#endif
void intra_mode_dc_64x64_av1_avx2_intrin(
    EbBool        isLeftAvailble,
    EbBool        isAboveAvailble,
    const uint32_t   size,                       //input parameter, denotes the size of the current PU
    uint8_t         *ref_samples,                 //input parameter, pointer to the reference samples
    uint8_t         *dst,              //output parameter, pointer to the prediction
    const uint32_t   prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
    const EbBool  skip)                       //skip half rows
{
    //printf("HEREEE");
    (void)skip;

    uint32_t leftOffset = 0;
    uint32_t topOffset = (size << 1) + 1;
    //uint32_t rowStride = skip ? 2 : 1;

    if (isLeftAvailble && !isAboveAvailble) {

        aom_dc_left_predictor_64x64_avx2(
            dst, prediction_buffer_stride,
            ref_samples + topOffset, ref_samples + leftOffset);

    }
    else if (isAboveAvailble && !isLeftAvailble) {

        aom_dc_top_predictor_64x64_avx2(
            dst, prediction_buffer_stride,
            ref_samples + topOffset, ref_samples + leftOffset);
    }
    else {
        aom_dc_predictor_64x64_avx2(
            dst, prediction_buffer_stride,
            ref_samples + topOffset, ref_samples + leftOffset);
    }
}
#if INTRA_ASM
void aom_dc_predictor_32x32_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above, const uint8_t *left) {
    const __m256i sum_above = dc_sum_32(above);
    __m256i sum_left = dc_sum_32(left);
    sum_left = _mm256_add_epi16(sum_left, sum_above);
    const __m256i thirtytwo = _mm256_set1_epi16(32);
    sum_left = _mm256_add_epi16(sum_left, thirtytwo);
    sum_left = _mm256_srai_epi16(sum_left, 6);
    const __m256i zero = _mm256_setzero_si256();
    __m256i row = _mm256_shuffle_epi8(sum_left, zero);
    row_store_32xh(&row, 32, dst, stride);
}
#endif
/***************************************************************************************************************************************************************************/
/***************************************************************************************intra_mode_angular_2_avx2_intrin***************************************************************************************/
void intra_mode_angular_2_avx2_intrin(
    const uint32_t      size,                       //input parameter, denotes the size of the current PU
    uint8_t            *ref_samples,                 //input parameter, pointer to the reference samples
    uint8_t            *prediction_ptr,              //output parameter, pointer to the prediction
    const uint32_t      prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
    const EbBool     skip)                       //skip one row 
{
    uint32_t pStride = prediction_buffer_stride;
    uint32_t leftOffset = 0;

    if (!skip) {

        if (size == 32) {
            uint32_t count;
            for (count = 0; count < 8; ++count) {
                _mm256_storeu_si256((__m256i *)prediction_ptr, _mm256_loadu_si256((__m256i *)(ref_samples + leftOffset + 1)));
                _mm256_storeu_si256((__m256i *)(prediction_ptr + pStride), _mm256_loadu_si256((__m256i *)(ref_samples + leftOffset + 2)));
                _mm256_storeu_si256((__m256i *)(prediction_ptr + 2 * pStride), _mm256_loadu_si256((__m256i *)(ref_samples + leftOffset + 3)));
                _mm256_storeu_si256((__m256i *)(prediction_ptr + 3 * pStride), _mm256_loadu_si256((__m256i *)(ref_samples + leftOffset + 4)));

                ref_samples += 4;
                prediction_ptr += (pStride << 2);
            }
        }
        else if (size == 16) {
            _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 1)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 2)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 3)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 4)));
            prediction_ptr += (pStride << 2);
            _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 5)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 6)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 7)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 8)));
            prediction_ptr += (pStride << 2);
            _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 9)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 10)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 11)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 12)));
            prediction_ptr += (pStride << 2);
            _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 13)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 14)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 15)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 16)));
        }
        else if (size == 8) {
            _mm_storel_epi64((__m128i *)prediction_ptr, _mm_loadl_epi64((__m128i *)(ref_samples + leftOffset + 1)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), _mm_loadl_epi64((__m128i *)(ref_samples + leftOffset + 2)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + leftOffset + 3)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + leftOffset + 4)));
            prediction_ptr += (pStride << 2);
            _mm_storel_epi64((__m128i *)prediction_ptr, _mm_loadl_epi64((__m128i *)(ref_samples + leftOffset + 5)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), _mm_loadl_epi64((__m128i *)(ref_samples + leftOffset + 6)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + leftOffset + 7)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + leftOffset + 8)));
        }
        else {
            *(uint32_t *)prediction_ptr = *(uint32_t *)(ref_samples + leftOffset + 1);
            *(uint32_t *)(prediction_ptr + pStride) = *(uint32_t *)(ref_samples + leftOffset + 2);
            *(uint32_t *)(prediction_ptr + 2 * pStride) = *(uint32_t *)(ref_samples + leftOffset + 3);
            *(uint32_t *)(prediction_ptr + 3 * pStride) = *(uint32_t *)(ref_samples + leftOffset + 4);
        }
    }
    else {
        if (size != 4) {
            pStride <<= 1;

            if (size == 32) {
                uint32_t count;

                for (count = 0; count < 4; ++count) {

                    _mm256_storeu_si256((__m256i *)prediction_ptr, _mm256_loadu_si256((__m256i *)(ref_samples + leftOffset + 1)));
                    _mm256_storeu_si256((__m256i *)(prediction_ptr + pStride), _mm256_loadu_si256((__m256i *)(ref_samples + leftOffset + 3)));
                    _mm256_storeu_si256((__m256i *)(prediction_ptr + 2 * pStride), _mm256_loadu_si256((__m256i *)(ref_samples + leftOffset + 5)));
                    _mm256_storeu_si256((__m256i *)(prediction_ptr + 3 * pStride), _mm256_loadu_si256((__m256i *)(ref_samples + leftOffset + 7)));
                    ref_samples += 8;
                    prediction_ptr += (pStride << 2);
                }
            }
            else if (size == 16) {

                _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 1)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 3)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 5)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 7)));
                prediction_ptr += (pStride << 2);
                _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 9)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 11)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 13)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + leftOffset + 15)));
            }
            else {
                _mm_storel_epi64((__m128i *)prediction_ptr, _mm_loadl_epi64((__m128i *)(ref_samples + leftOffset + 1)));
                _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), _mm_loadl_epi64((__m128i *)(ref_samples + leftOffset + 3)));
                _mm_storel_epi64((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + leftOffset + 5)));
                _mm_storel_epi64((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + leftOffset + 7)));
            }
        }
        else {
            *(uint32_t*)prediction_ptr = *(uint32_t*)(ref_samples + leftOffset + 1);
            *(uint32_t*)(prediction_ptr + 2 * pStride) = *(uint32_t*)(ref_samples + leftOffset + 3);
        }
    }
}

#define MIDRANGE_VALUE_8BIT    128

uint32_t update_neighbor_dc_intra_pred_avx2_intrin(
    uint8_t                           *yIntraReferenceArrayReverse,
    uint16_t                           input_height,
    uint16_t                           stride_y,
    EbByte                          buffer_y,
    uint16_t                           originY,
    uint16_t                           originX,
    uint32_t                           src_origin_x,
    uint32_t                           src_origin_y,
    uint32_t                           block_size,
    EbAsm                             asmType)
{

    uint32_t idx;
    uint8_t  *src_ptr;
    uint8_t  *dst_ptr;
    uint8_t  *readPtr;

    uint32_t count;

    uint8_t *yBorderReverse = yIntraReferenceArrayReverse;
    uint32_t height = input_height;
    uint32_t blockSizeHalf = block_size << 1;
    uint32_t topOffset = (block_size << 1) + 1;
    uint32_t leftOffset = 0;
    uint32_t    stride = stride_y;
    __m128i xmm0 = _mm_setzero_si128();
    __m256i xmm1 = _mm256_setzero_si256();
    __m256i ymm0;

    __m128i xmm_sad = _mm_setzero_si128();

    // Adjust the Source ptr to start at the origin of the block being updated
    src_ptr = buffer_y + (((src_origin_y + originY) * stride) + (src_origin_x + originX));

    // Adjust the Destination ptr to start at the origin of the Intra reference array
    dst_ptr = yBorderReverse;

    //CHKn here we need ref on Top+Left only. and memset is done only for border CUs

    //Initialise the Luma Intra Reference Array to the mid range value 128 (for CUs at the picture boundaries)
    memset(dst_ptr, MIDRANGE_VALUE_8BIT, (block_size << 2) + 1);

    // Get the left-column
    count = blockSizeHalf;

    if (src_origin_x != 0) {

        readPtr = src_ptr - 1;
        count = ((src_origin_y + count) > height) ? count - ((src_origin_y + count) - height) : count;

        for (idx = 0; idx < count; ++idx) {

            *dst_ptr = *readPtr;
            readPtr += stride;
            dst_ptr++;
        }

        dst_ptr += (blockSizeHalf - count);

    }
    else {

        dst_ptr += count;
    }



    if (block_size != 32) {

        __m128i xmm_mask1 = _mm_slli_si128(_mm_set1_epi8((signed char)0xFF), 1);
        __m128i xmm_mask2 = _mm_srli_si128(xmm_mask1, 15);
        __m128i xmm_C2 = _mm_set1_epi16(0x0002);


        if (block_size == 16) {
            __m128i xmm_predictionDcValue, xmm_top, xmm_left, xmm_sum, xmm_prediction_ptr_0;
            __m128i xmm_top_lo, xmm_top_hi, xmm_left_lo, xmm_left_hi, xmm_predictionDcValue_16, xmm_predictionDcValue_16_x2, xmm_predictionDcValue_16_x3;
            if (src_origin_y != 0)
            {
                xmm_top = _mm_loadu_si128((__m128i *)(src_ptr - stride));
            }
            else
            {
                xmm_top = _mm_loadu_si128((__m128i *)(yBorderReverse + topOffset));//_mm_set1_epi8(128);
            }
            xmm_left = _mm_loadu_si128((__m128i *)(yBorderReverse + leftOffset));
            xmm_top_lo = _mm_unpacklo_epi8(xmm_top, xmm0);
            xmm_top_hi = _mm_unpackhi_epi8(xmm_top, xmm0);
            xmm_left_lo = _mm_unpacklo_epi8(xmm_left, xmm0);
            xmm_left_hi = _mm_unpackhi_epi8(xmm_left, xmm0);

            xmm_sum = _mm_add_epi32(_mm_sad_epu8(xmm_top, xmm0), _mm_sad_epu8(xmm_left, xmm0));

            xmm_predictionDcValue = _mm_srli_epi32(_mm_add_epi32(_mm_add_epi32(_mm_srli_si128(xmm_sum, 8), xmm_sum), _mm_cvtsi32_si128(16)), 5);
            xmm_predictionDcValue = _mm_unpacklo_epi8(xmm_predictionDcValue, xmm_predictionDcValue);
            xmm_predictionDcValue = _mm_unpacklo_epi16(xmm_predictionDcValue, xmm_predictionDcValue);
            xmm_predictionDcValue = _mm_unpacklo_epi32(xmm_predictionDcValue, xmm_predictionDcValue);
            xmm_predictionDcValue = _mm_unpacklo_epi64(xmm_predictionDcValue, xmm_predictionDcValue);

            xmm_predictionDcValue_16 = _mm_srli_epi16(xmm_predictionDcValue, 8);
            xmm_predictionDcValue = _mm_and_si128(xmm_predictionDcValue, xmm_mask1);
            xmm_predictionDcValue_16_x2 = _mm_add_epi16(xmm_predictionDcValue_16, xmm_predictionDcValue_16);
            xmm_predictionDcValue_16_x3 = _mm_add_epi16(xmm_predictionDcValue_16_x2, xmm_predictionDcValue_16);

            xmm_top = _mm_packus_epi16(_mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(xmm_top_lo, xmm_predictionDcValue_16_x3), xmm_C2), 2),
                _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(xmm_top_hi, xmm_predictionDcValue_16_x3), xmm_C2), 2));

            xmm_left = _mm_srli_si128(_mm_packus_epi16(_mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(xmm_left_lo, xmm_predictionDcValue_16_x3), xmm_C2), 2),
                _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(xmm_left_hi, xmm_predictionDcValue_16_x3), xmm_C2), 2)), 1);

            xmm_prediction_ptr_0 = _mm_or_si128(_mm_and_si128(_mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(_mm_add_epi16(xmm_top_lo, xmm_left_lo), xmm_predictionDcValue_16_x2), xmm_C2), 2), xmm_mask2), _mm_and_si128(xmm_top, xmm_mask1));


            xmm_sad = _mm_add_epi32(xmm_sad, _mm_sad_epu8(_mm_loadu_si128((__m128i*)(src_ptr)), xmm_prediction_ptr_0));
            xmm_sad = _mm_add_epi32(xmm_sad, _mm_sad_epu8(_mm_loadu_si128((__m128i*)(src_ptr + stride)), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue)));
            xmm_left = _mm_srli_si128(xmm_left, 1);
            xmm_sad = _mm_add_epi32(xmm_sad, _mm_sad_epu8(_mm_loadu_si128((__m128i*)(src_ptr + (stride << 1))), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue)));
            xmm_left = _mm_srli_si128(xmm_left, 1);
            xmm_sad = _mm_add_epi32(xmm_sad, _mm_sad_epu8(_mm_loadu_si128((__m128i*)(src_ptr + 3 * stride)), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue)));
            xmm_left = _mm_srli_si128(xmm_left, 1);

            src_ptr += (stride << 2);

            for (idx = 4; idx < block_size; idx += 4) {
                xmm_sad = _mm_add_epi32(xmm_sad, _mm_sad_epu8(_mm_loadu_si128((__m128i*)(src_ptr)), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue)));
                xmm_left = _mm_srli_si128(xmm_left, 1);
                xmm_sad = _mm_add_epi32(xmm_sad, _mm_sad_epu8(_mm_loadu_si128((__m128i*)(src_ptr + stride)), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue)));
                xmm_left = _mm_srli_si128(xmm_left, 1);
                xmm_sad = _mm_add_epi32(xmm_sad, _mm_sad_epu8(_mm_loadu_si128((__m128i*)(src_ptr + (stride << 1))), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue)));
                xmm_left = _mm_srli_si128(xmm_left, 1);
                xmm_sad = _mm_add_epi32(xmm_sad, _mm_sad_epu8(_mm_loadu_si128((__m128i*)(src_ptr + 3 * stride)), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue)));
                xmm_left = _mm_srli_si128(xmm_left, 1);

                src_ptr += (stride << 2);

            }

            xmm_sad = _mm_add_epi32(xmm_sad, _mm_srli_si128(xmm_sad, 8));
            return _mm_cvtsi128_si32(xmm_sad);
        }

        else {


            __m128i xmm_left, xmm_top, xmm_top_lo, xmm_left_lo, xmm_predictionDcValue, xmm_predictionDcValue_16;
            __m128i xmm_predictionDcValue_16_x2, xmm_predictionDcValue_16_x3, xmm_prediction_ptr_0;
            if (src_origin_y != 0)
            {
                xmm_top = _mm_loadl_epi64((__m128i *)(src_ptr - stride));
            }
            else
            {
                xmm_top = _mm_loadl_epi64((__m128i *)(yBorderReverse + topOffset));//_mm_set1_epi8(128);//
            }
            xmm_left = _mm_loadl_epi64((__m128i *)(yBorderReverse + leftOffset));

            xmm_top_lo = _mm_unpacklo_epi8(xmm_top, xmm0);
            xmm_left_lo = _mm_unpacklo_epi8(xmm_left, xmm0);

            xmm_predictionDcValue = _mm_srli_epi32(_mm_add_epi32(_mm_add_epi32(_mm_sad_epu8(xmm_top, xmm0), _mm_sad_epu8(xmm_left, xmm0)), _mm_cvtsi32_si128(8)), 4);
            xmm_predictionDcValue = _mm_unpacklo_epi8(xmm_predictionDcValue, xmm_predictionDcValue);
            xmm_predictionDcValue = _mm_unpacklo_epi16(xmm_predictionDcValue, xmm_predictionDcValue);
            xmm_predictionDcValue = _mm_unpacklo_epi32(xmm_predictionDcValue, xmm_predictionDcValue);
            xmm_predictionDcValue = _mm_unpacklo_epi64(xmm_predictionDcValue, xmm_predictionDcValue);

            xmm_predictionDcValue_16 = _mm_srli_epi16(xmm_predictionDcValue, 8);
            xmm_predictionDcValue = _mm_and_si128(xmm_predictionDcValue, xmm_mask1);

            xmm_predictionDcValue_16_x2 = _mm_add_epi16(xmm_predictionDcValue_16, xmm_predictionDcValue_16);
            xmm_predictionDcValue_16_x3 = _mm_add_epi16(xmm_predictionDcValue_16_x2, xmm_predictionDcValue_16);

            xmm_top = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(xmm_top_lo, xmm_predictionDcValue_16_x3), xmm_C2), 2);
            xmm_left = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(xmm_left_lo, xmm_predictionDcValue_16_x3), xmm_C2), 2);

            xmm_left = _mm_srli_si128(_mm_packus_epi16(xmm_left, xmm_left), 1);
            xmm_top = _mm_and_si128(_mm_packus_epi16(xmm_top, xmm_top), xmm_mask1);

            xmm_prediction_ptr_0 = _mm_or_si128(_mm_and_si128(_mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(_mm_add_epi16(xmm_top_lo, xmm_left_lo), xmm_predictionDcValue_16_x2), xmm_C2), 2), xmm_mask2), xmm_top);



            xmm_sad = _mm_setzero_si128();

            xmm_sad = _mm_add_epi32(xmm_sad, _mm_sad_epu8(_mm_loadu_si128((__m128i*)(src_ptr)), xmm_prediction_ptr_0));
            xmm_sad = _mm_add_epi32(xmm_sad, _mm_sad_epu8(_mm_loadu_si128((__m128i*)(src_ptr + stride)), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue)));
            xmm_left = _mm_srli_si128(xmm_left, 1);
            xmm_sad = _mm_add_epi32(xmm_sad, _mm_sad_epu8(_mm_loadu_si128((__m128i*)(src_ptr + (stride << 1))), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue)));
            xmm_left = _mm_srli_si128(xmm_left, 1);
            xmm_sad = _mm_add_epi32(xmm_sad, _mm_sad_epu8(_mm_loadu_si128((__m128i*)(src_ptr + 3 * stride)), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue)));
            xmm_left = _mm_srli_si128(xmm_left, 1);

            src_ptr += (stride << 2);

            for (idx = 4; idx < block_size; idx += 4) {
                xmm_sad = _mm_add_epi32(xmm_sad, _mm_sad_epu8(_mm_loadu_si128((__m128i*)(src_ptr)), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue)));
                xmm_left = _mm_srli_si128(xmm_left, 1);
                xmm_sad = _mm_add_epi32(xmm_sad, _mm_sad_epu8(_mm_loadu_si128((__m128i*)(src_ptr + stride)), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue)));
                xmm_left = _mm_srli_si128(xmm_left, 1);
                xmm_sad = _mm_add_epi32(xmm_sad, _mm_sad_epu8(_mm_loadu_si128((__m128i*)(src_ptr + (stride << 1))), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue)));
                xmm_left = _mm_srli_si128(xmm_left, 1);
                xmm_sad = _mm_add_epi32(xmm_sad, _mm_sad_epu8(_mm_loadu_si128((__m128i*)(src_ptr + 3 * stride)), _mm_or_si128(_mm_and_si128(xmm_left, xmm_mask2), xmm_predictionDcValue)));
                xmm_left = _mm_srli_si128(xmm_left, 1);

                src_ptr += (stride << 2);

            }

            return _mm_cvtsi128_si32(xmm_sad);

        }
    }



    /*************************************************************************************************************************************************************************************************************/
    else {
        __m256i xmm_sum, xmm_sadleft, xmm_sadtop, xmm_toptmp, xmm_lefttmp, xmm_set, xmm_sum128_2, xmm_sum256, xmm_predictionDcValue;
        __m128i xmm_sumhi, xmm_sumlo, xmm_sum1, xmm_sum128, xmm_sumhitmp, xmm_sumlotmp, xmm_movelotmp, xmm_movehitmp;

        xmm_sumhi = xmm_sumlo = xmm_sum128 = xmm_sumhitmp = xmm_sumlotmp = _mm_setzero_si128();
        xmm_sum = xmm_sadleft = xmm_sadtop = xmm_toptmp = xmm_lefttmp = _mm256_setzero_si256();

        if (src_origin_y != 0)
        {
            xmm_toptmp = _mm256_sad_epu8(_mm256_set_m128i(_mm_loadu_si128((__m128i *)(src_ptr - stride + 16)), _mm_loadu_si128((__m128i *)(src_ptr - stride))), xmm1);

        }
        else
        {
            xmm_toptmp = _mm256_sad_epu8(_mm256_set_m128i(_mm_loadu_si128((__m128i *)(yBorderReverse + topOffset + 16)), _mm_loadu_si128((__m128i *)(yBorderReverse + topOffset))), xmm1);

        }
        xmm_lefttmp = _mm256_sad_epu8(_mm256_set_m128i(_mm_loadu_si128((__m128i *)(yBorderReverse + leftOffset + 16)), _mm_loadu_si128((__m128i *)(yBorderReverse + leftOffset))), xmm1);

        xmm_sum = _mm256_add_epi32(xmm_toptmp, xmm_lefttmp);
        xmm_sum = _mm256_hadd_epi32(xmm_sum, xmm_sum);

        xmm_sumlo = _mm256_extracti128_si256(xmm_sum, 0);
        xmm_sumhi = _mm256_extracti128_si256(xmm_sum, 1);

        xmm_movelotmp = _mm_move_epi64(xmm_sumlo);
        xmm_movehitmp = _mm_move_epi64(xmm_sumhi);

        xmm_sum1 = _mm_add_epi32(xmm_movelotmp, xmm_movehitmp);

        xmm_sum1 = _mm_hadd_epi32(xmm_sum1, xmm_sum1);

        xmm_sum256 = _mm256_castsi128_si256(xmm_sum1);


        xmm_set = _mm256_castsi128_si256(_mm_set1_epi32(32));

        xmm_sum128_2 = _mm256_add_epi32(xmm_sum256, xmm_set); // add offset
        xmm_predictionDcValue = _mm256_srli_epi32(xmm_sum128_2, 6); //_mm256_srli_epi32


        __m128i dc128 = _mm256_castsi256_si128(xmm_predictionDcValue);

        uint8_t dc = _mm_cvtsi128_si32(dc128);
        xmm_predictionDcValue = _mm256_set1_epi8(dc);//_mm_broadcastb_epi8


        // SAD
        ymm0 = _mm256_setzero_si256();
        for (idx = 0; idx < block_size; idx += 2) {
            ymm0 = _mm256_add_epi32(ymm0, _mm256_sad_epu8(_mm256_loadu_si256((__m256i*)src_ptr), xmm_predictionDcValue));
            xmm1 = _mm256_add_epi32(xmm1, _mm256_sad_epu8(_mm256_loadu_si256((__m256i*)(src_ptr + stride)), xmm_predictionDcValue));
            src_ptr += stride << 1;
        }
        ymm0 = _mm256_add_epi32(ymm0, xmm1);
        xmm0 = _mm_add_epi32(_mm256_extracti128_si256(ymm0, 0), _mm256_extracti128_si256(ymm0, 1));
        xmm0 = _mm_add_epi32(xmm0, _mm_srli_si128(xmm0, 8));
        return (uint32_t)_mm_cvtsi128_si32(xmm0);


    }
    (void)asmType;
}

/***********************************************************************************************************************************************************************************************
                                                                        intra_mode_angular_18_avx2_intrin
                                                                        ***********************************************************************************************************************************************************************************************/
void intra_mode_angular_18_avx2_intrin(
    const uint32_t      size,                       //input parameter, denotes the size of the current PU
    uint8_t            *ref_samples,                 //input parameter, pointer to the reference samples
    uint8_t            *prediction_ptr,              //output parameter, pointer to the prediction
    const uint32_t      prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
    const EbBool     skip)                    //skip one row 
{
    uint32_t pStride = prediction_buffer_stride;
    uint32_t topLeftOffset = (size << 1);

    if (!skip) {

        if (size == 32) {

            uint32_t count;
            for (count = 0; count < 8; ++count) {

                _mm256_storeu_si256((__m256i *)prediction_ptr, _mm256_loadu_si256((__m256i *)(ref_samples + topLeftOffset)));
                _mm256_storeu_si256((__m256i *)(prediction_ptr + pStride), _mm256_loadu_si256((__m256i *)(ref_samples + topLeftOffset - 1)));
                _mm256_storeu_si256((__m256i *)(prediction_ptr + 2 * pStride), _mm256_loadu_si256((__m256i *)(ref_samples + topLeftOffset - 2)));
                _mm256_storeu_si256((__m256i *)(prediction_ptr + 3 * pStride), _mm256_loadu_si256((__m256i *)(ref_samples + topLeftOffset - 3)));

                ref_samples -= 4;
                prediction_ptr += (pStride << 2);
            }
        }
        else if (size == 16) {
            _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 1)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 2)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 3)));
            prediction_ptr += (pStride << 2);
            _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 4)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 5)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 6)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 7)));
            prediction_ptr += (pStride << 2);
            _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 8)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 9)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 10)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 11)));
            prediction_ptr += (pStride << 2);
            _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 12)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 13)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 14)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 15)));
        }
        else if (size == 8) {
            _mm_storel_epi64((__m128i *)prediction_ptr, _mm_loadl_epi64((__m128i *)(ref_samples + topLeftOffset)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topLeftOffset - 1)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topLeftOffset - 2)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topLeftOffset - 3)));
            prediction_ptr += (pStride << 2);
            _mm_storel_epi64((__m128i *)prediction_ptr, _mm_loadl_epi64((__m128i *)(ref_samples + topLeftOffset - 4)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topLeftOffset - 5)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topLeftOffset - 6)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topLeftOffset - 7)));
        }
        else {
            *(uint32_t *)prediction_ptr = *(uint32_t *)(ref_samples + topLeftOffset);
            *(uint32_t *)(prediction_ptr + pStride) = *(uint32_t *)(ref_samples + topLeftOffset - 1);
            *(uint32_t *)(prediction_ptr + 2 * pStride) = *(uint32_t *)(ref_samples + topLeftOffset - 2);
            *(uint32_t *)(prediction_ptr + 3 * pStride) = *(uint32_t *)(ref_samples + topLeftOffset - 3);
        }
    }
    else {
        if (size != 4) {
            pStride <<= 1;

            if (size == 32) {
                uint32_t count;

                for (count = 0; count < 4; ++count) {

                    _mm256_storeu_si256((__m256i *)prediction_ptr, _mm256_loadu_si256((__m256i *)(ref_samples + topLeftOffset)));
                    _mm256_storeu_si256((__m256i *)(prediction_ptr + pStride), _mm256_loadu_si256((__m256i *)(ref_samples + topLeftOffset - 2)));
                    _mm256_storeu_si256((__m256i *)(prediction_ptr + 2 * pStride), _mm256_loadu_si256((__m256i *)(ref_samples + topLeftOffset - 4)));
                    _mm256_storeu_si256((__m256i *)(prediction_ptr + 3 * pStride), _mm256_loadu_si256((__m256i *)(ref_samples + topLeftOffset - 6)));

                    ref_samples -= 8;
                    prediction_ptr += (pStride << 2);
                }
            }
            else if (size == 16) {
                _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 2)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 4)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 6)));
                prediction_ptr += (pStride << 2);
                _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 8)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 10)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 12)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topLeftOffset - 14)));
            }
            else {
                _mm_storel_epi64((__m128i *)prediction_ptr, _mm_loadl_epi64((__m128i *)(ref_samples + topLeftOffset)));
                _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topLeftOffset - 2)));
                _mm_storel_epi64((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topLeftOffset - 4)));
                _mm_storel_epi64((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topLeftOffset - 6)));
            }
        }
        else {
            *(uint32_t*)prediction_ptr = *(uint32_t*)(ref_samples + topLeftOffset);
            *(uint32_t*)(prediction_ptr + 2 * pStride) = *(uint32_t*)(ref_samples + topLeftOffset - 2);
        }
    }
}

/*********************************************************************************************************************************************************************************************
                                                                intra_mode_angular_34_avx2_intrin
                                                                ***********************************************************************************************************************************************************************************************/

void intra_mode_angular_34_avx2_intrin(
    const uint32_t      size,                       //input parameter, denotes the size of the current PU
    uint8_t            *ref_samples,                 //input parameter, pointer to the reference samples
    uint8_t            *prediction_ptr,              //output parameter, pointer to the prediction
    const uint32_t      prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
    const EbBool     skip)                       //skip one row 
{
    uint32_t pStride = prediction_buffer_stride;
    uint32_t topOffset = ((size << 1) + 1);

    if (!skip) {

        if (size == 32) {

            uint32_t count;
            for (count = 0; count < 8; ++count) {

                _mm256_storeu_si256((__m256i *)prediction_ptr, _mm256_loadu_si256((__m256i *)(ref_samples + topOffset + 1)));
                _mm256_storeu_si256((__m256i *)(prediction_ptr + pStride), _mm256_loadu_si256((__m256i *)(ref_samples + topOffset + 2)));
                _mm256_storeu_si256((__m256i *)(prediction_ptr + 2 * pStride), _mm256_loadu_si256((__m256i *)(ref_samples + topOffset + 3)));
                _mm256_storeu_si256((__m256i *)(prediction_ptr + 3 * pStride), _mm256_loadu_si256((__m256i *)(ref_samples + topOffset + 4)));

                ref_samples += 4;
                prediction_ptr += (pStride << 2);
            }
        }
        else if (size == 16) {
            _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 1)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 2)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 3)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 4)));
            prediction_ptr += (pStride << 2);
            _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 5)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 6)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 7)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 8)));
            prediction_ptr += (pStride << 2);
            _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 9)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 10)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 11)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 12)));
            prediction_ptr += (pStride << 2);
            _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 13)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 14)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 15)));
            _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 16)));
        }
        else if (size == 8) {
            _mm_storel_epi64((__m128i *)prediction_ptr, _mm_loadl_epi64((__m128i *)(ref_samples + topOffset + 1)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topOffset + 2)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topOffset + 3)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topOffset + 4)));
            prediction_ptr += (pStride << 2);
            _mm_storel_epi64((__m128i *)prediction_ptr, _mm_loadl_epi64((__m128i *)(ref_samples + topOffset + 5)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topOffset + 6)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topOffset + 7)));
            _mm_storel_epi64((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topOffset + 8)));
        }
        else {
            *(uint32_t *)prediction_ptr = *(uint32_t *)(ref_samples + topOffset + 1);
            *(uint32_t *)(prediction_ptr + pStride) = *(uint32_t *)(ref_samples + topOffset + 2);
            *(uint32_t *)(prediction_ptr + 2 * pStride) = *(uint32_t *)(ref_samples + topOffset + 3);
            *(uint32_t *)(prediction_ptr + 3 * pStride) = *(uint32_t *)(ref_samples + topOffset + 4);
        }
    }
    else {
        if (size != 4) {
            pStride <<= 1;

            if (size == 32) {
                uint32_t count;

                for (count = 0; count < 4; ++count) {
                    _mm256_storeu_si256((__m256i *)prediction_ptr, _mm256_loadu_si256((__m256i *)(ref_samples + topOffset + 1)));
                    _mm256_storeu_si256((__m256i *)(prediction_ptr + pStride), _mm256_loadu_si256((__m256i *)(ref_samples + topOffset + 3)));
                    _mm256_storeu_si256((__m256i *)(prediction_ptr + 2 * pStride), _mm256_loadu_si256((__m256i *)(ref_samples + topOffset + 5)));
                    _mm256_storeu_si256((__m256i *)(prediction_ptr + 3 * pStride), _mm256_loadu_si256((__m256i *)(ref_samples + topOffset + 7)));

                    ref_samples += 8;
                    prediction_ptr += (pStride << 2);
                }
            }
            else if (size == 16) {

                _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 1)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 3)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 5)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 7)));
                prediction_ptr += (pStride << 2);
                _mm_storeu_si128((__m128i *)prediction_ptr, _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 9)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 11)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 13)));
                _mm_storeu_si128((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadu_si128((__m128i *)(ref_samples + topOffset + 15)));
            }
            else {
                _mm_storel_epi64((__m128i *)prediction_ptr, _mm_loadl_epi64((__m128i *)(ref_samples + topOffset + 1)));
                _mm_storel_epi64((__m128i *)(prediction_ptr + pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topOffset + 3)));
                _mm_storel_epi64((__m128i *)(prediction_ptr + 2 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topOffset + 5)));
                _mm_storel_epi64((__m128i *)(prediction_ptr + 3 * pStride), _mm_loadl_epi64((__m128i *)(ref_samples + topOffset + 7)));
            }
        }
        else {
            *(uint32_t*)prediction_ptr = *(uint32_t*)(ref_samples + topOffset + 1);
            *(uint32_t*)(prediction_ptr + 2 * pStride) = *(uint32_t*)(ref_samples + topOffset + 3);
        }
    }
}

// only define these intrinsics if immintrin.h doesn't have them 
#ifdef _WIN32
#if _MSC_VER < 1910
static inline int32_t _mm256_extract_epi32(__m256i a, const int32_t i)
{
    return a.m256i_i32[i & 7];
}

static inline __m256i _mm256_insert_epi32(__m256i a, int32_t b, const int32_t i)
{
    __m256i c = a;
    c.m256i_i32[i & 7] = b;
    return c;
}
#endif
#endif

#define PERM4x64(c0, c1, c2, c3) c0+(c1<<2)+(c2<<4)+(c3<<6)
#define PERM2x128(c0, c1) c0+(c1<<4)
void Predict_Z1_16bit_4x4(uint16_t *ref_samples, uint16_t *dst, const uint32_t prediction_buffer_stride, uint16_t dx, uint16_t bd)
{
    int32_t x;
    assert(dx > 0);

    const int32_t frac_bits = 6;
    const int32_t max_base_x = ((4 + 4) - 1);


    uint32_t aboveBy32[16];
    int_least32_t aboveDiff[16];

    // pre-filter above pixels
    // store in temp buffers:
    //   above[x] * 32 + 16
    //   above[x+1] - above[x]
    // final pixels will be caluculated as:
    //   (above[x] * 32 + 16 + (above[x+1] - above[x]) * shift) >> 5
    __m128i clip_bd = (bd == 8) ? _mm_set1_epi16(255) : (bd == 10) ? _mm_set1_epi16(1023) : _mm_set1_epi16(4095);
    __m256i a0, a1, diff, a32, a16;
    a16 = _mm256_set1_epi32(16);
    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(ref_samples)));
    int32_t a = _mm256_extract_epi32(a0, 4);
    a1 = _mm256_alignr_epi8(_mm256_set1_epi32(0), a0, 4);
    a1 = _mm256_insert_epi32(a1, a, 3);
    diff = _mm256_sub_epi32(a1, a0);
    diff = _mm256_insert_epi32(diff, 0, 7);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(aboveBy32), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff), diff);

    a32 = _mm256_set1_epi32(_mm256_extract_epi32(a0, 7));
    a32 = _mm256_slli_epi32(a32, 5);
    diff = _mm256_setzero_si256();
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 7), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 7), diff);

    x = dx;
    __m128i inc = _mm_set1_epi32(dx);
    __m128i shift = _mm_srli_epi32(_mm_and_si128(inc, _mm_set1_epi32(0x3f)), 1);
    for (int32_t r = 0; r < 4; r++, dst += prediction_buffer_stride) {
        __m128i a, b, res;

        int32_t base = x >> frac_bits;
        base = (base > max_base_x) ? max_base_x : base;

        a = _mm_loadu_si128((__m128i *)(aboveBy32 + base));
        b = _mm_loadu_si128((__m128i *)(aboveDiff + base));
        b = _mm_mullo_epi32(b, shift);
        res = _mm_add_epi32(a, b);
        res = _mm_srli_epi32(res, 5);

        x += dx;
        shift = _mm_srli_epi32(_mm_and_si128(_mm_set1_epi32(x), _mm_set1_epi32(0x3f)), 1);
        res = _mm_packus_epi32(res, res);
        res = _mm_min_epi16(clip_bd, res);
        _mm_storel_epi64((__m128i *)dst, res);
    }
}
void Predict_Z1_16bit_8x8(uint16_t *ref_samples, uint16_t *dst, const uint32_t prediction_buffer_stride, uint16_t dx, uint16_t bd)
{
    int32_t x;
    assert(dx > 0);

    const int32_t frac_bits = 6;
    const int32_t max_base_x = ((8 + 8) - 1);



    uint32_t aboveBy32[24];
    int_least32_t aboveDiff[24];

    // pre-filter above pixels
    // store in temp buffers:
    //   above[x] * 32 + 16
    //   above[x+1] - above[x]
    // final pixels will be caluculated as:
    //   (above[x] * 32 + 16 + (above[x+1] - above[x]) * shift) >> 5
    __m256i clip_bd = (bd == 8) ? _mm256_set1_epi16(255) : (bd == 10) ? _mm256_set1_epi16(1023) : _mm256_set1_epi16(4095);
    __m256i a0, a1, diff, a32, a16;
    a16 = _mm256_set1_epi32(16);
    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(ref_samples)));          // 01234567
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[1])));    // 89abcdef
    diff = _mm256_sub_epi32(a1, a0);                                        // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0, 5);                                            // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);                                        // a[x] * 32 + 16
    _mm256_storeu_si256((__m256i *)aboveBy32, a32);
    _mm256_storeu_si256((__m256i *)aboveDiff, diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[8])));
    int32_t a = _mm256_extract_epi32(a0, 4);
    a1 = _mm256_alignr_epi8(_mm256_set1_epi32(0), a0, 4);
    a1 = _mm256_insert_epi32(a1, a, 3);
    diff = _mm256_sub_epi32(a1, a0);
    diff = _mm256_insert_epi32(diff, 0, 7);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 8), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 8), diff);

    a32 = _mm256_set1_epi32(_mm256_extract_epi32(a0, 7));
    a32 = _mm256_slli_epi32(a32, 5);
    diff = _mm256_setzero_si256();
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 16), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 16), diff);

    x = dx;
    __m256i inc = _mm256_set1_epi32(dx);
    __m256i shift = _mm256_srli_epi32(_mm256_and_si256(inc, _mm256_set1_epi32(0x3f)), 1);
    for (int32_t r = 0; r < 8; r++, dst += prediction_buffer_stride) {
        __m256i a, b, res, res1;

        int32_t base = x >> frac_bits;
        base = (base > max_base_x) ? max_base_x : base;

        a = _mm256_loadu_si256((__m256i *)(aboveBy32 + base));
        b = _mm256_loadu_si256((__m256i *)(aboveDiff + base));
        b = _mm256_mullo_epi32(b, shift);
        res = _mm256_add_epi32(a, b);
        res = _mm256_srli_epi32(res, 5);

        x += dx;
        shift = _mm256_srli_epi32(_mm256_and_si256(_mm256_set1_epi32(x), _mm256_set1_epi32(0x3f)), 1);
        res1 = _mm256_permute2x128_si256(res, res, 0x01);
        res = _mm256_packus_epi32(res, res1);
        res = _mm256_min_epi16(clip_bd, res);
        _mm_storeu_si128((__m128i*)dst, _mm256_castsi256_si128(res));
    }

}
void Predict_Z1_16bit_16x16(uint16_t *ref_samples, uint16_t *dst, const uint32_t prediction_buffer_stride, uint16_t dx, uint16_t bd)
{
    int32_t x;
    assert(dx > 0);

    const int32_t frac_bits = 6;
    const int32_t max_base_x = ((16 + 16) - 1);


    uint32_t aboveBy32[48];
    int_least32_t aboveDiff[48];

    // pre-filter above pixels
    // store in temp buffers:
    //   above[x] * 32 + 16
    //   above[x+1] - above[x]
    // final pixels will be caluculated as:
    //   (above[x] * 32 + 16 + (above[x+1] - above[x]) * shift) >> 5
    __m256i clip_bd = (bd == 8) ? _mm256_set1_epi16(255) : (bd == 10) ? _mm256_set1_epi16(1023) : _mm256_set1_epi16(4095);
    __m256i a0, a1, diff, a32, a16;
    a16 = _mm256_set1_epi32(16);
    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(ref_samples)));          // 01234567
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[1])));    // 89abcdef
    diff = _mm256_sub_epi32(a1, a0);                                        // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0, 5);                                            // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);                                        // a[x] * 32 + 16
    _mm256_storeu_si256((__m256i *)aboveBy32, a32);
    _mm256_storeu_si256((__m256i *)aboveDiff, diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[8])));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[9])));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 8), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 8), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[16])));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[17])));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 16), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 16), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[24])));
    int32_t a = _mm256_extract_epi32(a0, 4);
    a1 = _mm256_alignr_epi8(_mm256_set1_epi32(0), a0, 4);
    a1 = _mm256_insert_epi32(a1, a, 3);
    diff = _mm256_sub_epi32(a1, a0);
    diff = _mm256_insert_epi32(diff, 0, 7);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 24), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 24), diff);

    a32 = _mm256_set1_epi32(_mm256_extract_epi32(a0, 7));
    a32 = _mm256_slli_epi32(a32, 5);
    diff = _mm256_setzero_si256();
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 31), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 31), diff);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 39), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 39), diff);


    x = dx;
    __m256i inc = _mm256_set1_epi32(dx);
    __m256i shift = _mm256_srli_epi32(_mm256_and_si256(inc, _mm256_set1_epi32(0x3f)), 1);
    for (int32_t r = 0; r < 16; r++, dst += prediction_buffer_stride) {
        __m256i a, b, res1, res2, res;

        int32_t base = x >> frac_bits;
        base = (base > max_base_x) ? max_base_x : base;

        a = _mm256_loadu_si256((__m256i *)(aboveBy32 + base));
        b = _mm256_loadu_si256((__m256i *)(aboveDiff + base));
        b = _mm256_mullo_epi32(b, shift);
        res1 = _mm256_add_epi32(a, b);
        res1 = _mm256_srli_epi32(res1, 5);

        a = _mm256_loadu_si256((__m256i *)(aboveBy32 + base + 8));
        b = _mm256_loadu_si256((__m256i *)(aboveDiff + base + 8));
        b = _mm256_mullo_epi32(b, shift);
        res2 = _mm256_add_epi32(a, b);
        res2 = _mm256_srli_epi32(res2, 5);

        x += dx;
        shift = _mm256_srli_epi32(_mm256_and_si256(_mm256_set1_epi32(x), _mm256_set1_epi32(0x3f)), 1);

        res = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2), PERM4x64(0, 2, 1, 3));
        res = _mm256_min_epi16(clip_bd, res);
        _mm256_storeu_si256((__m256i *)dst, res);
    }
}
void Predict_Z1_16bit_32x32(uint16_t *ref_samples, uint16_t *dst, const uint32_t prediction_buffer_stride, uint16_t dx, uint16_t bd)
{
    int32_t x;
    assert(dx > 0);
    const int32_t frac_bits = 6;
    const int32_t max_base_x = ((32 + 32) - 1);


    uint32_t aboveBy32[96];
    int_least32_t aboveDiff[96];

    // pre-filter above pixels
    // store in temp buffers:
    //   above[x] * 32 + 16
    //   above[x+1] - above[x]
    // final pixels will be caluculated as:
    //   (above[x] * 32 + 16 + (above[x+1] - above[x]) * shift) >> 5
    __m256i clip_bd = (bd == 8) ? _mm256_set1_epi16(255) : (bd == 10) ? _mm256_set1_epi16(1023) : _mm256_set1_epi16(4095);
    __m256i a0, a1, diff, a32, a16;
    a16 = _mm256_set1_epi32(16);
    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(ref_samples)));          // 01234567
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[1])));    // 89abcdef
    diff = _mm256_sub_epi32(a1, a0);                                        // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0, 5);                                            // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);                                        // a[x] * 32 + 16
    _mm256_storeu_si256((__m256i *)aboveBy32, a32);
    _mm256_storeu_si256((__m256i *)aboveDiff, diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[8])));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[9])));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 8), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 8), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[16])));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[17])));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 16), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 16), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[24])));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[25])));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 24), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 24), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[32])));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[33])));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 32), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 32), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[40])));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[41])));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 40), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 40), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[48])));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[49])));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 48), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 48), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[56])));
    int32_t a = _mm256_extract_epi32(a0, 4);
    a1 = _mm256_alignr_epi8(_mm256_set1_epi32(0), a0, 4);
    a1 = _mm256_insert_epi32(a1, a, 3);
    diff = _mm256_sub_epi32(a1, a0);
    diff = _mm256_insert_epi32(diff, 0, 7);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 56), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 56), diff);

    a32 = _mm256_set1_epi32(_mm256_extract_epi32(a0, 7));
    a32 = _mm256_slli_epi32(a32, 5);
    diff = _mm256_setzero_si256();
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 63), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 63), diff);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 71), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 71), diff);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 79), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 79), diff);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 87), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 87), diff);

    x = dx;
    __m256i inc = _mm256_set1_epi32(dx);
    __m256i shift = _mm256_srli_epi32(_mm256_and_si256(inc, _mm256_set1_epi32(0x3f)), 1);
    for (int32_t r = 0; r < 32; r++, dst += prediction_buffer_stride) {
        __m256i a, b, res1, res2, res3, res4, resLo, resHi;

        int32_t base = x >> frac_bits;
        base = (base > max_base_x) ? max_base_x : base;

        a = _mm256_loadu_si256((__m256i *)(aboveBy32 + base));
        b = _mm256_loadu_si256((__m256i *)(aboveDiff + base));
        b = _mm256_mullo_epi32(b, shift);
        res1 = _mm256_add_epi32(a, b);
        res1 = _mm256_srli_epi32(res1, 5);

        a = _mm256_loadu_si256((__m256i *)(aboveBy32 + base + 8));
        b = _mm256_loadu_si256((__m256i *)(aboveDiff + base + 8));
        b = _mm256_mullo_epi32(b, shift);
        res2 = _mm256_add_epi32(a, b);
        res2 = _mm256_srli_epi32(res2, 5);

        a = _mm256_loadu_si256((__m256i *)(aboveBy32 + base + 16));
        b = _mm256_loadu_si256((__m256i *)(aboveDiff + base + 16));
        b = _mm256_mullo_epi32(b, shift);
        res3 = _mm256_add_epi32(a, b);
        res3 = _mm256_srli_epi32(res3, 5);

        a = _mm256_loadu_si256((__m256i *)(aboveBy32 + base + 24));
        b = _mm256_loadu_si256((__m256i *)(aboveDiff + base + 24));
        b = _mm256_mullo_epi32(b, shift);
        res4 = _mm256_add_epi32(a, b);
        res4 = _mm256_srli_epi32(res4, 5);

        x += dx;
        shift = _mm256_srli_epi32(_mm256_and_si256(_mm256_set1_epi32(x), _mm256_set1_epi32(0x3f)), 1);

        resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2), PERM4x64(0, 2, 1, 3));
        resHi = _mm256_permute4x64_epi64(_mm256_packus_epi32(res3, res4), PERM4x64(0, 2, 1, 3));
        resLo = _mm256_min_epi16(clip_bd, resLo);
        _mm256_storeu_si256((__m256i *)dst, resLo);
        resHi = _mm256_min_epi16(clip_bd, resHi);
        _mm256_storeu_si256((__m256i *)(dst + 16), resHi);
    }
}
void Predict_Z1_16bit_64x64(uint16_t *ref_samples, uint16_t *dst, const uint32_t prediction_buffer_stride, uint16_t dx, uint16_t bd)
{
    int32_t x;
    assert(dx > 0);

    const int32_t frac_bits = 6;
    const int32_t max_base_x = ((64 + 64) - 1);

    uint32_t aboveBy32[192];
    int_least32_t aboveDiff[192];

    // pre-filter above pixels
    // store in temp buffers:
    //   above[x] * 32 + 16
    //   above[x+1] - above[x]
    // final pixels will be caluculated as:
    //   (above[x] * 32 + 16 + (above[x+1] - above[x]) * shift) >> 5
    __m256i clip_bd = (bd == 8) ? _mm256_set1_epi16(255) : (bd == 10) ? _mm256_set1_epi16(1023) : _mm256_set1_epi16(4095);
    __m256i a0, a1, diff, a32, a16;
    a16 = _mm256_set1_epi32(16);
    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(ref_samples)));          // 01234567
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[1])));    // 89abcdef
    diff = _mm256_sub_epi32(a1, a0);                                        // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0, 5);                                            // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);                                        // a[x] * 32 + 16
    _mm256_storeu_si256((__m256i *)aboveBy32, a32);
    _mm256_storeu_si256((__m256i *)aboveDiff, diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[8])));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[9])));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 8), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 8), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[16])));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[17])));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 16), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 16), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[24])));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[25])));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 24), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 24), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[32])));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[33])));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 32), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 32), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[40])));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[41])));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 40), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 40), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[48])));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[49])));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 48), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 48), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[56])));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[57])));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 56), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 56), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[64])));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[65])));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 64), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 64), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[72])));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[73])));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 72), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 72), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[80])));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[81])));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 80), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 80), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[88])));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[89])));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 88), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 88), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[96])));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[97])));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 96), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 96), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[104])));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[105])));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 104), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 104), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[112])));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[113])));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 112), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 112), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(&ref_samples[120])));
    int32_t a = _mm256_extract_epi32(a0, 4);
    a1 = _mm256_alignr_epi8(_mm256_set1_epi32(0), a0, 4);
    a1 = _mm256_insert_epi32(a1, a, 3);
    diff = _mm256_sub_epi32(a1, a0);
    diff = _mm256_insert_epi32(diff, 0, 7);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 120), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 120), diff);

    a32 = _mm256_set1_epi32(_mm256_extract_epi32(a0, 7));
    a32 = _mm256_slli_epi32(a32, 5);
    diff = _mm256_setzero_si256();
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 127), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 127), diff);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 135), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 135), diff);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 143), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 143), diff);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 151), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 151), diff);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 159), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 159), diff);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 167), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 167), diff);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 175), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 175), diff);
    _mm256_storeu_si256((__m256i *)(aboveBy32 + 183), a32);
    _mm256_storeu_si256((__m256i *)(aboveDiff + 183), diff);

    x = dx;
    __m256i inc = _mm256_set1_epi32(dx);
    __m256i shift = _mm256_srli_epi32(_mm256_and_si256(inc, _mm256_set1_epi32(0x3f)), 1);
    for (int32_t r = 0; r < 64; r++, dst += prediction_buffer_stride) {
        __m256i a, b, res1, res2, res3, res4, resLo, resHi;

        int32_t base = x >> frac_bits;
        base = (base > max_base_x) ? max_base_x : base;

        a = _mm256_loadu_si256((__m256i *)(aboveBy32 + base));
        b = _mm256_loadu_si256((__m256i *)(aboveDiff + base));
        b = _mm256_mullo_epi32(b, shift);
        res1 = _mm256_add_epi32(a, b);
        res1 = _mm256_srli_epi32(res1, 5);

        a = _mm256_loadu_si256((__m256i *)(aboveBy32 + base + 8));
        b = _mm256_loadu_si256((__m256i *)(aboveDiff + base + 8));
        b = _mm256_mullo_epi32(b, shift);
        res2 = _mm256_add_epi32(a, b);
        res2 = _mm256_srli_epi32(res2, 5);

        a = _mm256_loadu_si256((__m256i *)(aboveBy32 + base + 16));
        b = _mm256_loadu_si256((__m256i *)(aboveDiff + base + 16));
        b = _mm256_mullo_epi32(b, shift);
        res3 = _mm256_add_epi32(a, b);
        res3 = _mm256_srli_epi32(res3, 5);

        a = _mm256_loadu_si256((__m256i *)(aboveBy32 + base + 24));
        b = _mm256_loadu_si256((__m256i *)(aboveDiff + base + 24));
        b = _mm256_mullo_epi32(b, shift);
        res4 = _mm256_add_epi32(a, b);
        res4 = _mm256_srli_epi32(res4, 5);

        resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2), PERM4x64(0, 2, 1, 3));
        resHi = _mm256_permute4x64_epi64(_mm256_packus_epi32(res3, res4), PERM4x64(0, 2, 1, 3));


        a = _mm256_loadu_si256((__m256i *)(aboveBy32 + base + 32));
        b = _mm256_loadu_si256((__m256i *)(aboveDiff + base + 32));
        b = _mm256_mullo_epi32(b, shift);
        res1 = _mm256_add_epi32(a, b);
        res1 = _mm256_srli_epi32(res1, 5);

        a = _mm256_loadu_si256((__m256i *)(aboveBy32 + base + 40));
        b = _mm256_loadu_si256((__m256i *)(aboveDiff + base + 40));
        b = _mm256_mullo_epi32(b, shift);
        res2 = _mm256_add_epi32(a, b);
        res2 = _mm256_srli_epi32(res2, 5);

        a = _mm256_loadu_si256((__m256i *)(aboveBy32 + base + 48));
        b = _mm256_loadu_si256((__m256i *)(aboveDiff + base + 48));
        b = _mm256_mullo_epi32(b, shift);
        res3 = _mm256_add_epi32(a, b);
        res3 = _mm256_srli_epi32(res3, 5);

        a = _mm256_loadu_si256((__m256i *)(aboveBy32 + base + 56));
        b = _mm256_loadu_si256((__m256i *)(aboveDiff + base + 56));
        b = _mm256_mullo_epi32(b, shift);
        res4 = _mm256_add_epi32(a, b);
        res4 = _mm256_srli_epi32(res4, 5);

        x += dx;
        shift = _mm256_srli_epi32(_mm256_and_si256(_mm256_set1_epi32(x), _mm256_set1_epi32(0x3f)), 1);

        res1 = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2), PERM4x64(0, 2, 1, 3));
        res2 = _mm256_permute4x64_epi64(_mm256_packus_epi32(res3, res4), PERM4x64(0, 2, 1, 3));
        resLo = _mm256_min_epi16(clip_bd, resLo);
        _mm256_storeu_si256((__m256i *)dst, resLo);
        resHi = _mm256_min_epi16(clip_bd, resHi);
        _mm256_storeu_si256((__m256i *)(dst + 16), resHi);
        res1 = _mm256_min_epi16(clip_bd, res1);
        _mm256_storeu_si256((__m256i *)(dst + 32), res1);
        res2 = _mm256_min_epi16(clip_bd, res2);
        _mm256_storeu_si256((__m256i *)(dst + 48), res2);
    }
}

void intra_mode_angular_av1_z1_16bit_4x4_avx2(const uint32_t size, uint16_t *ref_samples, uint16_t *dst, const uint32_t prediction_buffer_stride, const EbBool skip, uint16_t dx, uint16_t dy, uint16_t bd)
{
    (void)dy;
    uint32_t toAboveOffset = (size << 1) + 1;
    uint32_t rowStride = skip ? 2 : 1;
    Predict_Z1_16bit_4x4(&ref_samples[toAboveOffset], dst, rowStride* prediction_buffer_stride, dx, bd);
}
void intra_mode_angular_av1_z1_16bit_8x8_avx2(const uint32_t size, uint16_t *ref_samples, uint16_t *dst, const uint32_t prediction_buffer_stride, const EbBool skip, uint16_t dx, uint16_t dy, uint16_t bd)
{
    (void)dy;
    uint32_t toAboveOffset = (size << 1) + 1;
    uint32_t rowStride = skip ? 2 : 1;
    Predict_Z1_16bit_8x8(&ref_samples[toAboveOffset], dst, rowStride* prediction_buffer_stride, dx, bd);
}
void intra_mode_angular_av1_z1_16bit_16x16_avx2(const uint32_t size, uint16_t *ref_samples, uint16_t *dst, const uint32_t prediction_buffer_stride, const EbBool skip, uint16_t dx, uint16_t dy, uint16_t bd)
{
    (void)dy;
    uint32_t toAboveOffset = (size << 1) + 1;
    uint32_t rowStride = skip ? 2 : 1;
    Predict_Z1_16bit_16x16(&ref_samples[toAboveOffset], dst, rowStride* prediction_buffer_stride, dx, bd);
}
void intra_mode_angular_av1_z1_16bit_32x32_avx2(const uint32_t size, uint16_t *ref_samples, uint16_t *dst, const uint32_t prediction_buffer_stride, const EbBool skip, uint16_t dx, uint16_t dy, uint16_t bd)
{
    (void)dy;
    uint32_t toAboveOffset = (size << 1) + 1;
    uint32_t rowStride = skip ? 2 : 1;
    Predict_Z1_16bit_32x32(&ref_samples[toAboveOffset], dst, rowStride* prediction_buffer_stride, dx, bd);
}
void intra_mode_angular_av1_z1_16bit_64x64_avx2(const uint32_t size, uint16_t *ref_samples, uint16_t *dst, const uint32_t prediction_buffer_stride, const EbBool skip, uint16_t dx, uint16_t dy, uint16_t bd)
{
    (void)dy;
    uint32_t toAboveOffset = (size << 1) + 1;
    uint32_t rowStride = skip ? 2 : 1;
    Predict_Z1_16bit_64x64(&ref_samples[toAboveOffset], dst, rowStride* prediction_buffer_stride, dx, bd);
}

void transpose_16bit_TX_4X4(const uint16_t *src, uint32_t srcStride, uint16_t *dst, uint32_t dstStride)
{
    assert(srcStride == 4);
    (void)srcStride;

    if (dstStride == 4)
    {
        __m128i s = _mm_loadu_si128((__m128i*)src);
        __m128i r1 = _mm_srli_si128(s, 8);
        __m128i r2 = _mm_loadu_si128((__m128i*)(src + 8));
        __m128i r3 = _mm_srli_si128(r2, 8);

        __m128i r0_Lo = _mm_unpacklo_epi16(s, r1);
        __m128i r2_Lo = _mm_unpacklo_epi16(r2, r3);
        __m128i r1_Lo = _mm_unpacklo_epi32(r0_Lo, r2_Lo);
        r0_Lo = _mm_unpackhi_epi32(r0_Lo, r2_Lo);

        _mm_storeu_si128((__m128i*)(dst + 0 * dstStride), r1_Lo);
        _mm_storeu_si128((__m128i*)(dst + 2 * dstStride), r0_Lo);
    }
    else
    {
        __m128i s = _mm_loadu_si128((__m128i*)src);
        __m128i r1 = _mm_srli_si128(s, 8);
        __m128i r2 = _mm_loadu_si128((__m128i*)(src + 8));
        __m128i r3 = _mm_srli_si128(r2, 8);

        __m128i r0_Lo = _mm_unpacklo_epi16(s, r1);
        __m128i r2_Lo = _mm_unpacklo_epi16(r2, r3);
        __m128i r1_Lo = _mm_unpacklo_epi32(r0_Lo, r2_Lo);
        r0_Lo = _mm_unpackhi_epi32(r0_Lo, r2_Lo);

        _mm_storel_epi64((__m128i*)(dst + 0 * dstStride), r1_Lo);
        _mm_storel_epi64((__m128i*)(dst + 1 * dstStride), _mm_srli_si128(r1_Lo, 8));
        _mm_storel_epi64((__m128i*)(dst + 2 * dstStride), r0_Lo);
        _mm_storel_epi64((__m128i*)(dst + 3 * dstStride), _mm_srli_si128(r0_Lo, 8));
    }
}
void transpose_16bit_TX_8X8(const uint16_t *src, uint32_t srcStride, uint16_t *dst, uint32_t dstStride)
{
    __m128i r0, r1, r2, r3, r4, r5, r6, r7, r0_Lo, r1_Lo, r2_Lo, r3_Lo, r4_Lo, r5_Lo, r6_Lo;
    r0 = _mm_loadu_si128((__m128i*)(src + 0 * srcStride));   // 07,06,05,04,03,02,01,00
    r1 = _mm_loadu_si128((__m128i*)(src + 1 * srcStride));   // 17,16,15,14,13,12,11,10
    r2 = _mm_loadu_si128((__m128i*)(src + 2 * srcStride));   // 27,26,25,24,23,22,21,20
    r3 = _mm_loadu_si128((__m128i*)(src + 3 * srcStride));   // 37,36,35,34,33,32,31,30
    r4 = _mm_loadu_si128((__m128i*)(src + 4 * srcStride));   // 47,46,45,44,43,42,41,40
    r5 = _mm_loadu_si128((__m128i*)(src + 5 * srcStride));   // 57,56,55,54,53,52,51,50
    r6 = _mm_loadu_si128((__m128i*)(src + 6 * srcStride));   // 67,66,65,64,63,62,61,60
    r7 = _mm_loadu_si128((__m128i*)(src + 7 * srcStride));   // 77,76,75,74,73,72,71,70

    r0_Lo = _mm_unpacklo_epi16(r0, r1);
    r2_Lo = _mm_unpacklo_epi16(r2, r3);
    r4_Lo = _mm_unpacklo_epi16(r4, r5);
    r6_Lo = _mm_unpacklo_epi16(r6, r7);

    r1_Lo = r0_Lo;
    r0_Lo = _mm_unpacklo_epi32(r0_Lo, r2_Lo);
    r1_Lo = _mm_unpackhi_epi32(r1_Lo, r2_Lo);
    r5_Lo = r4_Lo;
    r4_Lo = _mm_unpacklo_epi32(r4_Lo, r6_Lo);
    r5_Lo = _mm_unpackhi_epi32(r5_Lo, r6_Lo);
    r2_Lo = r0_Lo;
    r0_Lo = _mm_unpacklo_epi64(r0_Lo, r4_Lo); //64
    r2_Lo = _mm_unpackhi_epi64(r2_Lo, r4_Lo);
    r3_Lo = r1_Lo;
    r1_Lo = _mm_unpacklo_epi64(r1_Lo, r5_Lo);
    r3_Lo = _mm_unpackhi_epi64(r3_Lo, r5_Lo);

    _mm_storeu_si128((__m128i*)(dst + 0 * dstStride), r0_Lo);
    _mm_storeu_si128((__m128i*)(dst + 1 * dstStride), r2_Lo);
    _mm_storeu_si128((__m128i*)(dst + 2 * dstStride), r1_Lo);
    _mm_storeu_si128((__m128i*)(dst + 3 * dstStride), r3_Lo);

    r0 = _mm_unpackhi_epi16(r0, r1);
    r2 = _mm_unpackhi_epi16(r2, r3);
    r4 = _mm_unpackhi_epi16(r4, r5);
    r6 = _mm_unpackhi_epi16(r6, r7);

    r1 = r0;
    r0 = _mm_unpacklo_epi32(r0, r2);
    r1 = _mm_unpackhi_epi32(r1, r2);
    r5 = r4;
    r4 = _mm_unpacklo_epi32(r4, r6);
    r5 = _mm_unpackhi_epi32(r5, r6);
    r2 = r0;
    r0 = _mm_unpacklo_epi64(r0, r4);
    r2 = _mm_unpackhi_epi64(r2, r4);
    r3 = r1;
    r1 = _mm_unpacklo_epi64(r1, r5);
    r3 = _mm_unpackhi_epi64(r3, r5);

    _mm_storeu_si128((__m128i*)(dst + 4 * dstStride), r0);
    _mm_storeu_si128((__m128i*)(dst + 5 * dstStride), r2);
    _mm_storeu_si128((__m128i*)(dst + 6 * dstStride), r1);
    _mm_storeu_si128((__m128i*)(dst + 7 * dstStride), r3);
}
void transpose_16bit(const uint16_t *src, uint32_t srcStride, uint16_t *dst, uint32_t dstStride, int32_t width, int32_t height)
{
    for (int32_t j = 0; j < height; j += 8)
        for (int32_t i = 0; i < width; i += 8)
            transpose_16bit_TX_8X8(src + i * srcStride + j, srcStride, dst + j * dstStride + i, dstStride);
}
uint16_t z2BlendMaskabove_16bit[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };
uint16_t z2BlendMaskleft_16bit[] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
void Predict_Z2_16bit_4x4(uint16_t *refPel, uint16_t *dst, const uint32_t pitch, uint16_t d)
{
    uint32_t leftBy32[8];
    int_least32_t leftByDiff[8];
    int32_t c, x, y, base;
    const uint16_t* leftPtr = refPel - 1;
    const int32_t frac_bits = 6;
    __m128i a0, a1, diff, a32, a16;
    __m128i a, b, res1, res2, resLo;
    res2 = _mm_setzero_si128();

    a16 = _mm_set1_epi32(16);
    a0 = _mm_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));        // 01234567
    a1 = _mm_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 1)));    // 89abcdef
    diff = _mm_sub_epi32(a1, a0);                                        // a[x+1] - a[x]
    a32 = _mm_slli_epi32(a0, 5);                                        // a[x] * 32
    a32 = _mm_add_epi32(a32, a16);                                        // a[x] * 32 + 16
    _mm_storeu_si128((__m128i *)(leftBy32 + 4), a32);
    _mm_storeu_si128((__m128i *)(leftByDiff + 4), diff);

    for (c = 0; c < 4; ++c) {
        x = c + 1;
        y = -x * d;
        base = y >> frac_bits;

        __m128i inc = _mm_set1_epi32(y);
        __m128i shift = _mm_srli_epi32(_mm_and_si128(inc, _mm_set1_epi32(0x3f)), 1);

        base = base + 1;
        a = _mm_loadu_si128((__m128i *)(leftBy32 + base + 4));
        b = _mm_loadu_si128((__m128i *)(leftByDiff + base + 4));
        b = _mm_mullo_epi32(b, shift);
        res1 = _mm_add_epi32(a, b);
        res1 = _mm_srli_epi32(res1, 5);

        resLo = _mm_packus_epi32(res1, res2);
        _mm_storel_epi64((__m128i *)(dst + c * pitch), resLo);
    }
}
void Predict_Z2_16bit_4x4_Left(uint16_t *refPel, uint16_t *dst, const uint32_t pitch, uint16_t d, uint32_t toAboveLeftOffset)
{
    uint32_t leftBy32[8];
    int_least32_t leftByDiff[8];
    uint32_t leftBy32PlusOffset[8];
    int_least32_t leftByDiffPlusOffset[8];
    int32_t c, x, y, base, counter;
    const uint16_t* leftPtr = refPel - 1;
    const int32_t frac_bits = 6;
    __m128i a0, a1, diff, a32, a16;
    __m128i a, b, b1, res1, res2, resLo, leftPlusOffset;
    res2 = _mm_setzero_si128();

    a16 = _mm_set1_epi32(16);
    a0 = _mm_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));        // 01234567
    a1 = _mm_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 1)));    // 89abcdef
    diff = _mm_sub_epi32(a1, a0);                                        // a[x+1] - a[x]
    a32 = _mm_slli_epi32(a0, 5);                                        // a[x] * 32
    a32 = _mm_add_epi32(a32, a16);                                        // a[x] * 32 + 16
    _mm_storeu_si128((__m128i *)(leftBy32 + 4), a32);
    _mm_storeu_si128((__m128i *)(leftByDiff + 4), diff);

    leftPlusOffset = _mm_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 1 + toAboveLeftOffset)));
    diff = _mm_sub_epi32(a1, leftPlusOffset);
    a32 = _mm_slli_epi32(leftPlusOffset, 5);
    a32 = _mm_add_epi32(a32, a16);
    _mm_storeu_si128((__m128i *)(leftBy32PlusOffset + 4), a32);
    _mm_storeu_si128((__m128i *)(leftByDiffPlusOffset + 4), diff);

    for (c = 0; c < 4; ++c) {
        x = c + 1;
        y = -x * d;
        base = y >> frac_bits;

        __m128i inc = _mm_set1_epi32(y);
        __m128i shift = _mm_srli_epi32(_mm_and_si128(inc, _mm_set1_epi32(0x3f)), 1);

        if (base <= -1)
        {
            base = base + 1;
            a1 = _mm_loadu_si128((__m128i *)(leftBy32PlusOffset + base + 4));
            b1 = _mm_loadu_si128((__m128i *)(leftByDiffPlusOffset + base + 4));
            a = _mm_loadu_si128((__m128i *)(leftBy32 + base + 4));
            b = _mm_loadu_si128((__m128i *)(leftByDiff + base + 4));
            for (counter = 0; counter < (-base + 1); counter++)
            {
                if ((counter & 3) == 0)
                {
                    a = _mm_blend_epi32(a, a1, 1);
                    b = _mm_blend_epi32(b, b1, 1);
                }
                else if ((counter & 3) == 1)
                {
                    a = _mm_blend_epi32(a, a1, 2);
                    b = _mm_blend_epi32(b, b1, 2);
                }
                else if ((counter & 3) == 2)
                {
                    a = _mm_blend_epi32(a, a1, 4);
                    b = _mm_blend_epi32(b, b1, 4);
                }
                else
                {
                    a = _mm_blend_epi32(a, a1, 8);
                    b = _mm_blend_epi32(b, b1, 8);
                }

            }
            b = _mm_mullo_epi32(b, shift);
            res1 = _mm_add_epi32(a, b);
            res1 = _mm_srli_epi32(res1, 5);

            resLo = _mm_packus_epi32(res1, res2);
            _mm_storel_epi64((__m128i *)(dst + c * pitch), resLo);
        }
        else
        {
            base = base + 1;
            a = _mm_loadu_si128((__m128i *)(leftBy32 + base + 4));
            b = _mm_loadu_si128((__m128i *)(leftByDiff + base + 4));
            b = _mm_mullo_epi32(b, shift);
            res1 = _mm_add_epi32(a, b);
            res1 = _mm_srli_epi32(res1, 5);

            resLo = _mm_packus_epi32(res1, res2);
            _mm_storel_epi64((__m128i *)(dst + c * pitch), resLo);
        }
    }
}
void Predict_Z2_16bit_8x8(uint16_t *refPel, uint16_t *dst, const uint32_t pitch, uint16_t d)
{
    __m256i a0, a1, diff, a32, a16;
    __m256i a, b, res1, res2, resLo;
    int32_t c, x, y, base;
    const uint16_t* leftPtr = refPel - 1;
    const int32_t frac_bits = 6;

    res2 = _mm256_setzero_si256();
    uint32_t leftBy32[16];
    int_least32_t leftByDiff[16];

    a16 = _mm256_set1_epi32(16);
    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));        // 01234567
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 1)));    // 89abcdef
    diff = _mm256_sub_epi32(a1, a0);                                        // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0, 5);                                            // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);                                        // a[x] * 32 + 16
    _mm256_storeu_si256((__m256i *)(leftBy32 + 8), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 8), diff);

    for (c = 0; c < 8; ++c) {
        x = c + 1;
        y = -x * d;
        base = y >> frac_bits;

        __m256i inc = _mm256_set1_epi32(y);
        __m256i shift = _mm256_srli_epi32(_mm256_and_si256(inc, _mm256_set1_epi32(0x3f)), 1);

        base = base + 1;
        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 8));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 8));
        b = _mm256_mullo_epi32(b, shift);
        res1 = _mm256_add_epi32(a, b);
        res1 = _mm256_srli_epi32(res1, 5);

        resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2), PERM4x64(0, 2, 1, 3));
        _mm_storeu_si128((__m128i *)(dst + c * pitch), _mm256_castsi256_si128(resLo));
    }

}
void Predict_Z2_16bit_8x8_Left(uint16_t *refPel, uint16_t *dst, const uint32_t pitch, uint16_t d, uint32_t toAboveLeftOffset)
{
    __m256i a0, a1, diff, a32, a16;
    __m256i a, b, b1, res1, res2, resLo, leftPlusOffset;
    int32_t c, x, y, base, counter;
    const uint16_t* leftPtr = refPel - 1;
    const int32_t frac_bits = 6;

    res2 = _mm256_setzero_si256();
    uint32_t leftBy32[16];
    int_least32_t leftByDiff[16];
    uint32_t leftBy32PlusOffset[16];
    int_least32_t leftByDiffPlusOffset[16];

    a16 = _mm256_set1_epi32(16);
    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));        // 01234567
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 1)));    // 89abcdef
    diff = _mm256_sub_epi32(a1, a0);                                        // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0, 5);                                            // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);                                        // a[x] * 32 + 16
    _mm256_storeu_si256((__m256i *)(leftBy32 + 8), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 8), diff);

    leftPlusOffset = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 1 + toAboveLeftOffset)));
    diff = _mm256_sub_epi32(a1, leftPlusOffset);
    a32 = _mm256_slli_epi32(leftPlusOffset, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32PlusOffset + 8), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiffPlusOffset + 8), diff);

    for (c = 0; c < 8; ++c) {
        x = c + 1;
        y = -x * d;
        base = y >> frac_bits;

        __m256i inc = _mm256_set1_epi32(y);
        __m256i shift = _mm256_srli_epi32(_mm256_and_si256(inc, _mm256_set1_epi32(0x3f)), 1);

        if (base <= -1)
        {
            base = base + 1;
            a1 = _mm256_loadu_si256((__m256i *)(leftBy32PlusOffset + base + 8));
            b1 = _mm256_loadu_si256((__m256i *)(leftByDiffPlusOffset + base + 8));
            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 8));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 8));
            for (counter = 0; counter < (-base + 1); counter++)
            {
                if ((counter & 7) == 0)
                {
                    a = _mm256_blend_epi32(a, a1, 1);
                    b = _mm256_blend_epi32(b, b1, 1);
                }
                else if ((counter & 7) == 1)
                {
                    a = _mm256_blend_epi32(a, a1, 2);
                    b = _mm256_blend_epi32(b, b1, 2);
                }
                else if ((counter & 7) == 2)
                {
                    a = _mm256_blend_epi32(a, a1, 4);
                    b = _mm256_blend_epi32(b, b1, 4);
                }
                else if ((counter & 7) == 3)
                {
                    a = _mm256_blend_epi32(a, a1, 8);
                    b = _mm256_blend_epi32(b, b1, 8);
                }
                else if ((counter & 7) == 4)
                {
                    a = _mm256_blend_epi32(a, a1, 16);
                    b = _mm256_blend_epi32(b, b1, 16);
                }
                else if ((counter & 7) == 5)
                {
                    a = _mm256_blend_epi32(a, a1, 32);
                    b = _mm256_blend_epi32(b, b1, 32);
                }
                else if ((counter & 7) == 6)
                {
                    a = _mm256_blend_epi32(a, a1, 64);
                    b = _mm256_blend_epi32(b, b1, 64);
                }
                else
                {
                    a = _mm256_blend_epi32(a, a1, 128);
                    b = _mm256_blend_epi32(b, b1, 128);
                }

            }
            b = _mm256_mullo_epi32(b, shift);
            res1 = _mm256_add_epi32(a, b);
            res1 = _mm256_srli_epi32(res1, 5);

            resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2), PERM4x64(0, 2, 1, 3));
            _mm_storeu_si128((__m128i *)(dst + c * pitch), _mm256_castsi256_si128(resLo));
        }
        else
        {
            base = base + 1;
            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 8));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 8));
            b = _mm256_mullo_epi32(b, shift);
            res1 = _mm256_add_epi32(a, b);
            res1 = _mm256_srli_epi32(res1, 5);

            resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2), PERM4x64(0, 2, 1, 3));
            _mm_storeu_si128((__m128i *)(dst + c * pitch), _mm256_castsi256_si128(resLo));
        }
    }

}
void Predict_Z2_16bit_16x16(uint16_t *refPel, uint16_t *dst, const uint32_t pitch, uint16_t d)
{
    uint32_t leftBy32[32];
    int_least32_t leftByDiff[32];
    __m256i a0, a1, diff, a32, a16;
    int32_t c, x, y, base;
    const uint16_t* leftPtr = refPel - 1;
    const int32_t frac_bits = 6;

    a16 = _mm256_set1_epi32(16);
    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));        // 01234567
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 1)));    // 89abcdef
    diff = _mm256_sub_epi32(a1, a0);                                        // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0, 5);                                            // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);                                        // a[x] * 32 + 16
    _mm256_storeu_si256((__m256i *)(leftBy32 + 16), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 16), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 8)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 9)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 24), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 24), diff);

    for (c = 0; c < 16; ++c) {
        x = c + 1;
        y = -x * d;
        base = y >> frac_bits;

        __m256i inc = _mm256_set1_epi32(y);
        __m256i shift = _mm256_srli_epi32(_mm256_and_si256(inc, _mm256_set1_epi32(0x3f)), 1);

        __m256i a, b, res1, res2, resLo;
        base = base + 1;
        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 16));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 16));
        b = _mm256_mullo_epi32(b, shift);
        res1 = _mm256_add_epi32(a, b);
        res1 = _mm256_srli_epi32(res1, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 24));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 24));
        b = _mm256_mullo_epi32(b, shift);
        res2 = _mm256_add_epi32(a, b);
        res2 = _mm256_srli_epi32(res2, 5);

        resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2), PERM4x64(0, 2, 1, 3));
        _mm256_storeu_si256((__m256i *)(dst + c * pitch), resLo);
    }


}
void Predict_Z2_16bit_16x16_Left(uint16_t *refPel, uint16_t *dst, const uint32_t pitch, uint16_t d, uint32_t toAboveLeftOffset)
{
    uint32_t leftBy32[32];
    int_least32_t leftByDiff[32];
    uint32_t leftBy32PlusOffset[32];
    int_least32_t leftByDiffPlusOffset[32];
    __m256i a0, a1, diff, a32, a16, leftPlusOffset;
    int32_t c, x, y, base, counter;
    const uint16_t* leftPtr = refPel - 1;
    const int32_t frac_bits = 6;

    a16 = _mm256_set1_epi32(16);
    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));        // 01234567
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 1)));    // 89abcdef
    diff = _mm256_sub_epi32(a1, a0);                                        // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0, 5);                                            // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);                                        // a[x] * 32 + 16
    _mm256_storeu_si256((__m256i *)(leftBy32 + 16), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 16), diff);

    leftPlusOffset = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 1 + toAboveLeftOffset)));
    diff = _mm256_sub_epi32(a1, leftPlusOffset);
    a32 = _mm256_slli_epi32(leftPlusOffset, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32PlusOffset + 16), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiffPlusOffset + 16), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 8)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 9)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 24), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 24), diff);

    diff = _mm256_sub_epi32(a1, leftPlusOffset);
    a32 = _mm256_slli_epi32(leftPlusOffset, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32PlusOffset + 24), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiffPlusOffset + 24), diff);

    for (c = 0; c < 16; ++c) {
        x = c + 1;
        y = -x * d;
        base = y >> frac_bits;

        __m256i inc = _mm256_set1_epi32(y);
        __m256i shift = _mm256_srli_epi32(_mm256_and_si256(inc, _mm256_set1_epi32(0x3f)), 1);

        __m256i a, b, a1, b1, res1, res2, resLo;

        if (base <= -1)
        {
            base = base + 1;
            a1 = _mm256_loadu_si256((__m256i *)(leftBy32PlusOffset + base + 16));
            b1 = _mm256_loadu_si256((__m256i *)(leftByDiffPlusOffset + base + 16));
            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 16));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 16));
            for (counter = 0; counter < (-base + 1); counter++)
            {
                if ((counter & 7) == 0)
                {
                    a = _mm256_blend_epi32(a, a1, 1);
                    b = _mm256_blend_epi32(b, b1, 1);
                }
                else if ((counter & 7) == 1)
                {
                    a = _mm256_blend_epi32(a, a1, 2);
                    b = _mm256_blend_epi32(b, b1, 2);
                }
                else if ((counter & 7) == 2)
                {
                    a = _mm256_blend_epi32(a, a1, 4);
                    b = _mm256_blend_epi32(b, b1, 4);
                }
                else if ((counter & 7) == 3)
                {
                    a = _mm256_blend_epi32(a, a1, 8);
                    b = _mm256_blend_epi32(b, b1, 8);
                }
                else if ((counter & 7) == 4)
                {
                    a = _mm256_blend_epi32(a, a1, 16);
                    b = _mm256_blend_epi32(b, b1, 16);
                }
                else if ((counter & 7) == 5)
                {
                    a = _mm256_blend_epi32(a, a1, 32);
                    b = _mm256_blend_epi32(b, b1, 32);
                }
                else if ((counter & 7) == 6)
                {
                    a = _mm256_blend_epi32(a, a1, 64);
                    b = _mm256_blend_epi32(b, b1, 64);
                }
                else
                {
                    a = _mm256_blend_epi32(a, a1, 128);
                    b = _mm256_blend_epi32(b, b1, 128);
                }
            }
            b = _mm256_mullo_epi32(b, shift);
            res1 = _mm256_add_epi32(a, b);
            res1 = _mm256_srli_epi32(res1, 5);

            a1 = _mm256_loadu_si256((__m256i *)(leftBy32PlusOffset + base + 24));
            b1 = _mm256_loadu_si256((__m256i *)(leftByDiffPlusOffset + base + 24));
            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 24));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 24));
            for (counter = 8; counter < (-base + 1); counter++)
            {
                if ((counter & 7) == 0)
                {
                    a = _mm256_blend_epi32(a, a1, 1);
                    b = _mm256_blend_epi32(b, b1, 1);
                }
                else if ((counter & 7) == 1)
                {
                    a = _mm256_blend_epi32(a, a1, 2);
                    b = _mm256_blend_epi32(b, b1, 2);
                }
                else if ((counter & 7) == 2)
                {
                    a = _mm256_blend_epi32(a, a1, 4);
                    b = _mm256_blend_epi32(b, b1, 4);
                }
                else if ((counter & 7) == 3)
                {
                    a = _mm256_blend_epi32(a, a1, 8);
                    b = _mm256_blend_epi32(b, b1, 8);
                }
                else if ((counter & 7) == 4)
                {
                    a = _mm256_blend_epi32(a, a1, 16);
                    b = _mm256_blend_epi32(b, b1, 16);
                }
                else if ((counter & 7) == 5)
                {
                    a = _mm256_blend_epi32(a, a1, 32);
                    b = _mm256_blend_epi32(b, b1, 32);
                }
                else if ((counter & 7) == 6)
                {
                    a = _mm256_blend_epi32(a, a1, 64);
                    b = _mm256_blend_epi32(b, b1, 64);
                }
                else
                {
                    a = _mm256_blend_epi32(a, a1, 128);
                    b = _mm256_blend_epi32(b, b1, 128);
                }
            }
            b = _mm256_mullo_epi32(b, shift);
            res2 = _mm256_add_epi32(a, b);
            res2 = _mm256_srli_epi32(res2, 5);

            resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2), PERM4x64(0, 2, 1, 3));
            _mm256_storeu_si256((__m256i *)(dst + c * pitch), resLo);
        }
        else
        {
            base = base + 1;
            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 16));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 16));
            b = _mm256_mullo_epi32(b, shift);
            res1 = _mm256_add_epi32(a, b);
            res1 = _mm256_srli_epi32(res1, 5);

            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 24));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 24));
            b = _mm256_mullo_epi32(b, shift);
            res2 = _mm256_add_epi32(a, b);
            res2 = _mm256_srli_epi32(res2, 5);

            resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2), PERM4x64(0, 2, 1, 3));
            _mm256_storeu_si256((__m256i *)(dst + c * pitch), resLo);
        }
    }


}
void Predict_Z2_16bit_32x32(uint16_t *refPel, uint16_t *dst, const uint32_t pitch, uint16_t d)
{
    uint32_t leftBy32[72];
    int_least32_t leftByDiff[72];
    __m256i a0, a1, diff, a32, a16;
    int32_t c, x, y, base;
    const uint16_t* leftPtr = refPel - 1;
    const int32_t frac_bits = 6;

    a16 = _mm256_set1_epi32(16);
    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));      // 01234567
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 1)));    // 89abcdef
    diff = _mm256_sub_epi32(a1, a0);                                        // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0, 5);                                            // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);                                        // a[x] * 32 + 16
    _mm256_storeu_si256((__m256i *)(leftBy32 + 32), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 32), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 8)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 9)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 40), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 40), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 16)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 17)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 48), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 48), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 24)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 25)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 56), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 56), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 32)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 33)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 64), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 64), diff);


    for (c = 0; c < 32; ++c) {
        x = c + 1;
        y = -x * d;
        base = y >> frac_bits;

        __m256i inc = _mm256_set1_epi32(y);
        __m256i shift = _mm256_srli_epi32(_mm256_and_si256(inc, _mm256_set1_epi32(0x3f)), 1);

        __m256i a, b, res1, res2, res3, res4, resLo, resHi;

        base = base + 1;
        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 32));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 32));
        b = _mm256_mullo_epi32(b, shift);
        res1 = _mm256_add_epi32(a, b);
        res1 = _mm256_srli_epi32(res1, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 40));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 40));
        b = _mm256_mullo_epi32(b, shift);
        res2 = _mm256_add_epi32(a, b);
        res2 = _mm256_srli_epi32(res2, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 48));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 48));
        b = _mm256_mullo_epi32(b, shift);
        res3 = _mm256_add_epi32(a, b);
        res3 = _mm256_srli_epi32(res3, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 56));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 56));
        b = _mm256_mullo_epi32(b, shift);
        res4 = _mm256_add_epi32(a, b);
        res4 = _mm256_srli_epi32(res4, 5);

        resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2), PERM4x64(0, 2, 1, 3));
        resHi = _mm256_permute4x64_epi64(_mm256_packus_epi32(res3, res4), PERM4x64(0, 2, 1, 3));
        _mm256_storeu_si256((__m256i *)(dst + c * pitch), resLo);
        _mm256_storeu_si256((__m256i *)(dst + 16 + c * pitch), resHi);

    }


}
void Predict_Z2_16bit_32x32_Left(uint16_t *refPel, uint16_t *dst, const uint32_t pitch, uint16_t d, uint32_t toAboveLeftOffset)
{
    uint32_t leftBy32[64];
    int_least32_t leftByDiff[64];
    uint32_t leftBy32PlusOffset[64];
    int_least32_t leftByDiffPlusOffset[64];
    __m256i a0, a1, diff, a32, a16, leftPlusOffset;
    int32_t c, x, y, base;
    const uint16_t* leftPtr = refPel - 1;
    const int32_t frac_bits = 6;

    a16 = _mm256_set1_epi32(16);
    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));      // 01234567
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 1)));    // 89abcdef
    diff = _mm256_sub_epi32(a1, a0);                                        // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0, 5);                                            // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);                                        // a[x] * 32 + 16
    _mm256_storeu_si256((__m256i *)(leftBy32 + 32), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 32), diff);

    leftPlusOffset = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 1 + toAboveLeftOffset)));
    diff = _mm256_sub_epi32(a1, leftPlusOffset);
    a32 = _mm256_slli_epi32(leftPlusOffset, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32PlusOffset + 32), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiffPlusOffset + 32), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 8)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 9)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 40), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 40), diff);

    diff = _mm256_sub_epi32(a1, leftPlusOffset);
    a32 = _mm256_slli_epi32(leftPlusOffset, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32PlusOffset + 40), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiffPlusOffset + 40), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 16)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 17)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 48), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 48), diff);

    diff = _mm256_sub_epi32(a1, leftPlusOffset);
    a32 = _mm256_slli_epi32(leftPlusOffset, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32PlusOffset + 48), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiffPlusOffset + 48), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 24)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 25)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 56), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 56), diff);

    diff = _mm256_sub_epi32(a1, leftPlusOffset);
    a32 = _mm256_slli_epi32(leftPlusOffset, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32PlusOffset + 56), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiffPlusOffset + 56), diff);

    int32_t counter;
    for (c = 0; c < 32; ++c) {
        x = c + 1;
        y = -x * d;
        base = y >> frac_bits;

        __m256i inc = _mm256_set1_epi32(y);
        __m256i shift = _mm256_srli_epi32(_mm256_and_si256(inc, _mm256_set1_epi32(0x3f)), 1);

        __m256i a, b, a1, b1, res1, res2, res3, res4, resLo, resHi;

        if (base <= -1)
        {
            base = base + 1;
            a1 = _mm256_loadu_si256((__m256i *)(leftBy32PlusOffset + base + 32));
            b1 = _mm256_loadu_si256((__m256i *)(leftByDiffPlusOffset + base + 32));
            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 32));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 32));
            for (counter = 0; counter < (-base + 1); counter++)
            {
                if ((counter & 7) == 0)
                {
                    a = _mm256_blend_epi32(a, a1, 1);
                    b = _mm256_blend_epi32(b, b1, 1);
                }
                else if ((counter & 7) == 1)
                {
                    a = _mm256_blend_epi32(a, a1, 2);
                    b = _mm256_blend_epi32(b, b1, 2);
                }
                else if ((counter & 7) == 2)
                {
                    a = _mm256_blend_epi32(a, a1, 4);
                    b = _mm256_blend_epi32(b, b1, 4);
                }
                else if ((counter & 7) == 3)
                {
                    a = _mm256_blend_epi32(a, a1, 8);
                    b = _mm256_blend_epi32(b, b1, 8);
                }
                else if ((counter & 7) == 4)
                {
                    a = _mm256_blend_epi32(a, a1, 16);
                    b = _mm256_blend_epi32(b, b1, 16);
                }
                else if ((counter & 7) == 5)
                {
                    a = _mm256_blend_epi32(a, a1, 32);
                    b = _mm256_blend_epi32(b, b1, 32);
                }
                else if ((counter & 7) == 6)
                {
                    a = _mm256_blend_epi32(a, a1, 64);
                    b = _mm256_blend_epi32(b, b1, 64);
                }
                else
                {
                    a = _mm256_blend_epi32(a, a1, 128);
                    b = _mm256_blend_epi32(b, b1, 128);
                }
            }
            b = _mm256_mullo_epi32(b, shift);
            res1 = _mm256_add_epi32(a, b);
            res1 = _mm256_srli_epi32(res1, 5);

            a1 = _mm256_loadu_si256((__m256i *)(leftBy32PlusOffset + base + 40));
            b1 = _mm256_loadu_si256((__m256i *)(leftByDiffPlusOffset + base + 40));
            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 40));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 40));
            for (counter = 8; counter < (-base + 1); counter++)
            {
                if ((counter & 7) == 0)
                {
                    a = _mm256_blend_epi32(a, a1, 1);
                    b = _mm256_blend_epi32(b, b1, 1);
                }
                else if ((counter & 7) == 1)
                {
                    a = _mm256_blend_epi32(a, a1, 2);
                    b = _mm256_blend_epi32(b, b1, 2);
                }
                else if ((counter & 7) == 2)
                {
                    a = _mm256_blend_epi32(a, a1, 4);
                    b = _mm256_blend_epi32(b, b1, 4);
                }
                else if ((counter & 7) == 3)
                {
                    a = _mm256_blend_epi32(a, a1, 8);
                    b = _mm256_blend_epi32(b, b1, 8);
                }
                else if ((counter & 7) == 4)
                {
                    a = _mm256_blend_epi32(a, a1, 16);
                    b = _mm256_blend_epi32(b, b1, 16);
                }
                else if ((counter & 7) == 5)
                {
                    a = _mm256_blend_epi32(a, a1, 32);
                    b = _mm256_blend_epi32(b, b1, 32);
                }
                else if ((counter & 7) == 6)
                {
                    a = _mm256_blend_epi32(a, a1, 64);
                    b = _mm256_blend_epi32(b, b1, 64);
                }
                else
                {
                    a = _mm256_blend_epi32(a, a1, 128);
                    b = _mm256_blend_epi32(b, b1, 128);
                }
            }
            b = _mm256_mullo_epi32(b, shift);
            res2 = _mm256_add_epi32(a, b);
            res2 = _mm256_srli_epi32(res2, 5);


            a1 = _mm256_loadu_si256((__m256i *)(leftBy32PlusOffset + base + 48));
            b1 = _mm256_loadu_si256((__m256i *)(leftByDiffPlusOffset + base + 48));
            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 48));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 48));
            for (counter = 16; counter < (-base + 1); counter++)
            {
                if ((counter & 7) == 0)
                {
                    a = _mm256_blend_epi32(a, a1, 1);
                    b = _mm256_blend_epi32(b, b1, 1);
                }
                else if ((counter & 7) == 1)
                {
                    a = _mm256_blend_epi32(a, a1, 2);
                    b = _mm256_blend_epi32(b, b1, 2);
                }
                else if ((counter & 7) == 2)
                {
                    a = _mm256_blend_epi32(a, a1, 4);
                    b = _mm256_blend_epi32(b, b1, 4);
                }
                else if ((counter & 7) == 3)
                {
                    a = _mm256_blend_epi32(a, a1, 8);
                    b = _mm256_blend_epi32(b, b1, 8);
                }
                else if ((counter & 7) == 4)
                {
                    a = _mm256_blend_epi32(a, a1, 16);
                    b = _mm256_blend_epi32(b, b1, 16);
                }
                else if ((counter & 7) == 5)
                {
                    a = _mm256_blend_epi32(a, a1, 32);
                    b = _mm256_blend_epi32(b, b1, 32);
                }
                else if ((counter & 7) == 6)
                {
                    a = _mm256_blend_epi32(a, a1, 64);
                    b = _mm256_blend_epi32(b, b1, 64);
                }
                else
                {
                    a = _mm256_blend_epi32(a, a1, 128);
                    b = _mm256_blend_epi32(b, b1, 128);
                }
            }
            b = _mm256_mullo_epi32(b, shift);
            res3 = _mm256_add_epi32(a, b);
            res3 = _mm256_srli_epi32(res3, 5);

            a1 = _mm256_loadu_si256((__m256i *)(leftBy32PlusOffset + base + 56));
            b1 = _mm256_loadu_si256((__m256i *)(leftByDiffPlusOffset + base + 56));
            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 56));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 56));
            for (counter = 24; counter < (-base + 1); counter++)
            {
                if ((counter & 7) == 0)
                {
                    a = _mm256_blend_epi32(a, a1, 1);
                    b = _mm256_blend_epi32(b, b1, 1);
                }
                else if ((counter & 7) == 1)
                {
                    a = _mm256_blend_epi32(a, a1, 2);
                    b = _mm256_blend_epi32(b, b1, 2);
                }
                else if ((counter & 7) == 2)
                {
                    a = _mm256_blend_epi32(a, a1, 4);
                    b = _mm256_blend_epi32(b, b1, 4);
                }
                else if ((counter & 7) == 3)
                {
                    a = _mm256_blend_epi32(a, a1, 8);
                    b = _mm256_blend_epi32(b, b1, 8);
                }
                else if ((counter & 7) == 4)
                {
                    a = _mm256_blend_epi32(a, a1, 16);
                    b = _mm256_blend_epi32(b, b1, 16);
                }
                else if ((counter & 7) == 5)
                {
                    a = _mm256_blend_epi32(a, a1, 32);
                    b = _mm256_blend_epi32(b, b1, 32);
                }
                else if ((counter & 7) == 6)
                {
                    a = _mm256_blend_epi32(a, a1, 64);
                    b = _mm256_blend_epi32(b, b1, 64);
                }
                else
                {
                    a = _mm256_blend_epi32(a, a1, 128);
                    b = _mm256_blend_epi32(b, b1, 128);
                }
            }
            b = _mm256_mullo_epi32(b, shift);
            res4 = _mm256_add_epi32(a, b);
            res4 = _mm256_srli_epi32(res4, 5);

            resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2), PERM4x64(0, 2, 1, 3));
            resHi = _mm256_permute4x64_epi64(_mm256_packus_epi32(res3, res4), PERM4x64(0, 2, 1, 3));
            _mm256_storeu_si256((__m256i *)(dst + c * pitch), resLo);
            _mm256_storeu_si256((__m256i *)(dst + 16 + c * pitch), resHi);
        }
        else
        {
            base = base + 1;
            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 32));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 32));
            b = _mm256_mullo_epi32(b, shift);
            res1 = _mm256_add_epi32(a, b);
            res1 = _mm256_srli_epi32(res1, 5);

            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 40));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 40));
            b = _mm256_mullo_epi32(b, shift);
            res2 = _mm256_add_epi32(a, b);
            res2 = _mm256_srli_epi32(res2, 5);

            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 48));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 48));
            b = _mm256_mullo_epi32(b, shift);
            res3 = _mm256_add_epi32(a, b);
            res3 = _mm256_srli_epi32(res3, 5);

            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 56));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 56));
            b = _mm256_mullo_epi32(b, shift);
            res4 = _mm256_add_epi32(a, b);
            res4 = _mm256_srli_epi32(res4, 5);

            resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2), PERM4x64(0, 2, 1, 3));
            resHi = _mm256_permute4x64_epi64(_mm256_packus_epi32(res3, res4), PERM4x64(0, 2, 1, 3));
            _mm256_storeu_si256((__m256i *)(dst + c * pitch), resLo);
            _mm256_storeu_si256((__m256i *)(dst + 16 + c * pitch), resHi);
        }
    }


}
void Predict_Z2_16bit_64x64(uint16_t *refPel, uint16_t *dst, const uint32_t pitch, uint16_t d)
{
    uint32_t leftBy32[128];
    int_least32_t leftByDiff[128];
    __m256i a0, a1, diff, a32, a16;
    int32_t c, x, y, base;
    const uint16_t* leftPtr = refPel - 1;
    const int32_t frac_bits = 6;

    a16 = _mm256_set1_epi32(16);
    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));      // 01234567
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 1)));    // 89abcdef
    diff = _mm256_sub_epi32(a1, a0);                                        // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0, 5);                                            // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);                                        // a[x] * 32 + 16
    _mm256_storeu_si256((__m256i *)(leftBy32 + 64), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 64), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 8)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 9)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 72), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 72), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 16)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 17)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 80), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 80), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 24)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 25)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 88), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 88), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 32)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 33)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 96), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 96), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 40)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 41)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 104), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 104), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 48)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 49)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 112), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 112), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 56)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 57)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 120), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 120), diff);


    for (c = 0; c < 64; ++c) {
        x = c + 1;
        y = -x * d;
        base = y >> frac_bits;

        __m256i inc = _mm256_set1_epi32(y);
        __m256i shift = _mm256_srli_epi32(_mm256_and_si256(inc, _mm256_set1_epi32(0x3f)), 1);

        __m256i a, b, res1, res2, res3, res4, resLo, resHi;
        base = base + 1;
        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 64));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 64));
        b = _mm256_mullo_epi32(b, shift);
        res1 = _mm256_add_epi32(a, b);
        res1 = _mm256_srli_epi32(res1, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 72));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 72));
        b = _mm256_mullo_epi32(b, shift);
        res2 = _mm256_add_epi32(a, b);
        res2 = _mm256_srli_epi32(res2, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 80));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 80));
        b = _mm256_mullo_epi32(b, shift);
        res3 = _mm256_add_epi32(a, b);
        res3 = _mm256_srli_epi32(res3, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 88));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 88));
        b = _mm256_mullo_epi32(b, shift);
        res4 = _mm256_add_epi32(a, b);
        res4 = _mm256_srli_epi32(res4, 5);

        resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2), PERM4x64(0, 2, 1, 3));
        resHi = _mm256_permute4x64_epi64(_mm256_packus_epi32(res3, res4), PERM4x64(0, 2, 1, 3));
        _mm256_storeu_si256((__m256i *)(dst + c * pitch), resLo);
        _mm256_storeu_si256((__m256i *)(dst + 16 + c * pitch), resHi);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 96));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 96));
        b = _mm256_mullo_epi32(b, shift);
        res1 = _mm256_add_epi32(a, b);
        res1 = _mm256_srli_epi32(res1, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 104));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 104));
        b = _mm256_mullo_epi32(b, shift);
        res2 = _mm256_add_epi32(a, b);
        res2 = _mm256_srli_epi32(res2, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 112));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 112));
        b = _mm256_mullo_epi32(b, shift);
        res3 = _mm256_add_epi32(a, b);
        res3 = _mm256_srli_epi32(res3, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 120));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 120));
        b = _mm256_mullo_epi32(b, shift);
        res4 = _mm256_add_epi32(a, b);
        res4 = _mm256_srli_epi32(res4, 5);

        resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2), PERM4x64(0, 2, 1, 3));
        resHi = _mm256_permute4x64_epi64(_mm256_packus_epi32(res3, res4), PERM4x64(0, 2, 1, 3));
        _mm256_storeu_si256((__m256i *)(dst + 32 + c * pitch), resLo);
        _mm256_storeu_si256((__m256i *)(dst + 48 + c * pitch), resHi);
    }

}
void Predict_Z2_16bit_64x64_Left(uint16_t *refPel, uint16_t *dst, const uint32_t pitch, uint16_t d, uint32_t toAboveLeftOffset)
{
    uint32_t leftBy32[128];
    int_least32_t leftByDiff[128];
    uint32_t leftBy32PlusOffset[128];
    int_least32_t leftByDiffPlusOffset[128];
    __m256i a0, a1, diff, a32, a16, leftPlusOffset;
    int32_t c, x, y, base, counter;
    const uint16_t* leftPtr = refPel - 1;
    const int32_t frac_bits = 6;

    a16 = _mm256_set1_epi32(16);
    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));      // 01234567
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 1)));    // 89abcdef
    diff = _mm256_sub_epi32(a1, a0);                                        // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0, 5);                                            // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);                                        // a[x] * 32 + 16
    _mm256_storeu_si256((__m256i *)(leftBy32 + 64), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 64), diff);

    leftPlusOffset = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 1 + toAboveLeftOffset)));
    diff = _mm256_sub_epi32(a1, leftPlusOffset);
    a32 = _mm256_slli_epi32(leftPlusOffset, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32PlusOffset + 64), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiffPlusOffset + 64), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 8)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 9)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 72), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 72), diff);

    diff = _mm256_sub_epi32(a1, leftPlusOffset);
    a32 = _mm256_slli_epi32(leftPlusOffset, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32PlusOffset + 72), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiffPlusOffset + 72), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 16)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 17)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 80), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 80), diff);

    diff = _mm256_sub_epi32(a1, leftPlusOffset);
    a32 = _mm256_slli_epi32(leftPlusOffset, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32PlusOffset + 80), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiffPlusOffset + 80), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 24)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 25)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 88), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 88), diff);

    diff = _mm256_sub_epi32(a1, leftPlusOffset);
    a32 = _mm256_slli_epi32(leftPlusOffset, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32PlusOffset + 88), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiffPlusOffset + 88), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 32)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 33)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 96), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 96), diff);

    diff = _mm256_sub_epi32(a1, leftPlusOffset);
    a32 = _mm256_slli_epi32(leftPlusOffset, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32PlusOffset + 96), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiffPlusOffset + 96), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 40)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 41)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 104), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 104), diff);

    diff = _mm256_sub_epi32(a1, leftPlusOffset);
    a32 = _mm256_slli_epi32(leftPlusOffset, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32PlusOffset + 104), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiffPlusOffset + 104), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 48)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 49)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 112), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 112), diff);

    diff = _mm256_sub_epi32(a1, leftPlusOffset);
    a32 = _mm256_slli_epi32(leftPlusOffset, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32PlusOffset + 112), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiffPlusOffset + 112), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 56)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 57)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 120), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 120), diff);

    diff = _mm256_sub_epi32(a1, leftPlusOffset);
    a32 = _mm256_slli_epi32(leftPlusOffset, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32PlusOffset + 120), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiffPlusOffset + 120), diff);


    for (c = 0; c < 64; ++c) {
        x = c + 1;
        y = -x * d;
        base = y >> frac_bits;

        __m256i inc = _mm256_set1_epi32(y);
        __m256i shift = _mm256_srli_epi32(_mm256_and_si256(inc, _mm256_set1_epi32(0x3f)), 1);

        __m256i a, b, b1, res1, res2, res3, res4, resLo, resHi;

        if (base <= -1)
        {
            base = base + 1;
            a1 = _mm256_loadu_si256((__m256i *)(leftBy32PlusOffset + base + 64));
            b1 = _mm256_loadu_si256((__m256i *)(leftByDiffPlusOffset + base + 64));
            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 64));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 64));
            for (counter = 0; counter < (-base + 1); counter++)
            {
                if ((counter & 7) == 0)
                {
                    a = _mm256_blend_epi32(a, a1, 1);
                    b = _mm256_blend_epi32(b, b1, 1);
                }
                else if ((counter & 7) == 1)
                {
                    a = _mm256_blend_epi32(a, a1, 2);
                    b = _mm256_blend_epi32(b, b1, 2);
                }
                else if ((counter & 7) == 2)
                {
                    a = _mm256_blend_epi32(a, a1, 4);
                    b = _mm256_blend_epi32(b, b1, 4);
                }
                else if ((counter & 7) == 3)
                {
                    a = _mm256_blend_epi32(a, a1, 8);
                    b = _mm256_blend_epi32(b, b1, 8);
                }
                else if ((counter & 7) == 4)
                {
                    a = _mm256_blend_epi32(a, a1, 16);
                    b = _mm256_blend_epi32(b, b1, 16);
                }
                else if ((counter & 7) == 5)
                {
                    a = _mm256_blend_epi32(a, a1, 32);
                    b = _mm256_blend_epi32(b, b1, 32);
                }
                else if ((counter & 7) == 6)
                {
                    a = _mm256_blend_epi32(a, a1, 64);
                    b = _mm256_blend_epi32(b, b1, 64);
                }
                else
                {
                    a = _mm256_blend_epi32(a, a1, 128);
                    b = _mm256_blend_epi32(b, b1, 128);
                }
            }
            b = _mm256_mullo_epi32(b, shift);
            res1 = _mm256_add_epi32(a, b);
            res1 = _mm256_srli_epi32(res1, 5);

            a1 = _mm256_loadu_si256((__m256i *)(leftBy32PlusOffset + base + 72));
            b1 = _mm256_loadu_si256((__m256i *)(leftByDiffPlusOffset + base + 72));
            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 72));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 72));
            for (counter = 8; counter < (-base + 1); counter++)
            {
                if ((counter & 7) == 0)
                {
                    a = _mm256_blend_epi32(a, a1, 1);
                    b = _mm256_blend_epi32(b, b1, 1);
                }
                else if ((counter & 7) == 1)
                {
                    a = _mm256_blend_epi32(a, a1, 2);
                    b = _mm256_blend_epi32(b, b1, 2);
                }
                else if ((counter & 7) == 2)
                {
                    a = _mm256_blend_epi32(a, a1, 4);
                    b = _mm256_blend_epi32(b, b1, 4);
                }
                else if ((counter & 7) == 3)
                {
                    a = _mm256_blend_epi32(a, a1, 8);
                    b = _mm256_blend_epi32(b, b1, 8);
                }
                else if ((counter & 7) == 4)
                {
                    a = _mm256_blend_epi32(a, a1, 16);
                    b = _mm256_blend_epi32(b, b1, 16);
                }
                else if ((counter & 7) == 5)
                {
                    a = _mm256_blend_epi32(a, a1, 32);
                    b = _mm256_blend_epi32(b, b1, 32);
                }
                else if ((counter & 7) == 6)
                {
                    a = _mm256_blend_epi32(a, a1, 64);
                    b = _mm256_blend_epi32(b, b1, 64);
                }
                else
                {
                    a = _mm256_blend_epi32(a, a1, 128);
                    b = _mm256_blend_epi32(b, b1, 128);
                }
            }
            b = _mm256_mullo_epi32(b, shift);
            res2 = _mm256_add_epi32(a, b);
            res2 = _mm256_srli_epi32(res2, 5);


            a1 = _mm256_loadu_si256((__m256i *)(leftBy32PlusOffset + base + 80));
            b1 = _mm256_loadu_si256((__m256i *)(leftByDiffPlusOffset + base + 80));
            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 80));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 80));
            for (counter = 16; counter < (-base + 1); counter++)
            {
                if ((counter & 7) == 0)
                {
                    a = _mm256_blend_epi32(a, a1, 1);
                    b = _mm256_blend_epi32(b, b1, 1);
                }
                else if ((counter & 7) == 1)
                {
                    a = _mm256_blend_epi32(a, a1, 2);
                    b = _mm256_blend_epi32(b, b1, 2);
                }
                else if ((counter & 7) == 2)
                {
                    a = _mm256_blend_epi32(a, a1, 4);
                    b = _mm256_blend_epi32(b, b1, 4);
                }
                else if ((counter & 7) == 3)
                {
                    a = _mm256_blend_epi32(a, a1, 8);
                    b = _mm256_blend_epi32(b, b1, 8);
                }
                else if ((counter & 7) == 4)
                {
                    a = _mm256_blend_epi32(a, a1, 16);
                    b = _mm256_blend_epi32(b, b1, 16);
                }
                else if ((counter & 7) == 5)
                {
                    a = _mm256_blend_epi32(a, a1, 32);
                    b = _mm256_blend_epi32(b, b1, 32);
                }
                else if ((counter & 7) == 6)
                {
                    a = _mm256_blend_epi32(a, a1, 64);
                    b = _mm256_blend_epi32(b, b1, 64);
                }
                else
                {
                    a = _mm256_blend_epi32(a, a1, 128);
                    b = _mm256_blend_epi32(b, b1, 128);
                }
            }
            b = _mm256_mullo_epi32(b, shift);
            res3 = _mm256_add_epi32(a, b);
            res3 = _mm256_srli_epi32(res3, 5);

            a1 = _mm256_loadu_si256((__m256i *)(leftBy32PlusOffset + base + 88));
            b1 = _mm256_loadu_si256((__m256i *)(leftByDiffPlusOffset + base + 88));
            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 88));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 88));
            for (counter = 24; counter < (-base + 1); counter++)
            {
                if ((counter & 7) == 0)
                {
                    a = _mm256_blend_epi32(a, a1, 1);
                    b = _mm256_blend_epi32(b, b1, 1);
                }
                else if ((counter & 7) == 1)
                {
                    a = _mm256_blend_epi32(a, a1, 2);
                    b = _mm256_blend_epi32(b, b1, 2);
                }
                else if ((counter & 7) == 2)
                {
                    a = _mm256_blend_epi32(a, a1, 4);
                    b = _mm256_blend_epi32(b, b1, 4);
                }
                else if ((counter & 7) == 3)
                {
                    a = _mm256_blend_epi32(a, a1, 8);
                    b = _mm256_blend_epi32(b, b1, 8);
                }
                else if ((counter & 7) == 4)
                {
                    a = _mm256_blend_epi32(a, a1, 16);
                    b = _mm256_blend_epi32(b, b1, 16);
                }
                else if ((counter & 7) == 5)
                {
                    a = _mm256_blend_epi32(a, a1, 32);
                    b = _mm256_blend_epi32(b, b1, 32);
                }
                else if ((counter & 7) == 6)
                {
                    a = _mm256_blend_epi32(a, a1, 64);
                    b = _mm256_blend_epi32(b, b1, 64);
                }
                else
                {
                    a = _mm256_blend_epi32(a, a1, 128);
                    b = _mm256_blend_epi32(b, b1, 128);
                }
            }
            b = _mm256_mullo_epi32(b, shift);
            res4 = _mm256_add_epi32(a, b);
            res4 = _mm256_srli_epi32(res4, 5);

            resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2), PERM4x64(0, 2, 1, 3));
            resHi = _mm256_permute4x64_epi64(_mm256_packus_epi32(res3, res4), PERM4x64(0, 2, 1, 3));
            _mm256_storeu_si256((__m256i *)(dst + c * pitch), resLo);
            _mm256_storeu_si256((__m256i *)(dst + 16 + c * pitch), resHi);


            a1 = _mm256_loadu_si256((__m256i *)(leftBy32PlusOffset + base + 96));
            b1 = _mm256_loadu_si256((__m256i *)(leftByDiffPlusOffset + base + 96));
            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 96));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 96));
            for (counter = 32; counter < (-base + 1); counter++)
            {
                if ((counter & 7) == 0)
                {
                    a = _mm256_blend_epi32(a, a1, 1);
                    b = _mm256_blend_epi32(b, b1, 1);
                }
                else if ((counter & 7) == 1)
                {
                    a = _mm256_blend_epi32(a, a1, 2);
                    b = _mm256_blend_epi32(b, b1, 2);
                }
                else if ((counter & 7) == 2)
                {
                    a = _mm256_blend_epi32(a, a1, 4);
                    b = _mm256_blend_epi32(b, b1, 4);
                }
                else if ((counter & 7) == 3)
                {
                    a = _mm256_blend_epi32(a, a1, 8);
                    b = _mm256_blend_epi32(b, b1, 8);
                }
                else if ((counter & 7) == 4)
                {
                    a = _mm256_blend_epi32(a, a1, 16);
                    b = _mm256_blend_epi32(b, b1, 16);
                }
                else if ((counter & 7) == 5)
                {
                    a = _mm256_blend_epi32(a, a1, 32);
                    b = _mm256_blend_epi32(b, b1, 32);
                }
                else if ((counter & 7) == 6)
                {
                    a = _mm256_blend_epi32(a, a1, 64);
                    b = _mm256_blend_epi32(b, b1, 64);
                }
                else
                {
                    a = _mm256_blend_epi32(a, a1, 128);
                    b = _mm256_blend_epi32(b, b1, 128);
                }
            }
            b = _mm256_mullo_epi32(b, shift);
            res1 = _mm256_add_epi32(a, b);
            res1 = _mm256_srli_epi32(res1, 5);

            a1 = _mm256_loadu_si256((__m256i *)(leftBy32PlusOffset + base + 104));
            b1 = _mm256_loadu_si256((__m256i *)(leftByDiffPlusOffset + base + 104));
            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 104));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 104));
            for (counter = 40; counter < (-base + 1); counter++)
            {
                if ((counter & 7) == 0)
                {
                    a = _mm256_blend_epi32(a, a1, 1);
                    b = _mm256_blend_epi32(b, b1, 1);
                }
                else if ((counter & 7) == 1)
                {
                    a = _mm256_blend_epi32(a, a1, 2);
                    b = _mm256_blend_epi32(b, b1, 2);
                }
                else if ((counter & 7) == 2)
                {
                    a = _mm256_blend_epi32(a, a1, 4);
                    b = _mm256_blend_epi32(b, b1, 4);
                }
                else if ((counter & 7) == 3)
                {
                    a = _mm256_blend_epi32(a, a1, 8);
                    b = _mm256_blend_epi32(b, b1, 8);
                }
                else if ((counter & 7) == 4)
                {
                    a = _mm256_blend_epi32(a, a1, 16);
                    b = _mm256_blend_epi32(b, b1, 16);
                }
                else if ((counter & 7) == 5)
                {
                    a = _mm256_blend_epi32(a, a1, 32);
                    b = _mm256_blend_epi32(b, b1, 32);
                }
                else if ((counter & 7) == 6)
                {
                    a = _mm256_blend_epi32(a, a1, 64);
                    b = _mm256_blend_epi32(b, b1, 64);
                }
                else
                {
                    a = _mm256_blend_epi32(a, a1, 128);
                    b = _mm256_blend_epi32(b, b1, 128);
                }
            }
            b = _mm256_mullo_epi32(b, shift);
            res2 = _mm256_add_epi32(a, b);
            res2 = _mm256_srli_epi32(res2, 5);

            a1 = _mm256_loadu_si256((__m256i *)(leftBy32PlusOffset + base + 112));
            b1 = _mm256_loadu_si256((__m256i *)(leftByDiffPlusOffset + base + 112));
            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 112));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 112));
            for (counter = 48; counter < (-base + 1); counter++)
            {
                if ((counter & 7) == 0)
                {
                    a = _mm256_blend_epi32(a, a1, 1);
                    b = _mm256_blend_epi32(b, b1, 1);
                }
                else if ((counter & 7) == 1)
                {
                    a = _mm256_blend_epi32(a, a1, 2);
                    b = _mm256_blend_epi32(b, b1, 2);
                }
                else if ((counter & 7) == 2)
                {
                    a = _mm256_blend_epi32(a, a1, 4);
                    b = _mm256_blend_epi32(b, b1, 4);
                }
                else if ((counter & 7) == 3)
                {
                    a = _mm256_blend_epi32(a, a1, 8);
                    b = _mm256_blend_epi32(b, b1, 8);
                }
                else if ((counter & 7) == 4)
                {
                    a = _mm256_blend_epi32(a, a1, 16);
                    b = _mm256_blend_epi32(b, b1, 16);
                }
                else if ((counter & 7) == 5)
                {
                    a = _mm256_blend_epi32(a, a1, 32);
                    b = _mm256_blend_epi32(b, b1, 32);
                }
                else if ((counter & 7) == 6)
                {
                    a = _mm256_blend_epi32(a, a1, 64);
                    b = _mm256_blend_epi32(b, b1, 64);
                }
                else
                {
                    a = _mm256_blend_epi32(a, a1, 128);
                    b = _mm256_blend_epi32(b, b1, 128);
                }
            }
            b = _mm256_mullo_epi32(b, shift);
            res3 = _mm256_add_epi32(a, b);
            res3 = _mm256_srli_epi32(res3, 5);

            a1 = _mm256_loadu_si256((__m256i *)(leftBy32PlusOffset + base + 120));
            b1 = _mm256_loadu_si256((__m256i *)(leftByDiffPlusOffset + base + 120));
            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 120));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 120));
            for (counter = 56; counter < (-base + 1); counter++)
            {
                if ((counter & 7) == 0)
                {
                    a = _mm256_blend_epi32(a, a1, 1);
                    b = _mm256_blend_epi32(b, b1, 1);
                }
                else if ((counter & 7) == 1)
                {
                    a = _mm256_blend_epi32(a, a1, 2);
                    b = _mm256_blend_epi32(b, b1, 2);
                }
                else if ((counter & 7) == 2)
                {
                    a = _mm256_blend_epi32(a, a1, 4);
                    b = _mm256_blend_epi32(b, b1, 4);
                }
                else if ((counter & 7) == 3)
                {
                    a = _mm256_blend_epi32(a, a1, 8);
                    b = _mm256_blend_epi32(b, b1, 8);
                }
                else if ((counter & 7) == 4)
                {
                    a = _mm256_blend_epi32(a, a1, 16);
                    b = _mm256_blend_epi32(b, b1, 16);
                }
                else if ((counter & 7) == 5)
                {
                    a = _mm256_blend_epi32(a, a1, 32);
                    b = _mm256_blend_epi32(b, b1, 32);
                }
                else if ((counter & 7) == 6)
                {
                    a = _mm256_blend_epi32(a, a1, 64);
                    b = _mm256_blend_epi32(b, b1, 64);
                }
                else
                {
                    a = _mm256_blend_epi32(a, a1, 128);
                    b = _mm256_blend_epi32(b, b1, 128);
                }
            }
            b = _mm256_mullo_epi32(b, shift);
            res4 = _mm256_add_epi32(a, b);
            res4 = _mm256_srli_epi32(res4, 5);

            resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2), PERM4x64(0, 2, 1, 3));
            resHi = _mm256_permute4x64_epi64(_mm256_packus_epi32(res3, res4), PERM4x64(0, 2, 1, 3));
            _mm256_storeu_si256((__m256i *)(dst + 32 + c * pitch), resLo);
            _mm256_storeu_si256((__m256i *)(dst + 48 + c * pitch), resHi);
        }
        else
        {
            base = base + 1;
            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 64));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 64));
            b = _mm256_mullo_epi32(b, shift);
            res1 = _mm256_add_epi32(a, b);
            res1 = _mm256_srli_epi32(res1, 5);

            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 72));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 72));
            b = _mm256_mullo_epi32(b, shift);
            res2 = _mm256_add_epi32(a, b);
            res2 = _mm256_srli_epi32(res2, 5);

            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 80));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 80));
            b = _mm256_mullo_epi32(b, shift);
            res3 = _mm256_add_epi32(a, b);
            res3 = _mm256_srli_epi32(res3, 5);

            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 88));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 88));
            b = _mm256_mullo_epi32(b, shift);
            res4 = _mm256_add_epi32(a, b);
            res4 = _mm256_srli_epi32(res4, 5);

            resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2), PERM4x64(0, 2, 1, 3));
            resHi = _mm256_permute4x64_epi64(_mm256_packus_epi32(res3, res4), PERM4x64(0, 2, 1, 3));
            _mm256_storeu_si256((__m256i *)(dst + c * pitch), resLo);
            _mm256_storeu_si256((__m256i *)(dst + 16 + c * pitch), resHi);

            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 96));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 96));
            b = _mm256_mullo_epi32(b, shift);
            res1 = _mm256_add_epi32(a, b);
            res1 = _mm256_srli_epi32(res1, 5);

            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 104));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 104));
            b = _mm256_mullo_epi32(b, shift);
            res2 = _mm256_add_epi32(a, b);
            res2 = _mm256_srli_epi32(res2, 5);

            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 112));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 112));
            b = _mm256_mullo_epi32(b, shift);
            res3 = _mm256_add_epi32(a, b);
            res3 = _mm256_srli_epi32(res3, 5);

            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 120));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 120));
            b = _mm256_mullo_epi32(b, shift);
            res4 = _mm256_add_epi32(a, b);
            res4 = _mm256_srli_epi32(res4, 5);

            resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2), PERM4x64(0, 2, 1, 3));
            resHi = _mm256_permute4x64_epi64(_mm256_packus_epi32(res3, res4), PERM4x64(0, 2, 1, 3));
            _mm256_storeu_si256((__m256i *)(dst + 32 + c * pitch), resLo);
            _mm256_storeu_si256((__m256i *)(dst + 48 + c * pitch), resHi);
        }
    }


}

void intra_mode_angular_av1_z2_16bit_4x4_avx2(const uint32_t size, uint16_t *ref_samples, uint16_t *dst, const uint32_t prediction_buffer_stride, const EbBool skip, uint16_t dx, uint16_t dy, uint16_t bd)
{
    uint32_t toAboveOffset = (size << 1) + 1;
    uint32_t toAboveLeftOffset = (size << 1);
    uint32_t toLeftOffset = 0;
    uint32_t rowStride = skip ? 2 : 1;
    int32_t r, x, base, base1, c_pos;
    __m128i ma, ml, a, l, aboveReg, leftReg, result;
    const int32_t frac_bits_x = 6;

    uint16_t dstTmp[4 * 4] = { 0 };
    uint16_t dstTmpTransp[4 * 4] = { 0 };
    int32_t strideTmp = 4;
    __m128i clip_bd = (bd == 8) ? _mm_set1_epi16(255) : (bd == 10) ? _mm_set1_epi16(1023) : _mm_set1_epi16(4095);

    Predict_Z2_16bit_4x4(&ref_samples[toAboveOffset], dst, rowStride* prediction_buffer_stride, dx);
    Predict_Z2_16bit_4x4_Left(&ref_samples[toLeftOffset], dstTmp, strideTmp, dy, toAboveLeftOffset);
    transpose_16bit_TX_4X4(dstTmp, strideTmp, dstTmpTransp, strideTmp);

    for (r = 0; r < 4; r++) {
        x = -dx * (r + 1);
        base1 = x >> frac_bits_x;
        c_pos = -1 - base1;
        base = 0;
        if ((c_pos >= 0) && (c_pos < 64)) {
            base = 64 - c_pos;
        }

        ma = _mm_loadu_si128((__m128i*)(z2BlendMaskabove_16bit + base));
        ml = _mm_loadu_si128((__m128i*)(z2BlendMaskleft_16bit + base));
        a = _mm_loadu_si128((__m128i*)(dst + r * rowStride* prediction_buffer_stride));
        l = _mm_loadu_si128((__m128i*)(dstTmpTransp + r * strideTmp));
        aboveReg = _mm_and_si128(a, ma);
        leftReg = _mm_and_si128(l, ml);
        result = _mm_or_si128(aboveReg, leftReg);
        result = _mm_min_epi16(clip_bd, result);
        _mm_storel_epi64((__m128i*)(dst + r * rowStride* prediction_buffer_stride), result);

    }
}
void intra_mode_angular_av1_z2_16bit_8x8_avx2(const uint32_t size, uint16_t *ref_samples, uint16_t *dst, const uint32_t prediction_buffer_stride, const EbBool skip, uint16_t dx, uint16_t dy, uint16_t bd)
{
    uint32_t toAboveOffset = (size << 1) + 1;
    uint32_t toAboveLeftOffset = (size << 1);
    uint32_t toLeftOffset = 0;
    uint32_t rowStride = skip ? 2 : 1;
    int32_t r, x, base, base1, c_pos;
    const int32_t frac_bits_x = 6;
    __m128i ma, ml, a, l, aboveReg, leftReg, result;
    uint16_t dstTmp[8 * 8] = { 0 };
    uint16_t dstTmpTransp[8 * 8] = { 0 };
    int32_t strideTmp = 8;
    __m128i clip_bd = (bd == 8) ? _mm_set1_epi16(255) : (bd == 10) ? _mm_set1_epi16(1023) : _mm_set1_epi16(4095);

    Predict_Z2_16bit_8x8(&ref_samples[toAboveOffset], dst, rowStride* prediction_buffer_stride, dx);
    Predict_Z2_16bit_8x8_Left(&ref_samples[toLeftOffset], dstTmp, strideTmp, dy, toAboveLeftOffset);
    transpose_16bit(dstTmp, strideTmp, dstTmpTransp, strideTmp, 8, 8);

    for (r = 0; r < 8; r++) {
        x = -dx * (r + 1);
        base1 = x >> frac_bits_x;
        c_pos = -1 - base1;
        base = 0;
        if ((c_pos >= 0) && (c_pos < 64)) {
            base = 64 - c_pos;
        }

        ma = _mm_loadu_si128((__m128i*)(z2BlendMaskabove_16bit + base));
        ml = _mm_loadu_si128((__m128i*)(z2BlendMaskleft_16bit + base));
        a = _mm_loadu_si128((__m128i*)(dst + r * rowStride* prediction_buffer_stride));
        l = _mm_loadu_si128((__m128i*)(dstTmpTransp + r * strideTmp));
        aboveReg = _mm_and_si128(a, ma);
        leftReg = _mm_and_si128(l, ml);
        result = _mm_or_si128(aboveReg, leftReg);
        result = _mm_min_epi16(clip_bd, result);
        _mm_storeu_si128((__m128i*)(dst + r * rowStride* prediction_buffer_stride), result);
    }
}
void intra_mode_angular_av1_z2_16bit_16x16_avx2(const uint32_t size, uint16_t *ref_samples, uint16_t *dst, const uint32_t prediction_buffer_stride, const EbBool skip, uint16_t dx, uint16_t dy, uint16_t bd)
{
    uint32_t toAboveOffset = (size << 1) + 1;
    uint32_t toAboveLeftOffset = (size << 1);
    uint32_t toLeftOffset = 0;
    uint32_t rowStride = skip ? 2 : 1;
    int32_t r, x, base, base1, c_pos;
    const int32_t frac_bits_x = 6;
    __m256i ma, ml, a, l, aboveReg, leftReg, result;
    uint16_t dstTmp[16 * 16] = { 0 };
    uint16_t dstTmpTransp[16 * 16] = { 0 };
    int32_t strideTmp = 16;
    __m256i clip_bd = (bd == 8) ? _mm256_set1_epi16(255) : (bd == 10) ? _mm256_set1_epi16(1023) : _mm256_set1_epi16(4095);

    Predict_Z2_16bit_16x16(&ref_samples[toAboveOffset], dst, rowStride* prediction_buffer_stride, dx);
    Predict_Z2_16bit_16x16_Left(&ref_samples[toLeftOffset], dstTmp, strideTmp, dy, toAboveLeftOffset);
    transpose_16bit(dstTmp, strideTmp, dstTmpTransp, strideTmp, 16, 16);

    for (r = 0; r < 16; r++) {
        x = -dx * (r + 1);
        base1 = x >> frac_bits_x;
        c_pos = -1 - base1;
        base = 0;
        if (c_pos >= 0 && c_pos < 64) {
            base = 64 - c_pos;
        }

        ma = _mm256_loadu_si256((__m256i*)(z2BlendMaskabove_16bit + base));
        ml = _mm256_loadu_si256((__m256i*)(z2BlendMaskleft_16bit + base));
        a = _mm256_loadu_si256((__m256i*)(dst + r * rowStride* prediction_buffer_stride));
        l = _mm256_loadu_si256((__m256i*)(dstTmpTransp + r * strideTmp));
        aboveReg = _mm256_and_si256(a, ma);
        leftReg = _mm256_and_si256(l, ml);
        result = _mm256_or_si256(aboveReg, leftReg);
        result = _mm256_min_epi16(clip_bd, result);
        _mm256_storeu_si256((__m256i*)(dst + r * rowStride* prediction_buffer_stride), result);
    }
}
void intra_mode_angular_av1_z2_16bit_32x32_avx2(const uint32_t size, uint16_t *ref_samples, uint16_t *dst, const uint32_t prediction_buffer_stride, const EbBool skip, uint16_t dx, uint16_t dy, uint16_t bd)
{
    uint32_t toAboveOffset = (size << 1) + 1;
    uint32_t toLeftOffset = 0;
    uint32_t rowStride = skip ? 2 : 1;
    int32_t r, c, x, base, base1, c_pos;
    const int32_t frac_bits_x = 6;
    __m256i ma, ml, a, l, aboveReg, leftReg, result;
    uint16_t dstTmp[32 * 32] = { 0 };
    uint16_t dstTmpTransp[32 * 32] = { 0 };
    int32_t strideTmp = 32;
    __m256i clip_bd = (bd == 8) ? _mm256_set1_epi16(255) : (bd == 10) ? _mm256_set1_epi16(1023) : _mm256_set1_epi16(4095);

    Predict_Z2_16bit_32x32(&ref_samples[toAboveOffset], dst, rowStride* prediction_buffer_stride, dx);
    Predict_Z2_16bit_32x32_Left(&ref_samples[toLeftOffset], dstTmp, strideTmp, dy, (size << 1));
    transpose_16bit(dstTmp, strideTmp, dstTmpTransp, strideTmp, 32, 32);

    for (r = 0; r < 32; r++) {
        x = -dx * (r + 1);
        base1 = x >> frac_bits_x;
        c_pos = -1 - base1;
        base = 0;
        if (c_pos >= 0 && c_pos < 64) {
            base = 64 - c_pos;
        }

        for (c = 0; c < 32; c += 16) {
            ma = _mm256_loadu_si256((__m256i*)(z2BlendMaskabove_16bit + base + c));
            ml = _mm256_loadu_si256((__m256i*)(z2BlendMaskleft_16bit + base + c));
            a = _mm256_loadu_si256((__m256i*)(dst + r * rowStride* prediction_buffer_stride + c));
            l = _mm256_loadu_si256((__m256i*)(dstTmpTransp + r * strideTmp + c));
            aboveReg = _mm256_and_si256(a, ma);
            leftReg = _mm256_and_si256(l, ml);
            result = _mm256_or_si256(aboveReg, leftReg);
            result = _mm256_min_epi16(clip_bd, result);
            _mm256_storeu_si256((__m256i*)(dst + r * rowStride* prediction_buffer_stride + c), result);
        }
    }
}
void intra_mode_angular_av1_z2_16bit_64x64_avx2(const uint32_t size, uint16_t *ref_samples, uint16_t *dst, const uint32_t prediction_buffer_stride, const EbBool skip, uint16_t dx, uint16_t dy, uint16_t bd)
{
    uint32_t toAboveOffset = (size << 1) + 1;
    uint32_t toAboveLeftOffset = (size << 1);
    uint32_t toLeftOffset = 0;
    uint32_t rowStride = skip ? 2 : 1;
    int32_t r, c, x, base, base1, c_pos;
    const int32_t frac_bits_x = 6;
    __m256i ma, ml, a, l, aboveReg, leftReg, result;
    uint16_t dstTmp[64 * 64] = { 0 };
    uint16_t dstTmpTransp[64 * 64] = { 0 };
    int32_t strideTmp = 64;
    __m256i clip_bd = (bd == 8) ? _mm256_set1_epi16(255) : (bd == 10) ? _mm256_set1_epi16(1023) : _mm256_set1_epi16(4095);

    Predict_Z2_16bit_64x64(&ref_samples[toAboveOffset], dst, rowStride* prediction_buffer_stride, dx);
    Predict_Z2_16bit_64x64_Left(&ref_samples[toLeftOffset], dstTmp, strideTmp, dy, toAboveLeftOffset);
    transpose_16bit(dstTmp, strideTmp, dstTmpTransp, strideTmp, 64, 64);

    for (r = 0; r < 64; r++) {
        x = -dx * (r + 1);
        base1 = x >> frac_bits_x;
        c_pos = -1 - base1;
        base = 0;
        if (c_pos >= 0 && c_pos < 64) {
            base = 64 - c_pos;
        }

        for (c = 0; c < 64; c += 16) {
            ma = _mm256_loadu_si256((__m256i*)(z2BlendMaskabove_16bit + base + c));
            ml = _mm256_loadu_si256((__m256i*)(z2BlendMaskleft_16bit + base + c));
            a = _mm256_loadu_si256((__m256i*)(dst + r * rowStride* prediction_buffer_stride + c));
            l = _mm256_loadu_si256((__m256i*)(dstTmpTransp + r * strideTmp + c));
            aboveReg = _mm256_and_si256(a, ma);
            leftReg = _mm256_and_si256(l, ml);
            result = _mm256_or_si256(aboveReg, leftReg);
            result = _mm256_min_epi16(clip_bd, result);
            _mm256_storeu_si256((__m256i*)(dst + r * rowStride* prediction_buffer_stride + c), result);
        }
    }
}

void intra_mode_angular_av1_z3_16bit_4x4_avx2(const uint32_t size, uint16_t *ref_samples, uint16_t *dst, const uint32_t prediction_buffer_stride, const EbBool skip, uint16_t dx, uint16_t dy, uint16_t bd)
{
    (void)size;
    (void)dx;
    uint16_t dstT[4 * 4];
    uint32_t toLeftOffset = 0;
    uint32_t rowStride = skip ? 2 : 1;
    Predict_Z1_16bit_4x4(&ref_samples[toLeftOffset], dstT, 4, dy, bd);
    transpose_16bit_TX_4X4(dstT, 4, dst, rowStride*prediction_buffer_stride);
}
void intra_mode_angular_av1_z3_16bit_8x8_avx2(const uint32_t size, uint16_t *ref_samples, uint16_t *dst, const uint32_t prediction_buffer_stride, const EbBool skip, uint16_t dx, uint16_t dy, uint16_t bd)
{
    (void)size;
    (void)dx;
    uint16_t dstT[8 * 8];
    uint32_t toLeftOffset = 0;
    uint32_t rowStride = skip ? 2 : 1;
    Predict_Z1_16bit_8x8(&ref_samples[toLeftOffset], dstT, 8, dy, bd);
    transpose_16bit(dstT, 8, dst, rowStride*prediction_buffer_stride, 8, 8);
}
void intra_mode_angular_av1_z3_16bit_16x16_avx2(const uint32_t size, uint16_t *ref_samples, uint16_t *dst, const uint32_t prediction_buffer_stride, const EbBool skip, uint16_t dx, uint16_t dy, uint16_t bd)
{
    (void)size;
    (void)dx;
    uint16_t dstT[16 * 16];
    uint32_t toLeftOffset = 0;
    uint32_t rowStride = skip ? 2 : 1;
    Predict_Z1_16bit_16x16(&ref_samples[toLeftOffset], dstT, 16, dy, bd);
    transpose_16bit(dstT, 16, dst, rowStride*prediction_buffer_stride, 16, 16);
}
void intra_mode_angular_av1_z3_16bit_32x32_avx2(const uint32_t size, uint16_t *ref_samples, uint16_t *dst, const uint32_t prediction_buffer_stride, const EbBool skip, uint16_t dx, uint16_t dy, uint16_t bd)
{
    (void)size;
    (void)dx;
    uint16_t dstT[32 * 32];
    uint32_t toLeftOffset = 0;
    uint32_t rowStride = skip ? 2 : 1;
    Predict_Z1_16bit_32x32(&ref_samples[toLeftOffset], dstT, 32, dy, bd);
    transpose_16bit(dstT, 32, dst, rowStride*prediction_buffer_stride, 32, 32);
}
void intra_mode_angular_av1_z3_16bit_64x64_avx2(const uint32_t size, uint16_t *ref_samples, uint16_t *dst, const uint32_t prediction_buffer_stride, const EbBool skip, uint16_t dx, uint16_t dy, uint16_t bd)
{
    (void)size;
    (void)dx;
    uint16_t dstT[64 * 64];
    uint32_t toLeftOffset = 0;
    uint32_t rowStride = skip ? 2 : 1;
    Predict_Z1_16bit_64x64(&ref_samples[toLeftOffset], dstT, 64, dy, bd);
    transpose_16bit(dstT, 64, dst, rowStride*prediction_buffer_stride, 64, 64);
}
#if INTRAD_ASM
//#define PERM4x64(c0, c1, c2, c3) c0 + (c1 << 2) + (c2 << 4) + (c3 << 6)
//#define PERM2x128(c0, c1) c0 + (c1 << 4)


#define _mm256_set_m128i(/* __m128i */ hi, /* __m128i */ lo) \
    _mm256_insertf128_si256(_mm256_castsi128_si256(lo), (hi), 0x1)

// Low bit depth functions
static uint8_t BaseMask[33][32] = {
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0 },
    { 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff },
};

static AOM_FORCE_INLINE void dr_prediction_z1_4xN_internal_avx2(
    int32_t N, __m128i *dst, const uint8_t *above, int32_t upsample_above, int32_t dx) {
    const int32_t frac_bits = 6 - upsample_above;
    const int32_t max_base_x = ((N + 4) - 1) << upsample_above;
    int32_t x;
    // a assert(dx > 0);
    // pre-filter above pixels
    // store in temp buffers:
    //   above[x] * 32 + 16
    //   above[x+1] - above[x]
    // final pixels will be caluculated as:
    //   (above[x] * 32 + 16 + (above[x+1] - above[x]) * shift) >> 5
    __m256i a0, a1, a32, a16;
    __m256i diff;
    __m128i a_mbase_x;

    a16 = _mm256_set1_epi32(16);
    a_mbase_x = _mm_set1_epi8(above[max_base_x]);

    x = dx;
    for (int32_t r = 0; r < N; r++) {
        __m256i b, res, shift;
        __m128i res1, a0_128, a1_128;

        int32_t base = x >> frac_bits;
        int32_t base_max_diff = (max_base_x - base) >> upsample_above;
        if (base_max_diff <= 0) {
            for (int32_t i = r; i < N; ++i) {
                dst[i] = a_mbase_x;  // save 4 values
            }
            return;
        }
        if (base_max_diff > 4) base_max_diff = 4;
        a0_128 = _mm_loadu_si128((__m128i *)(above + base));
        a1_128 = _mm_srli_si128(a0_128, 1);
        a0 = _mm256_cvtepu8_epi32(a0_128);
        a1 = _mm256_cvtepu8_epi32(a1_128);

        if (upsample_above) {
            a0 = _mm256_permutevar8x32_epi32(
                a0, _mm256_set_epi32(7, 5, 3, 1, 6, 4, 2, 0));
            a1 = _mm256_castsi128_si256(_mm256_extracti128_si256(a0, 1));
            shift = _mm256_srli_epi32(
                _mm256_and_si256(
                _mm256_slli_epi32(_mm256_set1_epi32(x), upsample_above),
                _mm256_set1_epi32(0x3f)),
                1);
        }
        else {
            shift = _mm256_srli_epi32(
                _mm256_and_si256(_mm256_set1_epi32(x), _mm256_set1_epi32(0x3f)), 1);
        }

        diff = _mm256_sub_epi32(a1, a0);   // a[x+1] - a[x]
        a32 = _mm256_slli_epi32(a0, 5);    // a[x] * 32
        a32 = _mm256_add_epi32(a32, a16);  // a[x] * 32 + 16

        b = _mm256_mullo_epi32(diff, shift);
        res = _mm256_add_epi32(a32, b);
        res = _mm256_srli_epi32(res, 5);

        res1 = _mm256_castsi256_si128(res);
        res1 = _mm_packus_epi32(res1, res1);
        res1 = _mm_packus_epi16(res1, res1);

        dst[r] =
            _mm_blendv_epi8(a_mbase_x, res1, *(__m128i *)BaseMask[base_max_diff]);
        x += dx;
    }
}
static void dr_prediction_z1_4xN_avx2(int32_t N, uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above, int32_t upsample_above,
    int32_t dx) {
    __m128i dstvec[16];

    dr_prediction_z1_4xN_internal_avx2(N, dstvec, above, upsample_above, dx);
    for (int32_t i = 0; i < N; i++) {
        *(uint32_t *)(dst + stride * i) = _mm_cvtsi128_si32(dstvec[i]);
    }
}

static AOM_FORCE_INLINE void dr_prediction_z1_8xN_internal_avx2(
    int32_t N, __m128i *dst, const uint8_t *above, int32_t upsample_above, int32_t dx) {
    const int32_t frac_bits = 6 - upsample_above;
    const int32_t max_base_x = ((8 + N) - 1) << upsample_above;

    int32_t x;
    // pre-filter above pixels
    // store in temp buffers:
    //   above[x] * 32 + 16
    //   above[x+1] - above[x]
    // final pixels will be caluculated as:
    //   (above[x] * 32 + 16 + (above[x+1] - above[x]) * shift) >> 5
    __m256i a0, a1, a0_1, a1_1, a32, a16, diff;
    __m128i a_mbase_x;

    a16 = _mm256_set1_epi32(16);
    a_mbase_x = _mm_set1_epi8(above[max_base_x]);

    x = dx;
    for (int32_t r = 0; r < N; r++) {
        __m256i b, res, res1, shift;
        __m128i res128;

        int32_t base = x >> frac_bits;
        int32_t base_max_diff = (max_base_x - base) >> upsample_above;
        if (base_max_diff <= 0) {
            for (int32_t i = r; i < N; ++i) {
                dst[i] = a_mbase_x;  // save 16 values, 8 to be used furter
            }
            return;
        }
        if (base_max_diff > 8) base_max_diff = 8;

        a0 = _mm256_cvtepu8_epi32(_mm_loadu_si128((__m128i *)(above + base)));
        a1 = _mm256_cvtepu8_epi32(_mm_loadu_si128((__m128i *)(above + base + 1)));

        if (upsample_above) {
            a0 = _mm256_permutevar8x32_epi32(
                a0, _mm256_set_epi32(7, 5, 3, 1, 6, 4, 2, 0));
            a1 = _mm256_castsi128_si256(_mm256_extracti128_si256(a0, 1));

            a0_1 =
                _mm256_cvtepu8_epi32(_mm_loadu_si128((__m128i *)(above + base + 8)));
            a0_1 = _mm256_permutevar8x32_epi32(
                a0_1, _mm256_set_epi32(7, 5, 3, 1, 6, 4, 2, 0));
            a1_1 = _mm256_castsi128_si256(_mm256_extracti128_si256(a0_1, 1));

            a0 = _mm256_inserti128_si256(a0, _mm256_castsi256_si128(a0_1), 1);
            a1 = _mm256_inserti128_si256(a1, _mm256_castsi256_si128(a1_1), 1);

            shift = _mm256_srli_epi32(
                _mm256_and_si256(
                _mm256_slli_epi32(_mm256_set1_epi32(x), upsample_above),
                _mm256_set1_epi32(0x3f)),
                1);
        }
        else {
            shift = _mm256_srli_epi32(
                _mm256_and_si256(_mm256_set1_epi32(x), _mm256_set1_epi32(0x3f)), 1);
        }

        diff = _mm256_sub_epi32(a1, a0);   // a[x+1] - a[x]
        a32 = _mm256_slli_epi32(a0, 5);    // a[x] * 32
        a32 = _mm256_add_epi32(a32, a16);  // a[x] * 32 + 16

        b = _mm256_mullo_epi32(diff, shift);
        res = _mm256_add_epi32(a32, b);
        res = _mm256_srli_epi32(res, 5);

        res1 = _mm256_packus_epi32(
            res, _mm256_castsi128_si256(
            _mm256_extracti128_si256(res, 1)));  // goto 16 bit

        res128 = _mm_packus_epi16(_mm256_castsi256_si128(res1),
            _mm256_castsi256_si128(res1));  // goto 8 bit

        res128 =
            _mm_blendv_epi8(a_mbase_x, res128, *(__m128i *)BaseMask[base_max_diff]);
        dst[r] = res128;
        x += dx;
    }
}

static void dr_prediction_z1_8xN_avx2(int32_t N, uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above, int32_t upsample_above,
    int32_t dx) {
    __m128i dstvec[32];

    dr_prediction_z1_8xN_internal_avx2(N, dstvec, above, upsample_above, dx);
    for (int32_t i = 0; i < N; i++) {
        _mm_storel_epi64((__m128i *)(dst + stride * i), dstvec[i]);
    }
}

static AOM_FORCE_INLINE void dr_prediction_z1_16xN_internal_avx2(
    int32_t N, __m128i *dstvec, const uint8_t *above, int32_t upsample_above, int32_t dx) {
    int32_t x;
    // here upsample_above is 0 by design of av1_use_intra_edge_upsample
    (void)upsample_above;
    const int32_t frac_bits = 6;
    const int32_t max_base_x = ((16 + N) - 1);

    // pre-filter above pixels
    // store in temp buffers:
    //   above[x] * 32 + 16
    //   above[x+1] - above[x]
    // final pixels will be caluculated as:
    //   (above[x] * 32 + 16 + (above[x+1] - above[x]) * shift) >> 5
    __m256i a0, a0_1, a1, a1_1, diff, a32, a16;
    __m128i a_mbase_x;

    a16 = _mm256_set1_epi32(16);
    a_mbase_x = _mm_set1_epi8((uint8_t)above[max_base_x]);

    x = dx;
    for (int32_t r = 0; r < N; r++) {
        __m256i b, res[2];
        __m128i res128[2];
        int32_t base = x >> frac_bits;
        int32_t base_max_diff = (max_base_x - base);
        if (base_max_diff <= 0) {
            for (int32_t i = r; i < N; ++i) {
                dstvec[i] = a_mbase_x;  // save 16 values
            }
            return;
        }
        __m256i shift = _mm256_srli_epi32(
            _mm256_and_si256(_mm256_set1_epi32(x), _mm256_set1_epi32(0x3f)), 1);

        a0 = _mm256_cvtepu8_epi32(_mm_loadu_si128((__m128i *)(above + base)));
        a1 = _mm256_cvtepu8_epi32(_mm_loadu_si128((__m128i *)(above + base + 1)));

        diff = _mm256_sub_epi32(a1, a0);   // a[x+1] - a[x]
        a32 = _mm256_slli_epi32(a0, 5);    // a[x] * 32
        a32 = _mm256_add_epi32(a32, a16);  // a[x] * 32 + 16
        b = _mm256_mullo_epi32(diff, shift);

        res[0] = _mm256_add_epi32(a32, b);
        res[0] = _mm256_srli_epi32(res[0], 5);
        res[0] = _mm256_packus_epi32(
            res[0], _mm256_castsi128_si256(_mm256_extracti128_si256(res[0], 1)));
        res128[0] = _mm_packus_epi16(_mm256_castsi256_si128(res[0]),
            _mm256_castsi256_si128(res[0]));  // goto 8 bit

        if (base_max_diff > 8) {
            if (base_max_diff > 16) base_max_diff = 16;
            a0_1 =
                _mm256_cvtepu8_epi32(_mm_loadu_si128((__m128i *)(above + base + 8)));
            a1_1 =
                _mm256_cvtepu8_epi32(_mm_loadu_si128((__m128i *)(above + base + 9)));

            diff = _mm256_sub_epi32(a1_1, a0_1);  // a[x+1] - a[x]
            a32 = _mm256_slli_epi32(a0_1, 5);     // a[x] * 32
            a32 = _mm256_add_epi32(a32, a16);     // a[x] * 32 + 16
            b = _mm256_mullo_epi32(diff, shift);

            res[1] = _mm256_add_epi32(a32, b);
            res[1] = _mm256_srli_epi32(res[1], 5);
            res[1] = _mm256_packus_epi32(
                res[1], _mm256_castsi128_si256(_mm256_extracti128_si256(res[1], 1)));
            res128[1] =
                _mm_packus_epi16(_mm256_castsi256_si128(res[1]),
                _mm256_castsi256_si128(res[1]));  // goto 8 bit

        }
        else {
            res128[1] = a_mbase_x;
        }
        res128[0] = _mm_unpacklo_epi64(res128[0], res128[1]);  // 16 8bit values

        dstvec[r] = _mm_blendv_epi8(a_mbase_x, res128[0],
            *(__m128i *)BaseMask[base_max_diff]);
        x += dx;
    }
}
static void dr_prediction_z1_16xN_avx2(int32_t N, uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above, int32_t upsample_above,
    int32_t dx) {
    __m128i dstvec[64];

    dr_prediction_z1_16xN_internal_avx2(N, dstvec, above, upsample_above, dx);
    for (int32_t i = 0; i < N; i++) {
        _mm_storeu_si128((__m128i *)(dst + stride * i), dstvec[i]);
    }
}

static AOM_FORCE_INLINE void dr_prediction_z1_32xN_internal_avx2(
    int32_t N, __m256i *dstvec, const uint8_t *above, int32_t upsample_above, int32_t dx) {
    int32_t x;
    // here upsample_above is 0 by design of av1_use_intra_edge_upsample
    (void)upsample_above;
    const int32_t frac_bits = 6;
    const int32_t max_base_x = ((32 + N) - 1);

    // pre-filter above pixels
    // store in temp buffers:
    //   above[x] * 32 + 16
    //   above[x+1] - above[x]
    // final pixels will be caluculated as:
    //   (above[x] * 32 + 16 + (above[x+1] - above[x]) * shift) >> 5
    __m256i a0, a0_1, a1, a1_1, a32, a16;
    __m256i a_mbase_x, diff;

    a16 = _mm256_set1_epi32(16);
    a_mbase_x = _mm256_set1_epi8(above[max_base_x]);

    x = dx;
    for (int32_t r = 0; r < N; r++) {
        __m256i b, res[2], res16[2];

        int32_t base = x >> frac_bits;
        int32_t base_max_diff = (max_base_x - base);
        if (base_max_diff <= 0) {
            for (int32_t i = r; i < N; ++i) {
                dstvec[i] = a_mbase_x;  // save 32 values
            }
            return;
        }
        if (base_max_diff > 32) base_max_diff = 32;
        __m256i shift = _mm256_srli_epi32(
            _mm256_and_si256(_mm256_set1_epi32(x), _mm256_set1_epi32(0x3f)), 1);

        for (int32_t j = 0, jj = 0; j < 32; j += 16, jj++) {
            int32_t mdiff = base_max_diff - j;
            if (mdiff <= 0) {
                res16[jj] = a_mbase_x;
            }
            else {
                a0 = _mm256_cvtepu8_epi32(
                    _mm_loadu_si128((__m128i *)(above + base + j)));
                a1 = _mm256_cvtepu8_epi32(
                    _mm_loadu_si128((__m128i *)(above + base + 1 + j)));

                diff = _mm256_sub_epi32(a1, a0);   // a[x+1] - a[x]
                a32 = _mm256_slli_epi32(a0, 5);    // a[x] * 32
                a32 = _mm256_add_epi32(a32, a16);  // a[x] * 32 + 16
                b = _mm256_mullo_epi32(diff, shift);

                res[0] = _mm256_add_epi32(a32, b);
                res[0] = _mm256_srli_epi32(res[0], 5);
                res[0] = _mm256_packus_epi32(
                    res[0],
                    _mm256_castsi128_si256(_mm256_extracti128_si256(res[0], 1)));

                // goto 8 bit
                res[0] = _mm256_packus_epi16(res[0], res[0]);

                if (mdiff > 8) {
                    a0_1 = _mm256_cvtepu8_epi32(
                        _mm_loadu_si128((__m128i *)(above + base + 8 + j)));
                    a1_1 = _mm256_cvtepu8_epi32(
                        _mm_loadu_si128((__m128i *)(above + base + 9 + j)));

                    diff = _mm256_sub_epi32(a1_1, a0_1);  // a[x+1] - a[x]
                    a32 = _mm256_slli_epi32(a0_1, 5);     // a[x] * 32
                    a32 = _mm256_add_epi32(a32, a16);     // a[x] * 32 + 16
                    b = _mm256_mullo_epi32(diff, shift);

                    res[1] = _mm256_add_epi32(a32, b);
                    res[1] = _mm256_srli_epi32(res[1], 5);
                    res[1] = _mm256_packus_epi32(
                        res[1],
                        _mm256_castsi128_si256(_mm256_extracti128_si256(res[1], 1)));
                    res[1] = _mm256_packus_epi16(res[1], res[1]);
                    // goto 8 bit
                }
                else {
                    res[1] = a_mbase_x;
                }
                res16[jj] = _mm256_unpacklo_epi64(res[0], res[1]);  // 16 8bit values
            }
        }
        res16[1] =
            _mm256_inserti128_si256(res16[0], _mm256_castsi256_si128(res16[1]),
            1);  // 32 8bit values

        dstvec[r] = _mm256_blendv_epi8(
            a_mbase_x, res16[1],
            *(__m256i *)BaseMask[base_max_diff]);  // 32 8bit values
        x += dx;
    }
}

static void dr_prediction_z1_32xN_avx2(int32_t N, uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above, int32_t upsample_above,
    int32_t dx) {
    __m256i dstvec[64];
    dr_prediction_z1_32xN_internal_avx2(N, dstvec, above, upsample_above, dx);
    for (int32_t i = 0; i < N; i++) {
        _mm256_storeu_si256((__m256i *)(dst + stride * i), dstvec[i]);
    }
}

static void dr_prediction_z1_64xN_avx2(int32_t N, uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above, int32_t upsample_above,
    int32_t dx) {
    int32_t x;

    // here upsample_above is 0 by design of av1_use_intra_edge_upsample
    (void)upsample_above;
    const int32_t frac_bits = 6;
    const int32_t max_base_x = ((64 + N) - 1);

    // pre-filter above pixels
    // store in temp buffers:
    //   above[x] * 32 + 16
    //   above[x+1] - above[x]
    // final pixels will be caluculated as:
    //   (above[x] * 32 + 16 + (above[x+1] - above[x]) * shift) >> 5
    __m256i a0, a0_1, a1, a1_1, a32, a16;
    __m256i a_mbase_x, diff;
    __m128i max_base_x128, base_inc128, mask128;

    a16 = _mm256_set1_epi32(16);
    a_mbase_x = _mm256_set1_epi8(above[max_base_x]);
    max_base_x128 = _mm_set1_epi8(max_base_x);

    x = dx;
    for (int32_t r = 0; r < N; r++, dst += stride) {
        __m256i b, res[2];
        __m128i res1;

        int32_t base = x >> frac_bits;
        if (base >= max_base_x) {
            for (int32_t i = r; i < N; ++i) {
                _mm256_storeu_si256((__m256i *)dst, a_mbase_x);  // save 32 values
                _mm256_storeu_si256((__m256i *)(dst + 32), a_mbase_x);
                dst += stride;
            }
            return;
        }

        __m256i shift = _mm256_srli_epi32(
            _mm256_and_si256(_mm256_set1_epi32(x), _mm256_set1_epi32(0x3f)), 1);

        __m128i a0_128, a0_1_128, a1_128, a1_1_128;
        for (int32_t j = 0; j < 64; j += 16) {
            int32_t mdif = max_base_x - (base + j);
            if (mdif <= 0) {
                _mm_storeu_si128((__m128i *)(dst + j),
                    _mm256_castsi256_si128(a_mbase_x));
            }
            else {
                a0_128 = _mm_loadu_si128((__m128i *)(above + base + j));
                a1_128 = _mm_loadu_si128((__m128i *)(above + base + 1 + j));
                a0 = _mm256_cvtepu8_epi32(a0_128);
                a1 = _mm256_cvtepu8_epi32(a1_128);

                diff = _mm256_sub_epi32(a1, a0);   // a[x+1] - a[x]
                a32 = _mm256_slli_epi32(a0, 5);    // a[x] * 32
                a32 = _mm256_add_epi32(a32, a16);  // a[x] * 32 + 16
                b = _mm256_mullo_epi32(diff, shift);

                res[0] = _mm256_add_epi32(a32, b);
                res[0] = _mm256_srli_epi32(res[0], 5);
                res[0] = _mm256_packus_epi32(
                    res[0],
                    _mm256_castsi128_si256(_mm256_extracti128_si256(res[0], 1)));
                // goto 8 bit
                res[0] = _mm256_packus_epi16(res[0], res[0]);

                if (mdif > 8) {
                    a0_1_128 = _mm_loadu_si128((__m128i *)(above + base + 8 + j));
                    a1_1_128 = _mm_loadu_si128((__m128i *)(above + base + 9 + j));
                    a0_1 = _mm256_cvtepu8_epi32(a0_1_128);
                    a1_1 = _mm256_cvtepu8_epi32(a1_1_128);

                    diff = _mm256_sub_epi32(a1_1, a0_1);  // a[x+1] - a[x]
                    a32 = _mm256_slli_epi32(a0_1, 5);     // a[x] * 32
                    a32 = _mm256_add_epi32(a32, a16);     // a[x] * 32 + 16
                    b = _mm256_mullo_epi32(diff, shift);

                    res[1] = _mm256_add_epi32(a32, b);
                    res[1] = _mm256_srli_epi32(res[1], 5);
                    res[1] = _mm256_packus_epi32(
                        res[1],
                        _mm256_castsi128_si256(_mm256_extracti128_si256(res[1], 1)));
                    res[1] = _mm256_packus_epi16(res[1], res[1]);

                }
                else {
                    res[1] = a_mbase_x;
                }
                res1 = _mm_unpacklo_epi64(
                    _mm256_castsi256_si128(res[0]),
                    _mm256_castsi256_si128(res[1]));  // 16 8bit values

                base_inc128 = _mm_setr_epi8(
                    base + j, base + j + 1, base + j + 2, base + j + 3, base + j + 4,
                    base + j + 5, base + j + 6, base + j + 7, base + j + 8,
                    base + j + 9, base + j + 10, base + j + 11, base + j + 12,
                    base + j + 13, base + j + 14, base + j + 15);

                mask128 = _mm_cmpgt_epi8(_mm_subs_epu8(max_base_x128, base_inc128),
                    _mm_setzero_si128());
                res1 =
                    _mm_blendv_epi8(_mm256_castsi256_si128(a_mbase_x), res1, mask128);
                _mm_storeu_si128((__m128i *)(dst + j), res1);
            }
        }
        x += dx;
    }
}

// Directional prediction, zone 1: 0 < angle < 90
void av1_dr_prediction_z1_avx2(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh,
    const uint8_t *above, const uint8_t *left,
    int32_t upsample_above, int32_t dx, int32_t dy) {
    (void)left;
    (void)dy;
    switch (bw) {
    case 4:
        dr_prediction_z1_4xN_avx2(bh, dst, stride, above, upsample_above, dx);
        break;
    case 8:
        dr_prediction_z1_8xN_avx2(bh, dst, stride, above, upsample_above, dx);
        break;
    case 16:
        dr_prediction_z1_16xN_avx2(bh, dst, stride, above, upsample_above, dx);
        break;
    case 32:
        dr_prediction_z1_32xN_avx2(bh, dst, stride, above, upsample_above, dx);
        break;
    case 64:
        dr_prediction_z1_64xN_avx2(bh, dst, stride, above, upsample_above, dx);
        break;
    default: break;
    }
    return;
}

static AOM_FORCE_INLINE void highbd_dr_prediction_z1_4xN_internal_avx2(
    int32_t N, __m128i *dst, const uint16_t *above, int32_t upsample_above, int32_t dx) {
    const int32_t frac_bits = 6 - upsample_above;
    const int32_t max_base_x = ((N + 4) - 1) << upsample_above;
    int32_t x;
    // a assert(dx > 0);
    // pre-filter above pixels
    // store in temp buffers:
    //   above[x] * 32 + 16
    //   above[x+1] - above[x]
    // final pixels will be caluculated as:
    //   (above[x] * 32 + 16 + (above[x+1] - above[x]) * shift) >> 5
    __m256i a0, a1, a32, a16;
    __m256i diff;
    __m128i a_mbase_x, max_base_x128, base_inc128, mask128;

    a16 = _mm256_set1_epi32(16);
    a_mbase_x = _mm_set1_epi16(above[max_base_x]);
    max_base_x128 = _mm_set1_epi32(max_base_x);

    x = dx;
    for (int32_t r = 0; r < N; r++) {
        __m256i b, res, shift;
        __m128i res1;

        int32_t base = x >> frac_bits;
        if (base >= max_base_x) {
            for (int32_t i = r; i < N; ++i) {
                dst[i] = a_mbase_x;  // save 4 values
            }
            return;
        }

        a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(above + base)));
        a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(above + base + 1)));

        if (upsample_above) {
            a0 = _mm256_permutevar8x32_epi32(
                a0, _mm256_set_epi32(7, 5, 3, 1, 6, 4, 2, 0));
            a1 = _mm256_castsi128_si256(_mm256_extracti128_si256(a0, 1));
            base_inc128 = _mm_setr_epi32(base, base + 2, base + 4, base + 6);
            shift = _mm256_srli_epi32(
                _mm256_and_si256(
                _mm256_slli_epi32(_mm256_set1_epi32(x), upsample_above),
                _mm256_set1_epi32(0x3f)),
                1);
        }
        else {
            base_inc128 = _mm_setr_epi32(base, base + 1, base + 2, base + 3);
            shift = _mm256_srli_epi32(
                _mm256_and_si256(_mm256_set1_epi32(x), _mm256_set1_epi32(0x3f)), 1);
        }

        diff = _mm256_sub_epi32(a1, a0);   // a[x+1] - a[x]
        a32 = _mm256_slli_epi32(a0, 5);    // a[x] * 32
        a32 = _mm256_add_epi32(a32, a16);  // a[x] * 32 + 16

        b = _mm256_mullo_epi32(diff, shift);
        res = _mm256_add_epi32(a32, b);
        res = _mm256_srli_epi32(res, 5);

        res1 = _mm256_castsi256_si128(res);
        res1 = _mm_packus_epi32(res1, res1);

        mask128 = _mm_cmpgt_epi32(max_base_x128, base_inc128);
        mask128 = _mm_packs_epi32(mask128, mask128);  // goto 16 bit
        dst[r] = _mm_blendv_epi8(a_mbase_x, res1, mask128);
        x += dx;
    }
}

static void highbd_dr_prediction_z1_4xN_avx2(int32_t N, uint16_t *dst,
    ptrdiff_t stride,
    const uint16_t *above,
    int32_t upsample_above, int32_t dx) {
    __m128i dstvec[16];

    highbd_dr_prediction_z1_4xN_internal_avx2(N, dstvec, above, upsample_above,
        dx);
    for (int32_t i = 0; i < N; i++) {
        _mm_storel_epi64((__m128i *)(dst + stride * i), dstvec[i]);
    }
}

static AOM_FORCE_INLINE void highbd_dr_prediction_z1_8xN_internal_avx2(
    int32_t N, __m128i *dst, const uint16_t *above, int32_t upsample_above, int32_t dx) {
    const int32_t frac_bits = 6 - upsample_above;
    const int32_t max_base_x = ((8 + N) - 1) << upsample_above;

    int32_t x;
    // a assert(dx > 0);
    // pre-filter above pixels
    // store in temp buffers:
    //   above[x] * 32 + 16
    //   above[x+1] - above[x]
    // final pixels will be caluculated as:
    //   (above[x] * 32 + 16 + (above[x+1] - above[x]) * shift) >> 5
    __m256i a0, a1, a0_1, a1_1, a32, a16;
    __m256i a_mbase_x, diff, max_base_x256, base_inc256, mask256;

    a16 = _mm256_set1_epi32(16);
    a_mbase_x = _mm256_set1_epi16(above[max_base_x]);
    max_base_x256 = _mm256_set1_epi32(max_base_x);

    x = dx;
    for (int32_t r = 0; r < N; r++) {
        __m256i b, res, res1, shift;

        int32_t base = x >> frac_bits;
        if (base >= max_base_x) {
            for (int32_t i = r; i < N; ++i) {
                dst[i] = _mm256_castsi256_si128(a_mbase_x);  // save 8 values
            }
            return;
        }

        a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(above + base)));
        a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(above + base + 1)));

        if (upsample_above) {
            a0 = _mm256_permutevar8x32_epi32(
                a0, _mm256_set_epi32(7, 5, 3, 1, 6, 4, 2, 0));
            a1 = _mm256_castsi128_si256(_mm256_extracti128_si256(a0, 1));

            a0_1 =
                _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(above + base + 8)));
            a0_1 = _mm256_permutevar8x32_epi32(
                a0_1, _mm256_set_epi32(7, 5, 3, 1, 6, 4, 2, 0));
            a1_1 = _mm256_castsi128_si256(_mm256_extracti128_si256(a0_1, 1));

            a0 = _mm256_inserti128_si256(a0, _mm256_castsi256_si128(a0_1), 1);
            a1 = _mm256_inserti128_si256(a1, _mm256_castsi256_si128(a1_1), 1);
            base_inc256 =
                _mm256_setr_epi32(base, base + 2, base + 4, base + 6, base + 8,
                base + 10, base + 12, base + 14);
            shift = _mm256_srli_epi32(
                _mm256_and_si256(
                _mm256_slli_epi32(_mm256_set1_epi32(x), upsample_above),
                _mm256_set1_epi32(0x3f)),
                1);
        }
        else {
            base_inc256 = _mm256_setr_epi32(base, base + 1, base + 2, base + 3,
                base + 4, base + 5, base + 6, base + 7);
            shift = _mm256_srli_epi32(
                _mm256_and_si256(_mm256_set1_epi32(x), _mm256_set1_epi32(0x3f)), 1);
        }

        diff = _mm256_sub_epi32(a1, a0);   // a[x+1] - a[x]
        a32 = _mm256_slli_epi32(a0, 5);    // a[x] * 32
        a32 = _mm256_add_epi32(a32, a16);  // a[x] * 32 + 16

        b = _mm256_mullo_epi32(diff, shift);
        res = _mm256_add_epi32(a32, b);
        res = _mm256_srli_epi32(res, 5);

        res1 = _mm256_packus_epi32(
            res, _mm256_castsi128_si256(_mm256_extracti128_si256(res, 1)));

        mask256 = _mm256_cmpgt_epi32(max_base_x256, base_inc256);
        mask256 = _mm256_packs_epi32(
            mask256, _mm256_castsi128_si256(
            _mm256_extracti128_si256(mask256, 1)));  // goto 16 bit
        res1 = _mm256_blendv_epi8(a_mbase_x, res1, mask256);
        dst[r] = _mm256_castsi256_si128(res1);
        x += dx;
    }
}

static void highbd_dr_prediction_z1_8xN_avx2(int32_t N, uint16_t *dst,
    ptrdiff_t stride,
    const uint16_t *above,
    int32_t upsample_above, int32_t dx) {
    __m128i dstvec[32];

    highbd_dr_prediction_z1_8xN_internal_avx2(N, dstvec, above, upsample_above,
        dx);
    for (int32_t i = 0; i < N; i++) {
        _mm_storeu_si128((__m128i *)(dst + stride * i), dstvec[i]);
    }
}

static AOM_FORCE_INLINE void highbd_dr_prediction_z1_16xN_internal_avx2(
    int32_t N, __m256i *dstvec, const uint16_t *above, int32_t upsample_above, int32_t dx) {
    int32_t x;
    // here upsample_above is 0 by design of av1_use_intra_edge_upsample
    (void)upsample_above;
    const int32_t frac_bits = 6;
    const int32_t max_base_x = ((16 + N) - 1);

    // pre-filter above pixels
    // store in temp buffers:
    //   above[x] * 32 + 16
    //   above[x+1] - above[x]
    // final pixels will be caluculated as:
    //   (above[x] * 32 + 16 + (above[x+1] - above[x]) * shift) >> 5
    __m256i a0, a0_1, a1, a1_1, a32, a16;
    __m256i a_mbase_x, diff, max_base_x256, base_inc256, mask256;

    a16 = _mm256_set1_epi32(16);
    a_mbase_x = _mm256_set1_epi16(above[max_base_x]);
    max_base_x256 = _mm256_set1_epi16(max_base_x);

    x = dx;
    for (int32_t r = 0; r < N; r++) {
        __m256i b, res[2], res1;

        int32_t base = x >> frac_bits;
        if (base >= max_base_x) {
            for (int32_t i = r; i < N; ++i) {
                dstvec[i] = a_mbase_x;  // save 16 values
            }
            return;
        }
        __m256i shift = _mm256_srli_epi32(
            _mm256_and_si256(_mm256_set1_epi32(x), _mm256_set1_epi32(0x3f)), 1);

        a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(above + base)));
        a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(above + base + 1)));

        diff = _mm256_sub_epi32(a1, a0);   // a[x+1] - a[x]
        a32 = _mm256_slli_epi32(a0, 5);    // a[x] * 32
        a32 = _mm256_add_epi32(a32, a16);  // a[x] * 32 + 16
        b = _mm256_mullo_epi32(diff, shift);

        res[0] = _mm256_add_epi32(a32, b);
        res[0] = _mm256_srli_epi32(res[0], 5);
        res[0] = _mm256_packus_epi32(
            res[0], _mm256_castsi128_si256(_mm256_extracti128_si256(res[0], 1)));

        int32_t mdif = max_base_x - base;
        if (mdif > 8) {
            a0_1 =
                _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(above + base + 8)));
            a1_1 =
                _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(above + base + 9)));

            diff = _mm256_sub_epi32(a1_1, a0_1);  // a[x+1] - a[x]
            a32 = _mm256_slli_epi32(a0_1, 5);     // a[x] * 32
            a32 = _mm256_add_epi32(a32, a16);     // a[x] * 32 + 16
            b = _mm256_mullo_epi32(diff, shift);

            res[1] = _mm256_add_epi32(a32, b);
            res[1] = _mm256_srli_epi32(res[1], 5);
            res[1] = _mm256_packus_epi32(
                res[1], _mm256_castsi128_si256(_mm256_extracti128_si256(res[1], 1)));
        }
        else {
            res[1] = a_mbase_x;
        }
        res1 = _mm256_inserti128_si256(res[0], _mm256_castsi256_si128(res[1]),
            1);  // 16 16bit values

        base_inc256 = _mm256_setr_epi16(base, base + 1, base + 2, base + 3,
            base + 4, base + 5, base + 6, base + 7,
            base + 8, base + 9, base + 10, base + 11,
            base + 12, base + 13, base + 14, base + 15);
        mask256 = _mm256_cmpgt_epi16(max_base_x256, base_inc256);
        dstvec[r] = _mm256_blendv_epi8(a_mbase_x, res1, mask256);
        x += dx;
    }
}

static void highbd_dr_prediction_z1_16xN_avx2(int32_t N, uint16_t *dst,
    ptrdiff_t stride,
    const uint16_t *above,
    int32_t upsample_above, int32_t dx) {
    __m256i dstvec[64];
    highbd_dr_prediction_z1_16xN_internal_avx2(N, dstvec, above, upsample_above,
        dx);
    for (int32_t i = 0; i < N; i++) {
        _mm256_storeu_si256((__m256i *)(dst + stride * i), dstvec[i]);
    }
}

static AOM_FORCE_INLINE void highbd_dr_prediction_z1_32xN_internal_avx2(
    int32_t N, __m256i *dstvec, const uint16_t *above, int32_t upsample_above, int32_t dx) {
    int32_t x;
    // here upsample_above is 0 by design of av1_use_intra_edge_upsample
    (void)upsample_above;
    const int32_t frac_bits = 6;
    const int32_t max_base_x = ((32 + N) - 1);

    // pre-filter above pixels
    // store in temp buffers:
    //   above[x] * 32 + 16
    //   above[x+1] - above[x]
    // final pixels will be caluculated as:
    //   (above[x] * 32 + 16 + (above[x+1] - above[x]) * shift) >> 5
    __m256i a0, a0_1, a1, a1_1, a32, a16;
    __m256i a_mbase_x, diff, max_base_x256, base_inc256, mask256;

    a16 = _mm256_set1_epi32(16);
    a_mbase_x = _mm256_set1_epi16(above[max_base_x]);
    max_base_x256 = _mm256_set1_epi16(max_base_x);

    x = dx;
    for (int32_t r = 0; r < N; r++) {
        __m256i b, res[2], res1;

        int32_t base = x >> frac_bits;
        if (base >= max_base_x) {
            for (int32_t i = r; i < N; ++i) {
                dstvec[i] = a_mbase_x;  // save 32 values
                dstvec[i + N] = a_mbase_x;
            }
            return;
        }

        __m256i shift = _mm256_srli_epi32(
            _mm256_and_si256(_mm256_set1_epi32(x), _mm256_set1_epi32(0x3f)), 1);

        for (int32_t j = 0; j < 32; j += 16) {
            int32_t mdif = max_base_x - (base + j);
            if (mdif <= 0) {
                res1 = a_mbase_x;
            }
            else {
                a0 = _mm256_cvtepu16_epi32(
                    _mm_loadu_si128((__m128i *)(above + base + j)));
                a1 = _mm256_cvtepu16_epi32(
                    _mm_loadu_si128((__m128i *)(above + base + 1 + j)));

                diff = _mm256_sub_epi32(a1, a0);   // a[x+1] - a[x]
                a32 = _mm256_slli_epi32(a0, 5);    // a[x] * 32
                a32 = _mm256_add_epi32(a32, a16);  // a[x] * 32 + 16
                b = _mm256_mullo_epi32(diff, shift);

                res[0] = _mm256_add_epi32(a32, b);
                res[0] = _mm256_srli_epi32(res[0], 5);
                res[0] = _mm256_packus_epi32(
                    res[0],
                    _mm256_castsi128_si256(_mm256_extracti128_si256(res[0], 1)));
                if (mdif > 8) {
                    a0_1 = _mm256_cvtepu16_epi32(
                        _mm_loadu_si128((__m128i *)(above + base + 8 + j)));
                    a1_1 = _mm256_cvtepu16_epi32(
                        _mm_loadu_si128((__m128i *)(above + base + 9 + j)));

                    diff = _mm256_sub_epi32(a1_1, a0_1);  // a[x+1] - a[x]
                    a32 = _mm256_slli_epi32(a0_1, 5);     // a[x] * 32
                    a32 = _mm256_add_epi32(a32, a16);     // a[x] * 32 + 16
                    b = _mm256_mullo_epi32(diff, shift);

                    res[1] = _mm256_add_epi32(a32, b);
                    res[1] = _mm256_srli_epi32(res[1], 5);
                    res[1] = _mm256_packus_epi32(
                        res[1],
                        _mm256_castsi128_si256(_mm256_extracti128_si256(res[1], 1)));
                }
                else {
                    res[1] = a_mbase_x;
                }
                res1 = _mm256_inserti128_si256(res[0], _mm256_castsi256_si128(res[1]),
                    1);  // 16 16bit values
                base_inc256 = _mm256_setr_epi16(
                    base + j, base + j + 1, base + j + 2, base + j + 3, base + j + 4,
                    base + j + 5, base + j + 6, base + j + 7, base + j + 8,
                    base + j + 9, base + j + 10, base + j + 11, base + j + 12,
                    base + j + 13, base + j + 14, base + j + 15);

                mask256 = _mm256_cmpgt_epi16(max_base_x256, base_inc256);
                res1 = _mm256_blendv_epi8(a_mbase_x, res1, mask256);
            }
            if (!j)
                dstvec[r] = res1;
            else
                dstvec[r + N] = res1;
        }
        x += dx;
    }
}

static void highbd_dr_prediction_z1_32xN_avx2(int32_t N, uint16_t *dst,
    ptrdiff_t stride,
    const uint16_t *above,
    int32_t upsample_above, int32_t dx) {
    __m256i dstvec[128];

    highbd_dr_prediction_z1_32xN_internal_avx2(N, dstvec, above, upsample_above,
        dx);
    for (int32_t i = 0; i < N; i++) {
        _mm256_storeu_si256((__m256i *)(dst + stride * i), dstvec[i]);
        _mm256_storeu_si256((__m256i *)(dst + stride * i + 16), dstvec[i + N]);
    }
}

static void highbd_dr_prediction_z1_64xN_avx2(int32_t N, uint16_t *dst,
    ptrdiff_t stride,
    const uint16_t *above,
    int32_t upsample_above, int32_t dx) {
    int32_t x;

    // here upsample_above is 0 by design of av1_use_intra_edge_upsample
    (void)upsample_above;
    const int32_t frac_bits = 6;
    const int32_t max_base_x = ((64 + N) - 1);

    // pre-filter above pixels
    // store in temp buffers:
    //   above[x] * 32 + 16
    //   above[x+1] - above[x]
    // final pixels will be caluculated as:
    //   (above[x] * 32 + 16 + (above[x+1] - above[x]) * shift) >> 5
    __m256i a0, a0_1, a1, a1_1, a32, a16;
    __m256i a_mbase_x, diff, max_base_x256, base_inc256, mask256;

    a16 = _mm256_set1_epi32(16);
    a_mbase_x = _mm256_set1_epi16(above[max_base_x]);
    max_base_x256 = _mm256_set1_epi16(max_base_x);

    x = dx;
    for (int32_t r = 0; r < N; r++, dst += stride) {
        __m256i b, res[2], res1;

        int32_t base = x >> frac_bits;
        if (base >= max_base_x) {
            for (int32_t i = r; i < N; ++i) {
                _mm256_storeu_si256((__m256i *)dst, a_mbase_x);  // save 32 values
                _mm256_storeu_si256((__m256i *)(dst + 16), a_mbase_x);
                _mm256_storeu_si256((__m256i *)(dst + 32), a_mbase_x);
                _mm256_storeu_si256((__m256i *)(dst + 48), a_mbase_x);
                dst += stride;
            }
            return;
        }

        __m256i shift = _mm256_srli_epi32(
            _mm256_and_si256(_mm256_set1_epi32(x), _mm256_set1_epi32(0x3f)), 1);

        __m128i a0_128, a0_1_128, a1_128, a1_1_128;
        for (int32_t j = 0; j < 64; j += 16) {
            int32_t mdif = max_base_x - (base + j);
            if (mdif <= 0) {
                _mm256_storeu_si256((__m256i *)(dst + j), a_mbase_x);
            }
            else {
                a0_128 = _mm_loadu_si128((__m128i *)(above + base + j));
                a1_128 = _mm_loadu_si128((__m128i *)(above + base + 1 + j));
                a0 = _mm256_cvtepu16_epi32(a0_128);
                a1 = _mm256_cvtepu16_epi32(a1_128);

                diff = _mm256_sub_epi32(a1, a0);   // a[x+1] - a[x]
                a32 = _mm256_slli_epi32(a0, 5);    // a[x] * 32
                a32 = _mm256_add_epi32(a32, a16);  // a[x] * 32 + 16
                b = _mm256_mullo_epi32(diff, shift);

                res[0] = _mm256_add_epi32(a32, b);
                res[0] = _mm256_srli_epi32(res[0], 5);
                res[0] = _mm256_packus_epi32(
                    res[0],
                    _mm256_castsi128_si256(_mm256_extracti128_si256(res[0], 1)));
                if (mdif > 8) {
                    a0_1_128 = _mm_loadu_si128((__m128i *)(above + base + 8 + j));
                    a1_1_128 = _mm_loadu_si128((__m128i *)(above + base + 9 + j));
                    a0_1 = _mm256_cvtepu16_epi32(a0_1_128);
                    a1_1 = _mm256_cvtepu16_epi32(a1_1_128);

                    diff = _mm256_sub_epi32(a1_1, a0_1);  // a[x+1] - a[x]
                    a32 = _mm256_slli_epi32(a0_1, 5);     // a[x] * 32
                    a32 = _mm256_add_epi32(a32, a16);     // a[x] * 32 + 16
                    b = _mm256_mullo_epi32(diff, shift);

                    res[1] = _mm256_add_epi32(a32, b);
                    res[1] = _mm256_srli_epi32(res[1], 5);
                    res[1] = _mm256_packus_epi32(
                        res[1],
                        _mm256_castsi128_si256(_mm256_extracti128_si256(res[1], 1)));
                }
                else {
                    res[1] = a_mbase_x;
                }
                res1 = _mm256_inserti128_si256(res[0], _mm256_castsi256_si128(res[1]),
                    1);  // 16 16bit values
                base_inc256 = _mm256_setr_epi16(
                    base + j, base + j + 1, base + j + 2, base + j + 3, base + j + 4,
                    base + j + 5, base + j + 6, base + j + 7, base + j + 8,
                    base + j + 9, base + j + 10, base + j + 11, base + j + 12,
                    base + j + 13, base + j + 14, base + j + 15);

                mask256 = _mm256_cmpgt_epi16(max_base_x256, base_inc256);
                res1 = _mm256_blendv_epi8(a_mbase_x, res1, mask256);
                _mm256_storeu_si256((__m256i *)(dst + j), res1);
            }
        }
        x += dx;
    }
}

// Directional prediction, zone 1: 0 < angle < 90
void av1_highbd_dr_prediction_z1_avx2(uint16_t *dst, ptrdiff_t stride, int32_t bw,
    int32_t bh, const uint16_t *above,
    const uint16_t *left, int32_t upsample_above,
    int32_t dx, int32_t dy, int32_t bd) {
    (void)left;
    (void)dy;
    (void)bd;

    switch (bw) {
    case 4:
        highbd_dr_prediction_z1_4xN_avx2(bh, dst, stride, above, upsample_above,
            dx);
        break;
    case 8:
        highbd_dr_prediction_z1_8xN_avx2(bh, dst, stride, above, upsample_above,
            dx);
        break;
    case 16:
        highbd_dr_prediction_z1_16xN_avx2(bh, dst, stride, above, upsample_above,
            dx);
        break;
    case 32:
        highbd_dr_prediction_z1_32xN_avx2(bh, dst, stride, above, upsample_above,
            dx);
        break;
    case 64:
        highbd_dr_prediction_z1_64xN_avx2(bh, dst, stride, above, upsample_above,
            dx);
        break;
    default: break;
    }
    return;
}


static uint8_t LoadMaskx[8][16] = {
    { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 },
    { 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 },
    { 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 },
    { 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 },
    { 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 },
    { 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 },
    { 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8 },
};

static uint8_t EvenOddMaskx4[8][16] = {
    { 0, 2, 4, 6, 1, 3, 5, 7, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 1, 3, 5, 7, 2, 4, 6, 8, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 2, 4, 6, 8, 3, 5, 7, 9, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 3, 5, 7, 9, 4, 6, 8, 10, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 4, 6, 8, 10, 5, 7, 9, 11, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 5, 7, 9, 11, 6, 8, 10, 12, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 6, 8, 10, 12, 7, 9, 11, 13, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 7, 9, 11, 13, 8, 10, 12, 14, 0 }
};

static uint8_t EvenOddMaskx[8][16] = {
    { 0, 2, 4, 6, 8, 10, 12, 14, 1, 3, 5, 7, 0, 0, 0, 0 },
    { 0, 1, 3, 5, 7, 9, 11, 13, 15, 2, 4, 6, 8, 0, 0, 0 },
    { 0, 0, 2, 4, 6, 8, 10, 12, 14, 3, 5, 7, 9, 0, 0, 0 },
    { 0, 0, 0, 3, 5, 7, 9, 11, 13, 15, 4, 6, 8, 10, 0 },
    { 0, 0, 0, 0, 4, 6, 8, 10, 12, 14, 5, 7, 9, 11, 0, 0 },
    { 0, 0, 0, 0, 0, 5, 7, 9, 11, 13, 15, 6, 8, 10, 12, 0 },
    { 0, 0, 0, 0, 0, 0, 6, 8, 10, 12, 14, 7, 9, 11, 13, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 7, 9, 11, 13, 15, 8, 10, 12, 14 }
};

static void dr_prediction_z2_Nx4_avx2(int32_t N, uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above, const uint8_t *left,
    int32_t upsample_above, int32_t upsample_left,
    int32_t dx, int32_t dy) {
    const int32_t min_base_x = -(1 << upsample_above);
    const int32_t min_base_y = -(1 << upsample_left);
    const int32_t frac_bits_x = 6 - upsample_above;
    const int32_t frac_bits_y = 6 - upsample_left;

    // a assert(dx > 0);
    // pre-filter above pixels
    // store in temp buffers:
    //   above[x] * 32 + 16
    //   above[x+1] - above[x]
    // final pixels will be caluculated as:
    //   (above[x] * 32 + 16 + (above[x+1] - above[x]) * shift) >> 5
    __m256i a0_x, a1_x, a0_y, a1_y, a32, a16, diff;
    __m128i c3f, min_base_y128;

    a16 = _mm256_set1_epi32(16);
    c3f = _mm_set1_epi32(0x3f);
    min_base_y128 = _mm_set1_epi32(min_base_y);

    for (int32_t r = 0; r < N; r++) {
        __m256i b, res, shift;
        __m128i resx, resy, resxy;
        __m128i a0_x128, a1_x128;
        int32_t y = r + 1;
        int32_t base_x = (-y * dx) >> frac_bits_x;
        int32_t base_shift = 0;
        if (base_x < (min_base_x - 1)) {
            base_shift = (min_base_x - base_x - 1) >> upsample_above;
        }
        int32_t base_min_diff =
            (min_base_x - base_x + upsample_above) >> upsample_above;
        if (base_min_diff > 4) {
            base_min_diff = 4;
        }
        else {
            if (base_min_diff < 0) base_min_diff = 0;
        }

        if (base_shift > 3) {
            resx = _mm_setzero_si128();
        }
        else {
            a0_x128 = _mm_loadu_si128((__m128i *)(above + base_x + base_shift));
            if (upsample_above) {
                a0_x128 =
                    _mm_shuffle_epi8(a0_x128, *(__m128i *)EvenOddMaskx4[base_shift]);
                a1_x128 = _mm_srli_si128(a0_x128, 4);

                shift = _mm256_castsi128_si256(_mm_srli_epi32(
                    _mm_and_si128(
                    _mm_slli_epi32(
                    _mm_setr_epi32(-y * dx, (1 << 6) - y * dx,
                    (2 << 6) - y * dx, (3 << 6) - y * dx),
                    upsample_above),
                    c3f),
                    1));
            }
            else {
                a0_x128 = _mm_shuffle_epi8(a0_x128, *(__m128i *)LoadMaskx[base_shift]);
                a1_x128 = _mm_srli_si128(a0_x128, 1);
                shift = _mm256_castsi128_si256(_mm_srli_epi32(
                    _mm_and_si128(_mm_setr_epi32(-y * dx, (1 << 6) - y * dx,
                    (2 << 6) - y * dx, (3 << 6) - y * dx),
                    c3f),
                    1));
            }
            a0_x = _mm256_cvtepu8_epi32(a0_x128);
            a1_x = _mm256_cvtepu8_epi32(a1_x128);

            diff = _mm256_sub_epi32(a1_x, a0_x);  // a[x+1] - a[x]
            a32 = _mm256_slli_epi32(a0_x, 5);     // a[x] * 32
            a32 = _mm256_add_epi32(a32, a16);     // a[x] * 32 + 16

            b = _mm256_mullo_epi32(diff, shift);
            res = _mm256_add_epi32(a32, b);
            res = _mm256_srli_epi32(res, 5);

            resx = _mm256_castsi256_si128(res);
            resx = _mm_packus_epi32(resx, resx);
            resx = _mm_packus_epi16(resx, resx);
        }
        // y calc
        if (base_x < min_base_x) {
            DECLARE_ALIGNED(32, int32_t, base_y_c[4]);
            __m128i r6, c1234, dy128, y_c128, base_y_c128, mask128;
            r6 = _mm_set1_epi32(r << 6);
            dy128 = _mm_set1_epi32(dy);
            c1234 = _mm_setr_epi32(1, 2, 3, 4);
            y_c128 = _mm_sub_epi32(r6, _mm_mullo_epi32(c1234, dy128));
            base_y_c128 = _mm_srai_epi32(y_c128, frac_bits_y);
            mask128 = _mm_cmpgt_epi32(min_base_y128, base_y_c128);
            base_y_c128 = _mm_andnot_si128(mask128, base_y_c128);
            _mm_store_si128((__m128i *)base_y_c, base_y_c128);

            a0_y = _mm256_castsi128_si256(
                _mm_setr_epi32(left[base_y_c[0]], left[base_y_c[1]],
                left[base_y_c[2]], left[base_y_c[3]]));
            a1_y = _mm256_castsi128_si256(
                _mm_setr_epi32(left[base_y_c[0] + 1], left[base_y_c[1] + 1],
                left[base_y_c[2] + 1], left[base_y_c[3] + 1]));

            if (upsample_left) {
                shift = _mm256_castsi128_si256(_mm_srli_epi32(
                    _mm_and_si128(_mm_slli_epi32(y_c128, upsample_left), c3f), 1));
            }
            else {
                shift = _mm256_castsi128_si256(
                    _mm_srli_epi32(_mm_and_si128(y_c128, c3f), 1));
            }
            diff = _mm256_sub_epi32(a1_y, a0_y);  // a[x+1] - a[x]
            a32 = _mm256_slli_epi32(a0_y, 5);     // a[x] * 32
            a32 = _mm256_add_epi32(a32, a16);     // a[x] * 32 + 16

            b = _mm256_mullo_epi32(diff, shift);
            res = _mm256_add_epi32(a32, b);
            res = _mm256_srli_epi32(res, 5);

            resy = _mm256_castsi256_si128(res);
            resy = _mm_packus_epi32(resy, resy);
            resy = _mm_packus_epi16(resy, resy);
        }
        else {
            resy = resx;
        }
        resxy = _mm_blendv_epi8(resx, resy, *(__m128i *)BaseMask[base_min_diff]);
        *(uint32_t *)(dst) = _mm_cvtsi128_si32(resxy);
        dst += stride;
    }
}

static void dr_prediction_z2_Nx8_avx2(int32_t N, uint8_t *dst, ptrdiff_t stride,
    const uint8_t *above, const uint8_t *left,
    int32_t upsample_above, int32_t upsample_left,
    int32_t dx, int32_t dy) {
    const int32_t min_base_x = -(1 << upsample_above);
    const int32_t min_base_y = -(1 << upsample_left);
    const int32_t frac_bits_x = 6 - upsample_above;
    const int32_t frac_bits_y = 6 - upsample_left;

    // pre-filter above pixels
    // store in temp buffers:
    //   above[x] * 32 + 16
    //   above[x+1] - above[x]
    // final pixels will be caluculated as:
    //   (above[x] * 32 + 16 + (above[x+1] - above[x]) * shift) >> 5
    __m256i a0_x, a1_x, a0_y, a1_y, a32, a16, c3f;
    __m256i diff, min_base_y256;
    __m128i a0_x128, a1_x128;

    a16 = _mm256_set1_epi32(16);
    c3f = _mm256_set1_epi32(0x3f);
    min_base_y256 = _mm256_set1_epi32(min_base_y);

    for (int32_t r = 0; r < N; r++) {
        __m256i b, res, shift;
        __m128i resx, resy, resxy;

        int32_t y = r + 1;
        int32_t base_x = (-y * dx) >> frac_bits_x;
        int32_t base_shift = 0;
        if (base_x < (min_base_x - 1)) {
            base_shift = (min_base_x - base_x - 1) >> upsample_above;
        }
        int32_t base_min_diff =
            (min_base_x - base_x + upsample_above) >> upsample_above;
        if (base_min_diff > 8) {
            base_min_diff = 8;
        }
        else {
            if (base_min_diff < 0) base_min_diff = 0;
        }

        if (base_shift > 7) {
            resx = _mm_setzero_si128();
        }
        else {
            a0_x128 = _mm_loadu_si128((__m128i *)(above + base_x + base_shift));
            a1_x128 = _mm_loadu_si128((__m128i *)(above + base_x + 1 + base_shift));
            if (upsample_above) {
                a0_x128 =
                    _mm_shuffle_epi8(a0_x128, *(__m128i *)EvenOddMaskx[base_shift]);
                a1_x128 =
                    _mm_shuffle_epi8(a1_x128, *(__m128i *)EvenOddMaskx[base_shift]);

                shift = _mm256_srli_epi32(
                    _mm256_and_si256(
                    _mm256_slli_epi32(
                    _mm256_setr_epi32(-y * dx, (1 << 6) - y * dx,
                    (2 << 6) - y * dx, (3 << 6) - y * dx,
                    (4 << 6) - y * dx, (5 << 6) - y * dx,
                    (6 << 6) - y * dx, (7 << 6) - y * dx),
                    upsample_above),
                    c3f),
                    1);
            }
            else {
                a0_x128 = _mm_shuffle_epi8(a0_x128, *(__m128i *)LoadMaskx[base_shift]);
                a1_x128 = _mm_shuffle_epi8(a1_x128, *(__m128i *)LoadMaskx[base_shift]);

                shift = _mm256_srli_epi32(
                    _mm256_and_si256(
                    _mm256_setr_epi32(-y * dx, (1 << 6) - y * dx, (2 << 6) - y * dx,
                    (3 << 6) - y * dx, (4 << 6) - y * dx,
                    (5 << 6) - y * dx, (6 << 6) - y * dx,
                    (7 << 6) - y * dx),
                    c3f),
                    1);
            }
            a0_x = _mm256_cvtepu8_epi32(a0_x128);
            a1_x = _mm256_cvtepu8_epi32(a1_x128);

            diff = _mm256_sub_epi32(a1_x, a0_x);  // a[x+1] - a[x]
            a32 = _mm256_slli_epi32(a0_x, 5);     // a[x] * 32
            a32 = _mm256_add_epi32(a32, a16);     // a[x] * 32 + 16

            b = _mm256_mullo_epi32(diff, shift);
            res = _mm256_add_epi32(a32, b);
            res = _mm256_srli_epi32(res, 5);

            res = _mm256_packus_epi32(
                res, _mm256_castsi128_si256(_mm256_extracti128_si256(res, 1)));
            resx = _mm_packus_epi16(_mm256_castsi256_si128(res),
                _mm256_castsi256_si128(res));
        }
        // y calc
        if (base_x < min_base_x) {
            DECLARE_ALIGNED(32, int32_t, base_y_c[8]);
            __m256i r6, c256, dy256, y_c256, base_y_c256, mask256;
            r6 = _mm256_set1_epi32(r << 6);
            dy256 = _mm256_set1_epi32(dy);
            c256 = _mm256_setr_epi32(1, 2, 3, 4, 5, 6, 7, 8);
            y_c256 = _mm256_sub_epi32(r6, _mm256_mullo_epi32(c256, dy256));
            base_y_c256 = _mm256_srai_epi32(y_c256, frac_bits_y);
            mask256 = _mm256_cmpgt_epi32(min_base_y256, base_y_c256);
            base_y_c256 = _mm256_andnot_si256(mask256, base_y_c256);
            _mm256_store_si256((__m256i *)base_y_c, base_y_c256);

            a0_y = _mm256_setr_epi32(left[base_y_c[0]], left[base_y_c[1]],
                left[base_y_c[2]], left[base_y_c[3]],
                left[base_y_c[4]], left[base_y_c[5]],
                left[base_y_c[6]], left[base_y_c[7]]);
            a1_y = _mm256_setr_epi32(left[base_y_c[0] + 1], left[base_y_c[1] + 1],
                left[base_y_c[2] + 1], left[base_y_c[3] + 1],
                left[base_y_c[4] + 1], left[base_y_c[5] + 1],
                left[base_y_c[6] + 1], left[base_y_c[7] + 1]);

            if (upsample_left) {
                shift = _mm256_srli_epi32(
                    _mm256_and_si256(_mm256_slli_epi32(y_c256, upsample_left), c3f), 1);
            }
            else {
                shift = _mm256_srli_epi32(_mm256_and_si256(y_c256, c3f), 1);
            }
            diff = _mm256_sub_epi32(a1_y, a0_y);  // a[x+1] - a[x]
            a32 = _mm256_slli_epi32(a0_y, 5);     // a[x] * 32
            a32 = _mm256_add_epi32(a32, a16);     // a[x] * 32 + 16

            b = _mm256_mullo_epi32(diff, shift);
            res = _mm256_add_epi32(a32, b);
            res = _mm256_srli_epi32(res, 5);

            res = _mm256_packus_epi32(
                res, _mm256_castsi128_si256(_mm256_extracti128_si256(res, 1)));
            resy = _mm_packus_epi16(_mm256_castsi256_si128(res),
                _mm256_castsi256_si128(res));
        }
        else {
            resy = resx;
        }
        resxy = _mm_blendv_epi8(resx, resy, *(__m128i *)BaseMask[base_min_diff]);
        _mm_storel_epi64((__m128i *)(dst), resxy);
        dst += stride;
    }
}

static void dr_prediction_z2_HxW_avx2(int32_t H, int32_t W, uint8_t *dst,
    ptrdiff_t stride, const uint8_t *above,
    const uint8_t *left, int32_t upsample_above,
    int32_t upsample_left, int32_t dx, int32_t dy) {
    // here upsample_above and upsample_left are 0 by design of
    // av1_use_intra_edge_upsample
    const int32_t min_base_x = -1;
    const int32_t min_base_y = -1;
    (void)upsample_above;
    (void)upsample_left;
    const int32_t frac_bits_x = 6;
    const int32_t frac_bits_y = 6;

    // pre-filter above pixels
    // store in temp buffers:
    //   above[x] * 32 + 16
    //   above[x+1] - above[x]
    // final pixels will be caluculated as:
    //   (above[x] * 32 + 16 + (above[x+1] - above[x]) * shift) >> 5
    __m256i a0_x, a1_x, a0_y, a1_y, a32, a0_1_x, a1_1_x, a16;
    __m256i diff, min_base_y256, c3f;
    __m128i a0_x128, a1_x128, a0_1_x128, a1_1_x128;

    a16 = _mm256_set1_epi32(16);
    min_base_y256 = _mm256_set1_epi32(min_base_y);
    c3f = _mm256_set1_epi32(0x3f);

    for (int32_t r = 0; r < H; r++) {
        __m256i b, res, shift;
        __m128i resx[2], resy[2];
        __m128i resxy;
        for (int32_t j = 0; j < W; j += 16) {
            int32_t y = r + 1;
            int32_t base_x = (-y * dx) >> frac_bits_x;

            int32_t base_shift = 0;
            if ((base_x + j) < (min_base_x - 1)) {
                base_shift = (min_base_x - (base_x + j) - 1);
            }
            int32_t base_min_diff = (min_base_x - base_x - j);
            if (base_min_diff > 16) {
                base_min_diff = 16;
            }
            else {
                if (base_min_diff < 0) base_min_diff = 0;
            }
            if (base_shift > 7) {
                resx[0] = _mm_setzero_si128();
            }
            else {
                a0_x128 = _mm_loadu_si128((__m128i *)(above + base_x + base_shift + j));
                a1_x128 =
                    _mm_loadu_si128((__m128i *)(above + base_x + base_shift + 1 + j));
                a0_x128 = _mm_shuffle_epi8(a0_x128, *(__m128i *)LoadMaskx[base_shift]);
                a1_x128 = _mm_shuffle_epi8(a1_x128, *(__m128i *)LoadMaskx[base_shift]);

                a0_x = _mm256_cvtepu8_epi32(a0_x128);
                a1_x = _mm256_cvtepu8_epi32(a1_x128);

                shift = _mm256_srli_epi32(
                    _mm256_and_si256(
                    _mm256_setr_epi32(
                    ((0 + j) << 6) - y * dx, ((1 + j) << 6) - y * dx,
                    ((2 + j) << 6) - y * dx, ((3 + j) << 6) - y * dx,
                    ((4 + j) << 6) - y * dx, ((5 + j) << 6) - y * dx,
                    ((6 + j) << 6) - y * dx, ((7 + j) << 6) - y * dx),
                    c3f),
                    1);

                diff = _mm256_sub_epi32(a1_x, a0_x);  // a[x+1] - a[x]
                a32 = _mm256_slli_epi32(a0_x, 5);     // a[x] * 32
                a32 = _mm256_add_epi32(a32, a16);     // a[x] * 32 + 16

                b = _mm256_mullo_epi32(diff, shift);
                res = _mm256_add_epi32(a32, b);
                res = _mm256_srli_epi32(res, 5);

                res = _mm256_packus_epi32(
                    res, _mm256_castsi128_si256(_mm256_extracti128_si256(res, 1)));
                resx[0] = _mm_packus_epi16(_mm256_castsi256_si128(res),
                    _mm256_castsi256_si128(res));
            }
            base_shift = 0;
            if ((base_x + j + 8) < (min_base_x - 1)) {
                base_shift = (min_base_x - (base_x + j + 8) - 1);
            }
            if (base_shift > 7) {
                resx[1] = _mm_setzero_si128();
            }
            else {
                a0_1_x128 =
                    _mm_loadu_si128((__m128i *)(above + base_x + base_shift + 8 + j));
                a1_1_x128 =
                    _mm_loadu_si128((__m128i *)(above + base_x + base_shift + 9 + j));
                a0_1_x128 =
                    _mm_shuffle_epi8(a0_1_x128, *(__m128i *)LoadMaskx[base_shift]);
                a1_1_x128 =
                    _mm_shuffle_epi8(a1_1_x128, *(__m128i *)LoadMaskx[base_shift]);

                a0_1_x = _mm256_cvtepu8_epi32(a0_1_x128);
                a1_1_x = _mm256_cvtepu8_epi32(a1_1_x128);

                shift = _mm256_srli_epi32(
                    _mm256_and_si256(
                    _mm256_setr_epi32(
                    ((8 + j) << 6) - y * dx, ((9 + j) << 6) - y * dx,
                    ((10 + j) << 6) - y * dx, ((11 + j) << 6) - y * dx,
                    ((12 + j) << 6) - y * dx, ((13 + j) << 6) - y * dx,
                    ((14 + j) << 6) - y * dx, ((15 + j) << 6) - y * dx),
                    _mm256_set1_epi32(0x3f)),
                    1);

                diff = _mm256_sub_epi32(a1_1_x, a0_1_x);  // a[x+1] - a[x]
                a32 = _mm256_slli_epi32(a0_1_x, 5);       // a[x] * 32
                a32 = _mm256_add_epi32(a32, a16);         // a[x] * 32 + 16
                b = _mm256_mullo_epi32(diff, shift);

                res = _mm256_add_epi32(a32, b);
                res = _mm256_srli_epi32(res, 5);
                res = _mm256_packus_epi32(
                    res, _mm256_castsi128_si256(_mm256_extracti128_si256(res, 1)));
                resx[1] = _mm_packus_epi16(_mm256_castsi256_si128(res),
                    _mm256_castsi256_si128(res));
            }
            resx[0] = _mm_unpacklo_epi64(resx[0], resx[1]);

            // y calc
            if ((base_x < min_base_x)) {
                DECLARE_ALIGNED(32, int32_t, base_y_c[16]);
                __m256i r6, c256, dy256, y_c256, y_c_1_256, base_y_c256, mask256;
                r6 = _mm256_set1_epi32(r << 6);
                dy256 = _mm256_set1_epi32(dy);
                c256 = _mm256_setr_epi32(1 + j, 2 + j, 3 + j, 4 + j, 5 + j, 6 + j,
                    7 + j, 8 + j);
                y_c256 = _mm256_sub_epi32(r6, _mm256_mullo_epi32(c256, dy256));
                base_y_c256 = _mm256_srai_epi32(y_c256, frac_bits_y);
                mask256 = _mm256_cmpgt_epi32(min_base_y256, base_y_c256);
                base_y_c256 = _mm256_andnot_si256(mask256, base_y_c256);
                _mm256_store_si256((__m256i *)base_y_c, base_y_c256);
                c256 = _mm256_setr_epi32(9 + j, 10 + j, 11 + j, 12 + j, 13 + j, 14 + j,
                    15 + j, 16 + j);
                y_c_1_256 = _mm256_sub_epi32(r6, _mm256_mullo_epi32(c256, dy256));
                base_y_c256 = _mm256_srai_epi32(y_c_1_256, frac_bits_y);
                mask256 = _mm256_cmpgt_epi32(min_base_y256, base_y_c256);
                base_y_c256 = _mm256_andnot_si256(mask256, base_y_c256);
                _mm256_store_si256((__m256i *)(base_y_c + 8), base_y_c256);

                a0_y = _mm256_cvtepu8_epi32(_mm_setr_epi8(
                    left[base_y_c[0]], left[base_y_c[1]], left[base_y_c[2]],
                    left[base_y_c[3]], left[base_y_c[4]], left[base_y_c[5]],
                    left[base_y_c[6]], left[base_y_c[7]], 0, 0, 0, 0, 0, 0, 0, 0));
                a1_y = _mm256_cvtepu8_epi32(_mm_setr_epi8(
                    left[base_y_c[0] + 1], left[base_y_c[1] + 1], left[base_y_c[2] + 1],
                    left[base_y_c[3] + 1], left[base_y_c[4] + 1], left[base_y_c[5] + 1],
                    left[base_y_c[6] + 1], left[base_y_c[7] + 1], 0, 0, 0, 0, 0, 0, 0,
                    0));

                shift = _mm256_srli_epi32(_mm256_and_si256(y_c256, c3f), 1);

                diff = _mm256_sub_epi32(a1_y, a0_y);  // a[x+1] - a[x]
                a32 = _mm256_slli_epi32(a0_y, 5);     // a[x] * 32
                a32 = _mm256_add_epi32(a32, a16);     // a[x] * 32 + 16

                b = _mm256_mullo_epi32(diff, shift);
                res = _mm256_add_epi32(a32, b);
                res = _mm256_srli_epi32(res, 5);

                res = _mm256_packus_epi32(
                    res, _mm256_castsi128_si256(_mm256_extracti128_si256(res, 1)));
                resy[0] = _mm_packus_epi16(_mm256_castsi256_si128(res),
                    _mm256_castsi256_si128(res));

                a0_y = _mm256_cvtepu8_epi32(_mm_setr_epi8(
                    left[base_y_c[8]], left[base_y_c[9]], left[base_y_c[10]],
                    left[base_y_c[11]], left[base_y_c[12]], left[base_y_c[13]],
                    left[base_y_c[14]], left[base_y_c[15]], 0, 0, 0, 0, 0, 0, 0, 0));
                a1_y = _mm256_cvtepu8_epi32(
                    _mm_setr_epi8(left[base_y_c[8] + 1], left[base_y_c[9] + 1],
                    left[base_y_c[10] + 1], left[base_y_c[11] + 1],
                    left[base_y_c[12] + 1], left[base_y_c[13] + 1],
                    left[base_y_c[14] + 1], left[base_y_c[15] + 1], 0, 0,
                    0, 0, 0, 0, 0, 0));
                shift = _mm256_srli_epi32(_mm256_and_si256(y_c_1_256, c3f), 1);

                diff = _mm256_sub_epi32(a1_y, a0_y);  // a[x+1] - a[x]
                a32 = _mm256_slli_epi32(a0_y, 5);     // a[x] * 32
                a32 = _mm256_add_epi32(a32, a16);     // a[x] * 32 + 16

                b = _mm256_mullo_epi32(diff, shift);
                res = _mm256_add_epi32(a32, b);
                res = _mm256_srli_epi32(res, 5);

                res = _mm256_packus_epi32(
                    res, _mm256_castsi128_si256(_mm256_extracti128_si256(res, 1)));
                resy[1] =
                    _mm_packus_epi16(_mm256_castsi256_si128(res),
                    _mm256_castsi256_si128(res));  // 8 16bit values
                resy[0] = _mm_unpacklo_epi64(resy[0], resy[1]);
            }
            else {
                resy[0] = resx[0];
            }
            resxy = _mm_blendv_epi8(resx[0], resy[0],
                *(__m128i *)BaseMask[base_min_diff]);
            _mm_storeu_si128((__m128i *)(dst + j), resxy);
        }  // for j
        dst += stride;
    }
}

// Directional prediction, zone 2: 90 < angle < 180
void av1_dr_prediction_z2_avx2(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh,
    const uint8_t *above, const uint8_t *left,
    int32_t upsample_above, int32_t upsample_left, int32_t dx,
    int32_t dy) {
    assert(dx > 0);
    assert(dy > 0);
    switch (bw) {
    case 4:
        dr_prediction_z2_Nx4_avx2(bh, dst, stride, above, left, upsample_above,
            upsample_left, dx, dy);
        break;
    case 8:
        dr_prediction_z2_Nx8_avx2(bh, dst, stride, above, left, upsample_above,
            upsample_left, dx, dy);

        break;
    default:
        dr_prediction_z2_HxW_avx2(bh, bw, dst, stride, above, left,
            upsample_above, upsample_left, dx, dy);
        break;
    }
    return;
}

// z3 functions
static INLINE void transpose4x16_sse2(__m128i *x, __m128i *d) {
    __m128i w0, w1, w2, w3, ww0, ww1, ww2, ww3;
    w0 = _mm_unpacklo_epi8(x[0], x[1]);
    w1 = _mm_unpacklo_epi8(x[2], x[3]);
    w2 = _mm_unpackhi_epi8(x[0], x[1]);
    w3 = _mm_unpackhi_epi8(x[2], x[3]);

    ww0 = _mm_unpacklo_epi16(w0, w1);
    ww1 = _mm_unpacklo_epi16(w2, w3);
    ww2 = _mm_unpackhi_epi16(w0, w1);
    ww3 = _mm_unpackhi_epi16(w2, w3);

    w0 = _mm_unpacklo_epi32(ww0, ww1);
    w2 = _mm_unpacklo_epi32(ww2, ww3);
    w1 = _mm_unpackhi_epi32(ww0, ww1);
    w3 = _mm_unpackhi_epi32(ww2, ww3);

    d[0] = _mm_unpacklo_epi64(w0, w2);
    d[1] = _mm_unpackhi_epi64(w0, w2);
    d[2] = _mm_unpacklo_epi64(w1, w3);
    d[3] = _mm_unpackhi_epi64(w1, w3);

    d[4] = _mm_srli_si128(d[0], 8);
    d[5] = _mm_srli_si128(d[1], 8);
    d[6] = _mm_srli_si128(d[2], 8);
    d[7] = _mm_srli_si128(d[3], 8);

    d[8] = _mm_srli_si128(d[0], 4);
    d[9] = _mm_srli_si128(d[1], 4);
    d[10] = _mm_srli_si128(d[2], 4);
    d[11] = _mm_srli_si128(d[3], 4);

    d[12] = _mm_srli_si128(d[0], 12);
    d[13] = _mm_srli_si128(d[1], 12);
    d[14] = _mm_srli_si128(d[2], 12);
    d[15] = _mm_srli_si128(d[3], 12);
}

static INLINE void transpose16x32_avx2(__m256i *x, __m256i *d) {
    __m256i w0, w1, w2, w3, w4, w5, w6, w7, w8, w9;
    __m256i w10, w11, w12, w13, w14, w15;

    w0 = _mm256_unpacklo_epi8(x[0], x[1]);
    w1 = _mm256_unpacklo_epi8(x[2], x[3]);
    w2 = _mm256_unpacklo_epi8(x[4], x[5]);
    w3 = _mm256_unpacklo_epi8(x[6], x[7]);

    w8 = _mm256_unpacklo_epi8(x[8], x[9]);
    w9 = _mm256_unpacklo_epi8(x[10], x[11]);
    w10 = _mm256_unpacklo_epi8(x[12], x[13]);
    w11 = _mm256_unpacklo_epi8(x[14], x[15]);

    w4 = _mm256_unpacklo_epi16(w0, w1);
    w5 = _mm256_unpacklo_epi16(w2, w3);
    w12 = _mm256_unpacklo_epi16(w8, w9);
    w13 = _mm256_unpacklo_epi16(w10, w11);

    w6 = _mm256_unpacklo_epi32(w4, w5);
    w7 = _mm256_unpackhi_epi32(w4, w5);
    w14 = _mm256_unpacklo_epi32(w12, w13);
    w15 = _mm256_unpackhi_epi32(w12, w13);

    // Store first 4-line result
    d[0] = _mm256_unpacklo_epi64(w6, w14);
    d[1] = _mm256_unpackhi_epi64(w6, w14);
    d[2] = _mm256_unpacklo_epi64(w7, w15);
    d[3] = _mm256_unpackhi_epi64(w7, w15);

    w4 = _mm256_unpackhi_epi16(w0, w1);
    w5 = _mm256_unpackhi_epi16(w2, w3);
    w12 = _mm256_unpackhi_epi16(w8, w9);
    w13 = _mm256_unpackhi_epi16(w10, w11);

    w6 = _mm256_unpacklo_epi32(w4, w5);
    w7 = _mm256_unpackhi_epi32(w4, w5);
    w14 = _mm256_unpacklo_epi32(w12, w13);
    w15 = _mm256_unpackhi_epi32(w12, w13);

    // Store second 4-line result
    d[4] = _mm256_unpacklo_epi64(w6, w14);
    d[5] = _mm256_unpackhi_epi64(w6, w14);
    d[6] = _mm256_unpacklo_epi64(w7, w15);
    d[7] = _mm256_unpackhi_epi64(w7, w15);

    // upper half
    w0 = _mm256_unpackhi_epi8(x[0], x[1]);
    w1 = _mm256_unpackhi_epi8(x[2], x[3]);
    w2 = _mm256_unpackhi_epi8(x[4], x[5]);
    w3 = _mm256_unpackhi_epi8(x[6], x[7]);

    w8 = _mm256_unpackhi_epi8(x[8], x[9]);
    w9 = _mm256_unpackhi_epi8(x[10], x[11]);
    w10 = _mm256_unpackhi_epi8(x[12], x[13]);
    w11 = _mm256_unpackhi_epi8(x[14], x[15]);

    w4 = _mm256_unpacklo_epi16(w0, w1);
    w5 = _mm256_unpacklo_epi16(w2, w3);
    w12 = _mm256_unpacklo_epi16(w8, w9);
    w13 = _mm256_unpacklo_epi16(w10, w11);

    w6 = _mm256_unpacklo_epi32(w4, w5);
    w7 = _mm256_unpackhi_epi32(w4, w5);
    w14 = _mm256_unpacklo_epi32(w12, w13);
    w15 = _mm256_unpackhi_epi32(w12, w13);

    // Store first 4-line result
    d[8] = _mm256_unpacklo_epi64(w6, w14);
    d[9] = _mm256_unpackhi_epi64(w6, w14);
    d[10] = _mm256_unpacklo_epi64(w7, w15);
    d[11] = _mm256_unpackhi_epi64(w7, w15);

    w4 = _mm256_unpackhi_epi16(w0, w1);
    w5 = _mm256_unpackhi_epi16(w2, w3);
    w12 = _mm256_unpackhi_epi16(w8, w9);
    w13 = _mm256_unpackhi_epi16(w10, w11);

    w6 = _mm256_unpacklo_epi32(w4, w5);
    w7 = _mm256_unpackhi_epi32(w4, w5);
    w14 = _mm256_unpacklo_epi32(w12, w13);
    w15 = _mm256_unpackhi_epi32(w12, w13);

    // Store second 4-line result
    d[12] = _mm256_unpacklo_epi64(w6, w14);
    d[13] = _mm256_unpackhi_epi64(w6, w14);
    d[14] = _mm256_unpacklo_epi64(w7, w15);
    d[15] = _mm256_unpackhi_epi64(w7, w15);
}

static INLINE void transpose16x16_sse2(__m128i *x, __m128i *d) {
    __m128i w0, w1, w2, w3, w4, w5, w6, w7, w8, w9;
    __m128i w10, w11, w12, w13, w14, w15;

    w0 = _mm_unpacklo_epi8(x[0], x[1]);
    w1 = _mm_unpacklo_epi8(x[2], x[3]);
    w2 = _mm_unpacklo_epi8(x[4], x[5]);
    w3 = _mm_unpacklo_epi8(x[6], x[7]);

    w8 = _mm_unpacklo_epi8(x[8], x[9]);
    w9 = _mm_unpacklo_epi8(x[10], x[11]);
    w10 = _mm_unpacklo_epi8(x[12], x[13]);
    w11 = _mm_unpacklo_epi8(x[14], x[15]);

    w4 = _mm_unpacklo_epi16(w0, w1);
    w5 = _mm_unpacklo_epi16(w2, w3);
    w12 = _mm_unpacklo_epi16(w8, w9);
    w13 = _mm_unpacklo_epi16(w10, w11);

    w6 = _mm_unpacklo_epi32(w4, w5);
    w7 = _mm_unpackhi_epi32(w4, w5);
    w14 = _mm_unpacklo_epi32(w12, w13);
    w15 = _mm_unpackhi_epi32(w12, w13);

    // Store first 4-line result
    d[0] = _mm_unpacklo_epi64(w6, w14);
    d[1] = _mm_unpackhi_epi64(w6, w14);
    d[2] = _mm_unpacklo_epi64(w7, w15);
    d[3] = _mm_unpackhi_epi64(w7, w15);

    w4 = _mm_unpackhi_epi16(w0, w1);
    w5 = _mm_unpackhi_epi16(w2, w3);
    w12 = _mm_unpackhi_epi16(w8, w9);
    w13 = _mm_unpackhi_epi16(w10, w11);

    w6 = _mm_unpacklo_epi32(w4, w5);
    w7 = _mm_unpackhi_epi32(w4, w5);
    w14 = _mm_unpacklo_epi32(w12, w13);
    w15 = _mm_unpackhi_epi32(w12, w13);

    // Store second 4-line result
    d[4] = _mm_unpacklo_epi64(w6, w14);
    d[5] = _mm_unpackhi_epi64(w6, w14);
    d[6] = _mm_unpacklo_epi64(w7, w15);
    d[7] = _mm_unpackhi_epi64(w7, w15);

    // upper half
    w0 = _mm_unpackhi_epi8(x[0], x[1]);
    w1 = _mm_unpackhi_epi8(x[2], x[3]);
    w2 = _mm_unpackhi_epi8(x[4], x[5]);
    w3 = _mm_unpackhi_epi8(x[6], x[7]);

    w8 = _mm_unpackhi_epi8(x[8], x[9]);
    w9 = _mm_unpackhi_epi8(x[10], x[11]);
    w10 = _mm_unpackhi_epi8(x[12], x[13]);
    w11 = _mm_unpackhi_epi8(x[14], x[15]);

    w4 = _mm_unpacklo_epi16(w0, w1);
    w5 = _mm_unpacklo_epi16(w2, w3);
    w12 = _mm_unpacklo_epi16(w8, w9);
    w13 = _mm_unpacklo_epi16(w10, w11);

    w6 = _mm_unpacklo_epi32(w4, w5);
    w7 = _mm_unpackhi_epi32(w4, w5);
    w14 = _mm_unpacklo_epi32(w12, w13);
    w15 = _mm_unpackhi_epi32(w12, w13);

    // Store first 4-line result
    d[8] = _mm_unpacklo_epi64(w6, w14);
    d[9] = _mm_unpackhi_epi64(w6, w14);
    d[10] = _mm_unpacklo_epi64(w7, w15);
    d[11] = _mm_unpackhi_epi64(w7, w15);

    w4 = _mm_unpackhi_epi16(w0, w1);
    w5 = _mm_unpackhi_epi16(w2, w3);
    w12 = _mm_unpackhi_epi16(w8, w9);
    w13 = _mm_unpackhi_epi16(w10, w11);

    w6 = _mm_unpacklo_epi32(w4, w5);
    w7 = _mm_unpackhi_epi32(w4, w5);
    w14 = _mm_unpacklo_epi32(w12, w13);
    w15 = _mm_unpackhi_epi32(w12, w13);

    // Store second 4-line result
    d[12] = _mm_unpacklo_epi64(w6, w14);
    d[13] = _mm_unpackhi_epi64(w6, w14);
    d[14] = _mm_unpacklo_epi64(w7, w15);
    d[15] = _mm_unpackhi_epi64(w7, w15);
}

static void transpose_TX_8X8(const uint8_t *src, ptrdiff_t pitchSrc,
    uint8_t *dst, ptrdiff_t pitchDst) {
    __m128i r0, r1, r2, r3, r4, r5, r6, r7;
    __m128i d0d1, d2d3, d4d5, d6d7;
    r0 = _mm_loadl_epi64((__m128i *)(src + 0 * pitchSrc));
    r1 = _mm_loadl_epi64((__m128i *)(src + 1 * pitchSrc));
    r2 = _mm_loadl_epi64((__m128i *)(src + 2 * pitchSrc));
    r3 = _mm_loadl_epi64((__m128i *)(src + 3 * pitchSrc));
    r4 = _mm_loadl_epi64((__m128i *)(src + 4 * pitchSrc));
    r5 = _mm_loadl_epi64((__m128i *)(src + 5 * pitchSrc));
    r6 = _mm_loadl_epi64((__m128i *)(src + 6 * pitchSrc));
    r7 = _mm_loadl_epi64((__m128i *)(src + 7 * pitchSrc));

    transpose8x8_sse2(&r0, &r1, &r2, &r3, &r4, &r5, &r6, &r7, &d0d1, &d2d3, &d4d5,
        &d6d7);

    _mm_storel_epi64((__m128i *)(dst + 0 * pitchDst), d0d1);
    _mm_storel_epi64((__m128i *)(dst + 1 * pitchDst), _mm_srli_si128(d0d1, 8));
    _mm_storel_epi64((__m128i *)(dst + 2 * pitchDst), d2d3);
    _mm_storel_epi64((__m128i *)(dst + 3 * pitchDst), _mm_srli_si128(d2d3, 8));
    _mm_storel_epi64((__m128i *)(dst + 4 * pitchDst), d4d5);
    _mm_storel_epi64((__m128i *)(dst + 5 * pitchDst), _mm_srli_si128(d4d5, 8));
    _mm_storel_epi64((__m128i *)(dst + 6 * pitchDst), d6d7);
    _mm_storel_epi64((__m128i *)(dst + 7 * pitchDst), _mm_srli_si128(d6d7, 8));
}

static void transpose(const uint8_t *src, ptrdiff_t pitchSrc, uint8_t *dst,
    ptrdiff_t pitchDst, int32_t width, int32_t height) {
    for (int32_t j = 0; j < height; j += 8)
        for (int32_t i = 0; i < width; i += 8)
            transpose_TX_8X8(src + i * pitchSrc + j, pitchSrc, dst + j * pitchDst + i,
            pitchDst);
}

static void dr_prediction_z3_4x4_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *left, int32_t upsample_left,
    int32_t dy) {
    __m128i dstvec[4], d[4];

    dr_prediction_z1_4xN_internal_avx2(4, dstvec, left, upsample_left, dy);
    transpose4x8_8x4_low_sse2(&dstvec[0], &dstvec[1], &dstvec[2], &dstvec[3],
        &d[0], &d[1], &d[2], &d[3]);

    *(uint32_t *)(dst + stride * 0) = _mm_cvtsi128_si32(d[0]);
    *(uint32_t *)(dst + stride * 1) = _mm_cvtsi128_si32(d[1]);
    *(uint32_t *)(dst + stride * 2) = _mm_cvtsi128_si32(d[2]);
    *(uint32_t *)(dst + stride * 3) = _mm_cvtsi128_si32(d[3]);
    return;
}

static void dr_prediction_z3_8x8_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *left, int32_t upsample_left,
    int32_t dy) {
    __m128i dstvec[8], d[8];

    dr_prediction_z1_8xN_internal_avx2(8, dstvec, left, upsample_left, dy);
    transpose8x8_sse2(&dstvec[0], &dstvec[1], &dstvec[2], &dstvec[3], &dstvec[4],
        &dstvec[5], &dstvec[6], &dstvec[7], &d[0], &d[1], &d[2],
        &d[3]);

    _mm_storel_epi64((__m128i *)(dst + 0 * stride), d[0]);
    _mm_storel_epi64((__m128i *)(dst + 1 * stride), _mm_srli_si128(d[0], 8));
    _mm_storel_epi64((__m128i *)(dst + 2 * stride), d[1]);
    _mm_storel_epi64((__m128i *)(dst + 3 * stride), _mm_srli_si128(d[1], 8));
    _mm_storel_epi64((__m128i *)(dst + 4 * stride), d[2]);
    _mm_storel_epi64((__m128i *)(dst + 5 * stride), _mm_srli_si128(d[2], 8));
    _mm_storel_epi64((__m128i *)(dst + 6 * stride), d[3]);
    _mm_storel_epi64((__m128i *)(dst + 7 * stride), _mm_srli_si128(d[3], 8));
}

static void dr_prediction_z3_4x8_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *left, int32_t upsample_left,
    int32_t dy) {
    __m128i dstvec[4], d[8];

    dr_prediction_z1_8xN_internal_avx2(4, dstvec, left, upsample_left, dy);
    transpose4x8_8x4_sse2(&dstvec[0], &dstvec[1], &dstvec[2], &dstvec[3], &d[0],
        &d[1], &d[2], &d[3], &d[4], &d[5], &d[6], &d[7]);
    for (int32_t i = 0; i < 8; i++) {
        *(uint32_t *)(dst + stride * i) = _mm_cvtsi128_si32(d[i]);
    }
}

static void dr_prediction_z3_8x4_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *left, int32_t upsample_left,
    int32_t dy) {
    __m128i dstvec[8], d[4];

    dr_prediction_z1_4xN_internal_avx2(8, dstvec, left, upsample_left, dy);
    transpose8x8_low_sse2(&dstvec[0], &dstvec[1], &dstvec[2], &dstvec[3],
        &dstvec[4], &dstvec[5], &dstvec[6], &dstvec[7], &d[0],
        &d[1], &d[2], &d[3]);
    _mm_storel_epi64((__m128i *)(dst + 0 * stride), d[0]);
    _mm_storel_epi64((__m128i *)(dst + 1 * stride), d[1]);
    _mm_storel_epi64((__m128i *)(dst + 2 * stride), d[2]);
    _mm_storel_epi64((__m128i *)(dst + 3 * stride), d[3]);
}

static void dr_prediction_z3_8x16_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *left, int32_t upsample_left,
    int32_t dy) {
    __m128i dstvec[8], d[8];

    dr_prediction_z1_16xN_internal_avx2(8, dstvec, left, upsample_left, dy);
    transpose8x16_16x8_sse2(dstvec, dstvec + 1, dstvec + 2, dstvec + 3,
        dstvec + 4, dstvec + 5, dstvec + 6, dstvec + 7, d,
        d + 1, d + 2, d + 3, d + 4, d + 5, d + 6, d + 7);
    for (int32_t i = 0; i < 8; i++) {
        _mm_storel_epi64((__m128i *)(dst + i * stride), d[i]);
        _mm_storel_epi64((__m128i *)(dst + (i + 8) * stride),
            _mm_srli_si128(d[i], 8));
    }
}

static void dr_prediction_z3_16x8_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *left, int32_t upsample_left,
    int32_t dy) {
    __m128i dstvec[16], d[16];

    dr_prediction_z1_8xN_internal_avx2(16, dstvec, left, upsample_left, dy);
    transpose16x8_8x16_sse2(
        &dstvec[0], &dstvec[1], &dstvec[2], &dstvec[3], &dstvec[4], &dstvec[5],
        &dstvec[6], &dstvec[7], &dstvec[8], &dstvec[9], &dstvec[10], &dstvec[11],
        &dstvec[12], &dstvec[13], &dstvec[14], &dstvec[15], &d[0], &d[1], &d[2],
        &d[3], &d[4], &d[5], &d[6], &d[7]);

    for (int32_t i = 0; i < 8; i++) {
        _mm_storeu_si128((__m128i *)(dst + i * stride), d[i]);
    }
}

static void dr_prediction_z3_4x16_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *left, int32_t upsample_left,
    int32_t dy) {
    __m128i dstvec[4], d[16];

    dr_prediction_z1_16xN_internal_avx2(4, dstvec, left, upsample_left, dy);
    transpose4x16_sse2(dstvec, d);
    for (int32_t i = 0; i < 16; i++) {
        *(uint32_t *)(dst + stride * i) = _mm_cvtsi128_si32(d[i]);
    }
}

static void dr_prediction_z3_16x4_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *left, int32_t upsample_left,
    int32_t dy) {
    __m128i dstvec[16], d[8];

    dr_prediction_z1_4xN_internal_avx2(16, dstvec, left, upsample_left, dy);
    for (int32_t i = 4; i < 8; i++) {
        d[i] = _mm_setzero_si128();
    }
    transpose16x8_8x16_sse2(
        &dstvec[0], &dstvec[1], &dstvec[2], &dstvec[3], &dstvec[4], &dstvec[5],
        &dstvec[6], &dstvec[7], &dstvec[8], &dstvec[9], &dstvec[10], &dstvec[11],
        &dstvec[12], &dstvec[13], &dstvec[14], &dstvec[15], &d[0], &d[1], &d[2],
        &d[3], &d[4], &d[5], &d[6], &d[7]);

    for (int32_t i = 0; i < 4; i++) {
        _mm_storeu_si128((__m128i *)(dst + i * stride), d[i]);
    }
}

static void dr_prediction_z3_8x32_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *left, int32_t upsample_left,
    int32_t dy) {
    __m256i dstvec[16], d[16];

    dr_prediction_z1_32xN_internal_avx2(8, dstvec, left, upsample_left, dy);
    for (int32_t i = 8; i < 16; i++) {
        dstvec[i] = _mm256_setzero_si256();
    }
    transpose16x32_avx2(dstvec, d);

    for (int32_t i = 0; i < 16; i++) {
        _mm_storel_epi64((__m128i *)(dst + i * stride),
            _mm256_castsi256_si128(d[i]));
    }
    for (int32_t i = 0; i < 16; i++) {
        _mm_storel_epi64((__m128i *)(dst + (i + 16) * stride),
            _mm256_extracti128_si256(d[i], 1));
    }
}

static void dr_prediction_z3_32x8_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *left, int32_t upsample_left,
    int32_t dy) {
    __m128i dstvec[32], d[16];

    dr_prediction_z1_8xN_internal_avx2(32, dstvec, left, upsample_left, dy);

    transpose16x8_8x16_sse2(
        &dstvec[0], &dstvec[1], &dstvec[2], &dstvec[3], &dstvec[4], &dstvec[5],
        &dstvec[6], &dstvec[7], &dstvec[8], &dstvec[9], &dstvec[10], &dstvec[11],
        &dstvec[12], &dstvec[13], &dstvec[14], &dstvec[15], &d[0], &d[1], &d[2],
        &d[3], &d[4], &d[5], &d[6], &d[7]);
    transpose16x8_8x16_sse2(
        &dstvec[0 + 16], &dstvec[1 + 16], &dstvec[2 + 16], &dstvec[3 + 16],
        &dstvec[4 + 16], &dstvec[5 + 16], &dstvec[6 + 16], &dstvec[7 + 16],
        &dstvec[8 + 16], &dstvec[9 + 16], &dstvec[10 + 16], &dstvec[11 + 16],
        &dstvec[12 + 16], &dstvec[13 + 16], &dstvec[14 + 16], &dstvec[15 + 16],
        &d[0 + 8], &d[1 + 8], &d[2 + 8], &d[3 + 8], &d[4 + 8], &d[5 + 8],
        &d[6 + 8], &d[7 + 8]);

    for (int32_t i = 0; i < 8; i++) {
        _mm_storeu_si128((__m128i *)(dst + i * stride), d[i]);
        _mm_storeu_si128((__m128i *)(dst + i * stride + 16), d[i + 8]);
    }
}

static void dr_prediction_z3_16x16_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *left, int32_t upsample_left,
    int32_t dy) {
    __m128i dstvec[16], d[16];

    dr_prediction_z1_16xN_internal_avx2(16, dstvec, left, upsample_left, dy);
    transpose16x16_sse2(dstvec, d);

    for (int32_t i = 0; i < 16; i++) {
        _mm_storeu_si128((__m128i *)(dst + i * stride), d[i]);
    }
}

static void dr_prediction_z3_32x32_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *left, int32_t upsample_left,
    int32_t dy) {
    __m256i dstvec[32], d[32];

    dr_prediction_z1_32xN_internal_avx2(32, dstvec, left, upsample_left, dy);
    transpose16x32_avx2(dstvec, d);
    transpose16x32_avx2(dstvec + 16, d + 16);
    for (int32_t j = 0; j < 16; j++) {
        _mm_storeu_si128((__m128i *)(dst + j * stride),
            _mm256_castsi256_si128(d[j]));
        _mm_storeu_si128((__m128i *)(dst + j * stride + 16),
            _mm256_castsi256_si128(d[j + 16]));
    }
    for (int32_t j = 0; j < 16; j++) {
        _mm_storeu_si128((__m128i *)(dst + (j + 16) * stride),
            _mm256_extracti128_si256(d[j], 1));
        _mm_storeu_si128((__m128i *)(dst + (j + 16) * stride + 16),
            _mm256_extracti128_si256(d[j + 16], 1));
    }
}

static void dr_prediction_z3_64x64_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *left, int32_t upsample_left,
    int32_t dy) {
    DECLARE_ALIGNED(16, uint8_t, dstT[64 * 64]);
    dr_prediction_z1_64xN_avx2(64, dstT, 64, left, upsample_left, dy);
    transpose(dstT, 64, dst, stride, 64, 64);
}

static void dr_prediction_z3_16x32_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *left, int32_t upsample_left,
    int32_t dy) {
    __m256i dstvec[16], d[16];

    dr_prediction_z1_32xN_internal_avx2(16, dstvec, left, upsample_left, dy);
    transpose16x32_avx2(dstvec, d);
    // store
    for (int32_t j = 0; j < 16; j++) {
        _mm_storeu_si128((__m128i *)(dst + j * stride),
            _mm256_castsi256_si128(d[j]));
        _mm_storeu_si128((__m128i *)(dst + (j + 16) * stride),
            _mm256_extracti128_si256(d[j], 1));
    }
}

static void dr_prediction_z3_32x16_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *left, int32_t upsample_left,
    int32_t dy) {
    __m128i dstvec[32], d[16];

    dr_prediction_z1_16xN_internal_avx2(32, dstvec, left, upsample_left, dy);
    for (int32_t i = 0; i < 32; i += 16) {
        transpose16x16_sse2((dstvec + i), d);
        for (int32_t j = 0; j < 16; j++) {
            _mm_storeu_si128((__m128i *)(dst + j * stride + i), d[j]);
        }
    }
}

static void dr_prediction_z3_32x64_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *left, int32_t upsample_left,
    int32_t dy) {
    uint8_t dstT[64 * 32];
    dr_prediction_z1_64xN_avx2(32, dstT, 64, left, upsample_left, dy);
    transpose(dstT, 64, dst, stride, 32, 64);
}

static void dr_prediction_z3_64x32_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *left, int32_t upsample_left,
    int32_t dy) {
    uint8_t dstT[32 * 64];
    dr_prediction_z1_32xN_avx2(64, dstT, 32, left, upsample_left, dy);
    transpose(dstT, 32, dst, stride, 64, 32);
    return;
}

static void dr_prediction_z3_16x64_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *left, int32_t upsample_left,
    int32_t dy) {
    uint8_t dstT[64 * 16];
    dr_prediction_z1_64xN_avx2(16, dstT, 64, left, upsample_left, dy);
    transpose(dstT, 64, dst, stride, 16, 64);
}

static void dr_prediction_z3_64x16_avx2(uint8_t *dst, ptrdiff_t stride,
    const uint8_t *left, int32_t upsample_left,
    int32_t dy) {
    __m128i dstvec[64], d[16];

    dr_prediction_z1_16xN_internal_avx2(64, dstvec, left, upsample_left, dy);
    for (int32_t i = 0; i < 64; i += 16) {
        transpose16x16_sse2((dstvec + i), d);
        for (int32_t j = 0; j < 16; j++) {
            _mm_storeu_si128((__m128i *)(dst + j * stride + i), d[j]);
        }
    }
}

void av1_dr_prediction_z3_avx2(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh,
    const uint8_t *above, const uint8_t *left,
    int32_t upsample_left, int32_t dx, int32_t dy) {
    (void)above;
    (void)dx;
    assert(dx == 1);
    assert(dy > 0);

    if (bw == bh) {
        switch (bw) {
        case 4:
            dr_prediction_z3_4x4_avx2(dst, stride, left, upsample_left, dy);
            break;
        case 8:
            dr_prediction_z3_8x8_avx2(dst, stride, left, upsample_left, dy);
            break;
        case 16:
            dr_prediction_z3_16x16_avx2(dst, stride, left, upsample_left, dy);
            break;
        case 32:
            dr_prediction_z3_32x32_avx2(dst, stride, left, upsample_left, dy);
            break;
        case 64:
            dr_prediction_z3_64x64_avx2(dst, stride, left, upsample_left, dy);
            break;
        }
    }
    else {
        if (bw < bh) {
            if (bw + bw == bh) {
                switch (bw) {
                case 4:
                    dr_prediction_z3_4x8_avx2(dst, stride, left, upsample_left, dy);
                    break;
                case 8:
                    dr_prediction_z3_8x16_avx2(dst, stride, left, upsample_left, dy);
                    break;
                case 16:
                    dr_prediction_z3_16x32_avx2(dst, stride, left, upsample_left, dy);
                    break;
                case 32:
                    dr_prediction_z3_32x64_avx2(dst, stride, left, upsample_left, dy);
                    break;
                }
            }
            else {
                switch (bw) {
                case 4:
                    dr_prediction_z3_4x16_avx2(dst, stride, left, upsample_left, dy);
                    break;
                case 8:
                    dr_prediction_z3_8x32_avx2(dst, stride, left, upsample_left, dy);
                    break;
                case 16:
                    dr_prediction_z3_16x64_avx2(dst, stride, left, upsample_left, dy);
                    break;
                }
            }
        }
        else {
            if (bh + bh == bw) {
                switch (bh) {
                case 4:
                    dr_prediction_z3_8x4_avx2(dst, stride, left, upsample_left, dy);
                    break;
                case 8:
                    dr_prediction_z3_16x8_avx2(dst, stride, left, upsample_left, dy);
                    break;
                case 16:
                    dr_prediction_z3_32x16_avx2(dst, stride, left, upsample_left, dy);
                    break;
                case 32:
                    dr_prediction_z3_64x32_avx2(dst, stride, left, upsample_left, dy);
                    break;
                }
            }
            else {
                switch (bh) {
                case 4:
                    dr_prediction_z3_16x4_avx2(dst, stride, left, upsample_left, dy);
                    break;
                case 8:
                    dr_prediction_z3_32x8_avx2(dst, stride, left, upsample_left, dy);
                    break;
                case 16:
                    dr_prediction_z3_64x16_avx2(dst, stride, left, upsample_left, dy);
                    break;
                }
            }
        }
    }
    return;
}


uint16_t z2BlendMaskabove[] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF,
    0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF,
    0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF,
    0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF,
    0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF,
    0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF,
    0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF,
    0xFFFF, 0xFFFF
};
uint16_t z2BlendMaskleft[] = {
    0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF,
    0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF,
    0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF,
    0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF,
    0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF,
    0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF,
    0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF,
    0xFFFF, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0
};
void predict_64x64_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
    int32_t d, int32_t upsample) {
    uint32_t leftBy32[128];
    int_least32_t leftByDiff[128];
    __m256i a0, a1, diff, a32, a16;
    int32_t c, x, y, base;
    const uint16_t *leftPtr = refPel - 1;
    const int32_t frac_bits = 6 - upsample;

    a16 = _mm256_set1_epi32(16);
    a0 =
        _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
    a1 = _mm256_cvtepu16_epi32(
        _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
    diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
    _mm256_storeu_si256((__m256i *)(leftBy32 + 64), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 64), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 8)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 9)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 72), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 72), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 16)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 17)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 80), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 80), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 24)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 25)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 88), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 88), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 32)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 33)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 96), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 96), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 40)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 41)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 104), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 104), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 48)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 49)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 112), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 112), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 56)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 57)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 120), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 120), diff);

    for (c = 0; c < 64; ++c) {
        x = c + 1;
        y = -x * d;
        base = y >> frac_bits;

        __m256i inc = _mm256_set1_epi32(y);
        __m256i shift =
            _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
            _mm256_set1_epi32(0x3f)),
            1);

        __m256i a, b, res1, res2, res3, res4, resLo, resHi;
        base = base + 1;
        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 64));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 64));
        b = _mm256_mullo_epi32(b, shift);
        res1 = _mm256_add_epi32(a, b);
        res1 = _mm256_srli_epi32(res1, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 72));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 72));
        b = _mm256_mullo_epi32(b, shift);
        res2 = _mm256_add_epi32(a, b);
        res2 = _mm256_srli_epi32(res2, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 80));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 80));
        b = _mm256_mullo_epi32(b, shift);
        res3 = _mm256_add_epi32(a, b);
        res3 = _mm256_srli_epi32(res3, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 88));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 88));
        b = _mm256_mullo_epi32(b, shift);
        res4 = _mm256_add_epi32(a, b);
        res4 = _mm256_srli_epi32(res4, 5);

        resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
            PERM4x64(0, 2, 1, 3));
        resHi = _mm256_permute4x64_epi64(_mm256_packus_epi32(res3, res4),
            PERM4x64(0, 2, 1, 3));
        _mm256_storeu_si256((__m256i *)(dst + c * pitch), resLo);
        _mm256_storeu_si256((__m256i *)(dst + 16 + c * pitch), resHi);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 96));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 96));
        b = _mm256_mullo_epi32(b, shift);
        res1 = _mm256_add_epi32(a, b);
        res1 = _mm256_srli_epi32(res1, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 104));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 104));
        b = _mm256_mullo_epi32(b, shift);
        res2 = _mm256_add_epi32(a, b);
        res2 = _mm256_srli_epi32(res2, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 112));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 112));
        b = _mm256_mullo_epi32(b, shift);
        res3 = _mm256_add_epi32(a, b);
        res3 = _mm256_srli_epi32(res3, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 120));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 120));
        b = _mm256_mullo_epi32(b, shift);
        res4 = _mm256_add_epi32(a, b);
        res4 = _mm256_srli_epi32(res4, 5);

        resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
            PERM4x64(0, 2, 1, 3));
        resHi = _mm256_permute4x64_epi64(_mm256_packus_epi32(res3, res4),
            PERM4x64(0, 2, 1, 3));
        _mm256_storeu_si256((__m256i *)(dst + 32 + c * pitch), resLo);
        _mm256_storeu_si256((__m256i *)(dst + 48 + c * pitch), resHi);
    }
}
void predict_64x32_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
    int32_t d, int32_t upsample) {
    uint32_t leftBy32[128];
    int_least32_t leftByDiff[128];
    __m256i a0, a1, diff, a32, a16;
    int32_t c, x, y, base;
    const uint16_t *leftPtr = refPel - 1;
    const int32_t frac_bits = 6 - upsample;

    a16 = _mm256_set1_epi32(16);
    a0 =
        _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
    a1 = _mm256_cvtepu16_epi32(
        _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
    diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
    _mm256_storeu_si256((__m256i *)(leftBy32 + 64), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 64), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 8)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 9)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 72), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 72), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 16)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 17)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 80), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 80), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 24)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 25)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 88), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 88), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 32)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 33)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 96), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 96), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 40)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 41)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 104), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 104), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 48)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 49)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 112), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 112), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 56)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 57)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 120), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 120), diff);

    for (c = 0; c < 32; ++c) {
        x = c + 1;
        y = -x * d;
        base = y >> frac_bits;

        __m256i inc = _mm256_set1_epi32(y);
        __m256i shift =
            _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
            _mm256_set1_epi32(0x3f)),
            1);

        __m256i a, b, res1, res2, res3, res4, resLo, resHi;
        base = base + 1;
        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 64));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 64));
        b = _mm256_mullo_epi32(b, shift);
        res1 = _mm256_add_epi32(a, b);
        res1 = _mm256_srli_epi32(res1, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 72));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 72));
        b = _mm256_mullo_epi32(b, shift);
        res2 = _mm256_add_epi32(a, b);
        res2 = _mm256_srli_epi32(res2, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 80));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 80));
        b = _mm256_mullo_epi32(b, shift);
        res3 = _mm256_add_epi32(a, b);
        res3 = _mm256_srli_epi32(res3, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 88));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 88));
        b = _mm256_mullo_epi32(b, shift);
        res4 = _mm256_add_epi32(a, b);
        res4 = _mm256_srli_epi32(res4, 5);

        resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
            PERM4x64(0, 2, 1, 3));
        resHi = _mm256_permute4x64_epi64(_mm256_packus_epi32(res3, res4),
            PERM4x64(0, 2, 1, 3));
        _mm256_storeu_si256((__m256i *)(dst + c * pitch), resLo);
        _mm256_storeu_si256((__m256i *)(dst + 16 + c * pitch), resHi);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 96));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 96));
        b = _mm256_mullo_epi32(b, shift);
        res1 = _mm256_add_epi32(a, b);
        res1 = _mm256_srli_epi32(res1, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 104));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 104));
        b = _mm256_mullo_epi32(b, shift);
        res2 = _mm256_add_epi32(a, b);
        res2 = _mm256_srli_epi32(res2, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 112));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 112));
        b = _mm256_mullo_epi32(b, shift);
        res3 = _mm256_add_epi32(a, b);
        res3 = _mm256_srli_epi32(res3, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 120));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 120));
        b = _mm256_mullo_epi32(b, shift);
        res4 = _mm256_add_epi32(a, b);
        res4 = _mm256_srli_epi32(res4, 5);

        resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
            PERM4x64(0, 2, 1, 3));
        resHi = _mm256_permute4x64_epi64(_mm256_packus_epi32(res3, res4),
            PERM4x64(0, 2, 1, 3));
        _mm256_storeu_si256((__m256i *)(dst + 32 + c * pitch), resLo);
        _mm256_storeu_si256((__m256i *)(dst + 48 + c * pitch), resHi);
    }
}
void predict_64x16_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
    int32_t d, int32_t upsample) {
    uint32_t leftBy32[128];
    int_least32_t leftByDiff[128];
    __m256i a0, a1, diff, a32, a16;
    int32_t c, x, y, base;
    const uint16_t *leftPtr = refPel - 1;
    const int32_t frac_bits = 6 - upsample;

    a16 = _mm256_set1_epi32(16);
    a0 =
        _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
    a1 = _mm256_cvtepu16_epi32(
        _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
    diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
    _mm256_storeu_si256((__m256i *)(leftBy32 + 64), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 64), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 8)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 9)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 72), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 72), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 16)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 17)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 80), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 80), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 24)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 25)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 88), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 88), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 32)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 33)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 96), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 96), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 40)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 41)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 104), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 104), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 48)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 49)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 112), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 112), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 56)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 57)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 120), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 120), diff);

    for (c = 0; c < 16; ++c) {
        x = c + 1;
        y = -x * d;
        base = y >> frac_bits;

        __m256i inc = _mm256_set1_epi32(y);
        __m256i shift =
            _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
            _mm256_set1_epi32(0x3f)),
            1);

        __m256i a, b, res1, res2, res3, res4, resLo, resHi;
        base = base + 1;
        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 64));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 64));
        b = _mm256_mullo_epi32(b, shift);
        res1 = _mm256_add_epi32(a, b);
        res1 = _mm256_srli_epi32(res1, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 72));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 72));
        b = _mm256_mullo_epi32(b, shift);
        res2 = _mm256_add_epi32(a, b);
        res2 = _mm256_srli_epi32(res2, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 80));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 80));
        b = _mm256_mullo_epi32(b, shift);
        res3 = _mm256_add_epi32(a, b);
        res3 = _mm256_srli_epi32(res3, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 88));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 88));
        b = _mm256_mullo_epi32(b, shift);
        res4 = _mm256_add_epi32(a, b);
        res4 = _mm256_srli_epi32(res4, 5);

        resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
            PERM4x64(0, 2, 1, 3));
        resHi = _mm256_permute4x64_epi64(_mm256_packus_epi32(res3, res4),
            PERM4x64(0, 2, 1, 3));
        _mm256_storeu_si256((__m256i *)(dst + c * pitch), resLo);
        _mm256_storeu_si256((__m256i *)(dst + 16 + c * pitch), resHi);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 96));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 96));
        b = _mm256_mullo_epi32(b, shift);
        res1 = _mm256_add_epi32(a, b);
        res1 = _mm256_srli_epi32(res1, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 104));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 104));
        b = _mm256_mullo_epi32(b, shift);
        res2 = _mm256_add_epi32(a, b);
        res2 = _mm256_srli_epi32(res2, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 112));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 112));
        b = _mm256_mullo_epi32(b, shift);
        res3 = _mm256_add_epi32(a, b);
        res3 = _mm256_srli_epi32(res3, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 120));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 120));
        b = _mm256_mullo_epi32(b, shift);
        res4 = _mm256_add_epi32(a, b);
        res4 = _mm256_srli_epi32(res4, 5);

        resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
            PERM4x64(0, 2, 1, 3));
        resHi = _mm256_permute4x64_epi64(_mm256_packus_epi32(res3, res4),
            PERM4x64(0, 2, 1, 3));
        _mm256_storeu_si256((__m256i *)(dst + 32 + c * pitch), resLo);
        _mm256_storeu_si256((__m256i *)(dst + 48 + c * pitch), resHi);
    }
}
void predict_32x64_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
    int32_t d, int32_t upsample) {
    uint32_t leftBy32[64];
    int_least32_t leftByDiff[64];
    __m256i a0, a1, diff, a32, a16;
    int32_t c, x, y, base;
    const uint16_t *leftPtr = refPel - 1;
    const int32_t frac_bits = 6 - upsample;

    a16 = _mm256_set1_epi32(16);
    a0 =
        _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
    a1 = _mm256_cvtepu16_epi32(
        _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
    diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
    _mm256_storeu_si256((__m256i *)(leftBy32 + 32), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 32), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 8)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 9)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 40), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 40), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 16)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 17)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 48), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 48), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 24)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 25)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 56), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 56), diff);

    for (c = 0; c < 64; ++c) {
        x = c + 1;
        y = -x * d;
        base = y >> frac_bits;

        __m256i inc = _mm256_set1_epi32(y);
        __m256i shift =
            _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
            _mm256_set1_epi32(0x3f)),
            1);

        __m256i a, b, res1, res2, res3, res4, resLo, resHi;
        base = base + 1;
        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 32));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 32));
        b = _mm256_mullo_epi32(b, shift);
        res1 = _mm256_add_epi32(a, b);
        res1 = _mm256_srli_epi32(res1, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 40));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 40));
        b = _mm256_mullo_epi32(b, shift);
        res2 = _mm256_add_epi32(a, b);
        res2 = _mm256_srli_epi32(res2, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 48));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 48));
        b = _mm256_mullo_epi32(b, shift);
        res3 = _mm256_add_epi32(a, b);
        res3 = _mm256_srli_epi32(res3, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 56));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 56));
        b = _mm256_mullo_epi32(b, shift);
        res4 = _mm256_add_epi32(a, b);
        res4 = _mm256_srli_epi32(res4, 5);

        resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
            PERM4x64(0, 2, 1, 3));
        resHi = _mm256_permute4x64_epi64(_mm256_packus_epi32(res3, res4),
            PERM4x64(0, 2, 1, 3));
        _mm256_storeu_si256((__m256i *)(dst + c * pitch), resLo);
        _mm256_storeu_si256((__m256i *)(dst + 16 + c * pitch), resHi);
    }
}
void predict_32x32_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
    int32_t d, int32_t upsample) {
    uint32_t leftBy32[64];
    int_least32_t leftByDiff[64];
    __m256i a0, a1, diff, a32, a16;
    int32_t c, x, y, base;
    const uint16_t *leftPtr = refPel - 1;
    const int32_t frac_bits = 6 - upsample;

    a16 = _mm256_set1_epi32(16);
    a0 =
        _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
    a1 = _mm256_cvtepu16_epi32(
        _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
    diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
    _mm256_storeu_si256((__m256i *)(leftBy32 + 32), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 32), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 8)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 9)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 40), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 40), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 16)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 17)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 48), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 48), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 24)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 25)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 56), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 56), diff);

    for (c = 0; c < 32; ++c) {
        x = c + 1;
        y = -x * d;
        base = y >> frac_bits;

        __m256i inc = _mm256_set1_epi32(y);
        __m256i shift =
            _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
            _mm256_set1_epi32(0x3f)),
            1);

        __m256i a, b, res1, res2, res3, res4, resLo, resHi;
        base = base + 1;
        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 32));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 32));
        b = _mm256_mullo_epi32(b, shift);
        res1 = _mm256_add_epi32(a, b);
        res1 = _mm256_srli_epi32(res1, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 40));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 40));
        b = _mm256_mullo_epi32(b, shift);
        res2 = _mm256_add_epi32(a, b);
        res2 = _mm256_srli_epi32(res2, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 48));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 48));
        b = _mm256_mullo_epi32(b, shift);
        res3 = _mm256_add_epi32(a, b);
        res3 = _mm256_srli_epi32(res3, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 56));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 56));
        b = _mm256_mullo_epi32(b, shift);
        res4 = _mm256_add_epi32(a, b);
        res4 = _mm256_srli_epi32(res4, 5);

        resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
            PERM4x64(0, 2, 1, 3));
        resHi = _mm256_permute4x64_epi64(_mm256_packus_epi32(res3, res4),
            PERM4x64(0, 2, 1, 3));
        _mm256_storeu_si256((__m256i *)(dst + c * pitch), resLo);
        _mm256_storeu_si256((__m256i *)(dst + 16 + c * pitch), resHi);
    }
}
void predict_32x16_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
    int32_t d, int32_t upsample) {
    uint32_t leftBy32[64];
    int_least32_t leftByDiff[64];
    __m256i a0, a1, diff, a32, a16;
    int32_t c, x, y, base;
    const uint16_t *leftPtr = refPel - 1;
    const int32_t frac_bits = 6 - upsample;

    a16 = _mm256_set1_epi32(16);
    a0 =
        _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
    a1 = _mm256_cvtepu16_epi32(
        _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
    diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
    _mm256_storeu_si256((__m256i *)(leftBy32 + 32), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 32), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 8)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 9)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 40), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 40), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 16)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 17)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 48), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 48), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 24)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 25)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 56), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 56), diff);

    for (c = 0; c < 16; ++c) {
        x = c + 1;
        y = -x * d;
        base = y >> frac_bits;

        __m256i inc = _mm256_set1_epi32(y);
        __m256i shift =
            _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
            _mm256_set1_epi32(0x3f)),
            1);

        __m256i a, b, res1, res2, res3, res4, resLo, resHi;
        base = base + 1;
        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 32));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 32));
        b = _mm256_mullo_epi32(b, shift);
        res1 = _mm256_add_epi32(a, b);
        res1 = _mm256_srli_epi32(res1, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 40));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 40));
        b = _mm256_mullo_epi32(b, shift);
        res2 = _mm256_add_epi32(a, b);
        res2 = _mm256_srli_epi32(res2, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 48));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 48));
        b = _mm256_mullo_epi32(b, shift);
        res3 = _mm256_add_epi32(a, b);
        res3 = _mm256_srli_epi32(res3, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 56));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 56));
        b = _mm256_mullo_epi32(b, shift);
        res4 = _mm256_add_epi32(a, b);
        res4 = _mm256_srli_epi32(res4, 5);

        resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
            PERM4x64(0, 2, 1, 3));
        resHi = _mm256_permute4x64_epi64(_mm256_packus_epi32(res3, res4),
            PERM4x64(0, 2, 1, 3));
        _mm256_storeu_si256((__m256i *)(dst + c * pitch), resLo);
        _mm256_storeu_si256((__m256i *)(dst + 16 + c * pitch), resHi);
    }
}
void predict_32x8_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
    int32_t d, int32_t upsample) {
    uint32_t leftBy32[64];
    int_least32_t leftByDiff[64];
    __m256i a0, a1, diff, a32, a16;
    int32_t c, x, y, base;
    const uint16_t *leftPtr = refPel - 1;
    const int32_t frac_bits = 6 - upsample;

    a16 = _mm256_set1_epi32(16);
    a0 =
        _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
    a1 = _mm256_cvtepu16_epi32(
        _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
    diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
    _mm256_storeu_si256((__m256i *)(leftBy32 + 32), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 32), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 8)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 9)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 40), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 40), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 16)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 17)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 48), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 48), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 24)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 25)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 56), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 56), diff);

    for (c = 0; c < 8; ++c) {
        x = c + 1;
        y = -x * d;
        base = y >> frac_bits;

        __m256i inc = _mm256_set1_epi32(y);
        __m256i shift =
            _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
            _mm256_set1_epi32(0x3f)),
            1);

        __m256i a, b, res1, res2, res3, res4, resLo, resHi;
        base = base + 1;
        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 32));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 32));
        b = _mm256_mullo_epi32(b, shift);
        res1 = _mm256_add_epi32(a, b);
        res1 = _mm256_srli_epi32(res1, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 40));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 40));
        b = _mm256_mullo_epi32(b, shift);
        res2 = _mm256_add_epi32(a, b);
        res2 = _mm256_srli_epi32(res2, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 48));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 48));
        b = _mm256_mullo_epi32(b, shift);
        res3 = _mm256_add_epi32(a, b);
        res3 = _mm256_srli_epi32(res3, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 56));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 56));
        b = _mm256_mullo_epi32(b, shift);
        res4 = _mm256_add_epi32(a, b);
        res4 = _mm256_srli_epi32(res4, 5);

        resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
            PERM4x64(0, 2, 1, 3));
        resHi = _mm256_permute4x64_epi64(_mm256_packus_epi32(res3, res4),
            PERM4x64(0, 2, 1, 3));
        _mm256_storeu_si256((__m256i *)(dst + c * pitch), resLo);
        _mm256_storeu_si256((__m256i *)(dst + 16 + c * pitch), resHi);
    }
}
void predict_16x16_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
    int32_t d, int32_t upsample) {
    uint32_t leftBy32[32];
    int_least32_t leftByDiff[32];
    __m256i a0, a1, diff, a32, a16;
    int32_t c, x, y, base;
    const uint16_t *leftPtr = refPel - 1;
    const int32_t frac_bits = 6 - upsample;

    a16 = _mm256_set1_epi32(16);
    a0 =
        _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
    a1 = _mm256_cvtepu16_epi32(
        _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
    diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
    _mm256_storeu_si256((__m256i *)(leftBy32 + 16), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 16), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 8)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 9)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 24), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 24), diff);

    for (c = 0; c < 16; ++c) {
        x = c + 1;
        y = -x * d;
        base = y >> frac_bits;

        __m256i inc = _mm256_set1_epi32(y);
        __m256i shift =
            _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
            _mm256_set1_epi32(0x3f)),
            1);

        __m256i a, b, res1, res2, resLo;
        base = base + 1;
        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 16));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 16));
        b = _mm256_mullo_epi32(b, shift);
        res1 = _mm256_add_epi32(a, b);
        res1 = _mm256_srli_epi32(res1, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 24));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 24));
        b = _mm256_mullo_epi32(b, shift);
        res2 = _mm256_add_epi32(a, b);
        res2 = _mm256_srli_epi32(res2, 5);

        resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
            PERM4x64(0, 2, 1, 3));
        _mm256_storeu_si256((__m256i *)(dst + c * pitch), resLo);
    }
}
void predict_16x64_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
    int32_t d, int32_t upsample) {
    uint32_t leftBy32[32];
    int_least32_t leftByDiff[32];
    __m256i a0, a1, diff, a32, a16;
    int32_t c, x, y, base;
    const uint16_t *leftPtr = refPel - 1;
    const int32_t frac_bits = 6 - upsample;

    a16 = _mm256_set1_epi32(16);
    a0 =
        _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
    a1 = _mm256_cvtepu16_epi32(
        _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
    diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
    _mm256_storeu_si256((__m256i *)(leftBy32 + 16), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 16), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 8)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 9)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 24), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 24), diff);

    for (c = 0; c < 64; ++c) {
        x = c + 1;
        y = -x * d;
        base = y >> frac_bits;

        __m256i inc = _mm256_set1_epi32(y);
        __m256i shift =
            _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
            _mm256_set1_epi32(0x3f)),
            1);

        __m256i a, b, res1, res2, resLo;
        base = base + 1;
        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 16));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 16));
        b = _mm256_mullo_epi32(b, shift);
        res1 = _mm256_add_epi32(a, b);
        res1 = _mm256_srli_epi32(res1, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 24));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 24));
        b = _mm256_mullo_epi32(b, shift);
        res2 = _mm256_add_epi32(a, b);
        res2 = _mm256_srli_epi32(res2, 5);

        resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
            PERM4x64(0, 2, 1, 3));
        _mm256_storeu_si256((__m256i *)(dst + c * pitch), resLo);
    }
}
void predict_16x32_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
    int32_t d, int32_t upsample) {
    uint32_t leftBy32[32];
    int_least32_t leftByDiff[32];
    __m256i a0, a1, diff, a32, a16;
    int32_t c, x, y, base;
    const uint16_t *leftPtr = refPel - 1;
    const int32_t frac_bits = 6 - upsample;

    a16 = _mm256_set1_epi32(16);
    a0 =
        _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
    a1 = _mm256_cvtepu16_epi32(
        _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
    diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
    _mm256_storeu_si256((__m256i *)(leftBy32 + 16), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 16), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 8)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 9)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 24), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 24), diff);

    for (c = 0; c < 32; ++c) {
        x = c + 1;
        y = -x * d;
        base = y >> frac_bits;

        __m256i inc = _mm256_set1_epi32(y);
        __m256i shift =
            _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
            _mm256_set1_epi32(0x3f)),
            1);

        __m256i a, b, res1, res2, resLo;
        base = base + 1;
        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 16));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 16));
        b = _mm256_mullo_epi32(b, shift);
        res1 = _mm256_add_epi32(a, b);
        res1 = _mm256_srli_epi32(res1, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 24));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 24));
        b = _mm256_mullo_epi32(b, shift);
        res2 = _mm256_add_epi32(a, b);
        res2 = _mm256_srli_epi32(res2, 5);

        resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
            PERM4x64(0, 2, 1, 3));
        _mm256_storeu_si256((__m256i *)(dst + c * pitch), resLo);
    }
}
void predict_16x8_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
    int32_t d, int32_t upsample) {
    uint32_t leftBy32[32];
    int_least32_t leftByDiff[32];
    __m256i a0, a1, diff, a32, a16;
    int32_t c, x, y, base;
    const uint16_t *leftPtr = refPel - 1;
    const int32_t frac_bits = 6 - upsample;

    a16 = _mm256_set1_epi32(16);
    a0 =
        _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
    a1 = _mm256_cvtepu16_epi32(
        _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
    diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
    _mm256_storeu_si256((__m256i *)(leftBy32 + 16), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 16), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 8)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 9)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 24), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 24), diff);

    for (c = 0; c < 8; ++c) {
        x = c + 1;
        y = -x * d;
        base = y >> frac_bits;

        __m256i inc = _mm256_set1_epi32(y);
        __m256i shift =
            _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
            _mm256_set1_epi32(0x3f)),
            1);

        __m256i a, b, res1, res2, resLo;
        base = base + 1;
        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 16));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 16));
        b = _mm256_mullo_epi32(b, shift);
        res1 = _mm256_add_epi32(a, b);
        res1 = _mm256_srli_epi32(res1, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 24));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 24));
        b = _mm256_mullo_epi32(b, shift);
        res2 = _mm256_add_epi32(a, b);
        res2 = _mm256_srli_epi32(res2, 5);

        resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
            PERM4x64(0, 2, 1, 3));
        _mm256_storeu_si256((__m256i *)(dst + c * pitch), resLo);
    }
}
void predict_16x4_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
    int32_t d, int32_t upsample) {
    uint32_t leftBy32[32];
    int_least32_t leftByDiff[32];
    __m256i a0, a1, diff, a32, a16;
    int32_t c, x, y, base;
    const uint16_t *leftPtr = refPel - 1;
    const int32_t frac_bits = 6 - upsample;

    a16 = _mm256_set1_epi32(16);
    a0 =
        _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
    a1 = _mm256_cvtepu16_epi32(
        _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
    diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
    _mm256_storeu_si256((__m256i *)(leftBy32 + 16), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 16), diff);

    a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 8)));
    a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 9)));
    diff = _mm256_sub_epi32(a1, a0);
    a32 = _mm256_slli_epi32(a0, 5);
    a32 = _mm256_add_epi32(a32, a16);
    _mm256_storeu_si256((__m256i *)(leftBy32 + 24), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 24), diff);

    for (c = 0; c < 4; ++c) {
        x = c + 1;
        y = -x * d;
        base = y >> frac_bits;

        __m256i inc = _mm256_set1_epi32(y);
        __m256i shift =
            _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
            _mm256_set1_epi32(0x3f)),
            1);

        __m256i a, b, res1, res2, resLo;
        base = base + 1;
        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 16));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 16));
        b = _mm256_mullo_epi32(b, shift);
        res1 = _mm256_add_epi32(a, b);
        res1 = _mm256_srli_epi32(res1, 5);

        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 24));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 24));
        b = _mm256_mullo_epi32(b, shift);
        res2 = _mm256_add_epi32(a, b);
        res2 = _mm256_srli_epi32(res2, 5);

        resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
            PERM4x64(0, 2, 1, 3));
        _mm256_storeu_si256((__m256i *)(dst + c * pitch), resLo);
    }
}
void predict_8x32_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
    int32_t d, int32_t upsample) {
    __m256i a0, a1, diff, a32, a16;
    __m256i a, b, res1, res2, resLo;
    int32_t c, x, y, base;
    const uint16_t *leftPtr = refPel - 1;
    const int32_t frac_bits = 6 - upsample;

    res2 = _mm256_setzero_si256();
    uint32_t leftBy32[16];
    int_least32_t leftByDiff[16];

    a16 = _mm256_set1_epi32(16);
    a0 =
        _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
    a1 = _mm256_cvtepu16_epi32(
        _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
    diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
    _mm256_storeu_si256((__m256i *)(leftBy32 + 8), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 8), diff);

    for (c = 0; c < 32; ++c) {
        x = c + 1;
        y = -x * d;
        base = y >> frac_bits;

        __m256i inc = _mm256_set1_epi32(y);
        __m256i shift =
            _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
            _mm256_set1_epi32(0x3f)),
            1);

        base = base + 1;
        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 8));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 8));
        b = _mm256_mullo_epi32(b, shift);
        res1 = _mm256_add_epi32(a, b);
        res1 = _mm256_srli_epi32(res1, 5);

        resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
            PERM4x64(0, 2, 1, 3));
        _mm_storeu_si128((__m128i *)(dst + c * pitch),
            _mm256_castsi256_si128(resLo));
    }
}
void predict_8x16_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
    int32_t d, int32_t upsample) {
    __m256i a0, a1, diff, a32, a16;
    __m256i a, b, res1, res2, resLo;
    int32_t c, x, y, base;
    const uint16_t *leftPtr = refPel - 1;
    const int32_t frac_bits = 6 - upsample;

    res2 = _mm256_setzero_si256();
    uint32_t leftBy32[16];
    int_least32_t leftByDiff[16];

    a16 = _mm256_set1_epi32(16);
    a0 =
        _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
    a1 = _mm256_cvtepu16_epi32(
        _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
    diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
    a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
    a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
    _mm256_storeu_si256((__m256i *)(leftBy32 + 8), a32);
    _mm256_storeu_si256((__m256i *)(leftByDiff + 8), diff);

    for (c = 0; c < 16; ++c) {
        x = c + 1;
        y = -x * d;
        base = y >> frac_bits;

        __m256i inc = _mm256_set1_epi32(y);
        __m256i shift =
            _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
            _mm256_set1_epi32(0x3f)),
            1);

        base = base + 1;
        a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 8));
        b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 8));
        b = _mm256_mullo_epi32(b, shift);
        res1 = _mm256_add_epi32(a, b);
        res1 = _mm256_srli_epi32(res1, 5);

        resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
            PERM4x64(0, 2, 1, 3));
        _mm_storeu_si128((__m128i *)(dst + c * pitch),
            _mm256_castsi256_si128(resLo));
    }
}
void predict_8x8_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
    int32_t d, int32_t upsample) {
    __m256i a0, a1, diff, a32, a16;
    __m256i a, b, res1, res2, resLo;
    int32_t c, x, y, base;
    const uint16_t *leftPtr = refPel - 1;
    const int32_t frac_bits = 6 - upsample;

    if (upsample) {
        uint32_t leftBy32[32];
        int_least32_t leftByDiff[32];
        leftPtr = refPel - 2;
        a16 = _mm256_set1_epi32(16);
        a0 = _mm256_cvtepu16_epi32(
            _mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
        a1 = _mm256_cvtepu16_epi32(
            _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
        diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
        a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
        a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
        _mm256_storeu_si256((__m256i *)(leftBy32 + 16), a32);
        _mm256_storeu_si256((__m256i *)(leftByDiff + 16), diff);

        a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 8)));
        a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 9)));
        diff = _mm256_sub_epi32(a1, a0);
        a32 = _mm256_slli_epi32(a0, 5);
        a32 = _mm256_add_epi32(a32, a16);
        _mm256_storeu_si256((__m256i *)(leftBy32 + 24), a32);
        _mm256_storeu_si256((__m256i *)(leftByDiff + 24), diff);

        for (c = 0; c < 8; ++c) {
            x = c + 1;
            y = -x * d;
            base = y >> frac_bits;

            __m256i inc = _mm256_set1_epi32(y);
            __m256i shift =
                _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
                _mm256_set1_epi32(0x3f)),
                1);

            base = base + 2;
            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 16));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 16));
            b = _mm256_mullo_epi32(b, shift);
            res1 = _mm256_add_epi32(a, b);
            res1 = _mm256_srli_epi32(res1, 5);
            res1 = _mm256_permutevar8x32_epi32(
                res1, _mm256_set_epi32(14, 12, 10, 8, 6, 4, 2, 0));

            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 24));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 24));
            b = _mm256_mullo_epi32(b, shift);
            res2 = _mm256_add_epi32(a, b);
            res2 = _mm256_srli_epi32(res2, 5);
            res2 = _mm256_permutevar8x32_epi32(
                res2, _mm256_set_epi32(14, 12, 10, 8, 6, 4, 2, 0));
            resLo = _mm256_permute2x128_si256(res1, res2, 0x20);
            res1 = _mm256_permute2x128_si256(resLo, resLo, 0x01);
            res1 = _mm256_packus_epi32(resLo, res1);
            _mm_storeu_si128((__m128i *)(dst + c * pitch),
                _mm256_castsi256_si128(res1));
        }
    }
    else {
        res2 = _mm256_setzero_si256();
        uint32_t leftBy32[16];
        int_least32_t leftByDiff[16];

        a16 = _mm256_set1_epi32(16);
        a0 = _mm256_cvtepu16_epi32(
            _mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
        a1 = _mm256_cvtepu16_epi32(
            _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
        diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
        a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
        a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
        _mm256_storeu_si256((__m256i *)(leftBy32 + 8), a32);
        _mm256_storeu_si256((__m256i *)(leftByDiff + 8), diff);

        for (c = 0; c < 8; ++c) {
            x = c + 1;
            y = -x * d;
            base = y >> frac_bits;

            __m256i inc = _mm256_set1_epi32(y);
            __m256i shift =
                _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
                _mm256_set1_epi32(0x3f)),
                1);

            base = base + 1;
            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 8));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 8));
            b = _mm256_mullo_epi32(b, shift);
            res1 = _mm256_add_epi32(a, b);
            res1 = _mm256_srli_epi32(res1, 5);

            resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
                PERM4x64(0, 2, 1, 3));
            _mm_storeu_si128((__m128i *)(dst + c * pitch),
                _mm256_castsi256_si128(resLo));
        }
    }
}
void predict_8x4_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
    int32_t d, int32_t upsample) {
    __m256i a0, a1, diff, a32, a16;
    __m256i a, b, res1, res2, resLo;
    int32_t c, x, y, base;
    const uint16_t *leftPtr = refPel - 1;
    const int32_t frac_bits = 6 - upsample;

    if (upsample) {
        uint32_t leftBy32[32];
        int_least32_t leftByDiff[32];
        leftPtr = refPel - 2;
        a16 = _mm256_set1_epi32(16);
        a0 = _mm256_cvtepu16_epi32(
            _mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
        a1 = _mm256_cvtepu16_epi32(
            _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
        diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
        a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
        a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
        _mm256_storeu_si256((__m256i *)(leftBy32 + 16), a32);
        _mm256_storeu_si256((__m256i *)(leftByDiff + 16), diff);

        a0 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 8)));
        a1 = _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr + 9)));
        diff = _mm256_sub_epi32(a1, a0);
        a32 = _mm256_slli_epi32(a0, 5);
        a32 = _mm256_add_epi32(a32, a16);
        _mm256_storeu_si256((__m256i *)(leftBy32 + 24), a32);
        _mm256_storeu_si256((__m256i *)(leftByDiff + 24), diff);

        for (c = 0; c < 4; ++c) {
            x = c + 1;
            y = -x * d;
            base = y >> frac_bits;

            __m256i inc = _mm256_set1_epi32(y);
            __m256i shift =
                _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
                _mm256_set1_epi32(0x3f)),
                1);

            base = base + 2;
            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 16));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 16));
            b = _mm256_mullo_epi32(b, shift);
            res1 = _mm256_add_epi32(a, b);
            res1 = _mm256_srli_epi32(res1, 5);
            res1 = _mm256_permutevar8x32_epi32(
                res1, _mm256_set_epi32(14, 12, 10, 8, 6, 4, 2, 0));

            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 24));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 24));
            b = _mm256_mullo_epi32(b, shift);
            res2 = _mm256_add_epi32(a, b);
            res2 = _mm256_srli_epi32(res2, 5);
            res2 = _mm256_permutevar8x32_epi32(
                res2, _mm256_set_epi32(14, 12, 10, 8, 6, 4, 2, 0));
            resLo = _mm256_permute2x128_si256(res1, res2, 0x20);
            res1 = _mm256_permute2x128_si256(resLo, resLo, 0x01);
            res1 = _mm256_packus_epi32(resLo, res1);
            _mm_storeu_si128((__m128i *)(dst + c * pitch),
                _mm256_castsi256_si128(res1));
        }
    }
    else {
        res2 = _mm256_setzero_si256();
        uint32_t leftBy32[16];
        int_least32_t leftByDiff[16];

        a16 = _mm256_set1_epi32(16);
        a0 = _mm256_cvtepu16_epi32(
            _mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
        a1 = _mm256_cvtepu16_epi32(
            _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
        diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
        a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
        a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
        _mm256_storeu_si256((__m256i *)(leftBy32 + 8), a32);
        _mm256_storeu_si256((__m256i *)(leftByDiff + 8), diff);

        for (c = 0; c < 4; ++c) {
            x = c + 1;
            y = -x * d;
            base = y >> frac_bits;

            __m256i inc = _mm256_set1_epi32(y);
            __m256i shift =
                _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
                _mm256_set1_epi32(0x3f)),
                1);

            base = base + 1;
            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 8));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 8));
            b = _mm256_mullo_epi32(b, shift);
            res1 = _mm256_add_epi32(a, b);
            res1 = _mm256_srli_epi32(res1, 5);

            resLo = _mm256_permute4x64_epi64(_mm256_packus_epi32(res1, res2),
                PERM4x64(0, 2, 1, 3));
            _mm_storeu_si128((__m128i *)(dst + c * pitch),
                _mm256_castsi256_si128(resLo));
        }
    }
}
void predict_4x4_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
    int32_t d, int32_t upsample) {
    int32_t c, x, y, base;
    const uint16_t *leftPtr = refPel - 1;
    const int32_t frac_bits = 6 - upsample;

    if (upsample) {
        __m256i a0, a1, diff, a32, a16;
        __m256i a, b, res1;
        __m128i res;
        uint32_t leftBy32[16];
        int_least32_t leftByDiff[16];
        leftPtr = refPel - 2;

        a16 = _mm256_set1_epi32(16);
        a0 = _mm256_cvtepu16_epi32(
            _mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
        a1 = _mm256_cvtepu16_epi32(
            _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
        diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
        a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
        a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
        _mm256_storeu_si256((__m256i *)(leftBy32 + 8), a32);
        _mm256_storeu_si256((__m256i *)(leftByDiff + 8), diff);

        for (c = 0; c < 4; ++c) {
            x = c + 1;
            y = -x * d;
            base = y >> frac_bits;

            __m256i inc = _mm256_set1_epi32(y);
            __m256i shift =
                _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
                _mm256_set1_epi32(0x3f)),
                1);

            base = base + 2;
            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 8));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 8));
            b = _mm256_mullo_epi32(b, shift);
            res1 = _mm256_add_epi32(a, b);
            res1 = _mm256_srli_epi32(res1, 5);
            res1 = _mm256_permutevar8x32_epi32(
                res1, _mm256_set_epi32(14, 12, 10, 8, 6, 4, 2, 0));

            res = _mm256_extractf128_si256(res1, 0);
            _mm_storel_epi64((__m128i *)(dst + c * pitch),
                _mm_packus_epi32(res, res));
        }
    }
    else {
        uint32_t leftBy32[8];
        int_least32_t leftByDiff[8];
        __m128i a0, a1, diff, a32, a16;
        __m128i a, b, res1, res2, resLo;
        res2 = _mm_setzero_si128();

        a16 = _mm_set1_epi32(16);
        a0 = _mm_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
        a1 = _mm_cvtepu16_epi32(
            _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
        diff = _mm_sub_epi32(a1, a0);                    // a[x+1] - a[x]
        a32 = _mm_slli_epi32(a0, 5);                     // a[x] * 32
        a32 = _mm_add_epi32(a32, a16);                   // a[x] * 32 + 16
        _mm_storeu_si128((__m128i *)(leftBy32 + 4), a32);
        _mm_storeu_si128((__m128i *)(leftByDiff + 4), diff);

        for (c = 0; c < 4; ++c) {
            x = c + 1;
            y = -x * d;
            base = y >> frac_bits;

            __m128i inc = _mm_set1_epi32(y);
            __m128i shift = _mm_srli_epi32(
                _mm_and_si128(_mm_slli_epi32(inc, upsample), _mm_set1_epi32(0x3f)),
                1);

            base = base + 1;
            a = _mm_loadu_si128((__m128i *)(leftBy32 + base + 4));
            b = _mm_loadu_si128((__m128i *)(leftByDiff + base + 4));
            b = _mm_mullo_epi32(b, shift);
            res1 = _mm_add_epi32(a, b);
            res1 = _mm_srli_epi32(res1, 5);

            resLo = _mm_packus_epi32(res1, res2);
            _mm_storel_epi64((__m128i *)(dst + c * pitch), resLo);
        }
    }
}
void predict_4x8_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
    int32_t d, int32_t upsample) {
    int32_t c, x, y, base;
    const uint16_t *leftPtr = refPel - 1;
    const int32_t frac_bits = 6 - upsample;

    if (upsample) {
        __m256i a0, a1, diff, a32, a16;
        __m256i a, b, res1;
        __m128i res;
        uint32_t leftBy32[16];
        int_least32_t leftByDiff[16];
        leftPtr = refPel - 2;

        a16 = _mm256_set1_epi32(16);
        a0 = _mm256_cvtepu16_epi32(
            _mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
        a1 = _mm256_cvtepu16_epi32(
            _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
        diff = _mm256_sub_epi32(a1, a0);                 // a[x+1] - a[x]
        a32 = _mm256_slli_epi32(a0, 5);                  // a[x] * 32
        a32 = _mm256_add_epi32(a32, a16);                // a[x] * 32 + 16
        _mm256_storeu_si256((__m256i *)(leftBy32 + 8), a32);
        _mm256_storeu_si256((__m256i *)(leftByDiff + 8), diff);

        for (c = 0; c < 8; ++c) {
            x = c + 1;
            y = -x * d;
            base = y >> frac_bits;

            __m256i inc = _mm256_set1_epi32(y);
            __m256i shift =
                _mm256_srli_epi32(_mm256_and_si256(_mm256_slli_epi32(inc, upsample),
                _mm256_set1_epi32(0x3f)),
                1);

            base = base + 2;
            a = _mm256_loadu_si256((__m256i *)(leftBy32 + base + 8));
            b = _mm256_loadu_si256((__m256i *)(leftByDiff + base + 8));
            b = _mm256_mullo_epi32(b, shift);
            res1 = _mm256_add_epi32(a, b);
            res1 = _mm256_srli_epi32(res1, 5);
            res1 = _mm256_permutevar8x32_epi32(
                res1, _mm256_set_epi32(14, 12, 10, 8, 6, 4, 2, 0));

            res = _mm256_extractf128_si256(res1, 0);
            _mm_storel_epi64((__m128i *)(dst + c * pitch),
                _mm_packus_epi32(res, res));
        }
    }
    else {
        uint32_t leftBy32[8];
        int_least32_t leftByDiff[8];
        __m128i a0, a1, diff, a32, a16;
        __m128i a, b, res1, res2, resLo;
        res2 = _mm_setzero_si128();

        a16 = _mm_set1_epi32(16);
        a0 = _mm_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
        a1 = _mm_cvtepu16_epi32(
            _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
        diff = _mm_sub_epi32(a1, a0);                    // a[x+1] - a[x]
        a32 = _mm_slli_epi32(a0, 5);                     // a[x] * 32
        a32 = _mm_add_epi32(a32, a16);                   // a[x] * 32 + 16
        _mm_storeu_si128((__m128i *)(leftBy32 + 4), a32);
        _mm_storeu_si128((__m128i *)(leftByDiff + 4), diff);

        for (c = 0; c < 8; ++c) {
            x = c + 1;
            y = -x * d;
            base = y >> frac_bits;

            __m128i inc = _mm_set1_epi32(y);
            __m128i shift = _mm_srli_epi32(
                _mm_and_si128(_mm_slli_epi32(inc, upsample), _mm_set1_epi32(0x3f)),
                1);

            base = base + 1;
            a = _mm_loadu_si128((__m128i *)(leftBy32 + base + 4));
            b = _mm_loadu_si128((__m128i *)(leftByDiff + base + 4));
            b = _mm_mullo_epi32(b, shift);
            res1 = _mm_add_epi32(a, b);
            res1 = _mm_srli_epi32(res1, 5);

            resLo = _mm_packus_epi32(res1, res2);
            _mm_storel_epi64((__m128i *)(dst + c * pitch), resLo);
        }
    }
}
void predict_4x16_sse4(const uint16_t *refPel, uint16_t *dst, ptrdiff_t pitch,
    int32_t d, int32_t upsample) {
    uint32_t leftBy32[8];
    int_least32_t leftByDiff[8];
    __m128i a0, a1, diff, a32, a16;
    __m128i a, b, res1, res2, resLo;
    int32_t c, x, y, base;
    const uint16_t *leftPtr = refPel - 1;
    const int32_t frac_bits = 6 - upsample;

    res2 = _mm_setzero_si128();

    a16 = _mm_set1_epi32(16);
    a0 = _mm_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(leftPtr)));  // 01234567
    a1 = _mm_cvtepu16_epi32(
        _mm_loadu_si128((__m128i *)(leftPtr + 1)));  // 89abcdef
    diff = _mm_sub_epi32(a1, a0);                    // a[x+1] - a[x]
    a32 = _mm_slli_epi32(a0, 5);                     // a[x] * 32
    a32 = _mm_add_epi32(a32, a16);                   // a[x] * 32 + 16
    _mm_storeu_si128((__m128i *)(leftBy32 + 4), a32);
    _mm_storeu_si128((__m128i *)(leftByDiff + 4), diff);

    for (c = 0; c < 16; ++c) {
        x = c + 1;
        y = -x * d;
        base = y >> frac_bits;

        __m128i inc = _mm_set1_epi32(y);
        __m128i shift = _mm_srli_epi32(
            _mm_and_si128(_mm_slli_epi32(inc, upsample), _mm_set1_epi32(0x3f)), 1);

        base = base + 1;
        a = _mm_loadu_si128((__m128i *)(leftBy32 + base + 4));
        b = _mm_loadu_si128((__m128i *)(leftByDiff + base + 4));
        b = _mm_mullo_epi32(b, shift);
        res1 = _mm_add_epi32(a, b);
        res1 = _mm_srli_epi32(res1, 5);

        resLo = _mm_packus_epi32(res1, res2);
        _mm_storel_epi64((__m128i *)(dst + c * pitch), resLo);
    }
}

static void highbd_transpose_TX_8X8(const uint16_t *src, ptrdiff_t pitchSrc,
    uint16_t *dst, ptrdiff_t pitchDst) {
    __m128i r0, r1, r2, r3, r4, r5, r6, r7, r0_Lo, r1_Lo, r2_Lo, r3_Lo, r4_Lo,
        r5_Lo, r6_Lo;
    r0 = _mm_load_si128(
        (__m128i *)(src + 0 * pitchSrc));  // 07,06,05,04,03,02,01,00
    r1 = _mm_load_si128(
        (__m128i *)(src + 1 * pitchSrc));  // 17,16,15,14,13,12,11,10
    r2 = _mm_load_si128(
        (__m128i *)(src + 2 * pitchSrc));  // 27,26,25,24,23,22,21,20
    r3 = _mm_load_si128(
        (__m128i *)(src + 3 * pitchSrc));  // 37,36,35,34,33,32,31,30
    r4 = _mm_load_si128(
        (__m128i *)(src + 4 * pitchSrc));  // 47,46,45,44,43,42,41,40
    r5 = _mm_load_si128(
        (__m128i *)(src + 5 * pitchSrc));  // 57,56,55,54,53,52,51,50
    r6 = _mm_load_si128(
        (__m128i *)(src + 6 * pitchSrc));  // 67,66,65,64,63,62,61,60
    r7 = _mm_load_si128(
        (__m128i *)(src + 7 * pitchSrc));  // 77,76,75,74,73,72,71,70

    r0_Lo = _mm_unpacklo_epi16(r0, r1);
    r2_Lo = _mm_unpacklo_epi16(r2, r3);
    r4_Lo = _mm_unpacklo_epi16(r4, r5);
    r6_Lo = _mm_unpacklo_epi16(r6, r7);

    r1_Lo = r0_Lo;
    r0_Lo = _mm_unpacklo_epi32(r0_Lo, r2_Lo);
    r1_Lo = _mm_unpackhi_epi32(r1_Lo, r2_Lo);
    r5_Lo = r4_Lo;
    r4_Lo = _mm_unpacklo_epi32(r4_Lo, r6_Lo);
    r5_Lo = _mm_unpackhi_epi32(r5_Lo, r6_Lo);
    r2_Lo = r0_Lo;
    r0_Lo = _mm_unpacklo_epi64(r0_Lo, r4_Lo);  // 64
    r2_Lo = _mm_unpackhi_epi64(r2_Lo, r4_Lo);
    r3_Lo = r1_Lo;
    r1_Lo = _mm_unpacklo_epi64(r1_Lo, r5_Lo);
    r3_Lo = _mm_unpackhi_epi64(r3_Lo, r5_Lo);

    _mm_storeu_si128((__m128i *)(dst + 0 * pitchDst), r0_Lo);
    _mm_storeu_si128((__m128i *)(dst + 1 * pitchDst), r2_Lo);
    _mm_storeu_si128((__m128i *)(dst + 2 * pitchDst), r1_Lo);
    _mm_storeu_si128((__m128i *)(dst + 3 * pitchDst), r3_Lo);

    r0 = _mm_unpackhi_epi16(r0, r1);
    r2 = _mm_unpackhi_epi16(r2, r3);
    r4 = _mm_unpackhi_epi16(r4, r5);
    r6 = _mm_unpackhi_epi16(r6, r7);

    r1 = r0;
    r0 = _mm_unpacklo_epi32(r0, r2);
    r1 = _mm_unpackhi_epi32(r1, r2);
    r5 = r4;
    r4 = _mm_unpacklo_epi32(r4, r6);
    r5 = _mm_unpackhi_epi32(r5, r6);
    r2 = r0;
    r0 = _mm_unpacklo_epi64(r0, r4);
    r2 = _mm_unpackhi_epi64(r2, r4);
    r3 = r1;
    r1 = _mm_unpacklo_epi64(r1, r5);
    r3 = _mm_unpackhi_epi64(r3, r5);

    _mm_storeu_si128((__m128i *)(dst + 4 * pitchDst), r0);
    _mm_storeu_si128((__m128i *)(dst + 5 * pitchDst), r2);
    _mm_storeu_si128((__m128i *)(dst + 6 * pitchDst), r1);
    _mm_storeu_si128((__m128i *)(dst + 7 * pitchDst), r3);
}
static uint8_t HighbdLoadMaskx[8][16] = {
    { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 },
    { 0, 1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 },
    { 0, 1, 0, 1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 },
    { 0, 1, 0, 1, 0, 1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 },
    { 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 2, 3, 4, 5, 6, 7 },
    { 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 2, 3, 4, 5 },
    { 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 2, 3 },
    { 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1 },
};

static uint8_t HighbdEvenOddMaskx4[8][16] = {
    { 0, 1, 4, 5, 8, 9, 12, 13, 2, 3, 6, 7, 10, 11, 14,
    15 },  // 0=0,1, 1=2,3, 2=4,5, 3=6,7, 4=8,9, 5=10,11, 6=12,13, 7=14,15,
    // >7=0,1
    { 0, 1, 2, 3, 6, 7, 10, 11, 14, 15, 4, 5, 8, 9, 12, 13 },
    { 0, 1, 0, 1, 4, 5, 8, 9, 12, 13, 0, 1, 6, 7, 10, 11 },
    { 0, 1, 0, 1, 0, 1, 6, 7, 10, 11, 14, 15, 0, 1, 8, 9 },
    { 0, 1, 0, 1, 0, 1, 0, 1, 8, 9, 12, 13, 0, 1, 0, 1 },
    { 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 10, 11, 14, 15, 0, 1 },
    { 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 12, 13, 0, 1 },
    { 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 14, 15 }
};

static uint16_t HighbdEvenOddMaskx8_2[8][16] = {
    { 0, 2, 4, 6, 8, 10, 12, 14 }, { 2, 2, 4, 6, 8, 10, 12, 14 },
    { 4, 4, 4, 6, 8, 10, 12, 14 }, { 6, 6, 6, 6, 8, 10, 12, 14 },
    { 8, 8, 8, 8, 8, 10, 12, 14 }, { 10, 10, 10, 10, 10, 10, 12, 14 },
    { 12, 12, 12, 12, 12, 12, 12, 14 }, { 14, 14, 14, 14, 14, 14, 14, 14 },
};

static uint16_t HighbdBaseMask[17][16] = {
    {
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
    },
    { 0xffff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xffff, 0xffff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xffff, 0xffff, 0xffff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xffff, 0xffff, 0xffff, 0xffff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0 },
    { 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0, 0, 0, 0, 0, 0, 0,
    0, 0 },
    { 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0, 0, 0, 0,
    0, 0, 0, 0 },
    { 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0,
    0, 0, 0, 0, 0, 0 },
    { 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff,
    0xffff, 0, 0, 0, 0, 0, 0 },
    { 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff,
    0xffff, 0xffff, 0, 0, 0, 0, 0 },
    { 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff,
    0xffff, 0xffff, 0xffff, 0, 0, 0, 0 },
    { 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff,
    0xffff, 0xffff, 0xffff, 0xffff, 0, 0, 0 },
    { 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff,
    0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0, 0 },
    { 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff,
    0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0 },
    { 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff,
    0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff }
};

static void highbd_dr_prediction_z2_Nx4_avx2(
    int32_t N, uint16_t *dst, ptrdiff_t stride, const uint16_t *above,
    const uint16_t *left, int32_t upsample_above, int32_t upsample_left, int32_t dx,
    int32_t dy) {
    const int32_t min_base_x = -(1 << upsample_above);
    const int32_t min_base_y = -(1 << upsample_left);
    const int32_t frac_bits_x = 6 - upsample_above;
    const int32_t frac_bits_y = 6 - upsample_left;

    // a assert(dx > 0);
    // pre-filter above pixels
    // store in temp buffers:
    //   above[x] * 32 + 16
    //   above[x+1] - above[x]
    // final pixels will be caluculated as:
    //   (above[x] * 32 + 16 + (above[x+1] - above[x]) * shift) >> 5
    __m256i a0_x, a1_x, a0_y, a1_y, a32, a16;
    __m256i diff;
    __m128i c3f, min_base_y128;

    a16 = _mm256_set1_epi32(16);
    c3f = _mm_set1_epi32(0x3f);
    min_base_y128 = _mm_set1_epi32(min_base_y);

    for (int32_t r = 0; r < N; r++) {
        __m256i b, res, shift;
        __m128i resx, resy, resxy;
        __m128i a0_x128, a1_x128;
        int32_t y = r + 1;
        int32_t base_x = (-y * dx) >> frac_bits_x;
        int32_t base_shift = 0;
        if (base_x < (min_base_x - 1)) {
            base_shift = (min_base_x - base_x - 1) >> upsample_above;
        }
        int32_t base_min_diff =
            (min_base_x - base_x + upsample_above) >> upsample_above;
        if (base_min_diff > 4) {
            base_min_diff = 4;
        }
        else {
            if (base_min_diff < 0) base_min_diff = 0;
        }

        if (base_shift > 3) {
            resx = _mm_setzero_si128();
        }
        else {
            a0_x128 = _mm_loadu_si128((__m128i *)(above + base_x + base_shift));
            a1_x128 = _mm_loadu_si128((__m128i *)(above + base_x + 1 + base_shift));

            if (upsample_above) {
                a0_x128 = _mm_shuffle_epi8(a0_x128,
                    *(__m128i *)HighbdEvenOddMaskx4[base_shift]);
                a1_x128 = _mm_shuffle_epi8(a1_x128,
                    *(__m128i *)HighbdEvenOddMaskx4[base_shift]);
                shift = _mm256_castsi128_si256(_mm_srli_epi32(
                    _mm_and_si128(
                    _mm_slli_epi32(
                    _mm_setr_epi32(-y * dx, (1 << 6) - y * dx,
                    (2 << 6) - y * dx, (3 << 6) - y * dx),
                    upsample_above),
                    c3f),
                    1));
            }
            else {
                a0_x128 =
                    _mm_shuffle_epi8(a0_x128, *(__m128i *)HighbdLoadMaskx[base_shift]);
                a1_x128 =
                    _mm_shuffle_epi8(a1_x128, *(__m128i *)HighbdLoadMaskx[base_shift]);
                shift = _mm256_castsi128_si256(_mm_srli_epi32(
                    _mm_and_si128(_mm_setr_epi32(-y * dx, (1 << 6) - y * dx,
                    (2 << 6) - y * dx, (3 << 6) - y * dx),
                    c3f),
                    1));
            }
            a0_x = _mm256_cvtepu16_epi32(a0_x128);
            a1_x = _mm256_cvtepu16_epi32(a1_x128);

            diff = _mm256_sub_epi32(a1_x, a0_x);  // a[x+1] - a[x]
            a32 = _mm256_slli_epi32(a0_x, 5);     // a[x] * 32
            a32 = _mm256_add_epi32(a32, a16);     // a[x] * 32 + 16

            b = _mm256_mullo_epi32(diff, shift);
            res = _mm256_add_epi32(a32, b);
            res = _mm256_srli_epi32(res, 5);

            resx = _mm256_castsi256_si128(res);
            resx = _mm_packus_epi32(resx, resx);
        }
        // y calc
        if (base_x < min_base_x) {
            DECLARE_ALIGNED(32, int32_t, base_y_c[4]);
            __m128i r6, c1234, dy128, y_c128, base_y_c128, mask128;
            r6 = _mm_set1_epi32(r << 6);
            dy128 = _mm_set1_epi32(dy);
            c1234 = _mm_setr_epi32(1, 2, 3, 4);
            y_c128 = _mm_sub_epi32(r6, _mm_mullo_epi32(c1234, dy128));
            base_y_c128 = _mm_srai_epi32(y_c128, frac_bits_y);
            mask128 = _mm_cmpgt_epi32(min_base_y128, base_y_c128);
            base_y_c128 = _mm_andnot_si128(mask128, base_y_c128);
            _mm_store_si128((__m128i *)base_y_c, base_y_c128);

            a0_y = _mm256_castsi128_si256(
                _mm_setr_epi32(left[base_y_c[0]], left[base_y_c[1]],
                left[base_y_c[2]], left[base_y_c[3]]));
            a1_y = _mm256_castsi128_si256(
                _mm_setr_epi32(left[base_y_c[0] + 1], left[base_y_c[1] + 1],
                left[base_y_c[2] + 1], left[base_y_c[3] + 1]));

            if (upsample_left) {
                shift = _mm256_castsi128_si256(_mm_srli_epi32(
                    _mm_and_si128(_mm_slli_epi32(y_c128, upsample_left), c3f), 1));
            }
            else {
                shift = _mm256_castsi128_si256(
                    _mm_srli_epi32(_mm_and_si128(y_c128, c3f), 1));
            }
            diff = _mm256_sub_epi32(a1_y, a0_y);  // a[x+1] - a[x]
            a32 = _mm256_slli_epi32(a0_y, 5);     // a[x] * 32
            a32 = _mm256_add_epi32(a32, a16);     // a[x] * 32 + 16

            b = _mm256_mullo_epi32(diff, shift);
            res = _mm256_add_epi32(a32, b);
            res = _mm256_srli_epi32(res, 5);

            resy = _mm256_castsi256_si128(res);
            resy = _mm_packus_epi32(resy, resy);
        }
        else {
            resy = resx;
        }
        resxy =
            _mm_blendv_epi8(resx, resy, *(__m128i *)HighbdBaseMask[base_min_diff]);
        _mm_storel_epi64((__m128i *)(dst), resxy);
        dst += stride;
    }
}

static void highbd_dr_prediction_z2_Nx8_avx2(
    int32_t N, uint16_t *dst, ptrdiff_t stride, const uint16_t *above,
    const uint16_t *left, int32_t upsample_above, int32_t upsample_left, int32_t dx,
    int32_t dy) {
    const int32_t min_base_x = -(1 << upsample_above);
    const int32_t min_base_y = -(1 << upsample_left);
    const int32_t frac_bits_x = 6 - upsample_above;
    const int32_t frac_bits_y = 6 - upsample_left;

    // pre-filter above pixels
    // store in temp buffers:
    //   above[x] * 32 + 16
    //   above[x+1] - above[x]
    // final pixels will be caluculated as:
    //   (above[x] * 32 + 16 + (above[x+1] - above[x]) * shift) >> 5
    __m256i a0_x, a1_x, a0_y, a1_y, a32, a16, c3f, min_base_y256;
    __m256i diff;
    __m128i a0_x128, a1_x128;

    a16 = _mm256_set1_epi32(16);
    c3f = _mm256_set1_epi32(0x3f);
    min_base_y256 = _mm256_set1_epi32(min_base_y);

    for (int32_t r = 0; r < N; r++) {
        __m256i b, res, shift;
        __m128i resx, resy, resxy;
        int32_t y = r + 1;
        int32_t base_x = (-y * dx) >> frac_bits_x;
        int32_t base_shift = 0;
        if (base_x < (min_base_x - 1)) {
            base_shift = (min_base_x - base_x - 1) >> upsample_above;
        }
        int32_t base_min_diff =
            (min_base_x - base_x + upsample_above) >> upsample_above;
        if (base_min_diff > 8) {
            base_min_diff = 8;
        }
        else {
            if (base_min_diff < 0) base_min_diff = 0;
        }

        if (base_shift > 7) {
            resx = _mm_setzero_si128();
        }
        else {
            if (upsample_above) {
                a0_x128 = _mm_setr_epi16(
                    above[base_x + HighbdEvenOddMaskx8_2[base_shift][0]],
                    above[base_x + HighbdEvenOddMaskx8_2[base_shift][1]],
                    above[base_x + HighbdEvenOddMaskx8_2[base_shift][2]],
                    above[base_x + HighbdEvenOddMaskx8_2[base_shift][3]],
                    above[base_x + HighbdEvenOddMaskx8_2[base_shift][4]],
                    above[base_x + HighbdEvenOddMaskx8_2[base_shift][5]],
                    above[base_x + HighbdEvenOddMaskx8_2[base_shift][6]],
                    above[base_x + HighbdEvenOddMaskx8_2[base_shift][7]]);
                a1_x128 = _mm_setr_epi16(
                    above[base_x + 1 + HighbdEvenOddMaskx8_2[base_shift][0]],
                    above[base_x + 1 + HighbdEvenOddMaskx8_2[base_shift][1]],
                    above[base_x + 1 + HighbdEvenOddMaskx8_2[base_shift][2]],
                    above[base_x + 1 + HighbdEvenOddMaskx8_2[base_shift][3]],
                    above[base_x + 1 + HighbdEvenOddMaskx8_2[base_shift][4]],
                    above[base_x + 1 + HighbdEvenOddMaskx8_2[base_shift][5]],
                    above[base_x + 1 + HighbdEvenOddMaskx8_2[base_shift][6]],
                    above[base_x + 1 + HighbdEvenOddMaskx8_2[base_shift][7]]);
                shift = _mm256_srli_epi32(
                    _mm256_and_si256(
                    _mm256_slli_epi32(
                    _mm256_setr_epi32(-y * dx, (1 << 6) - y * dx,
                    (2 << 6) - y * dx, (3 << 6) - y * dx,
                    (4 << 6) - y * dx, (5 << 6) - y * dx,
                    (6 << 6) - y * dx, (7 << 6) - y * dx),
                    upsample_above),
                    c3f),
                    1);
            }
            else {
                a0_x128 = _mm_loadu_si128((__m128i *)(above + base_x + base_shift));
                a1_x128 = _mm_loadu_si128((__m128i *)(above + base_x + 1 + base_shift));
                a0_x128 =
                    _mm_shuffle_epi8(a0_x128, *(__m128i *)HighbdLoadMaskx[base_shift]);
                a1_x128 =
                    _mm_shuffle_epi8(a1_x128, *(__m128i *)HighbdLoadMaskx[base_shift]);

                shift = _mm256_srli_epi32(
                    _mm256_and_si256(
                    _mm256_setr_epi32(-y * dx, (1 << 6) - y * dx, (2 << 6) - y * dx,
                    (3 << 6) - y * dx, (4 << 6) - y * dx,
                    (5 << 6) - y * dx, (6 << 6) - y * dx,
                    (7 << 6) - y * dx),
                    c3f),
                    1);
            }

            a0_x = _mm256_cvtepu16_epi32(a0_x128);
            a1_x = _mm256_cvtepu16_epi32(a1_x128);

            diff = _mm256_sub_epi32(a1_x, a0_x);  // a[x+1] - a[x]
            a32 = _mm256_slli_epi32(a0_x, 5);     // a[x] * 32
            a32 = _mm256_add_epi32(a32, a16);     // a[x] * 32 + 16

            b = _mm256_mullo_epi32(diff, shift);
            res = _mm256_add_epi32(a32, b);
            res = _mm256_srli_epi32(res, 5);

            resx = _mm256_castsi256_si128(_mm256_packus_epi32(
                res, _mm256_castsi128_si256(_mm256_extracti128_si256(res, 1))));
        }
        // y calc
        if (base_x < min_base_x) {
            DECLARE_ALIGNED(32, int32_t, base_y_c[8]);
            __m256i r6, c256, dy256, y_c256, base_y_c256, mask256;
            r6 = _mm256_set1_epi32(r << 6);
            dy256 = _mm256_set1_epi32(dy);
            c256 = _mm256_setr_epi32(1, 2, 3, 4, 5, 6, 7, 8);
            y_c256 = _mm256_sub_epi32(r6, _mm256_mullo_epi32(c256, dy256));
            base_y_c256 = _mm256_srai_epi32(y_c256, frac_bits_y);
            mask256 = _mm256_cmpgt_epi32(min_base_y256, base_y_c256);
            base_y_c256 = _mm256_andnot_si256(mask256, base_y_c256);
            _mm256_store_si256((__m256i *)base_y_c, base_y_c256);

            a0_y = _mm256_cvtepu16_epi32(_mm_setr_epi16(
                left[base_y_c[0]], left[base_y_c[1]], left[base_y_c[2]],
                left[base_y_c[3]], left[base_y_c[4]], left[base_y_c[5]],
                left[base_y_c[6]], left[base_y_c[7]]));
            a1_y = _mm256_cvtepu16_epi32(_mm_setr_epi16(
                left[base_y_c[0] + 1], left[base_y_c[1] + 1], left[base_y_c[2] + 1],
                left[base_y_c[3] + 1], left[base_y_c[4] + 1], left[base_y_c[5] + 1],
                left[base_y_c[6] + 1], left[base_y_c[7] + 1]));

            if (upsample_left) {
                shift = _mm256_srli_epi32(
                    _mm256_and_si256(_mm256_slli_epi32((y_c256), upsample_left), c3f),
                    1);
            }
            else {
                shift = _mm256_srli_epi32(_mm256_and_si256(y_c256, c3f), 1);
            }
            diff = _mm256_sub_epi32(a1_y, a0_y);  // a[x+1] - a[x]
            a32 = _mm256_slli_epi32(a0_y, 5);     // a[x] * 32
            a32 = _mm256_add_epi32(a32, a16);     // a[x] * 32 + 16

            b = _mm256_mullo_epi32(diff, shift);
            res = _mm256_add_epi32(a32, b);
            res = _mm256_srli_epi32(res, 5);

            resy = _mm256_castsi256_si128(_mm256_packus_epi32(
                res, _mm256_castsi128_si256(_mm256_extracti128_si256(res, 1))));
        }
        else {
            resy = resx;
        }
        resxy =
            _mm_blendv_epi8(resx, resy, *(__m128i *)HighbdBaseMask[base_min_diff]);
        _mm_storeu_si128((__m128i *)(dst), resxy);
        dst += stride;
    }
}

static void highbd_dr_prediction_z2_HxW_avx2(
    int32_t H, int32_t W, uint16_t *dst, ptrdiff_t stride, const uint16_t *above,
    const uint16_t *left, int32_t upsample_above, int32_t upsample_left, int32_t dx,
    int32_t dy) {
    // here upsample_above and upsample_left are 0 by design of
    // av1_use_intra_edge_upsample
    const int32_t min_base_x = -1;
    const int32_t min_base_y = -1;
    (void)upsample_above;
    (void)upsample_left;
    const int32_t frac_bits_x = 6;
    const int32_t frac_bits_y = 6;

    // pre-filter above pixels
    // store in temp buffers:
    //   above[x] * 32 + 16
    //   above[x+1] - above[x]
    // final pixels will be caluculated as:
    //   (above[x] * 32 + 16 + (above[x+1] - above[x]) * shift) >> 5
    __m256i a0_x, a1_x, a0_y, a1_y, a32, a0_1_x, a1_1_x, a16;
    __m256i diff, min_base_y256, c3f;
    __m128i a0_x128, a1_x128, a0_1_x128, a1_1_x128;

    a16 = _mm256_set1_epi32(16);
    min_base_y256 = _mm256_set1_epi16(min_base_y);
    c3f = _mm256_set1_epi32(0x3f);
    for (int32_t r = 0; r < H; r++) {
        __m256i b, res, shift;
        __m256i resx[2], resy[2];
        __m256i resxy;
        for (int32_t j = 0; j < W; j += 16) {
            int32_t y = r + 1;
            int32_t base_x = (-y * dx) >> frac_bits_x;
            int32_t base_shift = 0;
            if ((base_x + j) < (min_base_x - 1)) {
                base_shift = (min_base_x - (base_x + j) - 1);
            }
            int32_t base_min_diff = (min_base_x - base_x - j);
            if (base_min_diff > 16) {
                base_min_diff = 16;
            }
            else {
                if (base_min_diff < 0) base_min_diff = 0;
            }

            if (base_shift > 7) {
                resx[0] = _mm256_setzero_si256();
            }
            else {
                a0_x128 = _mm_loadu_si128((__m128i *)(above + base_x + base_shift + j));
                a1_x128 =
                    _mm_loadu_si128((__m128i *)(above + base_x + base_shift + 1 + j));
                a0_x128 =
                    _mm_shuffle_epi8(a0_x128, *(__m128i *)HighbdLoadMaskx[base_shift]);
                a1_x128 =
                    _mm_shuffle_epi8(a1_x128, *(__m128i *)HighbdLoadMaskx[base_shift]);

                a0_x = _mm256_cvtepu16_epi32(a0_x128);
                a1_x = _mm256_cvtepu16_epi32(a1_x128);

                shift = _mm256_srli_epi32(
                    _mm256_and_si256(
                    _mm256_setr_epi32(
                    ((0 + j) << 6) - y * dx, ((1 + j) << 6) - y * dx,
                    ((2 + j) << 6) - y * dx, ((3 + j) << 6) - y * dx,
                    ((4 + j) << 6) - y * dx, ((5 + j) << 6) - y * dx,
                    ((6 + j) << 6) - y * dx, ((7 + j) << 6) - y * dx),
                    _mm256_set1_epi32(0x3f)),
                    1);

                diff = _mm256_sub_epi32(a1_x, a0_x);  // a[x+1] - a[x]
                a32 = _mm256_slli_epi32(a0_x, 5);     // a[x] * 32
                a32 = _mm256_add_epi32(a32, a16);     // a[x] * 32 + 16

                b = _mm256_mullo_epi32(diff, shift);
                res = _mm256_add_epi32(a32, b);
                res = _mm256_srli_epi32(res, 5);

                resx[0] = _mm256_packus_epi32(
                    res, _mm256_castsi128_si256(_mm256_extracti128_si256(res, 1)));
            }
            base_shift = 0;
            if ((base_x + j + 8) < (min_base_x - 1)) {
                base_shift = (min_base_x - (base_x + j + 8) - 1);
            }
            if (base_shift > 7) {
                resx[1] = _mm256_setzero_si256();
            }
            else {
                a0_1_x128 =
                    _mm_loadu_si128((__m128i *)(above + base_x + base_shift + 8 + j));
                a1_1_x128 =
                    _mm_loadu_si128((__m128i *)(above + base_x + base_shift + 9 + j));
                a0_1_x128 = _mm_shuffle_epi8(a0_1_x128,
                    *(__m128i *)HighbdLoadMaskx[base_shift]);
                a1_1_x128 = _mm_shuffle_epi8(a1_1_x128,
                    *(__m128i *)HighbdLoadMaskx[base_shift]);

                a0_1_x = _mm256_cvtepu16_epi32(a0_1_x128);
                a1_1_x = _mm256_cvtepu16_epi32(a1_1_x128);

                shift = _mm256_srli_epi32(
                    _mm256_and_si256(
                    _mm256_setr_epi32(
                    ((8 + j) << 6) - y * dx, ((9 + j) << 6) - y * dx,
                    ((10 + j) << 6) - y * dx, ((11 + j) << 6) - y * dx,
                    ((12 + j) << 6) - y * dx, ((13 + j) << 6) - y * dx,
                    ((14 + j) << 6) - y * dx, ((15 + j) << 6) - y * dx),
                    _mm256_set1_epi32(0x3f)),
                    1);

                diff = _mm256_sub_epi32(a1_1_x, a0_1_x);  // a[x+1] - a[x]
                a32 = _mm256_slli_epi32(a0_1_x, 5);       // a[x] * 32
                a32 = _mm256_add_epi32(a32, a16);         // a[x] * 32 + 16
                b = _mm256_mullo_epi32(diff, shift);

                resx[1] = _mm256_add_epi32(a32, b);
                resx[1] = _mm256_srli_epi32(resx[1], 5);
                resx[1] = _mm256_packus_epi32(
                    resx[1],
                    _mm256_castsi128_si256(_mm256_extracti128_si256(resx[1], 1)));
            }
            resx[0] =
                _mm256_inserti128_si256(resx[0], _mm256_castsi256_si128(resx[1]),
                1);  // 16 16bit values

            // y calc
            if ((base_x < min_base_x)) {
                DECLARE_ALIGNED(32, int32_t, base_y_c[16]);
                __m256i r6, c256, dy256, y_c256, y_c_1_256, base_y_c256, mask256;
                r6 = _mm256_set1_epi32(r << 6);
                dy256 = _mm256_set1_epi32(dy);
                c256 = _mm256_setr_epi32(1 + j, 2 + j, 3 + j, 4 + j, 5 + j, 6 + j,
                    7 + j, 8 + j);
                y_c256 = _mm256_sub_epi32(r6, _mm256_mullo_epi32(c256, dy256));
                base_y_c256 = _mm256_srai_epi32(y_c256, frac_bits_y);
                mask256 = _mm256_cmpgt_epi32(min_base_y256, base_y_c256);
                base_y_c256 = _mm256_andnot_si256(mask256, base_y_c256);
                _mm256_store_si256((__m256i *)base_y_c, base_y_c256);
                c256 = _mm256_setr_epi32(9 + j, 10 + j, 11 + j, 12 + j, 13 + j, 14 + j,
                    15 + j, 16 + j);
                y_c_1_256 = _mm256_sub_epi32(r6, _mm256_mullo_epi32(c256, dy256));
                base_y_c256 = _mm256_srai_epi32(y_c_1_256, frac_bits_y);
                mask256 = _mm256_cmpgt_epi32(min_base_y256, base_y_c256);
                base_y_c256 = _mm256_andnot_si256(mask256, base_y_c256);
                _mm256_store_si256((__m256i *)(base_y_c + 8), base_y_c256);

                a0_y = _mm256_cvtepu16_epi32(_mm_setr_epi16(
                    left[base_y_c[0]], left[base_y_c[1]], left[base_y_c[2]],
                    left[base_y_c[3]], left[base_y_c[4]], left[base_y_c[5]],
                    left[base_y_c[6]], left[base_y_c[7]]));
                a1_y = _mm256_cvtepu16_epi32(_mm_setr_epi16(
                    left[base_y_c[0] + 1], left[base_y_c[1] + 1], left[base_y_c[2] + 1],
                    left[base_y_c[3] + 1], left[base_y_c[4] + 1], left[base_y_c[5] + 1],
                    left[base_y_c[6] + 1], left[base_y_c[7] + 1]));

                shift = _mm256_srli_epi32(_mm256_and_si256(y_c256, c3f), 1);

                diff = _mm256_sub_epi32(a1_y, a0_y);  // a[x+1] - a[x]
                a32 = _mm256_slli_epi32(a0_y, 5);     // a[x] * 32
                a32 = _mm256_add_epi32(a32, a16);     // a[x] * 32 + 16

                b = _mm256_mullo_epi32(diff, shift);
                res = _mm256_add_epi32(a32, b);
                res = _mm256_srli_epi32(res, 5);

                resy[0] = _mm256_packus_epi32(
                    res, _mm256_castsi128_si256(_mm256_extracti128_si256(res, 1)));

                a0_y = _mm256_cvtepu16_epi32(_mm_setr_epi16(
                    left[base_y_c[8]], left[base_y_c[9]], left[base_y_c[10]],
                    left[base_y_c[11]], left[base_y_c[12]], left[base_y_c[13]],
                    left[base_y_c[14]], left[base_y_c[15]]));
                a1_y = _mm256_cvtepu16_epi32(
                    _mm_setr_epi16(left[base_y_c[8] + 1], left[base_y_c[9] + 1],
                    left[base_y_c[10] + 1], left[base_y_c[11] + 1],
                    left[base_y_c[12] + 1], left[base_y_c[13] + 1],
                    left[base_y_c[14] + 1], left[base_y_c[15] + 1]));
                shift = _mm256_srli_epi32(_mm256_and_si256(y_c_1_256, c3f), 1);

                diff = _mm256_sub_epi32(a1_y, a0_y);  // a[x+1] - a[x]
                a32 = _mm256_slli_epi32(a0_y, 5);     // a[x] * 32
                a32 = _mm256_add_epi32(a32, a16);     // a[x] * 32 + 16

                b = _mm256_mullo_epi32(diff, shift);
                res = _mm256_add_epi32(a32, b);
                res = _mm256_srli_epi32(res, 5);

                resy[1] = _mm256_packus_epi32(
                    res, _mm256_castsi128_si256(_mm256_extracti128_si256(res, 1)));

                resy[0] =
                    _mm256_inserti128_si256(resy[0], _mm256_castsi256_si128(resy[1]),
                    1);  // 16 16bit values
            }
            else {
                resy[0] = resx[0];
            }
            resxy = _mm256_blendv_epi8(resx[0], resy[0],
                *(__m256i *)HighbdBaseMask[base_min_diff]);
            _mm256_storeu_si256((__m256i *)(dst + j), resxy);
        }  // for j
        dst += stride;
    }
}

// Directional prediction, zone 2: 90 < angle < 180
void av1_highbd_dr_prediction_z2_avx2(uint16_t *dst, ptrdiff_t stride, int32_t bw,
    int32_t bh, const uint16_t *above, const uint16_t *left, int32_t upsample_above,
    int32_t upsample_left, int32_t dx, int32_t dy, int32_t bd) {
    (void)bd;
    assert(dx > 0);
    assert(dy > 0);
    switch (bw) {
    case 4:
        highbd_dr_prediction_z2_Nx4_avx2(bh, dst, stride, above, left,
            upsample_above, upsample_left, dx, dy);
        break;
    case 8:
        highbd_dr_prediction_z2_Nx8_avx2(bh, dst, stride, above, left,
            upsample_above, upsample_left, dx, dy);
        break;
    default:
        highbd_dr_prediction_z2_HxW_avx2(bh, bw, dst, stride, above, left,
            upsample_above, upsample_left, dx, dy);
        break;
    }
    return;
}
static void highbd_transpose(const uint16_t *src, ptrdiff_t pitchSrc,
    uint16_t *dst, ptrdiff_t pitchDst, int32_t width,
    int32_t height) {
    for (int32_t j = 0; j < height; j += 8)
        for (int32_t i = 0; i < width; i += 8)
            highbd_transpose_TX_8X8(src + i * pitchSrc + j, pitchSrc,
            dst + j * pitchDst + i, pitchDst);
}

static void highbd_dr_prediction_z3_4x4_avx2(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *left,
    int32_t upsample_left, int32_t dy) {
    __m128i dstvec[4], d[4];

    highbd_dr_prediction_z1_4xN_internal_avx2(4, dstvec, left, upsample_left, dy);
    highbd_transpose4x8_8x4_low_sse2(&dstvec[0], &dstvec[1], &dstvec[2],
        &dstvec[3], &d[0], &d[1], &d[2], &d[3]);
    _mm_storel_epi64((__m128i *)(dst + 0 * stride), d[0]);
    _mm_storel_epi64((__m128i *)(dst + 1 * stride), d[1]);
    _mm_storel_epi64((__m128i *)(dst + 2 * stride), d[2]);
    _mm_storel_epi64((__m128i *)(dst + 3 * stride), d[3]);
    return;
}

static void highbd_dr_prediction_z3_8x8_avx2(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *left,
    int32_t upsample_left, int32_t dy) {
    __m128i dstvec[8], d[8];

    highbd_dr_prediction_z1_8xN_internal_avx2(8, dstvec, left, upsample_left, dy);
    highbd_transpose8x8_sse2(&dstvec[0], &dstvec[1], &dstvec[2], &dstvec[3],
        &dstvec[4], &dstvec[5], &dstvec[6], &dstvec[7],
        &d[0], &d[1], &d[2], &d[3], &d[4], &d[5], &d[6],
        &d[7]);
    for (int32_t i = 0; i < 8; i++) {
        _mm_storeu_si128((__m128i *)(dst + i * stride), d[i]);
    }
}

static void highbd_dr_prediction_z3_4x8_avx2(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *left,
    int32_t upsample_left, int32_t dy) {
    __m128i dstvec[4], d[8];

    highbd_dr_prediction_z1_8xN_internal_avx2(4, dstvec, left, upsample_left, dy);
    highbd_transpose4x8_8x4_sse2(&dstvec[0], &dstvec[1], &dstvec[2], &dstvec[3],
        &d[0], &d[1], &d[2], &d[3], &d[4], &d[5], &d[6],
        &d[7]);
    for (int32_t i = 0; i < 8; i++) {
        _mm_storel_epi64((__m128i *)(dst + i * stride), d[i]);
    }
}

static void highbd_dr_prediction_z3_8x4_avx2(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *left,
    int32_t upsample_left, int32_t dy) {
    __m128i dstvec[8], d[4];

    highbd_dr_prediction_z1_4xN_internal_avx2(8, dstvec, left, upsample_left, dy);
    highbd_transpose8x8_low_sse2(&dstvec[0], &dstvec[1], &dstvec[2], &dstvec[3],
        &dstvec[4], &dstvec[5], &dstvec[6], &dstvec[7],
        &d[0], &d[1], &d[2], &d[3]);
    _mm_storeu_si128((__m128i *)(dst + 0 * stride), d[0]);
    _mm_storeu_si128((__m128i *)(dst + 1 * stride), d[1]);
    _mm_storeu_si128((__m128i *)(dst + 2 * stride), d[2]);
    _mm_storeu_si128((__m128i *)(dst + 3 * stride), d[3]);
}

static void highbd_dr_prediction_z3_8x16_avx2(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *left,
    int32_t upsample_left, int32_t dy) {
    __m256i dstvec[8], d[16];

    highbd_dr_prediction_z1_16xN_internal_avx2(8, dstvec, left, upsample_left,
        dy);
    highbd_transpose8x16_16x8_avx2(dstvec, d);
    for (int32_t i = 0; i < 8; i++) {
        _mm_storeu_si128((__m128i *)(dst + i * stride),
            _mm256_castsi256_si128(d[i]));
    }
    for (int32_t i = 8; i < 16; i++) {
        _mm_storeu_si128((__m128i *)(dst + i * stride),
            _mm256_extracti128_si256(d[i - 8], 1));
    }
}

static void highbd_dr_prediction_z3_16x8_avx2(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *left,
    int32_t upsample_left, int32_t dy) {
    __m128i dstvec[16], d[16];

    highbd_dr_prediction_z1_8xN_internal_avx2(16, dstvec, left, upsample_left,
        dy);
    for (int32_t i = 0; i < 16; i += 8) {
        highbd_transpose8x8_sse2(&dstvec[0 + i], &dstvec[1 + i], &dstvec[2 + i],
            &dstvec[3 + i], &dstvec[4 + i], &dstvec[5 + i],
            &dstvec[6 + i], &dstvec[7 + i], &d[0 + i],
            &d[1 + i], &d[2 + i], &d[3 + i], &d[4 + i],
            &d[5 + i], &d[6 + i], &d[7 + i]);
    }
    for (int32_t i = 0; i < 8; i++) {
        _mm_storeu_si128((__m128i *)(dst + i * stride), d[i]);
        _mm_storeu_si128((__m128i *)(dst + i * stride + 8), d[i + 8]);
    }
}

static void highbd_dr_prediction_z3_4x16_avx2(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *left,
    int32_t upsample_left, int32_t dy) {
    __m256i dstvec[4], d[4], d1;

    highbd_dr_prediction_z1_16xN_internal_avx2(4, dstvec, left, upsample_left,
        dy);
    highbd_transpose4x16_avx2(dstvec, d);
    for (int32_t i = 0; i < 4; i++) {
        _mm_storel_epi64((__m128i *)(dst + i * stride),
            _mm256_castsi256_si128(d[i]));
        d1 = _mm256_srli_si256(d[i], 8);
        _mm_storel_epi64((__m128i *)(dst + (i + 4) * stride),
            _mm256_castsi256_si128(d1));
        _mm_storel_epi64((__m128i *)(dst + (i + 8) * stride),
            _mm256_extracti128_si256(d[i], 1));
        _mm_storel_epi64((__m128i *)(dst + (i + 12) * stride),
            _mm256_extracti128_si256(d1, 1));
    }
}

static void highbd_dr_prediction_z3_16x4_avx2(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *left,
    int32_t upsample_left, int32_t dy) {
    __m128i dstvec[16], d[8];

    highbd_dr_prediction_z1_4xN_internal_avx2(16, dstvec, left, upsample_left,
        dy);
    highbd_transpose16x4_8x8_sse2(dstvec, d);

    _mm_storeu_si128((__m128i *)(dst + 0 * stride), d[0]);
    _mm_storeu_si128((__m128i *)(dst + 0 * stride + 8), d[1]);
    _mm_storeu_si128((__m128i *)(dst + 1 * stride), d[2]);
    _mm_storeu_si128((__m128i *)(dst + 1 * stride + 8), d[3]);
    _mm_storeu_si128((__m128i *)(dst + 2 * stride), d[4]);
    _mm_storeu_si128((__m128i *)(dst + 2 * stride + 8), d[5]);
    _mm_storeu_si128((__m128i *)(dst + 3 * stride), d[6]);
    _mm_storeu_si128((__m128i *)(dst + 3 * stride + 8), d[7]);
}

static void highbd_dr_prediction_z3_8x32_avx2(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *left,
    int32_t upsample_left, int32_t dy) {
    __m256i dstvec[16], d[16];

    highbd_dr_prediction_z1_32xN_internal_avx2(8, dstvec, left, upsample_left,
        dy);
    for (int32_t i = 0; i < 16; i += 8) {
        highbd_transpose8x16_16x8_avx2(dstvec + i, d + i);
    }

    for (int32_t i = 0; i < 8; i++) {
        _mm_storeu_si128((__m128i *)(dst + i * stride),
            _mm256_castsi256_si128(d[i]));
    }
    for (int32_t i = 0; i < 8; i++) {
        _mm_storeu_si128((__m128i *)(dst + (i + 8) * stride),
            _mm256_extracti128_si256(d[i], 1));
    }
    for (int32_t i = 8; i < 16; i++) {
        _mm_storeu_si128((__m128i *)(dst + (i + 8) * stride),
            _mm256_castsi256_si128(d[i]));
    }
    for (int32_t i = 8; i < 16; i++) {
        _mm_storeu_si128((__m128i *)(dst + (i + 16) * stride),
            _mm256_extracti128_si256(d[i], 1));
    }
}

static void highbd_dr_prediction_z3_32x8_avx2(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *left,
    int32_t upsample_left, int32_t dy) {
    __m128i dstvec[32], d[32];

    highbd_dr_prediction_z1_8xN_internal_avx2(32, dstvec, left, upsample_left,
        dy);
    for (int32_t i = 0; i < 32; i += 8) {
        highbd_transpose8x8_sse2(&dstvec[0 + i], &dstvec[1 + i], &dstvec[2 + i],
            &dstvec[3 + i], &dstvec[4 + i], &dstvec[5 + i],
            &dstvec[6 + i], &dstvec[7 + i], &d[0 + i],
            &d[1 + i], &d[2 + i], &d[3 + i], &d[4 + i],
            &d[5 + i], &d[6 + i], &d[7 + i]);
    }
    for (int32_t i = 0; i < 8; i++) {
        _mm_storeu_si128((__m128i *)(dst + i * stride), d[i]);
        _mm_storeu_si128((__m128i *)(dst + i * stride + 8), d[i + 8]);
        _mm_storeu_si128((__m128i *)(dst + i * stride + 16), d[i + 16]);
        _mm_storeu_si128((__m128i *)(dst + i * stride + 24), d[i + 24]);
    }
}

static void highbd_dr_prediction_z3_16x16_avx2(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *left,
    int32_t upsample_left, int32_t dy) {
    __m256i dstvec[16], d[16];

    highbd_dr_prediction_z1_16xN_internal_avx2(16, dstvec, left, upsample_left,
        dy);
    highbd_transpose16x16_avx2(dstvec, d);

    for (int32_t i = 0; i < 16; i++) {
        _mm256_storeu_si256((__m256i *)(dst + i * stride), d[i]);
    }
}

static void highbd_dr_prediction_z3_32x32_avx2(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *left,
    int32_t upsample_left, int32_t dy) {
    __m256i dstvec[64], d[16];

    highbd_dr_prediction_z1_32xN_internal_avx2(32, dstvec, left, upsample_left,
        dy);

    highbd_transpose16x16_avx2(dstvec, d);
    for (int32_t j = 0; j < 16; j++) {
        _mm256_storeu_si256((__m256i *)(dst + j * stride), d[j]);
    }
    highbd_transpose16x16_avx2(dstvec + 16, d);
    for (int32_t j = 0; j < 16; j++) {
        _mm256_storeu_si256((__m256i *)(dst + j * stride + 16), d[j]);
    }
    highbd_transpose16x16_avx2(dstvec + 32, d);
    for (int32_t j = 0; j < 16; j++) {
        _mm256_storeu_si256((__m256i *)(dst + (j + 16) * stride), d[j]);
    }
    highbd_transpose16x16_avx2(dstvec + 48, d);
    for (int32_t j = 0; j < 16; j++) {
        _mm256_storeu_si256((__m256i *)(dst + (j + 16) * stride + 16), d[j]);
    }
}

static void highbd_dr_prediction_z3_64x64_avx2(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *left,
    int32_t upsample_left, int32_t dy) {
    DECLARE_ALIGNED(16, uint16_t, dstT[64 * 64]);
    highbd_dr_prediction_z1_64xN_avx2(64, dstT, 64, left, upsample_left, dy);
    highbd_transpose(dstT, 64, dst, stride, 64, 64);
}

static void highbd_dr_prediction_z3_16x32_avx2(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *left,
    int32_t upsample_left, int32_t dy) {
    __m256i dstvec[32], d[32];

    highbd_dr_prediction_z1_32xN_internal_avx2(16, dstvec, left, upsample_left,
        dy);
    for (int32_t i = 0; i < 32; i += 8) {
        highbd_transpose8x16_16x8_avx2(dstvec + i, d + i);
    }
    // store
    for (int32_t j = 0; j < 32; j += 16) {
        for (int32_t i = 0; i < 8; i++) {
            _mm_storeu_si128((__m128i *)(dst + (i + j) * stride),
                _mm256_castsi256_si128(d[(i + j)]));
        }
        for (int32_t i = 0; i < 8; i++) {
            _mm_storeu_si128((__m128i *)(dst + (i + j) * stride + 8),
                _mm256_castsi256_si128(d[(i + j) + 8]));
        }
        for (int32_t i = 8; i < 16; i++) {
            _mm256_storeu_si256(
                (__m256i *)(dst + (i + j) * stride),
                _mm256_inserti128_si256(
                d[(i + j)], _mm256_extracti128_si256(d[(i + j) - 8], 1), 0));
        }
    }
}

static void highbd_dr_prediction_z3_32x16_avx2(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *left,
    int32_t upsample_left, int32_t dy) {
    __m256i dstvec[32], d[16];

    highbd_dr_prediction_z1_16xN_internal_avx2(32, dstvec, left, upsample_left,
        dy);
    for (int32_t i = 0; i < 32; i += 16) {
        highbd_transpose16x16_avx2((dstvec + i), d);
        for (int32_t j = 0; j < 16; j++) {
            _mm256_storeu_si256((__m256i *)(dst + j * stride + i), d[j]);
        }
    }
}

static void highbd_dr_prediction_z3_32x64_avx2(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *left,
    int32_t upsample_left, int32_t dy) {
    uint16_t dstT[64 * 32];
    highbd_dr_prediction_z1_64xN_avx2(32, dstT, 64, left, upsample_left, dy);
    highbd_transpose(dstT, 64, dst, stride, 32, 64);
}

static void highbd_dr_prediction_z3_64x32_avx2(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *left,
    int32_t upsample_left, int32_t dy) {
    DECLARE_ALIGNED(16, uint16_t, dstT[32 * 64]);
    highbd_dr_prediction_z1_32xN_avx2(64, dstT, 32, left, upsample_left, dy);
    highbd_transpose(dstT, 32, dst, stride, 64, 32);
    return;
}

static void highbd_dr_prediction_z3_16x64_avx2(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *left,
    int32_t upsample_left, int32_t dy) {
    DECLARE_ALIGNED(16, uint16_t, dstT[64 * 16]);
    highbd_dr_prediction_z1_64xN_avx2(16, dstT, 64, left, upsample_left, dy);
    highbd_transpose(dstT, 64, dst, stride, 16, 64);
}

static void highbd_dr_prediction_z3_64x16_avx2(uint16_t *dst, ptrdiff_t stride,
    const uint16_t *left,
    int32_t upsample_left, int32_t dy) {
    __m256i dstvec[64], d[16];

    highbd_dr_prediction_z1_16xN_internal_avx2(64, dstvec, left, upsample_left,
        dy);
    for (int32_t i = 0; i < 64; i += 16) {
        highbd_transpose16x16_avx2((dstvec + i), d);
        for (int32_t j = 0; j < 16; j++) {
            _mm256_storeu_si256((__m256i *)(dst + j * stride + i), d[j]);
        }
    }
}

void av1_highbd_dr_prediction_z3_avx2(uint16_t *dst, ptrdiff_t stride, int32_t bw,
    int32_t bh, const uint16_t *above,
    const uint16_t *left, int32_t upsample_left,
    int32_t dx, int32_t dy, int32_t bd) {
    (void)above;
    (void)dx;
    (void)bd;
    assert(dx == 1);
    assert(dy > 0);
    if (bw == bh) {
        switch (bw) {
        case 4:
            highbd_dr_prediction_z3_4x4_avx2(dst, stride, left, upsample_left, dy);
            break;
        case 8:
            highbd_dr_prediction_z3_8x8_avx2(dst, stride, left, upsample_left, dy);
            break;
        case 16:
            highbd_dr_prediction_z3_16x16_avx2(dst, stride, left, upsample_left,
                dy);
            break;
        case 32:
            highbd_dr_prediction_z3_32x32_avx2(dst, stride, left, upsample_left,
                dy);
            break;
        case 64:
            highbd_dr_prediction_z3_64x64_avx2(dst, stride, left, upsample_left,
                dy);
            break;
        }
    }
    else {
        if (bw < bh) {
            if (bw + bw == bh) {
                switch (bw) {
                case 4:
                    highbd_dr_prediction_z3_4x8_avx2(dst, stride, left, upsample_left,
                        dy);
                    break;
                case 8:
                    highbd_dr_prediction_z3_8x16_avx2(dst, stride, left, upsample_left,
                        dy);
                    break;
                case 16:
                    highbd_dr_prediction_z3_16x32_avx2(dst, stride, left, upsample_left,
                        dy);
                    break;
                case 32:
                    highbd_dr_prediction_z3_32x64_avx2(dst, stride, left, upsample_left,
                        dy);
                    break;
                }
            }
            else {
                switch (bw) {
                case 4:
                    highbd_dr_prediction_z3_4x16_avx2(dst, stride, left, upsample_left,
                        dy);
                    break;
                case 8:
                    highbd_dr_prediction_z3_8x32_avx2(dst, stride, left, upsample_left,
                        dy);
                    break;
                case 16:
                    highbd_dr_prediction_z3_16x64_avx2(dst, stride, left, upsample_left,
                        dy);
                    break;
                }
            }
        }
        else {
            if (bh + bh == bw) {
                switch (bh) {
                case 4:
                    highbd_dr_prediction_z3_8x4_avx2(dst, stride, left, upsample_left,
                        dy);
                    break;
                case 8:
                    highbd_dr_prediction_z3_16x8_avx2(dst, stride, left, upsample_left,
                        dy);
                    break;
                case 16:
                    highbd_dr_prediction_z3_32x16_avx2(dst, stride, left, upsample_left,
                        dy);
                    break;
                case 32:
                    highbd_dr_prediction_z3_64x32_avx2(dst, stride, left, upsample_left,
                        dy);
                    break;
                }
            }
            else {
                switch (bh) {
                case 4:
                    highbd_dr_prediction_z3_16x4_avx2(dst, stride, left, upsample_left,
                        dy);
                    break;
                case 8:
                    highbd_dr_prediction_z3_32x8_avx2(dst, stride, left, upsample_left,
                        dy);
                    break;
                case 16:
                    highbd_dr_prediction_z3_64x16_avx2(dst, stride, left, upsample_left,
                        dy);
                    break;
                }
            }
        }
    }
    return;
}

#endif


