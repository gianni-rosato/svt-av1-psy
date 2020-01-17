/*
 * Copyright (c) 2018, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef AOM_DSP_X86_TRANSPOSE_SSE2_H_
#define AOM_DSP_X86_TRANSPOSE_SSE2_H_

#include <emmintrin.h> // SSE2

// nclude "./aom_config.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __cplusplus
}
#endif

static INLINE __m128i transpose_8bit_4x4(const __m128i *const in) {
    // Unpack 16 bit elements. Goes from:
    // in[0]: 00 01 02 03
    // in[1]: 10 11 12 13
    // in[2]: 20 21 22 23
    // in[3]: 30 31 32 33
    // to:
    // a0:    00 10 01 11  02 12 03 13
    // a1:    20 30 21 31  22 32 23 33
    const __m128i a0 = _mm_unpacklo_epi8(in[0], in[1]);
    const __m128i a1 = _mm_unpacklo_epi8(in[2], in[3]);

    // Unpack 32 bit elements resulting in:
    // 00 10 20 30  01 11 21 31  02 12 22 32  03 13 23 33
    return _mm_unpacklo_epi16(a0, a1);
}

static INLINE void transpose_8bit_8x8(const __m128i *const in, __m128i *const out) {
    // Unpack 8 bit elements. Goes from:
    // in[0]: 00 01 02 03 04 05 06 07
    // in[1]: 10 11 12 13 14 15 16 17
    // in[2]: 20 21 22 23 24 25 26 27
    // in[3]: 30 31 32 33 34 35 36 37
    // in[4]: 40 41 42 43 44 45 46 47
    // in[5]: 50 51 52 53 54 55 56 57
    // in[6]: 60 61 62 63 64 65 66 67
    // in[7]: 70 71 72 73 74 75 76 77
    // to:
    // a0:    00 10 01 11 02 12 03 13  04 14 05 15 06 16 07 17
    // a1:    20 30 21 31 22 32 23 33  24 34 25 35 26 36 27 37
    // a2:    40 50 41 51 42 52 43 53  44 54 45 55 46 56 47 57
    // a3:    60 70 61 71 62 72 63 73  64 74 65 75 66 76 67 77
    const __m128i a0 = _mm_unpacklo_epi8(in[0], in[1]);
    const __m128i a1 = _mm_unpacklo_epi8(in[2], in[3]);
    const __m128i a2 = _mm_unpacklo_epi8(in[4], in[5]);
    const __m128i a3 = _mm_unpacklo_epi8(in[6], in[7]);

    // Unpack 16 bit elements resulting in:
    // b0: 00 10 20 30 01 11 21 31  02 12 22 32 03 13 23 33
    // b1: 40 50 60 70 41 51 61 71  42 52 62 72 43 53 63 73
    // b2: 04 14 24 34 05 15 25 35  06 16 26 36 07 17 27 37
    // b3: 44 54 64 74 45 55 65 75  46 56 66 76 47 57 67 77
    const __m128i b0 = _mm_unpacklo_epi16(a0, a1);
    const __m128i b1 = _mm_unpackhi_epi16(a0, a1);
    const __m128i b2 = _mm_unpacklo_epi16(a2, a3);
    const __m128i b3 = _mm_unpackhi_epi16(a2, a3);

    // Unpack 32 bit elements resulting in:
    // c0: 00 10 20 30 40 50 60 70  01 11 21 31 41 51 61 71
    // c1: 02 12 22 32 42 52 62 72  03 13 23 33 43 53 63 73
    // c2: 04 14 24 34 44 54 64 74  05 15 25 35 45 55 65 75
    // c3: 06 16 26 36 46 56 66 76  07 17 27 37 47 57 67 77
    const __m128i c0 = _mm_unpacklo_epi32(b0, b2);
    const __m128i c1 = _mm_unpackhi_epi32(b0, b2);
    const __m128i c2 = _mm_unpacklo_epi32(b1, b3);
    const __m128i c3 = _mm_unpackhi_epi32(b1, b3);

    // Unpack 64 bit elements resulting in:
    // out[0]: 00 10 20 30 40 50 60 70
    // out[1]: 01 11 21 31 41 51 61 71
    // out[2]: 02 12 22 32 42 52 62 72
    // out[3]: 03 13 23 33 43 53 63 73
    // out[4]: 04 14 24 34 44 54 64 74
    // out[5]: 05 15 25 35 45 55 65 75
    // out[6]: 06 16 26 36 46 56 66 76
    // out[7]: 07 17 27 37 47 57 67 77
    out[0] = c0;
    out[1] = _mm_srli_si128(c0, 8);
    out[2] = c1;
    out[3] = _mm_srli_si128(c1, 8);
    out[4] = c2;
    out[5] = _mm_srli_si128(c2, 8);
    out[6] = c3;
    out[7] = _mm_srli_si128(c3, 8);
}

static INLINE void partial_transpose_8bit_8x8(const __m128i *const in, __m128i *const out) {
    // Unpack 8 bit elements. Goes from:
    // in[0]: 00 01 02 03 04 05 06 07
    // in[1]: 10 11 12 13 14 15 16 17
    // in[2]: 20 21 22 23 24 25 26 27
    // in[3]: 30 31 32 33 34 35 36 37
    // in[4]: 40 41 42 43 44 45 46 47
    // in[5]: 50 51 52 53 54 55 56 57
    // in[6]: 60 61 62 63 64 65 66 67
    // in[7]: 70 71 72 73 74 75 76 77
    // to:
    // a0:    00 10 01 11 02 12 03 13  04 14 05 15 06 16 07 17
    // a1:    20 30 21 31 22 32 23 33  24 34 25 35 26 36 27 37
    // a2:    40 50 41 51 42 52 43 53  44 54 45 55 46 56 47 57
    // a3:    60 70 61 71 62 72 63 73  64 74 65 75 66 76 67 77
    const __m128i a0 = _mm_unpacklo_epi8(in[0], in[1]);
    const __m128i a1 = _mm_unpacklo_epi8(in[2], in[3]);
    const __m128i a2 = _mm_unpacklo_epi8(in[4], in[5]);
    const __m128i a3 = _mm_unpacklo_epi8(in[6], in[7]);

    // Unpack 16 bit elements resulting in:
    // b0: 00 10 20 30 01 11 21 31  02 12 22 32 03 13 23 33
    // b1: 40 50 60 70 41 51 61 71  42 52 62 72 43 53 63 73
    // b2: 04 14 24 34 05 15 25 35  06 16 26 36 07 17 27 37
    // b3: 44 54 64 74 45 55 65 75  46 56 66 76 47 57 67 77
    const __m128i b0 = _mm_unpacklo_epi16(a0, a1);
    const __m128i b1 = _mm_unpackhi_epi16(a0, a1);
    const __m128i b2 = _mm_unpacklo_epi16(a2, a3);
    const __m128i b3 = _mm_unpackhi_epi16(a2, a3);

    // Unpack 32 bit elements resulting in:
    // c0: 00 10 20 30 40 50 60 70  01 11 21 31 41 51 61 71
    // c1: 02 12 22 32 42 52 62 72  03 13 23 33 43 53 63 73
    // c2: 04 14 24 34 44 54 64 74  05 15 25 35 45 55 65 75
    // c3: 06 16 26 36 46 56 66 76  07 17 27 37 47 57 67 77
    out[0] = _mm_unpacklo_epi32(b0, b2);
    out[1] = _mm_unpackhi_epi32(b0, b2);
    out[2] = _mm_unpacklo_epi32(b1, b3);
    out[3] = _mm_unpackhi_epi32(b1, b3);
}

static INLINE void transpose_8bit_16x8(const __m128i *const in, __m128i *const out) {
    // Unpack 8 bit elements. Goes from:
    // in[0]: 00 01 02 03 04 05 06 07  08 09 0A 0B 0C 0D 0E 0F
    // in[1]: 10 11 12 13 14 15 16 17  18 19 1A 1B 1C 1D 1E 1F
    // in[2]: 20 21 22 23 24 25 26 27  28 29 2A 2B 2C 2D 2E 2F
    // in[3]: 30 31 32 33 34 35 36 37  38 39 3A 3B 3C 3D 3E 3F
    // in[4]: 40 41 42 43 44 45 46 47  48 49 4A 4B 4C 4D 4E 4F
    // in[5]: 50 51 52 53 54 55 56 57  58 59 5A 5B 5C 5D 5E 5F
    // in[6]: 60 61 62 63 64 65 66 67  68 69 6A 6B 6C 6D 6E 6F
    // in[7]: 70 71 72 73 74 75 76 77  78 79 7A 7B 7C 7D 7E 7F
    // to:
    // a0:    00 10 01 11 02 12 03 13  04 14 05 15 06 16 07 17
    // a1:    20 30 21 31 22 32 23 33  24 34 25 35 26 36 27 37
    // a2:    40 50 41 51 42 52 43 53  44 54 45 55 46 56 47 57
    // a3:    60 70 61 71 62 72 63 73  64 74 65 75 66 76 67 77
    // a4:    08 18 09 19 0A 1A 0B 1B  08 18 09 19 0A 1A 0B 1B
    // a5:    28 38 29 39 2A 3A 2B 3B  28 38 29 39 2A 3A 2B 3B
    // a6:    48 58 49 59 4A 5A 4B 5B  48 58 49 59 4A 5A 4B 5B
    // a7:    68 78 69 79 6A 7A 6B 7B  68 78 69 79 6A 7A 6B 7B
    const __m128i a0 = _mm_unpacklo_epi8(in[0], in[1]);
    const __m128i a1 = _mm_unpacklo_epi8(in[2], in[3]);
    const __m128i a2 = _mm_unpacklo_epi8(in[4], in[5]);
    const __m128i a3 = _mm_unpacklo_epi8(in[6], in[7]);
    const __m128i a4 = _mm_unpackhi_epi8(in[0], in[1]);
    const __m128i a5 = _mm_unpackhi_epi8(in[2], in[3]);
    const __m128i a6 = _mm_unpackhi_epi8(in[4], in[5]);
    const __m128i a7 = _mm_unpackhi_epi8(in[6], in[7]);

    // Unpack 16 bit elements resulting in:
    // b0: 00 10 20 30 01 11 21 31  02 12 22 32 03 13 23 33
    // b1: 40 50 60 70 41 51 61 71  42 52 62 72 43 53 63 73
    // b2: 04 14 24 34 05 15 25 35  06 16 26 36 07 17 27 37
    // b3: 44 54 64 74 45 55 65 75  46 56 66 76 47 57 67 77
    // b4: 08 18 28 38 09 19 29 39  0A 1A 2A 3A 0B 1B 2B 3B
    // b5: 48 58 68 78 49 59 69 79  4A 5A 6A 7A 4B 5B 6B 7B
    // b6: 0C 1C 2C 3C 0D 1D 2D 3D  0E 1E 2E 3E 0F 1F 2F 3F
    // b7: 4C 5C 6C 7C 4D 5D 6D 7D  4E 5E 6E 7E 4F 5F 6F 7F
    const __m128i b0 = _mm_unpacklo_epi16(a0, a1);
    const __m128i b1 = _mm_unpackhi_epi16(a0, a1);
    const __m128i b2 = _mm_unpacklo_epi16(a2, a3);
    const __m128i b3 = _mm_unpackhi_epi16(a2, a3);
    const __m128i b4 = _mm_unpacklo_epi16(a4, a5);
    const __m128i b5 = _mm_unpackhi_epi16(a4, a5);
    const __m128i b6 = _mm_unpacklo_epi16(a6, a7);
    const __m128i b7 = _mm_unpackhi_epi16(a6, a7);

    // Unpack 32 bit elements resulting in:
    // c0: 00 10 20 30 40 50 60 70  01 11 21 31 41 51 61 71
    // c1: 02 12 22 32 42 52 62 72  03 13 23 33 43 53 63 73
    // c2: 04 14 24 34 44 54 64 74  05 15 25 35 45 55 65 75
    // c3: 06 16 26 36 46 56 66 76  07 17 27 37 47 57 67 77
    // c4: 08 18 28 38 48 58 68 78  09 19 29 39 49 59 69 79
    // c5: 0A 1A 2A 3A 4A 5A 6A 7A  0B 1B 2B 3B 4B 5B 6B 7B
    // c6: 0C 1C 2C 3C 4C 5C 6C 7C  0D 1D 2D 3D 4D 5D 6D 7D
    // c7: 0E 1E 2E 3E 4E 5E 6E 7E  0F 1F 2F 3F 4F 5F 6F 7F
    out[0] = _mm_unpacklo_epi32(b0, b2);
    out[1] = _mm_unpackhi_epi32(b0, b2);
    out[2] = _mm_unpacklo_epi32(b1, b3);
    out[3] = _mm_unpackhi_epi32(b1, b3);
    out[4] = _mm_unpacklo_epi32(b4, b6);
    out[5] = _mm_unpackhi_epi32(b4, b6);
    out[6] = _mm_unpacklo_epi32(b5, b7);
    out[7] = _mm_unpackhi_epi32(b5, b7);
}

static INLINE void transpose_8bit_16x16_sse2(const __m128i *const in, __m128i *const out) {
    __m128i w0, w1, w2, w3, w4, w5, w6, w7, w8, w9;
    __m128i w10, w11, w12, w13, w14, w15;

    w0 = _mm_unpacklo_epi8(in[0], in[1]);
    w1 = _mm_unpacklo_epi8(in[2], in[3]);
    w2 = _mm_unpacklo_epi8(in[4], in[5]);
    w3 = _mm_unpacklo_epi8(in[6], in[7]);

    w8  = _mm_unpacklo_epi8(in[8], in[9]);
    w9  = _mm_unpacklo_epi8(in[10], in[11]);
    w10 = _mm_unpacklo_epi8(in[12], in[13]);
    w11 = _mm_unpacklo_epi8(in[14], in[15]);

    w4  = _mm_unpacklo_epi16(w0, w1);
    w5  = _mm_unpacklo_epi16(w2, w3);
    w12 = _mm_unpacklo_epi16(w8, w9);
    w13 = _mm_unpacklo_epi16(w10, w11);

    w6  = _mm_unpacklo_epi32(w4, w5);
    w7  = _mm_unpackhi_epi32(w4, w5);
    w14 = _mm_unpacklo_epi32(w12, w13);
    w15 = _mm_unpackhi_epi32(w12, w13);

    // Store first 4-line result
    out[0] = _mm_unpacklo_epi64(w6, w14);
    out[1] = _mm_unpackhi_epi64(w6, w14);
    out[2] = _mm_unpacklo_epi64(w7, w15);
    out[3] = _mm_unpackhi_epi64(w7, w15);

    w4  = _mm_unpackhi_epi16(w0, w1);
    w5  = _mm_unpackhi_epi16(w2, w3);
    w12 = _mm_unpackhi_epi16(w8, w9);
    w13 = _mm_unpackhi_epi16(w10, w11);

    w6  = _mm_unpacklo_epi32(w4, w5);
    w7  = _mm_unpackhi_epi32(w4, w5);
    w14 = _mm_unpacklo_epi32(w12, w13);
    w15 = _mm_unpackhi_epi32(w12, w13);

    // Store second 4-line result
    out[4] = _mm_unpacklo_epi64(w6, w14);
    out[5] = _mm_unpackhi_epi64(w6, w14);
    out[6] = _mm_unpacklo_epi64(w7, w15);
    out[7] = _mm_unpackhi_epi64(w7, w15);

    // upper half
    w0 = _mm_unpackhi_epi8(in[0], in[1]);
    w1 = _mm_unpackhi_epi8(in[2], in[3]);
    w2 = _mm_unpackhi_epi8(in[4], in[5]);
    w3 = _mm_unpackhi_epi8(in[6], in[7]);

    w8  = _mm_unpackhi_epi8(in[8], in[9]);
    w9  = _mm_unpackhi_epi8(in[10], in[11]);
    w10 = _mm_unpackhi_epi8(in[12], in[13]);
    w11 = _mm_unpackhi_epi8(in[14], in[15]);

    w4  = _mm_unpacklo_epi16(w0, w1);
    w5  = _mm_unpacklo_epi16(w2, w3);
    w12 = _mm_unpacklo_epi16(w8, w9);
    w13 = _mm_unpacklo_epi16(w10, w11);

    w6  = _mm_unpacklo_epi32(w4, w5);
    w7  = _mm_unpackhi_epi32(w4, w5);
    w14 = _mm_unpacklo_epi32(w12, w13);
    w15 = _mm_unpackhi_epi32(w12, w13);

    // Store first 4-line result
    out[8]  = _mm_unpacklo_epi64(w6, w14);
    out[9]  = _mm_unpackhi_epi64(w6, w14);
    out[10] = _mm_unpacklo_epi64(w7, w15);
    out[11] = _mm_unpackhi_epi64(w7, w15);

    w4  = _mm_unpackhi_epi16(w0, w1);
    w5  = _mm_unpackhi_epi16(w2, w3);
    w12 = _mm_unpackhi_epi16(w8, w9);
    w13 = _mm_unpackhi_epi16(w10, w11);

    w6  = _mm_unpacklo_epi32(w4, w5);
    w7  = _mm_unpackhi_epi32(w4, w5);
    w14 = _mm_unpacklo_epi32(w12, w13);
    w15 = _mm_unpackhi_epi32(w12, w13);

    // Store second 4-line result
    out[12] = _mm_unpacklo_epi64(w6, w14);
    out[13] = _mm_unpackhi_epi64(w6, w14);
    out[14] = _mm_unpacklo_epi64(w7, w15);
    out[15] = _mm_unpackhi_epi64(w7, w15);
}

static INLINE void transpose_16bit_4x4(const __m128i *const in, __m128i *const out) {
    // Unpack 16 bit elements. Goes from:
    // in[0]: 00 01 02 03  XX XX XX XX
    // in[1]: 10 11 12 13  XX XX XX XX
    // in[2]: 20 21 22 23  XX XX XX XX
    // in[3]: 30 31 32 33  XX XX XX XX
    // to:
    // a0:    00 10 01 11  02 12 03 13
    // a1:    20 30 21 31  22 32 23 33
    const __m128i a0 = _mm_unpacklo_epi16(in[0], in[1]);
    const __m128i a1 = _mm_unpacklo_epi16(in[2], in[3]);

    // Unpack 32 bit elements resulting in:
    // out[0]: 00 10 20 30
    // out[1]: 01 11 21 31
    // out[2]: 02 12 22 32
    // out[3]: 03 13 23 33
    out[0] = _mm_unpacklo_epi32(a0, a1);
    out[1] = _mm_srli_si128(out[0], 8);
    out[2] = _mm_unpackhi_epi32(a0, a1);
    out[3] = _mm_srli_si128(out[2], 8);
}

static INLINE void transpose_16bit_4x8(const __m128i *const in, __m128i *const out) {
    // Unpack 16 bit elements. Goes from:
    // in[0]: 00 01 02 03  XX XX XX XX
    // in[1]: 10 11 12 13  XX XX XX XX
    // in[2]: 20 21 22 23  XX XX XX XX
    // in[3]: 30 31 32 33  XX XX XX XX
    // in[4]: 40 41 42 43  XX XX XX XX
    // in[5]: 50 51 52 53  XX XX XX XX
    // in[6]: 60 61 62 63  XX XX XX XX
    // in[7]: 70 71 72 73  XX XX XX XX
    // to:
    // a0:    00 10 01 11  02 12 03 13
    // a1:    20 30 21 31  22 32 23 33
    // a2:    40 50 41 51  42 52 43 53
    // a3:    60 70 61 71  62 72 63 73
    const __m128i a0 = _mm_unpacklo_epi16(in[0], in[1]);
    const __m128i a1 = _mm_unpacklo_epi16(in[2], in[3]);
    const __m128i a2 = _mm_unpacklo_epi16(in[4], in[5]);
    const __m128i a3 = _mm_unpacklo_epi16(in[6], in[7]);

    // Unpack 32 bit elements resulting in:
    // b0: 00 10 20 30  01 11 21 31
    // b1: 40 50 60 70  41 51 61 71
    // b2: 02 12 22 32  03 13 23 33
    // b3: 42 52 62 72  43 53 63 73
    const __m128i b0 = _mm_unpacklo_epi32(a0, a1);
    const __m128i b1 = _mm_unpacklo_epi32(a2, a3);
    const __m128i b2 = _mm_unpackhi_epi32(a0, a1);
    const __m128i b3 = _mm_unpackhi_epi32(a2, a3);

    // Unpack 64 bit elements resulting in:
    // out[0]: 00 10 20 30  40 50 60 70
    // out[1]: 01 11 21 31  41 51 61 71
    // out[2]: 02 12 22 32  42 52 62 72
    // out[3]: 03 13 23 33  43 53 63 73
    out[0] = _mm_unpacklo_epi64(b0, b1);
    out[1] = _mm_unpackhi_epi64(b0, b1);
    out[2] = _mm_unpacklo_epi64(b2, b3);
    out[3] = _mm_unpackhi_epi64(b2, b3);
}

static INLINE void transpose_16bit_8x4(const __m128i *const in, __m128i *const out) {
    // Unpack 16 bit elements. Goes from:
    // in[0]: 00 01 02 03  04 05 06 07
    // in[1]: 10 11 12 13  14 15 16 17
    // in[2]: 20 21 22 23  24 25 26 27
    // in[3]: 30 31 32 33  34 35 36 37

    // to:
    // a0:    00 10 01 11  02 12 03 13
    // a1:    20 30 21 31  22 32 23 33
    // a4:    04 14 05 15  06 16 07 17
    // a5:    24 34 25 35  26 36 27 37
    const __m128i a0 = _mm_unpacklo_epi16(in[0], in[1]);
    const __m128i a1 = _mm_unpacklo_epi16(in[2], in[3]);
    const __m128i a4 = _mm_unpackhi_epi16(in[0], in[1]);
    const __m128i a5 = _mm_unpackhi_epi16(in[2], in[3]);

    // Unpack 32 bit elements resulting in:
    // b0: 00 10 20 30  01 11 21 31
    // b2: 04 14 24 34  05 15 25 35
    // b4: 02 12 22 32  03 13 23 33
    // b6: 06 16 26 36  07 17 27 37
    const __m128i b0 = _mm_unpacklo_epi32(a0, a1);
    const __m128i b2 = _mm_unpacklo_epi32(a4, a5);
    const __m128i b4 = _mm_unpackhi_epi32(a0, a1);
    const __m128i b6 = _mm_unpackhi_epi32(a4, a5);

    // Unpack 64 bit elements resulting in:
    // out[0]: 00 10 20 30  XX XX XX XX
    // out[1]: 01 11 21 31  XX XX XX XX
    // out[2]: 02 12 22 32  XX XX XX XX
    // out[3]: 03 13 23 33  XX XX XX XX
    // out[4]: 04 14 24 34  XX XX XX XX
    // out[5]: 05 15 25 35  XX XX XX XX
    // out[6]: 06 16 26 36  XX XX XX XX
    // out[7]: 07 17 27 37  XX XX XX XX
    const __m128i zeros = _mm_setzero_si128();
    out[0]              = _mm_unpacklo_epi64(b0, zeros);
    out[1]              = _mm_unpackhi_epi64(b0, zeros);
    out[2]              = _mm_unpacklo_epi64(b4, zeros);
    out[3]              = _mm_unpackhi_epi64(b4, zeros);
    out[4]              = _mm_unpacklo_epi64(b2, zeros);
    out[5]              = _mm_unpackhi_epi64(b2, zeros);
    out[6]              = _mm_unpacklo_epi64(b6, zeros);
    out[7]              = _mm_unpackhi_epi64(b6, zeros);
}

static INLINE void transpose_16bit_8x8(const __m128i *const in, __m128i *const out) {
    // Unpack 16 bit elements. Goes from:
    // in[0]: 00 01 02 03  04 05 06 07
    // in[1]: 10 11 12 13  14 15 16 17
    // in[2]: 20 21 22 23  24 25 26 27
    // in[3]: 30 31 32 33  34 35 36 37
    // in[4]: 40 41 42 43  44 45 46 47
    // in[5]: 50 51 52 53  54 55 56 57
    // in[6]: 60 61 62 63  64 65 66 67
    // in[7]: 70 71 72 73  74 75 76 77
    // to:
    // a0:    00 10 01 11  02 12 03 13
    // a1:    20 30 21 31  22 32 23 33
    // a2:    40 50 41 51  42 52 43 53
    // a3:    60 70 61 71  62 72 63 73
    // a4:    04 14 05 15  06 16 07 17
    // a5:    24 34 25 35  26 36 27 37
    // a6:    44 54 45 55  46 56 47 57
    // a7:    64 74 65 75  66 76 67 77
    const __m128i a0 = _mm_unpacklo_epi16(in[0], in[1]);
    const __m128i a1 = _mm_unpacklo_epi16(in[2], in[3]);
    const __m128i a2 = _mm_unpacklo_epi16(in[4], in[5]);
    const __m128i a3 = _mm_unpacklo_epi16(in[6], in[7]);
    const __m128i a4 = _mm_unpackhi_epi16(in[0], in[1]);
    const __m128i a5 = _mm_unpackhi_epi16(in[2], in[3]);
    const __m128i a6 = _mm_unpackhi_epi16(in[4], in[5]);
    const __m128i a7 = _mm_unpackhi_epi16(in[6], in[7]);

    // Unpack 32 bit elements resulting in:
    // b0: 00 10 20 30  01 11 21 31
    // b1: 40 50 60 70  41 51 61 71
    // b2: 04 14 24 34  05 15 25 35
    // b3: 44 54 64 74  45 55 65 75
    // b4: 02 12 22 32  03 13 23 33
    // b5: 42 52 62 72  43 53 63 73
    // b6: 06 16 26 36  07 17 27 37
    // b7: 46 56 66 76  47 57 67 77
    const __m128i b0 = _mm_unpacklo_epi32(a0, a1);
    const __m128i b1 = _mm_unpacklo_epi32(a2, a3);
    const __m128i b2 = _mm_unpacklo_epi32(a4, a5);
    const __m128i b3 = _mm_unpacklo_epi32(a6, a7);
    const __m128i b4 = _mm_unpackhi_epi32(a0, a1);
    const __m128i b5 = _mm_unpackhi_epi32(a2, a3);
    const __m128i b6 = _mm_unpackhi_epi32(a4, a5);
    const __m128i b7 = _mm_unpackhi_epi32(a6, a7);

    // Unpack 64 bit elements resulting in:
    // out[0]: 00 10 20 30  40 50 60 70
    // out[1]: 01 11 21 31  41 51 61 71
    // out[2]: 02 12 22 32  42 52 62 72
    // out[3]: 03 13 23 33  43 53 63 73
    // out[4]: 04 14 24 34  44 54 64 74
    // out[5]: 05 15 25 35  45 55 65 75
    // out[6]: 06 16 26 36  46 56 66 76
    // out[7]: 07 17 27 37  47 57 67 77
    out[0] = _mm_unpacklo_epi64(b0, b1);
    out[1] = _mm_unpackhi_epi64(b0, b1);
    out[2] = _mm_unpacklo_epi64(b4, b5);
    out[3] = _mm_unpackhi_epi64(b4, b5);
    out[4] = _mm_unpacklo_epi64(b2, b3);
    out[5] = _mm_unpackhi_epi64(b2, b3);
    out[6] = _mm_unpacklo_epi64(b6, b7);
    out[7] = _mm_unpackhi_epi64(b6, b7);
}

// Transpose in-place
static INLINE void transpose_16bit_16x16(__m128i *const left, __m128i *const right) {
    __m128i tbuf[8];
    transpose_16bit_8x8(left, left);
    transpose_16bit_8x8(right, tbuf);
    transpose_16bit_8x8(left + 8, right);
    transpose_16bit_8x8(right + 8, right + 8);

    left[8]  = tbuf[0];
    left[9]  = tbuf[1];
    left[10] = tbuf[2];
    left[11] = tbuf[3];
    left[12] = tbuf[4];
    left[13] = tbuf[5];
    left[14] = tbuf[6];
    left[15] = tbuf[7];
}

static INLINE void transpose_32bit_4x4(const __m128i *const in, __m128i *const out) {
    // Unpack 32 bit elements. Goes from:
    // in[0]: 00 01 02 03
    // in[1]: 10 11 12 13
    // in[2]: 20 21 22 23
    // in[3]: 30 31 32 33
    // to:
    // a0:    00 10 01 11
    // a1:    20 30 21 31
    // a2:    02 12 03 13
    // a3:    22 32 23 33

    const __m128i a0 = _mm_unpacklo_epi32(in[0], in[1]);
    const __m128i a1 = _mm_unpacklo_epi32(in[2], in[3]);
    const __m128i a2 = _mm_unpackhi_epi32(in[0], in[1]);
    const __m128i a3 = _mm_unpackhi_epi32(in[2], in[3]);

    // Unpack 64 bit elements resulting in:
    // out[0]: 00 10 20 30
    // out[1]: 01 11 21 31
    // out[2]: 02 12 22 32
    // out[3]: 03 13 23 33
    out[0] = _mm_unpacklo_epi64(a0, a1);
    out[1] = _mm_unpackhi_epi64(a0, a1);
    out[2] = _mm_unpacklo_epi64(a2, a3);
    out[3] = _mm_unpackhi_epi64(a2, a3);
}

static INLINE void transpose_32bit_4x4x2(const __m128i *const in, __m128i *const out) {
    // Unpack 32 bit elements. Goes from:
    // in[0]: 00 01 02 03
    // in[1]: 10 11 12 13
    // in[2]: 20 21 22 23
    // in[3]: 30 31 32 33
    // in[4]: 04 05 06 07
    // in[5]: 14 15 16 17
    // in[6]: 24 25 26 27
    // in[7]: 34 35 36 37
    // to:
    // a0:    00 10 01 11
    // a1:    20 30 21 31
    // a2:    02 12 03 13
    // a3:    22 32 23 33
    // a4:    04 14 05 15
    // a5:    24 34 25 35
    // a6:    06 16 07 17
    // a7:    26 36 27 37
    const __m128i a0 = _mm_unpacklo_epi32(in[0], in[1]);
    const __m128i a1 = _mm_unpacklo_epi32(in[2], in[3]);
    const __m128i a2 = _mm_unpackhi_epi32(in[0], in[1]);
    const __m128i a3 = _mm_unpackhi_epi32(in[2], in[3]);
    const __m128i a4 = _mm_unpacklo_epi32(in[4], in[5]);
    const __m128i a5 = _mm_unpacklo_epi32(in[6], in[7]);
    const __m128i a6 = _mm_unpackhi_epi32(in[4], in[5]);
    const __m128i a7 = _mm_unpackhi_epi32(in[6], in[7]);

    // Unpack 64 bit elements resulting in:
    // out[0]: 00 10 20 30
    // out[1]: 01 11 21 31
    // out[2]: 02 12 22 32
    // out[3]: 03 13 23 33
    // out[4]: 04 14 24 34
    // out[5]: 05 15 25 35
    // out[6]: 06 16 26 36
    // out[7]: 07 17 27 37
    out[0] = _mm_unpacklo_epi64(a0, a1);
    out[1] = _mm_unpackhi_epi64(a0, a1);
    out[2] = _mm_unpacklo_epi64(a2, a3);
    out[3] = _mm_unpackhi_epi64(a2, a3);
    out[4] = _mm_unpacklo_epi64(a4, a5);
    out[5] = _mm_unpackhi_epi64(a4, a5);
    out[6] = _mm_unpacklo_epi64(a6, a7);
    out[7] = _mm_unpackhi_epi64(a6, a7);
}

static INLINE void transpose_32bit_8x4(const __m128i *const in, __m128i *const out) {
    // Unpack 32 bit elements. Goes from:
    // in[0]: 00 01 02 03
    // in[1]: 04 05 06 07
    // in[2]: 10 11 12 13
    // in[3]: 14 15 16 17
    // in[4]: 20 21 22 23
    // in[5]: 24 25 26 27
    // in[6]: 30 31 32 33
    // in[7]: 34 35 36 37
    // to:
    // a0: 00 10 01 11
    // a1: 20 30 21 31
    // a2: 02 12 03 13
    // a3: 22 32 23 33
    // a4: 04 14 05 15
    // a5: 24 34 25 35
    // a6: 06 16 07 17
    // a7: 26 36 27 37
    const __m128i a0 = _mm_unpacklo_epi32(in[0], in[2]);
    const __m128i a1 = _mm_unpacklo_epi32(in[4], in[6]);
    const __m128i a2 = _mm_unpackhi_epi32(in[0], in[2]);
    const __m128i a3 = _mm_unpackhi_epi32(in[4], in[6]);
    const __m128i a4 = _mm_unpacklo_epi32(in[1], in[3]);
    const __m128i a5 = _mm_unpacklo_epi32(in[5], in[7]);
    const __m128i a6 = _mm_unpackhi_epi32(in[1], in[3]);
    const __m128i a7 = _mm_unpackhi_epi32(in[5], in[7]);

    // Unpack 64 bit elements resulting in:
    // out[0]: 00 10 20 30
    // out[1]: 01 11 21 31
    // out[2]: 02 12 22 32
    // out[3]: 03 13 23 33
    // out[4]: 04 14 24 34
    // out[5]: 05 15 25 35
    // out[6]: 06 16 26 36
    // out[7]: 07 17 27 37
    out[0] = _mm_unpacklo_epi64(a0, a1);
    out[1] = _mm_unpackhi_epi64(a0, a1);
    out[2] = _mm_unpacklo_epi64(a2, a3);
    out[3] = _mm_unpackhi_epi64(a2, a3);
    out[4] = _mm_unpacklo_epi64(a4, a5);
    out[5] = _mm_unpackhi_epi64(a4, a5);
    out[6] = _mm_unpacklo_epi64(a6, a7);
    out[7] = _mm_unpackhi_epi64(a6, a7);
}

#endif // AOM_DSP_X86_TRANSPOSE_SSE2_H_
