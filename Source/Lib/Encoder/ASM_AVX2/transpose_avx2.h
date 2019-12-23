/*
 * Copyright(c) 2019 Intel Corporation
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

#ifndef AOM_DSP_X86_TRANSPOSE_AVX2_H_
#define AOM_DSP_X86_TRANSPOSE_AVX2_H_

#include <immintrin.h> // AVX2
#include "EbDefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif

void transpose_8bit_16x16_reg128bit_instance_avx2(const __m128i *const in, __m128i *const out);

void transpose_32bit_8x8_reg256bit_instance_avx2(const __m256i *const in, __m256i *const out);

void transpose_64bit_4x4_reg256bit_instance_avx2(const __m256i *const in, __m256i *const out);

void transpose_64bit_4x6_reg256bit_instance_avx2(const __m256i *const in, __m256i *const out);

void transpose_64bit_4x8_reg256bit_instance_avx2(const __m256i *const in, __m256i *const out);

#ifdef __cplusplus
}
#endif

static INLINE __m256i _mm256_unpacklo_epi128(const __m256i in0, const __m256i in1) {
    return _mm256_inserti128_si256(in0, _mm256_extracti128_si256(in1, 0), 1);
}

static INLINE __m256i _mm256_unpackhi_epi128(const __m256i in0, const __m256i in1) {
    return _mm256_inserti128_si256(in1, _mm256_extracti128_si256(in0, 1), 0);
}

static INLINE void transpose_8bit_16x16_reg128bit_avx2(const __m128i *const in,
    __m128i *const out) {
    // Combine to 256 bit registers. Goes from:
    // in[ 0]: 00 01 02 03 04 05 06 07  08 09 0A 0B 0C 0D 0E 0F
    // in[ 1]: 10 11 12 13 14 15 16 17  18 19 1A 1B 1C 1D 1E 1F
    // in[ 2]: 20 21 22 23 24 25 26 27  28 29 2A 2B 2C 2D 2E 2F
    // in[ 3]: 30 31 32 33 34 35 36 37  38 39 3A 3B 3C 3D 3E 3F
    // in[ 4]: 40 41 42 43 44 45 46 47  48 49 4A 4B 4C 4D 4E 4F
    // in[ 5]: 50 51 52 53 54 55 56 57  58 59 5A 5B 5C 5D 5E 5F
    // in[ 6]: 60 61 62 63 64 65 66 67  68 69 6A 6B 6C 6D 6E 6F
    // in[ 7]: 70 71 72 73 74 75 76 77  78 79 7A 7B 7C 7D 7E 7F
    // in[ 8]: 80 81 82 83 84 85 86 87  88 89 8A 8B 8C 8D 8E 8F
    // in[ 9]: 90 91 92 93 94 95 96 97  98 99 9A 9B 9C 9D 9E 9F
    // in[10]: a0 A1 A2 A3 A4 A5 A6 A7  A8 A9 AA AB AC AD AE AF
    // in[11]: b0 B1 B2 B3 B4 B5 B6 B7  B8 B9 BA BB BC BD BE BF
    // in[12]: c0 C1 C2 C3 C4 C5 C6 C7  C8 C9 CA CB CC CD CE CF
    // in[13]: d0 D1 D2 D3 D4 D5 D6 D7  D8 D9 DA DB DC DD DE DF
    // in[14]: E0 E1 E2 E3 E4 E5 E6 E7  E8 E9 EA EB EC ED EE EF
    // in[15]: F0 F1 F2 F3 F4 F5 F6 F7  F8 F9 FA FB FC FD FE FF
    // to:
    // a0: 00 01 02 03 04 05 06 07  08 09 0A 0B 0C 0D 0E 0F   80 81 82 83 84 85 86 87  88 89 8A 8B 8C 8D 8E 8F
    // a1: 10 11 12 13 14 15 16 17  18 19 1A 1B 1C 1D 1E 1F   90 91 92 93 94 95 96 97  98 99 9A 9B 9C 9D 9E 9F
    // a2: 20 21 22 23 24 25 26 27  28 29 2A 2B 2C 2D 2E 2F   a0 A1 A2 A3 A4 A5 A6 A7  A8 A9 AA AB AC AD AE AF
    // a3: 30 31 32 33 34 35 36 37  38 39 3A 3B 3C 3D 3E 3F   b0 B1 B2 B3 B4 B5 B6 B7  B8 B9 BA BB BC BD BE BF
    // a4: 40 41 42 43 44 45 46 47  48 49 4A 4B 4C 4D 4E 4F   c0 C1 C2 C3 C4 C5 C6 C7  C8 C9 CA CB CC CD CE CF
    // a5: 50 51 52 53 54 55 56 57  58 59 5A 5B 5C 5D 5E 5F   d0 D1 D2 D3 D4 D5 D6 D7  D8 D9 DA DB DC DD DE DF
    // a6: 60 61 62 63 64 65 66 67  68 69 6A 6B 6C 6D 6E 6F   E0 E1 E2 E3 E4 E5 E6 E7  E8 E9 EA EB EC ED EE EF
    // a7: 70 71 72 73 74 75 76 77  78 79 7A 7B 7C 7D 7E 7F   F0 F1 F2 F3 F4 F5 F6 F7  F8 F9 FA FB FC FD FE FF
    const __m256i a0 = _mm256_inserti128_si256(_mm256_castsi128_si256(in[0]), in[8], 1);
    const __m256i a1 = _mm256_inserti128_si256(_mm256_castsi128_si256(in[1]), in[9], 1);
    const __m256i a2 = _mm256_inserti128_si256(_mm256_castsi128_si256(in[2]), in[10], 1);
    const __m256i a3 = _mm256_inserti128_si256(_mm256_castsi128_si256(in[3]), in[11], 1);
    const __m256i a4 = _mm256_inserti128_si256(_mm256_castsi128_si256(in[4]), in[12], 1);
    const __m256i a5 = _mm256_inserti128_si256(_mm256_castsi128_si256(in[5]), in[13], 1);
    const __m256i a6 = _mm256_inserti128_si256(_mm256_castsi128_si256(in[6]), in[14], 1);
    const __m256i a7 = _mm256_inserti128_si256(_mm256_castsi128_si256(in[7]), in[15], 1);

    // Unpack 8 bit elements resulting in:
    // b0: 00 10 01 11 02 12 03 13  04 14 05 15 06 16 07 17   80 90 81 91 82 92 83 93  84 94 85 95 86 96 87 97
    // b1: 20 30 21 31 22 32 23 33  24 34 25 35 26 36 27 37   a0 b0 A1 B1 A2 B2 A3 B3  A4 B4 A5 B5 A6 B6 A7 B7
    // b2: 40 50 41 51 42 52 43 53  44 54 45 55 46 56 47 57   c0 d0 C1 D1 C2 D2 C3 D3  C4 D4 C5 D5 C6 D6 C7 D7
    // b3: 60 70 61 71 62 72 63 73  64 74 65 75 66 76 67 77   E0 F0 E1 F1 E2 F2 E3 F3  E4 F4 E5 F5 E6 F6 E7 F7
    // b4: 08 18 09 19 0A 1A 0B 1B  0C 1C 0D 1D 0E 1E 0F 1F   88 98 89 99 8A 9A 8B 9B  8C 9C 8D 9D 8E 9E 8F 9F
    // b5: 28 38 29 39 2A 3A 2B 3B  2C 3C 2D 3D 2E 3E 2F 3F   A8 B8 A9 B9 AA BA AB BB  AC BC AD BD AE BE AF BF
    // b6: 48 58 49 59 4A 5A 4B 5B  4C 5C 4D 5D 4E 5E 4F 5F   C8 D8 C9 D9 CA DA CB DB  CC DC CD DD CE DE CF DF
    // b7: 68 78 69 79 6A 7A 6B 7B  6C 7C 6D 7D 6E 7E 6F 7F   E8 F8 E9 F9 EA FA EB FB  EC FC ED FD EE FE EF FF
    const __m256i b0 = _mm256_unpacklo_epi8(a0, a1);
    const __m256i b1 = _mm256_unpacklo_epi8(a2, a3);
    const __m256i b2 = _mm256_unpacklo_epi8(a4, a5);
    const __m256i b3 = _mm256_unpacklo_epi8(a6, a7);
    const __m256i b4 = _mm256_unpackhi_epi8(a0, a1);
    const __m256i b5 = _mm256_unpackhi_epi8(a2, a3);
    const __m256i b6 = _mm256_unpackhi_epi8(a4, a5);
    const __m256i b7 = _mm256_unpackhi_epi8(a6, a7);

    // Unpack 16 bit elements resulting in:
    // c0: 00 10 20 30 01 11 21 31  02 12 22 32 03 13 23 33   80 90 a0 b0 81 91 A1 B1  82 92 A2 B2 83 93 A3 B3
    // c1: 40 50 60 70 41 51 61 71  42 52 62 72 43 53 63 73   c0 d0 E0 F0 C1 D1 E1 F1  C2 D2 E2 F2 C3 D3 E3 F3
    // c2: 04 14 24 34 05 15 25 35  06 16 26 36 07 17 27 37   84 94 A4 B4 85 95 A5 B5  86 96 A6 B6 87 97 A7 B7
    // c3: 44 54 64 74 45 55 65 75  46 56 66 76 47 57 67 77   C4 D4 E4 F4 C5 D5 E5 F5  C6 D6 E6 F6 C7 D7 E7 F7
    // c4: 08 18 28 38 09 19 29 39  0A 1A 2A 3A 0B 1B 2B 3B   88 98 A8 B8 89 99 A9 B9  8A 9A AA BA 8B 9B AB BB
    // c5: 48 58 68 78 49 59 69 79  4A 5A 6A 7A 4B 5B 6B 7B   C8 D8 E8 F8 C9 D9 E9 F9  CA DA EA FA CB DB EB FB
    // c6: 0C 1C 2C 3C 0D 1D 2D 3D  0E 1E 2E 3E 0F 1F 2F 3F   8C 9C AC BC 8D 9D AD BD  8E 9E AE BE 8F 9F AF BF
    // c7: 4C 5C 6C 7C 4D 5D 6D 7D  4E 5E 6E 7E 4F 5F 6F 7F   CC DC EC FC CD DD ED FD  CE DE EE FE CF DF EF FF
    const __m256i c0 = _mm256_unpacklo_epi16(b0, b1);
    const __m256i c1 = _mm256_unpacklo_epi16(b2, b3);
    const __m256i c2 = _mm256_unpackhi_epi16(b0, b1);
    const __m256i c3 = _mm256_unpackhi_epi16(b2, b3);
    const __m256i c4 = _mm256_unpacklo_epi16(b4, b5);
    const __m256i c5 = _mm256_unpacklo_epi16(b6, b7);
    const __m256i c6 = _mm256_unpackhi_epi16(b4, b5);
    const __m256i c7 = _mm256_unpackhi_epi16(b6, b7);

    // Unpack 32 bit elements resulting in:
    // d0: 00 10 20 30 40 50 60 70  01 11 21 31 41 51 61 71   80 90 a0 b0 c0 d0 E0 F0  91 81 A1 B1 C1 D1 E1 F1
    // d1: 02 12 22 32 42 52 62 72  03 13 23 33 43 53 63 73   82 92 A2 B2 C2 D2 E2 F2  93 83 A3 B3 C3 D3 E3 F3
    // d2: 04 14 24 34 44 54 64 74  05 15 25 35 45 55 65 75   84 94 A4 B4 C4 D4 E4 F4  95 85 A5 B5 C5 D5 E5 F5
    // d3: 06 16 26 36 46 56 66 76  07 17 27 37 47 57 67 77   86 96 A6 B6 C6 D6 E6 F6  97 87 A7 B7 C7 D7 E7 F7
    // d4: 08 18 28 38 48 58 68 78  09 19 29 39 49 59 69 79   88 98 A8 B8 C8 D8 E8 F8  89 99 A9 B9 C9 D9 E9 F9
    // d5: 0A 1A 2A 3A 4A 5A 6A 7A  0B 1B 2B 3B 4B 5B 6B 7B   8A 9A AA BA CA DA EA FA  8B 9B AB BB CB DB EB FB
    // d6: 0C 1C 2C 3C 4C 5C 6C 7C  0D 1D 2D 3D 4D 5D 6D 7D   8C 9C AC BC CC DC EC FC  8D 9D AD BD CD DD ED FD
    // d7: 0E 1E 2E 3E 4E 5E 6E 7E  0F 1F 2F 3F 4F 5F 6F 7F   8E 9E AE BE CE DE EE FE  8F 9F AF BF CF DF EF FF
    const __m256i d0 = _mm256_unpacklo_epi32(c0, c1);
    const __m256i d1 = _mm256_unpackhi_epi32(c0, c1);
    const __m256i d2 = _mm256_unpacklo_epi32(c2, c3);
    const __m256i d3 = _mm256_unpackhi_epi32(c2, c3);
    const __m256i d4 = _mm256_unpacklo_epi32(c4, c5);
    const __m256i d5 = _mm256_unpackhi_epi32(c4, c5);
    const __m256i d6 = _mm256_unpacklo_epi32(c6, c7);
    const __m256i d7 = _mm256_unpackhi_epi32(c6, c7);

    // Permute 64 bit elements resulting in:
    // e0: 00 10 20 30 40 50 60 70  80 90 a0 b0 80 90 a0 b0   01 11 21 31 41 51 61 71  C1 D1 E1 F1 C1 D1 E1 F1
    // e1: 02 12 22 32 42 52 62 72  82 92 A2 B2 82 92 A2 B2   03 13 23 33 43 53 63 73  C3 D3 E3 F3 C3 D3 E3 F3
    // e2: 04 14 24 34 44 54 64 74  84 94 A4 B4 84 94 A4 B4   05 15 25 35 45 55 65 75  C5 D5 E5 F5 C5 D5 E5 F5
    // e3: 06 16 26 36 46 56 66 76  86 96 A6 B6 86 96 A6 B6   07 17 27 37 47 57 67 77  C7 D7 E7 F7 C7 D7 E7 F7
    // e4: 08 18 28 38 48 58 68 78  88 98 A8 B8 C8 D8 E8 F8   09 19 29 39 49 59 69 79  89 99 A9 B9 C9 D9 E9 F9
    // e5: 0A 1A 2A 3A 4A 5A 6A 7A  8A 9A AA BA CA DA EA FA   0B 1B 2B 3B 4B 5B 6B 7B  8B 9B AB BB CB DB EB FB
    // e6: 0C 1C 2C 3C 4C 5C 6C 7C  8C 9C AC BC CC DC EC FC   0D 1D 2D 3D 4D 5D 6D 7D  8D 9D AD BD CD DD ED FD
    // e7: 0E 1E 2E 3E 4E 5E 6E 7E  8E 9E AE BE CE DE EE FE   0F 1F 2F 3F 4F 5F 6F 7F  8F 9F AF BF CF DF EF FF
    const __m256i e0 = _mm256_permute4x64_epi64(d0, 0xD8);
    const __m256i e1 = _mm256_permute4x64_epi64(d1, 0xD8);
    const __m256i e2 = _mm256_permute4x64_epi64(d2, 0xD8);
    const __m256i e3 = _mm256_permute4x64_epi64(d3, 0xD8);
    const __m256i e4 = _mm256_permute4x64_epi64(d4, 0xD8);
    const __m256i e5 = _mm256_permute4x64_epi64(d5, 0xD8);
    const __m256i e6 = _mm256_permute4x64_epi64(d6, 0xD8);
    const __m256i e7 = _mm256_permute4x64_epi64(d7, 0xD8);

    out[0] = _mm256_castsi256_si128(e0);
    out[1] = _mm256_extracti128_si256(e0, 1);
    out[2] = _mm256_castsi256_si128(e1);
    out[3] = _mm256_extracti128_si256(e1, 1);
    out[4] = _mm256_castsi256_si128(e2);
    out[5] = _mm256_extracti128_si256(e2, 1);
    out[6] = _mm256_castsi256_si128(e3);
    out[7] = _mm256_extracti128_si256(e3, 1);
    out[8] = _mm256_castsi256_si128(e4);
    out[9] = _mm256_extracti128_si256(e4, 1);
    out[10] = _mm256_castsi256_si128(e5);
    out[11] = _mm256_extracti128_si256(e5, 1);
    out[12] = _mm256_castsi256_si128(e6);
    out[13] = _mm256_extracti128_si256(e6, 1);
    out[14] = _mm256_castsi256_si128(e7);
    out[15] = _mm256_extracti128_si256(e7, 1);
}
/* clang-format on */

static INLINE void transpose_32bit_8x8_avx2(const __m256i *const in, __m256i *const out) {
    // Unpack 32 bit elements. Goes from:
    // in[0]: 00 01 02 03  04 05 06 07
    // in[1]: 10 11 12 13  14 15 16 17
    // in[2]: 20 21 22 23  24 25 26 27
    // in[3]: 30 31 32 33  34 35 36 37
    // in[4]: 40 41 42 43  44 45 46 47
    // in[5]: 50 51 52 53  54 55 56 57
    // in[6]: 60 61 62 63  64 65 66 67
    // in[7]: 70 71 72 73  74 75 76 77
    // to:
    // a0:    00 10 01 11  04 14 05 15
    // a1:    20 30 21 31  24 34 25 35
    // a2:    40 50 41 51  44 54 45 55
    // a3:    60 70 61 71  64 74 65 75
    // a4:    02 12 03 13  06 16 07 17
    // a5:    22 32 23 33  26 36 27 37
    // a6:    42 52 43 53  46 56 47 57
    // a7:    62 72 63 73  66 76 67 77
    const __m256i a0 = _mm256_unpacklo_epi32(in[0], in[1]);
    const __m256i a1 = _mm256_unpacklo_epi32(in[2], in[3]);
    const __m256i a2 = _mm256_unpacklo_epi32(in[4], in[5]);
    const __m256i a3 = _mm256_unpacklo_epi32(in[6], in[7]);
    const __m256i a4 = _mm256_unpackhi_epi32(in[0], in[1]);
    const __m256i a5 = _mm256_unpackhi_epi32(in[2], in[3]);
    const __m256i a6 = _mm256_unpackhi_epi32(in[4], in[5]);
    const __m256i a7 = _mm256_unpackhi_epi32(in[6], in[7]);

    // Unpack 64 bit elements resulting in:
    // b0: 00 10 20 30  04 14 24 34
    // b1: 40 50 60 70  44 54 64 74
    // b2: 01 11 21 31  05 15 25 35
    // b3: 41 51 61 71  45 55 65 75
    // b4: 02 12 22 32  06 16 26 36
    // b5: 42 52 62 72  46 56 66 76
    // b6: 03 13 23 33  07 17 27 37
    // b7: 43 53 63 73  47 57 67 77
    const __m256i b0 = _mm256_unpacklo_epi64(a0, a1);
    const __m256i b1 = _mm256_unpacklo_epi64(a2, a3);
    const __m256i b2 = _mm256_unpackhi_epi64(a0, a1);
    const __m256i b3 = _mm256_unpackhi_epi64(a2, a3);
    const __m256i b4 = _mm256_unpacklo_epi64(a4, a5);
    const __m256i b5 = _mm256_unpacklo_epi64(a6, a7);
    const __m256i b6 = _mm256_unpackhi_epi64(a4, a5);
    const __m256i b7 = _mm256_unpackhi_epi64(a6, a7);

    // Unpack 128 bit elements resulting in:
    // out[0]: 00 10 20 30  40 50 60 70
    // out[1]: 01 11 21 31  41 51 61 71
    // out[2]: 02 12 22 32  42 52 62 72
    // out[3]: 03 13 23 33  43 53 63 73
    // out[4]: 04 14 24 34  44 54 64 74
    // out[5]: 05 15 25 35  45 55 65 75
    // out[6]: 06 16 26 36  46 56 66 76
    // out[7]: 07 17 27 37  47 57 67 77
    out[0] = _mm256_unpacklo_epi128(b0, b1);
    out[1] = _mm256_unpacklo_epi128(b2, b3);
    out[2] = _mm256_unpacklo_epi128(b4, b5);
    out[3] = _mm256_unpacklo_epi128(b6, b7);
    out[4] = _mm256_unpackhi_epi128(b0, b1);
    out[5] = _mm256_unpackhi_epi128(b2, b3);
    out[6] = _mm256_unpackhi_epi128(b4, b5);
    out[7] = _mm256_unpackhi_epi128(b6, b7);
}

static INLINE void transpose_64bit_4x4_avx2(const __m256i *const in, __m256i *const out) {
    // Unpack 32 bit elements. Goes from:
    // in[0]: 00 01 02 03
    // in[1]: 10 11 12 13
    // in[2]: 20 21 22 23
    // in[3]: 30 31 32 33
    // to:
    // a0:    00 10 02 12
    // a1:    20 30 22 32
    // a2:    01 11 03 13
    // a3:    21 31 23 33
    const __m256i a0 = _mm256_unpacklo_epi64(in[0], in[1]);
    const __m256i a1 = _mm256_unpacklo_epi64(in[2], in[3]);
    const __m256i a2 = _mm256_unpackhi_epi64(in[0], in[1]);
    const __m256i a3 = _mm256_unpackhi_epi64(in[2], in[3]);

    // Unpack 64 bit elements resulting in:
    // out[0]: 00 10 20 30
    // out[1]: 01 11 21 31
    // out[2]: 02 12 22 32
    // out[3]: 03 13 23 33
    out[0] = _mm256_inserti128_si256(a0, _mm256_extracti128_si256(a1, 0), 1);
    out[1] = _mm256_inserti128_si256(a2, _mm256_extracti128_si256(a3, 0), 1);
    out[2] = _mm256_inserti128_si256(a1, _mm256_extracti128_si256(a0, 1), 0);
    out[3] = _mm256_inserti128_si256(a3, _mm256_extracti128_si256(a2, 1), 0);
}

static INLINE void transpose_64bit_4x6_avx2(const __m256i *const in, __m256i *const out) {
    // Unpack 64 bit elements. Goes from:
    // in[0]: 00 01  02 03
    // in[1]: 10 11  12 13
    // in[2]: 20 21  22 23
    // in[3]: 30 31  32 33
    // in[4]: 40 41  42 43
    // in[5]: 50 51  52 53
    // to:
    // a0:    00 10  02 12
    // a1:    20 30  22 32
    // a2:    40 50  42 52
    // a4:    01 11  03 13
    // a5:    21 31  23 33
    // a6:    41 51  43 53
    const __m256i a0 = _mm256_unpacklo_epi64(in[0], in[1]);
    const __m256i a1 = _mm256_unpacklo_epi64(in[2], in[3]);
    const __m256i a2 = _mm256_unpacklo_epi64(in[4], in[5]);
    const __m256i a4 = _mm256_unpackhi_epi64(in[0], in[1]);
    const __m256i a5 = _mm256_unpackhi_epi64(in[2], in[3]);
    const __m256i a6 = _mm256_unpackhi_epi64(in[4], in[5]);

    // Unpack 128 bit elements resulting in:
    // b0: 00 10  20 30
    // b1: 40 50  40 50
    // b2: 01 11  21 31
    // b3: 41 51  41 51
    // b4: 02 12  22 32
    // b5: 42 52  42 52
    // b6: 03 13  23 33
    // b7: 43 53  43 53
    out[0] = _mm256_unpacklo_epi128(a0, a1);
    out[1] = _mm256_unpacklo_epi128(a2, a2);
    out[2] = _mm256_unpacklo_epi128(a4, a5);
    out[3] = _mm256_unpacklo_epi128(a6, a6);
    out[4] = _mm256_unpackhi_epi128(a0, a1);
    out[5] = _mm256_unpackhi_epi128(a2, a2);
    out[6] = _mm256_unpackhi_epi128(a4, a5);
    out[7] = _mm256_unpackhi_epi128(a6, a6);
}

static INLINE void transpose_64bit_4x8_avx2(const __m256i *const in, __m256i *const out) {
    // Unpack 64 bit elements. Goes from:
    // in[0]: 00 01  02 03
    // in[1]: 10 11  12 13
    // in[2]: 20 21  22 23
    // in[3]: 30 31  32 33
    // in[4]: 40 41  42 43
    // in[5]: 50 51  52 53
    // in[6]: 60 61  62 63
    // in[7]: 70 71  72 73
    // to:
    // a0:    00 10  02 12
    // a1:    20 30  22 32
    // a2:    40 50  42 52
    // a3:    60 70  62 72
    // a4:    01 11  03 13
    // a5:    21 31  23 33
    // a6:    41 51  43 53
    // a7:    61 71  63 73
    const __m256i a0 = _mm256_unpacklo_epi64(in[0], in[1]);
    const __m256i a1 = _mm256_unpacklo_epi64(in[2], in[3]);
    const __m256i a2 = _mm256_unpacklo_epi64(in[4], in[5]);
    const __m256i a3 = _mm256_unpacklo_epi64(in[6], in[7]);
    const __m256i a4 = _mm256_unpackhi_epi64(in[0], in[1]);
    const __m256i a5 = _mm256_unpackhi_epi64(in[2], in[3]);
    const __m256i a6 = _mm256_unpackhi_epi64(in[4], in[5]);
    const __m256i a7 = _mm256_unpackhi_epi64(in[6], in[7]);

    // Unpack 128 bit elements resulting in:
    // b0: 00 10  20 30
    // b1: 40 50  60 70
    // b2: 01 11  21 31
    // b3: 41 51  61 71
    // b4: 02 12  22 32
    // b5: 42 52  62 72
    // b6: 03 13  23 33
    // b7: 43 53  63 73
    out[0] = _mm256_unpacklo_epi128(a0, a1);
    out[1] = _mm256_unpacklo_epi128(a2, a3);
    out[2] = _mm256_unpacklo_epi128(a4, a5);
    out[3] = _mm256_unpacklo_epi128(a6, a7);
    out[4] = _mm256_unpackhi_epi128(a0, a1);
    out[5] = _mm256_unpackhi_epi128(a2, a3);
    out[6] = _mm256_unpackhi_epi128(a4, a5);
    out[7] = _mm256_unpackhi_epi128(a6, a7);
}

#endif // AOM_DSP_X86_TRANSPOSE_AVX2_H_
