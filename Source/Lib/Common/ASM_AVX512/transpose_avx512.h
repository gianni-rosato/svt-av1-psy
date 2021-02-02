/*
* Copyright(c) 2019 Intel Corporation
* Copyright (c) 2019, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#ifndef AOM_DSP_X86_TRANSPOSE_AVX512_H_
#define AOM_DSP_X86_TRANSPOSE_AVX512_H_

#include "EbDefinitions.h"
#if OPT_AVX512
#include <immintrin.h> // AVX2

#define ZERO (uint8_t)0U
#define ONE (uint8_t)1U
#define TWO (uint8_t)2U
#define THREE (uint8_t)3U

#define TRANSPOSE_4X4_AVX512(x0, x1, x2, x3, y0, y1, y2, y3) \
    do {                                                     \
        __m512i u0, u1, u2, u3;                              \
        u0 = _mm512_unpacklo_epi32(x0, x1);                  \
        u1 = _mm512_unpackhi_epi32(x0, x1);                  \
        u2 = _mm512_unpacklo_epi32(x2, x3);                  \
        u3 = _mm512_unpackhi_epi32(x2, x3);                  \
        y0 = _mm512_unpacklo_epi64(u0, u2);                  \
        y1 = _mm512_unpackhi_epi64(u0, u2);                  \
        y2 = _mm512_unpacklo_epi64(u1, u3);                  \
        y3 = _mm512_unpackhi_epi64(u1, u3);                  \
    } while (0)

static INLINE void transpose_16x16_avx512(int32_t stride, const __m512i *in, __m512i *out) {
    __m512i out1[16];
    TRANSPOSE_4X4_AVX512(in[0 * stride],
                         in[1 * stride],
                         in[2 * stride],
                         in[3 * stride],
                         out1[0],
                         out1[1],
                         out1[2],
                         out1[3]);
    TRANSPOSE_4X4_AVX512(in[4 * stride],
                         in[5 * stride],
                         in[6 * stride],
                         in[7 * stride],
                         out1[4],
                         out1[5],
                         out1[6],
                         out1[7]);
    TRANSPOSE_4X4_AVX512(in[8 * stride],
                         in[9 * stride],
                         in[10 * stride],
                         in[11 * stride],
                         out1[8],
                         out1[9],
                         out1[10],
                         out1[11]);
    TRANSPOSE_4X4_AVX512(in[12 * stride],
                         in[13 * stride],
                         in[14 * stride],
                         in[15 * stride],
                         out1[12],
                         out1[13],
                         out1[14],
                         out1[15]);

    __m128i *outptr = (__m128i *)(out + 0 * stride);

    //will get first row of transpose matrix from corresponding 4 vectors in out1
    outptr[0] = _mm512_extracti32x4_epi32(out1[0], ZERO);
    outptr[1] = _mm512_extracti32x4_epi32(out1[4], ZERO);
    outptr[2] = _mm512_extracti32x4_epi32(out1[8], ZERO);
    outptr[3] = _mm512_extracti32x4_epi32(out1[12], ZERO);

    //will get second row of transpose matrix from corresponding 4 vectors in out1
    outptr    = (__m128i *)(out + 1 * stride);
    outptr[0] = _mm512_extracti32x4_epi32(out1[1], ZERO);
    outptr[1] = _mm512_extracti32x4_epi32(out1[5], ZERO);
    outptr[2] = _mm512_extracti32x4_epi32(out1[9], ZERO);
    outptr[3] = _mm512_extracti32x4_epi32(out1[13], ZERO);

    //will get 3rd row of transpose matrix from corresponding 4 vectors in out1
    outptr    = (__m128i *)(out + 2 * stride);
    outptr[0] = _mm512_extracti32x4_epi32(out1[2], ZERO);
    outptr[1] = _mm512_extracti32x4_epi32(out1[6], ZERO);
    outptr[2] = _mm512_extracti32x4_epi32(out1[10], ZERO);
    outptr[3] = _mm512_extracti32x4_epi32(out1[14], ZERO);

    //will get 4th row of transpose matrix from corresponding 4 vectors in out1
    outptr    = (__m128i *)(out + 3 * stride);
    outptr[0] = _mm512_extracti32x4_epi32(out1[3], ZERO);
    outptr[1] = _mm512_extracti32x4_epi32(out1[7], ZERO);
    outptr[2] = _mm512_extracti32x4_epi32(out1[11], ZERO);
    outptr[3] = _mm512_extracti32x4_epi32(out1[15], ZERO);

    //will get 5th row of transpose matrix from corresponding 4 vectors in out1
    outptr    = (__m128i *)(out + 4 * stride);
    outptr[0] = _mm512_extracti32x4_epi32(out1[0], ONE);
    outptr[1] = _mm512_extracti32x4_epi32(out1[4], ONE);
    outptr[2] = _mm512_extracti32x4_epi32(out1[8], ONE);
    outptr[3] = _mm512_extracti32x4_epi32(out1[12], ONE);

    //will get 6th row of transpose matrix from corresponding 4 vectors in out1
    outptr    = (__m128i *)(out + 5 * stride);
    outptr[0] = _mm512_extracti32x4_epi32(out1[1], ONE);
    outptr[1] = _mm512_extracti32x4_epi32(out1[5], ONE);
    outptr[2] = _mm512_extracti32x4_epi32(out1[9], ONE);
    outptr[3] = _mm512_extracti32x4_epi32(out1[13], ONE);

    //will get 7th row of transpose matrix from corresponding 4 vectors in out1
    outptr    = (__m128i *)(out + 6 * stride);
    outptr[0] = _mm512_extracti32x4_epi32(out1[2], ONE);
    outptr[1] = _mm512_extracti32x4_epi32(out1[6], ONE);
    outptr[2] = _mm512_extracti32x4_epi32(out1[10], ONE);
    outptr[3] = _mm512_extracti32x4_epi32(out1[14], ONE);

    //will get 8th row of transpose matrix from corresponding 4 vectors in out1
    outptr    = (__m128i *)(out + 7 * stride);
    outptr[0] = _mm512_extracti32x4_epi32(out1[3], ONE);
    outptr[1] = _mm512_extracti32x4_epi32(out1[7], ONE);
    outptr[2] = _mm512_extracti32x4_epi32(out1[11], ONE);
    outptr[3] = _mm512_extracti32x4_epi32(out1[15], ONE);

    //will get 9th row of transpose matrix from corresponding 4 vectors in out1
    outptr    = (__m128i *)(out + 8 * stride);
    outptr[0] = _mm512_extracti32x4_epi32(out1[0], TWO);
    outptr[1] = _mm512_extracti32x4_epi32(out1[4], TWO);
    outptr[2] = _mm512_extracti32x4_epi32(out1[8], TWO);
    outptr[3] = _mm512_extracti32x4_epi32(out1[12], TWO);

    //will get 10th row of transpose matrix from corresponding 4 vectors in out1
    outptr    = (__m128i *)(out + 9 * stride);
    outptr[0] = _mm512_extracti32x4_epi32(out1[1], TWO);
    outptr[1] = _mm512_extracti32x4_epi32(out1[5], TWO);
    outptr[2] = _mm512_extracti32x4_epi32(out1[9], TWO);
    outptr[3] = _mm512_extracti32x4_epi32(out1[13], TWO);

    //will get 11th row of transpose matrix from corresponding 4 vectors in out1
    outptr    = (__m128i *)(out + 10 * stride);
    outptr[0] = _mm512_extracti32x4_epi32(out1[2], TWO);
    outptr[1] = _mm512_extracti32x4_epi32(out1[6], TWO);
    outptr[2] = _mm512_extracti32x4_epi32(out1[10], TWO);
    outptr[3] = _mm512_extracti32x4_epi32(out1[14], TWO);

    //will get 12th row of transpose matrix from corresponding 4 vectors in out1
    outptr    = (__m128i *)(out + 11 * stride);
    outptr[0] = _mm512_extracti32x4_epi32(out1[3], TWO);
    outptr[1] = _mm512_extracti32x4_epi32(out1[7], TWO);
    outptr[2] = _mm512_extracti32x4_epi32(out1[11], TWO);
    outptr[3] = _mm512_extracti32x4_epi32(out1[15], TWO);

    //will get 13th row of transpose matrix from corresponding 4 vectors in out1
    outptr    = (__m128i *)(out + 12 * stride);
    outptr[0] = _mm512_extracti32x4_epi32(out1[0], THREE);
    outptr[1] = _mm512_extracti32x4_epi32(out1[4], THREE);
    outptr[2] = _mm512_extracti32x4_epi32(out1[8], THREE);
    outptr[3] = _mm512_extracti32x4_epi32(out1[12], THREE);

    //will get 14th row of transpose matrix from corresponding 4 vectors in out1
    outptr    = (__m128i *)(out + 13 * stride);
    outptr[0] = _mm512_extracti32x4_epi32(out1[1], THREE);
    outptr[1] = _mm512_extracti32x4_epi32(out1[5], THREE);
    outptr[2] = _mm512_extracti32x4_epi32(out1[9], THREE);
    outptr[3] = _mm512_extracti32x4_epi32(out1[13], THREE);

    //will get 15th row of transpose matrix from corresponding 4 vectors in out1
    outptr    = (__m128i *)(out + 14 * stride);
    outptr[0] = _mm512_extracti32x4_epi32(out1[2], THREE);
    outptr[1] = _mm512_extracti32x4_epi32(out1[6], THREE);
    outptr[2] = _mm512_extracti32x4_epi32(out1[10], THREE);
    outptr[3] = _mm512_extracti32x4_epi32(out1[14], THREE);

    //will get 16th row of transpose matrix from corresponding 4 vectors in out1
    outptr    = (__m128i *)(out + 15 * stride);
    outptr[0] = _mm512_extracti32x4_epi32(out1[3], THREE);
    outptr[1] = _mm512_extracti32x4_epi32(out1[7], THREE);
    outptr[2] = _mm512_extracti32x4_epi32(out1[11], THREE);
    outptr[3] = _mm512_extracti32x4_epi32(out1[15], THREE);
}

static INLINE void transpose_16nx16n_avx512(int32_t txfm_size, const __m512i *input,
                                            __m512i *output) {
    const int32_t num_per_512 = 16;
    const int32_t row_size    = txfm_size;
    const int32_t col_size    = txfm_size / num_per_512;
    int32_t       r, c;

    // transpose each 16x16 block internally
    for (r = 0; r < row_size; r += 16) {
        for (c = 0; c < col_size; c++) {
            transpose_16x16_avx512(
                col_size, &input[r * col_size + c], &output[c * 16 * col_size + r / 16]);
        }
    }
}

static INLINE void transpose_16nx16m_inv_avx512(const __m512i *in, __m512i *out,
                                                const int32_t width, const int32_t height) {
    const int32_t numcol = height >> 4;
    const int32_t numrow = width >> 4;

    __m512i out1[16];

    for (int32_t j = 0; j < numrow; j++) {
        for (int32_t i = 0; i < numcol; i++) {
            TRANSPOSE_4X4_AVX512(in[i * width + j + (numrow * 0)],
                                 in[i * width + j + (numrow * 1)],
                                 in[i * width + j + (numrow * 2)],
                                 in[i * width + j + (numrow * 3)],
                                 out1[0],
                                 out1[1],
                                 out1[2],
                                 out1[3]);
            TRANSPOSE_4X4_AVX512(in[i * width + j + (numrow * 4)],
                                 in[i * width + j + (numrow * 5)],
                                 in[i * width + j + (numrow * 6)],
                                 in[i * width + j + (numrow * 7)],
                                 out1[4],
                                 out1[5],
                                 out1[6],
                                 out1[7]);
            TRANSPOSE_4X4_AVX512(in[i * width + j + (numrow * 8)],
                                 in[i * width + j + (numrow * 9)],
                                 in[i * width + j + (numrow * 10)],
                                 in[i * width + j + (numrow * 11)],
                                 out1[8],
                                 out1[9],
                                 out1[10],
                                 out1[11]);
            TRANSPOSE_4X4_AVX512(in[i * width + j + (numrow * 12)],
                                 in[i * width + j + (numrow * 13)],
                                 in[i * width + j + (numrow * 14)],
                                 in[i * width + j + (numrow * 15)],
                                 out1[12],
                                 out1[13],
                                 out1[14],
                                 out1[15]);

            __m128i *outptr = (__m128i *)(out + (j * height + i + (numcol * 0)));

            //will get first row of transpose matrix from corresponding 4 vectors in out1
            outptr[0] = _mm512_extracti32x4_epi32(out1[0], ZERO);
            outptr[1] = _mm512_extracti32x4_epi32(out1[4], ZERO);
            outptr[2] = _mm512_extracti32x4_epi32(out1[8], ZERO);
            outptr[3] = _mm512_extracti32x4_epi32(out1[12], ZERO);

            //will get second row of transpose matrix from corresponding 4 vectors in out1
            outptr    = (__m128i *)(out + (j * height + i + (numcol * 1)));
            outptr[0] = _mm512_extracti32x4_epi32(out1[1], ZERO);
            outptr[1] = _mm512_extracti32x4_epi32(out1[5], ZERO);
            outptr[2] = _mm512_extracti32x4_epi32(out1[9], ZERO);
            outptr[3] = _mm512_extracti32x4_epi32(out1[13], ZERO);

            //will get 3rd row of transpose matrix from corresponding 4 vectors in out1
            outptr    = (__m128i *)(out + (j * height + i + (numcol * 2)));
            outptr[0] = _mm512_extracti32x4_epi32(out1[2], ZERO);
            outptr[1] = _mm512_extracti32x4_epi32(out1[6], ZERO);
            outptr[2] = _mm512_extracti32x4_epi32(out1[10], ZERO);
            outptr[3] = _mm512_extracti32x4_epi32(out1[14], ZERO);

            //will get 4th row of transpose matrix from corresponding 4 vectors in out1
            outptr    = (__m128i *)(out + (j * height + i + (numcol * 3)));
            outptr[0] = _mm512_extracti32x4_epi32(out1[3], ZERO);
            outptr[1] = _mm512_extracti32x4_epi32(out1[7], ZERO);
            outptr[2] = _mm512_extracti32x4_epi32(out1[11], ZERO);
            outptr[3] = _mm512_extracti32x4_epi32(out1[15], ZERO);

            //will get 5th row of transpose matrix from corresponding 4 vectors in out1
            outptr    = (__m128i *)(out + (j * height + i + (numcol * 4)));
            outptr[0] = _mm512_extracti32x4_epi32(out1[0], ONE);
            outptr[1] = _mm512_extracti32x4_epi32(out1[4], ONE);
            outptr[2] = _mm512_extracti32x4_epi32(out1[8], ONE);
            outptr[3] = _mm512_extracti32x4_epi32(out1[12], ONE);

            //will get 6th row of transpose matrix from corresponding 4 vectors in out1
            outptr    = (__m128i *)(out + (j * height + i + (numcol * 5)));
            outptr[0] = _mm512_extracti32x4_epi32(out1[1], ONE);
            outptr[1] = _mm512_extracti32x4_epi32(out1[5], ONE);
            outptr[2] = _mm512_extracti32x4_epi32(out1[9], ONE);
            outptr[3] = _mm512_extracti32x4_epi32(out1[13], ONE);

            //will get 7th row of transpose matrix from corresponding 4 vectors in out1
            outptr    = (__m128i *)(out + (j * height + i + (numcol * 6)));
            outptr[0] = _mm512_extracti32x4_epi32(out1[2], ONE);
            outptr[1] = _mm512_extracti32x4_epi32(out1[6], ONE);
            outptr[2] = _mm512_extracti32x4_epi32(out1[10], ONE);
            outptr[3] = _mm512_extracti32x4_epi32(out1[14], ONE);

            //will get 8th row of transpose matrix from corresponding 4 vectors in out1
            outptr    = (__m128i *)(out + (j * height + i + (numcol * 7)));
            outptr[0] = _mm512_extracti32x4_epi32(out1[3], ONE);
            outptr[1] = _mm512_extracti32x4_epi32(out1[7], ONE);
            outptr[2] = _mm512_extracti32x4_epi32(out1[11], ONE);
            outptr[3] = _mm512_extracti32x4_epi32(out1[15], ONE);

            //will get 9th row of transpose matrix from corresponding 4 vectors in out1
            outptr    = (__m128i *)(out + (j * height + i + (numcol * 8)));
            outptr[0] = _mm512_extracti32x4_epi32(out1[0], TWO);
            outptr[1] = _mm512_extracti32x4_epi32(out1[4], TWO);
            outptr[2] = _mm512_extracti32x4_epi32(out1[8], TWO);
            outptr[3] = _mm512_extracti32x4_epi32(out1[12], TWO);

            //will get 10th row of transpose matrix from corresponding 4 vectors in out1
            outptr    = (__m128i *)(out + (j * height + i + (numcol * 9)));
            outptr[0] = _mm512_extracti32x4_epi32(out1[1], TWO);
            outptr[1] = _mm512_extracti32x4_epi32(out1[5], TWO);
            outptr[2] = _mm512_extracti32x4_epi32(out1[9], TWO);
            outptr[3] = _mm512_extracti32x4_epi32(out1[13], TWO);

            //will get 11th row of transpose matrix from corresponding 4 vectors in out1
            outptr    = (__m128i *)(out + (j * height + i + (numcol * 10)));
            outptr[0] = _mm512_extracti32x4_epi32(out1[2], TWO);
            outptr[1] = _mm512_extracti32x4_epi32(out1[6], TWO);
            outptr[2] = _mm512_extracti32x4_epi32(out1[10], TWO);
            outptr[3] = _mm512_extracti32x4_epi32(out1[14], TWO);

            //will get 12th row of transpose matrix from corresponding 4 vectors in out1
            outptr    = (__m128i *)(out + (j * height + i + (numcol * 11)));
            outptr[0] = _mm512_extracti32x4_epi32(out1[3], TWO);
            outptr[1] = _mm512_extracti32x4_epi32(out1[7], TWO);
            outptr[2] = _mm512_extracti32x4_epi32(out1[11], TWO);
            outptr[3] = _mm512_extracti32x4_epi32(out1[15], TWO);

            //will get 13th row of transpose matrix from corresponding 4 vectors in out1
            outptr    = (__m128i *)(out + (j * height + i + (numcol * 12)));
            outptr[0] = _mm512_extracti32x4_epi32(out1[0], THREE);
            outptr[1] = _mm512_extracti32x4_epi32(out1[4], THREE);
            outptr[2] = _mm512_extracti32x4_epi32(out1[8], THREE);
            outptr[3] = _mm512_extracti32x4_epi32(out1[12], THREE);

            //will get 14th row of transpose matrix from corresponding 4 vectors in out1
            outptr    = (__m128i *)(out + (j * height + i + (numcol * 13)));
            outptr[0] = _mm512_extracti32x4_epi32(out1[1], THREE);
            outptr[1] = _mm512_extracti32x4_epi32(out1[5], THREE);
            outptr[2] = _mm512_extracti32x4_epi32(out1[9], THREE);
            outptr[3] = _mm512_extracti32x4_epi32(out1[13], THREE);

            //will get 15th row of transpose matrix from corresponding 4 vectors in out1
            outptr    = (__m128i *)(out + (j * height + i + (numcol * 14)));
            outptr[0] = _mm512_extracti32x4_epi32(out1[2], THREE);
            outptr[1] = _mm512_extracti32x4_epi32(out1[6], THREE);
            outptr[2] = _mm512_extracti32x4_epi32(out1[10], THREE);
            outptr[3] = _mm512_extracti32x4_epi32(out1[14], THREE);

            //will get 16th row of transpose matrix from corresponding 4 vectors in out1
            outptr    = (__m128i *)(out + (j * height + i + (numcol * 15)));
            outptr[0] = _mm512_extracti32x4_epi32(out1[3], THREE);
            outptr[1] = _mm512_extracti32x4_epi32(out1[7], THREE);
            outptr[2] = _mm512_extracti32x4_epi32(out1[11], THREE);
            outptr[3] = _mm512_extracti32x4_epi32(out1[15], THREE);
        }
    }
}

#endif /*OPT_AVX512*/
#endif // AOM_DSP_X86_TRANSPOSE_AVX512_H_
