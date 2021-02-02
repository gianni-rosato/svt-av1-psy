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

#ifndef AOM_DSP_X86_TRANSPOSE_ENCODER_AVX512_H_
#define AOM_DSP_X86_TRANSPOSE_ENCODER_AVX512_H_

#include "transpose_avx512.h"
#if OPT_AVX512

static INLINE void transpose_16nx16m_avx512(const __m512i *in, __m512i *out, const int32_t width,
                                            const int32_t height) {
    const int32_t numcol = height >> 4;
    const int32_t numrow = width >> 4;
    __m512i       out1[16];
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

static AOM_FORCE_INLINE void transpose_16nx16n_N2_half_avx512(int32_t        txfm_size,
                                                              const __m512i *input,
                                                              __m512i *      output) {
    const int32_t num_per_512 = 16;
    const int32_t row_size    = txfm_size;
    const int32_t col_size    = txfm_size / num_per_512;
    int32_t       r, c;

    int32_t calc_numrow = row_size >> 1;
    if (!calc_numrow) {
        calc_numrow = 1;
    }

    // transpose each 16x16 block internally
    for (r = 0; r < calc_numrow; r += 16) {
        for (c = 0; c < col_size; c++) {
            transpose_16x16_avx512(
                col_size, &input[r * col_size + c], &output[c * 16 * col_size + r / 16]);
        }
    }
}

static AOM_FORCE_INLINE void transpose_16nx16n_N2_quad_avx512(int32_t        txfm_size,
                                                              const __m512i *input,
                                                              __m512i *      output) {
    const int32_t num_per_512 = 16;
    const int32_t row_size    = txfm_size;
    const int32_t col_size    = txfm_size / num_per_512;
    int32_t       r, c;

    int32_t calc_numrow = row_size >> 1;
    if (!calc_numrow) {
        calc_numrow = 1;
    }
    int32_t calc_numcol = col_size >> 1;
    if (!calc_numcol) {
        calc_numcol = 1;
    }

    // transpose each 16x16 block internally
    for (r = 0; r < calc_numrow; r += 16) {
        for (c = 0; c < calc_numcol; c++) {
            transpose_16x16_avx512(
                col_size, &input[r * col_size + c], &output[c * 16 * col_size + r / 16]);
        }
    }
}

static AOM_FORCE_INLINE void transpose_16nx16m_N2_half_avx512(const __m512i *in, __m512i *out,
                                                              const int32_t width,
                                                              const int32_t height) {
    const int32_t numcol = height >> 4;
    const int32_t numrow = width >> 4;
    __m512i       out1[16];

    int32_t calc_numcol = numcol >> 1;
    if (!calc_numcol) {
        calc_numcol = 1;
    }

    for (int32_t j = 0; j < numrow; j++) {
        for (int32_t i = 0; i < calc_numcol; i++) {
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

static AOM_FORCE_INLINE void transpose_16nx16m_N2_quad_avx512(const __m512i *in, __m512i *out,
                                                              const int32_t width,
                                                              const int32_t height) {
    const int32_t numcol = height >> 4;
    const int32_t numrow = width >> 4;

    int32_t calc_numcol = numcol >> 1;
    if (!calc_numcol) {
        calc_numcol = 1;
    }
    int32_t calc_numrow = numrow >> 1;
    if (!calc_numrow) {
        calc_numrow = 1;
    }

    __m512i out1[16];
    for (int32_t j = 0; j < calc_numrow; j++) {
        for (int32_t i = 0; i < calc_numcol; i++) {
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

static AOM_FORCE_INLINE void transpose_16nx16n_N4_half_avx512(int32_t        txfm_size,
                                                              const __m512i *input,
                                                              __m512i *      output) {
    const int32_t num_per_512 = 16;
    const int32_t row_size    = txfm_size;
    const int32_t col_size    = txfm_size / num_per_512;
    int32_t       r, c;

    int32_t calc_numrow = row_size >> 2;
    if (!calc_numrow) {
        calc_numrow = 1;
    }

    // transpose each 16x16 block internally
    for (r = 0; r < calc_numrow; r += 16) {
        for (c = 0; c < col_size; c++) {
            transpose_16x16_avx512(
                col_size, &input[r * col_size + c], &output[c * 16 * col_size + r / 16]);
        }
    }
}

static AOM_FORCE_INLINE void transpose_16nx16n_N4_quad_avx512(int32_t        txfm_size,
                                                              const __m512i *input,
                                                              __m512i *      output) {
    const int32_t num_per_512 = 16;
    const int32_t row_size    = txfm_size;
    const int32_t col_size    = txfm_size / num_per_512;
    int32_t       r, c;

    int32_t calc_numrow = row_size >> 2;
    if (!calc_numrow) {
        calc_numrow = 1;
    }
    int32_t calc_numcol = col_size >> 2;
    if (!calc_numcol) {
        calc_numcol = 1;
    }

    // transpose each 16x16 block internally
    for (r = 0; r < calc_numrow; r += 16) {
        for (c = 0; c < calc_numcol; c++) {
            transpose_16x16_avx512(
                col_size, &input[r * col_size + c], &output[c * 16 * col_size + r / 16]);
        }
    }
}
#endif /*OPT_AVX512*/
#endif // AOM_DSP_X86_TRANSPOSE_ENCODER_AVX512_H_
