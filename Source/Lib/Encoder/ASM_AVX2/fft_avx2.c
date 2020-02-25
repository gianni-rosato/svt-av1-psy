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

#include <immintrin.h>

#include "EbDefinitions.h"
#include "fft_common.h"

extern void eb_aom_transpose_float_sse2(const float *A, float *b, int32_t n);
extern void eb_aom_fft_unpack_2d_output_sse2(const float *col_fft, float *output, int32_t n);

// Generate the 1d forward transforms for float using _mm256
GEN_FFT_8(static INLINE void, avx2, float, __m256, _mm256_loadu_ps, _mm256_storeu_ps, _mm256_set1_ps,
          _mm256_add_ps, _mm256_sub_ps, _mm256_mul_ps);
GEN_FFT_16(static INLINE void, avx2, float, __m256, _mm256_loadu_ps, _mm256_storeu_ps, _mm256_set1_ps,
           _mm256_add_ps, _mm256_sub_ps, _mm256_mul_ps);
GEN_FFT_32(static INLINE void, avx2, float, __m256, _mm256_loadu_ps, _mm256_storeu_ps, _mm256_set1_ps,
           _mm256_add_ps, _mm256_sub_ps, _mm256_mul_ps);

void eb_aom_fft8x8_float_avx2(const float *input, float *temp, float *output) {
    eb_aom_fft_2d_gen(input,
                      temp,
                      output,
                      8,
                      eb_aom_fft1d_8_avx2,
                      eb_aom_transpose_float_sse2,
                      eb_aom_fft_unpack_2d_output_sse2,
                      8);
}

void eb_aom_fft16x16_float_avx2(const float *input, float *temp, float *output) {
    eb_aom_fft_2d_gen(input,
                      temp,
                      output,
                      16,
                      eb_aom_fft1d_16_avx2,
                      eb_aom_transpose_float_sse2,
                      eb_aom_fft_unpack_2d_output_sse2,
                      8);
}

void eb_aom_fft32x32_float_avx2(const float *input, float *temp, float *output) {
    eb_aom_fft_2d_gen(input,
                      temp,
                      output,
                      32,
                      eb_aom_fft1d_32_avx2,
                      eb_aom_transpose_float_sse2,
                      eb_aom_fft_unpack_2d_output_sse2,
                      8);
}

// Generate the 1d inverse transforms for float using _mm256
GEN_IFFT_8(static INLINE void, avx2, float, __m256, _mm256_loadu_ps, _mm256_storeu_ps, _mm256_set1_ps,
           _mm256_add_ps, _mm256_sub_ps, _mm256_mul_ps);
GEN_IFFT_16(static INLINE void, avx2, float, __m256, _mm256_loadu_ps, _mm256_storeu_ps,
            _mm256_set1_ps, _mm256_add_ps, _mm256_sub_ps, _mm256_mul_ps);
GEN_IFFT_32(static INLINE void, avx2, float, __m256, _mm256_loadu_ps, _mm256_storeu_ps,
            _mm256_set1_ps, _mm256_add_ps, _mm256_sub_ps, _mm256_mul_ps);

void eb_aom_ifft8x8_float_avx2(const float *input, float *temp, float *output) {
    eb_aom_ifft_2d_gen(input,
                       temp,
                       output,
                       8,
                       eb_aom_fft1d_8_float,
                       eb_aom_fft1d_8_avx2,
                       eb_aom_ifft1d_8_avx2,
                       eb_aom_transpose_float_sse2,
                       8);
}

void eb_aom_ifft16x16_float_avx2(const float *input, float *temp, float *output) {
    eb_aom_ifft_2d_gen(input,
                       temp,
                       output,
                       16,
                       eb_aom_fft1d_16_float,
                       eb_aom_fft1d_16_avx2,
                       eb_aom_ifft1d_16_avx2,
                       eb_aom_transpose_float_sse2,
                       8);
}

void eb_aom_ifft32x32_float_avx2(const float *input, float *temp, float *output) {
    eb_aom_ifft_2d_gen(input,
                       temp,
                       output,
                       32,
                       eb_aom_fft1d_32_float,
                       eb_aom_fft1d_32_avx2,
                       eb_aom_ifft1d_32_avx2,
                       eb_aom_transpose_float_sse2,
                       8);
}
