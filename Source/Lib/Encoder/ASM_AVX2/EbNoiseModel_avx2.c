/*
* Copyright(c) 2022 Intel Corporation
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#include <immintrin.h>
#include "EbDefinitions.h"
#include "common_dsp_rtcd.h"

#if FG_LOSSLES_OPT
void svt_av1_add_block_observations_internal_avx2(uint32_t n, const double val,
                                                  const double recp_sqr_norm, double *buffer,
                                                  double *buffer_norm, double *b, double *A) {
    uint32_t      i;
    const __m256d recp_spr_pd = _mm256_set1_pd(recp_sqr_norm);
    const __m256d val_pd      = _mm256_set1_pd(val);
    __m256d       buffer_pd, tmp_pd;

    for (i = 0; i + 8 - 1 < n; i += 8) {
        buffer_pd = _mm256_loadu_pd(buffer + i);
        buffer_pd = _mm256_mul_pd(buffer_pd, recp_spr_pd);
        _mm256_storeu_pd(buffer_norm + i, buffer_pd);
        buffer_pd = _mm256_mul_pd(buffer_pd, val_pd);
        tmp_pd    = _mm256_loadu_pd(b + i);
        tmp_pd    = _mm256_add_pd(tmp_pd, buffer_pd);
        _mm256_storeu_pd(b + i, tmp_pd);

        buffer_pd = _mm256_loadu_pd(buffer + i + 4);
        buffer_pd = _mm256_mul_pd(buffer_pd, recp_spr_pd);
        _mm256_storeu_pd(buffer_norm + i + 4, buffer_pd);
        buffer_pd = _mm256_mul_pd(buffer_pd, val_pd);
        tmp_pd    = _mm256_loadu_pd(b + i + 4);
        tmp_pd    = _mm256_add_pd(tmp_pd, buffer_pd);
        _mm256_storeu_pd(b + i + 4, tmp_pd);
    }
    for (; i < n; ++i) {
        buffer_norm[i] = buffer[i] * recp_sqr_norm;
        b[i] += buffer_norm[i] * val;
    }

    for (i = 0; i < n; ++i) {
        uint32_t      j              = 0;
        const double  buffer_norm_i  = buffer_norm[i];
        const __m256d buffer_norm_pd = _mm256_set1_pd(buffer_norm_i);

        for (j = 0; j + 8 - 1 < n; j += 8) {
            buffer_pd = _mm256_loadu_pd(buffer + j);
            tmp_pd    = _mm256_loadu_pd(A + i * n + j);
            buffer_pd = _mm256_mul_pd(buffer_pd, buffer_norm_pd);
            tmp_pd    = _mm256_add_pd(tmp_pd, buffer_pd);
            _mm256_storeu_pd(A + i * n + j, tmp_pd);

            buffer_pd = _mm256_loadu_pd(buffer + j + 4);
            tmp_pd    = _mm256_loadu_pd(A + i * n + j + 4);
            buffer_pd = _mm256_mul_pd(buffer_pd, buffer_norm_pd);
            tmp_pd    = _mm256_add_pd(tmp_pd, buffer_pd);
            _mm256_storeu_pd(A + i * n + j + 4, tmp_pd);
        }
        for (; j < n; ++j) { A[i * n + j] += (buffer_norm_i * buffer[j]); }
    }
}

void svt_av1_pointwise_multiply_avx2(const float *a, float *b, float *c, double *b_d, double *c_d,
                                     int32_t n) {
    int32_t i = 0;
    __m256   a_ps, b_ps, c_ps;
    __m128   tmp_ps1, tmp_ps2;

    for (; i + 8 - 1 < n; i += 8) {
        a_ps = _mm256_loadu_ps(a + i);

        tmp_ps1 = _mm256_cvtpd_ps(_mm256_loadu_pd(b_d + i));
        tmp_ps2 = _mm256_cvtpd_ps(_mm256_loadu_pd(b_d + i + 4));
        b_ps    = _mm256_insertf128_ps(_mm256_castps128_ps256(tmp_ps1), tmp_ps2, 0x1);
        tmp_ps1 = _mm256_cvtpd_ps(_mm256_loadu_pd(c_d + i));
        tmp_ps2 = _mm256_cvtpd_ps(_mm256_loadu_pd(c_d + i + 4));
        c_ps    = _mm256_insertf128_ps(_mm256_castps128_ps256(tmp_ps1), tmp_ps2, 0x1);

        _mm256_storeu_ps(b + i, _mm256_mul_ps(a_ps, b_ps));
        _mm256_storeu_ps(c + i, _mm256_mul_ps(a_ps, c_ps));
    }
    for (; i < n; i++) {
        b[i] = a[i] * (float)b_d[i];
        c[i] = a[i] * (float)c_d[i];
    }
}

void svt_av1_apply_window_function_to_plane_avx2(int32_t y_size, int32_t x_size, float *result_ptr,
                                                 uint32_t result_stride, float *block, float *plane,
                                                 const float *window_function) {
    __m256 block_ps, plane_ps, wnd_funct_ps, res_ps;
    for (int32_t y = 0; y < y_size; ++y) {
        int32_t x = 0;
        for (; x + 8 - 1 < x_size; x += 8) {
            block_ps     = _mm256_loadu_ps(block + y * x_size + x);
            plane_ps     = _mm256_loadu_ps(plane + y * x_size + x);
            wnd_funct_ps = _mm256_loadu_ps(window_function + y * x_size + x);
            res_ps       = _mm256_loadu_ps(result_ptr + y * result_stride + x);

            block_ps = _mm256_add_ps(block_ps, plane_ps);
            block_ps = _mm256_mul_ps(block_ps, wnd_funct_ps);
            block_ps = _mm256_add_ps(block_ps, res_ps);
            _mm256_storeu_ps(result_ptr + y * result_stride + x, block_ps);
        }
        for (; x < x_size; ++x) {
            result_ptr[y * result_stride + x] += (block[y * x_size + x] + plane[y * x_size + x]) *
                window_function[y * x_size + x];
        }
    }
}

#endif
