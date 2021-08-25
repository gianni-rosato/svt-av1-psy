/*
 * Copyright (c) 2019, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
 */

#include <assert.h>
#include <immintrin.h> /* AVX2 */

#include "EbDefinitions.h"
#include "./../ASM_SSE4_1/EbTemporalFiltering_constants.h"
#if FIXED_POINTS_PLANEWISE
#include "EbUtility.h"
#endif

#if FIXED_POINTS_PLANEWISE
/*value [i:0-15] (sqrt((float)i)*65536.0*/
static const uint32_t sqrt_array_fp16[16] = {
    0,
    65536,
    92681,
    113511,
    131072,
    146542,
    160529,
    173391,
    185363,
    196608,
    207243,
    217358,
    227023,
    236293,
    245213,
    253819
};

/*Calc sqrt linear max error 10%*/
static uint32_t sqrt_fast(uint32_t x)
{
    if (x > 15) {
        const int log2_half = svt_log2f(x)>>1;
        const int mul2 = log2_half<<1;
        int base = x >> (mul2 - 2);
        assert(base < 16);
        return sqrt_array_fp16[base] >> (17 - log2_half);
    }
    return sqrt_array_fp16[x] >> 16;
}

// T[X] =  exp(-(X)/16)  for x in [0..7], step 1/16 values in Fixed Points shift 16
static const int32_t expf_tab_fp16[] = {
65536,
61565,
57835,
54331,
51039,
47947,
45042,
42313,
39749,
37341,
35078,
32953,
30957,
29081,
27319,
25664,
24109,
22648,
21276,
19987,
18776,
17638,
16570,
15566,
14623,
13737,
12904,
12122,
11388,
10698,
10050,
9441,
8869,
8331,
7827,
7352,
6907,
6488,
6095,
5726,
5379,
5053,
4747,
4459,
4189,
3935,
3697,
3473,
3262,
3065,
2879,
2704,
2541,
2387,
2242,
2106,
1979,
1859,
1746,
1640,
1541,
1447,
1360,
1277,
1200,
1127,
1059,
995,
934,
878,
824,
774,
728,
683,
642,
603,
566,
532,
500,
470,
441,
414,
389,
366,
343,
323,
303,
285,
267,
251,
236,
222,
208,
195,
184,
172,
162,
152,
143,
134,
126,
118,
111,
104,
98,
92,
86,
81,
76,
72,
67,
63,
59,
56,
52,
49,
46,
43,
41,
38,
36,
34,
31,
30,
28,
26,
24,
23,
21
};
#endif /*FIXED_POINTS_PLANEWISE*/

#define SSE_STRIDE (BW + 2)

DECLARE_ALIGNED(32, static const uint32_t, sse_bytemask[4][8]) = {
    {0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0, 0, 0},
    {0, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0, 0},
    {0, 0, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0},
    {0, 0, 0, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF}};

DECLARE_ALIGNED(32, static const uint8_t, shufflemask_16b[2][16]) = {
    {0, 1, 0, 1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11},
    {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 10, 11, 10, 11}};

DECLARE_ALIGNED(32, static const uint8_t, shufflemask_32b[2][16]) = {
    {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 6, 7},
    {0, 1, 2, 3, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7}};

static AOM_FORCE_INLINE void get_squared_error_16x16_avx2(
    const uint8_t *frame1, const unsigned int stride, const uint8_t *frame2,
    const unsigned int stride2, const int block_width, const int block_height, uint16_t *frame_sse,
    const unsigned int sse_stride) {
    (void)block_width;
    const uint8_t *src1 = frame1;
    const uint8_t *src2 = frame2;
    uint16_t *     dst  = frame_sse;
    for (int i = 0; i < block_height; i++) {
        __m128i vf1_128, vf2_128;
        __m256i vf1, vf2, vdiff1, vsqdiff1;

        vf1_128  = _mm_loadu_si128((__m128i *)(src1));
        vf2_128  = _mm_loadu_si128((__m128i *)(src2));
        vf1      = _mm256_cvtepu8_epi16(vf1_128);
        vf2      = _mm256_cvtepu8_epi16(vf2_128);
        vdiff1   = _mm256_sub_epi16(vf1, vf2);
        vsqdiff1 = _mm256_mullo_epi16(vdiff1, vdiff1);

        _mm256_storeu_si256((__m256i *)(dst), vsqdiff1);
        // Set zero to unitialized memory to avoid uninitialized loads later
        *(uint32_t *)(dst + 16) = _mm_cvtsi128_si32(_mm_setzero_si128());

        src1 += stride, src2 += stride2;
        dst += sse_stride;
    }
}

static AOM_FORCE_INLINE void get_squared_error_32x32_avx2(
    const uint8_t *frame1, const unsigned int stride, const uint8_t *frame2,
    const unsigned int stride2, const int block_width, const int block_height, uint16_t *frame_sse,
    const unsigned int sse_stride) {
    (void)block_width;
    const uint8_t *src1 = frame1;
    const uint8_t *src2 = frame2;
    uint16_t *     dst  = frame_sse;
    for (int i = 0; i < block_height; i++) {
        __m256i vsrc1, vsrc2, vmin, vmax, vdiff, vdiff1, vdiff2, vres1, vres2;

        vsrc1 = _mm256_loadu_si256((__m256i *)src1);
        vsrc2 = _mm256_loadu_si256((__m256i *)src2);
        vmax  = _mm256_max_epu8(vsrc1, vsrc2);
        vmin  = _mm256_min_epu8(vsrc1, vsrc2);
        vdiff = _mm256_subs_epu8(vmax, vmin);

        __m128i vtmp1 = _mm256_castsi256_si128(vdiff);
        __m128i vtmp2 = _mm256_extracti128_si256(vdiff, 1);
        vdiff1        = _mm256_cvtepu8_epi16(vtmp1);
        vdiff2        = _mm256_cvtepu8_epi16(vtmp2);

        vres1 = _mm256_mullo_epi16(vdiff1, vdiff1);
        vres2 = _mm256_mullo_epi16(vdiff2, vdiff2);
        _mm256_storeu_si256((__m256i *)(dst), vres1);
        _mm256_storeu_si256((__m256i *)(dst + 16), vres2);
        // Set zero to unitialized memory to avoid uninitialized loads later
        *(uint32_t *)(dst + 32) = _mm_cvtsi128_si32(_mm_setzero_si128());

        src1 += stride;
        src2 += stride2;
        dst += sse_stride;
    }
}

static AOM_FORCE_INLINE __m256i xx_load_and_pad(uint16_t *src, int col, int block_width) {
    __m128i v128tmp = _mm_loadu_si128((__m128i *)(src));
    if (col == 0) {
        // For the first column, replicate the first element twice to the left
        v128tmp = _mm_shuffle_epi8(v128tmp, *(__m128i *)shufflemask_16b[0]);
    }
    if (col == block_width - 4) {
        // For the last column, replicate the last element twice to the right
        v128tmp = _mm_shuffle_epi8(v128tmp, *(__m128i *)shufflemask_16b[1]);
    }
    return _mm256_cvtepu16_epi32(v128tmp);
}

static AOM_FORCE_INLINE int32_t xx_mask_and_hadd(__m256i vsum, int i) {
    // Mask the required 5 values inside the vector
    __m256i vtmp = _mm256_and_si256(vsum, *(__m256i *)sse_bytemask[i]);
    __m128i v128a, v128b;
    // Extract 256b as two 128b registers A and B
    v128a = _mm256_castsi256_si128(vtmp);
    v128b = _mm256_extracti128_si256(vtmp, 1);
    // A = [A0+B0, A1+B1, A2+B2, A3+B3]
    v128a = _mm_add_epi32(v128a, v128b);
    // B = [A2+B2, A3+B3, 0, 0]
    v128b = _mm_srli_si128(v128a, 8);
    // A = [A0+B0+A2+B2, A1+B1+A3+B3, X, X]
    v128a = _mm_add_epi32(v128a, v128b);
    // B = [A1+B1+A3+B3, 0, 0, 0]
    v128b = _mm_srli_si128(v128a, 4);
    // A = [A0+B0+A2+B2+A1+B1+A3+B3, X, X, X]
    v128a = _mm_add_epi32(v128a, v128b);
    return _mm_extract_epi32(v128a, 0);
}
#if SIMD_APPROX_EXPF
static AOM_FORCE_INLINE __m256 exp_256_ps(__m256 _x) {

    __m256 t, f, p, r;
    __m256i i, j;

    const __m256 l2e = _mm256_set1_ps(1.442695041f); /* log2(e) */
    const __m256 l2h = _mm256_set1_ps(-6.93145752e-1f); /* -log(2)_hi */
    const __m256 l2l = _mm256_set1_ps(-1.42860677e-6f); /* -log(2)_lo */
    /* coefficients for core approximation to exp() in [-log(2)/2, log(2)/2] */
    const __m256 c0 = _mm256_set1_ps(0.008301110f);
    const __m256 c1 = _mm256_set1_ps(0.041906696f);
    const __m256 c2 = _mm256_set1_ps(0.166674897f);
    const __m256 c3 = _mm256_set1_ps(0.499990642f);
    const __m256 c4 = _mm256_set1_ps(0.999999762f);
    const __m256 c5 = _mm256_set1_ps(1.000000000f);

    /* exp(x) = 2^i * e^f; i = rint (log2(e) * x), f = x - log(2) * i */
    t = _mm256_mul_ps(_x, l2e);      /* t = log2(e) * x */
    r = _mm256_round_ps(t, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC); /* r = rint (t) */

    p = _mm256_mul_ps(r, l2h);      /* log(2)_hi * r */
    f = _mm256_add_ps(_x, p);        /* x - log(2)_hi * r */
    p = _mm256_mul_ps(r, l2l);      /* log(2)_lo * r */
    f = _mm256_add_ps(f, p);        /* f = x - log(2)_hi * r - log(2)_lo * r */

    i = _mm256_cvtps_epi32(t);       /* i = (int)rint(t) */

    /* p ~= exp (f), -log(2)/2 <= f <= log(2)/2 */
    p = c0;                          /* c0 */

    p = _mm256_mul_ps(p, f);        /* c0*f */
    p = _mm256_add_ps(p, c1);       /* c0*f+c1 */
    p = _mm256_mul_ps(p, f);        /* (c0*f+c1)*f */
    p = _mm256_add_ps(p, c2);       /* (c0*f+c1)*f+c2 */
    p = _mm256_mul_ps(p, f);        /* ((c0*f+c1)*f+c2)*f */
    p = _mm256_add_ps(p, c3);       /* ((c0*f+c1)*f+c2)*f+c3 */
    p = _mm256_mul_ps(p, f);        /* (((c0*f+c1)*f+c2)*f+c3)*f */
    p = _mm256_add_ps(p, c4);       /* (((c0*f+c1)*f+c2)*f+c3)*f+c4 ~= exp(f) */
    p = _mm256_mul_ps(p, f);        /* (((c0*f+c1)*f+c2)*f+c3)*f */
    p = _mm256_add_ps(p, c5);       /* (((c0*f+c1)*f+c2)*f+c3)*f+c4 ~= exp(f) */

    /* exp(x) = 2^i * p */
    j = _mm256_slli_epi32(i, 23); /* i << 23 */
    r = _mm256_castsi256_ps(_mm256_add_epi32(j, _mm256_castps_si256(p))); /* r = p * 2^i */

    return r;


}
#else /*SIMD_APPROX_EXPF*/
static AOM_FORCE_INLINE __m256 exp_256_ps(__m256 _x) {
    float buff[8];
    _mm256_storeu_ps(buff, _x);

    buff[0] = expf(buff[0]);
    buff[1] = expf(buff[1]);
    buff[2] = expf(buff[2]);
    buff[3] = expf(buff[3]);
    buff[4] = expf(buff[4]);
    buff[5] = expf(buff[5]);
    buff[6] = expf(buff[6]);
    buff[7] = expf(buff[7]);

    __m256 res = _mm256_loadu_ps(buff);
    return res;
}
#endif /*SIMD_APPROX_EXPF*/

static void apply_temporal_filter_planewise(struct MeContext *context_ptr, const uint8_t *frame1,
                                            const unsigned int stride, const uint8_t *frame2,
                                            const unsigned int stride2, const int block_width,
                                            const int block_height,
#if !FTR_TF_STRENGTH_PER_QP
                                            const double sigma,
                                            const int decay_control,
#endif
                                            unsigned int *accumulator,
                                            uint16_t *count, uint16_t *luma_sq_error,
                                            uint16_t *chroma_sq_error, int plane, int ss_x_shift,
                                            int ss_y_shift) {
    assert(TF_PLANEWISE_FILTER_WINDOW_LENGTH == 5);
    assert(((block_width == 32) && (block_height == 32)) ||
           ((block_width == 16) && (block_height == 16)));
    if (plane > PLANE_TYPE_Y)
        assert(chroma_sq_error != NULL);

    uint32_t acc_5x5_sse[BH][BW];
    // Larger noise -> larger filtering weight.

#if !FIXED_POINTS_PLANEWISE
#if FTR_TF_STRENGTH_PER_QP
    double n_decay_qr_inv = 1.0 / context_ptr->tf_decay_factor[plane];
#else
    const double n_decay                = (double)decay_control * (0.7 + log1p(sigma));
    const double n_decay_qr_inv         = 1.0 / (2 * n_decay * n_decay);
#endif
    const double block_balacne_inv      = 1.0 / (TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1);
    const double distance_threshold_inv = 1.0 /
        (double)AOMMAX(context_ptr->min_frame_size * TF_SEARCH_DISTANCE_THRESHOLD, 1);
#endif /*!FIXED_POINTS_PLANEWISE*/

    uint16_t *frame_sse = (plane == PLANE_TYPE_Y) ? luma_sq_error : chroma_sq_error;

    if (block_width == 32) {
        get_squared_error_32x32_avx2(
            frame1, stride, frame2, stride2, block_width, block_height, frame_sse, SSE_STRIDE);
    } else {
        get_squared_error_16x16_avx2(
            frame1, stride, frame2, stride2, block_width, block_height, frame_sse, SSE_STRIDE);
    }

    __m256i vsrc[5];

    // Traverse 4 columns at a time
    // First and last columns will require padding
    for (int col = 0; col < block_width; col += 4) {
        uint16_t *src = (col) ? frame_sse + col - 2 : frame_sse;

        // Load and pad(for first and last col) 3 rows from the top
        for (int i = 2; i < 5; i++) {
            vsrc[i] = xx_load_and_pad(src, col, block_width);
            src += SSE_STRIDE;
        }

        // Copy first row to first 2 vectors
        vsrc[0] = vsrc[2];
        vsrc[1] = vsrc[2];

        for (int row = 0; row < block_height; row++) {
            __m256i vsum = _mm256_setzero_si256();

            // Add 5 consecutive rows
            for (int i = 0; i < 5; i++) { vsum = _mm256_add_epi32(vsum, vsrc[i]); }

            // Push all elements by one element to the top
            for (int i = 0; i < 4; i++) { vsrc[i] = vsrc[i + 1]; }

            // Load next row to the last element
            if (row <= block_width - 4) {
                vsrc[4] = xx_load_and_pad(src, col, block_width);
                src += SSE_STRIDE;
            } else {
                vsrc[4] = vsrc[3];
            }

            // Accumulate the sum horizontally
            for (int i = 0; i < 4; i++) { acc_5x5_sse[row][col + i] = xx_mask_and_hadd(vsum, i); }
        }
    }

    int32_t idx_32x32      = context_ptr->tf_block_col + context_ptr->tf_block_row * 2;
#if FIXED_POINTS_PLANEWISE
    uint32_t num_ref_pixels = TF_PLANEWISE_FILTER_WINDOW_LENGTH * TF_PLANEWISE_FILTER_WINDOW_LENGTH;
#else
    int num_ref_pixels = TF_PLANEWISE_FILTER_WINDOW_LENGTH * TF_PLANEWISE_FILTER_WINDOW_LENGTH;
#endif

    if (plane != PLANE_TYPE_Y) {
        num_ref_pixels += (1 << ss_y_shift) * (1 << ss_x_shift);

#if FIXED_POINTS_PLANEWISE
        if (context_ptr->tf_ctrls.use_fixed_point && (ss_x_shift == 1) && (ss_y_shift == 1)) {
            for (int i = 0; i < block_height; i++) {
                int yy = (i << 1);
                for (int j = 0; j < block_width; j += 4) {
                    int     xx   = (j << 1);
                    __m256i line = _mm256_cvtepu16_epi32(
                        _mm_loadu_si128((__m128i *)(luma_sq_error + yy * SSE_STRIDE + xx)));
                    __m256i line2 = _mm256_cvtepu16_epi32(
                        _mm_loadu_si128((__m128i *)(luma_sq_error + (yy + 1) * SSE_STRIDE + xx)));
                    __m256i sum = _mm256_hadd_epi32(line, line2);
                    sum         = _mm256_permute4x64_epi64(sum, 0xd8);

                    __m128i sum_128 = _mm_add_epi32(_mm256_castsi256_si128(sum),
                                                    _mm256_extracti128_si256(sum, 1));

                    __m128i diff_sse = _mm_loadu_si128((__m128i *)(acc_5x5_sse[i] + j));
                    diff_sse         = _mm_add_epi32(diff_sse, sum_128);
                    _mm_storeu_si128((__m128i *)(acc_5x5_sse[i] + j), diff_sse);
                }
            }
        } else {
            for (int i = 0; i < block_height; i++) {
                for (int j = 0; j < block_width; j++) {
                    for (int ii = 0; ii < (1 << ss_y_shift); ++ii) {
                        for (int jj = 0; jj < (1 << ss_x_shift); ++jj) {
                            const int yy = (i << ss_y_shift) + ii; // Y-coord on Y-plane.
                            const int xx = (j << ss_x_shift) + jj; // X-coord on Y-plane.
                            acc_5x5_sse[i][j] += luma_sq_error[yy * SSE_STRIDE + xx];
                        }
                    }
                }
            }
        }
#else
        for (int i = 0; i < block_height; i++) {
            for (int j = 0; j < block_width; j++) {
                for (int ii = 0; ii < (1 << ss_y_shift); ++ii) {
                    for (int jj = 0; jj < (1 << ss_x_shift); ++jj) {
                        const int yy = (i << ss_y_shift) + ii; // Y-coord on Y-plane.
                        const int xx = (j << ss_x_shift) + jj; // X-coord on Y-plane.
                        acc_5x5_sse[i][j] += luma_sq_error[yy * SSE_STRIDE + xx];
                    }
                }
            }
        }
#endif
    }

    assert(!(block_width & 7));
#if FIXED_POINTS_PLANEWISE
    __m256i scaled_diff_B_fp4_arr[4], scaled_diff_A1_factor_fp14_arr[4], scaled_diff_A2_factor_fp10_arr[4];
    __m256i seven_fp   = _mm256_set1_epi32(7 * 16);
    __m256i tf_weight_scale_fp   = _mm256_set1_epi32(TF_WEIGHT_SCALE);
    int32_t move_arr[4];

    __m256d num_ref_pix            = _mm256_set1_pd((double)num_ref_pixels);
    __m256d blk_balacne_inv        = _mm256_set1_pd(1.0 / (TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1));
    __m256d wnd_blk_balacne_weight = _mm256_set1_pd((double)TF_WINDOW_BLOCK_BALANCE_WEIGHT);
    __m256  seven                  = _mm256_set1_ps(7.0f);
    __m256  zero                   = _mm256_set1_ps(0.0f);
    __m256  tf_weight_scale        = _mm256_set1_ps((float)TF_WEIGHT_SCALE);
    __m256d blk_errors[4];
    __m256d d_factor_mul_n_decay_qr_invs[4];
#else /*FIXED_POINTS_PLANEWISE*/
    __m256d num_ref_pix            = _mm256_set1_pd((double)num_ref_pixels);
    __m256d blk_balacne_inv        = _mm256_set1_pd(block_balacne_inv);
    __m256d wnd_blk_balacne_weight = _mm256_set1_pd((double)TF_WINDOW_BLOCK_BALANCE_WEIGHT);
    __m256  seven                  = _mm256_set1_ps(7.0f);
    __m256  zero                   = _mm256_set1_ps(0.0f);
    __m256  tf_weight_scale        = _mm256_set1_ps((float)TF_WEIGHT_SCALE);
    __m256d blk_errors[4];
    __m256d d_factor_mul_n_decay_qr_invs[4];
#endif /*FIXED_POINTS_PLANEWISE*/

    if (context_ptr->tf_32x32_block_split_flag[idx_32x32]) {
        for (int i = 0; i < 4; i++) {
#if FIXED_POINTS_PLANEWISE
            if (context_ptr->tf_ctrls.use_fixed_point) {
                uint32_t block_error_fp8 = (uint32_t)
                                              context_ptr->tf_16x16_block_error[idx_32x32 * 4 + i];
                int32_t  col                     = context_ptr->tf_16x16_mv_x[idx_32x32 * 4 + i];
                int32_t  row                     = context_ptr->tf_16x16_mv_y[idx_32x32 * 4 + i];
                uint32_t distance_fp4            = sqrt_fast(((uint32_t)(col * col + row * row)) << 8);
                uint32_t  distance_threshold_fp16 = AOMMAX(
                    (context_ptr->min_frame_size << 16) / 10, 1 << 16);

                int32_t d_factor_fp12 = AOMMAX(
                    (int32_t)(((((int64_t)distance_fp4) << 24) + (distance_threshold_fp16 >> 1)) /
                              (distance_threshold_fp16)), 1 << 12);

                int32_t scaled_diff_B_fp4 = (int32_t)(
                    (((int64_t)d_factor_fp12) *
                     ((block_error_fp8 + ((TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1) / 2)) /
                      (TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1))) /
                    ((context_ptr->tf_decay_factor_fp16[plane])));

                int32_t scaled_diff_A1_factor_fp14 = (int32_t)(
                    ((((int64_t)d_factor_fp12) << 16) * TF_WINDOW_BLOCK_BALANCE_WEIGHT +
                     (context_ptr->tf_decay_factor_fp16[plane] >> 3)) /
                    (context_ptr->tf_decay_factor_fp16[plane] >> 2));

                int32_t scaled_diff_A2_factor_fp10 =
                    ((1 << 10) + ((TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1) * num_ref_pixels) / 2) /
                    ((TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1) * num_ref_pixels);

                int move = 0;
                if (scaled_diff_A1_factor_fp14 > (1 << 20)) {
                    move = 18;
                } else if (scaled_diff_A1_factor_fp14 > (1 << 16)) {
                    move = 14;
                } else if (scaled_diff_A1_factor_fp14 > (1 << 12)) {
                    move = 10;
                } else if (scaled_diff_A1_factor_fp14 > (1 << 8)) {
                    move = 6;
                } else if (scaled_diff_A1_factor_fp14 > (1 << 4)) {
                    move = 2;
                } else {
                    move = 0;
                }

                scaled_diff_B_fp4_arr[i]          = _mm256_set1_epi32(scaled_diff_B_fp4);
                scaled_diff_A1_factor_fp14_arr[i] = _mm256_set1_epi32(scaled_diff_A1_factor_fp14);
                scaled_diff_A2_factor_fp10_arr[i] = _mm256_set1_epi32(scaled_diff_A2_factor_fp10);
                move_arr[i]                       = move;
            } else {
                blk_errors[i] = _mm256_set1_pd(
                    (double)context_ptr->tf_16x16_block_error[idx_32x32 * 4 + i] / 256.0);

                int16_t col      = context_ptr->tf_16x16_mv_x[idx_32x32 * 4 + i];
                int16_t row      = context_ptr->tf_16x16_mv_y[idx_32x32 * 4 + i];
                float   distance = sqrtf((float)(row * row + col * col));
                double  d_factor = AOMMAX(distance * (1.0 / (double)AOMMAX(context_ptr->min_frame_size *
                    TF_SEARCH_DISTANCE_THRESHOLD, 1)), 1);
                d_factor_mul_n_decay_qr_invs[i] = _mm256_set1_pd(
                    d_factor * (1.0 / context_ptr->tf_decay_factor[plane]));
            }
#else /*FIXED_POINTS_PLANEWISE*/
            blk_errors[i] = _mm256_set1_pd(
                (double)context_ptr->tf_16x16_block_error[idx_32x32 * 4 + i] / 256.0);

            int16_t col                     = context_ptr->tf_16x16_mv_x[idx_32x32 * 4 + i];
            int16_t row                     = context_ptr->tf_16x16_mv_y[idx_32x32 * 4 + i];
            float   distance                = sqrtf((float)(row * row + col * col));
            double  d_factor                = AOMMAX(distance * distance_threshold_inv, 1);
            d_factor_mul_n_decay_qr_invs[i] = _mm256_set1_pd(d_factor * n_decay_qr_inv);
#endif /*FIXED_POINTS_PLANEWISE*/
        }
    } else {
#if FIXED_POINTS_PLANEWISE
        if (context_ptr->tf_ctrls.use_fixed_point) {
            uint32_t block_error_fp8 = (uint32_t)context_ptr->tf_32x32_block_error[idx_32x32] >> 2;
            int32_t  col             = context_ptr->tf_32x32_mv_x[idx_32x32];
            int32_t  row             = context_ptr->tf_32x32_mv_y[idx_32x32];
            uint32_t distance_fp4    = sqrt_fast(((uint32_t)(col * col + row * row)) << 8);
            uint32_t distance_threshold_fp16 = AOMMAX((context_ptr->min_frame_size << 16) / 10,
                                                     1 << 16);

            int32_t d_factor_fp12 = AOMMAX(
                (int32_t)(((((int64_t)distance_fp4) << 24) + (distance_threshold_fp16 >> 1)) /
                          (distance_threshold_fp16)), 1 << 12);

            int32_t scaled_diff_B_fp4 = (int32_t)(
                (((int64_t)d_factor_fp12) *
                 ((block_error_fp8 + ((TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1) / 2)) /
                  (TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1))) /
                ((context_ptr->tf_decay_factor_fp16[plane])));

            int32_t scaled_diff_A1_factor_fp14 = (int32_t)(
                ((((int64_t)d_factor_fp12) << 16) * TF_WINDOW_BLOCK_BALANCE_WEIGHT +
                 (context_ptr->tf_decay_factor_fp16[plane] >> 3)) /
                (context_ptr->tf_decay_factor_fp16[plane] >> 2));

            int32_t scaled_diff_A2_factor_fp10 =
                ((1 << 10) + ((TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1) * num_ref_pixels) / 2) /
                ((TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1) * num_ref_pixels);

            int move = 0;
            if (scaled_diff_A1_factor_fp14 > (1 << 20)) {
                move = 18;
            } else if (scaled_diff_A1_factor_fp14 > (1 << 16)) {
                move = 14;
            } else if (scaled_diff_A1_factor_fp14 > (1 << 12)) {
                move = 10;
            } else if (scaled_diff_A1_factor_fp14 > (1 << 8)) {
                move = 6;
            } else if (scaled_diff_A1_factor_fp14 > (1 << 4)) {
                move = 2;
            } else {
                move = 0;
            }

            for (int i = 0; i < 4; i++) {
                scaled_diff_B_fp4_arr[i]          = _mm256_set1_epi32(scaled_diff_B_fp4);
                scaled_diff_A1_factor_fp14_arr[i] = _mm256_set1_epi32(scaled_diff_A1_factor_fp14);
                scaled_diff_A2_factor_fp10_arr[i] = _mm256_set1_epi32(scaled_diff_A2_factor_fp10);
                move_arr[i]                       = move;
            }
        } else {
            double      block_error = (double)context_ptr->tf_32x32_block_error[idx_32x32] / 1024.0;
            int16_t     col         = context_ptr->tf_32x32_mv_x[idx_32x32];
            int16_t     row         = context_ptr->tf_32x32_mv_y[idx_32x32];
            const float distance    = sqrtf((float)(row * row + col * col));
            const double d_factor   = AOMMAX(
                distance *
                    (1.0 /
                     (double)AOMMAX(context_ptr->min_frame_size * TF_SEARCH_DISTANCE_THRESHOLD, 1)),
                1);

            const double d_factor_mul_n_decay_qr_inv = d_factor *
                (1.0 / context_ptr->tf_decay_factor[plane]);

            for (int i = 0; i < 4; i++) {
                blk_errors[i]                   = _mm256_set1_pd(block_error);
                d_factor_mul_n_decay_qr_invs[i] = _mm256_set1_pd(d_factor_mul_n_decay_qr_inv);
            }
        }
#else /*FIXED_POINTS_PLANEWISE*/
        double       block_error = (double)context_ptr->tf_32x32_block_error[idx_32x32] / 1024.0;
        int16_t      col         = context_ptr->tf_32x32_mv_x[idx_32x32];
        int16_t      row         = context_ptr->tf_32x32_mv_y[idx_32x32];
        const float  distance    = sqrtf((float)(row * row + col * col));
        const double d_factor    = AOMMAX(distance * distance_threshold_inv, 1);

        const double d_factor_mul_n_decay_qr_inv = d_factor * n_decay_qr_inv;

        for (int i = 0; i < 4; i++) {
            blk_errors[i]                   = _mm256_set1_pd(block_error);
            d_factor_mul_n_decay_qr_invs[i] = _mm256_set1_pd(d_factor_mul_n_decay_qr_inv);
        }
#endif /*FIXED_POINTS_PLANEWISE*/
    }

    for (int i = 0; i < block_height; i++) {
        const int subblock_idx_h = (i >= block_height / 2) * 2;
        for (int j = 0; j < block_width; j += 8) {
#if FIXED_POINTS_PLANEWISE
            __m256i adjusted_weight1;
            if (context_ptr->tf_ctrls.use_fixed_point) {
                const int subblock_idx = subblock_idx_h + (j >= block_width / 2);

                __m256i diff_sse = _mm256_loadu_si256((__m256i *)(acc_5x5_sse[i] + j));

                //int32_t sum_square_part_fp10 = (((int32_t)sum_square_diff) * scaled_diff_A2_factor_fp10);
                __m256i sum_square_part_fp10 = _mm256_mullo_epi32(
                    diff_sse, scaled_diff_A2_factor_fp10_arr[subblock_idx]);

                //int32_t scaled_diff_A_fp4 = (int32_t)(((sum_square_part_fp10) * ((scaled_diff_A1_factor_fp14)>>move) + (1<<(19-move)))>>(20-move));
                __m256i scaled_diff_A1_factor_fp14_avx2 = _mm256_srli_epi32(
                    scaled_diff_A1_factor_fp14_arr[subblock_idx], move_arr[subblock_idx]);

                __m256i scaled_diff_A_fp4 = _mm256_mullo_epi32(sum_square_part_fp10,
                                                               scaled_diff_A1_factor_fp14_avx2);
                scaled_diff_A_fp4         = _mm256_add_epi32(
                    scaled_diff_A_fp4, _mm256_set1_epi32(1 << (19 - move_arr[subblock_idx])));
                scaled_diff_A_fp4 = _mm256_srli_epi32(scaled_diff_A_fp4,
                                                      20 - move_arr[subblock_idx]);

                //int32_t scaled_diff_fp4 = AOMMIN(scaled_diff_A_fp4 + scaled_diff_B_fp4, 7*16);
                __m256i scaled_diff_fp4 = _mm256_add_epi32(scaled_diff_A_fp4,
                                                           scaled_diff_B_fp4_arr[subblock_idx]);
                scaled_diff_fp4         = _mm256_min_epi32(scaled_diff_fp4, seven_fp);

                //int32_t adjusted_weight = (expf_tab_fp16[scaled_diff_fp4] * TF_WEIGHT_SCALE) >> 16;
                adjusted_weight1 = _mm256_i32gather_epi32(
                    expf_tab_fp16, scaled_diff_fp4, 4);
                adjusted_weight1 = _mm256_mullo_epi32(adjusted_weight1, tf_weight_scale_fp);
                adjusted_weight1 = _mm256_srli_epi32(adjusted_weight1, 16); //>> 16
            } else {
                //int diff_sse = acc_5x5_sse[i][j];
                __m128i diff_sse1 = _mm_loadu_si128((__m128i *)(acc_5x5_sse[i] + j));
                __m128i diff_sse2 = _mm_loadu_si128((__m128i *)(acc_5x5_sse[i] + j + 4));

                //double    window_error = (double)diff_sse / num_ref_pixels;
                __m256d diff_sse_pd1  = _mm256_cvtepi32_pd(diff_sse1);
                __m256d diff_sse_pd2  = _mm256_cvtepi32_pd(diff_sse2);
                __m256d window_error1 = _mm256_div_pd(diff_sse_pd1, num_ref_pix);
                __m256d window_error2 = _mm256_div_pd(diff_sse_pd2, num_ref_pix);

                //const int subblock_idx = subblock_idx_h + (j >= block_width / 2);
                __m256d blk_error1 = blk_errors[subblock_idx_h + (j >= block_width / 2)];
                __m256d blk_error2 = blk_errors[subblock_idx_h + ((j + 4) >= block_width / 2)];

                //double combined_error = (TF_WINDOW_BLOCK_BALANCE_WEIGHT * window_error + block_error) * block_balacne_inv;
                __m256d combined_error1 = _mm256_mul_pd(window_error1, wnd_blk_balacne_weight);
                combined_error1         = _mm256_add_pd(combined_error1, blk_error1);
                combined_error1         = _mm256_mul_pd(combined_error1, blk_balacne_inv);

                __m256d combined_error2 = _mm256_mul_pd(window_error2, wnd_blk_balacne_weight);
                combined_error2         = _mm256_add_pd(combined_error2, blk_error2);
                combined_error2         = _mm256_mul_pd(combined_error2, blk_balacne_inv);

                //double scaled_diff     = AOMMIN(combined_error * d_factor_mul_n_decay_qr_inv, 7);
                //float   exp_mul_tf_scale = expf((float)(-scaled_diff)) * TF_WEIGHT_SCALE;
                //int     adjusted_weight  = (int)(exp_mul_tf_scale);

                __m256d d_fact_mul_n_decay1 =
                    d_factor_mul_n_decay_qr_invs[subblock_idx_h + (j >= block_width / 2)];
                __m256d d_fact_mul_n_decay2 =
                    d_factor_mul_n_decay_qr_invs[subblock_idx_h + ((j + 4) >= block_width / 2)];

                combined_error1 = _mm256_mul_pd(combined_error1, d_fact_mul_n_decay1);
                combined_error2 = _mm256_mul_pd(combined_error2, d_fact_mul_n_decay2);

                __m128 combined_error_ps1 = _mm256_cvtpd_ps(combined_error1);
                __m128 combined_error_ps2 = _mm256_cvtpd_ps(combined_error2);

                __m256 combined_error_ps = _mm256_insertf128_ps(
                    _mm256_castps128_ps256(combined_error_ps1), combined_error_ps2, 0x1);

                __m256 scaled_diff_ps = _mm256_min_ps(combined_error_ps, seven);
                scaled_diff_ps        = _mm256_sub_ps(zero, scaled_diff_ps); //-scaled_diff

                //int    adjusted_weight = (int)(expf((float)(-scaled_diff)) * TF_WEIGHT_SCALE);
                /*
                ** Semi-lossless optimization:
                **     Write SIMD friendly approximating algorithm for exp, this will require changes in "C" kernel too
                */
                scaled_diff_ps = exp_256_ps(scaled_diff_ps);
                scaled_diff_ps = _mm256_mul_ps(scaled_diff_ps, tf_weight_scale);

                adjusted_weight1 = _mm256_cvttps_epi32(scaled_diff_ps);
            }
#else /*FIXED_POINTS_PLANEWISE*/
            //int diff_sse = acc_5x5_sse[i][j];
            __m128i diff_sse1 = _mm_loadu_si128((__m128i *)(acc_5x5_sse[i] + j));
            __m128i diff_sse2 = _mm_loadu_si128((__m128i *)(acc_5x5_sse[i] + j + 4));

            //double    window_error = (double)diff_sse / num_ref_pixels;
            __m256d diff_sse_pd1  = _mm256_cvtepi32_pd(diff_sse1);
            __m256d diff_sse_pd2  = _mm256_cvtepi32_pd(diff_sse2);
            __m256d window_error1 = _mm256_div_pd(diff_sse_pd1, num_ref_pix);
            __m256d window_error2 = _mm256_div_pd(diff_sse_pd2, num_ref_pix);

            //const int subblock_idx = subblock_idx_h + (j >= block_width / 2);
            __m256d blk_error1 = blk_errors[subblock_idx_h + (j >= block_width / 2)];
            __m256d blk_error2 = blk_errors[subblock_idx_h + ((j + 4) >= block_width / 2)];

            //double combined_error = (TF_WINDOW_BLOCK_BALANCE_WEIGHT * window_error + block_error) * block_balacne_inv;
            __m256d combined_error1 = _mm256_mul_pd(window_error1, wnd_blk_balacne_weight);
            combined_error1         = _mm256_add_pd(combined_error1, blk_error1);
            combined_error1         = _mm256_mul_pd(combined_error1, blk_balacne_inv);

            __m256d combined_error2 = _mm256_mul_pd(window_error2, wnd_blk_balacne_weight);
            combined_error2         = _mm256_add_pd(combined_error2, blk_error2);
            combined_error2         = _mm256_mul_pd(combined_error2, blk_balacne_inv);

            //double scaled_diff     = AOMMIN(combined_error * d_factor_mul_n_decay_qr_inv, 7);
            //float   exp_mul_tf_scale = expf((float)(-scaled_diff)) * TF_WEIGHT_SCALE;
            //int     adjusted_weight  = (int)(exp_mul_tf_scale);

            __m256d d_fact_mul_n_decay1 =
                d_factor_mul_n_decay_qr_invs[subblock_idx_h + (j >= block_width / 2)];
            __m256d d_fact_mul_n_decay2 =
                d_factor_mul_n_decay_qr_invs[subblock_idx_h + ((j + 4) >= block_width / 2)];

            combined_error1 = _mm256_mul_pd(combined_error1, d_fact_mul_n_decay1);
            combined_error2 = _mm256_mul_pd(combined_error2, d_fact_mul_n_decay2);

            __m128 combined_error_ps1 = _mm256_cvtpd_ps(combined_error1);
            __m128 combined_error_ps2 = _mm256_cvtpd_ps(combined_error2);

            __m256 combined_error_ps = _mm256_insertf128_ps(
                _mm256_castps128_ps256(combined_error_ps1), combined_error_ps2, 0x1);

            __m256 scaled_diff_ps = _mm256_min_ps(combined_error_ps, seven);
            scaled_diff_ps        = _mm256_sub_ps(zero, scaled_diff_ps); //-scaled_diff

            //int    adjusted_weight = (int)(expf((float)(-scaled_diff)) * TF_WEIGHT_SCALE);
            /*
            ** Semi-lossless optimization:
            **     Write SIMD friendly approximating algorithm for exp, this will require changes in "C" kernel too
            */
            scaled_diff_ps = exp_256_ps(scaled_diff_ps);
            scaled_diff_ps = _mm256_mul_ps(scaled_diff_ps, tf_weight_scale);

            __m256i adjusted_weight1      = _mm256_cvttps_epi32(scaled_diff_ps);
#endif /*FIXED_POINTS_PLANEWISE*/
            __m128i adjusted_weight_int16 = _mm_packus_epi32(
                _mm256_castsi256_si128(adjusted_weight1),
                _mm256_extractf128_si256(adjusted_weight1, 0x1));

            //count[i * stride2 + j] += adjusted_weight;
            __m128i count_array = _mm_loadu_si128((__m128i *)(count + i * stride2 + j));
#if FIX_TEMPORAL_FILTER_PLANEWISE
            count_array         = _mm_add_epi16(count_array, adjusted_weight_int16);
#else
            count_array         = _mm_adds_epu16(count_array, adjusted_weight_int16);
#endif
            _mm_storeu_si128((__m128i *)(count + i * stride2 + j), count_array);

            //accumulator[i * stride2 + j] += adjusted_weight * frame2[i * stride2 + j];
            __m256i accumulator_array = _mm256_loadu_si256(
                (__m256i *)(accumulator + i * stride2 + j));
            __m128i frame2_array     = _mm_loadl_epi64((__m128i *)(frame2 + i * stride2 + j));
            frame2_array             = _mm_cvtepu8_epi16(frame2_array);
            __m256i frame2_array_u32 = _mm256_cvtepi16_epi32(frame2_array);
            frame2_array_u32         = _mm256_mullo_epi32(frame2_array_u32, adjusted_weight1);

            accumulator_array = _mm256_add_epi32(accumulator_array, frame2_array_u32);
            _mm256_storeu_si256((__m256i *)(accumulator + i * stride2 + j), accumulator_array);
        }
    }
}

void svt_av1_apply_temporal_filter_planewise_avx2(
    struct MeContext *context_ptr, const uint8_t *y_src, int y_src_stride, const uint8_t *y_pre,
    int y_pre_stride, const uint8_t *u_src, const uint8_t *v_src, int uv_src_stride,
    const uint8_t *u_pre, const uint8_t *v_pre, int uv_pre_stride, unsigned int block_width,
    unsigned int block_height, int ss_x, int ss_y,
#if !FTR_TF_STRENGTH_PER_QP
    const double *noise_levels,
    const int decay_control,
#endif
    uint32_t *y_accum, uint16_t *y_count, uint32_t *u_accum,
    uint16_t *u_count, uint32_t *v_accum, uint16_t *v_count) {
    // Loop variables
    assert(block_width <= BW && "block width too large");
    assert(block_height <= BH && "block height too large");
    assert(block_width % 16 == 0 && "block width must be multiple of 16");
    assert(block_height % 2 == 0 && "block height must be even");
    assert((ss_x == 0 || ss_x == 1) && (ss_y == 0 || ss_y == 1) && "invalid chroma subsampling");
    const int num_planes = context_ptr->tf_chroma ? 3 : 1;
    uint16_t  luma_sq_error[SSE_STRIDE * BH];
    uint16_t  chroma_sq_error[SSE_STRIDE * BH];

    for (int plane = 0; plane < num_planes; ++plane) {
        const uint32_t plane_h    = plane ? (block_height >> ss_y) : block_height;
        const uint32_t plane_w    = plane ? (block_width >> ss_x) : block_width;
        const uint32_t src_stride = plane ? uv_src_stride : y_src_stride;
        const uint32_t pre_stride = plane ? uv_pre_stride : y_pre_stride;
        const int      ss_x_shift = plane ? ss_x : 0;
        const int      ss_y_shift = plane ? ss_y : 0;

        const uint8_t *ref  = plane == 0 ? y_src : plane == 1 ? u_src : v_src;
        const uint8_t *pred = plane == 0 ? y_pre : plane == 1 ? u_pre : v_pre;

        uint32_t *accum = plane == 0 ? y_accum : plane == 1 ? u_accum : v_accum;
        uint16_t *count = plane == 0 ? y_count : plane == 1 ? u_count : v_count;

        apply_temporal_filter_planewise(context_ptr,
                                        ref,
                                        src_stride,
                                        pred,
                                        pre_stride,
                                        plane_w,
                                        plane_h,
#if !FTR_TF_STRENGTH_PER_QP
                                        noise_levels[plane],
                                        decay_control,
#endif
                                        accum,
                                        count,
                                        luma_sq_error,
                                        chroma_sq_error,
                                        plane,
                                        ss_x_shift,
                                        ss_y_shift);
    }
}

static AOM_FORCE_INLINE void get_squared_error_16x16_hbd_avx2(
    const uint16_t *frame1, const unsigned int stride, const uint16_t *frame2,
    const unsigned int stride2, const int block_width, const int block_height, uint32_t *frame_sse,
    const unsigned int sse_stride) {
    (void)block_width;
    const uint16_t *src1 = frame1;
    const uint16_t *src2 = frame2;
    uint32_t *      dst  = frame_sse;
    for (int i = 0; i < block_height; i++) {
        __m128i vf1_128, vf2_128;
        __m256i vf1, vf2, vdiff1, vsqdiff1;

        vf1_128  = _mm_loadu_si128((__m128i *)(src1));
        vf2_128  = _mm_loadu_si128((__m128i *)(src2));
        vf1      = _mm256_cvtepu16_epi32(vf1_128);
        vf2      = _mm256_cvtepu16_epi32(vf2_128);
        vdiff1   = _mm256_sub_epi32(vf1, vf2);
        vsqdiff1 = _mm256_mullo_epi32(vdiff1, vdiff1);
        _mm256_storeu_si256((__m256i *)(dst), vsqdiff1);

        vf1_128  = _mm_loadu_si128((__m128i *)(src1 + 8));
        vf2_128  = _mm_loadu_si128((__m128i *)(src2 + 8));
        vf1      = _mm256_cvtepu16_epi32(vf1_128);
        vf2      = _mm256_cvtepu16_epi32(vf2_128);
        vdiff1   = _mm256_sub_epi32(vf1, vf2);
        vsqdiff1 = _mm256_mullo_epi32(vdiff1, vdiff1);
        _mm256_storeu_si256((__m256i *)(dst + 8), vsqdiff1);

        *(uint32_t *)(dst + 16) = _mm_cvtsi128_si32(_mm_setzero_si128());

        src1 += stride, src2 += stride2;
        dst += sse_stride;
    }
}

static AOM_FORCE_INLINE void get_squared_error_32x32_hbd_avx2(
    const uint16_t *frame1, const unsigned int stride, const uint16_t *frame2,
    const unsigned int stride2, const int block_width, const int block_height, uint32_t *frame_sse,
    const unsigned int sse_stride) {
    (void)block_width;
    const uint16_t *src1 = frame1;
    const uint16_t *src2 = frame2;
    uint32_t *      dst  = frame_sse;
    for (int i = 0; i < block_height; i++) {
        __m128i vf1_128, vf2_128;
        __m256i vf1, vf2, vdiff1, vsqdiff1;

        vf1_128  = _mm_loadu_si128((__m128i *)(src1));
        vf2_128  = _mm_loadu_si128((__m128i *)(src2));
        vf1      = _mm256_cvtepu16_epi32(vf1_128);
        vf2      = _mm256_cvtepu16_epi32(vf2_128);
        vdiff1   = _mm256_sub_epi32(vf1, vf2);
        vsqdiff1 = _mm256_mullo_epi32(vdiff1, vdiff1);
        _mm256_storeu_si256((__m256i *)(dst), vsqdiff1);

        vf1_128  = _mm_loadu_si128((__m128i *)(src1 + 8));
        vf2_128  = _mm_loadu_si128((__m128i *)(src2 + 8));
        vf1      = _mm256_cvtepu16_epi32(vf1_128);
        vf2      = _mm256_cvtepu16_epi32(vf2_128);
        vdiff1   = _mm256_sub_epi32(vf1, vf2);
        vsqdiff1 = _mm256_mullo_epi32(vdiff1, vdiff1);
        _mm256_storeu_si256((__m256i *)(dst + 8), vsqdiff1);

        vf1_128  = _mm_loadu_si128((__m128i *)(src1 + 16));
        vf2_128  = _mm_loadu_si128((__m128i *)(src2 + 16));
        vf1      = _mm256_cvtepu16_epi32(vf1_128);
        vf2      = _mm256_cvtepu16_epi32(vf2_128);
        vdiff1   = _mm256_sub_epi32(vf1, vf2);
        vsqdiff1 = _mm256_mullo_epi32(vdiff1, vdiff1);
        _mm256_storeu_si256((__m256i *)(dst + 16), vsqdiff1);

        vf1_128  = _mm_loadu_si128((__m128i *)(src1 + 24));
        vf2_128  = _mm_loadu_si128((__m128i *)(src2 + 24));
        vf1      = _mm256_cvtepu16_epi32(vf1_128);
        vf2      = _mm256_cvtepu16_epi32(vf2_128);
        vdiff1   = _mm256_sub_epi32(vf1, vf2);
        vsqdiff1 = _mm256_mullo_epi32(vdiff1, vdiff1);
        _mm256_storeu_si256((__m256i *)(dst + 24), vsqdiff1);

        *(uint32_t *)(dst + 32) = _mm_cvtsi128_si32(_mm_setzero_si128());

        src1 += stride, src2 += stride2;
        dst += sse_stride;
    }
}

static AOM_FORCE_INLINE __m256i xx_load_and_pad_hbd(uint32_t *src, int col, int block_width) {
    __m256i v256tmp = _mm256_loadu_si256((__m256i *)(src));
    if (col == 0) {
        __m128i v128tmp1 = _mm256_castsi256_si128(v256tmp);
        __m128i v128tmp2 = _mm_loadu_si128((__m128i *)(src + 2));
        v128tmp1         = _mm_shuffle_epi8(v128tmp1, *(__m128i *)shufflemask_32b[0]);
        v256tmp          = _mm256_inserti128_si256(_mm256_castsi128_si256(v128tmp1), v128tmp2, 1);
    }
    if (col == block_width - 4) {
        __m128i v128tmp1 = _mm256_extracti128_si256(v256tmp, 1);
        v128tmp1         = _mm_shuffle_epi8(v128tmp1, *(__m128i *)shufflemask_32b[1]);
        v256tmp          = _mm256_inserti128_si256(v256tmp, v128tmp1, 1);
    }
    return v256tmp;
}

static void apply_temporal_filter_planewise_hbd(
    struct MeContext *context_ptr, const uint16_t *frame1, const unsigned int stride,
    const uint16_t *frame2, const unsigned int stride2, const int block_width,
    const int block_height,
#if !FTR_TF_STRENGTH_PER_QP
    const double sigma, const int decay_control,
#endif
    unsigned int *accumulator,
    uint16_t *count, uint32_t *luma_sq_error, uint32_t *chroma_sq_error, int plane, int ss_x_shift,
    int ss_y_shift, uint32_t encoder_bit_depth) {
    assert(TF_PLANEWISE_FILTER_WINDOW_LENGTH == 5);
    assert(((block_width == 32) && (block_height == 32)) ||
           ((block_width == 16) && (block_height == 16)));
    if (plane > PLANE_TYPE_Y)
        assert(chroma_sq_error != NULL);

    uint32_t acc_5x5_sse[BH][BW];
    // Larger noise -> larger filtering weight.
#if !FIXED_POINTS_PLANEWISE
#if FTR_TF_STRENGTH_PER_QP
    double n_decay_qr_inv = 1.0 / context_ptr->tf_decay_factor[plane];
#else
    const double n_decay                = (double)decay_control * (0.7 + log1p((double)sigma));
    const double n_decay_qr_inv         = 1.0 / (2 * n_decay * n_decay);
#endif
    const double block_balacne_inv      = 1.0 / (TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1);
    const double distance_threshold_inv = 1.0 /
        (double)AOMMAX(context_ptr->min_frame_size * TF_SEARCH_DISTANCE_THRESHOLD, 1);
#endif /*!FIXED_POINTS_PLANEWISE*/
    uint32_t *frame_sse = (plane == PLANE_TYPE_Y) ? luma_sq_error : chroma_sq_error;

    if (block_width == 32) {
        get_squared_error_32x32_hbd_avx2(
            frame1, stride, frame2, stride2, block_width, block_height, frame_sse, SSE_STRIDE);
    } else {
        get_squared_error_16x16_hbd_avx2(
            frame1, stride, frame2, stride2, block_width, block_height, frame_sse, SSE_STRIDE);
    }

    __m256i vsrc[5];

    // Traverse 4 columns at a time
    // First and last columns will require padding
    for (int col = 0; col < block_width; col += 4) {
        uint32_t *src = (col) ? frame_sse + col - 2 : frame_sse;

        // Load and pad(for first and last col) 3 rows from the top
        for (int i = 2; i < 5; i++) {
            vsrc[i] = xx_load_and_pad_hbd(src, col, block_width);
            src += SSE_STRIDE;
        }

        // Copy first row to first 2 vectors
        vsrc[0] = vsrc[2];
        vsrc[1] = vsrc[2];

        for (int row = 0; row < block_height; row++) {
            __m256i vsum = _mm256_setzero_si256();

            // Add 5 consecutive rows
            for (int i = 0; i < 5; i++) { vsum = _mm256_add_epi32(vsum, vsrc[i]); }

            // Push all elements by one element to the top
            for (int i = 0; i < 4; i++) { vsrc[i] = vsrc[i + 1]; }

            // Load next row to the last element
            if (row <= block_width - 4) {
                vsrc[4] = xx_load_and_pad_hbd(src, col, block_width);
                src += SSE_STRIDE;
            } else {
                vsrc[4] = vsrc[3];
            }

            // Accumulate the sum horizontally
            for (int i = 0; i < 4; i++) { acc_5x5_sse[row][col + i] = xx_mask_and_hadd(vsum, i); }
        }
    }

    int32_t idx_32x32      = context_ptr->tf_block_col + context_ptr->tf_block_row * 2;
#if FIXED_POINTS_PLANEWISE
    uint32_t num_ref_pixels = TF_PLANEWISE_FILTER_WINDOW_LENGTH * TF_PLANEWISE_FILTER_WINDOW_LENGTH;
#else
    int num_ref_pixels = TF_PLANEWISE_FILTER_WINDOW_LENGTH * TF_PLANEWISE_FILTER_WINDOW_LENGTH;
#endif
    int shift_factor   = ((encoder_bit_depth - 8) * 2);

    if (plane != PLANE_TYPE_Y) {
        num_ref_pixels += (1 << ss_y_shift) * (1 << ss_x_shift);

#if FIXED_POINTS_PLANEWISE
        if (context_ptr->tf_ctrls.use_fixed_point && (ss_x_shift == 1) && (ss_y_shift == 1)) {
            for (int i = 0; i < block_height; i++) {
                int yy = (i << 1);
                for (int j = 0; j < block_width; j += 4) {
                    int xx = (j << 1);
                    __m256i line = _mm256_loadu_si256(
                        (__m256i *)(luma_sq_error + yy * SSE_STRIDE + xx));
                    __m256i line2 = _mm256_loadu_si256(
                        (__m256i *)(luma_sq_error + (yy + 1) * SSE_STRIDE + xx));
                    __m256i sum = _mm256_hadd_epi32(line, line2);
                    sum         = _mm256_permute4x64_epi64(sum, 0xd8);

                    __m128i sum_128 = _mm_add_epi32(_mm256_castsi256_si128(sum),
                                                    _mm256_extracti128_si256(sum, 1));

                    __m128i diff_sse = _mm_loadu_si128((__m128i *)(acc_5x5_sse[i] + j));
                    diff_sse         = _mm_add_epi32(diff_sse, sum_128);
                    _mm_storeu_si128((__m128i *)(acc_5x5_sse[i] + j), diff_sse);
                }
            }
        } else {
            for (int i = 0; i < block_height; i++) {
                for (int j = 0; j < block_width; j++) {
                    for (int ii = 0; ii < (1 << ss_y_shift); ++ii) {
                        for (int jj = 0; jj < (1 << ss_x_shift); ++jj) {
                            const int yy = (i << ss_y_shift) + ii; // Y-coord on Y-plane.
                            const int xx = (j << ss_x_shift) + jj; // X-coord on Y-plane.
                            acc_5x5_sse[i][j] += luma_sq_error[yy * SSE_STRIDE + xx];
                        }
                    }
                }
            }
        }
#else
        for (int i = 0; i < block_height; i++) {
            for (int j = 0; j < block_width; j++) {
                for (int ii = 0; ii < (1 << ss_y_shift); ++ii) {
                    for (int jj = 0; jj < (1 << ss_x_shift); ++jj) {
                        const int yy = (i << ss_y_shift) + ii; // Y-coord on Y-plane.
                        const int xx = (j << ss_x_shift) + jj; // X-coord on Y-plane.
                        acc_5x5_sse[i][j] += luma_sq_error[yy * SSE_STRIDE + xx];
                    }
                }
            }
        }
#endif /*FIXED_POINTS_PLANEWISE*/
    }

    assert(!(block_width & 7));
#if FIXED_POINTS_PLANEWISE
    __m256i scaled_diff_B_fp4_arr[4], scaled_diff_A1_factor_fp14_arr[4], scaled_diff_A2_factor_fp10_arr[4];
    __m256i seven_fp   = _mm256_set1_epi32(7 * 16);
    __m256i tf_weight_scale_fp   = _mm256_set1_epi32(TF_WEIGHT_SCALE);
    int32_t move_arr[4];

    __m256d num_ref_pix            = _mm256_set1_pd((double)num_ref_pixels);
    __m256d blk_balacne_inv        = _mm256_set1_pd(1.0 / (TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1));
    __m256d wnd_blk_balacne_weight = _mm256_set1_pd((double)TF_WINDOW_BLOCK_BALANCE_WEIGHT);
    __m256  seven                  = _mm256_set1_ps(7.0f);
    __m256  zero                   = _mm256_set1_ps(0.0f);
    __m256  tf_weight_scale        = _mm256_set1_ps((float)TF_WEIGHT_SCALE);
    __m256d blk_errors[4];
    __m256d d_factor_mul_n_decay_qr_invs[4];
#else /*FIXED_POINTS_PLANEWISE*/
    __m256d num_ref_pix            = _mm256_set1_pd((double)num_ref_pixels);
    __m256d blk_balacne_inv        = _mm256_set1_pd(block_balacne_inv);
    __m256d wnd_blk_balacne_weight = _mm256_set1_pd((double)TF_WINDOW_BLOCK_BALANCE_WEIGHT);
    __m256  seven                  = _mm256_set1_ps(7.0f);
    __m256  zero                   = _mm256_set1_ps(0.0f);
    __m256  tf_weight_scale        = _mm256_set1_ps((float)TF_WEIGHT_SCALE);
    __m256d blk_errors[4];
    __m256d d_factor_mul_n_decay_qr_invs[4];
#endif /*FIXED_POINTS_PLANEWISE*/

    if (context_ptr->tf_32x32_block_split_flag[idx_32x32]) {
        for (int i = 0; i < 4; i++) {
#if FIXED_POINTS_PLANEWISE
            if (context_ptr->tf_ctrls.use_fixed_point) {
                uint32_t block_error_fp8 = (uint32_t)(
                    context_ptr->tf_16x16_block_error[idx_32x32 * 4 + i] >> 4);
                int32_t  col                     = context_ptr->tf_16x16_mv_x[idx_32x32 * 4 + i];
                int32_t  row                     = context_ptr->tf_16x16_mv_y[idx_32x32 * 4 + i];
                uint32_t distance_fp4            = sqrt_fast(((uint32_t)(col * col + row * row)) << 8);
                uint32_t  distance_threshold_fp16 = AOMMAX(
                    (context_ptr->min_frame_size << 16) / 10, 1 << 16);

                int32_t d_factor_fp12 = AOMMAX(
                    (int32_t)(((((int64_t)distance_fp4) << 24) + (distance_threshold_fp16 >> 1)) /
                              (distance_threshold_fp16)),
                    1 << 12);

                int32_t scaled_diff_B_fp4 = (int32_t)(
                    (((int64_t)d_factor_fp12) *
                     ((block_error_fp8 + ((TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1) / 2)) /
                      (TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1))) /
                    ((context_ptr->tf_decay_factor_fp16[plane])));

                int32_t scaled_diff_A1_factor_fp14 = (int32_t)(
                    ((((int64_t)d_factor_fp12) << 16) * TF_WINDOW_BLOCK_BALANCE_WEIGHT +
                     (context_ptr->tf_decay_factor_fp16[plane] >> 3)) /
                    (context_ptr->tf_decay_factor_fp16[plane] >> 2));

                int32_t scaled_diff_A2_factor_fp10 =
                    ((1 << 10) + ((TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1) * num_ref_pixels) / 2) /
                    ((TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1) * num_ref_pixels);

                int move = 0;
                if (scaled_diff_A1_factor_fp14 > (1 << 20)) {
                    move = 18;
                } else if (scaled_diff_A1_factor_fp14 > (1 << 16)) {
                    move = 14;
                } else if (scaled_diff_A1_factor_fp14 > (1 << 12)) {
                    move = 10;
                } else if (scaled_diff_A1_factor_fp14 > (1 << 8)) {
                    move = 6;
                } else if (scaled_diff_A1_factor_fp14 > (1 << 4)) {
                    move = 2;
                } else {
                    move = 0;
                }

                scaled_diff_B_fp4_arr[i]          = _mm256_set1_epi32(scaled_diff_B_fp4);
                scaled_diff_A1_factor_fp14_arr[i] = _mm256_set1_epi32(scaled_diff_A1_factor_fp14);
                scaled_diff_A2_factor_fp10_arr[i] = _mm256_set1_epi32(scaled_diff_A2_factor_fp10);
                move_arr[i]                       = move;
            } else {
                blk_errors[i] = _mm256_set1_pd(
                    (double)(context_ptr->tf_16x16_block_error[idx_32x32 * 4 + i] >> 4) / 256.0);

                int16_t col      = context_ptr->tf_16x16_mv_x[idx_32x32 * 4 + i];
                int16_t row      = context_ptr->tf_16x16_mv_y[idx_32x32 * 4 + i];
                float   distance = sqrtf((float)(row * row + col * col));
                double  d_factor = AOMMAX(distance *(1.0 / (double)AOMMAX(context_ptr->min_frame_size *
                    TF_SEARCH_DISTANCE_THRESHOLD, 1)), 1);
                d_factor_mul_n_decay_qr_invs[i] = _mm256_set1_pd(
                    d_factor * (1.0 / context_ptr->tf_decay_factor[plane]));
            }
#else /*FIXED_POINTS_PLANEWISE*/
            blk_errors[i] = _mm256_set1_pd(
                (double)(context_ptr->tf_16x16_block_error[idx_32x32 * 4 + i] >> 4) / 256.0);

            int16_t col                     = context_ptr->tf_16x16_mv_x[idx_32x32 * 4 + i];
            int16_t row                     = context_ptr->tf_16x16_mv_y[idx_32x32 * 4 + i];
            float   distance                = sqrtf((float)(row * row + col * col));
            double  d_factor                = AOMMAX(distance * distance_threshold_inv, 1);
            d_factor_mul_n_decay_qr_invs[i] = _mm256_set1_pd(d_factor * n_decay_qr_inv);
#endif /*FIXED_POINTS_PLANEWISE*/
        }
    } else {
#if FIXED_POINTS_PLANEWISE
        if (context_ptr->tf_ctrls.use_fixed_point) {
            uint32_t block_error_fp8 = (uint32_t)(context_ptr->tf_32x32_block_error[idx_32x32] >>
                                                  6);
            int32_t  col             = context_ptr->tf_32x32_mv_x[idx_32x32];
            int32_t  row             = context_ptr->tf_32x32_mv_y[idx_32x32];
            uint32_t distance_fp4    = sqrt_fast(((uint32_t)(col * col + row * row)) << 8);
            uint32_t distance_threshold_fp16 = AOMMAX((context_ptr->min_frame_size << 16) / 10,
                                                     1 << 16);

            int32_t d_factor_fp12 = AOMMAX(
                (int32_t)(((((int64_t)distance_fp4) << 24) + (distance_threshold_fp16 >> 1)) /
                          (distance_threshold_fp16)),
                1 << 12);

            int32_t scaled_diff_B_fp4 = (int32_t)(
                (((int64_t)d_factor_fp12) *
                 ((block_error_fp8 + ((TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1) / 2)) /
                  (TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1))) /
                ((context_ptr->tf_decay_factor_fp16[plane])));

            int32_t scaled_diff_A1_factor_fp14 = (int32_t)(
                ((((int64_t)d_factor_fp12) << 16) * TF_WINDOW_BLOCK_BALANCE_WEIGHT +
                 (context_ptr->tf_decay_factor_fp16[plane] >> 3)) /
                (context_ptr->tf_decay_factor_fp16[plane] >> 2));

            int32_t scaled_diff_A2_factor_fp10 =
                ((1 << 10) + ((TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1) * num_ref_pixels) / 2) /
                ((TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1) * num_ref_pixels);

            int move = 0;
            if (scaled_diff_A1_factor_fp14 > (1 << 20)) {
                move = 18;
            } else if (scaled_diff_A1_factor_fp14 > (1 << 16)) {
                move = 14;
            } else if (scaled_diff_A1_factor_fp14 > (1 << 12)) {
                move = 10;
            } else if (scaled_diff_A1_factor_fp14 > (1 << 8)) {
                move = 6;
            } else if (scaled_diff_A1_factor_fp14 > (1 << 4)) {
                move = 2;
            } else {
                move = 0;
            }

            for (int i = 0; i < 4; i++) {
                scaled_diff_B_fp4_arr[i]          = _mm256_set1_epi32(scaled_diff_B_fp4);
                scaled_diff_A1_factor_fp14_arr[i] = _mm256_set1_epi32(scaled_diff_A1_factor_fp14);
                scaled_diff_A2_factor_fp10_arr[i] = _mm256_set1_epi32(scaled_diff_A2_factor_fp10);
                move_arr[i]                       = move;
            }
        } else {
            double block_error = (double)(context_ptr->tf_32x32_block_error[idx_32x32] >> 4) /
                1024.0;
            int16_t      col      = context_ptr->tf_32x32_mv_x[idx_32x32];
            int16_t      row      = context_ptr->tf_32x32_mv_y[idx_32x32];
            const float  distance = sqrtf((float)(row * row + col * col));
                double  d_factor = AOMMAX(distance *(1.0 / (double)AOMMAX(context_ptr->min_frame_size *
                    TF_SEARCH_DISTANCE_THRESHOLD, 1)), 1);

            const double d_factor_mul_n_decay_qr_inv = d_factor * (1.0 / context_ptr->tf_decay_factor[plane]);

            for (int i = 0; i < 4; i++) {
                blk_errors[i]                   = _mm256_set1_pd(block_error);
                d_factor_mul_n_decay_qr_invs[i] = _mm256_set1_pd(d_factor_mul_n_decay_qr_inv);
            }
        }
#else /*FIXED_POINTS_PLANEWISE*/
        double  block_error = (double)(context_ptr->tf_32x32_block_error[idx_32x32] >> 4) / 1024.0;
        int16_t col         = context_ptr->tf_32x32_mv_x[idx_32x32];
        int16_t row         = context_ptr->tf_32x32_mv_y[idx_32x32];
        const float  distance = sqrtf((float)(row * row + col * col));
        const double d_factor = AOMMAX(distance * distance_threshold_inv, 1);

        const double d_factor_mul_n_decay_qr_inv = d_factor * n_decay_qr_inv;

        for (int i = 0; i < 4; i++) {
            blk_errors[i]                   = _mm256_set1_pd(block_error);
            d_factor_mul_n_decay_qr_invs[i] = _mm256_set1_pd(d_factor_mul_n_decay_qr_inv);
        }
#endif /*FIXED_POINTS_PLANEWISE*/
    }

    for (int i = 0; i < block_height; i++) {
        const int subblock_idx_h = (i >= block_height / 2) * 2;
        for (int j = 0; j < block_width; j += 8) {
            //int diff_sse = acc_5x5_sse[i][j];
#if FIXED_POINTS_PLANEWISE
            __m256i adjusted_weight1;
            if (context_ptr->tf_ctrls.use_fixed_point) {
                const int subblock_idx = subblock_idx_h + (j >= block_width / 2);

                __m256i diff_sse = _mm256_loadu_si256((__m256i *)(acc_5x5_sse[i] + j));
                diff_sse         = _mm256_srli_epi32(diff_sse, shift_factor);

                //int32_t sum_square_part_fp10 = (((int32_t)sum_square_diff) * scaled_diff_A2_factor_fp10);
                __m256i sum_square_part_fp10 = _mm256_mullo_epi32(
                    diff_sse, scaled_diff_A2_factor_fp10_arr[subblock_idx]);

                //int32_t scaled_diff_A_fp4 = (int32_t)(((sum_square_part_fp10) * ((scaled_diff_A1_factor_fp14)>>move) + (1<<(19-move)))>>(20-move));
                __m256i scaled_diff_A1_factor_fp14_avx2 = _mm256_srli_epi32(
                    scaled_diff_A1_factor_fp14_arr[subblock_idx], move_arr[subblock_idx]);

                __m256i scaled_diff_A_fp4 = _mm256_mullo_epi32(sum_square_part_fp10,
                                                               scaled_diff_A1_factor_fp14_avx2);
                scaled_diff_A_fp4         = _mm256_add_epi32(
                    scaled_diff_A_fp4, _mm256_set1_epi32(1 << (19 - move_arr[subblock_idx])));
                scaled_diff_A_fp4 = _mm256_srli_epi32(scaled_diff_A_fp4,
                                                      20 - move_arr[subblock_idx]);

                //int32_t scaled_diff_fp4 = AOMMIN(scaled_diff_A_fp4 + scaled_diff_B_fp4, 7*16);
                __m256i scaled_diff_fp4 = _mm256_add_epi32(scaled_diff_A_fp4,
                                                           scaled_diff_B_fp4_arr[subblock_idx]);
                scaled_diff_fp4         = _mm256_min_epi32(scaled_diff_fp4, seven_fp);

                //int32_t adjusted_weight = (expf_tab_fp16[scaled_diff_fp4] * TF_WEIGHT_SCALE) >> 16;
                adjusted_weight1 = _mm256_i32gather_epi32(expf_tab_fp16, scaled_diff_fp4, 4);
                adjusted_weight1 = _mm256_mullo_epi32(adjusted_weight1, tf_weight_scale_fp);
                adjusted_weight1 = _mm256_srli_epi32(adjusted_weight1, 16); //>> 16
            } else {
                __m128i diff_sse1 = _mm_loadu_si128((__m128i *)(acc_5x5_sse[i] + j));
                __m128i diff_sse2 = _mm_loadu_si128((__m128i *)(acc_5x5_sse[i] + j + 4));

                //diff_sse >>= ((encoder_bit_depth - 8) * 2);
                diff_sse1 = _mm_srli_epi32(diff_sse1, shift_factor);
                diff_sse2 = _mm_srli_epi32(diff_sse2, shift_factor);

                //double    window_error = (double)diff_sse / num_ref_pixels;
                __m256d diff_sse_pd1  = _mm256_cvtepi32_pd(diff_sse1);
                __m256d diff_sse_pd2  = _mm256_cvtepi32_pd(diff_sse2);
                __m256d window_error1 = _mm256_div_pd(diff_sse_pd1, num_ref_pix);
                __m256d window_error2 = _mm256_div_pd(diff_sse_pd2, num_ref_pix);

                //const int subblock_idx = subblock_idx_h + (j >= block_width / 2);
                __m256d blk_error1 = blk_errors[subblock_idx_h + (j >= block_width / 2)];
                __m256d blk_error2 = blk_errors[subblock_idx_h + ((j + 4) >= block_width / 2)];

                //double combined_error = (TF_WINDOW_BLOCK_BALANCE_WEIGHT * window_error + block_error) * block_balacne_inv;
                __m256d combined_error1 = _mm256_mul_pd(window_error1, wnd_blk_balacne_weight);
                combined_error1         = _mm256_add_pd(combined_error1, blk_error1);
                combined_error1         = _mm256_mul_pd(combined_error1, blk_balacne_inv);

                __m256d combined_error2 = _mm256_mul_pd(window_error2, wnd_blk_balacne_weight);
                combined_error2         = _mm256_add_pd(combined_error2, blk_error2);
                combined_error2         = _mm256_mul_pd(combined_error2, blk_balacne_inv);

                //double scaled_diff     = AOMMIN(combined_error * d_factor_mul_n_decay_qr_inv, 7);
                //float   exp_mul_tf_scale = expf((float)(-scaled_diff)) * TF_WEIGHT_SCALE;
                //int     adjusted_weight  = (int)(exp_mul_tf_scale);

                __m256d d_fact_mul_n_decay1 =
                    d_factor_mul_n_decay_qr_invs[subblock_idx_h + (j >= block_width / 2)];
                __m256d d_fact_mul_n_decay2 =
                    d_factor_mul_n_decay_qr_invs[subblock_idx_h + ((j + 4) >= block_width / 2)];

                combined_error1 = _mm256_mul_pd(combined_error1, d_fact_mul_n_decay1);
                combined_error2 = _mm256_mul_pd(combined_error2, d_fact_mul_n_decay2);

                __m128 combined_error_ps1 = _mm256_cvtpd_ps(combined_error1);
                __m128 combined_error_ps2 = _mm256_cvtpd_ps(combined_error2);

                __m256 combined_error_ps = _mm256_insertf128_ps(
                    _mm256_castps128_ps256(combined_error_ps1), combined_error_ps2, 0x1);

                __m256 scaled_diff_ps = _mm256_min_ps(combined_error_ps, seven);
                scaled_diff_ps        = _mm256_sub_ps(zero, scaled_diff_ps); //-scaled_diff

                //int    adjusted_weight = (int)(expf((float)(-scaled_diff)) * TF_WEIGHT_SCALE);
                /*
            ** Semi-lossless optimization:
            **     Write SIMD friendly approximating algorithm for exp, this will require changes in "C" kernel too
            */
                scaled_diff_ps = exp_256_ps(scaled_diff_ps);
                scaled_diff_ps = _mm256_mul_ps(scaled_diff_ps, tf_weight_scale);

                adjusted_weight1 = _mm256_cvttps_epi32(scaled_diff_ps);
            }
#else /*FIXED_POINTS_PLANEWISE*/
            __m128i diff_sse1 = _mm_loadu_si128((__m128i *)(acc_5x5_sse[i] + j));
            __m128i diff_sse2 = _mm_loadu_si128((__m128i *)(acc_5x5_sse[i] + j + 4));

            //diff_sse >>= ((encoder_bit_depth - 8) * 2);
            diff_sse1 = _mm_srli_epi32(diff_sse1, shift_factor);
            diff_sse2 = _mm_srli_epi32(diff_sse2, shift_factor);

            //double    window_error = (double)diff_sse / num_ref_pixels;
            __m256d diff_sse_pd1  = _mm256_cvtepi32_pd(diff_sse1);
            __m256d diff_sse_pd2  = _mm256_cvtepi32_pd(diff_sse2);
            __m256d window_error1 = _mm256_div_pd(diff_sse_pd1, num_ref_pix);
            __m256d window_error2 = _mm256_div_pd(diff_sse_pd2, num_ref_pix);

            //const int subblock_idx = subblock_idx_h + (j >= block_width / 2);
            __m256d blk_error1 = blk_errors[subblock_idx_h + (j >= block_width / 2)];
            __m256d blk_error2 = blk_errors[subblock_idx_h + ((j + 4) >= block_width / 2)];

            //double combined_error = (TF_WINDOW_BLOCK_BALANCE_WEIGHT * window_error + block_error) * block_balacne_inv;
            __m256d combined_error1 = _mm256_mul_pd(window_error1, wnd_blk_balacne_weight);
            combined_error1         = _mm256_add_pd(combined_error1, blk_error1);
            combined_error1         = _mm256_mul_pd(combined_error1, blk_balacne_inv);

            __m256d combined_error2 = _mm256_mul_pd(window_error2, wnd_blk_balacne_weight);
            combined_error2         = _mm256_add_pd(combined_error2, blk_error2);
            combined_error2         = _mm256_mul_pd(combined_error2, blk_balacne_inv);

            //double scaled_diff     = AOMMIN(combined_error * d_factor_mul_n_decay_qr_inv, 7);
            //float   exp_mul_tf_scale = expf((float)(-scaled_diff)) * TF_WEIGHT_SCALE;
            //int     adjusted_weight  = (int)(exp_mul_tf_scale);

            __m256d d_fact_mul_n_decay1 =
                d_factor_mul_n_decay_qr_invs[subblock_idx_h + (j >= block_width / 2)];
            __m256d d_fact_mul_n_decay2 =
                d_factor_mul_n_decay_qr_invs[subblock_idx_h + ((j + 4) >= block_width / 2)];

            combined_error1 = _mm256_mul_pd(combined_error1, d_fact_mul_n_decay1);
            combined_error2 = _mm256_mul_pd(combined_error2, d_fact_mul_n_decay2);

            __m128 combined_error_ps1 = _mm256_cvtpd_ps(combined_error1);
            __m128 combined_error_ps2 = _mm256_cvtpd_ps(combined_error2);

            __m256 combined_error_ps = _mm256_insertf128_ps(
                _mm256_castps128_ps256(combined_error_ps1), combined_error_ps2, 0x1);

            __m256 scaled_diff_ps = _mm256_min_ps(combined_error_ps, seven);
            scaled_diff_ps        = _mm256_sub_ps(zero, scaled_diff_ps); //-scaled_diff

            //int    adjusted_weight = (int)(expf((float)(-scaled_diff)) * TF_WEIGHT_SCALE);
            /*
            ** Semi-lossless optimization:
            **     Write SIMD friendly approximating algorithm for exp, this will require changes in "C" kernel too
            */
            scaled_diff_ps = exp_256_ps(scaled_diff_ps);
            scaled_diff_ps = _mm256_mul_ps(scaled_diff_ps, tf_weight_scale);

            __m256i adjusted_weight1      = _mm256_cvttps_epi32(scaled_diff_ps);
#endif /*FIXED_POINTS_PLANEWISE*/
            __m128i adjusted_weight_int16 = _mm_packus_epi32(
                _mm256_castsi256_si128(adjusted_weight1),
                _mm256_extractf128_si256(adjusted_weight1, 0x1));

            //count[i * stride2 + j] += adjusted_weight;
            __m128i count_array = _mm_loadu_si128((__m128i *)(count + i * stride2 + j));
#if FIX_TEMPORAL_FILTER_PLANEWISE
            count_array         = _mm_add_epi16(count_array, adjusted_weight_int16);
#else
            count_array         = _mm_adds_epu16(count_array, adjusted_weight_int16);
#endif

            _mm_storeu_si128((__m128i *)(count + i * stride2 + j), count_array);

            //accumulator[i * stride2 + j] += adjusted_weight * frame2[i * stride2 + j];
            __m256i accumulator_array = _mm256_loadu_si256(
                (__m256i *)(accumulator + i * stride2 + j));
            __m128i frame2_array     = _mm_loadu_si128((__m128i *)(frame2 + i * stride2 + j));
            __m256i frame2_array_u32 = _mm256_cvtepi16_epi32(frame2_array);
            frame2_array_u32         = _mm256_mullo_epi32(frame2_array_u32, adjusted_weight1);

            accumulator_array = _mm256_add_epi32(accumulator_array, frame2_array_u32);
            _mm256_storeu_si256((__m256i *)(accumulator + i * stride2 + j), accumulator_array);
        }
    }
}

void svt_av1_apply_temporal_filter_planewise_hbd_avx2(
    struct MeContext *context_ptr, const uint16_t *y_src, int y_src_stride, const uint16_t *y_pre,
    int y_pre_stride, const uint16_t *u_src, const uint16_t *v_src, int uv_src_stride,
    const uint16_t *u_pre, const uint16_t *v_pre, int uv_pre_stride, unsigned int block_width,
    unsigned int block_height, int ss_x, int ss_y,
#if !FTR_TF_STRENGTH_PER_QP
    const double *noise_levels,
    const int decay_control,
#endif
    uint32_t *y_accum, uint16_t *y_count, uint32_t *u_accum,
    uint16_t *u_count, uint32_t *v_accum, uint16_t *v_count, uint32_t encoder_bit_depth) {
    // Loop variables
    assert(block_width <= BW && "block width too large");
    assert(block_height <= BH && "block height too large");
    assert(block_width % 16 == 0 && "block width must be multiple of 16");
    assert(block_height % 2 == 0 && "block height must be even");
    assert((ss_x == 0 || ss_x == 1) && (ss_y == 0 || ss_y == 1) && "invalid chroma subsampling");

    const int num_planes = context_ptr->tf_chroma ? 3 : 1;
    uint32_t  luma_sq_error[SSE_STRIDE * BH];
    uint32_t  chroma_sq_error[SSE_STRIDE * BH];

    for (int plane = 0; plane < num_planes; ++plane) {
        const uint32_t plane_h    = plane ? (block_height >> ss_y) : block_height;
        const uint32_t plane_w    = plane ? (block_width >> ss_x) : block_width;
        const uint32_t src_stride = plane ? uv_src_stride : y_src_stride;
        const uint32_t pre_stride = plane ? uv_pre_stride : y_pre_stride;
        const int      ss_x_shift = plane ? ss_x : 0;
        const int      ss_y_shift = plane ? ss_y : 0;

        const uint16_t *ref  = plane == 0 ? y_src : plane == 1 ? u_src : v_src;
        const uint16_t *pred = plane == 0 ? y_pre : plane == 1 ? u_pre : v_pre;

        uint32_t *accum = plane == 0 ? y_accum : plane == 1 ? u_accum : v_accum;
        uint16_t *count = plane == 0 ? y_count : plane == 1 ? u_count : v_count;
        apply_temporal_filter_planewise_hbd(context_ptr,
                                            ref,
                                            src_stride,
                                            pred,
                                            pre_stride,
                                            plane_w,
                                            plane_h,
#if !FTR_TF_STRENGTH_PER_QP
                                            noise_levels[plane],
                                            decay_control,
#endif
                                            accum,
                                            count,
                                            luma_sq_error,
                                            chroma_sq_error,
                                            plane,
                                            ss_x_shift,
                                            ss_y_shift,
                                            encoder_bit_depth);
    }
}
#if OPT_TFILTER
uint32_t calculate_squared_errors_sum_avx2(const uint8_t *s, int s_stride, const uint8_t *p,
    int p_stride, unsigned int w, unsigned int h) {
    assert(w % 16 == 0 && "block width must be multiple of 16");
    unsigned int i, j;

    __m256i sum = _mm256_setzero_si256();
    __m256i s_16, p_16, dif;

    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j += 16) {
            s_16 = _mm256_cvtepu8_epi16(_mm_loadu_si128((const __m128i *)(s + i * s_stride + j)));
            p_16 = _mm256_cvtepu8_epi16(_mm_loadu_si128((const __m128i *)(p + i * p_stride + j)));

            dif = _mm256_sub_epi16(s_16, p_16);
            sum = _mm256_add_epi32(sum, _mm256_madd_epi16(dif, dif));
        }
    }
    __m128i sum_128 = _mm_add_epi32(_mm256_castsi256_si128(sum),
        _mm256_extracti128_si256(sum, 0x1));
    sum_128 = _mm_hadd_epi32(sum_128, sum_128);
    sum_128 = _mm_hadd_epi32(sum_128, sum_128);

#if FIX_TEMPORAL_FILTER_PLANEWISE
    return _mm_cvtsi128_si32(sum_128);
#else
    return _mm_cvtsi128_si32(sum_128) / (w * h);
#endif
}

uint32_t  calculate_squared_errors_sum_highbd_avx2(const uint16_t *s, int s_stride,
    const uint16_t *p, int p_stride,
#if FIX_TEMPORAL_FILTER_PLANEWISE
    unsigned int w, unsigned int h, int shift_factor) {
#else
    unsigned int w, unsigned int h) {
#endif
    assert(w % 16 == 0 && "block width must be multiple of 16");
    unsigned int i, j;

    __m256i sum = _mm256_setzero_si256();
    __m256i s_16, p_16, dif;

    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j += 16) {
            s_16 = _mm256_loadu_si256((const __m256i *)(s + i * s_stride + j));
            p_16 = _mm256_loadu_si256((const __m256i *)(p + i * p_stride + j));

            dif = _mm256_sub_epi16(s_16, p_16);
            sum = _mm256_add_epi32(sum, _mm256_madd_epi16(dif, dif));
        }
    }
    __m128i sum_128 = _mm_add_epi32(_mm256_castsi256_si128(sum),
        _mm256_extracti128_si256(sum, 0x1));
    sum_128 = _mm_hadd_epi32(sum_128, sum_128);
    sum_128 = _mm_hadd_epi32(sum_128, sum_128);

#if FIX_TEMPORAL_FILTER_PLANEWISE
    return (_mm_cvtsi128_si32(sum_128) >> shift_factor);
#else
    return _mm_cvtsi128_si32(sum_128) / (w * h);
#endif
}

//exp(-x) for x in [0..7]
static double expf_tab[] = {
            1,
0.904837,
0.818731,
0.740818,
0.67032 ,
0.606531,
0.548812,
0.496585,
0.449329,
0.40657 ,
0.367879,
0.332871,
0.301194,
0.272532,
0.246597,
0.22313 ,
0.201896,
0.182683,
0.165299,
0.149569,
0.135335,
0.122456,
0.110803,
0.100259,
0.090718,
0.082085,
0.074274,
0.067206,
0.06081 ,
0.055023,
0.049787,
0.045049,
0.040762,
0.036883,
0.033373,
0.030197,
0.027324,
0.024724,
0.022371,
0.020242,
0.018316,
0.016573,
0.014996,
0.013569,
0.012277,
0.011109,
0.010052,
0.009095,
0.00823 ,
0.007447,
0.006738,
0.006097,
0.005517,
0.004992,
0.004517,
0.004087,
0.003698,
0.003346,
0.003028,
0.002739,
0.002479,
0.002243,
0.002029,
0.001836,
0.001662,
0.001503,
0.00136 ,
0.001231,
0.001114,
0.001008,
0.000912,
0.000825,
0.000747,
0.000676,
0.000611,
0.000553,
0.0005  ,
0.000453,
0.00041 ,
0.000371,
0.000335

};

void svt_av1_apply_temporal_filter_planewise_fast_avx2(
        struct MeContext *context_ptr, const uint8_t *y_src, int y_src_stride, const uint8_t *y_pre,
        int y_pre_stride, unsigned int block_width,
        unsigned int block_height, uint32_t *y_accum, uint16_t *y_count)
{
#if FIX_TEMPORAL_FILTER_PLANEWISE
    //to add chroma if need be
    uint32_t avg_err = calculate_squared_errors_sum_avx2(
        y_src, y_src_stride, y_pre, y_pre_stride, block_width, block_height) / (block_width * block_height);
#else
    //to add chroma if need be
    uint32_t avg_err = calculate_squared_errors_sum_avx2(
        y_src, y_src_stride, y_pre, y_pre_stride, block_width, block_height);
#endif

#if FIXED_POINTS_PLANEWISE
    int adjusted_weight;
    if (context_ptr->tf_ctrls.use_fixed_point) {
        //16*avg_err/context_ptr->tf_decay_factor[0];
        uint32_t scaled_diff_fp4 = AOMMIN((avg_err<<10) / (context_ptr->tf_decay_factor_fp16[0]>>10), 7*16);
        adjusted_weight = (expf_tab_fp16[scaled_diff_fp4] * TF_WEIGHT_SCALE) >> 16;
    } else {
#endif /*FIXED_POINTS_PLANEWISE*/
#if FTR_TF_STRENGTH_PER_QP
    double scaled_diff = AOMMIN(avg_err / context_ptr->tf_decay_factor[0], 7);
#else
    double scaled_diff = AOMMIN(avg_err / context_ptr->tf_decay_factor, 7);
#endif
#if FIXED_POINTS_PLANEWISE
        adjusted_weight = (int)(expf_tab[(int)(scaled_diff * 10)] * TF_WEIGHT_SCALE);
    }
#else
    int    adjusted_weight = (int)(expf_tab[(int)(scaled_diff * 10)] * TF_WEIGHT_SCALE);
#endif /*FIXED_POINTS_PLANEWISE*/
    const __m128i adjusted_weight_int16 = _mm_set1_epi16((int16_t)(adjusted_weight));
    const __m256i adjusted_weight_int32 = _mm256_set1_epi32((int32_t)(adjusted_weight));
    if (adjusted_weight) {
        unsigned int       i, j, k;
        for (i = 0; i < block_height; i++) {
            for (j = 0; j < block_width; j += 8) {
                k = i * y_pre_stride + j;

                //y_count[k] += adjusted_weight;
                __m128i count_array = _mm_loadu_si128((__m128i *)(y_count + k));
#if FIX_TEMPORAL_FILTER_PLANEWISE
                count_array = _mm_add_epi16(count_array, adjusted_weight_int16);
#else
                count_array = _mm_adds_epu16(count_array, adjusted_weight_int16);
#endif
                _mm_storeu_si128((__m128i *)(y_count + k), count_array);

                //y_accum[k] += adjusted_weight * pixel_value;
                __m256i accumulator_array = _mm256_loadu_si256((__m256i *)(y_accum + k));
                __m128i frame2_array = _mm_loadl_epi64((__m128i *)(y_pre + k));
                frame2_array = _mm_cvtepu8_epi16(frame2_array);
                __m256i frame2_array_u32 = _mm256_cvtepi16_epi32(frame2_array);
                frame2_array_u32 = _mm256_mullo_epi32(frame2_array_u32, adjusted_weight_int32);

                accumulator_array = _mm256_add_epi32(accumulator_array, frame2_array_u32);
                _mm256_storeu_si256((__m256i *)(y_accum + k), accumulator_array);
            }
        }
    }
}

void svt_av1_apply_temporal_filter_planewise_fast_hbd_avx2(
    struct MeContext *context_ptr, const uint16_t *y_src, int y_src_stride, const uint16_t *y_pre,
    int y_pre_stride, unsigned int block_width,
#if FIX_TEMPORAL_FILTER_PLANEWISE
    unsigned int block_height, uint32_t *y_accum, uint16_t *y_count, uint32_t encoder_bit_depth)
#else
    unsigned int block_height, uint32_t *y_accum, uint16_t *y_count)
#endif
{
#if FIX_TEMPORAL_FILTER_PLANEWISE
    int shift_factor  = ((encoder_bit_depth - 8) * 2);
    //to add chroma if need be
    uint32_t avg_err = calculate_squared_errors_sum_highbd_avx2(
        y_src, y_src_stride, y_pre, y_pre_stride, block_width, block_height, shift_factor)
            / (block_width * block_height);
#else
    //to add chroma if need be
    uint32_t avg_err = calculate_squared_errors_sum_highbd_avx2(
        y_src, y_src_stride, y_pre, y_pre_stride, block_width, block_height);
#endif

#if FIXED_POINTS_PLANEWISE
    int adjusted_weight;
    if (context_ptr->tf_ctrls.use_fixed_point) {
        //16*avg_err/context_ptr->tf_decay_factor[0];
        uint32_t scaled_diff_fp4 = AOMMIN(
            (avg_err << 10) / (context_ptr->tf_decay_factor_fp16[0] >> 10), 7 * 16);
        adjusted_weight = (expf_tab_fp16[scaled_diff_fp4] * TF_WEIGHT_SCALE) >> 16;
    } else {
#endif
#if FTR_TF_STRENGTH_PER_QP
    double scaled_diff = AOMMIN(avg_err / context_ptr->tf_decay_factor[0], 7);
#else
    double scaled_diff = AOMMIN(avg_err / context_ptr->tf_decay_factor, 7);
#endif
#if FIXED_POINTS_PLANEWISE
    adjusted_weight = (int)(expf_tab[(int)(scaled_diff * 10)] * TF_WEIGHT_SCALE);
    }
#else
    int adjusted_weight = (int)(expf_tab[(int)(scaled_diff * 10)] * TF_WEIGHT_SCALE);
#endif
    if (adjusted_weight) {
        unsigned int       i, j, k;
        const __m128i adjusted_weight_int16 = _mm_set1_epi16((int16_t)(adjusted_weight));
        const __m256i adjusted_weight_int32 = _mm256_set1_epi32((int32_t)(adjusted_weight));

        for (i = 0; i < block_height; i++) {
            for (j = 0; j < block_width; j += 8) {
                k = i * y_pre_stride + j;

                //y_count[k] += adjusted_weight;
                __m128i count_array = _mm_loadu_si128((__m128i *)(y_count + k));
#if FIX_TEMPORAL_FILTER_PLANEWISE
                count_array = _mm_add_epi16(count_array, adjusted_weight_int16);
#else
                count_array = _mm_adds_epu16(count_array, adjusted_weight_int16);
#endif
                _mm_storeu_si128((__m128i *)(y_count + k), count_array);

                //y_accum[k] += adjusted_weight * pixel_value;
                __m256i accumulator_array = _mm256_loadu_si256((__m256i *)(y_accum + k));
                __m128i frame2_array = _mm_loadu_si128((__m128i *)(y_pre + k));
                __m256i frame2_array_u32 = _mm256_cvtepi16_epi32(frame2_array);
                frame2_array_u32 = _mm256_mullo_epi32(frame2_array_u32, adjusted_weight_int32);

                accumulator_array = _mm256_add_epi32(accumulator_array, frame2_array_u32);
                _mm256_storeu_si256((__m256i *)(y_accum + k), accumulator_array);
            }
        }
    }
}

#if FIXED_POINTS_PLANEWISE
uint32_t calculate_squared_errors_sum_no_div_avx2(const uint8_t *s, int s_stride, const uint8_t *p,
    int p_stride, unsigned int w, unsigned int h) {
    assert(w % 16 == 0 && "block width must be multiple of 16");
    unsigned int i, j;

    __m256i sum = _mm256_setzero_si256();
    __m256i s_16, p_16, dif;

    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j += 16) {
            s_16 = _mm256_cvtepu8_epi16(_mm_loadu_si128((const __m128i *)(s + i * s_stride + j)));
            p_16 = _mm256_cvtepu8_epi16(_mm_loadu_si128((const __m128i *)(p + i * p_stride + j)));

            dif = _mm256_sub_epi16(s_16, p_16);
            sum = _mm256_add_epi32(sum, _mm256_madd_epi16(dif, dif));
        }
    }
    __m128i sum_128 = _mm_add_epi32(_mm256_castsi256_si128(sum),
        _mm256_extracti128_si256(sum, 0x1));
    sum_128 = _mm_hadd_epi32(sum_128, sum_128);
    sum_128 = _mm_hadd_epi32(sum_128, sum_128);

    return _mm_cvtsi128_si32(sum_128);
}
/*This function return 2 separate squared errors for two block 8xh, return value is stored in output array*/
void calculate_squared_errors_sum_2x8xh_no_div_avx2(const uint8_t *s, int s_stride, const uint8_t *p,
    int p_stride, unsigned int h, uint32_t *output ) {
    unsigned int i;

    __m256i sum = _mm256_setzero_si256();
    __m256i s_16, p_16, dif;

    for (i = 0; i < h; i++) {
        s_16 = _mm256_cvtepu8_epi16(_mm_loadu_si128((const __m128i *)(s + i * s_stride)));
        p_16 = _mm256_cvtepu8_epi16(_mm_loadu_si128((const __m128i *)(p + i * p_stride)));

        dif = _mm256_sub_epi16(s_16, p_16);
        sum = _mm256_add_epi32(sum, _mm256_madd_epi16(dif, dif));
    }
    sum = _mm256_hadd_epi32(sum, sum);
    sum = _mm256_hadd_epi32(sum, sum);

    output[0] = _mm_cvtsi128_si32(_mm256_castsi256_si128(sum));
    output[1] = _mm_cvtsi128_si32(_mm256_extracti128_si256(sum, 0x1));
}

uint32_t calculate_squared_errors_sum_no_div_highbd_avx2(const uint16_t *s, int s_stride,
                                                         const uint16_t *p, int p_stride,
                                                         unsigned int w, unsigned int h,
                                                         int shift_factor) {
    assert(w % 16 == 0 && "block width must be multiple of 16");
    unsigned int i, j;

    __m256i sum = _mm256_setzero_si256();
    __m256i s_16, p_16, dif;

    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j += 16) {
            s_16 = _mm256_loadu_si256((const __m256i *)(s + i * s_stride + j));
            p_16 = _mm256_loadu_si256((const __m256i *)(p + i * p_stride + j));

            dif = _mm256_sub_epi16(s_16, p_16);
            sum = _mm256_add_epi32(sum, _mm256_madd_epi16(dif, dif));
        }
    }
    __m128i sum_128 = _mm_add_epi32(_mm256_castsi256_si128(sum),
                                    _mm256_extracti128_si256(sum, 0x1));
    sum_128         = _mm_hadd_epi32(sum_128, sum_128);
    sum_128         = _mm_hadd_epi32(sum_128, sum_128);

    return _mm_cvtsi128_si32(sum_128) >> shift_factor;
}

/*This function return 2 separate squared errors for two block 8xh, return value is stored in output array*/
void calculate_squared_errors_sum_2x8xh_no_div_highbd_avx2(const uint16_t *s, int s_stride,
                                                            const uint16_t *p, int p_stride,
                                                            unsigned int h, int shift_factor,
                                                            uint32_t *output) {
    unsigned int i;

    __m256i sum = _mm256_setzero_si256();
    __m256i s_16, p_16, dif;

    for (i = 0; i < h; i++) {
        s_16 = _mm256_loadu_si256((const __m256i *)(s + i * s_stride));
        p_16 = _mm256_loadu_si256((const __m256i *)(p + i * p_stride));

        dif = _mm256_sub_epi16(s_16, p_16);
        sum = _mm256_add_epi32(sum, _mm256_madd_epi16(dif, dif));
    }
    sum = _mm256_hadd_epi32(sum, sum);
    sum = _mm256_hadd_epi32(sum, sum);

    output[0] = _mm_cvtsi128_si32(_mm256_castsi256_si128(sum)) >> shift_factor;
    output[1] = _mm_cvtsi128_si32(_mm256_extracti128_si256(sum, 0x1)) >> shift_factor;
}

static void svt_av1_apply_temporal_filter_planewise_medium_partial_avx2(
    struct MeContext *context_ptr, const uint8_t *y_src, int y_src_stride, const uint8_t *y_pre,
    int y_pre_stride, unsigned int block_width, unsigned int block_height, uint32_t *y_accum, uint16_t *y_count,
    const uint32_t tf_decay_factor,uint32_t luma_window_error_quad_fp8[4], int is_chroma) {
    unsigned int i, j, k, subblock_idx;

    int32_t idx_32x32 = context_ptr->tf_block_col + context_ptr->tf_block_row * 2;
    uint32_t distance_threshold_fp16 = AOMMAX((context_ptr->min_frame_size<<16) / 10, 1<<16);

    //Calculation for every quarter
    uint32_t  d_factor_fp8[4];
    uint32_t block_error_fp8[4];
    uint32_t  chroma_window_error_quad_fp8[4];
    uint32_t *window_error_quad_fp8 = is_chroma?chroma_window_error_quad_fp8:luma_window_error_quad_fp8;

    if (context_ptr->tf_32x32_block_split_flag[idx_32x32]) {
        for (i = 0; i < 4; ++i) {
            int32_t col = context_ptr->tf_16x16_mv_x[idx_32x32 * 4 + i];
            int32_t row = context_ptr->tf_16x16_mv_y[idx_32x32 * 4 + i];
            //const float  distance = sqrtf((float)col*col + row*row);
            uint32_t distance_fp4 = sqrt_fast(((uint32_t)(col*col + row*row))<<8);
            d_factor_fp8[i] = AOMMAX((distance_fp4<<12) / (distance_threshold_fp16>>8), 1<<8);
            FP_ASSERT(context_ptr->tf_16x16_block_error[idx_32x32 * 4 + i] < ((uint64_t)1<<31));
            //block_error[i] = (double)context_ptr->tf_16x16_block_error[idx_32x32 * 4 + i] / 256;
            block_error_fp8[i] = (uint32_t)(context_ptr->tf_16x16_block_error[idx_32x32 * 4 + i]);
        }
    }
    else {
        int32_t col = context_ptr->tf_32x32_mv_x[idx_32x32];
        int32_t row = context_ptr->tf_32x32_mv_y[idx_32x32];

        uint32_t distance_fp4 = sqrt_fast(((uint32_t)(col*col + row*row))<<8);
        //d_factor[0] = d_factor[1] = d_factor[2] = d_factor[3] = AOMMAX(distance / distance_threshold, 1);
        d_factor_fp8[0] = d_factor_fp8[1] = d_factor_fp8[2] = d_factor_fp8[3]
            = AOMMAX((distance_fp4<<12) / (distance_threshold_fp16>>8), 1<<8);
        FP_ASSERT(context_ptr->tf_32x32_block_error[idx_32x32] < ((uint64_t)1<<30));
        //block_error[0] = block_error[1] = block_error[2] = block_error[3] = (double)context_ptr->tf_32x32_block_error[idx_32x32] / 1024;
        block_error_fp8[0] = block_error_fp8[1] = block_error_fp8[2] = block_error_fp8[3]
            = (uint32_t)(context_ptr->tf_32x32_block_error[idx_32x32]>>2);
    }

    if (block_width == 32) {
            ////(((sum<<4) /bw_half)<<4)/bh_half;
            ////(((sum<<4) /16)<<4)/16;
            ////sum;
            window_error_quad_fp8[0] = calculate_squared_errors_sum_no_div_avx2(y_src, y_src_stride,
                                                                y_pre, y_pre_stride,
                                                                16, 16);
            window_error_quad_fp8[1] = calculate_squared_errors_sum_no_div_avx2(y_src + 16, y_src_stride,
                                                                y_pre + 16, y_pre_stride,
                                                                16, 16);
            window_error_quad_fp8[2] = calculate_squared_errors_sum_no_div_avx2(y_src + y_src_stride * 16,
                                                                y_src_stride, y_pre + y_pre_stride * 16,
                                                                y_pre_stride, 16, 16);
            window_error_quad_fp8[3] = calculate_squared_errors_sum_no_div_avx2(y_src + y_src_stride * 16 + 16,
                                                                y_src_stride, y_pre + y_pre_stride * 16 + 16,
                                                                y_pre_stride, 16, 16);
    }
    else {//block_width == 16
        //(((sum<<4) /bw_half)<<4)/bh_half;
        //(((sum<<4) /8)<<4)/8;
        //sum<<2;
        calculate_squared_errors_sum_2x8xh_no_div_avx2(
            y_src, y_src_stride, y_pre, y_pre_stride, 8, window_error_quad_fp8);

        calculate_squared_errors_sum_2x8xh_no_div_avx2(y_src + y_src_stride * 8,
                                                        y_src_stride,
                                                        y_pre + y_pre_stride * 8,
                                                        y_pre_stride,
                                                        8,
                                                        &window_error_quad_fp8[2]);
        window_error_quad_fp8[0] <<= 2;
        window_error_quad_fp8[1] <<= 2;
        window_error_quad_fp8[2] <<= 2;
        window_error_quad_fp8[3] <<= 2;
    }

    if (is_chroma) {
        for (i = 0; i < 4; ++i) {
            FP_ASSERT(((int64_t)window_error_quad_fp8[i] * 5 + luma_window_error_quad_fp8[i]) < ((int64_t)1<<31));
            window_error_quad_fp8[i] = (window_error_quad_fp8[i] * 5 + luma_window_error_quad_fp8[i])/6;
        }
    }

    __m128i adjusted_weight_int16[4];
    __m256i adjusted_weight_int32[4];

    for (subblock_idx = 0; subblock_idx < 4; subblock_idx++) {
        uint32_t combined_error_fp8 = (window_error_quad_fp8[subblock_idx] * TF_WINDOW_BLOCK_BALANCE_WEIGHT + block_error_fp8[subblock_idx]) /
                                (TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1);

        uint32_t avg_err_fp10 = ((combined_error_fp8>>3) * (d_factor_fp8[subblock_idx]>>3));
        FP_ASSERT((((int64_t)combined_error_fp8>>3) * (d_factor_fp8[subblock_idx]>>3)) < ((int64_t)1<<31));
        uint32_t scaled_diff16 = AOMMIN(/*((16*avg_err)<<8)*/(avg_err_fp10) / (tf_decay_factor>>10), 7*16);
        uint32_t adjusted_weight = (expf_tab_fp16[scaled_diff16] * TF_WEIGHT_SCALE) >> 16;

        adjusted_weight_int16[subblock_idx] = _mm_set1_epi16((int16_t)(adjusted_weight));
        adjusted_weight_int32[subblock_idx] = _mm256_set1_epi32((int32_t)(adjusted_weight));
    }

    for (i = 0; i < block_height; i++) {
        const int subblock_idx_h = (i >= block_height / 2) * 2;
        for (j = 0; j < block_width; j += 8) {
            k  = i * y_pre_stride + j;

                //y_count[k] += adjusted_weight;
                __m128i count_array = _mm_loadu_si128((__m128i *)(y_count + k));
                count_array = _mm_add_epi16(count_array, adjusted_weight_int16[subblock_idx_h + (j >= block_width / 2)]);
                _mm_storeu_si128((__m128i *)(y_count + k), count_array);

                //y_accum[k] += adjusted_weight * pixel_value;
                __m256i accumulator_array = _mm256_loadu_si256((__m256i *)(y_accum + k));
                __m128i frame2_array = _mm_loadl_epi64((__m128i *)(y_pre + k));
                frame2_array = _mm_cvtepu8_epi16(frame2_array);
                __m256i frame2_array_u32 = _mm256_cvtepi16_epi32(frame2_array);
                frame2_array_u32 = _mm256_mullo_epi32(frame2_array_u32, adjusted_weight_int32[subblock_idx_h + (j >= block_width / 2)]);

                accumulator_array = _mm256_add_epi32(accumulator_array, frame2_array_u32);
                _mm256_storeu_si256((__m256i *)(y_accum + k), accumulator_array);

        }
    }
}

void svt_av1_apply_temporal_filter_planewise_medium_avx2(
    struct MeContext *context_ptr, const uint8_t *y_src, int y_src_stride, const uint8_t *y_pre,
    int y_pre_stride, const uint8_t *u_src, const uint8_t *v_src, int uv_src_stride,
    const uint8_t *u_pre, const uint8_t *v_pre, int uv_pre_stride, unsigned int block_width,
    unsigned int block_height, int ss_x, int ss_y, uint32_t *y_accum, uint16_t *y_count, uint32_t *u_accum,
    uint16_t *u_count, uint32_t *v_accum, uint16_t *v_count)
{
    uint32_t luma_window_error_quad_fp8[4];

    svt_av1_apply_temporal_filter_planewise_medium_partial_avx2(context_ptr,
        y_src,
        y_src_stride,
        y_pre,
        y_pre_stride,
        (unsigned int)block_width,
        (unsigned int)block_height,
        y_accum,
        y_count,
        context_ptr->tf_decay_factor_fp16[C_Y],
        luma_window_error_quad_fp8,
        0);

    if (context_ptr->tf_chroma) {
        svt_av1_apply_temporal_filter_planewise_medium_partial_avx2(context_ptr,
            u_src,
            uv_src_stride,
            u_pre,
            uv_pre_stride,
            (unsigned int)block_width >> ss_x,
            (unsigned int)block_height >> ss_y,
            u_accum,
            u_count,
            context_ptr->tf_decay_factor_fp16[C_U],
            luma_window_error_quad_fp8,
            1);

        svt_av1_apply_temporal_filter_planewise_medium_partial_avx2(context_ptr,
            v_src,
            uv_src_stride,
            v_pre,
            uv_pre_stride,
            (unsigned int)block_width >> ss_x,
            (unsigned int)block_height >> ss_y,
            v_accum,
            v_count,
            context_ptr->tf_decay_factor_fp16[C_V],
            luma_window_error_quad_fp8,
            1);
    }
}

static void svt_av1_apply_temporal_filter_planewise_medium_hbd_partial_avx2(
    struct MeContext *context_ptr, const uint16_t *y_src, int y_src_stride, const uint16_t *y_pre,
    int y_pre_stride, unsigned int block_width,
    unsigned int block_height, uint32_t *y_accum, uint16_t *y_count, const uint32_t tf_decay_factor,
    uint32_t luma_window_error_quad_fp8[4], int is_chroma, uint32_t encoder_bit_depth) {
    unsigned int i, j, k, subblock_idx;

    int32_t idx_32x32 = context_ptr->tf_block_col + context_ptr->tf_block_row * 2;
    int shift_factor   = ((encoder_bit_depth - 8) * 2);
    uint32_t distance_threshold_fp16 = AOMMAX((context_ptr->min_frame_size<<16) / 10, 1<<16); //TODO Change to FP8

    //Calculation for every quarter
    uint32_t  d_factor_fp8[4];
    uint32_t block_error_fp8[4];
    uint32_t  chroma_window_error_quad_fp8[4];
    uint32_t *window_error_quad_fp8 = is_chroma?chroma_window_error_quad_fp8:luma_window_error_quad_fp8;

    if (context_ptr->tf_32x32_block_split_flag[idx_32x32]) {
        for (i = 0; i < 4; ++i) {
            int32_t col = context_ptr->tf_16x16_mv_x[idx_32x32 * 4 + i];
            int32_t row = context_ptr->tf_16x16_mv_y[idx_32x32 * 4 + i];
            //const float  distance = sqrtf((float)col*col + row*row);
            uint32_t distance_fp4 = sqrt_fast(((uint32_t)(col*col + row*row))<<8);
            d_factor_fp8[i] = AOMMAX((distance_fp4<<12) / (distance_threshold_fp16>>8), 1<<8);
            FP_ASSERT(context_ptr->tf_16x16_block_error[idx_32x32 * 4 + i] < ((uint64_t)1<<35));
            //block_error[i] = (double)context_ptr->tf_16x16_block_error[idx_32x32 * 4 + i] / 256;
            block_error_fp8[i] = (uint32_t)(context_ptr->tf_16x16_block_error[idx_32x32 * 4 + i] >>4);
        }
    }
    else {
        int32_t col = context_ptr->tf_32x32_mv_x[idx_32x32];
        int32_t row = context_ptr->tf_32x32_mv_y[idx_32x32];

        uint32_t distance_fp4 = sqrt_fast(((uint32_t)(col*col + row*row))<<8);
        //d_factor[0] = d_factor[1] = d_factor[2] = d_factor[3] = AOMMAX(distance / distance_threshold, 1);
        d_factor_fp8[0] = d_factor_fp8[1] = d_factor_fp8[2] = d_factor_fp8[3]
            = AOMMAX((distance_fp4<<12) / (distance_threshold_fp16>>8), 1<<8);
        FP_ASSERT(context_ptr->tf_32x32_block_error[idx_32x32] < ((uint64_t)1<<35));
        //block_error[0] = block_error[1] = block_error[2] = block_error[3] = (double)context_ptr->tf_32x32_block_error[idx_32x32] / 1024;
        block_error_fp8[0] = block_error_fp8[1] = block_error_fp8[2] = block_error_fp8[3]
            = (uint32_t)(context_ptr->tf_32x32_block_error[idx_32x32]>> 6);
    }

    if (block_width == 32) {
            ////(((sum<<4) /bw_half)<<4)/bh_half;
            ////(((sum<<4) /16)<<4)/16;
            ////sum;
            window_error_quad_fp8[0] = calculate_squared_errors_sum_no_div_highbd_avx2(y_src, y_src_stride,
                                                                y_pre, y_pre_stride,
                                                                16, 16, shift_factor);
            window_error_quad_fp8[1] = calculate_squared_errors_sum_no_div_highbd_avx2(y_src + 16, y_src_stride,
                                                                y_pre + 16, y_pre_stride,
                                                                16, 16, shift_factor);
            window_error_quad_fp8[2] = calculate_squared_errors_sum_no_div_highbd_avx2(y_src + y_src_stride * 16,
                                                                y_src_stride, y_pre + y_pre_stride * 16,
                                                                y_pre_stride, 16, 16, shift_factor);
            window_error_quad_fp8[3] = calculate_squared_errors_sum_no_div_highbd_avx2(y_src + y_src_stride * 16 + 16,
                                                                y_src_stride, y_pre + y_pre_stride * 16 + 16,
                                                                y_pre_stride, 16, 16, shift_factor);

    }
    else {//block_width == 16
        //(((sum<<4) /bw_half)<<4)/bh_half;
        //(((sum<<4) /8)<<4)/8;
        //sum<<2;
        calculate_squared_errors_sum_2x8xh_no_div_highbd_avx2(
            y_src, y_src_stride, y_pre, y_pre_stride, 8, shift_factor, window_error_quad_fp8);

        calculate_squared_errors_sum_2x8xh_no_div_highbd_avx2(y_src + y_src_stride * 8,
                                                        y_src_stride,
                                                        y_pre + y_pre_stride * 8,
                                                        y_pre_stride,
                                                        8, shift_factor,
                                                        &window_error_quad_fp8[2]);
        window_error_quad_fp8[0] <<= 2;
        window_error_quad_fp8[1] <<= 2;
        window_error_quad_fp8[2] <<= 2;
        window_error_quad_fp8[3] <<= 2;
    }

    if (is_chroma) {
        for (i = 0; i < 4; ++i) {
            FP_ASSERT(((int64_t)window_error_quad_fp8[i] * 5 + luma_window_error_quad_fp8[i]) < ((int64_t)1<<31));
            window_error_quad_fp8[i] = (window_error_quad_fp8[i] * 5 + luma_window_error_quad_fp8[i])/6;
        }
    }

    __m128i adjusted_weight_int16[4];
    __m256i adjusted_weight_int32[4];

    for (subblock_idx = 0; subblock_idx < 4; subblock_idx++) {
        uint32_t combined_error_fp8 = (window_error_quad_fp8[subblock_idx] * TF_WINDOW_BLOCK_BALANCE_WEIGHT + block_error_fp8[subblock_idx]) /
                                (TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1);

        uint32_t avg_err_fp10 = ((combined_error_fp8>>3) * (d_factor_fp8[subblock_idx]>>3));
        FP_ASSERT((((int64_t)combined_error_fp8>>3) * (d_factor_fp8[subblock_idx]>>3)) < ((int64_t)1<<31));
        uint32_t scaled_diff16 = AOMMIN(/*((16*avg_err)<<8)*/(avg_err_fp10) / (tf_decay_factor>>10), 7*16);
        int adjusted_weight = (expf_tab_fp16[scaled_diff16] * TF_WEIGHT_SCALE) >> 16;

        adjusted_weight_int16[subblock_idx] = _mm_set1_epi16((int16_t)(adjusted_weight));
        adjusted_weight_int32[subblock_idx] = _mm256_set1_epi32((int32_t)(adjusted_weight));
    }

    for (i = 0; i < block_height; i++) {
        const int subblock_idx_h = (i >= block_height / 2) * 2;
        for (j = 0; j < block_width; j += 8) {
            k = i * y_pre_stride + j;

            //y_count[k] += adjusted_weight;
            __m128i count_array = _mm_loadu_si128((__m128i *)(y_count + k));
            count_array = _mm_add_epi16(count_array, adjusted_weight_int16[subblock_idx_h + (j >= block_width / 2)]);
            _mm_storeu_si128((__m128i *)(y_count + k), count_array);

            //y_accum[k] += adjusted_weight * pixel_value;
            __m256i accumulator_array = _mm256_loadu_si256((__m256i *)(y_accum + k));
            __m128i frame2_array      = _mm_loadu_si128((__m128i *)(y_pre + k));
            __m256i frame2_array_u32  = _mm256_cvtepi16_epi32(frame2_array);
            frame2_array_u32          = _mm256_mullo_epi32(frame2_array_u32, adjusted_weight_int32[subblock_idx_h + (j >= block_width / 2)]);

            accumulator_array = _mm256_add_epi32(accumulator_array, frame2_array_u32);
            _mm256_storeu_si256((__m256i *)(y_accum + k), accumulator_array);
        }
    }
}

void svt_av1_apply_temporal_filter_planewise_medium_hbd_avx2(
    struct MeContext *context_ptr, const uint16_t *y_src, int y_src_stride, const uint16_t *y_pre,
    int y_pre_stride, const uint16_t *u_src, const uint16_t *v_src, int uv_src_stride,
    const uint16_t *u_pre, const uint16_t *v_pre, int uv_pre_stride, unsigned int block_width,
    unsigned int block_height, int ss_x, int ss_y, uint32_t *y_accum, uint16_t *y_count, uint32_t *u_accum,
    uint16_t *u_count, uint32_t *v_accum, uint16_t *v_count, uint32_t encoder_bit_depth)
{
    uint32_t luma_window_error_quad_fp8[4];

    svt_av1_apply_temporal_filter_planewise_medium_hbd_partial_avx2(context_ptr,
                                                        y_src,
                                                        y_src_stride,
                                                        y_pre,
                                                        y_pre_stride,
                                                        (unsigned int)block_width,
                                                        (unsigned int)block_height,
                                                        y_accum,
                                                        y_count,
                                                        context_ptr->tf_decay_factor_fp16[C_Y],
                                                        luma_window_error_quad_fp8,
                                                        0,
                                                        encoder_bit_depth);
    if (context_ptr->tf_chroma) {
        svt_av1_apply_temporal_filter_planewise_medium_hbd_partial_avx2(context_ptr,
                                                            u_src,
                                                            uv_src_stride,
                                                            u_pre,
                                                            uv_pre_stride,
                                                            (unsigned int)block_width >> ss_x,
                                                            (unsigned int)block_height >> ss_y,
                                                            u_accum,
                                                            u_count,
                                                            context_ptr->tf_decay_factor_fp16[C_U],
                                                            luma_window_error_quad_fp8,
                                                            1,
                                                            encoder_bit_depth);

        svt_av1_apply_temporal_filter_planewise_medium_hbd_partial_avx2(context_ptr,
                                                            v_src,
                                                            uv_src_stride,
                                                            v_pre,
                                                            uv_pre_stride,
                                                            (unsigned int)block_width >> ss_x,
                                                            (unsigned int)block_height >> ss_y,
                                                            v_accum,
                                                            v_count,
                                                            context_ptr->tf_decay_factor_fp16[C_V],
                                                            luma_window_error_quad_fp8,
                                                            1,
                                                            encoder_bit_depth);
    }
}
#endif /*FIXED_POINTS_PLANEWISE*/
#endif /*OPT_TFILTER*/

static INLINE __m256i __m256i_div_epi32(const __m256i *a, const __m256i *b) {
    __m256 d_f = _mm256_div_ps(_mm256_cvtepi32_ps(*a), _mm256_cvtepi32_ps(*b));
    return _mm256_cvtps_epi32(_mm256_floor_ps(d_f));
}

static void process_block_lbd_avx2(int h, int w, uint8_t *buff_lbd_start, uint32_t *accum,
                                   uint16_t *count, uint32_t stride) {
    int i, j, k;
    int pos = 0;
    for (i = 0, k = 0; i < h; i++) {
        for (j = 0; j < w; j += 16, k += 16) {
            //buff_lbd_start[pos] = (uint8_t)OD_DIVU(accum[k] + (count[k] >> 1), count[k]);
            //buff_lbd_start[pos] = (uint8_t)((accum[k] + (count[k] >> 1))/ count[k]);
            __m256i accum_a = _mm256_loadu_si256((__m256i *)(accum + k));
            __m256i accum_b = _mm256_loadu_si256((__m256i *)(accum + k + 8));
            __m256i count_a = _mm256_cvtepi16_epi32(_mm_loadu_si128((__m128i *)(count + k)));
            __m256i count_b = _mm256_cvtepi16_epi32(_mm_loadu_si128((__m128i *)(count + k + 8)));

            //accum[k] + (count[k] >> 1)
            __m256i tmp_a = _mm256_add_epi32(accum_a, _mm256_srli_epi32(count_a, 1));
            __m256i tmp_b = _mm256_add_epi32(accum_b, _mm256_srli_epi32(count_b, 1));

            //accum[k] + (count[k] >> 1))/ count[k]
            tmp_a          = __m256i_div_epi32(&tmp_a, &count_a);
            tmp_b          = __m256i_div_epi32(&tmp_b, &count_b);
            __m256i tmp_ab = _mm256_packs_epi32(tmp_a, tmp_b);

            tmp_ab = _mm256_permutevar8x32_epi32(tmp_ab, _mm256_setr_epi32(0, 1, 4, 5, 2, 3, 6, 7));

            __m128i tmp_ab_128 = _mm_packus_epi16(_mm256_castsi256_si128(tmp_ab),
                                                  _mm256_extracti128_si256(tmp_ab, 0x1));

            _mm_storeu_si128((__m128i *)(buff_lbd_start + pos), tmp_ab_128);

            pos += 16;
        }
        pos += stride;
    }
}

static void process_block_hbd_avx2(int h, int w, uint16_t *buff_hbd_start, uint32_t *accum,
                                   uint16_t *count, uint32_t stride) {
    int i, j, k;
    int pos = 0;
    for (i = 0, k = 0; i < h; i++) {
        for (j = 0; j < w; j += 16, k += 16) {
            //buff_hbd_start[pos] = (uint8_t)OD_DIVU(accum[k] + (count[k] >> 1), count[k]);
            //buff_hbd_start[pos] = (uint8_t)((accum[k] + (count[k] >> 1))/ count[k]);
            __m256i accum_a = _mm256_loadu_si256((__m256i *)(accum + k));
            __m256i accum_b = _mm256_loadu_si256((__m256i *)(accum + k + 8));
            __m256i count_a = _mm256_cvtepi16_epi32(_mm_loadu_si128((__m128i *)(count + k)));
            __m256i count_b = _mm256_cvtepi16_epi32(_mm_loadu_si128((__m128i *)(count + k + 8)));

            //accum[k] + (count[k] >> 1)
            __m256i tmp_a = _mm256_add_epi32(accum_a, _mm256_srli_epi32(count_a, 1));
            __m256i tmp_b = _mm256_add_epi32(accum_b, _mm256_srli_epi32(count_b, 1));

            //accum[k] + (count[k] >> 1))/ count[k]
            tmp_a          = __m256i_div_epi32(&tmp_a, &count_a);
            tmp_b          = __m256i_div_epi32(&tmp_b, &count_b);
            __m256i tmp_ab = _mm256_packs_epi32(tmp_a, tmp_b);

            tmp_ab = _mm256_permutevar8x32_epi32(tmp_ab, _mm256_setr_epi32(0, 1, 4, 5, 2, 3, 6, 7));

            _mm256_storeu_si256((__m256i *)(buff_hbd_start + pos), tmp_ab);

            pos += 16;
        }
        pos += stride;
    }
}

void get_final_filtered_pixels_avx2(MeContext *context_ptr, EbByte *src_center_ptr_start,
                                    uint16_t **altref_buffer_highbd_start, uint32_t **accum,
                                    uint16_t **count, const uint32_t *stride,
                                    int blk_y_src_offset, int blk_ch_src_offset,
                                    uint16_t blk_width_ch, uint16_t blk_height_ch,
                                    EbBool is_highbd) {
    assert(blk_width_ch % 16 == 0);
    assert(BW % 16 == 0);

    if (!is_highbd) {
         //Process luma
        process_block_lbd_avx2(BH, BW, &src_center_ptr_start[C_Y][blk_y_src_offset], accum[C_Y], count[C_Y], stride[C_Y]- BW);
        // Process chroma
        if (context_ptr->tf_chroma) {
            process_block_lbd_avx2(blk_height_ch, blk_width_ch, &src_center_ptr_start[C_U][blk_ch_src_offset], accum[C_U], count[C_U], stride[C_U]- blk_width_ch);
            process_block_lbd_avx2(blk_height_ch, blk_width_ch, &src_center_ptr_start[C_V][blk_ch_src_offset], accum[C_V], count[C_V], stride[C_V]- blk_width_ch);
        }
    } else {
        // Process luma
        process_block_hbd_avx2(BH, BW, &altref_buffer_highbd_start[C_Y][blk_y_src_offset], accum[C_Y], count[C_Y], stride[C_Y] - BW);
        // Process chroma
        if (context_ptr->tf_chroma) {
            process_block_hbd_avx2(blk_height_ch, blk_width_ch, &altref_buffer_highbd_start[C_U][blk_ch_src_offset], accum[C_U], count[C_U], stride[C_U]- blk_width_ch);
            process_block_hbd_avx2(blk_height_ch, blk_width_ch, &altref_buffer_highbd_start[C_V][blk_ch_src_offset], accum[C_V], count[C_V], stride[C_V]- blk_width_ch);
        }
    }
}

static void apply_filtering_central_loop_lbd(uint16_t w, uint16_t h, uint8_t *src,
                                             uint16_t src_stride, uint32_t *accum,
                                             uint16_t *count) {
    assert(w % 8 == 0);

    __m256i modifier       = _mm256_set1_epi32(TF_PLANEWISE_FILTER_WEIGHT_SCALE);
    __m128i modifier_epi16 = _mm_set1_epi16(TF_PLANEWISE_FILTER_WEIGHT_SCALE);

    for (uint16_t k = 0, i = 0; i < h; i++) {
        for (uint16_t j = 0; j < w; j += 8) {
            __m256i src_ = _mm256_cvtepu16_epi32(
                _mm_cvtepu8_epi16(_mm_loadl_epi64((__m128i *)(src + i * src_stride + j))));
            _mm256_storeu_si256((__m256i *)(accum + k), _mm256_mullo_epi32(modifier, src_));
            _mm_storeu_si128((__m128i *)(count + k), modifier_epi16);
            k += 8;
        }
    }
}

static void apply_filtering_central_loop_hbd(uint16_t w, uint16_t h, uint16_t *src,
                                             uint16_t src_stride, uint32_t *accum,
                                             uint16_t *count) {
    assert(w % 8 == 0);

    __m256i modifier       = _mm256_set1_epi32(TF_PLANEWISE_FILTER_WEIGHT_SCALE);
    __m128i modifier_epi16 = _mm_set1_epi16(TF_PLANEWISE_FILTER_WEIGHT_SCALE);

    for (uint16_t k = 0, i = 0; i < h; i++) {
        for (uint16_t j = 0; j < w; j += 8) {
            __m256i src_ = _mm256_cvtepu16_epi32(
                _mm_loadu_si128((__m128i *)(src + i * src_stride + j)));
            _mm256_storeu_si256((__m256i *)(accum + k), _mm256_mullo_epi32(modifier, src_));
            _mm_storeu_si128((__m128i *)(count + k), modifier_epi16);
            k += 8;
        }
    }
}

// Apply filtering to the central picture
void apply_filtering_central_avx2(MeContext *          context_ptr,
                                  EbPictureBufferDesc *input_picture_ptr_central,
                                  EbByte *src, uint32_t **accum, uint16_t **count,
                                  uint16_t blk_width, uint16_t blk_height, uint32_t ss_x,
                                  uint32_t ss_y) {
    uint16_t src_stride_y = input_picture_ptr_central->stride_y;

    // Luma
    apply_filtering_central_loop_lbd(
        blk_width, blk_height, src[C_Y], src_stride_y, accum[C_Y], count[C_Y]);

    // Chroma
    if (context_ptr->tf_chroma) {
        uint16_t blk_height_ch = blk_height >> ss_y;
        uint16_t blk_width_ch  = blk_width >> ss_x;
        uint16_t src_stride_ch = src_stride_y >> ss_x;
        apply_filtering_central_loop_lbd(
            blk_width_ch, blk_height_ch, src[C_U], src_stride_ch, accum[C_U], count[C_U]);
        apply_filtering_central_loop_lbd(
            blk_width_ch, blk_height_ch, src[C_V], src_stride_ch, accum[C_V], count[C_V]);
    }
}

// Apply filtering to the central picture
void apply_filtering_central_highbd_avx2(MeContext *          context_ptr,
                                         EbPictureBufferDesc *input_picture_ptr_central,
                                         uint16_t **src_16bit, uint32_t **accum,
                                         uint16_t **count, uint16_t blk_width,
                                         uint16_t blk_height, uint32_t ss_x, uint32_t ss_y) {
    uint16_t src_stride_y = input_picture_ptr_central->stride_y;

    // Luma
    apply_filtering_central_loop_hbd(
        blk_width, blk_height, src_16bit[C_Y], src_stride_y, accum[C_Y], count[C_Y]);

    // Chroma
    if (context_ptr->tf_chroma) {
        uint16_t blk_height_ch = blk_height >> ss_y;
        uint16_t blk_width_ch  = blk_width >> ss_x;
        uint16_t src_stride_ch = src_stride_y >> ss_x;
        apply_filtering_central_loop_hbd(
            blk_width_ch, blk_height_ch, src_16bit[C_U], src_stride_ch, accum[C_U], count[C_U]);
        apply_filtering_central_loop_hbd(
            blk_width_ch, blk_height_ch, src_16bit[C_V], src_stride_ch, accum[C_V], count[C_V]);
    }
}
