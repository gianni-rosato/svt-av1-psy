/*
 * Copyright (c) 2019, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <assert.h>
#include <immintrin.h> /* AVX2 */

#include "EbDefinitions.h"
#include "./../ASM_SSE4_1/EbTemporalFiltering_constants.h"


#define SSE_STRIDE (BW + 2)

DECLARE_ALIGNED(32, static const uint32_t, sse_bytemask[4][8]) = {
  { 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0, 0, 0 },
  { 0, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0, 0 },
  { 0, 0, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0 },
  { 0, 0, 0, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF }
};

DECLARE_ALIGNED(32, static const uint8_t, shufflemask_16b[2][16]) = {
  { 0, 1, 0, 1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 },
  { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 10, 11, 10, 11 }
};

static AOM_FORCE_INLINE void get_squared_error_16x16_avx2(
    const uint8_t *frame1, const unsigned int stride, const uint8_t *frame2,
    const unsigned int stride2, const int block_width, const int block_height,
    uint16_t *frame_sse, const unsigned int sse_stride) {
    (void)block_width;
    const uint8_t *src1 = frame1;
    const uint8_t *src2 = frame2;
    uint16_t *dst = frame_sse;
    for (int i = 0; i < block_height; i++) {
        __m128i vf1_128, vf2_128;
        __m256i vf1, vf2, vdiff1, vsqdiff1;

        vf1_128 = _mm_loadu_si128((__m128i *)(src1));
        vf2_128 = _mm_loadu_si128((__m128i *)(src2));
        vf1 = _mm256_cvtepu8_epi16(vf1_128);
        vf2 = _mm256_cvtepu8_epi16(vf2_128);
        vdiff1 = _mm256_sub_epi16(vf1, vf2);
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
    const unsigned int stride2, const int block_width, const int block_height,
    uint16_t *frame_sse, const unsigned int sse_stride) {
    (void)block_width;
    const uint8_t *src1 = frame1;
    const uint8_t *src2 = frame2;
    uint16_t *dst = frame_sse;
    for (int i = 0; i < block_height; i++) {
        __m256i vsrc1, vsrc2, vmin, vmax, vdiff, vdiff1, vdiff2, vres1, vres2;

        vsrc1 = _mm256_loadu_si256((__m256i *)src1);
        vsrc2 = _mm256_loadu_si256((__m256i *)src2);
        vmax = _mm256_max_epu8(vsrc1, vsrc2);
        vmin = _mm256_min_epu8(vsrc1, vsrc2);
        vdiff = _mm256_subs_epu8(vmax, vmin);

        __m128i vtmp1 = _mm256_castsi256_si128(vdiff);
        __m128i vtmp2 = _mm256_extracti128_si256(vdiff, 1);
        vdiff1 = _mm256_cvtepu8_epi16(vtmp1);
        vdiff2 = _mm256_cvtepu8_epi16(vtmp2);

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

static AOM_FORCE_INLINE __m256i xx_load_and_pad(uint16_t *src, int col,
    int block_width) {
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

static void apply_temporal_filter_planewise(
    const uint8_t *frame1, const unsigned int stride, const uint8_t *frame2,
    const unsigned int stride2, const int block_width, const int block_height,
    const double sigma, const int decay_control, unsigned int *accumulator,
    uint16_t *count, uint16_t *luma_sq_error, uint16_t *chroma_sq_error,
    int plane, int ss_x_shift, int ss_y_shift) {
    assert(TF_PLANEWISE_FILTER_WINDOW_LENGTH == 5);
    assert(((block_width == 32) && (block_height == 32)) ||
        ((block_width == 16) && (block_height == 16)));
    if (plane > PLANE_TYPE_Y) assert(chroma_sq_error != NULL);

    uint32_t acc_5x5_sse[BH][BW];
    const double h = decay_control * (0.7 + log(sigma + 1.0));
    uint16_t *frame_sse =
        (plane == PLANE_TYPE_Y) ? luma_sq_error : chroma_sq_error;

    if (block_width == 32) {
        get_squared_error_32x32_avx2(frame1, stride, frame2, stride2, block_width,
            block_height, frame_sse, SSE_STRIDE);
    }
    else {
        get_squared_error_16x16_avx2(frame1, stride, frame2, stride2, block_width,
            block_height, frame_sse, SSE_STRIDE);
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
            for (int i = 0; i < 5; i++) {
                vsum = _mm256_add_epi32(vsum, vsrc[i]);
            }

            // Push all elements by one element to the top
            for (int i = 0; i < 4; i++) {
                vsrc[i] = vsrc[i + 1];
            }

            // Load next row to the last element
            if (row <= block_width - 4) {
                vsrc[4] = xx_load_and_pad(src, col, block_width);
                src += SSE_STRIDE;
            }
            else {
                vsrc[4] = vsrc[3];
            }

            // Accumulate the sum horizontally
            for (int i = 0; i < 4; i++) {
                acc_5x5_sse[row][col + i] = xx_mask_and_hadd(vsum, i);
            }
        }
    }

    for (int i = 0; i < block_height; i++) {
        for (int j = 0; j < block_width; j++) {
            const int pixel_value = frame2[i * stride2 + j];

            int diff_sse = acc_5x5_sse[i][j];
            int num_ref_pixels =
                TF_PLANEWISE_FILTER_WINDOW_LENGTH * TF_PLANEWISE_FILTER_WINDOW_LENGTH;

             //Filter U-plane and V-plane using Y-plane. This is because motion
             //search is only done on Y-plane, so the information from Y-plane will
             //be more accurate.
            if (plane != PLANE_TYPE_Y) {
                for (int ii = 0; ii < (1 << ss_y_shift); ++ii) {
                    for (int jj = 0; jj < (1 << ss_x_shift); ++jj) {
                        const int yy = (i << ss_y_shift) + ii;  // Y-coord on Y-plane.
                        const int xx = (j << ss_x_shift) + jj;  // X-coord on Y-plane.
                        diff_sse += luma_sq_error[yy * SSE_STRIDE + xx];
                        ++num_ref_pixels;
                    }
                }
            }
            const double scaled_diff =
                AOMMAX(-(double)(diff_sse / num_ref_pixels) / (2 * h * h), -15.0);
            const int adjusted_weight =
                (int)(exp(scaled_diff) * TF_PLANEWISE_FILTER_WEIGHT_SCALE);
            // updated the index
            count[i * stride2 + j] += adjusted_weight;
            accumulator[i * stride2 + j] += adjusted_weight * pixel_value;
        }
    }
}
void svt_av1_apply_temporal_filter_planewise_avx2(
    const uint8_t *y_src, int y_src_stride, const uint8_t *y_pre, int y_pre_stride,
    const uint8_t *u_src, const uint8_t *v_src, int uv_src_stride, const uint8_t *u_pre,
    const uint8_t *v_pre, int uv_pre_stride, unsigned int block_width, unsigned int block_height,
    int ss_x, int ss_y, const double *noise_levels, const int decay_control, uint32_t *y_accum,
    uint16_t *y_count, uint32_t *u_accum, uint16_t *u_count, uint32_t *v_accum, uint16_t *v_count) {
    // Loop variables
    assert(block_width <= BW && "block width too large");
    assert(block_height <= BH && "block height too large");
    assert(block_width % 16 == 0 && "block width must be multiple of 16");
    assert(block_height % 2 == 0 && "block height must be even");
    assert((ss_x == 0 || ss_x == 1) && (ss_y == 0 || ss_y == 1) && "invalid chroma subsampling");

    const int num_planes = 3;
    uint16_t luma_sq_error[SSE_STRIDE * BH];
    uint16_t chroma_sq_error[SSE_STRIDE * BH];

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

        apply_temporal_filter_planewise(
            ref, src_stride, pred, pre_stride, plane_w, plane_h,
            noise_levels[plane], decay_control, accum,
            count, luma_sq_error, chroma_sq_error, plane,
            ss_x_shift, ss_y_shift);
    }
}
