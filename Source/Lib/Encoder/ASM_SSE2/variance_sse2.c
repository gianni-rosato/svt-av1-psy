/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
 */
#include "EbDefinitions.h"
#include <assert.h>
#include <emmintrin.h> // SSE2
#include "aom_dsp_rtcd.h"
#include "EbVariance_SSE2.h"
#include "synonyms.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __cplusplus
}
#endif

// Can handle 128 pixels' diff sum (such as 8x16 or 16x8)
// Slightly faster than variance_final_256_pel_sse2()
// diff sum of 128 pixels can still fit in 16bit integer
static INLINE void variance_final_128_pel_sse2(__m128i vsse, __m128i vsum, unsigned int *const sse,
                                               int *const sum) {
    *sse = add32x4_sse2(vsse);

    vsum = _mm_add_epi16(vsum, _mm_srli_si128(vsum, 8));
    vsum = _mm_add_epi16(vsum, _mm_srli_si128(vsum, 4));
    vsum = _mm_add_epi16(vsum, _mm_srli_si128(vsum, 2));
    *sum = (int16_t)_mm_extract_epi16(vsum, 0);
}

// Can handle 256 pixels' diff sum (such as 16x16)
static INLINE void variance_final_256_pel_sse2(__m128i vsse, __m128i vsum, unsigned int *const sse,
                                               int *const sum) {
    *sse = add32x4_sse2(vsse);

    vsum = _mm_add_epi16(vsum, _mm_srli_si128(vsum, 8));
    vsum = _mm_add_epi16(vsum, _mm_srli_si128(vsum, 4));
    *sum = (int16_t)_mm_extract_epi16(vsum, 0);
    *sum += (int16_t)_mm_extract_epi16(vsum, 1);
}

static INLINE void variance_kernel_sse2(const __m128i src, const __m128i ref, __m128i *const sse,
                                        __m128i *const sum) {
    const __m128i diff = _mm_sub_epi16(src, ref);
    *sse               = _mm_add_epi32(*sse, _mm_madd_epi16(diff, diff));
    *sum               = _mm_add_epi16(*sum, diff);
}

static INLINE void variance4_sse2(const uint8_t *src, const int src_stride, const uint8_t *ref,
                                  const int ref_stride, const int h, __m128i *const sse,
                                  __m128i *const sum) {
    assert(h <= 256); // May overflow for larger height.
    *sum = _mm_setzero_si128();

    for (int i = 0; i < h; i += 2) {
        const __m128i s = load4x2_sse2(src, src_stride);
        const __m128i r = load4x2_sse2(ref, ref_stride);

        variance_kernel_sse2(s, r, sse, sum);
        src += 2 * src_stride;
        ref += 2 * ref_stride;
    }
}

static INLINE void variance8_sse2(const uint8_t *src, const int src_stride, const uint8_t *ref,
                                  const int ref_stride, const int h, __m128i *const sse,
                                  __m128i *const sum) {
    assert(h <= 128); // May overflow for larger height.
    *sum = _mm_setzero_si128();
    for (int i = 0; i < h; i++) {
        const __m128i s = load8_8to16_sse2(src);
        const __m128i r = load8_8to16_sse2(ref);

        variance_kernel_sse2(s, r, sse, sum);
        src += src_stride;
        ref += ref_stride;
    }
}

#define AOM_VAR_NO_LOOP_SSE2(bw, bh, bits, max_pixels)                           \
    unsigned int eb_aom_variance##bw##x##bh##_sse2(const uint8_t *src,           \
                                                   int            src_stride,    \
                                                   const uint8_t *ref,           \
                                                   int            ref_stride,    \
                                                   unsigned int * sse) {          \
        __m128i vsse = _mm_setzero_si128();                                      \
        __m128i vsum;                                                            \
        int     sum = 0;                                                         \
        variance##bw##_sse2(src, src_stride, ref, ref_stride, bh, &vsse, &vsum); \
        variance_final_##max_pixels##_pel_sse2(vsse, vsum, sse, &sum);           \
        assert(sum <= 255 * bw * bh);                                            \
        assert(sum >= -255 * bw * bh);                                           \
        return *sse - (uint32_t)(((int64_t)sum * sum) >> bits);                  \
    }

AOM_VAR_NO_LOOP_SSE2(4, 4, 4, 128);
AOM_VAR_NO_LOOP_SSE2(4, 8, 5, 128);
AOM_VAR_NO_LOOP_SSE2(4, 16, 6, 128);

AOM_VAR_NO_LOOP_SSE2(8, 4, 5, 128);
AOM_VAR_NO_LOOP_SSE2(8, 8, 6, 128);
AOM_VAR_NO_LOOP_SSE2(8, 16, 7, 128);
AOM_VAR_NO_LOOP_SSE2(8, 32, 8, 256);

static INLINE const int16_t *av1_get_interp_filter_subpel_kernel(
    const InterpFilterParams filter_params, const int32_t subpel) {
    return filter_params.filter_ptr + filter_params.taps * subpel;
}
DECLARE_ALIGNED(256, static const InterpKernel, av1_bilinear_filters[SUBPEL_SHIFTS]) = {
    {0, 0, 0, 128, 0, 0, 0, 0},
    {0, 0, 0, 120, 8, 0, 0, 0},
    {0, 0, 0, 112, 16, 0, 0, 0},
    {0, 0, 0, 104, 24, 0, 0, 0},
    {0, 0, 0, 96, 32, 0, 0, 0},
    {0, 0, 0, 88, 40, 0, 0, 0},
    {0, 0, 0, 80, 48, 0, 0, 0},
    {0, 0, 0, 72, 56, 0, 0, 0},
    {0, 0, 0, 64, 64, 0, 0, 0},
    {0, 0, 0, 56, 72, 0, 0, 0},
    {0, 0, 0, 48, 80, 0, 0, 0},
    {0, 0, 0, 40, 88, 0, 0, 0},
    {0, 0, 0, 32, 96, 0, 0, 0},
    {0, 0, 0, 24, 104, 0, 0, 0},
    {0, 0, 0, 16, 112, 0, 0, 0},
    {0, 0, 0, 8, 120, 0, 0, 0}};

DECLARE_ALIGNED(256, static const InterpKernel, av1_sub_pel_filters_4[SUBPEL_SHIFTS]) = {
    {0, 0, 0, 128, 0, 0, 0, 0},
    {0, 0, -4, 126, 8, -2, 0, 0},
    {0, 0, -8, 122, 18, -4, 0, 0},
    {0, 0, -10, 116, 28, -6, 0, 0},
    {0, 0, -12, 110, 38, -8, 0, 0},
    {0, 0, -12, 102, 48, -10, 0, 0},
    {0, 0, -14, 94, 58, -10, 0, 0},
    {0, 0, -12, 84, 66, -10, 0, 0},
    {0, 0, -12, 76, 76, -12, 0, 0},
    {0, 0, -10, 66, 84, -12, 0, 0},
    {0, 0, -10, 58, 94, -14, 0, 0},
    {0, 0, -10, 48, 102, -12, 0, 0},
    {0, 0, -8, 38, 110, -12, 0, 0},
    {0, 0, -6, 28, 116, -10, 0, 0},
    {0, 0, -4, 18, 122, -8, 0, 0},
    {0, 0, -2, 8, 126, -4, 0, 0}};
DECLARE_ALIGNED(256, static const InterpKernel, av1_sub_pel_filters_4smooth[SUBPEL_SHIFTS]) = {
    {0, 0, 0, 128, 0, 0, 0, 0},
    {0, 0, 30, 62, 34, 2, 0, 0},
    {0, 0, 26, 62, 36, 4, 0, 0},
    {0, 0, 22, 62, 40, 4, 0, 0},
    {0, 0, 20, 60, 42, 6, 0, 0},
    {0, 0, 18, 58, 44, 8, 0, 0},
    {0, 0, 16, 56, 46, 10, 0, 0},
    {0, 0, 14, 54, 48, 12, 0, 0},
    {0, 0, 12, 52, 52, 12, 0, 0},
    {0, 0, 12, 48, 54, 14, 0, 0},
    {0, 0, 10, 46, 56, 16, 0, 0},
    {0, 0, 8, 44, 58, 18, 0, 0},
    {0, 0, 6, 42, 60, 20, 0, 0},
    {0, 0, 4, 40, 62, 22, 0, 0},
    {0, 0, 4, 36, 62, 26, 0, 0},
    {0, 0, 2, 34, 62, 30, 0, 0}};
DECLARE_ALIGNED(256, static const InterpKernel, av1_sub_pel_filters_8[SUBPEL_SHIFTS]) = {
    {0, 0, 0, 128, 0, 0, 0, 0},
    {0, 2, -6, 126, 8, -2, 0, 0},
    {0, 2, -10, 122, 18, -4, 0, 0},
    {0, 2, -12, 116, 28, -8, 2, 0},
    {0, 2, -14, 110, 38, -10, 2, 0},
    {0, 2, -14, 102, 48, -12, 2, 0},
    {0, 2, -16, 94, 58, -12, 2, 0},
    {0, 2, -14, 84, 66, -12, 2, 0},
    {0, 2, -14, 76, 76, -14, 2, 0},
    {0, 2, -12, 66, 84, -14, 2, 0},
    {0, 2, -12, 58, 94, -16, 2, 0},
    {0, 2, -12, 48, 102, -14, 2, 0},
    {0, 2, -10, 38, 110, -14, 2, 0},
    {0, 2, -8, 28, 116, -12, 2, 0},
    {0, 0, -4, 18, 122, -10, 2, 0},
    {0, 0, -2, 8, 126, -6, 2, 0}};

DECLARE_ALIGNED(256, static const InterpKernel, av1_sub_pel_filters_8sharp[SUBPEL_SHIFTS]) = {
    {0, 0, 0, 128, 0, 0, 0, 0},
    {-2, 2, -6, 126, 8, -2, 2, 0},
    {-2, 6, -12, 124, 16, -6, 4, -2},
    {-2, 8, -18, 120, 26, -10, 6, -2},
    {-4, 10, -22, 116, 38, -14, 6, -2},
    {-4, 10, -22, 108, 48, -18, 8, -2},
    {-4, 10, -24, 100, 60, -20, 8, -2},
    {-4, 10, -24, 90, 70, -22, 10, -2},
    {-4, 12, -24, 80, 80, -24, 12, -4},
    {-2, 10, -22, 70, 90, -24, 10, -4},
    {-2, 8, -20, 60, 100, -24, 10, -4},
    {-2, 8, -18, 48, 108, -22, 10, -4},
    {-2, 6, -14, 38, 116, -22, 10, -4},
    {-2, 6, -10, 26, 120, -18, 8, -2},
    {-2, 4, -6, 16, 124, -12, 6, -2},
    {0, 2, -2, 8, 126, -6, 2, -2}};

DECLARE_ALIGNED(256, static const InterpKernel, av1_sub_pel_filters_8smooth[SUBPEL_SHIFTS]) = {
    {0, 0, 0, 128, 0, 0, 0, 0},
    {0, 2, 28, 62, 34, 2, 0, 0},
    {0, 0, 26, 62, 36, 4, 0, 0},
    {0, 0, 22, 62, 40, 4, 0, 0},
    {0, 0, 20, 60, 42, 6, 0, 0},
    {0, 0, 18, 58, 44, 8, 0, 0},
    {0, 0, 16, 56, 46, 10, 0, 0},
    {0, -2, 16, 54, 48, 12, 0, 0},
    {0, -2, 14, 52, 52, 14, -2, 0},
    {0, 0, 12, 48, 54, 16, -2, 0},
    {0, 0, 10, 46, 56, 16, 0, 0},
    {0, 0, 8, 44, 58, 18, 0, 0},
    {0, 0, 6, 42, 60, 20, 0, 0},
    {0, 0, 4, 40, 62, 22, 0, 0},
    {0, 0, 4, 36, 62, 26, 0, 0},
    {0, 0, 2, 34, 62, 28, 2, 0}};
// For w<=4, MULTITAP_SHARP is the same as EIGHTTAP_REGULAR
static const InterpFilterParams av1_interp_4tap[SWITCHABLE_FILTERS + 1] = {
    {(const int16_t *)av1_sub_pel_filters_4, SUBPEL_TAPS, SUBPEL_SHIFTS, EIGHTTAP_REGULAR},
    {(const int16_t *)av1_sub_pel_filters_4smooth, SUBPEL_TAPS, SUBPEL_SHIFTS, EIGHTTAP_SMOOTH},
    {(const int16_t *)av1_sub_pel_filters_4, SUBPEL_TAPS, SUBPEL_SHIFTS, EIGHTTAP_REGULAR},
    {(const int16_t *)av1_bilinear_filters, SUBPEL_TAPS, SUBPEL_SHIFTS, BILINEAR},
};
static const InterpFilterParams av1_interp_filter_params_list[SWITCHABLE_FILTERS + 1] = {
    {(const int16_t *)av1_sub_pel_filters_8, SUBPEL_TAPS, SUBPEL_SHIFTS, EIGHTTAP_REGULAR},
    {(const int16_t *)av1_sub_pel_filters_8smooth, SUBPEL_TAPS, SUBPEL_SHIFTS, EIGHTTAP_SMOOTH},
    {(const int16_t *)av1_sub_pel_filters_8sharp, SUBPEL_TAPS, SUBPEL_SHIFTS, MULTITAP_SHARP},
    {(const int16_t *)av1_bilinear_filters, SUBPEL_TAPS, SUBPEL_SHIFTS, BILINEAR}};
static INLINE const InterpFilterParams *get_4tap_interp_filter_params(
    const InterpFilter interp_filter) {
    return &av1_interp_4tap[interp_filter];
}
static INLINE const InterpFilterParams *av1_get_filter(int subpel_search) {
    assert(subpel_search >= USE_2_TAPS);

    switch (subpel_search) {
    case USE_2_TAPS: return get_4tap_interp_filter_params(BILINEAR);
    case USE_4_TAPS: return get_4tap_interp_filter_params(EIGHTTAP_REGULAR);
    case USE_8_TAPS: return &av1_interp_filter_params_list[EIGHTTAP_REGULAR];
    default: assert(0); return NULL;
    }
}

void svt_aom_upsampled_pred_sse2(MacroBlockD *xd, const struct AV1Common *const cm, int mi_row,
                                 int mi_col, const MV *const mv, uint8_t *comp_pred, int width,
                                 int height, int subpel_x_q3, int subpel_y_q3, const uint8_t *ref,
                                 int ref_stride, int subpel_search) {
    (void)xd;
    (void)cm;
    (void)mi_row;
    (void)mi_col;
    (void)mv;
    const InterpFilterParams *filter = av1_get_filter(subpel_search);
    assert(filter!=NULL);
    int filter_taps = (subpel_search <= USE_4_TAPS) ? 4 : SUBPEL_TAPS;

    if (!subpel_x_q3 && !subpel_y_q3) {
        if (width >= 16) {
            int i;
            assert(!(width & 15));
            /*Read 16 pixels one row at a time.*/
            for (i = 0; i < height; i++) {
                int j;
                for (j = 0; j < width; j += 16) {
                    xx_storeu_128(comp_pred, xx_loadu_128(ref));
                    comp_pred += 16;
                    ref += 16;
                }
                ref += ref_stride - width;
            }
        } else if (width >= 8) {
            int i;
            assert(!(width & 7));
            assert(!(height & 1));
            /*Read 8 pixels two rows at a time.*/
            for (i = 0; i < height; i += 2) {
                __m128i s0 = xx_loadl_64(ref + 0 * ref_stride);
                __m128i s1 = xx_loadl_64(ref + 1 * ref_stride);
                xx_storeu_128(comp_pred, _mm_unpacklo_epi64(s0, s1));
                comp_pred += 16;
                ref += 2 * ref_stride;
            }
        } else {
            int i;
            assert(!(width & 3));
            assert(!(height & 3));
            /*Read 4 pixels four rows at a time.*/
            for (i = 0; i < height; i+=4) {
                const __m128i row0 = xx_loadl_64(ref + 0 * ref_stride);
                const __m128i row1 = xx_loadl_64(ref + 1 * ref_stride);
                const __m128i row2 = xx_loadl_64(ref + 2 * ref_stride);
                const __m128i row3 = xx_loadl_64(ref + 3 * ref_stride);
                const __m128i reg  = _mm_unpacklo_epi64(_mm_unpacklo_epi32(row0, row1),
                                                       _mm_unpacklo_epi32(row2, row3));
                xx_storeu_128(comp_pred, reg);
                comp_pred += 16;
                ref += 4 * ref_stride;
            }
        }
    } else if (!subpel_y_q3) {
        const int16_t *const kernel =
            av1_get_interp_filter_subpel_kernel(*filter, subpel_x_q3 << 1);
        eb_aom_convolve8_horiz(ref, ref_stride, comp_pred, width, kernel, 16, NULL, -1, width, height);
    } else if (!subpel_x_q3) {
        const int16_t *const kernel =
            av1_get_interp_filter_subpel_kernel(*filter, subpel_y_q3 << 1);
        eb_aom_convolve8_vert(ref, ref_stride, comp_pred, width, NULL, -1, kernel, 16, width, height);
    } else {
        DECLARE_ALIGNED(16, uint8_t, temp[((MAX_SB_SIZE * 2 + 16) + 16) * MAX_SB_SIZE]);
        const int16_t *const kernel_x =
            av1_get_interp_filter_subpel_kernel(*filter, subpel_x_q3 << 1);
        const int16_t *const kernel_y =
            av1_get_interp_filter_subpel_kernel(*filter, subpel_y_q3 << 1);
        const uint8_t *ref_start = ref - ref_stride * ((filter_taps >> 1) - 1);
        uint8_t *      temp_start_horiz =
            (subpel_search <= USE_4_TAPS) ? temp + (filter_taps >> 1) * MAX_SB_SIZE : temp;
        uint8_t *temp_start_vert     = temp + MAX_SB_SIZE * ((filter->taps >> 1) - 1);
        int      intermediate_height = (((height - 1) * 8 + subpel_y_q3) >> 3) + filter_taps;
        assert(intermediate_height <= (MAX_SB_SIZE * 2 + 16) + 16);
        eb_aom_convolve8_horiz(ref_start,
                               ref_stride,
                               temp_start_horiz,
                               MAX_SB_SIZE,
                               kernel_x,
                               16,
                               NULL,
                               -1,
                               width,
                               intermediate_height);
        eb_aom_convolve8_vert(
            temp_start_vert, MAX_SB_SIZE, comp_pred, width, NULL, -1, kernel_y, 16, width, height);
    }
}
