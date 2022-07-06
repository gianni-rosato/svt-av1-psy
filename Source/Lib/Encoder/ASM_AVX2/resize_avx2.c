/*
 * Copyright (c) 2017, Alliance for Open Media. All rights reserved
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
#include "EbMemory_AVX2.h"
#include "common_dsp_rtcd.h"
#include "EbInterPrediction.h"
#include "EbResize.h"

static INLINE __m256i RightShiftWithRounding_S32(const __m256i v_val_d, int bits) {
    const __m256i v_bias_d = _mm256_set1_epi32((1 << bits) >> 1);
    const __m256i v_tmp_d  = _mm256_add_epi32(v_val_d, v_bias_d);
    return _mm256_srai_epi32(v_tmp_d, bits);
}
static INLINE void mm256_clamp_epi32(__m256i *x, __m256i min, __m256i max) {
    *x = _mm256_min_epi32(*x, max);
    *x = _mm256_max_epi32(*x, min);
}
static INLINE void mm256_clamp_epi16(__m256i *x, __m256i min, __m256i max) {
    *x = _mm256_min_epi16(*x, max);
    *x = _mm256_max_epi16(*x, min);
}
static INLINE void mm_blend_load_lo(const uint8_t *const input, int blend_len, __m128i *dst) {
    const int      bits_to_mask = blend_len * 8;
    const uint64_t mask_c       = (1LL << bits_to_mask) - 1;
    const __m128i  blend_mask   = _mm_set1_epi64x(mask_c);
    const __m128i  src_0        = _mm_set1_epi8(input[0]);
    *dst = _mm_loadl_epi64((const __m128i *)(input)); // lo 64-bit: 0 1 2 3 4 5 6 7 8
    *dst = _mm_slli_epi64(*dst, bits_to_mask); // lo 64-bit: x x x 0 1 2 3 4 5
    *dst = _mm_blendv_epi8(*dst, src_0, blend_mask);
    *dst = _mm_cvtepu8_epi16(*dst);
    return;
}
static INLINE void mm_blend_load_hi(const uint8_t *const input, int blend_len, __m128i *dst) {
    const int      bits_to_mask = blend_len * 8;
    const uint64_t mask_c       = ~((1LL << (64 - bits_to_mask)) - 1);
    const __m128i  blend_mask   = _mm_set1_epi64x(mask_c);
    const __m128i  src_last     = _mm_set1_epi8(input[8 - blend_len - 1]);
    *dst = _mm_loadl_epi64((const __m128i *)(input - blend_len)); // lo 64-bit: x x x x 0 1 2 3
    *dst = _mm_srli_epi64(*dst, bits_to_mask); // lo 64-bit: 0 1 2 3 x x x x
    *dst = _mm_blendv_epi8(*dst, src_last, blend_mask);
    *dst = _mm_cvtepu8_epi16(*dst);
    return;
}
static INLINE void interpolate_core_w16_avx2(__m128i short_src[16], __m128i short_filter[16],
                                             uint8_t *optr) {
    const __m256i min = _mm256_set1_epi16(0);
    const __m256i max = _mm256_set1_epi16(255);

    __m256i sum[8];

    for (int j = 0; j < 4; ++j) {
        __m256i src_ab    = _mm256_set_m128i(short_src[j + 4], short_src[j]);
        __m256i filter_ab = _mm256_set_m128i(short_filter[j + 4], short_filter[j]);

        // two pixel per vector
        sum[j] = _mm256_madd_epi16(src_ab, filter_ab);
        // sum[0]: P00ab P00cd P00ef P00gh P04ab P04cd P04ef P04gh
        // sum[1]: P01ab P01cd P01ef P01gh P05ab P05cd P05ef P05gh
        // sum[2]: P02ab P02cd P02ef P02gh P06ab P06cd P06ef P06gh
        // sum[3]: P03ab P03cd P03ef P03gh P07ab P07cd P07ef P07gh
    }
    for (int j = 4; j < 8; ++j) {
        __m256i src_ab    = _mm256_set_m128i(short_src[j + 8], short_src[j + 4]);
        __m256i filter_ab = _mm256_set_m128i(short_filter[j + 8], short_filter[j + 4]);

        // two pixel per vector
        sum[j] = _mm256_madd_epi16(src_ab, filter_ab);
        // sum[0]: P00ab P00cd P00ef P00gh P04ab P04cd P04ef P04gh
        // sum[1]: P01ab P01cd P01ef P01gh P05ab P05cd P05ef P05gh
        // sum[2]: P02ab P02cd P02ef P02gh P06ab P06cd P06ef P06gh
        // sum[3]: P03ab P03cd P03ef P03gh P07ab P07cd P07ef P07gh
    }

    __m256i sum0to7, sum8to15;
    {
        // P00abcd P00efgh P01abcd P01efgh P04abcd P04efgh P05abcd P05efgh
        __m256i sum01 = _mm256_hadd_epi32(sum[0], sum[1]);
        // P02abcd P02efgh P03abcd P03efgh P06abcd P06efgh P07abcd P07efgh
        __m256i sum23 = _mm256_hadd_epi32(sum[2], sum[3]);
        // P00 P01 P02 P03 P04 P05 P06 P07
        __m256i sum0123 = _mm256_hadd_epi32(sum01, sum23);
        sum0to7         = RightShiftWithRounding_S32(sum0123, FILTER_BITS);
    }
    {
        __m256i sum01   = _mm256_hadd_epi32(sum[4], sum[5]);
        __m256i sum23   = _mm256_hadd_epi32(sum[6], sum[7]);
        __m256i sum0123 = _mm256_hadd_epi32(sum01, sum23);
        sum8to15        = RightShiftWithRounding_S32(sum0123, FILTER_BITS);
    }

    // 16* 16-bit: P00 P01 P02 P03 P08 P09 P10 P11 P04 P05 P06 P07 P12 P13 P14 P15
    __m256i sum0to15 = _mm256_packs_epi32(sum0to7, sum8to15);
    sum0to15         = _mm256_permute4x64_epi64(sum0to15, _MM_SHUFFLE(3, 1, 2, 0));
    // P00 000 P01 000 P02 000 P03 000 P04 000 P05 000 P06 000 P07 000
    // P08 000 P09 000 P10 000 P11 000 P12 000 P13 000 P14 000 P15 000
    mm256_clamp_epi16(&sum0to15, min, max);

    __m128i lo_lane    = _mm256_castsi256_si128(sum0to15);
    __m128i hi_lane    = _mm256_extracti128_si256(sum0to15, 1);
    __m128i m128_0to15 = _mm_packus_epi16(lo_lane, hi_lane);

    _mm_storeu_si128((__m128i *)optr, m128_0to15);
}
static INLINE void highbd_interpolate_core_w8_avx2(__m256i src[8], __m256i filter[8],
                                                   uint16_t *optr, __m256i max) {
    const int     steps = 8;
    const __m256i min = _mm256_set1_epi32(0);

    __m256i sum[8];
    for (int i = 0; i < steps; ++i) {
        sum[i] = _mm256_mullo_epi32(src[i], filter[i]);
    }

    __m256i sum0to3, sum4to7;
    {
        __m256i sum01   = _mm256_hadd_epi32(sum[0], sum[1]);
        sum01           = _mm256_permute4x64_epi64(sum01, _MM_SHUFFLE(3, 1, 2, 0));
        __m256i sum23   = _mm256_hadd_epi32(sum[2], sum[3]);
        sum23           = _mm256_permute4x64_epi64(sum23, _MM_SHUFFLE(3, 1, 2, 0));
        sum0to3         = _mm256_hadd_epi32(sum01, sum23);
        sum0to3         = _mm256_permute4x64_epi64(sum0to3, _MM_SHUFFLE(3, 1, 2, 0));
    }
    {
        __m256i sum45 = _mm256_hadd_epi32(sum[4], sum[5]);
        sum45         = _mm256_permute4x64_epi64(sum45, _MM_SHUFFLE(3, 1, 2, 0));
        __m256i sum67 = _mm256_hadd_epi32(sum[6], sum[7]);
        sum67         = _mm256_permute4x64_epi64(sum67, _MM_SHUFFLE(3, 1, 2, 0));
        sum4to7       = _mm256_hadd_epi32(sum45, sum67);
        sum4to7       = _mm256_permute4x64_epi64(sum4to7, _MM_SHUFFLE(3, 1, 2, 0));
    }

    __m256i sum8 = _mm256_hadd_epi32(sum0to3, sum4to7);
    sum8         = _mm256_permute4x64_epi64(sum8, _MM_SHUFFLE(3, 1, 2, 0));
    sum8         = RightShiftWithRounding_S32(sum8, FILTER_BITS);
    mm256_clamp_epi32(&sum8, min, max);

    __m128i lo_lane   = _mm256_castsi256_si128(sum8);
    __m128i hi_lane   = _mm256_extracti128_si256(sum8, 1);
    __m128i m128_sum8 = _mm_packus_epi32(lo_lane, hi_lane); // 8x 16-bit

    _mm_storeu_si128((__m128i *)optr, m128_sum8);
}
static INLINE void highbd_interpolate_core_w8_mid_part_avx2(const uint16_t *const input,
                                                            uint16_t            **output,
                                                            const int16_t *interp_filters, int *py,
                                                            const int delta, int length,
                                                            __m256i max) {
    const int steps = 8;
    uint16_t *optr = *output;
    int      y    = *py;

    __m256i filter[8];
    __m256i src[8];
    for (int i = 0; i < length; i += steps) {
        for (int j = 0; j < steps; ++j) {
            int            int_pel  = y >> RS_SCALE_SUBPEL_BITS;
            int            sub_pel  = (y >> RS_SCALE_EXTRA_BITS) & RS_SUBPEL_MASK; // 0~63
            const int16_t *filter_c = &interp_filters[sub_pel * SUBPEL_TAPS];
            const uint16_t *in_ptr   = &input[int_pel - SUBPEL_TAPS / 2 + 1];
            __m128i        src_128  = _mm_lddqu_si128((const __m128i *)in_ptr);
            src[j]                  = _mm256_cvtepu16_epi32(src_128);
            __m128i filter_128      = _mm_lddqu_si128((const __m128i *)filter_c);
            filter[j]               = _mm256_cvtepi16_epi32(filter_128);
            y += delta;
        }

        highbd_interpolate_core_w8_avx2(src, filter, optr, max);
        optr += steps;
    }
    *output = optr;
    *py     = y;
}
static INLINE void mm256_blend_load_lo(const uint16_t *const input, int blend_len, __m256i *dst) {
    uint16_t tmp[8];
    for (int i = 0; i < blend_len; ++i) {
        tmp[i] = input[0];
    }
    memcpy(tmp + blend_len, input, (8 - blend_len) * sizeof(uint16_t));

    __m128i src_128 = _mm_lddqu_si128((__m128i *)tmp);
    *dst            = _mm256_cvtepu16_epi32(src_128);
    return;
}
static INLINE void mm256_blend_load_hi(const uint16_t *const input, int blend_len, __m256i *dst) {
    uint16_t tmp[8];
    memcpy(tmp, input, (8 - blend_len) * sizeof(uint16_t));
    for (int i = 0; i < blend_len; ++i) {
        tmp[8 - blend_len + i] = input[8 - blend_len-1];
    }

    __m128i src_128 = _mm_lddqu_si128((__m128i *)tmp);
    *dst            = _mm256_cvtepu16_epi32(src_128);
    return;
}
static INLINE void highbd_interpolate_core_w8_init_part_avx2(const uint16_t *const input,
                                                             uint16_t            **output,
                                                             const int16_t *interp_filters, int *py,
                                                             const int delta, __m256i max) {
    const int steps = 8;
    uint16_t  *optr  = *output;
    int       y     = *py;

    __m256i filter[8];
    __m256i src[8];

    for (int j = 0; j < steps; ++j) {
        int            int_pel  = y >> RS_SCALE_SUBPEL_BITS;
        int            sub_pel  = (y >> RS_SCALE_EXTRA_BITS) & RS_SUBPEL_MASK; // 0~63
        const int16_t *filter_c = &interp_filters[sub_pel * SUBPEL_TAPS];
        int            in_offset = int_pel - SUBPEL_TAPS / 2 + 1;
        if (in_offset < 0) {
            mm256_blend_load_lo(input, SUBPEL_TAPS / 2 - int_pel - 1, &src[j]);
        } else {
            const uint16_t *in_ptr  = &input[in_offset];
            __m128i src_128 = _mm_lddqu_si128((const __m128i *)in_ptr);
            src[j]          = _mm256_cvtepu16_epi32(src_128);
        }
        __m128i filter_128      = _mm_lddqu_si128((const __m128i *)filter_c);
        filter[j]               = _mm256_cvtepi16_epi32(filter_128);
        y += delta;
    }

    highbd_interpolate_core_w8_avx2(src, filter, optr, max);
    optr += steps;

    *output = optr;
    *py     = y;
}
static INLINE void highbd_interpolate_core_w8_end_part_avx2(const uint16_t *const input,
                                                            int in_length, uint16_t **output,
                                                            const int16_t *interp_filters, int *py,
                                                            const int delta, __m256i max) {
    const int steps = 8;
    uint16_t  *optr  = *output;
    int       y     = *py;

    __m256i filter[8];
    __m256i src[8];

    for (int j = 0; j < steps; ++j) {
        int            int_pel   = y >> RS_SCALE_SUBPEL_BITS;
        int            sub_pel   = (y >> RS_SCALE_EXTRA_BITS) & RS_SUBPEL_MASK; // 0~63
        const int16_t *filter_c  = &interp_filters[sub_pel * SUBPEL_TAPS];
        int            in_offset = int_pel - SUBPEL_TAPS / 2 + 1;
        if (int_pel - SUBPEL_TAPS / 2 + 1 + 8 > in_length) {
            int blend_length = int_pel - SUBPEL_TAPS / 2 + 1 + 8 - in_length;
            mm256_blend_load_hi(&input[int_pel - SUBPEL_TAPS / 2 + 1], blend_length, &src[j]);
        } else {
            const uint16_t *in_ptr  = &input[in_offset];
            __m128i        src_128 = _mm_lddqu_si128((const __m128i *)in_ptr);
            src[j]                 = _mm256_cvtepu16_epi32(src_128);
        }
        __m128i filter_128 = _mm_lddqu_si128((const __m128i *)filter_c);
        filter[j]          = _mm256_cvtepi16_epi32(filter_128);
        y += delta;
    }

    highbd_interpolate_core_w8_avx2(src, filter, optr, max);
    optr += steps;

    *output = optr;
    *py     = y;
}
static INLINE void interpolate_core_w16_mid_part_avx2(const uint8_t *const input, uint8_t **output,
                                             const int16_t *interp_filters, int *py,
                                             const int delta, int length) {
    uint8_t *optr = *output;
    int      y    = *py;

    __m128i  short_filter[16];
    __m128i  short_src[16];
    for (int i = 0; i < length; i += 16) {
        for (int j = 0; j < 8; ++j) {
            int            int_pel = y >> RS_SCALE_SUBPEL_BITS;
            int            sub_pel = (y >> RS_SCALE_EXTRA_BITS) & RS_SUBPEL_MASK; // 0~63
            const int16_t *filter  = &interp_filters[sub_pel * SUBPEL_TAPS];
            const uint8_t *in_ptr  = &input[int_pel - SUBPEL_TAPS / 2 + 1];
            short_src[2 * j]       = _mm_loadl_epi64((const __m128i *)in_ptr);
            short_src[2 * j]       = _mm_cvtepu8_epi16(short_src[2 * j]);
            short_filter[2 * j]    = _mm_loadu_si128((const __m128i *)filter);
            y += delta;

            int_pel                 = y >> RS_SCALE_SUBPEL_BITS;
            sub_pel                 = (y >> RS_SCALE_EXTRA_BITS) & RS_SUBPEL_MASK; // 0~63
            filter                  = &interp_filters[sub_pel * SUBPEL_TAPS];
            in_ptr                  = &input[int_pel - SUBPEL_TAPS / 2 + 1];
            short_src[2 * j + 1]    = _mm_loadl_epi64((const __m128i *)in_ptr);
            short_src[2 * j + 1]    = _mm_cvtepu8_epi16(short_src[2 * j + 1]);
            short_filter[2 * j + 1] = _mm_loadu_si128((const __m128i *)filter);
            y += delta;
        }

        interpolate_core_w16_avx2(short_src, short_filter, optr);
        optr += 16;
    }
    *output = optr;
    *py     = y;
}
static INLINE void interpolate_core_w16_init_part_avx2(const uint8_t *const input, uint8_t **output,
                                                       const int16_t *interp_filters, int *py,
                                                       const int delta) {
    uint8_t *optr = *output;
    int      y    = *py;

    __m128i short_filter[16];
    __m128i short_src[16];

    for (int j = 0; j < 8; ++j) {
        int            int_pel = y >> RS_SCALE_SUBPEL_BITS;
        int            sub_pel = (y >> RS_SCALE_EXTRA_BITS) & RS_SUBPEL_MASK; // 0~63
        const int16_t *filter  = &interp_filters[sub_pel * SUBPEL_TAPS];
        if (int_pel - SUBPEL_TAPS / 2 + 1 < 0) {
            mm_blend_load_lo(input, SUBPEL_TAPS / 2 - int_pel - 1, &short_src[2 * j]);
        } else {
            const uint8_t *in_ptr = &input[int_pel - SUBPEL_TAPS / 2 + 1];
            short_src[2 * j]      = _mm_loadl_epi64((const __m128i *)in_ptr);
            short_src[2 * j]      = _mm_cvtepu8_epi16(short_src[2 * j]);
        }
        short_filter[2 * j] = _mm_loadu_si128((const __m128i *)filter);
        y += delta;

        int_pel = y >> RS_SCALE_SUBPEL_BITS;
        sub_pel = (y >> RS_SCALE_EXTRA_BITS) & RS_SUBPEL_MASK; // 0~63
        filter  = &interp_filters[sub_pel * SUBPEL_TAPS];
        if (int_pel - SUBPEL_TAPS / 2 + 1 < 0) {
            mm_blend_load_lo(input, SUBPEL_TAPS / 2 - int_pel - 1, &short_src[2 * j + 1]);
        } else {
            const uint8_t *in_ptr = &input[int_pel - SUBPEL_TAPS / 2 + 1];
            short_src[2 * j + 1]  = _mm_loadl_epi64((const __m128i *)in_ptr);
            short_src[2 * j + 1]  = _mm_cvtepu8_epi16(short_src[2 * j + 1]);
        }
        short_filter[2 * j + 1] = _mm_loadu_si128((const __m128i *)filter);
        y += delta;
    }

    interpolate_core_w16_avx2(short_src, short_filter, optr);
    optr += 16;

    *output = optr;
    *py     = y;
}
static INLINE void interpolate_core_w16_end_part_avx2(const uint8_t *const input, int in_length,
                                                      uint8_t      **output,
                                                      const int16_t *interp_filters, int *py,
                                                      const int delta) {
    uint8_t *optr = *output;
    int      y    = *py;

    __m128i short_filter[16];
    __m128i short_src[16];

    for (int j = 0; j < 8; ++j) {
        int            int_pel = y >> RS_SCALE_SUBPEL_BITS;
        int            sub_pel = (y >> RS_SCALE_EXTRA_BITS) & RS_SUBPEL_MASK; // 0~63
        const int16_t *filter  = &interp_filters[sub_pel * SUBPEL_TAPS];
        if (int_pel - SUBPEL_TAPS / 2 + 1 + 8 > in_length) {
            int blend_length = int_pel - SUBPEL_TAPS / 2 + 1 + 8 - in_length;
            mm_blend_load_hi(
                &input[int_pel - SUBPEL_TAPS / 2 + 1], blend_length, &short_src[2 * j]);
        } else {
            const uint8_t *in_ptr = &input[int_pel - SUBPEL_TAPS / 2 + 1];
            short_src[2 * j]      = _mm_loadl_epi64((const __m128i *)in_ptr);
            short_src[2 * j]      = _mm_cvtepu8_epi16(short_src[2 * j]);
        }
        short_filter[2 * j] = _mm_loadu_si128((const __m128i *)filter);
        y += delta;

        int_pel = y >> RS_SCALE_SUBPEL_BITS;
        sub_pel = (y >> RS_SCALE_EXTRA_BITS) & RS_SUBPEL_MASK; // 0~63
        filter  = &interp_filters[sub_pel * SUBPEL_TAPS];
        if (int_pel - SUBPEL_TAPS / 2 + 1 + 8 > in_length) {
            int blend_length = int_pel - SUBPEL_TAPS / 2 + 1 + 8 - in_length;
            mm_blend_load_hi(
                &input[int_pel - SUBPEL_TAPS / 2 + 1], blend_length, &short_src[2 * j + 1]);
        } else {
            const uint8_t *in_ptr = &input[int_pel - SUBPEL_TAPS / 2 + 1];
            short_src[2 * j + 1]  = _mm_loadl_epi64((const __m128i *)in_ptr);
            short_src[2 * j + 1]  = _mm_cvtepu8_epi16(short_src[2 * j + 1]);
        }
        short_filter[2 * j + 1] = _mm_loadu_si128((const __m128i *)filter);
        y += delta;
    }

    interpolate_core_w16_avx2(short_src, short_filter, optr);
    optr += 16;

    *output = optr;
    *py     = y;
}

static INLINE __m256i mm256_load2_m128i(__m128i hi, __m128i lo) {
    __m256i c = _mm256_castsi128_si256(lo);
    return _mm256_inserti128_si256(c, hi, 1);
}
static INLINE void down2_prepare_4vec(const __m128i *src, __m256i *dst) {
    __m128i tmp = _mm_alignr_epi8(src[1], src[0], 4);
    dst[0]      = mm256_load2_m128i(tmp, src[0]);
    dst[1]      = mm256_load2_m128i(_mm_alignr_epi8(src[1], src[0], 12),
                               _mm_alignr_epi8(src[1], src[0], 8));
    dst[2]      = mm256_load2_m128i(_mm_alignr_epi8(src[2], src[1], 4), src[1]);
    dst[3]      = mm256_load2_m128i(_mm_alignr_epi8(src[2], src[1], 12),
                               _mm_alignr_epi8(src[2], src[1], 8));
}
static INLINE __m256i down2_get_8sum(__m256i *vec, const __m256i *filter_4x, const __m256i *subtrahend) {
    vec[0]       = _mm256_shufflelo_epi16(vec[0], _MM_SHUFFLE(0, 1, 2, 3));
    vec[0]       = _mm256_shufflehi_epi16(vec[0], _MM_SHUFFLE(0, 1, 2, 3));
    __m256i sum0 = _mm256_add_epi16(vec[0], vec[1]);
    // sum0: P00a P00b P00c P00d P02a P02b P02c P02d   P01a P01b P01c P01d P03a P03b P03c P03d
    sum0 = _mm256_mullo_epi16(sum0, *filter_4x);
    sum0 = _mm256_subs_epi16(sum0, *subtrahend);  // to avoid exceeding 16-bit after P00a+P00b
    // sum0: P00a P00b P00c P00d P01a P01b P01c P01d   P02a P02b P02c P02d P03a P03b P03c P03d
    sum0 = _mm256_permute4x64_epi64(sum0, _MM_SHUFFLE(3, 1, 2, 0));

    vec[2]       = _mm256_shufflelo_epi16(vec[2], _MM_SHUFFLE(0, 1, 2, 3));
    vec[2]       = _mm256_shufflehi_epi16(vec[2], _MM_SHUFFLE(0, 1, 2, 3));
    __m256i sum1 = _mm256_add_epi16(vec[2], vec[3]);
    // sum1: P04a P04b P04c P04d P06a P06b P06c P06d   P05a P05b P05c P05d P07a P07b P07c P07d
    sum1 = _mm256_mullo_epi16(sum1, *filter_4x);
    sum1 = _mm256_subs_epi16(sum1, *subtrahend); // to avoid exceeding 16-bit after P00a+P00b
    // sum1: P04a P04b P04c P04d P05a P05b P05c P05d   P06a P06b P06c P06d P07a P07b P07c P07d
    sum1 = _mm256_permute4x64_epi64(sum1, _MM_SHUFFLE(3, 1, 2, 0));

    // P00ab P00cd P01ab P01cd P04ab P04cd P05ab P05cd  P02ab P02cd P03ab P03cd P06ab P06cd P07ab P07cd
    __m256i sum01 = _mm256_hadd_epi16(sum0, sum1);
    // P00ab P00cd P01ab P01cd P02ab P02cd P03ab P03cd  P04ab P04cd P05ab P05cd P06ab P06cd P07ab P07cd
    sum01 = _mm256_permute4x64_epi64(sum01, _MM_SHUFFLE(3, 1, 2, 0));
    return sum01;
}
static INLINE void down2_symeven_w16_avx2(const __m128i *load, const __m256i *filter_4x,
                                          uint8_t *optr) {
    const __m256i min         = _mm256_set1_epi16(0);
    const __m256i max         = _mm256_set1_epi16(255);
    const __m256i base_sum    = _mm256_set1_epi16(64);
    const __m256i subtrahend  = _mm256_set1_epi16(1 << 11);
    const __m256i recover_sub = _mm256_set1_epi16(1 << (13 - FILTER_BITS));

    __m256i vec[8];
    down2_prepare_4vec(load, vec);
    down2_prepare_4vec(load + 2, vec + 4);

    // obtain: P00ab P00cd P01ab P01cd P02ab P02cd P03ab P03cd  P04ab P04cd P05ab P05cd P06ab P06cd P07ab P07cd
    __m256i sum01 = down2_get_8sum(&vec[0], filter_4x, &subtrahend);
    __m256i sum23 = down2_get_8sum(&vec[4], filter_4x, &subtrahend);

    // P00 P01 P02 P03 P08 P09 P10 P11  P04 P05 P06 P07 P12 P13 P14 P15
    __m256i sum0123 = _mm256_hadd_epi16(sum01, sum23);
    // P00 P01 P02 P03 P04 P05 P06 P07  P08 P09 P10 P11 P12 P13 P14 P15
    sum0123 = _mm256_permute4x64_epi64(sum0123, _MM_SHUFFLE(3, 1, 2, 0));
    sum0123 = _mm256_adds_epi16(sum0123, base_sum);

    sum0123 = _mm256_srai_epi16(sum0123, FILTER_BITS);
    sum0123 = _mm256_adds_epi16(sum0123, recover_sub); // recover subtracted 2^13/2^7
    mm256_clamp_epi16(&sum0123, min, max);

    const __m128i lo_lane    = _mm256_castsi256_si128(sum0123);
    const __m128i hi_lane    = _mm256_extracti128_si256(sum0123, 1);
    __m128i m128_0to15 = _mm_packus_epi16(lo_lane, hi_lane);
    _mm_storeu_si128((__m128i *)optr, m128_0to15);
}
static INLINE void down2_symeven_w16_mid_part_avx2(const uint8_t *const input, int length,
                                                   uint8_t **output, const __m256i *filter_4x) {
    assert((length % 32) == 0);
    const uint8_t *in   = input - 3;
    uint8_t *optr = *output;
    __m128i  load[5];

    int i = 0;
    load[4] = _mm_loadl_epi64((const __m128i *)in);
    load[4] = _mm_cvtepu8_epi16(load[4]);
    do {
        load[0] = load[4];
        for (int n = 1; n < 5; ++n) {
            load[n] = _mm_loadl_epi64((const __m128i *)(in + n * 8));
            load[n] = _mm_cvtepu8_epi16(load[n]);
        }

        down2_symeven_w16_avx2(load, filter_4x, optr);

        optr += 16;
        i += 32;
        in += 32;
    } while (i < length);
    *output = optr;
}
static INLINE void down2_symeven_w16_init_part_avx2(const uint8_t *const input, uint8_t **output,
                                                    const __m256i *filter_4x) {
    const uint8_t *in   = input - 3;
    uint8_t       *optr = *output;
    __m128i        load[5];
    mm_blend_load_lo(input, 3, &load[0]);
    for (int n = 1; n < 5; ++n) {
        load[n] = _mm_loadl_epi64((const __m128i *)(in + n * 8));
        load[n] = _mm_cvtepu8_epi16(load[n]);
    }

    down2_symeven_w16_avx2(load, filter_4x, optr);
    optr += 16;

    *output = optr;
}
static INLINE void down2_symeven_w16_end_part_avx2(const uint8_t *const input, uint8_t **output,
                                                   const __m256i *filter_4x, int len_to_end) {
    const uint8_t *in   = input - 3;
    const uint8_t *end  = input + len_to_end;
    uint8_t       *optr = *output;
    __m128i        load[5];

    for (int n = 0; n < 4; ++n) {
        load[n] = _mm_loadl_epi64((const __m128i *)(in + n * 8));
        load[n] = _mm_cvtepu8_epi16(load[n]);
    }

    mm_blend_load_hi(in+32, (int)(in+32+8-end), &load[4]);

    down2_symeven_w16_avx2(load, filter_4x, optr);
    optr += 16;

    *output = optr;
}
static INLINE void highbd_down2_symeven_prepare_4pt(const __m128i load[2], __m256i filter_2x,
                                                    __m256i *vec_4pt) {
    // -3 -2 -1 00 01 02 03 04
    __m256i vec256_tmp_0 = _mm256_cvtepu16_epi32(load[0]);
    // -1 00 01 02 03 04 05 06
    __m128i vec128_tmp_1 = _mm_alignr_epi8(
        load[1], load[0], 4); // cat and shift right 2x 16-bit (4 bytes)
    __m256i vec256_tmp_1 = _mm256_cvtepu16_epi32(vec128_tmp_1);
    // -3 -2 -1 00 -1 00 01 02
    __m256i vec0 = _mm256_permute2x128_si256(vec256_tmp_0, vec256_tmp_1, 0x20);
    // 00 -1 -2 -3 02 01 00 -1
    vec0 = _mm256_shuffle_epi32(vec0, _MM_SHUFFLE(0, 1, 2, 3));

    // 01 02 03 04 03 04 05 06
    __m256i vec1 = _mm256_permute2x128_si256(vec256_tmp_0, vec256_tmp_1, 0x31);

    __m256i p01 = _mm256_add_epi32(vec0, vec1);
    // P00a P00b P00c P00d P01a P01b P01c P01d
    p01 = _mm256_mullo_epi32(p01, filter_2x);

    // 04 03 02 01 06 05 04 03
    vec0 = _mm256_shuffle_epi32(vec1, _MM_SHUFFLE(0, 1, 2, 3));
    // 05 06 07 08 09 10 11 12
    vec256_tmp_0 = _mm256_cvtepu16_epi32(load[1]);
    // 07 08 09 10 11 12 -3 -2
    vec128_tmp_1 = _mm_alignr_epi8(load[0], load[1], 4);
    vec256_tmp_1 = _mm256_cvtepu16_epi32(vec128_tmp_1);
    // 05 06 07 08 07 08 09 10
    vec1         = _mm256_permute2x128_si256(vec256_tmp_0, vec256_tmp_1, 0x20);

    __m256i p23 = _mm256_add_epi32(vec0, vec1);
    // P02a P02b P02c P02d P03a P03b P03c P03d
    p23 = _mm256_mullo_epi32(p23, filter_2x);

    // P00ab P00cd P02ab P02cd P01ab P01cd P03ab P03cd
    *vec_4pt = _mm256_hadd_epi32(p01, p23);
    // P00ab P00cd P01ab P01cd P02ab P02cd P03ab P03cd
    *vec_4pt = _mm256_permute4x64_epi64(*vec_4pt, _MM_SHUFFLE(3, 1, 2, 0));
}
static INLINE void highbd_down2_symeven_w8_mid_part_avx2(const uint16_t *const input, int length,
                                                          uint16_t **output, __m256i filter_2x,
                                                          __m256i max) {
    const int     steps    = 8; // output 8 pt by each loop
    const __m256i min      = _mm256_set1_epi32(0);
    const __m256i base_sum = _mm256_set1_epi32(64);

    const uint16_t *in   = input - 3;
    uint16_t       *optr = *output;
    int             i    = 0;

    __m128i load[2];
    // -3 -2 -1 00 01 02 03 04
    load[1] = _mm_lddqu_si128((__m128i *)in);

    do {
        load[0] = load[1];
        // 05 06 07 08 09 10 11 12
        load[1] = _mm_lddqu_si128((const __m128i *)(in + 8));

        {
            __m256i vec0123, vec4567;
            highbd_down2_symeven_prepare_4pt(load, filter_2x, &vec0123);

            load[0] = load[1];
            load[1] = _mm_lddqu_si128((const __m128i *)(in + 16));
            highbd_down2_symeven_prepare_4pt(load, filter_2x, &vec4567);

            // P00ab P00cd P01ab P01cd P02ab P02cd P03ab P03cd
            // P04ab P04cd P05ab P05cd P06ab P06cd P07ab P07cd

            // P00 P01 P04 P05 P02 P03 P06 P07
            __m256i vec_8pt = _mm256_hadd_epi32(vec0123, vec4567);
            vec_8pt         = _mm256_permute4x64_epi64(vec_8pt, _MM_SHUFFLE(3, 1, 2, 0));

            vec_8pt = _mm256_add_epi32(vec_8pt, base_sum);
            vec_8pt = _mm256_srai_epi32(vec_8pt, FILTER_BITS);
            mm256_clamp_epi32(&vec_8pt, min, max);

            // 8x 32-bit => 8x 16-bit
            __m128i lo_lane   = _mm256_castsi256_si128(vec_8pt);
            __m128i hi_lane   = _mm256_extracti128_si256(vec_8pt, 1);
            __m128i m128_sum8 = _mm_packus_epi32(lo_lane, hi_lane); // 8x 16-bit

            _mm_storeu_si128((__m128i *)optr, m128_sum8);
        }

        // DEBUG: compare with c
        /*{
            const int16_t *filter = av1_down2_symeven_half_filter;
            const int      filter_len_half = 4;
            const int      bd              = 10;
            uint16_t      *input              = in + 3;
            uint16_t       output[8];
            uint16_t      *optr_c = output;

            for (int i = 0; i < 16; i += 2) {
                int sum = (1 << (FILTER_BITS - 1));
                for (int j = 0; j < filter_len_half; ++j) {
                    sum += (input[i - j] + input[AOMMIN(i + 1 + j, length - 1)]) * filter[j];
                }
                sum >>= FILTER_BITS;
                *optr_c++ = clip_pixel_highbd(sum, bd);
            }

            if (memcmp(output, optr, 8 * sizeof(uint16_t)) != 0) {
                fprintf(stderr, "mismatch\n");
            }
        }*/

        optr += steps;
        i += steps * 2;
        in += steps * 2;
    } while (i < length);

    *output = optr;
}
static INLINE void highbd_down2_symeven_w8_init_part_avx2(const uint16_t *const input,
                                                           uint16_t **output, __m256i filter_2x,
                                                           __m256i max) {
    const int       steps = 8;
    const __m256i   min      = _mm256_set1_epi32(0);
    const __m256i   base_sum = _mm256_set1_epi32(64);

    const uint16_t *in   = input - 3;
    uint16_t       *optr = *output;

    __m128i load[2];
    {
        int      blend_len = 3;
        uint16_t tmp[8];
        for (int i = 0; i < blend_len; ++i) { tmp[i] = input[0]; }
        memcpy(tmp + blend_len, input, (8 - blend_len) * sizeof(uint16_t));

        // -3 -2 -1 00 01 02 03 04
        load[0] = _mm_lddqu_si128((__m128i *)tmp);
    }
    // 05 06 07 08 09 10 11 12
    load[1] = _mm_lddqu_si128((const __m128i *)(in + 8));

    __m256i vec0123, vec4567;
    highbd_down2_symeven_prepare_4pt(load, filter_2x, &vec0123);

    load[0] = load[1];
    load[1] = _mm_lddqu_si128((const __m128i *)(in + 16));
    highbd_down2_symeven_prepare_4pt(load, filter_2x, &vec4567);

    // P00ab P00cd P01ab P01cd P02ab P02cd P03ab P03cd
    // P04ab P04cd P05ab P05cd P06ab P06cd P07ab P07cd

    // P00 P01 P04 P05 P02 P03 P06 P07
    __m256i vec_8pt = _mm256_hadd_epi32(vec0123, vec4567);
    vec_8pt         = _mm256_permute4x64_epi64(vec_8pt, _MM_SHUFFLE(3, 1, 2, 0));

    vec_8pt = _mm256_add_epi32(vec_8pt, base_sum);
    vec_8pt = _mm256_srai_epi32(vec_8pt, FILTER_BITS);
    mm256_clamp_epi32(&vec_8pt, min, max);

    // 8x 32-bit => 8x 16-bit
    __m128i lo_lane   = _mm256_castsi256_si128(vec_8pt);
    __m128i hi_lane   = _mm256_extracti128_si256(vec_8pt, 1);
    __m128i m128_sum8 = _mm_packus_epi32(lo_lane, hi_lane); // 8x 16-bit

    _mm_storeu_si128((__m128i *)optr, m128_sum8);
    optr += steps;

    *output = optr;
}
static INLINE void highbd_down2_symeven_w8_end_part_avx2(const uint16_t *const input,
                                                          uint16_t            **output,
                                                          const __m256i filter_2x, int len_to_end,
                                                          __m256i max) {
    const int     steps    = 8;
    const __m256i min      = _mm256_set1_epi32(0);
    const __m256i base_sum = _mm256_set1_epi32(64);

    const uint16_t *in   = input - 3;
    const uint16_t *end  = input + len_to_end;
    uint16_t       *optr = *output;

    __m128i load[2];

    load[0] = _mm_lddqu_si128((const __m128i *)in);
    {
        if (in + 16 > end) {
            int      blend_len = (int)(in + 16 - end);
            uint16_t tmp[8];
            memcpy(tmp, in + 8, (8 - blend_len) * sizeof(uint16_t));
            for (int i = 0; i < blend_len; ++i) { tmp[8 - blend_len + i] = *(end - 1); }

            load[1] = _mm_lddqu_si128((__m128i *)tmp);
        } else {
            load[1] = _mm_lddqu_si128((const __m128i *)(in+8));
        }
    }

    {
        __m256i vec0123, vec4567;
        highbd_down2_symeven_prepare_4pt(load, filter_2x, &vec0123);

        load[0] = load[1];
        if (in + 24 > end) {
            int      blend_len = (int)(in + 24 - end);
            uint16_t tmp[8];
            memcpy(tmp, in + 16, (8 - blend_len) * sizeof(uint16_t));
            for (int i = 0; i < blend_len; ++i) { tmp[8 - blend_len + i] = *(end - 1); }

            load[1] = _mm_lddqu_si128((__m128i *)tmp);
        } else {
            load[1] = _mm_lddqu_si128((const __m128i *)(in + 16));
        }
        highbd_down2_symeven_prepare_4pt(load, filter_2x, &vec4567);

        // P00ab P00cd P01ab P01cd P02ab P02cd P03ab P03cd
        // P04ab P04cd P05ab P05cd P06ab P06cd P07ab P07cd

        // P00 P01 P04 P05 P02 P03 P06 P07
        __m256i vec_8pt = _mm256_hadd_epi32(vec0123, vec4567);
        vec_8pt         = _mm256_permute4x64_epi64(vec_8pt, _MM_SHUFFLE(3, 1, 2, 0));

        vec_8pt = _mm256_add_epi32(vec_8pt, base_sum);
        vec_8pt = _mm256_srai_epi32(vec_8pt, FILTER_BITS);
        mm256_clamp_epi32(&vec_8pt, min, max);

        // 8x 32-bit => 8x 16-bit
        __m128i lo_lane   = _mm256_castsi256_si128(vec_8pt);
        __m128i hi_lane   = _mm256_extracti128_si256(vec_8pt, 1);
        __m128i m128_sum8 = _mm_packus_epi32(lo_lane, hi_lane); // 8x 16-bit

        _mm_storeu_si128((__m128i *)optr, m128_sum8);
    }
    optr += steps;

    *output = optr;
}

void svt_av1_interpolate_core_avx2(const uint8_t *const input, int in_length, uint8_t *output,
                                  int out_length, const int16_t *interp_filters) {
    const int32_t steps = 16;
    const int32_t delta = (((uint32_t)in_length << RS_SCALE_SUBPEL_BITS) + out_length / 2) /
        out_length;
    const int32_t offset = in_length > out_length
               ? (((int32_t)(in_length - out_length) << (RS_SCALE_SUBPEL_BITS - 1)) + out_length / 2) /
            out_length
               : -(((int32_t)(out_length - in_length) << (RS_SCALE_SUBPEL_BITS - 1)) + out_length / 2) /
            out_length;
    uint8_t *optr = output;
    int32_t  x, x1, x2;
    int32_t  y;

    x = 0;
    y = offset + RS_SCALE_EXTRA_OFF;
    while ((y >> RS_SCALE_SUBPEL_BITS) < (SUBPEL_TAPS / 2 - 1)) { // (y >> 14) < 3 ?
        x++;
        y += delta;
    }
    x1 = x;
    x  = out_length - 1;
    y  = delta * x + offset + RS_SCALE_EXTRA_OFF;
    while ((y >> RS_SCALE_SUBPEL_BITS) + (int32_t)(SUBPEL_TAPS / 2) >= in_length) {
        x--;
        y -= delta;
    }
    x2 = x;
    assert(x1 <= x2);
    (void)x1;

    // Initial part.
    x = 0;
    y = offset + RS_SCALE_EXTRA_OFF;
    interpolate_core_w16_init_part_avx2(input, &optr, interp_filters, &y, delta);
    x += steps;

    // Middle part.
    int start = x;
    int mid  = (x2 - start + 1) & (~(steps - 1));
    interpolate_core_w16_mid_part_avx2(input, &optr, interp_filters, &y, delta, mid);
    x += mid;

    // Middle part + End part.
    if (out_length - x >= steps) {
        interpolate_core_w16_end_part_avx2(input, in_length, &optr, interp_filters, &y, delta);
        x += steps;
    }

    // End part.
    for (; x < out_length; ++x, y += delta) {
        int32_t        int_pel = y >> RS_SCALE_SUBPEL_BITS;
        int32_t        sub_pel = (y >> RS_SCALE_EXTRA_BITS) & RS_SUBPEL_MASK;
        const int16_t *filter  = &interp_filters[sub_pel * SUBPEL_TAPS];
        int32_t        sum     = 0;
        for (int32_t k = 0; k < SUBPEL_TAPS; ++k)
            sum += filter[k] *
                input[AOMMIN(int_pel - SUBPEL_TAPS / 2 + 1 + k,
                             in_length - 1)]; // clamp input offset to (any, in_length - 1)
        *optr++ = clip_pixel(ROUND_POWER_OF_TWO(sum, FILTER_BITS));
    }
}

void svt_av1_down2_symeven_avx2(const uint8_t *const input, int length, uint8_t *output) {
    const int16_t *filter          = av1_down2_symeven_half_filter;
    const int      filter_len_half = sizeof(av1_down2_symeven_half_filter) / 2;
    const int      steps           = 32;

    const __m128i  filter_1x = _mm_loadl_epi64((const __m128i *)filter);
    const __m256i  filter_4x = _mm256_broadcastq_epi64(filter_1x);

    int            i, j;
    uint8_t       *optr = output;
    const int      l1   = steps;
    int            l2   = (length - filter_len_half);
    l2 += (l2 & 1);

    // Initial part.
    down2_symeven_w16_init_part_avx2(input, &optr, &filter_4x);
    i = l1;

    // Middle part.
    int mid = (l2 - l1 + 1) & (~(steps-1));
    if (mid > 0) {
        down2_symeven_w16_mid_part_avx2(input + i, mid, &optr, &filter_4x);
        i += mid;
    }

    // Middle part + End part.
    if (length - i >= steps) {
        int len_to_end = length - i;
        down2_symeven_w16_end_part_avx2(input + i, &optr, &filter_4x, len_to_end);
        i += steps;
    }

    // End part.
    for (; i < length; i += 2) {
        int sum = (1 << (FILTER_BITS - 1));
        for (j = 0; j < filter_len_half; ++j) {
            sum += (input[i - j] + input[AOMMIN(i + 1 + j, length - 1)]) * filter[j];
        }
        sum >>= FILTER_BITS;
        *optr++ = clip_pixel(sum);
    }
}

void svt_av1_highbd_interpolate_core_avx2(const uint16_t *const input, int in_length,
                                          uint16_t *output, int out_length, int bd,
                                          const int16_t *interp_filters) {
    const int32_t steps = 8;
    __m256i       max;
    if (bd == 10) {
        max = _mm256_set1_epi32(1023);
    } else {
        max = _mm256_set1_epi32(4095);
    }
    const int32_t delta = (((uint32_t)in_length << RS_SCALE_SUBPEL_BITS) + out_length / 2) /
        out_length;
    const int32_t offset = in_length > out_length
        ? (((int32_t)(in_length - out_length) << (RS_SCALE_SUBPEL_BITS - 1)) + out_length / 2) /
            out_length
        : -(((int32_t)(out_length - in_length) << (RS_SCALE_SUBPEL_BITS - 1)) + out_length / 2) /
            out_length;
    uint16_t *optr = output;
    int32_t   x, x1, x2;
    int32_t   y;

    x = 0;
    y = offset + RS_SCALE_EXTRA_OFF;
    while ((y >> RS_SCALE_SUBPEL_BITS) < (SUBPEL_TAPS / 2 - 1)) { // (y >> 14) < 3 ?
        x++;
        y += delta;
    }
    x1 = x;
    x  = out_length - 1;
    y  = delta * x + offset + RS_SCALE_EXTRA_OFF;
    while ((y >> RS_SCALE_SUBPEL_BITS) + (int32_t)(SUBPEL_TAPS / 2) >= in_length) {
        x--;
        y -= delta;
    }
    x2 = x;
    assert(x1 <= x2);
    (void)x1;

    // Initial part.
    x = 0;
    y = offset + RS_SCALE_EXTRA_OFF;
    highbd_interpolate_core_w8_init_part_avx2(input, &optr, interp_filters, &y, delta, max);
    x += steps;

    // Middle part.
    int start = x;
    int mid2  = (x2 - start + 1) & (~(steps - 1));
    highbd_interpolate_core_w8_mid_part_avx2(input, &optr, interp_filters, &y, delta, mid2, max);
    x += mid2;

    // Middle part + End part.
    if (out_length - x >= steps) {
        highbd_interpolate_core_w8_end_part_avx2(
            input, in_length, &optr, interp_filters, &y, delta, max);
        x += steps;
    }

    // End part.
    for (; x < out_length; ++x, y += delta) {
        int32_t        int_pel = y >> RS_SCALE_SUBPEL_BITS;
        int32_t        sub_pel = (y >> RS_SCALE_EXTRA_BITS) & RS_SUBPEL_MASK;
        const int16_t *filter  = &interp_filters[sub_pel * SUBPEL_TAPS];
        int32_t        sum     = 0;
        for (int32_t k = 0; k < SUBPEL_TAPS; ++k)
            sum += filter[k] *
                input[AOMMIN(int_pel - SUBPEL_TAPS / 2 + 1 + k,
                             in_length - 1)]; // clamp input offset to (any, in_length - 1)
        *optr++ = clip_pixel_highbd(ROUND_POWER_OF_TWO(sum, FILTER_BITS), bd);
    }
}

void svt_av1_highbd_down2_symeven_avx2(const uint16_t *const input, int length, uint16_t *output,
                                       int bd) {
    const int16_t *filter          = av1_down2_symeven_half_filter;
    const int      filter_len_half = sizeof(av1_down2_symeven_half_filter) / 2;
    const int      steps           = 16;
    __m256i        max;
    if (bd == 10) {
        max = _mm256_set1_epi32(1023);
    } else {
        max = _mm256_set1_epi32(4095);
    }

    const __m128i filter_2x_128i = _mm_broadcastq_epi64(_mm_loadl_epi64((const __m128i *)filter));
    const __m256i filter_2x      = _mm256_cvtepi16_epi32(filter_2x_128i);

    int       i, j;
    uint16_t  *optr = output;
    const int l1   = steps;
    int       l2   = (length - filter_len_half);
    l2 += (l2 & 1);

    // Initial part.
    highbd_down2_symeven_w8_init_part_avx2(input, &optr, filter_2x, max);
    i = l1;

    // Middle part.
    int mid = (l2 - l1 + 1) & (~(steps - 1));
    if (mid > 0) {
        highbd_down2_symeven_w8_mid_part_avx2(input + i, mid, &optr, filter_2x, max);
        i += mid;
    }

    // Middle part + End part.
    if (length - i >= steps) {
        int len_to_end = length - i;
        highbd_down2_symeven_w8_end_part_avx2(input + i, &optr, filter_2x, len_to_end, max);
        i += steps;
    }

    // End part.
    for (; i < length; i += 2) {
        int sum = (1 << (FILTER_BITS - 1));
        for (j = 0; j < filter_len_half; ++j) {
            sum += (input[i - j] + input[AOMMIN(i + 1 + j, length - 1)]) * filter[j];
        }
        sum >>= FILTER_BITS;
        *optr++ = clip_pixel_highbd(sum, bd);
    }
}
