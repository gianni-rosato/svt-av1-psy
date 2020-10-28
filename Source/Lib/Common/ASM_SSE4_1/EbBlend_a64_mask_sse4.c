/*
 * Copyright(c) 2019 Intel Corporation
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
 */

#include <assert.h>
#include <smmintrin.h>
#include "EbDefinitions.h"
#include "EbBlend_sse4.h"
#include "common_dsp_rtcd.h"

//////////////////////////////////////////////////////////////////////////////
// No sub-sampling
//////////////////////////////////////////////////////////////////////////////

static void blend_a64_mask_w4_sse4_1(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0,
                                     uint32_t src0_stride, const uint8_t *src1,
                                     uint32_t src1_stride, const uint8_t *mask,
                                     uint32_t mask_stride, int w, int h) {
    (void)w;
    const __m128i v_maxval_b = _mm_set1_epi8(AOM_BLEND_A64_MAX_ALPHA);
    const __m128i _r         = _mm_set1_epi16(1 << (15 - AOM_BLEND_A64_ROUND_BITS));
    do {
        const __m128i v_m0_b  = xx_loadl_32(mask);
        const __m128i v_m1_b  = _mm_sub_epi8(v_maxval_b, v_m0_b);
        const __m128i v_res_b = blend_4_u8(src0, src1, &v_m0_b, &v_m1_b, &_r);
        xx_storel_32(dst, v_res_b);

        dst += dst_stride;
        src0 += src0_stride;
        src1 += src1_stride;
        mask += mask_stride;
    } while (--h);
}

static void blend_a64_mask_w8_sse4_1(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0,
                                     uint32_t src0_stride, const uint8_t *src1,
                                     uint32_t src1_stride, const uint8_t *mask,
                                     uint32_t mask_stride, int w, int h) {
    (void)w;
    const __m128i v_maxval_b = _mm_set1_epi8(AOM_BLEND_A64_MAX_ALPHA);
    const __m128i _r         = _mm_set1_epi16(1 << (15 - AOM_BLEND_A64_ROUND_BITS));
    do {
        const __m128i v_m0_b  = xx_loadl_64(mask);
        const __m128i v_m1_b  = _mm_sub_epi8(v_maxval_b, v_m0_b);
        const __m128i v_res_b = blend_8_u8(src0, src1, &v_m0_b, &v_m1_b, &_r);
        xx_storel_64(dst, v_res_b);

        dst += dst_stride;
        src0 += src0_stride;
        src1 += src1_stride;
        mask += mask_stride;
    } while (--h);
}

static void blend_a64_mask_w16n_sse4_1(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0,
                                       uint32_t src0_stride, const uint8_t *src1,
                                       uint32_t src1_stride, const uint8_t *mask,
                                       uint32_t mask_stride, int w, int h) {
    const __m128i v_maxval_b = _mm_set1_epi8(AOM_BLEND_A64_MAX_ALPHA);
    const __m128i _r         = _mm_set1_epi16(1 << (15 - AOM_BLEND_A64_ROUND_BITS));

    do {
        int c;
        for (c = 0; c < w; c += 16) {
            const __m128i v_m0_b = xx_loadu_128(mask + c);
            const __m128i v_m1_b = _mm_sub_epi8(v_maxval_b, v_m0_b);

            const __m128i v_res_b = blend_16_u8(src0 + c, src1 + c, &v_m0_b, &v_m1_b, &_r);

            xx_storeu_128(dst + c, v_res_b);
        }
        dst += dst_stride;
        src0 += src0_stride;
        src1 += src1_stride;
        mask += mask_stride;
    } while (--h);
}

//////////////////////////////////////////////////////////////////////////////
// Horizontal sub-sampling
//////////////////////////////////////////////////////////////////////////////

static void blend_a64_mask_sx_w4_sse4_1(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0,
                                        uint32_t src0_stride, const uint8_t *src1,
                                        uint32_t src1_stride, const uint8_t *mask,
                                        uint32_t mask_stride, int w, int h) {
    (void)w;

    const __m128i v_shuffle_b = xx_loadu_128(g_blend_a64_mask_shuffle);
    const __m128i v_maxval_b  = _mm_set1_epi8(AOM_BLEND_A64_MAX_ALPHA);
    const __m128i _r          = _mm_set1_epi16(1 << (15 - AOM_BLEND_A64_ROUND_BITS));
    do {
        const __m128i v_r_b    = xx_loadl_64(mask);
        const __m128i v_r0_s_b = _mm_shuffle_epi8(v_r_b, v_shuffle_b);
        const __m128i v_r_lo_b = _mm_unpacklo_epi64(v_r0_s_b, v_r0_s_b);
        const __m128i v_r_hi_b = _mm_unpackhi_epi64(v_r0_s_b, v_r0_s_b);
        const __m128i v_m0_b   = _mm_avg_epu8(v_r_lo_b, v_r_hi_b);
        const __m128i v_m1_b   = _mm_sub_epi8(v_maxval_b, v_m0_b);

        const __m128i v_res_b = blend_4_u8(src0, src1, &v_m0_b, &v_m1_b, &_r);
        xx_storel_32(dst, v_res_b);

        dst += dst_stride;
        src0 += src0_stride;
        src1 += src1_stride;
        mask += mask_stride;
    } while (--h);
}

static void blend_a64_mask_sx_w8_sse4_1(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0,
                                        uint32_t src0_stride, const uint8_t *src1,
                                        uint32_t src1_stride, const uint8_t *mask,
                                        uint32_t mask_stride, int w, int h) {
    (void)w;

    const __m128i v_shuffle_b = xx_loadu_128(g_blend_a64_mask_shuffle);
    const __m128i v_maxval_b  = _mm_set1_epi8(AOM_BLEND_A64_MAX_ALPHA);
    const __m128i _r          = _mm_set1_epi16(1 << (15 - AOM_BLEND_A64_ROUND_BITS));
    do {
        const __m128i v_r_b    = xx_loadu_128(mask);
        const __m128i v_r0_s_b = _mm_shuffle_epi8(v_r_b, v_shuffle_b);
        const __m128i v_r_lo_b = _mm_unpacklo_epi64(v_r0_s_b, v_r0_s_b);
        const __m128i v_r_hi_b = _mm_unpackhi_epi64(v_r0_s_b, v_r0_s_b);
        const __m128i v_m0_b   = _mm_avg_epu8(v_r_lo_b, v_r_hi_b);
        const __m128i v_m1_b   = _mm_sub_epi8(v_maxval_b, v_m0_b);

        const __m128i v_res_b = blend_8_u8(src0, src1, &v_m0_b, &v_m1_b, &_r);

        xx_storel_64(dst, v_res_b);

        dst += dst_stride;
        src0 += src0_stride;
        src1 += src1_stride;
        mask += mask_stride;
    } while (--h);
}

static void blend_a64_mask_sx_w16n_sse4_1(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0,
                                          uint32_t src0_stride, const uint8_t *src1,
                                          uint32_t src1_stride, const uint8_t *mask,
                                          uint32_t mask_stride, int w, int h) {
    const __m128i v_shuffle_b = xx_loadu_128(g_blend_a64_mask_shuffle);
    const __m128i v_maxval_b  = _mm_set1_epi8(AOM_BLEND_A64_MAX_ALPHA);
    const __m128i _r          = _mm_set1_epi16(1 << (15 - AOM_BLEND_A64_ROUND_BITS));

    do {
        int c;
        for (c = 0; c < w; c += 16) {
            const __m128i v_r0_b   = xx_loadu_128(mask + 2 * c);
            const __m128i v_r1_b   = xx_loadu_128(mask + 2 * c + 16);
            const __m128i v_r0_s_b = _mm_shuffle_epi8(v_r0_b, v_shuffle_b);
            const __m128i v_r1_s_b = _mm_shuffle_epi8(v_r1_b, v_shuffle_b);
            const __m128i v_r_lo_b = _mm_unpacklo_epi64(v_r0_s_b, v_r1_s_b);
            const __m128i v_r_hi_b = _mm_unpackhi_epi64(v_r0_s_b, v_r1_s_b);
            const __m128i v_m0_b   = _mm_avg_epu8(v_r_lo_b, v_r_hi_b);
            const __m128i v_m1_b   = _mm_sub_epi8(v_maxval_b, v_m0_b);

            const __m128i v_res_b = blend_16_u8(src0 + c, src1 + c, &v_m0_b, &v_m1_b, &_r);

            xx_storeu_128(dst + c, v_res_b);
        }
        dst += dst_stride;
        src0 += src0_stride;
        src1 += src1_stride;
        mask += mask_stride;
    } while (--h);
}

//////////////////////////////////////////////////////////////////////////////
// Vertical sub-sampling
//////////////////////////////////////////////////////////////////////////////

static void blend_a64_mask_sy_w4_sse4_1(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0,
                                        uint32_t src0_stride, const uint8_t *src1,
                                        uint32_t src1_stride, const uint8_t *mask,
                                        uint32_t mask_stride, int w, int h) {
    (void)w;

    const __m128i v_maxval_b = _mm_set1_epi8(AOM_BLEND_A64_MAX_ALPHA);
    const __m128i _r         = _mm_set1_epi16(1 << (15 - AOM_BLEND_A64_ROUND_BITS));

    do {
        const __m128i v_ra_b = xx_loadl_32(mask);
        const __m128i v_rb_b = xx_loadl_32(mask + mask_stride);
        const __m128i v_m0_b = _mm_avg_epu8(v_ra_b, v_rb_b);
        const __m128i v_m1_b = _mm_sub_epi8(v_maxval_b, v_m0_b);

        const __m128i v_res_b = blend_4_u8(src0, src1, &v_m0_b, &v_m1_b, &_r);

        xx_storel_32(dst, v_res_b);

        dst += dst_stride;
        src0 += src0_stride;
        src1 += src1_stride;
        mask += 2 * mask_stride;
    } while (--h);
}

static void blend_a64_mask_sy_w8_sse4_1(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0,
                                        uint32_t src0_stride, const uint8_t *src1,
                                        uint32_t src1_stride, const uint8_t *mask,
                                        uint32_t mask_stride, int w, int h) {
    (void)w;

    const __m128i v_maxval_b = _mm_set1_epi8(AOM_BLEND_A64_MAX_ALPHA);
    const __m128i _r         = _mm_set1_epi16(1 << (15 - AOM_BLEND_A64_ROUND_BITS));
    do {
        const __m128i v_ra_b  = xx_loadl_64(mask);
        const __m128i v_rb_b  = xx_loadl_64(mask + mask_stride);
        const __m128i v_m0_b  = _mm_avg_epu8(v_ra_b, v_rb_b);
        const __m128i v_m1_b  = _mm_sub_epi8(v_maxval_b, v_m0_b);
        const __m128i v_res_b = blend_8_u8(src0, src1, &v_m0_b, &v_m1_b, &_r);

        xx_storel_64(dst, v_res_b);

        dst += dst_stride;
        src0 += src0_stride;
        src1 += src1_stride;
        mask += 2 * mask_stride;
    } while (--h);
}

static void blend_a64_mask_sy_w16n_sse4_1(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0,
                                          uint32_t src0_stride, const uint8_t *src1,
                                          uint32_t src1_stride, const uint8_t *mask,
                                          uint32_t mask_stride, int w, int h) {
    const __m128i v_maxval_b = _mm_set1_epi8(AOM_BLEND_A64_MAX_ALPHA);
    const __m128i _r         = _mm_set1_epi16(1 << (15 - AOM_BLEND_A64_ROUND_BITS));
    do {
        int c;
        for (c = 0; c < w; c += 16) {
            const __m128i v_ra_b = xx_loadu_128(mask + c);
            const __m128i v_rb_b = xx_loadu_128(mask + c + mask_stride);
            const __m128i v_m0_b = _mm_avg_epu8(v_ra_b, v_rb_b);
            const __m128i v_m1_b = _mm_sub_epi8(v_maxval_b, v_m0_b);

            const __m128i v_res_b = blend_16_u8(src0 + c, src1 + c, &v_m0_b, &v_m1_b, &_r);

            xx_storeu_128(dst + c, v_res_b);
        }
        dst += dst_stride;
        src0 += src0_stride;
        src1 += src1_stride;
        mask += 2 * mask_stride;
    } while (--h);
}

//////////////////////////////////////////////////////////////////////////////
// Horizontal and Vertical sub-sampling
//////////////////////////////////////////////////////////////////////////////

static void blend_a64_mask_sx_sy_w4_sse4_1(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0,
                                           uint32_t src0_stride, const uint8_t *src1,
                                           uint32_t src1_stride, const uint8_t *mask,
                                           uint32_t mask_stride, int w, int h) {
    const __m128i v_shuffle_b = xx_loadu_128(g_blend_a64_mask_shuffle);
    const __m128i v_maxval_b  = _mm_set1_epi8(AOM_BLEND_A64_MAX_ALPHA);
    const __m128i _r          = _mm_set1_epi16(1 << (15 - AOM_BLEND_A64_ROUND_BITS));
    (void)w;

    do {
        const __m128i v_ra_b   = xx_loadl_64(mask);
        const __m128i v_rb_b   = xx_loadl_64(mask + mask_stride);
        const __m128i v_rvs_b  = _mm_add_epi8(v_ra_b, v_rb_b);
        const __m128i v_r_s_b  = _mm_shuffle_epi8(v_rvs_b, v_shuffle_b);
        const __m128i v_r0_s_w = _mm_cvtepu8_epi16(v_r_s_b);
        const __m128i v_r1_s_w = _mm_cvtepu8_epi16(_mm_srli_si128(v_r_s_b, 8));
        const __m128i v_rs_w   = _mm_add_epi16(v_r0_s_w, v_r1_s_w);
        const __m128i v_m0_w   = xx_roundn_epu16(v_rs_w, 2);
        const __m128i v_m0_b   = _mm_packus_epi16(v_m0_w, v_m0_w);
        const __m128i v_m1_b   = _mm_sub_epi8(v_maxval_b, v_m0_b);

        const __m128i v_res_b = blend_4_u8(src0, src1, &v_m0_b, &v_m1_b, &_r);

        xx_storel_32(dst, v_res_b);

        dst += dst_stride;
        src0 += src0_stride;
        src1 += src1_stride;
        mask += 2 * mask_stride;
    } while (--h);
}

static void blend_a64_mask_sx_sy_w8_sse4_1(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0,
                                           uint32_t src0_stride, const uint8_t *src1,
                                           uint32_t src1_stride, const uint8_t *mask,
                                           uint32_t mask_stride, int w, int h) {
    const __m128i v_shuffle_b = xx_loadu_128(g_blend_a64_mask_shuffle);
    const __m128i v_maxval_b  = _mm_set1_epi8(AOM_BLEND_A64_MAX_ALPHA);
    const __m128i _r          = _mm_set1_epi16(1 << (15 - AOM_BLEND_A64_ROUND_BITS));
    (void)w;

    do {
        const __m128i v_ra_b = xx_loadu_128(mask);
        const __m128i v_rb_b = xx_loadu_128(mask + mask_stride);

        const __m128i v_rvs_b  = _mm_add_epi8(v_ra_b, v_rb_b);
        const __m128i v_r_s_b  = _mm_shuffle_epi8(v_rvs_b, v_shuffle_b);
        const __m128i v_r0_s_w = _mm_cvtepu8_epi16(v_r_s_b);
        const __m128i v_r1_s_w = _mm_cvtepu8_epi16(_mm_srli_si128(v_r_s_b, 8));
        const __m128i v_rs_w   = _mm_add_epi16(v_r0_s_w, v_r1_s_w);
        const __m128i v_m0_w   = xx_roundn_epu16(v_rs_w, 2);
        const __m128i v_m0_b   = _mm_packus_epi16(v_m0_w, v_m0_w);
        const __m128i v_m1_b   = _mm_sub_epi8(v_maxval_b, v_m0_b);

        const __m128i v_res_b = blend_8_u8(src0, src1, &v_m0_b, &v_m1_b, &_r);

        xx_storel_64(dst, v_res_b);

        dst += dst_stride;
        src0 += src0_stride;
        src1 += src1_stride;
        mask += 2 * mask_stride;
    } while (--h);
}

static void blend_a64_mask_sx_sy_w16n_sse4_1(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0,
                                             uint32_t src0_stride, const uint8_t *src1,
                                             uint32_t src1_stride, const uint8_t *mask,
                                             uint32_t mask_stride, int w, int h) {
    const __m128i v_zmask_b =
        _mm_set_epi8(0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff);
    const __m128i v_maxval_b = _mm_set1_epi8(AOM_BLEND_A64_MAX_ALPHA);
    const __m128i _r         = _mm_set1_epi16(1 << (15 - AOM_BLEND_A64_ROUND_BITS));
    do {
        int c;
        for (c = 0; c < w; c += 16) {
            const __m128i v_ral_b   = xx_loadu_128(mask + 2 * c);
            const __m128i v_rah_b   = xx_loadu_128(mask + 2 * c + 16);
            const __m128i v_rbl_b   = xx_loadu_128(mask + mask_stride + 2 * c);
            const __m128i v_rbh_b   = xx_loadu_128(mask + mask_stride + 2 * c + 16);
            const __m128i v_rvsl_b  = _mm_add_epi8(v_ral_b, v_rbl_b);
            const __m128i v_rvsh_b  = _mm_add_epi8(v_rah_b, v_rbh_b);
            const __m128i v_rvsal_w = _mm_and_si128(v_rvsl_b, v_zmask_b);
            const __m128i v_rvsah_w = _mm_and_si128(v_rvsh_b, v_zmask_b);
            const __m128i v_rvsbl_w = _mm_and_si128(_mm_srli_si128(v_rvsl_b, 1), v_zmask_b);
            const __m128i v_rvsbh_w = _mm_and_si128(_mm_srli_si128(v_rvsh_b, 1), v_zmask_b);
            const __m128i v_rsl_w   = _mm_add_epi16(v_rvsal_w, v_rvsbl_w);
            const __m128i v_rsh_w   = _mm_add_epi16(v_rvsah_w, v_rvsbh_w);

            const __m128i v_m0l_w = xx_roundn_epu16(v_rsl_w, 2);
            const __m128i v_m0h_w = xx_roundn_epu16(v_rsh_w, 2);
            const __m128i v_m0_b  = _mm_packus_epi16(v_m0l_w, v_m0h_w);
            const __m128i v_m1_b  = _mm_sub_epi8(v_maxval_b, v_m0_b);

            const __m128i v_res_b = blend_16_u8(src0 + c, src1 + c, &v_m0_b, &v_m1_b, &_r);

            xx_storeu_128(dst + c, v_res_b);
        }
        dst += dst_stride;
        src0 += src0_stride;
        src1 += src1_stride;
        mask += 2 * mask_stride;
    } while (--h);
}

//////////////////////////////////////////////////////////////////////////////
// Dispatch
//////////////////////////////////////////////////////////////////////////////

void svt_aom_blend_a64_mask_sse4_1(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0,
                                   uint32_t src0_stride, const uint8_t *src1, uint32_t src1_stride,
                                   const uint8_t *mask, uint32_t mask_stride, int w, int h, int subx,
                                   int suby) {
    typedef void (*BlendFn)(uint8_t * dst,
                             uint32_t       dst_stride,
                             const uint8_t *src0,
                             uint32_t       src0_stride,
                             const uint8_t *src1,
                             uint32_t       src1_stride,
                             const uint8_t *mask,
                             uint32_t       mask_stride,
                             int            w,
                             int            h);
    assert(IMPLIES(src0 == dst, src0_stride == dst_stride));
    assert(IMPLIES(src1 == dst, src1_stride == dst_stride));

    assert(h >= 1);
    assert(w >= 1);
    assert(IS_POWER_OF_TWO(h));
    assert(IS_POWER_OF_TWO(w));

    if (UNLIKELY((h | w) & 3)) { // if (w <= 2 || h <= 2)
        svt_aom_blend_a64_mask_c(dst,
                                 dst_stride,
                                 src0,
                                 src0_stride,
                                 src1,
                                 src1_stride,
                                 mask,
                                 mask_stride,
                                 w,
                                 h,
                                 subx,
                                 suby);
    } else {
        // Dimensions are: width_index X subx X suby
        static const BlendFn blend[3][2][2] = {
            {// w % 16 == 0
             {blend_a64_mask_w16n_sse4_1, blend_a64_mask_sy_w16n_sse4_1},
             {blend_a64_mask_sx_w16n_sse4_1, blend_a64_mask_sx_sy_w16n_sse4_1}},
            {// w == 4
             {blend_a64_mask_w4_sse4_1, blend_a64_mask_sy_w4_sse4_1},
             {blend_a64_mask_sx_w4_sse4_1, blend_a64_mask_sx_sy_w4_sse4_1}},
            {// w == 8
             {blend_a64_mask_w8_sse4_1, blend_a64_mask_sy_w8_sse4_1},
             {blend_a64_mask_sx_w8_sse4_1, blend_a64_mask_sx_sy_w8_sse4_1}}};
        assert(((w >> 2) & 3) < 3);
        blend[(w >> 2) & 3][subx != 0][suby != 0](
            dst, dst_stride, src0, src0_stride, src1, src1_stride, mask, mask_stride, w, h);
    }
}

//////////////////////////////////////////////////////////////////////////////
// No sub-sampling
//////////////////////////////////////////////////////////////////////////////

static INLINE void blend_a64_mask_bn_w4_sse4_1(uint16_t *dst, uint32_t dst_stride,
                                               const uint16_t *src0, uint32_t src0_stride,
                                               const uint16_t *src1, uint32_t src1_stride,
                                               const uint8_t *mask, uint32_t mask_stride, int h,
                                               BlendUnitFn blend) {
    const __m128i v_maxval_w = _mm_set1_epi16(AOM_BLEND_A64_MAX_ALPHA);

    do {
        const __m128i v_m0_b = xx_loadl_32(mask);
        const __m128i v_m0_w = _mm_cvtepu8_epi16(v_m0_b);
        const __m128i v_m1_w = _mm_sub_epi16(v_maxval_w, v_m0_w);

        const __m128i v_res_w = blend(src0, src1, v_m0_w, v_m1_w);

        xx_storel_64(dst, v_res_w);

        dst += dst_stride;
        src0 += src0_stride;
        src1 += src1_stride;
        mask += mask_stride;
    } while (--h);
}

static void blend_a64_mask_b10_w4_sse4_1(uint16_t *dst, uint32_t dst_stride, const uint16_t *src0,
                                         uint32_t src0_stride, const uint16_t *src1,
                                         uint32_t src1_stride, const uint8_t *mask,
                                         uint32_t mask_stride, int w, int h) {
    (void)w;
    blend_a64_mask_bn_w4_sse4_1(
        dst, dst_stride, src0, src0_stride, src1, src1_stride, mask, mask_stride, h, blend_4_b10);
}

static void blend_a64_mask_b12_w4_sse4_1(uint16_t *dst, uint32_t dst_stride, const uint16_t *src0,
                                         uint32_t src0_stride, const uint16_t *src1,
                                         uint32_t src1_stride, const uint8_t *mask,
                                         uint32_t mask_stride, int w, int h) {
    (void)w;
    blend_a64_mask_bn_w4_sse4_1(
        dst, dst_stride, src0, src0_stride, src1, src1_stride, mask, mask_stride, h, blend_4_b12);
}

static INLINE void blend_a64_mask_bn_w8n_sse4_1(uint16_t *dst, uint32_t dst_stride,
                                                const uint16_t *src0, uint32_t src0_stride,
                                                const uint16_t *src1, uint32_t src1_stride,
                                                const uint8_t *mask, uint32_t mask_stride, int w,
                                                int h, BlendUnitFn blend) {
    const __m128i v_maxval_w = _mm_set1_epi16(AOM_BLEND_A64_MAX_ALPHA);

    do {
        int c;
        for (c = 0; c < w; c += 8) {
            const __m128i v_m0_b = xx_loadl_64(mask + c);
            const __m128i v_m0_w = _mm_cvtepu8_epi16(v_m0_b);
            const __m128i v_m1_w = _mm_sub_epi16(v_maxval_w, v_m0_w);

            const __m128i v_res_w = blend(src0 + c, src1 + c, v_m0_w, v_m1_w);

            xx_storeu_128(dst + c, v_res_w);
        }
        dst += dst_stride;
        src0 += src0_stride;
        src1 += src1_stride;
        mask += mask_stride;
    } while (--h);
}

static void blend_a64_mask_b10_w8n_sse4_1(uint16_t *dst, uint32_t dst_stride, const uint16_t *src0,
                                          uint32_t src0_stride, const uint16_t *src1,
                                          uint32_t src1_stride, const uint8_t *mask,
                                          uint32_t mask_stride, int w, int h) {
    blend_a64_mask_bn_w8n_sse4_1(dst,
                                 dst_stride,
                                 src0,
                                 src0_stride,
                                 src1,
                                 src1_stride,
                                 mask,
                                 mask_stride,
                                 w,
                                 h,
                                 blend_8_b10);
}

static void blend_a64_mask_b12_w8n_sse4_1(uint16_t *dst, uint32_t dst_stride, const uint16_t *src0,
                                          uint32_t src0_stride, const uint16_t *src1,
                                          uint32_t src1_stride, const uint8_t *mask,
                                          uint32_t mask_stride, int w, int h) {
    blend_a64_mask_bn_w8n_sse4_1(dst,
                                 dst_stride,
                                 src0,
                                 src0_stride,
                                 src1,
                                 src1_stride,
                                 mask,
                                 mask_stride,
                                 w,
                                 h,
                                 blend_8_b12);
}

//////////////////////////////////////////////////////////////////////////////
// Horizontal sub-sampling
//////////////////////////////////////////////////////////////////////////////

static INLINE void blend_a64_mask_bn_sx_w4_sse4_1(uint16_t *dst, uint32_t dst_stride,
                                                  const uint16_t *src0, uint32_t src0_stride,
                                                  const uint16_t *src1, uint32_t src1_stride,
                                                  const uint8_t *mask, uint32_t mask_stride, int h,
                                                  BlendUnitFn blend) {
    const __m128i v_zmask_b =
        _mm_set_epi8(0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff);
    const __m128i v_maxval_w = _mm_set1_epi16(AOM_BLEND_A64_MAX_ALPHA);

    do {
        const __m128i v_r_b = xx_loadl_64(mask);
        const __m128i v_a_b = _mm_avg_epu8(v_r_b, _mm_srli_si128(v_r_b, 1));

        const __m128i v_m0_w = _mm_and_si128(v_a_b, v_zmask_b);
        const __m128i v_m1_w = _mm_sub_epi16(v_maxval_w, v_m0_w);

        const __m128i v_res_w = blend(src0, src1, v_m0_w, v_m1_w);

        xx_storel_64(dst, v_res_w);

        dst += dst_stride;
        src0 += src0_stride;
        src1 += src1_stride;
        mask += mask_stride;
    } while (--h);
}

static void blend_a64_mask_b10_sx_w4_sse4_1(uint16_t *dst, uint32_t dst_stride,
                                            const uint16_t *src0, uint32_t src0_stride,
                                            const uint16_t *src1, uint32_t src1_stride,
                                            const uint8_t *mask, uint32_t mask_stride, int w,
                                            int h) {
    (void)w;
    blend_a64_mask_bn_sx_w4_sse4_1(
        dst, dst_stride, src0, src0_stride, src1, src1_stride, mask, mask_stride, h, blend_4_b10);
}

static void blend_a64_mask_b12_sx_w4_sse4_1(uint16_t *dst, uint32_t dst_stride,
                                            const uint16_t *src0, uint32_t src0_stride,
                                            const uint16_t *src1, uint32_t src1_stride,
                                            const uint8_t *mask, uint32_t mask_stride, int w,
                                            int h) {
    (void)w;
    blend_a64_mask_bn_sx_w4_sse4_1(
        dst, dst_stride, src0, src0_stride, src1, src1_stride, mask, mask_stride, h, blend_4_b12);
}

static INLINE void blend_a64_mask_bn_sx_w8n_sse4_1(uint16_t *dst, uint32_t dst_stride,
                                                   const uint16_t *src0, uint32_t src0_stride,
                                                   const uint16_t *src1, uint32_t src1_stride,
                                                   const uint8_t *mask, uint32_t mask_stride, int w,
                                                   int h, BlendUnitFn blend) {
    const __m128i v_zmask_b =
        _mm_set_epi8(0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff);
    const __m128i v_maxval_w = _mm_set1_epi16(AOM_BLEND_A64_MAX_ALPHA);

    do {
        int c;
        for (c = 0; c < w; c += 8) {
            const __m128i v_r_b = xx_loadu_128(mask + 2 * c);
            const __m128i v_a_b = _mm_avg_epu8(v_r_b, _mm_srli_si128(v_r_b, 1));

            const __m128i v_m0_w = _mm_and_si128(v_a_b, v_zmask_b);
            const __m128i v_m1_w = _mm_sub_epi16(v_maxval_w, v_m0_w);

            const __m128i v_res_w = blend(src0 + c, src1 + c, v_m0_w, v_m1_w);

            xx_storeu_128(dst + c, v_res_w);
        }
        dst += dst_stride;
        src0 += src0_stride;
        src1 += src1_stride;
        mask += mask_stride;
    } while (--h);
}

static void blend_a64_mask_b10_sx_w8n_sse4_1(uint16_t *dst, uint32_t dst_stride,
                                             const uint16_t *src0, uint32_t src0_stride,
                                             const uint16_t *src1, uint32_t src1_stride,
                                             const uint8_t *mask, uint32_t mask_stride, int w,
                                             int h) {
    blend_a64_mask_bn_sx_w8n_sse4_1(dst,
                                    dst_stride,
                                    src0,
                                    src0_stride,
                                    src1,
                                    src1_stride,
                                    mask,
                                    mask_stride,
                                    w,
                                    h,
                                    blend_8_b10);
}

static void blend_a64_mask_b12_sx_w8n_sse4_1(uint16_t *dst, uint32_t dst_stride,
                                             const uint16_t *src0, uint32_t src0_stride,
                                             const uint16_t *src1, uint32_t src1_stride,
                                             const uint8_t *mask, uint32_t mask_stride, int w,
                                             int h) {
    blend_a64_mask_bn_sx_w8n_sse4_1(dst,
                                    dst_stride,
                                    src0,
                                    src0_stride,
                                    src1,
                                    src1_stride,
                                    mask,
                                    mask_stride,
                                    w,
                                    h,
                                    blend_8_b12);
}

//////////////////////////////////////////////////////////////////////////////
// Vertical sub-sampling
//////////////////////////////////////////////////////////////////////////////

static INLINE void blend_a64_mask_bn_sy_w4_sse4_1(uint16_t *dst, uint32_t dst_stride,
                                                  const uint16_t *src0, uint32_t src0_stride,
                                                  const uint16_t *src1, uint32_t src1_stride,
                                                  const uint8_t *mask, uint32_t mask_stride, int h,
                                                  BlendUnitFn blend) {
    const __m128i v_maxval_w = _mm_set1_epi16(AOM_BLEND_A64_MAX_ALPHA);

    do {
        const __m128i v_ra_b = xx_loadl_32(mask);
        const __m128i v_rb_b = xx_loadl_32(mask + mask_stride);
        const __m128i v_a_b  = _mm_avg_epu8(v_ra_b, v_rb_b);

        const __m128i v_m0_w = _mm_cvtepu8_epi16(v_a_b);
        const __m128i v_m1_w = _mm_sub_epi16(v_maxval_w, v_m0_w);

        const __m128i v_res_w = blend(src0, src1, v_m0_w, v_m1_w);

        xx_storel_64(dst, v_res_w);

        dst += dst_stride;
        src0 += src0_stride;
        src1 += src1_stride;
        mask += 2 * mask_stride;
    } while (--h);
}

static void blend_a64_mask_b10_sy_w4_sse4_1(uint16_t *dst, uint32_t dst_stride,
                                            const uint16_t *src0, uint32_t src0_stride,
                                            const uint16_t *src1, uint32_t src1_stride,
                                            const uint8_t *mask, uint32_t mask_stride, int w,
                                            int h) {
    (void)w;
    blend_a64_mask_bn_sy_w4_sse4_1(
        dst, dst_stride, src0, src0_stride, src1, src1_stride, mask, mask_stride, h, blend_4_b10);
}

static void blend_a64_mask_b12_sy_w4_sse4_1(uint16_t *dst, uint32_t dst_stride,
                                            const uint16_t *src0, uint32_t src0_stride,
                                            const uint16_t *src1, uint32_t src1_stride,
                                            const uint8_t *mask, uint32_t mask_stride, int w,
                                            int h) {
    (void)w;
    blend_a64_mask_bn_sy_w4_sse4_1(
        dst, dst_stride, src0, src0_stride, src1, src1_stride, mask, mask_stride, h, blend_4_b12);
}

static INLINE void blend_a64_mask_bn_sy_w8n_sse4_1(uint16_t *dst, uint32_t dst_stride,
                                                   const uint16_t *src0, uint32_t src0_stride,
                                                   const uint16_t *src1, uint32_t src1_stride,
                                                   const uint8_t *mask, uint32_t mask_stride, int w,
                                                   int h, BlendUnitFn blend) {
    const __m128i v_maxval_w = _mm_set1_epi16(AOM_BLEND_A64_MAX_ALPHA);

    do {
        int c;
        for (c = 0; c < w; c += 8) {
            const __m128i v_ra_b = xx_loadl_64(mask + c);
            const __m128i v_rb_b = xx_loadl_64(mask + c + mask_stride);
            const __m128i v_a_b  = _mm_avg_epu8(v_ra_b, v_rb_b);

            const __m128i v_m0_w = _mm_cvtepu8_epi16(v_a_b);
            const __m128i v_m1_w = _mm_sub_epi16(v_maxval_w, v_m0_w);

            const __m128i v_res_w = blend(src0 + c, src1 + c, v_m0_w, v_m1_w);

            xx_storeu_128(dst + c, v_res_w);
        }
        dst += dst_stride;
        src0 += src0_stride;
        src1 += src1_stride;
        mask += 2 * mask_stride;
    } while (--h);
}

static void blend_a64_mask_b10_sy_w8n_sse4_1(uint16_t *dst, uint32_t dst_stride,
                                             const uint16_t *src0, uint32_t src0_stride,
                                             const uint16_t *src1, uint32_t src1_stride,
                                             const uint8_t *mask, uint32_t mask_stride, int w,
                                             int h) {
    blend_a64_mask_bn_sy_w8n_sse4_1(dst,
                                    dst_stride,
                                    src0,
                                    src0_stride,
                                    src1,
                                    src1_stride,
                                    mask,
                                    mask_stride,
                                    w,
                                    h,
                                    blend_8_b10);
}

static void blend_a64_mask_b12_sy_w8n_sse4_1(uint16_t *dst, uint32_t dst_stride,
                                             const uint16_t *src0, uint32_t src0_stride,
                                             const uint16_t *src1, uint32_t src1_stride,
                                             const uint8_t *mask, uint32_t mask_stride, int w,
                                             int h) {
    blend_a64_mask_bn_sy_w8n_sse4_1(dst,
                                    dst_stride,
                                    src0,
                                    src0_stride,
                                    src1,
                                    src1_stride,
                                    mask,
                                    mask_stride,
                                    w,
                                    h,
                                    blend_8_b12);
}

//////////////////////////////////////////////////////////////////////////////
// Horizontal and Vertical sub-sampling
//////////////////////////////////////////////////////////////////////////////

static INLINE void blend_a64_mask_bn_sx_sy_w4_sse4_1(uint16_t *dst, uint32_t dst_stride,
                                                     const uint16_t *src0, uint32_t src0_stride,
                                                     const uint16_t *src1, uint32_t src1_stride,
                                                     const uint8_t *mask, uint32_t mask_stride,
                                                     int h, BlendUnitFn blend) {
    const __m128i v_zmask_b =
        _mm_set_epi8(0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff);
    const __m128i v_maxval_w = _mm_set1_epi16(AOM_BLEND_A64_MAX_ALPHA);

    do {
        const __m128i v_ra_b   = xx_loadl_64(mask);
        const __m128i v_rb_b   = xx_loadl_64(mask + mask_stride);
        const __m128i v_rvs_b  = _mm_add_epi8(v_ra_b, v_rb_b);
        const __m128i v_rvsa_w = _mm_and_si128(v_rvs_b, v_zmask_b);
        const __m128i v_rvsb_w = _mm_and_si128(_mm_srli_si128(v_rvs_b, 1), v_zmask_b);
        const __m128i v_rs_w   = _mm_add_epi16(v_rvsa_w, v_rvsb_w);

        const __m128i v_m0_w = xx_roundn_epu16(v_rs_w, 2);
        const __m128i v_m1_w = _mm_sub_epi16(v_maxval_w, v_m0_w);

        const __m128i v_res_w = blend(src0, src1, v_m0_w, v_m1_w);

        xx_storel_64(dst, v_res_w);

        dst += dst_stride;
        src0 += src0_stride;
        src1 += src1_stride;
        mask += 2 * mask_stride;
    } while (--h);
}

static void blend_a64_mask_b10_sx_sy_w4_sse4_1(uint16_t *dst, uint32_t dst_stride,
                                               const uint16_t *src0, uint32_t src0_stride,
                                               const uint16_t *src1, uint32_t src1_stride,
                                               const uint8_t *mask, uint32_t mask_stride, int w,
                                               int h) {
    (void)w;
    blend_a64_mask_bn_sx_sy_w4_sse4_1(
        dst, dst_stride, src0, src0_stride, src1, src1_stride, mask, mask_stride, h, blend_4_b10);
}

static void blend_a64_mask_b12_sx_sy_w4_sse4_1(uint16_t *dst, uint32_t dst_stride,
                                               const uint16_t *src0, uint32_t src0_stride,
                                               const uint16_t *src1, uint32_t src1_stride,
                                               const uint8_t *mask, uint32_t mask_stride, int w,
                                               int h) {
    (void)w;
    blend_a64_mask_bn_sx_sy_w4_sse4_1(
        dst, dst_stride, src0, src0_stride, src1, src1_stride, mask, mask_stride, h, blend_4_b12);
}

static INLINE void blend_a64_mask_bn_sx_sy_w8n_sse4_1(uint16_t *dst, uint32_t dst_stride,
                                                      const uint16_t *src0, uint32_t src0_stride,
                                                      const uint16_t *src1, uint32_t src1_stride,
                                                      const uint8_t *mask, uint32_t mask_stride,
                                                      int w, int h, BlendUnitFn blend) {
    const __m128i v_zmask_b =
        _mm_set_epi8(0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff);
    const __m128i v_maxval_w = _mm_set1_epi16(AOM_BLEND_A64_MAX_ALPHA);

    do {
        int c;
        for (c = 0; c < w; c += 8) {
            const __m128i v_ra_b   = xx_loadu_128(mask + 2 * c);
            const __m128i v_rb_b   = xx_loadu_128(mask + 2 * c + mask_stride);
            const __m128i v_rvs_b  = _mm_add_epi8(v_ra_b, v_rb_b);
            const __m128i v_rvsa_w = _mm_and_si128(v_rvs_b, v_zmask_b);
            const __m128i v_rvsb_w = _mm_and_si128(_mm_srli_si128(v_rvs_b, 1), v_zmask_b);
            const __m128i v_rs_w   = _mm_add_epi16(v_rvsa_w, v_rvsb_w);

            const __m128i v_m0_w = xx_roundn_epu16(v_rs_w, 2);
            const __m128i v_m1_w = _mm_sub_epi16(v_maxval_w, v_m0_w);

            const __m128i v_res_w = blend(src0 + c, src1 + c, v_m0_w, v_m1_w);

            xx_storeu_128(dst + c, v_res_w);
        }
        dst += dst_stride;
        src0 += src0_stride;
        src1 += src1_stride;
        mask += 2 * mask_stride;
    } while (--h);
}

static void blend_a64_mask_b10_sx_sy_w8n_sse4_1(uint16_t *dst, uint32_t dst_stride,
                                                const uint16_t *src0, uint32_t src0_stride,
                                                const uint16_t *src1, uint32_t src1_stride,
                                                const uint8_t *mask, uint32_t mask_stride, int w,
                                                int h) {
    blend_a64_mask_bn_sx_sy_w8n_sse4_1(dst,
                                       dst_stride,
                                       src0,
                                       src0_stride,
                                       src1,
                                       src1_stride,
                                       mask,
                                       mask_stride,
                                       w,
                                       h,
                                       blend_8_b10);
}

static void blend_a64_mask_b12_sx_sy_w8n_sse4_1(uint16_t *dst, uint32_t dst_stride,
                                                const uint16_t *src0, uint32_t src0_stride,
                                                const uint16_t *src1, uint32_t src1_stride,
                                                const uint8_t *mask, uint32_t mask_stride, int w,
                                                int h) {
    blend_a64_mask_bn_sx_sy_w8n_sse4_1(dst,
                                       dst_stride,
                                       src0,
                                       src0_stride,
                                       src1,
                                       src1_stride,
                                       mask,
                                       mask_stride,
                                       w,
                                       h,
                                       blend_8_b12);
}

//////////////////////////////////////////////////////////////////////////////
// Dispatch
//////////////////////////////////////////////////////////////////////////////

void svt_aom_highbd_blend_a64_mask_8bit_sse4_1(uint8_t *dst_8, uint32_t dst_stride,
                                               const uint8_t *src0_8,
                                               uint32_t src0_stride, const uint8_t *src1_8,
                                               uint32_t src1_stride, const uint8_t *mask,
                                               uint32_t mask_stride, int w, int h, int subx, int suby,
                                               int bd) {
    typedef void (*BlendFn)(uint16_t * dst,
                             uint32_t        dst_stride,
                             const uint16_t *src0,
                             uint32_t        src0_stride,
                             const uint16_t *src1,
                             uint32_t        src1_stride,
                             const uint8_t * mask,
                             uint32_t        mask_stride,
                             int             w,
                             int             h);

    assert(IMPLIES(src0_8 == dst_8, src0_stride == dst_stride));
    assert(IMPLIES(src1_8 == dst_8, src1_stride == dst_stride));

    assert(h >= 1);
    assert(w >= 1);
    assert(IS_POWER_OF_TWO(h));
    assert(IS_POWER_OF_TWO(w));

    assert(bd == 8 || bd == 10 || bd == 12);
    if (UNLIKELY((h | w) & 3)) { // if (w <= 2 || h <= 2)
        svt_aom_highbd_blend_a64_mask_c(dst_8,
                                        dst_stride,
                                        src0_8,
                                        src0_stride,
                                        src1_8,
                                        src1_stride,
                                        mask,
                                        mask_stride,
                                        w,
                                        h,
                                        subx,
                                        suby,
                                        bd);
    } else {
        uint16_t *const       dst  = (uint16_t *)dst_8;
        const uint16_t *const src0 = (uint16_t *)src0_8;
        const uint16_t *const src1 = (uint16_t *)src1_8;
        // Dimensions are: bd_index X width_index X subx X suby
        static const BlendFn blend[2][2][2][2] = {
            {// bd == 8 or 10
             {// w % 8 == 0
              {blend_a64_mask_b10_w8n_sse4_1, blend_a64_mask_b10_sy_w8n_sse4_1},
              {blend_a64_mask_b10_sx_w8n_sse4_1, blend_a64_mask_b10_sx_sy_w8n_sse4_1}},
             {// w == 4
              {blend_a64_mask_b10_w4_sse4_1, blend_a64_mask_b10_sy_w4_sse4_1},
              {blend_a64_mask_b10_sx_w4_sse4_1, blend_a64_mask_b10_sx_sy_w4_sse4_1}}},
            {// bd == 12
             {// w % 8 == 0
              {blend_a64_mask_b12_w8n_sse4_1, blend_a64_mask_b12_sy_w8n_sse4_1},
              {blend_a64_mask_b12_sx_w8n_sse4_1, blend_a64_mask_b12_sx_sy_w8n_sse4_1}},
             {// w == 4
              {blend_a64_mask_b12_w4_sse4_1, blend_a64_mask_b12_sy_w4_sse4_1},
              {blend_a64_mask_b12_sx_w4_sse4_1, blend_a64_mask_b12_sx_sy_w4_sse4_1}}}};
        blend[bd == 12][(w >> 2) & 1][subx != 0][suby != 0](
            dst, dst_stride, src0, src0_stride, src1, src1_stride, mask, mask_stride, w, h);
    }
}

/*Vertical mask related blend functions*/
static void blend_a64_vmask_w4_sse4_1(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0,
                                      uint32_t src0_stride, const uint8_t *src1,
                                      uint32_t src1_stride, const uint8_t *mask, int w, int h) {
    const __m128i v_maxval_w = _mm_set1_epi16(AOM_BLEND_A64_MAX_ALPHA);

    (void)w;

    do {
        const __m128i v_m0_w = _mm_set1_epi16(*mask);
        const __m128i v_m1_w = _mm_sub_epi16(v_maxval_w, v_m0_w);

        const __m128i v_res_w = blend_4(src0, src1, &v_m0_w, &v_m1_w);

        const __m128i v_res_b = _mm_packus_epi16(v_res_w, v_res_w);

        xx_storel_32(dst, v_res_b);

        dst += dst_stride;
        src0 += src0_stride;
        src1 += src1_stride;
        mask += 1;
    } while (--h);
}

static void blend_a64_vmask_w8_sse4_1(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0,
                                      uint32_t src0_stride, const uint8_t *src1,
                                      uint32_t src1_stride, const uint8_t *mask, int w, int h) {
    const __m128i v_maxval_w = _mm_set1_epi16(AOM_BLEND_A64_MAX_ALPHA);

    (void)w;

    do {
        const __m128i v_m0_w = _mm_set1_epi16(*mask);
        const __m128i v_m1_w = _mm_sub_epi16(v_maxval_w, v_m0_w);

        const __m128i v_res_w = blend_8(src0, src1, &v_m0_w, &v_m1_w);

        const __m128i v_res_b = _mm_packus_epi16(v_res_w, v_res_w);

        xx_storel_64(dst, v_res_b);

        dst += dst_stride;
        src0 += src0_stride;
        src1 += src1_stride;
        mask += 1;
    } while (--h);
}

static void blend_a64_vmask_w16n_sse4_1(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0,
                                        uint32_t src0_stride, const uint8_t *src1,
                                        uint32_t src1_stride, const uint8_t *mask, int w, int h) {
    const __m128i v_maxval_w = _mm_set1_epi16(AOM_BLEND_A64_MAX_ALPHA);

    do {
        int           c;
        const __m128i v_m0_w = _mm_set1_epi16(*mask);
        const __m128i v_m1_w = _mm_sub_epi16(v_maxval_w, v_m0_w);
        for (c = 0; c < w; c += 16) {
            const __m128i v_resl_w = blend_8(src0 + c, src1 + c, &v_m0_w, &v_m1_w);
            const __m128i v_resh_w = blend_8(src0 + c + 8, src1 + c + 8, &v_m0_w, &v_m1_w);

            const __m128i v_res_b = _mm_packus_epi16(v_resl_w, v_resh_w);

            xx_storeu_128(dst + c, v_res_b);
        }
        dst += dst_stride;
        src0 += src0_stride;
        src1 += src1_stride;
        mask += 1;
    } while (--h);
}

//////////////////////////////////////////////////////////////////////////////
// Dispatch
//////////////////////////////////////////////////////////////////////////////

void svt_aom_blend_a64_vmask_sse4_1(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0,
                                    uint32_t src0_stride, const uint8_t *src1, uint32_t src1_stride,
                                    const uint8_t *mask, int w, int h) {
    typedef void (*BlendFn)(uint8_t * dst,
                             uint32_t       dst_stride,
                             const uint8_t *src0,
                             uint32_t       src0_stride,
                             const uint8_t *src1,
                             uint32_t       src1_stride,
                             const uint8_t *mask,
                             int            w,
                             int            h);

    // Dimension: width_index
    static const BlendFn blend[9] = {
        blend_a64_vmask_w16n_sse4_1, // w % 16 == 0
        svt_aom_blend_a64_vmask_c, // w == 1
        svt_aom_blend_a64_vmask_c, // w == 2
        NULL, // INVALID
        blend_a64_vmask_w4_sse4_1, // w == 4
        NULL, // INVALID
        NULL, // INVALID
        NULL, // INVALID
        blend_a64_vmask_w8_sse4_1, // w == 8
    };

    assert(IMPLIES(src0 == dst, src0_stride == dst_stride));
    assert(IMPLIES(src1 == dst, src1_stride == dst_stride));

    assert(h >= 1);
    assert(w >= 1);
    assert(IS_POWER_OF_TWO(h));
    assert(IS_POWER_OF_TWO(w));

    blend[w & 0xf](dst, dst_stride, src0, src0_stride, src1, src1_stride, mask, w, h);
}

//////////////////////////////////////////////////////////////////////////////
// Implementation - No sub-sampling
//////////////////////////////////////////////////////////////////////////////

static INLINE void blend_a64_vmask_bn_w4_sse4_1(uint16_t *dst, uint32_t dst_stride,
                                                const uint16_t *src0, uint32_t src0_stride,
                                                const uint16_t *src1, uint32_t src1_stride,
                                                const uint8_t *mask, int h, BlendUnitFn blend) {
    const __m128i v_maxval_w = _mm_set1_epi16(AOM_BLEND_A64_MAX_ALPHA);

    do {
        const __m128i v_m0_w = _mm_set1_epi16(*mask);
        const __m128i v_m1_w = _mm_sub_epi16(v_maxval_w, v_m0_w);

        const __m128i v_res_w = blend(src0, src1, v_m0_w, v_m1_w);

        xx_storel_64(dst, v_res_w);

        dst += dst_stride;
        src0 += src0_stride;
        src1 += src1_stride;
        mask += 1;
    } while (--h);
}

static void blend_a64_vmask_b10_w4_sse4_1(uint16_t *dst, uint32_t dst_stride, const uint16_t *src0,
                                          uint32_t src0_stride, const uint16_t *src1,
                                          uint32_t src1_stride, const uint8_t *mask, int w, int h) {
    (void)w;
    blend_a64_vmask_bn_w4_sse4_1(
        dst, dst_stride, src0, src0_stride, src1, src1_stride, mask, h, blend_4_b10);
}

static void blend_a64_vmask_b12_w4_sse4_1(uint16_t *dst, uint32_t dst_stride, const uint16_t *src0,
                                          uint32_t src0_stride, const uint16_t *src1,
                                          uint32_t src1_stride, const uint8_t *mask, int w, int h) {
    (void)w;
    blend_a64_vmask_bn_w4_sse4_1(
        dst, dst_stride, src0, src0_stride, src1, src1_stride, mask, h, blend_4_b12);
}

static INLINE void blend_a64_vmask_bn_w8n_sse4_1(uint16_t *dst, uint32_t dst_stride,
                                                 const uint16_t *src0, uint32_t src0_stride,
                                                 const uint16_t *src1, uint32_t src1_stride,
                                                 const uint8_t *mask, int w, int h,
                                                 BlendUnitFn blend) {
    const __m128i v_maxval_w = _mm_set1_epi16(AOM_BLEND_A64_MAX_ALPHA);

    do {
        int           c;
        const __m128i v_m0_w = _mm_set1_epi16(*mask);
        const __m128i v_m1_w = _mm_sub_epi16(v_maxval_w, v_m0_w);
        for (c = 0; c < w; c += 8) {
            const __m128i v_res_w = blend(src0 + c, src1 + c, v_m0_w, v_m1_w);

            xx_storeu_128(dst + c, v_res_w);
        }
        dst += dst_stride;
        src0 += src0_stride;
        src1 += src1_stride;
        mask += 1;
    } while (--h);
}

static void blend_a64_vmask_b10_w8n_sse4_1(uint16_t *dst, uint32_t dst_stride, const uint16_t *src0,
                                           uint32_t src0_stride, const uint16_t *src1,
                                           uint32_t src1_stride, const uint8_t *mask, int w,
                                           int h) {
    blend_a64_vmask_bn_w8n_sse4_1(
        dst, dst_stride, src0, src0_stride, src1, src1_stride, mask, w, h, blend_8_b10);
}

static void blend_a64_vmask_b12_w8n_sse4_1(uint16_t *dst, uint32_t dst_stride, const uint16_t *src0,
                                           uint32_t src0_stride, const uint16_t *src1,
                                           uint32_t src1_stride, const uint8_t *mask, int w,
                                           int h) {
    blend_a64_vmask_bn_w8n_sse4_1(
        dst, dst_stride, src0, src0_stride, src1, src1_stride, mask, w, h, blend_8_b12);
}

//////////////////////////////////////////////////////////////////////////////
// Dispatch
//////////////////////////////////////////////////////////////////////////////

void svt_aom_highbd_blend_a64_vmask_8bit_sse4_1(uint8_t *dst_8, uint32_t dst_stride,
                                                const uint8_t *src0_8,
                                                uint32_t src0_stride, const uint8_t *src1_8,
                                                uint32_t src1_stride, const uint8_t *mask, int w, int h,
                                                int bd) {
    typedef void (*BlendFn)(uint16_t * dst,
                             uint32_t        dst_stride,
                             const uint16_t *src0,
                             uint32_t        src0_stride,
                             const uint16_t *src1,
                             uint32_t        src1_stride,
                             const uint8_t * mask,
                             int             w,
                             int             h);

    assert(IMPLIES(src0_8 == dst_8, src0_stride == dst_stride));
    assert(IMPLIES(src1_8 == dst_8, src1_stride == dst_stride));

    assert(h >= 1);
    assert(w >= 1);
    assert(IS_POWER_OF_TWO(h));
    assert(IS_POWER_OF_TWO(w));

    assert(bd == 8 || bd == 10 || bd == 12);

    if (UNLIKELY((h | w) & 3)) { // if (w <= 2 || h <= 2)
        svt_aom_highbd_blend_a64_vmask_8bit_c(
            dst_8, dst_stride, src0_8, src0_stride, src1_8, src1_stride, mask, w, h, bd);
    } else {
        uint16_t *const       dst  = (uint16_t *)(dst_8); // CONVERT_TO_SHORTPTR(dst_8);
        const uint16_t *const src0 = (uint16_t *)(src0_8); //CONVERT_TO_SHORTPTR(src0_8);
        const uint16_t *const src1 = (uint16_t *)(src1_8); //CONVERT_TO_SHORTPTR(src1_8);
        // Dimensions are: bd_index X width_index
        static const BlendFn blend[2][2] = {{
                                                // bd == 8 or 10
                                                blend_a64_vmask_b10_w8n_sse4_1, // w % 8 == 0
                                                blend_a64_vmask_b10_w4_sse4_1, // w == 4
                                            },
                                            {
                                                // bd == 12
                                                blend_a64_vmask_b12_w8n_sse4_1, // w % 8 == 0
                                                blend_a64_vmask_b12_w4_sse4_1, // w == 4
                                            }};
        blend[bd == 12][(w >> 2) & 1](
            dst, dst_stride, src0, src0_stride, src1, src1_stride, mask, w, h);
    }
}

/*Horizontal related blend functions*/

// To start out, just dispatch to the function using the 2D mask and
// pass mask stride as 0. This can be improved upon if necessary.

void svt_aom_blend_a64_hmask_sse4_1(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0,
                                    uint32_t src0_stride, const uint8_t *src1, uint32_t src1_stride,
                                    const uint8_t *mask, int w, int h) {
    svt_aom_blend_a64_mask_sse4_1(
        dst, dst_stride, src0, src0_stride, src1, src1_stride, mask, 0, w, h, 0, 0);
}

void svt_aom_highbd_blend_a64_hmask_8bit_sse4_1(uint8_t *dst_8, uint32_t dst_stride,
                                                const uint8_t *src0_8,
                                                uint32_t src0_stride, const uint8_t *src1_8,
                                                uint32_t src1_stride, const uint8_t *mask, int w, int h,
                                                int bd) {
    svt_aom_highbd_blend_a64_mask_8bit_sse4_1(
        dst_8, dst_stride, src0_8, src0_stride, src1_8, src1_stride, mask, 0, w, h, 0, 0, bd);
}

void svt_aom_highbd_blend_a64_mask_16bit_sse4_1(uint16_t *dst, uint32_t dst_stride,
                                                const uint16_t *src0,
                                                uint32_t src0_stride, const uint16_t *src1,
                                                uint32_t src1_stride, const uint8_t *mask,
                                                uint32_t mask_stride, int w, int h, int subw, int subh,
                                                int bd) {
    typedef void (*BlendFn)(uint16_t * dst,
                             uint32_t        dst_stride,
                             const uint16_t *src0,
                             uint32_t        src0_stride,
                             const uint16_t *src1,
                             uint32_t        src1_stride,
                             const uint8_t * mask,
                             uint32_t        mask_stride,
                             int             w,
                             int             h);

    assert(IMPLIES(src0 == dst, src0_stride == dst_stride));
    assert(IMPLIES(src1 == dst, src1_stride == dst_stride));

    assert(h >= 1);
    assert(w >= 1);
    assert(IS_POWER_OF_TWO(h));
    assert(IS_POWER_OF_TWO(w));

    assert(bd == 8 || bd == 10 || bd == 12);
    if (UNLIKELY((h | w) & 3)) { // if (w <= 2 || h <= 2)
        svt_aom_highbd_blend_a64_mask_c((uint8_t *)dst,
                                        dst_stride,
                                        (uint8_t *)src0,
                                        src0_stride,
                                        (uint8_t *)src1,
                                        src1_stride,
                                        mask,
                                        mask_stride,
                                        w,
                                        h,
                                        subw,
                                        subh,
                                        bd);
    } else {
        //uint16_t *const dst = CONVERT_TO_SHORTPTR(dst_8);
        //const uint16_t *const src0 = CONVERT_TO_SHORTPTR(src0_8);
        //const uint16_t *const src1 = CONVERT_TO_SHORTPTR(src1_8);
        // Dimensions are: bd_index X width_index X subw X subh
        static const BlendFn blend[2][2][2][2] = {
            {// bd == 8 or 10
             {// w % 8 == 0
              {blend_a64_mask_b10_w8n_sse4_1, blend_a64_mask_b10_sy_w8n_sse4_1},
              {blend_a64_mask_b10_sx_w8n_sse4_1, blend_a64_mask_b10_sx_sy_w8n_sse4_1}},
             {// w == 4
              {blend_a64_mask_b10_w4_sse4_1, blend_a64_mask_b10_sy_w4_sse4_1},
              {blend_a64_mask_b10_sx_w4_sse4_1, blend_a64_mask_b10_sx_sy_w4_sse4_1}}},
            {// bd == 12
             {// w % 8 == 0
              {blend_a64_mask_b12_w8n_sse4_1, blend_a64_mask_b12_sy_w8n_sse4_1},
              {blend_a64_mask_b12_sx_w8n_sse4_1, blend_a64_mask_b12_sx_sy_w8n_sse4_1}},
             {// w == 4
              {blend_a64_mask_b12_w4_sse4_1, blend_a64_mask_b12_sy_w4_sse4_1},
              {blend_a64_mask_b12_sx_w4_sse4_1, blend_a64_mask_b12_sx_sy_w4_sse4_1}}}};
        blend[bd == 12][(w >> 2) & 1][subw != 0][subh != 0](
            dst, dst_stride, src0, src0_stride, src1, src1_stride, mask, mask_stride, w, h);
    }
}
void svt_aom_highbd_blend_a64_hmask_16bit_sse4_1(uint16_t *dst, uint32_t dst_stride,
                                                 const uint16_t *src0,
                                                 uint32_t src0_stride, const uint16_t *src1,
                                                 uint32_t src1_stride, const uint8_t *mask, int w, int h,
                                                 int bd) {
    svt_aom_highbd_blend_a64_mask_16bit_sse4_1(
        dst, dst_stride, src0, src0_stride, src1, src1_stride, mask, 0, w, h, 0, 0, bd);
}

void svt_aom_highbd_blend_a64_vmask_16bit_sse4_1(uint16_t *dst, uint32_t dst_stride,
                                                 const uint16_t *src0,
                                                 uint32_t src0_stride, const uint16_t *src1,
                                                 uint32_t src1_stride, const uint8_t *mask, int w, int h,
                                                 int bd) {
    typedef void (*BlendFn)(uint16_t * dst,
                             uint32_t        dst_stride,
                             const uint16_t *src0,
                             uint32_t        src0_stride,
                             const uint16_t *src1,
                             uint32_t        src1_stride,
                             const uint8_t * mask,
                             int             w,
                             int             h);

    assert(IMPLIES(src0 == dst, src0_stride == dst_stride));
    assert(IMPLIES(src1 == dst, src1_stride == dst_stride));

    assert(h >= 1);
    assert(w >= 1);
    assert(IS_POWER_OF_TWO(h));
    assert(IS_POWER_OF_TWO(w));

    assert(bd == 8 || bd == 10 || bd == 12);

    if (UNLIKELY((h | w) & 3)) { // if (w <= 2 || h <= 2)
        svt_aom_highbd_blend_a64_vmask_16bit_c(
            dst, dst_stride, src0, src0_stride, src1, src1_stride, mask, w, h, bd);
    } else {
        //uint16_t *const dst = CONVERT_TO_SHORTPTR(dst_8);
        //const uint16_t *const src0 = CONVERT_TO_SHORTPTR(src0_8);
        //const uint16_t *const src1 = CONVERT_TO_SHORTPTR(src1_8);
        // Dimensions are: bd_index X width_index
        static const BlendFn blend[2][2] = {{
                                                // bd == 8 or 10
                                                blend_a64_vmask_b10_w8n_sse4_1, // w % 8 == 0
                                                blend_a64_vmask_b10_w4_sse4_1, // w == 4
                                            },
                                            {
                                                // bd == 12
                                                blend_a64_vmask_b12_w8n_sse4_1, // w % 8 == 0
                                                blend_a64_vmask_b12_w4_sse4_1, // w == 4
                                            }};
        blend[bd == 12][(w >> 2) & 1](
            dst, dst_stride, src0, src0_stride, src1, src1_stride, mask, w, h);
    }
}
