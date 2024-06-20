/*
 * Copyright (c) 2023, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <arm_neon.h>
#include <assert.h>
#include "common_dsp_rtcd.h"
#include "mem_neon.h"
#include "definitions.h"

static INLINE uint8x16_t alpha_blend_a64_u8x16(uint8x16_t m, uint8x16_t a, uint8x16_t b) {
    uint16x8_t       blend_u16_lo, blend_u16_hi;
    uint8x8_t        blend_u8_lo, blend_u8_hi;
    const uint8x16_t m_inv = vsubq_u8(vdupq_n_u8(AOM_BLEND_A64_MAX_ALPHA), m);

    blend_u16_lo = vmull_u8(vget_low_u8(m), vget_low_u8(a));
    blend_u16_hi = vmull_u8(vget_high_u8(m), vget_high_u8(a));

    blend_u16_lo = vmlal_u8(blend_u16_lo, vget_low_u8(m_inv), vget_low_u8(b));
    blend_u16_hi = vmlal_u8(blend_u16_hi, vget_high_u8(m_inv), vget_high_u8(b));

    blend_u8_lo = vrshrn_n_u16(blend_u16_lo, AOM_BLEND_A64_ROUND_BITS);
    blend_u8_hi = vrshrn_n_u16(blend_u16_hi, AOM_BLEND_A64_ROUND_BITS);

    return vcombine_u8(blend_u8_lo, blend_u8_hi);
}

static INLINE uint8x8_t alpha_blend_a64_u8x8(uint8x8_t m, uint8x8_t a, uint8x8_t b) {
    uint16x8_t      blend_u16 = vmull_u8(m, a);
    const uint8x8_t m_inv     = vsub_u8(vdup_n_u8(AOM_BLEND_A64_MAX_ALPHA), m);

    blend_u16 = vmlal_u8(blend_u16, m_inv, b);

    return vrshrn_n_u16(blend_u16, AOM_BLEND_A64_ROUND_BITS);
}

static INLINE uint8x8_t avg_blend_u8x8(uint8x8_t a, uint8x8_t b) { return vrhadd_u8(a, b); }

static INLINE uint8x16_t avg_blend_u8x16(uint8x16_t a, uint8x16_t b) { return vrhaddq_u8(a, b); }

static INLINE uint8x8_t avg_blend_pairwise_u8x8(uint8x8_t a, uint8x8_t b) { return vrshr_n_u8(vpadd_u8(a, b), 1); }

static INLINE uint8x16_t avg_blend_pairwise_u8x16(uint8x16_t a, uint8x16_t b) {
    return vrshrq_n_u8(vpaddq_u8(a, b), 1);
}

static INLINE uint8x8_t avg_blend_pairwise_u8x8_4(uint8x8_t a, uint8x8_t b, uint8x8_t c, uint8x8_t d) {
    uint8x8_t a_c = vpadd_u8(a, c);
    uint8x8_t b_d = vpadd_u8(b, d);
    return vrshr_n_u8(vqadd_u8(a_c, b_d), 2);
}

static INLINE uint8x16_t avg_blend_pairwise_u8x16_4(uint8x16_t a, uint8x16_t b, uint8x16_t c, uint8x16_t d) {
    uint8x16_t a_c = vpaddq_u8(a, c);
    uint8x16_t b_d = vpaddq_u8(b, d);
    return vrshrq_n_u8(vqaddq_u8(a_c, b_d), 2);
}

void svt_aom_blend_a64_hmask_neon(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0, uint32_t src0_stride,
                                  const uint8_t *src1, uint32_t src1_stride, const uint8_t *mask, int w, int h) {
    assert(IMPLIES(src0 == dst, src0_stride == dst_stride));
    assert(IMPLIES(src1 == dst, src1_stride == dst_stride));

    assert(h >= 2);
    assert(w >= 2);
    assert(IS_POWER_OF_TWO(h));
    assert(IS_POWER_OF_TWO(w));
    uint8x8_t       tmp0, tmp1;
    uint8x16_t      res_q;
    uint16x8_t      res, res_low, res_high;
    int             i, j;
    const uint8x8_t vdup_64 = vdup_n_u8((uint8_t)64);

    if (w >= 16) {
        const uint8x16_t vdup_64_q = vdupq_n_u8((uint8_t)64);
        for (i = 0; i < h; ++i) {
            for (j = 0; j < w; j += 16) {
                __builtin_prefetch(src0);
                __builtin_prefetch(src1);
                const uint8x16_t tmp0_q        = vld1q_u8(src0);
                const uint8x16_t tmp1_q        = vld1q_u8(src1);
                const uint8x16_t m_q           = vld1q_u8(mask);
                const uint8x16_t max_minus_m_q = vsubq_u8(vdup_64_q, m_q);
                res_low                        = vmull_u8(vget_low_u8(m_q), vget_low_u8(tmp0_q));
                res_low                        = vmlal_u8(res_low, vget_low_u8(max_minus_m_q), vget_low_u8(tmp1_q));
                res_high                       = vmull_u8(vget_high_u8(m_q), vget_high_u8(tmp0_q));
                res_high                       = vmlal_u8(res_high, vget_high_u8(max_minus_m_q), vget_high_u8(tmp1_q));
                res_q                          = vcombine_u8(vrshrn_n_u16(res_low, AOM_BLEND_A64_ROUND_BITS),
                                    vrshrn_n_u16(res_high, AOM_BLEND_A64_ROUND_BITS));
                vst1q_u8(dst, res_q);
                src0 += 16;
                src1 += 16;
                dst += 16;
                mask += 16;
            }
            src0 += src0_stride - w;
            src1 += src1_stride - w;
            dst += dst_stride - w;
            mask -= w;
        }
    } else if (w == 8) {
        const uint8x8_t m           = vld1_u8(mask);
        const uint8x8_t max_minus_m = vsub_u8(vdup_64, m);
        for (i = 0; i < h; ++i) {
            __builtin_prefetch(src0);
            __builtin_prefetch(src1);
            tmp0 = vld1_u8(src0);
            tmp1 = vld1_u8(src1);
            res  = vmull_u8(m, tmp0);
            res  = vmlal_u8(res, max_minus_m, tmp1);
            vst1_u8(dst, vrshrn_n_u16(res, AOM_BLEND_A64_ROUND_BITS));
            src0 += src0_stride;
            src1 += src1_stride;
            dst += dst_stride;
        }
    } else if (w == 4) {
        assert(((uintptr_t)mask & 3) == 0);
        const uint8x8_t m           = vreinterpret_u8_u32(vld1_dup_u32((uint32_t *)mask));
        const uint8x8_t max_minus_m = vsub_u8(vdup_64, m);
        for (i = 0; i < h; i += 2) {
            __builtin_prefetch(src0 + 0 * src0_stride);
            __builtin_prefetch(src0 + 1 * src0_stride);
            __builtin_prefetch(src1 + 0 * src1_stride);
            __builtin_prefetch(src1 + 1 * src1_stride);
            tmp0                   = load_unaligned_u8_4x2(src0, src0_stride);
            tmp1                   = load_unaligned_u8_4x2(src1, src1_stride);
            res                    = vmull_u8(m, tmp0);
            res                    = vmlal_u8(res, max_minus_m, tmp1);
            const uint8x8_t result = vrshrn_n_u16(res, AOM_BLEND_A64_ROUND_BITS);
            store_unaligned_u8_4x1(dst + 0 * dst_stride, result, 0);
            store_unaligned_u8_4x1(dst + 1 * dst_stride, result, 1);
            src0 += (2 * src0_stride);
            src1 += (2 * src1_stride);
            dst += (2 * dst_stride);
        }
    } else if (w == 2) {
        assert(((uintptr_t)mask & 1) == 0);
        const uint8x8_t m           = vreinterpret_u8_u16(vld1_dup_u16((uint16_t *)mask));
        const uint8x8_t max_minus_m = vsub_u8(vdup_64, m);
        for (i = 0; i < h; i += 2) {
            __builtin_prefetch(src0 + 0 * src0_stride);
            __builtin_prefetch(src0 + 1 * src0_stride);
            __builtin_prefetch(src1 + 0 * src1_stride);
            __builtin_prefetch(src1 + 1 * src1_stride);
            tmp0                   = load_unaligned_u8_2x2(src0, src0_stride);
            tmp1                   = load_unaligned_u8_2x2(src1, src1_stride);
            res                    = vmull_u8(m, tmp0);
            res                    = vmlal_u8(res, max_minus_m, tmp1);
            const uint8x8_t result = vrshrn_n_u16(res, AOM_BLEND_A64_ROUND_BITS);
            store_unaligned_u8_2x1(dst + 0 * dst_stride, result, 0);
            store_unaligned_u8_2x1(dst + 1 * dst_stride, result, 1);
            src0 += (2 * src0_stride);
            src1 += (2 * src1_stride);
            dst += (2 * dst_stride);
        }
    }
}

void svt_aom_blend_a64_vmask_neon(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0, uint32_t src0_stride,
                                  const uint8_t *src1, uint32_t src1_stride, const uint8_t *mask, int w, int h) {
    uint8x8_t  tmp0, tmp1;
    uint8x16_t tmp0_q, tmp1_q, res_q;
    uint16x8_t res, res_low, res_high;
    int        i, j;
    assert(IMPLIES(src0 == dst, src0_stride == dst_stride));
    assert(IMPLIES(src1 == dst, src1_stride == dst_stride));

    assert(h >= 2);
    assert(w >= 2);
    assert(IS_POWER_OF_TWO(h));
    assert(IS_POWER_OF_TWO(w));

    if (w >= 16) {
        for (i = 0; i < h; ++i) {
            const uint8x8_t m           = vdup_n_u8((uint8_t)mask[i]);
            const uint8x8_t max_minus_m = vdup_n_u8(64 - (uint8_t)mask[i]);
            for (j = 0; j < w; j += 16) {
                __builtin_prefetch(src0);
                __builtin_prefetch(src1);
                tmp0_q   = vld1q_u8(src0);
                tmp1_q   = vld1q_u8(src1);
                res_low  = vmull_u8(m, vget_low_u8(tmp0_q));
                res_low  = vmlal_u8(res_low, max_minus_m, vget_low_u8(tmp1_q));
                res_high = vmull_u8(m, vget_high_u8(tmp0_q));
                res_high = vmlal_u8(res_high, max_minus_m, vget_high_u8(tmp1_q));
                res_q    = vcombine_u8(vrshrn_n_u16(res_low, AOM_BLEND_A64_ROUND_BITS),
                                    vrshrn_n_u16(res_high, AOM_BLEND_A64_ROUND_BITS));
                vst1q_u8(dst, res_q);
                src0 += 16;
                src1 += 16;
                dst += 16;
            }
            src0 += src0_stride - w;
            src1 += src1_stride - w;
            dst += dst_stride - w;
        }
    } else if (w == 8) {
        for (i = 0; i < h; ++i) {
            __builtin_prefetch(src0);
            __builtin_prefetch(src1);
            const uint8x8_t m           = vdup_n_u8((uint8_t)mask[i]);
            const uint8x8_t max_minus_m = vdup_n_u8(64 - (uint8_t)mask[i]);
            tmp0                        = vld1_u8(src0);
            tmp1                        = vld1_u8(src1);
            res                         = vmull_u8(m, tmp0);
            res                         = vmlal_u8(res, max_minus_m, tmp1);
            vst1_u8(dst, vrshrn_n_u16(res, AOM_BLEND_A64_ROUND_BITS));
            src0 += src0_stride;
            src1 += src1_stride;
            dst += dst_stride;
        }
    } else if (w == 4) {
        for (i = 0; i < h; i += 2) {
            __builtin_prefetch(src0 + 0 * src0_stride);
            __builtin_prefetch(src0 + 1 * src0_stride);
            __builtin_prefetch(src1 + 0 * src1_stride);
            __builtin_prefetch(src1 + 1 * src1_stride);
            const uint16x4_t m1           = vdup_n_u16((uint16_t)mask[i]);
            const uint16x4_t m2           = vdup_n_u16((uint16_t)mask[i + 1]);
            const uint8x8_t  m            = vmovn_u16(vcombine_u16(m1, m2));
            const uint16x4_t max_minus_m1 = vdup_n_u16(64 - (uint16_t)mask[i]);
            const uint16x4_t max_minus_m2 = vdup_n_u16(64 - (uint16_t)mask[i + 1]);
            const uint8x8_t  max_minus_m  = vmovn_u16(vcombine_u16(max_minus_m1, max_minus_m2));
            tmp0                          = load_unaligned_u8_4x2(src0, src0_stride);
            tmp1                          = load_unaligned_u8_4x2(src1, src1_stride);
            res                           = vmull_u8(m, tmp0);
            res                           = vmlal_u8(res, max_minus_m, tmp1);
            const uint8x8_t result        = vrshrn_n_u16(res, AOM_BLEND_A64_ROUND_BITS);
            store_unaligned_u8_4x1(dst + 0 * dst_stride, result, 0);
            store_unaligned_u8_4x1(dst + 1 * dst_stride, result, 1);
            src0 += (2 * src0_stride);
            src1 += (2 * src1_stride);
            dst += (2 * dst_stride);
        }
    } else if (w == 2) {
        for (i = 0; i < h; i += 2) {
            __builtin_prefetch(src0 + 0 * src0_stride);
            __builtin_prefetch(src0 + 1 * src0_stride);
            __builtin_prefetch(src1 + 0 * src1_stride);
            __builtin_prefetch(src1 + 1 * src1_stride);
            const uint8x8_t    m1              = vdup_n_u8(mask[i]);
            const uint8x8_t    m2              = vdup_n_u8(mask[i + 1]);
            const uint16x4x2_t m_trn           = vtrn_u16(vreinterpret_u16_u8(m1), vreinterpret_u16_u8(m2));
            const uint8x8_t    m               = vreinterpret_u8_u16(m_trn.val[0]);
            const uint8x8_t    max_minus_m1    = vdup_n_u8(64 - mask[i]);
            const uint8x8_t    max_minus_m2    = vdup_n_u8(64 - mask[i + 1]);
            const uint16x4x2_t max_minus_m_trn = vtrn_u16(vreinterpret_u16_u8(max_minus_m1),
                                                          vreinterpret_u16_u8(max_minus_m2));
            const uint8x8_t    max_minus_m     = vreinterpret_u8_u16(max_minus_m_trn.val[0]);
            tmp0                               = load_unaligned_u8_2x2(src0, src0_stride);
            tmp1                               = load_unaligned_u8_2x2(src1, src1_stride);
            res                                = vmull_u8(m, tmp0);
            res                                = vmlal_u8(res, max_minus_m, tmp1);
            const uint8x8_t result             = vrshrn_n_u16(res, AOM_BLEND_A64_ROUND_BITS);
            store_unaligned_u8_2x1(dst + 0 * dst_stride, result, 0);
            store_unaligned_u8_2x1(dst + 1 * dst_stride, result, 1);
            src0 += (2 * src0_stride);
            src1 += (2 * src1_stride);
            dst += (2 * dst_stride);
        }
    }
}

static INLINE void blend8x1(int16x8_t mask, int16x8_t src_0, int16x8_t src_1, const int16x8_t v_maxval,
                            int16x8_t *res) {
    int32x4_t       im_res_low, im_res_high;
    const int16x8_t max_minus_mask = vsubq_s16(v_maxval, mask);

    im_res_low = vmull_s16(vget_low_s16(mask), vget_low_s16(src_0));
    im_res_low = vmlal_s16(im_res_low, vget_low_s16(max_minus_mask), vget_low_s16(src_1));

    im_res_high = vmull_s16(vget_high_s16(mask), vget_high_s16(src_0));
    im_res_high = vmlal_s16(im_res_high, vget_high_s16(max_minus_mask), vget_high_s16(src_1));

    *res = vcombine_s16(vshrn_n_s32(im_res_low, AOM_BLEND_A64_ROUND_BITS),
                        vshrn_n_s32(im_res_high, AOM_BLEND_A64_ROUND_BITS));
}

static INLINE void blend_8x4(uint8_t *dst, uint32_t dst_stride, const CONV_BUF_TYPE *src0, uint32_t src0_stride,
                             const CONV_BUF_TYPE *src1, uint32_t src1_stride, int16x8_t mask0, int16x8_t mask1,
                             int16x8_t mask2, int16x8_t mask3, const int16x8_t v_maxval,
                             const uint16x8_t vec_round_offset, const int16x8_t vec_round_bits) {
    int16x8_t src0_0, src0_1, src0_2, src0_3;
    int16x8_t src1_0, src1_1, src1_2, src1_3;
    int16x8_t im_res_0, im_res_1, im_res_2, im_res_3;

    load_s16_8x4((int16_t *)src0, (int32_t)src0_stride, &src0_0, &src0_1, &src0_2, &src0_3);
    load_s16_8x4((int16_t *)src1, (int32_t)src1_stride, &src1_0, &src1_1, &src1_2, &src1_3);

    blend8x1(mask0, src0_0, src1_0, v_maxval, &im_res_0);
    blend8x1(mask1, src0_1, src1_1, v_maxval, &im_res_1);
    blend8x1(mask2, src0_2, src1_2, v_maxval, &im_res_2);
    blend8x1(mask3, src0_3, src1_3, v_maxval, &im_res_3);

    uint16x8_t im_res1_0 = vqsubq_u16(vreinterpretq_u16_s16(im_res_0), vec_round_offset);
    uint16x8_t im_res1_1 = vqsubq_u16(vreinterpretq_u16_s16(im_res_1), vec_round_offset);
    uint16x8_t im_res1_2 = vqsubq_u16(vreinterpretq_u16_s16(im_res_2), vec_round_offset);
    uint16x8_t im_res1_3 = vqsubq_u16(vreinterpretq_u16_s16(im_res_3), vec_round_offset);

    im_res_0 = vshlq_s16(vreinterpretq_s16_u16(im_res1_0), vec_round_bits);
    im_res_1 = vshlq_s16(vreinterpretq_s16_u16(im_res1_1), vec_round_bits);
    im_res_2 = vshlq_s16(vreinterpretq_s16_u16(im_res1_2), vec_round_bits);
    im_res_3 = vshlq_s16(vreinterpretq_s16_u16(im_res1_3), vec_round_bits);

    vst1_u8((dst + 0 * dst_stride), vqmovun_s16(im_res_0));
    vst1_u8((dst + 1 * dst_stride), vqmovun_s16(im_res_1));
    vst1_u8((dst + 2 * dst_stride), vqmovun_s16(im_res_2));
    vst1_u8((dst + 3 * dst_stride), vqmovun_s16(im_res_3));
}

static INLINE void blend_4x4(uint8_t *dst, uint32_t dst_stride, const CONV_BUF_TYPE *src0, uint32_t src0_stride,
                             const CONV_BUF_TYPE *src1, uint32_t src1_stride, int16x4_t mask0, int16x4_t mask1,
                             int16x4_t mask2, int16x4_t mask3, const int16x8_t v_maxval,
                             const uint16x8_t vec_round_offset, const int16x8_t vec_round_bits) {
    int16x8_t  src0_0, src0_1;
    int16x8_t  src1_0, src1_1;
    uint16x8_t tu0 = vdupq_n_u16(0);
    uint16x8_t tu1 = vdupq_n_u16(0);
    uint16x8_t tu2 = vdupq_n_u16(0);
    uint16x8_t tu3 = vdupq_n_u16(0);
    int16x8_t  mask0_1, mask2_3;
    int16x8_t  res0, res1;

    load_unaligned_u16_4x4(src0, src0_stride, &tu0, &tu1);
    load_unaligned_u16_4x4(src1, src1_stride, &tu2, &tu3);

    src0_0 = vreinterpretq_s16_u16(tu0);
    src0_1 = vreinterpretq_s16_u16(tu1);

    src1_0 = vreinterpretq_s16_u16(tu2);
    src1_1 = vreinterpretq_s16_u16(tu3);

    mask0_1 = vcombine_s16(mask0, mask1);
    mask2_3 = vcombine_s16(mask2, mask3);

    blend8x1(mask0_1, src0_0, src1_0, v_maxval, &res0);
    blend8x1(mask2_3, src0_1, src1_1, v_maxval, &res1);

    uint16x8_t im_res_0 = vqsubq_u16(vreinterpretq_u16_s16(res0), vec_round_offset);
    uint16x8_t im_res_1 = vqsubq_u16(vreinterpretq_u16_s16(res1), vec_round_offset);

    src0_0 = vshlq_s16(vreinterpretq_s16_u16(im_res_0), vec_round_bits);
    src0_1 = vshlq_s16(vreinterpretq_s16_u16(im_res_1), vec_round_bits);

    uint8x8_t res_0 = vqmovun_s16(src0_0);
    uint8x8_t res_1 = vqmovun_s16(src0_1);

    store_unaligned_u8_4x1(dst + 0 * dst_stride, res_0, 0);
    store_unaligned_u8_4x1(dst + 1 * dst_stride, res_0, 1);
    store_unaligned_u8_4x1(dst + 2 * dst_stride, res_1, 0);
    store_unaligned_u8_4x1(dst + 3 * dst_stride, res_1, 1);
}

void svt_aom_lowbd_blend_a64_d16_mask_neon(uint8_t *dst, uint32_t dst_stride, const CONV_BUF_TYPE *src0,
                                           uint32_t src0_stride, const CONV_BUF_TYPE *src1, uint32_t src1_stride,
                                           const uint8_t *mask, uint32_t mask_stride, int w, int h, int subw, int subh,
                                           ConvolveParams *conv_params) {
    int                  i        = 0;
    const int            bd       = 8;
    int                  w_tmp    = w;
    const uint8_t       *mask_tmp = mask;
    const CONV_BUF_TYPE *src0_tmp = src0;
    const CONV_BUF_TYPE *src1_tmp = src1;
    uint8_t             *dst_tmp  = dst;

    const int offset_bits  = bd + 2 * FILTER_BITS - conv_params->round_0;
    const int round_offset = (1 << (offset_bits - conv_params->round_1)) +
        (1 << (offset_bits - conv_params->round_1 - 1));
    const int round_bits = 2 * FILTER_BITS - conv_params->round_0 - conv_params->round_1;

    assert(IMPLIES((void *)src0 == dst, src0_stride == dst_stride));
    assert(IMPLIES((void *)src1 == dst, src1_stride == dst_stride));

    assert(h >= 4);
    assert(w >= 4);
    assert(IS_POWER_OF_TWO(h));
    assert(IS_POWER_OF_TWO(w));

    uint8x8_t        s0 = vdup_n_u8(0);
    uint8x8_t        s1 = vdup_n_u8(0);
    uint8x8_t        s2 = vdup_n_u8(0);
    uint8x8_t        s3 = vdup_n_u8(0);
    uint8x16_t       t0, t1, t2, t3, t4, t5, t6, t7;
    int16x8_t        mask0, mask1, mask2, mask3;
    int16x8_t        mask4, mask5, mask6, mask7;
    int32x4_t        m0_32, m1_32, m2_32, m3_32;
    int32x4_t        m4_32, m5_32, m6_32, m7_32;
    uint8x8_t        mask0_l, mask1_l, mask2_l, mask3_l;
    uint8x8_t        mask4_l, mask5_l, mask6_l, mask7_l;
    int16x4_t        mask0_low, mask1_low, mask2_low, mask3_low;
    const uint16x4_t vec_zero       = vdup_n_u16(0);
    const uint16_t   offset         = round_offset - (1 << (round_bits - 1));
    const int16x8_t  v_maxval       = vdupq_n_s16(AOM_BLEND_A64_MAX_ALPHA);
    const int16x8_t  vec_round_bits = vdupq_n_s16(-round_bits);
    const uint16x8_t vec_offset     = vdupq_n_u16(offset);

    if (subw == 0 && subh == 0) {
        if (w_tmp > 7) {
            do {
                w_tmp = w;
                do {
                    load_u8_8x4(mask_tmp, mask_stride, &s0, &s1, &s2, &s3);

                    mask0 = vmovl_s8(vreinterpret_s8_u8(s0));
                    mask1 = vmovl_s8(vreinterpret_s8_u8(s1));
                    mask2 = vmovl_s8(vreinterpret_s8_u8(s2));
                    mask3 = vmovl_s8(vreinterpret_s8_u8(s3));

                    blend_8x4(dst_tmp,
                              dst_stride,
                              src0_tmp,
                              src0_stride,
                              src1_tmp,
                              src1_stride,
                              mask0,
                              mask1,
                              mask2,
                              mask3,
                              v_maxval,
                              vec_offset,
                              vec_round_bits);

                    w_tmp -= 8;
                    mask_tmp += 8;
                    dst_tmp += 8;
                    src0_tmp += 8;
                    src1_tmp += 8;
                } while (w_tmp > 7);
                i += 4;
                mask_tmp += (4 * mask_stride) - w;
                dst_tmp += (4 * dst_stride) - w;
                src0_tmp += (4 * src0_stride) - w;
                src1_tmp += (4 * src1_stride) - w;
            } while (i < h);
        } else {
            do {
                load_unaligned_u8_4x4(mask_tmp, mask_stride, &s0, &s1);

                mask0 = vreinterpretq_s16_u16(vmovl_u8(s0));
                mask1 = vreinterpretq_s16_u16(vmovl_u8(s1));

                mask0_low = vget_low_s16(mask0);
                mask1_low = vget_high_s16(mask0);
                mask2_low = vget_low_s16(mask1);
                mask3_low = vget_high_s16(mask1);

                blend_4x4(dst_tmp,
                          dst_stride,
                          src0_tmp,
                          src0_stride,
                          src1_tmp,
                          src1_stride,
                          mask0_low,
                          mask1_low,
                          mask2_low,
                          mask3_low,
                          v_maxval,
                          vec_offset,
                          vec_round_bits);

                i += 4;
                mask_tmp += (4 * mask_stride);
                dst_tmp += (4 * dst_stride);
                src0_tmp += (4 * src0_stride);
                src1_tmp += (4 * src1_stride);
            } while (i < h);
        }
    } else if (subw == 1 && subh == 1) {
        if (w_tmp > 7) {
            do {
                w_tmp = w;
                do {
                    load_u8_16x8(mask_tmp, mask_stride, &t0, &t1, &t2, &t3, &t4, &t5, &t6, &t7);

                    mask0 = vreinterpretq_s16_u16(vaddl_u8(vget_low_u8(t0), vget_low_u8(t1)));
                    mask1 = vreinterpretq_s16_u16(vaddl_u8(vget_low_u8(t2), vget_low_u8(t3)));
                    mask2 = vreinterpretq_s16_u16(vaddl_u8(vget_low_u8(t4), vget_low_u8(t5)));
                    mask3 = vreinterpretq_s16_u16(vaddl_u8(vget_low_u8(t6), vget_low_u8(t7)));

                    mask4 = vreinterpretq_s16_u16(vaddl_u8(vget_high_u8(t0), vget_high_u8(t1)));
                    mask5 = vreinterpretq_s16_u16(vaddl_u8(vget_high_u8(t2), vget_high_u8(t3)));
                    mask6 = vreinterpretq_s16_u16(vaddl_u8(vget_high_u8(t4), vget_high_u8(t5)));
                    mask7 = vreinterpretq_s16_u16(vaddl_u8(vget_high_u8(t6), vget_high_u8(t7)));

                    m0_32 = vpaddlq_s16(mask0);
                    m1_32 = vpaddlq_s16(mask1);
                    m2_32 = vpaddlq_s16(mask2);
                    m3_32 = vpaddlq_s16(mask3);

                    m4_32 = vpaddlq_s16(mask4);
                    m5_32 = vpaddlq_s16(mask5);
                    m6_32 = vpaddlq_s16(mask6);
                    m7_32 = vpaddlq_s16(mask7);

                    mask0 = vcombine_s16(vqrshrn_n_s32(m0_32, 2), vqrshrn_n_s32(m4_32, 2));
                    mask1 = vcombine_s16(vqrshrn_n_s32(m1_32, 2), vqrshrn_n_s32(m5_32, 2));
                    mask2 = vcombine_s16(vqrshrn_n_s32(m2_32, 2), vqrshrn_n_s32(m6_32, 2));
                    mask3 = vcombine_s16(vqrshrn_n_s32(m3_32, 2), vqrshrn_n_s32(m7_32, 2));

                    blend_8x4(dst_tmp,
                              dst_stride,
                              src0_tmp,
                              src0_stride,
                              src1_tmp,
                              src1_stride,
                              mask0,
                              mask1,
                              mask2,
                              mask3,
                              v_maxval,
                              vec_offset,
                              vec_round_bits);

                    w_tmp -= 8;
                    mask_tmp += 16;
                    dst_tmp += 8;
                    src0_tmp += 8;
                    src1_tmp += 8;
                } while (w_tmp > 7);
                i += 4;
                mask_tmp += (8 * mask_stride) - (2 * w);
                dst_tmp += (4 * dst_stride) - w;
                src0_tmp += (4 * src0_stride) - w;
                src1_tmp += (4 * src1_stride) - w;
            } while (i < h);
        } else {
            do {
                load_u8_8x8(mask_tmp,
                            mask_stride,
                            &mask0_l,
                            &mask1_l,
                            &mask2_l,
                            &mask3_l,
                            &mask4_l,
                            &mask5_l,
                            &mask6_l,
                            &mask7_l);

                mask0 = vreinterpretq_s16_u16(vaddl_u8(mask0_l, mask1_l));
                mask1 = vreinterpretq_s16_u16(vaddl_u8(mask2_l, mask3_l));
                mask2 = vreinterpretq_s16_u16(vaddl_u8(mask4_l, mask5_l));
                mask3 = vreinterpretq_s16_u16(vaddl_u8(mask6_l, mask7_l));

                m0_32 = vpaddlq_s16(mask0);
                m1_32 = vpaddlq_s16(mask1);
                m2_32 = vpaddlq_s16(mask2);
                m3_32 = vpaddlq_s16(mask3);

                mask0_low = vqrshrn_n_s32(m0_32, 2);
                mask1_low = vqrshrn_n_s32(m1_32, 2);
                mask2_low = vqrshrn_n_s32(m2_32, 2);
                mask3_low = vqrshrn_n_s32(m3_32, 2);

                blend_4x4(dst_tmp,
                          dst_stride,
                          src0_tmp,
                          src0_stride,
                          src1_tmp,
                          src1_stride,
                          mask0_low,
                          mask1_low,
                          mask2_low,
                          mask3_low,
                          v_maxval,
                          vec_offset,
                          vec_round_bits);

                i += 4;
                mask_tmp += (8 * mask_stride);
                dst_tmp += (4 * dst_stride);
                src0_tmp += (4 * src0_stride);
                src1_tmp += (4 * src1_stride);
            } while (i < h);
        }
    } else if (subw == 1 && subh == 0) {
        if (w_tmp > 7) {
            do {
                w_tmp = w;
                do {
                    load_u8_16x4(mask_tmp, mask_stride, &t0, &t1, &t2, &t3);

                    mask0 = vreinterpretq_s16_u16(
                        vcombine_u16(vpaddl_u8(vget_low_u8(t0)), vpaddl_u8(vget_high_u8(t0))));
                    mask1 = vreinterpretq_s16_u16(
                        vcombine_u16(vpaddl_u8(vget_low_u8(t1)), vpaddl_u8(vget_high_u8(t1))));
                    mask2 = vreinterpretq_s16_u16(
                        vcombine_u16(vpaddl_u8(vget_low_u8(t2)), vpaddl_u8(vget_high_u8(t2))));
                    mask3 = vreinterpretq_s16_u16(
                        vcombine_u16(vpaddl_u8(vget_low_u8(t3)), vpaddl_u8(vget_high_u8(t3))));

                    mask0 = vmovl_s8(vqrshrn_n_s16(mask0, 1));
                    mask1 = vmovl_s8(vqrshrn_n_s16(mask1, 1));
                    mask2 = vmovl_s8(vqrshrn_n_s16(mask2, 1));
                    mask3 = vmovl_s8(vqrshrn_n_s16(mask3, 1));

                    blend_8x4(dst_tmp,
                              dst_stride,
                              src0_tmp,
                              src0_stride,
                              src1_tmp,
                              src1_stride,
                              mask0,
                              mask1,
                              mask2,
                              mask3,
                              v_maxval,
                              vec_offset,
                              vec_round_bits);
                    w_tmp -= 8;
                    mask_tmp += 16;
                    dst_tmp += 8;
                    src0_tmp += 8;
                    src1_tmp += 8;
                } while (w_tmp > 7);
                i += 4;
                mask_tmp += (4 * mask_stride) - (2 * w);
                dst_tmp += (4 * dst_stride) - w;
                src0_tmp += (4 * src0_stride) - w;
                src1_tmp += (4 * src1_stride) - w;
            } while (i < h);
        } else {
            do {
                load_u8_8x4(mask_tmp, mask_stride, &mask0_l, &mask1_l, &mask2_l, &mask3_l);

                mask0 = vreinterpretq_s16_u16(vcombine_u16(vpaddl_u8(mask0_l), vec_zero));
                mask1 = vreinterpretq_s16_u16(vcombine_u16(vpaddl_u8(mask1_l), vec_zero));
                mask2 = vreinterpretq_s16_u16(vcombine_u16(vpaddl_u8(mask2_l), vec_zero));
                mask3 = vreinterpretq_s16_u16(vcombine_u16(vpaddl_u8(mask3_l), vec_zero));

                mask0_low = vget_low_s16(vmovl_s8(vqrshrn_n_s16(mask0, 1)));
                mask1_low = vget_low_s16(vmovl_s8(vqrshrn_n_s16(mask1, 1)));
                mask2_low = vget_low_s16(vmovl_s8(vqrshrn_n_s16(mask2, 1)));
                mask3_low = vget_low_s16(vmovl_s8(vqrshrn_n_s16(mask3, 1)));

                blend_4x4(dst_tmp,
                          dst_stride,
                          src0_tmp,
                          src0_stride,
                          src1_tmp,
                          src1_stride,
                          mask0_low,
                          mask1_low,
                          mask2_low,
                          mask3_low,
                          v_maxval,
                          vec_offset,
                          vec_round_bits);

                i += 4;
                mask_tmp += (4 * mask_stride);
                dst_tmp += (4 * dst_stride);
                src0_tmp += (4 * src0_stride);
                src1_tmp += (4 * src1_stride);
            } while (i < h);
        }
    } else {
        if (w_tmp > 7) {
            do {
                w_tmp = w;
                do {
                    load_u8_8x8(mask_tmp,
                                mask_stride,
                                &mask0_l,
                                &mask1_l,
                                &mask2_l,
                                &mask3_l,
                                &mask4_l,
                                &mask5_l,
                                &mask6_l,
                                &mask7_l);

                    mask0 = vreinterpretq_s16_u16(vaddl_u8(mask0_l, mask1_l));
                    mask1 = vreinterpretq_s16_u16(vaddl_u8(mask2_l, mask3_l));
                    mask2 = vreinterpretq_s16_u16(vaddl_u8(mask4_l, mask5_l));
                    mask3 = vreinterpretq_s16_u16(vaddl_u8(mask6_l, mask7_l));

                    mask0 = vmovl_s8(vqrshrn_n_s16(mask0, 1));
                    mask1 = vmovl_s8(vqrshrn_n_s16(mask1, 1));
                    mask2 = vmovl_s8(vqrshrn_n_s16(mask2, 1));
                    mask3 = vmovl_s8(vqrshrn_n_s16(mask3, 1));

                    blend_8x4(dst_tmp,
                              dst_stride,
                              src0_tmp,
                              src0_stride,
                              src1_tmp,
                              src1_stride,
                              mask0,
                              mask1,
                              mask2,
                              mask3,
                              v_maxval,
                              vec_offset,
                              vec_round_bits);

                    w_tmp -= 8;
                    mask_tmp += 8;
                    dst_tmp += 8;
                    src0_tmp += 8;
                    src1_tmp += 8;
                } while (w_tmp > 7);
                i += 4;
                mask_tmp += (8 * mask_stride) - w;
                dst_tmp += (4 * dst_stride) - w;
                src0_tmp += (4 * src0_stride) - w;
                src1_tmp += (4 * src1_stride) - w;
            } while (i < h);
        } else {
            do {
                load_unaligned_u8_4x4(mask_tmp, 2 * mask_stride, &s0, &s1);
                load_unaligned_u8_4x4(mask_tmp + mask_stride, 2 * mask_stride, &s2, &s3);

                mask0 = vreinterpretq_s16_u16(vaddl_u8(s0, s2));
                mask1 = vreinterpretq_s16_u16(vaddl_u8(s1, s3));

                mask0 = vmovl_s8(vqrshrn_n_s16(mask0, 1));
                mask1 = vmovl_s8(vqrshrn_n_s16(mask1, 1));

                mask0_low = vget_low_s16(mask0);
                mask1_low = vget_high_s16(mask0);
                mask2_low = vget_low_s16(mask1);
                mask3_low = vget_high_s16(mask1);

                blend_4x4(dst_tmp,
                          dst_stride,
                          src0_tmp,
                          src0_stride,
                          src1_tmp,
                          src1_stride,
                          mask0_low,
                          mask1_low,
                          mask2_low,
                          mask3_low,
                          v_maxval,
                          vec_offset,
                          vec_round_bits);

                i += 4;
                mask_tmp += (8 * mask_stride);
                dst_tmp += (4 * dst_stride);
                src0_tmp += (4 * src0_stride);
                src1_tmp += (4 * src1_stride);
            } while (i < h);
        }
    }
}

void svt_aom_blend_a64_mask_neon(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0, uint32_t src0_stride,
                                 const uint8_t *src1, uint32_t src1_stride, const uint8_t *mask, uint32_t mask_stride,
                                 int w, int h, int subw, int subh) {
    int i;

    assert(IMPLIES(src0 == dst, src0_stride == dst_stride));
    assert(IMPLIES(src1 == dst, src1_stride == dst_stride));

    assert(h >= 1);
    assert(w >= 1);
    assert(IS_POWER_OF_TWO(h));
    assert(IS_POWER_OF_TWO(w));

    if ((subw | subh) == 0) {
        if (w > 8) {
            do {
                i = 0;
                do {
                    const uint8x16_t m0 = vld1q_u8(mask + i);
                    const uint8x16_t s0 = vld1q_u8(src0 + i);
                    const uint8x16_t s1 = vld1q_u8(src1 + i);

                    const uint8x16_t blend = alpha_blend_a64_u8x16(m0, s0, s1);

                    vst1q_u8(dst + i, blend);
                    i += 16;
                } while (i < w);

                mask += mask_stride;
                src0 += src0_stride;
                src1 += src1_stride;
                dst += dst_stride;
            } while (--h != 0);
        } else if (w == 8) {
            do {
                const uint8x8_t m0 = vld1_u8(mask);
                const uint8x8_t s0 = vld1_u8(src0);
                const uint8x8_t s1 = vld1_u8(src1);

                const uint8x8_t blend = alpha_blend_a64_u8x8(m0, s0, s1);

                vst1_u8(dst, blend);

                mask += mask_stride;
                src0 += src0_stride;
                src1 += src1_stride;
                dst += dst_stride;
            } while (--h != 0);
        } else {
            do {
                const uint8x8_t m0 = load_unaligned_u8_4x2(mask, mask_stride);
                const uint8x8_t s0 = load_unaligned_u8_4x2(src0, src0_stride);
                const uint8x8_t s1 = load_unaligned_u8_4x2(src1, src1_stride);

                const uint8x8_t blend = alpha_blend_a64_u8x8(m0, s0, s1);

                store_unaligned_u8_4x2(dst, dst_stride, blend);

                mask += 2 * mask_stride;
                src0 += 2 * src0_stride;
                src1 += 2 * src1_stride;
                dst += 2 * dst_stride;
                h -= 2;
            } while (h != 0);
        }
    } else if ((subw & subh) == 1) {
        if (w > 8) {
            do {
                i = 0;
                do {
                    const uint8x16_t m0 = vld1q_u8(mask + 0 * mask_stride + 2 * i);
                    const uint8x16_t m1 = vld1q_u8(mask + 1 * mask_stride + 2 * i);
                    const uint8x16_t m2 = vld1q_u8(mask + 0 * mask_stride + 2 * i + 16);
                    const uint8x16_t m3 = vld1q_u8(mask + 1 * mask_stride + 2 * i + 16);
                    const uint8x16_t s0 = vld1q_u8(src0 + i);
                    const uint8x16_t s1 = vld1q_u8(src1 + i);

                    const uint8x16_t m_avg = avg_blend_pairwise_u8x16_4(m0, m1, m2, m3);
                    const uint8x16_t blend = alpha_blend_a64_u8x16(m_avg, s0, s1);

                    vst1q_u8(dst + i, blend);

                    i += 16;
                } while (i < w);

                mask += 2 * mask_stride;
                src0 += src0_stride;
                src1 += src1_stride;
                dst += dst_stride;
            } while (--h != 0);
        } else if (w == 8) {
            do {
                const uint8x8_t m0 = vld1_u8(mask + 0 * mask_stride);
                const uint8x8_t m1 = vld1_u8(mask + 1 * mask_stride);
                const uint8x8_t m2 = vld1_u8(mask + 0 * mask_stride + 8);
                const uint8x8_t m3 = vld1_u8(mask + 1 * mask_stride + 8);
                const uint8x8_t s0 = vld1_u8(src0);
                const uint8x8_t s1 = vld1_u8(src1);

                const uint8x8_t m_avg = avg_blend_pairwise_u8x8_4(m0, m1, m2, m3);
                const uint8x8_t blend = alpha_blend_a64_u8x8(m_avg, s0, s1);

                vst1_u8(dst, blend);

                mask += 2 * mask_stride;
                src0 += src0_stride;
                src1 += src1_stride;
                dst += dst_stride;
            } while (--h != 0);
        } else {
            do {
                const uint8x8_t m0 = vld1_u8(mask + 0 * mask_stride);
                const uint8x8_t m1 = vld1_u8(mask + 1 * mask_stride);
                const uint8x8_t m2 = vld1_u8(mask + 2 * mask_stride);
                const uint8x8_t m3 = vld1_u8(mask + 3 * mask_stride);
                const uint8x8_t s0 = load_unaligned_u8_4x2(src0, src0_stride);
                const uint8x8_t s1 = load_unaligned_u8_4x2(src1, src1_stride);

                const uint8x8_t m_avg = avg_blend_pairwise_u8x8_4(m0, m1, m2, m3);
                const uint8x8_t blend = alpha_blend_a64_u8x8(m_avg, s0, s1);

                store_unaligned_u8_4x2(dst, dst_stride, blend);

                mask += 4 * mask_stride;
                src0 += 2 * src0_stride;
                src1 += 2 * src1_stride;
                dst += 2 * dst_stride;
                h -= 2;
            } while (h != 0);
        }
    } else if (subw == 1 && subh == 0) {
        if (w > 8) {
            do {
                i = 0;

                do {
                    const uint8x16_t m0 = vld1q_u8(mask + 2 * i);
                    const uint8x16_t m1 = vld1q_u8(mask + 2 * i + 16);
                    const uint8x16_t s0 = vld1q_u8(src0 + i);
                    const uint8x16_t s1 = vld1q_u8(src1 + i);

                    const uint8x16_t m_avg = avg_blend_pairwise_u8x16(m0, m1);
                    const uint8x16_t blend = alpha_blend_a64_u8x16(m_avg, s0, s1);

                    vst1q_u8(dst + i, blend);

                    i += 16;
                } while (i < w);

                mask += mask_stride;
                src0 += src0_stride;
                src1 += src1_stride;
                dst += dst_stride;
            } while (--h != 0);
        } else if (w == 8) {
            do {
                const uint8x8_t m0 = vld1_u8(mask);
                const uint8x8_t m1 = vld1_u8(mask + 8);
                const uint8x8_t s0 = vld1_u8(src0);
                const uint8x8_t s1 = vld1_u8(src1);

                const uint8x8_t m_avg = avg_blend_pairwise_u8x8(m0, m1);
                const uint8x8_t blend = alpha_blend_a64_u8x8(m_avg, s0, s1);

                vst1_u8(dst, blend);

                mask += mask_stride;
                src0 += src0_stride;
                src1 += src1_stride;
                dst += dst_stride;
            } while (--h != 0);
        } else {
            do {
                const uint8x8_t m0 = vld1_u8(mask + 0 * mask_stride);
                const uint8x8_t m1 = vld1_u8(mask + 1 * mask_stride);
                const uint8x8_t s0 = load_unaligned_u8_4x2(src0, src0_stride);
                const uint8x8_t s1 = load_unaligned_u8_4x2(src1, src1_stride);

                const uint8x8_t m_avg = avg_blend_pairwise_u8x8(m0, m1);
                const uint8x8_t blend = alpha_blend_a64_u8x8(m_avg, s0, s1);

                store_unaligned_u8_4x2(dst, dst_stride, blend);

                mask += 2 * mask_stride;
                src0 += 2 * src0_stride;
                src1 += 2 * src1_stride;
                dst += 2 * dst_stride;
                h -= 2;
            } while (h != 0);
        }
    } else {
        if (w > 8) {
            do {
                i = 0;
                do {
                    const uint8x16_t m0 = vld1q_u8(mask + 0 * mask_stride + i);
                    const uint8x16_t m1 = vld1q_u8(mask + 1 * mask_stride + i);
                    const uint8x16_t s0 = vld1q_u8(src0 + i);
                    const uint8x16_t s1 = vld1q_u8(src1 + i);

                    const uint8x16_t m_avg = avg_blend_u8x16(m0, m1);
                    const uint8x16_t blend = alpha_blend_a64_u8x16(m_avg, s0, s1);

                    vst1q_u8(dst + i, blend);

                    i += 16;
                } while (i < w);

                mask += 2 * mask_stride;
                src0 += src0_stride;
                src1 += src1_stride;
                dst += dst_stride;
            } while (--h != 0);
        } else if (w == 8) {
            do {
                const uint8x8_t m0 = vld1_u8(mask + 0 * mask_stride);
                const uint8x8_t m1 = vld1_u8(mask + 1 * mask_stride);
                const uint8x8_t s0 = vld1_u8(src0);
                const uint8x8_t s1 = vld1_u8(src1);

                const uint8x8_t m_avg = avg_blend_u8x8(m0, m1);
                const uint8x8_t blend = alpha_blend_a64_u8x8(m_avg, s0, s1);

                vst1_u8(dst, blend);

                mask += 2 * mask_stride;
                src0 += src0_stride;
                src1 += src1_stride;
                dst += dst_stride;
            } while (--h != 0);
        } else {
            do {
                const uint8x8_t m0_2 = load_unaligned_u8_4x2(mask + 0 * mask_stride, 2 * mask_stride);
                const uint8x8_t m1_3 = load_unaligned_u8_4x2(mask + 1 * mask_stride, 2 * mask_stride);
                const uint8x8_t s0   = load_unaligned_u8_4x2(src0, src0_stride);
                const uint8x8_t s1   = load_unaligned_u8_4x2(src1, src1_stride);

                const uint8x8_t m_avg = avg_blend_u8x8(m0_2, m1_3);
                const uint8x8_t blend = alpha_blend_a64_u8x8(m_avg, s0, s1);

                store_unaligned_u8_4x2(dst, dst_stride, blend);

                mask += 4 * mask_stride;
                src0 += 2 * src0_stride;
                src1 += 2 * src1_stride;
                dst += 2 * dst_stride;
                h -= 2;
            } while (h != 0);
        }
    }
}
