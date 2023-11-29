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
#include "EbDefinitions.h"
#include "common_dsp_rtcd.h"
#include "mem_neon.h"

/* ---------------------P R E D I C T I O N   Z 1--------------------------- */

static DECLARE_ALIGNED(16, uint8_t, EvenOddMaskx[8][16]) = {{0, 2, 4, 6, 8, 10, 12, 14, 1, 3, 5, 7, 9, 11, 13, 15},
                                                            {0, 1, 3, 5, 7, 9, 11, 13, 0, 2, 4, 6, 8, 10, 12, 14},
                                                            {0, 0, 2, 4, 6, 8, 10, 12, 0, 0, 3, 5, 7, 9, 11, 13},
                                                            {0, 0, 0, 3, 5, 7, 9, 11, 0, 0, 0, 4, 6, 8, 10, 12},
                                                            {0, 0, 0, 0, 4, 6, 8, 10, 0, 0, 0, 0, 5, 7, 9, 11},
                                                            {0, 0, 0, 0, 0, 5, 7, 9, 0, 0, 0, 0, 0, 6, 8, 10},
                                                            {0, 0, 0, 0, 0, 0, 6, 8, 0, 0, 0, 0, 0, 0, 7, 9},
                                                            {0, 0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 0, 0, 0, 0, 8}};

// Low bit depth functions
static DECLARE_ALIGNED(32, uint8_t, BaseMask[33][32]) = {
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0,    0,    0,    0,    0,    0,    0,    0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0,
     0,    0,    0,    0,    0,    0,    0,    0,    0, 0, 0, 0, 0, 0, 0, 0},
    {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0,
     0,    0,    0,    0,    0,    0,    0,    0,    0,    0, 0, 0, 0, 0, 0, 0},
    {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0,
     0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0, 0, 0, 0, 0, 0},
    {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0,
     0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0, 0, 0, 0, 0},
    {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0,
     0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0, 0, 0, 0},
    {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0,
     0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0, 0, 0},
    {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0,
     0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0, 0},
    {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0,
     0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0},
    {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
     0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0},
    {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
     0xff, 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0},
    {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
     0xff, 0xff, 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0},
    {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
     0xff, 0xff, 0xff, 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0},
    {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
     0xff, 0xff, 0xff, 0xff, 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0},
    {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
     0xff, 0xff, 0xff, 0xff, 0xff, 0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0},
    {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
     0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0,    0,    0,    0,    0,    0,    0,    0,    0,    0},
    {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
     0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0,    0,    0,    0,    0,    0,    0,    0,    0},
    {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
     0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0,    0,    0,    0,    0,    0,    0,    0},
    {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
     0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0,    0,    0,    0,    0,    0,    0},
    {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
     0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0,    0,    0,    0,    0,    0},
    {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
     0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0,    0,    0,    0,    0},
    {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
     0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0,    0,    0,    0},
    {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
     0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0,    0,    0},
    {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
     0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0,    0},
    {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
     0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0},
    {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
     0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff},
};

static INLINE void dr_prediction_z1_HxW_internal_neon_64(int H, int W, uint8x8_t *dst, const uint8_t *above,
                                                         int upsample_above, int dx) {
    const int frac_bits  = 6 - upsample_above;
    const int max_base_x = ((W + H) - 1) << upsample_above;

    assert(dx > 0);
    // pre-filter above pixels
    // store in temp buffers:
    //   above[x] * 32 + 16
    //   above[x+1] - above[x]
    // final pixels will be calculated as:
    //   (above[x] * 32 + 16 + (above[x+1] - above[x]) * shift) >> 5

    const uint16x8_t a16       = vdupq_n_u16(16);
    const uint8x8_t  a_mbase_x = vdup_n_u8(above[max_base_x]);
    const uint8x8_t  v_32      = vdup_n_u8(32);

    int x = dx;
    for (int r = 0; r < W; r++) {
        int base          = x >> frac_bits;
        int base_max_diff = (max_base_x - base) >> upsample_above;
        if (base_max_diff <= 0) {
            for (int i = r; i < W; ++i) {
                dst[i] = a_mbase_x; // save 4 values
            }
            return;
        }

        if (base_max_diff > H)
            base_max_diff = H;

        uint8x8x2_t a01_128;
        uint16x8_t  shift;
        if (upsample_above) {
            a01_128 = vld2_u8(above + base);
            shift   = vdupq_n_u16(((x << upsample_above) & 0x3f) >> 1);
        } else {
            a01_128.val[0] = vld1_u8(above + base);
            a01_128.val[1] = vld1_u8(above + base + 1);
            shift          = vdupq_n_u16((x & 0x3f) >> 1);
        }
        uint16x8_t diff = vsubl_u8(a01_128.val[1], a01_128.val[0]);
        uint16x8_t a32  = vmlal_u8(a16, a01_128.val[0], v_32);
        uint16x8_t res  = vmlaq_u16(a32, diff, shift);

        uint8x8_t mask = vld1_u8(BaseMask[base_max_diff]);
        dst[r]         = vbsl_u8(mask, vshrn_n_u16(res, 5), a_mbase_x);

        x += dx;
    }
}

static INLINE void dr_prediction_z1_4xN_neon(int N, uint8_t *dst, ptrdiff_t stride, const uint8_t *above,
                                             int upsample_above, int dx) {
    uint8x8_t dstvec[16];

    dr_prediction_z1_HxW_internal_neon_64(4, N, dstvec, above, upsample_above, dx);
    for (int i = 0; i < N; i++) { vst1_lane_u32((uint32_t *)(dst + stride * i), vreinterpret_u32_u8(dstvec[i]), 0); }
}

static INLINE void dr_prediction_z1_8xN_neon(int N, uint8_t *dst, ptrdiff_t stride, const uint8_t *above,
                                             int upsample_above, int dx) {
    uint8x8_t dstvec[32];

    dr_prediction_z1_HxW_internal_neon_64(8, N, dstvec, above, upsample_above, dx);
    for (int i = 0; i < N; i++) { vst1_u8(dst + stride * i, dstvec[i]); }
}

static INLINE void dr_prediction_z1_HxW_internal_neon(int H, int W, uint8x16_t *dst, const uint8_t *above,
                                                      int upsample_above, int dx) {
    const int frac_bits  = 6 - upsample_above;
    const int max_base_x = ((W + H) - 1) << upsample_above;

    assert(dx > 0);
    // pre-filter above pixels
    // store in temp buffers:
    //   above[x] * 32 + 16
    //   above[x+1] - above[x]
    // final pixels will be calculated as:
    //   (above[x] * 32 + 16 + (above[x+1] - above[x]) * shift) >> 5

    const uint16x8_t a16       = vdupq_n_u16(16);
    const uint8x16_t a_mbase_x = vdupq_n_u8(above[max_base_x]);
    const uint8x8_t  v_32      = vdup_n_u8(32);
    const uint8x16_t v_zero    = vdupq_n_u8(0);

    int x = dx;
    for (int r = 0; r < W; r++) {
        uint16x8x2_t res;
        uint16x8_t   shift;
        uint8x16_t   a0_128, a1_128;

        int base          = x >> frac_bits;
        int base_max_diff = (max_base_x - base) >> upsample_above;
        if (base_max_diff <= 0) {
            for (int i = r; i < W; ++i) {
                dst[i] = a_mbase_x; // save 4 values
            }
            return;
        }

        if (base_max_diff > H)
            base_max_diff = H;

        if (upsample_above) {
            uint8x8x2_t v_tmp_a0_128 = vld2_u8(above + base);
            a0_128                   = vcombine_u8(v_tmp_a0_128.val[0], v_tmp_a0_128.val[1]);
            a1_128                   = vextq_u8(a0_128, v_zero, 8);
            shift                    = vdupq_n_u16(((x << upsample_above) & 0x3f) >> 1);
        } else {
            a0_128 = vld1q_u8(above + base);
            a1_128 = vld1q_u8(above + base + 1);
            shift  = vdupq_n_u16((x & 0x3f) >> 1);
        }
        uint16x8x2_t diff, a32;
        diff.val[0]       = vsubl_u8(vget_low_u8(a1_128), vget_low_u8(a0_128));
        diff.val[1]       = vsubl_u8(vget_high_u8(a1_128), vget_high_u8(a0_128));
        a32.val[0]        = vmlal_u8(a16, vget_low_u8(a0_128), v_32);
        a32.val[1]        = vmlal_u8(a16, vget_high_u8(a0_128), v_32);
        res.val[0]        = vmlaq_u16(a32.val[0], diff.val[0], shift);
        res.val[1]        = vmlaq_u16(a32.val[1], diff.val[1], shift);
        uint8x16_t v_temp = vcombine_u8(vshrn_n_u16(res.val[0], 5), vshrn_n_u16(res.val[1], 5));

        uint8x16_t mask = vld1q_u8(BaseMask[base_max_diff]);
        dst[r]          = vbslq_u8(mask, v_temp, a_mbase_x);

        x += dx;
    }
}

static INLINE void dr_prediction_z1_16xN_neon(int N, uint8_t *dst, ptrdiff_t stride, const uint8_t *above,
                                              int upsample_above, int dx) {
    uint8x16_t dstvec[64];

    dr_prediction_z1_HxW_internal_neon(16, N, dstvec, above, upsample_above, dx);
    for (int i = 0; i < N; i++) { vst1q_u8(dst + stride * i, dstvec[i]); }
}

static INLINE void dr_prediction_z1_32xN_internal_neon(int N, uint8x16x2_t *dstvec, const uint8_t *above, int dx) {
    // here upsample_above is 0 by design of av1_use_intra_edge_upsample
    const int frac_bits  = 6;
    const int max_base_x = ((32 + N) - 1);

    // pre-filter above pixels
    // store in temp buffers:
    //   above[x] * 32 + 16
    //   above[x+1] - above[x]
    // final pixels will be calculated as:
    //   (above[x] * 32 + 16 + (above[x+1] - above[x]) * shift) >> 5

    const uint8x16_t a_mbase_x = vdupq_n_u8(above[max_base_x]);
    const uint16x8_t a16       = vdupq_n_u16(16);
    const uint8x8_t  v_32      = vdup_n_u8(32);

    int x = dx;
    for (int r = 0; r < N; r++) {
        uint8x16_t res16[2];

        int base          = x >> frac_bits;
        int base_max_diff = (max_base_x - base);
        if (base_max_diff <= 0) {
            for (int i = r; i < N; ++i) {
                dstvec[i].val[0] = a_mbase_x; // save 32 values
                dstvec[i].val[1] = a_mbase_x;
            }
            return;
        }
        if (base_max_diff > 32)
            base_max_diff = 32;

        uint16x8_t shift = vdupq_n_u16((x & 0x3f) >> 1);

        for (int j = 0, jj = 0; j < 32; j += 16, jj++) {
            int mdiff = base_max_diff - j;
            if (mdiff <= 0) {
                res16[jj] = a_mbase_x;
            } else {
                uint16x8x2_t a32, diff, res;
                uint8x16_t   a0_128, a1_128;
                a0_128      = vld1q_u8(above + base + j);
                a1_128      = vld1q_u8(above + base + j + 1);
                diff.val[0] = vsubl_u8(vget_low_u8(a1_128), vget_low_u8(a0_128));
                diff.val[1] = vsubl_u8(vget_high_u8(a1_128), vget_high_u8(a0_128));
                a32.val[0]  = vmlal_u8(a16, vget_low_u8(a0_128), v_32);
                a32.val[1]  = vmlal_u8(a16, vget_high_u8(a0_128), v_32);
                res.val[0]  = vmlaq_u16(a32.val[0], diff.val[0], shift);
                res.val[1]  = vmlaq_u16(a32.val[1], diff.val[1], shift);

                res16[jj] = vcombine_u8(vshrn_n_u16(res.val[0], 5), vshrn_n_u16(res.val[1], 5));
            }
        }

        uint8x16x2_t mask;

        mask.val[0]      = vld1q_u8(BaseMask[base_max_diff]);
        mask.val[1]      = vld1q_u8(BaseMask[base_max_diff] + 16);
        dstvec[r].val[0] = vbslq_u8(mask.val[0], res16[0], a_mbase_x);
        dstvec[r].val[1] = vbslq_u8(mask.val[1], res16[1], a_mbase_x);
        x += dx;
    }
}

static INLINE void dr_prediction_z1_32xN_neon(int N, uint8_t *dst, ptrdiff_t stride, const uint8_t *above, int dx) {
    uint8x16x2_t dstvec[64];
    dr_prediction_z1_32xN_internal_neon(N, dstvec, above, dx);
    for (int i = 0; i < N; i++) {
        vst1q_u8(dst + stride * i, dstvec[i].val[0]);
        vst1q_u8(dst + stride * i + 16, dstvec[i].val[1]);
    }
}

static INLINE void dr_prediction_z1_64xN_neon(int N, uint8_t *dst, ptrdiff_t stride, const uint8_t *above, int dx) {
    // here upsample_above is 0 by design of av1_use_intra_edge_upsample
    const int frac_bits  = 6;
    const int max_base_x = ((64 + N) - 1);

    // pre-filter above pixels
    // store in temp buffers:
    //   above[x] * 32 + 16
    //   above[x+1] - above[x]
    // final pixels will be calculated as:
    //   (above[x] * 32 + 16 + (above[x+1] - above[x]) * shift) >> 5

    const uint16x8_t a16           = vdupq_n_u16(16);
    const uint8x16_t a_mbase_x     = vdupq_n_u8(above[max_base_x]);
    const uint8x16_t max_base_x128 = vdupq_n_u8(max_base_x);
    const uint8x8_t  v_32          = vdup_n_u8(32);
    const uint8x16_t v_zero        = vdupq_n_u8(0);
    const uint8x16_t step          = vdupq_n_u8(16);

    int x = dx;
    for (int r = 0; r < N; r++, dst += stride) {
        int base = x >> frac_bits;
        if (base >= max_base_x) {
            for (int i = r; i < N; ++i) {
                vst1q_u8(dst, a_mbase_x);
                vst1q_u8(dst + 16, a_mbase_x);
                vst1q_u8(dst + 32, a_mbase_x);
                vst1q_u8(dst + 48, a_mbase_x);
                dst += stride;
            }
            return;
        }

        uint16x8_t shift       = vdupq_n_u16((x & 0x3f) >> 1);
        uint8x16_t base_inc128 = vaddq_u8(vdupq_n_u8(base),
                                          vcombine_u8(vcreate_u8(0x0706050403020100), vcreate_u8(0x0F0E0D0C0B0A0908)));

        for (int j = 0; j < 64; j += 16) {
            int mdif = max_base_x - (base + j);
            if (mdif <= 0) {
                vst1q_u8(dst + j, a_mbase_x);
            } else {
                uint16x8x2_t a32, diff, res;
                uint8x16_t   a0_128, a1_128, mask128, res128;
                a0_128            = vld1q_u8(above + base + j);
                a1_128            = vld1q_u8(above + base + 1 + j);
                diff.val[0]       = vsubl_u8(vget_low_u8(a1_128), vget_low_u8(a0_128));
                diff.val[1]       = vsubl_u8(vget_high_u8(a1_128), vget_high_u8(a0_128));
                a32.val[0]        = vmlal_u8(a16, vget_low_u8(a0_128), v_32);
                a32.val[1]        = vmlal_u8(a16, vget_high_u8(a0_128), v_32);
                res.val[0]        = vmlaq_u16(a32.val[0], diff.val[0], shift);
                res.val[1]        = vmlaq_u16(a32.val[1], diff.val[1], shift);
                uint8x16_t v_temp = vcombine_u8(vshrn_n_u16(res.val[0], 5), vshrn_n_u16(res.val[1], 5));

                mask128 = vcgtq_u8(vqsubq_u8(max_base_x128, base_inc128), v_zero);
                res128  = vbslq_u8(mask128, v_temp, a_mbase_x);
                vst1q_u8(dst + j, res128);

                base_inc128 = vaddq_u8(base_inc128, step);
            }
        }
        x += dx;
    }
}

// Directional prediction, zone 1: 0 < angle < 90
void svt_av1_dr_prediction_z1_neon(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint8_t *above,
                                   const uint8_t *left, int32_t upsample_above, int32_t dx, int32_t dy) {
    (void)left;
    (void)dy;

    switch (bw) {
    case 4: {
        dr_prediction_z1_4xN_neon(bh, dst, stride, above, upsample_above, dx);
        break;
    }
    case 8: {
        dr_prediction_z1_8xN_neon(bh, dst, stride, above, upsample_above, dx);
        break;
    }
    case 16: {
        dr_prediction_z1_16xN_neon(bh, dst, stride, above, upsample_above, dx);
        break;
    }
    case 32: {
        dr_prediction_z1_32xN_neon(bh, dst, stride, above, dx);
        break;
    }
    case 64: {
        dr_prediction_z1_64xN_neon(bh, dst, stride, above, dx);
        break;
    }
    }
}

/* ---------------------P R E D I C T I O N   Z 2--------------------------- */

static DECLARE_ALIGNED(16, uint8_t, LoadMaskz2[4][16]) = {
    {0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0},
    {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0},
    {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff}};

static INLINE void vector_shift_x4(uint8x8_t *vec, uint8x8_t *v_zero, int shift_value) {
    switch (shift_value) {
    case 1: {
        *vec = vext_u8(*v_zero, *vec, 7);
        break;
    }
    case 2: {
        *vec = vext_u8(*v_zero, *vec, 6);
        break;
    }
    case 3: {
        *vec = vext_u8(*v_zero, *vec, 5);
        break;
    }
    }
}

static INLINE void dr_prediction_z2_Nx4_neon(int N, uint8_t *dst, ptrdiff_t stride, const uint8_t *above,
                                             const uint8_t *left, int upsample_above, int upsample_left, int dx,
                                             int dy) {
    const int min_base_x  = -(1 << upsample_above);
    const int min_base_y  = -(1 << upsample_left);
    const int frac_bits_x = 6 - upsample_above;
    const int frac_bits_y = 6 - upsample_left;

    assert(dx > 0);
    // pre-filter above pixels
    // store in temp buffers:
    //   above[x] * 32 + 16
    //   above[x+1] - above[x]
    // final pixels will be calculated as:
    //   (above[x] * 32 + 16 + (above[x+1] - above[x]) * shift) >> 5
    uint16x8_t a0_x, a1_x, a32, diff;
    uint16x8_t v_32   = vdupq_n_u16(32);
    uint16x8_t v_zero = vdupq_n_u16(0);
    uint16x8_t a16    = vdupq_n_u16(16);

    uint8x8_t  v_zero_u8        = vdup_n_u8(0);
    uint16x4_t v_c3f            = vdup_n_u16(0x3f);
    uint16x4_t r6               = vcreate_u16(0x00C0008000400000);
    int16x4_t  v_upsample_left  = vdup_n_s16(upsample_left);
    int16x4_t  v_upsample_above = vdup_n_s16(upsample_above);
    int16x4_t  v_1234           = vcreate_s16(0x0004000300020001);
    int16x4_t  dy64             = vdup_n_s16(dy);
    int16x4_t  v_frac_bits_y    = vdup_n_s16(-frac_bits_y);
    int16x4_t  min_base_y64     = vdup_n_s16(min_base_y);

    // Use ext rather than loading left + 14 directly to avoid over-read.
    const uint8x16_t   left_m2   = vld1q_u8(left - 2);
    const uint8x16_t   left_0    = vld1q_u8(left);
    const uint8x16_t   left_14   = vextq_u8(left_0, left_0, 14);
    const uint8x16x2_t left_vals = {{left_m2, left_14}};

    for (int r = 0; r < N; r++) {
        uint16x8_t   res, shift;
        uint8x8_t    resx, resy;
        uint16x4x2_t v_shift;
        v_shift.val[1] = vdup_n_u16(0);
        int y          = r + 1;
        int base_x     = (-y * dx) >> frac_bits_x;
        int base_shift = 0;
        if (base_x < (min_base_x - 1)) {
            base_shift = (min_base_x - base_x - 1) >> upsample_above;
        }
        int base_min_diff = (min_base_x - base_x + upsample_above) >> upsample_above;
        if (base_min_diff > 4) {
            base_min_diff = 4;
        } else {
            if (base_min_diff < 0)
                base_min_diff = 0;
        }

        if (base_shift > 3) {
            a0_x           = v_zero;
            a1_x           = v_zero;
            v_shift.val[0] = vreinterpret_u16_u8(v_zero_u8);
            v_shift.val[1] = vreinterpret_u16_u8(v_zero_u8);
        } else {
            uint16x4_t ydx = vdup_n_u16(y * dx);

            if (upsample_above) {
                uint8x8x2_t v_tmp;
                v_tmp.val[0]           = vld1_u8(above + base_x + base_shift);
                v_tmp.val[1]           = vld1_u8(above + base_x + base_shift + 8);
                uint8x8_t v_index_low  = vld1_u8(EvenOddMaskx[base_shift]);
                uint8x8_t v_index_high = vld1_u8(EvenOddMaskx[base_shift] + 8);
                a0_x                   = vmovl_u8(vtbl2_u8(v_tmp, v_index_low));
                a1_x                   = vmovl_u8(vtbl2_u8(v_tmp, v_index_high));
                v_shift.val[0]         = vshr_n_u16(vand_u16(vshl_u16(vsub_u16(r6, ydx), v_upsample_above), v_c3f), 1);
            } else {
                uint8x8_t v_a0_x64 = vld1_u8(above + base_x + base_shift);
                vector_shift_x4(&v_a0_x64, &v_zero_u8, base_shift);
                uint8x8_t v_a1_x64 = vext_u8(v_a0_x64, v_zero_u8, 1);
                v_shift.val[0]     = vshr_n_u16(vand_u16(vsub_u16(r6, ydx), v_c3f), 1);
                a0_x               = vmovl_u8(v_a0_x64);
                a1_x               = vmovl_u8(v_a1_x64);
            }
        }

        // y calc
        if (base_x < min_base_x) {
            int16x4_t  v_r6       = vdup_n_s16(r << 6);
            int16x4_t  y_c64      = vmls_s16(v_r6, v_1234, dy64);
            int16x4_t  base_y_c64 = vshl_s16(y_c64, v_frac_bits_y);
            uint16x4_t mask64     = vcgt_s16(min_base_y64, base_y_c64);

            // Values in base_y_c64 range from -2 through 14 inclusive.
            base_y_c64 = vbic_s16(base_y_c64, vreinterpret_s16_u16(mask64));

            uint8x8_t left_idx0 = vreinterpret_u8_s16(base_y_c64 + 2); // [0, 16]
            uint8x8_t left_idx1 = vreinterpret_u8_s16(base_y_c64 + 3); // [1, 17]

            uint8x8_t a0_y = vtrn1_u8(vqtbl2_u8(left_vals, left_idx0), v_zero_u8);
            uint8x8_t a1_y = vtrn1_u8(vqtbl2_u8(left_vals, left_idx1), v_zero_u8);

            if (upsample_left) {
                v_shift.val[1] = vshr_n_u16(vand_u16(vshl_u16(vreinterpret_u16_s16(y_c64), v_upsample_left), v_c3f), 1);
            } else {
                v_shift.val[1] = vshr_n_u16(vand_u16(vreinterpret_u16_s16(y_c64), v_c3f), 1);
            }

            a0_x = vcombine_u16(vget_low_u16(a0_x), vreinterpret_u16_u8(a0_y));
            a1_x = vcombine_u16(vget_low_u16(a1_x), vreinterpret_u16_u8(a1_y));
        }
        shift = vcombine_u16(v_shift.val[0], v_shift.val[1]);
        diff  = vsubq_u16(a1_x, a0_x); // a[x+1] - a[x]
        a32   = vmlaq_u16(a16, a0_x, v_32); // a[x] * 32 + 16
        res   = vmlaq_u16(a32, diff, shift);
        resx  = vshrn_n_u16(res, 5);
        resy  = vext_u8(resx, v_zero_u8, 4);

        uint8x8_t mask    = vld1_u8(BaseMask[base_min_diff]);
        uint8x8_t v_resxy = vbsl_u8(mask, resy, resx);
        vst1_lane_u32((uint32_t *)dst, vreinterpret_u32_u8(v_resxy), 0);

        dst += stride;
    }
}

static INLINE void vector_shuffle(uint8x16_t *vec, uint8x16_t *vzero, int shift_value) {
    switch (shift_value) {
    case 1: {
        *vec = vextq_u8(*vzero, *vec, 15);
        break;
    }
    case 2: {
        *vec = vextq_u8(*vzero, *vec, 14);
        break;
    }
    case 3: {
        *vec = vextq_u8(*vzero, *vec, 13);
        break;
    }
    case 4: {
        *vec = vextq_u8(*vzero, *vec, 12);
        break;
    }
    case 5: {
        *vec = vextq_u8(*vzero, *vec, 11);
        break;
    }
    case 6: {
        *vec = vextq_u8(*vzero, *vec, 10);
        break;
    }
    case 7: {
        *vec = vextq_u8(*vzero, *vec, 9);
        break;
    }
    case 8: {
        *vec = vextq_u8(*vzero, *vec, 8);
        break;
    }
    case 9: {
        *vec = vextq_u8(*vzero, *vec, 7);
        break;
    }
    case 10: {
        *vec = vextq_u8(*vzero, *vec, 6);
        break;
    }
    case 11: {
        *vec = vextq_u8(*vzero, *vec, 5);
        break;
    }
    case 12: {
        *vec = vextq_u8(*vzero, *vec, 4);
        break;
    }
    case 13: {
        *vec = vextq_u8(*vzero, *vec, 3);
        break;
    }
    case 14: {
        *vec = vextq_u8(*vzero, *vec, 2);
        break;
    }
    case 15: {
        *vec = vextq_u8(*vzero, *vec, 1);
        break;
    }
    }
}

static INLINE void dr_prediction_z2_Nx8_neon(int N, uint8_t *dst, ptrdiff_t stride, const uint8_t *above,
                                             const uint8_t *left, int upsample_above, int upsample_left, int dx,
                                             int dy) {
    const int min_base_x  = -(1 << upsample_above);
    const int min_base_y  = -(1 << upsample_left);
    const int frac_bits_x = 6 - upsample_above;
    const int frac_bits_y = 6 - upsample_left;

    // pre-filter above pixels
    // store in temp buffers:
    //   above[x] * 32 + 16
    //   above[x+1] - above[x]
    // final pixels will be calculated as:
    //   (above[x] * 32 + 16 + (above[x+1] - above[x]) * shift) >> 5
    uint16x8x2_t diff, a32;
    uint8x16_t   v_zero           = vdupq_n_u8(0);
    int16x8_t    v_upsample_left  = vdupq_n_s16(upsample_left);
    int16x8_t    v_upsample_above = vdupq_n_s16(upsample_above);
    int16x8_t    v_frac_bits_y    = vdupq_n_s16(-frac_bits_y);

    uint16x8_t a16           = vdupq_n_u16(16);
    uint16x8_t c3f           = vdupq_n_u16(0x3f);
    int16x8_t  min_base_y128 = vdupq_n_s16(min_base_y);
    int16x8_t  dy128         = vdupq_n_s16(dy);
    uint16x8_t c1234         = vcombine_u16(vcreate_u16(0x0004000300020001), vcreate_u16(0x0008000700060005));

    // Use ext rather than loading left + 30 directly to avoid over-read.
    const uint8x16_t   left_m2   = vld1q_u8(left - 2);
    const uint8x16_t   left_0    = vld1q_u8(left + 0);
    const uint8x16_t   left_16   = vld1q_u8(left + 16);
    const uint8x16_t   left_14   = vextq_u8(left_0, left_16, 14);
    const uint8x16_t   left_30   = vextq_u8(left_16, left_16, 14);
    const uint8x16x3_t left_vals = {{left_m2, left_14, left_30}};

    for (int r = 0; r < N; r++) {
        uint8x8_t    resx, resy, resxy;
        uint16x8x2_t res, shift;
        shift.val[1] = vdupq_n_u16(0);

        int y          = r + 1;
        int base_x     = (-y * dx) >> frac_bits_x;
        int base_shift = 0;
        if (base_x < (min_base_x - 1)) {
            base_shift = (min_base_x - base_x - 1) >> upsample_above;
        }
        int base_min_diff = (min_base_x - base_x + upsample_above) >> upsample_above;
        if (base_min_diff > 8) {
            base_min_diff = 8;
        } else {
            if (base_min_diff < 0)
                base_min_diff = 0;
        }

        uint8x8_t a0_x0, a1_x0;
        if (base_shift > 7) {
            a0_x0        = vdup_n_u8(0);
            a1_x0        = vdup_n_u8(0);
            shift.val[0] = vreinterpretq_u16_u8(v_zero);
            shift.val[1] = vreinterpretq_u16_u8(v_zero);
        } else {
            uint16x8_t ydx = vdupq_n_u16(y * dx);
            uint16x8_t r6  = vshlq_n_u16(vextq_u16(c1234, vreinterpretq_u16_u8(v_zero), 2), 6);

            if (upsample_above) {
                uint8x8x2_t v_tmp;
                v_tmp.val[0]           = vld1_u8(above + base_x + base_shift);
                v_tmp.val[1]           = vld1_u8(above + base_x + base_shift + 8);
                uint8x8_t v_index_low  = vld1_u8(EvenOddMaskx[base_shift]);
                uint8x8_t v_index_high = vld1_u8(EvenOddMaskx[base_shift] + 8);
                shift.val[0] = vshrq_n_u16(vandq_u16(vshlq_u16(vsubq_u16(r6, ydx), v_upsample_above), c3f), 1);
                a0_x0        = vtbl2_u8(v_tmp, v_index_low);
                a1_x0        = vtbl2_u8(v_tmp, v_index_high);
            } else {
                uint8x16_t a0_x128, a1_x128;
                a0_x128 = vld1q_u8(above + base_x + base_shift);
                a1_x128 = vextq_u8(a0_x128, v_zero, 1);
                vector_shuffle(&a0_x128, &v_zero, base_shift);
                vector_shuffle(&a1_x128, &v_zero, base_shift);
                shift.val[0] = vshrq_n_u16(vandq_u16(vsubq_u16(r6, ydx), c3f), 1);
                a0_x0        = vget_low_u8(a0_x128);
                a1_x0        = vget_low_u8(a1_x128);
            }
        }

        diff.val[0] = vsubl_u8(a1_x0, a0_x0); // a[x+1] - a[x]
        a32.val[0]  = vmlal_u8(a16, a0_x0, vdup_n_u8(32)); // a[x] * 32 + 16
        res.val[0]  = vmlaq_u16(a32.val[0], diff.val[0], shift.val[0]);
        resx        = vshrn_n_u16(res.val[0], 5);

        // y calc
        if (base_x < min_base_x) {
            int16x8_t  y_c128, base_y_c128;
            uint16x8_t mask128;
            int16x8_t  v_r6 = vdupq_n_s16(r << 6);

            y_c128      = vmlsq_s16(v_r6, vreinterpretq_s16_u16(c1234), dy128);
            base_y_c128 = vshlq_s16(y_c128, v_frac_bits_y);
            mask128     = vcgtq_s16(min_base_y128, base_y_c128);

            // Values in base_y_c128 range from -2 through 31 inclusive.
            base_y_c128 = vbicq_s16(base_y_c128, vreinterpretq_s16_u16(mask128));

            uint8x16_t left_idx0  = vreinterpretq_u8_s16(base_y_c128 + 2); // [0, 33]
            uint8x16_t left_idx1  = vreinterpretq_u8_s16(base_y_c128 + 3); // [1, 34]
            uint8x16_t left_idx01 = vuzp1q_u8(left_idx0, left_idx1);

            uint8x16_t a01_x = vqtbl3q_u8(left_vals, left_idx01);
            uint8x8_t  a0_x1 = vget_low_u8(a01_x);
            uint8x8_t  a1_x1 = vget_high_u8(a01_x);

            if (upsample_left) {
                shift.val[1] = vshrq_n_u16(vandq_u16(vshlq_u16(vreinterpretq_u16_s16(y_c128), v_upsample_left), c3f),
                                           1);
            } else {
                shift.val[1] = vshrq_n_u16(vandq_u16(vreinterpretq_u16_s16(y_c128), c3f), 1);
            }

            diff.val[1]    = vsubl_u8(a1_x1, a0_x1);
            a32.val[1]     = vmlal_u8(a16, a0_x1, vdup_n_u8(32));
            res.val[1]     = vmlaq_u16(a32.val[1], diff.val[1], shift.val[1]);
            resy           = vshrn_n_u16(res.val[1], 5);
            uint8x8_t mask = vld1_u8(BaseMask[base_min_diff]);
            resxy          = vbsl_u8(mask, resy, resx);
            vst1_u8(dst, resxy);
        } else {
            vst1_u8(dst, resx);
        }

        dst += stride;
    }
}

static INLINE void dr_prediction_z2_HxW_neon(int H, int W, uint8_t *dst, ptrdiff_t stride, const uint8_t *above,
                                             const uint8_t *left, int dx, int dy) {
    // here upsample_above and upsample_left are 0 by design of
    // av1_use_intra_edge_upsample
    const int min_base_x  = -1;
    const int min_base_y  = -1;
    const int frac_bits_x = 6;
    const int frac_bits_y = 6;

    uint16x8x2_t a32, c0123, c1234, diff, shifty;
    uint8x16x2_t a0_x, a1_x;
    uint16x8_t   v_32          = vdupq_n_u16(32);
    uint8x16_t   v_zero        = vdupq_n_u8(0);
    int16x8_t    v_frac_bits_y = vdupq_n_s16(-frac_bits_y);

    uint16x8_t a16           = vdupq_n_u16(16);
    uint16x8_t c1            = vshrq_n_u16(a16, 4);
    int16x8_t  min_base_y256 = vdupq_n_s16(min_base_y);
    uint16x8_t c3f           = vdupq_n_u16(0x3f);
    int16x8_t  dy256         = vdupq_n_s16(dy);
    c0123.val[0]             = vcombine_u16(vcreate_u16(0x0003000200010000), vcreate_u16(0x0007000600050004));
    c0123.val[1]             = vcombine_u16(vcreate_u16(0x000B000A00090008), vcreate_u16(0x000F000E000D000C));
    c1234.val[0]             = vaddq_u16(c0123.val[0], c1);
    c1234.val[1]             = vaddq_u16(c0123.val[1], c1);

    const uint8x16_t   left_m1    = vld1q_u8(left - 1);
    const uint8x16_t   left_0     = vld1q_u8(left + 0);
    const uint8x16_t   left_16    = vld1q_u8(left + 16);
    const uint8x16_t   left_32    = vld1q_u8(left + 32);
    const uint8x16_t   left_48    = vld1q_u8(left + 48);
    const uint8x16_t   left_15    = vextq_u8(left_0, left_16, 15);
    const uint8x16_t   left_31    = vextq_u8(left_16, left_32, 15);
    const uint8x16_t   left_47    = vextq_u8(left_32, left_48, 15);
    const uint8x16x4_t left_vals0 = {{left_m1, left_15, left_31, left_47}};
    const uint8x16x4_t left_vals1 = {{left_0, left_16, left_32, left_48}};

    for (int r = 0; r < H; r++) {
        uint16x8x2_t res, r6, shift;
        uint16x8_t   j256;
        uint8x16_t   resx, resy, resxy;
        int          y   = r + 1;
        uint16x8_t   ydx = vdupq_n_u16((uint16_t)(y * dx));

        int base_x = (-y * dx) >> frac_bits_x;
        for (int j = 0; j < W; j += 16) {
            j256 = vdupq_n_u16(j);

            int base_shift = 0;
            if ((base_x + j) < (min_base_x - 1)) {
                base_shift = (min_base_x - (base_x + j) - 1);
            }
            int base_min_diff = (min_base_x - base_x - j);
            if (base_min_diff > 16) {
                base_min_diff = 16;
            } else {
                if (base_min_diff < 0)
                    base_min_diff = 0;
            }

            if (base_shift < 16) {
                uint8x16_t a0_x128, a1_x128;
                a0_x128 = vld1q_u8(above + base_x + base_shift + j);
                a1_x128 = vld1q_u8(above + base_x + base_shift + 1 + j);
                vector_shuffle(&a0_x128, &v_zero, base_shift);
                vector_shuffle(&a1_x128, &v_zero, base_shift);
                a0_x         = vzipq_u8(a0_x128, v_zero);
                a1_x         = vzipq_u8(a1_x128, v_zero);
                r6.val[0]    = vshlq_n_u16(vaddq_u16(c0123.val[0], j256), 6);
                r6.val[1]    = vshlq_n_u16(vaddq_u16(c0123.val[1], j256), 6);
                shift.val[0] = vshrq_n_u16(vandq_u16(vsubq_u16(r6.val[0], ydx), c3f), 1);
                shift.val[1] = vshrq_n_u16(vandq_u16(vsubq_u16(r6.val[1], ydx), c3f), 1);
                diff.val[0]  = vsubq_u16(vreinterpretq_u16_u8(a1_x.val[0]),
                                        vreinterpretq_u16_u8(a0_x.val[0])); // a[x+1] - a[x]
                diff.val[1]  = vsubq_u16(vreinterpretq_u16_u8(a1_x.val[1]),
                                        vreinterpretq_u16_u8(a0_x.val[1])); // a[x+1] - a[x]
                a32.val[0]   = vmlaq_u16(a16, vreinterpretq_u16_u8(a0_x.val[0]),
                                       v_32); // a[x] * 32 + 16
                a32.val[1]   = vmlaq_u16(a16, vreinterpretq_u16_u8(a0_x.val[1]),
                                       v_32); // a[x] * 32 + 16
                res.val[0]   = vmlaq_u16(a32.val[0], diff.val[0], shift.val[0]);
                res.val[1]   = vmlaq_u16(a32.val[1], diff.val[1], shift.val[1]);
                resx         = vcombine_u8(vshrn_n_u16(res.val[0], 5), vshrn_n_u16(res.val[1], 5));
            } else {
                resx = v_zero;
            }

            // y calc
            if (base_x < min_base_x) {
                uint16x8x2_t mask256;
                int16x8x2_t  c256, y_c256, base_y_c256, mul16;
                int16x8_t    v_r6 = vdupq_n_s16(r << 6);

                c256.val[0]   = vaddq_s16(vreinterpretq_s16_u16(j256), vreinterpretq_s16_u16(c1234.val[0]));
                c256.val[1]   = vaddq_s16(vreinterpretq_s16_u16(j256), vreinterpretq_s16_u16(c1234.val[1]));
                mul16.val[0]  = vreinterpretq_s16_u16(vminq_u16(vreinterpretq_u16_s16(vmulq_s16(c256.val[0], dy256)),
                                                               vshrq_n_u16(vreinterpretq_u16_s16(min_base_y256), 1)));
                mul16.val[1]  = vreinterpretq_s16_u16(vminq_u16(vreinterpretq_u16_s16(vmulq_s16(c256.val[1], dy256)),
                                                               vshrq_n_u16(vreinterpretq_u16_s16(min_base_y256), 1)));
                y_c256.val[0] = vsubq_s16(v_r6, mul16.val[0]);
                y_c256.val[1] = vsubq_s16(v_r6, mul16.val[1]);

                base_y_c256.val[0] = vshlq_s16(y_c256.val[0], v_frac_bits_y);
                base_y_c256.val[1] = vshlq_s16(y_c256.val[1], v_frac_bits_y);
                mask256.val[0]     = vcgtq_s16(min_base_y256, base_y_c256.val[0]);
                mask256.val[1]     = vcgtq_s16(min_base_y256, base_y_c256.val[1]);

                base_y_c256.val[0] = vbslq_s16(mask256.val[0], min_base_y256, base_y_c256.val[0]);
                base_y_c256.val[1] = vbslq_s16(mask256.val[1], min_base_y256, base_y_c256.val[1]);

                int16_t min_y       = vgetq_lane_s16(base_y_c256.val[1], 7);
                int16_t max_y       = vgetq_lane_s16(base_y_c256.val[0], 0);
                int16_t offset_diff = max_y - min_y;

                uint8x8_t a0_y0;
                uint8x8_t a0_y1;
                uint8x8_t a1_y0;
                uint8x8_t a1_y1;

                if (offset_diff < 16) {
                    assert(offset_diff >= 0);
                    int16x8_t min_y256 = vdupq_lane_s16(vget_high_s16(base_y_c256.val[1]), 3);

                    int16x8x2_t base_y_offset;
                    base_y_offset.val[0] = vsubq_s16(base_y_c256.val[0], min_y256);
                    base_y_offset.val[1] = vsubq_s16(base_y_c256.val[1], min_y256);

                    int8x16_t base_y_offset128 = vcombine_s8(vqmovn_s16(base_y_offset.val[0]),
                                                             vqmovn_s16(base_y_offset.val[1]));

                    uint8x16_t a0_y128, a1_y128;
                    uint8x16_t v_loadmaskz2 = vld1q_u8(LoadMaskz2[offset_diff / 4]);
                    a0_y128                 = vld1q_u8(left + min_y);
                    a0_y128                 = vandq_u8(a0_y128, v_loadmaskz2);
                    a1_y128                 = vld1q_u8(left + min_y + 1);
                    a1_y128                 = vandq_u8(a1_y128, v_loadmaskz2);
                    a0_y128                 = vqtbl1q_u8(a0_y128, vreinterpretq_u8_s8(base_y_offset128));
                    a1_y128                 = vqtbl1q_u8(a1_y128, vreinterpretq_u8_s8(base_y_offset128));

                    a0_y0 = vget_low_u8(a0_y128);
                    a0_y1 = vget_high_u8(a0_y128);
                    a1_y0 = vget_low_u8(a1_y128);
                    a1_y1 = vget_high_u8(a1_y128);
                } else {
                    // Values in base_y_c256 range from -1 through 62 inclusive.
                    base_y_c256.val[0] = vbicq_s16(base_y_c256.val[0], vreinterpretq_s16_u16(mask256.val[0]));
                    base_y_c256.val[1] = vbicq_s16(base_y_c256.val[1], vreinterpretq_s16_u16(mask256.val[1]));

                    // Values in left_idx{0,1} range from 0 through 63 inclusive.
                    uint8x16_t left_idx0 = vreinterpretq_u8_s16(base_y_c256.val[0] + 1);
                    uint8x16_t left_idx1 = vreinterpretq_u8_s16(base_y_c256.val[1] + 1);

                    uint8x16_t left_idx01 = vuzp1q_u8(left_idx0, left_idx1);

                    uint8x16_t a0_y01 = vqtbl4q_u8(left_vals0, left_idx01);
                    uint8x16_t a1_y01 = vqtbl4q_u8(left_vals1, left_idx01);

                    a0_y0 = vget_low_u8(a0_y01);
                    a0_y1 = vget_high_u8(a0_y01);
                    a1_y0 = vget_low_u8(a1_y01);
                    a1_y1 = vget_high_u8(a1_y01);
                }

                shifty.val[0] = vshrq_n_u16(vandq_u16(vreinterpretq_u16_s16(y_c256.val[0]), c3f), 1);
                shifty.val[1] = vshrq_n_u16(vandq_u16(vreinterpretq_u16_s16(y_c256.val[1]), c3f), 1);
                diff.val[0]   = vsubl_u8(a1_y0, a0_y0); // a[x+1] - a[x]
                diff.val[1]   = vsubl_u8(a1_y1, a0_y1); // a[x+1] - a[x]
                a32.val[0]    = vmlal_u8(a16, a0_y0, vdup_n_u8(32)); // a[x] * 32 + 16
                a32.val[1]    = vmlal_u8(a16, a0_y1, vdup_n_u8(32)); // a[x] * 32 + 16
                res.val[0]    = vmlaq_u16(a32.val[0], diff.val[0], shifty.val[0]);
                res.val[1]    = vmlaq_u16(a32.val[1], diff.val[1], shifty.val[1]);

                resy = vcombine_u8(vshrn_n_u16(res.val[0], 5), vshrn_n_u16(res.val[1], 5));
            } else {
                resy = v_zero;
            }
            uint8x16_t mask = vld1q_u8(BaseMask[base_min_diff]);
            resxy           = vbslq_u8(mask, resy, resx);
            vst1q_u8(dst + j, resxy);
        } // for j
        dst += stride;
    }
}

// Directional prediction, zone 2: 90 < angle < 180
void svt_av1_dr_prediction_z2_neon(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint8_t *above,
                                   const uint8_t *left, int32_t upsample_above, int32_t upsample_left, int32_t dx,
                                   int32_t dy) {
    assert(dx > 0);
    assert(dy > 0);

    switch (bw) {
    case 4: {
        dr_prediction_z2_Nx4_neon(bh, dst, stride, above, left, upsample_above, upsample_left, dx, dy);
        break;
    }
    case 8: {
        dr_prediction_z2_Nx8_neon(bh, dst, stride, above, left, upsample_above, upsample_left, dx, dy);
        break;
    }
    default: {
        dr_prediction_z2_HxW_neon(bh, bw, dst, stride, above, left, dx, dy);
        break;
    }
    }
}

/* ---------------------P R E D I C T I O N   Z 3--------------------------- */

static INLINE void transpose4x16_neon(uint8x16_t *x, uint16x8x2_t *d) {
    uint8x16x2_t w0, w1;

    w0 = vzipq_u8(x[0], x[1]);
    w1 = vzipq_u8(x[2], x[3]);

    d[0] = vzipq_u16(vreinterpretq_u16_u8(w0.val[0]), vreinterpretq_u16_u8(w1.val[0]));
    d[1] = vzipq_u16(vreinterpretq_u16_u8(w0.val[1]), vreinterpretq_u16_u8(w1.val[1]));
}

static INLINE void transpose4x8_8x4_low_neon(uint8x8_t *x, uint16x4x2_t *d) {
    uint8x8x2_t w0, w1;

    w0 = vzip_u8(x[0], x[1]);
    w1 = vzip_u8(x[2], x[3]);

    *d = vzip_u16(vreinterpret_u16_u8(w0.val[0]), vreinterpret_u16_u8(w1.val[0]));
}

static INLINE void transpose4x8_8x4_neon(uint8x8_t *x, uint16x4x2_t *d) {
    uint8x8x2_t w0, w1;

    w0 = vzip_u8(x[0], x[1]);
    w1 = vzip_u8(x[2], x[3]);

    d[0] = vzip_u16(vreinterpret_u16_u8(w0.val[0]), vreinterpret_u16_u8(w1.val[0]));
    d[1] = vzip_u16(vreinterpret_u16_u8(w0.val[1]), vreinterpret_u16_u8(w1.val[1]));
}

static INLINE void transpose8x8_low_neon(uint8x8_t *x, uint32x2x2_t *d) {
    uint8x8x2_t  w0, w1, w2, w3;
    uint16x4x2_t w4, w5;

    w0 = vzip_u8(x[0], x[1]);
    w1 = vzip_u8(x[2], x[3]);
    w2 = vzip_u8(x[4], x[5]);
    w3 = vzip_u8(x[6], x[7]);

    w4 = vzip_u16(vreinterpret_u16_u8(w0.val[0]), vreinterpret_u16_u8(w1.val[0]));
    w5 = vzip_u16(vreinterpret_u16_u8(w2.val[0]), vreinterpret_u16_u8(w3.val[0]));

    d[0] = vzip_u32(vreinterpret_u32_u16(w4.val[0]), vreinterpret_u32_u16(w5.val[0]));
    d[1] = vzip_u32(vreinterpret_u32_u16(w4.val[1]), vreinterpret_u32_u16(w5.val[1]));
}

static INLINE void transpose8x8_neon(uint8x8_t *x, uint32x2x2_t *d) {
    uint8x8x2_t  w0, w1, w2, w3;
    uint16x4x2_t w4, w5, w6, w7;

    w0 = vzip_u8(x[0], x[1]);
    w1 = vzip_u8(x[2], x[3]);
    w2 = vzip_u8(x[4], x[5]);
    w3 = vzip_u8(x[6], x[7]);

    w4 = vzip_u16(vreinterpret_u16_u8(w0.val[0]), vreinterpret_u16_u8(w1.val[0]));
    w5 = vzip_u16(vreinterpret_u16_u8(w2.val[0]), vreinterpret_u16_u8(w3.val[0]));

    d[0] = vzip_u32(vreinterpret_u32_u16(w4.val[0]), vreinterpret_u32_u16(w5.val[0]));
    d[1] = vzip_u32(vreinterpret_u32_u16(w4.val[1]), vreinterpret_u32_u16(w5.val[1]));

    w6 = vzip_u16(vreinterpret_u16_u8(w0.val[1]), vreinterpret_u16_u8(w1.val[1]));
    w7 = vzip_u16(vreinterpret_u16_u8(w2.val[1]), vreinterpret_u16_u8(w3.val[1]));

    d[2] = vzip_u32(vreinterpret_u32_u16(w6.val[0]), vreinterpret_u32_u16(w7.val[0]));
    d[3] = vzip_u32(vreinterpret_u32_u16(w6.val[1]), vreinterpret_u32_u16(w7.val[1]));
}

static INLINE void transpose16x8_8x16_neon(uint8x8_t *x, uint64x2_t *d) {
    uint8x8x2_t  w0, w1, w2, w3, w8, w9, w10, w11;
    uint16x4x2_t w4, w5, w12, w13;
    uint32x2x2_t w6, w7, w14, w15;

    w0 = vzip_u8(x[0], x[1]);
    w1 = vzip_u8(x[2], x[3]);
    w2 = vzip_u8(x[4], x[5]);
    w3 = vzip_u8(x[6], x[7]);

    w8  = vzip_u8(x[8], x[9]);
    w9  = vzip_u8(x[10], x[11]);
    w10 = vzip_u8(x[12], x[13]);
    w11 = vzip_u8(x[14], x[15]);

    w4  = vzip_u16(vreinterpret_u16_u8(w0.val[0]), vreinterpret_u16_u8(w1.val[0]));
    w5  = vzip_u16(vreinterpret_u16_u8(w2.val[0]), vreinterpret_u16_u8(w3.val[0]));
    w12 = vzip_u16(vreinterpret_u16_u8(w8.val[0]), vreinterpret_u16_u8(w9.val[0]));
    w13 = vzip_u16(vreinterpret_u16_u8(w10.val[0]), vreinterpret_u16_u8(w11.val[0]));

    w6  = vzip_u32(vreinterpret_u32_u16(w4.val[0]), vreinterpret_u32_u16(w5.val[0]));
    w7  = vzip_u32(vreinterpret_u32_u16(w4.val[1]), vreinterpret_u32_u16(w5.val[1]));
    w14 = vzip_u32(vreinterpret_u32_u16(w12.val[0]), vreinterpret_u32_u16(w13.val[0]));
    w15 = vzip_u32(vreinterpret_u32_u16(w12.val[1]), vreinterpret_u32_u16(w13.val[1]));

    // Store first 4-line result
    d[0] = vcombine_u64(vreinterpret_u64_u32(w6.val[0]), vreinterpret_u64_u32(w14.val[0]));
    d[1] = vcombine_u64(vreinterpret_u64_u32(w6.val[1]), vreinterpret_u64_u32(w14.val[1]));
    d[2] = vcombine_u64(vreinterpret_u64_u32(w7.val[0]), vreinterpret_u64_u32(w15.val[0]));
    d[3] = vcombine_u64(vreinterpret_u64_u32(w7.val[1]), vreinterpret_u64_u32(w15.val[1]));

    w4  = vzip_u16(vreinterpret_u16_u8(w0.val[1]), vreinterpret_u16_u8(w1.val[1]));
    w5  = vzip_u16(vreinterpret_u16_u8(w2.val[1]), vreinterpret_u16_u8(w3.val[1]));
    w12 = vzip_u16(vreinterpret_u16_u8(w8.val[1]), vreinterpret_u16_u8(w9.val[1]));
    w13 = vzip_u16(vreinterpret_u16_u8(w10.val[1]), vreinterpret_u16_u8(w11.val[1]));

    w6  = vzip_u32(vreinterpret_u32_u16(w4.val[0]), vreinterpret_u32_u16(w5.val[0]));
    w7  = vzip_u32(vreinterpret_u32_u16(w4.val[1]), vreinterpret_u32_u16(w5.val[1]));
    w14 = vzip_u32(vreinterpret_u32_u16(w12.val[0]), vreinterpret_u32_u16(w13.val[0]));
    w15 = vzip_u32(vreinterpret_u32_u16(w12.val[1]), vreinterpret_u32_u16(w13.val[1]));

    // Store second 4-line result
    d[4] = vcombine_u64(vreinterpret_u64_u32(w6.val[0]), vreinterpret_u64_u32(w14.val[0]));
    d[5] = vcombine_u64(vreinterpret_u64_u32(w6.val[1]), vreinterpret_u64_u32(w14.val[1]));
    d[6] = vcombine_u64(vreinterpret_u64_u32(w7.val[0]), vreinterpret_u64_u32(w15.val[0]));
    d[7] = vcombine_u64(vreinterpret_u64_u32(w7.val[1]), vreinterpret_u64_u32(w15.val[1]));
}

static INLINE void transpose8x16_16x8_neon(uint8x16_t *x, uint64x2_t *d) {
    uint8x16x2_t w0, w1, w2, w3;
    uint16x8x2_t w4, w5, w6, w7;
    uint32x4x2_t w8, w9, w10, w11;

    w0 = vzipq_u8(x[0], x[1]);
    w1 = vzipq_u8(x[2], x[3]);
    w2 = vzipq_u8(x[4], x[5]);
    w3 = vzipq_u8(x[6], x[7]);

    w4 = vzipq_u16(vreinterpretq_u16_u8(w0.val[0]), vreinterpretq_u16_u8(w1.val[0]));
    w5 = vzipq_u16(vreinterpretq_u16_u8(w2.val[0]), vreinterpretq_u16_u8(w3.val[0]));
    w6 = vzipq_u16(vreinterpretq_u16_u8(w0.val[1]), vreinterpretq_u16_u8(w1.val[1]));
    w7 = vzipq_u16(vreinterpretq_u16_u8(w2.val[1]), vreinterpretq_u16_u8(w3.val[1]));

    w8  = vzipq_u32(vreinterpretq_u32_u16(w4.val[0]), vreinterpretq_u32_u16(w5.val[0]));
    w9  = vzipq_u32(vreinterpretq_u32_u16(w6.val[0]), vreinterpretq_u32_u16(w7.val[0]));
    w10 = vzipq_u32(vreinterpretq_u32_u16(w4.val[1]), vreinterpretq_u32_u16(w5.val[1]));
    w11 = vzipq_u32(vreinterpretq_u32_u16(w6.val[1]), vreinterpretq_u32_u16(w7.val[1]));

    d[0] = vzip1q_u64(vreinterpretq_u64_u32(w8.val[0]), vreinterpretq_u64_u32(w9.val[0]));
    d[1] = vzip2q_u64(vreinterpretq_u64_u32(w8.val[0]), vreinterpretq_u64_u32(w9.val[0]));
    d[2] = vzip1q_u64(vreinterpretq_u64_u32(w8.val[1]), vreinterpretq_u64_u32(w9.val[1]));
    d[3] = vzip2q_u64(vreinterpretq_u64_u32(w8.val[1]), vreinterpretq_u64_u32(w9.val[1]));
    d[4] = vzip1q_u64(vreinterpretq_u64_u32(w10.val[0]), vreinterpretq_u64_u32(w11.val[0]));
    d[5] = vzip2q_u64(vreinterpretq_u64_u32(w10.val[0]), vreinterpretq_u64_u32(w11.val[0]));
    d[6] = vzip1q_u64(vreinterpretq_u64_u32(w10.val[1]), vreinterpretq_u64_u32(w11.val[1]));
    d[7] = vzip2q_u64(vreinterpretq_u64_u32(w10.val[1]), vreinterpretq_u64_u32(w11.val[1]));
}

static INLINE void transpose16x16_neon(uint8x16_t *x, uint64x2_t *d) {
    uint8x16x2_t w0, w1, w2, w3, w4, w5, w6, w7;
    uint16x8x2_t w8, w9, w10, w11;
    uint32x4x2_t w12, w13, w14, w15;

    w0 = vzipq_u8(x[0], x[1]);
    w1 = vzipq_u8(x[2], x[3]);
    w2 = vzipq_u8(x[4], x[5]);
    w3 = vzipq_u8(x[6], x[7]);

    w4 = vzipq_u8(x[8], x[9]);
    w5 = vzipq_u8(x[10], x[11]);
    w6 = vzipq_u8(x[12], x[13]);
    w7 = vzipq_u8(x[14], x[15]);

    w8  = vzipq_u16(vreinterpretq_u16_u8(w0.val[0]), vreinterpretq_u16_u8(w1.val[0]));
    w9  = vzipq_u16(vreinterpretq_u16_u8(w2.val[0]), vreinterpretq_u16_u8(w3.val[0]));
    w10 = vzipq_u16(vreinterpretq_u16_u8(w4.val[0]), vreinterpretq_u16_u8(w5.val[0]));
    w11 = vzipq_u16(vreinterpretq_u16_u8(w6.val[0]), vreinterpretq_u16_u8(w7.val[0]));

    w12 = vzipq_u32(vreinterpretq_u32_u16(w8.val[0]), vreinterpretq_u32_u16(w9.val[0]));
    w13 = vzipq_u32(vreinterpretq_u32_u16(w10.val[0]), vreinterpretq_u32_u16(w11.val[0]));
    w14 = vzipq_u32(vreinterpretq_u32_u16(w8.val[1]), vreinterpretq_u32_u16(w9.val[1]));
    w15 = vzipq_u32(vreinterpretq_u32_u16(w10.val[1]), vreinterpretq_u32_u16(w11.val[1]));

    d[0] = vzip1q_u64(vreinterpretq_u64_u32(w12.val[0]), vreinterpretq_u64_u32(w13.val[0]));
    d[1] = vzip2q_u64(vreinterpretq_u64_u32(w12.val[0]), vreinterpretq_u64_u32(w13.val[0]));
    d[2] = vzip1q_u64(vreinterpretq_u64_u32(w12.val[1]), vreinterpretq_u64_u32(w13.val[1]));
    d[3] = vzip2q_u64(vreinterpretq_u64_u32(w12.val[1]), vreinterpretq_u64_u32(w13.val[1]));
    d[4] = vzip1q_u64(vreinterpretq_u64_u32(w14.val[0]), vreinterpretq_u64_u32(w15.val[0]));
    d[5] = vzip2q_u64(vreinterpretq_u64_u32(w14.val[0]), vreinterpretq_u64_u32(w15.val[0]));
    d[6] = vzip1q_u64(vreinterpretq_u64_u32(w14.val[1]), vreinterpretq_u64_u32(w15.val[1]));
    d[7] = vzip2q_u64(vreinterpretq_u64_u32(w14.val[1]), vreinterpretq_u64_u32(w15.val[1]));

    // upper half
    w8  = vzipq_u16(vreinterpretq_u16_u8(w0.val[1]), vreinterpretq_u16_u8(w1.val[1]));
    w9  = vzipq_u16(vreinterpretq_u16_u8(w2.val[1]), vreinterpretq_u16_u8(w3.val[1]));
    w10 = vzipq_u16(vreinterpretq_u16_u8(w4.val[1]), vreinterpretq_u16_u8(w5.val[1]));
    w11 = vzipq_u16(vreinterpretq_u16_u8(w6.val[1]), vreinterpretq_u16_u8(w7.val[1]));

    w12 = vzipq_u32(vreinterpretq_u32_u16(w8.val[0]), vreinterpretq_u32_u16(w9.val[0]));
    w13 = vzipq_u32(vreinterpretq_u32_u16(w10.val[0]), vreinterpretq_u32_u16(w11.val[0]));
    w14 = vzipq_u32(vreinterpretq_u32_u16(w8.val[1]), vreinterpretq_u32_u16(w9.val[1]));
    w15 = vzipq_u32(vreinterpretq_u32_u16(w10.val[1]), vreinterpretq_u32_u16(w11.val[1]));

    d[8]  = vzip1q_u64(vreinterpretq_u64_u32(w12.val[0]), vreinterpretq_u64_u32(w13.val[0]));
    d[9]  = vzip2q_u64(vreinterpretq_u64_u32(w12.val[0]), vreinterpretq_u64_u32(w13.val[0]));
    d[10] = vzip1q_u64(vreinterpretq_u64_u32(w12.val[1]), vreinterpretq_u64_u32(w13.val[1]));
    d[11] = vzip2q_u64(vreinterpretq_u64_u32(w12.val[1]), vreinterpretq_u64_u32(w13.val[1]));
    d[12] = vzip1q_u64(vreinterpretq_u64_u32(w14.val[0]), vreinterpretq_u64_u32(w15.val[0]));
    d[13] = vzip2q_u64(vreinterpretq_u64_u32(w14.val[0]), vreinterpretq_u64_u32(w15.val[0]));
    d[14] = vzip1q_u64(vreinterpretq_u64_u32(w14.val[1]), vreinterpretq_u64_u32(w15.val[1]));
    d[15] = vzip2q_u64(vreinterpretq_u64_u32(w14.val[1]), vreinterpretq_u64_u32(w15.val[1]));
}

static INLINE void transpose16x32_neon(uint8x16x2_t *x, uint64x2x2_t *d) {
    uint8x16x2_t w0, w1, w2, w3, w8, w9, w10, w11;
    uint16x8x2_t w4, w5, w12, w13;
    uint32x4x2_t w6, w7, w14, w15;

    w0 = vzipq_u8(x[0].val[0], x[1].val[0]);
    w1 = vzipq_u8(x[2].val[0], x[3].val[0]);
    w2 = vzipq_u8(x[4].val[0], x[5].val[0]);
    w3 = vzipq_u8(x[6].val[0], x[7].val[0]);

    w8  = vzipq_u8(x[8].val[0], x[9].val[0]);
    w9  = vzipq_u8(x[10].val[0], x[11].val[0]);
    w10 = vzipq_u8(x[12].val[0], x[13].val[0]);
    w11 = vzipq_u8(x[14].val[0], x[15].val[0]);

    w4  = vzipq_u16(vreinterpretq_u16_u8(w0.val[0]), vreinterpretq_u16_u8(w1.val[0]));
    w5  = vzipq_u16(vreinterpretq_u16_u8(w2.val[0]), vreinterpretq_u16_u8(w3.val[0]));
    w12 = vzipq_u16(vreinterpretq_u16_u8(w8.val[0]), vreinterpretq_u16_u8(w9.val[0]));
    w13 = vzipq_u16(vreinterpretq_u16_u8(w10.val[0]), vreinterpretq_u16_u8(w11.val[0]));

    w6  = vzipq_u32(vreinterpretq_u32_u16(w4.val[0]), vreinterpretq_u32_u16(w5.val[0]));
    w7  = vzipq_u32(vreinterpretq_u32_u16(w4.val[1]), vreinterpretq_u32_u16(w5.val[1]));
    w14 = vzipq_u32(vreinterpretq_u32_u16(w12.val[0]), vreinterpretq_u32_u16(w13.val[0]));
    w15 = vzipq_u32(vreinterpretq_u32_u16(w12.val[1]), vreinterpretq_u32_u16(w13.val[1]));

    // Store first 4-line result

    d[0].val[0] = vzip1q_u64(vreinterpretq_u64_u32(w6.val[0]), vreinterpretq_u64_u32(w14.val[0]));
    d[0].val[1] = vzip2q_u64(vreinterpretq_u64_u32(w6.val[0]), vreinterpretq_u64_u32(w14.val[0]));
    d[1].val[0] = vzip1q_u64(vreinterpretq_u64_u32(w6.val[1]), vreinterpretq_u64_u32(w14.val[1]));
    d[1].val[1] = vzip2q_u64(vreinterpretq_u64_u32(w6.val[1]), vreinterpretq_u64_u32(w14.val[1]));
    d[2].val[0] = vzip1q_u64(vreinterpretq_u64_u32(w7.val[0]), vreinterpretq_u64_u32(w15.val[0]));
    d[2].val[1] = vzip2q_u64(vreinterpretq_u64_u32(w7.val[0]), vreinterpretq_u64_u32(w15.val[0]));
    d[3].val[0] = vzip1q_u64(vreinterpretq_u64_u32(w7.val[1]), vreinterpretq_u64_u32(w15.val[1]));
    d[3].val[1] = vzip2q_u64(vreinterpretq_u64_u32(w7.val[1]), vreinterpretq_u64_u32(w15.val[1]));

    w4  = vzipq_u16(vreinterpretq_u16_u8(w0.val[1]), vreinterpretq_u16_u8(w1.val[1]));
    w5  = vzipq_u16(vreinterpretq_u16_u8(w2.val[1]), vreinterpretq_u16_u8(w3.val[1]));
    w12 = vzipq_u16(vreinterpretq_u16_u8(w8.val[1]), vreinterpretq_u16_u8(w9.val[1]));
    w13 = vzipq_u16(vreinterpretq_u16_u8(w10.val[1]), vreinterpretq_u16_u8(w11.val[1]));

    w6  = vzipq_u32(vreinterpretq_u32_u16(w4.val[0]), vreinterpretq_u32_u16(w5.val[0]));
    w7  = vzipq_u32(vreinterpretq_u32_u16(w4.val[1]), vreinterpretq_u32_u16(w5.val[1]));
    w14 = vzipq_u32(vreinterpretq_u32_u16(w12.val[0]), vreinterpretq_u32_u16(w13.val[0]));
    w15 = vzipq_u32(vreinterpretq_u32_u16(w12.val[1]), vreinterpretq_u32_u16(w13.val[1]));

    // Store second 4-line result

    d[4].val[0] = vzip1q_u64(vreinterpretq_u64_u32(w6.val[0]), vreinterpretq_u64_u32(w14.val[0]));
    d[4].val[1] = vzip2q_u64(vreinterpretq_u64_u32(w6.val[0]), vreinterpretq_u64_u32(w14.val[0]));
    d[5].val[0] = vzip1q_u64(vreinterpretq_u64_u32(w6.val[1]), vreinterpretq_u64_u32(w14.val[1]));
    d[5].val[1] = vzip2q_u64(vreinterpretq_u64_u32(w6.val[1]), vreinterpretq_u64_u32(w14.val[1]));
    d[6].val[0] = vzip1q_u64(vreinterpretq_u64_u32(w7.val[0]), vreinterpretq_u64_u32(w15.val[0]));
    d[6].val[1] = vzip2q_u64(vreinterpretq_u64_u32(w7.val[0]), vreinterpretq_u64_u32(w15.val[0]));
    d[7].val[0] = vzip1q_u64(vreinterpretq_u64_u32(w7.val[1]), vreinterpretq_u64_u32(w15.val[1]));
    d[7].val[1] = vzip2q_u64(vreinterpretq_u64_u32(w7.val[1]), vreinterpretq_u64_u32(w15.val[1]));

    // upper half
    w0 = vzipq_u8(x[0].val[1], x[1].val[1]);
    w1 = vzipq_u8(x[2].val[1], x[3].val[1]);
    w2 = vzipq_u8(x[4].val[1], x[5].val[1]);
    w3 = vzipq_u8(x[6].val[1], x[7].val[1]);

    w8  = vzipq_u8(x[8].val[1], x[9].val[1]);
    w9  = vzipq_u8(x[10].val[1], x[11].val[1]);
    w10 = vzipq_u8(x[12].val[1], x[13].val[1]);
    w11 = vzipq_u8(x[14].val[1], x[15].val[1]);

    w4  = vzipq_u16(vreinterpretq_u16_u8(w0.val[0]), vreinterpretq_u16_u8(w1.val[0]));
    w5  = vzipq_u16(vreinterpretq_u16_u8(w2.val[0]), vreinterpretq_u16_u8(w3.val[0]));
    w12 = vzipq_u16(vreinterpretq_u16_u8(w8.val[0]), vreinterpretq_u16_u8(w9.val[0]));
    w13 = vzipq_u16(vreinterpretq_u16_u8(w10.val[0]), vreinterpretq_u16_u8(w11.val[0]));

    w6  = vzipq_u32(vreinterpretq_u32_u16(w4.val[0]), vreinterpretq_u32_u16(w5.val[0]));
    w7  = vzipq_u32(vreinterpretq_u32_u16(w4.val[1]), vreinterpretq_u32_u16(w5.val[1]));
    w14 = vzipq_u32(vreinterpretq_u32_u16(w12.val[0]), vreinterpretq_u32_u16(w13.val[0]));
    w15 = vzipq_u32(vreinterpretq_u32_u16(w12.val[1]), vreinterpretq_u32_u16(w13.val[1]));

    // Store first 4-line result

    d[8].val[0]  = vzip1q_u64(vreinterpretq_u64_u32(w6.val[0]), vreinterpretq_u64_u32(w14.val[0]));
    d[8].val[1]  = vzip2q_u64(vreinterpretq_u64_u32(w6.val[0]), vreinterpretq_u64_u32(w14.val[0]));
    d[9].val[0]  = vzip1q_u64(vreinterpretq_u64_u32(w6.val[1]), vreinterpretq_u64_u32(w14.val[1]));
    d[9].val[1]  = vzip2q_u64(vreinterpretq_u64_u32(w6.val[1]), vreinterpretq_u64_u32(w14.val[1]));
    d[10].val[0] = vzip1q_u64(vreinterpretq_u64_u32(w7.val[0]), vreinterpretq_u64_u32(w15.val[0]));
    d[10].val[1] = vzip2q_u64(vreinterpretq_u64_u32(w7.val[0]), vreinterpretq_u64_u32(w15.val[0]));
    d[11].val[0] = vzip1q_u64(vreinterpretq_u64_u32(w7.val[1]), vreinterpretq_u64_u32(w15.val[1]));
    d[11].val[1] = vzip2q_u64(vreinterpretq_u64_u32(w7.val[1]), vreinterpretq_u64_u32(w15.val[1]));

    w4  = vzipq_u16(vreinterpretq_u16_u8(w0.val[1]), vreinterpretq_u16_u8(w1.val[1]));
    w5  = vzipq_u16(vreinterpretq_u16_u8(w2.val[1]), vreinterpretq_u16_u8(w3.val[1]));
    w12 = vzipq_u16(vreinterpretq_u16_u8(w8.val[1]), vreinterpretq_u16_u8(w9.val[1]));
    w13 = vzipq_u16(vreinterpretq_u16_u8(w10.val[1]), vreinterpretq_u16_u8(w11.val[1]));

    w6  = vzipq_u32(vreinterpretq_u32_u16(w4.val[0]), vreinterpretq_u32_u16(w5.val[0]));
    w7  = vzipq_u32(vreinterpretq_u32_u16(w4.val[1]), vreinterpretq_u32_u16(w5.val[1]));
    w14 = vzipq_u32(vreinterpretq_u32_u16(w12.val[0]), vreinterpretq_u32_u16(w13.val[0]));
    w15 = vzipq_u32(vreinterpretq_u32_u16(w12.val[1]), vreinterpretq_u32_u16(w13.val[1]));

    // Store second 4-line result

    d[12].val[0] = vzip1q_u64(vreinterpretq_u64_u32(w6.val[0]), vreinterpretq_u64_u32(w14.val[0]));
    d[12].val[1] = vzip2q_u64(vreinterpretq_u64_u32(w6.val[0]), vreinterpretq_u64_u32(w14.val[0]));
    d[13].val[0] = vzip1q_u64(vreinterpretq_u64_u32(w6.val[1]), vreinterpretq_u64_u32(w14.val[1]));
    d[13].val[1] = vzip2q_u64(vreinterpretq_u64_u32(w6.val[1]), vreinterpretq_u64_u32(w14.val[1]));
    d[14].val[0] = vzip1q_u64(vreinterpretq_u64_u32(w7.val[0]), vreinterpretq_u64_u32(w15.val[0]));
    d[14].val[1] = vzip2q_u64(vreinterpretq_u64_u32(w7.val[0]), vreinterpretq_u64_u32(w15.val[0]));
    d[15].val[0] = vzip1q_u64(vreinterpretq_u64_u32(w7.val[1]), vreinterpretq_u64_u32(w15.val[1]));
    d[15].val[1] = vzip2q_u64(vreinterpretq_u64_u32(w7.val[1]), vreinterpretq_u64_u32(w15.val[1]));
}

static INLINE void transpose_TX_16X16(const uint8_t *src, ptrdiff_t pitchSrc, uint8_t *dst, ptrdiff_t pitchDst) {
    uint8x16_t r[16];
    uint64x2_t d[16];
    for (int i = 0; i < 16; i++) { r[i] = vld1q_u8(src + i * pitchSrc); }
    transpose16x16_neon(r, d);
    for (int i = 0; i < 16; i++) { vst1q_u8(dst + i * pitchDst, vreinterpretq_u8_u64(d[i])); }
}

static INLINE void transpose(const uint8_t *src, ptrdiff_t pitchSrc, uint8_t *dst, ptrdiff_t pitchDst, int width,
                             int height) {
    for (int j = 0; j < height; j += 16) {
        for (int i = 0; i < width; i += 16) {
            transpose_TX_16X16(src + i * pitchSrc + j, pitchSrc, dst + j * pitchDst + i, pitchDst);
        }
    }
}

static INLINE void dr_prediction_z3_4x4_neon(uint8_t *dst, ptrdiff_t stride, const uint8_t *left, int upsample_left,
                                             int dy) {
    uint8x8_t    dstvec[4];
    uint16x4x2_t dest;

    dr_prediction_z1_HxW_internal_neon_64(4, 4, dstvec, left, upsample_left, dy);
    transpose4x8_8x4_low_neon(dstvec, &dest);
    vst1_lane_u32((uint32_t *)(dst + stride * 0), vreinterpret_u32_u16(dest.val[0]), 0);
    vst1_lane_u32((uint32_t *)(dst + stride * 1), vreinterpret_u32_u16(dest.val[0]), 1);
    vst1_lane_u32((uint32_t *)(dst + stride * 2), vreinterpret_u32_u16(dest.val[1]), 0);
    vst1_lane_u32((uint32_t *)(dst + stride * 3), vreinterpret_u32_u16(dest.val[1]), 1);
}

static INLINE void dr_prediction_z3_8x8_neon(uint8_t *dst, ptrdiff_t stride, const uint8_t *left, int upsample_left,
                                             int dy) {
    uint8x8_t    dstvec[8];
    uint32x2x2_t d[4];

    dr_prediction_z1_HxW_internal_neon_64(8, 8, dstvec, left, upsample_left, dy);
    transpose8x8_neon(dstvec, d);
    vst1_u32((uint32_t *)(dst + 0 * stride), d[0].val[0]);
    vst1_u32((uint32_t *)(dst + 1 * stride), d[0].val[1]);
    vst1_u32((uint32_t *)(dst + 2 * stride), d[1].val[0]);
    vst1_u32((uint32_t *)(dst + 3 * stride), d[1].val[1]);
    vst1_u32((uint32_t *)(dst + 4 * stride), d[2].val[0]);
    vst1_u32((uint32_t *)(dst + 5 * stride), d[2].val[1]);
    vst1_u32((uint32_t *)(dst + 6 * stride), d[3].val[0]);
    vst1_u32((uint32_t *)(dst + 7 * stride), d[3].val[1]);
}

static INLINE void dr_prediction_z3_4x8_neon(uint8_t *dst, ptrdiff_t stride, const uint8_t *left, int upsample_left,
                                             int dy) {
    uint8x8_t    dstvec[4];
    uint16x4x2_t d[2];

    dr_prediction_z1_HxW_internal_neon_64(8, 4, dstvec, left, upsample_left, dy);
    transpose4x8_8x4_neon(dstvec, d);
    vst1_lane_u32((uint32_t *)(dst + stride * 0), vreinterpret_u32_u16(d[0].val[0]), 0);
    vst1_lane_u32((uint32_t *)(dst + stride * 1), vreinterpret_u32_u16(d[0].val[0]), 1);
    vst1_lane_u32((uint32_t *)(dst + stride * 2), vreinterpret_u32_u16(d[0].val[1]), 0);
    vst1_lane_u32((uint32_t *)(dst + stride * 3), vreinterpret_u32_u16(d[0].val[1]), 1);
    vst1_lane_u32((uint32_t *)(dst + stride * 4), vreinterpret_u32_u16(d[1].val[0]), 0);
    vst1_lane_u32((uint32_t *)(dst + stride * 5), vreinterpret_u32_u16(d[1].val[0]), 1);
    vst1_lane_u32((uint32_t *)(dst + stride * 6), vreinterpret_u32_u16(d[1].val[1]), 0);
    vst1_lane_u32((uint32_t *)(dst + stride * 7), vreinterpret_u32_u16(d[1].val[1]), 1);
}

static INLINE void dr_prediction_z3_8x4_neon(uint8_t *dst, ptrdiff_t stride, const uint8_t *left, int upsample_left,
                                             int dy) {
    uint8x8_t    dstvec[8];
    uint32x2x2_t d[2];

    dr_prediction_z1_HxW_internal_neon_64(4, 8, dstvec, left, upsample_left, dy);
    transpose8x8_low_neon(dstvec, d);
    vst1_u32((uint32_t *)(dst + 0 * stride), d[0].val[0]);
    vst1_u32((uint32_t *)(dst + 1 * stride), d[0].val[1]);
    vst1_u32((uint32_t *)(dst + 2 * stride), d[1].val[0]);
    vst1_u32((uint32_t *)(dst + 3 * stride), d[1].val[1]);
}

static INLINE void dr_prediction_z3_8x16_neon(uint8_t *dst, ptrdiff_t stride, const uint8_t *left, int upsample_left,
                                              int dy) {
    uint8x16_t dstvec[8];
    uint64x2_t d[8];

    dr_prediction_z1_HxW_internal_neon(16, 8, dstvec, left, upsample_left, dy);
    transpose8x16_16x8_neon(dstvec, d);
    for (int i = 0; i < 8; i++) {
        vst1_u8(dst + i * stride, vreinterpret_u8_u64(vget_low_u64(d[i])));
        vst1_u8(dst + (i + 8) * stride, vreinterpret_u8_u64(vget_high_u64(d[i])));
    }
}

static INLINE void dr_prediction_z3_16x8_neon(uint8_t *dst, ptrdiff_t stride, const uint8_t *left, int upsample_left,
                                              int dy) {
    uint8x8_t  dstvec[16];
    uint64x2_t d[8];

    dr_prediction_z1_HxW_internal_neon_64(8, 16, dstvec, left, upsample_left, dy);
    transpose16x8_8x16_neon(dstvec, d);
    for (int i = 0; i < 8; i++) { vst1q_u8(dst + i * stride, vreinterpretq_u8_u64(d[i])); }
}

static INLINE void dr_prediction_z3_4x16_neon(uint8_t *dst, ptrdiff_t stride, const uint8_t *left, int upsample_left,
                                              int dy) {
    uint8x16_t   dstvec[4];
    uint16x8x2_t d[2];

    dr_prediction_z1_HxW_internal_neon(16, 4, dstvec, left, upsample_left, dy);
    transpose4x16_neon(dstvec, d);
    vst1q_lane_u32((uint32_t *)(dst + stride * 0), vreinterpretq_u32_u16(d[0].val[0]), 0);
    vst1q_lane_u32((uint32_t *)(dst + stride * 1), vreinterpretq_u32_u16(d[0].val[0]), 1);
    vst1q_lane_u32((uint32_t *)(dst + stride * 2), vreinterpretq_u32_u16(d[0].val[0]), 2);
    vst1q_lane_u32((uint32_t *)(dst + stride * 3), vreinterpretq_u32_u16(d[0].val[0]), 3);

    vst1q_lane_u32((uint32_t *)(dst + stride * 4), vreinterpretq_u32_u16(d[0].val[1]), 0);
    vst1q_lane_u32((uint32_t *)(dst + stride * 5), vreinterpretq_u32_u16(d[0].val[1]), 1);
    vst1q_lane_u32((uint32_t *)(dst + stride * 6), vreinterpretq_u32_u16(d[0].val[1]), 2);
    vst1q_lane_u32((uint32_t *)(dst + stride * 7), vreinterpretq_u32_u16(d[0].val[1]), 3);

    vst1q_lane_u32((uint32_t *)(dst + stride * 8), vreinterpretq_u32_u16(d[1].val[0]), 0);
    vst1q_lane_u32((uint32_t *)(dst + stride * 9), vreinterpretq_u32_u16(d[1].val[0]), 1);
    vst1q_lane_u32((uint32_t *)(dst + stride * 10), vreinterpretq_u32_u16(d[1].val[0]), 2);
    vst1q_lane_u32((uint32_t *)(dst + stride * 11), vreinterpretq_u32_u16(d[1].val[0]), 3);

    vst1q_lane_u32((uint32_t *)(dst + stride * 12), vreinterpretq_u32_u16(d[1].val[1]), 0);
    vst1q_lane_u32((uint32_t *)(dst + stride * 13), vreinterpretq_u32_u16(d[1].val[1]), 1);
    vst1q_lane_u32((uint32_t *)(dst + stride * 14), vreinterpretq_u32_u16(d[1].val[1]), 2);
    vst1q_lane_u32((uint32_t *)(dst + stride * 15), vreinterpretq_u32_u16(d[1].val[1]), 3);
}

static INLINE void dr_prediction_z3_16x4_neon(uint8_t *dst, ptrdiff_t stride, const uint8_t *left, int upsample_left,
                                              int dy) {
    uint8x8_t  dstvec[16];
    uint64x2_t d[8];

    dr_prediction_z1_HxW_internal_neon_64(4, 16, dstvec, left, upsample_left, dy);
    transpose16x8_8x16_neon(dstvec, d);
    for (int i = 0; i < 4; i++) { vst1q_u8(dst + i * stride, vreinterpretq_u8_u64(d[i])); }
}

static INLINE void dr_prediction_z3_8x32_neon(uint8_t *dst, ptrdiff_t stride, const uint8_t *left, int dy) {
    uint8x16x2_t dstvec[16];
    uint64x2x2_t d[16];
    uint8x16_t   v_zero = vdupq_n_u8(0);
    dr_prediction_z1_32xN_internal_neon(8, dstvec, left, dy);
    for (int i = 8; i < 16; i++) {
        dstvec[i].val[0] = v_zero;
        dstvec[i].val[1] = v_zero;
    }
    transpose16x32_neon(dstvec, d);
    for (int i = 0; i < 16; i++) {
        vst1_u8(dst + 2 * i * stride, vreinterpret_u8_u64(vget_low_u64(d[i].val[0])));
        vst1_u8(dst + (2 * i + 1) * stride, vreinterpret_u8_u64(vget_low_u64(d[i].val[1])));
    }
}

static INLINE void dr_prediction_z3_32x8_neon(uint8_t *dst, ptrdiff_t stride, const uint8_t *left, int upsample_left,
                                              int dy) {
    uint8x8_t  dstvec[32];
    uint64x2_t d[16];

    dr_prediction_z1_HxW_internal_neon_64(8, 32, dstvec, left, upsample_left, dy);
    transpose16x8_8x16_neon(dstvec, d);
    transpose16x8_8x16_neon(dstvec + 16, d + 8);
    for (int i = 0; i < 8; i++) {
        vst1q_u8(dst + i * stride, vreinterpretq_u8_u64(d[i]));
        vst1q_u8(dst + i * stride + 16, vreinterpretq_u8_u64(d[i + 8]));
    }
}

static INLINE void dr_prediction_z3_16x16_neon(uint8_t *dst, ptrdiff_t stride, const uint8_t *left, int upsample_left,
                                               int dy) {
    uint8x16_t dstvec[16];
    uint64x2_t d[16];

    dr_prediction_z1_HxW_internal_neon(16, 16, dstvec, left, upsample_left, dy);
    transpose16x16_neon(dstvec, d);
    for (int i = 0; i < 16; i++) { vst1q_u8(dst + i * stride, vreinterpretq_u8_u64(d[i])); }
}

static INLINE void dr_prediction_z3_32x32_neon(uint8_t *dst, ptrdiff_t stride, const uint8_t *left, int dy) {
    uint8x16x2_t dstvec[32];
    uint64x2x2_t d[32];
    dr_prediction_z1_32xN_internal_neon(32, dstvec, left, dy);
    transpose16x32_neon(dstvec, d);
    transpose16x32_neon(dstvec + 16, d + 16);
    for (int i = 0; i < 16; i++) {
        vst1q_u8(dst + 2 * i * stride, vreinterpretq_u8_u64(d[i].val[0]));
        vst1q_u8(dst + 2 * i * stride + 16, vreinterpretq_u8_u64(d[i + 16].val[0]));
        vst1q_u8(dst + (2 * i + 1) * stride, vreinterpretq_u8_u64(d[i].val[1]));
        vst1q_u8(dst + (2 * i + 1) * stride + 16, vreinterpretq_u8_u64(d[i + 16].val[1]));
    }
}

static INLINE void dr_prediction_z3_64x64_neon(uint8_t *dst, ptrdiff_t stride, const uint8_t *left, int dy) {
    DECLARE_ALIGNED(16, uint8_t, dstT[64 * 64]);

    dr_prediction_z1_64xN_neon(64, dstT, 64, left, dy);
    transpose(dstT, 64, dst, stride, 64, 64);
}

static INLINE void dr_prediction_z3_16x32_neon(uint8_t *dst, ptrdiff_t stride, const uint8_t *left, int dy) {
    uint8x16x2_t dstvec[16];
    uint64x2x2_t d[16];
    dr_prediction_z1_32xN_internal_neon(16, dstvec, left, dy);
    transpose16x32_neon(dstvec, d);
    for (int i = 0; i < 16; i++) {
        vst1q_u8(dst + 2 * i * stride, vreinterpretq_u8_u64(d[i].val[0]));
        vst1q_u8(dst + (2 * i + 1) * stride, vreinterpretq_u8_u64(d[i].val[1]));
    }
}

static INLINE void dr_prediction_z3_32x16_neon(uint8_t *dst, ptrdiff_t stride, const uint8_t *left, int upsample_left,
                                               int dy) {
    uint8x16_t dstvec[32];
    uint64x2_t d[16];

    dr_prediction_z1_HxW_internal_neon(16, 32, dstvec, left, upsample_left, dy);
    for (int i = 0; i < 32; i += 16) {
        transpose16x16_neon((dstvec + i), d);
        for (int j = 0; j < 16; j++) { vst1q_u8(dst + j * stride + i, vreinterpretq_u8_u64(d[j])); }
    }
}

static INLINE void dr_prediction_z3_32x64_neon(uint8_t *dst, ptrdiff_t stride, const uint8_t *left, int dy) {
    uint8_t dstT[64 * 32];

    dr_prediction_z1_64xN_neon(32, dstT, 64, left, dy);
    transpose(dstT, 64, dst, stride, 32, 64);
}

static INLINE void dr_prediction_z3_64x32_neon(uint8_t *dst, ptrdiff_t stride, const uint8_t *left, int dy) {
    uint8_t dstT[32 * 64];

    dr_prediction_z1_32xN_neon(64, dstT, 32, left, dy);
    transpose(dstT, 32, dst, stride, 64, 32);
}

static INLINE void dr_prediction_z3_16x64_neon(uint8_t *dst, ptrdiff_t stride, const uint8_t *left, int dy) {
    uint8_t dstT[64 * 16];

    dr_prediction_z1_64xN_neon(16, dstT, 64, left, dy);
    transpose(dstT, 64, dst, stride, 16, 64);
}

static INLINE void dr_prediction_z3_64x16_neon(uint8_t *dst, ptrdiff_t stride, const uint8_t *left, int upsample_left,
                                               int dy) {
    uint8x16_t dstvec[64];
    uint64x2_t d[16];

    dr_prediction_z1_HxW_internal_neon(16, 64, dstvec, left, upsample_left, dy);
    for (int i = 0; i < 64; i += 16) {
        transpose16x16_neon((dstvec + i), d);
        for (int j = 0; j < 16; j++) { vst1q_u8(dst + j * stride + i, vreinterpretq_u8_u64(d[j])); }
    }
}

void svt_av1_dr_prediction_z3_neon(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh, const uint8_t *above,
                                   const uint8_t *left, int32_t upsample_left, int32_t dx, int32_t dy) {
    (void)above;
    (void)dx;
    assert(dx == 1);
    assert(dy > 0);

    if (bw == bh) {
        switch (bw) {
        case 4: {
            dr_prediction_z3_4x4_neon(dst, stride, left, upsample_left, dy);
            break;
        }
        case 8: {
            dr_prediction_z3_8x8_neon(dst, stride, left, upsample_left, dy);
            break;
        }
        case 16: {
            dr_prediction_z3_16x16_neon(dst, stride, left, upsample_left, dy);
            break;
        }
        case 32: {
            dr_prediction_z3_32x32_neon(dst, stride, left, dy);
            break;
        }
        case 64: {
            dr_prediction_z3_64x64_neon(dst, stride, left, dy);
            break;
        }
        }
    } else {
        if (bw < bh) {
            if (bw + bw == bh) {
                switch (bw) {
                case 4: {
                    dr_prediction_z3_4x8_neon(dst, stride, left, upsample_left, dy);
                    break;
                }
                case 8: {
                    dr_prediction_z3_8x16_neon(dst, stride, left, upsample_left, dy);
                    break;
                }
                case 16: {
                    dr_prediction_z3_16x32_neon(dst, stride, left, dy);
                    break;
                }
                case 32: {
                    dr_prediction_z3_32x64_neon(dst, stride, left, dy);
                    break;
                }
                }
            } else {
                switch (bw) {
                case 4: {
                    dr_prediction_z3_4x16_neon(dst, stride, left, upsample_left, dy);
                    break;
                }
                case 8: {
                    dr_prediction_z3_8x32_neon(dst, stride, left, dy);
                    break;
                }
                case 16: {
                    dr_prediction_z3_16x64_neon(dst, stride, left, dy);
                    break;
                }
                }
            }
        } else {
            if (bh + bh == bw) {
                switch (bh) {
                case 4: {
                    dr_prediction_z3_8x4_neon(dst, stride, left, upsample_left, dy);
                    break;
                }
                case 8: {
                    dr_prediction_z3_16x8_neon(dst, stride, left, upsample_left, dy);
                    break;
                }
                case 16: {
                    dr_prediction_z3_32x16_neon(dst, stride, left, upsample_left, dy);
                    break;
                }
                case 32: {
                    dr_prediction_z3_64x32_neon(dst, stride, left, dy);
                    break;
                }
                }
            } else {
                switch (bh) {
                case 4: {
                    dr_prediction_z3_16x4_neon(dst, stride, left, upsample_left, dy);
                    break;
                }
                case 8: {
                    dr_prediction_z3_32x8_neon(dst, stride, left, upsample_left, dy);
                    break;
                }
                case 16: {
                    dr_prediction_z3_64x16_neon(dst, stride, left, upsample_left, dy);
                    break;
                }
                }
            }
        }
    }
}
