/*
* Copyright (c) 2024, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at www.aomedia.org/license/software. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at www.aomedia.org/license/patent.
*/

#include <arm_neon.h>

#include "common_dsp_rtcd.h"
#include "definitions.h"
#include "filter.h"
#include "inter_prediction.h"
#include "mem_neon.h"
#include "neon_sve_bridge.h"
#include "neon_sve2_bridge.h"
#include "sum_neon.h"
#include "utility.h"

// clang-format off
DECLARE_ALIGNED(16, static const uint16_t, kDotProdMergeBlockTbl[24]) = {
    // Shift left and insert new last column in transposed 4x4 block.
    1, 2, 3, 0, 5, 6, 7, 4,
    // Shift left and insert two new columns in transposed 4x4 block.
    2, 3, 0, 1, 6, 7, 4, 5,
    // Shift left and insert three new columns in transposed 4x4 block.
    3, 0, 1, 2, 7, 4, 5, 6,
};
// clang-format on

static INLINE void transpose_concat_4x4(int16x4_t s0, int16x4_t s1, int16x4_t s2, int16x4_t s3, int16x8_t res[2]) {
    // Transpose 16-bit elements and concatenate result rows as follows:
    // s0: 00, 01, 02, 03
    // s1: 10, 11, 12, 13
    // s2: 20, 21, 22, 23
    // s3: 30, 31, 32, 33
    //
    // res[0]: 00 10 20 30 01 11 21 31
    // res[1]: 02 12 22 32 03 13 23 33

    int16x8_t s0q = vcombine_s16(s0, vdup_n_s16(0));
    int16x8_t s1q = vcombine_s16(s1, vdup_n_s16(0));
    int16x8_t s2q = vcombine_s16(s2, vdup_n_s16(0));
    int16x8_t s3q = vcombine_s16(s3, vdup_n_s16(0));

    int16x8_t s02 = vzip1q_s16(s0q, s2q);
    int16x8_t s13 = vzip1q_s16(s1q, s3q);

    int16x8x2_t s0123 = vzipq_s16(s02, s13);

    res[0] = s0123.val[0];
    res[1] = s0123.val[1];
}

static INLINE void transpose_concat_8x4(int16x8_t s0, int16x8_t s1, int16x8_t s2, int16x8_t s3, int16x8_t res[4]) {
    // Transpose 16-bit elements and concatenate result rows as follows:
    // s0: 00, 01, 02, 03, 04, 05, 06, 07
    // s1: 10, 11, 12, 13, 14, 15, 16, 17
    // s2: 20, 21, 22, 23, 24, 25, 26, 27
    // s3: 30, 31, 32, 33, 34, 35, 36, 37
    //
    // res[0]: 00 10 20 30 01 11 21 31
    // res[1]: 02 12 22 32 03 13 23 33
    // res[2]: 04 14 24 34 05 15 25 35
    // res[3]: 06 16 26 36 07 17 27 37

    int16x8x2_t s02      = vzipq_s16(s0, s2);
    int16x8x2_t s13      = vzipq_s16(s1, s3);
    int16x8x2_t s0123_lo = vzipq_s16(s02.val[0], s13.val[0]);
    int16x8x2_t s0123_hi = vzipq_s16(s02.val[1], s13.val[1]);

    res[0] = s0123_lo.val[0];
    res[1] = s0123_lo.val[1];
    res[2] = s0123_hi.val[0];
    res[3] = s0123_hi.val[1];
}

static INLINE void svt_tbl2x4_s16(int16x8_t t0[4], int16x8_t t1[4], uint16x8_t tbl, int16x8_t res[4]) {
    res[0] = svt_tbl2_s16(t0[0], t1[0], tbl);
    res[1] = svt_tbl2_s16(t0[1], t1[1], tbl);
    res[2] = svt_tbl2_s16(t0[2], t1[2], tbl);
    res[3] = svt_tbl2_s16(t0[3], t1[3], tbl);
}

static INLINE void svt_tbl2x2_s16(int16x8_t t0[2], int16x8_t t1[2], uint16x8_t tbl, int16x8_t res[2]) {
    res[0] = svt_tbl2_s16(t0[0], t1[0], tbl);
    res[1] = svt_tbl2_s16(t0[1], t1[1], tbl);
}

static INLINE uint16x4_t highbd_convolve8_4_y(int16x8_t samples_lo[2], int16x8_t samples_hi[2], int16x8_t filter,
                                              uint16x4_t max) {
    int64x2_t sum01 = svt_svdot_lane_s16(vdupq_n_s64(0), samples_lo[0], filter, 0);
    sum01           = svt_svdot_lane_s16(sum01, samples_hi[0], filter, 1);

    int64x2_t sum23 = svt_svdot_lane_s16(vdupq_n_s64(0), samples_lo[1], filter, 0);
    sum23           = svt_svdot_lane_s16(sum23, samples_hi[1], filter, 1);

    int32x4_t  sum0123 = vcombine_s32(vmovn_s64(sum01), vmovn_s64(sum23));
    uint16x4_t res     = vqrshrun_n_s32(sum0123, FILTER_BITS);
    return vmin_u16(res, max);
}

static INLINE uint16x8_t highbd_convolve8_8_y(int16x8_t samples_lo[4], int16x8_t samples_hi[4], int16x8_t filter,
                                              uint16x8_t max) {
    int64x2_t sum01 = svt_svdot_lane_s16(vdupq_n_s64(0), samples_lo[0], filter, 0);
    sum01           = svt_svdot_lane_s16(sum01, samples_hi[0], filter, 1);

    int64x2_t sum23 = svt_svdot_lane_s16(vdupq_n_s64(0), samples_lo[1], filter, 0);
    sum23           = svt_svdot_lane_s16(sum23, samples_hi[1], filter, 1);

    int64x2_t sum45 = svt_svdot_lane_s16(vdupq_n_s64(0), samples_lo[2], filter, 0);
    sum45           = svt_svdot_lane_s16(sum45, samples_hi[2], filter, 1);

    int64x2_t sum67 = svt_svdot_lane_s16(vdupq_n_s64(0), samples_lo[3], filter, 0);
    sum67           = svt_svdot_lane_s16(sum67, samples_hi[3], filter, 1);

    int32x4_t  sum0123 = vcombine_s32(vmovn_s64(sum01), vmovn_s64(sum23));
    int32x4_t  sum4567 = vcombine_s32(vmovn_s64(sum45), vmovn_s64(sum67));
    uint16x8_t res     = vcombine_u16(vqrshrun_n_s32(sum0123, FILTER_BITS), vqrshrun_n_s32(sum4567, FILTER_BITS));
    return vminq_u16(res, max);
}

static void highbd_convolve_y_sr_8tap_sve2(const uint16_t *src, ptrdiff_t src_stride, uint16_t *dst,
                                           ptrdiff_t dst_stride, int width, int height, const int16_t *filter_y,
                                           int bd) {
    assert(width >= 4 && height >= 4);

    const int16x8_t y_filter = vld1q_s16(filter_y);

    uint16x8x3_t merge_block_tbl = vld1q_u16_x3(kDotProdMergeBlockTbl);
    // Scale indices by size of the true vector length to avoid reading from an
    // 'undefined' portion of a vector on a system with SVE vectors > 128-bit.
    uint16x8_t correction0 = vreinterpretq_u16_u64(vdupq_n_u64(svcnth() * 0x0001000000000000ULL));
    merge_block_tbl.val[0] = vaddq_u16(merge_block_tbl.val[0], correction0);

    uint16x8_t correction1 = vreinterpretq_u16_u64(vdupq_n_u64(svcnth() * 0x0001000100000000ULL));
    merge_block_tbl.val[1] = vaddq_u16(merge_block_tbl.val[1], correction1);

    uint16x8_t correction2 = vreinterpretq_u16_u64(vdupq_n_u64(svcnth() * 0x0001000100010000ULL));
    merge_block_tbl.val[2] = vaddq_u16(merge_block_tbl.val[2], correction2);

    if (width == 4) {
        const uint16x4_t max = vdup_n_u16((1 << bd) - 1);
        int16_t         *s   = (int16_t *)src;

        int16x4_t s0, s1, s2, s3, s4, s5, s6;
        load_s16_4x7(s, src_stride, &s0, &s1, &s2, &s3, &s4, &s5, &s6);
        s += 7 * src_stride;

        // This operation combines a conventional transpose and the sample permute
        // required before computing the dot product.
        int16x8_t s0123[2], s1234[2], s2345[2], s3456[2];
        transpose_concat_4x4(s0, s1, s2, s3, s0123);
        transpose_concat_4x4(s1, s2, s3, s4, s1234);
        transpose_concat_4x4(s2, s3, s4, s5, s2345);
        transpose_concat_4x4(s3, s4, s5, s6, s3456);

        do {
            int16x4_t s7, s8, s9, s10;
            load_s16_4x4(s, src_stride, &s7, &s8, &s9, &s10);

            int16x8_t s4567[2], s5678[2], s6789[2], s789A[2];
            // Transpose and shuffle the 4 lines that were loaded.
            transpose_concat_4x4(s7, s8, s9, s10, s789A);

            // Merge new data into block from previous iteration.
            svt_tbl2x2_s16(s3456, s789A, merge_block_tbl.val[0], s4567);
            svt_tbl2x2_s16(s3456, s789A, merge_block_tbl.val[1], s5678);
            svt_tbl2x2_s16(s3456, s789A, merge_block_tbl.val[2], s6789);

            uint16x4_t d0 = highbd_convolve8_4_y(s0123, s4567, y_filter, max);
            uint16x4_t d1 = highbd_convolve8_4_y(s1234, s5678, y_filter, max);
            uint16x4_t d2 = highbd_convolve8_4_y(s2345, s6789, y_filter, max);
            uint16x4_t d3 = highbd_convolve8_4_y(s3456, s789A, y_filter, max);

            store_u16_4x4(dst, dst_stride, d0, d1, d2, d3);

            // Prepare block for next iteration - re-using as much as possible.
            // Shuffle everything up four rows.
            s0123[0] = s4567[0];
            s0123[1] = s4567[1];
            s1234[0] = s5678[0];
            s1234[1] = s5678[1];
            s2345[0] = s6789[0];
            s2345[1] = s6789[1];
            s3456[0] = s789A[0];
            s3456[1] = s789A[1];
            s += 4 * src_stride;
            dst += 4 * dst_stride;
            height -= 4;
        } while (height != 0);
    } else {
        const uint16x8_t max = vdupq_n_u16((1 << bd) - 1);

        do {
            int       h = height;
            int16_t  *s = (int16_t *)src;
            uint16_t *d = dst;

            int16x8_t s0, s1, s2, s3, s4, s5, s6;
            load_s16_8x7(s, src_stride, &s0, &s1, &s2, &s3, &s4, &s5, &s6);
            s += 7 * src_stride;

            // This operation combines a conventional transpose and the sample permute
            // required before computing the dot product.
            int16x8_t s0123[4], s1234[4], s2345[4], s3456[4];
            transpose_concat_8x4(s0, s1, s2, s3, s0123);
            transpose_concat_8x4(s1, s2, s3, s4, s1234);
            transpose_concat_8x4(s2, s3, s4, s5, s2345);
            transpose_concat_8x4(s3, s4, s5, s6, s3456);

            do {
                int16x8_t s7, s8, s9, s10;
                load_s16_8x4(s, src_stride, &s7, &s8, &s9, &s10);

                int16x8_t s4567[4], s5678[4], s6789[4], s789A[4];
                // Transpose and shuffle the 4 lines that were loaded.
                transpose_concat_8x4(s7, s8, s9, s10, s789A);

                // Merge new data into block from previous iteration.
                svt_tbl2x4_s16(s3456, s789A, merge_block_tbl.val[0], s4567);
                svt_tbl2x4_s16(s3456, s789A, merge_block_tbl.val[1], s5678);
                svt_tbl2x4_s16(s3456, s789A, merge_block_tbl.val[2], s6789);

                uint16x8_t d0 = highbd_convolve8_8_y(s0123, s4567, y_filter, max);
                uint16x8_t d1 = highbd_convolve8_8_y(s1234, s5678, y_filter, max);
                uint16x8_t d2 = highbd_convolve8_8_y(s2345, s6789, y_filter, max);
                uint16x8_t d3 = highbd_convolve8_8_y(s3456, s789A, y_filter, max);

                store_u16_8x4(d, dst_stride, d0, d1, d2, d3);

                // Prepare block for next iteration - re-using as much as possible.
                // Shuffle everything up four rows.
                s0123[0] = s4567[0];
                s0123[1] = s4567[1];
                s0123[2] = s4567[2];
                s0123[3] = s4567[3];
                s1234[0] = s5678[0];
                s1234[1] = s5678[1];
                s1234[2] = s5678[2];
                s1234[3] = s5678[3];
                s2345[0] = s6789[0];
                s2345[1] = s6789[1];
                s2345[2] = s6789[2];
                s2345[3] = s6789[3];
                s3456[0] = s789A[0];
                s3456[1] = s789A[1];
                s3456[2] = s789A[2];
                s3456[3] = s789A[3];

                s += 4 * src_stride;
                d += 4 * dst_stride;
                h -= 4;
            } while (h != 0);
            src += 8;
            dst += 8;
            width -= 8;
        } while (width != 0);
    }
}

static INLINE uint16x4_t highbd_convolve4_4_y(int16x8_t samples[2], int16x8_t filter, uint16x4_t max) {
    int64x2_t sum01 = svt_svdot_lane_s16(vdupq_n_s64(0), samples[0], filter, 0);
    int64x2_t sum23 = svt_svdot_lane_s16(vdupq_n_s64(0), samples[1], filter, 0);

    int32x4_t  sum0123 = vcombine_s32(vmovn_s64(sum01), vmovn_s64(sum23));
    uint16x4_t res     = vqrshrun_n_s32(sum0123, FILTER_BITS);
    return vmin_u16(res, max);
}

static INLINE uint16x8_t highbd_convolve4_8_y(int16x8_t samples[4], int16x8_t filter, uint16x8_t max) {
    int64x2_t sum01 = svt_svdot_lane_s16(vdupq_n_s64(0), samples[0], filter, 0);
    int64x2_t sum23 = svt_svdot_lane_s16(vdupq_n_s64(0), samples[1], filter, 0);
    int64x2_t sum45 = svt_svdot_lane_s16(vdupq_n_s64(0), samples[2], filter, 0);
    int64x2_t sum67 = svt_svdot_lane_s16(vdupq_n_s64(0), samples[3], filter, 0);

    int32x4_t  sum0123 = vcombine_s32(vmovn_s64(sum01), vmovn_s64(sum23));
    int32x4_t  sum4567 = vcombine_s32(vmovn_s64(sum45), vmovn_s64(sum67));
    uint16x8_t res     = vcombine_u16(vqrshrun_n_s32(sum0123, FILTER_BITS), vqrshrun_n_s32(sum4567, FILTER_BITS));
    return vminq_u16(res, max);
}

static void highbd_convolve_y_sr_4tap_sve2(const uint16_t *src, ptrdiff_t src_stride, uint16_t *dst,
                                           ptrdiff_t dst_stride, int width, int height, const int16_t *filter_y,
                                           int bd) {
    assert(width >= 4 && height >= 4);

    const int16x8_t y_filter = vcombine_s16(vld1_s16(filter_y + 2), vdup_n_s16(0));

    if (width == 4) {
        const uint16x4_t max = vdup_n_u16((1 << bd) - 1);
        int16_t         *s   = (int16_t *)src;

        int16x4_t s0, s1, s2;
        load_s16_4x3(s, src_stride, &s0, &s1, &s2);
        s += 3 * src_stride;

        do {
            int16x4_t s3, s4, s5, s6;
            load_s16_4x4(s, src_stride, &s3, &s4, &s5, &s6);

            // This operation combines a conventional transpose and the sample permute
            // required before computing the dot product.
            int16x8_t s0123[2], s1234[2], s2345[2], s3456[2];
            transpose_concat_4x4(s0, s1, s2, s3, s0123);
            transpose_concat_4x4(s1, s2, s3, s4, s1234);
            transpose_concat_4x4(s2, s3, s4, s5, s2345);
            transpose_concat_4x4(s3, s4, s5, s6, s3456);

            uint16x4_t d0 = highbd_convolve4_4_y(s0123, y_filter, max);
            uint16x4_t d1 = highbd_convolve4_4_y(s1234, y_filter, max);
            uint16x4_t d2 = highbd_convolve4_4_y(s2345, y_filter, max);
            uint16x4_t d3 = highbd_convolve4_4_y(s3456, y_filter, max);

            store_u16_4x4(dst, dst_stride, d0, d1, d2, d3);

            // Shuffle everything up four rows.
            s0 = s4;
            s1 = s5;
            s2 = s6;

            s += 4 * src_stride;
            dst += 4 * dst_stride;
            height -= 4;
        } while (height != 0);
    } else {
        const uint16x8_t max = vdupq_n_u16((1 << bd) - 1);

        do {
            int       h = height;
            int16_t  *s = (int16_t *)src;
            uint16_t *d = dst;

            int16x8_t s0, s1, s2;
            load_s16_8x3(s, src_stride, &s0, &s1, &s2);
            s += 3 * src_stride;

            do {
                int16x8_t s3, s4, s5, s6;
                load_s16_8x4(s, src_stride, &s3, &s4, &s5, &s6);

                // This operation combines a conventional transpose and the sample
                // permute required before computing the dot product.
                int16x8_t s0123[4], s1234[4], s2345[4], s3456[4];
                transpose_concat_8x4(s0, s1, s2, s3, s0123);
                transpose_concat_8x4(s1, s2, s3, s4, s1234);
                transpose_concat_8x4(s2, s3, s4, s5, s2345);
                transpose_concat_8x4(s3, s4, s5, s6, s3456);

                uint16x8_t d0 = highbd_convolve4_8_y(s0123, y_filter, max);
                uint16x8_t d1 = highbd_convolve4_8_y(s1234, y_filter, max);
                uint16x8_t d2 = highbd_convolve4_8_y(s2345, y_filter, max);
                uint16x8_t d3 = highbd_convolve4_8_y(s3456, y_filter, max);

                store_u16_8x4(d, dst_stride, d0, d1, d2, d3);

                // Shuffle everything up four rows.
                s0 = s4;
                s1 = s5;
                s2 = s6;

                s += 4 * src_stride;
                d += 4 * dst_stride;
                h -= 4;
            } while (h != 0);
            src += 8;
            dst += 8;
            width -= 8;
        } while (width != 0);
    }
}

void svt_av1_highbd_convolve_y_sr_sve2(const uint16_t *src, int src_stride, uint16_t *dst, int dst_stride, int w, int h,
                                       const InterpFilterParams *filter_params_x,
                                       const InterpFilterParams *filter_params_y, const int subpel_x_qn,
                                       const int subpel_y_qn, ConvolveParams *conv_params, int bd) {
    if (w == 2 || h == 2) {
        svt_av1_highbd_convolve_y_sr_c(src,
                                       src_stride,
                                       dst,
                                       dst_stride,
                                       w,
                                       h,
                                       filter_params_x,
                                       filter_params_y,
                                       subpel_x_qn,
                                       subpel_y_qn,
                                       conv_params,
                                       bd);
        return;
    }
    const int y_filter_taps = get_filter_tap(filter_params_y, subpel_y_qn);

    if (y_filter_taps == 6) {
        svt_av1_highbd_convolve_y_sr_neon(src,
                                          src_stride,
                                          dst,
                                          dst_stride,
                                          w,
                                          h,
                                          filter_params_x,
                                          filter_params_y,
                                          subpel_x_qn,
                                          subpel_y_qn,
                                          conv_params,
                                          bd);
        return;
    }

    const int      vert_offset  = filter_params_y->taps / 2 - 1;
    const int16_t *y_filter_ptr = av1_get_interp_filter_subpel_kernel(*filter_params_y, subpel_y_qn & SUBPEL_MASK);

    src -= vert_offset * src_stride;

    if (y_filter_taps == 4) {
        highbd_convolve_y_sr_4tap_sve2(src + 2 * src_stride, src_stride, dst, dst_stride, w, h, y_filter_ptr, bd);
        return;
    }

    highbd_convolve_y_sr_8tap_sve2(src, src_stride, dst, dst_stride, w, h, y_filter_ptr, bd);
}
