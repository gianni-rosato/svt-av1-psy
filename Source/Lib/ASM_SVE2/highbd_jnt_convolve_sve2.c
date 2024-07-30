/*
 * Copyright (c) 2024, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
 */

#include <arm_neon.h>

#include "common_dsp_rtcd.h"
#include "definitions.h"
#include "highbd_jnt_convolve_neon.h"
#include "highbd_convolve_sve2.h"
#include "mem_neon.h"
#include "neon_sve_bridge.h"

static INLINE uint16x4_t highbd_convolve8_4_y(int16x8_t samples_lo[2], int16x8_t samples_hi[2], int16x8_t filter,
                                              int64x2_t offset) {
    int64x2_t sum01 = svt_svdot_lane_s16(offset, samples_lo[0], filter, 0);
    sum01           = svt_svdot_lane_s16(sum01, samples_hi[0], filter, 1);

    int64x2_t sum23 = svt_svdot_lane_s16(offset, samples_lo[1], filter, 0);
    sum23           = svt_svdot_lane_s16(sum23, samples_hi[1], filter, 1);

    int32x4_t sum0123 = vcombine_s32(vmovn_s64(sum01), vmovn_s64(sum23));

    return vqrshrun_n_s32(sum0123, ROUND0_BITS);
}

static INLINE uint16x8_t highbd_convolve8_8_y(int16x8_t samples_lo[4], int16x8_t samples_hi[4], int16x8_t filter,
                                              int64x2_t offset) {
    int64x2_t sum01 = svt_svdot_lane_s16(offset, samples_lo[0], filter, 0);
    sum01           = svt_svdot_lane_s16(sum01, samples_hi[0], filter, 1);

    int64x2_t sum23 = svt_svdot_lane_s16(offset, samples_lo[1], filter, 0);
    sum23           = svt_svdot_lane_s16(sum23, samples_hi[1], filter, 1);

    int64x2_t sum45 = svt_svdot_lane_s16(offset, samples_lo[2], filter, 0);
    sum45           = svt_svdot_lane_s16(sum45, samples_hi[2], filter, 1);

    int64x2_t sum67 = svt_svdot_lane_s16(offset, samples_lo[3], filter, 0);
    sum67           = svt_svdot_lane_s16(sum67, samples_hi[3], filter, 1);

    int32x4_t sum0123 = vcombine_s32(vmovn_s64(sum01), vmovn_s64(sum23));
    int32x4_t sum4567 = vcombine_s32(vmovn_s64(sum45), vmovn_s64(sum67));

    return vcombine_u16(vqrshrun_n_s32(sum0123, ROUND0_BITS), vqrshrun_n_s32(sum4567, ROUND0_BITS));
}

static INLINE void highbd_dist_wtd_convolve_y_8tap_sve2(const uint16_t *src, int src_stride, uint16_t *dst,
                                                        int dst_stride, int width, int height,
                                                        const int16_t *y_filter_ptr, const int bd) {
    const int64x2_t offset   = vdupq_n_s64((1 << (bd + FILTER_BITS)) + (1 << (bd + FILTER_BITS - 1)));
    const int16x8_t y_filter = vld1q_s16(y_filter_ptr);

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
        int16_t  *s = (int16_t *)src;
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

            uint16x4_t d0 = highbd_convolve8_4_y(s0123, s4567, y_filter, offset);
            uint16x4_t d1 = highbd_convolve8_4_y(s1234, s5678, y_filter, offset);
            uint16x4_t d2 = highbd_convolve8_4_y(s2345, s6789, y_filter, offset);
            uint16x4_t d3 = highbd_convolve8_4_y(s3456, s789A, y_filter, offset);

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

                uint16x8_t d0 = highbd_convolve8_8_y(s0123, s4567, y_filter, offset);
                uint16x8_t d1 = highbd_convolve8_8_y(s1234, s5678, y_filter, offset);
                uint16x8_t d2 = highbd_convolve8_8_y(s2345, s6789, y_filter, offset);
                uint16x8_t d3 = highbd_convolve8_8_y(s3456, s789A, y_filter, offset);

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

void svt_av1_highbd_jnt_convolve_y_sve2(const uint16_t *src, int src_stride, uint16_t *dst, int dst_stride, int w,
                                        int h, const InterpFilterParams *filter_params_x,
                                        const InterpFilterParams *filter_params_y, const int subpel_x_qn,
                                        const int subpel_y_qn, ConvolveParams *conv_params, int bd) {
    if (h == 2 || w == 2) {
        svt_av1_highbd_jnt_convolve_y_c(src,
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

    DECLARE_ALIGNED(16, uint16_t, im_block[(MAX_SB_SIZE + MAX_FILTER_TAP) * MAX_SB_SIZE]);
    CONV_BUF_TYPE *dst16         = conv_params->dst;
    const int      y_filter_taps = get_filter_tap(filter_params_y, subpel_y_qn);

    if (y_filter_taps != 8) {
        svt_av1_highbd_jnt_convolve_y_neon(src,
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

    int       dst16_stride = conv_params->dst_stride;
    const int im_stride    = MAX_SB_SIZE;
    const int vert_offset  = filter_params_y->taps / 2 - 1;
    assert(FILTER_BITS == COMPOUND_ROUND1_BITS);

    const int16_t *y_filter_ptr = av1_get_interp_filter_subpel_kernel(*filter_params_y, subpel_y_qn & SUBPEL_MASK);

    src -= vert_offset * src_stride;

    if (conv_params->do_average) {
        highbd_dist_wtd_convolve_y_8tap_sve2(src, src_stride, im_block, im_stride, w, h, y_filter_ptr, bd);
        if (conv_params->use_jnt_comp_avg) {
            highbd_jnt_comp_avg_neon(im_block, im_stride, dst, dst_stride, w, h, conv_params, bd);
        } else {
            highbd_comp_avg_neon(im_block, im_stride, dst, dst_stride, w, h, conv_params, bd);
        }
    } else {
        highbd_dist_wtd_convolve_y_8tap_sve2(src, src_stride, dst16, dst16_stride, w, h, y_filter_ptr, bd);
    }
}
