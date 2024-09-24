/*
 * Copyright (c) 2024, Alliance for Open Media. All rights reserved.
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <arm_neon.h>
#include <arm_sve.h>

#include "aom_dsp_rtcd.h"
#include "definitions.h"
#include "mem_neon.h"
#include "neon_sve_bridge.h"
#include "pickrst_neon.h"
#include "pickrst_sve.h"
#include "restoration.h"
#include "restoration_pick.h"
#include "sum_neon.h"
#include "transpose_neon.h"
#include "utility.h"

static inline uint8_t find_average_sve(const uint8_t *src, int src_stride, int width, int height) {
    uint32x4_t avg_u32 = vdupq_n_u32(0);
    uint8x16_t ones    = vdupq_n_u8(1);

    // Use a predicate to compute the last columns.
    svbool_t pattern = svwhilelt_b8_u32(0, width % 16);

    int h = height;
    do {
        int            j       = width;
        const uint8_t *src_ptr = src;
        while (j >= 16) {
            uint8x16_t s = vld1q_u8(src_ptr);
            avg_u32      = vdotq_u32(avg_u32, s, ones);

            j -= 16;
            src_ptr += 16;
        }
        uint8x16_t s_end = svget_neonq_u8(svld1_u8(pattern, src_ptr));
        avg_u32          = vdotq_u32(avg_u32, s_end, ones);

        src += src_stride;
    } while (--h != 0);
    return (uint8_t)(vaddlvq_u32(avg_u32) / (width * height));
}

static inline void compute_sub_avg(const uint8_t *buf, int buf_stride, int avg, int16_t *buf_avg, int buf_avg_stride,
                                   int width, int height) {
    uint8x8_t avg_u8 = vdup_n_u8(avg);

    // Use a predicate to compute the last columns.
    svbool_t pattern = svwhilelt_b8_u32(0, width % 8);

    uint8x8_t avg_end = vget_low_u8(svget_neonq_u8(svdup_n_u8_z(pattern, avg)));

    do {
        int            j           = width;
        const uint8_t *buf_ptr     = buf;
        int16_t       *buf_avg_ptr = buf_avg;
        while (j >= 8) {
            uint8x8_t d = vld1_u8(buf_ptr);
            vst1q_s16(buf_avg_ptr, vreinterpretq_s16_u16(vsubl_u8(d, avg_u8)));

            j -= 8;
            buf_ptr += 8;
            buf_avg_ptr += 8;
        }
        uint8x8_t d_end = vget_low_u8(svget_neonq_u8(svld1_u8(pattern, buf_ptr)));
        vst1q_s16(buf_avg_ptr, vreinterpretq_s16_u16(vsubl_u8(d_end, avg_end)));

        buf += buf_stride;
        buf_avg += buf_avg_stride;
    } while (--height > 0);
}

void svt_av1_compute_stats_sve(int32_t wiener_win, const uint8_t *dgd, const uint8_t *src, int32_t h_start,
                               int32_t h_end, int32_t v_start, int32_t v_end, int32_t dgd_stride, int32_t src_stride,
                               int64_t *M, int64_t *H) {
    const int32_t wiener_win2    = wiener_win * wiener_win;
    const int32_t wiener_halfwin = (wiener_win >> 1);
    const int32_t width          = h_end - h_start;
    const int32_t height         = v_end - v_start;
    const int32_t d_stride       = (width + 2 * wiener_halfwin + 15) & ~15;
    const int32_t s_stride       = (width + 15) & ~15;
    int16_t      *d, *s;

    const uint8_t *dgd_start = dgd + h_start + v_start * dgd_stride;
    const uint8_t *src_start = src + h_start + v_start * src_stride;
    const uint16_t avg       = find_average_sve(dgd_start, dgd_stride, width, height);

    // The maximum input size is width * height, which is
    // (9 / 4) * RESTORATION_UNITSIZE_MAX * RESTORATION_UNITSIZE_MAX. Enlarge to
    // 3 * RESTORATION_UNITSIZE_MAX * RESTORATION_UNITSIZE_MAX considering
    // paddings.
    d = svt_aom_memalign(32, sizeof(*d) * 6 * RESTORATION_UNITSIZE_MAX * RESTORATION_UNITSIZE_MAX);
    s = d + 3 * RESTORATION_UNITSIZE_MAX * RESTORATION_UNITSIZE_MAX;

    compute_sub_avg(src_start, src_stride, avg, s, s_stride, width, height);
    compute_sub_avg(dgd + (v_start - wiener_halfwin) * dgd_stride + h_start - wiener_halfwin,
                    dgd_stride,
                    avg,
                    d,
                    d_stride,
                    width + 2 * wiener_halfwin,
                    height + 2 * wiener_halfwin);

    if (wiener_win == WIENER_WIN) {
        compute_stats_win7_sve(d, d_stride, s, s_stride, width, height, M, H);
    } else if (wiener_win == WIENER_WIN_CHROMA) {
        compute_stats_win5_sve(d, d_stride, s, s_stride, width, height, M, H);
    } else {
        assert(wiener_win == WIENER_WIN_3TAP);
        compute_stats_win3_sve(d, d_stride, s, s_stride, width, height, M, H);
    }

    // H is a symmetric matrix, so we only need to fill out the upper triangle.
    // We can copy it down to the lower triangle outside the (i, j) loops.
    diagonal_copy_stats_neon(wiener_win2, H);

    svt_aom_free(d);
}
