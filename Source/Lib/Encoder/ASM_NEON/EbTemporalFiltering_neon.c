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

#include <assert.h>
#include <arm_neon.h>
#include "EbDefinitions.h"
#include "EbTemporalFiltering_constants.h"
#include "EbUtility.h"

/* value [i:0-15] (sqrt((float)i)*65536.0 */
static const uint32_t sqrt_array_fp16[16] = {0,
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
                                             253819};

/* Calc sqrt linear max error 10% */
static uint32_t sqrt_fast(uint32_t x) {
    if (x > 15) {
        const int log2_half = svt_log2f(x) >> 1;
        const int mul2      = log2_half << 1;
        int       base      = x >> (mul2 - 2);
        assert(base < 16);
        return sqrt_array_fp16[base] >> (17 - log2_half);
    }
    return sqrt_array_fp16[x] >> 16;
}

#define SSE_STRIDE (BW + 2)

static uint32_t calculate_squared_errors_sum_no_div_16x16_neon(const uint8_t *s, int s_stride, const uint8_t *p,
                                                               int p_stride) {
    int32x4_t sum_hi = vdupq_n_s32(0);
    int32x4_t sum_lo = vdupq_n_s32(0);

    const int32x2_t stride = vzip1_s32(vdup_n_s32(s_stride), vdup_n_s32(p_stride));

    for (unsigned int i = 0; i < 16; i++) {
        const int32x2_t offset = vmul_s32(stride, vdup_n_s32(i));

        const uint8x16_t s_8 = vld1q_u8(s + vget_lane_s32(offset, 0));
        const uint8x16_t p_8 = vld1q_u8(p + vget_lane_s32(offset, 1));

        const int16x8_t s_8_hi = vreinterpretq_s16_u16(vmovl_u8(vget_high_u8(s_8)));
        const int16x8_t p_8_hi = vreinterpretq_s16_u16(vmovl_u8(vget_high_u8(p_8)));
        const int16x8_t s_8_lo = vreinterpretq_s16_u16(vmovl_u8(vget_low_u8(s_8)));
        const int16x8_t p_8_lo = vreinterpretq_s16_u16(vmovl_u8(vget_low_u8(p_8)));

        const int16x8_t dif_hi = vsubq_s16(s_8_hi, p_8_hi);
        const int16x8_t dif_lo = vsubq_s16(s_8_lo, p_8_lo);

        const int32x4_t pl_hi = vmull_s16(vget_low_s16(dif_hi), vget_low_s16(dif_hi));
        const int32x4_t ph_hi = vmull_high_s16(dif_hi, dif_hi);
        const int32x4_t pl_lo = vmull_s16(vget_low_s16(dif_lo), vget_low_s16(dif_lo));
        const int32x4_t ph_lo = vmull_high_s16(dif_lo, dif_lo);

        sum_hi = vaddq_s32(sum_hi, vpaddq_s32(pl_hi, ph_hi));
        sum_lo = vaddq_s32(sum_lo, vpaddq_s32(pl_lo, ph_lo));
    }

    int32x4_t sum = vpaddq_s32(sum_hi, sum_lo);

    sum = vpaddq_s32(sum, sum);
    sum = vpaddq_s32(sum, sum);

    return vgetq_lane_s32(sum, 0);
}

// T[X] =  exp(-(X)/16)  for x in [0..7], step 1/16 values in Fixed Points shift 16
static const int32_t expf_tab_fp16[] = {
    65536, 61565, 57835, 54331, 51039, 47947, 45042, 42313, 39749, 37341, 35078, 32953, 30957, 29081, 27319,
    25664, 24109, 22648, 21276, 19987, 18776, 17638, 16570, 15566, 14623, 13737, 12904, 12122, 11388, 10698,
    10050, 9441,  8869,  8331,  7827,  7352,  6907,  6488,  6095,  5726,  5379,  5053,  4747,  4459,  4189,
    3935,  3697,  3473,  3262,  3065,  2879,  2704,  2541,  2387,  2242,  2106,  1979,  1859,  1746,  1640,
    1541,  1447,  1360,  1277,  1200,  1127,  1059,  995,   934,   878,   824,   774,   728,   683,   642,
    603,   566,   532,   500,   470,   441,   414,   389,   366,   343,   323,   303,   285,   267,   251,
    236,   222,   208,   195,   184,   172,   162,   152,   143,   134,   126,   118,   111,   104,   98,
    92,    86,    81,    76,    72,    67,    63,    59,    56,    52,    49,    46,    43,    41,    38,
    36,    34,    31,    30,    28,    26,    24,    23,    21};

static void calculate_squared_errors_sum_2x8x8_no_div_neon(const uint8_t *s, int s_stride, const uint8_t *p,
                                                           int p_stride, uint32_t *output) {
    int32x4_t sum_lo = vdupq_n_s32(0);
    int32x4_t sum_hi = vdupq_n_s32(0);

    const int32x2_t stride = vzip1_s32(vdup_n_s32(s_stride), vdup_n_s32(p_stride));

    for (unsigned int i = 0; i < 8; i++) {
        const int32x2_t offset = vmul_s32(stride, vdup_n_s32(i));

        const uint8x16_t s_8 = vld1q_u8(s + vget_lane_s32(offset, 0));
        const uint8x16_t p_8 = vld1q_u8(p + vget_lane_s32(offset, 1));

        const int16x8_t s_8_lo = vreinterpretq_s16_u16(vmovl_u8(vget_low_u8(s_8)));
        const int16x8_t p_8_lo = vreinterpretq_s16_u16(vmovl_u8(vget_low_u8(p_8)));

        const int16x8_t s_8_hi = vreinterpretq_s16_u16(vmovl_u8(vget_high_u8(s_8)));
        const int16x8_t p_8_hi = vreinterpretq_s16_u16(vmovl_u8(vget_high_u8(p_8)));

        const int16x8_t dif_lo = vsubq_s16(s_8_lo, p_8_lo);
        const int16x8_t dif_hi = vsubq_s16(s_8_hi, p_8_hi);

        const int32x4_t pl_lo = vmull_s16(vget_low_s16(dif_lo), vget_low_s16(dif_lo));
        const int32x4_t ph_lo = vmull_high_s16(dif_lo, dif_lo);
        sum_lo                = vaddq_s32(sum_lo, vpaddq_s32(pl_lo, ph_lo));

        const int32x4_t pl_hi = vmull_s16(vget_low_s16(dif_hi), vget_low_s16(dif_hi));
        const int32x4_t ph_hi = vmull_high_s16(dif_hi, dif_hi);
        sum_hi                = vaddq_s32(sum_hi, vpaddq_s32(pl_hi, ph_hi));
    }
    sum_lo = vpaddq_s32(sum_lo, sum_lo);
    sum_lo = vpaddq_s32(sum_lo, sum_lo);

    sum_hi = vpaddq_s32(sum_hi, sum_hi);
    sum_hi = vpaddq_s32(sum_hi, sum_hi);

    output[0] = vgetq_lane_s32(sum_lo, 0);
    output[1] = vgetq_lane_s32(sum_hi, 0);
}

static void svt_av1_apply_temporal_filter_planewise_medium_partial_neon(
    struct MeContext *me_ctx, const uint8_t *y_src, int y_src_stride, const uint8_t *y_pre, int y_pre_stride,
    unsigned int block_width, unsigned int block_height, uint32_t *y_accum, uint16_t *y_count, uint32_t tf_decay_factor,
    uint32_t luma_window_error_quad_fp8[4], int is_chroma) {
    unsigned int i, j, k, subblock_idx;

    int32_t  idx_32x32               = me_ctx->tf_block_col + me_ctx->tf_block_row * 2;
    uint32_t distance_threshold_fp16 = AOMMAX((me_ctx->tf_mv_dist_th << 16) / 10, 1 << 16);
    // Calculation for every quarter
    uint32_t  d_factor_fp8[4];
    uint32_t  block_error_fp8[4];
    uint32_t  chroma_window_error_quad_fp8[4];
    uint32_t *window_error_quad_fp8 = is_chroma ? chroma_window_error_quad_fp8 : luma_window_error_quad_fp8;

    if (me_ctx->tf_32x32_block_split_flag[idx_32x32]) {
        const int32x4_t col = vmovl_s16(vld1_s16((int16_t *)&me_ctx->tf_16x16_mv_x[idx_32x32 * 4]));
        const int32x4_t row = vmovl_s16(vld1_s16((int16_t *)&me_ctx->tf_16x16_mv_y[idx_32x32 * 4]));

        const uint32x4_t hyp     = vreinterpretq_u32_s32(vaddq_s32(vmulq_s32(col, col), vmulq_s32(row, row)));
        const uint32x4_t hyp_256 = vshlq_n_u32(hyp, 8);

        uint32x4_t distance_fp4 = vcombine_u32(vzip1_u32(vdup_n_u32(sqrt_fast(vgetq_lane_u32(hyp_256, 0))),
                                                         vdup_n_u32(sqrt_fast(vgetq_lane_u32(hyp_256, 1)))),
                                               vzip1_u32(vdup_n_u32(sqrt_fast(vgetq_lane_u32(hyp_256, 2))),
                                                         vdup_n_u32(sqrt_fast(vgetq_lane_u32(hyp_256, 3)))));

        uint32x4_t d_factor_fp8_v = vcvtq_u32_f32(
            vdivq_f32(vcvtq_f32_u32(vshlq_n_u32(distance_fp4, 12)),
                      vcvtq_f32_u32(vshrq_n_u32(vdupq_n_u32(distance_threshold_fp16), 8))));
        d_factor_fp8_v = vmaxq_u32(d_factor_fp8_v, vdupq_n_u32(1 << 8));
        vst1q_u32(d_factor_fp8, d_factor_fp8_v);

        // ignore odd elements, since those are the higher 32 bits of every 64 bit entry
        uint32x4x2_t aux = vld2q_u32((uint32_t *)&me_ctx->tf_16x16_block_error[idx_32x32 * 4 + 0]);
        vst1q_u32(block_error_fp8, aux.val[0]);

    } else {
        tf_decay_factor <<= 1;
        int32_t col = me_ctx->tf_32x32_mv_x[idx_32x32];
        int32_t row = me_ctx->tf_32x32_mv_y[idx_32x32];

        uint32_t distance_fp4 = sqrt_fast(((uint32_t)(col * col + row * row)) << 8);
        d_factor_fp8[0] = d_factor_fp8[1] = d_factor_fp8[2] = d_factor_fp8[3] = AOMMAX(
            (distance_fp4 << 12) / (distance_threshold_fp16 >> 8), 1 << 8);
        FP_ASSERT(me_ctx->tf_32x32_block_error[idx_32x32] < ((uint64_t)1 << 30));
        block_error_fp8[0] = block_error_fp8[1] = block_error_fp8[2] = block_error_fp8[3] =
            (uint32_t)(me_ctx->tf_32x32_block_error[idx_32x32] >> 2);
    }

    if (block_width == 32) {
        window_error_quad_fp8[0] = calculate_squared_errors_sum_no_div_16x16_neon(
            y_src, y_src_stride, y_pre, y_pre_stride);
        window_error_quad_fp8[1] = calculate_squared_errors_sum_no_div_16x16_neon(
            y_src + 16, y_src_stride, y_pre + 16, y_pre_stride);
        window_error_quad_fp8[2] = calculate_squared_errors_sum_no_div_16x16_neon(
            y_src + y_src_stride * 16, y_src_stride, y_pre + y_pre_stride * 16, y_pre_stride);
        window_error_quad_fp8[3] = calculate_squared_errors_sum_no_div_16x16_neon(
            y_src + y_src_stride * 16 + 16, y_src_stride, y_pre + y_pre_stride * 16 + 16, y_pre_stride);
    } else {
        calculate_squared_errors_sum_2x8x8_no_div_neon(y_src, y_src_stride, y_pre, y_pre_stride, window_error_quad_fp8);

        calculate_squared_errors_sum_2x8x8_no_div_neon(
            y_src + y_src_stride * 8, y_src_stride, y_pre + y_pre_stride * 8, y_pre_stride, &window_error_quad_fp8[2]);

        uint32x4_t window_error_quad_fp8_v = vld1q_u32(window_error_quad_fp8);
        window_error_quad_fp8_v            = vshlq_n_u32(window_error_quad_fp8_v, 2);

        vst1q_u32(window_error_quad_fp8, window_error_quad_fp8_v);
    }

    if (is_chroma) {
        for (i = 0; i < 4; ++i) {
            FP_ASSERT(((int64_t)window_error_quad_fp8[i] * 5 + luma_window_error_quad_fp8[i]) < ((int64_t)1 << 31));
        }

        uint32x4_t window_error_quad_fp8_v      = vld1q_u32(window_error_quad_fp8);
        uint32x4_t luma_window_error_quad_fp8_v = vld1q_u32(luma_window_error_quad_fp8);

        window_error_quad_fp8_v = vmlaq_u32(luma_window_error_quad_fp8_v, window_error_quad_fp8_v, vdupq_n_u32(5));
        window_error_quad_fp8_v = vcvtq_u32_f32(vdivq_f32(vcvtq_f32_u32(window_error_quad_fp8_v), vdupq_n_f32(6.0f)));

        vst1q_u32(window_error_quad_fp8, window_error_quad_fp8_v);
    }

    int16x8_t adjusted_weight_int16[4];
    int32x4_t adjusted_weight_int32[4];

    for (subblock_idx = 0; subblock_idx < 4; subblock_idx++) {
        uint32_t combined_error_fp8 = (window_error_quad_fp8[subblock_idx] * TF_WINDOW_BLOCK_BALANCE_WEIGHT +
                                       block_error_fp8[subblock_idx]) /
            (TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1);

        uint64_t avg_err_fp10    = ((combined_error_fp8 >> 3) * (d_factor_fp8[subblock_idx] >> 3));
        uint32_t scaled_diff16   = (uint32_t)AOMMIN((avg_err_fp10) / AOMMAX((tf_decay_factor >> 10), 1), 7 * 16);
        uint32_t adjusted_weight = (expf_tab_fp16[scaled_diff16] * TF_WEIGHT_SCALE) >> 16;

        adjusted_weight_int16[subblock_idx] = vdupq_n_s16((int16_t)(adjusted_weight));
        adjusted_weight_int32[subblock_idx] = vdupq_n_s32((int32_t)(adjusted_weight));
    }

    for (i = 0; i < block_height; i++) {
        const int subblock_idx_h = (i >= block_height / 2) * 2;
        for (j = 0; j < block_width; j += 8) {
            k = i * y_pre_stride + j;

            uint16x8_t count_array = vld1q_u16(y_count + k);

            count_array = vaddq_u16(
                count_array, vreinterpretq_u16_s16(adjusted_weight_int16[subblock_idx_h + (j >= block_width / 2)]));

            vst1q_u16(y_count + k, count_array);
            uint32x4_t accumulator_array1 = vld1q_u32(y_accum + k);
            uint32x4_t accumulator_array2 = vld1q_u32(y_accum + k + 4);

            uint16x8_t frame2_array       = vmovl_u8(vld1_u8(y_pre + k));
            uint32x4_t frame2_array_u32_1 = vmovl_u16(vget_low_u16(frame2_array));
            uint32x4_t frame2_array_u32_2 = vmovl_u16(vget_high_u16(frame2_array));

            uint32x4_t adj_weight = vreinterpretq_u32_s32(
                adjusted_weight_int32[subblock_idx_h + (j >= block_width / 2)]);

            frame2_array_u32_1 = vmulq_u32(frame2_array_u32_1, adj_weight);
            frame2_array_u32_2 = vmulq_u32(frame2_array_u32_2, adj_weight);

            accumulator_array1 = vaddq_u32(accumulator_array1, frame2_array_u32_1);
            accumulator_array2 = vaddq_u32(accumulator_array2, frame2_array_u32_2);

            vst1q_u32(y_accum + k, accumulator_array1);
            vst1q_u32(y_accum + k + 4, accumulator_array2);
        }
    }
}

void svt_av1_apply_temporal_filter_planewise_medium_neon(
    struct MeContext *me_ctx, const uint8_t *y_src, int y_src_stride, const uint8_t *y_pre, int y_pre_stride,
    const uint8_t *u_src, const uint8_t *v_src, int uv_src_stride, const uint8_t *u_pre, const uint8_t *v_pre,
    int uv_pre_stride, unsigned int block_width, unsigned int block_height, int ss_x, int ss_y, uint32_t *y_accum,
    uint16_t *y_count, uint32_t *u_accum, uint16_t *u_count, uint32_t *v_accum, uint16_t *v_count) {
    uint32_t luma_window_error_quad_fp8[4];

    svt_av1_apply_temporal_filter_planewise_medium_partial_neon(me_ctx,
                                                                y_src,
                                                                y_src_stride,
                                                                y_pre,
                                                                y_pre_stride,
                                                                (unsigned int)block_width,
                                                                (unsigned int)block_height,
                                                                y_accum,
                                                                y_count,
                                                                me_ctx->tf_decay_factor_fp16[C_Y],
                                                                luma_window_error_quad_fp8,
                                                                0);

    if (me_ctx->tf_chroma) {
        svt_av1_apply_temporal_filter_planewise_medium_partial_neon(me_ctx,
                                                                    u_src,
                                                                    uv_src_stride,
                                                                    u_pre,
                                                                    uv_pre_stride,
                                                                    (unsigned int)block_width >> ss_x,
                                                                    (unsigned int)block_height >> ss_y,
                                                                    u_accum,
                                                                    u_count,
                                                                    me_ctx->tf_decay_factor_fp16[C_U],
                                                                    luma_window_error_quad_fp8,
                                                                    1);

        svt_av1_apply_temporal_filter_planewise_medium_partial_neon(me_ctx,
                                                                    v_src,
                                                                    uv_src_stride,
                                                                    v_pre,
                                                                    uv_pre_stride,
                                                                    (unsigned int)block_width >> ss_x,
                                                                    (unsigned int)block_height >> ss_y,
                                                                    v_accum,
                                                                    v_count,
                                                                    me_ctx->tf_decay_factor_fp16[C_V],
                                                                    luma_window_error_quad_fp8,
                                                                    1);
    }
}

// Divide two int32x4 vectors
static uint32x4_t div_u32(const uint32x4_t *a, const uint32x4_t *b) {
    uint32x4_t result = vdupq_n_u32(0);
    result            = vsetq_lane_u32(vdups_laneq_u32(*a, 0) / vdups_laneq_u32(*b, 0), result, 0);
    result            = vsetq_lane_u32(vdups_laneq_u32(*a, 1) / vdups_laneq_u32(*b, 1), result, 1);
    result            = vsetq_lane_u32(vdups_laneq_u32(*a, 2) / vdups_laneq_u32(*b, 2), result, 2);
    result            = vsetq_lane_u32(vdups_laneq_u32(*a, 3) / vdups_laneq_u32(*b, 3), result, 3);
    return result;
}

static void process_block_hbd_neon(int h, int w, uint16_t *buff_hbd_start, uint32_t *accum, uint16_t *count,
                                   uint32_t stride) {
    int i, j;
    int pos = 0;
    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j += 8) {
            //buff_lbd_start[pos] = (uint8_t)OD_DIVU(accum[k] + (count[k] >> 1), count[k]);
            //buff_lbd_start[pos] = (uint8_t)((accum[k] + (count[k] >> 1))/ count[k]);
            uint32x4_t accum_a = vld1q_u32(accum);
            uint32x4_t accum_b = vld1q_u32(&accum[4]);
            accum += 8;

            uint32x4_t count_a = vmovl_u16(vld1_u16(count));
            uint32x4_t count_b = vmovl_u16(vld1_u16(&count[4]));
            count += 8;

            //accum[k] + (count[k] >> 1)
            uint32x4_t tmp_a = vaddq_u32(accum_a, vshrq_n_u32(count_a, 1));
            uint32x4_t tmp_b = vaddq_u32(accum_b, vshrq_n_u32(count_b, 1));

            //accum[k] + (count[k] >> 1))/ count[k]
            tmp_a             = div_u32(&tmp_a, &count_a);
            tmp_b             = div_u32(&tmp_b, &count_b);
            uint16x8_t tmp_ab = vqmovn_high_u32(vqmovn_u32(tmp_a), tmp_b);

            vst1q_u16(buff_hbd_start + pos, tmp_ab);

            pos += 8;
        }
        pos += stride;
    }
}

static void process_block_lbd_neon(int h, int w, uint8_t *buff_lbd_start, uint32_t *accum, uint16_t *count,
                                   uint32_t stride) {
    int i, j;
    int pos = 0;
    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j += 8) {
            //buff_lbd_start[pos] = (uint8_t)OD_DIVU(accum[k] + (count[k] >> 1), count[k]);
            //buff_lbd_start[pos] = (uint8_t)((accum[k] + (count[k] >> 1))/ count[k]);
            uint32x4_t accum_a = vld1q_u32(accum);
            uint32x4_t accum_b = vld1q_u32(&accum[4]);
            accum += 8;

            uint32x4_t count_a = vmovl_u16(vld1_u16(count));
            uint32x4_t count_b = vmovl_u16(vld1_u16(&count[4]));
            count += 8;

            //accum[k] + (count[k] >> 1)
            uint32x4_t tmp_a = vaddq_u32(accum_a, vshrq_n_u32(count_a, 1));
            uint32x4_t tmp_b = vaddq_u32(accum_b, vshrq_n_u32(count_b, 1));

            //accum[k] + (count[k] >> 1))/ count[k]
            tmp_a              = div_u32(&tmp_a, &count_a);
            tmp_b              = div_u32(&tmp_b, &count_b);
            uint16x8_t tmp_ab1 = vqmovn_high_u32(vqmovn_u32(tmp_a), tmp_b);

            vst1_u8(buff_lbd_start + pos, vqmovn_u16(tmp_ab1));

            pos += 8;
        }
        pos += stride;
    }
}

void svt_aom_get_final_filtered_pixels_neon(MeContext *me_ctx, EbByte *src_center_ptr_start,
                                            uint16_t **altref_buffer_highbd_start, uint32_t **accum, uint16_t **count,
                                            const uint32_t *stride, int blk_y_src_offset, int blk_ch_src_offset,
                                            uint16_t blk_width_ch, uint16_t blk_height_ch, Bool is_highbd) {
    assert(blk_width_ch % 16 == 0);
    assert(BW % 16 == 0);

    if (!is_highbd) {
        //Process luma
        process_block_lbd_neon(
            BH, BW, &src_center_ptr_start[C_Y][blk_y_src_offset], accum[C_Y], count[C_Y], stride[C_Y] - BW);
        // Process chroma
        if (me_ctx->tf_chroma) {
            process_block_lbd_neon(blk_height_ch,
                                   blk_width_ch,
                                   &src_center_ptr_start[C_U][blk_ch_src_offset],
                                   accum[C_U],
                                   count[C_U],
                                   stride[C_U] - blk_width_ch);
            process_block_lbd_neon(blk_height_ch,
                                   blk_width_ch,
                                   &src_center_ptr_start[C_V][blk_ch_src_offset],
                                   accum[C_V],
                                   count[C_V],
                                   stride[C_V] - blk_width_ch);
        }
    } else {
        // Process luma
        process_block_hbd_neon(
            BH, BW, &altref_buffer_highbd_start[C_Y][blk_y_src_offset], accum[C_Y], count[C_Y], stride[C_Y] - BW);
        // Process chroma
        if (me_ctx->tf_chroma) {
            process_block_hbd_neon(blk_height_ch,
                                   blk_width_ch,
                                   &altref_buffer_highbd_start[C_U][blk_ch_src_offset],
                                   accum[C_U],
                                   count[C_U],
                                   stride[C_U] - blk_width_ch);
            process_block_hbd_neon(blk_height_ch,
                                   blk_width_ch,
                                   &altref_buffer_highbd_start[C_V][blk_ch_src_offset],
                                   accum[C_V],
                                   count[C_V],
                                   stride[C_V] - blk_width_ch);
        }
    }
}
