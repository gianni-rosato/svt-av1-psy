/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

/*
* Copyright (c) 2016, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at www.aomedia.org/license/software. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at www.aomedia.org/license/patent.
*/

#include <stdlib.h>
#include <string.h>

#include "EbIntraPrediction.h"
#include "EbUtility.h"
#include "EbModeDecision.h"
#include "EbCodingUnit.h"
#include "EbModeDecisionProcess.h"
#include "EbDefinitions.h"
#include "EbPictureControlSet.h"
#include "EbEncDecProcess.h"
#include "aom_dsp_rtcd.h"

void *eb_aom_memset16(void *dest, int32_t val, size_t length);

int32_t is_inter_block(const BlockModeInfo *mbmi);

// Some basic checks on weights for smooth predictor.
#define sm_weights_sanity_checks(weights_w, weights_h, weights_scale, \
                                 pred_scale)                          \
  assert(weights_w[0] < weights_scale);                               \
  assert(weights_h[0] < weights_scale);                               \
  assert(weights_scale - weights_w[bw - 1] < weights_scale);          \
  assert(weights_scale - weights_h[bh - 1] < weights_scale);          \
  assert(pred_scale < 31)  // ensures no overflow when calculating predictor.
#define MIDRANGE_VALUE_8BIT    128
#define MIDRANGE_VALUE_10BIT   512
static PartitionType from_shape_to_part[] = {
    PARTITION_NONE,
    PARTITION_HORZ,
    PARTITION_VERT,
    PARTITION_HORZ_A,
    PARTITION_HORZ_B,
    PARTITION_VERT_A,
    PARTITION_VERT_B,
    PARTITION_HORZ_4,
    PARTITION_VERT_4,
    PARTITION_SPLIT
};
int is_smooth(const BlockModeInfo *block_mi, int plane) {
    if (plane == 0) {
        const PredictionMode mode = block_mi->mode;
        return (mode == SMOOTH_PRED || mode == SMOOTH_V_PRED ||
            mode == SMOOTH_H_PRED);
    }
    else {
        // uv_mode is not set for inter blocks, so need to explicitly
        // detect that case.
        if (is_inter_block(block_mi)) return 0;

        const UvPredictionMode uv_mode = block_mi->uv_mode;
        return (uv_mode == UV_SMOOTH_PRED || uv_mode == UV_SMOOTH_V_PRED ||
            uv_mode == UV_SMOOTH_H_PRED);
    }
}

static int get_filt_type(const MacroBlockD *xd, int plane) {
    int ab_sm, le_sm;

    if (plane == 0) {
        const MbModeInfo *ab = xd->above_mbmi;
        const MbModeInfo *le = xd->left_mbmi;
        ab_sm = ab ? is_smooth(&ab->block_mi, plane) : 0;
        le_sm = le ? is_smooth(&le->block_mi, plane) : 0;
    }
    else {
        const MbModeInfo *ab = xd->chroma_above_mbmi;
        const MbModeInfo *le = xd->chroma_left_mbmi;
        ab_sm = ab ? is_smooth(&ab->block_mi, plane) : 0;
        le_sm = le ? is_smooth(&le->block_mi, plane) : 0;
    }

    return (ab_sm || le_sm) ? 1 : 0;
}

int32_t use_intra_edge_upsample(int32_t bs0, int32_t bs1, int32_t delta, int32_t type) {
    const int32_t d = abs(delta);
    const int32_t blk_wh = bs0 + bs1;
    if (d <= 0 || d >= 40) return 0;
    return type ? (blk_wh <= 8) : (blk_wh <= 16);
}

#define INTRA_EDGE_FILT 3
#define INTRA_EDGE_TAPS 5
void eb_av1_filter_intra_edge_high_c_old(uint8_t *p, int32_t sz, int32_t strength) {
    if (!strength) return;

    const int32_t kernel[INTRA_EDGE_FILT][INTRA_EDGE_TAPS] = {
      { 0, 4, 8, 4, 0 }, { 0, 5, 6, 5, 0 }, { 2, 4, 4, 4, 2 }
    };
    const int32_t filt = strength - 1;
    uint8_t edge[129];

    memcpy(edge, p, sz * sizeof(*p));
    for (int32_t i = 1; i < sz; i++) {
        int32_t s = 0;
        for (int32_t j = 0; j < INTRA_EDGE_TAPS; j++) {
            int32_t k = i - 2 + j;
            k = (k < 0) ? 0 : k;
            k = (k > sz - 1) ? sz - 1 : k;
            s += edge[k] * kernel[filt][j];
        }
        s = (s + 8) >> 4;
        p[i] = (uint8_t)s;
    }
}
void eb_av1_filter_intra_edge_high_c_left(uint8_t *p, int32_t sz, int32_t strength, uint32_t                      size) {
    if (!strength) return;
    const uint32_t          leftBlockEnd = 2 * (size);

    const int32_t kernel[INTRA_EDGE_FILT][INTRA_EDGE_TAPS] = {
      { 0, 4, 8, 4, 0 }, { 0, 5, 6, 5, 0 }, { 2, 4, 4, 4, 2 }
    };
    const int32_t filt = strength - 1;
    uint8_t edge[129];
    edge[0] = p[0];
    p = p - leftBlockEnd;

    memcpy(edge + 1, p, (sz - 1) * sizeof(*p));

    for (int32_t i = 1; i < sz; i++) {
        int32_t s = 0;
        for (int32_t j = 0; j < INTRA_EDGE_TAPS; j++) {
            int32_t k = i - 2 + j;
            k = (k < 0) ? 0 : k;
            k = (k > sz - 1) ? sz - 1 : k;
            s += edge[k] * kernel[filt][j];
        }
        s = (s + 8) >> 4;
        p[i - 1] = (uint8_t)s;
    }
}

int32_t intra_edge_filter_strength(int32_t bs0, int32_t bs1, int32_t delta, int32_t type) {
    const int32_t d = abs(delta);
    int32_t strength = 0;

    const int32_t blk_wh = bs0 + bs1;
    if (type == 0) {
        if (blk_wh <= 8) {
            if (d >= 56)
                strength = 1;
        }
        else if (blk_wh <= 12) {
            if (d >= 40)
                strength = 1;
        }
        else if (blk_wh <= 16) {
            if (d >= 40)
                strength = 1;
        }
        else if (blk_wh <= 24) {
            if (d >= 8)
                strength = 1;
            if (d >= 16)
                strength = 2;
            if (d >= 32)
                strength = 3;
        }
        else if (blk_wh <= 32) {
            if (d >= 1)
                strength = 1;
            if (d >= 4)
                strength = 2;
            if (d >= 32)
                strength = 3;
        }
        else {
            if (d >= 1)
                strength = 3;
        }
    }
    else {
        if (blk_wh <= 8) {
            if (d >= 40)
                strength = 1;
            if (d >= 64)
                strength = 2;
        }
        else if (blk_wh <= 16) {
            if (d >= 20)
                strength = 1;
            if (d >= 48)
                strength = 2;
        }
        else if (blk_wh <= 24) {
            if (d >= 4)
                strength = 3;
        }
        else {
            if (d >= 1)
                strength = 3;
        }
    }
    return strength;
}

#define MAX_UPSAMPLE_SZ 16
void av1_upsample_intra_edge_high_c_old(uint8_t *p, int32_t sz, int32_t bd) {
    // interpolate half-sample positions
    assert(sz <= MAX_UPSAMPLE_SZ);

    uint8_t in[MAX_UPSAMPLE_SZ + 3];
    // copy p[-1..(sz-1)] and extend first and last samples
    in[0] = p[-1];
    in[1] = p[-1];
    for (int32_t i = 0; i < sz; i++)
        in[i + 2] = p[i];
    in[sz + 2] = p[sz - 1];

    // interpolate half-sample edge positions
    p[-2] = in[0];
    for (int32_t i = 0; i < sz; i++) {
        int32_t s = -in[i] + (9 * in[i + 1]) + (9 * in[i + 2]) - in[i + 3];
        s = (s + 8) >> 4;
        s = clip_pixel_highbd(s, bd);
        p[2 * i - 1] = (uint8_t)s;
        p[2 * i] = in[i + 2];
    }
}

const uint16_t eb_dr_intra_derivative[90] = {
    // More evenly spread out angles and limited to 10-bit
    // Values that are 0 will never be used
    //                    Approx angle
    0,    0, 0,        //
    1023, 0, 0,        // 3, ...
    547,  0, 0,        // 6, ...
    372,  0, 0, 0, 0,  // 9, ...
    273,  0, 0,        // 14, ...
    215,  0, 0,        // 17, ...
    178,  0, 0,        // 20, ...
    151,  0, 0,        // 23, ... (113 & 203 are base angles)
    132,  0, 0,        // 26, ...
    116,  0, 0,        // 29, ...
    102,  0, 0, 0,     // 32, ...
    90,   0, 0,        // 36, ...
    80,   0, 0,        // 39, ...
    71,   0, 0,        // 42, ...
    64,   0, 0,        // 45, ... (45 & 135 are base angles)
    57,   0, 0,        // 48, ...
    51,   0, 0,        // 51, ...
    45,   0, 0, 0,     // 54, ...
    40,   0, 0,        // 58, ...
    35,   0, 0,        // 61, ...
    31,   0, 0,        // 64, ...
    27,   0, 0,        // 67, ... (67 & 157 are base angles)
    23,   0, 0,        // 70, ...
    19,   0, 0,        // 73, ...
    15,   0, 0, 0, 0,  // 76, ...
    11,   0, 0,        // 81, ...
    7,    0, 0,        // 84, ...
    3,    0, 0,        // 87, ...
};

// Get the shift (up-scaled by 256) in Y w.r.t a unit change in X.
// If angle > 0 && angle < 90, dy = 1;
// If angle > 90 && angle < 180, dy = (int32_t)(256 * t);
// If angle > 180 && angle < 270, dy = -((int32_t)(256 * t));

#define divide_round(value, bits) (((value) + (1 << ((bits)-1))) >> (bits))

static INLINE uint16_t get_dy(int32_t angle) {
    if (angle > 90 && angle < 180)
        return eb_dr_intra_derivative[angle - 90];
    else if (angle > 180 && angle < 270)
        return eb_dr_intra_derivative[270 - angle];
    else {
        // In this case, we are not really going to use dy. We may return any value.
        return 1;
    }
}
// Get the shift (up-scaled by 256) in X w.r.t a unit change in Y.
// If angle > 0 && angle < 90, dx = -((int32_t)(256 / t));
// If angle > 90 && angle < 180, dx = (int32_t)(256 / t);
// If angle > 180 && angle < 270, dx = 1;
static INLINE uint16_t get_dx(int32_t angle) {
    if (angle > 0 && angle < 90)
        return eb_dr_intra_derivative[angle];
    else if (angle > 90 && angle < 180)
        return eb_dr_intra_derivative[180 - angle];
    else {
        // In this case, we are not really going to use dx. We may return any value.
        return 1;
    }
}

// Directional prediction, zone 3: 180 < angle < 270
void eb_av1_dr_prediction_z3_c(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh,
    const uint8_t *above, const uint8_t *left,
    int32_t upsample_left, int32_t dx, int32_t dy) {
    int32_t r, c, y, base, shift, val;

    (void)above;
    (void)dx;

    assert(dx == 1);
    assert(dy > 0);

    const int32_t max_base_y = (bw + bh - 1) << upsample_left;
    const int32_t frac_bits = 6 - upsample_left;
    const int32_t base_inc = 1 << upsample_left;
    y = dy;
    for (c = 0; c < bw; ++c, y += dy) {
        base = y >> frac_bits;
        shift = ((y << upsample_left) & 0x3F) >> 1;

        for (r = 0; r < bh; ++r, base += base_inc) {
            if (base < max_base_y) {
                val = left[base] * (32 - shift) + left[base + 1] * shift;
                val = ROUND_POWER_OF_TWO(val, 5);
                dst[r * stride + c] = (uint8_t)clip_pixel_highbd(val, 8);
            }
            else {
                for (; r < bh; ++r) dst[r * stride + c] = left[max_base_y];
                break;
            }
        }
    }
}
void eb_av1_dr_prediction_z1_c(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh,
    const uint8_t *above, const uint8_t *left,
    int32_t upsample_above, int32_t dx, int32_t dy) {
    int32_t r, c, x, base, shift, val;

    (void)left;
    (void)dy;
    assert(dy == 1);
    assert(dx > 0);

    const int32_t max_base_x = ((bw + bh) - 1) << upsample_above;
    const int32_t frac_bits = 6 - upsample_above;
    const int32_t base_inc = 1 << upsample_above;
    x = dx;
    for (r = 0; r < bh; ++r, dst += stride, x += dx) {
        base = x >> frac_bits;
        shift = ((x << upsample_above) & 0x3F) >> 1;

        if (base >= max_base_x) {
            for (int32_t i = r; i < bh; ++i) {
                memset(dst, above[max_base_x], bw * sizeof(dst[0]));
                dst += stride;
            }
            return;
        }

        for (c = 0; c < bw; ++c, base += base_inc) {
            if (base < max_base_x) {
                val = above[base] * (32 - shift) + above[base + 1] * shift;
                val = ROUND_POWER_OF_TWO(val, 5);
                dst[c] = (uint8_t)clip_pixel_highbd(val, 8);
            }
            else
                dst[c] = above[max_base_x];
        }
    }
}

// Directional prediction, zone 2: 90 < angle < 180
void eb_av1_dr_prediction_z2_c(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh,
    const uint8_t *above, const uint8_t *left,
    int32_t upsample_above, int32_t upsample_left, int32_t dx,
    int32_t dy) {
    int32_t r, c, x, y, shift1, shift2, val, base1, base2;

    assert(dx > 0);
    assert(dy > 0);

    const int32_t min_base_x = -(1 << upsample_above);
    const int32_t frac_bits_x = 6 - upsample_above;
    const int32_t frac_bits_y = 6 - upsample_left;
    const int32_t base_inc_x = 1 << upsample_above;
    x = -dx;
    for (r = 0; r < bh; ++r, x -= dx, dst += stride) {
        base1 = x >> frac_bits_x;
        y = (r << 6) - dy;
        for (c = 0; c < bw; ++c, base1 += base_inc_x, y -= dy) {
            if (base1 >= min_base_x) {
                shift1 = ((x * (1 << upsample_above)) & 0x3F) >> 1;
                val = above[base1] * (32 - shift1) + above[base1 + 1] * shift1;
                val = ROUND_POWER_OF_TWO(val, 5);
            }
            else {
                base2 = y >> frac_bits_y;
                assert(base2 >= -(1 << upsample_left));
                shift2 = ((y * (1 << upsample_left)) & 0x3F) >> 1;
                val = left[base2] * (32 - shift2) + left[base2 + 1] * shift2;
                val = ROUND_POWER_OF_TWO(val, 5);
            }
            dst[c] = (uint8_t)clip_pixel_highbd(val, 8);
        }
    }
}

/* clang-format on */
void intra_mode_planar(
    const uint32_t   size,                       //input parameter, denotes the size of the current PU
    uint8_t         *ref_samples,                 //input parameter, pointer to the reference samples
    uint8_t         *dst,              //output parameter, pointer to the prediction
    const uint32_t   prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
    const EbBool  skip)                       //skip half rows
{
    uint32_t rowStride = skip ? 2 : 1;
    uint32_t leftOffset = 0;
    uint32_t topOffset = (size << 1) + 1;

    const uint8_t below_pred = ref_samples[leftOffset + size - 1];   // estimated by bottom-left pixel
    const uint8_t right_pred = ref_samples[topOffset + size - 1];  // estimated by top-right pixel

    const uint8_t *const sm_weights_w = sm_weight_arrays + size;
    const uint8_t *const sm_weights_h = sm_weight_arrays + size;
    // scale = 2 * 2^sm_weight_log2_scale
    const int32_t log2_scale = 1 + sm_weight_log2_scale;
    const uint16_t scale = (1 << sm_weight_log2_scale);
    // sm_weights_sanity_checks(sm_weights_w, sm_weights_h, scale,  log2_scale + sizeof(*dst));
    uint32_t r;
    for (r = 0; r < size; ++r) {
        uint32_t c;
        for (c = 0; c < size; ++c) {
            const uint8_t pixels[] = { ref_samples[topOffset + c], below_pred,  ref_samples[leftOffset + r], right_pred };

            const uint8_t weights[] = { sm_weights_h[r], (uint8_t)(scale - sm_weights_h[r]),
                sm_weights_w[c], (uint8_t)(scale - sm_weights_w[c]) };
            uint32_t this_pred = 0;
            int32_t i;
            assert(scale >= sm_weights_h[r] && scale >= sm_weights_w[c]);
            for (i = 0; i < 4; ++i)
                this_pred += weights[i] * pixels[i];
            dst[c] = (uint8_t)divide_round(this_pred, log2_scale);
        }
        dst += rowStride * prediction_buffer_stride;
    }

    return;
}

/*static INLINE*/ void smooth_v_predictor_c(uint8_t *dst, ptrdiff_t stride, int32_t bw,
    int32_t bh, const uint8_t *above,
    const uint8_t *left) {
    const uint8_t below_pred = left[bh - 1];  // estimated by bottom-left pixel
    const uint8_t *const sm_weights = sm_weight_arrays + bh;
    // scale = 2^sm_weight_log2_scale
    const int32_t log2_scale = sm_weight_log2_scale;
    const uint16_t scale = (1 << sm_weight_log2_scale);
    sm_weights_sanity_checks(sm_weights, sm_weights, scale,
        log2_scale + sizeof(*dst));

    int32_t r;
    for (r = 0; r < bh; r++) {
        int32_t c;
        for (c = 0; c < bw; ++c) {
            const uint8_t pixels[] = { above[c], below_pred };
            const uint8_t weights[] = { (uint8_t)(sm_weights[r]), (uint8_t)(scale - sm_weights[r]) };
            uint32_t this_pred = 0;
            assert(scale >= sm_weights[r]);
            int32_t i;
            for (i = 0; i < 2; ++i)
                this_pred += weights[i] * pixels[i];
            dst[c] = (uint8_t)divide_round(this_pred, log2_scale);
        }
        dst += stride;
    }
}

/*static INLINE */void smooth_h_predictor_c(uint8_t *dst, ptrdiff_t stride, int32_t bw,
    int32_t bh, const uint8_t *above,
    const uint8_t *left) {
    const uint8_t right_pred = above[bw - 1];  // estimated by top-right pixel
    const uint8_t *const sm_weights = sm_weight_arrays + bw;
    // scale = 2^sm_weight_log2_scale
    const int32_t log2_scale = sm_weight_log2_scale;
    const uint16_t scale = (1 << sm_weight_log2_scale);
    sm_weights_sanity_checks(sm_weights, sm_weights, scale,
        log2_scale + sizeof(*dst));

    int32_t r;
    for (r = 0; r < bh; r++) {
        int32_t c;
        for (c = 0; c < bw; ++c) {
            const uint8_t pixels[] = { left[r], right_pred };
            const uint8_t weights[] = { (uint8_t)(sm_weights[c]), (uint8_t)(scale - sm_weights[c]) };
            uint32_t this_pred = 0;
            assert(scale >= sm_weights[c]);
            int32_t i;
            for (i = 0; i < 2; ++i)
                this_pred += weights[i] * pixels[i];
            dst[c] = (uint8_t)divide_round(this_pred, log2_scale);
        }
        dst += stride;
    }
}
/* clang-format on */
void ebav1_smooth_v_predictor(
    const uint32_t   size,                       //input parameter, denotes the size of the current PU
    uint8_t         *ref_samples,                 //input parameter, pointer to the reference samples
    uint8_t         *dst,              //output parameter, pointer to the prediction
    const uint32_t   prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
    const EbBool  skip)                       //skip half rows
{
    (void)skip;

    const uint32_t          bottomLeftEnd = 1 * (size);
    const uint32_t          topLeftBlockEnd = 2 * (size)+1;

    const uint8_t below_pred = ref_samples[bottomLeftEnd - 1];//left[size - 1];  // estimated by bottom-left pixel
    const uint8_t *const sm_weights = sm_weight_arrays + size;
    // scale = 2^sm_weight_log2_scale
    const int32_t log2_scale = sm_weight_log2_scale;
    const uint16_t scale = (1 << sm_weight_log2_scale);
    //sm_weights_sanity_checks(sm_weights, sm_weights, scale,
    //                        log2_scale + sizeof(*dst));

    uint16_t r;
    for (r = 0; r < size; r++) {
        uint16_t c;
        for (c = 0; c < size; ++c) {
            const uint8_t pixels[] = { ref_samples[topLeftBlockEnd + c], below_pred };
            const uint8_t weights[] = { (uint8_t)(sm_weights[r]), (uint8_t)(scale - sm_weights[r]) };
            uint32_t this_pred = 0;
            assert(scale >= sm_weights[r]);
            int32_t i;
            for (i = 0; i < 2; ++i)
                this_pred += weights[i] * pixels[i];
            dst[c] = (uint8_t)divide_round(this_pred, log2_scale);
        }
        dst += prediction_buffer_stride;
    }
    return;
}
void ebav1_smooth_h_predictor(
    const uint32_t   size,                       //input parameter, denotes the size of the current PU
    uint8_t         *ref_samples,                 //input parameter, pointer to the reference samples
    uint8_t         *dst,              //output parameter, pointer to the prediction
    const uint32_t   prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
    const EbBool  skip)
{
    (void)skip;
    const uint32_t          topLeftBlockEnd = 2 * (size)+1;

    const uint8_t right_pred = ref_samples[topLeftBlockEnd + size - 1];  // estimated by top-right pixel
    const uint8_t *const sm_weights = sm_weight_arrays + size;
    // scale = 2^sm_weight_log2_scale
    const int32_t log2_scale = sm_weight_log2_scale;
    const uint16_t scale = (1 << sm_weight_log2_scale);
    //sm_weights_sanity_checks(sm_weights, sm_weights, scale,
    //                         log2_scale + sizeof(*dst));

    uint16_t r;
    for (r = 0; r < size; r++) {
        uint16_t c;
        for (c = 0; c < size; ++c) {
            const uint8_t pixels[] = { ref_samples[r], right_pred };
            const uint8_t weights[] = { (uint8_t)(sm_weights[c]), (uint8_t)(scale - scale - sm_weights[c]) };
            uint32_t this_pred = 0;
            assert(scale >= sm_weights[c]);
            int32_t i;
            for (i = 0; i < 2; ++i)
                this_pred += weights[i] * pixels[i];
            dst[c] = (uint8_t)divide_round(this_pred, log2_scale);
        }
        dst += prediction_buffer_stride;
    }
    return;
}

void ebav1_v_predictor(uint8_t *dst, const uint32_t stride, int32_t bw, int32_t bh,
    const uint8_t *ref_samples) {
    int32_t r;
    int32_t c;

    for (r = 0; r < bh; r++) {
        for (c = 0; c < bh; c++) {
            dst[c + r * stride] = ref_samples[bw + bh + 1 + c];
            //EB_MEMSET(dst, ref_samples[bw + bh + 1], bw);
           //dst += stride;
        }
    }
}

void ebav1_h_predictor(uint8_t *dst, const uint32_t stride, int32_t bw, int32_t bh,
    const uint8_t *ref_samples) {
    int32_t r;
    //(void)above;

    for (r = 0; r < bh; r++) {
        EB_MEMSET(dst, ref_samples[r], bw);
        dst += stride;
    }

    return;
}

void IntraModeAngular_AV1_Z1(
    const uint32_t   size,                       //input parameter, denotes the size of the current PU
    uint8_t         *ref_samples,                 //input parameter, pointer to the reference samples
    uint8_t         *dst,              //output parameter, pointer to the prediction
    const uint32_t   prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
    const EbBool  skip,
    uint16_t          dx,              //output parameter, pointer to the prediction
    uint16_t          dy              //output parameter, pointer to the prediction
)

{
    (void)dy;

    //    uint32_t row_index, colIndex;
    uint32_t toAboveOffset = (size << 1) + 1;
    uint32_t rowStride = skip ? 2 : 1;

    uint32_t r, c;
    int32_t  x, base, shift, val;
    const int32_t max_base_x = ((size + size) - 1);
    const int32_t frac_bits = 6;
    const int32_t base_inc = 1;

    x = dx;

    for (r = 0; r < size; ++r, dst += (rowStride* prediction_buffer_stride), x += dx) {
        base = x >> frac_bits;
        shift = ((x) & 0x3F) >> 1;

        if (base >= max_base_x) {
            for (uint32_t i = r; i < size; ++i) {
                EB_MEMSET(dst, ref_samples[toAboveOffset + max_base_x], size);
                dst += rowStride * prediction_buffer_stride;
            }
            return;
        }

        for (c = 0; c < size; ++c, base += base_inc) {
            if (base < max_base_x) {
                val = ref_samples[toAboveOffset + base] * (32 - shift) + ref_samples[toAboveOffset + base + 1] * shift;
                val = ROUND_POWER_OF_TWO(val, 5);

                dst[c] = (uint8_t)clip_pixel_highbd(val, 8);
            }
            else
                dst[c] = ref_samples[toAboveOffset + max_base_x];
        }
    }

    return;
}
void IntraModeAngular_AV1_Z2(
    const uint32_t   size,                       //input parameter, denotes the size of the current PU
    uint8_t         *ref_samples,                 //input parameter, pointer to the reference samples
    uint8_t         *dst,              //output parameter, pointer to the prediction
    const uint32_t   prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
    const EbBool  skip,
    uint16_t          dx,              //output parameter, pointer to the prediction
    uint16_t          dy              //output parameter, pointer to the prediction
)
{
    //    uint32_t row_index, colIndex;
    uint32_t toAboveOffset = (size << 1) + 1;
    uint32_t toAboveLeftOffset = (size << 1);
    uint32_t shiftToAboveLeft = 0;
    uint32_t toLeftOffset = 0;
    uint32_t rowStride = skip ? 2 : 1;

    uint32_t r, c;
    int32_t x, y, shift, val, base;

    const int32_t min_base_x = -(1);
    const int32_t frac_bits_x = 6;
    const int32_t frac_bits_y = 6;
    for (r = 0; r < size; ++r) {
        for (c = 0; c < size; ++c) {
            y = r + 1;
            x = (c << 6) - y * dx;
            base = x >> frac_bits_x;
            if (base >= min_base_x) {
                shift = ((x) & 0x3F) >> 1;
                shiftToAboveLeft = (base <= -1) ? -1 - base : 0;
                val = ref_samples[toAboveOffset + base + shiftToAboveLeft] * (32 - shift) + ref_samples[toAboveOffset + base + 1] * shift;
                val = ROUND_POWER_OF_TWO(val, 5);
            }
            else {
                x = c + 1;
                y = (r << 6) - x * dy;
                base = y >> frac_bits_y;
                shiftToAboveLeft = (base <= -1) ? toAboveLeftOffset - base : 0;
                shift = ((y) & 0x3F) >> 1;
                val = ref_samples[toLeftOffset + base + shiftToAboveLeft] * (32 - shift) + ref_samples[toLeftOffset + base + 1] * shift;
                val = ROUND_POWER_OF_TWO(val, 5);
            }

            dst[c] = (uint8_t)clip_pixel_highbd(val, 8);
        }
        dst += (rowStride* prediction_buffer_stride);
    }

    return;
}
void IntraModeAngular_AV1_Z3(
    const uint32_t   size,                        //input parameter, denotes the size of the current PU
    uint8_t         *ref_samples,                  //input parameter, pointer to the reference samples
    uint8_t         *dst,                         //output parameter, pointer to the prediction
    const uint32_t   prediction_buffer_stride,      //input parameter, denotes the stride for the prediction ptr
    const EbBool  skip,
    uint16_t          dx,                          //output parameter, pointer to the prediction
    uint16_t          dy                           //output parameter, pointer to the prediction
)
{
    //    uint32_t row_index, colIndex;
    //    uint32_t toAboveOffset = (size << 1) + 1 ;
    uint32_t toLeftOffset = 0;
    uint32_t rowStride = skip ? 2 : 1;

    uint32_t r, c;
    int32_t y, base, shift, val;

    (void)dx;
    assert(dx == 1);
    assert(dy > 0);

    const int32_t max_base_y = (size + size - 1);
    const int32_t frac_bits = 6;
    const int32_t base_inc = 1;
    y = dy;
    for (c = 0; c < size; ++c, y += dy) {
        base = y >> frac_bits;
        shift = ((y) & 0x3F) >> 1;

        for (r = 0; r < size; ++r, base += base_inc) {
            if (base < max_base_y) {
                val = ref_samples[toLeftOffset + base] * (32 - shift) + ref_samples[toLeftOffset + base + 1] * shift;
                val = ROUND_POWER_OF_TWO(val, 5);

                dst[r * (rowStride* prediction_buffer_stride) + c] = (uint8_t)clip_pixel_highbd(val, 8);
            }
            else {
                for (; r < size; ++r) dst[r * (rowStride* prediction_buffer_stride) + c] = ref_samples[toLeftOffset + max_base_y];
                break;
            }
        }
    }

    return;
}

void highbd_dc_predictor_16bit(
    EbBool        is_left_availble,
    EbBool        is_above_availble,
    const uint32_t   size,                       //input parameter, denotes the size of the current PU
    uint16_t         *ref_samples,                 //input parameter, pointer to the reference samples
    uint16_t         *dst,              //output parameter, pointer to the prediction
    const uint32_t   prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
    const EbBool  skip)                       //skip half rows
{
    //uint32_t sum = 0;
//    uint32_t index;
    uint32_t columnIndex, row_index;
    uint32_t writeIndex;
    uint32_t leftOffset = 0;
    uint32_t topOffset = (size << 1) + 1;
    uint32_t rowStride = skip ? 2 : 1;

    uint32_t i;
    uint32_t  expected_dc, sum = 0;
    const int32_t count = size + size;

    if (is_left_availble && !is_above_availble) {
        for (i = 0; i < size; i++)
            sum += ref_samples[leftOffset + i];
        expected_dc = (sum + (size >> 1)) / size;
    }
    else if (is_above_availble && !is_left_availble) {
        for (i = 0; i < size; i++)
            sum += ref_samples[topOffset + i];
        expected_dc = (sum + (size >> 1)) / size;
    }
    else {
        for (i = 0; i < size; i++)
            sum += ref_samples[topOffset + i];
        for (i = 0; i < size; i++)
            sum += ref_samples[leftOffset + i];
        expected_dc = (sum + (count >> 1)) / count;
    }

    // expected_dc = (sum + (count >> 1)) / count;

    /*for (r = 0; r < size; r++) {
    writeIndex = row_index * prediction_buffer_stride;

    EB_MEMSET( dst, expected_dc, size);

    dst += rowStride* prediction_buffer_stride;
    }*/

    // Generate the prediction
    for (row_index = 0; row_index < size; row_index += rowStride) {
        writeIndex = row_index * prediction_buffer_stride;
        for (columnIndex = 0; columnIndex < size; ++columnIndex) {
            dst[writeIndex] = (uint16_t)expected_dc;
            ++writeIndex;
        }
    }

    return;
}
/* clang-format on */
void IntraModePlanar_16bit(
    const uint32_t   size,                       //input parameter, denotes the size of the current PU
    uint16_t         *ref_samples,                 //input parameter, pointer to the reference samples
    uint16_t         *dst,              //output parameter, pointer to the prediction
    const uint32_t   prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
    const EbBool  skip)                       //skip half rows
{
    uint32_t rowStride = skip ? 2 : 1;
    uint32_t leftOffset = 0;
    uint32_t topOffset = (size << 1) + 1;

    const uint16_t below_pred = ref_samples[leftOffset + size - 1];   // estimated by bottom-left pixel
    const uint16_t right_pred = ref_samples[topOffset + size - 1];  // estimated by top-right pixel
    const uint8_t *const sm_weights_w = sm_weight_arrays + size;
    const uint8_t *const sm_weights_h = sm_weight_arrays + size;
    // scale = 2 * 2^sm_weight_log2_scale
    const int32_t log2_scale = 1 + sm_weight_log2_scale;
    const uint16_t scale = (1 << sm_weight_log2_scale);
    // sm_weights_sanity_checks(sm_weights_w, sm_weights_h, scale,  log2_scale + sizeof(*dst));
    uint32_t r;
    for (r = 0; r < size; ++r) {
        uint32_t c;
        for (c = 0; c < size; ++c) {
            const uint16_t pixels[] = { ref_samples[topOffset + c], below_pred, ref_samples[leftOffset + r], right_pred };
            const uint8_t weights[] = { sm_weights_h[r], (uint8_t)(scale - sm_weights_h[r]),
                sm_weights_w[c], (uint8_t)(scale - sm_weights_w[c]) };
            uint32_t this_pred = 0;
            int32_t i;
            assert(scale >= sm_weights_h[r] && scale >= sm_weights_w[c]);
            for (i = 0; i < 4; ++i)
                this_pred += weights[i] * pixels[i];
            dst[c] = (uint16_t)divide_round(this_pred, log2_scale);
        }
        dst += rowStride * prediction_buffer_stride;
    }

    return;
}

void v_predictor_16bit(uint16_t *dst, const uint32_t stride, int32_t bw, int32_t bh,
    const uint16_t *ref_samples) {
    int32_t r;
    int32_t c;

    for (r = 0; r < bh; r++) {
        for (c = 0; c < bh; c++) {
            dst[c + r * stride] = ref_samples[bw + bh + 1 + c];
            //EB_MEMSET(dst, ref_samples[bw + bh + 1], bw);
            //dst += stride;
        }
    }
}

void h_predictor_16bit(uint16_t *dst, const uint32_t stride, int32_t bw, int32_t bh,
    const uint16_t *ref_samples) {
    int32_t r;
    //(void)above;

    for (r = 0; r < bh; r++) {
        memset16bit(dst, ref_samples[r], bw);
        dst += stride;
    }

    return;
}

void intra_mode_angular_av1_z1_16bit(
    const uint32_t   size,                    //input parameter, denotes the size of the current PU
    uint16_t         *ref_samples,             //input parameter, pointer to the reference samples
    uint16_t         *dst,                    //output parameter, pointer to the prediction
    const uint32_t   prediction_buffer_stride,  //input parameter, denotes the stride for the prediction ptr
    const EbBool  skip,
    uint16_t          dx,                     //output parameter, pointer to the prediction
    uint16_t          dy,                      //output parameter, pointer to the prediction
    uint16_t          bd)

{
    (void)dy;
    //    uint32_t row_index, colIndex;
    uint32_t toAboveOffset = (size << 1) + 1;
    uint32_t rowStride = skip ? 2 : 1;

    uint32_t r, c;
    int32_t x, base, shift, val;
    const int32_t max_base_x = ((size + size) - 1);
    const int32_t frac_bits = 6;
    const int32_t base_inc = 1;

    x = dx;

    for (r = 0; r < size; ++r, dst += (rowStride* prediction_buffer_stride), x += dx) {
        base = x >> frac_bits;
        shift = ((x) & 0x3F) >> 1;

        if (base >= max_base_x) {
            for (uint32_t i = r; i < size; ++i) {
                memset16bit(dst, ref_samples[toAboveOffset + max_base_x], size);
                dst += rowStride * prediction_buffer_stride;
            }
            return;
        }

        for (c = 0; c < size; ++c, base += base_inc) {
            if (base < max_base_x) {
                val = ref_samples[toAboveOffset + base] * (32 - shift) + ref_samples[toAboveOffset + base + 1] * shift;
                val = ROUND_POWER_OF_TWO(val, 5);
                dst[c] = clip_pixel_highbd(val, bd);
            }
            else
                dst[c] = ref_samples[toAboveOffset + max_base_x];
        }
    }

    return;
}
void intra_mode_angular_av1_z2_16bit(
    const uint32_t   size,                       //input parameter, denotes the size of the current PU
    uint16_t         *ref_samples,                 //input parameter, pointer to the reference samples
    uint16_t         *dst,              //output parameter, pointer to the prediction
    const uint32_t   prediction_buffer_stride,     //input parameter, denotes the stride for the prediction ptr
    const EbBool  skip,
    uint16_t          dx,              //output parameter, pointer to the prediction
    uint16_t          dy,              //output parameter, pointer to the prediction
    uint16_t          bd)
{
    //    uint32_t row_index, colIndex;
    uint32_t toAboveOffset = (size << 1) + 1;
    uint32_t toAboveLeftOffset = (size << 1);
    uint32_t shiftToAboveLeft = 0;
    uint32_t toLeftOffset = 0;
    uint32_t rowStride = skip ? 2 : 1;

    uint32_t r, c;
    int32_t x, y, shift, val, base;

    const int32_t min_base_x = -(1);
    const int32_t frac_bits_x = 6;
    const int32_t frac_bits_y = 6;

    for (r = 0; r < size; ++r) {
        for (c = 0; c < size; ++c) {
            y = r + 1;
            x = (c << 6) - y * dx;
            base = x >> frac_bits_x;
            if (base >= min_base_x) {
                shift = ((x) & 0x3F) >> 1;
                shiftToAboveLeft = (base <= -1) ? -1 - base : 0;
                val = ref_samples[toAboveOffset + base + shiftToAboveLeft] * (32 - shift) + ref_samples[toAboveOffset + base + 1] * shift;
                val = ROUND_POWER_OF_TWO(val, 5);
            }
            else {
                x = c + 1;
                y = (r << 6) - x * dy;
                base = y >> frac_bits_y;
                shiftToAboveLeft = (base <= -1) ? toAboveLeftOffset - base : 0;
                shift = ((y) & 0x3F) >> 1;
                val = ref_samples[toLeftOffset + base + shiftToAboveLeft] * (32 - shift) + ref_samples[toLeftOffset + base + 1] * shift;
                val = ROUND_POWER_OF_TWO(val, 5);
            }
            dst[c] = clip_pixel_highbd(val, bd);
        }
        dst += (rowStride* prediction_buffer_stride);
    }

    return;
}
void intra_mode_angular_av1_z3_16bit(
    const uint32_t   size,                        //input parameter, denotes the size of the current PU
    uint16_t         *ref_samples,                  //input parameter, pointer to the reference samples
    uint16_t         *dst,                         //output parameter, pointer to the prediction
    const uint32_t   prediction_buffer_stride,      //input parameter, denotes the stride for the prediction ptr
    const EbBool  skip,
    uint16_t          dx,                          //output parameter, pointer to the prediction
    uint16_t          dy,                          //output parameter, pointer to the prediction
    uint16_t          bd)
{
    //    uint32_t row_index, colIndex;
    //    uint32_t toAboveOffset = (size << 1) + 1;
    uint32_t toLeftOffset = 0;
    uint32_t rowStride = skip ? 2 : 1;

    uint32_t r, c;
    int32_t  y, base, shift, val;

    (void)dx;
    assert(dx == 1);
    assert(dy > 0);

    const int32_t max_base_y = (size + size - 1);
    const int32_t frac_bits = 6;
    const int32_t base_inc = 1;
    y = dy;
    for (c = 0; c < size; ++c, y += dy) {
        base = y >> frac_bits;
        shift = ((y) & 0x3F) >> 1;

        for (r = 0; r < size; ++r, base += base_inc) {
            if (base < max_base_y) {
                val = ref_samples[toLeftOffset + base] * (32 - shift) + ref_samples[toLeftOffset + base + 1] * shift;
                val = ROUND_POWER_OF_TWO(val, 5);
                dst[r * (rowStride* prediction_buffer_stride) + c] = clip_pixel_highbd(val, bd);
            }
            else {
                for (; r < size; ++r) dst[r * (rowStride* prediction_buffer_stride) + c] = ref_samples[toLeftOffset + max_base_y];
                break;
            }
        }
    }

    return;
}

void intra_open_loop_reference_samples_dctor(EbPtr p)
{
    IntraReferenceSamplesOpenLoop *obj = (IntraReferenceSamplesOpenLoop*)p;
    obj->y_intra_reference_array_reverse--;
    EB_FREE(obj->y_intra_reference_array_reverse);
    EB_FREE(obj->y_intra_reference_array);
}
/**********************************************
 * Intra Reference Samples Ctor
 **********************************************/
EbErrorType intra_open_loop_reference_samples_ctor(
    IntraReferenceSamplesOpenLoop *context_ptr)
{
    context_ptr->dctor = intra_open_loop_reference_samples_dctor;
    EB_MALLOC(context_ptr->y_intra_reference_array, (4 * BLOCK_SIZE_64 + 1));

    EB_MALLOC(context_ptr->y_intra_reference_array_reverse, (4 * BLOCK_SIZE_64 + 2));

    context_ptr->y_intra_reference_array_reverse++;

    return EB_ErrorNone;
}

void cfl_luma_subsampling_420_lbd_c(
    const uint8_t *input,
    int32_t input_stride, int16_t *output_q3,
    int32_t width, int32_t height)
{
    for (int32_t j = 0; j < height; j += 2) {
        for (int32_t i = 0; i < width; i += 2) {
            const int32_t bot = i + input_stride;
            output_q3[i >> 1] =
                (input[i] + input[i + 1] + input[bot] + input[bot + 1]) << 1;
        }
        input += input_stride << 1;
        output_q3 += CFL_BUF_LINE;
    }
}
void cfl_luma_subsampling_420_hbd_c(
    const uint16_t *input,
    int32_t input_stride, int16_t *output_q3,
    int32_t width, int32_t height)
{
    for (int32_t j = 0; j < height; j += 2) {
        for (int32_t i = 0; i < width; i += 2) {
            const int32_t bot = i + input_stride;
            output_q3[i >> 1] =
                (input[i] + input[i + 1] + input[bot] + input[bot + 1]) << 1;
        }
        input += input_stride << 1;
        output_q3 += CFL_BUF_LINE;
    }
}
void eb_subtract_average_c(
    int16_t *pred_buf_q3,
    int32_t width,
    int32_t height,
    int32_t round_offset,
    int32_t num_pel_log2) {
    int32_t sum_q3 = 0;
    int16_t *pred_buf = pred_buf_q3;
    for (int32_t j = 0; j < height; j++) {
        // assert(pred_buf_q3 + tx_width <= cfl->pred_buf_q3 + CFL_BUF_SQUARE);
        for (int32_t i = 0; i < width; i++)
            sum_q3 += pred_buf[i];
        pred_buf += CFL_BUF_LINE;
    }
    const int32_t avg_q3 = (sum_q3 + round_offset) >> num_pel_log2;
    // Loss is never more than 1/2 (in Q3)
    // assert(abs((avg_q3 * (1 << num_pel_log2)) - sum_q3) <= 1 << num_pel_log2 >>
    //       1);
    for (int32_t j = 0; j < height; j++) {
        for (int32_t i = 0; i < width; i++)
            pred_buf_q3[i] -= (int16_t)(avg_q3);
        pred_buf_q3 += CFL_BUF_LINE;
    }
}

CFL_SUB_AVG_FN(c)

void eb_cfl_predict_lbd_c(
    const int16_t *pred_buf_q3,
    uint8_t *pred,// AMIR ADDED
    int32_t pred_stride,
    uint8_t *dst,// AMIR changed to 8 bit
    int32_t dst_stride,
    int32_t alpha_q3,
    int32_t bit_depth,
    int32_t width,
    int32_t height) {
    for (int32_t j = 0; j < height; j++) {
        for (int32_t i = 0; i < width; i++) {
            dst[i] = (uint8_t)clip_pixel_highbd(
                get_scaled_luma_q0(alpha_q3, pred_buf_q3[i]) + (int16_t)pred[i], bit_depth);
        }
        dst += dst_stride;
        pred += pred_stride;
        pred_buf_q3 += CFL_BUF_LINE;
    }
}
void eb_cfl_predict_hbd_c(
    const int16_t *pred_buf_q3,
    uint16_t *pred,// AMIR ADDED
    int32_t pred_stride,
    uint16_t *dst,// AMIR changed to 8 bit
    int32_t dst_stride,
    int32_t alpha_q3,
    int32_t bit_depth,
    int32_t width,
    int32_t height) {
    for (int32_t j = 0; j < height; j++) {
        for (int32_t i = 0; i < width; i++) {
            dst[i] = clip_pixel_highbd(
                get_scaled_luma_q0(alpha_q3, pred_buf_q3[i]) + (int16_t)pred[i], bit_depth);
        }
        dst += dst_stride;
        pred += pred_stride;
        pred_buf_q3 += CFL_BUF_LINE;
    }
}

const uint8_t extend_modes[INTRA_MODES] = {
    NEED_ABOVE | NEED_LEFT,                   // DC
    NEED_ABOVE,                               // V
    NEED_LEFT,                                // H
    NEED_ABOVE | NEED_ABOVERIGHT,             // D45
    NEED_LEFT | NEED_ABOVE | NEED_ABOVELEFT,  // D135
    NEED_LEFT | NEED_ABOVE | NEED_ABOVELEFT,  // D113
    NEED_LEFT | NEED_ABOVE | NEED_ABOVELEFT,  // D157
    NEED_LEFT | NEED_BOTTOMLEFT,              // D203
    NEED_ABOVE | NEED_ABOVERIGHT,             // D67
    NEED_LEFT | NEED_ABOVE,                   // SMOOTH
    NEED_LEFT | NEED_ABOVE,                   // SMOOTH_V
    NEED_LEFT | NEED_ABOVE,                   // SMOOTH_H
    NEED_LEFT | NEED_ABOVE | NEED_ABOVELEFT,  // PAETH
};

// Tables to store if the top-right reference pixels are available. The flags
// are represented with bits, packed into 8-bit integers. E.g., for the 32x32
// blocks in a 128x128 superblock, the index of the "o" block is 10 (in raster
// order), so its flag is stored at the 3rd bit of the 2nd entry in the table,
// i.e. (table[10 / 8] >> (10 % 8)) & 1.
//       . . . .
//       . . . .
//       . . o .
//       . . . .
static uint8_t has_tr_4x4[128] = {
    255, 255, 255, 255, 85, 85, 85, 85, 119, 119, 119, 119, 85, 85, 85, 85,
    127, 127, 127, 127, 85, 85, 85, 85, 119, 119, 119, 119, 85, 85, 85, 85,
    255, 127, 255, 127, 85, 85, 85, 85, 119, 119, 119, 119, 85, 85, 85, 85,
    127, 127, 127, 127, 85, 85, 85, 85, 119, 119, 119, 119, 85, 85, 85, 85,
    255, 255, 255, 127, 85, 85, 85, 85, 119, 119, 119, 119, 85, 85, 85, 85,
    127, 127, 127, 127, 85, 85, 85, 85, 119, 119, 119, 119, 85, 85, 85, 85,
    255, 127, 255, 127, 85, 85, 85, 85, 119, 119, 119, 119, 85, 85, 85, 85,
    127, 127, 127, 127, 85, 85, 85, 85, 119, 119, 119, 119, 85, 85, 85, 85,
};
static uint8_t has_tr_4x8[64] = {
    255, 255, 255, 255, 119, 119, 119, 119, 127, 127, 127, 127, 119,
    119, 119, 119, 255, 127, 255, 127, 119, 119, 119, 119, 127, 127,
    127, 127, 119, 119, 119, 119, 255, 255, 255, 127, 119, 119, 119,
    119, 127, 127, 127, 127, 119, 119, 119, 119, 255, 127, 255, 127,
    119, 119, 119, 119, 127, 127, 127, 127, 119, 119, 119, 119,
};
static uint8_t has_tr_8x4[64] = {
    255, 255, 0, 0, 85, 85, 0, 0, 119, 119, 0, 0, 85, 85, 0, 0,
    127, 127, 0, 0, 85, 85, 0, 0, 119, 119, 0, 0, 85, 85, 0, 0,
    255, 127, 0, 0, 85, 85, 0, 0, 119, 119, 0, 0, 85, 85, 0, 0,
    127, 127, 0, 0, 85, 85, 0, 0, 119, 119, 0, 0, 85, 85, 0, 0,
};
static uint8_t has_tr_8x8[32] = {
    255, 255, 85, 85, 119, 119, 85, 85, 127, 127, 85, 85, 119, 119, 85, 85,
    255, 127, 85, 85, 119, 119, 85, 85, 127, 127, 85, 85, 119, 119, 85, 85,
};
static uint8_t has_tr_8x16[16] = {
    255, 255, 119, 119, 127, 127, 119, 119,
    255, 127, 119, 119, 127, 127, 119, 119,
};
static uint8_t has_tr_16x8[16] = {
    255, 0, 85, 0, 119, 0, 85, 0, 127, 0, 85, 0, 119, 0, 85, 0,
};
static uint8_t has_tr_16x16[8] = {
    255, 85, 119, 85, 127, 85, 119, 85,
};
static uint8_t has_tr_16x32[4] = { 255, 119, 127, 119 };
static uint8_t has_tr_32x16[4] = { 15, 5, 7, 5 };
static uint8_t has_tr_32x32[2] = { 95, 87 };
static uint8_t has_tr_32x64[1] = { 127 };
static uint8_t has_tr_64x32[1] = { 19 };
static uint8_t has_tr_64x64[1] = { 7 };
static uint8_t has_tr_64x128[1] = { 3 };
static uint8_t has_tr_128x64[1] = { 1 };
static uint8_t has_tr_128x128[1] = { 1 };
static uint8_t has_tr_4x16[32] = {
    255, 255, 255, 255, 127, 127, 127, 127, 255, 127, 255,
    127, 127, 127, 127, 127, 255, 255, 255, 127, 127, 127,
    127, 127, 255, 127, 255, 127, 127, 127, 127, 127,
};
static uint8_t has_tr_16x4[32] = {
    255, 0, 0, 0, 85, 0, 0, 0, 119, 0, 0, 0, 85, 0, 0, 0,
    127, 0, 0, 0, 85, 0, 0, 0, 119, 0, 0, 0, 85, 0, 0, 0,
};
static uint8_t has_tr_8x32[8] = {
    255, 255, 127, 127, 255, 127, 127, 127,
};
static uint8_t has_tr_32x8[8] = {
    15, 0, 5, 0, 7, 0, 5, 0,
};
static uint8_t has_tr_16x64[2] = { 255, 127 };
static uint8_t has_tr_64x16[2] = { 3, 1 };

const uint8_t *const has_tr_tables[BlockSizeS_ALL] = {
    // 4X4
    has_tr_4x4,
    // 4X8,       8X4,            8X8
    has_tr_4x8, has_tr_8x4, has_tr_8x8,
    // 8X16,      16X8,           16X16
    has_tr_8x16, has_tr_16x8, has_tr_16x16,
    // 16X32,     32X16,          32X32
    has_tr_16x32, has_tr_32x16, has_tr_32x32,
    // 32X64,     64X32,          64X64
    has_tr_32x64, has_tr_64x32, has_tr_64x64,
    // 64x128,    128x64,         128x128
    has_tr_64x128, has_tr_128x64, has_tr_128x128,
    // 4x16,      16x4,            8x32
    has_tr_4x16, has_tr_16x4, has_tr_8x32,
    // 32x8,      16x64,           64x16
    has_tr_32x8, has_tr_16x64, has_tr_64x16
};

static uint8_t has_tr_vert_8x8[32] = {
    255, 255, 0, 0, 119, 119, 0, 0, 127, 127, 0, 0, 119, 119, 0, 0,
    255, 127, 0, 0, 119, 119, 0, 0, 127, 127, 0, 0, 119, 119, 0, 0,
};
static uint8_t has_tr_vert_16x16[8] = {
    255, 0, 119, 0, 127, 0, 119, 0,
};
static uint8_t has_tr_vert_32x32[2] = { 15, 7 };
static uint8_t has_tr_vert_64x64[1] = { 3 };

// The _vert_* tables are like the ordinary tables above, but describe the
// order we visit square blocks when doing a PARTITION_VERT_A or
// PARTITION_VERT_B. This is the same order as normal except for on the last
// split where we go vertically (TL, BL, TR, BR). We treat the rectangular block
// as a pair of squares, which means that these tables work correctly for both
// mixed vertical partition types.
//
// There are tables for each of the square sizes. Vertical rectangles (like
// BLOCK_16X32) use their respective "non-vert" table
const uint8_t *const has_tr_vert_tables[BlockSizeS] = {
    // 4X4
    NULL,
    // 4X8,      8X4,         8X8
    has_tr_4x8, NULL, has_tr_vert_8x8,
    // 8X16,     16X8,        16X16
    has_tr_8x16, NULL, has_tr_vert_16x16,
    // 16X32,    32X16,       32X32
    has_tr_16x32, NULL, has_tr_vert_32x32,
    // 32X64,    64X32,       64X64
    has_tr_32x64, NULL, has_tr_vert_64x64,
    // 64x128,   128x64,      128x128
    has_tr_64x128, NULL, has_tr_128x128
};

static const uint8_t *get_has_tr_table(PartitionType partition,
    BlockSize bsize) {
    const uint8_t *ret = NULL;
    // If this is a mixed vertical partition, look up bsize in orders_vert.
    if (partition == PARTITION_VERT_A || partition == PARTITION_VERT_B) {
        assert(bsize < BlockSizeS);
        ret = has_tr_vert_tables[bsize];
    }
    else
        ret = has_tr_tables[bsize];
    assert(ret);
    return ret;
}

int32_t intra_has_top_right(BlockSize   sb_size, BlockSize bsize, int32_t mi_row,
    int32_t mi_col, int32_t top_available, int32_t right_available,
    PartitionType partition, TxSize txsz, int32_t row_off,
    int32_t col_off, int32_t ss_x, int32_t ss_y) {
    if (!top_available || !right_available) return 0;

    const int32_t bw_unit = block_size_wide[bsize] >> tx_size_wide_log2[0];
    const int32_t plane_bw_unit = AOMMAX(bw_unit >> ss_x, 1);
    const int32_t top_right_count_unit = tx_size_wide_unit[txsz];

    if (row_off > 0) {  // Just need to check if enough pixels on the right.
        if (block_size_wide[bsize] > block_size_wide[BLOCK_64X64]) {
            // Special case: For 128x128 blocks, the transform unit whose
            // top-right corner is at the center of the block does in fact have
            // pixels available at its top-right corner.
            if (row_off == mi_size_high[BLOCK_64X64] >> ss_y &&
                col_off + top_right_count_unit == mi_size_wide[BLOCK_64X64] >> ss_x) {
                return 1;
            }
            const int32_t plane_bw_unit_64 = mi_size_wide[BLOCK_64X64] >> ss_x;
            const int32_t col_off_64 = col_off % plane_bw_unit_64;
            return col_off_64 + top_right_count_unit < plane_bw_unit_64;
        }
        return col_off + top_right_count_unit < plane_bw_unit;
    }
    else {
        // All top-right pixels are in the block above, which is already available.
        if (col_off + top_right_count_unit < plane_bw_unit) return 1;

        const int32_t bw_in_mi_log2 = mi_size_wide_log2[bsize];
        const int32_t bh_in_mi_log2 = mi_size_high_log2[bsize];
        const int32_t sb_mi_size = mi_size_high[sb_size];
        const int32_t blk_row_in_sb = (mi_row & (sb_mi_size - 1)) >> bh_in_mi_log2;
        const int32_t blk_col_in_sb = (mi_col & (sb_mi_size - 1)) >> bw_in_mi_log2;

        // Top row of superblock: so top-right pixels are in the top and/or
        // top-right superblocks, both of which are already available.
        if (blk_row_in_sb == 0) return 1;

        // Rightmost column of superblock (and not the top row): so top-right pixels
        // fall in the right superblock, which is not available yet.
        if (((blk_col_in_sb + 1) << bw_in_mi_log2) >= sb_mi_size)
            return 0;
        // General case (neither top row nor rightmost column): check if the
        // top-right block is coded before the current block.
        const int32_t this_blk_index =
            ((blk_row_in_sb + 0) << (MAX_MIB_SIZE_LOG2 - bw_in_mi_log2)) +
            blk_col_in_sb + 0;
        const int32_t idx1 = this_blk_index / 8;
        const int32_t idx2 = this_blk_index % 8;
        const uint8_t *has_tr_table = get_has_tr_table(partition, bsize);
        return (has_tr_table[idx1] >> idx2) & 1;
    }
}

// Similar to the has_tr_* tables, but store if the bottom-left reference
// pixels are available.
static uint8_t has_bl_4x4[128] = {
    84, 85, 85, 85, 16, 17, 17, 17, 84, 85, 85, 85, 0, 1, 1, 1, 84, 85, 85,
    85, 16, 17, 17, 17, 84, 85, 85, 85, 0, 0, 1, 0, 84, 85, 85, 85, 16, 17,
    17, 17, 84, 85, 85, 85, 0, 1, 1, 1, 84, 85, 85, 85, 16, 17, 17, 17, 84,
    85, 85, 85, 0, 0, 0, 0, 84, 85, 85, 85, 16, 17, 17, 17, 84, 85, 85, 85,
    0, 1, 1, 1, 84, 85, 85, 85, 16, 17, 17, 17, 84, 85, 85, 85, 0, 0, 1,
    0, 84, 85, 85, 85, 16, 17, 17, 17, 84, 85, 85, 85, 0, 1, 1, 1, 84, 85,
    85, 85, 16, 17, 17, 17, 84, 85, 85, 85, 0, 0, 0, 0,
};
static uint8_t has_bl_4x8[64] = {
    16, 17, 17, 17, 0, 1, 1, 1, 16, 17, 17, 17, 0, 0, 1, 0,
    16, 17, 17, 17, 0, 1, 1, 1, 16, 17, 17, 17, 0, 0, 0, 0,
    16, 17, 17, 17, 0, 1, 1, 1, 16, 17, 17, 17, 0, 0, 1, 0,
    16, 17, 17, 17, 0, 1, 1, 1, 16, 17, 17, 17, 0, 0, 0, 0,
};
static uint8_t has_bl_8x4[64] = {
    254, 255, 84, 85, 254, 255, 16, 17, 254, 255, 84, 85, 254, 255, 0, 1,
    254, 255, 84, 85, 254, 255, 16, 17, 254, 255, 84, 85, 254, 255, 0, 0,
    254, 255, 84, 85, 254, 255, 16, 17, 254, 255, 84, 85, 254, 255, 0, 1,
    254, 255, 84, 85, 254, 255, 16, 17, 254, 255, 84, 85, 254, 255, 0, 0,
};
static uint8_t has_bl_8x8[32] = {
    84, 85, 16, 17, 84, 85, 0, 1, 84, 85, 16, 17, 84, 85, 0, 0,
    84, 85, 16, 17, 84, 85, 0, 1, 84, 85, 16, 17, 84, 85, 0, 0,
};
static uint8_t has_bl_8x16[16] = {
    16, 17, 0, 1, 16, 17, 0, 0, 16, 17, 0, 1, 16, 17, 0, 0,
};
static uint8_t has_bl_16x8[16] = {
    254, 84, 254, 16, 254, 84, 254, 0, 254, 84, 254, 16, 254, 84, 254, 0,
};
static uint8_t has_bl_16x16[8] = {
    84, 16, 84, 0, 84, 16, 84, 0,
};
static uint8_t has_bl_16x32[4] = { 16, 0, 16, 0 };
static uint8_t has_bl_32x16[4] = { 78, 14, 78, 14 };
static uint8_t has_bl_32x32[2] = { 4, 4 };
static uint8_t has_bl_32x64[1] = { 0 };
static uint8_t has_bl_64x32[1] = { 34 };
static uint8_t has_bl_64x64[1] = { 0 };
static uint8_t has_bl_64x128[1] = { 0 };
static uint8_t has_bl_128x64[1] = { 0 };
static uint8_t has_bl_128x128[1] = { 0 };
static uint8_t has_bl_4x16[32] = {
    0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0,
    0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0,
};
static uint8_t has_bl_16x4[32] = {
    254, 254, 254, 84, 254, 254, 254, 16, 254, 254, 254, 84, 254, 254, 254, 0,
    254, 254, 254, 84, 254, 254, 254, 16, 254, 254, 254, 84, 254, 254, 254, 0,
};
static uint8_t has_bl_8x32[8] = {
    0, 1, 0, 0, 0, 1, 0, 0,
};
static uint8_t has_bl_32x8[8] = {
    238, 78, 238, 14, 238, 78, 238, 14,
};
static uint8_t has_bl_16x64[2] = { 0, 0 };
static uint8_t has_bl_64x16[2] = { 42, 42 };

const uint8_t *const has_bl_tables[BlockSizeS_ALL] = {
    // 4X4
    has_bl_4x4,
    // 4X8,         8X4,         8X8
    has_bl_4x8, has_bl_8x4, has_bl_8x8,
    // 8X16,        16X8,        16X16
    has_bl_8x16, has_bl_16x8, has_bl_16x16,
    // 16X32,       32X16,       32X32
    has_bl_16x32, has_bl_32x16, has_bl_32x32,
    // 32X64,       64X32,       64X64
    has_bl_32x64, has_bl_64x32, has_bl_64x64,
    // 64x128,      128x64,      128x128
    has_bl_64x128, has_bl_128x64, has_bl_128x128,
    // 4x16,        16x4,        8x32
    has_bl_4x16, has_bl_16x4, has_bl_8x32,
    // 32x8,        16x64,       64x16
    has_bl_32x8, has_bl_16x64, has_bl_64x16
};

static uint8_t has_bl_vert_8x8[32] = {
    254, 255, 16, 17, 254, 255, 0, 1, 254, 255, 16, 17, 254, 255, 0, 0,
    254, 255, 16, 17, 254, 255, 0, 1, 254, 255, 16, 17, 254, 255, 0, 0,
};
static uint8_t has_bl_vert_16x16[8] = {
    254, 16, 254, 0, 254, 16, 254, 0,
};
static uint8_t has_bl_vert_32x32[2] = { 14, 14 };
static uint8_t has_bl_vert_64x64[1] = { 2 };

// The _vert_* tables are like the ordinary tables above, but describe the
// order we visit square blocks when doing a PARTITION_VERT_A or
// PARTITION_VERT_B. This is the same order as normal except for on the last
// split where we go vertically (TL, BL, TR, BR). We treat the rectangular block
// as a pair of squares, which means that these tables work correctly for both
// mixed vertical partition types.
//
// There are tables for each of the square sizes. Vertical rectangles (like
// BLOCK_16X32) use their respective "non-vert" table
const uint8_t *const has_bl_vert_tables[BlockSizeS] = {
    // 4X4
    NULL,
    // 4X8,     8X4,         8X8
    has_bl_4x8, NULL, has_bl_vert_8x8,
    // 8X16,    16X8,        16X16
    has_bl_8x16, NULL, has_bl_vert_16x16,
    // 16X32,   32X16,       32X32
    has_bl_16x32, NULL, has_bl_vert_32x32,
    // 32X64,   64X32,       64X64
    has_bl_32x64, NULL, has_bl_vert_64x64,
    // 64x128,  128x64,      128x128
    has_bl_64x128, NULL, has_bl_128x128
};

static const uint8_t *get_has_bl_table(PartitionType partition,
    BlockSize bsize) {
    const uint8_t *ret = NULL;
    // If this is a mixed vertical partition, look up bsize in orders_vert.
    if (partition == PARTITION_VERT_A || partition == PARTITION_VERT_B) {
        assert(bsize < BlockSizeS);
        ret = has_bl_vert_tables[bsize];
    }
    else
        ret = has_bl_tables[bsize];
    assert(ret);
    return ret;
}

int32_t intra_has_bottom_left(BlockSize sb_size, BlockSize bsize, int32_t mi_row,
    int32_t mi_col, int32_t bottom_available, int32_t left_available,
    PartitionType partition, TxSize txsz, int32_t row_off,
    int32_t col_off, int32_t ss_x, int32_t ss_y) {
    if (!bottom_available || !left_available) return 0;

    // Special case for 128x* blocks, when col_off is half the block width.
    // This is needed because 128x* superblocks are divided into 64x* blocks in
    // raster order
    if (block_size_wide[bsize] > block_size_wide[BLOCK_64X64] && col_off > 0) {
        const int32_t plane_bw_unit_64 = mi_size_wide[BLOCK_64X64] >> ss_x;
        const int32_t col_off_64 = col_off % plane_bw_unit_64;
        if (col_off_64 == 0) {
            // We are at the left edge of top-right or bottom-right 64x* block.
            const int32_t plane_bh_unit_64 = mi_size_high[BLOCK_64X64] >> ss_y;
            const int32_t row_off_64 = row_off % plane_bh_unit_64;
            const int32_t plane_bh_unit =
                AOMMIN(mi_size_high[bsize] >> ss_y, plane_bh_unit_64);
            // Check if all bottom-left pixels are in the left 64x* block (which is
            // already coded).
            return row_off_64 + tx_size_high_unit[txsz] < plane_bh_unit;
        }
    }

    if (col_off > 0) {
        // Bottom-left pixels are in the bottom-left block, which is not available.
        return 0;
    }
    else {
        const int32_t bh_unit = block_size_high[bsize] >> tx_size_high_log2[0];
        const int32_t plane_bh_unit = AOMMAX(bh_unit >> ss_y, 1);
        const int32_t bottom_left_count_unit = tx_size_high_unit[txsz];

        // All bottom-left pixels are in the left block, which is already available.
        if (row_off + bottom_left_count_unit < plane_bh_unit) return 1;

        const int32_t bw_in_mi_log2 = mi_size_wide_log2[bsize];
        const int32_t bh_in_mi_log2 = mi_size_high_log2[bsize];
        const int32_t sb_mi_size = mi_size_high[sb_size];
        const int32_t blk_row_in_sb = (mi_row & (sb_mi_size - 1)) >> bh_in_mi_log2;
        const int32_t blk_col_in_sb = (mi_col & (sb_mi_size - 1)) >> bw_in_mi_log2;

        // Leftmost column of superblock: so bottom-left pixels maybe in the left
        // and/or bottom-left superblocks. But only the left superblock is
        // available, so check if all required pixels fall in that superblock.
        if (blk_col_in_sb == 0) {
            const int32_t blk_start_row_off = blk_row_in_sb
                << (bh_in_mi_log2 + MI_SIZE_LOG2 -
                    tx_size_wide_log2[0]) >>
                ss_y;
            const int32_t row_off_in_sb = blk_start_row_off + row_off;
            const int32_t sb_height_unit = sb_mi_size >> ss_y;
            return row_off_in_sb + bottom_left_count_unit < sb_height_unit;
        }

        // Bottom row of superblock (and not the leftmost column): so bottom-left
        // pixels fall in the bottom superblock, which is not available yet.
        if (((blk_row_in_sb + 1) << bh_in_mi_log2) >= sb_mi_size) return 0;

        // General case (neither leftmost column nor bottom row): check if the
        // bottom-left block is coded before the current block.
        const int32_t this_blk_index =
            ((blk_row_in_sb + 0) << (MAX_MIB_SIZE_LOG2 - bw_in_mi_log2)) +
            blk_col_in_sb + 0;
        const int32_t idx1 = this_blk_index / 8;
        const int32_t idx2 = this_blk_index % 8;
        const uint8_t *has_bl_table = get_has_bl_table(partition, bsize);
        return (has_bl_table[idx1] >> idx2) & 1;
    }
}

IntraPredFn pred[INTRA_MODES][TX_SIZES_ALL];
IntraPredFn dc_pred[2][2][TX_SIZES_ALL];

IntraHighPredFn pred_high[INTRA_MODES][TX_SIZES_ALL];
IntraHighPredFn dc_pred_high[2][2][TX_SIZES_ALL];

static INLINE void dc_128_predictor(uint8_t *dst, ptrdiff_t stride, int32_t bw,
    int32_t bh, const uint8_t *above,
    const uint8_t *left) {
    int32_t r;
    (void)above;
    (void)left;

    for (r = 0; r < bh; r++) {
        memset(dst, 128, bw);
        dst += stride;
    }
}

static INLINE void dc_left_predictor(uint8_t *dst, ptrdiff_t stride, int32_t bw,
    int32_t bh, const uint8_t *above,
    const uint8_t *left) {
    int32_t i, r, expected_dc, sum = 0;
    (void)above;

    for (i = 0; i < bh; i++) sum += left[i];
    expected_dc = (sum + (bh >> 1)) / bh;

    for (r = 0; r < bh; r++) {
        memset(dst, expected_dc, bw);
        dst += stride;
    }
}
static INLINE void dc_top_predictor(uint8_t *dst, ptrdiff_t stride, int32_t bw,
    int32_t bh, const uint8_t *above,
    const uint8_t *left) {
    int32_t i, r, expected_dc, sum = 0;
    (void)left;

    for (i = 0; i < bw; i++) sum += above[i];
    expected_dc = (sum + (bw >> 1)) / bw;

    for (r = 0; r < bh; r++) {
        memset(dst, expected_dc, bw);
        dst += stride;
    }
}
static INLINE void dc_predictor(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh,
    const uint8_t *above, const uint8_t *left) {
    int32_t i, r, expected_dc, sum = 0;
    const int32_t count = bw + bh;

    for (i = 0; i < bw; i++)
        sum += above[i];
    for (i = 0; i < bh; i++)
        sum += left[i];
    expected_dc = (sum + (count >> 1)) / count;

    for (r = 0; r < bh; r++) {
        memset(dst, expected_dc, bw);
        dst += stride;
    }
}
static INLINE void v_predictor(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh,
    const uint8_t *above, const uint8_t *left) {
    int32_t r;
    (void)left;

    for (r = 0; r < bh; r++) {
        memcpy(dst, above, bw);
        dst += stride;
    }
}

static INLINE void h_predictor(uint8_t *dst, ptrdiff_t stride, int32_t bw, int32_t bh,
    const uint8_t *above, const uint8_t *left) {
    int32_t r;
    (void)above;

    for (r = 0; r < bh; r++) {
        memset(dst, left[r], bw);
        dst += stride;
    }
}

static INLINE void smooth_predictor(uint8_t *dst, ptrdiff_t stride, int32_t bw,
    int32_t bh, const uint8_t *above,
    const uint8_t *left) {
    const uint8_t below_pred = left[bh - 1];   // estimated by bottom-left pixel
    const uint8_t right_pred = above[bw - 1];  // estimated by top-right pixel
    const uint8_t *const sm_weights_w = sm_weight_arrays + bw;
    const uint8_t *const sm_weights_h = sm_weight_arrays + bh;
    // scale = 2 * 2^sm_weight_log2_scale
    const int32_t log2_scale = 1 + sm_weight_log2_scale;
    const uint16_t scale = (1 << sm_weight_log2_scale);
    sm_weights_sanity_checks(sm_weights_w, sm_weights_h, scale,
        log2_scale + sizeof(*dst));
    int32_t r;
    for (r = 0; r < bh; ++r) {
        int32_t c;
        for (c = 0; c < bw; ++c) {
            const uint8_t pixels[] = { above[c], below_pred, left[r], right_pred };
            const uint8_t weights[] = { sm_weights_h[r], (uint8_t)(scale - sm_weights_h[r]),
                sm_weights_w[c], (uint8_t)(scale - sm_weights_w[c]) };
            uint32_t this_pred = 0;
            int32_t i;
            assert(scale >= sm_weights_h[r] && scale >= sm_weights_w[c]);
            for (i = 0; i < 4; ++i)
                this_pred += weights[i] * pixels[i];
            dst[c] = (uint8_t)divide_round(this_pred, log2_scale);
        }
        dst += stride;
    }
}

static INLINE void smooth_v_predictor(uint8_t *dst, ptrdiff_t stride, int32_t bw,
    int32_t bh, const uint8_t *above,
    const uint8_t *left) {
    const uint8_t below_pred = left[bh - 1];  // estimated by bottom-left pixel
    const uint8_t *const sm_weights = sm_weight_arrays + bh;
    // scale = 2^sm_weight_log2_scale
    const int32_t log2_scale = sm_weight_log2_scale;
    const uint16_t scale = (1 << sm_weight_log2_scale);
    sm_weights_sanity_checks(sm_weights, sm_weights, scale,
        log2_scale + sizeof(*dst));

    int32_t r;
    for (r = 0; r < bh; r++) {
        int32_t c;
        for (c = 0; c < bw; ++c) {
            const uint8_t pixels[] = { above[c], below_pred };
            const uint8_t weights[] = { sm_weights[r], (uint8_t)(scale - sm_weights[r]) };
            uint32_t this_pred = 0;
            assert(scale >= sm_weights[r]);
            int32_t i;
            for (i = 0; i < 2; ++i)
                this_pred += weights[i] * pixels[i];
            dst[c] = (uint8_t)divide_round(this_pred, log2_scale);
        }
        dst += stride;
    }
}

static INLINE void smooth_h_predictor(uint8_t *dst, ptrdiff_t stride, int32_t bw,
    int32_t bh, const uint8_t *above,
    const uint8_t *left) {
    const uint8_t right_pred = above[bw - 1];  // estimated by top-right pixel
    const uint8_t *const sm_weights = sm_weight_arrays + bw;
    // scale = 2^sm_weight_log2_scale
    const int32_t log2_scale = sm_weight_log2_scale;
    const uint16_t scale = (1 << sm_weight_log2_scale);
    sm_weights_sanity_checks(sm_weights, sm_weights, scale,
        log2_scale + sizeof(*dst));

    int32_t r;
    for (r = 0; r < bh; r++) {
        int32_t c;
        for (c = 0; c < bw; ++c) {
            const uint8_t pixels[] = { left[r], right_pred };
            const uint8_t weights[] = { sm_weights[c], (uint8_t)(scale - sm_weights[c]) };
            uint32_t this_pred = 0;
            assert(scale >= sm_weights[c]);
            int32_t i;
            for (i = 0; i < 2; ++i)
                this_pred += weights[i] * pixels[i];
            dst[c] = (uint8_t)divide_round(this_pred, log2_scale);
        }
        dst += stride;
    }
}
#undef DC_MULTIPLIER_1X2
#undef DC_MULTIPLIER_1X4

static INLINE void highbd_v_predictor(uint16_t *dst, ptrdiff_t stride, int32_t bw,
    int32_t bh, const uint16_t *above,
    const uint16_t *left, int32_t bd) {
    int32_t r;
    (void)left;
    (void)bd;
    for (r = 0; r < bh; r++) {
        memcpy(dst, above, bw * sizeof(uint16_t));
        dst += stride;
    }
}

static INLINE void highbd_h_predictor(uint16_t *dst, ptrdiff_t stride, int32_t bw,
    int32_t bh, const uint16_t *above,
    const uint16_t *left, int32_t bd) {
    int32_t r;
    (void)above;
    (void)bd;
    for (r = 0; r < bh; r++) {
        eb_aom_memset16(dst, left[r], bw);
        dst += stride;
    }
}
static INLINE int abs_diff(int a, int b) { return (a > b) ? a - b : b - a; }
static INLINE uint16_t paeth_predictor_single(uint16_t left, uint16_t top,
                                              uint16_t top_left) {
  const int base = top + left - top_left;
  const int p_left = abs_diff(base, left);
  const int p_top = abs_diff(base, top);
  const int p_top_left = abs_diff(base, top_left);

  // Return nearest to base of left, top and top_left.
  return (p_left <= p_top && p_left <= p_top_left)
             ? left
             : (p_top <= p_top_left) ? top : top_left;
}

static INLINE void paeth_predictor(uint8_t *dst, ptrdiff_t stride, int bw,
                                   int bh, const uint8_t *above,
                                   const uint8_t *left) {
  int r, c;
  const uint8_t ytop_left = above[-1];

  for (r = 0; r < bh; r++) {
    for (c = 0; c < bw; c++)
      dst[c] = (uint8_t)paeth_predictor_single(left[r], above[c], ytop_left);
    dst += stride;
  }
}

static INLINE void highbd_paeth_predictor(uint16_t *dst, ptrdiff_t stride,
                                          int bw, int bh, const uint16_t *above,
                                          const uint16_t *left, int bd) {
  int r, c;
  const uint16_t ytop_left = above[-1];
  (void)bd;

  for (r = 0; r < bh; r++) {
    for (c = 0; c < bw; c++)
      dst[c] = paeth_predictor_single(left[r], above[c], ytop_left);
    dst += stride;
  }
}
//static INLINE void highbd_paeth_predictor(uint16_t *dst, ptrdiff_t stride,
//    int32_t bw, int32_t bh, const uint16_t *above,
//    const uint16_t *left, int32_t bd) {
//    int32_t r, c;
//    const uint16_t ytop_left = above[-1];
//    (void)bd;
//
//    for (r = 0; r < bh; r++) {
//        for (c = 0; c < bw; c++)
//            dst[c] = paeth_predictor_single(left[r], above[c], ytop_left);
//        dst += stride;
//    }
//}

static INLINE void highbd_smooth_predictor(uint16_t *dst, ptrdiff_t stride,
    int32_t bw, int32_t bh,
    const uint16_t *above,
    const uint16_t *left, int32_t bd) {
    (void)bd;
    const uint16_t below_pred = left[bh - 1];   // estimated by bottom-left pixel
    const uint16_t right_pred = above[bw - 1];  // estimated by top-right pixel
    const uint8_t *const sm_weights_w = sm_weight_arrays + bw;
    const uint8_t *const sm_weights_h = sm_weight_arrays + bh;
    // scale = 2 * 2^sm_weight_log2_scale
    const int32_t log2_scale = 1 + sm_weight_log2_scale;
    const uint16_t scale = (1 << sm_weight_log2_scale);
    sm_weights_sanity_checks(sm_weights_w, sm_weights_h, scale,
        log2_scale + sizeof(*dst));
    int32_t r;
    for (r = 0; r < bh; ++r) {
        int32_t c;
        for (c = 0; c < bw; ++c) {
            const uint16_t pixels[] = { above[c], below_pred, left[r], right_pred };
            const uint8_t weights[] = { sm_weights_h[r], (uint8_t)(scale - sm_weights_h[r]),
                sm_weights_w[c], (uint8_t)(scale - sm_weights_w[c]) };
            uint32_t this_pred = 0;
            int32_t i;
            assert(scale >= sm_weights_h[r] && scale >= sm_weights_w[c]);
            for (i = 0; i < 4; ++i)
                this_pred += weights[i] * pixels[i];
            dst[c] = (uint16_t)divide_round(this_pred, log2_scale);
        }
        dst += stride;
    }
}

static INLINE void highbd_smooth_v_predictor(uint16_t *dst, ptrdiff_t stride,
    int32_t bw, int32_t bh,
    const uint16_t *above,
    const uint16_t *left, int32_t bd) {
    (void)bd;
    const uint16_t below_pred = left[bh - 1];  // estimated by bottom-left pixel
    const uint8_t *const sm_weights = sm_weight_arrays + bh;
    // scale = 2^sm_weight_log2_scale
    const int32_t log2_scale = sm_weight_log2_scale;
    const uint16_t scale = (1 << sm_weight_log2_scale);
    sm_weights_sanity_checks(sm_weights, sm_weights, scale,
        log2_scale + sizeof(*dst));

    int32_t r;
    for (r = 0; r < bh; r++) {
        int32_t c;
        for (c = 0; c < bw; ++c) {
            const uint16_t pixels[] = { above[c], below_pred };
            const uint8_t weights[] = { sm_weights[r], (uint8_t)(scale - sm_weights[r]) };
            uint32_t this_pred = 0;
            assert(scale >= sm_weights[r]);
            int32_t i;
            for (i = 0; i < 2; ++i)
                this_pred += weights[i] * pixels[i];
            dst[c] = (uint16_t)divide_round(this_pred, log2_scale);
        }
        dst += stride;
    }
}

static INLINE void highbd_smooth_h_predictor(uint16_t *dst, ptrdiff_t stride,
    int32_t bw, int32_t bh,
    const uint16_t *above,
    const uint16_t *left, int32_t bd) {
    (void)bd;
    const uint16_t right_pred = above[bw - 1];  // estimated by top-right pixel
    const uint8_t *const sm_weights = sm_weight_arrays + bw;
    // scale = 2^sm_weight_log2_scale
    const int32_t log2_scale = sm_weight_log2_scale;
    const uint16_t scale = (1 << sm_weight_log2_scale);
    sm_weights_sanity_checks(sm_weights, sm_weights, scale,
        log2_scale + sizeof(*dst));

    int32_t r;
    for (r = 0; r < bh; r++) {
        int32_t c;
        for (c = 0; c < bw; ++c) {
            const uint16_t pixels[] = { left[r], right_pred };
            const uint8_t weights[] = { sm_weights[c], (uint8_t)(scale - sm_weights[c]) };
            uint32_t this_pred = 0;
            assert(scale >= sm_weights[c]);
            int32_t i;
            for (i = 0; i < 2; ++i)
                this_pred += weights[i] * pixels[i];
            dst[c] = (uint16_t)divide_round(this_pred, log2_scale);
        }
        dst += stride;
    }
}

static INLINE void highbd_dc_128_predictor(uint16_t *dst, ptrdiff_t stride,
    int32_t bw, int32_t bh,
    const uint16_t *above,
    const uint16_t *left, int32_t bd) {
    int32_t r;
    (void)above;
    (void)left;

    for (r = 0; r < bh; r++) {
        eb_aom_memset16(dst, 128 << (bd - 8), bw);
        dst += stride;
    }
}

static INLINE void highbd_dc_left_predictor(uint16_t *dst, ptrdiff_t stride,
    int32_t bw, int32_t bh,
    const uint16_t *above,
    const uint16_t *left, int32_t bd) {
    int32_t i, r, expected_dc, sum = 0;
    (void)above;
    (void)bd;

    for (i = 0; i < bh; i++) sum += left[i];
    expected_dc = (sum + (bh >> 1)) / bh;

    for (r = 0; r < bh; r++) {
        eb_aom_memset16(dst, expected_dc, bw);
        dst += stride;
    }
}

static INLINE void highbd_dc_top_predictor(uint16_t *dst, ptrdiff_t stride,
    int32_t bw, int32_t bh,
    const uint16_t *above,
    const uint16_t *left, int32_t bd) {
    int32_t i, r, expected_dc, sum = 0;
    (void)left;
    (void)bd;

    for (i = 0; i < bw; i++) sum += above[i];
    expected_dc = (sum + (bw >> 1)) / bw;

    for (r = 0; r < bh; r++) {
        eb_aom_memset16(dst, expected_dc, bw);
        dst += stride;
    }
}

static INLINE void highbd_dc_predictor(uint16_t *dst, ptrdiff_t stride, int32_t bw,
    int32_t bh, const uint16_t *above,
    const uint16_t *left, int32_t bd) {
    int32_t i, r, expected_dc, sum = 0;
    const int32_t count = bw + bh;
    (void)bd;

    for (i = 0; i < bw; i++)
        sum += above[i];
    for (i = 0; i < bh; i++)
        sum += left[i];
    expected_dc = (sum + (count >> 1)) / count;

    for (r = 0; r < bh; r++) {
        eb_aom_memset16(dst, expected_dc, bw);
        dst += stride;
    }
}

//static INLINE void highbd_dc_predictor_rect(uint16_t *dst, ptrdiff_t stride,
//    int32_t bw, int32_t bh,
//    const uint16_t *above,
//    const uint16_t *left, int32_t bd,
//    int32_t shift1, uint32_t multiplier) {
//    int32_t sum = 0;
//    (void)bd;
//
//    for (int32_t i = 0; i < bw; i++) {
//        sum += above[i];
//    }
//    for (int32_t i = 0; i < bh; i++) {
//        sum += left[i];
//    }
//
//    const int32_t expected_dc = divide_using_multiply_shift(
//        sum + ((bw + bh) >> 1), shift1, multiplier, HIGHBD_DC_SHIFT2);
//    assert(expected_dc < (1 << bd));
//
//    for (int32_t r = 0; r < bh; r++) {
//        eb_aom_memset16(dst, expected_dc, bw);
//        dst += stride;
//    }
//}

//#undef HIGHBD_DC_SHIFT2
//
//void eb_aom_highbd_dc_predictor_4x8_c(uint16_t *dst, ptrdiff_t stride,
//    const uint16_t *above, const uint16_t *left,
//    int32_t bd) {
//    highbd_dc_predictor_rect(dst, stride, 4, 8, above, left, bd, 2,
//        HIGHBD_DC_MULTIPLIER_1X2);
//}
//
//void eb_aom_highbd_dc_predictor_8x4_c(uint16_t *dst, ptrdiff_t stride,
//    const uint16_t *above, const uint16_t *left,
//    int32_t bd) {
//    highbd_dc_predictor_rect(dst, stride, 8, 4, above, left, bd, 2,
//        HIGHBD_DC_MULTIPLIER_1X2);
//}
//
//void eb_aom_highbd_dc_predictor_4x16_c(uint16_t *dst, ptrdiff_t stride,
//    const uint16_t *above, const uint16_t *left,
//    int32_t bd) {
//    highbd_dc_predictor_rect(dst, stride, 4, 16, above, left, bd, 2,
//        HIGHBD_DC_MULTIPLIER_1X4);
//}
//
//void eb_aom_highbd_dc_predictor_16x4_c(uint16_t *dst, ptrdiff_t stride,
//    const uint16_t *above, const uint16_t *left,
//    int32_t bd) {
//    highbd_dc_predictor_rect(dst, stride, 16, 4, above, left, bd, 2,
//        HIGHBD_DC_MULTIPLIER_1X4);
//}
//
//void eb_aom_highbd_dc_predictor_8x16_c(uint16_t *dst, ptrdiff_t stride,
//    const uint16_t *above, const uint16_t *left,
//    int32_t bd) {
//    highbd_dc_predictor_rect(dst, stride, 8, 16, above, left, bd, 3,
//        HIGHBD_DC_MULTIPLIER_1X2);
//}
//
//void eb_aom_highbd_dc_predictor_16x8_c(uint16_t *dst, ptrdiff_t stride,
//    const uint16_t *above, const uint16_t *left,
//    int32_t bd) {
//    highbd_dc_predictor_rect(dst, stride, 16, 8, above, left, bd, 3,
//        HIGHBD_DC_MULTIPLIER_1X2);
//}
//
//void eb_aom_highbd_dc_predictor_8x32_c(uint16_t *dst, ptrdiff_t stride,
//    const uint16_t *above, const uint16_t *left,
//    int32_t bd) {
//    highbd_dc_predictor_rect(dst, stride, 8, 32, above, left, bd, 3,
//        HIGHBD_DC_MULTIPLIER_1X4);
//}
//
//void eb_aom_highbd_dc_predictor_32x8_c(uint16_t *dst, ptrdiff_t stride,
//    const uint16_t *above, const uint16_t *left,
//    int32_t bd) {
//    highbd_dc_predictor_rect(dst, stride, 32, 8, above, left, bd, 3,
//        HIGHBD_DC_MULTIPLIER_1X4);
//}
//
//void eb_aom_highbd_dc_predictor_16x32_c(uint16_t *dst, ptrdiff_t stride,
//    const uint16_t *above,
//    const uint16_t *left, int32_t bd) {
//    highbd_dc_predictor_rect(dst, stride, 16, 32, above, left, bd, 4,
//        HIGHBD_DC_MULTIPLIER_1X2);
//}
//
//void eb_aom_highbd_dc_predictor_32x16_c(uint16_t *dst, ptrdiff_t stride,
//    const uint16_t *above,
//    const uint16_t *left, int32_t bd) {
//    highbd_dc_predictor_rect(dst, stride, 32, 16, above, left, bd, 4,
//        HIGHBD_DC_MULTIPLIER_1X2);
//}
//
//void eb_aom_highbd_dc_predictor_16x64_c(uint16_t *dst, ptrdiff_t stride,
//    const uint16_t *above,
//    const uint16_t *left, int32_t bd) {
//    highbd_dc_predictor_rect(dst, stride, 16, 64, above, left, bd, 4,
//        HIGHBD_DC_MULTIPLIER_1X4);
//}
//
//void eb_aom_highbd_dc_predictor_64x16_c(uint16_t *dst, ptrdiff_t stride,
//    const uint16_t *above,
//    const uint16_t *left, int32_t bd) {
//    highbd_dc_predictor_rect(dst, stride, 64, 16, above, left, bd, 4,
//        HIGHBD_DC_MULTIPLIER_1X4);
//}
//
//void eb_aom_highbd_dc_predictor_32x64_c(uint16_t *dst, ptrdiff_t stride,
//    const uint16_t *above,
//    const uint16_t *left, int32_t bd) {
//    highbd_dc_predictor_rect(dst, stride, 32, 64, above, left, bd, 5,
//        HIGHBD_DC_MULTIPLIER_1X2);
//}
//
//void eb_aom_highbd_dc_predictor_64x32_c(uint16_t *dst, ptrdiff_t stride,
//    const uint16_t *above,
//    const uint16_t *left, int32_t bd) {
//    highbd_dc_predictor_rect(dst, stride, 64, 32, above, left, bd, 5,
//        HIGHBD_DC_MULTIPLIER_1X2);
//}
//
//#undef HIGHBD_DC_MULTIPLIER_1X2
//#undef HIGHBD_DC_MULTIPLIER_1X4

#define intra_pred_sized(type, width, height)                  \
  void eb_aom_##type##_predictor_##width##x##height##_c(          \
      uint8_t *dst, ptrdiff_t stride, const uint8_t *above,    \
      const uint8_t *left) {                                   \
    type##_predictor(dst, stride, width, height, above, left); \
  }

intra_pred_sized(dc, 2, 2)
intra_pred_sized(dc, 4, 4)
intra_pred_sized(dc, 8, 8)
intra_pred_sized(dc, 16, 16)
intra_pred_sized(dc, 32, 32)
intra_pred_sized(dc, 64, 64)
intra_pred_sized(dc, 4, 8)
intra_pred_sized(dc, 4, 16)
intra_pred_sized(dc, 8, 4)
intra_pred_sized(dc, 8, 16)
intra_pred_sized(dc, 8, 32)
intra_pred_sized(dc, 16, 4)
intra_pred_sized(dc, 16, 8)
intra_pred_sized(dc, 16, 32)
intra_pred_sized(dc, 16, 64)
intra_pred_sized(dc, 32, 8)
intra_pred_sized(dc, 32, 16)
intra_pred_sized(dc, 32, 64)
intra_pred_sized(dc, 64, 16)
intra_pred_sized(dc, 64, 32)

intra_pred_sized(dc_128, 2, 2)
intra_pred_sized(dc_128, 4, 4)
intra_pred_sized(dc_128, 8, 8)
intra_pred_sized(dc_128, 16, 16)
intra_pred_sized(dc_128, 32, 32)
intra_pred_sized(dc_128, 64, 64)
intra_pred_sized(dc_128, 4, 8)
intra_pred_sized(dc_128, 4, 16)
intra_pred_sized(dc_128, 8, 4)
intra_pred_sized(dc_128, 8, 16)
intra_pred_sized(dc_128, 8, 32)
intra_pred_sized(dc_128, 16, 4)
intra_pred_sized(dc_128, 16, 8)
intra_pred_sized(dc_128, 16, 32)
intra_pred_sized(dc_128, 16, 64)
intra_pred_sized(dc_128, 32, 8)
intra_pred_sized(dc_128, 32, 16)
intra_pred_sized(dc_128, 32, 64)
intra_pred_sized(dc_128, 64, 16)
intra_pred_sized(dc_128, 64, 32)

intra_pred_sized(dc_left, 2, 2)
intra_pred_sized(dc_left, 4, 4)
intra_pred_sized(dc_left, 8, 8)
intra_pred_sized(dc_left, 16, 16)
intra_pred_sized(dc_left, 32, 32)
intra_pred_sized(dc_left, 64, 64)
intra_pred_sized(dc_left, 4, 8)
intra_pred_sized(dc_left, 4, 16)
intra_pred_sized(dc_left, 8, 4)
intra_pred_sized(dc_left, 8, 16)
intra_pred_sized(dc_left, 8, 32)
intra_pred_sized(dc_left, 16, 4)
intra_pred_sized(dc_left, 16, 8)
intra_pred_sized(dc_left, 16, 32)
intra_pred_sized(dc_left, 16, 64)
intra_pred_sized(dc_left, 32, 8)
intra_pred_sized(dc_left, 32, 16)
intra_pred_sized(dc_left, 32, 64)
intra_pred_sized(dc_left, 64, 16)
intra_pred_sized(dc_left, 64, 32)

intra_pred_sized(dc_top, 2, 2)
intra_pred_sized(dc_top, 4, 4)
intra_pred_sized(dc_top, 8, 8)
intra_pred_sized(dc_top, 16, 16)
intra_pred_sized(dc_top, 32, 32)
intra_pred_sized(dc_top, 64, 64)
intra_pred_sized(dc_top, 4, 8)
intra_pred_sized(dc_top, 4, 16)
intra_pred_sized(dc_top, 8, 4)
intra_pred_sized(dc_top, 8, 16)
intra_pred_sized(dc_top, 8, 32)
intra_pred_sized(dc_top, 16, 4)
intra_pred_sized(dc_top, 16, 8)
intra_pred_sized(dc_top, 16, 32)
intra_pred_sized(dc_top, 16, 64)
intra_pred_sized(dc_top, 32, 8)
intra_pred_sized(dc_top, 32, 16)
intra_pred_sized(dc_top, 32, 64)
intra_pred_sized(dc_top, 64, 16)
intra_pred_sized(dc_top, 64, 32)
intra_pred_sized(v, 2, 2)
intra_pred_sized(v, 4, 4)
intra_pred_sized(v, 8, 8)
intra_pred_sized(v, 16, 16)
intra_pred_sized(v, 32, 32)
intra_pred_sized(v, 64, 64)
intra_pred_sized(v, 4, 8)
intra_pred_sized(v, 4, 16)
intra_pred_sized(v, 8, 4)
intra_pred_sized(v, 8, 16)
intra_pred_sized(v, 8, 32)
intra_pred_sized(v, 16, 4)
intra_pred_sized(v, 16, 8)
intra_pred_sized(v, 16, 32)
intra_pred_sized(v, 16, 64)
intra_pred_sized(v, 32, 8)
intra_pred_sized(v, 32, 16)
intra_pred_sized(v, 32, 64)
intra_pred_sized(v, 64, 16)
intra_pred_sized(v, 64, 32)
intra_pred_sized(h, 2, 2)
intra_pred_sized(h, 4, 4)
intra_pred_sized(h, 8, 8)
intra_pred_sized(h, 16, 16)
intra_pred_sized(h, 32, 32)
intra_pred_sized(h, 64, 64)
intra_pred_sized(h, 4, 8)
intra_pred_sized(h, 4, 16)
intra_pred_sized(h, 8, 4)
intra_pred_sized(h, 8, 16)
intra_pred_sized(h, 8, 32)
intra_pred_sized(h, 16, 4)
intra_pred_sized(h, 16, 8)
intra_pred_sized(h, 16, 32)
intra_pred_sized(h, 16, 64)
intra_pred_sized(h, 32, 8)
intra_pred_sized(h, 32, 16)
intra_pred_sized(h, 32, 64)
intra_pred_sized(h, 64, 16)
intra_pred_sized(h, 64, 32)
intra_pred_sized(smooth, 2, 2)
intra_pred_sized(smooth, 4, 4)
intra_pred_sized(smooth, 8, 8)
intra_pred_sized(smooth, 16, 16)
intra_pred_sized(smooth, 32, 32)
intra_pred_sized(smooth, 64, 64)
intra_pred_sized(smooth, 4, 8)
intra_pred_sized(smooth, 4, 16)
intra_pred_sized(smooth, 8, 4)
intra_pred_sized(smooth, 8, 16)
intra_pred_sized(smooth, 8, 32)
intra_pred_sized(smooth, 16, 4)
intra_pred_sized(smooth, 16, 8)
intra_pred_sized(smooth, 16, 32)
intra_pred_sized(smooth, 16, 64)
intra_pred_sized(smooth, 32, 8)
intra_pred_sized(smooth, 32, 16)
intra_pred_sized(smooth, 32, 64)
intra_pred_sized(smooth, 64, 16)
intra_pred_sized(smooth, 64, 32)
intra_pred_sized(smooth_h, 2, 2)
intra_pred_sized(smooth_h, 4, 4)
intra_pred_sized(smooth_h, 8, 8)
intra_pred_sized(smooth_h, 16, 16)
intra_pred_sized(smooth_h, 32, 32)
intra_pred_sized(smooth_h, 64, 64)
intra_pred_sized(smooth_h, 4, 8)
intra_pred_sized(smooth_h, 4, 16)
intra_pred_sized(smooth_h, 8, 4)
intra_pred_sized(smooth_h, 8, 16)
intra_pred_sized(smooth_h, 8, 32)
intra_pred_sized(smooth_h, 16, 4)
intra_pred_sized(smooth_h, 16, 8)
intra_pred_sized(smooth_h, 16, 32)
intra_pred_sized(smooth_h, 16, 64)
intra_pred_sized(smooth_h, 32, 8)
intra_pred_sized(smooth_h, 32, 16)
intra_pred_sized(smooth_h, 32, 64)
intra_pred_sized(smooth_h, 64, 16)
intra_pred_sized(smooth_h, 64, 32)
intra_pred_sized(smooth_v, 2, 2)
intra_pred_sized(smooth_v, 4, 4)
intra_pred_sized(smooth_v, 8, 8)
intra_pred_sized(smooth_v, 16, 16)
intra_pred_sized(smooth_v, 32, 32)
intra_pred_sized(smooth_v, 64, 64)
intra_pred_sized(smooth_v, 4, 8)
intra_pred_sized(smooth_v, 4, 16)
intra_pred_sized(smooth_v, 8, 4)
intra_pred_sized(smooth_v, 8, 16)
intra_pred_sized(smooth_v, 8, 32)
intra_pred_sized(smooth_v, 16, 4)
intra_pred_sized(smooth_v, 16, 8)
intra_pred_sized(smooth_v, 16, 32)
intra_pred_sized(smooth_v, 16, 64)
intra_pred_sized(smooth_v, 32, 8)
intra_pred_sized(smooth_v, 32, 16)
intra_pred_sized(smooth_v, 32, 64)
intra_pred_sized(smooth_v, 64, 16)
intra_pred_sized(smooth_v, 64, 32)
intra_pred_sized(paeth, 2, 2)
intra_pred_sized(paeth, 4, 4)
intra_pred_sized(paeth, 8, 8)
intra_pred_sized(paeth, 16, 16)
intra_pred_sized(paeth, 32, 32)
intra_pred_sized(paeth, 64, 64)
intra_pred_sized(paeth, 4, 8)
intra_pred_sized(paeth, 4, 16)
intra_pred_sized(paeth, 8, 4)
intra_pred_sized(paeth, 8, 16)
intra_pred_sized(paeth, 8, 32)
intra_pred_sized(paeth, 16, 4)
intra_pred_sized(paeth, 16, 8)
intra_pred_sized(paeth, 16, 32)
intra_pred_sized(paeth, 16, 64)
intra_pred_sized(paeth, 32, 8)
intra_pred_sized(paeth, 32, 16)
intra_pred_sized(paeth, 32, 64)
intra_pred_sized(paeth, 64, 16)
intra_pred_sized(paeth, 64, 32)
#define intra_pred_highbd_sized(type, width, height)                        \
  void eb_aom_highbd_##type##_predictor_##width##x##height##_c(                \
      uint16_t *dst, ptrdiff_t stride, const uint16_t *above,               \
      const uint16_t *left, int32_t bd) {                                       \
    highbd_##type##_predictor(dst, stride, width, height, above, left, bd); \
  }

intra_pred_highbd_sized(dc, 2, 2)
intra_pred_highbd_sized(dc, 4, 4)
intra_pred_highbd_sized(dc, 8, 8)
intra_pred_highbd_sized(dc, 16, 16)
intra_pred_highbd_sized(dc, 32, 32)
intra_pred_highbd_sized(dc, 64, 64)
intra_pred_highbd_sized(dc, 4, 8)
intra_pred_highbd_sized(dc, 4, 16)
intra_pred_highbd_sized(dc, 8, 4)
intra_pred_highbd_sized(dc, 8, 16)
intra_pred_highbd_sized(dc, 8, 32)
intra_pred_highbd_sized(dc, 16, 4)
intra_pred_highbd_sized(dc, 16, 8)
intra_pred_highbd_sized(dc, 16, 32)
intra_pred_highbd_sized(dc, 16, 64)
intra_pred_highbd_sized(dc, 32, 8)
intra_pred_highbd_sized(dc, 32, 16)
intra_pred_highbd_sized(dc, 32, 64)
intra_pred_highbd_sized(dc, 64, 16)
intra_pred_highbd_sized(dc, 64, 32)

intra_pred_highbd_sized(dc_128, 2, 2)
intra_pred_highbd_sized(dc_128, 4, 4)
intra_pred_highbd_sized(dc_128, 8, 8)
intra_pred_highbd_sized(dc_128, 16, 16)
intra_pred_highbd_sized(dc_128, 32, 32)
intra_pred_highbd_sized(dc_128, 64, 64)

intra_pred_highbd_sized(dc_128, 4, 8)
intra_pred_highbd_sized(dc_128, 4, 16)
intra_pred_highbd_sized(dc_128, 8, 4)
intra_pred_highbd_sized(dc_128, 8, 16)
intra_pred_highbd_sized(dc_128, 8, 32)
intra_pred_highbd_sized(dc_128, 16, 4)
intra_pred_highbd_sized(dc_128, 16, 8)
intra_pred_highbd_sized(dc_128, 16, 32)
intra_pred_highbd_sized(dc_128, 16, 64)
intra_pred_highbd_sized(dc_128, 32, 8)
intra_pred_highbd_sized(dc_128, 32, 16)
intra_pred_highbd_sized(dc_128, 32, 64)
intra_pred_highbd_sized(dc_128, 64, 16)
intra_pred_highbd_sized(dc_128, 64, 32)

intra_pred_highbd_sized(dc_left, 2, 2)
intra_pred_highbd_sized(dc_left, 4, 4)
intra_pred_highbd_sized(dc_left, 8, 8)
intra_pred_highbd_sized(dc_left, 16, 16)
intra_pred_highbd_sized(dc_left, 32, 32)
intra_pred_highbd_sized(dc_left, 64, 64)

intra_pred_highbd_sized(dc_left, 4, 8)
intra_pred_highbd_sized(dc_left, 4, 16)
intra_pred_highbd_sized(dc_left, 8, 4)
intra_pred_highbd_sized(dc_left, 8, 16)
intra_pred_highbd_sized(dc_left, 8, 32)
intra_pred_highbd_sized(dc_left, 16, 4)
intra_pred_highbd_sized(dc_left, 16, 8)
intra_pred_highbd_sized(dc_left, 16, 32)
intra_pred_highbd_sized(dc_left, 16, 64)
intra_pred_highbd_sized(dc_left, 32, 8)
intra_pred_highbd_sized(dc_left, 32, 16)
intra_pred_highbd_sized(dc_left, 32, 64)
intra_pred_highbd_sized(dc_left, 64, 16)
intra_pred_highbd_sized(dc_left, 64, 32)

intra_pred_highbd_sized(dc_top, 2, 2)
intra_pred_highbd_sized(dc_top, 4, 4)
intra_pred_highbd_sized(dc_top, 8, 8)
intra_pred_highbd_sized(dc_top, 16, 16)
intra_pred_highbd_sized(dc_top, 32, 32)
intra_pred_highbd_sized(dc_top, 64, 64)

intra_pred_highbd_sized(dc_top, 4, 8)
intra_pred_highbd_sized(dc_top, 4, 16)
intra_pred_highbd_sized(dc_top, 8, 4)
intra_pred_highbd_sized(dc_top, 8, 16)
intra_pred_highbd_sized(dc_top, 8, 32)
intra_pred_highbd_sized(dc_top, 16, 4)
intra_pred_highbd_sized(dc_top, 16, 8)
intra_pred_highbd_sized(dc_top, 16, 32)
intra_pred_highbd_sized(dc_top, 16, 64)
intra_pred_highbd_sized(dc_top, 32, 8)
intra_pred_highbd_sized(dc_top, 32, 16)
intra_pred_highbd_sized(dc_top, 32, 64)
intra_pred_highbd_sized(dc_top, 64, 16)
intra_pred_highbd_sized(dc_top, 64, 32)

intra_pred_highbd_sized(v, 2, 2)
intra_pred_highbd_sized(v, 4, 4)
intra_pred_highbd_sized(v, 8, 8)
intra_pred_highbd_sized(v, 16, 16)
intra_pred_highbd_sized(v, 32, 32)
intra_pred_highbd_sized(v, 64, 64)

intra_pred_highbd_sized(v, 4, 8)
intra_pred_highbd_sized(v, 4, 16)
intra_pred_highbd_sized(v, 8, 4)
intra_pred_highbd_sized(v, 8, 16)
intra_pred_highbd_sized(v, 8, 32)
intra_pred_highbd_sized(v, 16, 4)
intra_pred_highbd_sized(v, 16, 8)
intra_pred_highbd_sized(v, 16, 32)
intra_pred_highbd_sized(v, 16, 64)
intra_pred_highbd_sized(v, 32, 8)
intra_pred_highbd_sized(v, 32, 16)
intra_pred_highbd_sized(v, 32, 64)
intra_pred_highbd_sized(v, 64, 16)
intra_pred_highbd_sized(v, 64, 32)

intra_pred_highbd_sized(h, 2, 2)
intra_pred_highbd_sized(h, 4, 4)
intra_pred_highbd_sized(h, 8, 8)
intra_pred_highbd_sized(h, 16, 16)
intra_pred_highbd_sized(h, 32, 32)
intra_pred_highbd_sized(h, 64, 64)

intra_pred_highbd_sized(h, 4, 8)
intra_pred_highbd_sized(h, 4, 16)
intra_pred_highbd_sized(h, 8, 4)
intra_pred_highbd_sized(h, 8, 16)
intra_pred_highbd_sized(h, 8, 32)
intra_pred_highbd_sized(h, 16, 4)
intra_pred_highbd_sized(h, 16, 8)
intra_pred_highbd_sized(h, 16, 32)
intra_pred_highbd_sized(h, 16, 64)
intra_pred_highbd_sized(h, 32, 8)
intra_pred_highbd_sized(h, 32, 16)
intra_pred_highbd_sized(h, 32, 64)
intra_pred_highbd_sized(h, 64, 16)
intra_pred_highbd_sized(h, 64, 32)

intra_pred_highbd_sized(smooth, 2, 2)
intra_pred_highbd_sized(smooth, 4, 4)
intra_pred_highbd_sized(smooth, 8, 8)
intra_pred_highbd_sized(smooth, 16, 16)
intra_pred_highbd_sized(smooth, 32, 32)
intra_pred_highbd_sized(smooth, 64, 64)

intra_pred_highbd_sized(smooth, 4, 8)
intra_pred_highbd_sized(smooth, 4, 16)
intra_pred_highbd_sized(smooth, 8, 4)
intra_pred_highbd_sized(smooth, 8, 16)
intra_pred_highbd_sized(smooth, 8, 32)
intra_pred_highbd_sized(smooth, 16, 4)
intra_pred_highbd_sized(smooth, 16, 8)
intra_pred_highbd_sized(smooth, 16, 32)
intra_pred_highbd_sized(smooth, 16, 64)
intra_pred_highbd_sized(smooth, 32, 8)
intra_pred_highbd_sized(smooth, 32, 16)
intra_pred_highbd_sized(smooth, 32, 64)
intra_pred_highbd_sized(smooth, 64, 16)
intra_pred_highbd_sized(smooth, 64, 32)

intra_pred_highbd_sized(smooth_h, 2, 2)
intra_pred_highbd_sized(smooth_h, 4, 4)
intra_pred_highbd_sized(smooth_h, 8, 8)
intra_pred_highbd_sized(smooth_h, 16, 16)
intra_pred_highbd_sized(smooth_h, 32, 32)
intra_pred_highbd_sized(smooth_h, 64, 64)

intra_pred_highbd_sized(smooth_h, 4, 8)
intra_pred_highbd_sized(smooth_h, 4, 16)
intra_pred_highbd_sized(smooth_h, 8, 4)
intra_pred_highbd_sized(smooth_h, 8, 16)
intra_pred_highbd_sized(smooth_h, 8, 32)
intra_pred_highbd_sized(smooth_h, 16, 4)
intra_pred_highbd_sized(smooth_h, 16, 8)
intra_pred_highbd_sized(smooth_h, 16, 32)
intra_pred_highbd_sized(smooth_h, 16, 64)
intra_pred_highbd_sized(smooth_h, 32, 8)
intra_pred_highbd_sized(smooth_h, 32, 16)
intra_pred_highbd_sized(smooth_h, 32, 64)
intra_pred_highbd_sized(smooth_h, 64, 16)
intra_pred_highbd_sized(smooth_h, 64, 32)

intra_pred_highbd_sized(smooth_v, 2, 2)
intra_pred_highbd_sized(smooth_v, 4, 4)
intra_pred_highbd_sized(smooth_v, 8, 8)
intra_pred_highbd_sized(smooth_v, 16, 16)
intra_pred_highbd_sized(smooth_v, 32, 32)
intra_pred_highbd_sized(smooth_v, 64, 64)

intra_pred_highbd_sized(smooth_v, 4, 8)
intra_pred_highbd_sized(smooth_v, 4, 16)
intra_pred_highbd_sized(smooth_v, 8, 4)
intra_pred_highbd_sized(smooth_v, 8, 16)
intra_pred_highbd_sized(smooth_v, 8, 32)
intra_pred_highbd_sized(smooth_v, 16, 4)
intra_pred_highbd_sized(smooth_v, 16, 8)
intra_pred_highbd_sized(smooth_v, 16, 32)
intra_pred_highbd_sized(smooth_v, 16, 64)
intra_pred_highbd_sized(smooth_v, 32, 8)
intra_pred_highbd_sized(smooth_v, 32, 16)
intra_pred_highbd_sized(smooth_v, 32, 64)
intra_pred_highbd_sized(smooth_v, 64, 16)
intra_pred_highbd_sized(smooth_v, 64, 32)

intra_pred_highbd_sized(paeth, 2, 2)
intra_pred_highbd_sized(paeth, 4, 4)
intra_pred_highbd_sized(paeth, 8, 8)
intra_pred_highbd_sized(paeth, 16, 16)
intra_pred_highbd_sized(paeth, 32, 32)
intra_pred_highbd_sized(paeth, 64, 64)
intra_pred_highbd_sized(paeth, 4, 8)
intra_pred_highbd_sized(paeth, 4, 16)
intra_pred_highbd_sized(paeth, 8, 4)
intra_pred_highbd_sized(paeth, 8, 16)
intra_pred_highbd_sized(paeth, 8, 32)
intra_pred_highbd_sized(paeth, 16, 4)
intra_pred_highbd_sized(paeth, 16, 8)
intra_pred_highbd_sized(paeth, 16, 32)
intra_pred_highbd_sized(paeth, 16, 64)
intra_pred_highbd_sized(paeth, 32, 8)
intra_pred_highbd_sized(paeth, 32, 16)
intra_pred_highbd_sized(paeth, 32, 64)
intra_pred_highbd_sized(paeth, 64, 16)
intra_pred_highbd_sized(paeth, 64, 32)

IntraPredFnC  dc_pred_c[2][2];
IntraHighBdPredFnC  highbd_dc_pred_c[2][2];
void init_intra_dc_predictors_c_internal(void)
{
    dc_pred_c[0][0] = dc_128_predictor;
    dc_pred_c[0][1] = dc_top_predictor;
    dc_pred_c[1][0] = dc_left_predictor;
    dc_pred_c[1][1] = dc_predictor;

    highbd_dc_pred_c[0][0] = highbd_dc_128_predictor;
    highbd_dc_pred_c[0][1] = highbd_dc_top_predictor;
    highbd_dc_pred_c[1][0] = highbd_dc_left_predictor;
    highbd_dc_pred_c[1][1] = highbd_dc_predictor;
}

/*static*/ void init_intra_predictors_internal(void) {
    pred[V_PRED][TX_4X4] = eb_aom_v_predictor_4x4;
    pred[V_PRED][TX_8X8] = eb_aom_v_predictor_8x8;
    pred[V_PRED][TX_16X16] = eb_aom_v_predictor_16x16;
    pred[V_PRED][TX_32X32] = eb_aom_v_predictor_32x32;
    pred[V_PRED][TX_64X64] = eb_aom_v_predictor_64x64;
    pred[V_PRED][TX_4X8] = eb_aom_v_predictor_4x8;
    pred[V_PRED][TX_4X16] = eb_aom_v_predictor_4x16;

    pred[V_PRED][TX_8X4] = eb_aom_v_predictor_8x4;
    pred[V_PRED][TX_8X16] = eb_aom_v_predictor_8x16;
    pred[V_PRED][TX_8X32] = eb_aom_v_predictor_8x32;

    pred[V_PRED][TX_16X4] = eb_aom_v_predictor_16x4;
    pred[V_PRED][TX_16X8] = eb_aom_v_predictor_16x8;
    pred[V_PRED][TX_16X32] = eb_aom_v_predictor_16x32;
    pred[V_PRED][TX_16X64] = eb_aom_v_predictor_16x64;

    pred[V_PRED][TX_32X8] = eb_aom_v_predictor_32x8;
    pred[V_PRED][TX_32X16] = eb_aom_v_predictor_32x16;
    pred[V_PRED][TX_32X64] = eb_aom_v_predictor_32x64;

    pred[V_PRED][TX_64X16] = eb_aom_v_predictor_64x16;
    pred[V_PRED][TX_64X32] = eb_aom_v_predictor_64x32;

    pred[H_PRED][TX_4X4] = eb_aom_h_predictor_4x4;
    pred[H_PRED][TX_8X8] = eb_aom_h_predictor_8x8;
    pred[H_PRED][TX_16X16] = eb_aom_h_predictor_16x16;
    pred[H_PRED][TX_32X32] = eb_aom_h_predictor_32x32;
    pred[H_PRED][TX_64X64] = eb_aom_h_predictor_64x64;

    pred[H_PRED][TX_4X8] = eb_aom_h_predictor_4x8;
    pred[H_PRED][TX_4X16] = eb_aom_h_predictor_4x16;

    pred[H_PRED][TX_8X4] = eb_aom_h_predictor_8x4;
    pred[H_PRED][TX_8X16] = eb_aom_h_predictor_8x16;
    pred[H_PRED][TX_8X32] = eb_aom_h_predictor_8x32;

    pred[H_PRED][TX_16X4] = eb_aom_h_predictor_16x4;
    pred[H_PRED][TX_16X8] = eb_aom_h_predictor_16x8;
    pred[H_PRED][TX_16X32] = eb_aom_h_predictor_16x32;
    pred[H_PRED][TX_16X64] = eb_aom_h_predictor_16x64;

    pred[H_PRED][TX_32X8] = eb_aom_h_predictor_32x8;
    pred[H_PRED][TX_32X16] = eb_aom_h_predictor_32x16;
    pred[H_PRED][TX_32X64] = eb_aom_h_predictor_32x64;

    pred[H_PRED][TX_64X16] = eb_aom_h_predictor_64x16;
    pred[H_PRED][TX_64X32] = eb_aom_h_predictor_64x32;

    pred[SMOOTH_PRED][TX_4X4] = eb_aom_smooth_predictor_4x4;
    pred[SMOOTH_PRED][TX_8X8] = eb_aom_smooth_predictor_8x8;
    pred[SMOOTH_PRED][TX_16X16] = eb_aom_smooth_predictor_16x16;
    pred[SMOOTH_PRED][TX_32X32] = eb_aom_smooth_predictor_32x32;
    pred[SMOOTH_PRED][TX_64X64] = eb_aom_smooth_predictor_64x64;

    pred[SMOOTH_PRED][TX_4X8] = eb_aom_smooth_predictor_4x8;
    pred[SMOOTH_PRED][TX_4X16] = eb_aom_smooth_predictor_4x16;

    pred[SMOOTH_PRED][TX_8X4] = eb_aom_smooth_predictor_8x4;
    pred[SMOOTH_PRED][TX_8X16] = eb_aom_smooth_predictor_8x16;
    pred[SMOOTH_PRED][TX_8X32] = eb_aom_smooth_predictor_8x32;

    pred[SMOOTH_PRED][TX_16X4] = eb_aom_smooth_predictor_16x4;
    pred[SMOOTH_PRED][TX_16X8] = eb_aom_smooth_predictor_16x8;
    pred[SMOOTH_PRED][TX_16X32] = eb_aom_smooth_predictor_16x32;
    pred[SMOOTH_PRED][TX_16X64] = eb_aom_smooth_predictor_16x64;

    pred[SMOOTH_PRED][TX_32X8] = eb_aom_smooth_predictor_32x8;
    pred[SMOOTH_PRED][TX_32X16] = eb_aom_smooth_predictor_32x16;
    pred[SMOOTH_PRED][TX_32X64] = eb_aom_smooth_predictor_32x64;

    pred[SMOOTH_PRED][TX_64X16] = eb_aom_smooth_predictor_64x16;
    pred[SMOOTH_PRED][TX_64X32] = eb_aom_smooth_predictor_64x32;

    pred[SMOOTH_V_PRED][TX_4X4] = eb_aom_smooth_v_predictor_4x4;
    pred[SMOOTH_V_PRED][TX_8X8] = eb_aom_smooth_v_predictor_8x8;
    pred[SMOOTH_V_PRED][TX_16X16] = eb_aom_smooth_v_predictor_16x16;
    pred[SMOOTH_V_PRED][TX_32X32] = eb_aom_smooth_v_predictor_32x32;
    pred[SMOOTH_V_PRED][TX_64X64] = eb_aom_smooth_v_predictor_64x64;

    pred[SMOOTH_V_PRED][TX_4X8] = eb_aom_smooth_v_predictor_4x8;
    pred[SMOOTH_V_PRED][TX_4X16] = eb_aom_smooth_v_predictor_4x16;

    pred[SMOOTH_V_PRED][TX_8X4] = eb_aom_smooth_v_predictor_8x4;
    pred[SMOOTH_V_PRED][TX_8X16] = eb_aom_smooth_v_predictor_8x16;
    pred[SMOOTH_V_PRED][TX_8X32] = eb_aom_smooth_v_predictor_8x32;

    pred[SMOOTH_V_PRED][TX_16X4] = eb_aom_smooth_v_predictor_16x4;
    pred[SMOOTH_V_PRED][TX_16X8] = eb_aom_smooth_v_predictor_16x8;
    pred[SMOOTH_V_PRED][TX_16X32] = eb_aom_smooth_v_predictor_16x32;
    pred[SMOOTH_V_PRED][TX_16X64] = eb_aom_smooth_v_predictor_16x64;

    pred[SMOOTH_V_PRED][TX_32X8] = eb_aom_smooth_v_predictor_32x8;
    pred[SMOOTH_V_PRED][TX_32X16] = eb_aom_smooth_v_predictor_32x16;
    pred[SMOOTH_V_PRED][TX_32X64] = eb_aom_smooth_v_predictor_32x64;

    pred[SMOOTH_V_PRED][TX_64X16] = eb_aom_smooth_v_predictor_64x16;
    pred[SMOOTH_V_PRED][TX_64X32] = eb_aom_smooth_v_predictor_64x32;

    pred[SMOOTH_H_PRED][TX_4X4] = eb_aom_smooth_h_predictor_4x4;
    pred[SMOOTH_H_PRED][TX_8X8] = eb_aom_smooth_h_predictor_8x8;
    pred[SMOOTH_H_PRED][TX_16X16] = eb_aom_smooth_h_predictor_16x16;
    pred[SMOOTH_H_PRED][TX_32X32] = eb_aom_smooth_h_predictor_32x32;
    pred[SMOOTH_H_PRED][TX_64X64] = eb_aom_smooth_h_predictor_64x64;

    pred[SMOOTH_H_PRED][TX_4X8] = eb_aom_smooth_h_predictor_4x8;
    pred[SMOOTH_H_PRED][TX_4X16] = eb_aom_smooth_h_predictor_4x16;

    pred[SMOOTH_H_PRED][TX_8X4] = eb_aom_smooth_h_predictor_8x4;
    pred[SMOOTH_H_PRED][TX_8X16] = eb_aom_smooth_h_predictor_8x16;
    pred[SMOOTH_H_PRED][TX_8X32] = eb_aom_smooth_h_predictor_8x32;

    pred[SMOOTH_H_PRED][TX_16X4] = eb_aom_smooth_h_predictor_16x4;
    pred[SMOOTH_H_PRED][TX_16X8] = eb_aom_smooth_h_predictor_16x8;
    pred[SMOOTH_H_PRED][TX_16X32] = eb_aom_smooth_h_predictor_16x32;
    pred[SMOOTH_H_PRED][TX_16X64] = eb_aom_smooth_h_predictor_16x64;

    pred[SMOOTH_H_PRED][TX_32X8] = eb_aom_smooth_h_predictor_32x8;
    pred[SMOOTH_H_PRED][TX_32X16] = eb_aom_smooth_h_predictor_32x16;
    pred[SMOOTH_H_PRED][TX_32X64] = eb_aom_smooth_h_predictor_32x64;

    pred[SMOOTH_H_PRED][TX_64X16] = eb_aom_smooth_h_predictor_64x16;
    pred[SMOOTH_H_PRED][TX_64X32] = eb_aom_smooth_h_predictor_64x32;

    pred[PAETH_PRED][TX_4X4] = eb_aom_paeth_predictor_4x4;
    pred[PAETH_PRED][TX_8X8] = eb_aom_paeth_predictor_8x8;
    pred[PAETH_PRED][TX_16X16] = eb_aom_paeth_predictor_16x16;
    pred[PAETH_PRED][TX_32X32] = eb_aom_paeth_predictor_32x32;
    pred[PAETH_PRED][TX_64X64] = eb_aom_paeth_predictor_64x64;

    pred[PAETH_PRED][TX_4X8] = eb_aom_paeth_predictor_4x8;
    pred[PAETH_PRED][TX_4X16] = eb_aom_paeth_predictor_4x16;

    pred[PAETH_PRED][TX_8X4] = eb_aom_paeth_predictor_8x4;
    pred[PAETH_PRED][TX_8X16] = eb_aom_paeth_predictor_8x16;
    pred[PAETH_PRED][TX_8X32] = eb_aom_paeth_predictor_8x32;

    pred[PAETH_PRED][TX_16X4] = eb_aom_paeth_predictor_16x4;
    pred[PAETH_PRED][TX_16X8] = eb_aom_paeth_predictor_16x8;
    pred[PAETH_PRED][TX_16X32] = eb_aom_paeth_predictor_16x32;
    pred[PAETH_PRED][TX_16X64] = eb_aom_paeth_predictor_16x64;

    pred[PAETH_PRED][TX_32X8] = eb_aom_paeth_predictor_32x8;
    pred[PAETH_PRED][TX_32X16] = eb_aom_paeth_predictor_32x16;
    pred[PAETH_PRED][TX_32X64] = eb_aom_paeth_predictor_32x64;

    pred[PAETH_PRED][TX_64X16] = eb_aom_paeth_predictor_64x16;
    pred[PAETH_PRED][TX_64X32] = eb_aom_paeth_predictor_64x32;
    dc_pred[0][0][TX_4X4] = eb_aom_dc_128_predictor_4x4;
    dc_pred[0][0][TX_8X8] = eb_aom_dc_128_predictor_8x8;
    dc_pred[0][0][TX_16X16] = eb_aom_dc_128_predictor_16x16;
    dc_pred[0][0][TX_32X32] = eb_aom_dc_128_predictor_32x32;
    dc_pred[0][0][TX_64X64] = eb_aom_dc_128_predictor_64x64;

    dc_pred[0][0][TX_4X8] = eb_aom_dc_128_predictor_4x8;
    dc_pred[0][0][TX_4X16] = eb_aom_dc_128_predictor_4x16;

    dc_pred[0][0][TX_8X4] = eb_aom_dc_128_predictor_8x4;
    dc_pred[0][0][TX_8X16] = eb_aom_dc_128_predictor_8x16;
    dc_pred[0][0][TX_8X32] = eb_aom_dc_128_predictor_8x32;

    dc_pred[0][0][TX_16X4] = eb_aom_dc_128_predictor_16x4;
    dc_pred[0][0][TX_16X8] = eb_aom_dc_128_predictor_16x8;
    dc_pred[0][0][TX_16X32] = eb_aom_dc_128_predictor_16x32;
    dc_pred[0][0][TX_16X64] = eb_aom_dc_128_predictor_16x64;

    dc_pred[0][0][TX_32X8] = eb_aom_dc_128_predictor_32x8;
    dc_pred[0][0][TX_32X16] = eb_aom_dc_128_predictor_32x16;
    dc_pred[0][0][TX_32X64] = eb_aom_dc_128_predictor_32x64;

    dc_pred[0][0][TX_64X16] = eb_aom_dc_128_predictor_64x16;
    dc_pred[0][0][TX_64X32] = eb_aom_dc_128_predictor_64x32;

    dc_pred[0][1][TX_4X4] = eb_aom_dc_top_predictor_4x4;
    dc_pred[0][1][TX_8X8] = eb_aom_dc_top_predictor_8x8;
    dc_pred[0][1][TX_16X16] = eb_aom_dc_top_predictor_16x16;
    dc_pred[0][1][TX_32X32] = eb_aom_dc_top_predictor_32x32;
    dc_pred[0][1][TX_64X64] = eb_aom_dc_top_predictor_64x64;

    dc_pred[0][1][TX_4X8] = eb_aom_dc_top_predictor_4x8;
    dc_pred[0][1][TX_4X16] = eb_aom_dc_top_predictor_4x16;

    dc_pred[0][1][TX_8X4] = eb_aom_dc_top_predictor_8x4;
    dc_pred[0][1][TX_8X16] = eb_aom_dc_top_predictor_8x16;
    dc_pred[0][1][TX_8X32] = eb_aom_dc_top_predictor_8x32;

    dc_pred[0][1][TX_16X4] = eb_aom_dc_top_predictor_16x4;
    dc_pred[0][1][TX_16X8] = eb_aom_dc_top_predictor_16x8;
    dc_pred[0][1][TX_16X32] = eb_aom_dc_top_predictor_16x32;
    dc_pred[0][1][TX_16X64] = eb_aom_dc_top_predictor_16x64;

    dc_pred[0][1][TX_32X8] = eb_aom_dc_top_predictor_32x8;
    dc_pred[0][1][TX_32X16] = eb_aom_dc_top_predictor_32x16;
    dc_pred[0][1][TX_32X64] = eb_aom_dc_top_predictor_32x64;

    dc_pred[0][1][TX_64X16] = eb_aom_dc_top_predictor_64x16;
    dc_pred[0][1][TX_64X32] = eb_aom_dc_top_predictor_64x32;

    dc_pred[1][0][TX_4X4] = eb_aom_dc_left_predictor_4x4;
    dc_pred[1][0][TX_8X8] = eb_aom_dc_left_predictor_8x8;
    dc_pred[1][0][TX_16X16] = eb_aom_dc_left_predictor_16x16;
    dc_pred[1][0][TX_32X32] = eb_aom_dc_left_predictor_32x32;
    dc_pred[1][0][TX_64X64] = eb_aom_dc_left_predictor_64x64;
    dc_pred[1][0][TX_4X8] = eb_aom_dc_left_predictor_4x8;
    dc_pred[1][0][TX_4X16] = eb_aom_dc_left_predictor_4x16;

    dc_pred[1][0][TX_8X4] = eb_aom_dc_left_predictor_8x4;
    dc_pred[1][0][TX_8X16] = eb_aom_dc_left_predictor_8x16;
    dc_pred[1][0][TX_8X32] = eb_aom_dc_left_predictor_8x32;

    dc_pred[1][0][TX_16X4] = eb_aom_dc_left_predictor_16x4;
    dc_pred[1][0][TX_16X8] = eb_aom_dc_left_predictor_16x8;
    dc_pred[1][0][TX_16X32] = eb_aom_dc_left_predictor_16x32;
    dc_pred[1][0][TX_16X64] = eb_aom_dc_left_predictor_16x64;

    dc_pred[1][0][TX_32X8] = eb_aom_dc_left_predictor_32x8;
    dc_pred[1][0][TX_32X16] = eb_aom_dc_left_predictor_32x16;
    dc_pred[1][0][TX_32X64] = eb_aom_dc_left_predictor_32x64;

    dc_pred[1][0][TX_64X16] = eb_aom_dc_left_predictor_64x16;
    dc_pred[1][0][TX_64X32] = eb_aom_dc_left_predictor_64x32;

    dc_pred[1][1][TX_4X4] = eb_aom_dc_predictor_4x4;
    dc_pred[1][1][TX_8X8] = eb_aom_dc_predictor_8x8;
    dc_pred[1][1][TX_16X16] = eb_aom_dc_predictor_16x16;
    dc_pred[1][1][TX_32X32] = eb_aom_dc_predictor_32x32;
    dc_pred[1][1][TX_64X64] = eb_aom_dc_predictor_64x64;
    dc_pred[1][1][TX_4X8] = eb_aom_dc_predictor_4x8;
    dc_pred[1][1][TX_4X16] = eb_aom_dc_predictor_4x16;

    dc_pred[1][1][TX_8X4] = eb_aom_dc_predictor_8x4;
    dc_pred[1][1][TX_8X16] = eb_aom_dc_predictor_8x16;
    dc_pred[1][1][TX_8X32] = eb_aom_dc_predictor_8x32;

    dc_pred[1][1][TX_16X4] = eb_aom_dc_predictor_16x4;
    dc_pred[1][1][TX_16X8] = eb_aom_dc_predictor_16x8;
    dc_pred[1][1][TX_16X32] = eb_aom_dc_predictor_16x32;
    dc_pred[1][1][TX_16X64] = eb_aom_dc_predictor_16x64;

    dc_pred[1][1][TX_32X8] = eb_aom_dc_predictor_32x8;
    dc_pred[1][1][TX_32X16] = eb_aom_dc_predictor_32x16;
    dc_pred[1][1][TX_32X64] = eb_aom_dc_predictor_32x64;

    dc_pred[1][1][TX_64X16] = eb_aom_dc_predictor_64x16;
    dc_pred[1][1][TX_64X32] = eb_aom_dc_predictor_64x32;

    pred_high[V_PRED][TX_4X4] = eb_aom_highbd_v_predictor_4x4;
    pred_high[V_PRED][TX_8X8] = eb_aom_highbd_v_predictor_8x8;
    pred_high[V_PRED][TX_16X16] = eb_aom_highbd_v_predictor_16x16;
    pred_high[V_PRED][TX_32X32] = eb_aom_highbd_v_predictor_32x32;
    pred_high[V_PRED][TX_64X64] = eb_aom_highbd_v_predictor_64x64;

    pred_high[V_PRED][TX_4X8] = eb_aom_highbd_v_predictor_4x8;
    pred_high[V_PRED][TX_4X16] = eb_aom_highbd_v_predictor_4x16;

    pred_high[V_PRED][TX_8X4] = eb_aom_highbd_v_predictor_8x4;
    pred_high[V_PRED][TX_8X16] = eb_aom_highbd_v_predictor_8x16;
    pred_high[V_PRED][TX_8X32] = eb_aom_highbd_v_predictor_8x32;

    pred_high[V_PRED][TX_16X4] = eb_aom_highbd_v_predictor_16x4;
    pred_high[V_PRED][TX_16X8] = eb_aom_highbd_v_predictor_16x8;
    pred_high[V_PRED][TX_16X32] = eb_aom_highbd_v_predictor_16x32;
    pred_high[V_PRED][TX_16X64] = eb_aom_highbd_v_predictor_16x64;

    pred_high[V_PRED][TX_32X8] = eb_aom_highbd_v_predictor_32x8;
    pred_high[V_PRED][TX_32X16] = eb_aom_highbd_v_predictor_32x16;
    pred_high[V_PRED][TX_32X64] = eb_aom_highbd_v_predictor_32x64;

    pred_high[V_PRED][TX_64X16] = eb_aom_highbd_v_predictor_64x16;
    pred_high[V_PRED][TX_64X32] = eb_aom_highbd_v_predictor_64x32;

    pred_high[H_PRED][TX_4X4] = eb_aom_highbd_h_predictor_4x4;
    pred_high[H_PRED][TX_8X8] = eb_aom_highbd_h_predictor_8x8;
    pred_high[H_PRED][TX_16X16] = eb_aom_highbd_h_predictor_16x16;
    pred_high[H_PRED][TX_32X32] = eb_aom_highbd_h_predictor_32x32;
    pred_high[H_PRED][TX_64X64] = eb_aom_highbd_h_predictor_64x64;

    pred_high[H_PRED][TX_4X8] = eb_aom_highbd_h_predictor_4x8;
    pred_high[H_PRED][TX_4X16] = eb_aom_highbd_h_predictor_4x16;

    pred_high[H_PRED][TX_8X4] = eb_aom_highbd_h_predictor_8x4;
    pred_high[H_PRED][TX_8X16] = eb_aom_highbd_h_predictor_8x16;
    pred_high[H_PRED][TX_8X32] = eb_aom_highbd_h_predictor_8x32;

    pred_high[H_PRED][TX_16X4] = eb_aom_highbd_h_predictor_16x4;
    pred_high[H_PRED][TX_16X8] = eb_aom_highbd_h_predictor_16x8;
    pred_high[H_PRED][TX_16X32] = eb_aom_highbd_h_predictor_16x32;
    pred_high[H_PRED][TX_16X64] = eb_aom_highbd_h_predictor_16x64;

    pred_high[H_PRED][TX_32X8] = eb_aom_highbd_h_predictor_32x8;
    pred_high[H_PRED][TX_32X16] = eb_aom_highbd_h_predictor_32x16;
    pred_high[H_PRED][TX_32X64] = eb_aom_highbd_h_predictor_32x64;

    pred_high[H_PRED][TX_64X16] = eb_aom_highbd_h_predictor_64x16;
    pred_high[H_PRED][TX_64X32] = eb_aom_highbd_h_predictor_64x32;

    pred_high[SMOOTH_PRED][TX_4X4] = eb_aom_highbd_smooth_predictor_4x4;
    pred_high[SMOOTH_PRED][TX_8X8] = eb_aom_highbd_smooth_predictor_8x8;
    pred_high[SMOOTH_PRED][TX_16X16] = eb_aom_highbd_smooth_predictor_16x16;
    pred_high[SMOOTH_PRED][TX_32X32] = eb_aom_highbd_smooth_predictor_32x32;
    pred_high[SMOOTH_PRED][TX_64X64] = eb_aom_highbd_smooth_predictor_64x64;

    pred_high[SMOOTH_PRED][TX_4X8] = eb_aom_highbd_smooth_predictor_4x8;
    pred_high[SMOOTH_PRED][TX_4X16] = eb_aom_highbd_smooth_predictor_4x16;

    pred_high[SMOOTH_PRED][TX_8X4] = eb_aom_highbd_smooth_predictor_8x4;
    pred_high[SMOOTH_PRED][TX_8X16] = eb_aom_highbd_smooth_predictor_8x16;
    pred_high[SMOOTH_PRED][TX_8X32] = eb_aom_highbd_smooth_predictor_8x32;

    pred_high[SMOOTH_PRED][TX_16X4] = eb_aom_highbd_smooth_predictor_16x4;
    pred_high[SMOOTH_PRED][TX_16X8] = eb_aom_highbd_smooth_predictor_16x8;
    pred_high[SMOOTH_PRED][TX_16X32] = eb_aom_highbd_smooth_predictor_16x32;
    pred_high[SMOOTH_PRED][TX_16X64] = eb_aom_highbd_smooth_predictor_16x64;

    pred_high[SMOOTH_PRED][TX_32X8] = eb_aom_highbd_smooth_predictor_32x8;
    pred_high[SMOOTH_PRED][TX_32X16] = eb_aom_highbd_smooth_predictor_32x16;
    pred_high[SMOOTH_PRED][TX_32X64] = eb_aom_highbd_smooth_predictor_32x64;

    pred_high[SMOOTH_PRED][TX_64X16] = eb_aom_highbd_smooth_predictor_64x16;
    pred_high[SMOOTH_PRED][TX_64X32] = eb_aom_highbd_smooth_predictor_64x32;

    pred_high[SMOOTH_V_PRED][TX_4X4] = eb_aom_highbd_smooth_v_predictor_4x4;
    pred_high[SMOOTH_V_PRED][TX_8X8] = eb_aom_highbd_smooth_v_predictor_8x8;
    pred_high[SMOOTH_V_PRED][TX_16X16] = eb_aom_highbd_smooth_v_predictor_16x16;
    pred_high[SMOOTH_V_PRED][TX_32X32] = eb_aom_highbd_smooth_v_predictor_32x32;
    pred_high[SMOOTH_V_PRED][TX_64X64] = eb_aom_highbd_smooth_v_predictor_64x64;

    pred_high[SMOOTH_V_PRED][TX_4X8] = eb_aom_highbd_smooth_v_predictor_4x8;
    pred_high[SMOOTH_V_PRED][TX_4X16] = eb_aom_highbd_smooth_v_predictor_4x16;

    pred_high[SMOOTH_V_PRED][TX_8X4] = eb_aom_highbd_smooth_v_predictor_8x4;
    pred_high[SMOOTH_V_PRED][TX_8X16] = eb_aom_highbd_smooth_v_predictor_8x16;
    pred_high[SMOOTH_V_PRED][TX_8X32] = eb_aom_highbd_smooth_v_predictor_8x32;

    pred_high[SMOOTH_V_PRED][TX_16X4] = eb_aom_highbd_smooth_v_predictor_16x4;
    pred_high[SMOOTH_V_PRED][TX_16X8] = eb_aom_highbd_smooth_v_predictor_16x8;
    pred_high[SMOOTH_V_PRED][TX_16X32] = eb_aom_highbd_smooth_v_predictor_16x32;
    pred_high[SMOOTH_V_PRED][TX_16X64] = eb_aom_highbd_smooth_v_predictor_16x64;

    pred_high[SMOOTH_V_PRED][TX_32X8] = eb_aom_highbd_smooth_v_predictor_32x8;
    pred_high[SMOOTH_V_PRED][TX_32X16] = eb_aom_highbd_smooth_v_predictor_32x16;
    pred_high[SMOOTH_V_PRED][TX_32X64] = eb_aom_highbd_smooth_v_predictor_32x64;

    pred_high[SMOOTH_V_PRED][TX_64X16] = eb_aom_highbd_smooth_v_predictor_64x16;
    pred_high[SMOOTH_V_PRED][TX_64X32] = eb_aom_highbd_smooth_v_predictor_64x32;

    pred_high[SMOOTH_H_PRED][TX_4X4] = eb_aom_highbd_smooth_h_predictor_4x4;
    pred_high[SMOOTH_H_PRED][TX_8X8] = eb_aom_highbd_smooth_h_predictor_8x8;
    pred_high[SMOOTH_H_PRED][TX_16X16] = eb_aom_highbd_smooth_h_predictor_16x16;
    pred_high[SMOOTH_H_PRED][TX_32X32] = eb_aom_highbd_smooth_h_predictor_32x32;
    pred_high[SMOOTH_H_PRED][TX_64X64] = eb_aom_highbd_smooth_h_predictor_64x64;

    pred_high[SMOOTH_H_PRED][TX_4X8] = eb_aom_highbd_smooth_h_predictor_4x8;
    pred_high[SMOOTH_H_PRED][TX_4X16] = eb_aom_highbd_smooth_h_predictor_4x16;

    pred_high[SMOOTH_H_PRED][TX_8X4] = eb_aom_highbd_smooth_h_predictor_8x4;
    pred_high[SMOOTH_H_PRED][TX_8X16] = eb_aom_highbd_smooth_h_predictor_8x16;
    pred_high[SMOOTH_H_PRED][TX_8X32] = eb_aom_highbd_smooth_h_predictor_8x32;

    pred_high[SMOOTH_H_PRED][TX_16X4] = eb_aom_highbd_smooth_h_predictor_16x4;
    pred_high[SMOOTH_H_PRED][TX_16X8] = eb_aom_highbd_smooth_h_predictor_16x8;
    pred_high[SMOOTH_H_PRED][TX_16X32] = eb_aom_highbd_smooth_h_predictor_16x32;
    pred_high[SMOOTH_H_PRED][TX_16X64] = eb_aom_highbd_smooth_h_predictor_16x64;

    pred_high[SMOOTH_H_PRED][TX_32X8] = eb_aom_highbd_smooth_h_predictor_32x8;
    pred_high[SMOOTH_H_PRED][TX_32X16] = eb_aom_highbd_smooth_h_predictor_32x16;
    pred_high[SMOOTH_H_PRED][TX_32X64] = eb_aom_highbd_smooth_h_predictor_32x64;

    pred_high[SMOOTH_H_PRED][TX_64X16] = eb_aom_highbd_smooth_h_predictor_64x16;
    pred_high[SMOOTH_H_PRED][TX_64X32] = eb_aom_highbd_smooth_h_predictor_64x32;

    pred_high[PAETH_PRED][TX_4X4] = eb_aom_highbd_paeth_predictor_4x4;
    pred_high[PAETH_PRED][TX_8X8] = eb_aom_highbd_paeth_predictor_8x8;
    pred_high[PAETH_PRED][TX_16X16] = eb_aom_highbd_paeth_predictor_16x16;
    pred_high[PAETH_PRED][TX_32X32] = eb_aom_highbd_paeth_predictor_32x32;
    pred_high[PAETH_PRED][TX_64X64] = eb_aom_highbd_paeth_predictor_64x64;

    pred_high[PAETH_PRED][TX_4X8] = eb_aom_highbd_paeth_predictor_4x8;
    pred_high[PAETH_PRED][TX_4X16] = eb_aom_highbd_paeth_predictor_4x16;

    pred_high[PAETH_PRED][TX_8X4] = eb_aom_highbd_paeth_predictor_8x4;
    pred_high[PAETH_PRED][TX_8X16] = eb_aom_highbd_paeth_predictor_8x16;
    pred_high[PAETH_PRED][TX_8X32] = eb_aom_highbd_paeth_predictor_8x32;

    pred_high[PAETH_PRED][TX_16X4] = eb_aom_highbd_paeth_predictor_16x4;
    pred_high[PAETH_PRED][TX_16X8] = eb_aom_highbd_paeth_predictor_16x8;
    pred_high[PAETH_PRED][TX_16X32] = eb_aom_highbd_paeth_predictor_16x32;
    pred_high[PAETH_PRED][TX_16X64] = eb_aom_highbd_paeth_predictor_16x64;

    pred_high[PAETH_PRED][TX_32X8] = eb_aom_highbd_paeth_predictor_32x8;
    pred_high[PAETH_PRED][TX_32X16] = eb_aom_highbd_paeth_predictor_32x16;
    pred_high[PAETH_PRED][TX_32X64] = eb_aom_highbd_paeth_predictor_32x64;

    pred_high[PAETH_PRED][TX_64X16] = eb_aom_highbd_paeth_predictor_64x16;
    pred_high[PAETH_PRED][TX_64X32] = eb_aom_highbd_paeth_predictor_64x32;
    dc_pred_high[0][0][TX_4X4] = eb_aom_highbd_dc_128_predictor_4x4;
    dc_pred_high[0][0][TX_8X8] = eb_aom_highbd_dc_128_predictor_8x8;
    dc_pred_high[0][0][TX_16X16] = eb_aom_highbd_dc_128_predictor_16x16;
    dc_pred_high[0][0][TX_32X32] = eb_aom_highbd_dc_128_predictor_32x32;
    dc_pred_high[0][0][TX_64X64] = eb_aom_highbd_dc_128_predictor_64x64;

    dc_pred_high[0][0][TX_4X8] = eb_aom_highbd_dc_128_predictor_4x8;
    dc_pred_high[0][0][TX_4X16] = eb_aom_highbd_dc_128_predictor_4x16;

    dc_pred_high[0][0][TX_8X4] = eb_aom_highbd_dc_128_predictor_8x4;
    dc_pred_high[0][0][TX_8X16] = eb_aom_highbd_dc_128_predictor_8x16;
    dc_pred_high[0][0][TX_8X32] = eb_aom_highbd_dc_128_predictor_8x32;

    dc_pred_high[0][0][TX_16X4] = eb_aom_highbd_dc_128_predictor_16x4;
    dc_pred_high[0][0][TX_16X8] = eb_aom_highbd_dc_128_predictor_16x8;
    dc_pred_high[0][0][TX_16X32] = eb_aom_highbd_dc_128_predictor_16x32;
    dc_pred_high[0][0][TX_16X64] = eb_aom_highbd_dc_128_predictor_16x64;

    dc_pred_high[0][0][TX_32X8] = eb_aom_highbd_dc_128_predictor_32x8;
    dc_pred_high[0][0][TX_32X16] = eb_aom_highbd_dc_128_predictor_32x16;
    dc_pred_high[0][0][TX_32X64] = eb_aom_highbd_dc_128_predictor_32x64;

    dc_pred_high[0][0][TX_64X16] = eb_aom_highbd_dc_128_predictor_64x16;
    dc_pred_high[0][0][TX_64X32] = eb_aom_highbd_dc_128_predictor_64x32;

    dc_pred_high[0][1][TX_4X4] = eb_aom_highbd_dc_top_predictor_4x4;
    dc_pred_high[0][1][TX_8X8] = eb_aom_highbd_dc_top_predictor_8x8;
    dc_pred_high[0][1][TX_16X16] = eb_aom_highbd_dc_top_predictor_16x16;
    dc_pred_high[0][1][TX_32X32] = eb_aom_highbd_dc_top_predictor_32x32;
    dc_pred_high[0][1][TX_64X64] = eb_aom_highbd_dc_top_predictor_64x64;

    dc_pred_high[0][1][TX_4X8] = eb_aom_highbd_dc_top_predictor_4x8;
    dc_pred_high[0][1][TX_4X16] = eb_aom_highbd_dc_top_predictor_4x16;

    dc_pred_high[0][1][TX_8X4] = eb_aom_highbd_dc_top_predictor_8x4;
    dc_pred_high[0][1][TX_8X16] = eb_aom_highbd_dc_top_predictor_8x16;
    dc_pred_high[0][1][TX_8X32] = eb_aom_highbd_dc_top_predictor_8x32;

    dc_pred_high[0][1][TX_16X4] = eb_aom_highbd_dc_top_predictor_16x4;
    dc_pred_high[0][1][TX_16X8] = eb_aom_highbd_dc_top_predictor_16x8;
    dc_pred_high[0][1][TX_16X32] = eb_aom_highbd_dc_top_predictor_16x32;
    dc_pred_high[0][1][TX_16X64] = eb_aom_highbd_dc_top_predictor_16x64;

    dc_pred_high[0][1][TX_32X8] = eb_aom_highbd_dc_top_predictor_32x8;
    dc_pred_high[0][1][TX_32X16] = eb_aom_highbd_dc_top_predictor_32x16;
    dc_pred_high[0][1][TX_32X64] = eb_aom_highbd_dc_top_predictor_32x64;

    dc_pred_high[0][1][TX_64X16] = eb_aom_highbd_dc_top_predictor_64x16;
    dc_pred_high[0][1][TX_64X32] = eb_aom_highbd_dc_top_predictor_64x32;

    dc_pred_high[1][0][TX_4X4] = eb_aom_highbd_dc_left_predictor_4x4;
    dc_pred_high[1][0][TX_8X8] = eb_aom_highbd_dc_left_predictor_8x8;
    dc_pred_high[1][0][TX_16X16] = eb_aom_highbd_dc_left_predictor_16x16;
    dc_pred_high[1][0][TX_32X32] = eb_aom_highbd_dc_left_predictor_32x32;
    dc_pred_high[1][0][TX_64X64] = eb_aom_highbd_dc_left_predictor_64x64;

    dc_pred_high[1][0][TX_4X8] = eb_aom_highbd_dc_left_predictor_4x8;
    dc_pred_high[1][0][TX_4X16] = eb_aom_highbd_dc_left_predictor_4x16;

    dc_pred_high[1][0][TX_8X4] = eb_aom_highbd_dc_left_predictor_8x4;
    dc_pred_high[1][0][TX_8X16] = eb_aom_highbd_dc_left_predictor_8x16;
    dc_pred_high[1][0][TX_8X32] = eb_aom_highbd_dc_left_predictor_8x32;

    dc_pred_high[1][0][TX_16X4] = eb_aom_highbd_dc_left_predictor_16x4;
    dc_pred_high[1][0][TX_16X8] = eb_aom_highbd_dc_left_predictor_16x8;
    dc_pred_high[1][0][TX_16X32] = eb_aom_highbd_dc_left_predictor_16x32;
    dc_pred_high[1][0][TX_16X64] = eb_aom_highbd_dc_left_predictor_16x64;

    dc_pred_high[1][0][TX_32X8] = eb_aom_highbd_dc_left_predictor_32x8;
    dc_pred_high[1][0][TX_32X16] = eb_aom_highbd_dc_left_predictor_32x16;
    dc_pred_high[1][0][TX_32X64] = eb_aom_highbd_dc_left_predictor_32x64;

    dc_pred_high[1][0][TX_64X16] = eb_aom_highbd_dc_left_predictor_64x16;
    dc_pred_high[1][0][TX_64X32] = eb_aom_highbd_dc_left_predictor_64x32;

    dc_pred_high[1][1][TX_4X4] = eb_aom_highbd_dc_predictor_4x4;
    dc_pred_high[1][1][TX_8X8] = eb_aom_highbd_dc_predictor_8x8;
    dc_pred_high[1][1][TX_16X16] = eb_aom_highbd_dc_predictor_16x16;
    dc_pred_high[1][1][TX_32X32] = eb_aom_highbd_dc_predictor_32x32;
    dc_pred_high[1][1][TX_64X64] = eb_aom_highbd_dc_predictor_64x64;

    dc_pred_high[1][1][TX_4X8] = eb_aom_highbd_dc_predictor_4x8;
    dc_pred_high[1][1][TX_4X16] = eb_aom_highbd_dc_predictor_4x16;

    dc_pred_high[1][1][TX_8X4] = eb_aom_highbd_dc_predictor_8x4;
    dc_pred_high[1][1][TX_8X16] = eb_aom_highbd_dc_predictor_8x16;
    dc_pred_high[1][1][TX_8X32] = eb_aom_highbd_dc_predictor_8x32;

    dc_pred_high[1][1][TX_16X4] = eb_aom_highbd_dc_predictor_16x4;
    dc_pred_high[1][1][TX_16X8] = eb_aom_highbd_dc_predictor_16x8;
    dc_pred_high[1][1][TX_16X32] = eb_aom_highbd_dc_predictor_16x32;
    dc_pred_high[1][1][TX_16X64] = eb_aom_highbd_dc_predictor_16x64;

    dc_pred_high[1][1][TX_32X8] = eb_aom_highbd_dc_predictor_32x8;
    dc_pred_high[1][1][TX_32X16] = eb_aom_highbd_dc_predictor_32x16;
    dc_pred_high[1][1][TX_32X64] = eb_aom_highbd_dc_predictor_32x64;

    dc_pred_high[1][1][TX_64X16] = eb_aom_highbd_dc_predictor_64x16;
    dc_pred_high[1][1][TX_64X32] = eb_aom_highbd_dc_predictor_64x32;
}
void dr_predictor(uint8_t *dst, ptrdiff_t stride, TxSize tx_size,
    const uint8_t *above, const uint8_t *left,
    int32_t upsample_above, int32_t upsample_left, int32_t angle) {
    const int32_t dx = get_dx(angle);
    const int32_t dy = get_dy(angle);
    const int32_t bw = tx_size_wide[tx_size];
    const int32_t bh = tx_size_high[tx_size];
    assert(angle > 0 && angle < 270);

    if (angle > 0 && angle < 90) {
        eb_av1_dr_prediction_z1(dst, stride, bw, bh, above, left, upsample_above, dx,
            dy);
    }
    else if (angle > 90 && angle < 180) {
        eb_av1_dr_prediction_z2(dst, stride, bw, bh, above, left, upsample_above,
            upsample_left, dx, dy);
    }
    else if (angle > 180 && angle < 270) {
        eb_av1_dr_prediction_z3(dst, stride, bw, bh, above, left, upsample_left, dx,
            dy);
    }
    else if (angle == 90)
        pred[V_PRED][tx_size](dst, stride, above, left);
    else if (angle == 180)
        pred[H_PRED][tx_size](dst, stride, above, left);
}

void filter_intra_edge_corner(uint8_t *p_above, uint8_t *p_left) {
    const int32_t kernel[3] = { 5, 6, 5 };

    int32_t s = (p_left[0] * kernel[0]) + (p_above[-1] * kernel[1]) +
        (p_above[0] * kernel[2]);
    s = (s + 8) >> 4;
    p_above[-1] = (uint8_t)s;
    p_left[-1] = (uint8_t)s;
}

// Directional prediction, zone 1: 0 < angle < 90
void eb_av1_highbd_dr_prediction_z1_c(uint16_t *dst, ptrdiff_t stride, int32_t bw,
    int32_t bh, const uint16_t *above,
    const uint16_t *left, int32_t upsample_above,
    int32_t dx, int32_t dy, int32_t bd) {
    int32_t r, c, x, base, shift, val;

    (void)left;
    (void)dy;
    (void)bd;
    assert(dy == 1);
    assert(dx > 0);

    const int32_t max_base_x = ((bw + bh) - 1) << upsample_above;
    const int32_t frac_bits = 6 - upsample_above;
    const int32_t base_inc = 1 << upsample_above;
    x = dx;
    for (r = 0; r < bh; ++r, dst += stride, x += dx) {
        base = x >> frac_bits;
        shift = ((x << upsample_above) & 0x3F) >> 1;

        if (base >= max_base_x) {
            for (int32_t i = r; i < bh; ++i) {
                eb_aom_memset16(dst, above[max_base_x], bw);
                dst += stride;
            }
            return;
        }

        for (c = 0; c < bw; ++c, base += base_inc) {
            if (base < max_base_x) {
                val = above[base] * (32 - shift) + above[base + 1] * shift;
                val = ROUND_POWER_OF_TWO(val, 5);
                dst[c] = (uint16_t)clip_pixel_highbd(val, bd);
            }
            else
                dst[c] = above[max_base_x];
        }
    }
}

// Directional prediction, zone 2: 90 < angle < 180
void eb_av1_highbd_dr_prediction_z2_c(uint16_t *dst, ptrdiff_t stride, int32_t bw,
    int32_t bh, const uint16_t *above,
    const uint16_t *left, int32_t upsample_above,
    int32_t upsample_left, int32_t dx, int32_t dy, int32_t bd) {
    int32_t r, c, x, y, shift, val, base;

    (void)bd;
    assert(dx > 0);
    assert(dy > 0);

    const int32_t min_base_x = -(1 << upsample_above);
    const int32_t frac_bits_x = 6 - upsample_above;
    const int32_t frac_bits_y = 6 - upsample_left;
    for (r = 0; r < bh; ++r) {
        for (c = 0; c < bw; ++c) {
            y = r + 1;
            x = (c << 6) - y * dx;
            base = x >> frac_bits_x;
            if (base >= min_base_x) {
                shift = ((x * (1 << upsample_above)) & 0x3F) >> 1;
                val = above[base] * (32 - shift) + above[base + 1] * shift;
                val = ROUND_POWER_OF_TWO(val, 5);
            }
            else {
                x = c + 1;
                y = (r << 6) - x * dy;
                base = y >> frac_bits_y;
                shift = ((y * (1 << upsample_left)) & 0x3F) >> 1;
                val = left[base] * (32 - shift) + left[base + 1] * shift;
                val = ROUND_POWER_OF_TWO(val, 5);
            }
            dst[c] = (uint16_t)clip_pixel_highbd(val, bd);
        }
        dst += stride;
    }
}

// Directional prediction, zone 3: 180 < angle < 270
void eb_av1_highbd_dr_prediction_z3_c(uint16_t *dst, ptrdiff_t stride, int32_t bw,
    int32_t bh, const uint16_t *above,
    const uint16_t *left, int32_t upsample_left,
    int32_t dx, int32_t dy, int32_t bd) {
    int32_t r, c, y, base, shift, val;

    (void)above;
    (void)dx;
    (void)bd;
    assert(dx == 1);
    assert(dy > 0);

    const int32_t max_base_y = (bw + bh - 1) << upsample_left;
    const int32_t frac_bits = 6 - upsample_left;
    const int32_t base_inc = 1 << upsample_left;
    y = dy;
    for (c = 0; c < bw; ++c, y += dy) {
        base = y >> frac_bits;
        shift = ((y << upsample_left) & 0x3F) >> 1;

        for (r = 0; r < bh; ++r, base += base_inc) {
            if (base < max_base_y) {
                val = left[base] * (32 - shift) + left[base + 1] * shift;
                val = ROUND_POWER_OF_TWO(val, 5);
                dst[r * stride + c] = (uint16_t)clip_pixel_highbd(val, bd);
            }
            else {
                for (; r < bh; ++r) dst[r * stride + c] = left[max_base_y];
                break;
            }
        }
    }
}

void highbd_dr_predictor(uint16_t *dst, ptrdiff_t stride,
    TxSize tx_size, const uint16_t *above,
    const uint16_t *left, int32_t upsample_above,
    int32_t upsample_left, int32_t angle, int32_t bd) {
    const int32_t dx = get_dx(angle);
    const int32_t dy = get_dy(angle);
    const int32_t bw = tx_size_wide[tx_size];
    const int32_t bh = tx_size_high[tx_size];
    assert(angle > 0 && angle < 270);

    if (angle > 0 && angle < 90) {
        eb_av1_highbd_dr_prediction_z1(dst, stride, bw, bh, above, left,
            upsample_above, dx, dy, bd);
    }
    else if (angle > 90 && angle < 180) {
        eb_av1_highbd_dr_prediction_z2(dst, stride, bw, bh, above, left,
            upsample_above, upsample_left, dx, dy, bd);
    }
    else if (angle > 180 && angle < 270) {
        eb_av1_highbd_dr_prediction_z3(dst, stride, bw, bh, above, left, upsample_left,
            dx, dy, bd);
    }
    else if (angle == 90)
        pred_high[V_PRED][tx_size](dst, stride, above, left, bd);
    else if (angle == 180)
        pred_high[H_PRED][tx_size](dst, stride, above, left, bd);
}

void eb_av1_filter_intra_edge_high_c(uint16_t *p, int32_t sz, int32_t strength) {
    if (!strength) return;

    const int32_t kernel[INTRA_EDGE_FILT][INTRA_EDGE_TAPS] = {
        { 0, 4, 8, 4, 0 }, { 0, 5, 6, 5, 0 }, { 2, 4, 4, 4, 2 }
    };
    const int32_t filt = strength - 1;
    uint16_t edge[129];

    memcpy(edge, p, sz * sizeof(*p));
    for (int32_t i = 1; i < sz; i++) {
        int32_t s = 0;
        for (int32_t j = 0; j < INTRA_EDGE_TAPS; j++) {
            int32_t k = i - 2 + j;
            k = (k < 0) ? 0 : k;
            k = (k > sz - 1) ? sz - 1 : k;
            s += edge[k] * kernel[filt][j];
        }
        s = (s + 8) >> 4;
        p[i] = (uint16_t)s;
    }
}

void filter_intra_edge_corner_high(uint16_t *p_above, uint16_t *p_left) {
    const int32_t kernel[3] = { 5, 6, 5 };

    int32_t s = (p_left[0] * kernel[0]) + (p_above[-1] * kernel[1]) +
        (p_above[0] * kernel[2]);
    s = (s + 8) >> 4;
    p_above[-1] = (uint16_t)s;
    p_left[-1] = (uint16_t)s;
}

void eb_av1_upsample_intra_edge_high_c(uint16_t *p, int32_t sz, int32_t bd) {
    // interpolate half-sample positions
    assert(sz <= MAX_UPSAMPLE_SZ);

    uint16_t in[MAX_UPSAMPLE_SZ + 3];
    // copy p[-1..(sz-1)] and extend first and last samples
    in[0] = p[-1];
    in[1] = p[-1];
    for (int32_t i = 0; i < sz; i++)
        in[i + 2] = p[i];
    in[sz + 2] = p[sz - 1];

    // interpolate half-sample edge positions
    p[-2] = in[0];
    for (int32_t i = 0; i < sz; i++) {
        int32_t s = -in[i] + (9 * in[i + 1]) + (9 * in[i + 2]) - in[i + 3];
        s = (s + 8) >> 4;
        s = clip_pixel_highbd(s, bd);
        p[2 * i - 1] = (uint16_t)s;
        p[2 * i] = in[i + 2];
    }
}

void eb_av1_upsample_intra_edge_c(uint8_t *p, int32_t sz) {
    // interpolate half-sample positions
    assert(sz <= MAX_UPSAMPLE_SZ);

    uint8_t in[MAX_UPSAMPLE_SZ + 3];
    // copy p[-1..(sz-1)] and extend first and last samples
    in[0] = p[-1];
    in[1] = p[-1];
    for (int32_t i = 0; i < sz; i++)
        in[i + 2] = p[i];
    in[sz + 2] = p[sz - 1];

    // interpolate half-sample edge positions
    p[-2] = in[0];
    for (int32_t i = 0; i < sz; i++) {
        int32_t s = -in[i] + (9 * in[i + 1]) + (9 * in[i + 2]) - in[i + 3];
        s = clip_pixel((s + 8) >> 4);
        p[2 * i - 1] = (uint8_t)s;
        p[2 * i] = in[i + 2];
    }
}
/*static INLINE*/ BlockSize scale_chroma_bsize(BlockSize bsize, int32_t subsampling_x,
    int32_t subsampling_y) {
    BlockSize bs = bsize;
    switch (bsize) {
    case BLOCK_4X4:
        if (subsampling_x == 1 && subsampling_y == 1)
            bs = BLOCK_8X8;
        else if (subsampling_x == 1)
            bs = BLOCK_8X4;
        else if (subsampling_y == 1)
            bs = BLOCK_4X8;
        break;
    case BLOCK_4X8:
        if (subsampling_x == 1 && subsampling_y == 1)
            bs = BLOCK_8X8;
        else if (subsampling_x == 1)
            bs = BLOCK_8X8;
        else if (subsampling_y == 1)
            bs = BLOCK_4X8;
        break;
    case BLOCK_8X4:
        if (subsampling_x == 1 && subsampling_y == 1)
            bs = BLOCK_8X8;
        else if (subsampling_x == 1)
            bs = BLOCK_8X4;
        else if (subsampling_y == 1)
            bs = BLOCK_8X8;
        break;
    case BLOCK_4X16:
        if (subsampling_x == 1 && subsampling_y == 1)
            bs = BLOCK_8X16;
        else if (subsampling_x == 1)
            bs = BLOCK_8X16;
        else if (subsampling_y == 1)
            bs = BLOCK_4X16;
        break;
    case BLOCK_16X4:
        if (subsampling_x == 1 && subsampling_y == 1)
            bs = BLOCK_16X8;
        else if (subsampling_x == 1)
            bs = BLOCK_16X4;
        else if (subsampling_y == 1)
            bs = BLOCK_16X8;
        break;
    default: break;
    }
    return bs;
}

////////////########...........Recurssive intra prediction starting...........#########

DECLARE_ALIGNED(16, const int8_t,
                eb_av1_filter_intra_taps[FILTER_INTRA_MODES][8][8]) = {
  {
      { -6, 10, 0, 0, 0, 12, 0, 0 },
      { -5, 2, 10, 0, 0, 9, 0, 0 },
      { -3, 1, 1, 10, 0, 7, 0, 0 },
      { -3, 1, 1, 2, 10, 5, 0, 0 },
      { -4, 6, 0, 0, 0, 2, 12, 0 },
      { -3, 2, 6, 0, 0, 2, 9, 0 },
      { -3, 2, 2, 6, 0, 2, 7, 0 },
      { -3, 1, 2, 2, 6, 3, 5, 0 },
  },
  {
      { -10, 16, 0, 0, 0, 10, 0, 0 },
      { -6, 0, 16, 0, 0, 6, 0, 0 },
      { -4, 0, 0, 16, 0, 4, 0, 0 },
      { -2, 0, 0, 0, 16, 2, 0, 0 },
      { -10, 16, 0, 0, 0, 0, 10, 0 },
      { -6, 0, 16, 0, 0, 0, 6, 0 },
      { -4, 0, 0, 16, 0, 0, 4, 0 },
      { -2, 0, 0, 0, 16, 0, 2, 0 },
  },
  {
      { -8, 8, 0, 0, 0, 16, 0, 0 },
      { -8, 0, 8, 0, 0, 16, 0, 0 },
      { -8, 0, 0, 8, 0, 16, 0, 0 },
      { -8, 0, 0, 0, 8, 16, 0, 0 },
      { -4, 4, 0, 0, 0, 0, 16, 0 },
      { -4, 0, 4, 0, 0, 0, 16, 0 },
      { -4, 0, 0, 4, 0, 0, 16, 0 },
      { -4, 0, 0, 0, 4, 0, 16, 0 },
  },
  {
      { -2, 8, 0, 0, 0, 10, 0, 0 },
      { -1, 3, 8, 0, 0, 6, 0, 0 },
      { -1, 2, 3, 8, 0, 4, 0, 0 },
      { 0, 1, 2, 3, 8, 2, 0, 0 },
      { -1, 4, 0, 0, 0, 3, 10, 0 },
      { -1, 3, 4, 0, 0, 4, 6, 0 },
      { -1, 2, 3, 4, 0, 4, 4, 0 },
      { -1, 2, 2, 3, 4, 3, 3, 0 },
  },
  {
      { -12, 14, 0, 0, 0, 14, 0, 0 },
      { -10, 0, 14, 0, 0, 12, 0, 0 },
      { -9, 0, 0, 14, 0, 11, 0, 0 },
      { -8, 0, 0, 0, 14, 10, 0, 0 },
      { -10, 12, 0, 0, 0, 0, 14, 0 },
      { -9, 1, 12, 0, 0, 0, 12, 0 },
      { -8, 0, 0, 12, 0, 1, 11, 0 },
      { -7, 0, 0, 1, 12, 1, 9, 0 },
  },
};

void eb_av1_filter_intra_predictor_c(uint8_t *dst, ptrdiff_t stride,
                                  TxSize tx_size,
                                  const uint8_t *above,
                                  const uint8_t *left, int32_t mode) {
  int r, c;
  uint8_t buffer[33][33];
  const int bw = tx_size_wide[tx_size];
  const int bh = tx_size_high[tx_size];

  assert(bw <= 32 && bh <= 32);

  // The initialization is just for silencing Jenkins static analysis warnings
  for (r = 0; r < bh + 1; ++r)
    memset(buffer[r], 0, (bw + 1) * sizeof(buffer[0][0]));

  for (r = 0; r < bh; ++r) buffer[r + 1][0] = left[r];
  memcpy(buffer[0], &above[-1], (bw + 1) * sizeof(uint8_t));

  for (r = 1; r < bh + 1; r += 2)
    for (c = 1; c < bw + 1; c += 4) {
      const uint8_t p0 = buffer[r - 1][c - 1];
      const uint8_t p1 = buffer[r - 1][c];
      const uint8_t p2 = buffer[r - 1][c + 1];
      const uint8_t p3 = buffer[r - 1][c + 2];
      const uint8_t p4 = buffer[r - 1][c + 3];
      const uint8_t p5 = buffer[r][c - 1];
      const uint8_t p6 = buffer[r + 1][c - 1];
      for (int k = 0; k < 8; ++k) {
        int r_offset = k >> 2;
        int c_offset = k & 0x03;
        buffer[r + r_offset][c + c_offset] =
            clip_pixel(ROUND_POWER_OF_TWO_SIGNED(
                eb_av1_filter_intra_taps[mode][k][0] * p0 +
                    eb_av1_filter_intra_taps[mode][k][1] * p1 +
                    eb_av1_filter_intra_taps[mode][k][2] * p2 +
                    eb_av1_filter_intra_taps[mode][k][3] * p3 +
                    eb_av1_filter_intra_taps[mode][k][4] * p4 +
                    eb_av1_filter_intra_taps[mode][k][5] * p5 +
                    eb_av1_filter_intra_taps[mode][k][6] * p6,
                FILTER_INTRA_SCALE_BITS));
      }
    }

  for (r = 0; r < bh; ++r) {
    memcpy(dst, &buffer[r + 1][1], bw * sizeof(uint8_t));
    dst += stride;
  }
}

 void highbd_filter_intra_predictor(uint16_t *dst, ptrdiff_t stride,
                                          TxSize tx_size,
                                          const uint16_t *above,
                                          const uint16_t *left, int mode,
                                          int bd) {
    int r, c;
    uint16_t buffer[33][33];
    const int bw = tx_size_wide[tx_size];
    const int bh = tx_size_high[tx_size];

    assert(bw <= 32 && bh <= 32);

    // The initialization is just for silencing Jenkins static analysis warnings
    for (r = 0; r < bh + 1; ++r)
        memset(buffer[r], 0, (bw + 1) * sizeof(buffer[0][0]));

    for (r = 0; r < bh; ++r) buffer[r + 1][0] = left[r];
    memcpy(buffer[0], &above[-1], (bw + 1) * sizeof(buffer[0][0]));

    for (r = 1; r < bh + 1; r += 2)
        for (c = 1; c < bw + 1; c += 4) {
            const uint16_t p0 = buffer[r - 1][c - 1];
            const uint16_t p1 = buffer[r - 1][c];
            const uint16_t p2 = buffer[r - 1][c + 1];
            const uint16_t p3 = buffer[r - 1][c + 2];
            const uint16_t p4 = buffer[r - 1][c + 3];
            const uint16_t p5 = buffer[r][c - 1];
            const uint16_t p6 = buffer[r + 1][c - 1];
            for (int k = 0; k < 8; ++k) {
                int r_offset = k >> 2;
                int c_offset = k & 0x03;
                buffer[r + r_offset][c + c_offset] =
                    clip_pixel_highbd(ROUND_POWER_OF_TWO_SIGNED(
                        eb_av1_filter_intra_taps[mode][k][0] * p0 +
                        eb_av1_filter_intra_taps[mode][k][1] * p1 +
                        eb_av1_filter_intra_taps[mode][k][2] * p2 +
                        eb_av1_filter_intra_taps[mode][k][3] * p3 +
                        eb_av1_filter_intra_taps[mode][k][4] * p4 +
                        eb_av1_filter_intra_taps[mode][k][5] * p5 +
                        eb_av1_filter_intra_taps[mode][k][6] * p6,
                        FILTER_INTRA_SCALE_BITS),
                        bd);
            }
        }

    for (r = 0; r < bh; ++r) {
        memcpy(dst, &buffer[r + 1][1], bw * sizeof(dst[0]));
        dst += stride;
    }
}

////////////#####################...........Recurssive intra prediction ending...........#####################////////////

static void build_intra_predictors(
    const MacroBlockD *xd,
    uint8_t* topNeighArray,
    uint8_t* leftNeighArray,
    // const uint8_t *ref,    int32_t ref_stride,
    uint8_t *dst, int32_t dst_stride,
    PredictionMode mode, int32_t angle_delta,
    FilterIntraMode filter_intra_mode,
    TxSize tx_size, int32_t disable_edge_filter,
    int32_t n_top_px, int32_t n_topright_px,
    int32_t n_left_px, int32_t n_bottomleft_px,
    int32_t plane)
{
    (void)xd;
    int32_t i;

    int32_t ref_stride = 1;
    const uint8_t *above_ref = topNeighArray;//CHKN ref - ref_stride;
    const uint8_t *left_ref = leftNeighArray;//CHKN ref - 1;
    DECLARE_ALIGNED(32, uint8_t, left_data[MAX_TX_SIZE * 2 + 48]);
    DECLARE_ALIGNED(32, uint8_t, above_data[MAX_TX_SIZE * 2 + 48]);
    uint8_t *const above_row = above_data + 32;
    uint8_t *const left_col = left_data + 32;

    const int32_t txwpx = tx_size_wide[tx_size];
    const int32_t txhpx = tx_size_high[tx_size];
    int32_t need_left = extend_modes[mode] & NEED_LEFT;
    int32_t need_above = extend_modes[mode] & NEED_ABOVE;
    int32_t need_above_left = extend_modes[mode] & NEED_ABOVELEFT;
    int32_t p_angle = 0;
    const int32_t is_dr_mode = av1_is_directional_mode(mode);
    const int32_t use_filter_intra = filter_intra_mode != FILTER_INTRA_MODES;

    if (is_dr_mode) {
        p_angle = mode_to_angle_map[mode] + angle_delta * ANGLE_STEP;
        if (p_angle <= 90)
            need_above = 1, need_left = 0, need_above_left = 1;
        else if (p_angle < 180)
            need_above = 1, need_left = 1, need_above_left = 1;
        else
            need_above = 0, need_left = 1, need_above_left = 1;
    }
    if (use_filter_intra) need_left = need_above = need_above_left = 1;

    assert(n_top_px >= 0);
    assert(n_topright_px >= 0);
    assert(n_left_px >= 0);
    assert(n_bottomleft_px >= 0);

    if ((!need_above && n_left_px == 0) || (!need_left && n_top_px == 0)) {
        int32_t val;
        if (need_left)
            val = (n_top_px > 0) ? above_ref[0] : 129;
        else
            val = (n_left_px > 0) ? left_ref[0] : 127;
        for (i = 0; i < txhpx; ++i) {
            memset(dst, val, txwpx);
            dst += dst_stride;
        }
        return;
    }

    // NEED_LEFT
    if (need_left) {
        int32_t need_bottom = !!(extend_modes[mode] & NEED_BOTTOMLEFT);
        if (use_filter_intra) need_bottom = 0;
        if (is_dr_mode) need_bottom = p_angle > 180;
        const int32_t num_left_pixels_needed = txhpx + (need_bottom ? txwpx : 0);
        i = 0;
        if (n_left_px > 0) {
            for (; i < n_left_px; i++) left_col[i] = left_ref[i * ref_stride];
            if (need_bottom && n_bottomleft_px > 0) {
                assert(i == txhpx);
                for (; i < txhpx + n_bottomleft_px; i++)
                    left_col[i] = left_ref[i * ref_stride];
            }
            if (i < num_left_pixels_needed)
                memset(&left_col[i], left_col[i - 1], num_left_pixels_needed - i);
        }
        else {
            if (n_top_px > 0)
                memset(left_col, above_ref[0], num_left_pixels_needed);
            else
                memset(left_col, 129, num_left_pixels_needed);
        }
    }

    // NEED_ABOVE
    if (need_above) {
        int32_t need_right = !!(extend_modes[mode] & NEED_ABOVERIGHT);
        if (use_filter_intra) need_right = 0;
        if (is_dr_mode) need_right = p_angle < 90;
        const int32_t num_top_pixels_needed = txwpx + (need_right ? txhpx : 0);
        if (n_top_px > 0) {
            memcpy(above_row, above_ref, n_top_px);
            i = n_top_px;
            if (need_right && n_topright_px > 0) {
                assert(n_top_px == txwpx);
                memcpy(above_row + txwpx, above_ref + txwpx, n_topright_px);
                i += n_topright_px;
            }
            if (i < num_top_pixels_needed)
                memset(&above_row[i], above_row[i - 1], num_top_pixels_needed - i);
        }
        else {
            if (n_left_px > 0)
                memset(above_row, left_ref[0], num_top_pixels_needed);
            else
                memset(above_row, 127, num_top_pixels_needed);
        }
    }

    if (need_above_left) {
        if (n_top_px > 0 && n_left_px > 0)
            above_row[-1] = above_ref[-1];
        else if (n_top_px > 0)
            above_row[-1] = above_ref[0];
        else if (n_left_px > 0)
            above_row[-1] = left_ref[0];
        else
            above_row[-1] = 128;
        left_col[-1] = above_row[-1];
    }
#if FILTER_INTRA_FLAG
  if (use_filter_intra) {
    eb_av1_filter_intra_predictor(dst, dst_stride, tx_size, above_row, left_col,
                               filter_intra_mode);
    return;
  }
#else
    //    if (use_filter_intra) {
    ////        eb_av1_filter_intra_predictor(dst, dst_stride, tx_size, above_row, left_col,
    ////CHKN            filter_intra_mode);
    //        return;
    //    }
#endif
    if (is_dr_mode) {
        int32_t upsample_above = 0;
        int32_t upsample_left = 0;
        if (!disable_edge_filter) {
            const int32_t need_right = p_angle < 90;
            const int32_t need_bottom = p_angle > 180;
            const int32_t filt_type = get_filt_type(xd, plane);

            if (p_angle != 90 && p_angle != 180) {
                const int32_t ab_le = need_above_left ? 1 : 0;
                if (need_above && need_left && (txwpx + txhpx >= 24))
                    filter_intra_edge_corner(above_row, left_col);
                if (need_above && n_top_px > 0) {
                    const int32_t strength =
                        intra_edge_filter_strength(txwpx, txhpx, p_angle - 90, filt_type);
                    const int32_t n_px = n_top_px + ab_le + (need_right ? txhpx : 0);
                    eb_av1_filter_intra_edge(above_row - ab_le, n_px, strength);
                }
                if (need_left && n_left_px > 0) {
                    const int32_t strength = intra_edge_filter_strength(
                        txhpx, txwpx, p_angle - 180, filt_type);
                    const int32_t n_px = n_left_px + ab_le + (need_bottom ? txwpx : 0);
                    eb_av1_filter_intra_edge(left_col - ab_le, n_px, strength);
                }
            }
            upsample_above =
                use_intra_edge_upsample(txwpx, txhpx, p_angle - 90, filt_type);
            if (need_above && upsample_above) {
                const int32_t n_px = txwpx + (need_right ? txhpx : 0);
                eb_av1_upsample_intra_edge(above_row, n_px);
            }
            upsample_left =
                use_intra_edge_upsample(txhpx, txwpx, p_angle - 180, filt_type);
            if (need_left && upsample_left) {
                const int32_t n_px = txhpx + (need_bottom ? txwpx : 0);
                eb_av1_upsample_intra_edge(left_col, n_px);
            }
        }
        dr_predictor(dst, dst_stride, tx_size, above_row, left_col, upsample_above,
            upsample_left, p_angle);
        return;
    }

    // predict
    if (mode == DC_PRED) {
        dc_pred[n_left_px > 0][n_top_px > 0][tx_size](dst, dst_stride, above_row,
            left_col);
    }
    else
        pred[mode][tx_size](dst, dst_stride, above_row, left_col);
}
static void build_intra_predictors_high(
    const MacroBlockD *xd,
    uint16_t* topNeighArray, // int8_t
    uint16_t* leftNeighArray, // int8_t
    //const uint8_t *ref8, int32_t ref_stride,
    uint16_t *dst,//uint8_t *dst8
    int32_t dst_stride, PredictionMode mode, int32_t angle_delta,
    FilterIntraMode filter_intra_mode, TxSize tx_size,
    int32_t disable_edge_filter, int32_t n_top_px, int32_t n_topright_px, int32_t n_left_px,
    int32_t n_bottomleft_px, int32_t plane, int32_t bd) {
    (void)xd;
    int32_t i;
    //uint16_t *dst = CONVERT_TO_SHORTPTR(dst8);
    //uint16_t *ref = CONVERT_TO_SHORTPTR(ref8);

    DECLARE_ALIGNED(16, uint16_t, left_data[MAX_TX_SIZE * 2 + 32]);
    DECLARE_ALIGNED(16, uint16_t, above_data[MAX_TX_SIZE * 2 + 32]);
    uint16_t *const above_row = above_data + 16;
    uint16_t *const left_col = left_data + 16;
    const int32_t txwpx = tx_size_wide[tx_size];
    const int32_t txhpx = tx_size_high[tx_size];
    int32_t need_left = extend_modes[mode] & NEED_LEFT;
    int32_t need_above = extend_modes[mode] & NEED_ABOVE;
    int32_t need_above_left = extend_modes[mode] & NEED_ABOVELEFT;

    int32_t ref_stride = 1;
    const uint16_t *above_ref = topNeighArray;
    const uint16_t *left_ref = leftNeighArray;
    //const uint16_t *above_ref = ref - ref_stride;
    //const uint16_t *left_ref = ref - 1;
    int32_t p_angle = 0;
    const int32_t is_dr_mode = av1_is_directional_mode(mode);
    const int32_t use_filter_intra = filter_intra_mode != FILTER_INTRA_MODES;
    int32_t base = 128 << (bd - 8);

    // The default values if ref pixels are not available:
    // base-1 base-1 base-1 .. base-1 base-1 base-1 base-1 base-1 base-1
    // base+1   A      B  ..     Y      Z
    // base+1   C      D  ..     W      X
    // base+1   E      F  ..     U      V
    // base+1   G      H  ..     S      T      T      T      T      T

    if (is_dr_mode) {
        p_angle = mode_to_angle_map[mode] + angle_delta * ANGLE_STEP;
        if (p_angle <= 90)
            need_above = 1, need_left = 0, need_above_left = 1;
        else if (p_angle < 180)
            need_above = 1, need_left = 1, need_above_left = 1;
        else
            need_above = 0, need_left = 1, need_above_left = 1;
    }
    if (use_filter_intra) need_left = need_above = need_above_left = 1;

    assert(n_top_px >= 0);
    assert(n_topright_px >= 0);
    assert(n_left_px >= 0);
    assert(n_bottomleft_px >= 0);

    if ((!need_above && n_left_px == 0) || (!need_left && n_top_px == 0)) {
        int32_t val;
        if (need_left)
            val = (n_top_px > 0) ? above_ref[0] : base + 1;
        else
            val = (n_left_px > 0) ? left_ref[0] : base - 1;
        for (i = 0; i < txhpx; ++i) {
            eb_aom_memset16(dst, val, txwpx);
            dst += dst_stride;
        }
        return;
    }

    // NEED_LEFT
    if (need_left) {
        int32_t need_bottom = !!(extend_modes[mode] & NEED_BOTTOMLEFT);
        if (use_filter_intra) need_bottom = 0;
        if (is_dr_mode) need_bottom = p_angle > 180;
        const int32_t num_left_pixels_needed = txhpx + (need_bottom ? txwpx : 0);
        i = 0;
        if (n_left_px > 0) {
            for (; i < n_left_px; i++) left_col[i] = left_ref[i * ref_stride];
            if (need_bottom && n_bottomleft_px > 0) {
                assert(i == txhpx);
                for (; i < txhpx + n_bottomleft_px; i++)
                    left_col[i] = left_ref[i * ref_stride];
            }
            if (i < num_left_pixels_needed)
                eb_aom_memset16(&left_col[i], left_col[i - 1], num_left_pixels_needed - i);
        }
        else {
            if (n_top_px > 0)
                eb_aom_memset16(left_col, above_ref[0], num_left_pixels_needed);
            else
                eb_aom_memset16(left_col, base + 1, num_left_pixels_needed);
        }
    }

    // NEED_ABOVE
    if (need_above) {
        int32_t need_right = !!(extend_modes[mode] & NEED_ABOVERIGHT);
        if (use_filter_intra) need_right = 0;
        if (is_dr_mode) need_right = p_angle < 90;
        const int32_t num_top_pixels_needed = txwpx + (need_right ? txhpx : 0);
        if (n_top_px > 0) {
            memcpy(above_row, above_ref, n_top_px * sizeof(above_ref[0]));
            i = n_top_px;
            if (need_right && n_topright_px > 0) {
                assert(n_top_px == txwpx);
                memcpy(above_row + txwpx, above_ref + txwpx,
                    n_topright_px * sizeof(above_ref[0]));
                i += n_topright_px;
            }
            if (i < num_top_pixels_needed)
                eb_aom_memset16(&above_row[i], above_row[i - 1],
                    num_top_pixels_needed - i);
        }
        else {
            if (n_left_px > 0)
                eb_aom_memset16(above_row, left_ref[0], num_top_pixels_needed);
            else
                eb_aom_memset16(above_row, base - 1, num_top_pixels_needed);
        }
    }

    if (need_above_left) {
        if (n_top_px > 0 && n_left_px > 0)
            above_row[-1] = above_ref[-1];
        else if (n_top_px > 0)
            above_row[-1] = above_ref[0];
        else if (n_left_px > 0)
            above_row[-1] = left_ref[0];
        else
            above_row[-1] = (uint16_t)base;
        left_col[-1] = above_row[-1];
    }
#if FILTER_INTRA_FLAG
if (use_filter_intra) {
    highbd_filter_intra_predictor(dst, dst_stride, tx_size, above_row, left_col,
                               filter_intra_mode,10);
    return;
  }
#else
    // not added yet
    //if (use_filter_intra) {
    //    highbd_filter_intra_predictor(dst, dst_stride, tx_size, above_row, left_col,
    //        filter_intra_mode, xd->bd);
    //    return;
    //}
#endif
    if (is_dr_mode) {
        int32_t upsample_above = 0;
        int32_t upsample_left = 0;
        if (!disable_edge_filter) {
            const int32_t need_right = p_angle < 90;
            const int32_t need_bottom = p_angle > 180;
            //const int32_t filt_type = get_filt_type(xd, plane);
            const int32_t filt_type = get_filt_type(xd, plane);
            if (p_angle != 90 && p_angle != 180) {
                const int32_t ab_le = need_above_left ? 1 : 0;
                if (need_above && need_left && (txwpx + txhpx >= 24))
                    filter_intra_edge_corner_high(above_row, left_col);
                if (need_above && n_top_px > 0) {
                    const int32_t strength =
                        intra_edge_filter_strength(txwpx, txhpx, p_angle - 90, filt_type);
                    const int32_t n_px = n_top_px + ab_le + (need_right ? txhpx : 0);
                    eb_av1_filter_intra_edge_high(above_row - ab_le, n_px, strength);
                }
                if (need_left && n_left_px > 0) {
                    const int32_t strength = intra_edge_filter_strength(
                        txhpx, txwpx, p_angle - 180, filt_type);
                    const int32_t n_px = n_left_px + ab_le + (need_bottom ? txwpx : 0);

                    eb_av1_filter_intra_edge_high(left_col - ab_le, n_px, strength);
                }
            }
            upsample_above =
                use_intra_edge_upsample(txwpx, txhpx, p_angle - 90, filt_type);
            if (need_above && upsample_above) {
                const int32_t n_px = txwpx + (need_right ? txhpx : 0);
                //av1_upsample_intra_edge_high(above_row, n_px, bd);// AMIR : to be replaced by optimized code
                eb_av1_upsample_intra_edge_high_c(above_row, n_px, bd);
            }
            upsample_left =
                use_intra_edge_upsample(txhpx, txwpx, p_angle - 180, filt_type);
            if (need_left && upsample_left) {
                const int32_t n_px = txhpx + (need_bottom ? txwpx : 0);
                //av1_upsample_intra_edge_high(left_col, n_px, bd);// AMIR: to be replaced by optimized code
                eb_av1_upsample_intra_edge_high_c(left_col, n_px, bd);
            }
        }
        highbd_dr_predictor(dst, dst_stride, tx_size, above_row, left_col,
            upsample_above, upsample_left, p_angle, bd);
        return;
    }

    // predict
    if (mode == DC_PRED) {
        dc_pred_high[n_left_px > 0][n_top_px > 0][tx_size](
            dst, dst_stride, above_row, left_col, bd);
    }
    else
        pred_high[mode][tx_size](dst, dst_stride, above_row, left_col, bd);
}

void eb_av1_predict_intra_block(
    TileInfo * tile,
    STAGE       stage,
    const BlockGeom * blk_geom,
    const Av1Common *cm,
    int32_t wpx,
    int32_t hpx,
    TxSize tx_size,
    PredictionMode mode,
    int32_t angle_delta,
    int32_t use_palette,
#if PAL_SUP
    PaletteInfo  *palette_info,
#endif
    FilterIntraMode filter_intra_mode,
    uint8_t* topNeighArray,
    uint8_t* leftNeighArray,
    EbPictureBufferDesc  *recon_buffer,
    int32_t col_off,
    int32_t row_off,
    int32_t plane,
    BlockSize bsize,
    uint32_t tu_org_x_pict,
    uint32_t tu_org_y_pict,
    uint32_t bl_org_x_pict,
    uint32_t bl_org_y_pict,
    uint32_t bl_org_x_mb,
    uint32_t bl_org_y_mb)
{
    (void)use_palette;
    MacroBlockD xdS;
    MacroBlockD *xd = &xdS;

    uint32_t  pred_buf_x_offest;
    uint32_t  pred_buf_y_offest;

    if (stage == ED_STAGE) { // EncDec
        pred_buf_x_offest = plane ? ((bl_org_x_pict >> 3) << 3) >> 1 : tu_org_x_pict;
        pred_buf_y_offest = plane ? ((bl_org_y_pict >> 3) << 3) >> 1 : tu_org_y_pict;
    }
    else { // MD
        pred_buf_x_offest = bl_org_x_mb;
        pred_buf_y_offest = bl_org_y_mb;
    }

    // Adjust mirow , micol ;
    // All plane have the same values

    int32_t mirow = bl_org_y_pict >> 2;
    int32_t micol = bl_org_x_pict >> 2;
    xd->up_available   = (mirow > tile->mi_row_start);
    xd->left_available = (micol > tile->mi_col_start);
    const int32_t bw = mi_size_wide[bsize];
    const int32_t bh = mi_size_high[bsize];

    xd->mb_to_top_edge = -((mirow * MI_SIZE) * 8);
    xd->mb_to_bottom_edge = ((cm->mi_rows - bh - mirow) * MI_SIZE) * 8;
    xd->mb_to_left_edge = -((micol * MI_SIZE) * 8);
    xd->mb_to_right_edge = ((cm->mi_cols - bw - micol) * MI_SIZE) * 8;
    xd->tile.mi_col_start = tile->mi_col_start;
    xd->tile.mi_col_end = tile->mi_col_end;
    xd->tile.mi_row_start = tile->mi_row_start;
    xd->tile.mi_row_end = tile->mi_row_end;
    xd->n8_h = bh;
    xd->n8_w = bw;
    xd->is_sec_rect = 0;
    if (xd->n8_w < xd->n8_h) {
        // Only mark is_sec_rect as 1 for the last block.
        // For PARTITION_VERT_4, it would be (0, 0, 0, 1);
        // For other partitions, it would be (0, 1).
        if (!((micol + xd->n8_w) & (xd->n8_h - 1))) xd->is_sec_rect = 1;
    }

    if (xd->n8_w > xd->n8_h)
        if (mirow & (xd->n8_w - 1)) xd->is_sec_rect = 1;
    uint8_t  *dst;
    int32_t dst_stride;
    if (plane == 0) {
        dst = recon_buffer->buffer_y + pred_buf_x_offest + recon_buffer->origin_x + (pred_buf_y_offest + recon_buffer->origin_y)*recon_buffer->stride_y;
        dst_stride = recon_buffer->stride_y;
    }
    else if (plane == 1) {
        dst = recon_buffer->buffer_cb + (pred_buf_x_offest + recon_buffer->origin_x / 2 + (pred_buf_y_offest + recon_buffer->origin_y / 2)*recon_buffer->stride_cb);
        dst_stride = recon_buffer->stride_cb;
    }
    else {
        dst = recon_buffer->buffer_cr + (pred_buf_x_offest + recon_buffer->origin_x / 2 + (pred_buf_y_offest + recon_buffer->origin_y / 2)*recon_buffer->stride_cr);
        dst_stride = recon_buffer->stride_cr;
    }

    int32_t chroma_up_available = xd->up_available;
    int32_t chroma_left_available = xd->left_available;
    const int32_t ss_x = plane == 0 ? 0 : 1; //CHKN
    const int32_t ss_y = plane == 0 ? 0 : 1;

    if (ss_x && bw < mi_size_wide[BLOCK_8X8])
        chroma_left_available = (micol - 1) > tile->mi_col_start;
    if (ss_y && bh < mi_size_high[BLOCK_8X8])
        chroma_up_available = (mirow - 1) > tile->mi_row_start;

    int mi_stride = cm->mi_stride;
    const int32_t offset = mirow * mi_stride + micol;
    xd->mi = cm->pcs_ptr->mi_grid_base + offset;
    ModeInfo *miPtr = *xd->mi;

    if (xd->up_available) {
       // xd->above_mbmi = xd->mi[-xd->mi_stride].mbmi;
        xd->above_mbmi = &miPtr[-mi_stride].mbmi;
    }
    else
        xd->above_mbmi = NULL;
    if (xd->left_available) {
        //xd->left_mbmi = xd->mi[-1].mbmi;
        xd->left_mbmi = &miPtr[-1].mbmi;
    }
    else
        xd->left_mbmi = NULL;
    const int chroma_ref = ((mirow & 0x01) || !(bh & 0x01) || !ss_y) &&
        ((micol & 0x01) || !(bw & 0x01) || !ss_x);
    if (chroma_ref) {
        // To help calculate the "above" and "left" chroma blocks, note that the
        // current block may cover multiple luma blocks (eg, if partitioned into
        // 4x4 luma blocks).
        // First, find the top-left-most luma block covered by this chroma block

        ModeInfo *miPtr = xd->mi[-(mirow & ss_y) * mi_stride - (micol & ss_x)];

        // Then, we consider the luma region covered by the left or above 4x4 chroma
        // prediction. We want to point to the chroma reference block in that
        // region, which is the bottom-right-most mi unit.
        // This leads to the following offsets:
        MbModeInfo *chroma_above_mi =
            chroma_up_available ? &miPtr[-mi_stride + ss_x].mbmi : NULL;
        xd->chroma_above_mbmi = chroma_above_mi;

        MbModeInfo *chroma_left_mi =
            chroma_left_available ? &miPtr[ss_y * mi_stride - 1].mbmi : NULL;
        xd->chroma_left_mbmi = chroma_left_mi;
    }

    //CHKN  const MbModeInfo *const mbmi = xd->mi[0];
    const int32_t txwpx = tx_size_wide[tx_size];
    const int32_t txhpx = tx_size_high[tx_size];
    const int32_t x = col_off << tx_size_wide_log2[0];
    const int32_t y = row_off << tx_size_high_log2[0];

#if PAL_SUP
    if (use_palette) {
        int32_t r, c;

         const uint8_t *const map = palette_info->color_idx_map;
         const uint16_t *const palette =
             palette_info->pmi.palette_colors + plane * PALETTE_MAX_SIZE;
        for (r = 0; r < txhpx; ++r) {
            for (c = 0; c < txwpx; ++c) {
            dst[r * dst_stride + c] =
                (uint8_t)palette[map[(r + y) * wpx + c + x]];
            }
        }
        return;
    }
#else

    //if (use_palette) {
    //  int32_t r, c;
    //  const uint8_t *const map = xd->plane[plane != 0].color_index_map;
    //  const uint16_t *const palette =
    //      mbmi->palette_mode_info.palette_colors + plane * PALETTE_MAX_SIZE;
    //  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    //    uint16_t *dst16 = CONVERT_TO_SHORTPTR(dst);
    //    for (r = 0; r < txhpx; ++r) {
    //      for (c = 0; c < txwpx; ++c) {
    //        dst16[r * dst_stride + c] = palette[map[(r + y) * wpx + c + x]];
    //      }
    //    }
    //  } else {
    //    for (r = 0; r < txhpx; ++r) {
    //      for (c = 0; c < txwpx; ++c) {
    //        dst[r * dst_stride + c] =
    //            (uint8_t)palette[map[(r + y) * wpx + c + x]];
    //      }
    //    }
    //  }
    //  return;
    //}
#endif
    //CHKN BlockSize bsize = mbmi->sb_type;
    struct MacroblockdPlane  pd_s;
    struct MacroblockdPlane * pd = &pd_s;
    if (plane == 0)
        pd->subsampling_x = pd->subsampling_y = 0;
    else
        pd->subsampling_x = pd->subsampling_y = 1;
    const int32_t txw = tx_size_wide_unit[tx_size];
    const int32_t txh = tx_size_high_unit[tx_size];
    const int32_t have_top = row_off || (pd->subsampling_y ? /*xd->*/chroma_up_available
        : xd->up_available);
    const int32_t have_left =
        col_off ||
        (pd->subsampling_x ? /*xd->*/chroma_left_available : xd->left_available);
    const int32_t mi_row = -xd->mb_to_top_edge >> (3 + MI_SIZE_LOG2);
    const int32_t mi_col = -xd->mb_to_left_edge >> (3 + MI_SIZE_LOG2);
    const int32_t xr_chr_offset = 0;
    const int32_t yd_chr_offset = 0;

    // Distance between the right edge of this prediction block to
    // the frame right edge
    const int32_t xr = (xd->mb_to_right_edge >> (3 + pd->subsampling_x)) +
        (wpx - x - txwpx) - xr_chr_offset;
    // Distance between the bottom edge of this prediction block to
    // the frame bottom edge
    const int32_t yd = (xd->mb_to_bottom_edge >> (3 + pd->subsampling_y)) +
        (hpx - y - txhpx) - yd_chr_offset;
    const int32_t right_available =
        mi_col + ((col_off + txw) << pd->subsampling_x) < xd->tile.mi_col_end;
    const int32_t bottom_available =
        (yd > 0) &&
        (mi_row + ((row_off + txh) << pd->subsampling_y) < xd->tile.mi_row_end);

    const PartitionType partition = from_shape_to_part[blk_geom->shape]; //cu_ptr->part;// PARTITION_NONE;//CHKN this is good enough as the avail functions need to know if VERT part is used or not mbmi->partition;

    // force 4x4 chroma component block size.
    bsize = scale_chroma_bsize(bsize, pd->subsampling_x, pd->subsampling_y);

    const int32_t have_top_right = intra_has_top_right(
        cm->p_pcs_ptr->sequence_control_set_ptr->seq_header.sb_size, bsize,
        mi_row, mi_col, have_top, right_available, partition, tx_size,
        row_off, col_off, pd->subsampling_x, pd->subsampling_y);
    const int32_t have_bottom_left = intra_has_bottom_left(
        cm->p_pcs_ptr->sequence_control_set_ptr->seq_header.sb_size, bsize,
        mi_row, mi_col, bottom_available, have_left, partition,
        tx_size, row_off, col_off, pd->subsampling_x, pd->subsampling_y);

    const int32_t disable_edge_filter = 0;//CHKN !cm->seq_params.enable_intra_edge_filter;

    //if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    //  build_intra_predictors_high(
    //      xd, ref, ref_stride, dst, dst_stride, mode, angle_delta,
    //      filter_intra_mode, tx_size, disable_edge_filter,
    //      have_top ? AOMMIN(txwpx, xr + txwpx) : 0,
    //      have_top_right ? AOMMIN(txwpx, xr) : 0,
    //      have_left ? AOMMIN(txhpx, yd + txhpx) : 0,
    //      have_bottom_left ? AOMMIN(txhpx, yd) : 0, plane);
    //  return;
    //}

    build_intra_predictors(
        xd,
        topNeighArray,
        leftNeighArray,
        // ref, ref_stride,
        dst, dst_stride, mode,
        angle_delta, filter_intra_mode, tx_size,
        disable_edge_filter,
        have_top ? AOMMIN(txwpx, xr + txwpx) : 0,
        have_top_right ? AOMMIN(txwpx, xr) : 0,
        have_left ? AOMMIN(txhpx, yd + txhpx) : 0,
        have_bottom_left ? AOMMIN(txhpx, yd) : 0, plane);
}

void eb_av1_predict_intra_block_16bit(
    TileInfo * tile,
    STAGE       stage,
    const BlockGeom * blk_geom,
    const Av1Common *cm,
    int32_t wpx,
    int32_t hpx,
    TxSize tx_size,
    PredictionMode mode,
    int32_t angle_delta,
    int32_t use_palette,
#if PAL_SUP
    PaletteInfo  *palette_info,
#endif
    FilterIntraMode filter_intra_mode,
    uint16_t* topNeighArray,
    uint16_t* leftNeighArray,
    EbPictureBufferDesc  *recon_buffer,
    int32_t col_off,
    int32_t row_off,
    int32_t plane,
    BlockSize bsize,
    uint32_t tu_org_x_pict,
    uint32_t tu_org_y_pict,
    uint32_t bl_org_x_pict,
    uint32_t bl_org_y_pict,
    uint32_t bl_org_x_mb,
    uint32_t bl_org_y_mb)
{
    (void)use_palette;
    MacroBlockD xdS;
    MacroBlockD *xd = &xdS;

    uint32_t  pred_buf_x_offest;
    uint32_t  pred_buf_y_offest;

    if (stage == ED_STAGE) { // EncDec
        pred_buf_x_offest = plane ? ((bl_org_x_pict >> 3) << 3) >> 1 : tu_org_x_pict;
        pred_buf_y_offest = plane ? ((bl_org_y_pict >> 3) << 3) >> 1 : tu_org_y_pict;
    } else { // MD
        pred_buf_x_offest = bl_org_x_mb;
        pred_buf_y_offest = bl_org_y_mb;
    }

    int32_t mirow = bl_org_y_pict >> 2;
    int32_t micol = bl_org_x_pict >> 2;
    xd->up_available = (mirow > tile->mi_row_start);
    xd->left_available = (micol > tile->mi_col_start);

    const int32_t bw = mi_size_wide[bsize];
    const int32_t bh = mi_size_high[bsize];

    xd->mb_to_top_edge = -((mirow * MI_SIZE) * 8);
    xd->mb_to_bottom_edge = ((cm->mi_rows - bh - mirow) * MI_SIZE) * 8;
    xd->mb_to_left_edge = -((micol * MI_SIZE) * 8);
    xd->mb_to_right_edge = ((cm->mi_cols - bw - micol) * MI_SIZE) * 8;
    xd->tile.mi_col_start = tile->mi_col_start;
    xd->tile.mi_col_end = tile->mi_col_end;
    xd->tile.mi_row_start = tile->mi_row_start;
    xd->tile.mi_row_end = tile->mi_row_end;
    xd->n8_h = bh;
    xd->n8_w = bw;
    xd->is_sec_rect = 0;
    if (xd->n8_w < xd->n8_h) {
        // Only mark is_sec_rect as 1 for the last block.
        // For PARTITION_VERT_4, it would be (0, 0, 0, 1);
        // For other partitions, it would be (0, 1).
        if (!((micol + xd->n8_w) & (xd->n8_h - 1))) xd->is_sec_rect = 1;
    }

    if (xd->n8_w > xd->n8_h)
        if (mirow & (xd->n8_w - 1)) xd->is_sec_rect = 1;

    // Adjust prediction pointers
    uint16_t *dst;
    int32_t dst_stride;
    if (plane == 0) {
        dst = (uint16_t*)(recon_buffer->buffer_y) + pred_buf_x_offest + recon_buffer->origin_x + (pred_buf_y_offest + recon_buffer->origin_y)*recon_buffer->stride_y;
        dst_stride = recon_buffer->stride_y;
    }
    else if (plane == 1) {
        dst = (uint16_t*)(recon_buffer->buffer_cb) + (pred_buf_x_offest + recon_buffer->origin_x / 2 + (pred_buf_y_offest + recon_buffer->origin_y / 2)*recon_buffer->stride_cb);
        dst_stride = recon_buffer->stride_cb;
    }
    else {
        dst = (uint16_t*)(recon_buffer->buffer_cr) + (pred_buf_x_offest + recon_buffer->origin_x / 2 + (pred_buf_y_offest + recon_buffer->origin_y / 2)*recon_buffer->stride_cr);
        dst_stride = recon_buffer->stride_cr;
    }

    int32_t chroma_up_available = xd->up_available;
    int32_t chroma_left_available = xd->left_available;

    const int32_t ss_x = plane == 0 ? 0 : 1;
    const int32_t ss_y = plane == 0 ? 0 : 1;

    if (ss_x && bw < mi_size_wide[BLOCK_8X8])
        chroma_left_available = (micol - 1) > tile->mi_col_start;
    if (ss_y && bh < mi_size_high[BLOCK_8X8])
        chroma_up_available = (mirow - 1) > tile->mi_row_start;

    int mi_stride = cm->mi_stride;
    const int32_t offset = mirow * mi_stride + micol;
    xd->mi = cm->pcs_ptr->mi_grid_base + offset;
    ModeInfo *miPtr = *xd->mi;

    if (xd->up_available) {
        // xd->above_mbmi = xd->mi[-xd->mi_stride].mbmi;
        xd->above_mbmi = &miPtr[-mi_stride].mbmi;
    }
    else
        xd->above_mbmi = NULL;
    if (xd->left_available) {
        //xd->left_mbmi = xd->mi[-1].mbmi;
        xd->left_mbmi = &miPtr[-1].mbmi;
    }
    else
        xd->left_mbmi = NULL;
    const int chroma_ref = ((mirow & 0x01) || !(bh & 0x01) || !ss_y) &&
        ((micol & 0x01) || !(bw & 0x01) || !ss_x);
    if (chroma_ref) {
        // To help calculate the "above" and "left" chroma blocks, note that the
        // current block may cover multiple luma blocks (eg, if partitioned into
        // 4x4 luma blocks).
        // First, find the top-left-most luma block covered by this chroma block

        ModeInfo *miPtr = xd->mi[-(mirow & ss_y) * mi_stride - (micol & ss_x)];

        // Then, we consider the luma region covered by the left or above 4x4 chroma
        // prediction. We want to point to the chroma reference block in that
        // region, which is the bottom-right-most mi unit.
        // This leads to the following offsets:
        MbModeInfo *chroma_above_mi =
            chroma_up_available ? &miPtr[-mi_stride + ss_x].mbmi : NULL;
        xd->chroma_above_mbmi = chroma_above_mi;

        MbModeInfo *chroma_left_mi =
            chroma_left_available ? &miPtr[ss_y * mi_stride - 1].mbmi : NULL;
        xd->chroma_left_mbmi = chroma_left_mi;
    }

    //CHKN  const MbModeInfo *const mbmi = xd->mi[0];
    const int32_t txwpx = tx_size_wide[tx_size];
    const int32_t txhpx = tx_size_high[tx_size];
    const int32_t x = col_off << tx_size_wide_log2[0];
    const int32_t y = row_off << tx_size_high_log2[0];
#if PAL_SUP
    if (use_palette) {
        int32_t r, c;
        const uint8_t *const map = palette_info->color_idx_map;
        const uint16_t *const palette =
            palette_info->pmi.palette_colors + plane * PALETTE_MAX_SIZE;
        for (r = 0; r < txhpx; ++r) {
            for (c = 0; c < txwpx; ++c) {
                dst[r * dst_stride + c] = palette[map[(r + y) * wpx + c + x]];
            }
        }
        return;
    }
#else
    //if (use_palette) {
    //  int32_t r, c;
    //  const uint8_t *const map = xd->plane[plane != 0].color_index_map;
    //  const uint16_t *const palette =
    //      mbmi->palette_mode_info.palette_colors + plane * PALETTE_MAX_SIZE;
    //  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    //    uint16_t *dst16 = CONVERT_TO_SHORTPTR(dst);
    //    for (r = 0; r < txhpx; ++r) {
    //      for (c = 0; c < txwpx; ++c) {
    //        dst16[r * dst_stride + c] = palette[map[(r + y) * wpx + c + x]];
    //      }
    //    }
    //  } else {
    //    for (r = 0; r < txhpx; ++r) {
    //      for (c = 0; c < txwpx; ++c) {
    //        dst[r * dst_stride + c] =
    //            (uint8_t)palette[map[(r + y) * wpx + c + x]];
    //      }
    //    }
    //  }
    //  return;
    //}
#endif
    //CHKN BlockSize bsize = mbmi->sb_type;

    struct MacroblockdPlane  pd_s;
    struct MacroblockdPlane * pd = &pd_s;
    if (plane == 0)
        pd->subsampling_x = pd->subsampling_y = 0;
    else
        pd->subsampling_x = pd->subsampling_y = 1;
    const int32_t txw = tx_size_wide_unit[tx_size];
    const int32_t txh = tx_size_high_unit[tx_size];
    const int32_t have_top = row_off || (pd->subsampling_y ? /*xd->*/chroma_up_available
        : xd->up_available);
    const int32_t have_left =
        col_off ||
        (pd->subsampling_x ? /*xd->*/chroma_left_available : xd->left_available);
    const int32_t mi_row = -xd->mb_to_top_edge >> (3 + MI_SIZE_LOG2);
    const int32_t mi_col = -xd->mb_to_left_edge >> (3 + MI_SIZE_LOG2);
    const int32_t xr_chr_offset = 0;
    const int32_t yd_chr_offset = 0;

    // Distance between the right edge of this prediction block to
    // the frame right edge
    const int32_t xr = (xd->mb_to_right_edge >> (3 + pd->subsampling_x)) +
        (wpx - x - txwpx) - xr_chr_offset;
    // Distance between the bottom edge of this prediction block to
    // the frame bottom edge
    const int32_t yd = (xd->mb_to_bottom_edge >> (3 + pd->subsampling_y)) +
        (hpx - y - txhpx) - yd_chr_offset;
    const int32_t right_available =
        mi_col + ((col_off + txw) << pd->subsampling_x) < xd->tile.mi_col_end;
    const int32_t bottom_available =
        (yd > 0) &&
        (mi_row + ((row_off + txh) << pd->subsampling_y) < xd->tile.mi_row_end);

    const PartitionType partition = from_shape_to_part[blk_geom->shape]; //cu_ptr->part;// PARTITION_NONE;//CHKN this is good enough as the avail functions need to know if VERT part is used or not mbmi->partition;

    // force 4x4 chroma component block size.
    bsize = scale_chroma_bsize(bsize, pd->subsampling_x, pd->subsampling_y);

    const int32_t have_top_right = intra_has_top_right(
        cm->p_pcs_ptr->sequence_control_set_ptr->seq_header.sb_size, bsize,
        mi_row, mi_col, have_top, right_available, partition, tx_size,
        row_off, col_off, pd->subsampling_x, pd->subsampling_y);
    const int32_t have_bottom_left = intra_has_bottom_left(
        cm->p_pcs_ptr->sequence_control_set_ptr->seq_header.sb_size, bsize,
        mi_row, mi_col, bottom_available, have_left, partition,
        tx_size, row_off, col_off, pd->subsampling_x, pd->subsampling_y);

    const int32_t disable_edge_filter = 0;//CHKN !cm->seq_params.enable_intra_edge_filter;

    build_intra_predictors_high(
        xd,
        topNeighArray,
        leftNeighArray,
        // ref, ref_stride,
        dst, dst_stride, mode,
        angle_delta, filter_intra_mode, tx_size,
        disable_edge_filter,
        have_top ? AOMMIN(txwpx, xr + txwpx) : 0,
        have_top_right ? AOMMIN(txwpx, xr) : 0,
        have_left ? AOMMIN(txhpx, yd + txhpx) : 0,
        have_bottom_left ? AOMMIN(txhpx, yd) : 0, plane, EB_10BIT);
}

/** IntraPrediction()
is the main function to compute intra prediction for a PU
*/
EbErrorType eb_av1_intra_prediction_cl(
    ModeDecisionContext                  *md_context_ptr,
    PictureControlSet                    *picture_control_set_ptr,
    ModeDecisionCandidateBuffer           *candidate_buffer_ptr)
{
    EbErrorType return_error = EB_ErrorNone;

    uint32_t modeTypeLeftNeighborIndex = get_neighbor_array_unit_left_index(
        md_context_ptr->mode_type_neighbor_array,
        md_context_ptr->cu_origin_y);
    uint32_t modeTypeTopNeighborIndex = get_neighbor_array_unit_top_index(
        md_context_ptr->mode_type_neighbor_array,
        md_context_ptr->cu_origin_x);
    uint32_t intraLumaModeLeftNeighborIndex = get_neighbor_array_unit_left_index(
        md_context_ptr->intra_luma_mode_neighbor_array,
        md_context_ptr->cu_origin_y);
    uint32_t intraLumaModeTopNeighborIndex = get_neighbor_array_unit_top_index(
        md_context_ptr->intra_luma_mode_neighbor_array,
        md_context_ptr->cu_origin_x);

    uint32_t intraChromaModeLeftNeighborIndex = get_neighbor_array_unit_left_index(
        md_context_ptr->intra_chroma_mode_neighbor_array,
        md_context_ptr->round_origin_y >> 1);
    uint32_t intraChromaModeTopNeighborIndex = get_neighbor_array_unit_top_index(
        md_context_ptr->intra_chroma_mode_neighbor_array,
        md_context_ptr->round_origin_x >> 1);

    md_context_ptr->intra_luma_left_mode = (uint32_t)(
        (md_context_ptr->mode_type_neighbor_array->left_array[modeTypeLeftNeighborIndex] != INTRA_MODE) ? DC_PRED/*EB_INTRA_DC*/ :
        (uint32_t)md_context_ptr->intra_luma_mode_neighbor_array->left_array[intraLumaModeLeftNeighborIndex]);

    md_context_ptr->intra_luma_top_mode = (uint32_t)(
        (md_context_ptr->mode_type_neighbor_array->top_array[modeTypeTopNeighborIndex] != INTRA_MODE) ? DC_PRED/*EB_INTRA_DC*/ :
        (uint32_t)md_context_ptr->intra_luma_mode_neighbor_array->top_array[intraLumaModeTopNeighborIndex]);       //   use DC. This seems like we could use a LCU-width

    md_context_ptr->intra_chroma_left_mode = md_context_ptr->intra_luma_left_mode;
    md_context_ptr->intra_chroma_top_mode = md_context_ptr->intra_luma_top_mode;

    md_context_ptr->intra_chroma_left_mode = (uint32_t)(
        (md_context_ptr->mode_type_neighbor_array->left_array[modeTypeLeftNeighborIndex] != INTRA_MODE) ? UV_DC_PRED :
        (uint32_t)md_context_ptr->intra_chroma_mode_neighbor_array->left_array[intraChromaModeLeftNeighborIndex]);

    md_context_ptr->intra_chroma_top_mode = (uint32_t)(
        (md_context_ptr->mode_type_neighbor_array->top_array[modeTypeTopNeighborIndex] != INTRA_MODE) ? UV_DC_PRED :
        (uint32_t)md_context_ptr->intra_chroma_mode_neighbor_array->top_array[intraChromaModeTopNeighborIndex]);       //   use DC. This seems like we could use a LCU-width
    TxSize  tx_size = md_context_ptr->blk_geom->txsize[candidate_buffer_ptr->candidate_ptr->tx_depth][0]; // Nader - Intra 128x128 not supported
    TxSize  tx_size_Chroma = md_context_ptr->blk_geom->txsize_uv[candidate_buffer_ptr->candidate_ptr->tx_depth][0]; //Nader - Intra 128x128 not supported

    if(!md_context_ptr->hbd_mode_decision) {
        uint8_t    topNeighArray[64 * 2 + 1];
        uint8_t    leftNeighArray[64 * 2 + 1];
        PredictionMode mode;
        // Hsan: plane should be derived @ an earlier stage (e.g. @ the call of perform_fast_loop())
        int32_t start_plane = (md_context_ptr->uv_search_path) ? 1 : 0;
        int32_t end_plane = (md_context_ptr->blk_geom->has_uv && md_context_ptr->chroma_level <= CHROMA_MODE_1) ? (int)MAX_MB_PLANE : 1;

        for (int32_t plane = start_plane; plane < end_plane; ++plane) {
            if (plane == 0) {
                if (md_context_ptr->cu_origin_y != 0)
                    memcpy(topNeighArray + 1, md_context_ptr->luma_recon_neighbor_array->top_array + md_context_ptr->cu_origin_x, md_context_ptr->blk_geom->bwidth * 2);
                if (md_context_ptr->cu_origin_x != 0)
                    memcpy(leftNeighArray + 1, md_context_ptr->luma_recon_neighbor_array->left_array + md_context_ptr->cu_origin_y, md_context_ptr->blk_geom->bheight * 2);
                if (md_context_ptr->cu_origin_y != 0 && md_context_ptr->cu_origin_x != 0)
                    topNeighArray[0] = leftNeighArray[0] = md_context_ptr->luma_recon_neighbor_array->top_left_array[MAX_PICTURE_HEIGHT_SIZE + md_context_ptr->cu_origin_x - md_context_ptr->cu_origin_y];
            }

            else if (plane == 1) {
                if (md_context_ptr->round_origin_y != 0)
                    memcpy(topNeighArray + 1, md_context_ptr->cb_recon_neighbor_array->top_array + md_context_ptr->round_origin_x / 2, md_context_ptr->blk_geom->bwidth_uv * 2);

                if (md_context_ptr->round_origin_x != 0)
                    memcpy(leftNeighArray + 1, md_context_ptr->cb_recon_neighbor_array->left_array + md_context_ptr->round_origin_y / 2, md_context_ptr->blk_geom->bheight_uv * 2);

                if (md_context_ptr->round_origin_y != 0 && md_context_ptr->round_origin_x != 0)
                    topNeighArray[0] = leftNeighArray[0] = md_context_ptr->cb_recon_neighbor_array->top_left_array[MAX_PICTURE_HEIGHT_SIZE / 2 + md_context_ptr->round_origin_x / 2 - md_context_ptr->round_origin_y / 2];
            }
            else {
                if (md_context_ptr->round_origin_y != 0)
                    memcpy(topNeighArray + 1, md_context_ptr->cr_recon_neighbor_array->top_array + md_context_ptr->round_origin_x / 2, md_context_ptr->blk_geom->bwidth_uv * 2);

                if (md_context_ptr->round_origin_x != 0)
                    memcpy(leftNeighArray + 1, md_context_ptr->cr_recon_neighbor_array->left_array + md_context_ptr->round_origin_y / 2, md_context_ptr->blk_geom->bheight_uv * 2);

                if (md_context_ptr->round_origin_y != 0 && md_context_ptr->round_origin_x != 0)
                    topNeighArray[0] = leftNeighArray[0] = md_context_ptr->cr_recon_neighbor_array->top_left_array[MAX_PICTURE_HEIGHT_SIZE / 2 + md_context_ptr->round_origin_x / 2 - md_context_ptr->round_origin_y / 2];
            }

            if (plane)
                mode = (candidate_buffer_ptr->candidate_ptr->intra_chroma_mode == UV_CFL_PRED) ? (PredictionMode) UV_DC_PRED : (PredictionMode) candidate_buffer_ptr->candidate_ptr->intra_chroma_mode;
            else
                mode = candidate_buffer_ptr->candidate_ptr->pred_mode;

            eb_av1_predict_intra_block(
                &md_context_ptr->sb_ptr->tile_info,
                !ED_STAGE,
                md_context_ptr->blk_geom,
                picture_control_set_ptr->parent_pcs_ptr->av1_cm,                                      //const Av1Common *cm,
                plane ? md_context_ptr->blk_geom->bwidth_uv : md_context_ptr->blk_geom->bwidth,          //int32_t wpx,
                plane ? md_context_ptr->blk_geom->bheight_uv : md_context_ptr->blk_geom->bheight,          //int32_t hpx,
                plane ? tx_size_Chroma : tx_size,                                               //TxSize tx_size,
                mode,                                                                           //PredictionMode mode,
                plane ? candidate_buffer_ptr->candidate_ptr->angle_delta[PLANE_TYPE_UV] : candidate_buffer_ptr->candidate_ptr->angle_delta[PLANE_TYPE_Y],
#if PAL_SUP
                plane==0 ? (candidate_buffer_ptr->candidate_ptr->palette_info.pmi.palette_size[0]>0) : 0,
                plane==0 ? &candidate_buffer_ptr->candidate_ptr->palette_info : NULL,    //MD
#else
                0,                                                                              //int32_t use_palette,
#endif
#if FILTER_INTRA_FLAG
                plane ? FILTER_INTRA_MODES : candidate_buffer_ptr->candidate_ptr->filter_intra_mode,
#else
                FILTER_INTRA_MODES,                                                             //CHKN FilterIntraMode filter_intra_mode,
#endif
                topNeighArray + 1,
                leftNeighArray + 1,
                candidate_buffer_ptr->prediction_ptr,                                              //uint8_t *dst,
                                                                                                //int32_t dst_stride,
                0,                                                                              //int32_t col_off,
                0,                                                                              //int32_t row_off,
                plane,                                                                          //int32_t plane,
                md_context_ptr->blk_geom->bsize,       //uint32_t puSize,
                md_context_ptr->cu_origin_x,
                md_context_ptr->cu_origin_y,
                md_context_ptr->cu_origin_x,                  //uint32_t cuOrgX,
                md_context_ptr->cu_origin_y,                  //uint32_t cuOrgY
                plane ? ((md_context_ptr->blk_geom->origin_x >> 3) << 3) / 2 : md_context_ptr->blk_geom->origin_x,  //uint32_t cuOrgX used only for prediction Ptr
                plane ? ((md_context_ptr->blk_geom->origin_y >> 3) << 3) / 2 : md_context_ptr->blk_geom->origin_y   //uint32_t cuOrgY used only for prediction Ptr
            );
        }
    } else {
        uint16_t    topNeighArray[64 * 2 + 1];
        uint16_t    leftNeighArray[64 * 2 + 1];
        PredictionMode mode;
        // Hsan: plane should be derived @ an earlier stage (e.g. @ the call of perform_fast_loop())
        int32_t start_plane = (md_context_ptr->uv_search_path) ? 1 : 0;
        int32_t end_plane = (md_context_ptr->blk_geom->has_uv && md_context_ptr->chroma_level <= CHROMA_MODE_1) ? (int)MAX_MB_PLANE : 1;
        for (int32_t plane = start_plane; plane < end_plane; ++plane) {
            if (plane == 0) {
                if (md_context_ptr->cu_origin_y != 0)
                    memcpy(topNeighArray + 1, (uint16_t*)(md_context_ptr->luma_recon_neighbor_array16bit->top_array) + md_context_ptr->cu_origin_x, md_context_ptr->blk_geom->bwidth * 2 * sizeof(uint16_t));

                if (md_context_ptr->cu_origin_x != 0)
                    memcpy(leftNeighArray + 1, (uint16_t*)(md_context_ptr->luma_recon_neighbor_array16bit->left_array) + md_context_ptr->cu_origin_y, md_context_ptr->blk_geom->bheight * 2 * sizeof(uint16_t));

                if (md_context_ptr->cu_origin_y != 0 && md_context_ptr->cu_origin_x != 0)
                    topNeighArray[0] = leftNeighArray[0] = ((uint16_t*)(md_context_ptr->luma_recon_neighbor_array16bit->top_left_array) + MAX_PICTURE_HEIGHT_SIZE + md_context_ptr->cu_origin_x - md_context_ptr->cu_origin_y)[0];
            }
            else if (plane == 1) {
                if (md_context_ptr->round_origin_y != 0)
                    memcpy(topNeighArray + 1, (uint16_t*)(md_context_ptr->cb_recon_neighbor_array16bit->top_array) + md_context_ptr->round_origin_x / 2, md_context_ptr->blk_geom->bwidth_uv * 2 * sizeof(uint16_t));

                if (md_context_ptr->round_origin_x != 0)
                    memcpy(leftNeighArray + 1, (uint16_t*)(md_context_ptr->cb_recon_neighbor_array16bit->left_array) + md_context_ptr->round_origin_y / 2, md_context_ptr->blk_geom->bheight_uv * 2 * sizeof(uint16_t));

                if (md_context_ptr->round_origin_y != 0 && md_context_ptr->round_origin_x != 0)
                    topNeighArray[0] = leftNeighArray[0] = ((uint16_t*) (md_context_ptr->cb_recon_neighbor_array16bit->top_left_array) + MAX_PICTURE_HEIGHT_SIZE / 2 + md_context_ptr->round_origin_x / 2 - md_context_ptr->round_origin_y / 2)[0];
            }
            else {
                if (md_context_ptr->round_origin_y != 0)
                    memcpy(topNeighArray + 1, (uint16_t*)(md_context_ptr->cr_recon_neighbor_array16bit->top_array) + md_context_ptr->round_origin_x / 2, md_context_ptr->blk_geom->bwidth_uv * 2 * sizeof(uint16_t));

                if (md_context_ptr->round_origin_x != 0)
                    memcpy(leftNeighArray + 1, (uint16_t*)(md_context_ptr->cr_recon_neighbor_array16bit->left_array) + md_context_ptr->round_origin_y / 2, md_context_ptr->blk_geom->bheight_uv * 2 * sizeof(uint16_t));

                if (md_context_ptr->round_origin_y != 0 && md_context_ptr->round_origin_x != 0)
                    topNeighArray[0] = leftNeighArray[0] = ((uint16_t*) (md_context_ptr->cr_recon_neighbor_array16bit->top_left_array) + MAX_PICTURE_HEIGHT_SIZE / 2 + md_context_ptr->round_origin_x / 2 - md_context_ptr->round_origin_y / 2)[0];
            }

            if (plane)
                mode = (candidate_buffer_ptr->candidate_ptr->intra_chroma_mode == UV_CFL_PRED) ? (PredictionMode) UV_DC_PRED : (PredictionMode) candidate_buffer_ptr->candidate_ptr->intra_chroma_mode;
            else
                mode = candidate_buffer_ptr->candidate_ptr->pred_mode;

            eb_av1_predict_intra_block_16bit(
                &md_context_ptr->sb_ptr->tile_info,
                !ED_STAGE,
                md_context_ptr->blk_geom,
                picture_control_set_ptr->parent_pcs_ptr->av1_cm,                                      //const Av1Common *cm,
                plane ? md_context_ptr->blk_geom->bwidth_uv : md_context_ptr->blk_geom->bwidth,          //int32_t wpx,
                plane ? md_context_ptr->blk_geom->bheight_uv : md_context_ptr->blk_geom->bheight,          //int32_t hpx,
                plane ? tx_size_Chroma : tx_size,                                               //TxSize tx_size,
                mode,                                                                           //PredictionMode mode,
                plane ? candidate_buffer_ptr->candidate_ptr->angle_delta[PLANE_TYPE_UV] : candidate_buffer_ptr->candidate_ptr->angle_delta[PLANE_TYPE_Y],
#if PAL_SUP
                plane == 0 ? (candidate_buffer_ptr->candidate_ptr->palette_info.pmi.palette_size[0] > 0) : 0,
                plane == 0 ? &candidate_buffer_ptr->candidate_ptr->palette_info : NULL,    //MD
#else
                0,                                                                              //int32_t use_palette,
#endif
#if FILTER_INTRA_FLAG
                plane ? FILTER_INTRA_MODES : candidate_buffer_ptr->candidate_ptr->filter_intra_mode,
#else
                FILTER_INTRA_MODES,                                                             //CHKN FilterIntraMode filter_intra_mode,
#endif
                topNeighArray + 1,
                leftNeighArray + 1,
                candidate_buffer_ptr->prediction_ptr,                                              //uint8_t *dst,
                0,                                                                              //int32_t col_off,
                0,                                                                              //int32_t row_off,
                plane,                                                                          //int32_t plane,
                md_context_ptr->blk_geom->bsize,       //uint32_t puSize,
                md_context_ptr->cu_origin_x,
                md_context_ptr->cu_origin_y,
                md_context_ptr->cu_origin_x,                  //uint32_t cuOrgX,
                md_context_ptr->cu_origin_y,                  //uint32_t cuOrgY
                plane ? ((md_context_ptr->blk_geom->origin_x >> 3) << 3) / 2 : md_context_ptr->blk_geom->origin_x,  //uint32_t cuOrgX used only for prediction Ptr
                plane ? ((md_context_ptr->blk_geom->origin_y >> 3) << 3) / 2 : md_context_ptr->blk_geom->origin_y   //uint32_t cuOrgY used only for prediction Ptr
            );
        }
    }

    return return_error;
}

#if II_COMP_FLAG
EbErrorType  intra_luma_prediction_for_interintra(
    ModeDecisionContext         *md_context_ptr,
    PictureControlSet           *picture_control_set_ptr,
    INTERINTRA_MODE              interintra_mode,
    EbPictureBufferDesc         *prediction_ptr)
{
    EbErrorType return_error = EB_ErrorNone;

    uint32_t mode_type_left_neighbor_index = get_neighbor_array_unit_left_index(
        md_context_ptr->mode_type_neighbor_array,
        md_context_ptr->cu_origin_y);
    uint32_t mode_type_top_neighbor_index = get_neighbor_array_unit_top_index(
        md_context_ptr->mode_type_neighbor_array,
        md_context_ptr->cu_origin_x);
    uint32_t intra_luma_mode_left_neighbor_index = get_neighbor_array_unit_left_index(
        md_context_ptr->intra_luma_mode_neighbor_array,
        md_context_ptr->cu_origin_y);
    uint32_t intra_luma_mode_top_neighbor_index = get_neighbor_array_unit_top_index(
        md_context_ptr->intra_luma_mode_neighbor_array,
        md_context_ptr->cu_origin_x);

    md_context_ptr->intra_luma_left_mode = (uint32_t)(
        (md_context_ptr->mode_type_neighbor_array->left_array[mode_type_left_neighbor_index] != INTRA_MODE) ? DC_PRED:
        (uint32_t)md_context_ptr->intra_luma_mode_neighbor_array->left_array[intra_luma_mode_left_neighbor_index]);

    md_context_ptr->intra_luma_top_mode = (uint32_t)(
        (md_context_ptr->mode_type_neighbor_array->top_array[mode_type_top_neighbor_index] != INTRA_MODE) ? DC_PRED:
        (uint32_t)md_context_ptr->intra_luma_mode_neighbor_array->top_array[intra_luma_mode_top_neighbor_index]);       //   use DC. This seems like we could use a LCU-width

    TxSize  tx_size = md_context_ptr->blk_geom->txsize[0][0];  //CHKN  TOcheck
    PredictionMode mode = interintra_to_intra_mode[interintra_mode];

    if (!md_context_ptr->hbd_mode_decision) {
        uint8_t    top_neigh_array[64 * 2 + 1];
        uint8_t    left_neigh_array[64 * 2 + 1];

        if (md_context_ptr->cu_origin_y != 0)
            memcpy(top_neigh_array + 1, md_context_ptr->luma_recon_neighbor_array->top_array + md_context_ptr->cu_origin_x, md_context_ptr->blk_geom->bwidth * 2);
        if (md_context_ptr->cu_origin_x != 0)
            memcpy(left_neigh_array + 1, md_context_ptr->luma_recon_neighbor_array->left_array + md_context_ptr->cu_origin_y, md_context_ptr->blk_geom->bheight * 2);
        if (md_context_ptr->cu_origin_y != 0 && md_context_ptr->cu_origin_x != 0)
            top_neigh_array[0] = left_neigh_array[0] = md_context_ptr->luma_recon_neighbor_array->top_left_array[MAX_PICTURE_HEIGHT_SIZE + md_context_ptr->cu_origin_x - md_context_ptr->cu_origin_y];

        eb_av1_predict_intra_block(
            &md_context_ptr->sb_ptr->tile_info,
            !ED_STAGE,
            md_context_ptr->blk_geom,
            picture_control_set_ptr->parent_pcs_ptr->av1_cm,        //const Av1Common *cm,
            md_context_ptr->blk_geom->bwidth,                       //int32_t wpx,
            md_context_ptr->blk_geom->bheight,                      //int32_t hpx,
            tx_size,                                                //TxSize tx_size,
            mode,                                                   //PredictionMode mode,
            0,                                                      //candidate_buffer_ptr->candidate_ptr->angle_delta[PLANE_TYPE_Y],
            0,                                                      //int32_t use_palette,
#if PAL_SUP
            NULL,  //Inter-Intra
#endif
            FILTER_INTRA_MODES,                                     //CHKN FilterIntraMode filter_intra_mode,
            top_neigh_array + 1,
            left_neigh_array + 1,
            prediction_ptr,                                         //uint8_t *dst,
            md_context_ptr->blk_geom->tx_boff_x[0][0] >> 2,         //int32_t col_off,
            md_context_ptr->blk_geom->tx_boff_y[0][0] >> 2,         //int32_t row_off,
            PLANE_TYPE_Y,                                           //int32_t plane,
            md_context_ptr->blk_geom->bsize,                        //uint32_t puSize,
            md_context_ptr->cu_origin_x,
            md_context_ptr->cu_origin_y,
            md_context_ptr->cu_origin_x,                            //uint32_t cuOrgX,
            md_context_ptr->cu_origin_y,                            //uint32_t cuOrgY
            0,                                                      //cuOrgX used only for prediction Ptr
            0                                                       //cuOrgY used only for prediction Ptr
        );
    } else {
        uint16_t top_neigh_array[64 * 2 + 1];
        uint16_t left_neigh_array[64 * 2 + 1];

        if (md_context_ptr->cu_origin_y != 0)
            memcpy(top_neigh_array + 1, (uint16_t*)(md_context_ptr->luma_recon_neighbor_array16bit->top_array) + md_context_ptr->cu_origin_x, md_context_ptr->blk_geom->bwidth * 2 * sizeof(uint16_t));
        if (md_context_ptr->cu_origin_x != 0)
            memcpy(left_neigh_array + 1, (uint16_t*)(md_context_ptr->luma_recon_neighbor_array16bit->left_array) + md_context_ptr->cu_origin_y, md_context_ptr->blk_geom->bheight * 2 * sizeof(uint16_t));
        if (md_context_ptr->cu_origin_y != 0 && md_context_ptr->cu_origin_x != 0)
            top_neigh_array[0] = left_neigh_array[0] = ((uint16_t*)(md_context_ptr->luma_recon_neighbor_array16bit->top_left_array) + MAX_PICTURE_HEIGHT_SIZE + md_context_ptr->cu_origin_x - md_context_ptr->cu_origin_y)[0];

        eb_av1_predict_intra_block_16bit(
            &md_context_ptr->sb_ptr->tile_info,
            !ED_STAGE,
            md_context_ptr->blk_geom,
            picture_control_set_ptr->parent_pcs_ptr->av1_cm,        //const Av1Common *cm,
            md_context_ptr->blk_geom->bwidth,                       //int32_t wpx,
            md_context_ptr->blk_geom->bheight,                      //int32_t hpx,
            tx_size,                                                //TxSize tx_size,
            mode,                                                   //PredictionMode mode,
            0,                                                      //candidate_buffer_ptr->candidate_ptr->angle_delta[PLANE_TYPE_Y],
            0,                                                      //int32_t use_palette,
#if PAL_SUP
            NULL,  //Inter-Intra
#endif
            FILTER_INTRA_MODES,                                     //CHKN FilterIntraMode filter_intra_mode,
            top_neigh_array + 1,
            left_neigh_array + 1,
            prediction_ptr,                                         //uint8_t *dst,
            md_context_ptr->blk_geom->tx_boff_x[0][0] >> 2,         //int32_t col_off,
            md_context_ptr->blk_geom->tx_boff_y[0][0] >> 2,         //int32_t row_off,
            PLANE_TYPE_Y,                                           //int32_t plane,
            md_context_ptr->blk_geom->bsize,                        //uint32_t puSize,
            md_context_ptr->cu_origin_x,
            md_context_ptr->cu_origin_y,
            md_context_ptr->cu_origin_x,                            //uint32_t cuOrgX,
            md_context_ptr->cu_origin_y,                            //uint32_t cuOrgY
            0,                                                      //cuOrgX used only for prediction Ptr
            0                                                       //cuOrgY used only for prediction Ptr
        );
    }

    return return_error;
}
#endif

EbErrorType update_neighbor_samples_array_open_loop(
        uint8_t                           *above_ref,
        uint8_t                            *left_ref,
        EbPictureBufferDesc              *input_ptr,
        uint32_t                            stride,
        uint32_t                            src_origin_x,
        uint32_t                            src_origin_y,
        uint8_t                             bwidth,
        uint8_t                             bheight)
{
    EbErrorType    return_error = EB_ErrorNone;

    uint32_t idx;
    uint8_t  *src_ptr;
    uint8_t  *read_ptr;
    uint32_t count;

    uint32_t width = input_ptr->width;
    uint32_t height = input_ptr->height;
    uint32_t block_size_half = bwidth << 1;

    // Adjust the Source ptr to start at the origin of the block being updated
    src_ptr = input_ptr->buffer_y + (((src_origin_y + input_ptr->origin_y) * stride) + (src_origin_x + input_ptr->origin_x));

    //Initialise the Luma Intra Reference Array to the mid range value 128 (for CUs at the picture boundaries)
    EB_MEMSET(above_ref, 127, (bwidth << 1) + 1);
    EB_MEMSET(left_ref, 129, (bheight << 1) + 1);

    // Get the upper left sample
    if (src_origin_x != 0 && src_origin_y != 0) {
        read_ptr = src_ptr - stride - 1;
        *above_ref = *read_ptr;
        *left_ref = *read_ptr;
        left_ref++;
        above_ref++;
    }else {
        *above_ref = *left_ref = 128;
        left_ref++;
        above_ref++;
    }
    // Get the left-column
    count = block_size_half;
    if (src_origin_x != 0) {
        read_ptr = src_ptr - 1;
        count = ((src_origin_y + count) > height) ? count - ((src_origin_y + count) - height) : count;
        for (idx = 0; idx < count; ++idx) {
            *left_ref = *read_ptr;
            read_ptr += stride;
            left_ref++;
        }
        left_ref += (block_size_half - count);
    }else
        left_ref += count;

    // Get the top-row
    count = block_size_half;
    if (src_origin_y != 0) {
        read_ptr = src_ptr - stride;
        count = ((src_origin_x + count) > width) ? count - ((src_origin_x + count) - width) : count;
        EB_MEMCPY(above_ref, read_ptr, count);
        above_ref += (block_size_half - count);
    }else
        above_ref += count;

    return return_error;
}
/** intra_prediction_open_loop()
        performs Open-loop Intra candidate Search for a CU
 */
EbErrorType intra_prediction_open_loop(
     int32_t  p_angle ,
        uint8_t                          ois_intra_mode,
        uint32_t                         src_origin_x,
        uint32_t                         src_origin_y,
        TxSize                          tx_size,
        uint8_t                         *above_row,
        uint8_t                         *left_col,
        MotionEstimationContext_t       *context_ptr)                  // input parameter, ME context

{
    EbErrorType                return_error = EB_ErrorNone;
    PredictionMode mode = ois_intra_mode;
    const int32_t is_dr_mode = av1_is_directional_mode(mode);
    uint8_t *dst = (&(context_ptr->me_context_ptr->sb_buffer[0]));
    uint32_t dst_stride = context_ptr->me_context_ptr->sb_buffer_stride;

    if (is_dr_mode)
        dr_predictor(dst, dst_stride, tx_size, above_row, left_col, 0, 0, p_angle);
    else {
        // predict
        if (mode == DC_PRED) {
            dc_pred[src_origin_x > 0][src_origin_y > 0][tx_size](dst, dst_stride, above_row, left_col);
        } else
            pred[mode][tx_size](dst, dst_stride, above_row, left_col);
    }
    return return_error;
}
