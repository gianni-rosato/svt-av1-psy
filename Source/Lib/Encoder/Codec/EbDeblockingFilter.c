// clang-format off
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

#include <string.h>

#include "EbDeblockingFilter.h"
#include "EbDefinitions.h"
#include "EbUtility.h"
#include "EbPictureControlSet.h"
#include "EbCodingUnit.h"
#include "EbSequenceControlSet.h"
#include "EbReferenceObject.h"
#include "EbCommonUtils.h"
#include "EbLog.h"

#define   convertToChromaQp(iQpY)  ( ((iQpY) < 0) ? (iQpY) : (((iQpY) > 57) ? ((iQpY)-6) : (int32_t)(map_chroma_qp((uint32_t)iQpY))) )

static const int delta_lf_id_lut[MAX_MB_PLANE][2] = { { 0, 1 },
                                                      { 2, 2 },
                                                      { 3, 3 } };

static const SEG_LVL_FEATURES seg_lvl_lf_lut[MAX_MB_PLANE][2] = {
  { SEG_LVL_ALT_LF_Y_V, SEG_LVL_ALT_LF_Y_H },
  { SEG_LVL_ALT_LF_U, SEG_LVL_ALT_LF_U },
  { SEG_LVL_ALT_LF_V, SEG_LVL_ALT_LF_V }
};

/** setQpArrayBasedOnCU()
is used to set qp in the qp_array on a CU basis.
*/
void set_qp_array_based_on_cu(
    PictureControlSet *picture_control_set_ptr,          //input parameter
    uint32_t               cuPos_x,                       //input parameter, sample-based horizontal picture-wise locatin of the CU
    uint32_t               cuPos_y,                       //input parameter, sample-based vertical picture-wise locatin of the CU
    uint32_t               cu_size_in_min_cu_size,             //input parameter
    uint32_t               cu_qp)                          //input parameter, Qp of the CU
{
    uint32_t verticalIdx;
    uint32_t qpArrayIdx = (cuPos_y / MIN_BLOCK_SIZE) * picture_control_set_ptr->qp_array_stride + (cuPos_x / MIN_BLOCK_SIZE);

    for (verticalIdx = 0; verticalIdx < cu_size_in_min_cu_size; ++verticalIdx) {
        EB_MEMSET(picture_control_set_ptr->qp_array + qpArrayIdx + verticalIdx * picture_control_set_ptr->qp_array_stride,
            cu_qp, sizeof(uint8_t)*cu_size_in_min_cu_size);
    }

    return;
}

static INLINE int8_t signed_char_clamp(int32_t t) {
    return (int8_t)clamp(t, -128, 127);
}

static INLINE int16_t signed_char_clamp_high(int32_t t, int32_t bd) {
    switch (bd) {
    case 10: return (int16_t)clamp(t, -128 * 4, 128 * 4 - 1);
    case 12: return (int16_t)clamp(t, -128 * 16, 128 * 16 - 1);
    case 8:
    default: return (int16_t)clamp(t, -128, 128 - 1);
    }
}

// should we apply any filter at all: 11111111 yes, 00000000 no
static INLINE int8_t filter_mask2(uint8_t limit, uint8_t blimit, uint8_t p1,
    uint8_t p0, uint8_t q0, uint8_t q1) {
    int8_t mask = 0;
    mask |= (abs(p1 - p0) > limit) * -1;
    mask |= (abs(q1 - q0) > limit) * -1;
    mask |= (abs(p0 - q0) * 2 + abs(p1 - q1) / 2 > blimit) * -1;
    return ~mask;
}

static INLINE int8_t filter_mask(uint8_t limit, uint8_t blimit, uint8_t p3,
    uint8_t p2, uint8_t p1, uint8_t p0, uint8_t q0,
    uint8_t q1, uint8_t q2, uint8_t q3) {
    int8_t mask = 0;
    mask |= (abs(p3 - p2) > limit) * -1;
    mask |= (abs(p2 - p1) > limit) * -1;
    mask |= (abs(p1 - p0) > limit) * -1;
    mask |= (abs(q1 - q0) > limit) * -1;
    mask |= (abs(q2 - q1) > limit) * -1;
    mask |= (abs(q3 - q2) > limit) * -1;
    mask |= (abs(p0 - q0) * 2 + abs(p1 - q1) / 2 > blimit) * -1;
    return ~mask;
}

static INLINE int8_t filter_mask3_chroma(uint8_t limit, uint8_t blimit,
    uint8_t p2, uint8_t p1, uint8_t p0,
    uint8_t q0, uint8_t q1, uint8_t q2) {
    int8_t mask = 0;
    mask |= (abs(p2 - p1) > limit) * -1;
    mask |= (abs(p1 - p0) > limit) * -1;
    mask |= (abs(q1 - q0) > limit) * -1;
    mask |= (abs(q2 - q1) > limit) * -1;
    mask |= (abs(p0 - q0) * 2 + abs(p1 - q1) / 2 > blimit) * -1;
    return ~mask;
}

static INLINE int8_t flat_mask3_chroma(uint8_t thresh, uint8_t p2, uint8_t p1,
    uint8_t p0, uint8_t q0, uint8_t q1,
    uint8_t q2) {
    int8_t mask = 0;
    mask |= (abs(p1 - p0) > thresh) * -1;
    mask |= (abs(q1 - q0) > thresh) * -1;
    mask |= (abs(p2 - p0) > thresh) * -1;
    mask |= (abs(q2 - q0) > thresh) * -1;
    return ~mask;
}

static INLINE int8_t flat_mask4(uint8_t thresh, uint8_t p3, uint8_t p2,
    uint8_t p1, uint8_t p0, uint8_t q0, uint8_t q1,
    uint8_t q2, uint8_t q3) {
    int8_t mask = 0;
    mask |= (abs(p1 - p0) > thresh) * -1;
    mask |= (abs(q1 - q0) > thresh) * -1;
    mask |= (abs(p2 - p0) > thresh) * -1;
    mask |= (abs(q2 - q0) > thresh) * -1;
    mask |= (abs(p3 - p0) > thresh) * -1;
    mask |= (abs(q3 - q0) > thresh) * -1;
    return ~mask;
}

// is there high edge variance internal edge: 11111111 yes, 00000000 no
static INLINE int8_t hev_mask(uint8_t thresh, uint8_t p1, uint8_t p0,
    uint8_t q0, uint8_t q1) {
    int8_t hev = 0;
    hev |= (abs(p1 - p0) > thresh) * -1;
    hev |= (abs(q1 - q0) > thresh) * -1;
    return hev;
}

static INLINE void filter4(int8_t mask, uint8_t thresh, uint8_t *op1,
    uint8_t *op0, uint8_t *oq0, uint8_t *oq1) {
    int8_t filter1, filter2;

    const int8_t ps1 = (int8_t)*op1 ^ 0x80;
    const int8_t ps0 = (int8_t)*op0 ^ 0x80;
    const int8_t qs0 = (int8_t)*oq0 ^ 0x80;
    const int8_t qs1 = (int8_t)*oq1 ^ 0x80;
    const uint8_t hev = hev_mask(thresh, *op1, *op0, *oq0, *oq1);

    // add outer taps if we have high edge variance
    int8_t filter = signed_char_clamp(ps1 - qs1) & hev;

    // inner taps
    filter = signed_char_clamp(filter + 3 * (qs0 - ps0)) & mask;

    // save bottom 3 bits so that we round one side +4 and the other +3
    // if it equals 4 we'll set to adjust by -1 to account for the fact
    // we'd round 3 the other way
    filter1 = signed_char_clamp(filter + 4) >> 3;
    filter2 = signed_char_clamp(filter + 3) >> 3;

    *oq0 = signed_char_clamp(qs0 - filter1) ^ 0x80;
    *op0 = signed_char_clamp(ps0 + filter2) ^ 0x80;

    // outer tap adjustments
    filter = ROUND_POWER_OF_TWO(filter1, 1) & ~hev;

    *oq1 = signed_char_clamp(qs1 - filter) ^ 0x80;
    *op1 = signed_char_clamp(ps1 + filter) ^ 0x80;
}

void aom_lpf_horizontal_4_c(uint8_t *s, int32_t p /* pitch */,
    const uint8_t *blimit, const uint8_t *limit,
    const uint8_t *thresh) {
    int32_t i;
    int32_t count = 4;

    // loop filter designed to work using chars so that we can make maximum use
    // of 8 bit simd instructions.
    for (i = 0; i < count; ++i) {
        const uint8_t p1 = s[-2 * p], p0 = s[-p];
        const uint8_t q0 = s[0 * p], q1 = s[1 * p];
        const int8_t mask = filter_mask2(*limit, *blimit, p1, p0, q0, q1);
        filter4(mask, *thresh, s - 2 * p, s - 1 * p, s, s + 1 * p);
        ++s;
    }
}

void aom_lpf_vertical_4_c(uint8_t *s, int32_t pitch, const uint8_t *blimit,
    const uint8_t *limit, const uint8_t *thresh) {
    int32_t i;
    int32_t count = 4;

    // loop filter designed to work using chars so that we can make maximum use
    // of 8 bit simd instructions.
    for (i = 0; i < count; ++i) {
        const uint8_t p1 = s[-2], p0 = s[-1];
        const uint8_t q0 = s[0], q1 = s[1];
        const int8_t mask = filter_mask2(*limit, *blimit, p1, p0, q0, q1);
        filter4(mask, *thresh, s - 2, s - 1, s, s + 1);
        s += pitch;
    }
}


static INLINE void filter6(int8_t mask, uint8_t thresh, int8_t flat,
    uint8_t *op2, uint8_t *op1, uint8_t *op0,
    uint8_t *oq0, uint8_t *oq1, uint8_t *oq2) {
    if (flat && mask) {
        const uint8_t p2 = *op2, p1 = *op1, p0 = *op0;
        const uint8_t q0 = *oq0, q1 = *oq1, q2 = *oq2;

        // 5-tap filter [1, 2, 2, 2, 1]
        *op1 = ROUND_POWER_OF_TWO(p2 * 3 + p1 * 2 + p0 * 2 + q0, 3);
        *op0 = ROUND_POWER_OF_TWO(p2 + p1 * 2 + p0 * 2 + q0 * 2 + q1, 3);
        *oq0 = ROUND_POWER_OF_TWO(p1 + p0 * 2 + q0 * 2 + q1 * 2 + q2, 3);
        *oq1 = ROUND_POWER_OF_TWO(p0 + q0 * 2 + q1 * 2 + q2 * 3, 3);
    }
    else
        filter4(mask, thresh, op1, op0, oq0, oq1);
}

static INLINE void filter8(int8_t mask, uint8_t thresh, int8_t flat,
    uint8_t *op3, uint8_t *op2, uint8_t *op1,
    uint8_t *op0, uint8_t *oq0, uint8_t *oq1,
    uint8_t *oq2, uint8_t *oq3) {
    if (flat && mask) {
        const uint8_t p3 = *op3, p2 = *op2, p1 = *op1, p0 = *op0;
        const uint8_t q0 = *oq0, q1 = *oq1, q2 = *oq2, q3 = *oq3;

        // 7-tap filter [1, 1, 1, 2, 1, 1, 1]
        *op2 = ROUND_POWER_OF_TWO(p3 + p3 + p3 + 2 * p2 + p1 + p0 + q0, 3);
        *op1 = ROUND_POWER_OF_TWO(p3 + p3 + p2 + 2 * p1 + p0 + q0 + q1, 3);
        *op0 = ROUND_POWER_OF_TWO(p3 + p2 + p1 + 2 * p0 + q0 + q1 + q2, 3);
        *oq0 = ROUND_POWER_OF_TWO(p2 + p1 + p0 + 2 * q0 + q1 + q2 + q3, 3);
        *oq1 = ROUND_POWER_OF_TWO(p1 + p0 + q0 + 2 * q1 + q2 + q3 + q3, 3);
        *oq2 = ROUND_POWER_OF_TWO(p0 + q0 + q1 + 2 * q2 + q3 + q3 + q3, 3);
    }
    else
        filter4(mask, thresh, op1, op0, oq0, oq1);
}

void aom_lpf_horizontal_6_c(uint8_t *s, int32_t p, const uint8_t *blimit,
    const uint8_t *limit, const uint8_t *thresh) {
    int32_t i;
    int32_t count = 4;

    // loop filter designed to work using chars so that we can make maximum use
    // of 8 bit simd instructions.
    for (i = 0; i < count; ++i) {
        const uint8_t p2 = s[-3 * p], p1 = s[-2 * p], p0 = s[-p];
        const uint8_t q0 = s[0 * p], q1 = s[1 * p], q2 = s[2 * p];

        const int8_t mask =
            filter_mask3_chroma(*limit, *blimit, p2, p1, p0, q0, q1, q2);
        const int8_t flat = flat_mask3_chroma(1, p2, p1, p0, q0, q1, q2);
        filter6(mask, *thresh, flat, s - 3 * p, s - 2 * p, s - 1 * p, s, s + 1 * p,
            s + 2 * p);
        ++s;
    }
}

void aom_lpf_horizontal_8_c(uint8_t *s, int32_t p, const uint8_t *blimit,
    const uint8_t *limit, const uint8_t *thresh) {
    int32_t i;
    int32_t count = 4;

    // loop filter designed to work using chars so that we can make maximum use
    // of 8 bit simd instructions.
    for (i = 0; i < count; ++i) {
        const uint8_t p3 = s[-4 * p], p2 = s[-3 * p], p1 = s[-2 * p], p0 = s[-p];
        const uint8_t q0 = s[0 * p], q1 = s[1 * p], q2 = s[2 * p], q3 = s[3 * p];

        const int8_t mask =
            filter_mask(*limit, *blimit, p3, p2, p1, p0, q0, q1, q2, q3);
        const int8_t flat = flat_mask4(1, p3, p2, p1, p0, q0, q1, q2, q3);
        filter8(mask, *thresh, flat, s - 4 * p, s - 3 * p, s - 2 * p, s - 1 * p, s,
            s + 1 * p, s + 2 * p, s + 3 * p);
        ++s;
    }
}

void aom_lpf_vertical_8_c(uint8_t *s, int32_t pitch, const uint8_t *blimit,
    const uint8_t *limit, const uint8_t *thresh) {
    int32_t i;
    int32_t count = 4;

    for (i = 0; i < count; ++i) {
        const uint8_t p3 = s[-4], p2 = s[-3], p1 = s[-2], p0 = s[-1];
        const uint8_t q0 = s[0], q1 = s[1], q2 = s[2], q3 = s[3];
        const int8_t mask =
            filter_mask(*limit, *blimit, p3, p2, p1, p0, q0, q1, q2, q3);
        const int8_t flat = flat_mask4(1, p3, p2, p1, p0, q0, q1, q2, q3);
        filter8(mask, *thresh, flat, s - 4, s - 3, s - 2, s - 1, s, s + 1, s + 2,
            s + 3);
        s += pitch;
    }
}

// Should we apply any filter at all: 11111111 yes, 00000000 no ?
static INLINE int8_t highbd_filter_mask2(uint8_t limit, uint8_t blimit,
    uint16_t p1, uint16_t p0, uint16_t q0,
    uint16_t q1, int32_t bd) {
    int8_t mask = 0;
    int16_t limit16 = (uint16_t)limit << (bd - 8);
    int16_t blimit16 = (uint16_t)blimit << (bd - 8);
    mask |= (abs(p1 - p0) > limit16) * -1;
    mask |= (abs(q1 - q0) > limit16) * -1;
    mask |= (abs(p0 - q0) * 2 + abs(p1 - q1) / 2 > blimit16) * -1;
    return ~mask;
}

// Should we apply any filter at all: 11111111 yes, 00000000 no ?
static INLINE int8_t highbd_filter_mask(uint8_t limit, uint8_t blimit,
    uint16_t p3, uint16_t p2, uint16_t p1,
    uint16_t p0, uint16_t q0, uint16_t q1,
    uint16_t q2, uint16_t q3, int32_t bd) {
    int8_t mask = 0;
    int16_t limit16 = (uint16_t)limit << (bd - 8);
    int16_t blimit16 = (uint16_t)blimit << (bd - 8);
    mask |= (abs(p3 - p2) > limit16) * -1;
    mask |= (abs(p2 - p1) > limit16) * -1;
    mask |= (abs(p1 - p0) > limit16) * -1;
    mask |= (abs(q1 - q0) > limit16) * -1;
    mask |= (abs(q2 - q1) > limit16) * -1;
    mask |= (abs(q3 - q2) > limit16) * -1;
    mask |= (abs(p0 - q0) * 2 + abs(p1 - q1) / 2 > blimit16) * -1;
    return ~mask;
}

static INLINE int8_t highbd_flat_mask4(uint8_t thresh, uint16_t p3, uint16_t p2,
    uint16_t p1, uint16_t p0, uint16_t q0,
    uint16_t q1, uint16_t q2, uint16_t q3,
    int32_t bd) {
    int8_t mask = 0;
    int16_t thresh16 = (uint16_t)thresh << (bd - 8);
    mask |= (abs(p1 - p0) > thresh16) * -1;
    mask |= (abs(q1 - q0) > thresh16) * -1;
    mask |= (abs(p2 - p0) > thresh16) * -1;
    mask |= (abs(q2 - q0) > thresh16) * -1;
    mask |= (abs(p3 - p0) > thresh16) * -1;
    mask |= (abs(q3 - q0) > thresh16) * -1;
    return ~mask;
}

// Is there high edge variance internal edge:
// 11111111_11111111 yes, 00000000_00000000 no ?
static INLINE int16_t highbd_hev_mask(uint8_t thresh, uint16_t p1, uint16_t p0,
    uint16_t q0, uint16_t q1, int32_t bd) {
    int16_t hev = 0;
    int16_t thresh16 = (uint16_t)thresh << (bd - 8);
    hev |= (abs(p1 - p0) > thresh16) * -1;
    hev |= (abs(q1 - q0) > thresh16) * -1;
    return hev;
}

static INLINE void highbd_filter4(int8_t mask, uint8_t thresh, uint16_t *op1,
    uint16_t *op0, uint16_t *oq0, uint16_t *oq1,
    int32_t bd) {
    int16_t filter1, filter2;
    // ^0x80 equivalent to subtracting 0x80 from the values to turn them
    // into -128 to +127 instead of 0 to 255.
    int32_t shift = bd - 8;
    const int16_t ps1 = (int16_t)*op1 - (0x80 << shift);
    const int16_t ps0 = (int16_t)*op0 - (0x80 << shift);
    const int16_t qs0 = (int16_t)*oq0 - (0x80 << shift);
    const int16_t qs1 = (int16_t)*oq1 - (0x80 << shift);
    const uint16_t hev = highbd_hev_mask(thresh, *op1, *op0, *oq0, *oq1, bd);

    // Add outer taps if we have high edge variance.
    int16_t filter = signed_char_clamp_high(ps1 - qs1, bd) & hev;

    // Inner taps.
    filter = signed_char_clamp_high(filter + 3 * (qs0 - ps0), bd) & mask;

    // Save bottom 3 bits so that we round one side +4 and the other +3
    // if it equals 4 we'll set to adjust by -1 to account for the fact
    // we'd round 3 the other way.
    filter1 = signed_char_clamp_high(filter + 4, bd) >> 3;
    filter2 = signed_char_clamp_high(filter + 3, bd) >> 3;

    *oq0 = signed_char_clamp_high(qs0 - filter1, bd) + (0x80 << shift);
    *op0 = signed_char_clamp_high(ps0 + filter2, bd) + (0x80 << shift);

    // Outer tap adjustments.
    filter = ROUND_POWER_OF_TWO(filter1, 1) & ~hev;

    *oq1 = signed_char_clamp_high(qs1 - filter, bd) + (0x80 << shift);
    *op1 = signed_char_clamp_high(ps1 + filter, bd) + (0x80 << shift);
}

void aom_highbd_lpf_horizontal_4_c(uint16_t *s, int32_t p /* pitch */,
    const uint8_t *blimit, const uint8_t *limit,
    const uint8_t *thresh, int32_t bd) {
    int32_t i;
    int32_t count = 4;

    // loop filter designed to work using chars so that we can make maximum use
    // of 8 bit simd instructions.
    for (i = 0; i < count; ++i) {
        const uint16_t p1 = s[-2 * p];
        const uint16_t p0 = s[-p];
        const uint16_t q0 = s[0 * p];
        const uint16_t q1 = s[1 * p];
        const int8_t mask =
            highbd_filter_mask2(*limit, *blimit, p1, p0, q0, q1, bd);
        highbd_filter4(mask, *thresh, s - 2 * p, s - 1 * p, s, s + 1 * p, bd);
        ++s;
    }
}

void aom_highbd_lpf_vertical_4_c(uint16_t *s, int32_t pitch, const uint8_t *blimit,
    const uint8_t *limit, const uint8_t *thresh,
    int32_t bd) {
    int32_t i;
    int32_t count = 4;

    // loop filter designed to work using chars so that we can make maximum use
    // of 8 bit simd instructions.
    for (i = 0; i < count; ++i) {
        const uint16_t p1 = s[-2], p0 = s[-1];
        const uint16_t q0 = s[0], q1 = s[1];
        const int8_t mask =
            highbd_filter_mask2(*limit, *blimit, p1, p0, q0, q1, bd);
        highbd_filter4(mask, *thresh, s - 2, s - 1, s, s + 1, bd);
        s += pitch;
    }
}

static INLINE void highbd_filter8(int8_t mask, uint8_t thresh, int8_t flat,
    uint16_t *op3, uint16_t *op2, uint16_t *op1,
    uint16_t *op0, uint16_t *oq0, uint16_t *oq1,
    uint16_t *oq2, uint16_t *oq3, int32_t bd) {
    if (flat && mask) {
        const uint16_t p3 = *op3, p2 = *op2, p1 = *op1, p0 = *op0;
        const uint16_t q0 = *oq0, q1 = *oq1, q2 = *oq2, q3 = *oq3;

        // 7-tap filter [1, 1, 1, 2, 1, 1, 1]
        *op2 = ROUND_POWER_OF_TWO(p3 + p3 + p3 + 2 * p2 + p1 + p0 + q0, 3);
        *op1 = ROUND_POWER_OF_TWO(p3 + p3 + p2 + 2 * p1 + p0 + q0 + q1, 3);
        *op0 = ROUND_POWER_OF_TWO(p3 + p2 + p1 + 2 * p0 + q0 + q1 + q2, 3);
        *oq0 = ROUND_POWER_OF_TWO(p2 + p1 + p0 + 2 * q0 + q1 + q2 + q3, 3);
        *oq1 = ROUND_POWER_OF_TWO(p1 + p0 + q0 + 2 * q1 + q2 + q3 + q3, 3);
        *oq2 = ROUND_POWER_OF_TWO(p0 + q0 + q1 + 2 * q2 + q3 + q3 + q3, 3);
    }
    else
        highbd_filter4(mask, thresh, op1, op0, oq0, oq1, bd);
}

void aom_highbd_lpf_horizontal_8_c(uint16_t *s, int32_t p, const uint8_t *blimit,
    const uint8_t *limit, const uint8_t *thresh,
    int32_t bd) {
    int32_t i;
    int32_t count = 4;

    // loop filter designed to work using chars so that we can make maximum use
    // of 8 bit simd instructions.
    for (i = 0; i < count; ++i) {
        const uint16_t p3 = s[-4 * p], p2 = s[-3 * p], p1 = s[-2 * p], p0 = s[-p];
        const uint16_t q0 = s[0 * p], q1 = s[1 * p], q2 = s[2 * p], q3 = s[3 * p];

        const int8_t mask =
            highbd_filter_mask(*limit, *blimit, p3, p2, p1, p0, q0, q1, q2, q3, bd);
        const int8_t flat =
            highbd_flat_mask4(1, p3, p2, p1, p0, q0, q1, q2, q3, bd);
        highbd_filter8(mask, *thresh, flat, s - 4 * p, s - 3 * p, s - 2 * p,
            s - 1 * p, s, s + 1 * p, s + 2 * p, s + 3 * p, bd);
        ++s;
    }
}

void aom_highbd_lpf_vertical_8_c(uint16_t *s, int32_t pitch, const uint8_t *blimit,
    const uint8_t *limit, const uint8_t *thresh,
    int32_t bd) {
    int32_t i;
    int32_t count = 4;

    for (i = 0; i < count; ++i) {
        const uint16_t p3 = s[-4], p2 = s[-3], p1 = s[-2], p0 = s[-1];
        const uint16_t q0 = s[0], q1 = s[1], q2 = s[2], q3 = s[3];
        const int8_t mask =
            highbd_filter_mask(*limit, *blimit, p3, p2, p1, p0, q0, q1, q2, q3, bd);
        const int8_t flat =
            highbd_flat_mask4(1, p3, p2, p1, p0, q0, q1, q2, q3, bd);
        highbd_filter8(mask, *thresh, flat, s - 4, s - 3, s - 2, s - 1, s, s + 1,
            s + 2, s + 3, bd);
        s += pitch;
    }
}

//**********************************************************************************************************************//

//static const SEG_LVL_FEATURES seg_lvl_lf_lut[MAX_MB_PLANE][2] = {
//    { SEG_LVL_ALT_LF_Y_V, SEG_LVL_ALT_LF_Y_H },
//    { SEG_LVL_ALT_LF_U, SEG_LVL_ALT_LF_U },
//    { SEG_LVL_ALT_LF_V, SEG_LVL_ALT_LF_V }
//};

const int32_t mode_lf_lut[] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // INTRA_MODES
    1, 1, 0, 1,                             // INTER_MODES (GLOBALMV == 0)
    1, 1, 1, 1, 1, 1, 0, 1  // INTER_COMPOUND_MODES (GLOBAL_GLOBALMV == 0)
};

void update_sharpness(LoopFilterInfoN *lfi, int32_t sharpness_lvl) {
    int32_t lvl;

    // For each possible value for the loop filter fill out limits
    for (lvl = 0; lvl <= MAX_LOOP_FILTER; lvl++) {
        // Set loop filter parameters that control sharpness.
        int32_t block_inside_limit = lvl >> ((sharpness_lvl > 0) + (sharpness_lvl > 4));

        if (sharpness_lvl > 0) {
            if (block_inside_limit > (9 - sharpness_lvl))
                block_inside_limit = (9 - sharpness_lvl);
        }

        if (block_inside_limit < 1) block_inside_limit = 1;

        memset(lfi->lfthr[lvl].lim, block_inside_limit, SIMD_WIDTH);
        memset(lfi->lfthr[lvl].mblim, (2 * (lvl + 2) + block_inside_limit),
            SIMD_WIDTH);
    }
}

static int seg_feature_active(SegmentationParams *seg, int segment_id,
    SEG_LVL_FEATURES feature_id)
{
    return seg->segmentation_enabled && seg->feature_enabled[segment_id][feature_id];
}

uint8_t get_filter_level(
    FrameHeader* frm_hdr,
    const LoopFilterInfoN *lfi_n,
    const int32_t dir_idx, int32_t plane,
    int32_t *sb_delta_lf, uint8_t seg_id,
    PredictionMode pred_mode, MvReferenceFrame ref_frame_0)
{
    const int32_t segment_id = seg_id; /* const int32_t segment_id =  0; might cause encoder problem */
    PredictionMode mode; // Added to address 4x4 problem
    mode = (pred_mode == INTRA_MODE_4x4) ? DC_PRED : pred_mode;
    if (frm_hdr->delta_lf_params.delta_lf_present) {
        SVT_LOG("ERROR[AN]: delta_lf_present not supported yet\n");
        int32_t delta_lf = -1;
        if (frm_hdr->delta_lf_params.delta_lf_multi) {
            const int32_t delta_lf_idx = delta_lf_id_lut[plane][dir_idx];
            delta_lf = sb_delta_lf[delta_lf_idx];
        }
        else {
            delta_lf = sb_delta_lf[0];
        }
        int32_t base_level;
        if (plane == 0)
            base_level = frm_hdr->loop_filter_params.filter_level[dir_idx];
        else if (plane == 1)
            base_level = frm_hdr->loop_filter_params.filter_level_u;
        else
            base_level = frm_hdr->loop_filter_params.filter_level_v;
        int32_t lvl_seg = clamp(delta_lf + base_level, 0, MAX_LOOP_FILTER);
        assert(plane >= 0 && plane <= 2);
        const int32_t seg_lf_feature_id = seg_lvl_lf_lut[plane][dir_idx];
        if (seg_feature_active(&frm_hdr->segmentation_params, segment_id,
            seg_lf_feature_id))
        {
            const int32_t data = get_segdata(&frm_hdr->segmentation_params,
                                             segment_id, seg_lf_feature_id);
            lvl_seg = clamp(lvl_seg + data, 0, MAX_LOOP_FILTER);
        }

        if (frm_hdr->loop_filter_params.mode_ref_delta_enabled) {
            const int32_t scale = 1 << (lvl_seg >> 5);
            lvl_seg += frm_hdr->loop_filter_params.ref_deltas[ref_frame_0] * scale;
            if (ref_frame_0 > INTRA_FRAME)
                lvl_seg += frm_hdr->loop_filter_params.
                    mode_deltas[mode_lf_lut[mode]] * scale;
            lvl_seg = clamp(lvl_seg, 0, MAX_LOOP_FILTER);
        }
        return lvl_seg;
    }
    else {
        ASSERT(mode < MB_MODE_COUNT);
        return lfi_n->lvl[plane][segment_id][dir_idx][ref_frame_0]
            [mode_lf_lut[mode]];
    }
}

void eb_av1_loop_filter_init(PictureControlSet *pcs_ptr) {
    //assert(MB_MODE_COUNT == n_elements(mode_lf_lut));
    LoopFilterInfoN *lfi = &pcs_ptr->parent_pcs_ptr->lf_info;
    struct LoopFilter *lf = &pcs_ptr->parent_pcs_ptr->frm_hdr.loop_filter_params;
    int32_t lvl;

    lf->combine_vert_horz_lf = 1;

    // init limits for given sharpness
    update_sharpness(lfi, lf->sharpness_level);

    // init hev threshold const vectors
    for (lvl = 0; lvl <= MAX_LOOP_FILTER; lvl++)
        memset(lfi->lfthr[lvl].hev_thr, (lvl >> 4), SIMD_WIDTH);
}

// Update the loop filter for the current frame.
// This should be called before loop_filter_rows(),
// eb_av1_loop_filter_frame() calls this function directly.
void eb_av1_loop_filter_frame_init(FrameHeader *frm_hdr,
    LoopFilterInfoN *lfi, int32_t plane_start, int32_t plane_end)
{
    int32_t filt_lvl[MAX_MB_PLANE], filt_lvl_r[MAX_MB_PLANE];
    int32_t plane;
    int32_t seg_id;
    // n_shift is the multiplier for lf_deltas
    // the multiplier is 1 for when filter_lvl is between 0 and 31;
    // 2 when filter_lvl is between 32 and 63

    struct LoopFilter *const lf = &frm_hdr->loop_filter_params;
    // const struct segmentation *const seg = &pcs_ptr->parent_pcs_ptr->seg;

     // update sharpness limits
    update_sharpness(lfi, lf->sharpness_level);

    filt_lvl[0] = frm_hdr->loop_filter_params.filter_level[0];
    filt_lvl[1] = frm_hdr->loop_filter_params.filter_level_u;
    filt_lvl[2] = frm_hdr->loop_filter_params.filter_level_v;

    filt_lvl_r[0] = frm_hdr->loop_filter_params.filter_level[1];
    filt_lvl_r[1] = frm_hdr->loop_filter_params.filter_level_u;
    filt_lvl_r[2] = frm_hdr->loop_filter_params.filter_level_v;

    for (plane = plane_start; plane < plane_end; plane++) {
        if (plane == 0 && !filt_lvl[0] && !filt_lvl_r[0])
            break;
        else if (plane == 1 && !filt_lvl[1])
            continue;
        else if (plane == 2 && !filt_lvl[2])
            continue;

        for (seg_id = 0; seg_id < MAX_SEGMENTS; seg_id++) {
            for (int32_t dir = 0; dir < 2; ++dir) {
                int32_t lvl_seg = (dir == 0) ? filt_lvl[plane] : filt_lvl_r[plane];
                assert(plane >= 0 && plane <= 2);
                const int32_t seg_lf_feature_id = seg_lvl_lf_lut[plane][dir];
                if (seg_feature_active(&frm_hdr->segmentation_params, seg_id,
                    seg_lf_feature_id))
                {
                    const int32_t data = get_segdata(&frm_hdr->segmentation_params,
                                                     seg_id, seg_lf_feature_id);
                    lvl_seg = clamp(lvl_seg + data, 0, MAX_LOOP_FILTER);
                }

                if (!lf->mode_ref_delta_enabled) {
                    // we could get rid of this if we assume that deltas are set to
                    // zero when not in use; encoder always uses deltas
                    memset(lfi->lvl[plane][seg_id][dir], lvl_seg,
                        sizeof(lfi->lvl[plane][seg_id][dir]));
                }
                else {
                    int32_t ref, mode;
                    const int32_t scale = 1 << (lvl_seg >> 5);
                    const int32_t intra_lvl = lvl_seg + lf->ref_deltas[INTRA_FRAME] * scale;
                    lfi->lvl[plane][seg_id][dir][INTRA_FRAME][0] =
                        (uint8_t)clamp(intra_lvl, 0, MAX_LOOP_FILTER);

                    for (ref = LAST_FRAME; ref < REF_FRAMES; ++ref) {
                        for (mode = 0; mode < MAX_MODE_LF_DELTAS; ++mode) {
                            const int32_t inter_lvl = lvl_seg + lf->ref_deltas[ref] * scale +
                                lf->mode_deltas[mode] * scale;
                            lfi->lvl[plane][seg_id][dir][ref][mode] =
                                (uint8_t)clamp(inter_lvl, 0, MAX_LOOP_FILTER);
                        }
                    }
                }
            }
        }
    }
}
//***************************************************************************************************//

static INLINE int32_t scaled_buffer_offset(int32_t x_offset, int32_t y_offset, int32_t stride/*,
    const struct scale_factors *sf*/) {
    const int32_t x =
        /*sf ? sf->scale_value_x(x_offset, sf) >> SCALE_EXTRA_BITS :*/ x_offset;
    const int32_t y =
        /*sf ? sf->scale_value_y(y_offset, sf) >> SCALE_EXTRA_BITS :*/ y_offset;
    return y * stride + x;
}
static INLINE void setup_pred_plane(struct Buf2d *dst, BlockSize bsize,
    uint8_t *src, int32_t width, int32_t height,
    int32_t stride, int32_t mi_row, int32_t mi_col,
    /*const struct scale_factors *scale,*/
    int32_t subsampling_x, int32_t subsampling_y,
    int32_t is16Bit) {
    // Offset the buffer pointer
    if (subsampling_y && (mi_row & 0x01) && (mi_size_high[bsize] == 1))
        mi_row -= 1;
    if (subsampling_x && (mi_col & 0x01) && (mi_size_wide[bsize] == 1))
        mi_col -= 1;

    const int32_t x = (MI_SIZE * mi_col) >> subsampling_x;
    const int32_t y = (MI_SIZE * mi_row) >> subsampling_y;
    dst->buf = src + (scaled_buffer_offset(x, y, stride/*, scale*/) << is16Bit);
    dst->buf0 = src;
    dst->width = width;
    dst->height = height;
    dst->stride = stride;
}
void eb_av1_setup_dst_planes(struct MacroblockdPlane *planes, BlockSize bsize,
    //const Yv12BufferConfig *src,
    const EbPictureBufferDesc *src,
    int32_t mi_row, int32_t mi_col,
    const int32_t plane_start, const int32_t plane_end) {
    // We use AOMMIN(num_planes, MAX_MB_PLANE) instead of num_planes to quiet
    // the static analysis warnings.
    //for (int32_t i = plane_start; i < AOMMIN(plane_end, MAX_MB_PLANE); ++i) {
    //    struct MacroblockdPlane *const pd = &planes[i];
    //    const int32_t is_uv = i > 0;
    //    setup_pred_plane(&pd->dst, bsize, src->buffers[i], src->crop_widths[is_uv],
    //        src->crop_heights[is_uv], src->strides[is_uv], mi_row,
    //        mi_col, NULL, pd->subsampling_x, pd->subsampling_y);
    //}
    for (int32_t i = plane_start; i < AOMMIN(plane_end, 3); ++i) {
        if (i == 0) {
            struct MacroblockdPlane *const pd = &planes[0];
            setup_pred_plane(&pd->dst, bsize, &src->buffer_y[(src->origin_x + src->origin_y*src->stride_y) << pd->is16Bit], src->width,
                src->height, src->stride_y, mi_row,
                mi_col, /*NULL,*/ pd->subsampling_x, pd->subsampling_y, pd->is16Bit); //AMIR: Updated to point to the right location
        }
        else if (i == 1) {
            struct MacroblockdPlane *const pd = &planes[1];
            setup_pred_plane(&pd->dst, bsize, &src->buffer_cb[((src->origin_x + src->origin_y*src->stride_cb) << pd->is16Bit) / 2], src->width / 2,
                src->height / 2, src->stride_cb, mi_row,
                mi_col, /*NULL,*/ pd->subsampling_x, pd->subsampling_y, pd->is16Bit);
        }
        else if (i == 2) {
            struct MacroblockdPlane *const pd = &planes[2];
            setup_pred_plane(&pd->dst, bsize, &src->buffer_cr[((src->origin_x + src->origin_y*src->stride_cr) << pd->is16Bit) / 2], src->width / 2,
                src->height / 2, src->stride_cr, mi_row,
                mi_col,/* NULL,*/ pd->subsampling_x, pd->subsampling_y, pd->is16Bit);
        }
    }
}

#define INTER_TX_SIZE_BUF_LEN 16

//***************************************************************************************************//

static TxSize get_transform_size(const MacroBlockD *const xd,
    const MbModeInfo *const mbmi,
    const EDGE_DIR edge_dir, const int32_t mi_row,
    const int32_t mi_col, const int32_t plane,
    const struct MacroblockdPlane *plane_ptr) {
    assert(mbmi != NULL);
    (void)mi_row;
    (void)mi_col;
    (void)xd;
    //if (xd->lossless[mbmi->segment_id]) return TX_4X4;

    TxSize tx_size = (plane == COMPONENT_LUMA)
        ? (is_inter_block_no_intrabc(mbmi->block_mi.ref_frame[0])
        ? tx_depth_to_tx_size[0][mbmi->block_mi.sb_type]
        : tx_depth_to_tx_size[mbmi->tx_depth][mbmi->block_mi.sb_type]) // use max_tx_size
        : av1_get_max_uv_txsize(mbmi->block_mi.sb_type,
            plane_ptr->subsampling_x, plane_ptr->subsampling_y);
    assert(tx_size < TX_SIZES_ALL);
    if (((plane == COMPONENT_LUMA) &&
        is_inter_block_no_intrabc(mbmi->block_mi.ref_frame[0]) &&
        !mbmi->block_mi.skip)) {  // if split tx is used

        const TxSize mb_tx_size =
            tx_depth_to_tx_size[mbmi->tx_depth][mbmi->block_mi.sb_type]; // tx_size
        assert(mb_tx_size < TX_SIZES_ALL);
        tx_size = mb_tx_size;
    }
    // since in case of chrominance or non-square transorm need to convert
    // transform size into transform size in particular direction.
    // for vertical edge, filter direction is horizontal, for horizontal
    // edge, filter direction is vertical.
    tx_size = (VERT_EDGE == edge_dir) ? txsize_horz_map[tx_size]
        : txsize_vert_map[tx_size];
    return tx_size;
}

// Return TxSize from get_transform_size(), so it is plane and direction
// awared
static TxSize set_lpf_parameters(
    AV1_DEBLOCKING_PARAMETERS *const params, const uint64_t mode_step,
    const PictureControlSet *const  pcs_ptr, const MacroBlockD *const xd,
    const EDGE_DIR edge_dir, const uint32_t x, const uint32_t y,
    const int32_t plane, const struct MacroblockdPlane *const plane_ptr) {
    // reset to initial values
    params->filter_length = 0;

    // no deblocking is required
    const uint32_t width = plane_ptr->dst.width;
    const uint32_t height = plane_ptr->dst.height;
    if ((width <= x) || (height <= y)) {
        // just return the smallest transform unit size
        return TX_4X4;
    }

    const uint32_t scale_horz = plane_ptr->subsampling_x;
    const uint32_t scale_vert = plane_ptr->subsampling_y;
    // for sub8x8 block, chroma prediction mode is obtained from the bottom/right
    // mi structure of the co-located 8x8 luma block. so for chroma plane, mi_row
    // and mi_col should map to the bottom/right mi structure, i.e, both mi_row
    // and mi_col should be odd number for chroma plane.

    const int32_t mi_row = scale_vert | ((y << scale_vert) >> MI_SIZE_LOG2);
    const int32_t mi_col = scale_horz | ((x << scale_horz) >> MI_SIZE_LOG2);
    uint32_t mi_stride = pcs_ptr->mi_stride;
    const int32_t offset = mi_row * mi_stride + mi_col;
    ModeInfo **mi = (pcs_ptr->mi_grid_base + offset);
    //MbModeInfo **mi = cm->mi_grid_visible + mi_row * cm->mi_stride + mi_col;
    const MbModeInfo *mbmi = &mi[0]->mbmi;

    // If current mbmi is not correctly setup, return an invalid value to stop
    // filtering. One example is that if this tile is not coded, then its mbmi
    // it not set up.
    if (mbmi == NULL) return TX_INVALID;

    const TxSize ts =
        get_transform_size(xd, mbmi/*mi[0]*/, edge_dir, mi_row, mi_col, plane, plane_ptr);
    assert(ts < TX_SIZES_ALL);

    {
        const uint32_t coord = (VERT_EDGE == edge_dir) ? (x) : (y);
        const uint32_t transform_masks =
            edge_dir == VERT_EDGE ? tx_size_wide[ts] - 1 : tx_size_high[ts] - 1;
        const int32_t tu_edge = (coord & transform_masks) ? (0) : (1);

        if (!tu_edge) return ts;

        // prepare outer edge parameters. deblock the edge if it's an edge of a TU
        {
            const uint32_t curr_level =
                get_filter_level(&pcs_ptr->parent_pcs_ptr->frm_hdr,
                    &pcs_ptr->parent_pcs_ptr->lf_info, edge_dir, plane,
                    pcs_ptr->parent_pcs_ptr->curr_delta_lf, 0 /*segment_id*/,
                    mbmi->block_mi.mode, mbmi->block_mi.ref_frame[0]);

            const int32_t curr_skipped = mbmi->block_mi.skip &&
                is_inter_block_no_intrabc(mbmi->block_mi.ref_frame[0]);
            uint32_t level = curr_level;
            if (coord) {
                {
                    //const ModeInfo *const mi_prev = *(mi - mode_step);
                    const ModeInfo *const mi_prevTemp = *(mi - mode_step);
                    const MbModeInfo *const mi_prev = &mi_prevTemp[0].mbmi;
                    //
                    if (mi_prev == NULL) return TX_INVALID;
                    const int32_t pv_row =
                        (VERT_EDGE == edge_dir) ? (mi_row) : (mi_row - (1 << scale_vert));
                    const int32_t pv_col =
                        (VERT_EDGE == edge_dir) ? (mi_col - (1 << scale_horz)) : (mi_col);
                    const TxSize pv_ts = get_transform_size(
                        xd, mi_prev, edge_dir, pv_row, pv_col, plane, plane_ptr);
                    const uint32_t pv_lvl =
                        get_filter_level(&pcs_ptr->parent_pcs_ptr->frm_hdr,
                            &pcs_ptr->parent_pcs_ptr->lf_info, edge_dir, plane,
                            pcs_ptr->parent_pcs_ptr->curr_delta_lf, 0 /*segment_id*/,
                            mi_prev->block_mi.mode, mi_prev->block_mi.ref_frame[0]);

                    const int32_t pv_skip = mi_prev->block_mi.skip &&
                        is_inter_block_no_intrabc(mi_prev->block_mi.ref_frame[0]);

                    const BlockSize bsize =
                        get_plane_block_size(mbmi->block_mi.sb_type, plane_ptr->subsampling_x, plane_ptr->subsampling_y);
                    assert(bsize < BlockSizeS_ALL);
                    const int32_t prediction_masks = edge_dir == VERT_EDGE
                        ? block_size_wide[bsize] - 1
                        : block_size_high[bsize] - 1;
                    const int32_t pu_edge = !(coord & prediction_masks);
                    // if the current and the previous blocks are skipped,
                    // deblock the edge if the edge belongs to a PU's edge only.
                    if ((curr_level || pv_lvl) &&
                        (!pv_skip || !curr_skipped || pu_edge)) {
                        const TxSize min_ts = AOMMIN(ts, pv_ts);
                        if (TX_4X4 >= min_ts)
                            params->filter_length = 4;
                        else if (TX_8X8 == min_ts) {
                            if (plane != 0)
                                params->filter_length = 6;
                            else
                                params->filter_length = 8;
                        }
                        else {
                            params->filter_length = 14;
                            // No wide filtering for chroma plane
                            if (plane != 0)
                                params->filter_length = 6;
                        }

                        // update the level if the current block is skipped,
                        // but the previous one is not
                        level = (curr_level) ? (curr_level) : (pv_lvl);
                    }
                }
            }
            // prepare common parameters
            if (params->filter_length) {
                const LoopFilterThresh *const limits = pcs_ptr->parent_pcs_ptr->lf_info.lfthr + level;
                params->lim = limits->lim;
                params->mblim = limits->mblim;
                params->hev_thr = limits->hev_thr;
            }
        }
    }

    return ts;
}

void eb_av1_filter_block_plane_vert(
    const PictureControlSet *const  pcs_ptr,
    const MacroBlockD *const xd, const int32_t plane,
    const MacroblockdPlane *const plane_ptr,
    const uint32_t mi_row, const uint32_t mi_col) {
    SequenceControlSet *scs_ptr = (SequenceControlSet*)pcs_ptr->parent_pcs_ptr->sequence_control_set_wrapper_ptr->object_ptr;
    EbBool is16bit = scs_ptr->static_config.encoder_bit_depth > 8;
    const int32_t row_step = MI_SIZE >> MI_SIZE_LOG2;
    const uint32_t scale_horz = plane_ptr->subsampling_x;
    const uint32_t scale_vert = plane_ptr->subsampling_y;
    uint8_t *const dst_ptr = plane_ptr->dst.buf;
    const int32_t dst_stride = plane_ptr->dst.stride;
    const int32_t y_range = scs_ptr->seq_header.sb_size == BLOCK_128X128 ? (MAX_MIB_SIZE >> scale_vert) : (SB64_MIB_SIZE >> scale_vert);
    const int32_t x_range = scs_ptr->seq_header.sb_size == BLOCK_128X128 ? (MAX_MIB_SIZE >> scale_horz) : (SB64_MIB_SIZE >> scale_horz);
    for (int32_t y = 0; y < y_range; y += row_step) {
        uint8_t *p = dst_ptr + ((y * MI_SIZE * dst_stride) << plane_ptr->is16Bit);
        for (int32_t x = 0; x < x_range;) {
            // inner loop always filter vertical edges in a MI block. If MI size
            // is 8x8, it will filter the vertical edge aligned with a 8x8 block.
            // If 4x4 trasnform is used, it will then filter the internal edge
            //  aligned with a 4x4 block
            const uint32_t curr_x = ((mi_col * MI_SIZE) >> scale_horz) + x * MI_SIZE;
            const uint32_t curr_y = ((mi_row * MI_SIZE) >> scale_vert) + y * MI_SIZE;
            uint32_t advance_units;
            TxSize tx_size;
            AV1_DEBLOCKING_PARAMETERS params;
            memset(&params, 0, sizeof(params));

            tx_size =
                set_lpf_parameters(&params, ((uint64_t)1 << scale_horz), pcs_ptr, xd,
                    VERT_EDGE, curr_x, curr_y, plane, plane_ptr);
            if (tx_size == TX_INVALID) {
                params.filter_length = 0;
                tx_size = TX_4X4;
            }

            switch (params.filter_length) {
                // apply 4-tap filtering
            case 4:
                if (is16bit)
                    aom_highbd_lpf_vertical_4(
                    (uint16_t*)(p),//CONVERT_TO_SHORTPTR(p),
                        dst_stride,
                        params.mblim,
                        params.lim,
                        params.hev_thr,
                        scs_ptr->static_config.encoder_bit_depth);
                else
                    aom_lpf_vertical_4(
                        p,
                        dst_stride,
                        params.mblim,
                        params.lim,
                        params.hev_thr);
                break;
            case 6:  // apply 6-tap filter for chroma plane only
                assert(plane != 0);
                if (is16bit)
                    aom_highbd_lpf_vertical_6(
                    (uint16_t*)(p),//CONVERT_TO_SHORTPTR(p),
                        dst_stride,
                        params.mblim,
                        params.lim,
                        params.hev_thr,
                        scs_ptr->static_config.encoder_bit_depth);
                else
                    aom_lpf_vertical_6(
                        p,
                        dst_stride,
                        params.mblim,
                        params.lim,
                        params.hev_thr);
                break;
                // apply 8-tap filtering
            case 8:
                if (is16bit)
                    aom_highbd_lpf_vertical_8(
                    (uint16_t*)(p),//CONVERT_TO_SHORTPTR(p),
                        dst_stride,
                        params.mblim,
                        params.lim,
                        params.hev_thr,
                        scs_ptr->static_config.encoder_bit_depth);
                else
                    aom_lpf_vertical_8(
                        p,
                        dst_stride,
                        params.mblim,
                        params.lim,
                        params.hev_thr);
                break;
                // apply 14-tap filtering
            case 14:
                if (is16bit)
                    aom_highbd_lpf_vertical_14(
                    (uint16_t*)(p),//CONVERT_TO_SHORTPTR(p),
                        dst_stride,
                        params.mblim,
                        params.lim,
                        params.hev_thr,
                        scs_ptr->static_config.encoder_bit_depth);
                else
                    aom_lpf_vertical_14(
                        p,
                        dst_stride,
                        params.mblim,
                        params.lim,
                        params.hev_thr);
                break;
                // no filtering
            default: break;
            }
            // advance the destination pointer
            assert(tx_size < TX_SIZES_ALL);
            advance_units = tx_size_wide_unit[tx_size];
            x += advance_units;
            p += ((advance_units * MI_SIZE) << plane_ptr->is16Bit);
        }
    }
}

void eb_av1_filter_block_plane_horz(
    const PictureControlSet *const  pcs_ptr,
    const MacroBlockD *const xd, const int32_t plane,
    const MacroblockdPlane *const plane_ptr,
    const uint32_t mi_row, const uint32_t mi_col) {
    SequenceControlSet *scs_ptr = (SequenceControlSet*)pcs_ptr->parent_pcs_ptr->sequence_control_set_wrapper_ptr->object_ptr;
    EbBool is16bit = scs_ptr->static_config.encoder_bit_depth > 8;
    const int32_t col_step = MI_SIZE >> MI_SIZE_LOG2;
    const uint32_t scale_horz = plane_ptr->subsampling_x;
    const uint32_t scale_vert = plane_ptr->subsampling_y;
    uint8_t *const dst_ptr = plane_ptr->dst.buf;
    const int32_t dst_stride = plane_ptr->dst.stride;
    const int32_t y_range = scs_ptr->seq_header.sb_size == BLOCK_128X128 ? (MAX_MIB_SIZE >> scale_vert) : (SB64_MIB_SIZE >> scale_vert);
    const int32_t x_range = scs_ptr->seq_header.sb_size == BLOCK_128X128 ? (MAX_MIB_SIZE >> scale_horz) : (SB64_MIB_SIZE >> scale_horz);
    uint32_t mi_stride = pcs_ptr->mi_stride;
    for (int32_t x = 0; x < x_range; x += col_step) {
        uint8_t *p = dst_ptr + ((x * MI_SIZE) << plane_ptr->is16Bit);
        for (int32_t y = 0; y < y_range;) {
            // inner loop always filter vertical edges in a MI block. If MI size
            // is 8x8, it will first filter the vertical edge aligned with a 8x8
            // block. If 4x4 trasnform is used, it will then filter the internal
            // edge aligned with a 4x4 block
            const uint32_t curr_x = ((mi_col * MI_SIZE) >> scale_horz) + x * MI_SIZE;
            const uint32_t curr_y = ((mi_row * MI_SIZE) >> scale_vert) + y * MI_SIZE;
            uint32_t advance_units;
            TxSize tx_size;
            AV1_DEBLOCKING_PARAMETERS params;
            memset(&params, 0, sizeof(params));

            tx_size =
                set_lpf_parameters(
                    &params,
                    //(pcs_ptr->parent_pcs_ptr->av1_cm->mi_stride << scale_vert),
                    (mi_stride << scale_vert),
                    pcs_ptr,
                    xd,
                    HORZ_EDGE,
                    curr_x,
                    curr_y,
                    plane,
                    plane_ptr);
            if (tx_size == TX_INVALID) {
                params.filter_length = 0;
                tx_size = TX_4X4;
            }

            switch (params.filter_length) {
                // apply 4-tap filtering
            case 4:
                if (is16bit)
                    aom_highbd_lpf_horizontal_4(
                    (uint16_t*)(p),//CONVERT_TO_SHORTPTR(p),
                        dst_stride,
                        params.mblim,
                        params.lim,
                        params.hev_thr,
                        scs_ptr->static_config.encoder_bit_depth);
                else
                    aom_lpf_horizontal_4(
                        p,
                        dst_stride,
                        params.mblim,
                        params.lim,
                        params.hev_thr);
                break;
                // apply 6-tap filtering
            case 6:
                assert(plane != 0);
                if (is16bit)
                    aom_highbd_lpf_horizontal_6(
                    (uint16_t*)(p),//CONVERT_TO_SHORTPTR(p),
                        dst_stride,
                        params.mblim,
                        params.lim,
                        params.hev_thr,
                        scs_ptr->static_config.encoder_bit_depth);
                else
                    aom_lpf_horizontal_6(
                        p,
                        dst_stride,
                        params.mblim,
                        params.lim,
                        params.hev_thr);
                break;
                // apply 8-tap filtering
            case 8:
                if (is16bit)
                    aom_highbd_lpf_horizontal_8(
                    (uint16_t*)(p),//CONVERT_TO_SHORTPTR(p),
                        dst_stride,
                        params.mblim,
                        params.lim,
                        params.hev_thr,
                        scs_ptr->static_config.encoder_bit_depth);
                else
                    aom_lpf_horizontal_8(
                        p,
                        dst_stride,
                        params.mblim,
                        params.lim,
                        params.hev_thr);
                break;
                // apply 14-tap filtering
            case 14:
                if (is16bit)
                    aom_highbd_lpf_horizontal_14(
                    (uint16_t*)(p),//CONVERT_TO_SHORTPTR(p),
                        dst_stride,
                        params.mblim,
                        params.lim,
                        params.hev_thr,
                        scs_ptr->static_config.encoder_bit_depth);
                else
                    aom_lpf_horizontal_14(
                        p,
                        dst_stride,
                        params.mblim,
                        params.lim,
                        params.hev_thr);
                break;
                // no filtering
            default: break;
            }

            // advance the destination pointer
            assert(tx_size < TX_SIZES_ALL);
            advance_units = tx_size_high_unit[tx_size];
            y += advance_units;
            p += ((advance_units * dst_stride * MI_SIZE) << plane_ptr->is16Bit);
        }
    }
}

// New function to filter each sb (64x64)
void loop_filter_sb(
    EbPictureBufferDesc *frame_buffer,//reconpicture,
    //Yv12BufferConfig *frame_buffer,
    PictureControlSet *pcs_ptr,
    MacroBlockD *xd, int32_t mi_row, int32_t mi_col,
    int32_t plane_start, int32_t plane_end,
    uint8_t LastCol) {
    FrameHeader *frm_hdr = &pcs_ptr->parent_pcs_ptr->frm_hdr;
    struct MacroblockdPlane pd[3];
    int32_t plane;

    pd[0].subsampling_x = 0;
    pd[0].subsampling_y = 0;
    pd[0].plane_type = PLANE_TYPE_Y;
    pd[0].is16Bit = frame_buffer->bit_depth > 8;
    pd[1].subsampling_x = 1;
    pd[1].subsampling_y = 1;
    pd[1].plane_type = PLANE_TYPE_UV;
    pd[1].is16Bit = frame_buffer->bit_depth > 8;
    pd[2].subsampling_x = 1;
    pd[2].subsampling_y = 1;
    pd[2].plane_type = PLANE_TYPE_UV;
    pd[2].is16Bit = frame_buffer->bit_depth > 8;

    for (plane = plane_start; plane < plane_end; plane++) {
        if (plane == 0 && !(frm_hdr->loop_filter_params.filter_level[0]) && !(frm_hdr->loop_filter_params.filter_level[1]))
            break;
        else if (plane == 1 && !(frm_hdr->loop_filter_params.filter_level_u))
            continue;
        else if (plane == 2 && !(frm_hdr->loop_filter_params.filter_level_v))
            continue;

        if (frm_hdr->loop_filter_params.combine_vert_horz_lf) {
            // filter all vertical and horizontal edges in every 64x64 super block
            // filter vertical edges
            eb_av1_setup_dst_planes(pd, pcs_ptr->parent_pcs_ptr->sequence_control_set_ptr->seq_header.sb_size, frame_buffer, mi_row,
                mi_col, plane, plane + 1);
            eb_av1_filter_block_plane_vert(pcs_ptr, xd, plane, &pd[plane], mi_row,
                mi_col);
            // filter horizontal edges
            int32_t max_mib_size = pcs_ptr->parent_pcs_ptr->sequence_control_set_ptr->seq_header.sb_size == BLOCK_128X128 ? MAX_MIB_SIZE : SB64_MIB_SIZE;

            if (mi_col - max_mib_size >= 0) {
                eb_av1_setup_dst_planes(pd, pcs_ptr->parent_pcs_ptr->sequence_control_set_ptr->seq_header.sb_size, frame_buffer,
                    mi_row, mi_col - max_mib_size, plane,
                    plane + 1);
                eb_av1_filter_block_plane_horz(pcs_ptr, xd, plane, &pd[plane], mi_row,
                    mi_col - max_mib_size);
            }
            // Filter the horizontal edges of the last lcu in each row
            if (LastCol) {
                eb_av1_setup_dst_planes(pd, pcs_ptr->parent_pcs_ptr->sequence_control_set_ptr->seq_header.sb_size, frame_buffer,
                    mi_row, mi_col, plane,
                    plane + 1);
                eb_av1_filter_block_plane_horz(pcs_ptr, xd, plane, &pd[plane], mi_row,
                    mi_col);
            }
        }
        else {
            // filter all vertical edges in every 64x64 super block
            eb_av1_setup_dst_planes(pd, pcs_ptr->parent_pcs_ptr->sequence_control_set_ptr->seq_header.sb_size, frame_buffer, mi_row,
                mi_col, plane, plane + 1);

            eb_av1_filter_block_plane_vert(pcs_ptr, xd, plane, &pd[plane], mi_row,
                mi_col);

            // filter all horizontal edges in every 64x64 super block
            eb_av1_setup_dst_planes(pd, pcs_ptr->parent_pcs_ptr->sequence_control_set_ptr->seq_header.sb_size, frame_buffer, mi_row,
                mi_col, plane, plane + 1);
            eb_av1_filter_block_plane_horz(pcs_ptr, xd, plane, &pd[plane], mi_row,
                mi_col);
        }
    }
}

void eb_av1_loop_filter_frame(
    EbPictureBufferDesc *frame_buffer,
    PictureControlSet *picture_control_set_ptr,
    int32_t plane_start, int32_t plane_end) {
    SequenceControlSet *scs_ptr = (SequenceControlSet*)picture_control_set_ptr->parent_pcs_ptr->sequence_control_set_wrapper_ptr->object_ptr;
    //SuperBlock                     *sb_ptr;
    //uint16_t                                   sb_index;
    uint8_t                                   sb_size_Log2 = (uint8_t)Log2f(scs_ptr->sb_size_pix);
    uint32_t                                   x_lcu_index;
    uint32_t                                   y_lcu_index;
    uint32_t                                   sb_origin_x;
    uint32_t                                   sb_origin_y;
    EbBool                                  endOfRowFlag;

    uint32_t picture_width_in_sb = (scs_ptr->seq_header.max_frame_width + scs_ptr->sb_size_pix - 1) / scs_ptr->sb_size_pix;
    uint32_t picture_height_in_sb = (scs_ptr->seq_header.max_frame_height + scs_ptr->sb_size_pix - 1) / scs_ptr->sb_size_pix;

    eb_av1_loop_filter_frame_init(&picture_control_set_ptr->parent_pcs_ptr->frm_hdr,
        &picture_control_set_ptr->parent_pcs_ptr->lf_info, plane_start, plane_end);

    for (y_lcu_index = 0; y_lcu_index < picture_height_in_sb; ++y_lcu_index) {
        for (x_lcu_index = 0; x_lcu_index < picture_width_in_sb; ++x_lcu_index) {
            //sb_index        = (uint16_t)(y_lcu_index * picture_width_in_sb + x_lcu_index);
            //sb_ptr          = picture_control_set_ptr->sb_ptr_array[sb_index];
            sb_origin_x = x_lcu_index << sb_size_Log2;
            sb_origin_y = y_lcu_index << sb_size_Log2;
            endOfRowFlag = (x_lcu_index == picture_width_in_sb - 1) ? EB_TRUE : EB_FALSE;

            loop_filter_sb(
                frame_buffer,
                picture_control_set_ptr,
                NULL,
                sb_origin_y >> 2,
                sb_origin_x >> 2,
                plane_start,
                plane_end,
                endOfRowFlag);
        }
    }
}
extern int16_t eb_av1_ac_quant_Q3(int32_t qindex, int32_t delta, AomBitDepth bit_depth);

void EbCopyBuffer(
    EbPictureBufferDesc  *srcBuffer,
    EbPictureBufferDesc  *dstBuffer,
    PictureControlSet    *pcs_ptr,
    uint8_t                   plane) {
    EbBool is16bit = (EbBool)(pcs_ptr->parent_pcs_ptr->sequence_control_set_ptr->static_config.encoder_bit_depth > EB_8BIT);
    dstBuffer->origin_x = srcBuffer->origin_x;
    dstBuffer->origin_y = srcBuffer->origin_y;
    dstBuffer->width = srcBuffer->width;
    dstBuffer->height = srcBuffer->height;
    dstBuffer->max_width = srcBuffer->max_width;
    dstBuffer->max_height = srcBuffer->max_height;
    dstBuffer->bit_depth = srcBuffer->bit_depth;
    dstBuffer->luma_size = srcBuffer->luma_size;
    dstBuffer->chroma_size = srcBuffer->chroma_size;
    dstBuffer->packedFlag = srcBuffer->packedFlag;

    uint32_t   lumaBufferOffset = (srcBuffer->origin_x + srcBuffer->origin_y*srcBuffer->stride_y) << is16bit;
    uint16_t   luma_width = (uint16_t)(srcBuffer->width - pcs_ptr->parent_pcs_ptr->sequence_control_set_ptr->pad_right) << is16bit;
    uint16_t   luma_height = (uint16_t)(srcBuffer->height - pcs_ptr->parent_pcs_ptr->sequence_control_set_ptr->pad_bottom);
    uint16_t   chroma_width = (luma_width >> 1);
    if (plane == 0) {
        uint16_t stride_y = srcBuffer->stride_y << is16bit;

        dstBuffer->stride_y = srcBuffer->stride_y;
        dstBuffer->stride_bit_inc_y = srcBuffer->stride_bit_inc_y;

        for (int32_t inputRowIndex = 0; inputRowIndex < luma_height; inputRowIndex++) {
            EB_MEMCPY((dstBuffer->buffer_y + lumaBufferOffset + stride_y * inputRowIndex),
                (srcBuffer->buffer_y + lumaBufferOffset + stride_y * inputRowIndex),
                luma_width);
        }
    }
    else if (plane == 1) {
        uint16_t stride_cb = srcBuffer->stride_cb << is16bit;
        dstBuffer->stride_cb = srcBuffer->stride_cb;
        dstBuffer->stride_bit_inc_cb = srcBuffer->stride_bit_inc_cb;

        uint32_t   chromaBufferOffset = (srcBuffer->origin_x / 2 + srcBuffer->origin_y / 2 * srcBuffer->stride_cb) << is16bit;

        for (int32_t inputRowIndex = 0; inputRowIndex < luma_height >> 1; inputRowIndex++) {
            EB_MEMCPY((dstBuffer->buffer_cb + chromaBufferOffset + stride_cb * inputRowIndex),
                (srcBuffer->buffer_cb + chromaBufferOffset + stride_cb * inputRowIndex),
                chroma_width);
        }
    }
    else if (plane == 2) {
        uint16_t stride_cr = srcBuffer->stride_cr << is16bit;

        dstBuffer->stride_cr = srcBuffer->stride_cr;
        dstBuffer->stride_bit_inc_cr = srcBuffer->stride_bit_inc_cr;

        uint32_t   chromaBufferOffset = (srcBuffer->origin_x / 2 + srcBuffer->origin_y / 2 * srcBuffer->stride_cr) << is16bit;

        for (int32_t inputRowIndex = 0; inputRowIndex < luma_height >> 1; inputRowIndex++) {
            EB_MEMCPY((dstBuffer->buffer_cr + chromaBufferOffset + stride_cr * inputRowIndex),
                (srcBuffer->buffer_cr + chromaBufferOffset + stride_cr * inputRowIndex),
                chroma_width);
        }
    }
}

//int32_t av1_get_max_filter_level(const Av1Comp *cpi) {
//    if (cpi->oxcf.pass == 2) {
//        return cpi->twopass.section_intra_rating > 8 ? MAX_LOOP_FILTER * 3 / 4
//            : MAX_LOOP_FILTER;
//    }
//    else {
//        return MAX_LOOP_FILTER;
//    }
//}

uint64_t PictureSseCalculations(
    PictureControlSet    *picture_control_set_ptr,
    EbPictureBufferDesc *recon_ptr,
    int32_t plane)

{
    SequenceControlSet   *sequence_control_set_ptr = picture_control_set_ptr->parent_pcs_ptr->sequence_control_set_ptr;
    EbBool is16bit = (sequence_control_set_ptr->static_config.encoder_bit_depth > EB_8BIT);

    if (!is16bit) {
        EbPictureBufferDesc *input_picture_ptr = (EbPictureBufferDesc*)picture_control_set_ptr->parent_pcs_ptr->enhanced_picture_ptr;

        uint32_t   columnIndex;
        uint32_t   row_index = 0;
        uint64_t   residualDistortion = 0;
        EbByte  inputBuffer;
        EbByte  reconCoeffBuffer;
        if (plane == 0) {
            reconCoeffBuffer = &((recon_ptr->buffer_y)[recon_ptr->origin_x + recon_ptr->origin_y * recon_ptr->stride_y]);
            inputBuffer = &((input_picture_ptr->buffer_y)[input_picture_ptr->origin_x + input_picture_ptr->origin_y * input_picture_ptr->stride_y]);

            residualDistortion = 0;

            while (row_index < sequence_control_set_ptr->seq_header.max_frame_height) {
                columnIndex = 0;
                while (columnIndex < sequence_control_set_ptr->seq_header.max_frame_width) {
                    residualDistortion += (int64_t)SQR((int64_t)(inputBuffer[columnIndex]) - (reconCoeffBuffer[columnIndex]));
                    ++columnIndex;
                }
                inputBuffer += input_picture_ptr->stride_y;
                reconCoeffBuffer += recon_ptr->stride_y;
                ++row_index;
            }

            return residualDistortion;
        }

        else if (plane == 1) {
            reconCoeffBuffer = &((recon_ptr->buffer_cb)[recon_ptr->origin_x / 2 + recon_ptr->origin_y / 2 * recon_ptr->stride_cb]);
            inputBuffer = &((input_picture_ptr->buffer_cb)[input_picture_ptr->origin_x / 2 + input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cb]);

            residualDistortion = 0;
            row_index = 0;
            while (row_index < sequence_control_set_ptr->chroma_height) {
                columnIndex = 0;
                while (columnIndex < sequence_control_set_ptr->chroma_width) {
                    residualDistortion += (int64_t)SQR((int64_t)(inputBuffer[columnIndex]) - (reconCoeffBuffer[columnIndex]));
                    ++columnIndex;
                }

                inputBuffer += input_picture_ptr->stride_cb;
                reconCoeffBuffer += recon_ptr->stride_cb;
                ++row_index;
            }

            return residualDistortion;
        }
        else if (plane == 2) {
            reconCoeffBuffer = &((recon_ptr->buffer_cr)[recon_ptr->origin_x / 2 + recon_ptr->origin_y / 2 * recon_ptr->stride_cr]);
            inputBuffer = &((input_picture_ptr->buffer_cr)[input_picture_ptr->origin_x / 2 + input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cr]);
            residualDistortion = 0;
            row_index = 0;

            while (row_index < sequence_control_set_ptr->chroma_height) {
                columnIndex = 0;
                while (columnIndex < sequence_control_set_ptr->chroma_width) {
                    residualDistortion += (int64_t)SQR((int64_t)(inputBuffer[columnIndex]) - (reconCoeffBuffer[columnIndex]));
                    ++columnIndex;
                }

                inputBuffer += input_picture_ptr->stride_cr;
                reconCoeffBuffer += recon_ptr->stride_cr;
                ++row_index;
            }

            return residualDistortion;
        }
        return 0;
    }
    else {
        EbPictureBufferDesc *input_picture_ptr = (EbPictureBufferDesc*)picture_control_set_ptr->input_frame16bit;

        uint32_t   columnIndex;
        uint32_t   row_index = 0;
        uint64_t   residualDistortion = 0;
        uint16_t*  inputBuffer;
        uint16_t*  reconCoeffBuffer;
        if (plane == 0) {
            reconCoeffBuffer = (uint16_t*)&((recon_ptr->buffer_y)[(recon_ptr->origin_x + recon_ptr->origin_y * recon_ptr->stride_y) << is16bit]);
            inputBuffer = (uint16_t*)&((input_picture_ptr->buffer_y)[(input_picture_ptr->origin_x + input_picture_ptr->origin_y * input_picture_ptr->stride_y) << is16bit]);

            residualDistortion = 0;

            while (row_index < sequence_control_set_ptr->seq_header.max_frame_height) {
                columnIndex = 0;
                while (columnIndex < sequence_control_set_ptr->seq_header.max_frame_width) {
                    residualDistortion += (int64_t)SQR(((int64_t)inputBuffer[columnIndex]) - (int64_t)(reconCoeffBuffer[columnIndex]));
                    ++columnIndex;
                }

                inputBuffer += input_picture_ptr->stride_y;
                reconCoeffBuffer += recon_ptr->stride_y;
                ++row_index;
            }

            return residualDistortion;
        }

        else if (plane == 1) {
            reconCoeffBuffer = (uint16_t*)&((recon_ptr->buffer_cb)[(recon_ptr->origin_x / 2 + recon_ptr->origin_y / 2 * recon_ptr->stride_cb) << is16bit]);
            inputBuffer = (uint16_t*)&((input_picture_ptr->buffer_cb)[(input_picture_ptr->origin_x / 2 + input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cb) << is16bit]);

            residualDistortion = 0;
            row_index = 0;
            while (row_index < sequence_control_set_ptr->chroma_height) {
                columnIndex = 0;
                while (columnIndex < sequence_control_set_ptr->chroma_width) {
                    residualDistortion += (int64_t)SQR(((int64_t)inputBuffer[columnIndex]) - (int64_t)(reconCoeffBuffer[columnIndex]));
                    ++columnIndex;
                }

                inputBuffer += input_picture_ptr->stride_cb;
                reconCoeffBuffer += recon_ptr->stride_cb;
                ++row_index;
            }

            return residualDistortion;
        }
        else if (plane == 2) {
            reconCoeffBuffer = (uint16_t*)&((recon_ptr->buffer_cr)[(recon_ptr->origin_x / 2 + recon_ptr->origin_y / 2 * recon_ptr->stride_cr) << is16bit]);
            inputBuffer = (uint16_t*)&((input_picture_ptr->buffer_cr)[(input_picture_ptr->origin_x / 2 + input_picture_ptr->origin_y / 2 * input_picture_ptr->stride_cr) << is16bit]);
            residualDistortion = 0;
            row_index = 0;

            while (row_index < sequence_control_set_ptr->chroma_height) {
                columnIndex = 0;
                while (columnIndex < sequence_control_set_ptr->chroma_width) {
                    residualDistortion += (int64_t)SQR(((int64_t)inputBuffer[columnIndex]) - (int64_t)(reconCoeffBuffer[columnIndex]));
                    ++columnIndex;
                }

                inputBuffer += input_picture_ptr->stride_cr;
                reconCoeffBuffer += recon_ptr->stride_cr;
                ++row_index;
            }

            return residualDistortion;
        }

        return 0;
    }
}

static int64_t try_filter_frame(
    //const Yv12BufferConfig *sd,
    //Av1Comp *const cpi,
    const EbPictureBufferDesc *sd,
    EbPictureBufferDesc  *tempLfReconBuffer,
    PictureControlSet *pcs_ptr,
    int32_t filt_level,
    int32_t partial_frame, int32_t plane, int32_t dir) {
    (void)sd;
    (void)partial_frame;
    (void)sd;
    int64_t filt_err;
    FrameHeader *frm_hdr = &pcs_ptr->parent_pcs_ptr->frm_hdr;
    assert(plane >= 0 && plane <= 2);
    int32_t filter_level[2] = { filt_level, filt_level };
    if (plane == 0 && dir == 0) filter_level[1] = frm_hdr->loop_filter_params.filter_level[1];
    if (plane == 0 && dir == 1) filter_level[0] = frm_hdr->loop_filter_params.filter_level[0];

    EbBool is16bit = (EbBool)(pcs_ptr->parent_pcs_ptr->sequence_control_set_ptr->static_config.encoder_bit_depth > EB_8BIT);
    EbPictureBufferDesc  *recon_buffer = is16bit ? pcs_ptr->recon_picture16bit_ptr : pcs_ptr->recon_picture_ptr;
    if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE) {
        //get the 16bit form of the input LCU
        if (is16bit)
            recon_buffer = ((EbReferenceObject*)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->reference_picture16bit;
        else
            recon_buffer = ((EbReferenceObject*)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->reference_picture;
    }
    else { // non ref pictures
        recon_buffer = is16bit ? pcs_ptr->recon_picture16bit_ptr : pcs_ptr->recon_picture_ptr;
    }

    // set base filters for use of get_filter_level when in DELTA_Q_LF mode
    switch (plane) {
    case 0:
        frm_hdr->loop_filter_params.filter_level[0] = filter_level[0];
        frm_hdr->loop_filter_params.filter_level[1] = filter_level[1];
        break;
    case 1: frm_hdr->loop_filter_params.filter_level_u = filter_level[0]; break;
    case 2: frm_hdr->loop_filter_params.filter_level_v = filter_level[0]; break;
    }

    eb_av1_loop_filter_frame(recon_buffer, pcs_ptr, plane, plane + 1);

    filt_err = PictureSseCalculations(pcs_ptr, recon_buffer, plane);

    // Re-instate the unfiltered frame
    EbCopyBuffer(tempLfReconBuffer/*cpi->last_frame_uf*/, recon_buffer /*cm->frame_to_show*/, pcs_ptr, (uint8_t)plane);

    return filt_err;
}
static int32_t search_filter_level(
    //const Yv12BufferConfig *sd, Av1Comp *cpi,
    EbPictureBufferDesc *sd, // source
    EbPictureBufferDesc  *tempLfReconBuffer,
    PictureControlSet *pcs_ptr,
    int32_t partial_frame,
    const int32_t *last_frame_filter_level,
    double *best_cost_ret, int32_t plane, int32_t dir) {
    const int32_t min_filter_level = 0;
    const int32_t max_filter_level = MAX_LOOP_FILTER;// av1_get_max_filter_level(cpi);
    int32_t filt_direction = 0;
    int64_t best_err;
    int32_t filt_best;
    FrameHeader *frm_hdr = &pcs_ptr->parent_pcs_ptr->frm_hdr;
    //Macroblock *x = &cpi->td.mb;

    // Start the search at the previous frame filter level unless it is now out of
    // range.
    int32_t lvl;
    switch (plane) {
    case 0: lvl = last_frame_filter_level[dir]; break;
    case 1: lvl = last_frame_filter_level[2]; break;
    case 2: lvl = last_frame_filter_level[3]; break;
    default: assert(plane >= 0 && plane <= 2); return 0;
    }
    int32_t filt_mid = clamp(lvl, min_filter_level, max_filter_level);
    int32_t filter_step = filt_mid < 16 ? 4 : filt_mid / 4;

    EbBool is16bit = (EbBool)(pcs_ptr->parent_pcs_ptr->sequence_control_set_ptr->static_config.encoder_bit_depth > EB_8BIT);
    EbPictureBufferDesc  *recon_buffer = is16bit ? pcs_ptr->recon_picture16bit_ptr : pcs_ptr->recon_picture_ptr;

    if (pcs_ptr->parent_pcs_ptr->is_used_as_reference_flag == EB_TRUE) {
        //get the 16bit form of the input LCU
        if (is16bit)
            recon_buffer = ((EbReferenceObject*)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->reference_picture16bit;
        else
            recon_buffer = ((EbReferenceObject*)pcs_ptr->parent_pcs_ptr->reference_picture_wrapper_ptr->object_ptr)->reference_picture;
    }
    else { // non ref pictures
        recon_buffer = is16bit ? pcs_ptr->recon_picture16bit_ptr : pcs_ptr->recon_picture_ptr;
    }
    // Sum squared error at each filter level
    int64_t ss_err[MAX_LOOP_FILTER + 1];

    // Set each entry to -1
    memset(ss_err, 0xFF, sizeof(ss_err));
    // make a copy of recon_buffer
    EbCopyBuffer(recon_buffer/*cm->frame_to_show*/, tempLfReconBuffer/*&cpi->last_frame_uf*/, pcs_ptr, (uint8_t)plane);

    best_err = try_filter_frame(sd, tempLfReconBuffer, pcs_ptr, filt_mid, partial_frame, plane, dir);
    filt_best = filt_mid;
    ss_err[filt_mid] = best_err;

    if (pcs_ptr->parent_pcs_ptr->loop_filter_mode <= 2) {
        filter_step = 2;
        const int32_t filt_high = AOMMIN(filt_mid + filter_step, max_filter_level);
        const int32_t filt_low = AOMMAX(filt_mid - filter_step, min_filter_level);

        // Bias against raising loop filter in favor of lowering it.
        int64_t bias = (best_err >> (15 - (filt_mid / 8))) * filter_step;

        //if ((cpi->oxcf.pass == 2) && (cpi->twopass.section_intra_rating < 20))
        //    bias = (bias * cpi->twopass.section_intra_rating) / 20;

        // yx, bias less for large block size
        if (frm_hdr->tx_mode != ONLY_4X4) bias >>= 1;

        if (filt_direction <= 0 && filt_low != filt_mid) {
            // Get Low filter error score
            if (ss_err[filt_low] < 0) {
                ss_err[filt_low] =
                    try_filter_frame(sd, tempLfReconBuffer, pcs_ptr, filt_low, partial_frame, plane, dir);
            }
            // If value is close to the best so far then bias towards a lower loop
            // filter value.
            if (ss_err[filt_low] < (best_err + bias)) {
                // Was it actually better than the previous best?
                if (ss_err[filt_low] < best_err)
                    best_err = ss_err[filt_low];
                filt_best = filt_low;
            }
        }

        // Now look at filt_high
        if (filt_direction >= 0 && filt_high != filt_mid) {
            if (ss_err[filt_high] < 0) {
                ss_err[filt_high] =
                    try_filter_frame(sd, tempLfReconBuffer, pcs_ptr, filt_high, partial_frame, plane, dir);
            }
            // If value is significantly better than previous best, bias added against
            // raising filter value
            if (ss_err[filt_high] < (best_err - bias)) {
                best_err = ss_err[filt_high];
                filt_best = filt_high;
            }
        }
    }
    else {
        while (filter_step > 0) {
            const int32_t filt_high = AOMMIN(filt_mid + filter_step, max_filter_level);
            const int32_t filt_low = AOMMAX(filt_mid - filter_step, min_filter_level);

            // Bias against raising loop filter in favor of lowering it.
            int64_t bias = (best_err >> (15 - (filt_mid / 8))) * filter_step;

            //if ((cpi->oxcf.pass == 2) && (cpi->twopass.section_intra_rating < 20))
            //    bias = (bias * cpi->twopass.section_intra_rating) / 20;

            // yx, bias less for large block size
            if (frm_hdr->tx_mode != ONLY_4X4) bias >>= 1;

            if (filt_direction <= 0 && filt_low != filt_mid) {
                // Get Low filter error score
                if (ss_err[filt_low] < 0) {
                    ss_err[filt_low] =
                        try_filter_frame(sd, tempLfReconBuffer, pcs_ptr, filt_low, partial_frame, plane, dir);
                }
                // If value is close to the best so far then bias towards a lower loop
                // filter value.
                if (ss_err[filt_low] < (best_err + bias)) {
                    // Was it actually better than the previous best?
                    if (ss_err[filt_low] < best_err)
                        best_err = ss_err[filt_low];
                    filt_best = filt_low;
                }
            }

            // Now look at filt_high
            if (filt_direction >= 0 && filt_high != filt_mid) {
                if (ss_err[filt_high] < 0) {
                    ss_err[filt_high] =
                        try_filter_frame(sd, tempLfReconBuffer, pcs_ptr, filt_high, partial_frame, plane, dir);
                }
                // If value is significantly better than previous best, bias added against
                // raising filter value
                if (ss_err[filt_high] < (best_err - bias)) {
                    best_err = ss_err[filt_high];
                    filt_best = filt_high;
                }
            }

            // Half the step distance if the best filter value was the same as last time
            if (filt_best == filt_mid) {
                filter_step /= 2;
                filt_direction = 0;
            }
            else {
                filt_direction = (filt_best < filt_mid) ? -1 : 1;
                filt_mid = filt_best;
            }
        }
    }
    // Update best error
    best_err = ss_err[filt_best];

    if (best_cost_ret) *best_cost_ret = (double)best_err;//RDCOST_DBL(x->rdmult, 0, best_err);
    return filt_best;
}

void eb_av1_pick_filter_level(
    DlfContext            *context_ptr,
    EbPictureBufferDesc   *srcBuffer, // source input
    PictureControlSet     *pcs_ptr,
    LpfPickMethod          method) {
    SequenceControlSet *scs_ptr = (SequenceControlSet*)pcs_ptr->parent_pcs_ptr->sequence_control_set_wrapper_ptr->object_ptr;
    FrameHeader *frm_hdr = &pcs_ptr->parent_pcs_ptr->frm_hdr;

    const int32_t num_planes = 3;
    (void)srcBuffer;
    struct LoopFilter *const lf = &frm_hdr->loop_filter_params;
    lf->sharpness_level = frm_hdr->frame_type == KEY_FRAME ? 0 : LF_SHARPNESS;

    if (method == LPF_PICK_MINIMAL_LPF) {
        lf->filter_level[0] = 0;
        lf->filter_level[1] = 0;
    }
    else if (method >= LPF_PICK_FROM_Q) {
        const int32_t min_filter_level = 0;
        const int32_t max_filter_level = MAX_LOOP_FILTER;// av1_get_max_filter_level(cpi);
        const int32_t q = eb_av1_ac_quant_Q3(frm_hdr->quantization_params.base_q_idx, 0, (AomBitDepth)scs_ptr->static_config.encoder_bit_depth);
        // These values were determined by linear fitting the result of the
        // searched level for 8 bit depth:
        // Keyframes: filt_guess = q * 0.06699 - 1.60817
        // Other frames: filt_guess = q * 0.02295 + 2.48225
        //
        // And high bit depth separately:
        // filt_guess = q * 0.316206 + 3.87252
        int32_t filt_guess;
        switch (scs_ptr->static_config.encoder_bit_depth) {
        case EB_8BIT:
            filt_guess = (frm_hdr->frame_type == KEY_FRAME)
                ? ROUND_POWER_OF_TWO(q * 17563 - 421574, 18)
                : ROUND_POWER_OF_TWO(q * 6017 + 650707, 18);
            break;
        case EB_10BIT:
            filt_guess = ROUND_POWER_OF_TWO(q * 20723 + 4060632, 20);
            break;
        case EB_12BIT:
            filt_guess = ROUND_POWER_OF_TWO(q * 20723 + 16242526, 22);
            break;
        default:
            assert(0 &&
                "bit_depth should be AOM_BITS_8, AOM_BITS_10 "
                "or AOM_BITS_12");
            return;
        }
        if (scs_ptr->static_config.encoder_bit_depth != EB_8BIT && frm_hdr->frame_type == KEY_FRAME)
            filt_guess -= 4;

        filt_guess = filt_guess > 2 ? filt_guess - 2 : filt_guess > 1 ? filt_guess - 1 : filt_guess;
        int32_t filt_guess_chroma = filt_guess > 1 ? filt_guess / 2 : filt_guess;

        // TODO(chengchen): retrain the model for Y, U, V filter levels
        lf->filter_level[0] = clamp(filt_guess, min_filter_level, max_filter_level);
        lf->filter_level[1] = clamp(filt_guess, min_filter_level, max_filter_level);
        lf->filter_level_u = clamp(filt_guess_chroma, min_filter_level, max_filter_level);
        lf->filter_level_v = clamp(filt_guess_chroma, min_filter_level, max_filter_level);
    }
    else {
        const int32_t last_frame_filter_level[4] = { lf->filter_level[0],
            lf->filter_level[1],
            lf->filter_level_u,
            lf->filter_level_v };
        EbPictureBufferDesc  *tempLfReconBuffer = (scs_ptr->static_config.encoder_bit_depth != EB_8BIT) ? context_ptr->temp_lf_recon_picture16bit_ptr : context_ptr->temp_lf_recon_picture_ptr;

        lf->filter_level[0] = lf->filter_level[1] =
            search_filter_level(srcBuffer, tempLfReconBuffer, pcs_ptr, method == LPF_PICK_FROM_SUBIMAGE,
                last_frame_filter_level, NULL, 0, 2);

        if (num_planes > 1) {
            lf->filter_level_u =
                search_filter_level(srcBuffer, tempLfReconBuffer, pcs_ptr, method == LPF_PICK_FROM_SUBIMAGE,
                    last_frame_filter_level, NULL, 1, 0);
            lf->filter_level_v =
                search_filter_level(srcBuffer, tempLfReconBuffer, pcs_ptr, method == LPF_PICK_FROM_SUBIMAGE,
                    last_frame_filter_level, NULL, 2, 0);
        }
    }
}
// clang-format on
