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

#ifndef AOM_AV1_COMMON_ARM_CONVOLVE_NEON_H_
#define AOM_AV1_COMMON_ARM_CONVOLVE_NEON_H_

#include "inter_prediction.h"

static INLINE Bool is_convolve_2tap(const int16_t *const filter) {
    return (const void *)filter == (const void *)bilinear_filters;
}

static INLINE Bool is_convolve_4tap(const int16_t *const filter) {
    return (const void *)filter == (const void *)sub_pel_filters_4 ||
        (const void *)filter == (const void *)sub_pel_filters_4smooth;
}

static INLINE Bool is_convolve_6tap(const int16_t *const filter) {
    return (const void *)filter == (const void *)sub_pel_filters_8 ||
        (const void *)filter == (const void *)sub_pel_filters_8smooth;
}

static INLINE int32_t get_convolve_tap(const int16_t *const filter) {
    if (is_convolve_2tap(filter))
        return 2;
    else if (is_convolve_4tap(filter))
        return 4;
    else if (is_convolve_6tap(filter))
        return 6;
    else
        return 8;
}

#endif // AOM_AV1_COMMON_ARM_CONVOLVE_NEON_H_
