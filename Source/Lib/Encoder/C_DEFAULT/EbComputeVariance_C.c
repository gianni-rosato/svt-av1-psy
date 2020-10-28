/*
* Copyright(c) 2019 Intel Corporation
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#include <stdint.h>

static void variance_c(const uint8_t *a, int a_stride, const uint8_t *b, int b_stride, int w, int h,
                       uint32_t *sse, int *sum) {
    int i, j;

    *sum = 0;
    *sse = 0;

    for (i = 0; i < h; ++i) {
        for (j = 0; j < w; ++j) {
            const int diff = a[j] - b[j];
            *sum += diff;
            *sse += diff * diff;
        }

        a += a_stride;
        b += b_stride;
    }
}

// TODO: use or implement a simd version of this
uint32_t variance_highbd_c(const uint16_t *a, int a_stride, const uint16_t *b, int b_stride, int w,
                           int h, uint32_t *sse) {
    int i, j;

    int sad = 0;
    *sse    = 0;

    for (i = 0; i < h; ++i) {
        for (j = 0; j < w; ++j) {
            const int diff = a[j] - b[j];
            sad += diff;
            *sse += diff * diff;
        }

        a += a_stride;
        b += b_stride;
    }

    return *sse - (sad * sad) / (w * h);
}

#define VAR(W, H)                                                                        \
    uint32_t svt_aom_variance##W##x##H##_c(                                              \
        const uint8_t *a, int a_stride, const uint8_t *b, int b_stride, uint32_t *sse) { \
        int sum;                                                                         \
        variance_c(a, a_stride, b, b_stride, W, H, sse, &sum);                           \
        return *sse - (uint32_t)(((int64_t)sum * sum) / (W * H));                        \
    }

VAR(4, 4)
VAR(4, 8)
VAR(4, 16)
VAR(8, 4)
VAR(8, 8)
VAR(8, 16)
VAR(8, 32)
VAR(16, 4)
VAR(16, 8)
VAR(16, 16)
VAR(16, 32)
VAR(16, 64)
VAR(32, 8)
VAR(32, 16)
VAR(32, 32)
VAR(32, 64)
VAR(64, 16)
VAR(64, 32)
VAR(64, 64)
VAR(64, 128)
VAR(128, 64)
VAR(128, 128)
