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

#include "EbComputeSAD_C.h"
#include "EbUtility.h"
/*******************************************
*   returns NxM Sum of Absolute Differences
Note: moved from picture operators.
keep this function here for profiling
issues.
*******************************************/
uint32_t svt_fast_loop_nxm_sad_kernel(const uint8_t *src, // input parameter, source samples Ptr
                                      uint32_t       src_stride, // input parameter, source stride
                                      const uint8_t *ref, // input parameter, reference samples Ptr
                                      uint32_t       ref_stride, // input parameter, reference stride
                                      uint32_t       height, // input parameter, block height (M)
                                      uint32_t       width) // input parameter, block width (N)
{
    uint32_t x, y;
    uint32_t sad = 0;

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) sad += EB_ABS_DIFF(src[x], ref[x]);
        src += src_stride;
        ref += ref_stride;
    }

    return sad;
}

uint32_t sad_16b_kernel_c(uint16_t *src, // input parameter, source samples Ptr
                          uint32_t  src_stride, // input parameter, source stride
                          uint16_t *ref, // input parameter, reference samples Ptr
                          uint32_t  ref_stride, // input parameter, reference stride
                          uint32_t  height, // input parameter, block height (M)
                          uint32_t  width) // input parameter, block width (N)
{
    uint32_t x, y;
    uint32_t sad = 0;

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) sad += EB_ABS_DIFF(src[x], ref[x]);
        src += src_stride;
        ref += ref_stride;
    }

    return sad;
}

void svt_sad_loop_kernel_c(uint8_t * src, // input parameter, source samples Ptr
                           uint32_t  src_stride, // input parameter, source stride
                           uint8_t * ref, // input parameter, reference samples Ptr
                           uint32_t  ref_stride, // input parameter, reference stride
                           uint32_t  block_height, // input parameter, block height (M)
                           uint32_t  block_width, // input parameter, block width (N)
                           uint64_t *best_sad, int16_t *x_search_center, int16_t *y_search_center,
                           uint32_t src_stride_raw, // input parameter, source stride (no line skipping)
                           int16_t search_area_width, int16_t search_area_height) {
    int16_t x_search_index;
    int16_t y_search_index;

    *best_sad = 0xffffff;

    for (y_search_index = 0; y_search_index < search_area_height; y_search_index++) {
        for (x_search_index = 0; x_search_index < search_area_width; x_search_index++) {
            uint32_t x, y;
            uint32_t sad = 0;

            for (y = 0; y < block_height; y++) {
                for (x = 0; x < block_width; x++)
                    sad += EB_ABS_DIFF(src[y * src_stride + x],
                                       ref[x_search_index + y * ref_stride + x]);
            }

            // Update results
            if (sad < *best_sad) {
                *best_sad        = sad;
                *x_search_center = x_search_index;
                *y_search_center = y_search_index;
            }
        }

        ref += src_stride_raw;
    }

    return;
}

/* Sum the difference between every corresponding element of the buffers. */
static INLINE uint32_t sad_inline_c(const uint8_t *a, int a_stride, const uint8_t *b, int b_stride,
                                    int width, int height) {
    int          y, x;
    unsigned int sad = 0;

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) sad += EB_ABS_DIFF(a[x], b[x]);
        a += a_stride;
        b += b_stride;
    }
    return sad;
}

#define sadMxN(m, n)                                                              \
    uint32_t eb_aom_sad##m##x##n##_c(                                             \
        const uint8_t *src, int src_stride, const uint8_t *ref, int ref_stride) { \
        return sad_inline_c(src, src_stride, ref, ref_stride, m, n);              \
    }

// Calculate sad against 4 reference locations and store each in sad_array
#define sadMxNx4D(m, n)                                                                        \
    void eb_aom_sad##m##x##n##x4d_c(const uint8_t *      src,                                  \
                                    int                  src_stride,                           \
                                    const uint8_t *const ref_array[],                          \
                                    int                  ref_stride,                           \
                                    uint32_t *           sad_array) {                                     \
        int i;                                                                                 \
        for (i = 0; i < 4; ++i) {                                                              \
            sad_array[i] = eb_aom_sad##m##x##n##_c(src, src_stride, ref_array[i], ref_stride); \
        }                                                                                      \
    }

// 128x128
sadMxN(128, 128);
sadMxNx4D(128, 128);

// 128x64
sadMxN(128, 64);
sadMxNx4D(128, 64);

// 64x128
sadMxN(64, 128);
sadMxNx4D(64, 128);

// 64x64
sadMxN(64, 64);
sadMxNx4D(64, 64);

// 64x32
sadMxN(64, 32);
sadMxNx4D(64, 32);

// 32x64
sadMxN(32, 64);
sadMxNx4D(32, 64);

// 32x32
sadMxN(32, 32);
sadMxNx4D(32, 32);

// 32x16
sadMxN(32, 16);
sadMxNx4D(32, 16);

// 16x32
sadMxN(16, 32);
sadMxNx4D(16, 32);

// 16x16
sadMxN(16, 16);
sadMxNx4D(16, 16);

// 16x8
sadMxN(16, 8);
sadMxNx4D(16, 8);

// 8x16
sadMxN(8, 16);
sadMxNx4D(8, 16);

// 8x8
sadMxN(8, 8);
sadMxNx4D(8, 8);

// 8x4
sadMxN(8, 4);
sadMxNx4D(8, 4);

// 4x8
sadMxN(4, 8);
sadMxNx4D(4, 8);

// 4x4
sadMxN(4, 4);
sadMxNx4D(4, 4);

sadMxN(4, 16);
sadMxNx4D(4, 16);
sadMxN(16, 4);
sadMxNx4D(16, 4);
sadMxN(8, 32);
sadMxNx4D(8, 32);
sadMxN(32, 8);
sadMxNx4D(32, 8);
sadMxN(16, 64);
sadMxNx4D(16, 64);
sadMxN(64, 16);
sadMxNx4D(64, 16);

uint32_t svt_nxm_sad_kernel_helper_c(const uint8_t *src, uint32_t src_stride, const uint8_t *ref,
                                 uint32_t ref_stride, uint32_t height, uint32_t width) {

    return svt_fast_loop_nxm_sad_kernel(src, src_stride, ref, ref_stride, height, width);
};
