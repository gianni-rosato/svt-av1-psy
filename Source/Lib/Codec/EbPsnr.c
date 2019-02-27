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

#include "EbDefinitions.h"
#include <assert.h>
#include <math.h>
#include "EbPsnr.h"
#include "EbPictureBufferDesc.h"
#include "aom_dsp_rtcd.h"

double aom_sse_to_psnr(double samples, double peak, double sse) {
    if (sse > 0.0) {
        const double psnr = 10.0 * log10(samples * peak * peak / sse);
        return psnr > MAX_PSNR ? MAX_PSNR : psnr;
    }
    else {
        return MAX_PSNR;
    }
}

/* TODO(yaowu): The block_variance calls the unoptimized versions of variance()
 * and highbd_8_variance(). It should not.
 */
static int32_t encoder_variance(const uint8_t *a, int32_t a_stride, const uint8_t *b,
    int32_t b_stride, int32_t w, int32_t h) {
    int32_t i, j;
    int32_t sse = 0;

    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j++) {
            const int32_t diff = a[j] - b[j];
            sse += diff * diff;
        }

        a += a_stride;
        b += b_stride;
    }
    return sse;
}

static int64_t encoder_highbd_variance64(const uint8_t *a8, int32_t a_stride,
    const uint8_t *b8, int32_t b_stride, int32_t w, int32_t h) {
    const uint16_t *a = CONVERT_TO_SHORTPTR(a8);
    const uint16_t *b = CONVERT_TO_SHORTPTR(b8);
    int64_t sse = 0;
    for (int32_t i = 0; i < h; ++i) {
        for (int32_t j = 0; j < w; ++j) {
            const int32_t diff = a[j] - b[j];
            sse += (uint32_t)(diff * diff);
        }
        a += a_stride;
        b += b_stride;
    }
    return sse;
}

static void encoder_highbd_8_variance(const uint8_t *a8, int32_t a_stride,
    const uint8_t *b8, int32_t b_stride, int32_t w,
    int32_t h, uint32_t *sse) {
    *sse = (uint32_t)encoder_highbd_variance64(a8, a_stride, b8, b_stride, w, h);
}

static void variance(const uint8_t *a, int32_t a_stride, const uint8_t *b,
    int b_stride, int w, int h, uint32_t *sse, int32_t *sum) {
    int i, j;

    *sum = 0;
    *sse = 0;

    for (i = 0; i < h; ++i) {
        for (j = 0; j < w; ++j) {
            const int32_t diff = a[j] - b[j];
            *sum += diff;
            *sse += diff * diff;
        }

        a += a_stride;
        b += b_stride;
    }
}

uint32_t aom_mse16x16_c(const uint8_t *src_ptr, int32_t  source_stride,
    const uint8_t *ref_ptr, int32_t  recon_stride, uint32_t *sse){
    int32_t sum;
    variance(src_ptr, source_stride, ref_ptr, recon_stride,16, 16, sse, &sum);
    return *sse - (uint32_t)(((int64_t)sum * sum) / (16 * 16));
}


static int64_t get_sse(const uint8_t *a, int32_t a_stride, const uint8_t *b,
    int32_t b_stride, int32_t width, int32_t height) {
    const int32_t dw = width % 16;
    const int32_t dh = height % 16;
    int64_t total_sse = 0;
    uint32_t sse = 0;
    int32_t x, y;

    if (dw > 0) {
        sse = encoder_variance(&a[width - dw], a_stride, &b[width - dw], b_stride, dw,
            height);
        total_sse += sse;
    }

    if (dh > 0) {
        sse = encoder_variance(&a[(height - dh) * a_stride], a_stride,
            &b[(height - dh) * b_stride], b_stride, width - dw, dh);
        total_sse += sse;
    }

    for (y = 0; y < height / 16; ++y) {
        const uint8_t *pa = a;
        const uint8_t *pb = b;
        for (x = 0; x < width / 16; ++x) {

            aom_mse16x16(pa, a_stride, pb, b_stride, &sse);

            total_sse += sse;

            pa += 16;
            pb += 16;
        }

        a += 16 * a_stride;
        b += 16 * b_stride;
    }

    return total_sse;
}

static int64_t highbd_get_sse(const uint8_t *a, int32_t a_stride, const uint8_t *b,
    int32_t b_stride, int32_t width, int32_t height) {
    int64_t total_sse = 0;
    int32_t x, y;
    const int32_t dw = width % 16;
    const int32_t dh = height % 16;
    uint32_t sse = 0;
    if (dw > 0) {
        encoder_highbd_8_variance(&a[width - dw], a_stride, &b[width - dw],
            b_stride, dw, height, &sse);
        total_sse += sse;
    }
    if (dh > 0) {
        encoder_highbd_8_variance(&a[(height - dh) * a_stride], a_stride,
            &b[(height - dh) * b_stride], b_stride,
            width - dw, dh, &sse);
        total_sse += sse;
    }
    for (y = 0; y < height / 16; ++y) {
        const uint8_t *pa = a;
        const uint8_t *pb = b;
        for (x = 0; x < width / 16; ++x) {
#if RS_10BIT_FIX
            aom_highbd_8_mse16x16(pa, a_stride, pb, b_stride, &sse);
#else
            //aom_highbd_8_mse16x16(pa, a_stride, pb, b_stride, &sse);
#endif
            total_sse += sse;
            pa += 16;
            pb += 16;
        }
        a += 16 * a_stride;
        b += 16 * b_stride;
    }
    return total_sse;
}

int64_t aom_get_y_sse_part(const Yv12BufferConfig *a,
    const Yv12BufferConfig *b, int32_t hstart, int32_t width,
    int32_t vstart, int32_t height) {
    return get_sse(a->y_buffer + vstart * a->y_stride + hstart, a->y_stride,
        b->y_buffer + vstart * b->y_stride + hstart, b->y_stride,
        width, height);
}

int64_t aom_get_y_sse(const Yv12BufferConfig *a,
    const Yv12BufferConfig *b) {
    assert(a->y_crop_width == b->y_crop_width);
    assert(a->y_crop_height == b->y_crop_height);

    return get_sse(a->y_buffer, a->y_stride, b->y_buffer, b->y_stride,
        a->y_crop_width, a->y_crop_height);
}

int64_t aom_get_u_sse_part(const Yv12BufferConfig *a,
    const Yv12BufferConfig *b, int32_t hstart, int32_t width,
    int32_t vstart, int32_t height) {
    return get_sse(a->u_buffer + vstart * a->uv_stride + hstart, a->uv_stride,
        b->u_buffer + vstart * b->uv_stride + hstart, b->uv_stride,
        width, height);
}

int64_t aom_get_u_sse(const Yv12BufferConfig *a,
    const Yv12BufferConfig *b) {
    assert(a->uv_crop_width == b->uv_crop_width);
    assert(a->uv_crop_height == b->uv_crop_height);

    return get_sse(a->u_buffer, a->uv_stride, b->u_buffer, b->uv_stride,
        a->uv_crop_width, a->uv_crop_height);
}

int64_t aom_get_v_sse_part(const Yv12BufferConfig *a,
    const Yv12BufferConfig *b, int32_t hstart, int32_t width,
    int32_t vstart, int32_t height) {
    return get_sse(a->v_buffer + vstart * a->uv_stride + hstart, a->uv_stride,
        b->v_buffer + vstart * b->uv_stride + hstart, b->uv_stride,
        width, height);
}

int64_t aom_get_v_sse(const Yv12BufferConfig *a,
    const Yv12BufferConfig *b) {
    assert(a->uv_crop_width == b->uv_crop_width);
    assert(a->uv_crop_height == b->uv_crop_height);

    return get_sse(a->v_buffer, a->uv_stride, b->v_buffer, b->uv_stride,
        a->uv_crop_width, a->uv_crop_height);
}

int64_t aom_highbd_get_y_sse_part(const Yv12BufferConfig *a,
    const Yv12BufferConfig *b, int32_t hstart,
    int32_t width, int32_t vstart, int32_t height) {
    return highbd_get_sse(
        a->y_buffer + vstart * a->y_stride + hstart, a->y_stride,
        b->y_buffer + vstart * b->y_stride + hstart, b->y_stride, width, height);
}

int64_t aom_highbd_get_y_sse(const Yv12BufferConfig *a,
    const Yv12BufferConfig *b) {
    assert(a->y_crop_width == b->y_crop_width);
    assert(a->y_crop_height == b->y_crop_height);
    assert((a->flags & YV12_FLAG_HIGHBITDEPTH) != 0);
    assert((b->flags & YV12_FLAG_HIGHBITDEPTH) != 0);

    return highbd_get_sse(a->y_buffer, a->y_stride, b->y_buffer, b->y_stride,
        a->y_crop_width, a->y_crop_height);
}

int64_t aom_highbd_get_u_sse_part(const Yv12BufferConfig *a,
    const Yv12BufferConfig *b, int32_t hstart,
    int32_t width, int32_t vstart, int32_t height) {
    return highbd_get_sse(a->u_buffer + vstart * a->uv_stride + hstart,
        a->uv_stride,
        b->u_buffer + vstart * b->uv_stride + hstart,
        b->uv_stride, width, height);
}

int64_t aom_highbd_get_u_sse(const Yv12BufferConfig *a,
    const Yv12BufferConfig *b) {
    assert(a->uv_crop_width == b->uv_crop_width);
    assert(a->uv_crop_height == b->uv_crop_height);
    assert((a->flags & YV12_FLAG_HIGHBITDEPTH) != 0);
    assert((b->flags & YV12_FLAG_HIGHBITDEPTH) != 0);

    return highbd_get_sse(a->u_buffer, a->uv_stride, b->u_buffer, b->uv_stride,
        a->uv_crop_width, a->uv_crop_height);
}

int64_t aom_highbd_get_v_sse_part(const Yv12BufferConfig *a,
    const Yv12BufferConfig *b, int32_t hstart,
    int32_t width, int32_t vstart, int32_t height) {
    return highbd_get_sse(a->v_buffer + vstart * a->uv_stride + hstart,
        a->uv_stride,
        b->v_buffer + vstart * b->uv_stride + hstart,
        b->uv_stride, width, height);
}

int64_t aom_highbd_get_v_sse(const Yv12BufferConfig *a,
    const Yv12BufferConfig *b) {
    assert(a->uv_crop_width == b->uv_crop_width);
    assert(a->uv_crop_height == b->uv_crop_height);
    assert((a->flags & YV12_FLAG_HIGHBITDEPTH) != 0);
    assert((b->flags & YV12_FLAG_HIGHBITDEPTH) != 0);

    return highbd_get_sse(a->v_buffer, a->uv_stride, b->v_buffer, b->uv_stride,
        a->uv_crop_width, a->uv_crop_height);
}




