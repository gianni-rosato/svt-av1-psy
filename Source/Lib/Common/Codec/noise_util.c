/*
 * Copyright (c) 2017, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <math.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "EbDefinitions.h"
#include "noise_util.h"
#include "aom_dsp_rtcd.h"

void *aom_memalign(size_t align, size_t size);
void aom_free(void *memblk);

float aom_noise_psd_get_default_value(int32_t block_size, float factor) {
    return (factor * factor / 10000) * block_size * block_size / 8;
}

// Internal representation of noise transform. It keeps track of the
// transformed data and a temporary working buffer to use during the
// transform.
struct aom_noise_tx_t {
    DECLARE_ALIGNED(32, float, *tx_block);
    DECLARE_ALIGNED(32, float, *temp);
    int32_t block_size;
    void(*fft)(const float *, float *, float *);
    void(*ifft)(const float *, float *, float *);
};

struct aom_noise_tx_t *aom_noise_tx_malloc(int32_t block_size) {
    struct aom_noise_tx_t *noise_tx =
        (struct aom_noise_tx_t *)malloc(sizeof(struct aom_noise_tx_t));
    if (!noise_tx) return NULL;
    memset(noise_tx, 0, sizeof(*noise_tx));
    switch (block_size) {
    case 2:
        noise_tx->fft = aom_fft2x2_float;
        noise_tx->ifft = aom_ifft2x2_float;
        break;
    case 4:
        noise_tx->fft = aom_fft4x4_float;
        noise_tx->ifft = aom_ifft4x4_float;
        break;
    case 8:
        noise_tx->fft = aom_fft8x8_float;
        noise_tx->ifft = aom_ifft8x8_float;
        break;
    case 16:
        noise_tx->fft = aom_fft16x16_float;
        noise_tx->ifft = aom_ifft16x16_float;
        break;
    case 32:
        noise_tx->fft = aom_fft32x32_float;
        noise_tx->ifft = aom_ifft32x32_float;
        break;
    default:
        free(noise_tx);
        fprintf(stderr, "Unsupported block size %d\n", block_size);
        return NULL;
    }
    noise_tx->block_size = block_size;
    noise_tx->tx_block = (float *)aom_memalign(
        32, 2 * sizeof(*noise_tx->tx_block) * block_size * block_size);
    noise_tx->temp = (float *)aom_memalign(
        32, 2 * sizeof(*noise_tx->temp) * block_size * block_size);
    if (!noise_tx->tx_block || !noise_tx->temp) {
        aom_noise_tx_free(noise_tx);
        return NULL;
    }
    // Clear the buffers up front. Some outputs of the forward transform are
    // real only (the imaginary component will never be touched)
    memset(noise_tx->tx_block, 0,
        2 * sizeof(*noise_tx->tx_block) * block_size * block_size);
    memset(noise_tx->temp, 0,
        2 * sizeof(*noise_tx->temp) * block_size * block_size);
    return noise_tx;
}

void aom_noise_tx_forward(struct aom_noise_tx_t *noise_tx, const float *data) {
    noise_tx->fft(data, noise_tx->temp, noise_tx->tx_block);
}

void aom_noise_tx_filter(struct aom_noise_tx_t *noise_tx, const float *psd) {
    const int32_t block_size = noise_tx->block_size;
    const float kBeta = 1.1f;
    const float kEps = 1e-6f;
    for (int32_t y = 0; y < block_size; ++y) {
        for (int32_t x = 0; x < block_size; ++x) {
            int32_t i = y * block_size + x;
            float *c = noise_tx->tx_block + 2 * i;
            const float p = c[0] * c[0] + c[1] * c[1];
            if (p > kBeta * psd[i] && p > 1e-6) {
                noise_tx->tx_block[2 * i + 0] *= (p - psd[i]) / AOMMAX(p, kEps);
                noise_tx->tx_block[2 * i + 1] *= (p - psd[i]) / AOMMAX(p, kEps);
            }
            else {
                noise_tx->tx_block[2 * i + 0] *= (kBeta - 1.0f) / kBeta;
                noise_tx->tx_block[2 * i + 1] *= (kBeta - 1.0f) / kBeta;
            }
        }
    }
}

void aom_noise_tx_inverse(struct aom_noise_tx_t *noise_tx, float *data) {
    const int32_t n = noise_tx->block_size * noise_tx->block_size;
    noise_tx->ifft(noise_tx->tx_block, noise_tx->temp, data);
    for (int32_t i = 0; i < n; ++i) {
        data[i] /= n;
    }
}

void aom_noise_tx_free(struct aom_noise_tx_t *noise_tx) {
    if (!noise_tx) return;
    aom_free(noise_tx->tx_block);
    aom_free(noise_tx->temp);
    free(noise_tx);
}

double aom_normalized_cross_correlation(const double *a, const double *b,
    int32_t n) {
    double c = 0;
    double a_len = 0;
    double b_len = 0;
    for (int32_t i = 0; i < n; ++i) {
        a_len += a[i] * a[i];
        b_len += b[i] * b[i];
        c += a[i] * b[i];
    }
    return c / (sqrt(a_len) * sqrt(b_len));
}
