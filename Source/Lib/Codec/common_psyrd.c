// clang-format off
/*
* Copyright(c) 2019 Intel Corporation
* Copyright (c) 2016, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#include "definitions.h"
#include "common_psyrd.h"

#define HADAMARD4(d0, d1, d2, d3, s0, s1, s2, s3) { \
        uint32_t t0 = s0 + s1; \
        uint32_t t1 = s0 - s1; \
        uint32_t t2 = s2 + s3; \
        uint32_t t3 = s2 - s3; \
        d0 = t0 + t2; \
        d2 = t0 - t2; \
        d1 = t1 + t3; \
        d3 = t1 - t3; \
}

static int64_t satd_4x4_c(const TranLow *src, const TranLow *zero_buf);
static int64_t sa8d_8x8_c(const TranLow *src, const TranLow *zero_buf);
static int64_t sad_block_c(const TranLow *src, const TranLow *zero_buf, int size);

// Helper functions for AC energy calculation
static int64_t satd_4x4_c(const TranLow *src, const TranLow *zero_buf) {
    uint32_t temp[16], diff[16];
    int i;

    // Get differences
    for (i = 0; i < 16; i++) {
        diff[i] = src[i] - zero_buf[i];
    }

    // Horizontal transform
    for (i = 0; i < 16; i += 4) {
        HADAMARD4(temp[i], temp[i+1], temp[i+2], temp[i+3],
                 diff[i], diff[i+1], diff[i+2], diff[i+3]);
    }

    // Vertical transform
    for (i = 0; i < 4; i++) {
        HADAMARD4(diff[i], diff[i+4], diff[i+8], diff[i+12],
                 temp[i], temp[i+4], temp[i+8], temp[i+12]);
    }

    // Sum absolute values
    int64_t satd = 0;
    for (i = 0; i < 16; i++) {
        satd += diff[i];
    }

    return (satd + 2) >> 2;
}

static int64_t sa8d_8x8_c(const TranLow *src, const TranLow *zero_buf) {
    uint32_t temp[64], diff[64];
    int i, j;

    // Get differences
    for (i = 0; i < 64; i++) {
        diff[i] = src[i] - zero_buf[i];
    }

    // Transform rows - First stage
    for (i = 0; i < 64; i += 8) {
        HADAMARD4(temp[i], temp[i+1], temp[i+2], temp[i+3],
                 diff[i], diff[i+1], diff[i+2], diff[i+3]);
        HADAMARD4(temp[i+4], temp[i+5], temp[i+6], temp[i+7],
                 diff[i+4], diff[i+5], diff[i+6], diff[i+7]);
    }

    // Transform columns - First stage
    for (i = 0; i < 8; i++) {
        for (j = 0; j < 8; j += 4) {
            HADAMARD4(diff[i+j*8], diff[i+(j+1)*8], diff[i+(j+2)*8], diff[i+(j+3)*8],
                     temp[i+j*8], temp[i+(j+1)*8], temp[i+(j+2)*8], temp[i+(j+3)*8]);
        }
    }

    // Sum
    int sa8d = 0;
    for (i = 0; i < 64; i++) {
        sa8d += diff[i];
    }

    return (int64_t)((sa8d + 2) >> 2);
}

static int64_t sad_block_c(const TranLow *src, const TranLow *zero_buf, int size) {
    int64_t sad = 0;
    for (int i = 0; i < size * size; i++) {
        sad += abs(src[i] - zero_buf[i]);
    }
    return sad;
}

// Main psy-rd cost calculation function
int64_t calculate_psy_cost_c(const TranLow *orig_coeff, const TranLow *recon_coeff,
                                int width, int height) {
    static TranLow zero_buf[64] = {0}; // For 8x8 blocks
    int64_t psy_cost = 0;

    if (width >= 8 && height >= 8) {
        // Use SA8D for blocks 8x8 and larger
        for (int i = 0; i < height; i += 8) {
            for (int j = 0; j < width; j += 8) {
                // Calculate AC energy using SA8D - SAD method
                int source_energy = sa8d_8x8_c(&orig_coeff[i * width + j], zero_buf) -
                                  (sad_block_c(&orig_coeff[i * width + j], zero_buf, 8) >> 2);
                int recon_energy = sa8d_8x8_c(&recon_coeff[i * width + j], zero_buf) -
                                 (sad_block_c(&recon_coeff[i * width + j], zero_buf, 8) >> 2);
                psy_cost += abs(source_energy - recon_energy);
            }
        }
    } else if (width == 4 && height == 4) {
        // Use SATD for 4x4 blocks
        int source_energy = satd_4x4_c(orig_coeff, zero_buf) -
                           (sad_block_c(orig_coeff, zero_buf, 4) >> 2);
        int recon_energy = satd_4x4_c(recon_coeff, zero_buf) -
                          (sad_block_c(recon_coeff, zero_buf, 4) >> 2);
        psy_cost += abs(source_energy - recon_energy);
    } else {
        return 0;
    }

    return psy_cost;
}
