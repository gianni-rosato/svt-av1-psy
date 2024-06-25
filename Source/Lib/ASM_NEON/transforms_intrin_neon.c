/*
 * Copyright(c) 2024 Intel Corporation
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
 */

#include "definitions.h"
#include "motion_estimation.h"

#include <arm_neon.h>

static INLINE void energy_computation_kernel_neon(const int32_t *const in, uint64x2_t *const sum256) {
    const int32x4_t  input     = vld1q_s32(in);
    const int32x2_t  in_lo     = vget_low_s32(input);
    const int32x2_t  in_hi     = vget_high_s32(input);
    const uint64x2_t energy_lo = vreinterpretq_u64_s64(vmull_s32(in_lo, in_lo));
    const uint64x2_t energy_hi = vreinterpretq_u64_s64(vmull_s32(in_hi, in_hi));
    *sum256                    = vaddq_u64(*sum256, vaddq_u64(energy_lo, energy_hi));
}

static INLINE uint64_t hadd64_neon(const uint64x2_t sum256) {
    const uint64x2_t partial_sum = vpaddq_u64(sum256, sum256);
    return vgetq_lane_u64(partial_sum, 0);
}

static INLINE uint64_t energy_computation_neon(const int32_t *const in, const uint32_t size) {
    uint64x2_t sum = vdupq_n_u64(0);
    uint32_t   i   = 0;

    do {
        energy_computation_kernel_neon(in + i + 0, &sum);
        energy_computation_kernel_neon(in + i + 4, &sum);
        i += 8;
    } while (i < size);

    return hadd64_neon(sum);
}

static INLINE uint64_t energy_computation_64_neon(const int32_t *in, const uint32_t height) {
    uint64x2_t sum = vdupq_n_u64(0);
    uint32_t   i   = height;

    do {
        energy_computation_kernel_neon(in + 0 * 8 + 0, &sum);
        energy_computation_kernel_neon(in + 0 * 8 + 4, &sum);
        energy_computation_kernel_neon(in + 1 * 8 + 0, &sum);
        energy_computation_kernel_neon(in + 1 * 8 + 4, &sum);
        energy_computation_kernel_neon(in + 2 * 8 + 0, &sum);
        energy_computation_kernel_neon(in + 2 * 8 + 4, &sum);
        energy_computation_kernel_neon(in + 3 * 8 + 0, &sum);
        energy_computation_kernel_neon(in + 3 * 8 + 4, &sum);
        in += 64;
    } while (--i);

    return hadd64_neon(sum);
}

static INLINE void copy_32_bytes_neon(const int32_t *src, int32_t *dst) {
    const int32x4x2_t val = vld2q_s32(src + 0 * 8);
    vst2q_s32(dst + 0 * 8, val);
}

static INLINE void copy_256x_bytes_neon(const int32_t *src, int32_t *dst, const uint32_t height) {
    uint32_t h = height;

    do {
        copy_32_bytes_neon(src + 0 * 8, dst + 0 * 8);
        copy_32_bytes_neon(src + 1 * 8, dst + 1 * 8);
        copy_32_bytes_neon(src + 2 * 8, dst + 2 * 8);
        copy_32_bytes_neon(src + 3 * 8, dst + 3 * 8);
        src += 64;
        dst += 32;
    } while (--h);
}

uint64_t svt_handle_transform16x64_neon(int32_t *output) {
    //bottom 16x32 area.
    const uint64_t three_quad_energy = energy_computation_neon(output + 16 * 32, 16 * 32);
    return three_quad_energy;
}

uint64_t svt_handle_transform32x64_neon(int32_t *output) {
    //bottom 32x32 area.
    const uint64_t three_quad_energy = energy_computation_neon(output + 32 * 32, 32 * 32);
    return three_quad_energy;
}

uint64_t svt_handle_transform64x16_neon(int32_t *output) {
    // top - right 32x16 area.
    const uint64_t three_quad_energy = energy_computation_64_neon(output + 32, 16);
    // Re-pack non-zero coeffs in the first 32x16 indices.
    copy_256x_bytes_neon(output + 64, output + 32, 15);

    return three_quad_energy;
}

uint64_t svt_handle_transform64x32_neon(int32_t *output) {
    // top - right 32x32 area.
    const uint64_t three_quad_energy = energy_computation_64_neon(output + 32, 32);
    // Re-pack non-zero coeffs in the first 32x32 indices.
    copy_256x_bytes_neon(output + 64, output + 32, 31);

    return three_quad_energy;
}

uint64_t svt_handle_transform64x64_neon(int32_t *output) {
    uint64_t three_quad_energy;

    // top - right 32x32 area.
    three_quad_energy = energy_computation_64_neon(output + 32, 32);
    //bottom 64x32 area.
    three_quad_energy += energy_computation_neon(output + 32 * 64, 64 * 32);
    // Re-pack non-zero coeffs in the first 32x32 indices.
    copy_256x_bytes_neon(output + 64, output + 32, 31);

    return three_quad_energy;
}

uint64_t svt_handle_transform16x64_N2_N4_neon(int32_t *output) {
    (void)output;
    return 0;
}

uint64_t svt_handle_transform32x64_N2_N4_neon(int32_t *output) {
    (void)output;
    return 0;
}

uint64_t svt_handle_transform64x16_N2_N4_neon(int32_t *output) {
    // Re-pack non-zero coeffs in the first 32x16 indices.
    copy_256x_bytes_neon(output + 64, output + 32, 15);
    return 0;
}

uint64_t svt_handle_transform64x32_N2_N4_neon(int32_t *output) {
    // Re-pack non-zero coeffs in the first 32x32 indices.
    copy_256x_bytes_neon(output + 64, output + 32, 31);
    return 0;
}

uint64_t svt_handle_transform64x64_N2_N4_neon(int32_t *output) {
    // Re-pack non-zero coeffs in the first 32x32 indices.
    copy_256x_bytes_neon(output + 64, output + 32, 31);
    return 0;
}
