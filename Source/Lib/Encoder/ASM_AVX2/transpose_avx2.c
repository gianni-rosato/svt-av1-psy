/*
 * Copyright(c) 2019 Intel Corporation
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

// C instance for unit test.

#include "transpose_avx2.h"

void transpose_8bit_16x16_reg128bit_instance_avx2(const __m128i *const in, __m128i *const out) {
    transpose_8bit_16x16_reg128bit_avx2(in, out);
}

void transpose_32bit_8x8_reg256bit_instance_avx2(const __m256i *const in, __m256i *const out) {
    transpose_32bit_8x8_avx2(in, out);
}

void transpose_64bit_4x4_reg256bit_instance_avx2(const __m256i *const in, __m256i *const out) {
    transpose_64bit_4x4_avx2(in, out);
}

void transpose_64bit_4x6_reg256bit_instance_avx2(const __m256i *const in, __m256i *const out) {
    transpose_64bit_4x6_avx2(in, out);
}

void transpose_64bit_4x8_reg256bit_instance_avx2(const __m256i *const in, __m256i *const out) {
    transpose_64bit_4x8_avx2(in, out);
}
