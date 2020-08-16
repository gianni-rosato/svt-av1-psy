/*
* Copyright(c) 2019 Intel Corporation
* Copyright (c) 2019, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
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
