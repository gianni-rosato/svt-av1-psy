/*
 * Copyright(c) 2019 Intel Corporation
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

// C instance for unit test.

#include "EbDefinitions.h"

#include "transpose_sse2.h"

void transpose_8bit_4x4_reg128bit_instance_sse2(const __m128i *const in,
                                                __m128i *const out) {
    __m128i out_tmp[4] = {0};
    out_tmp[0] = transpose_8bit_4x4(in);
    for (size_t i = 0; i < 4; i++) {
        out[i] = _mm_loadu_si128((__m128i *)(&((uint8_t *)(&out_tmp))[4 * i]));
    }
}

void transpose_8bit_8x8_reg128bit_instance_sse2(const __m128i *const in,
                                                __m128i *const out) {
    transpose_8bit_8x8(in, out);
}

void transpose_8bit_16x8_reg128bit_instance_sse2(const __m128i *const in,
                                                 __m128i *const out) {
    transpose_8bit_16x8(in, out);
}

void transpose_8bit_16x16_reg128bit_instance_sse2(const __m128i *const in,
                                                  __m128i *const out) {
    transpose_8bit_16x16_sse2(in, out);
}

void partial_transpose_8bit_8x8_reg128bit_instance_sse2(const __m128i *const in,
                                                        __m128i *const out) {
    partial_transpose_8bit_8x8(in, out);
}

void transpose_16bit_4x4_reg128bit_instance_sse2(const __m128i *const in,
                                                 __m128i *const out) {
    transpose_16bit_4x4(in, out);
}

void transpose_16bit_4x8_reg128bit_instance_sse2(const __m128i *const in,
                                                 __m128i *const out) {
    transpose_16bit_4x8(in, out);
}

void transpose_16bit_8x4_reg128bit_instance_sse2(const __m128i *const in,
                                                 __m128i *const out) {
    transpose_16bit_8x4(in, out);
}

void transpose_16bit_8x8_reg128bit_instance_sse2(const __m128i *const in,
                                                 __m128i *const out) {
    transpose_16bit_8x8(in, out);
}

void transpose_16bit_16x16_reg128bit_instance_sse2(const __m128i *const in,
                                                   __m128i *const out) {
    memcpy(out, in, 16 * 16 * sizeof(__m128i));
    transpose_16bit_16x16(out, &out[16]);
}

void transpose_32bit_4x4_reg128bit_instance_sse2(const __m128i *const in,
                                                 __m128i *const out) {
    transpose_32bit_4x4(in, out);
}

void transpose_32bit_4x4x2_reg128bit_instance_sse2(const __m128i *const in,
                                                   __m128i *const out) {
    transpose_32bit_4x4x2(in, out);
}

void transpose_32bit_8x4_reg128bit_instance_sse2(const __m128i *const in,
                                                 __m128i *const out) {
    transpose_32bit_8x4(in, out);
}
