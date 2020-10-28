/*
 * Copyright(c) 2020 Intel Corporation
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at https://www.aomedia.org/license/software-license. If the
 * Alliance for Open Media Patent License 1.0 was not distributed with this
 * source code in the PATENTS file, you can obtain it at
 * https://www.aomedia.org/license/patent-license.
 */

#include <sstream>
#include "gtest/gtest.h"
#include "random.h"
#include "aom_dsp_rtcd.h"

namespace {

using svt_av1_test_tool::SVTRandom;

static const int block_size = 8 * 8;
static const int test_times = 10000;

static const uint16_t* prepare_data_8x8_10b(uint16_t* data, SVTRandom* rnd) {
    for (size_t i = 0; i < block_size; i++) {
        data[i] = (uint16_t)rnd->random();
    }
    return data;
}

static const uint8_t* prepare_data_8x8(uint8_t* data, SVTRandom* rnd) {
    for (size_t i = 0; i < block_size; i++) {
        data[i] = (uint8_t)rnd->random();
    }
    return data;
}
TEST(haar_ac_sad_8x8_uint8_input_test, 8bit) {
    SVTRandom rnd = SVTRandom(8, false);

    uint8_t input_data[block_size];

    for (int i = 0; i < test_times; i++) {
        prepare_data_8x8(input_data, &rnd);

        int output_tst = svt_av1_haar_ac_sad_8x8_uint8_input_avx2(input_data, 8, 0);
        int output_c_ref = svt_av1_haar_ac_sad_8x8_uint8_input_c(input_data, 8, 0);

        // compare results
        ASSERT_EQ(output_tst, output_c_ref)
            << "test"
            << "[" << i << "] "
            << "compute av1_haar_ac_sad_8x8_uint8_input_avx2 failed!\n";
    }
}

TEST(haar_ac_sad_8x8_uint8_input_test, 10bit) {
    SVTRandom rnd = SVTRandom(10, false);

    uint16_t input_data[block_size];

    for (int i = 0; i < test_times; i++) {
        prepare_data_8x8_10b(input_data, &rnd);

        uint8_t* inp_8bit_ptr = CONVERT_TO_BYTEPTR(input_data);

        int output_tst =
            svt_av1_haar_ac_sad_8x8_uint8_input_avx2(inp_8bit_ptr, 8, 1);
        int output_c_ref =
            svt_av1_haar_ac_sad_8x8_uint8_input_c(inp_8bit_ptr, 8, 1);

        // compare results
        ASSERT_EQ(output_tst, output_c_ref)
            << "test"
            << "[" << i << "] "
            << "compute av1_haar_ac_sad_8x8_uint8_input_avx2 failed!\n";
    }
}
}  // namespace
