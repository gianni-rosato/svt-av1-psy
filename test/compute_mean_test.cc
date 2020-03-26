/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file compute_mean_test.cc
 *
 * @brief Unit test for compute mean function:
 * - compute_mean8x8_sse2_intrin
 * - compute_mean_of_squared_values8x8_sse2_intrin
 * - compute_sub_mean8x8_sse2_intrin
 * - compute_subd_mean_of_squared_values8x8_sse2_intrin
 * - compute_mean8x8_avx2_intrin
 * - compute_interm_var_four8x8_avx2_intrin
 *
 * @author Cidana-Edmond,Cidana-Ivy
 *
 ******************************************************************************/

#include <sstream>
#include "gtest/gtest.h"
#include "EbComputeMean.h"
#include "random.h"
#include "aom_dsp_rtcd.h"
/**
 * @brief Unit test for compute mean function:
 * - compute_mean8x8_sse2_intrin
 * - compute_mean_of_squared_values8x8_sse2_intrin
 * - compute_sub_mean8x8_sse2_intrin
 * - compute_subd_mean_of_squared_values8x8_sse2_intrin
 * - compute_mean8x8_avx2_intrin
 *
 * Test strategy:
 * Verify this assembly code by comparing with reference c implementation.
 * Feed the same data and check test output and reference output.
 *
 * Expected result:
 * Output from assemble functions should be the same with output from c.
 *
 * Test coverage:
 * Test cases:
 * data buffer:
 * Buffer is filled with test data. The values of data in normal test are the
 * random 8-bit integer, and in boundary test are the large random integer
 * between 0xE0 and 0xFF
 */

namespace {

using svt_av1_test_tool::SVTRandom;

static const int block_size = 8 * 8;
static const int test_times = 10000;
static const std::string test_name[2] = {"Noraml Test:\n", "Boundary Test:\n"};

static const std::string print_data(const uint8_t* data, const int width,
                                    const int height) {
    std::string print_str;
    std::stringstream ss(print_str);
    ss << "test input data dump:\n";
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            ss << std::to_string(data[j * width + i]) << "\t";
        }
        ss << "\n";
    }
    return ss.str();
}

static const uint8_t* prepare_data_8x8(uint8_t* data, SVTRandom* rnd) {
    for (size_t i = 0; i < block_size; i++) {
        data[i] = (uint8_t)rnd->random();
    }
    return data;
}

TEST(ComputeMeanTest, run_compute_mean_test) {
    SVTRandom rnd[2] = {
        SVTRandom(8, false),  /**< random generator of normal test vector */
        SVTRandom(0xE0, 0xFF) /**< random generator of boundary test vector */
    };
    uint8_t input_data[block_size];

    for (size_t vi = 0; vi < 2; vi++) {
        for (int i = 0; i < test_times; i++) {
            // prepare data
            prepare_data_8x8(input_data, &rnd[vi]);

            // compute mean
            uint64_t output_sse2_tst =
                compute_mean8x8_sse2_intrin(input_data, 8, 8, 8);
            uint64_t output_avx2_tst =
                compute_mean8x8_avx2_intrin(input_data, 8, 8, 8);
            uint64_t output_c_ref = compute_mean_c(input_data, 8, 8, 8);

            // compare results
            ASSERT_EQ(output_sse2_tst, output_c_ref)
                << test_name[vi] << "[" << i << "] "
                << "compute mean with asm SSE2 failed!\n"
                << print_data(input_data, 8, 8);
            ASSERT_EQ(output_avx2_tst, output_c_ref)
                << test_name[vi] << "[" << i << "] "
                << "compute mean with asm AVX2 failed!\n"
                << print_data(input_data, 8, 8);
        }
    }
}

TEST(ComputeMeanTest, run_compute_mean_squared_values_test) {
    SVTRandom rnd[2] = {
        SVTRandom(8, false),  /**< random generator of normal test vector */
        SVTRandom(0xE0, 0xFF) /**< random generator of boundary test vector */
    };
    uint8_t input_data[block_size];

    for (size_t vi = 0; vi < 2; vi++) {
        for (int i = 0; i < test_times; i++) {
            // prepare data
            prepare_data_8x8(input_data, &rnd[vi]);

            // compute mean
            uint64_t output_sse2_tst =
                compute_mean_of_squared_values8x8_sse2_intrin(
                    input_data, 8, 8, 8);
            uint64_t output_c_ref =
                compute_mean_squared_values_c(input_data, 8, 8, 8);

            // compare results
            ASSERT_EQ(output_sse2_tst, output_c_ref)
                << test_name[vi] << "[" << i << "] "
                << "compute mean of squared values with asm SSE2 failed!\n"
                << print_data(input_data, 8, 8);
        }
    }
}

TEST(ComputeMeanTest, run_compute_sub_mean_test) {
    SVTRandom rnd[2] = {
        SVTRandom(8, false),  /**< random generator of normal test vector */
        SVTRandom(0xE0, 0xFF) /**< random generator of boundary test vector */
    };
    uint8_t input_data[block_size];

    for (size_t vi = 0; vi < 2; vi++) {
        for (int i = 0; i < test_times; i++) {
            // prepare data
            prepare_data_8x8(input_data, &rnd[vi]);

            // compute mean
            uint64_t output_sse2_tst =
                compute_sub_mean8x8_sse2_intrin(input_data, 8);
            uint64_t output_c_ref = compute_sub_mean_8x8_c(input_data, 8);

            // compare results
            ASSERT_EQ(output_sse2_tst, output_c_ref)
                << test_name[vi] << "[" << i << "] "
                << "compute sub mean with asm SSE2 failed!\n"
                << print_data(input_data, 8, 8);
        }
    }
}

TEST(ComputeMeanTest, run_compute_sub_mean_squared_values_test) {
    SVTRandom rnd[2] = {
        SVTRandom(8, false),  /**< random generator of normal test vector */
        SVTRandom(0xE0, 0xFF) /**< random generator of boundary test vector */
    };
    uint8_t input_data[block_size];

    for (size_t vi = 0; vi < 2; vi++) {
        for (int i = 0; i < test_times; i++) {
            // prepare data
            prepare_data_8x8(input_data, &rnd[vi]);

            // compute mean
            uint64_t output_sse2_tst =
                compute_subd_mean_of_squared_values8x8_sse2_intrin(input_data,
                                                                   8);
            uint64_t output_c_ref =
                compute_sub_mean_squared_values_c(input_data, 8, 8, 8);

            // compare results
            ASSERT_EQ(output_sse2_tst, output_c_ref)
                << test_name[vi] << "[" << i << "] "
                << "compute sub mean of squared values with asm SSE2 failed!\n"
                << print_data(input_data, 8, 8);
        }
    }
}

TEST(ComputeMeanTest, run_compute_mean_avx2_test) {
    SVTRandom rnd[2] = {
        SVTRandom(8, false),  /**< random generator of normal test vector */
        SVTRandom(0xE0, 0xFF) /**< random generator of boundary test vector */
    };

    uint8_t input_data[block_size];

    for (int i = 0; i < 2; i++) {
        // prepare data
        prepare_data_8x8(input_data, &rnd[i]);

        // compute mean
        uint64_t output_sse2_squared_tst =
            compute_subd_mean_of_squared_values8x8_sse2_intrin(input_data, 8);
        uint64_t output_sse2_tst =
            compute_sub_mean8x8_sse2_intrin(input_data, 8);

        uint64_t output_avx2_tst[4] = {0};
        uint64_t output_avx2_squared_tst[4] = {0};
        compute_interm_var_four8x8_avx2_intrin(
            input_data, 8, output_avx2_tst, output_avx2_squared_tst);

        // compare results
        EXPECT_EQ(output_avx2_tst[0], output_sse2_tst)
            << "compare mean of 8x8 block error"
            << print_data(input_data, 8, 8);
        EXPECT_EQ(output_avx2_squared_tst[0], output_sse2_squared_tst)
            << "compare mean of 8x8 squared block error"
            << print_data(input_data, 8, 8);
    }
}

}  // namespace
