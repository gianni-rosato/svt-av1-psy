/*
* Copyright(c) 2019 Netflix, Inc.
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

/******************************************************************************
 * @file ResidualTest.cc
 *
 * @brief Unit test for Residual related functions:
 * - residual_kernel_avx2
 * - svt_residual_kernel16bit_sse2_intrin
 * - residual_kernel_sub_sampled{w}x{h}_sse_intrin
 *
 * @author Cidana-Ivy, Cidana-Wenyao
 *
 ******************************************************************************/
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <limits.h>
#include <new>
// workaround to eliminate the compiling warning on linux
// The macro will conflict with definition in gtest.h
#ifdef __USE_GNU
#undef __USE_GNU  // defined in EbThreads.h
#endif
#ifdef _GNU_SOURCE
#undef _GNU_SOURCE  // defined in EbThreads.h
#endif

#include "EbDefinitions.h"
#include "EbPictureOperators.h"
#include "EbEncIntraPrediction.h"
#include "EbUtility.h"
#include "EbUnitTestUtility.h"
#include "aom_dsp_rtcd.h"
#include "random.h"
#include "util.h"

using svt_av1_test_tool::SVTRandom;  // to generate the random

namespace {
typedef enum { VAL_MIN, VAL_MAX, VAL_RANDOM } TestPattern;
TestPattern TEST_PATTERNS[] = {VAL_MIN, VAL_MAX, VAL_RANDOM};

/**
 * @Brief Base class for Residual related test, it will create input buffer,
 * prediction buffer and two residual buffers to store the results.
 */
class ResidualTestBase : public ::testing::Test {
  public:
    ResidualTestBase(TestPattern test_pattern) {
        test_pattern_ = test_pattern;
    }

    void SetUp() override {
        input_ = (uint8_t *)eb_aom_memalign(8, test_size_);
        pred_ = (uint8_t *)eb_aom_memalign(8, test_size_);
        residual1_ = (int16_t *)eb_aom_memalign(16, 2 * test_size_);
        residual2_ = (int16_t *)eb_aom_memalign(16, 2 * test_size_);
    }

    void TearDown() override {
        if (input_)
            eb_aom_free(input_);
        if (pred_)
            eb_aom_free(pred_);
        if (residual1_)
            eb_aom_free(residual1_);
        if (residual2_)
            eb_aom_free(residual2_);
    }

  protected:
    virtual void prepare_data() {
        const uint8_t mask = (1 << 8) - 1;
        switch (test_pattern_) {
        case VAL_MIN: {
            memset(input_, 0, test_size_ * sizeof(input_[0]));
            memset(pred_, 0, test_size_ * sizeof(pred_[0]));
            break;
        }
        case VAL_MAX: {
            memset(input_, mask, test_size_ * sizeof(input_[0]));
            memset(pred_, 0, test_size_ * sizeof(pred_[0]));
            break;
        }
        case VAL_RANDOM: {
            eb_buf_random_u8(input_, test_size_);
            eb_buf_random_u8(pred_, test_size_);
            break;
        }
        default: break;
        }
    }

    void check_residuals(uint32_t width, uint32_t height) {
        int mismatched_count_ = 0;
        for (uint32_t j = 0; j < height; j++) {
            for (uint32_t k = 0; k < width; k++) {
                if (residual1_[k + j * residual_stride_] !=
                    residual2_[k + j * residual_stride_])
                    ++mismatched_count_;
            }
        }
        EXPECT_EQ(0, mismatched_count_)
            << "compare residual result error"
            << "in test area for " << mismatched_count_ << " times";
    }

    TestPattern test_pattern_;
    uint32_t area_width_, area_height_;
    uint32_t input_stride_, pred_stride_, residual_stride_;
    uint8_t *input_, *pred_;
    int16_t *residual1_, *residual2_;
    uint32_t test_size_;
};

typedef void (*svt_residual_kernel8bit_func)(uint8_t *input, uint32_t input_stride,
                                    uint8_t *pred, uint32_t pred_stride,
                                    int16_t *residual, uint32_t residual_stride,
                                    uint32_t area_width, uint32_t area_height);

static const svt_residual_kernel8bit_func residual_kernel8bit_func_table[] = {
    svt_residual_kernel8bit_avx2,
#ifndef NON_AVX512_SUPPORT
    svt_residual_kernel8bit_avx512
#endif
};

typedef std::tuple<uint32_t, uint32_t> AreaSize;
AreaSize TEST_AREA_SIZES[] = {
    AreaSize(4, 4),    AreaSize(4, 8),   AreaSize(8, 4),   AreaSize(8, 8),
    AreaSize(16, 16),  AreaSize(4, 16),  AreaSize(16, 4),  AreaSize(16, 8),
    AreaSize(8, 16),   AreaSize(32, 32), AreaSize(32, 8),  AreaSize(16, 32),
    AreaSize(8, 32),   AreaSize(32, 16), AreaSize(16, 64), AreaSize(64, 16),
    AreaSize(64, 64),  AreaSize(64, 32), AreaSize(32, 64), AreaSize(128, 128),
    AreaSize(64, 128), AreaSize(128, 64)};
typedef std::tuple<AreaSize, TestPattern> TestResidualParam;

class ResidualKernelTest
    : public ResidualTestBase,
      public ::testing::WithParamInterface<TestResidualParam> {
  public:
    ResidualKernelTest() : ResidualTestBase(TEST_GET_PARAM(1)) {
        area_width_ = std::get<0>(TEST_GET_PARAM(0));
        area_height_ = std::get<1>(TEST_GET_PARAM(0));
        input_stride_ = pred_stride_ = residual_stride_ = MAX_SB_SIZE;
        test_size_ = MAX_SB_SQUARE;
        input16bit_ = (uint16_t *)eb_aom_memalign(16, 2 * test_size_);
        pred16bit_ = (uint16_t *)eb_aom_memalign(16, 2 * test_size_);
    }

    ~ResidualKernelTest() {
        if (input16bit_)
            eb_aom_free(input16bit_);
        if (pred16bit_)
            eb_aom_free(pred16bit_);
    }

  protected:
    void prepare_16bit_data() {
        const uint16_t mask = (1 << 16) - 1;
        switch (test_pattern_) {
        case VAL_MIN: {
            memset(input16bit_, 0, test_size_ * sizeof(input16bit_[0]));
            memset(pred16bit_, 0, test_size_ * sizeof(pred16bit_[0]));
            break;
        }
        case VAL_MAX: {
            for (uint32_t i = 0; i < test_size_; ++i)
                input16bit_[i] = mask;
            memset(pred16bit_, 0, test_size_ * sizeof(pred16bit_[0]));
            break;
        }
        case VAL_RANDOM: {
            eb_buf_random_u16(input16bit_, test_size_);
            eb_buf_random_u16(pred16bit_, test_size_);
            break;
        }
        default: break;
        }
    }

    void run_test() {
        prepare_data();

        svt_residual_kernel8bit_c(input_,
                                  input_stride_,
                                  pred_,
                                  pred_stride_,
                                  residual1_,
                                  residual_stride_,
                                  area_width_,
                                  area_height_);

        for (int i = 0; i < (int) (sizeof(residual_kernel8bit_func_table) /
                                sizeof(*residual_kernel8bit_func_table));
             i++) {
            eb_buf_random_s16(residual2_, test_size_);
            residual_kernel8bit_func_table[i](input_,
                                              input_stride_,
                                              pred_,
                                              pred_stride_,
                                              residual2_,
                                              residual_stride_,
                                              area_width_,
                                              area_height_);
            check_residuals(area_width_, area_height_);
        }
    }

    void speed_test() {
        double time_c, time_o;
        uint64_t start_time_seconds, start_time_useconds;
        uint64_t middle_time_seconds, middle_time_useconds;
        uint64_t finish_time_seconds, finish_time_useconds;

        const uint64_t num_loop = 10000000000 / (area_width_ * area_height_);

        prepare_data();

        svt_av1_get_time(&start_time_seconds, &start_time_useconds);

        for (uint64_t i = 0; i < num_loop; i++) {
            svt_residual_kernel8bit_c(input_,
                                      input_stride_,
                                      pred_,
                                      pred_stride_,
                                      residual1_,
                                      residual_stride_,
                                      area_width_,
                                      area_height_);
        }

        svt_av1_get_time(&middle_time_seconds, &middle_time_useconds);
        time_c = svt_av1_compute_overall_elapsed_time_ms(start_time_seconds,
                                                         start_time_useconds,
                                                         middle_time_seconds,
                                                         middle_time_useconds);

        for (int i = 0; i < (int) (sizeof(residual_kernel8bit_func_table) /
                                sizeof(*residual_kernel8bit_func_table));
             i++) {
            eb_buf_random_s16(residual2_, test_size_);

            svt_av1_get_time(&middle_time_seconds, &middle_time_useconds);

            for (uint64_t j = 0; j < num_loop; j++) {
                residual_kernel8bit_func_table[i](input_,
                                             input_stride_,
                                             pred_,
                                             pred_stride_,
                                             residual2_,
                                             residual_stride_,
                                             area_width_,
                                             area_height_);
            }
            check_residuals(area_width_, area_height_);

            svt_av1_get_time(&finish_time_seconds, &finish_time_useconds);
            time_o =
                svt_av1_compute_overall_elapsed_time_ms(middle_time_seconds,
                                                        middle_time_useconds,
                                                        finish_time_seconds,
                                                        finish_time_useconds);

            printf("svt_residual_kernel8bit(%3dx%3d): %6.2f\n",
                   area_width_,
                   area_height_,
                   time_c / time_o);
        }
    }

    void run_16bit_test() {
        prepare_16bit_data();

        svt_residual_kernel16bit_sse2_intrin(input16bit_,
                                             input_stride_,
                                             pred16bit_,
                                             pred_stride_,
                                             residual1_,
                                             residual_stride_,
                                             area_width_,
                                             area_height_);
        svt_residual_kernel16bit_c(input16bit_,
                                   input_stride_,
                                   pred16bit_,
                                   pred_stride_,
                                   residual2_,
                                   residual_stride_,
                                   area_width_,
                                   area_height_);

        check_residuals(area_width_, area_height_);
    }

    uint16_t *input16bit_, *pred16bit_;
};

TEST_P(ResidualKernelTest, MatchTest) {
    run_test();
};

TEST_P(ResidualKernelTest, DISABLED_SpeedTest) {
    speed_test();
};

TEST_P(ResidualKernelTest, 16bitMatchTest) {
    run_16bit_test();
};

INSTANTIATE_TEST_CASE_P(ResidualUtil, ResidualKernelTest,
                        ::testing::Combine(::testing::ValuesIn(TEST_AREA_SIZES),
                                           ::testing::ValuesIn(TEST_PATTERNS)));

typedef void (*SUB_SAMPLED_TEST_FUNC)(uint8_t *input, uint32_t input_stride,
                                      uint8_t *pred, uint32_t pred_stride,
                                      int16_t *residual,
                                      uint32_t residual_stride,
                                      uint32_t area_width, uint32_t area_height,
                                      uint8_t last_line);

typedef struct SubSampledParam {
    AreaSize area_size;
    SUB_SAMPLED_TEST_FUNC test_func;
} SubSampledParam;
typedef std::tuple<TestPattern, SubSampledParam> TestSubParam;

class ResidualSubSampledTest
    : public ResidualTestBase,
      public ::testing::WithParamInterface<TestSubParam> {
  public:
    ResidualSubSampledTest() : ResidualTestBase(TEST_GET_PARAM(0)) {
        area_width_ = std::get<0>(TEST_GET_PARAM(1).area_size);
        area_height_ = std::get<1>(TEST_GET_PARAM(1).area_size);
        sub_sampled_test_func_ = TEST_GET_PARAM(1).test_func;
        input_stride_ = pred_stride_ = residual_stride_ = MAX_PU_SIZE;
        test_size_ = MAX_PU_SIZE * MAX_PU_SIZE;
    }

  protected:
    void run_sub_sample_test() {
        prepare_data();

        sub_sampled_test_func_(input_,
                               input_stride_,
                               pred_,
                               pred_stride_,
                               residual1_,
                               residual_stride_,
                               area_width_,
                               area_height_,
                               false);

        check_residuals(area_width_, area_height_);
    }

    SUB_SAMPLED_TEST_FUNC sub_sampled_test_func_;
};

TEST_P(ResidualSubSampledTest, MatchTest) {
    run_sub_sample_test();
};

typedef std::tuple<int, TestPattern> TestSumParam;

class ResidualSumTest : public ::testing::Test,
                        public ::testing::WithParamInterface<TestSumParam> {
  public:
    ResidualSumTest() {
        test_pattern_ = TEST_GET_PARAM(1);
        size_ = TEST_GET_PARAM(0);
        residual_stride_ = MAX_SB_SIZE;
        residual_ = (int16_t *)eb_aom_memalign(16, 2 * MAX_SB_SQUARE);
    }

    ~ResidualSumTest() {
        if (residual_)
            eb_aom_free(residual_);
    }

  protected:
    void prepare_data() {
        // There is an assumption in sum_residual8bit_avx2_intrin, which is
        // 9bit or 11bit residual data.
        const int16_t mask = (1 << 10) - 1;
        SVTRandom rnd(-(1 << 10), mask);
        switch (test_pattern_) {
        case VAL_MIN: {
            memset(residual_, 0, MAX_SB_SQUARE * sizeof(residual_[0]));
            break;
        }
        case VAL_MAX: {
            for (uint32_t i = 0; i < MAX_SB_SQUARE; ++i)
                residual_[i] = mask;
            break;
        }
        case VAL_RANDOM: {
            for (uint32_t i = 0; i < MAX_SB_SQUARE; ++i)
                residual_[i] = rnd.random();
            break;
        }
        default: break;
        }
    }
    uint32_t size_;
    TestPattern test_pattern_;
    int16_t *residual_;
    uint32_t residual_stride_;
};
INSTANTIATE_TEST_CASE_P(ResidualUtil, ResidualSumTest,
                        ::testing::Combine(::testing::Values(4, 8, 16, 32, 64),
                                           ::testing::ValuesIn(TEST_PATTERNS)));
}  // namespace
