/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file PictureOperatorTest.cc
 *
 * @brief Unit test for PictureOperatorTest functions:
 * - picture_average_kernel_sse2_intrin
 * - picture_average_kernel1_line_sse2_intrin
 * - picture_copy_kernel_sse2
 *
 * @author Cidana-Ivy
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

#include "EbMcp_SSE2.h"
#include "EbEncIntraPrediction.h"
#include "EbDefinitions.h"
#include "random.h"
#include "util.h"
#include "common_dsp_rtcd.h"
using svt_av1_test_tool::SVTRandom;  // to generate the random

namespace {

typedef enum { REF_MAX, SRC_MAX, RANDOM } TestPattern;
TestPattern TEST_PATTERNS[] = {REF_MAX, SRC_MAX, RANDOM};

typedef std::tuple<uint32_t, uint32_t> PUSize;
PUSize TEST_PU_SIZES[] = {
    PUSize(64, 64), PUSize(64, 32), PUSize(32, 64), PUSize(32, 32),
    PUSize(32, 16), PUSize(16, 32), PUSize(16, 16), PUSize(16, 8),
    PUSize(8, 16),  PUSize(8, 8),   PUSize(8, 4),   PUSize(4, 4),
    PUSize(4, 8),   PUSize(4, 16),  PUSize(16, 4),  PUSize(8, 32),
    PUSize(32, 8),  PUSize(16, 64), PUSize(64, 16), PUSize(24, 24),
    PUSize(24, 16), PUSize(16, 24), PUSize(24, 8),  PUSize(8, 24),
    PUSize(64, 24), PUSize(48, 24), PUSize(32, 24), PUSize(24, 32),
    PUSize(48, 48), PUSize(48, 16), PUSize(48, 32), PUSize(16, 48),
    PUSize(32, 48), PUSize(48, 64), PUSize(64, 48)};

typedef std::tuple<PUSize, TestPattern> TestParam;

/**
 * @brief Unit test for Picture functions in PU size include:
 *  - picture_average_kernel_sse2_intrin
 *  - picture_average_kernel1_line_sse2_intrin
 *  - picture_copy_kernel_sse2
 *
 *
 * Test strategy:
 *  This test case combine different pu width{4-64} x search area
 * height{4-64} and different test pattern(REF_MAX, SRC_MAX, RANDOM). Check
 * the result by compare result from reference function and sse2 function.
 *
 *
 * Expect result:
 *  Results from reference function and sse2 function are
 * equal.
 *
 *
 * Test cases:
 *  Width {4, 8, 16, 24, 32, 48, 64} x height{ 4, 8, 16, 24, 32, 48, 64)
 *  Test vector pattern {REF_MAX, SRC_MAX, RANDOM}
 *
 */
class PictureOperatorTest : public ::testing::Test,
                            public ::testing::WithParamInterface<TestParam> {
  public:
    PictureOperatorTest()
        : pu_width_(std::get<0>(TEST_GET_PARAM(0))),
          pu_height_(std::get<1>(TEST_GET_PARAM(0))),
          test_pattern_(TEST_GET_PARAM(1)) {
        tst_size = 2 * MAX_PU_SIZE * MAX_PU_SIZE;
        tst_stride_ = 2 * MAX_PU_SIZE;
        tst1_aligned_ = (uint8_t *)eb_aom_memalign(8, tst_size);
        tst2_aligned_ = (uint8_t *)eb_aom_memalign(8, tst_size);
        dst1_aligned_ = (uint8_t *)eb_aom_memalign(8, tst_size);
        dst2_aligned_ = (uint8_t *)eb_aom_memalign(8, tst_size);
        memset(dst1_aligned_, 0, tst_size * sizeof(dst1_aligned_[0]));
        memset(dst2_aligned_, 0, tst_size * sizeof(dst2_aligned_[0]));
    }

    void TearDown() override {
        if (tst1_aligned_)
            eb_aom_free(tst1_aligned_);
        if (tst2_aligned_)
            eb_aom_free(tst2_aligned_);
        if (dst1_aligned_)
            eb_aom_free(dst1_aligned_);
        if (dst2_aligned_)
            eb_aom_free(dst2_aligned_);
    }

  protected:
    void prepare_data() {
        const int32_t mask = (1 << 8) - 1;
        SVTRandom rnd(0, mask);
        switch (test_pattern_) {
        case REF_MAX: {
            memset(tst1_aligned_, 0, tst_size * sizeof(tst1_aligned_[0]));
            memset(tst2_aligned_, 0, tst_size * sizeof(tst2_aligned_[0]));
            break;
        }
        case SRC_MAX: {
            memset(tst1_aligned_, mask, tst_size * sizeof(tst1_aligned_[0]));
            memset(tst2_aligned_, mask, tst_size * sizeof(tst2_aligned_[0]));
            break;
        }
        case RANDOM: {
            for (int i = 0; i < tst_size; i++) {
                tst1_aligned_[i] = rnd.random();
                tst2_aligned_[i] = rnd.random();
            }
            break;
        }
        default: break;
        }
    }

    void run_avg_test() {
        prepare_data();
        picture_average_kernel_sse2_intrin(tst1_aligned_,
                                           tst_stride_,
                                           tst2_aligned_,
                                           tst_stride_,
                                           dst1_aligned_,
                                           tst_stride_,
                                           pu_width_,
                                           pu_height_);
        picture_average_kernel_c(tst1_aligned_,
                               tst_stride_,
                               tst2_aligned_,
                               tst_stride_,
                               dst2_aligned_,
                               tst_stride_,
                               pu_width_,
                               pu_height_);

        int fail_pixel_count = 0;
        for (uint16_t j = 0; j < pu_height_; j++) {
            for (uint16_t k = 0; k < pu_width_; k++) {
                if (dst1_aligned_[k + j * tst_stride_] !=
                    dst2_aligned_[k + j * tst_stride_])
                    fail_pixel_count++;
            }
        }
        EXPECT_EQ(0, fail_pixel_count)
            << "compare picture average result error"
            << "in pu for " << fail_pixel_count << "times,"
            << " at func: [picture average] ";

        picture_average_kernel1_line_sse2_intrin(
            tst1_aligned_, tst2_aligned_, dst1_aligned_, pu_width_);
        picture_average_kernel1_line_c(
            tst1_aligned_, tst2_aligned_, dst2_aligned_, pu_width_);

        fail_pixel_count = 0;
        for (uint16_t k = 0; k < pu_width_; k++) {
            if (dst1_aligned_[k] != dst2_aligned_[k])
                fail_pixel_count++;
        }
        EXPECT_EQ(0, fail_pixel_count)
            << "compare picture average 1line result error"
            << "in pu for " << fail_pixel_count << "times,"
            << " at func: [picture average 1line] ";
    }

    void run_copy_test() {
        prepare_data();
        picture_copy_kernel_sse2(tst1_aligned_,
                                 tst_stride_,
                                 dst1_aligned_,
                                 tst_stride_,
                                 pu_width_,
                                 pu_height_);
        picture_copy_kernel(tst1_aligned_,
                            tst_stride_,
                            dst2_aligned_,
                            tst_stride_,
                            pu_width_,
                            pu_height_,
                            1);

        int fail_pixel_count = 0;
        for (uint16_t j = 0; j < pu_height_; j++) {
            for (uint16_t k = 0; k < pu_width_; k++) {
                if (dst1_aligned_[k + j * tst_stride_] !=
                    dst2_aligned_[k + j * tst_stride_])
                    fail_pixel_count++;
            }
        }
        EXPECT_EQ(0, fail_pixel_count)
            << "compare picture copy result error"
            << "in pu for " << fail_pixel_count << "times,"
            << " at func: [picture copy] ";
    }

    int tst_size;
    uint32_t pu_width_, pu_height_;
    TestPattern test_pattern_;
    uint32_t tst_stride_;
    uint8_t *tst1_aligned_, *tst2_aligned_;
    uint8_t *dst1_aligned_, *dst2_aligned_;
};

TEST_P(PictureOperatorTest, AvgTest) {
    run_avg_test();
};

TEST_P(PictureOperatorTest, CopyTest) {
    run_copy_test();
};

INSTANTIATE_TEST_CASE_P(PictureOperator, PictureOperatorTest,
                        ::testing::Combine(::testing::ValuesIn(TEST_PU_SIZES),
                                           ::testing::ValuesIn(TEST_PATTERNS)));

}  // namespace
