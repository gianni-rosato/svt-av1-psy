/*
 * Copyright(c) 2019 Netflix, Inc.
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at https://www.aomedia.org/license/software-license. If the
 * Alliance for Open Media Patent License 1.0 was not distributed with this
 * source code in the PATENTS file, you can obtain it at
 * https://www.aomedia.org/license/patent-license.
 */

/******************************************************************************
 * @file PictureOperatorTest.cc
 *
 * @brief Unit test for PictureOperatorTest functions:
 * - svt_picture_average_kernel_sse2_intrin
 * - svt_picture_average_kernel1_line_sse2_intrin
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

#include "mcp_sse2.h"
#include "enc_intra_prediction.h"
#include "definitions.h"
#include "random.h"
#include "util.h"
#include "common_dsp_rtcd.h"
#include "aom_dsp_rtcd.h"
#include "motion_estimation.h"
using svt_av1_test_tool::SVTRandom;  // to generate the random

namespace {

typedef enum { REF_MAX, SRC_MAX, RANDOM } TestPattern;

typedef std::tuple<uint32_t, uint32_t> PUSize;

#ifdef ARCH_X86_64
TestPattern TEST_PATTERNS[] = {REF_MAX, SRC_MAX, RANDOM};

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
#endif  // ARCH_X86_64

using PictureOperatorFunc = void (*)(EbByte src0, uint32_t src0_stride,
                                     EbByte src1, uint32_t src1_stride,
                                     EbByte dst, uint32_t dst_stride,
                                     uint32_t area_width, uint32_t area_height);

using PictureOperatorFunc1Line = void (*)(EbByte src0, EbByte src1, EbByte dst,
                                          uint32_t area_width);

typedef std::tuple<PUSize, TestPattern, PictureOperatorFunc,
                   PictureOperatorFunc1Line>
    TestParam;

/**
 * @brief Unit test for Picture functions in PU size include:
 *  - svt_picture_average_kernel_sse2_intrin
 *  - svt_picture_average_kernel1_line_sse2_intrin
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
          test_pattern_(TEST_GET_PARAM(1)),
          test_func_(TEST_GET_PARAM(2)),
          test_func_1line_(TEST_GET_PARAM(3)) {
        tst_size = 2 * 64 * 64;
        tst_stride_ = 2 * 64;
        tst1_aligned_ = (uint8_t *)svt_aom_memalign(8, tst_size);
        tst2_aligned_ = (uint8_t *)svt_aom_memalign(8, tst_size);
        dst1_aligned_ = (uint8_t *)svt_aom_memalign(8, tst_size);
        dst2_aligned_ = (uint8_t *)svt_aom_memalign(8, tst_size);
        memset(dst1_aligned_, 0, tst_size * sizeof(dst1_aligned_[0]));
        memset(dst2_aligned_, 0, tst_size * sizeof(dst2_aligned_[0]));
    }

    void TearDown() override {
        if (tst1_aligned_)
            svt_aom_free(tst1_aligned_);
        if (tst2_aligned_)
            svt_aom_free(tst2_aligned_);
        if (dst1_aligned_)
            svt_aom_free(dst1_aligned_);
        if (dst2_aligned_)
            svt_aom_free(dst2_aligned_);
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
        test_func_(tst1_aligned_,
                   tst_stride_,
                   tst2_aligned_,
                   tst_stride_,
                   dst1_aligned_,
                   tst_stride_,
                   pu_width_,
                   pu_height_);
        svt_picture_average_kernel_c(tst1_aligned_,
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

        svt_picture_average_kernel1_line_c(
            tst1_aligned_, tst2_aligned_, dst2_aligned_, pu_width_);
        test_func_1line_(
            tst1_aligned_, tst2_aligned_, dst1_aligned_, pu_width_);

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

    int tst_size;
    uint32_t pu_width_, pu_height_;
    TestPattern test_pattern_;
    PictureOperatorFunc test_func_;
    PictureOperatorFunc1Line test_func_1line_;
    uint32_t tst_stride_;
    uint8_t *tst1_aligned_, *tst2_aligned_;
    uint8_t *dst1_aligned_, *dst2_aligned_;
};
GTEST_ALLOW_UNINSTANTIATED_PARAMETERIZED_TEST(PictureOperatorTest);

TEST_P(PictureOperatorTest, AvgTest) {
    run_avg_test();
};

#if defined(ARCH_X86_64)
INSTANTIATE_TEST_SUITE_P(
    SSE2, PictureOperatorTest,
    ::testing::Combine(
        ::testing::ValuesIn(TEST_PU_SIZES), ::testing::ValuesIn(TEST_PATTERNS),
        ::testing::Values(svt_picture_average_kernel_sse2_intrin),
        ::testing::Values(svt_picture_average_kernel1_line_sse2_intrin)));

#endif  // ARCH_X86_64

typedef void (*downsample_2d_fn)(uint8_t *input_samples, uint32_t input_stride,
                                 uint32_t input_area_width,
                                 uint32_t input_area_height,
                                 uint8_t *decim_samples, uint32_t decim_stride,
                                 uint32_t decim_step);
uint32_t DECIM_STEPS[] = {2, 4, 8};
PUSize DOWNSAMPLE_SIZES[] = {
    PUSize(1920, 1080), PUSize(960, 540), PUSize(176, 144), PUSize(88, 72)};

typedef std::tuple<PUSize, uint32_t, downsample_2d_fn> Downsample2DParam;

class Downsample2DTest
    : public ::testing::Test,
      public ::testing::WithParamInterface<Downsample2DParam> {
  public:
    Downsample2DTest()
        : pu_width(std::get<0>(TEST_GET_PARAM(0))),
          pu_height(std::get<1>(TEST_GET_PARAM(0))),
          decim_step(TEST_GET_PARAM(1)),
          fn_ptr(TEST_GET_PARAM(2)) {
        max_size = sizeof(uint8_t) * (1920 + 3) * (1080 + 3);
        stride = pu_width + 3;
        decim_stride = (pu_width / decim_step) + 3;
        src_ptr = (uint8_t *)malloc(max_size);
        dst_ref_ptr = (uint8_t *)malloc(max_size);
        dst_tst_ptr = (uint8_t *)malloc(max_size);
    }

    void TearDown() override {
        if (src_ptr)
            free(src_ptr);
        if (dst_ref_ptr)
            free(dst_ref_ptr);
        if (dst_tst_ptr)
            free(dst_tst_ptr);
    }

  protected:
    void prepare_data() {
        const int32_t mask = (1 << 8) - 1;
        SVTRandom rnd(0, mask);
        for (int i = 0; i < max_size; i++) {
            src_ptr[i] = rnd.random();
        }
        uint8_t val = rnd.random();
        memset(dst_ref_ptr, val, max_size);
        memset(dst_tst_ptr, val, max_size);
    }

    void run_test() {
        prepare_data();
        svt_aom_downsample_2d_c(src_ptr,
                                stride,
                                pu_width,
                                pu_height,
                                dst_ref_ptr,
                                decim_stride,
                                decim_step);

        fn_ptr(src_ptr,
               stride,
               pu_width,
               pu_height,
               dst_tst_ptr,
               decim_stride,
               decim_step);

        EXPECT_EQ(memcmp(dst_ref_ptr, dst_tst_ptr, max_size), 0);
    }

    int max_size;
    uint32_t pu_width, pu_height;
    uint32_t decim_step;
    uint32_t decim_stride;
    uint32_t stride;
    uint8_t *src_ptr;
    uint8_t *dst_ref_ptr, *dst_tst_ptr;
    downsample_2d_fn fn_ptr;
};

TEST_P(Downsample2DTest, test) {
    for (int i = 0; i < 20; i++) {
        run_test();
    }
};

#if defined(ARCH_X86_64)

INSTANTIATE_TEST_SUITE_P(
    SSE4_1, Downsample2DTest,
    ::testing::Combine(::testing::ValuesIn(DOWNSAMPLE_SIZES),
                       ::testing::ValuesIn(DECIM_STEPS),
                       ::testing::Values(svt_aom_downsample_2d_sse4_1)));

INSTANTIATE_TEST_SUITE_P(
    AVX2, Downsample2DTest,
    ::testing::Combine(::testing::ValuesIn(DOWNSAMPLE_SIZES),
                       ::testing::ValuesIn(DECIM_STEPS),
                       ::testing::Values(svt_aom_downsample_2d_avx2)));

#endif

#if defined(ARCH_AARCH64)

INSTANTIATE_TEST_SUITE_P(
    NEON, Downsample2DTest,
    ::testing::Combine(::testing::ValuesIn(DOWNSAMPLE_SIZES),
                       ::testing::ValuesIn(DECIM_STEPS),
                       ::testing::Values(svt_aom_downsample_2d_neon)));

#endif

}  // namespace
