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
 * @file sad_Test.cc
 *
 * @brief Unit test for SAD functions:
 * - svt_nxm_sad_kernel_sub_sampled_func
 * - svt_nxm_sad_kernel_func
 * - nxm_sad_averaging_kernel_func
 * - nxm_sad_loop_kernel_sparse_func
 * - nxm_sad_loop_kernel_sparse_func
 * - get_eight_horizontal_search_point_results_8x8_16x16_func
 * - get_eight_horizontal_search_point_results_32x32_64x64_func
 * - svt_ext_ext_all_sad_calculation_8x8_16x16_func
 * - svt_ext_ext_eight_sad_calculation_32x32_64x64_func
 * - svt_ext_eigth_sad_calculation_nsq_func
 * - Extsad_Calculation_8x8_16x16_func
 * - Extsad_Calculation_32x32_64x64_func
 *
 * @author Cidana-Ryan, Cidana-Wenyao, Cidana-Ivy
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
#include "gtest/gtest.h"
#include "aom_dsp_rtcd.h"
#include "EbComputeSAD.h"
#include "EbMeSadCalculation.h"
#include "EbMotionEstimation.h"
#include "EbMotionEstimationContext.h"
#include "EbTime.h"
#include "random.h"
#include "util.h"

using svt_av1_test_tool::SVTRandom;  // to generate the random
extern "C" void svt_ext_all_sad_calculation_8x8_16x16_c(
    uint8_t *src, uint32_t src_stride, uint8_t *ref, uint32_t ref_stride,
    uint32_t mv, uint32_t *p_best_sad_8x8, uint32_t *p_best_sad_16x16,
    uint32_t *p_best_mv8x8, uint32_t *p_best_mv16x16,
    uint32_t p_eight_sad16x16[16][8], uint32_t p_eight_sad8x8[64][8], EbBool sub_sad);
extern "C" void svt_ext_eigth_sad_calculation_nsq_c(
    uint32_t p_sad8x8[64][8], uint32_t p_sad16x16[16][8],
    uint32_t p_sad32x32[4][8], uint32_t *p_best_sad_64x32,
    uint32_t *p_best_mv64x32, uint32_t *p_best_sad_32x16,
    uint32_t *p_best_mv32x16, uint32_t *p_best_sad_16x8,
    uint32_t *p_best_mv16x8, uint32_t *p_best_sad_32x64,
    uint32_t *p_best_mv32x64, uint32_t *p_best_sad_16x32,
    uint32_t *p_best_mv16x32, uint32_t *p_best_sad_8x16,
    uint32_t *p_best_mv8x16, uint32_t *p_best_sad_32x8, uint32_t *p_best_mv32x8,
    uint32_t *p_best_sad_8x32, uint32_t *p_best_mv8x32,
    uint32_t *p_best_sad_64x16, uint32_t *p_best_mv64x16,
    uint32_t *p_best_sad_16x64, uint32_t *p_best_mv16x64, uint32_t mv);
extern "C" void svt_ext_eight_sad_calculation_32x32_64x64_c(
    uint32_t p_sad16x16[16][8], uint32_t *p_best_sad_32x32,
    uint32_t *p_best_sad_64x64, uint32_t *p_best_mv32x32,
    uint32_t *p_best_mv64x64, uint32_t mv, uint32_t p_sad32x32[4][8]);

namespace {
/**
 * @Brief Test param definition.
 * - TEST_BLOCK_SIZES:All of sad_ calculation funcs in test cover block size
 *  Width {4, 8, 16, 24, 32, 48, 64} x height{ 4, 8, 16, 24, 32, 48, 64).
 */
#define MAX_BLOCK_SIZE (MAX_SB_SIZE * MAX_SB_SIZE)
#define MAX_REF_BLOCK_SIZE \
    ((MAX_SB_SIZE + PAD_VALUE) * (MAX_SB_SIZE + PAD_VALUE))
typedef std::tuple<int, int> BlkSize;
typedef enum { REF_MAX, SRC_MAX, RANDOM, UNALIGN } TestPattern;
typedef enum { BUF_MAX, BUF_MIN, BUF_SMALL, BUF_RANDOM } SADPattern;
BlkSize TEST_BLOCK_SIZES[] = {
    BlkSize(16, 10), BlkSize(16, 5),
    BlkSize(32, 10), BlkSize(32, 20), BlkSize(64, 20),
    BlkSize(64, 64), BlkSize(64, 32), BlkSize(32, 64), BlkSize(32, 32),
    BlkSize(32, 16), BlkSize(16, 32), BlkSize(16, 16), BlkSize(16, 8),
    BlkSize(8, 16),  BlkSize(8, 8),   BlkSize(8, 4),   BlkSize(4, 4),
    BlkSize(4, 8),   BlkSize(4, 16),  BlkSize(16, 4),  BlkSize(8, 32),
    BlkSize(32, 8),  BlkSize(16, 64), BlkSize(64, 16), BlkSize(24, 24),
    BlkSize(24, 16), BlkSize(16, 24), BlkSize(24, 8),  BlkSize(8, 24),
    BlkSize(64, 24), BlkSize(48, 24), BlkSize(32, 24), BlkSize(24, 32),
    BlkSize(48, 48), BlkSize(48, 16), BlkSize(48, 32), BlkSize(16, 48),
    BlkSize(32, 48), BlkSize(48, 64), BlkSize(64, 48), BlkSize(56, 32),
    BlkSize(40, 32)};

BlkSize TEST_BLOCK_SIZES_SMALL[] = {
    BlkSize(6, 2),   BlkSize(6, 4),   BlkSize(6, 8),   BlkSize(6, 16),
    BlkSize(6, 32),  BlkSize(12, 2),  BlkSize(12, 4),  BlkSize(12, 8),
    BlkSize(12, 16), BlkSize(12, 32), BlkSize(31, 1),  BlkSize(31, 2),
    BlkSize(31, 3),  BlkSize(15, 10), BlkSize(6, 10),  BlkSize(7, 10),
    BlkSize(5, 11),  BlkSize(6, 11),  BlkSize(4, 11),  BlkSize(8, 11),
    BlkSize(7, 11),  BlkSize(5, 10),  BlkSize(16, 10), BlkSize(17, 10),
    BlkSize(15, 11), BlkSize(16, 11), BlkSize(17, 11), BlkSize(18, 11),
    BlkSize(19, 11), BlkSize(31, 8),  BlkSize(16, 5),  BlkSize(16, 17),
    BlkSize(16, 31), BlkSize(16, 33), BlkSize(12, 5),  BlkSize(12, 17),
    BlkSize(12, 31), BlkSize(12, 33), BlkSize(31, 7),  BlkSize(31, 6),
    BlkSize(31, 5),  BlkSize(31, 4),  BlkSize(39, 1),  BlkSize(3, 40),
    BlkSize(43, 4),  BlkSize(43, 5),  BlkSize(41, 5),  BlkSize(55, 3),
    BlkSize(37, 37), BlkSize(41, 21), BlkSize(51, 21), BlkSize(63, 21),
    BlkSize(63, 27), BlkSize(63, 33), BlkSize(63, 32), BlkSize(4, 2),
    BlkSize(4, 3),   BlkSize(4, 4),   BlkSize(4, 8),   BlkSize(4, 9),
    BlkSize(4, 16),  BlkSize(4, 17),  BlkSize(4, 32),  BlkSize(4, 33),
    BlkSize(6, 3),   BlkSize(6, 9),   BlkSize(6, 17),  BlkSize(6, 33),
    BlkSize(8, 2),   BlkSize(8, 3),   BlkSize(8, 4),   BlkSize(8, 8),
    BlkSize(8, 9),   BlkSize(8, 15),  BlkSize(8, 16),  BlkSize(8, 31),
    BlkSize(8, 32),  BlkSize(12, 3),  BlkSize(12, 9),  BlkSize(12, 15),
    BlkSize(12, 31), BlkSize(16, 2),  BlkSize(16, 3),  BlkSize(16, 4),
    BlkSize(16, 8),  BlkSize(16, 9),  BlkSize(16, 15), BlkSize(16, 16),
    BlkSize(16, 31), BlkSize(16, 32), BlkSize(24, 2),  BlkSize(24, 3),
    BlkSize(24, 4),  BlkSize(24, 8),  BlkSize(24, 9),  BlkSize(24, 15),
    BlkSize(24, 16), BlkSize(24, 31), BlkSize(24, 32), BlkSize(32, 2),
    BlkSize(32, 3),  BlkSize(32, 4),  BlkSize(32, 8),  BlkSize(32, 9),
    BlkSize(32, 15), BlkSize(32, 16), BlkSize(32, 31), BlkSize(32, 32),
    BlkSize(32, 33),
};
TestPattern TEST_PATTERNS[] = {REF_MAX, SRC_MAX, RANDOM, UNALIGN};
SADPattern TEST_SAD_PATTERNS[] = {BUF_MAX, BUF_MIN, BUF_SMALL, BUF_RANDOM};
typedef std::tuple<TestPattern, BlkSize> Testsad_Param;

/**
 * @Brief Base class for SAD test. SADTestBase handle test vector in memory,
 * provide SAD and SAD avg reference function
 */
class SADTestBase : public ::testing::Test {
  public:
    SADTestBase(const int width, const int height, TestPattern test_pattern) {
        width_ = width;
        height_ = height;
        src_stride_ = MAX_SB_SIZE;
        ref1_stride_ = ref2_stride_ = MAX_SB_SIZE;
        test_pattern_ = test_pattern;
        src_aligned_ = nullptr;
        ref1_aligned_ = nullptr;
        ref2_aligned_ = nullptr;
    }

    SADTestBase(const int width, const int height, TestPattern test_pattern,
                SADPattern test_sad_pattern) {
        width_ = width;
        height_ = height;
        src_stride_ = MAX_SB_SIZE;
        ref1_stride_ = ref2_stride_ = MAX_SB_SIZE;
        test_pattern_ = test_pattern;
        test_sad_pattern_ = test_sad_pattern;
    }

    SADTestBase(const int width, const int height, TestPattern test_pattern,
                const int search_area_width, const int search_area_height) {
        width_ = width;
        height_ = height;
        src_stride_ = MAX_SB_SIZE;
        ref1_stride_ = ref2_stride_ = MAX_SB_SIZE;
        test_pattern_ = test_pattern;
        search_area_width_ = search_area_width;
        search_area_height_ = search_area_height;
    }

    void SetUp() override {
        src_aligned_ = (uint8_t *)svt_aom_memalign(32, MAX_BLOCK_SIZE);
        ref1_aligned_ = (uint8_t *)svt_aom_memalign(32, MAX_REF_BLOCK_SIZE);
        ref2_aligned_ = (uint8_t *)svt_aom_memalign(32, MAX_REF_BLOCK_SIZE);
        ASSERT_NE(src_aligned_, nullptr);
        ASSERT_NE(ref1_aligned_, nullptr);
        ASSERT_NE(ref2_aligned_, nullptr);
    }

    void TearDown() override {
        if (src_aligned_)
            svt_aom_free(src_aligned_);
        if (ref1_aligned_)
            svt_aom_free(ref1_aligned_);
        if (ref2_aligned_)
            svt_aom_free(ref2_aligned_);
    }

    void prepare_data() {
        const int32_t mask = (1 << 8) - 1;
        SVTRandom rnd(0, mask);
        switch (test_pattern_) {
        case REF_MAX: {
            for (int i = 0; i < MAX_BLOCK_SIZE; i++)
                src_aligned_[i] = 0;

            for (int i = 0; i < MAX_REF_BLOCK_SIZE; i++)
                ref1_aligned_[i] = ref2_aligned_[i] = mask;

            break;
        }
        case SRC_MAX: {
            for (int i = 0; i < MAX_BLOCK_SIZE; i++)
                src_aligned_[i] = mask;

            for (int i = 0; i < MAX_REF_BLOCK_SIZE; i++)
                ref1_aligned_[i] = ref2_aligned_[i] = 0;

            break;
        }
        case RANDOM: {
            for (int i = 0; i < MAX_BLOCK_SIZE; i++)
                src_aligned_[i] = rnd.random();

            for (int i = 0; i < MAX_REF_BLOCK_SIZE; ++i) {
                ref1_aligned_[i] = rnd.random();
                ref2_aligned_[i] = rnd.random();
            }
            break;
        };
        case UNALIGN: {
            for (int i = 0; i < MAX_BLOCK_SIZE; i++)
                src_aligned_[i] = rnd.random();

            for (int i = 0; i < MAX_REF_BLOCK_SIZE; ++i) {
                ref1_aligned_[i] = rnd.random();
                ref2_aligned_[i] = rnd.random();
            }
            ref1_stride_ -= 1;
            ref2_stride_ -= 1;
            break;
        }
        default: break;
        }
    }

    void fill_buf_with_value_16b(uint16_t *buf, int num, uint32_t value) {
        for (int i = 0; i < num; ++i)
            buf[i] = value;
    }

    void fill_buf_with_value(uint32_t *buf, int num, uint32_t value) {
        for (int i = 0; i < num; ++i)
            buf[i] = value;
    }

    void prepare_sad_data_16b(uint32_t best_sad32x32[4]) {
        const int32_t max = (1 << 16) - 1;
        SVTRandom rnd1(0, 4 * max);

        for (int i = 0; i < 4; i++)
            best_sad32x32[i] = rnd1.random();

        switch (test_sad_pattern_) {
        case BUF_MAX: {
            fill_buf_with_value_16b(&sad16x16_16b[0][0], 16 * 8, max);
            break;
        }
        case BUF_MIN: {
            fill_buf_with_value_16b(&sad16x16_16b[0][0], 16 * 8, 0);
            break;
        }
        case BUF_SMALL: {
            const int32_t mask = 256;
            SVTRandom rnd_small(0, mask);
            for (int i = 0; i < 16; i++)
                for (int j = 0; j < 8; j++)
                    sad16x16_16b[i][j] = rnd_small.random();
            break;
        }
        case BUF_RANDOM: {
            SVTRandom rnd(0, max);
            for (int i = 0; i < 16; i++)
                for (int j = 0; j < 8; j++)
                    sad16x16_16b[i][j] = rnd.random();
            break;
        }
        default: break;
        }
    }

    void prepare_sad_data_32b() {
        const int32_t mask = (1 << 8) - 1;
        SVTRandom rnd(0, mask);
        switch (test_sad_pattern_) {
        case BUF_MAX: {
            fill_buf_with_value(&sad16x16_32b[0][0], 16 * 8, mask);
            break;
        }
        case BUF_MIN: {
            fill_buf_with_value(&sad16x16_32b[0][0], 16 * 8, 0);
            break;
        }
        case BUF_SMALL: {
            SVTRandom rnd_small(0, 256);
            for (int i = 0; i < 16; i++)
                for (int j = 0; j < 8; j++)
                    sad16x16_32b[i][j] = rnd_small.random();
            break;
        }
        case BUF_RANDOM: {
            for (int i = 0; i < 16; i++)
                for (int j = 0; j < 8; j++)
                    sad16x16_32b[i][j] = rnd.random();
            break;
        }
        default: break;
        }
    }

    void prepare_nsq_sad_data() {
        const int32_t mask = (1 << 8) - 1;
        SVTRandom rnd(0, mask);
        switch (test_sad_pattern_) {
        case BUF_MAX: {
            fill_buf_with_value(&sad8x8[0][0], 64 * 8, mask);
            fill_buf_with_value(&sad16x16_32b[0][0], 16 * 8, mask);
            fill_buf_with_value(&sad32x32[0][0], 4 * 8, mask);
            break;
        }
        case BUF_MIN: {
            fill_buf_with_value(&sad8x8[0][0], 64 * 8, 0);
            fill_buf_with_value(&sad16x16_32b[0][0], 16 * 8, 0);
            fill_buf_with_value(&sad32x32[0][0], 4 * 8, 0);
            break;
        }
        case BUF_SMALL: {
            SVTRandom rnd_small(0, 256);
            for (int i = 0; i < 64; i++)
                for (int j = 0; j < 8; j++)
                    sad8x8[i][j] = rnd_small.random();

            for (int i = 0; i < 16; i++)
                for (int j = 0; j < 8; j++)
                    sad16x16_32b[i][j] = rnd_small.random();

            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 8; j++)
                    sad32x32[i][j] = rnd_small.random();
            break;
        }
        case BUF_RANDOM: {
            for (int i = 0; i < 64; i++)
                for (int j = 0; j < 8; j++)
                    sad8x8[i][j] = rnd.random();

            for (int i = 0; i < 16; i++)
                for (int j = 0; j < 8; j++)
                    sad16x16_32b[i][j] = rnd.random();

            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 8; j++)
                    sad32x32[i][j] = rnd.random();
            break;
        }
        default: break;
        }
    }

    uint32_t reference_sad() {
        unsigned int sad = 0;
        for (int h = 0; h < height_; ++h) {
            for (int w = 0; w < width_; ++w) {
                sad += abs(src_aligned_[h * src_stride_ + w] -
                           ref1_aligned_[h * ref1_stride_ + w]);
            }
        }
        return sad;
    }

    unsigned int reference_sad_avg() {
        unsigned int sad = 0;
        for (int h = 0; h < height_; ++h) {
            for (int w = 0; w < width_; ++w) {
                const int tmp = ref2_aligned_[h * ref2_stride_ + w] +
                                ref1_aligned_[h * ref1_stride_ + w];
                const uint8_t comp_pred = ROUND_POWER_OF_TWO(tmp, 1);
                sad += abs(src_aligned_[h * src_stride_ + w] - comp_pred);
            }
        }
        return sad;
    }

  protected:
    int width_, height_;
    int src_stride_;
    int ref1_stride_;
    int ref2_stride_;
    int16_t search_area_height_, search_area_width_;
    TestPattern test_pattern_;
    SADPattern test_sad_pattern_;
    uint8_t *src_aligned_;
    uint8_t *ref1_aligned_;
    uint8_t *ref2_aligned_;
    uint16_t sad16x16_16b[16][8];
    uint32_t sad8x8[64][8];
    uint32_t sad16x16_32b[16][8];
    uint32_t sad32x32[4][8];
};

/**
 * @brief Unit test for SAD sub smaple functions include:
 *  - svt_nxm_sad_kernel_helper_c
 *  - svt_nxm_sad_kernel_sub_sampled_helper_avx2
 *
 * Test strategy:
 *  This test case combine different width{4-64} x height{4-64} and different
 * test pattern(REF_MAX, SRC_MAX, RANDOM, UNALIGN). Check the result by compare
 *  result from reference function, non_avx2 function and avx2 function.
 *
 *
 * Expect result:
 *  Results from reference functon, non_avx2 function and avx2 funtion are
 * equal.
 *
 * Test coverage:
 *  All functions inside svt_nxm_sad_kernel_helper_c and
 * svt_nxm_sad_kernel_sub_sampled_helper_avx2.
 *
 * Test cases:
 *  Width {4, 8, 16, 24, 32, 48, 64} x height{ 4, 8, 16, 24, 32, 48, 64)
 *  Test vector pattern {REF_MAX, SRC_MAX, RANDOM, UNALIGN}
 *
 */
class SADTestSubSample : public ::testing::WithParamInterface<Testsad_Param>,
                         public SADTestBase {
  public:
    SADTestSubSample()
        : SADTestBase(std::get<0>(TEST_GET_PARAM(1)),
                      std::get<1>(TEST_GET_PARAM(1)), TEST_GET_PARAM(0)) {
    }

  protected:
    void check_sad() {
        uint32_t ref_sad = 0;
        uint32_t non_avx2_sad = 0;
        uint32_t avx2_sad = 0;

        prepare_data();

        ref_sad = reference_sad();
        non_avx2_sad = svt_nxm_sad_kernel_helper_c(src_aligned_,
                                               src_stride_,
                                               ref1_aligned_,
                                               ref1_stride_,
                                               height_,
                                               width_);
        avx2_sad = svt_nxm_sad_kernel_sub_sampled_helper_avx2(src_aligned_,
                                                          src_stride_,
                                                          ref1_aligned_,
                                                          ref1_stride_,
                                                          height_,
                                                          width_);
        EXPECT_EQ(non_avx2_sad, avx2_sad)
            << "compare non_avx2_sad(" << non_avx2_sad << ") and avx2_sad("
            << avx2_sad << ") error, ref: " << ref_sad;
    }
};

TEST_P(SADTestSubSample, SADTestSubSample) {
    check_sad();
}

INSTANTIATE_TEST_CASE_P(
    SAD, SADTestSubSample,
    ::testing::Combine(::testing::ValuesIn(TEST_PATTERNS),
                       ::testing::ValuesIn(TEST_BLOCK_SIZES)));
/**
 * @brief Unit test for SAD functions include:
 *  - svt_nxm_sad_kernel_helper_c
 *  - svt_nxm_sad_kernel_helper_avx2
 *
 * Test strategy:
 *  This test case combine different wight{4-64} x height{4-64}, different test
 *  vecotr pattern {REF_MAX, SRC_MAX, RANDOM, UNALIGN} to generate test vector.
 *  Check the result by compare result from reference function, non_avx2
 * function and avx2 function.
 *
 *
 * Expect result:
 *  Results from reference functon, non_avx2 function and avx2 funtion are
 *  equal.
 *
 * Test coverage:
 *  All functions inside svt_nxm_sad_kernel_helper_c and
 *  svt_nxm_sad_kernel_helper_avx2.
 *
 * Test cases:
 *  Width {4, 8, 16, 24, 32, 48, 64} x height{ 4, 8, 16, 24, 32, 48, 64)
 *  Test vector pattern {REF_MAX, SRC_MAX, RANDOM, UNALIGN}.
 *
 */
class SADTest : public ::testing::WithParamInterface<Testsad_Param>,
                public SADTestBase {
  public:
    SADTest()
        : SADTestBase(std::get<0>(TEST_GET_PARAM(1)),
                      std::get<1>(TEST_GET_PARAM(1)), TEST_GET_PARAM(0)) {
    }

  protected:
    void check_sad() {
        uint32_t ref_sad = 0;
        uint32_t non_avx2_sad = 0;
        uint32_t avx2_sad = 0;

        prepare_data();

        ref_sad = reference_sad();
        non_avx2_sad = svt_nxm_sad_kernel_helper_c(src_aligned_,
                                                   src_stride_,
                                                   ref1_aligned_,
                                                   ref1_stride_,
                                                   height_,
                                                   width_);
        avx2_sad = svt_nxm_sad_kernel_helper_avx2(src_aligned_,
                                                  src_stride_,
                                                  ref1_aligned_,
                                                  ref1_stride_,
                                                  height_,
                                                  width_);
        EXPECT_EQ(non_avx2_sad, avx2_sad)
            << "compare non_avx2_sad(" << non_avx2_sad << ") and avx2_sad("
            << avx2_sad << ") error, ref: " << ref_sad;
    }
};

TEST_P(SADTest, SADTest) {
    check_sad();
}

INSTANTIATE_TEST_CASE_P(
    SAD, SADTest,
    ::testing::Combine(::testing::ValuesIn(TEST_PATTERNS),
                       ::testing::ValuesIn(TEST_BLOCK_SIZES)));

typedef std::tuple<int16_t, int16_t> SearchArea;
#if !REMOVE_ME_SUBPEL_CODE
SearchArea TEST_AREAS[] = {
    SearchArea(64, 125),  SearchArea(192, 75),  SearchArea(128, 50),
    SearchArea(64, 25),   SearchArea(240, 200), SearchArea(144, 120),
    SearchArea(96, 80),   SearchArea(48, 40),   SearchArea(240, 120),
    SearchArea(144, 72),  SearchArea(96, 48),   SearchArea(48, 24),
    SearchArea(560, 320), SearchArea(336, 192), SearchArea(224, 128),
    SearchArea(112, 64),  SearchArea(640, 400), SearchArea(384, 240),
    SearchArea(256, 160), SearchArea(128, 80),  SearchArea(480, 120),
    SearchArea(288, 72),  SearchArea(192, 48),  SearchArea(96, 24),
    SearchArea(160, 60),  SearchArea(96, 36),   SearchArea(64, 24),
    SearchArea(32, 12)};
#endif
SearchArea TEST_LOOP_AREAS[] = {
    SearchArea(8, 15),    SearchArea(16, 31),   SearchArea(12, 31),
    SearchArea(64, 125),  SearchArea(192, 75),  SearchArea(128, 50),
    SearchArea(64, 25),   SearchArea(240, 200), SearchArea(144, 120),
    SearchArea(96, 80),   SearchArea(48, 40),   SearchArea(240, 120),
    SearchArea(144, 72),  SearchArea(96, 48),   SearchArea(48, 24),
    SearchArea(560, 320), SearchArea(336, 192), SearchArea(224, 128),
    SearchArea(112, 64),  SearchArea(640, 400), SearchArea(384, 240),
    SearchArea(256, 160), SearchArea(128, 80),  SearchArea(480, 120),
    SearchArea(288, 72),  SearchArea(192, 48),  SearchArea(96, 24),
    SearchArea(160, 60),  SearchArea(96, 36),   SearchArea(64, 24),
    SearchArea(32, 12),   SearchArea(15, 6)};

typedef void (*Ebsad_LoopKernelNxMType)(
    uint8_t *src,         // input parameter, source samples Ptr
    uint32_t src_stride,  // input parameter, source stride
    uint8_t *ref,         // input parameter, reference samples Ptr
    uint32_t ref_stride,  // input parameter, reference stride
    uint32_t height,      // input parameter, block height (M)
    uint32_t width,       // input parameter, block width (N)
    uint64_t *best_sad, int16_t *x_search_center, int16_t *y_search_center,
    uint32_t
        src_stride_raw,  // input parameter, source stride (no line skipping)
    int16_t search_area_width, int16_t search_area_height);

typedef std::tuple<Ebsad_LoopKernelNxMType, Ebsad_LoopKernelNxMType> FuncPair;

FuncPair TEST_FUNC_PAIRS[] = {
    FuncPair(svt_sad_loop_kernel_c, svt_sad_loop_kernel_sse4_1_intrin),
    FuncPair(svt_sad_loop_kernel_c, svt_sad_loop_kernel_avx2_intrin),
#ifndef NON_AVX512_SUPPORT
    FuncPair(svt_sad_loop_kernel_c, svt_sad_loop_kernel_avx512_intrin),
#endif
};

FuncPair TEST_FUNC_PAIRS_SMALL[] = {
    FuncPair(svt_sad_loop_kernel_c, svt_sad_loop_kernel_sse4_1_intrin),
    FuncPair(svt_sad_loop_kernel_c, svt_sad_loop_kernel_avx2_intrin),
#ifndef NON_AVX512_SUPPORT
    FuncPair(svt_sad_loop_kernel_c, svt_sad_loop_kernel_avx512_intrin),
#endif
};

typedef std::tuple<TestPattern, BlkSize, SearchArea, FuncPair>
    sad_LoopTestParam;

/**
 * @brief Unit test for SAD loop (sparse, hme) functions include:
 *  - svt_sad_loop_kernel_{sse4_1,avx2,avx512}
 *  - svt_sad_loop_kernel_sparse_{sse4_1,avx2}_intrin
 *  - svt_sad_loop_kernel_{sse4_1,avx2}_hme_l0_intrin
 *
 * Test strategy:
 *  This test case combine different wight(4-64) x height(4-64), different test
 * vecotr pattern(MaxRef, MaxSrc, Random, Unalign) to generate test vector.
 * Run func with test vector, compare result between  non_avx2 function and avx2
 * function.
 *
 *
 * Expect result:
 *  Results come from  non_avx2 function and avx2 funtion are
 * equal.
 *
 * Test coverage:
 *
 * Test cases:
 *
 */
class sad_LoopTest : public ::testing::WithParamInterface<sad_LoopTestParam>,
                     public SADTestBase {
  public:
    sad_LoopTest()
        : SADTestBase(std::get<0>(TEST_GET_PARAM(1)),
                      std::get<1>(TEST_GET_PARAM(1)), TEST_GET_PARAM(0),
                      std::get<0>(TEST_GET_PARAM(2)),
                      std::get<1>(TEST_GET_PARAM(2))),
          func_c_(std::get<0>(TEST_GET_PARAM(3))),
          func_o_(std::get<1>(TEST_GET_PARAM(3))) {
    }

  protected:
    Ebsad_LoopKernelNxMType func_c_;
    Ebsad_LoopKernelNxMType func_o_;

    void check_sad_loop() {
        prepare_data();

        uint64_t best_sad0 = UINT64_MAX;
        int16_t x_search_center0 = 0;
        int16_t y_search_center0 = 0;
        func_c_(src_aligned_,
                src_stride_,
                ref1_aligned_,
                ref1_stride_,
                height_,
                width_,
                &best_sad0,
                &x_search_center0,
                &y_search_center0,
                ref1_stride_,
                search_area_width_,
                search_area_height_);

        uint64_t best_sad1 = UINT64_MAX;
        int16_t x_search_center1 = 0;
        int16_t y_search_center1 = 0;
        func_o_(src_aligned_,
                src_stride_,
                ref1_aligned_,
                ref1_stride_,
                height_,
                width_,
                &best_sad1,
                &x_search_center1,
                &y_search_center1,
                ref1_stride_,
                search_area_width_,
                search_area_height_);

        EXPECT_EQ(best_sad0, best_sad1)
            << "compare best_sad error"
            << " block dim: [" << width_ << " x " << height_ << "] "
            << "search area [" << search_area_width_ << " x "
            << search_area_height_ << "]";
        EXPECT_EQ(x_search_center0, x_search_center1)
            << "compare x_search_center error"
            << " block dim: [" << width_ << " x " << height_ << "] "
            << "search area [" << search_area_width_ << " x "
            << search_area_height_ << "]";
        EXPECT_EQ(y_search_center0, y_search_center1)
            << "compare y_search_center error"
            << " block dim: [" << width_ << " x " << height_ << "] "
            << "search area [" << search_area_width_ << " x "
            << search_area_height_ << "]";
    }

    void speed_sad_loop() {
        const uint64_t num_loop = 100000;
        double time_c, time_o;
        uint64_t start_time_seconds, start_time_useconds;
        uint64_t middle_time_seconds, middle_time_useconds;
        uint64_t finish_time_seconds, finish_time_useconds;

        prepare_data();

        uint64_t best_sad0 = UINT64_MAX;
        uint64_t best_sad1 = UINT64_MAX;
        int16_t x_search_center0 = 0;
        int16_t x_search_center1 = 0;
        int16_t y_search_center0 = 0;
        int16_t y_search_center1 = 0;

        svt_av1_get_time(&start_time_seconds, &start_time_useconds);

        for (uint64_t i = 0; i < num_loop; i++) {
            func_c_(src_aligned_,
                    src_stride_,
                    ref1_aligned_,
                    ref1_stride_,
                    height_,
                    width_,
                    &best_sad0,
                    &x_search_center0,
                    &y_search_center0,
                    ref1_stride_,
                    search_area_width_,
                    search_area_height_);
        }

        svt_av1_get_time(&middle_time_seconds, &middle_time_useconds);

        for (uint64_t i = 0; i < num_loop; i++) {
            func_o_(src_aligned_,
                    src_stride_,
                    ref1_aligned_,
                    ref1_stride_,
                    height_,
                    width_,
                    &best_sad1,
                    &x_search_center1,
                    &y_search_center1,
                    ref1_stride_,
                    search_area_width_,
                    search_area_height_);
        }

        svt_av1_get_time(&finish_time_seconds, &finish_time_useconds);

        EXPECT_EQ(best_sad0, best_sad1)
            << "compare best_sad error"
            << " block dim: [" << width_ << " x " << height_ << "] "
            << "search area [" << search_area_width_ << " x "
            << search_area_height_ << "]";
        EXPECT_EQ(x_search_center0, x_search_center1)
            << "compare x_search_center error"
            << " block dim: [" << width_ << " x " << height_ << "] "
            << "search area [" << search_area_width_ << " x "
            << search_area_height_ << "]";
        EXPECT_EQ(y_search_center0, y_search_center1)
            << "compare y_search_center error"
            << " block dim: [" << width_ << " x " << height_ << "] "
            << "search area [" << search_area_width_ << " x "
            << search_area_height_ << "]";

        time_c = svt_av1_compute_overall_elapsed_time_ms(start_time_seconds,
                                                         start_time_useconds,
                                                         middle_time_seconds,
                                                         middle_time_useconds);
        time_o = svt_av1_compute_overall_elapsed_time_ms(middle_time_seconds,
                                                         middle_time_useconds,
                                                         finish_time_seconds,
                                                         finish_time_useconds);

        printf("    svt_sad_loop_kernel(%dx%d) search area[%dx%d]: %5.2fx)\n",
               width_,
               height_,
               search_area_width_,
               search_area_height_,
               time_c / time_o);
    }
};

TEST_P(sad_LoopTest, sad_LoopTest) {
    check_sad_loop();
}

TEST_P(sad_LoopTest, DISABLED_sad_LoopSpeedTest) {
    speed_sad_loop();
}

INSTANTIATE_TEST_CASE_P(
    LOOPSAD, sad_LoopTest,
    ::testing::Combine(::testing::ValuesIn(TEST_PATTERNS),
                       ::testing::ValuesIn(TEST_BLOCK_SIZES),
                       ::testing::ValuesIn(TEST_LOOP_AREAS),
                       ::testing::ValuesIn(TEST_FUNC_PAIRS)));

INSTANTIATE_TEST_CASE_P(
    LOOPSAD_SMALL, sad_LoopTest,
    ::testing::Combine(::testing::ValuesIn(TEST_PATTERNS),
                       ::testing::ValuesIn(TEST_BLOCK_SIZES_SMALL),
                       ::testing::ValuesIn(TEST_LOOP_AREAS),
                       ::testing::ValuesIn(TEST_FUNC_PAIRS_SMALL)));

/**
 * best_sadmxn in GetEightsad_Test,Allsad_CalculationTest and
 * Extsad_CalculationTest must be less than 0x7FFFFFFF because signed comparison
 * is used in test functions, which is as follow:
 *   - get_eight_horizontal_search_point_results_8x8_16x16_pu_avx2_intrin
 *   - get_eight_horizontal_search_point_results_32x32_64x64_pu_avx2_intrin
 *   - svt_ext_all_sad_calculation_8x8_16x16_avx2
 *   - svt_ext_sad_calculation_8x8_16x16_avx2_intrin
 *   - svt_ext_sad_calculation_32x32_64x64_sse4_intrin
 */
#define BEST_SAD_MAX 0x7FFFFFFF

typedef void (*get_eight_sad_8_16_func)(uint8_t *src, uint32_t src_stride,
                                        uint8_t *ref, uint32_t ref_stride,
                                        uint32_t *p_best_sad_8x8,
                                        uint32_t *p_best_mv8x8,
                                        uint32_t *p_best_sad_16x16,
                                        uint32_t *p_best_mv16x16, uint32_t mv,
                                        uint16_t *p_sad16x16, EbBool sub_sad);

typedef void (*get_eight_sad_32_64_func)(uint16_t *p_sad16x16,
                                         uint32_t *p_best_sad_32x32,
                                         uint32_t *p_best_sad_64x64,
                                         uint32_t *p_best_mv32x32,
                                         uint32_t *p_best_mv64x64, uint32_t mv);
typedef std::tuple<TestPattern, SADPattern> sad_CalTestParam;

/**
 * @brief Unit test for Allsad_Calculation Test functions include:
 *  -
 * svt_ext_all_sad_calculation_8x8_16x16_avx2
 * svt_ext_eight_sad_calculation_32x32_64x64_avx2
 * svt_ext_eigth_sad_calculation_nsq_avx2
 *
 *
 * Test strategy:
 *  This test use different test pattern {REF_MAX, SRC_MAX, RANDOM, UNALIGN}
 *  to generate test vector,sad pattern {BUF_MAX, BUF_MIN, BUF_RANDOM} to
 *  generate test sad16x16.Check the result by compare result from reference
 *  function, non_avx2 function and avx2 function.
 *
 *
 * Expect result:
 *  Results come from  non_avx2 function and avx2 funtion are
 * equal.
 *
 * Test coverage:
 *
 * Test cases:
 **/

class Allsad_CalculationTest
    : public ::testing::WithParamInterface<sad_CalTestParam>,
      public SADTestBase {
  public:
    Allsad_CalculationTest()
        : SADTestBase(16, 16, TEST_GET_PARAM(0), TEST_GET_PARAM(1)) {
        src_stride_ = ref1_stride_ = MAX_SB_SIZE;
    }

  protected:
    void check_get_8x8_sad() {
        uint32_t best_sad8x8[2][64];
        uint32_t best_mv8x8[2][64] = {{0}};
        uint32_t best_sad16x16[2][16];
        uint32_t best_mv16x16[2][16] = {{0}};
        uint32_t eight_sad16x16[2][16][8];
        uint32_t eight_sad8x8[2][64][8];
        EbBool sub_sad = false;
        fill_buf_with_value(&best_sad8x8[0][0], 2 * 64, BEST_SAD_MAX);
        fill_buf_with_value(&best_sad16x16[0][0], 2 * 16, UINT_MAX);
        fill_buf_with_value(&eight_sad16x16[0][0][0], 2 * 16 * 8, UINT_MAX);
        fill_buf_with_value(&eight_sad8x8[0][0][0], 2 * 64 * 8, UINT_MAX);

        prepare_data();

        svt_ext_all_sad_calculation_8x8_16x16_c(src_aligned_,
                                                src_stride_,
                                                ref1_aligned_,
                                                ref1_stride_,
                                                0,
                                                best_sad8x8[0],
                                                best_sad16x16[0],
                                                best_mv8x8[0],
                                                best_mv16x16[0],
                                                eight_sad16x16[0],
                                                eight_sad8x8[0],
                                                sub_sad);

        svt_ext_all_sad_calculation_8x8_16x16_avx2(src_aligned_,
                                                   src_stride_,
                                                   ref1_aligned_,
                                                   ref1_stride_,
                                                   0,
                                                   best_sad8x8[1],
                                                   best_sad16x16[1],
                                                   best_mv8x8[1],
                                                   best_mv16x16[1],
                                                   eight_sad16x16[1],
                                                   eight_sad8x8[1],
                                                   sub_sad);

        EXPECT_EQ(
            0, memcmp(best_sad8x8[0], best_sad8x8[1], sizeof(best_sad8x8[0])))
            << "compare best_sad8x8 error sub_sad false";
        EXPECT_EQ(0,
                  memcmp(best_mv8x8[0], best_mv8x8[1], sizeof(best_mv8x8[0])))
            << "compare best_mv8x8 error sub_sad false";
        EXPECT_EQ(
            0,
            memcmp(
                best_sad16x16[0], best_sad16x16[1], sizeof(best_sad16x16[0])))
            << "compare best_sad16x16 error sub_sad false";
        EXPECT_EQ(
            0,
            memcmp(best_mv16x16[0], best_mv16x16[1], sizeof(best_mv16x16[0])))
            << "compare best_mv16x16 error sub_sad false";
        EXPECT_EQ(
            0,
            memcmp(eight_sad8x8[0], eight_sad8x8[1], sizeof(eight_sad8x8[0])))
            << "compare eight_sad8x8 error sub_sad false";
        EXPECT_EQ(0,
                  memcmp(eight_sad16x16[0],
                         eight_sad16x16[1],
                         sizeof(eight_sad16x16[0])))
            << "compare eight_sad16x16 error sub_sad false";

        sub_sad = true;
        fill_buf_with_value(&best_sad8x8[0][0], 2 * 64, BEST_SAD_MAX);
        fill_buf_with_value(&best_sad16x16[0][0], 2 * 16, UINT_MAX);
        fill_buf_with_value(&eight_sad16x16[0][0][0], 2 * 16 * 8, UINT_MAX);
        fill_buf_with_value(&eight_sad8x8[0][0][0], 2 * 64 * 8, UINT_MAX);

        prepare_data();

        svt_ext_all_sad_calculation_8x8_16x16_c(src_aligned_,
            src_stride_,
            ref1_aligned_,
            ref1_stride_,
            0,
            best_sad8x8[0],
            best_sad16x16[0],
            best_mv8x8[0],
            best_mv16x16[0],
            eight_sad16x16[0],
            eight_sad8x8[0],
            sub_sad);

        svt_ext_all_sad_calculation_8x8_16x16_avx2(src_aligned_,
            src_stride_,
            ref1_aligned_,
            ref1_stride_,
            0,
            best_sad8x8[1],
            best_sad16x16[1],
            best_mv8x8[1],
            best_mv16x16[1],
            eight_sad16x16[1],
            eight_sad8x8[1],
            sub_sad);

        EXPECT_EQ(
            0, memcmp(best_sad8x8[0], best_sad8x8[1], sizeof(best_sad8x8[0])))
            << "compare best_sad8x8 error sub_sad true";
        EXPECT_EQ(0,
            memcmp(best_mv8x8[0], best_mv8x8[1], sizeof(best_mv8x8[0])))
            << "compare best_mv8x8 error sub_sad true";
        EXPECT_EQ(
            0,
            memcmp(
                best_sad16x16[0], best_sad16x16[1], sizeof(best_sad16x16[0])))
            << "compare best_sad16x16 error sub_sad true";
        EXPECT_EQ(
            0,
            memcmp(best_mv16x16[0], best_mv16x16[1], sizeof(best_mv16x16[0])))
            << "compare best_mv16x16 error sub_sad true";
        EXPECT_EQ(
            0,
            memcmp(eight_sad8x8[0], eight_sad8x8[1], sizeof(eight_sad8x8[0])))
            << "compare eight_sad8x8 error sub_sad true";
        EXPECT_EQ(0,
            memcmp(eight_sad16x16[0],
                eight_sad16x16[1],
                sizeof(eight_sad16x16[0])))
            << "compare eight_sad16x16 error sub_sad true";
    }

    void check_get_32x32_sad() {
        uint32_t best_sad32x32[2][4];
        uint32_t best_sad64x64[2];
        uint32_t best_mv32x32[2][4] = {{0}};
        uint32_t best_mv64x64[2] = {0};
        uint32_t sad32x32[2][4][8];
        fill_buf_with_value(&best_sad32x32[0][0], 2 * 4, UINT_MAX);
        fill_buf_with_value(&best_sad64x64[0], 2, UINT_MAX);
        fill_buf_with_value(&sad32x32[0][0][0], 2 * 4 * 8, UINT_MAX);

        prepare_sad_data_32b();

        svt_ext_eight_sad_calculation_32x32_64x64_c(sad16x16_32b,
                                                    best_sad32x32[0],
                                                    &best_sad64x64[0],
                                                    best_mv32x32[0],
                                                    &best_mv64x64[0],
                                                    0,
                                                    sad32x32[0]);

        svt_ext_eight_sad_calculation_32x32_64x64_avx2(sad16x16_32b,
                                                       best_sad32x32[1],
                                                       &best_sad64x64[1],
                                                       best_mv32x32[1],
                                                       &best_mv64x64[1],
                                                       0,
                                                       sad32x32[1]);

        EXPECT_EQ(
            0,
            memcmp(
                best_sad32x32[0], best_sad32x32[1], sizeof(best_sad32x32[0])))
            << "compare best_sad32x32 error";
        EXPECT_EQ(
            0,
            memcmp(best_mv32x32[0], best_mv32x32[1], sizeof(best_mv32x32[0])))
            << "compare best_mv32x32 error";
        EXPECT_EQ(best_sad64x64[0], best_sad64x64[1])
            << "compare best_sad64x64 error";
        EXPECT_EQ(best_mv64x64[0], best_mv64x64[1])
            << "compare best_mv64x64 error";
        EXPECT_EQ(0, memcmp(sad32x32[0], sad32x32[1], sizeof(sad32x32[0])))
            << "compare sad32x32 error";
    }
};
TEST_P(Allsad_CalculationTest, 8x8_16x16_Test) {
    check_get_8x8_sad();
}

TEST_P(Allsad_CalculationTest, 32x32_64x64_Test) {
    check_get_32x32_sad();
}
INSTANTIATE_TEST_CASE_P(
    ALLSAD, Allsad_CalculationTest,
    ::testing::Combine(::testing::ValuesIn(TEST_PATTERNS),
                       ::testing::ValuesIn(TEST_SAD_PATTERNS)));
/**
 * @brief Unit test for Extsad_Calculation Test functions include:
 *  -
 * svt_ext_sad_calculation_8x8_16x16_avx2_intrin
 * svt_ext_sad_calculation_32x32_64x64_sse4_intrin
 *
 * Test strategy:
 *  This test use different test pattern {REF_MAX, SRC_MAX, RANDOM, UNALIGN}
 *  to generate test vector,sad pattern {BUF_MAX, BUF_MIN, BUF_RANDOM} to
 *  generate test sad16x16.Check the result by compare result from reference
 *  function, non_avx2 function and avx2 function.
 *
 *
 * Expect result:
 *  Results come from  non_avx2 function and avx2 funtion are
 * equal.
 *
 * Test coverage:
 *
 * Test cases:
 **/

class Extsad_CalculationTest
    : public ::testing::WithParamInterface<sad_CalTestParam>,
      public SADTestBase {
  public:
    Extsad_CalculationTest()
        : SADTestBase(16, 16, TEST_GET_PARAM(0), TEST_GET_PARAM(1)) {
        src_stride_ = ref1_stride_ = MAX_SB_SIZE;
    }

  protected:
    void check_get_8x8_sad() {
        uint32_t best_sad8x8[2][4];
        uint32_t best_mv8x8[2][4] = {{0}};
        uint32_t best_sad16x16[2], best_mv16x16[2] = {0};
        uint32_t sad16x16[2];
        uint32_t sad8x8[2][4];
        EbBool sub_sad = false;
        fill_buf_with_value(&best_sad8x8[0][0], 2 * 4, BEST_SAD_MAX);
        fill_buf_with_value(&best_sad16x16[0], 2, UINT_MAX);
        fill_buf_with_value(&sad16x16[0], 2, UINT_MAX);
        fill_buf_with_value(&sad8x8[0][0], 2 * 4, UINT_MAX);

        prepare_data();

        svt_ext_sad_calculation_8x8_16x16_c(src_aligned_,
                                            src_stride_,
                                            ref1_aligned_,
                                            ref1_stride_,
                                            best_sad8x8[0],
                                            &best_sad16x16[0],
                                            best_mv8x8[0],
                                            &best_mv16x16[0],
                                            0,
                                            &sad16x16[0],
                                            sad8x8[0],
                                            sub_sad);

        svt_ext_sad_calculation_8x8_16x16_avx2_intrin(src_aligned_,
                                                      src_stride_,
                                                      ref1_aligned_,
                                                      ref1_stride_,
                                                      best_sad8x8[1],
                                                      &best_sad16x16[1],
                                                      best_mv8x8[1],
                                                      &best_mv16x16[1],
                                                      0,
                                                      &sad16x16[1],
                                                      sad8x8[1],
                                                      sub_sad);

        EXPECT_EQ(
            0, memcmp(best_sad8x8[0], best_sad8x8[1], sizeof(best_sad8x8[0])))
            << "compare best_sad8x8 error sub_sad false";
        EXPECT_EQ(0,
                  memcmp(best_mv8x8[0], best_mv8x8[1], sizeof(best_mv8x8[0])))
            << "compare best_mv8x8 error sub_sad false";
        EXPECT_EQ(best_sad16x16[0], best_sad16x16[1])
            << "compare best_sad16x16 error sub_sad false";
        EXPECT_EQ(best_mv16x16[0], best_mv16x16[1])
            << "compare best_mv16x16 error sub_sad false";
        EXPECT_EQ(0, memcmp(sad8x8[0], sad8x8[1], sizeof(sad8x8[0])))
            << "compare sad8x8 error sub_sad false";
        EXPECT_EQ(sad16x16[0], sad16x16[1]) << "compare sad16x16 error sub_sad false";

        sub_sad = true;
        fill_buf_with_value(&best_sad8x8[0][0], 2 * 4, BEST_SAD_MAX);
        fill_buf_with_value(&best_sad16x16[0], 2, UINT_MAX);
        fill_buf_with_value(&sad16x16[0], 2, UINT_MAX);
        fill_buf_with_value(&sad8x8[0][0], 2 * 4, UINT_MAX);

        prepare_data();

        svt_ext_sad_calculation_8x8_16x16_c(src_aligned_,
            src_stride_,
            ref1_aligned_,
            ref1_stride_,
            best_sad8x8[0],
            &best_sad16x16[0],
            best_mv8x8[0],
            &best_mv16x16[0],
            0,
            &sad16x16[0],
            sad8x8[0],
            sub_sad);

        svt_ext_sad_calculation_8x8_16x16_avx2_intrin(src_aligned_,
            src_stride_,
            ref1_aligned_,
            ref1_stride_,
            best_sad8x8[1],
            &best_sad16x16[1],
            best_mv8x8[1],
            &best_mv16x16[1],
            0,
            &sad16x16[1],
            sad8x8[1],
            sub_sad);

        EXPECT_EQ(
            0, memcmp(best_sad8x8[0], best_sad8x8[1], sizeof(best_sad8x8[0])))
            << "compare best_sad8x8 error sub_sad true";
        EXPECT_EQ(0,
            memcmp(best_mv8x8[0], best_mv8x8[1], sizeof(best_mv8x8[0])))
            << "compare best_mv8x8 error sub_sad true";
        EXPECT_EQ(best_sad16x16[0], best_sad16x16[1])
            << "compare best_sad16x16 error sub_sad true";
        EXPECT_EQ(best_mv16x16[0], best_mv16x16[1])
            << "compare best_mv16x16 error sub_sad true";
        EXPECT_EQ(0, memcmp(sad8x8[0], sad8x8[1], sizeof(sad8x8[0])))
            << "compare sad8x8 error sub_sad true";
        EXPECT_EQ(sad16x16[0], sad16x16[1]) << "compare sad16x16 error sub_sad true";
    }

    void check_get_32x32_sad() {
        uint32_t best_sad32x32[2][4];
        uint32_t best_mv32x32[2][4] = {{0}};
        uint32_t best_sad64x64[2], best_mv64x64[2] = {0};
        uint32_t sad32x32[2][4];
        fill_buf_with_value(&best_sad32x32[0][0], 2 * 4, BEST_SAD_MAX);
        fill_buf_with_value(&best_sad64x64[0], 2, UINT_MAX);
        fill_buf_with_value(&sad32x32[0][0], 2 * 4, UINT_MAX);

        prepare_sad_data_32b();

        svt_ext_sad_calculation_32x32_64x64_c(*sad16x16_32b,
                                              best_sad32x32[0],
                                              &best_sad64x64[0],
                                              best_mv32x32[0],
                                              &best_mv64x64[0],
                                              0,
                                              sad32x32[0]);

        svt_ext_sad_calculation_32x32_64x64_sse4_intrin(*sad16x16_32b,
                                                        best_sad32x32[1],
                                                        &best_sad64x64[1],
                                                        best_mv32x32[1],
                                                        &best_mv64x64[1],
                                                        0,
                                                        sad32x32[1]);

        EXPECT_EQ(
            0,
            memcmp(
                best_sad32x32[0], best_sad32x32[1], sizeof(best_sad32x32[0])))
            << "compare best_sad32x32 error";
        EXPECT_EQ(
            0,
            memcmp(best_mv32x32[0], best_mv32x32[1], sizeof(best_mv32x32[0])))
            << "compare best_mv32x32 error";
        EXPECT_EQ(best_sad64x64[0], best_sad64x64[1])
            << "compare best_sad64x64 error";
        EXPECT_EQ(best_mv64x64[0], best_mv64x64[1])
            << "compare best_mv64x64 error";
        EXPECT_EQ(0, memcmp(sad32x32[0], sad32x32[1], sizeof(sad32x32[0])))
            << "compare sad32x32 error";
    }
};

TEST_P(Extsad_CalculationTest, Extsad_8x8Test) {
    check_get_8x8_sad();
}

TEST_P(Extsad_CalculationTest, Extsad_32x32Test) {
    check_get_32x32_sad();
}

INSTANTIATE_TEST_CASE_P(
    EXTSAD, Extsad_CalculationTest,
    ::testing::Combine(::testing::ValuesIn(TEST_PATTERNS),
                       ::testing::ValuesIn(TEST_SAD_PATTERNS)));
using InitializeBuffer_param_t = ::testing::tuple<uint32_t, uint32_t>;
#define MAX_BUFFER_SIZE 100  // const value to simplify
class InitializeBuffer32
    : public ::testing::TestWithParam<InitializeBuffer_param_t> {
  public:
    InitializeBuffer32()
        : count128(TEST_GET_PARAM(0)),
          count32(TEST_GET_PARAM(1)),
          rnd_(0, (1 << 30) - 1) {
        value = rnd_.random();
        _ref_ = (uint32_t *)svt_aom_memalign(32, MAX_BUFFER_SIZE);
        _test_ = (uint32_t *)svt_aom_memalign(32, MAX_BUFFER_SIZE);
        memset(_ref_, 0, MAX_BUFFER_SIZE);
        memset(_test_, 0, MAX_BUFFER_SIZE);
    }

    ~InitializeBuffer32() {
        if (_ref_)
            svt_aom_free(_ref_);
        if (_test_)
            svt_aom_free(_test_);
    }

  protected:
    void checkWithSize() {
        svt_initialize_buffer_32bits_c(_ref_, count128, count32, value);
        svt_initialize_buffer_32bits_sse2_intrin(_test_, count128, count32, value);

        int cmpResult = memcmp(_ref_, _test_, MAX_BUFFER_SIZE);
        EXPECT_EQ(cmpResult, 0);
    }

  private:
    uint32_t *_ref_;
    uint32_t *_test_;
    uint32_t count128;
    uint32_t count32;
    uint32_t value;
    SVTRandom rnd_;
};

TEST_P(InitializeBuffer32, InitializeBuffer) {
    checkWithSize();
}

INSTANTIATE_TEST_CASE_P(InitializeBuffer32, InitializeBuffer32,
                        ::testing::Combine(::testing::Values(2, 3, 4),
                                           ::testing::Values(1, 2, 3)));
/**
 * @Brief Base class for SAD test. SADTestBasesad_16Bit handle test vector in
 * memory, provide SAD and SAD avg reference function
 */
class SADTestBase16bit : public ::testing::Test {
  public:
    SADTestBase16bit(const int width, const int height,
                     TestPattern test_pattern) {
        width_ = width;
        height_ = height;
        src_stride_ = MAX_SB_SIZE;
        ref_stride_ = MAX_SB_SIZE / 2;
        test_pattern_ = test_pattern;
        src_ = nullptr;
        ref_ = nullptr;
    }

    void SetUp() override {
        src_ = (uint16_t *)svt_aom_memalign(32, MAX_BLOCK_SIZE * sizeof(*src_));
        ref_ = (uint16_t *)svt_aom_memalign(32, MAX_BLOCK_SIZE * sizeof(*ref_));
        ASSERT_NE(src_, nullptr);
        ASSERT_NE(ref_, nullptr);
    }

    void TearDown() override {
        if (src_)
            svt_aom_free(src_);
        if (ref_)
            svt_aom_free(ref_);
    }

    void prepare_data() {
        const int32_t mask = (1 << 16) - 1;
        SVTRandom rnd(0, mask);
        switch (test_pattern_) {
        case REF_MAX: {
            for (int i = 0; i < MAX_BLOCK_SIZE; i++)
                src_[i] = 0;

            for (int i = 0; i < MAX_BLOCK_SIZE; i++)
                ref_[i] = mask;

            break;
        }
        case SRC_MAX: {
            for (int i = 0; i < MAX_BLOCK_SIZE; i++)
                src_[i] = mask;

            for (int i = 0; i < MAX_SB_SIZE; i++)
                ref_[i] = 0;

            break;
        }
        case RANDOM: {
            for (int i = 0; i < MAX_BLOCK_SIZE; i++)
                src_[i] = rnd.random();

            for (int i = 0; i < MAX_BLOCK_SIZE; ++i) {
                ref_[i] = rnd.random();
            }
            break;
        };
        case UNALIGN: {
            for (int i = 0; i < MAX_BLOCK_SIZE; i++)
                src_[i] = rnd.random();

            for (int i = 0; i < MAX_BLOCK_SIZE; ++i) {
                ref_[i] = rnd.random();
            }
            ref_stride_ -= 1;
            break;
        }
        default: break;
        }
    }

  protected:
    uint32_t width_, height_;
    uint32_t src_stride_;
    uint32_t ref_stride_;
    TestPattern test_pattern_;
    SADPattern test_sad_pattern_;
    uint16_t *src_;
    uint16_t *ref_;
};

/**
 * @brief Unit test for SAD sub smaple functions include:
 *  - sad_16b_kernel_c
 *  - sad_16bit_kernel_avx2
 *
 * Test strategy:
 *  This test case combine different width{4-64} x height{4-64} and different
 * test pattern(REF_MAX, SRC_MAX, RANDOM, UNALIGN). Check the result by compare
 *  result from reference function, non_avx2 function and avx2 function.
 *
 *
 * Expect result:
 *  Results from reference functon, non_avx2 function and avx2 funtion are
 * equal.
 *
 * Test coverage:
 *  All functions inside sad_16b_kernel_c and
 * sad_16bit_kernel_avx2.
 *
 * Test cases:
 *  Width {4, 8, 16, 24, 32, 48, 64, 128} x height{ 4, 8, 16, 24, 32, 48, 64,
 * 128) Test vector pattern {REF_MAX, SRC_MAX, RANDOM, UNALIGN}
 *
 */
class SADTestSubSample16bit
    : public ::testing::WithParamInterface<Testsad_Param>,
      public SADTestBase16bit {
  public:
    SADTestSubSample16bit()
        : SADTestBase16bit(std::get<0>(TEST_GET_PARAM(1)),
                           std::get<1>(TEST_GET_PARAM(1)), TEST_GET_PARAM(0)) {
    }

  protected:
    void check_sad() {
        uint32_t repeat = 1;
        if (test_pattern_ == RANDOM) {
            repeat = 30;
        }

        for (uint32_t i = 0; i < repeat; ++i) {
            uint32_t sad_c = 0;
            uint32_t sad_avx2 = 0;

            prepare_data();

            sad_c = sad_16b_kernel_c(
                src_, src_stride_, ref_, ref_stride_, height_, width_);

            sad_avx2 = sad_16bit_kernel_avx2(
                src_, src_stride_, ref_, ref_stride_, height_, width_);

            EXPECT_EQ(sad_c, sad_avx2)
                << "compare sad_16b_kernel_c and sad_16bit_kernel_avx2 error, "
                   "repeat: "
                << i;
        }
    }

    void RunSpeedTest() {
        uint32_t sad_c = 0;
        uint32_t sad_avx2 = 0;

        double time_c, time_o;
        uint64_t start_time_seconds, start_time_useconds;
        uint64_t middle_time_seconds, middle_time_useconds;
        uint64_t finish_time_seconds, finish_time_useconds;

        prepare_data();

        for (uint32_t area_width = 4; area_width <= 128; area_width += 4) {
            const uint32_t area_height = area_width;
            const int num_loops = 1000000000 / (area_width * area_height);
            svt_av1_get_time(&start_time_seconds, &start_time_useconds);

            for (int i = 0; i < num_loops; ++i) {
                sad_c = sad_16b_kernel_c(
                    src_, src_stride_, ref_, ref_stride_, height_, width_);
            }

            svt_av1_get_time(&middle_time_seconds, &middle_time_useconds);

            for (int i = 0; i < num_loops; ++i) {
                sad_avx2 = sad_16bit_kernel_avx2(
                    src_, src_stride_, ref_, ref_stride_, height_, width_);
            }
            svt_av1_get_time(&finish_time_seconds, &finish_time_useconds);

            EXPECT_EQ(sad_c, sad_avx2) << area_width << "x" << area_height;

            time_c =
                svt_av1_compute_overall_elapsed_time_ms(start_time_seconds,
                                                        start_time_useconds,
                                                        middle_time_seconds,
                                                        middle_time_useconds);
            time_o =
                svt_av1_compute_overall_elapsed_time_ms(middle_time_seconds,
                                                        middle_time_useconds,
                                                        finish_time_seconds,
                                                        finish_time_useconds);
            printf("Average Nanoseconds per Function Call\n");
            printf("    sad_16b_kernel_c  (%dx%d) : %6.2f\n",
                   area_width,
                   area_height,
                   1000000 * time_c / num_loops);
            printf(
                "    sad_16bit_kernel_avx2(%dx%d) : %6.2f   "
                "(Comparison: %5.2fx)\n",
                area_width,
                area_height,
                1000000 * time_o / num_loops,
                time_c / time_o);
        }
    }
};

BlkSize TEST_BLOCK_SAD_SIZES[] = {
    BlkSize(16, 10), BlkSize(16, 5),
    BlkSize(32, 10), BlkSize(32, 20), BlkSize(64, 20),
    BlkSize(64, 64),  BlkSize(64, 32), BlkSize(32, 64),  BlkSize(32, 32),
    BlkSize(32, 16),  BlkSize(16, 32), BlkSize(16, 16),  BlkSize(16, 8),
    BlkSize(8, 16),   BlkSize(8, 8),   BlkSize(8, 4),    BlkSize(4, 4),
    BlkSize(4, 8),    BlkSize(4, 16),  BlkSize(16, 4),   BlkSize(8, 32),
    BlkSize(32, 8),   BlkSize(16, 64), BlkSize(16, 128), BlkSize(128, 128),
    BlkSize(64, 16),  BlkSize(24, 24), BlkSize(24, 16),  BlkSize(16, 24),
    BlkSize(24, 8),   BlkSize(8, 24),  BlkSize(64, 24),  BlkSize(48, 24),
    BlkSize(32, 24),  BlkSize(24, 32), BlkSize(48, 48),  BlkSize(48, 16),
    BlkSize(48, 32),  BlkSize(16, 48), BlkSize(32, 48),  BlkSize(48, 64),
    BlkSize(64, 48),  BlkSize(64, 48), BlkSize(128, 64), BlkSize(64, 128),
    BlkSize(128, 128)};

TEST_P(SADTestSubSample16bit, SADTestSubSample16bit) {
    check_sad();
}

INSTANTIATE_TEST_CASE_P(
    SAD, SADTestSubSample16bit,
    ::testing::Combine(::testing::ValuesIn(TEST_PATTERNS),
                       ::testing::ValuesIn(TEST_BLOCK_SAD_SIZES)));

TEST_P(SADTestSubSample16bit, DISABLED_Speed) {
    RunSpeedTest();
}

}  // namespace
