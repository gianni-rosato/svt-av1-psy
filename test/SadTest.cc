/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file SadTest.cc
 *
 * @brief Unit test for SAD functions:
 * - nxm_sad_kernel_sub_sampled_func
 * - nxm_sad_kernel_func
 * - nxm_sad_averaging_kernel_func
 * - nxm_sad_loop_kernel_sparse_func
 * - nxm_sad_loop_kernel_sparse_func
 * - get_eight_horizontal_search_point_results_8x8_16x16_func
 * - get_eight_horizontal_search_point_results_32x32_64x64_func
 * - Ext_ext_all_sad_calculation_8x8_16x16_func
 * - Ext_ext_eight_sad_calculation_32x32_64x64_func
 * - Ext_eigth_sad_calculation_nsq_func
 * - ExtSadCalculation_8x8_16x16_func
 * - ExtSadCalculation_32x32_64x64_func
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
#include "EbMeSadCalculation_SSE2.h"
#include "EbMotionEstimation.h"
#include "EbMotionEstimationContext.h"
#include "EbTime.h"
#include "random.h"
#include "util.h"

using svt_av1_test_tool::SVTRandom;  // to generate the random
extern "C" void ext_all_sad_calculation_8x8_16x16_c(
    uint8_t *src, uint32_t src_stride, uint8_t *ref, uint32_t ref_stride,
    uint32_t mv, uint32_t *p_best_sad8x8, uint32_t *p_best_sad16x16,
    uint32_t *p_best_mv8x8, uint32_t *p_best_mv16x16,
    uint32_t p_eight_sad16x16[16][8], uint32_t p_eight_sad8x8[64][8]);
extern "C" void ext_eigth_sad_calculation_nsq_c(
    uint32_t p_sad8x8[64][8], uint32_t p_sad16x16[16][8],
    uint32_t p_sad32x32[4][8], uint32_t *p_best_sad64x32,
    uint32_t *p_best_mv64x32, uint32_t *p_best_sad32x16,
    uint32_t *p_best_mv32x16, uint32_t *p_best_sad16x8, uint32_t *p_best_mv16x8,
    uint32_t *p_best_sad32x64, uint32_t *p_best_mv32x64,
    uint32_t *p_best_sad16x32, uint32_t *p_best_mv16x32,
    uint32_t *p_best_sad8x16, uint32_t *p_best_mv8x16, uint32_t *p_best_sad32x8,
    uint32_t *p_best_mv32x8, uint32_t *p_best_sad8x32, uint32_t *p_best_mv8x32,
    uint32_t *p_best_sad64x16, uint32_t *p_best_mv64x16,
    uint32_t *p_best_sad16x64, uint32_t *p_best_mv16x64, uint32_t mv);
extern "C" void ext_eight_sad_calculation_32x32_64x64_c(
    uint32_t p_sad16x16[16][8], uint32_t *p_best_sad32x32,
    uint32_t *p_best_sad64x64, uint32_t *p_best_mv32x32,
    uint32_t *p_best_mv64x64, uint32_t mv, uint32_t p_sad32x32[4][8]);

namespace {
/**
 * @Brief Test param definition.
 * - TEST_BLOCK_SIZES:All of Sad calculation funcs in test cover block size
 *  Width {4, 8, 16, 24, 32, 48, 64} x height{ 4, 8, 16, 24, 32, 48, 64).
 */
#define MAX_BLOCK_SIZE (MAX_SB_SIZE * MAX_SB_SIZE)
#define MAX_REF_BLOCK_SIZE \
    ((MAX_SEARCH_AREA_WIDTH_CH) * (MAX_SEARCH_AREA_HEIGHT_CH))
typedef std::tuple<int, int> BlkSize;
typedef enum { REF_MAX, SRC_MAX, RANDOM, UNALIGN } TestPattern;
typedef enum { BUF_MAX, BUF_MIN, BUF_SMALL, BUF_RANDOM } SADPattern;
BlkSize TEST_BLOCK_SIZES[] = {
    BlkSize(64, 64), BlkSize(64, 32), BlkSize(32, 64), BlkSize(32, 32),
    BlkSize(32, 16), BlkSize(16, 32), BlkSize(16, 16), BlkSize(16, 8),
    BlkSize(8, 16),  BlkSize(8, 8),   BlkSize(8, 4),   BlkSize(4, 4),
    BlkSize(4, 8),   BlkSize(4, 16),  BlkSize(16, 4),  BlkSize(8, 32),
    BlkSize(32, 8),  BlkSize(16, 64), BlkSize(64, 16), BlkSize(24, 24),
    BlkSize(24, 16), BlkSize(16, 24), BlkSize(24, 8),  BlkSize(8, 24),
    BlkSize(64, 24), BlkSize(48, 24), BlkSize(32, 24), BlkSize(24, 32),
    BlkSize(48, 48), BlkSize(48, 16), BlkSize(48, 32), BlkSize(16, 48),
    BlkSize(32, 48), BlkSize(48, 64), BlkSize(64, 48)};
TestPattern TEST_PATTERNS[] = {REF_MAX, SRC_MAX, RANDOM, UNALIGN};
SADPattern TEST_SAD_PATTERNS[] = {BUF_MAX, BUF_MIN, BUF_SMALL, BUF_RANDOM};
typedef std::tuple<TestPattern, BlkSize> TestSadParam;

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
        ref1_stride_ = ref2_stride_ = MAX_SEARCH_AREA_WIDTH_CH;
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
        ref1_stride_ = ref2_stride_ = MAX_SEARCH_AREA_WIDTH_CH;
        test_pattern_ = test_pattern;
        test_sad_pattern_ = test_sad_pattern;
    }

    SADTestBase(const int width, const int height, TestPattern test_pattern,
                const int search_area_width, const int search_area_height) {
        width_ = width;
        height_ = height;
        src_stride_ = MAX_SB_SIZE;
        ref1_stride_ = ref2_stride_ = MAX_SEARCH_AREA_WIDTH_CH;
        test_pattern_ = test_pattern;
        search_area_width_ = search_area_width;
        search_area_height_ = search_area_height;
    }

    void SetUp() override {
        src_aligned_ = (uint8_t *)eb_aom_memalign(32, MAX_BLOCK_SIZE);
        ref1_aligned_ = (uint8_t *)eb_aom_memalign(32, MAX_REF_BLOCK_SIZE);
        ref2_aligned_ = (uint8_t *)eb_aom_memalign(32, MAX_REF_BLOCK_SIZE);
        ASSERT_NE(src_aligned_, nullptr);
        ASSERT_NE(ref1_aligned_, nullptr);
        ASSERT_NE(ref2_aligned_, nullptr);
    }

    void TearDown() override {
        if (src_aligned_)
            eb_aom_free(src_aligned_);
        if (ref1_aligned_)
            eb_aom_free(ref1_aligned_);
        if (ref2_aligned_)
            eb_aom_free(ref2_aligned_);
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
 *  - nxm_sad_kernel_helper_c
 *  - nxm_sad_kernel_sub_sampled_helper_avx2
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
 *  All functions inside nxm_sad_kernel_helper_c and
 * nxm_sad_kernel_sub_sampled_helper_avx2.
 *
 * Test cases:
 *  Width {4, 8, 16, 24, 32, 48, 64} x height{ 4, 8, 16, 24, 32, 48, 64)
 *  Test vector pattern {REF_MAX, SRC_MAX, RANDOM, UNALIGN}
 *
 */
class SADTestSubSample : public ::testing::WithParamInterface<TestSadParam>,
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
        non_avx2_sad = nxm_sad_kernel_helper_c(
                                         src_aligned_,
                                         src_stride_,
                                         ref1_aligned_,
                                         ref1_stride_,
                                         height_,
                                         width_);

        avx2_sad = nxm_sad_kernel_sub_sampled_helper_avx2(
                                 src_aligned_,
                                 src_stride_,
                                 ref1_aligned_,
                                 ref1_stride_,
                                 height_,
                                 width_);

        EXPECT_EQ(non_avx2_sad, avx2_sad)
            << "compare non_avx2 and non_avx2 error";

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
 *  - nxm_sad_kernel_helper_c
 *  - nxm_sad_kernel_helper_avx2
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
 *  All functions inside nxm_sad_kernel_helper_c and
 *  nxm_sad_kernel_helper_avx2.
 *
 * Test cases:
 *  Width {4, 8, 16, 24, 32, 48, 64} x height{ 4, 8, 16, 24, 32, 48, 64)
 *  Test vector pattern {REF_MAX, SRC_MAX, RANDOM, UNALIGN}.
 *
 */
class SADTest : public ::testing::WithParamInterface<TestSadParam>,
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
       
        non_avx2_sad = nxm_sad_kernel_helper_c(src_aligned_,
                                        src_stride_,
                                        ref1_aligned_,
                                        ref1_stride_,
                                        height_,
                                        width_);
    
        avx2_sad = nxm_sad_kernel_helper_avx2(src_aligned_,
                                        src_stride_,
                                        ref1_aligned_,
                                        ref1_stride_,
                                        height_,
                                        width_);

        EXPECT_EQ(non_avx2_sad, avx2_sad)
                << "compare non_avx2 and non_avx2 error";
    }
};

TEST_P(SADTest, SADTest) {
    check_sad();
}

INSTANTIATE_TEST_CASE_P(
    SAD, SADTest,
    ::testing::Combine(::testing::ValuesIn(TEST_PATTERNS),
                       ::testing::ValuesIn(TEST_BLOCK_SIZES)));
/**
 * @brief Unit test for SAD Avg functions include:
 *  - nxm_sad_avg_kernel_helper_c
 *  - nxm_sad_avg_kernel_helper_avx2
 *
 * Test strategy:
 *  This test case combine different width{4-64} x height{4-64} and different
 * test pattern(REF_MAX, SRC_MAX, RANDOM, UNALIGN). Check the result by compare
 *  result from reference function, non_avx2 function and avx2 function.
 *
 * Expect result:
 *  Results come from reference functon, non_avx2 function and avx2 funtion are
 * equal.
 *
 * Test coverage:
 *  All functions inside nxm_sad_avg_kernel_helper_c and
 * nxm_sad_avg_kernel_helper_avx2.
 *
 * Test cases:
 *  Width {4, 8, 16, 24, 32, 48, 64} x height {4, 8, 16, 24, 32, 48, 64)
 *  Test vector pattern {REF_MAX, SRC_MAX, RANDOM, UNALIGN}.
 *
 */
class SADAvgTest : public ::testing::WithParamInterface<TestSadParam>,
                   public SADTestBase {
  public:
    SADAvgTest()
        : SADTestBase(std::get<0>(TEST_GET_PARAM(1)),
                      std::get<1>(TEST_GET_PARAM(1)), TEST_GET_PARAM(0)) {
    }

  protected:
    void check_sad_avg() {
        uint32_t ref_sad = 0;
        uint32_t non_avx2_sad = 0;
        uint32_t avx2_sad = 0;

        prepare_data();

        ref_sad = reference_sad_avg();
        non_avx2_sad = nxm_sad_avg_kernel_helper_c(src_aligned_,
                                        src_stride_,
                                        ref1_aligned_,
                                        ref1_stride_,
                                        ref2_aligned_,
                                        ref2_stride_,
                                        height_,
                                        width_);

        avx2_sad = nxm_sad_avg_kernel_helper_avx2(src_aligned_,
                                 src_stride_,
                                 ref1_aligned_,
                                 ref1_stride_,
                                 ref2_aligned_,
                                 ref2_stride_,
                                 height_,
                                 width_);

        EXPECT_EQ(non_avx2_sad, avx2_sad)
                << "compare non_avx2 and non_avx2 error";

    }
};

TEST_P(SADAvgTest, SADAvgTest) {
    check_sad_avg();
}

INSTANTIATE_TEST_CASE_P(
    SAD, SADAvgTest,
    ::testing::Combine(::testing::ValuesIn(TEST_PATTERNS),
                       ::testing::ValuesIn(TEST_BLOCK_SIZES)));

typedef std::tuple<int16_t, int16_t> SearchArea;

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

SearchArea TEST_LOOP_AREAS[] = {
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

typedef void (*EbSadLoopKernelNxMType)(
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

typedef std::tuple<EbSadLoopKernelNxMType, EbSadLoopKernelNxMType> FuncPair;

FuncPair TEST_FUNC_PAIRS[] = {
    FuncPair(sad_loop_kernel_c, sad_loop_kernel_sse4_1_intrin),
    FuncPair(sad_loop_kernel_c, sad_loop_kernel_avx2_intrin),
    FuncPair(sad_loop_kernel_sparse_c, sad_loop_kernel_sparse_sse4_1_intrin),
    FuncPair(sad_loop_kernel_sparse_c, sad_loop_kernel_sparse_avx2_intrin),
#ifndef NON_AVX512_SUPPORT
    FuncPair(sad_loop_kernel_c, sad_loop_kernel_avx512_intrin),
#endif
};

FuncPair TEST_HME_FUNC_PAIRS[] = {
    FuncPair(sad_loop_kernel_c, sad_loop_kernel_sse4_1_hme_l0_intrin),
    FuncPair(sad_loop_kernel_c, sad_loop_kernel_avx2_hme_l0_intrin)};

typedef std::tuple<TestPattern, BlkSize, SearchArea, FuncPair> SadLoopTestParam;

/**
 * @brief Unit test for SAD loop (sparse, hme) functions include:
 *  - sad_loop_kernel_{sse4_1,avx2,avx512}
 *  - sad_loop_kernel_sparse_{sse4_1,avx2}_intrin
 *  - sad_loop_kernel_{sse4_1,avx2}_hme_l0_intrin
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
class SadLoopTest : public ::testing::WithParamInterface<SadLoopTestParam>,
                    public SADTestBase {
  public:
    SadLoopTest()
        : SADTestBase(std::get<0>(TEST_GET_PARAM(1)),
                      std::get<1>(TEST_GET_PARAM(1)), TEST_GET_PARAM(0),
                      std::get<0>(TEST_GET_PARAM(2)),
                      std::get<1>(TEST_GET_PARAM(2))),
          func_c_(std::get<0>(TEST_GET_PARAM(3))),
          func_o_(std::get<1>(TEST_GET_PARAM(3))) {
    }

  protected:
    EbSadLoopKernelNxMType func_c_;
    EbSadLoopKernelNxMType func_o_;

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

        EbStartTime(&start_time_seconds, &start_time_useconds);

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

        EbStartTime(&middle_time_seconds, &middle_time_useconds);

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

        EbStartTime(&finish_time_seconds, &finish_time_useconds);

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

        EbComputeOverallElapsedTimeMs(start_time_seconds,
                                      start_time_useconds,
                                      middle_time_seconds,
                                      middle_time_useconds,
                                      &time_c);
        EbComputeOverallElapsedTimeMs(middle_time_seconds,
                                      middle_time_useconds,
                                      finish_time_seconds,
                                      finish_time_useconds,
                                      &time_o);

        printf("    sad_loop_kernel(%dx%d) search area[%dx%d]: %5.2fx)\n",
               width_,
               height_,
               search_area_width_,
               search_area_height_,
               time_c / time_o);
    }
};

TEST_P(SadLoopTest, SadLoopTest) {
    check_sad_loop();
}

TEST_P(SadLoopTest, DISABLED_SadLoopSpeedTest) {
    speed_sad_loop();
}

INSTANTIATE_TEST_CASE_P(
    LOOPSAD, SadLoopTest,
    ::testing::Combine(::testing::ValuesIn(TEST_PATTERNS),
                       ::testing::ValuesIn(TEST_BLOCK_SIZES),
                       ::testing::ValuesIn(TEST_LOOP_AREAS),
                       ::testing::ValuesIn(TEST_FUNC_PAIRS)));

INSTANTIATE_TEST_CASE_P(
    HMESAD, SadLoopTest,
    ::testing::Combine(::testing::ValuesIn(TEST_PATTERNS),
                       ::testing::ValuesIn(TEST_BLOCK_SIZES),
                       ::testing::ValuesIn(TEST_AREAS),
                       ::testing::ValuesIn(TEST_HME_FUNC_PAIRS)));

/**
 * best_sadmxn in GetEightSadTest,AllSadCalculationTest and
 * ExtSadCalculationTest must be less than 0x7FFFFFFF because signed comparison
 * is used in test functions, which is as follow:
 *   - get_eight_horizontal_search_point_results_8x8_16x16_pu_avx2_intrin
 *   - get_eight_horizontal_search_point_results_32x32_64x64_pu_avx2_intrin
 *   - ext_all_sad_calculation_8x8_16x16_avx2
 *   - ext_sad_calculation_8x8_16x16_avx2_intrin
 *   - ext_sad_calculation_32x32_64x64_sse4_intrin
 */
#define BEST_SAD_MAX 0x7FFFFFFF

typedef void (*get_eight_sad_8_16_func)(uint8_t *src, uint32_t src_stride,
                                        uint8_t *ref, uint32_t ref_stride,
                                        uint32_t *p_best_sad8x8,
                                        uint32_t *p_best_mv8x8,
                                        uint32_t *p_best_sad16x16,
                                        uint32_t *p_best_mv16x16, uint32_t mv,
                                        uint16_t *p_sad16x16, EbBool sub_sad);

typedef void (*get_eight_sad_32_64_func)(uint16_t *p_sad16x16,
                                         uint32_t *p_best_sad32x32,
                                         uint32_t *p_best_sad64x64,
                                         uint32_t *p_best_mv32x32,
                                         uint32_t *p_best_mv64x64, uint32_t mv);

static const get_eight_sad_8_16_func get_eight_sad_8_16_func_table[] = {
    get_eight_horizontal_search_point_results_8x8_16x16_pu_sse41_intrin,
    get_eight_horizontal_search_point_results_8x8_16x16_pu_avx2_intrin};

static const get_eight_sad_32_64_func get_eight_sad_32_64_func_table[] = {
    get_eight_horizontal_search_point_results_32x32_64x64_pu_sse41_intrin,
    get_eight_horizontal_search_point_results_32x32_64x64_pu_avx2_intrin};

typedef std::tuple<TestPattern, SADPattern> SadCalTestParam;

/**
 * @brief Unit test for GetEightSadTest functions include:
 *
 *  - get_eight_horizontal_search_point_results_8x8_16x16_pu_{sse41,avx2}_intrin
 *  -
 * get_eight_horizontal_search_point_results_32x32_64x64_pu_{sse41,avx2}_intrin
 *
 * Test strategy:
 *  This test use different test pattern {REF_MAX, SRC_MAX, RANDOM, UNALIGN}
 *  to generate test vector,sad pattern {BUF_MAX, BUF_MIN, BUF_RANDOM} to
 * generate test sad16x16.Check the result by compare result from reference
 *  function, non_avx2 function and avx2 function.
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
class GetEightSadTest : public ::testing::WithParamInterface<SadCalTestParam>,
                        public SADTestBase {
  public:
    GetEightSadTest()
        : SADTestBase(16, 16, TEST_GET_PARAM(0), TEST_GET_PARAM(1)) {
        src_stride_ = ref1_stride_ = MAX_SB_SIZE;
    }

  protected:
    void check_get_eight_8_16() {
        for (int i = 0; i < 10; i++) {
            uint32_t best_sad8x8_1[4] = {BEST_SAD_MAX, 0, BEST_SAD_MAX, 0};
            uint32_t best_mv8x8_1[4] = {
                0x00830147, 0x0093FFD4, 0xFF371257, 0xF082F7DA};
            uint32_t best_sad16x16_1 = UINT_MAX, best_mv16x16_1 = 0x00ACFFBD;
            uint16_t sad16x16_1[8] = {0, 1, 2, 3, 4, 5, 6, 7};

            prepare_data();

            get_eight_horizontal_search_point_results_8x8_16x16_pu_c(
                src_aligned_,
                src_stride_,
                ref1_aligned_,
                ref1_stride_,
                best_sad8x8_1,
                best_mv8x8_1,
                &best_sad16x16_1,
                &best_mv16x16_1,
                0,
                sad16x16_1,
                false);

            for (int j = 0; j < sizeof(get_eight_sad_8_16_func_table) /
                                    sizeof(*get_eight_sad_8_16_func_table);
                 j++) {
                uint32_t best_sad8x8_2[4] = {BEST_SAD_MAX, 0, BEST_SAD_MAX, 0};
                uint32_t best_mv8x8_2[4] = {
                    0x00830147, 0x0093FFD4, 0xFF371257, 0xF082F7DA};
                uint32_t best_sad16x16_2 = UINT_MAX,
                         best_mv16x16_2 = 0x00ACFFBD;
                uint16_t sad16x16_2[8] = {8, 9, 10, 11, 12, 13, 14, 15};
                get_eight_sad_8_16_func_table[j](src_aligned_,
                                                 src_stride_,
                                                 ref1_aligned_,
                                                 ref1_stride_,
                                                 best_sad8x8_2,
                                                 best_mv8x8_2,
                                                 &best_sad16x16_2,
                                                 &best_mv16x16_2,
                                                 0,
                                                 sad16x16_2,
                                                 false);

                EXPECT_EQ(
                    0,
                    memcmp(best_sad8x8_1, best_sad8x8_2, sizeof(best_sad8x8_1)))
                    << "compare best_sad8x8 error";
                EXPECT_EQ(
                    0, memcmp(best_mv8x8_1, best_mv8x8_2, sizeof(best_mv8x8_1)))
                    << "compare best_mv8x8 error";
                EXPECT_EQ(best_sad16x16_1, best_sad16x16_2)
                    << "compare best_sad16x16 error";
                EXPECT_EQ(best_mv16x16_1, best_mv16x16_2)
                    << "compare best_mv16x16 error";
                EXPECT_EQ(0, memcmp(sad16x16_1, sad16x16_2, sizeof(sad16x16_1)))
                    << "compare sad16x16 error";
            }
        }
    }

    void speed_get_eight_8_16() {
        uint32_t best_sad8x8_1[4] = {BEST_SAD_MAX, 0, BEST_SAD_MAX, 0};
        uint32_t best_mv8x8_1[4] = {
            0x00830147, 0x0093FFD4, 0xFF371257, 0xF082F7DA};
        uint32_t best_sad16x16_1 = UINT_MAX, best_mv16x16_1 = 0x00ACFFBD;
        uint16_t sad16x16_1[8] = {0, 1, 2, 3, 4, 5, 6, 7};
        const uint64_t num_loop = 1000000;
        double time_c, time_o;
        uint64_t start_time_seconds, start_time_useconds;
        uint64_t middle_time_seconds, middle_time_useconds;
        uint64_t finish_time_seconds, finish_time_useconds;

        prepare_data();

        EbStartTime(&start_time_seconds, &start_time_useconds);

        for (uint64_t i = 0; i < num_loop; i++) {
            get_eight_horizontal_search_point_results_8x8_16x16_pu_c(
                src_aligned_,
                src_stride_,
                ref1_aligned_,
                ref1_stride_,
                best_sad8x8_1,
                best_mv8x8_1,
                &best_sad16x16_1,
                &best_mv16x16_1,
                0,
                sad16x16_1,
                false);
        }

        EbStartTime(&middle_time_seconds, &middle_time_useconds);
        EbComputeOverallElapsedTimeMs(start_time_seconds,
                                      start_time_useconds,
                                      middle_time_seconds,
                                      middle_time_useconds,
                                      &time_c);

        for (int i = 0; i < sizeof(get_eight_sad_8_16_func_table) /
                                sizeof(*get_eight_sad_8_16_func_table);
             i++) {
            uint32_t best_sad8x8_2[4] = {BEST_SAD_MAX, 0, BEST_SAD_MAX, 0};
            uint32_t best_mv8x8_2[4] = {
                0x00830147, 0x0093FFD4, 0xFF371257, 0xF082F7DA};
            uint32_t best_sad16x16_2 = UINT_MAX, best_mv16x16_2 = 0x00ACFFBD;
            uint16_t sad16x16_2[8] = {8, 9, 10, 11, 12, 13, 14, 15};

            EbStartTime(&middle_time_seconds, &middle_time_useconds);

            for (uint64_t j = 0; j < num_loop; j++) {
                get_eight_sad_8_16_func_table[i](src_aligned_,
                                                 src_stride_,
                                                 ref1_aligned_,
                                                 ref1_stride_,
                                                 best_sad8x8_2,
                                                 best_mv8x8_2,
                                                 &best_sad16x16_2,
                                                 &best_mv16x16_2,
                                                 0,
                                                 sad16x16_2,
                                                 false);
            }

            EbStartTime(&finish_time_seconds, &finish_time_useconds);

            EXPECT_EQ(
                0, memcmp(best_sad8x8_1, best_sad8x8_2, sizeof(best_sad8x8_1)))
                << "compare best_sad8x8 error";
            EXPECT_EQ(0,
                      memcmp(best_mv8x8_1, best_mv8x8_2, sizeof(best_mv8x8_1)))
                << "compare best_mv8x8 error";
            EXPECT_EQ(best_sad16x16_1, best_sad16x16_2)
                << "compare best_sad16x16 error";
            EXPECT_EQ(best_mv16x16_1, best_mv16x16_2)
                << "compare best_mv16x16 error";
            EXPECT_EQ(0, memcmp(sad16x16_1, sad16x16_2, sizeof(sad16x16_1)))
                << "compare sad16x16 error";

            EbComputeOverallElapsedTimeMs(middle_time_seconds,
                                          middle_time_useconds,
                                          finish_time_seconds,
                                          finish_time_useconds,
                                          &time_o);

            printf(
                "get_eight_horizontal_search_point_results_8x8_16x16_pu(%d): "
                "%5.2fx)\n",
                i,
                time_c / time_o);
        }
    }

    void check_get_eight_32_64() {
        for (int i = 0; i < 10000; i++) {
            uint32_t best_sad32x32[4], best_sad32x32_1[4];
            uint32_t best_mv32x32_1[4] = {
                0x00010002, 0x0003FFF4, 0xFFF70008, 0xFFF9FFF1};
            uint32_t best_sad64x64_1 = UINT_MAX, best_mv64x64_1 = 0x0078FF94;
            const uint32_t mv = 0xFF8300DB;

            prepare_sad_data_16b(best_sad32x32);
            memcpy(best_sad32x32_1, best_sad32x32, sizeof(best_sad32x32));

            get_eight_horizontal_search_point_results_32x32_64x64_pu_c(
                *sad16x16_16b,
                best_sad32x32_1,
                &best_sad64x64_1,
                best_mv32x32_1,
                &best_mv64x64_1,
                mv);

            for (int j = 0; j < sizeof(get_eight_sad_32_64_func_table) /
                                    sizeof(*get_eight_sad_32_64_func_table);
                 j++) {
                uint32_t best_sad32x32_2[4];
                uint32_t best_mv32x32_2[4] = {
                    0x00010002, 0x0003FFF4, 0xFFF70008, 0xFFF9FFF1};
                uint32_t best_sad64x64_2 = UINT_MAX,
                         best_mv64x64_2 = 0x0078FF94;

                memcpy(best_sad32x32_2, best_sad32x32, sizeof(best_sad32x32));
                get_eight_sad_32_64_func_table[j](*sad16x16_16b,
                                                  best_sad32x32_2,
                                                  &best_sad64x64_2,
                                                  best_mv32x32_2,
                                                  &best_mv64x64_2,
                                                  mv);

                EXPECT_EQ(0,
                          memcmp(best_sad32x32_1,
                                 best_sad32x32_2,
                                 sizeof(best_sad32x32_1)))
                    << "compare best_sad32x32 error";
                EXPECT_EQ(best_mv32x32_1[0], best_mv32x32_2[0])
                    << "compare best_mv32x32[0] error";
                EXPECT_EQ(best_mv32x32_1[1], best_mv32x32_2[1])
                    << "compare best_mv32x32[1] error";
                EXPECT_EQ(best_mv32x32_1[2], best_mv32x32_2[2])
                    << "compare best_mv32x32[2] error";
                EXPECT_EQ(best_mv32x32_1[3], best_mv32x32_2[3])
                    << "compare best_mv32x32[3] error";
                EXPECT_EQ(best_sad64x64_1, best_sad64x64_2)
                    << "compare best_sad64x64 error";
                EXPECT_EQ(best_mv64x64_1, best_mv64x64_2)
                    << "compare best_mv64x64 error";
            }
        }
    }

    void speed_get_eight_32_64() {
        uint32_t best_sad32x32_1[4] = {BEST_SAD_MAX, 0, BEST_SAD_MAX, 0};
        uint32_t best_sad32x32_2[4] = {BEST_SAD_MAX, 0, BEST_SAD_MAX, 0};
        uint32_t best_mv32x32_1[4] = {
            0x00010002, 0x0003FFF4, 0xFFF70008, 0xFFF9FFF1};
        uint32_t best_sad64x64_1 = UINT_MAX, best_mv64x64_1 = 0x0078FF94;
        const uint64_t num_loop = 100000000;
        double time_c, time_o;
        uint64_t start_time_seconds, start_time_useconds;
        uint64_t middle_time_seconds, middle_time_useconds;
        uint64_t finish_time_seconds, finish_time_useconds;

        prepare_sad_data_16b(best_sad32x32_1);
        memcpy(best_sad32x32_2, best_sad32x32_1, sizeof(best_sad32x32_1));

        EbStartTime(&start_time_seconds, &start_time_useconds);

        for (uint64_t i = 0; i < num_loop; i++) {
            best_sad64x64_1 = UINT_MAX;
            get_eight_horizontal_search_point_results_32x32_64x64_pu_c(
                *sad16x16_16b,
                best_sad32x32_1,
                &best_sad64x64_1,
                best_mv32x32_1,
                &best_mv64x64_1,
                0);
        }

        EbStartTime(&middle_time_seconds, &middle_time_useconds);
        EbComputeOverallElapsedTimeMs(start_time_seconds,
                                      start_time_useconds,
                                      middle_time_seconds,
                                      middle_time_useconds,
                                      &time_c);

        for (int i = 0; i < sizeof(get_eight_sad_32_64_func_table) /
                                sizeof(*get_eight_sad_32_64_func_table);
             i++) {
            uint32_t best_mv32x32_2[4] = {
                0x00010002, 0x0003FFF4, 0xFFF70008, 0xFFF9FFF1};
            uint32_t best_sad64x64_2 = UINT_MAX, best_mv64x64_2 = 0x0078FF94;

            EbStartTime(&middle_time_seconds, &middle_time_useconds);

            for (uint64_t j = 0; j < num_loop; j++) {
                best_sad64x64_2 = UINT_MAX;
                get_eight_sad_32_64_func_table[i](*sad16x16_16b,
                                                  best_sad32x32_2,
                                                  &best_sad64x64_2,
                                                  best_mv32x32_2,
                                                  &best_mv64x64_2,
                                                  0);
            }

            EbStartTime(&finish_time_seconds, &finish_time_useconds);

            EXPECT_EQ(
                0,
                memcmp(
                    best_sad32x32_1, best_sad32x32_2, sizeof(best_sad32x32_1)))
                << "compare best_sad32x32 error";
            EXPECT_EQ(
                0,
                memcmp(best_mv32x32_1, best_mv32x32_2, sizeof(best_mv32x32_1)))
                << "compare best_mv32x32 error";
            EXPECT_EQ(best_sad64x64_1, best_sad64x64_2)
                << "compare best_sad64x64 error";
            EXPECT_EQ(best_mv64x64_1, best_mv64x64_2)
                << "compare best_mv64x64 error";

            EbComputeOverallElapsedTimeMs(middle_time_seconds,
                                          middle_time_useconds,
                                          finish_time_seconds,
                                          finish_time_useconds,
                                          &time_o);

            printf(
                "get_eight_horizontal_search_point_results_32x32_64x64_pu(%d): "
                "%5.2fx)\n",
                i,
                time_c / time_o);
        }
    }
};

TEST_P(GetEightSadTest, EightSadTest_8_16) {
    check_get_eight_8_16();
}

TEST_P(GetEightSadTest, DISABLED_EightSadSpeedTest_8_16) {
    speed_get_eight_8_16();
}

TEST_P(GetEightSadTest, EightSadTest_32_64) {
    check_get_eight_32_64();
}

TEST_P(GetEightSadTest, DISABLED_EightSadSpeedTest_32_64) {
    speed_get_eight_32_64();
}

INSTANTIATE_TEST_CASE_P(
    EIGHTSAD, GetEightSadTest,
    ::testing::Combine(::testing::ValuesIn(TEST_PATTERNS),
                       ::testing::ValuesIn(TEST_SAD_PATTERNS)));

/**
 * @brief Unit test for AllSadCalculation Test functions include:
 *  -
 * ext_all_sad_calculation_8x8_16x16_avx2
 * ext_eight_sad_calculation_32x32_64x64_avx2
 * ext_eigth_sad_calculation_nsq_avx2
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

class AllSadCalculationTest
    : public ::testing::WithParamInterface<SadCalTestParam>,
      public SADTestBase {
  public:
    AllSadCalculationTest()
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

        fill_buf_with_value(&best_sad8x8[0][0], 2 * 64, BEST_SAD_MAX);
        fill_buf_with_value(&best_sad16x16[0][0], 2 * 16, UINT_MAX);
        fill_buf_with_value(&eight_sad16x16[0][0][0], 2 * 16 * 8, UINT_MAX);
        fill_buf_with_value(&eight_sad8x8[0][0][0], 2 * 64 * 8, UINT_MAX);

        prepare_data();

        ext_all_sad_calculation_8x8_16x16_c(src_aligned_,
                                            src_stride_,
                                            ref1_aligned_,
                                            ref1_stride_,
                                            0,
                                            best_sad8x8[0],
                                            best_sad16x16[0],
                                            best_mv8x8[0],
                                            best_mv16x16[0],
                                            eight_sad16x16[0],
                                            eight_sad8x8[0]);

        ext_all_sad_calculation_8x8_16x16_avx2(src_aligned_,
                                               src_stride_,
                                               ref1_aligned_,
                                               ref1_stride_,
                                               0,
                                               best_sad8x8[1],
                                               best_sad16x16[1],
                                               best_mv8x8[1],
                                               best_mv16x16[1],
                                               eight_sad16x16[1],
                                               eight_sad8x8[1]);

        EXPECT_EQ(
            0, memcmp(best_sad8x8[0], best_sad8x8[1], sizeof(best_sad8x8[0])))
            << "compare best_sad8x8 error";
        EXPECT_EQ(0,
                  memcmp(best_mv8x8[0], best_mv8x8[1], sizeof(best_mv8x8[0])))
            << "compare best_mv8x8 error";
        EXPECT_EQ(
            0,
            memcmp(
                best_sad16x16[0], best_sad16x16[1], sizeof(best_sad16x16[0])))
            << "compare best_sad16x16 error";
        EXPECT_EQ(
            0,
            memcmp(best_mv16x16[0], best_mv16x16[1], sizeof(best_mv16x16[0])))
            << "compare best_mv16x16 error";
        EXPECT_EQ(
            0,
            memcmp(eight_sad8x8[0], eight_sad8x8[1], sizeof(eight_sad8x8[0])))
            << "compare eight_sad8x8 error";
        EXPECT_EQ(0,
                  memcmp(eight_sad16x16[0],
                         eight_sad16x16[1],
                         sizeof(eight_sad16x16[0])))
            << "compare eight_sad16x16 error";
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

        ext_eight_sad_calculation_32x32_64x64_c(sad16x16_32b,
                                                best_sad32x32[0],
                                                &best_sad64x64[0],
                                                best_mv32x32[0],
                                                &best_mv64x64[0],
                                                0,
                                                sad32x32[0]);

        ext_eight_sad_calculation_32x32_64x64_avx2(sad16x16_32b,
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

    void check_get_nsq_sad() {
        uint32_t best_sad64x32[2][2];
        uint32_t best_sad32x64[2][2];
        uint32_t best_sad32x16[2][8];
        uint32_t best_sad16x32[2][8];
        uint32_t best_sad16x8[2][32];
        uint32_t best_sad8x16[2][32];
        uint32_t best_sad32x8[2][16];
        uint32_t best_sad8x32[2][16];
        uint32_t best_sad64x16[2][4];
        uint32_t best_sad16x64[2][4];

        uint32_t best_mv64x32[2][2] = {{0}};
        uint32_t best_mv32x64[2][2] = {{0}};
        uint32_t best_mv32x16[2][8] = {{0}};
        uint32_t best_mv16x32[2][8] = {{0}};
        uint32_t best_mv16x8[2][32] = {{0}};
        uint32_t best_mv8x16[2][32] = {{0}};
        uint32_t best_mv32x8[2][16] = {{0}};
        uint32_t best_mv8x32[2][16] = {{0}};
        uint32_t best_mv16x64[2][4] = {{0}};
        uint32_t best_mv64x16[2][4] = {{0}};

        fill_buf_with_value(&best_sad64x32[0][0], 2 * 2, UINT_MAX);
        fill_buf_with_value(&best_sad32x64[0][0], 2 * 2, UINT_MAX);
        fill_buf_with_value(&best_sad32x16[0][0], 2 * 8, UINT_MAX);
        fill_buf_with_value(&best_sad16x32[0][0], 2 * 8, UINT_MAX);
        fill_buf_with_value(&best_sad16x8[0][0], 2 * 32, UINT_MAX);
        fill_buf_with_value(&best_sad8x16[0][0], 2 * 32, UINT_MAX);
        fill_buf_with_value(&best_sad32x8[0][0], 2 * 16, UINT_MAX);
        fill_buf_with_value(&best_sad8x32[0][0], 2 * 16, UINT_MAX);
        fill_buf_with_value(&best_sad64x16[0][0], 2 * 4, UINT_MAX);
        fill_buf_with_value(&best_sad16x64[0][0], 2 * 4, UINT_MAX);

        prepare_nsq_sad_data();

        ext_eigth_sad_calculation_nsq_c(sad8x8,
                                        sad16x16_32b,
                                        sad32x32,
                                        best_sad64x32[0],
                                        best_mv64x32[0],
                                        best_sad32x16[0],
                                        best_mv32x16[0],
                                        best_sad16x8[0],
                                        best_mv16x8[0],
                                        best_sad32x64[0],
                                        best_mv32x64[0],
                                        best_sad16x32[0],
                                        best_mv16x32[0],
                                        best_sad8x16[0],
                                        best_mv8x16[0],
                                        best_sad32x8[0],
                                        best_mv32x8[0],
                                        best_sad8x32[0],
                                        best_mv8x32[0],
                                        best_sad64x16[0],
                                        best_mv64x16[0],
                                        best_sad16x64[0],
                                        best_mv16x64[0],
                                        0);

        ext_eigth_sad_calculation_nsq_avx2(sad8x8,
                                           sad16x16_32b,
                                           sad32x32,
                                           best_sad64x32[1],
                                           best_mv64x32[1],
                                           best_sad32x16[1],
                                           best_mv32x16[1],
                                           best_sad16x8[1],
                                           best_mv16x8[1],
                                           best_sad32x64[1],
                                           best_mv32x64[1],
                                           best_sad16x32[1],
                                           best_mv16x32[1],
                                           best_sad8x16[1],
                                           best_mv8x16[1],
                                           best_sad32x8[1],
                                           best_mv32x8[1],
                                           best_sad8x32[1],
                                           best_mv8x32[1],
                                           best_sad64x16[1],
                                           best_mv64x16[1],
                                           best_sad16x64[1],
                                           best_mv16x64[1],
                                           0);

        EXPECT_EQ(
            0,
            memcmp(
                best_sad64x32[0], best_sad64x32[1], sizeof(best_sad64x32[0])))
            << "compare best_sad64x32 error";
        EXPECT_EQ(
            0,
            memcmp(best_mv64x32[0], best_mv64x32[1], sizeof(best_mv64x32[0])))
            << "compare best_mv64x32 error";
        EXPECT_EQ(
            0,
            memcmp(
                best_sad32x16[0], best_sad32x16[1], sizeof(best_sad32x16[0])))
            << "compare best_sad32x16 error";
        EXPECT_EQ(
            0,
            memcmp(best_mv32x16[0], best_mv32x16[1], sizeof(best_mv32x16[0])))
            << "compare best_mv32x16 error";
        EXPECT_EQ(
            0,
            memcmp(best_sad16x8[0], best_sad16x8[1], sizeof(best_sad16x8[0])))
            << "compare best_sad16x8 error";
        EXPECT_EQ(
            0, memcmp(best_mv16x8[0], best_mv16x8[1], sizeof(best_mv16x8[0])))
            << "compare best_mv16x8 error";
        EXPECT_EQ(
            0,
            memcmp(
                best_sad32x64[0], best_sad32x64[1], sizeof(best_sad32x64[0])))
            << "compare best_sad32x64 error";
        EXPECT_EQ(
            0,
            memcmp(best_mv32x64[0], best_mv32x64[1], sizeof(best_mv32x64[0])))
            << "compare best_mv32x64 error";
        EXPECT_EQ(
            0,
            memcmp(best_mv16x32[0], best_mv16x32[1], sizeof(best_mv16x32[0])))
            << "compare best_mv16x32 error";
        EXPECT_EQ(
            0,
            memcmp(
                best_sad16x32[0], best_sad16x32[1], sizeof(best_sad16x32[0])))
            << "compare best_sad16x32 error";
        EXPECT_EQ(
            0,
            memcmp(best_sad8x16[0], best_sad8x16[1], sizeof(best_sad8x16[0])))
            << "compare best_sad8x16 error";
        EXPECT_EQ(
            0, memcmp(best_mv8x16[0], best_mv8x16[1], sizeof(best_mv8x16[0])))
            << "compare best_mv8x16 error";
        EXPECT_EQ(
            0,
            memcmp(best_sad32x8[0], best_sad32x8[1], sizeof(best_sad32x8[0])))
            << "compare best_sad32x8 error";
        EXPECT_EQ(
            0, memcmp(best_mv32x8[0], best_mv32x8[1], sizeof(best_mv32x8[0])))
            << "compare best_mv32x8 error";
        EXPECT_EQ(
            0,
            memcmp(best_sad8x32[0], best_sad8x32[1], sizeof(best_sad8x32[0])))
            << "compare best_sad8x32 error";
        EXPECT_EQ(
            0, memcmp(best_mv8x32[0], best_mv8x32[1], sizeof(best_mv8x32[0])))
            << "compare best_mv8x32 error";
        EXPECT_EQ(
            0,
            memcmp(
                best_sad64x16[0], best_sad64x16[1], sizeof(best_sad64x16[0])))
            << "compare best_sad64x16 error";
        EXPECT_EQ(
            0,
            memcmp(best_mv64x16[0], best_mv64x16[1], sizeof(best_mv64x16[0])))
            << "compare best_mv64x16 error";
        EXPECT_EQ(
            0,
            memcmp(
                best_sad16x64[0], best_sad16x64[1], sizeof(best_sad16x64[0])))
            << "compare best_sad16x64 error";
        EXPECT_EQ(
            0,
            memcmp(best_mv16x64[0], best_mv16x64[1], sizeof(best_mv16x64[0])))
            << "compare best_mv16x64 error";
    }
};

TEST_P(AllSadCalculationTest, 8x8_16x16_Test) {
    check_get_8x8_sad();
}

TEST_P(AllSadCalculationTest, 32x32_64x64_Test) {
    check_get_32x32_sad();
}

TEST_P(AllSadCalculationTest, nsq_sad_Test) {
    check_get_nsq_sad();
}

INSTANTIATE_TEST_CASE_P(
    ALLSAD, AllSadCalculationTest,
    ::testing::Combine(::testing::ValuesIn(TEST_PATTERNS),
                       ::testing::ValuesIn(TEST_SAD_PATTERNS)));
/**
 * @brief Unit test for ExtSadCalculation Test functions include:
 *  -
 * ext_sad_calculation_8x8_16x16_avx2_intrin
 * ext_sad_calculation_32x32_64x64_sse4_intrin
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

class ExtSadCalculationTest
    : public ::testing::WithParamInterface<SadCalTestParam>,
      public SADTestBase {
  public:
    ExtSadCalculationTest()
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
        fill_buf_with_value(&best_sad8x8[0][0], 2 * 4, BEST_SAD_MAX);
        fill_buf_with_value(&best_sad16x16[0], 2, UINT_MAX);
        fill_buf_with_value(&sad16x16[0], 2, UINT_MAX);
        fill_buf_with_value(&sad8x8[0][0], 2 * 4, UINT_MAX);

        prepare_data();

        ext_sad_calculation_8x8_16x16_c(src_aligned_,
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
                                        false);

        ext_sad_calculation_8x8_16x16_avx2_intrin(src_aligned_,
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
                                                  false);

        EXPECT_EQ(
            0, memcmp(best_sad8x8[0], best_sad8x8[1], sizeof(best_sad8x8[0])))
            << "compare best_sad8x8 error";
        EXPECT_EQ(0,
                  memcmp(best_mv8x8[0], best_mv8x8[1], sizeof(best_mv8x8[0])))
            << "compare best_mv8x8 error";
        EXPECT_EQ(best_sad16x16[0], best_sad16x16[1])
            << "compare best_sad16x16 error";
        EXPECT_EQ(best_mv16x16[0], best_mv16x16[1])
            << "compare best_mv16x16 error";
        EXPECT_EQ(0, memcmp(sad8x8[0], sad8x8[1], sizeof(sad8x8[0])))
            << "compare sad8x8 error";
        EXPECT_EQ(sad16x16[0], sad16x16[1]) << "compare sad16x16 error";
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

        ext_sad_calculation_32x32_64x64_c(*sad16x16_32b,
                                          best_sad32x32[0],
                                          &best_sad64x64[0],
                                          best_mv32x32[0],
                                          &best_mv64x64[0],
                                          0,
                                          sad32x32[0]);

        ext_sad_calculation_32x32_64x64_sse4_intrin(*sad16x16_32b,
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

TEST_P(ExtSadCalculationTest, ExtSad8x8Test) {
    check_get_8x8_sad();
}

TEST_P(ExtSadCalculationTest, ExtSad32x32Test) {
    check_get_32x32_sad();
}

INSTANTIATE_TEST_CASE_P(
    EXTSAD, ExtSadCalculationTest,
    ::testing::Combine(::testing::ValuesIn(TEST_PATTERNS),
                       ::testing::ValuesIn(TEST_SAD_PATTERNS)));

/**
 * @brief Unit test for SadCalculation Test functions include:
 *  -
 * sad_calculation_8x8_16x16_sse2_intrin
 * sad_calculation_32x32_64x64_sse2_intrin
 *
 * Test strategy:
 *  This test use different test pattern {REF_MAX, SRC_MAX, RANDOM, UNALIGN}
 *  to generate test vector,sad pattern {BUF_MAX, BUF_MIN, BUF_RANDOM} to
 *  generate test sad16x16.Check the result by compare result from reference
 *  function, non_SSE2 function and SSE2 function.
 *
 *
 * Expect result:
 *  Results come from  non_SSE2 function and SSE2 funtion are
 * equal.
 *
 * Test coverage:
 *
 * Test cases:
 **/

class SadCalculationTest
    : public ::testing::WithParamInterface<SadCalTestParam>,
      public SADTestBase {
  public:
    SadCalculationTest()
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
        fill_buf_with_value(&best_sad8x8[0][0], 2 * 4, BEST_SAD_MAX);
        fill_buf_with_value(&best_sad16x16[0], 2, UINT_MAX);
        fill_buf_with_value(&sad16x16[0], 2, UINT_MAX);
        fill_buf_with_value(&sad8x8[0][0], 2 * 4, UINT_MAX);

        prepare_data();

        sad_calculation_8x8_16x16_sse2_intrin(src_aligned_,
                                              src_stride_,
                                              ref1_aligned_,
                                              ref1_stride_,
                                              best_sad8x8[0],
                                              &best_sad16x16[0],
                                              best_mv8x8[0],
                                              &best_mv16x16[0],
                                              0,
                                              &sad16x16[0],
                                              false);

        sad_calculation_8x8_16x16_c(src_aligned_,
                                    src_stride_,
                                    ref1_aligned_,
                                    ref1_stride_,
                                    best_sad8x8[1],
                                    &best_sad16x16[1],
                                    best_mv8x8[1],
                                    &best_mv16x16[1],
                                    0,
                                    &sad16x16[1],
                                    false);

        EXPECT_EQ(
            0, memcmp(best_sad8x8[0], best_sad8x8[1], sizeof(best_sad8x8[0])))
            << "compare best_sad8x8 error";
        EXPECT_EQ(0,
                  memcmp(best_mv8x8[0], best_mv8x8[1], sizeof(best_mv8x8[0])))
            << "compare best_mv8x8 error";
        EXPECT_EQ(best_sad16x16[0], best_sad16x16[1])
            << "compare best_sad16x16 error";
        EXPECT_EQ(best_mv16x16[0], best_mv16x16[1])
            << "compare best_mv16x16 error";
        EXPECT_EQ(0, memcmp(sad8x8[0], sad8x8[1], sizeof(sad8x8[0])))
            << "compare sad8x8 error";
        EXPECT_EQ(sad16x16[0], sad16x16[1]) << "compare sad16x16 error";
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

        sad_calculation_32x32_64x64_sse2_intrin(*sad16x16_32b,
                                                best_sad32x32[0],
                                                &best_sad64x64[0],
                                                best_mv32x32[0],
                                                &best_mv64x64[0],
                                                0);

        sad_calculation_32x32_64x64_c(*sad16x16_32b,
                                      best_sad32x32[1],
                                      &best_sad64x64[1],
                                      best_mv32x32[1],
                                      &best_mv64x64[1],
                                      0);

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
    }
};

TEST_P(SadCalculationTest, Sad8x8Test) {
    check_get_8x8_sad();
}

TEST_P(SadCalculationTest, Sad32x32Test) {
    check_get_32x32_sad();
}

INSTANTIATE_TEST_CASE_P(
    CALSAD, SadCalculationTest,
    ::testing::Combine(::testing::ValuesIn(TEST_PATTERNS),
                       ::testing::ValuesIn(TEST_SAD_PATTERNS)));

/**
 * @brief Unit test for combined averaging ssd functions include:
 *  - combined_averaging_ssd_{c,avx2}
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

typedef uint32_t (*combined_averaging_ssd_func)(
    uint8_t *src, ptrdiff_t src_stride, uint8_t *ref1, ptrdiff_t ref1_stride,
    uint8_t *ref2, ptrdiff_t ref2_stride, uint32_t height, uint32_t width);

static const combined_averaging_ssd_func combined_averaging_ssd_func_table[] = {
    combined_averaging_ssd_avx2,
#ifndef NON_AVX512_SUPPORT
    combined_averaging_ssd_avx512
#endif
};

class SSDAvgTest : public ::testing::WithParamInterface<TestSadParam>,
                   public SADTestBase {
  public:
    SSDAvgTest()
        : SADTestBase(std::get<0>(TEST_GET_PARAM(1)),
                      std::get<1>(TEST_GET_PARAM(1)), TEST_GET_PARAM(0)) {
    }

  protected:
    void check_ssd_loop() {
        prepare_data();

        const uint32_t sum0_ssd = combined_averaging_ssd_c(src_aligned_,
                                                           src_stride_,
                                                           ref1_aligned_,
                                                           ref1_stride_,
                                                           ref2_aligned_,
                                                           ref2_stride_,
                                                           height_,
                                                           width_);

        for (int i = 0; i < sizeof(combined_averaging_ssd_func_table) /
                                sizeof(*combined_averaging_ssd_func_table);
             i++) {
            const uint32_t sum1_ssd =
                combined_averaging_ssd_func_table[i](src_aligned_,
                                                     src_stride_,
                                                     ref1_aligned_,
                                                     ref1_stride_,
                                                     ref2_aligned_,
                                                     ref2_stride_,
                                                     height_,
                                                     width_);

            EXPECT_EQ(sum0_ssd, sum1_ssd)
                << "compare sum combined averaging ssd error"
                << " block dim: [" << width_ << " x " << height_ << "] ";
        }
    }

    void check_ssd_speed() {
        uint32_t sum0_ssd;
        double time_c, time_o;
        uint64_t start_time_seconds, start_time_useconds;
        uint64_t middle_time_seconds, middle_time_useconds;
        uint64_t finish_time_seconds, finish_time_useconds;

        const uint64_t num_loop = 1000000000 / (width_ * height_);

        prepare_data();

        EbStartTime(&start_time_seconds, &start_time_useconds);

        for (uint64_t i = 0; i < num_loop; i++) {
            sum0_ssd = combined_averaging_ssd_c(src_aligned_,
                                                src_stride_,
                                                ref1_aligned_,
                                                ref1_stride_,
                                                ref2_aligned_,
                                                ref2_stride_,
                                                height_,
                                                width_);
        }

        EbStartTime(&middle_time_seconds, &middle_time_useconds);
        EbComputeOverallElapsedTimeMs(start_time_seconds,
                                      start_time_useconds,
                                      middle_time_seconds,
                                      middle_time_useconds,
                                      &time_c);

        for (int i = 0; i < sizeof(combined_averaging_ssd_func_table) /
                                sizeof(*combined_averaging_ssd_func_table);
             i++) {
            uint32_t sum1_ssd;

            EbStartTime(&middle_time_seconds, &middle_time_useconds);

            for (uint64_t j = 0; j < num_loop; j++) {
                sum1_ssd = combined_averaging_ssd_func_table[i](src_aligned_,
                                                                src_stride_,
                                                                ref1_aligned_,
                                                                ref1_stride_,
                                                                ref2_aligned_,
                                                                ref2_stride_,
                                                                height_,
                                                                width_);
            }

            EbStartTime(&finish_time_seconds, &finish_time_useconds);
            EbComputeOverallElapsedTimeMs(middle_time_seconds,
                                          middle_time_useconds,
                                          finish_time_seconds,
                                          finish_time_useconds,
                                          &time_o);

            EXPECT_EQ(sum0_ssd, sum1_ssd)
                << "compare sum combined averaging ssd error"
                << " block dim: [" << width_ << " x " << height_ << "] ";

            printf("combined_averaging_ssd(%3dx%3d): %6.2f\n",
                   width_,
                   height_,
                   time_c / time_o);
        }
    }
};

TEST_P(SSDAvgTest, SSDTest) {
    check_ssd_loop();
}

TEST_P(SSDAvgTest, DISABLED_SSDSpeedTest) {
    check_ssd_speed();
}

INSTANTIATE_TEST_CASE_P(
    SSDAvg, SSDAvgTest,
    ::testing::Combine(::testing::ValuesIn(TEST_PATTERNS),
                       ::testing::ValuesIn(TEST_BLOCK_SIZES)));
}  // namespace
