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
#include "EbComputeSAD.h"
#include "EbComputeSAD_SSE4_1.h"
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
#define MAX_BLOCK_SIZE (128 * 256)
typedef std::tuple<int, int> SadSize;
typedef enum { REF_MAX, SRC_MAX, RANDOM, UNALIGN } TestPattern;
typedef enum { BUF_MAX, BUF_MIN, BUF_RANDOM } SADPattern;
SadSize TEST_SAD_SIZES[] = {SadSize(64, 64),
                            SadSize(64, 32),
                            SadSize(32, 64),
                            SadSize(32, 32),
                            SadSize(32, 16),
                            SadSize(16, 32),
                            SadSize(16, 16),
                            SadSize(16, 8),
                            SadSize(8, 16),
                            SadSize(8, 8),
                            SadSize(8, 4),
                            SadSize(4, 8),
                            SadSize(4, 16),
                            SadSize(16, 4),
                            SadSize(8, 32),
                            SadSize(32, 8),
                            SadSize(16, 64),
                            SadSize(64, 16)};
TestPattern TEST_PATTERNS[] = {REF_MAX, SRC_MAX, RANDOM, UNALIGN};
SADPattern TEST_SAD_PATTERNS[] = {BUF_MAX, BUF_MIN, BUF_RANDOM};

typedef std::tuple<TestPattern, SadSize> TestSadParam;

/**
 * @Brief Base class for SAD test. SADTestBase handle test vector in memory,
 * provide SAD and SAD avg reference function
 */
class SADTestBase : public ::testing::Test {
  public:
    SADTestBase(const int width, const int height, TestPattern test_pattern) {
        width_ = width;
        height_ = height;
        src_stride_ = ref1_stride_ = ref2_stride_ = width_ * 2;
        test_pattern_ = test_pattern;
        src_aligned_ = nullptr;
        ref1_aligned_ = nullptr;
        ref2_aligned_ = nullptr;
    }

    SADTestBase(const int width, const int height, TestPattern test_pattern,
                SADPattern test_sad_pattern) {
        width_ = width;
        height_ = height;
        src_stride_ = ref1_stride_ = ref2_stride_ = width_ * 2;
        test_pattern_ = test_pattern;
        test_sad_pattern_ = test_sad_pattern;
    }

    SADTestBase(const int width, const int height, TestPattern test_pattern,
                const int search_area_width, const int search_area_height) {
        width_ = width;
        height_ = height;
        src_stride_ = ref1_stride_ = ref2_stride_ = width_ * 2;
        test_pattern_ = test_pattern;
        search_area_width_ = search_area_width;
        search_area_height_ = search_area_height;
    }

    void SetUp() override {
        src_aligned_ = (uint8_t * )eb_aom_memalign(32, MAX_BLOCK_SIZE);
        ref1_aligned_ = (uint8_t * )eb_aom_memalign(32, MAX_BLOCK_SIZE);
        ref2_aligned_ = (uint8_t * )eb_aom_memalign(32, MAX_BLOCK_SIZE);
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
            for (int i = 0; i < MAX_BLOCK_SIZE; i++) {
                src_aligned_[i] = 0;
                ref1_aligned_[i] = ref2_aligned_[i] = mask;
            }
            break;
        }
        case SRC_MAX: {
            for (int i = 0; i < MAX_BLOCK_SIZE; i++) {
                src_aligned_[i] = mask;
                ref1_aligned_[i] = ref2_aligned_[i] = 0;
            }
            break;
        }
        case RANDOM: {
            for (int i = 0; i < MAX_BLOCK_SIZE; i++) {
                src_aligned_[i] = rnd.random();
                ref1_aligned_[i] = rnd.random();
                ref2_aligned_[i] = rnd.random();
            }
            break;
        };
        case UNALIGN: {
            for (int i = 0; i < MAX_BLOCK_SIZE; i++) {
                src_aligned_[i] = rnd.random();
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

    void fill_buf_with_value(uint32_t *buf, int num, uint32_t value) {
        for (int i = 0; i < num; ++i)
            buf[i] = value;
    }

    void prepare_sad_data() {
        const int32_t mask = (1 << 8) - 1;
        SVTRandom rnd(0, mask);
        switch (test_sad_pattern_) {
        case BUF_MAX: {
            fill_buf_with_value(&sad16x16[0][0], 16 * 8, mask);
            break;
        }
        case BUF_MIN: {
            fill_buf_with_value(&sad16x16[0][0], 16 * 8, 0);
            break;
        }
        case BUF_RANDOM: {
            for (int i = 0; i < 16; i++)
                for (int j = 0; j < 8; j++)
                    sad16x16[i][j] = rnd.random();
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
            fill_buf_with_value(&sad16x16[0][0], 16 * 8, mask);
            fill_buf_with_value(&sad32x32[0][0], 4 * 8, mask);
            break;
        }
        case BUF_MIN: {
            fill_buf_with_value(&sad8x8[0][0], 64 * 8, 0);
            fill_buf_with_value(&sad16x16[0][0], 16 * 8, 0);
            fill_buf_with_value(&sad32x32[0][0], 4 * 8, 0);
            break;
        }
        case BUF_RANDOM: {
            for (int i = 0; i < 64; i++)
                for (int j = 0; j < 8; j++)
                    sad8x8[i][j] = rnd.random();

            for (int i = 0; i < 16; i++)
                for (int j = 0; j < 8; j++)
                    sad16x16[i][j] = rnd.random();

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
    uint32_t sad8x8[64][8];
    uint32_t sad16x16[16][8];
    uint32_t sad32x32[4][8];
};

/**
 * @brief Unit test for SAD sub smaple functions include:
 *  - nxm_sad_kernel_sub_sampled_func_ptr_array
 *  - nxm_sad_kernel_sub_sampled_func_ptr_array
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
 *  All functions inside nxm_sad_kernel_sub_sampled_func_ptr_array and
 * nxm_sad_kernel_sub_sampled_func_ptr_array.
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
        EbSadKernelNxMType non_avx2_func =
            nxm_sad_kernel_sub_sampled_func_ptr_array[ASM_NON_AVX2]
                                                     [width_ >> 3];
        EbSadKernelNxMType avx2_func =
            nxm_sad_kernel_sub_sampled_func_ptr_array[ASM_AVX2][width_ >> 3];

        uint32_t ref_sad = 0;
        uint32_t non_avx2_sad = 0;
        uint32_t avx2_sad = 0;

        prepare_data();

        ref_sad = reference_sad();
        if (non_avx2_func != nullptr)
            non_avx2_sad = non_avx2_func(src_aligned_,
                                         src_stride_,
                                         ref1_aligned_,
                                         ref1_stride_,
                                         height_,
                                         width_);
        if (avx2_func != nullptr)
            avx2_sad = avx2_func(src_aligned_,
                                 src_stride_,
                                 ref1_aligned_,
                                 ref1_stride_,
                                 height_,
                                 width_);

        if (non_avx2_func != nullptr)
            EXPECT_EQ(ref_sad, non_avx2_sad)
                << "compare ref and non_avx2 error";

        if (avx2_func != nullptr)
            EXPECT_EQ(ref_sad, avx2_sad) << "compare ref and avx2 error";

        if (non_avx2_func != nullptr && avx2_func != nullptr)
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
                       ::testing::ValuesIn(TEST_SAD_SIZES)));
/**
 * @brief Unit test for SAD functions include:
 *  - nxm_sad_kernel_func_ptr_array
 *  - nxm_sad_kernel_func_ptr_array
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
 *  All functions inside nxm_sad_kernel_func_ptr_array and
 *  nxm_sad_kernel_func_ptr_array.
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
        EbSadKernelNxMType non_avx2_func =
            nxm_sad_kernel_func_ptr_array[ASM_NON_AVX2][width_ >> 3];
        EbSadKernelNxMType avx2_func =
            nxm_sad_kernel_func_ptr_array[ASM_AVX2][width_ >> 3];

        uint32_t ref_sad = 0;
        uint32_t non_avx2_sad = 0;
        uint32_t avx2_sad = 0;

        prepare_data();

        ref_sad = reference_sad();
        if (non_avx2_func != nullptr)
            non_avx2_sad = non_avx2_func(src_aligned_,
                                         src_stride_,
                                         ref1_aligned_,
                                         ref1_stride_,
                                         height_,
                                         width_);
        if (avx2_func != nullptr)
            avx2_sad = avx2_func(src_aligned_,
                                 src_stride_,
                                 ref1_aligned_,
                                 ref1_stride_,
                                 height_,
                                 width_);

        if (non_avx2_func != nullptr)
            EXPECT_EQ(ref_sad, non_avx2_sad)
                << "compare ref and non_avx2 error";
        if (avx2_func != nullptr)
            EXPECT_EQ(ref_sad, avx2_sad) << "compare ref and avx2 error";
        if (non_avx2_func != nullptr && avx2_func != nullptr)
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
                       ::testing::ValuesIn(TEST_SAD_SIZES)));
/**
 * @brief Unit test for SAD Avg functions include:
 *  - nxm_sad_averaging_kernel_func_ptr_array
 *  - nxm_sad_averaging_kernel_func_ptr_array
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
 *  All functions inside nxm_sad_averaging_kernel_func_ptr_array and
 * nxm_sad_averaging_kernel_func_ptr_array.
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
        EbSadAvgKernelNxMType non_avx2_func =
            nxm_sad_averaging_kernel_func_ptr_array[ASM_NON_AVX2][width_ >> 3];
        EbSadAvgKernelNxMType avx2_func =
            nxm_sad_averaging_kernel_func_ptr_array[ASM_AVX2][width_ >> 3];

        uint32_t ref_sad = 0;
        uint32_t non_avx2_sad = 0;
        uint32_t avx2_sad = 0;

        prepare_data();

        ref_sad = reference_sad_avg();
        if (non_avx2_func != nullptr)
            non_avx2_sad = non_avx2_func(src_aligned_,
                                         src_stride_,
                                         ref1_aligned_,
                                         ref1_stride_,
                                         ref2_aligned_,
                                         ref2_stride_,
                                         height_,
                                         width_);
        if (avx2_func != nullptr)
            avx2_sad = avx2_func(src_aligned_,
                                 src_stride_,
                                 ref1_aligned_,
                                 ref1_stride_,
                                 ref2_aligned_,
                                 ref2_stride_,
                                 height_,
                                 width_);

        if (non_avx2_func != nullptr)
            EXPECT_EQ(ref_sad, non_avx2_sad)
                << "compare ref and non_avx2 error";
        if (avx2_func != nullptr)
            EXPECT_EQ(ref_sad, avx2_sad) << "compare ref and avx2 error";
        if (non_avx2_func != nullptr && avx2_func != nullptr)
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
                       ::testing::ValuesIn(TEST_SAD_SIZES)));

typedef std::tuple<int16_t, int16_t> SearchArea;
SadSize TEST_LOOP_SIZES[] = {
    SadSize(64, 64), SadSize(64, 32), SadSize(32, 64), SadSize(32, 32),
    SadSize(32, 16), SadSize(16, 32), SadSize(16, 16), SadSize(16, 8),
    SadSize(8, 16),  SadSize(8, 8),   SadSize(8, 4),   SadSize(4, 4),
    SadSize(4, 8),   SadSize(4, 16),  SadSize(16, 4),  SadSize(8, 32),
    SadSize(32, 8),  SadSize(16, 64), SadSize(64, 16), SadSize(24, 24),
    SadSize(24, 16), SadSize(16, 24), SadSize(24, 8),  SadSize(8, 24),
    SadSize(64, 24), SadSize(48, 24), SadSize(32, 24), SadSize(24, 32),
    SadSize(48, 48), SadSize(48, 16), SadSize(48, 32), SadSize(16, 48),
    SadSize(32, 48), SadSize(48, 64), SadSize(64, 48)};
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
typedef std::tuple<TestPattern, SadSize, SearchArea> SadLoopTestParam;

/**
 * @brief Unit test for SAD loop (sparse) functions include:
 *  - sad_loop_kernel_sparse_{sse4_1,avx2}_intrin
 *  - sad_loop_kernel_{sse4_1,avx2}_intrin
 *  - sad_loop_kernel_sse4_1_hme_l0_intrin
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
                      std::get<1>(TEST_GET_PARAM(2))) {
        src_stride_ = width_ * 2;
        ref1_stride_ = ref2_stride_ = 128;
    };

  protected:
    void check_sad_loop_sparse() {
        uint64_t best_sad1 = UINT64_MAX;
        int16_t x_search_center1 = 0;
        int16_t y_search_center1 = 0;

        prepare_data();

        sad_loop_kernel_sparse_sse4_1_intrin(src_aligned_,
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
        uint64_t best_sad2 = UINT64_MAX;
        int16_t x_search_center2 = 0;
        int16_t y_search_center2 = 0;
        sad_loop_kernel_sparse_avx2_intrin(src_aligned_,
                                           src_stride_,
                                           ref1_aligned_,
                                           ref1_stride_,
                                           height_,
                                           width_,
                                           &best_sad2,
                                           &x_search_center2,
                                           &y_search_center2,
                                           ref1_stride_,
                                           search_area_width_,
                                           search_area_height_);
        uint64_t best_sad3 = UINT64_MAX;
        int16_t x_search_center3 = 0;
        int16_t y_search_center3 = 0;
        sad_loop_kernel_sparse(src_aligned_,
                               src_stride_,
                               ref1_aligned_,
                               ref1_stride_,
                               height_,
                               width_,
                               &best_sad3,
                               &x_search_center3,
                               &y_search_center3,
                               ref1_stride_,
                               search_area_width_,
                               search_area_height_);

        EXPECT_EQ(best_sad1, best_sad3)
            << "compare best_sad error"
            << " block dim: [" << width_ << " x " << height_ << "] "
            << "search area [" << search_area_width_ << " x "
            << search_area_height_ << "]";
        EXPECT_EQ(x_search_center1, x_search_center3)
            << "compare x_search_center error"
            << " block dim: [" << width_ << " x " << height_ << "] "
            << "search area [" << search_area_width_ << " x "
            << search_area_height_ << "]";
        EXPECT_EQ(y_search_center1, y_search_center3)
            << "compare y_search_center error"
            << " block dim: [" << width_ << " x " << height_ << "] "
            << "search area [" << search_area_width_ << " x "
            << search_area_height_ << "]";
        EXPECT_EQ(best_sad2, best_sad3)
            << "compare best_sad error"
            << " block dim: [" << width_ << " x " << height_ << "] "
            << "search area [" << search_area_width_ << " x "
            << search_area_height_ << "]";
        EXPECT_EQ(x_search_center2, x_search_center3)
            << "compare x_search_center error"
            << " block dim: [" << width_ << " x " << height_ << "] "
            << "search area [" << search_area_width_ << " x "
            << search_area_height_ << "]";
        EXPECT_EQ(y_search_center2, y_search_center3)
            << "compare y_search_center error"
            << " block dim: [" << width_ << " x " << height_ << "] "
            << "search area [" << search_area_width_ << " x "
            << search_area_height_ << "]";
    }

    void check_sad_loop() {
        uint64_t best_sad1 = UINT64_MAX;
        int16_t x_search_center1 = 0;
        int16_t y_search_center1 = 0;

        prepare_data();

        sad_loop_kernel_sse4_1_intrin(src_aligned_,
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
        uint64_t best_sad2 = UINT64_MAX;
        int16_t x_search_center2 = 0;
        int16_t y_search_center2 = 0;
        sad_loop_kernel_avx2_intrin(src_aligned_,
                                    src_stride_,
                                    ref1_aligned_,
                                    ref1_stride_,
                                    height_,
                                    width_,
                                    &best_sad2,
                                    &x_search_center2,
                                    &y_search_center2,
                                    ref1_stride_,
                                    search_area_width_,
                                    search_area_height_);

        uint64_t best_sad3 = UINT64_MAX;
        int16_t x_search_center3 = 0;
        int16_t y_search_center3 = 0;
        sad_loop_kernel(src_aligned_,
                        src_stride_,
                        ref1_aligned_,
                        ref1_stride_,
                        height_,
                        width_,
                        &best_sad3,
                        &x_search_center3,
                        &y_search_center3,
                        ref1_stride_,
                        search_area_width_,
                        search_area_height_);

        EXPECT_EQ(best_sad1, best_sad3)
            << "compare bast_sad error"
            << " block dim: [" << width_ << " x " << height_ << "] "
            << "search area [" << search_area_width_ << " x "
            << search_area_height_ << "]";
        EXPECT_EQ(best_sad2, best_sad3)
            << "compare bast_sad error"
            << " block dim: [" << width_ << " x " << height_ << "] "
            << "search area [" << search_area_width_ << " x "
            << search_area_height_ << "]";
        EXPECT_EQ(x_search_center1, x_search_center3)
            << "compare x_search_center error"
            << " block dim: [" << width_ << " x " << height_ << "] "
            << "search area [" << search_area_width_ << " x "
            << search_area_height_ << "]";
        EXPECT_EQ(x_search_center2, x_search_center3)
            << "compare x_search_center error"
            << " block dim: [" << width_ << " x " << height_ << "] "
            << "search area [" << search_area_width_ << " x "
            << search_area_height_ << "]";
        EXPECT_EQ(y_search_center1, y_search_center3)
            << "compare y_search_center error"
            << " block dim: [" << width_ << " x " << height_ << "] "
            << "search area [" << search_area_width_ << " x "
            << search_area_height_ << "]";
        EXPECT_EQ(y_search_center2, y_search_center3)
            << "compare y_search_center error"
            << " block dim: [" << width_ << " x " << height_ << "] "
            << "search area [" << search_area_width_ << " x "
            << search_area_height_ << "]";
    }

    void check_hme_loop() {
        uint64_t best_sad1 = UINT64_MAX;
        int16_t x_search_center1 = 0;
        int16_t y_search_center1 = 0;

        prepare_data();

        sad_loop_kernel(src_aligned_,
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
        uint64_t best_sad2 = UINT64_MAX;
        int16_t x_search_center2 = 0;
        int16_t y_search_center2 = 0;
        sad_loop_kernel_sse4_1_hme_l0_intrin(src_aligned_,
                                             src_stride_,
                                             ref1_aligned_,
                                             ref1_stride_,
                                             height_,
                                             width_,
                                             &best_sad2,
                                             &x_search_center2,
                                             &y_search_center2,
                                             ref1_stride_,
                                             search_area_width_,
                                             search_area_height_);
        EXPECT_EQ(best_sad1, best_sad2)
            << "compare hme best_sad error"
            << " block dim: [" << width_ << " x " << height_ << "] "
            << "search area [" << search_area_width_ << " x "
            << search_area_height_ << "]";
        EXPECT_EQ(x_search_center1, x_search_center2)
            << "compare hme x_search_center error"
            << " block dim: [" << width_ << " x " << height_ << "] "
            << "search area [" << search_area_width_ << " x "
            << search_area_height_ << "]";
        EXPECT_EQ(y_search_center1, y_search_center2)
            << "compare hme y_search_center error"
            << " block dim: [" << width_ << " x " << height_ << "] "
            << "search area [" << search_area_width_ << " x "
            << search_area_height_ << "]";
    }
};

TEST_P(SadLoopTest, SadLoopSparseTest) {
    check_sad_loop_sparse();
}

TEST_P(SadLoopTest, SadLoopTest) {
    check_sad_loop();
}

TEST_P(SadLoopTest, HmeTest) {
    check_hme_loop();
}

INSTANTIATE_TEST_CASE_P(LOOPSAD, SadLoopTest,
                        ::testing::Combine(::testing::ValuesIn(TEST_PATTERNS),
                                           ::testing::ValuesIn(TEST_LOOP_SIZES),
                                           ::testing::ValuesIn(TEST_AREAS)));

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
        src_stride_ = ref1_stride_ = 256;
    }

  protected:
    void check_get_eight() {
        uint32_t best_sad8x8_1[4] = {UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX};
        uint32_t best_mv8x8_1[4] = {0};
        uint32_t best_sad16x16_1 = UINT_MAX, best_mv16x16_1 = 0;
        uint16_t sad16x16_1[8] = {0};

        prepare_data();

        get_eight_horizontal_search_point_results_8x8_16x16_pu_sse41_intrin(
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

        uint32_t best_sad8x8_2[4] = {UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX};
        uint32_t best_mv8x8_2[4] = {0};
        uint32_t best_sad16x16_2 = UINT_MAX, best_mv16x16_2 = 0;
        uint16_t sad16x16_2[8] = {0};
        get_eight_horizontal_search_point_results_8x8_16x16_pu_avx2_intrin(
            src_aligned_,
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

        EXPECT_EQ(0,
                  memcmp(best_sad8x8_1, best_sad8x8_2, sizeof(best_sad8x8_1)))
            << "compare best_sad8x8 error";
        EXPECT_EQ(0, memcmp(best_mv8x8_1, best_mv8x8_2, sizeof(best_mv8x8_1)))
            << "compare best_mv8x8 error";
        EXPECT_EQ(best_sad16x16_1, best_sad16x16_2)
            << "compare best_sad16x16 error";
        EXPECT_EQ(best_mv16x16_1, best_mv16x16_2)
            << "compare best_mv16x16 error";
        EXPECT_EQ(0, memcmp(sad16x16_1, sad16x16_2, sizeof(sad16x16_1)))
            << "compare sad16x16 error";
    }

    void check_get_eight_32_64() {
        uint32_t best_sad32x32_1[4] = {
            BEST_SAD_MAX, BEST_SAD_MAX, BEST_SAD_MAX, BEST_SAD_MAX};
        uint32_t best_mv32x32_1[4] = {0};
        uint32_t best_sad64x64_1 = UINT_MAX, best_mv64x64_1 = 0;

        prepare_sad_data();

        get_eight_horizontal_search_point_results_32x32_64x64_pu_sse41_intrin(
            (uint16_t *)sad16x16,
            best_sad32x32_1,
            &best_sad64x64_1,
            best_mv32x32_1,
            &best_mv64x64_1,
            0);

        uint32_t best_sad32x32_2[4] = {
            BEST_SAD_MAX, BEST_SAD_MAX, BEST_SAD_MAX, BEST_SAD_MAX};
        uint32_t best_mv32x32_2[4] = {0};
        uint32_t best_sad64x64_2 = UINT_MAX, best_mv64x64_2 = 0;
        get_eight_horizontal_search_point_results_32x32_64x64_pu_avx2_intrin(
            (uint16_t *)sad16x16,
            best_sad32x32_2,
            &best_sad64x64_2,
            best_mv32x32_2,
            &best_mv64x64_2,
            0);

        EXPECT_EQ(
            0,
            memcmp(best_sad32x32_1, best_sad32x32_2, sizeof(best_sad32x32_1)))
            << "compare best_sad32x32 error";
        EXPECT_EQ(
            0, memcmp(best_mv32x32_1, best_mv32x32_2, sizeof(best_mv32x32_1)))
            << "compare best_mv32x32 error";
        EXPECT_EQ(best_sad64x64_1, best_sad64x64_2)
            << "compare best_sad64x64 error";
        EXPECT_EQ(best_mv64x64_1, best_mv64x64_2)
            << "compare best_mv64x64 error";
    }
};

TEST_P(GetEightSadTest, GetEightSadTest) {
    check_get_eight();
}

TEST_P(GetEightSadTest, SadTest_32_64) {
    check_get_eight_32_64();
}

INSTANTIATE_TEST_CASE_P(
    EIGHTSAD, GetEightSadTest,
    ::testing::Combine(::testing::ValuesIn(TEST_PATTERNS),
                       ::testing::ValuesIn(TEST_SAD_PATTERNS)));

/**
 * @brief Unit test for AllSadCalculation Test functions include:
 *  -
 * ext_all_sad_calculation_8x8_16x16_{c,avx2}
 * ext_eight_sad_calculation_32x32_64x64_{c,avx2}
 * ext_eigth_sad_calculation_nsq_{c£¬avx2}
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
        src_stride_ = ref1_stride_ = 256;
    };

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

        prepare_sad_data();

        ext_eight_sad_calculation_32x32_64x64_c(sad16x16,
                                                best_sad32x32[0],
                                                &best_sad64x64[0],
                                                best_mv32x32[0],
                                                &best_mv64x64[0],
                                                0,
                                                sad32x32[0]);

        ext_eight_sad_calculation_32x32_64x64_avx2(sad16x16,
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
                                        sad16x16,
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
                                           sad16x16,
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
 * ext_sad_calculation_8x8_16x16
 * ext_sad_calculation_8x8_16x16_avx2_intrin
 * ext_sad_calculation_32x32_64x64
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
        src_stride_ = ref1_stride_ = 256;
    };

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

        ext_sad_calculation_8x8_16x16(src_aligned_,
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

        prepare_sad_data();

        ext_sad_calculation_32x32_64x64(*sad16x16,
                                        best_sad32x32[0],
                                        &best_sad64x64[0],
                                        best_mv32x32[0],
                                        &best_mv64x64[0],
                                        0,
                                        sad32x32[0]);

        ext_sad_calculation_32x32_64x64_sse4_intrin(*sad16x16,
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
}  // namespace
