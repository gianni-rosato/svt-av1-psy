/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file VarianceTest.cc
 *
 * @brief Unit test for variance, mse, sum square functions:
 * - aom_variance{4-128}x{4-128}_{c,sse2,avx2}
 * - aom_get_mb_ss_sse2
 * - aom_mse16x16_{c,avx2}
 *
 * @author Cidana-Ryan
 *
 ******************************************************************************/
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <new>
// workaround to eliminate the compiling warning on linux
// The macro will conflict with definition in gtest.h
#ifdef __USE_GNU
#undef __USE_GNU  // defined in EbThreads.h
#endif
#ifdef _GNU_SOURCE
#undef _GNU_SOURCE  // defined in EbThreads.h
#endif
#include "aom_dsp_rtcd.h"
#include "random.h"
#include "util.h"
#include "gtest/gtest.h"
#include "EbUtility.h"

using svt_av1_test_tool::SVTRandom;  // to generate the random
namespace {
#define MAX_BLOCK_SIZE (128 * 128)

// MSE test
static uint32_t mse_ref(const uint8_t *src, const uint8_t *ref, int width,
                        int height, int src_stride, int ref_stride,
                        uint32_t *sse_ptr) {
    int64_t se = 0;
    uint64_t sse = 0;
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int diff;
            diff = src[y * src_stride + x] - ref[y * ref_stride + x];
            se += diff;
            sse += diff * diff;
        }
    }
    *sse_ptr = static_cast<uint32_t>(sse);
    return static_cast<uint32_t>(sse - ((se * se) / (width * height)));
}

using MSE_NXM_FUNC = uint32_t (*)(const uint8_t *src_ptr, int32_t src_stride,
                                  const uint8_t *ref_ptr, int32_t ref_stride,
                                  uint32_t *sse);

typedef std::tuple<uint32_t, uint32_t, MSE_NXM_FUNC> TestMseParam;
/**
 * @brief Unit test for mse functions, target functions include:
 *  - aom_mse16x16_{c,avx2}
 *
 * Test strategy:
 *  This test case use random source, max source, zero source as test
 * pattern.
 *
 *
 * Expect result:
 *  Results come from reference functon and target function are
 * equal.
 *
 * Test coverage:
 *
 * Test cases:
 *
 */
class MseTest : public ::testing::TestWithParam<TestMseParam> {
  public:
    MseTest()
        : width_(TEST_GET_PARAM(0)),
          height_(TEST_GET_PARAM(1)),
          mse_tst_(TEST_GET_PARAM(2)) {
        src_data_ =
            reinterpret_cast<uint8_t *>(aom_memalign(32, MAX_BLOCK_SIZE));
        ref_data_ =
            reinterpret_cast<uint8_t *>(aom_memalign(32, MAX_BLOCK_SIZE));
    }

    ~MseTest() {
        aom_free(src_data_);
        src_data_ = nullptr;
        aom_free(ref_data_);
        ref_data_ = nullptr;
    }

    void run_match_test() {
        const int32_t mask = (1 << 8) - 1;
        SVTRandom rnd(0, mask);
        for (int i = 0; i < 10; ++i) {
            for (int j = 0; j < MAX_BLOCK_SIZE; ++j) {
                src_data_[j] = rnd.random();
                ref_data_[j] = rnd.random();
            }

            uint32_t sse_tst, sse_ref;
            mse_tst_(src_data_, width_, ref_data_, width_, &sse_tst);
            mse_ref(src_data_,
                    ref_data_,
                    width_,
                    height_,
                    width_,
                    width_,
                    &sse_ref);
            ASSERT_EQ(sse_tst, sse_ref) << "SSE Error at index: " << i;
        }
    }

    void run_max_test() {
        memset(src_data_, 255, MAX_BLOCK_SIZE);
        memset(ref_data_, 0, MAX_BLOCK_SIZE);
        uint32_t sse_tst;
        mse_tst_(src_data_, width_, ref_data_, width_, &sse_tst);
        const uint32_t expected = width_ * height_ * 255 * 255;
        ASSERT_EQ(sse_tst, expected) << "Error at MSE maximum test ";
    }

  private:
    uint8_t *src_data_;
    uint8_t *ref_data_;
    uint32_t width_;
    uint32_t height_;
    MSE_NXM_FUNC mse_tst_;
};

TEST_P(MseTest, MatchTest) {
    run_match_test();
};

TEST_P(MseTest, MaxTest) {
    run_max_test();
};

INSTANTIATE_TEST_CASE_P(Variance, MseTest,
                        ::testing::Values(TestMseParam(16, 16,
                                                       &aom_mse16x16_avx2)));

// sum of squares test
static uint32_t mb_ss_ref(const int16_t *src) {
    uint32_t res = 0;
    for (int i = 0; i < 256; ++i)
        res += src[i] * src[i];
    return res;
}

using SUM_SQUARE_FUNC = uint32_t (*)(const int16_t *src);
/**
 * @brief Unit test for sum square functions, target functions include:
 *  - aom_get_mb_ss_sse2
 *
 * Test strategy:
 *  This test case use random value, const value (0-255) as test
 * pattern.
 *
 *
 * Expect result:
 *  Results come from reference functon and target function are
 * equal.
 *
 * Test coverage:
 *
 * Test cases:
 *
 */
class SumSquareTest : public ::testing::TestWithParam<SUM_SQUARE_FUNC> {
  protected:
    void run_const_test() {
        int16_t mem[256];
        uint32_t res;
        for (int16_t v = 0; v < 256; ++v) {
            for (int i = 0; i < 256; ++i)
                mem[i] = v;

            SUM_SQUARE_FUNC sum_square_func = GetParam();
            res = sum_square_func(mem);
            ASSERT_EQ(res, 256u * (v * v))
                << "Error at sum test in const value: " << v;
        }
    }

    void run_match_test() {
        int32_t mask = (1 << 8) - 1;
        SVTRandom rnd(0, mask);
        int16_t mem[256];
        for (int i = 0; i < 100; ++i) {
            for (int j = 0; j < 256; ++j)
                mem[j] = rnd.random() - rnd.random();

            const uint32_t expected = mb_ss_ref(mem);
            SUM_SQUARE_FUNC sum_square_func = GetParam();
            uint32_t res = sum_square_func(mem);
            ASSERT_EQ(res, expected) << "Error at sum test in index: " << i;
        }
    }
};

TEST_P(SumSquareTest, ConstTest) {
    run_const_test();
};

TEST_P(SumSquareTest, MatchTest) {
    run_match_test();
};

extern "C" uint32_t aom_get_mb_ss_sse2(const int16_t *src);
INSTANTIATE_TEST_CASE_P(Variance, SumSquareTest,
                        ::testing::Values(aom_get_mb_ss_sse2));

// Variance test
using VARIANCE_NXM_FUNC = uint32_t (*)(const uint8_t *src_ptr,
                                       int32_t src_stride,
                                       const uint8_t *ref_ptr,
                                       int32_t recon_stride, uint32_t *sse);

using VarianceParam =
    std::tuple<uint32_t, uint32_t, VARIANCE_NXM_FUNC, VARIANCE_NXM_FUNC>;

/**
 * @brief Unit test for variance functions, target functions include:
 *  - - aom_variance{4-128}x{4-128}_{c,avx2}
 *
 * Test strategy:
 *  This test case contains zero test, random value test, one quarter test as
 * test pattern.
 *
 *
 * Expect result:
 *  Results come from reference functon and target function are
 * equal.
 *
 * Test coverage:
 *
 * Test cases:
 *
 */
class VarianceTest : public ::testing::TestWithParam<VarianceParam> {
  public:
    VarianceTest()
        : width_(TEST_GET_PARAM(0)),
          height_(TEST_GET_PARAM(1)),
          func_c_(TEST_GET_PARAM(2)),
          func_asm_(TEST_GET_PARAM(3)) {
        src_data_ =
            reinterpret_cast<uint8_t *>(aom_memalign(32, MAX_BLOCK_SIZE));
        ref_data_ =
            reinterpret_cast<uint8_t *>(aom_memalign(32, MAX_BLOCK_SIZE));
    }

    ~VarianceTest() {
        aom_free(src_data_);
        src_data_ = nullptr;
        aom_free(ref_data_);
        ref_data_ = nullptr;
    }

    void run_zero_test() {
        for (int i = 0; i < 256; ++i) {
            memset(src_data_, i, MAX_BLOCK_SIZE);
            for (int j = 0; j < 256; ++j) {
                memset(ref_data_, j, MAX_BLOCK_SIZE);

                uint32_t sse, var;
                var = func_asm_(src_data_, width_, ref_data_, width_, &sse);
                ASSERT_EQ(var, 0u)
                    << "Variance is mismatched, src values: " << i
                    << " ref values: " << j;
            }
        }
    }

    void run_match_test() {
        int32_t mask = (1 << 8) - 1;
        SVTRandom rnd(0, mask);
        for (int i = 0; i < 10; ++i) {
            for (int j = 0; j < MAX_BLOCK_SIZE; j++) {
                src_data_[j] = rnd.random();
                ref_data_[j] = rnd.random();
            }
            uint32_t sse_c, sse_asm, var_c, var_asm;

            var_c = func_c_(src_data_, width_, ref_data_, width_, &sse_c);
            var_asm = func_asm_(src_data_, width_, ref_data_, width_, &sse_asm);
            ASSERT_EQ(sse_c, sse_asm) << "Error at test index: " << i;
            ASSERT_EQ(var_c, var_asm) << "Error at test index: " << i;
        }
    }

    void run_one_quarter_test() {
        const int half = width_ * height_ / 2;
        memset(src_data_, 255, MAX_BLOCK_SIZE);
        memset(ref_data_, 255, half);
        memset(ref_data_ + half, 0, half);
        uint32_t sse_asm, var_asm, expected;
        var_asm = func_asm_(src_data_, width_, ref_data_, width_, &sse_asm);
        expected = width_ * height_ * 255 * 255 / 4;
        ASSERT_EQ(var_asm, expected);
    }

  private:
    uint8_t *src_data_;
    uint8_t *ref_data_;
    uint32_t width_;
    uint32_t height_;
    VARIANCE_NXM_FUNC func_c_;
    VARIANCE_NXM_FUNC func_asm_;
};

TEST_P(VarianceTest, ZeroTest) {
    run_zero_test();
};

TEST_P(VarianceTest, MatchTest) {
    run_match_test();
};

TEST_P(VarianceTest, OneQuarterTest) {
    run_one_quarter_test();
};

INSTANTIATE_TEST_CASE_P(
    Variance, VarianceTest,
    ::testing::Values(
        VarianceParam(4, 4, &aom_variance4x4_c, &aom_variance4x4_sse2),
        VarianceParam(4, 8, &aom_variance4x8_c, &aom_variance4x8_sse2),
        VarianceParam(4, 16, &aom_variance4x16_c, &aom_variance4x16_sse2),
        VarianceParam(8, 4, &aom_variance8x4_c, &aom_variance8x4_sse2),
        VarianceParam(8, 8, &aom_variance8x8_c, &aom_variance8x8_sse2),
        VarianceParam(8, 16, &aom_variance8x16_c, &aom_variance8x16_sse2),
        VarianceParam(8, 32, &aom_variance8x32_c, &aom_variance8x32_sse2),
        VarianceParam(16, 4, &aom_variance16x4_c, &aom_variance16x4_avx2),
        VarianceParam(16, 8, &aom_variance16x8_c, &aom_variance16x8_avx2),
        VarianceParam(16, 16, &aom_variance16x16_c, &aom_variance16x16_avx2),
        VarianceParam(16, 32, &aom_variance16x32_c, &aom_variance16x32_avx2),
        VarianceParam(16, 64, &aom_variance16x64_c, &aom_variance16x64_avx2),
        VarianceParam(32, 8, &aom_variance32x8_c, &aom_variance32x8_avx2),
        VarianceParam(32, 16, &aom_variance32x16_c, &aom_variance32x16_avx2),
        VarianceParam(32, 32, &aom_variance32x32_c, &aom_variance32x32_avx2),
        VarianceParam(32, 64, &aom_variance32x64_c, &aom_variance32x64_avx2),
        VarianceParam(64, 16, &aom_variance64x16_c, &aom_variance64x16_avx2),
        VarianceParam(64, 32, &aom_variance64x32_c, &aom_variance64x32_avx2),
        VarianceParam(64, 64, &aom_variance64x64_c, &aom_variance64x64_avx2),
        VarianceParam(64, 128, &aom_variance64x128_c, &aom_variance64x128_avx2),
        VarianceParam(128, 64, &aom_variance128x64_c, &aom_variance128x64_avx2),
        VarianceParam(128, 128, &aom_variance128x128_c,
                      &aom_variance128x128_avx2)));

}  // namespace
