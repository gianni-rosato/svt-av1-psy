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
 * @file HbdVarianceTest.cc
 *
 * @brief Unit test for HBD variance
 * functions:
 * - svt_aom_highbd_BD{8,10,12}_varianceW{8,16,32,64}xH{4,8,16,32,64}_sse2
 * - svt_aom_highbd_BD{8,10,12}_getS{8,16}xS{8,16}var_sse2
 *
 * @author  Cidana-Wenyao, Cidana-Edmond
 *
 ******************************************************************************/
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <new>
#include "gtest/gtest.h"
// workaround to eliminate the compiling warning on linux
// The macro will conflict with definition in gtest.h
#ifdef __USE_GNU
#undef __USE_GNU  // defined in EbThreads.h
#endif
#ifdef _GNU_SOURCE
#undef _GNU_SOURCE  // defined in EbThreads.h
#endif
#include "aom_dsp_rtcd.h"
#include "EbDefinitions.h"
#include "EbUtility.h"
#include "random.h"
#include "util.h"

using svt_av1_test_tool::SVTRandom;  // to generate the random
namespace {
#define MAX_BLOCK_SIZE (128 * 128)

using HighBdGetVarianceFunc = void (*)(const uint8_t* src8, int32_t src_stride,
                                       const uint8_t* ref8, int32_t ref_stride,
                                       uint32_t* sse, int32_t* sum);

using HighBdVarianceFunc = uint32_t (*)(const uint8_t* src8, int32_t src_stride,
                                        const uint8_t* ref8, int32_t ref_stride,
                                        uint32_t* sse);

// Truncate high bit depth results by downshifting (with rounding) by:
// 2 * (bit_depth - 8) for sse
// (bit_depth - 8) for se
static void round_hbd(const uint8_t bd, int64_t* se, uint64_t* sse) {
    switch (bd) {
    case 12:
        *sse = (*sse + 128) >> 8;
        *se = (*se + 8) >> 4;
        break;
    case 10:
        *sse = (*sse + 8) >> 4;
        *se = (*se + 2) >> 2;
        break;
    case 8:
    default: break;
    }
}

static void hbd_get_variance_ref(const uint32_t width, const uint32_t height,
                                 const uint8_t bd, const uint8_t* src8,
                                 int32_t src_stride, const uint8_t* ref8,
                                 int32_t ref_stride, uint32_t* sse,
                                 int32_t* sum) {
    int64_t sum_tmp = 0;
    uint64_t sse_tmp = 0;
    *sse = 0;
    for (uint32_t y = 0; y < height; y++) {
        for (uint32_t x = 0; x < width; x++) {
            int diff = CONVERT_TO_SHORTPTR(src8)[y * src_stride + x] -
                       CONVERT_TO_SHORTPTR(ref8)[y * ref_stride + x];
            sum_tmp += diff;
            sse_tmp += diff * diff;
        }
    }
    round_hbd(bd, &sum_tmp, &sse_tmp);
    *sse = static_cast<uint32_t>(sse_tmp);
    *sum = static_cast<int32_t>(sum_tmp);
}

static uint32_t hbd_variance_ref(const uint32_t width, const uint32_t height,
                                 const uint8_t bd, const uint8_t* src8,
                                 int32_t src_stride, const uint8_t* ref8,
                                 int32_t ref_stride, uint32_t* sse) {
    int32_t sum = 0;
    hbd_get_variance_ref(
        width, height, bd, src8, src_stride, ref8, ref_stride, sse, &sum);
    return static_cast<uint32_t>(
        *sse - ((((int64_t)sum * sum)) / ((int64_t)width * height)));
}

// High bit-depth variance test
using HbdVarianceParam = std::tuple<uint32_t,            /**< width */
                                    uint32_t,            /**< height */
                                    uint32_t,            /**< bit-depth */
                                    HighBdVarianceFunc>; /**< test function */

extern "C" {
unsigned int svt_aom_highbd_8_variance8x8_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_8_variance8x16_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_8_variance8x32_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_8_variance16x4_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_8_variance16x8_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_8_variance16x16_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_8_variance16x32_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_8_variance16x64_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_8_variance32x8_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_8_variance32x16_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_8_variance32x32_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_8_variance32x64_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_8_variance64x16_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_8_variance64x32_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_8_variance64x64_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_10_variance8x8_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_10_variance8x16_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_10_variance8x32_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_10_variance16x4_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_10_variance16x8_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_10_variance16x16_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_10_variance16x32_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_10_variance16x64_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_10_variance32x8_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_10_variance32x16_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_10_variance32x32_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_10_variance32x64_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_10_variance64x16_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_10_variance64x32_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_10_variance64x64_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_10_variance8x8_avx2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_10_variance8x16_avx2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_10_variance8x32_avx2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_10_variance16x8_avx2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_10_variance16x16_avx2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_10_variance16x32_avx2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_10_variance16x64_avx2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_10_variance32x8_avx2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_10_variance32x16_avx2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_10_variance32x32_avx2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_10_variance32x64_avx2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_10_variance64x16_avx2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_10_variance64x32_avx2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_10_variance64x64_avx2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_10_variance64x128_avx2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_10_variance128x64_avx2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_10_variance128x128_avx2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_12_variance8x8_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_12_variance8x16_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_12_variance8x32_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_12_variance16x4_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_12_variance16x8_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_12_variance16x16_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_12_variance16x32_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_12_variance16x64_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_12_variance32x8_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_12_variance32x16_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_12_variance32x32_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_12_variance32x64_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_12_variance64x16_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_12_variance64x32_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
unsigned int svt_aom_highbd_12_variance64x64_sse2(const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, int ref_stride, unsigned int *sse);
};

/**
 * @brief Unit test for HBD variance
 * functions:
 * - svt_aom_highbd_BD{8,10,12}_varianceW{8,16,32,64}xH{4,8,16,32,64}_sse2
 *
 * Test strategy:
 *  This test case use random source, max source, zero source as test
 * pattern.
 *
 *
 * Expected result:
 *  Results come from reference function and target function are
 * equal.
 *
 * Test cases:
 * - ZeroTest
 * - MaximumTest
 * - MatchTest
 */
class HbdVarianceTest : public ::testing::TestWithParam<HbdVarianceParam> {
  public:
    HbdVarianceTest()
        : rnd_(16, false),
          width_(TEST_GET_PARAM(0)),
          height_(TEST_GET_PARAM(1)),
          bd_(TEST_GET_PARAM(2)),
          tst_func_(TEST_GET_PARAM(3)) {
        src_data_ = reinterpret_cast<uint16_t*>(
            svt_aom_memalign(32, 2 * MAX_BLOCK_SIZE));
        ref_data_ = reinterpret_cast<uint16_t*>(
            svt_aom_memalign(32, 2 * MAX_BLOCK_SIZE));
    }

    ~HbdVarianceTest() {
        svt_aom_free(src_data_);
        src_data_ = nullptr;
        svt_aom_free(ref_data_);
        ref_data_ = nullptr;
    }

    void run_zero_test(int times) {
        for (int i = 0; i < times; ++i) {
            for (int j = 0; j < MAX_BLOCK_SIZE; ++j) {
                src_data_[j] = rnd_.random() & ((1 << bd_) - 1);
                ref_data_[j] = src_data_[j];
            }
            uint32_t sse_tst = 0;
            uint32_t var_tst = tst_func_(CONVERT_TO_BYTEPTR(src_data_),
                                         width_,
                                         CONVERT_TO_BYTEPTR(ref_data_),
                                         width_,
                                         &sse_tst);
            ASSERT_EQ(var_tst, 0u) << "Expect 0 variance, got: " << var_tst;
        }
    }

    void run_maximum_test() {
        for (int j = 0; j < MAX_BLOCK_SIZE; ++j) {
            src_data_[j] = 0;
            ref_data_[j] = (1 << bd_) - 1;
        }
        uint32_t sse_tst = 0, sse_ref = 0;
        uint32_t var_tst = tst_func_(CONVERT_TO_BYTEPTR(src_data_),
                                     width_,
                                     CONVERT_TO_BYTEPTR(ref_data_),
                                     width_,
                                     &sse_tst);
        uint32_t var_ref = hbd_variance_ref(width_,
                                            height_,
                                            bd_,
                                            CONVERT_TO_BYTEPTR(src_data_),
                                            width_,
                                            CONVERT_TO_BYTEPTR(ref_data_),
                                            width_,
                                            &sse_ref);
        ASSERT_EQ(var_tst, var_ref)
            << "Expect var " << var_ref << " got " << var_tst
            << " size: "<< width_ << "x" << height_;
        ASSERT_EQ(sse_tst, sse_ref)
            << "Expect sse " << sse_ref << " got " << sse_tst
            << " size: "<< width_ << "x" << height_;
    }

    void run_match_test(int times) {
        for (int i = 0; i < times; ++i) {
            for (int j = 0; j < MAX_BLOCK_SIZE; ++j) {
                src_data_[j] = rnd_.random() & ((1 << bd_) - 1);
                ref_data_[j] = rnd_.random() & ((1 << bd_) - 1);
            }
            uint32_t sse_tst = 0, sse_ref = 0;
            uint32_t var_tst = tst_func_(CONVERT_TO_BYTEPTR(src_data_),
                                         width_,
                                         CONVERT_TO_BYTEPTR(ref_data_),
                                         width_,
                                         &sse_tst);
            uint32_t var_ref = hbd_variance_ref(width_,
                                                height_,
                                                bd_,
                                                CONVERT_TO_BYTEPTR(src_data_),
                                                width_,
                                                CONVERT_TO_BYTEPTR(ref_data_),
                                                width_,
                                                &sse_ref);
            ASSERT_EQ(var_tst, var_ref)
                << "Error at variance test index: " << i << " size: "<< width_ << "x" << height_;
            ASSERT_EQ(sse_tst, sse_ref) << "Error at sse test index: " << i
                << " size: "<< width_ << "x" << height_;
        }
    }

  private:
    SVTRandom rnd_;
    uint16_t* src_data_;
    uint16_t* ref_data_;
    uint32_t width_;
    uint32_t height_;
    uint32_t bd_;
    HighBdVarianceFunc tst_func_;
};

TEST_P(HbdVarianceTest, ZeroTest) {
    run_zero_test(10);
};

TEST_P(HbdVarianceTest, MaximumTest) {
    run_maximum_test();
};

TEST_P(HbdVarianceTest, MatchTest) {
    run_match_test(10);
};


static const HbdVarianceParam HbdTestVector[] = {
    HbdVarianceParam(  8,  8, 8, svt_aom_highbd_8_variance8x8_sse2),
    HbdVarianceParam(  8, 16, 8, svt_aom_highbd_8_variance8x16_sse2),
    HbdVarianceParam(  8, 32, 8, svt_aom_highbd_8_variance8x32_sse2),
    HbdVarianceParam( 16,  4, 8, svt_aom_highbd_8_variance16x4_sse2),
    HbdVarianceParam( 16,  8, 8, svt_aom_highbd_8_variance16x8_sse2),
    HbdVarianceParam( 16, 16, 8, svt_aom_highbd_8_variance16x16_sse2),
    HbdVarianceParam( 16, 32, 8, svt_aom_highbd_8_variance16x32_sse2),
    HbdVarianceParam( 16, 64, 8, svt_aom_highbd_8_variance16x64_sse2),
    HbdVarianceParam( 32,  8, 8, svt_aom_highbd_8_variance32x8_sse2),
    HbdVarianceParam( 32, 16, 8, svt_aom_highbd_8_variance32x16_sse2),
    HbdVarianceParam( 32, 32, 8, svt_aom_highbd_8_variance32x32_sse2),
    HbdVarianceParam( 32, 64, 8, svt_aom_highbd_8_variance32x64_sse2),
    HbdVarianceParam( 64, 16, 8, svt_aom_highbd_8_variance64x16_sse2),
    HbdVarianceParam( 64, 32, 8, svt_aom_highbd_8_variance64x32_sse2),
    HbdVarianceParam( 64, 64, 8, svt_aom_highbd_8_variance64x64_sse2),
    HbdVarianceParam(  8,  8, 10, svt_aom_highbd_10_variance8x8_sse2),
    HbdVarianceParam(  8, 16, 10, svt_aom_highbd_10_variance8x16_sse2),
    HbdVarianceParam(  8, 32, 10, svt_aom_highbd_10_variance8x32_sse2),
    HbdVarianceParam( 16,  4, 10, svt_aom_highbd_10_variance16x4_sse2),
    HbdVarianceParam( 16,  8, 10, svt_aom_highbd_10_variance16x8_sse2),
    HbdVarianceParam( 16, 16, 10, svt_aom_highbd_10_variance16x16_sse2),
    HbdVarianceParam( 16, 32, 10, svt_aom_highbd_10_variance16x32_sse2),
    HbdVarianceParam( 16, 64, 10, svt_aom_highbd_10_variance16x64_sse2),
    HbdVarianceParam( 32,  8, 10, svt_aom_highbd_10_variance32x8_sse2),
    HbdVarianceParam( 32, 16, 10, svt_aom_highbd_10_variance32x16_sse2),
    HbdVarianceParam( 32, 32, 10, svt_aom_highbd_10_variance32x32_sse2),
    HbdVarianceParam( 32, 64, 10, svt_aom_highbd_10_variance32x64_sse2),
    HbdVarianceParam( 64, 16, 10, svt_aom_highbd_10_variance64x16_sse2),
    HbdVarianceParam( 64, 32, 10, svt_aom_highbd_10_variance64x32_sse2),
    HbdVarianceParam( 64, 64, 10, svt_aom_highbd_10_variance64x64_sse2),
    HbdVarianceParam(  8,  8, 10, svt_aom_highbd_10_variance8x8_avx2),
    HbdVarianceParam(  8, 16, 10, svt_aom_highbd_10_variance8x16_avx2),
    HbdVarianceParam(  8, 32, 10, svt_aom_highbd_10_variance8x32_avx2),
    HbdVarianceParam( 16,  8, 10, svt_aom_highbd_10_variance16x8_avx2),
    HbdVarianceParam( 16, 16, 10, svt_aom_highbd_10_variance16x16_avx2),
    HbdVarianceParam( 16, 32, 10, svt_aom_highbd_10_variance16x32_avx2),
    HbdVarianceParam( 16, 64, 10, svt_aom_highbd_10_variance16x64_avx2),
    HbdVarianceParam( 32,  8, 10, svt_aom_highbd_10_variance32x8_avx2),
    HbdVarianceParam( 32, 16, 10, svt_aom_highbd_10_variance32x16_avx2),
    HbdVarianceParam( 32, 32, 10, svt_aom_highbd_10_variance32x32_avx2),
    HbdVarianceParam( 32, 64, 10, svt_aom_highbd_10_variance32x64_avx2),
    HbdVarianceParam( 64, 16, 10, svt_aom_highbd_10_variance64x16_avx2),
    HbdVarianceParam( 64, 32, 10, svt_aom_highbd_10_variance64x32_avx2),
    HbdVarianceParam( 64, 64, 10, svt_aom_highbd_10_variance64x64_avx2),
    HbdVarianceParam( 64,128, 10, svt_aom_highbd_10_variance64x128_avx2),
    HbdVarianceParam(128, 64, 10, svt_aom_highbd_10_variance128x64_avx2),
    HbdVarianceParam(128,128, 10, svt_aom_highbd_10_variance128x128_avx2),
    HbdVarianceParam(  8,  8, 12, svt_aom_highbd_12_variance8x8_sse2),
    HbdVarianceParam(  8, 16, 12, svt_aom_highbd_12_variance8x16_sse2),
    HbdVarianceParam(  8, 32, 12, svt_aom_highbd_12_variance8x32_sse2),
    HbdVarianceParam( 16,  4, 12, svt_aom_highbd_12_variance16x4_sse2),
    HbdVarianceParam( 16,  8, 12, svt_aom_highbd_12_variance16x8_sse2),
    HbdVarianceParam( 16, 16, 12, svt_aom_highbd_12_variance16x16_sse2),
    HbdVarianceParam( 16, 32, 12, svt_aom_highbd_12_variance16x32_sse2),
    HbdVarianceParam( 16, 64, 12, svt_aom_highbd_12_variance16x64_sse2),
    HbdVarianceParam( 32,  8, 12, svt_aom_highbd_12_variance32x8_sse2),
    HbdVarianceParam( 32, 16, 12, svt_aom_highbd_12_variance32x16_sse2),
    HbdVarianceParam( 32, 32, 12, svt_aom_highbd_12_variance32x32_sse2),
    HbdVarianceParam( 32, 64, 12, svt_aom_highbd_12_variance32x64_sse2),
    HbdVarianceParam( 64, 16, 12, svt_aom_highbd_12_variance64x16_sse2),
    HbdVarianceParam( 64, 32, 12, svt_aom_highbd_12_variance64x32_sse2),
    HbdVarianceParam( 64, 64, 12, svt_aom_highbd_12_variance64x64_sse2),
};

INSTANTIATE_TEST_CASE_P(Variance, HbdVarianceTest,
                        ::testing::ValuesIn(HbdTestVector));
#define svt_aom_highbd_8_get8x8var_sse2 svt_aom_highbd_get8x8var_sse2
#define svt_aom_highbd_8_get16x16var_sse2 svt_aom_highbd_get16x16var_sse2
#define GEN_GET_VAR_FUNC_(S, bd) svt_aom_highbd_##bd##_get##S##x##S##var_sse2
#define GEN_SQUARE_VAR_PARAM_(S, bd) \
    HbdSquareVarianceParam(S, bd, GEN_GET_VAR_FUNC_(S, bd))
#define GEN_HBD_SQUARE_VAR_PARAM(S)                            \
    GEN_SQUARE_VAR_PARAM_(S, 8), GEN_SQUARE_VAR_PARAM_(S, 10), \
        GEN_SQUARE_VAR_PARAM_(S, 12)
#define EXTERN_HBD_SQUARE_VAR_FUNC(S)                        \
    extern "C" void GEN_GET_VAR_FUNC_(S, 8)(const uint8_t*,  \
                                            int32_t,         \
                                            const uint8_t*,  \
                                            int32_t,         \
                                            uint32_t*,       \
                                            int32_t*);       \
    extern "C" void GEN_GET_VAR_FUNC_(S, 10)(const uint8_t*, \
                                             int32_t,        \
                                             const uint8_t*, \
                                             int32_t,        \
                                             uint32_t*,      \
                                             int32_t*);      \
    extern "C" void GEN_GET_VAR_FUNC_(S, 12)(const uint8_t*, \
                                             int32_t,        \
                                             const uint8_t*, \
                                             int32_t,        \
                                             uint32_t*,      \
                                             int32_t*);

EXTERN_HBD_SQUARE_VAR_FUNC(8);
EXTERN_HBD_SQUARE_VAR_FUNC(16);

// High bit-depth square variance test
using HbdSquareVarianceParam =
    std::tuple<uint32_t,               /**< square length */
               uint32_t,               /**< bit-depth */
               HighBdGetVarianceFunc>; /**< test function */

/**
 * @brief Unit test for HBD variance
 * functions:
 * - svt_aom_highbd_BD{8,10,12}_getS{8,16}xS{8,16}var_sse2
 *
 * Test strategy:
 *  This test case use random source, max source, zero source as test
 * pattern.
 *
 *
 * Expected result:
 *  Results come from reference function and target function are
 * equal.
 *
 * Test cases:
 * - ZeroTest
 * - MaximumTest
 * - MatchTest
 */
class HbdSquareVarianceTest
    : public ::testing::TestWithParam<HbdSquareVarianceParam> {
  public:
    HbdSquareVarianceTest()
        : rnd_(16, false),
          length_(TEST_GET_PARAM(0)),
          bd_(TEST_GET_PARAM(1)),
          tst_func_(TEST_GET_PARAM(2)) {
        src_data_ = reinterpret_cast<uint16_t*>(
            svt_aom_memalign(32, 2 * MAX_BLOCK_SIZE));
        ref_data_ = reinterpret_cast<uint16_t*>(
            svt_aom_memalign(32, 2 * MAX_BLOCK_SIZE));
    }

    ~HbdSquareVarianceTest() {
        svt_aom_free(src_data_);
        src_data_ = nullptr;
        svt_aom_free(ref_data_);
        ref_data_ = nullptr;
    }

    void run_zero_test(int times) {
        for (int i = 0; i < times; ++i) {
            for (int j = 0; j < MAX_BLOCK_SIZE; ++j) {
                src_data_[j] = rnd_.random() & ((1 << bd_) - 1);
                ref_data_[j] = src_data_[j];
            }
            int32_t sum_tst = 0;
            uint32_t sse_tst = 0;
            tst_func_(CONVERT_TO_BYTEPTR(src_data_),
                      length_,
                      CONVERT_TO_BYTEPTR(ref_data_),
                      length_,
                      &sse_tst,
                      &sum_tst);
            ASSERT_EQ(sse_tst, 0u) << "Expect 0 sse, got: " << sse_tst;
            ASSERT_EQ(sum_tst, 0) << "Expect 0 sum, got: " << sum_tst;
        }
    }

    void run_maximum_test() {
        for (int j = 0; j < MAX_BLOCK_SIZE; ++j) {
            src_data_[j] = 0;
            ref_data_[j] = (1 << bd_) - 1;
        }
        uint32_t sse_tst = 0, sse_ref = 0;
        int32_t sum_tst = 0, sum_ref = 0;
        tst_func_(CONVERT_TO_BYTEPTR(src_data_),
                  length_,
                  CONVERT_TO_BYTEPTR(ref_data_),
                  length_,
                  &sse_tst,
                  &sum_tst);
        hbd_get_variance_ref(length_,
                             length_,
                             bd_,
                             CONVERT_TO_BYTEPTR(src_data_),
                             length_,
                             CONVERT_TO_BYTEPTR(ref_data_),
                             length_,
                             &sse_ref,
                             &sum_ref);
        ASSERT_EQ(sum_tst, sum_ref) << "Error at see in variance test";
        ASSERT_EQ(sse_tst, sse_ref) << "Error at error sum in variance test";
    }

    void run_match_test(int times) {
        for (int i = 0; i < times; ++i) {
            for (int j = 0; j < MAX_BLOCK_SIZE; ++j) {
                src_data_[j] = rnd_.random() & ((1 << bd_) - 1);
                ref_data_[j] = rnd_.random() & ((1 << bd_) - 1);
            }
            uint32_t sse_tst = 0, sse_ref = 0;
            int32_t sum_tst = 0, sum_ref = 0;
            tst_func_(CONVERT_TO_BYTEPTR(src_data_),
                      length_,
                      CONVERT_TO_BYTEPTR(ref_data_),
                      length_,
                      &sse_tst,
                      &sum_tst);
            hbd_get_variance_ref(length_,
                                 length_,
                                 bd_,
                                 CONVERT_TO_BYTEPTR(src_data_),
                                 length_,
                                 CONVERT_TO_BYTEPTR(ref_data_),
                                 length_,
                                 &sse_ref,
                                 &sum_ref);
            ASSERT_EQ(sum_tst, sum_ref) << "Error at error sum index: " << i;
            ASSERT_EQ(sse_tst, sse_ref) << "Error at sse index: " << i;
        }
    }

  private:
    SVTRandom rnd_;
    uint16_t* src_data_;
    uint16_t* ref_data_;
    uint32_t length_;
    uint32_t bd_;
    HighBdGetVarianceFunc tst_func_;
};

TEST_P(HbdSquareVarianceTest, ZeroTest) {
    run_zero_test(10);
};

TEST_P(HbdSquareVarianceTest, MaximumTest) {
    run_maximum_test();
};

TEST_P(HbdSquareVarianceTest, MatchTest) {
    run_match_test(10);
};

static const HbdSquareVarianceParam HbdSquareVarTestVector[] = {
    GEN_HBD_SQUARE_VAR_PARAM(8),
    GEN_HBD_SQUARE_VAR_PARAM(16),
};

INSTANTIATE_TEST_CASE_P(Variance, HbdSquareVarianceTest,
                        ::testing::ValuesIn(HbdSquareVarTestVector));


/**
 * @brief Unit test for different implementation of HBD variance with size 16x16 and 32x32
 * functions:
 * - variance_highbd_c
 * - variance_highbd_avx2
 *
 * Test strategy:
 *  This test case use random source, max source, zero source as test
 * pattern.
 *
 * Expected result:
 *  Results come from reference function and target function are
 * equal.
 *
 * Test cases:
 * - ZeroTest
 * - MaximumTest
 * - MatchTest
 *
 * @author  intel tszumski
 *
 */
using HighBdVarianceNoRoundFunc = uint32_t (*)(const uint16_t* src,
                                               int32_t src_stride,
                                               const uint16_t* ref,
                                               int32_t ref_stride, int w,
                                               int h, uint32_t* sse);

using HbdSquareVarianceNoRoundParam =
    std::tuple<uint32_t,                   /**< square length */
               HighBdVarianceNoRoundFunc>; /**< test function */

class HbdSquareVarianceNoRoundTest
    : public ::testing::TestWithParam<HbdSquareVarianceNoRoundParam> {
  public:
    HbdSquareVarianceNoRoundTest()
        : rnd_(16, false),
          bd_(10),
          length_(TEST_GET_PARAM(0)),
          tst_func_(TEST_GET_PARAM(1)) {
        src_data_ = reinterpret_cast<uint16_t*>(
            svt_aom_memalign(32, 2 * MAX_BLOCK_SIZE));
        ref_data_ = reinterpret_cast<uint16_t*>(
            svt_aom_memalign(32, 2 * MAX_BLOCK_SIZE));
    }

    ~HbdSquareVarianceNoRoundTest() {
        svt_aom_free(src_data_);
        src_data_ = nullptr;
        svt_aom_free(ref_data_);
        ref_data_ = nullptr;
    }

    void run_zero_test(int times) {
        for (int i = 0; i < times; ++i) {
            for (int j = 0; j < MAX_BLOCK_SIZE; ++j) {
                src_data_[j] = rnd_.random() & ((1 << bd_) - 1);
                ref_data_[j] = src_data_[j];
            }
            int32_t distortion_tst = 0;
            uint32_t sse_tst = 0;
            distortion_tst = tst_func_(src_data_,
                      length_,
                      ref_data_,
                      length_,
                      length_,
                      length_,
                      &sse_tst);
            ASSERT_EQ(sse_tst, 0u) << "Expect 0 sse, got: " << sse_tst;
            ASSERT_EQ(distortion_tst, 0)
                << "Expect 0 distortion, got: " << distortion_tst;
        }
    }

    void run_maximum_test() {
        for (int j = 0; j < MAX_BLOCK_SIZE; ++j) {
            src_data_[j] = 0;
            ref_data_[j] = (1 << bd_) - 1;
        }

        int32_t distortion_tst = 0, distortion_ref = 0;
        uint32_t sse_tst = 0, sse_ref = 0;

        distortion_tst = tst_func_(
            src_data_, length_, ref_data_, length_, length_, length_, &sse_tst);

        distortion_ref = variance_highbd_c(
            src_data_, length_, ref_data_, length_, length_, length_, &sse_ref);
        ASSERT_EQ(distortion_tst, distortion_ref)
            << "Error at distortion in variance test";
        ASSERT_EQ(sse_tst, sse_ref) << "Error at error sse in variance test";
    }

    void run_match_test(int times) {
        for (int i = 0; i < times; ++i) {
            for (int j = 0; j < MAX_BLOCK_SIZE; ++j) {
                src_data_[j] = rnd_.random() & ((1 << bd_) - 1);
                ref_data_[j] = rnd_.random() & ((1 << bd_) - 1);
            }
            int32_t distortion_tst = 0, distortion_ref = 0;
            uint32_t sse_tst = 0, sse_ref = 0;

            distortion_tst = tst_func_(src_data_,
                                       length_,
                                       ref_data_,
                                       length_,
                                       length_,
                                       length_,
                                       &sse_tst);

            distortion_ref = variance_highbd_c(src_data_,
                                               length_,
                                               ref_data_,
                                               length_,
                                               length_,
                                               length_,
                                               &sse_ref);
            ASSERT_EQ(distortion_tst, distortion_ref)
                << "Error at distortion in variance test";
            ASSERT_EQ(sse_tst, sse_ref)
                << "Error at error sse in variance test";
        }
    }

  private:
    SVTRandom rnd_;
    uint32_t bd_;
    uint32_t length_;
    HighBdVarianceNoRoundFunc tst_func_;
    uint16_t* src_data_;
    uint16_t* ref_data_;
};

TEST_P(HbdSquareVarianceNoRoundTest, ZeroTest) {
    run_zero_test(10);
};

TEST_P(HbdSquareVarianceNoRoundTest, MaximumTest) {
    run_maximum_test();
};

TEST_P(HbdSquareVarianceNoRoundTest, MatchTest) {
    run_match_test(10);
};

INSTANTIATE_TEST_CASE_P(
    Variance, HbdSquareVarianceNoRoundTest,
    ::testing::Combine(::testing::Values(16, 32),
                       ::testing::Values(variance_highbd_avx2)));

}  // namespace
