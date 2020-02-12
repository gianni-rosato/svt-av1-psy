/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file subtract_avg_cfl_test.cc
 *
 * @brief Unit test for cfl subtract average function:
 * - subtract_average_avx2
 *
 * @author Cidana-Edmond
 *
 ******************************************************************************/

#include <sstream>
// workaround to eliminate the compiling warning on linux
// The macro will conflict with definition in gtest.h
#ifdef __USE_GNU
#undef __USE_GNU  // defined in EbThreads.h
#endif
#ifdef _GNU_SOURCE
#undef _GNU_SOURCE  // defined in EbThreads.h
#endif
#include "gtest/gtest.h"
#include "EbEncIntraPrediction.h"
#include "aom_dsp_rtcd.h"
#include "random.h"
#include "util.h"

/**
 * @brief Unit test for cfl subtract average function:
 * - subtract_average_avx2
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
 * data buffer: Fill with random values
 * TxSize: all the TxSize.
 */

namespace {

using svt_av1_test_tool::SVTRandom;

/** Defines of test parameters */
static TxSize TEST_PARAMS[] = {TX_4X4,
                               TX_4X8,
                               TX_4X16,
                               TX_8X4,
                               TX_8X8,
                               TX_8X16,
                               TX_8X32,
                               TX_16X4,
                               TX_16X8,
                               TX_16X16,
                               TX_16X32,
                               TX_32X8,
                               TX_32X16,
                               TX_32X32};

typedef CflSubtractAverageFn (*get_sub_avg_fn)(TxSize tx_size);

/** CFL_SUB_AVG_FN is a wrapper for subtract_average_avx2 to setup function
 * arguments easier, defines it to enable AVX2 subtract average functions */
CFL_SUB_AVG_FN(avx2)

class CflSubAvgTest : public ::testing::TestWithParam<TxSize> {
  public:
    CflSubAvgTest() : tx_size_(GetParam()) {
        width_ = tx_size_wide[tx_size_];
        height_ = tx_size_high[tx_size_];
        memset(data_tst_, 0, sizeof(data_tst_));
        memset(data_ref_, 0, sizeof(data_ref_));

        /** get test and reference function */
        sub_avg_tst_ = eb_get_subtract_average_fn_avx2(tx_size_);
        sub_avg_ref_ = eb_get_subtract_average_fn_c(tx_size_);
    }

    virtual ~CflSubAvgTest() {
        aom_clear_system_state();
    }

    void prepare_data() {
        SVTRandom rnd(15, false);
        for (uint32_t j = 0; j < height_; j++) {
            for (uint32_t i = 0; i < width_; i++) {
                int pos = j * CFL_BUF_LINE + i;
                data_tst_[pos] = data_ref_[pos] = data_origin_[pos] =
                    (int16_t)rnd.random();
            }
        }
    }

    void compare_data() {
        for (uint32_t j = 0; j < height_; j++) {
            for (uint32_t i = 0; i < width_; i++) {
                ASSERT_EQ(data_tst_[j * CFL_BUF_LINE + i],
                          data_ref_[j * CFL_BUF_LINE + i])
                    << print_data();
            }
        }
    }

    const std::string print_data() {
        std::string print_str;
        std::stringstream ss(print_str);
        ss << "test input data dump:\n";
        for (uint32_t j = 0; j < height_; j++) {
            for (uint32_t i = 0; i < width_; i++) {
                ss << data_origin_[j * CFL_BUF_LINE + i] << "\t";
            }
            ss << "\n";
        }
        return ss.str();
    }

    /** Verify this assembly code by comparing with reference c implementation.
     * Feed the same data and check test output and reference output.
     */
    void run_asm_compare_test(uint32_t test_times) {
        for (uint32_t it = 0; it < test_times; it++) {
            prepare_data();
            if (sub_avg_tst_)
                sub_avg_tst_(data_tst_);
            if (sub_avg_ref_)
                sub_avg_ref_(data_ref_);
            compare_data();
        }
    }

  protected:
    int16_t data_tst_[CFL_BUF_SQUARE]; /**< data buffer for test function*/
    int16_t data_ref_[CFL_BUF_SQUARE]; /**< data buffer for reference function*/
    int16_t data_origin_[CFL_BUF_SQUARE]; /**< data buffer for original data*/
    const TxSize tx_size_;                /**< tranform type of this test */
    uint32_t width_;                      /**< width of test data */
    uint32_t height_;                     /**< height of test data */
    CflSubtractAverageFn
        sub_avg_tst_; /**< asm subtract average function for test*/
    CflSubtractAverageFn
        sub_avg_ref_; /**< c subtract average function for reference*/
};

TEST_P(CflSubAvgTest, subtract_average_asm_test) {
    run_asm_compare_test(1000);
}

INSTANTIATE_TEST_CASE_P(AVX2, CflSubAvgTest, ::testing::ValuesIn(TEST_PARAMS));

}  // namespace
