/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file FilterIntraPredTest.cc
 *
 * @brief Unit test for filter intra predictor functions:
 * - eb_av1_filter_intra_predictor_sse4_1
 *
 * @author Cidana-Edmond
 *
 ******************************************************************************/
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
#include "random.h"
#include "util.h"
#include "EbUtility.h"

using std::tuple;
using svt_av1_test_tool::SVTRandom;  // to generate the random

namespace {
using PredParams = tuple<FilterIntraMode, TxSize>;

class FilterIntraPredTest : public ::testing::TestWithParam<PredParams> {
  public:
    FilterIntraPredTest()
        : pred_mode_(TEST_GET_PARAM(0)),
          tx_size_(TEST_GET_PARAM(1)),
          rnd_(8, false) {
        input_ = reinterpret_cast<uint8_t*>(
            eb_aom_memalign(32, 2 * MAX_TX_SIZE + 1));
        pred_ref_ =
            reinterpret_cast<uint8_t*>(eb_aom_memalign(32, MAX_TX_SQUARE));
        pred_tst_ =
            reinterpret_cast<uint8_t*>(eb_aom_memalign(32, MAX_TX_SQUARE));
    }
    virtual ~FilterIntraPredTest() {
        eb_aom_free(input_);
        eb_aom_free(pred_ref_);
        eb_aom_free(pred_tst_);
        aom_clear_system_state();
    }

  protected:
    void run_test(size_t test_num) {
        const ptrdiff_t stride = tx_size_wide[tx_size_];
        const uint8_t* left = input_;
        const uint8_t* above = input_ + MAX_TX_SIZE;
        for (size_t i = 0; i < test_num; i++) {
            prepare_data();
            eb_av1_filter_intra_predictor_c(
                pred_ref_, stride, tx_size_, &above[1], left, pred_mode_);
            eb_av1_filter_intra_predictor_sse4_1(
                pred_tst_, stride, tx_size_, &above[1], left, pred_mode_);
            check_output(i);
        }
    }

    void prepare_data() {
        for (size_t i = 0; i < 2 * MAX_TX_SIZE + 1; i++)
            input_[i] = rnd_.random();
    }

    void check_output(size_t test_num) {
        for (int32_t i = 0; i < tx_size_wide[tx_size_] * tx_size_high[tx_size_];
             i++) {
            ASSERT_EQ(pred_tst_[i], pred_ref_[i])
                << "Error at position: " << i << " "
                << "Tx size: " << tx_size_wide[tx_size_] << "x"
                << tx_size_high[tx_size_] << " "
                << "Test number: " << test_num;
        }
    }

  protected:
    const FilterIntraMode pred_mode_;
    const TxSize tx_size_;
    uint8_t* input_;
    uint8_t* pred_tst_;
    uint8_t* pred_ref_;
    SVTRandom rnd_;
};

TEST_P(FilterIntraPredTest, RunCheckOutput) {
    run_test(1000);
}

const FilterIntraMode PRED_MODE_TABLE[] = {FILTER_DC_PRED,
                                           FILTER_V_PRED,
                                           FILTER_H_PRED,
                                           FILTER_D157_PRED,
                                           FILTER_PAETH_PRED};

const TxSize TX_SIZE_TABLE[] = {TX_4X4,
                                TX_8X8,
                                TX_16X16,
                                TX_32X32,
                                TX_4X8,
                                TX_8X4,
                                TX_8X16,
                                TX_16X8,
                                TX_16X32,
                                TX_32X16,
                                TX_4X16,
                                TX_16X4,
                                TX_8X32,
                                TX_32X8};

INSTANTIATE_TEST_CASE_P(AV1, FilterIntraPredTest,
                        ::testing::Combine(::testing::ValuesIn(PRED_MODE_TABLE),
                                           ::testing::ValuesIn(TX_SIZE_TABLE)));
}  // namespace
