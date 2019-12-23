/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file QuantAsmTest.c
 *
 * @brief Unit test for quantize avx2 functions:
 * - eb_aom_highbd_quantize_b_avx2
 * - eb_aom_highbd_quantize_b_32x32_avx2
 * - eb_aom_highbd_quantize_b_64x64_avx2
 *
 * @author Cidana-Zhengwen
 *
 ******************************************************************************/

#include <random>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "gtest/gtest.h"

// Workaround to eliminate the compiling warning on linux
// The macro will conflict with definition in gtest.h
#ifdef __USE_GNU
#undef __USE_GNU  // defined in EbThreads.h
#endif
#ifdef _GNU_SOURCE
#undef _GNU_SOURCE  // defined in EbThreads.h
#endif

#include "EbDefinitions.h"
#include "EbTransforms.h"
#include "EbPictureControlSet.h"
#include "aom_dsp_rtcd.h"
#include "util.h"
#include "random.h"

namespace QuantizeAsmTest {
extern "C" void eb_av1_build_quantizer(AomBitDepth bit_depth, int32_t y_dc_delta_q,
                                    int32_t u_dc_delta_q, int32_t u_ac_delta_q,
                                    int32_t v_dc_delta_q, int32_t v_ac_delta_q,
                                    Quants *const quants, Dequants *const deq);

using QuantizeFunc = void (*)(const TranLow *coeff_ptr, intptr_t n_coeffs,
                              int32_t skip_block, const int16_t *zbin_ptr,
                              const int16_t *round_ptr,
                              const int16_t *quant_ptr,
                              const int16_t *quant_shift_ptr,
                              TranLow *qcoeff_ptr, TranLow *dqcoeff_ptr,
                              const int16_t *dequant_ptr, uint16_t *eob_ptr,
                              const int16_t *scan, const int16_t *iscan);

using QuantizeParam = std::tuple<int, int>;

using svt_av1_test_tool::SVTRandom;  // to generate the random
/**
 * @brief Unit test for quantize avx2 functions:
 * - eb_aom_highbd_quantize_b_avx2
 * - eb_aom_highbd_quantize_b_32x32_avx2
 * - eb_aom_highbd_quantize_b_64x64_avx2
 *
 * Test strategy:
 * These tests use quantize C function as reference, input the same data and
 * compare C function output with avx2 function output, so as to check
 * quant/dequant/eob output consistancy of quantize avx2 implementation.
 *
 * Expect result:
 * avx2 output quant/dequant/eob should be exactly same as C output.
 *
 * Test coverage:
 * All combinations of all tx_size with 8bit/10bit.
 *
 * Test cases:
 * - AVX2/QuantizeBTest.input_zero_all
 * - AVX2/QuantizeBTest.input_dcac_minmax_q_n
 * - AVX2/QuantizeBTest.input_random_dc_only
 * - AVX2/QuantizeBTest.input_random_all_q_all
 */
class QuantizeBTest : public ::testing::TestWithParam<QuantizeParam> {
  protected:
    QuantizeBTest()
        : tx_size_(static_cast<TxSize>(TEST_GET_PARAM(0))),
          bd_(static_cast<AomBitDepth>(TEST_GET_PARAM(1))) {
        n_coeffs_ = av1_get_max_eob(tx_size_);
        coeff_min_ = -(1 << (7 + bd_));
        coeff_max_ = (1 << (7 + bd_)) - 1;
        rnd_ = new SVTRandom(coeff_min_, coeff_max_);

        eb_av1_build_quantizer(bd_, 0, 0, 0, 0, 0, &qtab_quants_, &qtab_deq_);
        setup_func_ptrs();
    }

    virtual ~QuantizeBTest() {
        delete rnd_;
        aom_clear_system_state();
    }

    void SetUp() override {
        coeff_in_ = reinterpret_cast<TranLow *>(
            eb_aom_memalign(32, MAX_TX_SQUARE * sizeof(TranLow)));
        qcoeff_ref_ = reinterpret_cast<TranLow *>(
            eb_aom_memalign(32, MAX_TX_SQUARE * sizeof(TranLow)));
        dqcoeff_ref_ = reinterpret_cast<TranLow *>(
            eb_aom_memalign(32, MAX_TX_SQUARE * sizeof(TranLow)));
        qcoeff_test_ = reinterpret_cast<TranLow *>(
            eb_aom_memalign(32, MAX_TX_SQUARE * sizeof(TranLow)));
        dqcoeff_test_ = reinterpret_cast<TranLow *>(
            eb_aom_memalign(32, MAX_TX_SQUARE * sizeof(TranLow)));
    }

    void TearDown() override {
        eb_aom_free(coeff_in_);
        eb_aom_free(qcoeff_ref_);
        eb_aom_free(dqcoeff_ref_);
        eb_aom_free(qcoeff_test_);
        eb_aom_free(dqcoeff_test_);
        aom_clear_system_state();
    }

    /*
     * @brief setup reference and target function ptrs
     * @see setup_rtcd_internal() in aom_dsp_rtcd.h
     */
    void setup_func_ptrs() {
        if (bd_ == AOM_BITS_8) {
            if (tx_size_ == TX_32X32) {
                quant_ref_ = eb_aom_quantize_b_32x32_c_ii;
                quant_test_ = eb_aom_quantize_b_32x32_avx2;
            } else if (tx_size_ == TX_64X64) {
                quant_ref_ = eb_aom_quantize_b_64x64_c_ii;
                quant_test_ = eb_aom_quantize_b_64x64_avx2;
            } else {
                quant_ref_ = eb_aom_quantize_b_c_ii;
                quant_test_ = eb_aom_quantize_b_avx2;
            }
        } else {
            if (tx_size_ == TX_32X32) {
                quant_ref_ = eb_aom_highbd_quantize_b_32x32_c;
                quant_test_ = eb_aom_highbd_quantize_b_32x32_avx2;
            } else if (tx_size_ == TX_64X64) {
                quant_ref_ = eb_aom_highbd_quantize_b_64x64_c;
                quant_test_ = eb_aom_highbd_quantize_b_64x64_avx2;
            } else {
                quant_ref_ = eb_aom_highbd_quantize_b_c;
                quant_test_ = eb_aom_highbd_quantize_b_avx2;
            }
        }
    }

    /*
     * @brief run quantize with the same input data, then check whether output
     * quant/dequant/eob are exactly same between C output and avx2 outptu
     */
    void run_quantize(int q) {
        const int32_t skip_block = 0;
        const ScanOrder *const sc = &av1_scan_orders[tx_size_][DCT_DCT];
        const int16_t *zbin = qtab_quants_.y_zbin[q];
        const int16_t *round = qtab_quants_.y_round[q];
        const int16_t *quant = qtab_quants_.y_quant[q];
        const int16_t *quant_shift = qtab_quants_.y_quant_shift[q];
        const int16_t *dequant = qtab_deq_.y_dequant_qtx[q];

        memset(qcoeff_ref_, 0, MAX_TX_SQUARE * sizeof(TranLow));
        memset(dqcoeff_ref_, 0, MAX_TX_SQUARE * sizeof(TranLow));
        memset(qcoeff_test_, 0, MAX_TX_SQUARE * sizeof(TranLow));
        memset(dqcoeff_test_, 0, MAX_TX_SQUARE * sizeof(TranLow));

        quant_ref_(coeff_in_,
                   n_coeffs_,
                   skip_block,
                   zbin,
                   round,
                   quant,
                   quant_shift,
                   qcoeff_ref_,
                   dqcoeff_ref_,
                   dequant,
                   &eob_ref_,
                   sc->scan,
                   sc->iscan);

        quant_test_(coeff_in_,
                    n_coeffs_,
                    skip_block,
                    zbin,
                    round,
                    quant,
                    quant_shift,
                    qcoeff_test_,
                    dqcoeff_test_,
                    dequant,
                    &eob_test_,
                    sc->scan,
                    sc->iscan);

        for (int j = 0; j < n_coeffs_; ++j) {
            ASSERT_EQ(qcoeff_ref_[j], qcoeff_test_[j])
                << "Q mismatch at position: " << j << ", Q: " << q
                << " coeff: " << coeff_in_[j];
        }

        for (int j = 0; j < n_coeffs_; ++j) {
            ASSERT_EQ(dqcoeff_ref_[j], dqcoeff_test_[j])
                << "Dq mismatch at position: " << j << ", Q: " << q
                << " coeff: " << coeff_in_[j];
        }

        ASSERT_EQ(eob_ref_, eob_test_) << "eobs mismatch, Q: " << q;
    }

    void fill_coeff_const(int i_begin, int i_end, TranLow c) {
        for (int i = i_begin; i < i_end; ++i) {
            coeff_in_[i] = c;
        }
    }

    void fill_coeff_random(int i_begin, int i_end) {
        for (int i = i_begin; i < i_end; ++i) {
            coeff_in_[i] = rnd_->random();
        }
    }

    SVTRandom *rnd_;          /**< random int for 8bit and 10bit coeffs */
    Quants qtab_quants_;      /**< quant table */
    Dequants qtab_deq_;       /**< dequant table */
    TranLow coeff_min_;       /**< min input coeff value */
    TranLow coeff_max_;       /**< max input coeff value */
    QuantizeFunc quant_ref_;  /**< reference quantize function */
    QuantizeFunc quant_test_; /**< test target quantize function */
    const TxSize tx_size_;    /**< input param tx_size */
    const AomBitDepth bd_;    /**< input param 8bit or 10bit */
    int n_coeffs_;            /**< coeff number */
    uint16_t eob_ref_;        /**< output ref eob */
    uint16_t eob_test_;       /**< output test eob */

    TranLow *coeff_in_;
    TranLow *qcoeff_ref_;
    TranLow *dqcoeff_ref_;
    TranLow *qcoeff_test_;
    TranLow *dqcoeff_test_;
};

/**
 * @brief AVX2/QuantizeBTest.input_zero_all
 *
 * test output data consistency of quantize C and avx2 functions with
 * input coef: all 0
 * q_index: 0
 */
TEST_P(QuantizeBTest, input_zero_all) {
    fill_coeff_const(0, n_coeffs_, 0);
    run_quantize(0);
}

/**
 * @brief AVX2/QuantizeBTest.input_dcac_minmax_q_n
 *
 * test output data consistency of quantize C and avx2 functions with
 * input coef: combine dc min/max with 1st ac min/max, other ac all 0
 * q_index: 0 ~ QINDEX_RANGE, step 25
 */
TEST_P(QuantizeBTest, input_dcac_minmax_q_n) {
    fill_coeff_const(2, n_coeffs_, 0);
    for (int q = 0; q < QINDEX_RANGE; q += 25) {
        fill_coeff_const(0, 1, coeff_min_);
        fill_coeff_const(1, 2, coeff_min_);
        run_quantize(q);
        fill_coeff_const(0, 1, coeff_min_);
        fill_coeff_const(1, 2, coeff_max_);
        run_quantize(q);
        fill_coeff_const(0, 1, coeff_max_);
        fill_coeff_const(1, 2, coeff_min_);
        run_quantize(q);
        fill_coeff_const(0, 1, coeff_max_);
        fill_coeff_const(1, 2, coeff_max_);
        run_quantize(q);
    }
}

/**
 * @brief AVX2/QuantizeBTest.input_random_dc_only
 *
 * test output data consistency of quantize C and avx2 functions with
 * input coef: dc random and ac all 0
 * q_index: 0
 */
TEST_P(QuantizeBTest, input_random_dc_only) {
    fill_coeff_const(1, n_coeffs_, 0);
    fill_coeff_random(0, 1);
    run_quantize(0);
}

/**
 * @brief AVX2/QuantizeBTest.input_random_all_q_all
 *
 * loop test output data consistency of quantize C and avx2 functions with
 * input coef: dc random and ac all random
 * q_index: from 0 to QINDEX_RANGE-1
 */
TEST_P(QuantizeBTest, input_random_all_q_all) {
    const int num_tests = 10;
    for (int q = 0; q < QINDEX_RANGE; ++q) {
        for (int i = 0; i < num_tests; ++i) {
            fill_coeff_random(0, n_coeffs_);
            run_quantize(q);
        }
    }
}

#ifndef FULL_UNIT_TEST
INSTANTIATE_TEST_CASE_P(
    Quant, QuantizeBTest,
    ::testing::Combine(::testing::Values(static_cast<int>(TX_16X16),
                                         static_cast<int>(TX_32X32),
                                         static_cast<int>(TX_64X64)),
                       ::testing::Values(static_cast<int>(AOM_BITS_8),
                                         static_cast<int>(AOM_BITS_10))));
#else
INSTANTIATE_TEST_CASE_P(
    Quant, QuantizeBTest,
    ::testing::Combine(::testing::Range(static_cast<int>(TX_4X4),
                                        static_cast<int>(TX_SIZES_ALL), 1),
                       ::testing::Values(static_cast<int>(AOM_BITS_8),
                                         static_cast<int>(AOM_BITS_10))));
#endif  // FULL_UNIT_TEST

}  // namespace QuantizeAsmTest
