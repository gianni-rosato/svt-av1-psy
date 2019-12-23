/*
 * Copyright (c) 2017, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

#include "gtest/gtest.h"

#include "aom_dsp_rtcd.h"
#include "random.h"
#include "util.h"

// workaround to eliminate the compiling warning on linux
// The macro will conflict with definition in gtest.h
#ifdef __USE_GNU
#undef __USE_GNU  // defined in EbThreads.h
#endif
#ifdef _GNU_SOURCE
#undef _GNU_SOURCE  // defined in EbThreads.h
#endif
#include "EbDefinitions.h"
#include "EbPictureControlSet.h"
#include "EbTransforms.h"
#include "EbUnitTestUtility.h"

namespace {
using svt_av1_test_tool::SVTRandom;

extern "C" void eb_av1_build_quantizer(
    AomBitDepth bit_depth, int32_t y_dc_delta_q, int32_t u_dc_delta_q,
    int32_t u_ac_delta_q, int32_t v_dc_delta_q, int32_t v_ac_delta_q,
    Quants *const quants, Dequants *const deq);

#define QUAN_PARAM_LIST                                                      \
    const TranLow *coeff_ptr, intptr_t n_coeffs, const int16_t *zbin_ptr,    \
        const int16_t *round_ptr, const int16_t *quant_ptr,                  \
        const int16_t *quant_shift_ptr, TranLow *qcoeff_ptr,                 \
        TranLow *dqcoeff_ptr, const int16_t *dequant_ptr, uint16_t *eob_ptr, \
        const int16_t *scan, const int16_t *iscan
#define QUAN_HBD_PARAM int16_t log_scale

typedef void (*QuantizeFunc)(QUAN_PARAM_LIST);
typedef void (*QuantizeHbdFunc)(QUAN_PARAM_LIST, QUAN_HBD_PARAM);

enum QuantType { TYPE_B, TYPE_DC, TYPE_FP };

using std::tuple;
typedef tuple<QuantizeFunc, QuantizeFunc, TxSize, QuantType, AomBitDepth>
    QuantizeParam;
typedef tuple<QuantizeHbdFunc, QuantizeHbdFunc, TxSize, QuantType, AomBitDepth>
    QuantizeHbdParam;

typedef struct {
    Quants quant;
    Dequants dequant;
} QuanTable;

const int kTestNum = 1000;

template <typename ParamType, typename FuncType>
class QuantizeTest : public ::testing::TestWithParam<ParamType> {
  protected:
    QuantizeTest()
        : rnd_(32, true),
          quant_ref_(nullptr),
          quant_(nullptr),
          tx_size_(TX_4X4),
          type_(TYPE_FP),
          bd_(AOM_BITS_8) {
    }

    virtual ~QuantizeTest() {
    }

    virtual void SetUp() {
        qtab_ =
            reinterpret_cast<QuanTable *>(eb_aom_memalign(32, sizeof(*qtab_)));
        const int n_coeffs = coeff_num();
        coeff_ = reinterpret_cast<TranLow *>(
            eb_aom_memalign(32, 6 * n_coeffs * sizeof(TranLow)));
        InitQuantizer();
    }

    virtual void TearDown() {
        eb_aom_free(qtab_);
        qtab_ = NULL;
        eb_aom_free(coeff_);
        coeff_ = NULL;
    }

    void InitQuantizer() {
        eb_av1_build_quantizer(
            bd_, 0, 0, 0, 0, 0, &qtab_->quant, &qtab_->dequant);
    }

    virtual void QuantizeRun(bool is_loop, int q = 0, int test_num = 1) = 0;

    void CompareResults(const TranLow *buf_ref, const TranLow *buf, int size,
                        const char *text, int q, int number) {
        int i;
        for (i = 0; i < size; ++i) {
            ASSERT_EQ(buf_ref[i], buf[i])
                << text << " mismatch on test: " << number
                << " at position: " << i << " Q: " << q;
        }
    }

    int coeff_num() const {
        return av1_get_max_eob(tx_size_);
    }

    void FillCoeff(TranLow c) {
        const int n_coeffs = coeff_num();
        for (int i = 0; i < n_coeffs; ++i) {
            coeff_[i] = c;
        }
    }

    void FillCoeffRandom() {
        const int n_coeffs = coeff_num();
        FillCoeffZero();
        int num = rnd_.random() % n_coeffs;
        for (int i = 0; i < num; ++i) {
            coeff_[i] = GetRandomCoeff();
        }
    }

    void FillCoeffRandomRows(int num) {
        FillCoeffZero();
        for (int i = 0; i < num; ++i) {
            coeff_[i] = GetRandomCoeff();
        }
    }

    void FillCoeffZero() {
        FillCoeff(0);
    }

    void FillCoeffConstant() {
        TranLow c = GetRandomCoeff();
        FillCoeff(c);
    }

    void FillDcOnly() {
        FillCoeffZero();
        coeff_[0] = GetRandomCoeff();
    }

    void FillDcLargeNegative() {
        FillCoeffZero();
        // Generate a qcoeff which contains 512/-512 (0x0100/0xFE00) to catch
        // issues like BUG=883 where the constant being compared was incorrectly
        // initialized.
        coeff_[0] = -8191;
    }

    TranLow GetRandomCoeff() {
        TranLow coeff;
        if (bd_ == AOM_BITS_8) {
            coeff = clamp(
                static_cast<int16_t>(rnd_.random()), INT16_MIN + 1, INT16_MAX);
        } else {
            TranLow min = -(1 << (7 + bd_));
            TranLow max = -min - 1;
            coeff = clamp(static_cast<TranLow>(rnd_.random()), min, max);
        }
        return coeff;
    }

    SVTRandom rnd_;
    QuanTable *qtab_;
    TranLow *coeff_;
    FuncType quant_ref_;
    FuncType quant_;
    TxSize tx_size_;
    QuantType type_;
    AomBitDepth bd_;
};

class QuantizeLbdTest : public QuantizeTest<QuantizeParam, QuantizeFunc> {
  protected:
    QuantizeLbdTest() {
        quant_ref_ = TEST_GET_PARAM(0);
        quant_ = TEST_GET_PARAM(1);
        tx_size_ = TEST_GET_PARAM(2);
        type_ = TEST_GET_PARAM(3);
        bd_ = TEST_GET_PARAM(4);
    }

    void QuantizeRun(bool is_loop, int q = 0, int test_num = 1) override {
        TranLow *coeff_ptr = coeff_;
        const intptr_t n_coeffs = coeff_num();

        TranLow *qcoeff_ref = coeff_ptr + n_coeffs;
        TranLow *dqcoeff_ref = qcoeff_ref + n_coeffs;

        TranLow *qcoeff = dqcoeff_ref + n_coeffs;
        TranLow *dqcoeff = qcoeff + n_coeffs;
        uint16_t *eob = (uint16_t *)(dqcoeff + n_coeffs);

        // Testing uses 2-D DCT scan order table
        const ScanOrder *const sc = &av1_scan_orders[tx_size_][DCT_DCT];

        // Testing uses luminance quantization table
        const int16_t *zbin = qtab_->quant.y_zbin[q];

        const int16_t *round = 0;
        const int16_t *quant = 0;
        if (type_ == TYPE_B) {
            round = qtab_->quant.y_round[q];
            quant = qtab_->quant.y_quant[q];
        } else if (type_ == TYPE_FP) {
            round = qtab_->quant.y_round_fp[q];
            quant = qtab_->quant.y_quant_fp[q];
        }

        const int16_t *quant_shift = qtab_->quant.y_quant_shift[q];
        const int16_t *dequant = qtab_->dequant.y_dequant_qtx[q];

        for (int i = 0; i < test_num; ++i) {
            if (is_loop)
                FillCoeffRandom();

            memset(qcoeff_ref, 0, 5 * n_coeffs * sizeof(*qcoeff_ref));

            quant_ref_(coeff_ptr,
                       n_coeffs,
                       zbin,
                       round,
                       quant,
                       quant_shift,
                       qcoeff_ref,
                       dqcoeff_ref,
                       dequant,
                       &eob[0],
                       sc->scan,
                       sc->iscan);

            quant_(coeff_ptr,
                   n_coeffs,
                   zbin,
                   round,
                   quant,
                   quant_shift,
                   qcoeff,
                   dqcoeff,
                   dequant,
                   &eob[1],
                   sc->scan,
                   sc->iscan);

            for (int j = 0; j < n_coeffs; ++j) {
                ASSERT_EQ(qcoeff_ref[j], qcoeff[j])
                    << "Q mismatch on test: " << i << " at position: " << j
                    << " Q: " << q << " coeff: " << coeff_ptr[j];
            }

            for (int j = 0; j < n_coeffs; ++j) {
                ASSERT_EQ(dqcoeff_ref[j], dqcoeff[j])
                    << "Dq mismatch on test: " << i << " at position: " << j
                    << " Q: " << q << " coeff: " << coeff_ptr[j];
            }

            ASSERT_EQ(eob[0], eob[1])
                << "eobs mismatch on test: " << i << " Q: " << q;
        }
    }
};

TEST_P(QuantizeLbdTest, ZeroInput) {
    FillCoeffZero();
    QuantizeRun(false);
}

TEST_P(QuantizeLbdTest, LargeNegativeInput) {
    FillDcLargeNegative();
    QuantizeRun(false, 0, 1);
}

TEST_P(QuantizeLbdTest, DcOnlyInput) {
    FillDcOnly();
    QuantizeRun(false, 0, 1);
}

TEST_P(QuantizeLbdTest, RandomInput) {
    QuantizeRun(true, 0, kTestNum);
}

TEST_P(QuantizeLbdTest, MultipleQ) {
    for (int q = 0; q < QINDEX_RANGE; ++q) {
        QuantizeRun(true, q, kTestNum);
    }
}

// Force the coeff to be half the value of the dequant.  This exposes a
// mismatch found in av1_quantize_fp_sse2().
TEST_P(QuantizeLbdTest, CoeffHalfDequant) {
    FillCoeff(16);
    QuantizeRun(false, 25, 1);
}

TEST_P(QuantizeLbdTest, DISABLED_Speed) {
    TranLow *coeff_ptr = coeff_;
    const intptr_t n_coeffs = coeff_num();

    TranLow *qcoeff_ref = coeff_ptr + n_coeffs;
    TranLow *dqcoeff_ref = qcoeff_ref + n_coeffs;

    TranLow *qcoeff = dqcoeff_ref + n_coeffs;
    TranLow *dqcoeff = qcoeff + n_coeffs;
    uint16_t *eob = (uint16_t *)(dqcoeff + n_coeffs);

    // Testing uses 2-D DCT scan order table
    const ScanOrder *const sc = &av1_scan_orders[tx_size_][DCT_DCT];

    // Testing uses luminance quantization table
    const int q = 22;
    const int16_t *zbin = qtab_->quant.y_zbin[q];
    const int16_t *round_fp = qtab_->quant.y_round_fp[q];
    const int16_t *quant_fp = qtab_->quant.y_quant_fp[q];
    const int16_t *quant_shift = qtab_->quant.y_quant_shift[q];
    const int16_t *dequant = qtab_->dequant.y_dequant_qtx[q];
    const int kNumTests = 5000000;
    double time_c, time_o;
    uint64_t start_time_seconds, start_time_useconds;
    uint64_t middle_time_seconds, middle_time_useconds;
    uint64_t finish_time_seconds, finish_time_useconds;
    int rows = tx_size_high[tx_size_];
    int cols = tx_size_wide[tx_size_];
    rows = AOMMIN(32, rows);
    cols = AOMMIN(32, cols);
    for (int cnt = 0; cnt <= rows; cnt++) {
        FillCoeffRandomRows(cnt * cols);

        eb_start_time(&start_time_seconds, &start_time_useconds);
        for (int n = 0; n < kNumTests; ++n) {
            quant_ref_(coeff_ptr,
                       n_coeffs,
                       zbin,
                       round_fp,
                       quant_fp,
                       quant_shift,
                       qcoeff,
                       dqcoeff,
                       dequant,
                       eob,
                       sc->scan,
                       sc->iscan);
        }

        eb_start_time(&middle_time_seconds, &middle_time_useconds);

        for (int n = 0; n < kNumTests; ++n) {
            quant_(coeff_ptr,
                   n_coeffs,
                   zbin,
                   round_fp,
                   quant_fp,
                   quant_shift,
                   qcoeff,
                   dqcoeff,
                   dequant,
                   eob,
                   sc->scan,
                   sc->iscan);
        }
        eb_start_time(&finish_time_seconds, &finish_time_useconds);
        eb_compute_overall_elapsed_time_ms(start_time_seconds,
                                      start_time_useconds,
                                      middle_time_seconds,
                                      middle_time_useconds,
                                      &time_c);
        eb_compute_overall_elapsed_time_ms(middle_time_seconds,
                                      middle_time_useconds,
                                      finish_time_seconds,
                                      finish_time_useconds,
                                      &time_o);

        printf("c_time = %f \t simd_time = %f \t Gain = %f \n",
               time_c,
               time_o,
               (time_c / time_o));
    }
}

class QuantizeHbdTest : public QuantizeTest<QuantizeHbdParam, QuantizeHbdFunc> {
  protected:
    QuantizeHbdTest() {
        quant_ref_ = TEST_GET_PARAM(0);
        quant_ = TEST_GET_PARAM(1);
        tx_size_ = TEST_GET_PARAM(2);
        type_ = TEST_GET_PARAM(3);
        bd_ = TEST_GET_PARAM(4);
    }

    void QuantizeRun(bool is_loop, int q = 0, int test_num = 1) override {
        TranLow *coeff_ptr = coeff_;
        const intptr_t n_coeffs = coeff_num();

        TranLow *qcoeff_ref = coeff_ptr + n_coeffs;
        TranLow *dqcoeff_ref = qcoeff_ref + n_coeffs;

        TranLow *qcoeff = dqcoeff_ref + n_coeffs;
        TranLow *dqcoeff = qcoeff + n_coeffs;
        uint16_t *eob = (uint16_t *)(dqcoeff + n_coeffs);

        // Testing uses 2-D DCT scan order table
        const ScanOrder *const sc = &av1_scan_orders[tx_size_][DCT_DCT];

        // Testing uses luminance quantization table
        const int16_t *zbin = qtab_->quant.y_zbin[q];

        const int16_t *round = 0;
        const int16_t *quant = 0;
        if (type_ == TYPE_B) {
            round = qtab_->quant.y_round[q];
            quant = qtab_->quant.y_quant[q];
        } else if (type_ == TYPE_FP) {
            round = qtab_->quant.y_round_fp[q];
            quant = qtab_->quant.y_quant_fp[q];
        }

        const int16_t *quant_shift = qtab_->quant.y_quant_shift[q];
        const int16_t *dequant = qtab_->dequant.y_dequant_qtx[q];

        for (int i = 0; i < test_num; ++i) {
            if (is_loop)
                FillCoeffRandom();

            memset(qcoeff_ref, 0, 5 * n_coeffs * sizeof(*qcoeff_ref));

            quant_ref_(coeff_ptr,
                       n_coeffs,
                       zbin,
                       round,
                       quant,
                       quant_shift,
                       qcoeff_ref,
                       dqcoeff_ref,
                       dequant,
                       &eob[0],
                       sc->scan,
                       sc->iscan,
                       av1_get_tx_scale(tx_size_));

            quant_(coeff_ptr,
                   n_coeffs,
                   zbin,
                   round,
                   quant,
                   quant_shift,
                   qcoeff,
                   dqcoeff,
                   dequant,
                   &eob[1],
                   sc->scan,
                   sc->iscan,
                   av1_get_tx_scale(tx_size_));

            for (int j = 0; j < n_coeffs; ++j) {
                ASSERT_EQ(qcoeff_ref[j], qcoeff[j])
                    << "Q mismatch on test: " << i << " at position: " << j
                    << " Q: " << q << " coeff: " << coeff_ptr[j];
            }

            for (int j = 0; j < n_coeffs; ++j) {
                ASSERT_EQ(dqcoeff_ref[j], dqcoeff[j])
                    << "Dq mismatch on test: " << i << " at position: " << j
                    << " Q: " << q << " coeff: " << coeff_ptr[j];
            }

            ASSERT_EQ(eob[0], eob[1])
                << "eobs mismatch on test: " << i << " Q: " << q;
        }
    }
};

TEST_P(QuantizeHbdTest, ZeroInput) {
    FillCoeffZero();
    QuantizeRun(false);
}

TEST_P(QuantizeHbdTest, LargeNegativeInput) {
    FillDcLargeNegative();
    QuantizeRun(false, 0, 1);
}

TEST_P(QuantizeHbdTest, DcOnlyInput) {
    FillDcOnly();
    QuantizeRun(false, 0, 1);
}

TEST_P(QuantizeHbdTest, RandomInput) {
    QuantizeRun(true, 0, kTestNum);
}

TEST_P(QuantizeHbdTest, MultipleQ) {
    for (int q = 0; q < QINDEX_RANGE; ++q) {
        QuantizeRun(true, q, kTestNum);
    }
}

// Force the coeff to be half the value of the dequant.  This exposes a
// mismatch found in av1_quantize_fp_sse2().
TEST_P(QuantizeHbdTest, CoeffHalfDequant) {
    FillCoeff(16);
    QuantizeRun(false, 25, 1);
}

TEST_P(QuantizeHbdTest, DISABLED_Speed) {
    TranLow *coeff_ptr = coeff_;
    const intptr_t n_coeffs = coeff_num();

    TranLow *qcoeff_ref = coeff_ptr + n_coeffs;
    TranLow *dqcoeff_ref = qcoeff_ref + n_coeffs;

    TranLow *qcoeff = dqcoeff_ref + n_coeffs;
    TranLow *dqcoeff = qcoeff + n_coeffs;
    uint16_t *eob = (uint16_t *)(dqcoeff + n_coeffs);

    // Testing uses 2-D DCT scan order table
    const ScanOrder *const sc = &av1_scan_orders[tx_size_][DCT_DCT];

    // Testing uses luminance quantization table
    const int q = 22;
    const int16_t *zbin = qtab_->quant.y_zbin[q];
    const int16_t *round_fp = qtab_->quant.y_round_fp[q];
    const int16_t *quant_fp = qtab_->quant.y_quant_fp[q];
    const int16_t *quant_shift = qtab_->quant.y_quant_shift[q];
    const int16_t *dequant = qtab_->dequant.y_dequant_qtx[q];
    const int kNumTests = 5000000;
    double time_c, time_o;
    uint64_t start_time_seconds, start_time_useconds;
    uint64_t middle_time_seconds, middle_time_useconds;
    uint64_t finish_time_seconds, finish_time_useconds;
    int rows = tx_size_high[tx_size_];
    int cols = tx_size_wide[tx_size_];
    rows = AOMMIN(32, rows);
    cols = AOMMIN(32, cols);
    for (int cnt = 0; cnt <= rows; cnt++) {
        FillCoeffRandomRows(cnt * cols);

        eb_start_time(&start_time_seconds, &start_time_useconds);
        for (int n = 0; n < kNumTests; ++n) {
            quant_ref_(coeff_ptr,
                       n_coeffs,
                       zbin,
                       round_fp,
                       quant_fp,
                       quant_shift,
                       qcoeff,
                       dqcoeff,
                       dequant,
                       eob,
                       sc->scan,
                       sc->iscan,
                       av1_get_tx_scale(tx_size_));
        }

        eb_start_time(&middle_time_seconds, &middle_time_useconds);

        for (int n = 0; n < kNumTests; ++n) {
            quant_(coeff_ptr,
                   n_coeffs,
                   zbin,
                   round_fp,
                   quant_fp,
                   quant_shift,
                   qcoeff,
                   dqcoeff,
                   dequant,
                   eob,
                   sc->scan,
                   sc->iscan,
                   av1_get_tx_scale(tx_size_));
        }
        eb_start_time(&finish_time_seconds, &finish_time_useconds);
        eb_compute_overall_elapsed_time_ms(start_time_seconds,
                                      start_time_useconds,
                                      middle_time_seconds,
                                      middle_time_useconds,
                                      &time_c);
        eb_compute_overall_elapsed_time_ms(middle_time_seconds,
                                      middle_time_useconds,
                                      finish_time_seconds,
                                      finish_time_useconds,
                                      &time_o);

        printf("c_time = %f \t simd_time = %f \t Gain = %f \n",
               time_c,
               time_o,
               (time_c / time_o));
    }
}

using std::make_tuple;

#if HAS_AVX2
const QuantizeParam kQParamArrayAvx2[] = {
    make_tuple(&eb_av1_quantize_fp_c, &eb_av1_quantize_fp_avx2,
               static_cast<TxSize>(TX_16X16), TYPE_FP, AOM_BITS_8),
    make_tuple(&eb_av1_quantize_fp_c, &eb_av1_quantize_fp_avx2,
               static_cast<TxSize>(TX_4X16), TYPE_FP, AOM_BITS_8),
    make_tuple(&eb_av1_quantize_fp_c, &eb_av1_quantize_fp_avx2,
               static_cast<TxSize>(TX_16X4), TYPE_FP, AOM_BITS_8),
    make_tuple(&eb_av1_quantize_fp_c, &eb_av1_quantize_fp_avx2,
               static_cast<TxSize>(TX_32X8), TYPE_FP, AOM_BITS_8),
    make_tuple(&eb_av1_quantize_fp_c, &eb_av1_quantize_fp_avx2,
               static_cast<TxSize>(TX_8X32), TYPE_FP, AOM_BITS_8),
    make_tuple(&eb_av1_quantize_fp_32x32_c, &eb_av1_quantize_fp_32x32_avx2,
               static_cast<TxSize>(TX_32X32), TYPE_FP, AOM_BITS_8),
    make_tuple(&eb_av1_quantize_fp_32x32_c, &eb_av1_quantize_fp_32x32_avx2,
               static_cast<TxSize>(TX_16X64), TYPE_FP, AOM_BITS_8),
    make_tuple(&eb_av1_quantize_fp_32x32_c, &eb_av1_quantize_fp_32x32_avx2,
               static_cast<TxSize>(TX_64X16), TYPE_FP, AOM_BITS_8),
    make_tuple(&eb_av1_quantize_fp_64x64_c, &eb_av1_quantize_fp_64x64_avx2,
               static_cast<TxSize>(TX_64X64), TYPE_FP, AOM_BITS_8)};

const QuantizeHbdParam kQHbdParamArrayAvx2[] = {
    make_tuple(&eb_av1_highbd_quantize_fp_c, &eb_av1_highbd_quantize_fp_avx2,
               static_cast<TxSize>(TX_16X16), TYPE_FP, AOM_BITS_8),
    make_tuple(&eb_av1_highbd_quantize_fp_c, &eb_av1_highbd_quantize_fp_avx2,
               static_cast<TxSize>(TX_4X16), TYPE_FP, AOM_BITS_8),
    make_tuple(&eb_av1_highbd_quantize_fp_c, &eb_av1_highbd_quantize_fp_avx2,
               static_cast<TxSize>(TX_16X4), TYPE_FP, AOM_BITS_8),
    make_tuple(&eb_av1_highbd_quantize_fp_c, &eb_av1_highbd_quantize_fp_avx2,
               static_cast<TxSize>(TX_32X8), TYPE_FP, AOM_BITS_8),
    make_tuple(&eb_av1_highbd_quantize_fp_c, &eb_av1_highbd_quantize_fp_avx2,
               static_cast<TxSize>(TX_8X32), TYPE_FP, AOM_BITS_8),
    make_tuple(&eb_av1_highbd_quantize_fp_c, &eb_av1_highbd_quantize_fp_avx2,
               static_cast<TxSize>(TX_16X16), TYPE_FP, AOM_BITS_10),
    make_tuple(&eb_av1_highbd_quantize_fp_c, &eb_av1_highbd_quantize_fp_avx2,
               static_cast<TxSize>(TX_4X16), TYPE_FP, AOM_BITS_10),
    make_tuple(&eb_av1_highbd_quantize_fp_c, &eb_av1_highbd_quantize_fp_avx2,
               static_cast<TxSize>(TX_16X4), TYPE_FP, AOM_BITS_10),
    make_tuple(&eb_av1_highbd_quantize_fp_c, &eb_av1_highbd_quantize_fp_avx2,
               static_cast<TxSize>(TX_32X8), TYPE_FP, AOM_BITS_10),
    make_tuple(&eb_av1_highbd_quantize_fp_c, &eb_av1_highbd_quantize_fp_avx2,
               static_cast<TxSize>(TX_8X32), TYPE_FP, AOM_BITS_10),
    make_tuple(&eb_av1_highbd_quantize_fp_c, &eb_av1_highbd_quantize_fp_avx2,
               static_cast<TxSize>(TX_16X16), TYPE_FP, AOM_BITS_12),
    make_tuple(&eb_av1_highbd_quantize_fp_c, &eb_av1_highbd_quantize_fp_avx2,
               static_cast<TxSize>(TX_4X16), TYPE_FP, AOM_BITS_12),
    make_tuple(&eb_av1_highbd_quantize_fp_c, &eb_av1_highbd_quantize_fp_avx2,
               static_cast<TxSize>(TX_16X4), TYPE_FP, AOM_BITS_12),
    make_tuple(&eb_av1_highbd_quantize_fp_c, &eb_av1_highbd_quantize_fp_avx2,
               static_cast<TxSize>(TX_32X8), TYPE_FP, AOM_BITS_12),
    make_tuple(&eb_av1_highbd_quantize_fp_c, &eb_av1_highbd_quantize_fp_avx2,
               static_cast<TxSize>(TX_8X32), TYPE_FP, AOM_BITS_12)};

INSTANTIATE_TEST_CASE_P(AVX2, QuantizeLbdTest,
                        ::testing::ValuesIn(kQParamArrayAvx2));
INSTANTIATE_TEST_CASE_P(AVX2, QuantizeHbdTest,
                        ::testing::ValuesIn(kQHbdParamArrayAvx2));
#endif  // HAS_AVX2
}  // namespace
