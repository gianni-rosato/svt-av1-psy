/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file intrapred_edge_filter_test.cc
 *
 * @brief Unit test for chroma from luma prediction:
 * - eb_cfl_predict_hbd_avx2
 * - eb_cfl_predict_lbd_avx2
 *
 * @author Cidana-Wenyao
 *
 ******************************************************************************/

#include "gtest/gtest.h"
#include "aom_dsp_rtcd.h"
#include "EbDefinitions.h"
#include "random.h"
namespace {
using svt_av1_test_tool::SVTRandom;

using CFL_PRED_HBD = void (*)(const int16_t *pred_buf_q3, uint16_t *pred,
                              int32_t pred_stride, uint16_t *dst,
                              int32_t dst_stride, int32_t alpha_q3,
                              int32_t bit_depth, int32_t width, int32_t height);
using CFL_PRED_LBD = void (*)(const int16_t *pred_buf_q3, uint8_t *pred,
                              int32_t pred_stride, uint8_t *dst,
                              int32_t dst_stride, int32_t alpha_q3,
                              int32_t bit_depth, int32_t width, int32_t height);
/**
 * @brief Unit test for chroma from luma prediction:
 * - eb_cfl_predict_hbd_avx2
 * - eb_cfl_predict_lbd_avx2
 *
 * Test strategy:
 * Verify this assembly code by comparing with reference c implementation.
 * Feed the same data and check test output and reference output.
 * Define a templete class to handle the common process, and
 * declare sub class to handle different bitdepth and function types.
 *
 * Expect result:
 * Output from assemble functions should be the same with output from c.
 *
 * Test coverage:
 * Test cases:
 * pred buffer and dst buffer: Fill with random values
 * TxSize: all the TxSize.
 * alpha_q3: [-16, 16]
 * BitDepth: 8bit and 10bit
 */
template <typename Sample, typename FuncType>
class CflPredTest {
  public:
    CflPredTest() {
        bd_ = 8;
        ref_func_ = nullptr;
        tst_func_ = nullptr;
        common_init();
    }

    virtual ~CflPredTest() {
    }

    void RunAllTest() {
        // for pred_buf, after sampling and subtracted from average
        SVTRandom pred_rnd(bd_ + 3 + 1, true);
        SVTRandom dst_rnd(8, false);
        for (int tx = TX_4X4; tx < TX_SIZES_ALL; ++tx) {
            const int c_w = tx_size_wide[tx] >> 1;
            const int c_h = tx_size_high[tx] >> 1;
            const int c_stride = CFL_BUF_LINE;
            memset(pred_buf_q3, 0, sizeof(pred_buf_q3));
            memset(dst_buf_ref_data_, 0, sizeof(dst_buf_ref_data_));
            memset(dst_buf_tst_data_, 0, sizeof(dst_buf_tst_data_));

            for (int alpha_q3 = -16; alpha_q3 <= 16; ++alpha_q3) {
                // prepare data
                // Implicit assumption: The dst_buf is supposed to be populated
                // by dc prediction before cfl prediction.
                const Sample rnd_dc = dst_rnd.random();
                for (int y = 0; y < c_h; ++y) {
                    for (int x = 0; x < c_w; ++x) {
                        pred_buf_q3[y * c_stride + x] =
                            (Sample)pred_rnd.random();
                        dst_buf_ref_[y * c_stride + x] =
                            dst_buf_tst_[y * c_stride + x] = rnd_dc;
                    }
                }

                ref_func_(pred_buf_q3,
                          dst_buf_ref_,
                          CFL_BUF_LINE,
                          dst_buf_ref_,
                          CFL_BUF_LINE,
                          alpha_q3,
                          bd_,
                          c_w,
                          c_h);
                tst_func_(pred_buf_q3,
                          dst_buf_tst_,
                          c_stride,
                          dst_buf_tst_,
                          c_stride,
                          alpha_q3,
                          bd_,
                          c_w,
                          c_h);

                for (int y = 0; y < c_h; ++y) {
                    for (int x = 0; x < c_w; ++x) {
                        ASSERT_EQ(dst_buf_ref_[y * c_stride + x],
                                  dst_buf_tst_[y * c_stride + x])
                            << "tx_size: " << tx << " alpha_q3 " << alpha_q3
                            << " expect " << dst_buf_ref_[y * c_stride + x]
                            << " got " << dst_buf_tst_[y * c_stride + x]
                            << " at [ " << x << " x " << y << " ]";
                    }
                }
            }
        }
    }

  protected:
    void common_init() {
        dst_buf_ref_ = reinterpret_cast<Sample *>(
            ((intptr_t)(dst_buf_ref_data_) + alignment - 1) & ~(alignment - 1));
        dst_buf_tst_ = reinterpret_cast<Sample *>(
            ((intptr_t)(dst_buf_tst_data_) + alignment - 1) & ~(alignment - 1));
    }

    static const int alignment = 32;
    int16_t pred_buf_q3[CFL_BUF_SQUARE];
    Sample dst_buf_ref_data_[CFL_BUF_SQUARE + alignment - 1];
    Sample dst_buf_tst_data_[CFL_BUF_SQUARE + alignment - 1];
    Sample *dst_buf_ref_;
    Sample *dst_buf_tst_;
    FuncType ref_func_;
    FuncType tst_func_;
    int bd_;
};

class LbdCflPredTest : public CflPredTest<uint8_t, CFL_PRED_LBD> {
  public:
    LbdCflPredTest() {
        bd_ = 8;
        ref_func_ = eb_cfl_predict_lbd_c;
        tst_func_ = eb_cfl_predict_lbd_avx2;
        common_init();
    }
};

class HbdCflPredTest : public CflPredTest<uint16_t, CFL_PRED_HBD> {
  public:
    HbdCflPredTest() {
        bd_ = 10;
        ref_func_ = eb_cfl_predict_hbd_c;
        tst_func_ = eb_cfl_predict_hbd_avx2;
        common_init();
    }
};

#define TEST_CLASS(tc_name, type_name)     \
    TEST(tc_name, match_test) {            \
        type_name *test = new type_name(); \
        test->RunAllTest();                \
        delete test;                       \
    }

TEST_CLASS(LbdCflPredMatchTest, LbdCflPredTest)
TEST_CLASS(HbdCflPredMatchTest, HbdCflPredTest)
}  // namespace
