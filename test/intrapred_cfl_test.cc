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
 * - cfl_luma_subsampling_420_lbd_avx2
 * - cfl_luma_subsampling_420_hbd_avx2
 *
 * @author Cidana-Wenyao
 *
 ******************************************************************************/

#include "gtest/gtest.h"
#include "aom_dsp_rtcd.h"
#include "EbDefinitions.h"
#include "random.h"
#include "util.h"
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

typedef void (*AomUpsampledPredFunc)(MacroBlockD *,
                                     const struct AV1Common *const, int, int,
                                     const MV *const, uint8_t *, int, int, int,
                                     int, const uint8_t *, int, int);
typedef ::testing::tuple<BlockSize, AomUpsampledPredFunc, int, int, int> AomUpsampledPredParam;

class AomUpsampledPredTest
    : public ::testing::TestWithParam<AomUpsampledPredParam> {
  public:
    AomUpsampledPredTest() : rnd_(0, 255){};
    virtual ~AomUpsampledPredTest() {
    }

    void TearDown() override {
        aom_clear_system_state();
    }

    void run_test() {
        const int block_size = TEST_GET_PARAM(0);
        AomUpsampledPredFunc test_impl = TEST_GET_PARAM(1);
        int subpel_search = TEST_GET_PARAM(2);
        int subpel_x_q3 = TEST_GET_PARAM(3);
        int subpel_y_q3 = TEST_GET_PARAM(4);
        const int width = block_size_wide[block_size];
        const int height = block_size_high[block_size];
        DECLARE_ALIGNED(16, uint8_t, ref_[MAX_SB_SQUARE * 2]);
        DECLARE_ALIGNED(16, uint8_t, comp_pred_ref_[MAX_SB_SQUARE * 2]);
        DECLARE_ALIGNED(16, uint8_t, comp_pred_tst_[MAX_SB_SQUARE * 2]);

        memset(comp_pred_ref_, 1, sizeof(comp_pred_ref_));
        memset(comp_pred_tst_, 1, sizeof(comp_pred_tst_));

        //Function aom_upsampled_pred_sse2 call inside function pointer which have to be set properly
        // by setup_common_rtcd_internal(), we want to test intrinsic version of it, so AVX2 flag is necessary
        setup_common_rtcd_internal(CPU_FLAGS_AVX2);

        const int run_times = 100;
        for (int i = 0; i < run_times; ++i) {
            memset(ref_, 1, sizeof(ref_));
            for (int j = 0; j < width * height+ 3 * width; j++) {
                ref_[j] = rnd_.random();
            }

            aom_upsampled_pred_c(NULL,
                                 NULL,
                                 0,
                                 0,
                                 NULL,
                                 comp_pred_ref_,
                                 width,
                                 height,
                                 subpel_x_q3,
                                 subpel_y_q3,
                                 ref_ + 3 * width,
                                 width,
                                 subpel_search);
            test_impl(NULL,
                      NULL,
                      0,
                      0,
                      NULL,
                      comp_pred_tst_,
                      width,
                      height,
                      subpel_x_q3,
                      subpel_y_q3,
                      ref_ + 3 * width,
                      width,
                      subpel_search);

            ASSERT_EQ(
                0,
                memcmp(comp_pred_ref_, comp_pred_tst_, sizeof(comp_pred_ref_)));
        }
    }

  private:
    SVTRandom rnd_;
};

TEST_P(AomUpsampledPredTest, MatchTest) {
    run_test();
}

INSTANTIATE_TEST_CASE_P(
    UPSAMPLED_PRED_TEST, AomUpsampledPredTest,
    ::testing::Combine(::testing::Range(BLOCK_4X4, BlockSizeS_ALL),
                       ::testing::Values(aom_upsampled_pred_sse2),
                       ::testing::Values(USE_2_TAPS, USE_4_TAPS, USE_8_TAPS),
                       ::testing::Values(0, 1, 2),
                       ::testing::Values(0, 1, 2)));


typedef void (*CflLumaSubsamplingLbdFunc)(const uint8_t *, int32_t, int16_t *, int32_t,
                                   int32_t);
typedef ::testing::tuple<BlockSize, CflLumaSubsamplingLbdFunc>
    CflLumaSubsamplingLbdParam;

class CflLumaSubsamplingLbdTest
    : public ::testing::TestWithParam<CflLumaSubsamplingLbdParam> {
  public:
    CflLumaSubsamplingLbdTest() : rnd_(0, 255){};
    virtual ~CflLumaSubsamplingLbdTest() {
    }

    void TearDown() override {
        aom_clear_system_state();
    }

    void run_test() {
        const int block_size = TEST_GET_PARAM(0);
        CflLumaSubsamplingLbdFunc test_impl = TEST_GET_PARAM(1);
        const int width = block_size_wide[block_size];
        // Output width is defined by CFL_BUF_LINE(32),
        // it lead to assumption that input width cannot be larger than 64,
        // otherwise computation will overwrite line "n" by line "n+1"
        if (width > 64)
            return;
        const int height = block_size_high[block_size];
        DECLARE_ALIGNED(16, uint8_t, input[MAX_SB_SQUARE]);
        DECLARE_ALIGNED(16, int16_t, output_q3_ref_[MAX_SB_SQUARE]);
        DECLARE_ALIGNED(16, int16_t, output_q3_tst_[MAX_SB_SQUARE]);

        memset(output_q3_ref_, 1, sizeof(output_q3_ref_));
        memset(output_q3_tst_, 1, sizeof(output_q3_tst_));

        const int run_times = 100;
        for (int i = 0; i < run_times; ++i) {

            memset(input, 1, sizeof(input));
            for (int j = 0; j < MAX_SB_SQUARE; j++) {
                input[j] = rnd_.random();
            }

            cfl_luma_subsampling_420_lbd_c(
                input, width, output_q3_ref_, width, height);

            test_impl(input, width, output_q3_tst_, width, height);

            ASSERT_EQ(
                0,
                memcmp(output_q3_ref_, output_q3_tst_, sizeof(output_q3_ref_)));

        }
    }

  private:
    SVTRandom rnd_;
};

TEST_P(CflLumaSubsamplingLbdTest, MatchTest) {
    run_test();
}

INSTANTIATE_TEST_CASE_P(
    CFL_LUMA_SUBSAMPLING_LBD, CflLumaSubsamplingLbdTest,
    ::testing::Combine(::testing::Range(BLOCK_4X4, BlockSizeS_ALL),
                       ::testing::Values(cfl_luma_subsampling_420_lbd_avx2)));


typedef void (*CflLumaSubsamplingHbdFunc)(const uint16_t *, int32_t, int16_t *,
                                          int32_t, int32_t);
typedef ::testing::tuple<BlockSize, CflLumaSubsamplingHbdFunc>
    CflLumaSubsamplingHbdParam;

class CflLumaSubsamplingHbdTest
    : public ::testing::TestWithParam<CflLumaSubsamplingHbdParam> {
  public:
    CflLumaSubsamplingHbdTest() : rnd_(0, 1023){};
    virtual ~CflLumaSubsamplingHbdTest() {
    }

    void TearDown() override {
        aom_clear_system_state();
    }

    void run_test() {
        const int block_size = TEST_GET_PARAM(0);
        CflLumaSubsamplingHbdFunc test_impl = TEST_GET_PARAM(1);
        const int width = block_size_wide[block_size];
        // Output width is defined by CFL_BUF_LINE(32),
        // it lead to assumption that input width cannot be larger than 64,
        // otherwise computation will overwrite line "n" by line "n+1"
        if (width > 64)
            return;
        const int height = block_size_high[block_size];
        DECLARE_ALIGNED(16, uint16_t, input[MAX_SB_SQUARE]);
        DECLARE_ALIGNED(16, int16_t, output_q3_ref_[MAX_SB_SQUARE]);
        DECLARE_ALIGNED(16, int16_t, output_q3_tst_[MAX_SB_SQUARE]);

        memset(output_q3_ref_, 1, sizeof(output_q3_ref_));
        memset(output_q3_tst_, 1, sizeof(output_q3_tst_));

        const int run_times = 100;
        for (int i = 0; i < run_times; ++i) {
            memset(input, 1, sizeof(input));
            for (int j = 0; j < MAX_SB_SQUARE; j++) {
                input[j] = rnd_.random();
            }

            cfl_luma_subsampling_420_hbd_c(
                input, width, output_q3_ref_, width, height);

            test_impl(input, width, output_q3_tst_, width, height);

            ASSERT_EQ(
                0,
                memcmp(output_q3_ref_, output_q3_tst_, sizeof(output_q3_ref_)));
        }
    }

  private:
    SVTRandom rnd_;
};

TEST_P(CflLumaSubsamplingHbdTest, MatchTest) {
    run_test();
}

INSTANTIATE_TEST_CASE_P(
    CFL_LUMA_SUBSAMPLING_HBD, CflLumaSubsamplingHbdTest,
    ::testing::Combine(::testing::Range(BLOCK_4X4, BlockSizeS_ALL),
                       ::testing::Values(cfl_luma_subsampling_420_hbd_avx2)));

}  // namespace
