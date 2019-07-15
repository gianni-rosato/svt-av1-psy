/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file weiner_convolve_test.cc
 *
 * @brief Unit test of weiner convolbe add source:
 * - av1_wiener_convolve_add_src_avx2
 * - av1_highbd_wiener_convolve_add_src_avx2
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
#include "EbDefinitions.h"
#include "EbRestoration.h"
#include "convolve.h"
#include "random.h"
#include "util.h"

/**
 * @brief Unit test of weiner convolbe add source:
 * - eb_av1_wiener_convolve_add_src_avx2
 * - eb_av1_highbd_wiener_convolve_add_src_avx2
 *
 * Test strategy:
 * Verify this assembly code by comparing with reference c implementation.
 * Feed the same data and check test output and reference output.
 * Define a templete class to handle the common process, and
 * declare sub class to handle different bitdepth and function types.
 *
 * Expected result:
 * Output from assemble functions should be the same with output from c.
 *
 * Test coverage:
 * Test cases:
 * BitDepth: 8bit, 10bit and 12bit
 */

namespace {
using std::make_tuple;
using svt_av1_test_tool::SVTRandom;

static const int w = 128, h = 128;

typedef void (*WienerConvolveFunc)(const uint8_t* src, ptrdiff_t src_stride,
                                   uint8_t* dst, ptrdiff_t dst_stride,
                                   const int16_t* filter_x, int x_step_q4,
                                   const int16_t* filter_y, int y_step_q4,
                                   int w, int h,
                                   const ConvolveParams* conv_params);
typedef void (*HbdWienerConvolveFunc)(const uint8_t* src, ptrdiff_t src_stride,
                                      uint8_t* dst, ptrdiff_t dst_stride,
                                      const int16_t* filter_x, int x_step_q4,
                                      const int16_t* filter_y, int y_step_q4,
                                      int w, int h,
                                      const ConvolveParams* conv_params,
                                      int32_t bps);
typedef std::tuple<int, int, WienerConvolveFunc> WienerConvolveParam;
typedef std::tuple<int, int, HbdWienerConvolveFunc, int> HbdWienerConvolveParam;

template <typename Sample, typename FuncType, typename ParamType>
class AV1WienerConvolveTest : public ::testing::TestWithParam<ParamType> {
  public:
    AV1WienerConvolveTest() : rnd_(16, false) {
        input_ = nullptr;
        output_tst_ = nullptr;
        output_ref_ = nullptr;
        out_w_ = 0;
        out_h_ = 0;
        func_tst_ = nullptr;
        bd_ = 8;
    }

    virtual ~AV1WienerConvolveTest() {
        if (input_) {
            delete[] input_;
            input_ = nullptr;
        }
        if (output_tst_) {
            delete[] output_tst_;
            output_tst_ = nullptr;
        }
        if (output_ref_) {
            delete[] output_ref_;
            output_ref_ = nullptr;
        }
        aom_clear_system_state();
    }

    void SetUp() override {
        rnd_.reset();
        malloc_data();
    }

  protected:
    void malloc_data() {
        input_ = new Sample[w * h];
        ASSERT_NE(input_, nullptr) << "create input buffer failed!";

        // The AVX2 convolve functions always write rows with widths that are
        // multiples of 16. So to avoid a buffer overflow, we may need to pad
        // rows to a multiple of 16.
        int output_n = ALIGN_POWER_OF_TWO(out_w_, 4) * out_h_;
        output_tst_ = new Sample[output_n];
        ASSERT_NE(output_tst_, nullptr) << "create test output buffer failed!";
        output_ref_ = new Sample[output_n];
        ASSERT_NE(output_ref_, nullptr)
            << "create refernece output buffer failed!";
    }

    void prepare_random_data() {
        for (int i = 0; i < h; ++i)
            for (int j = 0; j < w; ++j)
                input_[i * w + j] = (Sample)rnd_.random() & ((1 << bd_) - 1);
    }

    // Generate a random pair of filter kernels, using the ranges
    // of possible values from the loop-restoration experiment
    void generate_kernels(InterpKernel hkernel, InterpKernel vkernel,
                          int kernel_type = 2) {
        if (kernel_type == 0) {
            // Low possible values for filter coefficients
            hkernel[0] = hkernel[6] = vkernel[0] = vkernel[6] =
                WIENER_FILT_TAP0_MINV;
            hkernel[1] = hkernel[5] = vkernel[1] = vkernel[5] =
                WIENER_FILT_TAP1_MINV;
            hkernel[2] = hkernel[4] = vkernel[2] = vkernel[4] =
                WIENER_FILT_TAP2_MINV;
            hkernel[3] = vkernel[3] =
                -2 * (hkernel[0] + hkernel[1] + hkernel[2]);
            hkernel[7] = vkernel[7] = 0;
        } else if (kernel_type == 1) {
            // Max possible values for filter coefficients
            hkernel[0] = hkernel[6] = vkernel[0] = vkernel[6] =
                WIENER_FILT_TAP0_MAXV;
            hkernel[1] = hkernel[5] = vkernel[1] = vkernel[5] =
                WIENER_FILT_TAP1_MAXV;
            hkernel[2] = hkernel[4] = vkernel[2] = vkernel[4] =
                WIENER_FILT_TAP2_MAXV;
            hkernel[3] = vkernel[3] =
                -2 * (hkernel[0] + hkernel[1] + hkernel[2]);
            hkernel[7] = vkernel[7] = 0;
        } else {
            // Randomly generated values for filter coefficients
            hkernel[0] = hkernel[6] = WIENER_FILT_TAP0_MINV +
                                      pseudo_uniform(WIENER_FILT_TAP0_MAXV + 1 -
                                                     WIENER_FILT_TAP0_MINV);
            hkernel[1] = hkernel[5] = WIENER_FILT_TAP1_MINV +
                                      pseudo_uniform(WIENER_FILT_TAP1_MAXV + 1 -
                                                     WIENER_FILT_TAP1_MINV);
            hkernel[2] = hkernel[4] = WIENER_FILT_TAP2_MINV +
                                      pseudo_uniform(WIENER_FILT_TAP2_MAXV + 1 -
                                                     WIENER_FILT_TAP2_MINV);
            hkernel[3] = -2 * (hkernel[0] + hkernel[1] + hkernel[2]);
            hkernel[7] = 0;

            vkernel[0] = vkernel[6] = WIENER_FILT_TAP0_MINV +
                                      pseudo_uniform(WIENER_FILT_TAP0_MAXV + 2 -
                                                     WIENER_FILT_TAP0_MINV);
            vkernel[1] = vkernel[5] = WIENER_FILT_TAP1_MINV +
                                      pseudo_uniform(WIENER_FILT_TAP1_MAXV + 2 -
                                                     WIENER_FILT_TAP1_MINV);
            vkernel[2] = vkernel[4] = WIENER_FILT_TAP2_MINV +
                                      pseudo_uniform(WIENER_FILT_TAP2_MAXV + 2 -
                                                     WIENER_FILT_TAP2_MINV);
            vkernel[3] = -2 * (vkernel[0] + vkernel[1] + vkernel[2]);
            vkernel[7] = 0;
        }
    }

    int pseudo_uniform(const int range) {
        return rnd_.random() % (range + 1);
    }

    virtual void run_and_check(const InterpKernel& hkernel,
                               const InterpKernel& vkernel,
                               const ConvolveParams& params) = 0;

    virtual void run_random_test(const int test_times) {
        // Generate random filter kernels
        DECLARE_ALIGNED(16, InterpKernel, hkernel);
        DECLARE_ALIGNED(16, InterpKernel, vkernel);

        const ConvolveParams conv_params = get_conv_params_wiener(bd_);
        for (int kernel_type = 0; kernel_type < 3; kernel_type++) {
            generate_kernels(hkernel, vkernel, kernel_type);
            for (int i = 0;
                 i < test_times && !::testing::Test::HasFatalFailure();
                 ++i) {
                prepare_random_data();
                run_and_check(hkernel, vkernel, conv_params);
            }
        }
    }

  protected:
    FuncType func_tst_;  /**< function to test */
    uint32_t out_w_;     /**< width of output */
    uint32_t out_h_;     /**< height of output */
    uint32_t bd_;        /**< bit-depth */
    Sample* input_;      /**< buffer of input data */
    Sample* output_tst_; /**< buffer of test output data */
    Sample* output_ref_; /**< buffer of reference output data */
    SVTRandom rnd_;      /**< random tools*/
};

class AV1WienerConvolveLbdTest
    : public AV1WienerConvolveTest<uint8_t, WienerConvolveFunc,
                                   WienerConvolveParam> {
  protected:
    void SetUp() override {
        out_w_ = TEST_GET_PARAM(0);
        out_h_ = TEST_GET_PARAM(1);
        func_tst_ = TEST_GET_PARAM(2);
        AV1WienerConvolveTest::SetUp();
    }

    void run_and_check(const InterpKernel& hkernel, const InterpKernel& vkernel,
                       const ConvolveParams& params) override {
        uint8_t* input = input_;
        // Choose random locations within the source block
        int offset_r = 3 + pseudo_uniform(h - out_h_ - 7);
        int offset_c = 3 + pseudo_uniform(w - out_w_ - 7);
        eb_av1_wiener_convolve_add_src_c(input + offset_r * w + offset_c,
                                      w,
                                      output_ref_,
                                      out_w_,
                                      hkernel,
                                      16,
                                      vkernel,
                                      16,
                                      out_w_,
                                      out_h_,
                                      &params);
        func_tst_(input + offset_r * w + offset_c,
                  w,
                  output_tst_,
                  out_w_,
                  hkernel,
                  16,
                  vkernel,
                  16,
                  out_w_,
                  out_h_,
                  &params);

        for (uint32_t j = 0; j < out_w_ * out_h_; ++j)
            ASSERT_EQ(output_tst_[j], output_ref_[j])
                << "Pixel mismatch at index " << j << " = (" << (j % out_w_)
                << ", " << (j / out_w_) << ")";
    }

  public:
    static ::testing::internal::ParamGenerator<WienerConvolveParam> BuildParams(
        WienerConvolveFunc func) {
        const WienerConvolveParam params[] = {make_tuple(8, 8, func),
                                              make_tuple(8, 4, func),
                                              make_tuple(64, 24, func),
                                              make_tuple(64, 64, func),
                                              make_tuple(64, 56, func),
                                              make_tuple(32, 8, func),
                                              make_tuple(32, 28, func),
                                              make_tuple(32, 32, func),
                                              make_tuple(16, 34, func),
                                              make_tuple(32, 34, func),
                                              make_tuple(64, 34, func),
                                              make_tuple(8, 17, func),
                                              make_tuple(16, 17, func),
                                              make_tuple(32, 17, func)};
        return ::testing::ValuesIn(params);
    }
};

class AV1WienerConvolveHbdTest
    : public AV1WienerConvolveTest<uint16_t, HbdWienerConvolveFunc,
                                   HbdWienerConvolveParam> {
  protected:
    void SetUp() override {
        out_w_ = TEST_GET_PARAM(0);
        out_h_ = TEST_GET_PARAM(1);
        func_tst_ = TEST_GET_PARAM(2);
        bd_ = TEST_GET_PARAM(3);
        AV1WienerConvolveTest::SetUp();
    }

    void run_and_check(const InterpKernel& hkernel, const InterpKernel& vkernel,
                       const ConvolveParams& params) override {
        uint8_t* input = CONVERT_TO_BYTEPTR(input_);
        uint8_t* out_tst = CONVERT_TO_BYTEPTR(output_tst_);
        uint8_t* out_ref = CONVERT_TO_BYTEPTR(output_ref_);
        // Choose random locations within the source block
        int offset_r = 3 + pseudo_uniform(h - out_h_ - 7);
        int offset_c = 3 + pseudo_uniform(w - out_w_ - 7);
        eb_av1_highbd_wiener_convolve_add_src_c(input + offset_r * w + offset_c,
                                             w,
                                             out_ref,
                                             out_w_,
                                             hkernel,
                                             16,
                                             vkernel,
                                             16,
                                             out_w_,
                                             out_h_,
                                             &params,
                                             bd_);
        func_tst_(input + offset_r * w + offset_c,
                  w,
                  out_tst,
                  out_w_,
                  hkernel,
                  16,
                  vkernel,
                  16,
                  out_w_,
                  out_h_,
                  &params,
                  bd_);

        for (uint32_t j = 0; j < out_w_ * out_h_; ++j)
            ASSERT_EQ(output_tst_[j], output_ref_[j])
                << "Pixel mismatch at index " << j << " = (" << (j % out_w_)
                << ", " << (j / out_w_) << ")";
    }

  public:
    static ::testing::internal::ParamGenerator<HbdWienerConvolveParam>
    BuildParams(HbdWienerConvolveFunc func) {
        const HbdWienerConvolveParam params[] = {make_tuple(8, 8, func, 8),
                                                 make_tuple(64, 64, func, 8),
                                                 make_tuple(32, 8, func, 8),
                                                 make_tuple(8, 8, func, 10),
                                                 make_tuple(64, 64, func, 10),
                                                 make_tuple(32, 8, func, 10),
                                                 make_tuple(8, 8, func, 12),
                                                 make_tuple(64, 64, func, 12),
                                                 make_tuple(32, 8, func, 12)};
        return ::testing::ValuesIn(params);
    }
};

TEST_P(AV1WienerConvolveLbdTest, run_random_test) {
    run_random_test(1000);
}

INSTANTIATE_TEST_CASE_P(
    AV1, AV1WienerConvolveLbdTest,
    AV1WienerConvolveLbdTest::BuildParams(eb_av1_wiener_convolve_add_src_avx2));

TEST_P(AV1WienerConvolveHbdTest, run_random_test) {
    run_random_test(1000);
}

INSTANTIATE_TEST_CASE_P(AV1, AV1WienerConvolveHbdTest,
                        AV1WienerConvolveHbdTest::BuildParams(
                            eb_av1_highbd_wiener_convolve_add_src_avx2));

}  // namespace
