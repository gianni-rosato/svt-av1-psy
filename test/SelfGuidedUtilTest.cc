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
 * @file pixel_proj_err_test.cc
 *
 * @brief Unit test of project-related test in selfguided filter:
 *
 * - av1_lowbd_pixel_proj_error_avx2
 * - av1_highbd_pixel_proj_error_avx2
 * - svt_get_proj_subspace_avx2
 *
 * @author Cidana-Edmond, Cidana-Wenyao
 *
 ******************************************************************************/

#include "gtest/gtest.h"
#include "aom_dsp_rtcd.h"
#include "EbDefinitions.h"
#include "EbUtility.h"
#include "EbRestoration.h"
#include "EbUnitTestUtility.h"
#include "random.h"
#include "util.h"

#define MAX_DATA_BLOCK 384

namespace {
using std::make_tuple;
using svt_av1_test_tool::SVTRandom;

static const int min_test_times = 10;

typedef int64_t (*PixelProjFunc)(const uint8_t *src8, int32_t width,
                                 int32_t height, int32_t src_stride,
                                 const uint8_t *dat8, int32_t dat_stride,
                                 int32_t *flt0, int32_t flt0_stride,
                                 int32_t *flt1, int32_t flt1_stride,
                                 int32_t xq[2], const SgrParamsType *params);

typedef std::tuple<const PixelProjFunc, const PixelProjFunc>
    PixelProjErrorTestParam;

/**
 * @brief Unit test for pixel projection error:
 * - av1_lowbd_pixel_proj_error_avx2
 * - av1_highbd_pixel_proj_error_avx2
 *
 * Test strategy:
 * Verify this assembly code by comparing with reference c implementation.
 * Feed the same random data and check test output and reference output.
 * Define a template class to handle the common process, and
 * declare sub class to handle different bitdepth.
 *
 * Expected result:
 * Output from assemble functions should be the same with output from c.
 *
 * Test coverage:
 * Test cases:
 * input value: Fill with random values
 * test mode: fixed block size, random block size and extreme data check
 *
 */
template <typename Sample>
class PixelProjErrorTest
    : public ::testing::TestWithParam<PixelProjErrorTestParam> {
  public:
    PixelProjErrorTest()
        : rnd8_(8, false),
          rnd16_(16, false),
          rnd15s_(15, true),
          rnd_blk_size_(1, MAX_DATA_BLOCK) {
        tst_func_ = TEST_GET_PARAM(0);
        ref_func_ = TEST_GET_PARAM(1);
        src_ = nullptr;
        dgd_ = nullptr;
        flt0_ = nullptr;
        flt1_ = nullptr;
    }

    virtual void SetUp() {
        src_ = (Sample *)(svt_aom_malloc(MAX_DATA_BLOCK * MAX_DATA_BLOCK *
                                        sizeof(*src_)));
        ASSERT_NE(src_, nullptr);
        dgd_ = (Sample *)(svt_aom_malloc(MAX_DATA_BLOCK * MAX_DATA_BLOCK *
                                        sizeof(*dgd_)));
        ASSERT_NE(dgd_, nullptr);
        flt0_ = (int32_t *)(svt_aom_malloc(MAX_DATA_BLOCK * MAX_DATA_BLOCK *
                                          sizeof(*flt0_)));
        ASSERT_NE(flt0_, nullptr);
        flt1_ = (int32_t *)(svt_aom_malloc(MAX_DATA_BLOCK * MAX_DATA_BLOCK *
                                          sizeof(*flt1_)));
        ASSERT_NE(flt1_, nullptr);
    }

    virtual void TearDown() {
        svt_aom_free(src_);
        svt_aom_free(dgd_);
        svt_aom_free(flt0_);
        svt_aom_free(flt1_);
    }

    virtual void prepare_random_data() = 0;
    virtual void prepare_extreme_data() = 0;
    void run_and_check_data(const int index, const int fixed_size) {
        const int dgd_stride = MAX_DATA_BLOCK;
        const int src_stride = MAX_DATA_BLOCK;
        const int flt0_stride = MAX_DATA_BLOCK;
        const int flt1_stride = MAX_DATA_BLOCK;
        int h_end = fixed_size;
        int v_end = fixed_size;
        bool is_fixed_size = true;
        if (fixed_size == 0) {
            h_end = rnd_blk_size_.random();
            v_end = rnd_blk_size_.random();
            is_fixed_size = false;
        }

        int xq[2];
        xq[0] = rnd8_.random() % (1 << SGRPROJ_PRJ_BITS);
        xq[1] = rnd8_.random() % (1 << SGRPROJ_PRJ_BITS);
        SgrParamsType params;
        params.r[0] =
            !is_fixed_size ? (rnd8_.random() % MAX_RADIUS) : (index % 2);
        params.r[1] =
            !is_fixed_size ? (rnd8_.random() % MAX_RADIUS) : (index / 2);
        params.s[0] =
            !is_fixed_size ? (rnd8_.random() % MAX_RADIUS) : (index % 2);
        params.s[1] =
            !is_fixed_size ? (rnd8_.random() % MAX_RADIUS) : (index / 2);
        uint8_t *dgd =
            (sizeof(*dgd_) == 2) ? (CONVERT_TO_BYTEPTR(dgd_)) : (uint8_t *)dgd_;
        uint8_t *src =
            (sizeof(*src_) == 2) ? (CONVERT_TO_BYTEPTR(src_)) : (uint8_t *)src_;

        int64_t err_ref = ref_func_(src,
                                    h_end,
                                    v_end,
                                    src_stride,
                                    dgd,
                                    dgd_stride,
                                    flt0_,
                                    flt0_stride,
                                    flt1_,
                                    flt1_stride,
                                    xq,
                                    &params);
        int64_t err_test = tst_func_(src,
                                     h_end,
                                     v_end,
                                     src_stride,
                                     dgd,
                                     dgd_stride,
                                     flt0_,
                                     flt0_stride,
                                     flt1_,
                                     flt1_stride,
                                     xq,
                                     &params);
        ASSERT_EQ(err_ref, err_test);
    }

    virtual void run_random_test(const int run_times,
                                 const bool is_fixed_size) {
        const int iters = AOMMAX(run_times, min_test_times);
        for (int iter = 0; iter < iters && !HasFatalFailure(); ++iter) {
            prepare_random_data();
            run_and_check_data(iter, is_fixed_size ? 128 : 0);
        }
    }

    virtual void run_speed_test(const int size, const int r0, const int r1) {
        prepare_random_data();

        const int dgd_stride = MAX_DATA_BLOCK;
        const int src_stride = MAX_DATA_BLOCK;
        const int flt0_stride = MAX_DATA_BLOCK;
        const int flt1_stride = MAX_DATA_BLOCK;
        int h_end = size;
        int v_end = size;
        int xq[2];
        int64_t err_ref, err_test;
        double time_c, time_o;
        uint64_t start_time_seconds, start_time_useconds;
        uint64_t middle_time_seconds, middle_time_useconds;
        uint64_t finish_time_seconds, finish_time_useconds;

        xq[0] = rnd8_.random() % (1 << SGRPROJ_PRJ_BITS);
        xq[1] = rnd8_.random() % (1 << SGRPROJ_PRJ_BITS);
        SgrParamsType params;
        params.r[0] = r0;
        params.r[1] = r1;
        uint8_t *dgd =
            (sizeof(*dgd_) == 2) ? (CONVERT_TO_BYTEPTR(dgd_)) : (uint8_t *)dgd_;
        uint8_t *src =
            (sizeof(*src_) == 2) ? (CONVERT_TO_BYTEPTR(src_)) : (uint8_t *)src_;

        const uint64_t num_loop = 1000000;

        svt_av1_get_time(&start_time_seconds, &start_time_useconds);

        for (uint64_t i = 0; i < num_loop; i++) {
            err_ref = ref_func_(src,
                                h_end,
                                v_end,
                                src_stride,
                                dgd,
                                dgd_stride,
                                flt0_,
                                flt0_stride,
                                flt1_,
                                flt1_stride,
                                xq,
                                &params);
        }

        svt_av1_get_time(&middle_time_seconds, &middle_time_useconds);

        for (uint64_t i = 0; i < num_loop; i++) {
            err_test = tst_func_(src,
                                 h_end,
                                 v_end,
                                 src_stride,
                                 dgd,
                                 dgd_stride,
                                 flt0_,
                                 flt0_stride,
                                 flt1_,
                                 flt1_stride,
                                 xq,
                                 &params);
        }

        svt_av1_get_time(&finish_time_seconds, &finish_time_useconds);

        ASSERT_EQ(err_ref, err_test);

        time_c = svt_av1_compute_overall_elapsed_time_ms(start_time_seconds,
                                                         start_time_useconds,
                                                         middle_time_seconds,
                                                         middle_time_useconds);
        time_o = svt_av1_compute_overall_elapsed_time_ms(middle_time_seconds,
                                                         middle_time_useconds,
                                                         finish_time_seconds,
                                                         finish_time_useconds);

        printf("Average Nanoseconds per Function Call\n");
        printf(
            "   svt_av1_lowbd_pixel_proj_error_c(size: %d, r0: %d, r1: %d)   : "
            "%6.2f\n",
            size,
            r0,
            r1,
            1000000 * time_c / num_loop);
        printf(
            "   svt_av1_lowbd_pixel_proj_error_optsize: %d, r0: %d, r1: %d) : "
            "%6.2f   (Comparison: "
            "%5.2fx)\n",
            size,
            r0,
            r1,
            1000000 * time_o / num_loop,
            time_c / time_o);
    }

    virtual void run_extreme_test() {
        const int iters = min_test_times;
        for (int iter = 0; iter < iters && !HasFatalFailure(); ++iter) {
            prepare_extreme_data();
            run_and_check_data(iter, 192);
        }
    }

  protected:
    PixelProjFunc tst_func_;
    PixelProjFunc ref_func_;
    Sample *src_;
    Sample *dgd_;
    int32_t *flt0_;
    int32_t *flt1_;

    SVTRandom rnd8_;
    SVTRandom rnd16_;
    SVTRandom rnd15s_;
    SVTRandom rnd_blk_size_;
};

class PixelProjErrorLbdTest : public PixelProjErrorTest<uint8_t> {
    void prepare_random_data() override {
        for (int i = 0; i < MAX_DATA_BLOCK * MAX_DATA_BLOCK; ++i) {
            dgd_[i] = rnd8_.random();
            src_[i] = rnd8_.random();
            flt0_[i] = rnd15s_.random();
            flt1_[i] = rnd15s_.random();
        }
    }

    void prepare_extreme_data() override {
        for (int i = 0; i < MAX_DATA_BLOCK * MAX_DATA_BLOCK; ++i) {
            dgd_[i] = 0;
            src_[i] = 255;
            flt0_[i] = rnd15s_.random();
            flt1_[i] = rnd15s_.random();
        }
    }
};

TEST_P(PixelProjErrorLbdTest, MatchTestWithRandomValue) {
    run_random_test(50, true);
}
TEST_P(PixelProjErrorLbdTest, MatchTestWithRandomSizeAndValue) {
    run_random_test(50, false);
}
TEST_P(PixelProjErrorLbdTest, MatchTestWithExtremeValue) {
    run_extreme_test();
}
TEST_P(PixelProjErrorLbdTest, DISABLED_SpeedTest) {
    run_speed_test(256, 1, 1);
    run_speed_test(256, 1, 0);
    run_speed_test(256, 0, 0);
}

static const PixelProjErrorTestParam lbd_test_vector[] = {
    make_tuple(svt_av1_lowbd_pixel_proj_error_avx2,
               svt_av1_lowbd_pixel_proj_error_c),
#ifndef NON_AVX512_SUPPORT
    make_tuple(svt_av1_lowbd_pixel_proj_error_avx512,
               svt_av1_lowbd_pixel_proj_error_c)
#endif
};

INSTANTIATE_TEST_CASE_P(RST, PixelProjErrorLbdTest,
                        ::testing::ValuesIn(lbd_test_vector));

class PixelProjErrorHbdTest : public PixelProjErrorTest<uint16_t> {
  protected:
    PixelProjErrorHbdTest() : rnd12_(12, false) {
    }

    void prepare_random_data() override {
        for (int i = 0; i < MAX_DATA_BLOCK * MAX_DATA_BLOCK; ++i) {
            dgd_[i] = rnd12_.random();
            src_[i] = rnd12_.random();
            flt0_[i] = rnd15s_.random();
            flt1_[i] = rnd15s_.random();
        }
    }

    void prepare_extreme_data() override {
        for (int i = 0; i < MAX_DATA_BLOCK * MAX_DATA_BLOCK; ++i) {
            dgd_[i] = 0;
            src_[i] = (1 << 12) - 1;
            flt0_[i] = rnd15s_.random();
            flt1_[i] = rnd15s_.random();
        }
    }

  private:
    SVTRandom rnd12_;
};

TEST_P(PixelProjErrorHbdTest, MatchTestWithRandomValue) {
    run_random_test(50, true);
}
TEST_P(PixelProjErrorHbdTest, MatchTestWithRandomSizeAndValue) {
    run_random_test(50, false);
}
TEST_P(PixelProjErrorHbdTest, MatchTestWithExtremeValue) {
    run_extreme_test();
}

static const PixelProjErrorTestParam hbd_test_vector[] = {make_tuple(
    svt_av1_highbd_pixel_proj_error_avx2, svt_av1_highbd_pixel_proj_error_c)};

INSTANTIATE_TEST_CASE_P(RST, PixelProjErrorHbdTest,
                        ::testing::ValuesIn(hbd_test_vector));

// test svt_get_proj_subspace
TEST(SelfGuidedToolsTest, GetProjSubspaceMatchTest) {
    const int32_t pu_width = RESTORATION_PROC_UNIT_SIZE;
    const int32_t pu_height = RESTORATION_PROC_UNIT_SIZE;
    const int32_t width = 256, height = 256, stride = 288, out_stride = 288;
    const int NUM_ITERS = 2000;
    int i, j, k;

    uint8_t *input_ = (uint8_t *)svt_aom_memalign(
        32, stride * (height + 32) * sizeof(uint8_t));
    uint8_t *output_ = (uint8_t *)svt_aom_memalign(
        32, out_stride * (height + 32) * sizeof(uint8_t));
    int32_t *tmpbuf = (int32_t *)svt_aom_memalign(32, RESTORATION_TMPBUF_SIZE);
    uint8_t *input = input_ + stride * 16 + 16;
    uint8_t *output = output_ + out_stride * 16 + 16;
    int32_t *flt0 = tmpbuf;
    int32_t *flt1 = flt0 + RESTORATION_UNITPELS_MAX;
    int32_t flt_stride = ((width + 7) & ~7) + 8;

    // check all the sg params
    SVTRandom rnd(8, false);
    for (int iter = 0; iter < NUM_ITERS; ++iter) {
        // prepare src data and recon data
        for (i = -16; i < height + 16; ++i) {
            for (j = -16; j < width + 16; ++j) {
                input[i * stride + j] = rnd.random();
                if (iter == 0)
                    output[i * stride + j] = input[i * stride + j];
                else
                    output[i * stride + j] = rnd.random();
            }
        }

        for (int32_t ep = 0; ep < SGRPROJ_PARAMS; ++ep) {
            // apply selfguided filter to get A and b
            for (k = 0; k < height; k += pu_height) {
                for (j = 0; j < width; j += pu_width) {
                    int32_t w = AOMMIN(pu_width, width - j);
                    int32_t h = AOMMIN(pu_height, height - k);
                    uint8_t *output_p = output + k * out_stride + j;
                    int32_t *flt0_p = flt0 + k * flt_stride + j;
                    int32_t *flt1_p = flt1 + k * flt_stride + j;
                    assert(w * h <= RESTORATION_UNITPELS_MAX);

                    svt_av1_selfguided_restoration_avx2(output_p,
                                                        w,
                                                        h,
                                                        out_stride,
                                                        flt0_p,
                                                        flt1_p,
                                                        flt_stride,
                                                        ep,
                                                        8,
                                                        0);
                }
            }

            aom_clear_system_state();
            int32_t xqd_c[2] = {0};
            int32_t xqd_asm[2] = {0};
            const SgrParamsType *const params = &eb_sgr_params[ep];
            svt_get_proj_subspace_c(input,
                                    width,
                                    height,
                                    stride,
                                    output,
                                    out_stride,
                                    0,
                                    flt0,
                                    flt_stride,
                                    flt1,
                                    flt_stride,
                                    xqd_c,
                                    params);
            svt_get_proj_subspace_avx2(input,
                                       width,
                                       height,
                                       stride,
                                       output,
                                       out_stride,
                                       0,
                                       flt0,
                                       flt_stride,
                                       flt1,
                                       flt_stride,
                                       xqd_asm,
                                       params);
            ASSERT_EQ(xqd_c[0], xqd_asm[0])
                << "xqd_asm[0] does not match with xqd_asm[0] with iter "
                << iter << " ep " << ep;
            ASSERT_EQ(xqd_c[1], xqd_asm[1])
                << "xqd_asm[1] does not match with xqd_asm[1] with iter "
                << iter << " ep " << ep;
        }
    }

    svt_aom_free(input_);
    svt_aom_free(output_);
    svt_aom_free(tmpbuf);
}
}  // namespace
