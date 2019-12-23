/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <ctime>

#include "gtest/gtest.h"
#include "acm_random.h"
#include "aom_dsp_rtcd.h"
#include "EbDefinitions.h"
#include "EbRestoration.h"
#include "EbRestorationPick.h"
#include "EbUnitTestUtility.h"
#include "EbUtility.h"
#include "util.h"

namespace {

using libaom_test::ACMRandom;
using std::make_tuple;
using std::tuple;

typedef void (*SgrFunc)(const uint8_t *dat8, int32_t width, int32_t height,
                        int32_t stride, int32_t eps, const int32_t *xqd,
                        uint8_t *dst8, int32_t dst_stride, int32_t *tmpbuf,
                        int32_t bit_depth, int32_t highbd);

// Test parameter list:
//  <tst_fun_>
typedef tuple<SgrFunc> FilterTestParam;

class AV1SelfguidedFilterTest
    : public ::testing::TestWithParam<FilterTestParam> {
  public:
    virtual ~AV1SelfguidedFilterTest() {
    }
    virtual void SetUp() {
    }

    virtual void TearDown() {
        aom_clear_system_state();
    }

  protected:
    void RunSpeedTest() {
        tst_fun_ = TEST_GET_PARAM(0);
        const int32_t pu_width = RESTORATION_PROC_UNIT_SIZE;
        const int32_t pu_height = RESTORATION_PROC_UNIT_SIZE;
        const int32_t width = 256, height = 256, stride = 288, out_stride = 288;
        const int NUM_ITERS = 2000;
        int i, j, k;

        uint8_t *input_ = (uint8_t *)eb_aom_memalign(
            32, stride * (height + 32) * sizeof(uint8_t));
        uint8_t *output_ = (uint8_t *)eb_aom_memalign(
            32, out_stride * (height + 32) * sizeof(uint8_t));
        int32_t *tmpbuf = (int32_t *)eb_aom_memalign(32, RESTORATION_TMPBUF_SIZE);
        uint8_t *input = input_ + stride * 16 + 16;
        uint8_t *output = output_ + out_stride * 16 + 16;

        ACMRandom rnd(ACMRandom::DeterministicSeed());

        for (i = -16; i < height + 16; ++i)
            for (j = -16; j < width + 16; ++j)
                input[i * stride + j] = rnd.Rand16() & 0xFF;

        int32_t xqd[2] = {
            SGRPROJ_PRJ_MIN0 +
                rnd.PseudoUniform(SGRPROJ_PRJ_MAX0 + 1 - SGRPROJ_PRJ_MIN0),
            SGRPROJ_PRJ_MIN1 +
                rnd.PseudoUniform(SGRPROJ_PRJ_MAX1 + 1 - SGRPROJ_PRJ_MIN1)};
        // Fix a parameter set, since the speed depends slightly on r.
        // Change this to test different combinations of values of r.
        int32_t eps = 15;
        double ref_time, tst_time;
        uint64_t start_time_seconds, start_time_useconds;
        uint64_t middle_time_seconds, middle_time_useconds;
        uint64_t finish_time_seconds, finish_time_useconds;

        eb_start_time(&start_time_seconds, &start_time_useconds);
        for (i = 0; i < NUM_ITERS; ++i) {
            for (k = 0; k < height; k += pu_height)
                for (j = 0; j < width; j += pu_width) {
                    int32_t w = AOMMIN(pu_width, width - j);
                    int32_t h = AOMMIN(pu_height, height - k);
                    uint8_t *input_p = input + k * stride + j;
                    uint8_t *output_p = output + k * out_stride + j;
                    eb_apply_selfguided_restoration_c(input_p,
                                                   w,
                                                   h,
                                                   stride,
                                                   eps,
                                                   xqd,
                                                   output_p,
                                                   out_stride,
                                                   tmpbuf,
                                                   8,
                                                   0);
                }
        }
        eb_start_time(&middle_time_seconds, &middle_time_useconds);

        for (i = 0; i < NUM_ITERS; ++i) {
            for (k = 0; k < height; k += pu_height)
                for (j = 0; j < width; j += pu_width) {
                    int32_t w = AOMMIN(pu_width, width - j);
                    int32_t h = AOMMIN(pu_height, height - k);
                    uint8_t *input_p = input + k * stride + j;
                    uint8_t *output_p = output + k * out_stride + j;
                    tst_fun_(input_p,
                             w,
                             h,
                             stride,
                             eps,
                             xqd,
                             output_p,
                             out_stride,
                             tmpbuf,
                             8,
                             0);
                }
        }

        eb_start_time(&finish_time_seconds, &finish_time_useconds);
        eb_compute_overall_elapsed_time_ms(start_time_seconds,
                                      start_time_useconds,
                                      middle_time_seconds,
                                      middle_time_useconds,
                                      &ref_time);
        eb_compute_overall_elapsed_time_ms(middle_time_seconds,
                                      middle_time_useconds,
                                      finish_time_seconds,
                                      finish_time_useconds,
                                      &tst_time);

        std::cout << "[          ] C time = " << ref_time / 1000
                  << " ms, SIMD time = " << tst_time / 1000 << " ms\n";

        EXPECT_GT(ref_time, tst_time)
            << "Error: AV1SelfguidedFilterTest.SpeedTest, SIMD slower than C.\n"
            << "C time: " << ref_time << " us\n"
            << "SIMD time: " << tst_time << " us\n";

        eb_aom_free(input_);
        eb_aom_free(output_);
        eb_aom_free(tmpbuf);
    }

    void RunCorrectnessTest() {
        tst_fun_ = TEST_GET_PARAM(0);
        const int32_t pu_width = RESTORATION_PROC_UNIT_SIZE;
        const int32_t pu_height = RESTORATION_PROC_UNIT_SIZE;
        // Set the maximum width/height to test here. We actually test a small
        // range of sizes *up to* this size, so that we can check, eg.,
        // the behaviour on tiles which are not a multiple of 4 wide.
        const int32_t max_w = 260, max_h = 260, stride = 672, out_stride = 672;
        const int NUM_ITERS = 81;
        int i, j, k;

        uint8_t *input_ = (uint8_t *)eb_aom_memalign(
            32, stride * (max_h + 32) * sizeof(uint8_t));
        uint8_t *output_ = (uint8_t *)eb_aom_memalign(
            32, out_stride * (max_h + 32) * sizeof(uint8_t));
        uint8_t *output2_ = (uint8_t *)eb_aom_memalign(
            32, out_stride * (max_h + 32) * sizeof(uint8_t));
        int32_t *tmpbuf = (int32_t *)eb_aom_memalign(32, RESTORATION_TMPBUF_SIZE);

        uint8_t *input = input_ + stride * 16 + 16;
        uint8_t *output = output_ + out_stride * 16 + 16;
        uint8_t *output2 = output2_ + out_stride * 16 + 16;

        ACMRandom rnd(ACMRandom::DeterministicSeed());

        for (i = 0; i < NUM_ITERS; ++i) {
            for (j = -16; j < max_h + 16; ++j)
                for (k = -16; k < max_w + 16; ++k)
                    input[j * stride + k] = rnd.Rand16() & 0xFF;

            int32_t xqd[2] = {
                SGRPROJ_PRJ_MIN0 +
                    rnd.PseudoUniform(SGRPROJ_PRJ_MAX0 + 1 - SGRPROJ_PRJ_MIN0),
                SGRPROJ_PRJ_MIN1 +
                    rnd.PseudoUniform(SGRPROJ_PRJ_MAX1 + 1 - SGRPROJ_PRJ_MIN1)};
            int32_t eps = rnd.PseudoUniform(1 << SGRPROJ_PARAMS_BITS);

            // Test various tile sizes around 256x256
            int32_t test_w = max_w - (i / 9);
            int32_t test_h = max_h - (i % 9);

            for (k = 0; k < test_h; k += pu_height)
                for (j = 0; j < test_w; j += pu_width) {
                    int32_t w = AOMMIN(pu_width, test_w - j);
                    int32_t h = AOMMIN(pu_height, test_h - k);
                    uint8_t *input_p = input + k * stride + j;
                    uint8_t *output_p = output + k * out_stride + j;
                    uint8_t *output2_p = output2 + k * out_stride + j;
                    tst_fun_(input_p,
                             w,
                             h,
                             stride,
                             eps,
                             xqd,
                             output_p,
                             out_stride,
                             tmpbuf,
                             8,
                             0);
                    eb_apply_selfguided_restoration_c(input_p,
                                                   w,
                                                   h,
                                                   stride,
                                                   eps,
                                                   xqd,
                                                   output2_p,
                                                   out_stride,
                                                   tmpbuf,
                                                   8,
                                                   0);
                }

            for (j = 0; j < test_h; ++j)
                for (k = 0; k < test_w; ++k) {
                    ASSERT_EQ(output[j * out_stride + k],
                              output2[j * out_stride + k]);
                }
        }

        eb_aom_free(input_);
        eb_aom_free(output_);
        eb_aom_free(output2_);
        eb_aom_free(tmpbuf);
    }

  private:
    SgrFunc tst_fun_;
};

TEST_P(AV1SelfguidedFilterTest, DISABLED_SpeedTest) {
    RunSpeedTest();
}
TEST_P(AV1SelfguidedFilterTest, CorrectnessTest) {
    RunCorrectnessTest();
}

INSTANTIATE_TEST_CASE_P(
    AVX2, AV1SelfguidedFilterTest,
    ::testing::Values(make_tuple(eb_apply_selfguided_restoration_avx2)));

// Test parameter list:
//  <tst_fun_, bit_depth>
typedef tuple<SgrFunc, int32_t> HighbdFilterTestParam;

class AV1HighbdSelfguidedFilterTest
    : public ::testing::TestWithParam<HighbdFilterTestParam> {
  public:
    virtual ~AV1HighbdSelfguidedFilterTest() {
    }
    virtual void SetUp() {
    }

    virtual void TearDown() {
        aom_clear_system_state();
    }

  protected:
    void RunSpeedTest() {
        tst_fun_ = TEST_GET_PARAM(0);
        const int32_t pu_width = RESTORATION_PROC_UNIT_SIZE;
        const int32_t pu_height = RESTORATION_PROC_UNIT_SIZE;
        const int32_t width = 256, height = 256, stride = 288, out_stride = 288;
        const int NUM_ITERS = 2000;
        int i, j, k;
        int32_t bit_depth = TEST_GET_PARAM(1);
        int mask = (1 << bit_depth) - 1;

        uint16_t *input_ = (uint16_t *)eb_aom_memalign(
            32, stride * (height + 32) * sizeof(uint16_t));
        uint16_t *output_ = (uint16_t *)eb_aom_memalign(
            32, out_stride * (height + 32) * sizeof(uint16_t));
        int32_t *tmpbuf = (int32_t *)eb_aom_memalign(32, RESTORATION_TMPBUF_SIZE);
        uint16_t *input = input_ + stride * 16 + 16;
        uint16_t *output = output_ + out_stride * 16 + 16;

        ACMRandom rnd(ACMRandom::DeterministicSeed());

        for (i = -16; i < height + 16; ++i)
            for (j = -16; j < width + 16; ++j)
                input[i * stride + j] = rnd.Rand16() & mask;

        int32_t xqd[2] = {
            SGRPROJ_PRJ_MIN0 +
                rnd.PseudoUniform(SGRPROJ_PRJ_MAX0 + 1 - SGRPROJ_PRJ_MIN0),
            SGRPROJ_PRJ_MIN1 +
                rnd.PseudoUniform(SGRPROJ_PRJ_MAX1 + 1 - SGRPROJ_PRJ_MIN1)};
        // Fix a parameter set, since the speed depends slightly on r.
        // Change this to test different combinations of values of r.
        int32_t eps = 15;
        double ref_time, tst_time;
        uint64_t start_time_seconds, start_time_useconds;
        uint64_t middle_time_seconds, middle_time_useconds;
        uint64_t finish_time_seconds, finish_time_useconds;

        eb_start_time(&start_time_seconds, &start_time_useconds);
        for (i = 0; i < NUM_ITERS; ++i) {
            for (k = 0; k < height; k += pu_height)
                for (j = 0; j < width; j += pu_width) {
                    int32_t w = AOMMIN(pu_width, width - j);
                    int32_t h = AOMMIN(pu_height, height - k);
                    uint16_t *input_p = input + k * stride + j;
                    uint16_t *output_p = output + k * out_stride + j;
                    eb_apply_selfguided_restoration_c(CONVERT_TO_BYTEPTR(input_p),
                                                   w,
                                                   h,
                                                   stride,
                                                   eps,
                                                   xqd,
                                                   CONVERT_TO_BYTEPTR(output_p),
                                                   out_stride,
                                                   tmpbuf,
                                                   bit_depth,
                                                   1);
                }
        }
        eb_start_time(&middle_time_seconds, &middle_time_useconds);

        for (i = 0; i < NUM_ITERS; ++i) {
            for (k = 0; k < height; k += pu_height)
                for (j = 0; j < width; j += pu_width) {
                    int32_t w = AOMMIN(pu_width, width - j);
                    int32_t h = AOMMIN(pu_height, height - k);
                    uint16_t *input_p = input + k * stride + j;
                    uint16_t *output_p = output + k * out_stride + j;
                    tst_fun_(CONVERT_TO_BYTEPTR(input_p),
                             w,
                             h,
                             stride,
                             eps,
                             xqd,
                             CONVERT_TO_BYTEPTR(output_p),
                             out_stride,
                             tmpbuf,
                             bit_depth,
                             1);
                }
        }

        eb_start_time(&finish_time_seconds, &finish_time_useconds);
        eb_compute_overall_elapsed_time_ms(start_time_seconds,
                                      start_time_useconds,
                                      middle_time_seconds,
                                      middle_time_useconds,
                                      &ref_time);
        eb_compute_overall_elapsed_time_ms(middle_time_seconds,
                                      middle_time_useconds,
                                      finish_time_seconds,
                                      finish_time_useconds,
                                      &tst_time);

        std::cout << "[          ] C time = " << ref_time / 1000
                  << " ms, SIMD time = " << tst_time / 1000 << " ms\n";

        EXPECT_GT(ref_time, tst_time)
            << "Error: AV1HighbdSelfguidedFilterTest.SpeedTest, SIMD slower "
               "than "
               "C.\n"
            << "C time: " << ref_time << " us\n"
            << "SIMD time: " << tst_time << " us\n";

        eb_aom_free(input_);
        eb_aom_free(output_);
        eb_aom_free(tmpbuf);
    }

    void RunCorrectnessTest() {
        tst_fun_ = TEST_GET_PARAM(0);
        const int32_t pu_width = RESTORATION_PROC_UNIT_SIZE;
        const int32_t pu_height = RESTORATION_PROC_UNIT_SIZE;
        // Set the maximum width/height to test here. We actually test a small
        // range of sizes *up to* this size, so that we can check, eg.,
        // the behaviour on tiles which are not a multiple of 4 wide.
        const int32_t max_w = 260, max_h = 260, stride = 672, out_stride = 672;
        const int NUM_ITERS = 81;
        int i, j, k;
        int32_t bit_depth = TEST_GET_PARAM(1);
        int mask = (1 << bit_depth) - 1;

        uint16_t *input_ = (uint16_t *)eb_aom_memalign(
            32, stride * (max_h + 32) * sizeof(uint16_t));
        uint16_t *output_ = (uint16_t *)eb_aom_memalign(
            32, out_stride * (max_h + 32) * sizeof(uint16_t));
        uint16_t *output2_ = (uint16_t *)eb_aom_memalign(
            32, out_stride * (max_h + 32) * sizeof(uint16_t));
        int32_t *tmpbuf = (int32_t *)eb_aom_memalign(32, RESTORATION_TMPBUF_SIZE);

        uint16_t *input = input_ + stride * 16 + 16;
        uint16_t *output = output_ + out_stride * 16 + 16;
        uint16_t *output2 = output2_ + out_stride * 16 + 16;

        ACMRandom rnd(ACMRandom::DeterministicSeed());

        for (i = 0; i < NUM_ITERS; ++i) {
            for (j = -16; j < max_h + 16; ++j)
                for (k = -16; k < max_w + 16; ++k)
                    input[j * stride + k] = rnd.Rand16() & mask;

            int32_t xqd[2] = {
                SGRPROJ_PRJ_MIN0 +
                    rnd.PseudoUniform(SGRPROJ_PRJ_MAX0 + 1 - SGRPROJ_PRJ_MIN0),
                SGRPROJ_PRJ_MIN1 +
                    rnd.PseudoUniform(SGRPROJ_PRJ_MAX1 + 1 - SGRPROJ_PRJ_MIN1)};
            int32_t eps = rnd.PseudoUniform(1 << SGRPROJ_PARAMS_BITS);

            // Test various tile sizes around 256x256
            int32_t test_w = max_w - (i / 9);
            int32_t test_h = max_h - (i % 9);

            for (k = 0; k < test_h; k += pu_height)
                for (j = 0; j < test_w; j += pu_width) {
                    int32_t w = AOMMIN(pu_width, test_w - j);
                    int32_t h = AOMMIN(pu_height, test_h - k);
                    uint16_t *input_p = input + k * stride + j;
                    uint16_t *output_p = output + k * out_stride + j;
                    uint16_t *output2_p = output2 + k * out_stride + j;
                    tst_fun_(CONVERT_TO_BYTEPTR(input_p),
                             w,
                             h,
                             stride,
                             eps,
                             xqd,
                             CONVERT_TO_BYTEPTR(output_p),
                             out_stride,
                             tmpbuf,
                             bit_depth,
                             1);
                    eb_apply_selfguided_restoration_c(
                        CONVERT_TO_BYTEPTR(input_p),
                        w,
                        h,
                        stride,
                        eps,
                        xqd,
                        CONVERT_TO_BYTEPTR(output2_p),
                        out_stride,
                        tmpbuf,
                        bit_depth,
                        1);
                }

            for (j = 0; j < test_h; ++j)
                for (k = 0; k < test_w; ++k)
                    ASSERT_EQ(output[j * out_stride + k],
                              output2[j * out_stride + k]);
        }

        eb_aom_free(input_);
        eb_aom_free(output_);
        eb_aom_free(output2_);
        eb_aom_free(tmpbuf);
    }

  private:
    SgrFunc tst_fun_;
};

TEST_P(AV1HighbdSelfguidedFilterTest, DISABLED_SpeedTest) {
    RunSpeedTest();
}
TEST_P(AV1HighbdSelfguidedFilterTest, CorrectnessTest) {
    RunCorrectnessTest();
}

const int32_t highbd_params_avx2[] = {8, 10, 12};
INSTANTIATE_TEST_CASE_P(
    AVX2, AV1HighbdSelfguidedFilterTest,
    ::testing::Combine(::testing::Values(eb_apply_selfguided_restoration_avx2),
                       ::testing::ValuesIn(highbd_params_avx2)));

#if 0
// To test integral_images() and integral_images_highbd(), make them not static,
// and add declarations to header file.
static void init_data_integral_images(uint8_t **src8, uint16_t **src16,
                                      int32_t *src_stride) {
    *src_stride =
        eb_create_random_aligned_stride(2 * RESTORATION_UNITSIZE_MAX, 64);
    *src8 = (uint8_t *)malloc(sizeof(**src8) * 2 * RESTORATION_UNITSIZE_MAX *
                              *src_stride);
    *src16 = (uint16_t *)malloc(sizeof(**src16) * 2 * RESTORATION_UNITSIZE_MAX *
                                *src_stride);
    eb_buf_random_u8(*src8, 2 * RESTORATION_UNITSIZE_MAX * *src_stride);
    eb_buf_random_u16_with_bd(
        *src16, 2 * RESTORATION_UNITSIZE_MAX * *src_stride, 12);
}

static void uninit_data_integral_images(uint8_t *src8, uint16_t *src16) {
    free(src8);
    free(src16);
}

static void integral_images_c(const uint8_t *src, int32_t src_stride,
                              int32_t width, int32_t height, int32_t *C,
                              int32_t *D, int32_t buf_stride) {
    for (int x = 0; x < width; x++) {
        C[x + 1] = 0;
        D[x + 1] = 0;
    }

    for (int y = 0; y < height; y++) {
        C[(y + 1) * buf_stride] = 0;
        D[(y + 1) * buf_stride] = 0;
        for (int x = 0; x < width; x++) {
            C[(y + 1) * buf_stride + x + 1] =
                C[(y + 1) * buf_stride + x] +
                src[y * src_stride + x] * src[y * src_stride + x];
            D[(y + 1) * buf_stride + x + 1] =
                D[(y + 1) * buf_stride + x] + src[y * src_stride + x];
        }
    }

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            C[(y + 1) * buf_stride + x + 1] += C[y * buf_stride + x + 1];
            D[(y + 1) * buf_stride + x + 1] += D[y * buf_stride + x + 1];
        }
    }
}

static void integral_images_highbd_c(const uint16_t *src, int32_t src_stride,
                                     int32_t width, int32_t height, int32_t *C,
                                     int32_t *D, int32_t buf_stride) {
    for (int x = 0; x < width; x++) {
        C[x + 1] = 0;
        D[x + 1] = 0;
    }

    for (int y = 0; y < height; y++) {
        C[(y + 1) * buf_stride] = 0;
        D[(y + 1) * buf_stride] = 0;
        for (int x = 0; x < width; x++) {
            C[(y + 1) * buf_stride + x + 1] =
                C[(y + 1) * buf_stride + x] +
                src[y * src_stride + x] * src[y * src_stride + x];
            D[(y + 1) * buf_stride + x + 1] =
                D[(y + 1) * buf_stride + x] + src[y * src_stride + x];
        }
    }

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            C[(y + 1) * buf_stride + x + 1] += C[y * buf_stride + x + 1];
            D[(y + 1) * buf_stride + x + 1] += D[y * buf_stride + x + 1];
        }
    }
}

TEST(IntegralImagesTest, integral_images) {
    uint8_t *src8;
    uint16_t *src16;
    int32_t src_stride;

    // The ALIGN_POWER_OF_TWO macro here ensures that column 1 of atl, btl,
    // ctl and dtl is 32-byte aligned.
    const int32_t buf_elts = ALIGN_POWER_OF_TWO(RESTORATION_PROC_UNIT_PELS, 3);
    int32_t *buf_c, *buf_o;

    buf_c = (int32_t *)eb_aom_memalign(32, sizeof(*buf_c) * 4 * buf_elts);
    buf_o = (int32_t *)eb_aom_memalign(32, sizeof(*buf_o) * 4 * buf_elts);

    for (int i = 0; i < 10; i++) {
        init_data_integral_images(&src8, &src16, &src_stride);

        for (int height = 2; height < 32; height += 2) {
            for (int width = 2; width < 32; width += 2) {
                // From av1_selfguided_restoration_avx2()
                const int32_t width_ext = width + 2 * SGRPROJ_BORDER_HORZ;
                const int32_t height_ext = height + 2 * SGRPROJ_BORDER_VERT;

                // Adjusting the stride of A and b here appears to avoid bad
                // cache effects, leading to a significant speed
                // improvement. We also align the stride to a multiple of 32
                // bytes for efficiency.
                int32_t buf_stride = ALIGN_POWER_OF_TWO(width_ext + 16, 3);

                // The "tl" pointers point at the top-left of the
                // initialised data for the array.
                int32_t *Atl_c = buf_c + 0 * buf_elts + 7;
                int32_t *Btl_c = buf_c + 1 * buf_elts + 7;
                int32_t *Ctl_c = buf_c + 2 * buf_elts + 7;
                int32_t *Dtl_c = buf_c + 3 * buf_elts + 7;
                int32_t *Atl_o = buf_o + 0 * buf_elts + 7;
                int32_t *Btl_o = buf_o + 1 * buf_elts + 7;
                int32_t *Ctl_o = buf_o + 2 * buf_elts + 7;
                int32_t *Dtl_o = buf_o + 3 * buf_elts + 7;

                // The "0" pointers are (- SGRPROJ_BORDER_VERT,
                // -SGRPROJ_BORDER_HORZ). Note there's a zero row and column
                // in A, b (integral images), so we move down and right one
                // for them.
                const int32_t buf_diag_border =
                    SGRPROJ_BORDER_HORZ + buf_stride * SGRPROJ_BORDER_VERT;

                int32_t *A0_c = Atl_c + 1 + buf_stride;
                int32_t *B0_c = Btl_c + 1 + buf_stride;
                int32_t *C0_c = Ctl_c + 1 + buf_stride;
                int32_t *D0_c = Dtl_c + 1 + buf_stride;
                int32_t *A0_o = Atl_o + 1 + buf_stride;
                int32_t *B0_o = Btl_o + 1 + buf_stride;
                int32_t *C0_o = Ctl_o + 1 + buf_stride;
                int32_t *D0_o = Dtl_o + 1 + buf_stride;

                // Finally, A, b, C, D point at position (0, 0).
                int32_t *A_c = A0_c + buf_diag_border;
                int32_t *B_c = B0_c + buf_diag_border;
                int32_t *C_c = C0_c + buf_diag_border;
                int32_t *D_c = D0_c + buf_diag_border;
                int32_t *A_o = A0_o + buf_diag_border;
                int32_t *B_o = B0_o + buf_diag_border;
                int32_t *C_o = C0_o + buf_diag_border;
                int32_t *D_o = D0_o + buf_diag_border;

                (void)A_c;
                (void)A_o;
                (void)B_c;
                (void)B_o;
                (void)C_c;
                (void)C_o;
                (void)D_c;
                (void)D_o;

                const int32_t dgd_diag_border =
                    SGRPROJ_BORDER_HORZ + src_stride * SGRPROJ_BORDER_VERT;
                // Done from av1_selfguided_restoration_avx2()

                integral_images_c(src8,
                                  src_stride,
                                  width_ext,
                                  height_ext,
                                  Ctl_c,
                                  Dtl_c,
                                  buf_stride);

                integral_images(src8,
                                src_stride,
                                width_ext,
                                height_ext,
                                Ctl_o,
                                Dtl_o,
                                buf_stride);

                for (int y = 0; y < height_ext; y++) {
                    for (int x = 0; x < width_ext; x++) {
                        const int idx = y * buf_stride + x;
                        EXPECT_EQ(C0_c[idx], C0_o[idx]);
                        EXPECT_EQ(D0_c[idx], D0_o[idx]);
                    }
                }

                integral_images_highbd_c(src16,
                                         src_stride,
                                         width_ext,
                                         height_ext,
                                         Ctl_c,
                                         Dtl_c,
                                         buf_stride);

                integral_images_highbd(src16,
                                       src_stride,
                                       width_ext,
                                       height_ext,
                                       Ctl_o,
                                       Dtl_o,
                                       buf_stride);

                for (int y = 0; y < height_ext; y++) {
                    for (int x = 0; x < width_ext; x++) {
                        const int idx = y * buf_stride + x;
                        EXPECT_EQ(C0_c[idx], C0_o[idx]);
                        EXPECT_EQ(D0_c[idx], D0_o[idx]);
                    }
                }
            }
        }

        uninit_data_integral_images(src8, src16);
    }

    eb_aom_free(buf_c);
    eb_aom_free(buf_o);
}

TEST(IntegralImagesTest, DISABLED_integral_images_speed) {
    uint8_t *src8;
    uint16_t *src16;
    int32_t src_stride;
    double time_c, time_o;
    uint64_t start_time_seconds, start_time_useconds;
    uint64_t middle_time_seconds, middle_time_useconds;
    uint64_t finish_time_seconds, finish_time_useconds;

    // The ALIGN_POWER_OF_TWO macro here ensures that column 1 of atl, btl,
    // ctl and dtl is 32-byte aligned.
    const int32_t buf_elts = ALIGN_POWER_OF_TWO(RESTORATION_PROC_UNIT_PELS, 3);
    int32_t *buf_c, *buf_o;

    buf_c = (int32_t *)eb_aom_memalign(32, sizeof(*buf_c) * 4 * buf_elts);
    buf_o = (int32_t *)eb_aom_memalign(32, sizeof(*buf_o) * 4 * buf_elts);

    init_data_integral_images(&src8, &src16, &src_stride);
    const int32_t width = 32;
    const int32_t height = 32;

    // From av1_selfguided_restoration_avx2()
    const int32_t width_ext = width + 2 * SGRPROJ_BORDER_HORZ;
    const int32_t height_ext = height + 2 * SGRPROJ_BORDER_VERT;

    // Adjusting the stride of A and b here appears to avoid bad
    // cache effects, leading to a significant speed
    // improvement. We also align the stride to a multiple of 32
    // bytes for efficiency.
    int32_t buf_stride = ALIGN_POWER_OF_TWO(width_ext + 16, 3);

    // The "tl" pointers point at the top-left of the
    // initialised data for the array.
    int32_t *Atl_c = buf_c + 0 * buf_elts + 7;
    int32_t *Btl_c = buf_c + 1 * buf_elts + 7;
    int32_t *Ctl_c = buf_c + 2 * buf_elts + 7;
    int32_t *Dtl_c = buf_c + 3 * buf_elts + 7;
    int32_t *Atl_o = buf_o + 0 * buf_elts + 7;
    int32_t *Btl_o = buf_o + 1 * buf_elts + 7;
    int32_t *Ctl_o = buf_o + 2 * buf_elts + 7;
    int32_t *Dtl_o = buf_o + 3 * buf_elts + 7;

    // The "0" pointers are (- SGRPROJ_BORDER_VERT,
    // -SGRPROJ_BORDER_HORZ). Note there's a zero row and column
    // in A, b (integral images), so we move down and right one
    // for them.
    const int32_t buf_diag_border =
        SGRPROJ_BORDER_HORZ + buf_stride * SGRPROJ_BORDER_VERT;

    int32_t *A0_c = Atl_c + 1 + buf_stride;
    int32_t *B0_c = Btl_c + 1 + buf_stride;
    int32_t *C0_c = Ctl_c + 1 + buf_stride;
    int32_t *D0_c = Dtl_c + 1 + buf_stride;
    int32_t *A0_o = Atl_o + 1 + buf_stride;
    int32_t *B0_o = Btl_o + 1 + buf_stride;
    int32_t *C0_o = Ctl_o + 1 + buf_stride;
    int32_t *D0_o = Dtl_o + 1 + buf_stride;

    // Finally, A, b, C, D point at position (0, 0).
    int32_t *A_c = A0_c + buf_diag_border;
    int32_t *B_c = B0_c + buf_diag_border;
    int32_t *C_c = C0_c + buf_diag_border;
    int32_t *D_c = D0_c + buf_diag_border;
    int32_t *A_o = A0_o + buf_diag_border;
    int32_t *B_o = B0_o + buf_diag_border;
    int32_t *C_o = C0_o + buf_diag_border;
    int32_t *D_o = D0_o + buf_diag_border;

    (void)A_c;
    (void)A_o;
    (void)B_c;
    (void)B_o;
    (void)C_c;
    (void)C_o;
    (void)D_c;
    (void)D_o;

    const int32_t dgd_diag_border =
        SGRPROJ_BORDER_HORZ + src_stride * SGRPROJ_BORDER_VERT;
    // Done from av1_selfguided_restoration_avx2()

    const uint64_t num_loop = 1000000;

    eb_start_time(&start_time_seconds, &start_time_useconds);

    for (uint64_t i = 0; i < num_loop; i++) {
        integral_images_c(
            src8, src_stride, width_ext, height_ext, Ctl_c, Dtl_c, buf_stride);
    }

    eb_start_time(&middle_time_seconds, &middle_time_useconds);

    for (uint64_t i = 0; i < num_loop; i++) {
        integral_images(
            src8, src_stride, width_ext, height_ext, Ctl_o, Dtl_o, buf_stride);
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

    for (int y = 0; y < height_ext; y++) {
        for (int x = 0; x < width_ext; x++) {
            const int idx = y * buf_stride + x;
            EXPECT_EQ(C0_c[idx], C0_o[idx]);
            EXPECT_EQ(D0_c[idx], D0_o[idx]);
        }
    }

    printf("Average Nanoseconds per Function Call\n");
    printf("    integral_images_c    : %6.2f\n", 1000000 * time_c / num_loop);
    printf("    integral_images_avx2 : %6.2f   (Comparison: %5.2fx)\n",
           1000000 * time_o / num_loop,
           time_c / time_o);

    eb_start_time(&start_time_seconds, &start_time_useconds);

    for (uint64_t i = 0; i < num_loop; i++) {
        integral_images_highbd_c(
            src16, src_stride, width_ext, height_ext, Ctl_c, Dtl_c, buf_stride);
    }

    eb_start_time(&middle_time_seconds, &middle_time_useconds);

    for (uint64_t i = 0; i < num_loop; i++) {
        integral_images_highbd(
            src16, src_stride, width_ext, height_ext, Ctl_o, Dtl_o, buf_stride);
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

    for (int y = 0; y < height_ext; y++) {
        for (int x = 0; x < width_ext; x++) {
            const int idx = y * buf_stride + x;
            EXPECT_EQ(C0_c[idx], C0_o[idx]);
            EXPECT_EQ(D0_c[idx], D0_o[idx]);
        }
    }

    printf("Average Nanoseconds per Function Call\n");
    printf("    integral_images_highbd_c    : %6.2f\n",
           1000000 * time_c / num_loop);
    printf("    integral_images_highbd_avx2 : %6.2f   (Comparison: %5.2fx)\n",
           1000000 * time_o / num_loop,
           time_c / time_o);

    uninit_data_integral_images(src8, src16);
    eb_aom_free(buf_c);
    eb_aom_free(buf_o);
}
#endif

}  // namespace
