/*
 * Copyright (c) 2019, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "gtest/gtest.h"
#include "aom_dsp_rtcd.h"
#include "EbUnitTestUtility.h"
#include "EbUtility.h"
#include "random.h"
#include "util.h"

namespace {
typedef int64_t (*frame_error_func)(const uint8_t *const ref, int stride,
                                    const uint8_t *const dst, int p_width,
                                    int p_height, int p_stride);
const int kBlockWidth[] = {
    832,
    834,
    640,
    1280,
    1920,
};
const int kBlockHeight[] = {
    480,
    482,
    360,
    720,
    1080,
};
typedef std::tuple<frame_error_func, int, int> FrameErrorParam;

class AV1FrameErrorTest : public ::testing::TestWithParam<FrameErrorParam> {
  public:
    virtual ~AV1FrameErrorTest() {
    }
    virtual void SetUp() {
        rnd_ = new svt_av1_test_tool::SVTRandom(0, (1 << 8) - 1);
    }
    virtual void TearDown() {
        delete rnd_;
        aom_clear_system_state();
    }

  protected:
    void RandomValues(frame_error_func test_impl, int width, int height);
    void ExtremeValues(frame_error_func test_impl, int width, int height);
    void RunSpeedTest(frame_error_func test_impl, int width, int height);
    svt_av1_test_tool::SVTRandom *rnd_;
};

void AV1FrameErrorTest::RandomValues(frame_error_func test_impl, int width,
                                     int height) {
    const int stride = (((width * 3) / 2) + 15) & ~15;
    const int max_blk_size = stride * height;
    uint8_t *const dst =
        static_cast<uint8_t *>(eb_aom_memalign(16, max_blk_size * sizeof(*dst)));
    uint8_t *const ref =
        static_cast<uint8_t *>(eb_aom_memalign(16, max_blk_size * sizeof(*ref)));
    ASSERT_TRUE(dst != NULL);
    ASSERT_TRUE(ref != NULL);
    for (int i = 0; i < max_blk_size; ++i) {
        dst[i] = rnd_->Rand8();
        ref[i] = rnd_->Rand8();
    }
    const int64_t ref_error =
        eb_av1_calc_frame_error_c(ref, stride, dst, width, height, stride);
    const int64_t test_error =
        test_impl(ref, stride, dst, width, height, stride);
    ASSERT_EQ(test_error, ref_error) << width << "x" << height;
    eb_aom_free(dst);
    eb_aom_free(ref);
}

void AV1FrameErrorTest::ExtremeValues(frame_error_func test_impl, int width,
                                      int height) {
    const int stride = (((width * 3) / 2) + 15) & ~15;
    const int max_blk_size = stride * height;
    uint8_t *const dst =
        static_cast<uint8_t *>(eb_aom_memalign(16, max_blk_size * sizeof(*dst)));
    uint8_t *const ref =
        static_cast<uint8_t *>(eb_aom_memalign(16, max_blk_size * sizeof(*ref)));
    ASSERT_TRUE(dst != NULL);
    ASSERT_TRUE(ref != NULL);
    for (int r = 0; r < 2; r++) {
        if (r == 0) {
            memset(dst, 0, max_blk_size);
            memset(ref, 255, max_blk_size);
        } else if (r == 1) {
            memset(dst, 255, max_blk_size);
            memset(ref, 0, max_blk_size);
        }
        const int64_t ref_error =
            eb_av1_calc_frame_error_c(ref, stride, dst, width, height, stride);
        const int64_t test_error =
            test_impl(ref, stride, dst, width, height, stride);
        ASSERT_EQ(test_error, ref_error) << width << "x" << height;
    }
    eb_aom_free(dst);
    eb_aom_free(ref);
}

void AV1FrameErrorTest::RunSpeedTest(frame_error_func test_impl, int width,
                                     int height) {
    const int stride = (((width * 3) / 2) + 15) & ~15;
    const int max_blk_size = stride * height;
    uint8_t *const dst =
        static_cast<uint8_t *>(eb_aom_memalign(16, max_blk_size * sizeof(*dst)));
    uint8_t *const ref =
        static_cast<uint8_t *>(eb_aom_memalign(16, max_blk_size * sizeof(*ref)));
    ASSERT_TRUE(dst != NULL);
    ASSERT_TRUE(ref != NULL);
    for (int i = 0; i < max_blk_size; ++i) {
        dst[i] = ref[i] = rnd_->Rand8();
    }
    const int num_loops = 10000000 / (width + height);
    frame_error_func funcs[2] = {eb_av1_calc_frame_error_c, test_impl};
    double elapsed_time[2] = {0};
    for (int i = 0; i < 2; ++i) {
        double time;
        uint64_t start_time_seconds, start_time_useconds;
        uint64_t finish_time_seconds, finish_time_useconds;

        svt_av1_get_time(&start_time_seconds, &start_time_useconds);
        frame_error_func func = funcs[i];
        for (int j = 0; j < num_loops; ++j) {
            func(ref, stride, dst, width, height, stride);
        }
        svt_av1_get_time(&finish_time_seconds, &finish_time_useconds);
        time = svt_av1_compute_overall_elapsed_time_ms(start_time_seconds,
                                                       start_time_useconds,
                                                       finish_time_seconds,
                                                       finish_time_useconds);
        elapsed_time[i] = 1000.0 * time / num_loops;
    }
    eb_aom_free(dst);
    eb_aom_free(ref);
    printf("av1_calc_frame_error %3dx%-3d: %7.2f/%7.2fns",
           width,
           height,
           elapsed_time[0],
           elapsed_time[1]);
    printf("(%3.2f)\n", elapsed_time[0] / elapsed_time[1]);
}

TEST_P(AV1FrameErrorTest, CheckOutput) {
    RandomValues(TEST_GET_PARAM(0), TEST_GET_PARAM(1), TEST_GET_PARAM(2));
    ExtremeValues(TEST_GET_PARAM(0), TEST_GET_PARAM(1), TEST_GET_PARAM(2));
}

TEST_P(AV1FrameErrorTest, DISABLED_Speed) {
    RunSpeedTest(TEST_GET_PARAM(0), TEST_GET_PARAM(1), TEST_GET_PARAM(2));
}

INSTANTIATE_TEST_CASE_P(
    AVX2, AV1FrameErrorTest,
    ::testing::Combine(::testing::Values(&eb_av1_calc_frame_error_avx2),
                       ::testing::ValuesIn(kBlockWidth),
                       ::testing::ValuesIn(kBlockHeight)));

}  // namespace
