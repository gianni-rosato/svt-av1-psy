/*
 * Copyright(c) 2020 Intel Corporation
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

#include "aom_dsp_rtcd.h"
#include "gtest/gtest.h"
#include "EbTemporalFiltering.h"
#include "util.h"
#include "random.h"
#include "EbUtility.h"
#include "EbUnitTestUtility.h"


using svt_av1_test_tool::SVTRandom;

typedef void (*TemporalFilterFunc)(
    const uint8_t *y_src, int y_src_stride, const uint8_t *y_pre,
    int y_pre_stride, const uint8_t *u_src, const uint8_t *v_src,
    int uv_src_stride, const uint8_t *u_pre, const uint8_t *v_pre,
    int uv_pre_stride, unsigned int block_width, unsigned int block_height,
    int ss_x, int ss_y, const double *noise_levels, const int decay_control,
    uint32_t *y_accum, uint16_t *y_count, uint32_t *u_accum, uint16_t *u_count,
    uint32_t *v_accum, uint16_t *v_count);

#define MAX_STRIDE 256

typedef std::tuple<TemporalFilterFunc, TemporalFilterFunc>
    TemporalFilterWithParam;

class TemporalFilterTestPlanewise
    : public ::testing::TestWithParam<TemporalFilterWithParam> {
  public:
    TemporalFilterTestPlanewise() : rnd_(0, 255){};
    ~TemporalFilterTestPlanewise() {
    }

    void SetUp() {
        ref_func = TEST_GET_PARAM(0);
        tst_func = TEST_GET_PARAM(1);

        for (int color_channel = 0; color_channel < COLOR_CHANNELS; color_channel++) {
            src_ptr[color_channel] = reinterpret_cast<uint8_t *>(
                eb_aom_memalign(8, MAX_STRIDE * MAX_STRIDE));
            pred_ptr[color_channel] = reinterpret_cast<uint8_t *>(
                eb_aom_memalign(8, MAX_STRIDE * MAX_STRIDE));

            accum_ref_ptr[color_channel] =
                reinterpret_cast<uint32_t *>(
                eb_aom_memalign(8, MAX_STRIDE * MAX_STRIDE * sizeof(uint32_t)));
            count_ref_ptr[color_channel] =
                reinterpret_cast<uint16_t *>(
                eb_aom_memalign(8, MAX_STRIDE * MAX_STRIDE * sizeof(uint16_t)));
            accum_tst_ptr[color_channel] =
                reinterpret_cast<uint32_t *>(
                eb_aom_memalign(8, MAX_STRIDE * MAX_STRIDE * sizeof(uint32_t)));
            count_tst_ptr[color_channel] =
                reinterpret_cast<uint16_t *>(
                eb_aom_memalign(8, MAX_STRIDE * MAX_STRIDE * sizeof(uint16_t)));

            memset(accum_ref_ptr[color_channel],
                   0,
                   MAX_STRIDE * MAX_STRIDE * sizeof(uint32_t));
            memset(count_ref_ptr[color_channel],
                   0,
                   MAX_STRIDE * MAX_STRIDE * sizeof(uint16_t));
            memset(accum_tst_ptr[color_channel],
                   0,
                   MAX_STRIDE * MAX_STRIDE * sizeof(uint32_t));
            memset(count_tst_ptr[color_channel],
                   0,
                   MAX_STRIDE * MAX_STRIDE * sizeof(uint16_t));

            stride[color_channel] = MAX_STRIDE;
            stride_pred[color_channel] = MAX_STRIDE;

            noise_levels[color_channel] = 2.1002103677063437;
        }

        decay_control = 7;
    }

    void TearDown() {

        for (int color_channel = 0; color_channel < COLOR_CHANNELS; color_channel++) {
            eb_aom_free(src_ptr[color_channel]);
            eb_aom_free(pred_ptr[color_channel]);

            eb_aom_free(accum_ref_ptr[color_channel]);
            eb_aom_free(count_ref_ptr[color_channel]);
            eb_aom_free(accum_tst_ptr[color_channel]);
            eb_aom_free(count_tst_ptr[color_channel]);
        }

    }
    void RunTest(int width, int height, int run_times);

    void GenRandomData(int width, int height, int stride, int stride2) {
        for (int ii = 0; ii < height; ii++) {
            for (int jj = 0; jj < width; jj++) {
                src_ptr[C_Y][ii * stride + jj] = rnd_.Rand8();
                src_ptr[C_U][ii * stride + jj] = rnd_.Rand8();
                src_ptr[C_V][ii * stride + jj] = rnd_.Rand8();

                pred_ptr[C_Y][ii * stride2 + jj] = rnd_.Rand8();
                pred_ptr[C_U][ii * stride2 + jj] = rnd_.Rand8();
                pred_ptr[C_V][ii * stride2 + jj] = rnd_.Rand8();
            }
        }
    }

  private:
    TemporalFilterFunc ref_func;
    TemporalFilterFunc tst_func;
    SVTRandom rnd_;
    uint8_t *src_ptr[COLOR_CHANNELS];
    uint8_t *pred_ptr[COLOR_CHANNELS];
    uint32_t *accum_ref_ptr[COLOR_CHANNELS];
    uint16_t *count_ref_ptr[COLOR_CHANNELS];

    uint32_t *accum_tst_ptr[COLOR_CHANNELS];
    uint16_t *count_tst_ptr[COLOR_CHANNELS];

    uint32_t stride[COLOR_CHANNELS];
    uint32_t stride_pred[COLOR_CHANNELS];
    double noise_levels[COLOR_CHANNELS];
    int decay_control;
};

void TemporalFilterTestPlanewise::RunTest(int width, int height,
                                          int run_times) {
    if (run_times <= 100) {
        for (int j = 0; j < run_times; j++) {
            GenRandomData(width, height, MAX_STRIDE, MAX_STRIDE);
            ref_func(src_ptr[C_Y],
                     stride[C_Y],
                     pred_ptr[C_Y],
                     stride_pred[C_Y],
                     src_ptr[C_U],
                     src_ptr[C_V],
                     stride[C_U],
                     pred_ptr[C_U],
                     pred_ptr[C_V],
                     stride_pred[C_U],
                     width,
                     height,
                     1,  // subsampling
                     1,  // subsampling
                     noise_levels,
                     decay_control,
                     accum_ref_ptr[C_Y],
                     count_ref_ptr[C_Y],
                     accum_ref_ptr[C_U],
                     count_ref_ptr[C_U],
                     accum_ref_ptr[C_V],
                     count_ref_ptr[C_V]);

            tst_func(src_ptr[C_Y],
                     stride[C_Y],
                     pred_ptr[C_Y],
                     stride_pred[C_Y],
                     src_ptr[C_U],
                     src_ptr[C_V],
                     stride[C_U],
                     pred_ptr[C_U],
                     pred_ptr[C_V],
                     stride_pred[C_U],
                     width,
                     height,
                     1,  // subsampling
                     1,  // subsampling
                     noise_levels,
                     decay_control,
                     accum_tst_ptr[C_Y],
                     count_tst_ptr[C_Y],
                     accum_tst_ptr[C_U],
                     count_tst_ptr[C_U],
                     accum_tst_ptr[C_V],
                     count_tst_ptr[C_V]);

            for (int color_channel = 0; color_channel < COLOR_CHANNELS;
                 color_channel++) {
                EXPECT_EQ(memcmp(accum_ref_ptr[color_channel],
                                 accum_tst_ptr[color_channel],
                                 MAX_STRIDE * MAX_STRIDE * sizeof(uint32_t)),
                          0);

                EXPECT_EQ(memcmp(count_ref_ptr[color_channel],
                                 count_tst_ptr[color_channel],
                                 MAX_STRIDE * MAX_STRIDE * sizeof(uint16_t)),
                          0);
            }
        }
    }
    else{
        uint64_t ref_timer_seconds, ref_timer_useconds;
        uint64_t middle_timer_seconds, middle_timer_useconds;
        uint64_t test_timer_seconds, test_timer_useconds;
        double ref_time, tst_time;

        eb_start_time(&ref_timer_seconds, &ref_timer_useconds);
        for (int j = 0; j < run_times; j++) {
            ref_func(src_ptr[C_Y],
                     stride[C_Y],
                     pred_ptr[C_Y],
                     stride_pred[C_Y],
                     src_ptr[C_U],
                     src_ptr[C_V],
                     stride[C_U],
                     pred_ptr[C_U],
                     pred_ptr[C_V],
                     stride_pred[C_U],
                     width,
                     height,
                     1,  // subsampling
                     1,  // subsampling
                     noise_levels,
                     decay_control,
                     accum_ref_ptr[C_Y],
                     count_ref_ptr[C_Y],
                     accum_ref_ptr[C_U],
                     count_ref_ptr[C_U],
                     accum_ref_ptr[C_V],
                     count_ref_ptr[C_V]);
        }
        eb_start_time(&middle_timer_seconds, &middle_timer_useconds);

        for (int j = 0; j < run_times; j++) {
            tst_func(src_ptr[C_Y],
                     stride[C_Y],
                     pred_ptr[C_Y],
                     stride_pred[C_Y],
                     src_ptr[C_U],
                     src_ptr[C_V],
                     stride[C_U],
                     pred_ptr[C_U],
                     pred_ptr[C_V],
                     stride_pred[C_U],
                     width,
                     height,
                     1,  // subsampling
                     1,  // subsampling
                     noise_levels,
                     decay_control,
                     accum_tst_ptr[C_Y],
                     count_tst_ptr[C_Y],
                     accum_tst_ptr[C_U],
                     count_tst_ptr[C_U],
                     accum_tst_ptr[C_V],
                     count_tst_ptr[C_V]);
        }
        eb_start_time(&test_timer_seconds, &test_timer_useconds);

        eb_compute_overall_elapsed_time_ms(ref_timer_seconds,
                                           ref_timer_useconds,
                                           middle_timer_seconds,
                                           middle_timer_useconds,
                                           &ref_time);

        eb_compute_overall_elapsed_time_ms(middle_timer_seconds,
                                           middle_timer_useconds,
                                           test_timer_seconds,
                                           test_timer_useconds,
                                           &tst_time);

        printf(
            "c_time=%lf \t simd_time=%lf \t "
            "gain=%lf\t width=%d\t height=%d \n",
            ref_time / 1000,
            tst_time / 1000,
            ref_time / tst_time,
            width,
            height);
    }
}

TEST_P(TemporalFilterTestPlanewise, OperationCheck) {
    for (int height = 32; height <= 32; height = height * 2) {
        RunTest(height, height, 100);
    }
}

TEST_P(TemporalFilterTestPlanewise, DISABLED_Speed) {
    for (int height = 32; height <= 32; height = height * 2) {
        RunTest(height, height, 100000);
    }
}


INSTANTIATE_TEST_CASE_P(
    AVX2, TemporalFilterTestPlanewise,
    ::testing::Combine(::testing::Values(svt_av1_apply_temporal_filter_planewise_c),
                       ::testing::Values(svt_av1_apply_temporal_filter_planewise_avx2)));

