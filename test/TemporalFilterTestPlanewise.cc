/*
* Copyright(c) 2020 Intel Corporation
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
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
    struct MeContext *context_ptr, const uint8_t *y_src, int y_src_stride,
    const uint8_t *y_pre, int y_pre_stride, const uint8_t *u_src,
    const uint8_t *v_src, int uv_src_stride, const uint8_t *u_pre,
    const uint8_t *v_pre, int uv_pre_stride, unsigned int block_width,
    unsigned int block_height, int ss_x, int ss_y, const double *noise_levels,
    const int decay_control, uint32_t *y_accum, uint16_t *y_count,
    uint32_t *u_accum, uint16_t *u_count, uint32_t *v_accum, uint16_t *v_count);

typedef void (*TemporalFilterFuncHbd)(
    struct MeContext *context_ptr, const uint16_t *y_src, int y_src_stride,
    const uint16_t *y_pre, int y_pre_stride, const uint16_t *u_src,
    const uint16_t *v_src, int uv_src_stride, const uint16_t *u_pre,
    const uint16_t *v_pre, int uv_pre_stride, unsigned int block_width,
    unsigned int block_height, int ss_x, int ss_y, const double *noise_levels,
    const int decay_control, uint32_t *y_accum, uint16_t *y_count,
    uint32_t *u_accum, uint16_t *u_count, uint32_t *v_accum, uint16_t *v_count, uint32_t encoder_bit_depth);

#define MAX_STRIDE 256

typedef std::tuple<TemporalFilterFunc, TemporalFilterFunc>
    TemporalFilterWithParam;

typedef std::tuple<TemporalFilterFuncHbd, TemporalFilterFuncHbd>
    TemporalFilterWithParamHbd;


static void TemporalFilterFillMeContexts(MeContext *cnt1, MeContext *cnt2) {
    //Prepare two MeContext to cover all paths in code.
    cnt1->tf_32x32_block_error[0] = 0;
    cnt1->tf_32x32_block_error[1] = 156756;
    cnt1->tf_32x32_block_error[2] = 2;
    cnt1->tf_32x32_block_error[3] = 0;
    cnt1->tf_block_col = 0;  //or 1
    cnt1->tf_block_row = 0;
    cnt1->tf_32x32_block_split_flag[0] = 1;
    cnt1->tf_32x32_block_split_flag[1] = 0;
    cnt1->tf_32x32_block_split_flag[2] = 1;
    cnt1->tf_32x32_block_split_flag[3] = 0;
    cnt1->tf_16x16_block_error[0] = 29087;
    cnt1->tf_16x16_block_error[1] = 28482;
    cnt1->tf_16x16_block_error[2] = 44508;
    cnt1->tf_16x16_block_error[3] = 42356;
    cnt1->tf_16x16_block_error[4] = 35319;
    cnt1->tf_16x16_block_error[5] = 37670;
    cnt1->tf_16x16_block_error[6] = 38551;
    cnt1->tf_16x16_block_error[7] = 39320;
    cnt1->tf_16x16_block_error[8] = 28360;
    cnt1->tf_16x16_block_error[9] = 32753;
    cnt1->tf_16x16_block_error[10] = 38335;
    cnt1->tf_16x16_block_error[11] = 50679;
    cnt1->tf_16x16_block_error[12] = 44700;
    cnt1->tf_16x16_block_error[13] = 49620;
    cnt1->tf_16x16_block_error[14] = 47070;
    cnt1->tf_16x16_block_error[15] = 41171;
    cnt1->tf_32x32_block_error[0] = 156756;
    cnt1->tf_32x32_block_error[1] = 168763;
    cnt1->tf_32x32_block_error[2] = 143351;
    cnt1->tf_32x32_block_error[3] = 189005;
    cnt1->tf_16x16_mv_x[0] = 16;
    cnt1->tf_16x16_mv_x[1] = 13;
    cnt1->tf_16x16_mv_x[2] = 13;
    cnt1->tf_16x16_mv_x[3] = 13;
    cnt1->tf_16x16_mv_x[4] = 44;
    cnt1->tf_16x16_mv_x[5] = 12;
    cnt1->tf_16x16_mv_x[6] = 1;
    cnt1->tf_16x16_mv_x[7] = 13;
    cnt1->tf_16x16_mv_x[8] = 12;
    cnt1->tf_16x16_mv_x[9] = 52;
    cnt1->tf_16x16_mv_x[10] = 11;
    cnt1->tf_16x16_mv_x[11] = 52;
    cnt1->tf_16x16_mv_x[12] = -20;
    cnt1->tf_16x16_mv_x[13] = 13;
    cnt1->tf_16x16_mv_x[14] = 11;
    cnt1->tf_16x16_mv_x[15] = -4;
    cnt1->tf_16x16_mv_y[0] = -39;
    cnt1->tf_16x16_mv_y[1] = -35;
    cnt1->tf_16x16_mv_y[2] = -38;
    cnt1->tf_16x16_mv_y[3] = -60;
    cnt1->tf_16x16_mv_y[4] = -33;
    cnt1->tf_16x16_mv_y[5] = -19;
    cnt1->tf_16x16_mv_y[6] = -35;
    cnt1->tf_16x16_mv_y[7] = -34;
    cnt1->tf_16x16_mv_y[8] = -39;
    cnt1->tf_16x16_mv_y[9] = -19;
    cnt1->tf_16x16_mv_y[10] = -39;
    cnt1->tf_16x16_mv_y[11] = 21;
    cnt1->tf_16x16_mv_y[12] = -52;
    cnt1->tf_16x16_mv_y[13] = -31;
    cnt1->tf_16x16_mv_y[14] = -12;
    cnt1->tf_16x16_mv_y[15] = -36;
    cnt1->tf_32x32_mv_x[0] = -196;
    cnt1->tf_32x32_mv_x[1] = -188;
    cnt1->tf_32x32_mv_x[2] = -228;
    cnt1->tf_32x32_mv_x[3] = -220;
    cnt1->tf_32x32_mv_y[0] = -436;
    cnt1->tf_32x32_mv_y[1] = -436;
    cnt1->tf_32x32_mv_y[2] = -436;
    cnt1->tf_32x32_mv_y[3] = -420;
    cnt1->min_frame_size = 1080;
    *cnt2 = *cnt1;
    cnt2->tf_block_col = 1;
}

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

    struct MeContext context1, context2, *context_ptr;
    TemporalFilterFillMeContexts(&context1, &context2);

    if (run_times <= 100) {
        for (int j = 0; j < run_times; j++) {
            GenRandomData(width, height, MAX_STRIDE, MAX_STRIDE);
            if (j % 2 == 0) {
                context_ptr = &context1;
            } else {
                context_ptr = &context2;
            }
            ref_func(
                     context_ptr,
                     src_ptr[C_Y],
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

            tst_func(
                     context_ptr,
                     src_ptr[C_Y],
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

        svt_av1_get_time(&ref_timer_seconds, &ref_timer_useconds);
        for (int j = 0; j < run_times; j++) {
            if (j % 2 == 0) {
                context_ptr = &context1;
            } else {
                context_ptr = &context2;
            }
            ref_func(
                     context_ptr,
                     src_ptr[C_Y],
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
        svt_av1_get_time(&middle_timer_seconds, &middle_timer_useconds);

        for (int j = 0; j < run_times; j++) {
            if (j % 2 == 0) {
                context_ptr = &context1;
            } else {
                context_ptr = &context2;
            }
            tst_func(
                     context_ptr,
                     src_ptr[C_Y],
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
        svt_av1_get_time(&test_timer_seconds, &test_timer_useconds);

        ref_time =
            svt_av1_compute_overall_elapsed_time_ms(ref_timer_seconds,
                                                    ref_timer_useconds,
                                                    middle_timer_seconds,
                                                    middle_timer_useconds);

        tst_time =
            svt_av1_compute_overall_elapsed_time_ms(middle_timer_seconds,
                                                    middle_timer_useconds,
                                                    test_timer_seconds,
                                                    test_timer_useconds);

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

// MSVC fails because avx2 kernel does not exist, temporal fix by assigning C
// kernel instread of AVX2
INSTANTIATE_TEST_CASE_P(
    AVX2, TemporalFilterTestPlanewise,
    ::testing::Combine(::testing::Values(svt_av1_apply_temporal_filter_planewise_c),
                       ::testing::Values(svt_av1_apply_temporal_filter_planewise_avx2)));


class TemporalFilterTestPlanewiseHbd
    : public ::testing::TestWithParam<TemporalFilterWithParamHbd> {
  public:
    TemporalFilterTestPlanewiseHbd() : rnd_(0, (1 << 10) - 1){};
    ~TemporalFilterTestPlanewiseHbd() {
    }

    void SetUp() {
        ref_func = TEST_GET_PARAM(0);
        tst_func = TEST_GET_PARAM(1);

        for (int color_channel = 0; color_channel < COLOR_CHANNELS;
             color_channel++) {
            src_ptr[color_channel] = reinterpret_cast<uint16_t *>(
                eb_aom_memalign(8, MAX_STRIDE * MAX_STRIDE));
            pred_ptr[color_channel] = reinterpret_cast<uint16_t *>(
                eb_aom_memalign(8, MAX_STRIDE * MAX_STRIDE));

            accum_ref_ptr[color_channel] = reinterpret_cast<uint32_t *>(
                eb_aom_memalign(8, MAX_STRIDE * MAX_STRIDE * sizeof(uint32_t)));
            count_ref_ptr[color_channel] = reinterpret_cast<uint16_t *>(
                eb_aom_memalign(8, MAX_STRIDE * MAX_STRIDE * sizeof(uint16_t)));
            accum_tst_ptr[color_channel] = reinterpret_cast<uint32_t *>(
                eb_aom_memalign(8, MAX_STRIDE * MAX_STRIDE * sizeof(uint32_t)));
            count_tst_ptr[color_channel] = reinterpret_cast<uint16_t *>(
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
        for (int color_channel = 0; color_channel < COLOR_CHANNELS;
             color_channel++) {
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
                src_ptr[C_Y][ii * stride + jj] = rnd_.random();
                src_ptr[C_U][ii * stride + jj] = rnd_.random();
                src_ptr[C_V][ii * stride + jj] = rnd_.random();

                pred_ptr[C_Y][ii * stride2 + jj] = rnd_.random();
                pred_ptr[C_U][ii * stride2 + jj] = rnd_.random();
                pred_ptr[C_V][ii * stride2 + jj] = rnd_.random();
            }
        }
    }

  private:
    TemporalFilterFuncHbd ref_func;
    TemporalFilterFuncHbd tst_func;
    SVTRandom rnd_;
    uint16_t *src_ptr[COLOR_CHANNELS];
    uint16_t *pred_ptr[COLOR_CHANNELS];
    uint32_t *accum_ref_ptr[COLOR_CHANNELS];
    uint16_t *count_ref_ptr[COLOR_CHANNELS];

    uint32_t *accum_tst_ptr[COLOR_CHANNELS];
    uint16_t *count_tst_ptr[COLOR_CHANNELS];

    uint32_t stride[COLOR_CHANNELS];
    uint32_t stride_pred[COLOR_CHANNELS];
    double noise_levels[COLOR_CHANNELS];
    int decay_control;
    uint32_t encoder_bit_depth;
};

void TemporalFilterTestPlanewiseHbd::RunTest(int width, int height,
                                             int run_times) {

    struct MeContext context1, context2, *context_ptr;
    TemporalFilterFillMeContexts(&context1, &context2);

    if (run_times <= 100) {
        for (int j = 0; j < run_times; j++) {
            GenRandomData(width, height, MAX_STRIDE, MAX_STRIDE);
            if(j%2 == 0) {
                context_ptr = &context1;
            } else {
                context_ptr = &context2;
            }
            encoder_bit_depth = 10;
            ref_func(
                     context_ptr,
                     src_ptr[C_Y],
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
                     count_ref_ptr[C_V],
                     encoder_bit_depth);

            tst_func(
                     context_ptr,
                     src_ptr[C_Y],
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
                     count_tst_ptr[C_V],
                     encoder_bit_depth);

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
            encoder_bit_depth = 12;
            ref_func(
                context_ptr,
                src_ptr[C_Y],
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
                count_ref_ptr[C_V],
                encoder_bit_depth);

            tst_func(
                context_ptr,
                src_ptr[C_Y],
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
                count_tst_ptr[C_V],
                encoder_bit_depth);

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
    } else {
        uint64_t ref_timer_seconds, ref_timer_useconds;
        uint64_t middle_timer_seconds, middle_timer_useconds;
        uint64_t test_timer_seconds, test_timer_useconds;
        double ref_time, tst_time;

        encoder_bit_depth = 10;
        svt_av1_get_time(&ref_timer_seconds, &ref_timer_useconds);
        for (int j = 0; j < run_times; j++) {
            if (j % 2 == 0) {
                context_ptr = &context1;
            } else {
                context_ptr = &context2;
            }
            ref_func(
                     context_ptr,
                     src_ptr[C_Y],
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
                     count_ref_ptr[C_V],
                     encoder_bit_depth);
        }
        svt_av1_get_time(&middle_timer_seconds, &middle_timer_useconds);

        for (int j = 0; j < run_times; j++) {
            tst_func(
                     context_ptr,
                     src_ptr[C_Y],
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
                     count_tst_ptr[C_V],
                     encoder_bit_depth);
        }
        svt_av1_get_time(&test_timer_seconds, &test_timer_useconds);

        ref_time =
            svt_av1_compute_overall_elapsed_time_ms(ref_timer_seconds,
                                                    ref_timer_useconds,
                                                    middle_timer_seconds,
                                                    middle_timer_useconds);

        tst_time =
            svt_av1_compute_overall_elapsed_time_ms(middle_timer_seconds,
                                                    middle_timer_useconds,
                                                    test_timer_seconds,
                                                    test_timer_useconds);

        printf(
            "c_time=%lf \t simd_time=%lf \t "
            "gain=%lf\t width=%d\t height=%d \t encoder_bit_depth=%d \n",
            ref_time / 1000,
            tst_time / 1000,
            ref_time / tst_time,
            width,
            height,
            encoder_bit_depth);

        encoder_bit_depth = 12;
        svt_av1_get_time(&ref_timer_seconds, &ref_timer_useconds);
        for (int j = 0; j < run_times; j++) {
            if (j % 2 == 0) {
                context_ptr = &context1;
            }
            else {
                context_ptr = &context2;
            }
            ref_func(
                context_ptr,
                src_ptr[C_Y],
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
                count_ref_ptr[C_V],
                encoder_bit_depth);
        }
        svt_av1_get_time(&middle_timer_seconds, &middle_timer_useconds);

        for (int j = 0; j < run_times; j++) {
            tst_func(
                context_ptr,
                src_ptr[C_Y],
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
                count_tst_ptr[C_V],
                encoder_bit_depth);
        }
        svt_av1_get_time(&test_timer_seconds, &test_timer_useconds);

        ref_time =
            svt_av1_compute_overall_elapsed_time_ms(ref_timer_seconds,
                                                    ref_timer_useconds,
                                                    middle_timer_seconds,
                                                    middle_timer_useconds);

        tst_time =
            svt_av1_compute_overall_elapsed_time_ms(middle_timer_seconds,
                                                    middle_timer_useconds,
                                                    test_timer_seconds,
                                                    test_timer_useconds);

        printf(
            "c_time=%lf \t simd_time=%lf \t "
            "gain=%lf\t width=%d\t height=%d \t encoder_bit_depth=%d \n",
            ref_time / 1000,
            tst_time / 1000,
            ref_time / tst_time,
            width,
            height,
            encoder_bit_depth);
    }
}

TEST_P(TemporalFilterTestPlanewiseHbd, OperationCheck) {
    for (int height = 32; height <= 32; height = height * 2) {
        RunTest(height, height, 100);
    }
}

TEST_P(TemporalFilterTestPlanewiseHbd, DISABLED_Speed) {
    for (int height = 32; height <= 32; height = height * 2) {
        RunTest(height, height, 100000);
    }
}

INSTANTIATE_TEST_CASE_P(
    AVX2, TemporalFilterTestPlanewiseHbd,
    ::testing::Combine(
        ::testing::Values(svt_av1_apply_temporal_filter_planewise_hbd_c),
        ::testing::Values(svt_av1_apply_temporal_filter_planewise_hbd_avx2)));
