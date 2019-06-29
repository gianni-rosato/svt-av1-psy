/*
 * Copyright(c) 2019 Intel Corporation
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "aom_dsp_rtcd.h"
#include "EbDefinitions.h"
#include "EbPictureOperators_AVX2.h"
#include "EbTransforms.h"
#include "EbUnitTestUtility.h"
#include "random.h"
#include "util.h"

namespace {

class SpatialFullDistortionTest : public ::testing::Test {
  public:
    ~SpatialFullDistortionTest();

    void SetUp() {
        input_stride_ = eb_create_random_aligned_stride(MAX_SB_SIZE, 64);
        recon_stride_ = eb_create_random_aligned_stride(MAX_SB_SIZE, 64);
        input_ =
            (uint8_t *)malloc(sizeof(*input_) * MAX_SB_SIZE * input_stride_);
        recon_ =
            (uint8_t *)malloc(sizeof(*recon_) * MAX_SB_SIZE * recon_stride_);
    }
    void TearDown() {
        free(recon_);
        free(input_);
        aom_clear_system_state();
    }

  protected:
    void RunCheckOutput();
    void RunSpeedTest();

    void init_data() {
        eb_buf_random_u8(input_, MAX_SB_SIZE * input_stride_);
        eb_buf_random_u8(recon_, MAX_SB_SIZE * recon_stride_);
    }

    uint8_t *input_;
    uint8_t *recon_;
    uint32_t input_stride_;
    uint32_t recon_stride_;
};

SpatialFullDistortionTest::~SpatialFullDistortionTest() {
}

void SpatialFullDistortionTest::RunCheckOutput() {
    for (int i = 0; i < 10; i++) {
        init_data();
        for (uint32_t area_width = 4; area_width <= 128; area_width += 4) {
            for (uint32_t area_height = 4; area_height <= 32;
                 area_height += 4) {
                const uint64_t dist_org =
                    spatial_full_distortion_kernel_c(input_,
                                                     input_stride_,
                                                     recon_,
                                                     recon_stride_,
                                                     area_width,
                                                     area_height);
                const uint64_t dist_opt =
                    spatial_full_distortion_kernel_avx2(input_,
                                                        input_stride_,
                                                        recon_,
                                                        recon_stride_,
                                                        area_width,
                                                        area_height);

                EXPECT_EQ(dist_org, dist_opt)
                    << area_width << "x" << area_height;
            }
        }
    }
}

void SpatialFullDistortionTest::RunSpeedTest() {
    uint64_t dist_org, dist_opt;
    double time_c, time_o;
    uint64_t start_time_seconds, start_time_useconds;
    uint64_t middle_time_seconds, middle_time_useconds;
    uint64_t finish_time_seconds, finish_time_useconds;

    init_data();

    for (uint32_t area_width = 4; area_width <= 128; area_width += 4) {
        const uint32_t area_height = area_width;
        const int num_loops = 1000000000 / (area_width * area_height);
        EbStartTime(&start_time_seconds, &start_time_useconds);

        for (int i = 0; i < num_loops; ++i) {
            dist_org = spatial_full_distortion_kernel_c(input_,
                                                        input_stride_,
                                                        recon_,
                                                        recon_stride_,
                                                        area_width,
                                                        area_height);
        }

        EbStartTime(&middle_time_seconds, &middle_time_useconds);

        for (int i = 0; i < num_loops; ++i) {
            dist_opt = spatial_full_distortion_kernel_avx2(input_,
                                                           input_stride_,
                                                           recon_,
                                                           recon_stride_,
                                                           area_width,
                                                           area_height);
        }
        EbStartTime(&finish_time_seconds, &finish_time_useconds);

        EXPECT_EQ(dist_org, dist_opt) << area_width << "x" << area_height;

        EbComputeOverallElapsedTimeMs(start_time_seconds,
                                      start_time_useconds,
                                      middle_time_seconds,
                                      middle_time_useconds,
                                      &time_c);
        EbComputeOverallElapsedTimeMs(middle_time_seconds,
                                      middle_time_useconds,
                                      finish_time_seconds,
                                      finish_time_useconds,
                                      &time_o);
        printf("Average Nanoseconds per Function Call\n");
        printf("    spatial_full_distortion_kernel_c   (%dx%d) : %6.2f\n",
               area_width,
               area_height,
               1000000 * time_c / num_loops);
        printf(
            "    spatial_full_distortion_kernel_avx2(%dx%d) : %6.2f   "
            "(Comparison: %5.2fx)\n",
            area_width,
            area_height,
            1000000 * time_o / num_loops,
            time_c / time_o);
    }
}

TEST_F(SpatialFullDistortionTest, CheckOutput) {
    RunCheckOutput();
}

TEST_F(SpatialFullDistortionTest, DISABLED_Speed) {
    RunSpeedTest();
}

};  // namespace
