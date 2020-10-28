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
 * @file TemporalFilterTest.cc
 *
 * @brief Unit test for Temporal Filter functions:
 * - svt_av1_apply_temporal_filter_sse4_1
 * - svt_av1_highbd_apply_temporal_filter_sse4_1
 *
 * @author Cidana-Ivy
 *
 ******************************************************************************/
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <limits.h>

#include "EbPictureOperators.h"
#include "EbEncIntraPrediction.h"
#include "EbTemporalFiltering.h"
#include "random.h"
#include "util.h"
#include "aom_dsp_rtcd.h"
extern "C" {

void svt_av1_highbd_apply_temporal_filter_sse4_1(
    const uint16_t *y_src, int y_src_stride, const uint16_t *y_pre,
    int y_pre_stride, const uint16_t *u_src, const uint16_t *v_src,
    int uv_src_stride, const uint16_t *u_pre, const uint16_t *v_pre,
    int uv_pre_stride, unsigned int block_width, unsigned int block_height,
    int ss_x, int ss_y, int strength, const int *blk_fw, int use_whole_blk,
    uint32_t *y_accum, uint16_t *y_count, uint32_t *u_accum, uint16_t *u_count,
    uint32_t *v_accum, uint16_t *v_count);
}

using svt_av1_test_tool::SVTRandom;  // to generate the random

namespace {
typedef enum { FW_MIN, FW_MAX, FW_MID, FW_RANDOM } FilterWeightPattern;
FilterWeightPattern FW_PATTERNS[] = {FW_MIN, FW_MAX, FW_MID, FW_RANDOM};

typedef std::tuple<int16_t, int16_t> BlockIndex;
BlockIndex TEST_BLOCK[] = {
    BlockIndex(0, 0), BlockIndex(0, 1), BlockIndex(1, 0), BlockIndex(1, 1)};

static const uint32_t index_16x16_from_subindexes[4][4] = {
    {0, 1, 4, 5}, {2, 3, 6, 7}, {8, 9, 12, 13}, {10, 11, 14, 15}};

const int stride_pred_[COLOR_CHANNELS] = {BW, BW >> 1, BW >> 1};
const int ALTREF_STRENGTH[] = {0, 1, 2, 3, 4, 5, 6};

typedef std::tuple<BlockIndex, FilterWeightPattern, int> TestParam;

/**
 * @brief Unit test for Temporal Filter function include:
 *  - svt_av1_apply_temporal_filter_sse4_1
 *  - svt_av1_highbd_apply_temporal_filter_sse4_1
 *
 * Test strategy:
 * This test case combines different blk index,alt-ref strength,filter weight
 * and different test pattern. compare the accum and count from c function and
 * sse function.
 *
 * Expect result:
 * accum and count from c function and sse function are
 * equal.
 *
 * Test cases:
 *  block_index:row{0,1} x col{0,1}
 *  alt-ref strength{0,1,2,3,4,5,6}
 *  filter weight{0,1,2}
 *  Test vector pattern: radom pixel value of pred and src
 *
 */
class TemporalFilterTest : public ::testing::Test,
                           public ::testing::WithParamInterface<TestParam> {
  public:
    TemporalFilterTest()
        : fw_pattern_(TEST_GET_PARAM(1)),
          block_col_(std::get<0>(TEST_GET_PARAM(0))),
          block_row_(std::get<1>(TEST_GET_PARAM(0))),
          altref_strength_(TEST_GET_PARAM(2)) {
        ss_x_ = ss_y_ = 1;
        test_size_ = BLK_PELS;
        block_width_ = block_height_ = 32;

        width_ = height_ = 64;  // only 64x64 area used in test func
        stride_[0] = width_;
        stride_[1] = width_ >> 1;
        stride_[2] = width_ >> 1;

        source = reinterpret_cast<uint8_t *>(
            svt_aom_memalign(32, BLK_PELS * COLOR_CHANNELS));
        predictor = reinterpret_cast<uint8_t *>(
            svt_aom_memalign(32, BLK_PELS * COLOR_CHANNELS));
        hbd_source = reinterpret_cast<uint16_t *>(
            svt_aom_memalign(32, sizeof(uint16_t) * BLK_PELS * COLOR_CHANNELS));
        hbd_predictor = reinterpret_cast<uint16_t *>(
            svt_aom_memalign(32, sizeof(uint16_t) * BLK_PELS * COLOR_CHANNELS));
        accumulator_1_ = reinterpret_cast<uint32_t *>(
            svt_aom_memalign(16, sizeof(uint32_t) * BLK_PELS * COLOR_CHANNELS));
        accumulator_2_ = reinterpret_cast<uint32_t *>(
            svt_aom_memalign(16, sizeof(uint32_t) * BLK_PELS * COLOR_CHANNELS));
        accumulator_3_ = reinterpret_cast<uint32_t *>(
            svt_aom_memalign(16, sizeof(uint32_t) * BLK_PELS * COLOR_CHANNELS));
        counter_1_ = reinterpret_cast<uint16_t *>(
            svt_aom_memalign(16, sizeof(uint16_t) * BLK_PELS * COLOR_CHANNELS));
        counter_2_ = reinterpret_cast<uint16_t *>(
            svt_aom_memalign(16, sizeof(uint16_t) * BLK_PELS * COLOR_CHANNELS));
        counter_3_ = reinterpret_cast<uint16_t *>(
            svt_aom_memalign(16, sizeof(uint16_t) * BLK_PELS * COLOR_CHANNELS));
        memset(accumulator_1_,
               0,
               BLK_PELS * COLOR_CHANNELS * sizeof(accumulator_1_[0]));
        memset(
            counter_1_, 0, BLK_PELS * COLOR_CHANNELS * sizeof(counter_1_[0]));
        memset(accumulator_2_,
               0,
               BLK_PELS * COLOR_CHANNELS * sizeof(accumulator_2_[0]));
        memset(
            counter_2_, 0, BLK_PELS * COLOR_CHANNELS * sizeof(counter_2_[0]));
        memset(accumulator_3_,
               0,
               BLK_PELS * COLOR_CHANNELS * sizeof(accumulator_3_[0]));
        memset(
            counter_3_, 0, BLK_PELS * COLOR_CHANNELS * sizeof(counter_3_[0]));
        accum_1_[0] = accumulator_1_;
        accum_1_[1] = accumulator_1_ + BLK_PELS;
        accum_1_[2] = accumulator_1_ + (BLK_PELS << 1);
        accum_2_[0] = accumulator_2_;
        accum_2_[1] = accumulator_2_ + BLK_PELS;
        accum_2_[2] = accumulator_2_ + (BLK_PELS << 1);
        accum_3_[0] = accumulator_3_;
        accum_3_[1] = accumulator_3_ + BLK_PELS;
        accum_3_[2] = accumulator_3_ + (BLK_PELS << 1);
        count_1_[0] = counter_1_;
        count_1_[1] = counter_1_ + BLK_PELS;
        count_1_[2] = counter_1_ + (BLK_PELS << 1);
        count_2_[0] = counter_2_;
        count_2_[1] = counter_2_ + BLK_PELS;
        count_2_[2] = counter_2_ + (BLK_PELS << 1);
        count_3_[0] = counter_3_;
        count_3_[1] = counter_3_ + BLK_PELS;
        count_3_[2] = counter_3_ + (BLK_PELS << 1);
        pred_[0] = predictor;
        pred_[1] = predictor + BLK_PELS;
        pred_[2] = predictor + (BLK_PELS << 1);
        src_[0] = source;
        src_[1] = source + BLK_PELS;
        src_[2] = source + (BLK_PELS << 1);
        hbd_pred_[0] = hbd_predictor;
        hbd_pred_[1] = hbd_predictor + BLK_PELS;
        hbd_pred_[2] = hbd_predictor + (BLK_PELS << 1);
        hbd_src_[0] = hbd_source;
        hbd_src_[1] = hbd_source + BLK_PELS;
        hbd_src_[2] = hbd_source + (BLK_PELS << 1);
    }

    ~TemporalFilterTest() {
        if (source)
            svt_aom_free(source);
        if (predictor)
            svt_aom_free(predictor);
        if (hbd_source)
            svt_aom_free(hbd_source);
        if (hbd_predictor)
            svt_aom_free(hbd_predictor);
        if (accumulator_1_)
            svt_aom_free(accumulator_1_);
        if (accumulator_2_)
            svt_aom_free(accumulator_2_);
        if (accumulator_3_)
            svt_aom_free(accumulator_3_);
        if (counter_1_)
            svt_aom_free(counter_1_);
        if (counter_2_)
            svt_aom_free(counter_2_);
        if (counter_3_)
            svt_aom_free(counter_3_);
    }

  protected:
    void init_param() {
        uint8_t offset_src_buffer_y =
            block_row_ * (BH >> 1) * stride_[C_Y] + block_col_ * (BW >> 1);
        uint8_t offset_src_buffer_u =
            block_row_ * (BH >> 2) * stride_[C_U] + block_col_ * (BW >> 2);
        uint8_t offset_src_buffer_v =
            block_row_ * (BH >> 2) * stride_[C_V] + block_col_ * (BW >> 2);

        uint8_t offset_block_buffer_y =
            block_row_ * 32 * stride_pred_[C_Y] + block_col_ * 32;
        uint8_t offset_block_buffer_u =
            block_row_ * 16 * stride_pred_[C_U] + block_col_ * 16;
        uint8_t offset_block_buffer_v =
            block_row_ * 16 * stride_pred_[C_V] + block_col_ * 16;

        int idx_32x32 = block_row_ * 2 + block_col_;

        for (int ifw = 0; ifw < 4; ifw++) {
            int ifw_index = index_16x16_from_subindexes[idx_32x32][ifw];

            blk_fw_32x32_[ifw] = blk_fw_[ifw_index];
        }

        src_ptr_[C_Y] = src_[C_Y] + offset_src_buffer_y;
        src_ptr_[C_U] = src_[C_U] + offset_src_buffer_u;
        src_ptr_[C_V] = src_[C_V] + offset_src_buffer_v;

        pred_ptr_[C_Y] = pred_[C_Y] + offset_block_buffer_y;
        pred_ptr_[C_U] = pred_[C_U] + offset_block_buffer_u;
        pred_ptr_[C_V] = pred_[C_V] + offset_block_buffer_v;

        hbd_src_ptr_[C_Y] = hbd_src_[C_Y] + offset_src_buffer_y;
        hbd_src_ptr_[C_U] = hbd_src_[C_U] + offset_src_buffer_u;
        hbd_src_ptr_[C_V] = hbd_src_[C_V] + offset_src_buffer_v;

        hbd_pred_ptr_[C_Y] = hbd_pred_[C_Y] + offset_block_buffer_y;
        hbd_pred_ptr_[C_U] = hbd_pred_[C_U] + offset_block_buffer_u;
        hbd_pred_ptr_[C_V] = hbd_pred_[C_V] + offset_block_buffer_v;

        accum_ptr1_[C_Y] = accum_1_[C_Y] + offset_block_buffer_y;
        accum_ptr1_[C_U] = accum_1_[C_U] + offset_block_buffer_u;
        accum_ptr1_[C_V] = accum_1_[C_V] + offset_block_buffer_v;
        accum_ptr2_[C_Y] = accum_2_[C_Y] + offset_block_buffer_y;
        accum_ptr2_[C_U] = accum_2_[C_U] + offset_block_buffer_u;
        accum_ptr2_[C_V] = accum_2_[C_V] + offset_block_buffer_v;
        accum_ptr3_[C_Y] = accum_3_[C_Y] + offset_block_buffer_y;
        accum_ptr3_[C_U] = accum_3_[C_U] + offset_block_buffer_u;
        accum_ptr3_[C_V] = accum_3_[C_V] + offset_block_buffer_v;
        count_ptr1_[C_Y] = count_1_[C_Y] + offset_block_buffer_y;
        count_ptr1_[C_U] = count_1_[C_U] + offset_block_buffer_u;
        count_ptr1_[C_V] = count_1_[C_V] + offset_block_buffer_v;
        count_ptr2_[C_Y] = count_2_[C_Y] + offset_block_buffer_y;
        count_ptr2_[C_U] = count_2_[C_U] + offset_block_buffer_u;
        count_ptr2_[C_V] = count_2_[C_V] + offset_block_buffer_v;
        count_ptr3_[C_Y] = count_3_[C_Y] + offset_block_buffer_y;
        count_ptr3_[C_U] = count_3_[C_U] + offset_block_buffer_u;
        count_ptr3_[C_V] = count_3_[C_V] + offset_block_buffer_v;
    }

    void populate_list_with_value(int *list, int nelements, const int value) {
        for (int i = 0; i < nelements; i++)
            list[i] = value;
    }

    void prepare_data() {
        SVTRandom rnd_weight(0, 2);
        switch (fw_pattern_) {
        case FW_MIN: {
            populate_list_with_value(blk_fw_, 16, 0);
            break;
        }
        case FW_MAX: {
            populate_list_with_value(blk_fw_, 16, 2);
            break;
        }
        case FW_MID: {
            populate_list_with_value(blk_fw_, 16, 1);
            break;
        }
        case FW_RANDOM: {
            populate_list_with_value(blk_fw_, 16, rnd_weight.random());
            break;
        }
        default: break;
        }
    }

    void check_accum_output() {
        int fail_count = 0;
        for (uint32_t i = 0; i < BLK_PELS * COLOR_CHANNELS; i++) {
            if ((accumulator_1_[i] != accumulator_2_[i]) ||
                (accumulator_1_[i] != accumulator_3_[i]))
                fail_count++;
        }

        EXPECT_EQ(0, fail_count)
            << "compare result error"
            << "in test area for " << fail_count << "times";
    }

    void check_count_output() {
        int fail_count = 0;
        for (uint32_t i = 0; i < BLK_PELS; i++) {
            if ((counter_1_[i] != counter_2_[i]) ||
                (counter_1_[i] != counter_3_[i]))
                fail_count++;
        }
        EXPECT_EQ(0, fail_count)
            << "compare result error"
            << "in test area for " << fail_count << "times";
    }

    void run_test() {
        const int loops = 100;
        const uint8_t mask = (1 << 8) - 1;
        const uint16_t hbd_mask = (1 << 10) - 1;
        SVTRandom rnd_uint8(0, mask);
        SVTRandom rnd_uint16(0, hbd_mask);

        for (int j = 0; j < loops; j++) {
            for (uint32_t i = 0; i < BLK_PELS * COLOR_CHANNELS; i++) {
                predictor[i] = rnd_uint8.random();
                source[i] = rnd_uint8.random();
                hbd_predictor[i] = rnd_uint16.random();
                hbd_source[i] = rnd_uint16.random();
            }
            prepare_data();
            init_param();
            svt_av1_apply_filtering_c(src_ptr_[C_Y],
                                      stride_[C_Y],
                                      pred_ptr_[C_Y],
                                      stride_pred_[C_Y],
                                      src_ptr_[C_U],
                                      src_ptr_[C_V],
                                      stride_[C_U],
                                      pred_ptr_[C_U],
                                      pred_ptr_[C_V],
                                      stride_pred_[C_U],
                                      block_width_,
                                      block_height_,
                                      ss_x_,
                                      ss_y_,
                                      altref_strength_,
                                      blk_fw_32x32_,
                                      0,  // use_32x32
                                      accum_ptr1_[C_Y],
                                      count_ptr1_[C_Y],
                                      accum_ptr1_[C_U],
                                      count_ptr1_[C_U],
                                      accum_ptr1_[C_V],
                                      count_ptr1_[C_V]);
            svt_av1_apply_temporal_filter_sse4_1(src_ptr_[C_Y],
                                                 stride_[C_Y],
                                                 pred_ptr_[C_Y],
                                                 stride_pred_[C_Y],
                                                 src_ptr_[C_U],
                                                 src_ptr_[C_V],
                                                 stride_[C_U],
                                                 pred_ptr_[C_U],
                                                 pred_ptr_[C_V],
                                                 stride_pred_[C_U],
                                                 block_width_,
                                                 block_height_,
                                                 ss_x_,
                                                 ss_y_,
                                                 altref_strength_,
                                                 blk_fw_32x32_,
                                                 0,  // use_32x32
                                                 accum_ptr2_[C_Y],
                                                 count_ptr2_[C_Y],
                                                 accum_ptr2_[C_U],
                                                 count_ptr2_[C_U],
                                                 accum_ptr2_[C_V],
                                                 count_ptr2_[C_V]);
            svt_av1_highbd_apply_temporal_filter_sse4_1(hbd_src_ptr_[C_Y],
                                                        stride_[C_Y],
                                                        hbd_pred_ptr_[C_Y],
                                                        stride_pred_[C_Y],
                                                        hbd_src_ptr_[C_U],
                                                        hbd_src_ptr_[C_V],
                                                        stride_[C_U],
                                                        hbd_pred_ptr_[C_U],
                                                        hbd_pred_ptr_[C_V],
                                                        stride_pred_[C_U],
                                                        block_width_,
                                                        block_height_,
                                                        ss_x_,
                                                        ss_y_,
                                                        altref_strength_,
                                                        blk_fw_32x32_,
                                                        0,  // use_32x32
                                                        accum_ptr3_[C_Y],
                                                        count_ptr3_[C_Y],
                                                        accum_ptr3_[C_U],
                                                        count_ptr3_[C_U],
                                                        accum_ptr3_[C_V],
                                                        count_ptr3_[C_V]);
            check_accum_output();
            check_count_output();
            EXPECT_FALSE(HasFailure());
        }
    }

    FilterWeightPattern fw_pattern_;
    uint32_t block_width_, block_height_;
    int16_t block_col_, block_row_;
    uint8_t *source;
    uint8_t *predictor;
    uint16_t *hbd_source;
    uint16_t *hbd_predictor;
    uint32_t *accumulator_1_;
    uint16_t *counter_1_;
    uint32_t *accumulator_2_;
    uint16_t *counter_2_;
    uint32_t *accumulator_3_;
    uint16_t *counter_3_;
    uint8_t *src_[COLOR_CHANNELS];
    uint8_t *pred_[COLOR_CHANNELS];
    uint16_t *hbd_src_[COLOR_CHANNELS];
    uint16_t *hbd_pred_[COLOR_CHANNELS];
    uint32_t *accum_1_[COLOR_CHANNELS];
    uint16_t *count_1_[COLOR_CHANNELS];
    uint32_t *accum_2_[COLOR_CHANNELS];
    uint16_t *count_2_[COLOR_CHANNELS];
    uint32_t *accum_3_[COLOR_CHANNELS];
    uint16_t *count_3_[COLOR_CHANNELS];
    uint8_t *src_ptr_[COLOR_CHANNELS];
    uint8_t *pred_ptr_[COLOR_CHANNELS];
    uint16_t *hbd_src_ptr_[COLOR_CHANNELS];
    uint16_t *hbd_pred_ptr_[COLOR_CHANNELS];
    uint32_t *accum_ptr1_[COLOR_CHANNELS], *accum_ptr2_[COLOR_CHANNELS],
        *accum_ptr3_[COLOR_CHANNELS];
    uint16_t *count_ptr1_[COLOR_CHANNELS], *count_ptr2_[COLOR_CHANNELS],
        *count_ptr3_[COLOR_CHANNELS];
    uint8_t *src_altref_index_[COLOR_CHANNELS];
    int stride_[COLOR_CHANNELS];
    int16_t ss_x_, ss_y_;
    int altref_strength_;
    int blk_fw_[16];
    int blk_fw_32x32_[4];
    int16_t test_size_;
    uint32_t width_;
    uint32_t height_;
};

TEST_P(TemporalFilterTest, TemporalFilterTest) {
    run_test();
};

INSTANTIATE_TEST_CASE_P(
    TemporalFilter, TemporalFilterTest,
    ::testing::Combine(::testing::ValuesIn(TEST_BLOCK),
                       ::testing::ValuesIn(FW_PATTERNS),
                       ::testing::ValuesIn(ALTREF_STRENGTH)));

}  // namespace
