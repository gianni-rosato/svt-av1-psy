/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file SvtAv1E2EParamsTest.cc
 *
 * @brief Implementation of encoder parameter coverage test in E2E test
 *
 * @author Cidana-Edmond
 *
 ******************************************************************************/

#include <map>
#include <cmath>
#include "gtest/gtest.h"
#include "SvtAv1E2EFramework.h"
#include "../api_test/params.h"
#include "RefDecoder.h"

/**
 * @brief SVT-AV1 encoder parameter coverage E2E test
 *
 * Test strategy:
 * Config SVT-AV1 encoder with individual parameter, run the
 * conformance test and analyze the Bitstream to check if the params
 * take effect.
 *
 * Expected result:
 * No error is reported in encoding progress. The reconstructed frame
 * data is same as the output frame from reference decoder.
 *
 * Test coverage:
 * Almost all the encoder parameters except frame_rate_numerator and
 * frame_rate_denominator.
 */

using namespace svt_av1_e2e_test;
using namespace svt_av1_e2e_test_vector;
using std::map;
using std::stoul;
using std::string;
using std::to_string;

/** copied from EbRateControlProcess.c */
static const uint8_t quantizer_to_qindex[] = {
    0,   4,   8,   12,  16,  20,  24,  28,  32,  36,  40,  44,  48,
    52,  56,  60,  64,  68,  72,  76,  80,  84,  88,  92,  96,  100,
    104, 108, 112, 116, 120, 124, 128, 132, 136, 140, 144, 148, 152,
    156, 160, 164, 168, 172, 176, 180, 184, 188, 192, 196, 200, 204,
    208, 212, 216, 220, 224, 228, 232, 236, 240, 244, 249, 255,
};

/** get qp value with the given qindex */
static uint32_t get_qp(const int16_t qindex) {
    if (qindex > 255) {
        printf("qindex is larger than 255!\n");
        return 63;
    }

    uint32_t qp = 0;
    for (const uint8_t index : quantizer_to_qindex) {
        if (index == qindex)
            return qp;
        else if (index > qindex) {
            if ((index - qindex) > (qindex - quantizer_to_qindex[qp - 1]))
                return qp - 1;
            break;
        }
        qp++;
    }
    return qp;
}

/* clang-format off */
static const std::vector<EncTestSetting> default_enc_settings = {
    // test intra period length
    {"IntraPeriodTest1", {{"IntraPeriod", "3"}}, default_test_vectors},

    // test different qp
    {"QpTest1",
     {{"RateControlMode", "0"}, {"QP", "20"}}, default_test_vectors},
    {"QpTest2",
     {{"RateControlMode", "0"}, {"QP", "32"}}, default_test_vectors},
    {"QpTest3",
     {{"RateControlMode", "0"}, {"QP", "44"}}, default_test_vectors},

    // test rc mode {2, 3}, with {1Mbps, 0.75Mbps, 0.5Mbps} setting with
    // 480p
    {"RcTest1",
     {{"RateControlMode", "2"}, {"TargetBitRate", "1000000"}},
     default_test_vectors},
    {"RcTest2",
     {{"RateControlMode", "2"}, {"TargetBitRate", "750000"}},
     res_480p_test_vectors},
    {"RcTest3",
     {{"RateControlMode", "2"}, {"TargetBitRate", "500000"}},
     res_480p_test_vectors},
    {"RcTest4",
     {{"RateControlMode", "1"}, {"TargetBitRate", "1000000"}},
     res_480p_test_vectors},
    {"RcTest5",
     {{"RateControlMode", "1"}, {"TargetBitRate", "750000"}},
     res_480p_test_vectors},
    {"RcTest6",
     {{"RateControlMode", "1"}, {"TargetBitRate", "500000"}},
     res_480p_test_vectors},

    // test high bitrate with big min_qp, or low bitrate with small max_qp
    {"RcQpTest1",
     {{"RateControlMode", "2"}, {"TargetBitRate", "1000000"}, {"MinQpAllowed", "20"}},
     res_480p_test_vectors},
    {"RcQpTest2",
     {{"RateControlMode", "2"}, {"TargetBitRate", "500000"}, {"MaxQpAllowed", "50"}},
     res_480p_test_vectors},
    {"RcQpTest3",
     {{"RateControlMode", "2"}, {"TargetBitRate", "750000"}, {"MaxQpAllowed", "50"}, {"MinQpAllowed", "20"}},
     res_480p_test_vectors},
    {"RcQpTest4",
     {{"RateControlMode", "1"}, {"TargetBitRate", "1000000"}, {"MinQpAllowed", "20"}},
     res_480p_test_vectors},
    {"RcQpTest5",
     {{"RateControlMode", "1"}, {"TargetBitRate", "500000"}, {"MaxQpAllowed", "50"}},
     res_480p_test_vectors},
    {"RcQpTest6",
     {{"RateControlMode", "1"}, {"TargetBitRate", "750000"}, {"MaxQpAllowed", "50"}, {"MinQpAllowed", "20"}},
     res_480p_test_vectors},
};
/* clang-format on */

class CodingOptionTest : public SvtAv1E2ETestFramework {
  public:
    void config_test() override {
        enable_recon = true;
        enable_decoder = true;
        enable_analyzer = true;
        enable_config = true;
        SvtAv1E2ETestFramework::config_test();
    }

    void post_process() override {
        if (refer_dec_) {
            RefDecoder::StreamInfo *stream_info = refer_dec_->get_stream_info();
            validate_enc_setting(stream_info);
        }
    }

  protected:
    void validate_enc_setting(RefDecoder::StreamInfo *stream_info) {
        EbSvtAv1EncConfiguration *config = &av1enc_ctx_.enc_params;

        // check profile, level and tier
        EXPECT_EQ(config->profile, stream_info->profile)
            << "config profile: " << config->profile << "got "
            << stream_info->profile;

        // verify the superblock size
        EXPECT_EQ(config->super_block_size, stream_info->sb_size)
            << "config sb size: " << config->super_block_size << " got "
            << stream_info->sb_size;

        // Verify bit depth
        EXPECT_EQ(config->encoder_bit_depth, stream_info->bit_depth)
            << "config bitdepth: " << config->encoder_bit_depth << " got "
            << stream_info->bit_depth;

        // verify the color format
        EXPECT_EQ(config->encoder_color_format,
                  setup_video_format(stream_info->format))
            << "color format is mismatch";

        if (config->intra_period_length > 0) {
            EXPECT_EQ(config->intra_period_length,
                      stream_info->max_intra_period)
                << "config intra period " << config->intra_period_length
                << " got " << stream_info->max_intra_period;
        }

        // verify QP Setting
        uint32_t actual_min_qp = get_qp(stream_info->min_qindex);
        uint32_t actual_max_qp = get_qp(stream_info->max_qindex);
        EXPECT_LE(config->min_qp_allowed, actual_min_qp)
            << "Min qp allowd " << config->min_qp_allowed << " actual "
            << actual_min_qp;
        EXPECT_GE(config->max_qp_allowed, actual_max_qp)
            << "Max qp allowd " << config->max_qp_allowed << " actual "
            << actual_max_qp;
        if (config->rate_control_mode == 0) {
            EXPECT_EQ(actual_min_qp, actual_max_qp)
                << "QP fluctuate in const qp mode";
        }

        // verify the bitrate
        if (config->rate_control_mode == 3) {
            uint32_t avg_bit_rate =
                (config->frame_rate > 1000 ? config->frame_rate >> 16
                                           : config->frame_rate) *
                stream_info->frame_bit_rate;
            printf("%d--%d\n", config->target_bit_rate, avg_bit_rate);
            EXPECT_GE(config->target_bit_rate, avg_bit_rate)
                << "target bit-rate is less than actual: "
                << config->target_bit_rate << "--" << avg_bit_rate;
        }

        // verify tile row and tile column
        uint32_t expect_cols =
            (uint32_t)((video_src_->get_width_with_padding() >> 2) /
                       (1 << config->tile_columns));
        uint32_t expect_rows =
            (uint32_t)((video_src_->get_height_with_padding() >> 2) /
                       (1 << config->tile_rows));
        printf("expect_cols %d, expect_rows %d\n", expect_cols, expect_rows);
        printf("tile_cols %d, tile_rows %d\n",
               stream_info->tile_cols,
               stream_info->tile_rows);
        EXPECT_EQ(expect_cols, stream_info->tile_cols)
            << "Tile columns " << stream_info->tile_cols << " actual"
            << expect_cols;
        EXPECT_EQ(expect_rows, stream_info->tile_rows)
            << "Tile rows " << stream_info->tile_rows << " actual"
            << expect_rows;

        //
        // Verify the coding tools by checking the sps header
        //
        EXPECT_EQ(stream_info->enable_warped_motion,
                  config->enable_warped_motion);
    }

    bool is_valid_profile_setting() {
        /** check the color format according to spec 6.4.1 */
        if (av1enc_ctx_.enc_params.profile == 0) {
            /** main profile requires YUV420 or YUV400 Annex A */
            if (av1enc_ctx_.enc_params.encoder_bit_depth == 12)
                return false;
            if (av1enc_ctx_.enc_params.encoder_color_format != EB_YUV420 &&
                av1enc_ctx_.enc_params.encoder_color_format != EB_YUV400) {
                return false;
            }
        } else if (av1enc_ctx_.enc_params.profile == 1) {
            /** high profile requires 8bit/10bit YUV444 */
            if (av1enc_ctx_.enc_params.encoder_bit_depth == 12)
                return false;
            if (av1enc_ctx_.enc_params.encoder_color_format != EB_YUV444)
                return false;
        } else if (av1enc_ctx_.enc_params.profile == 2) {
            /** professional profile requires 8-bit/10-bit YUV422 or 12-bit
             *  YUV400, YUV420, YUV422 and YUV444
             */
            if (av1enc_ctx_.enc_params.encoder_bit_depth != 12 &&
                av1enc_ctx_.enc_params.encoder_color_format != EB_YUV422) {
                return false;
            }
        }
        return true;
    }
};

TEST_P(CodingOptionTest, CheckEncOptionsUsingBitstream) {
    run_death_test();
}

INSTANTIATE_TEST_CASE_P(SvtAv1, CodingOptionTest,
                        ::testing::ValuesIn(default_enc_settings),
                        EncTestSetting::GetSettingName);
