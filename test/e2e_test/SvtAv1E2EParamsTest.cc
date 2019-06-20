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
#include "ConfigEncoder.h"
/**
 * @brief SVT-AV1 encoder parameter coverage E2E test
 *
 * Test strategy:
 * Config SVT-AV1 encoder with individual parameter, run the
 * conformance test and analyze the bitstream to check if the params
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
static uint32_t get_qp(const uint8_t qindex) {
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

static const std::vector<EncTestSetting> default_enc_settings = {
    {"IntraPeriodTest1",
     {{"EncoderMode", "7"}, {"IntraPeriod", "3"}},
     default_test_vectors},
    {"EncModeTest1", {{"EncoderMode", "8"}}, default_test_vectors},
    {"SpeedControlTest1", {{"SpeedControlFlag", "1"}}, default_test_vectors}};

class CodingOptionTest : public SvtAv1E2ETestFramework {
  public:
    CodingOptionTest() {
        enc_config_ = create_enc_config();
    }

    virtual ~CodingOptionTest() {
        release_enc_config(enc_config_);
    }

    void config_test() override {
        enable_recon = true;
        enable_decoder = true;
        enable_analyzer = true;
        // iterate the mappings and update config
        for (auto &x : enc_setting.setting) {
            set_enc_config(enc_config_, x.first.c_str(), x.second.c_str());
            printf("EncSetting: %s = %s\n", x.first.c_str(), x.second.c_str());
        }
    }

    void update_enc_setting() override {
        copy_enc_param(&av1enc_ctx_.enc_params, enc_config_);
        setup_src_param(video_src_, av1enc_ctx_.enc_params);
        if (recon_queue_)
            av1enc_ctx_.enc_params.recon_enabled = 1;
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
        ASSERT_EQ(config->profile, stream_info->profile)
            << "config profile: " << config->profile << "got "
            << stream_info->profile;

        // verify the superblock size
        ASSERT_EQ(config->super_block_size, stream_info->sb_size)
            << "config sb size: " << config->super_block_size << " got "
            << stream_info->sb_size;

        // Verify bit depth
        ASSERT_EQ(config->encoder_bit_depth, stream_info->bit_depth)
            << "config bitdepth: " << config->encoder_bit_depth << " got "
            << stream_info->bit_depth;
        // TODO: verify the color format

        if (config->intra_period_length > 0) {
            ASSERT_EQ(config->intra_period_length,
                      stream_info->max_intra_period)
                << "config intra period " << config->intra_period_length
                << " got " << stream_info->max_intra_period;
        }

        // verify QP Setting
        int actual_min_qp = get_qp(stream_info->min_qindex);
        int actual_max_qp = get_qp(stream_info->max_qindex);
        ASSERT_LE(config->min_qp_allowed, actual_min_qp)
            << "Min qp allowd " << config->min_qp_allowed << " actual "
            << actual_min_qp;
        ASSERT_GE(config->max_qp_allowed, actual_max_qp)
            << "Max qp allowd " << config->max_qp_allowed << " actual "
            << actual_max_qp;
        if (config->rate_control_mode == 0)
            ASSERT_EQ(actual_min_qp, actual_max_qp)
                << "QP fluctuate in const qp mode";

        // TODO: verify the bitrate

        // TODO: verify tileRow and TileCol

        //
        // Verify the coding tools.
        //
        // Verify ext_block_flag
        // It fails if non-square partitions is not enabled in
        // encoder but found in bitstream. However if non-square partitions is
        // enabled, it can not gurantee that the bitstream must have non-square
        // partitions.
        ASSERT_GE(config->ext_block_flag, stream_info->ext_block_flag);

        // TODO: Verify enable_warped_motion, LoopFilterDisable
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

  protected:
    void *enc_config_;
};

TEST_P(CodingOptionTest, CheckEncOptionsUsingBitstream) {
    run_death_test();
}

INSTANTIATE_TEST_CASE_P(SvtAv1, CodingOptionTest,
                        ::testing::ValuesIn(default_enc_settings),
                        GetSettingName);
