/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file SvtAv1E2ETest.cc
 *
 * @brief Impelmentation of SVT-AV1 encoder E2E test
 *
 * @author Cidana-Edmond, Cidana-Wenyao
 *
 ******************************************************************************/

#include "EbSvtAv1Enc.h"
#include "gtest/gtest.h"
#include "SvtAv1E2EFramework.h"

using namespace svt_av1_e2e_test;
using namespace svt_av1_e2e_test_vector;
using std::string;

// generate all available enc mode settings
const std::vector<EncTestSetting> generate_enc_mode_settings() {
    string test = "EncModeTest";
    std::vector<EncTestSetting> settings;
    for (int i = 0; i <= MAX_ENC_PRESET; ++i) {
        string idx = std::to_string(i);
        string name = test + idx;
        EncTestSetting setting{
            name, {{"EncoderMode", idx}}, default_test_vectors};
        settings.push_back(setting);
    }

    return settings;
}

/**
 * @brief SVT-AV1 encoder simple E2E test
 *
 * Test strategy:
 * Setup SVT-AV1 encoder with default parameter, and encode the input YUV data
 * frames.
 *
 * Expected result:
 * No crash should occur in encoding progress. The output compressed data
 * is complete.
 *
 * Test coverage:
 * All test vectors
 */
class CrashDeathTest : public SvtAv1E2ETestFramework {
  protected:
    void config_test() override {
        enable_stat = true;
        enable_config = true;
        SvtAv1E2ETestFramework::config_test();
    }
};

TEST_P(CrashDeathTest, NotCrashTest) {
    run_death_test();
}

INSTANTIATE_TEST_CASE_P(SvtAv1, CrashDeathTest,
                        ::testing::ValuesIn(generate_enc_mode_settings()),
                        EncTestSetting::GetSettingName);

/**
 * @brief SVT-AV1 encoder E2E test with comparing the reconstructed frame with
 * output frame from decoder buffer list
 *
 * Test strategy:
 * Setup SVT-AV1 encoder with different parameter, and encode the input YUV data
 * frames. Collect the reconstructed frames and compared them with reference
 * decoder output.
 *
 * Expected result:
 * No error is reported in encoding progress. The reconstructed frame
 * data is same as the output frame from reference decoder.
 *
 * Test coverage:
 * All test vectors
 */
class ConformanceDeathTest : public SvtAv1E2ETestFramework {
  protected:
    void config_test() override {
        enable_decoder = true;
        enable_recon = true;
        enable_stat = true;
        enable_config = true;
        SvtAv1E2ETestFramework::config_test();
    }
};

TEST_P(ConformanceDeathTest, DefaultSettingTest) {
    run_death_test();
}

/* clang-format off */
std::vector<TestVideoVector> parkjoy = {
    std::make_tuple("park_joy_90p_8_420.y4m", Y4M_VIDEO_FILE, IMG_FMT_420, 160,
                    90, 8, 0, 0, 0),
};

static const std::vector<EncTestSetting> default_enc_settings = {
    {"EncModeTest1", {{"EncoderMode", "0"}}, default_test_vectors},
    {"EncModeTest2", {{"EncoderMode", "3"}}, default_test_vectors},
    {"EncModeTest3", {{"EncoderMode", "5"}}, default_test_vectors},
    {"EncModeTest4", {{"EncoderMode", "8"}}, default_test_vectors},

    // test intra period, default is -2;
    {"IntraPeriodTest1", {{"IntraPeriod", "-1"}}, default_test_vectors},
    {"IntraPeriodTest2", {{"IntraPeriod", "10"}}, default_test_vectors},

    // TODO: add intra_refresh_type, hierarchical_levels, pred_structure

    // test qps, default is 50
    {"QpTest1", {{"RateControlMode", "0"}, {"QP", "0"}}, default_test_vectors},
    {"QpTest2", {{"RateControlMode", "0"}, {"QP", "10"}}, default_test_vectors},
    {"QpTest3", {{"RateControlMode", "0"}, {"QP", "20"}}, default_test_vectors},
    {"QpTest4", {{"RateControlMode", "0"}, {"QP", "32"}}, default_test_vectors},
    {"QpTest5", {{"RateControlMode", "0"}, {"QP", "44"}}, default_test_vectors},
    {"QpTest6", {{"RateControlMode", "0"}, {"QP", "63"}}, default_test_vectors},

    // test disable_dlf_flag, default is 0
    {"DlfTest1", {{"LoopFilterDisable", "1"}}, default_test_vectors},

    // test film_grain_denoise_strength, default is 0
    {"FilmGrainTest1", {{"FilmGrain", "1"}}, default_test_vectors},

    // Skip enable_denoise_flag, enable_warped_motion, in_loop_me_flag
    // partition_depth and ext_block_flag, since they are
    // not used in encoder;

    // test constrained intra, default is 0
    {"ConstrainIntraTest1", {{"ConstrainedIntra", "1"}}, default_test_vectors},

    // test rate control modes, default is 0, 1 and 2 is supported
    {"RateControlTest1", {{"RateControlMode", "2"}}, default_test_vectors},
    {"RateControlTest2", {{"RateControlMode", "1"}}, default_test_vectors},

    // test scene change detection, default is 1
    {"SCDTest1", {{"SceneChangeDetection", "0"}}, default_test_vectors},

    // test look_ahead_distance, try other values than default value,
    // the default value should be 2 * 4 + 1 = 9
    {"LookAheadTest", {{"RateControlMode", "0"}, {"LookAheadDistance", "6"}, {"IntraPeriod", "4"}},
     default_test_vectors},

    // test ScreenContentMode, default 2 auto detection mode;
    {"ScreenToolTest1", {{"ScreenContentMode", "0"}}, default_test_vectors},
    {"ScreenToolTest2", {{"ScreenContentMode", "1"}, {"EncoderMode", "1"}}, screen_test_vectors},

    // test enable_adaptive_quantization, default is 0
    {"AdapQTest1", {{"AdaptiveQuantization", "1"}}, default_test_vectors},

    // test enable_altrefs, defalt is 1;
    {"AltrefTest1", {{"EnableAltRefs", "0"}}, default_test_vectors},

    // test tile settings
    {"TileTest1", {{"TileRow", "1"}}, default_test_vectors},
    {"TileTest2", {{"TileCol", "1"}}, default_test_vectors},
    {"TileTest3", {{"TileCol", "1"}, {"TileRow", "1"}}, default_test_vectors},

    {"SpeedControlTest1", {{"speed_control_flag", "1"}}, default_test_vectors},

    // Validate by setting a low bitrate and MaxQpAllowed, push the encoder to producing
    // large partitions.
    {"IncompleteSbTest1",
     {{"RateControlMode", "2"}, {"TargetBitRate", "100000"}, {"MaxQpAllowed", "40"}},
     parkjoy},
    // Validate by setting a high bitrate and MinQpAllowed, push the encoder to producing
    // small partitions.
    {"IncompleteSbTest2",
     {{"RateControlMode", "2"}, {"TargetBitRate", "1000000"}, {"MixQpAllowed", "10"}},
     parkjoy},

    // test pallete mode
    {"PaletteModeTest1",
     {{"PaletteMode", "1"}, {"ScreenContentMode", "1"}, {"EncoderMode", "1"}},
     screen_test_vectors},
    {"PaletteModeTest2",
     {{"PaletteMode", "2"}, {"ScreenContentMode", "1"}, {"EncoderMode", "1"}},
     screen_test_vectors},
    {"PaletteModeTest3",
     {{"PaletteMode", "3"}, {"ScreenContentMode", "1"}, {"EncoderMode", "1"}},
     screen_test_vectors},
    {"PaletteModeTest4",
     {{"PaletteMode", "4"}, {"ScreenContentMode", "1"}, {"EncoderMode", "1"}},
     screen_test_vectors},
    {"PaletteModeTest5",
     {{"PaletteMode", "5"}, {"ScreenContentMode", "1"}, {"EncoderMode", "1"}},
     screen_test_vectors},
    {"PaletteModeTest6",
     {{"PaletteMode", "6"}, {"ScreenContentMode", "1"}, {"EncoderMode", "1"}},
     screen_test_vectors},
    {"PaletteModeTest7",
     {{"PaletteMode", "0"}, {"ScreenContentMode", "1"}, {"EncoderMode", "1"}},
     screen_test_vectors},

    // test by using a dummy source of color bar
    {"DummySrcTest1", {{"EncoderMode", "8"}}, dummy_test_vectors},

    // only 420 input is supported
    //{"DummySrcTest2", {{"EncoderMode", "8"}, {"Profile", "2"}}, dummy_422_test_vectors},
    //{"DummySrcTest3", {{"EncoderMode", "8"}, {"Profile", "1"}}, dummy_444_test_vectors},
};

/* clang-format on */
INSTANTIATE_TEST_CASE_P(SvtAv1, ConformanceDeathTest,
                        ::testing::ValuesIn(default_enc_settings),
                        EncTestSetting::GetSettingName);

class LongtimeConformanceTest : public ConformanceDeathTest {};

TEST_P(LongtimeConformanceTest, DISABLED_LongtimeTest) {
    run_death_test();
}

INSTANTIATE_TEST_CASE_P(SvtAv1, LongtimeConformanceTest,
                        ::testing::ValuesIn(generate_vector_from_config(
                            "longtime_comformance_test.cfg")),
                        EncTestSetting::GetSettingName);
/**
 * @brief SVT-AV1 encoder E2E test with comparing the reconstructed frame with
 * output frame from decoder buffer list when the tile is inverted to prove
 * tile independence.
 *
 * Test strategy:
 * Setup SVT-AV1 encoder with different tile parameter, and encode the input YUV
 * data frames. Do the decode in inverted tile ordering and collect the
 * reconstructed frames and compared them with reference decoder output.
 *
 * Expected result:
 * No error is reported in encoding progress. The reconstructed frame data is
 * same as the output frame from reference decoder,which proved tiles are
 * considered independent and the test passes.
 *
 * Test coverage:
 * All test vectors of 640*480 */

class TileIndependenceTest : public SvtAv1E2ETestFramework {
  protected:
    void config_test() override {
        enable_decoder = true;
        enable_recon = true;
        enable_stat = true;
        enable_config = true;
        enable_invert_tile_decoding = true;
        SvtAv1E2ETestFramework::config_test();
    }
};

TEST_P(TileIndependenceTest, TileTest) {
    run_death_test();
}

static const std::vector<EncTestSetting> tile_settings = {
    {"TileTest1", {{"TileCol", "0"}, {"TileRow", "0"}}, default_test_vectors},
    {"TileTest2", {{"TileCol", "0"}, {"TileRow", "1"}}, default_test_vectors},
    {"TileTest3", {{"TileCol", "1"}, {"TileRow", "0"}}, default_test_vectors},
    {"TileTest4", {{"TileCol", "1"}, {"TileRow", "1"}}, default_test_vectors}};

INSTANTIATE_TEST_CASE_P(TILETEST, TileIndependenceTest,
                        ::testing::ValuesIn(tile_settings),
                        EncTestSetting::GetSettingName);
