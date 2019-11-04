/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file SvtAv1EncParamsTest.cc
 *
 * @brief SVT-AV1 encoder parameter configuration test
 *
 * @author Cidana-Edmond
 *
 ******************************************************************************/

#include "EbSvtAv1Enc.h"
#include "gtest/gtest.h"
#include "params.h"
#include "SvtAv1EncApiTest.h"

using namespace svt_av1_test_params;
using namespace svt_av1_test;

namespace {
/**
 * @brief SVT-AV1 encoder parameter configuration test
 *
 * Test strategy:
 * Feed default values, vaild values and invalid values of individual param
 * to the encoder and check if encoder return correctly.
 *
 * Expected result:
 * For default value and valid values, encoder should return EB_ErrorNone.
 * For invalid value, encoder should return EB_ErrorBadParameter.
 *
 * Test coverage:
 * Almost all the encoder parameters except frame_rate_numerator and
 * frame_rate_denominator
 */

/** @breif EncParamTestBase is the basic class framework to test individual
 * parameter
 */
class EncParamTestBase : public ::testing::Test {
  public:
    EncParamTestBase() {
        memset(&ctxt_, 0, sizeof(ctxt_));
        param_name_str_ = "";
    }
    EncParamTestBase(std::string param_name) {
        memset(&ctxt_, 0, sizeof(ctxt_));
        param_name_str_ = param_name;
    }
    virtual ~EncParamTestBase() {
    }

  public:
    /** Interfaces */
    /** Interface of running parameter check with default value*/
    virtual void run_default_param_check() = 0;
    /** Interface of running parameter check with valid value*/
    virtual void run_valid_param_check() = 0;
    /** Interface of running parameter check with invalid value*/
    virtual void run_invalid_param_check() = 0;
    /** Interface of running parameter check with special value*/
    virtual void run_special_param_check() = 0;

  protected:
    // Sets up the test fixture.
    virtual void SetUp() override {
        // initialize encoder and get handle
        ASSERT_EQ(EB_ErrorNone,
                  eb_init_handle(&ctxt_.enc_handle, &ctxt_, &ctxt_.enc_params))
            << "eb_init_handle failed";
        // setup encoder parameters with all default
        ASSERT_NE(nullptr, ctxt_.enc_handle) << "enc_handle is invalid";
        // setup source width/height with default if not in test source_width or
        // source_height
        if (param_name_str_.compare("source_width")) {
            const int width = 1280;
            ctxt_.enc_params.source_width = width;
        }
        if (param_name_str_.compare("source_height")) {
            const int height = 720;
            ctxt_.enc_params.source_height = height;
        }
    }

    // Tears down the test fixture.
    virtual void TearDown() override {
        // TODO: eb_deinit_encoder should not be called here, for this test does
        // not call eb_init_encoder, but there is huge memory leak if only calls
        // eb_deinit_handle. please remmove it after we pass
        // EncApiTest-->repeat_normal_setup
        ASSERT_EQ(EB_ErrorNone, eb_deinit_encoder(ctxt_.enc_handle))
            << "eb_deinit_encoder failed";
        // destory encoder
        ASSERT_EQ(EB_ErrorNone, eb_deinit_handle(ctxt_.enc_handle))
            << "eb_deinit_handle failed";
    }

    /** setup some of the params with related params modified before set
     * to encoder */
    void config_enc_param() {
        // special cases for parameter
        if (!param_name_str_.compare("max_qp_allowed")) {
            ctxt_.enc_params.rate_control_mode = 1;
            ctxt_.enc_params.min_qp_allowed = MIN_QP_VALUE;
        } else if (!param_name_str_.compare("min_qp_allowed")) {
            ctxt_.enc_params.rate_control_mode = 1;
            ctxt_.enc_params.max_qp_allowed = MAX_QP_VALUE;
        } else if (!param_name_str_.compare("profile")) {
            if (ctxt_.enc_params.profile == 0) {
                /** profile(0) requires YUV420 */
                ctxt_.enc_params.encoder_color_format = EB_YUV420;
            } else if (ctxt_.enc_params.profile == 1) {
                /** profile(1) requires 8-bit YUV444 */
                ctxt_.enc_params.encoder_bit_depth = 8;
                ctxt_.enc_params.encoder_color_format = EB_YUV444;
            } else if (ctxt_.enc_params.profile == 2) {
                /** profile(2) requires 8-bit/10-bit YUV422 */
                ctxt_.enc_params.encoder_bit_depth = 8;
                ctxt_.enc_params.encoder_color_format = EB_YUV422;
            }
        } else if (!param_name_str_.compare("target_bit_rate")) {
            ctxt_.enc_params.rate_control_mode = 1;
        } else if (!param_name_str_.compare("injector_frame_rate")) {
            ctxt_.enc_params.speed_control_flag = 1;
        } else if (!param_name_str_.compare("altref_strength") ||
                   !param_name_str_.compare("altref_nframes")) {
            ctxt_.enc_params.enable_altrefs = EB_TRUE;
        }
    }

  protected:
    SvtAv1Context ctxt_;         /**< context of encoder */
    std::string param_name_str_; /**< name of parameter for test */
};

/** Marcro defininition of printing parameter name when in failed */
#define PRINT_PARAM_FATAL(p) \
    << "eb_svt_enc_set_parameter " << #p << ": " << (int)(p) << " failed"

/** Marcro defininition of printing 2 parameters name when in failed */
#define PRINT_2PARAM_FATAL(p1, p2)                                             \
    << "eb_svt_enc_set_parameter " << #p1 << ": " << (int)(p1) << " + " << #p2 \
    << ": " << (int)(p2) << " failed"

/** Marcro defininition of batch processing check for default, valid, invalid
 * and special parameter check*/
#define PARAM_TEST(param_test)               \
    TEST_F(param_test, run_paramter_check) { \
        run_default_param_check();           \
        run_valid_param_check();             \
        run_invalid_param_check();           \
        run_special_param_check();           \
    }

/** @breif This class is a template based on EncParamTestBase to test each
 * parameter
 */
#define DEFINE_PARAM_TEST_CLASS(test_name, param_name)                        \
    class test_name : public EncParamTestBase {                               \
      public:                                                                 \
        test_name() : EncParamTestBase(#param_name) {                         \
        }                                                                     \
        virtual void run_default_param_check() override {                     \
            EncParamTestBase::SetUp();                                        \
            ASSERT_EQ(ctxt_.enc_params.param_name,                            \
                      GET_DEFAULT_PARAM(param_name));                         \
            EncParamTestBase::TearDown();                                     \
        }                                                                     \
        virtual void run_valid_param_check() override {                       \
            for (size_t i = 0; i < SIZE_VALID_PARAM(param_name); ++i) {       \
                EncParamTestBase::SetUp();                                    \
                ctxt_.enc_params.param_name = GET_VALID_PARAM(param_name, i); \
                config_enc_param();                                           \
                EXPECT_EQ(EB_ErrorNone,                                       \
                          eb_svt_enc_set_parameter(ctxt_.enc_handle,          \
                                                   &ctxt_.enc_params))        \
                PRINT_PARAM_FATAL(ctxt_.enc_params.param_name);               \
                EncParamTestBase::TearDown();                                 \
            }                                                                 \
        }                                                                     \
        virtual void run_invalid_param_check() override {                     \
            for (size_t i = 0; i < SIZE_INVALID_PARAM(param_name); ++i) {     \
                EncParamTestBase::SetUp();                                    \
                ctxt_.enc_params.param_name =                                 \
                    GET_INVALID_PARAM(param_name, i);                         \
                config_enc_param();                                           \
                EXPECT_EQ(EB_ErrorBadParameter,                               \
                          eb_svt_enc_set_parameter(ctxt_.enc_handle,          \
                                                   &ctxt_.enc_params))        \
                PRINT_PARAM_FATAL(ctxt_.enc_params.param_name);               \
                EncParamTestBase::TearDown();                                 \
            }                                                                 \
        }                                                                     \
        virtual void run_special_param_check() override {                     \
            /*do nothing for special cases*/                                  \
        }                                                                     \
                                                                              \
      protected:                                                              \
        virtual void SetUp() override {                                       \
            /* skip EncParamTestBase::SetUp() */                              \
        }                                                                     \
        virtual void TearDown() override {                                    \
            /* skip EncParamTestBase::TearDown() */                           \
        }                                                                     \
    };

/** Test case for enc_mode*/
DEFINE_PARAM_TEST_CLASS(EncParamEncModeTest, enc_mode);
PARAM_TEST(EncParamEncModeTest);

/** Test case for intra_period_length*/
DEFINE_PARAM_TEST_CLASS(EncParamIntraPeridLenTest, intra_period_length);
PARAM_TEST(EncParamIntraPeridLenTest);

/** Test case for intra_refresh_type*/
DEFINE_PARAM_TEST_CLASS(EncParamIntraRefreshTypeTest, intra_refresh_type);
PARAM_TEST(EncParamIntraRefreshTypeTest);

/** Test case for hierarchical_levels*/
DEFINE_PARAM_TEST_CLASS(EncParamHierarchicalLvlTest, hierarchical_levels);
PARAM_TEST(EncParamHierarchicalLvlTest);

/** Test case for pred_structure*/
DEFINE_PARAM_TEST_CLASS(EncParamPredStructTest, pred_structure);
PARAM_TEST(EncParamPredStructTest);

/** Test case for base_layer_switch_mode*/
DEFINE_PARAM_TEST_CLASS(EncParamBaseLayerSwitchModeTest,
                        base_layer_switch_mode);
PARAM_TEST(EncParamBaseLayerSwitchModeTest);

/** Test case for source_width*/
DEFINE_PARAM_TEST_CLASS(EncParamSrcWidthTest, source_width);
PARAM_TEST(EncParamSrcWidthTest);

/** Test case for source_height*/
DEFINE_PARAM_TEST_CLASS(EncParamSrcHeightTest, source_height);
PARAM_TEST(EncParamSrcHeightTest);

/** Test case for frame_rate*/
DEFINE_PARAM_TEST_CLASS(EncParamFrameRateTest, frame_rate);
PARAM_TEST(EncParamFrameRateTest);

/** Test case for encoder_bit_depth*/
DEFINE_PARAM_TEST_CLASS(EncParamEncBitDepthTest, encoder_bit_depth);
PARAM_TEST(EncParamEncBitDepthTest);

/** Test case for compressed_ten_bit_format*/
DEFINE_PARAM_TEST_CLASS(EncParamCompr10BitFmtTest, compressed_ten_bit_format);
PARAM_TEST(EncParamCompr10BitFmtTest);

/** Test case for frames_to_be_encoded*/
DEFINE_PARAM_TEST_CLASS(EncParamFrame2EncTest, frames_to_be_encoded);
PARAM_TEST(EncParamFrame2EncTest);

/** Test case for sb_sz*/
DEFINE_PARAM_TEST_CLASS(EncParamSbSizeTest, sb_sz);
PARAM_TEST(EncParamSbSizeTest);

/** Test case for super_block_size*/
DEFINE_PARAM_TEST_CLASS(EncParamSuperBlockSizeTest, super_block_size);
PARAM_TEST(EncParamSuperBlockSizeTest);

/** Test case for partition_depth*/
DEFINE_PARAM_TEST_CLASS(EncParamPartitionDepthTest, partition_depth);
PARAM_TEST(EncParamPartitionDepthTest);

/** Test case for qp*/
DEFINE_PARAM_TEST_CLASS(EncParamQPTest, qp);
PARAM_TEST(EncParamQPTest);

/** Test case for use_qp_file*/
DEFINE_PARAM_TEST_CLASS(EncParamUseQPFileTest, use_qp_file);
PARAM_TEST(EncParamUseQPFileTest);

/** Test case for enable_qp_scaling_flag*/
DEFINE_PARAM_TEST_CLASS(EncParamEnableQPScaleTest, enable_qp_scaling_flag);
PARAM_TEST(EncParamEnableQPScaleTest);

/** Test case for disable_dlf_flag*/
DEFINE_PARAM_TEST_CLASS(EncParamDisableDlfTest, disable_dlf_flag);
PARAM_TEST(EncParamDisableDlfTest);

/** Test case for enable_denoise_flag*/
DEFINE_PARAM_TEST_CLASS(EncParamEnableDenoiseTest, enable_denoise_flag);
PARAM_TEST(EncParamEnableDenoiseTest);

/** Test case for film_grain_denoise_strength*/
DEFINE_PARAM_TEST_CLASS(EncParamFilmGrainDenoiseStrTest,
                        film_grain_denoise_strength);
PARAM_TEST(EncParamFilmGrainDenoiseStrTest);

/** Test case for enable_warped_motion*/
DEFINE_PARAM_TEST_CLASS(EncParamEnableWarpedMotionTest, enable_warped_motion);
PARAM_TEST(EncParamEnableWarpedMotionTest);

/** Test case for enable_global_motion*/
DEFINE_PARAM_TEST_CLASS(EncParamEnableGlobalMotionTest, enable_global_motion);
PARAM_TEST(EncParamEnableGlobalMotionTest);

/** Test case for use_default_me_hme*/
DEFINE_PARAM_TEST_CLASS(EncParamUseDefaultMeHmeTest, use_default_me_hme);
PARAM_TEST(EncParamUseDefaultMeHmeTest);

/** Test case for enable_hme_flag*/
DEFINE_PARAM_TEST_CLASS(EncParamEnableHmeTest, enable_hme_flag);
PARAM_TEST(EncParamEnableHmeTest);

/** Test case for ext_block_flag*/
DEFINE_PARAM_TEST_CLASS(EncParamExtBlockTest, ext_block_flag);
PARAM_TEST(EncParamExtBlockTest);

/** Test case for in_loop_me_flag*/
DEFINE_PARAM_TEST_CLASS(EncParamInLoopMeTest, in_loop_me_flag);
PARAM_TEST(EncParamInLoopMeTest);

/** Test case for search_area_width*/
DEFINE_PARAM_TEST_CLASS(EncParamSearchAreaWidthTest, search_area_width);
PARAM_TEST(EncParamSearchAreaWidthTest);

/** Test case for search_area_height*/
DEFINE_PARAM_TEST_CLASS(EncParamSearchAreaHeightTest, search_area_height);
PARAM_TEST(EncParamSearchAreaHeightTest);

/** Test case for enable_palette*/
DEFINE_PARAM_TEST_CLASS(EncParamEnablePaletteTest, enable_palette);
PARAM_TEST(EncParamEnablePaletteTest);

/** Test case for constrained_intra*/
DEFINE_PARAM_TEST_CLASS(EncParamConstrainedIntraTest, constrained_intra);
PARAM_TEST(EncParamConstrainedIntraTest);

/** Test case for rate_control_mode*/
DEFINE_PARAM_TEST_CLASS(EncParamRateCtrlModeTest, rate_control_mode);
PARAM_TEST(EncParamRateCtrlModeTest);

/** Test case for scene_change_detection*/
DEFINE_PARAM_TEST_CLASS(EncParamSceneChangeDectTest, scene_change_detection);
PARAM_TEST(EncParamSceneChangeDectTest);

/** Test case for look_ahead_distance*/
DEFINE_PARAM_TEST_CLASS(EncParamLookAheadDistanceTest, look_ahead_distance);
PARAM_TEST(EncParamLookAheadDistanceTest);

/** Test case for target_bit_rate*/
DEFINE_PARAM_TEST_CLASS(EncParamTargetBitRateTest, target_bit_rate);
PARAM_TEST(EncParamTargetBitRateTest);

/** Test case for max_qp_allowed*/
DEFINE_PARAM_TEST_CLASS(EncParamMaxQPAllowTest, max_qp_allowed);
PARAM_TEST(EncParamMaxQPAllowTest);

/** Test case for min_qp_allowed*/
DEFINE_PARAM_TEST_CLASS(EncParamMinQPAllowTest, min_qp_allowed);
PARAM_TEST(EncParamMinQPAllowTest);

/** Test case for high_dynamic_range_input*/
DEFINE_PARAM_TEST_CLASS(EncParamHighDynamicRangeInputTest,
                        high_dynamic_range_input);
PARAM_TEST(EncParamHighDynamicRangeInputTest);

/** Test case for profile*/
DEFINE_PARAM_TEST_CLASS(EncParamProfileTest, profile);
PARAM_TEST(EncParamProfileTest);

/** Test case for tier*/
DEFINE_PARAM_TEST_CLASS(EncParamTierTest, tier);
PARAM_TEST(EncParamTierTest);

/** Test case for level*/
DEFINE_PARAM_TEST_CLASS(EncParamLevelTest, level);
PARAM_TEST(EncParamLevelTest);

/** Test case for use_cpu_flags*/
DEFINE_PARAM_TEST_CLASS(EncParamOplLevelTest, use_cpu_flags);
PARAM_TEST(EncParamOplLevelTest);

/** Test case for channel_id*/
DEFINE_PARAM_TEST_CLASS(EncParamChIdTest, channel_id);
PARAM_TEST(EncParamChIdTest);

/** Test case for active_channel_count*/
DEFINE_PARAM_TEST_CLASS(EncParamActiveChCountTest, active_channel_count);
PARAM_TEST(EncParamActiveChCountTest);

/** Test case for speed_control_flag*/
DEFINE_PARAM_TEST_CLASS(EncParamSpeedCtrlTest, speed_control_flag);
PARAM_TEST(EncParamSpeedCtrlTest);

/** Test case for injector_frame_rate*/
DEFINE_PARAM_TEST_CLASS(EncParamInjectorFrameRateTest, injector_frame_rate);
PARAM_TEST(EncParamInjectorFrameRateTest);

/** Test case for logical_processors*/
DEFINE_PARAM_TEST_CLASS(EncParamLogicalProcessorsTest, logical_processors);
PARAM_TEST(EncParamLogicalProcessorsTest);

/** Test case for target_socket*/
DEFINE_PARAM_TEST_CLASS(EncParamTargetSocketTest, target_socket);
PARAM_TEST(EncParamTargetSocketTest);

/** Test case for recon_enabled*/
DEFINE_PARAM_TEST_CLASS(EncParamReconEnabledTest, recon_enabled);
PARAM_TEST(EncParamReconEnabledTest);

#if TILES
/** Test case for tile_columns*/
DEFINE_PARAM_TEST_CLASS(EncParamTileColsTest, tile_columns);
PARAM_TEST(EncParamTileColsTest);

/** Test case for tile_rows*/
DEFINE_PARAM_TEST_CLASS(EncParamTileRowsTest, tile_rows);
PARAM_TEST(EncParamTileRowsTest);
#endif

/** Test case for screen_content_mode*/
DEFINE_PARAM_TEST_CLASS(EncParamScreenContentModeTest, screen_content_mode);
PARAM_TEST(EncParamScreenContentModeTest);

/** Test case for enable_altrefs*/
DEFINE_PARAM_TEST_CLASS(EncParamEnableAltRefsTest, enable_altrefs);
PARAM_TEST(EncParamEnableAltRefsTest);

/** Test case for altref_strength*/
DEFINE_PARAM_TEST_CLASS(EncParamAltRefsStrengthTest, altref_strength);
PARAM_TEST(EncParamAltRefsStrengthTest);

/** Test case for altref_nframes*/
DEFINE_PARAM_TEST_CLASS(EncParamAltRefsFramesNumTest, altref_nframes);
PARAM_TEST(EncParamAltRefsFramesNumTest);

/** Test case for enable_overlays*/
DEFINE_PARAM_TEST_CLASS(EncParamEnableOverlaysTest, enable_overlays);
PARAM_TEST(EncParamEnableOverlaysTest);

}  // namespace
