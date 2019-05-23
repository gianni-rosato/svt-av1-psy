/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file SvtAv1E2EParamsTest.cc
 *
 * @brief Impelmentation of encoder parameter coverage test in E2E test
 *
 * @author Cidana-Edmond
 *
 ******************************************************************************/

#include "gtest/gtest.h"
#include "SvtAv1E2EFramework.h"
#include "../api_test/params.h"

/**
 * @brief SVT-AV1 encoder parameter coverage E2E test
 *
 * Test strategy:
 * Setup SVT-AV1 encoder with individual parameter in vaild value and run the
 * conformance test progress. Check whether the result can match the output of
 * refence decoder.
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

/** Marcro defininition of batch processing check for default, valid, invalid
 * and special parameter check*/
#define PARAM_TEST(param_test)                          \
    TEST_P(param_test, run_paramter_conformance_test) { \
        run_conformance_test();                         \
    }                                                   \
    INSTANTIATE_TEST_CASE_P(                            \
        SVT_AV1,                                        \
        param_test,                                     \
        ::testing::ValuesIn(generate_vector_from_config("smoking_test.cfg")));

/** @breif This class is a template based on EncParamTestBase to test each
 * parameter
 */
#define DEFINE_PARAM_TEST_CLASS(test_name, param_name)                      \
    class test_name : public SvtAv1E2ETestFramework {                       \
      public:                                                               \
        test_name() {                                                       \
            param_name_str_ = #param_name;                                  \
            param_value_idx_ = 0;                                           \
        }                                                                   \
        /** initialization for test */                                      \
        void init_test() override {                                         \
            collect_ = new PerformanceCollect(typeid(this).name());         \
            av1enc_ctx_.enc_params.param_name =                             \
                GET_VALID_PARAM(param_name, param_value_idx_);              \
            /** create recon sink before setup parameter of encoder */      \
            VideoFrameParam param;                                          \
            memset(&param, 0, sizeof(param));                               \
            param.format = video_src_->get_image_format();                  \
            param.width = video_src_->get_width_with_padding();             \
            param.height = video_src_->get_height_with_padding();           \
            recon_queue_ = create_frame_queue(param);                       \
            ASSERT_NE(recon_queue_, nullptr)                                \
                << "can not create recon queue!!";                          \
            if (recon_queue_)                                               \
                av1enc_ctx_.enc_params.recon_enabled = 1;                   \
                                                                            \
            /** create reference decoder*/                                  \
            refer_dec_ = create_reference_decoder();                        \
            ASSERT_NE(refer_dec_, nullptr)                                  \
                << "can not create reference decoder!!";                    \
                                                                            \
            SvtAv1E2ETestFramework::init_test();                            \
        }                                                                   \
        /** close for test */                                               \
        void close_test() override {                                        \
            SvtAv1E2ETestFramework::close_test();                           \
            if (collect_) {                                                 \
                delete collect_;                                            \
                collect_ = nullptr;                                         \
            }                                                               \
        }                                                                   \
        /** run for the conformance test */                                 \
        void run_conformance_test() {                                       \
            for (param_value_idx_ = 0;                                      \
                 param_value_idx_ < SIZE_VALID_PARAM(param_name);           \
                 ++param_value_idx_) {                                      \
                SvtAv1E2ETestFramework::SetUp();                            \
                run_encode_process();                                       \
                SvtAv1E2ETestFramework::TearDown();                         \
            }                                                               \
        }                                                                   \
                                                                            \
      protected:                                                            \
        void SetUp() override {                                             \
            /* skip SvtAv1E2ETestFramework::SetUp() */                      \
        }                                                                   \
        void TearDown() override {                                          \
            /* skip SvtAv1E2ETestFramework::TearDown() */                   \
        }                                                                   \
                                                                            \
      protected:                                                            \
        std::string param_name_str_; /**< name of parameter for test */     \
        size_t param_value_idx_; /**< index of parameter value in vector */ \
    };

/** Test case for enc_mode*/
DEFINE_PARAM_TEST_CLASS(SvtAv1E2EParamEncModeTest, enc_mode);
PARAM_TEST(SvtAv1E2EParamEncModeTest);
