/*
 * Copyright(c) 2019 Intel Corporation
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */


 /******************************************************************************
  * @EbNNPredictTests.cc
  *
  * @brief Unit test for avi_nn_predict c vs sse3
  * @author JackZhouVCD
  *
  ******************************************************************************/

#include "gtest/gtest.h"
#include "EbDefinitions.h"
#include "EbUnitTest.h"
#include "EbUnitTestUtility.h"
#include <immintrin.h>

#include "ml.h"
#include "partition_model_weights.h"

#include "aom_dsp_rtcd.h"

#define FEATURE_SIZE_MAX_MIN_PART_PRED 13
#define MAX_NUM_CLASSES_MAX_MIN_PART_PRED 4

using namespace std;

namespace {

    typedef void (*av1_nn_predict_func)(
        const float *input_nodes,
        const NN_CONFIG *const nn_config,
        int reduce_prec,
        float *const output);

    class NNPredTest : public ::testing::TestWithParam<av1_nn_predict_func> {

    public:
        NNPredTest() : func_(GetParam()){}
        ~NNPredTest();

        void SetUp() override {
        }

        void TearDown() override {
        }
    protected:

        void run_test();

        void init_data() {

            srand(unsigned(time(NULL)));

            nn_config = &av1_max_part_pred_nn_config;
            memset(features1, 0, FEATURE_SIZE_MAX_MIN_PART_PRED * sizeof(float));
            memset(features2, 1, FEATURE_SIZE_MAX_MIN_PART_PRED*sizeof(float));

            for (int i = 0; i < FEATURE_SIZE_MAX_MIN_PART_PRED; ++i)
                features3[i] = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / 100));
        }

        void check_output(float* ref, float* test) {

            for (int i= 0; i < MAX_NUM_CLASSES_MAX_MIN_PART_PRED; i++)
                ASSERT_EQ(ref[i],test[i])
                    << " Mismatch at index " << i
                    << " C value: " << ref[i]
                    << " SSE3 value: " << test[i];
        }

        av1_nn_predict_func func_;

        const NN_CONFIG *nn_config;

        float   features1[FEATURE_SIZE_MAX_MIN_PART_PRED],
                features2[FEATURE_SIZE_MAX_MIN_PART_PRED],
                features3[FEATURE_SIZE_MAX_MIN_PART_PRED],
                ref_output[MAX_NUM_CLASSES_MAX_MIN_PART_PRED],
                test_output[MAX_NUM_CLASSES_MAX_MIN_PART_PRED];
    };


    NNPredTest::~NNPredTest(){}

    void NNPredTest::run_test() {
        init_data();

        // ALL zeroes
        av1_nn_predict_c(features1, nn_config, 1, ref_output);
        func_(features1,nn_config,1,test_output);
        check_output(ref_output, test_output);

        // ALL ones
        av1_nn_predict_c(features2, nn_config, 1, ref_output);
        func_(features2, nn_config, 1, test_output);
        check_output(ref_output, test_output);

        // Random input
        av1_nn_predict_c(features3, nn_config, 1, ref_output);
        func_(features3, nn_config, 1, test_output);
        check_output(ref_output, test_output);

    }

    TEST_P(NNPredTest, NNPredTest) {
        run_test();
    };

    INSTANTIATE_TEST_CASE_P(NNPRED, NNPredTest,
                            ::testing::Values(av1_nn_predict_sse3));

}  // namespace
