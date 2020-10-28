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
 * @file OBMCVarianceTest.cc
 *
 * @brief Unit test for obmc variance functions:
 * - svt_aom_obmc_variance{4-128}x{4-128}_{c, avx2}
 * - svt_aom_obmc_sub_pixel_variance{4-128}x{4-128}_{c, sse4_1}
 *
 * @author Cidana-Edmond
 *
 ******************************************************************************/
#include "gtest/gtest.h"
// workaround to eliminate the compiling warning on linux
// The macro will conflict with definition in gtest.h
#ifdef __USE_GNU
#undef __USE_GNU  // defined in EbThreads.h
#endif
#ifdef _GNU_SOURCE
#undef _GNU_SOURCE  // defined in EbThreads.h
#endif
#include "aom_dsp_rtcd.h"
#include "random.h"
#include "util.h"
#include "EbUtility.h"
#include "filter.h"

using std::tuple;
using svt_av1_test_tool::SVTRandom;  // to generate the random

namespace {
static const int MaskMax = 64;

using ObmcVarFunc = unsigned int (*)(const uint8_t *pre, int pre_stride,
                                     const int32_t *wsrc, const int32_t *mask,
                                     unsigned int *sse);
using ObmcVarParam = tuple<ObmcVarFunc, ObmcVarFunc>;

class OBMCVarianceTest : public ::testing::TestWithParam<ObmcVarParam> {
  public:
    OBMCVarianceTest()
        : rnd_(8, false),
          rnd_msk_(0, MaskMax * MaskMax + 1),
          func_ref_(TEST_GET_PARAM(0)),
          func_tst_(TEST_GET_PARAM(1)) {
        pre_ = reinterpret_cast<uint8_t *>(svt_aom_memalign(32, MAX_SB_SQUARE));
        wsrc_buf_ = reinterpret_cast<int32_t *>(
            svt_aom_memalign(32, MAX_SB_SQUARE * sizeof(int32_t)));
        mask_buf_ = reinterpret_cast<int32_t *>(
            svt_aom_memalign(32, MAX_SB_SQUARE * sizeof(int32_t)));
    }

    ~OBMCVarianceTest() {
        if (pre_)
            svt_aom_free(pre_);
        if (wsrc_buf_)
            svt_aom_free(wsrc_buf_);
        if (mask_buf_)
            svt_aom_free(mask_buf_);
    }

  protected:
    void run_test(size_t test_num) {
        for (size_t i = 0; i < test_num; i++) {
            for (size_t j = 0; j < MAX_SB_SQUARE; j++) {
                pre_[j] = rnd_.random();
                wsrc_buf_[j] = rnd_.random() * rnd_msk_.random();
                mask_buf_[j] = rnd_msk_.random();
            }

            unsigned int sse_ref = 0, sse_tst = 0;
            uint32_t var_ref =
                func_ref_(pre_, MAX_SB_SIZE, wsrc_buf_, mask_buf_, &sse_ref);
            uint32_t var_tst =
                func_tst_(pre_, MAX_SB_SIZE, wsrc_buf_, mask_buf_, &sse_tst);

            ASSERT_EQ(var_tst, var_ref) << "compare var error";
            ASSERT_EQ(sse_tst, sse_ref) << "compare sse error";
        }
    }

  protected:
    SVTRandom rnd_;
    SVTRandom rnd_msk_;
    ObmcVarFunc func_ref_;
    ObmcVarFunc func_tst_;
    uint8_t *pre_;
    int32_t *wsrc_buf_;
    int32_t *mask_buf_;
};

TEST_P(OBMCVarianceTest, RunCheckOutput) {
    run_test(1000);
};

#define OBMC_VAR_FUNC_C(W, H) svt_aom_obmc_variance##W##x##H##_c
#define OBMC_VAR_FUNC_AVX2(W, H) svt_aom_obmc_variance##W##x##H##_avx2
#define GEN_OBMC_VAR_TEST_PARAM(W, H) \
    ObmcVarParam(OBMC_VAR_FUNC_C(W, H), OBMC_VAR_FUNC_AVX2(W, H))
#define GEN_TEST_PARAMS(GEN_PARAM)                                          \
    GEN_PARAM(128, 128), GEN_PARAM(128, 64), GEN_PARAM(64, 128),            \
        GEN_PARAM(64, 64), GEN_PARAM(64, 32), GEN_PARAM(32, 64),            \
        GEN_PARAM(32, 32), GEN_PARAM(32, 16), GEN_PARAM(16, 32),            \
        GEN_PARAM(16, 16), GEN_PARAM(16, 8), GEN_PARAM(8, 16),              \
        GEN_PARAM(8, 8), GEN_PARAM(8, 4), GEN_PARAM(4, 8), GEN_PARAM(4, 4), \
        GEN_PARAM(4, 16), GEN_PARAM(16, 4), GEN_PARAM(8, 32),               \
        GEN_PARAM(32, 8), GEN_PARAM(16, 64), GEN_PARAM(64, 16)

static const ObmcVarParam obmc_var_test_params[] = {
    GEN_TEST_PARAMS(GEN_OBMC_VAR_TEST_PARAM)};

INSTANTIATE_TEST_CASE_P(OBMC, OBMCVarianceTest,
                        ::testing::ValuesIn(obmc_var_test_params));

using ObmcSubPixVarFunc = unsigned int (*)(const uint8_t *pre, int pre_stride,
                                           int xoffset, int yoffset,
                                           const int32_t *wsrc,
                                           const int32_t *mask,
                                           unsigned int *sse);
using ObmcSubPixVarParam = tuple<ObmcSubPixVarFunc, ObmcSubPixVarFunc>;

class OBMCSubPixelVarianceTest
    : public ::testing::TestWithParam<ObmcSubPixVarParam> {
  public:
    OBMCSubPixelVarianceTest()
        : rnd_(8, false),
          rnd_msk_(0, MaskMax * MaskMax + 1),
          rnd_offset_(0, BIL_SUBPEL_SHIFTS - 1),
          func_ref_(TEST_GET_PARAM(0)),
          func_tst_(TEST_GET_PARAM(1)) {
        pre_ = reinterpret_cast<uint8_t *>(svt_aom_memalign(32, MAX_SB_SQUARE));
        wsrc_buf_ = reinterpret_cast<int32_t *>(
            svt_aom_memalign(32, MAX_SB_SQUARE * sizeof(int32_t)));
        mask_buf_ = reinterpret_cast<int32_t *>(
            svt_aom_memalign(32, MAX_SB_SQUARE * sizeof(int32_t)));
    }

    ~OBMCSubPixelVarianceTest() {
        if (pre_)
            svt_aom_free(pre_);
        if (wsrc_buf_)
            svt_aom_free(wsrc_buf_);
        if (mask_buf_)
            svt_aom_free(mask_buf_);
    }

  protected:
    void run_test(size_t test_num) {
        for (size_t i = 0; i < test_num; i++) {
            for (size_t j = 0; j < MAX_SB_SQUARE; j++) {
                pre_[j] = rnd_.random();
                wsrc_buf_[j] = rnd_.random() * rnd_msk_.random();
                mask_buf_[j] = rnd_msk_.random();
            }

            int offset_x = rnd_offset_.random();
            int offset_y = rnd_offset_.random();
            unsigned int sse_ref = 0, sse_tst = 0;
            uint32_t var_ref = func_ref_(pre_,
                                         MAX_SB_SIZE,
                                         offset_x,
                                         offset_y,
                                         wsrc_buf_,
                                         mask_buf_,
                                         &sse_ref);
            uint32_t var_tst = func_tst_(pre_,
                                         MAX_SB_SIZE,
                                         offset_x,
                                         offset_y,
                                         wsrc_buf_,
                                         mask_buf_,
                                         &sse_tst);

            ASSERT_EQ(var_tst, var_ref)
                << "compare var error at offset x=" << offset_x
                << " y=" << offset_y;
            ASSERT_EQ(sse_tst, sse_ref)
                << "compare sse error at offset x=" << offset_x
                << " y=" << offset_y;
        }
    }

  protected:
    SVTRandom rnd_;
    SVTRandom rnd_msk_;
    SVTRandom rnd_offset_;
    ObmcSubPixVarFunc func_ref_;
    ObmcSubPixVarFunc func_tst_;
    uint8_t *pre_;
    int32_t *wsrc_buf_;
    int32_t *mask_buf_;
};

TEST_P(OBMCSubPixelVarianceTest, RunCheckOutput) {
    run_test(1000);
};

#define OBMC_SUB_PIX_VAR_FUNC_C(W, H) svt_aom_obmc_sub_pixel_variance##W##x##H##_c
#define OBMC_SUB_PIX_VAR_FUNC_SSE41(W, H) \
    svt_aom_obmc_sub_pixel_variance##W##x##H##_sse4_1
#define GEN_OBMC_SUB_PIX_VAR_TEST_PARAM(W, H)         \
    ObmcSubPixVarParam(OBMC_SUB_PIX_VAR_FUNC_C(W, H), \
                       OBMC_SUB_PIX_VAR_FUNC_SSE41(W, H))

static const ObmcSubPixVarParam obmc_sub_pix_var_test_params[] = {
    GEN_TEST_PARAMS(GEN_OBMC_SUB_PIX_VAR_TEST_PARAM)};

INSTANTIATE_TEST_CASE_P(OBMC, OBMCSubPixelVarianceTest,
                        ::testing::ValuesIn(obmc_sub_pix_var_test_params));

}  // namespace
