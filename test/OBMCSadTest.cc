/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file OBMCsad_Test.cc
 *
 * @brief Unit test for obmc sad functions:
 * - obmc_sad_w4_avx2
 * - obmc_sad_w8n_avx2
 *
 * @author Cidana-Ivy
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

using std::tuple;
using svt_av1_test_tool::SVTRandom;  // to generate the random

namespace {
static const int MaskMax = 64;

using Obmcsad_Func = uint32_t (*)(const uint8_t* pre, int pre_stride,
                                 const int32_t* wsrc, const int32_t* mask);
using Obmcsad_Param = tuple<Obmcsad_Func, Obmcsad_Func>;

class OBMCsad_Test : public ::testing::TestWithParam<Obmcsad_Param> {
  public:
    OBMCsad_Test()
        : rnd_(8, false),
          rnd_msk_(0, MaskMax * MaskMax + 1),
          func_ref_(TEST_GET_PARAM(0)),
          func_tst_(TEST_GET_PARAM(1)) {
        pre_ = reinterpret_cast<uint8_t*>(eb_aom_memalign(32, MAX_SB_SQUARE));
        wsrc_buf_ = reinterpret_cast<int32_t*>(
            eb_aom_memalign(32, MAX_SB_SQUARE * sizeof(int32_t)));
        mask_buf_ = reinterpret_cast<int32_t*>(
            eb_aom_memalign(32, MAX_SB_SQUARE * sizeof(int32_t)));
    }

    ~OBMCsad_Test() {
        if (pre_)
            eb_aom_free(pre_);
        if (wsrc_buf_)
            eb_aom_free(wsrc_buf_);
        if (mask_buf_)
            eb_aom_free(mask_buf_);
    }

  protected:
    void run_test(size_t test_num) {
        for (size_t i = 0; i < test_num; i++) {
            for (size_t j = 0; j < MAX_SB_SQUARE; j++) {
                pre_[j] = rnd_.random();
                wsrc_buf_[j] = rnd_.random() * rnd_msk_.random();
                mask_buf_[j] = rnd_msk_.random();
            }

            uint32_t sad_ref =
                func_ref_(pre_, MAX_SB_SIZE, wsrc_buf_, mask_buf_);
            uint32_t sad_tst =
                func_tst_(pre_, MAX_SB_SIZE, wsrc_buf_, mask_buf_);

            ASSERT_EQ(sad_tst, sad_ref) << "compare SAD error";
        }
    }

  protected:
    SVTRandom rnd_;
    SVTRandom rnd_msk_;
    Obmcsad_Func func_ref_;
    Obmcsad_Func func_tst_;
    uint8_t* pre_;
    int32_t* wsrc_buf_;
    int32_t* mask_buf_;
};

TEST_P(OBMCsad_Test, RunCheckOutput) {
    run_test(1000);
};

#define OBMC_SAD_FUNC_C(W, H) aom_obmc_sad##W##x##H##_c
#define OBMC_SAD_FUNC_AVX2(W, H) aom_obmc_sad##W##x##H##_avx2
#define GEN_OBMC_SAD_TEST_PARAM(W, H) \
    Obmcsad_Param(OBMC_SAD_FUNC_C(W, H), OBMC_SAD_FUNC_AVX2(W, H))
#define GEN_TEST_PARAMS(GEN_PARAM)                                          \
    GEN_PARAM(128, 128), GEN_PARAM(128, 64), GEN_PARAM(64, 128),            \
        GEN_PARAM(64, 64), GEN_PARAM(64, 32), GEN_PARAM(32, 64),            \
        GEN_PARAM(32, 32), GEN_PARAM(32, 16), GEN_PARAM(16, 32),            \
        GEN_PARAM(16, 16), GEN_PARAM(16, 8), GEN_PARAM(8, 16),              \
        GEN_PARAM(8, 8), GEN_PARAM(8, 4), GEN_PARAM(4, 8), GEN_PARAM(4, 4), \
        GEN_PARAM(4, 16), GEN_PARAM(16, 4), GEN_PARAM(8, 32),               \
        GEN_PARAM(32, 8), GEN_PARAM(16, 64), GEN_PARAM(64, 16)

static const Obmcsad_Param obmc_sad_test_params[] = {
    GEN_TEST_PARAMS(GEN_OBMC_SAD_TEST_PARAM)};

INSTANTIATE_TEST_CASE_P(OBMC, OBMCsad_Test,
                        ::testing::ValuesIn(obmc_sad_test_params));

}  // namespace
