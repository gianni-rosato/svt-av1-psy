/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file CompoundUtilTest.cc
 *
 * @brief Unit test for util functions in wedge prediction:
 * - av1_build_compound_diffwtd_mask_{, d16}avx2
 * - aom_blend_a64_mask_avx2/aom_lowbd_blend_a64_d16_mask_avx2
 *
 * @author Cidana-Wenyao
 *
 ******************************************************************************/
#include "gtest/gtest.h"
#include "EbDefinitions.h"
#include "EbUtility.h"
#include "aom_dsp_rtcd.h"
#include "random.h"
#include "convolve.h"
#include "util.h"

using svt_av1_test_tool::SVTRandom;
namespace {
extern "C" void aom_blend_a64_mask_c(uint8_t *dst, uint32_t dst_stride,
                                     const uint8_t *src0, uint32_t src0_stride,
                                     const uint8_t *src1, uint32_t src1_stride,
                                     const uint8_t *mask, uint32_t mask_stride,
                                     int w, int h, int subw, int subh);

template <typename SrcSample, typename DstSample>
class CompBlendTest : public ::testing::Test {
  public:
    void SetUp() override {
        ref_dst_ =
            (DstSample *)eb_aom_memalign(32, MAX_SB_SQUARE * sizeof(DstSample));
        tst_dst_ =
            (DstSample *)eb_aom_memalign(32, MAX_SB_SQUARE * sizeof(DstSample));
        src0_ =
            (SrcSample *)eb_aom_memalign(32, MAX_SB_SQUARE * sizeof(SrcSample));
        src1_ =
            (SrcSample *)eb_aom_memalign(32, MAX_SB_SQUARE * sizeof(SrcSample));
        mask_ = (uint8_t *)eb_aom_memalign(32, MAX_SB_SQUARE);
    }

    void TearDown() override {
        eb_aom_free(ref_dst_);
        eb_aom_free(tst_dst_);
        eb_aom_free(src0_);
        eb_aom_free(src1_);
        eb_aom_free(mask_);
        aom_clear_system_state();
    }

    void run_test() {
        const int loops = 10000;
        SVTRandom rnd(0, (1 << bd_) - 1);
        SVTRandom mask_rnd(0, 64);

        // generate random mask
        for (int i = 0; i < MAX_SB_SQUARE; ++i)
            mask_[i] = mask_rnd.random();

        for (int k = 0; k < loops; ++k) {
            for (int block_size = BLOCK_4X4; block_size < BlockSizeS_ALL;
                 block_size += 1) {
                w_ = block_size_wide[block_size];
                h_ = block_size_high[block_size];

                // initial the input data
                for (int i = 0; i < h_; ++i) {
                    for (int j = 0; j < w_; ++j) {
                        src0_[i * src_stride_ + j] = rnd.random();
                        src1_[i * src_stride_ + j] = rnd.random();
                    }
                }

                for (int subh = 0; subh < 2; ++subh) {
                    for (int subw = 0; subw < 2; ++subw) {
                        memset(ref_dst_, 0, MAX_SB_SQUARE * sizeof(DstSample));
                        memset(tst_dst_, 0, MAX_SB_SQUARE * sizeof(DstSample));

                        run_blend(subw, subh);

                        // check output
                        for (int i = 0; i < h_; ++i) {
                            for (int j = 0; j < w_; ++j) {
                                ASSERT_EQ(tst_dst_[i * dst_stride_ + j],
                                          ref_dst_[i * dst_stride_ + j]);
                            }
                        }
                    }
                }
            }
        }
    }

    virtual void run_blend(int subw, int subh) = 0;

  protected:
    static const int src_stride_ = MAX_SB_SIZE;
    static const int dst_stride_ = MAX_SB_SIZE;
    static const int mask_stride_ = MAX_SB_SIZE;
    DstSample *ref_dst_;
    DstSample *tst_dst_;
    SrcSample *src0_;
    SrcSample *src1_;
    uint8_t *mask_;
    int w_, h_;
    int bd_;
};

class LbdCompBlendTest : public CompBlendTest<uint8_t, uint8_t> {
  public:
    LbdCompBlendTest() {
        bd_ = 8;
    }

    void run_blend(int subw, int subh) override {
        aom_blend_a64_mask_c(ref_dst_,
                             dst_stride_,
                             src0_,
                             src_stride_,
                             src1_,
                             src_stride_,
                             mask_,
                             mask_stride_,
                             w_,
                             h_,
                             subw,
                             subh);
        aom_blend_a64_mask_avx2(tst_dst_,
                                dst_stride_,
                                src0_,
                                src_stride_,
                                src1_,
                                src_stride_,
                                mask_,
                                mask_stride_,
                                w_,
                                h_,
                                subw,
                                subh);
    }
};

TEST_F(LbdCompBlendTest, BlendA64Mask) {
    run_test();
}

class LbdCompBlendD16Test : public CompBlendTest<uint16_t, uint8_t> {
  public:
    LbdCompBlendD16Test() {
        bd_ = 10;
    }

    void run_blend(int subw, int subh) override {
        ConvolveParams conv_params;
        conv_params.round_0 = ROUND0_BITS;
        conv_params.round_1 = COMPOUND_ROUND1_BITS;
        aom_lowbd_blend_a64_d16_mask_avx2(tst_dst_,
                                          dst_stride_,
                                          src0_,
                                          src_stride_,
                                          src1_,
                                          src_stride_,
                                          mask_,
                                          mask_stride_,
                                          w_,
                                          h_,
                                          subw,
                                          subh,
                                          &conv_params);
        aom_lowbd_blend_a64_d16_mask_c(ref_dst_,
                                       dst_stride_,
                                       src0_,
                                       src_stride_,
                                       src1_,
                                       src_stride_,
                                       mask_,
                                       mask_stride_,
                                       w_,
                                       h_,
                                       subw,
                                       subh,
                                       &conv_params);
    }
};

TEST_F(LbdCompBlendD16Test, BlendA64MaskD16) {
    run_test();
}

typedef void (*BuildCompDiffwtdMaskedFunc)(uint8_t *mask,
                                           DIFFWTD_MASK_TYPE mask_type,
                                           const uint8_t *src0, int src0_stride,
                                           const uint8_t *src1, int src1_stride,
                                           int h, int w);
typedef ::testing::tuple<BlockSize, BuildCompDiffwtdMaskedFunc>
    BuildCompDiffwtdMaskParam;

class BuildCompDiffwtdMaskTest
    : public ::testing::TestWithParam<BuildCompDiffwtdMaskParam> {
  public:
    BuildCompDiffwtdMaskTest() : rnd_(0, 255){};
    virtual ~BuildCompDiffwtdMaskTest() {
    }

    void TearDown() override {
        aom_clear_system_state();
    }

    void run_test(const DIFFWTD_MASK_TYPE type) {
        const int block_size = TEST_GET_PARAM(0);
        BuildCompDiffwtdMaskedFunc test_impl = TEST_GET_PARAM(1);
        const int width = block_size_wide[block_size];
        const int height = block_size_high[block_size];
        DECLARE_ALIGNED(16, uint8_t, mask_ref[MAX_SB_SQUARE]);
        DECLARE_ALIGNED(16, uint8_t, mask_test[MAX_SB_SQUARE]);
        DECLARE_ALIGNED(16, uint8_t, src0[MAX_SB_SQUARE]);
        DECLARE_ALIGNED(16, uint8_t, src1[MAX_SB_SQUARE]);
        const int run_times = 100;
        for (int i = 0; i < run_times; ++i) {
            for (int i = 0; i < width * height; i++) {
                src0[i] = rnd_.random();
                src1[i] = rnd_.random();
            }

            av1_build_compound_diffwtd_mask_c(
                mask_ref, type, src0, width, src1, width, height, width);

            test_impl(mask_test, type, src0, width, src1, width, height, width);
            for (int r = 0; r < height; ++r) {
                for (int c = 0; c < width; ++c) {
                    ASSERT_EQ(mask_ref[c + r * width], mask_test[c + r * width])
                        << "[" << r << "," << c << "] " << i << " @ " << width
                        << "x" << height << " inv " << type;
                }
            }
        }
    }

  private:
    SVTRandom rnd_;
};

TEST_P(BuildCompDiffwtdMaskTest, MatchTest) {
    run_test(DIFFWTD_38);
    run_test(DIFFWTD_38_INV);
}

INSTANTIATE_TEST_CASE_P(
    CompUtilTest, BuildCompDiffwtdMaskTest,
    ::testing::Combine(
        ::testing::Range(BLOCK_4X4, BlockSizeS_ALL),
        ::testing::Values(av1_build_compound_diffwtd_mask_avx2)));

// test av1_build_compound_diffwtd_mask_d16_avx2
typedef void (*BuildCompDiffwtdMaskD16Func)(
    uint8_t *mask, DIFFWTD_MASK_TYPE mask_type, const CONV_BUF_TYPE *src0,
    int src0_stride, const CONV_BUF_TYPE *src1, int src1_stride, int h, int w,
    ConvolveParams *conv_params, int bd);

typedef ::testing::tuple<int, BuildCompDiffwtdMaskD16Func, BlockSize>
    BuildCompDiffwtdMaskD16Param;

class BuildCompDiffwtdMaskD16Test
    : public ::testing::TestWithParam<BuildCompDiffwtdMaskD16Param> {
  public:
    BuildCompDiffwtdMaskD16Test() : rnd_(16, false) {
    }

    ~BuildCompDiffwtdMaskD16Test() {
    }

    void TearDown() override {
        aom_clear_system_state();
    }

  protected:
    void run_test() {
        BuildCompDiffwtdMaskD16Func tst_func = TEST_GET_PARAM(1);
        const int block_size = TEST_GET_PARAM(2);
        const int bd = TEST_GET_PARAM(0);
        const int width = block_size_wide[block_size];
        const int height = block_size_high[block_size];
        DECLARE_ALIGNED(16, uint8_t, mask_ref[2 * MAX_SB_SQUARE]);
        DECLARE_ALIGNED(16, uint8_t, mask_test[2 * MAX_SB_SQUARE]);
        DECLARE_ALIGNED(32, uint16_t, src0[MAX_SB_SQUARE]);
        DECLARE_ALIGNED(32, uint16_t, src1[MAX_SB_SQUARE]);

        ConvolveParams conv_params =
            get_conv_params_no_round(0 /*unused*/, 0, 0, NULL, 0, 1, bd);

        int in_precision = bd + 2 * FILTER_BITS - conv_params.round_0 -
                           conv_params.round_1 + 2;

        for (int i = 0; i < MAX_SB_SQUARE; i++) {
            src0[i] = rnd_.random() & ((1 << in_precision) - 1);
            src1[i] = rnd_.random() & ((1 << in_precision) - 1);
        }

        for (int mask_type = 0; mask_type < DIFFWTD_MASK_TYPES; mask_type++) {
            av1_build_compound_diffwtd_mask_d16_c(mask_ref,
                                                  (DIFFWTD_MASK_TYPE)mask_type,
                                                  src0,
                                                  width,
                                                  src1,
                                                  width,
                                                  height,
                                                  width,
                                                  &conv_params,
                                                  bd);

            tst_func(mask_test,
                     (DIFFWTD_MASK_TYPE)mask_type,
                     src0,
                     width,
                     src1,
                     width,
                     height,
                     width,
                     &conv_params,
                     bd);

            for (int r = 0; r < height; ++r) {
                for (int c = 0; c < width; ++c) {
                    ASSERT_EQ(mask_ref[c + r * width], mask_test[c + r * width])
                        << "Mismatch at unit tests for "
                           "BuildCompDiffwtdMaskD16Test\n"
                        << " Pixel mismatch at index "
                        << "[" << r << "," << c << "] "
                        << " @ " << width << "x" << height << " inv "
                        << mask_type;
                }
            }
        }
    }
    SVTRandom rnd_;
};  // class BuildCompDiffwtdMaskD16Test

TEST_P(BuildCompDiffwtdMaskD16Test, MatchTest) {
    run_test();
}

INSTANTIATE_TEST_CASE_P(
    SSE4_1, BuildCompDiffwtdMaskD16Test,
    ::testing::Combine(
        ::testing::Range(8, 13, 2),
        ::testing::Values(av1_build_compound_diffwtd_mask_d16_avx2),
        ::testing::Range(BLOCK_4X4, BlockSizeS_ALL)));
}  // namespace
