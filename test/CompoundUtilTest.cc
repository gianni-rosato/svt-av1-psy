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
 * - aom_blend_a64_mask_sse4_1
 * - aom_highbd_blend_a64_mask_sse4_1/aom_highbd_blend_a64_d16_mask_avx2
 * - aom_blend_a64_hmask_sse4_1/aom_blend_a64_vmask_sse4_1
 * - aom_highbd_blend_a64_hmask_sse4_1/aom_highbd_blend_a64_vmask_sse4_1
 *
 * @author Cidana-Wenyao
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
#include "EbDefinitions.h"
#include "EbUtility.h"
#include "aom_dsp_rtcd.h"
#include "random.h"
#include "convolve.h"
#include "util.h"

using std::make_tuple;
using svt_av1_test_tool::SVTRandom;

#define MAX_MASK_SQUARE (4 * MAX_SB_SQUARE)
#define MAKE_PARAM(func_type) std::tuple<func_type, func_type, const char *>

namespace {
template <typename SrcSample, typename DstSample, typename BlendFunc,
          typename BlendTestParam>
class CompBlendTest : public ::testing::TestWithParam<BlendTestParam> {
  public:
    CompBlendTest() {
        tst_fn_name = "";
        func_ref_ = nullptr;
        func_tst_ = nullptr;
        no_sub_ = false;
    }

    void SetUp() override {
        ref_dst_ =
            (DstSample *)eb_aom_memalign(32, MAX_SB_SQUARE * sizeof(DstSample));
        tst_dst_ =
            (DstSample *)eb_aom_memalign(32, MAX_SB_SQUARE * sizeof(DstSample));
        src0_ =
            (SrcSample *)eb_aom_memalign(32, MAX_SB_SQUARE * sizeof(SrcSample));
        src1_ =
            (SrcSample *)eb_aom_memalign(32, MAX_SB_SQUARE * sizeof(SrcSample));
        mask_ = (uint8_t *)eb_aom_memalign(32, MAX_MASK_SQUARE);
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
        const int iterations = 1000;
        SVTRandom rnd(0, (1 << bd_) - 1);
        SVTRandom mask_rnd(0, 64);

        // generate random mask
        for (int i = 0; i < MAX_MASK_SQUARE; ++i)
            mask_[i] = mask_rnd.random();

        for (int k = 0; k < iterations; ++k) {
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

                int submax = 2;
                if (no_sub_)
                    submax = 1;
                for (int subh = 0; subh < submax; ++subh) {
                    for (int subw = 0; subw < submax; ++subw) {
                        memset(ref_dst_, 0, MAX_SB_SQUARE * sizeof(DstSample));
                        memset(tst_dst_, 0, MAX_SB_SQUARE * sizeof(DstSample));

                        run_blend(subw, subh);

                        // check output
                        for (int i = 0; i < h_; ++i) {
                            for (int j = 0; j < w_; ++j) {
                                ASSERT_EQ(tst_dst_[i * dst_stride_ + j],
                                          ref_dst_[i * dst_stride_ + j])
                                    << "Mismatch at unit tests for "
                                    << tst_fn_name << "\n"
                                    << " Pixel mismatch at index "
                                    << "[" << j << "x" << i
                                    << "]\nblock size: " << block_size
                                    << "\nbit-depth: " << bd_
                                    << "\niterator: " << k;
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
    static const int mask_stride_ = 2 * MAX_SB_SIZE;
    DstSample *ref_dst_;
    DstSample *tst_dst_;
    SrcSample *src0_;
    SrcSample *src1_;
    BlendFunc func_ref_;
    BlendFunc func_tst_;
    uint8_t *mask_;
    int w_, h_;
    int bd_;
    const char *tst_fn_name;  // test function name
    bool no_sub_;
};

using LbdBlendA64MaskFunc = void (*)(uint8_t *, uint32_t, const uint8_t *,
                                     uint32_t, const uint8_t *, uint32_t,
                                     const uint8_t *, uint32_t, int, int, int,
                                     int);

class LbdCompBlendTest
    : public CompBlendTest<uint8_t, uint8_t, LbdBlendA64MaskFunc,
                           MAKE_PARAM(LbdBlendA64MaskFunc)> {
  public:
    LbdCompBlendTest() {
        bd_ = 8;
        func_ref_ = TEST_GET_PARAM(0);
        func_tst_ = TEST_GET_PARAM(1);
        tst_fn_name = TEST_GET_PARAM(2);
    }

    void run_blend(int subw, int subh) override {
        func_ref_(ref_dst_,
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
        func_tst_(tst_dst_,
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

TEST_P(LbdCompBlendTest, BlendA64Mask) {
    run_test();
}

INSTANTIATE_TEST_CASE_P(
    BLEND, LbdCompBlendTest,
    ::testing::ValuesIn({
        make_tuple(aom_blend_a64_mask_c, aom_blend_a64_mask_avx2,
                   "aom_blend_a64_mask_avx2"),
        make_tuple(aom_blend_a64_mask_c, aom_blend_a64_mask_sse4_1,
                   "aom_blend_a64_mask_sse4_1"),
    }));

using LbdBlendA64D16MaskFunc = void (*)(uint8_t *, uint32_t,
                                        const CONV_BUF_TYPE *, uint32_t,
                                        const CONV_BUF_TYPE *, uint32_t,
                                        const uint8_t *, uint32_t, int, int,
                                        int, int, ConvolveParams *);

class LbdCompBlendD16Test
    : public CompBlendTest<uint16_t, uint8_t, LbdBlendA64D16MaskFunc,
                           MAKE_PARAM(LbdBlendA64D16MaskFunc)> {
  public:
    LbdCompBlendD16Test() {
        bd_ = 10;
        func_ref_ = TEST_GET_PARAM(0);
        func_tst_ = TEST_GET_PARAM(1);
        tst_fn_name = TEST_GET_PARAM(2);
    }

    void run_blend(int subw, int subh) override {
        ConvolveParams conv_params;
        conv_params.round_0 = ROUND0_BITS;
        conv_params.round_1 = COMPOUND_ROUND1_BITS;
        func_tst_(tst_dst_,
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
        func_ref_(ref_dst_,
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

TEST_P(LbdCompBlendD16Test, BlendA64MaskD16) {
    run_test();
}

INSTANTIATE_TEST_CASE_P(
    BLEND, LbdCompBlendD16Test,
    ::testing::ValuesIn({make_tuple(aom_lowbd_blend_a64_d16_mask_c,
                                    aom_lowbd_blend_a64_d16_mask_avx2,
                                    "aom_lowbd_blend_a64_d16_mask_avx2")}));

using LbdBlendA64HMaskFunc = void (*)(uint8_t *, uint32_t, const uint8_t *,
                                      uint32_t, const uint8_t *, uint32_t,
                                      const uint8_t *, int, int);

class LbdCompBlendHMaskTest
    : public CompBlendTest<uint8_t, uint8_t, LbdBlendA64HMaskFunc,
                           MAKE_PARAM(LbdBlendA64HMaskFunc)> {
  public:
    LbdCompBlendHMaskTest() {
        bd_ = 8;
        func_ref_ = TEST_GET_PARAM(0);
        func_tst_ = TEST_GET_PARAM(1);
        tst_fn_name = TEST_GET_PARAM(2);
        no_sub_ = true;
    }

    void run_blend(int subw, int subh) override {
        (void)subw;
        (void)subh;
        func_ref_(ref_dst_,
                  dst_stride_,
                  src0_,
                  src_stride_,
                  src1_,
                  src_stride_,
                  mask_,
                  w_,
                  h_);
        func_tst_(tst_dst_,
                  dst_stride_,
                  src0_,
                  src_stride_,
                  src1_,
                  src_stride_,
                  mask_,
                  w_,
                  h_);
    }
};

TEST_P(LbdCompBlendHMaskTest, BlendA64Mask) {
    run_test();
}

INSTANTIATE_TEST_CASE_P(BLEND, LbdCompBlendHMaskTest,
                        ::testing::ValuesIn({make_tuple(
                            aom_blend_a64_hmask_c, aom_blend_a64_hmask_sse4_1,
                            "aom_blend_a64_hmask_sse4_1")}));

using LbdBlendA64VMaskFunc = void (*)(uint8_t *, uint32_t, const uint8_t *,
                                      uint32_t, const uint8_t *, uint32_t,
                                      const uint8_t *, int, int);

class LbdCompBlendVMaskTest
    : public CompBlendTest<uint8_t, uint8_t, LbdBlendA64VMaskFunc,
                           MAKE_PARAM(LbdBlendA64VMaskFunc)> {
  public:
    LbdCompBlendVMaskTest() {
        bd_ = 8;
        func_ref_ = TEST_GET_PARAM(0);
        func_tst_ = TEST_GET_PARAM(1);
        tst_fn_name = TEST_GET_PARAM(2);
        no_sub_ = true;
    }

    void run_blend(int subw, int subh) override {
        (void)subw;
        (void)subh;
        func_ref_(ref_dst_,
                  dst_stride_,
                  src0_,
                  src_stride_,
                  src1_,
                  src_stride_,
                  mask_,
                  w_,
                  h_);
        func_tst_(tst_dst_,
                  dst_stride_,
                  src0_,
                  src_stride_,
                  src1_,
                  src_stride_,
                  mask_,
                  w_,
                  h_);
    }
};

TEST_P(LbdCompBlendVMaskTest, BlendA64Mask) {
    run_test();
}

INSTANTIATE_TEST_CASE_P(BLEND, LbdCompBlendVMaskTest,
                        ::testing::ValuesIn({make_tuple(
                            aom_blend_a64_vmask_c, aom_blend_a64_vmask_sse4_1,
                            "aom_blend_a64_vmask_sse4_1")}));

using HbdBlendA64MaskFunc = void (*)(uint8_t *, uint32_t, const uint8_t *,
                                     uint32_t, const uint8_t *, uint32_t,
                                     const uint8_t *, uint32_t, int, int, int,
                                     int, int);

class HbdCompBlendTest
    : public CompBlendTest<uint16_t, uint16_t, HbdBlendA64MaskFunc,
                           MAKE_PARAM(HbdBlendA64MaskFunc)> {
  public:
    HbdCompBlendTest() {
        bd_ = 8;
        func_ref_ = TEST_GET_PARAM(0);
        func_tst_ = TEST_GET_PARAM(1);
        tst_fn_name = TEST_GET_PARAM(2);
    }

    void run_hbd_test(uint8_t bd) {
        bd_ = bd;
        run_test();
    }

    void run_blend(int subw, int subh) override {
        func_ref_((uint8_t *)ref_dst_,
                  dst_stride_,
                  (uint8_t *)src0_,
                  src_stride_,
                  (uint8_t *)src1_,
                  src_stride_,
                  mask_,
                  mask_stride_,
                  w_,
                  h_,
                  subw,
                  subh,
                  bd_);
        func_tst_((uint8_t *)tst_dst_,
                  dst_stride_,
                  (uint8_t *)src0_,
                  src_stride_,
                  (uint8_t *)src1_,
                  src_stride_,
                  mask_,
                  mask_stride_,
                  w_,
                  h_,
                  subw,
                  subh,
                  bd_);
    }
};

TEST_P(HbdCompBlendTest, BlendA64Mask) {
    run_hbd_test(8);
    run_hbd_test(10);
    run_hbd_test(12);
}

INSTANTIATE_TEST_CASE_P(
    BLEND, HbdCompBlendTest,
    ::testing::ValuesIn({make_tuple(aom_highbd_blend_a64_mask_c,
                                    aom_highbd_blend_a64_mask_sse4_1,
                                    "aom_highbd_blend_a64_mask_sse4_1")}));

using HbdBlendA64D16MaskFunc = void (*)(uint8_t *, uint32_t,
                                        const CONV_BUF_TYPE *, uint32_t,
                                        const CONV_BUF_TYPE *, uint32_t,
                                        const uint8_t *, uint32_t, int, int,
                                        int, int, ConvolveParams *, const int);

class HbdCompBlendD16Test
    : public CompBlendTest<uint16_t, uint16_t, HbdBlendA64D16MaskFunc,
                           MAKE_PARAM(HbdBlendA64D16MaskFunc)> {
  public:
    HbdCompBlendD16Test() {
        bd_ = 10;
        func_ref_ = TEST_GET_PARAM(0);
        func_tst_ = TEST_GET_PARAM(1);
        tst_fn_name = TEST_GET_PARAM(2);
    }

    void run_hbd_test(uint8_t bd) {
        bd_ = bd;
        run_test();
    }

    void run_blend(int subw, int subh) override {
        ConvolveParams conv_params;
        conv_params.round_0 = ROUND0_BITS;
        conv_params.round_1 = COMPOUND_ROUND1_BITS;
        func_tst_((uint8_t *)tst_dst_,
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
                  &conv_params,
                  bd_);
        func_ref_((uint8_t *)ref_dst_,
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
                  &conv_params,
                  bd_);
    }
};

TEST_P(HbdCompBlendD16Test, BlendA64MaskD16) {
    run_hbd_test(8);
    run_hbd_test(10);
    run_hbd_test(12);
}

INSTANTIATE_TEST_CASE_P(
    BLEND, HbdCompBlendD16Test,
    ::testing::ValuesIn({make_tuple(aom_highbd_blend_a64_d16_mask_c,
                                    aom_highbd_blend_a64_d16_mask_avx2,
                                    "aom_highbd_blend_a64_d16_mask_avx2")}));

using HbdBlendA64HMaskFunc = void (*)(uint8_t *, uint32_t, const uint8_t *,
                                      uint32_t, const uint8_t *, uint32_t,
                                      const uint8_t *, int, int, int);

class HbdCompBlendHMaskTest
    : public CompBlendTest<uint16_t, uint16_t, HbdBlendA64HMaskFunc,
                           MAKE_PARAM(HbdBlendA64HMaskFunc)> {
  public:
    HbdCompBlendHMaskTest() {
        bd_ = 8;
        func_ref_ = TEST_GET_PARAM(0);
        func_tst_ = TEST_GET_PARAM(1);
        tst_fn_name = TEST_GET_PARAM(2);
        no_sub_ = true;
    }

    void run_hbd_test(uint8_t bd) {
        bd_ = bd;
        run_test();
    }

    void run_blend(int subw, int subh) override {
        (void)subw;
        (void)subh;
        func_ref_((uint8_t *)ref_dst_,
                  dst_stride_,
                  (uint8_t *)src0_,
                  src_stride_,
                  (uint8_t *)src1_,
                  src_stride_,
                  mask_,
                  w_,
                  h_,
                  bd_);
        func_tst_((uint8_t *)tst_dst_,
                  dst_stride_,
                  (uint8_t *)src0_,
                  src_stride_,
                  (uint8_t *)src1_,
                  src_stride_,
                  mask_,
                  w_,
                  h_,
                  bd_);
    }
};

TEST_P(HbdCompBlendHMaskTest, BlendA64Mask) {
    run_hbd_test(8);
    run_hbd_test(10);
    run_hbd_test(12);
}

INSTANTIATE_TEST_CASE_P(
    BLEND, HbdCompBlendHMaskTest,
    ::testing::ValuesIn({make_tuple(aom_highbd_blend_a64_hmask_c,
                                    aom_highbd_blend_a64_hmask_sse4_1,
                                    "aom_highbd_blend_a64_hmask_sse4_1")}));

using HbdBlendA64VMaskFunc = void (*)(uint8_t *, uint32_t, const uint8_t *,
                                      uint32_t, const uint8_t *, uint32_t,
                                      const uint8_t *, int, int, int);

class HbdCompBlendVMaskTest
    : public CompBlendTest<uint16_t, uint16_t, HbdBlendA64VMaskFunc,
                           MAKE_PARAM(HbdBlendA64VMaskFunc)> {
  public:
    HbdCompBlendVMaskTest() {
        bd_ = 8;
        func_ref_ = TEST_GET_PARAM(0);
        func_tst_ = TEST_GET_PARAM(1);
        tst_fn_name = TEST_GET_PARAM(2);
        no_sub_ = true;
    }

    void run_hbd_test(uint8_t bd) {
        bd_ = bd;
        run_test();
    }

    void run_blend(int subw, int subh) override {
        (void)subw;
        (void)subh;
        func_ref_((uint8_t *)ref_dst_,
                  dst_stride_,
                  (uint8_t *)src0_,
                  src_stride_,
                  (uint8_t *)src1_,
                  src_stride_,
                  mask_,
                  w_,
                  h_,
                  bd_);
        func_tst_((uint8_t *)tst_dst_,
                  dst_stride_,
                  (uint8_t *)src0_,
                  src_stride_,
                  (uint8_t *)src1_,
                  src_stride_,
                  mask_,
                  w_,
                  h_,
                  bd_);
    }
};

TEST_P(HbdCompBlendVMaskTest, BlendA64Mask) {
    run_hbd_test(8);
    run_hbd_test(10);
    run_hbd_test(12);
}

INSTANTIATE_TEST_CASE_P(
    BLEND, HbdCompBlendVMaskTest,
    ::testing::ValuesIn({make_tuple(aom_highbd_blend_a64_vmask_c,
                                    aom_highbd_blend_a64_vmask_sse4_1,
                                    "aom_highbd_blend_a64_vmask_sse4_1")}));

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
