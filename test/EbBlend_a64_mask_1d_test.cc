/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <stdlib.h>
#include <string.h>

#include "gtest/gtest.h"
#include "random.h"
#include "acm_random.h"

#include "aom_dsp_rtcd.h"

using libaom_test::ACMRandom;

namespace {

static void blend_a64_hmask_ref(uint8_t *dst, uint32_t dst_stride,
                                const uint8_t *src0, uint32_t src0_stride,
                                const uint8_t *src1, uint32_t src1_stride,
                                const uint8_t *mask, int w, int h);
static void blend_a64_vmask_ref(uint8_t *dst, uint32_t dst_stride,
                                const uint8_t *src0, uint32_t src0_stride,
                                const uint8_t *src1, uint32_t src1_stride,
                                const uint8_t *mask, int w, int h);
static void highbd_blend_a64_hmask_ref(
    uint8_t *dst, uint32_t dst_stride, const uint8_t *src0,
    uint32_t src0_stride, const uint8_t *src1, uint32_t src1_stride,
    const uint8_t *mask, int w, int h, int bd);
static void highbd_blend_a64_vmask_ref(
    uint8_t *dst, uint32_t dst_stride, const uint8_t *src0,
    uint32_t src0_stride, const uint8_t *src1, uint32_t src1_stride,
    const uint8_t *mask, int w, int h, int bd);

template <typename F, typename T>
class BlendA64Mask1DTest : public ACMRandom {
  public:
    static const int kIterations = 10000;
    static const int kMaxWidth =
        MAX_SB_SIZE * 5;  // * 5 to cover longer strides
    static const int kMaxHeight = MAX_SB_SIZE;
    static const int kBufSize = kMaxWidth * kMaxHeight;
    static const int kMaxMaskWidth = 2 * MAX_SB_SIZE;
    static const int kMaxMaskSize = kMaxMaskWidth;

    virtual ~BlendA64Mask1DTest() {
    }

    virtual void Execute(const T *, const T *) {};
    virtual void prepare_data(int) {};

    void Common(int type) {
        prepare_data(type);

        w_ = 2 << this->PseudoUniform(MAX_SB_SIZE_LOG2);
        h_ = 2 << this->PseudoUniform(MAX_SB_SIZE_LOG2);

        dst_offset_ = this->PseudoUniform(33);
        dst_stride_ = this->PseudoUniform(kMaxWidth + 1 - w_) + w_;

        src0_offset_ = this->PseudoUniform(33);
        src0_stride_ = this->PseudoUniform(kMaxWidth + 1 - w_) + w_;

        src1_offset_ = this->PseudoUniform(33);
        src1_stride_ = this->PseudoUniform(kMaxWidth + 1 - w_) + w_;

        T *p_src0;
        T *p_src1;

        switch (this->PseudoUniform(3)) {
        case 0:  // Separate sources
            p_src0 = src0_;
            p_src1 = src1_;
            break;
        case 1:  // src0 == dst
            p_src0 = dst_tst_;
            src0_stride_ = dst_stride_;
            src0_offset_ = dst_offset_;
            p_src1 = src1_;
            break;
        case 2:  // src1 == dst
            p_src0 = src0_;
            p_src1 = dst_tst_;
            src1_stride_ = dst_stride_;
            src1_offset_ = dst_offset_;
            break;
        default: FAIL();
        }

        Execute(p_src0, p_src1);

        for (int r = 0; r < h_; ++r) {
            for (int c = 0; c < w_; ++c) {
                ASSERT_EQ(dst_ref_[dst_offset_ + r * dst_stride_ + c],
                          dst_tst_[dst_offset_ + r * dst_stride_ + c]);
            }
        }
    }
    void RunTest() {
        for (int iter = 0;
             iter < kIterations && !::testing::Test::HasFatalFailure();
             iter++) {
            Common(0);  /*RandomValues*/
            Common(1);  /*ExtremeValues*/
        }
    }

    T dst_ref_[kBufSize];
    T dst_tst_[kBufSize];
    uint32_t dst_stride_;
    uint32_t dst_offset_;

    T src0_[kBufSize];
    uint32_t src0_stride_;
    uint32_t src0_offset_;

    T src1_[kBufSize];
    uint32_t src1_stride_;
    uint32_t src1_offset_;

    uint8_t mask_[kMaxMaskSize];

    int w_;
    int h_;
    int bd_;

    F ref_func_;
    F tst_func_;
};

//////////////////////////////////////////////////////////////////////////////
// 8 bit version
//////////////////////////////////////////////////////////////////////////////

typedef void (*Blend8B)(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0,
                        uint32_t src0_stride, const uint8_t *src1,
                        uint32_t src1_stride, const uint8_t *mask, int w,
                        int h);

class BlendA64Mask1DTest8B : public BlendA64Mask1DTest<Blend8B, uint8_t> {
  public:
    BlendA64Mask1DTest8B(Blend8B ref, Blend8B tst) {
        ref_func_ = ref;
        tst_func_ = tst;
        bd_ = 8;
    }
    void prepare_data(int type) override {
        if (!type) {
            for (int i = 0; i < kBufSize; ++i) {
                dst_ref_[i] = this->Rand8();
                dst_tst_[i] = this->Rand8();

                src0_[i] = this->Rand8();
                src1_[i] = this->Rand8();
            }

            for (int i = 0; i < kMaxMaskSize; ++i)
                mask_[i] = this->PseudoUniform(AOM_BLEND_A64_MAX_ALPHA + 1);
        } else {
            for (int i = 0; i < kBufSize; ++i) {
                dst_ref_[i] = this->PseudoUniform(2) + 254;
                dst_tst_[i] = this->PseudoUniform(2) + 254;
                src0_[i] = this->PseudoUniform(2) + 254;
                src1_[i] = this->PseudoUniform(2) + 254;
            }

            for (int i = 0; i < kMaxMaskSize; ++i)
                mask_[i] = this->PseudoUniform(2) + AOM_BLEND_A64_MAX_ALPHA - 1;
        }
    }
    void Execute(const uint8_t *p_src0, const uint8_t *p_src1) override {
        ref_func_(dst_ref_ + dst_offset_, dst_stride_, p_src0 + src0_offset_,
                  src0_stride_, p_src1 + src1_offset_, src1_stride_,
                  mask_, w_, h_);
        tst_func_(dst_tst_ + dst_offset_, dst_stride_, p_src0 + src0_offset_,
                  src0_stride_, p_src1 + src1_offset_, src1_stride_,
                  mask_, w_, h_);
    }
};

static void blend_a64_hmask_ref(uint8_t *dst, uint32_t dst_stride,
                                const uint8_t *src0, uint32_t src0_stride,
                                const uint8_t *src1, uint32_t src1_stride,
                                const uint8_t *mask, int w, int h) {
    uint8_t mask2d[BlendA64Mask1DTest8B::kMaxMaskSize]
                  [BlendA64Mask1DTest8B::kMaxMaskSize];

    for (int row = 0; row < h; ++row)
        for (int col = 0; col < w; ++col)
            mask2d[row][col] = mask[col];

    aom_blend_a64_mask_c(dst, dst_stride, src0, src0_stride, src1, src1_stride,
                         &mask2d[0][0], BlendA64Mask1DTest8B::kMaxMaskSize,
                         w, h, 0, 0);
}

static void blend_a64_vmask_ref(uint8_t *dst, uint32_t dst_stride,
                                const uint8_t *src0, uint32_t src0_stride,
                                const uint8_t *src1, uint32_t src1_stride,
                                const uint8_t *mask, int w, int h) {
    uint8_t mask2d[BlendA64Mask1DTest8B::kMaxMaskSize]
                  [BlendA64Mask1DTest8B::kMaxMaskSize];

    for (int row = 0; row < h; ++row)
        for (int col = 0; col < w; ++col)
            mask2d[row][col] = mask[row];

    aom_blend_a64_mask_c(dst, dst_stride, src0, src0_stride, src1, src1_stride,
                         &mask2d[0][0], BlendA64Mask1DTest8B::kMaxMaskSize,
                         w, h, 0, 0);
}

#define TEST_CLASS(type_name, ref, tst, match_test)     \
    TEST(type_name, match_test) {            \
        type_name *test = new type_name(ref, tst); \
        test->RunTest();                   \
        delete test;                       \
    }
// C
TEST_CLASS(BlendA64Mask1DTest8B, blend_a64_hmask_ref,
          aom_blend_a64_hmask_c, Horz_Blend_C)
TEST_CLASS(BlendA64Mask1DTest8B, blend_a64_vmask_ref,
           aom_blend_a64_vmask_c, Vert_Blend_C)
// Intrinsic
TEST_CLASS(BlendA64Mask1DTest8B, blend_a64_hmask_ref,
           aom_blend_a64_hmask_sse4_1, Horz_Blend_SSE4_1)
TEST_CLASS(BlendA64Mask1DTest8B, blend_a64_vmask_ref,
           aom_blend_a64_vmask_sse4_1, Vert_Blend_SSE4_1)

//////////////////////////////////////////////////////////////////////////////
// HBD version
//////////////////////////////////////////////////////////////////////////////
typedef void (*BlendHBD)(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0,
                         uint32_t src0_stride, const uint8_t *src1,
                         uint32_t src1_stride, const uint8_t *mask, int w,
                         int h, int bd);

class BlendA64Mask1DTestHBD : public BlendA64Mask1DTest<BlendHBD, uint16_t> {
  public:
    BlendA64Mask1DTestHBD(BlendHBD ref, BlendHBD tst) {
        ref_func_ = ref;
        tst_func_ = tst;
        bd_ = 10;
    }
    void prepare_data(int type) override {
        switch (this->PseudoUniform(3)) {
            case 0: bd_ = 8; break;
            case 1: bd_ = 10; break;
            default: bd_ = 12; break;
        }
        if (!type) {
            const int hi = 1 << bd_;

            for (int i = 0; i < kBufSize; ++i) {
              dst_ref_[i] = this->PseudoUniform(hi);
              dst_tst_[i] = this->PseudoUniform(hi);
              src0_[i] = this->PseudoUniform(hi);
              src1_[i] = this->PseudoUniform(hi);
            }

            for (int i = 0; i < kMaxMaskSize; ++i)
              mask_[i] = this->PseudoUniform(AOM_BLEND_A64_MAX_ALPHA + 1);

        } else {
            const int hi = 1 << bd_;
            const int lo = hi - 2;

            for (int i = 0; i < kBufSize; ++i) {
                dst_ref_[i] = this->PseudoUniform(hi - lo) + lo;
                dst_tst_[i] = this->PseudoUniform(hi - lo) + lo;
                src0_[i] = this->PseudoUniform(hi - lo) + lo;
                src1_[i] = this->PseudoUniform(hi - lo) + lo;
            }

            for (int i = 0; i < kMaxMaskSize; ++i)
                mask_[i] = this->PseudoUniform(2) + AOM_BLEND_A64_MAX_ALPHA - 1;
        }
    }
    void Execute(const uint16_t *p_src0, const uint16_t *p_src1) override {
        ref_func_((uint8_t *)(dst_ref_ + dst_offset_), dst_stride_,
                  (uint8_t *)(p_src0 + src0_offset_), src0_stride_,
                  (uint8_t *)(p_src1 + src1_offset_), src1_stride_,
                  mask_, w_, h_, bd_);
        tst_func_((uint8_t *)(dst_tst_ + dst_offset_), dst_stride_,
                  (uint8_t *)(p_src0 + src0_offset_), src0_stride_,
                  (uint8_t *)(p_src1 + src1_offset_), src1_stride_,
                  mask_, w_, h_, bd_);
    }
};

static void highbd_blend_a64_hmask_ref(
    uint8_t *dst, uint32_t dst_stride, const uint8_t *src0,
    uint32_t src0_stride, const uint8_t *src1, uint32_t src1_stride,
    const uint8_t *mask, int w, int h, int bd) {
    uint8_t mask2d[BlendA64Mask1DTestHBD::kMaxMaskSize]
                  [BlendA64Mask1DTestHBD::kMaxMaskSize];

    for (int row = 0; row < h; ++row)
        for (int col = 0; col < w; ++col)
            mask2d[row][col] = mask[col];

    aom_highbd_blend_a64_mask_c(dst, dst_stride, src0, src0_stride, src1,
        src1_stride, &mask2d[0][0], BlendA64Mask1DTestHBD::kMaxMaskSize,
        w, h, 0, 0, bd);
}

static void highbd_blend_a64_vmask_ref(
    uint8_t *dst, uint32_t dst_stride, const uint8_t *src0,
    uint32_t src0_stride, const uint8_t *src1, uint32_t src1_stride,
    const uint8_t *mask, int w, int h, int bd) {
    uint8_t mask2d[BlendA64Mask1DTestHBD::kMaxMaskSize]
                  [BlendA64Mask1DTestHBD::kMaxMaskSize];

    for (int row = 0; row < h; ++row)
        for (int col = 0; col < w; ++col)
            mask2d[row][col] = mask[row];

    aom_highbd_blend_a64_mask_c(dst, dst_stride, src0, src0_stride, src1,
        src1_stride, &mask2d[0][0], BlendA64Mask1DTestHBD::kMaxMaskSize,
        w, h, 0, 0,bd);
}

// C
TEST_CLASS(BlendA64Mask1DTestHBD, highbd_blend_a64_hmask_ref,
           aom_highbd_blend_a64_hmask_c, Horz_Blend_Hbd_C)
TEST_CLASS(BlendA64Mask1DTestHBD, highbd_blend_a64_vmask_ref,
           aom_highbd_blend_a64_vmask_c, Vert_Blend_Hbd_C)
// Intrinsic
TEST_CLASS(BlendA64Mask1DTestHBD, highbd_blend_a64_hmask_ref,
           aom_highbd_blend_a64_hmask_sse4_1, Horz_Blend_Hbd_SSE4_1)
TEST_CLASS(BlendA64Mask1DTestHBD, highbd_blend_a64_vmask_ref,
           aom_highbd_blend_a64_vmask_sse4_1, Vert_Blend_Hbd_SSE4_1)

}; // namespace
