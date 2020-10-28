/*
* Copyright(c) 2019 Netflix, Inc.
* Copyright (c) 2016, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#include <stdlib.h>
#include <string.h>

#include "gtest/gtest.h"
#include "random.h"
#include "acm_random.h"
#include "EbTransforms.h"

#include "aom_dsp_rtcd.h"

using libaom_test::ACMRandom;

namespace {

template <typename BlendA64Func, typename SrcPixel, typename DstPixel>
class BlendA64MaskTest : public ACMRandom {
  public:
    static const int kIterations = 10000;
    static const int kMaxWidth =
        MAX_SB_SIZE * 5;  // * 5 to cover longer strides
    static const int kMaxHeight = MAX_SB_SIZE;
    static const int kBufSize = kMaxWidth * kMaxHeight;
    static const int kMaxMaskWidth = 2 * MAX_SB_SIZE;
    static const int kMaxMaskSize = kMaxMaskWidth * kMaxMaskWidth;

    virtual ~BlendA64MaskTest() {
    }
    virtual void prepare_data(int type) = 0;

    virtual void Execute(const SrcPixel *p_src0, const SrcPixel *p_src1,
                         int run_times) = 0;

    template <typename Pixel>
    void GetSources(Pixel **src0, Pixel **src1, Pixel * /*dst*/,
                    int run_times) {
        if (run_times > 1) {
            *src0 = src0_;
            *src1 = src1_;
            return;
        }
        switch (this->PseudoUniform(3)) {
        case 0:  // Separate sources
            *src0 = src0_;
            *src1 = src1_;
            break;
        case 1:  // src0 == dst
            *src0 = dst_tst_;
            src0_stride_ = dst_stride_;
            src0_offset_ = dst_offset_;
            *src1 = src1_;
            break;
        case 2:  // src1 == dst
            *src0 = src0_;
            *src1 = dst_tst_;
            src1_stride_ = dst_stride_;
            src1_offset_ = dst_offset_;
            break;
        default: FAIL();
        }
    }

    void GetSources(uint16_t **src0, uint16_t **src1, uint8_t * /*dst*/,
                    int /*run_times*/) {
        *src0 = src0_;
        *src1 = src1_;
    }

    uint8_t Rand1() {
        return this->Rand8() & 1;
    }

    void RunOneTest(int block_size, int subx, int suby, int run_times) {
        w_ = block_size_wide[block_size];
        h_ = block_size_high[block_size];
        run_times = run_times > 1 ? run_times / w_ : 1;
        ASSERT_GT(run_times, 0);
        subx_ = subx;
        suby_ = suby;

        dst_offset_ = this->PseudoUniform(33);
        dst_stride_ = this->PseudoUniform(kMaxWidth + 1 - w_) + w_;

        src0_offset_ = this->PseudoUniform(33);
        src0_stride_ = this->PseudoUniform(kMaxWidth + 1 - w_) + w_;

        src1_offset_ = this->PseudoUniform(33);
        src1_stride_ = this->PseudoUniform(kMaxWidth + 1 - w_) + w_;

        mask_stride_ =
            this->PseudoUniform(kMaxWidth + 1 - w_ * (subx_ ? 2 : 1)) +
            w_ * (subx_ ? 2 : 1);

        SrcPixel *p_src0;
        SrcPixel *p_src1;

        p_src0 = src0_;
        p_src1 = src1_;

        GetSources(&p_src0, &p_src1, &dst_ref_[0], run_times);

        Execute(p_src0, p_src1, run_times);

        for (int r = 0; r < h_; ++r) {
            for (int c = 0; c < w_; ++c) {
                ASSERT_EQ(dst_ref_[dst_offset_ + r * dst_stride_ + c],
                    dst_tst_[dst_offset_ + r * dst_stride_ + c])
                    << w_ << "x" << h_ << " subx " << subx_ << " suby "
                    << suby_ << " r: " << r << " c: " << c;
            }
        }
    }

    void RunTest(int type) {
        int run_times = 1;
        for (int iter = 0;
             iter < kIterations && !::testing::Test::HasFatalFailure();
             ++iter) {
            prepare_data(type);
            int block_size = this->Rand8() % BLOCK_SIZES_ALL;
            subx_ = Rand1();
            suby_ = Rand1();
            RunOneTest(block_size, subx_, suby_, run_times);
        }
    }

    void Run() {
        RunTest(0);  // RandomValues
        RunTest(1);  // ExtremeValues
    }

    DstPixel dst_ref_[kBufSize];
    DstPixel dst_tst_[kBufSize];
    uint32_t dst_stride_;
    uint32_t dst_offset_;

    SrcPixel src0_[kBufSize];
    uint32_t src0_stride_;
    uint32_t src0_offset_;

    SrcPixel src1_[kBufSize];
    uint32_t src1_stride_;
    uint32_t src1_offset_;

    uint8_t mask_[kMaxMaskSize];
    size_t mask_stride_;

    int w_;
    int h_;

    int suby_;
    int subx_;

    BlendA64Func ref_func_;
    BlendA64Func tst_func_;
};

//////////////////////////////////////////////////////////////////////////////
// 8 bit version
//////////////////////////////////////////////////////////////////////////////

typedef void (*F8B)(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0,
                    uint32_t src0_stride, const uint8_t *src1,
                    uint32_t src1_stride, const uint8_t *mask,
                    uint32_t mask_stride, int w, int h, int subx, int suby);

class BlendA64MaskTest8B : public BlendA64MaskTest<F8B, uint8_t, uint8_t> {
  public:
    BlendA64MaskTest8B(F8B ref, F8B tst) {
        ref_func_ = ref;
        tst_func_ = tst;
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
    void Execute(const uint8_t *p_src0, const uint8_t *p_src1,
                 int run_times) override {
        for (int i = 0; i < run_times; ++i) {
            ref_func_(dst_ref_ + dst_offset_, dst_stride_, p_src0 + src0_offset_,
                      src0_stride_, p_src1 + src1_offset_, src1_stride_,
                      mask_, kMaxMaskWidth, w_, h_, subx_, suby_);
        }
        for (int i = 0; i < run_times; ++i) {
            tst_func_(dst_tst_ + dst_offset_, dst_stride_, p_src0 + src0_offset_,
                      src0_stride_, p_src1 + src1_offset_, src1_stride_,
                      mask_, kMaxMaskWidth, w_, h_, subx_, suby_);
        }
    }
};

#define TEST_CLASS(type_name, ref, tst, match_test) \
    TEST(type_name, match_test) {                   \
        type_name *test = new type_name(ref, tst);  \
        test->Run();                                \
        delete test;                                \
    }

TEST_CLASS(BlendA64MaskTest8B, svt_aom_blend_a64_mask_c, svt_aom_blend_a64_mask_sse4_1,
           Mask_Blend_SSE4_1)
TEST_CLASS(BlendA64MaskTest8B, svt_aom_blend_a64_mask_sse4_1,
           svt_aom_blend_a64_mask_avx2, Mask_Blend_AVX2)

//////////////////////////////////////////////////////////////////////////////
// 8 bit _d16 version
//////////////////////////////////////////////////////////////////////////////

typedef void (*F8B_D16)(uint8_t *dst, uint32_t dst_stride, const uint16_t *src0,
                        uint32_t src0_stride, const uint16_t *src1,
                        uint32_t src1_stride, const uint8_t *mask,
                        uint32_t mask_stride, int w, int h, int subx, int suby,
                        ConvolveParams *conv_params);

class BlendA64MaskTest8B_d16
    : public BlendA64MaskTest<F8B_D16, uint16_t, uint8_t> {
  // max number of bits used by the source
  static const int kSrcMaxBitsMask = 0x3fff;

public:
  BlendA64MaskTest8B_d16(F8B_D16 ref, F8B_D16 tst) {
      ref_func_ = ref;
      tst_func_ = tst;
  }

  void prepare_data(int type) override {
      UNUSED(type);
      if (!type) {
          for (int i = 0; i < kBufSize; ++i) {
              dst_ref_[i] = this->Rand8();
              dst_tst_[i] = this->Rand8();
              src0_[i] = this->Rand16() & kSrcMaxBitsMask;
              src1_[i] = this->Rand16() & kSrcMaxBitsMask;
          }
          for (int i = 0; i < kMaxMaskSize; ++i)
              mask_[i] = this->PseudoUniform(AOM_BLEND_A64_MAX_ALPHA + 1);
      }
      else {
          for (int i = 0; i < kBufSize; ++i) {
              dst_ref_[i] = (uint16_t)255;
              dst_tst_[i] = (uint16_t)255;
              src0_[i] = (uint16_t)kSrcMaxBitsMask;
              src1_[i] = (uint16_t)kSrcMaxBitsMask;
          }
          for (int i = 0; i < kMaxMaskSize; ++i)
              mask_[i] = AOM_BLEND_A64_MAX_ALPHA - 1;
      }
  }

  void Execute(const uint16_t *p_src0, const uint16_t *p_src1,
      int run_times) override
  {
    ConvolveParams conv_params;
    conv_params.round_0 = ROUND0_BITS;
    conv_params.round_1 = COMPOUND_ROUND1_BITS;
    for (int i = 0; i < run_times; ++i) {
      ref_func_(dst_ref_ + dst_offset_, dst_stride_,
                       p_src0 + src0_offset_, src0_stride_,
                       p_src1 + src1_offset_, src1_stride_, mask_,
                       kMaxMaskWidth, w_, h_, subx_, suby_, &conv_params);
    }
    for (int i = 0; i < run_times; ++i) {
      tst_func_(dst_tst_ + dst_offset_, dst_stride_,
                       p_src0 + src0_offset_, src0_stride_,
                       p_src1 + src1_offset_, src1_stride_, mask_,
                       kMaxMaskWidth, w_, h_, subx_, suby_, &conv_params);
    }
  }
};

TEST_CLASS(BlendA64MaskTest8B_d16, svt_aom_lowbd_blend_a64_d16_mask_c,
           svt_aom_lowbd_blend_a64_d16_mask_avx2, Mask_Blend_d16_AVX2)

//////////////////////////////////////////////////////////////////////////////
// High bit-depth version
//////////////////////////////////////////////////////////////////////////////

typedef void (*FHBD)(uint8_t *dst, uint32_t dst_stride, const uint8_t *src0,
                     uint32_t src0_stride, const uint8_t *src1,
                     uint32_t src1_stride, const uint8_t *mask,
                     uint32_t mask_stride, int w, int h, int subx, int suby,
                     int bd);

class BlendA64MaskTestHBD : public BlendA64MaskTest<FHBD, uint16_t, uint16_t> {
  public:
    BlendA64MaskTestHBD(FHBD ref, FHBD tst) {
        ref_func_ = ref;
        tst_func_ = tst;
    }
    void prepare_data(int type) override {
        if (!type) {
            switch (this->PseudoUniform(3)) {
            case 0: bit_depth_ = 8; break;
            case 1: bit_depth_ = 10; break;
            default: bit_depth_ = 12; break;
            }
            const int hi = 1 << bit_depth_;
            for (int i = 0; i < kBufSize; ++i) {
                dst_ref_[i] = this->PseudoUniform(hi);
                dst_tst_[i] = this->PseudoUniform(hi);
                src0_[i] = this->PseudoUniform(hi);
                src1_[i] = this->PseudoUniform(hi);
            }
            for (int i = 0; i < kMaxMaskSize; ++i)
                mask_[i] = this->PseudoUniform(AOM_BLEND_A64_MAX_ALPHA + 1);
        } else {
            switch (this->PseudoUniform(3)) {
            case 0: bit_depth_ = 8; break;
            case 1: bit_depth_ = 10; break;
            default: bit_depth_ = 12; break;
            }

            const int hi = 1 << bit_depth_;
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
    void Execute(const uint16_t *p_src0, const uint16_t *p_src1,
                 int run_times) override {
        for (int i = 0; i < run_times; ++i) {
            ref_func_((uint8_t *)(dst_ref_ + dst_offset_), dst_stride_,
                      (uint8_t *)(p_src0 + src0_offset_), src0_stride_,
                      (uint8_t *)(p_src1 + src1_offset_), src1_stride_,
                      mask_, kMaxMaskWidth, w_,  h_, subx_, suby_, bit_depth_);
        }
        for (int i = 0; i < run_times; ++i) {
            tst_func_((uint8_t *)(dst_tst_ + dst_offset_), dst_stride_,
                      (uint8_t *)(p_src0 + src0_offset_), src0_stride_,
                      (uint8_t *)(p_src1 + src1_offset_), src1_stride_,
                      mask_, kMaxMaskWidth, w_, h_, subx_, suby_, bit_depth_);
        }
    }

    int bit_depth_;
};

TEST_CLASS(BlendA64MaskTestHBD, svt_aom_highbd_blend_a64_mask_c,
           svt_aom_highbd_blend_a64_mask_8bit_sse4_1, Mask_Blend_Hbd_SSE4_1)

//////////////////////////////////////////////////////////////////////////////
// HBD _d16 version
//////////////////////////////////////////////////////////////////////////////

typedef void (*FHBD_D16)(uint8_t *dst, uint32_t dst_stride,
                         const CONV_BUF_TYPE *src0, uint32_t src0_stride,
                         const CONV_BUF_TYPE *src1, uint32_t src1_stride,
                         const uint8_t *mask, uint32_t mask_stride, int w,
                         int h, int subx, int suby, ConvolveParams *conv_params,
                         const int bd);

class BlendA64MaskTestHBD_d16
    : public BlendA64MaskTest<FHBD_D16, uint16_t, uint16_t> {
 public:
  // max number of bits used by the source
  static const int kSrcMaxBitsMask = (1 << 14) - 1;
  static const int kSrcMaxBitsMaskHBD = (1 << 16) - 1;
  BlendA64MaskTestHBD_d16(FHBD_D16 ref, FHBD_D16 tst) {
      ref_func_ = ref;
      tst_func_ = tst;
  }
  void prepare_data(int type) override {
      if (!type) {
          switch (this->PseudoUniform(3)) {
          case 0: bit_depth_ = 8; break;
          case 1: bit_depth_ = 10; break;
          default: bit_depth_ = 12; break;
          }
          src_max_bits_mask_ =
              (bit_depth_ == 8) ? kSrcMaxBitsMask : kSrcMaxBitsMaskHBD;

          for (int i = 0; i < kBufSize; ++i) {
              dst_ref_[i] = this->Rand8();
              dst_tst_[i] = this->Rand8();
              src0_[i] = this->Rand16() & src_max_bits_mask_;
              src1_[i] = this->Rand16() & src_max_bits_mask_;
          }
          for (int i = 0; i < kMaxMaskSize; ++i)
              mask_[i] = this->PseudoUniform(AOM_BLEND_A64_MAX_ALPHA + 1);
      } else {
          switch (this->PseudoUniform(3)) {
          case 0: bit_depth_ = 8; break;
          case 1: bit_depth_ = 10; break;
          default: bit_depth_ = 12; break;
          }

          const int hi = 1 << bit_depth_;
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
  void Execute(const uint16_t *p_src0, const uint16_t *p_src1, int run_times) override
  {
    ASSERT_GT(run_times, 0) << "Cannot run 0 iterations of the test.";
    ConvolveParams conv_params;
    conv_params.round_0 = (bit_depth_ == 12) ? ROUND0_BITS + 2 : ROUND0_BITS;
    conv_params.round_1 = COMPOUND_ROUND1_BITS;
    for (int i = 0; i < run_times; ++i) {
      ref_func_((uint8_t *)(dst_ref_ + dst_offset_), dst_stride_,
                p_src0 + src0_offset_, src0_stride_,
                p_src1 + src1_offset_, src1_stride_, mask_,
                kMaxMaskWidth, w_, h_, subx_, suby_, &conv_params,
                bit_depth_);
    }
    if (tst_func_) {
      for (int i = 0; i < run_times; ++i) {
        tst_func_((uint8_t *)(dst_tst_ + dst_offset_),
                  dst_stride_, p_src0 + src0_offset_, src0_stride_,
                  p_src1 + src1_offset_, src1_stride_, mask_,
                  kMaxMaskWidth, w_, h_, subx_, suby_, &conv_params,
                  bit_depth_);
      }
    }
  }

  int bit_depth_;
  int src_max_bits_mask_;
};

TEST_CLASS(BlendA64MaskTestHBD_d16, svt_aom_highbd_blend_a64_d16_mask_c,
           svt_aom_highbd_blend_a64_d16_mask_avx2, _Mask_Blend_Hbd_d16_AVX2)
}; // namespace
