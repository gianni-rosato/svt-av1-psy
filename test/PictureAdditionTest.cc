/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file PictureAdditionTest.cc
 *
 * @brief Unit test for Residual functions:
 * - picture_addition_kernel{m}x{n}_{sse,sse2}_intrin
 * - picture_addition_kernel{m}x{n}_av1_sse2_intrin
 *
 * @author Cidana-Ivy
 *
 ******************************************************************************/
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <limits.h>
#include <new>
// workaround to eliminate the compiling warning on linux
// The macro will conflict with definition in gtest.h
#ifdef __USE_GNU
#undef __USE_GNU  // defined in EbThreads.h
#endif
#ifdef _GNU_SOURCE
#undef _GNU_SOURCE  // defined in EbThreads.h
#endif

#include "EbPictureOperators.h"
#include "EbEncIntraPrediction.h"
#include "random.h"
#include "util.h"

#if !REMOVE_UNUSED_CODE
using svt_av1_test_tool::SVTRandom;  // to generate the random

namespace {
typedef enum { PRED_MIN, PRED_MAX, ALL_MIN, ALL_MAX, ALL_RANDOM } TestPattern;

typedef std::tuple<uint32_t, uint32_t> AreaSize;

typedef void (*LbdPictureAdditionFunc)(uint8_t *pred_ptr, uint32_t pred_stride,
                                       int16_t *residual_ptr,
                                       uint32_t residual_stride,
                                       uint8_t *recon_ptr,
                                       uint32_t recon_stride, uint32_t width,
                                       uint32_t height);
typedef void (*HbdPictureAdditionFunc)(uint8_t *pred_ptr, uint32_t pred_stride,
                                       int32_t *residual_ptr,
                                       uint32_t residual_stride,
                                       uint8_t *recon_ptr,
                                       uint32_t recon_stride, uint32_t width,
                                       uint32_t height, int32_t bd);

typedef struct PictureAdditionTestParam {
    AreaSize area_size;
    LbdPictureAdditionFunc lbd_test_func;
    HbdPictureAdditionFunc hbd_test_func;
} PictureAdditionTestParam;
typedef std::tuple<PictureAdditionTestParam, TestPattern> TestParam;
/**
 * @brief Unit test for Picture Addition functions include:
 *  - picture_addition_kernel{m}x{n}_{sse,sse2}_intrin
 *  - picture_addition_kernel{m}x{n}_av1_sse2_intrin
 *
 * Test strategy:
 * This test case combines different AreaSize and
 * different test pattern. Check the results of lbd function and hbd
 * function.
 *
 * Expect result:
 * Results from lbd function and hbd function are
 * equal.
 *
 * Test cases:
 *  Width {4, 8, 16, 32, 64} x height{ 4, 8, 16, 32, 64}
 *  Test vector pattern: different combination of pred and residual
 *
 */
class PictureAdditionTest : public ::testing::Test,
                            public ::testing::WithParamInterface<TestParam> {
  public:
    PictureAdditionTest()
        : test_pattern_(TEST_GET_PARAM(1)),
          area_width_(std::get<0>(TEST_GET_PARAM(0).area_size)),
          area_height_(std::get<1>(TEST_GET_PARAM(0).area_size)),
          lbd_test_func_(TEST_GET_PARAM(0).lbd_test_func),
          hbd_test_func_(TEST_GET_PARAM(0).hbd_test_func) {
        residual_ = nullptr;
        residual_int32_ = nullptr;
        pred_ = nullptr;
        pred16bit_ = nullptr;
        recon_1_ = recon_2_ = recon_c_ = nullptr;
        recon1_16bit_ = recon2_16bit_ = nullptr;
        recon_stride_ = pred_stride_ = residual_stride_ = MAX_PU_SIZE;
        test_size_ = MAX_PU_SIZE * MAX_PU_SIZE;
    }

    void SetUp() override {
        residual_ = reinterpret_cast<int16_t *>(
            eb_aom_memalign(32, sizeof(int16_t) * test_size_));
        residual_int32_ = reinterpret_cast<int32_t *>(
            eb_aom_memalign(32, sizeof(int32_t) * test_size_));
        pred_ = reinterpret_cast<uint8_t *>(eb_aom_memalign(32, test_size_));
        recon_1_ = reinterpret_cast<uint8_t *>(eb_aom_memalign(32, test_size_));
        recon_2_ = reinterpret_cast<uint8_t *>(eb_aom_memalign(32, test_size_));
        recon_c_ = reinterpret_cast<uint8_t *>(eb_aom_memalign(32, test_size_));
        pred16bit_ = reinterpret_cast<uint16_t *>(
            eb_aom_memalign(32, sizeof(uint16_t) * test_size_));
        recon1_16bit_ = reinterpret_cast<uint16_t *>(
            eb_aom_memalign(32, sizeof(uint16_t) * test_size_));
        recon2_16bit_ = reinterpret_cast<uint16_t *>(
            eb_aom_memalign(32, sizeof(uint16_t) * test_size_));
        memset(residual_int32_, 0, test_size_ * sizeof(residual_int32_[0]));
        memset(recon_1_, 0, test_size_ * sizeof(recon_1_[0]));
        memset(recon_2_, 0, test_size_ * sizeof(recon_2_[0]));
        memset(recon_c_, 0, test_size_ * sizeof(recon_c_[0]));
        memset(recon1_16bit_, 0, test_size_ * sizeof(recon1_16bit_[0]));
        memset(recon2_16bit_, 0, test_size_ * sizeof(recon2_16bit_[0]));
    }

    void TearDown() override {
        if (residual_)
            eb_aom_free(residual_);
        if (residual_int32_)
            eb_aom_free(residual_int32_);
        if (pred_)
            eb_aom_free(pred_);
        if (recon1_16bit_)
            eb_aom_free(recon1_16bit_);
        if (recon2_16bit_)
            eb_aom_free(recon2_16bit_);
        if (pred16bit_)
            eb_aom_free(pred16bit_);
        if (recon_1_)
            eb_aom_free(recon_1_);
        if (recon_2_)
            eb_aom_free(recon_2_);
        if (recon_c_)
            eb_aom_free(recon_c_);
    }

  protected:
    void prepare_data() {
        // Max bitdepth is 12, set the mask accordingly.
        const uint8_t mask_pred = (1 << 8) - 1;
        const int16_t mask_res = (1 << 12) - 1;
        SVTRandom rnd_uint8(0, mask_pred);
        SVTRandom rnd_int16(0, mask_res);
        switch (test_pattern_) {
        case PRED_MIN: {
            for (uint32_t i = 0; i < test_size_; i++) {
                pred_[i] = 0;
                residual_[i] = residual_int32_[i] = mask_res;
            }
            break;
        }
        case PRED_MAX: {
            for (uint32_t i = 0; i < test_size_; i++) {
                pred_[i] = mask_pred;
                residual_[i] = residual_int32_[i] = 0;
            }
            break;
        }
        case ALL_MIN: {
            for (uint32_t i = 0; i < test_size_; i++) {
                pred_[i] = 0;
                residual_[i] = residual_int32_[i] = 0;
            }
            break;
        }
        case ALL_MAX: {
            for (uint32_t i = 0; i < test_size_; i++) {
                pred_[i] = mask_pred;
                residual_[i] = residual_int32_[i] = mask_res;
            }
            break;
        }
        case ALL_RANDOM: {
            for (uint32_t i = 0; i < test_size_; i++) {
                pred_[i] = rnd_uint8.random();
                residual_[i] = residual_int32_[i] = rnd_int16.random();
            }
            break;
        }
        default: break;
        }
    }

    void prepare_16bit_data() {
        // Max bitdepth is 12, set the mask accordingly.
        const uint16_t mask_pred = (1 << 12) - 1;
        const int16_t mask_res = (1 << 15) - 1;
        SVTRandom rnd_uint16(0, mask_pred);
        SVTRandom rnd_int16(0, mask_res);
        switch (test_pattern_) {
        case PRED_MIN: {
            for (uint32_t i = 0; i < test_size_; i++) {
                pred16bit_[i] = 0;
                residual_[i] = residual_int32_[i] = mask_res;
            }
            break;
        }
        case PRED_MAX: {
            for (uint32_t i = 0; i < test_size_; i++) {
                pred16bit_[i] = mask_pred;
                residual_[i] = 0;
                residual_int32_[i] = 0;
            }
            break;
        }
        case ALL_MIN: {
            for (uint32_t i = 0; i < test_size_; i++) {
                pred16bit_[i] = 0;
                residual_[i] = 0;
                residual_int32_[i] = 0;
            }
            break;
        }
        case ALL_MAX: {
            for (uint32_t i = 0; i < test_size_; i++) {
                pred16bit_[i] = mask_pred;
                residual_[i] = residual_int32_[i] = mask_res;
            }
            break;
        }
        case ALL_RANDOM: {
            for (uint32_t i = 0; i < test_size_; i++) {
                pred16bit_[i] = rnd_uint16.random();
                residual_[i] = residual_int32_[i] = rnd_int16.random();
            }
            break;
        }
        default: break;
        }
    }

    void check_output(uint32_t width, uint32_t height, uint8_t *out_1,
                      uint8_t *out_2) {
        int fail_count = 0;
        for (uint32_t j = 0; j < height; j++) {
            for (uint32_t k = 0; k < width; k++) {
                if (out_1[k + j * recon_stride_] !=
                    out_2[k + j * recon_stride_])
                    fail_count++;
            }
        }
        EXPECT_EQ(0, fail_count)
            << "compare picture addition result error"
            << "in test area for " << fail_count << "times";
    }

    void check_16bit_output(uint32_t width, uint32_t height) {
        int fail_count = 0;
        for (uint32_t j = 0; j < height; j++) {
            for (uint32_t k = 0; k < width; k++) {
                if (recon1_16bit_[k + j * recon_stride_] !=
                    recon2_16bit_[k + j * recon_stride_])
                    fail_count++;
            }
        }
        EXPECT_EQ(0, fail_count)
            << "compare picture addition result error"
            << "in test area for " << fail_count << "times";
    }

    void run_test() {
        prepare_data();
        lbd_test_func_(pred_,
                       pred_stride_,
                       residual_,
                       residual_stride_,
                       recon_1_,
                       recon_stride_,
                       area_width_,
                       area_height_);
        hbd_test_func_(pred_,
                       pred_stride_,
                       residual_int32_,
                       residual_stride_,
                       recon_2_,
                       recon_stride_,
                       area_width_,
                       area_height_,
                       8);
        picture_addition_kernel(pred_,
                                pred_stride_,
                                residual_int32_,
                                residual_stride_,
                                recon_c_,
                                recon_stride_,
                                area_width_,
                                area_height_,
                                8);
        check_output(area_width_, area_height_, recon_1_, recon_c_);
        check_output(area_width_, area_height_, recon_2_, recon_c_);
        EXPECT_FALSE(HasFailure());
    }

    void run_16bit_test() {
        prepare_16bit_data();
        picture_addition_kernel16_bit(pred16bit_,
                                      pred_stride_,
                                      residual_int32_,
                                      residual_stride_,
                                      recon2_16bit_,
                                      recon_stride_,
                                      area_width_,
                                      area_height_,
                                      10);
        check_16bit_output(area_width_, area_height_);
        EXPECT_FALSE(HasFailure());
    }

    TestPattern test_pattern_;
    uint32_t area_width_, area_height_;
    uint32_t recon_stride_, pred_stride_, residual_stride_;
    uint8_t *pred_, *recon_1_, *recon_2_, *recon_c_;
    uint16_t *recon1_16bit_, *recon2_16bit_, *pred16bit_;
    int16_t *residual_;
    int32_t *residual_int32_;
    uint32_t test_size_;
    LbdPictureAdditionFunc lbd_test_func_;
    HbdPictureAdditionFunc hbd_test_func_;
};

TEST_P(PictureAdditionTest, PictureAdditionTest) {
    run_test();
};

TEST_P(PictureAdditionTest, PictureAddition16bitTest) {
    run_16bit_test();
};

}  // namespace
#endif
