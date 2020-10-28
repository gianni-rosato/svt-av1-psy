/******************************************************************************
 * @file EbSatdTest.cc
 *

 *
 * @author tszumski@intel.com
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

#include "EbEncIntraPrediction.h"
#include "EbDefinitions.h"
#include "random.h"
#include "util.h"


using svt_av1_test_tool::SVTRandom;  // to generate the random

namespace {

typedef enum { MAX, MIN, RANDOM } TestPattern;

class Compute8x8SatdTest : public ::testing::Test {
  public:
    Compute8x8SatdTest() : _width_(8), _height_(8), _dcSize_(10) {
        _stride_ = _width_ + 4;
        _dcValue1_ = (uint64_t *)svt_aom_malloc(_dcSize_ * sizeof(uint64_t));
        _dcValue2_ = (uint64_t *)svt_aom_malloc(_dcSize_ * sizeof(uint64_t));
        _src_ = (uint8_t *)svt_aom_malloc(_height_ * _stride_ * sizeof(uint8_t));
    }

    void TearDown() override {
        if (_dcValue1_)
            svt_aom_free(_dcValue1_);
        if (_dcValue2_)
            svt_aom_free(_dcValue2_);
        if (_src_)
            svt_aom_free(_src_);
    }

  protected:
    void prepare_data(TestPattern pattern) {
        const int32_t mask = (1 << 8) - 1;
        SVTRandom rnd(0, mask);
        switch (pattern) {
        case MAX: {
            memset(_src_, mask, _height_ * _stride_ * sizeof(uint8_t));
            for (uint32_t i = 0; i < _dcSize_; i++) {
                _dcValue1_[i] = _dcValue2_[i] = mask;
            }
            break;
        }
        case MIN: {
            memset(_src_, 0, _stride_ * sizeof(uint8_t));
            memset(_dcValue1_, 0, _dcSize_ * sizeof(uint64_t));
            memset(_dcValue2_, 0, _dcSize_ * sizeof(uint64_t));
            break;
        }
        case RANDOM: {
            for (uint32_t i = 0; i < _height_ * _stride_; i++) {
                _src_[i] = rnd.random();
            }

            for (uint32_t i = 0; i < _dcSize_; i++) {
                _dcValue1_[i] = _dcValue2_[i] = rnd.random();
            }

            break;
        }
        default: break;
        }
    }

    void run_min_test() {
        prepare_data(MIN);
        int fail_pixel_count = 0;
        for (uint16_t j = 0; j < _dcSize_; j++) {
            if (_dcValue1_[j] != _dcValue2_[j])
                fail_pixel_count++;
        }
        EXPECT_EQ(0, fail_pixel_count)
            << "error in dcValue for " << fail_pixel_count << "times";
    }

    void run_max_test() {
        prepare_data(MAX);
        int fail_pixel_count = 0;
        for (uint16_t j = 0; j < _dcSize_; j++) {
            if (_dcValue1_[j] != _dcValue2_[j])
                fail_pixel_count++;
        }
        EXPECT_EQ(0, fail_pixel_count)
            << "error in dcValue for " << fail_pixel_count << "times";
    }

    void run_rnd_test() {
        for (int i = 0; i < 500; i++) {
            prepare_data(RANDOM);
            int fail_pixel_count = 0;
            for (uint16_t j = 0; j < _dcSize_; j++) {
                if (_dcValue1_[j] != _dcValue2_[j])
                    fail_pixel_count++;
            }
            EXPECT_EQ(0, fail_pixel_count)
                << "error in dcValue for " << fail_pixel_count
                << "times";
        }
    }

    uint32_t _width_, _height_;
    uint32_t _stride_;
    uint32_t _dcSize_;
    uint64_t *_dcValue1_, *_dcValue2_;
    uint8_t *_src_;
};

TEST_F(Compute8x8SatdTest, minTest) {
    run_min_test();
};

TEST_F(Compute8x8SatdTest, maxTest) {
    run_max_test();
};

TEST_F(Compute8x8SatdTest, rndTest) {
    run_rnd_test();
};

}  // namespace
