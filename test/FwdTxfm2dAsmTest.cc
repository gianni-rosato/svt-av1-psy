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
 * @file FwdTxfm2dAsmTest.c
 *
 * @brief Unit test for forward 2d transform functions written in assembly code:
 * - svt_av1_fwd_txfm2d_{4, 8, 16, 32, 64}x{4, 8, 16, 32, 64}_avx2
 *
 * @author Cidana-Wenyao
 *
 ******************************************************************************/
#include "gtest/gtest.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <new>
#include <algorithm>

// workaround to eliminate the compiling warning on linux
// The macro will conflict with definition in gtest.h
#ifdef __USE_GNU
#undef __USE_GNU  // defined in EbThreads.h
#endif
#ifdef _GNU_SOURCE
#undef _GNU_SOURCE  // defined in EbThreads.h
#endif
#include "EbDefinitions.h"
#include "EbTransforms.h"
#include "random.h"
#include "util.h"
#include "aom_dsp_rtcd.h"
#include "EbTransforms.h"
#include "TxfmCommon.h"

using svt_av1_test_tool::SVTRandom;
namespace {
using FwdTxfm2dAsmParam = std::tuple<int, int>;
static const FwdTxfm2dFunc fwd_txfm_2d_asm_func[TX_SIZES_ALL] = {
    svt_av1_fwd_txfm2d_4x4_sse4_1, svt_av1_fwd_txfm2d_8x8_avx2,
    svt_av1_fwd_txfm2d_16x16_avx2, svt_av1_fwd_txfm2d_32x32_avx2,
    svt_av1_fwd_txfm2d_64x64_avx2, svt_av1_fwd_txfm2d_4x8_avx2,
    svt_av1_fwd_txfm2d_8x4_avx2,   svt_av1_fwd_txfm2d_8x16_avx2,
    svt_av1_fwd_txfm2d_16x8_avx2,  svt_av1_fwd_txfm2d_16x32_avx2,
    svt_av1_fwd_txfm2d_32x16_avx2, svt_av1_fwd_txfm2d_32x64_avx2,
    svt_av1_fwd_txfm2d_64x32_avx2, svt_av1_fwd_txfm2d_4x16_avx2,
    svt_av1_fwd_txfm2d_16x4_avx2,  svt_av1_fwd_txfm2d_8x32_avx2,
    svt_av1_fwd_txfm2d_32x8_avx2,  svt_av1_fwd_txfm2d_16x64_avx2,
    svt_av1_fwd_txfm2d_64x16_avx2,
};

/**
 * @brief Unit test for fwd tx 2d avx2 functions:
 * - svt_av1_fwd_txfm2d_{4, 8, 16, 32, 64}x{4, 8, 16, 32, 64}_avx2
 *
 * Test strategy:
 * Verify this assembly code by comparing with reference c implementation.
 * Feed the same data and check test output and reference output.
 * The test output and reference output are different at the beginning.
 *
 * Expect result:
 * Output from assemble function should be exactly same as output from c.
 *
 * Test coverage:
 * Test cases:
 * Input buffer: Fill with random values
 * TxSize: all the valid TxSize and TxType allowed.
 * BitDepth: 8bit and 10bit.
 *
 */
class FwdTxfm2dAsmTest : public ::testing::TestWithParam<FwdTxfm2dAsmParam> {
  public:
    FwdTxfm2dAsmTest()
        : tx_size_(static_cast<TxSize>(TEST_GET_PARAM(0))),
          bd_(TEST_GET_PARAM(1)) {
        // input are signed value with bitdepth + 1 bits
        rnd_ = new SVTRandom(-(1 << bd_) + 1, (1 << bd_) - 1);

        width_ = tx_size_wide[tx_size_];
        height_ = tx_size_high[tx_size_];
        // set default value to 0
        memset(output_test_buf_, 0, sizeof(output_test_buf_));
        // set default value to -1
        memset(output_ref_buf_, 255, sizeof(output_ref_buf_));
        input_ = ALIGNED_ADDR(int16_t, ALIGNMENT, input_buf_);
        output_test_ = ALIGNED_ADDR(int32_t, ALIGNMENT, output_test_buf_);
        output_ref_ = ALIGNED_ADDR(int32_t, ALIGNMENT, output_ref_buf_);
    }

    ~FwdTxfm2dAsmTest() {
        delete rnd_;
        aom_clear_system_state();
    }

    void run_match_test() {
        FwdTxfm2dFunc test_func = fwd_txfm_2d_asm_func[tx_size_];
        FwdTxfm2dFunc ref_func = fwd_txfm_2d_c_func[tx_size_];
        if (ref_func == nullptr || test_func == nullptr)
            return;

        ASSERT_NE(rnd_, nullptr) << "Failed to create random generator";
        for (int tx_type = 0; tx_type < TX_TYPES; ++tx_type) {
            TxType type = static_cast<TxType>(tx_type);
            // tx_type and tx_size are not compatible in the av1-spec.
            // like the max size of adst transform is 16, and max size of
            // identity transform is 32.
            if (is_txfm_allowed(type, tx_size_) == false)
                continue;

            const int loops = 100;
            for (int k = 0; k < loops; k++) {
                populate_with_random();

                ref_func(input_, output_ref_, stride_, type, (uint8_t)bd_);
                test_func(input_, output_test_, stride_, type, (uint8_t)bd_);

                for (int i = 0; i < height_; i++)
                    for (int j = 0; j < width_; j++)
                        ASSERT_EQ(output_ref_[i * width_ + j],
                                  output_test_[i * width_ + j])
                            << "loop: " << k << " tx_type: " << tx_type
                            << " tx_size: " << tx_size_ << " Mismatch at (" << j
                            << " x " << i << ")";
            }
        }
    }

  private:
    void populate_with_random() {
        for (int i = 0; i < height_; i++) {
            for (int j = 0; j < width_; j++) {
                input_[i * stride_ + j] = (int16_t)rnd_->random();
            }
        }

        return;
    }

  private:
    const TxSize tx_size_; /**< input param tx_size */
    const int bd_;         /**< input param 8bit or 10bit */
    int width_;
    int height_;
    SVTRandom *rnd_;
    static const int stride_ = MAX_TX_SIZE;
    uint8_t input_buf_[MAX_TX_SQUARE * sizeof(int16_t) + ALIGNMENT - 1];
    uint8_t output_test_buf_[MAX_TX_SQUARE * sizeof(int32_t) + ALIGNMENT - 1];
    uint8_t output_ref_buf_[MAX_TX_SQUARE * sizeof(int32_t) + ALIGNMENT - 1];
    int16_t *input_;       /**< aligned address for input */
    int32_t *output_test_; /**< aligned address for output test */
    int32_t *output_ref_;  /**< aligned address for output ref */
};

TEST_P(FwdTxfm2dAsmTest, match_test) {
    run_match_test();
}

INSTANTIATE_TEST_CASE_P(
    TX, FwdTxfm2dAsmTest,
    ::testing::Combine(::testing::Range(static_cast<int>(TX_4X4),
                                        static_cast<int>(TX_SIZES_ALL), 1),
                       ::testing::Values(static_cast<int>(AOM_BITS_8),
                                         static_cast<int>(AOM_BITS_10))));
}  // namespace
