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
#include "EbTime.h"

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
#define TEST_OFFSET 10

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

static const FwdTxfm2dFunc fwd_txfm_2d_N2_asm_func[TX_SIZES_ALL] = {
   svt_av1_fwd_txfm2d_4x4_N2_sse4_1, svt_av1_fwd_txfm2d_8x8_N2_avx2,
    svt_av1_fwd_txfm2d_16x16_N2_avx2, svt_av1_fwd_txfm2d_32x32_N2_avx2,
    svt_av1_fwd_txfm2d_64x64_N2_avx2, svt_av1_fwd_txfm2d_4x8_N2_avx2,
    svt_av1_fwd_txfm2d_8x4_N2_avx2,   svt_av1_fwd_txfm2d_8x16_N2_avx2,
    svt_av1_fwd_txfm2d_16x8_N2_avx2,  svt_av1_fwd_txfm2d_16x32_N2_avx2,
    svt_av1_fwd_txfm2d_32x16_N2_avx2, svt_av1_fwd_txfm2d_32x64_N2_avx2,
    svt_av1_fwd_txfm2d_64x32_N2_avx2, svt_av1_fwd_txfm2d_4x16_N2_avx2,
    svt_av1_fwd_txfm2d_16x4_N2_avx2,  svt_av1_fwd_txfm2d_8x32_N2_avx2,
    svt_av1_fwd_txfm2d_32x8_N2_avx2,  svt_av1_fwd_txfm2d_16x64_N2_avx2,
    svt_av1_fwd_txfm2d_64x16_N2_avx2,
};

static const FwdTxfm2dFunc fwd_txfm_2d_N4_asm_func[TX_SIZES_ALL] = {
    svt_av1_fwd_txfm2d_4x4_N4_sse4_1, svt_av1_fwd_txfm2d_8x8_N4_avx2,
    svt_av1_fwd_txfm2d_16x16_N4_avx2, svt_av1_fwd_txfm2d_32x32_N4_avx2,
    svt_av1_fwd_txfm2d_64x64_N4_avx2, svt_av1_fwd_txfm2d_4x8_N4_avx2,
    svt_av1_fwd_txfm2d_8x4_N4_avx2,   svt_av1_fwd_txfm2d_8x16_N4_avx2,
    svt_av1_fwd_txfm2d_16x8_N4_avx2,  svt_av1_fwd_txfm2d_16x32_N4_avx2,
    svt_av1_fwd_txfm2d_32x16_N4_avx2, svt_av1_fwd_txfm2d_32x64_N4_avx2,
    svt_av1_fwd_txfm2d_64x32_N4_avx2, svt_av1_fwd_txfm2d_4x16_N4_avx2,
    svt_av1_fwd_txfm2d_16x4_N4_avx2,  svt_av1_fwd_txfm2d_8x32_N4_avx2,
    svt_av1_fwd_txfm2d_32x8_N4_avx2,  svt_av1_fwd_txfm2d_16x64_N4_avx2,
    svt_av1_fwd_txfm2d_64x16_N4_avx2,
};

static const FwdTxfm2dFunc fwd_txfm_2d_N2_c_func[TX_SIZES_ALL] = {
    av1_transform_two_d_4x4_N2_c,   av1_transform_two_d_8x8_N2_c,
    av1_transform_two_d_16x16_N2_c, av1_transform_two_d_32x32_N2_c,
    av1_transform_two_d_64x64_N2_c, svt_av1_fwd_txfm2d_4x8_N2_c,
    svt_av1_fwd_txfm2d_8x4_N2_c,     svt_av1_fwd_txfm2d_8x16_N2_c,
    svt_av1_fwd_txfm2d_16x8_N2_c,    svt_av1_fwd_txfm2d_16x32_N2_c,
    svt_av1_fwd_txfm2d_32x16_N2_c,   svt_av1_fwd_txfm2d_32x64_N2_c,
    svt_av1_fwd_txfm2d_64x32_N2_c,   svt_av1_fwd_txfm2d_4x16_N2_c,
    svt_av1_fwd_txfm2d_16x4_N2_c,    svt_av1_fwd_txfm2d_8x32_N2_c,
    svt_av1_fwd_txfm2d_32x8_N2_c,    svt_av1_fwd_txfm2d_16x64_N2_c,
    svt_av1_fwd_txfm2d_64x16_N2_c,
};

static const FwdTxfm2dFunc fwd_txfm_2d_N4_c_func[TX_SIZES_ALL] = {
    av1_transform_two_d_4x4_N4_c,   av1_transform_two_d_8x8_N4_c,
    av1_transform_two_d_16x16_N4_c, av1_transform_two_d_32x32_N4_c,
    av1_transform_two_d_64x64_N4_c, svt_av1_fwd_txfm2d_4x8_N4_c,
    svt_av1_fwd_txfm2d_8x4_N4_c,     svt_av1_fwd_txfm2d_8x16_N4_c,
    svt_av1_fwd_txfm2d_16x8_N4_c,    svt_av1_fwd_txfm2d_16x32_N4_c,
    svt_av1_fwd_txfm2d_32x16_N4_c,   svt_av1_fwd_txfm2d_32x64_N4_c,
    svt_av1_fwd_txfm2d_64x32_N4_c,   svt_av1_fwd_txfm2d_4x16_N4_c,
    svt_av1_fwd_txfm2d_16x4_N4_c,    svt_av1_fwd_txfm2d_8x32_N4_c,
    svt_av1_fwd_txfm2d_32x8_N4_c,    svt_av1_fwd_txfm2d_16x64_N4_c,
    svt_av1_fwd_txfm2d_64x16_N4_c,
};

#ifndef NON_AVX512_SUPPORT
static const FwdTxfm2dFunc fwd_txfm_2d_N2_asm512_func[TX_SIZES_ALL] = {
    NULL,
    NULL,
    NULL,
    av1_fwd_txfm2d_32x32_N2_avx512,
    av1_fwd_txfm2d_64x64_N2_avx512,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    av1_fwd_txfm2d_32x64_N2_avx512,
    av1_fwd_txfm2d_64x32_N2_avx512,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
};

static const FwdTxfm2dFunc fwd_txfm_2d_N4_asm512_func[TX_SIZES_ALL] = {
    NULL, NULL, NULL, NULL, av1_fwd_txfm2d_64x64_N4_avx512,
    NULL, NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL,
};
#endif
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
        int over_buffer = sizeof(output_test_buf_) - width_ * height_;
        memset(output_test_buf_ + width_ * height_, 0xcd, over_buffer);
        memset(output_ref_buf_ + width_ * height_, 0xcd, over_buffer);
    }

    ~FwdTxfm2dAsmTest() {
        delete rnd_;
        aom_clear_system_state();
    }

    void run_match_test_default() {
        FwdTxfm2dFunc test_func = fwd_txfm_2d_asm_func[tx_size_];
        FwdTxfm2dFunc ref_func = fwd_txfm_2d_c_func[tx_size_];
        execute_test(test_func, ref_func, DEFAULT_SHAPE);
    }

    void run_match_test_N2() {
        FwdTxfm2dFunc test_func_asm = fwd_txfm_2d_N2_asm_func[tx_size_];
        FwdTxfm2dFunc test_func_c = fwd_txfm_2d_N2_c_func[tx_size_];
        FwdTxfm2dFunc ref_func = fwd_txfm_2d_c_func[tx_size_];
        execute_test(test_func_asm, ref_func, N2_SHAPE);
        execute_test(test_func_c, ref_func, N2_SHAPE);
    }

    void run_match_test_N4() {
        FwdTxfm2dFunc test_func_asm = fwd_txfm_2d_N4_asm_func[tx_size_];
        FwdTxfm2dFunc test_func_c = fwd_txfm_2d_N4_c_func[tx_size_];
        FwdTxfm2dFunc ref_func = fwd_txfm_2d_c_func[tx_size_];
        execute_test(test_func_asm, ref_func, N4_SHAPE);
        execute_test(test_func_c, ref_func, N4_SHAPE);
    }

    void speed_test() {
        FwdTxfm2dFunc test_func = fwd_txfm_2d_asm_func[tx_size_];
        FwdTxfm2dFunc ref_func = fwd_txfm_2d_c_func[tx_size_];
        run_speed_test("C   AND ASM", test_func, ref_func);
        run_speed_test(
            "ASM AND N2 ", fwd_txfm_2d_N2_asm_func[tx_size_], test_func);
        run_speed_test(
            "ASM AND N4 ", fwd_txfm_2d_N4_asm_func[tx_size_], test_func);
    }

#ifndef NON_AVX512_SUPPORT
    void run_match_test_N2_512() {
        FwdTxfm2dFunc test_func_asm = fwd_txfm_2d_N2_asm512_func[tx_size_];
        FwdTxfm2dFunc ref_func = fwd_txfm_2d_c_func[tx_size_];
        execute_test(test_func_asm, ref_func, N2_SHAPE);
    }

    void run_match_test_N4_512() {
        FwdTxfm2dFunc test_func_asm = fwd_txfm_2d_N4_asm512_func[tx_size_];
        FwdTxfm2dFunc ref_func = fwd_txfm_2d_c_func[tx_size_];
        execute_test(test_func_asm, ref_func, N4_SHAPE);
    }

    void speed_test_512() {
        FwdTxfm2dFunc test_func = fwd_txfm_2d_N2_asm512_func[tx_size_];
        FwdTxfm2dFunc ref_func = fwd_txfm_2d_N2_asm_func[tx_size_];
        if (test_func && ref_func) {
            run_speed_test("N2 AVX512 AVX2 ", test_func, ref_func);
        }

        test_func = fwd_txfm_2d_N4_asm512_func[tx_size_];
        ref_func = fwd_txfm_2d_N4_asm_func[tx_size_];
        if (test_func && ref_func) {
            run_speed_test("N4 AVX512 AVX2", test_func, ref_func);
        }
    }
#endif

  private:

       void execute_test(FwdTxfm2dFunc test_func, FwdTxfm2dFunc ref_func,
                      EB_TRANS_COEFF_SHAPE shape) {
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
                if (shape == N2_SHAPE) {
                    for (int i = 0;
                         i < (tx_size_wide[tx_size_] * tx_size_high[tx_size_]);
                         i++) {
                        if (i % tx_size_wide[tx_size_] >=
                                (tx_size_wide[tx_size_] >> 1) ||
                            i / tx_size_wide[tx_size_] >=
                                (tx_size_high[tx_size_] >> 1)) {
                            output_ref_[i] = 0;
                        }
                    }
                } else if (shape == N4_SHAPE) {
                    for (int i = 0;
                         i < (tx_size_wide[tx_size_] * tx_size_high[tx_size_]);
                         i++) {
                        if (i % tx_size_wide[tx_size_] >=
                                (tx_size_wide[tx_size_] >> 2) ||
                            i / tx_size_wide[tx_size_] >=
                                (tx_size_high[tx_size_] >> 2)) {
                            output_ref_[i] = 0;
                        }
                    }
                }
                test_func(input_, output_test_, stride_, type, (uint8_t)bd_);

                if (0 != memcmp(output_test_, output_ref_,
                           MAX_TX_SQUARE * sizeof(int32_t) + TEST_OFFSET)) {
                    for (int i = 0; i < height_; i++)
                        for (int j = 0; j < width_; j++) {
                            if (output_ref_[i * width_ + j] !=
                                output_test_[i * width_ + j]) {
                                printf("error in important part\n");
                            }

                            ASSERT_EQ(output_ref_[i * width_ + j],
                                      output_test_[i * width_ + j])
                                << "loop: " << k << " tx_type: " << tx_type
                                << " tx_size: " << tx_size_ << " Mismatch at ("
                                << j << " x " << i << ")";
                        }

                    ASSERT_EQ(1, 0);
                }
            }
        }
    }

    void run_speed_test(char *name_cmp, FwdTxfm2dFunc test_func,
                        FwdTxfm2dFunc ref_func) {
        double time_c, time_o;
        uint64_t start_time_seconds, start_time_useconds;
        uint64_t middle_time_seconds, middle_time_useconds;
        uint64_t finish_time_seconds, finish_time_useconds;
        const char *tx_type_name[] = {"DCT_DCT",
                                      "ADST_DCT",
                                      "DCT_ADST",
                                      "ADST_ADST",
                                      "FLIPADST_DCT",
                                      "DCT_FLIPADST",
                                      "FLIPADST_FLIPADST",
                                      "ADST_FLIPADST",
                                      "FLIPADST_ADST",
                                      "IDTX",
                                      "V_DCT",
                                      "H_DCT",
                                      "V_ADST",
                                      "H_ADST",
                                      "V_FLIPADST",
                                      "H_FLIPADST",
                                      "TX_TYPES"};

        if (ref_func == nullptr || test_func == nullptr)
            return;

        ASSERT_NE(rnd_, nullptr) << "Failed to create random generator";
        for (int tx_type = 0; tx_type < TX_TYPES; ++tx_type) {
            TxType type = static_cast<TxType>(tx_type);
            populate_with_random();
            // tx_type and tx_size are not compatible in the av1-spec.
            // like the max size of adst transform is 16, and max size of
            // identity transform is 32.
            if (is_txfm_allowed(type, tx_size_) == false)
                continue;

            const int loops = 500000;
            svt_av1_get_time(&start_time_seconds, &start_time_useconds);
            for (int k = 0; k < loops; k++) {
                ref_func(input_, output_ref_, stride_, type, (uint8_t)bd_);
            }
            svt_av1_get_time(&middle_time_seconds, &middle_time_useconds);
            for (int k = 0; k < loops; k++) {
                test_func(input_, output_test_, stride_, type, (uint8_t)bd_);
            }
            svt_av1_get_time(&finish_time_seconds, &finish_time_useconds);

            time_c = svt_av1_compute_overall_elapsed_time_ms(start_time_seconds,
                                                             start_time_useconds,
                                                             middle_time_seconds,
                                                             middle_time_useconds);

            time_o = svt_av1_compute_overall_elapsed_time_ms(middle_time_seconds,
                                                             middle_time_useconds,
                                                             finish_time_seconds,
                                                             finish_time_useconds);

            printf("[%s]; Transform: ;%02ix%02i; %17s; Speed compare: ;%5.2fx\n",
                    name_cmp,
                    tx_size_wide[tx_size_],
                    tx_size_high[tx_size_],
                    tx_type_name[tx_type],
                    time_c / time_o);
        }
    }
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
    uint8_t output_test_buf_[MAX_TX_SQUARE * sizeof(int32_t) + ALIGNMENT - 1 +
                             TEST_OFFSET];
    uint8_t output_ref_buf_[MAX_TX_SQUARE * sizeof(int32_t) + ALIGNMENT - 1 +
                            TEST_OFFSET];
    int16_t *input_;       /**< aligned address for input */
    int32_t *output_test_; /**< aligned address for output test */
    int32_t *output_ref_;  /**< aligned address for output ref */
};

TEST_P(FwdTxfm2dAsmTest, match_test) {
    run_match_test_default();
}

TEST_P(FwdTxfm2dAsmTest, match_test_N2) {
    run_match_test_N2();
}

TEST_P(FwdTxfm2dAsmTest, match_test_N4) {
    run_match_test_N4();
}

TEST_P(FwdTxfm2dAsmTest, DISABLED_speed_test) {
    speed_test();
}

#ifndef NON_AVX512_SUPPORT
TEST_P(FwdTxfm2dAsmTest, match_test_N2_512) {
    if (CPU_FLAGS_AVX512F & get_cpu_flags_to_use()) {
        run_match_test_N2_512();
    }
}

TEST_P(FwdTxfm2dAsmTest, match_test_N4_512) {
    if (CPU_FLAGS_AVX512F & get_cpu_flags_to_use()) {
        run_match_test_N4_512();
    }
}

TEST_P(FwdTxfm2dAsmTest, DISABLED_speed_test_512) {
    if (CPU_FLAGS_AVX512F & get_cpu_flags_to_use()) {
        speed_test_512();
    }
}
#endif
INSTANTIATE_TEST_CASE_P(
    TX, FwdTxfm2dAsmTest,
    ::testing::Combine(::testing::Range(static_cast<int>(TX_4X4),
                                        static_cast<int>(TX_SIZES_ALL), 1),
                       ::testing::Values(static_cast<int>(AOM_BITS_8),
                                         static_cast<int>(AOM_BITS_10))));
}  // namespace
