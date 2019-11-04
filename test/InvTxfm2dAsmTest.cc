/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file InvTxfm2dAsmTest.c
 *
 * @brief Unit test for forward 2d transform functions:
 * - Av1TransformTwoD_{4x4, 8x8, 16x16, 32x32, 64x64}
 * - eb_av1_fwd_txfm2d_{rectangle}
 *
 * @author Cidana-Wenyao
 *
 ******************************************************************************/
#include "gtest/gtest.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
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
#include "EbUnitTestUtility.h"
#include "TxfmCommon.h"
#include "av1_inv_txfm_ssse3.h"

using svt_av1_test_tool::SVTRandom;  // to generate the random
namespace {

using InvSqrTxfm2dFun = void (*)(const int32_t *input,
                                 uint16_t *output_r, int32_t stride_r,
                                 uint16_t *output_w, int32_t stride_w,
                                 TxType tx_type, int32_t bd);
using InvRectTxfm2dType1Func = void (*)(const int32_t *input,
                                        uint16_t *output_r, int32_t stride_r,
                                        uint16_t *output_w, int32_t stride_w,
                                        TxType tx_type, TxSize tx_size,
                                        int32_t eob, int32_t bd);
using InvRectTxfm2dType2Func = void (*)(const int32_t *input,
                                        uint16_t *output_r, int32_t stride_r,
                                        uint16_t *output_w, int32_t stride_w,
                                        TxType tx_type, TxSize tx_size,
                                        int32_t bd);
using LowbdInvTxfm2dFunc = void (*)(const int32_t *input,
                                    uint8_t *output_r, int32_t stride_r,
                                    uint8_t *output_w, int32_t stride_w,
                                    TxType tx_type, TxSize tx_size,
                                    int32_t eob);
typedef struct {
    const char *name;
    InvSqrTxfm2dFun ref_func;
    InvSqrTxfm2dFun test_func;
    IsTxTypeImpFunc check_imp_func;
} InvSqrTxfmFuncPair;

typedef struct {
    InvRectTxfm2dType2Func ref_func;
    InvRectTxfm2dType2Func test_func;
} InvRectType2TxfmFuncPair;

#define SQR_FUNC_PAIRS(name, type, is_tx_type_imp)                \
    {                                                             \
        #name, reinterpret_cast < InvSqrTxfm2dFun > (name##_c),   \
            reinterpret_cast < InvSqrTxfm2dFun > (name##_##type), \
            is_tx_type_imp                                        \
    }

#define EMPTY_FUNC_PAIRS(name) \
    { #name, nullptr, nullptr, nullptr }

static bool is_tx_type_imp_32x32_avx2(const TxType tx_type) {
    switch (tx_type) {
    case DCT_DCT:
    case IDTX: return true;
    default: return false;
    }
}

static bool is_tx_type_imp_64x64_sse4(const TxType tx_type) {
    if (tx_type == DCT_DCT)
        return true;
    return false;
}

static const InvSqrTxfmFuncPair inv_txfm_c_avx2_func_pairs[TX_64X64 + 1] = {
    SQR_FUNC_PAIRS(eb_av1_inv_txfm2d_add_4x4, avx2, all_txtype_imp),
    SQR_FUNC_PAIRS(eb_av1_inv_txfm2d_add_8x8, avx2, all_txtype_imp),
    SQR_FUNC_PAIRS(eb_av1_inv_txfm2d_add_16x16, avx2, all_txtype_imp),
    SQR_FUNC_PAIRS(eb_av1_inv_txfm2d_add_32x32, avx2, is_tx_type_imp_32x32_avx2),
    EMPTY_FUNC_PAIRS(eb_av1_inv_txfm2d_add_64x64),
};

static const InvSqrTxfmFuncPair inv_txfm_c_sse4_1_func_pairs[TX_64X64 + 1] = {
    SQR_FUNC_PAIRS(eb_av1_inv_txfm2d_add_4x4, sse4_1, dct_adst_combine_imp),
    SQR_FUNC_PAIRS(eb_av1_inv_txfm2d_add_8x8, sse4_1, dct_adst_combine_imp),
    SQR_FUNC_PAIRS(eb_av1_inv_txfm2d_add_16x16, sse4_1, dct_adst_combine_imp),
    EMPTY_FUNC_PAIRS(eb_av1_inv_txfm2d_add_32x32),
    SQR_FUNC_PAIRS(eb_av1_inv_txfm2d_add_64x64, sse4_1, is_tx_type_imp_64x64_sse4),
};

// from TX_4X8 to TX_SIZES_ALL
static const InvRectTxfm2dType1Func rect_type1_ref_funcs[TX_SIZES_ALL] = {
    // square transform
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,  // 4x8 and 8x4
    eb_av1_inv_txfm2d_add_8x16_c,
    eb_av1_inv_txfm2d_add_16x8_c,
    eb_av1_inv_txfm2d_add_16x32_c,
    eb_av1_inv_txfm2d_add_32x16_c,
    eb_av1_inv_txfm2d_add_32x64_c,
    eb_av1_inv_txfm2d_add_64x32_c,
    nullptr,
    nullptr,  // 4x16 and 16x4
    eb_av1_inv_txfm2d_add_8x32_c,
    eb_av1_inv_txfm2d_add_32x8_c,
    eb_av1_inv_txfm2d_add_16x64_c,
    eb_av1_inv_txfm2d_add_64x16_c};

static const InvRectType2TxfmFuncPair inv_4x8{eb_av1_inv_txfm2d_add_4x8_c,
                                              eb_av1_inv_txfm2d_add_4x8_sse4_1};
static const InvRectType2TxfmFuncPair inv_8x4{eb_av1_inv_txfm2d_add_8x4_c,
                                              eb_av1_inv_txfm2d_add_8x4_sse4_1};
static const InvRectType2TxfmFuncPair inv_4x16{eb_av1_inv_txfm2d_add_4x16_c,
                                               eb_av1_inv_txfm2d_add_4x16_sse4_1};
static const InvRectType2TxfmFuncPair inv_16x4{eb_av1_inv_txfm2d_add_16x4_c,
                                               eb_av1_inv_txfm2d_add_16x4_sse4_1};
static const InvRectType2TxfmFuncPair *get_rect_type2_func_pair(
    const TxSize tx_size) {
    switch (tx_size) {
    case TX_4X8: return &inv_4x8;
    case TX_8X4: return &inv_8x4;
    case TX_4X16: return &inv_4x16;
    case TX_16X4: return &inv_16x4;
    default: return nullptr;
    }
}

/**
 * @brief Unit test for inverse tx 2d avx2/sse4_1 functions:
 * - eb_av1_inv_txfm2d_{4, 8, 16, 32, 64}x{4, 8, 16, 32, 64}_avx2
 *
 * Test strategy:
 * Verify this assembly code by comparing with reference c implementation.
 * Feed the same data and check test output and reference output. Four tests
 * are required since there are three different function signatures and one
 * set of function for lowbd functions.
 *
 * Expect result:
 * Output from assemble function should be exactly same as output from c.
 *
 * Test coverage:
 * Test cases:
 * Input buffer: Fill with random values
 * TxSize: all the valid TxSize and TxType allowed.
 * BitDepth: 8bit and 10bit
 * AssembleType: avx2 and sse4_1
 *
 */

using InvTxfm2dParam = std::tuple<LowbdInvTxfm2dFunc, int>;

class InvTxfm2dAsmTest : public ::testing::TestWithParam<InvTxfm2dParam> {
  public:
    InvTxfm2dAsmTest()
        : bd_(TEST_GET_PARAM(1)), target_func_(TEST_GET_PARAM(0)) {
        // unsigned bd_ bits random
        u_bd_rnd_ = new SVTRandom(0, (1 << bd_) - 1);
        s_bd_rnd_ = new SVTRandom(-(1 << bd_) + 1, (1 << bd_) - 1);
    }

    ~InvTxfm2dAsmTest() {
        delete u_bd_rnd_;
        delete s_bd_rnd_;
        aom_clear_system_state();
    }

    void SetUp() override {
        pixel_input_ = reinterpret_cast<int16_t *>(
            eb_aom_memalign(32, MAX_TX_SQUARE * sizeof(int16_t)));
        input_ = reinterpret_cast<int32_t *>(
            eb_aom_memalign(32, MAX_TX_SQUARE * sizeof(int32_t)));
        output_test_ = reinterpret_cast<uint16_t *>(
            eb_aom_memalign(32, MAX_TX_SQUARE * sizeof(uint16_t)));
        output_ref_ = reinterpret_cast<uint16_t *>(
            eb_aom_memalign(32, MAX_TX_SQUARE * sizeof(uint16_t)));
        lowbd_output_test_ = reinterpret_cast<uint8_t *>(
            eb_aom_memalign(32, MAX_TX_SQUARE * sizeof(uint8_t)));
    }

    void TearDown() override {
        eb_aom_free(pixel_input_);
        eb_aom_free(input_);
        eb_aom_free(output_test_);
        eb_aom_free(output_ref_);
        eb_aom_free(lowbd_output_test_);
        aom_clear_system_state();
    }

    void run_sqr_txfm_match_test(const TxSize tx_size, int is_asm_kernel) {
        const int width = tx_size_wide[tx_size];
        const int height = tx_size_high[tx_size];
        InvSqrTxfmFuncPair pair = (is_asm_kernel == 0)
                                      ? inv_txfm_c_avx2_func_pairs[tx_size]
                                      : inv_txfm_c_sse4_1_func_pairs[tx_size];
        if (pair.ref_func == nullptr || pair.test_func == nullptr)
            return;
        for (int tx_type = DCT_DCT; tx_type < TX_TYPES; ++tx_type) {
            TxType type = static_cast<TxType>(tx_type);
            const IsTxTypeImpFunc is_tx_type_imp = pair.check_imp_func;

            if (is_txfm_allowed(type, tx_size) == false)
                continue;

            // Some tx_type is not implemented yet, so we will skip this;
            if (is_tx_type_imp(type) == false)
                continue;

            const int loops = 100;
            for (int k = 0; k < loops; k++) {
                populate_with_random(width, height, type, tx_size);

                pair.ref_func(input_,
                              output_ref_, stride_,
                              output_ref_, stride_,
                              type, bd_);
                pair.test_func(input_,
                               output_test_, stride_,
                               output_test_, stride_,
                               type,
                               bd_);

                EXPECT_EQ(0,
                          memcmp(output_ref_,
                                 output_test_,
                                 height * stride_ * sizeof(output_test_[0])))
                    << "loop: " << k << " tx_type: " << tx_type
                    << " tx_size: " << tx_size
                    << " is_asm_kernel: " << is_asm_kernel;
            }
        }
    }

    void run_rect_type1_txfm_match_test(const TxSize tx_size) {
        const int width = tx_size_wide[tx_size];
        const int height = tx_size_high[tx_size];
        const int max_eob = av1_get_max_eob(tx_size);

        const InvRectTxfm2dType1Func test_func = eb_av1_highbd_inv_txfm_add_avx2;
        const InvRectTxfm2dType1Func ref_func = rect_type1_ref_funcs[tx_size];
        if (ref_func == nullptr)
            return;

        for (int tx_type = DCT_DCT; tx_type < TX_TYPES; ++tx_type) {
            TxType type = static_cast<TxType>(tx_type);

            if (is_txfm_allowed(type, tx_size) == false)
                continue;

            const int loops = 10 * max_eob;
            SVTRandom eob_rnd(1, max_eob - 1);
            for (int k = 0; k < loops; k++) {
                int eob = k < max_eob - 1 ? k + 1 : eob_rnd.random();
                // prepare data by forward transform and then
                // clear the values between eob and max_eob
                populate_with_random(width, height, type, tx_size);
                clear_high_freq_coeffs(tx_size, type, eob, max_eob);

                ref_func(input_,
                    output_ref_, stride_, output_ref_, stride_,
                    type, tx_size, eob, bd_);
                test_func(input_,
                          output_test_, stride_, output_test_,
                          stride_, type, tx_size, eob, bd_);

                ASSERT_EQ(0,
                          memcmp(output_ref_,
                                 output_test_,
                                 height * stride_ * sizeof(output_test_[0])))
                    << "loop: " << k << " tx_type: " << tx_type
                    << " tx_size: " << tx_size << " eob: " << eob;
            }
        }
    }

    void run_rect_type2_txfm_match_test(const TxSize tx_size) {
        const int width = tx_size_wide[tx_size];
        const int height = tx_size_high[tx_size];
        const InvRectType2TxfmFuncPair *test_pair =
            get_rect_type2_func_pair(tx_size);
        if (test_pair == nullptr)
            return;

        for (int tx_type = DCT_DCT; tx_type < TX_TYPES; ++tx_type) {
            TxType type = static_cast<TxType>(tx_type);

            if (is_txfm_allowed(type, tx_size) == false)
                continue;

            const int loops = 100;
            for (int k = 0; k < loops; k++) {
                populate_with_random(width, height, type, tx_size);

                test_pair->ref_func(input_,
                                    output_ref_, stride_,
                                    output_ref_, stride_,
                                    type, tx_size, bd_);
                test_pair->test_func(input_,
                                     output_test_, stride_,
                                     output_test_, stride_,
                                     type, tx_size, bd_);

                ASSERT_EQ(0,
                          memcmp(output_ref_,
                                 output_test_,
                                 height * stride_ * sizeof(output_test_[0])))
                    << "loop: " << k << " tx_type: " << tx_type
                    << " tx_size: " << tx_size;
            }
        }
    }

    void run_lowbd_txfm_match_test(const TxSize tx_size) {
        if (bd_ > 8)
            return;
        const int width = tx_size_wide[tx_size];
        const int height = tx_size_high[tx_size];
        const int max_eob = av1_get_max_eob(tx_size);
        using LowbdInvRectTxfmRefFunc = void (*)(const int32_t *input,
                                                 uint16_t *output_r,
                                                 int32_t stride_r,
                                                 uint16_t *output_w,
                                                 int32_t stride_w,
                                                 TxType tx_type,
                                                 TxSize tx_size,
                                                 int32_t eob,
                                                 int32_t bd);
        using LowbdInvSqrTxfmRefFunc = void (*)(const int32_t *input,
                                                uint16_t *output_r,
                                                int32_t stride_r,
                                                uint16_t *output_w,
                                                int32_t stride_w,
                                                TxType tx_type,
                                                int32_t bd);
        const LowbdInvSqrTxfmRefFunc lowbd_sqr_ref_funcs[TX_SIZES] = {
            eb_av1_inv_txfm2d_add_4x4_c,
            eb_av1_inv_txfm2d_add_8x8_c,
            eb_av1_inv_txfm2d_add_16x16_c,
            eb_av1_inv_txfm2d_add_32x32_c,
            eb_av1_inv_txfm2d_add_64x64_c};
        const LowbdInvRectTxfmRefFunc lowbd_rect_ref_funcs[TX_SIZES_ALL] = {
            nullptr,
            nullptr,
            nullptr,
            nullptr,
            nullptr,
            nullptr,
            nullptr,
            eb_av1_inv_txfm2d_add_8x16_c,
            eb_av1_inv_txfm2d_add_16x8_c,
            eb_av1_inv_txfm2d_add_16x32_c,
            eb_av1_inv_txfm2d_add_32x16_c,
            eb_av1_inv_txfm2d_add_32x64_c,
            eb_av1_inv_txfm2d_add_64x32_c,
            nullptr,
            nullptr,
            eb_av1_inv_txfm2d_add_8x32_c,
            eb_av1_inv_txfm2d_add_32x8_c,
            eb_av1_inv_txfm2d_add_16x64_c,
            eb_av1_inv_txfm2d_add_64x16_c};

        if (tx_size >= TX_SIZES && lowbd_rect_ref_funcs[tx_size] == nullptr)
            return;

        for (int tx_type = DCT_DCT; tx_type < TX_TYPES; ++tx_type) {
            TxType type = static_cast<TxType>(tx_type);

            if (is_txfm_allowed(type, tx_size) == false)
                continue;

            const int loops = 10 * max_eob;
            SVTRandom eob_rnd(1, max_eob - 1);
            for (int k = 0; k < loops; k++) {
                int eob = k < max_eob - 1 ? k + 1 : eob_rnd.random();
                // prepare data by forward transform and then
                // clear the values between eob and max_eob
                populate_with_random(width, height, type, tx_size);
                clear_high_freq_coeffs(tx_size, type, eob, max_eob);
                // copy to lowbd output buffer from short buffer
                for (int i = 0; i < height; i++) {
                    for (int j = 0; j < width; j++)
                        lowbd_output_test_[i * stride_ + j] =
                            static_cast<uint8_t>(output_test_[i * stride_ + j]);
                }

                target_func_(input_,
                             lowbd_output_test_, stride_,
                             lowbd_output_test_, stride_,
                             type,
                             tx_size,
                             eob);
                if (tx_size >= TX_SIZES)
                    lowbd_rect_ref_funcs[tx_size](
                        input_,
                        output_ref_, stride_,
                        output_ref_, stride_,
                        type, tx_size, eob, bd_);
                else
                    lowbd_sqr_ref_funcs[tx_size](
                        input_,
                        output_ref_, stride_,
                        output_ref_, stride_,
                        type, bd_);

                // compare, note the output buffer has stride.
                for (int i = 0; i < height; i++) {
                    for (int j = 0; j < width; j++) {
                        uint8_t ref =
                            static_cast<uint8_t>(output_ref_[i * stride_ + j]);
                        ASSERT_EQ(ref, lowbd_output_test_[i * stride_ + j])
                            << "loop: " << k << " tx_type: " << tx_type
                            << " tx_size: " << tx_size << " eob: " << eob << " "
                            << j << " x " << i;
                    }
                }
            }
        }
    }

    void run_HandleTransform_match_test() {
        using HandleTxfmFunc = uint64_t (*)(int32_t * output);
        const int num_htf_sizes = 5;
        const TxSize htf_tx_size[num_htf_sizes] = {
            TX_16X64, TX_32X64, TX_64X16, TX_64X32, TX_64X64};
        const HandleTxfmFunc htf_ref_funcs[num_htf_sizes] = {
            HandleTransform16x64_c,
            HandleTransform32x64_c,
            HandleTransform64x16_c,
            HandleTransform64x32_c,
            HandleTransform64x64_c};
        const HandleTxfmFunc htf_asm_funcs[num_htf_sizes] = {
            HandleTransform16x64_avx2,
            HandleTransform32x64_avx2,
            HandleTransform64x16_avx2,
            HandleTransform64x32_avx2,
            HandleTransform64x64_avx2};
        DECLARE_ALIGNED(32, int32_t, input[MAX_TX_SQUARE]);

        for (int idx = 0; idx < num_htf_sizes; ++idx) {
            const TxSize tx_size = htf_tx_size[idx];

            eb_buf_random_s32(input_, MAX_TX_SQUARE);
            memcpy(input, input_, MAX_TX_SQUARE * sizeof(int32_t));

            const uint64_t energy_ref = htf_ref_funcs[idx](input_);
            const uint64_t energy_asm = htf_asm_funcs[idx](input);

            ASSERT_EQ(energy_ref, energy_asm);

            for (int i = 0; i < MAX_TX_SIZE; i++) {
                for (int j = 0; j < MAX_TX_SIZE; j++) {
                    ASSERT_EQ(input_[i * MAX_TX_SIZE + j],
                              input[i * MAX_TX_SIZE + j])
                        << " tx_size: " << tx_size << " " << j << " x " << i;
                }
            }
        }
    }

    void run_HandleTransform_speed_test() {
        using HandleTxfmFunc = uint64_t (*)(int32_t * output);
        const int num_htf_sizes = 5;
        const TxSize htf_tx_size[num_htf_sizes] = {
            TX_16X64, TX_32X64, TX_64X16, TX_64X32, TX_64X64};
        const int widths[num_htf_sizes] = {16, 32, 64, 64, 64};
        const int heights[num_htf_sizes] = {64, 64, 16, 32, 64};
        const HandleTxfmFunc htf_ref_funcs[num_htf_sizes] = {
            HandleTransform16x64_c,
            HandleTransform32x64_c,
            HandleTransform64x16_c,
            HandleTransform64x32_c,
            HandleTransform64x64_c};
        const HandleTxfmFunc htf_asm_funcs[num_htf_sizes] = {
            HandleTransform16x64_avx2,
            HandleTransform32x64_avx2,
            HandleTransform64x16_avx2,
            HandleTransform64x32_avx2,
            HandleTransform64x64_avx2};
        DECLARE_ALIGNED(32, int32_t, input[MAX_TX_SQUARE]);
        double time_c, time_o;
        uint64_t start_time_seconds, start_time_useconds;
        uint64_t middle_time_seconds, middle_time_useconds;
        uint64_t finish_time_seconds, finish_time_useconds;

        for (int idx = 0; idx < num_htf_sizes; ++idx) {
            const TxSize tx_size = htf_tx_size[idx];
            const uint64_t num_loop = 10000000;
            uint64_t energy_ref, energy_asm;

            eb_buf_random_s32(input_, MAX_TX_SQUARE);
            memcpy(input, input_, MAX_TX_SQUARE * sizeof(int32_t));

            EbStartTime(&start_time_seconds, &start_time_useconds);

            for (uint64_t i = 0; i < num_loop; i++)
                energy_ref = htf_ref_funcs[idx](input_);

            EbStartTime(&middle_time_seconds, &middle_time_useconds);

            for (uint64_t i = 0; i < num_loop; i++)
                energy_asm = htf_asm_funcs[idx](input);

            EbStartTime(&finish_time_seconds, &finish_time_useconds);
            EbComputeOverallElapsedTimeMs(start_time_seconds,
                                          start_time_useconds,
                                          middle_time_seconds,
                                          middle_time_useconds,
                                          &time_c);
            EbComputeOverallElapsedTimeMs(middle_time_seconds,
                                          middle_time_useconds,
                                          finish_time_seconds,
                                          finish_time_useconds,
                                          &time_o);

            ASSERT_EQ(energy_ref, energy_asm);

            for (int i = 0; i < MAX_TX_SIZE; i++) {
                for (int j = 0; j < MAX_TX_SIZE; j++) {
                    ASSERT_EQ(input_[i * MAX_TX_SIZE + j],
                              input[i * MAX_TX_SIZE + j])
                        << " tx_size: " << tx_size << " " << j << " x " << i;
                }
            }

            printf("Average Nanoseconds per Function Call\n");
            printf("    HandleTransform%dx%d_c     : %6.2f\n",
                   widths[idx],
                   heights[idx],
                   1000000 * time_c / num_loop);
            printf(
                "    HandleTransform%dx%d_avx2) : %6.2f   (Comparison: "
                "%5.2fx)\n",
                widths[idx],
                heights[idx],
                1000000 * time_o / num_loop,
                time_c / time_o);
        }
    }

  private:
    // clear the coeffs according to eob position, note the coeffs are
    // linear.
    void clear_high_freq_coeffs(const TxSize tx_size, const TxType tx_type,
                                const int eob, const int max_eob) {
        const ScanOrder *scan_order = &av1_scan_orders[tx_size][tx_type];
        const int16_t *scan = scan_order->scan;

        for (int i = eob; i < max_eob; ++i) {
            input_[scan[i]] = 0;
        }
    }

    // fill the pixel_input with random data and do forward transform,
    // Note that the forward transform do not re-pack the coefficients,
    // so we have to re-pack the coefficients after transform for
    // some tx_size;
    void populate_with_random(const int width, const int height,
                              const TxType tx_type, const TxSize tx_size) {
        using FwdTxfm2dFunc = void (*)(int16_t * input,
                                       int32_t * output,
                                       uint32_t stride,
                                       TxType tx_type,
                                       uint8_t bd);

        const FwdTxfm2dFunc fwd_txfm_func[TX_SIZES_ALL] = {
            Av1TransformTwoD_4x4_c,   Av1TransformTwoD_8x8_c,
            Av1TransformTwoD_16x16_c, Av1TransformTwoD_32x32_c,
            Av1TransformTwoD_64x64_c, eb_av1_fwd_txfm2d_4x8_c,
            eb_av1_fwd_txfm2d_8x4_c,     eb_av1_fwd_txfm2d_8x16_c,
            eb_av1_fwd_txfm2d_16x8_c,    eb_av1_fwd_txfm2d_16x32_c,
            eb_av1_fwd_txfm2d_32x16_c,   eb_av1_fwd_txfm2d_32x64_c,
            eb_av1_fwd_txfm2d_64x32_c,   eb_av1_fwd_txfm2d_4x16_c,
            eb_av1_fwd_txfm2d_16x4_c,    eb_av1_fwd_txfm2d_8x32_c,
            eb_av1_fwd_txfm2d_32x8_c,    eb_av1_fwd_txfm2d_16x64_c,
            eb_av1_fwd_txfm2d_64x16_c,
        };

        memset(output_ref_, 0, MAX_TX_SQUARE * sizeof(uint16_t));
        memset(output_test_, 0, MAX_TX_SQUARE * sizeof(uint16_t));
        memset(input_, 0, MAX_TX_SQUARE * sizeof(int32_t));
        memset(pixel_input_, 0, MAX_TX_SQUARE * sizeof(int16_t));
        memset(lowbd_output_test_, 0, MAX_TX_SQUARE * sizeof(uint8_t));
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                pixel_input_[i * stride_ + j] =
                    static_cast<int16_t>(s_bd_rnd_->random());
                output_ref_[i * stride_ + j] = output_test_[i * stride_ + j] =
                    static_cast<uint16_t>(u_bd_rnd_->random());
            }
        }

        fwd_txfm_func[tx_size](
            pixel_input_, input_, stride_, tx_type, static_cast<uint8_t>(bd_));
        // post-process, re-pack the coeffcients
        switch (tx_size) {
        case TX_64X64: HandleTransform64x64_c(input_); break;
        case TX_64X32: HandleTransform64x32_c(input_); break;
        case TX_32X64: HandleTransform32x64_c(input_); break;
        case TX_64X16: HandleTransform64x16_c(input_); break;
        case TX_16X64: HandleTransform16x64_c(input_); break;
        default: break;
        }
        return;
    }

  private:
    SVTRandom *u_bd_rnd_;
    SVTRandom *s_bd_rnd_;

    const int bd_; /**< input param 8bit or 10bit */
    static const int stride_ = MAX_TX_SIZE;
    int16_t *pixel_input_;
    int32_t *input_;
    uint16_t *output_test_;
    uint16_t *output_ref_;
    uint8_t *lowbd_output_test_;
    LowbdInvTxfm2dFunc target_func_;
};

TEST_P(InvTxfm2dAsmTest, sqr_txfm_match_test) {
    for (int i = TX_4X4; i <= TX_64X64; i++) {
        const TxSize tx_size = static_cast<TxSize>(i);
        run_sqr_txfm_match_test(tx_size, 0);
        run_sqr_txfm_match_test(tx_size, 1);
    }
}

TEST_P(InvTxfm2dAsmTest, rect_type1_txfm_match_test) {
    for (int i = TX_4X8; i < TX_SIZES_ALL; i++) {
        const TxSize tx_size = static_cast<TxSize>(i);
        run_rect_type1_txfm_match_test(tx_size);
    }
}

TEST_P(InvTxfm2dAsmTest, rect_type2_txfm_match_test) {
    for (int i = TX_4X8; i < TX_SIZES_ALL; i++) {
        const TxSize tx_size = static_cast<TxSize>(i);
        run_rect_type2_txfm_match_test(tx_size);
    }
}

TEST_P(InvTxfm2dAsmTest, lowbd_txfm_match_test) {
    for (int i = TX_4X4; i < TX_SIZES_ALL; i++) {
        const TxSize tx_size = static_cast<TxSize>(i);
        run_lowbd_txfm_match_test(tx_size);
    }
}

TEST_P(InvTxfm2dAsmTest, HandleTransform_match_test) {
    run_HandleTransform_match_test();
}

TEST_P(InvTxfm2dAsmTest, DISABLED_HandleTransform_speed_test) {
    run_HandleTransform_speed_test();
}

INSTANTIATE_TEST_CASE_P(
    TX, InvTxfm2dAsmTest,
    ::testing::Combine(::testing::Values(eb_av1_lowbd_inv_txfm2d_add_ssse3),
                       ::testing::Values(static_cast<int>(AOM_BITS_8),
                                         static_cast<int>(AOM_BITS_10))));

extern "C" void eb_av1_lowbd_inv_txfm2d_add_avx2(
    const int32_t *input, uint8_t *output_r, int32_t stride_r,
    uint8_t *output_w, int32_t stride_w, TxType tx_type, TxSize tx_size,
    int32_t eob);

INSTANTIATE_TEST_CASE_P(
    TX_AVX2, InvTxfm2dAsmTest,
    ::testing::Combine(::testing::Values(eb_av1_lowbd_inv_txfm2d_add_avx2),
                       ::testing::Values(static_cast<int>(AOM_BITS_8),
                                         static_cast<int>(AOM_BITS_10))));

}  // namespace
