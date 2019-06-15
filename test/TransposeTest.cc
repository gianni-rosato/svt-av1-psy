/*
 * Copyright(c) 2019 Intel Corporation
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "aom_dsp_rtcd.h"
#include "EbDefinitions.h"
#include "EbUnitTestUtility.h"
#include "EbTransforms.h"
#include "random.h"
#include "util.h"
#include "transpose_avx2.h"
#include "transpose_sse2.h"

namespace {

static INLINE void transpose_8bit_16x16_c(const uint8_t *const in,
                                          uint8_t *const out) {
    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < 16; j++) {
            out[16 * j + i] = in[16 * i + j];
        }
    }
}

class TransposeTest : public ::testing::Test {
  public:
    ~TransposeTest();

    void SetUp() {
        rnd_ = new svt_av1_test_tool::SVTRandom(0, (1 << 16) - 1);
    }
    void TearDown() {
        aom_clear_system_state();
    }

  protected:
    void RunCheckOutput();
    void RunSpeedTest();

    uint8_t in_[256];
    uint8_t out_c_[256];
    uint8_t out_o_[256];
    svt_av1_test_tool::SVTRandom *rnd_;
};

TransposeTest::~TransposeTest() {
}

void TransposeTest::RunCheckOutput() {
    __m128i in[16], out[16];

    for (int i = 0; i < 256; ++i) {
        in_[i] = rnd_->Rand8();
    }

    transpose_8bit_16x16_c(in_, out_c_);

    for (int i = 0; i < 16; i++)
        in[i] = _mm_loadu_si128((__m128i *)(in_ + 16 * i));
    transpose_8bit_16x16_reg128bit_avx2(in, out);
    for (int i = 0; i < 16; i++)
        _mm_storeu_si128((__m128i *)(out_o_ + 16 * i), out[i]);
    for (int i = 0; i < 256; ++i)
        ASSERT_EQ(out_c_[i], out_o_[i]) << "[" << i << "]";
}

void TransposeTest::RunSpeedTest() {
    const int num_loops = 100000000;
    __m128i in[16], out_c[16], out_o[16];
    double time_c, time_o;
    uint64_t start_time_seconds, start_time_useconds;
    uint64_t middle_time_seconds, middle_time_useconds;
    uint64_t finish_time_seconds, finish_time_useconds;

    for (int i = 0; i < 256; ++i) {
        in_[i] = rnd_->Rand8();
    }

    for (int i = 0; i < 16; i++)
        in[i] = _mm_loadu_si128((__m128i *)(in_ + 16 * i));

    EbStartTime(&start_time_seconds, &start_time_useconds);

    for (int i = 0; i < num_loops; ++i) {
        transpose_8bit_16x16_sse2(in, out_c);
    }

    EbStartTime(&middle_time_seconds, &middle_time_useconds);

    for (int i = 0; i < num_loops; ++i) {
        transpose_8bit_16x16_reg128bit_avx2(in, out_o);
    }
    EbStartTime(&finish_time_seconds, &finish_time_useconds);

    for (int i = 0; i < 16; i++) {
        _mm_storeu_si128((__m128i *)(out_c_ + 16 * i), out_c[i]);
        _mm_storeu_si128((__m128i *)(out_o_ + 16 * i), out_o[i]);
    }

    for (int i = 0; i < 256; ++i)
        ASSERT_EQ(out_c_[i], out_o_[i]) << "[" << i << "]";

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
    printf("Average Nanoseconds per Function Call\n");
    printf("    transpose_8bit_16x16 (SSE2) : %6.2f\n",
           1000000 * time_c / num_loops);
    printf(
        "    transpose_8bit_16x16 (AVX2) : %6.2f   (Comparison: "
        "%5.2fx)\n",
        1000000 * time_o / num_loops,
        time_c / time_o);
}

TEST_F(TransposeTest, CheckOutput) {
    RunCheckOutput();
}

TEST_F(TransposeTest, DISABLED_Speed) {
    RunSpeedTest();
}

};  // namespace
