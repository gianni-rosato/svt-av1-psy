/*
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
#include "gtest/gtest.h"
#include "aom_dsp_rtcd.h"
#include "EbDefinitions.h"
#include "random.h"
#include "util.h"
#include "EbUnitTestUtility.h"
#include "acm_random.h"
#include "corner_match.h"

using libaom_test::ACMRandom;

namespace {

typedef double (*ComputeCrossCorrFunc)(unsigned char *im1, int stride1, int x1,
                                       int y1, unsigned char *im2, int stride2,
                                       int x2, int y2);

using ::testing::make_tuple;
using ::testing::tuple;
typedef tuple<int, ComputeCrossCorrFunc> CornerMatchParam;

class AV1CornerMatchTest : public ::testing::TestWithParam<CornerMatchParam> {
 public:
  virtual ~AV1CornerMatchTest();
  virtual void SetUp();

  virtual void TearDown();

 protected:
  void RunCheckOutput(int run_times);
  ComputeCrossCorrFunc target_func;

  libaom_test::ACMRandom rnd_;
};

AV1CornerMatchTest::~AV1CornerMatchTest() {}
void AV1CornerMatchTest::SetUp() {
  rnd_.Reset(ACMRandom::DeterministicSeed());
    target_func = TEST_GET_PARAM(1);
}
void AV1CornerMatchTest::TearDown() { }

void AV1CornerMatchTest::RunCheckOutput(int run_times) {
  const int w = 128, h = 128;
  const int num_iters = 10000;
  int i, j;
  uint64_t start_time_seconds, start_time_useconds;
  uint64_t middle_time_seconds, middle_time_useconds;
  uint64_t finish_time_seconds, finish_time_useconds;
  double time_c = 0, time_o = 0, time;

  uint8_t *input1 = new uint8_t[w * h];
  uint8_t *input2 = new uint8_t[w * h];

  // Test the two extreme cases:
  // i) Random data, should have correlation close to 0
  // ii) Linearly related data + noise, should have correlation close to 1
  int mode = TEST_GET_PARAM(0);
  if (mode == 0) {
    for (i = 0; i < h; ++i)
      for (j = 0; j < w; ++j) {
        input1[i * w + j] = rnd_.Rand8();
        input2[i * w + j] = rnd_.Rand8();
      }
  } else if (mode == 1) {
    for (i = 0; i < h; ++i)
      for (j = 0; j < w; ++j) {
        int v = rnd_.Rand8();
        input1[i * w + j] = v;
        input2[i * w + j] = (v / 2) + (rnd_.Rand8() & 15);
      }
  }

  for (i = 0; i < num_iters; ++i) {
    int x1 = MATCH_SZ_BY2 + rnd_.PseudoUniform(w - 2 * MATCH_SZ_BY2);
    int y1 = MATCH_SZ_BY2 + rnd_.PseudoUniform(h - 2 * MATCH_SZ_BY2);
    int x2 = MATCH_SZ_BY2 + rnd_.PseudoUniform(w - 2 * MATCH_SZ_BY2);
    int y2 = MATCH_SZ_BY2 + rnd_.PseudoUniform(h - 2 * MATCH_SZ_BY2);

    double res_c =
        eb_av1_compute_cross_correlation_c(input1, w, x1, y1, input2, w, x2, y2);
    double res_simd = target_func(input1, w, x1, y1, input2, w, x2, y2);

    if (run_times > 1) {
      svt_av1_get_time(&start_time_seconds, &start_time_useconds);
      for (j = 0; j < run_times; j++) {
          eb_av1_compute_cross_correlation_c(
              input1, w, x1, y1, input2, w, x2, y2);
      }
      svt_av1_get_time(&middle_time_seconds, &middle_time_useconds);

      for (j = 0; j < run_times; j++) {
        target_func(input1, w, x1, y1, input2, w, x2, y2);
      }

      svt_av1_get_time(&finish_time_seconds, &finish_time_useconds);

      time = svt_av1_compute_overall_elapsed_time_ms(start_time_seconds,
                                                     start_time_useconds,
                                                     middle_time_seconds,
                                                     middle_time_useconds);
      time_c += time;
      time = svt_av1_compute_overall_elapsed_time_ms(middle_time_seconds,
                                                     middle_time_useconds,
                                                     finish_time_seconds,
                                                     finish_time_useconds);
      time_o += time;



    } else {
      ASSERT_EQ(res_simd, res_c);
    }
  }

   if (run_times > 1) {
      printf("Average Nanoseconds per Function Call\n");
      printf("    eb_av1_compute_cross_correlation_c : %6.2f\n",
             1000000 * time_c / run_times * num_iters);
      printf(
          "    av1_compute_cross_correlation (AVX2) : %6.2f   (Comparison: "
          "%5.2fx)\n",
          1000000 * time_o / run_times * num_iters,
          time_c / time_o);

   }

  delete[] input1;
  delete[] input2;
}

TEST_P(AV1CornerMatchTest, CheckOutput) { RunCheckOutput(1); }
TEST_P(AV1CornerMatchTest, DISABLED_Speed) { RunCheckOutput(1000); }

INSTANTIATE_TEST_CASE_P(
    AV1CornerMatchTest, AV1CornerMatchTest,
    ::testing::Values(make_tuple(0, &eb_av1_compute_cross_correlation_avx2),
                      make_tuple(1, &eb_av1_compute_cross_correlation_avx2)));

}
