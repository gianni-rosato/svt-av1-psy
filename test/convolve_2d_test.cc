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
 * @file Convolve2dTest.cc
 *
 * @brief Unit test for interpolation in inter prediction:
 * - svt_av1_highbd_convolve_2d_copy_sr_avx2
 * - svt_av1_highbd_jnt_convolve_2d_copy_avx2
 * - svt_av1_highbd_convolve_x_sr_avx2
 * - svt_av1_highbd_convolve_y_sr_avx2
 * - svt_av1_highbd_convolve_2d_sr_avx2
 * - svt_av1_highbd_jnt_convolve_x_avx2
 * - svt_av1_highbd_jnt_convolve_y_avx2
 * - svt_av1_highbd_jnt_convolve_2d_avx2
 * - svt_av1_convolve_2d_copy_sr_avx2
 * - svt_av1_jnt_convolve_2d_copy_avx2
 * - svt_av1_convolve_x_sr_avx2
 * - svt_av1_convolve_y_sr_avx2
 * - svt_av1_convolve_2d_sr_avx2
 * - svt_av1_jnt_convolve_x_avx2
 * - svt_av1_jnt_convolve_y_avx2
 * - svt_av1_jnt_convolve_2d_avx2
 *
 * @author Cidana-Wenyao
 *
 ******************************************************************************/
#include <stdlib.h>
#include "gtest/gtest.h"
#include "aom_dsp_rtcd.h"
#include "EbDefinitions.h"
#include "random.h"
#include "util.h"
#include "EbTime.h"
#include "EbUtility.h"
#include "convolve.h"
#include "convolve_avx2.h"
#include "filter.h"
#if defined(_MSC_VER)
#pragma warning(suppress : 4324)
#endif

const int kMaxSize = 128 + 32;  // padding
using svt_av1_test_tool::SVTRandom;
namespace {
using highbd_convolve_func =
    void (*)(const uint16_t *src, int src_stride, uint16_t *dst, int dst_stride,
             int w, int h, const InterpFilterParams *filter_params_x,
             const InterpFilterParams *filter_params_y, const int subpel_x_qn,
             const int subpel_y_qn, ConvolveParams *conv_params, int bd);

using lowbd_convolve_func = void (*)(const uint8_t *src, int src_stride,
                                     uint8_t *dst, int dst_stride, int w, int h,
                                     InterpFilterParams *filter_params_x,
                                     InterpFilterParams *filter_params_y,
                                     const int subpel_x_qn,
                                     const int subpel_y_qn,
                                     ConvolveParams *conv_params);

static const lowbd_convolve_func lowbd_convolve_2d_sr_func_table[] = {
    svt_av1_convolve_2d_sr_avx2,
#ifndef NON_AVX512_SUPPORT
    svt_av1_convolve_2d_sr_avx512
#endif
};

static const lowbd_convolve_func lowbd_convolve_x_sr_func_table[] = {
    svt_av1_convolve_x_sr_avx2,
#ifndef NON_AVX512_SUPPORT
    svt_av1_convolve_x_sr_avx512
#endif
};

static const lowbd_convolve_func lowbd_convolve_y_sr_func_table[] = {
    svt_av1_convolve_y_sr_avx2,
#ifndef NON_AVX512_SUPPORT
    svt_av1_convolve_y_sr_avx512
#endif
};

static const lowbd_convolve_func lowbd_convolve_copy_sr_func_table[] = {
    svt_av1_convolve_2d_copy_sr_avx2,
#ifndef NON_AVX512_SUPPORT
    svt_av1_convolve_2d_copy_sr_avx512
#endif
};

static const lowbd_convolve_func lowbd_jnt_convolve_2d_func_table[] = {
    svt_av1_jnt_convolve_2d_avx2,
#ifndef NON_AVX512_SUPPORT
    svt_av1_jnt_convolve_2d_avx512
#endif
};

static const lowbd_convolve_func lowbd_jnt_convolve_x_func_table[] = {
    svt_av1_jnt_convolve_x_avx2,
#ifndef NON_AVX512_SUPPORT
    svt_av1_jnt_convolve_x_avx512
#endif
};

static const lowbd_convolve_func lowbd_jnt_convolve_y_func_table[] = {
    svt_av1_jnt_convolve_y_avx2,
#ifndef NON_AVX512_SUPPORT
    svt_av1_jnt_convolve_y_avx512
#endif
};

static const lowbd_convolve_func lowbd_jnt_convolve_copy_func_table[] = {
    svt_av1_jnt_convolve_2d_copy_avx2,
#ifndef NON_AVX512_SUPPORT
    svt_av1_jnt_convolve_2d_copy_avx512
#endif
};

/**
 * @brief Unit test for interpolation in inter prediction:
 * - av1_{highbd, }_{jnt, }_convolve_{x, y, 2d}_{sr, copy}_avx2
 *
 * Test strategy:
 * Verify this assembly code by comparing with reference c implementation.
 * Feed the same data and check test output and reference output.
 * Define a template class to handle the common process, and
 * declare sub class to handle different bitdepth and function types.
 *
 * Expect result:
 * Output from assemble functions should be the same with output from c.
 *
 * Test coverage:
 * Test cases:
 * input value: Fill with random values
 * modes: jnt, sr, copy, x, y modes
 * TxSize: all the TxSize.
 * BitDepth: 8bit, 10bit, 12bit
 *
 */
typedef ::testing::tuple<int, int, int, int, BlockSize> Convolve2DParam;
template <typename Sample, typename FuncType>
class AV1Convolve2DTest : public ::testing::TestWithParam<Convolve2DParam> {
  public:
    AV1Convolve2DTest() : bd_(TEST_GET_PARAM(0)) {
        is_jnt_ = 0;
    }

    virtual ~AV1Convolve2DTest() {
    }

    // make the address algined to 32.
    void SetUp() override {
        conv_buf_init_ = reinterpret_cast<ConvBufType *>(
            svt_aom_memalign(32, MAX_SB_SQUARE * sizeof(ConvBufType)));
        conv_buf_ref_ = reinterpret_cast<ConvBufType *>(
            svt_aom_memalign(32, MAX_SB_SQUARE * sizeof(ConvBufType)));
        conv_buf_tst_ = reinterpret_cast<ConvBufType *>(
            svt_aom_memalign(32, MAX_SB_SQUARE * sizeof(ConvBufType)));
        output_init_ = reinterpret_cast<Sample *>(
            svt_aom_memalign(32, MAX_SB_SQUARE * sizeof(Sample)));
        output_ref_ = reinterpret_cast<Sample *>(
            svt_aom_memalign(32, MAX_SB_SQUARE * sizeof(Sample)));
        output_tst_ = reinterpret_cast<Sample *>(
            svt_aom_memalign(32, MAX_SB_SQUARE * sizeof(Sample)));
    }

    void TearDown() override {
        svt_aom_free(conv_buf_init_);
        svt_aom_free(conv_buf_ref_);
        svt_aom_free(conv_buf_tst_);
        svt_aom_free(output_init_);
        svt_aom_free(output_ref_);
        svt_aom_free(output_tst_);
        aom_clear_system_state();
    }

    virtual void run_convolve(int offset_r, int offset_c, int src_stride,
                              int dst_stride, int w, int h,
                              InterpFilterParams *filter_params_x,
                              InterpFilterParams *filter_params_y,
                              const int32_t subpel_x_q4,
                              const int32_t subpel_y_q4,
                              ConvolveParams *conv_params1,
                              ConvolveParams *conv_params2, int bd) {
        (void)offset_r;
        (void)offset_c;
        (void)src_stride;
        (void)dst_stride;
        (void)w;
        (void)h;
        (void)filter_params_x;
        (void)filter_params_y;
        (void)subpel_x_q4;
        (void)subpel_y_q4;
        (void)conv_params1;
        (void)conv_params2;
        (void)bd;
    }

    virtual void speed_convolve(int offset_r, int offset_c, int src_stride,
                                int dst_stride, int w, int h,
                                InterpFilterParams *filter_params_x,
                                InterpFilterParams *filter_params_y,
                                const int32_t subpel_x_q4,
                                const int32_t subpel_y_q4,
                                ConvolveParams *conv_params1,
                                ConvolveParams *conv_params2, int bd) {
        (void)offset_r;
        (void)offset_c;
        (void)src_stride;
        (void)dst_stride;
        (void)w;
        (void)h;
        (void)filter_params_x;
        (void)filter_params_y;
        (void)subpel_x_q4;
        (void)subpel_y_q4;
        (void)conv_params1;
        (void)conv_params2;
        (void)bd;
    }

    void test_convolve(int has_subx, int has_suby, int src_stride,
                       int dst_stride, int output_w, int output_h,
                       InterpFilterParams *filter_params_x,
                       InterpFilterParams *filter_params_y,
                       ConvolveParams *conv_params_ref,
                       ConvolveParams *conv_params_tst, int bd) {
        const int subx_range = has_subx ? 16 : 1;
        const int suby_range = has_suby ? 16 : 1;
        int subx, suby;
        for (subx = 0; subx < subx_range; ++subx) {
            for (suby = 0; suby < suby_range; ++suby) {
                const int offset_r = 3;
                const int offset_c = 3;

                reset_output();

                run_convolve(offset_r,
                             offset_c,
                             src_stride,
                             dst_stride,
                             output_w,
                             output_h,
                             filter_params_x,
                             filter_params_y,
                             subx,
                             suby,
                             conv_params_ref,
                             conv_params_tst,
                             bd);

                if (memcmp(output_ref_,
                           output_tst_,
                           MAX_SB_SQUARE * sizeof(output_ref_[0]))) {
                    for (int i = 0; i < MAX_SB_SIZE; ++i) {
                        for (int j = 0; j < MAX_SB_SIZE; ++j) {
                            int idx = i * MAX_SB_SIZE + j;
                            ASSERT_EQ(output_ref_[idx], output_tst_[idx])
                                << output_w << "x" << output_h
                                << " Pixel mismatch at "
                                   "index "
                                << idx << " = (" << j << ", " << i
                                << "), sub pixel offset = (" << suby << ", "
                                << subx << ")"
                                << " tap = ("
                                << get_convolve_tap(filter_params_x->filter_ptr)
                                << "x"
                                << get_convolve_tap(filter_params_y->filter_ptr)
                                << ") do_average: "
                                << conv_params_tst->do_average
                                << " use_jnt_comp_avg: "
                                << conv_params_tst->use_jnt_comp_avg;
                        }
                    }
                }

                if (is_jnt_) {
                    if (memcmp(conv_buf_ref_,
                               conv_buf_tst_,
                               MAX_SB_SQUARE * sizeof(conv_buf_ref_[0]))) {
                        for (int i = 0; i < MAX_SB_SIZE; ++i) {
                            for (int j = 0; j < MAX_SB_SIZE; ++j) {
                                int idx = i * MAX_SB_SIZE + j;
                                ASSERT_EQ(conv_buf_ref_[idx],
                                          conv_buf_tst_[idx])
                                    << output_w << "x" << output_h
                                    << " Pixel mismatch at index " << idx
                                    << " = (" << j << ", " << i
                                    << "), sub pixel offset = (" << suby << ", "
                                    << subx << ")"
                                    << " tap = ("
                                    << get_convolve_tap(
                                           filter_params_x->filter_ptr)
                                    << "x"
                                    << get_convolve_tap(
                                           filter_params_y->filter_ptr)
                                    << ") do_average: "
                                    << conv_params_tst->do_average
                                    << " use_jnt_comp_avg: "
                                    << conv_params_tst->use_jnt_comp_avg;
                            }
                        }
                    }
                }
            }
        }
    }

    void test_speed(int has_subx, int has_suby, int src_stride, int dst_stride,
                    int output_w, int output_h,
                    InterpFilterParams *filter_params_x,
                    InterpFilterParams *filter_params_y,
                    ConvolveParams *conv_params_ref,
                    ConvolveParams *conv_params_tst, int bd) {
        const int subx_range = has_subx ? 16 : 1;
        const int suby_range = has_suby ? 16 : 1;
        int subx, suby;
        for (subx = 0; subx < subx_range; ++subx) {
            if (subx)
                continue;
            for (suby = 0; suby < suby_range; ++suby) {
                if (suby)
                    continue;
                const int offset_r = 3;
                const int offset_c = 3;

                speed_convolve(offset_r,
                               offset_c,
                               src_stride,
                               dst_stride,
                               output_w,
                               output_h,
                               filter_params_x,
                               filter_params_y,
                               subx,
                               suby,
                               conv_params_ref,
                               conv_params_tst,
                               bd);

                if (memcmp(output_ref_,
                           output_tst_,
                           MAX_SB_SQUARE * sizeof(output_ref_[0]))) {
                    for (int i = 0; i < output_h; ++i) {
                        for (int j = 0; j < output_w; ++j) {
                            int idx = i * MAX_SB_SIZE + j;
                            ASSERT_EQ(output_ref_[idx], output_tst_[idx])
                                << output_w << "x" << output_h
                                << " Pixel mismatch at "
                                   "index "
                                << idx << " = (" << j << ", " << i
                                << "), sub pixel offset = (" << suby << ", "
                                << subx << ") do_average: "
                                << conv_params_tst->do_average
                                << " use_jnt_comp_avg: "
                                << conv_params_tst->use_jnt_comp_avg;
                        }
                    }
                }

                if (is_jnt_) {
                    if (memcmp(conv_buf_ref_,
                               conv_buf_tst_,
                               MAX_SB_SQUARE * sizeof(conv_buf_ref_[0]))) {
                        for (int i = 0; i < output_h; ++i) {
                            for (int j = 0; j < output_w; ++j) {
                                int idx = i * MAX_SB_SIZE + j;
                                ASSERT_EQ(conv_buf_ref_[idx],
                                          conv_buf_tst_[idx])
                                    << output_w << "x" << output_h
                                    << " Pixel mismatch at index " << idx
                                    << " = (" << j << ", " << i
                                    << "), sub pixel offset = (" << suby << ", "
                                    << subx << ") do_average: "
                                    << conv_params_tst->do_average
                                    << " use_jnt_comp_avg: "
                                    << conv_params_tst->use_jnt_comp_avg;
                            }
                        }
                    }
                }
            }
        }
    }

  protected:
    void prepare_data(int w, int h) {
        SVTRandom rnd_(bd_, false);  // bd_-bits, unsigned
        SVTRandom rnd12_(12, false);

        for (int i = 0; i < h; ++i)
            for (int j = 0; j < w; ++j)
                input[i * w + j] = (Sample)rnd_.random();

        for (int i = 0; i < MAX_SB_SQUARE; ++i) {
            conv_buf_init_[i] = rnd12_.random();
            output_init_[i] = (Sample)rnd_.random();
        }
    }

    void reset_output() {
        memcpy(conv_buf_ref_,
               conv_buf_init_,
               MAX_SB_SQUARE * sizeof(*conv_buf_init_));
        memcpy(conv_buf_tst_,
               conv_buf_init_,
               MAX_SB_SQUARE * sizeof(*conv_buf_init_));
        memcpy(
            output_ref_, output_init_, MAX_SB_SQUARE * sizeof(*output_init_));
        memcpy(
            output_tst_, output_init_, MAX_SB_SQUARE * sizeof(*output_init_));
    }

    void run_test() {
        const int input_w = kMaxSize, input_h = kMaxSize;
        const int bd = TEST_GET_PARAM(0);
        const int has_subx = TEST_GET_PARAM(1);
        const int has_suby = TEST_GET_PARAM(2);
        const int block_idx = TEST_GET_PARAM(4);
        int hfilter, vfilter;

        // fill the input data with random
        prepare_data(input_w, input_h);

        setup_common_rtcd_internal(get_cpu_flags_to_use());

        // loop the filter type and subpixel position
        const int output_w = block_size_wide[block_idx];
        const int output_h = block_size_high[block_idx];
        const InterpFilter max_hfilter =
            has_subx ? INTERP_FILTERS_ALL : EIGHTTAP_SMOOTH;
        const InterpFilter max_vfilter =
            has_suby ? INTERP_FILTERS_ALL : EIGHTTAP_SMOOTH;
        for (int compIdx = 0; compIdx < 2; ++compIdx) {
            for (hfilter = EIGHTTAP_REGULAR; hfilter < max_hfilter; ++hfilter) {
                for (vfilter = EIGHTTAP_REGULAR; vfilter < max_vfilter;
                     ++vfilter) {
                    InterpFilterParams filter_params_x =
                        av1_get_interp_filter_params_with_block_size(
                            (InterpFilter)hfilter, output_w >> compIdx);
                    InterpFilterParams filter_params_y =
                        av1_get_interp_filter_params_with_block_size(
                            (InterpFilter)vfilter, output_h >> compIdx);
                    for (int do_average = 0; do_average < 1 + is_jnt_;
                         ++do_average) {
                        // setup convolveParams according to jnt or sr
                        ConvolveParams conv_params_ref;
                        ConvolveParams conv_params_tst;
                        if (is_jnt_) {
                            conv_params_ref =
                                get_conv_params_no_round(0,
                                                         do_average,
                                                         0,
                                                         conv_buf_ref_,
                                                         MAX_SB_SIZE,
                                                         1,
                                                         bd);
                            conv_params_tst =
                                get_conv_params_no_round(0,
                                                         do_average,
                                                         0,
                                                         conv_buf_tst_,
                                                         MAX_SB_SIZE,
                                                         1,
                                                         bd);
                            // Test special case where dist_wtd_comp_avg is not
                            // used
                            conv_params_ref.use_jnt_comp_avg = 0;
                            conv_params_tst.use_jnt_comp_avg = 0;

                        } else {
                            conv_params_ref = get_conv_params_no_round(
                                0, do_average, 0, nullptr, 0, 0, bd);
                            conv_params_tst = get_conv_params_no_round(
                                0, do_average, 0, nullptr, 0, 0, bd);
                        }

                        test_convolve(has_subx,
                                      has_suby,
                                      input_w,
                                      MAX_SB_SIZE,
                                      output_w >> compIdx,
                                      output_h >> compIdx,
                                      &filter_params_x,
                                      &filter_params_y,
                                      &conv_params_ref,
                                      &conv_params_tst,
                                      bd);

                        // AV1 standard won't have 32x4 case.
                        // This only favors some optimization feature which
                        // subsamples 32x8 to 32x4 and triggers 4-tap filter.
                        if ((is_jnt_ == 0) && has_suby &&
                            ((output_w >> compIdx) == 32) &&
                            ((output_h >> compIdx) == 8)) {
                            filter_params_y =
                                av1_get_interp_filter_params_with_block_size(
                                    (InterpFilter)vfilter, 4);
                            test_convolve(has_subx,
                                          has_suby,
                                          input_w,
                                          MAX_SB_SIZE,
                                          32,
                                          4,
                                          &filter_params_x,
                                          &filter_params_y,
                                          &conv_params_ref,
                                          &conv_params_tst,
                                          bd);
                        }

                        if (is_jnt_ == 0)
                            continue;

                        const int quant_dist_lookup_table[2][4][2] = {
                            {{9, 7}, {11, 5}, {12, 4}, {13, 3}},
                            {{7, 9}, {5, 11}, {4, 12}, {3, 13}},
                        };
                        // Test different combination of fwd and bck offset
                        // weights
                        for (int k = 0; k < 2; ++k) {
                            for (int l = 0; l < 4; ++l) {
                                conv_params_ref.use_jnt_comp_avg = 1;
                                conv_params_tst.use_jnt_comp_avg = 1;
                                conv_params_ref.fwd_offset =
                                    quant_dist_lookup_table[k][l][0];
                                conv_params_ref.bck_offset =
                                    quant_dist_lookup_table[k][l][1];
                                conv_params_tst.fwd_offset =
                                    quant_dist_lookup_table[k][l][0];
                                conv_params_tst.bck_offset =
                                    quant_dist_lookup_table[k][l][1];

                                test_convolve(has_subx,
                                              has_suby,
                                              input_w,
                                              MAX_SB_SIZE,
                                              output_w >> compIdx,
                                              output_h >> compIdx,
                                              &filter_params_x,
                                              &filter_params_y,
                                              &conv_params_ref,
                                              &conv_params_tst,
                                              bd);
                            }
                        }
                    }
                }
            }
        }
    }

    void speed_test() {
        const int input_w = kMaxSize, input_h = kMaxSize;
        const int bd = TEST_GET_PARAM(0);
        const int has_subx = TEST_GET_PARAM(1);
        const int has_suby = TEST_GET_PARAM(2);
        const int block_idx = TEST_GET_PARAM(4);
        int hfilter, vfilter;

        // fill the input data with random
        prepare_data(input_w, input_h);
        reset_output();

        // loop the filter type and subpixel position
        const int output_w = block_size_wide[block_idx];
        const int output_h = block_size_high[block_idx];
        const InterpFilter max_hfilter =
            has_subx ? INTERP_FILTERS_ALL : EIGHTTAP_SMOOTH;
        const InterpFilter max_vfilter =
            has_suby ? INTERP_FILTERS_ALL : EIGHTTAP_SMOOTH;
        for (hfilter = EIGHTTAP_REGULAR; hfilter < max_hfilter; ++hfilter) {
            if ((max_hfilter == INTERP_FILTERS_ALL) && (hfilter != BILINEAR))
                continue;
            for (vfilter = EIGHTTAP_REGULAR; vfilter < max_vfilter; ++vfilter) {
                if ((max_vfilter == INTERP_FILTERS_ALL) &&
                    (vfilter != BILINEAR))
                    continue;

                InterpFilterParams filter_params_x =
                    av1_get_interp_filter_params_with_block_size(
                        (InterpFilter)hfilter, output_w);
                InterpFilterParams filter_params_y =
                    av1_get_interp_filter_params_with_block_size(
                        (InterpFilter)vfilter, output_h);
                const int32_t h_tap =
                    get_convolve_tap(filter_params_x.filter_ptr);
                const int32_t v_tap =
                    get_convolve_tap(filter_params_y.filter_ptr);
                if (h_tap != v_tap)
                    continue;
                // setup convolveParams according to jnt or sr
                ConvolveParams conv_params_ref;
                ConvolveParams conv_params_tst;

                conv_params_ref = get_conv_params_no_round(
                    0, is_jnt_, 0, conv_buf_ref_, MAX_SB_SIZE, is_jnt_, bd);
                conv_params_tst = get_conv_params_no_round(
                    0, is_jnt_, 0, conv_buf_tst_, MAX_SB_SIZE, is_jnt_, bd);
                conv_params_ref.use_jnt_comp_avg = 0;
                conv_params_tst.use_jnt_comp_avg = 0;
                conv_params_ref.fwd_offset = conv_params_tst.fwd_offset = 9;
                conv_params_ref.bck_offset = conv_params_tst.bck_offset = 7;

                test_speed(has_subx,
                           has_suby,
                           input_w,
                           MAX_SB_SIZE,
                           output_w,
                           output_h,
                           &filter_params_x,
                           &filter_params_y,
                           &conv_params_ref,
                           &conv_params_tst,
                           bd);
            }
        }
    }

  protected:
    FuncType func_ref_;
    FuncType func_tst_;
    int bd_;
    int is_jnt_;  // jnt or single reference
    Sample input[kMaxSize * kMaxSize];
    ConvBufType *conv_buf_init_;  // aligned address
    ConvBufType *conv_buf_ref_;   // aligned address
    ConvBufType *conv_buf_tst_;   // aligned address
    Sample *output_init_;         // aligned address
    Sample *output_ref_;          // aligned address
    Sample *output_tst_;          // aligned address
};

::testing::internal::ParamGenerator<Convolve2DParam> BuildParams(int has_subx,
                                                                 int has_suby,
                                                                 int func_idx,
                                                                 int highbd) {
    if (highbd)
        return ::testing::Combine(::testing::Range(8, 13, 2),
                                  ::testing::Values(has_subx),
                                  ::testing::Values(has_suby),
                                  ::testing::Values(func_idx),
                                  ::testing::Range(BLOCK_4X4, BlockSizeS_ALL));
    else
        return ::testing::Combine(::testing::Values(8),
                                  ::testing::Values(has_subx),
                                  ::testing::Values(has_suby),
                                  ::testing::Values(func_idx),
                                  ::testing::Range(BLOCK_4X4, BlockSizeS_ALL));
}

class AV1LbdConvolve2DTest
    : public AV1Convolve2DTest<uint8_t, lowbd_convolve_func> {
  public:
    AV1LbdConvolve2DTest() {
    }
    virtual ~AV1LbdConvolve2DTest() {
    }
    void run_convolve(int offset_r, int offset_c, int src_stride,
                      int dst_stride, int output_w, int output_h,
                      InterpFilterParams *filter_params_x,
                      InterpFilterParams *filter_params_y,
                      const int32_t subpel_x_q4, const int32_t subpel_y_q4,
                      ConvolveParams *conv_params_ref,
                      ConvolveParams *conv_params_tst, int bd) override {
        (void)bd;
        func_ref_(input + offset_r * src_stride + offset_c,
                  src_stride,
                  output_ref_,
                  dst_stride,
                  output_w,
                  output_h,
                  filter_params_x,
                  filter_params_y,
                  subpel_x_q4,
                  subpel_y_q4,
                  conv_params_ref);
        func_tst_(input + offset_r * src_stride + offset_c,
                  src_stride,
                  output_tst_,
                  dst_stride,
                  output_w,
                  output_h,
                  filter_params_x,
                  filter_params_y,
                  subpel_x_q4,
                  subpel_y_q4,
                  conv_params_tst);
    }

    void speed_convolve(int offset_r, int offset_c, int src_stride,
                        int dst_stride, int output_w, int output_h,
                        InterpFilterParams *filter_params_x,
                        InterpFilterParams *filter_params_y,
                        const int32_t subpel_x_q4, const int32_t subpel_y_q4,
                        ConvolveParams *conv_params_ref,
                        ConvolveParams *conv_params_tst, int bd) override {
        (void)bd;
        double time_c, time_o;
        uint64_t start_time_seconds, start_time_useconds;
        uint64_t middle_time_seconds, middle_time_useconds;
        uint64_t finish_time_seconds, finish_time_useconds;

        const uint64_t num_loop = 10000000000 / (output_w * output_h);

        svt_av1_get_time(&start_time_seconds, &start_time_useconds);

        for (uint64_t i = 0; i < num_loop; i++) {
            func_ref_(input + offset_r * src_stride + offset_c,
                      src_stride,
                      output_ref_,
                      dst_stride,
                      output_w,
                      output_h,
                      filter_params_x,
                      filter_params_y,
                      subpel_x_q4,
                      subpel_y_q4,
                      conv_params_ref);
        }

        svt_av1_get_time(&middle_time_seconds, &middle_time_useconds);

        for (uint64_t i = 0; i < num_loop; i++) {
            func_tst_(input + offset_r * src_stride + offset_c,
                      src_stride,
                      output_tst_,
                      dst_stride,
                      output_w,
                      output_h,
                      filter_params_x,
                      filter_params_y,
                      subpel_x_q4,
                      subpel_y_q4,
                      conv_params_tst);
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

        printf("convolve(%3dx%3d, tap (%d, %d)): %6.2f\n",
               output_w,
               output_h,
               get_convolve_tap(filter_params_x->filter_ptr),
               get_convolve_tap(filter_params_y->filter_ptr),
               time_c / time_o);
    }
};

class AV1LbdJntConvolve2DTest : public AV1LbdConvolve2DTest {
  public:
    AV1LbdJntConvolve2DTest() {
        is_jnt_ = 1;
        func_ref_ = svt_av1_jnt_convolve_2d_c;
        const int has_subx = TEST_GET_PARAM(1);
        const int has_suby = TEST_GET_PARAM(2);
        const int func_idx = TEST_GET_PARAM(3);
        if (has_subx == 1 && has_suby == 1)
            func_tst_ = lowbd_jnt_convolve_2d_func_table[func_idx];
        else if (has_subx == 1)
            func_tst_ = lowbd_jnt_convolve_x_func_table[func_idx];
        else if (has_suby == 1)
            func_tst_ = lowbd_jnt_convolve_y_func_table[func_idx];
        else
            func_tst_ = lowbd_jnt_convolve_copy_func_table[func_idx];
        bd_ = TEST_GET_PARAM(0);
    }
    virtual ~AV1LbdJntConvolve2DTest() {
    }
};

TEST_P(AV1LbdJntConvolve2DTest, MatchTest) {
    run_test();
}

TEST_P(AV1LbdJntConvolve2DTest, DISABLED_SpeedTest) {
    speed_test();
}

INSTANTIATE_TEST_CASE_P(ConvolveTestCOPY, AV1LbdJntConvolve2DTest,
                        BuildParams(0, 0, 0, 0));
INSTANTIATE_TEST_CASE_P(ConvolveTestX, AV1LbdJntConvolve2DTest,
                        BuildParams(1, 0, 0, 0));
INSTANTIATE_TEST_CASE_P(ConvolveTestY, AV1LbdJntConvolve2DTest,
                        BuildParams(0, 1, 0, 0));
INSTANTIATE_TEST_CASE_P(ConvolveTest2D, AV1LbdJntConvolve2DTest,
                        BuildParams(1, 1, 0, 0));

#ifndef NON_AVX512_SUPPORT
INSTANTIATE_TEST_CASE_P(ConvolveTestCOPY_AVX512, AV1LbdJntConvolve2DTest,
                        BuildParams(0, 0, 1, 0));
INSTANTIATE_TEST_CASE_P(ConvolveTestX_AVX512, AV1LbdJntConvolve2DTest,
                        BuildParams(1, 0, 1, 0));
INSTANTIATE_TEST_CASE_P(ConvolveTestY_AVX512, AV1LbdJntConvolve2DTest,
                        BuildParams(0, 1, 1, 0));
INSTANTIATE_TEST_CASE_P(ConvolveTest2D_AVX512, AV1LbdJntConvolve2DTest,
                        BuildParams(1, 1, 1, 0));
#endif

class AV1LbdSrConvolve2DTest : public AV1LbdConvolve2DTest {
  public:
    AV1LbdSrConvolve2DTest() {
        is_jnt_ = 0;
        func_ref_ = svt_av1_convolve_2d_sr_c;
        const int has_subx = TEST_GET_PARAM(1);
        const int has_suby = TEST_GET_PARAM(2);
        const int func_idx = TEST_GET_PARAM(3);
        if (has_subx == 1 && has_suby == 1)
            func_tst_ = lowbd_convolve_2d_sr_func_table[func_idx];
        else if (has_subx == 1)
            func_tst_ = lowbd_convolve_x_sr_func_table[func_idx];
        else if (has_suby == 1)
            func_tst_ = lowbd_convolve_y_sr_func_table[func_idx];
        else
            func_tst_ = lowbd_convolve_copy_sr_func_table[func_idx];
        bd_ = TEST_GET_PARAM(0);
    }
    virtual ~AV1LbdSrConvolve2DTest() {
    }
};

TEST_P(AV1LbdSrConvolve2DTest, MatchTest) {
    run_test();
}

TEST_P(AV1LbdSrConvolve2DTest, DISABLED_SpeedTest) {
    speed_test();
}

INSTANTIATE_TEST_CASE_P(ConvolveTestCopy, AV1LbdSrConvolve2DTest,
                        BuildParams(0, 0, 0, 0));
INSTANTIATE_TEST_CASE_P(ConvolveTestX, AV1LbdSrConvolve2DTest,
                        BuildParams(1, 0, 0, 0));
INSTANTIATE_TEST_CASE_P(ConvolveTestY, AV1LbdSrConvolve2DTest,
                        BuildParams(0, 1, 0, 0));
INSTANTIATE_TEST_CASE_P(ConvolveTest2D, AV1LbdSrConvolve2DTest,
                        BuildParams(1, 1, 0, 0));

#ifndef NON_AVX512_SUPPORT
INSTANTIATE_TEST_CASE_P(ConvolveTestCopy_AVX512, AV1LbdSrConvolve2DTest,
                        BuildParams(0, 0, 1, 0));
INSTANTIATE_TEST_CASE_P(ConvolveTestX_AVX512, AV1LbdSrConvolve2DTest,
                        BuildParams(1, 0, 1, 0));
INSTANTIATE_TEST_CASE_P(ConvolveTestY_AVX512, AV1LbdSrConvolve2DTest,
                        BuildParams(0, 1, 1, 0));
INSTANTIATE_TEST_CASE_P(ConvolveTest2D_AVX512, AV1LbdSrConvolve2DTest,
                        BuildParams(1, 1, 1, 0));
#endif

class AV1HbdConvolve2DTest
    : public AV1Convolve2DTest<uint16_t, highbd_convolve_func> {
  public:
    AV1HbdConvolve2DTest() {
    }
    virtual ~AV1HbdConvolve2DTest() {
    }
    void run_convolve(int offset_r, int offset_c, int src_stride,
                      int dst_stride, int blk_w, int blk_h,
                      InterpFilterParams *filter_params_x,
                      InterpFilterParams *filter_params_y,
                      const int32_t subpel_x_q4, const int32_t subpel_y_q4,
                      ConvolveParams *conv_params_ref,
                      ConvolveParams *conv_params_tst, int bd) override {
        func_ref_(input + offset_r * src_stride + offset_c,
                  src_stride,
                  output_ref_,
                  dst_stride,
                  blk_w,
                  blk_h,
                  filter_params_x,
                  filter_params_y,
                  subpel_x_q4,
                  subpel_y_q4,
                  conv_params_ref,
                  bd);
        func_tst_(input + offset_r * src_stride + offset_c,
                  src_stride,
                  output_tst_,
                  dst_stride,
                  blk_w,
                  blk_h,
                  filter_params_x,
                  filter_params_y,
                  subpel_x_q4,
                  subpel_y_q4,
                  conv_params_tst,
                  bd);
    }

    void speed_convolve(int offset_r, int offset_c, int src_stride,
                        int dst_stride, int output_w, int output_h,
                        InterpFilterParams *filter_params_x,
                        InterpFilterParams *filter_params_y,
                        const int32_t subpel_x_q4, const int32_t subpel_y_q4,
                        ConvolveParams *conv_params_ref,
                        ConvolveParams *conv_params_tst, int bd) override {
        double time_c, time_o;
        uint64_t start_time_seconds, start_time_useconds;
        uint64_t middle_time_seconds, middle_time_useconds;
        uint64_t finish_time_seconds, finish_time_useconds;

        const uint64_t num_loop = 10000000 / (output_w * output_h);

        svt_av1_get_time(&start_time_seconds, &start_time_useconds);

        for (uint64_t i = 0; i < num_loop; i++) {
            func_ref_(input + offset_r * src_stride + offset_c,
                      src_stride,
                      output_ref_,
                      dst_stride,
                      output_w,
                      output_h,
                      filter_params_x,
                      filter_params_y,
                      subpel_x_q4,
                      subpel_y_q4,
                      conv_params_ref,
                      bd);
        }

        svt_av1_get_time(&middle_time_seconds, &middle_time_useconds);

        for (uint64_t i = 0; i < num_loop; i++) {
            func_tst_(input + offset_r * src_stride + offset_c,
                      src_stride,
                      output_tst_,
                      dst_stride,
                      output_w,
                      output_h,
                      filter_params_x,
                      filter_params_y,
                      subpel_x_q4,
                      subpel_y_q4,
                      conv_params_tst,
                      bd);
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

        printf(
            "convolve(%3dx%3d): %6.2f\n", output_w, output_h, time_c / time_o);
    }
};

class AV1HbdJntConvolve2DTest : public AV1HbdConvolve2DTest {
  public:
    AV1HbdJntConvolve2DTest() {
        is_jnt_ = 1;
        func_ref_ = svt_av1_highbd_jnt_convolve_2d_c;
        const int has_subx = TEST_GET_PARAM(1);
        const int has_suby = TEST_GET_PARAM(2);
        if (has_subx == 1 && has_suby == 1)
            func_tst_ = svt_av1_highbd_jnt_convolve_2d_avx2;
        else if (has_subx == 1)
            func_tst_ = svt_av1_highbd_jnt_convolve_x_avx2;
        else if (has_suby == 1)
            func_tst_ = svt_av1_highbd_jnt_convolve_y_avx2;
        else
            func_tst_ = svt_av1_highbd_jnt_convolve_2d_copy_avx2;

        bd_ = TEST_GET_PARAM(0);
    }
    virtual ~AV1HbdJntConvolve2DTest() {
    }
};

TEST_P(AV1HbdJntConvolve2DTest, MatchTest) {
    run_test();
}

TEST_P(AV1HbdJntConvolve2DTest, DISABLED_SpeedTest) {
    speed_test();
}

INSTANTIATE_TEST_CASE_P(AVX2_COPY, AV1HbdJntConvolve2DTest,
                        BuildParams(0, 0, 0, 1));
INSTANTIATE_TEST_CASE_P(ConvolveTest2D, AV1HbdJntConvolve2DTest,
                        BuildParams(1, 1, 0, 1));
INSTANTIATE_TEST_CASE_P(ConvolveTestX, AV1HbdJntConvolve2DTest,
                        BuildParams(1, 0, 0, 1));
INSTANTIATE_TEST_CASE_P(ConvolveTestY, AV1HbdJntConvolve2DTest,
                        BuildParams(0, 1, 0, 1));

class AV1HbdSrConvolve2DTest : public AV1HbdConvolve2DTest {
  public:
    AV1HbdSrConvolve2DTest() {
        is_jnt_ = 0;
        func_ref_ = svt_av1_highbd_convolve_2d_sr_c;
        const int has_subx = TEST_GET_PARAM(1);
        const int has_suby = TEST_GET_PARAM(2);
        if (has_subx == 1 && has_suby == 1)
            func_tst_ = svt_av1_highbd_convolve_2d_sr_avx2;
        else if (has_subx == 1)
            func_tst_ = svt_av1_highbd_convolve_x_sr_avx2;
        else if (has_suby == 1)
            func_tst_ = svt_av1_highbd_convolve_y_sr_avx2;
        else
            func_tst_ = svt_av1_highbd_convolve_2d_copy_sr_avx2;
        bd_ = TEST_GET_PARAM(0);
    }
    virtual ~AV1HbdSrConvolve2DTest() {
    }
};

TEST_P(AV1HbdSrConvolve2DTest, MatchTest) {
    run_test();
}

TEST_P(AV1HbdSrConvolve2DTest, DISABLED_SpeedTest) {
    speed_test();
}

INSTANTIATE_TEST_CASE_P(ConvolveTestX, AV1HbdSrConvolve2DTest,
                        BuildParams(1, 0, 0, 1));
INSTANTIATE_TEST_CASE_P(ConvolveTest2D, AV1HbdSrConvolve2DTest,
                        BuildParams(1, 1, 0, 1));
INSTANTIATE_TEST_CASE_P(ConvolveTestY, AV1HbdSrConvolve2DTest,
                        BuildParams(0, 1, 0, 1));
INSTANTIATE_TEST_CASE_P(ConvolveTestCopy, AV1HbdSrConvolve2DTest,
                        BuildParams(0, 0, 0, 1));
}  // namespace
