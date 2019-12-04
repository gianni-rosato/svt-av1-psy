/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file AVCStyleMcpTest.cc
 *
 * @brief Unit test for Mcp functions:
 * - avc_style_copy_sse2
 * - avc_style_luma_interpolation_filter_pose_ssse3
 * - avc_style_luma_interpolation_filter_posf_ssse3
 * - avc_style_luma_interpolation_filter_posg_ssse3
 * - avc_style_luma_interpolation_filter_posi_ssse3
 * - avc_style_luma_interpolation_filter_posj_ssse3
 * - avc_style_luma_interpolation_filter_posk_ssse3
 * - avc_style_luma_interpolation_filter_posp_ssse3
 * - avc_style_luma_interpolation_filter_posq_ssse3
 * - avc_style_luma_interpolation_filter_posr_ssse3
 * - avc_style_luma_interpolation_filter_vertical_ssse3_intrin
 * - avc_style_luma_interpolation_filter_horizontal_ssse3_intrin
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
#include "EbAvcStyleMcp.h"
#include "EbMcp_SSE2.h"
#include "EbAvcStyleMcp_SSSE3.h"
#include "EbAvcStyleMcp_SSE2.h"
#include "EbDefinitions.h"
#include "EbIntraPrediction.h"
#include "random.h"
#include "util.h"

using svt_av1_test_tool::SVTRandom;  // to generate the random

namespace {
#define MAX_SEARCH_REGION_SIZE \
    ((MAX_SEARCH_AREA_WIDTH_CH) * (MAX_SEARCH_AREA_HEIGHT_CH))
#define MAX_PIC_SIZE \
    (MAX_PICTURE_WIDTH_SIZE + 64 + 4 + 64 + 4) * MAX_PICTURE_HEIGHT_SIZE

typedef enum { REF_MIN, REF_MAX, REF_RANDOM } TestPattern;
TestPattern TEST_PATTERNS[] = {REF_MIN, REF_MAX, REF_RANDOM};

typedef std::tuple<uint32_t, uint32_t> PUSize;
PUSize TEST_PU_SIZES[] = {
    PUSize(64, 64), PUSize(64, 32), PUSize(32, 64), PUSize(32, 32),
    PUSize(32, 16), PUSize(16, 32), PUSize(16, 16), PUSize(16, 8),
    PUSize(8, 16),  PUSize(8, 8),   PUSize(8, 4),   PUSize(4, 4),
    PUSize(4, 8),   PUSize(4, 16),  PUSize(16, 4),  PUSize(8, 32),
    PUSize(32, 8),  PUSize(16, 64), PUSize(64, 16), PUSize(24, 24),
    PUSize(24, 16), PUSize(16, 24), PUSize(24, 8),  PUSize(8, 24),
    PUSize(64, 24), PUSize(48, 24), PUSize(32, 24), PUSize(24, 32),
    PUSize(48, 48), PUSize(48, 16), PUSize(48, 32), PUSize(16, 48),
    PUSize(32, 48), PUSize(48, 64), PUSize(64, 48)};

typedef std::tuple<uint32_t, uint32_t> SearchArea;
SearchArea TEST_AREAS[] = {SearchArea(64, 64),
                           SearchArea(64, 32),
                           SearchArea(48, 16),
                           SearchArea(16, 10),
                           SearchArea(112, 112),
                           SearchArea(128, 128),
                           SearchArea(1280, 1280),
                           SearchArea(640, 640),
                           SearchArea(288, 248),
                           SearchArea(208, 168),
                           SearchArea(168, 128),
                           SearchArea(128, 80),
                           SearchArea(64, 48),
                           SearchArea(80, 80)};

PUSize TEST_PU_HELPER[] = {PUSize(64, 64)};

typedef void (*MCP_REF_FUNC)(EbByte refPic, uint32_t srcStride, EbByte dst,
                             uint32_t dstStride, uint32_t puWidth,
                             uint32_t puHeight, EbByte tempBuf,
                             uint32_t fracPos);
typedef void (*MCP_TEST_FUNC)(EbByte ref_pic, uint32_t src_stride, EbByte dst,
                              uint32_t dst_stride, uint32_t pu_width,
                              uint32_t pu_height, EbByte temp_buf, EbBool skip,
                              uint32_t frac_pos);

typedef struct {
    const char *name;
    MCP_REF_FUNC ref_func;
    MCP_TEST_FUNC test_func;
} AVCStyleMcpFuncPair;
static const AVCStyleMcpFuncPair AVC_style_c_sse3_func_pairs[] = {
    {"posA", avc_style_copy_c, avc_style_copy_sse2},
    {"pose",
     avc_style_luma_interpolation_filter_pose_c,
     avc_style_luma_interpolation_filter_pose_ssse3},
    {"posf",
     avc_style_luma_interpolation_filter_posf_c,
     avc_style_luma_interpolation_filter_posf_ssse3},
    {"posg",
     avc_style_luma_interpolation_filter_posg_c,
     avc_style_luma_interpolation_filter_posg_ssse3},
    {"posi",
     avc_style_luma_interpolation_filter_posi_c,
     avc_style_luma_interpolation_filter_posi_ssse3},
    {"posj",
     avc_style_luma_interpolation_filter_posj_c,
     avc_style_luma_interpolation_filter_posj_ssse3},
    {"posk",
     avc_style_luma_interpolation_filter_posk_c,
     avc_style_luma_interpolation_filter_posk_ssse3},
    {"posp",
     avc_style_luma_interpolation_filter_posp_c,
     avc_style_luma_interpolation_filter_posp_ssse3},
    {"posq",
     avc_style_luma_interpolation_filter_posq_c,
     avc_style_luma_interpolation_filter_posq_ssse3},
    {"posr",
     avc_style_luma_interpolation_filter_posr_c,
     avc_style_luma_interpolation_filter_posr_ssse3},
    {"vertical",
     avc_style_luma_interpolation_filter_vertical_c,
     avc_style_luma_interpolation_filter_vertical_ssse3_intrin},
    {"horizontal",
     avc_style_luma_interpolation_filter_horizontal_c,
     avc_style_luma_interpolation_filter_horizontal_ssse3_intrin}
};

const int NUM_FUNCS = sizeof(AVC_style_c_sse3_func_pairs)
                    / sizeof(AVC_style_c_sse3_func_pairs[0]);

typedef std::tuple<PUSize, TestPattern> TestPUParam;
typedef std::tuple<SearchArea, TestPattern> TestSearchRegionParam;

/**
 * @brief Unit test for Mcp functions include:
 * aseembly functions in AVC_style_c_sse3_func_pairs
 *
 *
 * Test strategy:
 *  This test cases combine different pu size or search area with different
 *  test pattern(REF_MIN, REF_MAX, RANDOM). Check the result by compare
 *  with c implementation.
 *
 * Expect result:
 *  Results from assembly and c functions are equal.
 *
 * Test coverage:
 *  All functions inside avc_style_luma_interpolation_filter_helper_ssse3
 *
 * Test cases:
 *  PU Size: Width {4, 8, 16, 24, 32, 48, 64} x height{ 4, 8, 16, 24, 32, 48,
 * 64) Search Area: Width:{ 16, 48, 64, 80, 112, 128, 168, 208, 288, 640, 1280}
 * x height{ 9, 16, 32, 48, 64, 80, 112, 128, 168, 248, 640, 1280) Test
 * vector pattern: {REF_MIN, REF_MAX, RANDOM}
 *
 */
class AVCStyleMcpTestBase : public ::testing::Test {
  public:
    AVCStyleMcpTestBase(const uint32_t block_width, const uint32_t block_height,
                        TestPattern test_pattern, EbBool is_PU_block) {
        test_pattern_ = test_pattern;
        src_stride_ = MAX_PICTURE_WIDTH_SIZE + 64 + 4 + 64 + 4;
        is_PU_block_ = is_PU_block;
        if (is_PU_block) {
            // interpolate PU blocks
            block_width_ = block_width;
            block_height_ = block_height;
            block_size_ = MAX_PU_SIZE * MAX_PU_SIZE;
            dst_stride_ = MAX_PU_SIZE;
            src_offset_ = 2 + 2 * (MAX_PICTURE_WIDTH_SIZE + 64 + 4 + 64 + 4);
        } else {  // interpolate search areas
            // The width must be oversized by 2 to account for edge conditions,
            // and should be a multiple of 8 to align with the ASM kernel.
            block_width_ = ROUND_UP_MUL_8(block_width + BLOCK_SIZE_64 + 2);
            block_height_ = block_height + BLOCK_SIZE_64 + 4;
            block_size_ = block_width_ * block_height_;
            dst_stride_ = block_width_;
            src_offset_ = 2 + 2 * src_stride_;
        }
        // dst size should be the same with block size.
        dst_size_ = block_size_;
        // The interpolation will required more than triple times of block size,
        // so we allocate four times of block sizes.
        tmp_buf_size_ = 4 * block_size_;

        src_ = (uint8_t *)eb_aom_memalign(8, MAX_PIC_SIZE * sizeof(*src_));
        tmp_buf_ =
            (uint8_t *)eb_aom_memalign(8, tmp_buf_size_ * sizeof(*tmp_buf_));
        dst1_ = (uint8_t *)eb_aom_memalign(8, dst_size_ * sizeof(*dst1_));
        dst2_ = (uint8_t *)eb_aom_memalign(8, dst_size_ * sizeof(*dst2_));
        memset(tmp_buf_, 0, tmp_buf_size_ * sizeof(*tmp_buf_));
        memset(dst1_, 0, dst_size_ * sizeof(*dst1_));
        memset(dst2_, 0, dst_size_ * sizeof(*dst2_));
    }

    void TearDown() override {
        if (src_)
            eb_aom_free(src_);
        if (tmp_buf_)
            eb_aom_free(tmp_buf_);
        if (dst1_)
            eb_aom_free(dst1_);
        if (dst2_)
            eb_aom_free(dst2_);
    }

    virtual void run_Mcp_test() {
        prepare_data();

        // Skip the avc_style_copy_sse2 for search area test cases;
        // Since the function requires that the width should be
        // 4, 8, 12, 16, 24, 32, 48, 64 or 128, and it can not
        // be meet in the search area test cases.
        uint32_t i = is_PU_block_ ? 0 : 1;
        for (; i < NUM_FUNCS; i++) {
            AVC_style_c_sse3_func_pairs[i].ref_func(src_ + src_offset_,
                                                    src_stride_,
                                                    dst1_,
                                                    dst_stride_,
                                                    block_width_,
                                                    block_height_,
                                                    tmp_buf_,
                                                    2);

            AVC_style_c_sse3_func_pairs[i].test_func(src_ + src_offset_,
                                                     src_stride_,
                                                     dst2_,
                                                     dst_stride_,
                                                     block_width_,
                                                     block_height_,
                                                     tmp_buf_,
                                                     EB_FALSE,
                                                     2);

            int fail_pixel_Count = 0;
            for (uint32_t j = 0; j < block_height_; j++) {
                for (uint32_t k = 0; k < block_width_; k++) {
                    if (dst1_[k + j * dst_stride_] !=
                        dst2_[k + j * dst_stride_])
                        fail_pixel_Count++;
                }
            }
            EXPECT_EQ(0, fail_pixel_Count)
                << "Mcp test failed, "
                << " failed count: " << fail_pixel_Count << " at func: ["
                << AVC_style_c_sse3_func_pairs[i].name << "] ";
        }
    }

    virtual void run_Mcp_test_helper() {
        prepare_data();

        // Skip the avc_style_copy_sse2 for search area test cases;
        // Since the function requires that the width should be
        // 4, 8, 12, 16, 24, 32, 48, 64 or 128, and it can not
        // be meet in the search area test cases.
        uint32_t i = is_PU_block_ ? 0 : 1;
        for (; i < 16; i++) {
            avc_style_luma_interpolation_filter_helper_c(src_ + src_offset_,
                                                src_stride_,
                                                dst1_,
                                                dst_stride_,
                                                block_width_,
                                                block_height_,
                                                tmp_buf_,
                                                EB_FALSE,
                                                2,
                                                i);

            avc_style_luma_interpolation_filter_helper_ssse3(
                                                src_ + src_offset_,
                                                src_stride_,
                                                dst2_,
                                                dst_stride_,
                                                block_width_,
                                                block_height_,
                                                tmp_buf_,
                                                EB_FALSE,
                                                2,
                                                i);

            int fail_pixel_Count = 0;
            for (uint32_t j = 0; j < block_height_; j++) {
                for (uint32_t k = 0; k < block_width_; k++) {
                    if (dst1_[k + j * dst_stride_] !=
                        dst2_[k + j * dst_stride_])
                        fail_pixel_Count++;
                }
            }
            EXPECT_EQ(0, fail_pixel_Count)
                << "Mcp test failed, "
                << " failed count: " << fail_pixel_Count << " at func: ["
                << AVC_style_c_sse3_func_pairs[i].name << "] ";
        }
    }

  protected:
    virtual void prepare_data() {
        const int32_t mask = (1 << 8) - 1;
        SVTRandom rnd(0, mask);
        switch (test_pattern_) {
        case REF_MIN: {
            memset(src_, 0, MAX_PIC_SIZE * sizeof(*src_));
            break;
        }
        case REF_MAX: {
            memset(src_, mask, MAX_PIC_SIZE * sizeof(*src_));
            break;
        }
        case REF_RANDOM: {
            for (uint32_t i = 0; i < MAX_PIC_SIZE; i++)
                src_[i] = rnd.random();
            break;
        }
        default: break;
        }
    }

    TestPattern test_pattern_;

    // block dimensions, the original block could be a PU or search area;
    uint32_t block_width_, block_height_;
    // buffer sizes
    uint32_t block_size_;
    uint32_t tmp_buf_size_, dst_size_;
    // interpolation requires neighbor pixels
    uint32_t src_offset_;
    uint32_t src_stride_;
    uint32_t dst_stride_;
    // buffer to store the original full pixels
    uint8_t *src_;
    // buffer to store the intermedia interpolated result.
    uint8_t *tmp_buf_;
    // buffer to store the final interpolated values;
    uint8_t *dst1_, *dst2_;
    EbBool is_PU_block_;
};

class AVCStyleMcpPUTest : public AVCStyleMcpTestBase,
                          public ::testing::WithParamInterface<TestPUParam> {
  public:
    AVCStyleMcpPUTest()
        : AVCStyleMcpTestBase(std::get<0>(TEST_GET_PARAM(0)),
                              std::get<1>(TEST_GET_PARAM(0)), TEST_GET_PARAM(1),
                              true) {
    }
};

TEST_P(AVCStyleMcpPUTest, AVCStyleMcpPUTest) {
    run_Mcp_test();
};

INSTANTIATE_TEST_CASE_P(AVCMCPPU, AVCStyleMcpPUTest,
                        ::testing::Combine(::testing::ValuesIn(TEST_PU_SIZES),
                                           ::testing::ValuesIn(TEST_PATTERNS)));

class AVCStyleMcpSearchRegionTest
    : public AVCStyleMcpTestBase,
      public ::testing::WithParamInterface<TestSearchRegionParam> {
  public:
    AVCStyleMcpSearchRegionTest()
        : AVCStyleMcpTestBase(std::get<0>(TEST_GET_PARAM(0)),
                              std::get<1>(TEST_GET_PARAM(0)), TEST_GET_PARAM(1),
                              false) {
    }
};

TEST_P(AVCStyleMcpSearchRegionTest, AVCStyleMcpSearchRegionTest) {
    run_Mcp_test();
};

INSTANTIATE_TEST_CASE_P(AVCMCPSearchRegion, AVCStyleMcpSearchRegionTest,
                        ::testing::Combine(::testing::ValuesIn(TEST_AREAS),
                                           ::testing::ValuesIn(TEST_PATTERNS)));

class AVCStyleMcpPUTestHelper
    : public AVCStyleMcpTestBase,
      public ::testing::WithParamInterface<TestPUParam> {
  public:
    AVCStyleMcpPUTestHelper()
        : AVCStyleMcpTestBase(std::get<0>(TEST_GET_PARAM(0)),
                              std::get<1>(TEST_GET_PARAM(0)), TEST_GET_PARAM(1),
                              true) {
    }
};

TEST_P(AVCStyleMcpPUTestHelper, AVCStyleMcpPUTestHelper) {
    run_Mcp_test_helper();
};

INSTANTIATE_TEST_CASE_P(AVCMCPPU_HELPER, AVCStyleMcpPUTestHelper,
                        ::testing::Combine(::testing::ValuesIn(TEST_PU_HELPER),
                                           ::testing::ValuesIn(TEST_PATTERNS)));

class AVCStyleMcpSearchRegionTestHelper
    : public AVCStyleMcpTestBase,
      public ::testing::WithParamInterface<TestSearchRegionParam> {
  public:
    AVCStyleMcpSearchRegionTestHelper()
        : AVCStyleMcpTestBase(std::get<0>(TEST_GET_PARAM(0)),
                              std::get<1>(TEST_GET_PARAM(0)), TEST_GET_PARAM(1),
                              false) {
    }
};

TEST_P(AVCStyleMcpSearchRegionTestHelper, AVCStyleMcpSearchRegionTestHelper) {
    run_Mcp_test_helper();
};

INSTANTIATE_TEST_CASE_P(AVCMCPSearchRegionHelper,
                        AVCStyleMcpSearchRegionTestHelper,
                        ::testing::Combine(::testing::ValuesIn(TEST_PU_HELPER),
                                           ::testing::ValuesIn(TEST_PATTERNS)));

}  // namespace
