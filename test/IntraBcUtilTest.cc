/*
 * Copyright(c) 2019 Netflix, Inc.
 * SPDX - License - Identifier: BSD - 2 - Clause - Patent
 */

/******************************************************************************
 * @file IntraBcUtilTest.cc
 *
 * @brief Unit test of Intra BC utility:
 * - av1_is_dv_valid
 *
 * @author Cidana-Edmond
 *
 ******************************************************************************/

#include "gtest/gtest.h"
// workaround to eliminate the compiling warning on linux
// The macro will conflict with definition in gtest.h
#ifdef __USE_GNU
#undef __USE_GNU  // defined in EbThreads.h
#endif
#ifdef _GNU_SOURCE
#undef _GNU_SOURCE  // defined in EbThreads.h
#endif
#include "EbDefinitions.h"
#include "EbAdaptiveMotionVectorPrediction.h"
#include "util.h"

namespace {
using std::make_tuple;

typedef std::tuple<MV,            /**< DV */
                   int,           /**< offset of MI row */
                   int,           /**< offset of MI colum */
                   BlockSize, int /**< the result of dv validation function */
                   >
    DvValidationParam;

static const int kSubPelScale = 8;
static const int kTileMaxMibWidth = 8;
static const DvValidationParam dv_validation_params[] = {
    DvValidationParam({0, 0}, 0, 0, BLOCK_128X128, 0),
    DvValidationParam({0, 0}, 0, 0, BLOCK_64X64, 0),
    DvValidationParam({0, 0}, 0, 0, BLOCK_32X32, 0),
    DvValidationParam({0, 0}, 0, 0, BLOCK_16X16, 0),
    DvValidationParam({0, 0}, 0, 0, BLOCK_8X8, 0),
    DvValidationParam({0, 0}, 0, 0, BLOCK_4X4, 0),
    DvValidationParam({-MAX_SB_SIZE * kSubPelScale, -MAX_SB_SIZE* kSubPelScale},
                      MAX_SB_SIZE / MI_SIZE, MAX_SB_SIZE / MI_SIZE, BLOCK_16X16,
                      1),
    DvValidationParam({0, -MAX_SB_SIZE* kSubPelScale}, MAX_SB_SIZE / MI_SIZE,
                      MAX_SB_SIZE / MI_SIZE, BLOCK_16X16, 0),
    DvValidationParam({-MAX_SB_SIZE * kSubPelScale, 0}, MAX_SB_SIZE / MI_SIZE,
                      MAX_SB_SIZE / MI_SIZE, BLOCK_16X16, 1),
    DvValidationParam({MAX_SB_SIZE * kSubPelScale, 0}, MAX_SB_SIZE / MI_SIZE,
                      MAX_SB_SIZE / MI_SIZE, BLOCK_16X16, 0),
    DvValidationParam({0, MAX_SB_SIZE* kSubPelScale}, MAX_SB_SIZE / MI_SIZE,
                      MAX_SB_SIZE / MI_SIZE, BLOCK_16X16, 0),
    DvValidationParam({-32 * kSubPelScale, -32 * kSubPelScale},
                      MAX_SB_SIZE / MI_SIZE, MAX_SB_SIZE / MI_SIZE, BLOCK_32X32,
                      1),
    DvValidationParam({-32 * kSubPelScale, -32 * kSubPelScale}, 32 / MI_SIZE,
                      32 / MI_SIZE, BLOCK_32X32, 0),
    DvValidationParam(
        {-32 * kSubPelScale - kSubPelScale / 2, -32 * kSubPelScale},
        MAX_SB_SIZE / MI_SIZE, MAX_SB_SIZE / MI_SIZE, BLOCK_32X32, 0),
    DvValidationParam({-33 * kSubPelScale, -32 * kSubPelScale},
                      MAX_SB_SIZE / MI_SIZE, MAX_SB_SIZE / MI_SIZE, BLOCK_32X32,
                      1),
    DvValidationParam(
        {-32 * kSubPelScale, -32 * kSubPelScale - kSubPelScale / 2},
        MAX_SB_SIZE / MI_SIZE, MAX_SB_SIZE / MI_SIZE, BLOCK_32X32, 0),
    DvValidationParam({-32 * kSubPelScale, -33 * kSubPelScale},
                      MAX_SB_SIZE / MI_SIZE, MAX_SB_SIZE / MI_SIZE, BLOCK_32X32,
                      1),
    DvValidationParam({-MAX_SB_SIZE * kSubPelScale, -MAX_SB_SIZE* kSubPelScale},
                      MAX_SB_SIZE / MI_SIZE, MAX_SB_SIZE / MI_SIZE,
                      BLOCK_LARGEST, 1),
    DvValidationParam(
        {-(MAX_SB_SIZE + 1) * kSubPelScale, -MAX_SB_SIZE* kSubPelScale},
        MAX_SB_SIZE / MI_SIZE, MAX_SB_SIZE / MI_SIZE, BLOCK_LARGEST, 0),
    DvValidationParam(
        {-MAX_SB_SIZE * kSubPelScale, -(MAX_SB_SIZE + 1) * kSubPelScale},
        MAX_SB_SIZE / MI_SIZE, MAX_SB_SIZE / MI_SIZE, BLOCK_LARGEST, 0),
    DvValidationParam(
        {-(MAX_SB_SIZE - 1) * kSubPelScale, -MAX_SB_SIZE* kSubPelScale},
        MAX_SB_SIZE / MI_SIZE, MAX_SB_SIZE / MI_SIZE, BLOCK_LARGEST, 0),
    DvValidationParam(
        {-MAX_SB_SIZE * kSubPelScale, -(MAX_SB_SIZE - 1) * kSubPelScale},
        MAX_SB_SIZE / MI_SIZE, MAX_SB_SIZE / MI_SIZE, BLOCK_LARGEST, 1),
    DvValidationParam(
        {-(MAX_SB_SIZE - 1) * kSubPelScale, -(MAX_SB_SIZE - 1) * kSubPelScale},
        MAX_SB_SIZE / MI_SIZE, MAX_SB_SIZE / MI_SIZE, BLOCK_LARGEST, 0),
    DvValidationParam({-MAX_SB_SIZE * kSubPelScale, MAX_SB_SIZE* kSubPelScale},
                      MAX_SB_SIZE / MI_SIZE, MAX_SB_SIZE / MI_SIZE,
                      BLOCK_LARGEST, 0),
    DvValidationParam({-MAX_SB_SIZE * kSubPelScale,
                       (kTileMaxMibWidth - 2) * MAX_SB_SIZE* kSubPelScale},
                      MAX_SB_SIZE / MI_SIZE, MAX_SB_SIZE / MI_SIZE,
                      BLOCK_LARGEST, 0),
    DvValidationParam(
        {-MAX_SB_SIZE * kSubPelScale,
         ((kTileMaxMibWidth - 2) * MAX_SB_SIZE + 1) * kSubPelScale},
        MAX_SB_SIZE / MI_SIZE, MAX_SB_SIZE / MI_SIZE, BLOCK_LARGEST, 0),
};

class DvValiationTest : public ::testing::TestWithParam<DvValidationParam> {
  public:
    void is_dv_validate() {
        MacroBlockD xd;
        memset(&xd, 0, sizeof(xd));
        xd.tile.mi_row_start = 8 * MAX_MIB_SIZE;
        xd.tile.mi_row_end = 16 * MAX_MIB_SIZE;
        xd.tile.mi_col_start = 24 * MAX_MIB_SIZE;
        xd.tile.mi_col_end =
            xd.tile.mi_col_start + kTileMaxMibWidth * MAX_MIB_SIZE;
        xd.plane[1].subsampling_x = 1;
        xd.plane[1].subsampling_y = 1;
        xd.plane[2].subsampling_x = 1;
        xd.plane[2].subsampling_y = 1;

        ASSERT_EQ(TEST_GET_PARAM(4),
                  av1_is_dv_valid(TEST_GET_PARAM(0),
                                  &xd,
                                  xd.tile.mi_row_start + TEST_GET_PARAM(1),
                                  xd.tile.mi_col_start + TEST_GET_PARAM(2),
                                  TEST_GET_PARAM(3),
                                  MAX_MIB_SIZE_LOG2));
    }
};

TEST_P(DvValiationTest, IsDvValidate) {
    is_dv_validate();
}

INSTANTIATE_TEST_CASE_P(AV1, DvValiationTest,
                        ::testing::ValuesIn(dv_validation_params));

}  // namespace
