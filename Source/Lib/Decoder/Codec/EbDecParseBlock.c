/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

/*
* Copyright (c) 2016, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at www.aomedia.org/license/software. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at www.aomedia.org/license/patent.
*/

#include "EbDefinitions.h"
#include "EbPictureBufferDesc.h"

#include "EbSvtAv1Dec.h"
#include "EbDecHandle.h"
#include "EbDecBitReader.h"
#include "EbObuParse.h"

#include "EbDecParseHelper.h"
#include "EbTransforms.h"

#include "EbDecNbr.h"

#if ENABLE_ENTROPY_TRACE
FILE* temp_fp;
#endif

#define MV_BORDER (16 << 3)  // Allow 16 pels in 1/8th pel units
#define MAX_DIFFWTD_MASK_BITS 1
#define READ_REF_BIT(pname) \
  svt_read_symbol(r, get_pred_cdf_##pname(xd), 2, ACCT_STR)
#define SQR_BLOCK_SIZES 6

typedef AomCdfProb(*base_cdf_arr)[CDF_SIZE(4)];
typedef AomCdfProb(*br_cdf_arr)[CDF_SIZE(BR_CDF_SIZE)];

static const TxSize tx_mode_to_biggest_tx_size[TX_MODES] = {
  TX_4X4,    // ONLY_4X4
  TX_64X64,  // TX_MODE_LARGEST
  TX_64X64,  // TX_MODE_SELECT
};

typedef struct txb_ctx {
    int txb_skip_ctx;
    int dc_sign_ctx;
} TXB_CTX;

static const int16_t k_eob_group_start[12] = { 0,  1,  2,  3,   5,   9,
                                        17, 33, 65, 129, 257, 513 };
static const int16_t k_eob_offset_bits[12] = { 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

extern const int8_t av1_nz_map_ctx_offset_4x4[16];/* = {
  0, 1, 6, 6, 1, 6, 6, 21, 6, 6, 21, 21, 6, 21, 21, 21,
};*/

extern const int8_t av1_nz_map_ctx_offset_8x8[64];/* = {
  0,  1,  6,  6,  21, 21, 21, 21, 1,  6,  6,  21, 21, 21, 21, 21,
  6,  6,  21, 21, 21, 21, 21, 21, 6,  21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
};*/

extern const int8_t av1_nz_map_ctx_offset_16x16[256];/* = {
  0,  1,  6,  6,  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 1,  6,  6,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 6,  6,  21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 6,  21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21,
};*/

extern const int8_t av1_nz_map_ctx_offset_32x32[1024];/* = {
  0,  1,  6,  6,  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 1,  6,  6,  21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 6,  6,  21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 6,  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
};*/

extern const int8_t av1_nz_map_ctx_offset_8x4[32];/* = {
  0,  16, 6,  6,  21, 21, 21, 21, 16, 16, 6,  21, 21, 21, 21, 21,
  16, 16, 21, 21, 21, 21, 21, 21, 16, 16, 21, 21, 21, 21, 21, 21,
};*/

extern const int8_t av1_nz_map_ctx_offset_8x16[128];/* = {
  0,  11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 6,  6,  21,
  21, 21, 21, 21, 21, 6,  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
};*/

extern const int8_t av1_nz_map_ctx_offset_16x8[128];/* = {
  0,  16, 6,  6,  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 16, 16, 6,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 16, 16, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 16, 16, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 16, 16, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 16, 16, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 16, 16, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 16, 16,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
};*/

extern const int8_t av1_nz_map_ctx_offset_16x32[512];/* = {
  0,  11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11,
  11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 6,  6,  21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 6,  21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
};*/

extern const int8_t av1_nz_map_ctx_offset_32x16[512];/* = {
  0,  16, 6,  6,  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 16, 16, 6,  21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 16, 16, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 16, 16, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 16, 16, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 16, 16, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 16, 16, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 16, 16, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 16, 16, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 16, 16, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 16, 16, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 16, 16, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 16, 16, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 16, 16,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 16, 16, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 16, 16, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
};*/

extern const int8_t av1_nz_map_ctx_offset_32x64[1024];/* = {
  0,  11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11,
  11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11,
  11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11,
  11, 11, 11, 11, 11, 11, 11, 6,  6,  21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 6,  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
};*/

extern const int8_t av1_nz_map_ctx_offset_64x32[1024];/* = {
  0,  16, 6,  6,  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 16, 16, 6,  21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 16, 16, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 16, 16, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 16, 16, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 16, 16, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 16, 16, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 16, 16, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 16, 16, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 16, 16, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 16, 16, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 16, 16, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 16, 16, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 16, 16,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 16, 16, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 16, 16, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 16,
  16, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 16, 16, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 16, 16, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  16, 16, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 16, 16, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 16, 16, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 16, 16, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 16, 16, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 16, 16, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 16, 16, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 16, 16, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 16, 16, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 16, 16, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 16, 16, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 16, 16, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 16, 16, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
};*/

extern const int8_t av1_nz_map_ctx_offset_4x16[64];/* = {
  0,  11, 11, 11, 11, 11, 11, 11, 6,  6,  21, 21, 6,  21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
};*/

extern const int8_t av1_nz_map_ctx_offset_16x4[64];/* = {
  0,  16, 6,  6,  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  16, 16, 6,  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  16, 16, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  16, 16, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
};*/

extern const int8_t av1_nz_map_ctx_offset_8x32[256];/* = {
  0,  11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 6,  6,  21,
  21, 21, 21, 21, 21, 6,  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21,
};*/

extern const int8_t av1_nz_map_ctx_offset_32x8[256];/* = {
  0,  16, 6,  6,  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 16, 16, 6,  21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 16, 16, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 16, 16, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 16, 16, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 16, 16, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 16, 16, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 16, 16, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21,
};*/

static const int8_t *av1_nz_map_ctx_offset[19] = {
  av1_nz_map_ctx_offset_4x4,    // TX_4x4
  av1_nz_map_ctx_offset_8x8,    // TX_8x8
  av1_nz_map_ctx_offset_16x16,  // TX_16x16
  av1_nz_map_ctx_offset_32x32,  // TX_32x32
  av1_nz_map_ctx_offset_32x32,  // TX_32x32
  av1_nz_map_ctx_offset_4x16,   // TX_4x8
  av1_nz_map_ctx_offset_8x4,    // TX_8x4
  av1_nz_map_ctx_offset_8x32,   // TX_8x16
  av1_nz_map_ctx_offset_16x8,   // TX_16x8
  av1_nz_map_ctx_offset_16x32,  // TX_16x32
  av1_nz_map_ctx_offset_32x16,  // TX_32x16
  av1_nz_map_ctx_offset_32x64,  // TX_32x64
  av1_nz_map_ctx_offset_64x32,  // TX_64x32
  av1_nz_map_ctx_offset_4x16,   // TX_4x16
  av1_nz_map_ctx_offset_16x4,   // TX_16x4
  av1_nz_map_ctx_offset_8x32,   // TX_8x32
  av1_nz_map_ctx_offset_32x8,   // TX_32x8
  av1_nz_map_ctx_offset_16x32,  // TX_16x64
  av1_nz_map_ctx_offset_64x32,  // TX_64x16
};

#define NZ_MAP_CTX_0 SIG_COEF_CONTEXTS_2D
#define NZ_MAP_CTX_5 (NZ_MAP_CTX_0 + 5)
#define NZ_MAP_CTX_10 (NZ_MAP_CTX_0 + 10)

const int nz_map_ctx_offset_1d[32] = {
  NZ_MAP_CTX_0,  NZ_MAP_CTX_5,  NZ_MAP_CTX_10, NZ_MAP_CTX_10, NZ_MAP_CTX_10,
  NZ_MAP_CTX_10, NZ_MAP_CTX_10, NZ_MAP_CTX_10, NZ_MAP_CTX_10, NZ_MAP_CTX_10,
  NZ_MAP_CTX_10, NZ_MAP_CTX_10, NZ_MAP_CTX_10, NZ_MAP_CTX_10, NZ_MAP_CTX_10,
  NZ_MAP_CTX_10, NZ_MAP_CTX_10, NZ_MAP_CTX_10, NZ_MAP_CTX_10, NZ_MAP_CTX_10,
  NZ_MAP_CTX_10, NZ_MAP_CTX_10, NZ_MAP_CTX_10, NZ_MAP_CTX_10, NZ_MAP_CTX_10,
  NZ_MAP_CTX_10, NZ_MAP_CTX_10, NZ_MAP_CTX_10, NZ_MAP_CTX_10, NZ_MAP_CTX_10,
  NZ_MAP_CTX_10, NZ_MAP_CTX_10,
};

static const MV kZeroMv = { 0, 0 };

static INLINE int get_nz_mag(const uint8_t *const levels,
    const int bwl, const TxClass tx_class) {
    int mag;

#define CLIP_MAX3(x) x > 3 ? 3 : x

    // Note: AOMMIN(level, 3) is useless for decoder since level < 3.
    mag = CLIP_MAX3(levels[1]);                         // { 0, 1 }
    mag += CLIP_MAX3(levels[(1 << bwl) + TX_PAD_HOR]);  // { 1, 0 }

    if (tx_class == TX_CLASS_2D) {
        mag += CLIP_MAX3(levels[(1 << bwl) + TX_PAD_HOR + 1]);          // { 1, 1 }
        mag += CLIP_MAX3(levels[2]);                                    // { 0, 2 }
        mag += CLIP_MAX3(levels[(2 << bwl) + (2 << TX_PAD_HOR_LOG2)]);  // { 2, 0 }
    }
    else if (tx_class == TX_CLASS_VERT) {
        mag += CLIP_MAX3(levels[(2 << bwl) + (2 << TX_PAD_HOR_LOG2)]);  // { 2, 0 }
        mag += CLIP_MAX3(levels[(3 << bwl) + (3 << TX_PAD_HOR_LOG2)]);  // { 3, 0 }
        mag += CLIP_MAX3(levels[(4 << bwl) + (4 << TX_PAD_HOR_LOG2)]);  // { 4, 0 }
    }
    else {
        mag += CLIP_MAX3(levels[2]);  // { 0, 2 }
        mag += CLIP_MAX3(levels[3]);  // { 0, 3 }
        mag += CLIP_MAX3(levels[4]);  // { 0, 4 }
    }

    return mag;
}

static INLINE int get_nz_map_ctx_from_stats(
    const int stats,
    const int coeff_idx,  // raster order
    const int bwl, const TxSize tx_size, const TxClass tx_class) {
    if ((tx_class | coeff_idx) == 0) return 0;
    int ctx = (stats + 1) >> 1;
    ctx = AOMMIN(ctx, 4);
    switch (tx_class) {
    case TX_CLASS_2D: {
        return ctx + av1_nz_map_ctx_offset[tx_size][coeff_idx];
    }
    case TX_CLASS_HORIZ: {
        const int row = coeff_idx >> bwl;
        const int col = coeff_idx - (row << bwl);
        return ctx + nz_map_ctx_offset_1d[col];
    }
    case TX_CLASS_VERT: {
        const int row = coeff_idx >> bwl;
        return ctx + nz_map_ctx_offset_1d[row];
    }
    default: break;
    }
    return 0;
}

void palette_mode_info(/*EbDecHandle *dec_handle,
    int mi_row, int mi_col, SvtReader *r*/) {
    //TO-DO
    assert(0);
}

void filter_intra_mode_info(EbDecHandle *dec_handle,
    PartitionInfo_t *xd, SvtReader *r) {
    ModeInfo_t *const mbmi = xd->mi;
    FilterIntraModeInfo_t *filter_intra_mode_info =
        &mbmi->filter_intra_mode_info;
    ParseCtxt *parse_ctxt = (ParseCtxt *)dec_handle->pv_parse_ctxt;
    FRAME_CONTEXT *frm_ctx = &parse_ctxt->frm_ctx[0];

    if (filter_intra_allowed(dec_handle, mbmi)) {
        filter_intra_mode_info->use_filter_intra = svt_read_symbol(
            r, frm_ctx->filter_intra_cdfs[mbmi->sb_type], 2, ACCT_STR);
        if (filter_intra_mode_info->use_filter_intra) {
            filter_intra_mode_info->filter_intra_mode = svt_read_symbol(
                r, frm_ctx->filter_intra_mode_cdf, FILTER_INTRA_MODES, ACCT_STR);
        }
    }
    else
        filter_intra_mode_info->use_filter_intra = 0;
}

uint8_t read_cfl_alphas(FRAME_CONTEXT *ec_ctx, SvtReader *r,
    uint8_t *signs_out) {
    const int joint_sign =
        svt_read_symbol(r, ec_ctx->cfl_sign_cdf, CFL_JOINT_SIGNS, "cfl:signs");
    uint8_t idx = 0;
    // Magnitudes are only coded for nonzero values
    if (CFL_SIGN_U(joint_sign) != CFL_SIGN_ZERO) {
        AomCdfProb *cdf_u = ec_ctx->cfl_alpha_cdf[CFL_CONTEXT_U(joint_sign)];
        idx = svt_read_symbol(r, cdf_u, CFL_ALPHABET_SIZE, "cfl:alpha_u")
            << CFL_ALPHABET_SIZE_LOG2;
    }
    if (CFL_SIGN_V(joint_sign) != CFL_SIGN_ZERO) {
        AomCdfProb *cdf_v = ec_ctx->cfl_alpha_cdf[CFL_CONTEXT_V(joint_sign)];
        idx += svt_read_symbol(r, cdf_v, CFL_ALPHABET_SIZE, "cfl:alpha_v");
    }
    *signs_out = joint_sign;
    return idx;
}

void read_cdef(EbDecHandle *dec_handle, SvtReader *r, PartitionInfo_t *xd,
    int mi_col, int mi_row, int8_t *cdef_strength) {
    ModeInfo_t *const mbmi = xd->mi;

    if (dec_handle->frame_header.coded_lossless || dec_handle->frame_header.allow_intrabc
        || !dec_handle->seq_header.enable_cdef || mbmi->skip_mode)
        return;

    if (!(mi_col & (dec_handle->seq_header.sb_mi_size - 1)) &&
        !(mi_row & (dec_handle->seq_header.sb_mi_size - 1))) {  // Top left?
        cdef_strength[0] = cdef_strength[1] = cdef_strength[2] =
            cdef_strength[3] = -1;
    }
    // Read CDEF param at the first non-skip coding block
    const int mask = (1 << (6 - MI_SIZE_LOG2));
    const int index = dec_handle->seq_header.sb_size == BLOCK_128X128
        ? !!(mi_col & mask) + 2 * !!(mi_row & mask)
        : 0;
    cdef_strength[index] = cdef_strength[index] == -1 && !mbmi->skip
        ? svt_read_literal(r, dec_handle->frame_header.CDEF_params.cdef_bits, ACCT_STR)
        : cdef_strength[index];
}

int read_delta_qindex(EbDecHandle *dec_handle,
    SvtReader *r, ModeInfo_t *const mbmi,
    int mi_col, int mi_row) {
    int sign, abs, reduced_delta_qindex = 0;
    BlockSize bsize = mbmi->sb_type;
    const int b_col = mi_col & (dec_handle->seq_header.sb_mi_size - 1);
    const int b_row = mi_row & (dec_handle->seq_header.sb_mi_size - 1);
    const int read_delta_q_flag = (b_col == 0 && b_row == 0);

    if ((bsize != dec_handle->seq_header.sb_size || mbmi->skip == 0) &&
        read_delta_q_flag) {
        ParseCtxt *parse_ctxt = (ParseCtxt *)dec_handle->pv_parse_ctxt;
        abs = svt_read_symbol(r, parse_ctxt->frm_ctx[0].delta_q_cdf, DELTA_Q_PROBS + 1, ACCT_STR);
        const int smallval = (abs < DELTA_Q_SMALL);

        if (!smallval) {
            const int rem_bits = svt_read_literal(r, 3, ACCT_STR) + 1;
            const int thr = (1 << rem_bits) + 1;
            abs = svt_read_literal(r, rem_bits, ACCT_STR) + thr;
        }

        if (abs)
            sign = svt_read_bit(r, ACCT_STR);
        else
            sign = 1;
        reduced_delta_qindex = sign ? -abs : abs;
    }
    return reduced_delta_qindex;
}

int read_delta_lflevel(EbDecHandle *dec_handle, SvtReader *r,
    AomCdfProb *cdf,
    ModeInfo_t *mbmi, int mi_col,
    int mi_row) {
    int reduced_delta_lflevel = 0;
    const BlockSize bsize = mbmi->sb_type;
    const int b_col = mi_col & (dec_handle->seq_header.sb_mi_size - 1);
    const int b_row = mi_row & (dec_handle->seq_header.sb_mi_size - 1);
    const int read_delta_lf_flag = (b_col == 0 && b_row == 0);

    if ((bsize != dec_handle->seq_header.sb_size || mbmi->skip == 0) &&
        read_delta_lf_flag) {
        int abs = svt_read_symbol(r, cdf, DELTA_LF_PROBS + 1, ACCT_STR);
        const int smallval = (abs < DELTA_LF_SMALL);
        if (!smallval) {
            const int rem_bits = svt_read_literal(r, 3, ACCT_STR) + 1;
            const int thr = (1 << rem_bits) + 1;
            abs = svt_read_literal(r, rem_bits, ACCT_STR) + thr;
        }
        const int sign = abs ? svt_read_bit(r, ACCT_STR) : 1;
        reduced_delta_lflevel = sign ? -abs : abs;
    }
    return reduced_delta_lflevel;
}

int seg_feature_active(SegmentationParams *seg, int segment_id,
    SEG_LVL_FEATURES feature_id)
{
    return seg->segmentation_enabled && seg->feature_data[segment_id][feature_id];
}

int read_skip(EbDecHandle *dec_handle, PartitionInfo_t *xd, int segment_id,
    SvtReader *r)
{
    ParseCtxt *parse_ctxt = (ParseCtxt *)dec_handle->pv_parse_ctxt;
    uint8_t segIdPreSkip = dec_handle->frame_header.segmentation_params.seg_id_pre_skip;
    if (segIdPreSkip && seg_feature_active(&dec_handle->frame_header.segmentation_params,
        segment_id, SEG_LVL_SKIP))
        return 1;
    else {
        const int above_skip = xd->above_mbmi ? xd->above_mbmi->skip : 0;
        const int left_skip = xd->left_mbmi ? xd->left_mbmi->skip : 0;
        int ctx = above_skip + left_skip;
        return svt_read_symbol(r, parse_ctxt->frm_ctx[0].skip_cdfs[ctx], 2, ACCT_STR);
    }
}

int read_skip_mode(EbDecHandle *dec_handle, PartitionInfo_t *xd, int segment_id,
    SvtReader *r)
{
    if (seg_feature_active(&dec_handle->frame_header.segmentation_params, segment_id, SEG_LVL_SKIP) ||
        seg_feature_active(&dec_handle->frame_header.segmentation_params, segment_id, SEG_LVL_REF_FRAME) ||
        seg_feature_active(&dec_handle->frame_header.segmentation_params, segment_id, SEG_LVL_GLOBALMV) ||
        !dec_handle->frame_header.skip_mode_params.skip_mode_present ||
        block_size_wide[xd->mi->sb_type] < 8 ||
        block_size_high[xd->mi->sb_type] < 8) {
        return 0;
    }
    ParseCtxt *parse_ctxt = (ParseCtxt *)dec_handle->pv_parse_ctxt;
    int above_skip_mode = xd->above_mbmi ? xd->above_mbmi->skip_mode : 0;
    int left_skip_mode = xd->left_mbmi ? xd->left_mbmi->skip_mode : 0;
    int ctx = above_skip_mode + left_skip_mode;
    return svt_read_symbol(r, parse_ctxt->frm_ctx[0].skip_mode_cdfs[ctx], 2, ACCT_STR);
}

// If delta q is present, reads delta_q index.
// Also reads delta_q loop filter levels, if present.
static void read_delta_params(EbDecHandle *dec_handle, SvtReader *r,
                              PartitionInfo_t *xd,
                              const int mi_row, const int mi_col)
{
    ParseCtxt       *parse_ctxt     = (ParseCtxt *)dec_handle->pv_parse_ctxt;
    DeltaQParams    *delta_q_params = &dec_handle->frame_header.delta_q_params;
    DeltaLFParams   *delta_lf_params = &dec_handle->frame_header.delta_lf_params;
    SBInfo          *sb_info = xd->sb_info;

    if (delta_q_params->delta_q_present) {
        ModeInfo_t *const mbmi = &xd->mi[0];
        int current_qindex;
        int base_qindex = dec_handle->frame_header.quantization_params.base_q_idx;

        sb_info->sb_delta_q[0] = read_delta_qindex(dec_handle, r, mbmi, mi_col, mi_row) *
                                delta_q_params->delta_q_res;
        current_qindex = base_qindex + sb_info->sb_delta_q[0];
        /* Normative: Clamp to [1,MAXQ] to not interfere with lossless mode */
        /* TODO: Should add cur q idx and remove this repetitive? */
        sb_info->sb_delta_q[0] = clamp(current_qindex, 1, MAXQ) - base_qindex;

        FRAME_CONTEXT *const ec_ctx = &parse_ctxt->frm_ctx[0];

        if (delta_lf_params->delta_lf_present) {
            if (delta_lf_params->delta_lf_multi) {
                EbColorConfig *color_info = &dec_handle->seq_header.color_config;
                int num_planes = color_info->mono_chrome ? 1 : MAX_MB_PLANE;
                const int frame_lf_count =
                    num_planes > 1 ? FRAME_LF_COUNT : FRAME_LF_COUNT - 2;
                for (int lf_id = 0; lf_id < frame_lf_count; ++lf_id) {
                    int tmp_lvl;
                    int base_lvl = dec_handle->frame_header.loop_filter_params.loop_filter_level[lf_id];
                    sb_info->sb_delta_lf[lf_id] =
                        read_delta_lflevel(dec_handle, r, ec_ctx->delta_lf_multi_cdf[lf_id],
                                           mbmi, mi_col, mi_row) *
                                           delta_lf_params->delta_lf_res;
                    tmp_lvl = base_lvl + sb_info->sb_delta_lf[lf_id];
                    tmp_lvl = clamp(tmp_lvl, -MAX_LOOP_FILTER, MAX_LOOP_FILTER);
                    sb_info->sb_delta_lf[lf_id] = tmp_lvl - base_lvl;
                }
            } else {
                int tmp_lvl;
                int base_lvl = dec_handle->frame_header.loop_filter_params.loop_filter_level[0];
                sb_info->sb_delta_lf[0] =
                                read_delta_lflevel(dec_handle, r, ec_ctx->delta_lf_cdf,
                                                    mbmi, mi_col, mi_row) *
                                delta_lf_params->delta_lf_res;
                tmp_lvl = base_lvl + sb_info->sb_delta_lf[0];
                tmp_lvl = clamp(tmp_lvl, -MAX_LOOP_FILTER, MAX_LOOP_FILTER);
                sb_info->sb_delta_lf[0] = tmp_lvl - base_lvl;
            }
        }
    }
}

int is_directional_mode(PredictionMode mode) {
    return mode >= V_PRED && mode <= D67_PRED;
}

int intra_angle_info(SvtReader *r, AomCdfProb *cdf, PredictionMode mode, BlockSize bsize) {
    int angleDeltaY = 0;
    if (use_angle_delta(bsize) && is_directional_mode(mode))
    {
        const int sym = svt_read_symbol(r, cdf, 2 * MAX_ANGLE_DELTA + 1, ACCT_STR);
        angleDeltaY = sym - MAX_ANGLE_DELTA;
    }
    return angleDeltaY;
}

int get_segment_id(FrameHeader *frm_info, uint8_t *segment_ids,
    BlockSize bsize, uint32_t mi_row, uint32_t mi_col) {
    const int mi_offset = mi_row * frm_info->mi_cols + mi_col;
    const uint32_t bw = mi_size_wide[bsize];
    const uint32_t bh = mi_size_high[bsize];
    const int xmis = AOMMIN(frm_info->mi_cols - mi_col, bw);
    const int ymis = AOMMIN(frm_info->mi_rows - mi_row, bh);
    int x, y, segment_id = MAX_SEGMENTS;

    for (y = 0; y < ymis; ++y)
        for (x = 0; x < xmis; ++x)
            segment_id =
            AOMMIN(segment_id, segment_ids[mi_offset + y * frm_info->mi_cols + x]);

    assert(segment_id >= 0 && segment_id < MAX_SEGMENTS);
    return segment_id;
}

static int read_segment_id(EbDecHandle *dec_handle, PartitionInfo_t *xd, uint32_t mi_row,
    uint32_t mi_col, SvtReader *r, int skip)
{
    int cdf_num = 0;
    ParseCtxt *parse_ctxt = (ParseCtxt *)dec_handle->pv_parse_ctxt;

    int prev_ul = -1;  // top left segment_id
    int prev_l = -1;   // left segment_id
    int prev_u = -1;   // top segment_id
    int pred = -1;
    if ((xd->up_available) && (xd->left_available)) {
        prev_ul = get_segment_id(&dec_handle->frame_header, parse_ctxt->parse_nbr4x4_ctxt.segment_maps, BLOCK_4X4, mi_row - 1,
            mi_col - 1);
    }
    if (xd->up_available) {
        prev_u = get_segment_id(&dec_handle->frame_header, parse_ctxt->parse_nbr4x4_ctxt.segment_maps, BLOCK_4X4, mi_row - 1,
            mi_col);
    }
    if (xd->left_available) {
        prev_l = get_segment_id(&dec_handle->frame_header, parse_ctxt->parse_nbr4x4_ctxt.segment_maps, BLOCK_4X4, mi_row,
            mi_col - 1);
    }
    if ((prev_ul == prev_u) && (prev_ul == prev_l))
        cdf_num = 2;
    else if ((prev_ul == prev_u) || (prev_ul == prev_l) || (prev_u == prev_l))
        cdf_num = 1;

    // If 2 or more are identical returns that as predictor, otherwise prev_l.
    if (prev_u == -1)  // edge case
        pred = prev_l == -1 ? 0 : prev_l;
    if (prev_l == -1)  // edge case
        pred = prev_u;
    pred = (prev_ul == prev_u) ? prev_u : prev_l;
    return pred;

    if (skip) return pred;

    FRAME_CONTEXT *ec_ctx = &parse_ctxt->frm_ctx[0];
    SegmentationParams *seg = &(dec_handle->frame_header.segmentation_params);

    struct segmentation_probs *segp = &ec_ctx->seg;
    AomCdfProb *pred_cdf = segp->spatial_pred_seg_cdf[cdf_num];

    int coded_id = svt_read_symbol(r, pred_cdf, MAX_SEGMENTS, ACCT_STR);
    return neg_deinterleave(coded_id, pred, seg->last_active_seg_id + 1);
}

int intra_segment_id(EbDecHandle *dec_handle, PartitionInfo_t *xd, int mi_row, int mi_col,
    int bsize, SvtReader *r, int skip)
{
    SegmentationParams *seg = &dec_handle->frame_header.segmentation_params;
    int segment_id = 0;

    if (seg->segmentation_enabled) {
        const int mi_offset = mi_row * dec_handle->frame_header.mi_cols + mi_col;
        const int bw = mi_size_wide[bsize];
        const int bh = mi_size_high[bsize];
        const int x_mis = AOMMIN((int32_t)(dec_handle->frame_header.mi_cols - mi_col), bw);
        const int y_mis = AOMMIN((int32_t)(dec_handle->frame_header.mi_rows - mi_row), bh);
        segment_id = read_segment_id(dec_handle, xd, mi_row, mi_col, r, skip);
        set_segment_id(dec_handle, mi_offset, x_mis, y_mis, segment_id);
    }
    return segment_id;
}

void intra_frame_mode_info(EbDecHandle *dec_handle,
    PartitionInfo_t *xd, int mi_row,
    int mi_col, SvtReader *r, int8_t *cdef_strength)
{
    ParseCtxt *parse_ctxt = (ParseCtxt *)dec_handle->pv_parse_ctxt;
    ModeInfo_t *const mbmi = xd->mi;
    const ModeInfo_t *above_mi = xd->above_mbmi;
    const ModeInfo_t *left_mi = xd->left_mbmi;
    const BlockSize bsize = mbmi->sb_type;
    SegmentationParams *const seg = &dec_handle->frame_header.segmentation_params;
    EbColorConfig color_config = dec_handle->seq_header.color_config;
    uint8_t     *lossless_array = &dec_handle->frame_header.lossless_array[0];

    if (seg->seg_id_pre_skip)
        mbmi->segment_id =
        intra_segment_id(dec_handle, xd, mi_row, mi_col, bsize, r, 0);

    mbmi->skip = read_skip(dec_handle, xd, mbmi->segment_id, r);

    if (!seg->seg_id_pre_skip)
        mbmi->segment_id =
        intra_segment_id(dec_handle, xd, mi_row, mi_col, bsize, r, mbmi->skip);

    read_cdef(dec_handle, r, xd, mi_col, mi_row, cdef_strength);

    read_delta_params(dec_handle, r, xd, mi_row, mi_col);

    mbmi->ref_frame[0] = INTRA_FRAME;
    mbmi->ref_frame[1] = NONE_FRAME;

    mbmi->use_intrabc = 0;
    if (allow_intrabc(dec_handle))
        mbmi->use_intrabc = svt_read_symbol(r, parse_ctxt->frm_ctx[0].intrabc_cdf, 2, ACCT_STR);

    if (mbmi->use_intrabc) {
        mbmi->mode = DC_PRED;
        mbmi->uv_mode = UV_DC_PRED;
        mbmi->motion_mode = SIMPLE_TRANSLATION;
        mbmi->compound_mode = COMPOUND_AVERAGE;
        dec_handle->frame_header.interpolation_filter = BILINEAR;
        //find_mv_stack(0)
        //assign_mv(0)
    }
    else
    {
        mbmi->mode = read_intra_mode(r, get_y_mode_cdf(&parse_ctxt->frm_ctx[0], above_mi, left_mi));
        mbmi->angle_delta[PLANE_TYPE_Y] =
            intra_angle_info(r, &parse_ctxt->frm_ctx[0].angle_delta_cdf[mbmi->mode - V_PRED][0], mbmi->mode, bsize);
        const int has_chroma =
                dec_is_chroma_reference(mi_row, mi_col, bsize, color_config.subsampling_x,
                color_config.subsampling_y);
        if (has_chroma) {
            mbmi->uv_mode =
                read_intra_mode_uv(&parse_ctxt->frm_ctx[0], r, is_cfl_allowed(xd, &color_config, lossless_array), mbmi->mode);
            if (mbmi->uv_mode == UV_CFL_PRED)
                mbmi->cfl_alpha_idx = read_cfl_alphas(&parse_ctxt->frm_ctx[0], r, &mbmi->cfl_alpha_signs);
            mbmi->angle_delta[PLANE_TYPE_UV] =
            intra_angle_info(r, &parse_ctxt->frm_ctx[0].angle_delta_cdf[mbmi->uv_mode - V_PRED][0], dec_get_uv_mode(mbmi->uv_mode), bsize);
        }

        if (allow_palette(dec_handle->frame_header.allow_screen_content_tools, bsize))
            palette_mode_info(/*dec_handle, mi_row, mi_col, r*/);
        filter_intra_mode_info(dec_handle, xd, r);
    }
}

static INLINE int get_pred_context_seg_id(const PartitionInfo_t *xd) {
    const ModeInfo_t *const above_mi = xd->above_mbmi;
    const ModeInfo_t *const left_mi = xd->left_mbmi;
    const int above_sip = (above_mi != NULL) ? above_mi->seg_id_predicted : 0;
    const int left_sip = (left_mi != NULL) ? left_mi->seg_id_predicted : 0;

    return above_sip + left_sip;
}

void update_seg_ctx(ParseNbr4x4Ctxt *ngr_ctx, int blk_col,
    int w4, int h4, int seg_id_predicted) {
    uint8_t *const above_seg_ctx = ngr_ctx->above_seg_pred_ctx + blk_col;
    uint8_t *const left_seg_ctx = ngr_ctx->left_seg_pred_ctx;

    memset(above_seg_ctx, seg_id_predicted, w4);
    memset(left_seg_ctx, seg_id_predicted, h4);
}

int read_inter_segment_id(EbDecHandle *dec_handle, PartitionInfo_t *xd,
    uint32_t mi_row, uint32_t mi_col, int preskip,
    SvtReader *r) {
    SegmentationParams *seg = &dec_handle->frame_header.segmentation_params;
    ModeInfo_t *const mbmi = xd->mi;
    FrameHeader *frame_header = &dec_handle->frame_header;
    ParseCtxt *parse_ctxt = (ParseCtxt *)dec_handle->pv_parse_ctxt;
    const int mi_offset = mi_row * frame_header->mi_cols + mi_col;
    const uint32_t bw = mi_size_wide[mbmi->sb_type];
    const uint32_t bh = mi_size_high[mbmi->sb_type];

    const int x_mis = AOMMIN(frame_header->mi_cols - mi_col, bw);
    const int y_mis = AOMMIN(frame_header->mi_rows - mi_row, bh);

    //TODO - check neighbour access index
    int predictedSegmentId = get_segment_id(frame_header,
        parse_ctxt->parse_nbr4x4_ctxt.segment_maps, mbmi->sb_type, mi_row, mi_col);

    if (!seg->segmentation_enabled) return 0;  // Default for disabled segmentation

    if (!seg->segmentation_update_map)
        return predictedSegmentId;
    int segment_id;
    if (preskip) {
        if (!seg->seg_id_pre_skip)
          return 0;
    }
    else {
        if (mbmi->skip) {
            if (seg->segmentation_temporal_update)
                mbmi->seg_id_predicted = 0;
            update_seg_ctx(&parse_ctxt->parse_nbr4x4_ctxt, mi_col, bw, bh, mbmi->seg_id_predicted);
            segment_id = read_segment_id(dec_handle, xd, mi_row, mi_col, r, 1);
            set_segment_id(dec_handle, mi_offset, x_mis, y_mis, segment_id);
            return segment_id;
        }
    }

    if (seg->segmentation_temporal_update) {
        const int ctx = get_pred_context_seg_id(xd);
        ParseCtxt *parse_ctxt = (ParseCtxt *)dec_handle->pv_parse_ctxt;
        struct segmentation_probs *const segp = &parse_ctxt->frm_ctx[0].seg;
        AomCdfProb *pred_cdf = segp->pred_cdf[ctx];
        mbmi->seg_id_predicted = svt_read_symbol(r, pred_cdf, 2, ACCT_STR);
        if (mbmi->seg_id_predicted)
            segment_id = predictedSegmentId;
        else
            segment_id = read_segment_id(dec_handle, xd, mi_row, mi_col, r, 0);
        update_seg_ctx(&parse_ctxt->parse_nbr4x4_ctxt, mi_col, bw, bh, mbmi->seg_id_predicted);
    }
    else
        segment_id = read_segment_id(dec_handle, xd, mi_row, mi_col, r, 0);
    set_segment_id(dec_handle, mi_offset, x_mis, y_mis, segment_id);

    return segment_id;
}

void intra_block_mode_info(EbDecHandle *dec_handle, int mi_row,
    int mi_col, PartitionInfo_t* xd,
    ModeInfo_t *mbmi,
    SvtReader *r) {
    ParseCtxt *parse_ctxt = (ParseCtxt *)dec_handle->pv_parse_ctxt;
    const BlockSize bsize = mbmi->sb_type;
    mbmi->ref_frame[0] = INTRA_FRAME;
    mbmi->ref_frame[1] = NONE_FRAME;

    EbColorConfig *color_cfg = &dec_handle->seq_header.color_config;
    uint8_t     *lossless_array = &dec_handle->frame_header.lossless_array[0];

    mbmi->mode = read_intra_mode(r, parse_ctxt->frm_ctx[0].y_mode_cdf[size_group_lookup[bsize]]);

    mbmi->angle_delta[PLANE_TYPE_Y] =
        intra_angle_info(r, &parse_ctxt->frm_ctx[0].angle_delta_cdf[mbmi->mode - V_PRED][0], mbmi->mode, bsize);
    const int has_chroma =
        dec_is_chroma_reference(mi_row, mi_col, bsize, color_cfg->subsampling_x,
            color_cfg->subsampling_y);
    xd->has_chroma = has_chroma;
    if (has_chroma) {
        mbmi->uv_mode =
            read_intra_mode_uv(&parse_ctxt->frm_ctx[0], r, is_cfl_allowed(xd, color_cfg, lossless_array), mbmi->mode);
        if (mbmi->uv_mode == UV_CFL_PRED) {
            mbmi->cfl_alpha_idx =
                read_cfl_alphas(&parse_ctxt->frm_ctx[0], r, &mbmi->cfl_alpha_signs);
        }
        mbmi->angle_delta[PLANE_TYPE_UV] = intra_angle_info(r,
            &parse_ctxt->frm_ctx[0].angle_delta_cdf[mbmi->uv_mode - V_PRED][0], dec_get_uv_mode(mbmi->uv_mode), bsize);
    }

    if (allow_palette(dec_handle->seq_header.seq_force_screen_content_tools, bsize))
        palette_mode_info(/*dec_handle, mi_row, mi_col, r*/);

    filter_intra_mode_info(dec_handle, xd, r);
}

int read_is_inter(EbDecHandle* dec_handle, PartitionInfo_t * xd,
    int segment_id, SvtReader *r) {
    ParseCtxt *parse_ctxt = (ParseCtxt *)dec_handle->pv_parse_ctxt;
    int is_inter = 0;
    if (xd->mi->skip_mode)
        is_inter = 1;
    if (seg_feature_active(&dec_handle->frame_header.segmentation_params, segment_id, SEG_LVL_REF_FRAME))
        is_inter = get_segdata(&dec_handle->frame_header.segmentation_params, segment_id, SEG_LVL_REF_FRAME) != INTRA_FRAME;
    if (seg_feature_active(&dec_handle->frame_header.segmentation_params, segment_id, SEG_LVL_GLOBALMV))
        is_inter = 1;
    else {
        const int ctx = get_intra_inter_context(xd);
        is_inter = svt_read_symbol(r, parse_ctxt->frm_ctx[0].intra_inter_cdf[ctx], 2, ACCT_STR);
    }
    return is_inter;
}

void inter_frame_mode_info(EbDecHandle *dec_handle,
    PartitionInfo_t * xd, uint32_t mi_row,
    uint32_t mi_col, SvtReader *r, int8_t *cdef_strength) {
    ModeInfo_t *mbmi = xd->mi;
    int inter_block = 1;

    mbmi->mv[0].as_int = 0;
    mbmi->mv[1].as_int = 0;
    mbmi->segment_id = read_inter_segment_id(dec_handle, xd, mi_row, mi_col, 1, r);

    mbmi->skip_mode = read_skip_mode(dec_handle, xd, mbmi->segment_id, r);

    if (mbmi->skip_mode)
        mbmi->skip = 1;
    else
        mbmi->skip = read_skip(dec_handle, xd, mbmi->segment_id, r);

    if (!dec_handle->frame_header.segmentation_params.seg_id_pre_skip)
        mbmi->segment_id = read_inter_segment_id(dec_handle, xd, mi_row, mi_col, 0, r);

    read_cdef(dec_handle, r, xd, mi_col, mi_row, cdef_strength);

    /*TODO : Change to read_delta_params */
    assert(0);

    if (!mbmi->skip_mode)
        inter_block = read_is_inter(dec_handle, xd, mbmi->segment_id, r);

    if (inter_block);
        /*TO-DO fix for inter parse
        inter_block_mode_info(dec_handle, xd, mbmi, mi_row, mi_col, r);*/
    else
        intra_block_mode_info(dec_handle, mi_row, mi_col, xd, mbmi, r);
}

void mode_info(EbDecHandle *decHandle, PartitionInfo_t *part_info, uint32_t mi_row,
    uint32_t mi_col, SvtReader *r, int8_t *cdef_strength) {
    ModeInfo_t *mi = part_info->mi;
    mi->use_intrabc = 0;

    if (decHandle->frame_header.frame_type == KEY_FRAME || decHandle->frame_header.frame_type == INTRA_ONLY_FRAME)
        intra_frame_mode_info(decHandle, part_info, mi_row, mi_col, r, cdef_strength);
    else
        // TO-DO
        inter_frame_mode_info(decHandle, part_info, mi_row, mi_col, r, cdef_strength);
}

TxSize read_tx_size(EbDecHandle *dec_handle, PartitionInfo_t *xd,
                    int allow_select, SvtReader *r)
{
    ModeInfo_t *mbmi = xd->mi;
    const TxMode tx_mode = dec_handle->frame_header.tx_mode;
    const BlockSize bsize = xd->mi->sb_type;
    if (dec_handle->frame_header.lossless_array[mbmi->segment_id]) return TX_4X4;

    if (bsize > BLOCK_4X4 && allow_select && tx_mode == TX_MODE_SELECT) {
        const TxSize coded_tx_size = read_selected_tx_size(xd, r, dec_handle);
        return coded_tx_size;
    }
    assert(IMPLIES(tx_mode == ONLY_4X4, bsize == BLOCK_4X4));
    TxSize tx_size = max_txsize_rect_lookup[bsize];
    ParseCtxt *parse_ctxt = (ParseCtxt *)dec_handle->pv_parse_ctxt;
    update_tx_context(parse_ctxt, xd->mi_row, xd->mi_col, bsize, tx_size);
    return tx_size;
}

static INLINE void txfm_nbr_update(uint8_t *above_ctx,
    uint8_t *left_ctx,
    TxSize tx_size, TxSize txb_size) {
    BlockSize bsize = txsize_to_bsize[txb_size];
    int bh = mi_size_high[bsize];
    int bw = mi_size_wide[bsize];
    uint8_t txw = tx_size_wide[tx_size];
    uint8_t txh = tx_size_high[tx_size];
    memset(above_ctx, txw, bw);
    memset(left_ctx, txh, bh);
}

void read_var_tx_size(EbDecHandle *dec_handle, PartitionInfo_t *pi,
                               TxSize tx_size,
                               int blk_row, int blk_col) {
    ModeInfo_t *mbmi = pi->mi;
    const BlockSize bsize = mbmi->sb_type;
    (void)dec_handle;
    const int max_blocks_high = max_block_high(pi, bsize, 0);
    const int max_blocks_wide = max_block_wide(pi, bsize, 0);
    if (blk_row >= max_blocks_high || blk_col >= max_blocks_wide) return;
    assert(tx_size > TX_4X4);
    TxSize txs = max_txsize_rect_lookup[bsize];
    for (int level = 0; level < MAX_VARTX_DEPTH - 1; ++level)
        txs = sub_tx_size_map[txs];
    assert(0);
}

static INLINE TxSize av1_get_max_uv_txsize(BlockSize bsize, int subsampling_x,
    int subsampling_y) {
    const BlockSize plane_bsize =
        get_plane_block_size(bsize, subsampling_x, subsampling_y);
    assert(plane_bsize < BlockSizeS_ALL);
    const TxSize uv_tx = max_txsize_rect_lookup[plane_bsize];
    return av1_get_adjusted_tx_size(uv_tx);
}

/* Update Flat Transform Info for Intra Case! */
void update_flat_trans_info(EbDecHandle *dec_handle, PartitionInfo_t *part_info,
                            BlockSize bsize, TxSize tx_size)
{
    ParseCtxt *parse_ctx = (ParseCtxt*)dec_handle->pv_parse_ctxt;
    ModeInfo_t *mbmi = part_info->mi;
    SBInfo     *sb_info = part_info->sb_info;
    EbColorConfig color_config = dec_handle->seq_header.color_config;

    int num_luma_tus = 0;
    int num_chroma_tus = 0;

    int sx = color_config.subsampling_x;
    int sy = color_config.subsampling_x;

    TransformInfo_t *luma_trans_info = sb_info->sb_luma_trans_info +
                                    mbmi->first_luma_tu_offset;
    TransformInfo_t *chroma_trans_info = sb_info->sb_chroma_trans_info +
                                     mbmi->first_chroma_tu_offset;

    const TxSize max_tx_size = max_txsize_rect_lookup[bsize];
    const int bh = tx_size_high_unit[max_tx_size];
    const int bw = tx_size_wide_unit[max_tx_size];
    const int width = block_size_wide[bsize] >> tx_size_wide_log2[0];
    const int height = block_size_high[bsize] >> tx_size_high_log2[0];

    TxSize tx_size_uv = av1_get_max_uv_txsize(bsize, sx, sy);
    /* TODO: Make plane loop and avoid the unroll */
    for (int idy = 0; idy < height; idy += bh) {
        for (int idx = 0; idx < width; idx += bw) {
            /* Update Luma Transform Info */
            int stepr = tx_size_high_unit[tx_size];
            int stepc = tx_size_wide_unit[tx_size];

            /* TODO : Can cause prblm for incomplete SBs. Fix! */
            for (int blk_row = idy; blk_row < (idy + bh); blk_row += stepr) {
                for (int blk_col = idx; blk_col < (idx + bw); blk_col += stepc) {
                    luma_trans_info->tx_size = tx_size;
                    luma_trans_info++;
                    num_luma_tus++;
                }
            }

            if(color_config.mono_chrome)
                continue;

            /* Update Chroma Transform Info */
            if (!dec_is_chroma_reference(part_info->mi_row, part_info->mi_col,
                bsize, sx, sy))
                continue;

            stepr = tx_size_high_unit[tx_size_uv];
            stepc = tx_size_wide_unit[tx_size_uv];

            /* UV dim. of 4 special case! */
            const int unit_height = ROUND_POWER_OF_TWO(
                AOMMIN(bh + idy, bh), sy);
            const int unit_width = ROUND_POWER_OF_TWO(
                AOMMIN(bw + idx, bw), sx);
            /* TODO : Can cause prblm for incomplete SBs. Fix! */
            for (int blk_row = 0; blk_row < unit_height; blk_row += stepr) {
                for (int blk_col = 0; blk_col < unit_width; blk_col += stepc) {
                    chroma_trans_info->tx_size = tx_size_uv;
                    chroma_trans_info++;
                    num_chroma_tus++;
                }
            }
            for (int blk_row = 0; blk_row < unit_height; blk_row += stepr) {
                for (int blk_col = 0; blk_col < unit_width; blk_col += stepc) {
                    chroma_trans_info->tx_size = tx_size_uv;
                    chroma_trans_info++;
                }
            }
        }
    }

    mbmi->num_luma_tus = num_luma_tus;
    mbmi->num_chroma_tus = num_chroma_tus;

    parse_ctx->first_luma_tu_offset += num_luma_tus;
    parse_ctx->first_chroma_tu_offset += 2 * num_chroma_tus;
}

void read_block_tx_size(EbDecHandle *dec_handle, SvtReader *r, PartitionInfo_t *part_info,
                           BlockSize bsize) {
    ModeInfo_t *mbmi = part_info->mi;

    int inter_block_tx = dec_is_inter_block(mbmi);

    if (dec_handle->frame_header.tx_mode == TX_MODE_SELECT && bsize > BLOCK_4X4 &&
        !mbmi->skip && inter_block_tx &&
        !dec_handle->frame_header.lossless_array[mbmi->segment_id])
    {
        const TxSize max_tx_size = max_txsize_rect_lookup[bsize];
        const int bh = tx_size_high_unit[max_tx_size];
        const int bw = tx_size_wide_unit[max_tx_size];
        const int width = block_size_wide[bsize] >> tx_size_wide_log2[0];
        const int height = block_size_high[bsize] >> tx_size_high_log2[0];

        for (int idy = 0; idy < height; idy += bh)
            for (int idx = 0; idx < width; idx += bw)
                read_var_tx_size(dec_handle, part_info, max_tx_size, idx, idy);
    }
    else {
        TxSize tx_size = read_tx_size(dec_handle, part_info, !mbmi->skip || !inter_block_tx, r);
        /* Update Flat Transform Info */
        update_flat_trans_info(dec_handle, part_info, bsize, tx_size);
    }
}

int read_mv_component(SvtReader *r, NmvComponent *mvcomp,
    int use_subpel, int usehp) {
    int mag, d, fr, hp;
    const int sign = svt_read_symbol(r, mvcomp->sign_cdf, 2, ACCT_STR);
    const int mv_class =
        svt_read_symbol(r, mvcomp->classes_cdf, MV_CLASSES, ACCT_STR);
    const int class0 = mv_class == MV_CLASS_0;

    // Integer part
    if (class0) {
        d = svt_read_symbol(r, mvcomp->class0_cdf, CLASS0_SIZE, ACCT_STR);
        mag = 0;
    }
    else {
        const int n = mv_class + CLASS0_BITS - 1;  // number of bits
        d = 0;
        for (int i = 0; i < n; ++i)
            d |= svt_read_symbol(r, mvcomp->bits_cdf[i], 2, ACCT_STR) << i;
        mag = CLASS0_SIZE << (mv_class + 2);
    }

    fr = use_subpel ? svt_read_symbol(r, class0 ? mvcomp->class0_fp_cdf[d] :
        mvcomp->fp_cdf, MV_FP_SIZE, ACCT_STR) : 3;

    hp = usehp ? svt_read_symbol(r, class0 ? mvcomp->class0_hp_cdf :
        mvcomp->hp_cdf, 2, ACCT_STR) : 1;

    // Result
    mag += ((d << 3) | (fr << 1) | hp) + 1;
    return sign ? -mag : mag;
}

void read_mv(SvtReader *r, MV *mv, MV *ref,
    NmvContext *ctx, MvSubpelPrecision precision) {
    MV diff = kZeroMv;

    const MvJointType joint_type =
        (MvJointType)svt_read_symbol(r, ctx->joints_cdf, MV_JOINTS, ACCT_STR);

    if (mv_joint_vertical(joint_type))
        diff.row = read_mv_component(r, &ctx->comps[0], precision > MV_SUBPEL_NONE,
            precision > MV_SUBPEL_LOW_PRECISION);

    if (mv_joint_horizontal(joint_type))
        diff.col = read_mv_component(r, &ctx->comps[1], precision > MV_SUBPEL_NONE,
            precision > MV_SUBPEL_LOW_PRECISION);

    mv->row = ref->row + diff.row;
    mv->col = ref->col + diff.col;
}

BlockSize get_plane_residual_size(BlockSize bsize,
    int subsampling_x,
    int subsampling_y) {
    if (bsize == BLOCK_INVALID) return BLOCK_INVALID;
    return ss_size_lookup[bsize][subsampling_x][subsampling_y];
}

TxSetType get_tx_set(TxSize tx_size, int is_inter,
    int use_reduced_set) {
    return get_ext_tx_set_type(tx_size, is_inter, use_reduced_set);
}

void parse_transform_type(EbDecHandle *dec_handle, PartitionInfo_t *xd,
     TxSize tx_size, SvtReader *r, TransformInfo_t *trans_info) {
    ModeInfo_t *mbmi = xd->mi;

    TxType *tx_type = &trans_info->txk_type;
    *tx_type = DCT_DCT;

    ParseCtxt *parse_ctxt = (ParseCtxt *)dec_handle->pv_parse_ctxt;
    FRAME_CONTEXT *frm_ctx = &parse_ctxt->frm_ctx[0];

    // No need to read transform type if block is skipped.
    if (mbmi->skip || seg_feature_active(&dec_handle->frame_header.segmentation_params, mbmi->segment_id, SEG_LVL_SKIP))
        return;

    const int qindex = dec_handle->frame_header.quantization_params.base_q_idx;
    if (qindex == 0) return;

    const int inter_block = dec_is_inter_block(mbmi);
    const TxSetType tx_set_type =
        get_ext_tx_set_type(tx_size, inter_block, dec_handle->frame_header.reduced_tx_set);
    if (tx_set_type > 1) {
        const int eset =
            get_ext_tx_set(tx_size, inter_block, dec_handle->frame_header.reduced_tx_set);
        // eset == 0 should correspond to a set with only DCT_DCT and
        // there is no need to read the tx_type
        assert(eset != 0);

        const TxSize square_tx_size = txsize_sqr_map[tx_size];
        if (inter_block) {
            *tx_type = av1_ext_tx_inv[tx_set_type][svt_read_symbol(
                r, frm_ctx->inter_ext_tx_cdf[eset][square_tx_size],
                av1_num_ext_tx_set[tx_set_type], ACCT_STR)];
        }
        else {
            const PredictionMode intra_mode =
                mbmi->filter_intra_mode_info.use_filter_intra
                ? fimode_to_intradir[mbmi->filter_intra_mode_info
                .filter_intra_mode]
                : mbmi->mode;
            *tx_type = av1_ext_tx_inv[tx_set_type][svt_read_symbol(
                r, frm_ctx->intra_ext_tx_cdf[eset][square_tx_size][intra_mode],
                av1_num_ext_tx_set[tx_set_type], ACCT_STR)];
        }
    }
}

void clamp_mv_row_col()
{
    //TO-DO

    /*clamp_mv(mv, xd->mb_to_left_edge - bw * 8 - MV_BORDER,
        xd->mb_to_right_edge + bw * 8 + MV_BORDER,
        xd->mb_to_top_edge - bh * 8 - MV_BORDER,
        xd->mb_to_bottom_edge + bh * 8 + MV_BORDER);*/
}

const ScanOrder* get_scan(TxSize tx_size, TxType tx_type) {
    return &av1_scan_orders[tx_size][tx_type];
}

TxType compute_tx_type(PlaneType plane_type,
    const PartitionInfo_t *xd, uint32_t blk_row,
    uint32_t blk_col, TxSize tx_size,
    int reduced_tx_set, EbColorConfig *color_cfg,
    uint8_t *lossless_array, TransformInfo_t *trans_info)
{
    const ModeInfo_t *const mbmi = xd->mi;
    const TxSetType tx_set_type =
        get_ext_tx_set_type(tx_size, dec_is_inter_block(mbmi), reduced_tx_set);

    TxType tx_type;
    if (lossless_array[mbmi->segment_id] || txsize_sqr_up_map[tx_size] > TX_32X32)
        tx_type = DCT_DCT;
    else {
        if (plane_type == PLANE_TYPE_Y)
            tx_type = trans_info->txk_type;
        else if (dec_is_inter_block(mbmi)) {
            // scale back to y plane's coordinate
            blk_row <<= color_cfg->subsampling_y;
            blk_col <<= color_cfg->subsampling_x;

            assert(0);
        }
        else {
            // In intra mode, uv planes don't share the same prediction mode as y
            // plane, so the tx_type should not be shared
            tx_type = intra_mode_to_tx_type(mbmi, PLANE_TYPE_UV);
        }
    }
    assert(tx_type < TX_TYPES);
    if (!av1_ext_tx_used[tx_set_type][tx_type]) return DCT_DCT;
    return tx_type;
}

static INLINE int is_interintra_wedge_used(BlockSize sb_type) {
    return wedge_params_lookup[sb_type].bits > 0;
}

static INLINE int is_comp_ref_allowed(BlockSize bsize) {
    return AOMMIN(block_size_wide[bsize], block_size_high[bsize]) >= 8;
}

static INLINE int is_masked_compound_type(CompoundType type) {
    return (type == COMPOUND_WEDGE || type == COMPOUND_DIFFWTD);
}

static INLINE int is_interinter_compound_used(CompoundType type,
    BlockSize sb_type) {
    const int comp_allowed = is_comp_ref_allowed(sb_type);
    switch (type) {
    case COMPOUND_AVERAGE:
    case COMPOUND_DISTWTD:
    case COMPOUND_DIFFWTD: return comp_allowed;
    case COMPOUND_WEDGE:
        return comp_allowed && wedge_params_lookup[sb_type].bits > 0;
    default: assert(0); return 0;
    }
}

static INLINE int is_any_masked_compound_used(BlockSize sb_type) {
    CompoundType comp_type;
    int i;
    if (!is_comp_ref_allowed(sb_type)) return 0;
    for (i = 0; i < COMPOUND_TYPES; i++) {
        comp_type = (CompoundType)i;
        if (is_masked_compound_type(comp_type) &&
            is_interinter_compound_used(comp_type, sb_type))
            return 1;
    }
    return 0;
}

void read_interintra_mode(EbDecHandle *dec_handle, ModeInfo_t *mbmi,
    SvtReader *r, int size_group) {
    ParseCtxt *parse_ctxt = (ParseCtxt *)dec_handle->pv_parse_ctxt;
    FRAME_CONTEXT *frm_ctx = &parse_ctxt->frm_ctx[0];
    BlockSize bsize = mbmi->sb_type;
    if (dec_handle->seq_header.enable_interintra_compound && !mbmi->skip_mode &&
        is_interintra_allowed(mbmi)) {
        const int bsize_group = size_group_lookup[bsize];
        const int interintra =
            svt_read_symbol(r, frm_ctx->interintra_cdf[bsize_group], 2, ACCT_STR);
        assert(mbmi->ref_frame[1] == NONE_FRAME);
        if (interintra) {
            mbmi->is_inter_intra =
                (InterIntraMode)svt_read_symbol(
                    r, frm_ctx->interintra_mode_cdf[size_group], INTERINTRA_MODES,
                    ACCT_STR);
            mbmi->ref_frame[1] = INTRA_FRAME;
            mbmi->interintra_mode.interintra_mode = mbmi->is_inter_intra;
            mbmi->angle_delta[PLANE_TYPE_Y] = 0;
            mbmi->angle_delta[PLANE_TYPE_UV] = 0;
            mbmi->filter_intra_mode_info.use_filter_intra = 0;
            if (is_interintra_wedge_used(bsize)) {
                mbmi->interintra_mode.wedge_interintra = svt_read_symbol(
                    r, frm_ctx->wedge_interintra_cdf[bsize], 2, ACCT_STR);
                if (mbmi->interintra_mode.wedge_interintra) {
                    mbmi->interintra_mode.interintra_wedge_index =
                        svt_read_symbol(r, frm_ctx->wedge_idx_cdf[bsize], 16, ACCT_STR);
                }
            }
        }
    }
}

static INLINE int get_comp_group_idx_context(const PartitionInfo_t *xd) {
    const ModeInfo_t *const above_mi = xd->above_mbmi;
    const ModeInfo_t *const left_mi = xd->left_mbmi;
    int above_ctx = 0, left_ctx = 0;

    if (above_mi) {
        if (has_second_ref(above_mi))
            above_ctx = above_mi->inter_compound.comp_group_idx;
        else if (above_mi->ref_frame[0] == ALTREF_FRAME)
            above_ctx = 3;
    }
    if (left_mi) {
        if (has_second_ref(left_mi))
            left_ctx = left_mi->inter_compound.comp_group_idx;
        else if (left_mi->ref_frame[0] == ALTREF_FRAME)
            left_ctx = 3;
    }

    return AOMMIN(5, above_ctx + left_ctx);
}

void read_compound_type(EbDecHandle *dec_handle, PartitionInfo_t *xd, ModeInfo_t *mbmi,
    SvtReader *r) {
    BlockSize bsize = mbmi->sb_type;
    mbmi->inter_compound.comp_group_idx = 0;
    mbmi->inter_compound.compound_idx = 1;
    ParseCtxt *parse_ctxt = (ParseCtxt *)dec_handle->pv_parse_ctxt;
    FRAME_CONTEXT *frm_ctx = &parse_ctxt->frm_ctx[0];

    if(mbmi->skip_mode) mbmi->inter_compound.type = COMPOUND_AVERAGE;

    if (has_second_ref(mbmi)) {
        // Read idx to indicate current compound inter prediction mode group
        const int masked_compound_used = is_any_masked_compound_used(bsize) &&
            dec_handle->seq_header.enable_masked_compound;

        if (masked_compound_used) {
            const int ctx_comp_group_idx = get_comp_group_idx_context(xd);
            mbmi->inter_compound.comp_group_idx = svt_read_symbol(
                r, frm_ctx->comp_group_idx_cdf[ctx_comp_group_idx], 2, ACCT_STR);
        }

        if (mbmi->inter_compound.comp_group_idx == 0) {
            if (dec_handle->seq_header.order_hint_info.enable_jnt_comp) {
                const int comp_index_ctx = get_comp_index_context(dec_handle/*, xd*/);
                mbmi->inter_compound.compound_idx = svt_read_symbol(
                    r, frm_ctx->compound_index_cdf[comp_index_ctx], 2, ACCT_STR);
                mbmi->inter_compound.type =
                    mbmi->inter_compound.compound_idx ? COMPOUND_AVERAGE : COMPOUND_DISTWTD;
            }
            else {
                // Distance-weighted compound is disabled, so always use average
                mbmi->inter_compound.compound_idx = 1;
                mbmi->inter_compound.type = COMPOUND_AVERAGE;
            }
        }
        else {
            assert(dec_handle->frame_header.reference_mode != SINGLE_REFERENCE &&
                is_inter_compound_mode(mbmi->mode) &&
                mbmi->motion_mode == SIMPLE_TRANSLATION);
            assert(masked_compound_used);

            // compound_diffwtd, wedge
            if (is_interinter_compound_used(COMPOUND_WEDGE, bsize))
                mbmi->inter_compound.type =
                COMPOUND_WEDGE + svt_read_symbol(r,
                    frm_ctx->compound_type_cdf[bsize],
                    MASKED_COMPOUND_TYPES, ACCT_STR);
            else
                mbmi->inter_compound.type = COMPOUND_DIFFWTD;

            if (mbmi->inter_compound.type == COMPOUND_WEDGE) {
                assert(is_interinter_compound_used(COMPOUND_WEDGE, bsize));
                mbmi->inter_compound.wedge_index =
                    svt_read_symbol(r, frm_ctx->wedge_idx_cdf[bsize], 16, ACCT_STR);
                mbmi->inter_compound.wedge_sign = svt_read_bit(r, ACCT_STR);
            }
            else {
                assert(mbmi->inter_compound.type == COMPOUND_DIFFWTD);
                mbmi->inter_compound.mask_type =
                    svt_read_literal(r, MAX_DIFFWTD_MASK_BITS, ACCT_STR);
            }
        }
    }
    else {
        if (mbmi->is_inter_intra)
            mbmi->inter_compound.type = mbmi->interintra_mode.wedge_interintra ?
            COMPOUND_WEDGE : COMPOUND_INTRA;
        else
            mbmi->inter_compound.type = COMPOUND_AVERAGE;
    }
}

MotionMode read_motion_mode(EbDecHandle *dec_handle,
    ModeInfo_t *mbmi, SvtReader *r) {
    ParseCtxt *parse_ctxt = (ParseCtxt *)dec_handle->pv_parse_ctxt;
    FRAME_CONTEXT *frm_ctx = &parse_ctxt->frm_ctx[0];

    //if (dec_handle->switchable_motion_mode == 0) return SIMPLE_TRANSLATION;
    if (mbmi->skip_mode) return SIMPLE_TRANSLATION;

    MotionMode last_motion_mode_allowed = SIMPLE_TRANSLATION;
    int motion_mode;

    if (last_motion_mode_allowed == SIMPLE_TRANSLATION) return SIMPLE_TRANSLATION;

    if (last_motion_mode_allowed == OBMC_CAUSAL) {
        motion_mode =
            svt_read_symbol(r, frm_ctx->obmc_cdf[mbmi->sb_type], 2, ACCT_STR);
        return (MotionMode)(SIMPLE_TRANSLATION + motion_mode);
    }
    else {
        motion_mode =
            svt_read_symbol(r, frm_ctx->motion_mode_cdf[mbmi->sb_type],
                MOTION_MODES, ACCT_STR);
        return (MotionMode)(SIMPLE_TRANSLATION + motion_mode);
    }
}

PredictionMode getMode(PredictionMode yMode, int refList)
{
    PredictionMode compMode = GLOBALMV;
    if (yMode == NEW_NEWMV) compMode = NEWMV;
    if (refList == 0)
    {
        if (yMode < NEAREST_NEARESTMV)
            compMode = yMode;
        else if (yMode == NEW_NEARESTMV || yMode == NEW_NEARMV)
            compMode = NEWMV;
        else if (yMode == NEAREST_NEARESTMV || yMode == NEAREST_NEWMV)
            compMode = NEARESTMV;
        else if (yMode == NEAR_NEARMV || yMode == NEAR_NEWMV)
            compMode = NEARMV;
    }
    else
    {
        if (yMode == NEAREST_NEWMV || yMode == NEAR_NEWMV)
            compMode = NEWMV;
        else if (yMode == NEAREST_NEARESTMV || yMode == NEW_NEARESTMV)
            compMode = NEARESTMV;
        else if (yMode == NEAR_NEARMV || yMode == NEW_NEARMV)
            compMode = NEARMV;
    }
    return compMode;
}

static INLINE int is_mv_valid(const MV *mv) {
    return mv->row > MV_LOW && mv->row < MV_UPP && mv->col > MV_LOW &&
        mv->col < MV_UPP;
}

int assign_mv(EbDecHandle *dec_handle, PartitionInfo_t *xd,
    MvReferenceFrame ref_frame[2], IntMv mv[2],
    IntMv ref_mv[2], IntMv nearest_mv[2],
    IntMv near_mv[2], int mi_row, int mi_col,
    int is_compound, int allow_hp, SvtReader *r)
{
    ModeInfo_t *mbmi = xd->mi;
    BlockSize bsize = mbmi->sb_type;
    ParseCtxt *parse_ctxt = (ParseCtxt *)dec_handle->pv_parse_ctxt;
    FRAME_CONTEXT *frm_ctx = &parse_ctxt->frm_ctx[0];
    if (dec_handle->frame_header.force_integer_mv)
        allow_hp = MV_SUBPEL_NONE;
    PredictionMode mode;
    for (int i = 0; i < is_compound + 1; i++)
    {
        if (mbmi->use_intrabc)
            mode = NEWMV;
        else
            mode = getMode(mbmi->mode, i);

        switch (mode) {
        case NEWMV: {
            NmvContext *const nmvc = &frm_ctx->nmvc;
            read_mv(r, &mv[0].as_mv, &ref_mv[0].as_mv, nmvc, allow_hp);
            break;
        }
        case NEARESTMV: {
            mv[0].as_int = nearest_mv[0].as_int;
            break;
        }
        case NEARMV: {
            mv[0].as_int = near_mv[0].as_int;
            break;
        }
        case GLOBALMV: {
            mv[0].as_int =
                gm_get_motion_vector(&xd->global_motion[ref_frame[0]],
                    dec_handle->frame_header.allow_high_precision_mv, bsize, mi_col,
                    mi_row, dec_handle->frame_header.force_integer_mv)
                .as_int;
            break;
        }
        case NEW_NEWMV: {
            assert(is_compound);
            for (int i = 0; i < 2; ++i) {
                NmvContext *const nmvc = &frm_ctx->nmvc;
                read_mv(r, &mv[i].as_mv, &ref_mv[i].as_mv, nmvc, allow_hp);
            }
            break;
        }
        case NEAREST_NEARESTMV: {
            assert(is_compound);
            mv[0].as_int = nearest_mv[0].as_int;
            mv[1].as_int = nearest_mv[1].as_int;
            break;
        }
        case NEAR_NEARMV: {
            assert(is_compound);
            mv[0].as_int = near_mv[0].as_int;
            mv[1].as_int = near_mv[1].as_int;
            break;
        }
        case NEW_NEARESTMV: {
            NmvContext *const nmvc = &frm_ctx->nmvc;
            read_mv(r, &mv[0].as_mv, &ref_mv[0].as_mv, nmvc, allow_hp);
            assert(is_compound);
            mv[1].as_int = nearest_mv[1].as_int;
            break;
        }
        case NEAREST_NEWMV: {
            NmvContext *const nmvc = &frm_ctx->nmvc;
            mv[0].as_int = nearest_mv[0].as_int;
            read_mv(r, &mv[1].as_mv, &ref_mv[1].as_mv, nmvc, allow_hp);
            assert(is_compound);
            break;
        }
        case NEAR_NEWMV: {
            NmvContext *const nmvc = &frm_ctx->nmvc;
            mv[0].as_int = near_mv[0].as_int;
            read_mv(r, &mv[1].as_mv, &ref_mv[1].as_mv, nmvc, allow_hp);
            assert(is_compound);
            break;
        }
        case NEW_NEARMV: {
            NmvContext *const nmvc = &frm_ctx->nmvc;
            read_mv(r, &mv[0].as_mv, &ref_mv[0].as_mv, nmvc, allow_hp);
            assert(is_compound);
            mv[1].as_int = near_mv[1].as_int;
            break;
        }
        case GLOBAL_GLOBALMV: {
            assert(is_compound);
            mv[0].as_int =
                gm_get_motion_vector(&xd->global_motion[ref_frame[0]],
                    dec_handle->frame_header.allow_high_precision_mv, bsize, mi_col,
                    mi_row, dec_handle->frame_header.force_integer_mv)
                .as_int;
            mv[1].as_int =
                gm_get_motion_vector(&xd->global_motion[ref_frame[1]],
                    dec_handle->frame_header.allow_high_precision_mv, bsize, mi_col,
                    mi_row, dec_handle->frame_header.force_integer_mv)
                .as_int;
            break;
        }
        default: { return 0; }
        }
    }
        int ret = is_mv_valid(&mv[0].as_mv);
        if (is_compound)
            ret = ret && is_mv_valid(&mv[1].as_mv);
    return ret;
}

static AOM_FORCE_INLINE int get_br_ctx_eob(const int c,  // raster order
    const int bwl,
    const TxClass tx_class) {
    const int row = c >> bwl;
    const int col = c - (row << bwl);
    if (c == 0) return 0;
    if ((tx_class == TX_CLASS_2D && row < 2 && col < 2) ||
        (tx_class == TX_CLASS_HORIZ && col == 0) ||
        (tx_class == TX_CLASS_VERT && row == 0))
        return 7;
    return 14;
}

void update_coeff_ctx(EbDecHandle *dec_handle, int plane, PartitionInfo_t *pi, TxSize tx_size,
    uint32_t blk_row, uint32_t blk_col, int above_off, int left_off, int cul_level, int dc_val)
{
    ParseCtxt *parse_ctxt = (ParseCtxt *)dec_handle->pv_parse_ctxt;
    ParseNbr4x4Ctxt *ngr_ctx = &parse_ctxt->parse_nbr4x4_ctxt;

    uint8_t suby = plane ? dec_handle->seq_header.color_config.subsampling_y : 0;
    uint8_t subx = plane ? dec_handle->seq_header.color_config.subsampling_x : 0;

    uint8_t *const above_dc_ctx = ngr_ctx->above_dc_ctx[plane] + blk_col;
    uint8_t *const left_dc_ctx = ngr_ctx->left_dc_ctx[plane] + (blk_row - (parse_ctxt->sb_row_mi >> suby));

    uint8_t *const above_level_ctx = ngr_ctx->above_level_ctx[plane] + blk_col;
    uint8_t *const left_level_ctx = ngr_ctx->left_level_ctx[plane] + (blk_row - (parse_ctxt->sb_row_mi >> suby));

    const int txs_wide = tx_size_wide_unit[tx_size];
    const int txs_high = tx_size_high_unit[tx_size];

    if (pi->mb_to_right_edge < 0) {
        int plane_bsize = (pi->mi->sb_type == BLOCK_INVALID) ? BLOCK_INVALID :
            ss_size_lookup[pi->mi->sb_type][subx][suby];
        const int blocks_wide = max_block_wide(pi, plane_bsize, subx);
        const int above_contexts = AOMMIN(txs_wide, (blocks_wide - above_off));

        memset(above_dc_ctx, dc_val, above_contexts);
        memset(above_dc_ctx + above_contexts, 0, (txs_wide - above_contexts));
        memset(above_level_ctx, cul_level, above_contexts);
        memset(above_level_ctx + above_contexts, 0, (txs_wide - above_contexts));
    }
    else {
        memset(above_dc_ctx, dc_val, txs_wide);
        memset(above_level_ctx, cul_level, txs_wide);
    }

    if (pi->mb_to_bottom_edge < 0) {
        int plane_bsize = (pi->mi->sb_type == BLOCK_INVALID) ? BLOCK_INVALID :
            ss_size_lookup[pi->mi->sb_type][subx][suby];
        const int blocks_high = max_block_high(pi, plane_bsize, suby);
        const int left_contexts = AOMMIN(txs_high, (blocks_high - left_off));

        memset(left_dc_ctx, dc_val, left_contexts);
        memset(left_dc_ctx + left_contexts, 0, txs_high - left_contexts);
        memset(left_level_ctx, cul_level, left_contexts);
        memset(left_level_ctx + left_contexts, 0, txs_high - left_contexts);
    }
    else {
        memset(left_dc_ctx, dc_val, txs_high);
        memset(left_level_ctx, cul_level, txs_high);
    }
}

static INLINE int rec_eob_pos(const int eob_token, const int extra) {
    int eob = k_eob_group_start[eob_token];
    if (eob > 2)
        eob += extra;
    return eob;
}

static INLINE int get_lower_levels_ctx_2d(const uint8_t *levels, int coeff_idx,
    int bwl, TxSize tx_size) {
    assert(coeff_idx > 0);
    int mag;
    levels = levels + get_padded_idx(coeff_idx, bwl);
    mag = AOMMIN(levels[1], 3);                                     // { 0, 1 }
    mag += AOMMIN(levels[(1 << bwl) + TX_PAD_HOR], 3);              // { 1, 0 }
    mag += AOMMIN(levels[(1 << bwl) + TX_PAD_HOR + 1], 3);          // { 1, 1 }
    mag += AOMMIN(levels[2], 3);                                    // { 0, 2 }
    mag += AOMMIN(levels[(2 << bwl) + (2 << TX_PAD_HOR_LOG2)], 3);  // { 2, 0 }

    const int ctx = AOMMIN((mag + 1) >> 1, 4);
    return ctx + av1_nz_map_ctx_offset[tx_size][coeff_idx];
}
static INLINE int get_lower_levels_ctx(const uint8_t *levels,
    int coeff_idx, int bwl,
    TxSize tx_size,
    TxClass tx_class) {
    const int stats =
        get_nz_mag(levels + get_padded_idx(coeff_idx, bwl), bwl, tx_class);
    return get_nz_map_ctx_from_stats(stats, coeff_idx, bwl, tx_size, tx_class);
}

static int read_golomb(SvtReader *r) {
    int x = 1;
    int length = 0;
    int i = 0;

    while (!i) {
        i = svt_read_bit(r, ACCT_STR);
        ++length;
        if (length > 20) {
            printf("Invalid length in read_golomb");
            break;
        }
    }

    for (i = 0; i < length - 1; ++i) {
        x <<= 1;
        x += svt_read_bit(r, ACCT_STR);
    }

    return x - 1;
}

static INLINE int get_br_ctx(const uint8_t *const levels,
    const int c,  // raster order
    const int bwl, const TxClass tx_class) {
    const int row = c >> bwl;
    const int col = c - (row << bwl);
    const int stride = (1 << bwl) + TX_PAD_HOR;
    const int pos = row * stride + col;
    int mag = levels[pos + 1];
    mag += levels[pos + stride];
    switch (tx_class) {
    case TX_CLASS_2D:
        mag += levels[pos + stride + 1];
        mag = AOMMIN((mag + 1) >> 1, 6);
        if (c == 0) return mag;
        if ((row < 2) && (col < 2)) return mag + 7;
        break;
    case TX_CLASS_HORIZ:
        mag += levels[pos + 2];
        mag = AOMMIN((mag + 1) >> 1, 6);
        if (c == 0) return mag;
        if (col == 0) return mag + 7;
        break;
    case TX_CLASS_VERT:
        mag += levels[pos + (stride << 1)];
        mag = AOMMIN((mag + 1) >> 1, 6);
        if (c == 0) return mag;
        if (row == 0) return mag + 7;
        break;
    default: break;
    }

    return mag + 14;
}

static INLINE int get_br_ctx_2d(const uint8_t *const levels,
    const int c,  // raster order
    const int bwl) {
    assert(c > 0);
    const int row = c >> bwl;
    const int col = c - (row << bwl);
    const int stride = (1 << bwl) + TX_PAD_HOR;
    const int pos = row * stride + col;
    int mag = AOMMIN(levels[pos + 1], MAX_BASE_BR_RANGE) +
        AOMMIN(levels[pos + stride], MAX_BASE_BR_RANGE) +
        AOMMIN(levels[pos + 1 + stride], MAX_BASE_BR_RANGE);
    mag = AOMMIN((mag + 1) >> 1, 6);
    if ((row | col) < 2) return mag + 7;
    return mag + 14;
}

static INLINE void read_coeffs_reverse_2d(SvtReader *r, TxSize tx_size,
    int start_si, int end_si,
    const int16_t *scan, int bwl,
    uint8_t *levels,
    base_cdf_arr base_cdf,
    br_cdf_arr br_cdf) {
    for (int c = end_si; c >= start_si; --c) {
        const int pos = scan[c];
        const int coeff_ctx = get_lower_levels_ctx_2d(levels, pos, bwl, tx_size);
        const int nsymbs = 4;
        int level = svt_read_symbol(r, base_cdf[coeff_ctx], nsymbs, ACCT_STR);
        if (level > NUM_BASE_LEVELS) {
            const int br_ctx = get_br_ctx_2d(levels, pos, bwl);
            AomCdfProb *cdf = br_cdf[br_ctx];
            for (int idx = 0; idx < COEFF_BASE_RANGE; idx += BR_CDF_SIZE - 1) {
                const int k = svt_read_symbol(r, cdf, BR_CDF_SIZE, ACCT_STR);
                level += k;
                if (k < BR_CDF_SIZE - 1) break;
            }
        }
        levels[get_padded_idx(pos, bwl)] = level;
    }
}

static INLINE void read_coeffs_reverse(SvtReader *r, TxSize tx_size,
    TxClass tx_class, int start_si,
    int end_si, const int16_t *scan, int bwl,
    uint8_t *levels, base_cdf_arr base_cdf,
    br_cdf_arr br_cdf) {
    for (int c = end_si; c >= start_si; --c) {
        const int pos = scan[c];
        const int coeff_ctx =
            get_lower_levels_ctx(levels, pos, bwl, tx_size, tx_class);
        const int nsymbs = 4;
        int level = svt_read_symbol(r, base_cdf[coeff_ctx], nsymbs, ACCT_STR);
        if (level > NUM_BASE_LEVELS) {
            const int br_ctx = get_br_ctx(levels, pos, bwl, tx_class);
            AomCdfProb *cdf = br_cdf[br_ctx];
            for (int idx = 0; idx < COEFF_BASE_RANGE; idx += BR_CDF_SIZE - 1) {
                const int k = svt_read_symbol(r, cdf, BR_CDF_SIZE, ACCT_STR);
                level += k;
                if (k < BR_CDF_SIZE - 1) break;
            }
        }
        levels[get_padded_idx(pos, bwl)] = level;
    }
}

uint16_t parse_coeffs(EbDecHandle *dec_handle, PartitionInfo_t *xd, SvtReader *r,
    uint32_t blk_row, uint32_t blk_col, int above_off, int left_off, int plane,
    int txb_skip_ctx, int dc_sign_ctx,
    TxSize tx_size, int32_t *coeff_buf, TransformInfo_t *trans_info)
{
    const int width = get_txb_wide(tx_size);
    const int height = get_txb_high(tx_size);

    ParseCtxt *parse_ctxt = (ParseCtxt *)dec_handle->pv_parse_ctxt;
    FRAME_CONTEXT *frm_ctx = &parse_ctxt->frm_ctx[0];

    TxSize txs_ctx = (TxSize)((txsize_sqr_map[tx_size] + txsize_sqr_up_map[tx_size] + 1) >> 1);
    PlaneType plane_type = (plane == 0) ? PLANE_TYPE_Y : PLANE_TYPE_UV;
    int cul_level = 0;
    int dc_val = 0;
    uint8_t levels_buf[TX_PAD_2D];
    uint8_t *const levels = set_levels(levels_buf, width);
    EbColorConfig *color_cfg = &dec_handle->seq_header.color_config;

    const int all_zero = svt_read_symbol(
        r, frm_ctx->txb_skip_cdf[txs_ctx][txb_skip_ctx], 2, ACCT_STR);

    const int bwl = get_txb_bwl(tx_size);

    uint16_t eob = 0;
    uint16_t max_scan_line = 0;
    if (all_zero) {
        if (plane == 0) {
            trans_info->txk_type = DCT_DCT;
            trans_info->cbf      = 0;
        }

        update_coeff_ctx(dec_handle, plane, xd, tx_size, blk_row, blk_col,
            above_off, left_off, cul_level, dc_val);

        return 0;
    }

    if (plane == AOM_PLANE_Y)
        parse_transform_type(dec_handle, xd, tx_size, r, trans_info);

    uint8_t     *lossless_array = &dec_handle->frame_header.lossless_array[0];
    trans_info->txk_type = compute_tx_type(plane_type, xd, blk_row, blk_col,
        tx_size, dec_handle->frame_header.reduced_tx_set, color_cfg, lossless_array, trans_info);
    const ScanOrder *scan_order = get_scan(tx_size, trans_info->txk_type);
    const int16_t *const scan = scan_order->scan;
    const int eob_multi_size = txsize_log2_minus4[tx_size];

    const TxClass tx_class = tx_type_to_class[trans_info->txk_type];
    const int eob_multi_ctx = (tx_class == TX_CLASS_2D) ? 0 : 1;

    int eob_extra = 0;
    int eob_pt = 1;

    switch (eob_multi_size) {
    case 0:
        eob_pt =
            svt_read_symbol(r, frm_ctx->eob_flag_cdf16[plane_type][eob_multi_ctx],
                5, ACCT_STR) + 1;
        break;
    case 1:
        eob_pt =
            svt_read_symbol(r, frm_ctx->eob_flag_cdf32[plane_type][eob_multi_ctx],
                6, ACCT_STR) + 1;
        break;
    case 2:
        eob_pt =
            svt_read_symbol(r, frm_ctx->eob_flag_cdf64[plane_type][eob_multi_ctx],
                7, ACCT_STR) + 1;
        break;
    case 3:
        eob_pt =
            svt_read_symbol(r, frm_ctx->eob_flag_cdf128[plane_type][eob_multi_ctx],
                8, ACCT_STR) + 1;
        break;
    case 4:
        eob_pt =
            svt_read_symbol(r, frm_ctx->eob_flag_cdf256[plane_type][eob_multi_ctx],
                9, ACCT_STR) + 1;
        break;
    case 5:
        eob_pt =
            svt_read_symbol(r, frm_ctx->eob_flag_cdf512[plane_type][eob_multi_ctx],
                10, ACCT_STR) + 1;
        break;
    default:
        eob_pt = svt_read_symbol(
            r, frm_ctx->eob_flag_cdf1024[plane_type][eob_multi_ctx], 11,
            ACCT_STR) + 1;
        break;
    }

    int eob_shift = k_eob_offset_bits[eob_pt];
    if (eob_shift > 0) {
        const int eob_ctx = eob_pt;
        int bit = svt_read_symbol(
            r, frm_ctx->eob_extra_cdf[txs_ctx][plane_type][eob_ctx], 2, ACCT_STR);
        if (bit)
            eob_extra += (1 << (eob_shift - 1));
        for (int i = 1; i < eob_shift; i++) {
            bit = svt_read_bit(r, ACCT_STR);
            if (bit)
                eob_extra += (1 << (eob_shift - 1 - i));
        }
    }
    eob = rec_eob_pos(eob_pt, eob_extra);

    if (eob > 1) {
        memset(levels_buf, 0,
            sizeof(*levels_buf) *
            ((width + TX_PAD_HOR) * (height + TX_PAD_VER) + TX_PAD_END));
    }

    {
        int i = eob - 1;
        const int pos = scan[i];
        const int coeff_ctx = get_lower_levels_ctx_eob(bwl, height, i);
        const int nsymbs = 3;
        AomCdfProb *cdf =
            frm_ctx->coeff_base_eob_cdf[txs_ctx][plane_type][coeff_ctx];
        int level = svt_read_symbol(r, cdf, nsymbs, ACCT_STR) + 1;
        if (level > NUM_BASE_LEVELS) {
            const int br_ctx = get_br_ctx_eob(pos, bwl, tx_class);
            cdf = frm_ctx->coeff_br_cdf[AOMMIN(txs_ctx, TX_32X32)][plane_type][br_ctx];
            for (int idx = 0; idx < COEFF_BASE_RANGE/ (BR_CDF_SIZE - 1); idx ++) {
                int coeff_br = svt_read_symbol(r, cdf, BR_CDF_SIZE, ACCT_STR);
                level += coeff_br;
                if (coeff_br < BR_CDF_SIZE - 1) break;
            }
        }
        levels[get_padded_idx(pos, bwl)] = level;
    }

    if (eob > 1) {
        base_cdf_arr base_cdf = frm_ctx->coeff_base_cdf[txs_ctx][plane_type];
        br_cdf_arr br_cdf =
            frm_ctx->coeff_br_cdf[AOMMIN(txs_ctx, TX_32X32)][plane_type];
        if (tx_class == TX_CLASS_2D) {
            read_coeffs_reverse_2d(r, tx_size, 1, eob - 1 - 1, scan, bwl, levels,
                base_cdf, br_cdf);
            read_coeffs_reverse(r, tx_size, tx_class, 0, 0, scan, bwl, levels,
                base_cdf, br_cdf);
        }
        else {
            read_coeffs_reverse(r, tx_size, tx_class, 0, eob - 1 - 1, scan, bwl,
                levels, base_cdf, br_cdf);
        }
    }
#if SVT_DEC_COEFF_DEBUG
    {
        int16_t    *cur_coeff = (int16_t *)coeff_buf;
        /* cur_coeff[0] used for debug */
        cur_coeff[1] = eob;
    }
#else
    coeff_buf[0] = eob;
#endif

    for (int c = 0; c < eob; ++c) {
        const int pos = scan[c];
        uint8_t sign = 0;
        TranLow level = levels[get_padded_idx(pos, bwl)];
        if (level) {
            max_scan_line = AOMMAX(max_scan_line, pos);
            if (c == 0) {
                sign = svt_read_symbol(r, frm_ctx->dc_sign_cdf[plane_type][dc_sign_ctx],
                    2, ACCT_STR);
            }
            else
                sign = svt_read_bit(r, ACCT_STR);
            if (level >= MAX_BASE_BR_RANGE)
                level += read_golomb(r);
            if (c == 0) dc_val = sign ? 1 : 2;

        level &= 0xfffff;
        cul_level += level;
        }
        coeff_buf[c + 1] = sign ? -level : level;
    }

    cul_level = AOMMIN(COEFF_CONTEXT_MASK, cul_level);

    update_coeff_ctx(dec_handle, plane, xd, tx_size, blk_row, blk_col,
        above_off, left_off, cul_level, dc_val);

    trans_info->cbf = 1; assert(eob);

    return eob;
}

void palette_tokens(int plane) {
    assert(plane == 0 || plane == 1);

    // TO-DO
}

void get_palette_color_context()
{
    //TO-DO
}

#define CHECK_BACKWARD_REFS(ref_frame) \
  (((ref_frame) >= BWDREF_FRAME) && ((ref_frame) <= ALTREF_FRAME))
#define IS_BACKWARD_REF_FRAME(ref_frame) CHECK_BACKWARD_REFS(ref_frame)

int get_reference_mode_context(const PartitionInfo_t *xd) {
    int ctx;
    const ModeInfo_t *const above_mbmi = xd->above_mbmi;
    const ModeInfo_t *const left_mbmi = xd->left_mbmi;
    const int has_above = xd->up_available;
    const int has_left = xd->left_available;

    // Note:
    // The mode info data structure has a one element border above and to the
    // left of the entries corresponding to real macroblocks.
    // The prediction flags in these dummy entries are initialized to 0.
    if (has_above && has_left) {  // both edges available
        if (!has_second_ref(above_mbmi) && !has_second_ref(left_mbmi))
            // neither edge uses comp pred (0/1)
            ctx = IS_BACKWARD_REF_FRAME(above_mbmi->ref_frame[0]) ^
            IS_BACKWARD_REF_FRAME(left_mbmi->ref_frame[0]);
        else if (!has_second_ref(above_mbmi))
            // one of two edges uses comp pred (2/3)
            ctx = 2 + (IS_BACKWARD_REF_FRAME(above_mbmi->ref_frame[0]) ||
                !dec_is_inter_block(above_mbmi));
        else if (!has_second_ref(left_mbmi))
            // one of two edges uses comp pred (2/3)
            ctx = 2 + (IS_BACKWARD_REF_FRAME(left_mbmi->ref_frame[0]) ||
                !dec_is_inter_block(left_mbmi));
        else  // both edges use comp pred (4)
            ctx = 4;
    }
    else if (has_above || has_left) {  // one edge available
        const ModeInfo_t *edge_mbmi = has_above ? above_mbmi : left_mbmi;

        if (!has_second_ref(edge_mbmi))
            // edge does not use comp pred (0/1)
            ctx = IS_BACKWARD_REF_FRAME(edge_mbmi->ref_frame[0]);
        else
            // edge uses comp pred (3)
            ctx = 3;
    }
    else {  // no edges available (1)
        ctx = 1;
    }
    assert(ctx >= 0 && ctx < COMP_INTER_CONTEXTS);
    return ctx;
}

void read_ref_frames(EbDecHandle *dec_handle, PartitionInfo_t *const xd,
    SvtReader *r, int segment_id, MvReferenceFrame ref_frame[2])
{
    (void)r;
    //ParseCtxt *parse_ctxt = (ParseCtxt *)dec_handle->pv_parse_ctxt;

    if (xd->mi->skip_mode)
        return;
    if (seg_feature_active(&dec_handle->frame_header.segmentation_params, segment_id, SEG_LVL_REF_FRAME)) {
        ref_frame[0] = (MvReferenceFrame)get_segdata(&dec_handle->frame_header.segmentation_params, segment_id,
            SEG_LVL_REF_FRAME);
        ref_frame[1] = NONE_FRAME;
    }
    else if (seg_feature_active(&dec_handle->frame_header.segmentation_params, segment_id, SEG_LVL_SKIP) ||
        seg_feature_active(&dec_handle->frame_header.segmentation_params, segment_id, SEG_LVL_GLOBALMV)) {
        ref_frame[0] = LAST_FRAME;
        ref_frame[1] = NONE_FRAME;
    }
    else {
        assert(0);
        /* ReferenceMode mode = SINGLE_REFERENCE;
        int bw4 = mi_size_wide[xd->mi->sb_type];
        int bh4 = mi_size_high[xd->mi->sb_type];
        if (dec_handle->frame_header.reference_mode == REFERENCE_MODE_SELECT &&
            (AOMMIN(bw4, bh4) >= 2))
        {
            const int ctx = get_reference_mode_context(xd);
            const ReferenceMode mode = (ReferenceMode)svt_read_symbol(
                r, parse_ctxt->frm_ctx[0].comp_inter_cdf[ctx], 2, ACCT_STR);
        }*/

        //if (mode == COMPOUND_REFERENCE) {
        //    const int ctx = get_comp_reference_type_context(xd);
        //    const COMP_REFERENCE_TYPE comp_ref_type =
        //        (COMP_REFERENCE_TYPE)svt_read_symbol(
        //            r, &parse_ctxt->frm_ctx[0].comp_ref_type_cdf[ctx], 2, ACCT_STR);

        //    if (comp_ref_type == UNIDIR_COMP_REFERENCE) {
        //        const int bit = READ_REF_BIT(uni_comp_ref_p);
        //        if (bit) {
        //            ref_frame[0] = BWDREF_FRAME;
        //            ref_frame[1] = ALTREF_FRAME;
        //        }
        //        else {
        //            const int bit1 = READ_REF_BIT(uni_comp_ref_p1);
        //            if (bit1) {
        //                const int bit2 = READ_REF_BIT(uni_comp_ref_p2);
        //                if (bit2) {
        //                    ref_frame[0] = LAST_FRAME;
        //                    ref_frame[1] = GOLDEN_FRAME;
        //                }
        //                else {
        //                    ref_frame[0] = LAST_FRAME;
        //                    ref_frame[1] = LAST3_FRAME;
        //                }
        //            }
        //            else {
        //                ref_frame[0] = LAST_FRAME;
        //                ref_frame[1] = LAST2_FRAME;
        //            }
        //        }
        //        return;
        //    }

        //    assert(comp_ref_type == BIDIR_COMP_REFERENCE);

        //    const int idx = 1;
        //    const int bit = READ_REF_BIT(comp_ref_p);
        //    // Decode forward references.
        //    if (!bit) {
        //        const int bit1 = READ_REF_BIT(comp_ref_p1);
        //        ref_frame[!idx] = bit1 ? LAST2_FRAME : LAST_FRAME;
        //    }
        //    else {
        //        const int bit2 = READ_REF_BIT(comp_ref_p2);
        //        ref_frame[!idx] = bit2 ? GOLDEN_FRAME : LAST3_FRAME;
        //    }

        //    // Decode backward references.
        //    const int bit_bwd = READ_REF_BIT(comp_bwdref_p);
        //    if (!bit_bwd) {
        //        const int bit1_bwd = READ_REF_BIT(comp_bwdref_p1);
        //        ref_frame[idx] = bit1_bwd? ALTREF2_FRAME : BWDREF_FRAME;
        //    }
        //    else {
        //        ref_frame[idx] = ALTREF_FRAME;
        //    }
        //}
        //else if (mode == SINGLE_REFERENCE) {
        //    const int bit0 = READ_REF_BIT(single_ref_p1);
        //    if (bit0) {
        //        const int bit1 = READ_REF_BIT(single_ref_p2);
        //        if (!bit1) {
        //            const int bit5 = READ_REF_BIT(single_ref_p6);
        //            ref_frame[0] = bit5 ? ALTREF2_FRAME : BWDREF_FRAME;
        //        }
        //        else {
        //            ref_frame[0] = ALTREF_FRAME;
        //        }
        //    }
        //    else {
        //        const int bit2 = READ_REF_BIT(single_ref_p3);
        //        if (bit2) {
        //            const int bit4 = READ_REF_BIT(single_ref_p5);
        //            ref_frame[0] = bit4 ? GOLDEN_FRAME : LAST3_FRAME;
        //        }
        //        else {
        //            const int bit3 = READ_REF_BIT(single_ref_p4);
        //            ref_frame[0] = bit3 ? LAST2_FRAME : LAST_FRAME;
        //        }
        //    }

        //    ref_frame[1] = NONE_FRAME;
        //}
        //else {
        //    assert(0 && "Invalid prediction mode.");
        //}
    }
}

int partition_plane_context(int mi_row,
    int mi_col, BlockSize bsize, ParseCtxt *parse_ctxt) {
    const uint8_t *above_ctx = parse_ctxt->parse_nbr4x4_ctxt.above_part_wd + mi_col;
    const uint8_t *left_ctx =
        parse_ctxt->parse_nbr4x4_ctxt.left_part_ht + ((mi_row- parse_ctxt->sb_row_mi) & MAX_MIB_MASK);

    // Minimum partition point is 8x8. Offset the bsl accordingly.
    int bsl = mi_size_wide_log2[bsize] - mi_size_wide_log2[BLOCK_8X8];
    int above = (*above_ctx >> bsl) & 1, left = (*left_ctx >> bsl) & 1;

    assert(mi_size_wide_log2[bsize] == mi_size_high_log2[bsize]);
    assert(bsl >= 0);

    return (left * 2 + above) + bsl * PARTITION_PLOFFSET;
}

AomCdfProb cdf_element_prob(AomCdfProb *cdf,
    size_t element) {
    assert(cdf != NULL);
    return (element > 0 ? cdf[element - 1] : CDF_PROB_TOP) - cdf[element];
}

void partition_gather_horz_alike(AomCdfProb *out,
     AomCdfProb *in,
    BlockSize bsize) {
    (void)bsize;
    out[0] = CDF_PROB_TOP;
    out[0] -= cdf_element_prob(in, PARTITION_HORZ);
    out[0] -= cdf_element_prob(in, PARTITION_SPLIT);
    out[0] -= cdf_element_prob(in, PARTITION_HORZ_A);
    out[0] -= cdf_element_prob(in, PARTITION_HORZ_B);
    out[0] -= cdf_element_prob(in, PARTITION_VERT_A);
    if (bsize != BLOCK_128X128) out[0] -= cdf_element_prob(in, PARTITION_HORZ_4);
    out[0] = AOM_ICDF(out[0]);
    out[1] = AOM_ICDF(CDF_PROB_TOP);
}

void partition_gather_vert_alike(AomCdfProb *out,
     AomCdfProb *in,
    BlockSize bsize) {
    (void)bsize;
    out[0] = CDF_PROB_TOP;
    out[0] -= cdf_element_prob(in, PARTITION_VERT);
    out[0] -= cdf_element_prob(in, PARTITION_SPLIT);
    out[0] -= cdf_element_prob(in, PARTITION_HORZ_A);
    out[0] -= cdf_element_prob(in, PARTITION_VERT_A);
    out[0] -= cdf_element_prob(in, PARTITION_VERT_B);
    if (bsize != BLOCK_128X128) out[0] -= cdf_element_prob(in, PARTITION_VERT_4);
    out[0] = AOM_ICDF(out[0]);
    out[1] = AOM_ICDF(CDF_PROB_TOP);
}

PartitionType parse_partition_type(uint32_t blk_row, uint32_t blk_col, SvtReader *reader,
    BlockSize bsize, int has_rows, int has_cols, EbDecHandle *dec_handle)
{
    ParseCtxt *parse_ctxt = (ParseCtxt *)dec_handle->pv_parse_ctxt;

    int partition_cdf_length = bsize <= BLOCK_8X8 ? PARTITION_TYPES :
        (bsize == BLOCK_128X128 ? EXT_PARTITION_TYPES - 2 : EXT_PARTITION_TYPES);
    int ctx = partition_plane_context(blk_row, blk_col, bsize, parse_ctxt);

    if (bsize < BLOCK_8X8) return PARTITION_NONE;
    else if (has_rows && has_cols)
    {
        return (PartitionType)svt_read_symbol(
            reader, parse_ctxt->frm_ctx[0].partition_cdf[ctx], partition_cdf_length, ACCT_STR);
    }
    else if (has_cols)
    {
        assert(bsize > BLOCK_8X8);
        AomCdfProb cdf[2];
        partition_gather_vert_alike(cdf, parse_ctxt->frm_ctx[0].partition_cdf[ctx], bsize);
        assert(cdf[1] == AOM_ICDF(CDF_PROB_TOP));
        return svt_read_cdf(reader, cdf, 2, ACCT_STR) ? PARTITION_SPLIT : PARTITION_HORZ;
    }
    else if (has_rows)
    {
        assert(has_rows && !has_cols);
        assert(bsize > BLOCK_8X8);
        AomCdfProb cdf[2];
        partition_gather_horz_alike(cdf, parse_ctxt->frm_ctx[0].partition_cdf[ctx], bsize);
        assert(cdf[1] == AOM_ICDF(CDF_PROB_TOP));
        return svt_read_cdf(reader, cdf, 2, ACCT_STR) ? PARTITION_SPLIT : PARTITION_VERT;
    }
    return PARTITION_INVALID;
}

static INLINE void dec_get_txb_ctx(int plane_bsize, const TxSize tx_size, const int plane,
                               int blk_row, int blk_col,
                               EbDecHandle *dec_handle, TXB_CTX *const txb_ctx)
{
#define MAX_TX_SIZE_UNIT 16

    ParseCtxt *parse_ctx = (ParseCtxt*)dec_handle->pv_parse_ctxt;
    ParseNbr4x4Ctxt *nbr_ctx = &parse_ctx->parse_nbr4x4_ctxt;
    EbColorConfig *clr_cfg = &dec_handle->seq_header.color_config;
    int txb_w_unit = tx_size_wide_unit[tx_size];
    int txb_h_unit = tx_size_high_unit[tx_size];
    uint8_t suby   = plane ? clr_cfg->subsampling_y : 0;

    int dc_sign = 0;
    int k = 0;
    uint8_t *above_dc_ctx = nbr_ctx->above_dc_ctx[plane] + blk_col;
    uint8_t *left_dc_ctx = nbr_ctx->left_dc_ctx[plane] + (blk_row - (parse_ctx->sb_row_mi>> suby));

    do {
        const unsigned int sign = above_dc_ctx[k];
        if (sign == 1)
            dc_sign--;
        else if (sign == 2)
            dc_sign++;
    } while (++k < txb_w_unit);

    k = 0;
    do {
        const unsigned int sign = left_dc_ctx[k];
        if (sign == 1)
            dc_sign--;
        else if (sign == 2)
            dc_sign++;
    } while (++k < txb_h_unit);

    if (dc_sign < 0)
        txb_ctx->dc_sign_ctx = 1;
    else if (dc_sign > 0)
        txb_ctx->dc_sign_ctx = 2;
    else
        txb_ctx->dc_sign_ctx = 0;

    uint8_t *above_level_ctx = nbr_ctx->above_level_ctx[plane] + blk_col;
    uint8_t *left_level_ctx = nbr_ctx->left_level_ctx[plane] + (blk_row - (parse_ctx->sb_row_mi >> suby));

    if (plane == 0) {
        if (plane_bsize == txsize_to_bsize[tx_size])
            txb_ctx->txb_skip_ctx = 0;
        else {
            static const uint8_t skip_contexts[5][5] = { { 1, 2, 2, 2, 3 },
                                                         { 1, 4, 4, 4, 5 },
                                                         { 1, 4, 4, 4, 5 },
                                                         { 1, 4, 4, 4, 5 },
                                                         { 1, 4, 4, 4, 6 } };
            int top = 0;
            int left = 0;

            k = 0;
            do {
                top |= above_level_ctx[k];
            } while (++k < txb_w_unit);

            k = 0;
            do {
                left |= left_level_ctx[k];
            } while (++k < txb_h_unit);
            const int max = AOMMIN(top | left, 4);
            const int min = AOMMIN(AOMMIN(top, left), 4);

            txb_ctx->txb_skip_ctx = skip_contexts[min][max];
        }
    }
    else {
        //TODO check this - change to libAOM implementation
        int above = 0, left = 0;
        k = 0;
        do {
                above |= above_level_ctx[k];
                above |= above_dc_ctx[k];
        } while (++k < txb_w_unit);

        k = 0;
        do {
            left |= left_level_ctx[k];
            left |= left_dc_ctx[k];
        } while (++k < txb_h_unit);

        const int ctx_offset = (num_pels_log2_lookup[plane_bsize] >
            num_pels_log2_lookup[txsize_to_bsize[tx_size]])
            ? 10
            : 7;
        int ctx = (above != 0) + (left != 0);
        txb_ctx->txb_skip_ctx = ctx + ctx_offset;
    }
#undef MAX_TX_SIZE_UNIT
}

uint16_t parse_transform_block(EbDecHandle *dec_handle, PartitionInfo_t *pi, SvtReader *r,
                              int32_t *coeff, TransformInfo_t *trans_info, int plane,
                              int blk_col, int blk_row, BlockSize mi_row, BlockSize mi_col,
                              TxSize tx_size, int is_inter, int skip)
{
    // ParseCtxt *parse_ctxt = (ParseCtxt *)dec_handle->pv_parse_ctxt;
    (void)is_inter;
    uint16_t eob = 0 , sub_x, sub_y;
    BlockSize bsize = pi->mi->sb_type;

    sub_x = (plane > 0) ? dec_handle->seq_header.color_config.subsampling_x : 0;
    sub_y = (plane > 0) ? dec_handle->seq_header.color_config.subsampling_y : 0;

    uint32_t start_x = (mi_col >> sub_x) + blk_col;
    uint32_t start_y = (mi_row >> sub_y) + blk_row;

    if ( start_x >= (dec_handle->frame_header.mi_cols >> sub_x) ||
         start_y >= (dec_handle->frame_header.mi_rows >> sub_y) )
        return eob;

    if (!skip) {
        int plane_bsize = (bsize == BLOCK_INVALID) ? BLOCK_INVALID :
            ss_size_lookup[bsize][sub_x][sub_y];
        TXB_CTX txb_ctx;
        dec_get_txb_ctx(plane_bsize, tx_size, plane,
                    start_y, start_x,
                    dec_handle, &txb_ctx);

        eob = parse_coeffs(dec_handle, pi, r, start_y, start_x, blk_col, blk_row,
            plane, txb_ctx.txb_skip_ctx, txb_ctx.dc_sign_ctx, tx_size, coeff, trans_info);
    }
    else{
        update_coeff_ctx(dec_handle, plane, pi, tx_size,
            start_y, start_x, blk_col, blk_row, 0, 0);
    }
    return eob;
}

/* Gives the pointer to current block's transform info. Will work only for Intra */
static INLINE TransformInfo_t* get_cur_trans_info_intra(int plane, SBInfo  *sb_info, ModeInfo_t *mi){
    return (plane == 0) ? (sb_info->sb_luma_trans_info + mi->first_luma_tu_offset) :
                          (sb_info->sb_chroma_trans_info + mi->first_chroma_tu_offset + (plane-1));
}

/* TODO: Clean the logic! */
static INLINE TxSize av1_get_tx_size(int plane, PartitionInfo_t *pi, int sub_x, int sub_y) {
    /*const */ModeInfo_t *mbmi = pi->mi;
    if (plane == 0) return (get_cur_trans_info_intra(plane, pi->sb_info,mbmi)->tx_size);
    return av1_get_max_uv_txsize(mbmi->sb_type, sub_x, sub_y);
}

void parse_residual(EbDecHandle *dec_handle, PartitionInfo_t *pi, SvtReader *r,
                    int mi_row, int mi_col, BlockSize mi_size)
{
    ParseCtxt   *parse_ctx = (ParseCtxt*)dec_handle->pv_parse_ctxt;
    TxSize      tx_size;
    EbColorConfig *color_info = &dec_handle->seq_header.color_config;
    int num_planes = color_info->mono_chrome ? 1 : MAX_MB_PLANE;
    ModeInfo_t *mode = pi->mi;

    int is_inter = dec_is_inter_block(mode);
    int skip     = mode->skip;
    int lossless = (&dec_handle->frame_header.lossless_array[0])[mode->segment_id];

    if (is_inter)
        assert(0);
    else {
        const int max_blocks_wide = max_block_wide(pi, mi_size, 0);
        const int max_blocks_high = max_block_high(pi, mi_size, 0);
        const BlockSize max_unit_bsize = BLOCK_64X64;
        int mu_blocks_wide = block_size_wide[max_unit_bsize] >> tx_size_wide_log2[0];
        int mu_blocks_high = block_size_high[max_unit_bsize] >> tx_size_high_log2[0];
        mu_blocks_wide = AOMMIN(max_blocks_wide, mu_blocks_wide);
        mu_blocks_high = AOMMIN(max_blocks_high, mu_blocks_high);
        TransformInfo_t *trans_info = NULL, *next_trans_chroma = NULL, *next_trans_luma = NULL;

        for (int row = 0; row < max_blocks_high; row += mu_blocks_high) {
            for (int col = 0; col < max_blocks_wide; col += mu_blocks_wide) {
                for (int plane = 0; plane < num_planes; ++plane) {
                    int sub_x = (plane > 0) ? color_info->subsampling_x : 0;
                    int sub_y = (plane > 0) ? color_info->subsampling_y : 0;

                    if (!dec_is_chroma_reference(mi_row, mi_col, mi_size, sub_x, sub_y))
                        continue;

                    if (row == 0 && col == 0)
                    trans_info = get_cur_trans_info_intra(plane, pi->sb_info, mode);
                    else
                        trans_info = plane ? next_trans_chroma : next_trans_luma;

                    if (lossless)
                        tx_size = TX_4X4;
                    else
                        tx_size = av1_get_tx_size(plane, pi, sub_x, sub_y);

                    int step_x = tx_size_high[tx_size] >> 2;
                    int step_y = tx_size_wide[tx_size] >> 2;

                    const int unit_height = ROUND_POWER_OF_TWO(
                        AOMMIN(mu_blocks_high + row, max_blocks_high), sub_y);
                    const int unit_width = ROUND_POWER_OF_TWO(
                        AOMMIN(mu_blocks_wide + col, max_blocks_wide), sub_x);

                    for (int blk_row = row >> sub_y; blk_row < unit_height; blk_row += step_x) {
                        for (int blk_col = col >> sub_x; blk_col < unit_width; blk_col += step_y) {
                            int32_t *coeff = plane ? parse_ctx->cur_chroma_coeff_buf :
                                                     parse_ctx->cur_luma_coeff_buf;
#if SVT_DEC_COEFF_DEBUG
                            {
                            uint8_t  *cur_coeff = (uint8_t*)coeff;
                            uint8_t  cur_loc = (mi_row + blk_row) & 0xFF;
                            cur_coeff[0] = cur_loc;
                            cur_loc = (mi_col + blk_col) & 0xFF;
                            cur_coeff[1] = cur_loc;
                            }
#endif
                            int32_t eob = parse_transform_block(dec_handle, pi, r, coeff,
                                trans_info, plane,
                                blk_col, blk_row, mi_row, mi_col,
                                tx_size, is_inter, skip);

                            if (eob != 0) {
                                plane ? (parse_ctx->cur_chroma_coeff_buf += (eob + 1)) : (parse_ctx->cur_luma_coeff_buf += (eob + 1));
                                trans_info->cbf = 1;
                                }
                            else
                                trans_info->cbf = 0;

                            // increment transform pointer
                            trans_info++;
                        } //for blk_col
                    } //for blk_row

                    // Remembers trans_info for next transform block within a block of 128xH / Wx128
                    if (plane == 0)
                        next_trans_luma = trans_info;
                    else
                        next_trans_chroma = trans_info;
                }//for plane
            }
        }
    }//intra
}

void parse_block(EbDecHandle *dec_handle,
    uint32_t mi_row, uint32_t mi_col, SvtReader *r,
    BlockSize subsize, TileInfo *tile,
    SBInfo *sb_info,
    PartitionType partition)
{
    ParseCtxt *parse_ctx = (ParseCtxt*)dec_handle->pv_parse_ctxt;

    ModeInfo_t *mode = parse_ctx->cur_mode_info;

    int8_t      *cdef_strength = sb_info->sb_cdef_strength;

    int bw4 = mi_size_wide[subsize];
    int bh4 = mi_size_high[subsize];

    uint32_t mi_cols = (&dec_handle->frame_header)->mi_cols;
    uint32_t mi_rows = (&dec_handle->frame_header)->mi_rows;

    /* TODO: Can move to a common init fun for parse & decode */
    PartitionInfo_t part_info;
    part_info.mi        = mode;
    part_info.sb_info   = sb_info;
    part_info.mi_row = mi_row;
    part_info.mi_col = mi_col;

    mode->partition = partition;
    /* TU offset update from parse ctxt info of previous block */
    mode->first_luma_tu_offset  = parse_ctx->first_luma_tu_offset;
    mode->first_chroma_tu_offset= parse_ctx->first_chroma_tu_offset;
#if MODE_INFO_DBG
    mode->mi_row = mi_row;
    mode->mi_col = mi_col;
#endif
    EbColorConfig color_config = dec_handle->seq_header.color_config;
    if (bh4 == 1 && color_config.subsampling_y && (mi_row & 1) == 0)
        part_info.has_chroma = 0;
    else if (bw4 == 1 && color_config.subsampling_x && (mi_col & 1) == 0)
        part_info.has_chroma = 0;
    else
        part_info.has_chroma = color_config.mono_chrome ? 1: MAX_MB_PLANE;

    /* TODO : tile->tile_rows boundary condn check is wrong */
    part_info.up_available = ((int32_t)mi_row > tile->mi_row_start);
    part_info.left_available = ((int32_t)mi_col > tile->mi_col_start);
    part_info.chroma_up_available = part_info.up_available;
    part_info.chroma_left_available = part_info.left_available;

    part_info.mb_to_right_edge = ((int32_t)mi_cols - bw4 - mi_col) * MI_SIZE * 8;
    part_info.mb_to_bottom_edge = ((int32_t)mi_rows - bh4 - mi_row) * MI_SIZE * 8;

    if (part_info.has_chroma)
    {
        if (bh4 == 1 && color_config.subsampling_y)
            part_info.chroma_up_available = (int32_t)(mi_row - 1) > tile->mi_row_start;
        if (bw4 == 1 && color_config.subsampling_x)
            part_info.chroma_left_available = (int32_t)(mi_col - 1) > tile->mi_col_start;
    }
    else
    {
        part_info.chroma_up_available = 0;
        part_info.chroma_left_available = 0;
    }

    if (part_info.up_available)
        part_info.above_mbmi = get_top_mode_info(dec_handle, mi_row, mi_col, sb_info);
    else
        part_info.above_mbmi = NULL;
    if (part_info.left_available)
        part_info.left_mbmi = get_left_mode_info(dec_handle, mi_row, mi_col, sb_info);
    else
        part_info.left_mbmi = NULL;
    mode->sb_type = subsize;
    mode_info(dec_handle, &part_info, mi_row, mi_col, r, cdef_strength);

    /* Replicating same chroma mode for block pairs or 4x4 blks when chroma is present in last block*/
    if(0 == parse_ctx->prev_blk_has_chroma)
    {
        /* if the previous block does not have chroma info then         */
        /* current uv mode is stores in the previous ModeInfo structre. */
        /* this is done to simplify neighbour access for mode deviation */
        if (mode->sb_type != BLOCK_4X4)
        {
            //2 partition case
            mode[-1].uv_mode = mode->uv_mode;
            assert(part_info.has_chroma != 0);
        }
        else
        {
            /* wait unitl the last 4x4 block is parsed */
            if (part_info.has_chroma)
            {
                //4 partition case
                mode[-1].uv_mode = mode->uv_mode;
                mode[-2].uv_mode = mode->uv_mode;
                mode[-3].uv_mode = mode->uv_mode;
            }
        }
    }

    /* current block's has_chroma info is stored for useage in next block */
    parse_ctx->prev_blk_has_chroma = part_info.has_chroma;

    if (!dec_is_inter_block(mode)) {
        for (int plane = 0; plane < AOMMIN(2, dec_handle->seq_header.color_config.mono_chrome ? 1 : MAX_MB_PLANE); ++plane)
            palette_tokens(plane);
    }

    read_block_tx_size(dec_handle, r, &part_info, subsize);

    parse_residual(dec_handle, &part_info, r, mi_row, mi_col, subsize);

    /* Update block level MI map */
    update_block_nbrs(dec_handle, mi_row, mi_col, subsize);
    parse_ctx->cur_mode_info_cnt++;
    parse_ctx->cur_mode_info++;
}

static INLINE void update_partition_context(ParseCtxt *parse_ctx,
    int mi_row, int mi_col, BlockSize subsize,
    BlockSize bsize) {
    ParseNbr4x4Ctxt *ngr_ctx = &parse_ctx->parse_nbr4x4_ctxt;
    uint8_t *const above_ctx = ngr_ctx->above_part_wd + mi_col;
    uint8_t *const left_ctx = ngr_ctx->left_part_ht + ((mi_row- parse_ctx->sb_row_mi) & MAX_MIB_MASK);

    const int bw = mi_size_wide[bsize];
    const int bh = mi_size_high[bsize];
    memset(above_ctx, partition_context_lookup[subsize].above, bw);
    memset(left_ctx, partition_context_lookup[subsize].left, bh);
}

static INLINE void update_ext_partition_context(ParseCtxt *parse_ctx, int mi_row,
    int mi_col, BlockSize subsize,
    BlockSize bsize,
    PartitionType partition) {
    if (bsize >= BLOCK_8X8) {
        const int hbs = mi_size_wide[bsize] / 2;
        BlockSize bsize2 = Partition_Subsize[PARTITION_SPLIT][bsize];
        switch (partition) {
        case PARTITION_SPLIT:
            if (bsize != BLOCK_8X8) break;
        case PARTITION_NONE:
        case PARTITION_HORZ:
        case PARTITION_VERT:
        case PARTITION_HORZ_4:
        case PARTITION_VERT_4:
            update_partition_context(parse_ctx, mi_row, mi_col, subsize, bsize);
            break;
        case PARTITION_HORZ_A:
            update_partition_context(parse_ctx, mi_row, mi_col, bsize2, subsize);
            update_partition_context(parse_ctx, mi_row + hbs, mi_col, subsize, subsize);
            break;
        case PARTITION_HORZ_B:
            update_partition_context(parse_ctx, mi_row, mi_col, subsize, subsize);
            update_partition_context(parse_ctx, mi_row + hbs, mi_col, bsize2, subsize);
            break;
        case PARTITION_VERT_A:
            update_partition_context(parse_ctx, mi_row, mi_col, bsize2, subsize);
            update_partition_context(parse_ctx, mi_row, mi_col + hbs, subsize, subsize);
            break;
        case PARTITION_VERT_B:
            update_partition_context(parse_ctx, mi_row, mi_col, subsize, subsize);
            update_partition_context(parse_ctx, mi_row, mi_col + hbs, bsize2, subsize);
            break;
        default: assert(0 && "Invalid partition type");
        }
    }
}

void parse_partition(EbDecHandle *dec_handle, uint32_t blk_row,
    uint32_t blk_col, SvtReader *reader, BlockSize bsize, SBInfo *sb_info)
{
    ParseCtxt *parse_ctx = (ParseCtxt*)dec_handle->pv_parse_ctxt;

    if (blk_row >= dec_handle->frame_header.mi_rows ||
        blk_col >= dec_handle->frame_header.mi_cols)
        return;

    int num4x4 = mi_size_wide[bsize];
    int half_block_4x4 = num4x4 >> 1;
    int quarter_block_4x4 = half_block_4x4 >> 1;

    int has_rows = (blk_row + half_block_4x4) < dec_handle->frame_header.mi_rows;
    int has_cols = (blk_col + half_block_4x4) < dec_handle->frame_header.mi_cols;

    PartitionType partition;

    partition = (bsize < BLOCK_8X8) ? PARTITION_NONE
        : parse_partition_type(blk_row, blk_col, reader, bsize, has_rows, has_cols, dec_handle);
    int subSize = Partition_Subsize[(int)partition][bsize];
    int splitSize = Partition_Subsize[PARTITION_SPLIT][bsize];

#define PARSE_BLOCK(db_r, db_c, db_subsize)                 \
parse_block(dec_handle, db_r, db_c, reader, db_subsize,     \
    &parse_ctx->cur_tile_info, sb_info, partition);

#define PARSE_PARTITION(db_r, db_c, db_subsize)                 \
  parse_partition(dec_handle, (db_r), (db_c), reader,           \
                   (db_subsize), sb_info)

    switch ((int)partition) {
    case PARTITION_NONE: PARSE_BLOCK(blk_row, blk_col, subSize); break;
    case PARTITION_HORZ:
        PARSE_BLOCK(blk_row, blk_col, subSize);
        if (has_rows) PARSE_BLOCK(blk_row + half_block_4x4, blk_col, subSize);
        break;
    case PARTITION_VERT:
        PARSE_BLOCK(blk_row, blk_col, subSize);
        if (has_cols) PARSE_BLOCK(blk_row, blk_col + half_block_4x4, subSize);
        break;
    case PARTITION_SPLIT:
        PARSE_PARTITION(blk_row, blk_col, subSize);
        PARSE_PARTITION(blk_row, blk_col + half_block_4x4, subSize);
        PARSE_PARTITION(blk_row + half_block_4x4, blk_col, subSize);
        PARSE_PARTITION(blk_row + half_block_4x4, blk_col + half_block_4x4, subSize);
        break;
    case PARTITION_HORZ_A:
        PARSE_BLOCK(blk_row, blk_col, splitSize);
        PARSE_BLOCK(blk_row, blk_col + half_block_4x4, splitSize);
        PARSE_BLOCK(blk_row + half_block_4x4, blk_col, subSize);
        break;
    case PARTITION_HORZ_B:
        PARSE_BLOCK(blk_row, blk_col, subSize);
        PARSE_BLOCK(blk_row + half_block_4x4, blk_col, splitSize);
        PARSE_BLOCK(blk_row + half_block_4x4, blk_col + half_block_4x4, splitSize);
        break;
    case PARTITION_VERT_A:
        PARSE_BLOCK(blk_row, blk_col, splitSize);
        PARSE_BLOCK(blk_row + half_block_4x4, blk_col, splitSize);
        PARSE_BLOCK(blk_row, blk_col + half_block_4x4, subSize);
        break;
    case PARTITION_VERT_B:
        PARSE_BLOCK(blk_row, blk_col, subSize);
        PARSE_BLOCK(blk_row, blk_col + half_block_4x4, splitSize);
        PARSE_BLOCK(blk_row + half_block_4x4, blk_col + half_block_4x4, splitSize);
        break;
    case PARTITION_HORZ_4:
        for (int i = 0; i < 4; ++i) {
            uint32_t this_blk_row = blk_row + (uint32_t)(i * quarter_block_4x4);
            if (i > 0 && this_blk_row >= dec_handle->frame_header.mi_rows) break;
            PARSE_BLOCK(this_blk_row, blk_col, subSize);
        }
        break;
    case PARTITION_VERT_4:
        for (int i = 0; i < 4; ++i) {
            uint32_t this_blk_col = blk_col + (uint32_t)(i * quarter_block_4x4);
            if (i > 0 && this_blk_col >= dec_handle->frame_header.mi_cols) break;
            PARSE_BLOCK(blk_row, this_blk_col, subSize);
        }
        break;
    default: assert(0 && "Invalid partition type");
    }
    update_ext_partition_context(parse_ctx, blk_row, blk_col, subSize, bsize, partition);
}

void parse_super_block(EbDecHandle *dec_handle,
    uint32_t blk_row, uint32_t blk_col, SBInfo *sbInfo)
{
    ParseCtxt *parse_ctx = (ParseCtxt*)dec_handle->pv_parse_ctxt;
    SvtReader *reader = &parse_ctx->r;

    parse_partition(dec_handle, blk_row, blk_col, reader,
        dec_handle->seq_header.sb_size, sbInfo);
}
