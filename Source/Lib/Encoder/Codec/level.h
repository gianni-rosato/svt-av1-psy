/*
 * Copyright (c) 2019, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef AOM_AV1_ENCODER_LEVEL_H_
#define AOM_AV1_ENCODER_LEVEL_H_

#include "EbDefinitions.h"

#if TWOPASS_RC
enum {
    SEQ_LEVEL_2_0,
    SEQ_LEVEL_2_1,
    SEQ_LEVEL_2_2,
    SEQ_LEVEL_2_3,
    SEQ_LEVEL_3_0,
    SEQ_LEVEL_3_1,
    SEQ_LEVEL_3_2,
    SEQ_LEVEL_3_3,
    SEQ_LEVEL_4_0,
    SEQ_LEVEL_4_1,
    SEQ_LEVEL_4_2,
    SEQ_LEVEL_4_3,
    SEQ_LEVEL_5_0,
    SEQ_LEVEL_5_1,
    SEQ_LEVEL_5_2,
    SEQ_LEVEL_5_3,
    SEQ_LEVEL_6_0,
    SEQ_LEVEL_6_1,
    SEQ_LEVEL_6_2,
    SEQ_LEVEL_6_3,
    SEQ_LEVEL_7_0,
    SEQ_LEVEL_7_1,
    SEQ_LEVEL_7_2,
    SEQ_LEVEL_7_3,
    SEQ_LEVELS,
    SEQ_LEVEL_MAX = 31
} UENUM1BYTE(AV1_LEVEL);


typedef BitstreamProfile BITSTREAM_PROFILE;

// AV1 Level Specifications
typedef struct {
  AV1_LEVEL level;
  int max_picture_size;
  int max_h_size;
  int max_v_size;
  int max_header_rate;
  int max_tile_rate;
  int max_tiles;
  int max_tile_cols;
  int64_t max_display_rate;
  int64_t max_decode_rate;
  double main_mbps;
  double high_mbps;
  double main_cr;
  double high_cr;
} AV1LevelSpec;

typedef struct AV1LevelParams {
  // Specifies the level that the coded video sequence conforms to for each
  // operating point.
  //AV1_LEVEL target_seq_level_idx[MAX_NUM_OPERATING_POINTS];
  // Bit mask to indicate whether to keep level stats for corresponding
  // operating points.
  uint32_t keep_level_stats;
  // Level information for each operating point.
  //AV1LevelInfo *level_info[MAX_NUM_OPERATING_POINTS];
  // Count the number of OBU_FRAME and OBU_FRAME_HEADER for level calculation.
  int frame_header_count;
} AV1LevelParams;

#endif  // TWOPASS_RC
#endif  // AOM_AV1_ENCODER_LEVEL_H_
