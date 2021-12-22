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

#include "level.h"

#define UNDEFINED_LEVEL                                                                    \
    {                                                                                      \
        .level = SEQ_LEVEL_MAX, .max_picture_size = 0, .max_h_size = 0, .max_v_size = 0,   \
        .max_display_rate = 0, .max_decode_rate = 0, .max_header_rate = 0, .main_mbps = 0, \
        .high_mbps = 0, .main_cr = 0, .high_cr = 0, .max_tiles = 0, .max_tile_cols = 0     \
    }

typedef enum {
    LUMA_PIC_SIZE_TOO_LARGE,
    LUMA_PIC_H_SIZE_TOO_LARGE,
    LUMA_PIC_V_SIZE_TOO_LARGE,
    LUMA_PIC_H_SIZE_TOO_SMALL,
    LUMA_PIC_V_SIZE_TOO_SMALL,
    TOO_MANY_TILE_COLUMNS,
    TOO_MANY_TILES,
    TILE_RATE_TOO_HIGH,
    TILE_TOO_LARGE,
    SUPERRES_TILE_WIDTH_TOO_LARGE,
    CROPPED_TILE_WIDTH_TOO_SMALL,
    CROPPED_TILE_HEIGHT_TOO_SMALL,
    TILE_WIDTH_INVALID,
    FRAME_HEADER_RATE_TOO_HIGH,
    DISPLAY_RATE_TOO_HIGH,
    DECODE_RATE_TOO_HIGH,
    CR_TOO_SMALL,
    TILE_SIZE_HEADER_RATE_TOO_HIGH,
    BITRATE_TOO_HIGH,
    DECODER_MODEL_FAIL,

    TARGET_LEVEL_FAIL_IDS,
    TARGET_LEVEL_OK,
} TARGET_LEVEL_FAIL_ID;
