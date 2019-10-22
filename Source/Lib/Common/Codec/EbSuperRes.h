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

#ifndef EbSuperRes_h
#define EbSuperRes_h

#include "EbAv1Structs.h"

#ifdef __cplusplus
extern "C" {
#endif

    void av1_tile_set_col(TilesInfo *tile, TilesInfo *tiles_info, int32_t mi_cols,
        uint32_t *mi_col_start, uint32_t *mi_col_end, int col);

#ifdef __cplusplus
}
#endif
#endif // EbSuperRes_h
