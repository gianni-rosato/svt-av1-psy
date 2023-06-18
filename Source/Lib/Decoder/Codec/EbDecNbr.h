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

#ifndef EbDecNbr_h
#define EbDecNbr_h

#include "EbObuParse.h"
#include "EbDecParseFrame.h"
BlockModeInfo *svt_aom_get_cur_mode_info(void *pv_dec_handle, int mi_row, int mi_col, SBInfo *sb_info);

void svt_aom_update_block_nbrs(EbDecHandle *dec_handle, ParseCtxt *parse_ctx, int mi_row, int mi_col,
                               BlockSize subsize);

BlockModeInfo *svt_aom_get_left_mode_info(EbDecHandle *dec_handle, int mi_row, int mi_col, SBInfo *sb_info);

BlockModeInfo *svt_aom_get_top_mode_info(EbDecHandle *dec_handle, int mi_row, int mi_col, SBInfo *sb_info);

#endif //EbDecNbr_h
