/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbDecNbr_h
#define EbDecNbr_h

#include "EbObuParse.h"
#include "EbDecParseFrame.h"
BlockModeInfo *get_cur_mode_info(void *pv_dec_handle, int mi_row, int mi_col, SBInfo *sb_info);

void update_block_nbrs(EbDecHandle *dec_handle, ParseCtxt *parse_ctx, int mi_row, int mi_col,
                       BlockSize subsize);

BlockModeInfo *get_left_mode_info(EbDecHandle *dec_handle, int mi_row, int mi_col, SBInfo *sb_info);

BlockModeInfo *get_top_mode_info(EbDecHandle *dec_handle, int mi_row, int mi_col, SBInfo *sb_info);

#endif //EbDecNbr_h
