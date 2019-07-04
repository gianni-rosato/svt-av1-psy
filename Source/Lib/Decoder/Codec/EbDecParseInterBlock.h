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

#ifndef EbDecParseInterBlock_h
#define EbDecParseInterBlock_h

#ifdef __cplusplus
extern "C" {
#endif

#include "EbDefinitions.h"
#include "EbPictureBufferDesc.h"

#include "EbSvtAv1Dec.h"
#include "EbDecHandle.h"
#include "EbDecBitReader.h"
#include "EbObuParse.h"

#include "EbDecParseHelper.h"

#include "EbDecNbr.h"
#include "EbDecPicMgr.h"
#include "EbDecUtils.h"

#define MVREF_ROW_COLS 3
#define MV_BORDER (16 << 3)  // Allow 16 pels in 1/8th pel units
#define MAX_DIFFWTD_MASK_BITS 1
#define NELEMENTS(x) (int)(sizeof(x) / sizeof(x[0]))

static const MV kZeroMv = { 0, 0 };

extern  int8_t av1_ref_frame_type(const MvReferenceFrame *const rf);
extern void av1_set_ref_frame(MvReferenceFrame *rf, int8_t ref_frame_type);

void inter_block_mode_info(EbDecHandle *dec_handle, PartitionInfo_t* pi,
    int mi_row, int mi_col, SvtReader *r);

void av1_find_mv_refs(EbDecHandle *dec_handle, PartitionInfo_t *pi,
    MvReferenceFrame ref_frame, CandidateMv_dec ref_mv_stack[][MAX_REF_MV_STACK_SIZE],
    IntMv_dec mv_ref_list[][MAX_MV_REF_CANDIDATES], IntMv_dec global_mvs[2],
    int mi_row, int mi_col, int16_t *mode_context, MvCount *mv_cnt);

#ifdef __cplusplus
}
#endif
#endif // EbDecParseInterBlock_h
