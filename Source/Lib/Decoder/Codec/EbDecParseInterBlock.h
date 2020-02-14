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
#include "EbDecInterPrediction.h"

#define MVREF_ROW_COLS 3
// Set the upper limit of the motion vector component magnitude.
// This would make a motion vector fit in 26 bits. Plus 3 bits for the
// reference frame index. A tuple of motion vector can hence be stored within
// 32 bit range for efficient load/store operations.
#define REFMVS_LIMIT ((1 << 12) - 1)

static const MV k_zero_mv = {0, 0};

extern int8_t av1_ref_frame_type(const MvReferenceFrame *const rf);
extern void   av1_set_ref_frame(MvReferenceFrame *rf, int8_t ref_frame_type);

void inter_block_mode_info(EbDecHandle *dec_handle, ParseCtxt *parse_ctxt, PartitionInfo *pi);

void av1_find_mv_refs(EbDecHandle *dec_handle, PartitionInfo *pi, ParseCtxt *parse_ctx,
                      MvReferenceFrame ref_frame, CandidateMv ref_mv_stack[][MAX_REF_MV_STACK_SIZE],
                      IntMv mv_ref_list[][MAX_MV_REF_CANDIDATES], IntMv global_mvs[2],
                      int16_t *mode_context, MvCount *mv_cnt);
void assign_intrabc_mv(ParseCtxt *parse_ctxt, IntMv ref_mvs[INTRA_FRAME + 1][MAX_MV_REF_CANDIDATES],
                       PartitionInfo *pi);
void palette_tokens(EbDecHandle *dec_handle, ParseCtxt *parse_ctx, PartitionInfo *pi);
#ifdef __cplusplus
}
#endif
#endif // EbDecParseInterBlock_h
