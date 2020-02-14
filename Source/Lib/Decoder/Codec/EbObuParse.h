/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbDecObuParser_h
#define EbDecObuParser_h

#include "assert.h"
#include "EbCabacContextModel.h" //for ENTROPY_CONTEXT
#include "EbDecBitstream.h"
#include "EbDecBitReader.h"

#define HEADER_DUMP 0

#if HEADER_DUMP
#define PRINT_NL SVT_LOG("\n");
#define PRINT(name, val) SVT_LOG("\n%s :\t%X", name, val);
#define PRINT_NAME(name) SVT_LOG("\n%s :\t", name);
#define PRINT_FRAME(name, val) SVT_LOG("\n%s :\t%X", name, val);
#else
#define PRINT_NL
#define PRINT(name, val)
#define PRINT_NAME(name)
#define PRINT_FRAME(name, val)
#endif

#define MAX_NUM_TEMPORAL_LAYERS 8
#define MAX_NUM_SPATIAL_LAYERS 4
#ifdef MAX_NUM_OPERATING_POINTS
#undef MAX_NUM_OPERATING_POINTS
#endif // MAX_NUM_OPERATING_POINTS

#define MAX_NUM_OPERATING_POINTS MAX_NUM_TEMPORAL_LAYERS *MAX_NUM_SPATIAL_LAYERS
#define OP_POINTS_CNT_MINUS_1_BITS 5
#define OP_POINTS_IDC_BITS 12

#define SELECT_SCREEN_CONTENT_TOOLS 2
#define SELECT_INTEGER_MV 2
#define NUM_REF_FRAMES 8

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

int get_qindex(SegmentationParams *seg_params, int segment_id, int base_q_idx);
void svt_setup_motion_field(EbDecHandle *dec_handle, DecThreadCtxt *thread_ctxt);
EbErrorType decode_multiple_obu(EbDecHandle *dec_handle_ptr, uint8_t **data, size_t data_size,
                                uint32_t is_annexb);

static INLINE int allow_intrabc(const EbDecHandle *dec_handle) {
    return (dec_handle->frame_header.frame_type == KEY_FRAME ||
            dec_handle->frame_header.frame_type == INTRA_ONLY_FRAME) &&
           dec_handle->seq_header.seq_force_screen_content_tools &&
           dec_handle->frame_header.allow_intrabc;
}

#endif // EbDecObuParser_h
