/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbEntropyCodingProcess_h
#define EbEntropyCodingProcess_h

#include "EbDefinitions.h"
#include "EbSyntaxElements.h"

#include "EbSystemResourceManager.h"
#include "EbPictureBufferDesc.h"
#include "EbEncInterPrediction.h"
#include "EbEntropyCoding.h"
#include "EbTransQuantBuffers.h"
#include "EbCodingUnit.h"
#include "EbObject.h"

/**************************************
 * Enc Dec Context
 **************************************/
typedef struct EntropyCodingContext {
    EbDctor  dctor;
    EbFifo * enc_dec_input_fifo_ptr;
    EbFifo * entropy_coding_output_fifo_ptr; // to packetization
    EbFifo * rate_control_output_fifo_ptr; // feedback to rate control
    uint32_t sb_total_count;
    // Coding Unit Workspace---------------------------
    EbPictureBufferDesc *coeff_buffer_sb; //Used to hold quantized coeff for one TB in EncPass.

    //  Context Variables---------------------------------
    BlkStruct *blk_ptr;
    //const CodedBlockStats           *cu_stats;
    uint32_t        blk_index;
    uint8_t         cu_depth;
    uint32_t        cu_size;
    uint32_t        sb_sz;
    uint32_t        cu_size_log2;
    uint32_t        blk_origin_x;
    uint32_t        blk_origin_y;
    uint32_t        sb_origin_x;
    uint32_t        sb_origin_y;
    uint32_t        pu_itr;
    PredictionUnit *pu_ptr;
    uint32_t        pu_origin_x;
    uint32_t        pu_origin_y;
    uint32_t        pu_width;
    uint32_t        pu_height;
    MvUnit          mv_unit;

    uint32_t       txb_itr;
    uint32_t       txb_origin_x;
    uint32_t       txb_origin_y;
    uint32_t       txb_size;

    // MCP Context
    EbBool      is_16bit; //enable 10 bit encode in CL
    int32_t     coded_area_sb;
    int32_t     coded_area_sb_uv;
    TOKENEXTRA *tok;
} EntropyCodingContext;

/**************************************
 * Extern Function Declarations
 **************************************/
extern EbErrorType entropy_coding_context_ctor(EbThreadContext *  thread_context_ptr,
                                               const EbEncHandle *enc_handle_ptr, int index,
                                               int rate_control_index);

extern void *entropy_coding_kernel(void *input_ptr);

#endif // EbEntropyCodingProcess_h
