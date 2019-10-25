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
#include "EbModeDecision.h"
#include "EbInterPrediction.h"
#include "EbEntropyCoding.h"
#include "EbTransQuantBuffers.h"
#include "EbReferenceObject.h"
#include "EbNeighborArrays.h"
#include "EbCodingUnit.h"
#include "EbObject.h"

/**************************************
 * Enc Dec Context
 **************************************/
typedef struct EntropyCodingContext
{
    EbDctor                      dctor;
    EbFifo                       *enc_dec_input_fifo_ptr;
    EbFifo                       *entropy_coding_output_fifo_ptr;  // to packetization
    EbFifo                       *rate_control_output_fifo_ptr; // feedback to rate control

    uint32_t                        sb_total_count;
    // Lambda
#if ADD_DELTA_QP_SUPPORT
#else
    uint8_t                         qp;
    uint8_t                         chroma_qp;
#endif
    // Coding Unit Workspace---------------------------
    EbPictureBufferDesc           *coeff_buffer_sb;                              //Used to hold quantized coeff for one TB in EncPass.

    //  Context Variables---------------------------------
    CodingUnit                     *cu_ptr;
    const CodedUnitStats           *cu_stats;
    uint32_t                          cu_index;
    uint8_t                           cu_depth;
    uint32_t                          cu_size;
    uint32_t                          sb_sz;
    uint32_t                          cu_size_log2;
    uint32_t                          cu_origin_x;
    uint32_t                          cu_origin_y;
    uint32_t                          sb_origin_x;
    uint32_t                          sb_origin_y;
    uint32_t                          pu_itr;
    PredictionUnit                 *pu_ptr;
    const PredictionUnitStats      *pu_stats;
    uint32_t                          pu_origin_x;
    uint32_t                          pu_origin_y;
    uint32_t                          pu_width;
    uint32_t                          pu_height;
    MvUnit                          mv_unit;

    uint32_t                          txb_itr;
    TransformUnit                  *txb_ptr;
    uint32_t                          txb_origin_x;
    uint32_t                          txb_origin_y;
    uint32_t                          txb_size;

    // MCP Context
    EbBool                            is16bit; //enable 10 bit encode in CL
    int32_t                           coded_area_sb;
    int32_t                           coded_area_sb_uv;
#if PAL_SUP
    TOKENEXTRA *tok;
#endif
} EntropyCodingContext;

/**************************************
 * Extern Function Declarations
 **************************************/
extern EbErrorType entropy_coding_context_ctor(
    EntropyCodingContext  *context_ptr,
    EbFifo                *enc_dec_input_fifo_ptr,
    EbFifo                *packetization_output_fifo_ptr,
    EbFifo                *rate_control_output_fifo_ptr,
    EbBool                   is16bit);

extern void* entropy_coding_kernel(void *input_ptr);

#endif // EbEntropyCodingProcess_h
