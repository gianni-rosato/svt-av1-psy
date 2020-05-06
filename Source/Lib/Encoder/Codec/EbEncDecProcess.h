/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbEncDecProcess_h
#define EbEncDecProcess_h

#include "EbDefinitions.h"
#include "EbSyntaxElements.h"
#include "EbModeDecisionProcess.h"
#include "EbSystemResourceManager.h"
#include "EbPictureBufferDesc.h"
#include "EbModeDecision.h"
#include "EbEncInterPrediction.h"
#include "EbEntropyCoding.h"
#include "EbTransQuantBuffers.h"
#include "EbReferenceObject.h"
#include "EbNeighborArrays.h"
#include "EbCodingUnit.h"
#include "EbObject.h"

#ifdef __cplusplus
extern "C" {
#endif

#define PM_STRIDE 4

/**************************************
     * Enc Dec Context
     **************************************/
typedef struct EncDecContext {
    EbFifo *                 mode_decision_input_fifo_ptr;
    EbFifo *                 enc_dec_output_fifo_ptr;
    EbFifo *                 enc_dec_feedback_fifo_ptr;
    EbFifo *                 picture_demux_output_fifo_ptr; // to picture-manager
    MdRateEstimationContext *md_rate_estimation_ptr;
    EbBool                   is_md_rate_estimation_ptr_owner;
    ModeDecisionContext *    md_context;
    const BlockGeom *        blk_geom;
    // MCP Context
    MotionCompensationPredictionContext *mcp_context;

    // Coding Unit Workspace---------------------------
    EbPictureBufferDesc *residual_buffer;
    EbPictureBufferDesc *transform_buffer;
    EbPictureBufferDesc *input_samples;
    EbPictureBufferDesc *input_sample16bit_buffer;
    // temporary buffers for decision making of LF (LPF_PICK_FROM_FULL_IMAGE).
    // Since recon switches between reconPtr and referencePtr, the temporary buffers sizes used the referencePtr's which has padding,...
    EbPictureBufferDesc *inverse_quant_buffer;
    // Lambda
    uint16_t qp;
    uint8_t  chroma_qp;
    uint32_t fast_lambda;
    uint32_t full_lambda;
    uint32_t full_chroma_lambda_sao;

    //  Context Variables---------------------------------
    BlkStruct *blk_ptr;
    //const CodedBlockStats                *cu_stats;
    uint16_t      blk_origin_x; // within the picture
    uint16_t      blk_origin_y; // within the picture
    uint8_t       sb_sz;
    uint32_t      sb_index;
    MvUnit        mv_unit;
    uint8_t       txb_itr;
    EbBool        is_16bit; //enable 10 bit encode in CL
    uint32_t      bit_depth;
    EbColorFormat color_format;
    uint64_t      tot_intra_coded_area;
    uint8_t       intra_coded_area_sb
        [MAX_NUMBER_OF_TREEBLOCKS_PER_PICTURE]; //percentage of intra coded area 0-100%
    uint16_t qp_index;
    uint64_t three_quad_energy;

    // Needed for DC prediction
    uint8_t upsample_left;
    uint8_t upsample_above;
    uint16_t coded_area_sb;
    uint16_t coded_area_sb_uv;

    uint8_t is_inter;
    uint8_t reduced_tx_set_used;
    EbBool
            evaluate_cfl_ep; // 0: CFL is evaluated @ mode decision, 1: CFL is evaluated @ encode pass
    uint8_t md_skip_blk;

    uint16_t tile_group_index;
    uint16_t tile_index;
    uint32_t coded_sb_count;
} EncDecContext;

/**************************************
     * Extern Function Declarations
     **************************************/
extern EbErrorType enc_dec_context_ctor(EbThreadContext *  thread_context_ptr,
                                        const EbEncHandle *enc_handle_ptr, int index,
                                        int tasks_index, int demux_index);

extern void *enc_dec_kernel(void *input_ptr);

#ifdef __cplusplus
}
#endif
#endif // EbEncDecProcess_h
