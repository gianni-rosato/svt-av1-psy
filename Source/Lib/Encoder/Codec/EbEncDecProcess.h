/*
* Copyright(c) 2019 Intel Corporation
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#ifndef EbEncDecProcess_h
#define EbEncDecProcess_h

#include "EbDefinitions.h"
#include "EbModeDecisionProcess.h"
#include "EbSystemResourceManager.h"
#include "EbPictureBufferDesc.h"
#include "EbModeDecision.h"
#include "EbEncInterPrediction.h"
#include "EbEntropyCoding.h"
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
    EbFifo              *mode_decision_input_fifo_ptr;
    EbFifo              *enc_dec_output_fifo_ptr;
    EbFifo              *enc_dec_feedback_fifo_ptr;
    EbFifo              *picture_demux_output_fifo_ptr; // to picture-manager
    ModeDecisionContext *md_ctx;
    const BlockGeom     *blk_geom;
    // Coding Unit Workspace---------------------------
    EbPictureBufferDesc *residual_buffer;
    EbPictureBufferDesc *transform_buffer;
    EbPictureBufferDesc *input_samples;
    EbPictureBufferDesc *input_sample16bit_buffer;
    // temporary buffers for decision making of LF (LPF_PICK_FROM_FULL_IMAGE).
    // Since recon switches between reconPtr and referencePtr, the temporary buffers sizes used the referencePtr's which has padding,...
    EbPictureBufferDesc *inverse_quant_buffer;
    uint32_t             pic_fast_lambda[2];
    uint32_t             pic_full_lambda[2];

    //  Context Variables---------------------------------
    BlkStruct *blk_ptr;
    //const CodedBlockStats                *cu_stats;
    uint16_t      blk_org_x; // within the picture
    uint16_t      blk_org_y; // within the picture
    uint32_t      sb_index;
    MvUnit        mv_unit;
    uint8_t       txb_itr;
    Bool          is_16bit; //enable 10 bit encode in CL
    uint32_t      bit_depth;
    EbColorFormat color_format;
    uint64_t      tot_intra_coded_area;
    uint64_t      tot_skip_coded_area;
    uint64_t      three_quad_energy;

    // Needed for DC prediction
    uint8_t  upsample_left;
    uint8_t  upsample_above;
    uint16_t coded_area_sb;
    uint16_t coded_area_sb_uv;

    uint8_t is_inter;
    uint8_t reduced_tx_set_used;
    uint8_t md_skip_blk;

    uint16_t tile_group_index;
    uint16_t tile_index;
    uint32_t coded_sb_count;
} EncDecContext;

/**************************************
 * Extern Function Declarations
 **************************************/
extern EbErrorType svt_aom_enc_dec_context_ctor(EbThreadContext *thread_ctx, const EbEncHandle *enc_handle_ptr,
                                                int index, int tasks_index);

extern void *svt_aom_mode_decision_kernel(void *input_ptr);

#ifdef __cplusplus
}
#endif
#endif // EbEncDecProcess_h
