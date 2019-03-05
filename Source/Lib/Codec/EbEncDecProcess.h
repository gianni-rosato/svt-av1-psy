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
#include "EbInterPrediction.h"
#include "EbEntropyCoding.h"
#include "EbTransQuantBuffers.h"
#include "EbReferenceObject.h"
#include "EbNeighborArrays.h"
#include "EbCodingUnit.h"

#ifdef __cplusplus
extern "C" {
#endif

#define PM_STRIDE   4

    typedef struct EbPmCand_s
    {
        int16_t        tr_coeff[4 * 4];
        int16_t        qu_coeff[4 * 4];
        int16_t        iq_coeff[4 * 4];
        uint8_t        masking_level;
        uint64_t       cost;
        uint32_t       nz_coeff;
    } EbPmCand_t;

    /**************************************
     * Enc Dec Context
     **************************************/
    typedef struct EncDecContext_s
    {
        EbFifo_t                              *mode_decision_input_fifo_ptr;
        EbFifo_t                              *enc_dec_output_fifo_ptr;
        EbFifo_t                              *enc_dec_feedback_fifo_ptr;
        EbFifo_t                              *picture_demux_output_fifo_ptr;   // to picture-manager
        int16_t                               *transform_inner_array_ptr;
        MdRateEstimationContext_t             *md_rate_estimation_ptr;
        ModeDecisionContext_t                 *md_context;
        const BlockGeom                       *blk_geom;

        // TMVP
        EbReferenceObject_t                   *reference_object_write_ptr;

        // MCP Context
        MotionCompensationPredictionContext_t *mcp_context;
        SsMeContext_t                         *ss_mecontext;

        // Intra Reference Samples
        IntraReferenceSamples_t               *intra_ref_ptr;
        IntraReference16bitSamples_t          *intra_ref_ptr16;  //We need a different buffer for ENC pass then the MD one.
        
        // Coding Unit Workspace---------------------------
        EbPictureBufferDesc_t                 *residual_buffer;
        EbPictureBufferDesc_t                 *transform_buffer;
        EbPictureBufferDesc_t                 *input_samples;
        EbPictureBufferDesc_t                 *input_sample16bit_buffer;
#if !FILT_PROC
        EbPictureBufferDesc_t                 *trial_frame_rst;
#endif
        // temporary buffers for decision making of LF (LPF_PICK_FROM_FULL_IMAGE).
        // Since recon switches between reconPtr and referencePtr, the temporary buffers sizes used the referencePtr's which has padding,...
#if !FILT_PROC
        EbPictureBufferDesc_t                 *temp_lf_recon_picture_ptr;
        EbPictureBufferDesc_t                 *temp_lf_recon_picture16bit_ptr;
#endif
        EbPictureBufferDesc_t                 *inverse_quant_buffer;
        // Lambda
#if ADD_DELTA_QP_SUPPORT
        uint16_t                               qp;
#else
        uint8_t                                qp;
#endif
        uint8_t                                chroma_qp;
        uint32_t                               fast_lambda;
        uint32_t                               full_lambda;
        uint32_t                               fast_chroma_lambda;
        uint32_t                               full_chroma_lambda;
        uint32_t                               full_chroma_lambda_sao;

        //  Context Variables---------------------------------
        CodingUnit_t                          *cu_ptr;
        const CodedUnitStats_t                *cu_stats;
        uint16_t                               cu_origin_x; // within the picture
        uint16_t                               cu_origin_y; // within the picture
        uint8_t                                sb_sz;
        uint32_t                               sb_index;
        MvUnit_t                               mv_unit;
        int16_t                                x_mv_amvp_candidate_array_list0[MAX_NUM_OF_AMVP_CANDIDATES];
        uint8_t                                txb_itr;
        EbBool                                 is16bit; //enable 10 bit encode in CL
        uint64_t                               tot_intra_coded_area;
        uint8_t                                intra_coded_area_sb[MAX_NUMBER_OF_TREEBLOCKS_PER_PICTURE];//percentage of intra coded area 0-100%
        uint8_t                                pmp_masking_level_enc_dec;
        EbBool                                 skip_qpm_flag;
        int16_t                                min_delta_qp_weight;
        int16_t                                max_delta_qp_weight;
        int8_t                                 min_delta_qp[4];
        int8_t                                 max_delta_qp[4];
        int8_t                                 non_moving_delta_qp;
        EbBool                                 grass_enhancement_flag;
        EbBool                                 backgorund_enhancement;
#if ADD_DELTA_QP_SUPPORT                       
        uint16_t                               qpmQp;
#else                                          
        uint8_t                                qpmQp;
#endif
        EB_TRANS_COEFF_SHAPE                   trans_coeff_shape_luma;
        EB_TRANS_COEFF_SHAPE                   trans_coeff_shape_chroma;
        EbPmCand_t                             pm_cand_buffer[5];
        uint16_t                               qp_index;
        uint64_t                               three_quad_energy;

        // Needed for DC prediction
        EbBool                                 is_left_availble;
        EbBool                                 is_above_availble;
        uint8_t                                upsample_left;
        uint8_t                                upsample_above;
        uint8_t                                upsample_left_chroma;
        uint8_t                                upsample_above_chroma; 
#if !CHROMA_BLIND
        int16_t                                pred_buf_q3[CFL_BUF_SQUARE];
#endif
        uint16_t                               coded_area_sb;
        uint16_t                               coded_area_sb_uv;

#if ENCDEC_TX_SEARCH
        uint8_t                                is_inter;
        uint8_t                                reduced_tx_set_used;
#endif
#if CHROMA_BLIND
        EbBool                                 evaluate_cfl_ep; // 0: CFL is evaluated @ mode decision, 1: CFL is evaluated @ encode pass
#endif
    } EncDecContext_t;

    /**************************************
     * Extern Function Declarations
     **************************************/
    extern EbErrorType enc_dec_context_ctor(
        EncDecContext_t        **context_dbl_ptr,
        EbFifo_t                *mode_decision_configuration_input_fifo_ptr,
        EbFifo_t                *packetization_output_fifo_ptr,
        EbFifo_t                *feedback_fifo_ptr,
        EbFifo_t                *picture_demux_fifo_ptr,
        EbBool                   is16bit,
        uint32_t                 max_input_luma_width,
        uint32_t                 max_input_luma_height);

    extern void* EncDecKernel(void *input_ptr);

#ifdef __cplusplus
}
#endif
#endif // EbEncDecProcess_h