/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbModeDecision_h
#define EbModeDecision_h

#include "EbDefinitions.h"
#include "EbUtility.h"
#include "EbPictureControlSet.h"
#include "EbCodingUnit.h"
#include "EbPredictionUnit.h"
#include "EbSyntaxElements.h"
#include "EbPictureBufferDesc.h"
#include "EbAdaptiveMotionVectorPrediction.h"
#include "EbPictureOperators.h"
#include "EbNeighborArrays.h"
#ifdef __cplusplus
extern "C" {
#endif
#define ENABLE_AMVP_MV_FOR_RC_PU    0
#define MAX_MB_PLANE                3
#define MAX_MPM_CANDIDATES          3
#define MERGE_PENALTY                          10

    // Create incomplete struct definition for the following function pointer typedefs
    struct ModeDecisionCandidateBuffer_s;
    struct ModeDecisionContext_s;

    /**************************************
    * Function Ptrs Definitions
    **************************************/
    typedef EbErrorType(*EB_PREDICTION_FUNC)(
        struct ModeDecisionContext_s           *context_ptr,
        uint32_t                                component_mask,
        PictureControlSet_t                    *picture_control_set_ptr,
        struct ModeDecisionCandidateBuffer_s   *candidate_buffer_ptr,
        EbAsm                                   asm_type);

    typedef EbErrorType(*EB_FAST_COST_FUNC)(
        struct ModeDecisionContext_s           *context_ptr,
        CodingUnit_t                           *cu_ptr,
        struct ModeDecisionCandidateBuffer_s   *candidate_buffer_ptr,
        uint32_t                                qp,
        uint64_t                                luma_distortion,
        uint64_t                                chroma_distortion,
        uint64_t                                lambda,
        PictureControlSet_t                    *picture_control_set_ptr);

    typedef EbErrorType(*EB_FULL_COST_FUNC)(
        LargestCodingUnit_t                    *sb_ptr,
        CodingUnit_t                           *cu_ptr,
        uint32_t                                cu_size,
        uint32_t                                cu_size_log2,
        struct ModeDecisionCandidateBuffer_s   *candidate_buffer_ptr,
        uint32_t                                qp,
        uint64_t                               *y_distortion,
        uint64_t                               *cb_distortion,
        uint64_t                               *cr_distortion,
        uint64_t                                lambda,
        uint64_t                                lambda_chroma,
        uint64_t                               *y_coeff_bits,
        uint64_t                               *cb_coeff_bits,
        uint64_t                               *cr_coeff_bits,
        uint32_t                                transform_size,
        uint32_t                                transform_chroma_size,
        PictureControlSet_t                    *picture_control_set_ptr);

    typedef EbErrorType(*EB_AV1_FULL_COST_FUNC)(
        PictureControlSet_t                    *picture_control_set_ptr,
        struct ModeDecisionContext_s           *context_ptr,
        struct ModeDecisionCandidateBuffer_s   *candidate_buffer_ptr,
        CodingUnit_t                           *cu_ptr,
        uint64_t                               *y_distortion,
        uint64_t                               *cb_distortion,
        uint64_t                               *cr_distortion,
        uint64_t                                lambda,
        uint64_t                               *y_coeff_bits,
        uint64_t                               *cb_coeff_bits,
        uint64_t                               *cr_coeff_bits,
        BlockSize                               bsize);

    typedef EbErrorType(*EB_FULL_LUMA_COST_FUNC)(
        CodingUnit_t                           *cu_ptr,
        uint32_t                                cu_size,
        uint32_t                                cu_size_log2,
        struct ModeDecisionCandidateBuffer_s   *candidate_buffer_ptr,
        uint64_t                               *y_distortion,
        uint64_t                                lambda,
        uint64_t                               *y_coeff_bits,
        uint32_t                                transform_size);

    /**************************************
    * Mode Decision Candidate
    **************************************/
    typedef struct ModeDecisionCandidate_s
    {
        // *Warning - this struct has been organized to be cache efficient when being
        //    constructured in the function GenerateAmvpMergeInterIntraMdCandidatesCU.
        //    Changing the ordering could affect performance
        union {
            struct {
                unsigned                        me_distortion    : 20;
                unsigned                        distortion_ready : 1;
                unsigned : 3;
                unsigned                        intra_luma_mode  : 8; // HEVC mode, use pred_mode for AV1
            };
            uint32_t ois_results;
        };
        union {
            struct {
                union {
                    struct {
                        int16_t              motionVector_x_L0;  //Note: Do not change the order of these fields
                        int16_t              motionVector_y_L0;
                    };
                    uint32_t MVsL0;
                };
                union {
                    struct {
                        int16_t              motionVector_x_L1;  //Note: Do not change the order of these fields
                        int16_t              motionVector_y_L1;
                    };
                    uint32_t MVsL1;
                };
            };
            uint64_t MVs;
        };

        uint8_t                                skip_flag;
        EbBool                                 merge_flag;
        uint8_t                                merge_index;
        uint16_t                               count_non_zero_coeffs;
        EbBool                                 prediction_is_ready_luma;
        uint8_t                                type;
        EbBool                                 mpm_flag;

        // MD Rate Estimation Ptr
        MdRateEstimationContext_t             *md_rate_estimation_ptr; // 64 bits
        uint64_t                               fast_luma_rate;
        uint64_t                               fast_chroma_rate;
        uint64_t                               chroma_distortion;
        uint64_t                               chroma_distortion_inter_depth;
        uint32_t                               luma_distortion;
        uint32_t                               full_distortion;

        EbPtr                                 prediction_context_ptr;
        PictureControlSet_t                   *picture_control_set_ptr;
        EbPredDirection                        prediction_direction[MAX_NUM_OF_PU_PER_CU]; // 2 bits

        int16_t                                motion_vector_pred_x[MAX_NUM_OF_REF_PIC_LIST]; // 16 bits
        int16_t                                motion_vector_pred_y[MAX_NUM_OF_REF_PIC_LIST]; // 16 bits
        uint8_t                                motion_vector_pred_idx[MAX_NUM_OF_REF_PIC_LIST]; // 2 bits
        uint8_t                                block_has_coeff;             // ?? bit - determine empirically
        uint8_t                                u_has_coeff;               // ?? bit
        uint8_t                                v_has_coeff;               // ?? bit
        uint32_t                               y_has_coeff;                // Issue, should be less than 32

        PredictionMode                         pred_mode; // AV1 mode, no need to convert
        uint8_t                                drl_index;

        // Intra Mode
        int32_t                                angle_delta[PLANE_TYPES];
        EbBool                                 is_directional_mode_flag;
        EbBool                                 is_directional_chroma_mode_flag;
        EbBool                                 use_angle_delta;
        uint32_t                               intra_chroma_mode; // AV1 mode, no need to convert

        // Index of the alpha Cb and alpha Cr combination
        int32_t                                cfl_alpha_idx;
        // Joint sign of alpha Cb and alpha Cr
        int32_t                                cfl_alpha_signs;

        // Inter Mode
        PredictionMode                         inter_mode;
        EbBool                                 is_compound;
        uint32_t                               pred_mv_weight;
        uint8_t                                ref_frame_type;
        uint8_t                                ref_mv_index;
        EbBool                                 is_skip_mode_flag;
        EbBool                                 is_new_mv;
        EbBool                                 is_zero_mv;
        TxType                                 transform_type[PLANE_TYPES];
        MacroblockPlane                        candidate_plane[MAX_MB_PLANE];
        uint16_t                               eob[MAX_MB_PLANE][MAX_TXB_COUNT];
        int32_t                                quantized_dc[3];
        uint32_t                               interp_filters;
        uint8_t                                tu_width;
        uint8_t                                tu_height;
        MOTION_MODE                            motion_mode;
        uint16_t                               num_proj_ref;
        EbBool                                 local_warp_valid;
        EbWarpedMotionParams                   wm_params;
    } ModeDecisionCandidate_t;

    /**************************************
    * Mode Decision Candidate Buffer
    **************************************/
    typedef struct IntraChromaCandidateBuffer_s {
        uint32_t                                mode;
        uint64_t                                cost;
        uint64_t                                distortion;
        EbPictureBufferDesc_t                  *prediction_ptr;
        EbPictureBufferDesc_t                  *residual_ptr;
    } IntraChromaCandidateBuffer_t;

    /**************************************
    * Mode Decision Candidate Buffer
    **************************************/
    typedef struct ModeDecisionCandidateBuffer_s {
        // Candidate Ptr
        ModeDecisionCandidate_t                *candidate_ptr;

        // Video Buffers
        EbPictureBufferDesc_t                  *prediction_ptr;
        EbPictureBufferDesc_t                  *predictionPtrTemp;
        EbPictureBufferDesc_t                  *cflTempPredictionPtr;
        EbPictureBufferDesc_t                  *residualQuantCoeffPtr;// One buffer for residual and quantized coefficient
        EbPictureBufferDesc_t                  *reconCoeffPtr;
        EbPictureBufferDesc_t                  *residual_ptr;

        // *Note - We should be able to combine the reconCoeffPtr & reconPtr pictures (they aren't needed at the same time)
        EbPictureBufferDesc_t                  *reconPtr;

        // Distortion (SAD)
        uint64_t                                residual_luma_sad;
        uint64_t                                full_lambda_rate;
        uint64_t                                full_cost_luma;
                                               
        // Costs                               
        uint64_t                               *fast_cost_ptr;
        uint64_t                               *full_cost_ptr;
        uint64_t                               *full_cost_skip_ptr;
        uint64_t                               *full_cost_merge_ptr;
        //                                     
        uint64_t                                cb_coeff_bits;
        uint64_t                                cb_distortion[2];
        uint64_t                                cr_coeff_bits;
        uint64_t                                cr_distortion[2];
        EbBool                                  sub_sampled_pred; //do prediction every other line, duplicate the residual
        EbBool                                  sub_sampled_pred_chroma;
        uint64_t                                y_full_distortion[DIST_CALC_TOTAL];
        uint64_t                                y_coeff_bits;

    } ModeDecisionCandidateBuffer_t;

    /**************************************
    * Extern Function Declarations
    **************************************/
    extern EbErrorType mode_decision_candidate_buffer_ctor(
        ModeDecisionCandidateBuffer_t **buffer_dbl_ptr,
        uint16_t                        sb_max_size,
        EB_BITDEPTH                     max_bitdepth,
        uint64_t                       *fast_cost_ptr,
        uint64_t                       *full_cost_ptr,
        uint64_t                       *full_cost_skip_ptr,
        uint64_t                       *full_cost_merge_ptr
    );

    uint8_t product_full_mode_decision(
        struct ModeDecisionContext_s   *context_ptr,
        CodingUnit_t                   *cu_ptr,
        uint8_t                         bwidth,
        uint8_t                         bheight,
        ModeDecisionCandidateBuffer_t **buffer_ptr_array,
        uint32_t                        candidate_total_count,
        uint8_t                        *best_candidate_index_array,
        uint32_t                       *best_intra_mode);

    EbErrorType PreModeDecision(
        CodingUnit_t                   *cu_ptr,
        uint32_t                        buffer_total_count,
        ModeDecisionCandidateBuffer_t **buffer_ptr_array,
        uint32_t                       *full_candidate_total_count_ptr,
        uint8_t                        *best_candidate_index_array,
        uint8_t                        *disable_merge_index,
#if TX_SEARCH_LEVELS
        uint64_t                       *ref_fast_cost,
#endif
        EbBool                          same_fast_full_candidate);

    typedef EbErrorType(*EB_INTRA_4x4_FAST_LUMA_COST_FUNC)(
        struct ModeDecisionContext_s           *context_ptr,
        uint32_t                                pu_index,
        struct ModeDecisionCandidateBuffer_s   *candidate_buffer_ptr,
        uint64_t                                luma_distortion,
        uint64_t                                lambda);

    typedef EbErrorType(*EB_INTRA_4x4_FULL_LUMA_COST_FUNC)(
        struct ModeDecisionCandidateBuffer_s   *candidate_buffer_ptr,
        uint64_t                               *y_distortion,
        uint64_t                                lambda,
        uint64_t                               *y_coeff_bits,
        uint32_t                                transform_size);

    typedef EbErrorType(*EB_FULL_NXN_COST_FUNC)(
        PictureControlSet_t                    *picture_control_set_ptr,
        struct ModeDecisionCandidateBuffer_s   *candidate_buffer_ptr,
        uint32_t                                qp,
        uint64_t                               *y_distortion,
        uint64_t                               *cb_distortion,
        uint64_t                               *cr_distortion,
        uint64_t                                lambda,
        uint64_t                                lambda_chroma,
        uint64_t                               *y_coeff_bits,
        uint64_t                               *cb_coeff_bits,
        uint64_t                               *cr_coeff_bits,
        uint32_t                                transform_size);

    struct CodingLoopContext_s;
#ifdef __cplusplus
}
#endif
#endif // EbModeDecision_h
