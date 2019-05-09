/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbModeDecisionProcess_h
#define EbModeDecisionProcess_h

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

#ifdef __cplusplus
extern "C" {
#endif
    /**************************************
     * Defines
     **************************************/
#define IBC_CAND 2 //two intra bc candidates
#if CHECK_CAND
#if MRP_DUPLICATION_FIX
#if EIGTH_PEL_MV
#define MODE_DECISION_CANDIDATE_MAX_COUNT               (470+IBC_CAND )
#else
#define MODE_DECISION_CANDIDATE_MAX_COUNT               (440 +IBC_CAND) 
#endif
#else
#define MODE_DECISION_CANDIDATE_MAX_COUNT               (170 +IBC_CAND) 
#endif
#else
#define MODE_DECISION_CANDIDATE_MAX_COUNT               (124+IBC_CAND) /* 61 Intra & 18+2x8+2x8 Inter*/
#endif
#define DEPTH_ONE_STEP   21
#define DEPTH_TWO_STEP    5
#define DEPTH_THREE_STEP  1

     /**************************************
      * Macros
      **************************************/

#define GROUP_OF_4_8x8_BLOCKS(origin_x,origin_y) (((origin_x >> 3) & 0x1) && ((origin_y >> 3) & 0x1) ? EB_TRUE : EB_FALSE)
#define GROUP_OF_4_16x16_BLOCKS(origin_x,origin_y) (((((origin_x >> 3) & 0x2) == 0x2) && (((origin_y >> 3) & 0x2) == 0x2)) ? EB_TRUE : EB_FALSE)
#define GROUP_OF_4_32x32_BLOCKS(origin_x,origin_y) (((((origin_x >> 3) & 0x4) == 0x4) && (((origin_y >> 3) & 0x4) == 0x4)) ? EB_TRUE : EB_FALSE)

      /**************************************
       * Coding Loop Context
       **************************************/
    typedef struct MdEncPassCuData
    {
        uint64_t                    skip_cost;
        uint64_t                    merge_cost;
        uint64_t                    chroma_distortion;
        uint64_t                    y_full_distortion[DIST_CALC_TOTAL];
        uint64_t                    y_coeff_bits;
        uint32_t                    y_has_coeff;
        uint64_t                    fast_luma_rate;
        uint16_t                    y_count_non_zero_coeffs[4];// Store nonzero CoeffNum, per TU. If one TU, stored in 0, otherwise 4 tus stored in 0 to 3

    } MdEncPassCuData;


    typedef struct MdCodingUnit
    {
        unsigned                    tested_cu_flag                  : 1;   //tells whether this CU is tested in MD.
        unsigned                    mdc_array_index                 : 7;
        unsigned                    count_non_zero_coeffs           : 11;
#if M8_SKIP_BLK
        unsigned                    top_neighbor_depth              : 8;
        unsigned                    left_neighbor_depth             : 8;
#else
        unsigned                    top_neighbor_depth              : 3;
        unsigned                    left_neighbor_depth             : 3;
#endif
        unsigned                    top_neighbor_mode               : 2;
        unsigned                    left_neighbor_mode              : 2;
        unsigned                    full_distortion                 : 32;
        unsigned                    chroma_distortion               : 32;
        unsigned                    chroma_distortion_inter_depth   : 32;
        PartitionContextType           left_neighbor_partition;
        PartitionContextType           above_neighbor_partition;
        uint64_t                    cost;
        uint64_t                    cost_luma;
        CandidateMv ed_ref_mv_stack[MODE_CTX_REF_FRAMES][MAX_REF_MV_STACK_SIZE];//to be used in MD and EncDec
#if RED_CU
        uint8_t                     avail_blk_flag ;   //tells whether this CU is tested in MD and have a valid cu data
#endif

    } MdCodingUnit;


    typedef struct ModeDecisionContext
    {
        EbFifo                       *mode_decision_configuration_input_fifo_ptr;
        EbFifo                       *mode_decision_output_fifo_ptr;
        int16_t                        *transform_inner_array_ptr;

        ModeDecisionCandidate       **fast_candidate_ptr_array;
        ModeDecisionCandidate        *fast_candidate_array;
        ModeDecisionCandidateBuffer **candidate_buffer_ptr_array;
        MdRateEstimationContext      *md_rate_estimation_ptr;
        InterPredictionContext       *inter_prediction_context;
        MdCodingUnit                  md_local_cu_unit[BLOCK_MAX_COUNT_SB_128];
        CodingUnit                    md_cu_arr_nsq[BLOCK_MAX_COUNT_SB_128];

        NeighborArrayUnit            *intra_luma_mode_neighbor_array;
        NeighborArrayUnit            *intra_chroma_mode_neighbor_array;
        NeighborArrayUnit            *mv_neighbor_array;
        NeighborArrayUnit            *skip_flag_neighbor_array;
        NeighborArrayUnit            *mode_type_neighbor_array;
        NeighborArrayUnit            *leaf_depth_neighbor_array;
        NeighborArrayUnit            *luma_recon_neighbor_array;
        NeighborArrayUnit            *cb_recon_neighbor_array;
        NeighborArrayUnit            *cr_recon_neighbor_array;
#if !REMOVE_SKIP_COEFF_NEIGHBOR_ARRAY
        NeighborArrayUnit            *skip_coeff_neighbor_array;
#endif
        NeighborArrayUnit            *luma_dc_sign_level_coeff_neighbor_array; // Stored per 4x4. 8 bit: lower 6 bits (COEFF_CONTEXT_BITS), shows if there is at least one Coef. Top 2 bit store the sign of DC as follow: 0->0,1->-1,2-> 1
        NeighborArrayUnit            *cr_dc_sign_level_coeff_neighbor_array; // Stored per 4x4. 8 bit: lower 6 bits(COEFF_CONTEXT_BITS), shows if there is at least one Coef. Top 2 bit store the sign of DC as follow: 0->0,1->-1,2-> 1
        NeighborArrayUnit            *cb_dc_sign_level_coeff_neighbor_array; // Stored per 4x4. 8 bit: lower 6 bits(COEFF_CONTEXT_BITS), shows if there is at least one Coef. Top 2 bit store the sign of DC as follow: 0->0,1->-1,2-> 1
        NeighborArrayUnit            *inter_pred_dir_neighbor_array;
        NeighborArrayUnit            *ref_frame_type_neighbor_array;
        NeighborArrayUnit            *leaf_partition_neighbor_array;
        NeighborArrayUnit32          *interpolation_type_neighbor_array;
#if !OPT_LOSSLESS_0
        // TMVP
        EbReferenceObject            *reference_object_write_ptr;
#endif
        // Intra Reference Samples
        IntraReferenceSamples        *intra_ref_ptr;

        // Transform and Quantization Buffers
        EbTransQuantBuffers          *trans_quant_buffers_ptr;
        struct EncDecContext         *enc_dec_context_ptr;

        uint64_t                       *fast_cost_array;
        uint64_t                       *full_cost_array;
        uint64_t                       *full_cost_skip_ptr;
        uint64_t                       *full_cost_merge_ptr;
        // Lambda
#if ADD_DELTA_QP_SUPPORT
        uint16_t                        qp;
#else
        uint8_t                         qp;
#endif
        uint8_t                         chroma_qp;
        uint32_t                        fast_lambda;
        uint32_t                        full_lambda;
        uint32_t                        fast_chroma_lambda;
        uint32_t                        full_chroma_lambda;
        uint32_t                        full_chroma_lambda_sao;

        //  Context Variables---------------------------------
        LargestCodingUnit            *sb_ptr;
        TransformUnit                *txb_ptr;
        CodingUnit                   *cu_ptr;
        const BlockGeom                *blk_geom;
        PredictionUnit               *pu_ptr;
        const PredictionUnitStats    *pu_stats;
        MvUnit                        mv_unit;

        // Entropy Coder
        EntropyCoder                 *coeff_est_entropy_coder_ptr;
        MdEncPassCuData               md_ep_pipe_sb[BLOCK_MAX_COUNT_SB_128];
#if !OPT_LOSSLESS_0
        uint8_t                         group_of8x8_blocks_count;
        uint8_t                         group_of16x16_blocks_count;
#endif
        uint8_t                         pu_itr;
        uint8_t                         cu_size_log2;
        uint8_t                         best_candidate_index_array[MAX_NFL + 2];
        uint8_t                         sorted_candidate_index_array[MAX_NFL];
        uint16_t                        cu_origin_x;
        uint16_t                        cu_origin_y;
#if !OPT_LOSSLESS_0
        uint64_t                        chroma_weight;
#endif
        uint8_t                         sb_sz;
        uint32_t                        sb_origin_x;
        uint32_t                        sb_origin_y;
        uint32_t                        round_origin_x;
        uint32_t                        round_origin_y;
        uint16_t                        pu_origin_x;
        uint16_t                        pu_origin_y;
        uint16_t                        pu_width;
        uint16_t                        pu_height;
        EbPfMode                        pf_md_mode;
#if !OPT_LOSSLESS_0
        unsigned                        luma_intra_ref_samples_gen_done      : 2; // only 1 bit is needed, but used two for rounding
        unsigned                        chroma_intra_ref_samples_gen_done    : 2; // only 1 bit is needed, but used two for rounding
        unsigned                        generate_mvp                         : 2; // only 1 bit is needed, but used two for rounding
#endif
        uint32_t                        full_recon_search_count;
        EbBool                          cu_use_ref_src_flag;
        uint16_t                        qp_index;
        uint64_t                        three_quad_energy;
#if SEARCH_UV_MODE
        EbBool                          uv_search_path;
        UvPredictionMode                best_uv_mode    [UV_PAETH_PRED + 1][(MAX_ANGLE_DELTA << 1) + 1];
        int32_t                         best_uv_angle   [UV_PAETH_PRED + 1][(MAX_ANGLE_DELTA << 1) + 1];
        uint64_t                        best_uv_cost    [UV_PAETH_PRED + 1][(MAX_ANGLE_DELTA << 1) + 1];
        uint64_t                        fast_luma_rate  [UV_PAETH_PRED + 1][(MAX_ANGLE_DELTA << 1) + 1];
        uint64_t                        fast_chroma_rate[UV_PAETH_PRED + 1][(MAX_ANGLE_DELTA << 1) + 1];
#endif
        // Needed for DC prediction
        EbBool                          is_left_availble;
        EbBool                          is_above_availble;
        int32_t                         is_inter_ctx;
        uint8_t                         intra_luma_left_mode;
        uint8_t                         intra_luma_top_mode;
        uint8_t                         intra_chroma_left_mode;
        uint8_t                         intra_chroma_top_mode;
        int16_t                         pred_buf_q3[CFL_BUF_SQUARE]; // Hsan: both MD and EP to use pred_buf_q3 (kept 1, and removed the 2nd)
#if MRP_DUPLICATION_FIX
        uint8_t                           injected_ref_type_l0_array[MODE_DECISION_CANDIDATE_MAX_COUNT]; // used to do not inject existing MV
        uint8_t                           injected_ref_type_l1_array[MODE_DECISION_CANDIDATE_MAX_COUNT]; // used to do not inject existing MV
        uint8_t                           injected_ref_type_bipred_array[MODE_DECISION_CANDIDATE_MAX_COUNT]; // used to do not inject existing MV
#endif
        int16_t                           injected_mv_x_l0_array[MODE_DECISION_CANDIDATE_MAX_COUNT]; // used to do not inject existing MV
        int16_t                           injected_mv_y_l0_array[MODE_DECISION_CANDIDATE_MAX_COUNT]; // used to do not inject existing MV
        uint8_t                           injected_mv_count_l0;
                                          
        int16_t                           injected_mv_x_l1_array[MODE_DECISION_CANDIDATE_MAX_COUNT]; // used to do not inject existing MV
        int16_t                           injected_mv_y_l1_array[MODE_DECISION_CANDIDATE_MAX_COUNT]; // used to do not inject existing MV
        uint8_t                           injected_mv_count_l1;
                                          
        int16_t                           injected_mv_x_bipred_l0_array[MODE_DECISION_CANDIDATE_MAX_COUNT]; // used to do not inject existing MV
        int16_t                           injected_mv_y_bipred_l0_array[MODE_DECISION_CANDIDATE_MAX_COUNT]; // used to do not inject existing MV
        int16_t                           injected_mv_x_bipred_l1_array[MODE_DECISION_CANDIDATE_MAX_COUNT]; // used to do not inject existing MV
        int16_t                           injected_mv_y_bipred_l1_array[MODE_DECISION_CANDIDATE_MAX_COUNT]; // used to do not inject existing MV
        uint8_t                           injected_mv_count_bipred;
        uint32_t                          fast_candidate_intra_count;
        uint32_t                          fast_candidate_inter_count;
#if MEMORY_FOOTPRINT_OPT_ME_MV   
        uint32_t                          me_block_offset;
#endif
        // Multi-modes signal(s) 
        uint8_t                           nfl_level;
        uint8_t                           skip_interpolation_search;
        uint8_t                           parent_sq_type[MAX_PARENT_SQ];
        uint8_t                           parent_sq_has_coeff[MAX_PARENT_SQ];
        uint8_t                           parent_sq_pred_mode[MAX_PARENT_SQ];
        uint8_t                           chroma_level;
        PART                              nsq_table[NSQ_TAB_SIZE];
        uint8_t                           decoupled_fast_loop_search_method; 
        uint8_t                           decouple_intra_inter_fast_loop;
        uint8_t                           full_loop_escape;
        uint8_t                           global_mv_injection;
#if M9_NEAR_INJECTION
        uint8_t                           near_mv_injection;
#endif
        uint8_t                           warped_motion_injection;
        uint8_t                           unipred3x3_injection;
        uint8_t                           bipred3x3_injection;
        uint8_t                           interpolation_filter_search_blk_size;
        uint8_t                           redundant_blk;
#if CFL_FIX
        uint8_t                            cfl_temp_luma_recon[128 * 128];
#endif
#if SPATIAL_SSE
        EbBool                            spatial_sse_full_loop;
#endif
#if M9_INTER_SRC_SRC_FAST_LOOP
        uint8_t                           inter_fast_loop_src_src;
#endif
#if  BLK_SKIP_DECISION
        EbBool                            blk_skip_decision;
#endif  
#if OPT_QUANT_COEFF
        EbBool                            trellis_quant_coeff_optimization;
#endif
    } ModeDecisionContext;

    typedef void(*EbAv1LambdaAssignFunc)(
        uint32_t                    *fast_lambda,
        uint32_t                    *full_lambda,
        uint32_t                    *fast_chroma_lambda,
        uint32_t                    *full_chroma_lambda,
        uint8_t                      bit_depth,
        uint16_t                     qp_index);

    typedef void(*EbLambdaAssignFunc)(
        uint32_t                    *fast_lambda,
        uint32_t                    *full_lambda,
        uint32_t                    *fast_chroma_lambda,
        uint32_t                    *full_chroma_lambda,
        uint32_t                    *full_chroma_lambda_sao,
        uint8_t                      qp_hierarchical_position,
        uint8_t                      qp,
        uint8_t                      chroma_qp);

    /**************************************
     * Extern Function Declarations
     **************************************/
    extern EbErrorType mode_decision_context_ctor(
        ModeDecisionContext      **context_dbl_ptr,
        EbColorFormat              color_format,
        EbFifo                    *mode_decision_configuration_input_fifo_ptr,
        EbFifo                    *mode_decision_output_fifo_ptr);

    extern void reset_mode_decision_neighbor_arrays(
        PictureControlSet *picture_control_set_ptr);

    extern void lambda_assign_low_delay(
        uint32_t                    *fast_lambda,
        uint32_t                    *full_lambda,
        uint32_t                    *fast_chroma_lambda,
        uint32_t                    *full_chroma_lambda,
        uint32_t                    *full_chroma_lambda_sao,
        uint8_t                      qp_hierarchical_position,
        uint8_t                      qp,
        uint8_t                      chroma_qp);

    extern void lambda_assign_random_access(
        uint32_t                    *fast_lambda,
        uint32_t                    *full_lambda,
        uint32_t                    *fast_chroma_lambda,
        uint32_t                    *full_chroma_lambda,
        uint32_t                    *full_chroma_lambda_sao,
        uint8_t                      qp_hierarchical_position,
        uint8_t                      qp,
        uint8_t                      chroma_qp);

    extern const EbLambdaAssignFunc     lambda_assignment_function_table[4];
    extern const EbAv1LambdaAssignFunc av1_lambda_assignment_function_table[4];

    // Table that converts 0-63 Q-range values passed in outside to the Qindex
    // range used internally.
    static const uint8_t quantizer_to_qindex[] = {
        0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48,
        52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100,
        104, 108, 112, 116, 120, 124, 128, 132, 136, 140, 144, 148, 152,
        156, 160, 164, 168, 172, 176, 180, 184, 188, 192, 196, 200, 204,
        208, 212, 216, 220, 224, 228, 232, 236, 240, 244, 249, 255};

    extern void reset_mode_decision(
        ModeDecisionContext   *context_ptr,
        PictureControlSet     *picture_control_set_ptr,
        SequenceControlSet    *sequence_control_set_ptr,
        uint32_t                 segment_index);

    extern void mode_decision_configure_lcu(
        ModeDecisionContext   *context_ptr,
        LargestCodingUnit     *sb_ptr,
        PictureControlSet     *picture_control_set_ptr,
        SequenceControlSet    *sequence_control_set_ptr,
        uint8_t                  picture_qp,
        uint8_t                  sb_qp);


    extern void cfl_rd_pick_alpha(
        PictureControlSet             *picture_control_set_ptr,
        ModeDecisionCandidateBuffer   *candidateBuffer,
        LargestCodingUnit             *sb_ptr,
        ModeDecisionContext           *context_ptr,
        EbPictureBufferDesc           *input_picture_ptr,
        uint32_t                         inputCbOriginIndex,
        uint32_t                         cuChromaOriginIndex,
        EbAsm                            asm_type);


#ifdef __cplusplus
}
#endif
#endif // EbModeDecisionProcess_h
