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
#if IMPROVED_BIPRED_INJECTION || IMPROVED_UNIPRED_INJECTION
#define MODE_DECISION_CANDIDATE_MAX_COUNT               113// 61 Intra & 18+2x8+2x8 Inter
#else
#define MODE_DECISION_CANDIDATE_MAX_COUNT               79//35//20 // 61 Intra & 18 Inter
#endif
#if INC_NFL12
#define MODE_DECISION_CANDIDATE_BUFFER_MAX_COUNT        (MAX_NFL*6) //up to 6 depths
#else
#define MODE_DECISION_CANDIDATE_BUFFER_MAX_COUNT        33
#endif

#define INDEPENDENT_INTRA_CHROMA_MODE_TOTAL_COUNT       4       // Planar, Vertical, Horizontal, DC

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
    typedef struct MdEncPassCuData_s
    {
        uint64_t                    skip_cost;
        uint64_t                    merge_cost;
        uint64_t                    chroma_distortion;
        uint64_t                    y_full_distortion[DIST_CALC_TOTAL];
        uint64_t                    y_coeff_bits;
        uint32_t                    y_has_coeff;
        uint64_t                    fast_luma_rate;
        uint16_t                    y_count_non_zero_coeffs[4];// Store nonzero CoeffNum, per TU. If one TU, stored in 0, otherwise 4 tus stored in 0 to 3

    } MdEncPassCuData_t;


    typedef struct MdCodingUnit_s
    {
        unsigned                    tested_cu_flag                  : 1;   //tells whether this CU is tested in MD.
        unsigned                    mdc_array_index                 : 7;
        unsigned                    count_non_zero_coeffs           : 11;
        unsigned                    top_neighbor_depth              : 3;
        unsigned                    left_neighbor_depth             : 3;
        unsigned                    top_neighbor_mode               : 2;
        unsigned                    left_neighbor_mode              : 2;
        unsigned                    full_distortion                 : 32;
        unsigned                    chroma_distortion               : 32;
        unsigned                    chroma_distortion_inter_depth   : 32;
        PARTITION_CONTEXT           left_neighbor_partition;
        PARTITION_CONTEXT           above_neighbor_partition;
        uint64_t                    cost;
        uint64_t                    cost_luma;
        CandidateMv ed_ref_mv_stack[MODE_CTX_REF_FRAMES][MAX_REF_MV_STACK_SIZE];//to be used in MD and EncDec

    } MdCodingUnit_t;


    typedef struct ModeDecisionContext_s
    {
        EbFifo_t                       *mode_decision_configuration_input_fifo_ptr;
        EbFifo_t                       *mode_decision_output_fifo_ptr;
        int16_t                        *transform_inner_array_ptr;

        ModeDecisionCandidate_t       **fast_candidate_ptr_array;
        ModeDecisionCandidate_t        *fast_candidate_array;
        ModeDecisionCandidateBuffer_t **candidate_buffer_ptr_array;
        MdRateEstimationContext_t      *md_rate_estimation_ptr;
        InterPredictionContext_t       *inter_prediction_context;
        MdCodingUnit_t                  md_local_cu_unit[BLOCK_MAX_COUNT];
        CodingUnit_t                    md_cu_arr_nsq[BLOCK_MAX_COUNT];

        NeighborArrayUnit_t            *intra_luma_mode_neighbor_array;
        NeighborArrayUnit_t            *intra_chroma_mode_neighbor_array;
        NeighborArrayUnit_t            *mv_neighbor_array;
        NeighborArrayUnit_t            *skip_flag_neighbor_array;
        NeighborArrayUnit_t            *mode_type_neighbor_array;
        NeighborArrayUnit_t            *leaf_depth_neighbor_array;
        NeighborArrayUnit_t            *luma_recon_neighbor_array;
        NeighborArrayUnit_t            *cb_recon_neighbor_array;
        NeighborArrayUnit_t            *cr_recon_neighbor_array;
        NeighborArrayUnit_t            *skip_coeff_neighbor_array;
        NeighborArrayUnit_t            *luma_dc_sign_level_coeff_neighbor_array; // Stored per 4x4. 8 bit: lower 6 bits (COEFF_CONTEXT_BITS), shows if there is at least one Coef. Top 2 bit store the sign of DC as follow: 0->0,1->-1,2-> 1
        NeighborArrayUnit_t            *cr_dc_sign_level_coeff_neighbor_array; // Stored per 4x4. 8 bit: lower 6 bits(COEFF_CONTEXT_BITS), shows if there is at least one Coef. Top 2 bit store the sign of DC as follow: 0->0,1->-1,2-> 1
        NeighborArrayUnit_t            *cb_dc_sign_level_coeff_neighbor_array; // Stored per 4x4. 8 bit: lower 6 bits(COEFF_CONTEXT_BITS), shows if there is at least one Coef. Top 2 bit store the sign of DC as follow: 0->0,1->-1,2-> 1
        NeighborArrayUnit_t            *inter_pred_dir_neighbor_array;
        NeighborArrayUnit_t            *ref_frame_type_neighbor_array;
        NeighborArrayUnit_t            *leaf_partition_neighbor_array;
        NeighborArrayUnit32_t          *interpolation_type_neighbor_array;

        // TMVP
        EbReferenceObject_t            *reference_object_write_ptr;

        // Intra Reference Samples
        IntraReferenceSamples_t        *intra_ref_ptr;

        // Transform and Quantization Buffers
        EbTransQuantBuffers_t          *trans_quant_buffers_ptr;
        struct EncDecContext_s         *enc_dec_context_ptr;

        uint64_t                       *fast_cost_array;
        uint64_t                       *full_cost_array;
        uint64_t                       *full_cost_skip_ptr;
        uint64_t                       *full_cost_merge_ptr;

        // Fast loop buffers
        uint8_t                         buffer_depth_index_start[MAX_LEVEL_COUNT];
        uint8_t                         buffer_depth_index_width[MAX_LEVEL_COUNT];

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
        LargestCodingUnit_t            *sb_ptr;
        TransformUnit_t                *txb_ptr;
        CodingUnit_t                   *cu_ptr;
        const BlockGeom                *blk_geom;
        PredictionUnit_t               *pu_ptr;
        const PredictionUnitStats_t    *pu_stats;
        MvUnit_t                        mv_unit;

        // Entropy Coder
        EntropyCoder_t                 *coeff_est_entropy_coder_ptr;
        MdEncPassCuData_t               md_ep_pipe_sb[BLOCK_MAX_COUNT];

        uint8_t                         group_of8x8_blocks_count;
        uint8_t                         group_of16x16_blocks_count;
        uint8_t                         pu_itr;
        uint8_t                         cu_size_log2;
        uint8_t                         best_candidate_index_array[MAX_NFL];
        uint16_t                        cu_origin_x;
        uint16_t                        cu_origin_y;
        uint64_t                        chroma_weight;
        uint32_t                        use_chroma_information_in_fast_loop;
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
        unsigned                        luma_intra_ref_samples_gen_done      : 2; // only 1 bit is needed, but used two for rounding
        unsigned                        chroma_intra_ref_samples_gen_done    : 2; // only 1 bit is needed, but used two for rounding
        unsigned                        generate_mvp                         : 2; // only 1 bit is needed, but used two for rounding
        unsigned                        round_mv_to_integer                  : 2; // only 1 bit is needed, but used two for rounding
        uint32_t                        full_recon_search_count;
        EbBool                          cu_use_ref_src_flag;
        uint16_t                        qp_index;
        uint64_t                        three_quad_energy;

        // Needed for DC prediction
        EbBool                          is_left_availble;
        EbBool                          is_above_availble;
        int32_t                         is_inter_ctx;
        uint8_t                         intra_luma_left_mode;
        uint8_t                         intra_luma_top_mode;
        uint8_t                         intra_chroma_left_mode;
        uint8_t                         intra_chroma_top_mode;
        int16_t                         pred_buf_q3[CFL_BUF_SQUARE];
#if INTRA_CORE_OPT
        DECLARE_ALIGNED(16, uint8_t, left_data[MAX_MB_PLANE][MAX_TX_SIZE * 2 + 32]);
        DECLARE_ALIGNED(16, uint8_t, above_data[MAX_MB_PLANE][MAX_TX_SIZE * 2 + 32]);
        BlockSize  scaled_chroma_bsize;
#endif


        // Multi-modes signal(s) 
        uint8_t                           nfl_level;
#if INTERPOLATION_SEARCH_LEVELS
        uint8_t                           skip_interpolation_search;
#endif

    } ModeDecisionContext_t;

    typedef void(*EB_AV1_LAMBDA_ASSIGN_FUNC)(
        uint32_t                    *fast_lambda,
        uint32_t                    *full_lambda,
        uint32_t                    *fast_chroma_lambda,
        uint32_t                    *full_chroma_lambda,
        uint8_t                      bit_depth,
        uint16_t                     qp_index);

    typedef void(*EB_LAMBDA_ASSIGN_FUNC)(
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
        ModeDecisionContext_t      **context_dbl_ptr,
        EbFifo_t                    *mode_decision_configuration_input_fifo_ptr,
        EbFifo_t                    *mode_decision_output_fifo_ptr);

    extern void reset_mode_decision_neighbor_arrays(
        PictureControlSet_t *picture_control_set_ptr);

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

    extern const EB_LAMBDA_ASSIGN_FUNC     lambda_assignment_function_table[4];
    extern const EB_AV1_LAMBDA_ASSIGN_FUNC av1_lambda_assignment_function_table[4];

    // Table that converts 0-63 Q-range values passed in outside to the Qindex
    // range used internally.
    static const uint8_t quantizer_to_qindex[] = {
        0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48,
        52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100,
        104, 108, 112, 116, 120, 124, 128, 132, 136, 140, 144, 148, 152,
        156, 160, 164, 168, 172, 176, 180, 184, 188, 192, 196, 200, 204,
        208, 212, 216, 220, 224, 228, 232, 236, 240, 244, 249, 255};

    extern void reset_mode_decision(
        ModeDecisionContext_t   *context_ptr,
        PictureControlSet_t     *picture_control_set_ptr,
        SequenceControlSet_t    *sequence_control_set_ptr,
        uint32_t                 segment_index);

    extern void ModeDecisionConfigureLcu(
        ModeDecisionContext_t   *context_ptr,
        LargestCodingUnit_t     *sb_ptr,
        PictureControlSet_t     *picture_control_set_ptr,
        SequenceControlSet_t    *sequence_control_set_ptr,
        uint8_t                  picture_qp,
        uint8_t                  sb_qp);

#ifdef __cplusplus
}
#endif
#endif // EbModeDecisionProcess_h