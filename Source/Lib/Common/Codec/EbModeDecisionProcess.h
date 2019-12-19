// clang-format off
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
#include "EbObject.h"

#ifdef __cplusplus
extern "C" {
#endif
    /**************************************
     * Defines
     **************************************/
#define MODE_DECISION_CANDIDATE_MAX_COUNT               1855
#define DEPTH_ONE_STEP   21
#define DEPTH_TWO_STEP    5
#define DEPTH_THREE_STEP  1

#if MULTI_PASS_PD
#define PREDICTIVE_ME_MAX_MVP_CANIDATES  4
#define PREDICTIVE_ME_DEVIATION_TH      50
#define FULL_PEL_REF_WINDOW_WIDTH        7
#define FULL_PEL_REF_WINDOW_HEIGHT       5
#define HALF_PEL_REF_WINDOW              3
#define QUARTER_PEL_REF_WINDOW           3
#if M0_OPT
#define FULL_PEL_REF_WINDOW_WIDTH_EXTENDED  15
#define FULL_PEL_REF_WINDOW_HEIGHT_EXTENDED 15
#endif
#define EIGHT_PEL_REF_WINDOW 3
#endif

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

#if PAL_SUP
    typedef struct {
        uint8_t best_palette_color_map[MAX_PALETTE_SQUARE];
        int kmeans_data_buf[2 * MAX_PALETTE_SQUARE];
    } PALETTE_BUFFER;
#endif
    typedef struct MdCodingUnit
    {
        unsigned                    tested_cu_flag                  : 1;   //tells whether this CU is tested in MD.
        unsigned                    mdc_array_index                 : 7;
        unsigned                    count_non_zero_coeffs           : 11;
        unsigned                    top_neighbor_depth              : 8;
        unsigned                    left_neighbor_depth             : 8;
        unsigned                    top_neighbor_mode               : 2;
        unsigned                    left_neighbor_mode              : 2;
        unsigned                    full_distortion                 : 32;
        unsigned                    chroma_distortion               : 32;
        unsigned                    chroma_distortion_inter_depth   : 32;
        PartitionContextType           left_neighbor_partition;
        PartitionContextType           above_neighbor_partition;
        uint64_t                    cost;
        CandidateMv ed_ref_mv_stack[MODE_CTX_REF_FRAMES][MAX_REF_MV_STACK_SIZE];//to be used in MD and EncDec
        uint8_t                     avail_blk_flag ;   //tells whether this CU is tested in MD and have a valid cu data
    } MdCodingUnit;

    typedef struct ModeDecisionContext
    {
        EbDctor                      dctor;
        EbFifo                       *mode_decision_configuration_input_fifo_ptr;
        EbFifo                       *mode_decision_output_fifo_ptr;
        int16_t                        *transform_inner_array_ptr;

        ModeDecisionCandidate       **fast_candidate_ptr_array;
        ModeDecisionCandidate        *fast_candidate_array;
        ModeDecisionCandidateBuffer **candidate_buffer_ptr_array;
        ModeDecisionCandidateBuffer  *scratch_candidate_buffer;
        MdRateEstimationContext      *md_rate_estimation_ptr;
        EbBool                        is_md_rate_estimation_ptr_owner;
        InterPredictionContext       *inter_prediction_context;
        MdCodingUnit                  *md_local_cu_unit;
        CodingUnit                    *md_cu_arr_nsq;

        NeighborArrayUnit            *intra_luma_mode_neighbor_array;
        NeighborArrayUnit            *intra_chroma_mode_neighbor_array;
        NeighborArrayUnit            *mv_neighbor_array;
        NeighborArrayUnit            *skip_flag_neighbor_array;
        NeighborArrayUnit            *mode_type_neighbor_array;
        NeighborArrayUnit            *leaf_depth_neighbor_array;
        NeighborArrayUnit            *luma_recon_neighbor_array;
        NeighborArrayUnit            *cb_recon_neighbor_array;
        NeighborArrayUnit            *cr_recon_neighbor_array;
        NeighborArrayUnit            *tx_search_luma_recon_neighbor_array;
        NeighborArrayUnit            *luma_recon_neighbor_array16bit;
        NeighborArrayUnit            *cb_recon_neighbor_array16bit;
        NeighborArrayUnit            *cr_recon_neighbor_array16bit;
        NeighborArrayUnit            *tx_search_luma_recon_neighbor_array16bit;
        NeighborArrayUnit            *skip_coeff_neighbor_array;
        NeighborArrayUnit            *luma_dc_sign_level_coeff_neighbor_array; // Stored per 4x4. 8 bit: lower 6 bits (COEFF_CONTEXT_BITS), shows if there is at least one Coef. Top 2 bit store the sign of DC as follow: 0->0,1->-1,2-> 1
        NeighborArrayUnit            *full_loop_luma_dc_sign_level_coeff_neighbor_array; // Stored per 4x4. 8 bit: lower 6 bits (COEFF_CONTEXT_BITS), shows if there is at least one Coef. Top 2 bit store the sign of DC as follow: 0->0,1->-1,2-> 1
        NeighborArrayUnit            *tx_search_luma_dc_sign_level_coeff_neighbor_array; // Stored per 4x4. 8 bit: lower 6 bits (COEFF_CONTEXT_BITS), shows if there is at least one Coef. Top 2 bit store the sign of DC as follow: 0->0,1->-1,2-> 1
        NeighborArrayUnit            *cr_dc_sign_level_coeff_neighbor_array; // Stored per 4x4. 8 bit: lower 6 bits(COEFF_CONTEXT_BITS), shows if there is at least one Coef. Top 2 bit store the sign of DC as follow: 0->0,1->-1,2-> 1
        NeighborArrayUnit            *cb_dc_sign_level_coeff_neighbor_array; // Stored per 4x4. 8 bit: lower 6 bits(COEFF_CONTEXT_BITS), shows if there is at least one Coef. Top 2 bit store the sign of DC as follow: 0->0,1->-1,2-> 1
        NeighborArrayUnit            *txfm_context_array;
        NeighborArrayUnit            *inter_pred_dir_neighbor_array;
        NeighborArrayUnit            *ref_frame_type_neighbor_array;
        NeighborArrayUnit            *leaf_partition_neighbor_array;
        NeighborArrayUnit32          *interpolation_type_neighbor_array;

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
        SuperBlock                   *sb_ptr;
        TransformUnit                *txb_ptr;
        CodingUnit                   *cu_ptr;
        const BlockGeom                *blk_geom;
        PredictionUnit               *pu_ptr;
        MvUnit                        mv_unit;
#if PAL_SUP
        PALETTE_BUFFER            palette_buffer;
        PaletteInfo              palette_cand_array[MAX_PAL_CAND];
#endif
        // Entropy Coder
        EntropyCoder                 *coeff_est_entropy_coder_ptr;
        MdEncPassCuData               *md_ep_pipe_sb;
        uint8_t                         pu_itr;
        uint8_t                         cu_size_log2;
        uint32_t                         best_candidate_index_array[MAX_NFL_BUFF];
        uint32_t                         sorted_candidate_index_array[MAX_NFL];
        uint16_t                        cu_origin_x;
        uint16_t                        cu_origin_y;
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
        uint8_t                         hbd_mode_decision;
        uint16_t                        qp_index;
        uint64_t                        three_quad_energy;
        uint32_t                        txb_1d_offset;
        EbBool                          uv_search_path;
        UvPredictionMode                best_uv_mode    [UV_PAETH_PRED + 1][(MAX_ANGLE_DELTA << 1) + 1];
        int32_t                         best_uv_angle   [UV_PAETH_PRED + 1][(MAX_ANGLE_DELTA << 1) + 1];
        uint64_t                        best_uv_cost    [UV_PAETH_PRED + 1][(MAX_ANGLE_DELTA << 1) + 1];
        uint64_t                        fast_luma_rate  [UV_PAETH_PRED + 1][(MAX_ANGLE_DELTA << 1) + 1];
        uint64_t                        fast_chroma_rate[UV_PAETH_PRED + 1][(MAX_ANGLE_DELTA << 1) + 1];
        // Needed for DC prediction
        EbBool                          is_left_availble;
        EbBool                          is_above_availble;
        int32_t                         is_inter_ctx;
        uint8_t                         intra_luma_left_mode;
        uint8_t                         intra_luma_top_mode;
        uint8_t                         intra_chroma_left_mode;
        uint8_t                         intra_chroma_top_mode;
        EB_ALIGN(64) int16_t            pred_buf_q3[CFL_BUF_SQUARE]; // Hsan: both MD and EP to use pred_buf_q3 (kept 1, and removed the 2nd)
        uint8_t                         injected_ref_type_l0_array[MODE_DECISION_CANDIDATE_MAX_COUNT]; // used to do not inject existing MV
        uint8_t                         injected_ref_type_l1_array[MODE_DECISION_CANDIDATE_MAX_COUNT]; // used to do not inject existing MV
        uint8_t                         injected_ref_type_bipred_array[MODE_DECISION_CANDIDATE_MAX_COUNT]; // used to do not inject existing MV
        int16_t                         injected_mv_x_l0_array[MODE_DECISION_CANDIDATE_MAX_COUNT]; // used to do not inject existing MV
        int16_t                         injected_mv_y_l0_array[MODE_DECISION_CANDIDATE_MAX_COUNT]; // used to do not inject existing MV
        uint8_t                         injected_mv_count_l0;

        int16_t                         injected_mv_x_l1_array[MODE_DECISION_CANDIDATE_MAX_COUNT]; // used to do not inject existing MV
        int16_t                         injected_mv_y_l1_array[MODE_DECISION_CANDIDATE_MAX_COUNT]; // used to do not inject existing MV
        uint8_t                         injected_mv_count_l1;

        int16_t                         injected_mv_x_bipred_l0_array[MODE_DECISION_CANDIDATE_MAX_COUNT]; // used to do not inject existing MV
        int16_t                         injected_mv_y_bipred_l0_array[MODE_DECISION_CANDIDATE_MAX_COUNT]; // used to do not inject existing MV
        int16_t                         injected_mv_x_bipred_l1_array[MODE_DECISION_CANDIDATE_MAX_COUNT]; // used to do not inject existing MV
        int16_t                         injected_mv_y_bipred_l1_array[MODE_DECISION_CANDIDATE_MAX_COUNT]; // used to do not inject existing MV
        uint8_t                         injected_mv_count_bipred;
        uint32_t                        fast_candidate_inter_count;
        uint32_t                        me_block_offset;
        uint8_t                         tx_depth;
        uint8_t                         txb_itr;
        uint32_t                        me_sb_addr;
        uint32_t                        geom_offset_x;
        uint32_t                        geom_offset_y;
        int16_t                         luma_txb_skip_context;
        int16_t                         luma_dc_sign_context;
        int16_t                         cb_txb_skip_context;
        int16_t                         cb_dc_sign_context;
        int16_t                         cr_txb_skip_context;
        int16_t                         cr_dc_sign_context;
        // Multi-modes signal(s)
        uint8_t                         parent_sq_type[MAX_PARENT_SQ];
        uint8_t                         parent_sq_has_coeff[MAX_PARENT_SQ];
        uint8_t                         parent_sq_pred_mode[MAX_PARENT_SQ];
        uint8_t                         chroma_level;
        PART                            nsq_table[NSQ_TAB_SIZE];
        uint8_t                         full_loop_escape;
        uint8_t                         global_mv_injection;
        uint8_t                         nx4_4xn_parent_mv_injection;
#if MULTI_PASS_PD
        uint8_t                         new_nearest_injection;
#endif
        uint8_t                         new_nearest_near_comb_injection;
        uint8_t                         warped_motion_injection;
        uint8_t                         unipred3x3_injection;
        uint8_t                         bipred3x3_injection;
        uint8_t                         predictive_me_level;
        uint8_t                         interpolation_filter_search_blk_size;
        uint8_t                         redundant_blk;
        uint8_t                         *cfl_temp_luma_recon;
        uint16_t                        *cfl_temp_luma_recon16bit;
        EbBool                          spatial_sse_full_loop;
        EbBool                          blk_skip_decision;
        EbBool                          trellis_quant_coeff_optimization;
        int16_t                         best_spatial_pred_mv[2][4][2];
        int8_t                          valid_refined_mv[2][4];
        EbPictureBufferDesc            *input_sample16bit_buffer;
        DECLARE_ALIGNED(16, uint8_t, pred0[2 * MAX_SB_SQUARE]);
        DECLARE_ALIGNED(16, uint8_t, pred1[2 * MAX_SB_SQUARE]);
        DECLARE_ALIGNED(32, int16_t, residual1[MAX_SB_SQUARE]);
        DECLARE_ALIGNED(32, int16_t, diff10[MAX_SB_SQUARE]);
    unsigned int prediction_mse ;
    EbBool      variance_ready;

    uint32_t                            cand_buff_indices[CAND_CLASS_TOTAL][MAX_NFL_BUFF];
    uint8_t                             md_staging_mode;
#if MULTI_PASS_PD
    uint8_t                             md_staging_count_level;
#endif
    uint8_t                             bypass_md_stage_1[CAND_CLASS_TOTAL];

    uint32_t                            md_stage_0_count[CAND_CLASS_TOTAL];
    uint32_t                            md_stage_1_count[CAND_CLASS_TOTAL];
    uint32_t                            md_stage_2_count[CAND_CLASS_TOTAL];

    uint32_t                            md_stage_1_total_count;
    uint32_t                            md_stage_2_total_count;

    uint8_t                             combine_class12; // 1:class1 and 2 are combined.

    CAND_CLASS                          target_class;

    // fast_loop_core signals
    EbBool                              md_staging_use_bilinear;
    EbBool                              md_staging_skip_interpolation_search;
    EbBool                              md_staging_skip_inter_chroma_pred;

    // full_loop_core signals
    EbBool                              md_staging_skip_full_pred; // 0: perform luma & chroma prediction + interpolation search, 2: nothing (use information from previous stages)
    EbBool                              md_staging_skip_atb;
    EbBool                              md_staging_tx_search; // 0: skip, 1: use ref cost, 2: no shortcuts
    EbBool                              md_staging_skip_full_chroma;
    EbBool                              md_staging_skip_rdoq;

#if II_COMP_FLAG
    DECLARE_ALIGNED(16, uint8_t,      intrapred_buf[INTERINTRA_MODES][2 * 32 * 32]); //MAX block size for inter intra is 32x32
#endif
    uint64_t                           *ref_best_cost_sq_table;
    uint32_t                           *ref_best_ref_sq_table;
    uint8_t                             edge_based_skip_angle_intra;
    EbBool                              coeff_based_skip_atb;
    uint8_t                             prune_ref_frame_for_rec_partitions;
    unsigned int                        source_variance; // input block variance
    unsigned int                        inter_inter_wedge_variance_th;
    uint64_t                            md_exit_th;
    uint64_t                            md_stage_1_cand_prune_th;
    uint64_t                            md_stage_1_class_prune_th;
    uint64_t                            md_stage_2_cand_prune_th;
    uint64_t                            md_stage_2_class_prune_th;
    DECLARE_ALIGNED(16, uint8_t, obmc_buff_0[2* 2 * MAX_MB_PLANE * MAX_SB_SQUARE]);
    DECLARE_ALIGNED(16, uint8_t, obmc_buff_1[2* 2 * MAX_MB_PLANE * MAX_SB_SQUARE]);
    DECLARE_ALIGNED(16, uint8_t, obmc_buff_0_8b[2 * MAX_MB_PLANE * MAX_SB_SQUARE]);
    DECLARE_ALIGNED(16, uint8_t, obmc_buff_1_8b[2 * MAX_MB_PLANE * MAX_SB_SQUARE]);
    DECLARE_ALIGNED(16, int32_t, wsrc_buf[MAX_SB_SQUARE]);
    DECLARE_ALIGNED(16, int32_t, mask_buf[MAX_SB_SQUARE]);
    unsigned int pred_sse[REF_FRAMES];
    uint8_t                            *above_txfm_context;
    uint8_t                            *left_txfm_context;
#if COMBINE_MDC_NSQ_TABLE
    PART best_nsq_sahpe1;
    PART best_nsq_sahpe2;
    PART best_nsq_sahpe3;
    PART best_nsq_sahpe4;
    PART best_nsq_sahpe5;
    PART best_nsq_sahpe6;
    PART best_nsq_sahpe7;
    PART best_nsq_sahpe8;
#endif
    // square cost weighting for deciding if a/b shapes could be skipped
    uint32_t sq_weight;

    // signal for enabling shortcut to skip search depths
    uint8_t enable_auto_max_partition;


#if MULTI_PASS_PD
    MD_COMP_TYPE compound_types_to_try;
    uint8_t best_me_cand_only_flag;
    uint8_t dc_cand_only_flag;
    EbBool disable_angle_z2_intra_flag;
    uint8_t full_cost_shut_fast_rate_flag;
    EbBool coeff_based_nsq_cand_reduction;
    uint8_t tx_search_level;
    uint64_t tx_weight;
    uint8_t tx_search_reduced_set;
    uint8_t interpolation_search_level;
    EbBool rdoq_quantize_fp;
    uint8_t md_atb_mode;
    uint8_t md_pic_obmc_mode;
    uint8_t md_enable_inter_intra;
    uint8_t md_filter_intra_mode;
    uint8_t md_max_ref_count;
    EbBool md_skip_mvp_generation;
    int16_t full_pel_ref_window_width_th;
    int16_t full_pel_ref_window_height_th;

    // Signal to control initial and final pass PD setting(s)
    PD_PASS pd_pass;
#endif
    } ModeDecisionContext;

    typedef void(*EbAv1LambdaAssignFunc)(
        uint32_t                    *fast_lambda,
        uint32_t                    *full_lambda,
        uint32_t                    *fast_chroma_lambda,
        uint32_t                    *full_chroma_lambda,
        uint8_t                      bit_depth,
        uint16_t                     qp_index,
        EbBool                       hbd_mode_decision);

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
        ModeDecisionContext       *context_ptr,
        EbColorFormat              color_format,
        EbFifo                    *mode_decision_configuration_input_fifo_ptr,
        EbFifo                    *mode_decision_output_fifo_ptr,
        uint8_t                    enable_hbd_mode_decision
#if PAL_SUP
        ,uint8_t                 cfg_palette
#endif
    );

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
        SequenceControlSet    *sequence_control_set_ptr,
        ModeDecisionContext   *context_ptr,
        PictureControlSet     *picture_control_set_ptr,
        uint32_t                 segment_index);

    extern void mode_decision_configure_lcu(
        ModeDecisionContext   *context_ptr,
        PictureControlSet     *picture_control_set_ptr,
        uint8_t                  sb_qp);

    extern void cfl_rd_pick_alpha(
        PictureControlSet             *picture_control_set_ptr,
        ModeDecisionCandidateBuffer   *candidate_buffer,
        SuperBlock                    *sb_ptr,
        ModeDecisionContext           *context_ptr,
        EbPictureBufferDesc           *input_picture_ptr,
        uint32_t                         inputCbOriginIndex,
        uint32_t                         cuChromaOriginIndex);

#ifdef __cplusplus
}
#endif
#endif // EbModeDecisionProcess_h
// clang-format on
