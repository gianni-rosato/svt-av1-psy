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

#ifndef EbModeDecisionProcess_h
#define EbModeDecisionProcess_h

#include "EbDefinitions.h"
#include "EbModeDecision.h"
#include "EbSyntaxElements.h"
#include "EbSystemResourceManager.h"
#include "EbPictureBufferDesc.h"
#include "EbEntropyCoding.h"
#include "EbTransQuantBuffers.h"
#include "EbReferenceObject.h"
#include "EbNeighborArrays.h"
#include "EbObject.h"
#include "EbEncInterPrediction.h"

#ifdef __cplusplus
extern "C" {
#endif
/**************************************
     * Defines
     **************************************/
#define MODE_DECISION_CANDIDATE_MAX_COUNT_Y 1855
#define MODE_DECISION_CANDIDATE_MAX_COUNT \
    (MODE_DECISION_CANDIDATE_MAX_COUNT_Y + 84)
#define DEPTH_ONE_STEP 21
#define DEPTH_TWO_STEP 5
#define DEPTH_THREE_STEP 1
#define MAX_MVP_CANIDATES 4
/**************************************
      * Macros
      **************************************/

#define GROUP_OF_4_8x8_BLOCKS(origin_x, origin_y) \
    (((origin_x >> 3) & 0x1) && ((origin_y >> 3) & 0x1) ? EB_TRUE : EB_FALSE)
#define GROUP_OF_4_16x16_BLOCKS(origin_x, origin_y) \
    (((((origin_x >> 3) & 0x2) == 0x2) && (((origin_y >> 3) & 0x2) == 0x2)) ? EB_TRUE : EB_FALSE)
#define GROUP_OF_4_32x32_BLOCKS(origin_x, origin_y) \
    (((((origin_x >> 3) & 0x4) == 0x4) && (((origin_y >> 3) & 0x4) == 0x4)) ? EB_TRUE : EB_FALSE)

/**************************************
       * Coding Loop Context
       **************************************/
typedef struct MdEncPassCuData {
    uint64_t skip_cost;
    uint64_t merge_cost;
    uint64_t chroma_distortion;
} MdEncPassCuData;

typedef struct {
    uint8_t best_palette_color_map[MAX_PALETTE_SQUARE];
    int     kmeans_data_buf[2 * MAX_PALETTE_SQUARE];
} PALETTE_BUFFER;
typedef struct MdBlkStruct {
    unsigned             tested_blk_flag : 1; //tells whether this CU is tested in MD.
    unsigned             mdc_array_index : 7;
    unsigned             count_non_zero_coeffs : 11;
    unsigned             top_neighbor_depth : 8;
    unsigned             left_neighbor_depth : 8;
    unsigned             top_neighbor_mode : 2;
    unsigned             left_neighbor_mode : 2;
    unsigned             full_distortion : 32;
    unsigned             chroma_distortion : 32;
    unsigned             chroma_distortion_inter_depth : 32;
    PartitionContextType left_neighbor_partition;
    PartitionContextType above_neighbor_partition;
    uint64_t             cost;
    uint64_t             default_cost; // Similar to cost but does not get updated @ d1_non_square_block_decision() and d2_inter_depth_block_decision()
    CandidateMv          ed_ref_mv_stack[MODE_CTX_REF_FRAMES]
                               [MAX_REF_MV_STACK_SIZE]; //to be used in MD and EncDec
    uint8_t avail_blk_flag; //tells whether this CU is tested in MD and have a valid cu data
    IntMv ref_mvs[MODE_CTX_REF_FRAMES][MAX_MV_REF_CANDIDATES]; //used only for nonCompound modes.
    uint32_t best_d1_blk;
    uint8_t *neigh_left_recon[3]; //only for MD
    uint8_t *neigh_top_recon[3];
    uint16_t *neigh_left_recon_16bit[3];
    uint16_t *neigh_top_recon_16bit[3];
    uint8_t merge_flag;
    uint8_t sse_gradian_band[NUMBER_OF_SHAPES];
    int8_t pred_depth_refinement;
    int8_t pred_depth;
    // wm
    EbWarpedMotionParams wm_params_l0;
    EbWarpedMotionParams wm_params_l1;
    // ref frame
    int8_t ref_frame_index_l0;
    int8_t ref_frame_index_l1;
    // compound
    uint8_t compound_idx;
    InterInterCompoundData interinter_comp;
    // txb
    uint8_t u_has_coeff[TRANSFORM_UNIT_MAX_COUNT];
    uint8_t v_has_coeff[TRANSFORM_UNIT_MAX_COUNT];
    uint8_t y_has_coeff[TRANSFORM_UNIT_MAX_COUNT];
} MdBlkStruct;

struct ModeDecisionCandidate;
struct ModeDecisionCandidateBuffer;
struct InterPredictionContext;

typedef struct RefResults {
    uint8_t  list_i;   // list index of this ref
    uint8_t  ref_i;    // ref list index of this ref
    uint32_t dist;     // distortion
    uint8_t  do_ref;   // to process this ref  or not
    EbBool valid_ref;
} RefResults;
typedef enum InterCandGroup {
    PA_ME_GROUP         = 0,
    UNI_3x3_GROUP       = 1,
    BI_3x3_GROUP        = 2,
    NRST_NEW_NEAR_GROUP = 3,
    WARP_GROUP          = 4,
    NRST_NEAR_GROUP     = 5,
    PRED_ME_GROUP       = 6,
    GLOBAL_GROUP        = 7,
    TOT_INTER_GROUP     = 8,
} InterCandGroup;
typedef struct  ObmcControls {
    uint8_t enabled;
    uint8_t me_count;      //how many me candidates to consider injecting obmc
    uint8_t pme_best_ref;  //limit injection to best ref in pme
    uint8_t mvp_ref_count; //closest references allowed in mvp 0:4
    uint8_t near_count;    //how many near to consider injecting obmc 0..3
}ObmcControls;
typedef struct  AMdCycleRControls {
    uint8_t enabled; // On/Off feature control
    uint16_t skip_nsq_th;  // Threshold to bypass nsq <the higher th the higher speed>
    uint16_t switch_mode_th;
    uint8_t mode_offset;
}AMdCycleRControls;
typedef struct  TxtCycleRControls {
    uint8_t enabled;    // On/Off feature control
    uint16_t intra_th;  // Threshold to bypass intra TXT <the higher th the higher speed>
    uint16_t inter_th;  // Threshold to bypass inter TXT <the higher th the higher speed>
}TxtCycleRControls;
typedef struct  TxsCycleRControls {
    uint8_t enabled;    // On/Off feature control
    uint16_t intra_th;  // Threshold to bypass intra TXS <the higher th the higher speed>
    uint16_t inter_th;  // Threshold to bypass inter TXS <the higher th the higher speed>
}TxsCycleRControls;

typedef struct  InterCompoundControls {
    uint8_t enabled;
    uint8_t mrp_pruning_w_distortion;   // 0: OFF ; 1: Prune number of references based on ME/PME distortion
    uint8_t mrp_pruning_w_distance;     // 4: ALL ; 1: Prune number of references based on reference distance (Best 1)
    uint8_t wedge_search_mode;          // 0: Fast search: estimate Wedge sign 1: full search
}InterCompoundControls;

typedef struct RefPruningControls {
    uint8_t enabled; // 0: OFF; 1: use inter to inter distortion deviation to derive best_refs
    uint8_t best_refs[TOT_INTER_GROUP];     // 0: OFF; 1: limit the injection to the best references based on distortion
    uint8_t closest_refs[TOT_INTER_GROUP];  // 0: OFF; 1: limit the injection to the closest references based on distance (LAST/BWD)
}RefPruningControls;
typedef struct DepthRefinementCtrls {
    uint8_t enabled;

    int64_t sub_to_current_th; // decrease towards a more agressive level
    int64_t parent_to_current_th; // decrease towards a more agressive level

}DepthRefinementCtrls;
typedef struct MdNsqMotionSearchCtrls {
    uint8_t enabled;                    // 0: NSQ motion search @ MD OFF; 1: NSQ motion search @ MD ON
    uint8_t use_ssd;                    // 0: search using SAD; 1: search using SSD
    uint8_t full_pel_search_width;      // Full Pel search area width
    uint8_t full_pel_search_height;     // Full Pel search area height
}MdNsqMotionSearchCtrls;
typedef struct MdSqMotionSearchCtrls {
    uint8_t enabled;                    // 0: SQ motion search @ MD OFF; 1: SQ motion search @ MD ON
    uint8_t use_ssd;                    // 0: search using SAD; 1: search using SSD

    uint16_t pame_distortion_th;        // TH for pa_me distortion to determine whether to search (distortion per pixel)

    uint8_t  sprs_lev0_enabled;         // 0: OFF; 1: ON
    uint8_t  sprs_lev0_step;            // Sparse search step
    uint16_t sprs_lev0_w;               // Sparse search area width
    uint16_t sprs_lev0_h;               // Sparse search area height
    uint16_t max_sprs_lev0_w;           // Max Sparse search area width
    uint16_t max_sprs_lev0_h;           // Max Sparse search area height
    int16_t sprs_lev0_multiplier;       // search area multiplier (is a % -- 100 is no scaling)

    uint8_t  sprs_lev1_enabled;         // 0: OFF; 1: ON
    uint8_t  sprs_lev1_step;            // Sparse search step
    uint16_t sprs_lev1_w;               // Sparse search area width
    uint16_t sprs_lev1_h;               // Sparse search area height
    uint16_t max_sprs_lev1_w;           // Max Sparse search area width
    uint16_t max_sprs_lev1_h;           // Max Sparse search area height
    int16_t sprs_lev1_multiplier;       // search area multiplier (is a % -- 100 is no scaling)

    uint8_t  sprs_lev2_enabled;         // 0: OFF; 1: ON
    uint8_t  sprs_lev2_step;            // Sparse search step
    uint16_t sprs_lev2_w;               // Sparse search area width
    uint16_t sprs_lev2_h;               // Sparse search area height
}MdSqMotionSearchCtrls;
typedef struct MdPmeCtrls {
    uint8_t enabled;                    // 0: PME search @ MD OFF; 1: PME search @ MD ON
    uint8_t use_ssd;                    // 0: search using SAD; 1: search using SSD
    uint8_t full_pel_search_width;      // Full Pel search area width
    uint8_t full_pel_search_height;     // Full Pel search area height
    int pre_fp_pme_to_me_cost_th;   // If pre_fp_pme_to_me_cost higher than pre_fp_pme_to_me_cost_th then PME_MV = ME_MV and exit (decrease towards a faster level)
    int pre_fp_pme_to_me_mv_th;     // If pre_fp_pme_to_me_mv smaller than pre_fp_pme_to_me_mv_th then PME_MV = ME_MV and exit (increase towards a faster level)
    int post_fp_pme_to_me_cost_th;  // If post_fp_pme_to_me_cost higher than post_fp_pme_to_me_cost_th then PME_MV = ME_MV and exit (decrease towards a faster level)
    int post_fp_pme_to_me_mv_th;    // If post_fp_pme_to_me_mv smaller than post_fp_pme_to_me_mv_th then PME_MV = ME_MV and exit (increase towards a faster level)
}MdPmeCtrls;
typedef struct MdSubPelSearchCtrls {
    uint8_t enabled;                             // 0: subpel search @ MD OFF; 1: subpel search @ MD ON
    SUBPEL_SEARCH_TYPE subpel_search_type;       // USE_8_TAPS | USE_4_TAPS | USE_2_TAPS | USE_2_TAPS_ORIG (not supported)
    int subpel_iters_per_step;                   // Maximum number of steps in logarithmic subpel search before giving up.
    uint8_t eight_pel_search_enabled;            // 0: OFF; 1: ON
}MdSubPelSearchCtrls;
typedef struct CoeffBSwMdCtrls {
    uint8_t enabled;                // 0:  OFF; 1:  ON
    uint8_t mode_offset;            // Offset to the mode to switch to
    uint8_t skip_block;             // Allow skipping NSQ blocks
}CoeffBSwMdCtrls;
typedef struct ModeDecisionContext {
    EbDctor  dctor;
    EbFifo * mode_decision_configuration_input_fifo_ptr;
    EbFifo * mode_decision_output_fifo_ptr;
    int16_t *transform_inner_array_ptr;

    ModeDecisionCandidate **      fast_candidate_ptr_array;
    ModeDecisionCandidate *       fast_candidate_array;
    ModeDecisionCandidateBuffer **candidate_buffer_ptr_array;
    ModeDecisionCandidateBuffer *candidate_buffer_tx_depth_1;
    ModeDecisionCandidateBuffer *candidate_buffer_tx_depth_2;
    MdRateEstimationContext *     md_rate_estimation_ptr;
    EbBool                        is_md_rate_estimation_ptr_owner;
    struct MdRateEstimationContext rate_est_table;
    InterPredictionContext *      inter_prediction_context;
    MdBlkStruct *                md_local_blk_unit;
    BlkStruct *                  md_blk_arr_nsq;
    MdcSbData *mdc_sb_array;

    NeighborArrayUnit *intra_luma_mode_neighbor_array;
    NeighborArrayUnit *intra_chroma_mode_neighbor_array;
    NeighborArrayUnit *mv_neighbor_array;
    NeighborArrayUnit *skip_flag_neighbor_array;
    NeighborArrayUnit *mode_type_neighbor_array;
    NeighborArrayUnit *leaf_depth_neighbor_array;
    NeighborArrayUnit *luma_recon_neighbor_array;
    NeighborArrayUnit *cb_recon_neighbor_array;
    NeighborArrayUnit *cr_recon_neighbor_array;
    NeighborArrayUnit *tx_search_luma_recon_neighbor_array;
    NeighborArrayUnit *luma_recon_neighbor_array16bit;
    NeighborArrayUnit *cb_recon_neighbor_array16bit;
    NeighborArrayUnit *cr_recon_neighbor_array16bit;
    NeighborArrayUnit *tx_search_luma_recon_neighbor_array16bit;
    NeighborArrayUnit *skip_coeff_neighbor_array;
    NeighborArrayUnit *
        luma_dc_sign_level_coeff_neighbor_array; // Stored per 4x4. 8 bit: lower 6 bits (COEFF_CONTEXT_BITS), shows if there is at least one Coef. Top 2 bit store the sign of DC as follow: 0->0,1->-1,2-> 1
    NeighborArrayUnit *
        full_loop_luma_dc_sign_level_coeff_neighbor_array; // Stored per 4x4. 8 bit: lower 6 bits (COEFF_CONTEXT_BITS), shows if there is at least one Coef. Top 2 bit store the sign of DC as follow: 0->0,1->-1,2-> 1
    NeighborArrayUnit *
        tx_search_luma_dc_sign_level_coeff_neighbor_array; // Stored per 4x4. 8 bit: lower 6 bits (COEFF_CONTEXT_BITS), shows if there is at least one Coef. Top 2 bit store the sign of DC as follow: 0->0,1->-1,2-> 1
    NeighborArrayUnit *
        cr_dc_sign_level_coeff_neighbor_array; // Stored per 4x4. 8 bit: lower 6 bits(COEFF_CONTEXT_BITS), shows if there is at least one Coef. Top 2 bit store the sign of DC as follow: 0->0,1->-1,2-> 1
    NeighborArrayUnit *
                         cb_dc_sign_level_coeff_neighbor_array; // Stored per 4x4. 8 bit: lower 6 bits(COEFF_CONTEXT_BITS), shows if there is at least one Coef. Top 2 bit store the sign of DC as follow: 0->0,1->-1,2-> 1
    NeighborArrayUnit *  txfm_context_array;
    NeighborArrayUnit *  inter_pred_dir_neighbor_array;
    NeighborArrayUnit *  ref_frame_type_neighbor_array;
    NeighborArrayUnit *  leaf_partition_neighbor_array;
    NeighborArrayUnit32 *interpolation_type_neighbor_array;

    // Transform and Quantization Buffers
    EbTransQuantBuffers * trans_quant_buffers_ptr;
    struct EncDecContext *enc_dec_context_ptr;

    uint64_t *fast_cost_array;
    uint64_t *full_cost_array;
    uint64_t *full_cost_skip_ptr;
    uint64_t *full_cost_merge_ptr;
    // Lambda
    uint32_t fast_lambda_md[2];
    uint32_t full_lambda_md[2];
    uint32_t full_sb_lambda_md[2]; // for the case of lambda modulation (blk_lambda_tuning), full_lambda_md/fast_lambda_md corresponds
                                   // to block lambda and full_sb_lambda_md is the full lambda per sb
    EbBool       blk_lambda_tuning;
    //  Context Variables---------------------------------
    SuperBlock *     sb_ptr;
    BlkStruct *     blk_ptr;
    const BlockGeom *blk_geom;
    PredictionUnit * pu_ptr;
    MvUnit           mv_unit;
    PALETTE_BUFFER   palette_buffer;
    PaletteInfo      palette_cand_array[MAX_PAL_CAND];
    // Entropy Coder
    MdEncPassCuData *md_ep_pipe_sb;
    uint8_t          pu_itr;
    uint8_t          cu_size_log2;
    uint32_t         best_candidate_index_array[MAX_NFL_BUFF];
    uint32_t         sorted_candidate_index_array[MAX_NFL];
    uint16_t         blk_origin_x;
    uint16_t         blk_origin_y;
    uint8_t          sb_sz;
    uint32_t         sb_origin_x;
    uint32_t         sb_origin_y;
    uint32_t         round_origin_x;
    uint32_t         round_origin_y;
    uint16_t         pu_origin_x;
    uint16_t         pu_origin_y;
    uint16_t         pu_width;
    uint16_t         pu_height;
    EbPfMode         pf_md_mode;
    uint8_t          hbd_mode_decision;
    uint8_t         qp_index;
    uint64_t         three_quad_energy;
    uint32_t         txb_1d_offset;
    EbBool           uv_intra_comp_only;
    UvPredictionMode best_uv_mode[UV_PAETH_PRED + 1][(MAX_ANGLE_DELTA << 1) + 1];
    int32_t          best_uv_angle[UV_PAETH_PRED + 1][(MAX_ANGLE_DELTA << 1) + 1];
    uint64_t         best_uv_cost[UV_PAETH_PRED + 1][(MAX_ANGLE_DELTA << 1) + 1];
    uint64_t         fast_luma_rate[UV_PAETH_PRED + 1][(MAX_ANGLE_DELTA << 1) + 1];
    uint64_t         fast_chroma_rate[UV_PAETH_PRED + 1][(MAX_ANGLE_DELTA << 1) + 1];
    // Needed for DC prediction
    int32_t is_inter_ctx;
    uint8_t intra_luma_left_mode;
    uint8_t intra_luma_top_mode;
    uint8_t intra_chroma_left_mode;
    uint8_t intra_chroma_top_mode;
    EB_ALIGN(64)
    int16_t pred_buf_q3
        [CFL_BUF_SQUARE]; // Hsan: both MD and EP to use pred_buf_q3 (kept 1, and removed the 2nd)
    uint8_t injected_ref_type_l0_array
        [MODE_DECISION_CANDIDATE_MAX_COUNT]; // used to do not inject existing MV
    uint8_t injected_ref_type_l1_array
        [MODE_DECISION_CANDIDATE_MAX_COUNT]; // used to do not inject existing MV
    uint8_t injected_ref_type_bipred_array
        [MODE_DECISION_CANDIDATE_MAX_COUNT]; // used to do not inject existing MV
    int16_t injected_mv_x_l0_array
        [MODE_DECISION_CANDIDATE_MAX_COUNT]; // used to do not inject existing MV
    int16_t injected_mv_y_l0_array
        [MODE_DECISION_CANDIDATE_MAX_COUNT]; // used to do not inject existing MV
    uint8_t injected_mv_count_l0;

    int16_t injected_mv_x_l1_array
        [MODE_DECISION_CANDIDATE_MAX_COUNT]; // used to do not inject existing MV
    int16_t injected_mv_y_l1_array
        [MODE_DECISION_CANDIDATE_MAX_COUNT]; // used to do not inject existing MV
    uint8_t injected_mv_count_l1;

    int16_t injected_mv_x_bipred_l0_array
        [MODE_DECISION_CANDIDATE_MAX_COUNT]; // used to do not inject existing MV
    int16_t injected_mv_y_bipred_l0_array
        [MODE_DECISION_CANDIDATE_MAX_COUNT]; // used to do not inject existing MV
    int16_t injected_mv_x_bipred_l1_array
        [MODE_DECISION_CANDIDATE_MAX_COUNT]; // used to do not inject existing MV
    int16_t injected_mv_y_bipred_l1_array
        [MODE_DECISION_CANDIDATE_MAX_COUNT]; // used to do not inject existing MV
    uint8_t  injected_mv_count_bipred;
    uint32_t fast_candidate_inter_count;
    uint32_t me_block_offset;
    uint32_t me_cand_offset;
    EbPictureBufferDesc *cfl_temp_prediction_ptr;
    EbPictureBufferDesc *prediction_ptr_temp;
    EbPictureBufferDesc
        *residual_quant_coeff_ptr; // One buffer for residual and quantized coefficient
    uint8_t  tx_depth;
    uint8_t  txb_itr;
    uint32_t me_sb_addr;
    uint32_t geom_offset_x;
    uint32_t geom_offset_y;
    int16_t  luma_txb_skip_context;
    int16_t  luma_dc_sign_context;
    int16_t  cb_txb_skip_context;
    int16_t  cb_dc_sign_context;
    int16_t  cr_txb_skip_context;
    int16_t  cr_dc_sign_context;
    // Multi-modes signal(s)
    uint8_t              parent_sq_type[MAX_PARENT_SQ];
    uint8_t              parent_sq_has_coeff[MAX_PARENT_SQ];
    uint8_t              parent_sq_pred_mode[MAX_PARENT_SQ];
    uint8_t              chroma_level;
    uint8_t              chroma_at_last_md_stage;
    uint64_t             chroma_at_last_md_stage_intra_th;
    uint64_t             chroma_at_last_md_stage_cfl_th;
    uint8_t              global_mv_injection;
    uint8_t              new_nearest_injection;
    uint8_t              new_nearest_near_comb_injection;
    uint8_t              warped_motion_injection;
    uint8_t              unipred3x3_injection;
    uint8_t              bipred3x3_injection;
    uint8_t              redundant_blk;
    uint8_t              nic_level;
    uint8_t              similar_blk_avail;
    uint16_t             similar_blk_mds;
    uint8_t              inject_inter_candidates;
    uint8_t *            cfl_temp_luma_recon;
    uint16_t *           cfl_temp_luma_recon16bit;
    EbBool               spatial_sse_full_loop_level;
    EbBool               blk_skip_decision;
    int8_t               rdoq_level;
    int16_t              sb_me_mv[BLOCK_MAX_COUNT_SB_128][MAX_NUM_OF_REF_PIC_LIST][MAX_REF_IDX][2];
    MV                   fp_me_mv[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];
    MV                   sub_me_mv[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];
    uint32_t             post_subpel_me_mv_cost[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];
    int16_t              best_pme_mv[MAX_NUM_OF_REF_PIC_LIST][MAX_REF_IDX][2];
    int8_t               valid_pme_mv[MAX_NUM_OF_REF_PIC_LIST][MAX_REF_IDX];
    EbPictureBufferDesc *input_sample16bit_buffer;
    uint16_t             tile_index;
    DECLARE_ALIGNED(16, uint8_t, pred0[2 * MAX_SB_SQUARE]);
    DECLARE_ALIGNED(16, uint8_t, pred1[2 * MAX_SB_SQUARE]);
    DECLARE_ALIGNED(32, int16_t, residual1[MAX_SB_SQUARE]);
    DECLARE_ALIGNED(32, int16_t, diff10[MAX_SB_SQUARE]);
    unsigned int prediction_mse;
    MdStage md_stage;
    uint32_t     cand_buff_indices[CAND_CLASS_TOTAL][MAX_NFL_BUFF];
    uint8_t      md_staging_mode;
    uint8_t      md_staging_count_level;
    uint8_t      bypass_md_stage_1[CAND_CLASS_TOTAL];
    uint8_t bypass_md_stage_2[CAND_CLASS_TOTAL];
    uint32_t md_stage_0_count[CAND_CLASS_TOTAL];
    uint32_t md_stage_1_count[CAND_CLASS_TOTAL];
    uint32_t md_stage_2_count[CAND_CLASS_TOTAL];
    uint32_t md_stage_3_count[CAND_CLASS_TOTAL];
    uint32_t md_stage_1_total_count;
    uint32_t md_stage_2_total_count;
    uint32_t md_stage_3_total_count;
    uint32_t md_stage_3_total_intra_count;
    uint64_t best_intra_cost;
    uint64_t best_inter_cost;
    uint16_t skip_cfl_cost_dev_th;
    uint16_t mds3_intra_prune_th;
    CandClass target_class;

    // fast_loop_core signals
    EbBool md_staging_skip_interpolation_search;
    EbBool md_staging_skip_chroma_pred;
    // full_loop_core signals
    EbBool md_staging_perform_inter_pred; // 0: perform luma & chroma prediction + interpolation search, 2: nothing (use information from previous stages)
    EbBool md_staging_tx_size_mode; // 0: Tx Size recon only, 1:Tx Size search and recon
    EbBool md_staging_tx_search; // 0: skip, 1: use ref cost, 2: no shortcuts
    EbBool md_staging_skip_full_chroma;
    EbBool md_staging_skip_rdoq;
    EbBool md_staging_spatial_sse_full_loop_level;
    EbBool md_staging_perform_intra_chroma_pred;
    DECLARE_ALIGNED(
        16, uint8_t,
        intrapred_buf[INTERINTRA_MODES][2 * 32 * 32]); //MAX block size for inter intra is 32x32
    uint64_t *   ref_best_cost_sq_table;
    uint32_t *   ref_best_ref_sq_table;
    uint64_t     md_stage_1_cand_prune_th;
    uint64_t     md_stage_1_class_prune_th;
    uint64_t     md_stage_2_3_cand_prune_th;
    uint64_t     md_stage_2_3_class_prune_th;
    DECLARE_ALIGNED(16, uint8_t, obmc_buff_0[2 * 2 * MAX_MB_PLANE * MAX_SB_SQUARE]);
    DECLARE_ALIGNED(16, uint8_t, obmc_buff_1[2 * 2 * MAX_MB_PLANE * MAX_SB_SQUARE]);
    DECLARE_ALIGNED(16, uint8_t, obmc_buff_0_8b[2 * MAX_MB_PLANE * MAX_SB_SQUARE]);
    DECLARE_ALIGNED(16, uint8_t, obmc_buff_1_8b[2 * MAX_MB_PLANE * MAX_SB_SQUARE]);
    DECLARE_ALIGNED(16, int32_t, wsrc_buf[MAX_SB_SQUARE]);
    DECLARE_ALIGNED(16, int32_t, mask_buf[MAX_SB_SQUARE]);
    unsigned int pred_sse[REF_FRAMES];
    uint8_t *    above_txfm_context;
    uint8_t *    left_txfm_context;
    // square cost weighting for deciding if a/b shapes could be skipped
    uint32_t sq_weight;
    // signal for enabling shortcut to skip search depths
    MD_COMP_TYPE compound_types_to_try;
    uint8_t      dc_cand_only_flag;
    EbBool       disable_angle_z2_intra_flag;
    uint8_t      shut_fast_rate; // use coeff rate and slipt flag rate only (no MVP derivation)
    uint8_t      tx_search_level;
    uint8_t      interpolation_search_level;
    uint8_t      md_tx_size_search_mode;
    uint8_t      md_pic_obmc_level;
    uint8_t      md_enable_paeth;
    uint8_t      md_enable_smooth;
    uint8_t      md_inter_intra_level;
    uint8_t      md_filter_intra_level;
    uint8_t      md_intra_angle_delta;
    uint8_t      md_allow_intrabc;
    uint8_t      md_palette_level;
    uint8_t      dist_based_ref_pruning;
    uint8_t block_based_depth_refinement_level;
    DepthRefinementCtrls depth_refinement_ctrls;
    // Control signals for MD sparse search (used for increasing ME search for active clips)
    uint8_t md_sq_mv_search_level;
    MdSqMotionSearchCtrls md_sq_me_ctrls;
    uint8_t md_nsq_mv_search_level ;
    MdNsqMotionSearchCtrls md_nsq_motion_search_ctrls;
    uint8_t md_pme_level;
    MdPmeCtrls md_pme_ctrls;
    uint8_t md_subpel_me_level;
    MdSubPelSearchCtrls md_subpel_me_ctrls;
    uint8_t md_subpel_pme_level;
    MdSubPelSearchCtrls md_subpel_pme_ctrls;
    uint8_t      md_max_ref_count;
    RefResults    pme_res[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];
    ObmcControls obmc_ctrls;
    InterCompoundControls inter_comp_ctrls;
    RefResults ref_filtering_res[TOT_INTER_GROUP][MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];
    RefPruningControls ref_pruning_ctrls;
    // Signal to control initial and final pass PD setting(s)
    PdPass pd_pass;

    EbBool        md_disable_cfl;
    TxtCycleRControls txt_cycles_red_ctrls;
    TxsCycleRControls txs_cycles_red_ctrls;
    AMdCycleRControls admd_cycles_red_ctrls;
    uint8_t disallow_4x4;
    uint8_t       md_disallow_nsq;
    uint64_t best_nsq_default_cost;
    uint64_t default_cost_per_shape[NUMBER_OF_SHAPES];
    uint8_t enable_area_based_cycles_allocation;
    uint8_t sb_class;
    uint16_t coeff_area_based_bypass_nsq_th;
    uint8_t sb_size;

    EbPictureBufferDesc *recon_coeff_ptr[TX_TYPES];
    EbPictureBufferDesc *recon_ptr[TX_TYPES];
    uint32_t part_cnt[NUMBER_OF_SHAPES-1][FB_NUM][SSEG_NUM];
    uint16_t part_prob[NUMBER_OF_SHAPES-1][FB_NUM][SSEG_NUM];
    uint32_t pred_depth_count[DEPTH_DELTA_NUM][NUMBER_OF_SHAPES-1];
    uint32_t depth_prob[DEPTH_DELTA_NUM];
    uint32_t ad_md_prob[DEPTH_DELTA_NUM][NUMBER_OF_SHAPES-1];
    uint32_t txt_cnt[TXT_DEPTH_DELTA_NUM][TX_TYPES];
    uint32_t txt_prob[TXT_DEPTH_DELTA_NUM][TX_TYPES];
    uint8_t skip_intra;
    EbPictureBufferDesc* temp_residual_ptr;
    EbPictureBufferDesc* temp_recon_ptr;
    // Array for all nearest/near MVs for a block for single ref case
    MV mvp_array[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH][MAX_MVP_CANIDATES];
    // Count of all nearest/near MVs for a block for single ref case
    int8_t  mvp_count[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];
    // Start/end position for MD sparse search
    int16_t sprs_lev0_start_x;
    int16_t sprs_lev0_end_x;
    int16_t sprs_lev0_start_y;
    int16_t sprs_lev0_end_y;
    uint8_t txs_in_inter_classes;
    uint8_t nic_scaling_level;
    uint8_t inter_compound_mode;
    uint8_t switch_md_mode_based_on_sq_coeff;
    CoeffBSwMdCtrls cb_sw_md_ctrls;
    MV ref_mv;
} ModeDecisionContext;

typedef void (*EbAv1LambdaAssignFunc)(PictureControlSet* pcs_ptr, uint32_t *fast_lambda, uint32_t *full_lambda,
                                      uint8_t bit_depth, uint16_t qp_index,
                                      EbBool multiply_lambda);

/**************************************
     * Extern Function Declarations
     **************************************/
extern EbErrorType mode_decision_context_ctor(ModeDecisionContext *context_ptr,
                                              EbColorFormat        color_format,
                                              uint8_t sb_size,
                                              EbFifo *mode_decision_configuration_input_fifo_ptr,
                                              EbFifo *mode_decision_output_fifo_ptr,
                                              uint8_t enable_hbd_mode_decision,
                                              uint8_t cfg_palette);

extern const EbAv1LambdaAssignFunc av1_lambda_assignment_function_table[4];

// Table that converts 0-63 Q-range values passed in outside to the Qindex
// range used internally.
static const uint8_t quantizer_to_qindex[] = {
    0,   4,   8,   12,  16,  20,  24,  28,  32,  36,  40,  44,  48,  52,  56,  60,
    64,  68,  72,  76,  80,  84,  88,  92,  96,  100, 104, 108, 112, 116, 120, 124,
    128, 132, 136, 140, 144, 148, 152, 156, 160, 164, 168, 172, 176, 180, 184, 188,
    192, 196, 200, 204, 208, 212, 216, 220, 224, 228, 232, 236, 240, 244, 249, 255};

extern void reset_mode_decision(SequenceControlSet *scs_ptr, ModeDecisionContext *context_ptr,
                                PictureControlSet *pcs_ptr, uint16_t tile_row_idx, uint32_t segment_index);

extern void mode_decision_configure_sb(ModeDecisionContext *context_ptr, PictureControlSet *pcs_ptr,
                                       uint8_t sb_qp);
extern void md_cfl_rd_pick_alpha(PictureControlSet *          pcs_ptr,
                              ModeDecisionCandidateBuffer *candidate_buffer, SuperBlock *sb_ptr,
                              ModeDecisionContext *context_ptr,
                              EbPictureBufferDesc *input_picture_ptr,
                              uint32_t input_cb_origin_in_index, uint32_t blk_chroma_origin_index);

#ifdef __cplusplus
}
#endif
#endif // EbModeDecisionProcess_h
