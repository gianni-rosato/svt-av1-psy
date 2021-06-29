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
#define MODE_DECISION_CANDIDATE_MAX_COUNT (MODE_DECISION_CANDIDATE_MAX_COUNT_Y + 84)
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
#if !CLN_MOVE_SKIP_MODE_CHECK
    uint64_t skip_cost;
    uint64_t merge_cost;
#endif
    uint64_t chroma_distortion;
} MdEncPassCuData;

typedef struct {
    uint8_t best_palette_color_map[MAX_PALETTE_SQUARE];
    int     kmeans_data_buf[2 * MAX_PALETTE_SQUARE];
} PALETTE_BUFFER;
typedef struct MdBlkStruct {
    unsigned             mdc_array_index : 7;
    unsigned             count_non_zero_coeffs : 12;
    unsigned             top_neighbor_depth : 8;
    unsigned             left_neighbor_depth : 8;
    unsigned             full_distortion : 32;
    uint64_t             rec_dist_per_quadrant[4];
    PartitionContextType left_neighbor_partition;
    PartitionContextType above_neighbor_partition;
    uint64_t             cost;
    uint64_t
                default_cost; // Similar to cost but does not get updated @ d1_non_square_block_decision() and d2_inter_depth_block_decision()
    CandidateMv ed_ref_mv_stack[MODE_CTX_REF_FRAMES]
                               [MAX_REF_MV_STACK_SIZE]; //to be used in MD and EncDec
#if !SS_CLN_MVP_TABLE
    IntMv    ref_mvs[MODE_CTX_REF_FRAMES][MAX_MV_REF_CANDIDATES]; //used only for nonCompound modes.
#endif
#if !CLN_ENC_DEC
    uint32_t best_d1_blk;
#endif
    uint8_t *neigh_left_recon[3]; //only for MD
    uint8_t *neigh_top_recon[3];
    uint16_t *neigh_left_recon_16bit[3];
    uint16_t *neigh_top_recon_16bit[3];
    int32_t quantized_dc[3][MAX_TXB_COUNT];

#if !CLN_MOVE_SKIP_MODE_CHECK
    uint8_t   skip_mode_allowed;
#endif
    // wm
    EbWarpedMotionParams wm_params_l0;
    EbWarpedMotionParams wm_params_l1;
#if !SS_OPT_MD
    // compound
    uint8_t                compound_idx;
    InterInterCompoundData interinter_comp;
#endif
    // txb
    uint8_t u_has_coeff[TRANSFORM_UNIT_MAX_COUNT];
    uint8_t v_has_coeff[TRANSFORM_UNIT_MAX_COUNT];
    uint8_t y_has_coeff[TRANSFORM_UNIT_MAX_COUNT];
} MdBlkStruct;

struct ModeDecisionCandidate;
struct ModeDecisionCandidateBuffer;
struct InterPredictionContext;

typedef struct RefResults {
#if SS_CLN_REF_PRUNE
    uint8_t  do_ref; // to process this ref or not
#else
    uint8_t  list_i; // list index of this ref
    uint8_t  ref_i; // ref list index of this ref
    uint32_t dist; // distortion
    uint8_t  do_ref; // to process this ref  or not
    EbBool   valid_ref;
#endif
} RefResults;
#if SS_CLN_REF_PRUNE
typedef struct PmeResults {
    uint32_t dist; // distortion
} PmeResults;
#endif
typedef enum InterCandGroup {
    // elementary-groups
    PA_ME_GROUP,
    UNI_3x3_GROUP,
    BI_3x3_GROUP,
    NRST_NEW_NEAR_GROUP,
    NRST_NEAR_GROUP,
    PRED_ME_GROUP,
    GLOBAL_GROUP,
    // complex-groups
    WARP_GROUP,
    OBMC_GROUP,
    INTER_INTRA_GROUP,
    COMP_DIST,
    COMP_DIFF,
    COMP_WEDGE,
    TOT_INTER_GROUP
} InterCandGroup;
typedef struct InterCompCtrls {
#if OPT_COMP_MODE_CHECK
    uint8_t tot_comp_types;                        // total compound types to test; 0: OFF, 1: AVG, 2: AVG/DIST, 3: AVG/DIST/DIFF/WEDGE, 4: AVG/DIST/DIFF/WEDGE
#else
    uint8_t allowed_comp_types
        [MD_COMP_TYPES]; // Compound types to inject; AVG/DIST/DIFF/WEDGE (if a comp type is disallowed here, it will
        // override distance-based settings)
#endif
    uint8_t do_me;                                   // if true, test all compound types for me
    uint8_t do_pme;                                  // if true, test all compound types for pme
    uint8_t do_nearest_nearest;                      // if true, test all compound types for nearest_nearest
    uint8_t do_near_near;                            // if true, test all compound types for near_near
    uint8_t do_nearest_near_new;                     // if true, test all compound types for nearest_near_new
    uint8_t do_3x3_bi;                               // if true, test all compound types for 3x3_bipred

    uint8_t pred0_to_pred1_mult;                     // multiplier to the pred0_to_pred1_sad; 0: no pred0_to_pred1_sad-based pruning, >= 1: towards more inter-inter compound
    uint8_t use_rate;                                // if true, use rate @ compound params derivation
} InterCompCtrls;
typedef struct InterIntraCompCtrls {
    uint8_t enabled;
} InterIntraCompCtrls;
typedef struct ObmcControls {
    uint8_t enabled;
    EbBool  max_blk_size_16x16; // if true, cap the max block size that OBMC can be used to 16x16
} ObmcControls;
typedef struct TxtControls {
    uint8_t enabled;

    uint8_t txt_group_inter_lt_16x16; // group to use when inter and tx block < 16x16
    uint8_t txt_group_inter_gt_eq_16x16; // group to use when inter and tx block >= 16x16

    uint8_t txt_group_intra_lt_16x16; // group to use when intra and tx block < 16x16
    uint8_t txt_group_intra_gt_eq_16x16; // group to use when intra and tx block >= 16x16
#if TUNE_NEW_TXT_LVLS
    uint32_t early_exit_dist_th;  // Per unit distortion TH; if best TX distortion is below the TH, skip remaining TX types (0 OFF, higher = more aggressive)
    uint32_t early_exit_coeff_th; // If best TX number of coeffs is less than the TH, skip remaining TX types (0 OFF, higher = more aggressive)
#endif
} TxtControls;
typedef struct TxsCycleRControls {
    uint8_t  enabled; // On/Off feature control
    uint16_t intra_th; // Threshold to bypass intra TXS <the higher th the higher speed>
    uint16_t inter_th; // Threshold to bypass inter TXS <the higher th the higher speed>
} TxsCycleRControls;

typedef struct NearCountCtrls {
    uint8_t enabled;

    uint8_t near_count;       // max # of near to consider
    uint8_t near_near_count;    // max # of near_near to consider
}NearCountCtrls;

typedef struct RefPruningControls {
    uint8_t enabled; // 0: OFF; 1: use inter to inter distortion deviation to derive best_refs
    uint32_t max_dev_to_best[TOT_INTER_GROUP];     // 0: OFF; 1: limit the injection to the best references based on distortion
    uint32_t ref_idx_2_offset;
    uint32_t ref_idx_3_offset;
    uint8_t closest_refs
        [TOT_INTER_GROUP]; // 0: OFF; 1: limit the injection to the closest references based on distance (LAST/BWD)
} RefPruningControls;
#if !OPT_REFACTOR_DEPTH_REFINEMENT_CTRLS
typedef struct DepthRefinementCtrls {
    uint8_t enabled;

    int64_t sub_to_current_th; // decrease towards a more agressive level
    int64_t parent_to_current_th; // decrease towards a more agressive level
    uint8_t up_to_2_depth;                        // when 1, a maximum of 2 depth per block (PRED+Parent or PRED+Sub), 0: no restriction(s)
    uint8_t
        use_pred_block_cost; // add an offset to sub_to_current_th and parent_to_current_th on the cost range of the predicted block; use default ths for high cost(s) and more aggressive TH(s) for low cost(s)
} DepthRefinementCtrls;
#endif
typedef struct DepthRemovalCtrls {
    uint8_t enabled;
    uint8_t disallow_below_64x64;  // remove 32x32 blocks and below based on the sb_64x64 (me_distortion, variance)
    uint8_t disallow_below_32x32;  // remove 16x16 blocks and below based on the sb_64x64 (me_distortion, variance)
    uint8_t disallow_below_16x16;  // remove 8x8 blocks and below based on the sb_64x64 (me_distortion, variance)
}DepthRemovalCtrls;
typedef struct DepthCtrls {
    int8_t s_depth; // start depth; 0: consider no parent blocks; else number of parent blocks to consider, specified as a negative number (e.g. -2 means consider 2 parents)
    int8_t e_depth; // end depth; 0: consider no child blocks; else number of child blocks to consider, specified as a positive number (e.g. 2 means consider 2 children)
}DepthCtrls;
#define MAX_RANGE_CNT 8
typedef struct InDepthBlockSkipCtrls {
    uint16_t base_weight;                      // 0: in-depth-block-skip OFF; 1: in-depth-block-skip ON
                                               // higher towards more aggressive level(s)
                                               // 0: the estimated cost for the next children is not taken into account and the action will be lossless compared to no in - depth - block - skip
                                               // 100 : the normalized cost of next children is assumed to be equal to the normalized cost of past children

    uint8_t  cost_band_based_modulation;       // whether to amplify the base_weight based on the cost range of the parent block or not
    uint16_t max_cost_multiplier;              // the max cost beyond which the base_weight is zeroed out
    uint8_t  max_band_cnt;                     // the number of band(s)
    uint16_t weight_per_band[MAX_RANGE_CNT];   // the weight per band

    uint8_t  child_cnt_based_modulation;       // whether to modulate based on the child count
    uint16_t cnt_based_weight[3];              // to specify the weight per child cnt

} InDepthBlockSkipCtrls;
typedef struct LowerDepthBlockSkipCtrls {
    uint8_t enabled;
    float min_distortion_cost_ratio; // the distortion-to-cost ratio under wich the quad_deviation_th is zeroed out (feature is disabled)
    float quad_deviation_th;         // do not perform sub_depth if std_deviation of the 4 quadrants src-to-rec dist is less than std_deviation_th
    uint8_t skip_all;                // whether to skip all or only next depth; 0: skip only next depth; 1: skip all lower depths
}LowerDepthBlockSkipCtrls;
#if OPT_REFACTOR_DEPTH_REFINEMENT_CTRLS
#define MAX_RANGE_CNT 8
typedef struct DepthRefinementCtrls {
    uint8_t enabled;
    int64_t parent_to_current_th;                 // maximum allowed parent-to-current cost deviation beyond which the previous depth will not be added to PRED
    int64_t sub_to_current_th;                    // maximum allowed sub-to-current cost deviation beyond which the next depth will not be added to PRED
    uint8_t up_to_2_depth;                        // when 1, a maximum of 2 depth per block (PRED+Parent or PRED+Sub), 0: no restriction(s)
    uint8_t  cost_band_based_modulation;          // whether to decrement parent_to_current_th and sub_to_current_th based on the cost range of the parent block or not
    uint16_t max_cost_multiplier;                 // the max cost beyond which the decrement is ignored
    uint8_t  max_band_cnt;                        // the number of band(s)
    int64_t  decrement_per_band[MAX_RANGE_CNT];   // the offset per band
} DepthRefinementCtrls;
#endif
#if FTR_SUBSAMPLE_RESIDUAL
typedef struct SubresCtrls {
    uint8_t step;                                 // Residual sub-sampling step (0:OFF)
    uint8_t odd_to_even_deviation_th;             // Set step to 0 if the deviation between: (1) the pred-to-src SAD of even rows and (2) the pred-to-src SAD of odd rows
                                                  // of the 1st 64x64 block @ mds3 of PD0 is higher than odd_to_even_deviation_th
                                                  // when 0, the detection is OFF
} SubresCtrls;
#endif
typedef struct PfCtrls {
    EB_TRANS_COEFF_SHAPE pf_shape;
} PfCtrls;
typedef struct MdNsqMotionSearchCtrls {
    uint8_t enabled; // 0: NSQ motion search @ MD OFF; 1: NSQ motion search @ MD ON
    uint8_t use_ssd; // 0: search using SAD; 1: search using SSD
    uint8_t full_pel_search_width; // Full Pel search area width
    uint8_t full_pel_search_height; // Full Pel search area height
#if FTR_USE_PSAD
    uint8_t enable_psad; //Enable pSad
#endif
} MdNsqMotionSearchCtrls;
typedef struct MdSqMotionSearchCtrls {
    uint8_t enabled; // 0: SQ motion search @ MD OFF; 1: SQ motion search @ MD ON
    uint8_t use_ssd; // 0: search using SAD; 1: search using SSD

    uint16_t
        pame_distortion_th; // TH for pa_me distortion to determine whether to search (distortion per pixel)

    uint8_t  sprs_lev0_enabled; // 0: OFF; 1: ON
    uint8_t  sprs_lev0_step; // Sparse search step
    uint16_t sprs_lev0_w; // Sparse search area width
    uint16_t sprs_lev0_h; // Sparse search area height
    uint16_t max_sprs_lev0_w; // Max Sparse search area width
    uint16_t max_sprs_lev0_h; // Max Sparse search area height
    int16_t  sprs_lev0_multiplier; // search area multiplier (is a % -- 100 is no scaling)

    uint8_t  sprs_lev1_enabled; // 0: OFF; 1: ON
    uint8_t  sprs_lev1_step; // Sparse search step
    uint16_t sprs_lev1_w; // Sparse search area width
    uint16_t sprs_lev1_h; // Sparse search area height
    uint16_t max_sprs_lev1_w; // Max Sparse search area width
    uint16_t max_sprs_lev1_h; // Max Sparse search area height
    int16_t  sprs_lev1_multiplier; // search area multiplier (is a % -- 100 is no scaling)

    uint8_t  sprs_lev2_enabled; // 0: OFF; 1: ON
    uint8_t  sprs_lev2_step; // Sparse search step
    uint16_t sprs_lev2_w; // Sparse search area width
    uint16_t sprs_lev2_h; // Sparse search area height
#if FTR_USE_PSAD
    uint8_t enable_psad; // Enable pSad
#endif
} MdSqMotionSearchCtrls;
typedef struct MdPmeCtrls {
    uint8_t enabled; // 0: PME search @ MD OFF; 1: PME search @ MD ON
    uint8_t use_ssd; // 0: search using SAD; 1: search using SSD
#if FTR_IMPROVE_PME
    int     early_check_mv_th_multiplier; // Do not perform PME search for blocks that have a valid ME_MV unless the ME_MV has a different direction than all MVP(s) and the ME_MV mag is higher than MV_TH=f(early_check_mv_th_multiplier)
#endif
    uint8_t full_pel_search_width; // Full Pel search area width
    uint8_t full_pel_search_height; // Full Pel search area height
    int     pre_fp_pme_to_me_cost_th; // If pre_fp_pme_to_me_cost higher than pre_fp_pme_to_me_cost_th then PME_MV = ME_MV and exit (decrease towards a faster level)
    int     pre_fp_pme_to_me_mv_th; // If pre_fp_pme_to_me_mv smaller than pre_fp_pme_to_me_mv_th then PME_MV = ME_MV and exit (increase towards a faster level)
    int     post_fp_pme_to_me_cost_th; // If post_fp_pme_to_me_cost higher than post_fp_pme_to_me_cost_th then PME_MV = ME_MV and exit (decrease towards a faster level)
    int     post_fp_pme_to_me_mv_th; // If post_fp_pme_to_me_mv smaller than post_fp_pme_to_me_mv_th then PME_MV = ME_MV and exit (increase towards a faster level)
#if TUNE_BLOCK_SIZE
    uint8_t modulate_pme_for_blk_size_res; // only turn pme on for 32x32 blocks or 1080p 64x64 blocks
#endif
#if FTR_USE_PSAD
    uint8_t enable_psad; // Enable pSad
#endif
} MdPmeCtrls;
typedef struct MdSubPelSearchCtrls {
    uint8_t enabled; // 0: subpel search @ MD OFF; 1: subpel search @ MD ON
    SUBPEL_SEARCH_TYPE
        subpel_search_type; // USE_8_TAPS | USE_4_TAPS | USE_2_TAPS | USE_2_TAPS_ORIG (not supported)
    int subpel_iters_per_step; // Maximum number of steps in logarithmic subpel search before giving up.
    uint8_t eight_pel_search_enabled; // 0: OFF; 1: ON
#if TUNE_CTR_QUARTER_PEL
    SUBPEL_FORCE_STOP forced_stop_lt_eq_32; // When to stop subpel search if sq_size less than or equal than 32
    SUBPEL_FORCE_STOP forced_stop_gt_32;    // When to stop subpel search.if sq_size is greater than 32
#endif
    SUBPEL_SEARCH_METHODS subpel_search_method;   // Subpel_search_method can only be subpel_tree which does a subpixel
                                                  // logarithmic search that keeps stepping at 1/2 pixel units until
                                                  // you stop getting a gain, and then goes on to 1/4 and repeats
                                                  // the same process. Along the way it skips many diagonals.
#if OPT_M11_SUBPEL
    int pred_variance_th;
    uint8_t abs_th_mult;
    int round_dev_th;
#endif
#if OPT11_SUBPEL
    uint8_t skip_diag_refinement;
#endif
#if OPT_SKIP_SUBPEL_ZZ
    uint8_t skip_zz_mv;                          // don't perform subpel for (0,0) MVs
#endif
} MdSubPelSearchCtrls;
typedef struct ParentSqCoeffAreaBasedCyclesReductionCtrls {
    EbBool enabled;

    uint8_t high_freq_band1_th;         // cutoff for the highest coeff-area band [0-100]
    uint8_t high_freq_band1_level;      // level of action to use if luma coeff-area of parent SQ is >= high_freq_band1_th
    uint8_t high_freq_band2_th;         // cutoff for the second high coeff-area band [0-100]; should be less than high_freq_band1_th
    uint8_t high_freq_band2_level;      // level of action to use if luma coeff-area of parent SQ is >= high_freq_band2_th
    uint8_t high_freq_band3_th;         // cutoff for the third high coeff-area band [0-100]; should be less than high_freq_band2_th
    uint8_t high_freq_band3_level;      // level of action to use if luma coeff-area of parent SQ is >= high_freq_band3_th

    uint8_t enable_zero_coeff_action;   // enable for whether to apply action when parent SQ has 0 luma coefficients
    uint8_t zero_coeff_action;          // level of action to use if parent SQ has 0 luma coeffs
    uint8_t enable_one_coeff_action;    // enable for whether to apply action when parent SQ has 1 luma coefficients
    uint8_t one_coeff_action;           // level of action to use if parent SQ has 1 luma coeff

    uint8_t low_freq_band1_th;          // cutoff for the lowest coeff-area band [0-100]; should be less than high_freq_band2_th
    uint8_t low_freq_band1_level;       // level of action to use if luma coeff-area of parent SQ is < low_freq_band1_th
    uint8_t low_freq_band2_th;          // cutoff for the lowest coeff-area band [0-100]; should be less than high_freq_band2_th and larger than low_freq_band1_th
    uint8_t low_freq_band2_level;       // level of action to use if luma coeff-area of parent SQ is < low_freq_band2_th
}ParentSqCoeffAreaBasedCyclesReductionCtrls;
#if FTR_RDOQ_ONLY_DCT_DCT

typedef struct RdoqCtrls {
    uint8_t enabled;

    uint8_t eob_fast_y_inter;    // 0: do not use eob_fast  for luma inter; 1:  use eob_fast  for luma inter
    uint8_t eob_fast_y_intra;    // 0: do not use eob_fast  for luma intra; 1:  use eob_fast  for luma intra
    uint8_t eob_fast_uv_inter;   // 0: do not use eob_fast  for chroma inter; 1:  use eob_fast  for chroma inter
    uint8_t eob_fast_uv_intra;   // 0: do not use eob_fast  for chroma intra; 1:  use eob_fast  for chroma intra
    uint8_t fp_q_y;              // 0: use default quant for luma; 1: use fp_quant for luma
    uint8_t fp_q_uv;             // 0: use default quant for chroma; 1: use fp_quant for chroma
    uint8_t satd_factor;         // Do not perform rdoq if the tx satd > satd_factor
    uint8_t early_exit_th;       // Do not perform rdoq based on an early skip/non-skip cost
    uint8_t skip_uv;             // [MD Only] 0: RDOQ for both Luma & Chroma, 1: RDOQ for only Luma
#if OPT_RDOQ_EOB
    uint8_t dct_dct_only;        // [MD Only] 0: RDOQ for All txt type(s), 1: RDOQ for only DCT_DCT
    uint8_t eob_th;              // eob_th beyond which RDOQ is shut
#else
    uint8_t dct_dct_only;        // [MD Only] 0: RDOQ for both Luma & Chroma, 1: RDOQ for only Luma
#endif
#if FTR_RDO_OPT
    uint8_t eob_fast_th;
#endif
} RdoqCtrls;
#else
typedef struct RdoqCtrls {
    uint8_t enabled;

    uint8_t
        eob_fast_l_inter; // 0: do not use eob_fast  for luma inter; 1:  use eob_fast  for luma inter
    uint8_t
        eob_fast_l_intra; // 0: do not use eob_fast  for luma intra; 1:  use eob_fast  for luma intra
    uint8_t
        eob_fast_c_inter; // 0: do not use eob_fast  for chroma inter; 1:  use eob_fast  for chroma inter
    uint8_t
        eob_fast_c_intra; // 0: do not use eob_fast  for chroma intra; 1:  use eob_fast  for chroma intra
    uint8_t fp_q_l; // 0: use default quant for luma; 1: use fp_quant for luma
    uint8_t fp_q_c; // 0: use default quant for chroma; 1: use fp_quant for chroma
    uint8_t satd_factor; // do not perform rdoq if the tx satd > satd_factor
    uint8_t
        early_exit_th; // do not perform rdoq based on an early skip/non-skip cost, threshold for early exit is 5
    uint8_t disallow_md_rdoq_uv;
    uint8_t md_satd_factor;
} RdoqCtrls;
#endif
typedef struct NicCtrls {
    uint8_t stage1_scaling_num; // Scaling numerator for post-stage 0 NICS: <x>/16
    uint8_t stage2_scaling_num; // Scaling numerator for post-stage 1 NICS: <x>/16
    uint8_t stage3_scaling_num; // Scaling numerator for post-stage 2 NICS: <x>/16
} NicCtrls;
typedef struct NicPruningCtrls {

    // class pruning signal(s)
    // mdsx_class_th (for class removal); reduce cand if deviation to the best_cand is higher than mdsx_cand_th

    // All bands (except the last) are derived as follows:
    // For band_index=0 to band_index=(mdsx_band_cnt-2),
    //     band=[band_index*band_width, (band_index+1)*band_width]; band_width = mdsx_class_th/(band_cnt-1)
    //     multiplier= 1 / ((band_index+1)*2)
    // Last band is [mds1_class_th, +?] = kill (nic=0)

    // e.g. mds1_class_th=20 and mds1_band_cnt=3
    // band_index  |0         |1        | 2       |
    // band        |0 to 10   |10 to 20 | 20 to +?|
    // action      |nic * 1   |nic * 1/2| nic *  0|

   // Post mds0
    uint64_t mds1_class_th;
    uint8_t  mds1_band_cnt; // >=2

    // Post mds1
    uint64_t mds2_class_th;
    uint8_t  mds2_band_cnt; // >=2

    // Post mds2
    uint64_t mds3_class_th;
    uint8_t  mds3_band_cnt; // >=2

    // cand pruning signal(s)
    // mdsx_cand_th (for single cand removal per class); remove cand if deviation to the best_cand for @ the target class is higher than mdsx_cand_th
    // mdsx_cand_th = base_th + sq_offset_th + intra_class_offset_th

    // Post mds0
    uint64_t mds1_cand_base_th;               // base_th


    // Post mds1
    uint64_t mds2_cand_base_th;

    // Post mds2
    uint64_t mds3_cand_base_th;

} NicPruningCtrls;
typedef struct CandEliminationCtlrs {
    uint32_t enabled;
    uint8_t dc_only;
    uint8_t inject_new_me;
    uint8_t inject_new_pme;
    uint8_t inject_new_warp;
#if OPT_EARLY_ELIM_TH
    uint8_t th_multiplier;      // factor to scale base TH by for distortion check
#endif
}CandEliminationCtlrs;
#if OPT_TXS_SEARCH
typedef struct TxsControls {
    uint8_t enabled;
    uint8_t prev_depth_coeff_exit; // Skip current depth if previous depth has no coeff
    uint8_t intra_class_max_depth; // Max number of depth(s) for INTRA classes
    uint8_t inter_class_max_depth; // Max number of depth(s) for INTER classes
    int depth1_txt_group_offset;   // Offset to be subtracted from default txt-group to derive the txt-group of depth-1
    int depth2_txt_group_offset;   // Offset to be subtracted from default txt-group to derive the txt-group of depth-2
} TxsControls;
#endif
#if FTR_NEW_WM_LVL
typedef struct WmCtrls {
    uint8_t enabled;
    uint8_t use_wm_for_mvp;        // allow/disallow MW for MVP candidates
    uint8_t num_new_mv_refinement; // [0-12] number of refinement positions around NEW_MVs to use
} WmCtrls;
#endif
#if CHROMA_CLEANUP
typedef struct UvCtrls {
    uint8_t enabled;
    uint8_t uv_mode; // Indicates the chroma search mode
                     // CHROMA_MODE_0: Full chroma search @ MD
                     // CHROMA_MODE_1: Fast chroma search @ MD - No independent chroma mode search
                     // CHROMA_MODE_2: Chroma blind @ MD
    uint8_t nd_uv_serach_mode; // Non-direct chroma search 0: pre chroma search is used, 1: chroma search at last md_stage is used
    uint64_t uv_intra_th; // Threshold to skip  Non-direct chroma search.
    uint64_t uv_cfl_th;// Threshold to skip clf.
} UvCtrls;
#endif
#if TUNE_BLOCK_SIZE
typedef struct InterpolationSearchCtrls {
    IfsLevel interpolation_search_level; // set the interpolation search level for md stage
    uint8_t quarter_pel_only; // IFS only on for quarter-pel
    uint8_t modulate_filter_per_resolution; // use regular+smooth filters for >=720p, and regular+sharp for <=480p
#if OPT_IFS
    uint8_t early_skip; // Early skip the cheking of the filter based on cost of the syntax
    uint8_t subsampled_distortion; // use subsambped distortion for block_height greater than 16
    uint8_t skip_sse_rd_model; // bypass rd model
#endif
} InterpolationSearchCtrls;
#endif
#if TUNE_BLOCK_SIZE
typedef struct SpatialSSECtrls {
    EbBool spatial_sse_full_loop_level; // enable spatial-sse for each superblock
    uint8_t blk_size_res_based_modulation; // turns off spatial sse for blocks that are not 64x64 and have a resolution >480p, speed feature
} SpatialSSECtrls;
#endif
#if FTR_SKIP_MDS1
typedef struct SkipMDS1Ctrls {
    EbBool enabled;                 // enable skipping MDS1 in PD1
    uint32_t use_mds3_shortcuts_th; // if (best_mds0_distortion/QP < TH) use shortcuts for candidate at MDS3; 0: OFF, higher: more aggressive
    uint32_t force_1_cand_th;       // if (best_mds0_distortion/QP < TH) consider only the best candidate after MDS0; 0: OFF, higher: more aggressive.  Should be safer than use_mds3_shortcuts_th
} SkipMDS1Ctrls;
#endif
#if REMOVE_CLOSE_MVS
typedef struct RedundantCandCtrls {
    int score_th;
    int mag_th;
} RedundantCandCtrls;
#endif
#if OPT_EARLY_CAND_ELIM
typedef struct EarlyCandElimCtrls {
    uint8_t enabled;
    uint16_t mds0_distortion_th; // TH used to compare candidate distortion to best cost; higher is safer
} EarlyCandElimCtrls;
#endif
#if OPT_USE_INTRA_NEIGHBORING
typedef struct UseNeighbouringModeCtrls {
    uint8_t enabled;
} UseNeighbouringModeCtrls;
#endif
#if LIGHT_PD1
typedef struct BlockLocation {
    uint32_t input_origin_index;        //luma   block location in picture
    uint32_t input_cb_origin_in_index;  //chroma block location in picture
    uint32_t blk_origin_index;          //luma   block location in SB
    uint32_t blk_chroma_origin_index;   //chroma block location in SB
} BlockLocation;
#endif
#if FTR_LPD1_DETECTOR
typedef struct Lpd1Ctrls {
    EbBool enabled;             // Whether light-PD1 can be used for an SB
    EbBool use_light_pd1;       // Whether light-PD1 is set to be used for an SB (the detector may change this)
    EbBool use_lpd1_detector;   // Whether to use a detector; if use_light_pd1 is set to 1, the detector will protect tough SBs;
                                // if use_light_pd1 is set to 0, the detector will select easy SBs to use light-PD1 for (not yet implemented)
    uint16_t cost_th_dist;      // Distortion value used in cost TH for detector
    uint16_t coeff_th;          // Num non-zero coeffs used in detector
#if OPT_USE_MVS_LPD1_DETECTOR
    uint16_t max_mv_length;     // Max MV length TH used in the detector: 0 - (0,0) MV only; (uint16_t)~0 means no MV check (higher is more aggressive)
#endif
#if OPT_LIGHT_PD1_USE_ME_DIST_VAR
    uint32_t me_8x8_cost_variance_th; // me_8x8_cost_variance_th beyond which the PD1 is used (instead of light-PD1)
#endif
} Lpd1Ctrls;
#endif
typedef struct ModeDecisionContext {
    EbDctor  dctor;
    EbFifo * mode_decision_configuration_input_fifo_ptr;
    EbFifo * mode_decision_output_fifo_ptr;
    int16_t *transform_inner_array_ptr;

    ModeDecisionCandidate **       fast_candidate_ptr_array;
    ModeDecisionCandidate *        fast_candidate_array;
    ModeDecisionCandidateBuffer ** candidate_buffer_ptr_array;
    ModeDecisionCandidateBuffer *  candidate_buffer_tx_depth_1;
    ModeDecisionCandidateBuffer *  candidate_buffer_tx_depth_2;
    MdRateEstimationContext *      md_rate_estimation_ptr;
    EbBool                         is_md_rate_estimation_ptr_owner;
    struct MdRateEstimationContext rate_est_table;
    InterPredictionContext *       inter_prediction_context;
    MdBlkStruct *                  md_local_blk_unit;
    BlkStruct *                    md_blk_arr_nsq;
    uint8_t *                      avail_blk_flag;
    uint8_t* tested_blk_flag; //tells whether this CU is tested in MD.
    uint8_t* do_not_process_blk;
    MdcSbData *                    mdc_sb_array;

    NeighborArrayUnit *intra_luma_mode_neighbor_array;
#if !OPT_NA_SKIP
    NeighborArrayUnit *skip_flag_neighbor_array;
#endif
    NeighborArrayUnit *mode_type_neighbor_array;

    NeighborArrayUnit *luma_recon_neighbor_array;
    NeighborArrayUnit *cb_recon_neighbor_array;
    NeighborArrayUnit *cr_recon_neighbor_array;
    NeighborArrayUnit *tx_search_luma_recon_neighbor_array;
    NeighborArrayUnit *luma_recon_neighbor_array16bit;
    NeighborArrayUnit *cb_recon_neighbor_array16bit;
    NeighborArrayUnit *cr_recon_neighbor_array16bit;
    NeighborArrayUnit *tx_search_luma_recon_neighbor_array16bit;
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
#if !OPT_NA_IFS
    NeighborArrayUnit *  ref_frame_type_neighbor_array;
#endif
    NeighborArrayUnit *  leaf_partition_neighbor_array;
#if !OPT_NA_IFS
    NeighborArrayUnit32 *interpolation_type_neighbor_array;
#endif
    // Transform and Quantization Buffers
    EbTransQuantBuffers * trans_quant_buffers_ptr;
    struct EncDecContext *enc_dec_context_ptr;

    uint64_t *fast_cost_array;
    uint64_t *full_cost_array;
#if !CLN_MOVE_SKIP_MODE_CHECK
    uint64_t *full_cost_skip_ptr;
    uint64_t *full_cost_merge_ptr;
#endif
    // Lambda
    uint32_t fast_lambda_md[2];
    uint32_t full_lambda_md[2];
    uint32_t full_sb_lambda_md
        [2]; // for the case of lambda modulation (blk_lambda_tuning), full_lambda_md/fast_lambda_md corresponds
        // to block lambda and full_sb_lambda_md is the full lambda per sb
    EbBool blk_lambda_tuning;
    //  Context Variables---------------------------------
    SuperBlock *     sb_ptr;
    BlkStruct *      blk_ptr;
    const BlockGeom *blk_geom;
    PredictionUnit * pu_ptr;
#if !OPT_NA_SKIP
    MvUnit           mv_unit;
#endif
    PALETTE_BUFFER   palette_buffer;
    PaletteInfo      palette_cand_array[MAX_PAL_CAND];
    // Entropy Coder
    MdEncPassCuData *md_ep_pipe_sb;

    uint8_t         sb64_sq_no4xn_geom;   //simple geometry 64x64SB, Sq only, no 4xN
    uint8_t          pu_itr;
    uint32_t         *best_candidate_index_array;
    uint16_t         blk_origin_x;
    uint16_t         blk_origin_y;
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
    uint8_t          qp_index;
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
    uint8_t              injected_mv_count_bipred;
    uint32_t             fast_candidate_inter_count;
    uint32_t             me_block_offset;
    uint32_t             me_cand_offset;
#if OPT_INTER_PRED // rename
    // Pointer to a scratch buffer used by CFL & IFS
    EbPictureBufferDesc *scratch_prediction_ptr;
#else
    EbPictureBufferDesc *cfl_temp_prediction_ptr;
#endif
#if !REFCTR_ADD_QUANT_COEFF_BUFF_MD
    EbPictureBufferDesc
        *    residual_quant_coeff_ptr; // One buffer for residual and quantized coefficient
#endif
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
#if !SS_OPT_MD
    uint8_t              parent_sq_type[MAX_PARENT_SQ];
    uint8_t              parent_sq_pred_mode[MAX_PARENT_SQ];
#endif
#if !CHROMA_CLEANUP
    uint8_t              chroma_level;
    uint8_t              chroma_at_last_md_stage;
    uint64_t             chroma_at_last_md_stage_intra_th;
    uint64_t             chroma_at_last_md_stage_cfl_th;
#endif
    uint8_t              global_mv_injection;
    uint8_t              new_nearest_injection;
    uint8_t              new_nearest_near_comb_injection;
#if FTR_NEW_WM_LVL
    WmCtrls              wm_ctrls;
#else
    uint8_t              warped_motion_injection;
#endif
#if CHROMA_CLEANUP
    UvCtrls              uv_ctrls;
#endif
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
    MdStage      md_stage;
    uint32_t     *cand_buff_indices[CAND_CLASS_TOTAL];
    uint8_t      md_staging_mode;
#if SS_OPT_MD
    uint8_t      bypass_md_stage_1;
    uint8_t      bypass_md_stage_2;
#else
    uint8_t      bypass_md_stage_1[CAND_CLASS_TOTAL];
    uint8_t      bypass_md_stage_2[CAND_CLASS_TOTAL];
#endif
    uint32_t     md_stage_0_count[CAND_CLASS_TOTAL];
    uint32_t     md_stage_1_count[CAND_CLASS_TOTAL];
    uint32_t     md_stage_2_count[CAND_CLASS_TOTAL];
    uint32_t     md_stage_3_count[CAND_CLASS_TOTAL];
    uint32_t     md_stage_1_total_count;
    uint32_t     md_stage_2_total_count;
    uint32_t     md_stage_3_total_count;
    uint32_t     md_stage_3_total_intra_count;
    uint64_t     best_intra_cost;
    uint64_t     best_inter_cost;
    CandClass    target_class;
#if FTR_SKIP_MDS1
    uint8_t perform_mds1;
    uint8_t use_tx_shortcuts_mds3;
    SkipMDS1Ctrls skip_mds1_ctrls;
#endif
    // fast_loop_core signals
    EbBool md_staging_skip_interpolation_search;
    EbBool md_staging_skip_chroma_pred;
    // full_loop_core signals
    EbBool
           md_staging_perform_inter_pred; // 0: perform luma & chroma prediction + interpolation search, 2: nothing (use information from previous stages)
    EbBool md_staging_tx_size_mode; // 0: Tx Size recon only, 1:Tx Size search and recon
    EbBool md_staging_txt_level;
    EbBool md_staging_skip_full_chroma;
    EbBool md_staging_skip_rdoq;
    EbBool md_staging_spatial_sse_full_loop_level;
    EbBool md_staging_perform_intra_chroma_pred;
    DECLARE_ALIGNED(
        16, uint8_t,
        intrapred_buf[INTERINTRA_MODES][2 * 32 * 32]); //MAX block size for inter intra is 32x32
    uint64_t *ref_best_cost_sq_table;
    uint32_t *ref_best_ref_sq_table;
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
    uint32_t max_part0_to_part1_dev;
    // signal for enabling shortcut to skip search depths
    uint8_t dc_cand_only_flag;
    EbBool  disable_angle_z2_intra_flag;
    uint8_t shut_skip_ctx_dc_sign_update;
    uint8_t shut_fast_rate; // use coeff rate and slipt flag rate only (no MVP derivation)
    uint8_t fast_coeff_est_level; // estimate the rate of the first (eob/N) coeff(s) and last coeff only
#if FTR_ENHANCE_FAST_COEFF_RATE
    uint8_t md_staging_fast_coeff_est_level; // Control fast_coeff_est_level per mds
#endif
#if FTR_SUBSAMPLE_RESIDUAL
    uint8_t md_staging_subres_step; // Control subres_step per mds
#endif
#if !TUne_TX_TYPE_LEVELS
    uint8_t tx_search_level;
#endif
#if !TUNE_BLOCK_SIZE
    uint8_t              interpolation_search_level;
#endif
#if !OPT_TXS_SEARCH
    uint8_t              md_tx_size_search_mode;
#endif
    uint8_t              md_pic_obmc_level;
    uint8_t              md_enable_paeth;
    uint8_t              md_enable_smooth;
    uint8_t              md_inter_intra_level;
    uint8_t              md_filter_intra_level;
    uint8_t              md_intra_angle_delta;
    uint8_t              md_allow_intrabc;
    uint8_t              md_palette_level;
    uint8_t              dist_based_ref_pruning;
    DepthRemovalCtrls    depth_removal_ctrls;
    DepthCtrls           depth_ctrls; // control which depths can be considered in PD1
    InDepthBlockSkipCtrls in_depth_block_skip_ctrls;
    LowerDepthBlockSkipCtrls lower_depth_block_skip_ctrls;
    DepthRefinementCtrls depth_refinement_ctrls;
    int64_t parent_to_current_deviation;
    int64_t child_to_current_deviation;
#if !OPT_REDUCE_TX
    uint8_t              pf_level;
#endif
#if FTR_SUBSAMPLE_RESIDUAL
    SubresCtrls          subres_ctrls;
    uint8_t              is_subres_safe;
#endif
    PfCtrls              pf_ctrls;
    // Control signals for MD sparse search (used for increasing ME search for active clips)
    uint8_t                md_sq_mv_search_level;
    MdSqMotionSearchCtrls  md_sq_me_ctrls;
    uint8_t                md_nsq_mv_search_level;
    MdNsqMotionSearchCtrls md_nsq_motion_search_ctrls;
    uint8_t                md_pme_level;
    MdPmeCtrls             md_pme_ctrls;
    uint8_t                md_subpel_me_level;
    MdSubPelSearchCtrls    md_subpel_me_ctrls;
    uint8_t                md_subpel_pme_level;
    MdSubPelSearchCtrls    md_subpel_pme_ctrls;
#if SS_CLN_REF_PRUNE
    PmeResults             pme_res[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];
#else
    RefResults             pme_res[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];
#endif
    ObmcControls           obmc_ctrls;
    InterCompCtrls         inter_comp_ctrls;
    InterIntraCompCtrls    inter_intra_comp_ctrls;
    RefResults ref_filtering_res[TOT_INTER_GROUP][MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];
    RefPruningControls ref_pruning_ctrls;
    // Signal to control initial and final pass PD setting(s)
    PdPass pd_pass;

    EbBool            md_disable_cfl;
#if OPT_TXS_SEARCH
    TxsControls       txs_ctrls;
#endif
    TxtControls       txt_ctrls;
    NearCountCtrls near_count_ctrls;
    RdoqCtrls         rdoq_ctrls;
    uint8_t           disallow_4x4;
    uint8_t           md_disallow_nsq;
    uint64_t          best_nsq_default_cost;
    uint64_t          default_cost_per_shape[NUMBER_OF_SHAPES];
    ParentSqCoeffAreaBasedCyclesReductionCtrls parent_sq_coeff_area_based_cycles_reduction_ctrls;
    uint8_t           sb_size;
    EbPictureBufferDesc *recon_coeff_ptr[TX_TYPES];
    EbPictureBufferDesc *recon_ptr[TX_TYPES];
#if REFCTR_ADD_QUANT_COEFF_BUFF_MD
    EbPictureBufferDesc *quant_coeff_ptr[TX_TYPES];
#endif

    uint8_t              skip_intra;
    EbPictureBufferDesc *temp_residual_ptr;
    EbPictureBufferDesc *temp_recon_ptr;
    // Array for all nearest/near MVs for a block for single ref case
    MV mvp_array[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH][MAX_MVP_CANIDATES];
    // Count of all nearest/near MVs for a block for single ref case
    int8_t mvp_count[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];
    // Start/end position for MD sparse search
    int16_t sprs_lev0_start_x;
    int16_t sprs_lev0_end_x;
    int16_t sprs_lev0_start_y;
    int16_t sprs_lev0_end_y;
#if !OPT_TXS_SEARCH
    uint8_t         md_staging_tx_size_level;
#endif
    NicCtrls        nic_ctrls;
    NicPruningCtrls nic_pruning_ctrls;
    uint8_t         inter_compound_mode;
    MV              ref_mv;
#if !OPT_INTER_PRED
    uint8_t         ifs_is_regular_last; // If regular is last performed interp_filters @ IFS
#endif
#if !SS_OPT_MD
    uint8_t         use_prev_mds_res;
#endif
    uint16_t        sb_index;
#if OPT_EARLY_CAND_ELIM
    EarlyCandElimCtrls early_cand_elimination_ctrls;
#else
    uint8_t         early_cand_elimination;
#endif
    uint64_t        mds0_best_cost;
    uint8_t         mds0_best_class;
    uint8_t reduce_last_md_stage_candidate;
    uint32_t mds0_best_idx;
    CandClass mds0_best_class_it;
    uint32_t mds1_best_idx;
    CandClass mds1_best_class_it;
#if SS_OPT_MDS0
    uint8_t mds0_dist_type;
#else
    uint8_t use_var_in_mds0;
#endif
#if !CLN_MISC_CLEANUP
    uint32_t md_me_cost[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];
#endif
    uint32_t md_me_dist;
    uint8_t inject_new_me;
    uint8_t inject_new_pme;
    uint8_t inject_new_warp;
    uint8_t merge_inter_classes;
    uint8_t bypass_tx_search_when_zcoef;
#if OPT_ESTIMATE_REF_BITS
    uint64_t estimate_ref_frames_num_bits[MODE_CTX_REF_FRAMES]; // [TOTAL_REFS_PER_FRAME + 1]
#else
    uint64_t estimate_ref_frames_num_bits[MODE_CTX_REF_FRAMES][2]; // [TOTAL_REFS_PER_FRAME + 1][is_compound]
#endif
    CandEliminationCtlrs cand_elimination_ctrs;
#if !TUNE_NEW_TXT_LVLS
    uint32_t early_txt_search_exit_level; // should be moved to txt_ctrls
#endif
    uint8_t ep_use_md_skip_decision;
    uint32_t max_nics ; // Maximum number of candidates MD can support
    uint32_t max_nics_uv ; // Maximum number of candidates MD can support
#if TUNE_BLOCK_SIZE
    InterpolationSearchCtrls ifs_ctrls;
#endif
#if FTR_BYPASS_ENCDEC
    EbBool bypass_encdec; // If enabled, will bypass EncDec and copy recon/quant coeffs from MD; only supported for 8bit
    EbBool pred_depth_only; // Indicates whether only pred depth refinement is used in PD1 - not yet supported
    uint16_t coded_area_sb;
    uint16_t coded_area_sb_uv;
#endif
#if LIGHT_PD0
    EbBool use_light_pd0; // Use light PD0 path.  Assumes one class, no NSQ, no 4x4, TXT off, TXS off, PME off, etc.
#endif
#if LIGHT_PD1
#if FTR_LPD1_DETECTOR
    Lpd1Ctrls lpd1_ctrls;
#else
    EbBool use_light_pd1; // Use light PD1 path.
#endif
#else
    uint8_t use_best_mds0;
#endif
#if TUNE_BLOCK_SIZE
    SpatialSSECtrls      spatial_sse_ctrls;
#endif

#if CLN_GEOM
    uint16_t   init_max_block_cnt;
#endif
#if CHROMA_CLEANUP
    uint8_t end_plane;
#endif
#if FTR_LOW_AC_COST_EST
    int masked_compound_used;
    int ctx_comp_group_idx;
    int comp_index_ctx;
    uint8_t use_low_precision_cost_estimation;
#endif
#if FTR_USE_PSAD
    uint8_t enable_psad; // Enable pSad
#endif
#if FTR_PD_EARLY_EXIT
#if SS_CLN_LIGHT_PD0_PATH
    uint32_t pd0_early_exit_th;
    uint32_t pd0_inter_depth_bias;
#else
    uint64_t pd0_early_exit_th;
    uint64_t pd0_inter_depth_bias;
#endif
#endif
#if OPT_SUBPEL && !OPT_SUBPEL_SKIP_TH
    uint8_t resolution;
#endif
#if REMOVE_CLOSE_MVS
    RedundantCandCtrls redundant_cand_ctrls;
#endif
#if FTR_FASTER_CFL
    uint8_t cfl_itr_th;
#endif
#if FTR_REDUCE_UNI_PRED
    uint8_t reduce_unipred_candidates;
#endif
#if OPT_USE_INTRA_NEIGHBORING
    UseNeighbouringModeCtrls use_neighbouring_mode_ctrls;
    uint8_t is_intra_bordered;
    uint8_t updated_enable_pme;
#endif
#if FTR_SKIP_COST_ZERO_COEFF
    uint16_t lpd1_zero_y_coeff_exit; // skip cost calc and chroma TX/compensation if there are zero luma coeffs
#endif
#if CHROMA_DETECTOR
    COMPONENT_TYPE lpd1_chroma_comp; // chroma components to compensate at MDS3 of LPD1
#endif
#if FIX_DO_NOT_TEST_CORRUPTED_MVS
    uint8_t corrupted_mv_check;
#endif
} ModeDecisionContext;

typedef void (*EbAv1LambdaAssignFunc)(PictureControlSet *pcs_ptr, uint32_t *fast_lambda,
                                      uint32_t *full_lambda, uint8_t bit_depth, uint16_t qp_index,
                                      EbBool multiply_lambda);

/**************************************
     * Extern Function Declarations
     **************************************/
extern EbErrorType mode_decision_context_ctor(ModeDecisionContext *context_ptr,
                                              EbColorFormat color_format, uint8_t sb_size,
                                              uint8_t enc_mode,
#if CLN_GEOM
                                              uint16_t max_block_cnt,
#endif
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

#if FTR_6L_QPS
#define FIXED_QP_OFFSET_COUNT 6
static const int percents[2][FIXED_QP_OFFSET_COUNT] =
{
    { 75, 70, 60, 20, 15, 0 },
    { 76, 60, 30, 15,  8, 4 } // libaom offsets
};
#endif
extern void reset_mode_decision(SequenceControlSet *scs_ptr, ModeDecisionContext *context_ptr,
                                PictureControlSet *pcs_ptr, uint16_t tile_row_idx,
                                uint32_t segment_index);

extern void mode_decision_configure_sb(ModeDecisionContext *context_ptr, PictureControlSet *pcs_ptr,
                                       uint8_t sb_qp);
extern void md_cfl_rd_pick_alpha(PictureControlSet *          pcs_ptr,
#if OPT_CHROMA_PATH
                                 ModeDecisionCandidateBuffer *candidate_buffer,
#else
                                 ModeDecisionCandidateBuffer *candidate_buffer, SuperBlock *sb_ptr,
#endif
                                 ModeDecisionContext *context_ptr,
                                 EbPictureBufferDesc *input_picture_ptr,
                                 uint32_t             input_cb_origin_in_index,
                                 uint32_t             blk_chroma_origin_index);

#ifdef __cplusplus
}
#endif
#endif // EbModeDecisionProcess_h
