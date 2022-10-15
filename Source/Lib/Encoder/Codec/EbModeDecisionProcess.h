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
    (((origin_x >> 3) & 0x1) && ((origin_y >> 3) & 0x1) ? TRUE : FALSE)
#define GROUP_OF_4_16x16_BLOCKS(origin_x, origin_y) \
    (((((origin_x >> 3) & 0x2) == 0x2) && (((origin_y >> 3) & 0x2) == 0x2)) ? TRUE : FALSE)
#define GROUP_OF_4_32x32_BLOCKS(origin_x, origin_y) \
    (((((origin_x >> 3) & 0x4) == 0x4) && (((origin_y >> 3) & 0x4) == 0x4)) ? TRUE : FALSE)

/**************************************
       * Coding Loop Context
       **************************************/
typedef struct MdEncPassCuData {
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
    uint8_t  *neigh_left_recon[3]; //only for MD
    uint8_t  *neigh_top_recon[3];
    uint16_t *neigh_left_recon_16bit[3];
    uint16_t *neigh_top_recon_16bit[3];
    int32_t   quantized_dc[3][MAX_TXB_COUNT];

    // wm
    EbWarpedMotionParams wm_params_l0;
    EbWarpedMotionParams wm_params_l1;
    // txb
    uint8_t u_has_coeff[TRANSFORM_UNIT_MAX_COUNT];
    uint8_t v_has_coeff[TRANSFORM_UNIT_MAX_COUNT];
    uint8_t y_has_coeff[TRANSFORM_UNIT_MAX_COUNT];
} MdBlkStruct;

struct ModeDecisionCandidate;
struct ModeDecisionCandidateBuffer;
typedef struct RefResults {
    uint8_t do_ref; // to process this ref or not
} RefResults;
typedef struct PmeResults {
    uint32_t dist; // distortion
} PmeResults;
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
    uint8_t
            tot_comp_types; // total compound types to test; 0: OFF, 1: AVG, 2: AVG/DIST, 3: AVG/DIST/DIFF/WEDGE, 4: AVG/DIST/DIFF/WEDGE
    uint8_t do_me; // if true, test all compound types for me
    uint8_t do_pme; // if true, test all compound types for pme
    uint8_t do_nearest_nearest; // if true, test all compound types for nearest_nearest
    uint8_t do_near_near; // if true, test all compound types for near_near
    uint8_t do_nearest_near_new; // if true, test all compound types for nearest_near_new
    uint8_t do_3x3_bi; // if true, test all compound types for 3x3_bipred
    uint8_t
        skip_mvp_on_ref_info; // Skip MVP compound based on ref frame type and neighbour ref frame types
    uint8_t
            pred0_to_pred1_mult; // multiplier to the pred0_to_pred1_sad; 0: no pred0_to_pred1_sad-based pruning, >= 1: towards more inter-inter compound
    uint8_t use_rate; // if true, use rate @ compound params derivation
} InterCompCtrls;
typedef struct InterIntraCompCtrls {
    uint8_t enabled;
} InterIntraCompCtrls;
typedef struct ObmcControls {
    uint8_t enabled;
    Bool    max_blk_size_to_refine_16x16; // if true, cap the max block size to refine to 16x16
    Bool    max_blk_size_16x16; // if true, cap the max block size to test to 16x16
} ObmcControls;
typedef struct TxtControls {
    uint8_t enabled;

    uint8_t txt_group_inter_lt_16x16; // group to use when inter and tx block < 16x16
    uint8_t txt_group_inter_gt_eq_16x16; // group to use when inter and tx block >= 16x16

    uint8_t txt_group_intra_lt_16x16; // group to use when intra and tx block < 16x16
    uint8_t txt_group_intra_gt_eq_16x16; // group to use when intra and tx block >= 16x16
    uint32_t
        early_exit_dist_th; // Per unit distortion TH; if best TX distortion is below the TH, skip remaining TX types (0 OFF, higher = more aggressive)
    uint32_t
        early_exit_coeff_th; // If best TX number of coeffs is less than the TH, skip remaining TX types (0 OFF, higher = more aggressive)
    uint16_t
        satd_early_exit_th_intra; // If dev. of satd of current tx_type to best is greater than TH, exit early from testing current tx_type. 0: off; larger value is safer
    uint16_t
        satd_early_exit_th_inter; // If dev. of satd of current tx_type to best is greater than TH, exit early from testing current tx_type. 0: off; larger value is safer
} TxtControls;
typedef struct TxsCycleRControls {
    uint8_t  enabled; // On/Off feature control
    uint16_t intra_th; // Threshold to bypass intra TXS <the higher th the higher speed>
    uint16_t inter_th; // Threshold to bypass inter TXS <the higher th the higher speed>
} TxsCycleRControls;

typedef struct NearCountCtrls {
    uint8_t enabled;

    uint8_t near_count; // max # of near to consider
    uint8_t near_near_count; // max # of near_near to consider
} NearCountCtrls;

typedef struct RefPruningControls {
    uint8_t  enabled; // 0: OFF; 1: use inter to inter distortion deviation to derive best_refs
    uint32_t max_dev_to_best
        [TOT_INTER_GROUP]; // 0: OFF; 1: limit the injection to the best references based on distortion
    uint8_t use_tpl_info_offset;
    uint8_t check_closest_multiplier;
    uint8_t closest_refs
        [TOT_INTER_GROUP]; // 0: OFF; 1: limit the injection to the closest references based on distance (LAST/BWD)
} RefPruningControls;
typedef struct DepthRemovalCtrls {
    uint8_t enabled;
    uint8_t
        disallow_below_64x64; // remove 32x32 blocks and below based on the sb_64x64 (me_distortion, variance)
    uint8_t
        disallow_below_32x32; // remove 16x16 blocks and below based on the sb_64x64 (me_distortion, variance)
    uint8_t
        disallow_below_16x16; // remove 8x8 blocks and below based on the sb_64x64 (me_distortion, variance)
    uint8_t
        disallow_4x4; // remove 8x8 blocks and below based on the sb_64x64 (me_distortion, variance)
} DepthRemovalCtrls;
typedef struct DepthCtrls {
    int8_t
        s_depth; // start depth; 0: consider no parent blocks; else number of parent blocks to consider, specified as a negative number (e.g. -2 means consider 2 parents)
    int8_t
         e_depth; // end depth; 0: consider no child blocks; else number of child blocks to consider, specified as a positive number (e.g. 2 means consider 2 children)
    Bool allow_nsq_in_child_depths; // Allow NSQ in depths below the PD0-selected depths.
} DepthCtrls;
#define MAX_RANGE_CNT 8
#define MAX_RANGE_CNT 8
typedef struct DepthRefinementCtrls {
    uint8_t enabled;
    int64_t
        parent_to_current_th; // maximum allowed parent-to-current cost deviation beyond which the previous depth will not be added to PRED
    int64_t
        sub_to_current_th; // maximum allowed sub-to-current cost deviation beyond which the next depth will not be added to PRED
    uint8_t
        up_to_2_depth; // when 1, a maximum of 2 depth per block (PRED+Parent or PRED+Sub), 0: no restriction(s)
    uint8_t
             cost_band_based_modulation; // whether to decrement parent_to_current_th and sub_to_current_th based on the cost range of the parent block or not
    uint16_t max_cost_multiplier; // the max cost beyond which the decrement is ignored
    uint8_t  max_band_cnt; // the number of band(s)
    int64_t  decrement_per_band[MAX_RANGE_CNT]; // the offset per band
    int64_t
        limit_4x4_depth; // an offset (towards more aggressive pruning) when the child-depth is 4x4
    unsigned int
        sub_to_current_pd0_coeff_th; // If current depth has fewer than sub_to_current_pd0_coeff_th PD0 coeffs, subtract sub_to_current_pd0_coeff_offset from
    // sub_to_current_th towards skipping child depths (doesn't work with very light PD0)
    int sub_to_current_pd0_coeff_offset;
} DepthRefinementCtrls;
typedef struct SubresCtrls {
    uint8_t step; // Residual sub-sampling step (0:OFF)
    uint8_t
        odd_to_even_deviation_th; // Set step to 0 if the deviation between: (1) the pred-to-src SAD of even rows and (2) the pred-to-src SAD of odd rows
    // of the 1st 64x64 block @ mds3 of PD0 is higher than odd_to_even_deviation_th
    // when 0, the detection is OFF
} SubresCtrls;
typedef struct PfCtrls {
    EB_TRANS_COEFF_SHAPE pf_shape;
} PfCtrls;
typedef struct MdNsqMotionSearchCtrls {
    uint8_t enabled; // 0: NSQ motion search @ MD OFF; 1: NSQ motion search @ MD ON
    uint8_t use_ssd; // 0: search using SAD; 1: search using SSD
    uint8_t full_pel_search_width; // Full Pel search area width
    uint8_t full_pel_search_height; // Full Pel search area height
    uint8_t enable_psad; //Enable pSad
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
    uint8_t  enable_psad; // Enable pSad
} MdSqMotionSearchCtrls;
typedef struct MdPmeCtrls {
    uint8_t enabled; // 0: PME search @ MD OFF; 1: PME search @ MD ON
    uint8_t use_ssd; // 0: search using SAD; 1: search using SSD
    int     early_check_mv_th_multiplier; // Do not perform PME search for blocks that have a valid ME_MV unless the ME_MV has a different direction than all MVP(s) and the ME_MV mag is higher than MV_TH=f(early_check_mv_th_multiplier)
    uint8_t full_pel_search_width; // Full Pel search area width
    uint8_t full_pel_search_height; // Full Pel search area height
    int     pre_fp_pme_to_me_cost_th; // If pre_fp_pme_to_me_cost higher than pre_fp_pme_to_me_cost_th then PME_MV = ME_MV and exit (decrease towards a faster level)
    int     pre_fp_pme_to_me_mv_th; // If pre_fp_pme_to_me_mv smaller than pre_fp_pme_to_me_mv_th then PME_MV = ME_MV and exit (increase towards a faster level)
    int     post_fp_pme_to_me_cost_th; // If post_fp_pme_to_me_cost higher than post_fp_pme_to_me_cost_th then PME_MV = ME_MV and exit (decrease towards a faster level)
    int     post_fp_pme_to_me_mv_th; // If post_fp_pme_to_me_mv smaller than post_fp_pme_to_me_mv_th then PME_MV = ME_MV and exit (increase towards a faster level)
    uint8_t
        modulate_pme_for_blk_size_res; // only turn pme on for 32x32 blocks or 1080p 64x64 blocks
    uint8_t enable_psad; // Enable pSad
} MdPmeCtrls;
typedef struct MdSubPelSearchCtrls {
    uint8_t
        enabled; // Specifies whether the Sub-Pel search will be performed or not (0: OFF, 1: ON)
    SUBPEL_SEARCH_TYPE
    subpel_search_type; // Specifies the interpolation filter tap (1: 2-tap filter, 2: 4-tap filter, 3: 8-tap filter)
    SUBPEL_FORCE_STOP
    max_precision; // Specifies the refinement precision (or number of rounds) (0: 1/8-Pel (3 rounds), 1: 1/4-Pel (2 rounds), 2: 1/2-Pel (1 round), 3: Full-Pel-no refinement (0 round))
    SUBPEL_SEARCH_METHODS
    subpel_search_method; // Specifies whether pruning will be applied to 1/2-Pel position(s) or not (SUBPEL_TREE: No, SUBPEL_TREE_PRUNED: YES)
    int subpel_iters_per_step; // Specifies the maximum number of steps in logarithmic subpel search before giving up
    int pred_variance_th; // Specifies the Full-Pel prediction-block-variance threshold under which the Sub-Pel search is not performed; do not perform Sub-Pel if the variance of the Full-Pel prediction-block is low (where interpolation will unlikely modify the Full-Pel samples)
    uint8_t
            abs_th_mult; // Specifies the Full-Pel prediction-block-error-threshold below which the Sub-Pel search is not performed; do not perform Sub-Pel if the prediction-block-error is already low
    int     round_dev_th; // Specifies the prediction-block-error deviation threshold between round-(N-1) and round-(N-2) under which the refinement is paused; pause the refinement if the prediction-block-error is not getting better through the process (the check takes place at only the 2nd round (prior to the 1/4-Pel refinement) and the 3rd round (prior to the 1/8-Pel refinement).
    uint8_t skip_diag_refinement; // Specifies the refinement accuracy for diagonal position(s).
    uint8_t
        skip_zz_mv; // Specifies whether the Sub-Pel search will be performed for around (0,0) or not (0: OFF, 1: ON)
} MdSubPelSearchCtrls;
typedef struct ParentSqCoeffAreaBasedCyclesReductionCtrls {
    Bool enabled;

    uint8_t high_freq_band1_th; // cutoff for the highest coeff-area band [0-100]
    uint8_t
        high_freq_band1_level; // level of action to use if luma coeff-area of parent SQ is >= high_freq_band1_th
    uint8_t
        high_freq_band2_th; // cutoff for the second high coeff-area band [0-100]; should be less than high_freq_band1_th
    uint8_t
        high_freq_band2_level; // level of action to use if luma coeff-area of parent SQ is >= high_freq_band2_th
    uint8_t
        high_freq_band3_th; // cutoff for the third high coeff-area band [0-100]; should be less than high_freq_band2_th
    uint8_t
        high_freq_band3_level; // level of action to use if luma coeff-area of parent SQ is >= high_freq_band3_th

    uint8_t
            enable_zero_coeff_action; // enable for whether to apply action when parent SQ has 0 luma coefficients
    uint8_t zero_coeff_action; // level of action to use if parent SQ has 0 luma coeffs
    uint8_t
            enable_one_coeff_action; // enable for whether to apply action when parent SQ has 1 luma coefficients
    uint8_t one_coeff_action; // level of action to use if parent SQ has 1 luma coeff

    uint8_t
        low_freq_band1_th; // cutoff for the lowest coeff-area band [0-100]; should be less than high_freq_band2_th
    uint8_t
        low_freq_band1_level; // level of action to use if luma coeff-area of parent SQ is < low_freq_band1_th
    uint8_t
        low_freq_band2_th; // cutoff for the lowest coeff-area band [0-100]; should be less than high_freq_band2_th and larger than low_freq_band1_th
    uint8_t
        low_freq_band2_level; // level of action to use if luma coeff-area of parent SQ is < low_freq_band2_th
} ParentSqCoeffAreaBasedCyclesReductionCtrls;

typedef struct RdoqCtrls {
    uint8_t enabled;

    uint8_t
        eob_fast_y_inter; // 0: do not use eob_fast  for luma inter; 1:  use eob_fast  for luma inter
    uint8_t
        eob_fast_y_intra; // 0: do not use eob_fast  for luma intra; 1:  use eob_fast  for luma intra
    uint8_t
        eob_fast_uv_inter; // 0: do not use eob_fast  for chroma inter; 1:  use eob_fast  for chroma inter
    uint8_t
        eob_fast_uv_intra; // 0: do not use eob_fast  for chroma intra; 1:  use eob_fast  for chroma intra
    uint8_t fp_q_y; // 0: use default quant for luma; 1: use fp_quant for luma
    uint8_t fp_q_uv; // 0: use default quant for chroma; 1: use fp_quant for chroma
    uint8_t satd_factor; // Do not perform rdoq if the tx satd > satd_factor
    uint8_t early_exit_th; // Do not perform rdoq based on an early skip/non-skip cost
    uint8_t skip_uv; // [MD Only] 0: RDOQ for both Luma & Chroma, 1: RDOQ for only Luma
    uint8_t dct_dct_only; // [MD Only] 0: RDOQ for All txt type(s), 1: RDOQ for only DCT_DCT
    uint8_t eob_th; // eob_th beyond which RDOQ is shut
    uint8_t eob_fast_th;
} RdoqCtrls;
typedef struct NicScalingCtrls {
    uint8_t stage1_scaling_num; // Scaling numerator for post-stage 0 NICS: <x>/16
    uint8_t stage2_scaling_num; // Scaling numerator for post-stage 1 NICS: <x>/16
    uint8_t stage3_scaling_num; // Scaling numerator for post-stage 2 NICS: <x>/16
} NicScalingCtrls;
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
    uint64_t mds1_cand_base_th_intra; // Post-MDS0 candidate pruning TH for intra classes
    uint64_t mds1_cand_base_th_inter; // Post-MDS0 candidate pruning TH for inter classes

    uint16_t
        mds1_cand_th_rank_factor; // Use a more aggressive MDS0 candidate pruning TH based on how the candidate is ranked compared to others
    // (2nd best candidate will use a safer pruning TH than the 3rd best, and so on).
    // The new_th = mds1_cand_base_th / (mds1_cand_th_rank_factor * cand_count); where cand_count is the number of
    // candidates with better cost than the current candidate. 0 is off. Higher is more aggressive.

    // Post mds1
    uint64_t mds2_cand_base_th;
    uint16_t
        mds2_cand_th_rank_factor; // Use a more aggressive MDS1 candidate pruning TH based on how the candidate is ranked compared to others
    // (2nd best candidate will use a safer pruning TH than the 3rd best, and so on).
    // The new_th = mds2_cand_base_th / (mds2_cand_th_rank_factor * cand_count); where cand_count is the number of
    // candidates with better cost than the current candidate. 0 is off. Higher is more aggressive.
    uint16_t
        mds2_relative_dev_th; // Prune candidates based on the relative deviation compared to the previous candidate. Deviations are still computed
    // compared to the best candidate. If dev > previous_dev + th then prune the candidate.

    // Post mds2
    uint64_t mds3_cand_base_th;

    Bool enable_skipping_mds1; // enable skipping MDS1 in PD1 when there is only 1 cand post-mds0
    uint32_t
        force_1_cand_th; // if (best_mds0_distortion/QP < TH) consider only the best candidate after MDS0; 0: OFF, higher: more aggressive.
} NicPruningCtrls;
typedef struct NicCtrls {
    NicPruningCtrls pruning_ctrls;
    NicScalingCtrls scaling_ctrls;
    uint8_t         md_staging_mode; // to specify the number of md-stage(s)
} NicCtrls;
typedef struct CandEliminationCtlrs {
    uint32_t enabled;
    uint8_t  dc_only;
    uint8_t  inject_new_me;
    uint8_t  inject_new_pme;
    uint8_t  inject_new_warp;
    uint8_t  th_multiplier; // factor to scale base TH by for distortion check
} CandEliminationCtlrs;
typedef struct TxsControls {
    uint8_t  enabled;
    uint8_t  prev_depth_coeff_exit; // Skip current depth if previous depth has no coeff
    uint8_t  intra_class_max_depth_sq; // Max number of depth(s) for INTRA classes in SQ blocks
    uint8_t  intra_class_max_depth_nsq; // Max number of depth(s) for INTRA classes in NSQ blocks
    uint8_t  inter_class_max_depth_sq; // Max number of depth(s) for INTER classes in SQ blocks
    uint8_t  inter_class_max_depth_nsq; // Max number of depth(s) for INTER classes in NSQ blocks
    int      depth1_txt_group_offset; // Offset to be subtracted from default txt-group to derive the txt-group of depth-1
    int      depth2_txt_group_offset; // Offset to be subtracted from default txt-group to derive the txt-group of depth-2
    uint16_t min_sq_size; // Min. sq size to use TXS for
} TxsControls;
typedef struct WmCtrls {
    uint8_t enabled;
    uint8_t use_wm_for_mvp; // allow/disallow MW for MVP candidates
    uint8_t
        refinement_iterations; // Number of iterations to use in the refinement search; each iteration searches the cardinal neighbours around the best-so-far position; 0 is no refinement
} WmCtrls;
typedef struct UvCtrls {
    uint8_t enabled;
    uint8_t uv_mode; // Indicates the chroma search mode
        // CHROMA_MODE_0: Full chroma search @ MD
        // CHROMA_MODE_1: Fast chroma search @ MD - No independent chroma mode search
        // CHROMA_MODE_2: Chroma blind @ MD
    uint8_t
             nd_uv_serach_mode; // Non-direct chroma search 0: pre chroma search is used, 1: chroma search at last md_stage is used
    uint8_t  uv_nic_scaling_num; // Scaling numerator for independent chroma NICS: <x>/16
    uint32_t uv_intra_th; // Threshold to skip non-directional chroma search.
    uint32_t uv_cfl_th; // Threshold to skip cfl.
} UvCtrls;
typedef struct InterpolationSearchCtrls {
    IfsLevel
        level; // Specifies the MD Stage where the interpolation filter search will take place (IFS_MDS0, IFS_MDS1, IFS_MDS2, or IFS_MDS3 for respectively MD Stage 0, MD Stage 1, MD Stage 2, and MD Stage 3)
    uint8_t
        quarter_pel_only; // Specifies whether the interpolation filter search will use 1/8-Pel precision or 1/4-Pel precision (0: 1/8-Pel precision, 1: 1/4-Pel precision)
    uint8_t
        early_skip; // Specifies whether an early interpolation filter search exit could take place based on the cost of signaling a switchable filter type (0: OFF, 1: ON)
    uint8_t
        subsampled_distortion; // Specifies whether sub-sampled input/prediction will be used at the distortion computation (0: OFF, 1: ON, NA for block height 16 and lower)
    uint8_t
        skip_sse_rd_model; // Specifies whether a model wll be used for rate estimation or not (0: NO (assume rate is 0), 1: estimate rate from distortion)
} InterpolationSearchCtrls;
typedef struct SpatialSSECtrls {
    Bool spatial_sse_full_loop_level; // enable spatial-sse for each superblock
} SpatialSSECtrls;
typedef struct RedundantCandCtrls {
    int score_th;
    int mag_th;
} RedundantCandCtrls;
typedef struct UseNeighbouringModeCtrls {
    uint8_t enabled;
} UseNeighbouringModeCtrls;
typedef struct BlockLocation {
    uint32_t input_origin_index; //luma   block location in picture
    uint32_t input_cb_origin_in_index; //chroma block location in picture
    uint32_t blk_origin_index; //luma   block location in SB
    uint32_t blk_chroma_origin_index; //chroma block location in SB
} BlockLocation;

typedef struct Lpd0Ctrls {
    Pd0Level
         pd0_level; // Whether light-PD0 is set to be used for an SB (the detector may change this)
    Bool use_lpd0_detector
        [LPD0_LEVELS]; // Whether to use a detector; if use_light_pd0 is set to 1, the detector will protect tough SBs
    Bool use_ref_info
        [LPD0_LEVELS]; // Use info of ref frames - incl. colocated SBs - such as mode, coeffs, etc. in the detector
    uint32_t me_8x8_cost_variance_th
        [LPD0_LEVELS]; // me_8x8_cost_variance_th beyond which the PD0 is used (instead of light-PD0)
} Lpd0Ctrls;

typedef struct Lpd1Ctrls {
    int8_t
         pd1_level; // Whether light-PD1 is set to be used for an SB (the detector may change this)
    Bool use_lpd1_detector
        [LPD1_LEVELS]; // Whether to use a detector; if use_light_pd1 is set to 1, the detector will protect tough SBs
    Bool use_ref_info
        [LPD1_LEVELS]; // Use info of ref frames - incl. colocated SBs - such as mode, coeffs, etc. in the detector
    uint32_t cost_th_dist[LPD1_LEVELS]; // Distortion value used in cost TH for detector
    uint32_t coeff_th[LPD1_LEVELS]; // Num non-zero coeffs used in detector
    uint16_t max_mv_length
        [LPD1_LEVELS]; // Max MV length TH used in the detector: 0 - (0,0) MV only; (uint16_t)~0 means no MV check (higher is more aggressive)
    uint32_t me_8x8_cost_variance_th
        [LPD1_LEVELS]; // me_8x8_cost_variance_th beyond which the PD1 is used (instead of light-PD1)
    uint32_t skip_pd0_edge_dist_th
        [LPD1_LEVELS]; // ME_64x64_dist threshold used for edge SBs when PD0 is skipped
    uint16_t skip_pd0_me_shift
        [LPD1_LEVELS]; // Shift applied to ME dist and var of top and left SBs when PD0 is skipped
} Lpd1Ctrls;

typedef struct DetectHighFreqCtrls {
    int8_t enabled;
    uint16_t
        me_8x8_sad_var_th; // me-8x8 SADs deviation threshold beyond which the SB is not considered
    uint16_t high_satd_th; // 32x32 satd threshold beyond which the SB is tagged
    uint16_t
            satd_to_sad_dev_th; // me-SAD-to-SATD deviation of the 32x32 blocks threshold beyond which the SB is tagged (~2x is the fundamental deviation)
    uint8_t depth_removal_shift; // depth-removal level left-shift for the detected SB(s)
    uint8_t max_pic_lpd0_lvl; // maximum lpd0 level for the detected SB(s)
    uint8_t max_pic_lpd1_lvl; // maximum lpd1 level for the detected SB(s)
    uint8_t max_pd1_txt_lvl; // maximum pd1-txt level for the detected SB(s)
} DetectHighFreqCtrls;

typedef struct Lpd1TxCtrls {
    uint8_t
        zero_y_coeff_exit; // skip cost calc and chroma TX/compensation if there are zero luma coeffs
    uint8_t
        chroma_detector_level; // Control aggressiveness of chroma detector (used to skip chroma TX when luma has 0 coeffs): 0: OFF, 1: saftest, 2, 3: medium
    uint8_t
             skip_nrst_nrst_luma_tx; // Skip luma TX for NRST_NRST candidates if dist/QP is low and if top and left neighbours have no coeffs and are NRST_NRST
    uint16_t skip_tx_th; // if (skip_tx_th/QP < TH) skip TX at MDS3; 0: OFF, higher: more aggressive
    uint32_t
        use_mds3_shortcuts_th; // if (best_mds0_distortion/QP < TH) use shortcuts for candidate at MDS3; 0: OFF, higher: more aggressive
    uint8_t
        use_neighbour_info; // if true, use info from neighbouring blocks to use more aggressive THs/actions
    uint8_t
        use_uv_shortcuts_on_y_coeffs; // Apply shortcuts to the chroma TX path if luma has few coeffs
} Lpd1TxCtrls;
typedef struct CflCtrls {
    Bool    enabled;
    uint8_t itr_th; // Early exit to reduce the number of iterations to compute CFL parameters
} CflCtrls;
typedef struct MdRateEstCtrls {
    Bool
         update_skip_ctx_dc_sign_ctx; // If true, update skip context and dc_sign context (updates are done in the same func, so control together)
    Bool update_skip_coeff_ctx; // If true, update skip coeff context
    uint8_t
           coeff_rate_est_lvl; // 0: OFF (always use approx for coeff rate), 1: full; always compute coeff rate, 2: when num_coeff is low, use approximation for coeff rate
    int8_t lpd0_qp_offset; // Offset applied to qindex at quantization in LPD0
    uint8_t
        pd0_fast_coeff_est_level; // estimate the rate of the first (eob/N) coeff(s) and last coeff only
} MdRateEstCtrls;
typedef struct IntraCtrls {
    uint8_t enable_intra;
    uint8_t
        intra_mode_end; // the last intra prediciton mode generated starting from DC_PRED, min: DC_PRED, max: PAETH_PRED
    uint8_t
        angular_pred_level; // 0: angular off; 1: angular full; 2/3: limit num. angular candidates; 4: H + V only
} IntraCtrls;
typedef struct TxShortcutCtrls {
    uint8_t bypass_tx_when_zcoeff; // Skip TX at MDS3 if the MDS1 TX gave 0 coeffs
    uint8_t apply_pf_on_coeffs; // Apply pf based on the number of coeffs
    uint8_t
        chroma_detector_level; // Use a detector to protect chroma from aggressive actions based on luma info: 0: OFF, 1: saftest, 2, 3: medium
    uint32_t
        use_mds3_shortcuts_th; // if (best_mds0_distortion/QP < TH) use shortcuts for candidates at MDS3; 0: OFF, higher: more aggressive
    uint8_t
        use_neighbour_info; // if true, use info from neighbouring blocks to use more aggressive THs/actions
} TxShortcutCtrls;
typedef struct Mds0Ctrls {
    uint8_t mds0_dist_type; // Distortion metric to use MDS0: SSD, VAR, SAD
    uint8_t
        enable_cost_based_early_exit; // Skip cost computation if distortion is mds0_distortion_th % higher than best candidate cost (applies to reg. PD1 only)
    uint16_t
        mds0_distortion_th; // % TH used to compare candidate distortion to best cost; higher is safer (applies to reg. PD1 only)
} Mds0Ctrls;
typedef struct CandReductionCtrls {
    uint8_t            merge_inter_classes;
    RedundantCandCtrls redundant_cand_ctrls;
    NearCountCtrls     near_count_ctrls;
    uint8_t lpd1_mvp_best_me_list; // inject unipred MVP candidates only for the best ME list
    UseNeighbouringModeCtrls use_neighbouring_mode_ctrls;
    CandEliminationCtlrs     cand_elimination_ctrls;
    uint8_t                  reduce_unipred_candidates;
    uint8_t                  mds0_reduce_intra;
} CandReductionCtrls;
typedef struct Skip4x4DepthCtrls {
    uint8_t enabled;
    float
        min_distortion_cost_ratio; // the distortion-to-cost ratio under wich the quad_deviation_th is zeroed out (feature is disabled)
    float
        quad_deviation_th; // do not perform sub_depth if std_deviation of the 4 quadrants src-to-rec dist is less than std_deviation_th
    uint8_t
        skip_all; // whether to skip all or only next depth; 0: skip only next depth; 1: skip all lower depths
} Skip4x4DepthCtrls;
typedef struct ModeDecisionContext {
    EbDctor dctor;

    EbFifo                       *mode_decision_configuration_input_fifo_ptr;
    EbFifo                       *mode_decision_output_fifo_ptr;
    ModeDecisionCandidate       **fast_candidate_ptr_array;
    ModeDecisionCandidate        *fast_candidate_array;
    ModeDecisionCandidateBuffer **candidate_buffer_ptr_array;
    ModeDecisionCandidateBuffer  *candidate_buffer_tx_depth_1;
    ModeDecisionCandidateBuffer  *candidate_buffer_tx_depth_2;
    MdRateEstimationContext      *md_rate_estimation_ptr;
    MdRateEstimationContext      *rate_est_table;
    MdBlkStruct                  *md_local_blk_unit;
    BlkStruct                    *md_blk_arr_nsq;
    uint8_t                      *avail_blk_flag;
    uint8_t                      *tested_blk_flag; //tells whether this CU is tested in MD.
    uint8_t                      *do_not_process_blk;
    MdcSbData                    *mdc_sb_array;

    NeighborArrayUnit *intra_luma_mode_neighbor_array;
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
    NeighborArrayUnit                   *
        cb_dc_sign_level_coeff_neighbor_array; // Stored per 4x4. 8 bit: lower 6 bits(COEFF_CONTEXT_BITS), shows if there is at least one Coef. Top 2 bit store the sign of DC as follow: 0->0,1->-1,2-> 1
    NeighborArrayUnit    *txfm_context_array;
    NeighborArrayUnit    *leaf_partition_neighbor_array;
    NeighborArrayUnit    *skip_coeff_neighbor_array;
    struct EncDecContext *enc_dec_context_ptr;

    uint64_t *fast_cost_array;
    uint64_t *full_cost_array;
    // Lambda
    uint32_t fast_lambda_md[2];
    uint32_t full_lambda_md[2];
    uint32_t full_sb_lambda_md
        [2]; // for the case of lambda modulation (blk_lambda_tuning), full_lambda_md/fast_lambda_md corresponds
    // to block lambda and full_sb_lambda_md is the full lambda per sb
    Bool blk_lambda_tuning;
    //  Context Variables---------------------------------
    SuperBlock      *sb_ptr;
    BlkStruct       *blk_ptr;
    const BlockGeom *blk_geom;
    PredictionUnit  *pu_ptr;
    PALETTE_BUFFER   palette_buffer;
    PaletteInfo      palette_cand_array[MAX_PAL_CAND];
    // MD palette search
    uint8_t         *palette_size_array_0;
    uint8_t         *palette_size_array_1;
    uint8_t          sb64_sq_no4xn_geom; //simple geometry 64x64SB, Sq only, no 4xN
    uint8_t          pu_itr;
    uint32_t        *best_candidate_index_array;
    uint16_t         blk_origin_x;
    uint16_t         blk_origin_y;
    uint32_t         sb_origin_x;
    uint32_t         sb_origin_y;
    uint32_t         round_origin_x;
    uint32_t         round_origin_y;
    uint8_t          hbd_mode_decision;
    uint8_t          encoder_bit_depth;
    uint8_t          qp_index;
    uint8_t          me_q_index;
    uint64_t         three_quad_energy;
    uint32_t         txb_1d_offset;
    Bool             uv_intra_comp_only;
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
    Mv               **
        injected_mvs; // Track all MVs that are prepared for candidates prior to MDS0. Used to avoid MV duplication.
    MvReferenceFrame *injected_ref_types; // Track the reference types for each MV
    uint16_t          injected_mv_count;
    uint32_t          me_block_offset;
    uint32_t          me_cand_offset;
    // Pointer to a scratch buffer used by CFL & IFS
    EbPictureBufferDesc *scratch_prediction_ptr;
    uint8_t              tx_depth;
    uint8_t              txb_itr;
    uint32_t             me_sb_addr;
    uint32_t             geom_offset_x;
    uint32_t             geom_offset_y;
    int16_t              luma_txb_skip_context;
    int16_t              luma_dc_sign_context;
    int16_t              cb_txb_skip_context;
    int16_t              cb_dc_sign_context;
    int16_t              cr_txb_skip_context;
    int16_t              cr_dc_sign_context;
    // Multi-modes signal(s)
    uint8_t              global_mv_injection;
    uint8_t              new_nearest_injection;
    uint8_t              new_nearest_near_comb_injection;
    WmCtrls              wm_ctrls;
    UvCtrls              uv_ctrls;
    uint8_t              unipred3x3_injection;
    uint8_t              bipred3x3_injection;
    uint8_t              redundant_blk;
    uint8_t              nic_level;
    uint8_t              similar_blk_avail;
    uint16_t             similar_blk_mds;
    uint8_t              inject_inter_candidates;
    uint8_t             *cfl_temp_luma_recon;
    uint16_t            *cfl_temp_luma_recon16bit;
    Bool                 spatial_sse_full_loop_level;
    Bool                 blk_skip_decision;
    int8_t               rdoq_level;
    int16_t              sb_me_mv[BLOCK_MAX_COUNT_SB_128][MAX_NUM_OF_REF_PIC_LIST][MAX_REF_IDX][2];
    MV                   fp_me_mv[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];
    MV                   sub_me_mv[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];
    uint32_t             post_subpel_me_mv_cost[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];
    int16_t              best_pme_mv[MAX_NUM_OF_REF_PIC_LIST][MAX_REF_IDX][2];
    int8_t               valid_pme_mv[MAX_NUM_OF_REF_PIC_LIST][MAX_REF_IDX];
    EbPictureBufferDesc *input_sample16bit_buffer;
    uint8_t  hbd_pack_done; //set to 1 once the packing of 10bit source is done for each SB
    uint16_t tile_index;
    DECLARE_ALIGNED(16, uint8_t, pred0[2 * MAX_SB_SQUARE]);
    DECLARE_ALIGNED(16, uint8_t, pred1[2 * MAX_SB_SQUARE]);
    DECLARE_ALIGNED(32, int16_t, residual1[MAX_SB_SQUARE]);
    DECLARE_ALIGNED(32, int16_t, diff10[MAX_SB_SQUARE]);
    MdStage   md_stage;
    uint32_t *cand_buff_indices[CAND_CLASS_TOTAL];
    uint8_t   bypass_md_stage_1;
    uint8_t   bypass_md_stage_2;
    uint32_t  md_stage_0_count[CAND_CLASS_TOTAL];
    uint32_t  md_stage_1_count[CAND_CLASS_TOTAL];
    uint32_t  md_stage_2_count[CAND_CLASS_TOTAL];
    uint32_t  md_stage_3_count[CAND_CLASS_TOTAL];
    uint32_t  md_stage_1_total_count;
    uint32_t  md_stage_2_total_count;
    uint32_t  md_stage_3_total_count;
    uint32_t  md_stage_3_total_intra_count;
    uint64_t  best_intra_cost;
    uint64_t  best_inter_cost;
    CandClass target_class;
    uint8_t   perform_mds1;
    uint8_t   use_tx_shortcuts_mds3;
    uint8_t   lpd1_allow_skipping_tx;
    // fast_loop_core signals
    Bool md_staging_skip_interpolation_search;
    Bool md_staging_skip_chroma_pred;
    // full_loop_core signals
    Bool
         md_staging_perform_inter_pred; // 0: perform luma & chroma prediction + interpolation search, 2: nothing (use information from previous stages)
    Bool md_staging_tx_size_mode; // 0: Tx Size recon only, 1:Tx Size search and recon
    Bool md_staging_txt_level;
    Bool md_staging_skip_full_chroma;
    Bool md_staging_skip_rdoq;
    Bool md_staging_spatial_sse_full_loop_level;
    Bool md_staging_perform_intra_chroma_pred;
    DECLARE_ALIGNED(
        16, uint8_t,
        intrapred_buf[INTERINTRA_MODES][2 * 32 * 32]); //MAX block size for inter intra is 32x32
    DECLARE_ALIGNED(16, uint8_t, obmc_buff_0[2 * 2 * MAX_MB_PLANE * MAX_SB_SQUARE]);
    DECLARE_ALIGNED(16, uint8_t, obmc_buff_1[2 * 2 * MAX_MB_PLANE * MAX_SB_SQUARE]);
    DECLARE_ALIGNED(16, int32_t, wsrc_buf[MAX_SB_SQUARE]);
    DECLARE_ALIGNED(16, int32_t, mask_buf[MAX_SB_SQUARE]);
    unsigned int pred_sse[REF_FRAMES];
    uint8_t     *above_txfm_context;
    uint8_t     *left_txfm_context;
    // square cost weighting for deciding if a/b shapes could be skipped
    uint32_t sq_weight;
    uint32_t max_part0_to_part1_dev;
    uint32_t
        skip_hv4_on_best_part; // if true, skip H4/V4 shapes when best partition so far is not H/V
    uint8_t
        md_depth_early_exit_th; // Skip testing remaining blocks at the current depth if (curr_cost * 100 > pic_depth_early_exit_th * parent_cost)
    // [0-100], 0 is OFF, lower percentage is more aggressive
    IntraCtrls        intra_ctrls;
    MdRateEstCtrls    rate_est_ctrls;
    uint8_t           shut_fast_rate; // use coeff rate and slipt flag rate only (no MVP derivation)
    uint8_t           md_staging_fast_coeff_est_level; // Control fast_coeff_est_level per mds
    uint8_t           md_staging_subres_step; // Control subres_step per mds
    uint8_t           md_pic_obmc_level;
    uint8_t           md_inter_intra_level;
    uint8_t           md_filter_intra_level;
    uint8_t           md_allow_intrabc;
    uint8_t           md_palette_level;
    uint8_t           dist_based_ref_pruning;
    DepthRemovalCtrls depth_removal_ctrls;
    DepthCtrls        depth_ctrls; // control which depths can be considered in PD1
    DepthRefinementCtrls depth_refinement_ctrls;
    Skip4x4DepthCtrls    skip_4x4_depth_ctrls;
    int64_t              parent_to_current_deviation;
    int64_t              child_to_current_deviation;
    SubresCtrls          subres_ctrls;
    uint8_t              is_subres_safe;
    PfCtrls              pf_ctrls;
    // Control signals for MD sparse search (used for increasing ME search for active clips)
    MdSqMotionSearchCtrls  md_sq_me_ctrls;
    MdNsqMotionSearchCtrls md_nsq_motion_search_ctrls;
    MdPmeCtrls             md_pme_ctrls;
    uint8_t                md_subpel_me_level;
    MdSubPelSearchCtrls    md_subpel_me_ctrls;
    uint8_t                md_subpel_pme_level;
    MdSubPelSearchCtrls    md_subpel_pme_ctrls;
    PmeResults             pme_res[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];
    ObmcControls           obmc_ctrls;
    InterCompCtrls         inter_comp_ctrls;
    InterIntraCompCtrls    inter_intra_comp_ctrls;
    RefResults ref_filtering_res[TOT_INTER_GROUP][MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];
    RefPruningControls ref_pruning_ctrls;
    // Signal to control initial and final pass PD setting(s)
    PdPass                                     pd_pass;
    CflCtrls                                   cfl_ctrls;
    TxsControls                                txs_ctrls;
    TxtControls                                txt_ctrls;
    CandReductionCtrls                         cand_reduction_ctrls;
    RdoqCtrls                                  rdoq_ctrls;
    uint8_t                                    disallow_4x4;
    uint8_t                                    md_disallow_nsq;
    ParentSqCoeffAreaBasedCyclesReductionCtrls parent_sq_coeff_area_based_cycles_reduction_ctrls;
    uint8_t                                    sb_size;
    EbPictureBufferDesc                       *recon_coeff_ptr[TX_TYPES];
    EbPictureBufferDesc                       *recon_ptr[TX_TYPES];
    EbPictureBufferDesc                       *quant_coeff_ptr[TX_TYPES];
    EbPictureBufferDesc *tx_coeffs; // buffer used to store transformed coeffs during TX/Q/IQ.
        // TX'd coeffs are only needed temporarily, so no need to save for each TX type.

    uint8_t              skip_intra;
    EbPictureBufferDesc *temp_residual_ptr;
    EbPictureBufferDesc *temp_recon_ptr;
    // Array for all nearest/near MVs for a block for single ref case
    MV mvp_array[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH][MAX_MVP_CANIDATES];
    // Count of all nearest/near MVs for a block for single ref case
    int8_t mvp_count[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];
    // Start/end position for MD sparse search
    int16_t         sprs_lev0_start_x;
    int16_t         sprs_lev0_end_x;
    int16_t         sprs_lev0_start_y;
    int16_t         sprs_lev0_end_y;
    NicCtrls        nic_ctrls;
    uint8_t         inter_compound_mode;
    MV              ref_mv;
    uint16_t        sb_index;
    uint64_t        mds0_best_cost;
    uint8_t         mds0_best_class;
    uint32_t        mds0_best_idx;
    CandClass       mds0_best_class_it;
    uint32_t        mds1_best_idx;
    CandClass       mds1_best_class_it;
    Mds0Ctrls       mds0_ctrls;
    uint32_t        md_me_dist;
    uint8_t         inject_new_me;
    uint8_t         inject_new_pme;
    uint8_t         inject_new_warp;
    TxShortcutCtrls tx_shortcut_ctrls;
    uint64_t        estimate_ref_frames_num_bits[MODE_CTX_REF_FRAMES]; // [TOTAL_REFS_PER_FRAME + 1]
    uint32_t        max_nics; // Maximum number of candidates MD can support
    uint32_t        max_nics_uv; // Maximum number of candidates MD can support
    InterpolationSearchCtrls ifs_ctrls;
    Bool bypass_encdec; // If enabled, will bypass EncDec and copy recon/quant coeffs from MD
    Bool
        pred_depth_only; // Indicates whether only pred depth refinement is used in PD1 (set per SB)
    Bool
        pic_pred_depth_only; // Indicates whether only pred depth refinement is used in PD1 (set per frame)
    // Per frame is necessary because some shortcuts can only be taken if the whole frame
    // uses pred depth only
    uint16_t            coded_area_sb;
    uint16_t            coded_area_sb_uv;
    Lpd0Ctrls           lpd0_ctrls;
    DetectHighFreqCtrls detect_high_freq_ctrls;
    // 0 : Use regular PD0
    // 1 : Use light PD0 path.  Assumes one class, no NSQ, no 4x4, TXT off, TXS off, PME off, etc.
    // 2 : Use very light PD0 path: only mds0 (no transform path), no compensation(s) @ mds0
    //     (only umpired candidates, and read directly from reference buffer(s) to compute distortion(s) without
    //     accessing the pred buffer(s)), SSD as distortion metric(kept the use of full lambda @ in-depth-exit and @ inter-depth),
    //     only split flag for rate(coefficient(s) rate is assumed to be 0),
    //     use the 2nd PD1-level classifier (as the regular PD1 classifier uses the number of non-zero coefficient(s)).
    // 3: Skip pd0 if block size is equal to or greater than 32x32
    Lpd1Ctrls       lpd1_ctrls;
    SpatialSSECtrls spatial_sse_ctrls;

    uint16_t init_max_block_cnt;
    uint8_t  end_plane;
    uint8_t
        need_hbd_comp_mds3; // set to true if MDS3 needs to perform a full 10bit compensation in MDS3 (to make MDS3 conformant when using bypass_encdec)
    uint8_t
                approx_inter_rate; // use approximate rate for inter cost (set at pic-level b/c some pic-level initializations will be removed)
    uint8_t     enable_psad; // Enable pSad
    uint32_t    inter_depth_bias;
    uint8_t     bipred_available;
    uint8_t     is_intra_bordered;
    uint8_t     updated_enable_pme;
    Lpd1TxCtrls lpd1_tx_ctrls;
    uint8_t
        chroma_complexity; // Indicates which chroma components (if any) are complex, relative to luma. Chroma TX shortcuts based on luma should not be used when chroma is complex.
    uint8_t
        lpd1_skip_inter_tx_level; // Signal to skip INTER TX in LPD1; should only be used by M13 as this causes blocking artifacts.
    // 0: OFF, 1: Skip INTER TX if neighs have 0 coeffs, 2: skip all INTER TX
    COMPONENT_TYPE lpd1_chroma_comp; // chroma components to compensate at MDS3 of LPD1
    uint8_t        corrupted_mv_check;
    uint8_t        skip_pd0;
    uint8_t
        scale_palette; //   when MD is done on 8bit, scale  palette colors to 10bit (valid when bypass is 1)
    uint32_t b32_satd[4];
    uint8_t  high_freq_present;
} ModeDecisionContext;

typedef void (*EbAv1LambdaAssignFunc)(PictureControlSet *pcs_ptr, uint32_t *fast_lambda,
                                      uint32_t *full_lambda, uint8_t bit_depth, uint16_t qp_index,
                                      Bool multiply_lambda);

/**************************************
     * Extern Function Declarations
     **************************************/
extern EbErrorType mode_decision_context_ctor(
    ModeDecisionContext *context_ptr, EbColorFormat color_format, uint8_t sb_size, uint8_t enc_mode,
    uint16_t max_block_cnt, uint32_t encoder_bit_depth,
    EbFifo *mode_decision_configuration_input_fifo_ptr, EbFifo *mode_decision_output_fifo_ptr,
    uint8_t enable_hbd_mode_decision, uint8_t cfg_palette, uint32_t hierarchical_levels);

extern const EbAv1LambdaAssignFunc av1_lambda_assignment_function_table[4];

// Table that converts 0-63 Q-range values passed in outside to the Qindex
// range used internally.
static const uint8_t quantizer_to_qindex[] = {
    0,   4,   8,   12,  16,  20,  24,  28,  32,  36,  40,  44,  48,  52,  56,  60,
    64,  68,  72,  76,  80,  84,  88,  92,  96,  100, 104, 108, 112, 116, 120, 124,
    128, 132, 136, 140, 144, 148, 152, 156, 160, 164, 168, 172, 176, 180, 184, 188,
    192, 196, 200, 204, 208, 212, 216, 220, 224, 228, 232, 236, 240, 244, 249, 255};

#define FIXED_QP_OFFSET_COUNT 6
static const int percents[2][FIXED_QP_OFFSET_COUNT] = {
    {75, 70, 60, 20, 15, 0}, {76, 60, 30, 15, 8, 4} // libaom offsets
};
static const uint8_t uni_psy_bias[] = {
    85,  85,  85,  85,  85,  85,  85,  85,  85,  85,  85,  85,  85,  85,  85,  85,
    95,  95,  95,  95,  95,  95,  95,  95,  95,  95,  95,  95,  95,  95,  95,  95,
    95,  95,  95,  95,  95,  95,  95,  95,  95,  95,  95,  95,  95,  95,  95,  95,
    100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
};
extern void reset_mode_decision(SequenceControlSet *scs_ptr, ModeDecisionContext *context_ptr,
                                PictureControlSet *pcs_ptr, uint16_t tile_row_idx,
                                uint32_t segment_index);

extern void mode_decision_configure_sb(ModeDecisionContext *context_ptr, PictureControlSet *pcs_ptr,
                                       uint8_t sb_qp, uint8_t me_sb_qp);
extern void md_cfl_rd_pick_alpha(PictureControlSet           *pcs_ptr,
                                 ModeDecisionCandidateBuffer *candidate_buffer,
                                 ModeDecisionContext         *context_ptr,
                                 EbPictureBufferDesc         *input_picture_ptr,
                                 uint32_t                     input_cb_origin_in_index,
                                 uint32_t                     blk_chroma_origin_index);

#ifdef __cplusplus
}
#endif
#endif // EbModeDecisionProcess_h
