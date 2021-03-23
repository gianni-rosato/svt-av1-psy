/*
* Copyright(c) 2020 Intel Corporation
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/


/*
* This file contains only debug macros that are used during the development
* and are supposed to be cleaned up every tag cycle
* all macros must have the following format:
* - adding a new feature should be prefixed by FTR_
* - tuning a feature should be prefixed by TUNE_
* - enabling a feature should be prefixed by EN_
* - disabling a feature should be prefixed by DIS_
* - bug fixes should be prefixed by FIX_
* - code refactors should be prefixed by RFCTR_
* - code cleanups should be prefixed by CLN_
* - all macros must have a coherent comment explaining what the MACRO is doing
* - #if 0 / #if 1 are not to be used
*/


#ifndef EbDebugMacros_h
#define EbDebugMacros_h

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

    // ============= START SVT_04 =============
#define FIX_IME                       1 // Fix inloop me
#define FTR_TPL_TR                    1 // TPL trailing  support

#define RFCTR_MD_BLOCK_LOOP           1 // Refactor the loop that iterates over all blocks at MD
#define CLN_REMOVE_UNUSED_SIGNALS     1 // Remove mds3_intra_prune_th, and skip_cfl_cost_dev_th
#define FIX_Y_COEFF_FLAG_UPDATE       1 // Fix bug where y_has_coeff flag is overwritten by non-selected tx_types during tx_type_search
#define CLN_ENC_MODE_CHECK            1 // Make enc mode check conform to convention of using "<="
#define FTR_NEW_CYCLES_ALLOC          1 // Replace old cycles allocation with a cycles allocation algorithm that
                                        // does not depend on stats.  Merge cycles allocation and zero-sq-coeff feature.
#define TUNE_M8_TO_MATCH_M7                1
#define FTR_PD2_BLOCK_REDUCTION            1 // Reduce depth refinement based on the complexity of the SB.
#define FTR_REDUCE_MDS2_CAND               1 // Reduce mds3 candidates when mds0 and mds1 select the same best candidate.
#define FTR_DISABLE_ADAPTIVE_ME            1 // Disable adaptive ME
#define FTR_PD2_REDUCE_INTRA               1 // Reduce intra when the me distortion is low
#define CLN_CLEANUP_MDC_CTX                1 // Cleanup mdc context
#define FTR_USE_VAR_IN_FAST_LOOP           1 // Use var in fast loop
#define CLN_REMOVE_UNUSED_CODE             1 // Remove unused code related to nsq stat
#define FTR_PD2_REDUCE_MDS0                1 // Reduce the number of injected blocks
#define FTR_REDUCE_TXT_BASED_ON_DISTORTION 1 // Reduce the number of injected blocks
#define FTR_USE_VAR_IN_FAST_LOOP10BIT      1 // Use var in fast loop 10bit

#define FTR_EARLY_DEPTH_REMOVAL            1 // Create a dedicated ctrl for depth removal (i.e. independent from post-PD0 depth refinement), generate depth removal and PD0 depth refinement settings 1 time per SB (rather than per PD), use (me_distortion, variance) to reduce number of input depth(s) to PD0
#define FTR_NIC_PRUNING                    1 // Add nic_pruning_ctrls to group the loose the per mds pruning signals: mdsx_class_th, mdsx_cand_base_th .. , and to expose the hidden nic operations: mdsx_cand_sq_offset_th, mdsx_cand_intra_class_offset_th,..
#define CLN_FAST_COST                      1 // Clean-up fast-cost estimation kernels: remove  remove md_pass)
#define FTR_REF_BITS                       1 // Call estimate_ref_frame_type_bits() 1 time per block before md_stage_0() for all ref rather than per candidate.
#define FTR_BYPASS_IF_NO_FAST_RATE         1 // Bypass skip_flag_context and is_inter_ctx derivation if no fast rate estimation
#define CLN_MDC_CTX                        1 // Remove inter_pred_dir_neighbor_array, md_inter_pred_dir_neighbor_array, remove md_intra_chroma_mode_neighbor_array and and md_mv_neighbor_array
#define FTR_SCALE_FACTOR                   1 // Bypass scale_factor generation if not is_superres_none
#define CLN_INIT_OP                        1 // Bypass useless pme res and pme initializations
#define CLN_GET_LIST_GET_REF               1 // Use look-up tables instead of check(s) at get_list_idx(), get_ref_idx(), ..
#define FTR_EOB_CUL_LEVEL                  1 // Bypass useless eob memset(s) and early exit cul_level derivation if max COEFF_CONTEXT_MASK
#define FTR_TXT_SKIP_RATE_EST              1 // Do not perform rate estimation @ tx_type search if current tx_type dist is higher than best_cost

#define TUNE_PRESETS_CLEANUP               1 // Preset tuning for MRS, MR, M0 and M1

#define CLN_MD_CANDS                       1 // Cleanup MD candidate struct and candidate generation functions
#define CLN_NSQ_AND_STATS                  1 // Clean-up code for NSQ block processing and remove stats code
#define CLN_REMOVE_TX_BUFFS                1 // Remove unused TX buffers
#define RFCTR_INTRA_TX_INIT_FUNC           1 // refactor av1_get_tx_type()
#define CLN_REMOVE_ME_SSD_CALCS            1 // remove unused SSD calculation from ME
#define TUNE_M8_MAX_ME                     1 // Decrease the max ME area for M8
#define FTR_PF_PD0_DETECTOR                1 // Turn on PF in PD0 using a variance-based detector
#define FTR_ME_HME_PROTECT_CLOSEST_REF     1 // Protect the closest reference frames from ME/HME SAD-based pruning
#define TUNE_LOWER_PRESETS                 1 // Tune the lower presets M0, M1, M2, M3, and M4
#define TUNE_HBD_MODE_DECISION             1 // Push M0 level of HBD Mode Decision to M1

#define FTR_NEW_REF_PRUNING_CTRLS          1 // Replace old ref pruning features with new, unified controls; enabled MRP to M7
#define TUNE_HME_ME_SETTINGS               1 // Tune HME area for M1-M8 and create new levels of HME/ME based ref pruning
#define FTR_HME_REF_IDX_RESIZING           1 // Add ability to scale down HME size based on ref index
#define FTR_UPGRADE_COMP_LEVELS            1 // Add the ability to use distortion only and pred0-to-pred1 diff  @ inter-inter compound params generation,
                                             // add the ability to control inter-inter compound per elementary group, move MVP selection and MV tracking
                                             // before the compound_type loop as same MVs for all compound_type

#define FTR_CDEF_CHROMA_FOLLOWS_LUMA           1 // Bypass chroma V in the cdef search
#define FTR_GM_OPT_BASED_ON_ME                 1 // Bypass Global Motion on ME
#define FTR_OPTIMISE_TF                        1 // Bypass halfpel in temporal filtering

#define TUNE_M4_M8                             1 // Preset optimization: M4-M8

#define FTR_REDUCE_MDS3_COMPLEXITY             1 // Added 3 levels of tuning to reduce MDS3 complexity
#define TUNE_10BIT_MD_SETTINGS                 1 // New settings for hbd mode decision
#define TUNE_NEW_PRESETS_MR_M8                 1 // Tune presets MR to M8
#define TUNE_SB_SIZE                           1 // New SB size settings for M5 and M6
#define FTR_NSQ_RED_USING_RECON                1 // Use the deviation recon-to-src distortion to prune H/V
#define TUNE_IMPROVE_DEPTH_REFINEMENT          1 // Better M7 and better M8 depth refinement level(s)
#define TUNE_SC_SETTINGS                       1 // New SC settings for depth refinement, MRP, IBC, and depth removal

#define FTR_MODULATE_SRC_REC_TH                1 // Modulate the scr-to-recon TH
#define FTR_IMPROVE_DEPTH_REMOVAL              1 // Added stage-2 (pd2-based pd2 depth reduction), and more agressive levels for stage-1
#define TUNE_M9_LEVELS                         1 // Add M9 levels for stage1 depth reduction
#define TUNE_DEPTH_REMOVAL_PER_RESOLUTION      1 // Improve depth removal
#define TUNE_PRESETS_AND_PRUNING               1 // Shifted features in M4 and M5, and updated MD pruning levels for M6 and M9
#define TUNE_M9_FEATURES                       1 // Shifting features from M8 into M9 and remove I_SLICE check from disallow_4x4 in M9
#define FTR_M9_AGRESSIVE_LAST_MD_STAGE         1 // Agressive last md stage for M9
#define FTR_M9_AGRESSIVE_EARLY_CAN_ELIMINATION 1 // Agressive early candidate elimintation for M9
#define OPT_RDOQ_FOR_M9                        1 // RDOQ for M9
#define TUNE_M9_TF_LEVEL                       1 // create tf level 5 for m9 where temp_layer ==0 and res<=720p
#define TUNE_CFL_LEVEL_M8_M9                   1 // turn off cfl for non-base layer pictures for M8 and M9
#define TUNE_M9_OPT_DEPTH_REMOVAL              1 // Depth removal optimisation
#define TUNE_M9_TF_1080P                       1 // TF level for M9 for res=1080p
#define TUNE_M9_GM_INTER_COMPOUND              1 // Turn OFF GMV and inter compound mode for M9
#define TUNE_M9_ME_HME_TXT                     1 // Tune ME/HME and TxT for M9
#define TUNE_M9_PME                            1 // PME level for M9
#define FTR_M9_CDEF                            1 // CDEF optimisation
#define CLN_CANDIDATE_ELEMINATION_CTR          1 // Cleanup candidate elimination
#define TUNE_M9_ME_HME                         1 // Tune ME/HME for M9, making use of resolution checks to assign different values based on resolution
#define TUNE_M9_PF                             1 // create PF level for M9
#define OPT_TPL                                1 // Optimise tpl by using PF and SAD for intra_search
#define OPT_M9_TXT_PRED_DEPTH_PRUNING          1 // optimizations for m9 that include pruning, txt_level off for high-res, pred-depth
#define OPT_TX_TYPE_SEARCH                     1 // Early exit TXT search
#define TUNE_M9_SKIP_CTX_DC_SIGN               1 // Skip skip_ctx and dc_sign @ PD2:
#define TUNE_M9_IFS_SSE_ADAPT_ME_MV_NEAR_WM_TF 1 // Tune M9 for ifs, spatial sse, adaptive me, mv merge near cand, wm, temp f.
#define FTR_REDUCE_ME_INJECTION                1 // Reduce me candidates
#define TUNE_REMOVE_INTRA_STATS_TRACKING       1 // Remove INTRA stats tracking
#define TUNE_REMOVE_TXT_STATS                  1 // Remove TXT stats
#define TUNE_M9_TF_BASE                        1 // Remove temporal layer 1 checks in TF for 1080p only
#define FTR_SHUT_ENCDEC_MVP                    1 // Shut EncDec MVP
#define FTR_SHUT_ENCDEC_CBF_ZERO               1 // Remove EncDec Cbf Zero

#define TUNE_M9_DEPTH_REMOVAL                  1 // Add 32x32 vs 16x16 cost deviation for depth removal, apply to M9
#define FTR_TPL_REDUCE_NUMBER_OF_REF           1 // Reduce the number of references to search in tpl
#define TUNE_M7_M9                             1 // Adjust preset features for M7-M9
#define TUNE_M9_HME                            1 // Reduce HME search area for low-motion clips
#define TUNE_M9_GM_DETECTOR                    1 // Check the presence of stationnary block(s) per SB to avoid calling the complex GM detection for
                                                 // clip(s) with no global motion, and distance-based active_th
#define TUNE_FIX_TF                            1 // Remove noise_level feature from TF
#define TUNE_CDEF2                             0 // Optimize skip decision for CDEF filtering
#define FTR_USE_SKIP_MD                        1 // Use md skip decision to bypass residual generation, transform, quantization, inverse quant and
                                                 // inverse transform when it is skip
#define FTR_FAST_RATE_ESTIMATION               1 // Estimate the rate of the first (eob/N) coeff(s) and last coeff only
#define FTR_PRUNED_SUBPEL_TREE                 1 // Add pruned subpel tree function; this macro adds/removes files
#define TUNE_M6_FEATURES                       1 // Optimize M6 with features from M5 and M7
#define OPT_SB_CLASS                           1 // Remove SB class
#define OPT_REFINEMENT_SIGNALS                 1 // Remove useless block refinement signals
#define OPT_PME_RES_PREP                       1 // Optimize pme res
#define OPT_BLK_REFINEMENT_PREP                1 // Faster block refinement
#define FTR_SIMULATE_P_BASE                    1 // Use list0 only if BASE (mimik a P)
#define TUNE_M4_M5_DEC2                        1 // Tune for M4 and M5

#define TUNE_DEPTH_SKIP                        1 // Fix the problematic knob(s) of depth_skip
#define OPT_WM                                 1 // Optimize warp
#define OPT_LF                                 0 // Optimize lf
#define TUNE_TXT_M9                            1 // Tune txt for m9, add th_cost along with tuning coeff_th
#define TUNE_DEPTH_REMOVAL_M9                  1 // Tune depth_removal_level with base/nbase and new level 10
#define TUNE_M3_FEATURES                       1 // Optimize M3 with features from M2 and M4
#define TUNE_M2_FEATURES                       1 // Optimize M2 with features from M1 and M3
#define TUNE_M1_FEATURES                       1 // Optimize M1 with features from M0 and M2

#define FIX_DEPTH_REMOVAL_MODULATION           1 // Removed noise_level-based modulation,  added a 32x32 to 8x8 distortion deviation (besides the existing 32x32 to 16x16 distortion deviation check) before disallowing below 32x32, and tuning
#define FIX_DEPTH_REMOVAL_SETTINGS             1 // Upgrade depth_removal settings
#define TUNE_MR_M0_FEATURES                    1 // Optimize M0 with features from MR and M1, reduce M0 features by pushing them to MR
#define FTR_M10                                1 // Create M10
#define TUNE_M7_SC                             1 // Create a new ME area level for M7 SC
#define TUNE_M8_FEATURES                       1 // Optimize m8 features dec 15
#define TUNE_M5_FEATURES                       1 // Optimize m5 features dec 16
#define TUNE_M4_FEATURES                       1 // Optimize m4 featuers dec 16
#define TUNE_M6_M7_FEATURES                    1 // Optimize m6 and m7 features dec 16
#define CLN_NEAR_CTRLS                         1 // Add controsl for near nearest count
#define FTR_FINAL_M10                          1 // Finalize M10
#if FTR_FINAL_M10
#define TUNE_M10_TPL_TRAILING_FRAME_CNT        1
#define TUNE_M10_SHUT_NEAR_NEAR                1
#define TUNE_M10_USE_DC_ONLY                   1
#define TUNE_M10_P_BASE                        1
#define TUNE_M10_BYPASS_HME_LEVEL_1_2          1
#define TUNE_M10_NIC_PRUNING                   1
#define TUNE_M10_MD_EXIT                       1
#define TUNE_M10_MERGE_INTER_CLASSES           1
#define TUNE_M10_SUBPEL                        1
#endif

#define FTR_TPL_SEGMENTS                       1 // Add segments to TPL dispenser
#if FTR_TPL_SEGMENTS
#define TPL_KERNEL                  1 // Infrastrcture to add tpl dispenser kernel
#define TPL_ENABLE_TPL_KERNEL       1 // ENABLE tpl dispenser kernel
#define TPL_SEG                     1 // Add segments to tpl dispenser kernel
#endif

#define TUNE_MATCH_TR                          1 // Make trailing frames settings similar to non-trailing
#define TUNE_NEW_ME_HME                        1 // Incrase the ME size for multi-threaded modes
#define FIX_FE_CDF_UPDATE_CRASH_NBASE          1 // fix a crash when using frame_end_cdf_update for non-BASE pics only
#define DISABLE_FE_CDF_UPDATE_BASE             1 // use frame_end_cdf_update for non-BASE pics only for multi-threaded M9
#define TUNE_INTRA_PRED_MODE_MT                1 // Tune M9 setting of intra_pred_mode for multi-threaded M9
#define TUNE_M0_M3_BASE_NBASE                  1 // tune preset M0-M3 with base checks
#define DIS_TRAILING_PICTURES                  1 // turn off trailing pictures for all presets
#define TUNE_M4_BASE_NBASE                     1 // tune preset M4 with base checks

#define TUNE_TPL_USE_PRED_SAD                  1
#define TUNE_SUPER_BLOCK_SIZE_M4_M5            1 // change superblock size to 64x64 for select resolutions in M4 and M5
#define TUNE_UPDATE_CDF_LEVEL                  1 // tune update_cdf_level for presets M4-M8 using i-pic/base checks
#define FIX_R2R_10B_LAMBDA                     1 //tpl lambda/qp calc uses a  random hbd_mode_decision as not beeing set yet.

#define FTR_VBR_MT                             1 // Add supports for RC in multithreaded mode.
                                          // Move base_frame_target,this_frame_target and projected_frame_size to PCS
#if FTR_VBR_MT
#define FTR_VBR_MT_MINIGOP_FIX                 1 // Fix the non 5L minigops
#define FTR_VBR_MT_REMOVE_DEC_ORDER            1 // Remove the decode order constraint
#endif
#define TUNE_DEFAULT_RECODE_LOOP               1 // default recode-loop setting, reenc2 for M0~5 and reenc1 for M6~8
#define TUNE_FIRST_PASS_FRAME                  1 // make first_pass_frame segment based for improving cpu usage

#define TUNE_M0_REPOSITION                     1 // reposition M0 preset
#define TUNE_M1_REPOSITION                     1 // reposition M1 preset
#define TUNE_M3_REPOSITION                     1 // reposition M3 preset
#define TUNE_M4_REPOSITION                     1 // reposition M4 preset
#define TUNE_SHIFT_M2_M1                       1 // shift M2 differences into M1
#define TUNE_SHIFT_PRESETS_DOWN                1 // shift higher presets down all the way to M2
#define FIX_SC_MVP_TABLE_GEN                   1 // Fix an invalid input @ the call of generate_av1_mvp_table()
#define PR_1660                                1 // Fix R2R issue kernel mismatch in kernel svt_aom_highbd_10_variance16x4_avx2()

#define NEW_PRESETS                            1 // New preset changes, designed with fps considerations
#define USE_SB64_M3                            1 // Use SB 64x64 for M3+ in all resolutions
#define  FTR_LAD_MG                            1 // Add a look-ahead on top of the default re-ordering
#define  FTR_PRE_HME                           1 // preHME stage to catch vertical and horizontal high motion
#define  FTR_ALIGN_SC_DETECOR                  1 // Upgrade sc detector to three classes
#define CLN_MEM_REF                            1 // Clean-up non used src reference
#define FIX_KF_BOOST_CAP                       1 // Cap kf_boost
#define CLN_REST_FILTER                        1 // Remove fixed-cost operations when restoration filtering is on (not lossless)
#define FIX_ADD_MVP_MEMSET                     1 // Fix bitstream corruption for longer clips
#define FIX_COMPUTE_MEAN_8X8                   1 // Fix crash when block_mean_calc_prec is set to BLOCK_MEAN_PREC_FULL
#define FIX_UPDATE_DQPRESENT_FLAG              0 // Disable q_present flag when no qp offset is in ]-th, +th[ intervall
#define FIX_ADD_TPL_VALID                      1 // Set tpl to invalid state when mc_dep_cost_base is 0
#define FTR_BYPASS_RDOQ_CHROMA_QP_BASED        1 // Disable rdoq when qp offset is higher than certain th.
#define FTR_USE_LAD_TPL                        1 // Use extended tpl group

#define CLN_SB_DATA                            1 //  Memory optimizations in final_blk_arr (from svt-01)
#define CLN_DLF_RES_PROCESS                    1 //  Memory optimizations in DLF and restoration processes
#define CLN_MD_CAND_BUFF                       1 //  Memory optimizations in MD candidate buffers
#define FIX_TPL_4MG                            1 //  Fix the number of allocated pictures when extended tpl is used.

#define CLN_REMOVE_MEAN                        1 //  Remove mean calculation
#define TUNE_PICT_PARALLEL                     1 //  Tune picture parallelization
#define CLN_RES_PROCESS                        1 //  Memory optimizations in DLF and restoration processes
#define CLN_CTRL_INIT_MRP                      1 //  Control MRP settings
#define CLN_TPL_CONNECT_FLAG                   1 // Connect tpl_opt_flag to the remaining ofptimisation flags, lossless
#define PR_1650                                1 //Fix build AVX512
#if PR_1650
#define OPT_AVX512                             1 //Optimize some AVX2 and AVX512 transform, transpose, fdct16x16_avx512, fadst16x16_avx512, av1_fdct64_new_avx512, av1_fdct32_new_avx512, etc
#endif

#define CLN_PPCS                              1 // Remove 1first pass memory  allocation when first pass is not used
#define CLN_FA                                1 // Memory optimizations in MD final_blk_arr
#define CLN_REC                               1 // Update the padding to be function of super_block_size
#define CLN_BN                                1 // Move RestorationInfo from PPCS to CPCS
#define TUNE_PICT_PARALLEL_II                 1 // Tune picture parallelization phase 2
#define OPT_R0_FOR_LOW_MOTION                  1 // Adjust r0 for kf low motion clip

#define TUNE_M0_M8_MEGA_FEB                    1 // tuning for MEGA_FEB testset
#define TUNE_ME_M9_OPT                         1 // optimized me for faster m9
#define  FTR_NEW_RPS                          1 // New MRP struct keeping top layer as reference
#if FTR_NEW_RPS
//to support full MRP, need to enable all below macros
#define INC_TD       1 //lossles- needed to support new MRP structure
#define BASE_RPS     0
#define TOP_LAY_REF  0
#define ADD_16       0
#define TOPL           0
#define TOPL_PLUS    0
#define SIM_OLD_REF  1 //tune TPL to support new MRP referencing
#define LIMIT_TO_43  1 //limit the MRP changes to only 4,3 mode
#define IMP_4L       1 //improve 4L MRP
#endif




#define FTR_AOM_QPS_IF_TPL_OFF                 1 // Use identical QPS as aom when TPL is OFF

#define FIX_TXS_RATE_R2R    1 // fix R2R related tx size rate estimation
#define FIX_SCD                     1 // remove use of scd_delay in TPL
#define TUNE_NEAR_CTRLS                       1 // Tune near settings
#define FIX_USELESS_ENCDEC_CYCLE              1 // Remove useless EncDec operations
#define TUNE_M8_FAST                          1 // tune me, hme, set-me-sr-adjust for m8

#define      CLN_STRUCT                      1 // create new structure for enc-dec for pixel and coeffs
#define      CLN_REDUCE_TPL_RECON            1 // No recon needed in some levels of TPL
#define      CLN_DLF                         1 // No recon needed when frame DFL is OFF
#define      CLN_REST                        1 // No restoration buffers when it is OFF
#define      CLN_OIS                         1 // Move to in loop OIS (OIS is now done in TPL)

#define TUNE_TPL_VBR                          1 // Tune VBR based on the latest tpl changes
#define DIS_2PASS_CRF                         1 // Disable 2pass CRF
#define TUNE_VBR                              1 // Tune VBR to limit qp decrease in non base layer. Update the default rate tables
#define VARIANCE_AVX512                       1 // svt_aom_varianceWxH_avx512 for sizes 32x64, 32x32, 32x16, 32x8, 64x16, 64x32, 64x64, 64x128, 128x64, 128x128
#define LOG2_FN_PTR                           1 //Replace log2f_32() with svt_log2f()
#define OPT_IFS_FOR_SKIP_MODE                 1 // Check IFS for skip_mode candidates; if IFS selected, disables skip_mode; rename merge_flag to skip_mode_allowed (no macro)
#define FIX_UPDATE_COEFF_FLAGS                1 // If the cost of not sending coeffs is favourable, update the coeff flags to reflect the decision
#define FIX_CHROMA_QUANT_MODE                 1 // Fix the prediction mode passed to the quant/inv. quant function for chroma


#define OPT3_DECIMATION                       1 // Unify the use of the downsampling method whether decimation or filtering.
#define OPT1_REMOVE_FLAT_NOISE                1 // remove flat noise
#define OPT5_TPL_REDUCE_PAD                   1 // reduce pad in tpl search
#define OPT6_DEPTH_REFINEMENT                 1 // reduce unnecessary operations when nsq/4x4 are off
#define OPT7_BUILD_CAND_BLK                   1 // reduce unnecessary operations when nsq/4x4 are off
#define OPT10_MI_MAP                          1 // cleanup of mi update


#define TUNE_REDESIGN_TF_CTRLS                 1 // Improve TF ctrls, and TF performance

#define FIX_GM_R2R                             1 // Fix GM initialization

#if FTR_PRUNED_SUBPEL_TREE
#define SUB_PIXEL_VAR_AVX512                  1 // svt_aom_sub_pixel_varianceWxH_avx512 for sizes 32x16, 32x32, 32x64, 64x32, 64x64, 64x128, 128x64, 128x128
#endif
#define TUNE_UPDATE_SCD_DELAY                  1 // Update the scd_delay based on the the number of future frames @ BASE

#define  PIC_BASED_MFMV                       1 // Suport for picture based MFMV. Lossless.
#define RFCTR_MFMV                            1 // Lossless MFMV refactoring
#if RFCTR_MFMV
#define  OPT_MFMV       1
#define  OPT_MFMV_2     1
#define  OPT_MFMV_3     1
#define  OPT_MFMV_4     1
#endif
#define  PIC_BASED_MFMV_R0                    0 //Use r0  to control MFMV
#define  PIC_BASED_MFMV_SAD                   0 //Use sad to control MFMV

#define OPT_BUILD_CAND_BLK_2                  1 // Lossless optimization of build_cand_block_array()
#define OPT_MI_UPDATE                         1 // Lossless opt for updating mode info
#define OPT12_PD0                             1 // Bypass neighbor update when intra is not allowed in PD0
#define OPT13_PD0                             1 // Bypass all mds1, mds2 computations in PD0
#define OPT_INLINE_FILTER_FUNCS               1 // Make funcs inline to improve efficiency
#define CLN_OLD_RC                            1 // Clean up old useless rate control code
#define FIX_SCD_DELAY                         1 // Fix the run2run issue caused by scd_delay
// ============= END SVT_04 =============
//FOR DEBUGGING - Do not remove
#define NO_ENCDEC               0 // bypass encDec to test cmpliance of MD. complained achieved when skip_flag is OFF. Port sample code from VCI-SW_AV1_Candidate1 branch
#define DEBUG_TPL               0 // Prints to debug TPL
#define DETAILED_FRAME_OUTPUT   0 // Prints detailed frame output from the library for debugging
#define TUNE_CHROMA_SSIM        0 // Allows for Chroma and SSIM BDR-based Tuning
#define FTR_ENABLE_FIXED_QINDEX_OFFSETS 1

#define FIX_DDL                 1 // Fix deadlock issues
#define MIN_PIC_PARALLELIZATION 0 // Use the minimum amount of picture parallelization
#define  SRM_REPORT             0 // Report SRM status
#define  LAD_MG_PRINT             0 // Report LAD
#ifdef __cplusplus
}
#endif // __cplusplus

#endif // EbDebugMacros_h
