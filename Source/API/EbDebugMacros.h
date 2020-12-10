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
    // ============= END SVT_04 =============
//FOR DEBUGGING - Do not remove
#define NO_ENCDEC               0 // bypass encDec to test cmpliance of MD. complained achieved when skip_flag is OFF. Port sample code from VCI-SW_AV1_Candidate1 branch
#define DEBUG_TPL               0 // Prints to debug TPL
#define DETAILED_FRAME_OUTPUT   0 // Prints detailed frame output from the library for debugging
#define TUNE_CHROMA_SSIM        0 // Allows for Chroma and SSIM BDR-based Tuning
#define FTR_ENABLE_FIXED_QINDEX_OFFSETS 1

#define FIX_DDL                 1 // Fix deadlock issues
#ifdef __cplusplus
}
#endif // __cplusplus

#endif // EbDebugMacros_h
