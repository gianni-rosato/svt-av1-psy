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
* - enabling a feature should be prefixed by ENABLE_
* - disableing a feature should be prefixed by DISABLE_
* - tuning a feature should be prefixed by TUNE_
* - adding a new feature should be prefixed by FEATURE_
* - bug fixes should be prefixed by FIX_
* - all macros must have a coherent comment explaining what the MACRO is doing
* - #if 0 / #if 1 are not to be used
*/


#ifndef EbDebugMacros_h
#define EbDebugMacros_h

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

// undefining this macro would allow the AVX512 optimization to be enabled by default
#ifndef NON_AVX512_SUPPORT
#define NON_AVX512_SUPPORT
#endif

#define FIX_RC_BUG 1 // Fix the one pass QP assignment using frames_to_be_encoded
#define FIX_VBR_BUG 1 // Fix 1st pass bug (bug from rebasing the branch)
#define FIX_10BIT     1 // fix 1st pass for 10bit input
#define FIX_RC_TOKEN     1 // fix RC token check to include double dash

#define FEATURE_NEW_DELAY             1 // I frames with GOP resetting (aka IDR) are delayed in Picture Decision to wait for buffering of the next mini-gop frames
#define FEATURE_INL_ME                1 //Enable in-loop ME
#if FEATURE_INL_ME
#define TUNE_IME_REUSE_TPL_RESULT     1 // Reuse TPL results for iLoopME
#define TUNE_INL_TPL_ENHANCEMENT      1 // Refinement for TPL
#define TUNE_INL_ME_RECON_INPUT       1 // Perform ME GM, TPL on input/recon: 1 on input, 0 on recon
#if !TUNE_IME_REUSE_TPL_RESULT
#define TUNE_SIGNAL_TPL_ME_OQ         1 // A separate signal_xxx_oq for TPL ME
#endif
#endif

#define FEATURE_PA_ME                1 // The ability to do ME inloop or in PAME

#define FEATURE_TPL_SOP              1  // Move TPL to source based operation for ime == 0

#define FEATURE_IN_LOOP_TPL 1 // Moving TPL to in loop
#if FEATURE_IN_LOOP_TPL
#define ENABLE_TPL_ZERO_LAD     1 // Enable TPL in loop to work with zero LAD
#define TUNE_TPL                1   // Tuning TPL algorithm for QP assignment
#define ENABLE_TPL_TRAILING     1 //enable trailing pictures for TPL
#define TUNE_TPL_OPT            1  // Tune TPL for better BDR/speed , added signals
#define TUNE_TPL_LOSSLESS       1  // Algorithmic  TPL clean up
#define TUNE_TPL_OIS            1 // move ois to inloop TPL, can be done in me kernel with scs_ptr->in_loop_ois = 0
#define TUNE_TPL_RATE           1 // remove  uncessary rate calculation
#define TUNE_TPL_FWD_FRAME      1 // Tune TPL for FWD FRAME
#define TUNE_TPL_END_OF_GOP     1 // Tune TPL for end of Gop frames
#endif
#define FIX_GM_COMPUTATION      1  // Fix global motion computation for different modes
#define FEATURE_MDS2 1 // TXT @ MDS2 if CLASS_0_3, and TXS/RDOQ @ MDS2 if CLASS_1_2
#define FEATURE_NIC_SCALING_PER_STAGE            1 // Add ability to scale NICs per stage; improve current trade-offs
#define TUNE_NICS                                1 // Tune settings for NIC scaling/pruning/# of stages to improve trade-offs with new scaling
#define FEATURE_PARTIAL_FREQUENCY                        1 //Calculate partial frequency transforms N2 or N4
#define TUNE_SC_QPS_IMP                          1 // Improve QP assignment for SC
#define FEATURE_REMOVE_CIRCULAR                  1 // Remove circular actions from current NSQ feautres; replace them with non-circular levels
#define FEATURE_NEW_INTER_COMP_LEVELS            1 // Add new levels and controls for inter compound; remove old levels
#define FEATURE_NEW_OBMC_LEVELS                  1 // Add new levels and controls for OBMC
#define TUNE_CDF                                 1 // Update CDF Levels
#define TUNE_TX_TYPE_LEVELS                      1 // Add Tx Type Grouping Levels
#define TUNE_REMOVE_UNUSED_NEIG_ARRAY            1 // Removes unused neighbor array
#define TUNE_INIT_BLOCK_OPT                      1 // optimize block initialization
#define TUNE_ME_IDX_LUPT                         1 // get index using lookuptable
#define FEATURE_INTER_INTRA_LEVELS               1 // Cleanup and modify inter-intra levels
#define TUNE_QPS_QPM                             1 // Improve the QPS settings for Keyframe. Improve QPM for nonI base frames
#define TUNE_CDEF_FILTER                         1 // Added new fast search for CDEF
#define FIX_OPTIMIZE_BUILD_QUANTIZER                 1 // Optimize eb_av1_build_quantizer():  called for each single frame (while the generated data does not change per frame). Moved buffer to sps, and performed 1 @ 1st frame only.
#define FIX_REMOVE_UNUSED_CODE                       1 // Remove unused code
#define FEATURE_OPT_IFS                              1 // Reverse IFS search order; regular to be performed last since the most probable then bypass the last evaluation if regular is the winner. 1 chroma compensation could be avoided if we add the ability to do chroma only when calling inter_comp.
#define FIX_BYPASS_USELESS_OPERATIONS                1 // Bypass useless operations when fast rate is OFF
#define FIX_USE_MDS_CNT_INIT                         1 // Use the actual number of candidates @ the init of cand_buff_indices
#define FIX_MOVE_PME_RES_INIT_UNDER_PME              1 // Bypass useless pme_res init
#define FIX_REMOVE_MD_SKIP_COEFF_CIRCUITERY          1 // Remove MD skip_coeff_context circuitery
#define FIX_REMOVE_MVP_MEMSET                        1 // Remove MVP generation useless memset()
#define FIX_OPT_FAST_COST_INIT                       1 // Use the actual number of md_stage_0 candidates @ fast_cost_array init
#define FIX_TUNIFY_SORTING_ARRAY                     1 // Unify MD sorting arrays into 1
#define FIX_IFS                                      1 // Fix IFS to use the actual motion_mode and the actual is_inter_intra
#define FEATURE_COST_BASED_PRED_REFINEMENT           1 // Add an offset to sub_to_current_th and parent_to_current_th on the cost range of the predicted block; use default ths for high cost(s) and more aggressive TH(s) for low cost(s)
#define FEATURE_PD0_CUT_DEPTH                        1 // Shut 16x16 & lower depth(s) based on the 64x64 distortion if sb_64x64
#define FEATURE_PD0_SHUT_SKIP_DC_SIGN_UPDATE         1 // Skip dc_sign derivation/update, and bypass useless mi_info updates
#define FEATURE_OPT_RDOQ                             1 // Use default_quant for chroma and rdoq_bypass = f(satd)
                                                       // lossless, early exit rdo, disables last md search tools (rdoq, txtype search, interpolation search)
#define FEATURE_OPT_TF                               1 // Add the ability to perform luma only @ tf, control tf_16x16 using tf_32x32 pred error, apply tf @ base only
#define TUNE_CFL_REF_ONLY                            1 // CFL only @ REF
#define FEATURE_GM_OPT                               1 // GM @ REF, bipred only, rotzoom model omly
#define TUNE_HME_ME_TUNING                           1 // HME/ME:HME_L1=8x3 instead of 16x16, HME_L2=8x3 instead of 16x16, MAX_ME=64x32 instead 64x64
#define FEATURE_DC_ONLY_AT_NON_REF                       1 // use only intra dc at non reference frame
#define TUNE_PALETTE_LEVEL                       1 // palette level will only be 6 for temporal layer == 0, not encode preset <=M3
#define FEATURE_MDS0_ELIMINATE_CAND                  1 // Eliminate candidates based on the estimated cost of the distortion in mds0.
#define TUNE_TPL_TOWARD_CHROMA                       1 //Tune TPL for better chroma. Only for 240P
#define TUNE_LOW_DELAY                               1 // Tuning the 0B, 1B and 3B settings to support mingop 1, 2 and 4


//FOR DEBUGGING - Do not remove
#define NO_ENCDEC         0 // bypass encDec to test cmpliance of MD. complained achieved when skip_flag is OFF. Port sample code from VCI-SW_AV1_Candidate1 branch
#define DEBUG_TPL         0 // Prints to debug TPL
#ifdef __cplusplus
}
#endif // __cplusplus

#endif // EbDebugMacros_h
