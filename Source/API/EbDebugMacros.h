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
* - optimizations should be prefixed by OPT_
* - all macros must have a coherent comment explaining what the MACRO is doing
* - #if 0 / #if 1 are not to be used
*/

#ifndef EbDebugMacros_h
#define EbDebugMacros_h

// clang-format off

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

// svt-09 macros
#define OPT_NSQ_MV                                1 // Improve the refinement of the NSQ MV(s)
#define OPT_TF_PATH                               1 // Clean TF path
#define OPT_Q_PRUNE_TH_WEIGHT                     1 // Modulate mds-pruning th(s) using q and distortion
#define OPT_OBMC                                  1 // To do
#define OPT_MERGE_INTER_CANDS                     1 // Perform the inter-candidates merge @ block-basis using q and best PME/ME pred-error, instead of @ input-basis using coeff-level
#define OPT_TF_SP_EXIT                            1 // Apply TF subpel_early_exit feature to all bsizes with same TH for all
#define OPT_COEFF_LVL_NORM                        1 // Normalize the coeff_lvl feature
#define FIX_NSQ_CTRL                              1 // Break the nsq_ctrls into nsq_geom_ctrls and nsq_search_ctrls and use qp-banding for only the nsq_search_level derivation
#define CLN_SMALL_SIGS                            1 // Remove signals and levels having an insignificant impact on the behaviour
#define OPT_ME_SP_8TH_PEL                         1 // Optimize ME subpel
#define OPT_NSQ_HIGH_FREQ                         1 // Use more conservative NSQ settings in the presence of a high energy area 
#define OPT_1P_VBR                                1 // Optimized one-pass VBR
#define CLN_MVP_DIST_CALC                         1 // Move MVP distortion calc to one place to avoid recomputing
#define CLN_USE_BEST_PME_DIST                     1 // Use best PME dist in cand. reduction tool, instead of just best of ref_idx 0 cands
#define CLN_ADD_FIXED_PRED_SIG                    1 // Add signal to track if PDx has a fixed prediction structure
#define OPT_HME_L0_Q                              1 // Use Q @ the derivation of HME-L0 search-area (w,h)
#define CLN_NSQ_COPIES                            1 // Cleanup copying of neighbour arrays for NSQ shapes to avoid unnecessary copying
#define FIX_NSQ_SETTINGS                          1 // Fix when settings are reset in MD loop (settings can be modified by NSQ features)

//FOR DEBUGGING - Do not remove
#define OPT_LD_LATENCY2         0 // Latency optimization for low delay
#define LOG_ENC_DONE            0 // log encoder job one
#define NO_ENCDEC               0 // bypass encDec to test cmpliance of MD. complained achieved when skip_flag is OFF. Port sample code from VCI-SW_AV1_Candidate1 branch
#define DEBUG_TPL               0 // Prints to debug TPL
#define DETAILED_FRAME_OUTPUT   0 // Prints detailed frame output from the library for debugging
#define TUNE_CHROMA_SSIM        0 // Allows for Chroma and SSIM BDR-based Tuning
#define TUNE_CQP_CHROMA_SSIM    0 // Tune CQP qp scaling towards improved chroma and SSIM BDR

#define MIN_PIC_PARALLELIZATION 0 // Use the minimum amount of picture parallelization
#define SRM_REPORT              0 // Report SRM status
#define LAD_MG_PRINT            0 // Report LAD
#define RC_NO_R2R               0 // This is a debugging flag for RC and makes encoder to run with no R2R in RC mode
                                  // Note that the speed might impacted significantly
#if RC_NO_R2R
#define REMOVE_LP1_LPN_DIFF     1 // Disallow single-thread/multi-thread differences
#else
#define REMOVE_LP1_LPN_DIFF     0 // Disallow single-thread/multi-thread differences
#define FTR_KF_ON_FLY_SAMPLE      0 // Sample code to signal KF
#define FTR_RES_ON_FLY_SAMPLE     0 // Sample functions to change the resolution on the fly
#define FTR_RATE_ON_FLY_SAMPLE    0 // Sample functions to change bit rate
#endif
// Super-resolution debugging code
#define DEBUG_SCALING           0
#define DEBUG_TF                0
#define DEBUG_UPSCALING         0
#define DEBUG_SUPERRES_RECODE   0
#define DEBUG_SUPERRES_ENERGY   0
#define DEBUG_RC_CAP_LOG        0 // Prints for RC cap

// Switch frame debugging code
#define DEBUG_SFRAME            0
// Quantization matrices
#define DEBUG_QM_LEVEL          0
#define DEBUG_STARTUP_MG_SIZE   0
#define DEBUG_SEGMENT_QP        0
#define DEBUG_ROI               0
#ifdef __cplusplus
}
#endif // __cplusplus

// clang-format on

#endif // EbDebugMacros_h
