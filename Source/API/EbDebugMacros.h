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

// svt-10
#define OPT_DEPTHS_LVL                            1 // Optimize depth_level: levels and qp-banding
#define OPT_OBMC_DIST                             1 // Use scs-qp @ obmc pic early exit
#define OPT_TXS                                   1 // Optimize txs_level: levels and qp-banding
#define OPT_HIGH_FREQ                             1 // Opt high freq
#define OPT_NIC                                   1 // Optimize nic_level: levels and qp-banding
#define OPT_NEW_NSQ_LVLS                          1 // Add new NSQ levels and tune based on QP and coeff_level
#define OPT_TXT                                   1 // Optimize for txt levels: remove unnecessary levels and create new levels with better spacing in terms of speed
#define OPT_OBMC_LVLS                             1 // Optimize obmc_level
#define OPT_ME_HME_REF_PRUNE                      1 // Optimize levels for ME/HME ref pruning
#define OPT_WARP_LVLS                             1 // Optimize warp_level
#define OPT_CAND_RED                              1 // Optimize cand reduction level for M6
#define OPT_R0_NSQ                                1 // NSQ-offset = f(r0)
#define TUNE_SHUT_ME_CAND_PRUNE                   1 // Disable me cand pruning in M6

#define OPT_LIST1_REMOVAL                         1 // Optimize levels for List0-only @ Base

#define OPT_SUBPEL_LVLS                           1 // Optimize ME subpel levels

#define OPT_M7                                    1 // M7 tuning using M5 as a reference
#define OPT_DEPTH_REFINEMENT_LVLS                 1 // Use r0 instead of coeff-check/base-check @ depth-refinement lvl derivation
#define OPT_TPL_SUB_LVL                           1 // Optimize the tpl sub-lvl
#define OPT_DEPTH_REMOVAL                         1 // Optimize depth_removal
#define CLOSE_THE_SPEED_GAP_TO_1_8                1 // M5_cand_red + M5_txs
#define OPT_MED                                   1
#if OPT_MED
#define OPT_DIST_DEPTHS                           1 // med-banding for depths
#define OPT_DIST_BAND                             1 // med-banding for nsq
#define OPT_DIST_TXS                              1 // med-banding for depth-removal
#define OPT_STRENGTHEN_MED                        1 // conservative med
#endif
#define OPT_NSQ_CLASSIFIER_1                      1 // add/use sc_class3

#define TF_BRIGHTNESS_CHECK_V1                    1 // prune tf reference frame(s) using the brightness deviation to the central-frame
#define TF_HME_ME_                                1 // OPT HME of tf
#define OPT_SWAP_C1_C2                            1 // Swap CAND_CLASS_1 and CAND_CLASS_2 to process MVP candidates first
#define OPT_MFMV_FD                               1 // Use coeff_lvl for MFMV fast-decode levels instead of current ME_dist/QP
#define OPT_CDEF_0                                1 // Remove the fast-decode check from CDEF's lvl derivation
#define OPT_SG                                    1 // Tune SG settings
#define CLN_REMOVE_USELESS_CHECKS                 1 // Remove useless preset checks (lossless)
#define OPT_WM_BIS                                1 // Tune WM settings
#define CNL_INJ_NON_SIMPLE_MODES                  1 // Use common function to inject II/WM/OBMC
#define LAMBDA_SKIP_BIAS                          1 // Lambda tuning towards more skip selection
#define OPT_MFMV_FD_NEW                           1 // Use coeff_lvl for MFMV fast-decode levels instead of current ME_dist/QP
#define TUNE_M5                                   1 // Tune M5 features
#define FASTER_DEPTH_REMOVAL                      1 // Faster depth-removal
#define FASTER_NSQ_SEARCH                         1 // Faster nsq-search
#define FASTER_NSQ_GEOM                           1 // Faster nsq-geom
#define TUNE_M7                                   1 // Tune M7 features
#define SMOOTH_RIGHT                              1 // Smooth QP bands for high and medium QPs
#define SMOOTH_LEFT                               1 // Neutralize left qp-band for depth, txs, and nic
#define HIGHER_LAMBDA_EXTREME_RIGHT_BAND          1 // Neutralize left qp-band for depth, txs, and nic
#define TUNE_M8                                   1 // Tune M8 features
#define TUNE_M9                                   1 // Tune M9 features
#define TUNE_M10                                  1 // Tune M10 features
#define TUNE_M4                                   1 // Tune M4
#define TUNE_M0                                   1 // Tune M0 features
#define TUNE_M3                                   1 // Tune M3
#define TUNE_M11                                  1 // Tune M11
#define TUNE_M4_2                                 1 // Tune M4
#define TUNE_M13                                  1 // Tune M13
#define TUNE_SHIFT_MR_M0                          1 // Shift M0->M1, MR->M0, add new MR
#define TUNE_M0_M1                                1 // Tune M0, M1
#define FIX_TF_M13                                1 // Fix TF M13
#define FIX_LAMBDA                                1 // Fix lambda
#define USE_DEFAULT_M13_FOR_360p_480p             1 // Dont hard-code M12 for M13 360p/480p
#define FIX_RESIZE_R2R                            1 // Force SB 64x64 when resize mode is used to address long-standing r2r
#define ADD_LAMBDA_WEIGHT_SGN                     1 // Add a signal for lambda-weight
#define TUNE_M9_2                                 1 // Tune M9 features
#define TUNE_M10_2                                1 // Tune M10 features
#define TUNE_M8_2                                 1 // Tune M8 features
#define REMOVE_M13_WARN                           1 // Remove warnings for M13
#define TUNE_LIMIT_FD_M1                          1 // Limit fast-decode mode to M1-M10
#define OPT_LD                                    1 // Tune LD
#define CLN_REMOVE_UNUSED_FUNCS                   1 // Remove unused functions; removes a file

//FOR DEBUGGING - Do not remove
#define OPT_LD_LATENCY2         1 // Latency optimization for low delay - to keep the Macro for backwards testing until 3.0
#define OPT_FAST_DECODE_LVLS    0 // Add the ability to use different lvls for fast-deocde
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

// Variance boost debugging code
#define DEBUG_VAR_BOOST         0
#define DEBUG_VAR_BOOST_QP      0
#define DEBUG_VAR_BOOST_STATS   0

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
