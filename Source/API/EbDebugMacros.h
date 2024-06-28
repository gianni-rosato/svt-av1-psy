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

#define TUNE_M2                                   1 // Tune M0-M2
#define OPT_ME_ON_MV                              1 // Increase ME SA based on MV search centre
#define OPT_MV_DIFF_RATE                          1 // Include drl rate when choosing best MV pred
#define TUNE_COMP_LEVEL                           1 // Tune compound level for M1/2
#define OPT_NRST_NR_NEW                           1 // Enable NRST_NEW, etc. prediction types
#define TUNE_SB_SIZE_M2_RES                       1 // Use different QP th for 720p+ for SB size in M2
#define TUNE_UPDATE_CDF                           1 // Tune M2,3 level of update cdf
#define CLN_REMOVE_UNUSED_MACROS                  1 // remove unused/useless macros
#define CLN_REMOVE_UNUSED_SCS                     1 // Remove unused and useless signals under scs
#define OPT_LPD1                                  1
#if OPT_LPD1
#define OPTIMIZE_LPD1_PATH                        1 // Optimize the LPD1-settings
#define OPTIMIZE_LPD1_MDS0_DIST                   1 // Unify the mds0-distortion weighing between LPD1 and regular-PD1
#endif
#define OPT_LPD0_DEFS                             1 // Use seq-qp instead of pic-qp for LPD0-level modulation, and optimize LPD0 level-8
#define OPT_DLF                                   1 // Modulate the DLF-level using skip-information
#define TUNE_M5                                   1 // Tune M5
#define OPT_L0_ONLY_BASE                          1 // Use variance in list0_only_base classifier
#define TUNE_M5_2                                 1 // Tune M5
#define OPTIMIZE_LP_LEVELS                        1 // Optimize multithreaded settings
#define TUNE_M7                                   1 // Tune M7
#define SHUT_QP_SMOOTHING                         1 // Shut QP smoothing
#define OPT_DEPTH_LEVEL                           1 // Protect isolated cplx SB(s)
#define OPT_DLF_FD                                1 // Opt DLF for M7 fast-decode
#define TUNE_DR_M7_RES                            1 // Adopt 1080p depth-removal settings for 480p/720p in M7
#define TUNE_M5_TPL                               1 // Tune M5 tpl settings
#define TUNE_M0                                   1 // Tune M0
#define TUNE_MR                                   1 // Tune MR
#define TUNE_M1                                   1 // Tune M1
#define TUNE_M2_2                                 1 // Tune M2
#define OPT_DR                                    1 // Modulate depth-removal level using high-frequency
#define OPT_NSQ_PART0_PART1                       1 // Improve the modulation of the part0-to-part1 deviation-threshold: (1) the SQ-distortion to SQ-cost ratio:
                                                    // bypass part0-to-part1 if the SQ-distortion to SQ-cost ratio is lower than 50%, and
                                                    // (2) the SQ-pred_mode: apply a higher weight when pred-mode is NEW, and reduce the weight
                                                    // when pred-mode is Merge
#define OPT_NSQ_PSQ_CPLX                          1 // Apply a more aggressive psq_cplx_lvl to NSQ-level 16 and higher
#define TUNE_M3                                   1 // Tune M3
#define OPT_FILTER_INTRA                          1 // Optimize filter intra
#define TUNE_M4                                   1 // Tune M4
#define OPT_DEFAULT_M7                            1 // Tune M7
#define TUNE_DEFAULT_M5                           1 // Tune M5
#define OPT_DEFAULT_M8                            1 // Tune M8
#define TUNE_M13                                  1 // Tune M13
#define OPT_DEFAULT_M9                            1 // Tune M9
#define TUNE_M11                                  1 // Tune M11
#define OPT_DEFAULT_M10                           1 // Tune M10

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
