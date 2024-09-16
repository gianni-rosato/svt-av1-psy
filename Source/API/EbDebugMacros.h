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

#define OPT_ACC_PART_CTX_FEAT                     1 // Update the use_accurate_part_ctx feature to use the correct partition context
#define CLN_FD_USE_LOCAL                          1 // Use local variable for fast-decode to follow convention in signal derivation
#define UNIFY_DLF_LVL                             1 // Unified the DLF-level derivation for fd_3 and fd_1; used coeff_level-based modulation and skip_selection-based modulation for fd_3
#define FIX_PART_NEIGH_UPDATE                     1 // Fix partition neighbour updates
#define FIX_USE_SB_LAMBDA_DR                      1 // Use the SB lambda at depth-removal, instead of non-updated block-based lambda
#define OPT_FAST_DECODE_LVLS                      1 // Add the ability to use different lvls for fast-decode
#define OPT_FD3_CDEF                              1 // Optimize cost bias for cdef for fd3
#define OPT_SUBPEL                                1 // Bias towards fpel at the MD subpel-search: apply a penalty to the cost of fractional positions during the subpel-search each time we check against a full-pel MV
#define CLEAN_UP_FD_SIG                           1 // Clean-up fd signal(s); fd_2 to match fd_3, then remove fd_-1 and fd_3.
#define OPT_FD2_MFMV                              1 // Create a very-low coeff-level band (VLOW_LVL) and apply mfmv only to that band in fd2
#define TUNE_FD2                                  1 // Preset tuning for only fd2
#define TUNE_M4_M5_FD2                            1 // Tune fd2 for M4 and M5
#define TUNE_M5_M7                                1 // Tune M5 and M7 for fd0, fd1, and fd2
#define FIX_SKIP_INTRA_CHECK                      1 // Ensure intra is not skipped for islice
#define TUNE_M3_FD2                               1 // Tune M3 for fd2
#define TUNE_M2_FD2                               1 // Tune M2 for fd2
#define TUNE_M1_FD2                               1 // Tune M1 for fd2
#define ENABLE_FD_M0_MR                           1 // Enable fast decode modes for M0 and MR
#define TUNE_M0_MR_FD2                            1 // Tune M0 and MR for fd2
#define TUNE_M2_FD1                               1 // Tune M3 and below fd1
#define TUNE_M7_CDEF_FD1                          1 // Enable CDEF cost bias for FD1 in M7-M9
#define LIMIT_FD2_MR_M7                           1 // Limit fast-decode 2 to M7 and below presets
#define TUNE_LIMIT_FD                             1 // Enable fast-decode 2 up to M9; disable fast-decode in M10
#define CLN_UNUSED_CHECKS                         1 // Cleanup unused checks in the code
#define CLN_MRP_LVL                               1 // Cleanup how MRP level is set

//FOR DEBUGGING - Do not remove
#define OPT_LD_LATENCY2         1 // Latency optimization for low delay - to keep the Macro for backwards testing until 3.0
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
