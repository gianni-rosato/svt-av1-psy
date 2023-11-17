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

#define FTR_KF_ON_FLY             1 // Add key frames on the fly
#define FTR_RES_ON_FLY            1 // Add Resolution on the fly
#define FTR_RES_ON_FLY2           1 // move parameters from scs to enc_ctx
#define FTR_RES_ON_FLY3           1 // Release scs_wrapper
#define FTR_RES_ON_FLY4           1 // changing Resolution update functions
#define FTR_RES_ON_FLY5           1 // changing Resolution update function calls
#define FTR_RES_ON_FLY6           1 // changing Resolution update signalling functions
                                    // limit resolution on the fly to supported modes
#define FTR_RES_ON_FLY7           1 // changing Resolution fix TF
#define FTR_RES_ON_FLY_APP        1 // Application changes for change the resolution on the fly
#define FTR_RATE_ON_FLY           1 // changing bit rate

#define CLN_ME_OFFSET                             1 // Create unified function for deriving me offset in MD
#define NEW_BLK_GEOM 1
#if NEW_BLK_GEOM
#define OPT_REORDER_GEOM                          1 // reorder geom to have HVA/HVB partitions at the end
#define OPT_REMOVE_HVAB_GEOM                      1 // Add new block geom with HVA/HVB shapes removed
#endif
#define OPT_REST_SIZE_FULL                        1 // Use the largest restoration unit size for 240p content
#define OPT_BLK_EARLY_EXIT                        1 // Exit early from process_block() based on known settings

#define OPT_VBR_PACKAGE         1 // VBR optimization
#if OPT_VBR_PACKAGE
#define OPT_VBR2                1 // remove IPPP and use lookahead data in VBR
#define OPT_VBR3                1 // lossless changes to save cycle for 1P VBR
#define OPT_VBR4                1 // lossless changes clean up
#define OPT_VBR_2P              1 // remove DG from 2P
#define OPT_VBR6                1 // VBR improvement
#define OPT_2PVBR               1 // 2 Pass VBR improvement
#endif
#define OPT_TF_LVL3                               1 // Opt tf lvl3

#define OPT_NSQ_RATE_THS                          1 // Apply offset to nsq rate THs for small block sizes
#define CLN_REMOVE_COND0                          1 // Remove cond0 from skip_sub_depth_ctrls
#define TUNE_NSQ_LEVELS                           1 // Tune NSQ settings
#define OPT_TPL                                   1 // Optimize TPL
#define TUNE_M8_M12                               1 // Tune M8-M12 for fast-decode-crf-lp1
#define OPT_MRP_VBR                               1 // Optimize MRP in VBR
#define TUNE_ENABLE_ME_8X8                        1 // tune enable_me_8x8 for 1080p

#define USE_PRED_MODE                             1 // Optimize the INTRA-search classifier of PD0, and use the pred-mode @ the sub-depth(s) refinement
#if USE_PRED_MODE
#define OPT_M4_BEYOND_DEPTHS                      0 // Not yet ready
#define OPT_M3_BEYOND_DEPTHS                      1 // Shut sub-depth refinement if INTRA pred-block (speed)
#define OPT_M2_BELOW_DEPTHS                       1 // (-2,2) for MR, and (-2,2) for INTRA pred-block for M0 (BD-rate)
#endif
#define OPT_PRE_MDS0_SEARCH                       1 // Variance instead of SSD @ the PME and NSQ-MV-refinement MR (BD-rate) and M0 (BD-rate)

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
