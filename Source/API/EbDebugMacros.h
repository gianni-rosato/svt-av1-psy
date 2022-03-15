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

// clang-format off

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus
#define FTR_CBR                 1 // Adding CBR mode for LD mode
#define FIX_TF_LDB              1 // Fix TF parameters for LD
#define FTR_LDB_MT              1 // New Multi-Threading parameters for LD

#define VMAF_OPT                1 // Optimizing Y-PSNR/Y-SSIM/VMAF
#define TUNE_M13                1 // Optimizing M13 for better slope
#define TUNE_4L_M11             1 // Optimizing 4L compared to 5L for M11
#define OPT_REMOVE_TL_CHECKS    1 // Remove temporal layer checks in M7-M13, more removal in M6-M7
#define TUNE_4L_M10             1 // Optimizing 4L compared to 5L for M10
#define TUNE_4L_M9              1 // Optimizing 4L compared to 5L for M9
#define TUNE_4L_M8              1 // Optimizing 4L compared to 5L for M8
#define TUNE_4L_M7              1 // Optimizing 4L compared to 5L for M7

#define CLN_MD_CTX              1 // Cleanup MD ctx fields and memory alloc
#define OPT_UPDATE_CDF_MEM      1 // Only allocate memory for update_cdf rate estimation table if feature is on
#define OPT_CAND_BUFF_MEM       1 // Don't alloc memory for independent chroma search candidates when not used
#define OPT_MV_INJ_CHECK        1 // Optimize memory structures for tracking previously injected MVs
#define FIX_LPD1_W_OBMC         1 // Enable LPD1 to be used when OBMC is on at the frame level (previously caused corrupted bitstream)
#define OPT_LPD_LVL2            1 // Use a more conservative cost_th_dist
#define CLN_DLF_MEM_ALLOC       1 // Remove unnecessary check for DLF level during PCS memory allocation
#define CLN_SIG_DERIV           1 // Cleanup signal derivation functions
#define TUNE_M7                 1 // Slowdown M7 preset to improve spacing
#define TUNE_4L_M6              1 // Optimizing 4L compared to 5L for M6
#define TUNE_4L_M5              1 // Optimizing 4L compared to 5L for M5
#define TUNE_4L_M4              1 // Optimizing 4L compared to 5L for M4
#define TUNE_4L_M3              1 // Optimizing 4L compared to 5L for M3
#define OPT_VQ_MODE             1 // Improve VMAF while preserving the VQ gain
#define OPT_M10_SUBJ            1 // Speedup M10 for subjective mode
#define OPT_M8_SUBJ             1 // Speedup M8 for subjective mode
#define OPT_M7_SUBJ             1 // Speedup M7 for subjective mode

#define FIX_AQ_MODE             1 // Fix aq-mode 0
// RC refactoring
#define FRFCTR_RC_P1            1 // clean up extra variables
#define FRFCTR_RC_P2            1 // clean up overlay signals and num_stats_required_for_gfu_boost
#define FRFCTR_RC_P3            1 //remove KeyFrameCfg
#define FRFCTR_RC_P4            1 //remove GFConfig
#define FRFCTR_RC_P5            1 //remove rc min max interval
#define FRFCTR_RC_P6            1 //CBR and flat VBR clean up
#define FRFCTR_RC_P8            1 //remove EncodeFrameParams and CurrentFrame
#define FRFCTR_RC_P9            1 //clean up rc variables
#define FIX_TUNE_0              1
#if FIX_TUNE_0
#define FIX_VQ_PRED_STRUCT      1 // Undo the forcing of 4L if VQ mode
#define FIX_VQ_MODE_TF          1 // Ignore the rate component at RDOQ for negative delta - QP towards quant - coeff preservation
#define FIX_VQ_MODE_RDOQ        1 // Use a more conservative input variance threshold @ TF towards a weaker level selection
#define FIX_VQ_MODE_TPL         1 // Treat the 1st BASE after a transition as an I_SLICE at the delta - QP derivation(aggressive action when lowering the QP)
#endif


#define CLN_MD_CAND             1 // Cleanup MD candidate struct
#if CLN_MD_CAND
#define CLN_CAND_MV             1 // Cleanup candidate MV arrays
#define CLN_MOVE_COSTS          1 // Move info about candidate results (e.g. cost, distortion, etc.) to candidate_buffer
#define CLN_REMOVE_REDUND       1 // Remove redundant fields that can be derived from other fields (is_compound)
#define CLN_REMOVE_REDUND_2     1 // Remove redundant fields that can be derived from other fields (prediction_direction)
#define CLN_REMOVE_REDUND_3     1 // Remove redundant fields that can be derived from other fields (intra_luma_mode)
#define CLN_REMOVE_REDUND_4     1 // Remove redundant fields that can be derived from other fields (type)
#define CLN_REMOVE_REDUND_5     1 // Remove redundant fields that can be derived from other fields (local_warp_valid)
#define CLN_MOVE_COSTS_2        1 // Move info about candidate results (e.g. eob, quantized_dc, etc.) to candidate_buffer
#define CLN_REMOVE_REDUND_6     1 // Remove redundant fields that can be derived from other fields (is_direction_mode_flag)
#define CLN_CAND_TYPES          1 // Reorder entries in MD candidate struct and reduce size of some fields
#endif
#define CLN_SKIP_NAMING         1 // Rename skip_mode/skip coeff variables for clarity
#define CLN_DEFINITIONS         1 // Cleanup EbDefinitions.h
#define CLN_IS_VALID_MV_DIFF    1 // Cleanup arguments of is_valid_mv_diff()
#define TUNE_FAST_DECODE        1 // Tune existing fast-decode levels and add support for --fast-decode in M0-M4

//FOR DEBUGGING - Do not remove
#define LOG_ENC_DONE            0 // log encoder job one
#define NO_ENCDEC               0 // bypass encDec to test cmpliance of MD. complained achieved when skip_flag is OFF. Port sample code from VCI-SW_AV1_Candidate1 branch
#define DEBUG_TPL               0 // Prints to debug TPL
#define DETAILED_FRAME_OUTPUT   0 // Prints detailed frame output from the library for debugging
#define TUNE_CHROMA_SSIM        0 // Allows for Chroma and SSIM BDR-based Tuning

#define MIN_PIC_PARALLELIZATION 0 // Use the minimum amount of picture parallelization
#define SRM_REPORT              0 // Report SRM status
#define LAD_MG_PRINT            0 // Report LAD
#define RC_NO_R2R               0 // This is a debugging flag for RC and makes encoder to run with no R2R in RC mode
                                  // Note that the speed might impacted significantly
#if RC_NO_R2R
#define REMOVE_LP1_LPN_DIFF     1 // Disallow single-thread/multi-thread differences
#else
#define REMOVE_LP1_LPN_DIFF     0 // Disallow single-thread/multi-thread differences
#endif
#define FIX_1PVBR               1 // Derive initial qp based on target bitrate
// Super-resolution debugging code
#define DEBUG_SCALING           0
#define DEBUG_TF                0
#define DEBUG_UPSCALING         0
#define DEBUG_SUPERRES_RECODE   0
#define DEBUG_SUPERRES_ENERGY   0
#define DEBUG_RC_CAP_LOG        0 // Prints for RC cap

#ifdef __cplusplus
}
#endif // __cplusplus

// clang-format on

#endif // EbDebugMacros_h
