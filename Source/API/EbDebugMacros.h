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

#define FTR_STARTUP_MG_SIZE     1 // specify another mini-gop configuration for the first mini-gop after the key-frame
#define OPT_STARTUP_MG_SIZE     1 // Optimize startup mini-gop

#define OPT_LD                  1 // Optimize the performance of low delay mode

#if OPT_LD
#define OPT_LD_MRP              1 // Optimize MRP in low delay mode
#define OPT_LD_TF               1 // Optimize TF in low delay mode
#define OPT_LD_M13              1 // Optimize M13 in low delay mode
#define OPT_LD_M9               1 // Optimize M9 in low delay mode
#define OPT_LD_M10              1 // Optimize M10 in low delay mode
#define OPT_LD_M11              1 // Optimize M11 in low delay mode
#define OPT_LD_QPM              1 // Optimize QPM for low delay mode
#define OPT_LD_M12_13           1 // Optimize M12 and M13 in low delay mode
#define OPT_LD_LATENCY          1 // Optimize for low latency
#define OPT_LD_LATENCY_MD       1 // Optimize MD for low latency KF
#define OPT_LD_P2               1 // Optimize bdrate/speed trade off for low delay mode
#define OPT_LD_MRP2             1 // Optimize MRP in low delay mode by changing the references
#define OPT_CBR                 1 // Improve CBR by limiting the QP decrease between two base pictures
#endif
#define OPT_RPS_CONSTR             1 // Construct RPS in decode order; implement DPB at PD
#if OPT_RPS_CONSTR
#define OPT_RPS_CONSTR_2           1
#define OPT_RPS_CONSTR_3           1
#define FIX_INCOMP_MG_2            1 // Use 3L for incomplete MGs for all RA cases
#define CLN_REMOVE_REF_CNT         1 // Remove scs->reference_count
#define FIX_OVERLAY                1 // Fix overlay issue when OPT_RPS_CONSTR changes calling sequence
#define OPT_RPS_REFS               1 // Remove the duplicate references
#endif

#define FIX_LAYER_SIGNAL           1
#define FTR_PRED_STRUCT_CLASSIFIER 1
#define FTR_PRED_STRUCT_CLASSIFIER2 1

#define FIX_2009                   1 // fix for issue 2009, fixes mismatch between recon flag ON/OFF
#define EN_WARNING_FOR_MISMATCH    1 // create a new warning for mismatch that is expected for recon ON/OFF and stat-report ON/OFF for 10bit path
#define FIX_2042                   1 // fix for issue 2042, fixes corrupted bitstream when over boundary blocks are disabled

// LD Improvement
// applicable in all presets
#define OPT_LD_PD0                 1 // Optimize LPD0 for low delay mode
#define OPT_LD_MRP3                1 // Optimize MRP setting for low delay mode

#define OPT_LD2_M8                  1 // Optimize M8 in low delay mode
#define OPT_LD2_M9                  1 // Optimize M9 in low delay mode
#define OPT_LD2_M10                 1 // Optimize M10 in low delay mode
#define OPT_LD2_M11                 1 // Optimize M11 in low delay mode
#define OPT_LD2_M13                 1 // Optimize M13 in low delay mode
#define OPT_LD_TX_SHORT_CUT_OFF     1 // Fix TX short cut, enable bypassEncdec // to be tested in RA
#define OPT_LD_MRP5                 1 // Use 3 reference in list 0 --> to create a new level for ld only
#define OPT_LD_B_FIX                1 // Set the reference information of list 1 in low delay B mode
#define OPT_LD_PD1_1                1 // Optimize light PD1 setting in low delay mode
#define OPT_LD_DR                   1 // Optimize depth removal by adding the skip_pd0 option
#define OPT_LD_SKIPTX               1 // Optimize skip tx based on early skip estimation
#define OPT_LD_CDEF1                1 // Improve cdef for M12-13
#define OPT_LD_SKIP_TX_NEAREST_LPD1 1 // Enable lpd1_skip_inter_tx_level for low delay mode
#define OPT_LD_ME                   1 // loss less changes in ME, disable pre-HME for RTC M13
#define OPT_LD_RC3                  1 // Improve RC for RTC
#define OPT_LD_CAND_RED_LVL         1 // optimize candidate reduction controls for low-delay mode
#define OPT_LD_PALLET               1 // optimize pallet setting for low-delay mode
#define OPT_LD_SC_TF                1 // optimize TF setting for Screen content in the low-delay mode
#define OPT_LD_SPEED_M11_M12        1 // speed-up presets M11 and M12
#define OPT_LD2_SC_M12           1 // Optimize M12 for SC in low delay mode
#define OPT_LD_SC_PD0            1 // Optimize LPD0 for SC in low delay mode
#define OPT_LD_SC_HF             1 // Optimize High frequency setting for SC in low delay mode
#define OPT_LD_SC_ME             1 // Optimize PreHme/me settings for SC in low delay mode
#define OPT_LD2_SC_M11           1 // Optimize M11 for SC in low delay mode
#define OPT_LD2_SC_M10           1 // Optimize M10 for SC in low delay mode
#define OPT_LD2_SC_M9            1 // Optimize M9 for SC in low delay mode
#define OPT_LD2_SC_M8            1 // Optimize M8 for SC in low delay mode

#define OPT_LD_SC_MDS0           1 // Optimize intra cost syntax generation
#define OPT_LD_SC_RDOQ           1 // Disable RDOQ for RTC SC in M13

#define OPT_LD_CLEANUP           1 // clean up signal derivation functions with rtc checks
#define OPT_LD_CLEANUP_II       1 // Refactor rtc check to use block level signals (like fast_decode)

//FOR DEBUGGING - Do not remove
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

#ifdef __cplusplus
}
#endif // __cplusplus

// clang-format on

#endif // EbDebugMacros_h
