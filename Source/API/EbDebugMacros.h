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

#define FIX_TOGGLE_MG              1 // Fix condition for updating lay0_toggle and lay1_toggle to account for case where MG start is not 0
#define FIX_OVERLAY_6L             1 // Set refresh_frame_mask to 0 for overlay frames
#define FIX_SHORT_KEYINT_WARN      1 // Fix warning when running short intra-period with 1-pass (implies the default configuration is 2-pass)
#define OPT_REPLACE_DEP_CNT        1 // Update the live count of PA refs at runtime based on the actual references, not init-time info
#define OPT_REPLACE_DEP_CNT_CL     1 // Update the live count of closed-loop refs at runtime based on the actual references, not init-time info
#define REMOVE_MANUAL_PRED         1 // Remove support for specifying a manual prediction structure
#define CLN_ENC_CTX                1 // Remove unused entries from encode_context
#define REMOVE_DEP_PIC_LIST        1 // Remove the list of dependent pics in the pred struct entries
#define CLN_PRED_STRUCT_CTOR       1 // Cleanup prediction_structure_ctor
#define OPT_PD_REF_QUEUE           1 // Refactor the picture decision PA ref queue
#define OPT_PM_REF_QUEUE           1 // Refactor the picture manager reference queue
#define OPT_TPL_REF_BUFFERS        1 // Create TPL reference buffers, instead of using closed loop reference objects, towards memory savings
#define CLN_PIC_DEC_PROC           1 // Cleanup picture_decision_kernel
#define CLN_PIC_MGR_PROC           1 // Cleanup picture_manager_kernel
#define FIX_REF_SIGN_BIAS          1 // Fix ref_frame_sign_bias derivation
#define FTR_USE_3L_AT_MG_END       1 // Use 3L for incomplete minigops
#define FTR_USE_3_BASE_REF         1 // Use 3 base pics as ref in 6L
#define FIX_OVERLAY_LA_BUFFS       1 // Make the overlay buffers dependent on the lookahead
#define CLN_TPL_LAD_MG_SET         1 // Add comment to clarify tpl_lad_mg settings
#define FIX_ME_PRUNING_R2R         1 // Prune r2r when using prune_me_candidates_th

#define FIX_DEFAULT_PRED_STRUCT    1 // Fix the default hierarchical levels

#define OPT_6L_TF                  1 // Increase the tf window size if 6L
#define FIX_INTRA_SELECTION        1 // Fixed the intra-percentage at reference frame(s) calculation; do not consider the stats if the reference frame is INTRA
#define TUNE_TPL_QPM_LAMBDA        1 // Use intra-selection at ref for QPS of L1&higher
                                     // Generate/use r0 for L1 of 6L
                                     // r0-based q modulation for L1 of 6L
                                     // Tune lambda modulation

#define OPT_NSQ                    1 // Optimize nsq for M1-M5
#define TUNE_NSQ                   1 // Tune M1, M2, M5, M6 for nsq optimization

#define TUNE_6L_M5                 1 // Tune M5 for 6L
#define TUNE_6L_M4                 1 // Tune M4 for 6L
#define TUNE_6L_M1                 1 // Tune M1 for 6L
#define TUNE_6L_M3                 1 // Tune M3 for 6L
#define OPT_WM                     1 // Optimize wm_level

#define TUNE_M1                    1 // Adjust M1 position
#define TUNE_M2                    1 // Adjust M2 position
#define TUNE_M3                    1 // Adjust M3 position
#define TUNE_M4                    1 // Adjust M4 position
#define TUNE_M5                    1 // Adjust M5 position
#define TUNE_M6_M7                 1 // Adjust M6/M7 position

#define FIX_BYPASS_ENCDEC          1 // Fix encdec condition
#define CLN_NSQ                    1 // Reorganize NSQ signals into 1 control
#if CLN_NSQ
#define ADD_NSQ_ENABLE             1 // Replace disallow_nsq
#endif

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

#ifdef __cplusplus
}
#endif // __cplusplus

// clang-format on

#endif // EbDebugMacros_h
