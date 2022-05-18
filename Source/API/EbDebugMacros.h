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


#define TUNE_DEFAULT_M10_M11    1 // Tune default M10 M11 with consideration to decode speed
#define TUNE_DEFAULT_M9         1 // Tune default M9 with consideration to decode speed
#define TUNE_DEFAULT_M8         1 // Tune default M8 with consideration to decode speed
#define OPT_MAX_P0_P1_NSQ       1 // Optimize max_part0_to_part1_dev NSQ-targeting feature
#define TUNE_DEFAULT_M7         1 // Tune default M7 with consideration to decode speed
#define TUNE_DEFAULT_M6         1 // Tune default M6 with consideration to decode speed
#define TUNE_DEFAULT_M5         1 // Tune default M5 with consideration to decode speed
#define TUNE_DEFAULT_M3         1 // Tune default M3 with consideration to decode speed

#define CBR_OPT                 1 // Optimisation of the cbr mode towards better coding efficiency and behavior
#define CBR_QPM                 1 // QP modulation for the cbr mode based on Motion estimatiom distorsion (no tpl)

#define TUNE_SSIM_M13           1 // Tune default M13 for SSIM
#define TUNE_SSIM_M12           1 // Tune default M12 for SSIM
#define TUNE_SSIM_M11           1 // Tune default M11 for SSIM
#define TUNE_SSIM_M8            1 // Tune default M8 for SSIM
#define TUNE_SSIM_M5            1 // Tune default M5 for SSIM
#define TUNE_SSIM_M2            1 // Tune default M2 for SSIM
#define TUNE_SSIM_M1            1 // Tune default M1 for SSIM

#define FIX_M9_M10_DLF          1 // Adjust M9 M10 DLF level to fix decoder speed
#define TUNE_M7_M8_DLF          1 // Adjust M7 M8 DLF level to improve decoder speed
#define FIX_DEC_SPEED_M6        1 // Fix M6 decoder speed

#define TUNE_INTER_COMPOUND     1 // Push down M3 Inter Compound Mode to M0
#define FIX_DISALLOW_8x8_SC     1 // Align sc and nsc settings for disallow 8x8 (keep NSC memory footprint reduction)
#define FG_LOSSLES_OPT          1 // Film Graing Lossless optimization
#define FIX_DELTAQ_BOUNDARY     1 // Fix boundary issue of av1_get_deltaq_offset()
                                  // When base_qindex is 255 and beta > 1, the original
                                  //    function will return invalid deltaq which leads to
                                  //    invalid q_index.
#define FIX_RDMULT_OVERFLOW     1 // Avoid int overflow in rdmult calculation
#define FIX_RATE_EST_SIGN       1 // Count sign bit cost in tpl rate cost

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
// Super-resolution debugging code
#define DEBUG_SCALING           0
#define DEBUG_TF                0
#define DEBUG_UPSCALING         0
#define DEBUG_SUPERRES_RECODE   0
#define DEBUG_SUPERRES_ENERGY   0
#define DEBUG_RC_CAP_LOG        0 // Prints for RC cap

// Switch frame debugging code
#define DEBUG_SFRAME            0

#ifdef __cplusplus
}
#endif // __cplusplus

// clang-format on

#endif // EbDebugMacros_h
