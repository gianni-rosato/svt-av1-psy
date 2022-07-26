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
#define OPT_TXS                 1 // Cleanup TXS signalling and improve controls
#define OPT_TPL_QPS             1 // Optimize QPS when using TPL
#define FTR_TPL_SUBPEL          1 // Add subpel search to TPL path

#define FTR_USE_COEFF_LVL       1 // Use me_8x8_distortion and QP to predict the coeff level per frame
#define FTR_NO_TILE_FAST_DEC    1 // Disable use of tiles by default when fast-decode is used

#define FTR_RC_VBR_IMR          1 // optimize TPL for VBR, 1 Pass VBR improvement
#define TUNE_M1_M3_BDR          1 // Tune m1 to m3 to recover maximum bdr
#define CLN_TPL_OPT             1 // Remove tpl_opt_flag that acts as an enable switch for many TPL features, despite having their own signals
#define FTR_USE_TPL_INTRA       1 // Use TPL results to restrict tested INTRA modes in MD
#define FIX_GMV_LVL4            1 // GMV - Shut unipred-ref restriction (params generation)
#define FIX_GMV_DOWN            1 // GMV - Add the ability to specify the down-sampling method per input (instead of per preset), then use the average ME distortion to switch between Full and 1/4
#define TUNE_M1_M3              1 // Tune M1 and M3 features
#define TUNE_M4                 1 // Tune M4 features for SSIM
#define TUNE_M3_NSQ             1 // Use HV4 blocks in M3 (all frames).  Adopt new shortcuts to offset the speed cost.
#define TUNE_M5                 1 // Tune M5 features for SSIM
#define FTR_OPTIMIZE_OBMC       1 // Bypass the OBMC-refinement (both fullpel and subpel), and use the regular MV(s)  (instead of shutting OBMC) for 32x32&above block(s))
#define TUNE_M3                 1 // Tune m3 md inter intra level
#define TUNE_M5_MRP             1 // Tune m5 MRP level
#define TUNE_CHROMA_RDOQ        1 // Increase the lambda used for RDOQ of chroma towards rate savings
#define OPT_QPS_WEIGHT          1 // Optimize the weighting used in QPS to improve VMAF in low presets, SSIM in high presets
#define FIX_REST_SANITIZER      1 // Fix sanitizer failure in restoration
#define FIX_SCALE_SANITIZER     1 // Fix sanitizer failure in restoration when scaled
#define FIX_UV_QINDEX_OFFSET    1 // Make av1_quantize_inv_quantize() chroma qindex aware
                                  // Decouple use-fixed-qindex-offsets and chroma qp offset for non-RC configurations(s)

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
