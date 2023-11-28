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

// svt-08-temp macros
#define CLN_NSQ_LVLS                              1 // cleanup NSQ levels
#define OPT_LPD0                                  1 // optimize LPD0 levels, lpd0 enc_dec_deriv, and lpd0 detector thresholds
#define CLN_ME_OFFSET                             1 // Create unified function for deriving me offset in MD
#define NEW_BLK_GEOM                              1
#if NEW_BLK_GEOM
#define OPT_REORDER_GEOM                          1 // reorder geom to have HVA/HVB partitions at the end
#define OPT_REMOVE_HVAB_GEOM                      1 // Add new block geom with HVA/HVB shapes removed
#endif
#define OPT_REST_SIZE_FULL                        1 // Use the largest restoration unit size for 240p content
#define OPT_BLK_EARLY_EXIT                        1 // Exit early from process_block() based on known settings
#define OPT_VBR_PACKAGE                           1 // VBR optimization
#if OPT_VBR_PACKAGE
#define OPT_VBR2                                  1 // remove IPPP and use lookahead data in VBR
#define OPT_VBR3                                  1 // lossless changes to save cycle for 1P VBR
#define OPT_VBR4                                  1 // lossless changes clean up
#define OPT_VBR_2P                                1 // remove DG from 2P
#define OPT_VBR6                                  1 // VBR improvement
#define OPT_2PVBR                                 1 // 2 Pass VBR improvement
#define CLN_VBR                                   1 // VBR Clean up
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
#define OPT_M3_BEYOND_DEPTHS                      1 // Shut sub-depth refinement if INTRA pred-block (speed)
#define OPT_M2_BELOW_DEPTHS                       1 // (-2,2) for MR, and (-2,2) for INTRA pred-block for M0 (BD-rate)#endif
#endif
#define OPT_PRE_MDS0_SEARCH                       1 // Variance instead of SSD @ the PME and NSQ-MV-refinement MR (BD-rate) and M0 (BD-rate)
#define OPT_PD0_NSQ                               1 // Optimize the NSQ of PD 0 MR (speed) and M0 (BD-rate)
#define TUNE_M0_M1_M3_M5                          1 // Optimize m0, m1, m3, m5
#define OPT_UNIPRED_3x3                           1 // Opt unipred_3x3 for MR and M0
#define OPT_PME                                   1 // Opt PME for MR and M0
#define TUNE_M7_M13                               1 // Tune m7-m13 for lp8
#define CLN_MISC                                  1 // Clean redundancies
#define OPT_UV_RDOQ                               1 // Shut early uv quant-coeff shaving
#define OPT_EOB_RDOQ                              1 // Improve RDOQ EOB check
#define OPT_TF_FACTOR_LARGE_BLOCKS                1 // Double the decay factor when large blocks are selected
#define OPT_TF_8X8_BLOCKS                         1 // Include 8x8 blocks in TF motion search and prediction
#define OPT_TF_AHD_TH                             1 // Adjust the THs used to prune frames in TF based on AHD
#define FIX_TF_64X64_PRED                         1 // Fix improper block size in TF 64x64 prediction path
#define OPT_MR_M0                                 1 // optimize mr and m0 features
#define OPT_M1                                    1 // tune m1 features
#define OPT_INTRA_LVL                             1 // optimize intra levels for higher presets
#define OPT_WARP_REFINEMENT                       1  // Opt Warp refinment for MR/M0
#define OPT_OBMC_REFINEMENT                       1  // Opt OBMC refinment for MR/M0/M1
#define OPT_FAST_LAMBDA                           1  // Use full_lambda modulation process for fast_lambda derivation
#define FIX_ENC_MODE_TF                           1  // From unsigned to signed
#define TUNE_M2_M3_M4                             1 // tune m2-m4 for fast decode lp8
#define TUNE_M5_M6                                1 // tune m5-m6 for fast decode lp8
#define OPT_ME_8x8                                1 // optimize ME SA based on 8x8 var cost
#define OPT_M2                                    1 // optimize m2 features
#define TUNE_M9                                   1 // tune M9
#define OPT_BIPRED3x3                             1 // Use unipred info to cut the number of search position(s)
#define OPT_HP_MV                                 1 // Optimize hp; use the percentage of hp at ref frame(s)
#define OPT_M3                                    1 // optimize m3 features
#define OPT_GEOM_SB12B_B4                         1 // 128x128 --> 8x8 (H/V + H4/V4 only)
#define DIS_LPD0_128x128                          1 // Disable lpd0 when sb size is 128
#define OPT_SB_SIZE                               1 // optimize sb size feature
#define OPT_M4                                    1 // optimize m4 features
#define OPT_M5                                    1 // optimize m5 features
#define OPT_NSQ_LEVELS                            1 // optimize nsq levels
#define CLN_CMPOUND                               1 //Lossless clean up for compound - avoid extra useless computation
#define OPT_CMPOUND                               1 //Optimize compound level
#define CLN_IS_REF                                1 // cleanup is_ref flags
#define OPT_M0_M1                                 1 // opt m0 m1 features
#define USE_QP_COMP                               1 // Modulate the accuracy of the prediction tools using QP
#if USE_QP_COMP
#define OPT_NSQ_QP                                1
#define OPT_NIC_QP                                1
#define OPT_DEPTH_QP                              1 // temp m5-m7
#define OPT_COEFF_LVL                             1
#define OPT_LPD0_QP                               1
#define OPT_TX_BYPASS                             1
#define OPT_DLF_QP                                1
#endif
#define FIX_NSQ_MEM_ALLOC                         1 // Fix memory allocation for NSQ
#define OPT_SC_PAL_DEPTH                          1 // improve SC trade-offs for M3-M12, changes to depth_level and palette_level
#define TUNE_MRP                                  1 // Tune mrp m6-m7
#define OPT_Q_TXT                                 1 // Modulate txt-th(s) using q
#define OPT_ME_RP                                 1 // opt ME ref pruning
#define OPT_DEPTH_LEVELS                          1 // Optimize depth level
#define OPT_Q_PME                                 1 // Modulate PME (w,h) using q
#define OPT_DEPTH_REFINE                          1 // optimize block based depth refinement
#define TUNE_M0_MR                                1 // tune m0 / mr features
#define TUNE_TXS_LEVELS                           1 // tune TXS feature levels
#define TUNE_DEPTH_REMOVAL                        1 // tune fast decode depth removal
#define OPT_M1M2                                  1 // Opt m1 and m2
#define OPT_DEPTH_REFIN                           1 // Use split rate to skip testing parent depth in PD1
#define TUNE_M6_II                                1 // tune M6
#define OPT_TX_SIZE                               1 // optimize TXS feature levels
#define OPT_Q_WARP                                1 // Ctrl Warp using q and distortion
#define OPT_Q_GLOBAL                              1 // Ctrl GM using q
#define OPT_Q_TF                                  1 // Ctrl TF using q
#if OPT_Q_TF
#define OPT_Q_OFFSET_TF                           1
#define OPT_Q_ME_TF                               1
#define PPT_TF_PRESET                             1
#endif
#define OPT_Q_OBMC                                1 // Ctrl OBMC using q and distortion
#define TUNE_M1_M4                                1 // tune m1-m4 fast decode
#define OPT_FI                                    1 // Filter intra opt for M4
#define TUNE_M7_II                                1 // Tune M7
#define NEW_M0                                    1 // New preset candidate (between M0 / MR). Macro on + run MR for the new preset
#define TUNE_M5_II                                1 // tune m5 fast decode
#define OPT_CHILD_DEPTH_RATE                      1 // Use split rate info to skip child depth in PD1
#define TUNE_DEPTH_LVLS                           1 // Tune depth levels for M8-M10
#define TUNE_ILV                                  1 // Tune Intra level
#define OPT_DR_QP                                 1 // Use QP to modulate depth refinement thresholds
#define TUNE_M4_II                                1 // Tune m4
#define TUNE_M8_II                                1 // Tune M8
#define TUNE_CDEF_M56                             1 // Tune CDEF levels for M5/6
#define OPT_SG                                    1 // Opt sg filter
#define OPT_SC_M1_ABOVE                           1 // optimize presets M7 and above for FPS in VBR mode, optimize ME_SA for SC for presets M1-M6
#define TUNE_SHIFT_PRESETS                        1 //Shift presets to match Master BD-rate
#define FIX_FI_R2R                                1 // Fix to allow FI to change at picture level
#define TUNE_4X4                                  1 // Use base and islice modulation for 4x4 levels
#define FIX_EXPOSE_4X4_DR                         1 // Expose depth removal feature that targets 4x4
#define TUNE_M7_M8                                1 // tune new m7-m8
#define FIX_DEPTH_REMOVAL_OFF                     1 // When depth removal is set to off, don't allow high_freq_present to turn it on
#define TUNE_NEW_M1                               1 // tune new m1
#define TUNE_M10                                  1 // tune m10
#define OPT_Q_ME                                  1 // Ctrl ME using q
#define OPT_PME_LVL                               1 // Improve the PME lvl
#define TUNE_NEW_M9                               1 // tune new m9
#define TUNE_M11                                  1 // tune m11
#define OPT_MRP_VBR_FPS                           1 // create separate control structure for mrp for VBR
#define TUNE_M12                                  1 // tune m12
#define TUNE_M8_III                               1 // tune m8
#define TUNE_TPL                                  1 // Tune  tpl levels
#define TUNE_M9_M10                               1 // tune m9-m10
#define OPT_LPD1_DET                              1 // Optimize LPD1 detector
#define TUNE_M6_M7                                1 // tune m6-m7
#define TUNE_M9_LPD0_LVLS                         1 // Remove QP banding from M9 LPD0 levels
#define TUNE_M6_CDEF                              1 // New M6 CDEF level
#define TUNE_M2_M4                                1 // Tune m2 m4 fast decode
#define TUNE_M1_II                                1 // tune m1 default
#define TUNE_M13                                  1 // Opt M13 bd-rate
#define FIX_DEPTH_R2R                             1 // Fix r2r from using uninit'd data in depth refinement
#define FIX_BYPASS_ED                             1 // Various bug fixes related to bypassing EncDec
#if FIX_BYPASS_ED
#define FIX_BYPASS_ED_COEFF                       1 // Fix how coeffs are saved when encdec is bypassed and fix chroma bugs for 128xN blocks
#define TUNE_BYPASS_ED                            1 // Bypass encdec for M5/6 in 8bit
#define FIX_BYPASS_ED_10BIT                       1 // Add missing neighbour array copies for 10bit when NSQ and bypass_endec are used
#define FIX_NZ_COEFF_SKIP                         1 // Set the count of non-zero coeffs to 0 when skip is selected in svt_aom_full_cost
#define FIX_RECON_COPIES                          1 // Fix the 8bit recon copies when 8bit MD is used for 10bit content
#endif
#define TUNE_M5_III                               1 // tune m5 fast decode
#define TUNE_M7_III                               1 // tune M7
#define TUNE_M8_IIII                              1 // tune m8
#define TUNE_M10_II                               1 // tune M10
#define TUNE_M11_II                               1 // tune M11
#define TUNE_M13_II                               1 // tune M13
#define TUNE_M12_II                               1 // tune m12
#define DIS_UNSUPPORTED_MODES                     1 // disable unsupported modes
#define OPT_TF_SUBPEL                             1 // use 8-bit subpel-search for 10-bit input(s)
#define FIX_LD_CBR_MODE                           1 // fix regression in branch in order to match v1.7.0 performance for ld-cbr mode - lpd1
#define FIX_LD_CBR_MODE_II                        1 // fix regression in branch in order to match v1.7.0 performance for ld-cbr mode - rdoq
#define FIX_LD_CBR_MODE_III                       1 // fix regression in branch in order to match v1.7.0 performance for ld-cbr mode - lpd0
#define OPT_M11_VBR_SPEED                         1 // move M11 in 1/2-pass VBR mode in order to fix the speed spacing between M10-M11-M12

#define TUNE_NSQ_HIGH_RES                         1 // tune nsq for 720p+ content

#define CLN_MISC_II                               1 // CI cleanup
#define FIX_MEM_ALLOC_ON_THE_FLY                  1 // Ensure memory is allocated for features whose level depends on resolution/qp (which can be updated on the fly)

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
