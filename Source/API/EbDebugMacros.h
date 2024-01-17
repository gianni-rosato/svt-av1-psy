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

// svt-09 macros
#define OPT_NSQ_MV                                1 // Improve the refinement of the NSQ MV(s)
#define OPT_TF_PATH                               1 // Clean TF path
#define OPT_Q_PRUNE_TH_WEIGHT                     1 // Modulate mds-pruning th(s) using q and distortion
#define OPT_OBMC                                  1 // To do
#define OPT_MERGE_INTER_CANDS                     1 // Perform the inter-candidates merge @ block-basis using q and best PME/ME pred-error, instead of @ input-basis using coeff-level
#define OPT_TF_SP_EXIT                            1 // Apply TF subpel_early_exit feature to all bsizes with same TH for all
#define OPT_COEFF_LVL_NORM                        1 // Normalize the coeff_lvl feature
#define FIX_NSQ_CTRL                              1 // Break the nsq_ctrls into nsq_geom_ctrls and nsq_search_ctrls and use qp-banding for only the nsq_search_level derivation
#define CLN_SMALL_SIGS                            1 // Remove signals and levels having an insignificant impact on the behaviour
#define OPT_ME_SP_8TH_PEL                         1 // Optimize ME subpel
#define OPT_NSQ_HIGH_FREQ                         1 // Use more conservative NSQ settings in the presence of a high energy area
#define OPT_1P_VBR                                1 // Optimized one-pass VBR
#define CLN_MVP_DIST_CALC                         1 // Move MVP distortion calc to one place to avoid recomputing
#define CLN_USE_BEST_PME_DIST                     1 // Use best PME dist in cand. reduction tool, instead of just best of ref_idx 0 cands
#define CLN_ADD_FIXED_PRED_SIG                    1 // Add signal to track if PDx has a fixed prediction structure
#define OPT_HME_L0_Q                              1 // Use Q @ the derivation of HME-L0 search-area (w,h)
#define CLN_NSQ_COPIES                            1 // Cleanup copying of neighbour arrays for NSQ shapes to avoid unnecessary copying
#define FIX_NSQ_SETTINGS                          1 // Fix when settings are reset in MD loop (settings can be modified by NSQ features)

#define OPT_MPASS_VBR_PACKAGE                     1 // Optimized Multi-pass VBR
#if OPT_MPASS_VBR_PACKAGE
#define OPT_MPASS_VBR1                            1 // Remove the dependancy of IPP and middle pass
#define OPT_MPASS_VBR2                            1 // Lossless optimization for multiPASS VBR
#define OPT_MPASS_VBR3                            1 // Lossless optimization for multiPASS VBR, remove IPP pass
#define OPT_MPASS_VBR4                            1 // Lossless optimization for multiPASS VBR, remove ENC_FIRST_PASS
#define OPT_MPASS_VBR5                            1 // Lossless optimization for multiPASS VBR, remove FIRST PASSME
#define OPT_MPASS_VBR6                            1 // Lossless optimization for multiPASS VBR, Shift passes
#define OPT_MPASS_VBR7                            1 // Lossless optimization for multiPASS VBR, replace ENC_MIDDLE_PASS by ENC_FIRST_PASS
#define OPT_MPASS_VBR8                            1 // Lossless optimization for multiPASS VBR, refactor mid_pass
#define FIX_REMOVE_PASS3                          1 // Change bounds of --pass to exclude pass 3, which has been removed
#endif

#define OPT_NSQ_INCOMP                            1 // Optimize how NSQ is used for incomplete SBs
                                                    // Enable NSQ shapes to be used in LPD0 (all sizes) and LPD1 (8x8 and higher) for incomplete blocks when the SQ is not allowed.
                                                    // When nsq geom is enabled, incomplete blocks will always use NSQ shapes (when SQ is disallowed), even when md_disallow_nsq is true.
                                                    // TODO: change md_disallow_nsq to disallow search, or add enabled signal to nsq_search_ctrls.

#define TUNE_M11_M13                              1 // tune m11-m13 lp8
#define TUNE_M9_M10                               1 // tune m9-m10 lp8
#define TUNE_MR                                   1 // Temporary mr candidate
#define TUNE_M0_M1                                1 // tune m0-m1 lp8
#define TUNE_M7_M8                                1 // tune m7-m8 lp8
#define TUNE_M2                                   1 // tune m2 lp8
#define TUNE_M5_M6                                1 // tune m5-m6 lp8 fd
#define TUNE_M3_M4                                1 // tune m3-m4 lp8
#define TUNE_M5_M6_2                              1 // speed alignment m5-m6 lp8 fd 2

#define CLN_MD_LOOP                               1 // Cleanup how NSQ is processed in MD
#define OPT_NSQ_MEM                               1 // Optimize memory of NSQ blocks in MD search
                                                    // Merge MdBlkStruct with BlkStruct - no macros
#define CLN_BLK_STRUCT                            1 // Cleanup fields of BlkStruct
#define CLN_SB_ME_MV                              1 // Cleanup how MVs are store in MD ctx
#define CLN_MDC_ARRAY                             1 // Allocate MD scan array for the correct number of blocks
#define CLN_QUAD_REC                              1 // Move rec_dist_per_quadrant from BlkStruct to MD ctx. Fix redundant blocks using incorrect data for SQ blocks.
#define CLN_BLK_STRUCT_2                          1 // Remove prediction_unit_array from BlkStruct
#define CLN_BLK_STRUCT_3                          1 // Reorder fields of BlkStruct
#define CLN_MOVE_CFL_BUFF                         1 // Allocate temp cfl buffer on heap
#define CLN_MOVE_PAL_BUFF                         1 // Allocate PALETTE_BUFFER on heap
#define FIX_REDUND_PAL                            1 // Add missing copies for palette for redundant blocks
#define CLN_SEG_MASK                              1 // Remove seg_mask* from InterInterCompoundData and remove the now identical EcInterInterCompoundData
#define CLN_INTER_MODE_CTX                        1 // Move inter_mode_ctx array to MD ctx and save only the winning result in blk_ptr
#define CLN_BLK_STRUCT_4                          1 // Move min_mz_h/v to ctx
#define OPT_COEFF_LVL_NOISE                       1 // Use the input noise-level @ the derivation of the coeff-level low/high bands
#define OPT_BLOCK_SETTINGS                        1 // Adopt same nic scaling per bsize and remove sig_deriv_block
#define CLN_UNUSED_SETTING                        1 // Remove unused setting in update_md_settings
#define CLN_TX_DATA                               1 // Cleanup TX info in BlkStruct
#define CLN_EC_PAL_STRUCT                         1 // Remove EcPaletteInfo becuase it's the same as PaletteInfo
#define CLN_QUANT_ONE_BYTE                        1 // Update MD cand buff to use new Eob and QuantDc structs; save quant DC is single byte, since that's how it's used

#define TUNE_LPD1                                 1 // tune m10-m11 lpd1 level
#define CLN_QUANT_FUNC                            1 // Remove useless arguments from quantization function
#define CLN_EC_BLK_STRUCT                         1 // Cleanup fields in EcBlkStruct
#define CLN_MD_DISALLOW_NSQ                       1 // Cleanup md_disallow_nsq signal to refer only to the nsq search
#define OPT_CHECK_SKIP_NSQ_ONCE                   1 // Only call get_skip_processing_nsq_block and faster_md_settings_nsq for the first block in an NSQ shape
#define FIX_NSQ_HIGH_ENERGY                       1 // Move resetting of signals in MD to not overwrite NSQ settings when high_energy_weight is enabled
#define CLN_TXTYPE_WIN                            1 // Force TxType type to be 1 byte on windows to save memory in EcBlkStruct

#define OPT_SB128                                 1 // sb128 = f(resolution)
#define CLN_ENCDEC_FUNCS                          1 // Merge 8bit/16bit encdec functions
#define FIX_RECON_PADDING                         1 // Fix recon buffer padding (top padding is set to bottom and vice-versa)
#define TUNE_M12_II                               1 // tune m12 lp8
#define CLN_SKIP_PD0_SIG                          1 // Remove skip_pd0 signal because not used for skipping PD0

#define TUNE_M5_2                                 1 // tune m5 lp1
#define TUNE_M6                                   1 // tune m6 lp1
#define OPT_DEPTH_FEAT_128                        1 // enable depth removal and other depth shortcut relying on ME info when using 128x128 SB
#define OPT_USE_FAST_TX_PATH_B64                  1 // Use optimized DCT_DCT only path for block sizes <= 64 (instead of checking SB size)
#define OPT_DR_COEFF_LVL                          1 // Set depth removal level using coeff_level

#define DIS_DLF_SG_QP                             1 // remove fd dlf / sg qp bands
#define DIS_COMP_QP                               1 // remove get_th_qp bands
#define TUNE_TPL_LVL                              1 // Use tpl-group cplx @ the derivation of the tpl-params-lvl

#define TUNE_M3_M4_2                              1 // tune def m3-m4 lp1
#define FIX_PSQ_TXS_UPDATE                        1 // Fix updating psq_txs ctrls from high-frequency NSQ feature

#define TUNE_M7_M8_2                              1 // tune def m7-m8 lp1
#define TUNE_M9_M10_2                             1 // tune def m9-m10 lp1
#define TUNE_M0                                   1 // tune def M0 lp1
#define TUNE_M1                                   1 // tune def M1 lp1

#define TUNE_M5_M6_3                              1 // tune fd m5-m6 lp1
#define FIX_REDUND                                1 // Fix updates for redundant blocks
#define OPT_PALETTE_MEM                           1 // Only allocate palette memory for blocks where palette is allowed
#define TUNE_LD                                   1 // Tune LD NSC settings
#define CLN_REMOVE_UNUSED_ARGS                    1 // Removed unused args from RDOQ func
#define CLN_M13_CHECKS                            1 // Remove checks on <=M13 because they're useless

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
