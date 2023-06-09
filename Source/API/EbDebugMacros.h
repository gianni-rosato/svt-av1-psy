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

#define FTR_ROI                 1 // implement region of interest feature
#define FIX_ISSUE_2064          0
#define FIX_ISSUE_2064_ALT      1

#define OPT_SUBPEL_ME_LD        1 // optimize subpel_me for low delay mode
#define OPT_LPD1_M13            1 // optimize pic_lpd1_lvl for M13
#define OPT_LD_MDS0             1 // Optimize MDS for better mode selection
#define OPT_LD_DLF              1 // Bypass dlf application for non-ref
#define OPT_LD_PD0              1 // Opimize PD0 interdepth decision
#define OPT_LD_EN_BASE_FRAME    1 // change ld_enchanced_base_frame frequency

#define OPT_LD_TF               1 // TF optimization for low delay mode

#define OPT_LD_M12              1
#define OPT_LD_M11              1
#define OPT_LD_M10              1

#define OPT_LPD0_8BIT_ONLY                         1 // Make MDS0 8bit only
#define OPT_SKIP_DISALLOWED_NSQ                    1 // Skip testing NSQ shapes when SQ block is disallowed (as cost of SQ is set to 0 and will be selected anyway)
#define TUNE_SIG_DERIV_HL                          1 // Use picture-level hierarchical levels in signal derivation functions at MDC and encdec
#define FIX_DEPTH_EXIT_CHECK                       1 // Move check for depth early exit to account for all quadrants
#define OPT_DEPTH_REFIN                            1 // Opt depth refinement
#define CLN_ME_DO_REF                              1 // Lossless Clean-Up for setting do ref flag in ME
#define OPT_ZZ_REF_PRUNE                           1 // ZZ sad based ref pruning
#define TUNE_11_SEC_TOP                            1 // limit 1,1 to safe second last layer
#define TUNE_MRP_M8_TO_M12                         1 // enable mrp-6 for M8-M12 + Use Ref scheme 0 for mrp level 6
#define OPT_PHME_PRUNE                             1 // PreHme based ref pruning
#define CLN_PHME_DATA                              1 //lossles fix for pre-hme to consider checking validity of list0 data
                                                     //TODO: remove performed_phme
#define GM_REFINFO                                 1 // Use References to determine wether or not to do Global
#define TUNE_GM_LEVELS                             1 // Tune GM Levels

#define NSQ_PACKAGE                                1 // Optimizations towards NSQ for M7 & M8
#if NSQ_PACKAGE
#define OPT_NSQ_M7_M8                              1 // Turn ON NSQ for M7 and M8
#define FIX_NSQ_LPD1                               1 // LPD1=f(NSQ)
#define OPT_NSQ_SUPPORT                            1 // Optimize NSQ infrastructure set0
#define OPT_NSQ_PSQ_MODE                           1 // Skip NSQ based on the pred-mode of the SQ
#define OPT_32x16_16x32_GEOM                       1 // Shut 16x8/8x16, and optimixed thegeom; from 169 to 105 blocks
#define FIX_NSQ_MD_SETTINGS_SHIFT                  1 // Optimize NSQ infrastructure set1
#define OPT_NON_NORMATIVE_TXS                      1 // Optimize the non-normative SQ txs
#define OPT_SKIP_NSQ_PATH                          1 // Optimize NSQ infrastructure set2
#define OPT_NSQ_PSQ_MODE_BIS                       1 // Optimize NSQ infrastructure set2

#define FIX_NSQ_SETTINGS                           1 // Optimize ParentSqCmplxCtrls

#define OPT_NSQ_SETTINGS                           1 // Optimize ParentSqCmplxCtrls
#if OPT_NSQ_SETTINGS
#define OPT_SHUT_PME_ACTION                        1 // Optimize ParentSqCmplxCtrls
#define OPT_SHUT_MRP_ACTION                        1 // Optimize ParentSqCmplxCtrls
#endif
#define TAG_SUB_PRED_NSQ                           1 // Opt the NSQ(s) that belong to the child depth(s); relative to PRED (the output of PD0)
#endif

#define CHROMA_PACKAGE 1 // All chroma optimizations and cleanups
#if CHROMA_PACKAGE
#define REMOVE_UV_CFL_TH                           1 // Remove uv_cfl_th from chroma ctrls because it is useless
#define CLN_IND_UV_SEARCH                          1 // Cleanup how independent chroma data is used
#define FIX_FI_UV_RATE                             1 // Fix rate used for independent chroma filter-intra candidates
#define FIX_UV_INTRA_MODE_RATE                     1 // Fix chroma fast rate when CFL is used. Enabled independent chroma modes for palette since chroma palette is not supported
#define CLN_CFL_SIGNALLING                         1 // Clean how CFL testing is signalled and add new fast chroma search path
#if CLN_CFL_SIGNALLING
#define TEST_UV_DC_ALWAYS 1 // Always test chroma DC prediction in independent chroma search
#endif
#endif
#define NEW_TPL_LVL_3                              1 // New TPL level 3 based on existing features
#define CLN_PSQ_CPLX_LVL                           1 // Remove unused NSQ levels

#define OPT_TF_REF_PICS                            1 // Modulate BASE and L1 tf ref-pics using the filt_INTRA-to-unfilterd_INTRA distortion range
#if OPT_TF_REF_PICS
#define INCREASE_REF_PICS                          1 // Limit the ref-pics increase for only high filt_INTRA-to-unfilterd_INTRA distortion range
#define ENHANCE_KEY_FRAME                          1
#endif
#define OPT_BIAS_PARENT                            1 // Apply a bias @ d2-decision towards a larger block(s)

#define GM_PP                                      1 //GM detection pre-processor for TF pictures
#define OPT_GM_CORNER                              1

#define TUNE_MRP_LEVELS                            1 //re-tune new MRP levels
#define FTR_ME_COST_VAR                            1 // Modify ME search area based on the 8x8 SAD variance of the search centre

#define OPT_CMP                                    1 //optimize inter compound

#define TUNE_NSQ                                   1 // Tune NSQ THs
#define TUNE_M12                                   1 // Tune m12 for vmaf neg
#define TUNE_M11                                   1 // Tune m11 for vmaf neg
#define TUNE_M1                                    1 // Tune m1 for vmaf neg
#define TUNE_M2                                    1 // Tune m2 for vmaf neg
#define TUNE_M3                                    1 // Tune m3 for vmaf neg
#define TUNE_M0                                    1 // Tune m0 for vmaf neg
#define TUNE_M4                                    1 // Tune m4 for vmaf neg
#define TUNE_M5                                    1 // Tune m5 for vmaf neg
#define TUNE_M6                                    1 // Tune m6 for vmaf neg
#define TUNE_M7                                    1 // Tune m7 for vmaf neg
#define TUNE_M8                                    1 // Tune m8 for vmaf neg

#define OPT_NIC_PER_BSIZE                          1 // Add NIC scaling settings for block size

#define DEPTH_PACKAGE                              1
#if DEPTH_PACKAGE
#define OPT_DEPTH_MODULATION_WIDE                  1 // Shut depth-level modulation when the band is not wide enough
#define OPT_DEPTH_LVL_TUNING                       1 // Use a more agressive depth-level

#define OPT_SUB_DEPTH_ALL_EXIT                     1 // Cond0: use the nsq-to-sq cost deviation to skip sub-depth(s)
                                                     // Cond1: use the 4 quad(s) src-to-recon cost deviation to skip sub-depth(s)
#endif
#define CLN_TXS                                    1 // Lossless TXS cleanup
#define CLN_MVP_ARR                                1 // clean up for the MVP array
#define OPT_SPEL                                   1 // optimize subpel
#define FIX_TXT                                    1 // Fix TXT restriction that prevented some allowable TX types from being searched
#define OPT_NSQ_HIGH_COST_COMP                     1 // Skip NSQ shapes when the cost of the SQ is concentrated in either distortion or rate
#define CLN_FUNC_DECL                              1 // Move some hanging functions declarations to appropriate header files

#define RATE_FIXES 1
#if RATE_FIXES
#define FIX_SKIP_COST                              1 // Fix how skip cost is computed
#define FIX_II_RATE                                1 // Fix inter-intra rate in inter fast cost
#define FIX_COEF_RATE_UPDATE                       1 // Fix coeff rate updates at EncDec
#define FIX_EOB_RATE_UPDATE                        1 // Fix eob rate updates at EncDec
#define FIX_TXS_RATE_UPDATE                        1 // Fix TXS rate updates at EncDec
#define FIX_IFS_RATE_UPDATE                        1 // Fix IFS rate updates at EncDec
#endif

#define OPT_SG                                     1 // optimize SG filter
#define FIX_SKIP_MODE_RATE                         1 //Fix Skip mode rate estimation
#define FIX_SCT                                    1 //Fix frame header allow-sc-tool flag to be off when palette level is off
#define FIX_ME_HME_PRUNE_CTRLS                     1 // fix missing break in me_hme_prune_ctrls, re-adjust levels
#define OPT_INTRA_LEVELS                           1 // Optimize intra levels used for M0-M3

#define FIX_NSQ_COEFF_ACTION                      1 // fix onion ring for nsq coeff action
#define FIX_IFS_PATH                              1 // Allow IFS for <8x8 blocks and fix IFS rate calc

#define OPT_SQ_WEIGHT                             1 // Optimize sq-weight
#define FIX_ISSUE_2070                            1 // Fix issue #2070

#define TUNE_M0_2                                 1 // Tune m0 for vmaf neg
#define TUNE_M1_2                                 1 // Tune m1 for vmaf neg
#define TUNE_M5_2                                 1 // Tune m5 for vmaf neg
#define TUNE_M13_2                                1 // Tune m13 for vmaf neg
#define TUNE_M12_2                                1 // Tune m12 for vmaf neg
#define TUNE_M11_2                                1 // Tune m11 for vmaf neg
#define TUNE_M2_2                                 1 // Tune m1 for vmaf neg
#define TUNE_M3_2                                 1 // Tune m5 for vmaf neg
#define TUNE_M7_2                                 1 // Tune m7 for vmaf neg
#define TUNE_M4_2                                 1 // Tune M4 for vmaf neg
#define TUNE_M6_2                                 1 // Tune M6 for Vmaf neg
#define TUNE_M10_2                                1 // Tune m10 for vmaf neg
#define TUNE_M9_2                                 1 // Tune m9 for vmaf neg
#define TUNE_M8_2                                 1 // Tune m8 for vmaf neg

#define CLN_IFS_PATH                              1 // Cleanup the IFS path and unused signals
#define CLN_OBMC_WARP                             1 // Cleanup OBMC/Warp

#define FIX_INTER_CMP                             1 // Fix inter compound to not send side information when compound is off.
#define FIX_FILTER_INTRA                          1 // Fix filter intra to not send side information when filter intra is off.
#define CLN_REMOVE_NEIGH_ARRAYS                   1 // Remove unused neighbor arrays
#define CLN_REMOVE_NEIGH_ARRAYS_2                 1 // Remove unused neighbor arrays
#define FIX_SKIP_NEIGH_ARRAY                      1 // Fix how skip_coeff context is derived; remove neighbour array
#define CLN_MISC_CLEANUPS                         1 // Miscellaneous cleanups*
                                                    // *TODO: Involves removing a file: EbSyntaxElements.h
                                                    // *TODO: Renames sb_type --> bsize without a macro. Before renaming, in cdef_seg_search()
                                                    //        rename "int32_t bsize[3];" to "int32_t plane_bsize[3];" to avoid conflicts

#define TUNE_4K_8K                                1 // Tune memory sensitive settings for 4K+ resolutions . Note we can remove the debug macro: REMOVE_LP1_LPN_DIFF
#define CLN_TF                                    1 // Cleanup tf

#define OPT_OBMC_UV_ONLY                          1 // Allow OBMC prediction to be done for chroma only
#define OPT_WEDGE_UV_ONLY                         1 // Allow wedge compound prediction to be done for chroma only
#define OPT_WM_UV_ONLY                            1 // Allow WM prediction to be done for chroma only

#define TUNE_FAST_DEC_5L                          1 // Make fast decode use fixed 5L
#define TUNE_REST_PROCESS                         1 // Tune number of restoration processes for nlp
#define MEM_SG                                    1 // Memory optimization for self-guided

#define TUNE_M1_FD                                1 // Tune m1 fast decode
#define TUNE_M2_FD                                1 // Tune m2 fast decode
#define TUNE_M3_FD                                1 // Tune m3 fast decode
#define TUNE_M4_FD                                1 // Tune m4 fast decode
#define TUNE_M10_FD                               1 // Tune m10 fast decode
#define TUNE_M5_FD                                1 // Tune m5 fast decode
#define TUNE_M6_FD                                1 // Tune m6 fast decode
#define TUNE_M8_FD                                1 // Tune m8 fast decode
#define TUNE_M7_FD                                1 // Tune m7 fast decode

#define OPT_CHILD_PCS                             1 //memory optimization to child pcs final block data
#if OPT_CHILD_PCS
#define SHUT_TEMP_BUFF                            1
#define SHUT_TOTAL_RATE                           1
#define SHUT_PART                                 1
#define SHUT_BEST_D1_BLK                          1
#define SHUT_SEG_MASK                             1
#define CLN_MBMI                                  1
#define CLN_MBMI_2                                1 //skip_mode
#define CLN_MBMI_3                                1 //prediction mode flag
#define CLN_MBMI_4                                1 //pred mode
#define CLN_MBMI_5                                1 //interp filters
#define CLN_MBMI_6                                1 //compound_idx + compound group idx
#define CLN_MBMI_7                                1 //tx depth
#define CLN_MBMI_8                                1 //use intra bc
#define CLN_MBMI_9                                1 //intra chroma mode
#define CLN_MBMI_10                               1 //split flag
#define CLN_MBMI_11                               1 //mvs
#define CLN_MBMI_12                               1 //inter mode ctx
#define CLN_MBMI_13                               1 //ref frame

#endif

#define TUNE_MR                                   1 // Tune MR
#define OPT_STARTUP_MG                            1 // Optimize performance when a shorter MG is used at the start of a GOP

#define FIX_ISSUE_2078                            1 // Fix issue #2078
#define FIX_GM_PP                                 1 // fix non-lp R2R

#define CLN_SKIP_HV4                              1 // Clean up unused skip hv4 signal
#define ENABLE_LD_DLF_RA                          1 // Enable LD DLF changes in RA
#define CLN_FAST_COST_FUNCS                       1 // Cleanup fast cost functions
#define CLN_UNUSED_DEFNS                          1 // Remove unused definitions
#define FIX_GM_CI                                 1 // Clean Up GM code

#define OPT_LD_LATENCY2         0 // Latency optimization for low delay
#define FIX_DEADLOCK            1 // Deadlock Fix
#define OPT_LD_M12_SC           1
#define OPT_LD_M8_SC            1
#define OPT_LD_M9_SC            1
#define OPT_LD_SC_PRESETS       1
#define OPT_LD_SC_RC            1
#define OPT_LD_SC_ME            1

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
#define DEBUG_SEGMENT_QP        0
#define DEBUG_ROI               0
#ifdef __cplusplus
}
#endif // __cplusplus

// clang-format on

#endif // EbDebugMacros_h
