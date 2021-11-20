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

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus
#define FTR_RC_CAP             1 // Add support for Capped CRF
                                 // Add the support to re-encode when max bitrate is specified for CRF
                                 // Add the rate stat queue for the storing rate data
#define FTR_2PASS_CBR          1 // to support 2pass CBR
#define FTR_1PASS_CBR          1 // to support 1pass CBR
#define FTR_1PASS_CBR_FIX      1 // 1pass CBR fix for pcs_total_rate
#define FTR_1PASS_CBR_RT       1 // to support 1pass CBR real time mode
#define FTR_1PASS_CBR_RT_MT    1 // to improve quality for 1pass CBR real time multi thread
#define FTR_1PASS_CBR_RT_MT_TUNE 1 // improve quality for 1pass CBR real time
#define FTR_1PAS_VBR                1 // 1 PASS VBR with ability to use LAD stats
                                      // Calculate stat total in 1 pass VBR based on the look ahead
#define FTR_LAD_INPUT               1 // Add look ahead as an input to the encoder
#define FTR_2PASS_1PASS_UNIFICATION 1 // Unifty the VBR algorithms between 2 pass and 1 pass
                                      // Create a ME functions for the IPPP processing of 2PASS and 1 PASS VBR/CBR

#define FIX_LOW_DELAY               1 // Fix the 6L, 5L and 4L with LOW_DELAY_P
// Multi pass encode
#define FTR_NEW_MULTI_PASS          1 // New multipass VBR
#define FTR_MULTI_PASS_API          1 // the API needed to run multi pass (3 or 2) VBR

#define TUNE_MULTI_PASS             1 // Tune the multipass VBR. Add the option to remove IPPP
#define FTR_OPT_MPASS               1 // Optimize the middle pass for better speed

#define RFCTR_RC_P1                 1 // Rate control code refactoring Part 1
#define RFCTR_RC_P2                 1 // Rate control code refactoring Part 2

#define TUNE_RC                     1 // Tune RC setting for each preset
#if TUNE_MULTI_PASS
typedef enum MultiPassModes {
    SINGLE_PASS, //single pass mode
    TWO_PASS_IPP_FINAL, // two pass: IPP + final
    TWO_PASS_SAMEPRED_FINAL, // two pass: Same Pred + final
    THREE_PASS_IPP_SAMEPRED_FINAL, // three pass: IPP + Same Pred + final
} MultiPassModes;
#endif
#if FTR_OPT_MPASS
#define FTR_OPT_MPASS_CDEF 0
#define FTR_OPT_MPASS_RDOQ_OFF 0
#define FTR_OPT_MPASS_NEAR0 0
#define FTR_OPT_MPASS_MRP1REF 0
#define FTR_OPT_MPASS_WARP_OFF 0
#define FTR_OPT_MPASS_DLF_OFF 0
#define FTR_OP_TEST 0 // for debugging to run the optimization of middle pass in normal CRF mode
#define FTR_OPT_MPASS_DOWN_SAMPLE 1 // Down sample the input picture by 2 in each dimention to be used in middle pass

#define FTR_OPT_MPASS_BYPASS_FRAMES 1
#endif
#define FTR_OPT_IPP_DOWN_SAMPLE 1 // Down sample the input picture by 2 in each dimention to be used in IPP pass
#define DS_SC_FACT 23
#define SVT_05          1
#define PRIVATE_MACROS   1

#define FIX_PA_REF_RELEASE_HANG               1 // Fix a hang that occured for incomplete MGs when MRP is ON
#define FIX_QUANT_COEFF_BUFF                  1 // Fix how the quant coeff buffer is released, and cleanup the init

#if SVT_05 //SVT05_TOBE_PORTED
#define CLN_FTR_EARLY_DEPTH_REMOVAL        1 // Clean-up cppcheck issues introduced with the FTR_EARLY_DEPTH_REMOVAL macro
// ============= START SVT_05 =============
// depth-skip-ctrls
#define OPT_REFACTOR_IN_DEPTH_CTRLS            1 // Refactor in_depth_skip_ctrls
#if OPT_REFACTOR_IN_DEPTH_CTRLS
#define OPT_IN_DEPTH_SKIP_M9                   1
#endif
#define OPT_REFACTOR_DEPTH_REFINEMENT_CTRLS    1 // Refactor depth_refinement_ctrls
#if OPT_REFACTOR_DEPTH_REFINEMENT_CTRLS
#define OPT_DEPTH_REFINEMENT_M3_M4             1
#define TUNE_PRED_DEPTH_M9                     1 // turn on pred depth for m9
#endif
// depth-skip-ctrls
#define TUNE_INTRA_LEVELS                      1 // Create new intra levels and assign them to presets
#define FTR_MULTI_STAGE_CDEF                   1 // multi_stages CDEF
#if FTR_MULTI_STAGE_CDEF
#define MS_CDEF_OPT1                           1 // Optimised skip cdef condition
#define MS_CDEF_OPT2                           1 // Add the ability to reduce the number of filter strengths for chroma components
#define MS_CDEF_OPT3                           1 // Add the ability to use the best filter strength from previously coded frames
#endif
#define TUNE_SG_FILTER_LVLS                    1 // Fix OFF level for sg_filter_mode; re-order existing levels to follow conventions
#define FTR_NEW_WN_LVLS                        1 // Add new wiener filter levels
#define TUNE_DLF_LVLS                          1 // New DLF level for M5-M8
#define SIMD_APPROX_EXPF                       1 // Optimize exp @ TF weight(s) generation
#define TUNE_NEW_M9_LEVEL                      1 // Adopt M8 levels in M9 for certain features to gain bd-rate in M9
#define TUNE_M9_23_BDR                         1 // push features to M10 to achieve 23% bd-rate loss compared to M8
#define TUNE_RDOQ_LEVEL_M9                     1 // replace level 2/3 with modified level 1
#define TUNE_TXT_LEVEL_M9                      1 // tune level 6 to be used in M9 in base layer
#define CLN_REDUCE_CDEF                        1 // Reduce fixed cost of cdef list computation and skipping decision
#define OPT_INTER_PRED                         1 // Skip luma prediction @ mds3 when IFS is OFF, and skip luma prediction if regular @ md3 when IFS is ON
#define FIX_REMOVE_PD1                         1 // Remove PD_PASS_1
#define OPT_TXS_SEARCH                         1 // Improve TXS ctrls, and TXS performance; up to 1 tx_depth for INTRA, and decrement the txt group per an offset if tx_depth 1 or higher
#define TUNE_NEW_TXT_LVLS                      1 // Improve TXT levels
#define TUNE_CDEF_TF_LEVELS                    1 // tune both cdef_level and tf_level for M9
#define FTR_LIMIT_GM_REFINEMENT                1 // Add ability to limit the GM params refinement; use in M9
#define TUNE_IN_DEPTH_LIST0_M9_LEVEL           1 // tune in_depth_block_skip_level and list0_only_base for M9
#define FTR_NEW_WM_LVL                         1 // Add new intermediate level for WM
#define CLN_FP_SEARCH_KERNEL                   1 // Do not use md_sq_me_ctrls under md_full_pel_search()
#define FTR_IMPROVE_PME                        1 // Do not perform PME search for blocks that have a valid ME_MV unless the ME_MV has a different direction than all MVP(s) and the ME_MV mag is higher than MV_TH (not around(0,0))
#define CLN_NIC_PRUNING                        1 // clean up nic pruning levels
#define FTR_RDOQ_ONLY_DCT_DCT                  1 // RDOQ for only DCT-DCT (md_context->rdoq_ctrls.dct_dct_only)
#define TUNE_M10_INITIAL_PRESET                1 // create new M10 preset
#define TUNE_TXS_IFS_MFMV_DEPTH_M9             1 // tune M9/M10 levels for txs, ifs, mfmv, block-based-depth-refinement, tf
#define TUNE_HME_SR                            1 // tune HME and pre-HME search areas
#define FTR_ENHANCE_FAST_COEFF_RATE            1 // Use fast-coeff-rate estimation @ mds1 of PD1, and tune the PD0 fast-coeff-rate estimation;
#define CLN_ADD_LIST0_ONLY_CTRL                1 // add list0_only_base controls
#define TUNE_NEW_M11                           1 // create new M11 level

#define OPT_REDUCE_TX                          1 // extend SB-based PF detector to PD1, and use MDS1 coeff info to cut TXT path in MDS3
#define FTR_QP_BASED_DEPTH_REMOVAL             1 // Modulate depth_removal level @ BASE based on the qp_offset band
#define TUNE_M7_M10_PRESETS                    1 // tune M7-M10 presets
#define FTR_CDEF_SUBSAMPLING                   1 // Allow subsampling (i.e. skipping rows) in the CDEF search
#if FTR_CDEF_SUBSAMPLING
#define FTR_SKIP_LINES_CDEF_SEARCH             1
#define FTR_SKIP_LINES_CDEF_DIST               1
#endif

#define FTR_TPL_TX_SUBSAMPLE                  1 // Subsample TX in tpl
#define FTR_SUBSAMPLE_RESIDUAL                1 // Add the ability the subsample the residual @ mds3 of PD0
#if FTR_SUBSAMPLE_RESIDUAL
#define TUNE_DEPTH_REMOVAL_1080P              1
#define FTR_TX_SUB_PD0_INTRA                  1
#endif
#define PD1_SUBSAMPLE_RESIDUAL                1 // Add the ability the subsample the residual @ mds1 of PD1
#define TUNE_M8_M9_FEB24                      1 // tune m8 and m9
#define TUNE_M9_TPL_ME_HME                    1 // tune tpl, prehme, me, hme for m9
#define TUNE_CTR_QUARTER_PEL                  1 // Added the ability to signal when to stop subpel
#define TUNE_M0_M7_MEGA_FEB                   1 // tuned features for m0-m7
#define TUNE_MATCH_04_M8                      1 // synchronize svt-04 m8 with svt-05
#define TUNE_SHIFT_PRESETS_DOWN                1 // shift higher presets down all the way to M2

#define TUNE_M3_M6_MEM_OPT                    1 // tune M3-M6



#define CLN_MOVE_SKIP_MODE_CHECK              1 // Move check on skip_mode to MD from EncDec; semi-lossless
#define CLN_ENC_DEC                           1 // Clean up the EncDec function



// cleanup


#define OPT8_MDC                              1
// Reduce the number of context for partition.
#define OPT9_RATE_ESTIMATION                  1
// tpl optimisation for M9
#define OPT5_TPL_REDUCE_PIC                   1 // Please do not test this macro - feature under contruction.

#define TUNE_BLOCK_SIZE                       1 // optimize ifs for M9/M10
                                                // optimize spatial sse using block size modulation
                                                // add block size and resolution checks to most aggresive levels
#define OPT11_SUBPEL                          1 // Skip next refinement when cost is not improved after checking left-top-right-bottom

#define REFCTR_SEP_ENCDEC                     1 // Separate EncDec into two parts: true Encdec, and filtering/updates
#define REFCTR_ADD_QUANT_COEFF_BUFF_MD        1 // Move quant_coeff buffer to be per candidate to forward coeffs for best candidate
#define FTR_BYPASS_ENCDEC                     1 // if signal enabled, FWD quant_coeffs and recon from MD to EncDec and overwrite the EncDec buffers with the MD results; 8bit only
                                                // Disables features that cause recon/quant coeffs to mismatch between MD/EncDec; shut_skip_ctx_dc_sign_update
                                                // fundamentally different between MD/EncDec - must address this separately
                                                // Tune M9 settings to acheive best trade-offs when EncDec is skipped

#define TUNE_M9_M10_MAR                       1 // tune M9 and M10
#define OPT_REDUCE_CDEF_MEMSET                1 // Reduce memory initialization in cdef search
#define OPT14_TPL                             1 // tune M9 tpl
#define TUNE_M10_DEPTH_ME                     1 // create new depth removal level, update me area for m10 and push features out for better bd-rate
#define OPT_RECON_COPY                        1 // Optimize copying the recon when EncDec is bypassed
                                                // Skip copying recon for CFL when CFL is disabled
#define CLN_DEPTH_REFINEMENT                  1 // cleanup perform_pred_depth_refinement()
#define OPT_TF                                1 // TF optimisation
#define TUNE_M9_MARCH                         1 // tuning M9
#define OPT_INLINE_FUNCS                      1 // inline oft-used functions: av1_ref_frame_type() and av1_set_ref_frame()
                                                // inline get_blk_geom_mds() and removed unnecessary initializations for blk geom data
#define FIX_UPDATE_GM_SIGS                    1 // Set the GM controls properly for ME_PCS, used at ME
#define TUNE_MEGA_M9_M4                       1 // tuning m9-m4
/*
NOTE : PART OF LIGHT_PD0  code was committed to svt-04-final-rebased under OPT_INIT_XD
*/
#define LIGHT_PD0                             1 // Use an efficient MD path for PD0 that assumes most features are off
                                                // udpate intra context updates only if shut_fast_rate is OFF (inlcudes light PD0 path and regular path)
                                                // skip partition context update in light PD0 path - not lossless
                                                //LIMIT PD0 to 2 ME candidates. lossy
                                                //to add prefetch. lossless
                                                //assumes 2 ME candidates as the nunber of cands could go up to 10 and max nic=6. lossy due to sorting.
                                                // don't update lambda for PD0 - macro can be removed
                                                // remove unnecessary inits in light PD0
                                                // Move init of blk_ptr->avx1d to one place (affects PD0 + PD1)
                                                // Move parent_depth_idx_mds, d1_depth_offset, and ns_depth_offset to blk_geom
                                                // skip d1/d2 decision when unnecessary.  Affects PD0 and PD1; semi-lossless when light-PD0 is used.


/*
NOTE : PART OF LIGHT_PD0_2  code was committed to svt-04-final-rebased under OPT_INIT_XD_2 & OPT_INIT
*/
#define LIGHT_PD0_2                           1 // Use optimized parent depth comparison for skipping MD blocks. Merge with LIGHT_PD0
                                                // Remove unnecessary code/initializations for light-PD0 candidate generation. Merge with LIGHT_PD0
                                                // Reduce number of data initializations (PD0 + PD1)
                                                // Clean up the light-PD0 TX path. Merge with LIGHT_PD0
                                                //merge with LIGHT_PD0_2
                                                //merge with LIGHT_PD0_2
                                                //light rate estimatinon for PD0, one Coeff optimization for for both PD0/PD1

                                                //leight weight Quant
                                                //remove unessary zero out of transformed coeff for TX_64xn

                                                // Move init of blk_ptr->avx1d to one place (affects PD0 + PD1)
                                                // remove unused functions


#define TUNE_NEW_M10_M11                      1 // New M10 and M11

#define FTR_OPT_PF_PD1                        1 //Use PF_N4 under certain threshold
#define FIX_PRED_DEPTH_REFINEMNT_4X4          1 // limit pred depth refinement to go to 4x4 if 4x4 is disabled.
#define CLN_GEOM                              1 // Lossless. add geom types and enable geom85 when SB64 + NSQ OFF + 4x4 OFF
#define CLN_SSE                               1 // Cleanup picture_full_distortion32_bits(); each function call now does only one component
#define OPT_CHROMA_PATH                       1 // Merge full_loop_r() and cu_full_distortion_fast_txb_mode_r(); remove unnecessary distortion calcs. Includes fix for 128x128 coeff-rate est.
#define OPT_TX_PATH                           1 // Use PF size in distortion calc. To be tested.
                                                // Cleanup the PD1 tx path

#define FTR_PREHME_SUB                        1 // skips every other search region line in pre-Hme

#define CHROMA_CLEANUP                        1 // Clean up Chroma in MD
#define TUNE_M10_M7                           1 // tuning M10-M7; slowing down M8 and M9; speeding up M7 and M10


#define TUNE_FIRSTPASS_HME2                     1 // firstpass optimizations: turn OFF hme LVL2 // OMK DIFF


#define TUNE_TF                                 1 // SAD subsampling in TF

#define OPT_ME                                1 //opt construct me results

#define TPL_CLEANUP                           1 // Clean up tpl in MD
#define FTR_LOW_AC_COST_EST                   1 // low precision fast cost calculation
#define FTR_LOW_AC_SUBPEL                     1 // lower precision subpel
#define CLN_INTER_PRED                        1 // Cleanup av1_inter_prediction()



#define TUNE_HME_SUB                          1 // skips every other search region line in HmeLvl0

#define CLN_DLF_SIGNALS                       1 // cleanup DLF signals


#define FTR_DEPTH_REMOVAL_QP                  1 // new depth removal levels
#define TUNE_M11                              1 // tune m11

#define FTR_USE_PSAD                          1 // added the pSad path in all md searches - pSad is enabled in pme for M6-M9

#define TUNE_M7_M10_MT                        1 // tuned M7-M10 using multithreaded and lp1 together

#define OPT_PD0_PATH                          1 // A light-PD0 version of copying the neighbour arrays, since only recon is needed and only when intra is not skipped
                                                // Create a light-pd0 signal derivation function (SB-level)




#define OPT_TFILTER                           1 //fast temporal  filter
                                                //avoid 8x8 ME results storage for TF. Lossless all presets (to fix Unit test)
                                                //fast TF subpel. avoid 2D quarter pel positions and use of 2 tap filter

#define SS_FIX_TF_BUG                         1 //set a TF signl that might cause R2R. to link to OPT_TF if possible

#define FTR_SKIP_MDS1                         1 // Skip MDS1 if only one candidate is output from MDS0; apply MDS3 shortcuts if distortion/QP < TH
#define FIX_OPT_RECON_COPY                    1 // Bug fix for OPT_RECON_COPY - to be merged with original macro when merged to svt-05-temp

#define OPT_PD0_TXT                           1 // bypass coeff rate estimation when the number of coeff is small

#define SS_OPT_MOVE_SC_DETECTION              1 // Perform SC detection only for I_SLICE, since other frame types follow I_SLICE SC class
#define OPT_UV_USE_Y_COEFF                    1 // Use MDS3 luma coeff info for chroma search, instead of MDS1 coeff info

#define TUNE_M5_M6                            1 // tuning m5-m6

#define FTR_PD_EARLY_EXIT                     1 // Early exit from mds0 in pd0 and bias inter depth decision towards large blocks
#define OPT_SUBPEL                            1 // Early exit from md_subpelsearch

#define TUNE_M10_FASTER                       1 // tuning m10 to make it faster
#define SS_OPT_CDEF                           1 // lossless optimizations for CDEF path
#define OPT_CDEF_LEVELS                       1 // Create new CDEF levels
#define FTR_CDEF_SEARCH_BEST_REF              1 // Limit CDEF search to the best filter strengths of the nearest ref frames
#define FTR_ADJUST_SR_FOR_STILL               1 // Adjust me/hme for still area

#define TUNE_SC_M7_M10                        1 // tune M7-M10 for screen content, cfl, spatial_sse, prehme, hme_level2

#define OPT_NOISE_LEVEL                       1 //estimate TF noise level for only I pictures
#define TUNE_SKIP_INTRA_PD1                   1 // set skip_intra to 1 for non-ref
#define TUNE_PREHME_M10                       1 // tuning prehme level for m10 and above
#define FTR_BIAS_STAT                         1 // Bias me stat
#define SS_ITXF                               1 //avoid doing itx when skip intra, filtering is off
#define FTR_INTRA_DETECTOR                    1 // infrastructure for selective intra and dlf
#define FTR_SELECTIVE_DLF                     1 // Enable DLF for base pictures and disable it selectively for other pictures
#define OPT_RDOQ_EOB                          1 // Control RDOQ f(eob)
#define FIX_SUBRES                            1 // Lossless bug fix (to link to the original macro FTR_SUBSAMPLE_RESIDUAL)

#define TUNE_M10_M9_1                         1 // tuning m10

#define DOWNSAMPLE_2D_AVX2                    1 //avx2 implementation of downsample_2d()
#define SS_CLN_MVP_TABLE                      1 // Cleanup and optimize MVP table generation; store all MVPs under ed_ref_mv_stack
#define SS_CLN_ME_CAND_ARRAY                  1 // create mrp off path for ME candidate array generation
#define TUNE_TXT_LEVEL                        1 // adding new level for txt_level
#define FTR_SELECTIVE_MFMV                    1 // Selectively disable MFMV

#define REMOVE_CLOSE_MVS                      1 // faster redundant bipred removal
#define FTR_ADJUST_SR_USING_LIST0             1 // adjust me sr for all pics if first ref pic in list0 has a small sad
#define FTR_RDO_OPT                           1 // Reduce the eob to speedup rdoq.
#define CLN_REMOVE_CHECK_MV_VALIDITY          1 // removing calls to check_mv_validity
#define SS_OPT_SUBPEL_PATH                    1 // Lossless opts to subpel path
#define OPT_SUBPEL_SKIP_TH                    1 // change early_exit_th used in subpel to depend on block size, not resolution
#define SS_OPT_EST_REF_TYPE                   1 // lossless opt.
#define OPT_INTER_FAST_COST                   1 // Use light inter fast cost
#define FTR_SIMPLIFIED_DEPTH_REMOVAL          1 // Simplify depth_removal ctrls
#define TUNE_NEW_M11_2                        1 // new levels for m11
#define  SS_OPT_NA                            1 // Loosless optimization of neighbour array access
#if SS_OPT_NA
#define  OPT_NA_IFS                           1 //interpolation filter
#define  OPT_NA_SKIP                          1 //skip flag
#define  OPT_NA_MODE                          1 //intra luma mode
#define  OPT_NA_ISINTER                       1 //is inter mode
#define  OPT_NA_INTRA                         1 //intra recon samples
#endif
#define SS_OPT_MD                             1 // Opt MD path
#define FTR_SIMPLIFIED_MV_COST                1 // Simplified MV cost estimation
#define SS_OPT_PALETTE_COST                   1 // Optimize Palette fast cost computation
#define SS_OPT_INTRABC                        1 // Remove unnecessary hash table initializations

#define OPT_MVP_READ                          1 // MVP selection
#define OPT_COMP_MODE_CHECK                   1 // Compound type validity check
#define CLN_UPDATE_CDF                        1 // bypassing unnecessary mv copy, syntax update and coeff update when update_cdf is off.

#define OPT_TXT_COPY                          1 // Copied txt search data to cand_buffer for DCT  case, then bypassed txt_data_cpy if DCT is the winner

#define TUNE_M10_M3_1                         1 // tuning m10-m3
#define OPT_MEMORY                            1 // Optimize memory footprint
#if OPT_MEMORY
#define OPT_MEMORY_REST         1  // RestProcess
#define OPT_MEMORY_CDF          1  // CDF
#define OPT_MEMORY_HASH_TABLE   1  // hash_table
#define OPT_MEMORY_MIP          1 //mip
#endif

#define FIX_PRESET_TUNING       1

#define SS_CLN_REF_PRUNE                      1 // Cleanup/optimize MD reference pruning
#define SS_OPT_MDS0                           1 // Optimize MDS0 scratch buffer sorting
#define SS_OPT_TXT_OFF_PATH                   1 // Create special transform path for TXT/TXS off

#define OPT_ESTIMATE_REF_BITS                 1 // Moved the redundant context derivation(s) @ estimate_ref_frames_num_bits() to be done per block instead of per frame-type
#define OPT_IFS                               1 // Speed up interpolation filter search
#define SS_OPT_TPL                            1 // Loosless optimization of TPL
#define SS_OPT_SET_LAMDA                      1 // Optimize labmda update for when block tuning is used
#define OPT_EARLY_CAND_ELIM                   1 // Upgrade early_cand_elimination feature
#define SS_CLN_LIGHT_PD0_PATH                 1 // Remove unnecessary opts from light-pd0 path
#define FIX_LIGHT_PD0_BEST_COST               1 // Correctly set MDS0 best cost for light-PD0
#define SS_OPT_SORTING                        1 // bypass fast cost sorting if MDS1 has only 1 candidate (can merge with SS_OPT_MDS0)

#define SHUT_8x8_IF_NON_ISLICE                1 // shut 8x8(s) if non-ISLICE

#define SS_OPT_INIT                           1    // Loosless optimization of init time (there is an R2R for M3 and below)
#if SS_OPT_INIT
#define  SIM_OLD_TPL 0
#endif
#define  FIX_TPL_NON_VALID_REF                1  // exclude TPL non valid ref (no TPL recon data available)

#define TUNE_M11_2                            1 // tuning m11
#define FTR_COEFF_DETECTOR                    1 // use detector based on the skipped area % in the reference picture
#define FTR_PREHME_OPT                        1 // reduce the complexity of the prehme list1 search based on list0 results

#define OPT_UPGRADE_TF                        1 // Added a 64x64 compensation path @ TF
#define TUNE_M9_SLOW                          1 // slowing down m9
#define FTR_FASTER_CFL                        1 // faster cfl
#define FTR_CDEF_BIAS_ZERO_COST               1 // Bias the cost of the (0,0) filter strength for CDEF
#if OPT_UPGRADE_TF
#define OPT_EARLY_TF_ME_EXIT                  1 // Early ME_TF exit; exit after HME if low distortion, then call the 64x64 path
#endif
#define TUNE_SC_ME_SR_ADJUST                  1 // set_me_sr_adjustment_ctrls to most aggresive level for sc for all presets

#define SS_OPT_CDEF_APPL                      1 // Save CDEF directional search results from search to use in application; remove unneeded copies in application
#define TUNE_SUBRES_TX                        1 // New settings for M11 subres tx
#if OPT_UPGRADE_TF
#define OPT_TUNE_DECAY_UN_8_10_M11            1 // Tune the decay value for fast_tf_filtering, use a more agressive abs_th to exit search, unify TF for 8Bit and 10BIT for M11: both to use the same search and filtering methods
#endif
#define FTR_16X16_TPL_MAP                     1 // Use a 16x16 based TPL map
#define SS_MORE_OPT_TPL                       0 // Loosless optimization of TPL
#define FTR_HME_ME_EARLY_EXIT                 1 // Uss zz sad  to early exit prehme/hme/me.
#define OPT_M11_PME                           1 // M11 level for PME
#define TUNE_ADD_SUBRESS_FACTOR4              1 // Add the ability to use a factor of 4 @ subres
#define TUNE_M10_M0                           1 // tuning m10-m0
#define FTR_TUNE_PRUNING                      1

#define OPT_M11_SUBPEL                        1 // Do not perform 1/4-Pel if the 1/2-Pel-to-fp-Pel error deviation is not high (could be also used to skip 1/8-Pel; skip next round if current round-to-previous-round deviation is not high) and use the central error to skip sub-Pel search
#define OPT_NIC_PRUNING                       1 // Skip class pruning step for best class; allows class pruning TH to be 0
#define TUNE_MDS0_SAD                         1 // Tune MDS0 metric for high presets in light-PD0
#define FTR_LIMIT_ME_CANDS                    1 // Add ability to limit the number of ME cands when MRP is off
#define OPT_FAST_MDS0_PD0                     1 // Skip distortion calc in light-PD0 MDS0 if only one cand
#define FTR_REDUCE_UNI_PRED                   1 // Reduce uni pred candidate
#define OPT_SUPEL_VAR_CHECK                   1 // Exit subpel search if the variance of the full-pel predicted block is low (i.e. where likely interpolation will not modify the integer samples)
#if OPT_UPGRADE_TF
#define FIX_SVT_POSITION_CHECK_CPP            1 // Fix the cpp error caused by svt_check_position()
#endif

#define OPT_DEPTH_REMOVAL_I_SLICE             1 // 1st depth_removal for I_SLICE using the variance of 64x64(s)
#define OPT_USE_INTRA_NEIGHBORING             1 // Fast INTER @ PD1 if left and above are coded as INTRA
#define TUNE_TXS_M11                          1 // tuning txs for m11
#define OPT_CDEF                              1 // Use skip are percentage ofthe ref to disable cdef
#define TUNE_M7_11                            1 // Tuning M7 to 11 at the end of May
#define OPT_TPL_64X64_32X32                   1 // Add the ability to search 32x32 block(s) and 64x64 block(s) @ tpl_dispenser(), and the ability to subsample the residual by 4 @ tpl_dispenser().
#if OPT_TPL_64X64_32X32
#define OPT_TPL_ALL32X32            1
#define OPT_TPL_ALL64X64            0
#endif
#define OPT_COMBINE_TPL_FOR_LAD               1 // Harmonize tpl for lad_mg 0 and lad_mg 1

#define SS_CLN_INIT_IFS_MDS0                  1 // Move IFS init to MDS0
#define LIGHT_PD1_MACRO                             1 // Add a light-PD1 path
#define CLN_DECPL_TX_FEATS                    1 // Decouple mds1 skipping tx shortcuts and reduce_last_md_stage_candidate

#define FTR_PD0_OPT                           1 // optimize pd0

#define SS_MEM_VAR                            1 // Optimize the variance array
#define SS_MEM_HIS                            1 // Optimize the histogram array
#define SS_MEM_DLF                            1 // Run time FDL scratch buffer

#if LIGHT_PD1_MACRO
#define FTR_LPD1_DETECTOR           1 // Add a detector for using light-PD1.  Set signal before depth refinement so you can use the detector to force pred depth
#define OPT_LPD1_MRP                1 // Enable light-PD1 for BASE (when WM and MRP are on).  NB light-PD1 does not do reference pruning
#define OPT_LPD1_PME                0 // Add PME to light-PD1
#endif

#define SS_MEM_TPL                            1 // Optimize TPL when lad_mg is zero
#define TUNE_M8_M10                           1 // tuning m8-m10 using lad_mg = 0 and light_pd1 on
#define FIX_LAD_MG_0_HANG                     1 // fixing the hang in M11 when lad_mg = 0, the original fix was using FTR_USE_LAD_TPL as the macro
#define FTR_INCREASE_QP_BY_8                  1 // Increasing QP by 8
#define CLN_PME_BLK_SIZE                      1 // Clean up PME block size modulation
#define OPT_USE_MVS_LPD1_DETECTOR             1 // Use MV length in light-PD1 detector
#define OPT_REMOVE_TXT_LPD1                   1 // Remove TXT from light-PD1
#define TUNE_M7_M8                            1 // tuning m8-m7 to make them faster
#define FTR_SKIP_COST_ZERO_COEFF              1 // Skip cost computation if zero luma coeffs
#define TUNE_M11_SUBPEL                       1 // Tune M11 subpel level
#if FTR_SKIP_COST_ZERO_COEFF
#define CHROMA_DETECTOR                       1 // Add detector for chroma complexity when skipping chroma TX/compensation
#endif
#define TUNE_M9_M10                           1 // tuning m10-m9 to make them faster
#define CLN_MISC_CLEANUP                      1 // Do some miscellaneous cleanup
#define ADJUST_LAMBDA                         1 // Adjust lambda based on deltaqp
#define OPT_REDUCE_COPIES_LPD1                1 // Reduce copying recon in LPD1
#define CLN_READ_REFINE_MVS                   1 // Create read_refine_me_mvs() func for LPD1
#define FIX_LAD_MG_WHEN_NO_TPL                1 // Set lad_mg to 0 if no TPL
#define FIX_DO_NOT_TEST_CORRUPTED_MVS         1 // Ignore corrupted MV(s) @ MD
#define FTR_6L_QPS                            1 // Improve the QPS offsets for 6L (used when TPL is OFF)
#define FTR_TF_STRENGTH_PER_QP                1 // Use the picture QP @ the strength derivation of TF
#define OPT_SKIP_SUBPEL_ZZ                    1 // Skip subpel for (0,0) MVs
#define OPT_EARLY_ELIM_TH                     1 // Add ability to control TH used for cand elimination
#define OPT_LIGHT_PD1                         1 // Use a more aggressive light_pd1 classifier

#define OPT_PA_REF                            1 // Memory Optimization of PA REF - remove the storage redundancy with input luma
#define TUNE_REMOVE_CDEF_COST_BIAS            1 // set cost bias to zero for cdef
#define TUNE_MDS0                             1 // Tune max number of candidate
#define FIX_HME_ME_EARLY_EXIT                 1 // Fix hme_me_early_exit r2r
#define OPT_LIGHT_PD1_USE_ME_DIST_VAR         1 // Force the use of PD1 if the variance of the 8x8 distortions is high
#define TUNE_M8_M11_MT                        1 // tune M8-M11 using multi-threaded framework
#define TUNE_4K_M8_M11                        1 // tuning 4k for m8-m11
#define FTR_TPL_SYNTH                         1 // Make TPL sythn to support 32x32
#define ME_8X8                                1 // memory optimization for me data arrays when 8x8 blocks are disallowed
#define TUNE_M7_MT                            1 // tune m7 using multi-threaded framework
#define FIX_TEMPORAL_FILTER_PLANEWISE         1 //Fix calculation for hbd fast, and mismatch exp C and AVX2

// DYNAMIC GOP
#define FIX_2PASS_CRF                         1 // Fix 2pass crf
#define CLIP_BASED_DYNAMIC_MINIGOP            1 // Ajust minigop size based on the first pass statistics for each clip
#define GOP_BASED_DYNAMIC_MINIGOP             1 // Ajust minigop size based on the first pass statistics at gop level
#define OPT_FIRST_PASS                        1 // Reduce the compexity of the first pass
#define OPTIMIZE_L6                           1
#define FIX_DATA_RACE_2PASS                   1 // fixing data race issues with the 2pass

#define TUNE_M6_MT                            1 // tune m6 using multi-threaded framework
#define OPT_UPDATE_MI_MAP                     1 // opt update_mi_map() and some signal setting
#define SS_CLN_10BIT_TPL_BUFFER                  1 // memory optimization for 10-bit by using 8-bit tpl buffer for 10-bit
#define SS_CLN_NOISE_DENOISE_BUFFERS             1 // stopped allocating memory for noise and denoise buffers as they were unused
#define FTR_SKIP_TX_LPD1                      1 // Skip luma TX in LPD1 based on neighbour info, distortion/QP, and candidate type
#define FTR_MVP_BEST_ME_LIST                  1 // Inject MVP unipred cands for best PA_ME ref only in LPD1
#define TUNE_REDUCE_NR_LPD1                   1 // Opt near count level for LPD1
#define OPT_TX_SKIP                           1 // Add new level for LPD1 TX skipping
#define OPT_MD_SEARCHES_10BIT                 1 // use 8-bit path for 10bit content when performing MD searches

#define FTR_10BIT_MDS3_REG_PD1                1 // Add ability to bypass EncDec for 10bit content (incl. when using 8bit MD)
#define FTR_10BIT_MDS3_LPD1                   1 // Add ability to bypass EncDec for 10bit content when using 8bit MD (LPD1 path) - to merge with above macro

#define SS_OPT_TPL_REC                           1 // Lossless Optimize memory for TPL recon buffers
#define CLN_RTIME_MEM_ALLOC                   1 // put run time memory allocation calls into function prefixed with rtime_alloc_...

#define TUNE_M7_M8_2                          1 // tuning m7 and m8
#define TUNE_4K_M11                           1 // tuning M11 for 4k

#define FIX_PA_REF_RELEASE_HANG               1 // Fix a hang that occured for incomplete MGs when MRP is ON
#define OPT_IBC_HASH_SEARCH                   1 // optimize hash search memory and speed by turning off higher block sizes and 4x4 for higher presets

#define TUNE_M8_M10_4K_SUPER                  1 // tuning M8 to M10 for 4K, july2021

#define OPT_BYPASS_ED_10BIT                   0 // Remove unnecessary buffer when copying recon for 10bit bypass-encdec; to merge with FTR_10BIT_MDS3_REG_PD1
#define FIX_10BIT_R2R                         1 // Fix 10bit r2r for bypassing encdec - merge with FTR_10BIT_MDS3_LPD1
#define FIX_COST_CALC_CHECK                   1 // Fix check to skip cost calcs
#define FIX_QUANT_COEFF_BUFF                  1 // Fix how the quant coeff buffer is released, and cleanup the init

#define SANITIZER_FIX                         1 // Fix thread sanitizer: race:svt_memcpy_small and race:svt_memcpy_sse
#define FIX_SANITIZER_RACE_CONDS              1 // Fix race conditions in MD/EncDec

//svt-05-03 start

#define FIX_QPS_OPEN_GOP                      1 // Use the BASE QP-Offset for CRA (instead of using the IDR QP-Offset)
#define FIX_TF_OPEN_GOP                       1 // Use past frame(s) to tf CRA (instead of using future frame(s) only)
#define FIX_PD0_RDOQ                          1 // Enable RDOQ @ mds3 of PD0
#define FIX_TF_HME                            1 // Apply a penalty to the HME_MV cost(@ the post - HME(0, 0) vs.HME_MV distortion check) when the HME_MV distortion is high (towards more search ~ (0,0) if difficult 64x64
#define FIX_SKIP_COEFF_CONTEXT                1 // Use skip-coeff-context update @ full-cost derivation
#define TUNE_HME_LEVEL0_M0                    1 // Tune HME-Level0 for 4K
#define TUNE_PME_M0                           1 // Tune Full-Pel PME

#define OPT_MI_MAP_MEMORY                     1 // Reduce fields in ModeInfo struct.  Only allocate MI for 8x8 blocks when 4x4 is disallowed.  Sub-macros can be merged with this one.
#if OPT_MI_MAP_MEMORY
#define OPT_INTRA_MI_MEM 1
#define OPT_INTER_MI_MEM 1
#define OPT_MODE_MI_MEM  1
#define OPT_PALETTE_MEM  1
#define OPT_TX_MI_MEM    1
#endif

#define SS_2B_COMPRESS                       1 // compress 10bit pictures into 10 bits from 16 bits, lossless, memory saving
#if SS_2B_COMPRESS
#define INC_PAD68                            1  // for 10bit increase the pad of source from 68 to 72 to be mutliple of 8 to accomodate 2bit-compression flow
#endif
#define FIX_TPL_R2R_BUFF                      1 // Fix r2r from using uninit'd recon buffer for REF pics
#define FIX_RACE_CONDS                        1 // Fix race conditions, causing sanitizer thread test failure
#define OPT_TF_8BIT_SUBPEL                    1 // Use 8bit pic in TF subpel
#define SS_TF_OPTS                           1 // lossless tf optimizations, turning off chroma memory allocation in higher presets
#define TUNE_LPD1_DETECTOR                    1 // Adjust LPD1 detector values for M11
#define OPT_10BIT_MDS3                        1 // Use 8bit for MDS3 when bypassing EncDec if recon is not needed
#define OPT_COEFF_BIT_EST                     1 // Skip coeff bit estimation when have many coeffs
#define FIX_PREHME_ADD                        1 // fix the addition of the prehme candidate into the hme_level0 candidates
#define FTR_TX_NEIGH_INFO                     1 // Use neighbour blocks' info to apply TX shortcuts

#define TUNE_M7_M8_3                          1 // tune m7-m8
#define FIX_LPD1_DETECTOR                     1 // Fix left SB check for LPD1 detector when PD0 is skipped

#define FIX_COMPRESSED_10BIT                  1 // compressed 10bit files working (use with flag --compressed-ten-bit-format)
#define OPT_TX_SHORTS                         1 // Extend coeff bit est skipping to regular PD1 path; use TX shrotcuts when PD0 is skipped
#define FIX_TPL_BOOST                         1 // fix r2r between intra-period -1 and infinit intra-period.
#define FIX_ISSUE_45                          1 // Issue 45 fix

#define TUNE_M9_11_3                          1 // tune m8-11, new levels for deep_removal_level and md_pme_level. Small changes in M8 to improve spacing.

#define FIX_INTRA_PERIOD_2PASS                1 // fixing the issue of intra-period -1 being lossy against intra-period greater than the number of frames in 2pass

//svt-05-04 start

#define FIX_ISSUE_46                          1 // disbale 2pass when I period is less than 16
#define FIX_DO_NOT_TEST_CORRUPTED_MV_DIFF     1 // Ignore corrupted MV-to-MVP diffrence @ MD

#define FTR_MEM_OPT                           1 // Remove reference buffer duplication in 10bit
#define FTR_MEM_OPT_WM                        1 // Remove reference buffer duplication in 10bit for Warped motion

#define TUNE_SUBRES_LEVEL                     1 // New level for subres in M11 LPD0
#define OPT_MALLOC_TRIM                       1 // use malloc trim tool to release freed memory  back to the OS
#define OPT_PD0_PF_LEVEL                      1 // Improve detector for using PF in LPD0

#define FIX_ISSUE_45_M0_6                     1 // fix R2R due to enable_restoration checks to sc content
#define OPT_TXS_WM                            1 // Add new TXS level restricting TXS usage per sq size; tune M11 WM level

#define OPTIMIZE_COMPRESS_PACK_SB             1 //Optimize svt_compressed_packmsb_avx2_intrin() to work with any width
#define OPTIMIZE_SVT_UNPACK_2B                1 //AVX2 Implementation of kernel svt_unpack_and_2bcompress_c()
#if FTR_MEM_OPT_WM
#define FIX_UT_FTR_MEM_OPT_WM                 1 // Fix Unit Tests for macro FTR_MEM_OPT_WM
#endif


#if FIX_TEMPORAL_FILTER_PLANEWISE && FTR_TF_STRENGTH_PER_QP && SIMD_APPROX_EXPF
#define FIXED_POINTS_PLANEWISE                1 //Calculate Temporal Filter Planewise on Fixed Points
#else
#define FIXED_POINTS_PLANEWISE                0 //Calculate Temporal Filter Planewise on Fixed Points
#endif
#if FIXED_POINTS_PLANEWISE
#define ENABLE_FIXED_POINTS_PLANEWISE         1 //Enable Temporal Filter Planewise on Fixed Points. Need decide on what level should be used.
#define ENABLE_MEDIUM_PLANEWISE               0 //Enable Temporal Filter Planewise Medium on Fixed Points instead of normal filter. Need decide on what level should be used.
//#define FIXED_POINT_ASSERT_TEST             1 //Enable special FP_ASSERT () tests for all builds. No define macro, then test only for debug
#endif

#define FIX_VBR_R2R                           1 // Fixes the 1 PASS VBR run to run
#define TUNE_LPD1_DETECTOR_LVL                1 // Tune the LPD1 detector thresholds for M11

#define TUNE_MEDIUM_TFILTER                   1 // Tunes the new medium TF filters

#define OPT_MMAP_FILE                         1 // Add support for memory mapped files in linux
#define TUNE_MRP_MEM                          1 // Tune memory for  MRP level 3
#define TUNE_BYPASS_MEM                       1 // Tune memory for  BYPASS ENCDEC OFF

#define FIX_DG                                1 // fixing data race issues with the 2pass without impacting the bdrate and enable//CI
                                                // the 2pass crf with DG for up to M10
#define OPT_FIRST_PASS2                       1 // Reduce the compexity of the first pass//CI
#define OPT_FIRST_PASS3                       1 // Bypass the odd sbs in both directions.//CI
#define OPT_FIRST_PASS4                       1 // Skipintra comp and Downsampling.//CI

#define OPT_TPL_DATA                          1 // move the TPL data from the parent pcs to the me results buffer
#define CLN_MDCONTEXT                         1 // clean up pointers in the ModeDecisionContext struct
#define OPT_MEM_PALETTE                       1  // make palette_info runtime memory allocation

#define OPT_1P                                0 //New encoder pipeline for IPP based first pass //CI
#define FIX_ISSUE_50                          1
#define FIX_ISSUE_49                          1

#define ENBLE_SKIP_FRAME_IN_VBR_MODE          1 // Skip frame step of 8 for VBR.
#define ENBLE_DG_IN_VBR_MODE                  1 // Enable DG for VBR.

#define FIX_TXB_INIT_LEVELS                   1 // Fix EC to use AVX2 @ svt_av1_txb_init_levels()
#define FTR_VLPD1                             1 // Add a new very-light-PD1 path, incl. new levels for existing features, and tracking info from co-located SBs

#define FTR_MOD_DEPTH_REMOVAL_LVL             1 // Apply a more agressive depth_removal level for Layer0 based on the qp_offset band
#define FIX_UNPACK_10BIT_ED_BYPASS            1 // Fix how the 10bit recon pictures are unpacked when generating the 8bit recon when ED is bypassed

#define OPT_PREHME                            1 // Use TF motion to direct preHme search
#define TUNE_EXIT_TH                          1 // Skip ME L1 ref based on L0 ref.
#define TUNE_VLPD1_LVLS                       1 // Improve the features used in VLPD1
#define TUNE_REG_PD1                          1 // Tune feature levels used in reg. PD1
#define FIX_LPD1_FOR_ISLICE                   1 // Make LPD1 and VLPD1 compatible with I_SLICE

#define TUNE_M9_M11_OPTIMIZED_SUPER4KL        1 // Optimize 4K in M11 and slowing down M9 by around 1%
#define TUNE_M7_SLOWDOWN                      1 // Slowdown M7 for about 20% to get to a distance of around 82% from M6
#define TUNE_M8_SLOWDOWN                      1 // Slowdown M8 for about 29% to get to a distance of around 80% from the new M7 after slowdown
#define TUNE_M9_SLOWDOWN                      1 // Slowdown M9 for about 20% to get to a distance of around 73% from the new M8 after slowdown
#define TUNE_M10_SLOWDOWN                     1 // Slowdown M10 for about 16% to get to a distance of around 63% from the new M9 after slowdown
#define TUNE_M11_SLOWDOWN                     1 // Slowdown M11 for about 22% to get to a distance of around 51% from the new M10 after slowdown
#define TUNE_NEW_M12                          1 // Enabling M12 as the old M11

#define FIX_FTR_PREHME_SUB                    1 // Fix mismatch between svt_sad_loop_kernel_{c/sse4_1/avx2/avx512}_intrin introduced by FTR_PREHME_SUB macro

#define FTR_VLPD0                             1 // Very Light PD0

#define OPT_M12_SPEED                         1 // optimize subpel_me for higher resolution for speed and modify ftr_vlpd0 levels

#define FTR_LESS_BI                           0 // Avoid ME BiPred candidates if they have close enough MV to already injected Bi Pred.

#define FTR_M13                               1 // add M13 preset to the encoder, should match the cvh path
#define SS_CLN_CFL_CTRLS                      1 // Create ctrl struct for CFL
#define CLN_RATE_EST_CTRLS                    1 // Create ctrl struct for rate estimation signals
#define CLN_REMOVE_UNUSED_FEATS               1 // Remove unneeded features

#define FIX_INIT_ZZ_CAND                      1 // Fix initialized variables @ inject_zz_backup_candidate()
#define CLN_TF_CTRLS                          1
#define CLN_DEPTH_REMOVAL                     1
#define CLN_INDEPTH                           1
#define CLN_TF_LVL                            1
#define CLN_P_BASE_LVL                        1
#define CLN_HIGH_LVL_SYNTAX                   1
#define CLN_SUBPEL_CTRLS                      1
#define CLN_REF_PRUNING                       1
#define CLN_RDOQ_CTRLS                        1
#define FTR_VLPD0_INTER_DEPTH                 1 // Modulate the inter-depth bias based on the QP, and the temporal complexity of the SB towards more split for low QPs or/and complex SBs, and less split for high QPs or/and easy SBs (to compensate for the absence of the coeff rate).

#define CLN_RED_CAND                          1 // Remove some actions of reduce_last_md_stage_candidate (actions related to RDOQ, IFS, and class pruning)
#define TUNE_M10_M11                          1 // tuned m10 and m11 oct 6
#define SS_FIX_MOVE_SKIP_INTRA_PIC            1 // Change skip_intra to be set at the picture level
#define CLN_LPD1_LVLS                         1 // Make LPD1 have several levels, instead of enc_mode checks per feature
#define CLN_REG_PD1_DIFFS                     1 // Remove feature level diffs from reg. PD1 path for high presets, where LPD1 is used mostly
#define CLN_TPL_LEVEL_7                       1 // Remove TPL level 7 including modulate_depth_removal_level flag
#define CLN_LPD1_TX_CTRLS                     1 // Clean up TX shortcut controls used in LPD1
#define CLN_LPD0_CTRL                         1
#define CLN_SUBRES                            1
#define FTR_NEW_QPS                           1 // Added the ability to use the libaom QPS model (only CQP 6L)
#define OPT_PACK_HBD                          1 // Avoid redundunt packing of high bit depth source pixels between different PD passes. Lossless.
#define CLN_INTRA_CTRLS                       1 // Cleanup intra controls
#define TUNE_HBD_MD                           1 // Add checks for M9-13, hybrid 10bit/8bit MD for INTRA, 8bit MD for NON-INTRA.
#define CLN_MERGE_LPD0_VLPD0                  1
#define CLN_REG_PD1_TX_CTRLS                  1 // Cleanup controls for TX shortcuts in reg. PD1
#define CLN_M10_M12_DIFFS                     1 // Remove useless feature diffs in M10-M12
#define CLN_NIC_PRUNE_CTRLS                   1 // Merge mds1_skip_level with nic_pruning ctrls

#define CLN_REMOVE_CORRUPTED_MV_PRINTF        1
#define FIX_I64                               1 // Remove tile use from 1st pass
#define FIX_I51                               1 // Fix release of PA references for P pictures
#define FIX_R2R_TPL_IXX                       1 // Move TPL control from ress coord to pic decision where pcs->hierarchical_levels is finally set.
#define CLN_MDS0_CTRLS                        1 // Cleanup ctrls for MDS0
#define CLN_REF_AREA                          1 // Make ref_is_high_skip and ref_is_high_intra return a percentage (out of 100)
#define CLN_LIST0_ONLY_BASE_IFS_MFMV          1
#define CLN_NIC_SIGS                          1 // Merge NIC pruning/scaling controls under one signal
#define DIS_VBR_HL0                           1 // Block VBR in HL0
#define CLN_M6_M12_FEATURES                   1 // push all small speed features out to m11 vs m12 difference
#define FTR_16K                               0 // Support Encoding of up to 16K resolutions
#define CLN_FEAT_LEVEL                        1 // Remove useless feature levels/diffs

#define CLN_CAND_REDUCTION_CTRLS              1 // Merge near_count_level, lpd1_mvp_best_me_list, reduce_unipred_candidates, redundant_cand_level, eliminate_candidate_based_on_pme_me_results, use_neighbouring_mode, merge_inter_classes into one control.

#define CLN_USE_SELECTIVE_MFMV_TH             1
#define TUNE_VBR_OVERSHOOT                    1 // Tune VBR to reduce the overshoot
#define RFCT_ME8X8                            1 // refactor the enable_me_8x8 and disallow_below_16x16 code so that the two features are linked
#define CLN_UPDATE_CDF                        1 // Cleanup levels for update_cdf feature in MD

#define FIX_TF_FILTER_64x64_PATH              1
#define TUNE_SC_DETECTOR                      1 // Disable SC detector for 4k and higher resolutions

#define CLN_TF                                1
#define CLN_MD_STAGING_CTRLS                  1
#define CLN_SUBPEL_LVL                        1
#define IPP_CTRL                              1

#define TUNE_IMPROVE_M12                      1 // Pushing out M11 M12 differences to M13, slow down by 0.5%, gain 0.5% BDrate.

#define TUNE_4K_STABILITY                     1 // Change 4k levels of depth_removal_level and pd0_level level because of R2R
#define TUNE_M13_LPD1_LVL                     1 // Turn LPD1 off for BASE pics in M13

#define TUNE_IMPROVE_M11_M10                  1 // Improve spacing and slope for M11 and M10, M11 slowdown by 2% gain 1.5% BDrate, M10 speedup 3% lose .3% BDrate. tpl_lad_mg 1 to 0 is moved from M8 level to M11 level.
#define TUNE_M1_M8                            1 // tune M1 to M8 after cleanup
#define FIX_VBR_R2R                           1 // Fix run-to-run in 3 pass VBR
#define TUNE_LAD_MAX                          1 // Change the max LAD size to 120
#define FIX_SUBRES_R2R                        1 // Fix uninit'd subres signals

#define CLN_4K_CHECKS                         1

#define CLN_REG_PD_SIG                        1
#if CLN_REG_PD_SIG
#define CLN_REG_PD_SIG_SET_0                  1
#define CLN_REG_PD_SIG_SET_1                  1
#define CLN_REG_PD_SIG_SET_2                  1
#endif
#define TUNE_ED_BYPASS_8BIT                   1 // Tune 8bit bypass_encdec level
#define FIX_THREE_QUAD_ENERGY                 1 // Initialize three_quad_energy to 0 in perform_dct_dct_tx_light_pd1

#define CLN_SB_LEVEL_SIG                      1
#if CLN_SB_LEVEL_SIG
#define CLN_WM_SIG                            1 // centralize and clean up Warped motion signalling - Lossless
#define CLN_4X4_SIG                           1 // Move 4x4 signal from SB to Pic Level - Lossless
#define CLN_BYP_ED_SIG                        1 // Move bypas enc dec     signal from SB to Pic Level - Lossless
#define CLN_PD0_SIG                           1 // Move pd0_level signal from SB to Pic Level - Lossless
#define CLN_SKIP_PD0_SIG                      1 // Move skip_pd0 signal from SB to Pic Level - Lossless
#define CLN_DISALLOW_BELOW_16X16_SIG          1 // Move disallow_below_16x16 signal from SB to Pic Level - Lossless
#define CLN_DEPTH_REMOVAL_SIG                 1 // Move depth_removal_level signal from SB to Pic Level - Lossless
#define CLN_BLOCK_BASED_DEPTH_SIG             1 // Move block_based_depth_refinement_level signal from SB to Pic Level - Lossless
#define CLN_LPD1_LVL_SIG                      1 // Move lpd1_lvl signal from SB to Pic Level - Lossless


#define IMPROVE_PRESET_SPACING                1 // Improve spacing between the presets
#if IMPROVE_PRESET_SPACING
#define CLN_RES_CHECKS                        1
#define CLN_RES_CDEF                          1
#define CLN_RES_TXT                           1
#define CLN_RES_ME                            1
#define CLN_RES_IFS                           1
#define CLN_RES_DP                            1
#define CLN_RES_DISALLOW_B16                  1
#define CLN_RES_DEPTH_REMOVAL                 1
#define CLN_RES_ME_BIS                        1
#define CLN_RES_LPD1_BIS                      1
#endif

#endif

#define TUNE_M9_M13                           1 // tune both M9 and M13 for better slopes, will affect stream for M6-M8 also
#define TUNE_2PASS_SETTINGS                   1 // move 1 pass to M9 for crf
#define TUNE_PICT_PARALLEL                    1 // Tune picture parallelization
#define FIX_TPL_PORTS                         1 // Separate the TPL thread port numbers from the EncDec port numbers
#define FTR_MG_PARALELL                       1  // MG level picture based pralellism
#define FIX_ED_PORT                           1  // Fix unused Enc-Dec port FIFO causing a seg fault
#define FIX_RC_PORT                           1  // Fix RC port FIFO causing a seg fault
#define FIX_PMG_PORT                          1  // Fix Pic-Mgr port FIFO causing a seg fault

#define CLN_MATHUTIL                          1 //Cleanup linsolve(), and more
#define CLN_RANSAC                            1 //Cleanup ransac operations
#define OPT_FILM_GRAIN                        1 //lossless optimization of film grain feature
#define OPT_CORNER_MATCH                      1 //lossless optimization of corner match
#define OPT_FP_LOG                            1 //lossless optimization of first pass
#define OPT_CODE_LOG                          1 //lossless optimization of remove few log()
#define CLN_RATE_CONTROL                      1 //Cleanup RateControlProcess
#define CLN_MD_MEAN_CALC                      1 //Cleanup ModeDecision and RateControlProcess geometrin mean calc
#define OPT_SBO_CALC_FACTORS                  1 //lossless optimization of SourceBaseOperations
#define OPT_CALC_TPL_MC                       1 //lossless optimization of SourceBaseOperations
#define CLN_2PASS                             1 //lossless cleanup 1 and 2 PASS structures and calculations
#define TUNE_CAPPED_CRF                       1 // Improve the Capped CRF algorithm

#define TUNE_SC_SPACING                       1 // Tune SC features for presets
#define TUNE_10BIT_M5_M8                      1 // tune 10bit setting for m5-m8
#define FIX_TILES_COL                         1 // Fix the hard-coded Tiles setting
#define FIX_I87                               1 // Fix issue 87
#define FIX_I80                               1 // Fix issue 80
#define FIX_RDCOST_OVERFLOW                   1 //Fix overflow uint32_t when calculate RDCOST. Found by mismatch between GCC and CLANG. Fix issue 101
#define FIX_CLANG_GCC_MISMATCH                1 //Fix mismatch between GCC and Clang, issue 101

#define TUNE_MIDDLEP_VBR                      1 // Set the best middlepass preset for each encoder mode.
#define OPT_1P_VBR                            1 // Optimise 1pass VBR.
#define TUNE_RECODE_VBR                       1 // Set recode_loop level for VBR.#endif //----------------------------------- all svt-05 features should be place are above this line -------------------------

#define FIX_ISSUE_99                          1 // invalid parent-to-pred cost deviation @ the final Pred-Depth refinement stage (the up to 2 depth(s) only check) for incomplete SB(s) as not derived @ the regular refinement stage (weather to consider or not the parent depth); initialized to MAX before (and not inside) the regular refinement stage

#define TUNE_OVERSHOOT_I83                    1 // Tune VBR to decrease the overshoot
#endif //----------------------------------- all svt-05 features should be place are above this line -------------------------
#if !PRIVATE_MACROS

#undef LIGHT_PD0
#undef LIGHT_PD0_2
#undef  CLN_GEOM

#undef FTR_MULTI_STAGE_CDEF                   // multi_stages CDEF
#undef MS_CDEF_OPT1                           // Optimised skip cdef condition
#undef MS_CDEF_OPT2                           // Add the ability to reduce the number of filter strengths for chroma components
#undef MS_CDEF_OPT3

#undef CLN_REDUCE_CDEF

#undef TUNE_CDEF_TF_LEVELS

#undef FTR_CDEF_SUBSAMPLING                    // Allow subsampling (i.e. skipping rows) in the CDEF search
#undef FTR_SKIP_LINES_CDEF_SEARCH
#undef FTR_SKIP_LINES_CDEF_DIST

#undef OPT_REDUCE_CDEF_MEMSET
#undef SS_OPT_CDEF
#undef OPT_CDEF_LEVELS
#undef FTR_CDEF_SEARCH_BEST_REF

#undef SS_ITXF
#undef SS_OPT_MD

#undef SS_OPT_MDS0
#undef OPT_EARLY_CAND_ELIM
#undef SS_CLN_LIGHT_PD0_PATH
#undef FIX_LIGHT_PD0_BEST_COST
#undef SIM_OLD_TPL
#undef FIX_TPL_NON_VALID_REF
#undef TUNE_M11_2
#undef FTR_COEFF_DETECTOR
#undef OPT_UPGRADE_TF
#undef TUNE_M9_SLOW
#undef FTR_CDEF_BIAS_ZERO_COST
#undef OPT_EARLY_TF_ME_EXIT
#undef SS_OPT_CDEF_APPL
#undef FTR_16X16_TPL_MAP
#undef SS_MORE_OPT_TPL
#undef TUNE_ADD_SUBRESS_FACTOR4
#undef TUNE_M10_M0
#undef OPT_M11_SUBPEL
#undef TUNE_MDS0_SAD
#undef OPT_FAST_MDS0_PD0
#undef OPT_SUPEL_VAR_CHECK
#undef OPT_USE_INTRA_NEIGHBORING
#undef OPT_CDEF
#undef TUNE_M7_11
#undef OPT_TPL_64X64_32X32
#undef OPT_TPL_ALL32X32
#undef OPT_TPL_ALL64X64
#undef OPT_COMBINE_TPL_FOR_LAD
#undef LIGHT_PD1_MACRO
#undef FTR_PD0_OPT
#undef FTR_LPD1_DETECTOR
#undef OPT_LPD1_MRP
#undef OPT_LPD1_PME
#undef TUNE_M8_M10
#undef TUNE_M7_M8
#undef CLN_MISC_CLEANUP
#undef CLN_READ_REFINE_MVS
#undef FTR_TF_STRENGTH_PER_QP
#undef OPT_LIGHT_PD1
#undef TUNE_REMOVE_CDEF_COST_BIAS
#undef OPT_LIGHT_PD1_USE_ME_DIST_VAR
#undef TUNE_M8_M11_MT
#undef TUNE_4K_M8_M11
#undef FTR_TPL_SYNTH
#undef TUNE_M7_MT

#undef TUNE_M6_MT
#undef ME_8X8
#undef CLN_RTIME_MEM_ALLOC
#undef OPT_IBC_HASH_SEARCH
#undef FTR_10BIT_MDS3_REG_PD1
#undef FTR_10BIT_MDS3_LPD1
#undef TUNE_4K_M11
#undef TUNE_M8_M10_4K_SUPER
#undef FTR_TX_NEIGH_INFO
#undef TUNE_M7_M8_3
#undef TUNE_M9_11_3
#undef OPT_COEFF_BIT_EST
#undef OPT_TX_SHORTS
#undef FIXED_POINTS_PLANEWISE
#undef TUNE_MEDIUM_TFILTER
#undef FTR_MEM_OPT
#undef FTR_MEM_OPT_WM
#undef OPT_MEM_PALETTE
#undef TUNE_BYPASS_MEM
#undef OPT_TPL_DATA
#undef FTR_VLPD1
#undef FTR_VLPD0
#undef OPT_M12_SPEED
#undef SS_CLN_CFL_CTRLS
#undef CLN_RATE_EST_CTRLS
#undef CLN_TF_CTRLS
#undef CLN_SUBPEL_CTRLS
#undef CLN_LPD1_LVLS
#undef CLN_TPL_LEVEL_7
#endif
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

#endif // EbDebugMacros_h
