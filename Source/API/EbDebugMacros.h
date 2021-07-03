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
#define FTR_1PASS_CBR_RT       1 // to support 1pass CBR real time mode
#define FTR_1PAS_VBR                1 // 1 PASS VBR with ability to use LAD stats
                                      // Calculate stat total in 1 pass VBR based on the look ahead
#define FTR_LAD_INPUT               1 // Add look ahead as an input to the encoder
#define FTR_2PASS_1PASS_UNIFICATION 1 // Unifty the VBR algorithms between 2 pass and 1 pass
                                      // Create a ME functions for the IPPP processing of 2PASS and 1 PASS VBR/CBR

#define FIX_LOW_DELAY               1 // Fix the 6L, 5L and 4L with LOW_DELAY_P

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
#define  SIM_OLD_TPL 1
#endif
#define  FIX_TPL_NON_VALID_REF  0

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
#define LIGHT_PD1                             1 // Add a light-PD1 path
#define CLN_DECPL_TX_FEATS                    1 // Decouple mds1 skipping tx shortcuts and reduce_last_md_stage_candidate

#define FTR_PD0_OPT                           1 // optimize pd0

#define SS_MEM_VAR                            1 // Optimize the variance array
#define SS_MEM_HIS                            1 // Optimize the histogram array
#define SS_MEM_DLF                            1 // Run time FDL scratch buffer

#if LIGHT_PD1
#define FTR_LPD1_DETECTOR           1 // Add a detector for using light-PD1.  Set signal before depth refinement so you can use the detector to force pred depth
#define OPT_LPD1_MRP                1 // Enable light-PD1 for BASE (when WM and MRP are on).  NB light-PD1 does not do reference pruning
#define OPT_LPD1_PME                1 // Add PME to light-PD1
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
#endif

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
#undef LIGHT_PD1
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
#endif
//FOR DEBUGGING - Do not remove
#define NO_ENCDEC               0 // bypass encDec to test cmpliance of MD. complained achieved when skip_flag is OFF. Port sample code from VCI-SW_AV1_Candidate1 branch
#define DEBUG_TPL               0 // Prints to debug TPL
#define DETAILED_FRAME_OUTPUT   0 // Prints detailed frame output from the library for debugging
#define TUNE_CHROMA_SSIM        0 // Allows for Chroma and SSIM BDR-based Tuning

#define MIN_PIC_PARALLELIZATION 0 // Use the minimum amount of picture parallelization
#define TUNE_PICT_PARALLEL      0 //  Tune picture parallelization
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
