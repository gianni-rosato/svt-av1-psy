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

#define FTR_GOP_CONST_RC        1 // When enabled, the rate control matches the target rate for each GoP
                                  // Currently, only applicable for VBR and  when GoP size is greater than 119 frames.
#define CLN_B64_RENAMING        1 // Use b64 where appropriate rather than sb variant
#define EN_HL2                  1 // Enable and improve the performance of HL2
#define TUNE_HL2                1 // Tune warp, nics, supbel and cdef for HL2
#define FIX_LD_R2R              1 // Restrict the number of past frames to 1 in LD

#define FIX_PARTITION_COST      1 // Remove the extra Parition_None cost

#define FTR_BLOCK4x4_PER_SB     1 // skip 4x4-block(s) on the fly per SB at PD1 (when an 8x8 is selected) based on the variance of the 4 recon-to-src quadrants
#define OPT_DEPTH_REFINMENT     1 // use more conservative inter-depth deviation thresholds, and applied an additional offset (towards more aggressive pruning) when the child-depth is 4x4.
#define OPT_DEPTH_REMOVAL       1 // use the me-cost of the 8x8-depth, and the deviation of the 8x8 me-distortion(s) to control the 4x4-depth
#define OPT_PME                 1 // turn the problematic knobs and used more conservative settings for the safe knob(s)
#define OPT_INTRA_PD0           1 // use the me-cost of the 8x8-depth, and the 64x64 variance to control INTRA at PD0
#define OPT_M5                  1
#define OPT_MRP                 1 // Use TPL-stats at the pruning mechanism; derive the highest TPL reference-index per list and per SB, then apply an offset to the default pruning deviation threshold for the reference indices higher than the highest TPL reference-index
                                  // Shut MRP based on the block distortion(s) to the closest ref
                                  // Tune HME/ME, and MD pruning th(s)
#define FTR_MDS0_INTRA_OPT      1 // Process CLASS_0 candidates through iterations at mds0
#define OPT_GMV_M5              1 // Adaptively switch between GM_DOWN and GM_DOWN16 based on the average ME distortion, and the picture variance
#define OPT_TF_1                1 // tf using 64x64 MV(s)
#define OPT_TF_2                1 // shut subpel early exit
#define OPT_LAMBDA_MODULATION   1 // Improve lambda modulation
#define OPT_WN_FILTER           1 // Optimize wn filter in M5
#define OPT_SG_FILTER           1 // Optimize sg filter in m5

#define OPT_TF                  1 // Use medium filter for all resolutions in M0-M9
#define CLN_COEFF_BUFFER        1 // Cleanup the coeff buffer used in TX - should remove files
#define FTR_USE_SATD_TXT        1 // Use SATD to skip TX types in TXT search
#define FTR_MDS_PRUNE           1 // Scale MDS0 cand prune TH based on rank of current candidate compared to the best. Use RDOQ for INTRA candidates in MDS1
#define OPT_TF_SPLIT_TH         1 // Remove MAX/MIN 16x16 block distortion difference from TF split decision
#define OPT_REMOVE_RDOQ_FEAT    1 // Remove skip_rdoq_uv_qp_based_th
#define OPT_WM_INJ              1 // Add function for WM refinement
#define FIX_CHECK_OBMC_MVS      1 // Check that refined OBMC MVs are valid
#define OPT_MDS0_PRUNE          1 // Split MDS0 cand pruning TH into intra/inter THs
#define OPT_INTER_COMP          1 // Opt inter compound
#define OPT_SATD_PER_TYPE       1 // Make SATD early TXT exit inter/intra dependent
#define FTR_MDS1_PRUNE          1 // Optimize MDS1 candidate pruning
#define OPT_MDS3_TX_BYPASS      1 // Optimize bypassing TX in MDS3 when MDS1 has 0 coeffs
#define OPT_NSQ_M5              1 // Enable NSQ for I_SLICE in M5
#define TUNE_TPL_WEIGHTS        1 // Tune TPL weights; off for now. ~-0.25% BDR gain for Y-PSNR BDR; ~0.05% BDR loss for Y-SSIM BDR. To test for VBR as well
#define FIX_CAND_COUNT          1 // Increase mem_max_can_count for low presets
#define OPT_M4                  1 // Tune default M4 for SSIM
#define OPT_M3                  1 // Tune default M3 for SSIM
#define OPT_M2                  1 // Tune default M2 for SSIM
#define FTR_DEPTH_EARLY_EXIT    1 // Exit depth early if cost is > TH * higher_depth_cost
#define OPT_RDOQ_MDS1           1 // Perform RDOQ in MDS1 only if MDS2 is bypassed
#define OPT_M0_TXT              1 // Use SATD early exit TH of 20 for M0
#define FTR_USE_PD0_COEFF       1 // Use PD0 coeff info towards more aggressive depth refinement
#define OPT_M1                  1 // Tune default M1 for SSIM
#define OPT_M6                  1 // Tune default M6 for SSIM
#define REMOVE_PART_ASSERT      1 // For debugging: Temporarily remove an assert in svt_aom_partition_rate_cost() that is triggered when called from update_d1_data(). The
                                  // function should not be called for disallowed blocks, but currently there is no check (block is not evaluated at MD but still calls the func).
#define OPT_M7                  1 // Tune default M7 for SSIM
#define OPT_M5_2                1 // Tune default M5 for SSIM
#define OPT_M8                  1 // Tune default M8 for SSIM
#define OPT_M9                  1 // Tune default M9 for SSIM
#define OPT_M10                 1 // Tune default M10 for SSIM
#define OPT_M11                 1 // Tune default M11 for SSIM
#define OPT_M12                 1 // Tune default M12 for SSIM
#define OPT_M13                 1 // Tune default M12 for SSIM

#define FIX_ISSUE_1969          1 // Fix issue #1969
#define FIX_DISALLOW_CTRL       1 // Move get_disallow_4x4() and get_disallow_nsq() to a header

#define OPT_DAV1D_INV_TXFM      1 // Enable inverse transform functions from dav1d (written in AVX2 assembler)
#define OPT_DAV1D_BLEND_V_H     1 // Enable blend vertical/horizontal functions from dav1d (written in AVX2 assembler)

#define FIX_TYPO_IN_6L_RPS      1 // Fix a typo in the 6L RPS construction

#define OPT_COST_COEFFS_TXB     1 // Use coding trick to optimize cost calculation; branch misprediction reduction from ~30% to ~10%

#define FIX_MDS1_COUNT          1 // Fix mds1 count
#define FIX_CRA_R2R             1 // Fix r2r when using open-gop configuration caused by not initializing the cra_flag under pcs_ptr
#define FIX_SUPERRES_MEM_LEAK   1 // Fix memory leak when superres is used

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
