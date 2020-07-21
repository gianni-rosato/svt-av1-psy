/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbModeDecisionConfigurationProcess_h
#define EbModeDecisionConfigurationProcess_h

#include "EbDefinitions.h"
#include "EbModeDecision.h"
#include "EbSystemResourceManager.h"
#include "EbMdRateEstimation.h"
#include "EbRateControlProcess.h"
#include "EbSequenceControlSet.h"
#include "EbObject.h"
#include "EbInvTransforms.h"

#ifdef __cplusplus
extern "C" {
#endif

/**************************************
     * Defines
     **************************************/
static const uint8_t  depth_offset[4]   = {85, 21, 5, 1};
static const uint32_t ns_blk_offset[10] = {0, 1, 3, 25, 5, 8, 11, 14, 17, 21};
static const uint32_t ns_blk_num[10]    = {1, 2, 2, 4, 3, 3, 3, 3, 4, 4};

typedef struct MdcpLocalBlkStruct {
    uint64_t early_cost;
    EbBool   early_split_flag;
    uint32_t split_context;
    EbBool   selected_cu;
    EbBool   stop_split;
} MdcpLocalBlkStruct;

typedef struct ModeDecisionConfigurationContext {
    EbFifo *                 rate_control_input_fifo_ptr;
    EbFifo *                 mode_decision_configuration_output_fifo_ptr;
    MdRateEstimationContext *md_rate_estimation_ptr;
    EbBool                   is_md_rate_estimation_ptr_owner;
    uint8_t                  qp;

    // Adaptive Depth Partitioning
    uint32_t *sb_score_array;
    uint8_t   cost_depth_mode[SB_SQ_NON4_BLOCKS_DEPTH_MODE];
    uint8_t * sb_cost_array;
    uint32_t  predicted_cost;
    uint32_t  budget;
    int8_t    score_th[MAX_SUPPORTED_SEGMENTS];
    uint8_t   interval_cost[MAX_SUPPORTED_SEGMENTS];
    uint8_t   number_of_segments;
    uint32_t  sb_min_score;
    uint32_t  sb_max_score;
    uint32_t  sb_average_score;

    const BlockGeom *      blk_geom;
    ModeDecisionCandidate *mdc_candidate_ptr;
    CandidateMv *          mdc_ref_mv_stack;
    BlkStruct *           mdc_blk_ptr;
    uint8_t                qp_index;

    // Multi - Mode signal(s)
#if REMOVE_MR_MACRO
    EbEncMode adp_level;
#else
    uint8_t adp_level;
#endif
} ModeDecisionConfigurationContext;

/**************************************
     * Extern Function Declarations
     **************************************/
EbErrorType mode_decision_configuration_context_ctor(EbThreadContext *  thread_context_ptr,
                                                     const EbEncHandle *enc_handle_ptr,
                                                     int input_index, int output_index);

extern void *mode_decision_configuration_kernel(void *input_ptr);
#ifdef __cplusplus
}
#endif
#endif // EbModeDecisionConfigurationProcess_h
