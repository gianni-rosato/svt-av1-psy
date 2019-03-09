/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbModeDecisionConfigurationProcess_h
#define EbModeDecisionConfigurationProcess_h

#include "EbSystemResourceManager.h"
#include "EbMdRateEstimation.h"
#include "EbDefinitions.h"
#include "EbRateControlProcess.h"
#include "EbSequenceControlSet.h"
#if MDC_FIX_0
#include "EbModeDecision.h"
#endif
#ifdef __cplusplus
extern "C" {
#endif

    /**************************************
     * Defines
     **************************************/
    static const uint8_t DepthOffset[4] = { 85,21,5,1 };
    static const uint32_t ns_blk_offset[10] = { 0, 1, 3, 25, 5, 8, 11, 14 ,17, 21 };
    static const uint32_t ns_blk_num[10] = { 1, 2, 2, 4, 3, 3, 3, 3, 4, 4 };

    typedef struct MdcpLocalCodingUnit_s
    {
        uint64_t                          earlyCost;
        EbBool                         earlySplitFlag;
        uint32_t                          splitContext;
        EbBool                         slectedCu;
        EbBool                         stopSplit;
    } MdcpLocalCodingUnit_t;

    typedef struct ModeDecisionConfigurationContext_s
    {
        EbFifo_t                            *rateControlInputFifoPtr;
        EbFifo_t                            *modeDecisionConfigurationOutputFifoPtr;

        MdRateEstimationContext_t           *md_rate_estimation_ptr;

        uint8_t                               qp;
        uint64_t                              lambda;
        MdcpLocalCodingUnit_t               localCuArray[CU_MAX_COUNT];

        // Inter depth decision
        uint8_t                               group_of8x8_blocks_count;
        uint8_t                               group_of16x16_blocks_count;
        uint64_t                              interComplexityMinimum;
        uint64_t                              interComplexityMaximum;
        uint64_t                              interComplexityAverage;
        uint64_t                              intraComplexityMinimum;
        uint64_t                              intraComplexityMaximum;
        uint64_t                              intraComplexityAverage;
        int16_t                              min_delta_qp_weight;
        int16_t                              max_delta_qp_weight;
        int8_t                               min_delta_qp[4];
        int8_t                               max_delta_qp[4];


#if ADAPTIVE_DEPTH_PARTITIONING
        // Adaptive Depth Partitioning
        uint32_t                             *sb_score_array;
        uint8_t                               cost_depth_mode[SB_PRED_OPEN_LOOP_1_NFL_DEPTH_MODE];
        uint8_t                              *sb_cost_array;
        uint32_t                              predicted_cost;
        uint32_t                              budget;
        int8_t                                score_th[MAX_SUPPORTED_SEGMENTS];
        uint8_t                               interval_cost[MAX_SUPPORTED_SEGMENTS];
        uint8_t                               number_of_segments;
        uint32_t                              sb_min_score;
        uint32_t                              sb_max_score;
#else
        // Budgeting
        uint32_t                             *lcuScoreArray;

        uint8_t                                costDepthMode[LCU_PRED_OPEN_LOOP_1_NFL_DEPTH_MODE];

        uint8_t                              *lcuCostArray;
        uint32_t                              predictedCost;
        uint32_t                              budget;

        int8_t                               scoreTh[MAX_SUPPORTED_SEGMENTS];
        uint8_t                               intervalCost[MAX_SUPPORTED_SEGMENTS];
        uint8_t                               numberOfSegments;

        uint32_t                              lcuMinScore;
        uint32_t                              lcuMaxScore;
        EbBool                             depthSensitivePictureFlag;
        EbBool                             performRefinement;
#endif
#if MDC_FIX_0
        const BlockGeom                      *blk_geom;
        ModeDecisionCandidate_t              *mdc_candidate_ptr;
        CandidateMv                          *mdc_ref_mv_stack;
        CodingUnit_t                         *mdc_cu_ptr;
#endif
        uint8_t                               qp_index;

#if ADAPTIVE_DEPTH_PARTITIONING
        // Multi - Mode signal(s)
        uint8_t                               adp_level; // Hsan: to use
#endif
    } ModeDecisionConfigurationContext_t;



    /**************************************
     * Extern Function Declarations
     **************************************/
    extern EbErrorType ModeDecisionConfigurationContextCtor(
        ModeDecisionConfigurationContext_t **context_dbl_ptr,
        EbFifo_t                            *rateControlInputFifoPtr,

        EbFifo_t                            *modeDecisionConfigurationOutputFifoPtr,
        uint16_t                                 sb_total_count);


    extern void* ModeDecisionConfigurationKernel(void *input_ptr);
#ifdef __cplusplus
}
#endif
#endif // EbModeDecisionConfigurationProcess_h