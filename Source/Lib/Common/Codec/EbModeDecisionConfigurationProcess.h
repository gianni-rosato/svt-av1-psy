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
#include "EbModeDecision.h"
#if ADD_MDC_FULL_COST
#include "EbTransQuantBuffers.h"
#endif
#include "EbObject.h"
#ifdef __cplusplus
extern "C" {
#endif

    /**************************************
     * Defines
     **************************************/
    static const uint8_t depth_offset[4] = { 85,21,5,1 };
    static const uint32_t ns_blk_offset[10] = { 0, 1, 3, 25, 5, 8, 11, 14 ,17, 21 };
    static const uint32_t ns_blk_num[10] = { 1, 2, 2, 4, 3, 3, 3, 3, 4, 4 };

    typedef struct MdcpLocalCodingUnit
    {
        uint64_t early_cost;
        EbBool   early_split_flag;
        uint32_t split_context;
        EbBool   selected_cu;
        EbBool   stop_split;
#if PREDICT_NSQ_SHAPE
        PartitionType part;
        uint32_t best_d1_blk;
#endif
    } MdcpLocalCodingUnit;

    typedef struct ModeDecisionConfigurationContext
    {
        EbDctor                              dctor;
        EbFifo                              *rate_control_input_fifo_ptr;
        EbFifo                              *mode_decision_configuration_output_fifo_ptr;

        MdRateEstimationContext           *md_rate_estimation_ptr;
        EbBool                             is_md_rate_estimation_ptr_owner;

        uint8_t                              qp;
        uint64_t                             lambda;
#if PREDICT_NSQ_SHAPE
        MdcpLocalCodingUnit                   local_cu_array[BLOCK_MAX_COUNT_SB_128];
#else
        MdcpLocalCodingUnit                  local_cu_array[CU_MAX_COUNT];
#endif

        // Inter depth decision
        uint8_t                              group_of8x8_blocks_count;
        uint8_t                              group_of16x16_blocks_count;
        uint64_t                              inter_complexity_minimum;
        uint64_t                              inter_complexity_maximum;
        uint64_t                              inter_complexity_average;
        uint64_t                              intra_complexity_minimum;
        uint64_t                              intra_complexity_maximum;
        uint64_t                              intra_complexity_average;
        int16_t                              min_delta_qp_weight;
        int16_t                              max_delta_qp_weight;
        int8_t                               min_delta_qp[4];
        int8_t                               max_delta_qp[4];

        // Adaptive Depth Partitioning
        uint32_t                            *sb_score_array;
        uint8_t                              cost_depth_mode[SB_PRED_OPEN_LOOP_DEPTH_MODE];
        uint8_t                             *sb_cost_array;
        uint32_t                             predicted_cost;
        uint32_t                             budget;
        int8_t                               score_th[MAX_SUPPORTED_SEGMENTS];
        uint8_t                              interval_cost[MAX_SUPPORTED_SEGMENTS];
        uint8_t                              number_of_segments;
        uint32_t                             sb_min_score;
        uint32_t                             sb_max_score;
        uint32_t                             sb_average_score;

        const BlockGeom                     *blk_geom;
        ModeDecisionCandidate             *mdc_candidate_ptr;
        CandidateMv                         *mdc_ref_mv_stack;
        CodingUnit                        *mdc_cu_ptr;
        uint8_t                              qp_index;

        // Multi - Mode signal(s)
        uint8_t                             adp_level;
#if PREDICT_NSQ_SHAPE
        uint32_t                            mds_idx;
#endif
#if ADD_MDC_FULL_COST
        MdcpLocalCodingUnit                  *cu_ptr;
        ModeDecisionCandidateBuffer          *candidate_buffer;
        uint16_t                             sb_origin_y;
        uint16_t                             sb_origin_x;
        uint16_t                             cu_origin_y;
        uint16_t                             cu_origin_x;
        uint16_t                             round_origin_x;
        uint16_t                             round_origin_y;
        uint8_t                              cu_size_log2;
        uint64_t                             three_quad_energy;
        // Transform and Quantization Buffers
        EbTransQuantBuffers                  *trans_quant_buffers_ptr;
        // Trasform Scratch Memory
        int16_t                              *transform_inner_array_ptr;
        uint8_t                               spatial_sse_full_loop;
        EntropyCoder                         *coeff_est_entropy_coder_ptr;
        uint64_t                             full_lambda;
        ModeDecisionCandidate                **fast_candidate_ptr_array;
        ModeDecisionCandidate                *fast_candidate_array;
        ModeDecisionCandidateBuffer          **candidate_buffer_ptr_array;
#endif
    } ModeDecisionConfigurationContext;

    /**************************************
     * Extern Function Declarations
     **************************************/
    extern EbErrorType mode_decision_configuration_context_ctor(
        ModeDecisionConfigurationContext  *context_ptr,
        EbFifo                            *rate_control_input_fifo_ptr,
        EbFifo                            *mode_decision_configuration_output_fifo_ptr,
        uint16_t                           sb_total_count);

    extern void* mode_decision_configuration_kernel(void *input_ptr);
#ifdef __cplusplus
}
#endif
#endif // EbModeDecisionConfigurationProcess_h
