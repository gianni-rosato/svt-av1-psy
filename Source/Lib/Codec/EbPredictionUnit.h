/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPredictionUnit_h
#define EbPredictionUnit_h

#include "EbDefinitions.h"
#include "EbSyntaxElements.h"
#include "EbDefinitions.h"
#include "EbMotionVectorUnit.h"
#ifdef __cplusplus
extern "C" {
#endif
#pragma pack(push, 1)
    typedef struct PredictionUnit_s
    {
        Mv_t                         mv[MAX_NUM_OF_REF_PIC_LIST];  // 16-bytes
        Mvd_t                        mvd[MAX_NUM_OF_REF_PIC_LIST]; // 16-bytes
#if !INTRA_INTER_FAST_LOOP
        unsigned                     merge_index                : 5;
#endif
        unsigned                     merge_flag                 : 1;
        unsigned                     inter_pred_direction_index : 2;
        unsigned                     intra_luma_mode            : 6;
        unsigned                     intra_luma_left_mode       : 6;
        unsigned                     intra_luma_top_mode        : 6;

        uint8_t                      intra_chroma_left_mode;
        uint8_t                      intra_chroma_top_mode;
        // Intra Mode
        int32_t                      angle_delta[PLANE_TYPES];
        EbBool                       is_directional_mode_flag;
        EbBool                       is_directional_chroma_mode_flag;
        EbBool                       use_angle_delta;
        uint32_t                     intra_chroma_mode;

        // Inter Mode
        PredictionMode               inter_mode;
        EbBool                       is_compound;
        uint32_t                     pred_mv_weight;
        uint8_t                      ref_frame_type;
        uint8_t                      ref_mv_index;
#if !INTRA_INTER_FAST_LOOP
        EbBool                       is_skip_mode_flag;
#endif
        EbBool                       is_new_mv;
        EbBool                       is_zero_mv;

        MOTION_MODE                   motion_mode;
        uint16_t                      num_proj_ref;
        EbWarpedMotionParams          wm_params;
        uint32_t                      overlappable_neighbors[2];

        // Index of the alpha Cb and alpha Cr combination
        int32_t                      cfl_alpha_idx;
        // Joint sign of alpha Cb and alpha Cr
        int32_t                      cfl_alpha_signs;

    } PredictionUnit_t;
#pragma pack(pop)


#ifdef __cplusplus
}
#endif
#endif //EbPredictionUnit_h