/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbAdaptiveMotionVectorPrediction_h
#define EbAdaptiveMotionVectorPrediction_h

#include "EbUtility.h"
#include "EbPictureControlSet.h"
#include "EbCodingUnit.h"
#include "EbPredictionUnit.h"
#include "EbNeighborArrays.h"
#include "EbMvMerge.h"
#include "EbWarpedMotion.h"

#ifdef __cplusplus
extern "C" {
#endif

    struct ModeDecisionContext_s;
    struct InterPredictionContext_s;

    typedef enum TmvpPos {
        TmvpColocatedBottomRight = 0,
        TmvpColocatedCenter = 1
    } TmvpPos;

    // TMVP items corresponding to one LCU
    typedef struct TmvpUnit_s {
        Mv_t              mv[MAX_NUM_OF_REF_PIC_LIST][MAX_TMVP_CAND_PER_LCU];
        uint64_t            refPicPOC[MAX_NUM_OF_REF_PIC_LIST][MAX_TMVP_CAND_PER_LCU];
        EbPredDirection  prediction_direction[MAX_TMVP_CAND_PER_LCU];
        EbBool              availabilityFlag[MAX_TMVP_CAND_PER_LCU];

        //*Note- list 1 motion info will be added when B-slices are ready

    } TmvpUnit_t;

    extern EbErrorType clip_mv(
        uint32_t                   cu_origin_x,
        uint32_t                   cu_origin_y,
        int16_t                  *MVx,
        int16_t                  *MVy,
        uint32_t                   picture_width,
        uint32_t                   picture_height,
        uint32_t                   tbSize);

    void generate_av1_mvp_table(
        struct ModeDecisionContext_s            *context_ptr,
        CodingUnit_t                     *cu_ptr,
        const BlockGeom                   * blk_geom,
        uint16_t                            cu_origin_x,
        uint16_t                            cu_origin_y,
        MvReferenceFrame                *refFrames,
        uint32_t                            TotRefs,
        PictureControlSet_t              *picture_control_set_ptr);

    void get_av1_mv_pred_drl(
        struct ModeDecisionContext_s            *context_ptr,
        CodingUnit_t      *cu_ptr,
        MvReferenceFrame ref_frame,
        uint8_t              is_compound,
        PredictionMode    mode,
        uint8_t              drl_index,    //valid value of drl_index
        IntMv             nearestmv[2],
        IntMv             nearmv[2],
        IntMv             ref_mv[2]);

    void enc_pass_av1_mv_pred(
        struct ModeDecisionContext_s            *md_context_ptr,
        CodingUnit_t                     *cu_ptr,
        const BlockGeom                   * blk_geom,
        uint16_t                            cu_origin_x,
        uint16_t                            cu_origin_y,
        PictureControlSet_t              *picture_control_set_ptr,
        MvReferenceFrame                ref_frame,
        uint8_t                             is_compound,
        PredictionMode                   mode,
        IntMv                            ref_mv[2] //[OUT]
    );

    void update_mi_map(
        CodingUnit_t                   *cu_ptr,
        uint32_t                          cu_origin_x,
        uint32_t                          cu_origin_y,
        const BlockGeom               * blk_geom,
        const CodedUnitStats_t         *cu_stats,
        PictureControlSet_t            *picture_control_set_ptr);

    uint16_t wm_find_samples(
        CodingUnit_t                       *cu_ptr,
        const BlockGeom                    *blk_geom,
        uint16_t                            pu_origin_x,
        uint16_t                            pu_origin_y,
        PictureControlSet_t                *picture_control_set_ptr,
        int32_t                            *pts,
        int32_t                            *pts_inref);

    EbBool warped_motion_parameters(
        PictureControlSet_t              *picture_control_set_ptr,
        CodingUnit_t                     *cu_ptr,
        MvUnit_t                         *mv_unit,
        const BlockGeom                  *blk_geom,
        uint16_t                          pu_origin_x,
        uint16_t                          pu_origin_y,
        EbWarpedMotionParams             *wm_params,
        uint16_t                         *num_samples);

#ifdef __cplusplus
}
#endif
#endif // EbAdaptiveMotionVectorPrediction_h
