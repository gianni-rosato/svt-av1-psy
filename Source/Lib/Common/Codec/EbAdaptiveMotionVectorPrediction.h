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

    struct ModeDecisionContext;
    struct InterPredictionContext;

    typedef enum TmvpPos
    {
        TmvpColocatedBottomRight = 0,
        TmvpColocatedCenter = 1
    } TmvpPos;

    // TMVP items corresponding to one LCU
    typedef struct TmvpUnit
    {
        Mv              mv[MAX_NUM_OF_REF_PIC_LIST][MAX_TMVP_CAND_PER_LCU];
        uint64_t            ref_pic_poc[MAX_NUM_OF_REF_PIC_LIST][MAX_TMVP_CAND_PER_LCU];
        EbPredDirection  prediction_direction[MAX_TMVP_CAND_PER_LCU];
        EbBool              availability_flag[MAX_TMVP_CAND_PER_LCU];

        //*Note- list 1 motion info will be added when B-slices are ready
    } TmvpUnit;

    extern EbErrorType clip_mv(
        uint32_t  cu_origin_x,
        uint32_t  cu_origin_y,
        int16_t  *mv_x,
        int16_t  *mv_y,
        uint32_t  picture_width,
        uint32_t  picture_height,
        uint32_t  tb_size);

    void generate_av1_mvp_table(
        TileInfo                   *tile,
        struct ModeDecisionContext *context_ptr,
        CodingUnit                 *cu_ptr,
        const BlockGeom            *blk_geom,
        uint16_t                    cu_origin_x,
        uint16_t                    cu_origin_y,
        MvReferenceFrame           *ref_frames,
        uint32_t                    tot_refs,
        PictureControlSet          *picture_control_set_ptr);

    void get_av1_mv_pred_drl(
        struct ModeDecisionContext *context_ptr,
        CodingUnit                 *cu_ptr,
        MvReferenceFrame            ref_frame,
        uint8_t                     is_compound,
        PredictionMode              mode,
        uint8_t                     drl_index,
        IntMv                       nearestmv[2],
        IntMv                       nearmv[2],
        IntMv                       ref_mv[2]);

    void enc_pass_av1_mv_pred(
        TileInfo                   *tile,
        struct ModeDecisionContext *md_context_ptr,
        CodingUnit                 *cu_ptr,
        const BlockGeom            *blk_geom,
        uint16_t                    cu_origin_x,
        uint16_t                    cu_origin_y,
        PictureControlSet          *picture_control_set_ptr,
        MvReferenceFrame            ref_frame,
        uint8_t                     is_compound,
        PredictionMode              mode,
        IntMv                       ref_mv[2]);

    void update_mi_map(
        struct ModeDecisionContext   *context_ptr,
        CodingUnit                   *cu_ptr,
        uint32_t                          cu_origin_x,
        uint32_t                          cu_origin_y,
        const BlockGeom               * blk_geom,
        const CodedUnitStats         *cu_stats,
        PictureControlSet            *picture_control_set_ptr);

    uint16_t wm_find_samples(
        CodingUnit                       *cu_ptr,
        const BlockGeom                    *blk_geom,
        uint16_t                            cu_origin_x,
        uint16_t                            cu_origin_y,
        MvReferenceFrame                    rf0,
        PictureControlSet                *picture_control_set_ptr,
        int32_t                            *pts,
        int32_t                            *pts_inref);

    void wm_count_samples(
        CodingUnit                       *cu_ptr,
        const BlockGeom                    *blk_geom,
        uint16_t                            cu_origin_x,
        uint16_t                            cu_origin_y,
        uint8_t                             ref_frame_type,
        PictureControlSet                *picture_control_set_ptr,
        uint16_t                           *num_samples);

    EbBool warped_motion_parameters(
        PictureControlSet              *picture_control_set_ptr,
        CodingUnit                     *cu_ptr,
        MvUnit                         *mv_unit,
        const BlockGeom                  *blk_geom,
        uint16_t                          cu_origin_x,
        uint16_t                          cu_origin_y,
        uint8_t                           ref_frame_type,
        EbWarpedMotionParams             *wm_params,
        uint16_t                         *num_samples);

    static INLINE EbBool is_motion_variation_allowed_bsize(const BlockSize bsize)
    {
        return (block_size_wide[bsize] >= 8 && block_size_high[bsize] >= 8);
    }

    static INLINE int is_neighbor_overlappable(const MbModeInfo *mbmi)
    {
        return /*is_intrabc_block(mbmi) ||*/ mbmi->ref_frame[0] > INTRA_FRAME; // TODO: modify when add intra_bc
    }

    static INLINE EbBool has_overlappable_candidates(const CodingUnit *cu_ptr)
    {
        return (cu_ptr->prediction_unit_array[0].overlappable_neighbors[0] != 0
             || cu_ptr->prediction_unit_array[0].overlappable_neighbors[1] != 0);
    }

    void av1_count_overlappable_neighbors(
        const PictureControlSet        *picture_control_set_ptr,
        CodingUnit                     *cu_ptr,
        const BlockSize                   bsize,
        int32_t                           mi_row,
        int32_t                           mi_col);

    void av1_find_best_ref_mvs_from_stack(int allow_hp,
        CandidateMv ref_mv_stack[][MAX_REF_MV_STACK_SIZE],
        MacroBlockD * xd,
        MvReferenceFrame ref_frame,
        IntMv *nearest_mv, IntMv *near_mv,
        int is_integer);
    void av1_find_ref_dv(IntMv *ref_dv, const TileInfo *const tile,
        int mib_size, int mi_row, int mi_col);
    int av1_is_dv_valid(const MV dv,
        const MacroBlockD *xd, int mi_row, int mi_col,
        BlockSize bsize, int mib_size_log2);

#ifdef __cplusplus
}
#endif
#endif // EbAdaptiveMotionVectorPrediction_h
