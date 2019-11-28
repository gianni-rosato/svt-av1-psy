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

#define MV_BORDER (16 << 3)  // Allow 16 pels in 1/8th pel units

#define INTRABC_DELAY_PIXELS 256  //  Delay of 256 pixels
#define INTRABC_DELAY_SB64  (INTRABC_DELAY_PIXELS / 64)

#define NELEMENTS(x) (int)(sizeof(x) / sizeof(x[0]))

// Although we assign 32 bit integers, all the values are strictly under 14
// bits.
    static int div_mult[32] = { 0,    16384, 8192, 5461, 4096, 3276, 2730, 2340,
                                2048, 1820,  1638, 1489, 1365, 1260, 1170, 1092,
                                1024, 963,   910,  862,  819,  780,  744,  712,
                                682,  655,   630,  606,  585,  564,  546,  528 };

    struct ModeDecisionContext;
    struct InterPredictionContext;

    extern EbErrorType clip_mv(
        uint32_t  cu_origin_x,
        uint32_t  cu_origin_y,
        int16_t  *mv_x,
        int16_t  *mv_y,
        uint32_t  picture_width,
        uint32_t  picture_height,
        uint32_t  tb_size);

#if MULTI_PASS_PD
    void mvp_bypass_init(
        PictureControlSet          *picture_control_set_ptr,
        struct ModeDecisionContext *context_ptr);
#endif

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
        return /*is_intrabc_block(mbmi) ||*/ mbmi->block_mi.ref_frame[0] > INTRA_FRAME; // TODO: modify when add intra_bc
    }

    static INLINE EbBool has_overlappable_candidates(const CodingUnit *cu_ptr)
    {
        return (cu_ptr->prediction_unit_array[0].overlappable_neighbors[0] != 0
             || cu_ptr->prediction_unit_array[0].overlappable_neighbors[1] != 0);
    }

    static INLINE void integer_mv_precision(MV *mv) {
        int mod = (mv->row % 8);
        if (mod != 0) {
            mv->row -= mod;
            if (abs(mod) > 4) {
                if (mod > 0)
                    mv->row += 8;
                else
                    mv->row -= 8;
            }
        }

        mod = (mv->col % 8);
        if (mod != 0) {
            mv->col -= mod;
            if (abs(mod) > 4) {
                if (mod > 0)
                    mv->col += 8;
                else
                    mv->col -= 8;
            }
        }
    }

    static INLINE void lower_mv_precision(MV *mv, int allow_hp, int is_integer) {
        if (is_integer)
            integer_mv_precision(mv);
        else {
            if (!allow_hp) {
                if (mv->row & 1) mv->row += (mv->row > 0 ? -1 : 1);
                if (mv->col & 1) mv->col += (mv->col > 0 ? -1 : 1);
            }
        }
    }

    static INLINE void get_mv_projection(MV *output, MV ref, int num, int den) {
        den = AOMMIN(den, MAX_FRAME_DISTANCE);
        num = num > 0 ? AOMMIN(num, MAX_FRAME_DISTANCE)
            : AOMMAX(num, -MAX_FRAME_DISTANCE);
        const int mv_row =
            ROUND_POWER_OF_TWO_SIGNED(ref.row * num * div_mult[den], 14);
        const int mv_col =
            ROUND_POWER_OF_TWO_SIGNED(ref.col * num * div_mult[den], 14);
        const int clamp_max = MV_UPP - 1;
        const int clamp_min = MV_LOW + 1;
        output->row = (int16_t)clamp(mv_row, clamp_min, clamp_max);
        output->col = (int16_t)clamp(mv_col, clamp_min, clamp_max);
    }

    static INLINE PredictionMode compound_ref0_mode(PredictionMode mode) {
        static PredictionMode lut[] = {
          MB_MODE_COUNT,  // DC_PRED
          MB_MODE_COUNT,  // V_PRED
          MB_MODE_COUNT,  // H_PRED
          MB_MODE_COUNT,  // D45_PRED
          MB_MODE_COUNT,  // D135_PRED
          MB_MODE_COUNT,  // D113_PRED
          MB_MODE_COUNT,  // D157_PRED
          MB_MODE_COUNT,  // D203_PRED
          MB_MODE_COUNT,  // D67_PRED
          MB_MODE_COUNT,  // SMOOTH_PRED
          MB_MODE_COUNT,  // SMOOTH_V_PRED
          MB_MODE_COUNT,  // SMOOTH_H_PRED
          MB_MODE_COUNT,  // PAETH_PRED
          MB_MODE_COUNT,  // NEARESTMV
          MB_MODE_COUNT,  // NEARMV
          MB_MODE_COUNT,  // GLOBALMV
          MB_MODE_COUNT,  // NEWMV
          NEARESTMV,      // NEAREST_NEARESTMV
          NEARMV,         // NEAR_NEARMV
          NEARESTMV,      // NEAREST_NEWMV
          NEWMV,          // NEW_NEARESTMV
          NEARMV,         // NEAR_NEWMV
          NEWMV,          // NEW_NEARMV
          GLOBALMV,       // GLOBAL_GLOBALMV
          NEWMV,          // NEW_NEWMV
        };
        assert(NELEMENTS(lut) == MB_MODE_COUNT);
        assert(is_inter_compound_mode(mode));
        return lut[mode];
    }

    static INLINE PredictionMode compound_ref1_mode(PredictionMode mode) {
        static PredictionMode lut[] = {
          MB_MODE_COUNT,  // DC_PRED
          MB_MODE_COUNT,  // V_PRED
          MB_MODE_COUNT,  // H_PRED
          MB_MODE_COUNT,  // D45_PRED
          MB_MODE_COUNT,  // D135_PRED
          MB_MODE_COUNT,  // D113_PRED
          MB_MODE_COUNT,  // D157_PRED
          MB_MODE_COUNT,  // D203_PRED
          MB_MODE_COUNT,  // D67_PRED
          MB_MODE_COUNT,  // SMOOTH_PRED
          MB_MODE_COUNT,  // SMOOTH_V_PRED
          MB_MODE_COUNT,  // SMOOTH_H_PRED
          MB_MODE_COUNT,  // PAETH_PRED
          MB_MODE_COUNT,  // NEARESTMV
          MB_MODE_COUNT,  // NEARMV
          MB_MODE_COUNT,  // GLOBALMV
          MB_MODE_COUNT,  // NEWMV
          NEARESTMV,      // NEAREST_NEARESTMV
          NEARMV,         // NEAR_NEARMV
          NEWMV,          // NEAREST_NEWMV
          NEARESTMV,      // NEW_NEARESTMV
          NEWMV,          // NEAR_NEWMV
          NEARMV,         // NEW_NEARMV
          GLOBALMV,       // GLOBAL_GLOBALMV
          NEWMV,          // NEW_NEWMV
        };
        assert(NELEMENTS(lut) == MB_MODE_COUNT);
        assert(is_inter_compound_mode(mode));
        return lut[mode];
    }

    static INLINE int check_sb_border(const int mi_row, const int mi_col,
        const int row_offset, const int col_offset)
    {
        const int sb_mi_size = mi_size_wide[BLOCK_64X64];
        const int row = mi_row & (sb_mi_size - 1);
        const int col = mi_col & (sb_mi_size - 1);

        if (row + row_offset < 0 || row + row_offset >= sb_mi_size ||
            col + col_offset < 0 || col + col_offset >= sb_mi_size)
            return 0;

        return 1;
    }

    static INLINE int is_global_mv_block(
        const PredictionMode          mode,
        const BlockSize               bsize,
        TransformationType            type)
    {
        return (mode == GLOBALMV || mode == GLOBAL_GLOBALMV)
            && type > TRANSLATION
            && is_motion_variation_allowed_bsize(bsize);
    }

    static INLINE void clamp_mv(MV *mv, int32_t min_col, int32_t max_col,
        int32_t min_row, int32_t max_row)
    {
        mv->col = (int16_t)clamp(mv->col, min_col, max_col);
        mv->row = (int16_t)clamp(mv->row, min_row, max_row);
    }

    static INLINE int32_t is_mv_valid(const MV *mv) {
        return mv->row > MV_LOW && mv->row < MV_UPP && mv->col > MV_LOW &&
            mv->col < MV_UPP;
    }

    void eb_av1_count_overlappable_neighbors(
        const PictureControlSet        *picture_control_set_ptr,
        CodingUnit                     *cu_ptr,
        const BlockSize                   bsize,
        int32_t                           mi_row,
        int32_t                           mi_col);

    void eb_av1_find_best_ref_mvs_from_stack(int allow_hp,
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
    int is_inside_tile_boundary(TileInfo *tile, int16_t mvx, int16_t mvy, int mi_col, int mi_row, BlockSize bsize);

    IntMv gm_get_motion_vector_enc(
        const EbWarpedMotionParams *gm,
        int32_t allow_hp,
        BlockSize bsize,
        int32_t mi_col, int32_t mi_row,
        int32_t is_integer);

#ifdef __cplusplus
}
#endif
#endif // EbAdaptiveMotionVectorPrediction_h
