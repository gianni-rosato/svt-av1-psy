/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

/*
* Copyright (c) 2016, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at www.aomedia.org/license/software. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at www.aomedia.org/license/patent.
*/

#include <string.h>

#include "EbDefinitions.h"
#include "EbUtility.h"
#include "EbAdaptiveMotionVectorPrediction.h"
#include "EbSequenceControlSet.h"
#include "EbReferenceObject.h"
#include "EbErrorCodes.h"
#include "EbModeDecisionProcess.h"

#define UNUSED_FUNC
static PartitionType from_shape_to_part[] = {
    PARTITION_NONE,
    PARTITION_HORZ,
    PARTITION_VERT,
    PARTITION_HORZ_A,
    PARTITION_HORZ_B,
    PARTITION_VERT_A,
    PARTITION_VERT_B,
    PARTITION_HORZ_4,
    PARTITION_VERT_4,
    PARTITION_SPLIT
};

/** ScaleMV
        is used to scale the motion vector in AMVP process.
 */
/*
static inline void scale_mv(
    uint64_t    current_pic_poc,                // Iuput parameter, the POC of the current picture to be encoded.
    uint64_t    target_ref_pic_poc,              // Iuput parameter, the POC of the reference picture where the inter coding is searching for.
    uint64_t    col_pu_pic_poc,                  // Iuput parameter, the POC of picture where the co-located PU is.
    uint64_t    col_pu_ref_pic_poc,               // Iuput parameter, the POC of the reference picture where the MV of the co-located PU points to.
    int16_t    *mvx,                          // Output parameter,
    int16_t    *mvy)                          // Output parameter,
{
    int16_t td = (int16_t)(col_pu_pic_poc - col_pu_ref_pic_poc);
    int16_t tb = (int16_t)(current_pic_poc - target_ref_pic_poc);
    int16_t scale_factor;
    int16_t temp;

    if (td != tb) {
        tb = CLIP3(-128, 127, tb);
        td = CLIP3(-128, 127, td);
        temp = (int16_t)((0x4000 + ABS(td >> 1)) / td);
        scale_factor = CLIP3(-4096, 4095, (tb * temp + 32) >> 6);

        *mvx = CLIP3(-32768, 32767, (scale_factor * (*mvx) + 127 + (scale_factor * (*mvx) < 0)) >> 8);
        *mvy = CLIP3(-32768, 32767, (scale_factor * (*mvy) + 127 + (scale_factor * (*mvy) < 0)) >> 8);
    }

    return;
}
*/

EbErrorType clip_mv(
    uint32_t                   cu_origin_x,
    uint32_t                   cu_origin_y,
    int16_t                    *mvx,
    int16_t                    *mvy,
    uint32_t                   picture_width,
    uint32_t                   picture_height,
    uint32_t                   tb_size)
{
    EbErrorType return_error = EB_ErrorNone;

    // horizontal clipping
    (*mvx) = CLIP3(((int16_t)((1 - cu_origin_x - 8 - tb_size) << 2)), ((int16_t)((picture_width + 8 - cu_origin_x - 1) << 2)), (*mvx));
    // vertical clipping
    (*mvy) = CLIP3(((int16_t)((1 - cu_origin_y - 8 - tb_size) << 2)), ((int16_t)((picture_height + 8 - cu_origin_y - 1) << 2)), (*mvy));
#if AV1_UPGRADE
    const int32_t clamp_max = MV_UPP - 1;
    const int32_t clamp_min = MV_LOW + 1;
    // horizontal clipping
    (*mvx) = CLIP3(clamp_min, clamp_max, (*mvx));
    // vertical clipping
    (*mvy) = CLIP3(clamp_min, clamp_max, (*mvy));
#endif
    return return_error;
}

#define MVREF_ROWS 3
#define MVREF_COLS 3

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
    return lut[mode];
}


/*static INLINE*/ int32_t is_inter_block(const MbModeInfo *mbmi) {
    return /*is_intrabc_block(mbmi) ||*/ mbmi->ref_frame[0] > INTRA_FRAME;
}

static int32_t have_newmv_in_inter_mode(PredictionMode mode) {
    return (mode == NEWMV || mode == NEW_NEWMV || mode == NEAREST_NEWMV ||
        mode == NEW_NEARESTMV || mode == NEAR_NEWMV || mode == NEW_NEARMV);
}
#define NELEMENTS(x) (int32_t)(sizeof(x) / sizeof(x[0]))
MvReferenceFrame comp_ref0(int32_t ref_idx) {
    static const MvReferenceFrame lut[] = {
        LAST_FRAME,     // LAST_LAST2_FRAMES,
        LAST_FRAME,     // LAST_LAST3_FRAMES,
        LAST_FRAME,     // LAST_GOLDEN_FRAMES,
        BWDREF_FRAME,   // BWDREF_ALTREF_FRAMES,
        LAST2_FRAME,    // LAST2_LAST3_FRAMES
        LAST2_FRAME,    // LAST2_GOLDEN_FRAMES,
        LAST3_FRAME,    // LAST3_GOLDEN_FRAMES,
        BWDREF_FRAME,   // BWDREF_ALTREF2_FRAMES,
        ALTREF2_FRAME,  // ALTREF2_ALTREF_FRAMES,
    };
    assert(NELEMENTS(lut) == TOTAL_UNIDIR_COMP_REFS);
    return lut[ref_idx];
}

MvReferenceFrame comp_ref1(int32_t ref_idx) {
    static const MvReferenceFrame lut[] = {
        LAST2_FRAME,    // LAST_LAST2_FRAMES,
        LAST3_FRAME,    // LAST_LAST3_FRAMES,
        GOLDEN_FRAME,   // LAST_GOLDEN_FRAMES,
        ALTREF_FRAME,   // BWDREF_ALTREF_FRAMES,
        LAST3_FRAME,    // LAST2_LAST3_FRAMES
        GOLDEN_FRAME,   // LAST2_GOLDEN_FRAMES,
        GOLDEN_FRAME,   // LAST3_GOLDEN_FRAMES,
        ALTREF2_FRAME,  // BWDREF_ALTREF2_FRAMES,
        ALTREF_FRAME,   // ALTREF2_ALTREF_FRAMES,
    };
    assert(NELEMENTS(lut) == TOTAL_UNIDIR_COMP_REFS);
    return lut[ref_idx];
}

typedef struct position {
    int32_t row;
    int32_t col;
} Position;

#define MV_BORDER (16 << 3)  // Allow 16 pels in 1/8th pel units

static MvReferenceFrame ref_frame_map[TOTAL_COMP_REFS][2] = {
    { LAST_FRAME, BWDREF_FRAME },{ LAST2_FRAME, BWDREF_FRAME },
    { LAST3_FRAME, BWDREF_FRAME },{ GOLDEN_FRAME, BWDREF_FRAME },

    { LAST_FRAME, ALTREF2_FRAME },{ LAST2_FRAME, ALTREF2_FRAME },
    { LAST3_FRAME, ALTREF2_FRAME },{ GOLDEN_FRAME, ALTREF2_FRAME },

    { LAST_FRAME, ALTREF_FRAME },{ LAST2_FRAME, ALTREF_FRAME },
    { LAST3_FRAME, ALTREF_FRAME },{ GOLDEN_FRAME, ALTREF_FRAME },

    { LAST_FRAME, LAST2_FRAME },{ LAST_FRAME, LAST3_FRAME },
    { LAST_FRAME, GOLDEN_FRAME },{ BWDREF_FRAME, ALTREF_FRAME },

    // NOTE: Following reference frame pairs are not supported to be explicitly
    //       signalled, but they are possibly chosen by the use of skip_mode,
    //       which may use the most recent one-sided reference frame pair.
    { LAST2_FRAME, LAST3_FRAME },{ LAST2_FRAME, GOLDEN_FRAME },
    { LAST3_FRAME, GOLDEN_FRAME },{ BWDREF_FRAME, ALTREF2_FRAME },
    { ALTREF2_FRAME, ALTREF_FRAME }
};

static void clamp_mv(
    MV *mv,
    int32_t min_col,
    int32_t max_col,
    int32_t min_row,
    int32_t max_row) {
    mv->col = (int16_t)clamp(mv->col, min_col, max_col);
    mv->row = (int16_t)clamp(mv->row, min_row, max_row);
}

// clang-format on

void av1_set_ref_frame(MvReferenceFrame *rf,
    int8_t ref_frame_type) {
    if (ref_frame_type >= TOTAL_REFS_PER_FRAME) {
        rf[0] = ref_frame_map[ref_frame_type - TOTAL_REFS_PER_FRAME][0];
        rf[1] = ref_frame_map[ref_frame_type - TOTAL_REFS_PER_FRAME][1];
    }
    else {
        rf[0] = ref_frame_type;
        rf[1] = NONE_FRAME;
        // assert(ref_frame_type > NONE_FRAME); AMIR
    }
}
int8_t get_uni_comp_ref_idx(const MvReferenceFrame *const rf) {
    // Single ref pred
    if (rf[1] <= INTRA_FRAME) return -1;

    // Bi-directional comp ref pred
    if ((rf[0] < BWDREF_FRAME) && (rf[1] >= BWDREF_FRAME)) return -1;

    for (int8_t ref_idx = 0; ref_idx < TOTAL_UNIDIR_COMP_REFS; ++ref_idx) {
        if (rf[0] == comp_ref0(ref_idx) && rf[1] == comp_ref1(ref_idx))
            return ref_idx;
    }
    return -1;
}

extern INLINE int8_t av1_ref_frame_type(const MvReferenceFrame *const rf) {
    if (rf[1] > INTRA_FRAME) {
        const int8_t uni_comp_ref_idx = get_uni_comp_ref_idx(rf);
        if (uni_comp_ref_idx >= 0) {
            assert((TOTAL_REFS_PER_FRAME + FWD_REFS * BWD_REFS + uni_comp_ref_idx) <
                MODE_CTX_REF_FRAMES);
            return TOTAL_REFS_PER_FRAME + FWD_REFS * BWD_REFS + uni_comp_ref_idx;
        }
        else {
            return TOTAL_REFS_PER_FRAME + FWD_RF_OFFSET(rf[0]) +
                BWD_RF_OFFSET(rf[1]) * FWD_REFS;
        }
    }

    return rf[0];
}

static INLINE IntMv get_sub_block_mv(const ModeInfo *candidate, int32_t which_mv,
    int32_t search_col) {

    (void)search_col;
    return candidate->mbmi.mv[which_mv];
}
static INLINE int32_t is_inside(const TileInfo *const tile, int32_t mi_col, int32_t mi_row,
    int32_t mi_rows, const Position *mi_pos) {

    const int32_t dependent_horz_tile_flag = 0;
    if (dependent_horz_tile_flag && !tile->tg_horz_boundary) {
        return !(mi_row + mi_pos->row < 0 ||
            mi_col + mi_pos->col < tile->mi_col_start ||
            mi_row + mi_pos->row >= mi_rows ||
            mi_col + mi_pos->col >= tile->mi_col_end);
    }
    else {
        return !(mi_row + mi_pos->row < tile->mi_row_start ||
            mi_col + mi_pos->col < tile->mi_col_start ||
            mi_row + mi_pos->row >= tile->mi_row_end ||
            mi_col + mi_pos->col >= tile->mi_col_end);
    }
}

static INLINE void clamp_mv_ref(MV *mv, int32_t bw, int32_t bh, const MacroBlockD *xd) {
    clamp_mv(mv, xd->mb_to_left_edge - bw * 8 - MV_BORDER,
        xd->mb_to_right_edge + bw * 8 + MV_BORDER,
        xd->mb_to_top_edge - bh * 8 - MV_BORDER,
        xd->mb_to_bottom_edge + bh * 8 + MV_BORDER);
}

static void add_ref_mv_candidate(
    const ModeInfo *const candidate_mi, const MbModeInfo *const candidate,
    const MvReferenceFrame rf[2], uint8_t refmv_counts[MODE_CTX_REF_FRAMES],
    uint8_t ref_match_counts[MODE_CTX_REF_FRAMES],
    uint8_t newmv_counts[MODE_CTX_REF_FRAMES],
    CandidateMv ref_mv_stacks[][MAX_REF_MV_STACK_SIZE], int32_t len,
#if USE_CUR_GM_REFMV
    IntMv *gm_mv_candidates, const EbWarpedMotionParams *gm_params,
#endif  // USE_CUR_GM_REFMV
    int32_t col, int32_t weight)
{

    if (!is_inter_block(candidate)) return;  // for intrabc
    int32_t index = 0, ref;
    assert(weight % 2 == 0);

    if (rf[1] == NONE_FRAME) {
        uint8_t *refmv_count = &refmv_counts[rf[0]];
        uint8_t *ref_match_count = &ref_match_counts[rf[0]];
        uint8_t *newmv_count = &newmv_counts[rf[0]];
        CandidateMv *ref_mv_stack = ref_mv_stacks[rf[0]];
        (void)ref_match_count;

        // single reference frame
        for (ref = 0; ref < 2; ++ref) {
            if (candidate->ref_frame[ref] == rf[0]) {
                IntMv this_refmv;
#if USE_CUR_GM_REFMV    //CHKN this should not be used for TRANSLATION ME model
                if (is_global_mv_block(candidate_mi, gm_params[rf[0]].wmtype))
                    this_refmv = gm_mv_candidates[0];
                else
#endif  // USE_CUR_GM_REFMV
                    this_refmv = get_sub_block_mv(candidate_mi, ref, col);

                for (index = 0; index < *refmv_count; ++index)
                    if (ref_mv_stack[index].this_mv.as_int == this_refmv.as_int) break;

                if (index < *refmv_count) ref_mv_stack[index].weight += weight * len;

                // Add a new item to the list.
                if (index == *refmv_count && *refmv_count < MAX_REF_MV_STACK_SIZE) {
                    ref_mv_stack[index].this_mv = this_refmv;
                    ref_mv_stack[index].weight = weight * len;
                    ++(*refmv_count);
                }
                if (have_newmv_in_inter_mode(candidate->mode))++*newmv_count;
                ++*ref_match_count;
            }
        }
    }
    else {
        MvReferenceFrame ref_frame = av1_ref_frame_type(rf);
        uint8_t *refmv_count = &refmv_counts[ref_frame];
        uint8_t *ref_match_count = &ref_match_counts[ref_frame];
        uint8_t *newmv_count = &newmv_counts[ref_frame];
        CandidateMv *ref_mv_stack = ref_mv_stacks[ref_frame];
        (void)ref_match_count;

        // compound reference frame
        if (candidate->ref_frame[0] == rf[0] && candidate->ref_frame[1] == rf[1]) {
            IntMv this_refmv[2];

            for (ref = 0; ref < 2; ++ref) {
#if USE_CUR_GM_REFMV
                if (is_global_mv_block(candidate_mi, gm_params[rf[ref]].wmtype))
                    this_refmv[ref] = gm_mv_candidates[ref];
                else
#endif  // USE_CUR_GM_REFMV
                    this_refmv[ref] = get_sub_block_mv(candidate_mi, ref, col);
        }

            for (index = 0; index < *refmv_count; ++index)
                if ((ref_mv_stack[index].this_mv.as_int == this_refmv[0].as_int) &&
                    (ref_mv_stack[index].comp_mv.as_int == this_refmv[1].as_int))
                    break;

            if (index < *refmv_count) ref_mv_stack[index].weight += weight * len;

            // Add a new item to the list.
            if (index == *refmv_count && *refmv_count < MAX_REF_MV_STACK_SIZE) {
                ref_mv_stack[index].this_mv = this_refmv[0];
                ref_mv_stack[index].comp_mv = this_refmv[1];
                ref_mv_stack[index].weight = weight * len;
                ++(*refmv_count);
            }
            if (have_newmv_in_inter_mode(candidate->mode))++*newmv_count;
            ++*ref_match_count;
    }
}
}
static void scan_row_mbmi(const Av1Common *cm, const MacroBlockD *xd,
    int32_t mi_row, int32_t mi_col,
    const MvReferenceFrame rf[2], int32_t row_offset,
    CandidateMv ref_mv_stack[][MAX_REF_MV_STACK_SIZE],
    uint8_t refmv_count[MODE_CTX_REF_FRAMES],
    uint8_t ref_match_count[MODE_CTX_REF_FRAMES],
    uint8_t newmv_count[MODE_CTX_REF_FRAMES],
#if USE_CUR_GM_REFMV
    IntMv *gm_mv_candidates,
#endif  // USE_CUR_GM_REFMV
    int32_t max_row_offset, int32_t *processed_rows) {
    int32_t end_mi = AOMMIN(xd->n8_w, cm->mi_cols - mi_col);
    end_mi = AOMMIN(end_mi, mi_size_wide[BLOCK_64X64]);
    const int32_t n8_w_8 = mi_size_wide[BLOCK_8X8];
    const int32_t n8_w_16 = mi_size_wide[BLOCK_16X16];
    int32_t i;
    int32_t col_offset = 0;
    const int32_t shift = 0;
    // TODO(jingning): Revisit this part after cb4x4 is stable.
    if (abs(row_offset) > 1) {
        col_offset = 1;
        if (mi_col & 0x01 && xd->n8_w < n8_w_8) --col_offset;
    }
    const int32_t use_step_16 = (xd->n8_w >= 16);
    ModeInfo **const candidate_mi0 = xd->mi + row_offset * xd->mi_stride;
    (void)mi_row;

    for (i = 0; i < end_mi;) {
        const ModeInfo *const candidate_mi = candidate_mi0[col_offset + i];
        const MbModeInfo *const candidate = &candidate_mi->mbmi;
        const int32_t candidate_bsize = candidate->sb_type;
        ASSERT(candidate_bsize < BlockSizeS_ALL);
        const int32_t n8_w = mi_size_wide[candidate_bsize];
        int32_t len = AOMMIN(xd->n8_w, n8_w);
        if (use_step_16)
            len = AOMMAX(n8_w_16, len);
        else if (abs(row_offset) > 1)
            len = AOMMAX(len, n8_w_8);

        int32_t weight = 2;
        if (xd->n8_w >= n8_w_8 && xd->n8_w <= n8_w) {
            int32_t inc = AOMMIN(-max_row_offset + row_offset + 1,
                mi_size_high[candidate_bsize]);
            // Obtain range used in weight calculation.
            weight = AOMMAX(weight, (inc << shift));
            // Update processed rows.
            *processed_rows = inc - row_offset - 1;
        }

        add_ref_mv_candidate(candidate_mi, candidate, rf, refmv_count,
            ref_match_count, newmv_count, ref_mv_stack, len,
#if USE_CUR_GM_REFMV
            gm_mv_candidates, cm->global_motion,
#endif  // USE_CUR_GM_REFMV
            col_offset + i, weight);

        i += len;
    }
}

static void scan_col_mbmi(const Av1Common *cm, const MacroBlockD *xd,
    int32_t mi_row, int32_t mi_col,
    const MvReferenceFrame rf[2], int32_t col_offset,
    CandidateMv ref_mv_stack[][MAX_REF_MV_STACK_SIZE],
    uint8_t refmv_count[MODE_CTX_REF_FRAMES],
    uint8_t ref_match_count[MODE_CTX_REF_FRAMES],
    uint8_t newmv_count[MODE_CTX_REF_FRAMES],
#if USE_CUR_GM_REFMV
    IntMv *gm_mv_candidates,
#endif  // USE_CUR_GM_REFMV
    int32_t max_col_offset, int32_t *processed_cols) {
    int32_t end_mi = AOMMIN(xd->n8_h, cm->mi_rows - mi_row);
    end_mi = AOMMIN(end_mi, mi_size_high[BLOCK_64X64]);
    const int32_t n8_h_8 = mi_size_high[BLOCK_8X8];
    const int32_t n8_h_16 = mi_size_high[BLOCK_16X16];
    int32_t i;
    int32_t row_offset = 0;
    const int32_t shift = 0;
    if (abs(col_offset) > 1) {
        row_offset = 1;
        if (mi_row & 0x01 && xd->n8_h < n8_h_8) --row_offset;
    }
    const int32_t use_step_16 = (xd->n8_h >= 16);
    (void)mi_col;

    for (i = 0; i < end_mi;) {
        const ModeInfo *const candidate_mi =
            xd->mi[(row_offset + i) * xd->mi_stride + col_offset];
        const MbModeInfo *const candidate = &candidate_mi->mbmi;
        const int32_t candidate_bsize = candidate->sb_type;
        ASSERT(candidate_bsize < BlockSizeS_ALL);
        const int32_t n8_h = mi_size_high[candidate_bsize];
        int32_t len = AOMMIN(xd->n8_h, n8_h);
        if (use_step_16)
            len = AOMMAX(n8_h_16, len);
        else if (abs(col_offset) > 1)
            len = AOMMAX(len, n8_h_8);

        int32_t weight = 2;
        if (xd->n8_h >= n8_h_8 && xd->n8_h <= n8_h) {
            int32_t inc = AOMMIN(-max_col_offset + col_offset + 1,
                mi_size_wide[candidate_bsize]);
            // Obtain range used in weight calculation.
            weight = AOMMAX(weight, (inc << shift));
            // Update processed cols.
            *processed_cols = inc - col_offset - 1;
        }

        add_ref_mv_candidate(candidate_mi, candidate, rf, refmv_count,
            ref_match_count, newmv_count, ref_mv_stack, len,
#if USE_CUR_GM_REFMV
            gm_mv_candidates, cm->global_motion,
#endif  // USE_CUR_GM_REFMV
            col_offset, weight);

        i += len;
    }
}

static void scan_blk_mbmi(const Av1Common *cm, const MacroBlockD *xd,
    const int32_t mi_row, const int32_t mi_col,
    const MvReferenceFrame rf[2], int32_t row_offset,
    int32_t col_offset,
    CandidateMv ref_mv_stack[][MAX_REF_MV_STACK_SIZE],
    uint8_t ref_match_count[MODE_CTX_REF_FRAMES],
    uint8_t newmv_count[MODE_CTX_REF_FRAMES],
#if USE_CUR_GM_REFMV
    IntMv *gm_mv_candidates,
#endif  // USE_CUR_GM_REFMV
    uint8_t refmv_count[MODE_CTX_REF_FRAMES]) {
    const TileInfo *const tile = &xd->tile;
    Position mi_pos;

    mi_pos.row = row_offset;
    mi_pos.col = col_offset;

    if (is_inside(tile, mi_col, mi_row, cm->mi_rows, &mi_pos)) {
        const ModeInfo *const candidate_mi =
            xd->mi[mi_pos.row * xd->mi_stride + mi_pos.col];
        const MbModeInfo *const candidate = &candidate_mi->mbmi;
        const int32_t len = mi_size_wide[BLOCK_8X8];

        add_ref_mv_candidate(candidate_mi, candidate, rf, refmv_count,
            ref_match_count, newmv_count, ref_mv_stack, len,
#if USE_CUR_GM_REFMV
            gm_mv_candidates, cm->global_motion,
#endif  // USE_CUR_GM_REFMV
            mi_pos.col, 2);
    }  // Analyze a single 8x8 block motion information.
}

static int32_t has_top_right(const Av1Common *cm, const MacroBlockD *xd,
    int32_t mi_row, int32_t mi_col, int32_t bs) {

    (void)xd;
    (void)cm;
    const int32_t sb_mi_size = mi_size_wide[cm->p_pcs_ptr->sequence_control_set_ptr->sb_size];
    const int32_t mask_row = mi_row & (sb_mi_size - 1);
    const int32_t mask_col = mi_col & (sb_mi_size - 1);

    if (bs > mi_size_wide[BLOCK_64X64]) return 0;

    // In a split partition all apart from the bottom right has a top right
    int32_t has_tr = !((mask_row & bs) && (mask_col & bs));

    // bs > 0 and bs is a power of 2
    assert(bs > 0 && !(bs & (bs - 1)));

    // For each 4x4 group of blocks, when the bottom right is decoded the blocks
    // to the right have not been decoded therefore the bottom right does
    // not have a top right
    while (bs < sb_mi_size) {
        if (mask_col & bs) {
            if ((mask_col & (2 * bs)) && (mask_row & (2 * bs))) {
                has_tr = 0;
                break;
            }
        }
        else {
            break;
        }
        bs <<= 1;
    }

    // The left hand of two vertical rectangles always has a top right (as the
    // block above will have been decoded)
    if (xd->n8_w < xd->n8_h)
        if (!xd->is_sec_rect) has_tr = 1;

    // The bottom of two horizontal rectangles never has a top right (as the block
    // to the right won't have been decoded)
    if (xd->n8_w > xd->n8_h)
        if (xd->is_sec_rect) has_tr = 0;

    // The bottom left square of a Vertical A (in the old format) does
    // not have a top right as it is decoded before the right hand
    // rectangle of the partition
    if (xd->mi[0]->mbmi.partition == PARTITION_VERT_A) {
        if (xd->n8_w == xd->n8_h)
            if (mask_row & bs) has_tr = 0;
    }

    return has_tr;
}
static INLINE int32_t find_valid_row_offset(const TileInfo *const tile, int32_t mi_row,
    int32_t mi_rows, int32_t row_offset) {
    const int32_t dependent_horz_tile_flag = 0;
    if (dependent_horz_tile_flag && !tile->tg_horz_boundary)
        return clamp(row_offset, -mi_row, mi_rows - mi_row - 1);
    else
        return clamp(row_offset, tile->mi_row_start - mi_row,
            tile->mi_row_end - mi_row - 1);
}

static INLINE int32_t find_valid_col_offset(const TileInfo *const tile, int32_t mi_col,
    int32_t col_offset) {
    return clamp(col_offset, tile->mi_col_start - mi_col,
        tile->mi_col_end - mi_col - 1);
}

void setup_ref_mv_list(
    const Av1Common *cm, const MacroBlockD *xd, MvReferenceFrame ref_frame,
    uint8_t refmv_count[MODE_CTX_REF_FRAMES],
    CandidateMv ref_mv_stack[][MAX_REF_MV_STACK_SIZE],
    IntMv mv_ref_list[][MAX_MV_REF_CANDIDATES],
    IntMv *gm_mv_candidates,
    int32_t mi_row, int32_t mi_col, int16_t *mode_context)
{

    const int32_t bs = AOMMAX(xd->n8_w, xd->n8_h);
    const int32_t has_tr = has_top_right(cm, xd, mi_row, mi_col, bs);
    MvReferenceFrame rf[2];

    const TileInfo *const tile = &xd->tile;
    int32_t max_row_offset = 0, max_col_offset = 0;
    const int32_t row_adj = (xd->n8_h < mi_size_high[BLOCK_8X8]) && (mi_row & 0x01);
    const int32_t col_adj = (xd->n8_w < mi_size_wide[BLOCK_8X8]) && (mi_col & 0x01);
    int32_t processed_rows = 0;
    int32_t processed_cols = 0;

    av1_set_ref_frame(rf, ref_frame);
    mode_context[ref_frame] = 0;
    refmv_count[ref_frame] = 0;

    // Find valid maximum row/col offset.
    if (xd->up_available) {
        max_row_offset = -(MVREF_ROWS << 1) + row_adj;

        if (xd->n8_h < mi_size_high[BLOCK_8X8])
            max_row_offset = -(2 << 1) + row_adj;

        max_row_offset =
            find_valid_row_offset(tile, mi_row, cm->mi_rows, max_row_offset);
    }

    if (xd->left_available) {
        max_col_offset = -(MVREF_COLS << 1) + col_adj;

        if (xd->n8_w < mi_size_wide[BLOCK_8X8])
            max_col_offset = -(2 << 1) + col_adj;

        max_col_offset = find_valid_col_offset(tile, mi_col, max_col_offset);
    }

    uint8_t ref_match_count[MODE_CTX_REF_FRAMES] = { 0 };
    uint8_t col_match_count[MODE_CTX_REF_FRAMES] = { 0 };
    uint8_t row_match_count[MODE_CTX_REF_FRAMES] = { 0 };
    uint8_t newmv_count[MODE_CTX_REF_FRAMES] = { 0 };


    //CHKN-------------    ROW-1

    // Scan the first above row mode info. row_offset = -1;
    if (abs(max_row_offset) >= 1)
        scan_row_mbmi(cm, xd, mi_row, mi_col, rf, -1, ref_mv_stack, refmv_count,
            row_match_count, newmv_count,
#if USE_CUR_GM_REFMV
            gm_mv_candidates,
#endif  // USE_CUR_GM_REFMV
            max_row_offset, &processed_rows);


    //CHKN-------------    COL-1
    // Scan the first left column mode info. col_offset = -1;
    if (abs(max_col_offset) >= 1)
        scan_col_mbmi(cm, xd, mi_row, mi_col, rf, -1, ref_mv_stack, refmv_count,
            col_match_count, newmv_count,
#if USE_CUR_GM_REFMV
            gm_mv_candidates,
#endif  // USE_CUR_GM_REFMV
            max_col_offset, &processed_cols);



    //CHKN-------------    TOP-RIGHT

    // Check top-right boundary
    if (has_tr)
        scan_blk_mbmi(cm, xd, mi_row, mi_col, rf, -1, xd->n8_w, ref_mv_stack,
            row_match_count, newmv_count,
#if USE_CUR_GM_REFMV
            gm_mv_candidates,
#endif  // USE_CUR_GM_REFMV
            refmv_count);

    uint8_t nearest_match[MODE_CTX_REF_FRAMES];
    uint8_t nearest_refmv_count[MODE_CTX_REF_FRAMES];

    nearest_match[ref_frame] =
        (row_match_count[ref_frame] > 0) + (col_match_count[ref_frame] > 0);
    nearest_refmv_count[ref_frame] = refmv_count[ref_frame];

    // TODO(yunqing): for comp_search, do it for all 3 cases.
    for (int32_t idx = 0; idx < nearest_refmv_count[ref_frame]; ++idx)
        ref_mv_stack[ref_frame][idx].weight += REF_CAT_LEVEL;


    //-------------------------- TMVP --------------------------
    //CHKN  TMVP - get canididates from reference frames- orderHint has to be on,
    // in order to scale the vectors.
    //CHKN this checks all colocated block(s) + extra 3 positions. this changes
    // the mode context too.

    //-------------------------- TMVP --------------------------

    //CHKN------------- TOP-LEFT
    uint8_t dummy_newmv_count[MODE_CTX_REF_FRAMES] = { 0 };

    // Scan the second outer area.
    scan_blk_mbmi(cm, xd, mi_row, mi_col, rf, -1, -1, ref_mv_stack,
        row_match_count, dummy_newmv_count,
#if USE_CUR_GM_REFMV
        gm_mv_candidates,
#endif  // USE_CUR_GM_REFMV
        refmv_count);


    //CHKN-------------    ROW-3  COL-3     ROW-5   COL-5
    for (int32_t idx = 2; idx <= MVREF_ROWS; ++idx) {
        const int32_t row_offset = -(idx << 1) + 1 + row_adj;
        const int32_t col_offset = -(idx << 1) + 1 + col_adj;

        if (abs(row_offset) <= abs(max_row_offset) &&
            abs(row_offset) > processed_rows)
            scan_row_mbmi(cm, xd, mi_row, mi_col, rf, row_offset, ref_mv_stack,
                refmv_count, row_match_count, dummy_newmv_count,
#if USE_CUR_GM_REFMV
                gm_mv_candidates,
#endif  // USE_CUR_GM_REFMV
                max_row_offset, &processed_rows);

        if (abs(col_offset) <= abs(max_col_offset) &&
            abs(col_offset) > processed_cols)
            scan_col_mbmi(cm, xd, mi_row, mi_col, rf, col_offset, ref_mv_stack,
                refmv_count, col_match_count, dummy_newmv_count,
#if USE_CUR_GM_REFMV
                gm_mv_candidates,
#endif  // USE_CUR_GM_REFMV
                max_col_offset, &processed_cols);
    }

    //---------- Mode Context Derivation based on 3 counters -------------
    ref_match_count[ref_frame] =
        (row_match_count[ref_frame] > 0) + (col_match_count[ref_frame] > 0);

    switch (nearest_match[ref_frame]) {
    case 0:
        mode_context[ref_frame] |= 0;
        if (ref_match_count[ref_frame] >= 1) mode_context[ref_frame] |= 1;
        if (ref_match_count[ref_frame] == 1)
            mode_context[ref_frame] |= (1 << REFMV_OFFSET);
        else if (ref_match_count[ref_frame] >= 2)
            mode_context[ref_frame] |= (2 << REFMV_OFFSET);
        break;
    case 1:
        mode_context[ref_frame] |= (newmv_count[ref_frame] > 0) ? 2 : 3;
        if (ref_match_count[ref_frame] == 1)
            mode_context[ref_frame] |= (3 << REFMV_OFFSET);
        else if (ref_match_count[ref_frame] >= 2)
            mode_context[ref_frame] |= (4 << REFMV_OFFSET);
        break;
    case 2:
    default:
        if (newmv_count[ref_frame] >= 1)
            mode_context[ref_frame] |= 4;
        else
            mode_context[ref_frame] |= 5;

        mode_context[ref_frame] |= (5 << REFMV_OFFSET);
        break;
    }
    //---------- Mode Context Derivation based on 3 counters -------------

    // Rank the likelihood and assign nearest and near mvs.
    int32_t len = nearest_refmv_count[ref_frame];
    while (len > 0) {
        int32_t nr_len = 0;
        for (int32_t idx = 1; idx < len; ++idx) {
            if (ref_mv_stack[ref_frame][idx - 1].weight <
                ref_mv_stack[ref_frame][idx].weight) {
                CandidateMv tmp_mv = ref_mv_stack[ref_frame][idx - 1];
                ref_mv_stack[ref_frame][idx - 1] = ref_mv_stack[ref_frame][idx];
                ref_mv_stack[ref_frame][idx] = tmp_mv;
                nr_len = idx;
            }
        }
        len = nr_len;
    }

    len = refmv_count[ref_frame];
    while (len > nearest_refmv_count[ref_frame]) {
        int32_t nr_len = nearest_refmv_count[ref_frame];
        for (int32_t idx = nearest_refmv_count[ref_frame] + 1; idx < len; ++idx) {
            if (ref_mv_stack[ref_frame][idx - 1].weight <
                ref_mv_stack[ref_frame][idx].weight) {
                CandidateMv tmp_mv = ref_mv_stack[ref_frame][idx - 1];
                ref_mv_stack[ref_frame][idx - 1] = ref_mv_stack[ref_frame][idx];
                ref_mv_stack[ref_frame][idx] = tmp_mv;
                nr_len = idx;
            }
        }
        len = nr_len;
    }


    //CHKN finish the Tables.

    if (rf[1] > NONE_FRAME) {

        // TODO(jingning, yunqing): Refactor and consolidate the compound and
        // single reference frame modes. Reduce unnecessary redundancy.

        //CHKN we get here only when refMVCount=0 or 1

        if (refmv_count[ref_frame] < 2) {
            IntMv ref_id[2][2], ref_diff[2][2];
            int32_t ref_id_count[2] = { 0 }, ref_diff_count[2] = { 0 };

            int32_t mi_width = AOMMIN(mi_size_wide[BLOCK_64X64], xd->n8_w);
            mi_width = AOMMIN(mi_width, cm->mi_cols - mi_col);
            int32_t mi_height = AOMMIN(mi_size_high[BLOCK_64X64], xd->n8_h);
            mi_height = AOMMIN(mi_height, cm->mi_rows - mi_row);
            int32_t mi_size = AOMMIN(mi_width, mi_height);

            //CHKN  scan ROW=-1 again but with more relaxed constraints
            for (int32_t idx = 0; abs(max_row_offset) >= 1 && idx < mi_size;) {
                const ModeInfo *const candidate_mi = xd->mi[-xd->mi_stride + idx];
                const MbModeInfo *const candidate = &candidate_mi->mbmi;
                const int32_t candidate_bsize = candidate->sb_type;

                for (int32_t rf_idx = 0; rf_idx < 2; ++rf_idx) {
                    MvReferenceFrame can_rf = candidate->ref_frame[rf_idx];

                    for (int32_t cmp_idx = 0; cmp_idx < 2; ++cmp_idx) {
                        if (can_rf == rf[cmp_idx] && ref_id_count[cmp_idx] < 2) {
                            ref_id[cmp_idx][ref_id_count[cmp_idx]] = candidate->mv[rf_idx];
                            ++ref_id_count[cmp_idx];
                        }
                        else if (can_rf > INTRA_FRAME && ref_diff_count[cmp_idx] < 2) {
                            IntMv this_mv = candidate->mv[rf_idx];
                            if (cm->ref_frame_sign_bias[can_rf] !=
                                cm->ref_frame_sign_bias[rf[cmp_idx]]) {
                                this_mv.as_mv.row = -this_mv.as_mv.row;
                                this_mv.as_mv.col = -this_mv.as_mv.col;
                            }
                            ref_diff[cmp_idx][ref_diff_count[cmp_idx]] = this_mv;
                            ++ref_diff_count[cmp_idx];
                        }
                    }
                }
                idx += mi_size_wide[candidate_bsize];
            }

            //CHKN  scan COL=-1 again but with more relaxed constraints
            for (int32_t idx = 0; abs(max_col_offset) >= 1 && idx < mi_size;) {
                const ModeInfo *const candidate_mi = xd->mi[idx * xd->mi_stride - 1];
                const MbModeInfo *const candidate = &candidate_mi->mbmi;
                const int32_t candidate_bsize = candidate->sb_type;

                for (int32_t rf_idx = 0; rf_idx < 2; ++rf_idx) {
                    MvReferenceFrame can_rf = candidate->ref_frame[rf_idx];

                    for (int32_t cmp_idx = 0; cmp_idx < 2; ++cmp_idx) {
                        if (can_rf == rf[cmp_idx] && ref_id_count[cmp_idx] < 2) {
                            ref_id[cmp_idx][ref_id_count[cmp_idx]] = candidate->mv[rf_idx];
                            ++ref_id_count[cmp_idx];
                        }
                        else if (can_rf > INTRA_FRAME && ref_diff_count[cmp_idx] < 2) {
                            IntMv this_mv = candidate->mv[rf_idx];
                            if (cm->ref_frame_sign_bias[can_rf] !=
                                cm->ref_frame_sign_bias[rf[cmp_idx]]) {
                                this_mv.as_mv.row = -this_mv.as_mv.row;
                                this_mv.as_mv.col = -this_mv.as_mv.col;
                            }
                            ref_diff[cmp_idx][ref_diff_count[cmp_idx]] = this_mv;
                            ++ref_diff_count[cmp_idx];
                        }
                    }
                }
                idx += mi_size_high[candidate_bsize];
            }

            //CHKN I dont think we need the third

            // Build up the compound mv predictor
            IntMv comp_list[3][2];

            for (int32_t idx = 0; idx < 2; ++idx) {
                int32_t comp_idx = 0;
                for (int32_t list_idx = 0; list_idx < ref_id_count[idx]
                    && comp_idx < 3; ++list_idx, ++comp_idx)
                    comp_list[comp_idx][idx] = ref_id[idx][list_idx];

                for (int32_t list_idx = 0; list_idx < ref_diff_count[idx]
                    && comp_idx < 3; ++list_idx, ++comp_idx)
                    comp_list[comp_idx][idx] = ref_diff[idx][list_idx];

                for (; comp_idx < 3; ++comp_idx)
                    comp_list[comp_idx][idx] = gm_mv_candidates[idx];
            }

            //CHKN fill the stack, increment the counter
            if (refmv_count[ref_frame]) { //CHKN RefMvCount=1
                assert(refmv_count[ref_frame] == 1);
                if (comp_list[0][0].as_int == ref_mv_stack[ref_frame][0].this_mv.as_int
                    &&  comp_list[0][1].as_int == ref_mv_stack[ref_frame][0].comp_mv.as_int) {
                    ref_mv_stack[ref_frame][refmv_count[ref_frame]].this_mv = comp_list[1][0];
                    ref_mv_stack[ref_frame][refmv_count[ref_frame]].comp_mv = comp_list[1][1];
                }
                else {
                    ref_mv_stack[ref_frame][refmv_count[ref_frame]].this_mv = comp_list[0][0];
                    ref_mv_stack[ref_frame][refmv_count[ref_frame]].comp_mv = comp_list[0][1];
                }
                ref_mv_stack[ref_frame][refmv_count[ref_frame]].weight = 2;
                ++refmv_count[ref_frame];

            }
            else {//CHKN RefMvCount=0
                for (int32_t idx = 0; idx < MAX_MV_REF_CANDIDATES; ++idx) {
                    ref_mv_stack[ref_frame][refmv_count[ref_frame]].this_mv = comp_list[idx][0];
                    ref_mv_stack[ref_frame][refmv_count[ref_frame]].comp_mv = comp_list[idx][1];
                    ref_mv_stack[ref_frame][refmv_count[ref_frame]].weight = 2;
                    ++refmv_count[ref_frame];
                }
            }
        }

        assert(refmv_count[ref_frame] >= 2);

        for (int32_t idx = 0; idx < refmv_count[ref_frame]; ++idx) {
            clamp_mv_ref(&ref_mv_stack[ref_frame][idx].this_mv.as_mv,
                xd->n8_w << MI_SIZE_LOG2, xd->n8_h << MI_SIZE_LOG2, xd);
            clamp_mv_ref(&ref_mv_stack[ref_frame][idx].comp_mv.as_mv,
                xd->n8_w << MI_SIZE_LOG2, xd->n8_h << MI_SIZE_LOG2, xd);
        }
    }
    else {
        // Handle single reference frame extension
        int32_t mi_width = AOMMIN(mi_size_wide[BLOCK_64X64], xd->n8_w);
        mi_width = AOMMIN(mi_width, cm->mi_cols - mi_col);
        int32_t mi_height = AOMMIN(mi_size_high[BLOCK_64X64], xd->n8_h);
        mi_height = AOMMIN(mi_height, cm->mi_rows - mi_row);
        int32_t mi_size = AOMMIN(mi_width, mi_height);


        //CHKn if count is still < 2, re-scan ROW=-1 with less constraints.
        //     Order is already fixed. the added candidates are stored as we go at the bottom of the Stack.
        //CHKN TODO: confirm this could be avoided if we have already 2(DRL:OFF), or 4(DRL:ON) candidates
        for (int32_t idx = 0; abs(max_row_offset) >= 1 && idx < mi_size &&
            refmv_count[ref_frame] < MAX_MV_REF_CANDIDATES;) {

            const ModeInfo *const candidate_mi = xd->mi[-xd->mi_stride + idx];
            const MbModeInfo *const candidate = &candidate_mi->mbmi;
            const int32_t candidate_bsize = candidate->sb_type;

            // TODO(jingning): Refactor the following code.
            for (int32_t rf_idx = 0; rf_idx < 2; ++rf_idx) {
                if (candidate->ref_frame[rf_idx] > INTRA_FRAME) {
                    IntMv this_mv = candidate->mv[rf_idx];
                    if (cm->ref_frame_sign_bias[candidate->ref_frame[rf_idx]] !=
                        cm->ref_frame_sign_bias[ref_frame]) {
                        this_mv.as_mv.row = -this_mv.as_mv.row;
                        this_mv.as_mv.col = -this_mv.as_mv.col;
                    }
                    int32_t stack_idx;
                    for (stack_idx = 0; stack_idx < refmv_count[ref_frame]; ++stack_idx) {
                        IntMv stack_mv = ref_mv_stack[ref_frame][stack_idx].this_mv;
                        if (this_mv.as_int == stack_mv.as_int) break;
                    }

                    if (stack_idx == refmv_count[ref_frame]) {
                        ref_mv_stack[ref_frame][stack_idx].this_mv = this_mv;

                        // TODO(jingning): Set an arbitrary small number here. The weight
                        // doesn't matter as long as it is properly initialized.
                        ref_mv_stack[ref_frame][stack_idx].weight = 2;
                        ++refmv_count[ref_frame];
                    }
                }
            }
            idx += mi_size_wide[candidate_bsize];
        }

        //CHKn if count is still < 2, re-scan COL=-1 with less constraints. the added candidates are stored as we go at the bottom of the Stack.
        for (int32_t idx = 0; abs(max_col_offset) >= 1 && idx < mi_size &&
            refmv_count[ref_frame] < MAX_MV_REF_CANDIDATES;) {
            const ModeInfo *const candidate_mi = xd->mi[idx * xd->mi_stride - 1];
            const MbModeInfo *const candidate = &candidate_mi->mbmi;
            const int32_t candidate_bsize = candidate->sb_type;

            // TODO(jingning): Refactor the following code.
            for (int32_t rf_idx = 0; rf_idx < 2; ++rf_idx) {
                if (candidate->ref_frame[rf_idx] > INTRA_FRAME) {
                    IntMv this_mv = candidate->mv[rf_idx];
                    if (cm->ref_frame_sign_bias[candidate->ref_frame[rf_idx]] !=
                        cm->ref_frame_sign_bias[ref_frame]) {
                        this_mv.as_mv.row = -this_mv.as_mv.row;
                        this_mv.as_mv.col = -this_mv.as_mv.col;
                    }
                    int32_t stack_idx;
                    for (stack_idx = 0; stack_idx < refmv_count[ref_frame]; ++stack_idx) {
                        IntMv stack_mv = ref_mv_stack[ref_frame][stack_idx].this_mv;
                        if (this_mv.as_int == stack_mv.as_int) break;
                    }

                    if (stack_idx == refmv_count[ref_frame]) {
                        ref_mv_stack[ref_frame][stack_idx].this_mv = this_mv;

                        // TODO(jingning): Set an arbitrary small number here. The weight
                        // doesn't matter as long as it is properly initialized.
                        ref_mv_stack[ref_frame][stack_idx].weight = 2;
                        ++refmv_count[ref_frame];
                    }
                }
            }
            idx += mi_size_high[candidate_bsize];
        }


        //CHKN  THIS IS a Single Reference case

        //CHKN if the stack has at least 2 cand, then copy the top 2 to the final mvp Table,
        // if the stack has less than 2, use Gm to fill the mvpTable.
        //     the stack counter remains 0 or 1 in this case.

        for (int32_t idx = refmv_count[ref_frame]; idx < MAX_MV_REF_CANDIDATES; ++idx)
            mv_ref_list[rf[0]][idx].as_int = gm_mv_candidates[0].as_int;

        for (int32_t idx = 0; idx < refmv_count[ref_frame]; ++idx) {
            clamp_mv_ref(&ref_mv_stack[ref_frame][idx].this_mv.as_mv,
                xd->n8_w << MI_SIZE_LOG2, xd->n8_h << MI_SIZE_LOG2, xd);
        }

        for (int32_t idx = 0;
            idx < AOMMIN(MAX_MV_REF_CANDIDATES, refmv_count[ref_frame]); ++idx) {
            mv_ref_list[rf[0]][idx].as_int =
                ref_mv_stack[ref_frame][idx].this_mv.as_int;
        }
    }
    (void)nearest_match;
}

static INLINE void integer_mv_precision(MV *mv) {
    int32_t mod = (mv->row % 8);
    if (mod != 0) {
        mv->row -= (int16_t)mod;
        if (abs(mod) > 4) {
            if (mod > 0) {
                mv->row += 8;
            }
            else {
                mv->row -= 8;
            }
        }
    }

    mod = (mv->col % 8);
    if (mod != 0) {
        mv->col -= (int16_t)mod;
        if (abs(mod) > 4) {
            if (mod > 0) {
                mv->col += 8;
            }
            else {
                mv->col -= 8;
            }
        }
    }
}

static INLINE IntMv gm_get_motion_vector(
    const EbWarpedMotionParams *gm,
    int32_t allow_hp,
    BlockSize bsize,
    int32_t mi_col, int32_t mi_row,
    int32_t is_integer)

{
    IntMv res;

    res.as_int = 0;


    (void)bsize;
    (void)mi_col;
    (void)mi_row;
    (void)allow_hp;

    if (gm->wmtype <= TRANSLATION) {
        // All global motion vectors are stored with WARPEDMODEL_PREC_BITS (16)
        // bits of fractional precision. The offset for a translation is stored in
        // entries 0 and 1. For translations, all but the top three (two if
        // cm->allow_high_precision_mv is false) fractional bits are always zero.
        //
        // After the right shifts, there are 3 fractional bits of precision. If
        // allow_hp is false, the bottom bit is always zero (so we don't need a
        // call to convert_to_trans_prec here)
        res.as_mv.row = (int16_t)(gm->wmmat[0] >> GM_TRANS_ONLY_PREC_DIFF);
        res.as_mv.col = (int16_t)(gm->wmmat[1] >> GM_TRANS_ONLY_PREC_DIFF);
        assert(IMPLIES(1 & (res.as_mv.row | res.as_mv.col), allow_hp));

        if (is_integer) {
            integer_mv_precision(&res.as_mv);
        }

        return res;
    }
    /*else
        printf("ERROR - INVALID GLOBAL MV - GLOBAL ROTATION AND AFFINE ARE NOT SUPPORTED FOR NOW!!");*/

    return res;

}

void generate_av1_mvp_table(
    ModeDecisionContext_t            *context_ptr,
    CodingUnit_t                     *cu_ptr,
    const BlockGeom                  *blk_geom,
    uint16_t                          cu_origin_x,
    uint16_t                          cu_origin_y,
    MvReferenceFrame                 *ref_frames,
    uint32_t                          tot_refs,
    PictureControlSet_t              *picture_control_set_ptr)
{
    int32_t mi_row = cu_origin_y >> MI_SIZE_LOG2;
    int32_t mi_col = cu_origin_x >> MI_SIZE_LOG2;
    Av1Common  *cm = picture_control_set_ptr->parent_pcs_ptr->av1_cm;
    MacroBlockD  *xd = cu_ptr->av1xd;
    xd->n8_w = blk_geom->bwidth >> MI_SIZE_LOG2;
    xd->n8_h = blk_geom->bheight >> MI_SIZE_LOG2;
    xd->n4_w = blk_geom->bwidth >> MI_SIZE_LOG2;
    xd->n4_h = blk_geom->bheight >> MI_SIZE_LOG2;
    BlockSize bsize = blk_geom->bsize;
    const int32_t bw = mi_size_wide[bsize];
    const int32_t bh = mi_size_high[bsize];

    xd->mb_to_top_edge = -((mi_row * MI_SIZE) * 8);
    xd->mb_to_bottom_edge = ((cm->mi_rows - bh - mi_row) * MI_SIZE) * 8;
    xd->mb_to_left_edge = -((mi_col * MI_SIZE) * 8);
    xd->mb_to_right_edge = ((cm->mi_cols - bw - mi_col) * MI_SIZE) * 8;

    memset(xd->ref_mv_count, 0, sizeof(xd->ref_mv_count));
    memset(context_ptr->md_local_cu_unit[blk_geom->blkidx_mds].ed_ref_mv_stack, 0, sizeof(context_ptr->md_local_cu_unit[blk_geom->blkidx_mds].ed_ref_mv_stack));

    xd->up_available = (mi_row > 0);
    xd->left_available = (mi_col > 0);

    xd->n8_h = bh;
    xd->n8_w = bw;
    xd->is_sec_rect = 0;
    if (xd->n8_w < xd->n8_h) {
        // Only mark is_sec_rect as 1 for the last block.
        // For PARTITION_VERT_4, it would be (0, 0, 0, 1);
        // For other partitions, it would be (0, 1).
        if (!((mi_col + xd->n8_w) & (xd->n8_h - 1))) xd->is_sec_rect = 1;
    }

    if (xd->n8_w > xd->n8_h)
        if (mi_row & (xd->n8_w - 1)) xd->is_sec_rect = 1;


    xd->tile.mi_col_start = 0;
    xd->tile.mi_col_end = cm->mi_cols;
    xd->tile.mi_row_start = 0;
    xd->tile.mi_row_end = cm->mi_rows;

    //these could be done at init time
    xd->mi_stride = picture_control_set_ptr->parent_pcs_ptr->sequence_control_set_ptr->picture_width_in_sb*(BLOCK_SIZE_64 / 4);
    const int32_t offset = mi_row * xd->mi_stride + mi_col;
    xd->mi = picture_control_set_ptr->mi_grid_base + offset;

    xd->mi[0]->mbmi.partition = from_shape_to_part[blk_geom->shape];

    uint32_t refIt;
    for (refIt = 0; refIt < tot_refs; ++refIt) {

        MvReferenceFrame ref_frame = ref_frames[refIt];
        IntMv zeromv[2] = { {0}, {0} };

        MvReferenceFrame rf[2];
        av1_set_ref_frame(rf, ref_frame);

        if (ref_frame != INTRA_FRAME) {
            zeromv[0].as_int =
                gm_get_motion_vector(&picture_control_set_ptr->parent_pcs_ptr->global_motion[rf[0]],
                    picture_control_set_ptr->parent_pcs_ptr->allow_high_precision_mv, bsize, mi_col, mi_row,
                    picture_control_set_ptr->parent_pcs_ptr->cur_frame_force_integer_mv)
                .as_int;
            zeromv[1].as_int = (rf[1] != NONE_FRAME)
                ? gm_get_motion_vector(&picture_control_set_ptr->parent_pcs_ptr->global_motion[rf[1]],
                    picture_control_set_ptr->parent_pcs_ptr->allow_high_precision_mv,
                    bsize, mi_col, mi_row,
                    picture_control_set_ptr->parent_pcs_ptr->cur_frame_force_integer_mv)
                .as_int
                : 0;
        }
        else {
            zeromv[0].as_int = zeromv[1].as_int = 0;
        }
        setup_ref_mv_list(cm,
            xd,
            ref_frame,
            xd->ref_mv_count,
            context_ptr->md_local_cu_unit[blk_geom->blkidx_mds].ed_ref_mv_stack,
            cu_ptr->ref_mvs,
            zeromv,
            mi_row,
            mi_col,
            cu_ptr->inter_mode_ctx);

    }

}
void get_av1_mv_pred_drl(
    ModeDecisionContext_t            *context_ptr,
    CodingUnit_t      *cu_ptr,
    MvReferenceFrame   ref_frame,
    uint8_t            is_compound,
    PredictionMode     mode,
    uint8_t            drl_index,    //valid value of drl_index
    IntMv              nearestmv[2],
    IntMv              nearmv[2],
    IntMv              ref_mv[2])
{
    MacroBlockD*  xd = cu_ptr->av1xd;

    if (!is_compound &&  mode != GLOBALMV) {

        //av1_find_best_ref_mvs(allow_hp, ref_mvs[mbmi->ref_frame[0]], &nearestmv[0], &nearmv[0], cm->cur_frame_force_integer_mv);

        nearestmv[0] = cu_ptr->ref_mvs[ref_frame][0];
        nearmv[0] = cu_ptr->ref_mvs[ref_frame][1];
    }

    if (is_compound && mode != GLOBAL_GLOBALMV) {
        int32_t ref_mv_idx = drl_index + 1;
        nearestmv[0] = context_ptr->md_local_cu_unit[cu_ptr->mds_idx].ed_ref_mv_stack[ref_frame][0].this_mv;
        nearestmv[1] = context_ptr->md_local_cu_unit[cu_ptr->mds_idx].ed_ref_mv_stack[ref_frame][0].comp_mv;
        nearmv[0]    = context_ptr->md_local_cu_unit[cu_ptr->mds_idx].ed_ref_mv_stack[ref_frame][ref_mv_idx].this_mv;
        nearmv[1]    = context_ptr->md_local_cu_unit[cu_ptr->mds_idx].ed_ref_mv_stack[ref_frame][ref_mv_idx].comp_mv;
    }
    else if (drl_index > 0 && mode == NEARMV) {
        ASSERT((1 + drl_index) < MAX_REF_MV_STACK_SIZE);
        IntMv cur_mv = context_ptr->md_local_cu_unit[cu_ptr->mds_idx].ed_ref_mv_stack[ref_frame][1 + drl_index].this_mv;
        nearmv[0] = cur_mv;
    }

    ref_mv[0] = nearestmv[0];
    ref_mv[1] = nearestmv[1];

    if (is_compound) {
        int32_t ref_mv_idx = drl_index;
        // Special case: NEAR_NEWMV and NEW_NEARMV modes use
        // 1 + mbmi->ref_mv_idx (like NEARMV) instead of
        // mbmi->ref_mv_idx (like NEWMV)
        if (mode == NEAR_NEWMV || mode == NEW_NEARMV)
            ref_mv_idx = 1 + drl_index;


        // TODO(jingning, yunqing): Do we need a lower_mv_precision() call here?
        if (compound_ref0_mode(mode) == NEWMV)
            ref_mv[0] = context_ptr->md_local_cu_unit[cu_ptr->mds_idx].ed_ref_mv_stack[ref_frame][ref_mv_idx].this_mv;

        if (compound_ref1_mode(mode) == NEWMV)
            ref_mv[1] = context_ptr->md_local_cu_unit[cu_ptr->mds_idx].ed_ref_mv_stack[ref_frame][ref_mv_idx].comp_mv;

    }
    else {
        if (mode == NEWMV) {
            if (xd->ref_mv_count[ref_frame] > 1)
                ref_mv[0] = context_ptr->md_local_cu_unit[cu_ptr->mds_idx].ed_ref_mv_stack[ref_frame][drl_index].this_mv;
        }
    }


}

void enc_pass_av1_mv_pred(
    ModeDecisionContext_t            *md_context_ptr,
    CodingUnit_t                     *cu_ptr,
    const BlockGeom                  *blk_geom,
    uint16_t                          cu_origin_x,
    uint16_t                          cu_origin_y,
    PictureControlSet_t              *picture_control_set_ptr,
    MvReferenceFrame                  ref_frame,
    uint8_t                           is_compound,
    PredictionMode                    mode,
    IntMv                             ref_mv[2]){ //[OUT]

    (void)mode;
    IntMv    nearestmv[2], nearmv[2];

    generate_av1_mvp_table(
        md_context_ptr,
        cu_ptr,
        blk_geom,
        cu_origin_x,
        cu_origin_y,
        &ref_frame,
        1,
        picture_control_set_ptr);

    get_av1_mv_pred_drl(
        md_context_ptr,
        cu_ptr,
        ref_frame,
        is_compound,
        cu_ptr->pred_mode,
        cu_ptr->drl_index,
        nearestmv,
        nearmv,
        ref_mv);
}

void update_av1_mi_map(
    CodingUnit_t                   *cu_ptr,
    uint32_t                        cu_origin_x,
    uint32_t                        cu_origin_y,
    const BlockGeom                *blk_geom,
    PictureControlSet_t            *picture_control_set_ptr)
{
    uint32_t mi_stride = picture_control_set_ptr->parent_pcs_ptr->sequence_control_set_ptr->picture_width_in_sb*(BLOCK_SIZE_64 >> MI_SIZE_LOG2);
    int32_t mi_row = cu_origin_y >> MI_SIZE_LOG2;
    int32_t mi_col = cu_origin_x >> MI_SIZE_LOG2;

    const int32_t offset = mi_row * mi_stride + mi_col;
    ModeInfo *miPtr = *(picture_control_set_ptr->mi_grid_base + offset);
    MvReferenceFrame rf[2];
    av1_set_ref_frame(rf, cu_ptr->prediction_unit_array->ref_frame_type);
    uint8_t  miX, miY;

    for (miY = 0; miY < (blk_geom->bheight >> MI_SIZE_LOG2); miY++) {

        for (miX = 0; miX < (blk_geom->bwidth >> MI_SIZE_LOG2); miX++) {

            //these needed for mvPred
            {
                miPtr[miX + miY * mi_stride].mbmi.mode = cu_ptr->pred_mode;

                if (cu_ptr->prediction_mode_flag == INTRA_MODE && cu_ptr->pred_mode == INTRA_MODE_4x4) {
                    miPtr[miX + miY * mi_stride].mbmi.tx_size = 0;
                    miPtr[miX + miY * mi_stride].mbmi.sb_type = BLOCK_4X4;
                }
                else {
                    int32_t txb_itr;
                    for (txb_itr = 0; txb_itr < blk_geom->txb_count; txb_itr++) {
                        miPtr[miX + miY * mi_stride].mbmi.tx_size = blk_geom->txsize[txb_itr]; // Nader - TO_DO
                    }
                    miPtr[miX + miY * mi_stride].mbmi.sb_type = blk_geom->bsize;


                }

                miPtr[miX + miY * mi_stride].mbmi.ref_frame[0] = rf[0];
                miPtr[miX + miY * mi_stride].mbmi.ref_frame[1] = rf[1];
                if (cu_ptr->prediction_unit_array->inter_pred_direction_index == UNI_PRED_LIST_0) {
                    miPtr[miX + miY * mi_stride].mbmi.mv[0].as_mv.col = cu_ptr->prediction_unit_array->mv[0].x;
                    miPtr[miX + miY * mi_stride].mbmi.mv[0].as_mv.row = cu_ptr->prediction_unit_array->mv[0].y;
                }
                else if (cu_ptr->prediction_unit_array->inter_pred_direction_index == UNI_PRED_LIST_1) {
                    miPtr[miX + miY * mi_stride].mbmi.mv[0].as_mv.col = cu_ptr->prediction_unit_array->mv[1].x;
                    miPtr[miX + miY * mi_stride].mbmi.mv[0].as_mv.row = cu_ptr->prediction_unit_array->mv[1].y;
                }
                else {
                    miPtr[miX + miY * mi_stride].mbmi.mv[0].as_mv.col = cu_ptr->prediction_unit_array->mv[0].x;
                    miPtr[miX + miY * mi_stride].mbmi.mv[0].as_mv.row = cu_ptr->prediction_unit_array->mv[0].y;
                    miPtr[miX + miY * mi_stride].mbmi.mv[1].as_mv.col = cu_ptr->prediction_unit_array->mv[1].x;
                    miPtr[miX + miY * mi_stride].mbmi.mv[1].as_mv.row = cu_ptr->prediction_unit_array->mv[1].y;
                }

                miPtr[miX + miY * mi_stride].mbmi.partition = from_shape_to_part[blk_geom->shape];// cu_ptr->part;
            }


            //needed for CDEF
            miPtr[miX + miY * mi_stride].mbmi.skip = cu_ptr->block_has_coeff ? EB_FALSE : EB_TRUE;
        }
    }
}


void update_mi_map(
    CodingUnit_t                   *cu_ptr,
    uint32_t                        cu_origin_x,
    uint32_t                        cu_origin_y,
    const BlockGeom                *blk_geom,
    const CodedUnitStats_t         *cu_stats,
    PictureControlSet_t            *picture_control_set_ptr)
{
    UNUSED(cu_stats);
    uint32_t mi_stride = picture_control_set_ptr->parent_pcs_ptr->sequence_control_set_ptr->picture_width_in_sb*(BLOCK_SIZE_64 >> MI_SIZE_LOG2);
    int32_t mi_row = cu_origin_y >> MI_SIZE_LOG2;
    int32_t mi_col = cu_origin_x >> MI_SIZE_LOG2;

    const int32_t offset = mi_row * mi_stride + mi_col;
    ModeInfo *miPtr = *(picture_control_set_ptr->mi_grid_base + offset);
    MvReferenceFrame rf[2];
    av1_set_ref_frame(rf, cu_ptr->prediction_unit_array->ref_frame_type);

    uint8_t  miX, miY;
    for (miY = 0; miY < (blk_geom->bheight >> MI_SIZE_LOG2); miY++) {

        for (miX = 0; miX < (blk_geom->bwidth >> MI_SIZE_LOG2); miX++) {

            //these needed for mvPred
            {
                miPtr[miX + miY * mi_stride].mbmi.mode = cu_ptr->pred_mode;

                if (cu_ptr->prediction_mode_flag == INTRA_MODE && cu_ptr->pred_mode == INTRA_MODE_4x4) {
                    miPtr[miX + miY * mi_stride].mbmi.tx_size = 0;
                    miPtr[miX + miY * mi_stride].mbmi.sb_type = BLOCK_4X4;
                }
                else {

                    int32_t txb_itr;
                    for (txb_itr = 0; txb_itr < blk_geom->txb_count; txb_itr++) {
                        miPtr[miX + miY * mi_stride].mbmi.tx_size = blk_geom->txsize[txb_itr]; // Nader - TO_DO
                    }

                    miPtr[miX + miY * mi_stride].mbmi.sb_type = blk_geom->bsize;
                }

                miPtr[miX + miY * mi_stride].mbmi.ref_frame[0] = rf[0];
                miPtr[miX + miY * mi_stride].mbmi.ref_frame[1] = rf[1];
                if (cu_ptr->prediction_unit_array->inter_pred_direction_index == UNI_PRED_LIST_0) {
                    miPtr[miX + miY * mi_stride].mbmi.mv[0].as_mv.col = cu_ptr->prediction_unit_array->mv[0].x;
                    miPtr[miX + miY * mi_stride].mbmi.mv[0].as_mv.row = cu_ptr->prediction_unit_array->mv[0].y;
                }
                else if (cu_ptr->prediction_unit_array->inter_pred_direction_index == UNI_PRED_LIST_1) {
                    miPtr[miX + miY * mi_stride].mbmi.mv[0].as_mv.col = cu_ptr->prediction_unit_array->mv[1].x;
                    miPtr[miX + miY * mi_stride].mbmi.mv[0].as_mv.row = cu_ptr->prediction_unit_array->mv[1].y;
                }
                else {
                    miPtr[miX + miY * mi_stride].mbmi.mv[0].as_mv.col = cu_ptr->prediction_unit_array->mv[0].x;
                    miPtr[miX + miY * mi_stride].mbmi.mv[0].as_mv.row = cu_ptr->prediction_unit_array->mv[0].y;
                    miPtr[miX + miY * mi_stride].mbmi.mv[1].as_mv.col = cu_ptr->prediction_unit_array->mv[1].x;
                    miPtr[miX + miY * mi_stride].mbmi.mv[1].as_mv.row = cu_ptr->prediction_unit_array->mv[1].y;
                }

                miPtr[miX + miY * mi_stride].mbmi.partition = from_shape_to_part[blk_geom->shape];// cu_ptr->part;
            }

            if (blk_geom->has_uv)
                miPtr[miX + miY * mi_stride].mbmi.skip = (cu_ptr->transform_unit_array[0].y_has_coeff == 0 && cu_ptr->transform_unit_array[0].v_has_coeff == 0 && cu_ptr->transform_unit_array[0].u_has_coeff == 0) ? EB_TRUE : EB_FALSE;
            else
                miPtr[miX + miY * mi_stride].mbmi.skip = (cu_ptr->transform_unit_array[0].y_has_coeff == 0) ? EB_TRUE : EB_FALSE;

        }
    }
}


static INLINE void record_samples(
    MbModeInfo *mbmi,
    int *pts,
    int *pts_inref,
    int row_offset,
    int sign_r,
    int col_offset,
    int sign_c)
{
    uint8_t bw = block_size_wide[mbmi->sb_type];
    uint8_t bh = block_size_high[mbmi->sb_type];
    int x = col_offset * MI_SIZE + sign_c * AOMMAX(bw, MI_SIZE) / 2 - 1;
    int y = row_offset * MI_SIZE + sign_r * AOMMAX(bh, MI_SIZE) / 2 - 1;

    pts[0] = (x * 8);
    pts[1] = (y * 8);
    pts_inref[0] = (x * 8) + mbmi->mv[0].as_mv.col;
    pts_inref[1] = (y * 8) + mbmi->mv[0].as_mv.row;
}


// Select samples according to the motion vector difference.
int select_samples(
    MV *mv,
    int *pts,
    int *pts_inref,
    int len,
    BlockSize bsize)
{
  const uint8_t bw = block_size_wide[bsize];
  const uint8_t bh = block_size_high[bsize];
  const int thresh = clamp(AOMMAX(bw, bh), 16, 112);
  int pts_mvd[SAMPLES_ARRAY_SIZE] = { 0 };
  int i, j, k, l = len;
  int ret = 0;
  // assert(len <= LEAST_SQUARES_SAMPLES_MAX);

  // Obtain the motion vector difference.
  for (i = 0; i < len; ++i) {
    pts_mvd[i] = abs(pts_inref[2 * i] - pts[2 * i] - mv->col) +
                 abs(pts_inref[2 * i + 1] - pts[2 * i + 1] - mv->row);

    if (pts_mvd[i] > thresh)
      pts_mvd[i] = -1;
    else
      ret++;
  }

  // Keep at least 1 sample.
  if (!ret)
    return 1;

  i = 0;
  j = l - 1;
  for (k = 0; k < l - ret; k++) {
    while (pts_mvd[i] != -1) i++;
    while (pts_mvd[j] == -1) j--;
    assert(i != j);
    if (i > j)
        break;

    // Replace the discarded samples;
    pts_mvd[i] = pts_mvd[j];
    pts[2 * i] = pts[2 * j];
    pts[2 * i + 1] = pts[2 * j + 1];
    pts_inref[2 * i] = pts_inref[2 * j];
    pts_inref[2 * i + 1] = pts_inref[2 * j + 1];
    i++;
    j--;
  }

  return ret;
}


// Note: Samples returned are at 1/8-pel precision
// Sample are the neighbor block center point's coordinates relative to the
// left-top pixel of current block.
int av1_find_samples(
    const Av1Common *cm,
    MacroBlockD *xd,
    int mi_row,
    int mi_col,
    int *pts,
    int *pts_inref)
{
  MbModeInfo *const mbmi0 = &xd->mi[0]->mbmi;
  int ref_frame = mbmi0->ref_frame[0];
  int up_available = xd->up_available;
  int left_available = xd->left_available;
  int i, mi_step = 1, np = 0;

  const TileInfo *const tile = &xd->tile;
  int do_tl = 1;
  int do_tr = 1;

  // scan the nearest above rows
  if (up_available) {
    int mi_row_offset = -1;
    MbModeInfo *mbmi = &xd->mi[mi_row_offset * xd->mi_stride]->mbmi;
    uint8_t n4_w = mi_size_wide[mbmi->sb_type];

    if (xd->n4_w <= n4_w) {
      // Handle "current block width <= above block width" case.
      int col_offset = -mi_col % n4_w;

      if (col_offset < 0) do_tl = 0;
      if (col_offset + n4_w > xd->n4_w) do_tr = 0;

      if (mbmi->ref_frame[0] == ref_frame && mbmi->ref_frame[1] == NONE_FRAME) {
        record_samples(mbmi, pts, pts_inref, 0, -1, col_offset, 1);
        pts += 2;
        pts_inref += 2;
        np++;
        if (np >= LEAST_SQUARES_SAMPLES_MAX)
            return LEAST_SQUARES_SAMPLES_MAX;
      }
    } else {
      // Handle "current block width > above block width" case.
      for (i = 0; i < AOMMIN(xd->n4_w, cm->mi_cols - mi_col); i += mi_step) {
        int mi_col_offset = i;
        mbmi = &xd->mi[mi_col_offset + mi_row_offset * xd->mi_stride]->mbmi;
        n4_w = mi_size_wide[mbmi->sb_type];
        mi_step = AOMMIN(xd->n4_w, n4_w);

        if (mbmi->ref_frame[0] == ref_frame &&
            mbmi->ref_frame[1] == NONE_FRAME)
        {
          record_samples(mbmi, pts, pts_inref, 0, -1, i, 1);
          pts += 2;
          pts_inref += 2;
          np++;
          if (np >= LEAST_SQUARES_SAMPLES_MAX)
            return LEAST_SQUARES_SAMPLES_MAX;
        }
      }
    }
  }

  // scan the nearest left columns
  if (left_available) {
    int mi_col_offset = -1;

    MbModeInfo *mbmi = &xd->mi[mi_col_offset]->mbmi;
    uint8_t n4_h = mi_size_high[mbmi->sb_type];

    if (xd->n4_h <= n4_h) {
      // Handle "current block height <= above block height" case.
      int row_offset = -mi_row % n4_h;

      if (row_offset < 0) do_tl = 0;

      if (mbmi->ref_frame[0] == ref_frame && mbmi->ref_frame[1] == NONE_FRAME) {
        record_samples(mbmi, pts, pts_inref, row_offset, 1, 0, -1);
        pts += 2;
        pts_inref += 2;
        np++;
        if (np >= LEAST_SQUARES_SAMPLES_MAX)
            return LEAST_SQUARES_SAMPLES_MAX;
      }
    } else {
      // Handle "current block height > above block height" case.
      for (i = 0; i < AOMMIN(xd->n4_h, cm->mi_rows - mi_row); i += mi_step) {
        int mi_row_offset = i;
        mbmi = &xd->mi[mi_col_offset + mi_row_offset * xd->mi_stride]->mbmi;
        n4_h = mi_size_high[mbmi->sb_type];
        mi_step = AOMMIN(xd->n4_h, n4_h);

        if (mbmi->ref_frame[0] == ref_frame &&
            mbmi->ref_frame[1] == NONE_FRAME)
        {
          record_samples(mbmi, pts, pts_inref, i, 1, 0, -1);
          pts += 2;
          pts_inref += 2;
          np++;
          if (np >= LEAST_SQUARES_SAMPLES_MAX)
            return LEAST_SQUARES_SAMPLES_MAX;
        }
      }
    }
  }

  // Top-left block
  if (do_tl && left_available && up_available) {
    int mi_row_offset = -1;
    int mi_col_offset = -1;

    MbModeInfo *mbmi = &xd->mi[mi_col_offset + mi_row_offset * xd->mi_stride]->mbmi;

    if (mbmi->ref_frame[0] == ref_frame && mbmi->ref_frame[1] == NONE_FRAME) {
      record_samples(mbmi, pts, pts_inref, 0, -1, 0, -1);
      pts += 2;
      pts_inref += 2;
      np++;
      if (np >= LEAST_SQUARES_SAMPLES_MAX)
        return LEAST_SQUARES_SAMPLES_MAX;
    }
  }

  // Top-right block
  if (do_tr &&
      has_top_right(cm, xd, mi_row, mi_col, AOMMAX(xd->n4_w, xd->n4_h)))
  {
    Position trb_pos = { -1, xd->n4_w };

    if (is_inside(tile, mi_col, mi_row, cm->mi_rows, &trb_pos)) {
      int mi_row_offset = -1;
      int mi_col_offset = xd->n4_w;

      MbModeInfo *mbmi =
          &xd->mi[mi_col_offset + mi_row_offset * xd->mi_stride]->mbmi;

      if (mbmi->ref_frame[0] == ref_frame && mbmi->ref_frame[1] == NONE_FRAME) {
        record_samples(mbmi, pts, pts_inref, 0, -1, xd->n4_w, 1);
        np++;
        if (np >= LEAST_SQUARES_SAMPLES_MAX)
            return LEAST_SQUARES_SAMPLES_MAX;
      }
    }
  }

  return np;
}


uint16_t wm_find_samples(
    CodingUnit_t                       *cu_ptr,
    const BlockGeom                    *blk_geom,
    uint16_t                            pu_origin_x,
    uint16_t                            pu_origin_y,
    PictureControlSet_t                *picture_control_set_ptr,
    int32_t                            *pts,
    int32_t                            *pts_inref)
{
    Av1Common  *cm = picture_control_set_ptr->parent_pcs_ptr->av1_cm;
    MacroBlockD  *xd = cu_ptr->av1xd;

    int32_t mi_row = pu_origin_y >> MI_SIZE_LOG2;
    int32_t mi_col = pu_origin_x >> MI_SIZE_LOG2;

    xd->n4_w = blk_geom->bwidth  >> MI_SIZE_LOG2;
    xd->n4_h = blk_geom->bheight >> MI_SIZE_LOG2;

    return (uint16_t) av1_find_samples(cm, xd, mi_row, mi_col, pts, pts_inref);
}


EbBool warped_motion_parameters(
    PictureControlSet_t              *picture_control_set_ptr,
    CodingUnit_t                     *cu_ptr,
    MvUnit_t                         *mv_unit,
    const BlockGeom                  *blk_geom,
    uint16_t                          pu_origin_x,
    uint16_t                          pu_origin_y,
    EbWarpedMotionParams             *wm_params,
    uint16_t                         *num_samples)
{
    MacroBlockD  *xd = cu_ptr->av1xd;
    BlockSize bsize = blk_geom->bsize;
    EbBool apply_wm = EB_FALSE;

    int pts[SAMPLES_ARRAY_SIZE], pts_inref[SAMPLES_ARRAY_SIZE];
    int32_t mi_row = pu_origin_y >> MI_SIZE_LOG2;
    int32_t mi_col = pu_origin_x >> MI_SIZE_LOG2;
    xd->n4_w = blk_geom->bwidth  >> MI_SIZE_LOG2;
    xd->n4_h = blk_geom->bheight >> MI_SIZE_LOG2;

    *num_samples = 0;
    if (blk_geom->bwidth < 8 || blk_geom->bheight < 8)
        return apply_wm;

    uint16_t nsamples = wm_find_samples(
        cu_ptr,
        blk_geom,
        pu_origin_x,
        pu_origin_y,
        picture_control_set_ptr,
        pts,
        pts_inref);

    if(nsamples==0)
        return apply_wm;

    MV mv;
    mv.col = mv_unit->mv[REF_LIST_0].x;
    mv.row = mv_unit->mv[REF_LIST_0].y;

    if(nsamples > 1)
        nsamples = select_samples(&mv, pts, pts_inref, nsamples, bsize);
    *num_samples = nsamples;

    apply_wm = !find_projection(
        (int) nsamples,
        pts,
        pts_inref,
        bsize,
        mv.row,
        mv.col,
        wm_params,
        (int) mi_row,
        (int) mi_col);

    return apply_wm;
}

