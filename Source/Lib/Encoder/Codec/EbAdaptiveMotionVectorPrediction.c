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
#include "EbAdaptiveMotionVectorPrediction.h"
#include "EbSvtAv1.h"
#include "EbModeDecisionProcess.h"
#include "EbCommonUtils.h"
#include "EbEntropyCoding.h"

#define UNUSED_FUNC

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
static PartitionType from_shape_to_part[] = {PARTITION_NONE,
                                             PARTITION_HORZ,
                                             PARTITION_VERT,
                                             PARTITION_HORZ_A,
                                             PARTITION_HORZ_B,
                                             PARTITION_VERT_A,
                                             PARTITION_VERT_B,
                                             PARTITION_HORZ_4,
                                             PARTITION_VERT_4,
                                             PARTITION_SPLIT};

#define MVREF_ROWS 3
#define MVREF_COLS 3

static int32_t have_newmv_in_inter_mode(PredictionMode mode) {
    return (mode == NEWMV || mode == NEW_NEWMV || mode == NEAREST_NEWMV || mode == NEW_NEARESTMV ||
            mode == NEAR_NEWMV || mode == NEW_NEARMV);
}

typedef struct position {
    int32_t row;
    int32_t col;
} Position;

// clang-format on

static INLINE IntMv get_sub_block_mv(const ModeInfo *candidate, int32_t which_mv,
                                     int32_t search_col) {
    (void)search_col;
    return candidate->mbmi.block_mi.mv[which_mv];
}
static INLINE int32_t is_inside(const TileInfo *const tile, int32_t mi_col, int32_t mi_row,
                                int32_t mi_rows, const Position *mi_pos) {
    const int32_t dependent_horz_tile_flag = 0;
    if (dependent_horz_tile_flag && !tile->tg_horz_boundary) {
        return !(mi_row + mi_pos->row < 0 || mi_col + mi_pos->col < tile->mi_col_start ||
                 mi_row + mi_pos->row >= mi_rows || mi_col + mi_pos->col >= tile->mi_col_end);
    } else {
        return !(mi_row + mi_pos->row < tile->mi_row_start ||
                 mi_col + mi_pos->col < tile->mi_col_start ||
                 mi_row + mi_pos->row >= tile->mi_row_end ||
                 mi_col + mi_pos->col >= tile->mi_col_end);
    }
}

static INLINE void clamp_mv_ref(MV *mv, int32_t bw, int32_t bh, const MacroBlockD *xd) {
    clamp_mv(mv,
             xd->mb_to_left_edge - bw * 8 - MV_BORDER,
             xd->mb_to_right_edge + bw * 8 + MV_BORDER,
             xd->mb_to_top_edge - bh * 8 - MV_BORDER,
             xd->mb_to_bottom_edge + bh * 8 + MV_BORDER);
}

static void add_ref_mv_candidate(const ModeInfo *const   candidate_mi,
                                 const MbModeInfo *const candidate, const MvReferenceFrame rf[2],
                                 uint8_t     refmv_counts[MODE_CTX_REF_FRAMES],
                                 uint8_t     ref_match_counts[MODE_CTX_REF_FRAMES],
                                 uint8_t     newmv_counts[MODE_CTX_REF_FRAMES],
                                 CandidateMv ref_mv_stacks[][MAX_REF_MV_STACK_SIZE], int32_t len,
                                 IntMv *gm_mv_candidates, const EbWarpedMotionParams *gm_params,
                                 int32_t col, int32_t weight) {
    if (!is_inter_block(&candidate->block_mi)) return; // for intrabc
    int32_t index = 0, ref;
    assert(weight % 2 == 0);

    if (rf[1] == NONE_FRAME) {
        uint8_t *    refmv_count     = &refmv_counts[rf[0]];
        uint8_t *    ref_match_count = &ref_match_counts[rf[0]];
        uint8_t *    newmv_count     = &newmv_counts[rf[0]];
        CandidateMv *ref_mv_stack    = ref_mv_stacks[rf[0]];
        (void)ref_match_count;

        // single reference frame
        for (ref = 0; ref < 2; ++ref) {
            if (candidate->block_mi.ref_frame[ref] == rf[0]) {
                IntMv this_refmv;
                if (is_global_mv_block(candidate->block_mi.mode,
                                       candidate->block_mi.sb_type,
                                       gm_params[rf[0]].wmtype))
                    this_refmv = gm_mv_candidates[0];
                else
                    this_refmv = get_sub_block_mv(candidate_mi, ref, col);

                for (index = 0; index < *refmv_count; ++index)
                    if (ref_mv_stack[index].this_mv.as_int == this_refmv.as_int) break;

                if (index < *refmv_count) ref_mv_stack[index].weight += weight * len;

                // Add a new item to the list.
                if (index == *refmv_count && *refmv_count < MAX_REF_MV_STACK_SIZE) {
                    ref_mv_stack[index].this_mv = this_refmv;
                    ref_mv_stack[index].weight  = weight * len;
                    ++(*refmv_count);
                }
                if (have_newmv_in_inter_mode(candidate->block_mi.mode)) ++*newmv_count;
                ++*ref_match_count;
            }
        }
    } else {
        MvReferenceFrame ref_frame       = av1_ref_frame_type(rf);
        uint8_t *        refmv_count     = &refmv_counts[ref_frame];
        uint8_t *        ref_match_count = &ref_match_counts[ref_frame];
        uint8_t *        newmv_count     = &newmv_counts[ref_frame];
        CandidateMv *    ref_mv_stack    = ref_mv_stacks[ref_frame];
        (void)ref_match_count;

        // compound reference frame
        if (candidate->block_mi.ref_frame[0] == rf[0] &&
            candidate->block_mi.ref_frame[1] == rf[1]) {
            IntMv this_refmv[2];

            for (ref = 0; ref < 2; ++ref) {
                if (is_global_mv_block(candidate->block_mi.mode,
                                       candidate->block_mi.sb_type,
                                       gm_params[rf[ref]].wmtype))
                    this_refmv[ref] = gm_mv_candidates[ref];
                else
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
                ref_mv_stack[index].weight  = weight * len;
                ++(*refmv_count);
            }
            if (have_newmv_in_inter_mode(candidate->block_mi.mode)) ++*newmv_count;
            ++*ref_match_count;
        }
    }
}

static void scan_row_mbmi(const Av1Common *cm, const MacroBlockD *xd, int32_t mi_row,
                          int32_t mi_col, const MvReferenceFrame rf[2], int32_t row_offset,
                          CandidateMv ref_mv_stack[][MAX_REF_MV_STACK_SIZE],
                          uint8_t     refmv_count[MODE_CTX_REF_FRAMES],
                          uint8_t     ref_match_count[MODE_CTX_REF_FRAMES],
                          uint8_t newmv_count[MODE_CTX_REF_FRAMES], IntMv *gm_mv_candidates,
                          const EbWarpedMotionParams *gm_params, int32_t max_row_offset,
                          int32_t *processed_rows) {
    int32_t end_mi        = AOMMIN(xd->n8_w, cm->mi_cols - mi_col);
    end_mi                = AOMMIN(end_mi, mi_size_wide[BLOCK_64X64]);
    const int32_t n8_w_8  = mi_size_wide[BLOCK_8X8];
    const int32_t n8_w_16 = mi_size_wide[BLOCK_16X16];
    int32_t       i;
    int32_t       col_offset = 0;
    const int32_t shift      = 0;
    if (abs(row_offset) > 1) {
        col_offset = 1;
        if (mi_col & 0x01 && xd->n8_w < n8_w_8) --col_offset;
    }
    const int32_t    use_step_16   = (xd->n8_w >= 16);
    ModeInfo **const candidate_mi0 = xd->mi + row_offset * xd->mi_stride;
    (void)mi_row;

    for (i = 0; i < end_mi;) {
        const ModeInfo *const   candidate_mi    = candidate_mi0[col_offset + i];
        const MbModeInfo *const candidate       = &candidate_mi->mbmi;
        const int32_t           candidate_bsize = candidate->block_mi.sb_type;
        assert(candidate_bsize < BlockSizeS_ALL);
        const int32_t n8_w = mi_size_wide[candidate_bsize];
        int32_t       len  = AOMMIN(xd->n8_w, n8_w);
        if (use_step_16)
            len = AOMMAX(n8_w_16, len);
        else if (abs(row_offset) > 1)
            len = AOMMAX(len, n8_w_8);

        int32_t weight = 2;
        if (xd->n8_w >= n8_w_8 && xd->n8_w <= n8_w) {
            int32_t inc = AOMMIN(-max_row_offset + row_offset + 1, mi_size_high[candidate_bsize]);
            // Obtain range used in weight calculation.
            weight = AOMMAX(weight, (inc << shift));
            // Update processed rows.
            *processed_rows = inc - row_offset - 1;
        }

        add_ref_mv_candidate(candidate_mi,
                             candidate,
                             rf,
                             refmv_count,
                             ref_match_count,
                             newmv_count,
                             ref_mv_stack,
                             len,
                             gm_mv_candidates,
                             gm_params,
                             col_offset + i,
                             weight);

        i += len;
    }
}

static void scan_col_mbmi(const Av1Common *cm, const MacroBlockD *xd, int32_t mi_row,
                          int32_t mi_col, const MvReferenceFrame rf[2], int32_t col_offset,
                          CandidateMv ref_mv_stack[][MAX_REF_MV_STACK_SIZE],
                          uint8_t     refmv_count[MODE_CTX_REF_FRAMES],
                          uint8_t     ref_match_count[MODE_CTX_REF_FRAMES],
                          uint8_t newmv_count[MODE_CTX_REF_FRAMES], IntMv *gm_mv_candidates,
                          const EbWarpedMotionParams *gm_params, int32_t max_col_offset,
                          int32_t *processed_cols) {
    int32_t end_mi        = AOMMIN(xd->n8_h, cm->mi_rows - mi_row);
    end_mi                = AOMMIN(end_mi, mi_size_high[BLOCK_64X64]);
    const int32_t n8_h_8  = mi_size_high[BLOCK_8X8];
    const int32_t n8_h_16 = mi_size_high[BLOCK_16X16];
    int32_t       i;
    int32_t       row_offset = 0;
    const int32_t shift      = 0;
    if (abs(col_offset) > 1) {
        row_offset = 1;
        if (mi_row & 0x01 && xd->n8_h < n8_h_8) --row_offset;
    }
    const int32_t use_step_16 = (xd->n8_h >= 16);
    (void)mi_col;

    for (i = 0; i < end_mi;) {
        const ModeInfo *const candidate_mi = xd->mi[(row_offset + i) * xd->mi_stride + col_offset];
        const MbModeInfo *const candidate  = &candidate_mi->mbmi;
        const int32_t           candidate_bsize = candidate->block_mi.sb_type;
        assert(candidate_bsize < BlockSizeS_ALL);
        const int32_t n8_h = mi_size_high[candidate_bsize];
        int32_t       len  = AOMMIN(xd->n8_h, n8_h);
        if (use_step_16)
            len = AOMMAX(n8_h_16, len);
        else if (abs(col_offset) > 1)
            len = AOMMAX(len, n8_h_8);

        int32_t weight = 2;
        if (xd->n8_h >= n8_h_8 && xd->n8_h <= n8_h) {
            int32_t inc = AOMMIN(-max_col_offset + col_offset + 1, mi_size_wide[candidate_bsize]);
            // Obtain range used in weight calculation.
            weight = AOMMAX(weight, (inc << shift));
            // Update processed cols.
            *processed_cols = inc - col_offset - 1;
        }

        add_ref_mv_candidate(candidate_mi,
                             candidate,
                             rf,
                             refmv_count,
                             ref_match_count,
                             newmv_count,
                             ref_mv_stack,
                             len,
                             gm_mv_candidates,
                             gm_params,
                             col_offset,
                             weight);

        i += len;
    }
}

static void scan_blk_mbmi(const Av1Common *cm, const MacroBlockD *xd, const int32_t mi_row,
                          const int32_t mi_col, const MvReferenceFrame rf[2], int32_t row_offset,
                          int32_t col_offset, CandidateMv ref_mv_stack[][MAX_REF_MV_STACK_SIZE],
                          uint8_t ref_match_count[MODE_CTX_REF_FRAMES],
                          uint8_t newmv_count[MODE_CTX_REF_FRAMES], IntMv *gm_mv_candidates,
                          const EbWarpedMotionParams *gm_params,
                          uint8_t                     refmv_count[MODE_CTX_REF_FRAMES]) {
    const TileInfo *const tile = &xd->tile;
    Position              mi_pos;

    mi_pos.row = row_offset;
    mi_pos.col = col_offset;

    if (is_inside(tile, mi_col, mi_row, cm->mi_rows, &mi_pos)) {
        const ModeInfo *const   candidate_mi = xd->mi[mi_pos.row * xd->mi_stride + mi_pos.col];
        const MbModeInfo *const candidate    = &candidate_mi->mbmi;
        const int32_t           len          = mi_size_wide[BLOCK_8X8];

        add_ref_mv_candidate(candidate_mi,
                             candidate,
                             rf,
                             refmv_count,
                             ref_match_count,
                             newmv_count,
                             ref_mv_stack,
                             len,
                             gm_mv_candidates,
                             gm_params,
                             mi_pos.col,
                             2);
    } // Analyze a single 8x8 block motion information.
}

static int32_t has_top_right(const Av1Common *cm,  const BlockSize sb_size, const MacroBlockD *xd, int32_t mi_row,
                             int32_t mi_col, int32_t bs) {
    (void)xd;
    (void)cm;
    const int32_t sb_mi_size = mi_size_wide[sb_size];
    const int32_t mask_row   = mi_row & (sb_mi_size - 1);
    const int32_t mask_col   = mi_col & (sb_mi_size - 1);

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
        } else
            break;
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
    if (xd->mi[0]->mbmi.block_mi.partition == PARTITION_VERT_A) {
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
        return clamp(row_offset, tile->mi_row_start - mi_row, tile->mi_row_end - mi_row - 1);
}

static INLINE int32_t find_valid_col_offset(const TileInfo *const tile, int32_t mi_col,
                                            int32_t col_offset) {
    return clamp(col_offset, tile->mi_col_start - mi_col, tile->mi_col_end - mi_col - 1);
}
static INLINE int get_relative_dist(const OrderHintInfo *oh, int a, int b) {
    if (!oh->enable_order_hint) return 0;

    const int bits = oh->order_hint_bits;

    assert(bits >= 1);
    assert(a >= 0 && a < (1 << bits));
    assert(b >= 0 && b < (1 << bits));

    int       diff = a - b;
    const int m    = 1 << (bits - 1);
    diff           = (diff & (m - 1)) - (diff & m);
    return diff;
}
static int add_tpl_ref_mv(const Av1Common *cm, PictureControlSet *pcs_ptr, const MacroBlockD *xd,
                          int mi_row, int mi_col, MvReferenceFrame ref_frame, int blk_row,
                          int blk_col, IntMv *gm_mv_candidates, uint8_t *const refmv_count,
                          CandidateMv ref_mv_stack[MAX_REF_MV_STACK_SIZE], int16_t *mode_context) {
    Position mi_pos;
    mi_pos.row = (mi_row & 0x01) ? blk_row : blk_row + 1;
    mi_pos.col = (mi_col & 0x01) ? blk_col : blk_col + 1;

    if (!is_inside(&xd->tile, mi_col, mi_row, xd->tile.mi_row_end, &mi_pos)) return 0;

    const TPL_MV_REF *prev_frame_mvs = pcs_ptr->tpl_mvs +
                                       ((mi_row + mi_pos.row) >> 1) * (cm->mi_stride >> 1) +
                                       ((mi_col + mi_pos.col) >> 1);
    if (prev_frame_mvs->mfmv0.as_int == INVALID_MV) return 0;

    MvReferenceFrame rf[2];
    av1_set_ref_frame(rf, ref_frame);

    uint8_t list_idx0, list_idx1, ref_idx_l0, ref_idx_l1;
    list_idx0  = get_list_idx(rf[0]);
    ref_idx_l0 = get_ref_frame_idx(rf[0]);
    if (rf[1] == NONE_FRAME) {
        list_idx1  = get_list_idx(rf[0]);
        ref_idx_l1 = get_ref_frame_idx(rf[0]);
    } else {
        list_idx1  = get_list_idx(rf[1]);
        ref_idx_l1 = get_ref_frame_idx(rf[1]);
    }

    const uint16_t     weight_unit     = 1;
    const int          cur_frame_index = pcs_ptr->parent_pcs_ptr->cur_order_hint;
    EbReferenceObject *buf_0 =
        (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[list_idx0][ref_idx_l0]->object_ptr;

    const int frame0_index = buf_0->order_hint;
    const int cur_offset_0 =
        get_relative_dist(&pcs_ptr->parent_pcs_ptr->scs_ptr->seq_header.order_hint_info,
                          cur_frame_index,
                          frame0_index);
    int idx;

    IntMv this_refmv;
    get_mv_projection(&this_refmv.as_mv,
                      prev_frame_mvs->mfmv0.as_mv,
                      cur_offset_0,
                      prev_frame_mvs->ref_frame_offset);
    lower_mv_precision(
        &this_refmv.as_mv, pcs_ptr->parent_pcs_ptr->frm_hdr.allow_high_precision_mv, 0);
    if (rf[1] == NONE_FRAME) {
        if (blk_row == 0 && blk_col == 0) {
            if (abs(this_refmv.as_mv.row - gm_mv_candidates[0].as_mv.row) >= 16 ||
                abs(this_refmv.as_mv.col - gm_mv_candidates[0].as_mv.col) >= 16)
                mode_context[ref_frame] |= (1 << GLOBALMV_OFFSET);
        }

        for (idx = 0; idx < *refmv_count; ++idx)
            if (this_refmv.as_int == ref_mv_stack[idx].this_mv.as_int) break;

        if (idx < *refmv_count) ref_mv_stack[idx].weight += 2 * weight_unit;

        if (idx == *refmv_count && *refmv_count < MAX_REF_MV_STACK_SIZE) {
            ref_mv_stack[idx].this_mv.as_int = this_refmv.as_int;
            ref_mv_stack[idx].weight         = 2 * weight_unit;
            ++(*refmv_count);
        }
    } else {
        // Process compound inter mode
        EbReferenceObject *buf_1 =
            (EbReferenceObject *)pcs_ptr->ref_pic_ptr_array[list_idx1][ref_idx_l1]->object_ptr;

        const int frame1_index = buf_1->order_hint;
        const int cur_offset_1 =
            get_relative_dist(&pcs_ptr->parent_pcs_ptr->scs_ptr->seq_header.order_hint_info,
                              cur_frame_index,
                              frame1_index);
        IntMv comp_refmv;
        get_mv_projection(&comp_refmv.as_mv,
                          prev_frame_mvs->mfmv0.as_mv,
                          cur_offset_1,
                          prev_frame_mvs->ref_frame_offset);
        lower_mv_precision(
            &comp_refmv.as_mv, pcs_ptr->parent_pcs_ptr->frm_hdr.allow_high_precision_mv, 0);
        if (blk_row == 0 && blk_col == 0) {
            if (abs(this_refmv.as_mv.row - gm_mv_candidates[0].as_mv.row) >= 16 ||
                abs(this_refmv.as_mv.col - gm_mv_candidates[0].as_mv.col) >= 16 ||
                abs(comp_refmv.as_mv.row - gm_mv_candidates[1].as_mv.row) >= 16 ||
                abs(comp_refmv.as_mv.col - gm_mv_candidates[1].as_mv.col) >= 16)
                mode_context[ref_frame] |= (1 << GLOBALMV_OFFSET);
        }

        for (idx = 0; idx < *refmv_count; ++idx) {
            if (this_refmv.as_int == ref_mv_stack[idx].this_mv.as_int &&
                comp_refmv.as_int == ref_mv_stack[idx].comp_mv.as_int)
                break;
        }

        if (idx < *refmv_count) ref_mv_stack[idx].weight += 2 * weight_unit;

        if (idx == *refmv_count && *refmv_count < MAX_REF_MV_STACK_SIZE) {
            ref_mv_stack[idx].this_mv.as_int = this_refmv.as_int;
            ref_mv_stack[idx].comp_mv.as_int = comp_refmv.as_int;
            ref_mv_stack[idx].weight         = 2 * weight_unit;
            ++(*refmv_count);
        }
    }

    return 1;
}

void setup_ref_mv_list(PictureControlSet *pcs_ptr, const Av1Common *cm, const MacroBlockD *xd,
                       MvReferenceFrame ref_frame, uint8_t refmv_count[MODE_CTX_REF_FRAMES],
                       CandidateMv ref_mv_stack[][MAX_REF_MV_STACK_SIZE],
                       IntMv mv_ref_list[][MAX_MV_REF_CANDIDATES], IntMv *gm_mv_candidates,
                       const EbWarpedMotionParams *gm_params, int32_t mi_row, int32_t mi_col,
                       int16_t *mode_context) {
    const int32_t    bs     = AOMMAX(xd->n8_w, xd->n8_h);
    const int32_t    has_tr = has_top_right(cm,
            ((SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr)->seq_header.sb_size, xd, mi_row, mi_col, bs);
    MvReferenceFrame rf[2];

    const TileInfo *const tile           = &xd->tile;
    int32_t               max_row_offset = 0, max_col_offset = 0;
    const int32_t         row_adj        = (xd->n8_h < mi_size_high[BLOCK_8X8]) && (mi_row & 0x01);
    const int32_t         col_adj        = (xd->n8_w < mi_size_wide[BLOCK_8X8]) && (mi_col & 0x01);
    int32_t               processed_rows = 0;
    int32_t               processed_cols = 0;

    av1_set_ref_frame(rf, ref_frame);
    mode_context[ref_frame] = 0;
    refmv_count[ref_frame]  = 0;

    // Find valid maximum row/col offset.
    if (xd->up_available) {
        max_row_offset = -(MVREF_ROWS << 1) + row_adj;

        if (xd->n8_h < mi_size_high[BLOCK_8X8]) max_row_offset = -(2 << 1) + row_adj;

        max_row_offset = find_valid_row_offset(tile, mi_row, cm->mi_rows, max_row_offset);
    }

    if (xd->left_available) {
        max_col_offset = -(MVREF_COLS << 1) + col_adj;

        if (xd->n8_w < mi_size_wide[BLOCK_8X8]) max_col_offset = -(2 << 1) + col_adj;

        max_col_offset = find_valid_col_offset(tile, mi_col, max_col_offset);
    }

    uint8_t ref_match_count[MODE_CTX_REF_FRAMES] = {0};
    uint8_t col_match_count[MODE_CTX_REF_FRAMES] = {0};
    uint8_t row_match_count[MODE_CTX_REF_FRAMES] = {0};
    uint8_t newmv_count[MODE_CTX_REF_FRAMES]     = {0};

    //CHKN-------------    ROW-1

    // Scan the first above row mode info. row_offset = -1;
    if (abs(max_row_offset) >= 1)
        scan_row_mbmi(cm,
                      xd,
                      mi_row,
                      mi_col,
                      rf,
                      -1,
                      ref_mv_stack,
                      refmv_count,
                      row_match_count,
                      newmv_count,
                      gm_mv_candidates,
                      gm_params,
                      max_row_offset,
                      &processed_rows);

    //CHKN-------------    COL-1
    // Scan the first left column mode info. col_offset = -1;
    if (abs(max_col_offset) >= 1)
        scan_col_mbmi(cm,
                      xd,
                      mi_row,
                      mi_col,
                      rf,
                      -1,
                      ref_mv_stack,
                      refmv_count,
                      col_match_count,
                      newmv_count,
                      gm_mv_candidates,
                      gm_params,
                      max_col_offset,
                      &processed_cols);

    //CHKN-------------    TOP-RIGHT

    // Check top-right boundary
    if (has_tr)
        scan_blk_mbmi(cm,
                      xd,
                      mi_row,
                      mi_col,
                      rf,
                      -1,
                      xd->n8_w,
                      ref_mv_stack,
                      row_match_count,
                      newmv_count,
                      gm_mv_candidates,
                      gm_params,
                      refmv_count);

    uint8_t nearest_match[MODE_CTX_REF_FRAMES];
    uint8_t nearest_refmv_count[MODE_CTX_REF_FRAMES];

    nearest_match[ref_frame] = (row_match_count[ref_frame] > 0) + (col_match_count[ref_frame] > 0);
    nearest_refmv_count[ref_frame] = refmv_count[ref_frame];

    for (int32_t idx = 0; idx < nearest_refmv_count[ref_frame]; ++idx)
        ref_mv_stack[ref_frame][idx].weight += REF_CAT_LEVEL;

    //CHKN  MFMV - get canididates from reference frames- orderHint has to be on, in order to scale the vectors.
    if (pcs_ptr->parent_pcs_ptr->frm_hdr.use_ref_frame_mvs) {
        int       is_available = 0;
        const int voffset      = AOMMAX(mi_size_high[BLOCK_8X8], xd->n4_h);
        const int hoffset      = AOMMAX(mi_size_wide[BLOCK_8X8], xd->n4_w);
        const int blk_row_end  = AOMMIN(xd->n4_h, mi_size_high[BLOCK_64X64]);
        const int blk_col_end  = AOMMIN(xd->n4_w, mi_size_wide[BLOCK_64X64]);

        const int tpl_sample_pos[3][2] = {
            {voffset, -2},
            {voffset, hoffset},
            {voffset - 2, hoffset},
        };
        const int allow_extension =
            (xd->n4_h >= mi_size_high[BLOCK_8X8]) && (xd->n4_h < mi_size_high[BLOCK_64X64]) &&
            (xd->n4_w >= mi_size_wide[BLOCK_8X8]) && (xd->n4_w < mi_size_wide[BLOCK_64X64]);

        const int step_h = (xd->n4_h >= mi_size_high[BLOCK_64X64]) ? mi_size_high[BLOCK_16X16]
                                                                   : mi_size_high[BLOCK_8X8];
        const int step_w = (xd->n4_w >= mi_size_wide[BLOCK_64X64]) ? mi_size_wide[BLOCK_16X16]
                                                                   : mi_size_wide[BLOCK_8X8];

        for (int blk_row = 0; blk_row < blk_row_end; blk_row += step_h) {
            for (int blk_col = 0; blk_col < blk_col_end; blk_col += step_w) {
                int ret = add_tpl_ref_mv(cm,
                                         pcs_ptr,
                                         xd,
                                         mi_row,
                                         mi_col,
                                         ref_frame,
                                         blk_row,
                                         blk_col,
                                         gm_mv_candidates,
                                         &refmv_count[ref_frame],
                                         ref_mv_stack[ref_frame],
                                         mode_context);
                if (blk_row == 0 && blk_col == 0) is_available = ret;
            }
        }

        if (is_available == 0) mode_context[ref_frame] |= (1 << GLOBALMV_OFFSET);

        for (int i = 0; i < 3 && allow_extension; ++i) {
            const int blk_row = tpl_sample_pos[i][0];
            const int blk_col = tpl_sample_pos[i][1];

            if (!check_sb_border(mi_row, mi_col, blk_row, blk_col)) continue;
            add_tpl_ref_mv(cm,
                           pcs_ptr,
                           xd,
                           mi_row,
                           mi_col,
                           ref_frame,
                           blk_row,
                           blk_col,
                           gm_mv_candidates,
                           &refmv_count[ref_frame],
                           ref_mv_stack[ref_frame],
                           mode_context);
        }
    }

    //CHKN------------- TOP-LEFT
    uint8_t dummy_newmv_count[MODE_CTX_REF_FRAMES] = {0};

    // Scan the second outer area.
    scan_blk_mbmi(cm,
                  xd,
                  mi_row,
                  mi_col,
                  rf,
                  -1,
                  -1,
                  ref_mv_stack,
                  row_match_count,
                  dummy_newmv_count,
                  gm_mv_candidates,
                  gm_params,
                  refmv_count);

    //CHKN-------------    ROW-3  COL-3     ROW-5   COL-5
    for (int32_t idx = 2; idx <= MVREF_ROWS; ++idx) {
        const int32_t row_offset = -(idx << 1) + 1 + row_adj;
        const int32_t col_offset = -(idx << 1) + 1 + col_adj;

        if (abs(row_offset) <= abs(max_row_offset) && abs(row_offset) > processed_rows)
            scan_row_mbmi(cm,
                          xd,
                          mi_row,
                          mi_col,
                          rf,
                          row_offset,
                          ref_mv_stack,
                          refmv_count,
                          row_match_count,
                          dummy_newmv_count,
                          gm_mv_candidates,
                          gm_params,
                          max_row_offset,
                          &processed_rows);

        if (abs(col_offset) <= abs(max_col_offset) && abs(col_offset) > processed_cols)
            scan_col_mbmi(cm,
                          xd,
                          mi_row,
                          mi_col,
                          rf,
                          col_offset,
                          ref_mv_stack,
                          refmv_count,
                          col_match_count,
                          dummy_newmv_count,
                          gm_mv_candidates,
                          gm_params,
                          max_col_offset,
                          &processed_cols);
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
            if (ref_mv_stack[ref_frame][idx - 1].weight < ref_mv_stack[ref_frame][idx].weight) {
                CandidateMv tmp_mv               = ref_mv_stack[ref_frame][idx - 1];
                ref_mv_stack[ref_frame][idx - 1] = ref_mv_stack[ref_frame][idx];
                ref_mv_stack[ref_frame][idx]     = tmp_mv;
                nr_len                           = idx;
            }
        }
        len = nr_len;
    }

    len = refmv_count[ref_frame];
    while (len > nearest_refmv_count[ref_frame]) {
        int32_t nr_len = nearest_refmv_count[ref_frame];
        for (int32_t idx = nearest_refmv_count[ref_frame] + 1; idx < len; ++idx) {
            if (ref_mv_stack[ref_frame][idx - 1].weight < ref_mv_stack[ref_frame][idx].weight) {
                CandidateMv tmp_mv               = ref_mv_stack[ref_frame][idx - 1];
                ref_mv_stack[ref_frame][idx - 1] = ref_mv_stack[ref_frame][idx];
                ref_mv_stack[ref_frame][idx]     = tmp_mv;
                nr_len                           = idx;
            }
        }
        len = nr_len;
    }

    //CHKN finish the Tables.

    if (rf[1] > NONE_FRAME) {

        //CHKN we get here only when refMVCount=0 or 1

        if (refmv_count[ref_frame] < 2) {
            IntMv   ref_id[2][2], ref_diff[2][2];
            int32_t ref_id_count[2] = {0}, ref_diff_count[2] = {0};

            int32_t mi_width  = AOMMIN(mi_size_wide[BLOCK_64X64], xd->n8_w);
            mi_width          = AOMMIN(mi_width, cm->mi_cols - mi_col);
            int32_t mi_height = AOMMIN(mi_size_high[BLOCK_64X64], xd->n8_h);
            mi_height         = AOMMIN(mi_height, cm->mi_rows - mi_row);
            int32_t mi_size   = AOMMIN(mi_width, mi_height);

            //CHKN  scan ROW=-1 again but with more relaxed constraints
            for (int32_t idx = 0; abs(max_row_offset) >= 1 && idx < mi_size;) {
                const ModeInfo *const   candidate_mi    = xd->mi[-xd->mi_stride + idx];
                const MbModeInfo *const candidate       = &candidate_mi->mbmi;
                const int32_t           candidate_bsize = candidate->block_mi.sb_type;

                for (int32_t rf_idx = 0; rf_idx < 2; ++rf_idx) {
                    MvReferenceFrame can_rf = candidate->block_mi.ref_frame[rf_idx];

                    for (int32_t cmp_idx = 0; cmp_idx < 2; ++cmp_idx) {
                        if (can_rf == rf[cmp_idx] && ref_id_count[cmp_idx] < 2) {
                            ref_id[cmp_idx][ref_id_count[cmp_idx]] = candidate->block_mi.mv[rf_idx];
                            ++ref_id_count[cmp_idx];
                        } else if (can_rf > INTRA_FRAME && ref_diff_count[cmp_idx] < 2) {
                            IntMv this_mv = candidate->block_mi.mv[rf_idx];
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
                const ModeInfo *const   candidate_mi    = xd->mi[idx * xd->mi_stride - 1];
                const MbModeInfo *const candidate       = &candidate_mi->mbmi;
                const int32_t           candidate_bsize = candidate->block_mi.sb_type;

                for (int32_t rf_idx = 0; rf_idx < 2; ++rf_idx) {
                    MvReferenceFrame can_rf = candidate->block_mi.ref_frame[rf_idx];

                    for (int32_t cmp_idx = 0; cmp_idx < 2; ++cmp_idx) {
                        if (can_rf == rf[cmp_idx] && ref_id_count[cmp_idx] < 2) {
                            ref_id[cmp_idx][ref_id_count[cmp_idx]] = candidate->block_mi.mv[rf_idx];
                            ++ref_id_count[cmp_idx];
                        } else if (can_rf > INTRA_FRAME && ref_diff_count[cmp_idx] < 2) {
                            IntMv this_mv = candidate->block_mi.mv[rf_idx];
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

            // Build up the compound mv predictor
            IntMv comp_list[MAX_MV_REF_CANDIDATES][2];

            for (int32_t idx = 0; idx < 2; ++idx) {
                int32_t comp_idx = 0;
                for (int32_t list_idx = 0;
                     list_idx < ref_id_count[idx] && comp_idx < MAX_MV_REF_CANDIDATES;
                     ++list_idx, ++comp_idx)
                    comp_list[comp_idx][idx] = ref_id[idx][list_idx];
                for (int32_t list_idx = 0;
                     list_idx < ref_diff_count[idx] && comp_idx < MAX_MV_REF_CANDIDATES;
                     ++list_idx, ++comp_idx)
                    comp_list[comp_idx][idx] = ref_diff[idx][list_idx];
                for (; comp_idx < MAX_MV_REF_CANDIDATES; ++comp_idx)
                    comp_list[comp_idx][idx] = gm_mv_candidates[idx];
            }

            //CHKN fill the stack, increment the counter
            if (refmv_count[ref_frame]) { //CHKN RefMvCount=1
                assert(refmv_count[ref_frame] == 1);
                if (comp_list[0][0].as_int == ref_mv_stack[ref_frame][0].this_mv.as_int &&
                    comp_list[0][1].as_int == ref_mv_stack[ref_frame][0].comp_mv.as_int) {
                    ref_mv_stack[ref_frame][refmv_count[ref_frame]].this_mv = comp_list[1][0];
                    ref_mv_stack[ref_frame][refmv_count[ref_frame]].comp_mv = comp_list[1][1];
                } else {
                    ref_mv_stack[ref_frame][refmv_count[ref_frame]].this_mv = comp_list[0][0];
                    ref_mv_stack[ref_frame][refmv_count[ref_frame]].comp_mv = comp_list[0][1];
                }
                ref_mv_stack[ref_frame][refmv_count[ref_frame]].weight = 2;
                ++refmv_count[ref_frame];
            } else { //CHKN RefMvCount=0
                for (int32_t idx = 0; idx < MAX_MV_REF_CANDIDATES; ++idx) {
                    ref_mv_stack[ref_frame][refmv_count[ref_frame]].this_mv = comp_list[idx][0];
                    ref_mv_stack[ref_frame][refmv_count[ref_frame]].comp_mv = comp_list[idx][1];
                    ref_mv_stack[ref_frame][refmv_count[ref_frame]].weight  = 2;
                    ++refmv_count[ref_frame];
                }
            }
        }

        assert(refmv_count[ref_frame] >= 2);

        for (int32_t idx = 0; idx < refmv_count[ref_frame]; ++idx) {
            clamp_mv_ref(&ref_mv_stack[ref_frame][idx].this_mv.as_mv,
                         xd->n8_w << MI_SIZE_LOG2,
                         xd->n8_h << MI_SIZE_LOG2,
                         xd);
            clamp_mv_ref(&ref_mv_stack[ref_frame][idx].comp_mv.as_mv,
                         xd->n8_w << MI_SIZE_LOG2,
                         xd->n8_h << MI_SIZE_LOG2,
                         xd);
        }
    } else {
        // Handle single reference frame extension
        int32_t mi_width  = AOMMIN(mi_size_wide[BLOCK_64X64], xd->n8_w);
        mi_width          = AOMMIN(mi_width, cm->mi_cols - mi_col);
        int32_t mi_height = AOMMIN(mi_size_high[BLOCK_64X64], xd->n8_h);
        mi_height         = AOMMIN(mi_height, cm->mi_rows - mi_row);
        int32_t mi_size   = AOMMIN(mi_width, mi_height);

        //CHKn if count is still < 2, re-scan ROW=-1 with less constraints.
        //     Order is already fixed. the added candidates are stored as we go at the bottom of the Stack.
        //CHKN TODO: confirm this could be avoided if we have already 2(DRL:OFF), or 4(DRL:ON) candidates
        for (int32_t idx = 0; abs(max_row_offset) >= 1 && idx < mi_size &&
                              refmv_count[ref_frame] < MAX_MV_REF_CANDIDATES;) {
            const ModeInfo *const   candidate_mi    = xd->mi[-xd->mi_stride + idx];
            const MbModeInfo *const candidate       = &candidate_mi->mbmi;
            const int32_t           candidate_bsize = candidate->block_mi.sb_type;

            for (int32_t rf_idx = 0; rf_idx < 2; ++rf_idx) {
                if (candidate->block_mi.ref_frame[rf_idx] > INTRA_FRAME) {
                    IntMv this_mv = candidate->block_mi.mv[rf_idx];
                    if (cm->ref_frame_sign_bias[candidate->block_mi.ref_frame[rf_idx]] !=
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
            const ModeInfo *const   candidate_mi    = xd->mi[idx * xd->mi_stride - 1];
            const MbModeInfo *const candidate       = &candidate_mi->mbmi;
            const int32_t           candidate_bsize = candidate->block_mi.sb_type;

            for (int32_t rf_idx = 0; rf_idx < 2; ++rf_idx) {
                if (candidate->block_mi.ref_frame[rf_idx] > INTRA_FRAME) {
                    IntMv this_mv = candidate->block_mi.mv[rf_idx];
                    if (cm->ref_frame_sign_bias[candidate->block_mi.ref_frame[rf_idx]] !=
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
                         xd->n8_w << MI_SIZE_LOG2,
                         xd->n8_h << MI_SIZE_LOG2,
                         xd);
        }

        for (int32_t idx = 0; idx < AOMMIN(MAX_MV_REF_CANDIDATES, refmv_count[ref_frame]); ++idx) {
            mv_ref_list[rf[0]][idx].as_int = ref_mv_stack[ref_frame][idx].this_mv.as_int;
        }
    }
    (void)nearest_match;
}

static INLINE int convert_to_trans_prec(int allow_hp, int coor) {
    if (allow_hp)
        return ROUND_POWER_OF_TWO_SIGNED(coor, WARPEDMODEL_PREC_BITS - 3);
    else
        return ROUND_POWER_OF_TWO_SIGNED(coor, WARPEDMODEL_PREC_BITS - 2) * 2;
}

static INLINE int block_center_x(int mi_col, BlockSize bs) {
    const int bw = block_size_wide[bs];
    return mi_col * MI_SIZE + bw / 2 - 1;
}

static INLINE int block_center_y(int mi_row, BlockSize bs) {
    const int bh = block_size_high[bs];
    return mi_row * MI_SIZE + bh / 2 - 1;
}

IntMv gm_get_motion_vector_enc(const EbWarpedMotionParams *gm, int32_t allow_hp, BlockSize bsize,
                               int32_t mi_col, int32_t mi_row, int32_t is_integer) {
    IntMv res;

    if (gm->wmtype == IDENTITY) {
        res.as_int = 0;
        return res;
    }

    const int32_t *mat = gm->wmmat;
    int            x, y, tx, ty;

    if (gm->wmtype == TRANSLATION) {
        // All global motion vectors are stored with WARPEDMODEL_PREC_BITS (16)
        // bits of fractional precision. The offset for a translation is stored in
        // entries 0 and 1. For translations, all but the top three (two if
        // cm->allow_high_precision_mv is false) fractional bits are always zero.
        //
        // After the right shifts, there are 3 fractional bits of precision. If
        // allow_hp is false, the bottom bit is always zero (so we don't need a
        // call to convert_to_trans_prec here)
        res.as_mv.row = gm->wmmat[0] >> GM_TRANS_ONLY_PREC_DIFF;
        res.as_mv.col = gm->wmmat[1] >> GM_TRANS_ONLY_PREC_DIFF;
        assert(IMPLIES(1 & (res.as_mv.row | res.as_mv.col), allow_hp));
        if (is_integer) { integer_mv_precision(&res.as_mv); }
        return res;
    }

    x = block_center_x(mi_col, bsize);
    y = block_center_y(mi_row, bsize);

    if (gm->wmtype == ROTZOOM) {
        assert(gm->wmmat[5] == gm->wmmat[2]);
        assert(gm->wmmat[4] == -gm->wmmat[3]);
    }

    const int xc = (mat[2] - (1 << WARPEDMODEL_PREC_BITS)) * x + mat[3] * y + mat[0];
    const int yc = mat[4] * x + (mat[5] - (1 << WARPEDMODEL_PREC_BITS)) * y + mat[1];
    tx           = convert_to_trans_prec(allow_hp, xc);
    ty           = convert_to_trans_prec(allow_hp, yc);

    res.as_mv.row = ty;
    res.as_mv.col = tx;

    if (is_integer) { integer_mv_precision(&res.as_mv); }
    return res;
}
void mvp_bypass_init(PictureControlSet *pcs_ptr, ModeDecisionContext *context_ptr) {
    TileInfo *tile = &context_ptr->sb_ptr->tile_info;

    int32_t       mi_row = context_ptr->blk_origin_y >> MI_SIZE_LOG2;
    int32_t       mi_col = context_ptr->blk_origin_x >> MI_SIZE_LOG2;
    Av1Common *   cm     = pcs_ptr->parent_pcs_ptr->av1_cm;
    MacroBlockD * xd     = context_ptr->blk_ptr->av1xd;
    BlockSize     bsize  = context_ptr->blk_geom->bsize;
    const int32_t bw     = mi_size_wide[bsize];
    const int32_t bh     = mi_size_high[bsize];

    xd->n8_w = context_ptr->blk_geom->bwidth >> MI_SIZE_LOG2;
    xd->n8_h = context_ptr->blk_geom->bheight >> MI_SIZE_LOG2;
    xd->n4_w = context_ptr->blk_geom->bwidth >> MI_SIZE_LOG2;
    xd->n4_h = context_ptr->blk_geom->bheight >> MI_SIZE_LOG2;

    xd->mb_to_top_edge    = -((mi_row * MI_SIZE) * 8);
    xd->mb_to_bottom_edge = ((cm->mi_rows - bh - mi_row) * MI_SIZE) * 8;
    xd->mb_to_left_edge   = -((mi_col * MI_SIZE) * 8);
    xd->mb_to_right_edge  = ((cm->mi_cols - bw - mi_col) * MI_SIZE) * 8;

    xd->up_available   = (mi_row > tile->mi_row_start);
    xd->left_available = (mi_col > tile->mi_col_start);

    xd->n8_h        = bh;
    xd->n8_w        = bw;
    xd->is_sec_rect = 0;
    if (xd->n8_w < xd->n8_h) {
        // Only mark is_sec_rect as 1 for the last block.
        // For PARTITION_VERT_4, it would be (0, 0, 0, 1);
        // For other partitions, it would be (0, 1).
        if (!((mi_col + xd->n8_w) & (xd->n8_h - 1))) xd->is_sec_rect = 1;
    }

    if (xd->n8_w > xd->n8_h)
        if (mi_row & (xd->n8_w - 1)) xd->is_sec_rect = 1;

    xd->tile.mi_col_start = tile->mi_col_start;
    xd->tile.mi_col_end   = tile->mi_col_end;
    xd->tile.mi_row_start = tile->mi_row_start;
    xd->tile.mi_row_end   = tile->mi_row_end;

    xd->mi_stride        = pcs_ptr->mi_stride;
    const int32_t offset = mi_row * xd->mi_stride + mi_col;
    xd->mi               = pcs_ptr->mi_grid_base + offset;

    xd->mi[0]->mbmi.block_mi.partition = from_shape_to_part[context_ptr->blk_geom->shape];

    // Set to 0 the fields which would have been set by setup_ref_mv_list()
    memset(xd->ref_mv_count, 0, sizeof(uint8_t) * MODE_CTX_REF_FRAMES);
    memset(context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack,
           0,
           sizeof(CandidateMv) * MODE_CTX_REF_FRAMES * MAX_REF_MV_STACK_SIZE);
    memset(context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds].ref_mvs,
           0,
           sizeof(context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds].ref_mvs[0][0]) *
               MODE_CTX_REF_FRAMES * MAX_MV_REF_CANDIDATES);
    memset(context_ptr->blk_ptr->inter_mode_ctx, 0, sizeof(int16_t) * MODE_CTX_REF_FRAMES);
    memset(xd->ref_mv_count, 0, sizeof(int8_t) * MODE_CTX_REF_FRAMES);
}

void generate_av1_mvp_table(TileInfo *tile, ModeDecisionContext *context_ptr, BlkStruct *blk_ptr,
                            const BlockGeom *blk_geom, uint16_t blk_origin_x, uint16_t blk_origin_y,
                            MvReferenceFrame *ref_frames, uint32_t tot_refs,
                            PictureControlSet *pcs_ptr) {
    int32_t      mi_row  = blk_origin_y >> MI_SIZE_LOG2;
    int32_t      mi_col  = blk_origin_x >> MI_SIZE_LOG2;
    Av1Common *  cm      = pcs_ptr->parent_pcs_ptr->av1_cm;
    FrameHeader *frm_hdr = &pcs_ptr->parent_pcs_ptr->frm_hdr;
    MacroBlockD *xd      = blk_ptr->av1xd;
    xd->n8_w             = blk_geom->bwidth >> MI_SIZE_LOG2;
    xd->n8_h             = blk_geom->bheight >> MI_SIZE_LOG2;
    xd->n4_w             = blk_geom->bwidth >> MI_SIZE_LOG2;
    xd->n4_h             = blk_geom->bheight >> MI_SIZE_LOG2;
    BlockSize     bsize  = blk_geom->bsize;
    const int32_t bw     = mi_size_wide[bsize];
    const int32_t bh     = mi_size_high[bsize];

    xd->mb_to_top_edge    = -((mi_row * MI_SIZE) * 8);
    xd->mb_to_bottom_edge = ((cm->mi_rows - bh - mi_row) * MI_SIZE) * 8;
    xd->mb_to_left_edge   = -((mi_col * MI_SIZE) * 8);
    xd->mb_to_right_edge  = ((cm->mi_cols - bw - mi_col) * MI_SIZE) * 8;

    memset(xd->ref_mv_count, 0, sizeof(xd->ref_mv_count));
    memset(context_ptr->md_local_blk_unit[blk_geom->blkidx_mds].ed_ref_mv_stack,
           0,
           sizeof(context_ptr->md_local_blk_unit[blk_geom->blkidx_mds].ed_ref_mv_stack));

    xd->up_available   = (mi_row > tile->mi_row_start);
    xd->left_available = (mi_col > tile->mi_col_start);

    xd->n8_h        = bh;
    xd->n8_w        = bw;
    xd->is_sec_rect = 0;
    if (xd->n8_w < xd->n8_h) {
        // Only mark is_sec_rect as 1 for the last block.
        // For PARTITION_VERT_4, it would be (0, 0, 0, 1);
        // For other partitions, it would be (0, 1).
        if (!((mi_col + xd->n8_w) & (xd->n8_h - 1))) xd->is_sec_rect = 1;
    }

    if (xd->n8_w > xd->n8_h)
        if (mi_row & (xd->n8_w - 1)) xd->is_sec_rect = 1;

    xd->tile.mi_col_start = tile->mi_col_start;
    xd->tile.mi_col_end   = tile->mi_col_end;
    xd->tile.mi_row_start = tile->mi_row_start;
    xd->tile.mi_row_end   = tile->mi_row_end;

    xd->mi_stride        = pcs_ptr->mi_stride;
    const int32_t offset = mi_row * xd->mi_stride + mi_col;
    xd->mi               = pcs_ptr->mi_grid_base + offset;

    xd->mi[0]->mbmi.block_mi.partition = from_shape_to_part[blk_geom->shape];

    uint32_t ref_it;
    for (ref_it = 0; ref_it < tot_refs; ++ref_it) {
        MvReferenceFrame ref_frame = ref_frames[ref_it];

        MvReferenceFrame rf[2];
        av1_set_ref_frame(rf, ref_frame);

        IntMv gm_mv[2];

        if (ref_frame == INTRA_FRAME) {
            gm_mv[0].as_int = gm_mv[1].as_int = 0;
        } else {
            if (ref_frame < REF_FRAMES) {
                gm_mv[0] =
                    gm_get_motion_vector_enc(&pcs_ptr->parent_pcs_ptr->global_motion[ref_frame],
                                             frm_hdr->allow_high_precision_mv,
                                             bsize,
                                             mi_col,
                                             mi_row,
                                             frm_hdr->force_integer_mv);
                gm_mv[1].as_int = 0;
            } else {
                MvReferenceFrame rf[2];
                av1_set_ref_frame(rf, ref_frame);
                gm_mv[0] = gm_get_motion_vector_enc(&pcs_ptr->parent_pcs_ptr->global_motion[rf[0]],
                                                    frm_hdr->allow_high_precision_mv,
                                                    bsize,
                                                    mi_col,
                                                    mi_row,
                                                    frm_hdr->force_integer_mv);
                gm_mv[1] = gm_get_motion_vector_enc(&pcs_ptr->parent_pcs_ptr->global_motion[rf[1]],
                                                    frm_hdr->allow_high_precision_mv,
                                                    bsize,
                                                    mi_col,
                                                    mi_row,
                                                    frm_hdr->force_integer_mv);
            }
        }

        setup_ref_mv_list(pcs_ptr,
                          cm,
                          xd,
                          ref_frame,
                          xd->ref_mv_count,
                          context_ptr->md_local_blk_unit[blk_geom->blkidx_mds].ed_ref_mv_stack,
                          context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds].ref_mvs,
                          gm_mv,
                          pcs_ptr->parent_pcs_ptr->global_motion,
                          mi_row,
                          mi_col,
                          blk_ptr->inter_mode_ctx);
    }
}
void get_av1_mv_pred_drl(ModeDecisionContext *context_ptr, BlkStruct *blk_ptr,
                         MvReferenceFrame ref_frame, uint8_t is_compound, PredictionMode mode,
                         uint8_t drl_index, //valid value of drl_index
                         IntMv nearestmv[2], IntMv nearmv[2], IntMv ref_mv[2]) {
    MacroBlockD *xd = blk_ptr->av1xd;

    if (!is_compound && mode != GLOBALMV) {
        //av1_find_best_ref_mvs(allow_hp, ref_mvs[mbmi->ref_frame[0]], &nearestmv[0], &nearmv[0], cm->cur_frame_force_integer_mv);
        nearestmv[0] =
            context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds].ref_mvs[ref_frame][0];
        nearmv[0] =
            context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds].ref_mvs[ref_frame][1];
    }

    if (is_compound && mode != GLOBAL_GLOBALMV) {
        int32_t ref_mv_idx = drl_index + 1;
        nearestmv[0] =
            context_ptr->md_local_blk_unit[blk_ptr->mds_idx].ed_ref_mv_stack[ref_frame][0].this_mv;
        nearestmv[1] =
            context_ptr->md_local_blk_unit[blk_ptr->mds_idx].ed_ref_mv_stack[ref_frame][0].comp_mv;
        nearmv[0] = context_ptr->md_local_blk_unit[blk_ptr->mds_idx]
                        .ed_ref_mv_stack[ref_frame][ref_mv_idx]
                        .this_mv;
        nearmv[1] = context_ptr->md_local_blk_unit[blk_ptr->mds_idx]
                        .ed_ref_mv_stack[ref_frame][ref_mv_idx]
                        .comp_mv;
    } else if (drl_index > 0 && mode == NEARMV) {
        assert((1 + drl_index) < MAX_REF_MV_STACK_SIZE);
        IntMv cur_mv = context_ptr->md_local_blk_unit[blk_ptr->mds_idx]
                           .ed_ref_mv_stack[ref_frame][1 + drl_index]
                           .this_mv;
        nearmv[0] = cur_mv;
    }

    ref_mv[0] = nearestmv[0];
    ref_mv[1] = nearestmv[1];

    if (is_compound) {
        int32_t ref_mv_idx = drl_index;
        // Special case: NEAR_NEWMV and NEW_NEARMV modes use
        // 1 + mbmi->ref_mv_idx (like NEARMV) instead of
        // mbmi->ref_mv_idx (like NEWMV)
        if (mode == NEAR_NEWMV || mode == NEW_NEARMV) ref_mv_idx = 1 + drl_index;

        if (compound_ref0_mode(mode) == NEWMV)
            ref_mv[0] = context_ptr->md_local_blk_unit[blk_ptr->mds_idx]
                            .ed_ref_mv_stack[ref_frame][ref_mv_idx]
                            .this_mv;

        if (compound_ref1_mode(mode) == NEWMV)
            ref_mv[1] = context_ptr->md_local_blk_unit[blk_ptr->mds_idx]
                            .ed_ref_mv_stack[ref_frame][ref_mv_idx]
                            .comp_mv;
    } else {
        if (mode == NEWMV) {
            if (xd->ref_mv_count[ref_frame] > 1)
                ref_mv[0] = context_ptr->md_local_blk_unit[blk_ptr->mds_idx]
                                .ed_ref_mv_stack[ref_frame][drl_index]
                                .this_mv;
        }
    }
}

void enc_pass_av1_mv_pred(TileInfo *tile, ModeDecisionContext *md_context_ptr, BlkStruct *blk_ptr,
                          const BlockGeom *blk_geom, uint16_t blk_origin_x, uint16_t blk_origin_y,
                          PictureControlSet *pcs_ptr, MvReferenceFrame ref_frame,
                          uint8_t is_compound, PredictionMode mode,
                          IntMv ref_mv[2]) { //[OUT]

    (void)mode;
    IntMv nearestmv[2], nearmv[2];
    memset(nearestmv, 0, sizeof(nearestmv));

    generate_av1_mvp_table(tile,
                           md_context_ptr,
                           blk_ptr,
                           blk_geom,
                           blk_origin_x,
                           blk_origin_y,
                           &ref_frame,
                           1,
                           pcs_ptr);

    get_av1_mv_pred_drl(md_context_ptr,
                        blk_ptr,
                        ref_frame,
                        is_compound,
                        blk_ptr->pred_mode,
                        blk_ptr->drl_index,
                        nearestmv,
                        nearmv,
                        ref_mv);
}

void update_av1_mi_map(BlkStruct *blk_ptr, uint32_t blk_origin_x, uint32_t blk_origin_y,
                       const BlockGeom *blk_geom, PictureControlSet *pcs_ptr) {
    uint32_t mi_stride = pcs_ptr->mi_stride;
    int32_t  mi_row    = blk_origin_y >> MI_SIZE_LOG2;
    int32_t  mi_col    = blk_origin_x >> MI_SIZE_LOG2;

    const int32_t    offset = mi_row * mi_stride + mi_col;
    ModeInfo *       mi_ptr = *(pcs_ptr->mi_grid_base + offset);
    MvReferenceFrame rf[2];
    av1_set_ref_frame(rf, blk_ptr->prediction_unit_array->ref_frame_type);
    uint8_t mi_x, mi_y;

    for (mi_y = 0; mi_y < (blk_geom->bheight >> MI_SIZE_LOG2); mi_y++) {
        for (mi_x = 0; mi_x < (blk_geom->bwidth >> MI_SIZE_LOG2); mi_x++) {
            //these needed for mvPred
            {
                mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.mode = blk_ptr->pred_mode;
                mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.uv_mode =
                    blk_ptr->prediction_unit_array->intra_chroma_mode;
                if (blk_ptr->prediction_mode_flag == INTRA_MODE &&
                    blk_ptr->pred_mode == INTRA_MODE_4x4) {
                    mi_ptr[mi_x + mi_y * mi_stride].mbmi.tx_size          = 0;
                    mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.sb_type = BLOCK_4X4;
                    mi_ptr[mi_x + mi_y * mi_stride].mbmi.tx_depth         = blk_ptr->tx_depth;
                } else {
                    mi_ptr[mi_x + mi_y * mi_stride].mbmi.tx_size =
                        blk_geom->txsize[blk_ptr->tx_depth]
                                        [0]; // inherit tx_size from 1st transform block
                    mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.sb_type = blk_geom->bsize;
                    mi_ptr[mi_x + mi_y * mi_stride].mbmi.tx_depth         = blk_ptr->tx_depth;
                }
                mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.use_intrabc =
                    blk_ptr->av1xd->use_intrabc;
                mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.ref_frame[0] = rf[0];
                mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.ref_frame[1] =
                    blk_ptr->is_interintra_used ? INTRA_FRAME : rf[1];
                if (blk_ptr->prediction_unit_array->inter_pred_direction_index == UNI_PRED_LIST_0) {
                    mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.mv[0].as_mv.col =
                        blk_ptr->prediction_unit_array->mv[0].x;
                    mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.mv[0].as_mv.row =
                        blk_ptr->prediction_unit_array->mv[0].y;
                } else if (blk_ptr->prediction_unit_array->inter_pred_direction_index ==
                           UNI_PRED_LIST_1) {
                    mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.mv[0].as_mv.col =
                        blk_ptr->prediction_unit_array->mv[1].x;
                    mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.mv[0].as_mv.row =
                        blk_ptr->prediction_unit_array->mv[1].y;
                } else {
                    mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.mv[0].as_mv.col =
                        blk_ptr->prediction_unit_array->mv[0].x;
                    mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.mv[0].as_mv.row =
                        blk_ptr->prediction_unit_array->mv[0].y;
                    mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.mv[1].as_mv.col =
                        blk_ptr->prediction_unit_array->mv[1].x;
                    mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.mv[1].as_mv.row =
                        blk_ptr->prediction_unit_array->mv[1].y;
                }

                mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.partition =
                    from_shape_to_part[blk_geom->shape]; // blk_ptr->part;
            }
            mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.interp_filters = blk_ptr->interp_filters;
            mi_ptr[mi_x + mi_y * mi_stride].mbmi.comp_group_idx          = blk_ptr->comp_group_idx;
            mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.compound_idx   = blk_ptr->compound_idx;
            mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.interintra_mode_params.interintra_mode =
                blk_ptr->interintra_mode;
            mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.interintra_mode_params.wedge_interintra =
                blk_ptr->use_wedge_interintra;
            mi_ptr[mi_x + mi_y * mi_stride]
                .mbmi.block_mi.interintra_mode_params.interintra_wedge_index =
                blk_ptr->interintra_wedge_index; //in
            memcpy(&mi_ptr[mi_x + mi_y * mi_stride].mbmi.palette_mode_info,
                   &blk_ptr->palette_info.pmi,
                   sizeof(PaletteModeInfo));
            //needed for CDEF
            mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.skip =
                blk_ptr->block_has_coeff ? EB_FALSE : EB_TRUE;
            mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.skip_mode  = (int8_t)blk_ptr->skip_flag;
            mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.segment_id = blk_ptr->segment_id;
            mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.seg_id_predicted =
                blk_ptr->seg_id_predicted;
            mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.ref_mv_idx = blk_ptr->drl_index;
            mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.motion_mode =
                blk_ptr->prediction_unit_array[0].motion_mode;
            mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.cfl_alpha_idx =
                blk_ptr->prediction_unit_array->cfl_alpha_idx;
            mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.cfl_alpha_signs =
                blk_ptr->prediction_unit_array->cfl_alpha_signs;
            mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.angle_delta[PLANE_TYPE_Y] =
                blk_ptr->prediction_unit_array[0].angle_delta[PLANE_TYPE_Y];
            mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.angle_delta[PLANE_TYPE_UV] =
                blk_ptr->prediction_unit_array[0].angle_delta[PLANE_TYPE_UV];
        }
    }
}

void update_mi_map(struct ModeDecisionContext *context_ptr, BlkStruct *blk_ptr,
                   uint32_t blk_origin_x, uint32_t blk_origin_y, const BlockGeom *blk_geom,
                   uint8_t avail_blk_flag, PictureControlSet *pcs_ptr) {
    uint32_t mi_stride = pcs_ptr->mi_stride;
    int32_t  mi_row    = blk_origin_y >> MI_SIZE_LOG2;
    int32_t  mi_col    = blk_origin_x >> MI_SIZE_LOG2;

    const int32_t    offset = mi_row * mi_stride + mi_col;
    ModeInfo *       mi_ptr = *(pcs_ptr->mi_grid_base + offset);
    MvReferenceFrame rf[2] = {0,0};
    if (avail_blk_flag)
        av1_set_ref_frame(rf, blk_ptr->prediction_unit_array->ref_frame_type);

    uint8_t mi_x, mi_y;
    for (mi_y = 0; mi_y < (blk_geom->bheight >> MI_SIZE_LOG2); mi_y++) {
        for (mi_x = 0; mi_x < (blk_geom->bwidth >> MI_SIZE_LOG2); mi_x++) {
            //these needed for mvPred
            {
                mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.mode = blk_ptr->pred_mode;
                mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.uv_mode =
                    blk_ptr->prediction_unit_array->intra_chroma_mode;
                if (avail_blk_flag && blk_ptr->prediction_mode_flag == INTRA_MODE &&
                    blk_ptr->pred_mode == INTRA_MODE_4x4) {
                    mi_ptr[mi_x + mi_y * mi_stride].mbmi.tx_size          = 0;
                    mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.sb_type = BLOCK_4X4;
                    mi_ptr[mi_x + mi_y * mi_stride].mbmi.tx_depth         = blk_ptr->tx_depth;
                } else {
                    mi_ptr[mi_x + mi_y * mi_stride].mbmi.tx_size =
                        blk_geom->txsize[blk_ptr->tx_depth]
                                        [0]; // inherit tx_size from 1st transform block

                    mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.sb_type = blk_geom->bsize;
                    mi_ptr[mi_x + mi_y * mi_stride].mbmi.tx_depth         = blk_ptr->tx_depth;
                }
                mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.use_intrabc =
                    blk_ptr->av1xd->use_intrabc;
                mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.ref_frame[0] = rf[0];
                mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.ref_frame[1] =
                    avail_blk_flag && blk_ptr->is_interintra_used ?
                    INTRA_FRAME : rf[1];

                if (avail_blk_flag)
                {
                    if (blk_ptr->prediction_unit_array->inter_pred_direction_index == UNI_PRED_LIST_0) {
                        mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.mv[0].as_mv.col =
                            blk_ptr->prediction_unit_array->mv[0].x;
                        mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.mv[0].as_mv.row =
                            blk_ptr->prediction_unit_array->mv[0].y;
                    } else if (blk_ptr->prediction_unit_array->inter_pred_direction_index ==
                               UNI_PRED_LIST_1) {
                        mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.mv[0].as_mv.col =
                            blk_ptr->prediction_unit_array->mv[1].x;
                        mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.mv[0].as_mv.row =
                            blk_ptr->prediction_unit_array->mv[1].y;
                    } else {
                        mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.mv[0].as_mv.col =
                            blk_ptr->prediction_unit_array->mv[0].x;
                        mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.mv[0].as_mv.row =
                            blk_ptr->prediction_unit_array->mv[0].y;
                        mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.mv[1].as_mv.col =
                            blk_ptr->prediction_unit_array->mv[1].x;
                        mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.mv[1].as_mv.row =
                            blk_ptr->prediction_unit_array->mv[1].y;
                    }
                }

                mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.partition =
                    from_shape_to_part[blk_geom->shape]; // blk_ptr->part;
            }

            if (!avail_blk_flag)
                mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.skip = EB_TRUE;
            else {
                if (blk_geom->has_uv && context_ptr->chroma_level <= CHROMA_MODE_1)
                    mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.skip =
                        (context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                                 .y_has_coeff[0] == 0 &&
                         context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                                 .v_has_coeff[0] == 0 &&
                         context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                                 .u_has_coeff[0] == 0)
                            ? EB_TRUE
                            : EB_FALSE;
                else
                    mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.skip =
                        (context_ptr->md_local_blk_unit[context_ptr->blk_geom->blkidx_mds]
                             .y_has_coeff[0] == 0)
                            ? EB_TRUE
                            : EB_FALSE;
            }

            mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.interp_filters = blk_ptr->interp_filters;
            mi_ptr[mi_x + mi_y * mi_stride].mbmi.comp_group_idx          = blk_ptr->comp_group_idx;
            mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.compound_idx   = blk_ptr->compound_idx;
            mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.interintra_mode_params.interintra_mode =
                blk_ptr->interintra_mode;
            mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.interintra_mode_params.wedge_interintra =
                blk_ptr->use_wedge_interintra;
            mi_ptr[mi_x + mi_y * mi_stride]
                .mbmi.block_mi.interintra_mode_params.interintra_wedge_index =
                blk_ptr->interintra_wedge_index; //in
            memcpy(&mi_ptr[mi_x + mi_y * mi_stride].mbmi.palette_mode_info,
                   &blk_ptr->palette_info.pmi,
                   sizeof(PaletteModeInfo));
            mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.skip_mode  = (int8_t)blk_ptr->skip_flag;
            mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.segment_id = blk_ptr->segment_id;
            mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.seg_id_predicted =
                blk_ptr->seg_id_predicted;
            mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.ref_mv_idx = blk_ptr->drl_index;
            mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.motion_mode =
                blk_ptr->prediction_unit_array[0].motion_mode;
            mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.cfl_alpha_idx =
                blk_ptr->prediction_unit_array->cfl_alpha_idx;
            mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.cfl_alpha_signs =
                blk_ptr->prediction_unit_array->cfl_alpha_signs;
            mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.angle_delta[PLANE_TYPE_Y] =
                blk_ptr->prediction_unit_array[0].angle_delta[PLANE_TYPE_Y];
            mi_ptr[mi_x + mi_y * mi_stride].mbmi.block_mi.angle_delta[PLANE_TYPE_UV] =
                blk_ptr->prediction_unit_array[0].angle_delta[PLANE_TYPE_UV];
        }
    }
}

static INLINE void record_samples(MbModeInfo *mbmi, int *pts, int *pts_inref, int row_offset,
                                  int sign_r, int col_offset, int sign_c) {
    uint8_t bw = block_size_wide[mbmi->block_mi.sb_type];
    uint8_t bh = block_size_high[mbmi->block_mi.sb_type];
    int     x  = col_offset * MI_SIZE + sign_c * AOMMAX(bw, MI_SIZE) / 2 - 1;
    int     y  = row_offset * MI_SIZE + sign_r * AOMMAX(bh, MI_SIZE) / 2 - 1;

    pts[0]       = (x * 8);
    pts[1]       = (y * 8);
    pts_inref[0] = (x * 8) + mbmi->block_mi.mv[0].as_mv.col;
    pts_inref[1] = (y * 8) + mbmi->block_mi.mv[0].as_mv.row;
}

// Note: Samples returned are at 1/8-pel precision
// Sample are the neighbor block center point's coordinates relative to the
// left-top pixel of current block.
int av1_find_samples(const Av1Common *cm, const BlockSize sb_size, MacroBlockD *xd, int mi_row, int mi_col,
                     MvReferenceFrame rf0, int *pts, int *pts_inref) {
    int up_available   = xd->up_available;
    int left_available = xd->left_available;
    int i, mi_step = 1, np = 0;

    const TileInfo *const tile  = &xd->tile;
    int                   do_tl = 1;
    int                   do_tr = 1;

    // scan the nearest above rows
    if (up_available) {
        int         mi_row_offset = -1;
        MbModeInfo *mbmi          = &xd->mi[mi_row_offset * xd->mi_stride]->mbmi;
        uint8_t     n4_w          = mi_size_wide[mbmi->block_mi.sb_type];

        if (xd->n4_w <= n4_w) {
            // Handle "current block width <= above block width" case.
            int col_offset = -mi_col % n4_w;

            if (col_offset < 0) do_tl = 0;
            if (col_offset + n4_w > xd->n4_w) do_tr = 0;

            if (mbmi->block_mi.ref_frame[0] == rf0 && mbmi->block_mi.ref_frame[1] == NONE_FRAME) {
                record_samples(mbmi, pts, pts_inref, 0, -1, col_offset, 1);
                pts += 2;
                pts_inref += 2;
                np++;
                if (np >= LEAST_SQUARES_SAMPLES_MAX) return LEAST_SQUARES_SAMPLES_MAX;
            }
        } else {
            // Handle "current block width > above block width" case.
            for (i = 0; i < AOMMIN(xd->n4_w, cm->mi_cols - mi_col); i += mi_step) {
                int mi_col_offset = i;
                mbmi              = &xd->mi[mi_col_offset + mi_row_offset * xd->mi_stride]->mbmi;
                n4_w              = mi_size_wide[mbmi->block_mi.sb_type];
                mi_step           = AOMMIN(xd->n4_w, n4_w);

                if (mbmi->block_mi.ref_frame[0] == rf0 &&
                    mbmi->block_mi.ref_frame[1] == NONE_FRAME) {
                    record_samples(mbmi, pts, pts_inref, 0, -1, i, 1);
                    pts += 2;
                    pts_inref += 2;
                    np++;
                    if (np >= LEAST_SQUARES_SAMPLES_MAX) return LEAST_SQUARES_SAMPLES_MAX;
                }
            }
        }
    }

    // scan the nearest left columns
    if (left_available) {
        int         mi_col_offset = -1;
        MbModeInfo *mbmi          = &xd->mi[mi_col_offset]->mbmi;
        uint8_t     n4_h          = mi_size_high[mbmi->block_mi.sb_type];

        if (xd->n4_h <= n4_h) {
            // Handle "current block height <= above block height" case.
            int row_offset = -mi_row % n4_h;
            if (row_offset < 0) do_tl = 0;

            if (mbmi->block_mi.ref_frame[0] == rf0 && mbmi->block_mi.ref_frame[1] == NONE_FRAME) {
                record_samples(mbmi, pts, pts_inref, row_offset, 1, 0, -1);
                pts += 2;
                pts_inref += 2;
                np++;
                if (np >= LEAST_SQUARES_SAMPLES_MAX) return LEAST_SQUARES_SAMPLES_MAX;
            }
        } else {
            // Handle "current block height > above block height" case.
            for (i = 0; i < AOMMIN(xd->n4_h, cm->mi_rows - mi_row); i += mi_step) {
                int mi_row_offset = i;
                mbmi              = &xd->mi[mi_col_offset + mi_row_offset * xd->mi_stride]->mbmi;
                n4_h              = mi_size_high[mbmi->block_mi.sb_type];
                mi_step           = AOMMIN(xd->n4_h, n4_h);

                if (mbmi->block_mi.ref_frame[0] == rf0 &&
                    mbmi->block_mi.ref_frame[1] == NONE_FRAME) {
                    record_samples(mbmi, pts, pts_inref, i, 1, 0, -1);
                    pts += 2;
                    pts_inref += 2;
                    np++;
                    if (np >= LEAST_SQUARES_SAMPLES_MAX) return LEAST_SQUARES_SAMPLES_MAX;
                }
            }
        }
    }

    // Top-left block
    if (do_tl && left_available && up_available) {
        int         mi_row_offset = -1;
        int         mi_col_offset = -1;
        MbModeInfo *mbmi          = &xd->mi[mi_col_offset + mi_row_offset * xd->mi_stride]->mbmi;

        if (mbmi->block_mi.ref_frame[0] == rf0 && mbmi->block_mi.ref_frame[1] == NONE_FRAME) {
            record_samples(mbmi, pts, pts_inref, 0, -1, 0, -1);
            pts += 2;
            pts_inref += 2;
            np++;
            if (np >= LEAST_SQUARES_SAMPLES_MAX) return LEAST_SQUARES_SAMPLES_MAX;
        }
    }

    // Top-right block
    if (do_tr && has_top_right(cm, sb_size, xd, mi_row, mi_col, AOMMAX(xd->n4_w, xd->n4_h))) {
        Position trb_pos = {-1, xd->n4_w};

        if (is_inside(tile, mi_col, mi_row, cm->mi_rows, &trb_pos)) {
            int mi_row_offset = -1;
            int mi_col_offset = xd->n4_w;

            MbModeInfo *mbmi = &xd->mi[mi_col_offset + mi_row_offset * xd->mi_stride]->mbmi;

            if (mbmi->block_mi.ref_frame[0] == rf0 && mbmi->block_mi.ref_frame[1] == NONE_FRAME) {
                record_samples(mbmi, pts, pts_inref, 0, -1, xd->n4_w, 1);
                np++;
                if (np >= LEAST_SQUARES_SAMPLES_MAX) return LEAST_SQUARES_SAMPLES_MAX;
            }
        }
    }

    return np;
}

void wm_count_samples(BlkStruct *blk_ptr, const BlockSize sb_size, const BlockGeom *blk_geom, uint16_t blk_origin_x,
                      uint16_t blk_origin_y, uint8_t ref_frame_type, PictureControlSet *pcs_ptr,
                      uint16_t *num_samples) {
    Av1Common *  cm = pcs_ptr->parent_pcs_ptr->av1_cm;
    MacroBlockD *xd = blk_ptr->av1xd;

    int32_t mi_row = blk_origin_y >> MI_SIZE_LOG2;
    int32_t mi_col = blk_origin_x >> MI_SIZE_LOG2;

    xd->n4_w = blk_geom->bwidth >> MI_SIZE_LOG2;
    xd->n4_h = blk_geom->bheight >> MI_SIZE_LOG2;

    int up_available   = xd->up_available;
    int left_available = xd->left_available;
    int i, mi_step = 1, np = 0;

    const TileInfo *const tile  = &xd->tile;
    int                   do_tl = 1;
    int                   do_tr = 1;

    MvReferenceFrame rf[2];
    av1_set_ref_frame(rf, ref_frame_type);

    if (up_available) {
        int         mi_row_offset = -1;
        MbModeInfo *mbmi          = &xd->mi[mi_row_offset * xd->mi_stride]->mbmi;
        uint8_t     n4_w          = mi_size_wide[mbmi->block_mi.sb_type];

        if (xd->n4_w <= n4_w) {
            int col_offset = -mi_col % n4_w;
            if (col_offset < 0) do_tl = 0;
            if (col_offset + n4_w > xd->n4_w) do_tr = 0;

            if (mbmi->block_mi.ref_frame[0] == rf[0] && mbmi->block_mi.ref_frame[1] == NONE_FRAME) {
                np++;
                if (np >= LEAST_SQUARES_SAMPLES_MAX) {
                    *num_samples = LEAST_SQUARES_SAMPLES_MAX;
                    return;
                }
            }
        } else {
            for (i = 0; i < AOMMIN(xd->n4_w, cm->mi_cols - mi_col); i += mi_step) {
                int mi_col_offset = i;
                mbmi              = &xd->mi[mi_col_offset + mi_row_offset * xd->mi_stride]->mbmi;
                n4_w              = mi_size_wide[mbmi->block_mi.sb_type];
                mi_step           = AOMMIN(xd->n4_w, n4_w);

                if (mbmi->block_mi.ref_frame[0] == rf[0] &&
                    mbmi->block_mi.ref_frame[1] == NONE_FRAME) {
                    np++;
                    if (np >= LEAST_SQUARES_SAMPLES_MAX) {
                        *num_samples = LEAST_SQUARES_SAMPLES_MAX;
                        return;
                    }
                }
            }
        }
    }

    if (left_available) {
        int         mi_col_offset = -1;
        MbModeInfo *mbmi          = &xd->mi[mi_col_offset]->mbmi;
        uint8_t     n4_h          = mi_size_high[mbmi->block_mi.sb_type];
        if (xd->n4_h <= n4_h) {
            int row_offset = -mi_row % n4_h;
            if (row_offset < 0) do_tl = 0;
            if (mbmi->block_mi.ref_frame[0] == rf[0] && mbmi->block_mi.ref_frame[1] == NONE_FRAME) {
                np++;
                if (np >= LEAST_SQUARES_SAMPLES_MAX) {
                    *num_samples = LEAST_SQUARES_SAMPLES_MAX;
                    return;
                }
            }
        } else {
            for (i = 0; i < AOMMIN(xd->n4_h, cm->mi_rows - mi_row); i += mi_step) {
                int mi_row_offset = i;
                mbmi              = &xd->mi[mi_col_offset + mi_row_offset * xd->mi_stride]->mbmi;
                n4_h              = mi_size_high[mbmi->block_mi.sb_type];
                mi_step           = AOMMIN(xd->n4_h, n4_h);

                if (mbmi->block_mi.ref_frame[0] == rf[0] &&
                    mbmi->block_mi.ref_frame[1] == NONE_FRAME) {
                    np++;
                    if (np >= LEAST_SQUARES_SAMPLES_MAX) {
                        *num_samples = LEAST_SQUARES_SAMPLES_MAX;
                        return;
                    }
                }
            }
        }
    }

    if (do_tl && left_available && up_available) {
        int         mi_row_offset = -1;
        int         mi_col_offset = -1;
        MbModeInfo *mbmi          = &xd->mi[mi_col_offset + mi_row_offset * xd->mi_stride]->mbmi;
        if (mbmi->block_mi.ref_frame[0] == rf[0] && mbmi->block_mi.ref_frame[1] == NONE_FRAME) {
            np++;
            if (np >= LEAST_SQUARES_SAMPLES_MAX) {
                *num_samples = LEAST_SQUARES_SAMPLES_MAX;
                return;
            }
        }
    }

    if (do_tr && has_top_right(cm, sb_size, xd, mi_row, mi_col, AOMMAX(xd->n4_w, xd->n4_h))) {
        Position trb_pos = {-1, xd->n4_w};
        if (is_inside(tile, mi_col, mi_row, cm->mi_rows, &trb_pos)) {
            int         mi_row_offset = -1;
            int         mi_col_offset = xd->n4_w;
            MbModeInfo *mbmi = &xd->mi[mi_col_offset + mi_row_offset * xd->mi_stride]->mbmi;

            if (mbmi->block_mi.ref_frame[0] == rf[0] && mbmi->block_mi.ref_frame[1] == NONE_FRAME) {
                np++;
                if (np >= LEAST_SQUARES_SAMPLES_MAX) {
                    *num_samples = LEAST_SQUARES_SAMPLES_MAX;
                    return;
                }
            }
        }
    }
    *num_samples = np;
}

uint16_t wm_find_samples(BlkStruct *blk_ptr, const BlockGeom *blk_geom, uint16_t blk_origin_x,
                         uint16_t blk_origin_y, MvReferenceFrame rf0, PictureControlSet *pcs_ptr,
                         int32_t *pts, int32_t *pts_inref) {
    Av1Common *  cm = pcs_ptr->parent_pcs_ptr->av1_cm;
    MacroBlockD *xd = blk_ptr->av1xd;

    int32_t mi_row = blk_origin_y >> MI_SIZE_LOG2;
    int32_t mi_col = blk_origin_x >> MI_SIZE_LOG2;

    xd->n4_w = blk_geom->bwidth >> MI_SIZE_LOG2;
    xd->n4_h = blk_geom->bheight >> MI_SIZE_LOG2;

    return (uint16_t)av1_find_samples(cm,
            ((SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr)->seq_header.sb_size,
            xd, mi_row, mi_col, rf0, pts, pts_inref);
}

EbBool warped_motion_parameters(PictureControlSet *pcs_ptr, BlkStruct *blk_ptr, MvUnit *mv_unit,
                                const BlockGeom *blk_geom, uint16_t blk_origin_x,
                                uint16_t blk_origin_y, uint8_t ref_frame_type,
                                EbWarpedMotionParams *wm_params, uint16_t *num_samples) {
    MacroBlockD *xd       = blk_ptr->av1xd;
    BlockSize    bsize    = blk_geom->bsize;
    EbBool       apply_wm = EB_FALSE;

    int     pts[SAMPLES_ARRAY_SIZE], pts_inref[SAMPLES_ARRAY_SIZE];
    int32_t mi_row = blk_origin_y >> MI_SIZE_LOG2;
    int32_t mi_col = blk_origin_x >> MI_SIZE_LOG2;
    xd->n4_w       = blk_geom->bwidth >> MI_SIZE_LOG2;
    xd->n4_h       = blk_geom->bheight >> MI_SIZE_LOG2;

    *num_samples = 0;
    if (blk_geom->bwidth < 8 || blk_geom->bheight < 8) return apply_wm;

    MvReferenceFrame rf[2];
    av1_set_ref_frame(rf, ref_frame_type);

    uint16_t nsamples = wm_find_samples(
        blk_ptr, blk_geom, blk_origin_x, blk_origin_y, rf[0], pcs_ptr, pts, pts_inref);

    if (nsamples == 0) return apply_wm;

    MV mv;
    mv.col = mv_unit->mv[mv_unit->pred_direction % BI_PRED].x;
    mv.row = mv_unit->mv[mv_unit->pred_direction % BI_PRED].y;
    if (nsamples > 1) nsamples = select_samples(&mv, pts, pts_inref, nsamples, bsize);
    *num_samples = nsamples;

    apply_wm = !eb_find_projection(
        (int)nsamples, pts, pts_inref, bsize, mv.row, mv.col, wm_params, (int)mi_row, (int)mi_col);

    return apply_wm;
}

//foreach_overlappable_nb_above
int count_overlappable_nb_above(const Av1Common *cm, MacroBlockD *xd, int32_t mi_col, int nb_max) {
    int nb_count = 0;
    if (!xd->up_available) return nb_count;

    // prev_row_mi points into the mi array, starting at the beginning of the
    // previous row.
    ModeInfo **prev_row_mi = xd->mi - mi_col - 1 * xd->mi_stride;
    const int  end_col     = MIN(mi_col + xd->n4_w, cm->mi_cols);
    uint8_t    mi_step;

    for (int above_mi_col = mi_col; above_mi_col < end_col && nb_count < nb_max;
         above_mi_col += mi_step) {
        ModeInfo **above_mi = prev_row_mi + above_mi_col;
        mi_step = MIN(mi_size_wide[above_mi[0]->mbmi.block_mi.sb_type], mi_size_wide[BLOCK_64X64]);

        // If we're considering a block with width 4, it should be treated as
        // half of a pair of blocks with chroma information in the second. Move
        // above_mi_col back to the start of the pair if needed, set above_mbmi
        // to point at the block with chroma information, and set mi_step to 2 to
        // step over the entire pair at the end of the iteration.
        if (mi_step == 1) {
            above_mi_col &= ~1;
            above_mi = prev_row_mi + above_mi_col + 1;
            mi_step  = 2;
        }
        if (is_neighbor_overlappable(&(*above_mi)->mbmi)) ++nb_count;
    }

    return nb_count;
}

int count_overlappable_nb_left(const Av1Common *cm, MacroBlockD *xd, int32_t mi_row, int nb_max) {
    int nb_count = 0;
    if (!xd->left_available) return nb_count;

    // prev_col_mi points into the mi array, starting at the top of the
    // previous column
    ModeInfo **prev_col_mi = xd->mi - 1 - mi_row * xd->mi_stride;
    const int  end_row     = MIN(mi_row + xd->n4_h, cm->mi_rows);
    uint8_t    mi_step;

    for (int left_mi_row = mi_row; left_mi_row < end_row && nb_count < nb_max;
         left_mi_row += mi_step) {
        ModeInfo **left_mi = prev_col_mi + left_mi_row * xd->mi_stride;
        mi_step = MIN(mi_size_high[left_mi[0]->mbmi.block_mi.sb_type], mi_size_high[BLOCK_64X64]);
        if (mi_step == 1) {
            left_mi_row &= ~1;
            left_mi = prev_col_mi + (left_mi_row + 1) * xd->mi_stride;
            mi_step = 2;
        }

        if (is_neighbor_overlappable(&(*left_mi)->mbmi)) ++nb_count;
    }

    return nb_count;
}

void eb_av1_count_overlappable_neighbors(const PictureControlSet *pcs_ptr, BlkStruct *blk_ptr,
                                         const BlockSize bsize, int32_t mi_row, int32_t mi_col) {
    Av1Common *  cm                                             = pcs_ptr->parent_pcs_ptr->av1_cm;
    MacroBlockD *xd                                             = blk_ptr->av1xd;
    blk_ptr->prediction_unit_array[0].overlappable_neighbors[0] = 0;
    blk_ptr->prediction_unit_array[0].overlappable_neighbors[1] = 0;

    if (!is_motion_variation_allowed_bsize(bsize)) return;

    blk_ptr->prediction_unit_array[0].overlappable_neighbors[0] =
        count_overlappable_nb_above(cm, xd, mi_col, MAX_SIGNED_VALUE);

    blk_ptr->prediction_unit_array[0].overlappable_neighbors[1] =
        count_overlappable_nb_left(cm, xd, mi_row, MAX_SIGNED_VALUE);
}

int av1_is_dv_valid(const MV dv, const MacroBlockD *xd, int mi_row, int mi_col, BlockSize bsize,
                    int mib_size_log2) {
    const int bw             = block_size_wide[bsize];
    const int bh             = block_size_high[bsize];
    const int scale_px_to_mv = 8;
    // Disallow subpixel for now
    // SUBPEL_MASK is not the correct scale
    if (((dv.row & (scale_px_to_mv - 1)) || (dv.col & (scale_px_to_mv - 1)))) return 0;

    const TileInfo *const tile = &xd->tile;
    // Is the source top-left inside the current tile?
    const int src_top_edge  = mi_row * MI_SIZE * scale_px_to_mv + dv.row;
    const int tile_top_edge = tile->mi_row_start * MI_SIZE * scale_px_to_mv;
    if (src_top_edge < tile_top_edge) return 0;
    const int src_left_edge  = mi_col * MI_SIZE * scale_px_to_mv + dv.col;
    const int tile_left_edge = tile->mi_col_start * MI_SIZE * scale_px_to_mv;
    if (src_left_edge < tile_left_edge) return 0;
    // Is the bottom right inside the current tile?
    const int src_bottom_edge  = (mi_row * MI_SIZE + bh) * scale_px_to_mv + dv.row;
    const int tile_bottom_edge = tile->mi_row_end * MI_SIZE * scale_px_to_mv;
    if (src_bottom_edge > tile_bottom_edge) return 0;
    const int src_right_edge  = (mi_col * MI_SIZE + bw) * scale_px_to_mv + dv.col;
    const int tile_right_edge = tile->mi_col_end * MI_SIZE * scale_px_to_mv;
    if (src_right_edge > tile_right_edge) return 0;

    // Special case for sub 8x8 chroma cases, to prevent referring to chroma
    // pixels outside current tile.
    for (int plane = 1; plane < 3 /* av1_num_planes(cm)*/; ++plane) {
        //const struct MacroBlockDPlane *const pd = &xd->plane[plane];

        if (is_chroma_reference(mi_row, mi_col, bsize, 1, 1/* pd->subsampling_x,
            pd->subsampling_y*/)) {
            if (bw < 8 /*&& pd->subsampling_x*/)
                if (src_left_edge < tile_left_edge + 4 * scale_px_to_mv) return 0;
            if (bh < 8 /* && pd->subsampling_y*/)
                if (src_top_edge < tile_top_edge + 4 * scale_px_to_mv) return 0;
        }
    }

    // Is the bottom right within an already coded SB? Also consider additional
    // constraints to facilitate HW decoder.
    const int max_mib_size       = 1 << mib_size_log2;
    const int active_sb_row      = mi_row >> mib_size_log2;
    const int active_sb64_col    = (mi_col * MI_SIZE) >> 6;
    const int sb_size            = max_mib_size * MI_SIZE;
    const int src_sb_row         = ((src_bottom_edge >> 3) - 1) / sb_size;
    const int src_sb64_col       = ((src_right_edge >> 3) - 1) >> 6;
    const int total_sb64_per_row = ((tile->mi_col_end - tile->mi_col_start - 1) >> 4) + 1;
    const int active_sb64        = active_sb_row * total_sb64_per_row + active_sb64_col;
    const int src_sb64           = src_sb_row * total_sb64_per_row + src_sb64_col;
    if (src_sb64 >= active_sb64 - INTRABC_DELAY_SB64) return 0;

    // Wavefront constraint: use only top left area of frame for reference.
    const int gradient  = 1 + INTRABC_DELAY_SB64 + (sb_size > 64);
    const int wf_offset = gradient * (active_sb_row - src_sb_row);
    if (src_sb_row > active_sb_row ||
        src_sb64_col >= active_sb64_col - INTRABC_DELAY_SB64 + wf_offset)
        return 0;

    //add a SW-Wavefront constraint
    if (sb_size == 64) {
        if (src_sb64_col > active_sb64_col + (active_sb_row - src_sb_row)) return 0;
    } else {
        const int src_sb128_col    = ((src_right_edge >> 3) - 1) >> 7;
        const int active_sb128_col = (mi_col * MI_SIZE) >> 7;

        if (src_sb128_col > active_sb128_col + (active_sb_row - src_sb_row)) return 0;
    }

    return 1;
}

int is_inside_tile_boundary(TileInfo *tile, int16_t mvx, int16_t mvy, int mi_col, int mi_row,
                            BlockSize bsize) {
    const int bw             = block_size_wide[bsize];
    const int bh             = block_size_high[bsize];
    const int scale_px_to_mv = 8;
    const int src_top_edge   = mi_row * MI_SIZE * scale_px_to_mv + mvy;
    const int tile_top_edge  = tile->mi_row_start * MI_SIZE * scale_px_to_mv;
    if (src_top_edge < tile_top_edge) return 0;
    const int src_left_edge  = mi_col * MI_SIZE * scale_px_to_mv + mvx;
    const int tile_left_edge = tile->mi_col_start * MI_SIZE * scale_px_to_mv;
    if (src_left_edge < tile_left_edge) return 0;
    // Is the bottom right inside the current tile?
    const int src_bottom_edge  = (mi_row * MI_SIZE + bh) * scale_px_to_mv + mvy;
    const int tile_bottom_edge = tile->mi_row_end * MI_SIZE * scale_px_to_mv;
    if (src_bottom_edge > tile_bottom_edge) return 0;
    const int src_right_edge  = (mi_col * MI_SIZE + bw) * scale_px_to_mv + mvx;
    const int tile_right_edge = tile->mi_col_end * MI_SIZE * scale_px_to_mv;
    if (src_right_edge > tile_right_edge) return 0;

    return 1;
}

IntMv eb_av1_get_ref_mv_from_stack(int ref_idx, const MvReferenceFrame *ref_frame, int ref_mv_idx,
                                   CandidateMv  ref_mv_stack[][MAX_REF_MV_STACK_SIZE],
                                   MacroBlockD *xd
                                   /*const MB_MODE_INFO_EXT *mbmi_ext*/) {
    const int8_t       ref_frame_type = av1_ref_frame_type(ref_frame);
    const CandidateMv *curr_ref_mv_stack =
        /*mbmi_ext->*/ ref_mv_stack[ref_frame_type];
    IntMv ref_mv;
    ref_mv.as_int = INVALID_MV;

    if (ref_frame[1] > INTRA_FRAME) {
        if (ref_idx == 0)
            ref_mv = curr_ref_mv_stack[ref_mv_idx].this_mv;
        else {
            assert(ref_idx == 1);
            ref_mv = curr_ref_mv_stack[ref_mv_idx].comp_mv;
        }
    } else {
        assert(ref_idx == 0);
        if (ref_mv_idx < /*mbmi_ext->*/ xd->ref_mv_count[ref_frame_type])
            ref_mv = curr_ref_mv_stack[ref_mv_idx].this_mv;
        else {
            //CHKN got this from decoder read_intrabc_info global_mvs[ref_frame].as_int = INVALID_MV;
            ref_mv.as_int = INVALID_MV; // mbmi_ext->global_mvs[ref_frame_type];
        }
    }
    return ref_mv;
}

void eb_av1_find_best_ref_mvs_from_stack(int allow_hp,
                                         //const MB_MODE_INFO_EXT *mbmi_ext,
                                         CandidateMv  ref_mv_stack[][MAX_REF_MV_STACK_SIZE],
                                         MacroBlockD *xd, MvReferenceFrame ref_frame,
                                         IntMv *nearest_mv, IntMv *near_mv, int is_integer) {
    const int        ref_idx       = 0;
    MvReferenceFrame ref_frames[2] = {ref_frame, NONE_FRAME};
    *nearest_mv =
        eb_av1_get_ref_mv_from_stack(ref_idx, ref_frames, 0, ref_mv_stack /*mbmi_ext*/, xd);
    lower_mv_precision(&nearest_mv->as_mv, allow_hp, is_integer);
    *near_mv = eb_av1_get_ref_mv_from_stack(ref_idx, ref_frames, 1, ref_mv_stack /*mbmi_ext*/, xd);
    lower_mv_precision(&near_mv->as_mv, allow_hp, is_integer);
}
