/*
* Copyright(c) 2019 Netflix, Inc.
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

#include "EbDecParseInterBlock.h"
#include "EbCommonUtils.h"
#include "EbWarpedMotion.h"

typedef const int (*ColorCost)[PALETTE_SIZES][PALETTE_COLOR_INDEX_CONTEXTS][PALETTE_COLORS];
typedef AomCdfProb (*MapCdf)[PALETTE_SIZES][PALETTE_COLOR_INDEX_CONTEXTS];

typedef struct {
    int       rows;
    int       n_colors;
    int       plane_width;
    int       plane_height;
} Av1ColorMapParam;

#define MAX_COLOR_CONTEXT_HASH 8
#define NUM_PALETTE_NEIGHBORS 3 // left, top-left and top.
#define COLOR_MAP_STRIDE 128 // worst case

// Negative values are invalid
extern int palette_color_index_context_lookup[MAX_COLOR_CONTEXT_HASH + 1];

static uint16_t compound_mode_ctx_map[3][COMP_NEWMV_CTXS] = {
    {0, 1, 1, 1, 1},
    {1, 2, 3, 4, 4},
    {4, 4, 5, 6, 7},
};

static INLINE void svt_collect_neighbors_ref_counts(PartitionInfo *pi) {
    ZERO_ARRAY(&pi->neighbors_ref_counts[0], sizeof(pi->neighbors_ref_counts[0]) * REF_FRAMES);

    uint8_t *const ref_counts = pi->neighbors_ref_counts;

    const BlockModeInfo *const above_mbmi     = pi->above_mbmi;
    const BlockModeInfo *const left_mbmi      = pi->left_mbmi;
    const int                  above_in_image = pi->up_available;
    const int                  left_in_image  = pi->left_available;

    // Above neighbor
    if (above_in_image && is_inter_block(above_mbmi)) {
        ref_counts[above_mbmi->ref_frame[0]]++;
        if (has_second_ref(above_mbmi)) ref_counts[above_mbmi->ref_frame[1]]++;
    }

    // Left neighbor
    if (left_in_image && is_inter_block(left_mbmi)) {
        ref_counts[left_mbmi->ref_frame[0]]++;
        if (has_second_ref(left_mbmi)) ref_counts[left_mbmi->ref_frame[1]]++;
    }
}

static INLINE int is_inside(TileInfo *tile, int mi_col, int mi_row) {
    return (mi_col >= tile->mi_col_start && mi_col < tile->mi_col_end &&
            mi_row >= tile->mi_row_start && mi_row < tile->mi_row_end);
}

static int get_reference_mode_context(const PartitionInfo *xd) {
    int                        ctx;
    const BlockModeInfo *const above_mbmi = xd->above_mbmi;
    const BlockModeInfo *const left_mbmi  = xd->left_mbmi;
    const int                  has_above  = xd->up_available;
    const int                  has_left   = xd->left_available;

    // Note:
    // The mode info data structure has a one element border above and to the
    // left of the entries corresponding to real macroblocks.
    // The prediction flags in these dummy entries are initialized to 0.
    if (has_above && has_left) { // both edges available
        if (!has_second_ref(above_mbmi) && !has_second_ref(left_mbmi))
            // neither edge uses comp pred (0/1)
            ctx = IS_BACKWARD_REF_FRAME(above_mbmi->ref_frame[0]) ^
                  IS_BACKWARD_REF_FRAME(left_mbmi->ref_frame[0]);
        else if (!has_second_ref(above_mbmi))
            // one of two edges uses comp pred (2/3)
            ctx = 2 +
                  (IS_BACKWARD_REF_FRAME(above_mbmi->ref_frame[0]) || !is_inter_block(above_mbmi));
        else if (!has_second_ref(left_mbmi))
            // one of two edges uses comp pred (2/3)
            ctx =
                2 + (IS_BACKWARD_REF_FRAME(left_mbmi->ref_frame[0]) || !is_inter_block(left_mbmi));
        else // both edges use comp pred (4)
            ctx = 4;
    } else if (has_above || has_left) { // one edge available
        const BlockModeInfo *edge_mbmi = has_above ? above_mbmi : left_mbmi;

        if (!has_second_ref(edge_mbmi))
            // edge does not use comp pred (0/1)
            ctx = IS_BACKWARD_REF_FRAME(edge_mbmi->ref_frame[0]);
        else
            // edge uses comp pred (3)
            ctx = 3;
    } else { // no edges available (1)
        ctx = 1;
    }
    assert(ctx >= 0 && ctx < COMP_INTER_CONTEXTS);
    return ctx;
}

static int32_t get_pred_context_comp_ref_p(PartitionInfo *pi) {
    const uint8_t *const ref_counts = &pi->neighbors_ref_counts[0];

    // Count of LAST + LAST2
    const int32_t last_last2_count = ref_counts[LAST_FRAME] + ref_counts[LAST2_FRAME];
    // Count of LAST3 + GOLDEN
    const int32_t last3_gld_count = ref_counts[LAST3_FRAME] + ref_counts[GOLDEN_FRAME];

    const int32_t pred_context =
        (last_last2_count == last3_gld_count) ? 1 : ((last_last2_count < last3_gld_count) ? 0 : 2);

    assert(pred_context >= 0 && pred_context < REF_CONTEXTS);
    return pred_context;
}

static int32_t get_pred_context_comp_bwdref_p(PartitionInfo *pi) {
    const uint8_t *const ref_counts = &pi->neighbors_ref_counts[0];

    // Counts of BWDREF, ALTREF2, or ALTREF frames (b, A2, or A)
    const int32_t brfarf2_count = ref_counts[BWDREF_FRAME] + ref_counts[ALTREF2_FRAME];
    const int32_t arf_count     = ref_counts[ALTREF_FRAME];

    const int32_t pred_context =
        (brfarf2_count == arf_count) ? 1 : ((brfarf2_count < arf_count) ? 0 : 2);

    assert(pred_context >= 0 && pred_context < REF_CONTEXTS);
    return pred_context;
}

static int32_t get_pred_context_comp_bwdref_p1(PartitionInfo *pi) {
    const uint8_t *const ref_counts = &pi->neighbors_ref_counts[0];

    // Count of BWDREF frames (b)
    const int32_t brf_count = ref_counts[BWDREF_FRAME];
    // Count of ALTREF2 frames (A2)
    const int32_t arf2_count = ref_counts[ALTREF2_FRAME];

    const int32_t pred_context = (brf_count == arf2_count) ? 1 : ((brf_count < arf2_count) ? 0 : 2);

    assert(pred_context >= 0 && pred_context < REF_CONTEXTS);
    return pred_context;
}

static int get_pred_context_uni_comp_ref_p2(PartitionInfo *pi) {
    const uint8_t *const ref_counts = &pi->neighbors_ref_counts[0];

    // Count of LAST3
    const int last3_count = ref_counts[LAST3_FRAME];
    // Count of GOLDEN
    const int gld_count = ref_counts[GOLDEN_FRAME];

    const int pred_context = (last3_count == gld_count) ? 1 : ((last3_count < gld_count) ? 0 : 2);

    assert(pred_context >= 0 && pred_context < UNI_COMP_REF_CONTEXTS);
    return pred_context;
}

static int get_pred_context_uni_comp_ref_p1(PartitionInfo *pi) {
    const uint8_t *const ref_counts = &pi->neighbors_ref_counts[0];

    // Count of LAST2
    const int last2_count = ref_counts[LAST2_FRAME];
    // Count of LAST3 or GOLDEN
    const int last3_or_gld_count = ref_counts[LAST3_FRAME] + ref_counts[GOLDEN_FRAME];

    const int pred_context =
        (last2_count == last3_or_gld_count) ? 1 : ((last2_count < last3_or_gld_count) ? 0 : 2);

    assert(pred_context >= 0 && pred_context < UNI_COMP_REF_CONTEXTS);
    return pred_context;
}

static int get_pred_context_uni_comp_ref_p(PartitionInfo *pi) {
    const uint8_t *const ref_counts = &pi->neighbors_ref_counts[0];

    // Count of forward references (L, L2, L3, or G)
    const int frf_count = ref_counts[LAST_FRAME] + ref_counts[LAST2_FRAME] +
                          ref_counts[LAST3_FRAME] + ref_counts[GOLDEN_FRAME];
    // Count of backward references (b or A)
    const int brf_count =
        ref_counts[BWDREF_FRAME] + ref_counts[ALTREF2_FRAME] + ref_counts[ALTREF_FRAME];

    const int pred_context = (frf_count == brf_count) ? 1 : ((frf_count < brf_count) ? 0 : 2);

    assert(pred_context >= 0 && pred_context < UNI_COMP_REF_CONTEXTS);
    return pred_context;
}

static int32_t get_pred_context_single_ref_p1(PartitionInfo *pi) {
    const uint8_t *const ref_counts = &pi->neighbors_ref_counts[0];

    // Count of forward reference frames
    const int32_t fwd_count = ref_counts[LAST_FRAME] + ref_counts[LAST2_FRAME] +
                              ref_counts[LAST3_FRAME] + ref_counts[GOLDEN_FRAME];
    // Count of backward reference frames
    const int32_t bwd_count =
        ref_counts[BWDREF_FRAME] + ref_counts[ALTREF2_FRAME] + ref_counts[ALTREF_FRAME];

    const int32_t pred_context = (fwd_count == bwd_count) ? 1 : ((fwd_count < bwd_count) ? 0 : 2);

    assert(pred_context >= 0 && pred_context < REF_CONTEXTS);
    return pred_context;
}

static int32_t get_pred_context_single_ref_p4(PartitionInfo *pi) {
    const uint8_t *const ref_counts = &pi->neighbors_ref_counts[0];

    // Count of LAST
    const int32_t last_count = ref_counts[LAST_FRAME];
    // Count of LAST2
    const int32_t last2_count = ref_counts[LAST2_FRAME];

    const int32_t pred_context =
        (last_count == last2_count) ? 1 : ((last_count < last2_count) ? 0 : 2);

    assert(pred_context >= 0 && pred_context < REF_CONTEXTS);
    return pred_context;
}

static int32_t get_pred_context_last3_or_gld(PartitionInfo *pi) {
    const uint8_t *const ref_counts = &pi->neighbors_ref_counts[0];

    // Count of LAST3
    const int32_t last3_count = ref_counts[LAST3_FRAME];
    // Count of GOLDEN
    const int32_t gld_count = ref_counts[GOLDEN_FRAME];

    const int32_t pred_context =
        (last3_count == gld_count) ? 1 : ((last3_count < gld_count) ? 0 : 2);

    assert(pred_context >= 0 && pred_context < REF_CONTEXTS);
    return pred_context;
}

static void read_ref_frames(ParseCtxt *parse_ctxt, PartitionInfo *const pi) {
    SvtReader *         r          = &parse_ctxt->r;
    int                 segment_id = pi->mi->segment_id;
    MvReferenceFrame *  ref_frame  = pi->mi->ref_frame;
    AomCdfProb *        cdf;
    SegmentationParams *seg_params = &parse_ctxt->frame_header->segmentation_params;
    if (pi->mi->skip_mode) {
        ref_frame[0] =
            (MvReferenceFrame)(parse_ctxt->frame_header->skip_mode_params.ref_frame_idx_0);
        ref_frame[1] =
            (MvReferenceFrame)(parse_ctxt->frame_header->skip_mode_params.ref_frame_idx_1);
    } else if (seg_feature_active(seg_params, segment_id, SEG_LVL_REF_FRAME)) {
        ref_frame[0] = (MvReferenceFrame)get_segdata(seg_params, segment_id, SEG_LVL_REF_FRAME);
        ref_frame[1] = NONE_FRAME;
    } else if (seg_feature_active(seg_params, segment_id, SEG_LVL_SKIP) ||
               seg_feature_active(seg_params, segment_id, SEG_LVL_GLOBALMV)) {
        ref_frame[0] = LAST_FRAME;
        ref_frame[1] = NONE_FRAME;
    } else {
        ReferenceMode mode = SINGLE_REFERENCE;
        int           bw4  = mi_size_wide[pi->mi->sb_type];
        int           bh4  = mi_size_high[pi->mi->sb_type];
        if (parse_ctxt->frame_header->reference_mode == REFERENCE_MODE_SELECT &&
            (AOMMIN(bw4, bh4) >= 2)) {
            const int ctx = get_reference_mode_context(pi);
            mode          = (ReferenceMode)svt_read_symbol(
                r, parse_ctxt->cur_tile_ctx.comp_inter_cdf[ctx], 2, ACCT_STR);
        }

        if (mode == COMPOUND_REFERENCE) {
            int                     pred_context;
            const int               ctx           = get_comp_reference_type_context(pi);
            const CompReferenceType comp_ref_type = (CompReferenceType)svt_read_symbol(
                r, parse_ctxt->cur_tile_ctx.comp_ref_type_cdf[ctx], 2, ACCT_STR);

            if (comp_ref_type == UNIDIR_COMP_REFERENCE) {
                pred_context = get_pred_context_uni_comp_ref_p(pi);
                uint16_t bit = (uint16_t)svt_read_symbol(
                    r, parse_ctxt->cur_tile_ctx.uni_comp_ref_cdf[pred_context][0], 2, ACCT_STR);
                if (bit) {
                    ref_frame[0] = BWDREF_FRAME;
                    ref_frame[1] = ALTREF_FRAME;
                } else {
                    pred_context  = get_pred_context_uni_comp_ref_p1(pi);
                    uint16_t bit1 = (uint16_t)svt_read_symbol(
                        r, parse_ctxt->cur_tile_ctx.uni_comp_ref_cdf[pred_context][1], 2, ACCT_STR);
                    if (bit1) {
                        pred_context  = get_pred_context_uni_comp_ref_p2(pi);
                        uint16_t bit2 = (uint16_t)svt_read_symbol(
                            r,
                            parse_ctxt->cur_tile_ctx.uni_comp_ref_cdf[pred_context][2],
                            2,
                            ACCT_STR);
                        if (bit2) {
                            ref_frame[0] = LAST_FRAME;
                            ref_frame[1] = GOLDEN_FRAME;
                        } else {
                            ref_frame[0] = LAST_FRAME;
                            ref_frame[1] = LAST3_FRAME;
                        }
                    } else {
                        ref_frame[0] = LAST_FRAME;
                        ref_frame[1] = LAST2_FRAME;
                    }
                }
                return;
            }

            assert(comp_ref_type == BIDIR_COMP_REFERENCE);

            const int idx = 1;
            pred_context  = get_pred_context_comp_ref_p(pi);
            uint16_t bit  = (uint16_t)svt_read_symbol(
                r, parse_ctxt->cur_tile_ctx.comp_ref_cdf[pred_context][0], 2, ACCT_STR);
            // Decode forward references.
            if (!bit) {
                uint16_t bit1 = (uint16_t)svt_read_symbol(
                    r,
                    parse_ctxt->cur_tile_ctx.comp_ref_cdf[get_pred_context_single_ref_p4(pi)][1],
                    2,
                    ACCT_STR);
                ref_frame[!idx] = bit1 ? LAST2_FRAME : LAST_FRAME;
            } else {
                uint16_t bit2 = (uint16_t)svt_read_symbol(
                    r,
                    parse_ctxt->cur_tile_ctx.comp_ref_cdf[get_pred_context_last3_or_gld(pi)][2],
                    2,
                    ACCT_STR);
                ref_frame[!idx] = bit2 ? GOLDEN_FRAME : LAST3_FRAME;
            }

            // Decode backward references.
            pred_context     = get_pred_context_comp_bwdref_p(pi);
            uint16_t bit_bwd = (uint16_t)svt_read_symbol(
                r, parse_ctxt->cur_tile_ctx.comp_bwdref_cdf[pred_context][0], 2, ACCT_STR);
            if (!bit_bwd) {
                pred_context      = get_pred_context_comp_bwdref_p1(pi);
                uint16_t bit1_bwd = (uint16_t)svt_read_symbol(
                    r, parse_ctxt->cur_tile_ctx.comp_bwdref_cdf[pred_context][1], 2, ACCT_STR);
                ref_frame[idx] = bit1_bwd ? ALTREF2_FRAME : BWDREF_FRAME;
            } else {
                ref_frame[idx] = ALTREF_FRAME;
            }
        } else if (mode == SINGLE_REFERENCE) {
            cdf = parse_ctxt->cur_tile_ctx.single_ref_cdf[get_pred_context_single_ref_p1(pi)][0];
            const int32_t bit0 = svt_read_symbol(r, cdf, 2, ACCT_STR);

            if (bit0) {
                cdf =
                    parse_ctxt->cur_tile_ctx.single_ref_cdf[get_pred_context_comp_bwdref_p(pi)][1];
                const int32_t bit1 = svt_read_symbol(r, cdf, 2, ACCT_STR);
                if (!bit1) {
                    cdf = parse_ctxt->cur_tile_ctx
                              .single_ref_cdf[get_pred_context_comp_bwdref_p1(pi)][5];
                    const int32_t bit5 = svt_read_symbol(r, cdf, 2, ACCT_STR);
                    ref_frame[0]       = bit5 ? ALTREF2_FRAME : BWDREF_FRAME;
                } else {
                    ref_frame[0] = ALTREF_FRAME;
                }
            } else {
                cdf = parse_ctxt->cur_tile_ctx.single_ref_cdf[get_pred_context_comp_ref_p(pi)][2];
                const int32_t bit2 = svt_read_symbol(r, cdf, 2, ACCT_STR);
                if (bit2) {
                    cdf = parse_ctxt->cur_tile_ctx
                              .single_ref_cdf[get_pred_context_last3_or_gld(pi)][4];
                    const int32_t bit4 = svt_read_symbol(r, cdf, 2, ACCT_STR);
                    ref_frame[0]       = bit4 ? GOLDEN_FRAME : LAST3_FRAME;
                } else {
                    cdf = parse_ctxt->cur_tile_ctx
                              .single_ref_cdf[get_pred_context_single_ref_p4(pi)][3];
                    const int32_t bit3 = svt_read_symbol(r, cdf, 2, ACCT_STR);
                    ref_frame[0]       = bit3 ? LAST2_FRAME : LAST_FRAME;
                }
            }

            ref_frame[1] = NONE_FRAME;
        } else
            assert(0 && "Invalid prediction mode.");
    }
}

int has_newmv(PredictionMode mode) {
    return (mode == NEWMV || mode == NEW_NEWMV || mode == NEAR_NEWMV || mode == NEW_NEARMV ||
            mode == NEAREST_NEWMV || mode == NEW_NEARESTMV);
}

static void add_ref_mv_candidate(EbDecHandle *dec_handle, const BlockModeInfo *const candidate,
                                 const MvReferenceFrame rf[2], uint8_t *num_mv_found,
                                 uint8_t *found_match, uint8_t *newmv_count,
                                 CandidateMv *ref_mv_stack, IntMv *gm_mv_candidates, int weight) {
    if (!is_inter_block(candidate)) return; // for intrabc
    int index = 0, ref;
    assert(weight % 2 == 0);

    EbDecPicBuf *buf = dec_handle->cur_pic_buf[0];
    if (rf[1] == NONE_FRAME) {
        // single reference frame
        for (ref = 0; ref < 2; ++ref) {
            if (candidate->ref_frame[ref] == rf[0]) {
                IntMv this_refmv;
                if (is_global_mv_block(
                        candidate->mode, candidate->sb_type, buf->global_motion[rf[0]].gm_type)) {
                    this_refmv = gm_mv_candidates[0];
                } else
                    this_refmv = candidate->mv[ref];

                for (index = 0; index < *num_mv_found; ++index)
                    if (ref_mv_stack[index].this_mv.as_int == this_refmv.as_int) break;

                if (index < *num_mv_found) ref_mv_stack[index].weight += weight;

                // Add a new item to the list.
                if (index == *num_mv_found && *num_mv_found < MAX_REF_MV_STACK_SIZE) {
                    ref_mv_stack[index].this_mv.as_int = this_refmv.as_int;
                    ref_mv_stack[index].weight         = weight;
                    ++(*num_mv_found);
                }
                if (has_newmv(candidate->mode)) ++*newmv_count;
                ++*found_match;
            }
        }
    } else {
        // compound reference frame
        if (candidate->ref_frame[0] == rf[0] && candidate->ref_frame[1] == rf[1]) {
            IntMv this_refmv[2];
            for (ref = 0; ref < 2; ++ref) {
                if (is_global_mv_block(
                        candidate->mode, candidate->sb_type, buf->global_motion[rf[ref]].gm_type)) {
                    this_refmv[ref] = gm_mv_candidates[ref];
                } else
                    this_refmv[ref] = candidate->mv[ref];
            }

            //*found_match = 1;

            for (index = 0; index < *num_mv_found; ++index)
                if ((ref_mv_stack[index].this_mv.as_int == this_refmv[0].as_int) &&
                    (ref_mv_stack[index].comp_mv.as_int == this_refmv[1].as_int))
                    break;

            if (index < *num_mv_found) ref_mv_stack[index].weight += weight;

            // Add a new item to the list.
            if (index == *num_mv_found && *num_mv_found < MAX_REF_MV_STACK_SIZE) {
                ref_mv_stack[index].this_mv.as_int = this_refmv[0].as_int;
                ref_mv_stack[index].comp_mv.as_int = this_refmv[1].as_int;
                ref_mv_stack[index].weight         = weight;
                ++(*num_mv_found);
            }
            if (has_newmv(candidate->mode)) ++*newmv_count;
            ++*found_match;
        }
    }
}

static void scan_row_mbmi(EbDecHandle *dec_handle, ParseCtxt *parse_ctx, PartitionInfo *pi,
                          int delta_row, const MvReferenceFrame rf[2], CandidateMv *ref_mv_stack,
                          uint8_t *num_mv_found, uint8_t *found_match, uint8_t *newmv_count,
                          IntMv *gm_mv_candidates, int max_row_offset, int *processed_rows) {
    int mi_row = pi->mi_row;
    int mi_col = pi->mi_col;

    int          bw4         = mi_size_wide[pi->mi->sb_type];
    FrameHeader *frm_header  = &dec_handle->frame_header;
    int          end4        = AOMMIN(AOMMIN(bw4, (int)frm_header->mi_cols - mi_col), 16);
    int          delta_col   = 0;
    int          use_step_16 = (bw4 >= 16);
    const int    n8_w_8      = mi_size_wide[BLOCK_8X8];
    const int    n8_w_16     = mi_size_wide[BLOCK_16X16];

    if (abs(delta_row) > 1) {
        delta_col = 1;
        if ((mi_col & 0x01) && bw4 < n8_w_8) --delta_col;
    }

    for (int i = 0; i < end4;) {
        int mv_row = mi_row + delta_row;
        int mv_col = mi_col + delta_col + i;
        if (!is_inside(&parse_ctx->cur_tile_info, mv_col, mv_row)) break;
        BlockModeInfo *candidate = get_cur_mode_info(dec_handle, mv_row, mv_col, pi->sb_info);
        int            len       = AOMMIN(bw4, mi_size_wide[candidate->sb_type]);
        const int      n4_w      = mi_size_wide[candidate->sb_type];
        if (use_step_16)
            len = AOMMAX(n8_w_16, len);
        else if (abs(delta_row) > 1)
            len = AOMMAX(n8_w_8, len);

        int weight = 2;
        if (bw4 >= n8_w_8 && bw4 <= n4_w) {
            int inc = AOMMIN(-max_row_offset + delta_row + 1, mi_size_high[candidate->sb_type]);
            // Obtain range used in weight calculation.
            weight          = AOMMAX(weight, inc);
            *processed_rows = inc - delta_row - 1;
        }
        add_ref_mv_candidate(dec_handle,
                             candidate,
                             rf,
                             num_mv_found,
                             found_match,
                             newmv_count,
                             ref_mv_stack,
                             gm_mv_candidates,
                             len * weight);

        i += len;
    }
}

static void scan_col_mbmi(EbDecHandle *dec_handle, ParseCtxt *parse_ctx, PartitionInfo *pi,
                          int delta_col, const MvReferenceFrame rf[2], CandidateMv *ref_mv_stack,
                          uint8_t *num_mv_found, uint8_t *found_match, uint8_t *newmv_count,
                          IntMv *gm_mv_candidates, int max_col_offset, int *processed_cols) {
    int          mi_row      = pi->mi_row;
    int          mi_col      = pi->mi_col;
    int          bh4         = mi_size_high[pi->mi->sb_type];
    FrameHeader *frm_header  = &dec_handle->frame_header;
    int          end4        = AOMMIN(AOMMIN(bh4, (int)frm_header->mi_rows - mi_row), 16);
    int          delta_row   = 0;
    int          use_step_16 = (bh4 >= 16);
    const int    n8_h_8      = mi_size_high[BLOCK_8X8];

    if (abs(delta_col) > 1) {
        delta_row = 1;
        if ((mi_row & 0x01) && bh4 < n8_h_8) --delta_row;
    }

    for (int i = 0; i < end4;) {
        int mv_row = mi_row + delta_row + i;
        int mv_col = mi_col + delta_col;
        if (!is_inside(&parse_ctx->cur_tile_info, mv_col, mv_row)) break;
        BlockModeInfo *candidate = get_cur_mode_info(dec_handle, mv_row, mv_col, pi->sb_info);
        int            len       = AOMMIN(bh4, mi_size_high[candidate->sb_type]);
        const int      n4_h      = mi_size_high[candidate->sb_type];
        if (abs(delta_col) > 1) len = AOMMAX(2, len);
        if (use_step_16) len = AOMMAX(4, len);

        int weight = 2;
        if (bh4 >= n8_h_8 && bh4 <= n4_h) {
            int inc = AOMMIN(-max_col_offset + delta_col + 1, mi_size_wide[candidate->sb_type]);
            // Obtain range used in weight calculation.
            weight          = AOMMAX(weight, inc);
            *processed_cols = inc - delta_col - 1;
        }

        add_ref_mv_candidate(dec_handle,
                             candidate,
                             rf,
                             num_mv_found,
                             found_match,
                             newmv_count,
                             ref_mv_stack,
                             gm_mv_candidates,
                             len * weight);

        i += len;
    }
}

static void scan_blk_mbmi(EbDecHandle *dec_handle, ParseCtxt *parse_ctx, PartitionInfo *pi,
                          int delta_row, int delta_col, const MvReferenceFrame rf[2],
                          CandidateMv *ref_mv_stack, uint8_t *found_match, uint8_t *newmv_count,
                          IntMv *gm_mv_candidates, uint8_t num_mv_found[MODE_CTX_REF_FRAMES]) {
    int mv_row = pi->mi_row + delta_row;
    int mv_col = pi->mi_col + delta_col;
    int weight = 4;

    if (is_inside(&parse_ctx->cur_tile_info, mv_col, mv_row)) {
        BlockModeInfo *candidate = get_cur_mode_info(dec_handle, mv_row, mv_col, pi->sb_info);

        add_ref_mv_candidate(dec_handle,
                             candidate,
                             rf,
                             num_mv_found,
                             found_match,
                             newmv_count,
                             ref_mv_stack,
                             gm_mv_candidates,
                             weight);
    } // Analyze a single 8x8 block motion information.
}

/* TODO: Harmonize with Encoder. */
static int has_top_right(EbDecHandle *dec_handle, PartitionInfo *pi, int bs) {
    int       n4_w       = mi_size_wide[pi->mi->sb_type];
    int       n4_h       = mi_size_high[pi->mi->sb_type];
    const int sb_mi_size = mi_size_wide[dec_handle->seq_header.sb_size];
    const int mask_row   = pi->mi_row & (sb_mi_size - 1);
    const int mask_col   = pi->mi_col & (sb_mi_size - 1);

    if (bs > mi_size_wide[BLOCK_64X64]) return 0;
    int has_tr = !((mask_row & bs) && (mask_col & bs));

    assert(bs > 0 && !(bs & (bs - 1)));

    while (bs < sb_mi_size) {
        if (mask_col & bs) {
            if ((mask_col & (2 * bs)) && (mask_row & (2 * bs))) {
                has_tr = 0;
                break;
            }
        } else {
            break;
        }
        bs <<= 1;
    }

    if (n4_w < n4_h)
        if (!pi->is_sec_rect) has_tr = 1;
    if (n4_w > n4_h)
        if (pi->is_sec_rect) has_tr = 0;
    if (pi->mi->partition == PARTITION_VERT_A) {
        if (n4_w == n4_h)
            if (mask_row & bs) has_tr = 0;
    }
    return has_tr;
}

static int add_tpl_ref_mv(EbDecHandle *dec_handle, ParseCtxt *parse_ctx, int mi_row, int mi_col,
                          MvReferenceFrame ref_frame, int blk_row, int blk_col,
                          IntMv *gm_mv_candidates, uint8_t *num_mv_found,
                          CandidateMv ref_mv_stacks[][MAX_REF_MV_STACK_SIZE],
                          int16_t *   mode_context) {
    uint8_t      idx;
    FrameHeader *frm_header = &dec_handle->frame_header;
    int          mv_row     = (mi_row + blk_row) | 1;
    int          mv_col     = (mi_col + blk_col) | 1;

    if (!is_inside(&parse_ctx->cur_tile_info, mv_col, mv_row)) return 0;

    int x8 = mv_col >> 1;
    int y8 = mv_row >> 1;

    MvReferenceFrame rf[2];
    av1_set_ref_frame(rf, ref_frame);

    const TemporalMvRef *tpl_mvs =
        dec_handle->master_frame_buf.tpl_mvs + y8 * (frm_header->mi_stride >> 1) + x8;
    const IntMv prev_frame_mvs = tpl_mvs->mf_mv0;
    if (rf[1] == NONE_FRAME) {
        int                      cur_frame_index = dec_handle->cur_pic_buf[0]->order_hint;
        const EbDecPicBuf *const buf_0           = get_ref_frame_buf(dec_handle, rf[0]);
        int                      frame0_index    = buf_0->order_hint;
        int                      cur_offset_0    = get_relative_dist(
            &dec_handle->seq_header.order_hint_info, cur_frame_index, frame0_index);
        CandidateMv *ref_mv_stack = ref_mv_stacks[rf[0]];

        if (prev_frame_mvs.as_int == INVALID_MV) return 0;

        IntMv this_refmv;
        get_mv_projection(
            &this_refmv.as_mv, prev_frame_mvs.as_mv, cur_offset_0, tpl_mvs->ref_frame_offset);

        lower_mv_precision(
            &this_refmv.as_mv, frm_header->allow_high_precision_mv, frm_header->force_integer_mv);

        if (blk_row == 0 && blk_col == 0) {
            if (abs(this_refmv.as_mv.row - gm_mv_candidates[0].as_mv.row) >= 16 ||
                abs(this_refmv.as_mv.col - gm_mv_candidates[0].as_mv.col) >= 16)
                /*zero_mv_ctxt*/
                mode_context[ref_frame] |= (1 << GLOBALMV_OFFSET);
        }

        for (idx = 0; idx < *num_mv_found; ++idx)
            if (this_refmv.as_int == ref_mv_stack[idx].this_mv.as_int) break;

        if (idx < *num_mv_found)
            ref_mv_stack[idx].weight += 2;
        else if (*num_mv_found < MAX_REF_MV_STACK_SIZE) {
            ref_mv_stack[idx].this_mv.as_int = this_refmv.as_int;
            ref_mv_stack[idx].weight         = 2;
            ++(*num_mv_found);
        }
        return 1;
    } else {
        // Process compound inter mode
        int                      cur_frame_index = dec_handle->cur_pic_buf[0]->order_hint;
        const EbDecPicBuf *const buf_0           = get_ref_frame_buf(dec_handle, rf[0]);
        int                      frame0_index = buf_0->order_hint;

        int cur_offset_0 = get_relative_dist(
            &dec_handle->seq_header.order_hint_info, cur_frame_index, frame0_index);
        const EbDecPicBuf *const buf_1 = get_ref_frame_buf(dec_handle, rf[1]);
        int                      frame1_index = buf_1->order_hint;
        int                      cur_offset_1 = get_relative_dist(
            &dec_handle->seq_header.order_hint_info, cur_frame_index, frame1_index);
        CandidateMv *ref_mv_stack = ref_mv_stacks[ref_frame];

        if (prev_frame_mvs.as_int == INVALID_MV) return 0;

        IntMv this_refmv;
        IntMv comp_refmv;
        get_mv_projection(
            &this_refmv.as_mv, prev_frame_mvs.as_mv, cur_offset_0, tpl_mvs->ref_frame_offset);
        get_mv_projection(
            &comp_refmv.as_mv, prev_frame_mvs.as_mv, cur_offset_1, tpl_mvs->ref_frame_offset);

        lower_mv_precision(
            &this_refmv.as_mv, frm_header->allow_high_precision_mv, frm_header->force_integer_mv);
        lower_mv_precision(
            &comp_refmv.as_mv, frm_header->allow_high_precision_mv, frm_header->force_integer_mv);

        if (blk_row == 0 && blk_col == 0) {
            if (abs(this_refmv.as_mv.row - gm_mv_candidates[0].as_mv.row) >= 16 ||
                abs(this_refmv.as_mv.col - gm_mv_candidates[0].as_mv.col) >= 16 ||
                abs(comp_refmv.as_mv.row - gm_mv_candidates[1].as_mv.row) >= 16 ||
                abs(comp_refmv.as_mv.col - gm_mv_candidates[1].as_mv.col) >= 16)
                /*zero_mv_ctxt*/
                mode_context[ref_frame] |= (1 << GLOBALMV_OFFSET);
        }

        for (idx = 0; idx < *num_mv_found; ++idx)
            if (this_refmv.as_int == ref_mv_stack[idx].this_mv.as_int &&
                comp_refmv.as_int == ref_mv_stack[idx].comp_mv.as_int)
                break;

        if (idx < *num_mv_found)
            ref_mv_stack[idx].weight += 2;
        else if (*num_mv_found < MAX_REF_MV_STACK_SIZE) {
            ref_mv_stack[idx].this_mv.as_int = this_refmv.as_int;
            ref_mv_stack[idx].comp_mv.as_int = comp_refmv.as_int;
            ref_mv_stack[idx].weight         = 2;
            ++(*num_mv_found);
        }
    }
    return 1;
}

static void add_extra_mv_candidate(BlockModeInfo *candidate, EbDecHandle *dec_handle,
                                   MvReferenceFrame *rf, IntMv ref_id[2][2], int ref_id_count[2],
                                   IntMv ref_diff[2][2], int ref_diff_count[2]) {
    FrameHeader *frm_header = &dec_handle->frame_header;
    for (int rf_idx = 0; rf_idx < 2; ++rf_idx) {
        MvReferenceFrame can_rf = candidate->ref_frame[rf_idx];
        if (can_rf > INTRA_FRAME) {
            for (int cmp_idx = 0; cmp_idx < 2; ++cmp_idx) {
                if (can_rf == rf[cmp_idx] && ref_id_count[cmp_idx] < 2) {
                    ref_id[cmp_idx][ref_id_count[cmp_idx]] = candidate->mv[rf_idx];
                    ++ref_id_count[cmp_idx];
                } else if (ref_diff_count[cmp_idx] < 2) {
                    IntMv this_mv = candidate->mv[rf_idx];
                    if (frm_header->ref_frame_sign_bias[can_rf] !=
                        frm_header->ref_frame_sign_bias[rf[cmp_idx]]) {
                        this_mv.as_mv.row = -this_mv.as_mv.row;
                        this_mv.as_mv.col = -this_mv.as_mv.col;
                    }
                    ref_diff[cmp_idx][ref_diff_count[cmp_idx]] = this_mv;
                    ++ref_diff_count[cmp_idx];
                }
            }
        }
    }
}

static void process_single_ref_mv_candidate(BlockModeInfo *candidate, EbDecHandle *dec_handle,
                                            MvReferenceFrame ref_frame,
                                            uint8_t          refmv_count[MODE_CTX_REF_FRAMES],
                                            CandidateMv ref_mv_stack[][MAX_REF_MV_STACK_SIZE]) {
    FrameHeader *frm_header = &dec_handle->frame_header;
    for (int rf_idx = 0; rf_idx < 2; ++rf_idx) {
        if (candidate->ref_frame[rf_idx] > INTRA_FRAME) {
            IntMv this_mv = candidate->mv[rf_idx];
            if (frm_header->ref_frame_sign_bias[candidate->ref_frame[rf_idx]] !=
                frm_header->ref_frame_sign_bias[ref_frame]) {
                this_mv.as_mv.row = -this_mv.as_mv.row;
                this_mv.as_mv.col = -this_mv.as_mv.col;
            }
            int stack_idx;
            for (stack_idx = 0; stack_idx < refmv_count[ref_frame]; ++stack_idx) {
                const IntMv stack_mv = ref_mv_stack[ref_frame][stack_idx].this_mv;
                if (this_mv.as_int == stack_mv.as_int) break;
            }

            if (stack_idx == refmv_count[ref_frame]) {
                ref_mv_stack[ref_frame][stack_idx].this_mv = this_mv;
                ref_mv_stack[ref_frame][stack_idx].weight  = 2;
                ++refmv_count[ref_frame];
            }
        }
    }
}

static INLINE void clamp_mv_ref(MV *mv, int bw, int bh, PartitionInfo *pi) {
    clamp_mv(mv,
             pi->mb_to_left_edge - bw * 8 - MV_BORDER,
             pi->mb_to_right_edge + bw * 8 + MV_BORDER,
             pi->mb_to_top_edge - bh * 8 - MV_BORDER,
             pi->mb_to_bottom_edge + bh * 8 + MV_BORDER);
}

static void dec_setup_ref_mv_list(EbDecHandle *dec_handle, ParseCtxt *parse_ctx, PartitionInfo *pi,
                                  MvReferenceFrame ref_frame,
                                  CandidateMv      ref_mv_stack[][MAX_REF_MV_STACK_SIZE],
                                  IntMv            mv_ref_list[][MAX_MV_REF_CANDIDATES],
                                  IntMv *gm_mv_candidates, int16_t *mode_context, MvCount *mv_cnt) {
    int              n4_w = mi_size_wide[pi->mi->sb_type];
    int              n4_h = mi_size_high[pi->mi->sb_type];
    const int        bs   = AOMMAX(n4_w, n4_h);
    MvReferenceFrame rf[2];

    FrameHeader *         frame_info     = parse_ctx->frame_header;
    const TileInfo *const tile           = &parse_ctx->cur_tile_info;
    int                   max_row_offset = 0, max_col_offset = 0;
    int32_t               mi_row         = pi->mi_row;
    int32_t               mi_col         = pi->mi_col;
    const int             row_adj        = (n4_h < mi_size_high[BLOCK_8X8]) && (mi_row & 0x01);
    const int             col_adj        = (n4_w < mi_size_wide[BLOCK_8X8]) && (mi_col & 0x01);
    int                   processed_rows = 0;
    int                   processed_cols = 0;

    av1_set_ref_frame(rf, ref_frame);
    mode_context[ref_frame] = 0;

    // Find valid maximum row/col offset.
    if (pi->up_available) {
        max_row_offset = -(MVREF_ROW_COLS << 1) + row_adj;

        if (n4_h < mi_size_high[BLOCK_8X8]) max_row_offset = -(2 << 1) + row_adj;

        max_row_offset =
            clamp(max_row_offset, tile->mi_row_start - mi_row, tile->mi_row_end - mi_row - 1);
    }

    if (pi->left_available) {
        max_col_offset = -(MVREF_ROW_COLS << 1) + col_adj;

        if (n4_w < mi_size_wide[BLOCK_8X8]) max_col_offset = -(2 << 1) + col_adj;

        max_col_offset =
            clamp(max_col_offset, tile->mi_col_start - mi_col, tile->mi_col_end - mi_col - 1);
    }
    memset(mv_cnt, 0, sizeof(*mv_cnt));

    // Scan the first above row mode info. row_offset = -1;
    if (abs(max_row_offset) >= 1) {
        scan_row_mbmi(dec_handle,
                      parse_ctx,
                      pi,
                      -1,
                      rf,
                      ref_mv_stack[ref_frame],
                      &mv_cnt->num_mv_found[ref_frame],
                      &mv_cnt->found_above_match,
                      &mv_cnt->newmv_count,
                      gm_mv_candidates,
                      max_row_offset,
                      &processed_rows);
    }

    // Scan the first left column mode info. col_offset = -1;
    if (abs(max_col_offset) >= 1) {
        scan_col_mbmi(dec_handle,
                      parse_ctx,
                      pi,
                      -1,
                      rf,
                      ref_mv_stack[ref_frame],
                      &mv_cnt->num_mv_found[ref_frame],
                      &mv_cnt->found_left_match,
                      &mv_cnt->newmv_count,
                      gm_mv_candidates,
                      max_col_offset,
                      &processed_cols);
    }

    if (has_top_right(dec_handle, pi, bs)) {
        scan_blk_mbmi(dec_handle,
                      parse_ctx,
                      pi,
                      -1,
                      n4_w,
                      rf,
                      ref_mv_stack[ref_frame],
                      &mv_cnt->found_above_match,
                      &mv_cnt->newmv_count,
                      gm_mv_candidates,
                      &mv_cnt->num_mv_found[ref_frame]);
    }

    const uint8_t nearest_match = (mv_cnt->found_above_match > 0) + (mv_cnt->found_left_match > 0);
    const uint8_t num_nearest   = mv_cnt->num_mv_found[ref_frame];
    const uint8_t num_new       = mv_cnt->newmv_count;

    for (int idx = 0; idx < num_nearest; ++idx)
        ref_mv_stack[ref_frame][idx].weight += REF_CAT_LEVEL;

    if (frame_info->use_ref_frame_mvs) {
        int       is_available = 0;
        const int voffset      = AOMMAX(mi_size_high[BLOCK_8X8], n4_h);
        const int hoffset      = AOMMAX(mi_size_wide[BLOCK_8X8], n4_w);
        const int blk_row_end  = AOMMIN(n4_h, mi_size_high[BLOCK_64X64]);
        const int blk_col_end  = AOMMIN(n4_w, mi_size_wide[BLOCK_64X64]);

        const int tpl_sample_pos[3][2] = {
            {voffset, -2},
            {voffset, hoffset},
            {voffset - 2, hoffset},
        };
        const int allow_extension =
            (n4_h >= mi_size_high[BLOCK_8X8]) && (n4_h < mi_size_high[BLOCK_64X64]) &&
            (n4_w >= mi_size_wide[BLOCK_8X8]) && (n4_w < mi_size_wide[BLOCK_64X64]);

        const int step_h = (n4_h >= mi_size_high[BLOCK_64X64]) ? mi_size_high[BLOCK_16X16]
                                                               : mi_size_high[BLOCK_8X8];
        const int step_w = (n4_w >= mi_size_wide[BLOCK_64X64]) ? mi_size_wide[BLOCK_16X16]
                                                               : mi_size_wide[BLOCK_8X8];

        for (int blk_row = 0; blk_row < blk_row_end; blk_row += step_h) {
            for (int blk_col = 0; blk_col < blk_col_end; blk_col += step_w) {
                int ret = add_tpl_ref_mv(dec_handle,
                                         parse_ctx,
                                         mi_row,
                                         mi_col,
                                         ref_frame,
                                         blk_row,
                                         blk_col,
                                         gm_mv_candidates,
                                         &mv_cnt->num_mv_found[ref_frame],
                                         ref_mv_stack,
                                         mode_context);
                if (blk_row == 0 && blk_col == 0) is_available = ret;
            }
        }

        if (is_available == 0) mode_context[ref_frame] |= (1 << GLOBALMV_OFFSET);

        if (allow_extension) {
            for (int i = 0; i < 3; ++i) {
                const int blk_row = tpl_sample_pos[i][0];
                const int blk_col = tpl_sample_pos[i][1];

                if (check_sb_border(mi_row, mi_col, blk_row, blk_col)) {
                    add_tpl_ref_mv(dec_handle,
                                   parse_ctx,
                                   mi_row,
                                   mi_col,
                                   ref_frame,
                                   blk_row,
                                   blk_col,
                                   gm_mv_candidates,
                                   &mv_cnt->num_mv_found[ref_frame],
                                   ref_mv_stack,
                                   mode_context);
                }
            }
        }
    }

    // Scan the second outer area.
    scan_blk_mbmi(dec_handle,
                  parse_ctx,
                  pi,
                  -1,
                  -1,
                  rf,
                  ref_mv_stack[ref_frame],
                  &mv_cnt->found_above_match,
                  &mv_cnt->newmv_count,
                  gm_mv_candidates,
                  &mv_cnt->num_mv_found[ref_frame]);

    for (int idx = 2; idx <= MVREF_ROW_COLS; ++idx) {
        const int row_offset = -(idx << 1) + 1 + row_adj;
        const int col_offset = -(idx << 1) + 1 + col_adj;
        if (abs(row_offset) <= abs(max_row_offset) && abs(row_offset) > processed_rows) {
            scan_row_mbmi(dec_handle,
                          parse_ctx,
                          pi,
                          row_offset,
                          rf,
                          ref_mv_stack[ref_frame],
                          &mv_cnt->num_mv_found[ref_frame],
                          &mv_cnt->found_above_match,
                          &mv_cnt->newmv_count,
                          gm_mv_candidates,
                          max_row_offset,
                          &processed_rows);
        }

        if (abs(col_offset) <= abs(max_col_offset) && abs(col_offset) > processed_cols) {
            scan_col_mbmi(dec_handle,
                          parse_ctx,
                          pi,
                          col_offset,
                          rf,
                          ref_mv_stack[ref_frame],
                          &mv_cnt->num_mv_found[ref_frame],
                          &mv_cnt->found_left_match,
                          &mv_cnt->newmv_count,
                          gm_mv_candidates,
                          max_col_offset,
                          &processed_cols);
        }
    }

    /* sorting process*/
    int start = 0;
    int end   = num_nearest;
    while (end > start) {
        int new_end = start;
        for (int idx = start + 1; idx < end; ++idx) {
            if (ref_mv_stack[ref_frame][idx - 1].weight < ref_mv_stack[ref_frame][idx].weight) {
                CandidateMv tmp_mv               = ref_mv_stack[ref_frame][idx - 1];
                ref_mv_stack[ref_frame][idx - 1] = ref_mv_stack[ref_frame][idx];
                ref_mv_stack[ref_frame][idx]     = tmp_mv;
                new_end                          = idx;
            }
        }
        end = new_end;
    }

    start = num_nearest;
    end   = mv_cnt->num_mv_found[ref_frame];
    while (end > start) {
        int new_end = start;
        for (int idx = start + 1; idx < end; ++idx) {
            if (ref_mv_stack[ref_frame][idx - 1].weight < ref_mv_stack[ref_frame][idx].weight) {
                CandidateMv tmp_mv               = ref_mv_stack[ref_frame][idx - 1];
                ref_mv_stack[ref_frame][idx - 1] = ref_mv_stack[ref_frame][idx];
                ref_mv_stack[ref_frame][idx]     = tmp_mv;
                new_end                          = idx;
            }
        }
        end = new_end;
    }

    /* extra search process */
    if (mv_cnt->num_mv_found[ref_frame] < MAX_MV_REF_CANDIDATES) {
        IntMv ref_id[2][2], ref_diff[2][2];
        int   ref_id_count[2] = {0}, ref_diff_count[2] = {0};

        int mi_width  = AOMMIN(16, n4_w);
        mi_width      = AOMMIN(mi_width, (int)frame_info->mi_cols - mi_col);
        int mi_height = AOMMIN(16, n4_h);
        mi_height     = AOMMIN(mi_height, (int)frame_info->mi_rows - mi_row);
        int mi_size   = AOMMIN(mi_width, mi_height);

        for (int pass = 0; pass < 2; pass++) {
            int idx = 0;
            while (idx < mi_size && mv_cnt->num_mv_found[ref_frame] < MAX_MV_REF_CANDIDATES) {
                int mv_row, mv_col;
                if (pass == 0) {
                    mv_row = mi_row - 1;
                    mv_col = mi_col + idx;
                } else {
                    mv_row = mi_row + idx;
                    mv_col = mi_col - 1;
                }

                if (!is_inside(&parse_ctx->cur_tile_info, mv_col, mv_row)) break;

                BlockModeInfo *nbr = get_cur_mode_info(dec_handle, mv_row, mv_col, pi->sb_info);

                if (rf[1] != NONE_FRAME)
                    add_extra_mv_candidate(
                        nbr, dec_handle, rf, ref_id, ref_id_count, ref_diff, ref_diff_count);
                else
                    process_single_ref_mv_candidate(
                        nbr, dec_handle, ref_frame, mv_cnt->num_mv_found, ref_mv_stack);

                idx += pass ? mi_size_high[nbr->sb_type] : mi_size_wide[nbr->sb_type];
            }
        }

        if (rf[1] > NONE_FRAME) {
            IntMv comp_list[3][2];

            for (int idx = 0; idx < 2; ++idx) {
                int comp_idx = 0;
                for (int list_idx = 0; list_idx < ref_id_count[idx]; ++list_idx, comp_idx++) {
                    comp_list[comp_idx][idx] = ref_id[idx][list_idx];
                }

                for (int list_idx = 0; list_idx < ref_diff_count[idx] && comp_idx < 2;
                     ++list_idx, ++comp_idx) {
                    comp_list[comp_idx][idx] = ref_diff[idx][list_idx];
                }

                for (; comp_idx < 2; ++comp_idx) comp_list[comp_idx][idx] = gm_mv_candidates[idx];
            }

            if (mv_cnt->num_mv_found[ref_frame]) {
                assert(mv_cnt->num_mv_found[ref_frame] == 1);
                if (comp_list[0][0].as_int == ref_mv_stack[ref_frame][0].this_mv.as_int &&
                    comp_list[0][1].as_int == ref_mv_stack[ref_frame][0].comp_mv.as_int) {
                    ref_mv_stack[ref_frame][mv_cnt->num_mv_found[ref_frame]].this_mv =
                        comp_list[1][0];
                    ref_mv_stack[ref_frame][mv_cnt->num_mv_found[ref_frame]].comp_mv =
                        comp_list[1][1];
                } else {
                    ref_mv_stack[ref_frame][mv_cnt->num_mv_found[ref_frame]].this_mv =
                        comp_list[0][0];
                    ref_mv_stack[ref_frame][mv_cnt->num_mv_found[ref_frame]].comp_mv =
                        comp_list[0][1];
                }
                ref_mv_stack[ref_frame][mv_cnt->num_mv_found[ref_frame]].weight = 2;
                ++mv_cnt->num_mv_found[ref_frame];
            } else {
                for (int idx = 0; idx < MAX_MV_REF_CANDIDATES; ++idx) {
                    ref_mv_stack[ref_frame][mv_cnt->num_mv_found[ref_frame]].this_mv =
                        comp_list[idx][0];
                    ref_mv_stack[ref_frame][mv_cnt->num_mv_found[ref_frame]].comp_mv =
                        comp_list[idx][1];
                    ref_mv_stack[ref_frame][mv_cnt->num_mv_found[ref_frame]].weight = 2;
                    ++mv_cnt->num_mv_found[ref_frame];
                }
            }
        }
    }

    /* context and clamping process */
    //int num_lists = rf[1] > NONE_FRAME ? 2 : 1;
    //for (int list = 0; list < num_lists; list++) {
    //    for (int idx = 0; idx < mv_cnt->num_mv_found[ref_frame]; idx++) {
    //        IntMv refMv = ref_mv_stack[ref_frame][idx].this_mv;
    //        refMv.as_mv.row = clamp_mv_row(dec_handle, pi, mi_row,
    //            refMv.as_mv.row, MV_BORDER + n4_h * 8);
    //        refMv.as_mv.col = clamp_mv_col(dec_handle, pi, mi_col,
    //            refMv.as_mv.col, MV_BORDER + n4_w * 8);
    //        ref_mv_stack[ref_frame][idx].this_mv = refMv;
    //    }
    //}

    if (rf[1] > NONE_FRAME) {
        for (int idx = 0; idx < mv_cnt->num_mv_found[ref_frame]; ++idx) {
            clamp_mv_ref(&ref_mv_stack[ref_frame][idx].this_mv.as_mv,
                         n4_w << MI_SIZE_LOG2,
                         n4_h << MI_SIZE_LOG2,
                         pi);
            clamp_mv_ref(&ref_mv_stack[ref_frame][idx].comp_mv.as_mv,
                         n4_w << MI_SIZE_LOG2,
                         n4_h << MI_SIZE_LOG2,
                         pi);
        }
    } else {
        for (int idx = 0; idx < mv_cnt->num_mv_found[ref_frame]; ++idx) {
            clamp_mv_ref(&ref_mv_stack[ref_frame][idx].this_mv.as_mv,
                         n4_w << MI_SIZE_LOG2,
                         n4_h << MI_SIZE_LOG2,
                         pi);
        }
    }

    const uint8_t ref_match_count =
        (mv_cnt->found_above_match > 0) + (mv_cnt->found_left_match > 0);
    switch (nearest_match) {
    case 0:
        mode_context[ref_frame] |= 0;
        if (ref_match_count >= 1) mode_context[ref_frame] |= 1;
        if (ref_match_count == 1)
            mode_context[ref_frame] |= (1 << REFMV_OFFSET);
        else if (ref_match_count >= 2)
            mode_context[ref_frame] |= (2 << REFMV_OFFSET);
        break;
    case 1:
        mode_context[ref_frame] |= (num_new > 0) ? 2 : 3;
        if (ref_match_count == 1)
            mode_context[ref_frame] |= (3 << REFMV_OFFSET);
        else if (ref_match_count >= 2)
            mode_context[ref_frame] |= (4 << REFMV_OFFSET);
        break;
    case 2:
    default:
        if (num_new >= 1)
            mode_context[ref_frame] |= 4;
        else
            mode_context[ref_frame] |= 5;

        mode_context[ref_frame] |= (5 << REFMV_OFFSET);
        break;
    }

    if (rf[1] == NONE_FRAME && mv_ref_list != NULL) {
        for (int idx = mv_cnt->num_mv_found[ref_frame]; idx < MAX_MV_REF_CANDIDATES; ++idx)
            mv_ref_list[rf[0]][idx].as_int = gm_mv_candidates[0].as_int;

        for (int idx = 0; idx < AOMMIN(MAX_MV_REF_CANDIDATES, mv_cnt->num_mv_found[ref_frame]);
             ++idx) {
            mv_ref_list[rf[0]][idx].as_int = ref_mv_stack[ref_frame][idx].this_mv.as_int;
        }
    }
}

static INLINE int16_t svt_mode_context_analyzer(const int16_t *const          mode_context,
                                                const MvReferenceFrame *const rf) {
    const int8_t ref_frame = av1_ref_frame_type(rf);

    if (rf[1] <= INTRA_FRAME) return mode_context[ref_frame];

    const int16_t newmv_ctx = mode_context[ref_frame] & NEWMV_CTX_MASK;
    const int16_t refmv_ctx = (mode_context[ref_frame] >> REFMV_OFFSET) & REFMV_CTX_MASK;

    const int16_t comp_ctx =
        compound_mode_ctx_map[refmv_ctx >> 1][AOMMIN(newmv_ctx, COMP_NEWMV_CTXS - 1)];
    return comp_ctx;
}

void av1_find_mv_refs(EbDecHandle *dec_handle, PartitionInfo *pi, ParseCtxt *parse_ctx,
                      MvReferenceFrame ref_frame, CandidateMv ref_mv_stack[][MAX_REF_MV_STACK_SIZE],
                      IntMv mv_ref_list[][MAX_MV_REF_CANDIDATES], IntMv global_mvs[2],
                      int16_t *mode_context, MvCount *mv_cnt) {
    BlockSize        bsize = pi->mi->sb_type;
    MvReferenceFrame rf[2];
    av1_set_ref_frame(rf, ref_frame);

    /* setup global mv process */
    global_mvs[0].as_int = 0;
    global_mvs[1].as_int = 0;
    if (ref_frame != INTRA_FRAME) {
        EbDecPicBuf *buf = dec_handle->cur_pic_buf[0];
        global_mvs[0].as_int =
            gm_get_motion_vector(&buf->global_motion[rf[0]],
                                 dec_handle->frame_header.allow_high_precision_mv,
                                 bsize,
                                 pi->mi_col,
                                 pi->mi_row,
                                 dec_handle->frame_header.force_integer_mv)
                .as_int;

        global_mvs[1].as_int =
            (rf[1] != NONE_FRAME)
                ? gm_get_motion_vector(&buf->global_motion[rf[1]],
                                       dec_handle->frame_header.allow_high_precision_mv,
                                       bsize,
                                       pi->mi_col,
                                       pi->mi_row,
                                       dec_handle->frame_header.force_integer_mv)
                      .as_int
                : 0;
    }
    dec_setup_ref_mv_list(dec_handle,
                          parse_ctx,
                          pi,
                          ref_frame,
                          ref_mv_stack,
                          mv_ref_list,
                          global_mvs,
                          mode_context,
                          mv_cnt);
}

static PredictionMode read_inter_compound_mode(ParseCtxt *parse_ctxt, int16_t ctx) {
    SvtReader *r    = &parse_ctxt->r;
    const int  mode = svt_read_symbol(
        r, parse_ctxt->cur_tile_ctx.inter_compound_mode_cdf[ctx], INTER_COMPOUND_MODES, ACCT_STR);
    assert(is_inter_compound_mode(NEAREST_NEARESTMV + mode));
    return NEAREST_NEARESTMV + mode;
}

static INLINE int has_nearmv(PredictionMode mode) {
    return (mode == NEARMV || mode == NEAR_NEARMV || mode == NEAR_NEWMV || mode == NEW_NEARMV);
}

static INLINE uint8_t get_drl_ctx(const CandidateMv *ref_mv_stack, int ref_idx) {
    if (ref_mv_stack[ref_idx].weight >= REF_CAT_LEVEL &&
        ref_mv_stack[ref_idx + 1].weight < REF_CAT_LEVEL) {
        return 1;
    }

    if (ref_mv_stack[ref_idx].weight < REF_CAT_LEVEL &&
        ref_mv_stack[ref_idx + 1].weight < REF_CAT_LEVEL) {
        return 2;
    }

    return 0;
}

static void read_drl_idx(ParseCtxt *parse_ctxt, PartitionInfo *pi, BlockModeInfo *mbmi,
                         int num_mv_found) {
    SvtReader *r              = &parse_ctxt->r;
    uint8_t    ref_frame_type = av1_ref_frame_type(mbmi->ref_frame);
    mbmi->ref_mv_idx          = 0;
    if (mbmi->mode == NEWMV || mbmi->mode == NEW_NEWMV) {
        for (int idx = 0; idx < 2; ++idx) {
            if (num_mv_found > idx + 1) {
                uint8_t drl_ctx = get_drl_ctx(pi->ref_mv_stack[ref_frame_type], idx);
                int     drl_idx =
                    svt_read_symbol(r, parse_ctxt->cur_tile_ctx.drl_cdf[drl_ctx], 2, ACCT_STR);
                mbmi->ref_mv_idx = idx;
                if (!drl_idx) return;
                mbmi->ref_mv_idx = idx + 1;
            }
        }
    }
    if (have_nearmv_in_inter_mode(mbmi->mode)) {
        for (int idx = 1; idx < 3; ++idx) {
            if (num_mv_found > idx + 1) {
                uint8_t drl_ctx = get_drl_ctx(pi->ref_mv_stack[ref_frame_type], idx);
                int     drl_idx =
                    svt_read_symbol(r, parse_ctxt->cur_tile_ctx.drl_cdf[drl_ctx], 2, ACCT_STR);
                mbmi->ref_mv_idx = idx + drl_idx - 1;
                if (!drl_idx) return;
            }
        }
    }
}

/* TODO: Harmonize*/
static void svt_find_best_ref_mvs(int allow_hp, IntMv *mvlist, IntMv *nearest_mv, IntMv *near_mv,
                                  int is_integer) {
    int i;
    // Make sure all the candidates are properly clamped etc
    for (i = 0; i < MAX_MV_REF_CANDIDATES; ++i) {
        lower_mv_precision(&mvlist[i].as_mv, allow_hp, is_integer);
    }
    *nearest_mv = mvlist[0];
    *near_mv    = mvlist[1];
}

int read_mv_component(SvtReader *r, NmvComponent *mvcomp, int use_subpel, int usehp) {
    int       mag, d, fr, hp;
    const int sign     = svt_read_symbol(r, mvcomp->sign_cdf, 2, ACCT_STR);
    const int mv_class = svt_read_symbol(r, mvcomp->classes_cdf, MV_CLASSES, ACCT_STR);
    const int class0   = mv_class == MV_CLASS_0;

    // Integer part
    if (class0) {
        d   = svt_read_symbol(r, mvcomp->class0_cdf, CLASS0_SIZE, ACCT_STR);
        mag = 0;
    } else {
        d = 0;
        for (int i = 0; i < mv_class; ++i)
            d |= svt_read_symbol(r, mvcomp->bits_cdf[i], 2, ACCT_STR) << i;
        mag = CLASS0_SIZE << (mv_class + 2);
    }

    fr = use_subpel
             ? svt_read_symbol(
                   r, class0 ? mvcomp->class0_fp_cdf[d] : mvcomp->fp_cdf, MV_FP_SIZE, ACCT_STR)
             : 3;

    hp = usehp ? svt_read_symbol(r, class0 ? mvcomp->class0_hp_cdf : mvcomp->hp_cdf, 2, ACCT_STR)
               : 1;

    // Result
    mag += ((d << 3) | (fr << 1) | hp) + 1;
    return sign ? -mag : mag;
}

static INLINE void read_mv(SvtReader *r, MV *mv, MV *ref, NmvContext *ctx,
                           MvSubpelPrecision precision) {
    MV diff = k_zero_mv;

    const MvJointType joint_type =
        (MvJointType)svt_read_symbol(r, ctx->joints_cdf, MV_JOINTS, ACCT_STR);

    if (mv_joint_vertical(joint_type))
        diff.row = read_mv_component(
            r, &ctx->comps[0], precision > MV_SUBPEL_NONE, precision > MV_SUBPEL_LOW_PRECISION);

    if (mv_joint_horizontal(joint_type))
        diff.col = read_mv_component(
            r, &ctx->comps[1], precision > MV_SUBPEL_NONE, precision > MV_SUBPEL_LOW_PRECISION);

    mv->row = ref->row + diff.row;
    mv->col = ref->col + diff.col;
}

static INLINE int assign_mv(ParseCtxt *parse_ctxt, PartitionInfo *pi, IntMv mv[2],
                            IntMv *global_mvs, IntMv ref_mv[2], IntMv nearest_mv[2],
                            IntMv near_mv[2], int is_compound, int allow_hp) {
    SvtReader *    r    = &parse_ctxt->r;
    BlockModeInfo *mbmi = pi->mi;

    if (parse_ctxt->frame_header->force_integer_mv) allow_hp = MV_SUBPEL_NONE;

    switch (mbmi->mode) {
    case NEWMV: {
        NmvContext *const nmvc = &parse_ctxt->cur_tile_ctx.nmvc;
        read_mv(r, &mv[0].as_mv, &ref_mv[0].as_mv, nmvc, allow_hp);
        break;
    }
    case NEARESTMV: {
        mv[0].as_int = nearest_mv[0].as_int;
        break;
    }
    case NEARMV: {
        mv[0].as_int = near_mv[0].as_int;
        break;
    }
    case GLOBALMV: {
        mv[0].as_int = global_mvs[0].as_int;
        break;
    }
    case NEW_NEWMV: {
        assert(is_compound);
        for (int i = 0; i < 2; ++i) {
            NmvContext *const nmvc = &parse_ctxt->cur_tile_ctx.nmvc;
            read_mv(r, &mv[i].as_mv, &ref_mv[i].as_mv, nmvc, allow_hp);
        }
        break;
    }
    case NEAREST_NEARESTMV: {
        assert(is_compound);
        mv[0].as_int = nearest_mv[0].as_int;
        mv[1].as_int = nearest_mv[1].as_int;
        break;
    }
    case NEAR_NEARMV: {
        assert(is_compound);
        mv[0].as_int = near_mv[0].as_int;
        mv[1].as_int = near_mv[1].as_int;
        break;
    }
    case NEW_NEARESTMV: {
        NmvContext *const nmvc = &parse_ctxt->cur_tile_ctx.nmvc;
        read_mv(r, &mv[0].as_mv, &ref_mv[0].as_mv, nmvc, allow_hp);
        assert(is_compound);
        mv[1].as_int = nearest_mv[1].as_int;
        break;
    }
    case NEAREST_NEWMV: {
        mv[0].as_int           = nearest_mv[0].as_int;
        NmvContext *const nmvc = &parse_ctxt->cur_tile_ctx.nmvc;
        read_mv(r, &mv[1].as_mv, &ref_mv[1].as_mv, nmvc, allow_hp);
        assert(is_compound);
        break;
    }
    case NEAR_NEWMV: {
        mv[0].as_int           = near_mv[0].as_int;
        NmvContext *const nmvc = &parse_ctxt->cur_tile_ctx.nmvc;
        read_mv(r, &mv[1].as_mv, &ref_mv[1].as_mv, nmvc, allow_hp);
        assert(is_compound);
        break;
    }
    case NEW_NEARMV: {
        NmvContext *const nmvc = &parse_ctxt->cur_tile_ctx.nmvc;
        read_mv(r, &mv[0].as_mv, &ref_mv[0].as_mv, nmvc, allow_hp);
        assert(is_compound);
        mv[1].as_int = near_mv[1].as_int;
        break;
    }
    case GLOBAL_GLOBALMV: {
        assert(is_compound);
        mv[0].as_int = global_mvs[0].as_int;
        mv[1].as_int = global_mvs[1].as_int;
        break;
    }
    default: {
        return 0;
    }
    }

    int ret = is_mv_valid(&mv[0].as_mv);
    if (is_compound) ret = ret && is_mv_valid(&mv[1].as_mv);
    return ret;
}

static INLINE int is_dv_valid(MV dv, ParseCtxt *parse_ctx, PartitionInfo *pi) {
    int       mi_row         = pi->mi_row;
    int       mi_col         = pi->mi_col;
    int       mib_size_log2  = parse_ctx->seq_header->sb_size_log2;
    int       subsampling_x  = parse_ctx->seq_header->color_config.subsampling_x;
    int       subsampling_y  = parse_ctx->seq_header->color_config.subsampling_y;
    BlockSize bsize          = pi->mi->sb_type;
    const int bw             = block_size_wide[bsize];
    const int bh             = block_size_high[bsize];
    const int scale_px_to_mv = 8;
    if (((dv.row & (scale_px_to_mv - 1)) || (dv.col & (scale_px_to_mv - 1)))) return 0;

    TileInfo *tile          = &parse_ctx->cur_tile_info;
    const int src_top_edge  = mi_row * MI_SIZE * scale_px_to_mv + dv.row;
    const int tile_top_edge = tile->mi_row_start * MI_SIZE * scale_px_to_mv;
    if (src_top_edge < tile_top_edge) return 0;
    const int src_left_edge  = mi_col * MI_SIZE * scale_px_to_mv + dv.col;
    const int tile_left_edge = tile->mi_col_start * MI_SIZE * scale_px_to_mv;
    if (src_left_edge < tile_left_edge) return 0;
    const int src_bottom_edge  = (mi_row * MI_SIZE + bh) * scale_px_to_mv + dv.row;
    const int tile_bottom_edge = tile->mi_row_end * MI_SIZE * scale_px_to_mv;
    if (src_bottom_edge > tile_bottom_edge) return 0;
    const int src_right_edge  = (mi_col * MI_SIZE + bw) * scale_px_to_mv + dv.col;
    const int tile_right_edge = tile->mi_col_end * MI_SIZE * scale_px_to_mv;
    if (src_right_edge > tile_right_edge) return 0;

    // Special case for sub 8x8 chroma cases, to prevent referring to chroma
    // pixels outside current tile.
    int     num_planes    = parse_ctx->seq_header->color_config.mono_chrome ? 1 : MAX_MB_PLANE;
    int32_t is_chroma_ref = pi->is_chroma_ref;
    for (int plane = 1; plane < num_planes; ++plane) {
        if (is_chroma_ref) {
            if (bw < 8 && subsampling_x)
                if (src_left_edge < tile_left_edge + 4 * scale_px_to_mv) return 0;
            if (bh < 8 && subsampling_y)
                if (src_top_edge < tile_top_edge + 4 * scale_px_to_mv) return 0;
        }
    }

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

    return 1;
}

int dec_assign_dv(ParseCtxt *parse_ctxt, PartitionInfo *pi, IntMv *mv, IntMv *ref_mv) {
    SvtReader *    r       = &parse_ctxt->r;
    FRAME_CONTEXT *frm_ctx = &parse_ctxt->cur_tile_ctx;
    read_mv(r, &mv->as_mv, &ref_mv->as_mv, &frm_ctx->ndvc, MV_SUBPEL_NONE);
    // DV should not have sub-pel.
    assert((mv->as_mv.col & 7) == 0);
    assert((mv->as_mv.row & 7) == 0);
    mv->as_mv.col = (mv->as_mv.col >> 3) * 8;
    mv->as_mv.row = (mv->as_mv.row >> 3) * 8;
    int valid     = is_mv_valid(&mv->as_mv) && is_dv_valid(mv->as_mv, parse_ctxt, pi);
    return valid;
}

void assign_intrabc_mv(ParseCtxt *parse_ctxt, IntMv ref_mvs[INTRA_FRAME + 1][MAX_MV_REF_CANDIDATES],
                       PartitionInfo *pi) {
    BlockModeInfo *mbmi = pi->mi;
    IntMv          nearestmv, nearmv;
    svt_find_best_ref_mvs(0, ref_mvs[INTRA_FRAME], &nearestmv, &nearmv, 0);
    IntMv dv_ref = nearestmv.as_int == 0 ? nearmv : nearestmv;
    if (dv_ref.as_int == 0) {
        av1_find_ref_dv(&dv_ref,
                        &parse_ctxt->cur_tile_info,
                        parse_ctxt->seq_header->sb_mi_size,
                        pi->mi_row,
                        pi->mi_col);
    }
    // Ref DV should not have sub-pel.
    int valid_dv     = (dv_ref.as_mv.col & 7) == 0 && (dv_ref.as_mv.row & 7) == 0;
    dv_ref.as_mv.col = (dv_ref.as_mv.col >> 3) * 8;
    dv_ref.as_mv.row = (dv_ref.as_mv.row >> 3) * 8;
    valid_dv         = valid_dv && dec_assign_dv(parse_ctxt, pi, &mbmi->mv[0], &dv_ref);
}

void read_interintra_mode(ParseCtxt *parse_ctxt, BlockModeInfo *mbmi) {
    SvtReader *    r       = &parse_ctxt->r;
    FRAME_CONTEXT *frm_ctx = &parse_ctxt->cur_tile_ctx;
    BlockSize      bsize   = mbmi->sb_type;
    if (parse_ctxt->seq_header->enable_interintra_compound && !mbmi->skip_mode &&
        is_interintra_allowed(mbmi)) {
        const int bsize_group = size_group_lookup[bsize];
        mbmi->is_inter_intra =
            svt_read_symbol(r, frm_ctx->interintra_cdf[bsize_group], 2, ACCT_STR);
        assert(mbmi->ref_frame[1] == NONE_FRAME);
        if (mbmi->is_inter_intra) {
            mbmi->interintra_mode_params.interintra_mode = (InterIntraMode)svt_read_symbol(
                r, frm_ctx->interintra_mode_cdf[bsize_group], INTERINTRA_MODES, ACCT_STR);
            mbmi->ref_frame[1]                            = INTRA_FRAME;
            mbmi->angle_delta[PLANE_TYPE_Y]               = 0;
            mbmi->angle_delta[PLANE_TYPE_UV]              = 0;
            mbmi->filter_intra_mode_info.use_filter_intra = 0;
            if (is_interintra_wedge_used(bsize)) {
                mbmi->interintra_mode_params.wedge_interintra =
                    svt_read_symbol(r, frm_ctx->wedge_interintra_cdf[bsize], 2, ACCT_STR);
                if (mbmi->interintra_mode_params.wedge_interintra) {
                    mbmi->interintra_mode_params.interintra_wedge_index =
                        svt_read_symbol(r, frm_ctx->wedge_idx_cdf[bsize], 16, ACCT_STR);
                }
            }
        }
    }
}

static INLINE void add_samples(BlockModeInfo *mbmi, int *pts, int *pts_inref, int row_offset,
                               int sign_r, int col_offset, int sign_c) {
    int bw = block_size_wide[mbmi->sb_type];
    int bh = block_size_high[mbmi->sb_type];
    int x  = col_offset * MI_SIZE + sign_c * AOMMAX(bw, MI_SIZE) / 2 - 1;
    int y  = row_offset * MI_SIZE + sign_r * AOMMAX(bh, MI_SIZE) / 2 - 1;

    pts[0]       = (x * 8);
    pts[1]       = (y * 8);
    pts_inref[0] = (x * 8) + mbmi->mv[0].as_mv.col;
    pts_inref[1] = (y * 8) + mbmi->mv[0].as_mv.row;
}

int find_warp_samples(EbDecHandle *dec_handle, TileInfo *tile, PartitionInfo *pi, int *pts,
                      int *pts_inref) {
    int mi_row = pi->mi_row;
    int mi_col = pi->mi_col;

    BlockModeInfo *const mbmi0          = pi->mi;
    int                  ref_frame      = mbmi0->ref_frame[0];
    int                  up_available   = pi->up_available;
    int                  left_available = pi->left_available;
    int                  i, mi_step = 1, np = 0;

    int do_tl = 1;
    int do_tr = 1;
    int b4_w  = mi_size_wide[pi->mi->sb_type];
    int b4_h  = mi_size_high[pi->mi->sb_type];

    // scan the nearest above rows
    if (up_available) {
        BlockModeInfo *mbmi = get_cur_mode_info(dec_handle, mi_row - 1, mi_col, pi->sb_info);
        uint8_t        n4_w = mi_size_wide[mbmi->sb_type];

        if (b4_w <= n4_w) {
            // Handle "current block width <= above block width" case.
            int col_offset = -mi_col % n4_w;

            if (col_offset < 0) do_tl = 0;
            if (col_offset + n4_w > b4_w) do_tr = 0;

            if (mbmi->ref_frame[0] == ref_frame && mbmi->ref_frame[1] == NONE_FRAME) {
                add_samples(mbmi, pts, pts_inref, 0, -1, col_offset, 1);
                pts += 2;
                pts_inref += 2;
                np++;
                if (np >= LEAST_SQUARES_SAMPLES_MAX) return LEAST_SQUARES_SAMPLES_MAX;
            }
        } else {
            // Handle "current block width > above block width" case.
            for (i = 0; i < AOMMIN(b4_w, tile->mi_col_end - mi_col); i += mi_step) {
                mbmi    = get_cur_mode_info(dec_handle, mi_row - 1, mi_col + i, pi->sb_info);
                n4_w    = mi_size_wide[mbmi->sb_type];
                mi_step = AOMMIN(b4_w, n4_w);

                if (mbmi->ref_frame[0] == ref_frame && mbmi->ref_frame[1] == NONE_FRAME) {
                    add_samples(mbmi, pts, pts_inref, 0, -1, i, 1);
                    pts += 2;
                    pts_inref += 2;
                    np++;
                    if (np >= LEAST_SQUARES_SAMPLES_MAX) return LEAST_SQUARES_SAMPLES_MAX;
                }
            }
        }
    }
    assert(np <= LEAST_SQUARES_SAMPLES_MAX);

    // scan the nearest left columns
    if (left_available) {
        BlockModeInfo *mbmi = get_cur_mode_info(dec_handle, mi_row, mi_col - 1, pi->sb_info);
        uint8_t        n4_h = mi_size_high[mbmi->sb_type];

        if (b4_h <= n4_h) {
            // Handle "current block height <= above block height" case.
            int row_offset = -mi_row % n4_h;

            if (row_offset < 0) do_tl = 0;

            if (mbmi->ref_frame[0] == ref_frame && mbmi->ref_frame[1] == NONE_FRAME) {
                add_samples(mbmi, pts, pts_inref, row_offset, 1, 0, -1);
                pts += 2;
                pts_inref += 2;
                np++;
                if (np >= LEAST_SQUARES_SAMPLES_MAX) return LEAST_SQUARES_SAMPLES_MAX;
            }
        } else {
            // Handle "current block height > above block height" case.
            for (i = 0; i < AOMMIN(b4_h, tile->mi_row_end - mi_row); i += mi_step) {
                mbmi    = get_cur_mode_info(dec_handle, mi_row + i, mi_col - 1, pi->sb_info);
                n4_h    = mi_size_high[mbmi->sb_type];
                mi_step = AOMMIN(b4_h, n4_h);

                if (mbmi->ref_frame[0] == ref_frame && mbmi->ref_frame[1] == NONE_FRAME) {
                    add_samples(mbmi, pts, pts_inref, i, 1, 0, -1);
                    pts += 2;
                    pts_inref += 2;
                    np++;
                    if (np >= LEAST_SQUARES_SAMPLES_MAX) return LEAST_SQUARES_SAMPLES_MAX;
                }
            }
        }
    }
    assert(np <= LEAST_SQUARES_SAMPLES_MAX);

    // Top-left block
    if (do_tl && left_available && up_available) {
        BlockModeInfo *mbmi = get_cur_mode_info(dec_handle, mi_row - 1, mi_col - 1, pi->sb_info);

        if (mbmi->ref_frame[0] == ref_frame && mbmi->ref_frame[1] == NONE_FRAME) {
            add_samples(mbmi, pts, pts_inref, 0, -1, 0, -1);
            pts += 2;
            pts_inref += 2;
            np++;
            if (np >= LEAST_SQUARES_SAMPLES_MAX) return LEAST_SQUARES_SAMPLES_MAX;
        }
    }
    assert(np <= LEAST_SQUARES_SAMPLES_MAX);

    // Top-right block
    if (do_tr && has_top_right(dec_handle, pi, AOMMAX(b4_w, b4_h))) {
        int mv_row = mi_row - 1;
        int mv_col = mi_col + b4_w;

        if (is_inside(tile, mv_col, mv_row)) {
            BlockModeInfo *mbmi = get_cur_mode_info(dec_handle, mv_row, mv_col, pi->sb_info);

            if (mbmi->ref_frame[0] == ref_frame && mbmi->ref_frame[1] == NONE_FRAME) {
                add_samples(mbmi, pts, pts_inref, 0, -1, b4_w, 1);
                np++;
                if (np >= LEAST_SQUARES_SAMPLES_MAX) return LEAST_SQUARES_SAMPLES_MAX;
            }
        }
    }
    assert(np <= LEAST_SQUARES_SAMPLES_MAX);

    return np;
}

int has_overlappable_cand(EbDecHandle *dec_handle, ParseCtxt *parse_ctx, PartitionInfo *pi) {
    int                   mi_row = pi->mi_row;
    int                   mi_col = pi->mi_col;
    const TileInfo *const tile   = &parse_ctx->cur_tile_info;
    BlockModeInfo *       mbmi   = pi->mi;
    if (!is_motion_variation_allowed_bsize(mbmi->sb_type)) return 0;

    if (pi->up_available) {
        int w4 = mi_size_wide[mbmi->sb_type];
        int x4 = mi_col;
        while (x4 < AOMMIN(tile->mi_col_end, mi_col + w4)) {
            BlockModeInfo *top_nb_mode =
                get_cur_mode_info(dec_handle, mi_row - 1, x4 | 1, pi->sb_info);
            x4 += AOMMAX(2, mi_size_wide[top_nb_mode->sb_type] >> 2);
            if (is_inter_block(top_nb_mode)) return 1;
        }
    }
    if (pi->left_available) {
        int h4 = mi_size_high[mbmi->sb_type];
        int y4 = mi_row;
        while (y4 < AOMMIN(tile->mi_row_end, mi_row + h4)) {
            BlockModeInfo *left_nb_mode =
                get_cur_mode_info(dec_handle, y4 | 1, mi_col - 1, pi->sb_info);
            y4 += AOMMAX(2, mi_size_high[left_nb_mode->sb_type] >> 2);
            if (is_inter_block(left_nb_mode)) return 1;
        }
    }
    return 0;
}

static INLINE MotionMode is_motion_mode_allowed(EbDecHandle *dec_handle, ParseCtxt *parse_ctx,
                                                GlobalMotionParams *gm_params, PartitionInfo *pi,
                                                int allow_warped_motion) {
    BlockModeInfo *mbmi = pi->mi;
    if (dec_handle->frame_header.force_integer_mv == 0) {
        const TransformationType gm_type = gm_params[mbmi->ref_frame[0]].gm_type;
        if (is_global_mv_block(mbmi->mode, mbmi->sb_type, gm_type)) return SIMPLE_TRANSLATION;
    }
    if ((block_size_wide[mbmi->sb_type] >= 8 && block_size_high[mbmi->sb_type] >= 8) &&
        (mbmi->mode >= NEARESTMV && mbmi->mode < MB_MODE_COUNT) &&
        mbmi->ref_frame[1] != INTRA_FRAME && !has_second_ref(mbmi)) {
        if (!has_overlappable_cand(dec_handle, parse_ctx, pi)) return SIMPLE_TRANSLATION;
        assert(!has_second_ref(mbmi));

        if (pi->num_samples >= 1 && (allow_warped_motion && !av1_is_scaled(pi->block_ref_sf[0]))) {
            if (dec_handle->frame_header.force_integer_mv) { return OBMC_CAUSAL; }
            return WARPED_CAUSAL;
        }
        return OBMC_CAUSAL;
    } else {
        return SIMPLE_TRANSLATION;
    }
}

MotionMode read_motion_mode(EbDecHandle *dec_handle, ParseCtxt *parse_ctxt, PartitionInfo *pi) {
    SvtReader *    r                   = &parse_ctxt->r;
    FRAME_CONTEXT *frm_ctx             = &parse_ctxt->cur_tile_ctx;
    FrameHeader *  frame_info          = &dec_handle->frame_header;
    int            allow_warped_motion = frame_info->allow_warped_motion;
    BlockModeInfo *mbmi                = pi->mi;

    if (dec_handle->frame_header.is_motion_mode_switchable == 0) return SIMPLE_TRANSLATION;
    if (mbmi->skip_mode) return SIMPLE_TRANSLATION;

    const MotionMode last_motion_mode_allowed = is_motion_mode_allowed(
        dec_handle, parse_ctxt, dec_handle->cur_pic_buf[0]->global_motion, pi, allow_warped_motion);
    int motion_mode;

    if (last_motion_mode_allowed == SIMPLE_TRANSLATION) return SIMPLE_TRANSLATION;

    if (last_motion_mode_allowed == OBMC_CAUSAL) {
        motion_mode = svt_read_symbol(r, frm_ctx->obmc_cdf[mbmi->sb_type], 2, ACCT_STR);
        return (MotionMode)(motion_mode);
    } else {
        motion_mode =
            svt_read_symbol(r, frm_ctx->motion_mode_cdf[mbmi->sb_type], MOTION_MODES, ACCT_STR);
        return (MotionMode)(motion_mode);
    }
}

static INLINE int get_comp_group_idx_context(ParseCtxt *parse_ctxt, const PartitionInfo *xd) {
    const BlockModeInfo *const above_mi  = xd->above_mbmi;
    const BlockModeInfo *const left_mi   = xd->left_mbmi;
    int                        above_ctx = 0, left_ctx = 0;

    if (above_mi) {
        if (has_second_ref(above_mi)) {
            above_ctx =
                parse_ctxt->parse_above_nbr4x4_ctxt
                    ->above_comp_grp_idx[xd->mi_col - parse_ctxt->cur_tile_info.mi_col_start];
        } else if (above_mi->ref_frame[0] == ALTREF_FRAME)
            above_ctx = 3;
    }
    if (left_mi) {
        if (has_second_ref(left_mi)) {
            left_ctx = parse_ctxt->parse_left_nbr4x4_ctxt
                           ->left_comp_grp_idx[xd->mi_row - parse_ctxt->sb_row_mi];
        } else if (left_mi->ref_frame[0] == ALTREF_FRAME)
            left_ctx = 3;
    }

    return AOMMIN(5, above_ctx + left_ctx);
}

int get_comp_index_context(EbDecHandle *dec_handle, PartitionInfo *pi) {
    BlockModeInfo *mbmi       = pi->mi;
    SeqHeader *    seq_params = &dec_handle->seq_header;
    FrameHeader *  frm_header = &dec_handle->frame_header;

    int bck_frame_index = 0, fwd_frame_index = 0;
    int cur_frame_index = frm_header->order_hint;

    EbDecPicBuf *bck_buf = get_ref_frame_buf(dec_handle, mbmi->ref_frame[0]);
    EbDecPicBuf *fwd_buf = get_ref_frame_buf(dec_handle, mbmi->ref_frame[1]);

    if (bck_buf != NULL) bck_frame_index = bck_buf->order_hint;
    if (fwd_buf != NULL) fwd_frame_index = fwd_buf->order_hint;

    int fwd =
        abs(get_relative_dist(&seq_params->order_hint_info, fwd_frame_index, cur_frame_index));
    int bck =
        abs(get_relative_dist(&seq_params->order_hint_info, cur_frame_index, bck_frame_index));

    const BlockModeInfo *const above_mi = pi->above_mbmi;
    const BlockModeInfo *const left_mi  = pi->left_mbmi;

    int       above_ctx = 0, left_ctx = 0;
    const int offset = (fwd == bck);

    if (above_mi != NULL) {
        if (has_second_ref(above_mi))
            above_ctx = above_mi->compound_idx;
        else if (above_mi->ref_frame[0] == ALTREF_FRAME)
            above_ctx = 1;
    }

    if (left_mi != NULL) {
        if (has_second_ref(left_mi))
            left_ctx = left_mi->compound_idx;
        else if (left_mi->ref_frame[0] == ALTREF_FRAME)
            left_ctx = 1;
    }

    return above_ctx + left_ctx + 3 * offset;
}

void update_compound_ctx(ParseCtxt *parse_ctxt, PartitionInfo *pi, uint32_t blk_row,
                         uint32_t blk_col, uint32_t comp_grp_idx) {
    ParseAboveNbr4x4Ctxt *above_parse_ctx = parse_ctxt->parse_above_nbr4x4_ctxt;
    ParseLeftNbr4x4Ctxt * left_parse_ctx  = parse_ctxt->parse_left_nbr4x4_ctxt;

    const uint32_t bw = mi_size_wide[pi->mi->sb_type];
    const uint32_t bh = mi_size_high[pi->mi->sb_type];

    int8_t *above_ctx =
        above_parse_ctx->above_comp_grp_idx + blk_col - parse_ctxt->cur_tile_info.mi_col_start;
    int8_t *left_ctx =
        left_parse_ctx->left_comp_grp_idx + ((blk_row - parse_ctxt->sb_row_mi) & MAX_MIB_MASK);

    memset(above_ctx, comp_grp_idx, bw);
    memset(left_ctx, comp_grp_idx, bh);
}

void read_compound_type(EbDecHandle *dec_handle, ParseCtxt *parse_ctxt, PartitionInfo *pi) {
    SvtReader *    r              = &parse_ctxt->r;
    BlockModeInfo *mbmi           = pi->mi;
    BlockSize      bsize          = mbmi->sb_type;
    int32_t        comp_group_idx = 0;
    mbmi->compound_idx            = 1;
    FRAME_CONTEXT *frm_ctx        = &parse_ctxt->cur_tile_ctx;

    if (mbmi->skip_mode) mbmi->inter_inter_compound.type = COMPOUND_AVERAGE;

    if (has_second_ref(mbmi) && !mbmi->skip_mode) {
        // Read idx to indicate current compound inter prediction mode group
        const int masked_compound_used =
            is_any_masked_compound_used(bsize) && dec_handle->seq_header.enable_masked_compound;

        if (masked_compound_used) {
            const int ctx_comp_group_idx = get_comp_group_idx_context(parse_ctxt, pi);
            comp_group_idx =
                svt_read_symbol(r, frm_ctx->comp_group_idx_cdf[ctx_comp_group_idx], 2, ACCT_STR);
        }

        if (comp_group_idx == 0) {
            if (dec_handle->seq_header.order_hint_info.enable_jnt_comp) {
                const int comp_index_ctx = get_comp_index_context(dec_handle, pi);
                mbmi->compound_idx =
                    svt_read_symbol(r, frm_ctx->compound_index_cdf[comp_index_ctx], 2, ACCT_STR);
                mbmi->inter_inter_compound.type =
                    mbmi->compound_idx ? COMPOUND_AVERAGE : COMPOUND_DISTWTD;
            } else {
                // Distance-weighted compound is disabled, so always use average
                mbmi->compound_idx              = 1;
                mbmi->inter_inter_compound.type = COMPOUND_AVERAGE;
            }
        } else {
            assert(dec_handle->frame_header.reference_mode != SINGLE_REFERENCE &&
                   is_inter_compound_mode(mbmi->mode) && mbmi->motion_mode == SIMPLE_TRANSLATION);
            assert(masked_compound_used);

            // compound_diffwtd, wedge
            if (is_interinter_compound_used(COMPOUND_WEDGE, bsize))
                mbmi->inter_inter_compound.type =
                    COMPOUND_WEDGE +
                    svt_read_symbol(
                        r, frm_ctx->compound_type_cdf[bsize], MASKED_COMPOUND_TYPES, ACCT_STR);
            else
                mbmi->inter_inter_compound.type = COMPOUND_DIFFWTD;

            if (mbmi->inter_inter_compound.type == COMPOUND_WEDGE) {
                assert(is_interinter_compound_used(COMPOUND_WEDGE, bsize));
                mbmi->inter_inter_compound.wedge_index =
                    svt_read_symbol(r, frm_ctx->wedge_idx_cdf[bsize], 16, ACCT_STR);
                mbmi->inter_inter_compound.wedge_sign = svt_read_bit(r, ACCT_STR);
            } else {
                assert(mbmi->inter_inter_compound.type == COMPOUND_DIFFWTD);
                mbmi->inter_inter_compound.mask_type =
                    svt_read_literal(r, MAX_DIFFWTD_MASK_BITS, ACCT_STR);
            }
        }
    }

    update_compound_ctx(parse_ctxt, pi, pi->mi_row, pi->mi_col, comp_group_idx);
}

static INLINE int is_nontrans_global_motion(PartitionInfo *pi, GlobalMotionParams *gm_params) {
    int            ref;
    BlockModeInfo *mbmi = pi->mi;
    // First check if all modes are GLOBALMV
    if (mbmi->mode != GLOBALMV && mbmi->mode != GLOBAL_GLOBALMV) return 0;

    if (AOMMIN(mi_size_wide[mbmi->sb_type], mi_size_high[mbmi->sb_type]) < 2) return 0;

    // Now check if all global motion is non translational
    for (ref = 0; ref < 1 + has_second_ref(mbmi); ++ref) {
        if (gm_params[mbmi->ref_frame[ref]].gm_type == TRANSLATION) return 0;
    }
    return 1;
}

static INLINE int av1_is_interp_needed(PartitionInfo *pi, GlobalMotionParams *gm_params) {
    BlockModeInfo *mbmi = pi->mi;
    if (mbmi->skip_mode) return 0;
    if (mbmi->motion_mode == WARPED_CAUSAL) return 0;
    if (is_nontrans_global_motion(pi, gm_params)) return 0;
    return 1;
}

static InterpFilter get_ref_filter_type(const BlockModeInfo *ref_mbmi, int dir,
                                        MvReferenceFrame ref_frame) {
    return ((ref_mbmi->ref_frame[0] == ref_frame || ref_mbmi->ref_frame[1] == ref_frame)
                ? av1_extract_interp_filter(ref_mbmi->interp_filters, dir & 0x01)
                : SWITCHABLE_FILTERS);
}

int get_context_interp(PartitionInfo *pi, int dir) {
    const BlockModeInfo *const mbmi = pi->mi;
    const int ctx_offset            = (mbmi->ref_frame[1] > INTRA_FRAME) * INTER_FILTER_COMP_OFFSET;
    assert(dir == 0 || dir == 1);
    const MvReferenceFrame ref_frame = mbmi->ref_frame[0];

    int filter_type_ctx = ctx_offset + (dir & 0x01) * INTER_FILTER_DIR_OFFSET;
    int left_type       = SWITCHABLE_FILTERS;
    int above_type      = SWITCHABLE_FILTERS;

    if (pi->left_available) left_type = get_ref_filter_type(pi->left_mbmi, dir, ref_frame);

    if (pi->up_available) above_type = get_ref_filter_type(pi->above_mbmi, dir, ref_frame);

    if (left_type == above_type) {
        filter_type_ctx += left_type;
    } else if (left_type == SWITCHABLE_FILTERS) {
        assert(above_type != SWITCHABLE_FILTERS);
        filter_type_ctx += above_type;
    } else if (above_type == SWITCHABLE_FILTERS) {
        assert(left_type != SWITCHABLE_FILTERS);
        filter_type_ctx += left_type;
    } else {
        filter_type_ctx += SWITCHABLE_FILTERS;
    }

    return filter_type_ctx;
}

void inter_block_mode_info(EbDecHandle *dec_handle, ParseCtxt *parse_ctxt, PartitionInfo *pi) {
    BlockModeInfo *     mbmi     = pi->mi;
    SvtReader *         r        = &parse_ctxt->r;
    const int           allow_hp = dec_handle->frame_header.allow_high_precision_mv;
    IntMv               ref_mvs[MODE_CTX_REF_FRAMES][MAX_MV_REF_CANDIDATES] = {{{0}}};
    int16_t             inter_mode_ctx[MODE_CTX_REF_FRAMES];
    int                 pts[SAMPLES_ARRAY_SIZE], pts_inref[SAMPLES_ARRAY_SIZE];
    SegmentationParams *seg = &dec_handle->frame_header.segmentation_params;
    MvCount             mv_cnt;

    mbmi->palette_size[0] = 0;
    mbmi->palette_size[1] = 0;

    /* TO-DO initialize palette info */

    svt_collect_neighbors_ref_counts(pi);

    read_ref_frames(parse_ctxt, pi);
    /* if ((pi->mi->ref_frame[0] >= BWDREF_FRAME && pi->mi->ref_frame[0] <= ALTREF_FRAME) ||
        (pi->mi->ref_frame[1] >= BWDREF_FRAME && pi->mi->ref_frame[1] <= ALTREF_FRAME)) {
        SVT_LOG("ALTREF found - frame : %d\n", dec_handle->dec_cnt);
        exit(0);
    }*/
    const int is_compound = has_second_ref(mbmi);

    MvReferenceFrame ref_frame = av1_ref_frame_type(mbmi->ref_frame);
    IntMv            global_mvs[2];
    av1_find_mv_refs(dec_handle,
                     pi,
                     parse_ctxt,
                     ref_frame,
                     pi->ref_mv_stack,
                     ref_mvs,
                     global_mvs,
                     inter_mode_ctx,
                     &mv_cnt);

#if EXTRA_DUMP
    if (enable_dump) {
        SVT_LOG("\n mi_row: %d mi_col: %d\n", pi->mi_row, pi->mi_col);
        /*for (int i = 0; i < MODE_CTX_REF_FRAMES; i++)
            for (int j = 0; j < MAX_REF_MV_STACK_SIZE; j++)
                SVT_LOG("ref_mv_stack[%d][%d]=%d\t", i, j, pi->ref_mv_stack[i][j].this_mv.as_int);
        SVT_LOG("\n");*/
        fflush(stdout);
    }
#endif

    int mode_ctx     = svt_mode_context_analyzer(inter_mode_ctx, mbmi->ref_frame);
    mbmi->ref_mv_idx = 0;

    if (mbmi->skip_mode) {
        assert(is_compound);
        mbmi->mode = NEAREST_NEARESTMV;
    } else {
        if (seg_feature_active(seg, mbmi->segment_id, SEG_LVL_SKIP) ||
            seg_feature_active(seg, mbmi->segment_id, SEG_LVL_GLOBALMV))
            mbmi->mode = GLOBALMV;
        else {
            if (is_compound)
                mbmi->mode = read_inter_compound_mode(parse_ctxt, mode_ctx);
            else {
                int new_mv = svt_read_symbol(
                    r, parse_ctxt->cur_tile_ctx.newmv_cdf[mode_ctx & NEWMV_CTX_MASK], 2, ACCT_STR);
                if (new_mv) {
                    int zero_mv = svt_read_symbol(
                        r,
                        parse_ctxt->cur_tile_ctx
                            .zeromv_cdf[(mode_ctx >> GLOBALMV_OFFSET) & GLOBALMV_CTX_MASK],
                        2,
                        ACCT_STR);
                    if (zero_mv) {
                        int ref_mv = svt_read_symbol(
                            r,
                            parse_ctxt->cur_tile_ctx
                                .refmv_cdf[(mode_ctx >> REFMV_OFFSET) & REFMV_CTX_MASK],
                            2,
                            ACCT_STR);
                        mbmi->mode = ref_mv ? NEARMV : NEARESTMV;
                    } else
                        mbmi->mode = GLOBALMV;
                } else
                    mbmi->mode = NEWMV;
            }
            if (mbmi->mode == NEWMV || mbmi->mode == NEW_NEWMV || has_nearmv(mbmi->mode))
                read_drl_idx(parse_ctxt, pi, mbmi, mv_cnt.num_mv_found[ref_frame]);
        }
    }
    mbmi->uv_mode = UV_DC_PRED;

    IntMv ref_mv[2];
    IntMv nearestmv[2], nearmv[2];
    memset(nearestmv, 0, sizeof(nearestmv));
    if (!is_compound && mbmi->mode != GLOBALMV) {
        svt_find_best_ref_mvs(allow_hp,
                              ref_mvs[mbmi->ref_frame[0]],
                              &nearestmv[0],
                              &nearmv[0],
                              dec_handle->frame_header.force_integer_mv);
    }
    if (is_compound && mbmi->mode != GLOBAL_GLOBALMV) {
        int ref_mv_idx = mbmi->ref_mv_idx + 1;
        nearestmv[0]   = pi->ref_mv_stack[ref_frame][0].this_mv;
        nearestmv[1]   = pi->ref_mv_stack[ref_frame][0].comp_mv;
        nearmv[0]      = pi->ref_mv_stack[ref_frame][ref_mv_idx].this_mv;
        nearmv[1]      = pi->ref_mv_stack[ref_frame][ref_mv_idx].comp_mv;
        lower_mv_precision(
            &nearestmv[0].as_mv, allow_hp, dec_handle->frame_header.force_integer_mv);
        lower_mv_precision(
            &nearestmv[1].as_mv, allow_hp, dec_handle->frame_header.force_integer_mv);
        lower_mv_precision(&nearmv[0].as_mv, allow_hp, dec_handle->frame_header.force_integer_mv);
        lower_mv_precision(&nearmv[1].as_mv, allow_hp, dec_handle->frame_header.force_integer_mv);
    } else if (mbmi->ref_mv_idx > 0 && mbmi->mode == NEARMV) {
        IntMv cur_mv = pi->ref_mv_stack[mbmi->ref_frame[0]][1 + mbmi->ref_mv_idx].this_mv;
        nearmv[0]    = cur_mv;
    }

    ref_mv[0] = nearestmv[0];
    ref_mv[1] = nearestmv[1];

    if (is_compound) {
        int ref_mv_idx = mbmi->ref_mv_idx;
        // Special case: NEAR_NEWMV and NEW_NEARMV modes use
        // 1 + mbmi->ref_mv_idx (like NEARMV) instead of
        // mbmi->ref_mv_idx (like NEWMV)
        if (mbmi->mode == NEAR_NEWMV || mbmi->mode == NEW_NEARMV) ref_mv_idx = 1 + mbmi->ref_mv_idx;

        if (compound_ref0_mode(mbmi->mode) == NEWMV)
            ref_mv[0] = pi->ref_mv_stack[ref_frame][ref_mv_idx].this_mv;

        if (compound_ref1_mode(mbmi->mode) == NEWMV)
            ref_mv[1] = pi->ref_mv_stack[ref_frame][ref_mv_idx].comp_mv;
    } else {
        if (mbmi->mode == NEWMV) {
            if (mv_cnt.num_mv_found[ref_frame] > 1)
                ref_mv[0] = pi->ref_mv_stack[ref_frame][mbmi->ref_mv_idx].this_mv;
        }
    }

    assign_mv(
        parse_ctxt, pi, mbmi->mv, global_mvs, ref_mv, nearestmv, nearmv, is_compound, allow_hp);

#if EXTRA_DUMP
    if (enable_dump) {
        SVT_LOG("\n mode %d MV %d %d \n", mbmi->mode, mbmi->mv[0].as_mv.row, mbmi->mv[0].as_mv.col);
        fflush(stdout);
    }
#endif
    read_interintra_mode(parse_ctxt, mbmi);

    for (int ref = 0; ref < 1 + has_second_ref(mbmi); ++ref) {
        const MvReferenceFrame frame = mbmi->ref_frame[ref];
        pi->block_ref_sf[ref]        = get_ref_scale_factors(dec_handle, frame);
    }

    pi->num_samples = find_warp_samples(dec_handle, &parse_ctxt->cur_tile_info, pi, pts, pts_inref);

    mbmi->motion_mode = read_motion_mode(dec_handle, parse_ctxt, pi);

    read_compound_type(dec_handle, parse_ctxt, pi);

    if (!av1_is_interp_needed(pi, dec_handle->cur_pic_buf[0]->global_motion)) {
        mbmi->interp_filters = av1_broadcast_interp_filter(
            av1_unswitchable_filter(dec_handle->frame_header.interpolation_filter));
    } else {
        if (dec_handle->frame_header.interpolation_filter != SWITCHABLE) {
            mbmi->interp_filters =
                av1_broadcast_interp_filter(dec_handle->frame_header.interpolation_filter);
        } else {
            InterpFilter ref0_filter[2] = {EIGHTTAP_REGULAR, EIGHTTAP_REGULAR};
            for (int dir = 0; dir < 2; ++dir) {
                const int ctx    = get_context_interp(pi, dir);
                ref0_filter[dir] = (InterpFilter)svt_read_symbol(
                    r,
                    parse_ctxt->cur_tile_ctx.switchable_interp_cdf[ctx],
                    SWITCHABLE_FILTERS,
                    ACCT_STR);
                if (dec_handle->seq_header.enable_dual_filter == 0) {
                    ref0_filter[1] = ref0_filter[0];
                    break;
                }
            }
            mbmi->interp_filters = av1_make_interp_filters(ref0_filter[0], ref0_filter[1]);
        }
    }
}

int get_palette_color_context(uint8_t (*color_map)[COLOR_MAP_STRIDE][COLOR_MAP_STRIDE], int r,
                              int c, int palette_size, uint8_t *color_order) {
    // Get color indices of neighbors.
    int color_neighbors[NUM_PALETTE_NEIGHBORS];
    color_neighbors[0] = (c - 1 >= 0) ? (*color_map)[r][c - 1] : -1;
    color_neighbors[1] = (c - 1 >= 0 && r - 1 >= 0) ? (*color_map)[(r - 1)][c - 1] : -1;
    color_neighbors[2] = (r - 1 >= 0) ? (*color_map)[(r - 1)][c] : -1;

    int              scores[PALETTE_MAX_SIZE + 10] = {0};
    int              i;
    static const int weights[NUM_PALETTE_NEIGHBORS] = {2, 1, 2};
    for (i = 0; i < NUM_PALETTE_NEIGHBORS; ++i) {
        if (color_neighbors[i] >= 0) { scores[color_neighbors[i]] += weights[i]; }
    }

    for (i = 0; i < PALETTE_MAX_SIZE; ++i) { color_order[i] = i; }

    for (i = 0; i < NUM_PALETTE_NEIGHBORS; i++) {
        int max_score = scores[i];
        int max_id    = i;
        for (int j = i + 1; j < palette_size; j++) {
            if (scores[j] > max_score) {
                max_score = scores[j];
                max_id    = j;
            }
        }
        if (max_id != i) {
            max_score           = scores[max_id];
            int max_color_order = color_order[max_id];
            for (int k = max_id; k > i; k--) {
                scores[k]      = scores[k - 1];
                color_order[k] = color_order[k - 1];
            }
            scores[i]      = max_score;
            color_order[i] = max_color_order;
        }
    }
    int              color_index_ctx_hash                    = 0;
    static const int hash_multipliers[NUM_PALETTE_NEIGHBORS] = {1, 2, 2};
    for (int i = 0; i < NUM_PALETTE_NEIGHBORS; i++) {
        color_index_ctx_hash += scores[i] * hash_multipliers[i];
    }
    assert(color_index_ctx_hash > 0);
    assert(color_index_ctx_hash <= MAX_COLOR_CONTEXT_HASH);

    const int color_index_ctx = palette_color_index_context_lookup[color_index_ctx_hash];
    assert(color_index_ctx >= 0);
    assert(color_index_ctx < PALETTE_COLOR_INDEX_CONTEXTS);
    return color_index_ctx;
}

void palette_tokens(EbDecHandle *dec_handle, ParseCtxt *parse_ctx, PartitionInfo *pi) {
    int            mi_row           = pi->mi_row;
    int            mi_col           = pi->mi_col;
    BlockModeInfo *mbmi             = pi->mi;
    BlockSize      bsize            = mbmi->sb_type;
    FRAME_CONTEXT *frm_ctx          = &parse_ctx->cur_tile_ctx;
    SvtReader *    r                = &parse_ctx->r;
    int            block_height     = block_size_high[bsize];
    int            block_width      = block_size_wide[bsize];
    int            mi_cols          = (&dec_handle->frame_header)->mi_cols;
    int            mi_rows          = (&dec_handle->frame_header)->mi_rows;
    int            on_screen_height = MIN(block_height, (mi_rows - mi_row) * MI_SIZE);
    int            on_screen_width  = MIN(block_width, (mi_cols - mi_col) * MI_SIZE);

    int32_t is_chroma_ref = pi->is_chroma_ref;
    uint8_t color_order[PALETTE_MAX_SIZE];
    uint8_t color_map[COLOR_MAP_STRIDE][COLOR_MAP_STRIDE];
    int     sub_x, sub_y;
    for (int plane_itr = 0; plane_itr < MAX_MB_PLANE; plane_itr++) {
        uint8_t palette_size = mbmi->palette_size[plane_itr != 0];
        sub_x                = plane_itr ? dec_handle->seq_header.color_config.subsampling_x : 0;
        sub_y                = plane_itr ? dec_handle->seq_header.color_config.subsampling_y : 0;
        if (plane_itr < PLANE_TYPES && palette_size) {
            block_height     = block_height >> sub_y;
            block_width      = block_width >> sub_x;
            on_screen_height = on_screen_height >> sub_y;
            on_screen_width  = on_screen_width >> sub_x;

            if (plane_itr) {
                if (block_width < 4) {
                    block_width += 2;
                    on_screen_width += 2;
                }
                if (block_height < 4) {
                    block_height += 2;
                    on_screen_height += 2;
                }
            }

            if ((plane_itr ? is_chroma_ref : 1)) {
                int color_index_map = svt_read_ns_ae(r, palette_size, ACCT_STR);
                color_map[0][0]     = color_index_map;
                for (int i = 1; i < on_screen_height + on_screen_width - 1; i++) {
                    for (int j = MIN(i, on_screen_width - 1); j >= MAX(0, i - on_screen_height + 1);
                         j--) {
                        int color_ctx = get_palette_color_context(
                            &color_map, (i - j), j, palette_size, color_order);
                        int palette_color_idx = svt_read_symbol(
                            r,
                            plane_itr
                                ? frm_ctx->palette_uv_color_index_cdf[palette_size -
                                                                      PALETTE_MIN_SIZE][color_ctx]
                                : frm_ctx->palette_y_color_index_cdf[palette_size -
                                                                     PALETTE_MIN_SIZE][color_ctx],
                            palette_size,
                            ACCT_STR);
                        color_map[(i - j)][j] = color_order[palette_color_idx];
                    }
                }
                for (int i = 0; i < on_screen_height; i++) {
                    for (int j = on_screen_width; j < block_width; j++) {
                        color_map[i][j] = color_map[i][on_screen_width - 1];
                    }
                }
                for (int i = on_screen_height; i < block_height; i++) {
                    for (int j = 0; j < block_width; j++) {
                        color_map[i][j] = color_map[on_screen_height - 1][j];
                    }
                }
            }
        }

        if ((plane_itr ? is_chroma_ref : 1)) {
            if (palette_size) {
                /* Palette prediction process */
                void *               blk_recon_buf;
                int32_t              recon_stride;
                EbPictureBufferDesc *recon_picture_buf = dec_handle->cur_pic_buf[0]->ps_pic_buf;

                derive_blk_pointers(recon_picture_buf,
                                    plane_itr,
                                    (mi_col >> sub_x) * MI_SIZE,
                                    (mi_row >> sub_y) * MI_SIZE,
                                    &blk_recon_buf,
                                    &recon_stride,
                                    sub_x,
                                    sub_y);
                uint16_t *palette = parse_ctx->palette_colors[plane_itr];
                if (recon_picture_buf->bit_depth == EB_8BIT &&
                    !(dec_handle->is_16bit_pipeline))
                {
                    uint8_t *temp_buf = (uint8_t *)blk_recon_buf;
                    for (int i = 0; i < block_height; i++) {
                        for (int j = 0; j < block_width; j++) {
                            temp_buf[i * recon_stride + j] = (uint8_t)palette[color_map[i][j]];
                        }
                    }
                } else {
                    uint16_t *temp_buf = (uint16_t *)blk_recon_buf;
                    for (int i = 0; i < block_height; i++) {
                        for (int j = 0; j < block_width; j++) {
                            temp_buf[i * recon_stride + j] = palette[color_map[i][j]];
                        }
                    }
                }
            }
        }
    }
}
