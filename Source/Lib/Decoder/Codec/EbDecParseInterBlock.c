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

static uint16_t compound_mode_ctx_map[3][COMP_NEWMV_CTXS] = {
    { 0, 1, 1, 1, 1 },
    { 1, 2, 3, 4, 4 },
    { 4, 4, 5, 6, 7 },
};

static INLINE void svt_collect_neighbors_ref_counts(PartitionInfo_t *pi) {

    ZERO_ARRAY(&pi->neighbors_ref_counts[0],
        sizeof(pi->neighbors_ref_counts[0])*REF_FRAMES);

    uint8_t *const ref_counts = pi->neighbors_ref_counts;

    const ModeInfo_t *const above_mbmi = pi->above_mbmi;
    const ModeInfo_t *const left_mbmi = pi->left_mbmi;
    const int above_in_image = pi->up_available;
    const int left_in_image = pi->left_available;

    // Above neighbor
    if (above_in_image && dec_is_inter_block(above_mbmi)) {
        ref_counts[above_mbmi->ref_frame[0]]++;
        if (has_second_ref(above_mbmi))
            ref_counts[above_mbmi->ref_frame[1]]++;
    }

    // Left neighbor
    if (left_in_image && dec_is_inter_block(left_mbmi)) {
        ref_counts[left_mbmi->ref_frame[0]]++;
        if (has_second_ref(left_mbmi))
            ref_counts[left_mbmi->ref_frame[1]]++;
    }
}

#define CHECK_BACKWARD_REFS(ref_frame) \
  (((ref_frame) >= BWDREF_FRAME) && ((ref_frame) <= ALTREF_FRAME))
#define IS_BACKWARD_REF_FRAME(ref_frame) CHECK_BACKWARD_REFS(ref_frame)

static INLINE int is_inside(TileInfo *tile, int mi_col, int mi_row) {
    return (mi_col >= tile->mi_col_start && mi_col < tile->mi_col_end &&
        mi_row >= tile->mi_row_start && mi_row < tile->mi_row_end);
}

static int get_reference_mode_context(const PartitionInfo_t *xd) {
    int ctx;
    const ModeInfo_t *const above_mbmi = xd->above_mbmi;
    const ModeInfo_t *const left_mbmi = xd->left_mbmi;
    const int has_above = xd->up_available;
    const int has_left = xd->left_available;

    // Note:
    // The mode info data structure has a one element border above and to the
    // left of the entries corresponding to real macroblocks.
    // The prediction flags in these dummy entries are initialized to 0.
    if (has_above && has_left) {  // both edges available
        if (!has_second_ref(above_mbmi) && !has_second_ref(left_mbmi))
            // neither edge uses comp pred (0/1)
            ctx = IS_BACKWARD_REF_FRAME(above_mbmi->ref_frame[0]) ^
            IS_BACKWARD_REF_FRAME(left_mbmi->ref_frame[0]);
        else if (!has_second_ref(above_mbmi))
            // one of two edges uses comp pred (2/3)
            ctx = 2 + (IS_BACKWARD_REF_FRAME(above_mbmi->ref_frame[0]) ||
                !dec_is_inter_block(above_mbmi));
        else if (!has_second_ref(left_mbmi))
            // one of two edges uses comp pred (2/3)
            ctx = 2 + (IS_BACKWARD_REF_FRAME(left_mbmi->ref_frame[0]) ||
                !dec_is_inter_block(left_mbmi));
        else  // both edges use comp pred (4)
            ctx = 4;
    }
    else if (has_above || has_left) {  // one edge available
        const ModeInfo_t *edge_mbmi = has_above ? above_mbmi : left_mbmi;

        if (!has_second_ref(edge_mbmi))
            // edge does not use comp pred (0/1)
            ctx = IS_BACKWARD_REF_FRAME(edge_mbmi->ref_frame[0]);
        else
            // edge uses comp pred (3)
            ctx = 3;
    }
    else {  // no edges available (1)
        ctx = 1;
    }
    assert(ctx >= 0 && ctx < COMP_INTER_CONTEXTS);
    return ctx;
}

static int32_t get_pred_context_comp_ref_p(PartitionInfo_t *pi) {
    const uint8_t *const ref_counts = &pi->neighbors_ref_counts[0];

    // Count of LAST + LAST2
    const int32_t last_last2_count = ref_counts[LAST_FRAME] + ref_counts[LAST2_FRAME];
    // Count of LAST3 + GOLDEN
    const int32_t last3_gld_count =
        ref_counts[LAST3_FRAME] + ref_counts[GOLDEN_FRAME];

    const int32_t pred_context = (last_last2_count == last3_gld_count) ? 1 :
        ((last_last2_count < last3_gld_count) ? 0 : 2);

    assert(pred_context >= 0 && pred_context < REF_CONTEXTS);
    return pred_context;
}

static int32_t get_pred_context_comp_bwdref_p(PartitionInfo_t *pi) {
    const uint8_t *const ref_counts = &pi->neighbors_ref_counts[0];

    // Counts of BWDREF, ALTREF2, or ALTREF frames (B, A2, or A)
    const int32_t brfarf2_count =
        ref_counts[BWDREF_FRAME] + ref_counts[ALTREF2_FRAME];
    const int32_t arf_count = ref_counts[ALTREF_FRAME];

    const int32_t pred_context =
        (brfarf2_count == arf_count) ? 1 : ((brfarf2_count < arf_count) ? 0 : 2);

    assert(pred_context >= 0 && pred_context < REF_CONTEXTS);
    return pred_context;
}

static int32_t get_pred_context_comp_bwdref_p1(PartitionInfo_t *pi) {
    const uint8_t *const ref_counts = &pi->neighbors_ref_counts[0];

    // Count of BWDREF frames (B)
    const int32_t brf_count = ref_counts[BWDREF_FRAME];
    // Count of ALTREF2 frames (A2)
    const int32_t arf2_count = ref_counts[ALTREF2_FRAME];

    const int32_t pred_context =
        (brf_count == arf2_count) ? 1 : ((brf_count < arf2_count) ? 0 : 2);

    assert(pred_context >= 0 && pred_context < REF_CONTEXTS);
    return pred_context;
}

static int get_pred_context_uni_comp_ref_p2(PartitionInfo_t *pi) {
    const uint8_t *const ref_counts = &pi->neighbors_ref_counts[0];

    // Count of LAST3
    const int last3_count = ref_counts[LAST3_FRAME];
    // Count of GOLDEN
    const int gld_count = ref_counts[GOLDEN_FRAME];

    const int pred_context =
        (last3_count == gld_count) ? 1 : ((last3_count < gld_count) ? 0 : 2);

    assert(pred_context >= 0 && pred_context < UNI_COMP_REF_CONTEXTS);
    return pred_context;
}

static int get_pred_context_uni_comp_ref_p1(PartitionInfo_t *pi) {
    const uint8_t *const ref_counts = &pi->neighbors_ref_counts[0];

    // Count of LAST2
    const int last2_count = ref_counts[LAST2_FRAME];
    // Count of LAST3 or GOLDEN
    const int last3_or_gld_count =
        ref_counts[LAST3_FRAME] + ref_counts[GOLDEN_FRAME];

    const int pred_context = (last2_count == last3_or_gld_count)
        ? 1
        : ((last2_count < last3_or_gld_count) ? 0 : 2);

    assert(pred_context >= 0 && pred_context < UNI_COMP_REF_CONTEXTS);
    return pred_context;
}

static int get_pred_context_uni_comp_ref_p(PartitionInfo_t *pi) {
    const uint8_t *const ref_counts = &pi->neighbors_ref_counts[0];

    // Count of forward references (L, L2, L3, or G)
    const int frf_count = ref_counts[LAST_FRAME] + ref_counts[LAST2_FRAME] +
        ref_counts[LAST3_FRAME] + ref_counts[GOLDEN_FRAME];
    // Count of backward references (B or A)
    const int brf_count = ref_counts[BWDREF_FRAME] + ref_counts[ALTREF2_FRAME] +
        ref_counts[ALTREF_FRAME];

    const int pred_context =
        (frf_count == brf_count) ? 1 : ((frf_count < brf_count) ? 0 : 2);

    assert(pred_context >= 0 && pred_context < UNI_COMP_REF_CONTEXTS);
    return pred_context;
}

static int32_t get_pred_context_single_ref_p1(PartitionInfo_t *pi) {
    const uint8_t *const ref_counts = &pi->neighbors_ref_counts[0];

    // Count of forward reference frames
    const int32_t fwd_count = ref_counts[LAST_FRAME] + ref_counts[LAST2_FRAME] +
        ref_counts[LAST3_FRAME] + ref_counts[GOLDEN_FRAME];
    // Count of backward reference frames
    const int32_t bwd_count = ref_counts[BWDREF_FRAME] + ref_counts[ALTREF2_FRAME] +
        ref_counts[ALTREF_FRAME];

    const int32_t pred_context =
        (fwd_count == bwd_count) ? 1 : ((fwd_count < bwd_count) ? 0 : 2);

    assert(pred_context >= 0 && pred_context < REF_CONTEXTS);
    return pred_context;
}

static int32_t get_pred_context_single_ref_p4(PartitionInfo_t *pi) {
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

static int32_t get_pred_context_last3_or_gld(PartitionInfo_t *pi) {
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

static void read_ref_frames(EbDecHandle *dec_handle, PartitionInfo_t *const pi,
    SvtReader *r)
{
    int segment_id = pi->mi->segment_id;
    MvReferenceFrame *ref_frame = pi->mi->ref_frame;
    ParseCtxt *parse_ctxt = (ParseCtxt *)dec_handle->pv_parse_ctxt;
    AomCdfProb *cdf;
    SegmentationParams *seg_params = &dec_handle->frame_header.segmentation_params;
    if (pi->mi->skip_mode) {
        ref_frame[0] = (MvReferenceFrame)(LAST_FRAME +
            dec_handle->frame_header.skip_mode_params.ref_frame_idx_0);
        ref_frame[1] = (MvReferenceFrame)(LAST_FRAME +
            dec_handle->frame_header.skip_mode_params.ref_frame_idx_1);
    }
    else if (seg_feature_active(seg_params, segment_id, SEG_LVL_REF_FRAME)) {
        ref_frame[0] = (MvReferenceFrame)get_segdata(seg_params, segment_id,
            SEG_LVL_REF_FRAME);
        ref_frame[1] = NONE_FRAME;
    }
    else if (seg_feature_active(seg_params, segment_id, SEG_LVL_SKIP) ||
        seg_feature_active(seg_params, segment_id, SEG_LVL_GLOBALMV))
    {
        ref_frame[0] = LAST_FRAME;
        ref_frame[1] = NONE_FRAME;
    }
    else {
        ReferenceMode mode = SINGLE_REFERENCE;
        int bw4 = mi_size_wide[pi->mi->sb_type];
        int bh4 = mi_size_high[pi->mi->sb_type];
        if (dec_handle->frame_header.reference_mode == REFERENCE_MODE_SELECT
            && (AOMMIN(bw4, bh4) >= 2))
        {
            const int ctx = get_reference_mode_context(pi);
            mode = (ReferenceMode)svt_read_symbol(r,
                parse_ctxt->cur_tile_ctx.comp_inter_cdf[ctx], 2, ACCT_STR);
        }

        if (mode == COMPOUND_REFERENCE) {
            int pred_context;
            const int ctx = get_comp_reference_type_context(pi);
            const CompReferenceType comp_ref_type =
                (CompReferenceType)svt_read_symbol(
                    r, parse_ctxt->cur_tile_ctx.comp_ref_type_cdf[ctx], 2, ACCT_STR);

            if (comp_ref_type == UNIDIR_COMP_REFERENCE) {
                pred_context = get_pred_context_uni_comp_ref_p(pi);
                uint16_t bit = (uint16_t)svt_read_symbol(r,
                    parse_ctxt->cur_tile_ctx.uni_comp_ref_cdf[pred_context][0],
                    2, ACCT_STR);
                if (bit) {
                    ref_frame[0] = BWDREF_FRAME;
                    ref_frame[1] = ALTREF_FRAME;
                }
                else {
                    pred_context = get_pred_context_uni_comp_ref_p1(pi);
                    uint16_t bit1 = (uint16_t)svt_read_symbol(r,
                        parse_ctxt->cur_tile_ctx.uni_comp_ref_cdf[pred_context][1],
                        2, ACCT_STR);
                    if (bit1) {
                        pred_context = get_pred_context_uni_comp_ref_p2(pi);
                        uint16_t bit2 = (uint16_t)svt_read_symbol(r,
                            parse_ctxt->cur_tile_ctx.uni_comp_ref_cdf[pred_context][2],
                            2, ACCT_STR);
                        if (bit2) {
                            ref_frame[0] = LAST_FRAME;
                            ref_frame[1] = GOLDEN_FRAME;
                        }
                        else {
                            ref_frame[0] = LAST_FRAME;
                            ref_frame[1] = LAST3_FRAME;
                        }
                    }
                    else {
                        ref_frame[0] = LAST_FRAME;
                        ref_frame[1] = LAST2_FRAME;
                    }
                }
                return;
            }

            assert(comp_ref_type == BIDIR_COMP_REFERENCE);

            const int idx = 1;
            pred_context = get_pred_context_comp_ref_p(pi);
            uint16_t bit = (uint16_t)svt_read_symbol(r,
                parse_ctxt->cur_tile_ctx.comp_ref_cdf[pred_context][0], 2, ACCT_STR);
            // Decode forward references.
            if (!bit) {
                uint16_t bit1 = (uint16_t)svt_read_symbol(r, parse_ctxt->cur_tile_ctx.
                    comp_ref_cdf[get_pred_context_single_ref_p4(pi)][1], 2, ACCT_STR);
                ref_frame[!idx] = bit1 ? LAST2_FRAME : LAST_FRAME;
            }
            else {
                uint16_t bit2 = (uint16_t)svt_read_symbol(r, parse_ctxt->cur_tile_ctx.
                    comp_ref_cdf[get_pred_context_last3_or_gld(pi)][2], 2, ACCT_STR);
                ref_frame[!idx] = bit2 ? GOLDEN_FRAME : LAST3_FRAME;
            }

            // Decode backward references.
            pred_context = get_pred_context_comp_bwdref_p(pi);
            uint16_t bit_bwd = (uint16_t)svt_read_symbol(r, parse_ctxt->cur_tile_ctx.
                comp_bwdref_cdf[pred_context][0], 2, ACCT_STR);
            if (!bit_bwd) {
                pred_context = get_pred_context_comp_bwdref_p1(pi);
                uint16_t bit1_bwd = (uint16_t)svt_read_symbol(r, parse_ctxt->cur_tile_ctx.
                    comp_bwdref_cdf[pred_context][1], 2, ACCT_STR);
                ref_frame[idx] = bit1_bwd ? ALTREF2_FRAME : BWDREF_FRAME;
            }
            else {
                ref_frame[idx] = ALTREF_FRAME;
            }
        }
        else if (mode == SINGLE_REFERENCE) {

            cdf = parse_ctxt->cur_tile_ctx.
                single_ref_cdf[get_pred_context_single_ref_p1(pi)][0];
            const int32_t bit0 = svt_read_symbol(r, cdf, 2, ACCT_STR);

            if (bit0) {
                cdf = parse_ctxt->cur_tile_ctx.
                    single_ref_cdf[get_pred_context_comp_bwdref_p(pi)][1];
                const int32_t bit1 = svt_read_symbol(r, cdf, 2, ACCT_STR);
                if (!bit1) {
                    cdf = parse_ctxt->cur_tile_ctx.
                        single_ref_cdf[get_pred_context_comp_bwdref_p1(pi)][5];
                    const int32_t bit5 = svt_read_symbol(r, cdf, 2, ACCT_STR);
                    ref_frame[0] = bit5 ? ALTREF2_FRAME : BWDREF_FRAME;
                }
                else {
                    ref_frame[0] = ALTREF_FRAME;
                }
            }
            else {
                cdf = parse_ctxt->cur_tile_ctx.
                    single_ref_cdf[get_pred_context_comp_ref_p(pi)][2];
                const int32_t bit2 = svt_read_symbol(r, cdf, 2, ACCT_STR);
                if (bit2) {
                    cdf = parse_ctxt->cur_tile_ctx.
                        single_ref_cdf[get_pred_context_last3_or_gld(pi)][4];
                    const int32_t bit4 = svt_read_symbol(r, cdf, 2, ACCT_STR);
                    ref_frame[0] = bit4 ? GOLDEN_FRAME : LAST3_FRAME;
                }
                else {
                    cdf = parse_ctxt->cur_tile_ctx.
                        single_ref_cdf[get_pred_context_single_ref_p4(pi)][3];
                    const int32_t bit3 = svt_read_symbol(r, cdf, 2, ACCT_STR);
                    ref_frame[0] = bit3 ? LAST2_FRAME : LAST_FRAME;
                }
            }

            ref_frame[1] = NONE_FRAME;
        }
        else
            assert(0 && "Invalid prediction mode.");
    }
}

static INLINE int is_global_mv_block(const ModeInfo_t *const mbmi,
    TransformationType type)
{
    const PredictionMode mode = mbmi->mode;
    const BlockSize bsize = mbmi->sb_type;
    const int block_size_allowed =
        AOMMIN(block_size_wide[bsize], block_size_high[bsize]) >= 8;
    return (mode == GLOBALMV || mode == GLOBAL_GLOBALMV) && type > TRANSLATION &&
        block_size_allowed;
}

int has_newmv(PredictionMode mode) {
    return (mode == NEWMV ||
        mode == NEW_NEWMV ||
        mode == NEAR_NEWMV ||
        mode == NEW_NEARMV ||
        mode == NEAREST_NEWMV ||
        mode == NEW_NEARESTMV);
}

static void add_ref_mv_candidate(EbDecHandle *dec_handle,
    const ModeInfo_t *const candidate, const MvReferenceFrame rf[2],
    uint8_t *num_mv_found, uint8_t *found_match, uint8_t *newmv_count,
    CandidateMv_dec *ref_mv_stack, IntMv_dec *gm_mv_candidates, int weight)
{
    if (!dec_is_inter_block(candidate)) return;  // for intrabc
    int index = 0, ref;
    assert(weight % 2 == 0);

    EbDecPicBuf *buf = dec_handle->cur_pic_buf[0];
    if (rf[1] == NONE_FRAME) {
        // single reference frame
        for (ref = 0; ref < 2; ++ref) {
            if (candidate->ref_frame[ref] == rf[0]) {
                IntMv_dec this_refmv;
                if (is_global_mv_block(candidate, buf->global_motion[rf[0]].gm_type))
                    this_refmv = gm_mv_candidates[0];
                else
                    this_refmv = candidate->mv[ref];

                for (index = 0; index < *num_mv_found; ++index)
                    if (ref_mv_stack[index].this_mv.as_int == this_refmv.as_int) break;

                if (index < *num_mv_found) ref_mv_stack[index].weight += weight;

                // Add a new item to the list.
                if (index == *num_mv_found && *num_mv_found < MAX_REF_MV_STACK_SIZE) {
                    ref_mv_stack[index].this_mv.as_int = this_refmv.as_int;
                    ref_mv_stack[index].weight = weight;
                    ++(*num_mv_found);
                }
                if (has_newmv(candidate->mode))++*newmv_count;
                ++*found_match;
            }
        }
    }
    else {
        // compound reference frame
        if (candidate->ref_frame[0] == rf[0] && candidate->ref_frame[1] == rf[1]) {
            IntMv_dec this_refmv[2];
            for (ref = 0; ref < 2; ++ref) {
                if (is_global_mv_block(candidate, buf->global_motion[rf[ref]].gm_type))
                    this_refmv[ref] = gm_mv_candidates[ref];
                else
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
                ref_mv_stack[index].weight = weight;
                ++(*num_mv_found);
            }
            if (has_newmv(candidate->mode))++*newmv_count;
            ++*found_match;
        }
    }
}

static void scan_row_mbmi(EbDecHandle *dec_handle, PartitionInfo_t *pi,
    int delta_row, int mi_row, int mi_col, const MvReferenceFrame rf[2],
    CandidateMv_dec *ref_mv_stack, uint8_t *num_mv_found, uint8_t *found_match,
    uint8_t *newmv_count, IntMv_dec *gm_mv_candidates, int max_row_offset,
    int *processed_rows)
{
    int bw4 = mi_size_wide[pi->mi->sb_type];
    ParseCtxt *parse_ctx = (ParseCtxt*)dec_handle->pv_parse_ctxt;
    FrameHeader *frm_header = &dec_handle->frame_header;
    int end4 = AOMMIN(AOMMIN(bw4, (int)frm_header->mi_cols - mi_col), 16);
    int delta_col = 0;
    int use_step_16 = (bw4 >= 16);
    const int n8_w_8 = mi_size_wide[BLOCK_8X8];
    const int n8_w_16 = mi_size_wide[BLOCK_16X16];

    if (abs(delta_row) > 1) {
        delta_col = 1;
        if ((mi_col & 0x01) && bw4 < n8_w_8) --delta_col;
    }

    for (int i = 0; i < end4;) {
        int mv_row = mi_row + delta_row;
        int mv_col = mi_col + delta_col + i;
        if (!is_inside(&parse_ctx->cur_tile_info, mv_col, mv_row))
            break;
        ModeInfo_t *candidate = get_cur_mode_info(dec_handle,
            mv_row, mv_col, pi->sb_info);
        int len = AOMMIN(bw4, mi_size_wide[candidate->sb_type]);
        const int n4_w = mi_size_wide[candidate->sb_type];
        if (use_step_16)
            len = AOMMAX(n8_w_16, len);
        else if (abs(delta_row) > 1)
            len = AOMMAX(n8_w_8, len);

        int weight = 2;
        if (bw4 >= n8_w_8 && bw4 <= n4_w) {
            int inc = AOMMIN(-max_row_offset + delta_row + 1,
                mi_size_high[candidate->sb_type]);
            // Obtain range used in weight calculation.
            weight = AOMMAX(weight, inc);
            *processed_rows = inc - delta_row - 1;
        }
        add_ref_mv_candidate(dec_handle, candidate, rf, num_mv_found, found_match,
            newmv_count, ref_mv_stack, gm_mv_candidates, len * weight);

        i += len;
    }
}

static void scan_col_mbmi(EbDecHandle *dec_handle, PartitionInfo_t *pi,
    int delta_col, int mi_row, int mi_col, const MvReferenceFrame rf[2],
    CandidateMv_dec *ref_mv_stack, uint8_t *num_mv_found, uint8_t *found_match,
    uint8_t *newmv_count, IntMv_dec *gm_mv_candidates, int max_col_offset,
    int *processed_cols)
{
    int bh4 = mi_size_high[pi->mi->sb_type];
    ParseCtxt *parse_ctx = (ParseCtxt*)dec_handle->pv_parse_ctxt;
    FrameHeader *frm_header = &dec_handle->frame_header;
    int end4 = AOMMIN(AOMMIN(bh4, (int)frm_header->mi_rows - mi_row), 16);
    int delta_row = 0;
    int use_step_16 = (bh4 >= 16);
    const int n8_h_8 = mi_size_high[BLOCK_8X8];

    if (abs(delta_col) > 1) {
        delta_row = 1;
        if ((mi_row & 0x01) && bh4 < n8_h_8) --delta_row;
    }

    for (int i = 0; i < end4;) {
        int mv_row = mi_row + delta_row + i;
        int mv_col = mi_col + delta_col;
        if (!is_inside(&parse_ctx->cur_tile_info, mv_col, mv_row))
            break;
        ModeInfo_t *candidate = get_cur_mode_info(dec_handle,
            mv_row, mv_col, pi->sb_info);
        int len = AOMMIN(bh4, mi_size_high[candidate->sb_type]);
        const int n4_h = mi_size_high[candidate->sb_type];
        if (abs(delta_col) > 1)
            len = AOMMAX(2, len);
        if (use_step_16)
            len = AOMMAX(4, len);

        int weight = 2;
        if (bh4 >= n8_h_8 && bh4 <= n4_h) {
            int inc = AOMMIN(-max_col_offset + delta_col + 1,
                mi_size_wide[candidate->sb_type]);
            // Obtain range used in weight calculation.
            weight = AOMMAX(weight, inc);
            *processed_cols = inc - delta_col - 1;
        }

        add_ref_mv_candidate(dec_handle, candidate, rf, num_mv_found, found_match,
            newmv_count, ref_mv_stack, gm_mv_candidates, len * weight);

        i += len;
    }
}

static void scan_blk_mbmi(EbDecHandle *dec_handle, PartitionInfo_t *pi,
    int delta_row, int delta_col, const int mi_row, const int mi_col,
    const MvReferenceFrame rf[2], CandidateMv_dec *ref_mv_stack,
    uint8_t *found_match, uint8_t *newmv_count, IntMv_dec *gm_mv_candidates,
    uint8_t num_mv_found[MODE_CTX_REF_FRAMES])
{
    ParseCtxt *parse_ctx = (ParseCtxt*)dec_handle->pv_parse_ctxt;

    int mv_row = mi_row + delta_row;
    int mv_col = mi_col + delta_col;
    int weight = 4;

    if (is_inside(&parse_ctx->cur_tile_info, mv_col, mv_row)) {
        ModeInfo_t * candidate = get_cur_mode_info(dec_handle,
            mv_row, mv_col, pi->sb_info);

        add_ref_mv_candidate(dec_handle, candidate, rf, num_mv_found, found_match,
            newmv_count, ref_mv_stack, gm_mv_candidates, weight);
    }  // Analyze a single 8x8 block motion information.
}

/* TODO: Harmonize with Encoder. */
static int has_top_right(EbDecHandle *dec_handle, PartitionInfo_t *pi,
    int mi_row, int mi_col, int bs)
{
    int n4_w = mi_size_wide[pi->mi->sb_type];
    int n4_h = mi_size_high[pi->mi->sb_type];
    const int sb_mi_size = mi_size_wide[dec_handle->seq_header.sb_size];
    const int mask_row = mi_row & (sb_mi_size - 1);
    const int mask_col = mi_col & (sb_mi_size - 1);

    if (bs > mi_size_wide[BLOCK_64X64]) return 0;
    int has_tr = !((mask_row & bs) && (mask_col & bs));

    assert(bs > 0 && !(bs & (bs - 1)));

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

static INLINE void integer_mv_precision(MV *mv) {
    int mod = (mv->row % 8);
    if (mod != 0) {
        mv->row -= mod;
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
        mv->col -= mod;
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

static int add_tpl_ref_mv(EbDecHandle *dec_handle, int mi_row, int mi_col,
    MvReferenceFrame ref_frame, int blk_row, int blk_col,
    IntMv_dec *gm_mv_candidates, uint8_t *num_mv_found,
    CandidateMv_dec ref_mv_stacks[][MAX_REF_MV_STACK_SIZE], int16_t *mode_context)
{
    uint8_t idx;
    ParseCtxt *parse_ctx = (ParseCtxt*)dec_handle->pv_parse_ctxt;
    FrameHeader *frm_header = &dec_handle->frame_header;
    int mv_row = (mi_row + blk_row) | 1;
    int mv_col = (mi_col + blk_col) | 1;

    if (!is_inside(&parse_ctx->cur_tile_info, mv_col, mv_row)) return 0;

    int x8 = mv_col >> 1;
    int y8 = mv_row >> 1;

    MvReferenceFrame rf[2];
    av1_set_ref_frame(rf, ref_frame);

    const TemporalMvRef *tpl_mvs = dec_handle->master_frame_buf.tpl_mvs +
        y8 * (frm_header->mi_stride >> 1) + x8;
    const IntMv_dec prev_frame_mvs = tpl_mvs->mf_mv0;

    if (rf[1] == NONE_FRAME) {
        CandidateMv_dec *ref_mv_stack = ref_mv_stacks[rf[0]];

        if (prev_frame_mvs.as_int == INVALID_MV)
            return 0;

        IntMv_dec this_refmv;

        lower_mv_precision(&this_refmv.as_mv, frm_header->allow_high_precision_mv,
            frm_header->force_integer_mv);

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
            ref_mv_stack[idx].weight = 2;
            ++(*num_mv_found);
        }
        return 1;
    }
    else {
        // Process compound inter mode
        CandidateMv_dec *ref_mv_stack = ref_mv_stacks[ref_frame];
        if (prev_frame_mvs.as_int == INVALID_MV)
            return 0;

        IntMv_dec this_refmv;
        IntMv_dec comp_refmv;

        lower_mv_precision(&this_refmv.as_mv, frm_header->allow_high_precision_mv,
            frm_header->force_integer_mv);
        lower_mv_precision(&comp_refmv.as_mv, frm_header->allow_high_precision_mv,
            frm_header->force_integer_mv);

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
            ref_mv_stack[idx].weight = 2;
            ++(*num_mv_found);
        }
    }
    return 1;
}

static int check_sb_border(const int mi_row, const int mi_col,
    const int row_offset, const int col_offset) {
    const int sb_mi_size = mi_size_wide[BLOCK_64X64];
    const int row = (mi_row & (sb_mi_size - 1)) + row_offset;
    const int col = (mi_col & (sb_mi_size - 1)) + col_offset;

    return ( (row >= 0) && (row < sb_mi_size) &&
             (col >= 0) && (col < sb_mi_size) );
}


static void add_extra_mv_candidate(ModeInfo_t * candidate, EbDecHandle *dec_handle,
    MvReferenceFrame *rf, IntMv_dec ref_id[2][2], int ref_id_count[2],
    IntMv_dec ref_diff[2][2], int ref_diff_count[2])
{
    FrameHeader *frm_header = &dec_handle->frame_header;
    for (int rf_idx = 0; rf_idx < 2; ++rf_idx) {
        MvReferenceFrame can_rf = candidate->ref_frame[rf_idx];
        if (can_rf > INTRA_FRAME) {
            for (int cmp_idx = 0; cmp_idx < 2; ++cmp_idx) {
                if (can_rf == rf[cmp_idx] && ref_id_count[cmp_idx] < 2) {
                    ref_id[cmp_idx][ref_id_count[cmp_idx]] = candidate->mv[rf_idx];
                    ++ref_id_count[cmp_idx];
                }
                else if (ref_diff_count[cmp_idx] < 2) {
                    IntMv_dec this_mv = candidate->mv[rf_idx];
                    if (frm_header->ref_frame_sign_bias[can_rf] !=
                        frm_header->ref_frame_sign_bias[rf[cmp_idx]])
                    {
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

static void process_single_ref_mv_candidate(ModeInfo_t * candidate,
    EbDecHandle *dec_handle, MvReferenceFrame ref_frame,
    uint8_t refmv_count[MODE_CTX_REF_FRAMES],
    CandidateMv_dec ref_mv_stack[][MAX_REF_MV_STACK_SIZE])
{
    FrameHeader *frm_header = &dec_handle->frame_header;
    for (int rf_idx = 0; rf_idx < 2; ++rf_idx) {
        if (candidate->ref_frame[rf_idx] > INTRA_FRAME) {
            IntMv_dec this_mv = candidate->mv[rf_idx];
            if (frm_header->ref_frame_sign_bias[candidate->ref_frame[rf_idx]] !=
                frm_header->ref_frame_sign_bias[ref_frame])
            {
                this_mv.as_mv.row = -this_mv.as_mv.row;
                this_mv.as_mv.col = -this_mv.as_mv.col;
            }
            int stack_idx;
            for (stack_idx = 0; stack_idx < refmv_count[ref_frame]; ++stack_idx) {
                const IntMv_dec stack_mv = ref_mv_stack[ref_frame][stack_idx].this_mv;
                if (this_mv.as_int == stack_mv.as_int) break;
            }

            if (stack_idx == refmv_count[ref_frame]) {
                ref_mv_stack[ref_frame][stack_idx].this_mv = this_mv;
                ref_mv_stack[ref_frame][stack_idx].weight = 2;
                ++refmv_count[ref_frame];
            }
        }
    }
}

static INLINE void clamp_mv(MV *mv, int min_col, int max_col, int min_row,
    int max_row) {
    mv->col = clamp(mv->col, min_col, max_col);
    mv->row = clamp(mv->row, min_row, max_row);
}

static INLINE void clamp_mv_ref(MV *mv, int bw, int bh, PartitionInfo_t *pi) {
    clamp_mv(mv, pi->mb_to_left_edge - bw * 8 - MV_BORDER,
        pi->mb_to_right_edge + bw * 8 + MV_BORDER,
        pi->mb_to_top_edge - bh * 8 - MV_BORDER,
        pi->mb_to_bottom_edge + bh * 8 + MV_BORDER);
}

static void dec_setup_ref_mv_list(
    EbDecHandle *dec_handle, PartitionInfo_t *pi, MvReferenceFrame ref_frame,
    CandidateMv_dec ref_mv_stack[][MAX_REF_MV_STACK_SIZE],
    IntMv_dec mv_ref_list[][MAX_MV_REF_CANDIDATES], IntMv_dec *gm_mv_candidates,
    int mi_row, int mi_col, int16_t *mode_context, MvCount *mv_cnt)
{
    int n4_w = mi_size_wide[pi->mi->sb_type];
    int n4_h = mi_size_high[pi->mi->sb_type];
    const int bs = AOMMAX(n4_w, n4_h);
    MvReferenceFrame rf[2];

    ParseCtxt *parse_ctx = (ParseCtxt*)dec_handle->pv_parse_ctxt;
    FrameHeader *frame_info = &dec_handle->frame_header;
    const TileInfo *const tile = &parse_ctx->cur_tile_info;
    int max_row_offset = 0, max_col_offset = 0;
    const int row_adj = (n4_h < mi_size_high[BLOCK_8X8]) && (mi_row & 0x01);
    const int col_adj = (n4_w < mi_size_wide[BLOCK_8X8]) && (mi_col & 0x01);
    int processed_rows = 0;
    int processed_cols = 0;

    av1_set_ref_frame(rf, ref_frame);
    mode_context[ref_frame] = 0;

    // Find valid maximum row/col offset.
    if (pi->up_available) {
        max_row_offset = -(MVREF_ROW_COLS << 1) + row_adj;

        if (n4_h < mi_size_high[BLOCK_8X8])
            max_row_offset = -(2 << 1) + row_adj;

        max_row_offset = clamp(max_row_offset, tile->mi_row_start - mi_row,
            tile->mi_row_end - mi_row - 1);
    }

    if (pi->left_available) {
        max_col_offset = -(MVREF_ROW_COLS << 1) + col_adj;

        if (n4_w < mi_size_wide[BLOCK_8X8])
            max_col_offset = -(2 << 1) + col_adj;

        max_col_offset = clamp(max_col_offset, tile->mi_col_start - mi_col,
            tile->mi_col_end - mi_col - 1);
    }
    memset(mv_cnt, 0, sizeof(*mv_cnt));

    // Scan the first above row mode info. row_offset = -1;
    if (abs(max_row_offset) >= 1) {
        scan_row_mbmi(dec_handle, pi, -1, mi_row, mi_col, rf, ref_mv_stack[ref_frame],
            &mv_cnt->num_mv_found[ref_frame], &mv_cnt->found_above_match,
            &mv_cnt->newmv_count, gm_mv_candidates, max_row_offset, &processed_rows);
    }

    // Scan the first left column mode info. col_offset = -1;
    if (abs(max_col_offset) >= 1) {
        scan_col_mbmi(dec_handle, pi, -1, mi_row, mi_col, rf, ref_mv_stack[ref_frame],
            &mv_cnt->num_mv_found[ref_frame], &mv_cnt->found_left_match,
            &mv_cnt->newmv_count, gm_mv_candidates, max_col_offset, &processed_cols);
    }

    if (has_top_right(dec_handle, pi, mi_row, mi_col, bs)) {
        scan_blk_mbmi(dec_handle, pi, -1, n4_w, mi_row, mi_col, rf,
            ref_mv_stack[ref_frame], &mv_cnt->found_above_match, &mv_cnt->newmv_count,
            gm_mv_candidates, &mv_cnt->num_mv_found[ref_frame]);
    }

    const uint8_t nearest_match = (mv_cnt->found_above_match > 0) +
        (mv_cnt->found_left_match > 0);
    const uint8_t num_nearest = mv_cnt->num_mv_found[ref_frame];
    const uint8_t num_new = mv_cnt->newmv_count;

    for (int idx = 0; idx < num_nearest; ++idx)
        ref_mv_stack[ref_frame][idx].weight += REF_CAT_LEVEL;

    if (frame_info->use_ref_frame_mvs) {
        int is_available = 0;
        const int voffset = AOMMAX(mi_size_high[BLOCK_8X8], n4_h);
        const int hoffset = AOMMAX(mi_size_wide[BLOCK_8X8], n4_w);
        const int blk_row_end = AOMMIN(n4_h, mi_size_high[BLOCK_64X64]);
        const int blk_col_end = AOMMIN(n4_w, mi_size_wide[BLOCK_64X64]);

        const int tpl_sample_pos[3][2] = {
          { voffset, -2 },
          { voffset, hoffset },
          { voffset - 2, hoffset },
        };
        const int allow_extension = (n4_h >= mi_size_high[BLOCK_8X8]) &&
            (n4_h < mi_size_high[BLOCK_64X64]) &&
            (n4_w >= mi_size_wide[BLOCK_8X8]) &&
            (n4_w < mi_size_wide[BLOCK_64X64]);

        const int step_h = (n4_h >= mi_size_high[BLOCK_64X64])
            ? mi_size_high[BLOCK_16X16]
            : mi_size_high[BLOCK_8X8];
        const int step_w = (n4_w >= mi_size_wide[BLOCK_64X64])
            ? mi_size_wide[BLOCK_16X16]
            : mi_size_wide[BLOCK_8X8];

        for (int blk_row = 0; blk_row < blk_row_end; blk_row += step_h) {
            for (int blk_col = 0; blk_col < blk_col_end; blk_col += step_w) {

                int ret = add_tpl_ref_mv(dec_handle, mi_row, mi_col,
                    ref_frame, blk_row, blk_col, gm_mv_candidates,
                    &mv_cnt->num_mv_found[ref_frame], ref_mv_stack, mode_context);
                if (blk_row == 0 && blk_col == 0) is_available = ret;
            }
        }

        if (is_available == 0) mode_context[ref_frame] |= (1 << GLOBALMV_OFFSET);

        if (allow_extension) {
            for (int i = 0; i < 3; ++i) {
                const int blk_row = tpl_sample_pos[i][0];
                const int blk_col = tpl_sample_pos[i][1];

                if (check_sb_border(mi_row, mi_col, blk_row, blk_col)) {
                    add_tpl_ref_mv(dec_handle, mi_row, mi_col, ref_frame, blk_row,
                        blk_col, gm_mv_candidates, &mv_cnt->num_mv_found[ref_frame],
                        ref_mv_stack, mode_context);
                }
            }
        }
    }

    // Scan the second outer area.
    scan_blk_mbmi(dec_handle, pi, -1, -1, mi_row, mi_col, rf,
        ref_mv_stack[ref_frame], &mv_cnt->found_above_match, &mv_cnt->newmv_count,
        gm_mv_candidates, &mv_cnt->num_mv_found[ref_frame]);

    for (int idx = 2; idx <= MVREF_ROW_COLS; ++idx) {
        const int row_offset = -(idx << 1) + 1 + row_adj;
        const int col_offset = -(idx << 1) + 1 + col_adj;
        if (abs(row_offset) <= abs(max_row_offset) && abs(row_offset) > processed_rows) {
            scan_row_mbmi(dec_handle, pi, row_offset, mi_row, mi_col, rf,
                ref_mv_stack[ref_frame], &mv_cnt->num_mv_found[ref_frame],
                &mv_cnt->found_above_match, &mv_cnt->newmv_count,
                gm_mv_candidates, max_row_offset, &processed_rows);
        }

        if (abs(col_offset) <= abs(max_col_offset) && abs(col_offset) > processed_cols) {
            scan_col_mbmi(dec_handle, pi, col_offset, mi_row, mi_col, rf,
                ref_mv_stack[ref_frame], &mv_cnt->num_mv_found[ref_frame],
                &mv_cnt->found_left_match, &mv_cnt->newmv_count,
                gm_mv_candidates, max_col_offset, &processed_cols);
        }
    }

    /* sorting process*/
    int start = 0;
    int end = num_nearest;
    while (end > start) {
        int new_end = start;
        for (int idx = start + 1; idx < end; ++idx) {
            if (ref_mv_stack[ref_frame][idx - 1].weight <
                ref_mv_stack[ref_frame][idx].weight) {
                CandidateMv_dec tmp_mv = ref_mv_stack[ref_frame][idx - 1];
                ref_mv_stack[ref_frame][idx - 1] = ref_mv_stack[ref_frame][idx];
                ref_mv_stack[ref_frame][idx] = tmp_mv;
                new_end = idx;
            }
        }
        end = new_end;
    }

    start = num_nearest;
    end = mv_cnt->num_mv_found[ref_frame];
    while (end > start) {
        int new_end = start;
        for (int idx = start + 1; idx < end; ++idx) {
            if (ref_mv_stack[ref_frame][idx - 1].weight <
                ref_mv_stack[ref_frame][idx].weight) {
                CandidateMv_dec tmp_mv = ref_mv_stack[ref_frame][idx - 1];
                ref_mv_stack[ref_frame][idx - 1] = ref_mv_stack[ref_frame][idx];
                ref_mv_stack[ref_frame][idx] = tmp_mv;
                new_end = idx;
            }
        }
        end = new_end;
    }

    /* extra search process */
    if (mv_cnt->num_mv_found[ref_frame] < MAX_MV_REF_CANDIDATES) {
        IntMv_dec ref_id[2][2], ref_diff[2][2];
        int ref_id_count[2] = { 0 }, ref_diff_count[2] = { 0 };

        int mi_width = AOMMIN(16, n4_w);
        mi_width = AOMMIN(mi_width, (int)frame_info->mi_cols - mi_col);
        int mi_height = AOMMIN(16, n4_h);
        mi_height = AOMMIN(mi_height, (int)frame_info->mi_rows - mi_row);
        int mi_size = AOMMIN(mi_width, mi_height);

        for (int pass = 0; pass < 2; pass++) {
            int idx = 0;
            while (idx < mi_size &&
                mv_cnt->num_mv_found[ref_frame] < MAX_MV_REF_CANDIDATES)
            {
                int mv_row, mv_col;
                if (pass == 0) { mv_row = mi_row - 1; mv_col = mi_col + idx; }
                else { mv_row = mi_row + idx; mv_col = mi_col - 1; }

                if (!is_inside(&parse_ctx->cur_tile_info, mv_col, mv_row)) break;

                ModeInfo_t *nbr = get_cur_mode_info(dec_handle,
                    mv_row, mv_col, pi->sb_info);

                if (rf[1] != NONE_FRAME)
                    add_extra_mv_candidate(nbr, dec_handle, rf,
                        ref_id, ref_id_count, ref_diff, ref_diff_count);
                else
                    process_single_ref_mv_candidate(nbr, dec_handle, ref_frame,
                        mv_cnt->num_mv_found, ref_mv_stack);

                idx += pass ? mi_size_high[nbr->sb_type] : mi_size_wide[nbr->sb_type];
            }
        }

        if (rf[1] > NONE_FRAME) {
            IntMv_dec comp_list[3][2];

            for (int idx = 0; idx < 2; ++idx) {
                int comp_idx = 0;
                for (int list_idx = 0; list_idx < ref_id_count[idx];
                    ++list_idx, comp_idx++) {
                    comp_list[comp_idx][idx] = ref_id[idx][list_idx];
                }

                for (int list_idx = 0; list_idx < ref_diff_count[idx]
                    && comp_idx < 2; ++list_idx, ++comp_idx)
                {
                    comp_list[comp_idx][idx] = ref_diff[idx][list_idx];
                }

                for (; comp_idx < 2; ++comp_idx)
                    comp_list[comp_idx][idx] = gm_mv_candidates[idx];
            }

            if (mv_cnt->num_mv_found[ref_frame]) {
                assert(mv_cnt->num_mv_found[ref_frame] == 1);
                if (comp_list[0][0].as_int ==
                    ref_mv_stack[ref_frame][0].this_mv.as_int &&
                    comp_list[0][1].as_int ==
                    ref_mv_stack[ref_frame][0].comp_mv.as_int)
                {
                    ref_mv_stack[ref_frame][mv_cnt->num_mv_found[ref_frame]].this_mv =
                        comp_list[1][0];
                    ref_mv_stack[ref_frame][mv_cnt->num_mv_found[ref_frame]].comp_mv =
                        comp_list[1][1];
                }
                else {
                    ref_mv_stack[ref_frame][mv_cnt->num_mv_found[ref_frame]].this_mv =
                        comp_list[0][0];
                    ref_mv_stack[ref_frame][mv_cnt->num_mv_found[ref_frame]].comp_mv =
                        comp_list[0][1];
                }
                ref_mv_stack[ref_frame][mv_cnt->num_mv_found[ref_frame]].weight = 2;
                ++mv_cnt->num_mv_found[ref_frame];
            }
            else {
                for (int idx = 0; idx < MAX_MV_REF_CANDIDATES; ++idx) {
                    ref_mv_stack[ref_frame][mv_cnt->num_mv_found[ref_frame]]
                        .this_mv = comp_list[idx][0];
                    ref_mv_stack[ref_frame][mv_cnt->num_mv_found[ref_frame]]
                        .comp_mv = comp_list[idx][1];
                    ref_mv_stack[ref_frame][mv_cnt->num_mv_found[ref_frame]]
                        .weight = 2;
                    ++mv_cnt->num_mv_found[ref_frame];
                }
            }
        }
        /*else
        {
            for (int idx = mv_cnt->num_mv_found[ref_frame]; idx < 2; idx++) {
                ref_mv_stack[ref_frame][idx].this_mv = gm_mv_candidates[0];
            }
        }*/
    }

    /* context and clamping process */
    //int num_lists = rf[1] > NONE_FRAME ? 2 : 1;
    //for (int list = 0; list < num_lists; list++) {
    //    for (int idx = 0; idx < mv_cnt->num_mv_found[ref_frame]; idx++) {
    //        IntMv_dec refMv = ref_mv_stack[ref_frame][idx].this_mv;
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
                n4_w << MI_SIZE_LOG2, n4_h << MI_SIZE_LOG2, pi);
            clamp_mv_ref(&ref_mv_stack[ref_frame][idx].comp_mv.as_mv,
                n4_w << MI_SIZE_LOG2, n4_h << MI_SIZE_LOG2, pi);
        }
    }
    else {
        for (int idx = 0; idx < mv_cnt->num_mv_found[ref_frame]; ++idx) {
            clamp_mv_ref(&ref_mv_stack[ref_frame][idx].this_mv.as_mv,
                n4_w << MI_SIZE_LOG2, n4_h << MI_SIZE_LOG2, pi);
        }
    }


    const uint8_t ref_match_count = (mv_cnt->found_above_match > 0) +
        (mv_cnt->found_left_match > 0);
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
        for (int idx = mv_cnt->num_mv_found[ref_frame];
            idx < MAX_MV_REF_CANDIDATES; ++idx)
            mv_ref_list[rf[0]][idx].as_int = gm_mv_candidates[0].as_int;

        for (int idx = 0; idx < AOMMIN(MAX_MV_REF_CANDIDATES,
            mv_cnt->num_mv_found[ref_frame]); ++idx)
        {
            mv_ref_list[rf[0]][idx].as_int =
                ref_mv_stack[ref_frame][idx].this_mv.as_int;
        }
    }
}

static INLINE int16_t svt_mode_context_analyzer(
    const int16_t *const mode_context, const MvReferenceFrame *const rf)
{
    const int8_t ref_frame = av1_ref_frame_type(rf);

    if (rf[1] <= INTRA_FRAME) return mode_context[ref_frame];

    const int16_t newmv_ctx = mode_context[ref_frame] & NEWMV_CTX_MASK;
    const int16_t refmv_ctx =
        (mode_context[ref_frame] >> REFMV_OFFSET) & REFMV_CTX_MASK;

    const int16_t comp_ctx = compound_mode_ctx_map[refmv_ctx >> 1][AOMMIN(
        newmv_ctx, COMP_NEWMV_CTXS - 1)];
    return comp_ctx;
}

void av1_find_mv_refs(EbDecHandle *dec_handle, PartitionInfo_t *pi,
    MvReferenceFrame ref_frame, CandidateMv_dec ref_mv_stack[][MAX_REF_MV_STACK_SIZE],
    IntMv_dec mv_ref_list[][MAX_MV_REF_CANDIDATES], IntMv_dec global_mvs[2],
    int mi_row, int mi_col, int16_t *mode_context, MvCount *mv_cnt)
{
    BlockSize bsize = pi->mi->sb_type;
    MvReferenceFrame rf[2];
    av1_set_ref_frame(rf, ref_frame);

    /* setup global mv process */
    global_mvs[0].as_int = 0;
    global_mvs[1].as_int = 0;
    if (ref_frame != INTRA_FRAME) {
        EbDecPicBuf *buf = dec_handle->cur_pic_buf[0];
        global_mvs[0].as_int = gm_get_motion_vector(&buf->global_motion[rf[0]],
            dec_handle->frame_header.allow_high_precision_mv, bsize,
            mi_col, mi_row, dec_handle->frame_header.force_integer_mv).as_int;

        global_mvs[1].as_int = (rf[1] != NONE_FRAME) ?
            gm_get_motion_vector(&buf->global_motion[rf[1]],
                dec_handle->frame_header.allow_high_precision_mv, bsize,
                mi_col, mi_row, dec_handle->frame_header.force_integer_mv).as_int : 0;
    }
    dec_setup_ref_mv_list(dec_handle, pi, ref_frame, ref_mv_stack, mv_ref_list,
        global_mvs, mi_row, mi_col, mode_context, mv_cnt);
}

static PredictionMode read_inter_compound_mode(EbDecHandle *dec_handle,
    SvtReader *r, int16_t ctx)
{
    ParseCtxt *parse_ctxt = (ParseCtxt *)dec_handle->pv_parse_ctxt;
    const int mode =
        svt_read_symbol(r, parse_ctxt->cur_tile_ctx.inter_compound_mode_cdf[ctx],
            INTER_COMPOUND_MODES, ACCT_STR);
    assert(is_inter_compound_mode(NEAREST_NEARESTMV + mode));
    return NEAREST_NEARESTMV + mode;
}

static INLINE int has_nearmv(PredictionMode mode) {
    return (mode == NEARMV || mode == NEAR_NEARMV ||
        mode == NEAR_NEWMV || mode == NEW_NEARMV);
}

static INLINE uint8_t get_drl_ctx(const CandidateMv_dec *ref_mv_stack, int ref_idx) {

    if (ref_mv_stack[ref_idx].weight >= REF_CAT_LEVEL &&
        ref_mv_stack[ref_idx + 1].weight < REF_CAT_LEVEL)
    {
        return 1;
    }

    if (ref_mv_stack[ref_idx].weight < REF_CAT_LEVEL &&
        ref_mv_stack[ref_idx + 1].weight < REF_CAT_LEVEL)
    {
        return 2;
    }

    return 0;
}

static void read_drl_idx(EbDecHandle *dec_handle, PartitionInfo_t *pi,
    ModeInfo_t *mbmi, SvtReader *r, int num_mv_found)
{
    ParseCtxt *parse_ctxt = (ParseCtxt *)dec_handle->pv_parse_ctxt;
    uint8_t ref_frame_type = av1_ref_frame_type(mbmi->ref_frame);
    mbmi->ref_mv_idx = 0;
    if (mbmi->mode == NEWMV || mbmi->mode == NEW_NEWMV) {
        for (int idx = 0; idx < 2; ++idx) {
            if (num_mv_found > idx + 1) {
                uint8_t drl_ctx = get_drl_ctx(pi->ref_mv_stack[ref_frame_type], idx);
                int drl_idx = svt_read_symbol(r, parse_ctxt->cur_tile_ctx.
                    drl_cdf[drl_ctx], 2, ACCT_STR);
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
                int drl_idx = svt_read_symbol(r, parse_ctxt->cur_tile_ctx.
                    drl_cdf[drl_ctx], 2, ACCT_STR);
                mbmi->ref_mv_idx = idx + drl_idx - 1;
                if (!drl_idx) return;
            }
        }
    }
}

/* TODO: Harmonize*/
static void svt_find_best_ref_mvs(int allow_hp, IntMv_dec *mvlist,
    IntMv_dec *nearest_mv, IntMv_dec *near_mv, int is_integer)
{
    int i;
    // Make sure all the candidates are properly clamped etc
    for (i = 0; i < MAX_MV_REF_CANDIDATES; ++i) {
        lower_mv_precision(&mvlist[i].as_mv, allow_hp, is_integer);
    }
    *nearest_mv = mvlist[0];
    *near_mv = mvlist[1];
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

int read_mv_component(SvtReader *r, NmvComponent *mvcomp, int use_subpel, int usehp) {
    int mag, d, fr, hp;
    const int sign = svt_read_symbol(r, mvcomp->sign_cdf, 2, ACCT_STR);
    const int mv_class =
        svt_read_symbol(r, mvcomp->classes_cdf, MV_CLASSES, ACCT_STR);
    const int class0 = mv_class == MV_CLASS_0;

    // Integer part
    if (class0) {
        d = svt_read_symbol(r, mvcomp->class0_cdf, CLASS0_SIZE, ACCT_STR);
        mag = 0;
    }
    else {
        d = 0;
        for (int i = 0; i < mv_class; ++i)
            d |= svt_read_symbol(r, mvcomp->bits_cdf[i], 2, ACCT_STR) << i;
        mag = CLASS0_SIZE << (mv_class + 2);
    }

    fr = use_subpel ? svt_read_symbol(r, class0 ? mvcomp->class0_fp_cdf[d] :
        mvcomp->fp_cdf, MV_FP_SIZE, ACCT_STR) : 3;

    hp = usehp ? svt_read_symbol(r, class0 ? mvcomp->class0_hp_cdf :
        mvcomp->hp_cdf, 2, ACCT_STR) : 1;

    // Result
    mag += ((d << 3) | (fr << 1) | hp) + 1;
    return sign ? -mag : mag;
}

void read_mv(SvtReader *r, MV *mv, MV *ref, NmvContext *ctx,
    int intra_bc, MvSubpelPrecision precision) {
    MV diff = kZeroMv;

    const MvJointType joint_type =
        (MvJointType)svt_read_symbol(r, &ctx->joints_cdf[intra_bc], MV_JOINTS, ACCT_STR);

    if (mv_joint_vertical(joint_type))
        diff.row = read_mv_component(r, &ctx->comps[0], precision > MV_SUBPEL_NONE,
            precision > MV_SUBPEL_LOW_PRECISION);

    if (mv_joint_horizontal(joint_type))
        diff.col = read_mv_component(r, &ctx->comps[1], precision > MV_SUBPEL_NONE,
            precision > MV_SUBPEL_LOW_PRECISION);

    mv->row = ref->row + diff.row;
    mv->col = ref->col + diff.col;
}

static INLINE int is_mv_valid(MV *mv) {
    return mv->row > MV_LOW && mv->row < MV_UPP && mv->col > MV_LOW &&
        mv->col < MV_UPP;
}

static INLINE int assign_mv(EbDecHandle *dec_handle, PartitionInfo_t *pi,
    IntMv_dec mv[2], IntMv_dec *global_mvs, IntMv_dec ref_mv[2],
    IntMv_dec nearest_mv[2], IntMv_dec near_mv[2],
    int is_compound, int allow_hp, SvtReader *r)
{
    ModeInfo_t *mbmi = pi->mi;
    ParseCtxt *parse_ctxt = (ParseCtxt *)dec_handle->pv_parse_ctxt;

    if (dec_handle->frame_header.force_integer_mv)
        allow_hp = MV_SUBPEL_NONE;

    switch (mbmi->mode) {
    case NEWMV: {
        NmvContext *const nmvc = &parse_ctxt->cur_tile_ctx.nmvc;
        read_mv(r, &mv[0].as_mv, &ref_mv[0].as_mv, nmvc, mbmi->use_intrabc, allow_hp);
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
            read_mv(r, &mv[i].as_mv, &ref_mv[i].as_mv, nmvc, mbmi->use_intrabc, allow_hp);
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
        read_mv(r, &mv[0].as_mv, &ref_mv[0].as_mv, nmvc, mbmi->use_intrabc, allow_hp);
        assert(is_compound);
        mv[1].as_int = nearest_mv[1].as_int;
        break;
    }
    case NEAREST_NEWMV: {
        mv[0].as_int = nearest_mv[0].as_int;
        NmvContext *const nmvc = &parse_ctxt->cur_tile_ctx.nmvc;
        read_mv(r, &mv[1].as_mv, &ref_mv[1].as_mv, nmvc, mbmi->use_intrabc, allow_hp);
        assert(is_compound);
        break;
    }
    case NEAR_NEWMV: {
        mv[0].as_int = near_mv[0].as_int;
        NmvContext *const nmvc = &parse_ctxt->cur_tile_ctx.nmvc;
        read_mv(r, &mv[1].as_mv, &ref_mv[1].as_mv, nmvc, mbmi->use_intrabc, allow_hp);
        assert(is_compound);
        break;
    }
    case NEW_NEARMV: {
        NmvContext *const nmvc = &parse_ctxt->cur_tile_ctx.nmvc;
        read_mv(r, &mv[0].as_mv, &ref_mv[0].as_mv, nmvc, mbmi->use_intrabc, allow_hp);
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
    default: { return 0; }
    }

    int ret = is_mv_valid(&mv[0].as_mv);
    if (is_compound)
        ret = ret && is_mv_valid(&mv[1].as_mv);
    return ret;
}

static INLINE int is_interintra_wedge_used(BlockSize sb_type) {
    return wedge_params_lookup[sb_type].bits > 0;
}

void read_interintra_mode(EbDecHandle *dec_handle,
    ModeInfo_t *mbmi, SvtReader *r)
{
    ParseCtxt *parse_ctxt = (ParseCtxt *)dec_handle->pv_parse_ctxt;
    FRAME_CONTEXT *frm_ctx = &parse_ctxt->cur_tile_ctx;
    BlockSize bsize = mbmi->sb_type;
    if (dec_handle->seq_header.enable_interintra_compound
        && !mbmi->skip_mode && is_interintra_allowed(mbmi))
    {
        const int bsize_group = size_group_lookup[bsize];
        mbmi->is_inter_intra =
            svt_read_symbol(r, frm_ctx->interintra_cdf[bsize_group], 2, ACCT_STR);
        assert(mbmi->ref_frame[1] == NONE_FRAME);
        if (mbmi->is_inter_intra) {
            mbmi->interintra_mode.interintra_mode = (InterIntraMode)svt_read_symbol(
                r, frm_ctx->interintra_mode_cdf[bsize_group], INTERINTRA_MODES,
                ACCT_STR);
            mbmi->ref_frame[1] = INTRA_FRAME;
            mbmi->angle_delta[PLANE_TYPE_Y] = 0;
            mbmi->angle_delta[PLANE_TYPE_UV] = 0;
            mbmi->filter_intra_mode_info.use_filter_intra = 0;
            if (is_interintra_wedge_used(bsize)) {
                mbmi->interintra_mode.wedge_interintra = svt_read_symbol(
                    r, frm_ctx->wedge_interintra_cdf[bsize], 2, ACCT_STR);
                if (mbmi->interintra_mode.wedge_interintra) {
                    mbmi->interintra_mode.interintra_wedge_index =
                        svt_read_symbol(r, frm_ctx->wedge_idx_cdf[bsize], 16, ACCT_STR);
                }
            }
        }
    }
}

static INLINE void add_samples(ModeInfo_t *mbmi, int *pts, int *pts_inref,
    int row_offset, int sign_r, int col_offset, int sign_c)
{
    int bw = block_size_wide[mbmi->sb_type];
    int bh = block_size_high[mbmi->sb_type];
    int x = col_offset * MI_SIZE + sign_c * AOMMAX(bw, MI_SIZE) / 2 - 1;
    int y = row_offset * MI_SIZE + sign_r * AOMMAX(bh, MI_SIZE) / 2 - 1;

    pts[0] = (x * 8);
    pts[1] = (y * 8);
    pts_inref[0] = (x * 8) + mbmi->mv[0].as_mv.col;
    pts_inref[1] = (y * 8) + mbmi->mv[0].as_mv.row;
}

int find_warp_samples(EbDecHandle *dec_handle, PartitionInfo_t *pi,
    int mi_row, int mi_col, int *pts, int *pts_inref)
{
    ModeInfo_t *const mbmi0 = pi->mi;
    int ref_frame = mbmi0->ref_frame[0];
    int up_available = pi->up_available;
    int left_available = pi->left_available;
    int i, mi_step = 1, np = 0;

    ParseCtxt *parse_ctx = (ParseCtxt*)dec_handle->pv_parse_ctxt;
    TileInfo *tile = &parse_ctx->cur_tile_info;
    int do_tl = 1;
    int do_tr = 1;
    int b4_w = mi_size_wide[pi->mi->sb_type];
    int b4_h = mi_size_high[pi->mi->sb_type];

    // scan the nearest above rows
    if (up_available) {
        ModeInfo_t *mbmi = get_cur_mode_info(dec_handle,
            mi_row - 1, mi_col, pi->sb_info);
        uint8_t n4_w = mi_size_wide[mbmi->sb_type];

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
        }
        else {
            // Handle "current block width > above block width" case.
            for (i = 0; i < AOMMIN(b4_w, tile->mi_col_end - mi_col); i += mi_step) {
                mbmi = get_cur_mode_info(dec_handle,
                    mi_row - 1, mi_col + i, pi->sb_info);
                n4_w = mi_size_wide[mbmi->sb_type];
                mi_step = AOMMIN(b4_w, n4_w);

                if (mbmi->ref_frame[0] == ref_frame &&
                    mbmi->ref_frame[1] == NONE_FRAME) {
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

        ModeInfo_t *mbmi = get_cur_mode_info(dec_handle,
            mi_row, mi_col - 1, pi->sb_info);
        uint8_t n4_h = mi_size_high[mbmi->sb_type];

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
        }
        else {
            // Handle "current block height > above block height" case.
            for (i = 0; i < AOMMIN(b4_h, tile->mi_row_end - mi_row); i += mi_step) {
                mbmi = get_cur_mode_info(dec_handle,
                    mi_row + i, mi_col - 1, pi->sb_info);
                n4_h = mi_size_high[mbmi->sb_type];
                mi_step = AOMMIN(b4_h, n4_h);

                if (mbmi->ref_frame[0] == ref_frame &&
                    mbmi->ref_frame[1] == NONE_FRAME) {
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

        ModeInfo_t *mbmi = get_cur_mode_info(dec_handle,
            mi_row - 1, mi_col - 1, pi->sb_info);

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
    if (do_tr &&
        has_top_right(dec_handle, pi, mi_row, mi_col, AOMMAX(b4_w, b4_h))) {

        int mv_row = mi_row - 1;
        int mv_col = mi_col + b4_w;

        if (is_inside(tile, mv_col, mv_row)) {
            ModeInfo_t *mbmi =
                get_cur_mode_info(dec_handle,
                    mv_row, mv_col, pi->sb_info);

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

int has_overlappable_cand(EbDecHandle *dec_handle, PartitionInfo_t *pi,
    int mi_row, int mi_col)
{
    ParseCtxt *parse_ctx = (ParseCtxt*)dec_handle->pv_parse_ctxt;
    const TileInfo *const tile = &parse_ctx->cur_tile_info;
    ModeInfo_t *mbmi = pi->mi;
    if (!is_motion_variation_allowed_bsize(mbmi->sb_type)) return 0;

    if (pi->up_available) {
        int w4 = mi_size_wide[mbmi->sb_type];
        int x4 = mi_col;
        while (x4 < AOMMIN(tile->mi_col_end, mi_col + w4)) {
            ModeInfo_t *top_nb_mode = get_cur_mode_info(dec_handle,
                mi_row - 1, x4 | 1, pi->sb_info);
            x4 += AOMMAX(2, mi_size_wide[top_nb_mode->sb_type] >> 2);
            if (dec_is_inter_block(top_nb_mode))
                return 1;
        }
    }
    if (pi->left_available) {
        int h4 = mi_size_high[mbmi->sb_type];
        int y4 = mi_row;
        while (y4 < AOMMIN(tile->mi_row_end, mi_row + h4)) {
            ModeInfo_t *left_nb_mode = get_cur_mode_info(dec_handle,
                y4 | 1, mi_col - 1, pi->sb_info);
            y4 += AOMMAX(2, mi_size_high[left_nb_mode->sb_type] >> 2);
            if (dec_is_inter_block(left_nb_mode))
                return 1;
        }
    }
    return 0;
}

static INLINE MotionMode is_motion_mode_allowed(EbDecHandle *dec_handle, GlobalMotionParams *gm_params, PartitionInfo_t *pi, int mi_row,
    int mi_col, int allow_warped_motion)
{
    ModeInfo_t *mbmi = pi->mi;
    if (dec_handle->frame_header.force_integer_mv == 0) {
        const TransformationType gm_type = gm_params[mbmi->ref_frame[0]].gm_type;
        if (is_global_mv_block(mbmi, gm_type)) return SIMPLE_TRANSLATION;
    }
    if ((block_size_wide[mbmi->sb_type] >= 8 && block_size_high[mbmi->sb_type] >= 8) &&
        (mbmi->mode >= NEARESTMV && mbmi->mode < MB_MODE_COUNT)
        && mbmi->ref_frame[1] != INTRA_FRAME && !has_second_ref(mbmi)) {
        if (!has_overlappable_cand(dec_handle, pi, mi_row, mi_col))
            return SIMPLE_TRANSLATION;
        assert(!has_second_ref(mbmi));

        if (pi->num_samples >= 1 &&
            (allow_warped_motion /*&&
                !av1_is_scaled(xd->block_ref_scale_factors[0])*/)) {
            if (dec_handle->frame_header.force_integer_mv) {
                return OBMC_CAUSAL;
            }
            return WARPED_CAUSAL;
        }
        return OBMC_CAUSAL;
    }
    else {
        return SIMPLE_TRANSLATION;
    }
}

MotionMode read_motion_mode(EbDecHandle *dec_handle,
    PartitionInfo_t *pi, int mi_row, int mi_col, SvtReader *r)
{
    ParseCtxt *parse_ctxt = (ParseCtxt *)dec_handle->pv_parse_ctxt;
    FRAME_CONTEXT *frm_ctx = &parse_ctxt->cur_tile_ctx;
    FrameHeader *frame_info = &dec_handle->frame_header;
    int allow_warped_motion = frame_info->allow_warped_motion;
    ModeInfo_t *mbmi = pi->mi;

    if (dec_handle->frame_header.is_motion_mode_switchable == 0) return SIMPLE_TRANSLATION;
    if (mbmi->skip_mode) return SIMPLE_TRANSLATION;

    const MotionMode last_motion_mode_allowed =
        is_motion_mode_allowed(dec_handle,
            dec_handle->cur_pic_buf[0]->global_motion, pi,
            mi_row, mi_col, allow_warped_motion);
    int motion_mode;

    if (last_motion_mode_allowed == SIMPLE_TRANSLATION) return SIMPLE_TRANSLATION;

    if (last_motion_mode_allowed == OBMC_CAUSAL) {
        motion_mode =
            svt_read_symbol(r, frm_ctx->obmc_cdf[mbmi->sb_type], 2, ACCT_STR);
        return (MotionMode)(motion_mode);
    }
    else {
        motion_mode =
            svt_read_symbol(r, frm_ctx->motion_mode_cdf[mbmi->sb_type],
                MOTION_MODES, ACCT_STR);
        return (MotionMode)(motion_mode);
    }
}


static INLINE int is_comp_ref_allowed(BlockSize bsize) {
    return AOMMIN(block_size_wide[bsize], block_size_high[bsize]) >= 8;
}

static INLINE int is_masked_compound_type(CompoundType type) {
    return (type == COMPOUND_WEDGE || type == COMPOUND_DIFFWTD);
}

static INLINE int is_interinter_compound_used(CompoundType type, BlockSize sb_type) {
    const int comp_allowed = is_comp_ref_allowed(sb_type);
    switch (type) {
    case COMPOUND_AVERAGE:
    case COMPOUND_DISTWTD:
    case COMPOUND_DIFFWTD: return comp_allowed;
    case COMPOUND_WEDGE:
        return comp_allowed && wedge_params_lookup[sb_type].bits > 0;
    default: assert(0); return 0;
    }
}

static INLINE int is_any_masked_compound_used(BlockSize sb_type) {
    CompoundType comp_type;
    int i;
    if (!is_comp_ref_allowed(sb_type)) return 0;
    for (i = 0; i < COMPOUND_TYPES; i++) {
        comp_type = (CompoundType)i;
        if (is_masked_compound_type(comp_type) &&
            is_interinter_compound_used(comp_type, sb_type))
            return 1;
    }
    return 0;
}

static INLINE int get_comp_group_idx_context(const PartitionInfo_t *xd) {
    const ModeInfo_t *const above_mi = xd->above_mbmi;
    const ModeInfo_t *const left_mi = xd->left_mbmi;
    int above_ctx = 0, left_ctx = 0;

    if (above_mi) {
        if (has_second_ref(above_mi))
            above_ctx = above_mi->inter_compound.comp_group_idx;
        else if (above_mi->ref_frame[0] == ALTREF_FRAME)
            above_ctx = 3;
    }
    if (left_mi) {
        if (has_second_ref(left_mi))
            left_ctx = left_mi->inter_compound.comp_group_idx;
        else if (left_mi->ref_frame[0] == ALTREF_FRAME)
            left_ctx = 3;
    }

    return AOMMIN(5, above_ctx + left_ctx);
}

int get_comp_index_context(EbDecHandle *dec_handle, PartitionInfo_t *pi) {
    ModeInfo_t *mbmi = pi->mi;
    SeqHeader *seq_params = &dec_handle->seq_header;
    FrameHeader *frm_header = &dec_handle->frame_header;

    int bck_frame_index = 0, fwd_frame_index = 0;
    int cur_frame_index = frm_header->order_hint;

    EbDecPicBuf *bck_buf = get_ref_frame_buf(dec_handle, mbmi->ref_frame[0]);
    EbDecPicBuf *fwd_buf = get_ref_frame_buf(dec_handle, mbmi->ref_frame[1]);

    if (bck_buf != NULL) bck_frame_index = bck_buf->order_hint;
    if (fwd_buf != NULL) fwd_frame_index = fwd_buf->order_hint;

    int fwd = abs(get_relative_dist(&seq_params->order_hint_info,
        fwd_frame_index, cur_frame_index));
    int bck = abs(get_relative_dist(&seq_params->order_hint_info,
        cur_frame_index, bck_frame_index));

    const ModeInfo_t *const above_mi = pi->above_mbmi;
    const ModeInfo_t *const left_mi = pi->left_mbmi;

    int above_ctx = 0, left_ctx = 0;
    const int offset = (fwd == bck);

    if (above_mi != NULL) {
        if (has_second_ref(above_mi))
            above_ctx = above_mi->inter_compound.compound_idx;
        else if (above_mi->ref_frame[0] == ALTREF_FRAME)
            above_ctx = 1;
    }

    if (left_mi != NULL) {
        if (has_second_ref(left_mi))
            left_ctx = left_mi->inter_compound.compound_idx;
        else if (left_mi->ref_frame[0] == ALTREF_FRAME)
            left_ctx = 1;
    }

    return above_ctx + left_ctx + 3 * offset;
}

void read_compound_type(EbDecHandle *dec_handle, PartitionInfo_t *pi, SvtReader *r)
{
    ModeInfo_t *mbmi = pi->mi;
    BlockSize bsize = mbmi->sb_type;
    mbmi->inter_compound.comp_group_idx = 0;
    mbmi->inter_compound.compound_idx = 1;
    ParseCtxt *parse_ctxt = (ParseCtxt *)dec_handle->pv_parse_ctxt;
    FRAME_CONTEXT *frm_ctx = &parse_ctxt->cur_tile_ctx;

    if (mbmi->skip_mode) mbmi->inter_compound.type = COMPOUND_AVERAGE;

    if (has_second_ref(mbmi) && !mbmi->skip_mode) {
        // Read idx to indicate current compound inter prediction mode group
        const int masked_compound_used = is_any_masked_compound_used(bsize) &&
            dec_handle->seq_header.enable_masked_compound;

        if (masked_compound_used) {
            const int ctx_comp_group_idx = get_comp_group_idx_context(pi);
            mbmi->inter_compound.comp_group_idx = svt_read_symbol(
                r, frm_ctx->comp_group_idx_cdf[ctx_comp_group_idx], 2, ACCT_STR);
        }

        if (mbmi->inter_compound.comp_group_idx == 0) {
            if (dec_handle->seq_header.order_hint_info.enable_jnt_comp) {
                const int comp_index_ctx = get_comp_index_context(dec_handle, pi);
                mbmi->inter_compound.compound_idx = svt_read_symbol(
                    r, frm_ctx->compound_index_cdf[comp_index_ctx], 2, ACCT_STR);
                mbmi->inter_compound.type =
                    mbmi->inter_compound.compound_idx ? COMPOUND_AVERAGE : COMPOUND_DISTWTD;
            }
            else {
                // Distance-weighted compound is disabled, so always use average
                mbmi->inter_compound.compound_idx = 1;
                mbmi->inter_compound.type = COMPOUND_AVERAGE;
            }
        }
        else {
            assert(dec_handle->frame_header.reference_mode != SINGLE_REFERENCE &&
                is_inter_compound_mode(mbmi->mode) &&
                mbmi->motion_mode == SIMPLE_TRANSLATION);
            assert(masked_compound_used);

            // compound_diffwtd, wedge
            if (is_interinter_compound_used(COMPOUND_WEDGE, bsize))
                mbmi->inter_compound.type = COMPOUND_WEDGE +
                svt_read_symbol(r, frm_ctx->compound_type_cdf[bsize],
                    MASKED_COMPOUND_TYPES, ACCT_STR);
            else
                mbmi->inter_compound.type = COMPOUND_DIFFWTD;

            if (mbmi->inter_compound.type == COMPOUND_WEDGE) {
                assert(is_interinter_compound_used(COMPOUND_WEDGE, bsize));
                mbmi->inter_compound.wedge_index =
                    svt_read_symbol(r, frm_ctx->wedge_idx_cdf[bsize], 16, ACCT_STR);
                mbmi->inter_compound.wedge_sign = svt_read_bit(r, ACCT_STR);
            }
            else {
                assert(mbmi->inter_compound.type == COMPOUND_DIFFWTD);
                mbmi->inter_compound.mask_type =
                    svt_read_literal(r, MAX_DIFFWTD_MASK_BITS, ACCT_STR);
            }
        }
    }
   /* else {
        if (mbmi->is_inter_intra)
            mbmi->inter_compound.type = mbmi->interintra_mode.wedge_interintra ?
            COMPOUND_WEDGE : COMPOUND_INTRA;
        else
            mbmi->inter_compound.type = COMPOUND_AVERAGE;
    }*/
}

static INLINE int is_nontrans_global_motion(PartitionInfo_t *pi,
    GlobalMotionParams *gm_params)
{
    int ref;
    ModeInfo_t * mbmi = pi->mi;
    // First check if all modes are GLOBALMV
    if (mbmi->mode != GLOBALMV && mbmi->mode != GLOBAL_GLOBALMV) return 0;

    if (AOMMIN(mi_size_wide[mbmi->sb_type], mi_size_high[mbmi->sb_type]) < 2)
        return 0;

    // Now check if all global motion is non translational
    for (ref = 0; ref < 1 + has_second_ref(mbmi); ++ref) {
        if (gm_params[mbmi->ref_frame[ref]].gm_type == TRANSLATION) return 0;
    }
    return 1;
}

static INLINE int av1_is_interp_needed(PartitionInfo_t *pi,
    GlobalMotionParams *gm_params)
{
    ModeInfo_t * mbmi = pi->mi;
    if (mbmi->skip_mode) return 0;
    if (mbmi->motion_mode == WARPED_CAUSAL) return 0;
    if (is_nontrans_global_motion(pi, gm_params)) return 0;
    return 1;
}

static INLINE InterpFilters make_interp_filters(InterpFilter y_filter,
    InterpFilter x_filter) {
    uint16_t y16 = y_filter & 0xffff;
    uint16_t x16 = x_filter & 0xffff;
    return y16 | ((uint32_t)x16 << 16);
}

static INLINE InterpFilters broadcast_interp_filter(InterpFilter filter) {
    return make_interp_filters(filter, filter);
}

static INLINE InterpFilter unswitchable_filter(InterpFilter filter) {
    return filter == SWITCHABLE ? EIGHTTAP_REGULAR : filter;
}

static INLINE void set_default_interp_filters(
    ModeInfo_t *mbmi, InterpFilter frame_interp_filter) {
    mbmi->interp_filters =
        broadcast_interp_filter(unswitchable_filter(frame_interp_filter));
}

static InterpFilter get_ref_filter_type(const ModeInfo_t *ref_mbmi,
    int dir, MvReferenceFrame ref_frame)
{
    return ((ref_mbmi->ref_frame[0] == ref_frame ||
        ref_mbmi->ref_frame[1] == ref_frame)
        ? av1_extract_interp_filter(ref_mbmi->interp_filters, dir & 0x01)
        : SWITCHABLE_FILTERS);
}

int get_context_interp(PartitionInfo_t *pi, int dir) {
    const ModeInfo_t *const mbmi = pi->mi;
    const int ctx_offset =
        (mbmi->ref_frame[1] > INTRA_FRAME) * INTER_FILTER_COMP_OFFSET;
    assert(dir == 0 || dir == 1);
    const MvReferenceFrame ref_frame = mbmi->ref_frame[0];

    int filter_type_ctx = ctx_offset + (dir & 0x01) * INTER_FILTER_DIR_OFFSET;
    int left_type = SWITCHABLE_FILTERS;
    int above_type = SWITCHABLE_FILTERS;

    if (pi->left_available)
        left_type = get_ref_filter_type(pi->left_mbmi, dir, ref_frame);

    if (pi->up_available)
        above_type =
        get_ref_filter_type(pi->above_mbmi, dir, ref_frame);

    if (left_type == above_type) {
        filter_type_ctx += left_type;
    }
    else if (left_type == SWITCHABLE_FILTERS) {
        assert(above_type != SWITCHABLE_FILTERS);
        filter_type_ctx += above_type;
    }
    else if (above_type == SWITCHABLE_FILTERS) {
        assert(left_type != SWITCHABLE_FILTERS);
        filter_type_ctx += left_type;
    }
    else {
        filter_type_ctx += SWITCHABLE_FILTERS;
    }

    return filter_type_ctx;
}

void inter_block_mode_info(EbDecHandle *dec_handle, PartitionInfo_t* pi,
    int mi_row, int mi_col, SvtReader *r)
{
    ModeInfo_t *mbmi = pi->mi;
    const int allow_hp = dec_handle->frame_header.allow_high_precision_mv;
    IntMv_dec ref_mvs[MODE_CTX_REF_FRAMES][MAX_MV_REF_CANDIDATES] = { { { 0 } } };
    int16_t inter_mode_ctx[MODE_CTX_REF_FRAMES];
    int pts[SAMPLES_ARRAY_SIZE], pts_inref[SAMPLES_ARRAY_SIZE];
    SegmentationParams *seg = &dec_handle->frame_header.segmentation_params;
    ParseCtxt *parse_ctxt = (ParseCtxt *)dec_handle->pv_parse_ctxt;
    MvCount mv_cnt;

    /* TO-DO initialize palette info */

    svt_collect_neighbors_ref_counts(pi);

    read_ref_frames(dec_handle, pi, r);
   /* if ((pi->mi->ref_frame[0] >= BWDREF_FRAME && pi->mi->ref_frame[0] <= ALTREF_FRAME) ||
        (pi->mi->ref_frame[1] >= BWDREF_FRAME && pi->mi->ref_frame[1] <= ALTREF_FRAME)) {
        printf("ALTREF found - frame : %d\n", dec_handle->dec_cnt);
        exit(0);
    }*/
    const int is_compound = has_second_ref(mbmi);

    MvReferenceFrame ref_frame = av1_ref_frame_type(mbmi->ref_frame);
    IntMv_dec global_mvs[2];
    av1_find_mv_refs(dec_handle, pi, ref_frame, pi->ref_mv_stack,
        ref_mvs, global_mvs, mi_row, mi_col,
        inter_mode_ctx, &mv_cnt);

#if EXTRA_DUMP
    if (enable_dump) {
        printf("\n mi_row: %d mi_col: %d\n", mi_row, mi_col);
        /*for (int i = 0; i < MODE_CTX_REF_FRAMES; i++)
            for (int j = 0; j < MAX_REF_MV_STACK_SIZE; j++)
                printf("ref_mv_stack[%d][%d]=%d\t", i, j, pi->ref_mv_stack[i][j].this_mv.as_int);
        printf("\n");*/
        fflush(stdout);
    }
#endif

    int mode_ctx = svt_mode_context_analyzer(inter_mode_ctx, mbmi->ref_frame);
    mbmi->ref_mv_idx = 0;

    if (mbmi->skip_mode) {
        assert(is_compound);
        mbmi->mode = NEAREST_NEARESTMV;
    }
    else {
        if (seg_feature_active(seg, mbmi->segment_id, SEG_LVL_SKIP) ||
            seg_feature_active(seg, mbmi->segment_id, SEG_LVL_GLOBALMV))
            mbmi->mode = GLOBALMV;
        else {
            if (is_compound)
                mbmi->mode = read_inter_compound_mode(dec_handle, r, mode_ctx);
            else {
                int new_mv = svt_read_symbol(r, parse_ctxt->cur_tile_ctx.
                    newmv_cdf[mode_ctx & NEWMV_CTX_MASK], 2, ACCT_STR);
                if (new_mv) {
                    int zero_mv = svt_read_symbol(r,
                        parse_ctxt->cur_tile_ctx.zeromv_cdf
                        [(mode_ctx >> GLOBALMV_OFFSET) & GLOBALMV_CTX_MASK],
                        2, ACCT_STR);
                    if (zero_mv) {
                        int ref_mv = svt_read_symbol(r, parse_ctxt->cur_tile_ctx.
                            refmv_cdf[(mode_ctx >> REFMV_OFFSET) & REFMV_CTX_MASK],
                            2, ACCT_STR);
                        mbmi->mode = ref_mv ? NEARMV : NEARESTMV;
                    }
                    else
                        mbmi->mode = GLOBALMV;
                }
                else
                    mbmi->mode = NEWMV;
            }
            if (mbmi->mode == NEWMV || mbmi->mode == NEW_NEWMV ||
                has_nearmv(mbmi->mode))
                read_drl_idx(dec_handle, pi, mbmi, r, mv_cnt.num_mv_found[ref_frame]);
        }
    }
    mbmi->uv_mode = UV_DC_PRED;

    IntMv_dec ref_mv[2];
    IntMv_dec nearestmv[2], nearmv[2];
    if (!is_compound && mbmi->mode != GLOBALMV) {
        svt_find_best_ref_mvs(allow_hp, ref_mvs[mbmi->ref_frame[0]], &nearestmv[0],
            &nearmv[0], dec_handle->frame_header.force_integer_mv);
    }
    if (is_compound && mbmi->mode != GLOBAL_GLOBALMV) {
        int ref_mv_idx = mbmi->ref_mv_idx + 1;
        nearestmv[0] = pi->ref_mv_stack[ref_frame][0].this_mv;
        nearestmv[1] = pi->ref_mv_stack[ref_frame][0].comp_mv;
        nearmv[0] = pi->ref_mv_stack[ref_frame][ref_mv_idx].this_mv;
        nearmv[1] = pi->ref_mv_stack[ref_frame][ref_mv_idx].comp_mv;
        lower_mv_precision(&nearestmv[0].as_mv, allow_hp,
            dec_handle->frame_header.force_integer_mv);
        lower_mv_precision(&nearestmv[1].as_mv, allow_hp,
            dec_handle->frame_header.force_integer_mv);
        lower_mv_precision(&nearmv[0].as_mv, allow_hp,
            dec_handle->frame_header.force_integer_mv);
        lower_mv_precision(&nearmv[1].as_mv, allow_hp,
            dec_handle->frame_header.force_integer_mv);
    }
    else if (mbmi->ref_mv_idx > 0 && mbmi->mode == NEARMV) {
        IntMv_dec cur_mv =
            pi->ref_mv_stack[mbmi->ref_frame[0]][1 + mbmi->ref_mv_idx].this_mv;
        nearmv[0] = cur_mv;
    }

    ref_mv[0] = nearestmv[0];
    ref_mv[1] = nearestmv[1];

    if (is_compound) {
        int ref_mv_idx = mbmi->ref_mv_idx;
        // Special case: NEAR_NEWMV and NEW_NEARMV modes use
        // 1 + mbmi->ref_mv_idx (like NEARMV) instead of
        // mbmi->ref_mv_idx (like NEWMV)
        if (mbmi->mode == NEAR_NEWMV || mbmi->mode == NEW_NEARMV)
            ref_mv_idx = 1 + mbmi->ref_mv_idx;

        if (compound_ref0_mode(mbmi->mode) == NEWMV)
            ref_mv[0] = pi->ref_mv_stack[ref_frame][ref_mv_idx].this_mv;

        if (compound_ref1_mode(mbmi->mode) == NEWMV)
            ref_mv[1] = pi->ref_mv_stack[ref_frame][ref_mv_idx].comp_mv;
    }
    else {
        if (mbmi->mode == NEWMV) {
            if (mv_cnt.num_mv_found[ref_frame] > 1)
                ref_mv[0] = pi->ref_mv_stack[ref_frame][mbmi->ref_mv_idx].this_mv;
        }
    }

    assign_mv(dec_handle, pi, mbmi->mv, global_mvs,
        ref_mv, nearestmv, nearmv, is_compound, allow_hp, r);

#if EXTRA_DUMP
    if (enable_dump) {
        printf("\n mode %d MV %d %d \n", mbmi->mode, mbmi->mv[0].as_mv.row, mbmi->mv[0].as_mv.col);
        fflush(stdout);
    }
#endif
    read_interintra_mode(dec_handle, mbmi, r);

    pi->num_samples = find_warp_samples(dec_handle, pi, mi_row, mi_col, pts, pts_inref);

    mbmi->motion_mode = read_motion_mode(dec_handle, pi, mi_row, mi_col, r);

    read_compound_type(dec_handle, pi, r);

    if (!av1_is_interp_needed(pi, dec_handle->cur_pic_buf[0]->global_motion)) {
        set_default_interp_filters(mbmi,
            dec_handle->frame_header.interpolation_filter);
    }
    else {
        if (dec_handle->frame_header.interpolation_filter != SWITCHABLE) {
            mbmi->interp_filters =
                av1_broadcast_interp_filter(dec_handle->frame_header.interpolation_filter);
        }
        else {
            InterpFilter ref0_filter[2] = { EIGHTTAP_REGULAR, EIGHTTAP_REGULAR };
            for (int dir = 0; dir < 2; ++dir) {
                const int ctx = get_context_interp(pi, dir);
                ref0_filter[dir] = (InterpFilter)svt_read_symbol(
                    r, parse_ctxt->cur_tile_ctx.switchable_interp_cdf[ctx], SWITCHABLE_FILTERS, ACCT_STR);
                if (dec_handle->seq_header.enable_dual_filter == 0) {
                    ref0_filter[1] = ref0_filter[0];
                    break;
                }
            }
            mbmi->interp_filters =
                av1_make_interp_filters(ref0_filter[0], ref0_filter[1]);
        }
    }
}
