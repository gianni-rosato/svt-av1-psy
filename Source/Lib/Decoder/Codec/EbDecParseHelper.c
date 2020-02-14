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

#include "EbDefinitions.h"
#include "EbSvtAv1Dec.h"
#include "EbDecHandle.h"

#include "EbDecParseHelper.h"
#include "EbDecUtils.h"
#include "EbInvTransforms.h"
#include "EbInterPrediction.h"

int neg_deinterleave(const int diff, int ref, int max) {
    if (!ref) return diff;
    if (ref >= max - 1) return max - diff - 1;
    if (2 * ref < max) {
        if (diff <= 2 * ref) {
            if (diff & 1)
                return ref + ((diff + 1) >> 1);
            else
                return ref - (diff >> 1);
        }
        return diff;
    } else {
        if (diff <= 2 * (max - ref - 1)) {
            if (diff & 1)
                return ref + ((diff + 1) >> 1);
            else
                return ref - (diff >> 1);
        }
        return max - (diff + 1);
    }
}

void set_segment_id(EbDecHandle *dec_handle, int mi_offset, int x_mis, int y_mis, int segment_id) {
    assert(segment_id >= 0 && segment_id < MAX_SEGMENTS);
    FrameHeader *frm_header = &dec_handle->frame_header;

    for (int y = 0; y < y_mis; y++)
        for (int x = 0; x < x_mis; x++)
            dec_handle->cur_pic_buf[0]->segment_maps[mi_offset + y * frm_header->mi_cols + x] =
                segment_id;
}

static INLINE int get_tx_size_context(const PartitionInfo *xd, ParseCtxt *parse_ctxt) {
    const BlockModeInfo *      mbmi        = xd->mi;
    const BlockModeInfo *const above_mbmi  = xd->above_mbmi;
    const BlockModeInfo *const left_mbmi   = xd->left_mbmi;
    const TxSize               max_tx_size = max_txsize_rect_lookup[mbmi->sb_type];
    const uint8_t              max_tx_wide = tx_size_wide[max_tx_size];
    const uint8_t              max_tx_high = tx_size_high[max_tx_size];
    const int                  has_above   = xd->up_available;
    const int                  has_left    = xd->left_available;

    int above =
        parse_ctxt->parse_above_nbr4x4_ctxt
            ->above_tx_wd[xd->mi_col - parse_ctxt->cur_tile_info.mi_col_start] >= max_tx_wide;
    int left = parse_ctxt->parse_left_nbr4x4_ctxt->left_tx_ht[xd->mi_row - parse_ctxt->sb_row_mi] >=
               max_tx_high;

    if (has_above)
        if (is_inter_block(above_mbmi)) above = block_size_wide[above_mbmi->sb_type] >= max_tx_wide;

    if (has_left)
        if (is_inter_block(left_mbmi)) left = block_size_high[left_mbmi->sb_type] >= max_tx_high;

    if (has_above && has_left)
        return (above + left);
    else if (has_above)
        return above;
    else if (has_left)
        return left;
    else
        return 0;
}

static INLINE TxSize depth_to_tx_size(int depth, BlockSize bsize) {
    TxSize tx_size = max_txsize_rect_lookup[bsize];
    for (int d = 0; d < depth; ++d) tx_size = sub_tx_size_map[tx_size];
    return tx_size;
}

void update_tx_context(ParseCtxt *parse_ctxt, PartitionInfo *pi, BlockSize bsize, TxSize tx_size,
                       int blk_row, int blk_col) {
    int                   mi_row          = pi->mi_row;
    int                   mi_col          = pi->mi_col;
    ParseAboveNbr4x4Ctxt *above_parse_ctx = parse_ctxt->parse_above_nbr4x4_ctxt;
    ParseLeftNbr4x4Ctxt * left_parse_ctx  = parse_ctxt->parse_left_nbr4x4_ctxt;
    BlockSize             b_size          = bsize;
    if (is_inter_block(pi->mi)) b_size = txsize_to_bsize[tx_size];
    uint8_t *const above_ctx =
        above_parse_ctx->above_tx_wd + (mi_col - parse_ctxt->cur_tile_info.mi_col_start + blk_col);
    uint8_t *const left_ctx =
        left_parse_ctx->left_tx_ht + (mi_row - parse_ctxt->sb_row_mi + blk_row);

    const uint8_t tx_wide = tx_size_wide[tx_size];
    const uint8_t tx_high = tx_size_high[tx_size];

    const int bw = mi_size_wide[b_size];
    const int bh = mi_size_high[b_size];
    memset(above_ctx, tx_wide, bw);
    memset(left_ctx, tx_high, bh);
}

TxSize read_selected_tx_size(PartitionInfo *xd, ParseCtxt *parse_ctxt) {
    SvtReader *     r            = &parse_ctxt->r;
    const BlockSize bsize        = xd->mi->sb_type;
    const int32_t   tx_size_cat  = bsize_to_tx_size_cat(bsize);
    const int       max_tx_depth = bsize_to_max_depth(bsize);
    const int       ctx          = get_tx_size_context(xd, parse_ctxt);
    /*TODO : Change idx */
    const int depth = svt_read_symbol(
        r, parse_ctxt->cur_tile_ctx.tx_size_cdf[tx_size_cat][ctx], max_tx_depth + 1, ACCT_STR);
    assert(depth >= 0 && depth <= max_tx_depth);
    const TxSize tx_size = depth_to_tx_size(depth, bsize);
    return tx_size;
}

int get_intra_inter_context(PartitionInfo *xd) {
    const BlockModeInfo *const above_mbmi = xd->above_mbmi;
    const BlockModeInfo *const left_mbmi  = xd->left_mbmi;
    const int                  has_above  = xd->up_available;
    const int                  has_left   = xd->left_available;

    if (has_above && has_left) { // both edges available
        const int above_intra = !is_inter_block(above_mbmi);
        const int left_intra  = !is_inter_block(left_mbmi);
        return left_intra && above_intra ? 3 : left_intra || above_intra;
    } else if (has_above || has_left) { // one edge available
        return 2 * !is_inter_block(has_above ? above_mbmi : left_mbmi);
    } else
        return 0;
}

PredictionMode read_intra_mode(SvtReader *r, AomCdfProb *cdf) {
    return (PredictionMode)svt_read_symbol(r, cdf, INTRA_MODES, ACCT_STR);
}

UvPredictionMode read_intra_mode_uv(FRAME_CONTEXT *ec_ctx, SvtReader *r, CflAllowedType cfl_allowed,
                                    PredictionMode y_mode) {
    const UvPredictionMode uv_mode = svt_read_symbol(
        r, ec_ctx->uv_mode_cdf[cfl_allowed][y_mode], UV_INTRA_MODES - !cfl_allowed, ACCT_STR);
    return uv_mode;
}

static INLINE int block_center_x(int mi_col, BlockSize bs) {
    const int bw = block_size_wide[bs];
    return mi_col * MI_SIZE + bw / 2 - 1;
}

static INLINE int block_center_y(int mi_row, BlockSize bs) {
    const int bh = block_size_high[bs];
    return mi_row * MI_SIZE + bh / 2 - 1;
}

static INLINE int convert_to_trans_prec(int allow_hp, int coor) {
    if (allow_hp)
        return ROUND_POWER_OF_TWO_SIGNED(coor, WARPEDMODEL_PREC_BITS - 3);
    else
        return ROUND_POWER_OF_TWO_SIGNED(coor, WARPEDMODEL_PREC_BITS - 2) * 2;
}

IntMv gm_get_motion_vector(const GlobalMotionParams *gm, int allow_hp, BlockSize bsize, int mi_col,
                           int mi_row, int is_integer) {
    IntMv res;

    if (gm->gm_type == IDENTITY) {
        res.as_int = 0;
        return res;
    }

    const int32_t *mat = gm->gm_params;
    int            x = 0, y = 0, tx = 0, ty = 0;

    if (gm->gm_type == TRANSLATION) {
        res.as_mv.row = gm->gm_params[0] >> GM_TRANS_ONLY_PREC_DIFF;
        res.as_mv.col = gm->gm_params[1] >> GM_TRANS_ONLY_PREC_DIFF;
        assert(IMPLIES(1 & (res.as_mv.row | res.as_mv.col), allow_hp));
        if (is_integer) integer_mv_precision(&res.as_mv);
        return res;
    }

    x = block_center_x(mi_col, bsize);
    y = block_center_y(mi_row, bsize);

    if (gm->gm_type == ROTZOOM) {
        assert(gm->gm_params[5] == gm->gm_params[2]);
        assert(gm->gm_params[4] == -gm->gm_params[3]);
    }

    const int xc = (mat[2] - (1 << WARPEDMODEL_PREC_BITS)) * x + mat[3] * y + mat[0];
    const int yc = mat[4] * x + (mat[5] - (1 << WARPEDMODEL_PREC_BITS)) * y + mat[1];
    tx           = convert_to_trans_prec(allow_hp, xc);
    ty           = convert_to_trans_prec(allow_hp, yc);

    res.as_mv.row = ty;
    res.as_mv.col = tx;

    if (is_integer) integer_mv_precision(&res.as_mv);
    return res;
}

static INLINE int has_uni_comp_refs(const BlockModeInfo *mbmi) {
    return has_second_ref(mbmi) &&
           (!((mbmi->ref_frame[0] >= BWDREF_FRAME) ^ (mbmi->ref_frame[1] >= BWDREF_FRAME)));
}

int get_comp_reference_type_context(const PartitionInfo *xd) {
    int                        pred_context;
    const BlockModeInfo *const above_mbmi     = xd->above_mbmi;
    const BlockModeInfo *const left_mbmi      = xd->left_mbmi;
    const int                  above_in_image = xd->up_available;
    const int                  left_in_image  = xd->left_available;

    if (above_in_image && left_in_image) { // both edges available
        const int above_intra = !is_inter_block(above_mbmi);
        const int left_intra  = !is_inter_block(left_mbmi);

        if (above_intra && left_intra) { // intra/intra
            pred_context = 2;
        } else if (above_intra || left_intra) { // intra/inter
            const BlockModeInfo *inter_mbmi = above_intra ? left_mbmi : above_mbmi;

            if (!has_second_ref(inter_mbmi)) // single pred
                pred_context = 2;
            else // comp pred
                pred_context = 1 + 2 * has_uni_comp_refs(inter_mbmi);
        } else { // inter/inter
            const int              a_sg = !has_second_ref(above_mbmi);
            const int              l_sg = !has_second_ref(left_mbmi);
            const MvReferenceFrame frfa = above_mbmi->ref_frame[0];
            const MvReferenceFrame frfl = left_mbmi->ref_frame[0];

            if (a_sg && l_sg) { // single/single
                pred_context =
                    1 + 2 * (!(IS_BACKWARD_REF_FRAME(frfa) ^ IS_BACKWARD_REF_FRAME(frfl)));
            } else if (l_sg || a_sg) { // single/comp
                const int uni_rfc =
                    a_sg ? has_uni_comp_refs(left_mbmi) : has_uni_comp_refs(above_mbmi);

                if (!uni_rfc) // comp bidir
                    pred_context = 1;
                else // comp unidir
                    pred_context =
                        3 + (!(IS_BACKWARD_REF_FRAME(frfa) ^ IS_BACKWARD_REF_FRAME(frfl)));
            } else { // comp/comp
                const int a_uni_rfc = has_uni_comp_refs(above_mbmi);
                const int l_uni_rfc = has_uni_comp_refs(left_mbmi);

                if (!a_uni_rfc && !l_uni_rfc) // bidir/bidir
                    pred_context = 0;
                else if (!a_uni_rfc || !l_uni_rfc) // unidir/bidir
                    pred_context = 2;
                else // unidir/unidir
                    pred_context = 3 + (!((frfa == BWDREF_FRAME) ^ (frfl == BWDREF_FRAME)));
            }
        }
    } else if (above_in_image || left_in_image) { // one edge available
        const BlockModeInfo *edge_mbmi = above_in_image ? above_mbmi : left_mbmi;

        if (!is_inter_block(edge_mbmi)) { // intra
            pred_context = 2;
        } else { // inter
            if (!has_second_ref(edge_mbmi)) // single pred
                pred_context = 2;
            else // comp pred
                pred_context = 4 * has_uni_comp_refs(edge_mbmi);
        }
    } else { // no edges available
        pred_context = 2;
    }

    assert(pred_context >= 0 && pred_context < COMP_REF_TYPE_CONTEXTS);
    return pred_context;
}

int seg_feature_active(SegmentationParams *seg, int segment_id, SEG_LVL_FEATURES feature_id) {
    return seg->segmentation_enabled && seg->feature_enabled[segment_id][feature_id];
}
