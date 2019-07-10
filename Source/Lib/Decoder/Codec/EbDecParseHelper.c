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
#include "EbUtility.h"
#include "EbPictureBufferDesc.h"

#include "EbSvtAv1Dec.h"
#include "EbDecHandle.h"
#include "EbObuParse.h"

#include "EbDecParseHelper.h"
#include "EbDecUtils.h"
#include "EbTransforms.h"

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
    }
    else {
        if (diff <= 2 * (max - ref - 1)) {
            if (diff & 1)
                return ref + ((diff + 1) >> 1);
            else
                return ref - (diff >> 1);
        }
        return max - (diff + 1);
    }
}

void set_segment_id(EbDecHandle *dec_handle, int mi_offset,
    int x_mis, int y_mis, int segment_id)
{
    assert(segment_id >= 0 && segment_id < MAX_SEGMENTS);
    FrameHeader *frm_header = &dec_handle->frame_header;
    ParseCtxt *parse_ctxt = (ParseCtxt *)dec_handle->pv_parse_ctxt;

    for (int y = 0; y < y_mis; y++)
        for (int x = 0; x < x_mis; x++)
            parse_ctxt->parse_nbr4x4_ctxt.segment_maps[mi_offset +
            y * frm_header->mi_cols + x] = segment_id;
}

int bsize_to_max_depth(BlockSize bsize) {
    TxSize tx_size = max_txsize_rect_lookup[bsize];
    int depth = 0;
    while (depth < MAX_TX_DEPTH && tx_size != TX_4X4) {
        depth++;
        tx_size = sub_tx_size_map[tx_size];
    }
    return depth;
}

static INLINE int bsize_to_tx_size_cat(BlockSize bsize) {
    TxSize tx_size = max_txsize_rect_lookup[bsize];
    assert(tx_size != TX_4X4);
    int depth = 0;
    while (tx_size != TX_4X4) {
    depth++;
    tx_size = sub_tx_size_map[tx_size];
    assert(depth < 10);
    }
    assert(depth <= MAX_TX_CATS);
    return depth - 1;
}

int get_tx_size_context(const PartitionInfo_t *xd, ParseCtxt *parse_ctxt) {
    const ModeInfo_t *mbmi = xd->mi;
    const ModeInfo_t *const above_mbmi = xd->above_mbmi;
    const ModeInfo_t *const left_mbmi = xd->left_mbmi;
    const TxSize max_tx_size = max_txsize_rect_lookup[mbmi->sb_type];
    const uint8_t max_tx_wide = tx_size_wide[max_tx_size];
    const uint8_t max_tx_high = tx_size_high[max_tx_size];
    const int has_above = xd->up_available;
    const int has_left = xd->left_available;

    int above = parse_ctxt->parse_nbr4x4_ctxt.above_tx_wd[xd->mi_col] >= max_tx_wide;
    int left = parse_ctxt->parse_nbr4x4_ctxt.
        left_tx_ht[xd->mi_row - parse_ctxt->sb_row_mi] >= max_tx_high;

    if (has_above)
        if (dec_is_inter_block(above_mbmi))
            above = block_size_wide[above_mbmi->sb_type] >= max_tx_wide;

    if (has_left)
        if (dec_is_inter_block(left_mbmi))
            left = block_size_high[left_mbmi->sb_type] >= max_tx_high;

    if (has_above && has_left)
        return (above + left);
    else if (has_above)
        return above;
    else if (has_left)
        return left;
    else
        return 0;
}

TxSize depth_to_tx_size(int depth, BlockSize bsize) {
    TxSize tx_size = max_txsize_rect_lookup[bsize];
    for (int d = 0; d < depth; ++d) tx_size = sub_tx_size_map[tx_size];
    return tx_size;
}

void update_tx_context(ParseCtxt *parse_ctxt, PartitionInfo_t *pi,
    BlockSize bsize, TxSize txSize,
    int blk_row, int blk_col)
{
    int mi_row = pi->mi_row;
    int mi_col = pi->mi_col;
    ParseNbr4x4Ctxt *ngr_ctx = &parse_ctxt->parse_nbr4x4_ctxt;
    BlockSize b_size = bsize;
    if(dec_is_inter_block(pi->mi))
        b_size = txsize_to_bsize[txSize];
    uint8_t *const above_ctx = ngr_ctx->above_tx_wd + mi_col + blk_col;
    uint8_t *const left_ctx = ngr_ctx->left_tx_ht +
        (mi_row - parse_ctxt->sb_row_mi + blk_row);

    const uint8_t tx_wide = tx_size_wide[txSize];
    const uint8_t tx_high = tx_size_high[txSize];

    const int bw = mi_size_wide[b_size];
    const int bh = mi_size_high[b_size];
    memset(above_ctx, tx_wide, bw);
    memset(left_ctx, tx_high, bh);
}

TxSize read_selected_tx_size(PartitionInfo_t *xd, SvtReader *r,
    EbDecHandle *dec_handle)
{
    ParseCtxt *parse_ctxt = (ParseCtxt *)dec_handle->pv_parse_ctxt;
    const BlockSize bsize = xd->mi->sb_type;
    const int32_t tx_size_cat = bsize_to_tx_size_cat(bsize);
    const int maxTxDepth = bsize_to_max_depth(bsize);
    const int ctx = get_tx_size_context(xd, parse_ctxt);
    /*TODO : Change idx */
    const int depth = svt_read_symbol(r, parse_ctxt->cur_tile_ctx
        .tx_size_cdf[tx_size_cat][ctx],
        maxTxDepth + 1, ACCT_STR);
    assert(depth >= 0 && depth <= maxTxDepth);
    const TxSize tx_size = depth_to_tx_size(depth, bsize);
    return tx_size;
}

int is_intrabc_block(const ModeInfo_t *mbmi) {
    return mbmi->use_intrabc;
}

/* TODO : Harmonize with is_inter_block*/
int dec_is_inter_block(const ModeInfo_t *mbmi) {
    return is_intrabc_block(mbmi) || mbmi->ref_frame[0] > INTRA_FRAME;
}

int max_block_wide(PartitionInfo_t *part_info, int plane_bsize, int subx) {
    int max_blocks_wide = block_size_wide[plane_bsize];
    if (part_info->mb_to_right_edge < 0)
        max_blocks_wide += part_info->mb_to_right_edge >> (3 + subx);
    //Scale width in the transform block unit.
    return max_blocks_wide >> tx_size_wide_log2[0];
}

int max_block_high(PartitionInfo_t *part_info, int plane_bsize, int suby) {
    int max_blocks_high = block_size_high[plane_bsize];
    if (part_info->mb_to_bottom_edge < 0)
        max_blocks_high += part_info->mb_to_bottom_edge >> (3 + suby);
    // Scale the height in the transform block unit.
    return max_blocks_high >> tx_size_high_log2[0];
}

TxSize get_sqr_tx_size(int tx_dim) {
    switch (tx_dim) {
    case 128:
    case 64: return TX_64X64; break;
    case 32: return TX_32X32; break;
    case 16: return TX_16X16; break;
    case 8: return TX_8X8; break;
    default: return TX_4X4;
    }
}

int txfm_partition_context(TXFM_CONTEXT *above_ctx,
    TXFM_CONTEXT *left_ctx, BlockSize bsize, TxSize tx_size)
{
    const uint8_t txw = tx_size_wide[tx_size];
    const uint8_t txh = tx_size_high[tx_size];
    const int above = *above_ctx < txw;
    const int left = *left_ctx < txh;
    int category = TXFM_PARTITION_CONTEXTS;

    // dummy return, not used by others.
    if (tx_size <= TX_4X4) return 0;

    TxSize max_tx_size =
        get_sqr_tx_size(AOMMAX(block_size_wide[bsize], block_size_high[bsize]));

    if (max_tx_size >= TX_8X8) {
        category =
            (txsize_sqr_up_map[tx_size] != max_tx_size && max_tx_size > TX_8X8) +
            (TX_SIZES - 1 - max_tx_size) * 2;
    }
    assert(category != TXFM_PARTITION_CONTEXTS);
    return category * 3 + above + left;
}

int get_segdata(SegmentationParams *seg, int segment_id,
    SEG_LVL_FEATURES feature_id)
{
    return seg->feature_data[segment_id][feature_id];
}

int get_intra_inter_context(PartitionInfo_t *xd) {
    const ModeInfo_t *const above_mbmi = xd->above_mbmi;
    const ModeInfo_t *const left_mbmi = xd->left_mbmi;
    const int has_above = xd->up_available;
    const int has_left = xd->left_available;

    if (has_above && has_left) {  // both edges available
        const int above_intra = !dec_is_inter_block(above_mbmi);
        const int left_intra = !dec_is_inter_block(left_mbmi);
        return left_intra && above_intra ? 3 : left_intra || above_intra;
    }
    else if (has_above || has_left) {  // one edge available
        return 2 * !dec_is_inter_block(has_above ? above_mbmi : left_mbmi);
    }
    else
        return 0;
}

int use_angle_delta(BlockSize bsize) {
    return bsize >= BLOCK_8X8;
}

PredictionMode read_intra_mode(SvtReader *r, AomCdfProb *cdf) {
    return (PredictionMode)svt_read_symbol(r, cdf, INTRA_MODES, ACCT_STR);
}

/* TODO : Should reuse encode function */
int dec_is_chroma_reference(int mi_row, int mi_col, BlockSize bsize,
    int subsampling_x, int subsampling_y)
{
    const int bw = mi_size_wide[bsize];
    const int bh = mi_size_high[bsize];
    int ref_pos = ((mi_row & 0x01) || !(bh & 0x01) || !subsampling_y) &&
        ((mi_col & 0x01) || !(bw & 0x01) || !subsampling_x);
    return ref_pos;
}

UvPredictionMode read_intra_mode_uv(FRAME_CONTEXT *ec_ctx, SvtReader *r,
    CflAllowedType cfl_allowed, PredictionMode y_mode)
{
    const UvPredictionMode uv_mode =
        svt_read_symbol(r, ec_ctx->uv_mode_cdf[cfl_allowed][y_mode],
            UV_INTRA_MODES - !cfl_allowed, ACCT_STR);
    return uv_mode;
}

CflAllowedType is_cfl_allowed(PartitionInfo_t *xd,
    EbColorConfig* color_cfg, uint8_t *lossless_array)
{
    const ModeInfo_t *mbmi = xd->mi;
    const BlockSize bsize = mbmi->sb_type;
    assert(bsize < BlockSizeS_ALL);
    if (lossless_array[mbmi->segment_id]) {
        // In lossless, CfL is available when the partition size is equal to the
        // transform size.
        const int ssx = color_cfg->subsampling_x;
        const int ssy = color_cfg->subsampling_y;
        const int plane_bsize = get_plane_block_size(bsize, ssx, ssy);
        return (CflAllowedType)(plane_bsize == BLOCK_4X4);
    }
    // Spec: CfL is available to luma partitions lesser than or equal to 32x32
    return (CflAllowedType)(block_size_wide[bsize] <= 32 &&
        block_size_high[bsize] <= 32);
}

int allow_palette(int allow_screen_content_tools, BlockSize sb_type) {
    return allow_screen_content_tools && block_size_wide[sb_type] <= 64 &&
        block_size_high[sb_type] <= 64 && sb_type >= BLOCK_8X8;
}

int filter_intra_allowed_bsize(EbDecHandle *dec_handle, BlockSize bs) {
    if (!dec_handle->seq_header.enable_filter_intra || bs == BLOCK_INVALID)
        return 0;

    return block_size_wide[bs] <= 32 && block_size_high[bs] <= 32;
}

int filter_intra_allowed(EbDecHandle *dec_handle, const ModeInfo_t *mbmi) {
    return mbmi->mode == DC_PRED &&
        /* TO DO : Add when palette support comes */
        /*mbmi->palette_mode_info.palette_size[0] == 0 &&*/
        filter_intra_allowed_bsize(dec_handle, mbmi->sb_type);
}

int allow_intrabc(const EbDecHandle *dec_handle) {
    return  (dec_handle->frame_header.frame_type == KEY_FRAME
            || dec_handle->frame_header.frame_type == INTRA_ONLY_FRAME)
            && dec_handle->seq_header.seq_force_screen_content_tools
            && dec_handle->frame_header.allow_intrabc;
}

/*TODO: Move to common after segregating from encoder */
PredictionMode dec_get_uv_mode(UvPredictionMode mode) {
    assert(mode < UV_INTRA_MODES);
    static const PredictionMode uv2y[] = {
      DC_PRED,        // UV_DC_PRED
      V_PRED,         // UV_V_PRED
      H_PRED,         // UV_H_PRED
      D45_PRED,       // UV_D45_PRED
      D135_PRED,      // UV_D135_PRED
      D113_PRED,      // UV_D113_PRED
      D157_PRED,      // UV_D157_PRED
      D203_PRED,      // UV_D203_PRED
      D67_PRED,       // UV_D67_PRED
      SMOOTH_PRED,    // UV_SMOOTH_PRED
      SMOOTH_V_PRED,  // UV_SMOOTH_V_PRED
      SMOOTH_H_PRED,  // UV_SMOOTH_H_PRED
      PAETH_PRED,     // UV_PAETH_PRED
      DC_PRED,        // UV_CFL_PRED
      INTRA_INVALID,  // UV_INTRA_MODES
      INTRA_INVALID,  // UV_MODE_INVALID
    };
    return uv2y[mode];
}

TxType intra_mode_to_tx_type(const ModeInfo_t *mbmi, PlaneType plane_type) {
    static const TxType _intra_mode_to_tx_type[INTRA_MODES] = {
        DCT_DCT,    // DC
        ADST_DCT,   // V
        DCT_ADST,   // H
        DCT_DCT,    // D45
        ADST_ADST,  // D135
        ADST_DCT,   // D117
        DCT_ADST,   // D153
        DCT_ADST,   // D207
        ADST_DCT,   // D63
        ADST_ADST,  // SMOOTH
        ADST_DCT,   // SMOOTH_V
        DCT_ADST,   // SMOOTH_H
        ADST_ADST,  // PAETH
    };
    const PredictionMode mode =
        (plane_type == PLANE_TYPE_Y) ? mbmi->mode : dec_get_uv_mode(mbmi->uv_mode);
    assert(mode < INTRA_MODES);
    return _intra_mode_to_tx_type[mode];
}

int has_second_ref(const ModeInfo_t *mbmi) {
    return mbmi->ref_frame[1] > INTRA_FRAME;
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

IntMv_dec gm_get_motion_vector(const GlobalMotionParams *gm, int allow_hp,
    BlockSize bsize, int mi_col, int mi_row, int is_integer)
{
    IntMv_dec res;

    if (gm->gm_type == IDENTITY) {
        res.as_int = 0;
        return res;
    }

    const int32_t *mat = gm->gm_params;
    int x = 0, y = 0, tx = 0, ty = 0;

    if (gm->gm_type == TRANSLATION) {
        res.as_mv.row = gm->gm_params[0] >> GM_TRANS_ONLY_PREC_DIFF;
        res.as_mv.col = gm->gm_params[1] >> GM_TRANS_ONLY_PREC_DIFF;
        assert(IMPLIES(1 & (res.as_mv.row | res.as_mv.col), allow_hp));
        if (is_integer)
            integer_mv_precision(&res.as_mv);
        return res;
    }

    x = block_center_x(mi_col, bsize);
    y = block_center_y(mi_row, bsize);

    if (gm->gm_type == ROTZOOM) {
        assert(gm->gm_params[5] == gm->gm_params[2]);
        assert(gm->gm_params[4] == -gm->gm_params[3]);
    }

    const int xc =
        (mat[2] - (1 << WARPEDMODEL_PREC_BITS)) * x + mat[3] * y + mat[0];
    const int yc =
        mat[4] * x + (mat[5] - (1 << WARPEDMODEL_PREC_BITS)) * y + mat[1];
    tx = convert_to_trans_prec(allow_hp, xc);
    ty = convert_to_trans_prec(allow_hp, yc);

    res.as_mv.row = ty;
    res.as_mv.col = tx;

    if (is_integer)
        integer_mv_precision(&res.as_mv);
    return res;
}

int get_txb_wide(TxSize tx_size) {
    tx_size = av1_get_adjusted_tx_size(tx_size);
    return tx_size_wide[tx_size];
}

int get_txb_high(TxSize tx_size) {
    tx_size = av1_get_adjusted_tx_size(tx_size);
    return tx_size_high[tx_size];
}

int get_lower_levels_ctx_eob(int bwl, int height, int scan_idx) {
    if (scan_idx == 0) return 0;
    if (scan_idx <= (height << bwl) / 8) return 1;
    if (scan_idx <= (height << bwl) / 4) return 2;
    return 3;
}

uint8_t *set_levels(uint8_t *const levels_buf, const int width) {
    return levels_buf + TX_PAD_TOP * (width + TX_PAD_HOR);
}

int get_padded_idx(const int idx, const int bwl) {
    return idx + ((idx >> bwl) << TX_PAD_HOR_LOG2);
}

int get_txb_bwl(TxSize tx_size) {
    tx_size = av1_get_adjusted_tx_size(tx_size);
    return tx_size_wide_log2[tx_size];
}

static INLINE int has_uni_comp_refs(const ModeInfo_t *mbmi) {
  return has_second_ref(mbmi) && (!((mbmi->ref_frame[0] >= BWDREF_FRAME) ^
                                    (mbmi->ref_frame[1] >= BWDREF_FRAME)));
}

int get_comp_reference_type_context(const PartitionInfo_t *xd) {
#define CHECK_BACKWARD_REFS(ref_frame) \
  (((ref_frame) >= BWDREF_FRAME) && ((ref_frame) <= ALTREF_FRAME))
#define IS_BACKWARD_REF_FRAME(ref_frame) CHECK_BACKWARD_REFS(ref_frame)

    int pred_context;
    const ModeInfo_t *const above_mbmi = xd->above_mbmi;
    const ModeInfo_t *const left_mbmi = xd->left_mbmi;
    const int above_in_image = xd->up_available;
    const int left_in_image = xd->left_available;

    if (above_in_image && left_in_image) {  // both edges available
        const int above_intra = !dec_is_inter_block(above_mbmi);
        const int left_intra = !dec_is_inter_block(left_mbmi);

        if (above_intra && left_intra) {  // intra/intra
            pred_context = 2;
        }
        else if (above_intra || left_intra) {  // intra/inter
            const ModeInfo_t *inter_mbmi = above_intra ? left_mbmi : above_mbmi;

            if (!has_second_ref(inter_mbmi))  // single pred
                pred_context = 2;
            else  // comp pred
                pred_context = 1 + 2 * has_uni_comp_refs(inter_mbmi);
        }
        else {  // inter/inter
            const int a_sg = !has_second_ref(above_mbmi);
            const int l_sg = !has_second_ref(left_mbmi);
            const MvReferenceFrame frfa = above_mbmi->ref_frame[0];
            const MvReferenceFrame frfl = left_mbmi->ref_frame[0];

            if (a_sg && l_sg) {  // single/single
                pred_context = 1 + 2 * (!(IS_BACKWARD_REF_FRAME(frfa) ^
                    IS_BACKWARD_REF_FRAME(frfl)));
            }
            else if (l_sg || a_sg) {  // single/comp
                const int uni_rfc =
                    a_sg ? has_uni_comp_refs(left_mbmi) : has_uni_comp_refs(above_mbmi);

                if (!uni_rfc)  // comp bidir
                    pred_context = 1;
                else  // comp unidir
                    pred_context = 3 + (!(IS_BACKWARD_REF_FRAME(frfa) ^
                        IS_BACKWARD_REF_FRAME(frfl)));
            }
            else {  // comp/comp
                const int a_uni_rfc = has_uni_comp_refs(above_mbmi);
                const int l_uni_rfc = has_uni_comp_refs(left_mbmi);

                if (!a_uni_rfc && !l_uni_rfc)  // bidir/bidir
                    pred_context = 0;
                else if (!a_uni_rfc || !l_uni_rfc)  // unidir/bidir
                    pred_context = 2;
                else  // unidir/unidir
                    pred_context =
                    3 + (!((frfa == BWDREF_FRAME) ^ (frfl == BWDREF_FRAME)));
            }
        }
    }
    else if (above_in_image || left_in_image) {  // one edge available
        const ModeInfo_t *edge_mbmi = above_in_image ? above_mbmi : left_mbmi;

        if (!dec_is_inter_block(edge_mbmi)) {  // intra
            pred_context = 2;
        }
        else {                           // inter
            if (!has_second_ref(edge_mbmi))  // single pred
                pred_context = 2;
            else  // comp pred
                pred_context = 4 * has_uni_comp_refs(edge_mbmi);
        }
    }
    else {  // no edges available
        pred_context = 2;
    }

    assert(pred_context >= 0 && pred_context < COMP_REF_TYPE_CONTEXTS);
    return pred_context;
}

AomCdfProb *get_y_mode_cdf(FRAME_CONTEXT *tile_ctx,
    const ModeInfo_t *above_mi, const ModeInfo_t *left_mi)
{
    const PredictionMode above = above_mi ? above_mi->mode : DC_PRED;
    const PredictionMode left = left_mi ? left_mi->mode : DC_PRED;
    const int above_ctx = intra_mode_context[above];
    const int left_ctx = intra_mode_context[left];
    return tile_ctx->kf_y_cdf[above_ctx][left_ctx];
}

int is_interintra_allowed_bsize(const BlockSize bsize) {
    return (bsize >= BLOCK_8X8) && (bsize <= BLOCK_32X32);
}

int is_interintra_allowed_mode(const PredictionMode mode) {
    return (mode >= SINGLE_INTER_MODE_START) && (mode < SINGLE_INTER_MODE_END);
}

int is_interintra_allowed_ref(const MvReferenceFrame rf[2]) {
    return (rf[0] > INTRA_FRAME) && (rf[1] <= INTRA_FRAME);
}

int is_interintra_allowed(const ModeInfo_t *mbmi) {
    return is_interintra_allowed_bsize(mbmi->sb_type) &&
        is_interintra_allowed_mode(mbmi->mode) &&
        is_interintra_allowed_ref(mbmi->ref_frame);
}

MotionMode dec_motion_mode_allowed() {
    return SIMPLE_TRANSLATION;
}

int seg_feature_active(SegmentationParams *seg, int segment_id,
    SEG_LVL_FEATURES feature_id)
{
    return seg->segmentation_enabled && seg->feature_data[segment_id][feature_id];
}
