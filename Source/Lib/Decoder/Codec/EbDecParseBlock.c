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
#include "EbPictureBufferDesc.h"

#include "EbDecHandle.h"
#include "EbDecBitReader.h"
#include "EbObuParse.h"

#include "EbDecParseHelper.h"
#include "EbInvTransforms.h"

#include "EbDecNbr.h"
#include "EbDecPicMgr.h"
#include "EbDecUtils.h"

#include "EbDecParseInterBlock.h"
#include "EbDecProcessFrame.h"

#include "EbCommonUtils.h"
#include "EbCoefficients.h"
#include "EbIntraPrediction.h"
#include "EbCabacContextModel.h"

#include "EbLog.h"

#if ENABLE_ENTROPY_TRACE
FILE *temp_fp;
int   enable_dump;
#endif

#define READ_REF_BIT(pname) svt_read_symbol(r, get_pred_cdf_##pname(pi), 2, ACCT_STR)
#define SQR_BLOCK_SIZES 6

typedef AomCdfProb (*BaseCdfArr)[CDF_SIZE(4)];
typedef AomCdfProb (*BrCdfArr)[CDF_SIZE(BR_CDF_SIZE)];

typedef struct txb_ctx {
    int txb_skip_ctx;
    int dc_sign_ctx;
} TXB_CTX;

static INLINE int get_palette_mode_ctx(PartitionInfo *pi) {
    BlockModeInfo *above_mi = pi->above_mbmi;
    BlockModeInfo *left_mi  = pi->left_mbmi;
    int            ctx      = 0;
    if (above_mi) ctx += (above_mi->palette_size[0] > 0);
    if (left_mi) ctx += (left_mi->palette_size[0] > 0);
    return ctx;
}

static INLINE void palette_add_to_cache(uint16_t *cache, int *n, uint16_t val) {
    // Do not add an already existing value
    if (*n > 0 && val == cache[*n - 1]) return;

    cache[(*n)++] = val;
}

int av1_get_palette_cache(ParseCtxt *parse_ctx, PartitionInfo *pi, int plane, uint16_t *cache) {
    const int row = -pi->mb_to_top_edge >> 3;
    // Do not refer to above SB row when on SB boundary.
    BlockModeInfo *       above_mi  = (row % (1 << MIN_SB_SIZE_LOG2)) ? pi->above_mbmi : NULL;
    ParseAboveNbr4x4Ctxt *above_ctx = parse_ctx->parse_above_nbr4x4_ctxt;
    ParseLeftNbr4x4Ctxt * left_ctx  = parse_ctx->parse_left_nbr4x4_ctxt;
    BlockModeInfo *       left_mi   = pi->left_mbmi;
    uint8_t               above_n = 0, left_n = 0;
    if (above_mi) above_n = above_mi->palette_size[plane != 0];
    if (left_mi) left_n = left_mi->palette_size[plane != 0];
    if (above_n == 0 && left_n == 0) return 0;
    int above_idx = 0;
    int left_idx  = 0;
    int n         = 0;

    const uint16_t *above_colors =
        above_mi ? above_ctx->above_palette_colors[plane] +
                       (PALETTE_MAX_SIZE * ((pi->mi_col - parse_ctx->sb_col_mi) % 16))
                 : NULL;
    const uint16_t *left_colors = left_mi
                                      ? left_ctx->left_palette_colors[plane] +
                                            (PALETTE_MAX_SIZE * (pi->mi_row - parse_ctx->sb_row_mi))
                                      : NULL;
    // Merge the sorted lists of base colors from above and left to get
    // combined sorted color cache.
    while (above_n > 0 && left_n > 0) {
        uint16_t v_above = above_colors[above_idx];
        uint16_t v_left  = left_colors[left_idx];
        if (v_left < v_above) {
            palette_add_to_cache(cache, &n, v_left);
            ++left_idx, --left_n;
        } else {
            palette_add_to_cache(cache, &n, v_above);
            ++above_idx, --above_n;
            if (v_left == v_above) ++left_idx, --left_n;
        }
    }
    while (above_n-- > 0) {
        uint16_t val = above_colors[above_idx++];
        palette_add_to_cache(cache, &n, val);
    }
    while (left_n-- > 0) {
        uint16_t val = left_colors[left_idx++];
        palette_add_to_cache(cache, &n, val);
    }
    assert(n <= 2 * PALETTE_MAX_SIZE);
    return n;
}

// Merge the sorted list of cached colors(cached_colors[0...n_cached_colors-1])
// and the sorted list of transmitted colors(colors[n_cached_colors...n-1]) into
// one single sorted list(colors[...]).
static void merge_colors(uint16_t *colors, uint16_t *cached_colors, int n_colors,
                         int n_cached_colors) {
    if (n_cached_colors == 0) return;
    int cache_idx = 0, trans_idx = n_cached_colors;
    for (int i = 0; i < n_colors; ++i) {
        if (cache_idx < n_cached_colors &&
            (trans_idx >= n_colors || cached_colors[cache_idx] <= colors[trans_idx])) {
            colors[i] = cached_colors[cache_idx++];
        } else {
            assert(trans_idx < n_colors);
            colors[i] = colors[trans_idx++];
        }
    }
}

static void read_palette_colors_y(ParseCtxt *parse_ctx, PartitionInfo *pi, int bit_depth) {
    SvtReader *r = &parse_ctx->r;
    uint16_t   color_cache[2 * PALETTE_MAX_SIZE];
    uint16_t   cached_colors[PALETTE_MAX_SIZE];
    const int  n_cache = av1_get_palette_cache(parse_ctx, pi, 0, color_cache);
    const int  n       = pi->mi->palette_size[AOM_PLANE_Y];
    int        idx     = 0;
    for (int i = 0; i < n_cache && idx < n; ++i)
        if (svt_read_bit(r, ACCT_STR)) cached_colors[idx++] = color_cache[i];
    if (idx < n) {
        const int n_cached_colors                     = idx;
        parse_ctx->palette_colors[AOM_PLANE_Y][idx++] = svt_read_literal(r, bit_depth, ACCT_STR);
        if (idx < n) {
            const int min_bits = bit_depth - 3;
            int       bits     = min_bits + svt_read_literal(r, 2, ACCT_STR);
            int range = (1 << bit_depth) - parse_ctx->palette_colors[AOM_PLANE_Y][idx - 1] - 1;
            for (; idx < n; ++idx) {
                assert(range >= 0);
                const int delta = svt_read_literal(r, bits, ACCT_STR) + 1;
                parse_ctx->palette_colors[AOM_PLANE_Y][idx] =
                    clamp(parse_ctx->palette_colors[AOM_PLANE_Y][idx - 1] + delta,
                          0,
                          (1 << bit_depth) - 1);
                range -= (parse_ctx->palette_colors[AOM_PLANE_Y][idx] -
                          parse_ctx->palette_colors[AOM_PLANE_Y][idx - 1]);
                bits = AOMMIN(bits, av1_ceil_log2(range));
            }
        }
        merge_colors(parse_ctx->palette_colors[AOM_PLANE_Y], cached_colors, n, n_cached_colors);
    } else {
        eb_memcpy(parse_ctx->palette_colors[AOM_PLANE_Y], cached_colors, n * sizeof(cached_colors[0]));
    }
}

static void read_palette_colors_uv(ParseCtxt *parse_ctx, PartitionInfo *pi, int bit_depth) {
    SvtReader *r = &parse_ctx->r;
    const int  n = pi->mi->palette_size[AOM_PLANE_U];
    // U channel colors.
    uint16_t  color_cache[2 * PALETTE_MAX_SIZE];
    uint16_t  cached_colors[PALETTE_MAX_SIZE];
    const int n_cache = av1_get_palette_cache(parse_ctx, pi, 1, color_cache);
    int       idx     = 0;
    for (int i = 0; i < n_cache && idx < n; ++i)
        if (svt_read_bit(r, ACCT_STR)) cached_colors[idx++] = color_cache[i];
    if (idx < n) {
        const int n_cached_colors                     = idx;
        parse_ctx->palette_colors[AOM_PLANE_U][idx++] = svt_read_literal(r, bit_depth, ACCT_STR);
        if (idx < n) {
            const int min_bits = bit_depth - 3;
            int       bits     = min_bits + svt_read_literal(r, 2, ACCT_STR);
            int       range    = (1 << bit_depth) - parse_ctx->palette_colors[AOM_PLANE_U][idx - 1];
            for (; idx < n; ++idx) {
                assert(range >= 0);
                const int delta = svt_read_literal(r, bits, ACCT_STR);
                parse_ctx->palette_colors[AOM_PLANE_U][idx] =
                    clamp(parse_ctx->palette_colors[AOM_PLANE_U][idx - 1] + delta,
                          0,
                          (1 << bit_depth) - 1);
                range -= (parse_ctx->palette_colors[AOM_PLANE_U][idx] -
                          parse_ctx->palette_colors[AOM_PLANE_U][idx - 1]);
                bits = AOMMIN(bits, av1_ceil_log2(range));
            }
        }
        merge_colors(parse_ctx->palette_colors[AOM_PLANE_U], cached_colors, n, n_cached_colors);
    } else {
        eb_memcpy(parse_ctx->palette_colors[AOM_PLANE_U], cached_colors, n * sizeof(cached_colors[0]));
    }

    // V channel colors.
    if (svt_read_bit(r, ACCT_STR)) { // Delta encoding.
        const int min_bits_v                      = bit_depth - 4;
        const int max_val                         = 1 << bit_depth;
        int       bits                            = min_bits_v + svt_read_literal(r, 2, ACCT_STR);
        parse_ctx->palette_colors[AOM_PLANE_V][0] = svt_read_literal(r, bit_depth, ACCT_STR);
        for (int i = 1; i < n; ++i) {
            int delta = svt_read_literal(r, bits, ACCT_STR);
            if (delta && svt_read_bit(r, ACCT_STR)) delta = -delta;
            int val = (int)parse_ctx->palette_colors[AOM_PLANE_V][i - 1] + delta;
            if (val < 0) val += max_val;
            if (val >= max_val) val -= max_val;
            parse_ctx->palette_colors[AOM_PLANE_V][i] = val;
        }
    } else {
        for (int i = 0; i < n; ++i) {
            parse_ctx->palette_colors[AOM_PLANE_V][i] = svt_read_literal(r, bit_depth, ACCT_STR);
        }
    }
}

void palette_mode_info(ParseCtxt *parse_ctxt, PartitionInfo *pi) {
    SvtReader *          r          = &parse_ctxt->r;
    FRAME_CONTEXT *      frm_ctx    = &parse_ctxt->cur_tile_ctx;
    EbColorConfig *      color_info = &parse_ctxt->seq_header->color_config;
    const int            num_planes = color_info->mono_chrome ? 1 : MAX_MB_PLANE;
    BlockModeInfo *const mbmi       = pi->mi;
    const BlockSize      bsize      = mbmi->sb_type;
    assert(allow_palette(parse_ctxt->frame_header->allow_screen_content_tools, bsize));
    const int bsize_ctx = get_palette_bsize_ctx(bsize);

    if (mbmi->mode == DC_PRED) {
        const int palette_mode_ctx = get_palette_mode_ctx(pi);
        const int modev            = svt_read_symbol(
            r, frm_ctx->palette_y_mode_cdf[bsize_ctx][palette_mode_ctx], 2, ACCT_STR);
        if (modev) {
            mbmi->palette_size[0] =
                svt_read_symbol(
                    r, frm_ctx->palette_y_size_cdf[bsize_ctx], PALETTE_SIZES, ACCT_STR) +
                2;
            memset(parse_ctxt->palette_colors[AOM_PLANE_Y],
                   0,
                   mbmi->palette_size[0] * sizeof(uint16_t));
            read_palette_colors_y(parse_ctxt, pi, parse_ctxt->seq_header->color_config.bit_depth);
        }
    }

    if (num_planes > 1 && mbmi->uv_mode == UV_DC_PRED && pi->is_chroma_ref) {
        const int palette_uv_mode_ctx = (mbmi->palette_size[0] > 0);
        const int modev =
            svt_read_symbol(r, frm_ctx->palette_uv_mode_cdf[palette_uv_mode_ctx], 2, ACCT_STR);
        if (modev) {
            mbmi->palette_size[1] =
                svt_read_symbol(
                    r, frm_ctx->palette_uv_size_cdf[bsize_ctx], PALETTE_SIZES, ACCT_STR) +
                2;
            memset(parse_ctxt->palette_colors[AOM_PLANE_U],
                   0,
                   mbmi->palette_size[1] * sizeof(uint16_t));
            memset(parse_ctxt->palette_colors[AOM_PLANE_V],
                   0,
                   mbmi->palette_size[1] * sizeof(uint16_t));
            read_palette_colors_uv(parse_ctxt, pi, parse_ctxt->seq_header->color_config.bit_depth);
        }
    }
}

static INLINE int filter_intra_allowed_bsize(ParseCtxt *parse_ctxt, BlockSize bs) {
#if FILTER_INTRA_CLI
    if (!parse_ctxt->seq_header->filter_intra_level || bs == BLOCK_INVALID) return 0;
#else
    if (!parse_ctxt->seq_header->enable_filter_intra || bs == BLOCK_INVALID) return 0;
#endif

    return block_size_wide[bs] <= 32 && block_size_high[bs] <= 32;
}

static INLINE int filter_intra_allowed(ParseCtxt *parse_ctxt, const BlockModeInfo *mbmi) {
    return mbmi->mode == DC_PRED && mbmi->palette_size[0] == 0 &&
           filter_intra_allowed_bsize(parse_ctxt, mbmi->sb_type);
}

void filter_intra_mode_info(ParseCtxt *parse_ctxt, PartitionInfo *xd) {
    SvtReader *            r                      = &parse_ctxt->r;
    BlockModeInfo *const   mbmi                   = xd->mi;
    FilterIntraModeInfo_t *filter_intra_mode_info = &mbmi->filter_intra_mode_info;
    FRAME_CONTEXT *        frm_ctx                = &parse_ctxt->cur_tile_ctx;

    if (filter_intra_allowed(parse_ctxt, mbmi)) {
        filter_intra_mode_info->use_filter_intra =
            svt_read_symbol(r, frm_ctx->filter_intra_cdfs[mbmi->sb_type], 2, ACCT_STR);
        if (filter_intra_mode_info->use_filter_intra) {
            filter_intra_mode_info->filter_intra_mode =
                svt_read_symbol(r, frm_ctx->filter_intra_mode_cdf, FILTER_INTRA_MODES, ACCT_STR);
        }
    } else
        filter_intra_mode_info->use_filter_intra = 0;
}

uint8_t read_cfl_alphas(FRAME_CONTEXT *ec_ctx, SvtReader *r, uint8_t *signs_out) {
    const int joint_sign = svt_read_symbol(r, ec_ctx->cfl_sign_cdf, CFL_JOINT_SIGNS, "cfl:signs");
    uint8_t   idx        = 0;
    // Magnitudes are only coded for nonzero values
    if (CFL_SIGN_U(joint_sign) != CFL_SIGN_ZERO) {
        AomCdfProb *cdf_u = ec_ctx->cfl_alpha_cdf[CFL_CONTEXT_U(joint_sign)];
        idx = svt_read_symbol(r, cdf_u, CFL_ALPHABET_SIZE, "cfl:alpha_u") << CFL_ALPHABET_SIZE_LOG2;
    }
    if (CFL_SIGN_V(joint_sign) != CFL_SIGN_ZERO) {
        AomCdfProb *cdf_v = ec_ctx->cfl_alpha_cdf[CFL_CONTEXT_V(joint_sign)];
        idx += svt_read_symbol(r, cdf_v, CFL_ALPHABET_SIZE, "cfl:alpha_v");
    }
    *signs_out = joint_sign;
    return idx;
}

void read_cdef(ParseCtxt *parse_ctxt, PartitionInfo *xd) {
    SvtReader *          r    = &parse_ctxt->r;
    BlockModeInfo *const mbmi = xd->mi;
    if (mbmi->skip || parse_ctxt->frame_header->coded_lossless ||
        !parse_ctxt->seq_header->enable_cdef || parse_ctxt->frame_header->allow_intrabc) {
        return;
    }
    int       cdf_size = mi_size_wide[BLOCK_64X64];
    int       row      = xd->mi_row & cdf_size;
    int       col      = xd->mi_col & cdf_size;
    const int index = parse_ctxt->seq_header->sb_size == BLOCK_128X128 ? !!(col) + 2 * !!(row) : 0;
    int8_t *  cdef_strength = xd->cdef_strength;
    if (cdef_strength[index] == -1) {
        cdef_strength[index] =
            svt_read_literal(r, parse_ctxt->frame_header->cdef_params.cdef_bits, ACCT_STR);
        /* Populate to nearby 64x64s if needed based on h4 & w4 */
        if (parse_ctxt->seq_header->sb_size == BLOCK_128X128) {
            int w4 = mi_size_wide[mbmi->sb_type];
            int h4 = mi_size_high[mbmi->sb_type];
            for (int i = row; i < row + h4; i += cdf_size) {
                for (int j = col; j < col + w4; j += cdf_size) {
                    cdef_strength[!!(j & cdf_size) + 2 * !!(i & cdf_size)] = cdef_strength[index];
                }
            }
        }
    }
}

void read_delta_qindex(ParseCtxt *parse_ctxt, BlockModeInfo *const mbmi, int32_t *cur_qind,
                       int32_t *sb_delta_q) {
    SvtReader *   r = &parse_ctxt->r;
    int           sign, abs, reduced_delta_qindex = 0;
    BlockSize     bsize          = mbmi->sb_type;
    DeltaQParams *delta_q_params = &parse_ctxt->frame_header->delta_q_params;

    if ((bsize != parse_ctxt->seq_header->sb_size || mbmi->skip == 0)) {
        abs = svt_read_symbol(r, parse_ctxt->cur_tile_ctx.delta_q_cdf, DELTA_Q_PROBS + 1, ACCT_STR);

        if (abs == DELTA_Q_SMALL) {
            const int rem_bits = svt_read_literal(r, 3, ACCT_STR) + 1;
            const int thr      = (1 << rem_bits) + 1;
            abs                = svt_read_literal(r, rem_bits, ACCT_STR) + thr;
        }

        if (abs)
            sign = svt_read_bit(r, ACCT_STR);
        else
            sign = 1;
        reduced_delta_qindex = sign ? -abs : abs;
        reduced_delta_qindex =
            clamp(*cur_qind + (reduced_delta_qindex << delta_q_params->delta_q_res), 1, MAXQ);
        *sb_delta_q = *cur_qind = reduced_delta_qindex;
    }
}

int read_delta_lflevel(ParseCtxt *parse_ctxt, AomCdfProb *cdf, BlockModeInfo *mbmi,
                       int32_t delta_lf) {
    SvtReader *     r                     = &parse_ctxt->r;
    int             tmp_lvl               = 0;
    int             reduced_delta_lflevel = 0;
    const BlockSize bsize                 = mbmi->sb_type;

    if (bsize == parse_ctxt->seq_header->sb_size && mbmi->skip) return delta_lf;

    DeltaLfParams *delta_lf_params = &parse_ctxt->frame_header->delta_lf_params;

    int abs = svt_read_symbol(r, cdf, DELTA_LF_PROBS + 1, ACCT_STR);
    if (abs == DELTA_LF_SMALL) {
        const int rem_bits = svt_read_literal(r, 3, ACCT_STR) + 1;
        const int thr      = (1 << rem_bits) + 1;
        abs                = svt_read_literal(r, rem_bits, ACCT_STR) + thr;
    }
    const int sign        = abs ? svt_read_bit(r, ACCT_STR) : 1;
    reduced_delta_lflevel = sign ? -abs : abs;
    tmp_lvl = (clamp(delta_lf + (reduced_delta_lflevel << delta_lf_params->delta_lf_res),
                     -MAX_LOOP_FILTER,
                     MAX_LOOP_FILTER));
    return tmp_lvl;
}

int read_skip(ParseCtxt *parse_ctxt, PartitionInfo *xd, int segment_id) {
    SvtReader *r = &parse_ctxt->r;
    //uint8_t segIdPreSkip = dec_handle->frame_header.segmentation_params.seg_id_pre_skip;
    if (seg_feature_active(
            &parse_ctxt->frame_header->segmentation_params, segment_id, SEG_LVL_SKIP)) {
        return 1;
    } else {
        const int above_skip = xd->above_mbmi ? xd->above_mbmi->skip : 0;
        const int left_skip  = xd->left_mbmi ? xd->left_mbmi->skip : 0;
        int       ctx        = above_skip + left_skip;
        return svt_read_symbol(r, parse_ctxt->cur_tile_ctx.skip_cdfs[ctx], 2, ACCT_STR);
    }
}

int read_skip_mode(ParseCtxt *parse_ctxt, PartitionInfo *xd, int segment_id) {
    SvtReader *         r   = &parse_ctxt->r;
    SegmentationParams *seg = &parse_ctxt->frame_header->segmentation_params;
    assert(xd->mi->sb_type < BlockSizeS_ALL);
    if (seg_feature_active(seg, segment_id, SEG_LVL_SKIP) ||
        seg_feature_active(seg, segment_id, SEG_LVL_REF_FRAME) ||
        seg_feature_active(seg, segment_id, SEG_LVL_GLOBALMV) ||
        !parse_ctxt->frame_header->skip_mode_params.skip_mode_flag ||
        block_size_wide[xd->mi->sb_type] < 8 || block_size_high[xd->mi->sb_type] < 8) {
        return 0;
    }
    int above_skip_mode = xd->above_mbmi ? xd->above_mbmi->skip_mode : 0;
    int left_skip_mode  = xd->left_mbmi ? xd->left_mbmi->skip_mode : 0;
    int ctx             = above_skip_mode + left_skip_mode;
    return svt_read_symbol(r, parse_ctxt->cur_tile_ctx.skip_mode_cdfs[ctx], 2, ACCT_STR);
}

// If delta q is present, reads delta_q index.
// Also reads delta_q loop filter levels, if present.
static void read_delta_params(ParseCtxt *parse_ctxt, PartitionInfo *xd) {
    DeltaQParams *       delta_q_params  = &parse_ctxt->frame_header->delta_q_params;
    DeltaLfParams *      delta_lf_params = &parse_ctxt->frame_header->delta_lf_params;
    SBInfo *             sb_info         = xd->sb_info;
    BlockModeInfo *const mbmi            = &xd->mi[0];

    if (!parse_ctxt->read_deltas) return;

    if (delta_q_params->delta_q_present) {
        read_delta_qindex(parse_ctxt, mbmi, &parse_ctxt->cur_q_ind, &sb_info->sb_delta_q[0]);
    }

    FRAME_CONTEXT *const ec_ctx = &parse_ctxt->cur_tile_ctx;

    int frame_lf_count = 1;
    if (delta_lf_params->delta_lf_present) {
        if (delta_lf_params->delta_lf_multi) {
            EbColorConfig *color_info = &parse_ctxt->seq_header->color_config;
            int            num_planes = color_info->mono_chrome ? 1 : MAX_MB_PLANE;
            frame_lf_count            = num_planes > 1 ? FRAME_LF_COUNT : FRAME_LF_COUNT - 2;

            for (int lf_id = 0; lf_id < frame_lf_count; ++lf_id) {
                parse_ctxt->delta_lf[lf_id] = sb_info->sb_delta_lf[lf_id] =
                    read_delta_lflevel(parse_ctxt,
                                       ec_ctx->delta_lf_multi_cdf[lf_id],
                                       mbmi,
                                       parse_ctxt->delta_lf[lf_id]);
            }
        } else {
            parse_ctxt->delta_lf[0] = sb_info->sb_delta_lf[0] =
                read_delta_lflevel(parse_ctxt, ec_ctx->delta_lf_cdf, mbmi, parse_ctxt->delta_lf[0]);
        }
    }
}

int intra_angle_info(SvtReader *r, AomCdfProb *cdf, PredictionMode mode, BlockSize bsize) {
    int angle_delta_y = 0;
    if (av1_use_angle_delta(bsize, EB_TRUE) && av1_is_directional_mode(mode)) {
        const int sym = svt_read_symbol(r, cdf, 2 * MAX_ANGLE_DELTA + 1, ACCT_STR);
        angle_delta_y = sym - MAX_ANGLE_DELTA;
    }
    return angle_delta_y;
}

int get_segment_id(FrameHeader *frm_info, uint8_t *segment_ids, BlockSize bsize, uint32_t mi_row,
                   uint32_t mi_col) {
    const int      mi_offset = mi_row * frm_info->mi_cols + mi_col;
    const uint32_t bw        = mi_size_wide[bsize];
    const uint32_t bh        = mi_size_high[bsize];
    const int      xmis      = AOMMIN(frm_info->mi_cols - mi_col, bw);
    const int      ymis      = AOMMIN(frm_info->mi_rows - mi_row, bh);
    int            x, y, segment_id = MAX_SEGMENTS - 1;

    for (y = 0; y < ymis; ++y)
        for (x = 0; x < xmis; ++x)
            segment_id = AOMMIN(segment_id, segment_ids[mi_offset + y * frm_info->mi_cols + x]);

    assert(segment_id >= 0 && segment_id < MAX_SEGMENTS);
    return segment_id;
}

static int read_segment_id(EbDecHandle *dec_handle, ParseCtxt *parse_ctxt, PartitionInfo *xd,
                           int skip) {
    SvtReader *r       = &parse_ctxt->r;
    int        cdf_num = 0;

    int      prev_ul  = -1; // top left segment_id
    int      prev_l   = -1; // left segment_id
    int      prev_u   = -1; // top segment_id
    int      pred     = -1;
    uint32_t mi_row   = xd->mi_row;
    uint32_t mi_col   = xd->mi_col;
    uint8_t *seg_maps = dec_handle->cur_pic_buf[0]->segment_maps;

    if ((xd->up_available) && (xd->left_available)) {
        prev_ul =
            get_segment_id(parse_ctxt->frame_header, seg_maps, BLOCK_4X4, mi_row - 1, mi_col - 1);
    }
    if (xd->up_available) {
        prev_u = get_segment_id(parse_ctxt->frame_header, seg_maps, BLOCK_4X4, mi_row - 1, mi_col);
    }
    if (xd->left_available) {
        prev_l = get_segment_id(parse_ctxt->frame_header, seg_maps, BLOCK_4X4, mi_row, mi_col - 1);
    }

    // Pick CDF index based on number of matching/out-of-bounds segment IDs.
    if (prev_ul < 0) /* Edge cases */
        cdf_num = 0;
    else if ((prev_ul == prev_u) && (prev_ul == prev_l))
        cdf_num = 2;
    else if ((prev_ul == prev_u) || (prev_ul == prev_l) || (prev_u == prev_l))
        cdf_num = 1;

    // If 2 or more are identical returns that as predictor, otherwise prev_l.
    if (prev_u == -1) // edge case
        pred = prev_l == -1 ? 0 : prev_l;
    else if (prev_l == -1) // edge case
        pred = prev_u;
    else
        pred = (prev_ul == prev_u) ? prev_u : prev_l;

    if (skip) return pred;

    FRAME_CONTEXT *     ec_ctx = &parse_ctxt->cur_tile_ctx;
    SegmentationParams *seg    = &(parse_ctxt->frame_header->segmentation_params);

    struct segmentation_probs *segp     = &ec_ctx->seg;
    AomCdfProb *               pred_cdf = segp->spatial_pred_seg_cdf[cdf_num];

    int coded_id = svt_read_symbol(r, pred_cdf, MAX_SEGMENTS, ACCT_STR);
    return neg_deinterleave(coded_id, pred, seg->last_active_seg_id + 1);
}

int intra_segment_id(EbDecHandle *dec_handle, ParseCtxt *parse_ctxt, PartitionInfo *xd, int bsize,
                     int skip) {
    SegmentationParams *seg        = &parse_ctxt->frame_header->segmentation_params;
    int                 segment_id = 0;
    int                 mi_row     = xd->mi_row;
    int                 mi_col     = xd->mi_col;
    if (seg->segmentation_enabled) {
        const int mi_offset = mi_row * parse_ctxt->frame_header->mi_cols + mi_col;
        const int bw        = mi_size_wide[bsize];
        const int bh        = mi_size_high[bsize];
        const int x_mis     = AOMMIN((int32_t)(parse_ctxt->frame_header->mi_cols - mi_col), bw);
        const int y_mis     = AOMMIN((int32_t)(parse_ctxt->frame_header->mi_rows - mi_row), bh);
        segment_id          = read_segment_id(dec_handle, parse_ctxt, xd, skip);
        set_segment_id(dec_handle, mi_offset, x_mis, y_mis, segment_id);
    }
    return segment_id;
}

static INLINE void update_palette_context(ParseCtxt *parse_ctx, int mi_row, int mi_col,
                                          BlockModeInfo *mi) {
    BlockSize             bsize     = mi->sb_type;
    ParseAboveNbr4x4Ctxt *above_ctx = parse_ctx->parse_above_nbr4x4_ctxt;
    ParseLeftNbr4x4Ctxt * left_ctx  = parse_ctx->parse_left_nbr4x4_ctxt;
    const int             bw        = mi_size_wide[bsize];
    const int             bh        = mi_size_high[bsize];
    for (int plane = 0; plane < MAX_MB_PLANE; plane++) {
        uint16_t *above_pal_col = above_ctx->above_palette_colors[plane] +
                                  (PALETTE_MAX_SIZE * ((mi_col - parse_ctx->sb_col_mi) % 16));
        uint16_t *left_pal_col = left_ctx->left_palette_colors[plane] +
                                 (PALETTE_MAX_SIZE * (mi_row - parse_ctx->sb_row_mi));
        for (int i = 0; i < bw; i++) {
            eb_memcpy(above_pal_col,
                   &parse_ctx->palette_colors[plane][0],
                   mi->palette_size[plane != 0] * sizeof(uint16_t));
            above_pal_col += PALETTE_MAX_SIZE;
        }
        for (int i = 0; i < bh; i++) {
            eb_memcpy(left_pal_col,
                   &parse_ctx->palette_colors[plane][0],
                   mi->palette_size[plane != 0] * sizeof(uint16_t));
            left_pal_col += PALETTE_MAX_SIZE;
        }
    }
}

static INLINE AomCdfProb *get_y_mode_cdf(FRAME_CONTEXT *tile_ctx, const BlockModeInfo *above_mi,
                                         const BlockModeInfo *left_mi) {
    const PredictionMode above     = above_mi ? above_mi->mode : DC_PRED;
    const PredictionMode left      = left_mi ? left_mi->mode : DC_PRED;
    const int            above_ctx = intra_mode_context[above];
    const int            left_ctx  = intra_mode_context[left];
    return tile_ctx->kf_y_cdf[above_ctx][left_ctx];
}

void intra_frame_mode_info(EbDecHandle *dec_handle, ParseCtxt *parse_ctxt, PartitionInfo *xd) {
    SvtReader *               r              = &parse_ctxt->r;
    BlockModeInfo *const      mbmi           = xd->mi;
    const BlockModeInfo *     above_mi       = xd->above_mbmi;
    const BlockModeInfo *     left_mi        = xd->left_mbmi;
    const BlockSize           bsize          = mbmi->sb_type;
    SegmentationParams *const seg            = &parse_ctxt->frame_header->segmentation_params;
    EbColorConfig             color_config   = parse_ctxt->seq_header->color_config;
    uint8_t *                 lossless_array = &parse_ctxt->frame_header->lossless_array[0];
    IntMv                     ref_mvs[INTRA_FRAME + 1][MAX_MV_REF_CANDIDATES] = {{{0}}};
    int16_t                   inter_mode_ctx[MODE_CTX_REF_FRAMES];
    MvCount *                 mv_cnt = (MvCount *)malloc(sizeof(MvCount));

    if (seg->seg_id_pre_skip) {
        mbmi->segment_id = intra_segment_id(dec_handle, parse_ctxt, xd, bsize, 0);
    }

    mbmi->skip = read_skip(parse_ctxt, xd, mbmi->segment_id);

    if (!seg->seg_id_pre_skip) {
        mbmi->segment_id = intra_segment_id(dec_handle, parse_ctxt, xd, bsize, mbmi->skip);
    }

    read_cdef(parse_ctxt, xd);

    read_delta_params(parse_ctxt, xd);
    parse_ctxt->read_deltas = 0;

    mbmi->ref_frame[0] = INTRA_FRAME;
    mbmi->ref_frame[1] = NONE_FRAME;

    mbmi->palette_size[0] = 0;
    mbmi->palette_size[1] = 0;

    mbmi->use_intrabc = 0;
    if (allow_intrabc(dec_handle))
        mbmi->use_intrabc = svt_read_symbol(r, parse_ctxt->cur_tile_ctx.intrabc_cdf, 2, ACCT_STR);

    mbmi->inter_inter_compound.type = COMPOUND_AVERAGE;

    if (mbmi->use_intrabc) {
        mbmi->mode           = DC_PRED;
        mbmi->uv_mode        = UV_DC_PRED;
        mbmi->motion_mode    = SIMPLE_TRANSLATION;
        mbmi->compound_mode  = COMPOUND_AVERAGE;
        mbmi->interp_filters = av1_broadcast_interp_filter(BILINEAR);
        IntMv global_mvs[2];
        av1_find_mv_refs(dec_handle,
                         xd,
                         parse_ctxt,
                         INTRA_FRAME,
                         xd->ref_mv_stack,
                         ref_mvs,
                         global_mvs,
                         inter_mode_ctx,
                         mv_cnt);

        assign_intrabc_mv(parse_ctxt, ref_mvs, xd);
    } else {
        AomCdfProb *y_mode_cdf = get_y_mode_cdf(&parse_ctxt->cur_tile_ctx, above_mi, left_mi);
        mbmi->mode             = read_intra_mode(r, y_mode_cdf);
        mbmi->angle_delta[PLANE_TYPE_Y] =
            intra_angle_info(r,
                             &parse_ctxt->cur_tile_ctx.angle_delta_cdf[mbmi->mode - V_PRED][0],
                             mbmi->mode,
                             bsize);
        if (xd->is_chroma_ref && !color_config.mono_chrome) {
            mbmi->uv_mode = read_intra_mode_uv(&parse_ctxt->cur_tile_ctx,
                                               r,
                                               is_cfl_allowed(xd, &color_config, lossless_array),
                                               mbmi->mode);
            if (mbmi->uv_mode == UV_CFL_PRED) {
                mbmi->cfl_alpha_idx =
                    read_cfl_alphas(&parse_ctxt->cur_tile_ctx, r, &mbmi->cfl_alpha_signs);
            }
            mbmi->angle_delta[PLANE_TYPE_UV] = intra_angle_info(
                r,
                &parse_ctxt->cur_tile_ctx.angle_delta_cdf[mbmi->uv_mode - V_PRED][0],
                get_uv_mode(mbmi->uv_mode),
                bsize);
        } else
            mbmi->uv_mode = UV_DC_PRED;

        if (allow_palette(parse_ctxt->frame_header->allow_screen_content_tools, bsize)) {
            palette_mode_info(parse_ctxt, xd);
            update_palette_context(parse_ctxt, xd->mi_row, xd->mi_col, mbmi);
        }
        filter_intra_mode_info(parse_ctxt, xd);
    }
    free(mv_cnt);
}

static INLINE int get_pred_context_seg_id(const PartitionInfo *xd) {
    const BlockModeInfo *const above_mi  = xd->above_mbmi;
    const BlockModeInfo *const left_mi   = xd->left_mbmi;
    const int                  above_sip = (above_mi != NULL) ? above_mi->seg_id_predicted : 0;
    const int                  left_sip  = (left_mi != NULL) ? left_mi->seg_id_predicted : 0;

    return above_sip + left_sip;
}

static INLINE void update_seg_ctx(ParseCtxt *parse_ctxt, int blk_col, int w4, int h4,
                                  int seg_id_predicted) {
    ParseAboveNbr4x4Ctxt *above_ctx = parse_ctxt->parse_above_nbr4x4_ctxt;
    ParseLeftNbr4x4Ctxt * left_ctx  = parse_ctxt->parse_left_nbr4x4_ctxt;

    uint8_t *const above_seg_ctx =
        above_ctx->above_seg_pred_ctx + (blk_col - parse_ctxt->cur_tile_info.mi_col_start);
    uint8_t *const left_seg_ctx = left_ctx->left_seg_pred_ctx;

    memset(above_seg_ctx, seg_id_predicted, w4);
    memset(left_seg_ctx, seg_id_predicted, h4);
}

static INLINE void copy_segment_id(ParseCtxt *parse_ctxt, const uint8_t *last_segment_ids,
                                   uint8_t *current_segment_ids, int mi_offset, int x_mis,
                                   int y_mis) {
    FrameHeader *frame_header = parse_ctxt->frame_header;
    for (int y = 0; y < y_mis; y++)
        for (int x = 0; x < x_mis; x++) {
            current_segment_ids[mi_offset + y * frame_header->mi_cols + x] =
                last_segment_ids ? last_segment_ids[mi_offset + y * frame_header->mi_cols + x] : 0;
        }
}

int read_inter_segment_id(EbDecHandle *dec_handle, ParseCtxt *parse_ctxt, PartitionInfo *xd,
                          int preskip) {
    SvtReader *          r            = &parse_ctxt->r;
    FrameHeader *        frame_header = parse_ctxt->frame_header;
    SegmentationParams * seg          = &frame_header->segmentation_params;
    BlockModeInfo *const mbmi         = xd->mi;
    uint32_t             mi_row       = xd->mi_row;
    uint32_t             mi_col       = xd->mi_col;
    const int            mi_offset    = mi_row * frame_header->mi_cols + mi_col;
    const uint32_t       bw           = mi_size_wide[mbmi->sb_type];
    const uint32_t       bh           = mi_size_high[mbmi->sb_type];

    const int x_mis = AOMMIN(frame_header->mi_cols - mi_col, bw);
    const int y_mis = AOMMIN(frame_header->mi_rows - mi_row, bh);

    if (!seg->segmentation_enabled) return 0; // Default for disabled segmentation
    if (!seg->segmentation_update_map) {
        copy_segment_id(parse_ctxt,
                        dec_handle->cm
                            .last_frame_seg_map,
                        dec_handle->cur_pic_buf[0]->segment_maps,
                        mi_offset,
                        x_mis,
                        y_mis);
        return dec_handle->cm.last_frame_seg_map
                   ? get_segment_id(
                         frame_header, dec_handle->cm.last_frame_seg_map, mbmi->sb_type, mi_row, mi_col)
                   : 0;
    }

    int segment_id;
    if (preskip) {
        if (!seg->seg_id_pre_skip) return 0;
    } else {
        if (mbmi->skip) {
            mbmi->seg_id_predicted = 0;
            update_seg_ctx(parse_ctxt, mi_col, bw, bh, mbmi->seg_id_predicted);
            segment_id = read_segment_id(dec_handle, parse_ctxt, xd, 1);
            set_segment_id(dec_handle, mi_offset, x_mis, y_mis, segment_id);
            return segment_id;
        }
    }

    if (seg->segmentation_temporal_update) {
        const int                        ctx  = get_pred_context_seg_id(xd);
        struct segmentation_probs *const segp = &parse_ctxt->cur_tile_ctx.seg;
        mbmi->seg_id_predicted = svt_read_symbol(r, segp->pred_cdf[ctx], 2, ACCT_STR);
        if (mbmi->seg_id_predicted) {
            segment_id =
                dec_handle->cm.last_frame_seg_map
                    ? get_segment_id(
                          frame_header, dec_handle->cm.last_frame_seg_map, mbmi->sb_type, mi_row, mi_col)
                    : 0;
        } else
            segment_id = read_segment_id(dec_handle, parse_ctxt, xd, 0);
        update_seg_ctx(parse_ctxt, mi_col, bw, bh, mbmi->seg_id_predicted);
    } else
        segment_id = read_segment_id(dec_handle, parse_ctxt, xd, 0);
    set_segment_id(dec_handle, mi_offset, x_mis, y_mis, segment_id);

    return segment_id;
}

static int get_block_position(FrameHeader *frame_info, int *mi_r, int *mi_c, int blk_row,
                              int blk_col, MV mv, int sign_bias) {
    const int32_t base_blk_row = (blk_row >> 3) << 3;
    const int32_t base_blk_col = (blk_col >> 3) << 3;

    const int row_offset =
        (mv.row >= 0) ? (mv.row >> (4 + MI_SIZE_LOG2)) : -((-mv.row) >> (4 + MI_SIZE_LOG2));

    const int col_offset =
        (mv.col >= 0) ? (mv.col >> (4 + MI_SIZE_LOG2)) : -((-mv.col) >> (4 + MI_SIZE_LOG2));

    const int32_t row = (sign_bias == 1) ? blk_row - row_offset : blk_row + row_offset;
    const int32_t col = (sign_bias == 1) ? blk_col - col_offset : blk_col + col_offset;

    if (row < 0 || (uint32_t)row >= (frame_info->mi_rows >> 1) || col < 0 ||
        (uint32_t)col >= (frame_info->mi_cols >> 1))
        return 0;

    if (row < base_blk_row - (MAX_OFFSET_HEIGHT >> 3) ||
        row >= base_blk_row + 8 + (MAX_OFFSET_HEIGHT >> 3) ||
        col < base_blk_col - (MAX_OFFSET_WIDTH >> 3) ||
        col >= base_blk_col + 8 + (MAX_OFFSET_WIDTH >> 3))
        return 0;

    *mi_r = row;
    *mi_c = col;

    return 1;
}

int motion_field_projection_check(EbDecHandle *dec_handle, MvReferenceFrame start_frame) {
    FrameHeader *            frame_info      = &dec_handle->frame_header;
    const EbDecPicBuf *const start_frame_buf = get_ref_frame_buf(dec_handle, start_frame);

    if (start_frame_buf == NULL) return 0;

    if (start_frame_buf->frame_type == KEY_FRAME || start_frame_buf->frame_type == INTRA_ONLY_FRAME)
        return 0;

    uint32_t mi_cols = 2 * ((start_frame_buf->frame_width + 7) >> 3);
    uint32_t mi_rows = 2 * ((start_frame_buf->frame_height + 7) >> 3);
    if (mi_rows != frame_info->mi_rows || mi_cols != frame_info->mi_cols) return 0;
    return 1;
}

// Note: motion_filed_projection finds motion vectors of current frame's
// reference frame, and projects them to current frame. To make it clear,
// let's call current frame's reference frame as start frame.
// Call Start frame's reference frames as reference frames.
// Call ref_offset as frame distances between start frame and its reference
// frames.
int motion_field_projection_row(EbDecHandle *dec_handle, MvReferenceFrame start_frame, int sb_row,
                                int num_blk_mv_rows, int dir) {
    FrameHeader *frame_info = &dec_handle->frame_header;

    if (!motion_field_projection_check(dec_handle, start_frame)) return 0;

    const EbDecPicBuf *const start_frame_buf = get_ref_frame_buf(dec_handle, start_frame);

    TemporalMvRef *mv_ref_base = start_frame_buf->mvs;

    const int start_frame_order_hint = start_frame_buf->order_hint;

    const unsigned int *const ref_order_hints = &start_frame_buf->ref_order_hints[0];

    const int cur_order_hint = dec_handle->cur_pic_buf[0]->order_hint;

    int start_to_current_frame_offset = get_relative_dist(
        &dec_handle->seq_header.order_hint_info, start_frame_order_hint, cur_order_hint);

    if (dir == 2) start_to_current_frame_offset = -start_to_current_frame_offset;

    int       ref_offset[REF_FRAMES] = {0};
    const int mvs_cols               = (frame_info->mi_cols + 1) >> 1; //8x8 unit level

    TemporalMvRef *tpl_mvs_base = dec_handle->master_frame_buf.tpl_mvs;

    for (MvReferenceFrame rf = LAST_FRAME; rf <= INTER_REFS_PER_FRAME; ++rf) {
        ref_offset[rf] = get_relative_dist(&dec_handle->seq_header.order_hint_info,
                                           start_frame_order_hint,
                                           ref_order_hints[rf - LAST_FRAME]);
    }

    for (int blk_row = 0; blk_row < num_blk_mv_rows; ++blk_row) {
        for (int blk_col = 0; blk_col < mvs_cols; ++blk_col) {
            int blk_row_in_frm = (sb_row << 3) + blk_row;

            TemporalMvRef *mv_ref = &mv_ref_base[blk_row_in_frm * mvs_cols + blk_col];
            MV             fwd_mv = mv_ref->mf_mv0.as_mv;

            if (mv_ref->ref_frame_offset > INTRA_FRAME) {
                IntMv     this_mv;
                int       mi_r, mi_c;
                const int ref_frame_offset = ref_offset[mv_ref->ref_frame_offset];

                int pos_valid = abs(ref_frame_offset) <= MAX_FRAME_DISTANCE &&
                                ref_frame_offset > 0 &&
                                abs(start_to_current_frame_offset) <= MAX_FRAME_DISTANCE;

                if (pos_valid) {
                    get_mv_projection(
                        &this_mv.as_mv, fwd_mv, start_to_current_frame_offset, ref_frame_offset);

                    pos_valid = get_block_position(
                        frame_info, &mi_r, &mi_c, blk_row_in_frm, blk_col, this_mv.as_mv, dir >> 1);
                }

                if (pos_valid) {
                    const int mi_offset = mi_r * (frame_info->mi_stride >> 1) + mi_c;
                    tpl_mvs_base[mi_offset].mf_mv0.as_int = mv_ref->mf_mv0.as_int;/*Same as fwd_mv*/
                    tpl_mvs_base[mi_offset].ref_frame_offset = ref_frame_offset;
                }
            }
        }
    }
    return 1;
}

/* motion_field_projection_row for all reference frames */
void motion_field_projections_row(EbDecHandle *dec_handle, int sb_row, const EbDecPicBuf **ref_buf,
                                  int *ref_order_hint) {
    OrderHintInfo *order_hint_info = &dec_handle->seq_header.order_hint_info;
    TemporalMvRef *tpl_mvs_base    = dec_handle->master_frame_buf.tpl_mvs;

    const int cur_order_hint = dec_handle->cur_pic_buf[0]->order_hint;
    const int mvs_rows       = (dec_handle->frame_header.mi_rows + 1) >> 1; //8x8 unit level
    const int sb_mvs_rows    = (mvs_rows + 7) >> 3; //64x64 unit level

    int num_blk_mv_rows = 8; //8 8x8 rows in 64x64

    // if last row then calcaute the actual mv blk rows
    if (sb_row == (sb_mvs_rows - 1)) num_blk_mv_rows = mvs_rows - (sb_row << 3);

    int offset = (sb_row << 3) * (dec_handle->frame_header.mi_stride >> 1);

    int size = 8 * (dec_handle->frame_header.mi_stride >> 1); //8 8x8 rows in 64x64

    for (int idx = 0; idx < size; ++idx) {
        tpl_mvs_base[offset + idx].mf_mv0.as_int    = INVALID_MV;
        tpl_mvs_base[offset + idx].ref_frame_offset = 0;
    }

    int ref_stamp = MFMV_STACK_SIZE - 1;

    if (ref_buf[LAST_FRAME - LAST_FRAME] != NULL) {
        const int alt_of_lst_order_hint =
            ref_buf[LAST_FRAME - LAST_FRAME]->ref_order_hints[ALTREF_FRAME - LAST_FRAME];

        const int is_lst_overlay =
            (alt_of_lst_order_hint == ref_order_hint[GOLDEN_FRAME - LAST_FRAME]);
        if (!is_lst_overlay)
            motion_field_projection_row(dec_handle, LAST_FRAME, sb_row, num_blk_mv_rows, 2);
        --ref_stamp;
    }

    if (get_relative_dist(
            order_hint_info, ref_order_hint[BWDREF_FRAME - LAST_FRAME], cur_order_hint) > 0) {
        if (motion_field_projection_row(dec_handle, BWDREF_FRAME, sb_row, num_blk_mv_rows, 0))
            --ref_stamp;
    }

    if (get_relative_dist(
            order_hint_info, ref_order_hint[ALTREF2_FRAME - LAST_FRAME], cur_order_hint) > 0) {
        if (motion_field_projection_row(dec_handle, ALTREF2_FRAME, sb_row, num_blk_mv_rows, 0))
            --ref_stamp;
    }

    if (get_relative_dist(
            order_hint_info, ref_order_hint[ALTREF_FRAME - LAST_FRAME], cur_order_hint) > 0 &&
        ref_stamp >= 0)
        if (motion_field_projection_row(dec_handle, ALTREF_FRAME, sb_row, num_blk_mv_rows, 0))
            --ref_stamp;

    if (ref_stamp >= 0)
        motion_field_projection_row(dec_handle, LAST2_FRAME, sb_row, num_blk_mv_rows, 2);
}

void svt_setup_motion_field(EbDecHandle *dec_handle, DecThreadCtxt *thread_ctxt) {
    DecMtFrameData *dec_mt_frame_data =
        &dec_handle->master_frame_buf.cur_frame_bufs[0].dec_mt_frame_data;
    EbBool is_mt = dec_handle->dec_config.threads > 1;
    EbBool do_memset = EB_TRUE;

    if (is_mt) {
        volatile EbBool *start_motion_proj = &dec_mt_frame_data->start_motion_proj;

        while (*start_motion_proj != EB_TRUE)
            eb_block_on_semaphore(NULL == thread_ctxt ? dec_handle->thread_semaphore
                                                      : thread_ctxt->thread_semaphore);

        DecMtMotionProjInfo *motion_proj_info =
            &dec_mt_frame_data->motion_proj_info;
        do_memset = EB_FALSE;
        //lock mutex
        eb_block_on_mutex(motion_proj_info->motion_proj_mutex);
        //Check memset needed
        if (motion_proj_info->motion_proj_init_done == EB_FALSE) {
            do_memset = EB_TRUE;
            motion_proj_info->motion_proj_init_done = EB_TRUE;
        }
        //unlock mutex
        eb_release_mutex(motion_proj_info->motion_proj_mutex);
    }

    if (do_memset) {
        memset(dec_handle->master_frame_buf.ref_frame_side,
               0,
               sizeof(dec_handle->master_frame_buf.ref_frame_side));
    }

    EbBool no_proj_flag = (dec_handle->frame_header.show_existing_frame ||
                           (0 == dec_handle->frame_header.use_ref_frame_mvs));

    OrderHintInfo *order_hint_info = &dec_handle->seq_header.order_hint_info;

    if (!order_hint_info->enable_order_hint) return;

    const int cur_order_hint = dec_handle->cur_pic_buf[0]->order_hint;

    const EbDecPicBuf *ref_buf[INTER_REFS_PER_FRAME];
    int                ref_order_hint[INTER_REFS_PER_FRAME];

    for (int ref_frame = LAST_FRAME; ref_frame <= ALTREF_FRAME; ref_frame++) {
        const int                ref_idx    = ref_frame - LAST_FRAME;
        const EbDecPicBuf *const buf        = get_ref_frame_buf(dec_handle, ref_frame);
        int                      order_hint = 0;

        if (buf != NULL) order_hint = buf->order_hint;

        ref_buf[ref_idx]        = buf;
        ref_order_hint[ref_idx] = order_hint;

        if (get_relative_dist(order_hint_info, order_hint, cur_order_hint) > 0)
            dec_handle->master_frame_buf.ref_frame_side[ref_frame] = 1;
        else if (order_hint == cur_order_hint)
            dec_handle->master_frame_buf.ref_frame_side[ref_frame] = -1;
    }
    if (!no_proj_flag) {
        //branch of point for MT
        if (is_mt) {
            DecMtMotionProjInfo *motion_proj_info =
                &dec_handle->master_frame_buf.cur_frame_bufs[0].dec_mt_frame_data.motion_proj_info;

            while (1) {
                int32_t proj_row = -1;

                //lock mutex
                eb_block_on_mutex(motion_proj_info->motion_proj_mutex);

                //pick up a row and increment the sb row counter
                if (motion_proj_info->motion_proj_row_to_process !=
                    motion_proj_info->num_motion_proj_rows) {
                    proj_row = motion_proj_info->motion_proj_row_to_process;
                    motion_proj_info->motion_proj_row_to_process++;
                }

                //unlock mutex
                eb_release_mutex(motion_proj_info->motion_proj_mutex);

                //wait for parse
                if (-1 != proj_row)
                    motion_field_projections_row(dec_handle, proj_row, ref_buf, ref_order_hint);

                /*if all sb rows have been picked up for processing then break the while loop */
                if (motion_proj_info->motion_proj_row_to_process ==
                    motion_proj_info->num_motion_proj_rows) {
                    break;
                }
            }
        } else {
            const int mvs_rows    = (dec_handle->frame_header.mi_rows + 1) >> 1; //8x8 unit level
            const int sb_mvs_rows = (mvs_rows + 7) >> 3; //64x64 unit level
            for (int sb_row = 0; sb_row < sb_mvs_rows; sb_row++)
                motion_field_projections_row(dec_handle, sb_row, ref_buf, ref_order_hint);
        }
    }

    if (is_mt) {
        eb_block_on_mutex(dec_mt_frame_data->temp_mutex);
        dec_mt_frame_data->num_threads_header++;
        if (dec_handle->dec_config.threads == dec_mt_frame_data->num_threads_header) {
            dec_mt_frame_data->start_motion_proj = EB_FALSE;
        }
        eb_release_mutex(dec_mt_frame_data->temp_mutex);

        volatile uint32_t *num_threads_header = &dec_mt_frame_data->num_threads_header;
        while (*num_threads_header != dec_handle->dec_config.threads &&
              (EB_FALSE == dec_mt_frame_data->end_flag))
            ;
    }
}

void intra_block_mode_info(ParseCtxt *parse_ctxt, PartitionInfo *xd) {
    SvtReader *     r     = &parse_ctxt->r;
    BlockModeInfo * mbmi  = xd->mi;
    const BlockSize bsize = mbmi->sb_type;
    mbmi->ref_frame[0]    = INTRA_FRAME;
    mbmi->ref_frame[1]    = NONE_FRAME;

    EbColorConfig *color_cfg      = &parse_ctxt->seq_header->color_config;
    uint8_t *      lossless_array = &parse_ctxt->frame_header->lossless_array[0];

    mbmi->mode = read_intra_mode(r, parse_ctxt->cur_tile_ctx.y_mode_cdf[size_group_lookup[bsize]]);

    mbmi->angle_delta[PLANE_TYPE_Y] = intra_angle_info(
        r, &parse_ctxt->cur_tile_ctx.angle_delta_cdf[mbmi->mode - V_PRED][0], mbmi->mode, bsize);

    if (xd->is_chroma_ref && !color_cfg->mono_chrome) {
        mbmi->uv_mode = read_intra_mode_uv(&parse_ctxt->cur_tile_ctx,
                                           r,
                                           is_cfl_allowed(xd, color_cfg, lossless_array),
                                           mbmi->mode);
        if (mbmi->uv_mode == UV_CFL_PRED) {
            mbmi->cfl_alpha_idx =
                read_cfl_alphas(&parse_ctxt->cur_tile_ctx, r, &mbmi->cfl_alpha_signs);
        }
        mbmi->angle_delta[PLANE_TYPE_UV] =
            intra_angle_info(r,
                             &parse_ctxt->cur_tile_ctx.angle_delta_cdf[mbmi->uv_mode - V_PRED][0],
                             get_uv_mode(mbmi->uv_mode),
                             bsize);
    }

    mbmi->palette_size[0] = 0;
    mbmi->palette_size[1] = 0;

    if (allow_palette(parse_ctxt->frame_header->allow_screen_content_tools, bsize)) {
        palette_mode_info(parse_ctxt, xd);
        update_palette_context(parse_ctxt, xd->mi_row, xd->mi_col, mbmi);
    }

    filter_intra_mode_info(parse_ctxt, xd);
}

int read_is_inter(ParseCtxt *parse_ctxt, PartitionInfo *xd, int segment_id) {
    SvtReader *         r          = &parse_ctxt->r;
    int                 is_inter   = 0;
    SegmentationParams *seg_params = &parse_ctxt->frame_header->segmentation_params;
    if (seg_feature_active(seg_params, segment_id, SEG_LVL_REF_FRAME))
        is_inter = get_segdata(seg_params, segment_id, SEG_LVL_REF_FRAME) != INTRA_FRAME;
    else if (seg_feature_active(seg_params, segment_id, SEG_LVL_GLOBALMV))
        is_inter = 1;
    else {
        const int ctx = get_intra_inter_context(xd);
        is_inter = svt_read_symbol(r, parse_ctxt->cur_tile_ctx.intra_inter_cdf[ctx], 2, ACCT_STR);
    }
    return is_inter;
}

void inter_frame_mode_info(EbDecHandle *dec_handle, ParseCtxt *parse_ctxt, PartitionInfo *pi) {
    BlockModeInfo *mbmi = pi->mi;
    mbmi->use_intrabc   = 0;
    int inter_block     = 1;

    mbmi->mv[0].as_int = 0;
    mbmi->mv[1].as_int = 0;

    mbmi->inter_inter_compound.type = COMPOUND_AVERAGE;

    mbmi->segment_id = read_inter_segment_id(dec_handle, parse_ctxt, pi, 1);

    mbmi->skip_mode = read_skip_mode(parse_ctxt, pi, mbmi->segment_id);

    if (mbmi->skip_mode)
        mbmi->skip = 1;
    else
        mbmi->skip = read_skip(parse_ctxt, pi, mbmi->segment_id);

    if (!parse_ctxt->frame_header->segmentation_params.seg_id_pre_skip) {
        mbmi->segment_id = read_inter_segment_id(dec_handle, parse_ctxt, pi, 0);
    }

    read_cdef(parse_ctxt, pi);

    read_delta_params(parse_ctxt, pi);
    parse_ctxt->read_deltas = 0;

    if (!mbmi->skip_mode) inter_block = read_is_inter(parse_ctxt, pi, mbmi->segment_id);

    if (inter_block)
        inter_block_mode_info(dec_handle, parse_ctxt, pi);
    else
        intra_block_mode_info(parse_ctxt, pi);
}

static void intra_copy_frame_mvs(EbDecHandle *dec_handle, int mi_row, int mi_col, int x_mis,
                                 int y_mis) {
    FrameHeader *  frame_info       = &dec_handle->frame_header;
    const int      frame_mvs_stride = ROUND_POWER_OF_TWO(frame_info->mi_cols, 1);
    TemporalMvRef *frame_mvs =
        dec_handle->cur_pic_buf[0]->mvs + (mi_row >> 1) * frame_mvs_stride + (mi_col >> 1);
    x_mis = ROUND_POWER_OF_TWO(x_mis, 1);
    y_mis = ROUND_POWER_OF_TWO(y_mis, 1);

    for (int h = 0; h < y_mis; h++) {
        TemporalMvRef *mv = frame_mvs;
        for (int w = 0; w < x_mis; w++) {
            mv->ref_frame_offset = NONE_FRAME;
            mv++;
        }
        frame_mvs += frame_mvs_stride;
    }
}

void inter_copy_frame_mvs(EbDecHandle *dec_handle, BlockModeInfo *mi, int mi_row, int mi_col,
                          int x_mis, int y_mis) {
    FrameHeader *  frame_info       = &dec_handle->frame_header;
    const int      frame_mvs_stride = ROUND_POWER_OF_TWO(frame_info->mi_cols, 1);
    TemporalMvRef *frame_mvs =
        dec_handle->cur_pic_buf[0]->mvs + (mi_row >> 1) * frame_mvs_stride + (mi_col >> 1);
    x_mis = ROUND_POWER_OF_TWO(x_mis, 1);
    y_mis = ROUND_POWER_OF_TWO(y_mis, 1);

    TemporalMvRef *mv;
    TemporalMvRef  cur_mv;
    cur_mv.ref_frame_offset = NONE_FRAME;
    cur_mv.mf_mv0.as_int    = 0;

    for (int idx = 0; idx < 2; ++idx) {
        MvReferenceFrame ref_frame = mi->ref_frame[idx];
        if (ref_frame > INTRA_FRAME) {
            int8_t ref_idx = dec_handle->master_frame_buf.ref_frame_side[ref_frame];
            if (ref_idx) continue;
            if ((abs(mi->mv[idx].as_mv.row) > REFMVS_LIMIT) ||
                (abs(mi->mv[idx].as_mv.col) > REFMVS_LIMIT))
                continue;
            cur_mv.ref_frame_offset = ref_frame;
            cur_mv.mf_mv0.as_int    = mi->mv[idx].as_int;
        }
    }
    for (int h = 0; h < y_mis; h++) {
        mv = frame_mvs;
        for (int w = 0; w < x_mis; w++) {
            *mv = cur_mv;
            mv++;
        }
        frame_mvs += frame_mvs_stride;
    }
}

void mode_info(EbDecHandle *dec_handle, PartitionInfo *part_info, ParseCtxt *parse_ctxt) {
    BlockModeInfo *mi         = part_info->mi;
    FrameHeader *  frame_info = parse_ctxt->frame_header;
    //BlockSize bsize = mi->sb_type
    mi->use_intrabc = 0;
    mi->segment_id  = 0;

    uint16_t       mi_row = part_info->mi_row;
    uint16_t       mi_col = part_info->mi_col;
    const uint32_t bw     = mi_size_wide[mi->sb_type];
    const uint32_t bh     = mi_size_high[mi->sb_type];
    const int      x_mis  = AOMMIN(bw, frame_info->mi_cols - mi_col);
    const int      y_mis  = AOMMIN(bh, frame_info->mi_rows - mi_row);

    if (frame_info->frame_type == KEY_FRAME || frame_info->frame_type == INTRA_ONLY_FRAME) {
        intra_frame_mode_info(dec_handle, parse_ctxt, part_info);
        intra_copy_frame_mvs(dec_handle, mi_row, mi_col, x_mis, y_mis);
    } else {
        inter_frame_mode_info(dec_handle, parse_ctxt, part_info);
        inter_copy_frame_mvs(dec_handle, mi, mi_row, mi_col, x_mis, y_mis);
    }
}

TxSize read_tx_size(ParseCtxt *parse_ctxt, PartitionInfo *xd, int allow_select) {
    BlockModeInfo * mbmi    = xd->mi;
    const TxMode    tx_mode = parse_ctxt->frame_header->tx_mode;
    const BlockSize bsize   = xd->mi->sb_type;
    if (parse_ctxt->frame_header->lossless_array[mbmi->segment_id]) return TX_4X4;

    if (bsize > BLOCK_4X4 && allow_select && tx_mode == TX_MODE_SELECT) {
        const TxSize coded_tx_size = read_selected_tx_size(xd, parse_ctxt);
        return coded_tx_size;
    }
    assert(IMPLIES(tx_mode == ONLY_4X4, bsize == BLOCK_4X4));
    TxSize tx_size = max_txsize_rect_lookup[bsize];
    return tx_size;
}

/* Update Chroma Transform Info for Inter Case! */
void update_chroma_trans_info(ParseCtxt *parse_ctx, PartitionInfo *part_info, BlockSize bsize) {
    BlockModeInfo *mbmi         = part_info->mi;
    SBInfo *       sb_info      = part_info->sb_info;
    EbColorConfig  color_config = parse_ctx->seq_header->color_config;

    int              num_chroma_tus = 0, step_r, step_c, force_split_cnt = 0, total_chroma_tus = 0;
    int              sx = color_config.subsampling_x;
    int              sy = color_config.subsampling_y;
    TransformInfo_t *chroma_trans_info =
        sb_info->sb_trans_info[AOM_PLANE_U] + mbmi->first_txb_offset[AOM_PLANE_U];

    const int       max_blocks_wide = max_block_wide(part_info, bsize, 0);
    const int       max_blocks_high = max_block_high(part_info, bsize, 0);
    const BlockSize max_unit_bsize  = BLOCK_64X64;
    int             width           = block_size_wide[max_unit_bsize] >> tx_size_wide_log2[0];
    int             height          = block_size_high[max_unit_bsize] >> tx_size_high_log2[0];
    width                           = AOMMIN(width, max_blocks_wide);
    height                          = AOMMIN(height, max_blocks_high);

    TxSize tx_size_uv = av1_get_max_uv_txsize(bsize, sx, sy);
    assert(parse_ctx->frame_header->lossless_array[mbmi->segment_id] != 1);

    /* TODO: Make plane loop and avoid the unroll */
    for (int idy = 0; idy < max_blocks_high; idy += height) {
        for (int idx = 0; idx < max_blocks_wide; idx += width, force_split_cnt++) {
            num_chroma_tus = 0;

            if (color_config.mono_chrome) continue;

            /* Update Chroma Transform Info */
            if (!part_info->is_chroma_ref) continue;

            step_r = tx_size_high_unit[tx_size_uv];
            step_c = tx_size_wide_unit[tx_size_uv];

            /* UV dim. of 4 special case! */
            int unit_height = ROUND_POWER_OF_TWO(AOMMIN(height + idy, max_blocks_high), sy);
            int unit_width  = ROUND_POWER_OF_TWO(AOMMIN(width + idx, max_blocks_wide), sx);

            /* TODO : Can cause prblm for incomplete SBs. Fix! */
            for (int blk_row = idy >> sy; blk_row < unit_height; blk_row += step_r) {
                for (int blk_col = idx >> sx; blk_col < unit_width; blk_col += step_c) {
                    // Chroma Cb
                    chroma_trans_info->tx_size      = tx_size_uv;
                    chroma_trans_info->txb_x_offset = blk_col;
                    chroma_trans_info->txb_y_offset = blk_row;
                    chroma_trans_info++;
                    num_chroma_tus++;
                }
            }

            parse_ctx->num_tus[AOM_PLANE_U][force_split_cnt] = num_chroma_tus;
            parse_ctx->num_tus[AOM_PLANE_V][force_split_cnt] = num_chroma_tus;
        }
    }

    total_chroma_tus = parse_ctx->num_tus[AOM_PLANE_U][0] + parse_ctx->num_tus[AOM_PLANE_U][1] +
                       parse_ctx->num_tus[AOM_PLANE_U][2] + parse_ctx->num_tus[AOM_PLANE_U][3];
    /* Cr Transform Info Update from Cb */
    if (total_chroma_tus) {
        assert((chroma_trans_info - total_chroma_tus) ==
               sb_info->sb_trans_info[AOM_PLANE_U] + mbmi->first_txb_offset[AOM_PLANE_U]);
        eb_memcpy(chroma_trans_info,
               chroma_trans_info - total_chroma_tus,
               total_chroma_tus * sizeof(*chroma_trans_info));
    }
    mbmi->num_tus[AOM_PLANE_U] = total_chroma_tus;
    parse_ctx->first_txb_offset[AOM_PLANE_U] += 2 * total_chroma_tus;
}

static INLINE TxSize find_tx_size(int w, int h) {
    int tx_sz;
    for (tx_sz = 0; tx_sz < TX_SIZES_ALL; tx_sz++)
        if (tx_size_wide[tx_sz] == w && tx_size_high[tx_sz] == h) break;
    return tx_sz;
}

int get_txfm_split_ctx(PartitionInfo *pi, ParseCtxt *parse_ctx, TxSize tx_size, int blk_row,
                       int blk_col) {
    int above = parse_ctx->parse_above_nbr4x4_ctxt
                    ->above_tx_wd[pi->mi_col - parse_ctx->cur_tile_info.mi_col_start + blk_col] <
                tx_size_wide[tx_size];
    int left =
        parse_ctx->parse_left_nbr4x4_ctxt->left_tx_ht[pi->mi_row - parse_ctx->sb_row_mi + blk_row] <
        tx_size_high[tx_size];
    int size = MIN(64, MAX(block_size_wide[pi->mi->sb_type], block_size_high[pi->mi->sb_type]));

    TxSize max_tx_size  = find_tx_size(size, size);
    TxSize tx_sz_sqr_up = txsize_sqr_up_map[tx_size];
    return ((tx_sz_sqr_up != max_tx_size) * 3 + (TX_SIZES - 1 - max_tx_size) * 6 + above + left);
}

void read_var_tx_size(ParseCtxt *parse_ctx, PartitionInfo *pi, TxSize tx_size, int blk_row,
                      int blk_col, int depth, int *num_luma_tus) {
    SvtReader *     r               = &parse_ctx->r;
    BlockModeInfo * mbmi            = pi->mi;
    const BlockSize bsize           = mbmi->sb_type;
    const int       max_blocks_high = max_block_high(pi, bsize, 0);
    const int       max_blocks_wide = max_block_wide(pi, bsize, 0);
    int             i, j, txfm_split = 0;

    if (blk_row >= max_blocks_high || blk_col >= max_blocks_wide) return;

    if (tx_size == TX_4X4 || depth == MAX_VARTX_DEPTH)
        txfm_split = 0;
    else {
        int ctx = get_txfm_split_ctx(pi, parse_ctx, tx_size, blk_row, blk_col);
        txfm_split =
            svt_read_symbol(r, parse_ctx->cur_tile_ctx.txfm_partition_cdf[ctx], 2, ACCT_STR);
    }

    int w4 = tx_size_wide_unit[tx_size];
    int h4 = tx_size_high_unit[tx_size];
    if (txfm_split) {
        TxSize sub_tx_sz = sub_tx_size_map[tx_size];
        int    step_w    = tx_size_wide_unit[sub_tx_sz];
        int    step_h    = tx_size_high_unit[sub_tx_sz];

        for (i = 0; i < h4; i += step_h)
            for (j = 0; j < w4; j += step_w)
                read_var_tx_size(
                    parse_ctx, pi, sub_tx_sz, blk_row + i, blk_col + j, depth + 1, num_luma_tus);
    } else {
        parse_ctx->cur_luma_trans_info->tx_size      = tx_size;
        parse_ctx->cur_luma_trans_info->txb_x_offset = blk_col;
        parse_ctx->cur_luma_trans_info->txb_y_offset = blk_row;
        parse_ctx->cur_luma_trans_info++;
        parse_ctx->cur_blk_luma_count++;
        *num_luma_tus += 1;
        update_tx_context(parse_ctx, pi, bsize, tx_size, blk_row, blk_col);
    }
}

/* Update Flat Transform Info for Intra Case! */
void update_flat_trans_info(ParseCtxt *parse_ctx, PartitionInfo *part_info, BlockSize bsize,
                            TxSize tx_size) {
    BlockModeInfo *mbmi         = part_info->mi;
    SBInfo *       sb_info      = part_info->sb_info;
    EbColorConfig  color_config = parse_ctx->seq_header->color_config;

    int sx = color_config.subsampling_x;
    int sy = color_config.subsampling_y;
    int num_luma_tus, num_chroma_tus, force_split_cnt = 0, total_luma_tus = 0, total_chroma_tus = 0;
    TransformInfo_t *luma_trans_info =
        sb_info->sb_trans_info[AOM_PLANE_Y] + mbmi->first_txb_offset[AOM_PLANE_Y];
    TransformInfo_t *chroma_trans_info =
        sb_info->sb_trans_info[AOM_PLANE_U] + mbmi->first_txb_offset[AOM_PLANE_U];

    const int       max_blocks_wide = max_block_wide(part_info, bsize, 0);
    const int       max_blocks_high = max_block_high(part_info, bsize, 0);
    const BlockSize max_unit_bsize  = BLOCK_64X64;
    int             width           = block_size_wide[max_unit_bsize] >> tx_size_wide_log2[0];
    int             height          = block_size_high[max_unit_bsize] >> tx_size_high_log2[0];
    width                           = AOMMIN(width, max_blocks_wide);
    height                          = AOMMIN(height, max_blocks_high);

    TxSize tx_size_uv = parse_ctx->frame_header->lossless_array[mbmi->segment_id]
                            ? TX_4X4
                            : av1_get_max_uv_txsize(bsize, sx, sy);

    /* TODO: Make plane loop and avoid the unroll */
    for (int idy = 0; idy < max_blocks_high; idy += height) {
        for (int idx = 0; idx < max_blocks_wide; idx += width, force_split_cnt++) {
            num_luma_tus   = 0;
            num_chroma_tus = 0;

            /* Update Luma Transform Info */
            int step_r = tx_size_high_unit[tx_size];
            int step_c = tx_size_wide_unit[tx_size];

            /* Y dim. of 4 special case! */
            int unit_height = ROUND_POWER_OF_TWO(AOMMIN(height + idy, max_blocks_high), 0);
            int unit_width  = ROUND_POWER_OF_TWO(AOMMIN(width + idx, max_blocks_wide), 0);

            /* TODO : Can cause prblm for incomplete SBs. Fix! */
            for (int blk_row = idy; blk_row < unit_height; blk_row += step_r) {
                for (int blk_col = idx; blk_col < unit_width; blk_col += step_c) {
                    luma_trans_info->tx_size      = tx_size;
                    luma_trans_info->txb_x_offset = blk_col;
                    luma_trans_info->txb_y_offset = blk_row;
                    luma_trans_info++;
                    num_luma_tus++;
                    total_luma_tus++;
                }
            }
            parse_ctx->num_tus[AOM_PLANE_Y][force_split_cnt] = num_luma_tus;

            if (color_config.mono_chrome) continue;

            /* Update Chroma Transform Info */
            if (!part_info->is_chroma_ref) continue;

            step_r = tx_size_high_unit[tx_size_uv];
            step_c = tx_size_wide_unit[tx_size_uv];

            /* UV dim. of 4 special case! */
            unit_height = ROUND_POWER_OF_TWO(AOMMIN(height + idy, max_blocks_high), sy);
            unit_width  = ROUND_POWER_OF_TWO(AOMMIN(width + idx, max_blocks_wide), sx);
            /* TODO : Can cause prblm for incomplete SBs. Fix! */
            for (int blk_row = idy >> sy; blk_row < unit_height; blk_row += step_r) {
                for (int blk_col = idx >> sx; blk_col < unit_width; blk_col += step_c) {
                    chroma_trans_info->tx_size      = tx_size_uv;
                    chroma_trans_info->txb_x_offset = blk_col;
                    chroma_trans_info->txb_y_offset = blk_row;
                    chroma_trans_info++;
                    num_chroma_tus++;
                    total_chroma_tus++;
                }
            }

            parse_ctx->num_tus[AOM_PLANE_U][force_split_cnt] = num_chroma_tus;
            parse_ctx->num_tus[AOM_PLANE_V][force_split_cnt] = num_chroma_tus;
        }
    }

    /* Cr Transform Info Update from Cb */
    if (total_chroma_tus) {
        assert((chroma_trans_info - total_chroma_tus) ==
               sb_info->sb_trans_info[AOM_PLANE_U] + mbmi->first_txb_offset[AOM_PLANE_U]);
        eb_memcpy(chroma_trans_info,
               chroma_trans_info - total_chroma_tus,
               total_chroma_tus * sizeof(*chroma_trans_info));
    }

    mbmi->num_tus[AOM_PLANE_Y] = total_luma_tus;
    mbmi->num_tus[AOM_PLANE_U] = total_chroma_tus;

    parse_ctx->first_txb_offset[AOM_PLANE_Y] += total_luma_tus;
    parse_ctx->first_txb_offset[AOM_PLANE_U] += 2 * total_chroma_tus;
}

static INLINE void set_txfm_ctxs(ParseCtxt *parse_ctx, TxSize tx_size, int n4_w, int n4_h, int skip,
                                 const PartitionInfo *pi) {
    ParseAboveNbr4x4Ctxt *above_parse_ctx = parse_ctx->parse_above_nbr4x4_ctxt;
    ParseLeftNbr4x4Ctxt * left_parse_ctx  = parse_ctx->parse_left_nbr4x4_ctxt;
    int                   mi_row          = pi->mi_row;
    int                   mi_col          = pi->mi_col;
    uint8_t               tx_wide         = tx_size_wide[tx_size];
    uint8_t               tx_high         = tx_size_high[tx_size];
    uint8_t *const        above_ctx =
        above_parse_ctx->above_tx_wd + (mi_col - parse_ctx->cur_tile_info.mi_col_start);
    uint8_t *const left_ctx = left_parse_ctx->left_tx_ht + (mi_row - parse_ctx->sb_row_mi);
    if (skip) {
        tx_wide = n4_w * MI_SIZE;
        tx_high = n4_h * MI_SIZE;
    }
    memset(above_ctx, tx_wide, n4_w);
    memset(left_ctx, tx_high, n4_h);
}

void read_block_tx_size(ParseCtxt *parse_ctx, PartitionInfo *part_info, BlockSize bsize) {
    BlockModeInfo *mbmi           = part_info->mi;
    SBInfo *       sb_info        = part_info->sb_info;
    int            inter_block_tx = is_inter_block(mbmi);

    if (parse_ctx->frame_header->tx_mode == TX_MODE_SELECT && bsize > BLOCK_4X4 && !mbmi->skip &&
        inter_block_tx && !parse_ctx->frame_header->lossless_array[mbmi->segment_id]) {
        const TxSize max_tx_size     = max_txsize_rect_lookup[bsize];
        const int    bh              = tx_size_high_unit[max_tx_size];
        const int    bw              = tx_size_wide_unit[max_tx_size];
        const int    width           = block_size_wide[bsize] >> tx_size_wide_log2[0];
        const int    height          = block_size_high[bsize] >> tx_size_high_log2[0];
        int          force_split_cnt = 0, num_luma_tus = 0;

        // Current luma trans_info and offset initialization
        parse_ctx->cur_luma_trans_info =
            sb_info->sb_trans_info[AOM_PLANE_Y] + mbmi->first_txb_offset[AOM_PLANE_Y];
        parse_ctx->cur_blk_luma_count = 0;

        // Luma trans_info update
        for (int idy = 0; idy < height; idy += bh)
            for (int idx = 0; idx < width; idx += bw) {
                num_luma_tus = 0;
                read_var_tx_size(parse_ctx, part_info, max_tx_size, idy, idx, 0, &num_luma_tus);
                parse_ctx->num_tus[AOM_PLANE_Y][force_split_cnt] = num_luma_tus;
                force_split_cnt++;
            }

        // Chroma trans_info update
        update_chroma_trans_info(parse_ctx, part_info, bsize);
        mbmi->num_tus[AOM_PLANE_Y] = parse_ctx->cur_blk_luma_count;
        parse_ctx->first_txb_offset[AOM_PLANE_Y] += parse_ctx->cur_blk_luma_count;
    } else {
        TxSize tx_size = read_tx_size(parse_ctx, part_info, !mbmi->skip || !inter_block_tx);

        int b4_w = mi_size_wide[mbmi->sb_type];
        int b4_h = mi_size_high[mbmi->sb_type];

        set_txfm_ctxs(
            parse_ctx, tx_size, b4_w, b4_h, mbmi->skip && is_inter_block(mbmi), part_info);

        /* Update Flat Transform Info */
        update_flat_trans_info(parse_ctx, part_info, bsize, tx_size);
    }
}

void parse_transform_type(ParseCtxt *parse_ctxt, PartitionInfo *xd, TxSize tx_size,
                          TransformInfo_t *trans_info) {
    BlockModeInfo *mbmi    = xd->mi;
    SvtReader *    r       = &parse_ctxt->r;
    TxType *       tx_type = &trans_info->txk_type;
    *tx_type               = DCT_DCT;

    FRAME_CONTEXT *frm_ctx = &parse_ctxt->cur_tile_ctx;

    // No need to read transform type if block is skipped.
    if (mbmi->skip ||
        seg_feature_active(
            &parse_ctxt->frame_header->segmentation_params, mbmi->segment_id, SEG_LVL_SKIP))
        return;
    const int qindex = parse_ctxt->frame_header->quantization_params.
                                            qindex[mbmi->segment_id];
    if (qindex == 0) return;

    const int       inter_block = is_inter_block(mbmi);
    const TxSetType tx_set_type =
        get_ext_tx_set_type(tx_size, inter_block, parse_ctxt->frame_header->reduced_tx_set);
    if (av1_num_ext_tx_set[tx_set_type] > 1) {
        const int eset =
            get_ext_tx_set(tx_size, inter_block, parse_ctxt->frame_header->reduced_tx_set);
        // eset == 0 should correspond to a set with only DCT_DCT and
        // there is no need to read the tx_type
        assert(eset != 0);

        const TxSize square_tx_size = txsize_sqr_map[tx_size];
        if (inter_block) {
            *tx_type = av1_ext_tx_inv[tx_set_type][svt_read_symbol(
                r,
                frm_ctx->inter_ext_tx_cdf[eset][square_tx_size],
                av1_num_ext_tx_set[tx_set_type],
                ACCT_STR)];
        } else {
            const PredictionMode intra_mode =
                mbmi->filter_intra_mode_info.use_filter_intra
                    ? fimode_to_intradir[mbmi->filter_intra_mode_info.filter_intra_mode]
                    : mbmi->mode;
            *tx_type = av1_ext_tx_inv[tx_set_type][svt_read_symbol(
                r,
                frm_ctx->intra_ext_tx_cdf[eset][square_tx_size][intra_mode],
                av1_num_ext_tx_set[tx_set_type],
                ACCT_STR)];
        }
    }
}

TxType compute_tx_type(PlaneType plane_type, const PartitionInfo *xd, TxSize tx_size,
                       int reduced_tx_set, uint8_t *lossless_array, TransformInfo_t *trans_info) {
    const BlockModeInfo *const mbmi = xd->mi;
    const TxSetType            tx_set_type =
        get_ext_tx_set_type(tx_size, is_inter_block(mbmi), reduced_tx_set);

    TxType tx_type = DCT_DCT;
    if (lossless_array[mbmi->segment_id] || txsize_sqr_up_map[tx_size] > TX_32X32)
        tx_type = DCT_DCT;
    else {
        if (plane_type == PLANE_TYPE_Y || is_inter_block(mbmi)) {
            tx_type = trans_info->txk_type;
        } else {
            // In intra mode, uv planes don't share the same prediction mode as y
            // plane, so the tx_type should not be shared
            tx_type = intra_mode_to_tx_type(mbmi, PLANE_TYPE_UV);
        }
    }
    assert(tx_type < TX_TYPES);
    if (!av1_ext_tx_used[tx_set_type][tx_type]) return DCT_DCT;
    return tx_type;
}

void reset_skip_context(ParseCtxt *parse_ctxt, PartitionInfo *pi, BlockSize bsize,
                        const int num_planes) {
    int i, nplanes;
    nplanes = 1 + (num_planes - 1) * pi->is_chroma_ref;

    for (i = 0; i < nplanes; i++) {
        uint8_t         suby        = i ? parse_ctxt->seq_header->color_config.subsampling_y : 0;
        uint8_t         subx        = i ? parse_ctxt->seq_header->color_config.subsampling_x : 0;
        const BlockSize plane_bsize = get_plane_block_size(bsize, subx, suby);
        const int       txs_wide    = block_size_wide[plane_bsize] >> tx_size_wide_log2[0];
        const int       txs_high    = block_size_high[plane_bsize] >> tx_size_high_log2[0];

        uint8_t *const above_ctx = parse_ctxt->parse_above_nbr4x4_ctxt->above_ctx[i] +
                                   ((pi->mi_col - parse_ctxt->cur_tile_info.mi_col_start) >> subx);
        uint8_t *const left_ctx = parse_ctxt->parse_left_nbr4x4_ctxt->left_ctx[i] +
                                  ((pi->mi_row - parse_ctxt->sb_row_mi) >> suby);

        memset(above_ctx, 0, txs_wide);
        memset(left_ctx, 0, txs_high);
    }
}

void update_coeff_ctx(ParseCtxt *parse_ctxt, int plane, PartitionInfo *pi, TxSize tx_size,
                      uint32_t blk_row, uint32_t blk_col, int above_off, int left_off,
                      int cul_level) {
    ParseAboveNbr4x4Ctxt *above_parse_ctx = parse_ctxt->parse_above_nbr4x4_ctxt;
    ParseLeftNbr4x4Ctxt * left_parse_ctx  = parse_ctxt->parse_left_nbr4x4_ctxt;

    uint8_t suby = plane ? parse_ctxt->seq_header->color_config.subsampling_y : 0;
    uint8_t subx = plane ? parse_ctxt->seq_header->color_config.subsampling_x : 0;

    uint8_t *const above_ctx = above_parse_ctx->above_ctx[plane] + blk_col -
                               (parse_ctxt->cur_tile_info.mi_col_start >> subx);
    uint8_t *const left_ctx =
        left_parse_ctx->left_ctx[plane] + (blk_row - (parse_ctxt->sb_row_mi >> suby));

    const int txs_wide = tx_size_wide_unit[tx_size];
    const int txs_high = tx_size_high_unit[tx_size];

    if (pi->mb_to_right_edge < 0) {
        int plane_bsize = (pi->mi->sb_type == BLOCK_INVALID)
                              ? BLOCK_INVALID
                              : ss_size_lookup[pi->mi->sb_type][subx][suby];
        const int blocks_wide    = max_block_wide(pi, plane_bsize, subx);
        const int above_contexts = AOMMIN(txs_wide, (blocks_wide - above_off));

        memset(above_ctx, cul_level, above_contexts);
        memset(above_ctx + above_contexts, 0, (txs_wide - above_contexts));
    } else {
        memset(above_ctx, cul_level, txs_wide);
    }

    if (pi->mb_to_bottom_edge < 0) {
        int plane_bsize = (pi->mi->sb_type == BLOCK_INVALID)
                              ? BLOCK_INVALID
                              : ss_size_lookup[pi->mi->sb_type][subx][suby];
        const int blocks_high   = max_block_high(pi, plane_bsize, suby);
        const int left_contexts = AOMMIN(txs_high, (blocks_high - left_off));

        memset(left_ctx, cul_level, left_contexts);
        memset(left_ctx + left_contexts, 0, txs_high - left_contexts);
    } else {
        memset(left_ctx, cul_level, txs_high);
    }
}

static INLINE int rec_eob_pos(const int eob_token, const int extra) {
    int eob = eb_k_eob_group_start[eob_token];
    if (eob > 2) eob += extra;
    return eob;
}

static INLINE int get_lower_levels_ctx_2d(const uint8_t *levels, int coeff_idx, int bwl,
                                          TxSize tx_size) {
    assert(coeff_idx > 0);
    int mag;
    levels = levels + get_padded_idx(coeff_idx, bwl);
    mag    = AOMMIN(levels[1], 3); // { 0, 1 }
    mag += AOMMIN(levels[(1 << bwl) + TX_PAD_HOR], 3); // { 1, 0 }
    mag += AOMMIN(levels[(1 << bwl) + TX_PAD_HOR + 1], 3); // { 1, 1 }
    mag += AOMMIN(levels[2], 3); // { 0, 2 }
    mag += AOMMIN(levels[(2 << bwl) + (2 << TX_PAD_HOR_LOG2)], 3); // { 2, 0 }

    const int ctx = AOMMIN((mag + 1) >> 1, 4);
    return ctx + eb_av1_nz_map_ctx_offset[tx_size][coeff_idx];
}

static INLINE int read_golomb(SvtReader *r) {
    int x      = 1;
    int length = 0;
    int i      = 0;

    while (!i) {
        i = svt_read_bit(r, ACCT_STR);
        ++length;
        if (length > 20) {
            SVT_LOG("Invalid length in read_golomb");
            break;
        }
    }

    for (i = 0; i < length - 1; ++i) {
        x <<= 1;
        x += svt_read_bit(r, ACCT_STR);
    }

    return x - 1;
}

static INLINE int get_br_ctx_2d(const uint8_t *const levels,
                                const int            c, // raster order
                                const int            bwl) {
    assert(c > 0);
    const int row    = c >> bwl;
    const int col    = c - (row << bwl);
    const int stride = (1 << bwl) + TX_PAD_HOR;
    const int pos    = row * stride + col;
    int       mag    = AOMMIN(levels[pos + 1], MAX_BASE_BR_RANGE) +
              AOMMIN(levels[pos + stride], MAX_BASE_BR_RANGE) +
              AOMMIN(levels[pos + 1 + stride], MAX_BASE_BR_RANGE);
    mag = AOMMIN((mag + 1) >> 1, 6);
    if ((row | col) < 2) return mag + 7;
    return mag + 14;
}

static INLINE void read_coeffs_reverse_2d(SvtReader *r, TxSize tx_size, int start_si, int end_si,
                                          const int16_t *scan, int bwl, uint8_t *levels,
                                          BaseCdfArr base_cdf, BrCdfArr br_cdf) {
    for (int c = end_si; c >= start_si; --c) {
        const int pos       = scan[c];
        const int coeff_ctx = get_lower_levels_ctx_2d(levels, pos, bwl, tx_size);
        const int nsymbs    = 4;
        int       level     = svt_read_symbol(r, base_cdf[coeff_ctx], nsymbs, ACCT_STR);
        if (level > NUM_BASE_LEVELS) {
            const int   br_ctx = get_br_ctx_2d(levels, pos, bwl);
            AomCdfProb *cdf    = br_cdf[br_ctx];
            for (int idx = 0; idx < COEFF_BASE_RANGE; idx += BR_CDF_SIZE - 1) {
                const int k = svt_read_symbol(r, cdf, BR_CDF_SIZE, ACCT_STR);
                level += k;
                if (k < BR_CDF_SIZE - 1) break;
            }
        }
        levels[get_padded_idx(pos, bwl)] = level;
    }
}

static INLINE void read_coeffs_reverse(SvtReader *r, TxSize tx_size, TxType tx_type, int start_si,
                                       int end_si, const int16_t *scan, int bwl, uint8_t *levels,
                                       BaseCdfArr base_cdf, BrCdfArr br_cdf) {
    TxClass tx_class = tx_type_to_class[tx_type];
    for (int c = end_si; c >= start_si; --c) {
        const int pos       = scan[c];
        const int coeff_ctx = get_lower_levels_ctx(levels, pos, bwl, tx_size, tx_class);
        const int nsymbs    = 4;
        int       level     = svt_read_symbol(r, base_cdf[coeff_ctx], nsymbs, ACCT_STR);
        if (level > NUM_BASE_LEVELS) {
            const int br_ctx = get_br_ctx(levels, pos, bwl, tx_class);
            AomCdfProb *cdf    = br_cdf[br_ctx];
            for (int idx = 0; idx < COEFF_BASE_RANGE; idx += BR_CDF_SIZE - 1) {
                const int k = svt_read_symbol(r, cdf, BR_CDF_SIZE, ACCT_STR);
                level += k;
                if (k < BR_CDF_SIZE - 1) break;
            }
        }
        levels[get_padded_idx(pos, bwl)] = level;
    }
}

static INLINE void set_dc_sign(int *cul_level, int dc_val) {
    if (dc_val < 0)
        *cul_level |= 1 << COEFF_CONTEXT_BITS;
    else if (dc_val > 0)
        *cul_level += 2 << COEFF_CONTEXT_BITS;
}

uint16_t parse_coeffs(ParseCtxt *parse_ctxt, PartitionInfo *xd, uint32_t blk_row, uint32_t blk_col,
                      int above_off, int left_off, int plane, int txb_skip_ctx, int dc_sign_ctx,
                      TxSize tx_size, int32_t *coeff_buf, TransformInfo_t *trans_info) {
    SvtReader *r      = &parse_ctxt->r;
    const int  width  = get_txb_wide(tx_size);
    const int  height = get_txb_high(tx_size);

    FRAME_CONTEXT *frm_ctx = &parse_ctxt->cur_tile_ctx;

    TxSize    txs_ctx = (TxSize)((txsize_sqr_map[tx_size] + txsize_sqr_up_map[tx_size] + 1) >> 1);
    PlaneType plane_type = (plane == 0) ? PLANE_TYPE_Y : PLANE_TYPE_UV;
    int       cul_level  = 0;
    int       dc_val     = 0;
    uint8_t   levels_buf[TX_PAD_2D];
    uint8_t *const levels = set_levels(levels_buf, width);

    const int all_zero =
        svt_read_symbol(r, frm_ctx->txb_skip_cdf[txs_ctx][txb_skip_ctx], 2, ACCT_STR);

    const int bwl = get_txb_bwl(tx_size);

    uint16_t eob           = 0;
    uint16_t max_scan_line = 0;
    if (all_zero) {
        if (plane == 0) {
            trans_info->txk_type = DCT_DCT;
            trans_info->cbf      = 0;
        }

        update_coeff_ctx(
            parse_ctxt, plane, xd, tx_size, blk_row, blk_col, above_off, left_off, cul_level);

        return 0;
    }

    if (plane == AOM_PLANE_Y) parse_transform_type(parse_ctxt, xd, tx_size, trans_info);

    uint8_t *        lossless_array = &parse_ctxt->frame_header->lossless_array[0];
    TransformInfo_t *trans_buf =
        (is_inter_block(xd->mi) && plane) ? parse_ctxt->inter_trans_chroma : trans_info;
    trans_info->txk_type                = compute_tx_type(plane_type,
                                           xd,
                                           tx_size,
                                           parse_ctxt->frame_header->reduced_tx_set,
                                           lossless_array,
                                           trans_buf);
    const ScanOrder *    scan_order     = &av1_scan_orders[tx_size][trans_info->txk_type];
    const int16_t *const scan           = scan_order->scan;
    const int            eob_multi_size = txsize_log2_minus4[tx_size];

    const TxClass tx_class      = tx_type_to_class[trans_info->txk_type];
    const int     eob_multi_ctx = (tx_class == TX_CLASS_2D) ? 0 : 1;

    int eob_extra = 0;
    int eob_pt    = 1;

    switch (eob_multi_size) {
    case 0:
        eob_pt =
            svt_read_symbol(r, frm_ctx->eob_flag_cdf16[plane_type][eob_multi_ctx], 5, ACCT_STR) + 1;
        break;
    case 1:
        eob_pt =
            svt_read_symbol(r, frm_ctx->eob_flag_cdf32[plane_type][eob_multi_ctx], 6, ACCT_STR) + 1;
        break;
    case 2:
        eob_pt =
            svt_read_symbol(r, frm_ctx->eob_flag_cdf64[plane_type][eob_multi_ctx], 7, ACCT_STR) + 1;
        break;
    case 3:
        eob_pt =
            svt_read_symbol(r, frm_ctx->eob_flag_cdf128[plane_type][eob_multi_ctx], 8, ACCT_STR) +
            1;
        break;
    case 4:
        eob_pt =
            svt_read_symbol(r, frm_ctx->eob_flag_cdf256[plane_type][eob_multi_ctx], 9, ACCT_STR) +
            1;
        break;
    case 5:
        eob_pt =
            svt_read_symbol(r, frm_ctx->eob_flag_cdf512[plane_type][eob_multi_ctx], 10, ACCT_STR) +
            1;
        break;
    default:
        eob_pt =
            svt_read_symbol(r, frm_ctx->eob_flag_cdf1024[plane_type][eob_multi_ctx], 11, ACCT_STR) +
            1;
        break;
    }

    int eob_shift = eb_k_eob_offset_bits[eob_pt];
    if (eob_shift > 0) {
        const int eob_ctx = eob_pt;
        int       bit =
            svt_read_symbol(r, frm_ctx->eob_extra_cdf[txs_ctx][plane_type][eob_ctx], 2, ACCT_STR);
        if (bit) eob_extra += (1 << (eob_shift - 1));
        for (int i = 1; i < eob_shift; i++) {
            bit = svt_read_bit(r, ACCT_STR);
            if (bit) eob_extra += (1 << (eob_shift - 1 - i));
        }
    }
    eob = rec_eob_pos(eob_pt, eob_extra);

    if (eob > 1) {
        memset(levels_buf,
               0,
               sizeof(*levels_buf) * ((width + TX_PAD_HOR) * (height + TX_PAD_VER) + TX_PAD_END));
    }

    int         i         = eob - 1;
    const int   pos       = scan[i];
    const int   coeff_ctx = get_lower_levels_ctx_eob(bwl, height, i);
    const int   nsymbs    = 3;
    AomCdfProb *cdf       = frm_ctx->coeff_base_eob_cdf[txs_ctx][plane_type][coeff_ctx];
    int         level     = svt_read_symbol(r, cdf, nsymbs, ACCT_STR) + 1;
    if (level > NUM_BASE_LEVELS) {
        const int br_ctx = get_br_ctx_eob(pos, bwl, tx_class);
        cdf              = frm_ctx->coeff_br_cdf[AOMMIN(txs_ctx, TX_32X32)][plane_type][br_ctx];
        for (int idx = 0; idx < COEFF_BASE_RANGE / (BR_CDF_SIZE - 1); idx++) {
            int coeff_br = svt_read_symbol(r, cdf, BR_CDF_SIZE, ACCT_STR);
            level += coeff_br;
            if (coeff_br < BR_CDF_SIZE - 1) break;
        }
    }
    levels[get_padded_idx(pos, bwl)] = level;

    if (eob > 1) {
        BaseCdfArr base_cdf = frm_ctx->coeff_base_cdf[txs_ctx][plane_type];
        BrCdfArr   br_cdf   = frm_ctx->coeff_br_cdf[AOMMIN(txs_ctx, TX_32X32)][plane_type];
        if (tx_class == TX_CLASS_2D) {
            read_coeffs_reverse_2d(r, tx_size, 1, eob - 1 - 1, scan, bwl, levels, base_cdf, br_cdf);
            read_coeffs_reverse(
                r, tx_size, trans_info->txk_type, 0, 0, scan, bwl, levels, base_cdf, br_cdf);
        } else {
            read_coeffs_reverse(r,
                                tx_size,
                                trans_info->txk_type,
                                0,
                                eob - 1 - 1,
                                scan,
                                bwl,
                                levels,
                                base_cdf,
                                br_cdf);
        }
    }
#if SVT_DEC_COEFF_DEBUG
    {
        int16_t *cur_coeff = (int16_t *)coeff_buf;
        /* cur_coeff[0] used for debug */
        cur_coeff[1] = eob;
    }
#else
    coeff_buf[0] = eob;
#endif

    for (int c = 0; c < eob; ++c) {
        const int pos   = scan[c];
        uint8_t   sign  = 0;
        TranLow   level = levels[get_padded_idx(pos, bwl)];
        if (level) {
            max_scan_line = AOMMAX(max_scan_line, pos);
            if (c == 0) {
                sign =
                    svt_read_symbol(r, frm_ctx->dc_sign_cdf[plane_type][dc_sign_ctx], 2, ACCT_STR);
            } else
                sign = svt_read_bit(r, ACCT_STR);
            if (level >= MAX_BASE_BR_RANGE) level += read_golomb(r);
            if (c == 0) dc_val = sign ? -level : level;

            level &= 0xfffff;
            cul_level += level;
        }
        coeff_buf[c + 1] = sign ? -level : level;
    }

    cul_level = AOMMIN(COEFF_CONTEXT_MASK, cul_level);
    set_dc_sign(&cul_level, dc_val);

    update_coeff_ctx(
        parse_ctxt, plane, xd, tx_size, blk_row, blk_col, above_off, left_off, cul_level);

    trans_info->cbf = 1;
    assert(eob);

    return eob;
}

static INLINE int partition_plane_context(int mi_row, int mi_col, BlockSize bsize,
                                          ParseCtxt *parse_ctxt) {
    const uint8_t *above_ctx = parse_ctxt->parse_above_nbr4x4_ctxt->above_part_wd + mi_col -
                               parse_ctxt->cur_tile_info.mi_col_start;
    const uint8_t *left_ctx = parse_ctxt->parse_left_nbr4x4_ctxt->left_part_ht +
                              ((mi_row - parse_ctxt->sb_row_mi) & MAX_MIB_MASK);

    // Minimum partition point is 8x8. Offset the bsl accordingly.
    int bsl   = mi_size_wide_log2[bsize] - mi_size_wide_log2[BLOCK_8X8];
    int above = (*above_ctx >> bsl) & 1, left = (*left_ctx >> bsl) & 1;

    assert(mi_size_wide_log2[bsize] == mi_size_high_log2[bsize]);
    assert(bsl >= 0);

    return (left * 2 + above) + bsl * PARTITION_PLOFFSET;
}

PartitionType parse_partition_type(uint32_t blk_row, uint32_t blk_col, BlockSize bsize,
                                   int has_rows, int has_cols, ParseCtxt *parse_ctxt) {
    SvtReader *reader = &parse_ctxt->r;
    int        partition_cdf_length =
        bsize <= BLOCK_8X8
            ? PARTITION_TYPES
            : (bsize == BLOCK_128X128 ? EXT_PARTITION_TYPES - 2 : EXT_PARTITION_TYPES);
    int ctx = partition_plane_context(blk_row, blk_col, bsize, parse_ctxt);

    if (bsize < BLOCK_8X8)
        return PARTITION_NONE;
    else if (has_rows && has_cols) {
        return (PartitionType)svt_read_symbol(
            reader, parse_ctxt->cur_tile_ctx.partition_cdf[ctx], partition_cdf_length, ACCT_STR);
    } else if (has_cols) {
        assert(bsize > BLOCK_8X8);
        AomCdfProb cdf[3];
        partition_gather_vert_alike(cdf, parse_ctxt->cur_tile_ctx.partition_cdf[ctx], bsize);
        assert(cdf[1] == AOM_ICDF(CDF_PROB_TOP));
        return svt_read_cdf(reader, cdf, 2, ACCT_STR) ? PARTITION_SPLIT : PARTITION_HORZ;
    } else if (has_rows) {
        assert(has_rows && !has_cols);
        assert(bsize > BLOCK_8X8);
        AomCdfProb cdf[3];
        partition_gather_horz_alike(cdf, parse_ctxt->cur_tile_ctx.partition_cdf[ctx], bsize);
        assert(cdf[1] == AOM_ICDF(CDF_PROB_TOP));
        return svt_read_cdf(reader, cdf, 2, ACCT_STR) ? PARTITION_SPLIT : PARTITION_VERT;
    }
    return PARTITION_SPLIT;
}

static INLINE int combine_entropy_contexts(int8_t a, int8_t b) { return (a != 0) + (b != 0); }

static INLINE int get_entropy_context(TxSize tx_size, const uint8_t *a, const uint8_t *l) {
    uint8_t above_ec = 0, left_ec = 0;

    switch (tx_size) {
    case TX_4X4:
        above_ec = a[0] != 0;
        left_ec  = l[0] != 0;
        break;
    case TX_4X8:
        above_ec = a[0] != 0;
        left_ec  = !!*(const uint16_t *)l;
        break;
    case TX_8X4:
        above_ec = !!*(const uint16_t *)a;
        left_ec  = l[0] != 0;
        break;
    case TX_8X16:
        above_ec = !!*(const uint16_t *)a;
        left_ec  = !!*(const uint32_t *)l;
        break;
    case TX_16X8:
        above_ec = !!*(const uint32_t *)a;
        left_ec  = !!*(const uint16_t *)l;
        break;
    case TX_16X32:
        above_ec = !!*(const uint32_t *)a;
        left_ec  = !!*(const uint64_t *)l;
        break;
    case TX_32X16:
        above_ec = !!*(const uint64_t *)a;
        left_ec  = !!*(const uint32_t *)l;
        break;
    case TX_8X8:
        above_ec = !!*(const uint16_t *)a;
        left_ec  = !!*(const uint16_t *)l;
        break;
    case TX_16X16:
        above_ec = !!*(const uint32_t *)a;
        left_ec  = !!*(const uint32_t *)l;
        break;
    case TX_32X32:
        above_ec = !!*(const uint64_t *)a;
        left_ec  = !!*(const uint64_t *)l;
        break;
    case TX_64X64:
        above_ec = !!(*(const uint64_t *)a | *(const uint64_t *)(a + 8));
        left_ec  = !!(*(const uint64_t *)l | *(const uint64_t *)(l + 8));
        break;
    case TX_32X64:
        above_ec = !!*(const uint64_t *)a;
        left_ec  = !!(*(const uint64_t *)l | *(const uint64_t *)(l + 8));
        break;
    case TX_64X32:
        above_ec = !!(*(const uint64_t *)a | *(const uint64_t *)(a + 8));
        left_ec  = !!*(const uint64_t *)l;
        break;
    case TX_4X16:
        above_ec = a[0] != 0;
        left_ec  = !!*(const uint32_t *)l;
        break;
    case TX_16X4:
        above_ec = !!*(const uint32_t *)a;
        left_ec  = l[0] != 0;
        break;
    case TX_8X32:
        above_ec = !!*(const uint16_t *)a;
        left_ec  = !!*(const uint64_t *)l;
        break;
    case TX_32X8:
        above_ec = !!*(const uint64_t *)a;
        left_ec  = !!*(const uint16_t *)l;
        break;
    case TX_16X64:
        above_ec = !!*(const uint32_t *)a;
        left_ec  = !!(*(const uint64_t *)l | *(const uint64_t *)(l + 8));
        break;
    case TX_64X16:
        above_ec = !!(*(const uint64_t *)a | *(const uint64_t *)(a + 8));
        left_ec  = !!*(const uint32_t *)l;
        break;
    default: assert(0 && "Invalid transform size."); break;
    }
    return combine_entropy_contexts(above_ec, left_ec);
}

static INLINE void dec_get_txb_ctx(ParseCtxt *parse_ctx, const TxSize tx_size, const int plane,
                                   int plane_bsize, int txb_h_unit, int txb_w_unit, int blk_row,
                                   int blk_col, TXB_CTX *const txb_ctx) {
#define MAX_TX_SIZE_UNIT 16

    ParseAboveNbr4x4Ctxt *above_parse_ctx = parse_ctx->parse_above_nbr4x4_ctxt;
    ParseLeftNbr4x4Ctxt * left_parse_ctx  = parse_ctx->parse_left_nbr4x4_ctxt;
    EbColorConfig *       clr_cfg         = &parse_ctx->seq_header->color_config;

    int suby    = plane ? clr_cfg->subsampling_y : 0;
    int subx    = plane ? clr_cfg->subsampling_x : 0;
    int dc_sign = 0;
    int k       = 0;

    static const int8_t signs[3]                                   = {0, -1, 1};
    static const int8_t dc_sign_contexts[4 * MAX_TX_SIZE_UNIT + 1] = {
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
    uint8_t *above_ctx = above_parse_ctx->above_ctx[plane] + blk_col -
                         (parse_ctx->cur_tile_info.mi_col_start >> subx);
    uint8_t *left_ctx =
        left_parse_ctx->left_ctx[plane] + (blk_row - (parse_ctx->sb_row_mi >> suby));

    do {
        const unsigned int sign = ((uint8_t)above_ctx[k] >> COEFF_CONTEXT_BITS);
        assert(sign <= 2);
        dc_sign += signs[sign];
    } while (++k < txb_w_unit);

    k = 0;
    do {
        const unsigned int sign = ((uint8_t)left_ctx[k] >> COEFF_CONTEXT_BITS);
        assert(sign <= 2);
        dc_sign += signs[sign];
    } while (++k < txb_h_unit);

    txb_ctx->dc_sign_ctx = dc_sign_contexts[dc_sign + 2 * MAX_TX_SIZE_UNIT];

    if (plane == 0) {
        if (plane_bsize == txsize_to_bsize[tx_size])
            txb_ctx->txb_skip_ctx = 0;
        else {
            static const uint8_t skip_contexts[5][5] = {{1, 2, 2, 2, 3},
                                                        {1, 4, 4, 4, 5},
                                                        {1, 4, 4, 4, 5},
                                                        {1, 4, 4, 4, 5},
                                                        {1, 4, 4, 4, 6}};
            int                  top                 = 0;
            int                  left                = 0;

            k = 0;
            do { top |= above_ctx[k]; } while (++k < txb_w_unit);
            top &= COEFF_CONTEXT_MASK;

            k = 0;
            do { left |= left_ctx[k]; } while (++k < txb_h_unit);
            left &= COEFF_CONTEXT_MASK;

            const int max = AOMMIN(top | left, 4);
            const int min = AOMMIN(AOMMIN(top, left), 4);

            txb_ctx->txb_skip_ctx = skip_contexts[min][max];
        }
    } else {
        const int ctx_base = get_entropy_context(tx_size, above_ctx, left_ctx);
        const int ctx_offset =
            (num_pels_log2_lookup[plane_bsize] > num_pels_log2_lookup[txsize_to_bsize[tx_size]])
                ? 10
                : 7;
        txb_ctx->txb_skip_ctx = ctx_base + ctx_offset;
    }
#undef MAX_TX_SIZE_UNIT
}

uint16_t parse_transform_block(ParseCtxt *parse_ctx, PartitionInfo *pi, int32_t *coeff,
                               TransformInfo_t *trans_info, int plane, int blk_col, int blk_row,
                               int start_x, int start_y, TxSize tx_size, int sub_x, int sub_y) {
    uint16_t eob = 0;
    TXB_CTX  txb_ctx;

    BlockSize bsize = pi->mi->sb_type;
    int       plane_bsize =
        (bsize == BLOCK_INVALID) ? BLOCK_INVALID : ss_size_lookup[bsize][sub_x][sub_y];

    int txb_w_unit = tx_size_wide_unit[tx_size];
    int txb_h_unit = tx_size_high_unit[tx_size];

    if (pi->mb_to_right_edge < 0) {
        int plane_bsize = (pi->mi->sb_type == BLOCK_INVALID)
                              ? BLOCK_INVALID
                              : ss_size_lookup[pi->mi->sb_type][sub_x][sub_y];
        const int blocks_wide = max_block_wide(pi, plane_bsize, sub_x);
        txb_w_unit            = AOMMIN(txb_w_unit, (blocks_wide - blk_col));
    }

    if (pi->mb_to_bottom_edge < 0) {
        int plane_bsize = (pi->mi->sb_type == BLOCK_INVALID)
                              ? BLOCK_INVALID
                              : ss_size_lookup[pi->mi->sb_type][sub_x][sub_y];
        const int blocks_high = max_block_high(pi, plane_bsize, sub_y);
        txb_h_unit            = AOMMIN(txb_h_unit, (blocks_high - blk_row));
    }

    dec_get_txb_ctx(
        parse_ctx, tx_size, plane, plane_bsize, txb_h_unit, txb_w_unit, start_y, start_x, &txb_ctx);

    eob = parse_coeffs(parse_ctx,
                       pi,
                       start_y,
                       start_x,
                       blk_col,
                       blk_row,
                       plane,
                       txb_ctx.txb_skip_ctx,
                       txb_ctx.dc_sign_ctx,
                       tx_size,
                       coeff,
                       trans_info);

    return eob;
}

void parse_residual(ParseCtxt *parse_ctx, PartitionInfo *pi, BlockSize mi_size) {
    EbColorConfig *color_info      = &parse_ctx->seq_header->color_config;
    SBInfo *       sb_info         = pi->sb_info;
    int            num_planes      = color_info->mono_chrome ? 1 : MAX_MB_PLANE;
    BlockModeInfo *mode            = pi->mi;
    int            mi_row          = pi->mi_row;
    int            mi_col          = pi->mi_col;
    int            skip            = mode->skip;
    uint32_t       force_split_cnt = 0;
    uint32_t       num_tu, total_num_tu;

    const int       max_blocks_wide = max_block_wide(pi, mi_size, 0);
    const int       max_blocks_high = max_block_high(pi, mi_size, 0);
    const BlockSize max_unit_bsize  = BLOCK_64X64;
    int             mu_blocks_wide  = block_size_wide[max_unit_bsize] >> tx_size_wide_log2[0];
    int             mu_blocks_high  = block_size_high[max_unit_bsize] >> tx_size_high_log2[0];
    mu_blocks_wide                  = AOMMIN(max_blocks_wide, mu_blocks_wide);
    mu_blocks_high                  = AOMMIN(max_blocks_high, mu_blocks_high);

    TransformInfo_t *trans_info[MAX_MB_PLANE];
    uint8_t          lossless = parse_ctx->frame_header->lossless_array[pi->mi->segment_id];
    int lossless_block = (lossless && ((mi_size >= BLOCK_64X64) && (mi_size <= BLOCK_128X128)));

    int num_chroma_tus = lossless_block ? ((max_blocks_wide * max_blocks_high) >>
                                           (color_info->subsampling_x + color_info->subsampling_y))
                                        : mode->num_tus[AOM_PLANE_U];

    trans_info[AOM_PLANE_Y] =
        (sb_info->sb_trans_info[AOM_PLANE_Y] + mode->first_txb_offset[AOM_PLANE_Y]);
    trans_info[AOM_PLANE_U] =
        (sb_info->sb_trans_info[AOM_PLANE_U] + mode->first_txb_offset[AOM_PLANE_U]);
    trans_info[AOM_PLANE_V] =
        (sb_info->sb_trans_info[AOM_PLANE_U] + mode->first_txb_offset[AOM_PLANE_U]) +
        num_chroma_tus;

    for (int row = 0; row < max_blocks_high; row += mu_blocks_high) {
        for (int col = 0; col < max_blocks_wide; col += mu_blocks_wide) {
            for (int plane = 0; plane < num_planes; ++plane) {
                int sub_x = (plane > 0) ? color_info->subsampling_x : 0;
                int sub_y = (plane > 0) ? color_info->subsampling_y : 0;

                if (plane && !pi->is_chroma_ref) continue;

                if (is_inter_block(mode) && !plane)
                    parse_ctx->inter_trans_chroma = trans_info[plane];
                if (lossless_block) {
                    int unit_height =
                        ROUND_POWER_OF_TWO(AOMMIN(mu_blocks_high + row, max_blocks_high), 0);
                    int unit_width =
                        ROUND_POWER_OF_TWO(AOMMIN(mu_blocks_wide + col, max_blocks_wide), 0);
                    assert(trans_info[plane]->tx_size == TX_4X4);
                    num_tu = ((unit_width - col) * (unit_height - row)) >> (sub_x + sub_y);
                } else {
                    total_num_tu = mode->num_tus[!!plane];
                    num_tu       = parse_ctx->num_tus[plane][force_split_cnt];

                    assert(total_num_tu != 0);
                    assert(total_num_tu ==
                           (uint32_t)(parse_ctx->num_tus[plane][0] + parse_ctx->num_tus[plane][1] +
                                      parse_ctx->num_tus[plane][2] + parse_ctx->num_tus[plane][3]));
                }
                assert(num_tu != 0);

                (void)total_num_tu;

                for (uint32_t tu = 0; tu < num_tu; tu++) {
                    assert(trans_info[plane]->txb_x_offset <= max_blocks_wide);
                    assert(trans_info[plane]->txb_y_offset <= max_blocks_high);

                    int32_t *coeff = parse_ctx->cur_coeff_buf[plane];
#if SVT_DEC_COEFF_DEBUG
                    {
                        uint8_t *cur_coeff = (uint8_t *)coeff;
                        uint8_t  cur_loc   = (mi_row + trans_info[plane]->txb_y_offset) & 0xFF;
                        cur_coeff[0]       = cur_loc;
                        cur_loc            = (mi_col + trans_info[plane]->txb_x_offset) & 0xFF;
                        cur_coeff[1]       = cur_loc;
                    }
#endif
                    int32_t eob = 0;

                    int blk_col = trans_info[plane]->txb_x_offset;
                    int blk_row = trans_info[plane]->txb_y_offset;

                    uint32_t start_x = (mi_col >> sub_x) + blk_col;
                    uint32_t start_y = (mi_row >> sub_y) + blk_row;

                    if (start_x >= (parse_ctx->frame_header->mi_cols >> sub_x) ||
                        start_y >= (parse_ctx->frame_header->mi_rows >> sub_y))
                        return;

                    if (!skip) {
                        eob = parse_transform_block(parse_ctx,
                                                    pi,
                                                    coeff,
                                                    trans_info[plane],
                                                    plane,
                                                    blk_col,
                                                    blk_row,
                                                    start_x,
                                                    start_y,
                                                    trans_info[plane]->tx_size,
                                                    sub_x,
                                                    sub_y);
                    }

                    if (eob != 0) {
                        parse_ctx->cur_coeff_buf[plane] += (eob + 1);
                        trans_info[plane]->cbf = 1;
                    } else
                        trans_info[plane]->cbf = 0;

                    // increment transform pointer
                    trans_info[plane]++;
                } //for tu
            } //for plane
            force_split_cnt++;
        }
    }
}

void parse_block(EbDecHandle *dec_handle, ParseCtxt *parse_ctx, uint32_t mi_row, uint32_t mi_col,
                 BlockSize subsize, TileInfo *tile, SBInfo *sb_info, PartitionType partition) {
    BlockModeInfo *mode = parse_ctx->cur_mode_info;

    int bw4 = mi_size_wide[subsize];
    int bh4 = mi_size_high[subsize];

    uint32_t mi_cols = parse_ctx->frame_header->mi_cols;
    uint32_t mi_rows = parse_ctx->frame_header->mi_rows;

    int32_t sub_x         = dec_handle->seq_header.color_config.subsampling_x;
    int32_t sub_y         = dec_handle->seq_header.color_config.subsampling_y;
    int32_t is_chroma_ref = is_chroma_reference(mi_row, mi_col, subsize, sub_x, sub_y);

    /* TODO: Can move to a common init fun for parse & decode */
    PartitionInfo part_info;
    part_info.mi      = mode;
    part_info.sb_info = sb_info;
    part_info.mi_row  = mi_row;
    part_info.mi_col  = mi_col;

    part_info.cdef_strength = sb_info->sb_cdef_strength;
    part_info.is_chroma_ref = is_chroma_ref;

    mode->partition = partition;
    /* TU offset update from parse ctxt info of previous block */
    mode->first_txb_offset[AOM_PLANE_Y] = parse_ctx->first_txb_offset[AOM_PLANE_Y];
    mode->first_txb_offset[AOM_PLANE_U] = parse_ctx->first_txb_offset[AOM_PLANE_U];
#if MODE_INFO_DBG
    mode->mi_row = mi_row;
    mode->mi_col = mi_col;
#endif

    sb_info->num_block++;
    mode->mi_row_in_sb = mi_row - parse_ctx->sb_row_mi;
    mode->mi_col_in_sb = mi_col - parse_ctx->sb_col_mi;

    /* TODO : tile->tile_rows boundary condn check is wrong */
    part_info.up_available      = ((int32_t)mi_row > tile->mi_row_start);
    part_info.left_available    = ((int32_t)mi_col > tile->mi_col_start);
    part_info.mb_to_left_edge   = -(((int32_t)mi_col * MI_SIZE) * 8);
    part_info.mb_to_right_edge  = ((int32_t)mi_cols - bw4 - mi_col) * MI_SIZE * 8;
    part_info.mb_to_top_edge    = -(((int32_t)mi_row * MI_SIZE) * 8);
    part_info.mb_to_bottom_edge = ((int32_t)mi_rows - bh4 - mi_row) * MI_SIZE * 8;

    part_info.is_sec_rect = 0;
    if (bw4 < bh4) {
        if (!((mi_col + bw4) & (bh4 - 1))) part_info.is_sec_rect = 1;
    }

    if (bw4 > bh4)
        if (mi_row & (bw4 - 1)) part_info.is_sec_rect = 1;

    if (part_info.up_available)
        part_info.above_mbmi = get_top_mode_info(dec_handle, mi_row, mi_col, sb_info);
    else
        part_info.above_mbmi = NULL;
    if (part_info.left_available)
        part_info.left_mbmi = get_left_mode_info(dec_handle, mi_row, mi_col, sb_info);
    else
        part_info.left_mbmi = NULL;
    mode->sb_type = subsize;

    mode_info(dec_handle, &part_info, parse_ctx);

    /* Initialize block or force splt block tu count to 0*/
    ZERO_ARRAY(parse_ctx->num_tus[AOM_PLANE_Y], 4);
    ZERO_ARRAY(parse_ctx->num_tus[AOM_PLANE_U], 4);
    ZERO_ARRAY(parse_ctx->num_tus[AOM_PLANE_V], 4);

    if (!is_inter_block(mode)) palette_tokens(dec_handle, parse_ctx, &part_info);

    read_block_tx_size(parse_ctx, &part_info, subsize);

    if (mode->skip) {
        int num_planes = dec_handle->seq_header.color_config.mono_chrome ? 1 : MAX_MB_PLANE;
        reset_skip_context(parse_ctx, &part_info, mode->sb_type, num_planes);
    }
    parse_residual(parse_ctx, &part_info, subsize);

    /* Update block level MI map */
    update_block_nbrs(dec_handle, parse_ctx, mi_row, mi_col, subsize);
    parse_ctx->cur_mode_info_cnt++;
    parse_ctx->cur_mode_info++;
}

static INLINE void update_partition_context(ParseCtxt *parse_ctx, int mi_row, int mi_col,
                                            BlockSize subsize, BlockSize bsize) {
    ParseAboveNbr4x4Ctxt *above_parse_ctx = parse_ctx->parse_above_nbr4x4_ctxt;
    ParseLeftNbr4x4Ctxt * left_parse_ctx  = parse_ctx->parse_left_nbr4x4_ctxt;
    uint8_t *const        above_ctx =
        above_parse_ctx->above_part_wd + mi_col - parse_ctx->cur_tile_info.mi_col_start;
    uint8_t *const left_ctx =
        left_parse_ctx->left_part_ht + ((mi_row - parse_ctx->sb_row_mi) & MAX_MIB_MASK);

    const int bw = mi_size_wide[bsize];
    const int bh = mi_size_high[bsize];
    memset(above_ctx, partition_context_lookup[subsize].above, bw);
    memset(left_ctx, partition_context_lookup[subsize].left, bh);
}

static INLINE void update_ext_partition_context(ParseCtxt *parse_ctx, int mi_row, int mi_col,
                                                BlockSize subsize, BlockSize bsize,
                                                PartitionType partition) {
    if (bsize >= BLOCK_8X8) {
        const int hbs    = mi_size_wide[bsize] / 2;
        BlockSize bsize2 = partition_subsize[PARTITION_SPLIT][bsize];
        switch (partition) {
        case PARTITION_SPLIT:
            if (bsize != BLOCK_8X8) break;
            goto PARTITIONS;
        case PARTITION_NONE:
        case PARTITION_HORZ:
        case PARTITION_VERT:
        case PARTITION_HORZ_4:
        case PARTITION_VERT_4:
        PARTITIONS:
            update_partition_context(parse_ctx, mi_row, mi_col, subsize, bsize);
            break;
        case PARTITION_HORZ_A:
            update_partition_context(parse_ctx, mi_row, mi_col, bsize2, subsize);
            update_partition_context(parse_ctx, mi_row + hbs, mi_col, subsize, subsize);
            break;
        case PARTITION_HORZ_B:
            update_partition_context(parse_ctx, mi_row, mi_col, subsize, subsize);
            update_partition_context(parse_ctx, mi_row + hbs, mi_col, bsize2, subsize);
            break;
        case PARTITION_VERT_A:
            update_partition_context(parse_ctx, mi_row, mi_col, bsize2, subsize);
            update_partition_context(parse_ctx, mi_row, mi_col + hbs, subsize, subsize);
            break;
        case PARTITION_VERT_B:
            update_partition_context(parse_ctx, mi_row, mi_col, subsize, subsize);
            update_partition_context(parse_ctx, mi_row, mi_col + hbs, bsize2, subsize);
            break;
        default: assert(0 && "Invalid partition type");
        }
    }
}

void parse_partition(EbDecHandle *dec_handle, ParseCtxt *parse_ctx, uint32_t blk_row,
                     uint32_t blk_col, BlockSize bsize, SBInfo *sb_info) {
    if (blk_row >= parse_ctx->frame_header->mi_rows || blk_col >= parse_ctx->frame_header->mi_cols)
        return;

    int num4x4            = mi_size_wide[bsize];
    int half_block_4x4    = num4x4 >> 1;
    int quarter_block_4x4 = half_block_4x4 >> 1;

    int has_rows = (blk_row + half_block_4x4) < parse_ctx->frame_header->mi_rows;
    int has_cols = (blk_col + half_block_4x4) < parse_ctx->frame_header->mi_cols;

    PartitionType partition;

    partition = (bsize < BLOCK_8X8)
                    ? PARTITION_NONE
                    : parse_partition_type(blk_row, blk_col, bsize, has_rows, has_cols, parse_ctx);
    int sub_size   = partition_subsize[(int)partition][bsize];
    int split_size = partition_subsize[PARTITION_SPLIT][bsize];

#define PARSE_BLOCK(db_r, db_c, db_subsize) \
    parse_block(dec_handle,                 \
                parse_ctx,                  \
                db_r,                       \
                db_c,                       \
                db_subsize,                 \
                &parse_ctx->cur_tile_info,  \
                sb_info,                    \
                partition);

#define PARSE_PARTITION(db_r, db_c, db_subsize) \
    parse_partition(dec_handle, parse_ctx, (db_r), (db_c), (db_subsize), sb_info)

    switch ((int)partition) {
    case PARTITION_NONE: PARSE_BLOCK(blk_row, blk_col, sub_size); break;
    case PARTITION_HORZ:
        PARSE_BLOCK(blk_row, blk_col, sub_size);
        if (has_rows) PARSE_BLOCK(blk_row + half_block_4x4, blk_col, sub_size);
        break;
    case PARTITION_VERT:
        PARSE_BLOCK(blk_row, blk_col, sub_size);
        if (has_cols) PARSE_BLOCK(blk_row, blk_col + half_block_4x4, sub_size);
        break;
    case PARTITION_SPLIT:
        PARSE_PARTITION(blk_row, blk_col, sub_size);
        PARSE_PARTITION(blk_row, blk_col + half_block_4x4, sub_size);
        PARSE_PARTITION(blk_row + half_block_4x4, blk_col, sub_size);
        PARSE_PARTITION(blk_row + half_block_4x4, blk_col + half_block_4x4, sub_size);
        break;
    case PARTITION_HORZ_A:
        PARSE_BLOCK(blk_row, blk_col, split_size);
        PARSE_BLOCK(blk_row, blk_col + half_block_4x4, split_size);
        PARSE_BLOCK(blk_row + half_block_4x4, blk_col, sub_size);
        break;
    case PARTITION_HORZ_B:
        PARSE_BLOCK(blk_row, blk_col, sub_size);
        PARSE_BLOCK(blk_row + half_block_4x4, blk_col, split_size);
        PARSE_BLOCK(blk_row + half_block_4x4, blk_col + half_block_4x4, split_size);
        break;
    case PARTITION_VERT_A:
        PARSE_BLOCK(blk_row, blk_col, split_size);
        PARSE_BLOCK(blk_row + half_block_4x4, blk_col, split_size);
        PARSE_BLOCK(blk_row, blk_col + half_block_4x4, sub_size);
        break;
    case PARTITION_VERT_B:
        PARSE_BLOCK(blk_row, blk_col, sub_size);
        PARSE_BLOCK(blk_row, blk_col + half_block_4x4, split_size);
        PARSE_BLOCK(blk_row + half_block_4x4, blk_col + half_block_4x4, split_size);
        break;
    case PARTITION_HORZ_4:
        for (int i = 0; i < 4; ++i) {
            uint32_t this_blk_row = blk_row + (uint32_t)(i * quarter_block_4x4);
            if (i > 0 && this_blk_row >= parse_ctx->frame_header->mi_rows) break;
            PARSE_BLOCK(this_blk_row, blk_col, sub_size);
        }
        break;
    case PARTITION_VERT_4:
        for (int i = 0; i < 4; ++i) {
            uint32_t this_blk_col = blk_col + (uint32_t)(i * quarter_block_4x4);
            if (i > 0 && this_blk_col >= parse_ctx->frame_header->mi_cols) break;
            PARSE_BLOCK(blk_row, this_blk_col, sub_size);
        }
        break;
    default: assert(0 && "Invalid partition type");
    }
    update_ext_partition_context(parse_ctx, blk_row, blk_col, sub_size, bsize, partition);
}

static INLINE int count_units_in_frame(int unitSize, int frame_size) {
    return MAX((frame_size + (unitSize >> 1)) / unitSize, 1);
}

static INLINE int decode_subexp_bool(int num_syms, int k, SvtReader *reader) {
    int i = 0, mk = 0;
    while (1) {
        int b2 = i ? k + i - 1 : k;
        int a  = 1 << b2;
        if (num_syms <= mk + 3 * a) {
            return svt_read_ns_ae(reader, num_syms - mk, ACCT_STR) + mk;
        } else {
            if (svt_read_literal(reader, 1, ACCT_STR)) {
                i++;
                mk += a;
            } else {
                return svt_read_literal(reader, b2, ACCT_STR) + mk;
            }
        }
    }
}

static INLINE int decode_unsigned_subexp_with_ref_bool(int mx, int k, int r, SvtReader *reader) {
    int v = decode_subexp_bool(mx, k, reader);
    if ((r << 1) <= mx) return inverse_recenter(r, v);
    return mx - 1 - inverse_recenter(mx - 1 - r, v);
}

static INLINE int decode_signed_subexp_with_ref_bool(int low, int high, int k, int r,
                                                     SvtReader *reader) {
    int x = decode_unsigned_subexp_with_ref_bool(high - low, k, r - low, reader);
    return x + low;
}
void read_wiener_filter(int wiener_win, WienerInfo *wiener_info, WienerInfo *ref_wiener_info,
                        SvtReader *reader) {
    memset(wiener_info->vfilter, 0, sizeof(wiener_info->vfilter));
    memset(wiener_info->hfilter, 0, sizeof(wiener_info->hfilter));

    // vfilter[0] and vfilter[6]
    if (wiener_win == WIENER_WIN) {
        wiener_info->vfilter[0] = wiener_info->vfilter[WIENER_WIN - 1] =
            decode_signed_subexp_with_ref_bool(WIENER_FILT_TAP0_MINV,
                                               WIENER_FILT_TAP0_MAXV + 1,
                                               WIENER_FILT_TAP0_SUBEXP_K,
                                               ref_wiener_info->vfilter[0],
                                               reader);
    } else
        wiener_info->vfilter[0] = wiener_info->vfilter[WIENER_WIN - 1] = 0;

    // vfilter[1] and vfilter[5]
    wiener_info->vfilter[1] = wiener_info->vfilter[WIENER_WIN - 2] =
        decode_signed_subexp_with_ref_bool(WIENER_FILT_TAP1_MINV,
                                           WIENER_FILT_TAP1_MAXV + 1,
                                           WIENER_FILT_TAP1_SUBEXP_K,
                                           ref_wiener_info->vfilter[1],
                                           reader);

    // vfilter[2] and vfilter[4]
    wiener_info->vfilter[2] = wiener_info->vfilter[WIENER_WIN - 3] =
        decode_signed_subexp_with_ref_bool(WIENER_FILT_TAP2_MINV,
                                           WIENER_FILT_TAP2_MAXV + 1,
                                           WIENER_FILT_TAP2_SUBEXP_K,
                                           ref_wiener_info->vfilter[2],
                                           reader);

    // vfilter[3] - The central element has an implicit +WIENER_FILT_STEP
    wiener_info->vfilter[WIENER_HALFWIN] =
        -2 * (wiener_info->vfilter[0] + wiener_info->vfilter[1] + wiener_info->vfilter[2]);

    // hfilter[0] and hfilter[6]
    if (wiener_win == WIENER_WIN) {
        wiener_info->hfilter[0] = wiener_info->hfilter[WIENER_WIN - 1] =
            decode_signed_subexp_with_ref_bool(WIENER_FILT_TAP0_MINV,
                                               WIENER_FILT_TAP0_MAXV + 1,
                                               WIENER_FILT_TAP0_SUBEXP_K,
                                               ref_wiener_info->hfilter[0],
                                               reader);
    } else
        wiener_info->hfilter[0] = wiener_info->hfilter[WIENER_WIN - 1] = 0;

    // hfilter[1] and hfilter[5]
    wiener_info->hfilter[1] = wiener_info->hfilter[WIENER_WIN - 2] =
        decode_signed_subexp_with_ref_bool(WIENER_FILT_TAP1_MINV,
                                           WIENER_FILT_TAP1_MAXV + 1,
                                           WIENER_FILT_TAP1_SUBEXP_K,
                                           ref_wiener_info->hfilter[1],
                                           reader);

    // hfilter[2] and hfilter[4]
    wiener_info->hfilter[2] = wiener_info->hfilter[WIENER_WIN - 3] =
        decode_signed_subexp_with_ref_bool(WIENER_FILT_TAP2_MINV,
                                           WIENER_FILT_TAP2_MAXV + 1,
                                           WIENER_FILT_TAP2_SUBEXP_K,
                                           ref_wiener_info->hfilter[2],
                                           reader);

    // hfilter[3] - The central element has an implicit +WIENER_FILT_STEP
    wiener_info->hfilter[WIENER_HALFWIN] =
        -2 * (wiener_info->hfilter[0] + wiener_info->hfilter[1] + wiener_info->hfilter[2]);

    eb_memcpy(ref_wiener_info, wiener_info, sizeof(*wiener_info));
}

void read_sgrproj_filter(SgrprojInfo *sgrproj_info, SgrprojInfo *ref_sgrproj_info,
                         SvtReader *reader) {
    sgrproj_info->ep = svt_read_literal(reader, SGRPROJ_PARAMS_BITS, ACCT_STR);
    int *r           = (int *)&eb_sgr_params[sgrproj_info->ep];

    if (r[0] == 0) {
        sgrproj_info->xqd[0] = 0;
        sgrproj_info->xqd[1] = decode_signed_subexp_with_ref_bool(SGRPROJ_PRJ_MIN1,
                                                                  SGRPROJ_PRJ_MAX1 + 1,
                                                                  SGRPROJ_PRJ_SUBEXP_K,
                                                                  ref_sgrproj_info->xqd[1],
                                                                  reader);
    } else if (r[1] == 0) {
        sgrproj_info->xqd[0] = decode_signed_subexp_with_ref_bool(SGRPROJ_PRJ_MIN0,
                                                                  SGRPROJ_PRJ_MAX0 + 1,
                                                                  SGRPROJ_PRJ_SUBEXP_K,
                                                                  ref_sgrproj_info->xqd[0],
                                                                  reader);
        sgrproj_info->xqd[1] = clamp(
            (1 << SGRPROJ_PRJ_BITS) - sgrproj_info->xqd[0], SGRPROJ_PRJ_MIN1, SGRPROJ_PRJ_MAX1);
    } else {
        sgrproj_info->xqd[0] = decode_signed_subexp_with_ref_bool(SGRPROJ_PRJ_MIN0,
                                                                  SGRPROJ_PRJ_MAX0 + 1,
                                                                  SGRPROJ_PRJ_SUBEXP_K,
                                                                  ref_sgrproj_info->xqd[0],
                                                                  reader);
        sgrproj_info->xqd[1] = decode_signed_subexp_with_ref_bool(SGRPROJ_PRJ_MIN1,
                                                                  SGRPROJ_PRJ_MAX1 + 1,
                                                                  SGRPROJ_PRJ_SUBEXP_K,
                                                                  ref_sgrproj_info->xqd[1],
                                                                  reader);
    }

    eb_memcpy(ref_sgrproj_info, sgrproj_info, sizeof(*sgrproj_info));
}

void read_lr_unit(ParseCtxt *parse_ctxt, int32_t plane, RestorationUnitInfo *lr_unit) {
    SvtReader *     reader     = &parse_ctxt->r;
    FrameHeader *   frame_info = parse_ctxt->frame_header;
    const LrParams *lrp        = &frame_info->lr_params[plane];
    if (lrp->frame_restoration_type == RESTORE_NONE) return;

    lr_unit->restoration_type = RESTORE_NONE;
    if (lrp->frame_restoration_type == RESTORE_SWITCHABLE) {
        lr_unit->restoration_type = svt_read_symbol(reader,
                                                    parse_ctxt->cur_tile_ctx.switchable_restore_cdf,
                                                    RESTORE_SWITCHABLE_TYPES,
                                                    ACCT_STR);
    } else if (lrp->frame_restoration_type == RESTORE_WIENER) {
        if (svt_read_symbol(reader, parse_ctxt->cur_tile_ctx.wiener_restore_cdf, 2, ACCT_STR)) {
            lr_unit->restoration_type = RESTORE_WIENER;
        }
    } else if (lrp->frame_restoration_type == RESTORE_SGRPROJ) {
        if (svt_read_symbol(reader, parse_ctxt->cur_tile_ctx.sgrproj_restore_cdf, 2, ACCT_STR)) {
            lr_unit->restoration_type = RESTORE_SGRPROJ;
        }
    }

    RestorationUnitInfo *ref_lr_plane = &parse_ctxt->ref_lr_unit[plane];
    WienerInfo *         wiener_info  = &lr_unit->wiener_info;
    SgrprojInfo *        sgrproj_info = &lr_unit->sgrproj_info;
    const int            wiener_win   = (plane > 0) ? WIENER_WIN_CHROMA : WIENER_WIN;

    switch (lr_unit->restoration_type) {
    case RESTORE_WIENER:
        read_wiener_filter(wiener_win, wiener_info, &ref_lr_plane->wiener_info, reader);
        break;
    case RESTORE_SGRPROJ:
        read_sgrproj_filter(sgrproj_info, &ref_lr_plane->sgrproj_info, reader);
        break;
    default: assert(lr_unit->restoration_type == RESTORE_NONE); break;
    }
}

void read_lr(EbDecHandle *dec_handle, ParseCtxt *parse_ctx, int32_t row, int32_t col) {
    FrameHeader *  frame_info   = parse_ctx->frame_header;
    SeqHeader *    seq_header   = parse_ctx->seq_header;
    EbColorConfig *color_config = &seq_header->color_config;
    FrameSize *    frame_size   = &frame_info->frame_size;
    if (frame_info->allow_intrabc) return;

    int     width      = mi_size_wide[seq_header->sb_size];
    int     height     = mi_size_high[seq_header->sb_size];
    int     num_planes = color_config->mono_chrome ? 1 : MAX_MB_PLANE;
    LrCtxt *lr_ctxt    = (LrCtxt *)dec_handle->pv_lr_ctxt;

    for (int plane = 0; plane < num_planes; plane++) {
        if (frame_info->lr_params[plane].frame_restoration_type != RESTORE_NONE) {
            int sub_x     = (plane == 0) ? 0 : color_config->subsampling_x;
            int sub_y     = (plane == 0) ? 0 : color_config->subsampling_y;
            int unit_size = frame_info->lr_params[plane].loop_restoration_size;
            int unit_rows = count_units_in_frame(
                unit_size, ROUND_POWER_OF_TWO(frame_size->frame_height, sub_y));
            int unit_cols = count_units_in_frame(
                unit_size, ROUND_POWER_OF_TWO(frame_size->superres_upscaled_width, sub_x));
            int unit_row_start = (row * (MI_SIZE >> sub_y) + unit_size - 1) / unit_size;
            int unit_row_end =
                MIN(unit_rows, ((row + height) * (MI_SIZE >> sub_y) + unit_size - 1) / unit_size);
            int numerator = 0, denominator = 0;
            if (!(frame_size->frame_width == frame_size->superres_upscaled_width)) {
                numerator   = (MI_SIZE >> sub_x) * frame_size->superres_denominator;
                denominator = unit_size * SCALE_NUMERATOR;
            } else {
                numerator   = MI_SIZE >> sub_x;
                denominator = unit_size;
            }
            int unit_col_start = (col * numerator + denominator - 1) / denominator;
            int unit_col_end =
                MIN(unit_cols, ((col + width) * numerator + denominator - 1) / denominator);
            for (int unit_row = unit_row_start; unit_row < unit_row_end; unit_row++) {
                for (int unit_col = unit_col_start; unit_col < unit_col_end; unit_col++) {
                    RestorationUnitInfo *cur_lr =
                        lr_ctxt->lr_unit[plane] + (unit_row * lr_ctxt->lr_stride[plane]) + unit_col;
                    read_lr_unit(parse_ctx, plane, cur_lr);
                }
            }
        }
    }
}

void parse_super_block(EbDecHandle *dec_handle, ParseCtxt *parse_ctx, uint32_t blk_row,
                       uint32_t blk_col, SBInfo *sb_info) {
    parse_ctx->read_deltas = parse_ctx->frame_header->delta_q_params.delta_q_present;

    read_lr(dec_handle, parse_ctx, blk_row, blk_col);

    parse_partition(
        dec_handle, parse_ctx, blk_row, blk_col, parse_ctx->seq_header->sb_size, sb_info);
}
