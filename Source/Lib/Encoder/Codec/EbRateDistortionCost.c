/*
* Copyright(c) 2019 Intel Corporation
* Copyright (c) 2016, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

/***************************************
* Includes
***************************************/
#include "EbRateDistortionCost.h"
#include "EbCommonUtils.h"
#include "aom_dsp_rtcd.h"
#include "EbLog.h"
#include "EbEncInterPrediction.h"
#include "EbFullLoop.h"
#include "EbEntropyCoding.h"

#include <assert.h>

#define AV1_COST_PRECISION 0
#define MV_COST_WEIGHT 108
int svt_aom_get_reference_mode_context_new(const MacroBlockD *xd);
int svt_av1_get_pred_context_uni_comp_ref_p(const MacroBlockD *xd);
int svt_av1_get_pred_context_uni_comp_ref_p1(const MacroBlockD *xd);
int svt_av1_get_pred_context_uni_comp_ref_p2(const MacroBlockD *xd);
int svt_aom_get_comp_reference_type_context_new(const MacroBlockD *xd);

int  svt_aom_get_palette_bsize_ctx(BlockSize bsize);
int  svt_aom_get_palette_mode_ctx(const MacroBlockD *xd);
int  svt_aom_write_uniform_cost(int n, int v);
int  svt_get_palette_cache_y(const MacroBlockD *const xd, uint16_t *cache);
int  svt_av1_palette_color_cost_y(const PaletteModeInfo *const pmi, uint16_t *color_cache, const int palette_size,
                                  int n_cache, int bit_depth);
int  svt_av1_cost_color_map(ModeDecisionCandidate *cand, MdRateEstimationContext *rate_table,

                            BlkStruct *blk_ptr, int plane, BlockSize bsize, COLOR_MAP_TYPE type);
void svt_aom_get_block_dimensions(BlockSize bsize, int plane, const MacroBlockD *xd, int *width, int *height,
                                  int *rows_within_bounds, int *cols_within_bounds);
int  svt_aom_allow_palette(int allow_screen_content_tools, BlockSize bsize);
int  svt_aom_allow_intrabc(const FrameHeader *frm_hdr, SliceType slice_type);
/* Symbols for coding which components are zero jointly */
//#define MV_JOINTS 4
//typedef enum {
//    MV_JOINT_ZERO = 0,   /* zero vector */
//    MV_JOINT_HNZVZ = 1,  /* Vert zero, hor nonzero */
//    MV_JOINT_HZVNZ = 2,  /* Hor zero, vert nonzero */
//    MV_JOINT_HNZVNZ = 3, /* Both components nonzero */
//} MvJointType;

MvJointType svt_av1_get_mv_joint(const MV *mv) {
    if (mv->row == 0)
        return mv->col == 0 ? MV_JOINT_ZERO : MV_JOINT_HNZVZ;
    else
        return mv->col == 0 ? MV_JOINT_HZVNZ : MV_JOINT_HNZVNZ;
}
static int32_t mv_cost(const MV *mv, const int32_t *joint_cost, int32_t *const comp_cost[2]) {
    int32_t jn_c = svt_av1_get_mv_joint(mv);
    int32_t res  = joint_cost[jn_c] + comp_cost[0][CLIP3(MV_LOW, MV_UPP, mv->row)] +
        comp_cost[1][CLIP3(MV_LOW, MV_UPP, mv->col)];
    return res;
}
int32_t svt_av1_mv_bit_cost_light(const MV *mv, const MV *ref) {
    const uint32_t factor     = 50;
    const uint32_t absmvdiffx = ABS(mv->col - ref->col);
    const uint32_t absmvdiffy = ABS(mv->row - ref->row);
    const uint32_t mv_rate    = 1296 + (factor * (absmvdiffx + absmvdiffy));
    return mv_rate;
}
int32_t svt_av1_mv_bit_cost(const MV *mv, const MV *ref, const int32_t *mvjcost, int32_t *mvcost[2], int32_t weight) {
    // Restrict the size of the MV diff to be within the max AV1 range.  If the MV diff
    // is outside this range, the diff will index beyond the cost array, causing a seg fault.
    // Both the MVs and the MV diffs should be within the allowable range for accessing the MV cost
    // infrastructure.
    MV temp_diff  = {mv->row - ref->row, mv->col - ref->col};
    temp_diff.row = MAX(temp_diff.row, MV_LOW);
    temp_diff.row = MIN(temp_diff.row, MV_UPP);
    temp_diff.col = MAX(temp_diff.col, MV_LOW);
    temp_diff.col = MIN(temp_diff.col, MV_UPP);

    const MV diff = temp_diff;
    return ROUND_POWER_OF_TWO(mv_cost(&diff, mvjcost, mvcost) * weight, 7);
}

/////////////////////////////COEFFICIENT CALCULATION //////////////////////////////////////////////
static INLINE int32_t get_golomb_cost(int32_t abs_qc) {
    if (abs_qc >= 1 + NUM_BASE_LEVELS + COEFF_BASE_RANGE) {
        const int32_t r      = abs_qc - COEFF_BASE_RANGE - NUM_BASE_LEVELS;
        const int32_t length = get_msb(r) + 1;
        return av1_cost_literal(2 * length - 1);
    }
    return 0;
}

void svt_av1_txb_init_levels_c(const TranLow *const coeff, const int32_t width, const int32_t height,
                               uint8_t *const levels) {
    const int32_t stride = width + TX_PAD_HOR;
    uint8_t      *ls     = levels;

    memset(levels - TX_PAD_TOP * stride, 0, sizeof(*levels) * TX_PAD_TOP * stride);
    memset(levels + stride * height, 0, sizeof(*levels) * (TX_PAD_BOTTOM * stride + TX_PAD_END));

    for (int32_t i = 0; i < height; i++) {
        for (int32_t j = 0; j < width; j++) *ls++ = (uint8_t)clamp(abs(coeff[i * width + j]), 0, INT8_MAX);
        for (int32_t j = 0; j < TX_PAD_HOR; j++) *ls++ = 0;
    }
}

static int32_t av1_transform_type_rate_estimation(struct ModeDecisionContext *ctx, uint8_t allow_update_cdf,
                                                  FRAME_CONTEXT *fc, struct ModeDecisionCandidateBuffer *cand_bf,
                                                  Bool is_inter, TxSize transform_size, TxType transform_type,
                                                  Bool reduced_tx_set_used) {
    // const MbModeInfo *mbmi = &xd->mi[0]->mbmi;
    // const int32_t is_inter = is_inter_block(mbmi);

    if (get_ext_tx_types(transform_size, is_inter, reduced_tx_set_used) >
        1 /*&&    !xd->lossless[xd->mi[0]->mbmi.segment_id]  WE ARE NOT LOSSLESS*/) {
        const TxSize square_tx_size = txsize_sqr_map[transform_size];
        assert(square_tx_size < EXT_TX_SIZES);

        const int32_t ext_tx_set = get_ext_tx_set(transform_size, is_inter, reduced_tx_set_used);
        if (is_inter) {
            if (ext_tx_set > 0) {
                if (allow_update_cdf) {
                    const TxSetType tx_set_type = get_ext_tx_set_type(transform_size, is_inter, reduced_tx_set_used);

                    update_cdf(fc->inter_ext_tx_cdf[ext_tx_set][square_tx_size],
                               av1_ext_tx_ind[tx_set_type][transform_type],
                               av1_num_ext_tx_set[tx_set_type]);
                }
                return ctx->md_rate_est_ctx->inter_tx_type_fac_bits[ext_tx_set][square_tx_size][transform_type];
            }
        } else {
            if (ext_tx_set > 0) {
                PredictionMode intra_dir;
                if (cand_bf->cand->filter_intra_mode != FILTER_INTRA_MODES)
                    intra_dir = fimode_to_intradir[cand_bf->cand->filter_intra_mode];
                else
                    intra_dir = cand_bf->cand->pred_mode;
                assert(intra_dir < INTRA_MODES);
                const TxSetType tx_set_type = get_ext_tx_set_type(transform_size, is_inter, reduced_tx_set_used);

                if (allow_update_cdf) {
                    update_cdf(fc->intra_ext_tx_cdf[ext_tx_set][square_tx_size][intra_dir],
                               av1_ext_tx_ind[tx_set_type][transform_type],
                               av1_num_ext_tx_set[tx_set_type]);
                }
                return ctx->md_rate_est_ctx
                    ->intra_tx_type_fac_bits[ext_tx_set][square_tx_size][intra_dir][transform_type];
            }
        }
    }
    return 0;
}

static const int8_t eob_to_pos_small[33] = {
    0, 1, 2, // 0-2
    3, 3, // 3-4
    4, 4, 4, 4, // 5-8
    5, 5, 5, 5, 5, 5, 5, 5, // 9-16
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 // 17-32
};

static const int8_t eob_to_pos_large[17] = {
    6, // place holder
    7, // 33-64
    8,
    8, // 65-128
    9,
    9,
    9,
    9, // 129-256
    10,
    10,
    10,
    10,
    10,
    10,
    10,
    10, // 257-512
    11 // 513-
};

static INLINE int32_t get_eob_pos_token(const int32_t eob, int32_t *const extra) {
    int32_t t;

    if (eob < 33)
        t = eob_to_pos_small[eob];
    else {
        const int32_t e = AOMMIN((eob - 1) >> 5, 16);
        t               = eob_to_pos_large[e];
    }

    *extra = eob - eb_k_eob_group_start[t];

    return t;
}
#define TX_SIZE TxSize
void svt_av1_update_eob_context(int eob, TxSize tx_size, TxClass tx_class, PlaneType plane, FRAME_CONTEXT *ec_ctx,
                                uint8_t allow_update_cdf) {
    int          eob_extra;
    const int    eob_pt  = get_eob_pos_token(eob, &eob_extra);
    const TxSize txs_ctx = (TxSize)((txsize_sqr_map[tx_size] + txsize_sqr_up_map[tx_size] + 1) >> 1);
    assert(txs_ctx < TX_SIZES);
    const int eob_multi_size = txsize_log2_minus4[tx_size];
    const int eob_multi_ctx  = (tx_class == TX_CLASS_2D) ? 0 : 1;

    switch (eob_multi_size) {
    case 0:
#if CONFIG_ENTROPY_STATS
        ++counts->eob_multi16[cdf_idx][plane][eob_multi_ctx][eob_pt - 1];
#endif
        if (allow_update_cdf)
            update_cdf(ec_ctx->eob_flag_cdf16[plane][eob_multi_ctx], eob_pt - 1, 5);
        break;
    case 1:
#if CONFIG_ENTROPY_STATS
        ++counts->eob_multi32[cdf_idx][plane][eob_multi_ctx][eob_pt - 1];
#endif
        if (allow_update_cdf)
            update_cdf(ec_ctx->eob_flag_cdf32[plane][eob_multi_ctx], eob_pt - 1, 6);
        break;
    case 2:
#if CONFIG_ENTROPY_STATS
        ++counts->eob_multi64[cdf_idx][plane][eob_multi_ctx][eob_pt - 1];
#endif
        if (allow_update_cdf)
            update_cdf(ec_ctx->eob_flag_cdf64[plane][eob_multi_ctx], eob_pt - 1, 7);
        break;
    case 3:
#if CONFIG_ENTROPY_STATS
        ++counts->eob_multi128[cdf_idx][plane][eob_multi_ctx][eob_pt - 1];
#endif
        if (allow_update_cdf) {
            update_cdf(ec_ctx->eob_flag_cdf128[plane][eob_multi_ctx], eob_pt - 1, 8);
        }
        break;
    case 4:
#if CONFIG_ENTROPY_STATS
        ++counts->eob_multi256[cdf_idx][plane][eob_multi_ctx][eob_pt - 1];
#endif
        if (allow_update_cdf) {
            update_cdf(ec_ctx->eob_flag_cdf256[plane][eob_multi_ctx], eob_pt - 1, 9);
        }
        break;
    case 5:
#if CONFIG_ENTROPY_STATS
        ++counts->eob_multi512[cdf_idx][plane][eob_multi_ctx][eob_pt - 1];
#endif
        if (allow_update_cdf) {
            update_cdf(ec_ctx->eob_flag_cdf512[plane][eob_multi_ctx], eob_pt - 1, 10);
        }
        break;
    case 6:
    default:
#if CONFIG_ENTROPY_STATS
        ++counts->eob_multi1024[cdf_idx][plane][eob_multi_ctx][eob_pt - 1];
#endif
        if (allow_update_cdf) {
            update_cdf(ec_ctx->eob_flag_cdf1024[plane][eob_multi_ctx], eob_pt - 1, 11);
        }
        break;
    }

    const int16_t eob_offset_bits = eb_k_eob_offset_bits[eob_pt];
    if (eob_offset_bits > 0) {
        int eob_shift = eob_offset_bits - 1;
        int bit       = (eob_extra & (1 << eob_shift)) ? 1 : 0;
#if CONFIG_ENTROPY_STATS
        counts->eob_extra[cdf_idx][txs_ctx][plane][eob_pt][bit]++;
#endif // CONFIG_ENTROPY_STATS
        if (allow_update_cdf)
            update_cdf(ec_ctx->eob_extra_cdf[txs_ctx][plane][eob_pt], bit, 2);
    }
}
// Transform end of block bit estimation
static int get_eob_cost(int eob, const LvMapEobCost *txb_eob_costs, const LvMapCoeffCost *txb_costs, TxClass tx_class) {
    int       eob_extra;
    const int eob_pt        = get_eob_pos_token(eob, &eob_extra);
    int       eob_cost      = 0;
    const int eob_multi_ctx = (tx_class == TX_CLASS_2D) ? 0 : 1;
    eob_cost                = txb_eob_costs->eob_cost[eob_multi_ctx][eob_pt - 1];

    if (eb_k_eob_offset_bits[eob_pt] > 0) {
        const int eob_ctx   = eob_pt - 3;
        const int eob_shift = eb_k_eob_offset_bits[eob_pt] - 1;
        const int bit       = (eob_extra & (1 << eob_shift)) ? 1 : 0;
        eob_cost += txb_costs->eob_extra_cost[eob_ctx][bit];
        const int offset_bits = eb_k_eob_offset_bits[eob_pt];
        if (offset_bits > 1)
            eob_cost += av1_cost_literal(offset_bits - 1);
    }
    return eob_cost;
}
static INLINE int32_t av1_cost_skip_txb(struct ModeDecisionContext *ctx, uint8_t allow_update_cdf,
                                        FRAME_CONTEXT *ec_ctx, TxSize transform_size, PlaneType plane_type,
                                        int16_t txb_skip_ctx) {
    const TxSize txs_ctx = (TxSize)((txsize_sqr_map[transform_size] + txsize_sqr_up_map[transform_size] + 1) >> 1);
    assert(txs_ctx < TX_SIZES);
    const LvMapCoeffCost *const coeff_costs = &ctx->md_rate_est_ctx->coeff_fac_bits[txs_ctx][plane_type];
    if (allow_update_cdf)
        update_cdf(ec_ctx->txb_skip_cdf[txs_ctx][txb_skip_ctx], 1, 2);
    return coeff_costs->txb_skip_cost[txb_skip_ctx][1];
}

static INLINE int32_t av1_cost_coeffs_txb_loop_cost_one_eob(const TranLow *const qcoeff, int8_t *const coeff_contexts,
                                                            const LvMapCoeffCost *coeff_costs, int16_t dc_sign_ctx) {
    const TranLow v         = qcoeff[0];
    const int32_t level     = abs(v);
    const int32_t coeff_ctx = coeff_contexts[0];

    assert((AOMMIN(level, 3) - 1) >= 0);
    int32_t cost = coeff_costs->base_eob_cost[coeff_ctx][AOMMIN(level, 3) - 1];

    if (v != 0) {
        const int32_t sign = (v < 0) ? 1 : 0;
        // sign bit cost
        cost += coeff_costs->dc_sign_cost[dc_sign_ctx][sign];

        if (level > NUM_BASE_LEVELS) {
            const int32_t base_range = level - 1 - NUM_BASE_LEVELS;

            if (base_range < COEFF_BASE_RANGE)
                cost += coeff_costs->lps_cost[0][base_range];
            else
                cost += coeff_costs->lps_cost[0][COEFF_BASE_RANGE];

            if (level >= 1 + NUM_BASE_LEVELS + COEFF_BASE_RANGE)
                cost += get_golomb_cost(level);
        }
    }
    return cost;
}

static INLINE int32_t av1_cost_coeffs_txb_loop_cost_eob(struct ModeDecisionContext *md_ctx, uint16_t eob,
                                                        const int16_t *const scan, const TranLow *const qcoeff,
                                                        int8_t *const coeff_contexts, const LvMapCoeffCost *coeff_costs,
                                                        int16_t dc_sign_ctx, uint8_t *const levels, const int32_t bwl,
                                                        TxType transform_type) {
    const uint32_t cost_literal = av1_cost_literal(1);
    int32_t        cost         = 0;

    //Optimized/simplified function when eob is 1
    if (eob == 1)
        return av1_cost_coeffs_txb_loop_cost_one_eob(qcoeff, coeff_contexts, coeff_costs, dc_sign_ctx);

    //  first (eob - 1) index
    {
        const int32_t pos       = scan[eob - 1];
        const TranLow v         = qcoeff[pos];
        const int32_t level     = abs(v);
        const int32_t coeff_ctx = coeff_contexts[pos];

        assert((AOMMIN(level, 3) - 1) >= 0);
        cost += coeff_costs->base_eob_cost[coeff_ctx][AOMMIN(level, 3) - 1];

        if (v != 0) {
            cost += cost_literal;
            if (level > NUM_BASE_LEVELS) {
                int32_t       ctx        = get_br_ctx(levels, pos, bwl, tx_type_to_class[transform_type]);
                const int32_t base_range = level - 1 - NUM_BASE_LEVELS;

                if (base_range < COEFF_BASE_RANGE)
                    cost += coeff_costs->lps_cost[ctx][base_range];
                else
                    cost += coeff_costs->lps_cost[ctx][COEFF_BASE_RANGE];

                if (level >= 1 + NUM_BASE_LEVELS + COEFF_BASE_RANGE)
                    cost += get_golomb_cost(level);
            }
        }
    }
    // last (0) index
    {
        const TranLow v         = qcoeff[0];
        const int32_t level     = abs(v);
        const int32_t coeff_ctx = coeff_contexts[0];

        cost += coeff_costs->base_cost[coeff_ctx][AOMMIN(level, 3)];

        if (v != 0) {
            const int32_t sign = (v < 0) ? 1 : 0;
            // sign bit cost

            cost += coeff_costs->dc_sign_cost[dc_sign_ctx][sign];

            if (level > NUM_BASE_LEVELS) {
                int32_t       ctx        = get_br_ctx(levels, 0, bwl, tx_type_to_class[transform_type]);
                const int32_t base_range = level - 1 - NUM_BASE_LEVELS;

                if (base_range < COEFF_BASE_RANGE)
                    cost += coeff_costs->lps_cost[ctx][base_range];
                else
                    cost += coeff_costs->lps_cost[ctx][COEFF_BASE_RANGE];

                if (level >= 1 + NUM_BASE_LEVELS + COEFF_BASE_RANGE)
                    cost += get_golomb_cost(level);
            }
        }
    }
    int32_t c;
    /* Optimized Loop, omitted first (eob - 1) and last (0) index */
    // Estimate the rate of the first(eob / fast_coeff_est_level) coeff(s), DC and last coeff only
    int32_t  c_start = MIN(eob - 2, eob / MAX(1, (int)(md_ctx->mds_fast_coeff_est_level - md_ctx->mds_subres_step)));
    uint32_t cost_literal_cnt = 0;
    for (c = c_start; c >= 1; --c) {
        const int32_t pos = scan[c];
        cost_literal_cnt += !!(qcoeff[pos]);
        const int32_t level = abs(qcoeff[pos]);
        if (level > NUM_BASE_LEVELS) {
            int32_t       ctx        = get_br_ctx(levels, pos, bwl, tx_type_to_class[transform_type]);
            const int32_t base_range = level - 1 - NUM_BASE_LEVELS;

            cost += coeff_costs->base_cost[coeff_contexts[pos]][3];
            if (base_range < COEFF_BASE_RANGE) {
                cost += coeff_costs->lps_cost[ctx][base_range];
            } else {
                cost += get_golomb_cost(level) + coeff_costs->lps_cost[ctx][COEFF_BASE_RANGE];
            }
        } else {
            cost += coeff_costs->base_cost[coeff_contexts[pos]][level];
        }
    }
    cost += cost_literal_cnt * cost_literal;

    return cost;
}

// Note: don't call this function when eob is 0.
uint64_t svt_av1_cost_coeffs_txb(struct ModeDecisionContext *ctx, uint8_t allow_update_cdf, FRAME_CONTEXT *ec_ctx,
                                 struct ModeDecisionCandidateBuffer *cand_bf, const TranLow *const qcoeff, uint16_t eob,
                                 PlaneType plane_type, TxSize transform_size, TxType transform_type,
                                 int16_t txb_skip_ctx, int16_t dc_sign_ctx, Bool reduced_transform_set_flag)

{
    //Note: there is a different version of this function in AOM that seems to be efficient as its name is:
    //warehouse_efficients_txb

    const TxSize  txs_ctx  = (TxSize)((txsize_sqr_map[transform_size] + txsize_sqr_up_map[transform_size] + 1) >> 1);
    const TxClass tx_class = tx_type_to_class[transform_type];
    int32_t       cost;
    const int32_t bwl                 = get_txb_bwl_tab[transform_size];
    const int32_t width               = get_txb_wide_tab[transform_size];
    const int32_t height              = get_txb_high_tab[transform_size];
    const ScanOrder *const scan_order = &av1_scan_orders[transform_size][transform_type]; // get_scan(tx_size, tx_type);
    const int16_t *const   scan       = scan_order->scan;
    uint8_t                levels_buf[TX_PAD_2D];
    uint8_t *const         levels = set_levels(levels_buf, width);
    DECLARE_ALIGNED(16, int8_t, coeff_contexts[MAX_TX_SQUARE]);
    assert(txs_ctx < TX_SIZES);
    const LvMapCoeffCost *const coeff_costs = &ctx->md_rate_est_ctx->coeff_fac_bits[txs_ctx][plane_type];

    const int32_t             eob_multi_size = txsize_log2_minus4[transform_size];
    const LvMapEobCost *const eob_bits       = &ctx->md_rate_est_ctx->eob_frac_bits[eob_multi_size][plane_type];
    // eob must be greater than 0 here.
    assert(eob > 0);
    cost = coeff_costs->txb_skip_cost[txb_skip_ctx][0];

    if (allow_update_cdf)
        update_cdf(ec_ctx->txb_skip_cdf[txs_ctx][txb_skip_ctx], eob == 0, 2);

    if (eob > 1)
        svt_av1_txb_init_levels(qcoeff,
                                width,
                                height,
                                levels); // NM - Needs to be optimized - to be combined with the quantisation.
    const Bool is_inter = is_inter_mode(cand_bf->cand->pred_mode);
    // Transform type bit estimation
    cost += plane_type > PLANE_TYPE_Y ? 0
                                      : av1_transform_type_rate_estimation(ctx,
                                                                           allow_update_cdf,
                                                                           ec_ctx,
                                                                           cand_bf,
                                                                           is_inter,
                                                                           transform_size,
                                                                           transform_type,
                                                                           reduced_transform_set_flag);

    // Transform eob bit estimation
    int32_t eob_cost = get_eob_cost(eob, eob_bits, coeff_costs, tx_class);
    cost += eob_cost;
    if (allow_update_cdf)
        svt_av1_update_eob_context(eob, transform_size, tx_class, plane_type, ec_ctx, allow_update_cdf);
    // Transform non-zero coeff bit estimation
    svt_av1_get_nz_map_contexts(levels,
                                scan,
                                eob,
                                transform_size,
                                tx_class,
                                coeff_contexts); // NM - Assembly version is available in AOM
    assert(eob <= width * height);
    if (allow_update_cdf) {
        for (int c = eob - 1; c >= 0; --c) {
            const int     pos       = scan[c];
            const int     coeff_ctx = coeff_contexts[pos];
            const TranLow v         = qcoeff[pos];
            const TranLow level     = abs(v);
            if (c == eob - 1) {
                assert(coeff_ctx < 4);
                update_cdf(ec_ctx->coeff_base_eob_cdf[txs_ctx][plane_type][coeff_ctx], AOMMIN(level, 3) - 1, 3);
            } else {
                update_cdf(ec_ctx->coeff_base_cdf[txs_ctx][plane_type][coeff_ctx], AOMMIN(level, 3), 4);
            }

            {
                if (c == eob - 1) {
                    assert(coeff_ctx < 4);
#if CONFIG_ENTROPY_STATS
                    ++td->counts
                          ->coeff_base_eob_multi[cdf_idx][txsize_ctx][plane_type][coeff_ctx][AOMMIN(level, 3) - 1];
                } else {
                    ++td->counts->coeff_base_multi[cdf_idx][txsize_ctx][plane_type][coeff_ctx][AOMMIN(level, 3)];
#endif
                }
            }

            if (level > NUM_BASE_LEVELS) {
                const int base_range = level - 1 - NUM_BASE_LEVELS;
                int       br_ctx;
                if (eob == 1)
                    br_ctx = 0;
                else
                    br_ctx = get_br_ctx(levels, pos, bwl, tx_class);

                for (int idx = 0; idx < COEFF_BASE_RANGE; idx += BR_CDF_SIZE - 1) {
                    const int k = AOMMIN(base_range - idx, BR_CDF_SIZE - 1);
                    update_cdf(ec_ctx->coeff_br_cdf[AOMMIN(txs_ctx, TX_32X32)][plane_type][br_ctx], k, BR_CDF_SIZE);
                    for (int lps = 0; lps < BR_CDF_SIZE - 1; lps++) {
#if CONFIG_ENTROPY_STATS
                        ++td->counts->coeff_lps[AOMMIN(txsize_ctx, TX_32X32)][plane_type][lps][br_ctx][lps == k];
#endif // CONFIG_ENTROPY_STATS
                        if (lps == k)
                            break;
                    }
#if CONFIG_ENTROPY_STATS
                    ++td->counts->coeff_lps_multi[cdf_idx][AOMMIN(txsize_ctx, TX_32X32)][plane_type][br_ctx][k];
#endif
                    if (k < BR_CDF_SIZE - 1)
                        break;
                }
            }
        }

        if (qcoeff[0] != 0)
            update_cdf(ec_ctx->dc_sign_cdf[plane_type][dc_sign_ctx], qcoeff[0] < 0, 2);

        //TODO: CHKN  for 128x128 where we need more than one TXb, we need to update the txb_context(dc_sign+skip_ctx) in a Txb basis.

        return 0;
    }

    cost += av1_cost_coeffs_txb_loop_cost_eob(
        ctx, eob, scan, qcoeff, coeff_contexts, coeff_costs, dc_sign_ctx, levels, bwl, transform_type);
    return cost;
}
uint64_t svt_aom_get_intra_uv_fast_rate(PictureControlSet *pcs, struct ModeDecisionContext *ctx,
                                        ModeDecisionCandidateBuffer *cand_bf, bool use_accurate_cfl) {
    const BlockGeom *const blk_geom = ctx->blk_geom;
    ModeDecisionCandidate *cand     = cand_bf->cand;
    assert(blk_geom->has_uv);
    assert(!(svt_aom_allow_intrabc(&pcs->ppcs->frm_hdr, pcs->ppcs->slice_type) && cand->use_intrabc));
    MdRateEstimationContext *md_rate_est_ctx = ctx->md_rate_est_ctx;
    const uint8_t            is_cfl_allowed  = (blk_geom->bwidth <= 32 && blk_geom->bheight <= 32) ? 1 : 0;
    PredictionMode           intra_mode      = (PredictionMode)cand->pred_mode;
    // If CFL alphas are not known yet, calculate the chroma mode bits based on DC Mode. If CFL is selected the chroma mode bits must be updated later
    const UvPredictionMode chroma_mode = cand->intra_chroma_mode == UV_CFL_PRED && !use_accurate_cfl
        ? UV_DC_PRED
        : cand->intra_chroma_mode;
    const uint32_t         mi_row      = ctx->blk_org_y >> MI_SIZE_LOG2;
    const uint32_t         mi_col      = ctx->blk_org_x >> MI_SIZE_LOG2;
    // Subsampling assumes YUV 420 content
    const uint8_t ss_x = 1;
    const uint8_t ss_y = 1;

    uint64_t chroma_rate = 0;
    // Estimate chroma nominal intra mode bits
    chroma_rate += (uint64_t)md_rate_est_ctx->intra_uv_mode_fac_bits[is_cfl_allowed][intra_mode][chroma_mode];

    // Estimate chroma angular mode bits; angular offset only allow for bsize >= 8x8
    if (blk_geom->bsize >= BLOCK_8X8 && av1_is_directional_mode(get_uv_mode(chroma_mode))) {
        chroma_rate += md_rate_est_ctx->angle_delta_fac_bits[chroma_mode - V_PRED]
                                                            [MAX_ANGLE_DELTA + cand->angle_delta[PLANE_TYPE_UV]];
    }

    // Estimate CFL factor bits when CFL is used
    if (chroma_mode == UV_CFL_PRED) {
        chroma_rate += (uint64_t)md_rate_est_ctx
                           ->cfl_alpha_fac_bits[cand->cfl_alpha_signs][CFL_PRED_U][CFL_IDX_U(cand->cfl_alpha_idx)] +
            (uint64_t)
                md_rate_est_ctx->cfl_alpha_fac_bits[cand->cfl_alpha_signs][CFL_PRED_V][CFL_IDX_V(cand->cfl_alpha_idx)];
    }

    // Estimate chroma palette mode bits (currently not supported, so just cost of signalling off)
    if (chroma_mode == UV_DC_PRED &&
        svt_aom_allow_palette(pcs->ppcs->frm_hdr.allow_screen_content_tools, blk_geom->bsize) &&
        is_chroma_reference(mi_row, mi_col, blk_geom->bsize, ss_x, ss_y)) {
        const int use_palette_y  = cand->palette_info && (cand->palette_size[0] > 0);
        const int use_palette_uv = cand->palette_info && (cand->palette_size[1] > 0);
        chroma_rate += ctx->md_rate_est_ctx->palette_uv_mode_fac_bits[use_palette_y][use_palette_uv];
    }

    return chroma_rate;
}
uint64_t svt_aom_intra_fast_cost(PictureControlSet *pcs, struct ModeDecisionContext *ctx,
                                 ModeDecisionCandidateBuffer *cand_bf, uint64_t lambda, uint64_t luma_distortion,
                                 uint64_t chroma_distortion) {
    const BlockGeom       *blk_geom = ctx->blk_geom;
    BlkStruct             *blk_ptr  = ctx->blk_ptr;
    ModeDecisionCandidate *cand     = cand_bf->cand;
    if (svt_aom_allow_intrabc(&pcs->ppcs->frm_hdr, pcs->ppcs->slice_type) && cand->use_intrabc) {
        uint64_t rate = 0;

        RefList ref_list_idx = 0;
        int16_t pred_ref_x   = cand->pred_mv[ref_list_idx].x;
        int16_t pred_ref_y   = cand->pred_mv[ref_list_idx].y;
        int16_t mv_ref_x     = cand->mv[ref_list_idx].x;
        int16_t mv_ref_y     = cand->mv[ref_list_idx].y;
        MV      mv;
        mv.row = mv_ref_y;
        mv.col = mv_ref_x;
        MV ref_mv;
        ref_mv.row        = pred_ref_y;
        ref_mv.col        = pred_ref_x;
        int    *dvcost[2] = {(int *)&ctx->md_rate_est_ctx->dv_cost[0][MV_MAX],
                             (int *)&ctx->md_rate_est_ctx->dv_cost[1][MV_MAX]};
        int32_t mv_rate   = svt_av1_mv_bit_cost(
            &mv, &ref_mv, ctx->md_rate_est_ctx->dv_joint_cost, dvcost, MV_COST_WEIGHT_SUB);

        rate                      = mv_rate + ctx->md_rate_est_ctx->intrabc_fac_bits[cand->use_intrabc];
        cand_bf->fast_luma_rate   = RATE_WEIGHT * rate;
        cand_bf->fast_chroma_rate = 0;
        uint64_t luma_sad         = (LUMA_WEIGHT * luma_distortion) << AV1_COST_PRECISION;
        uint64_t chromasad_       = (CHROMA_WEIGHT * chroma_distortion) << AV1_COST_PRECISION;
        uint64_t total_distortion = luma_sad + chromasad_;

        return (RDCOST(lambda, rate, total_distortion));
    } else {
        // Number of bits for each synatax element
        uint64_t       intra_mode_bits_num          = 0;
        uint64_t       intra_luma_mode_bits_num     = 0;
        uint64_t       intra_luma_ang_mode_bits_num = 0;
        uint64_t       intra_filter_mode_bits_num   = 0;
        uint64_t       skip_mode_rate               = 0;
        const uint8_t  skip_mode_ctx                = ctx->skip_mode_ctx;
        PredictionMode intra_mode                   = (PredictionMode)cand->pred_mode;
        // Luma and chroma rate
        uint32_t rate;
        uint32_t luma_rate   = 0;
        uint32_t chroma_rate = 0;
        uint64_t luma_sad, chromasad_;
        assert(intra_mode < INTRA_MODES);
        // Luma and chroma distortion
        uint64_t total_distortion;
        intra_mode_bits_num = pcs->slice_type != I_SLICE
            ? (uint64_t)ctx->md_rate_est_ctx->mb_mode_fac_bits[size_group_lookup[blk_geom->bsize]][intra_mode]
            : ZERO_COST;

        skip_mode_rate = pcs->slice_type != I_SLICE && pcs->ppcs->frm_hdr.skip_mode_params.skip_mode_flag &&
                is_comp_ref_allowed(blk_geom->bsize)
            ? (uint64_t)ctx->md_rate_est_ctx->skip_mode_fac_bits[skip_mode_ctx][0]
            : ZERO_COST;
        // Estimate luma nominal intra mode bits for key frame
        intra_luma_mode_bits_num = pcs->slice_type == I_SLICE
            ? (uint64_t)
                  ctx->md_rate_est_ctx->y_mode_fac_bits[ctx->intra_luma_top_ctx][ctx->intra_luma_left_ctx][intra_mode]
            : ZERO_COST;
        // Estimate luma angular mode bits
        if (blk_geom->bsize >= BLOCK_8X8 && av1_is_directional_mode(cand->pred_mode)) {
            assert((intra_mode - V_PRED) < 8);
            assert((intra_mode - V_PRED) >= 0);
            intra_luma_ang_mode_bits_num =
                ctx->md_rate_est_ctx
                    ->angle_delta_fac_bits[intra_mode - V_PRED][MAX_ANGLE_DELTA + cand->angle_delta[PLANE_TYPE_Y]];
        }
        if (svt_aom_allow_palette(pcs->ppcs->frm_hdr.allow_screen_content_tools, blk_geom->bsize) &&
            intra_mode == DC_PRED) {
            const int use_palette = cand->palette_info ? (cand->palette_size[0] > 0) : 0;
            const int bsize_ctx   = svt_aom_get_palette_bsize_ctx(blk_geom->bsize);
            const int mode_ctx    = svt_aom_get_palette_mode_ctx(blk_ptr->av1xd);
            intra_luma_mode_bits_num += ctx->md_rate_est_ctx->palette_ymode_fac_bits[bsize_ctx][mode_ctx][use_palette];
            if (use_palette) {
                const uint8_t *const color_map = cand->palette_info->color_idx_map;
                int                  block_width, block_height, rows, cols;
                svt_aom_get_block_dimensions(
                    blk_geom->bsize, 0, blk_ptr->av1xd, &block_width, &block_height, &rows, &cols);
                const int plt_size = cand->palette_size[0];
                int       palette_mode_cost =
                    ctx->md_rate_est_ctx->palette_ysize_fac_bits[bsize_ctx][plt_size - PALETTE_MIN_SIZE] +
                    svt_aom_write_uniform_cost(plt_size, color_map[0]);
                uint16_t  color_cache[2 * PALETTE_MAX_SIZE];
                const int n_cache = svt_get_palette_cache_y(blk_ptr->av1xd, color_cache);
                palette_mode_cost += svt_av1_palette_color_cost_y(&cand->palette_info->pmi,
                                                                  color_cache,
                                                                  cand->palette_size[0],
                                                                  n_cache,
                                                                  pcs->ppcs->scs->encoder_bit_depth);
                if (!pcs->ppcs->palette_ctrls.reduce_palette_cost_precision)
                    palette_mode_cost += svt_av1_cost_color_map(
                        cand, ctx->md_rate_est_ctx, blk_ptr, 0, blk_geom->bsize, PALETTE_MAP);
                intra_luma_mode_bits_num += palette_mode_cost;
            }
        }

        if (svt_aom_filter_intra_allowed(pcs->ppcs->scs->seq_header.filter_intra_level,
                                         blk_geom->bsize,
                                         cand->palette_info ? cand->palette_size[0] : 0,
                                         intra_mode)) {
            intra_filter_mode_bits_num =
                ctx->md_rate_est_ctx
                    ->filter_intra_fac_bits[blk_geom->bsize][cand->filter_intra_mode != FILTER_INTRA_MODES];
            if (cand->filter_intra_mode != FILTER_INTRA_MODES) {
                intra_filter_mode_bits_num += ctx->md_rate_est_ctx->filter_intra_mode_fac_bits[cand->filter_intra_mode];
            }
        }
        if (blk_geom->has_uv) {
            // CFL info not known in fasta loop, so assume DC mode when CFL is allowed
            chroma_rate = (uint32_t)svt_aom_get_intra_uv_fast_rate(pcs, ctx, cand_bf, 0);
        }

        uint32_t is_inter_rate = pcs->slice_type != I_SLICE
            ? ctx->md_rate_est_ctx->intra_inter_fac_bits[ctx->is_inter_ctx][0]
            : 0;
        luma_rate              = (uint32_t)(intra_mode_bits_num + skip_mode_rate + intra_luma_mode_bits_num +
                               intra_luma_ang_mode_bits_num + is_inter_rate + intra_filter_mode_bits_num);
        if (svt_aom_allow_intrabc(&pcs->ppcs->frm_hdr, pcs->ppcs->slice_type)) {
            svt_aom_assert_err(cand->use_intrabc == 0, "this block ibc should be off\n");
            luma_rate += ctx->md_rate_est_ctx->intrabc_fac_bits[cand->use_intrabc];
        }
        // Keep the Fast Luma and Chroma rate for future use
        cand_bf->fast_luma_rate   = luma_rate * RATE_WEIGHT;
        cand_bf->fast_chroma_rate = chroma_rate * RATE_WEIGHT;
        luma_sad                  = (LUMA_WEIGHT * luma_distortion) << AV1_COST_PRECISION;
        chromasad_                = (CHROMA_WEIGHT * chroma_distortion) << AV1_COST_PRECISION;
        total_distortion          = luma_sad + chromasad_;

        rate = luma_rate + chroma_rate;

        // Assign fast cost
        return (RDCOST(lambda, rate, total_distortion));
    }
}
static INLINE int svt_aom_has_second_ref(const MbModeInfo *mbmi) { return mbmi->block_mi.ref_frame[1] > INTRA_FRAME; }

static INLINE int svt_aom_has_uni_comp_refs(const MbModeInfo *mbmi) {
    return svt_aom_has_second_ref(mbmi) &&
        (!((mbmi->block_mi.ref_frame[0] >= BWDREF_FRAME) ^ (mbmi->block_mi.ref_frame[1] >= BWDREF_FRAME)));
}

// This function encodes the reference frame
uint64_t estimate_ref_frame_type_bits(struct ModeDecisionContext *ctx, BlkStruct *blk_ptr, uint8_t ref_frame_type,
                                      Bool is_compound) {
    uint64_t ref_rate_bits = 0;

    // const MbModeInfo *const mbmi = &blk_ptr->av1xd->mi[0]->mbmi;
    MbModeInfo *const mbmi = &blk_ptr->av1xd->mi[0]->mbmi;
    MvReferenceFrame  ref_type[2];
    av1_set_ref_frame(ref_type, ref_frame_type);
    mbmi->block_mi.ref_frame[0] = ref_type[0];
    mbmi->block_mi.ref_frame[1] = ref_type[1];
    //const int is_compound = svt_aom_has_second_ref(mbmi);
    {
        if (is_compound) {
            const CompReferenceType comp_ref_type = svt_aom_has_uni_comp_refs(mbmi) ? UNIDIR_COMP_REFERENCE
                                                                                    : BIDIR_COMP_REFERENCE;

            ref_rate_bits += ctx->md_rate_est_ctx->comp_ref_type_fac_bits[svt_aom_get_comp_reference_type_context_new(
                blk_ptr->av1xd)][comp_ref_type];
            /*aom_write_symbol(w, comp_ref_type,
               svt_aom_get_comp_reference_type_cdf(blk_ptr->av1xd), 2);*/

            if (comp_ref_type == UNIDIR_COMP_REFERENCE) {
                // SVT_LOG("ERROR[AN]: UNIDIR_COMP_REFERENCE not supported\n");
                const int bit = mbmi->block_mi.ref_frame[0] == BWDREF_FRAME;

                ref_rate_bits += ctx->md_rate_est_ctx->uni_comp_ref_fac_bits[svt_av1_get_pred_context_uni_comp_ref_p(
                    blk_ptr->av1xd)][0][bit];
                // blk_ptr->av1xd->tile_ctx->uni_comp_ref_cdf[pred_context][0];
                // WRITE_REF_BIT(bit, uni_comp_ref_p);

                if (!bit) {
                    assert(mbmi->block_mi.ref_frame[0] == LAST_FRAME);
                    const int bit1 = mbmi->block_mi.ref_frame[1] == LAST3_FRAME ||
                        mbmi->block_mi.ref_frame[1] == GOLDEN_FRAME;
                    ref_rate_bits +=
                        ctx->md_rate_est_ctx
                            ->uni_comp_ref_fac_bits[svt_av1_get_pred_context_uni_comp_ref_p1(blk_ptr->av1xd)][1][bit1];
                    // ref_rate_d = blk_ptr->av1xd->tile_ctx->uni_comp_ref_cdf[pred_context][1];
                    // WRITE_REF_BIT(bit1, uni_comp_ref_p1);
                    if (bit1) {
                        const int bit2 = mbmi->block_mi.ref_frame[1] == GOLDEN_FRAME;
                        ref_rate_bits +=
                            ctx->md_rate_est_ctx->uni_comp_ref_fac_bits[svt_av1_get_pred_context_uni_comp_ref_p2(
                                blk_ptr->av1xd)][2][bit2];

                        // ref_rate_e = blk_ptr->av1xd->tile_ctx->uni_comp_ref_cdf[pred_context][2];
                        //WRITE_REF_BIT(bit2, uni_comp_ref_p2);
                    }
                }
                return ref_rate_bits;
                //return;
            }

            assert(comp_ref_type == BIDIR_COMP_REFERENCE);

            const int bit = (mbmi->block_mi.ref_frame[0] == GOLDEN_FRAME || mbmi->block_mi.ref_frame[0] == LAST3_FRAME);
            const int pred_ctx = svt_av1_get_pred_context_comp_ref_p(blk_ptr->av1xd);
            ref_rate_bits += ctx->md_rate_est_ctx->comp_ref_fac_bits[pred_ctx][0][bit];
            // ref_rate_f = blk_ptr->av1xd->tile_ctx->comp_ref_cdf[pred_ctx][0];
            // WRITE_REF_BIT(bit, comp_ref_p);

            if (!bit) {
                const int bit1 = mbmi->block_mi.ref_frame[0] == LAST2_FRAME;
                ref_rate_bits += ctx->md_rate_est_ctx
                                     ->comp_ref_fac_bits[svt_av1_get_pred_context_comp_ref_p1(blk_ptr->av1xd)][1][bit1];
                // ref_rate_g = blk_ptr->av1xd->tile_ctx->comp_ref_cdf[pred_context][1];
                // WRITE_REF_BIT(bit1, comp_ref_p1);
            } else {
                const int bit2 = mbmi->block_mi.ref_frame[0] == GOLDEN_FRAME;
                ref_rate_bits += ctx->md_rate_est_ctx
                                     ->comp_ref_fac_bits[svt_av1_get_pred_context_comp_ref_p2(blk_ptr->av1xd)][2][bit2];
                // ref_rate_h = blk_ptr->av1xd->tile_ctx->comp_ref_cdf[pred_context][2];
                // WRITE_REF_BIT(bit2, comp_ref_p2);
            }

            const int bit_bwd    = mbmi->block_mi.ref_frame[1] == ALTREF_FRAME;
            const int pred_ctx_2 = svt_av1_get_pred_context_comp_bwdref_p(blk_ptr->av1xd);
            ref_rate_bits += ctx->md_rate_est_ctx->comp_bwd_ref_fac_bits[pred_ctx_2][0][bit_bwd];
            // ref_rate_i = blk_ptr->av1xd->tile_ctx->comp_bwdref_cdf[pred_ctx_2][0];
            // WRITE_REF_BIT(bit_bwd, comp_bwdref_p);

            if (!bit_bwd) {
                ref_rate_bits += ctx->md_rate_est_ctx->comp_bwd_ref_fac_bits[svt_av1_get_pred_context_comp_bwdref_p1(
                    blk_ptr->av1xd)][1][ref_type[1] == ALTREF2_FRAME];
                // ref_rate_j = blk_ptr->av1xd->tile_ctx->comp_bwdref_cdf[pred_context][1];
                // WRITE_REF_BIT(mbmi->block_mi.ref_frame[1] == ALTREF2_FRAME, comp_bwdref_p1);
            }
        } else {
            const int bit0 = (mbmi->block_mi.ref_frame[0] <= ALTREF_FRAME &&
                              mbmi->block_mi.ref_frame[0] >= BWDREF_FRAME);
            ref_rate_bits += ctx->md_rate_est_ctx
                                 ->single_ref_fac_bits[svt_av1_get_pred_context_single_ref_p1(blk_ptr->av1xd)][0][bit0];
            // ref_rate_k =
            // blk_ptr->av1xd->tile_ctx->single_ref_cdf[svt_av1_get_pred_context_single_ref_p1(blk_ptr->av1xd)][0];
            // WRITE_REF_BIT(bit0, single_ref_p1);

            if (bit0) {
                const int bit1 = mbmi->block_mi.ref_frame[0] == ALTREF_FRAME;
                ref_rate_bits += ctx->md_rate_est_ctx->single_ref_fac_bits[svt_av1_get_pred_context_single_ref_p2(
                    blk_ptr->av1xd)][1][bit1];
                // ref_rate_l =
                // blk_ptr->av1xd->tile_ctx->single_ref_cdf[svt_av1_get_pred_context_single_ref_p2(blk_ptr->av1xd)][1];
                // WRITE_REF_BIT(bit1, single_ref_p2);
                if (!bit1) {
                    ref_rate_bits += ctx->md_rate_est_ctx->single_ref_fac_bits[svt_av1_get_pred_context_single_ref_p6(
                        blk_ptr->av1xd)][5][ref_frame_type == ALTREF2_FRAME];
                    // ref_rate_m =
                    // blk_ptr->av1xd->tile_ctx->single_ref_cdf[svt_av1_get_pred_context_single_ref_p6(blk_ptr->av1xd)][5];
                    // WRITE_REF_BIT(mbmi->block_mi.ref_frame[0] == ALTREF2_FRAME, single_ref_p6);
                }
            } else {
                const int bit2 = (mbmi->block_mi.ref_frame[0] == LAST3_FRAME ||
                                  mbmi->block_mi.ref_frame[0] == GOLDEN_FRAME);
                ref_rate_bits += ctx->md_rate_est_ctx->single_ref_fac_bits[svt_av1_get_pred_context_single_ref_p3(
                    blk_ptr->av1xd)][2][bit2];
                // ref_rate_n =
                // blk_ptr->av1xd->tile_ctx->single_ref_cdf[svt_av1_get_pred_context_single_ref_p3(blk_ptr->av1xd)][2];
                // WRITE_REF_BIT(bit2, single_ref_p3);
                if (!bit2) {
                    const int bit3 = mbmi->block_mi.ref_frame[0] != LAST_FRAME;
                    ref_rate_bits += ctx->md_rate_est_ctx->single_ref_fac_bits[svt_av1_get_pred_context_single_ref_p4(
                        blk_ptr->av1xd)][3][bit3];
                    // ref_rate_o =
                    // blk_ptr->av1xd->tile_ctx->single_ref_cdf[svt_av1_get_pred_context_single_ref_p4(blk_ptr->av1xd)][3];
                    // WRITE_REF_BIT(bit3, single_ref_p4);
                } else {
                    const int bit4 = mbmi->block_mi.ref_frame[0] != LAST3_FRAME;
                    ref_rate_bits += ctx->md_rate_est_ctx->single_ref_fac_bits[svt_av1_get_pred_context_single_ref_p5(
                        blk_ptr->av1xd)][4][bit4];
                    // ref_rate_p =
                    // blk_ptr->av1xd->tile_ctx->single_ref_cdf[svt_av1_get_pred_context_single_ref_p5(blk_ptr->av1xd)][4];
                    // WRITE_REF_BIT(bit4, single_ref_p5);
                }
            }
        }
    }
    return ref_rate_bits;
}

int                    svt_aom_get_comp_group_idx_context_enc(const MacroBlockD *xd);
int                    is_any_masked_compound_used(BlockSize bsize);
static INLINE uint32_t get_compound_mode_rate(struct ModeDecisionContext *ctx, ModeDecisionCandidate *cand,
                                              BlkStruct *blk_ptr, uint8_t ref_frame_type, BlockSize bsize,
                                              SequenceControlSet *scs, PictureControlSet *pcs) {
    uint32_t          comp_rate = 0;
    MbModeInfo *const mbmi      = &blk_ptr->av1xd->mi[0]->mbmi;
    MvReferenceFrame  rf[2];
    av1_set_ref_frame(rf, ref_frame_type);
    mbmi->block_mi.ref_frame[0] = rf[0];
    mbmi->block_mi.ref_frame[1] = rf[1];

    //NOTE  :  Make sure, any cuPtr data is already set before   usage

    if (svt_aom_has_second_ref(mbmi)) {
        const int masked_compound_used = is_any_masked_compound_used(bsize) && scs->seq_header.enable_masked_compound;

        if (masked_compound_used) {
            const int ctx_comp_group_idx = svt_aom_get_comp_group_idx_context_enc(blk_ptr->av1xd);
            comp_rate = ctx->md_rate_est_ctx->comp_group_idx_fac_bits[ctx_comp_group_idx][cand->comp_group_idx];
        } else {
            assert(cand->comp_group_idx == 0);
        }

        if (cand->comp_group_idx == 0) {
            if (cand->compound_idx)
                assert(cand->interinter_comp.type == COMPOUND_AVERAGE);

            if (scs->seq_header.order_hint_info.enable_jnt_comp) {
                const int comp_index_ctx = svt_aom_get_comp_index_context_enc(pcs->ppcs,
                                                                              pcs->ppcs->cur_order_hint,
                                                                              pcs->ppcs->ref_order_hint[rf[0] - 1],
                                                                              pcs->ppcs->ref_order_hint[rf[1] - 1],
                                                                              blk_ptr->av1xd);
                comp_rate += ctx->md_rate_est_ctx->comp_idx_fac_bits[comp_index_ctx][cand->compound_idx];
            } else {
                assert(cand->compound_idx == 1);
            }
        } else {
            assert(pcs->ppcs->frm_hdr.reference_mode != SINGLE_REFERENCE && is_inter_compound_mode(cand->pred_mode));
            assert(masked_compound_used);
            // compound_diffwtd, wedge
            assert(cand->interinter_comp.type == COMPOUND_WEDGE || cand->interinter_comp.type == COMPOUND_DIFFWTD);

            if (is_interinter_compound_used(COMPOUND_WEDGE, bsize))
                comp_rate +=
                    ctx->md_rate_est_ctx->compound_type_fac_bits[bsize][cand->interinter_comp.type - COMPOUND_WEDGE];

            if (cand->interinter_comp.type == COMPOUND_WEDGE) {
                assert(is_interinter_compound_used(COMPOUND_WEDGE, bsize));
                comp_rate += ctx->md_rate_est_ctx->wedge_idx_fac_bits[bsize][cand->interinter_comp.wedge_index];
                comp_rate += av1_cost_literal(1);
            } else {
                assert(cand->interinter_comp.type == COMPOUND_DIFFWTD);
                comp_rate += av1_cost_literal(1);
            }
        }
    }

    return comp_rate;
}
int             svt_aom_is_interintra_wedge_used(BlockSize bsize);
static uint64_t av1_inter_fast_cost_light(struct ModeDecisionContext *ctx, BlkStruct *blk_ptr,
                                          ModeDecisionCandidateBuffer *cand_bf, uint64_t luma_distortion,
                                          uint64_t chroma_distortion, uint64_t lambda, PictureControlSet *pcs,
                                          CandidateMv *ref_mv_stack) {
    ModeDecisionCandidate *cand = cand_bf->cand;
    // NM - fast inter cost estimation
    MdRateEstimationContext *r = ctx->md_rate_est_ctx;
    //_mm_prefetch(p, _MM_HINT_T2);
    // Luma rate
    uint32_t luma_rate   = 0;
    uint32_t chroma_rate = 0;
    uint64_t mv_rate     = 0;
    // Luma and chroma distortion
    uint64_t luma_sad;
    uint64_t chromasad_;
    uint64_t total_distortion;

    uint32_t             rate;
    const PredictionMode inter_mode          = (PredictionMode)cand->pred_mode;
    const uint8_t        have_nearmv         = have_nearmv_in_inter_mode(inter_mode);
    uint64_t             inter_mode_bits_num = 0;
    const uint8_t        skip_mode_ctx       = ctx->skip_mode_ctx;
    MvReferenceFrame     rf[2];
    av1_set_ref_frame(rf, cand->ref_frame_type);
    const uint8_t  is_compound  = is_inter_compound_mode(cand->pred_mode);
    const uint32_t mode_context = svt_aom_mode_context_analyzer(blk_ptr->inter_mode_ctx[cand->ref_frame_type], rf);
    uint64_t       reference_picture_bits_num = 0;
    reference_picture_bits_num                = ctx->estimate_ref_frames_num_bits[cand->ref_frame_type];
    if (is_compound) {
        assert(INTER_COMPOUND_OFFSET(inter_mode) < INTER_COMPOUND_MODES);
        inter_mode_bits_num += r->inter_compound_mode_fac_bits[mode_context][INTER_COMPOUND_OFFSET(inter_mode)];
    } else {
        int16_t newmv_ctx = mode_context & NEWMV_CTX_MASK;
        //aom_write_symbol(ec_writer, mode != NEWMV, frame_context->newmv_cdf[newmv_ctx], 2);
        inter_mode_bits_num += r->new_mv_mode_fac_bits[newmv_ctx][inter_mode != NEWMV];
        if (inter_mode != NEWMV) {
            const int16_t zero_mv_ctx = (mode_context >> GLOBALMV_OFFSET) & GLOBALMV_CTX_MASK;
            //aom_write_symbol(ec_writer, mode != GLOBALMV, frame_context->zeromv_cdf[zero_mv_ctx], 2);
            inter_mode_bits_num += r->zero_mv_mode_fac_bits[zero_mv_ctx][inter_mode != GLOBALMV];
            if (inter_mode != GLOBALMV) {
                int16_t ref_mv_ctx = (mode_context >> REFMV_OFFSET) & REFMV_CTX_MASK;
                /*aom_write_symbol(ec_writer, mode != NEARESTMV, frame_context->refmv_cdf[refmv_ctx], 2);*/
                inter_mode_bits_num += r->ref_mv_mode_fac_bits[ref_mv_ctx][inter_mode != NEARESTMV];
            }
        }
    }
    if (inter_mode == NEWMV || inter_mode == NEW_NEWMV || have_nearmv) {
        //drLIdex cost estimation
        const int32_t new_mv = inter_mode == NEWMV || inter_mode == NEW_NEWMV;
        if (new_mv) {
            int32_t idx;
            for (idx = 0; idx < 2; ++idx) {
                if (blk_ptr->av1xd->ref_mv_count[cand->ref_frame_type] > idx + 1) {
                    uint8_t drl_1_ctx = av1_drl_ctx(ref_mv_stack, idx);
                    inter_mode_bits_num += r->drl_mode_fac_bits[drl_1_ctx][cand->drl_index != idx];
                    if (cand->drl_index == idx)
                        break;
                }
            }
        }
        if (have_nearmv) {
            int32_t idx;
            for (idx = 1; idx < 3; ++idx) {
                if (blk_ptr->av1xd->ref_mv_count[cand->ref_frame_type] > idx + 1) {
                    uint8_t drl_ctx = av1_drl_ctx(ref_mv_stack, idx);
                    inter_mode_bits_num += r->drl_mode_fac_bits[drl_ctx][cand->drl_index != (idx - 1)];
                    if (cand->drl_index == (idx - 1))
                        break;
                }
            }
        }
    }
    if (svt_aom_have_newmv_in_inter_mode(inter_mode)) {
        const uint16_t factor = pcs->ppcs->frm_hdr.allow_screen_content_tools ? 20 : 50;
        if (is_compound) {
            mv_rate = 0;
            if (inter_mode == NEW_NEWMV) {
                for (RefList ref_list_idx = 0; ref_list_idx < 2; ++ref_list_idx) {
                    MV mv = {
                        .row = cand->mv[ref_list_idx].y,
                        .col = cand->mv[ref_list_idx].x,
                    };

                    MV ref_mv = {
                        .row = cand->pred_mv[ref_list_idx].y,
                        .col = cand->pred_mv[ref_list_idx].x,
                    };
                    const uint16_t absmvdiffx = ABS(mv.col - ref_mv.col);
                    const uint16_t absmvdiffy = ABS(mv.row - ref_mv.row);
                    mv_rate += 1296 + (factor * (absmvdiffx + absmvdiffy));
                }
            } else if (inter_mode == NEAREST_NEWMV || inter_mode == NEAR_NEWMV) {
                MV mv = {
                    .row = cand->mv[REF_LIST_1].y,
                    .col = cand->mv[REF_LIST_1].x,
                };

                MV ref_mv = {
                    .row = cand->pred_mv[REF_LIST_1].y,
                    .col = cand->pred_mv[REF_LIST_1].x,
                };
                const uint16_t absmvdiffx = ABS(mv.col - ref_mv.col);
                const uint16_t absmvdiffy = ABS(mv.row - ref_mv.row);
                mv_rate += 1296 + (factor * (absmvdiffx + absmvdiffy));
            } else {
                assert(inter_mode == NEW_NEARESTMV || inter_mode == NEW_NEARMV);
                MV mv = {
                    .row = cand->mv[REF_LIST_0].y,
                    .col = cand->mv[REF_LIST_0].x,
                };

                MV ref_mv = {
                    .row = cand->pred_mv[REF_LIST_0].y,
                    .col = cand->pred_mv[REF_LIST_0].x,
                };
                const uint16_t absmvdiffx = ABS(mv.col - ref_mv.col);
                const uint16_t absmvdiffy = ABS(mv.row - ref_mv.row);
                mv_rate += 1296 + (factor * (absmvdiffx + absmvdiffy));
            }
        } else {
            assert(!is_compound); // single ref inter prediction
            RefList ref_list_idx = get_list_idx(rf[0]);
            MV      mv           = {
                               .row = cand->mv[ref_list_idx].y,
                               .col = cand->mv[ref_list_idx].x,
            };

            MV ref_mv = {
                .row = cand->pred_mv[ref_list_idx].y,
                .col = cand->pred_mv[ref_list_idx].x,
            };
            const uint16_t absmvdiffx = ABS(mv.col - ref_mv.col);
            const uint16_t absmvdiffy = ABS(mv.row - ref_mv.row);
            mv_rate += 1296 + (factor * (absmvdiffx + absmvdiffy));
        }
    }
    // NM - To be added when the overlappable mode is adopted
    //    read_compound_type(is_compound)
    // NM - To be added when switchable filter is adopted
    //    if (interpolation_filter == SWITCHABLE) {
    //        for (dir = 0; dir < (enable_dual_filter ? 2 : 1); dir++) {
    //            if (needs_interp_filter()) {
    //            interp_filter[1] = interp_filter[0]
    //    }
    //    else {
    //        for (dir = 0; dir < 2; dir++)
    //            interp_filter[dir] = interpolation_filter
    //    }
    uint32_t is_inter_rate = r->intra_inter_fac_bits[ctx->is_inter_ctx][1];

    uint32_t skip_mode_rate = pcs->ppcs->frm_hdr.skip_mode_params.skip_mode_flag &&
            is_comp_ref_allowed(ctx->blk_geom->bsize)
        ? r->skip_mode_fac_bits[skip_mode_ctx][0]
        : 0;
    luma_rate = (uint32_t)(reference_picture_bits_num + skip_mode_rate + inter_mode_bits_num + mv_rate + is_inter_rate);
    //chroma_rate = intra_chroma_mode_bits_num + intra_chroma_ang_mode_bits_num;

    // Keep the Fast Luma and Chroma rate for future use
    cand_bf->fast_luma_rate   = luma_rate * RATE_WEIGHT;
    cand_bf->fast_chroma_rate = chroma_rate * RATE_WEIGHT;
    luma_sad                  = (LUMA_WEIGHT * luma_distortion) << AV1_COST_PRECISION;
    chromasad_                = (CHROMA_WEIGHT * chroma_distortion) << AV1_COST_PRECISION;
    total_distortion          = luma_sad + chromasad_;
    //if (blk_geom->has_uv == 0 && chromasad_ != 0)
    //    SVT_LOG("svt_aom_inter_fast_cost: Chroma error");
    rate = luma_rate + chroma_rate;
    // Assign fast cost
    if (cand->skip_mode_allowed) {
        skip_mode_rate = r->skip_mode_fac_bits[skip_mode_ctx][1];
        if (skip_mode_rate < rate)
            return (RDCOST(lambda, skip_mode_rate, total_distortion));
    }
    return (RDCOST(lambda, rate, total_distortion));
}
uint64_t svt_aom_inter_fast_cost(PictureControlSet *pcs, struct ModeDecisionContext *ctx,
                                 ModeDecisionCandidateBuffer *cand_bf, uint64_t lambda, uint64_t luma_distortion,
                                 uint64_t chroma_distortion) {
    const BlockGeom       *blk_geom     = ctx->blk_geom;
    BlkStruct             *blk_ptr      = ctx->blk_ptr;
    ModeDecisionCandidate *cand         = cand_bf->cand;
    CandidateMv           *ref_mv_stack = &(
        ctx->md_local_blk_unit[blk_geom->blkidx_mds].ed_ref_mv_stack[cand->ref_frame_type][0]);

    if (ctx->approx_inter_rate)
        return av1_inter_fast_cost_light(
            ctx, blk_ptr, cand_bf, luma_distortion, chroma_distortion, lambda, pcs, ref_mv_stack);
    FrameHeader *frm_hdr = &pcs->ppcs->frm_hdr;

    // Luma rate
    uint32_t luma_rate   = 0;
    uint32_t chroma_rate = 0;
    uint64_t mv_rate     = 0;
    // Luma and chroma distortion
    uint64_t luma_sad;
    uint64_t chromasad_;
    uint64_t total_distortion;

    uint32_t       rate;
    PredictionMode inter_mode = (PredictionMode)cand->pred_mode;

    uint64_t inter_mode_bits_num = 0;

    const uint8_t    skip_mode_ctx = ctx->skip_mode_ctx;
    MvReferenceFrame rf[2];
    av1_set_ref_frame(rf, cand->ref_frame_type);
    const uint8_t is_compound  = is_inter_compound_mode(cand->pred_mode);
    uint32_t      mode_context = svt_aom_mode_context_analyzer(blk_ptr->inter_mode_ctx[cand->ref_frame_type], rf);
    uint64_t      reference_picture_bits_num = 0;

    //Reference Type and Mode Bit estimation
    reference_picture_bits_num = ctx->estimate_ref_frames_num_bits[cand->ref_frame_type];
    if (is_compound) {
        assert(INTER_COMPOUND_OFFSET(inter_mode) < INTER_COMPOUND_MODES);
        inter_mode_bits_num +=
            ctx->md_rate_est_ctx->inter_compound_mode_fac_bits[mode_context][INTER_COMPOUND_OFFSET(inter_mode)];
    } else {
        // uint32_t newmv_ctx = mode_context & NEWMV_CTX_MASK;
        // inter_mode_bits_num = cand_bf->cand->md_rate_est_ctx->new_mv_mode_fac_bits[mode_ctx][0];

        int16_t newmv_ctx = mode_context & NEWMV_CTX_MASK;
        // aom_write_symbol(ec_writer, mode != NEWMV, frame_context->newmv_cdf[newmv_ctx], 2);
        inter_mode_bits_num += ctx->md_rate_est_ctx->new_mv_mode_fac_bits[newmv_ctx][inter_mode != NEWMV];
        if (inter_mode != NEWMV) {
            const int16_t zero_mv_ctx = (mode_context >> GLOBALMV_OFFSET) & GLOBALMV_CTX_MASK;
            // aom_write_symbol(ec_writer, mode != GLOBALMV, frame_context->zeromv_cdf[zero_mv_ctx],
            // 2);
            inter_mode_bits_num += ctx->md_rate_est_ctx->zero_mv_mode_fac_bits[zero_mv_ctx][inter_mode != GLOBALMV];
            if (inter_mode != GLOBALMV) {
                int16_t ref_mv_ctx = (mode_context >> REFMV_OFFSET) & REFMV_CTX_MASK;
                /*aom_write_symbol(ec_writer, mode != NEARESTMV,
                 * frame_context->refmv_cdf[refmv_ctx], 2);*/
                inter_mode_bits_num += ctx->md_rate_est_ctx->ref_mv_mode_fac_bits[ref_mv_ctx][inter_mode != NEARESTMV];
            }
        }
    }
    if (inter_mode == NEWMV || inter_mode == NEW_NEWMV || have_nearmv_in_inter_mode(inter_mode)) {
        //drLIdex cost estimation
        const int32_t new_mv = inter_mode == NEWMV || inter_mode == NEW_NEWMV;
        if (new_mv) {
            int32_t idx;
            for (idx = 0; idx < 2; ++idx) {
                if (blk_ptr->av1xd->ref_mv_count[cand->ref_frame_type] > idx + 1) {
                    uint8_t drl_1_ctx = av1_drl_ctx(ref_mv_stack, idx);
                    inter_mode_bits_num += ctx->md_rate_est_ctx->drl_mode_fac_bits[drl_1_ctx][cand->drl_index != idx];
                    if (cand->drl_index == idx)
                        break;
                }
            }
        }

        if (have_nearmv_in_inter_mode(inter_mode)) {
            int32_t idx;
            for (idx = 1; idx < 3; ++idx) {
                if (blk_ptr->av1xd->ref_mv_count[cand->ref_frame_type] > idx + 1) {
                    uint8_t drl_ctx = av1_drl_ctx(ref_mv_stack, idx);
                    inter_mode_bits_num +=
                        ctx->md_rate_est_ctx->drl_mode_fac_bits[drl_ctx][cand->drl_index != (idx - 1)];

                    if (cand->drl_index == (idx - 1))
                        break;
                }
            }
        }
    }

    if (svt_aom_have_newmv_in_inter_mode(inter_mode)) {
        if (is_compound) {
            mv_rate = 0;

            if (inter_mode == NEW_NEWMV) {
                for (RefList ref_list_idx = 0; ref_list_idx < 2; ++ref_list_idx) {
                    MV mv = {
                        .row = cand->mv[ref_list_idx].y,
                        .col = cand->mv[ref_list_idx].x,
                    };

                    MV ref_mv = {
                        .row = cand->pred_mv[ref_list_idx].y,
                        .col = cand->pred_mv[ref_list_idx].x,
                    };
                    mv_rate += svt_av1_mv_bit_cost(&mv,
                                                   &ref_mv,
                                                   ctx->md_rate_est_ctx->nmv_vec_cost,
                                                   ctx->md_rate_est_ctx->nmvcoststack,
                                                   MV_COST_WEIGHT);
                }
            } else if (inter_mode == NEAREST_NEWMV || inter_mode == NEAR_NEWMV) {
                MV mv = {
                    .row = cand->mv[REF_LIST_1].y,
                    .col = cand->mv[REF_LIST_1].x,
                };

                MV ref_mv = {
                    .row = cand->pred_mv[REF_LIST_1].y,
                    .col = cand->pred_mv[REF_LIST_1].x,
                };
                mv_rate += svt_av1_mv_bit_cost(&mv,
                                               &ref_mv,
                                               ctx->md_rate_est_ctx->nmv_vec_cost,
                                               ctx->md_rate_est_ctx->nmvcoststack,
                                               MV_COST_WEIGHT);
            } else {
                assert(inter_mode == NEW_NEARESTMV || inter_mode == NEW_NEARMV);
                MV mv = {
                    .row = cand->mv[REF_LIST_0].y,
                    .col = cand->mv[REF_LIST_0].x,
                };

                MV ref_mv = {
                    .row = cand->pred_mv[REF_LIST_0].y,
                    .col = cand->pred_mv[REF_LIST_0].x,
                };
                mv_rate += svt_av1_mv_bit_cost(&mv,
                                               &ref_mv,
                                               ctx->md_rate_est_ctx->nmv_vec_cost,
                                               ctx->md_rate_est_ctx->nmvcoststack,
                                               MV_COST_WEIGHT);
            }
        } else {
            assert(!is_compound); // single ref inter prediction
            RefList ref_list_idx = get_list_idx(rf[0]);
            MV      mv           = {
                               .row = cand->mv[ref_list_idx].y,
                               .col = cand->mv[ref_list_idx].x,
            };

            MV ref_mv = {
                .row = cand->pred_mv[ref_list_idx].y,
                .col = cand->pred_mv[ref_list_idx].x,
            };
            mv_rate = svt_av1_mv_bit_cost(
                &mv, &ref_mv, ctx->md_rate_est_ctx->nmv_vec_cost, ctx->md_rate_est_ctx->nmvcoststack, MV_COST_WEIGHT);
        }
    }
    // inter intra mode rate
    if (pcs->ppcs->frm_hdr.reference_mode != COMPOUND_REFERENCE &&
        pcs->ppcs->scs->seq_header.enable_interintra_compound &&
        /* Check if inter-intra is allowed for current block size / mode (even if the feature is off
        * for the current block, we still need to signal inter-intra off.
        */
        svt_is_interintra_allowed(TRUE, blk_geom->bsize, cand->pred_mode, rf)) {
        const int interintra  = cand->is_interintra_used;
        const int bsize_group = size_group_lookup[blk_geom->bsize];

        inter_mode_bits_num += ctx->md_rate_est_ctx->inter_intra_fac_bits[bsize_group][cand->is_interintra_used];

        if (interintra) {
            inter_mode_bits_num += ctx->md_rate_est_ctx->inter_intra_mode_fac_bits[bsize_group][cand->interintra_mode];

            if (svt_aom_is_interintra_wedge_used(blk_geom->bsize)) {
                inter_mode_bits_num +=
                    ctx->md_rate_est_ctx->wedge_inter_intra_fac_bits[blk_geom->bsize][cand->use_wedge_interintra];

                if (cand->use_wedge_interintra) {
                    inter_mode_bits_num +=
                        ctx->md_rate_est_ctx->wedge_idx_fac_bits[blk_geom->bsize][cand->interintra_wedge_index];
                }
            }
        }
    }
    Bool is_inter = inter_mode >= SINGLE_INTER_MODE_START && inter_mode < SINGLE_INTER_MODE_END;
    if (is_inter && frm_hdr->is_motion_mode_switchable && rf[1] != INTRA_FRAME) {
        MotionMode motion_mode_rd                      = cand->motion_mode;
        BlockSize  bsize                               = blk_geom->bsize;
        blk_ptr->prediction_unit_array[0].num_proj_ref = cand->num_proj_ref;
        MotionMode last_motion_mode_allowed            = svt_aom_motion_mode_allowed(
            pcs,
            blk_ptr->prediction_unit_array[0].num_proj_ref,
            blk_ptr->prediction_unit_array[0].overlappable_neighbors,
            bsize,
            rf[0],
            rf[1],
            inter_mode);
        switch (last_motion_mode_allowed) {
        case SIMPLE_TRANSLATION: break;
        case OBMC_CAUSAL:
            inter_mode_bits_num += ctx->md_rate_est_ctx->motion_mode_fac_bits1[bsize][motion_mode_rd == OBMC_CAUSAL];
            break;
        default: inter_mode_bits_num += ctx->md_rate_est_ctx->motion_mode_fac_bits[bsize][motion_mode_rd];
        }
    }
    // this func return 0 if masked=0 and distance=0
    inter_mode_bits_num += get_compound_mode_rate(
        ctx, cand, blk_ptr, cand->ref_frame_type, blk_geom->bsize, pcs->ppcs->scs, pcs);
    // NM - To be added when the overlappable mode is adopted
    //    read_compound_type(is_compound)
    // NM - To be added when switchable filter is adopted
    //    if (interpolation_filter == SWITCHABLE) {
    //        for (dir = 0; dir < (enable_dual_filter ? 2 : 1); dir++) {
    //            if (needs_interp_filter()) {
    //                interp_filter[dir]    S()
    //            }
    //            else {
    //                interp_filter[dir] = EIGHTTAP
    //            }
    //        }
    //        if (!enable_dual_filter)
    //            interp_filter[1] = interp_filter[0]
    //    }
    //    else {
    //        for (dir = 0; dir < 2; dir++)
    //            interp_filter[dir] = interpolation_filter
    //    }
    uint32_t is_inter_rate  = ctx->md_rate_est_ctx->intra_inter_fac_bits[ctx->is_inter_ctx][1];
    uint32_t skip_mode_rate = pcs->ppcs->frm_hdr.skip_mode_params.skip_mode_flag && is_comp_ref_allowed(blk_geom->bsize)
        ? ctx->md_rate_est_ctx->skip_mode_fac_bits[skip_mode_ctx][0]
        : 0;
    luma_rate = (uint32_t)(reference_picture_bits_num + skip_mode_rate + inter_mode_bits_num + mv_rate + is_inter_rate);

    // chroma_rate = intra_chroma_mode_bits_num + intra_chroma_ang_mode_bits_num;

    // Keep the Fast Luma and Chroma rate for future use
    cand_bf->fast_luma_rate   = luma_rate * RATE_WEIGHT;
    cand_bf->fast_chroma_rate = chroma_rate * RATE_WEIGHT;
    luma_sad                  = (LUMA_WEIGHT * luma_distortion) << AV1_COST_PRECISION;
    chromasad_                = (CHROMA_WEIGHT * chroma_distortion) << AV1_COST_PRECISION;
    total_distortion          = luma_sad + chromasad_;
    if (blk_geom->has_uv == 0 && chromasad_ != 0)
        SVT_ERROR("svt_aom_inter_fast_cost: Chroma error");
    rate = luma_rate + chroma_rate;
    // Assign fast cost
    if (cand->skip_mode_allowed) {
        skip_mode_rate = ctx->md_rate_est_ctx->skip_mode_fac_bits[skip_mode_ctx][1];
        if (skip_mode_rate < rate)
            return (RDCOST(lambda, skip_mode_rate, total_distortion));
    }
    return (RDCOST(lambda, rate, total_distortion));
}
/*
 */
EbErrorType svt_aom_txb_estimate_coeff_bits_light_pd0(struct ModeDecisionContext         *ctx,
                                                      struct ModeDecisionCandidateBuffer *cand_bf,
                                                      uint32_t txb_origin_index, EbPictureBufferDesc *coeff_buffer_sb,
                                                      uint32_t y_eob, uint64_t *y_txb_coeff_bits, TxSize txsize) {
    if (y_eob) {
        *y_txb_coeff_bits = svt_av1_cost_coeffs_txb(
            ctx,
            0,
            0,
            cand_bf,
            (int32_t *)&coeff_buffer_sb->buffer_y[txb_origin_index * sizeof(int32_t)],
            (uint16_t)y_eob,
            PLANE_TYPE_Y,
            txsize,
            DCT_DCT,
            0,
            0,
            0);

        *y_txb_coeff_bits = (*y_txb_coeff_bits) << ctx->mds_subres_step;

    } else {
        *y_txb_coeff_bits = av1_cost_skip_txb(ctx, 0, 0, txsize, PLANE_TYPE_Y, 0);
    }

    return EB_ErrorNone;
}
EbErrorType svt_aom_txb_estimate_coeff_bits(struct ModeDecisionContext *ctx, uint8_t allow_update_cdf,
                                            FRAME_CONTEXT *ec_ctx, PictureControlSet *pcs,
                                            struct ModeDecisionCandidateBuffer *cand_bf, uint32_t txb_origin_index,
                                            uint32_t txb_chroma_origin_index, EbPictureBufferDesc *coeff_buffer_sb,
                                            uint32_t y_eob, uint32_t cb_eob, uint32_t cr_eob,
                                            uint64_t *y_txb_coeff_bits, uint64_t *cb_txb_coeff_bits,
                                            uint64_t *cr_txb_coeff_bits, TxSize txsize, TxSize txsize_uv,
                                            TxType tx_type, TxType tx_type_uv, COMPONENT_TYPE component_type) {
    EbErrorType return_error = EB_ErrorNone;

    FrameHeader *frm_hdr = &pcs->ppcs->frm_hdr;

    int32_t *coeff_buffer;
    int16_t  luma_txb_skip_context = ctx->luma_txb_skip_context;
    int16_t  luma_dc_sign_context  = ctx->luma_dc_sign_context;
    int16_t  cb_txb_skip_context   = ctx->cb_txb_skip_context;
    int16_t  cb_dc_sign_context    = ctx->cb_dc_sign_context;
    int16_t  cr_txb_skip_context   = ctx->cr_txb_skip_context;
    int16_t  cr_dc_sign_context    = ctx->cr_dc_sign_context;

    Bool reduced_transform_set_flag = frm_hdr->reduced_tx_set ? TRUE : FALSE;

    //Estimate the rate of the transform type and coefficient for Luma

    if (component_type == COMPONENT_LUMA || component_type == COMPONENT_ALL) {
        if (y_eob) {
            coeff_buffer = (int32_t *)&coeff_buffer_sb->buffer_y[txb_origin_index * sizeof(int32_t)];

            *y_txb_coeff_bits = svt_av1_cost_coeffs_txb(ctx,
                                                        allow_update_cdf,
                                                        ec_ctx,
                                                        cand_bf,
                                                        coeff_buffer,
                                                        (uint16_t)y_eob,
                                                        PLANE_TYPE_Y,
                                                        txsize,
                                                        tx_type,
                                                        luma_txb_skip_context,
                                                        luma_dc_sign_context,
                                                        reduced_transform_set_flag);
            *y_txb_coeff_bits = (*y_txb_coeff_bits) << ctx->mds_subres_step;
        } else {
            *y_txb_coeff_bits = av1_cost_skip_txb(
                ctx, allow_update_cdf, ec_ctx, txsize, PLANE_TYPE_Y, luma_txb_skip_context);
        }
    }
    // Estimate the rate of the transform type and coefficient for chroma Cb

    if (component_type == COMPONENT_CHROMA_CB || component_type == COMPONENT_CHROMA ||
        component_type == COMPONENT_ALL) {
        if (cb_eob) {
            coeff_buffer = (int32_t *)&coeff_buffer_sb->buffer_cb[txb_chroma_origin_index * sizeof(int32_t)];

            *cb_txb_coeff_bits = svt_av1_cost_coeffs_txb(ctx,
                                                         allow_update_cdf,
                                                         ec_ctx,
                                                         cand_bf,
                                                         coeff_buffer,
                                                         (uint16_t)cb_eob,
                                                         PLANE_TYPE_UV,
                                                         txsize_uv,
                                                         tx_type_uv,
                                                         cb_txb_skip_context,
                                                         cb_dc_sign_context,
                                                         reduced_transform_set_flag);
        } else {
            *cb_txb_coeff_bits = av1_cost_skip_txb(
                ctx, allow_update_cdf, ec_ctx, txsize_uv, PLANE_TYPE_UV, cb_txb_skip_context);
        }
    }

    if (component_type == COMPONENT_CHROMA_CR || component_type == COMPONENT_CHROMA ||
        component_type == COMPONENT_ALL) {
        //Estimate the rate of the transform type and coefficient for chroma Cr
        if (cr_eob) {
            coeff_buffer = (int32_t *)&coeff_buffer_sb->buffer_cr[txb_chroma_origin_index * sizeof(int32_t)];

            *cr_txb_coeff_bits = svt_av1_cost_coeffs_txb(ctx,
                                                         allow_update_cdf,
                                                         ec_ctx,
                                                         cand_bf,
                                                         coeff_buffer,
                                                         (uint16_t)cr_eob,
                                                         PLANE_TYPE_UV,
                                                         txsize_uv,
                                                         tx_type_uv,
                                                         cr_txb_skip_context,
                                                         cr_dc_sign_context,
                                                         reduced_transform_set_flag);
        } else {
            *cr_txb_coeff_bits = av1_cost_skip_txb(
                ctx, allow_update_cdf, ec_ctx, txsize_uv, PLANE_TYPE_UV, cr_txb_skip_context);
        }
    }

    return return_error;
}

EbErrorType svt_aom_full_cost_light_pd0(ModeDecisionContext *ctx, struct ModeDecisionCandidateBuffer *cand_bf,
                                        uint64_t *y_distortion, uint64_t lambda, uint64_t *y_coeff_bits) {
    EbErrorType return_error = EB_ErrorNone;

    uint64_t coeff_rate = (*y_coeff_bits + (uint64_t)ctx->md_rate_est_ctx->skip_fac_bits[0][0]);

    // Assign full cost
    // Use context index 0 for the partition rate as an approximation to skip call to
    // av1_partition_rate_cost Partition cost is only needed for > 4x4 blocks, but light-PD0 assumes
    // 4x4 blocks are disallowed
    *(cand_bf->full_cost) = RDCOST(
        lambda, coeff_rate + ctx->md_rate_est_ctx->partition_fac_bits[0][PARTITION_NONE], y_distortion[0]);
    return return_error;
}
/*********************************************************************************
 * svt_aom_av1_full_cost function is used to estimate the cost of a candidate mode
 * for full mode decision module.
 **********************************************************************************/
void svt_aom_full_cost(PictureControlSet *pcs, ModeDecisionContext *ctx, struct ModeDecisionCandidateBuffer *cand_bf,
                       uint64_t lambda, uint64_t y_distortion[DIST_TOTAL][DIST_CALC_TOTAL],
                       uint64_t cb_distortion[DIST_TOTAL][DIST_CALC_TOTAL],
                       uint64_t cr_distortion[DIST_TOTAL][DIST_CALC_TOTAL], uint64_t *y_coeff_bits,
                       uint64_t *cb_coeff_bits, uint64_t *cr_coeff_bits) {
    const uint8_t skip_coeff_ctx        = ctx->skip_coeff_ctx;
    const bool    update_full_cost_ssim = (ctx->tune_ssim_level > SSIM_LVL_0 || pcs->scs->static_config.tune == 3) ? true : false;

    // Get the TX size rate for skip and non-skip block. Need both to make non-skip decision
    uint64_t non_skip_tx_size_bits = 0, skip_tx_size_bits = 0;
    if (!ctx->shut_fast_rate && pcs->ppcs->frm_hdr.tx_mode == TX_MODE_SELECT) {
        if (cand_bf->block_has_coeff) {
            non_skip_tx_size_bits = svt_aom_get_tx_size_bits(
                cand_bf, ctx, pcs, cand_bf->cand->tx_depth, /*cand_bf->block_has_coeff*/ 1);
        }

        skip_tx_size_bits = svt_aom_get_tx_size_bits(
            cand_bf, ctx, pcs, cand_bf->cand->tx_depth, /*cand_bf->block_has_coeff*/ 0);
    }

    assert(IMPLIES(is_inter_mode(cand_bf->cand->pred_mode), skip_tx_size_bits == 0));

    // Decide if block should be signalled as skip (send no coeffs)
    if (ctx->blk_skip_decision && cand_bf->block_has_coeff && is_inter_mode(cand_bf->cand->pred_mode)) {
        const uint64_t non_skip_cost = RDCOST(
            lambda,
            (*y_coeff_bits + *cb_coeff_bits + *cr_coeff_bits + non_skip_tx_size_bits +
             (uint64_t)ctx->md_rate_est_ctx->skip_fac_bits[skip_coeff_ctx][0]),
            (y_distortion[DIST_SSD][0] + cb_distortion[DIST_SSD][0] + cr_distortion[DIST_SSD][0]));

        const uint64_t skip_cost = RDCOST(
            lambda,
            ((uint64_t)ctx->md_rate_est_ctx->skip_fac_bits[skip_coeff_ctx][1]) + skip_tx_size_bits,
            (y_distortion[DIST_SSD][1] + cb_distortion[DIST_SSD][1] + cr_distortion[DIST_SSD][1]));

        // Update signals to correspond to skip_mode values (no coeffs, etc.)
        if (skip_cost < non_skip_cost) {
            y_distortion[DIST_SSD][0]  = y_distortion[DIST_SSD][1];
            cb_distortion[DIST_SSD][0] = cb_distortion[DIST_SSD][1];
            cr_distortion[DIST_SSD][0] = cr_distortion[DIST_SSD][1];

            y_distortion[DIST_SSIM][0]  = y_distortion[DIST_SSIM][1];
            cb_distortion[DIST_SSIM][0] = cb_distortion[DIST_SSIM][1];
            cr_distortion[DIST_SSIM][0] = cr_distortion[DIST_SSIM][1];
            cand_bf->block_has_coeff    = 0;
            cand_bf->y_has_coeff        = 0;
            cand_bf->u_has_coeff        = 0;
            cand_bf->v_has_coeff        = 0;
            cand_bf->cnt_nz_coeff       = 0;

            // For inter modes, signalling skip means no TX depth is used and the TX type will be DCT_DCT
            cand_bf->cand->tx_depth          = 0;
            cand_bf->cand->transform_type_uv = DCT_DCT;
            assert(DCT_DCT == 0);
            memset(cand_bf->cand->transform_type, DCT_DCT, 16 * sizeof(cand_bf->cand->transform_type[0]));
            memset(cand_bf->quantized_dc, 0, 3 * 16 * sizeof(cand_bf->quantized_dc[0][0]));
            memset(cand_bf->eob, 0, 3 * 16 * sizeof(cand_bf->eob[0][0]));
        }
    }

    uint64_t coeff_rate = 0;
    if (cand_bf->block_has_coeff) {
        coeff_rate = (*y_coeff_bits + *cb_coeff_bits + *cr_coeff_bits + non_skip_tx_size_bits +
                      (uint64_t)ctx->md_rate_est_ctx->skip_fac_bits[skip_coeff_ctx][0]);
    } else {
        coeff_rate = ctx->md_rate_est_ctx->skip_fac_bits[skip_coeff_ctx][1] + skip_tx_size_bits;
    }

    uint64_t mode_rate            = cand_bf->fast_luma_rate + cand_bf->fast_chroma_rate + coeff_rate;
    uint64_t mode_distortion      = y_distortion[DIST_SSD][0] + cb_distortion[DIST_SSD][0] + cr_distortion[DIST_SSD][0];
    uint64_t mode_ssim_distortion = update_full_cost_ssim
        ? y_distortion[DIST_SSIM][0] + cb_distortion[DIST_SSIM][0] + cr_distortion[DIST_SSIM][0]
        : 0;
    uint64_t mode_cost            = RDCOST(lambda, mode_rate, mode_distortion);

    // If skip_mode is allowed for this candidate, check cost of skip mode compared to regular cost
    if (cand_bf->cand->skip_mode_allowed == TRUE) {
        const uint8_t skip_mode_ctx = ctx->skip_mode_ctx;

        // Skip mode cost
        const uint64_t skip_mode_rate       = ctx->md_rate_est_ctx->skip_mode_fac_bits[skip_mode_ctx][1];
        const uint64_t skip_mode_distortion = y_distortion[DIST_SSD][1] + cb_distortion[DIST_SSD][1] +
            cr_distortion[DIST_SSD][1];
        const uint64_t skip_mode_ssim_distortion = update_full_cost_ssim
            ? y_distortion[DIST_SSIM][1] + cb_distortion[DIST_SSIM][1] + cr_distortion[DIST_SSIM][1]
            : 0;
        const uint64_t skip_mode_cost            = RDCOST(lambda, skip_mode_rate, skip_mode_distortion);

        cand_bf->cand->skip_mode = FALSE;
        if (skip_mode_cost <= mode_cost) {
            // Update candidate cost
            mode_cost                = skip_mode_cost;
            mode_rate                = skip_mode_rate;
            mode_distortion          = skip_mode_distortion;
            mode_ssim_distortion     = skip_mode_ssim_distortion;
            cand_bf->cand->skip_mode = TRUE;

            // Update signals to correspond to skip_mode values (no coeffs, etc.)
            cand_bf->block_has_coeff = 0;
            cand_bf->y_has_coeff     = 0;
            cand_bf->u_has_coeff     = 0;
            cand_bf->v_has_coeff     = 0;
            cand_bf->cnt_nz_coeff    = 0;
            cand_bf->cand->tx_depth  = 0;
            assert(DCT_DCT == 0);
            memset(cand_bf->cand->transform_type, DCT_DCT, 16 * sizeof(cand_bf->cand->transform_type[0]));
            cand_bf->cand->transform_type_uv = DCT_DCT;
            memset(cand_bf->quantized_dc, 0, 3 * 16 * sizeof(cand_bf->quantized_dc[0][0]));
            memset(cand_bf->eob, 0, 3 * 16 * sizeof(cand_bf->eob[0][0]));
        }
    }

    // Assign full cost
    *(cand_bf->full_cost) = mode_cost;
    cand_bf->total_rate   = mode_rate;
    cand_bf->full_dist    = (uint32_t)mode_distortion;
    if (update_full_cost_ssim) {
        assert(ctx->pd_pass == PD_PASS_1);
        assert(ctx->md_stage == MD_STAGE_3);
        *(cand_bf->full_cost_ssim) = RDCOST(lambda, mode_rate, mode_ssim_distortion);
    }
    return;
}

/************************************************************
 * Coding Loop Context Generation
 ************************************************************/
void svt_aom_coding_loop_context_generation(PictureControlSet *pcs, ModeDecisionContext *ctx) {
    BlkStruct         *blk_ptr                       = ctx->blk_ptr;
    uint32_t           blk_org_x                     = ctx->blk_org_x;
    uint32_t           blk_org_y                     = ctx->blk_org_y;
    NeighborArrayUnit *leaf_partition_na             = ctx->leaf_partition_na;
    uint32_t           partition_left_neighbor_index = get_neighbor_array_unit_left_index(leaf_partition_na, blk_org_y);
    uint32_t           partition_above_neighbor_index = get_neighbor_array_unit_top_index(leaf_partition_na, blk_org_x);
    MacroBlockD       *xd                             = blk_ptr->av1xd;
    if (!ctx->shut_fast_rate) {
        if (pcs->slice_type == I_SLICE) {
            svt_aom_get_kf_y_mode_ctx(xd, &ctx->intra_luma_top_ctx, &ctx->intra_luma_left_ctx);
        }
        ctx->is_inter_ctx  = svt_av1_get_intra_inter_context(xd);
        ctx->skip_mode_ctx = av1_get_skip_mode_context(xd);
    }
    // Generate Partition context
    ctx->md_local_blk_unit[blk_ptr->mds_idx].above_neighbor_partition =
        (((PartitionContext *)leaf_partition_na->top_array)[partition_above_neighbor_index].above ==
         (char)INVALID_NEIGHBOR_DATA)
        ? 0
        : ((PartitionContext *)leaf_partition_na->top_array)[partition_above_neighbor_index].above;

    ctx->md_local_blk_unit[blk_ptr->mds_idx].left_neighbor_partition =
        (((PartitionContext *)leaf_partition_na->left_array)[partition_left_neighbor_index].left ==
         (char)INVALID_NEIGHBOR_DATA)
        ? 0
        : ((PartitionContext *)leaf_partition_na->left_array)[partition_left_neighbor_index].left;

    // Collect Neighbor ref cout
    if (pcs->slice_type != I_SLICE || pcs->ppcs->frm_hdr.allow_intrabc)
        svt_aom_collect_neighbors_ref_counts_new(blk_ptr->av1xd);

    // Skip Coeff Context
    ctx->skip_coeff_ctx = ctx->rate_est_ctrls.update_skip_coeff_ctx ? av1_get_skip_context(xd) : 0;
}

static INLINE int block_signals_txsize(BlockSize bsize) { return bsize > BLOCK_4X4; }

static INLINE int get_vartx_max_txsize(/*const MbModeInfo *xd,*/ BlockSize bsize, int plane) {
    /* if (xd->lossless[xd->mi[0]->segment_id]) return TX_4X4;*/
    const TxSize max_txsize = max_txsize_rect_lookup[bsize];
    if (plane == 0)
        return max_txsize; // luma
    return av1_get_adjusted_tx_size(max_txsize); // chroma
}

static INLINE int max_block_wide(const MacroBlockD *xd, BlockSize bsize, int plane) {
    int max_blocks_wide = block_size_wide[bsize];

    if (xd->mb_to_right_edge < 0)
        max_blocks_wide += gcc_right_shift(xd->mb_to_right_edge, 3 + !!plane);

    // Scale the width in the transform block unit.
    return max_blocks_wide >> tx_size_wide_log2[0];
}

static INLINE int max_block_high(const MacroBlockD *xd, BlockSize bsize, int plane) {
    int max_blocks_high = block_size_high[bsize];

    if (xd->mb_to_bottom_edge < 0)
        max_blocks_high += gcc_right_shift(xd->mb_to_bottom_edge, 3 + !!plane);

    // Scale the height in the transform block unit.
    return max_blocks_high >> tx_size_high_log2[0];
}

static INLINE void txfm_partition_update(TXFM_CONTEXT *above_ctx, TXFM_CONTEXT *left_ctx, TxSize tx_size,
                                         TxSize txb_size) {
    BlockSize bsize = txsize_to_bsize[txb_size];
    assert(bsize < BlockSizeS_ALL);
    int     bh  = mi_size_high[bsize];
    int     bw  = mi_size_wide[bsize];
    uint8_t txw = tx_size_wide[tx_size];
    uint8_t txh = tx_size_high[tx_size];
    int     i;
    for (i = 0; i < bh; ++i) left_ctx[i] = txh;
    for (i = 0; i < bw; ++i) above_ctx[i] = txw;
}

static INLINE TxSize get_sqr_tx_size(int tx_dim) {
    switch (tx_dim) {
    case 128:
    case 64: return TX_64X64; break;
    case 32: return TX_32X32; break;
    case 16: return TX_16X16; break;
    case 8: return TX_8X8; break;
    default: return TX_4X4;
    }
}
static INLINE int txfm_partition_context(TXFM_CONTEXT *above_ctx, TXFM_CONTEXT *left_ctx, BlockSize bsize,
                                         TxSize tx_size) {
    const uint8_t txw      = tx_size_wide[tx_size];
    const uint8_t txh      = tx_size_high[tx_size];
    const int     above    = *above_ctx < txw;
    const int     left     = *left_ctx < txh;
    int           category = TXFM_PARTITION_CONTEXTS;

    // dummy return, not used by others.
    if (tx_size == TX_4X4)
        return 0;

    TxSize max_tx_size = get_sqr_tx_size(AOMMAX(block_size_wide[bsize], block_size_high[bsize]));

    if (max_tx_size >= TX_8X8) {
        category = (txsize_sqr_up_map[tx_size] != max_tx_size && max_tx_size > TX_8X8) +
            (TX_SIZES - 1 - max_tx_size) * 2;
    }
    assert(category != TXFM_PARTITION_CONTEXTS);
    return category * 3 + above + left;
}

static uint64_t cost_tx_size_vartx(MacroBlockD *xd, const MbModeInfo *mbmi, TxSize tx_size, int depth, int blk_row,
                                   int blk_col, MdRateEstimationContext *md_rate_est_ctx, FRAME_CONTEXT *ec_ctx,
                                   uint8_t allow_update_cdf) {
    uint64_t  bits            = 0;
    const int max_blocks_high = max_block_high(xd, mbmi->block_mi.bsize, 0);
    const int max_blocks_wide = max_block_wide(xd, mbmi->block_mi.bsize, 0);

    if (blk_row >= max_blocks_high || blk_col >= max_blocks_wide)
        return bits;

    if (depth == MAX_VARTX_DEPTH) {
        txfm_partition_update(xd->above_txfm_context + blk_col, xd->left_txfm_context + blk_row, tx_size, tx_size);

        return bits;
    }

    const int ctx = txfm_partition_context(
        xd->above_txfm_context + blk_col, xd->left_txfm_context + blk_row, mbmi->block_mi.bsize, tx_size);
    const int write_txfm_partition = (tx_size == tx_depth_to_tx_size[mbmi->block_mi.tx_depth][mbmi->block_mi.bsize]);
    if (write_txfm_partition) {
        bits += md_rate_est_ctx->txfm_partition_fac_bits[ctx][0];

        if (allow_update_cdf)
            update_cdf(ec_ctx->txfm_partition_cdf[ctx], 0, 2);

        txfm_partition_update(xd->above_txfm_context + blk_col, xd->left_txfm_context + blk_row, tx_size, tx_size);

    } else {
        assert(tx_size < TX_SIZES_ALL);
        const TxSize sub_txs = sub_tx_size_map[tx_size];
        const int    bsw     = tx_size_wide_unit[sub_txs];
        const int    bsh     = tx_size_high_unit[sub_txs];

        bits += md_rate_est_ctx->txfm_partition_fac_bits[ctx][1];

        if (allow_update_cdf)
            update_cdf(ec_ctx->txfm_partition_cdf[ctx], 1, 2);

        if (sub_txs == TX_4X4) {
            txfm_partition_update(xd->above_txfm_context + blk_col, xd->left_txfm_context + blk_row, sub_txs, tx_size);

            return bits;
        }

        assert(bsw > 0 && bsh > 0);
        for (int row = 0; row < tx_size_high_unit[tx_size]; row += bsh)
            for (int col = 0; col < tx_size_wide_unit[tx_size]; col += bsw) {
                int offsetr = blk_row + row;
                int offsetc = blk_col + col;
                bits += cost_tx_size_vartx(
                    xd, mbmi, sub_txs, depth + 1, offsetr, offsetc, md_rate_est_ctx, ec_ctx, allow_update_cdf);
            }
    }
    return bits;
}

static INLINE void set_txfm_ctx(TXFM_CONTEXT *txfm_ctx, uint8_t txs, int len) {
    int i;
    for (i = 0; i < len; ++i) txfm_ctx[i] = txs;
}

static INLINE void set_txfm_ctxs(TxSize tx_size, int n8_w, int n8_h, int skip, const MacroBlockD *xd) {
    uint8_t bw = tx_size_wide[tx_size];
    uint8_t bh = tx_size_high[tx_size];

    if (skip) {
        bw = n8_w * MI_SIZE;
        bh = n8_h * MI_SIZE;
    }

    set_txfm_ctx(xd->above_txfm_context, bw, n8_w);
    set_txfm_ctx(xd->left_txfm_context, bh, n8_h);
}

static INLINE int tx_size_to_depth(TxSize tx_size, BlockSize bsize) {
    TxSize ctx_size = max_txsize_rect_lookup[bsize];
    int    depth    = 0;
    while (tx_size != ctx_size) {
        depth++;
        ctx_size = sub_tx_size_map[ctx_size];
        assert(depth <= MAX_TX_DEPTH);
    }
    return depth;
}

// Returns a context number for the given MB prediction signal
// The mode info data structure has a one element border above and to the
// left of the entries corresponding to real blocks.
// The prediction flags in these dummy entries are initialized to 0.
static INLINE int get_tx_size_context(const MacroBlockD *xd) {
    const ModeInfo         *mi          = xd->mi[0];
    const MbModeInfo       *mbmi        = &mi->mbmi;
    const MbModeInfo *const above_mbmi  = xd->above_mbmi;
    const MbModeInfo *const left_mbmi   = xd->left_mbmi;
    const TxSize            max_tx_size = max_txsize_rect_lookup[mbmi->block_mi.bsize];
    const int               max_tx_wide = tx_size_wide[max_tx_size];
    const int               max_tx_high = tx_size_high[max_tx_size];
    const int               has_above   = xd->up_available;
    const int               has_left    = xd->left_available;

    int above = xd->above_txfm_context[0] >= max_tx_wide;
    int left  = xd->left_txfm_context[0] >= max_tx_high;

    if (has_above)
        if (is_inter_block(&above_mbmi->block_mi))
            above = block_size_wide[above_mbmi->block_mi.bsize] >= max_tx_wide;

    if (has_left)
        if (is_inter_block(&left_mbmi->block_mi))
            left = block_size_high[left_mbmi->block_mi.bsize] >= max_tx_high;

    if (has_above && has_left)
        return (above + left);
    else if (has_above)
        return above;
    else if (has_left)
        return left;
    else
        return 0;
}
static uint64_t cost_selected_tx_size(const MacroBlockD *xd, MdRateEstimationContext *md_rate_est_ctx, TxSize tx_size,
                                      FRAME_CONTEXT *ec_ctx, uint8_t allow_update_cdf) {
    const ModeInfo *const   mi    = xd->mi[0];
    const MbModeInfo *const mbmi  = &mi->mbmi;
    const BlockSize         bsize = mbmi->block_mi.bsize;
    uint64_t                bits  = 0;

    if (block_signals_txsize(bsize)) {
        const int tx_size_ctx = get_tx_size_context(xd);
        assert(bsize < BlockSizeS_ALL);
        const int     depth       = tx_size_to_depth(tx_size, bsize);
        const int32_t tx_size_cat = bsize_to_tx_size_cat(bsize);
        bits += md_rate_est_ctx->tx_size_fac_bits[tx_size_cat][tx_size_ctx][depth];

        if (allow_update_cdf) {
            const int max_depths = bsize_to_max_depth(bsize);
            assert(depth >= 0 && depth <= max_depths);
            assert(!is_inter_block(&mbmi->block_mi));
            assert(IMPLIES(is_rect_tx(tx_size), is_rect_tx_allowed(/*xd,*/ mbmi)));
            update_cdf(ec_ctx->tx_size_cdf[tx_size_cat][tx_size_ctx], depth, max_depths + 1);
        }
    }

    return bits;
}

/* Get the TXS rate and update the txfm context.  If allow_update_cdf is true, the TX size CDFs will
be updated. */
uint64_t svt_aom_tx_size_bits(MdRateEstimationContext *md_rate_est_ctx, MacroBlockD *xd, const MbModeInfo *mbmi,
                              TxSize tx_size, TxMode tx_mode, BlockSize bsize, uint8_t skip, FRAME_CONTEXT *ec_ctx,
                              uint8_t allow_update_cdf) {
    uint64_t bits        = 0;
    int      is_inter_tx = is_inter_block(&mbmi->block_mi);
    if (tx_mode == TX_MODE_SELECT && block_signals_txsize(bsize) &&
        !(is_inter_tx && skip) /*&& !xd->lossless[segment_id]*/) {
        if (is_inter_tx) { // This implies skip flag is 0.
            const TxSize max_tx_size = get_vartx_max_txsize(/*xd,*/ bsize, 0);
            const int    txbh        = tx_size_high_unit[max_tx_size];
            const int    txbw        = tx_size_wide_unit[max_tx_size];
            const int    width       = block_size_wide[bsize] >> tx_size_wide_log2[0];
            const int    height      = block_size_high[bsize] >> tx_size_high_log2[0];
            int          idx, idy;
            for (idy = 0; idy < height; idy += txbh)
                for (idx = 0; idx < width; idx += txbw)
                    bits += cost_tx_size_vartx(
                        xd, mbmi, max_tx_size, 0, idy, idx, md_rate_est_ctx, ec_ctx, allow_update_cdf);
        } else {
            bits += cost_selected_tx_size(xd, md_rate_est_ctx, tx_size, ec_ctx, allow_update_cdf);
            set_txfm_ctxs(tx_size, xd->n8_w, xd->n8_h, 0, xd);
        }
    } else {
        set_txfm_ctxs(tx_size, xd->n8_w, xd->n8_h, skip && is_inter_block(&mbmi->block_mi), xd);
    }

    return bits;
}

/* Get the TXS rate.  A dummy txfm context array will be used, so context updates will not be saved for
future blocks. */
uint64_t svt_aom_get_tx_size_bits(ModeDecisionCandidateBuffer *candidateBuffer, ModeDecisionContext *ctx,
                                  PictureControlSet *pcs, uint8_t tx_depth, Bool block_has_coeff) {
    NeighborArrayUnit *txfm_context_array      = ctx->txfm_context_array;
    uint32_t           txfm_context_left_index = get_neighbor_array_unit_left_index(txfm_context_array, ctx->blk_org_y);
    uint32_t           txfm_context_above_index = get_neighbor_array_unit_top_index(txfm_context_array, ctx->blk_org_x);

    TxMode       tx_mode = pcs->ppcs->frm_hdr.tx_mode;
    MacroBlockD *xd      = ctx->blk_ptr->av1xd;
    BlockSize    bsize   = ctx->blk_geom->bsize;
    MbModeInfo  *mbmi    = &xd->mi[0]->mbmi;

    svt_memcpy(ctx->above_txfm_context,
               &(txfm_context_array->top_array[txfm_context_above_index]),
               (ctx->blk_geom->bwidth >> MI_SIZE_LOG2) * sizeof(TXFM_CONTEXT));
    svt_memcpy(ctx->left_txfm_context,
               &(txfm_context_array->left_array[txfm_context_left_index]),
               (ctx->blk_geom->bheight >> MI_SIZE_LOG2) * sizeof(TXFM_CONTEXT));

    xd->above_txfm_context      = ctx->above_txfm_context;
    xd->left_txfm_context       = ctx->left_txfm_context;
    mbmi->block_mi.bsize        = ctx->blk_geom->bsize;
    mbmi->block_mi.use_intrabc  = candidateBuffer->cand->use_intrabc;
    mbmi->block_mi.ref_frame[0] = candidateBuffer->cand->ref_frame_type;
    mbmi->block_mi.tx_depth     = tx_depth;

    const uint64_t bits = svt_aom_tx_size_bits(
        ctx->md_rate_est_ctx, xd, mbmi, ctx->blk_geom->txsize[tx_depth], tx_mode, bsize, !block_has_coeff, NULL, 0);
    return bits;
}

/*
 * av1_partition_rate_cost function is used to generate the rate of signaling the
 * partition type for a given block.
 */
uint64_t svt_aom_partition_rate_cost(PictureParentControlSet *pcs, ModeDecisionContext *ctx, uint32_t blk_mds_idx,
                                     PartitionType p, uint64_t lambda, bool use_accurate_part_ctx,
                                     MdRateEstimationContext *md_rate_est_ctx) {
    const BlockGeom *blk_geom = get_blk_geom_mds(blk_mds_idx);
    const BlockSize  bsize    = blk_geom->bsize;
    assert(mi_size_wide_log2[bsize] == mi_size_high_log2[bsize]);
    assert(bsize < BlockSizeS_ALL);
    const Bool is_partition_point = (bsize >= BLOCK_8X8);

    if (!is_partition_point) {
        return 0;
    }

    const int blk_org_x = ctx->sb_origin_x + blk_geom->org_x;
    const int blk_org_y = ctx->sb_origin_y + blk_geom->org_y;
    const int hbs       = (mi_size_wide[bsize] << 2) >> 1;
    const int has_rows  = (blk_org_y + hbs) < pcs->aligned_height;
    const int has_cols  = (blk_org_x + hbs) < pcs->aligned_width;

    if (!has_rows && !has_cols) {
        return 0;
    }
    const PartitionContextType left_ctx  = ctx->md_local_blk_unit[blk_mds_idx].left_neighbor_partition ==
            (char)(INVALID_NEIGHBOR_DATA)
         ? 0
         : ctx->md_local_blk_unit[blk_mds_idx].left_neighbor_partition;
    const PartitionContextType above_ctx = ctx->md_local_blk_unit[blk_mds_idx].above_neighbor_partition ==
            (char)(INVALID_NEIGHBOR_DATA)
        ? 0
        : ctx->md_local_blk_unit[blk_mds_idx].above_neighbor_partition;
    const int                  bsl       = mi_size_wide_log2[bsize] - mi_size_wide_log2[BLOCK_8X8];
    assert(bsl >= 0);

    const int      above = (above_ctx >> bsl) & 1, left = (left_ctx >> bsl) & 1;
    const int      partitio_ploffset = use_accurate_part_ctx ? PARTITION_PLOFFSET : 0;
    const uint32_t context_index     = (left * 2 + above) + bsl * partitio_ploffset;

    uint64_t split_rate = 0;

    if (has_rows && has_cols) {
        split_rate = (uint64_t)md_rate_est_ctx->partition_fac_bits[context_index][p];
    } else if (!has_rows && has_cols) {
        split_rate = bsize == BLOCK_128X128
            ? (uint64_t)md_rate_est_ctx->partition_vert_alike_128x128_fac_bits[context_index][p == PARTITION_SPLIT]
            : (uint64_t)md_rate_est_ctx->partition_vert_alike_fac_bits[context_index][p == PARTITION_SPLIT];
    } else {
        split_rate = bsize == BLOCK_128X128
            ? (uint64_t)md_rate_est_ctx->partition_horz_alike_128x128_fac_bits[context_index][p == PARTITION_SPLIT]
            : (uint64_t)md_rate_est_ctx->partition_horz_alike_fac_bits[context_index][p == PARTITION_SPLIT];
    }

    return (RDCOST(lambda, split_rate, 0));
}
