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

/***************************************
* Includes
***************************************/
#include "EbRateDistortionCost.h"
#include "aom_dsp_rtcd.h"

#include <assert.h>

#define AV1_COST_PRECISION          0
#define MV_COST_WEIGHT              108
int av1_get_reference_mode_context_new(const MacroBlockD *xd);
int av1_get_pred_context_uni_comp_ref_p(const MacroBlockD *xd);
int av1_get_pred_context_uni_comp_ref_p1(const MacroBlockD *xd);
int av1_get_pred_context_uni_comp_ref_p2(const MacroBlockD *xd);
int av1_get_comp_reference_type_context_new(const MacroBlockD *xd);

BlockSize GetBlockSize(uint8_t cu_size) {
    return (cu_size == 64 ? BLOCK_64X64 : cu_size == 32 ? BLOCK_32X32 : cu_size == 16 ? BLOCK_16X16 : cu_size == 8 ? BLOCK_8X8 : BLOCK_4X4);
}

int av1_allow_intrabc(const Av1Common *const cm);
int32_t is_chroma_reference(int32_t mi_row, int32_t mi_col, BlockSize bsize,
    int32_t subsampling_x, int32_t subsampling_y) {
    const int32_t bw = mi_size_wide[bsize];
    const int32_t bh = mi_size_high[bsize];
    int32_t ref_pos = ((mi_row & 0x01) || !(bh & 0x01) || !subsampling_y) &&
        ((mi_col & 0x01) || !(bw & 0x01) || !subsampling_x);
    return ref_pos;
}

uint8_t av1_drl_ctx(const CandidateMv *ref_mv_stack,
    int32_t ref_idx) {
    if (ref_mv_stack[ref_idx].weight >= REF_CAT_LEVEL &&
        ref_mv_stack[ref_idx + 1].weight >= REF_CAT_LEVEL)
        return 0;

    if (ref_mv_stack[ref_idx].weight >= REF_CAT_LEVEL &&
        ref_mv_stack[ref_idx + 1].weight < REF_CAT_LEVEL)
        return 1;

    if (ref_mv_stack[ref_idx].weight < REF_CAT_LEVEL &&
        ref_mv_stack[ref_idx + 1].weight < REF_CAT_LEVEL)
        return 2;

    return 0;
}

/* Symbols for coding which components are zero jointly */
//#define MV_JOINTS 4
//typedef enum {
//    MV_JOINT_ZERO = 0,   /* Zero vector */
//    MV_JOINT_HNZVZ = 1,  /* Vert zero, hor nonzero */
//    MV_JOINT_HZVNZ = 2,  /* Hor zero, vert nonzero */
//    MV_JOINT_HNZVNZ = 3, /* Both components nonzero */
//} MvJointType;

MvJointType av1_get_mv_joint(const MV *mv) {
    if (mv->row == 0)
        return mv->col == 0 ? MV_JOINT_ZERO : MV_JOINT_HNZVZ;
    else
        return mv->col == 0 ? MV_JOINT_HZVNZ : MV_JOINT_HNZVNZ;
}
int32_t mv_cost(const MV *mv, const int32_t *joint_cost,
    int32_t *const comp_cost[2]) {
    int32_t jnC = av1_get_mv_joint(mv);
    int32_t res =
        joint_cost[jnC] + comp_cost[0][mv->row] +
        comp_cost[1][mv->col];

    return res;
}

int32_t av1_mv_bit_cost(const MV *mv, const MV *ref, const int32_t *mvjcost,
    int32_t *mvcost[2], int32_t weight) {
    const MV diff = { mv->row - ref->row, mv->col - ref->col };
    return ROUND_POWER_OF_TWO(mv_cost(&diff, mvjcost, mvcost) * weight, 7);
}

/////////////////////////////COEFFICIENT CALCULATION //////////////////////////////////////////////
static INLINE int32_t get_golomb_cost(int32_t abs_qc) {
    if (abs_qc >= 1 + NUM_BASE_LEVELS + COEFF_BASE_RANGE) {
        const int32_t r = abs_qc - COEFF_BASE_RANGE - NUM_BASE_LEVELS;
        const int32_t length = get_msb(r) + 1;
        return av1_cost_literal(2 * length - 1);
    }
    return 0;
}

static INLINE int32_t get_txb_bwl(TxSize tx_size) {
    tx_size = av1_get_adjusted_tx_size(tx_size);
    assert(tx_size < TX_SIZES_ALL);
    return tx_size_wide_log2[tx_size];
}

static INLINE int32_t get_txb_wide(TxSize tx_size) {
    tx_size = av1_get_adjusted_tx_size(tx_size);
    assert(tx_size < TX_SIZES_ALL);
    return tx_size_wide[tx_size];
}

static INLINE int32_t get_txb_high(TxSize tx_size) {
    tx_size = av1_get_adjusted_tx_size(tx_size);
    assert(tx_size < TX_SIZES_ALL);
    return tx_size_high[tx_size];
}
static INLINE uint8_t *set_levels(uint8_t *const levels_buf, const int32_t width) {
    return levels_buf + TX_PAD_TOP * (width + TX_PAD_HOR);
}

void av1_txb_init_levels_c(
    const TranLow *const coeff,
    const int32_t width,
    const int32_t height,
    uint8_t *const levels) {
    const int32_t stride = width + TX_PAD_HOR;
    uint8_t *ls = levels;

    memset(levels - TX_PAD_TOP * stride, 0,
        sizeof(*levels) * TX_PAD_TOP * stride);
    memset(levels + stride * height, 0,
        sizeof(*levels) * (TX_PAD_BOTTOM * stride + TX_PAD_END));

    for (int32_t i = 0; i < height; i++) {
        for (int32_t j = 0; j < width; j++)
            *ls++ = (uint8_t)clamp(abs(coeff[i * width + j]), 0, INT8_MAX);
        for (int32_t j = 0; j < TX_PAD_HOR; j++)
            *ls++ = 0;
    }
}

static const PredictionMode fimode_to_intradir[FILTER_INTRA_MODES] = {
    DC_PRED, V_PRED, H_PRED, D157_PRED, DC_PRED
};
// TODO(angiebird): use this function whenever it's possible
int32_t Av1TransformTypeRateEstimation(
    uint8_t        allow_update_cdf,
    FRAME_CONTEXT *fc,
    struct ModeDecisionCandidateBuffer    *candidate_buffer_ptr,
    EbBool                                  is_inter,
    EbBool                                  useFilterIntraFlag,
    TxSize                                  transform_size,
    TxType                                  transform_type,
    EbBool                                  reduced_tx_set_used)
{
    uint8_t filterIntraMode = 0; // AMIR to check// NM- hardcoded to zero for the moment until we support different intra filtering modes.

    //const MbModeInfo *mbmi = &xd->mi[0]->mbmi;
    //const int32_t is_inter = is_inter_block(mbmi);

    if (get_ext_tx_types(transform_size, is_inter, reduced_tx_set_used) > 1  /*&&    !xd->lossless[xd->mi[0]->mbmi.segment_id]  WE ARE NOT LOSSLESS*/) {
        const TxSize square_tx_size = txsize_sqr_map[transform_size];
        assert(square_tx_size < EXT_TX_SIZES);

        const int32_t ext_tx_set = get_ext_tx_set(transform_size, is_inter, reduced_tx_set_used);
        if (is_inter) {
            if (ext_tx_set > 0)
            {
                if (allow_update_cdf) {
                    const TxSetType tx_set_type =
                        get_ext_tx_set_type(transform_size, is_inter, reduced_tx_set_used);

                    update_cdf(fc->inter_ext_tx_cdf[ext_tx_set][square_tx_size],
                        av1_ext_tx_ind[tx_set_type][transform_type],
                        av1_num_ext_tx_set[tx_set_type]);
                }
                return candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->inter_tx_type_fac_bits[ext_tx_set][square_tx_size][transform_type];
            }
        }
        else {
            if (ext_tx_set > 0) {
                PredictionMode intra_dir;
                if (useFilterIntraFlag)
                    intra_dir = fimode_to_intradir[filterIntraMode];
                else
                    intra_dir = candidate_buffer_ptr->candidate_ptr->pred_mode;
                assert(intra_dir < INTRA_MODES);
                const TxSetType tx_set_type =
                    get_ext_tx_set_type(transform_size, is_inter, reduced_tx_set_used);

                if (allow_update_cdf) {
                    update_cdf(
                        fc->intra_ext_tx_cdf[ext_tx_set][square_tx_size][intra_dir],
                        av1_ext_tx_ind[tx_set_type][transform_type],
                        av1_num_ext_tx_set[tx_set_type]);
                }
                return candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->intra_tx_type_fac_bits[ext_tx_set][square_tx_size][intra_dir][transform_type];
            }
        }
    }
    return 0;
}

const int16_t k_eob_group_start[12] = { 0, 1, 2, 3, 5, 9, 17, 33, 65, 129, 257, 513 };
const int16_t k_eob_offset_bits[12] = { 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

static const int8_t eob_to_pos_small[33] = {
    0, 1, 2,                                        // 0-2
    3, 3,                                           // 3-4
    4, 4, 4, 4,                                     // 5-8
    5, 5, 5, 5, 5, 5, 5, 5,                         // 9-16
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6  // 17-32
};

static const int8_t eob_to_pos_large[17] = {
    6,                               // place holder
    7,                               // 33-64
    8, 8,                           // 65-128
    9, 9, 9, 9,                   // 129-256
    10, 10, 10, 10, 10, 10, 10, 10,  // 257-512
    11                               // 513-
};

static INLINE int32_t get_eob_pos_token(const int32_t eob, int32_t *const extra) {
    int32_t t;

    if (eob < 33)
        t = eob_to_pos_small[eob];
    else {
        const int32_t e = AOMMIN((eob - 1) >> 5, 16);
        t = eob_to_pos_large[e];
    }

    *extra = eob - k_eob_group_start[t];

    return t;
}
#define TX_SIZE TxSize
static INLINE TX_SIZE get_txsize_entropy_ctx(TX_SIZE txsize) {
    return (TX_SIZE)((txsize_sqr_map[txsize] + txsize_sqr_up_map[txsize] + 1) >>
        1);
}
void av1_update_eob_context(int eob, TX_SIZE tx_size, TxClass tx_class,
    PlaneType plane, FRAME_CONTEXT *ec_ctx,
    uint8_t allow_update_cdf) {
    int eob_extra;
    const int eob_pt = get_eob_pos_token(eob, &eob_extra);
    TX_SIZE txs_ctx = get_txsize_entropy_ctx(tx_size);
    assert(txs_ctx < TX_SIZES);
    const int eob_multi_size = txsize_log2_minus4[tx_size];
    const int eob_multi_ctx = (tx_class == TX_CLASS_2D) ? 0 : 1;

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
            update_cdf(ec_ctx->eob_flag_cdf128[plane][eob_multi_ctx], eob_pt - 1,
                8);
        }
        break;
    case 4:
#if CONFIG_ENTROPY_STATS
        ++counts->eob_multi256[cdf_idx][plane][eob_multi_ctx][eob_pt - 1];
#endif
        if (allow_update_cdf) {
            update_cdf(ec_ctx->eob_flag_cdf256[plane][eob_multi_ctx], eob_pt - 1,
                9);
        }
        break;
    case 5:
#if CONFIG_ENTROPY_STATS
        ++counts->eob_multi512[cdf_idx][plane][eob_multi_ctx][eob_pt - 1];
#endif
        if (allow_update_cdf) {
            update_cdf(ec_ctx->eob_flag_cdf512[plane][eob_multi_ctx], eob_pt - 1,
                10);
        }
        break;
    case 6:
    default:
#if CONFIG_ENTROPY_STATS
        ++counts->eob_multi1024[cdf_idx][plane][eob_multi_ctx][eob_pt - 1];
#endif
        if (allow_update_cdf) {
            update_cdf(ec_ctx->eob_flag_cdf1024[plane][eob_multi_ctx], eob_pt - 1,
                11);
        }
        break;
    }

    if (k_eob_offset_bits[eob_pt] > 0) {
        int eob_ctx = eob_pt - 3;
        int eob_shift = k_eob_offset_bits[eob_pt] - 1;
        int bit = (eob_extra & (1 << eob_shift)) ? 1 : 0;
#if CONFIG_ENTROPY_STATS
        counts->eob_extra[cdf_idx][txs_ctx][plane][eob_pt][bit]++;
#endif  // CONFIG_ENTROPY_STATS
        if (allow_update_cdf)
            update_cdf(ec_ctx->eob_extra_cdf[txs_ctx][plane][eob_ctx], bit, 2);
    }
}
static int32_t get_eob_cost(int32_t eob, const LvMapEobCost *txb_eob_costs,
    const LvMapCoeffCost *txb_costs, TxType tx_type) {
    int32_t eob_extra;
    const int32_t eob_pt = get_eob_pos_token(eob, &eob_extra);
    int32_t eob_cost = 0;
    const int32_t eob_multi_ctx = (tx_type_to_class[tx_type] == TX_CLASS_2D) ? 0 : 1;
    eob_cost = txb_eob_costs->eob_cost[eob_multi_ctx][eob_pt - 1];

    if (k_eob_offset_bits[eob_pt] > 0) {
        const int32_t eob_shift = k_eob_offset_bits[eob_pt] - 1;
        const int32_t bit = (eob_extra & (1 << eob_shift)) ? 1 : 0;
        eob_cost += txb_costs->eob_extra_cost[eob_pt][bit];
        const int32_t offset_bits = k_eob_offset_bits[eob_pt];
        if (offset_bits > 1) eob_cost += av1_cost_literal(offset_bits - 1);
    }
    return eob_cost;
}

static INLINE int32_t get_br_ctx(const uint8_t *const levels,
    const int32_t c,  // raster order
    const int32_t bwl, const TxType tx_type) {
    const int32_t row = c >> bwl;
    const int32_t col = c - (row << bwl);
    const int32_t stride = (1 << bwl) + TX_PAD_HOR;
    const TxClass tx_class = tx_type_to_class[tx_type];
    const int32_t pos = row * stride + col;
    int32_t mag = levels[pos + 1];
    mag += levels[pos + stride];
    switch (tx_class) {
    case TX_CLASS_2D:
        mag += levels[pos + stride + 1];
        mag = AOMMIN((mag + 1) >> 1, 6);
        if (c == 0) return mag;
        if ((row < 2) && (col < 2)) return mag + 7;
        break;
    case TX_CLASS_HORIZ:
        mag += levels[pos + 2];
        mag = AOMMIN((mag + 1) >> 1, 6);
        if (c == 0) return mag;
        if (col == 0) return mag + 7;
        break;
    case TX_CLASS_VERT:
        mag += levels[pos + (stride << 1)];
        mag = AOMMIN((mag + 1) >> 1, 6);
        if (c == 0) return mag;
        if (row == 0) return mag + 7;
        break;
    default: break;
    }

    return mag + 14;
}

static INLINE int32_t av1_cost_skip_txb(
    uint8_t        allow_update_cdf,
    FRAME_CONTEXT *ec_ctx,
    struct ModeDecisionCandidateBuffer    *candidate_buffer_ptr,
    TxSize                                  transform_size,
    PlaneType                               plane_type,
    int16_t                                   txb_skip_ctx)
{
    const TxSize txs_ctx = (TxSize)((txsize_sqr_map[transform_size] + txsize_sqr_up_map[transform_size] + 1) >> 1);
    assert(txs_ctx < TX_SIZES);
    const LvMapCoeffCost *const coeff_costs = &candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->coeff_fac_bits[txs_ctx][plane_type];
    if (allow_update_cdf)
        update_cdf(ec_ctx->txb_skip_cdf[txs_ctx][txb_skip_ctx], 1, 2);
    return coeff_costs->txb_skip_cost[txb_skip_ctx][1];
}
// Note: don't call this function when eob is 0.
uint64_t av1_cost_coeffs_txb(
    uint8_t                             allow_update_cdf,
    FRAME_CONTEXT                      *ec_ctx,
    struct ModeDecisionCandidateBuffer *candidate_buffer_ptr,
    const TranLow                      *const qcoeff,
    uint16_t                            eob,
    PlaneType                           plane_type,
    TxSize                              transform_size,
    TxType                              transform_type,
    int16_t                             txb_skip_ctx,
    int16_t                             dc_sign_ctx,
    EbBool                              reducedTransformSetFlag)

{
    //Note: there is a different version of this function in AOM that seems to be efficient as its name is:
    //warehouse_efficients_txb

    const TxSize txs_ctx = (TxSize)((txsize_sqr_map[transform_size] + txsize_sqr_up_map[transform_size] + 1) >> 1);
    const TxClass tx_class = tx_type_to_class[transform_type];
    int32_t c, cost;
    const int32_t bwl = get_txb_bwl(transform_size);
    const int32_t width = get_txb_wide(transform_size);
    const int32_t height = get_txb_high(transform_size);
    const ScanOrder *const scan_order = &av1_scan_orders[transform_size][transform_type]; // get_scan(tx_size, tx_type);
    const int16_t *const scan = scan_order->scan;
    uint8_t levels_buf[TX_PAD_2D];
    uint8_t *const levels = set_levels(levels_buf, width);
    DECLARE_ALIGNED(16, int8_t, coeff_contexts[MAX_TX_SQUARE]);
    assert(txs_ctx < TX_SIZES);
    const LvMapCoeffCost *const coeff_costs = &candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->coeff_fac_bits[txs_ctx][plane_type];

    const int32_t eob_multi_size = txsize_log2_minus4[transform_size];
    const LvMapEobCost *const eobBits = &candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->eob_frac_bits[eob_multi_size][plane_type];
    // eob must be greater than 0 here.
    assert(eob > 0);
    cost = coeff_costs->txb_skip_cost[txb_skip_ctx][0];

    if (allow_update_cdf)
        update_cdf(ec_ctx->txb_skip_cdf[txs_ctx][txb_skip_ctx], eob == 0, 2);
    av1_txb_init_levels(qcoeff, width, height, levels); // NM - Needs to be optimized - to be combined with the quantisation.

    // Transform type bit estimation
    cost += plane_type > PLANE_TYPE_Y ? 0 :
        Av1TransformTypeRateEstimation(
            allow_update_cdf,
            ec_ctx,
            candidate_buffer_ptr,
            candidate_buffer_ptr->candidate_ptr->type == INTER_MODE ? EB_TRUE : EB_FALSE,
            EB_FALSE, // NM - Hardcoded to false for the moment until we support the intra filtering
            transform_size,
            transform_type,
            reducedTransformSetFlag);

    // Transform ebo bit estimation
    int32_t eob_cost = get_eob_cost(eob, eobBits, coeff_costs, transform_type);
    cost += eob_cost;
    if (allow_update_cdf)
        av1_update_eob_context(eob, transform_size, tx_class,
            plane_type, ec_ctx, allow_update_cdf);
    // Transform non-zero coeff bit estimation
    av1_get_nz_map_contexts(
        levels,
        scan,
        eob,
        transform_size,
        tx_class,
        coeff_contexts); // NM - Assembly version is available in AOM

    if (allow_update_cdf)
    {
        for (int c = eob - 1; c >= 0; --c) {
            const int pos = scan[c];
            const int coeff_ctx = coeff_contexts[pos];
            const TranLow v = qcoeff[pos];
            const TranLow level = abs(v);

            if (allow_update_cdf) {
                if (c == eob - 1) {
                    assert(coeff_ctx < 4);
                    update_cdf(
                        ec_ctx->coeff_base_eob_cdf[txs_ctx][plane_type][coeff_ctx],
                        AOMMIN(level, 3) - 1, 3);
                }
                else {
                    update_cdf(ec_ctx->coeff_base_cdf[txs_ctx][plane_type][coeff_ctx],
                        AOMMIN(level, 3), 4);
                }
            }

            {
                if (c == eob - 1) {
                    assert(coeff_ctx < 4);
#if CONFIG_ENTROPY_STATS
                    ++td->counts->coeff_base_eob_multi[cdf_idx][txsize_ctx][plane_type]
                        [coeff_ctx][AOMMIN(level, 3) - 1];
                }
                else {
                    ++td->counts->coeff_base_multi[cdf_idx][txsize_ctx][plane_type]
                        [coeff_ctx][AOMMIN(level, 3)];
#endif
                }
            }

            if (level > NUM_BASE_LEVELS) {
                const int base_range = level - 1 - NUM_BASE_LEVELS;
                const int br_ctx = get_br_ctx(levels, pos, bwl, (const TxType)tx_class);

                for (int idx = 0; idx < COEFF_BASE_RANGE; idx += BR_CDF_SIZE - 1) {
                    const int k = AOMMIN(base_range - idx, BR_CDF_SIZE - 1);
                    if (allow_update_cdf) {
                        update_cdf(ec_ctx->coeff_br_cdf[AOMMIN(txs_ctx, TX_32X32)]
                            [plane_type][br_ctx],
                            k, BR_CDF_SIZE);
                    }
                    for (int lps = 0; lps < BR_CDF_SIZE - 1; lps++) {
#if CONFIG_ENTROPY_STATS
                        ++td->counts->coeff_lps[AOMMIN(txsize_ctx, TX_32X32)][plane_type][lps]
                            [br_ctx][lps == k];
#endif  // CONFIG_ENTROPY_STATS
                        if (lps == k) break;
                    }
#if CONFIG_ENTROPY_STATS
                    ++td->counts->coeff_lps_multi[cdf_idx][AOMMIN(txsize_ctx, TX_32X32)]
                        [plane_type][br_ctx][k];
#endif
                    if (k < BR_CDF_SIZE - 1) break;
                }
            }
        }

        if (qcoeff[0] != 0) {
            const int dc_sign = (qcoeff[0] < 0) ? 1 : 0;
            if (allow_update_cdf)
                update_cdf(ec_ctx->dc_sign_cdf[plane_type][dc_sign_ctx], dc_sign, 2);
        }

        //TODO: CHKN  for 128x128 where we need more than one TXb, we need to update the txb_context(dc_sign+skip_ctx) in a Txb basis.

        return 0;
    }

    for (c = eob - 1; c >= 0; --c) {
        const int32_t pos = scan[c];
        const TranLow v = qcoeff[pos];
        const int32_t is_nz = (v != 0);
        const int32_t level = abs(v);
        const int32_t coeff_ctx = coeff_contexts[pos];

        if (c == eob - 1) {
            assert((AOMMIN(level, 3) - 1) >= 0);
            cost += coeff_costs->base_eob_cost[coeff_ctx][AOMMIN(level, 3) - 1];
        }
        else
            cost += coeff_costs->base_cost[coeff_ctx][AOMMIN(level, 3)];
        if (is_nz) {
            const int32_t sign = (v < 0) ? 1 : 0;
            // sign bit cost
            if (c == 0)
                cost += coeff_costs->dc_sign_cost[dc_sign_ctx][sign];
            else
                cost += av1_cost_literal(1);
            if (level > NUM_BASE_LEVELS) {
                int32_t ctx;
                ctx = get_br_ctx(levels, pos, bwl, transform_type);

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
    return cost;
}
/*static*/ void model_rd_from_sse(
    BlockSize bsize,
    int16_t quantizer,
    //const Av1Comp *const cpi,
    //const MacroBlockD *const xd,
    //BlockSize bsize,
    //int32_t plane,
    uint64_t sse,
    uint32_t *rate,
    uint64_t *dist);

uint64_t av1_intra_fast_cost(
    CodingUnit            *cu_ptr,
    ModeDecisionCandidate *candidate_ptr,
    uint32_t                 qp,
    uint64_t                 luma_distortion,
    uint64_t                 chroma_distortion,
    uint64_t                 lambda,
    EbBool                   use_ssd,
    PictureControlSet     *picture_control_set_ptr,
    CandidateMv             *ref_mv_stack,
    const BlockGeom         *blk_geom,
    uint32_t                 miRow,
    uint32_t                 miCol,
    uint8_t                 md_pass,
    uint32_t                 left_neighbor_mode,
    uint32_t                 top_neighbor_mode)

{
    UNUSED(qp);
    UNUSED(ref_mv_stack);
    UNUSED(miRow);
    UNUSED(miCol);
    UNUSED(left_neighbor_mode);
    UNUSED(top_neighbor_mode);
    UNUSED(md_pass);

    FrameHeader *frm_hdr = &picture_control_set_ptr->parent_pcs_ptr->frm_hdr;
    if (av1_allow_intrabc(picture_control_set_ptr->parent_pcs_ptr->av1_cm) && candidate_ptr->use_intrabc) {
        uint64_t lumaSad = (LUMA_WEIGHT * luma_distortion) << AV1_COST_PRECISION;
        uint64_t chromaSad = chroma_distortion << AV1_COST_PRECISION;
        uint64_t totalDistortion = lumaSad + chromaSad;

        uint64_t rate = 0;

        EbReflist refListIdx = 0;
        int16_t predRefX = candidate_ptr->motion_vector_pred_x[refListIdx];
        int16_t predRefY = candidate_ptr->motion_vector_pred_y[refListIdx];
        int16_t mvRefX = candidate_ptr->motion_vector_xl0;
        int16_t mvRefY = candidate_ptr->motion_vector_yl0;
        MV mv;
        mv.row = mvRefY;
        mv.col = mvRefX;
        MV ref_mv;
        ref_mv.row = predRefY;
        ref_mv.col = predRefX;

        int *dvcost[2] = { (int *)&candidate_ptr->md_rate_estimation_ptr->dv_cost[0][MV_MAX],
                           (int *)&candidate_ptr->md_rate_estimation_ptr->dv_cost[1][MV_MAX] };

        int32_t mvRate = av1_mv_bit_cost(
            &mv,
            &ref_mv,
            candidate_ptr->md_rate_estimation_ptr->dv_joint_cost,
            dvcost, MV_COST_WEIGHT_SUB);

        rate = mvRate + candidate_ptr->md_rate_estimation_ptr->intrabc_fac_bits[candidate_ptr->use_intrabc];

        candidate_ptr->fast_luma_rate = rate;
        candidate_ptr->fast_chroma_rate = 0;

        lumaSad = (LUMA_WEIGHT * luma_distortion) << AV1_COST_PRECISION;
        chromaSad = chroma_distortion << AV1_COST_PRECISION;
        totalDistortion = lumaSad + chromaSad;

        return(RDCOST(lambda, rate, totalDistortion));
    }
    else {
    EbBool isMonochromeFlag = EB_FALSE; // NM - isMonochromeFlag is harcoded to false.
    EbBool isCflAllowed = (blk_geom->bwidth <= 32 && blk_geom->bheight <= 32) ? 1 : 0;

    uint8_t   subSamplingX = 1; // NM - subsampling_x is harcoded to 1 for 420 chroma sampling.
    uint8_t   subSamplingY = 1; // NM - subsampling_y is harcoded to 1 for 420 chroma sampling.
    // In fast loop CFL alphas are not know yet. The chroma mode bits are calculated based on DC Mode, and if CFL is the winner compared to CFL, ChromaBits are updated
    uint32_t chroma_mode = candidate_ptr->intra_chroma_mode == UV_CFL_PRED ? UV_DC_PRED : candidate_ptr->intra_chroma_mode;

    // Number of bits for each synatax element
    uint64_t intraModeBitsNum = 0;
    uint64_t intraLumaModeBitsNum = 0;
    uint64_t intraLumaAngModeBitsNum = 0;
    uint64_t intraChromaModeBitsNum = 0;
    uint64_t intraChromaAngModeBitsNum = 0;
    uint64_t skipModeRate = 0;
    uint8_t  skipModeCtx = cu_ptr->skip_flag_context; // NM - Harcoded to 1 until the skip_mode context is added.
    PredictionMode intra_mode = (PredictionMode)candidate_ptr->pred_mode;
    // Luma and chroma rate
    uint32_t rate;
    uint32_t lumaRate = 0;
    uint32_t chromaRate = 0;
    uint64_t lumaSad, chromaSad;

    // Luma and chroma distortion
    uint64_t totalDistortion;
    const int32_t AboveCtx = intra_mode_context[top_neighbor_mode];
    const int32_t LeftCtx = intra_mode_context[left_neighbor_mode];
    intraModeBitsNum = picture_control_set_ptr->slice_type != I_SLICE ? (uint64_t)candidate_ptr->md_rate_estimation_ptr->mb_mode_fac_bits[size_group_lookup[blk_geom->bsize]][intra_mode] : ZERO_COST;
    skipModeRate = picture_control_set_ptr->slice_type != I_SLICE ? (uint64_t)candidate_ptr->md_rate_estimation_ptr->skip_mode_fac_bits[skipModeCtx][0] : ZERO_COST;

    // Estimate luma nominal intra mode bits
    intraLumaModeBitsNum = picture_control_set_ptr->slice_type == I_SLICE ? (uint64_t)candidate_ptr->md_rate_estimation_ptr->y_mode_fac_bits[AboveCtx][LeftCtx][intra_mode] : ZERO_COST;
    // Estimate luma angular mode bits
    if (blk_geom->bsize >= BLOCK_8X8 && candidate_ptr->is_directional_mode_flag) {
        assert((intra_mode - V_PRED) < 8);
        assert((intra_mode - V_PRED) >= 0);
        intraLumaAngModeBitsNum = candidate_ptr->md_rate_estimation_ptr->angle_delta_fac_bits[intra_mode - V_PRED][MAX_ANGLE_DELTA + candidate_ptr->angle_delta[PLANE_TYPE_Y]];
    }

    // NM- Harcoded assuming luma mode is equal to chroma mode
    //if (!cm->seq_params.monochrome &&
    //    is_chroma_reference(mi_row, mi_col, bsize, xd->plane[1].subsampling_x,
    //    xd->plane[1].subsampling_y)) {
    //    mbmi->uv_mode =
    //        read_intra_mode_uv(ec_ctx, r, is_cfl_allowed(xd), mbmi->mode);
    //    if (mbmi->uv_mode == UV_CFL_PRED) {
    //        mbmi->cfl_alpha_idx =
    //            read_cfl_alphas(xd->tile_ctx, r, &mbmi->cfl_alpha_signs);
    //        xd->cfl.store_y = 1;
    //    }
    //    else {
    //        xd->cfl.store_y = 0;
    //    }
    //    mbmi->angle_delta[PLANE_TYPE_UV] =
    //        use_angle_delta && av1_is_directional_mode(get_uv_mode(mbmi->uv_mode))
    //        ? read_angle_delta(r,
    //        ec_ctx->angle_delta_cdf[mbmi->uv_mode - V_PRED])
    //        : 0;
    //}
    //else {
    //    // Avoid decoding angle_info if there is is no chroma prediction
    //    mbmi->uv_mode = UV_DC_PRED;
    //    xd->cfl.is_chroma_reference = 0;
    //    xd->cfl.store_y = 1;
    //}
    if (blk_geom->has_uv) {
        if (!isMonochromeFlag && is_chroma_reference(miRow, miCol, blk_geom->bsize, subSamplingX, subSamplingY)) {
            // Estimate luma nominal intra mode bits
            intraChromaModeBitsNum = (uint64_t)candidate_ptr->md_rate_estimation_ptr->intra_uv_mode_fac_bits[isCflAllowed][intra_mode][chroma_mode];
            // Estimate luma angular mode bits
            if (blk_geom->bsize >= BLOCK_8X8 && candidate_ptr->is_directional_chroma_mode_flag) {
                intraChromaAngModeBitsNum = candidate_ptr->md_rate_estimation_ptr->angle_delta_fac_bits[chroma_mode - V_PRED][MAX_ANGLE_DELTA + candidate_ptr->angle_delta[PLANE_TYPE_UV]];
            }
        }
    }

    uint32_t isInterRate = picture_control_set_ptr->slice_type != I_SLICE ? candidate_ptr->md_rate_estimation_ptr->intra_inter_fac_bits[cu_ptr->is_inter_ctx][0] : 0;
    lumaRate = (uint32_t)(intraModeBitsNum + skipModeRate + intraLumaModeBitsNum + intraLumaAngModeBitsNum + isInterRate);
    if (av1_allow_intrabc(picture_control_set_ptr->parent_pcs_ptr->av1_cm))
        lumaRate += candidate_ptr->md_rate_estimation_ptr->intrabc_fac_bits[candidate_ptr->use_intrabc];

    chromaRate = (uint32_t)(intraChromaModeBitsNum + intraChromaAngModeBitsNum);

    // Keep the Fast Luma and Chroma rate for future use
    candidate_ptr->fast_luma_rate = lumaRate;
    candidate_ptr->fast_chroma_rate = chromaRate;
    if (use_ssd) {
        int32_t current_q_index = frm_hdr->quantization_params.base_q_idx;
        Dequants *const dequants = &picture_control_set_ptr->parent_pcs_ptr->deq;

        int16_t quantizer = dequants->y_dequant_Q3[current_q_index][1];
        rate = 0;
        model_rd_from_sse(
            blk_geom->bsize,
            quantizer,
            luma_distortion,
            &rate,
            &lumaSad);
        lumaRate += rate;
        totalDistortion = lumaSad;

        rate = 0;
        model_rd_from_sse(
            blk_geom->bsize_uv,
            quantizer,
            chroma_distortion,
            &chromaRate,
            &chromaSad);
        chromaRate += rate;
        totalDistortion += chromaSad;

        rate = lumaRate + chromaRate;

        return(RDCOST(lambda, rate, totalDistortion));
    }
    else {
        lumaSad = (LUMA_WEIGHT * luma_distortion) << AV1_COST_PRECISION;
        chromaSad = chroma_distortion << AV1_COST_PRECISION;
        totalDistortion = lumaSad + chromaSad;

        rate = lumaRate + chromaRate;

        // Assign fast cost
        return(RDCOST(lambda, rate, totalDistortion));
    }
    }
}

//extern INLINE int32_t have_newmv_in_inter_mode(PredictionMode mode);
static INLINE int32_t have_newmv_in_inter_mode(PredictionMode mode) {
    return (mode == NEWMV || mode == NEW_NEWMV || mode == NEAREST_NEWMV ||
        mode == NEW_NEARESTMV || mode == NEAR_NEWMV || mode == NEW_NEARMV);
}

extern void av1_set_ref_frame(MvReferenceFrame *rf,
    int8_t ref_frame_type);

static INLINE int has_second_ref(const MbModeInfo *mbmi) {
    return mbmi->ref_frame[1] > INTRA_FRAME;
}

static INLINE int has_uni_comp_refs(const MbModeInfo *mbmi) {
    return has_second_ref(mbmi) && (!((mbmi->ref_frame[0] >= BWDREF_FRAME) ^
        (mbmi->ref_frame[1] >= BWDREF_FRAME)));
}

// This function encodes the reference frame
uint64_t EstimateRefFramesNumBits(
    PictureControlSet                    *picture_control_set_ptr,
    ModeDecisionCandidate                *candidate_ptr,
    CodingUnit                           *cu_ptr,
    uint32_t                                 bwidth,
    uint32_t                                 bheight,
    uint8_t                                  ref_frame_type,
    uint8_t                                   md_pass,
    EbBool                                is_compound)
{

    FrameHeader *frm_hdr = &picture_control_set_ptr->parent_pcs_ptr->frm_hdr;
    uint64_t refRateBits = 0;

    if (md_pass == 1) {
        uint64_t refRateA = 0;
        uint64_t refRateB = 0;
        uint64_t refRateC = 0;
        uint64_t refRateD = 0;
        uint64_t refRateE = 0;
        uint64_t refRateF = 0;
        uint64_t refRateG = 0;
        uint64_t refRateH = 0;
        uint64_t refRateI = 0;
        uint64_t refRateJ = 0;
        uint64_t refRateK = 0;
        uint64_t refRateL = 0;
        uint64_t refRateM = 0;
        uint64_t refRateN = 0;
        uint64_t refRateO = 0;
        uint64_t refRateP = 0;
        // const MbModeInfo *const mbmi = &cu_ptr->av1xd->mi[0]->mbmi;
        MbModeInfo *const mbmi = &cu_ptr->av1xd->mi[0]->mbmi;
        MvReferenceFrame refType[2];
        av1_set_ref_frame(refType, ref_frame_type);
        mbmi->ref_frame[0] = refType[0];
        mbmi->ref_frame[1] = refType[1];
        //const int is_compound = has_second_ref(mbmi);
        {
            // does the feature use compound prediction or not
            // (if not specified at the frame/segment level)
            if (frm_hdr->reference_mode == REFERENCE_MODE_SELECT) {
                if (MIN(bwidth, bheight) >= 8) {
                    //aom_write_symbol(w, is_compound, av1_get_reference_mode_cdf(cu_ptr->av1xd), 2);
                    int32_t context = av1_get_reference_mode_context_new(cu_ptr->av1xd);
                    refRateA = candidate_ptr->md_rate_estimation_ptr->comp_inter_fac_bits[context][is_compound];
                }
            }
            else {
                assert((!is_compound) ==
                    (frm_hdr->reference_mode == SINGLE_REFERENCE));
            }

            if (is_compound) {
                const CompReferenceType comp_ref_type = has_uni_comp_refs(mbmi)
                    ? UNIDIR_COMP_REFERENCE
                    : BIDIR_COMP_REFERENCE;

                const int pred_context = av1_get_comp_reference_type_context_new(cu_ptr->av1xd);
                refRateB = candidate_ptr->md_rate_estimation_ptr->comp_ref_type_fac_bits[pred_context][comp_ref_type];
                /*aom_write_symbol(w, comp_ref_type, av1_get_comp_reference_type_cdf(cu_ptr->av1xd),
                    2);*/

                if (comp_ref_type == UNIDIR_COMP_REFERENCE) {
                    //printf("ERROR[AN]: UNIDIR_COMP_REFERENCE not supported\n");
                    const int bit = mbmi->ref_frame[0] == BWDREF_FRAME;

                    const int pred_context = av1_get_pred_context_uni_comp_ref_p(cu_ptr->av1xd);
                    refRateC = candidate_ptr->md_rate_estimation_ptr->uni_comp_ref_fac_bits[pred_context][0][bit];
                    //cu_ptr->av1xd->tile_ctx->uni_comp_ref_cdf[pred_context][0];
                    //WRITE_REF_BIT(bit, uni_comp_ref_p);

                    if (!bit) {
                        assert(mbmi->ref_frame[0] == LAST_FRAME);
                        const int bit1 = mbmi->ref_frame[1] == LAST3_FRAME ||
                            mbmi->ref_frame[1] == GOLDEN_FRAME;
                        const int pred_context = av1_get_pred_context_uni_comp_ref_p1(cu_ptr->av1xd);
                        refRateD = candidate_ptr->md_rate_estimation_ptr->uni_comp_ref_fac_bits[pred_context][1][bit1];
                        //refRateD = cu_ptr->av1xd->tile_ctx->uni_comp_ref_cdf[pred_context][1];
                        //WRITE_REF_BIT(bit1, uni_comp_ref_p1);
                        if (bit1) {
                            const int bit2 = mbmi->ref_frame[1] == GOLDEN_FRAME;
                            const int pred_context = av1_get_pred_context_uni_comp_ref_p2(cu_ptr->av1xd);
                            refRateE = candidate_ptr->md_rate_estimation_ptr->uni_comp_ref_fac_bits[pred_context][2][bit2];

                            // refRateE = cu_ptr->av1xd->tile_ctx->uni_comp_ref_cdf[pred_context][2];
                             //WRITE_REF_BIT(bit2, uni_comp_ref_p2);
                        }
                    }
                    //else {
                    //    assert(mbmi->ref_frame[1] == ALTREF_FRAME);
                    //}
                    refRateBits = refRateA + refRateB + refRateC + refRateD + refRateE + refRateF + refRateG + refRateH + refRateI + refRateJ + refRateK + refRateL + refRateM;
                    return refRateBits;
                    //return;
                }

                assert(comp_ref_type == BIDIR_COMP_REFERENCE);

                const int bit = (mbmi->ref_frame[0] == GOLDEN_FRAME ||
                    mbmi->ref_frame[0] == LAST3_FRAME);
                const int pred_ctx = av1_get_pred_context_comp_ref_p(cu_ptr->av1xd);
                refRateF = candidate_ptr->md_rate_estimation_ptr->comp_ref_fac_bits[pred_ctx][0][bit];
                //refRateF = cu_ptr->av1xd->tile_ctx->comp_ref_cdf[pred_ctx][0];
                //WRITE_REF_BIT(bit, comp_ref_p);

                if (!bit) {
                    const int bit1 = mbmi->ref_frame[0] == LAST2_FRAME;
                    const int pred_context = av1_get_pred_context_comp_ref_p1(cu_ptr->av1xd);
                    refRateG = candidate_ptr->md_rate_estimation_ptr->comp_ref_fac_bits[pred_context][1][bit1];
                    //refRateG = cu_ptr->av1xd->tile_ctx->comp_ref_cdf[pred_context][1];
                    //WRITE_REF_BIT(bit1, comp_ref_p1);
                }
                else {
                    const int bit2 = mbmi->ref_frame[0] == GOLDEN_FRAME;
                    const int pred_context = av1_get_pred_context_comp_ref_p2(cu_ptr->av1xd);
                    refRateH = candidate_ptr->md_rate_estimation_ptr->comp_ref_fac_bits[pred_context][2][bit2];
                    //refRateH = cu_ptr->av1xd->tile_ctx->comp_ref_cdf[pred_context][2];
                    //WRITE_REF_BIT(bit2, comp_ref_p2);
                }

                const int bit_bwd = mbmi->ref_frame[1] == ALTREF_FRAME;
                const int pred_ctx_2 = av1_get_pred_context_comp_bwdref_p(cu_ptr->av1xd);
                refRateI = candidate_ptr->md_rate_estimation_ptr->comp_bwd_ref_fac_bits[pred_ctx_2][0][bit_bwd];
                //refRateI = cu_ptr->av1xd->tile_ctx->comp_bwdref_cdf[pred_ctx_2][0];
                //WRITE_REF_BIT(bit_bwd, comp_bwdref_p);

                if (!bit_bwd) {
                    const int pred_context = av1_get_pred_context_comp_bwdref_p1(cu_ptr->av1xd);
                    refRateJ = candidate_ptr->md_rate_estimation_ptr->comp_bwd_ref_fac_bits[pred_context][1][refType[1] == ALTREF2_FRAME];
                    //refRateJ = cu_ptr->av1xd->tile_ctx->comp_bwdref_cdf[pred_context][1];
                    //WRITE_REF_BIT(mbmi->ref_frame[1] == ALTREF2_FRAME, comp_bwdref_p1);
                }
            }
            else {
                const int bit0 = (mbmi->ref_frame[0] <= ALTREF_FRAME &&
                    mbmi->ref_frame[0] >= BWDREF_FRAME);
                refRateK = candidate_ptr->md_rate_estimation_ptr->single_ref_fac_bits[av1_get_pred_context_single_ref_p1(cu_ptr->av1xd)][0][bit0];
                //refRateK = cu_ptr->av1xd->tile_ctx->single_ref_cdf[av1_get_pred_context_single_ref_p1(cu_ptr->av1xd)][0];
                //WRITE_REF_BIT(bit0, single_ref_p1);

                if (bit0) {
                    const int bit1 = mbmi->ref_frame[0] == ALTREF_FRAME;
                    refRateL = candidate_ptr->md_rate_estimation_ptr->single_ref_fac_bits[av1_get_pred_context_single_ref_p2(cu_ptr->av1xd)][1][bit1];
                    //refRateL = cu_ptr->av1xd->tile_ctx->single_ref_cdf[av1_get_pred_context_single_ref_p2(cu_ptr->av1xd)][1];
                    //WRITE_REF_BIT(bit1, single_ref_p2);
                    if (!bit1) {
                        refRateM = candidate_ptr->md_rate_estimation_ptr->single_ref_fac_bits[av1_get_pred_context_single_ref_p6(cu_ptr->av1xd)][5][ref_frame_type == ALTREF2_FRAME];
                        //refRateM = cu_ptr->av1xd->tile_ctx->single_ref_cdf[av1_get_pred_context_single_ref_p6(cu_ptr->av1xd)][5];
                        //WRITE_REF_BIT(mbmi->ref_frame[0] == ALTREF2_FRAME, single_ref_p6);
                    }
                }
                else {
                    const int bit2 = (mbmi->ref_frame[0] == LAST3_FRAME ||
                        mbmi->ref_frame[0] == GOLDEN_FRAME);
                    refRateN = candidate_ptr->md_rate_estimation_ptr->single_ref_fac_bits[av1_get_pred_context_single_ref_p3(cu_ptr->av1xd)][2][bit2];
                    //refRateN = cu_ptr->av1xd->tile_ctx->single_ref_cdf[av1_get_pred_context_single_ref_p3(cu_ptr->av1xd)][2];
                    //WRITE_REF_BIT(bit2, single_ref_p3);
                    if (!bit2) {
                        const int bit3 = mbmi->ref_frame[0] != LAST_FRAME;
                        refRateO = candidate_ptr->md_rate_estimation_ptr->single_ref_fac_bits[av1_get_pred_context_single_ref_p4(cu_ptr->av1xd)][3][bit3];
                        //refRateO = cu_ptr->av1xd->tile_ctx->single_ref_cdf[av1_get_pred_context_single_ref_p4(cu_ptr->av1xd)][3];
                        //WRITE_REF_BIT(bit3, single_ref_p4);
                    }
                    else {
                        const int bit4 = mbmi->ref_frame[0] != LAST3_FRAME;
                        refRateP = candidate_ptr->md_rate_estimation_ptr->single_ref_fac_bits[av1_get_pred_context_single_ref_p5(cu_ptr->av1xd)][4][bit4];
                        //refRateP = cu_ptr->av1xd->tile_ctx->single_ref_cdf[av1_get_pred_context_single_ref_p5(cu_ptr->av1xd)][4];
                        //WRITE_REF_BIT(bit4, single_ref_p5);
                    }
                }
            }
        }
        refRateBits = refRateA + refRateB + refRateC + refRateD + refRateE + refRateF + refRateG + refRateH + refRateI +
            refRateJ + refRateK + refRateL + refRateM + refRateN + refRateO + refRateP;
    }
    else {
    uint64_t refRateA = 0;
    uint64_t refRateB = 0;
    uint64_t refRateC = 0;
    uint64_t refRateD = 0;
    uint64_t refRateE = 0;
    uint64_t refRateF = 0;
    uint64_t refRateG = 0;
    uint64_t refRateH = 0;
    uint64_t refRateI = 0;
    uint64_t refRateJ = 0;
    uint64_t refRateK = 0;
    uint64_t refRateL = 0;
    uint64_t refRateM = 0;

    // If segment level coding of this signal is disabled...
    // or the segment allows multiple reference frame options
    /*if (segfeature_active(&cm->seg, segment_id, SEG_LVL_REF_FRAME)) {
    assert(!is_compound);
    assert(mbmi->ref_frame[0] ==
    get_segdata(&cm->seg, segment_id, SEG_LVL_REF_FRAME));
    }
    else if (segfeature_active(&cm->seg, segment_id, SEG_LVL_SKIP) ||
    segfeature_active(&cm->seg, segment_id, SEG_LVL_GLOBALMV)) {
    assert(!is_compound);
    assert(mbmi->ref_frame[0] == LAST_FRAME);
    }
    else*/ {
    // does the feature use compound prediction or not
    // (if not specified at the frame/segment level)
        if (frm_hdr->reference_mode == REFERENCE_MODE_SELECT) {
            if (MIN(bwidth, bheight) >= 8) {
                int32_t context = 0;
                context = cu_ptr->reference_mode_context;
                assert(context >= 0 && context < 5);
                refRateA = candidate_ptr->md_rate_estimation_ptr->comp_inter_fac_bits[context][is_compound];
            }
        }
        else
            assert((!is_compound) == (frm_hdr->reference_mode == SINGLE_REFERENCE));
        int32_t context = 0;
        if (is_compound) {
            const CompReferenceType comp_ref_type = /*has_uni_comp_refs(mbmi)
                                                      ? UNIDIR_COMP_REFERENCE
                                                      : */BIDIR_COMP_REFERENCE;
            MvReferenceFrame refType[2];
            av1_set_ref_frame(refType, ref_frame_type);

            context = cu_ptr->compoud_reference_type_context;
            assert(context >= 0 && context < 5);
            refRateB = candidate_ptr->md_rate_estimation_ptr->comp_ref_type_fac_bits[context][comp_ref_type];

            if (comp_ref_type == UNIDIR_COMP_REFERENCE) {
                printf("ERROR[AN]: UNIDIR_COMP_REFERENCE not supported\n");
                //const int32_t bit = mbmi->ref_frame[0] == BWDREF_FRAME;
                //WRITE_REF_BIT(bit, uni_comp_ref_p);

                //if (!bit) {
                //    assert(mbmi->ref_frame[0] == LAST_FRAME);
                //    const int32_t bit1 = mbmi->ref_frame[1] == LAST3_FRAME ||
                //        mbmi->ref_frame[1] == GOLDEN_FRAME;
                //    WRITE_REF_BIT(bit1, uni_comp_ref_p1);
                //    if (bit1) {
                //        const int32_t bit2 = mbmi->ref_frame[1] == GOLDEN_FRAME;
                //        WRITE_REF_BIT(bit2, uni_comp_ref_p2);
                //    }
                //}
                //else {
                //    assert(mbmi->ref_frame[1] == ALTREF_FRAME);
                //}

                //return;
            }

            assert(comp_ref_type == BIDIR_COMP_REFERENCE);

            const int32_t bit = (refType[0] == GOLDEN_FRAME ||
                refType[0] == LAST3_FRAME);

            context = av1_get_pred_context_comp_ref_p(cu_ptr->av1xd);
            assert(context >= 0 && context < 3);
            refRateC = candidate_ptr->md_rate_estimation_ptr->comp_ref_fac_bits[context][0][bit];
            //            WRITE_REF_BIT(bit, comp_ref_p);

            if (!bit) {
                const int32_t bit1 = (refType[0] == LAST2_FRAME);
                context = av1_get_pred_context_comp_ref_p1(cu_ptr->av1xd);
                /*aom_write_symbol(ec_writer, bit1, frameContext->comp_ref_cdf[context][1],
                    2);*/
                assert(context >= 0 && context < 3);
                refRateD = candidate_ptr->md_rate_estimation_ptr->comp_ref_fac_bits[context][1][bit1];

                //WRITE_REF_BIT(bit1, comp_ref_p1);
            }
            else {
                const int32_t bit2 = (refType[0] == GOLDEN_FRAME);
                context = av1_get_pred_context_comp_ref_p2(cu_ptr->av1xd);
                /*aom_write_symbol(ec_writer, bit2, frameContext->comp_ref_cdf[context][2],
                    2);*/
                assert(context >= 0 && context < 3);
                refRateE = candidate_ptr->md_rate_estimation_ptr->comp_ref_fac_bits[context][2][bit2];

                //WRITE_REF_BIT(bit2, comp_ref_p2);
            }

            const int32_t bit_bwd = (refType[1] == ALTREF_FRAME);
            context = av1_get_pred_context_comp_bwdref_p(cu_ptr->av1xd);
            /*aom_write_symbol(ec_writer, bit_bwd, frameContext->comp_bwdref_cdf[context][0],
                2);*/
            assert(context >= 0 && context < 3);
            refRateF = candidate_ptr->md_rate_estimation_ptr->comp_bwd_ref_fac_bits[context][0][bit_bwd];
            //WRITE_REF_BIT(bit_bwd, comp_bwdref_p);

            if (!bit_bwd) {
                context = av1_get_pred_context_comp_bwdref_p1(cu_ptr->av1xd);
                /*aom_write_symbol(ec_writer, refType[1] == ALTREF2_FRAME, frameContext->comp_bwdref_cdf[context][1],
                    2);*/
                assert(context >= 0 && context < 3);
                refRateG = candidate_ptr->md_rate_estimation_ptr->comp_bwd_ref_fac_bits[context][1][refType[1] == ALTREF2_FRAME];
                //WRITE_REF_BIT(mbmi->ref_frame[1] == ALTREF2_FRAME, comp_bwdref_p1);
            }
        }
        else {
            const int32_t bit0 = (ref_frame_type <= ALTREF_FRAME &&
                ref_frame_type >= BWDREF_FRAME);//0

            context = av1_get_pred_context_single_ref_p1(cu_ptr->av1xd);
            /*aom_write_symbol(ec_writer, bit0, frameContext->single_ref_cdf[context][0],
                2);*/
            assert(context >= 0 && context < 3);
            refRateH = candidate_ptr->md_rate_estimation_ptr->single_ref_fac_bits[context][0][bit0];
            //WRITE_REF_BIT(bit0, single_ref_p1);

            if (bit0) {
                const int32_t bit1 = (ref_frame_type == ALTREF_FRAME);
                context = av1_get_pred_context_single_ref_p2(cu_ptr->av1xd);
                assert(context >= 0 && context < 3);
                /*aom_write_symbol(ec_writer, bit1, frameContext->single_ref_cdf[context][1],
                    2);*/
                refRateI = candidate_ptr->md_rate_estimation_ptr->single_ref_fac_bits[context][1][bit1];
                //WRITE_REF_BIT(bit1, single_ref_p2);

                if (!bit1) {
                    context = av1_get_pred_context_single_ref_p6(cu_ptr->av1xd);
                    /*aom_write_symbol(ec_writer, cu_ptr->prediction_unit_array[0].ref_frame_type == ALTREF2_FRAME, frameContext->single_ref_cdf[context][5],
                        2);*/
                    assert(context >= 0 && context < 3);
                    refRateJ = candidate_ptr->md_rate_estimation_ptr->single_ref_fac_bits[context][5][ref_frame_type == ALTREF2_FRAME];
                    //WRITE_REF_BIT(mbmi->ref_frame[0] == ALTREF2_FRAME, single_ref_p6);
                }
            }
            else {
                const int32_t bit2 = (ref_frame_type == LAST3_FRAME ||
                    ref_frame_type == GOLDEN_FRAME); //0
                context = av1_get_pred_context_single_ref_p3(cu_ptr->av1xd);
                /*aom_write_symbol(ec_writer, bit2, frameContext->single_ref_cdf[context][2],
                    2);*/
                assert(context >= 0 && context < 3);
                refRateK = candidate_ptr->md_rate_estimation_ptr->single_ref_fac_bits[context][2][bit2];
                //WRITE_REF_BIT(bit2, single_ref_p3);

                if (!bit2) {
                    const int32_t bit3 = (ref_frame_type != LAST_FRAME); //0;
                    context = av1_get_pred_context_single_ref_p4(cu_ptr->av1xd);
                    assert(context >= 0 && context < 3);
                    /*aom_write_symbol(ec_writer, bit3, frameContext->single_ref_cdf[context][3],
                        2);*/
                    refRateL = candidate_ptr->md_rate_estimation_ptr->single_ref_fac_bits[context][3][bit3];
                    //WRITE_REF_BIT(bit3, single_ref_p4);
                }
                else {
                    const int32_t bit4 = (ref_frame_type != LAST3_FRAME);
                    context = av1_get_pred_context_single_ref_p5(cu_ptr->av1xd);
                    /*aom_write_symbol(ec_writer, bit4, frameContext->single_ref_cdf[context][4],
                        2);*/
                    assert(context >= 0 && context < 3);
                    refRateM = candidate_ptr->md_rate_estimation_ptr->single_ref_fac_bits[context][4][bit4];
                    //WRITE_REF_BIT(bit4, single_ref_p5);
                }
            }
        }
    }

    refRateBits = refRateA + refRateB + refRateC + refRateD + refRateE + refRateF + refRateG + refRateH + refRateI + refRateJ + refRateK + refRateL + refRateM;

    }
    return refRateBits;
}
//extern INLINE int16_t Av1ModeContextAnalyzer(const int16_t *const mode_context, const MvReferenceFrame *const rf);

extern  int8_t av1_ref_frame_type(const MvReferenceFrame *const rf);
uint16_t compound_mode_ctx_map_2[3][COMP_NEWMV_CTXS] = {
   { 0, 1, 1, 1, 1 },
   { 1, 2, 3, 4, 4 },
   { 4, 4, 5, 6, 7 },
};
static INLINE int16_t Av1ModeContextAnalyzer(
    const int16_t *const mode_context, const MvReferenceFrame *const rf) {
    const int8_t ref_frame = av1_ref_frame_type(rf);

    if (rf[1] <= INTRA_FRAME) return mode_context[ref_frame];

    const int16_t newmv_ctx = mode_context[ref_frame] & NEWMV_CTX_MASK;
    const int16_t refmv_ctx =
        (mode_context[ref_frame] >> REFMV_OFFSET) & REFMV_CTX_MASK;
    assert((refmv_ctx >> 1) < 3);
    const int16_t comp_ctx = compound_mode_ctx_map_2[refmv_ctx >> 1][AOMMIN(
        newmv_ctx, COMP_NEWMV_CTXS - 1)];
    return comp_ctx;
}

uint64_t av1_inter_fast_cost(
    CodingUnit            *cu_ptr,
    ModeDecisionCandidate *candidate_ptr,
    uint32_t                 qp,
    uint64_t                 luma_distortion,
    uint64_t                 chroma_distortion,
    uint64_t                 lambda,
    EbBool                   use_ssd,
    PictureControlSet     *picture_control_set_ptr,
    CandidateMv             *ref_mv_stack,
    const BlockGeom         *blk_geom,
    uint32_t                 miRow,
    uint32_t                 miCol,
    uint8_t                 md_pass,
    uint32_t                 left_neighbor_mode,
    uint32_t                 top_neighbor_mode)

{
    UNUSED(top_neighbor_mode);
    UNUSED(left_neighbor_mode);
    UNUSED(miCol);
    UNUSED(miRow);

    FrameHeader *frm_hdr = &picture_control_set_ptr->parent_pcs_ptr->frm_hdr;

    // Luma rate
    uint32_t           lumaRate = 0;
    uint32_t           chromaRate = 0;
    uint64_t           mvRate = 0;
    uint64_t           skipModeRate;
    // Luma and chroma distortion
    uint64_t           lumaSad;
    uint64_t             chromaSad;
    uint64_t           totalDistortion;

    uint32_t           rate;

    int16_t           predRefX;
    int16_t           predRefY;
    int16_t           mvRefX;
    int16_t           mvRefY;

    EbReflist       refListIdx;

    (void)qp;

    PredictionMode inter_mode = (PredictionMode)candidate_ptr->pred_mode;

    uint64_t interModeBitsNum = 0;

    uint8_t skipModeCtx = cu_ptr->skip_flag_context;
    MvReferenceFrame rf[2];
    av1_set_ref_frame(rf, candidate_ptr->ref_frame_type);
    uint32_t modeCtx = Av1ModeContextAnalyzer(cu_ptr->inter_mode_ctx, rf);
    skipModeRate = candidate_ptr->md_rate_estimation_ptr->skip_mode_fac_bits[skipModeCtx][0];
    uint64_t referencePictureBitsNum = 0;

    //Reference Type and Mode Bit estimation

    referencePictureBitsNum = EstimateRefFramesNumBits(
        picture_control_set_ptr,
        candidate_ptr,
        cu_ptr,
        blk_geom->bwidth,
        blk_geom->bheight,
        candidate_ptr->ref_frame_type,
        md_pass,
        candidate_ptr->is_compound);

    if (candidate_ptr->is_compound)
        interModeBitsNum += candidate_ptr->md_rate_estimation_ptr->inter_compound_mode_fac_bits[modeCtx][INTER_COMPOUND_OFFSET(inter_mode)];
    else {
        //uint32_t newmv_ctx = modeCtx & NEWMV_CTX_MASK;
        //interModeBitsNum = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->new_mv_mode_fac_bits[mode_ctx][0];

        int16_t newmv_ctx = modeCtx & NEWMV_CTX_MASK;
        //aom_write_symbol(ec_writer, mode != NEWMV, frameContext->newmv_cdf[newmv_ctx], 2);
        interModeBitsNum += candidate_ptr->md_rate_estimation_ptr->new_mv_mode_fac_bits[newmv_ctx][inter_mode != NEWMV];
        if (inter_mode != NEWMV) {
            const int16_t zeromvCtx = (modeCtx >> GLOBALMV_OFFSET) & GLOBALMV_CTX_MASK;
            //aom_write_symbol(ec_writer, mode != GLOBALMV, frameContext->zeromv_cdf[zeromvCtx], 2);
            interModeBitsNum += candidate_ptr->md_rate_estimation_ptr->zero_mv_mode_fac_bits[zeromvCtx][inter_mode != GLOBALMV];
            if (inter_mode != GLOBALMV) {
                int16_t refmvCtx = (modeCtx >> REFMV_OFFSET) & REFMV_CTX_MASK;
                /*aom_write_symbol(ec_writer, mode != NEARESTMV, frameContext->refmv_cdf[refmv_ctx], 2);*/
                interModeBitsNum += candidate_ptr->md_rate_estimation_ptr->ref_mv_mode_fac_bits[refmvCtx][inter_mode != NEARESTMV];
            }
        }
    }
    if (inter_mode == NEWMV || inter_mode == NEW_NEWMV || have_nearmv_in_inter_mode(inter_mode)) {
        //drLIdex cost estimation
        const int32_t new_mv = inter_mode == NEWMV || inter_mode == NEW_NEWMV;
        if (new_mv) {
            int32_t idx;
            for (idx = 0; idx < 2; ++idx) {
                if (cu_ptr->av1xd->ref_mv_count[candidate_ptr->ref_frame_type] > idx + 1) {
                    uint8_t drl1Ctx =
                        av1_drl_ctx(ref_mv_stack, idx);
                    interModeBitsNum += candidate_ptr->md_rate_estimation_ptr->drl_mode_fac_bits[drl1Ctx][candidate_ptr->drl_index != idx];
                    if (candidate_ptr->drl_index == idx) break;
                }
            }
        }

        if (have_nearmv_in_inter_mode(inter_mode)) {
            int32_t idx;
            // TODO(jingning): Temporary solution to compensate the NEARESTMV offset.
            for (idx = 1; idx < 3; ++idx) {
                if (cu_ptr->av1xd->ref_mv_count[candidate_ptr->ref_frame_type] > idx + 1) {
                    uint8_t drl_ctx =
                        av1_drl_ctx(ref_mv_stack, idx);
                    interModeBitsNum += candidate_ptr->md_rate_estimation_ptr->drl_mode_fac_bits[drl_ctx][candidate_ptr->drl_index != (idx - 1)];

                    if (candidate_ptr->drl_index == (idx - 1)) break;
                }
            }
        }
    }

    if (have_newmv_in_inter_mode(inter_mode)) {
        if (candidate_ptr->is_compound) {
            mvRate = 0;

            if (inter_mode == NEW_NEWMV) {
                for (refListIdx = 0; refListIdx < 2; ++refListIdx) {
                    predRefX = candidate_ptr->motion_vector_pred_x[refListIdx];
                    predRefY = candidate_ptr->motion_vector_pred_y[refListIdx];
                    mvRefX = refListIdx == REF_LIST_1 ? candidate_ptr->motion_vector_xl1 : candidate_ptr->motion_vector_xl0;
                    mvRefY = refListIdx == REF_LIST_1 ? candidate_ptr->motion_vector_yl1 : candidate_ptr->motion_vector_yl0;

                    MV mv;
                    mv.row = mvRefY;
                    mv.col = mvRefX;

                    MV ref_mv;
                    ref_mv.row = predRefY;
                    ref_mv.col = predRefX;

                    mvRate += av1_mv_bit_cost(
                        &mv,
                        &ref_mv,
                        candidate_ptr->md_rate_estimation_ptr->nmv_vec_cost,
                        candidate_ptr->md_rate_estimation_ptr->nmvcoststack,
                        MV_COST_WEIGHT);
                }
            }
            else if (inter_mode == NEAREST_NEWMV || inter_mode == NEAR_NEWMV) {
                predRefX = candidate_ptr->motion_vector_pred_x[REF_LIST_1];
                predRefY = candidate_ptr->motion_vector_pred_y[REF_LIST_1];
                mvRefX = candidate_ptr->motion_vector_xl1;
                mvRefY = candidate_ptr->motion_vector_yl1;

                MV mv;
                mv.row = mvRefY;
                mv.col = mvRefX;

                MV ref_mv;
                ref_mv.row = predRefY;
                ref_mv.col = predRefX;

                mvRate += av1_mv_bit_cost(
                    &mv,
                    &ref_mv,
                    candidate_ptr->md_rate_estimation_ptr->nmv_vec_cost,
                    candidate_ptr->md_rate_estimation_ptr->nmvcoststack,
                    MV_COST_WEIGHT);
            }
            else {
                assert(inter_mode == NEW_NEARESTMV || inter_mode == NEW_NEARMV);

                predRefX = candidate_ptr->motion_vector_pred_x[REF_LIST_0];
                predRefY = candidate_ptr->motion_vector_pred_y[REF_LIST_0];
                mvRefX = candidate_ptr->motion_vector_xl0;
                mvRefY = candidate_ptr->motion_vector_yl0;

                MV mv;
                mv.row = mvRefY;
                mv.col = mvRefX;

                MV ref_mv;
                ref_mv.row = predRefY;
                ref_mv.col = predRefX;

                mvRate += av1_mv_bit_cost(
                    &mv,
                    &ref_mv,
                    candidate_ptr->md_rate_estimation_ptr->nmv_vec_cost,
                    candidate_ptr->md_rate_estimation_ptr->nmvcoststack,
                    MV_COST_WEIGHT);
            }
        }
        else {
            refListIdx = candidate_ptr->prediction_direction[0] == 0 ? 0 : 1;

            predRefX = candidate_ptr->motion_vector_pred_x[refListIdx];
            predRefY = candidate_ptr->motion_vector_pred_y[refListIdx];

            mvRefX = refListIdx == 0 ? candidate_ptr->motion_vector_xl0 : candidate_ptr->motion_vector_xl1;
            mvRefY = refListIdx == 0 ? candidate_ptr->motion_vector_yl0 : candidate_ptr->motion_vector_yl1;

            MV mv;
            mv.row = mvRefY;
            mv.col = mvRefX;

            MV ref_mv;
            ref_mv.row = predRefY;
            ref_mv.col = predRefX;

            mvRate = av1_mv_bit_cost(
                &mv,
                &ref_mv,
                candidate_ptr->md_rate_estimation_ptr->nmv_vec_cost,
                candidate_ptr->md_rate_estimation_ptr->nmvcoststack,
                MV_COST_WEIGHT);
        }
    }

    // NM - To be added when the intrainter mode is adopted
    //  read_interintra_mode(is_compound)

    EbBool is_inter = inter_mode >= SINGLE_INTER_MODE_START && inter_mode < SINGLE_INTER_MODE_END;
    if (is_inter
        && frm_hdr->is_motion_mode_switchable
        && rf[1] != INTRA_FRAME)
    {
        MotionMode motion_mode_rd = candidate_ptr->motion_mode;
        BlockSize bsize = blk_geom->bsize;
        cu_ptr->prediction_unit_array[0].num_proj_ref = candidate_ptr->num_proj_ref;
        MotionMode last_motion_mode_allowed = motion_mode_allowed(
            picture_control_set_ptr,
            cu_ptr,
            bsize,
            rf[0],
            rf[1],
            inter_mode);

        switch (last_motion_mode_allowed) {
        case SIMPLE_TRANSLATION: break;
        case OBMC_CAUSAL:
            assert(motion_mode_rd == SIMPLE_TRANSLATION); // TODO: remove when OBMC added
            interModeBitsNum += candidate_ptr->md_rate_estimation_ptr->motion_mode_fac_bits1[bsize][motion_mode_rd];
            break;
        default:
            interModeBitsNum += candidate_ptr->md_rate_estimation_ptr->motion_mode_fac_bits[bsize][motion_mode_rd];
        }
    }

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
    uint32_t isInterRate = candidate_ptr->md_rate_estimation_ptr->intra_inter_fac_bits[cu_ptr->is_inter_ctx][1];
    lumaRate = (uint32_t)(referencePictureBitsNum + skipModeRate + interModeBitsNum + mvRate + isInterRate);

    //chromaRate = intraChromaModeBitsNum + intraChromaAngModeBitsNum;

    // Keep the Fast Luma and Chroma rate for future use
    candidate_ptr->fast_luma_rate = lumaRate;
    candidate_ptr->fast_chroma_rate = chromaRate;

    if (use_ssd) {
        int32_t current_q_index = frm_hdr->quantization_params.base_q_idx;
        Dequants *const dequants = &picture_control_set_ptr->parent_pcs_ptr->deq;

        int16_t quantizer = dequants->y_dequant_Q3[current_q_index][1];
        rate = 0;
        model_rd_from_sse(
            blk_geom->bsize,
            quantizer,
            luma_distortion,
            &rate,
            &lumaSad);
        lumaRate += rate;
        totalDistortion = lumaSad;

        rate = 0;
        model_rd_from_sse(
            blk_geom->bsize_uv,
            quantizer,
            chroma_distortion,
            &chromaRate,
            &chromaSad);
        chromaRate += rate;
        totalDistortion += chromaSad;

        rate = lumaRate + chromaRate;

        if (candidate_ptr->merge_flag) {
            uint64_t skipModeRate = candidate_ptr->md_rate_estimation_ptr->skip_mode_fac_bits[skipModeCtx][1];
            if (skipModeRate < rate)
                return(RDCOST(lambda, skipModeRate, totalDistortion));
        }
        return(RDCOST(lambda, rate, totalDistortion));
    }
    else {
        lumaSad = (LUMA_WEIGHT * luma_distortion) << AV1_COST_PRECISION;
        chromaSad = chroma_distortion << AV1_COST_PRECISION;
        totalDistortion = lumaSad + chromaSad;
        if (blk_geom->has_uv == 0 && chromaSad != 0)
            printf("av1_inter_fast_cost: Chroma error");
        rate = lumaRate + chromaRate;

        // Assign fast cost
        if (candidate_ptr->merge_flag) {
            uint64_t skipModeRate = candidate_ptr->md_rate_estimation_ptr->skip_mode_fac_bits[skipModeCtx][1];
            if (skipModeRate < rate)
                return(RDCOST(lambda, skipModeRate, totalDistortion));
        }
        return(RDCOST(lambda, rate, totalDistortion));
    }
}

EbErrorType av1_tu_estimate_coeff_bits(
    struct ModeDecisionContext         *md_context,
    uint8_t                             allow_update_cdf,
    FRAME_CONTEXT                      *ec_ctx,
    PictureControlSet                  *picture_control_set_ptr,
    struct ModeDecisionCandidateBuffer *candidate_buffer_ptr,
    uint32_t                            tu_origin_index,
    uint32_t                            tu_chroma_origin_index,
    EntropyCoder                       *entropy_coder_ptr,
    EbPictureBufferDesc                *coeff_buffer_sb,
    uint32_t                            y_eob,
    uint32_t                            cb_eob,
    uint32_t                            cr_eob,
    uint64_t                           *y_tu_coeff_bits,
    uint64_t                           *cb_tu_coeff_bits,
    uint64_t                           *cr_tu_coeff_bits,
    TxSize                              txsize,
    TxSize                              txsize_uv,
    TxType                              tx_type,
    TxType                              tx_type_uv,
    COMPONENT_TYPE                      component_type,
    EbAsm                               asm_type)
{
    (void)asm_type;
    (void)entropy_coder_ptr;
    EbErrorType return_error = EB_ErrorNone;

    FrameHeader *frm_hdr = &picture_control_set_ptr->parent_pcs_ptr->frm_hdr;

    int32_t *coeff_buffer;
    int16_t  luma_txb_skip_context = md_context->luma_txb_skip_context;
    int16_t  luma_dc_sign_context = md_context->luma_dc_sign_context;
    int16_t  cb_txb_skip_context = md_context->cb_txb_skip_context;
    int16_t  cb_dc_sign_context = md_context->cb_dc_sign_context;
    int16_t  cr_txb_skip_context = md_context->cr_txb_skip_context;
    int16_t  cr_dc_sign_context = md_context->cr_dc_sign_context;

    EbBool reducedTransformSetFlag = frm_hdr->reduced_tx_set ? EB_TRUE : EB_FALSE;

    //Estimate the rate of the transform type and coefficient for Luma

    if (component_type == COMPONENT_LUMA || component_type == COMPONENT_ALL) {
        if (y_eob) {
            coeff_buffer = (int32_t*)&coeff_buffer_sb->buffer_y[tu_origin_index * sizeof(int32_t)];

            *y_tu_coeff_bits = av1_cost_coeffs_txb(
                allow_update_cdf,
                ec_ctx,
                candidate_buffer_ptr,
                coeff_buffer,
                (uint16_t)y_eob,
                PLANE_TYPE_Y,
                txsize,
                tx_type,
                luma_txb_skip_context,
                luma_dc_sign_context,
                reducedTransformSetFlag);
        }
        else {
            *y_tu_coeff_bits = av1_cost_skip_txb(
                allow_update_cdf,
                ec_ctx,
                candidate_buffer_ptr,
                txsize,
                PLANE_TYPE_Y,
                luma_txb_skip_context);
        }
    }
    //Estimate the rate of the transform type and coefficient for chroma Cb

    if (component_type == COMPONENT_CHROMA_CB || component_type == COMPONENT_CHROMA || component_type == COMPONENT_ALL) {
        if (cb_eob) {
            coeff_buffer = (int32_t*)&coeff_buffer_sb->buffer_cb[tu_chroma_origin_index * sizeof(int32_t)];

            *cb_tu_coeff_bits = av1_cost_coeffs_txb(
                allow_update_cdf,
                ec_ctx,
                candidate_buffer_ptr,
                coeff_buffer,
                (uint16_t)cb_eob,
                PLANE_TYPE_UV,
                txsize_uv,
                tx_type_uv,
                cb_txb_skip_context,
                cb_dc_sign_context,
                reducedTransformSetFlag);
        }
        else {
            *cb_tu_coeff_bits = av1_cost_skip_txb(
                allow_update_cdf,
                ec_ctx,
                candidate_buffer_ptr,
                txsize_uv,
                PLANE_TYPE_UV,
                cb_txb_skip_context);
        }
    }

    if (component_type == COMPONENT_CHROMA_CR || component_type == COMPONENT_CHROMA || component_type == COMPONENT_ALL) {
        //Estimate the rate of the transform type and coefficient for chroma Cr
        if (cr_eob) {
            coeff_buffer = (int32_t*)&coeff_buffer_sb->buffer_cr[tu_chroma_origin_index * sizeof(int32_t)];

            *cr_tu_coeff_bits = av1_cost_coeffs_txb(
                allow_update_cdf,
                ec_ctx,
                candidate_buffer_ptr,
                coeff_buffer,
                (uint16_t)cr_eob,
                PLANE_TYPE_UV,
                txsize_uv,
                tx_type_uv,
                cr_txb_skip_context,
                cr_dc_sign_context,
                reducedTransformSetFlag);
        }
        else {
            *cr_tu_coeff_bits = av1_cost_skip_txb(
                allow_update_cdf,
                ec_ctx,
                candidate_buffer_ptr,
                txsize_uv,
                PLANE_TYPE_UV,
                cr_txb_skip_context);
        }
    }

    return return_error;
}

uint64_t estimate_tx_size_bits(
    PictureControlSet       *pcsPtr,
    uint32_t                 cu_origin_x,
    uint32_t                 cu_origin_y,
    CodingUnit               *cu_ptr,
    const BlockGeom          *blk_geom,
    NeighborArrayUnit        *txfm_context_array,
    uint8_t                   tx_depth,
    MdRateEstimationContext  *md_rate_estimation_ptr);

/*********************************************************************************
* av1_intra_full_cost function is used to estimate the cost of an intra candidate mode
* for full mode decisoion module.
*
*   @param *cu_ptr(input)
*       cu_ptr is the pointer of the target CU.
*   @param *candidate_buffer_ptr(input)
*       chromaBufferPtr is the buffer pointer of the candidate luma mode.
*   @param qp(input)
*       qp is the quantizer parameter.
*   @param luma_distortion (input)
*       luma_distortion is the intra condidate luma distortion.
*   @param lambda(input)
*       lambda is the Lagrange multiplier
**********************************************************************************/
EbErrorType Av1FullCost(
    PictureControlSet                    *picture_control_set_ptr,
    ModeDecisionContext                  *context_ptr,
    struct ModeDecisionCandidateBuffer   *candidate_buffer_ptr,
    CodingUnit                           *cu_ptr,
    uint64_t                               *y_distortion,
    uint64_t                               *cb_distortion,
    uint64_t                               *cr_distortion,
    uint64_t                                lambda,
    uint64_t                               *y_coeff_bits,
    uint64_t                               *cb_coeff_bits,
    uint64_t                               *cr_coeff_bits,
    BlockSize                               bsize)
{
    UNUSED(picture_control_set_ptr);
    UNUSED(bsize);
    UNUSED(cu_ptr);
    EbErrorType return_error = EB_ErrorNone;

    // Luma and chroma rate
    uint64_t lumaRate = 0;
    uint64_t chromaRate = 0;
    uint64_t coeffRate = 0;

    // Luma and chroma SSE
    uint64_t luma_sse;
    uint64_t chromaSse;
    uint64_t totalDistortion;
    uint64_t rate;

    //Estimate the rate of the transform type and coefficient for Luma
    // Add fast rate to get the total rate of the subject mode
    lumaRate += candidate_buffer_ptr->candidate_ptr->fast_luma_rate;
    chromaRate += candidate_buffer_ptr->candidate_ptr->fast_chroma_rate;

    // For CFL, costs of alphas are not computed in fast loop, since they are computed in the full loop. The rate costs are added to the full loop.
    // In fast loop CFL alphas are not know yet. The chroma mode bits are calculated based on DC Mode, and if CFL is the winner compared to CFL, ChromaBits are updated in Full loop
    if (context_ptr->blk_geom->has_uv) {
        if (candidate_buffer_ptr->candidate_ptr->type == INTRA_MODE && candidate_buffer_ptr->candidate_ptr->intra_chroma_mode == UV_CFL_PRED) {
            EbBool isCflAllowed = (context_ptr->blk_geom->bwidth <= 32 &&
                context_ptr->blk_geom->bheight <= 32) ? 1 : 0;

            chromaRate += candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->cfl_alpha_fac_bits[candidate_buffer_ptr->candidate_ptr->cfl_alpha_signs][CFL_PRED_U][CFL_IDX_U(candidate_buffer_ptr->candidate_ptr->cfl_alpha_idx)] +
                candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->cfl_alpha_fac_bits[candidate_buffer_ptr->candidate_ptr->cfl_alpha_signs][CFL_PRED_V][CFL_IDX_V(candidate_buffer_ptr->candidate_ptr->cfl_alpha_idx)];

            chromaRate += (uint64_t)candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->intra_uv_mode_fac_bits[isCflAllowed][candidate_buffer_ptr->candidate_ptr->intra_luma_mode][UV_CFL_PRED];
            chromaRate -= (uint64_t)candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->intra_uv_mode_fac_bits[isCflAllowed][candidate_buffer_ptr->candidate_ptr->intra_luma_mode][UV_DC_PRED];
        }
    }

    // Coeff rate

    if (context_ptr->blk_skip_decision && candidate_buffer_ptr->candidate_ptr->type != INTRA_MODE) {
        uint64_t non_skip_cost = RDCOST(lambda, (*y_coeff_bits + *cb_coeff_bits + *cr_coeff_bits + (uint64_t)candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->skip_fac_bits[cu_ptr->skip_coeff_context][0]), (y_distortion[0] + cb_distortion[0] + cr_distortion[0]));
        uint64_t skip_cost = RDCOST(lambda, ((uint64_t)candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->skip_fac_bits[cu_ptr->skip_coeff_context][1]), (y_distortion[1] + cb_distortion[1] + cr_distortion[1]));
        if ((candidate_buffer_ptr->candidate_ptr->block_has_coeff == 0) || (skip_cost < non_skip_cost)) {
            y_distortion[0] = y_distortion[1];
            cb_distortion[0] = cb_distortion[1];
            cr_distortion[0] = cr_distortion[1];
            candidate_buffer_ptr->candidate_ptr->block_has_coeff = 0;
        }
        if (candidate_buffer_ptr->candidate_ptr->block_has_coeff)
            coeffRate = (*y_coeff_bits + *cb_coeff_bits + *cr_coeff_bits + (uint64_t)candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->skip_fac_bits[cu_ptr->skip_coeff_context][0]);
        else
            coeffRate = MIN((uint64_t)candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->skip_fac_bits[cu_ptr->skip_coeff_context][1],
            (*y_coeff_bits + *cb_coeff_bits + *cr_coeff_bits + (uint64_t)candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->skip_fac_bits[cu_ptr->skip_coeff_context][0]));
    }
    else
        coeffRate = (*y_coeff_bits + *cb_coeff_bits + *cr_coeff_bits + (uint64_t)candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->skip_fac_bits[cu_ptr->skip_coeff_context][0]);
    luma_sse = y_distortion[0];
    chromaSse = cb_distortion[0] + cr_distortion[0];
    totalDistortion = luma_sse + chromaSse;

    rate = lumaRate + chromaRate + coeffRate;

    if (candidate_buffer_ptr->candidate_ptr->block_has_coeff) {
        uint64_t tx_size_bits = estimate_tx_size_bits(
            picture_control_set_ptr,
            context_ptr->cu_origin_x,
            context_ptr->cu_origin_y,
            context_ptr->cu_ptr,
            context_ptr->blk_geom,
            context_ptr->txfm_context_array,
            candidate_buffer_ptr->candidate_ptr->tx_depth,
            context_ptr->md_rate_estimation_ptr);

        rate += tx_size_bits;
    }

    // Assign full cost
    *(candidate_buffer_ptr->full_cost_ptr) = RDCOST(lambda, rate, totalDistortion);

    candidate_buffer_ptr->full_lambda_rate = *candidate_buffer_ptr->full_cost_ptr - totalDistortion;
    coeffRate = *y_coeff_bits;
    candidate_buffer_ptr->full_cost_luma = RDCOST(lambda, lumaRate + *y_coeff_bits, luma_sse);

    return return_error;
}

/*********************************************************************************
* merge_skip_full_cost function is used to estimate the cost of an AMVPSkip candidate
* mode for full mode decisoion module.
*
*   @param *cu_ptr(input)
*       cu_ptr is the pointer of the target CU.
*   @param *candidate_buffer_ptr(input)
*       chromaBufferPtr is the buffer pointer of the candidate luma mode.
*   @param qp(input)
*       qp is the quantizer parameter.
*   @param luma_distortion (input)
*       luma_distortion is the inter condidate luma distortion.
*   @param lambda(input)
*       lambda is the Lagrange multiplier
**********************************************************************************/
EbErrorType  Av1MergeSkipFullCost(
    PictureControlSet                    *picture_control_set_ptr,
    ModeDecisionContext                  *context_ptr,
    struct ModeDecisionCandidateBuffer   *candidate_buffer_ptr,
    CodingUnit                           *cu_ptr,
    uint64_t                               *y_distortion,
    uint64_t                               *cb_distortion,
    uint64_t                               *cr_distortion,
    uint64_t                                lambda,
    uint64_t                               *y_coeff_bits,
    uint64_t                               *cb_coeff_bits,
    uint64_t                               *cr_coeff_bits,
    BlockSize                               bsize)
{
    UNUSED(bsize);
    UNUSED(context_ptr);
    UNUSED(picture_control_set_ptr);

    EbErrorType  return_error = EB_ErrorNone;
    uint64_t skipModeCtx = cu_ptr->skip_flag_context;
    uint64_t mergeRate = 0;
    uint64_t skipRate = 0;
    // Merge
    //uint64_t mergeChromaRate;
    uint64_t mergeDistortion;
    uint64_t merge_cost;
    //uint64_t mergeLumaCost;
    uint64_t mergeLumaSse;
    uint64_t mergeChromaSse;
    uint64_t coeffRate;
    //uint64_t lumaCoeffRate;

    // SKIP
    uint64_t skipDistortion;
    uint64_t skip_cost;
    //uint64_t skipLumaCost;

    // Luma and chroma transform size shift for the distortion
    uint64_t skipLumaSse;
    uint64_t skipChromaSse;

    uint64_t skipModeRate = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->skip_mode_fac_bits[skipModeCtx][1];

    // Coeff rate
    coeffRate = (*y_coeff_bits + *cb_coeff_bits + *cr_coeff_bits);

    // Compute Merge Cost
    mergeLumaSse = y_distortion[0] << AV1_COST_PRECISION;
    mergeChromaSse = (cb_distortion[0] + cr_distortion[0]) << AV1_COST_PRECISION;

    skipLumaSse = y_distortion[1] << AV1_COST_PRECISION;
    skipChromaSse = (cb_distortion[1] + cr_distortion[1]) << AV1_COST_PRECISION;

    // *Note - As in JCTVC-G1102, the JCT-VC uses the Mode Decision forumula where the chromaSse has been weighted
    //  CostMode = (luma_sse + wchroma * chromaSse) + lambdaSse * rateMode

    //if (picture_control_set_ptr->parent_pcs_ptr->pred_structure == EB_PRED_RANDOM_ACCESS) {
    //    // Random Access
    //    if (picture_control_set_ptr->temporal_layer_index == 0) {
    //        mergeChromaSse = (((mergeChromaSse * chroma_weight_factor_ra[qp]) + CHROMA_WEIGHT_OFFSET) >> CHROMA_WEIGHT_SHIFT);
    //    }
    //    else if (picture_control_set_ptr->temporal_layer_index < 3) {
    //        mergeChromaSse = (((mergeChromaSse * chroma_weight_factor_ra_qp_scaling_l1[qp]) + CHROMA_WEIGHT_OFFSET) >> CHROMA_WEIGHT_SHIFT);
    //    }
    //    else {
    //        mergeChromaSse = (((mergeChromaSse * chroma_weight_factor_ra_qp_scaling_l3[qp]) + CHROMA_WEIGHT_OFFSET) >> CHROMA_WEIGHT_SHIFT);
    //    }
    //}
    //else {
    //    // Low delay
    //    if (picture_control_set_ptr->temporal_layer_index == 0) {
    //        mergeChromaSse = (((mergeChromaSse * chroma_weight_factor_ld[qp]) + CHROMA_WEIGHT_OFFSET) >> CHROMA_WEIGHT_SHIFT);
    //    }
    //    else {
    //        mergeChromaSse = (((mergeChromaSse * chroma_weight_factor_ld_qp_scaling[qp]) + CHROMA_WEIGHT_OFFSET) >> CHROMA_WEIGHT_SHIFT);
    //    }
    //}

    // Add fast rate to get the total rate of the subject mode
    mergeRate += candidate_buffer_ptr->candidate_ptr->fast_luma_rate;
    mergeRate += candidate_buffer_ptr->candidate_ptr->fast_chroma_rate;

    mergeRate += coeffRate;

    if (candidate_buffer_ptr->candidate_ptr->block_has_coeff) {
        uint64_t tx_size_bits = estimate_tx_size_bits(
            picture_control_set_ptr,
            context_ptr->cu_origin_x,
            context_ptr->cu_origin_y,
            context_ptr->cu_ptr,
            context_ptr->blk_geom,
            context_ptr->txfm_context_array,
            candidate_buffer_ptr->candidate_ptr->tx_depth,
            context_ptr->md_rate_estimation_ptr);

        mergeRate += tx_size_bits;
    }

    mergeDistortion = (mergeLumaSse + mergeChromaSse);

    //merge_cost = mergeDistortion + (((lambda * coeffRate + lambda * mergeLumaRate + lambda_chroma * mergeChromaRate) + MD_OFFSET) >> MD_SHIFT);

    merge_cost = RDCOST(lambda, mergeRate, mergeDistortion);
    // mergeLumaCost = mergeLumaSse    + (((lambda * lumaCoeffRate + lambda * mergeLumaRate) + MD_OFFSET) >> MD_SHIFT);

    // *Note - As in JCTVC-G1102, the JCT-VC uses the Mode Decision forumula where the chromaSse has been weighted
    //  CostMode = (luma_sse + wchroma * chromaSse) + lambdaSse * rateMode

    //if (picture_control_set_ptr->parent_pcs_ptr->pred_structure == EB_PRED_RANDOM_ACCESS) {
    //    if (picture_control_set_ptr->temporal_layer_index == 0) {
    //        skipChromaSse = (((skipChromaSse * chroma_weight_factor_ra[qp]) + CHROMA_WEIGHT_OFFSET) >> CHROMA_WEIGHT_SHIFT);
    //    }
    //    else if (picture_control_set_ptr->temporal_layer_index < 3) {
    //        skipChromaSse = (((skipChromaSse * chroma_weight_factor_ra_qp_scaling_l1[qp]) + CHROMA_WEIGHT_OFFSET) >> CHROMA_WEIGHT_SHIFT);
    //    }
    //    else {
    //        skipChromaSse = (((skipChromaSse * chroma_weight_factor_ra_qp_scaling_l3[qp]) + CHROMA_WEIGHT_OFFSET) >> CHROMA_WEIGHT_SHIFT);
    //    }
    //}
    //else {
    //    // Low Delay
    //    if (picture_control_set_ptr->temporal_layer_index == 0) {
    //        skipChromaSse = (((skipChromaSse * chroma_weight_factor_ld[qp]) + CHROMA_WEIGHT_OFFSET) >> CHROMA_WEIGHT_SHIFT);
    //    }
    //    else {
    //        skipChromaSse = (((skipChromaSse * chroma_weight_factor_ld_qp_scaling[qp]) + CHROMA_WEIGHT_OFFSET) >> CHROMA_WEIGHT_SHIFT);
    //    }
    //}

    skipDistortion = skipLumaSse + skipChromaSse;
    skipRate = skipModeRate;
    skip_cost = RDCOST(lambda, skipRate, skipDistortion);

    // Assigne full cost
    *candidate_buffer_ptr->full_cost_ptr = (skip_cost <= merge_cost) ? skip_cost : merge_cost;

    uint64_t tempDistortion;
    tempDistortion = (skip_cost <= merge_cost) ? skipDistortion : mergeDistortion;
    candidate_buffer_ptr->full_lambda_rate = *candidate_buffer_ptr->full_cost_ptr - tempDistortion;
    *candidate_buffer_ptr->full_cost_merge_ptr = merge_cost;
    *candidate_buffer_ptr->full_cost_skip_ptr = skip_cost;
    // Assigne merge flag
    candidate_buffer_ptr->candidate_ptr->merge_flag = EB_TRUE;
    // Assigne skip flag

    candidate_buffer_ptr->candidate_ptr->skip_flag = (skip_cost <= merge_cost) ? EB_TRUE : EB_FALSE;

    //CHKN:  skip_flag context is not accurate as MD does not keep skip info in sync with EncDec.

    return return_error;
}
/*********************************************************************************
* av1_intra_full_cost function is used to estimate the cost of an intra candidate mode
* for full mode decisoion module.
*
*   @param *cu_ptr(input)
*       cu_ptr is the pointer of the target CU.
*   @param *candidate_buffer_ptr(input)
*       chromaBufferPtr is the buffer pointer of the candidate luma mode.
*   @param qp(input)
*       qp is the quantizer parameter.
*   @param luma_distortion (input)
*       luma_distortion is the intra condidate luma distortion.
*   @param lambda(input)
*       lambda is the Lagrange multiplier
**********************************************************************************/
EbErrorType av1_intra_full_cost(
    PictureControlSet                    *picture_control_set_ptr,
    ModeDecisionContext                  *context_ptr,
    struct ModeDecisionCandidateBuffer   *candidate_buffer_ptr,
    CodingUnit                           *cu_ptr,
    uint64_t                                 *y_distortion,
    uint64_t                                 *cb_distortion,
    uint64_t                                 *cr_distortion,
    uint64_t                                  lambda,
    uint64_t                                 *y_coeff_bits,
    uint64_t                                 *cb_coeff_bits,
    uint64_t                                 *cr_coeff_bits,
    BlockSize                              bsize)

{
    EbErrorType return_error = EB_ErrorNone;

    Av1FullCost(
        picture_control_set_ptr,
        context_ptr,
        candidate_buffer_ptr,
        cu_ptr,
        y_distortion,
        cb_distortion,
        cr_distortion,
        lambda,
        y_coeff_bits,
        cb_coeff_bits,
        cr_coeff_bits,
        bsize);

    return return_error;
}

/*********************************************************************************
* av1_inter_full_cost function is used to estimate the cost of an inter candidate mode
* for full mode decisoion module in inter frames.
*
*   @param *cu_ptr(input)
*       cu_ptr is the pointer of the target CU.
*   @param *candidate_buffer_ptr(input)
*       chromaBufferPtr is the buffer pointer of the candidate luma mode.
*   @param qp(input)
*       qp is the quantizer parameter.
*   @param luma_distortion (input)
*       luma_distortion is the inter condidate luma distortion.
*   @param lambda(input)
*       lambda is the Lagrange multiplier
**********************************************************************************/
EbErrorType av1_inter_full_cost(
    PictureControlSet                    *picture_control_set_ptr,
    ModeDecisionContext                  *context_ptr,
    struct ModeDecisionCandidateBuffer   *candidate_buffer_ptr,
    CodingUnit                           *cu_ptr,
    uint64_t                                 *y_distortion,
    uint64_t                                 *cb_distortion,
    uint64_t                                 *cr_distortion,
    uint64_t                                  lambda,
    uint64_t                                 *y_coeff_bits,
    uint64_t                                 *cb_coeff_bits,
    uint64_t                                 *cr_coeff_bits,
    BlockSize                              bsize
)
{
    EbErrorType  return_error = EB_ErrorNone;

    if (candidate_buffer_ptr->candidate_ptr->merge_flag == EB_TRUE) {
        Av1MergeSkipFullCost(
            picture_control_set_ptr,
            context_ptr,
            candidate_buffer_ptr,
            cu_ptr,
            y_distortion,
            cb_distortion,
            cr_distortion,
            lambda,
            y_coeff_bits,
            cb_coeff_bits,
            cr_coeff_bits,
            bsize);
    }
    else {
        Av1FullCost(
            picture_control_set_ptr,
            context_ptr,
            candidate_buffer_ptr,
            cu_ptr,
            y_distortion,
            cb_distortion,
            cr_distortion,
            lambda,
            y_coeff_bits,
            cb_coeff_bits,
            cr_coeff_bits,
            bsize);
    }
    return return_error;
}

/************************************************************
* Coding Loop Context Generation
************************************************************/
void coding_loop_context_generation(
    ModeDecisionContext      *context_ptr,
    CodingUnit               *cu_ptr,
    uint32_t                      cu_origin_x,
    uint32_t                      cu_origin_y,
    uint32_t                      sb_sz,
    NeighborArrayUnit        *skip_coeff_neighbor_array,
    NeighborArrayUnit        *inter_pred_dir_neighbor_array,
    NeighborArrayUnit        *ref_frame_type_neighbor_array,
    NeighborArrayUnit        *intra_luma_mode_neighbor_array,
    NeighborArrayUnit        *skip_flag_neighbor_array,
    NeighborArrayUnit        *mode_type_neighbor_array,
    NeighborArrayUnit        *leaf_depth_neighbor_array,
    NeighborArrayUnit       *leaf_partition_neighbor_array)
{
    (void)sb_sz;
    UNUSED(ref_frame_type_neighbor_array);
    uint32_t modeTypeLeftNeighborIndex = get_neighbor_array_unit_left_index(
        mode_type_neighbor_array,
        cu_origin_y);
    uint32_t modeTypeTopNeighborIndex = get_neighbor_array_unit_top_index(
        mode_type_neighbor_array,
        cu_origin_x);
    uint32_t leafDepthLeftNeighborIndex = get_neighbor_array_unit_left_index(
        leaf_depth_neighbor_array,
        cu_origin_y);
    uint32_t leafDepthTopNeighborIndex = get_neighbor_array_unit_top_index(
        leaf_depth_neighbor_array,
        cu_origin_x);
    uint32_t skipFlagLeftNeighborIndex = get_neighbor_array_unit_left_index(
        skip_flag_neighbor_array,
        cu_origin_y);
    uint32_t skipFlagTopNeighborIndex = get_neighbor_array_unit_top_index(
        skip_flag_neighbor_array,
        cu_origin_x);
    uint32_t intraLumaModeLeftNeighborIndex = get_neighbor_array_unit_left_index(
        intra_luma_mode_neighbor_array,
        cu_origin_y);
    uint32_t intraLumaModeTopNeighborIndex = get_neighbor_array_unit_top_index(
        intra_luma_mode_neighbor_array,
        cu_origin_x);

    uint32_t partition_left_neighbor_index = get_neighbor_array_unit_left_index(
        leaf_partition_neighbor_array,
        cu_origin_y);
    uint32_t partition_above_neighbor_index = get_neighbor_array_unit_top_index(
        leaf_partition_neighbor_array,
        cu_origin_x);

    // Intra Luma Neighbor Modes

    cu_ptr->prediction_unit_array->intra_luma_left_mode = (uint32_t)(
        (mode_type_neighbor_array->left_array[modeTypeLeftNeighborIndex] != INTRA_MODE) ? (uint32_t)DC_PRED :
        intra_luma_mode_neighbor_array->left_array[intraLumaModeLeftNeighborIndex]);

    cu_ptr->prediction_unit_array->intra_luma_top_mode = (uint32_t)(
        (mode_type_neighbor_array->top_array[modeTypeTopNeighborIndex] != INTRA_MODE) ? (uint32_t)DC_PRED :
        intra_luma_mode_neighbor_array->top_array[intraLumaModeTopNeighborIndex]);

    int32_t contextIndex;
    if (mode_type_neighbor_array->left_array[modeTypeLeftNeighborIndex] != (uint8_t)INVALID_MODE && mode_type_neighbor_array->top_array[modeTypeTopNeighborIndex] != (uint8_t)INVALID_MODE) {
        contextIndex = (mode_type_neighbor_array->left_array[modeTypeLeftNeighborIndex] == (uint8_t)INTRA_MODE && mode_type_neighbor_array->top_array[modeTypeTopNeighborIndex] == (uint8_t)INTRA_MODE) ? 3 :
            (mode_type_neighbor_array->left_array[modeTypeLeftNeighborIndex] == (uint8_t)INTRA_MODE || mode_type_neighbor_array->top_array[modeTypeTopNeighborIndex] == (uint8_t)INTRA_MODE) ? 1 : 0;
    }
    else  if (mode_type_neighbor_array->left_array[modeTypeLeftNeighborIndex] != (uint8_t)INVALID_MODE)
        contextIndex = (mode_type_neighbor_array->left_array[modeTypeLeftNeighborIndex] == (uint8_t)INTRA_MODE) ? 2 : 0;
    else if (mode_type_neighbor_array->top_array[modeTypeTopNeighborIndex] != (uint8_t)INVALID_MODE)
        contextIndex = (mode_type_neighbor_array->top_array[modeTypeTopNeighborIndex] == (uint8_t)INTRA_MODE) ? 2 : 0;
    else
        contextIndex = 0;
    cu_ptr->is_inter_ctx = contextIndex;
    //  if(cu_ptr->is_inter_ctx!=0) //
    //      printf("ctx:%i \n",cu_ptr->is_inter_ctx);

      //   Top Intra Mode Neighbor Array instead of a Full
      // Skip Flag Context
    cu_ptr->skip_flag_context =
        (mode_type_neighbor_array->left_array[modeTypeLeftNeighborIndex] == (uint8_t)INVALID_MODE) ? 0 :
        (skip_flag_neighbor_array->left_array[skipFlagLeftNeighborIndex] == EB_TRUE) ? 1 : 0;
    cu_ptr->skip_flag_context +=
        (mode_type_neighbor_array->top_array[modeTypeTopNeighborIndex] == (uint8_t)INVALID_MODE) ? 0 :
        (skip_flag_neighbor_array->top_array[skipFlagTopNeighborIndex] == EB_TRUE) ? 1 : 0;

    // Split Flag Context (neighbor info)
    context_ptr->md_local_cu_unit[cu_ptr->mds_idx].left_neighbor_mode = (uint32_t)(
        (mode_type_neighbor_array->left_array[modeTypeLeftNeighborIndex] != INTRA_MODE) ? (uint32_t)DC_PRED :
        intra_luma_mode_neighbor_array->left_array[intraLumaModeLeftNeighborIndex]);
    context_ptr->md_local_cu_unit[cu_ptr->mds_idx].left_neighbor_depth = leaf_depth_neighbor_array->left_array[leafDepthLeftNeighborIndex];
    context_ptr->md_local_cu_unit[cu_ptr->mds_idx].top_neighbor_mode = (uint32_t)(
        (mode_type_neighbor_array->top_array[modeTypeTopNeighborIndex] != INTRA_MODE) ? (uint32_t)DC_PRED :
        intra_luma_mode_neighbor_array->top_array[intraLumaModeTopNeighborIndex]);
    context_ptr->md_local_cu_unit[cu_ptr->mds_idx].top_neighbor_depth = leaf_depth_neighbor_array->top_array[leafDepthTopNeighborIndex];

    // Generate Partition context
    context_ptr->md_local_cu_unit[cu_ptr->mds_idx].above_neighbor_partition = (((PartitionContext*)leaf_partition_neighbor_array->top_array)[partition_above_neighbor_index].above == (int8_t)INVALID_NEIGHBOR_DATA) ?
        0 : ((PartitionContext*)leaf_partition_neighbor_array->top_array)[partition_above_neighbor_index].above;

    context_ptr->md_local_cu_unit[cu_ptr->mds_idx].left_neighbor_partition = (((PartitionContext*)leaf_partition_neighbor_array->left_array)[partition_left_neighbor_index].left == (int8_t)INVALID_NEIGHBOR_DATA) ?
        0 : ((PartitionContext*)leaf_partition_neighbor_array->left_array)[partition_left_neighbor_index].left;
    // Skip Coeff AV1 Context
    uint32_t skipCoeffLeftNeighborIndex = get_neighbor_array_unit_left_index(
        skip_coeff_neighbor_array,
        cu_origin_y);
    uint32_t skipCoeffTopNeighborIndex = get_neighbor_array_unit_top_index(
        skip_coeff_neighbor_array,
        cu_origin_x);

    cu_ptr->skip_coeff_context =
        (skip_coeff_neighbor_array->left_array[skipCoeffLeftNeighborIndex] == (uint8_t)INVALID_NEIGHBOR_DATA) ? 0 :
        (skip_coeff_neighbor_array->left_array[skipCoeffLeftNeighborIndex]) ? 1 : 0;

    cu_ptr->skip_coeff_context +=
        (skip_coeff_neighbor_array->top_array[skipCoeffTopNeighborIndex] == (uint8_t)INVALID_NEIGHBOR_DATA) ? 0 :
        (skip_coeff_neighbor_array->top_array[skipCoeffTopNeighborIndex]) ? 1 : 0;
    // Generate reference mode context

    cu_ptr->reference_mode_context = (uint8_t)av1_get_reference_mode_context(
        cu_origin_x,
        cu_origin_y,
        mode_type_neighbor_array,
        inter_pred_dir_neighbor_array);

    cu_ptr->compoud_reference_type_context = (uint8_t)av1_get_comp_reference_type_context(
        cu_origin_x,
        cu_origin_y,
        mode_type_neighbor_array,
        inter_pred_dir_neighbor_array);

    //Collect Neighbor ref cout
    av1_collect_neighbors_ref_counts_new(cu_ptr->av1xd);

    return;
}

/********************************************
* tu_calc_cost
*   computes TU Cost and generetes TU Cbf
********************************************/
EbErrorType av1_tu_calc_cost(
    ModeDecisionCandidate *candidate_ptr,                        // input parameter, prediction result Ptr
    int16_t                   txb_skip_ctx,
    uint32_t                   tu_index,                             // input parameter, TU index inside the CU
    uint32_t                   y_count_non_zero_coeffs,                 // input parameter, number of non zero Y quantized coefficients
    uint32_t                   cb_count_non_zero_coeffs,                // input parameter, number of non zero cb quantized coefficients
    uint32_t                   cr_count_non_zero_coeffs,                // input parameter, number of non zero cr quantized coefficients
    uint64_t                   y_tu_distortion[DIST_CALC_TOTAL],      // input parameter, Y distortion for both Normal and Cbf zero modes
    uint64_t                   cb_tu_distortion[DIST_CALC_TOTAL],     // input parameter, Cb distortion for both Normal and Cbf zero modes
    uint64_t                   cr_tu_distortion[DIST_CALC_TOTAL],     // input parameter, Cr distortion for both Normal and Cbf zero modes
    COMPONENT_TYPE           component_type,
    uint64_t                  *y_tu_coeff_bits,                        // input parameter, Y quantized coefficients rate
    uint64_t                  *cb_tu_coeff_bits,                       // input parameter, Cb quantized coefficients rate
    uint64_t                  *cr_tu_coeff_bits,                       // input parameter, Cr quantized coefficients rate
    TxSize                  txsize,
    uint64_t                   lambda)                              // input parameter, lambda for Luma

{
    (void)cr_tu_coeff_bits;
    (void)cb_tu_coeff_bits;
    (void)cr_tu_distortion;
    (void)cb_tu_distortion;
    EbErrorType return_error = EB_ErrorNone;
    // Non Zero coeff mode variables
    uint64_t y_nonzero_coeff_distortion = y_tu_distortion[DIST_CALC_RESIDUAL];
    uint64_t y_nonzero_coeff_rate;

    uint64_t y_nonzero_coeff_cost = 0;

    // Zero Cbf mode variables
    uint64_t y_zero_coeff_distortion = y_tu_distortion[DIST_CALC_PREDICTION];

    uint64_t y_zero_coeff_luma_flag_bits_num = 0;

    uint64_t y_zero_coeff_rate;

    uint64_t y_zero_coeff_cost = 0;
    if (component_type == COMPONENT_LUMA || component_type == COMPONENT_ALL) {
        // Non Zero Distortion
        // *Note - As of Oct 2011, the JCT-VC uses the PSNR forumula
        //  PSNR = (LUMA_WEIGHT * PSNRy + PSNRu + PSNRv) / (2+LUMA_WEIGHT)
        y_nonzero_coeff_distortion = LUMA_WEIGHT * (y_nonzero_coeff_distortion << AV1_COST_PRECISION);

        // Zero distortion
        // *Note - As of Oct 2011, the JCT-VC uses the PSNR forumula
        //  PSNR = (LUMA_WEIGHT * PSNRy + PSNRu + PSNRv) / (2+LUMA_WEIGHT)
        y_zero_coeff_distortion = LUMA_WEIGHT * (y_zero_coeff_distortion << AV1_COST_PRECISION);

        // **Compute Rate

        // Esimate Cbf's Bits

        const TxSize txs_ctx = (TxSize)((txsize_sqr_map[txsize] + txsize_sqr_up_map[txsize] + 1) >> 1);
        assert(txs_ctx < TX_SIZES);
        const LvMapCoeffCost *const coeff_costs = &candidate_ptr->md_rate_estimation_ptr->coeff_fac_bits[txs_ctx][0];

        y_zero_coeff_luma_flag_bits_num = coeff_costs->txb_skip_cost[txb_skip_ctx][1];

        y_nonzero_coeff_rate = *y_tu_coeff_bits; // yNonZeroCbfLumaFlagBitsNum is already calculated inside y_tu_coeff_bits

        y_zero_coeff_rate = y_zero_coeff_luma_flag_bits_num;

        if (1)
            y_zero_coeff_cost = 0xFFFFFFFFFFFFFFFFull;
        else
            y_zero_coeff_cost = RDCOST(lambda, y_zero_coeff_rate, y_zero_coeff_distortion);
        // **Compute Cost
        y_nonzero_coeff_cost = RDCOST(lambda, y_nonzero_coeff_rate, y_nonzero_coeff_distortion);

        candidate_ptr->y_has_coeff |= (((y_count_non_zero_coeffs != 0) && (y_nonzero_coeff_cost < y_zero_coeff_cost)) << tu_index);
        *y_tu_coeff_bits = (y_nonzero_coeff_cost < y_zero_coeff_cost) ? *y_tu_coeff_bits : 0;
        y_tu_distortion[DIST_CALC_RESIDUAL] = (y_nonzero_coeff_cost < y_zero_coeff_cost) ? y_tu_distortion[DIST_CALC_RESIDUAL] : y_tu_distortion[DIST_CALC_PREDICTION];
        }
    if (component_type == COMPONENT_CHROMA_CB || component_type == COMPONENT_CHROMA || component_type == COMPONENT_ALL)
        candidate_ptr->u_has_coeff |= ((cb_count_non_zero_coeffs != 0) << tu_index);
    if (component_type == COMPONENT_CHROMA_CR || component_type == COMPONENT_CHROMA || component_type == COMPONENT_ALL)
        candidate_ptr->v_has_coeff |= ((cr_count_non_zero_coeffs != 0) << tu_index);
    return return_error;
    }

/********************************************
* tu_calc_cost
*   computes TU Cost and generetes TU Cbf
********************************************/

EbErrorType av1_tu_calc_cost_luma(
    int16_t                   txb_skip_ctx,
    ModeDecisionCandidate *candidate_ptr,                        // input parameter, prediction result Ptr
    uint32_t                   tu_index,                             // input parameter, TU index inside the CU
    TxSize                  tx_size,
    uint32_t                   y_count_non_zero_coeffs,                 // input parameter, number of non zero Y quantized coefficients
    uint64_t                   y_tu_distortion[DIST_CALC_TOTAL],      // input parameter, Y distortion for both Normal and Cbf zero modes
    uint64_t                  *y_tu_coeff_bits,                        // input parameter, Y quantized coefficients rate
    uint64_t                  *y_full_cost,
    uint64_t                   lambda)                              // input parameter, lambda for Luma

{
    EbErrorType return_error = EB_ErrorNone;

    // Non Zero Cbf mode variables
    uint64_t yNonZeroCbfDistortion = y_tu_distortion[DIST_CALC_RESIDUAL];

    uint64_t yNonZeroCbfRate;

    uint64_t yNonZeroCbfCost = 0;

    // Zero Cbf mode variables
    uint64_t yZeroCbfDistortion = y_tu_distortion[DIST_CALC_PREDICTION];

    uint64_t yZeroCbfLumaFlagBitsNum = 0;

    uint64_t yZeroCbfRate;

    uint64_t yZeroCbfCost = 0;

    // **Compute distortion
    // Non Zero Distortion
    // *Note - As of Oct 2011, the JCT-VC uses the PSNR forumula
    //  PSNR = (LUMA_WEIGHT * PSNRy + PSNRu + PSNRv) / (2+LUMA_WEIGHT)
    yNonZeroCbfDistortion = LUMA_WEIGHT * (yNonZeroCbfDistortion << AV1_COST_PRECISION);

    // Zero distortion
    // *Note - As of Oct 2011, the JCT-VC uses the PSNR forumula
    //  PSNR = (LUMA_WEIGHT * PSNRy + PSNRu + PSNRv) / (2+LUMA_WEIGHT)
    yZeroCbfDistortion = LUMA_WEIGHT * (yZeroCbfDistortion << AV1_COST_PRECISION);

    // **Compute Rate

    // Esimate Cbf's Bits

    const TxSize txs_ctx = (TxSize)((txsize_sqr_map[tx_size] + txsize_sqr_up_map[tx_size] + 1) >> 1);
    assert(txs_ctx < TX_SIZES);
    const LvMapCoeffCost *const coeff_costs = &candidate_ptr->md_rate_estimation_ptr->coeff_fac_bits[txs_ctx][0];

    yZeroCbfLumaFlagBitsNum = coeff_costs->txb_skip_cost[txb_skip_ctx][1];

    yNonZeroCbfRate = *y_tu_coeff_bits; // yNonZeroCbfLumaFlagBitsNum is already calculated inside y_tu_coeff_bits

    yZeroCbfRate = yZeroCbfLumaFlagBitsNum;

    if (1)
        yZeroCbfCost = 0xFFFFFFFFFFFFFFFFull;
    else
        yZeroCbfCost = RDCOST(lambda, yZeroCbfRate, yZeroCbfDistortion);
    // **Compute Cost
    yNonZeroCbfCost = RDCOST(lambda, yNonZeroCbfRate, yNonZeroCbfDistortion);
    candidate_ptr->y_has_coeff |= ((y_count_non_zero_coeffs != 0) << tu_index);
    *y_tu_coeff_bits = (yNonZeroCbfCost < yZeroCbfCost) ? *y_tu_coeff_bits : 0;
    y_tu_distortion[DIST_CALC_RESIDUAL] = (yNonZeroCbfCost < yZeroCbfCost) ? y_tu_distortion[DIST_CALC_RESIDUAL] : y_tu_distortion[DIST_CALC_PREDICTION];

    *y_full_cost = MIN(yNonZeroCbfCost, yZeroCbfCost);

    return return_error;
    }

//static INLINE int32_t partition_plane_context(const MacroBlockD *xd, int32_t mi_row,
//    int32_t mi_col, BlockSize bsize) {
//    const PartitionContextType *above_ctx = xd->above_seg_context + mi_col;
//    const PartitionContextType *left_ctx =
//        xd->left_seg_context + (mi_row & MAX_MIB_MASK);
//    // Minimum partition point is 8x8. Offset the bsl accordingly.
//    const int32_t bsl = mi_size_wide_log2[bsize] - mi_size_wide_log2[BLOCK_8X8];
//    int32_t above = (*above_ctx >> bsl) & 1, left = (*left_ctx >> bsl) & 1;
//
//    assert(mi_size_wide_log2[bsize] == mi_size_high_log2[bsize]);
//    assert(bsl >= 0);
//
//    return (left * 2 + above) + bsl * PARTITION_PLOFFSET;
//}

/*********************************************************************************
* split_flag_rate function is used to generate the Split rate
*
*   @param *cu_ptr(input)
*       cu_ptr is the pointer of the target CU.
*   @param split_flag(input)
*       split_flag is the split flag value.
*   @param split_rate(output)
*       split_rate contains rate.
*   @param lambda(input)
*       lambda is the Lagrange multiplier
*   @param md_rate_estimation_ptr(input)
*       md_rate_estimation_ptr is pointer to MD rate Estimation Tables
**********************************************************************************/
EbErrorType av1_split_flag_rate(
    SequenceControlSet                  *sequence_control_set_ptr,
    ModeDecisionContext                  *context_ptr,
    CodingUnit                           *cu_ptr,
    uint32_t                                  leaf_index,
    PartitionType                          partitionType,
    uint64_t                                 *split_rate,
    uint64_t                                  lambda,
    MdRateEstimationContext              *md_rate_estimation_ptr,
    uint32_t                                  tb_max_depth)
{
    (void)tb_max_depth;
    (void)leaf_index;

    const BlockGeom          *blk_geom = get_blk_geom_mds(cu_ptr->mds_idx);
    EbErrorType return_error = EB_ErrorNone;

    uint32_t cu_origin_x = context_ptr->sb_origin_x + blk_geom->origin_x;
    uint32_t cu_origin_y = context_ptr->sb_origin_y + blk_geom->origin_y;

    PartitionType p = partitionType;

    uint32_t cu_depth = blk_geom->depth;
    UNUSED(cu_depth);
    BlockSize bsize = blk_geom->bsize;
    assert(bsize<BlockSizeS_ALL);
    const int32_t is_partition_point = blk_geom->bsize >= BLOCK_8X8;

    if (is_partition_point) {
        const int32_t hbs = (mi_size_wide[bsize] << 2) >> 1;
        const int32_t hasRows = (cu_origin_y + hbs) < sequence_control_set_ptr->seq_header.max_frame_height;
        const int32_t hasCols = (cu_origin_x + hbs) < sequence_control_set_ptr->seq_header.max_frame_width;

        uint32_t contextIndex = 0;

        const PartitionContextType left_ctx = context_ptr->md_local_cu_unit[cu_ptr->mds_idx].left_neighbor_partition == (int8_t)(INVALID_NEIGHBOR_DATA) ? 0 : context_ptr->md_local_cu_unit[cu_ptr->mds_idx].left_neighbor_partition;
        const PartitionContextType above_ctx = context_ptr->md_local_cu_unit[cu_ptr->mds_idx].above_neighbor_partition == (int8_t)(INVALID_NEIGHBOR_DATA) ? 0 : context_ptr->md_local_cu_unit[cu_ptr->mds_idx].above_neighbor_partition;

        const int32_t bsl = mi_size_wide_log2[bsize] - mi_size_wide_log2[BLOCK_8X8];

        int32_t above = (above_ctx >> bsl) & 1, left = (left_ctx >> bsl) & 1;

        assert(mi_size_wide_log2[bsize] == mi_size_high_log2[bsize]);
        assert(bsl >= 0);

        contextIndex = (left * 2 + above) + bsl * PARTITION_PLOFFSET;

        if (hasRows && hasCols) {
            *split_rate = (uint64_t)md_rate_estimation_ptr->partition_fac_bits[contextIndex][partitionType];

        }
        else if (!hasRows && hasCols) {
            *split_rate = (uint64_t)md_rate_estimation_ptr->partition_fac_bits[2][p == PARTITION_SPLIT];

        }
        else {
            *split_rate = (uint64_t)md_rate_estimation_ptr->partition_fac_bits[2][p == PARTITION_SPLIT];

        }
    }
    else
        *split_rate = (uint64_t)md_rate_estimation_ptr->partition_fac_bits[0][partitionType];
    *split_rate = RDCOST(lambda, *split_rate, 0);

    return return_error;
}

/********************************************
* tu_calc_cost
*   Computes TU Cost and generetes TU Cbf
*   at the level of the encode pass
********************************************/
EbErrorType av1_encode_tu_calc_cost(
    EncDecContext          *context_ptr,
    uint32_t                   *count_non_zero_coeffs,
    uint64_t                    y_tu_distortion[DIST_CALC_TOTAL],
    uint64_t                   *y_tu_coeff_bits,
    uint32_t                    component_mask
)
{
    CodingUnit              *cu_ptr = context_ptr->cu_ptr;
    uint32_t                     tu_index = context_ptr->txb_itr;
    MdRateEstimationContext *md_rate_estimation_ptr = context_ptr->md_rate_estimation_ptr;
    uint64_t                     lambda = context_ptr->full_lambda;
    uint32_t                     y_count_non_zero_coeffs = count_non_zero_coeffs[0];
    uint32_t                     cb_count_non_zero_coeffs = count_non_zero_coeffs[1];
    uint32_t                     cr_count_non_zero_coeffs = count_non_zero_coeffs[2];

    EbErrorType return_error = EB_ErrorNone;

    // Non Zero Cbf mode variables
    uint64_t yNonZeroCbfDistortion = y_tu_distortion[DIST_CALC_RESIDUAL];

    uint64_t yNonZeroCbfRate;

    uint64_t yNonZeroCbfCost = 0;

    // Zero Cbf mode variables
    uint64_t yZeroCbfDistortion = y_tu_distortion[DIST_CALC_PREDICTION];

    uint64_t yZeroCbfLumaFlagBitsNum = 0;

    uint64_t yZeroCbfRate;

    uint64_t yZeroCbfCost = 0;
    int16_t  txb_skip_ctx = context_ptr->md_context->luma_txb_skip_context;
    // **Compute distortion
    if (component_mask == PICTURE_BUFFER_DESC_LUMA_MASK || component_mask == PICTURE_BUFFER_DESC_FULL_MASK) {
        // Non Zero Distortion
        // *Note - As of Oct 2011, the JCT-VC uses the PSNR forumula
        //  PSNR = (LUMA_WEIGHT * PSNRy + PSNRu + PSNRv) / (2+LUMA_WEIGHT)
        yNonZeroCbfDistortion = LUMA_WEIGHT * (yNonZeroCbfDistortion << AV1_COST_PRECISION);

        // Zero distortion
        // *Note - As of Oct 2011, the JCT-VC uses the PSNR forumula
        //  PSNR = (LUMA_WEIGHT * PSNRy + PSNRu + PSNRv) / (2+LUMA_WEIGHT)
        yZeroCbfDistortion = LUMA_WEIGHT * (yZeroCbfDistortion << AV1_COST_PRECISION);
        TxSize    txSize = context_ptr->blk_geom->txsize[cu_ptr->tx_depth][context_ptr->txb_itr];
        assert(txSize < TX_SIZES_ALL);

        const TxSize txs_ctx = (TxSize)((txsize_sqr_map[txSize] + txsize_sqr_up_map[txSize] + 1) >> 1);
        assert(txs_ctx < TX_SIZES);
        const LvMapCoeffCost *const coeff_costs = &md_rate_estimation_ptr->coeff_fac_bits[txs_ctx][0];

        yZeroCbfLumaFlagBitsNum = coeff_costs->txb_skip_cost[txb_skip_ctx][1];

        yNonZeroCbfRate = *y_tu_coeff_bits; // yNonZeroCbfLumaFlagBitsNum is already calculated inside y_tu_coeff_bits

        yZeroCbfRate = yZeroCbfLumaFlagBitsNum;
        TransformUnit       *txb_ptr = &cu_ptr->transform_unit_array[context_ptr->txb_itr];
        if (txb_ptr->transform_type[PLANE_TYPE_Y] != DCT_DCT) {
            yZeroCbfCost = 0xFFFFFFFFFFFFFFFFull;
        }
        else
            yZeroCbfCost = RDCOST(lambda, yZeroCbfRate, yZeroCbfDistortion);
        // **Compute Cost
        yNonZeroCbfCost = RDCOST(lambda, yNonZeroCbfRate, yNonZeroCbfDistortion);
        cu_ptr->transform_unit_array[tu_index].y_has_coeff = ((y_count_non_zero_coeffs != 0) && (yNonZeroCbfCost < yZeroCbfCost)) ? EB_TRUE : EB_FALSE;
        *y_tu_coeff_bits = (yNonZeroCbfCost < yZeroCbfCost) ? *y_tu_coeff_bits : 0;
        y_tu_distortion[DIST_CALC_RESIDUAL] = (yNonZeroCbfCost < yZeroCbfCost) ? y_tu_distortion[DIST_CALC_RESIDUAL] : y_tu_distortion[DIST_CALC_PREDICTION];
        }
    else
        cu_ptr->transform_unit_array[tu_index].y_has_coeff = EB_FALSE;
    cu_ptr->transform_unit_array[tu_index].u_has_coeff = cb_count_non_zero_coeffs != 0 ? EB_TRUE : EB_FALSE;
    cu_ptr->transform_unit_array[tu_index].v_has_coeff = cr_count_non_zero_coeffs != 0 ? EB_TRUE : EB_FALSE;

    return return_error;
    }

uint64_t GetPMCost(
    uint64_t                   lambda,
    uint64_t                   tuDistortion,
    uint64_t                   y_tu_coeff_bits
)
{
    uint64_t yNonZeroCbfDistortion = LUMA_WEIGHT * (tuDistortion << COST_PRECISION);
    uint64_t yNonZeroCbfRate = (y_tu_coeff_bits);
    uint64_t yNonZeroCbfCost = yNonZeroCbfDistortion + (((lambda       * yNonZeroCbfRate) + MD_OFFSET) >> MD_SHIFT);

    return yNonZeroCbfCost;
}
