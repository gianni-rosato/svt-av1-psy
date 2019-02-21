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

BlockSize GetBlockSize(uint8_t cu_size) {
    return (cu_size == 64 ? BLOCK_64X64 : cu_size == 32 ? BLOCK_32X32 : cu_size == 16 ? BLOCK_16X16 : cu_size == 8 ? BLOCK_8X8 : BLOCK_4X4);
}

static INLINE int32_t is_chroma_reference(int32_t mi_row, int32_t mi_col, BlockSize bsize,
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
//} MV_JOINT_TYPE;



MV_JOINT_TYPE av1_get_mv_joint(const MV *mv) {
    if (mv->row == 0) {
        return mv->col == 0 ? MV_JOINT_ZERO : MV_JOINT_HNZVZ;
    }
    else {
        return mv->col == 0 ? MV_JOINT_HZVNZ : MV_JOINT_HNZVNZ;
    }
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
    return tx_size_wide_log2[tx_size];
}

static INLINE int32_t get_txb_wide(TxSize tx_size) {
    tx_size = av1_get_adjusted_tx_size(tx_size);
    return tx_size_wide[tx_size];
}

static INLINE int32_t get_txb_high(TxSize tx_size) {
    tx_size = av1_get_adjusted_tx_size(tx_size);
    return tx_size_high[tx_size];
}
static INLINE uint8_t *set_levels(uint8_t *const levels_buf, const int32_t width) {
    return levels_buf + TX_PAD_TOP * (width + TX_PAD_HOR);
}

void av1_txb_init_levels_c(
    const tran_low_t *const coeff,
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
        for (int32_t j = 0; j < width; j++) {

            *ls++ = (uint8_t)clamp(abs(coeff[i * width + j]), 0, INT8_MAX);

        }
        for (int32_t j = 0; j < TX_PAD_HOR; j++) {
            *ls++ = 0;
        }
    }
}


static const PredictionMode fimode_to_intradir[FILTER_INTRA_MODES] = {
    DC_PRED, V_PRED, H_PRED, D157_PRED, DC_PRED
};
// TODO(angiebird): use this function whenever it's possible
int32_t Av1TransformTypeRateEstimation(
    struct ModeDecisionCandidateBuffer_s    *candidate_buffer_ptr,
    EbBool                                  is_inter,
    EbBool                                  useFilterIntraFlag,
    TxSize                                  transform_size,
    TxType                                  transform_type,
    EbBool                                  reduced_tx_set_used)
{


    uint8_t filterIntraMode = 0; // NM- hardcoded to zero for the moment until we support different intra filtering modes.
    const TxSize square_tx_size = txsize_sqr_map[transform_size];

    //const MbModeInfo *mbmi = &xd->mi[0]->mbmi;
    //const int32_t is_inter = is_inter_block(mbmi);

    if (get_ext_tx_types(transform_size, is_inter, reduced_tx_set_used) > 1  /*&&    !xd->lossless[xd->mi[0]->mbmi.segment_id]  WE ARE NOT LOSSLESS*/) {

        const int32_t ext_tx_set = get_ext_tx_set(transform_size, is_inter, reduced_tx_set_used);
        if (is_inter) {
            if (ext_tx_set > 0)
                return candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->interTxTypeFacBits[ext_tx_set][square_tx_size][transform_type];
        }
        else {
            if (ext_tx_set > 0) {
                PredictionMode intra_dir;
                if (useFilterIntraFlag)
                    intra_dir = fimode_to_intradir[filterIntraMode];
                else
                    intra_dir = candidate_buffer_ptr->candidate_ptr->pred_mode;
                ASSERT(intra_dir < INTRA_MODES);
                return candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->intraTxTypeFacBits[ext_tx_set][square_tx_size][intra_dir][transform_type];
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

    if (eob < 33) {
        t = eob_to_pos_small[eob];
    }
    else {
        const int32_t e = AOMMIN((eob - 1) >> 5, 16);
        t = eob_to_pos_large[e];
    }

    *extra = eob - k_eob_group_start[t];

    return t;
}

static int32_t get_eob_cost(int32_t eob, const LV_MAP_EOB_COST *txb_eob_costs,
    const LV_MAP_COEFF_COST *txb_costs, TxType tx_type) {
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

// The ctx offset table when TX is TX_CLASS_2D.
// TX col and row indices are clamped to 4.
const int8_t av1_nz_map_ctx_offset[TX_SIZES_ALL][5][5] = {
    // TX_4X4
    { { 0, 1, 6, 6, 0 },
    { 1, 6, 6, 21, 0 },
    { 6, 6, 21, 21, 0 },
    { 6, 21, 21, 21, 0 },
    { 0, 0, 0, 0, 0 } },
    // TX_8X8
    { { 0, 1, 6, 6, 21 },
    { 1, 6, 6, 21, 21 },
    { 6, 6, 21, 21, 21 },
    { 6, 21, 21, 21, 21 },
    { 21, 21, 21, 21, 21 } },
    // TX_16X16
    { { 0, 1, 6, 6, 21 },
    { 1, 6, 6, 21, 21 },
    { 6, 6, 21, 21, 21 },
    { 6, 21, 21, 21, 21 },
    { 21, 21, 21, 21, 21 } },
    // TX_32X32
    { { 0, 1, 6, 6, 21 },
    { 1, 6, 6, 21, 21 },
    { 6, 6, 21, 21, 21 },
    { 6, 21, 21, 21, 21 },
    { 21, 21, 21, 21, 21 } },
    // TX_64X64
    { { 0, 1, 6, 6, 21 },
    { 1, 6, 6, 21, 21 },
    { 6, 6, 21, 21, 21 },
    { 6, 21, 21, 21, 21 },
    { 21, 21, 21, 21, 21 } },
    // TX_4X8
    { { 0, 11, 11, 11, 0 },
    { 11, 11, 11, 11, 0 },
    { 6, 6, 21, 21, 0 },
    { 6, 21, 21, 21, 0 },
    { 21, 21, 21, 21, 0 } },
    // TX_8X4
    { { 0, 16, 6, 6, 21 },
    { 16, 16, 6, 21, 21 },
    { 16, 16, 21, 21, 21 },
    { 16, 16, 21, 21, 21 },
    { 0, 0, 0, 0, 0 } },
    // TX_8X16
    { { 0, 11, 11, 11, 11 },
    { 11, 11, 11, 11, 11 },
    { 6, 6, 21, 21, 21 },
    { 6, 21, 21, 21, 21 },
    { 21, 21, 21, 21, 21 } },
    // TX_16X8
    { { 0, 16, 6, 6, 21 },
    { 16, 16, 6, 21, 21 },
    { 16, 16, 21, 21, 21 },
    { 16, 16, 21, 21, 21 },
    { 16, 16, 21, 21, 21 } },
    // TX_16X32
    { { 0, 11, 11, 11, 11 },
    { 11, 11, 11, 11, 11 },
    { 6, 6, 21, 21, 21 },
    { 6, 21, 21, 21, 21 },
    { 21, 21, 21, 21, 21 } },
    // TX_32X16
    { { 0, 16, 6, 6, 21 },
    { 16, 16, 6, 21, 21 },
    { 16, 16, 21, 21, 21 },
    { 16, 16, 21, 21, 21 },
    { 16, 16, 21, 21, 21 } },
    // TX_32X64
    { { 0, 11, 11, 11, 11 },
    { 11, 11, 11, 11, 11 },
    { 6, 6, 21, 21, 21 },
    { 6, 21, 21, 21, 21 },
    { 21, 21, 21, 21, 21 } },
    // TX_64X32
    { { 0, 16, 6, 6, 21 },
    { 16, 16, 6, 21, 21 },
    { 16, 16, 21, 21, 21 },
    { 16, 16, 21, 21, 21 },
    { 16, 16, 21, 21, 21 } },
    // TX_4X16
    { { 0, 11, 11, 11, 0 },
    { 11, 11, 11, 11, 0 },
    { 6, 6, 21, 21, 0 },
    { 6, 21, 21, 21, 0 },
    { 21, 21, 21, 21, 0 } },
    // TX_16X4
    { { 0, 16, 6, 6, 21 },
    { 16, 16, 6, 21, 21 },
    { 16, 16, 21, 21, 21 },
    { 16, 16, 21, 21, 21 },
    { 0, 0, 0, 0, 0 } },
    // TX_8X32
    { { 0, 11, 11, 11, 11 },
    { 11, 11, 11, 11, 11 },
    { 6, 6, 21, 21, 21 },
    { 6, 21, 21, 21, 21 },
    { 21, 21, 21, 21, 21 } },
    // TX_32X8
    { { 0, 16, 6, 6, 21 },
    { 16, 16, 6, 21, 21 },
    { 16, 16, 21, 21, 21 },
    { 16, 16, 21, 21, 21 },
    { 16, 16, 21, 21, 21 } },
    // TX_16X64
    { { 0, 11, 11, 11, 11 },
    { 11, 11, 11, 11, 11 },
    { 6, 6, 21, 21, 21 },
    { 6, 21, 21, 21, 21 },
    { 21, 21, 21, 21, 21 } },
    // TX_64X16
    { { 0, 16, 6, 6, 21 },
    { 16, 16, 6, 21, 21 },
    { 16, 16, 21, 21, 21 },
    { 16, 16, 21, 21, 21 },
    { 16, 16, 21, 21, 21 } }
};

static INLINE int32_t get_br_ctx(const uint8_t *const levels,
    const int32_t c,  // raster order
    const int32_t bwl, const TxType tx_type) {
    const int32_t row = c >> bwl;
    const int32_t col = c - (row << bwl);
    const int32_t stride = (1 << bwl) + TX_PAD_HOR;
    const TX_CLASS tx_class = tx_type_to_class[tx_type];
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
    struct ModeDecisionCandidateBuffer_s    *candidate_buffer_ptr,
    TxSize                                  transform_size,
    PLANE_TYPE                               plane_type,
    int16_t                                   txb_skip_ctx)
{
    const TxSize txs_ctx = (TxSize)((txsize_sqr_map[transform_size] + txsize_sqr_up_map[transform_size] + 1) >> 1);
    const LV_MAP_COEFF_COST *const coeff_costs = &candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->coeffFacBits[txs_ctx][plane_type];
    return coeff_costs->txb_skip_cost[txb_skip_ctx][1];
}
// Note: don't call this function when eob is 0.
uint64_t av1_cost_coeffs_txb(
    struct ModeDecisionCandidateBuffer_s    *candidate_buffer_ptr,
    const tran_low_t                        *const qcoeff,
    uint16_t                                   eob,
    PLANE_TYPE                               plane_type,
    TxSize                                  transform_size,
    /*const uint32_t                             areaSize,
    const uint32_t                             stride,*/
    int16_t                                   txb_skip_ctx,
    int16_t                                   dc_sign_ctx,
    EbBool                                  reducedTransformSetFlag)

{

    //Note: there is a different version of this function in AOM that seems to be efficient as its name is:
    //warehouse_efficients_txb

    const TxSize txs_ctx = (TxSize)((txsize_sqr_map[transform_size] + txsize_sqr_up_map[transform_size] + 1) >> 1);
    const TxType transform_type = candidate_buffer_ptr->candidate_ptr->transform_type[plane_type];
    const TX_CLASS tx_class = tx_type_to_class[transform_type];
    int32_t c, cost;
    const int32_t bwl = get_txb_bwl(transform_size);
    const int32_t width = get_txb_wide(transform_size);
    const int32_t height = get_txb_high(transform_size);
    const SCAN_ORDER *const scan_order = &av1_scan_orders[transform_size][transform_type]; // get_scan(tx_size, tx_type);
    const int16_t *const scan = scan_order->scan;
    uint8_t levels_buf[TX_PAD_2D];
    uint8_t *const levels = set_levels(levels_buf, width);
    DECLARE_ALIGNED(16, int8_t, coeff_contexts[MAX_TX_SQUARE]);
    const LV_MAP_COEFF_COST *const coeff_costs = &candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->coeffFacBits[txs_ctx][plane_type];

    const int32_t eob_multi_size = txsize_log2_minus4[transform_size];
    const LV_MAP_EOB_COST *const eobBits = &candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->eobFracBits[eob_multi_size][plane_type];
    // eob must be greater than 0 here.
    assert(eob > 0);
    cost = coeff_costs->txb_skip_cost[txb_skip_ctx][0];


    av1_txb_init_levels(qcoeff, width, height, levels); // NM - Needs to be optimized - to be combined with the quantisation.


    // Transform type bit estimation
    cost += plane_type > PLANE_TYPE_Y ? 0 :
        Av1TransformTypeRateEstimation(
            candidate_buffer_ptr,
            candidate_buffer_ptr->candidate_ptr->type == INTER_MODE ? EB_TRUE : EB_FALSE,
            EB_FALSE, // NM - Hardcoded to false for the moment until we support the intra filtering
            transform_size,
            transform_type,
            reducedTransformSetFlag);

    // Transform ebo bit estimation
    int32_t eob_cost = get_eob_cost(eob, eobBits, coeff_costs, transform_type);
    cost += eob_cost;

    // Transform non-zero coeff bit estimation
    av1_get_nz_map_contexts(
        levels,
        scan,
        eob,
        transform_size,
        tx_class,
        coeff_contexts); // NM - Assembly version is available in AOM

    for (c = eob - 1; c >= 0; --c) {

        const int32_t pos = scan[c];
        const tran_low_t v = qcoeff[pos];
        const int32_t is_nz = (v != 0);
        const int32_t level = abs(v);
        const int32_t coeff_ctx = coeff_contexts[pos];

        if (c == eob - 1) {
            ASSERT((AOMMIN(level, 3) - 1) >= 0);
            cost += coeff_costs->base_eob_cost[coeff_ctx][AOMMIN(level, 3) - 1];
        }
        else {
            cost += coeff_costs->base_cost[coeff_ctx][AOMMIN(level, 3)];
        }
        if (is_nz) {
            const int32_t sign = (v < 0) ? 1 : 0;
            // sign bit cost
            if (c == 0) {
                cost += coeff_costs->dc_sign_cost[dc_sign_ctx][sign];
            }
            else {
                cost += av1_cost_literal(1);
            }
            if (level > NUM_BASE_LEVELS) {
                int32_t ctx;
                ctx = get_br_ctx(levels, pos, bwl, transform_type);

                const int32_t base_range = level - 1 - NUM_BASE_LEVELS;
                if (base_range < COEFF_BASE_RANGE) {
                    cost += coeff_costs->lps_cost[ctx][base_range];
                }
                else {
                    cost += coeff_costs->lps_cost[ctx][COEFF_BASE_RANGE];
                }

                if (level >= 1 + NUM_BASE_LEVELS + COEFF_BASE_RANGE) {
                    cost += get_golomb_cost(level);
                }
            }
        }
    }
    return cost;
}


//static INLINE int32_t av1_get_skip_mode_context(const MacroBlockD *xd) {
//    const MbModeInfo *const above_mi = xd->above_mbmi;
//    const MbModeInfo *const left_mi = xd->left_mbmi;
//    const int32_t above_skip_mode = above_mi ? above_mi->skip_mode : 0;
//    const int32_t left_skip_mode = left_mi ? left_mi->skip_mode : 0;
//    return above_skip_mode + left_skip_mode;
//}
//
//static INLINE int32_t av1_get_skip_context(const MacroBlockD *xd) {
//    const MbModeInfo *const above_mi = xd->above_mbmi;
//    const MbModeInfo *const left_mi = xd->left_mbmi;
//    const int32_t above_skip = above_mi ? above_mi->skip : 0;
//    const int32_t left_skip = left_mi ? left_mi->skip : 0;
//    return above_skip + left_skip;
//}
/*********************************************************************************
* Av1IntraFastCost function is used to estimate the cost of an intra candidate mode
* for fast mode decisoion module in Intra or inter frame.
* Chroma cost is excluded from fast cost functions. Only the fast_chroma_rate is stored
* for future use in full loop
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
EbErrorType Av1IntraFastCost(
    struct ModeDecisionContext_s            *context_ptr,
    CodingUnit_t                            *cu_ptr,
    struct ModeDecisionCandidateBuffer_s    *candidate_buffer_ptr,
    uint32_t                                  qp,
    uint64_t                                  luma_distortion,
    uint64_t                                  chroma_distortion,
    uint64_t                                  lambda,
    PictureControlSet_t                     *picture_control_set_ptr)
{

    EbErrorType return_error = EB_ErrorNone;

    (void)qp;
    (void)picture_control_set_ptr;

    EbBool isMonochromeFlag = EB_FALSE; // NM - isMonochromeFlag is harcoded to false.

    EbBool isCflAllowed = (context_ptr->blk_geom->bwidth <= 32 &&
        context_ptr->blk_geom->bheight <= 32) ? 1 : 0;


    uint8_t   subSamplingX = 1; // NM - subsampling_x is harcoded to 1 for 420 chroma sampling.
    uint8_t   subSamplingY = 1; // NM - subsampling_y is harcoded to 1 for 420 chroma sampling.

    BlockSize cuSizeIndex = context_ptr->blk_geom->bsize;

    uint32_t miRow = context_ptr->cu_origin_y >> MI_SIZE_LOG2;
    uint32_t miCol = context_ptr->cu_origin_x >> MI_SIZE_LOG2;

    // In fast loop CFL alphas are not know yet. The chroma mode bits are calculated based on DC Mode, and if CFL is the winner compared to CFL, ChromaBits are updated
    uint32_t chroma_mode = candidate_buffer_ptr->candidate_ptr->intra_chroma_mode == UV_CFL_PRED ? UV_DC_PRED : candidate_buffer_ptr->candidate_ptr->intra_chroma_mode;

    // Number of bits for each synatax element
    uint64_t intraModeBitsNum = 0;
    uint64_t intraLumaModeBitsNum = 0;
    uint64_t intraLumaAngModeBitsNum = 0;
    uint64_t intraChromaModeBitsNum = 0;
    uint64_t intraChromaAngModeBitsNum = 0;
    uint64_t skipModeRate = 0;
    uint8_t  skipModeCtx = cu_ptr->skip_flag_context; // NM - Harcoded to 1 until the skip_mode context is added.

    PredictionMode intra_mode = (PredictionMode)candidate_buffer_ptr->candidate_ptr->pred_mode;

    // Luma and chroma rate
    uint64_t rate;
    uint64_t lumaRate = 0;
    uint64_t chromaRate = 0;
    uint64_t lumaSad, chromaSad;

    // Luma and chroma distortion
    uint64_t totalDistortion;
    uint32_t left_neighbor_mode = context_ptr->intra_luma_left_mode;
    uint32_t top_neighbor_mode = context_ptr->intra_luma_top_mode;

    const int32_t AboveCtx = intra_mode_context[top_neighbor_mode];
    const int32_t LeftCtx = intra_mode_context[left_neighbor_mode];

    intraModeBitsNum = picture_control_set_ptr->slice_type != I_SLICE ? (uint64_t)candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->mbModeFacBits[size_group_lookup[cuSizeIndex]][intra_mode] : ZERO_COST;

    skipModeRate = picture_control_set_ptr->slice_type != I_SLICE ? (uint64_t)candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->skipModeFacBits[skipModeCtx][0] : ZERO_COST;

    // Estimate luma nominal intra mode bits
    intraLumaModeBitsNum = picture_control_set_ptr->slice_type == I_SLICE ? (uint64_t)candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->yModeFacBits[AboveCtx][LeftCtx][intra_mode] : ZERO_COST;

    // Estimate luma angular mode bits
    if (candidate_buffer_ptr->candidate_ptr->is_directional_mode_flag && candidate_buffer_ptr->candidate_ptr->use_angle_delta) {
        ASSERT((intra_mode - V_PRED) < 8);
        ASSERT((intra_mode - V_PRED) >= 0);
        intraLumaAngModeBitsNum = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->angleDeltaFacBits[intra_mode - V_PRED][MAX_ANGLE_DELTA + candidate_buffer_ptr->candidate_ptr->angle_delta[PLANE_TYPE_Y]];
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

    if (context_ptr->blk_geom->has_uv) {

        if (!isMonochromeFlag && is_chroma_reference(miRow, miCol, cuSizeIndex, subSamplingX, subSamplingY)) {

            // Estimate luma nominal intra mode bits
            intraChromaModeBitsNum = (uint64_t)candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->intraUVmodeFacBits[isCflAllowed][intra_mode][chroma_mode];

            // Estimate luma angular mode bits
            if (candidate_buffer_ptr->candidate_ptr->is_directional_chroma_mode_flag && candidate_buffer_ptr->candidate_ptr->use_angle_delta) {
                intraChromaAngModeBitsNum = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->angleDeltaFacBits[chroma_mode - V_PRED][MAX_ANGLE_DELTA + candidate_buffer_ptr->candidate_ptr->angle_delta[PLANE_TYPE_UV]];


            }
        }
    }


    uint32_t isInterRate = picture_control_set_ptr->slice_type != I_SLICE ? candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->intraInterFacBits[cu_ptr->is_inter_ctx][0] : 0;
    lumaRate = intraModeBitsNum + skipModeRate + intraLumaModeBitsNum + intraLumaAngModeBitsNum + isInterRate;

    chromaRate = intraChromaModeBitsNum + intraChromaAngModeBitsNum;

    // Keep the Fast Luma and Chroma rate for future use
    candidate_buffer_ptr->candidate_ptr->fast_luma_rate = lumaRate;
    candidate_buffer_ptr->candidate_ptr->fast_chroma_rate = chromaRate;

    lumaSad = (LUMA_WEIGHT * luma_distortion) << AV1_COST_PRECISION;
    chromaSad = chroma_distortion << AV1_COST_PRECISION;
    totalDistortion = lumaSad + chromaSad;

    rate = lumaRate + chromaRate;

    // Assign fast cost
    *(candidate_buffer_ptr->fast_cost_ptr) = RDCOST(lambda, rate, totalDistortion);

    return return_error;
}

//extern INLINE int32_t have_newmv_in_inter_mode(PredictionMode mode);
static INLINE int32_t have_newmv_in_inter_mode(PredictionMode mode) {
    return (mode == NEWMV || mode == NEW_NEWMV || mode == NEAREST_NEWMV ||
        mode == NEW_NEARESTMV || mode == NEAR_NEWMV || mode == NEW_NEARMV);
}

extern void av1_set_ref_frame(MvReferenceFrame *rf,
    int8_t ref_frame_type);

// This function encodes the reference frame
uint64_t EstimateRefFramesNumBits(
    PictureControlSet_t                    *picture_control_set_ptr,
    ModeDecisionCandidateBuffer_t          *candidate_buffer_ptr,
    CodingUnit_t                           *cu_ptr,
    uint32_t                                 bwidth,
    uint32_t                                 bheight,
    uint8_t                                  ref_frame_type,
    EbBool                                is_compound)
{

    uint64_t refRateBits = 0;
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
        if (picture_control_set_ptr->parent_pcs_ptr->reference_mode == REFERENCE_MODE_SELECT) {
            if (MIN(bwidth, bheight) >= 8) {
                int32_t context = 0;
                context = cu_ptr->reference_mode_context;
                assert(context >= 0 && context < 5);
                refRateA = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->compInterFacBits[context][is_compound];

            }
        }
        else {
            assert((!is_compound) == (picture_control_set_ptr->parent_pcs_ptr->reference_mode == SINGLE_REFERENCE));
        }
        int32_t context = 0;
        if (is_compound) {
            const COMP_REFERENCE_TYPE comp_ref_type = /*has_uni_comp_refs(mbmi)
                                                      ? UNIDIR_COMP_REFERENCE
                                                      : */BIDIR_COMP_REFERENCE;
            MvReferenceFrame refType[2];
            av1_set_ref_frame(refType, ref_frame_type);


            context = cu_ptr->compoud_reference_type_context;
            assert(context >= 0 && context < 5);
            refRateB = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->compRefTypeFacBits[context][comp_ref_type];


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
            refRateC = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->compRefFacBits[context][0][bit];

            //            WRITE_REF_BIT(bit, comp_ref_p);

            if (!bit) {
                const int32_t bit1 = (refType[0] == LAST2_FRAME);
                context = av1_get_pred_context_comp_ref_p1(cu_ptr->av1xd);
                /*aom_write_symbol(ecWriter, bit1, frameContext->comp_ref_cdf[context][1],
                    2);*/
                assert(context >= 0 && context < 3);
                refRateD = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->compRefFacBits[context][1][bit1];


                //WRITE_REF_BIT(bit1, comp_ref_p1);
            }
            else {
                const int32_t bit2 = (refType[0] == GOLDEN_FRAME);
                context = av1_get_pred_context_comp_ref_p2(cu_ptr->av1xd);
                /*aom_write_symbol(ecWriter, bit2, frameContext->comp_ref_cdf[context][2],
                    2);*/
                assert(context >= 0 && context < 3);
                refRateE = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->compRefFacBits[context][2][bit2];

                //WRITE_REF_BIT(bit2, comp_ref_p2);
            }

            const int32_t bit_bwd = (refType[1] == ALTREF_FRAME);
            context = av1_get_pred_context_comp_bwdref_p(cu_ptr->av1xd);
            /*aom_write_symbol(ecWriter, bit_bwd, frameContext->comp_bwdref_cdf[context][0],
                2);*/
            assert(context >= 0 && context < 3);
            refRateF = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->compBwdRefFacBits[context][0][bit_bwd];
            //WRITE_REF_BIT(bit_bwd, comp_bwdref_p);


            if (!bit_bwd) {
                context = av1_get_pred_context_comp_bwdref_p1(cu_ptr->av1xd);
                /*aom_write_symbol(ecWriter, refType[1] == ALTREF2_FRAME, frameContext->comp_bwdref_cdf[context][1],
                    2);*/
                assert(context >= 0 && context < 3);
                refRateG = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->compBwdRefFacBits[context][1][refType[1] == ALTREF2_FRAME];
                //WRITE_REF_BIT(mbmi->ref_frame[1] == ALTREF2_FRAME, comp_bwdref_p1);

            }

        }
        else {
            const int32_t bit0 = (ref_frame_type <= ALTREF_FRAME &&
                ref_frame_type >= BWDREF_FRAME);//0

            context = av1_get_pred_context_single_ref_p1(cu_ptr->av1xd);
            /*aom_write_symbol(ecWriter, bit0, frameContext->single_ref_cdf[context][0],
                2);*/
            assert(context >= 0 && context < 3);
            refRateH = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->singleRefFacBits[context][0][bit0];
            //WRITE_REF_BIT(bit0, single_ref_p1);


            if (bit0) {
                const int32_t bit1 = (ref_frame_type == ALTREF_FRAME);
                context = av1_get_pred_context_single_ref_p2(cu_ptr->av1xd);
                assert(context >= 0 && context < 3);
                /*aom_write_symbol(ecWriter, bit1, frameContext->single_ref_cdf[context][1],
                    2);*/
                refRateI = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->singleRefFacBits[context][1][bit1];
                //WRITE_REF_BIT(bit1, single_ref_p2);


                if (!bit1) {
                    context = av1_get_pred_context_single_ref_p6(cu_ptr->av1xd);
                    /*aom_write_symbol(ecWriter, cu_ptr->prediction_unit_array[0].ref_frame_type == ALTREF2_FRAME, frameContext->single_ref_cdf[context][5],
                        2);*/
                    assert(context >= 0 && context < 3);
                    refRateJ = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->singleRefFacBits[context][5][ref_frame_type == ALTREF2_FRAME];
                    //WRITE_REF_BIT(mbmi->ref_frame[0] == ALTREF2_FRAME, single_ref_p6);

                }
            }
            else {
                const int32_t bit2 = (ref_frame_type == LAST3_FRAME ||
                    ref_frame_type == GOLDEN_FRAME); //0
                context = av1_get_pred_context_single_ref_p3(cu_ptr->av1xd);
                /*aom_write_symbol(ecWriter, bit2, frameContext->single_ref_cdf[context][2],
                    2);*/
                assert(context >= 0 && context < 3);
                refRateK = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->singleRefFacBits[context][2][bit2];
                //WRITE_REF_BIT(bit2, single_ref_p3);

                if (!bit2) {
                    const int32_t bit3 = (ref_frame_type != LAST_FRAME); //0;
                    context = av1_get_pred_context_single_ref_p4(cu_ptr->av1xd);
                    assert(context >= 0 && context < 3);
                    /*aom_write_symbol(ecWriter, bit3, frameContext->single_ref_cdf[context][3],
                        2);*/
                    refRateL = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->singleRefFacBits[context][3][bit3];
                    //WRITE_REF_BIT(bit3, single_ref_p4);

                }
                else {
                    const int32_t bit4 = (ref_frame_type != LAST3_FRAME);
                    context = av1_get_pred_context_single_ref_p5(cu_ptr->av1xd);
                    /*aom_write_symbol(ecWriter, bit4, frameContext->single_ref_cdf[context][4],
                        2);*/
                    assert(context >= 0 && context < 3);
                    refRateM = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->singleRefFacBits[context][4][bit4];
                    //WRITE_REF_BIT(bit4, single_ref_p5);

                }
            }
        }
    }

    refRateBits = refRateA + refRateB + refRateC + refRateD + refRateE + refRateF + refRateG + refRateH + refRateI + refRateJ + refRateK + refRateL + refRateM;
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
    ASSERT((refmv_ctx >> 1) < 3);
    const int16_t comp_ctx = compound_mode_ctx_map_2[refmv_ctx >> 1][AOMMIN(
        newmv_ctx, COMP_NEWMV_CTXS - 1)];
    return comp_ctx;
}

/*********************************************************************************
* Av1InterFastCost function is used to estimate the cost of an inter candidate mode
* for fast mode decisoion module in Inter frame.
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
EbErrorType Av1InterFastCost(
    struct ModeDecisionContext_s           *context_ptr,
    CodingUnit_t                           *cu_ptr,
    ModeDecisionCandidateBuffer_t          *candidate_buffer_ptr,
    uint32_t                                  qp,
    uint64_t                                  luma_distortion,
    uint64_t                                  chroma_distortion,
    uint64_t                                  lambda,
    PictureControlSet_t                    *picture_control_set_ptr)
{
    EbErrorType  return_error = EB_ErrorNone;
    ModeDecisionCandidate_t *candidate_ptr = candidate_buffer_ptr->candidate_ptr;
    // Luma rate
    uint64_t           lumaRate = 0;
    uint64_t           chromaRate = 0;
    uint64_t           mvRate = 0;
    uint64_t           skipModeRate;
    // Luma and chroma distortion
    uint64_t           lumaSad;
    uint64_t             chromaSad;
    uint64_t           totalDistortion;

    uint64_t           rate;

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
    skipModeRate = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->skipModeFacBits[skipModeCtx][0];

    uint64_t referencePictureBitsNum = 0;

    //Reference Type and Mode Bit estimation

    referencePictureBitsNum = EstimateRefFramesNumBits(
        picture_control_set_ptr,
        candidate_buffer_ptr,
        cu_ptr,
        context_ptr->blk_geom->bwidth,
        context_ptr->blk_geom->bheight,
        candidate_ptr->ref_frame_type,
        candidate_ptr->is_compound);


    if (candidate_ptr->is_compound) {
        interModeBitsNum += candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->interCompoundModeFacBits[modeCtx][INTER_COMPOUND_OFFSET(inter_mode)];
    }
    else {
        //uint32_t newmv_ctx = modeCtx & NEWMV_CTX_MASK;
        //interModeBitsNum = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->newMvModeFacBits[mode_ctx][0];

        int16_t newmv_ctx = modeCtx & NEWMV_CTX_MASK;
        //aom_write_symbol(ecWriter, mode != NEWMV, frameContext->newmv_cdf[newmv_ctx], 2);
        interModeBitsNum += candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->newMvModeFacBits[newmv_ctx][inter_mode != NEWMV];

        if (inter_mode != NEWMV) {
            const int16_t zeromvCtx = (modeCtx >> GLOBALMV_OFFSET) & GLOBALMV_CTX_MASK;
            //aom_write_symbol(ecWriter, mode != GLOBALMV, frameContext->zeromv_cdf[zeromvCtx], 2);
            interModeBitsNum += candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->zeroMvModeFacBits[zeromvCtx][inter_mode != GLOBALMV];

            if (inter_mode != GLOBALMV) {
                int16_t refmvCtx = (modeCtx >> REFMV_OFFSET) & REFMV_CTX_MASK;
                /*aom_write_symbol(ecWriter, mode != NEARESTMV, frameContext->refmv_cdf[refmv_ctx], 2);*/
                interModeBitsNum += candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->refMvModeFacBits[refmvCtx][inter_mode != NEARESTMV];
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
                        av1_drl_ctx(context_ptr->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack[candidate_ptr->ref_frame_type], idx);

                    interModeBitsNum += candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->drlModeFacBits[drl1Ctx][candidate_ptr->drl_index != idx];
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
                        av1_drl_ctx(context_ptr->md_local_cu_unit[context_ptr->blk_geom->blkidx_mds].ed_ref_mv_stack[candidate_ptr->ref_frame_type], idx);

                    interModeBitsNum += candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->drlModeFacBits[drl_ctx][candidate_ptr->drl_index != (idx - 1)];


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
                    mvRefX = refListIdx == REF_LIST_1 ? candidate_ptr->motionVector_x_L1 : candidate_ptr->motionVector_x_L0;
                    mvRefY = refListIdx == REF_LIST_1 ? candidate_ptr->motionVector_y_L1 : candidate_ptr->motionVector_y_L0;


                    MV mv;
                    mv.row = mvRefY;
                    mv.col = mvRefX;

                    MV ref_mv;
                    ref_mv.row = predRefY;
                    ref_mv.col = predRefX;

                    mvRate += av1_mv_bit_cost(
                        &mv,
                        &ref_mv,
                        candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->nmv_vec_cost,
                        candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->nmvcoststack,
                        MV_COST_WEIGHT);
                }

            }
            else if (inter_mode == NEAREST_NEWMV || inter_mode == NEAR_NEWMV) {

                predRefX = candidate_ptr->motion_vector_pred_x[REF_LIST_1];
                predRefY = candidate_ptr->motion_vector_pred_y[REF_LIST_1];
                mvRefX = candidate_ptr->motionVector_x_L1;
                mvRefY = candidate_ptr->motionVector_y_L1;


                MV mv;
                mv.row = mvRefY;
                mv.col = mvRefX;

                MV ref_mv;
                ref_mv.row = predRefY;
                ref_mv.col = predRefX;

                mvRate += av1_mv_bit_cost(
                    &mv,
                    &ref_mv,
                    candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->nmv_vec_cost,
                    candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->nmvcoststack,
                    MV_COST_WEIGHT);

            }
            else {
                assert(inter_mode == NEW_NEARESTMV || inter_mode == NEW_NEARMV);

                predRefX = candidate_ptr->motion_vector_pred_x[REF_LIST_0];
                predRefY = candidate_ptr->motion_vector_pred_y[REF_LIST_0];
                mvRefX = candidate_ptr->motionVector_x_L0;
                mvRefY = candidate_ptr->motionVector_y_L0;

                MV mv;
                mv.row = mvRefY;
                mv.col = mvRefX;

                MV ref_mv;
                ref_mv.row = predRefY;
                ref_mv.col = predRefX;

                mvRate += av1_mv_bit_cost(
                    &mv,
                    &ref_mv,
                    candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->nmv_vec_cost,
                    candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->nmvcoststack,
                    MV_COST_WEIGHT);


            }
        }
        else {

            refListIdx = candidate_ptr->prediction_direction[0] == 0 ? 0 : 1;

            predRefX = candidate_ptr->motion_vector_pred_x[refListIdx];
            predRefY = candidate_ptr->motion_vector_pred_y[refListIdx];

            mvRefX = refListIdx == 0 ? candidate_ptr->motionVector_x_L0 : candidate_ptr->motionVector_x_L1;
            mvRefY = refListIdx == 0 ? candidate_ptr->motionVector_y_L0 : candidate_ptr->motionVector_y_L1;

            MV mv;
            mv.row = mvRefY;
            mv.col = mvRefX;

            MV ref_mv;
            ref_mv.row = predRefY;
            ref_mv.col = predRefX;

            mvRate = av1_mv_bit_cost(
                &mv,
                &ref_mv,
                candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->nmv_vec_cost,
                candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->nmvcoststack,
                MV_COST_WEIGHT);
        }
    }

    // NM - To be added when the intrainter mode is adopted
    //  read_interintra_mode(is_compound)

    EbBool is_inter = inter_mode >= SINGLE_INTER_MODE_START && inter_mode < SINGLE_INTER_MODE_END;
    if (is_inter
        && picture_control_set_ptr->parent_pcs_ptr->switchable_motion_mode
        && rf[1] != INTRA_FRAME)
    {
        MOTION_MODE motion_mode_rd = candidate_buffer_ptr->candidate_ptr->motion_mode;
        BlockSize bsize = context_ptr->blk_geom->bsize;

        cu_ptr->prediction_unit_array[0].overlappable_neighbors[0] = 1;
        cu_ptr->prediction_unit_array[0].overlappable_neighbors[1] = 1;

        MOTION_MODE motion_allowed =  motion_mode_allowed(
            picture_control_set_ptr,
            cu_ptr,
            bsize,
            rf[0],
            rf[1],
            inter_mode);

        switch (motion_allowed) {
        case SIMPLE_TRANSLATION: break;
        case OBMC_CAUSAL:
            interModeBitsNum += candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->motionModeFacBits1[bsize][SIMPLE_TRANSLATION]; // TODO: modify when OBMC added
            break;
        default:
            interModeBitsNum += candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->motionModeFacBits[bsize][motion_mode_rd];
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

    uint32_t isInterRate = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->intraInterFacBits[cu_ptr->is_inter_ctx][1];
    lumaRate = referencePictureBitsNum + skipModeRate + interModeBitsNum + mvRate + isInterRate;


    //chromaRate = intraChromaModeBitsNum + intraChromaAngModeBitsNum;

    // Keep the Fast Luma and Chroma rate for future use
    candidate_buffer_ptr->candidate_ptr->fast_luma_rate = lumaRate;
    candidate_buffer_ptr->candidate_ptr->fast_chroma_rate = chromaRate;


    lumaSad = (LUMA_WEIGHT * luma_distortion) << AV1_COST_PRECISION;
    chromaSad = chroma_distortion << AV1_COST_PRECISION;
    totalDistortion = lumaSad + chromaSad;

    if (context_ptr->blk_geom->has_uv == 0 && chromaSad != 0) {
        printf("Av1InterFastCost: Chroma error");
    }


    rate = lumaRate + chromaRate;


    // Assign fast cost
    *(candidate_buffer_ptr->fast_cost_ptr) = RDCOST(lambda, rate, totalDistortion);

    if (candidate_buffer_ptr->candidate_ptr->merge_flag) {
        uint64_t skipModeRate = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->skipModeFacBits[skipModeCtx][1];
        if (skipModeRate < rate) {

            *(candidate_buffer_ptr->fast_cost_ptr) = RDCOST(lambda, skipModeRate, totalDistortion);
        }

    }


    return return_error;
}


EbErrorType Av1TuEstimateCoeffBits(
    PictureControlSet_t                    *picture_control_set_ptr,
    struct ModeDecisionCandidateBuffer_s   *candidate_buffer_ptr,
    CodingUnit_t                           *cu_ptr,
    uint32_t                                  tuOriginIndex,
    uint32_t                                  tuChromaOriginIndex,
    EntropyCoder_t                         *entropy_coder_ptr,
    EbPictureBufferDesc_t                  *coeff_buffer_sb,
    uint32_t                                 yEob,
    uint32_t                                 cbEob,
    uint32_t                                 crEob,
    uint64_t                                 *yTuCoeffBits,
    uint64_t                                 *cbTuCoeffBits,
    uint64_t                                 *crTuCoeffBits,
    TxSize                                 txsize,
    TxSize                                 txsize_uv,
    COMPONENT_TYPE                          componentType,
    EbAsm                                  asm_type)
{
    (void)asm_type;
    (void)entropy_coder_ptr;
    EbErrorType return_error = EB_ErrorNone;


    int32_t *coeffBuffer;


    int16_t  luma_txb_skip_context = cu_ptr->luma_txb_skip_context;
    int16_t  luma_dc_sign_context = cu_ptr->luma_dc_sign_context;
    int16_t  cb_txb_skip_context = cu_ptr->cb_txb_skip_context;
    int16_t  cb_dc_sign_context = cu_ptr->cb_dc_sign_context;
    int16_t  cr_txb_skip_context = cu_ptr->cr_txb_skip_context;
    int16_t  cr_dc_sign_context = cu_ptr->cr_dc_sign_context;


    EbBool reducedTransformSetFlag = picture_control_set_ptr->parent_pcs_ptr->reduced_tx_set_used ? EB_TRUE : EB_FALSE;

    //Estimate the rate of the transform type and coefficient for Luma

    if (componentType == COMPONENT_LUMA || componentType == COMPONENT_ALL) {
        if (yEob) {
            coeffBuffer = (int32_t*)&coeff_buffer_sb->bufferY[tuOriginIndex * sizeof(int32_t)];

            *yTuCoeffBits = av1_cost_coeffs_txb(
                candidate_buffer_ptr,
                coeffBuffer,
                (uint16_t)yEob,
                PLANE_TYPE_Y,
                txsize,
                luma_txb_skip_context,
                luma_dc_sign_context,
                reducedTransformSetFlag);
        }
        else {
            *yTuCoeffBits = av1_cost_skip_txb(
                candidate_buffer_ptr,
                txsize,
                PLANE_TYPE_Y,
                luma_txb_skip_context);
        }
    }
    //Estimate the rate of the transform type and coefficient for chroma Cb

    if (componentType == COMPONENT_CHROMA_CB || componentType == COMPONENT_CHROMA || componentType == COMPONENT_ALL) {

        if (cbEob) {

            coeffBuffer = (int32_t*)&coeff_buffer_sb->bufferCb[tuChromaOriginIndex * sizeof(int32_t)];


            *cbTuCoeffBits = av1_cost_coeffs_txb(
                candidate_buffer_ptr,
                coeffBuffer,
                (uint16_t)cbEob,
                PLANE_TYPE_UV,
                txsize_uv,
                cb_txb_skip_context,
                cb_dc_sign_context,
                reducedTransformSetFlag);

        }
        else {
            *cbTuCoeffBits = av1_cost_skip_txb(
                candidate_buffer_ptr,
                txsize_uv,
                PLANE_TYPE_UV,
                cb_txb_skip_context);
        }
    }

    if (componentType == COMPONENT_CHROMA_CR || componentType == COMPONENT_CHROMA || componentType == COMPONENT_ALL) {

        //Estimate the rate of the transform type and coefficient for chroma Cr
        if (crEob) {

            coeffBuffer = (int32_t*)&coeff_buffer_sb->bufferCr[tuChromaOriginIndex * sizeof(int32_t)];

            *crTuCoeffBits = av1_cost_coeffs_txb(
                candidate_buffer_ptr,
                coeffBuffer,
                (uint16_t)crEob,
                PLANE_TYPE_UV,
                txsize_uv,
                cr_txb_skip_context,
                cr_dc_sign_context,
                reducedTransformSetFlag);

        }
        else {
            *crTuCoeffBits = av1_cost_skip_txb(
                candidate_buffer_ptr,
                txsize_uv,
                PLANE_TYPE_UV,
                cr_txb_skip_context);
        }
    }


    return return_error;
}
/*********************************************************************************
* Av1IntraFullCost function is used to estimate the cost of an intra candidate mode
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
    PictureControlSet_t                    *picture_control_set_ptr,
    ModeDecisionContext_t                  *context_ptr,
    struct ModeDecisionCandidateBuffer_s   *candidate_buffer_ptr,
    CodingUnit_t                           *cu_ptr,
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

            chromaRate += candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->cflAlphaFacBits[candidate_buffer_ptr->candidate_ptr->cfl_alpha_signs][CFL_PRED_U][CFL_IDX_U(candidate_buffer_ptr->candidate_ptr->cfl_alpha_idx)] +
                candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->cflAlphaFacBits[candidate_buffer_ptr->candidate_ptr->cfl_alpha_signs][CFL_PRED_V][CFL_IDX_V(candidate_buffer_ptr->candidate_ptr->cfl_alpha_idx)];

            chromaRate += (uint64_t)candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->intraUVmodeFacBits[isCflAllowed][candidate_buffer_ptr->candidate_ptr->intra_luma_mode][UV_CFL_PRED];
            chromaRate -= (uint64_t)candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->intraUVmodeFacBits[isCflAllowed][candidate_buffer_ptr->candidate_ptr->intra_luma_mode][UV_DC_PRED];
        }
    }

    // Coeff rate
    coeffRate = (*y_coeff_bits + *cb_coeff_bits + *cr_coeff_bits);
    luma_sse = y_distortion[0];
    chromaSse = cb_distortion[0] + cr_distortion[0];

    // *Note - As of Oct 2011, the JCT-VC uses the PSNR forumula
    //  PSNR = (LUMA_WEIGHT * PSNRy + PSNRu + PSNRv) / (2+LUMA_WEIGHT)
    luma_sse = LUMA_WEIGHT * (luma_sse << AV1_COST_PRECISION);

    // *Note - As in JCTVC-G1102, the JCT-VC uses the Mode Decision forumula where the chromaSse has been weighted
    //  CostMode = (luma_sse + wchroma * chromaSse) + lambdaSse * rateMode
    //chromaSse = (((chromaSse * ChromaWeightFactorLd[qp]) + CHROMA_WEIGHT_OFFSET) >> CHROMA_WEIGHT_SHIFT); // Low delay and Random access have the same value of chroma weight

    chromaSse = (chromaSse << AV1_COST_PRECISION);

    totalDistortion = luma_sse + chromaSse;

    rate = lumaRate + chromaRate + coeffRate;
    // Assign full cost
    *(candidate_buffer_ptr->full_cost_ptr) = RDCOST(lambda, rate, totalDistortion);

    candidate_buffer_ptr->full_lambda_rate = *candidate_buffer_ptr->full_cost_ptr - totalDistortion;
    coeffRate = *y_coeff_bits;
    candidate_buffer_ptr->full_cost_luma = RDCOST(lambda, lumaRate + *y_coeff_bits, luma_sse);

    return return_error;
}

/*********************************************************************************
* MergeSkipFullCost function is used to estimate the cost of an AMVPSkip candidate
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
    PictureControlSet_t                    *picture_control_set_ptr,
    ModeDecisionContext_t                  *context_ptr,
    struct ModeDecisionCandidateBuffer_s   *candidate_buffer_ptr,
    CodingUnit_t                           *cu_ptr,
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

    uint64_t skipModeRate = candidate_buffer_ptr->candidate_ptr->md_rate_estimation_ptr->skipModeFacBits[skipModeCtx][1];

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
    //        mergeChromaSse = (((mergeChromaSse * ChromaWeightFactorRa[qp]) + CHROMA_WEIGHT_OFFSET) >> CHROMA_WEIGHT_SHIFT);
    //    }
    //    else if (picture_control_set_ptr->temporal_layer_index < 3) {
    //        mergeChromaSse = (((mergeChromaSse * ChromaWeightFactorRaQpScalingL1[qp]) + CHROMA_WEIGHT_OFFSET) >> CHROMA_WEIGHT_SHIFT);
    //    }
    //    else {
    //        mergeChromaSse = (((mergeChromaSse * ChromaWeightFactorRaQpScalingL3[qp]) + CHROMA_WEIGHT_OFFSET) >> CHROMA_WEIGHT_SHIFT);
    //    }
    //}
    //else {
    //    // Low delay
    //    if (picture_control_set_ptr->temporal_layer_index == 0) {
    //        mergeChromaSse = (((mergeChromaSse * ChromaWeightFactorLd[qp]) + CHROMA_WEIGHT_OFFSET) >> CHROMA_WEIGHT_SHIFT);
    //    }
    //    else {
    //        mergeChromaSse = (((mergeChromaSse * ChromaWeightFactorLdQpScaling[qp]) + CHROMA_WEIGHT_OFFSET) >> CHROMA_WEIGHT_SHIFT);
    //    }
    //}

    // Add fast rate to get the total rate of the subject mode
    mergeRate += candidate_buffer_ptr->candidate_ptr->fast_luma_rate;
    mergeRate += candidate_buffer_ptr->candidate_ptr->fast_chroma_rate;


    mergeRate += coeffRate;

    mergeDistortion = (mergeLumaSse + mergeChromaSse);

    //merge_cost = mergeDistortion + (((lambda * coeffRate + lambda * mergeLumaRate + lambda_chroma * mergeChromaRate) + MD_OFFSET) >> MD_SHIFT);

    merge_cost = RDCOST(lambda, mergeRate, mergeDistortion);
    // mergeLumaCost = mergeLumaSse    + (((lambda * lumaCoeffRate + lambda * mergeLumaRate) + MD_OFFSET) >> MD_SHIFT);


    // *Note - As in JCTVC-G1102, the JCT-VC uses the Mode Decision forumula where the chromaSse has been weighted
    //  CostMode = (luma_sse + wchroma * chromaSse) + lambdaSse * rateMode

    //if (picture_control_set_ptr->parent_pcs_ptr->pred_structure == EB_PRED_RANDOM_ACCESS) {

    //    if (picture_control_set_ptr->temporal_layer_index == 0) {
    //        skipChromaSse = (((skipChromaSse * ChromaWeightFactorRa[qp]) + CHROMA_WEIGHT_OFFSET) >> CHROMA_WEIGHT_SHIFT);
    //    }
    //    else if (picture_control_set_ptr->temporal_layer_index < 3) {
    //        skipChromaSse = (((skipChromaSse * ChromaWeightFactorRaQpScalingL1[qp]) + CHROMA_WEIGHT_OFFSET) >> CHROMA_WEIGHT_SHIFT);
    //    }
    //    else {
    //        skipChromaSse = (((skipChromaSse * ChromaWeightFactorRaQpScalingL3[qp]) + CHROMA_WEIGHT_OFFSET) >> CHROMA_WEIGHT_SHIFT);
    //    }
    //}
    //else {
    //    // Low Delay
    //    if (picture_control_set_ptr->temporal_layer_index == 0) {
    //        skipChromaSse = (((skipChromaSse * ChromaWeightFactorLd[qp]) + CHROMA_WEIGHT_OFFSET) >> CHROMA_WEIGHT_SHIFT);
    //    }
    //    else {
    //        skipChromaSse = (((skipChromaSse * ChromaWeightFactorLdQpScaling[qp]) + CHROMA_WEIGHT_OFFSET) >> CHROMA_WEIGHT_SHIFT);
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
* Av1IntraFullCost function is used to estimate the cost of an intra candidate mode
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
EbErrorType Av1IntraFullCost(
    PictureControlSet_t                    *picture_control_set_ptr,
    ModeDecisionContext_t                  *context_ptr,
    struct ModeDecisionCandidateBuffer_s   *candidate_buffer_ptr,
    CodingUnit_t                           *cu_ptr,
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
* Av1InterFullCost function is used to estimate the cost of an inter candidate mode
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
EbErrorType Av1InterFullCost(
    PictureControlSet_t                    *picture_control_set_ptr,
    ModeDecisionContext_t                  *context_ptr,
    struct ModeDecisionCandidateBuffer_s   *candidate_buffer_ptr,
    CodingUnit_t                           *cu_ptr,
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
void CodingLoopContextGeneration(
    ModeDecisionContext_t      *context_ptr,
    CodingUnit_t               *cu_ptr,
    uint32_t                      cu_origin_x,
    uint32_t                      cu_origin_y,
    uint32_t                      sb_sz,

    NeighborArrayUnit_t        *skip_coeff_neighbor_array,
    NeighborArrayUnit_t        *luma_dc_sign_level_coeff_neighbor_array,
    NeighborArrayUnit_t        *cb_dc_sign_level_coeff_neighbor_array,
    NeighborArrayUnit_t        *cr_dc_sign_level_coeff_neighbor_array,
    NeighborArrayUnit_t        *inter_pred_dir_neighbor_array,
    NeighborArrayUnit_t        *ref_frame_type_neighbor_array,

    NeighborArrayUnit_t        *intraLumaNeighborArray,
    NeighborArrayUnit_t        *skip_flag_neighbor_array,
    NeighborArrayUnit_t        *mode_type_neighbor_array,
    NeighborArrayUnit_t        *leaf_depth_neighbor_array,
    NeighborArrayUnit_t       *leaf_partition_neighbor_array)
{
    (void)sb_sz;
    uint32_t modeTypeLeftNeighborIndex = GetNeighborArrayUnitLeftIndex(
        mode_type_neighbor_array,
        cu_origin_y);
    uint32_t modeTypeTopNeighborIndex = GetNeighborArrayUnitTopIndex(
        mode_type_neighbor_array,
        cu_origin_x);
    uint32_t leafDepthLeftNeighborIndex = GetNeighborArrayUnitLeftIndex(
        leaf_depth_neighbor_array,
        cu_origin_y);
    uint32_t leafDepthTopNeighborIndex = GetNeighborArrayUnitTopIndex(
        leaf_depth_neighbor_array,
        cu_origin_x);
    uint32_t skipFlagLeftNeighborIndex = GetNeighborArrayUnitLeftIndex(
        skip_flag_neighbor_array,
        cu_origin_y);
    uint32_t skipFlagTopNeighborIndex = GetNeighborArrayUnitTopIndex(
        skip_flag_neighbor_array,
        cu_origin_x);
    uint32_t intraLumaModeLeftNeighborIndex = GetNeighborArrayUnitLeftIndex(
        intraLumaNeighborArray,
        cu_origin_y);
    uint32_t intraLumaModeTopNeighborIndex = GetNeighborArrayUnitTopIndex(
        intraLumaNeighborArray,
        cu_origin_x);

    uint32_t partition_left_neighbor_index = GetNeighborArrayUnitLeftIndex(
        leaf_partition_neighbor_array,
        cu_origin_y);
    uint32_t partition_above_neighbor_index = GetNeighborArrayUnitTopIndex(
        leaf_partition_neighbor_array,
        cu_origin_x);

    // Intra Luma Neighbor Modes

    cu_ptr->prediction_unit_array->intra_luma_left_mode = (uint32_t)(
        (mode_type_neighbor_array->leftArray[modeTypeLeftNeighborIndex] != INTRA_MODE) ? (uint32_t)DC_PRED :
        intraLumaNeighborArray->leftArray[intraLumaModeLeftNeighborIndex]);

    cu_ptr->prediction_unit_array->intra_luma_top_mode = (uint32_t)(
        (mode_type_neighbor_array->topArray[modeTypeTopNeighborIndex] != INTRA_MODE) ? (uint32_t)DC_PRED :
        intraLumaNeighborArray->topArray[intraLumaModeTopNeighborIndex]);

    int32_t contextIndex;
    if (mode_type_neighbor_array->leftArray[modeTypeLeftNeighborIndex] != (uint8_t)INVALID_MODE && mode_type_neighbor_array->topArray[modeTypeTopNeighborIndex] != (uint8_t)INVALID_MODE) {
        contextIndex = (mode_type_neighbor_array->leftArray[modeTypeLeftNeighborIndex] == (uint8_t)INTRA_MODE && mode_type_neighbor_array->topArray[modeTypeTopNeighborIndex] == (uint8_t)INTRA_MODE) ? 3 :
            (mode_type_neighbor_array->leftArray[modeTypeLeftNeighborIndex] == (uint8_t)INTRA_MODE || mode_type_neighbor_array->topArray[modeTypeTopNeighborIndex] == (uint8_t)INTRA_MODE) ? 1 : 0;

    }
    else  if (mode_type_neighbor_array->leftArray[modeTypeLeftNeighborIndex] != (uint8_t)INVALID_MODE) {
        contextIndex = (mode_type_neighbor_array->leftArray[modeTypeLeftNeighborIndex] == (uint8_t)INTRA_MODE) ? 2 : 0;
    }
    else if (mode_type_neighbor_array->topArray[modeTypeTopNeighborIndex] != (uint8_t)INVALID_MODE) {
        contextIndex = (mode_type_neighbor_array->topArray[modeTypeTopNeighborIndex] == (uint8_t)INTRA_MODE) ? 2 : 0;
    }
    else {
        contextIndex = 0;
    }

    cu_ptr->is_inter_ctx = contextIndex;
    //  if(cu_ptr->is_inter_ctx!=0) //
    //      printf("ctx:%i \n",cu_ptr->is_inter_ctx);

      //   Top Intra Mode Neighbor Array instead of a Full
      // Skip Flag Context
    cu_ptr->skip_flag_context =
        (mode_type_neighbor_array->leftArray[modeTypeLeftNeighborIndex] == (uint8_t)INVALID_MODE) ? 0 :
        (skip_flag_neighbor_array->leftArray[skipFlagLeftNeighborIndex] == EB_TRUE) ? 1 : 0;
    cu_ptr->skip_flag_context +=
        (mode_type_neighbor_array->topArray[modeTypeTopNeighborIndex] == (uint8_t)INVALID_MODE) ? 0 :
        (skip_flag_neighbor_array->topArray[skipFlagTopNeighborIndex] == EB_TRUE) ? 1 : 0;

    // Split Flag Context (neighbor info)
    context_ptr->md_local_cu_unit[cu_ptr->mds_idx].left_neighbor_mode = (uint32_t)(
        (mode_type_neighbor_array->leftArray[modeTypeLeftNeighborIndex] != INTRA_MODE) ? (uint32_t)DC_PRED :
        intraLumaNeighborArray->leftArray[intraLumaModeLeftNeighborIndex]);
    context_ptr->md_local_cu_unit[cu_ptr->mds_idx].left_neighbor_depth = leaf_depth_neighbor_array->leftArray[leafDepthLeftNeighborIndex];
    context_ptr->md_local_cu_unit[cu_ptr->mds_idx].top_neighbor_mode = (uint32_t)(
        (mode_type_neighbor_array->topArray[modeTypeTopNeighborIndex] != INTRA_MODE) ? (uint32_t)DC_PRED :
        intraLumaNeighborArray->topArray[intraLumaModeTopNeighborIndex]);
    context_ptr->md_local_cu_unit[cu_ptr->mds_idx].top_neighbor_depth = leaf_depth_neighbor_array->topArray[leafDepthTopNeighborIndex];


    // Generate Partition context
    context_ptr->md_local_cu_unit[cu_ptr->mds_idx].above_neighbor_partition = (((PartitionContext*)leaf_partition_neighbor_array->topArray)[partition_above_neighbor_index].above == (int8_t)INVALID_NEIGHBOR_DATA) ?
        0 : ((PartitionContext*)leaf_partition_neighbor_array->topArray)[partition_above_neighbor_index].above;

    context_ptr->md_local_cu_unit[cu_ptr->mds_idx].left_neighbor_partition = (((PartitionContext*)leaf_partition_neighbor_array->leftArray)[partition_left_neighbor_index].left == (int8_t)INVALID_NEIGHBOR_DATA) ?
        0 : ((PartitionContext*)leaf_partition_neighbor_array->leftArray)[partition_left_neighbor_index].left;

    // Skip Coeff AV1 Context
    uint32_t skipCoeffLeftNeighborIndex = GetNeighborArrayUnitLeftIndex(
        skip_coeff_neighbor_array,
        cu_origin_y);
    uint32_t skipCoeffTopNeighborIndex = GetNeighborArrayUnitTopIndex(
        skip_coeff_neighbor_array,
        cu_origin_x);


    cu_ptr->skip_coeff_context =
        (skip_coeff_neighbor_array->leftArray[skipCoeffLeftNeighborIndex] == (uint8_t)INVALID_NEIGHBOR_DATA) ? 0 :
        (skip_coeff_neighbor_array->leftArray[skipCoeffLeftNeighborIndex]) ? 1 : 0;


    cu_ptr->skip_coeff_context +=
        (skip_coeff_neighbor_array->topArray[skipCoeffTopNeighborIndex] == (uint8_t)INVALID_NEIGHBOR_DATA) ? 0 :
        (skip_coeff_neighbor_array->topArray[skipCoeffTopNeighborIndex]) ? 1 : 0;

    // Skip and Dc sign context generation

    BlockSize plane_bsize = context_ptr->blk_geom->bsize;

    cu_ptr->luma_txb_skip_context = 0;
    cu_ptr->luma_dc_sign_context = 0;
    cu_ptr->cb_txb_skip_context = 0;
    cu_ptr->cb_dc_sign_context = 0;
    cu_ptr->cr_txb_skip_context = 0;
    cu_ptr->cr_dc_sign_context = 0;

    int32_t txb_count = context_ptr->blk_geom->txb_count;
    int32_t txb_itr = 0;
    for (txb_itr = 0; txb_itr < txb_count; txb_itr++) {


        GetTxbCtx(                  //SB128_TODO move inside Full loop
            COMPONENT_LUMA,
            luma_dc_sign_level_coeff_neighbor_array,
            cu_origin_x,
            cu_origin_y,
            plane_bsize,
            context_ptr->blk_geom->txsize[txb_itr],
            &cu_ptr->luma_txb_skip_context,
            &cu_ptr->luma_dc_sign_context);


#if CHROMA_BLIND
        if (context_ptr->blk_geom->has_uv && context_ptr->chroma_level == CHROMA_MODE_0) {
#else
        if (context_ptr->blk_geom->has_uv) {
#endif
            GetTxbCtx(
                COMPONENT_CHROMA,
                cb_dc_sign_level_coeff_neighbor_array,
                context_ptr->round_origin_x >> 1,
                context_ptr->round_origin_y >> 1,
                context_ptr->blk_geom->bsize_uv,
                context_ptr->blk_geom->txsize_uv[txb_itr],
                &cu_ptr->cb_txb_skip_context,
                &cu_ptr->cb_dc_sign_context);
            GetTxbCtx(
                COMPONENT_CHROMA,
                cr_dc_sign_level_coeff_neighbor_array,
                context_ptr->round_origin_x >> 1,
                context_ptr->round_origin_y >> 1,
                context_ptr->blk_geom->bsize_uv,
                context_ptr->blk_geom->txsize_uv[txb_itr],
                &cu_ptr->cr_txb_skip_context,
                &cu_ptr->cr_dc_sign_context);
        }
    }

    // Generate reference mode context

    cu_ptr->reference_mode_context = (uint8_t)Av1GetReferenceModeContext(
        cu_origin_x,
        cu_origin_y,
        mode_type_neighbor_array,
        inter_pred_dir_neighbor_array);

    cu_ptr->compoud_reference_type_context = (uint8_t)Av1GetCompReferenceTypeContext(
        cu_origin_x,
        cu_origin_y,
        mode_type_neighbor_array,
        inter_pred_dir_neighbor_array);

    //Collect Neighbor ref cout
    Av1CollectNeighborsRefCounts(
        cu_ptr,
        cu_origin_x,
        cu_origin_y,
        mode_type_neighbor_array,
        inter_pred_dir_neighbor_array,
        ref_frame_type_neighbor_array);

    return;
}

/********************************************
* TuCalcCost
*   computes TU Cost and generetes TU Cbf
********************************************/
EbErrorType Av1TuCalcCost(
    ModeDecisionCandidate_t *candidate_ptr,                        // input parameter, prediction result Ptr
    int16_t                   txbSkipCtx,
    uint32_t                   tu_index,                             // input parameter, TU index inside the CU
    uint32_t                   y_count_non_zero_coeffs,                 // input parameter, number of non zero Y quantized coefficients
    uint32_t                   cbCountNonZeroCoeffs,                // input parameter, number of non zero cb quantized coefficients
    uint32_t                   crCountNonZeroCoeffs,                // input parameter, number of non zero cr quantized coefficients
    uint64_t                   yTuDistortion[DIST_CALC_TOTAL],      // input parameter, Y distortion for both Normal and Cbf zero modes
    uint64_t                   cbTuDistortion[DIST_CALC_TOTAL],     // input parameter, Cb distortion for both Normal and Cbf zero modes
    uint64_t                   crTuDistortion[DIST_CALC_TOTAL],     // input parameter, Cr distortion for both Normal and Cbf zero modes
    COMPONENT_TYPE           componentType,
    uint64_t                  *yTuCoeffBits,                        // input parameter, Y quantized coefficients rate
    uint64_t                  *cbTuCoeffBits,                       // input parameter, Cb quantized coefficients rate
    uint64_t                  *crTuCoeffBits,                       // input parameter, Cr quantized coefficients rate
    TxSize                  txsize,
    uint64_t                   lambda)                              // input parameter, lambda for Luma

{
    (void)crTuCoeffBits;
    (void)cbTuCoeffBits;
    (void)crTuDistortion;
    (void)cbTuDistortion;
    EbErrorType return_error = EB_ErrorNone;
    // Non Zero coeff mode variables
    uint64_t y_nonzero_coeff_distortion = yTuDistortion[DIST_CALC_RESIDUAL];
    uint64_t y_nonzero_coeff_rate;

    uint64_t y_nonzero_coeff_cost = 0;

    // Zero Cbf mode variables
    uint64_t y_zero_coeff_distortion = yTuDistortion[DIST_CALC_PREDICTION];

    uint64_t y_zero_coeff_luma_flag_bits_num = 0;

    uint64_t y_zero_coeff_rate;

    uint64_t y_zero_coeff_cost = 0;
    if (componentType == COMPONENT_LUMA || componentType == COMPONENT_ALL) {

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
        const LV_MAP_COEFF_COST *const coeff_costs = &candidate_ptr->md_rate_estimation_ptr->coeffFacBits[txs_ctx][0];

        y_zero_coeff_luma_flag_bits_num = coeff_costs->txb_skip_cost[txbSkipCtx][1];

        y_nonzero_coeff_rate = *yTuCoeffBits; // yNonZeroCbfLumaFlagBitsNum is already calculated inside yTuCoeffBits

        y_zero_coeff_rate = y_zero_coeff_luma_flag_bits_num;

#if CBF_ZERO_OFF
        if (1) {
#else
        if (candidate_ptr->type == INTRA_MODE) {
#endif
            y_zero_coeff_cost = 0xFFFFFFFFFFFFFFFFull;

        }
        else {

            y_zero_coeff_cost = RDCOST(lambda, y_zero_coeff_rate, y_zero_coeff_distortion);
        }

        // **Compute Cost
        y_nonzero_coeff_cost = RDCOST(lambda, y_nonzero_coeff_rate, y_nonzero_coeff_distortion);

        candidate_ptr->y_has_coeff |= (((y_count_non_zero_coeffs != 0) && (y_nonzero_coeff_cost < y_zero_coeff_cost)) << tu_index);
        *yTuCoeffBits = (y_nonzero_coeff_cost < y_zero_coeff_cost) ? *yTuCoeffBits : 0;
        yTuDistortion[DIST_CALC_RESIDUAL] = (y_nonzero_coeff_cost < y_zero_coeff_cost) ? yTuDistortion[DIST_CALC_RESIDUAL] : yTuDistortion[DIST_CALC_PREDICTION];
        }
    if (componentType == COMPONENT_CHROMA_CB || componentType == COMPONENT_CHROMA || componentType == COMPONENT_ALL) {

        candidate_ptr->u_has_coeff |= ((cbCountNonZeroCoeffs != 0) << tu_index);
    }
    if (componentType == COMPONENT_CHROMA_CR || componentType == COMPONENT_CHROMA || componentType == COMPONENT_ALL) {

        candidate_ptr->v_has_coeff |= ((crCountNonZeroCoeffs != 0) << tu_index);
    }

    return return_error;
    }

/********************************************
* TuCalcCost
*   computes TU Cost and generetes TU Cbf
********************************************/

EbErrorType Av1TuCalcCostLuma(

    int16_t                   txbSkipCtx,
    ModeDecisionCandidate_t *candidate_ptr,                        // input parameter, prediction result Ptr
    uint32_t                   tu_index,                             // input parameter, TU index inside the CU
    TxSize                  txSize,
    uint32_t                   y_count_non_zero_coeffs,                 // input parameter, number of non zero Y quantized coefficients
    uint64_t                   yTuDistortion[DIST_CALC_TOTAL],      // input parameter, Y distortion for both Normal and Cbf zero modes
    uint64_t                  *yTuCoeffBits,                        // input parameter, Y quantized coefficients rate
    uint64_t                  *yFullCost,
    uint64_t                   lambda)                              // input parameter, lambda for Luma

{

    EbErrorType return_error = EB_ErrorNone;

    // Non Zero Cbf mode variables
    uint64_t yNonZeroCbfDistortion = yTuDistortion[DIST_CALC_RESIDUAL];

    uint64_t yNonZeroCbfRate;

    uint64_t yNonZeroCbfCost = 0;

    // Zero Cbf mode variables
    uint64_t yZeroCbfDistortion = yTuDistortion[DIST_CALC_PREDICTION];

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

    const TxSize txs_ctx = (TxSize)((txsize_sqr_map[txSize] + txsize_sqr_up_map[txSize] + 1) >> 1);
    const LV_MAP_COEFF_COST *const coeff_costs = &candidate_ptr->md_rate_estimation_ptr->coeffFacBits[txs_ctx][0];

    yZeroCbfLumaFlagBitsNum = coeff_costs->txb_skip_cost[txbSkipCtx][1];

    yNonZeroCbfRate = *yTuCoeffBits; // yNonZeroCbfLumaFlagBitsNum is already calculated inside yTuCoeffBits

    yZeroCbfRate = yZeroCbfLumaFlagBitsNum;

#if CBF_ZERO_OFF
    if (1) {
#else
    if (candidate_ptr->type == INTRA_MODE) {
#endif
        yZeroCbfCost = 0xFFFFFFFFFFFFFFFFull;

    }
    else {

        yZeroCbfCost = RDCOST(lambda, yZeroCbfRate, yZeroCbfDistortion);
    }

    // **Compute Cost
    yNonZeroCbfCost = RDCOST(lambda, yNonZeroCbfRate, yNonZeroCbfDistortion);
    candidate_ptr->y_has_coeff |= ((y_count_non_zero_coeffs != 0) << tu_index);
    *yTuCoeffBits = (yNonZeroCbfCost < yZeroCbfCost) ? *yTuCoeffBits : 0;
    yTuDistortion[DIST_CALC_RESIDUAL] = (yNonZeroCbfCost < yZeroCbfCost) ? yTuDistortion[DIST_CALC_RESIDUAL] : yTuDistortion[DIST_CALC_PREDICTION];

    *yFullCost = MIN(yNonZeroCbfCost, yZeroCbfCost);

    return return_error;
    }

static INLINE int32_t partition_cdf_length(BlockSize bsize) {
    if (bsize <= BLOCK_8X8)
        return PARTITION_TYPES;
    else if (bsize == BLOCK_128X128)
        return EXT_PARTITION_TYPES - 2;
    else
        return EXT_PARTITION_TYPES;
}

static int32_t cdf_element_prob(const int32_t *cdf,
    size_t element) {
    assert(cdf != NULL);
    return (element > 0 ? cdf[element - 1] : CDF_PROB_TOP) - cdf[element];
}
static void partition_gather_horz_alike(int32_t *out,
    BlockSize bsize,
    const int32_t *const in) {
    out[0] = CDF_PROB_TOP;
    out[0] -= cdf_element_prob(in, PARTITION_HORZ);
    out[0] -= cdf_element_prob(in, PARTITION_SPLIT);
    out[0] -= cdf_element_prob(in, PARTITION_HORZ_A);
    out[0] -= cdf_element_prob(in, PARTITION_HORZ_B);
    out[0] -= cdf_element_prob(in, PARTITION_VERT_A);
    if (bsize != BLOCK_128X128) out[0] -= cdf_element_prob(in, PARTITION_HORZ_4);
    out[0] = AOM_ICDF(out[0]);
    out[1] = AOM_ICDF(CDF_PROB_TOP);
}

static void partition_gather_vert_alike(int32_t *out,
    BlockSize bsize,
    const int32_t *const in) {
    out[0] = CDF_PROB_TOP;
    out[0] -= cdf_element_prob(in, PARTITION_VERT);
    out[0] -= cdf_element_prob(in, PARTITION_SPLIT);
    out[0] -= cdf_element_prob(in, PARTITION_HORZ_A);
    out[0] -= cdf_element_prob(in, PARTITION_VERT_A);
    out[0] -= cdf_element_prob(in, PARTITION_VERT_B);
    if (bsize != BLOCK_128X128) out[0] -= cdf_element_prob(in, PARTITION_VERT_4);
    out[0] = AOM_ICDF(out[0]);
    out[1] = AOM_ICDF(CDF_PROB_TOP);
}

//static INLINE int32_t partition_plane_context(const MacroBlockD *xd, int32_t mi_row,
//    int32_t mi_col, BlockSize bsize) {
//    const PARTITION_CONTEXT *above_ctx = xd->above_seg_context + mi_col;
//    const PARTITION_CONTEXT *left_ctx =
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
* SplitFlagRate function is used to generate the Split rate
*
*   @param *cu_ptr(input)
*       cu_ptr is the pointer of the target CU.
*   @param split_flag(input)
*       split_flag is the split flag value.
*   @param splitRate(output)
*       splitRate contains rate.
*   @param lambda(input)
*       lambda is the Lagrange multiplier
*   @param md_rate_estimation_ptr(input)
*       md_rate_estimation_ptr is pointer to MD rate Estimation Tables
**********************************************************************************/
EbErrorType Av1SplitFlagRate(
    SequenceControlSet_t                  *sequence_control_set_ptr,
    ModeDecisionContext_t                  *context_ptr,
    CodingUnit_t                           *cu_ptr,
    uint32_t                                  leaf_index,
    PartitionType                          partitionType,
    uint64_t                                 *splitRate,
    uint64_t                                  lambda,
    MdRateEstimationContext_t              *md_rate_estimation_ptr,
    uint32_t                                  tbMaxDepth)
{
    (void)tbMaxDepth;
    (void)leaf_index;

    const BlockGeom          *blk_geom = Get_blk_geom_mds(cu_ptr->mds_idx);
    EbErrorType return_error = EB_ErrorNone;

    uint32_t cu_origin_x = context_ptr->sb_origin_x + blk_geom->origin_x;
    uint32_t cu_origin_y = context_ptr->sb_origin_y + blk_geom->origin_y;

    PartitionType p = partitionType;

    uint32_t cu_depth = blk_geom->depth;
    UNUSED(cu_depth);
    BlockSize bsize = blk_geom->bsize;
    ASSERT(bsize<BlockSizeS_ALL);
    const int32_t is_partition_point = blk_geom->bsize >= BLOCK_8X8;


    if (is_partition_point) {

        const int32_t hbs = (mi_size_wide[bsize] << 2) >> 1;
        const int32_t hasRows = (cu_origin_y + hbs) < sequence_control_set_ptr->luma_height;
        const int32_t hasCols = (cu_origin_x + hbs) < sequence_control_set_ptr->luma_width;

        uint32_t contextIndex = 0;

        const PARTITION_CONTEXT left_ctx = context_ptr->md_local_cu_unit[cu_ptr->mds_idx].left_neighbor_partition == (int8_t)(INVALID_NEIGHBOR_DATA) ? 0 : context_ptr->md_local_cu_unit[cu_ptr->mds_idx].left_neighbor_partition;
        const PARTITION_CONTEXT above_ctx = context_ptr->md_local_cu_unit[cu_ptr->mds_idx].above_neighbor_partition == (int8_t)(INVALID_NEIGHBOR_DATA) ? 0 : context_ptr->md_local_cu_unit[cu_ptr->mds_idx].above_neighbor_partition;

        const int32_t bsl = mi_size_wide_log2[bsize] - mi_size_wide_log2[BLOCK_8X8];

        int32_t above = (above_ctx >> bsl) & 1, left = (left_ctx >> bsl) & 1;

        assert(mi_size_wide_log2[bsize] == mi_size_high_log2[bsize]);
        assert(bsl >= 0);

        contextIndex = (left * 2 + above) + bsl * PARTITION_PLOFFSET;

        if (hasRows && hasCols) {

            *splitRate = (uint64_t)md_rate_estimation_ptr->partitionFacBits[partition_cdf_length(bsize)][partitionType];

        }
        else if (!hasRows && hasCols) {
            int32_t cdf[2];
            partition_gather_vert_alike(cdf, bsize, md_rate_estimation_ptr->partitionFacBits[contextIndex]);
            *splitRate = (uint64_t)md_rate_estimation_ptr->partitionFacBits[partition_cdf_length(bsize)][partitionType];

            *splitRate = (uint64_t)cdf[p == PARTITION_SPLIT];
        }
        else {
            int32_t cdf[2];
            partition_gather_horz_alike(cdf, bsize, md_rate_estimation_ptr->partitionFacBits[contextIndex]);
            *splitRate = (uint64_t)cdf[p == PARTITION_SPLIT];
        }
    }
    else {
        *splitRate = (uint64_t)md_rate_estimation_ptr->partitionFacBits[0][partitionType];
    }

    *splitRate = RDCOST(lambda, *splitRate, 0);

    return return_error;
}

/********************************************
* TuCalcCost
*   Computes TU Cost and generetes TU Cbf
*   at the level of the encode pass
********************************************/
EbErrorType Av1EncodeTuCalcCost(
    EncDecContext_t          *context_ptr,
    uint32_t                   *count_non_zero_coeffs,
    uint64_t                    yTuDistortion[DIST_CALC_TOTAL],
    uint64_t                   *yTuCoeffBits,
    uint32_t                    component_mask
)
{
    CodingUnit_t              *cu_ptr = context_ptr->cu_ptr;
    uint32_t                     tu_index = context_ptr->txb_itr;
    MdRateEstimationContext_t *md_rate_estimation_ptr = context_ptr->md_rate_estimation_ptr;
    uint64_t                     lambda = context_ptr->full_lambda;
    uint32_t                     y_count_non_zero_coeffs = count_non_zero_coeffs[0];
    uint32_t                     cbCountNonZeroCoeffs = count_non_zero_coeffs[1];
    uint32_t                     crCountNonZeroCoeffs = count_non_zero_coeffs[2];

    EbErrorType return_error = EB_ErrorNone;

    // Non Zero Cbf mode variables
    uint64_t yNonZeroCbfDistortion = yTuDistortion[DIST_CALC_RESIDUAL];

    uint64_t yNonZeroCbfRate;

    uint64_t yNonZeroCbfCost = 0;

    // Zero Cbf mode variables
    uint64_t yZeroCbfDistortion = yTuDistortion[DIST_CALC_PREDICTION];

    uint64_t yZeroCbfLumaFlagBitsNum = 0;

    uint64_t yZeroCbfRate;

    uint64_t yZeroCbfCost = 0;
    int16_t  txbSkipCtx = cu_ptr->luma_txb_skip_context;

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
        TxSize    txSize = context_ptr->blk_geom->txsize[context_ptr->txb_itr];

        const TxSize txs_ctx = (TxSize)((txsize_sqr_map[txSize] + txsize_sqr_up_map[txSize] + 1) >> 1);
        const LV_MAP_COEFF_COST *const coeff_costs = &md_rate_estimation_ptr->coeffFacBits[txs_ctx][0];

        yZeroCbfLumaFlagBitsNum = coeff_costs->txb_skip_cost[txbSkipCtx][1];

        yNonZeroCbfRate = *yTuCoeffBits; // yNonZeroCbfLumaFlagBitsNum is already calculated inside yTuCoeffBits

        yZeroCbfRate = yZeroCbfLumaFlagBitsNum;
#if CBF_ZERO_OFF || TX_TYPE_FIX
        if (1) {
#else
        if (cu_ptr->prediction_mode_flag == INTRA_MODE) {
#endif
            yZeroCbfCost = 0xFFFFFFFFFFFFFFFFull;

        }
        else {

            yZeroCbfCost = RDCOST(lambda, yZeroCbfRate, yZeroCbfDistortion);
        }

        // **Compute Cost
        yNonZeroCbfCost = RDCOST(lambda, yNonZeroCbfRate, yNonZeroCbfDistortion);
        cu_ptr->transform_unit_array[tu_index].y_has_coeff = ((y_count_non_zero_coeffs != 0) && (yNonZeroCbfCost < yZeroCbfCost)) ? EB_TRUE : EB_FALSE;
        *yTuCoeffBits = (yNonZeroCbfCost < yZeroCbfCost) ? *yTuCoeffBits : 0;
        yTuDistortion[DIST_CALC_RESIDUAL] = (yNonZeroCbfCost < yZeroCbfCost) ? yTuDistortion[DIST_CALC_RESIDUAL] : yTuDistortion[DIST_CALC_PREDICTION];

        }
    else {
        cu_ptr->transform_unit_array[tu_index].y_has_coeff = EB_FALSE;
    }
    cu_ptr->transform_unit_array[tu_index].u_has_coeff = cbCountNonZeroCoeffs != 0 ? EB_TRUE : EB_FALSE;
    cu_ptr->transform_unit_array[tu_index].v_has_coeff = crCountNonZeroCoeffs != 0 ? EB_TRUE : EB_FALSE;

    return return_error;
    }


uint64_t GetPMCost(
    uint64_t                   lambda,
    uint64_t                   tuDistortion,
    uint64_t                   yTuCoeffBits
)
{

    uint64_t yNonZeroCbfDistortion = LUMA_WEIGHT * (tuDistortion << COST_PRECISION);
    uint64_t yNonZeroCbfRate = (yTuCoeffBits);
    uint64_t yNonZeroCbfCost = yNonZeroCbfDistortion + (((lambda       * yNonZeroCbfRate) + MD_OFFSET) >> MD_SHIFT);

    return yNonZeroCbfCost;
}

