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

#include <stdlib.h>

#include "EbDefinitions.h"
#include "EbMdRateEstimation.h"

#include "EbBitstreamUnit.h"

static INLINE int32_t get_interinter_wedge_bits(BlockSize sb_type) {
    const int32_t wbits = wedge_params_lookup[sb_type].bits;
    return (wbits > 0) ? wbits + 1 : 0;
}

/**************************************************************
* AV1GetCostSymbold
* Calculate the cost of a symbol with
* probability p15 / 2^15
***************************************************************/
static INLINE int32_t av1_cost_symbol(AomCdfProb p15) {
    assert(0 < p15 && p15 < CDF_PROB_TOP);
    const int32_t shift = CDF_PROB_BITS - 1 - get_msb(p15);
    const int32_t prob = get_prob(p15 << shift, CDF_PROB_TOP);
    assert(prob >= 128);
    return av1_prob_cost[prob - 128] + av1_cost_literal(shift);
}

/*************************************************************
* av1_get_syntax_rate_from_cdf
**************************************************************/
void av1_get_syntax_rate_from_cdf(
    int32_t                      *costs,
    const AomCdfProb       *cdf,
    const int32_t                *inv_map)
{
    int32_t i;
    AomCdfProb prev_cdf = 0;
    for (i = 0;; ++i) {
        AomCdfProb p15 = AOM_ICDF(cdf[i]) - prev_cdf;
        p15 = (p15 < EC_MIN_PROB) ? EC_MIN_PROB : p15;
        prev_cdf = AOM_ICDF(cdf[i]);

        if (inv_map)
            costs[inv_map[i]] = av1_cost_symbol(p15);
        else
            costs[i] = av1_cost_symbol(p15);

        // Stop once we reach the end of the CDF
        if (cdf[i] == AOM_ICDF(CDF_PROB_TOP)) break;
    }
}

///tmp function to be removed once we have updated all syntax CDFs
void av1_estimate_syntax_rate___partial(
    MdRateEstimationContext  *md_rate_estimation_array,
    FRAME_CONTEXT              *fc)
{
    int32_t i, j;

    md_rate_estimation_array->initialized = 1;
#if CABAC_UP1
    for (i = 0; i < PARTITION_CONTEXTS; ++i)
        av1_get_syntax_rate_from_cdf(md_rate_estimation_array->partitionFacBits[i], fc->partition_cdf[i], NULL);
#endif

#if CABAC_UP2
    //if (cm->skip_mode_flag) { // NM - Hardcoded to true
    for (i = 0; i < SKIP_CONTEXTS; ++i)
        av1_get_syntax_rate_from_cdf(md_rate_estimation_array->skipModeFacBits[i], fc->skip_mode_cdfs[i], NULL);
    //}
#endif

    for (i = TX_4X4; i < EXT_TX_SIZES; ++i) {
        int32_t s;
        for (s = 1; s < EXT_TX_SETS_INTER; ++s) {
            if (use_inter_ext_tx_for_txsize[s][i])
                av1_get_syntax_rate_from_cdf(md_rate_estimation_array->inter_tx_type_fac_bits[s][i], fc->inter_ext_tx_cdf[s][i], av1_ext_tx_inv[av1_ext_tx_set_idx_to_type[1][s]]);
        }
        for (s = 1; s < EXT_TX_SETS_INTRA; ++s) {
            if (use_intra_ext_tx_for_txsize[s][i]) {
                for (j = 0; j < INTRA_MODES; ++j)
                    av1_get_syntax_rate_from_cdf(md_rate_estimation_array->intra_tx_type_fac_bits[s][i][j], fc->intra_ext_tx_cdf[s][i][j], av1_ext_tx_inv[av1_ext_tx_set_idx_to_type[0][s]]);
            }
        }
    }
}
#if FILTER_INTRA_FLAG
int av1_filter_intra_allowed_bsize(  uint8_t enable_filter_intra,  BlockSize bs);
#if !PAL_SUP
int av1_filter_intra_allowed(uint8_t   enable_filter_intra, BlockSize bsize, uint32_t  mode);
#endif
#endif
/*************************************************************
* av1_estimate_syntax_rate()
* Estimate the rate for each syntax elements and for
* all scenarios based on the frame CDF
**************************************************************/
void av1_estimate_syntax_rate(
    MdRateEstimationContext  *md_rate_estimation_array,
    EbBool                     is_i_slice,
    FRAME_CONTEXT              *fc)
{
    int32_t i, j;

    md_rate_estimation_array->initialized = 1;

    for (i = 0; i < PARTITION_CONTEXTS; ++i)
        av1_get_syntax_rate_from_cdf(md_rate_estimation_array->partition_fac_bits[i], fc->partition_cdf[i], NULL);

    //if (cm->skip_mode_flag) { // NM - Hardcoded to true
    for (i = 0; i < SKIP_CONTEXTS; ++i)
        av1_get_syntax_rate_from_cdf(md_rate_estimation_array->skip_mode_fac_bits[i], fc->skip_mode_cdfs[i], NULL);
    //}

    for (i = 0; i < SKIP_CONTEXTS; ++i)
        av1_get_syntax_rate_from_cdf(md_rate_estimation_array->skip_fac_bits[i], fc->skip_cdfs[i], NULL);
    for (i = 0; i < KF_MODE_CONTEXTS; ++i)
        for (j = 0; j < KF_MODE_CONTEXTS; ++j)
            av1_get_syntax_rate_from_cdf(md_rate_estimation_array->y_mode_fac_bits[i][j], fc->kf_y_cdf[i][j], NULL);

    for (i = 0; i < BlockSize_GROUPS; ++i)
        av1_get_syntax_rate_from_cdf(md_rate_estimation_array->mb_mode_fac_bits[i], fc->y_mode_cdf[i], NULL);

    for (i = 0; i < CFL_ALLOWED_TYPES; ++i) {
        for (j = 0; j < INTRA_MODES; ++j)
            av1_get_syntax_rate_from_cdf(md_rate_estimation_array->intra_uv_mode_fac_bits[i][j], fc->uv_mode_cdf[i][j], NULL);
    }

    av1_get_syntax_rate_from_cdf(md_rate_estimation_array->filter_intra_mode_fac_bits, fc->filter_intra_mode_cdf, NULL);
#if FILTER_INTRA_FLAG
    for (i = 0; i < BlockSizeS_ALL; ++i) {
        if (av1_filter_intra_allowed_bsize(1,i))
            av1_get_syntax_rate_from_cdf(md_rate_estimation_array->filter_intra_fac_bits[i], fc->filter_intra_cdfs[i], NULL);
    }
#else
    // NM - To be added when intra filtering is adopted
    /*for (i = 0; i < BlockSizeS_ALL; ++i) {
        if (av1_filter_intra_allowed_bsize(cm, i))
            av1_FacBits_tokens_from_cdf(md_rate_estimation_array->filter_intra_fac_bits[i],
            fc->filter_intra_cdfs[i], NULL);
    }*/

    // NM - To be added when inter filtering is adopted
#endif
    for (i = 0; i < SWITCHABLE_FILTER_CONTEXTS; ++i)
        av1_get_syntax_rate_from_cdf(md_rate_estimation_array->switchable_interp_fac_bitss[i], fc->switchable_interp_cdf[i], NULL);

    for (i = 0; i < PALATTE_BSIZE_CTXS; ++i) {
        av1_get_syntax_rate_from_cdf(md_rate_estimation_array->palette_ysize_fac_bits[i], fc->palette_y_size_cdf[i], NULL);
        av1_get_syntax_rate_from_cdf(md_rate_estimation_array->palette_uv_size_fac_bits[i], fc->palette_uv_size_cdf[i], NULL);
        for (j = 0; j < PALETTE_Y_MODE_CONTEXTS; ++j)
            av1_get_syntax_rate_from_cdf(md_rate_estimation_array->palette_ymode_fac_bits[i][j], fc->palette_y_mode_cdf[i][j], NULL);
    }

    for (i = 0; i < PALETTE_UV_MODE_CONTEXTS; ++i)
        av1_get_syntax_rate_from_cdf(md_rate_estimation_array->palette_uv_mode_fac_bits[i], fc->palette_uv_mode_cdf[i], NULL);
    for (i = 0; i < PALETTE_SIZES; ++i) {
        for (j = 0; j < PALETTE_COLOR_INDEX_CONTEXTS; ++j) {
            av1_get_syntax_rate_from_cdf(md_rate_estimation_array->palette_ycolor_fac_bitss[i][j], fc->palette_y_color_index_cdf[i][j], NULL);
            av1_get_syntax_rate_from_cdf(md_rate_estimation_array->palette_uv_color_fac_bits[i][j], fc->palette_uv_color_index_cdf[i][j], NULL);
        }
    }

    int32_t sign_FacBits[CFL_JOINT_SIGNS];
    av1_get_syntax_rate_from_cdf(sign_FacBits, fc->cfl_sign_cdf, NULL);
    for (int32_t joint_sign = 0; joint_sign < CFL_JOINT_SIGNS; joint_sign++) {
        int32_t *FacBits_u = md_rate_estimation_array->cfl_alpha_fac_bits[joint_sign][CFL_PRED_U];
        int32_t *FacBits_v = md_rate_estimation_array->cfl_alpha_fac_bits[joint_sign][CFL_PRED_V];
        if (CFL_SIGN_U(joint_sign) == CFL_SIGN_ZERO)
            memset(FacBits_u, 0, CFL_ALPHABET_SIZE * sizeof(*FacBits_u));
        else {
            const AomCdfProb *cdf_u = fc->cfl_alpha_cdf[CFL_CONTEXT_U(joint_sign)];
            av1_get_syntax_rate_from_cdf(FacBits_u, cdf_u, NULL);
        }
        if (CFL_SIGN_V(joint_sign) == CFL_SIGN_ZERO)
            memset(FacBits_v, 0, CFL_ALPHABET_SIZE * sizeof(*FacBits_v));
        else {
            int32_t cdf_index = CFL_CONTEXT_V(joint_sign);
            if ((cdf_index < CFL_ALPHA_CONTEXTS) && (cdf_index >= 0)) {
                const AomCdfProb *cdf_v = fc->cfl_alpha_cdf[cdf_index];
                av1_get_syntax_rate_from_cdf(FacBits_v, cdf_v, NULL);
            }
        }
        for (int32_t u = 0; u < CFL_ALPHABET_SIZE; u++)
            FacBits_u[u] += sign_FacBits[joint_sign];
    }

    for (i = 0; i < MAX_TX_CATS; ++i)
        for (j = 0; j < TX_SIZE_CONTEXTS; ++j)
            av1_get_syntax_rate_from_cdf(md_rate_estimation_array->tx_size_fac_bits[i][j], fc->tx_size_cdf[i][j],
                NULL);

    for (i = 0; i < TXFM_PARTITION_CONTEXTS; ++i) {
        av1_get_syntax_rate_from_cdf(md_rate_estimation_array->txfm_partition_fac_bits[i],
            fc->txfm_partition_cdf[i], NULL);
    }

    for (i = TX_4X4; i < EXT_TX_SIZES; ++i) {
        int32_t s;
        for (s = 1; s < EXT_TX_SETS_INTER; ++s) {
            if (use_inter_ext_tx_for_txsize[s][i])
                av1_get_syntax_rate_from_cdf(md_rate_estimation_array->inter_tx_type_fac_bits[s][i], fc->inter_ext_tx_cdf[s][i], av1_ext_tx_inv[av1_ext_tx_set_idx_to_type[1][s]]);
        }
        for (s = 1; s < EXT_TX_SETS_INTRA; ++s) {
            if (use_intra_ext_tx_for_txsize[s][i]) {
                for (j = 0; j < INTRA_MODES; ++j)
                    av1_get_syntax_rate_from_cdf(md_rate_estimation_array->intra_tx_type_fac_bits[s][i][j], fc->intra_ext_tx_cdf[s][i][j], av1_ext_tx_inv[av1_ext_tx_set_idx_to_type[0][s]]);
            }
        }
    }
    for (i = 0; i < DIRECTIONAL_MODES; ++i)
        av1_get_syntax_rate_from_cdf(md_rate_estimation_array->angle_delta_fac_bits[i], fc->angle_delta_cdf[i], NULL);
    av1_get_syntax_rate_from_cdf(md_rate_estimation_array->switchable_restore_fac_bits, fc->switchable_restore_cdf, NULL);
    av1_get_syntax_rate_from_cdf(md_rate_estimation_array->wiener_restore_fac_bits, fc->wiener_restore_cdf, NULL);
    av1_get_syntax_rate_from_cdf(md_rate_estimation_array->sgrproj_restore_fac_bits, fc->sgrproj_restore_cdf, NULL);
    av1_get_syntax_rate_from_cdf(md_rate_estimation_array->intrabc_fac_bits, fc->intrabc_cdf, NULL);

    if (!is_i_slice) { // NM - Hardcoded to true
    //if (1){
        for (i = 0; i < COMP_INTER_CONTEXTS; ++i)
            av1_get_syntax_rate_from_cdf(md_rate_estimation_array->comp_inter_fac_bits[i], fc->comp_inter_cdf[i], NULL);
        for (i = 0; i < REF_CONTEXTS; ++i) {
            for (j = 0; j < SINGLE_REFS - 1; ++j)
                av1_get_syntax_rate_from_cdf(md_rate_estimation_array->single_ref_fac_bits[i][j], fc->single_ref_cdf[i][j], NULL);
        }

        for (i = 0; i < COMP_REF_TYPE_CONTEXTS; ++i)
            av1_get_syntax_rate_from_cdf(md_rate_estimation_array->comp_ref_type_fac_bits[i], fc->comp_ref_type_cdf[i], NULL);
        for (i = 0; i < UNI_COMP_REF_CONTEXTS; ++i) {
            for (j = 0; j < UNIDIR_COMP_REFS - 1; ++j)
                av1_get_syntax_rate_from_cdf(md_rate_estimation_array->uni_comp_ref_fac_bits[i][j], fc->uni_comp_ref_cdf[i][j], NULL);
        }

        for (i = 0; i < REF_CONTEXTS; ++i) {
            for (j = 0; j < FWD_REFS - 1; ++j)
                av1_get_syntax_rate_from_cdf(md_rate_estimation_array->comp_ref_fac_bits[i][j], fc->comp_ref_cdf[i][j], NULL);
        }

        for (i = 0; i < REF_CONTEXTS; ++i) {
            for (j = 0; j < BWD_REFS - 1; ++j)
                av1_get_syntax_rate_from_cdf(md_rate_estimation_array->comp_bwd_ref_fac_bits[i][j], fc->comp_bwdref_cdf[i][j], NULL);
        }

        for (i = 0; i < INTRA_INTER_CONTEXTS; ++i)
            av1_get_syntax_rate_from_cdf(md_rate_estimation_array->intra_inter_fac_bits[i], fc->intra_inter_cdf[i], NULL);
        for (i = 0; i < NEWMV_MODE_CONTEXTS; ++i)
            av1_get_syntax_rate_from_cdf(md_rate_estimation_array->new_mv_mode_fac_bits[i], fc->newmv_cdf[i], NULL);
        for (i = 0; i < GLOBALMV_MODE_CONTEXTS; ++i)
            av1_get_syntax_rate_from_cdf(md_rate_estimation_array->zero_mv_mode_fac_bits[i], fc->zeromv_cdf[i], NULL);
        for (i = 0; i < REFMV_MODE_CONTEXTS; ++i)
            av1_get_syntax_rate_from_cdf(md_rate_estimation_array->ref_mv_mode_fac_bits[i], fc->refmv_cdf[i], NULL);
        for (i = 0; i < DRL_MODE_CONTEXTS; ++i)
            av1_get_syntax_rate_from_cdf(md_rate_estimation_array->drl_mode_fac_bits[i], fc->drl_cdf[i], NULL);
        for (i = 0; i < INTER_MODE_CONTEXTS; ++i)
            av1_get_syntax_rate_from_cdf(md_rate_estimation_array->inter_compound_mode_fac_bits[i], fc->inter_compound_mode_cdf[i], NULL);
        for (i = 0; i < BlockSizeS_ALL; ++i)
            av1_get_syntax_rate_from_cdf(md_rate_estimation_array->compound_type_fac_bits[i], fc->compound_type_cdf[i], NULL);
        for (i = 0; i < BlockSizeS_ALL; ++i) {
            if (get_interinter_wedge_bits((BlockSize)i))
                av1_get_syntax_rate_from_cdf(md_rate_estimation_array->wedge_idx_fac_bits[i], fc->wedge_idx_cdf[i], NULL);
        }
        for (i = 0; i < BlockSize_GROUPS; ++i) {
            av1_get_syntax_rate_from_cdf(md_rate_estimation_array->inter_intra_fac_bits[i], fc->interintra_cdf[i], NULL);
            av1_get_syntax_rate_from_cdf(md_rate_estimation_array->inter_intra_mode_fac_bits[i], fc->interintra_mode_cdf[i], NULL);
        }
        for (i = 0; i < BlockSizeS_ALL; ++i)
            av1_get_syntax_rate_from_cdf(md_rate_estimation_array->wedge_inter_intra_fac_bits[i], fc->wedge_interintra_cdf[i], NULL);
        for (i = BLOCK_8X8; i < BlockSizeS_ALL; i++)
            av1_get_syntax_rate_from_cdf(md_rate_estimation_array->motion_mode_fac_bits[i], fc->motion_mode_cdf[i], NULL);
        for (i = BLOCK_8X8; i < BlockSizeS_ALL; i++)
            av1_get_syntax_rate_from_cdf(md_rate_estimation_array->motion_mode_fac_bits1[i], fc->obmc_cdf[i], NULL);
        for (i = 0; i < COMP_INDEX_CONTEXTS; ++i)
            av1_get_syntax_rate_from_cdf(md_rate_estimation_array->comp_idx_fac_bits[i], fc->compound_index_cdf[i], NULL);
        for (i = 0; i < COMP_GROUP_IDX_CONTEXTS; ++i)
            av1_get_syntax_rate_from_cdf(md_rate_estimation_array->comp_group_idx_fac_bits[i], fc->comp_group_idx_cdf[i], NULL);
    }
}

static const uint8_t log_in_base_2[] = {
    0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10
};

static INLINE int32_t mv_class_base(MvClassType c) {
    return c ? CLASS0_SIZE << (c + 2) : 0;
}

MvClassType av1_get_mv_class(int32_t z, int32_t *offset) {
    const MvClassType c = (z >= CLASS0_SIZE * 4096)
        ? MV_CLASS_10
        : (MvClassType)log_in_base_2[z >> 3];
    if (offset) *offset = z - mv_class_base(c);
    return c;
}

//void eb_av1_build_nmv_cost_table(int32_t *mvjoint, int32_t *mvcost[2],
//    const NmvContext *ctx,
//    MvSubpelPrecision precision)

void eb_av1_build_nmv_cost_table(int32_t *mvjoint, int32_t *mvcost[2],
    const NmvContext *ctx,
    MvSubpelPrecision precision);

/**************************************************************************
* av1_estimate_mv_rate()
* Estimate the rate of motion vectors
* based on the frame CDF
***************************************************************************/
void av1_estimate_mv_rate(
    PictureControlSet     *picture_control_set_ptr,
    MdRateEstimationContext  *md_rate_estimation_array,
    NmvContext                *nmv_ctx)
{
    int32_t *nmvcost[2];
    int32_t *nmvcost_hp[2];
    FrameHeader *frm_hdr = &picture_control_set_ptr->parent_pcs_ptr->frm_hdr;

    nmvcost[0] = &md_rate_estimation_array->nmv_costs[0][MV_MAX];
    nmvcost[1] = &md_rate_estimation_array->nmv_costs[1][MV_MAX];
    nmvcost_hp[0] = &md_rate_estimation_array->nmv_costs_hp[0][MV_MAX];
    nmvcost_hp[1] = &md_rate_estimation_array->nmv_costs_hp[1][MV_MAX];

    eb_av1_build_nmv_cost_table(
        md_rate_estimation_array->nmv_vec_cost,//out
        frm_hdr->allow_high_precision_mv ? nmvcost_hp : nmvcost, //out
        nmv_ctx,
        frm_hdr->allow_high_precision_mv);
#if EIGHT_PEL_FIX
    md_rate_estimation_array->nmvcoststack[0] = frm_hdr->allow_high_precision_mv ?
        &md_rate_estimation_array->nmv_costs_hp[0][MV_MAX] : &md_rate_estimation_array->nmv_costs[0][MV_MAX];
    md_rate_estimation_array->nmvcoststack[1] = frm_hdr->allow_high_precision_mv ?
        &md_rate_estimation_array->nmv_costs_hp[1][MV_MAX] : &md_rate_estimation_array->nmv_costs[1][MV_MAX];
#else
    md_rate_estimation_array->nmvcoststack[0] = &md_rate_estimation_array->nmv_costs[0][MV_MAX];
    md_rate_estimation_array->nmvcoststack[1] = &md_rate_estimation_array->nmv_costs[1][MV_MAX];
#endif
    if (frm_hdr->allow_intrabc) {
        int32_t *dvcost[2] = { &md_rate_estimation_array->dv_cost[0][MV_MAX], &md_rate_estimation_array->dv_cost[1][MV_MAX] };
        eb_av1_build_nmv_cost_table(md_rate_estimation_array->dv_joint_cost, dvcost, &picture_control_set_ptr->coeff_est_entropy_coder_ptr->fc->ndvc,
            MV_SUBPEL_NONE);
    }
}
/**************************************************************************
* av1_estimate_coefficients_rate()
* Estimate the rate of the quantised coefficient
* based on the frame CDF
***************************************************************************/
void av1_estimate_coefficients_rate(
    MdRateEstimationContext  *md_rate_estimation_array,
    FRAME_CONTEXT              *fc)
{
    int32_t num_planes = 3; // NM - Hardcoded to 3
    const int32_t nplanes = AOMMIN(num_planes, PLANE_TYPES);
    int32_t eob_multi_size = 0;
    int32_t plane = 0;
    int32_t ctx = 0;
    int32_t tx_size = 0;

    for (eob_multi_size = 0; eob_multi_size < 7; ++eob_multi_size) {
        for (plane = 0; plane < nplanes; ++plane) {
            LvMapEobCost *pcost = &md_rate_estimation_array->eob_frac_bits[eob_multi_size][plane];
            for (ctx = 0; ctx < 2; ++ctx) {
                AomCdfProb *pcdf;
                switch (eob_multi_size) {
                case 0: pcdf = fc->eob_flag_cdf16[plane][ctx]; break;
                case 1: pcdf = fc->eob_flag_cdf32[plane][ctx]; break;
                case 2: pcdf = fc->eob_flag_cdf64[plane][ctx]; break;
                case 3: pcdf = fc->eob_flag_cdf128[plane][ctx]; break;
                case 4: pcdf = fc->eob_flag_cdf256[plane][ctx]; break;
                case 5: pcdf = fc->eob_flag_cdf512[plane][ctx]; break;
                case 6:
                default: pcdf = fc->eob_flag_cdf1024[plane][ctx]; break;
                }
                av1_get_syntax_rate_from_cdf(pcost->eob_cost[ctx], pcdf, NULL);
            }
        }
    }
    for (tx_size = 0; tx_size < TX_SIZES; ++tx_size) {
        for (plane = 0; plane < nplanes; ++plane) {
            LvMapCoeffCost *pcost = &md_rate_estimation_array->coeff_fac_bits[tx_size][plane];

            for (ctx = 0; ctx < TXB_SKIP_CONTEXTS; ++ctx)
                av1_get_syntax_rate_from_cdf(pcost->txb_skip_cost[ctx],
                    fc->txb_skip_cdf[tx_size][ctx], NULL);

            for (ctx = 0; ctx < SIG_COEF_CONTEXTS_EOB; ++ctx)
                av1_get_syntax_rate_from_cdf(pcost->base_eob_cost[ctx],
                    fc->coeff_base_eob_cdf[tx_size][plane][ctx],
                    NULL);
            for (ctx = 0; ctx < SIG_COEF_CONTEXTS; ++ctx)
                av1_get_syntax_rate_from_cdf(pcost->base_cost[ctx],
                    fc->coeff_base_cdf[tx_size][plane][ctx], NULL);
            for (int ctx = 0; ctx < SIG_COEF_CONTEXTS; ++ctx) {
                pcost->base_cost[ctx][4] = 0;
                pcost->base_cost[ctx][5] = pcost->base_cost[ctx][1] +
                    av1_cost_literal(1) -
                    pcost->base_cost[ctx][0];
                pcost->base_cost[ctx][6] =
                    pcost->base_cost[ctx][2] - pcost->base_cost[ctx][1];
                pcost->base_cost[ctx][7] =
                    pcost->base_cost[ctx][3] - pcost->base_cost[ctx][2];
            }
            for (ctx = 0; ctx < EOB_COEF_CONTEXTS; ++ctx)
                av1_get_syntax_rate_from_cdf(pcost->eob_extra_cost[ctx],
                    fc->eob_extra_cdf[tx_size][plane][ctx], NULL);

            for (ctx = 0; ctx < DC_SIGN_CONTEXTS; ++ctx)
                av1_get_syntax_rate_from_cdf(pcost->dc_sign_cost[ctx],
                    fc->dc_sign_cdf[plane][ctx], NULL);

            for (ctx = 0; ctx < LEVEL_CONTEXTS; ++ctx) {
                int32_t br_rate[BR_CDF_SIZE];
                int32_t prev_cost = 0;
                int32_t i, j;
                av1_get_syntax_rate_from_cdf(br_rate, fc->coeff_br_cdf[tx_size][plane][ctx], NULL);
                // printf("br_rate: ");
                // for(j = 0; j < BR_CDF_SIZE; j++)
                //  printf("%4d ", br_rate[j]);
                // printf("\n");
                for (i = 0; i < COEFF_BASE_RANGE; i += BR_CDF_SIZE - 1) {
                    for (j = 0; j < BR_CDF_SIZE - 1; j++)
                        pcost->lps_cost[ctx][i + j] = prev_cost + br_rate[j];
                    prev_cost += br_rate[j];
                }
                pcost->lps_cost[ctx][i] = prev_cost;
                // printf("lps_cost: %d %d %2d : ", tx_size, plane, ctx);
                // for (i = 0; i <= COEFF_BASE_RANGE; i++)
                //  printf("%5d ", pcost->lps_cost[ctx][i]);
                // printf("\n");
            }
            for (int ctx = 0; ctx < LEVEL_CONTEXTS; ++ctx) {
                pcost->lps_cost[ctx][0 + COEFF_BASE_RANGE + 1] =
                    pcost->lps_cost[ctx][0];
                for (int i = 1; i <= COEFF_BASE_RANGE; ++i) {
                    pcost->lps_cost[ctx][i + COEFF_BASE_RANGE + 1] =
                        pcost->lps_cost[ctx][i] - pcost->lps_cost[ctx][i - 1];
                }
            }
        }
    }
}
