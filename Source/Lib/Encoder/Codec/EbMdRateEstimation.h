// clang-format off
/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbMdRateEstimation_h
#define EbMdRateEstimation_h

#include "EbCabacContextModel.h"
#include "EbPictureControlSet.h"
#include "EbCodingUnit.h"
#ifdef __cplusplus
extern "C" {
#endif
    /**************************************
     * MD Rate Estimation Defines
     **************************************/
#define MV_COST_WEIGHT_SUB 120

#define TOTAL_NUMBER_OF_MD_RATE_ESTIMATION_CASE_BUFFERS (TOTAL_NUMBER_OF_QP_VALUES * TOTAL_NUMBER_OF_SLICE_TYPES)
#define NUMBER_OF_SPLIT_FLAG_CASES                            6       // number of cases for bit estimation for split flag
#define NUMBER_OF_MVD_CASES                                  12       // number of cases for bit estimation for motion vector difference
     // Set to (1 << 5) if the 32-ary codebooks are used for any bock size
#define MAX_WEDGE_TYPES                                      (1 << 4)
     // The factor to scale from cost in bits to cost in av1_prob_cost units.
#define AV1_PROB_COST_SHIFT 9
     // Cost of coding an n bit literal, using 128 (i.e. 50%) probability for each bit.
#define av1_cost_literal(n) ((n) * (1 << AV1_PROB_COST_SHIFT))

    typedef struct LvMapEobCost {
        int32_t eob_cost[2][11];
    } LvMapEobCost;

    typedef struct LvMapCoeffCost
    {
        int32_t txb_skip_cost[TXB_SKIP_CONTEXTS][2];
        int32_t base_eob_cost[SIG_COEF_CONTEXTS_EOB][3];
        int32_t base_cost[SIG_COEF_CONTEXTS][8];
        int32_t eob_extra_cost[EOB_COEF_CONTEXTS][2];
        int32_t dc_sign_cost[DC_SIGN_CONTEXTS][2];
        int32_t lps_cost[LEVEL_CONTEXTS][COEFF_BASE_RANGE + 1 + COEFF_BASE_RANGE + 1];
    } LvMapCoeffCost;

    /**************************************
     * The EbBitFraction is used to define the bit fraction numbers
     **************************************/
    typedef uint32_t EbBitFraction;

    /**************************************
     * MD Rate Estimation Structure
     **************************************/
    typedef struct MdRateEstimationContext
    {
        // Partition
        int32_t partition_fac_bits[PARTITION_CONTEXTS][EXT_PARTITION_TYPES];

        // MV Mode
        int32_t skip_mode_fac_bits[SKIP_CONTEXTS][2];
        int32_t new_mv_mode_fac_bits[NEWMV_MODE_CONTEXTS][2];
        int32_t zero_mv_mode_fac_bits[GLOBALMV_MODE_CONTEXTS][2];
        int32_t ref_mv_mode_fac_bits[REFMV_MODE_CONTEXTS][2];
        int32_t drl_mode_fac_bits[DRL_MODE_CONTEXTS][2];

        int32_t nmv_vec_cost[MV_JOINTS];
        int32_t nmv_costs[2][MV_VALS];
        int32_t nmv_costs_hp[2][MV_VALS];
        int32_t *nmvcoststack[2];
        int dv_cost[2][MV_VALS];
        int dv_joint_cost[MV_JOINTS];

        // Compouned Mode
        int32_t inter_compound_mode_fac_bits[INTER_MODE_CONTEXTS][INTER_COMPOUND_MODES];
        int32_t compound_type_fac_bits[BlockSizeS_ALL][MASKED_COMPOUND_TYPES];
        int32_t single_ref_fac_bits[REF_CONTEXTS][SINGLE_REFS - 1][2];
        int32_t comp_ref_type_fac_bits[COMP_REF_TYPE_CONTEXTS][2];
        int32_t uni_comp_ref_fac_bits[UNI_COMP_REF_CONTEXTS][UNIDIR_COMP_REFS - 1][2];
        int32_t comp_ref_fac_bits[REF_CONTEXTS][FWD_REFS - 1][2];
        int32_t comp_bwd_ref_fac_bits[REF_CONTEXTS][BWD_REFS - 1][2];
        int32_t comp_idx_fac_bits[COMP_INDEX_CONTEXTS][2];
        int32_t comp_group_idx_fac_bits[COMP_GROUP_IDX_CONTEXTS][2];
        int32_t comp_inter_fac_bits[COMP_INTER_CONTEXTS][2];

        // Wedge Mode
        int32_t wedge_idx_fac_bits[BlockSizeS_ALL][16];
        int32_t inter_intra_fac_bits[BlockSize_GROUPS][2];
        int32_t wedge_inter_intra_fac_bits[BlockSizeS_ALL][2];
        int32_t inter_intra_mode_fac_bits[BlockSize_GROUPS][INTERINTRA_MODES];
        int32_t motion_mode_fac_bits[BlockSizeS_ALL][MOTION_MODES];
        int32_t motion_mode_fac_bits1[BlockSizeS_ALL][2];

        // Intra Mode
        int32_t intrabc_fac_bits[2];
        int32_t intra_inter_fac_bits[INTRA_INTER_CONTEXTS][2];
        int32_t filter_intra_fac_bits[BlockSizeS_ALL][2];
        int32_t filter_intra_mode_fac_bits[FILTER_INTRA_MODES];
        int32_t switchable_restore_fac_bits[RESTORE_SWITCHABLE_TYPES];
        int32_t wiener_restore_fac_bits[2];
        int32_t sgrproj_restore_fac_bits[2];
        int32_t y_mode_fac_bits[KF_MODE_CONTEXTS][KF_MODE_CONTEXTS][INTRA_MODES];
        int32_t mb_mode_fac_bits[BlockSize_GROUPS][INTRA_MODES];
        int32_t intra_uv_mode_fac_bits[CFL_ALLOWED_TYPES][INTRA_MODES][UV_INTRA_MODES];
        int32_t angle_delta_fac_bits[DIRECTIONAL_MODES][2 * MAX_ANGLE_DELTA + 1];
        int32_t cfl_alpha_fac_bits[CFL_JOINT_SIGNS][CFL_PRED_PLANES][CFL_ALPHABET_SIZE];

        // Palette Mode
        int32_t palette_ysize_fac_bits[PALATTE_BSIZE_CTXS][PALETTE_SIZES];
        int32_t palette_uv_size_fac_bits[PALATTE_BSIZE_CTXS][PALETTE_SIZES];
        int32_t palette_ycolor_fac_bitss[PALETTE_SIZES][PALETTE_COLOR_INDEX_CONTEXTS][PALETTE_COLORS];
        int32_t palette_uv_color_fac_bits[PALETTE_SIZES][PALETTE_COLOR_INDEX_CONTEXTS][PALETTE_COLORS];
        int32_t palette_ymode_fac_bits[PALATTE_BSIZE_CTXS][PALETTE_Y_MODE_CONTEXTS][2];
        int32_t palette_uv_mode_fac_bits[PALETTE_UV_MODE_CONTEXTS][2];

        // Tx and Coeff Rate Estimation
        LvMapCoeffCost coeff_fac_bits[TX_SIZES][PLANE_TYPES];
        LvMapEobCost eob_frac_bits[7][2];
        int32_t txfm_partition_fac_bits[TXFM_PARTITION_CONTEXTS][2];
        int32_t skip_fac_bits[SKIP_CONTEXTS][2];
        int32_t tx_size_fac_bits[MAX_TX_CATS][TX_SIZE_CONTEXTS][MAX_TX_DEPTH + 1];
        int32_t intra_tx_type_fac_bits[EXT_TX_SETS_INTRA][EXT_TX_SIZES][INTRA_MODES][TX_TYPES];
        int32_t inter_tx_type_fac_bits[EXT_TX_SETS_INTER][EXT_TX_SIZES][TX_TYPES];
        int32_t switchable_interp_fac_bitss[SWITCHABLE_FILTER_CONTEXTS][SWITCHABLE_FILTERS];

        int32_t initialized;
    } MdRateEstimationContext;
    /***************************************************************************
    * AV1 Probability table
    * // round(-log2(i/256.) * (1 << AV1_PROB_COST_SHIFT)); i = 128~255.
    ***************************************************************************/
    static const uint16_t av1_prob_cost[128] = {
        512, 506, 501, 495, 489, 484, 478, 473, 467, 462, 456, 451, 446, 441, 435,
        430, 425, 420, 415, 410, 405, 400, 395, 390, 385, 380, 375, 371, 366, 361,
        356, 352, 347, 343, 338, 333, 329, 324, 320, 316, 311, 307, 302, 298, 294,
        289, 285, 281, 277, 273, 268, 264, 260, 256, 252, 248, 244, 240, 236, 232,
        228, 224, 220, 216, 212, 209, 205, 201, 197, 194, 190, 186, 182, 179, 175,
        171, 168, 164, 161, 157, 153, 150, 146, 143, 139, 136, 132, 129, 125, 122,
        119, 115, 112, 109, 105, 102, 99, 95, 92, 89, 86, 82, 79, 76, 73,
        70, 66, 63, 60, 57, 54, 51, 48, 45, 42, 38, 35, 32, 29, 26,
        23, 20, 18, 15, 12, 9, 6, 3,
    };
    static const int use_inter_ext_tx_for_txsize[EXT_TX_SETS_INTER]
        [EXT_TX_SIZES] = {
          { 1, 1, 1, 1 },  // unused
          { 1, 1, 0, 0 },
          { 0, 0, 1, 0 },
          { 0, 1, 1, 1 },
    };
    static const int32_t use_intra_ext_tx_for_txsize[EXT_TX_SETS_INTRA][EXT_TX_SIZES] =
    {
        { 1, 1, 1, 1 },  // unused
        { 1, 1, 0, 0 },
        { 0, 0, 1, 0 },
    };

    static const int32_t av1_ext_tx_set_idx_to_type[2][AOMMAX(EXT_TX_SETS_INTRA, EXT_TX_SETS_INTER)] = {
            {
            // Intra
            EXT_TX_SET_DCTONLY,
            EXT_TX_SET_DTT4_IDTX_1DDCT,
            EXT_TX_SET_DTT4_IDTX,
        },
        {
            // Inter
            EXT_TX_SET_DCTONLY,
            EXT_TX_SET_ALL16,
            EXT_TX_SET_DTT9_IDTX_1DDCT,
            EXT_TX_SET_DCT_IDTX,
        },
    };

    /***************************************************************************
    * av1_get_syntax_rate_from_cdf
    ***************************************************************************/
    extern void av1_get_syntax_rate_from_cdf(
        int32_t                        *costs,
        const AomCdfProb             *cdf,
        const int32_t                  *inv_map);
    /**************************************************************************
    * Estimate the rate for each syntax elements and for
    * all scenarios based on the frame CDF
    ***************************************************************************/
    extern void av1_estimate_syntax_rate(
        MdRateEstimationContext      *md_rate_estimation_array,
        EbBool                          is_i_slice,
        FRAME_CONTEXT                  *fc);
    /**************************************************************************
    * Estimate the rate of the quantised coefficient
    * based on the frame CDF
    ***************************************************************************/
    extern void av1_estimate_coefficients_rate(
        MdRateEstimationContext  *md_rate_estimation_array,
        FRAME_CONTEXT              *fc);
    /**************************************************************************
    * av1_estimate_mv_rate()
    * Estimate the rate of motion vectors
    * based on the frame CDF
    ***************************************************************************/
extern void av1_estimate_mv_rate(
        struct PictureControlSet *pcs_ptr,
        MdRateEstimationContext  *md_rate_estimation_array,
        FRAME_CONTEXT            *fc);
#define AVG_CDF_WEIGHT_LEFT      3
#define AVG_CDF_WEIGHT_TOP       1

static AOM_INLINE void avg_cdf_symbol(AomCdfProb *cdf_ptr_left,
    AomCdfProb *cdf_ptr_tr, int num_cdfs,
    int cdf_stride, int nsymbs, int wt_left,
    int wt_tr) {
    for (int i = 0; i < num_cdfs; i++) {
        for (int j = 0; j <= nsymbs; j++) {
            cdf_ptr_left[i * cdf_stride + j] =
                (AomCdfProb)(((int)cdf_ptr_left[i * cdf_stride + j] * wt_left +
                (int)cdf_ptr_tr[i * cdf_stride + j] * wt_tr +
                    ((wt_left + wt_tr) / 2)) /
                    (wt_left + wt_tr));
            assert(cdf_ptr_left[i * cdf_stride + j] >= 0 &&
                cdf_ptr_left[i * cdf_stride + j] < CDF_PROB_TOP);
        }
    }
}

#define AVERAGE_CDF(cname_left, cname_tr, nsymbs) \
  AVG_CDF_STRIDE(cname_left, cname_tr, nsymbs, CDF_SIZE(nsymbs))

#define AVG_CDF_STRIDE(cname_left, cname_tr, nsymbs, cdf_stride)           \
  do {                                                                     \
    AomCdfProb *cdf_ptr_left = (AomCdfProb *)cname_left;               \
    AomCdfProb *cdf_ptr_tr = (AomCdfProb *)cname_tr;                   \
    int array_size = (int)sizeof(cname_left) / sizeof(AomCdfProb);       \
    int num_cdfs = array_size / cdf_stride;                                \
    avg_cdf_symbol(cdf_ptr_left, cdf_ptr_tr, num_cdfs, cdf_stride, nsymbs, \
                   wt_left, wt_tr);                                        \
  } while (0)

static AOM_INLINE void avg_nmv(NmvContext *nmv_left, NmvContext *nmv_tr,
    int wt_left, int wt_tr) {
    AVERAGE_CDF(nmv_left->joints_cdf, nmv_tr->joints_cdf, 4);
    for (int i = 0; i < 2; i++) {
        AVERAGE_CDF(nmv_left->comps[i].classes_cdf, nmv_tr->comps[i].classes_cdf,
            MV_CLASSES);
        AVERAGE_CDF(nmv_left->comps[i].class0_fp_cdf,
            nmv_tr->comps[i].class0_fp_cdf, MV_FP_SIZE);
        AVERAGE_CDF(nmv_left->comps[i].fp_cdf, nmv_tr->comps[i].fp_cdf, MV_FP_SIZE);
        AVERAGE_CDF(nmv_left->comps[i].sign_cdf, nmv_tr->comps[i].sign_cdf, 2);
        AVERAGE_CDF(nmv_left->comps[i].class0_hp_cdf,
            nmv_tr->comps[i].class0_hp_cdf, 2);
        AVERAGE_CDF(nmv_left->comps[i].hp_cdf, nmv_tr->comps[i].hp_cdf, 2);
        AVERAGE_CDF(nmv_left->comps[i].class0_cdf, nmv_tr->comps[i].class0_cdf,
            CLASS0_SIZE);
        AVERAGE_CDF(nmv_left->comps[i].bits_cdf, nmv_tr->comps[i].bits_cdf, 2);
    }
}

// Since the top and left SBs are completed, we can average the top SB's CDFs and
// the left SB's CDFs and use it for current SB's encoding to
// improve the performance. This function facilitates the averaging
// of CDF.
static AOM_INLINE void avg_cdf_symbols(FRAME_CONTEXT *ctx_left,
    FRAME_CONTEXT *ctx_tr, int wt_left,
    int wt_tr) {
    AVERAGE_CDF(ctx_left->txb_skip_cdf, ctx_tr->txb_skip_cdf, 2);
    AVERAGE_CDF(ctx_left->eob_extra_cdf, ctx_tr->eob_extra_cdf, 2);
    AVERAGE_CDF(ctx_left->dc_sign_cdf, ctx_tr->dc_sign_cdf, 2);
    AVERAGE_CDF(ctx_left->eob_flag_cdf16, ctx_tr->eob_flag_cdf16, 5);
    AVERAGE_CDF(ctx_left->eob_flag_cdf32, ctx_tr->eob_flag_cdf32, 6);
    AVERAGE_CDF(ctx_left->eob_flag_cdf64, ctx_tr->eob_flag_cdf64, 7);
    AVERAGE_CDF(ctx_left->eob_flag_cdf128, ctx_tr->eob_flag_cdf128, 8);
    AVERAGE_CDF(ctx_left->eob_flag_cdf256, ctx_tr->eob_flag_cdf256, 9);
    AVERAGE_CDF(ctx_left->eob_flag_cdf512, ctx_tr->eob_flag_cdf512, 10);
    AVERAGE_CDF(ctx_left->eob_flag_cdf1024, ctx_tr->eob_flag_cdf1024, 11);
    AVERAGE_CDF(ctx_left->coeff_base_eob_cdf, ctx_tr->coeff_base_eob_cdf, 3);
    AVERAGE_CDF(ctx_left->coeff_base_cdf, ctx_tr->coeff_base_cdf, 4);
    AVERAGE_CDF(ctx_left->coeff_br_cdf, ctx_tr->coeff_br_cdf, BR_CDF_SIZE);
    AVERAGE_CDF(ctx_left->newmv_cdf, ctx_tr->newmv_cdf, 2);
    AVERAGE_CDF(ctx_left->zeromv_cdf, ctx_tr->zeromv_cdf, 2);
    AVERAGE_CDF(ctx_left->refmv_cdf, ctx_tr->refmv_cdf, 2);
    AVERAGE_CDF(ctx_left->drl_cdf, ctx_tr->drl_cdf, 2);
    AVERAGE_CDF(ctx_left->inter_compound_mode_cdf,
        ctx_tr->inter_compound_mode_cdf, INTER_COMPOUND_MODES);
    AVERAGE_CDF(ctx_left->compound_type_cdf, ctx_tr->compound_type_cdf,
        MASKED_COMPOUND_TYPES);
    AVERAGE_CDF(ctx_left->wedge_idx_cdf, ctx_tr->wedge_idx_cdf, 16);
    AVERAGE_CDF(ctx_left->interintra_cdf, ctx_tr->interintra_cdf, 2);
    AVERAGE_CDF(ctx_left->wedge_interintra_cdf, ctx_tr->wedge_interintra_cdf, 2);
    AVERAGE_CDF(ctx_left->interintra_mode_cdf, ctx_tr->interintra_mode_cdf,
        INTERINTRA_MODES);
    AVERAGE_CDF(ctx_left->motion_mode_cdf, ctx_tr->motion_mode_cdf, MOTION_MODES);
    AVERAGE_CDF(ctx_left->obmc_cdf, ctx_tr->obmc_cdf, 2);
    AVERAGE_CDF(ctx_left->palette_y_size_cdf, ctx_tr->palette_y_size_cdf,
        PALETTE_SIZES);
    AVERAGE_CDF(ctx_left->palette_uv_size_cdf, ctx_tr->palette_uv_size_cdf,
        PALETTE_SIZES);
    for (int j = 0; j < PALETTE_SIZES; j++) {
        int nsymbs = j + PALETTE_MIN_SIZE;
        AVG_CDF_STRIDE(ctx_left->palette_y_color_index_cdf[j],
            ctx_tr->palette_y_color_index_cdf[j], nsymbs,
            CDF_SIZE(PALETTE_COLORS));
        AVG_CDF_STRIDE(ctx_left->palette_uv_color_index_cdf[j],
            ctx_tr->palette_uv_color_index_cdf[j], nsymbs,
            CDF_SIZE(PALETTE_COLORS));
    }
    AVERAGE_CDF(ctx_left->palette_y_mode_cdf, ctx_tr->palette_y_mode_cdf, 2);
    AVERAGE_CDF(ctx_left->palette_uv_mode_cdf, ctx_tr->palette_uv_mode_cdf, 2);
    AVERAGE_CDF(ctx_left->comp_inter_cdf, ctx_tr->comp_inter_cdf, 2);
    AVERAGE_CDF(ctx_left->single_ref_cdf, ctx_tr->single_ref_cdf, 2);
    AVERAGE_CDF(ctx_left->comp_ref_type_cdf, ctx_tr->comp_ref_type_cdf, 2);
    AVERAGE_CDF(ctx_left->uni_comp_ref_cdf, ctx_tr->uni_comp_ref_cdf, 2);
    AVERAGE_CDF(ctx_left->comp_ref_cdf, ctx_tr->comp_ref_cdf, 2);
    AVERAGE_CDF(ctx_left->comp_bwdref_cdf, ctx_tr->comp_bwdref_cdf, 2);
    AVERAGE_CDF(ctx_left->txfm_partition_cdf, ctx_tr->txfm_partition_cdf, 2);
    AVERAGE_CDF(ctx_left->compound_index_cdf, ctx_tr->compound_index_cdf, 2);
    AVERAGE_CDF(ctx_left->comp_group_idx_cdf, ctx_tr->comp_group_idx_cdf, 2);
    AVERAGE_CDF(ctx_left->skip_mode_cdfs, ctx_tr->skip_mode_cdfs, 2);
    AVERAGE_CDF(ctx_left->skip_cdfs, ctx_tr->skip_cdfs, 2);
    AVERAGE_CDF(ctx_left->intra_inter_cdf, ctx_tr->intra_inter_cdf, 2);
    avg_nmv(&ctx_left->nmvc, &ctx_tr->nmvc, wt_left, wt_tr);
    avg_nmv(&ctx_left->ndvc, &ctx_tr->ndvc, wt_left, wt_tr);
    AVERAGE_CDF(ctx_left->intrabc_cdf, ctx_tr->intrabc_cdf, 2);
    AVERAGE_CDF(ctx_left->seg.tree_cdf, ctx_tr->seg.tree_cdf, MAX_SEGMENTS);
    AVERAGE_CDF(ctx_left->seg.pred_cdf, ctx_tr->seg.pred_cdf, 2);
    AVERAGE_CDF(ctx_left->seg.spatial_pred_seg_cdf,
        ctx_tr->seg.spatial_pred_seg_cdf, MAX_SEGMENTS);
    AVERAGE_CDF(ctx_left->filter_intra_cdfs, ctx_tr->filter_intra_cdfs, 2);
    AVERAGE_CDF(ctx_left->filter_intra_mode_cdf, ctx_tr->filter_intra_mode_cdf,
        FILTER_INTRA_MODES);
    AVERAGE_CDF(ctx_left->switchable_restore_cdf, ctx_tr->switchable_restore_cdf,
        RESTORE_SWITCHABLE_TYPES);
    AVERAGE_CDF(ctx_left->wiener_restore_cdf, ctx_tr->wiener_restore_cdf, 2);
    AVERAGE_CDF(ctx_left->sgrproj_restore_cdf, ctx_tr->sgrproj_restore_cdf, 2);
    AVERAGE_CDF(ctx_left->y_mode_cdf, ctx_tr->y_mode_cdf, INTRA_MODES);
    AVG_CDF_STRIDE(ctx_left->uv_mode_cdf[0], ctx_tr->uv_mode_cdf[0],
        UV_INTRA_MODES - 1, CDF_SIZE(UV_INTRA_MODES));
    AVERAGE_CDF(ctx_left->uv_mode_cdf[1], ctx_tr->uv_mode_cdf[1], UV_INTRA_MODES);
    for (int i = 0; i < PARTITION_CONTEXTS; i++) {
        if (i < 4) {
            AVG_CDF_STRIDE(ctx_left->partition_cdf[i], ctx_tr->partition_cdf[i], 4,
                CDF_SIZE(10));
        }
        else if (i < 16) {
            AVERAGE_CDF(ctx_left->partition_cdf[i], ctx_tr->partition_cdf[i], 10);
        }
        else {
            AVG_CDF_STRIDE(ctx_left->partition_cdf[i], ctx_tr->partition_cdf[i], 8,
                CDF_SIZE(10));
        }
    }
    AVERAGE_CDF(ctx_left->switchable_interp_cdf, ctx_tr->switchable_interp_cdf,
        SWITCHABLE_FILTERS);
    AVERAGE_CDF(ctx_left->kf_y_cdf, ctx_tr->kf_y_cdf, INTRA_MODES);
    AVERAGE_CDF(ctx_left->angle_delta_cdf, ctx_tr->angle_delta_cdf,
        2 * MAX_ANGLE_DELTA + 1);
    AVG_CDF_STRIDE(ctx_left->tx_size_cdf[0], ctx_tr->tx_size_cdf[0], MAX_TX_DEPTH,
        CDF_SIZE(MAX_TX_DEPTH + 1));
    AVERAGE_CDF(ctx_left->tx_size_cdf[1], ctx_tr->tx_size_cdf[1],
        MAX_TX_DEPTH + 1);
    AVERAGE_CDF(ctx_left->tx_size_cdf[2], ctx_tr->tx_size_cdf[2],
        MAX_TX_DEPTH + 1);
    AVERAGE_CDF(ctx_left->tx_size_cdf[3], ctx_tr->tx_size_cdf[3],
        MAX_TX_DEPTH + 1);
    AVERAGE_CDF(ctx_left->delta_q_cdf, ctx_tr->delta_q_cdf, DELTA_Q_PROBS + 1);
    AVERAGE_CDF(ctx_left->delta_lf_cdf, ctx_tr->delta_lf_cdf, DELTA_LF_PROBS + 1);
    for (int i = 0; i < FRAME_LF_COUNT; i++) {
        AVERAGE_CDF(ctx_left->delta_lf_multi_cdf[i], ctx_tr->delta_lf_multi_cdf[i],
            DELTA_LF_PROBS + 1);
    }
    AVG_CDF_STRIDE(ctx_left->intra_ext_tx_cdf[1], ctx_tr->intra_ext_tx_cdf[1], 7,
        CDF_SIZE(TX_TYPES));
    AVG_CDF_STRIDE(ctx_left->intra_ext_tx_cdf[2], ctx_tr->intra_ext_tx_cdf[2], 5,
        CDF_SIZE(TX_TYPES));
    AVG_CDF_STRIDE(ctx_left->inter_ext_tx_cdf[1], ctx_tr->inter_ext_tx_cdf[1], 16,
        CDF_SIZE(TX_TYPES));
    AVG_CDF_STRIDE(ctx_left->inter_ext_tx_cdf[2], ctx_tr->inter_ext_tx_cdf[2], 12,
        CDF_SIZE(TX_TYPES));
    AVG_CDF_STRIDE(ctx_left->inter_ext_tx_cdf[3], ctx_tr->inter_ext_tx_cdf[3], 2,
        CDF_SIZE(TX_TYPES));
    AVERAGE_CDF(ctx_left->cfl_sign_cdf, ctx_tr->cfl_sign_cdf, CFL_JOINT_SIGNS);
    AVERAGE_CDF(ctx_left->cfl_alpha_cdf, ctx_tr->cfl_alpha_cdf,
        CFL_ALPHABET_SIZE);
}
/*******************************************************************************
 * Updates all the syntax stats/CDF for the current block
 ******************************************************************************/
void update_stats(
    struct PictureControlSet   *pcs_ptr,
    struct BlkStruct          *blk_ptr,
    int                         mi_row,
    int                         mi_col);
/*******************************************************************************
 * Updates the partition stats/CDF for the current block
 ******************************************************************************/
void update_part_stats(
    struct PictureControlSet   *pcs_ptr,
    struct BlkStruct          *blk_ptr,
    uint16_t                    tile_idx,
    int                         mi_row,
    int                         mi_col);

#ifdef __cplusplus
}
#endif

#endif //EbMdRateEstimationTables_h
// clang-format on
