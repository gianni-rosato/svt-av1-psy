/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbMdRateEstimation_h
#define EbMdRateEstimation_h

#include "EbDefinitions.h"
#include "EbCabacContextModel.h"
#include "EbPictureControlSet.h"
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
        EbBitFraction  split_flag_bits[NUMBER_OF_SPLIT_FLAG_CASES];
        EbBitFraction  mvd_bits[NUMBER_OF_MVD_CASES];
        // Partition
        int32_t partition_fac_bits[PARTITION_CONTEXTS][CDF_SIZE(EXT_PARTITION_TYPES)];

        // MV Mode
        int32_t skip_mode_fac_bits[SKIP_CONTEXTS][CDF_SIZE(2)];
        int32_t new_mv_mode_fac_bits[NEWMV_MODE_CONTEXTS][CDF_SIZE(2)];
        int32_t zero_mv_mode_fac_bits[GLOBALMV_MODE_CONTEXTS][CDF_SIZE(2)];
        int32_t ref_mv_mode_fac_bits[REFMV_MODE_CONTEXTS][CDF_SIZE(2)];
        int32_t drl_mode_fac_bits[DRL_MODE_CONTEXTS][CDF_SIZE(2)];

        int32_t nmv_vec_cost[MV_JOINTS];
        int32_t nmv_costs[2][MV_VALS];
        int32_t nmv_costs_hp[2][MV_VALS];
        int32_t *nmvcoststack[2];
        int dv_cost[2][MV_VALS];
        int dv_joint_cost[MV_JOINTS];

        // Compouned Mode
        int32_t inter_compound_mode_fac_bits[INTER_MODE_CONTEXTS][CDF_SIZE(INTER_COMPOUND_MODES)];
        int32_t compound_type_fac_bits[BlockSizeS_ALL][CDF_SIZE(MASKED_COMPOUND_TYPES)];
        int32_t single_ref_fac_bits[REF_CONTEXTS][SINGLE_REFS - 1][CDF_SIZE(2)];
        int32_t comp_ref_type_fac_bits[COMP_REF_TYPE_CONTEXTS][CDF_SIZE(2)];
        int32_t uni_comp_ref_fac_bits[UNI_COMP_REF_CONTEXTS][UNIDIR_COMP_REFS - 1][CDF_SIZE(2)];
        int32_t comp_ref_fac_bits[REF_CONTEXTS][FWD_REFS - 1][CDF_SIZE(2)];
        int32_t comp_bwd_ref_fac_bits[REF_CONTEXTS][BWD_REFS - 1][CDF_SIZE(2)];
        int32_t comp_idx_fac_bits[COMP_INDEX_CONTEXTS][CDF_SIZE(2)];
        int32_t comp_group_idx_fac_bits[COMP_GROUP_IDX_CONTEXTS][CDF_SIZE(2)];
        int32_t comp_inter_fac_bits[COMP_INTER_CONTEXTS][CDF_SIZE(2)];

        // Wedge Mode
        int32_t wedge_idx_fac_bits[BlockSizeS_ALL][CDF_SIZE(16)];
        int32_t inter_intra_fac_bits[BlockSize_GROUPS][CDF_SIZE(2)];
        int32_t wedge_inter_intra_fac_bits[BlockSizeS_ALL][CDF_SIZE(2)];
        int32_t inter_intra_mode_fac_bits[BlockSize_GROUPS][CDF_SIZE(INTERINTRA_MODES)];
        int32_t motion_mode_fac_bits[BlockSizeS_ALL][MOTION_MODES];
        int32_t motion_mode_fac_bits1[BlockSizeS_ALL][2];

        // Intra Mode
        int32_t intrabc_fac_bits[CDF_SIZE(2)];
        int32_t intra_inter_fac_bits[INTRA_INTER_CONTEXTS][2];
        int32_t filter_intra_fac_bits[BlockSizeS_ALL][CDF_SIZE(2)];
        int32_t filter_intra_mode_fac_bits[CDF_SIZE(FILTER_INTRA_MODES)];
        int32_t switchable_restore_fac_bits[CDF_SIZE(RESTORE_SWITCHABLE_TYPES)];
        int32_t wiener_restore_fac_bits[CDF_SIZE(2)];
        int32_t sgrproj_restore_fac_bits[CDF_SIZE(2)];
        int32_t y_mode_fac_bits[KF_MODE_CONTEXTS][KF_MODE_CONTEXTS][CDF_SIZE(INTRA_MODES)];
        int32_t mb_mode_fac_bits[BlockSize_GROUPS][CDF_SIZE(INTRA_MODES)];
        int32_t intra_uv_mode_fac_bits[CFL_ALLOWED_TYPES][INTRA_MODES][CDF_SIZE(UV_INTRA_MODES)];
        int32_t angle_delta_fac_bits[DIRECTIONAL_MODES][CDF_SIZE(2 * MAX_ANGLE_DELTA + 1)];
        int32_t cfl_alpha_fac_bits[CFL_JOINT_SIGNS][CFL_PRED_PLANES][CFL_ALPHABET_SIZE];

        // Palette Mode
        int32_t palette_ysize_fac_bits[PALATTE_BSIZE_CTXS][CDF_SIZE(PALETTE_SIZES)];
        int32_t palette_uv_size_fac_bits[PALATTE_BSIZE_CTXS][CDF_SIZE(PALETTE_SIZES)];
#if PAL_SUP
        int32_t palette_ycolor_fac_bitss[PALETTE_SIZES][PALETTE_COLOR_INDEX_CONTEXTS][PALETTE_COLORS];
        int32_t palette_uv_color_fac_bits[PALETTE_SIZES][PALETTE_COLOR_INDEX_CONTEXTS][PALETTE_COLORS];
#else
        int32_t palette_ycolor_fac_bitss[PALETTE_SIZES][PALETTE_COLOR_INDEX_CONTEXTS][CDF_SIZE(PALETTE_COLORS)];
        int32_t palette_uv_color_fac_bits[PALETTE_SIZES][PALETTE_COLOR_INDEX_CONTEXTS][CDF_SIZE(PALETTE_COLORS)];
#endif
        int32_t palette_ymode_fac_bits[PALATTE_BSIZE_CTXS][PALETTE_Y_MODE_CONTEXTS][CDF_SIZE(2)];
        int32_t palette_uv_mode_fac_bits[PALETTE_UV_MODE_CONTEXTS][CDF_SIZE(2)];

        // Tx and Coeff Rate Estimation
        LvMapCoeffCost coeff_fac_bits[TX_SIZES][PLANE_TYPES];
        LvMapEobCost eob_frac_bits[7][2];
        int32_t txfm_partition_fac_bits[TXFM_PARTITION_CONTEXTS][CDF_SIZE(2)];
        int32_t skip_fac_bits[SKIP_CONTEXTS][CDF_SIZE(2)];
        int32_t tx_size_fac_bits[MAX_TX_CATS][TX_SIZE_CONTEXTS][CDF_SIZE(MAX_TX_DEPTH + 1)];
        int32_t intra_tx_type_fac_bits[EXT_TX_SETS_INTRA][EXT_TX_SIZES][INTRA_MODES][CDF_SIZE(TX_TYPES)];
        int32_t inter_tx_type_fac_bits[EXT_TX_SETS_INTER][EXT_TX_SIZES][CDF_SIZE(TX_TYPES)];
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

    static const int32_t use_inter_ext_tx_for_txsize[EXT_TX_SETS_INTER][EXT_TX_SIZES] =
    {
        { 1, 1, 1, 1 },  // unused
        { 1, 1, 0, 0 },
        { 0, 0, 1, 0 },
        { 0, 0, 0, 1 },
    };

    static const int32_t use_intra_ext_tx_for_txsize[EXT_TX_SETS_INTRA][EXT_TX_SIZES] =
    {
        { 1, 1, 1, 1 },  // unused
        { 1, 1, 0, 0 },
        { 0, 0, 1, 0 },
    };

    DECLARE_ALIGNED(16, static uint8_t, wedge_signflip_lookup[BlockSizeS_ALL][MAX_WEDGE_TYPES]) = {
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },  // not used
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },  // not used
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },  // not used
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, },
        { 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, },
        { 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, },
        { 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, },
        { 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },  // not used
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },  // not used
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },  // not used
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },  // not used
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },  // not used
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },  // not used
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },  // not used
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },  // not used
        { 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, },
        { 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },  // not used
        { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },  // not used
    };

    typedef uint8_t *WedgeMasksType[MAX_WEDGE_TYPES];
    static WedgeMasksType wedge_masks[BlockSizeS_ALL][2];
    // Angles are with respect to horizontal anti-clockwise
    typedef enum WedgeDirectionType
    {
        WEDGE_HORIZONTAL = 0,
        WEDGE_VERTICAL = 1,
        WEDGE_OBLIQUE27 = 2,
        WEDGE_OBLIQUE63 = 3,
        WEDGE_OBLIQUE117 = 4,
        WEDGE_OBLIQUE153 = 5,
        WEDGE_DIRECTIONS
    } WedgeDirectionType;
    // 3-tuple: {direction, x_offset, y_offset}
    typedef struct WedgeCodeType
    {
        WedgeDirectionType direction;
        int32_t x_offset;
        int32_t y_offset;
    } WedgeCodeType;

    typedef struct WedgeParamsType
    {
        int32_t bits;
        const WedgeCodeType *codebook;
        uint8_t *signflip;
        WedgeMasksType *masks;
    } WedgeParamsType;

    static const WedgeCodeType wedge_codebook_16_hgtw[16] = {
        { WEDGE_OBLIQUE27, 4, 4 }, { WEDGE_OBLIQUE63, 4, 4 },
        { WEDGE_OBLIQUE117, 4, 4 }, { WEDGE_OBLIQUE153, 4, 4 },
        { WEDGE_HORIZONTAL, 4, 2 }, { WEDGE_HORIZONTAL, 4, 4 },
        { WEDGE_HORIZONTAL, 4, 6 }, { WEDGE_VERTICAL, 4, 4 },
        { WEDGE_OBLIQUE27, 4, 2 }, { WEDGE_OBLIQUE27, 4, 6 },
        { WEDGE_OBLIQUE153, 4, 2 }, { WEDGE_OBLIQUE153, 4, 6 },
        { WEDGE_OBLIQUE63, 2, 4 }, { WEDGE_OBLIQUE63, 6, 4 },
        { WEDGE_OBLIQUE117, 2, 4 }, { WEDGE_OBLIQUE117, 6, 4 },
    };

    static const WedgeCodeType wedge_codebook_16_hltw[16] = {
        { WEDGE_OBLIQUE27, 4, 4 }, { WEDGE_OBLIQUE63, 4, 4 },
        { WEDGE_OBLIQUE117, 4, 4 }, { WEDGE_OBLIQUE153, 4, 4 },
        { WEDGE_VERTICAL, 2, 4 }, { WEDGE_VERTICAL, 4, 4 },
        { WEDGE_VERTICAL, 6, 4 }, { WEDGE_HORIZONTAL, 4, 4 },
        { WEDGE_OBLIQUE27, 4, 2 }, { WEDGE_OBLIQUE27, 4, 6 },
        { WEDGE_OBLIQUE153, 4, 2 }, { WEDGE_OBLIQUE153, 4, 6 },
        { WEDGE_OBLIQUE63, 2, 4 }, { WEDGE_OBLIQUE63, 6, 4 },
        { WEDGE_OBLIQUE117, 2, 4 }, { WEDGE_OBLIQUE117, 6, 4 },
    };

    static const WedgeCodeType wedge_codebook_16_heqw[16] = {
        { WEDGE_OBLIQUE27, 4, 4 }, { WEDGE_OBLIQUE63, 4, 4 },
        { WEDGE_OBLIQUE117, 4, 4 }, { WEDGE_OBLIQUE153, 4, 4 },
        { WEDGE_HORIZONTAL, 4, 2 }, { WEDGE_HORIZONTAL, 4, 6 },
        { WEDGE_VERTICAL, 2, 4 }, { WEDGE_VERTICAL, 6, 4 },
        { WEDGE_OBLIQUE27, 4, 2 }, { WEDGE_OBLIQUE27, 4, 6 },
        { WEDGE_OBLIQUE153, 4, 2 }, { WEDGE_OBLIQUE153, 4, 6 },
        { WEDGE_OBLIQUE63, 2, 4 }, { WEDGE_OBLIQUE63, 6, 4 },
        { WEDGE_OBLIQUE117, 2, 4 }, { WEDGE_OBLIQUE117, 6, 4 },
    };

    static const WedgeParamsType wedge_params_lookup[BlockSizeS_ALL] = {
        { 0, NULL, NULL, NULL },
        { 0, NULL, NULL, NULL },
        { 0, NULL, NULL, NULL },
        { 4, wedge_codebook_16_heqw, wedge_signflip_lookup[BLOCK_8X8],
        wedge_masks[BLOCK_8X8] },
        { 4, wedge_codebook_16_hgtw, wedge_signflip_lookup[BLOCK_8X16],
        wedge_masks[BLOCK_8X16] },
        { 4, wedge_codebook_16_hltw, wedge_signflip_lookup[BLOCK_16X8],
        wedge_masks[BLOCK_16X8] },
        { 4, wedge_codebook_16_heqw, wedge_signflip_lookup[BLOCK_16X16],
        wedge_masks[BLOCK_16X16] },
        { 4, wedge_codebook_16_hgtw, wedge_signflip_lookup[BLOCK_16X32],
        wedge_masks[BLOCK_16X32] },
        { 4, wedge_codebook_16_hltw, wedge_signflip_lookup[BLOCK_32X16],
        wedge_masks[BLOCK_32X16] },
        { 4, wedge_codebook_16_heqw, wedge_signflip_lookup[BLOCK_32X32],
        wedge_masks[BLOCK_32X32] },
        { 0, NULL, NULL, NULL },
        { 0, NULL, NULL, NULL },
        { 0, NULL, NULL, NULL },
        { 0, NULL, NULL, NULL },
        { 0, NULL, NULL, NULL },
        { 0, NULL, NULL, NULL },
        { 0, NULL, NULL, NULL },
        { 0, NULL, NULL, NULL },
        { 4, wedge_codebook_16_hgtw, wedge_signflip_lookup[BLOCK_8X32],
        wedge_masks[BLOCK_8X32] },
        { 4, wedge_codebook_16_hltw, wedge_signflip_lookup[BLOCK_32X8],
        wedge_masks[BLOCK_32X8] },
        { 0, NULL, NULL, NULL },
        { 0, NULL, NULL, NULL },
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
        struct PictureControlSet     *picture_control_set_ptr,
        MdRateEstimationContext  *md_rate_estimation_array,
        NmvContext                *nmv_ctx);

#ifdef __cplusplus
}
#endif

#endif //EbMdRateEstimationTables_h
