/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbMdRateEstimation_h
#define EbMdRateEstimation_h

#include "EbDefinitions.h"
#include "EbCabacContextModel.h"
#if ICOPY
#include "EbPictureControlSet.h"
#endif
#ifdef __cplusplus
extern "C" {
#endif
    /**************************************
     * MD Rate Estimation Defines
     **************************************/
#if ICOPY
#define MV_COST_WEIGHT_SUB 120
#endif

#define TOTAL_NUMBER_OF_MD_RATE_ESTIMATION_CASE_BUFFERS (TOTAL_NUMBER_OF_QP_VALUES * TOTAL_NUMBER_OF_SLICE_TYPES)
#define NUMBER_OF_SPLIT_FLAG_CASES                            6       // number of cases for bit estimation for split flag
#define NUMBER_OF_MVD_CASES                                  12       // number of cases for bit estimation for motion vector difference
     // Set to (1 << 5) if the 32-ary codebooks are used for any bock size
#define MAX_WEDGE_TYPES                                      (1 << 4)
     // The factor to scale from cost in bits to cost in av1_prob_cost units.
#define AV1_PROB_COST_SHIFT 9
     // Cost of coding an n bit literal, using 128 (i.e. 50%) probability for each bit.
#define av1_cost_literal(n) ((n) * (1 << AV1_PROB_COST_SHIFT))

    typedef struct {
        int32_t eob_cost[2][11];
    } LV_MAP_EOB_COST;

    typedef struct {
        int32_t txb_skip_cost[TXB_SKIP_CONTEXTS][2];
        int32_t base_eob_cost[SIG_COEF_CONTEXTS_EOB][3];
        int32_t base_cost[SIG_COEF_CONTEXTS][4];
        int32_t eob_extra_cost[EOB_COEF_CONTEXTS][2];
        int32_t dc_sign_cost[DC_SIGN_CONTEXTS][2];
        int32_t lps_cost[LEVEL_CONTEXTS][COEFF_BASE_RANGE + 1];
    } LV_MAP_COEFF_COST;


    /**************************************
     * The EB_BitFraction is used to define the bit fraction numbers
     **************************************/
    typedef uint32_t EB_BitFraction;

    /**************************************
     * MD Rate Estimation Structure
     **************************************/
    typedef struct MdRateEstimationContext_s {
        EB_BitFraction  splitFlagBits[NUMBER_OF_SPLIT_FLAG_CASES];
        EB_BitFraction  mvdBits[NUMBER_OF_MVD_CASES];
        // Partition
        int32_t partitionFacBits[PARTITION_CONTEXTS][CDF_SIZE(EXT_PARTITION_TYPES)];

        // MV Mode
        int32_t skipModeFacBits[SKIP_CONTEXTS][CDF_SIZE(2)];
        int32_t newMvModeFacBits[NEWMV_MODE_CONTEXTS][CDF_SIZE(2)];
        int32_t zeroMvModeFacBits[GLOBALMV_MODE_CONTEXTS][CDF_SIZE(2)];
        int32_t refMvModeFacBits[REFMV_MODE_CONTEXTS][CDF_SIZE(2)];
        int32_t drlModeFacBits[DRL_MODE_CONTEXTS][CDF_SIZE(2)];

        int32_t nmv_vec_cost[MV_JOINTS];
        int32_t nmv_costs[2][MV_VALS];
        int32_t nmv_costs_hp[2][MV_VALS];
        int32_t *nmvcoststack[2];
#if ICOPY
        int dv_cost[2][MV_VALS];
        int dv_joint_cost[MV_JOINTS];
#endif

        // Compouned Mode
        int32_t interCompoundModeFacBits[INTER_MODE_CONTEXTS][CDF_SIZE(INTER_COMPOUND_MODES)];
        int32_t compoundTypeFacBits[BlockSizeS_ALL][CDF_SIZE(COMPOUND_TYPES - 1)];
        int32_t singleRefFacBits[REF_CONTEXTS][SINGLE_REFS - 1][CDF_SIZE(2)];
        int32_t compRefTypeFacBits[COMP_REF_TYPE_CONTEXTS][CDF_SIZE(2)];
        int32_t uniCompRefFacBits[UNI_COMP_REF_CONTEXTS][UNIDIR_COMP_REFS - 1][CDF_SIZE(2)];
        int32_t compRefFacBits[REF_CONTEXTS][FWD_REFS - 1][CDF_SIZE(2)];
        int32_t compBwdRefFacBits[REF_CONTEXTS][BWD_REFS - 1][CDF_SIZE(2)];
        int32_t compIdxFacBits[COMP_INDEX_CONTEXTS][CDF_SIZE(2)];
        int32_t compGroupIdxFacBits[COMP_GROUP_IDX_CONTEXTS][CDF_SIZE(2)];
        int32_t compInterFacBits[COMP_INTER_CONTEXTS][CDF_SIZE(2)];

        // Wedge Mode
        int32_t wedgeIdxFacBits[BlockSizeS_ALL][CDF_SIZE(16)];
        int32_t interIntraFacBits[BlockSize_GROUPS][CDF_SIZE(2)];
        int32_t wedgeInterIntraFacBits[BlockSizeS_ALL][CDF_SIZE(2)];
        int32_t interIntraModeFacBits[BlockSize_GROUPS][CDF_SIZE(INTERINTRA_MODES)];
        int32_t motionModeFacBits[BlockSizeS_ALL][MOTION_MODES];
        int32_t motionModeFacBits1[BlockSizeS_ALL][2];

        // Intra Mode
        int32_t intrabcFacBits[CDF_SIZE(2)];
        int32_t intraInterFacBits[INTRA_INTER_CONTEXTS][2];
        int32_t filter_intra_FacBits[BlockSizeS_ALL][CDF_SIZE(2)];
        int32_t filterIntraModeFacBits[CDF_SIZE(FILTER_INTRA_MODES)];
        int32_t switchableRestoreFacBits[CDF_SIZE(RESTORE_SWITCHABLE_TYPES)];
        int32_t wienerRestoreFacBits[CDF_SIZE(2)];
        int32_t sgrprojRestoreFacBits[CDF_SIZE(2)];
        int32_t yModeFacBits[KF_MODE_CONTEXTS][KF_MODE_CONTEXTS][CDF_SIZE(INTRA_MODES)];
        int32_t mbModeFacBits[BlockSize_GROUPS][CDF_SIZE(INTRA_MODES)];
        int32_t intraUVmodeFacBits[CFL_ALLOWED_TYPES][INTRA_MODES][CDF_SIZE(UV_INTRA_MODES)];
        int32_t angleDeltaFacBits[DIRECTIONAL_MODES][CDF_SIZE(2 * MAX_ANGLE_DELTA + 1)];
        int32_t cflAlphaFacBits[CFL_JOINT_SIGNS][CFL_PRED_PLANES][CFL_ALPHABET_SIZE];

        // Palette Mode
        int32_t paletteYsizeFacBits[PALATTE_BSIZE_CTXS][CDF_SIZE(PALETTE_SIZES)];
        int32_t paletteUVsizeFacBits[PALATTE_BSIZE_CTXS][CDF_SIZE(PALETTE_SIZES)];
        int32_t paletteYcolorFacBitss[PALETTE_SIZES][PALETTE_COLOR_INDEX_CONTEXTS][CDF_SIZE(PALETTE_COLORS)];
        int32_t paletteUVcolorFacBits[PALETTE_SIZES][PALETTE_COLOR_INDEX_CONTEXTS][CDF_SIZE(PALETTE_COLORS)];
        int32_t paletteYmodeFacBits[PALATTE_BSIZE_CTXS][PALETTE_Y_MODE_CONTEXTS][CDF_SIZE(2)];
        int32_t paletteUVmodeFacBits[PALETTE_UV_MODE_CONTEXTS][CDF_SIZE(2)];

        // Tx and Coeff Rate Estimation
        LV_MAP_COEFF_COST coeffFacBits[TX_SIZES][PLANE_TYPES];
        LV_MAP_EOB_COST eobFracBits[7][2];
        int32_t txfmPartitionFacBits[TXFM_PARTITION_CONTEXTS][CDF_SIZE(2)];
        int32_t skipFacBits[SKIP_CONTEXTS][CDF_SIZE(2)];
        int32_t txSizeFacBits[MAX_TX_CATS][TX_SIZE_CONTEXTS][CDF_SIZE(MAX_TX_DEPTH + 1)];
        int32_t intraTxTypeFacBits[EXT_TX_SETS_INTRA][EXT_TX_SIZES][INTRA_MODES][CDF_SIZE(TX_TYPES)];
        int32_t interTxTypeFacBits[EXT_TX_SETS_INTER][EXT_TX_SIZES][CDF_SIZE(TX_TYPES)];
        int32_t switchable_interp_FacBitss[SWITCHABLE_FILTER_CONTEXTS][SWITCHABLE_FILTERS];
        int32_t initialized;

    } MdRateEstimationContext_t;

    /**************************************
    * Extern Function Declarations
    **************************************/
    extern EbErrorType MdRateEstimationContextCtor(MdRateEstimationContext_t *md_rate_estimation_array);
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

    typedef uint8_t *wedge_masks_type[MAX_WEDGE_TYPES];
    static wedge_masks_type wedge_masks[BlockSizeS_ALL][2];
    // Angles are with respect to horizontal anti-clockwise
    typedef enum {
        WEDGE_HORIZONTAL = 0,
        WEDGE_VERTICAL = 1,
        WEDGE_OBLIQUE27 = 2,
        WEDGE_OBLIQUE63 = 3,
        WEDGE_OBLIQUE117 = 4,
        WEDGE_OBLIQUE153 = 5,
        WEDGE_DIRECTIONS
    } WedgeDirectionType;
    // 3-tuple: {direction, x_offset, y_offset}
    typedef struct {
        WedgeDirectionType direction;
        int32_t x_offset;
        int32_t y_offset;
    } wedge_code_type;

    typedef struct {
        int32_t bits;
        const wedge_code_type *codebook;
        uint8_t *signflip;
        wedge_masks_type *masks;
    } wedge_params_type;

    static const wedge_code_type wedge_codebook_16_hgtw[16] = {
        { WEDGE_OBLIQUE27, 4, 4 }, { WEDGE_OBLIQUE63, 4, 4 },
        { WEDGE_OBLIQUE117, 4, 4 }, { WEDGE_OBLIQUE153, 4, 4 },
        { WEDGE_HORIZONTAL, 4, 2 }, { WEDGE_HORIZONTAL, 4, 4 },
        { WEDGE_HORIZONTAL, 4, 6 }, { WEDGE_VERTICAL, 4, 4 },
        { WEDGE_OBLIQUE27, 4, 2 }, { WEDGE_OBLIQUE27, 4, 6 },
        { WEDGE_OBLIQUE153, 4, 2 }, { WEDGE_OBLIQUE153, 4, 6 },
        { WEDGE_OBLIQUE63, 2, 4 }, { WEDGE_OBLIQUE63, 6, 4 },
        { WEDGE_OBLIQUE117, 2, 4 }, { WEDGE_OBLIQUE117, 6, 4 },
    };

    static const wedge_code_type wedge_codebook_16_hltw[16] = {
        { WEDGE_OBLIQUE27, 4, 4 }, { WEDGE_OBLIQUE63, 4, 4 },
        { WEDGE_OBLIQUE117, 4, 4 }, { WEDGE_OBLIQUE153, 4, 4 },
        { WEDGE_VERTICAL, 2, 4 }, { WEDGE_VERTICAL, 4, 4 },
        { WEDGE_VERTICAL, 6, 4 }, { WEDGE_HORIZONTAL, 4, 4 },
        { WEDGE_OBLIQUE27, 4, 2 }, { WEDGE_OBLIQUE27, 4, 6 },
        { WEDGE_OBLIQUE153, 4, 2 }, { WEDGE_OBLIQUE153, 4, 6 },
        { WEDGE_OBLIQUE63, 2, 4 }, { WEDGE_OBLIQUE63, 6, 4 },
        { WEDGE_OBLIQUE117, 2, 4 }, { WEDGE_OBLIQUE117, 6, 4 },
    };

    static const wedge_code_type wedge_codebook_16_heqw[16] = {
        { WEDGE_OBLIQUE27, 4, 4 }, { WEDGE_OBLIQUE63, 4, 4 },
        { WEDGE_OBLIQUE117, 4, 4 }, { WEDGE_OBLIQUE153, 4, 4 },
        { WEDGE_HORIZONTAL, 4, 2 }, { WEDGE_HORIZONTAL, 4, 6 },
        { WEDGE_VERTICAL, 2, 4 }, { WEDGE_VERTICAL, 6, 4 },
        { WEDGE_OBLIQUE27, 4, 2 }, { WEDGE_OBLIQUE27, 4, 6 },
        { WEDGE_OBLIQUE153, 4, 2 }, { WEDGE_OBLIQUE153, 4, 6 },
        { WEDGE_OBLIQUE63, 2, 4 }, { WEDGE_OBLIQUE63, 6, 4 },
        { WEDGE_OBLIQUE117, 2, 4 }, { WEDGE_OBLIQUE117, 6, 4 },
    };

    static const wedge_params_type wedge_params_lookup[BlockSizeS_ALL] = {
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
        const aom_cdf_prob             *cdf,
        const int32_t                  *inv_map);
    /**************************************************************************
    * Estimate the rate for each syntax elements and for
    * all scenarios based on the frame CDF
    ***************************************************************************/
    extern void av1_estimate_syntax_rate(
        MdRateEstimationContext_t      *md_rate_estimation_array,
        EbBool                          is_i_slice,
        FRAME_CONTEXT                  *fc);
    /**************************************************************************
    * Estimate the rate of the quantised coefficient
    * based on the frame CDF
    ***************************************************************************/
    extern void av1_estimate_coefficients_rate(
        MdRateEstimationContext_t  *md_rate_estimation_array,
        FRAME_CONTEXT              *fc);
    /**************************************************************************
    * av1_estimate_mv_rate()
    * Estimate the rate of motion vectors
    * based on the frame CDF
    ***************************************************************************/
    extern void av1_estimate_mv_rate(
#if ICOPY
        struct PictureControlSet_s     *picture_control_set_ptr,
#endif
        MdRateEstimationContext_t  *md_rate_estimation_array,
        nmv_context                *nmv_ctx);


#ifdef __cplusplus
}
#endif

#endif //EbMdRateEstimationTables_h