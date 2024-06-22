/*
* Copyright(c) 2019 Intel Corporation
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#ifndef EbCodingUnit_h
#define EbCodingUnit_h

#include "me_sb_results.h"
#include "pic_buffer_desc.h"
#include "block_structures.h"
#include "cabac_context_model.h"
#include "hash.h"
#include "definitions.h"
#include "motion_vector_unit.h"

#ifdef __cplusplus
extern "C" {
#endif
/*
    Requirements:
    -Must have enough CodingUnits for every single CU pattern
    -Easy to expand/insert CU
    -Easy to collapse a CU
    -Easy to replace CUs
    -Statically Allocated
    -Contains the leaf count
    */

// Macros for deblocking filter
#define MAX_CU_COST (0xFFFFFFFFFFFFFFFFull >> 1)
#define MAX_MODE_COST (13754408443200 * 8) // RDCOST(6544618, 128 * 128 * 255 * 255, 128 * 128 * 255 * 255) * 8;

static const uint32_t intra_luma_to_chroma[INTRA_MODES] = // EB_INTRA_PLANAR
    {
        UV_DC_PRED, // Average of above and left pixels
        UV_V_PRED, // Vertical
        UV_H_PRED, // Horizontal
        UV_D45_PRED, // Directional 45  degree
        UV_D135_PRED, // Directional 135 degree
        UV_D113_PRED, // Directional 113 degree
        UV_D157_PRED, // Directional 157 degree
        UV_D203_PRED, // Directional 203 degree
        UV_D67_PRED, // Directional 67  degree
        UV_SMOOTH_PRED, // Combination of horizontal and vertical interpolation
        UV_SMOOTH_V_PRED, // Vertical interpolation
        UV_SMOOTH_H_PRED, // Horizontal interpolation
        UV_PAETH_PRED, // Predict from the direction of smallest gradient
};

typedef struct {
    IntMv   mfmv0;
    uint8_t ref_frame_offset;
} TPL_MV_REF;
typedef struct {
    IntMv            mv;
    MvReferenceFrame ref_frame;
} MV_REF;

typedef struct ModeInfo {
    MbModeInfo mbmi;
} ModeInfo;

typedef struct MacroBlockDPlane {
    int          subsampling_x;
    int          subsampling_y;
    struct Buf2D dst;
    struct Buf2D pre[2];
    // block size in pixels
    uint8_t width, height;
} MacroBlockDPlane;

typedef struct MacroBlockPlane {
    struct Buf2D src;
    /*
        DECLARE_ALIGNED(16, int16_t, src_diff[MAX_SB_SQUARE]);
        TranLow *qcoeff;
        TranLow *coeff;
        uint16_t *eobs;
        uint8_t *txb_entropy_ctx;

        // Quantizer setings
        // These are used/accessed only in the quantization process
        // RDO does not / must not depend on any of these values
        // All values below share the coefficient scale/shift used in TX
        const int16_t *quant_fp_qtx;
        const int16_t *round_fp_qtx;
        const int16_t *quant_qtx;
        const int16_t *quant_shift_qtx;
        const int16_t *zbin_qtx;
        const int16_t *round_qtx;
        const int16_t *dequant_qtx;
*/
} MacroBlockPlane;

typedef struct MacroBlockD {
    // block dimension in the unit of mode_info.
    uint8_t    n8_w, n8_h;
    uint8_t    n4_w, n4_h; // for warped motion
    uint8_t    ref_mv_count[MODE_CTX_REF_FRAMES];
    uint8_t    is_sec_rect;
    int8_t     up_available;
    int8_t     left_available;
    int8_t     chroma_up_available;
    int8_t     chroma_left_available;
    TileInfo   tile;
    int32_t    mi_stride;
    ModeInfo **mi;

    /* Distance of MB away from frame edges in subpixels (1/8th pixel)  */
    int32_t        mb_to_left_edge;
    int32_t        mb_to_right_edge;
    int32_t        mb_to_top_edge;
    int32_t        mb_to_bottom_edge;
    int            mi_row; // Row position in mi units
    int            mi_col; // Column position in mi units
    uint8_t        neighbors_ref_counts[TOTAL_REFS_PER_FRAME];
    MbModeInfo    *above_mbmi;
    MbModeInfo    *left_mbmi;
    MbModeInfo    *chroma_above_mbmi;
    MbModeInfo    *chroma_left_mbmi;
    FRAME_CONTEXT *tile_ctx;
    TXFM_CONTEXT  *above_txfm_context;
    TXFM_CONTEXT  *left_txfm_context;
    BlockSize      bsize;
} MacroBlockD;

typedef struct Macroblock {
    int32_t rdmult;
    int32_t switchable_restore_cost[RESTORE_SWITCHABLE_TYPES];
    int32_t wiener_restore_cost[2];
    int32_t sgrproj_restore_cost[2];
} Macroblock;

typedef struct IntraBcContext {
    int32_t                 rdmult;
    struct MacroBlockDPlane xdplane[MAX_MB_PLANE];
    struct MacroBlockPlane  plane[MAX_MB_PLANE];
    MvLimits                mv_limits;
    // The equivalend SAD error of one (whole) bit at the current quantizer
    // for large blocks.
    int sadperbit16;
    // The equivalent error at the current rdmult of one whole bit (not one
    // bitcost unit).
    int errorperbit;
    // Store the best motion vector during motion search
    IntMv best_mv;
    // Store the second best motion vector during full-pixel motion search
    IntMv        second_best_mv;
    MacroBlockD *xd;
    int         *nmv_vec_cost;
    int        **mv_cost_stack;
    // buffer for hash value calculation of a block
    // used only in svt_av1_get_block_hash_value()
    // [first hash/second hash]
    // [two buffers used ping-pong]
    uint32_t      *hash_value_buffer[2][2];
    uint8_t        is_exhaustive_allowed;
    CRC_CALCULATOR crc_calculator1;
    CRC_CALCULATOR crc_calculator2;
    // use approximate rate for inter cost (set at pic-level b/c some pic-level initializations will
    // be removed)
    uint8_t approx_inter_rate;
} IntraBcContext;
typedef struct EobData {
    uint16_t y[MAX_TXB_COUNT];
    uint16_t u[MAX_TXB_COUNT_UV];
    uint16_t v[MAX_TXB_COUNT_UV];
} EobData;

typedef struct QuantDcData {
    uint8_t y[MAX_TXB_COUNT];
    uint8_t u[MAX_TXB_COUNT_UV];
    uint8_t v[MAX_TXB_COUNT_UV];
} QuantDcData;

typedef struct BlkStruct {
    MacroBlockD *av1xd;
    // only for MD
    uint8_t  *neigh_left_recon[3];
    uint8_t  *neigh_top_recon[3];
    uint16_t *neigh_left_recon_16bit[3];
    uint16_t *neigh_top_recon_16bit[3];
    // buffer to store quantized coeffs from MD for the final mode of each block
    // Used when encdec is bypassed
    EbPictureBufferDesc *coeff_tmp;
    // buffer to store recon from MD for the final mode of each block
    // Used when encdec is bypassed
    EbPictureBufferDesc *recon_tmp;
    uint64_t             cost;
    // Similar to cost but does not get updated @ svt_aom_d1_non_square_block_decision() and
    // svt_aom_d2_inter_depth_block_decision()
    uint64_t     default_cost;
    uint64_t     total_rate;
    uint32_t     full_dist;
    QuantDcData  quant_dc;
    EobData      eob;
    TxType       tx_type[MAX_TXB_COUNT];
    TxType       tx_type_uv;
    uint16_t     y_has_coeff;
    uint8_t      u_has_coeff;
    uint8_t      v_has_coeff;
    PaletteInfo *palette_info;
    uint8_t      palette_mem; // status of palette info alloc
    uint8_t      palette_size[2];

    Mv      mv[MAX_NUM_OF_REF_PIC_LIST]; // 16-bytes
    IntMv   predmv[MAX_NUM_OF_REF_PIC_LIST];
    uint8_t inter_pred_direction_index;

    // Intra Mode
    int8_t           angle_delta[PLANE_TYPES];
    UvPredictionMode intra_chroma_mode;
    // Inter Mode
    uint8_t    ref_frame_type;
    MotionMode motion_mode;
    uint16_t   num_proj_ref;
    uint32_t   overlappable_neighbors;
    uint8_t    cfl_alpha_idx; // Index of the alpha Cb and alpha Cr combination
    uint8_t    cfl_alpha_signs; // Joint sign of alpha Cb and alpha Cr

    InterInterCompoundData interinter_comp;
    uint32_t               interp_filters;
    uint8_t                interintra_wedge_index;
    int16_t                inter_mode_ctx;
    // equivalent of leaf_index in the nscu context. we will keep both for now and use the right one
    // on a case by case basis.
    uint16_t mds_idx;
    // txb
    uint8_t tx_depth;
    uint8_t compound_idx;
    uint8_t comp_group_idx;

    uint8_t qindex;
    uint8_t split_flag;
    uint8_t skip_mode; // skips mode_info + coeff. as defined in section 6.10.10 of the av1 text
    uint8_t drl_index;
    // Store the drl ctx in coding loop to avoid storing final_ref_mv_stack and ref_mv_count for EC
    int8_t drl_ctx[2];
    // Store the drl ctx in coding loop to avoid storing final_ref_mv_stack and ref_mv_count for EC
    int8_t         drl_ctx_near[2];
    PredictionMode pred_mode;

    uint8_t segment_id;

    PartitionType  part;
    uint16_t       best_d1_blk;
    InterIntraMode interintra_mode;
    uint8_t        is_interintra_used;
    uint8_t        use_wedge_interintra;
    uint8_t        filter_intra_mode;
    uint8_t        use_intrabc;

    PartitionContextType left_neighbor_partition;
    PartitionContextType above_neighbor_partition;

    // wm
    EbWarpedMotionParams wm_params_l0;
    EbWarpedMotionParams wm_params_l1;

    unsigned cnt_nz_coeff : 12;
    unsigned prediction_mode_flag : 2;
    // ec; skip coeff only. as defined in section 6.10.11 of the av1 text
    unsigned block_has_coeff : 1;
} BlkStruct;
typedef struct EcBlkStruct {
    MacroBlockD *av1xd;
    EobData      eob;
    TxType       tx_type[MAX_TXB_COUNT];
    TxType       tx_type_uv;

    PaletteInfo           *palette_info;
    uint8_t                palette_size[2];
    IntMv                  predmv[2];
    InterInterCompoundData interinter_comp;

    // Intra Mode
    int8_t angle_delta[PLANE_TYPES];
    // Inter Mode
    MotionMode motion_mode;
    uint16_t   num_proj_ref;
    uint32_t   overlappable_neighbors;
    uint8_t    cfl_alpha_idx; // Index of the alpha Cb and alpha Cr combination
    uint8_t    cfl_alpha_signs; // Joint sign of alpha Cb and alpha Cr
    uint8_t    interintra_wedge_index;

    int16_t inter_mode_ctx;
    // equivalent of leaf_index in the nscu context. we will keep both for now and use the right one
    // on a case by case basis.
    uint16_t mds_idx;

    uint8_t qindex;

    uint8_t drl_index;
    // Store the drl ctx in coding loop to avoid storing final_ref_mv_stack and ref_mv_count for EC
    int8_t drl_ctx[2];
    // Store the drl ctx in coding loop to avoid storing final_ref_mv_stack and ref_mv_count for EC
    int8_t drl_ctx_near[2];

    InterIntraMode interintra_mode;
    uint8_t        is_interintra_used;
    uint8_t        use_wedge_interintra;
    uint8_t        filter_intra_mode;
} EcBlkStruct;

typedef struct TplStats {
    int64_t  srcrf_dist;
    int64_t  recrf_dist;
    int64_t  srcrf_rate;
    int64_t  recrf_rate;
    int64_t  mc_dep_rate;
    int64_t  mc_dep_dist;
    MV       mv;
    uint64_t ref_frame_poc;
} TplStats;

typedef struct TplSrcStats {
    int64_t        srcrf_dist;
    int64_t        srcrf_rate;
    uint64_t       ref_frame_poc;
    MV             mv;
    uint8_t        best_mode;
    int32_t        best_rf_idx;
    PredictionMode best_intra_mode;
} TplSrcStats;
typedef struct SuperBlock {
    EbDctor                   dctor;
    struct PictureControlSet *pcs;
    EcBlkStruct              *final_blk_arr;
    //for memory free only
    MacroBlockD   *av1xd;
    PartitionType *cu_partition_array;
    unsigned       index : 32;
    unsigned       org_x : 32;
    unsigned       org_y : 32;
    uint8_t        qindex;
    TileInfo       tile_info;
    uint16_t       final_blk_cnt; // number of block(s) posted from EncDec to EC
} SuperBlock;

extern EbErrorType svt_aom_largest_coding_unit_ctor(SuperBlock *larget_coding_unit_ptr, uint8_t sb_size,
                                                    uint16_t sb_origin_x, uint16_t sb_origin_y, uint16_t sb_index,
                                                    EncMode enc_mode, uint16_t max_block_cnt,

                                                    struct PictureControlSet *picture_control_set);

#ifdef __cplusplus
}
#endif
#endif // EbCodingUnit_h
