/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbCodingUnit_h
#define EbCodingUnit_h

#include "EbDefinitions.h"
#include "EbSyntaxElements.h"
#include "EbMotionEstimationLcuResults.h"
#include "EbDefinitions.h"
#include "EbPictureBufferDesc.h"
#include "EbPredictionUnit.h"
#include "EbTransformUnit.h"
#include "EbCabacContextModel.h"

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
#define MAX_LCU_SIZE_IN_4X4BLK                                  (BLOCK_SIZE_64 >> 2)
#define VERTICAL_EDGE_BS_ARRAY_SIZE                             (MAX_LCU_SIZE_IN_4X4BLK * MAX_LCU_SIZE_IN_4X4BLK)
#define HORIZONTAL_EDGE_BS_ARRAY_SIZE                           (MAX_LCU_SIZE_IN_4X4BLK * MAX_LCU_SIZE_IN_4X4BLK)

#define MAX_NUMBER_OF_BS_EDGES_PER_TREEBLOCK    128
#define MAX_NUMBER_OF_LEAFS_PER_TREEBLOCK       64
#define MAX_NUMBER_OF_4x4_TUs_IN_8x8_LEAF       4
#define MAX_CU_COST (0xFFFFFFFFFFFFFFFFull >> 1)
#define MAX_MODE_COST ( 13616969489728 * 8) // RDCOST(6544618, 128 * 128 * 255 * 255, 128 * 128 * 255 * 255) * 8;
#define INVALID_FAST_CANDIDATE_INDEX    ~0
#define MAX_OPEN_LOOP_INTRA_CANDIDATES  18//18

    static const uint32_t intra_hev_cmode_to_intra_av1_mode[35] = {
        /*SMOOTH_PRED   */  SMOOTH_PRED,                                                        // EB_INTRA_PLANAR
        /*DC_PRED       */  DC_PRED,                                                            // EB_INTRA_DC
        /*D203_PRED     */  D203_PRED, D203_PRED, D203_PRED, D203_PRED, D203_PRED, D203_PRED,   // EB_INTRA_MODE_2 -> EB_INTRA_MODE_7
        /*H_PRED        */  H_PRED,  H_PRED, H_PRED, H_PRED, H_PRED,                            // EB_INTRA_MODE_8 -> EB_INTRA_HORIZONTAL -> EB_INTRA_MODE_12
        /*D157_PRED     */  D157_PRED, D157_PRED, D157_PRED,                                    // EB_INTRA_MODE_13 -> EB_INTRA_MODE_15
        /*D135_PRED     */  D135_PRED, D135_PRED, D135_PRED, D135_PRED, D135_PRED,              // EB_INTRA_MODE_16 -> EB_INTRA_MODE_20
        /*D113_PRED     */  D113_PRED, D113_PRED, D113_PRED,                                    // EB_INTRA_MODE_21 -> EB_INTRA_MODE_23
        /*V_PRED        */  V_PRED, V_PRED, V_PRED,  V_PRED,                                    // EB_INTRA_MODE_24 -> EB_INTRA_VERTICAL -> EB_INTRA_MODE_27
        /*D67_PRED      */  D67_PRED, D67_PRED, D67_PRED,                                       // EB_INTRA_MODE_28 -> EB_INTRA_MODE_30
        /*D45_PRED      */  D45_PRED, D45_PRED, D45_PRED, D45_PRED,                             // EB_INTRA_MODE_31 -> EB_INTRA_MODE_34
    };
    static const int8_t hevcMode_to_angle_delta_map[35] = {
        /*SMOOTH_PRED   */   0,                                                                 // EB_INTRA_PLANAR ->
        /*DC_PRED       */   0,                                                                 // EB_INTRA_DC
        /*D203_PRED     */   3, 2, 1, 0, -1, -2,                                                // EB_INTRA_MODE_2 -> EB_INTRA_MODE_7
        /*H_PRED        */   3, 2, 0, -2, -3,                                                   // EB_INTRA_MODE_8 -> EB_INTRA_HORIZONTAL -> EB_INTRA_MODE_12
        /*D157_PRED     */   2, 0, -2,                                                          // EB_INTRA_MODE_13 -> EB_INTRA_MODE_15
        /*D135_PRED     */   3, 2, 0, -2, -3,                                                   // EB_INTRA_MODE_16 -> EB_INTRA_MODE_20
        /*D113_PRED     */   2, 0, -2,                                                          // EB_INTRA_MODE_21 -> EB_INTRA_MODE_23
        /*V_PRED        */   3, 2, 0, -2,                                                       // EB_INTRA_MODE_24 -> EB_INTRA_VERTICAL -> EB_INTRA_MODE_27
        /*D67_PRED      */   2, 0, -2,                                                          // EB_INTRA_MODE_28 -> EB_INTRA_MODE_30
        /*D45_PRED      */   3, 2, 0, -2,                                                       // EB_INTRA_MODE_31 -> EB_INTRA_MODE_34
    };
#if IMPROVE_CHROMA_MODE
    static const uint32_t intra_luma_to_chroma[INTRA_MODES] = {                                                                            // EB_INTRA_PLANAR
       UV_DC_PRED,        // Average of above and left pixels
       UV_V_PRED,         // Vertical
       UV_H_PRED,         // Horizontal
       UV_D45_PRED,       // Directional 45  degree
       UV_D135_PRED,      // Directional 135 degree
       UV_D113_PRED,      // Directional 113 degree
       UV_D157_PRED,      // Directional 157 degree
       UV_D203_PRED,      // Directional 203 degree
       UV_D67_PRED,       // Directional 67  degree
       UV_SMOOTH_PRED,    // Combination of horizontal and vertical interpolation
       UV_SMOOTH_V_PRED,  // Vertical interpolation
       UV_SMOOTH_H_PRED,  // Horizontal interpolation
       UV_PAETH_PRED,     // Predict from the direction of smallest gradient
    };
#else
    static const uint32_t intra_luma_to_chroma[INTRA_MODES] = {                                                                            // EB_INTRA_PLANAR
        /*DC_PRED       */  UV_DC_PRED,
        /*V_PRED        */  UV_SMOOTH_PRED,
        /*H_PRED        */  UV_SMOOTH_PRED,
        /*D45_PRED      */  UV_D45_PRED,
        /*D135_PRED     */  UV_D135_PRED,
        /*D113_PRED     */  UV_D135_PRED,
        /*D157_PRED     */  UV_D135_PRED,
        /*D203_PRED     */  UV_SMOOTH_PRED,
        /*D67_PRED      */  UV_D45_PRED,
        /*SMOOTH_PRED   */  UV_SMOOTH_PRED,
        /*SMOOTH_V_PRED */  UV_SMOOTH_PRED,
        /*SMOOTH_H_PRED */  UV_SMOOTH_PRED,
        /*PAETH_PRED    */  UV_PAETH_PRED,
    };
#endif
    static const TxType chroma_transform_type[14] = {
        /*UV_DC_PRED,          */   DCT_DCT   ,
        /*UV_V_PRED,           */   ADST_DCT  ,
        /*UV_H_PRED,           */   DCT_ADST  ,
        /*UV_D45_PRED,         */   DCT_DCT   ,
        /*UV_D135_PRED,        */   ADST_ADST ,
        /*UV_D113_PRED,        */   ADST_DCT  ,
        /*UV_D157_PRED,        */   DCT_ADST  ,
        /*UV_D203_PRED,        */   DCT_ADST  ,
        /*UV_D67_PRED,         */   ADST_DCT  ,
        /*UV_SMOOTH_PRED,      */   ADST_ADST ,
        /*UV_SMOOTH_V_PRED,    */   ADST_DCT  ,
        /*UV_SMOOTH_H_PRED,    */   DCT_ADST  ,
        /*UV_PAETH_PRED,       */   ADST_ADST ,
        /*UV_CFL_PRED,          */  DCT_DCT,
    };
    static const uint8_t av1_is_directional_chroma[UV_INTRA_MODES] = {
        /*UV_DC_PRED,        */  0,
        /*UV_V_PRED,         */  0,
        /*UV_H_PRED,         */  1,
        /*UV_D45_PRED,       */  1,
        /*UV_D135_PRED,      */  1,
        /*UV_D113_PRED,      */  1,
        /*UV_D157_PRED,      */  1,
        /*UV_D203_PRED,      */  1,
        /*UV_D67_PRED,       */  1,
        /*UV_SMOOTH_PRED,    */  0,
        /* UV_SMOOTH_PRED,   */  0,
        /* UV_SMOOTH_V_PRED, */  0,
        /* UV_SMOOTH_H_PRED, */  0,
        /* UV_PAETH_PRED,    */  0,
    };
    static const uint8_t av1_is_directional[35] = {
        0,                           // EB_INTRA_PLANAR
        0,                           // EB_INTRA_DC
        1, 1, 1, 1, 1, 1, 1, 1,      // EB_INTRA_MODE_2 -> EB_INTRA_MODE_9
        1,                           // EB_INTRA_HORIZONTAL
        1, 1, 1, 1, 1,               // EB_INTRA_MODE_11 -> EB_INTRA_MODE_15
        1, 1, 1, 1, 1,               // EB_INTRA_MODE_16 -> EB_INTRA_MODE_20
        1, 1, 1, 1, 1,               // EB_INTRA_MODE_21 -> EB_INTRA_MODE_25
        1,                           // EB_INTRA_VERTICAL
        1, 1, 1, 1,                  // EB_INTRA_MODE_27 -> EB_INTRA_MODE_30
        1, 1, 1, 1,                  // EB_INTRA_MODE_31 -> EB_INTRA_MODE_34
    };
    static const uint32_t mode_to_angle_map[] = {
        0, 90, 180, 45, 135, 113, 157, 203, 67, 0, 0, 0, 0,
    };
    static INLINE int32_t av1_is_directional_mode(PredictionMode mode) {
        return mode >= V_PRED && mode <= D67_PRED;
    }
    struct PictureControlSet_s;

    typedef struct MV {
        int16_t row;
        int16_t col;
    } MV;
    typedef union  IntMv {
        uint32_t as_int;
        MV as_mv;
    } IntMv; /* facilitates faster equality tests and copies */
    typedef struct CandidateMv {
        IntMv this_mv;
        IntMv comp_mv;
        int32_t weight;
    } CandidateMv;

    typedef struct MbModeInfo {
        // Common for both INTER and INTRA blocks
        block_size sb_type;
        PredictionMode mode;
        //TxSize tx_size;
        //uint8_t inter_tx_size[INTER_TX_SIZE_BUF_LEN];
        //int8_t skip;
        //int8_t skip_mode;
        //int8_t segment_id;
        //int8_t seg_id_predicted;  // valid only when temporal_update is enabled
        // Only for INTRA blocks
        UV_PredictionMode uv_mode;
        //PALETTE_MODE_INFO palette_mode_info;
        //uint8_t use_intrabc;
        // Only for INTER blocks
        //InterpFilters interp_filters;
        MvReferenceFrame ref_frame[2];
        //TxType txk_type[TXK_TYPE_BUF_LEN];
        //FILTER_INTRA_MODE_INFO filter_intra_mode_info;
        // The actual prediction angle is the base angle + (angle_delta * step).
        //int8_t angle_delta[PLANE_TYPES];
        // interintra members
        //INTERINTRA_MODE interintra_mode;
        // TODO(debargha): Consolidate these flags
        //int32_t use_wedge_interintra;
        //int32_t interintra_wedge_index;
        //int32_t interintra_wedge_sign;
        // interinter members
        //COMPOUND_TYPE interinter_compound_type;
        //int32_t wedge_index;
        //int32_t wedge_sign;
        //SEG_MASK_TYPE mask_type;
        //MOTION_MODE motion_mode;
        //int32_t overlappable_neighbors[2];
        IntMv mv[2];
        //IntMv pred_mv[2];
        //uint8_t ref_mv_idx;
        PartitionType partition;
        /* deringing gain *per-superblock* */
        //int8_t cdef_strength;
        //int32_t current_q_index;
        //int32_t current_delta_lf_from_base;
        //int32_t curr_delta_lf[FRAME_LF_COUNT];
#if CONFIG_RD_DEBUG
        RD_STATS rd_stats;
        int32_t mi_row;
        int32_t mi_col;
#endif
        EbWarpedMotionParams wm_params;
        // Index of the alpha Cb and alpha Cr combination
        // int32_t cfl_alpha_idx;
        // Joint sign of alpha Cb and alpha Cr
        // int32_t cfl_alpha_signs;
        //int32_t compound_idx;
        //int32_t comp_group_idx;
        int8_t skip;
        int8_t cdef_strength;
        TxSize tx_size;
    } MbModeInfo;
    typedef struct ModeInfo {
        MbModeInfo mbmi;
    } ModeInfo;
    typedef struct TileInfo {
        int32_t mi_row_start, mi_row_end;
        int32_t mi_col_start, mi_col_end;
        int32_t tg_horz_boundary;
#if TILES
        int32_t tile_row;
        int32_t tile_col;
#endif
    } TileInfo;
    typedef struct MacroBlockD {
        // block dimension in the unit of mode_info.
        uint8_t n8_w, n8_h;
        uint8_t n4_w, n4_h;  // for warped motion
        uint8_t ref_mv_count[MODE_CTX_REF_FRAMES];
        CandidateMv final_ref_mv_stack[MAX_REF_MV_STACK_SIZE];       
        uint8_t is_sec_rect;
        int32_t up_available;
        int32_t left_available;
        TileInfo tile;
        int32_t mi_stride;
        ModeInfo **mi;

        /* Distance of MB away from frame edges in subpixels (1/8th pixel)  */
        int32_t mb_to_left_edge;
        int32_t mb_to_right_edge;
        int32_t mb_to_top_edge;
        int32_t mb_to_bottom_edge;
        uint8_t neighbors_ref_counts[TOTAL_REFS_PER_FRAME];
    } MacroBlockD;
    typedef struct Macroblock {
        int32_t rdmult;
        int32_t switchable_restore_cost[RESTORE_SWITCHABLE_TYPES];
        int32_t wiener_restore_cost[2];
        int32_t sgrproj_restore_cost[2];
    } Macroblock;
    typedef struct CodingUnit_s
    {
        TransformUnit_t             transform_unit_array[TRANSFORM_UNIT_MAX_COUNT]; // 2-bytes * 21 = 42-bytes
        PredictionUnit_t            prediction_unit_array[MAX_NUM_OF_PU_PER_CU];    // 35-bytes * 4 = 140 bytes

        unsigned                    skip_flag_context       : 2;
        unsigned                    prediction_mode_flag    : 2;
        unsigned                    block_has_coeff         : 1;
        unsigned                    split_flag_context      : 2;
#if !ADD_DELTA_QP_SUPPORT
        unsigned                    qp                      : 6;
        unsigned                    ref_qp                  : 6;
        signed                      delta_qp                : 8; // can be signed 8bits
        signed                      org_delta_qp            : 8;
#else
        uint16_t                    qp;
        uint16_t                    ref_qp;
        int16_t                     delta_qp; // can be signed 8bits
        int16_t                     org_delta_qp;
#endif
        // Coded Tree
        struct {
            unsigned                leaf_index           : 8;
            unsigned                split_flag           : 1;
            unsigned                skip_flag            : 1;
#if FIX_INTER_DEPTH
            unsigned                mdc_split_flag      : 1;
#endif
        };
#if NO_ENCDEC
        EbPictureBufferDesc_t      *quant_tmp;
        EbPictureBufferDesc_t      *coeff_tmp;
        EbPictureBufferDesc_t      *recon_tmp;
        uint32_t                    cand_buff_index;
#endif
        MacroBlockD                *av1xd;
        // uint8_t ref_mv_count[MODE_CTX_REF_FRAMES];
        int16_t                     inter_mode_ctx[MODE_CTX_REF_FRAMES];
        IntMv                       ref_mvs[MODE_CTX_REF_FRAMES][MAX_MV_REF_CANDIDATES]; //used only for nonCompound modes.
        uint8_t                     drl_index;
        PredictionMode              pred_mode;
        IntMv                       predmv[2];
        uint8_t                     skip_coeff_context;
        int16_t                     luma_txb_skip_context;
        int16_t                     luma_dc_sign_context;
        int16_t                     cb_txb_skip_context;
        int16_t                     cb_dc_sign_context;
        int16_t                     cr_txb_skip_context;
        int16_t                     cr_dc_sign_context;
        uint8_t                     reference_mode_context;
        uint8_t                     compoud_reference_type_context;
        int32_t                     quantized_dc[3];
        uint32_t                    is_inter_ctx;
        uint32_t                    interp_filters;
        PartitionType               part;
        PART                        shape;
        uint16_t                    mds_idx;     //equivalent of leaf_index in the nscu context. we will keep both for now and use the right one on a case by case basis.
        uint8_t                    *neigh_left_recon[3];  //only for MD
        uint8_t                    *neigh_top_recon[3];
        uint32_t                    best_d1_blk;
    } CodingUnit_t;
    typedef struct OisCandidate_s {
        union {
            struct {
                unsigned distortion : 20;
                unsigned valid_distortion : 1;
                unsigned : 3;
                unsigned intra_mode : 8;
            };
            uint32_t ois_results;
        };
    } OisCandidate_t;
    typedef struct OisLcuResults_s
    {
        uint8_t           total_intra_luma_mode[CU_MAX_COUNT];
        OisCandidate_t    sorted_ois_candidate[CU_MAX_COUNT][MAX_OPEN_LOOP_INTRA_CANDIDATES];
    } OisLcuResults_t;
    typedef struct OisCu32Cu16Results_s
    {
        uint8_t            total_intra_luma_mode[21];
        OisCandidate_t*    sorted_ois_candidate[21];
    } OisCu32Cu16Results_t;
    typedef struct OisCu8Results_s
    {
        uint8_t            total_intra_luma_mode[64];
        OisCandidate_t*    sorted_ois_candidate[64];
    } OisCu8Results_t;
    typedef struct QpmLcuResults_s {
        uint8_t  cu_qp;
        uint8_t  cu_intra_qp;
        uint8_t  cu_inter_qp;
        int8_t   delta_qp;
        int8_t   inner_sb_cu_delta_qp;
    } QpmLcuResults_t; // to be cleaned up
    typedef struct EdgeLcuResults_s {
        uint8_t  edge_block_num;
        uint8_t  isolated_high_intensity_sb; // to be cleanedup
    } EdgeLcuResults_t;
    typedef struct LargestCodingUnit_s {
        struct PictureControlSet_s     *picture_control_set_ptr;
        CodingUnit_t                   *final_cu_arr;
        uint32_t                        tot_final_cu;
        PartitionType                  *cu_partition_array;

        // Coding Units
        EbAuraStatus                    aura_status_iii; // aura status for Gold 4K only, used in testing more depths
#if !ADD_DELTA_QP_SUPPORT
        unsigned                        qp                      : 8;
#endif                                                          
        unsigned                        picture_left_edge_flag  : 1;
        unsigned                        picture_top_edge_flag   : 1;
        unsigned                        picture_right_edge_flag : 1;
        unsigned                        pred64                  : 2;
        unsigned                        index                   : 12;
        unsigned                        origin_x                : 12;
        unsigned                        origin_y                : 12;
#if ADD_DELTA_QP_SUPPORT
        uint16_t                        qp;
        int16_t                         delta_qp;
        int16_t                         org_delta_qp;
#endif
        //Bits only used for quantized coeffs
        uint32_t                        quantized_coeffs_bits;
        uint32_t                        total_bits;

        // Quantized Coefficients
        EbPictureBufferDesc_t          *quantized_coeff;
#if TILES
        TileInfo tile_info;
#endif

    } LargestCodingUnit_t;

    extern EbErrorType largest_coding_unit_ctor(
        LargestCodingUnit_t          **larget_coding_unit_dbl_ptr,
        uint8_t                        sb_sz,
        uint16_t                       sb_origin_x,
        uint16_t                       sb_origin_y,
        uint16_t                       sb_index,
        struct PictureControlSet_s    *picture_control_set);

#ifdef __cplusplus
}
#endif
#endif // EbCodingUnit_h
