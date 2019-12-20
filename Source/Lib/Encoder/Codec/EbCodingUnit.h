// clang-format off
/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbCodingUnit_h
#define EbCodingUnit_h

#include "EbMotionEstimationLcuResults.h"
#include "EbPictureBufferDesc.h"
#include "EbPredictionUnit.h"
#include "EbTransformUnit.h"
#include "EbCabacContextModel.h"
#include "hash.h"


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
#define MAX_OIS_CANDIDATES  61  //18//18

    static const uint32_t intra_hev_cmode_to_intra_av1_mode[35] =
    {
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
    static const int8_t hevcMode_to_angle_delta_map[35] =
    {
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
    static const uint32_t intra_luma_to_chroma[INTRA_MODES] = // EB_INTRA_PLANAR
    {
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

    static const TxType chroma_transform_type[14] =
    {
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
    static const uint8_t av1_is_directional_chroma[UV_INTRA_MODES] =
    {
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
    static const uint8_t av1_is_directional[35] =
    {
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

    static INLINE EbBool av1_use_angle_delta(BlockSize bsize) {
        return bsize >= BLOCK_8X8;
    }

    struct PictureControlSet;

    typedef struct MV
    {
        int16_t row;
        int16_t col;
    } MV;

    typedef union  IntMv
    {
        uint32_t as_int;
        MV as_mv;
    } IntMv; /* facilitates faster equality tests and copies */

    typedef struct mv32 {
        int32_t row;
        int32_t col;
    } MV32;

    typedef struct CandidateMv
    {
        IntMv this_mv;
        IntMv comp_mv;
        int32_t weight;
    } CandidateMv;

#define INTER_TX_SIZE_BUF_LEN 16
#define TXK_TYPE_BUF_LEN 64

    typedef struct FilterIntraModeInfo {
        /*!< Specifies the type of intra filtering, and can represent any of the following:
         * FILTER_DC_PRED, FILTER_V_PRED, FILTER_H_PRED, FILTER_D157_PRED, FILTER_PAETH_PRED */
        FilterIntraMode filter_intra_mode;

        /*!< This bit specifies whether or not intra filtering can be used. */
        uint8_t use_filter_intra;
    } FilterIntraModeInfo_t;

    typedef struct InterIntraModeParams {
        /*!< Specifies the type of intra prediction to be used */
        InterIntraMode interintra_mode;

        /*!< equal to 1 specifies that wedge blending should be used.
            * wedge_interintra equal to 0 specifies that intra blending should be used. */
        uint8_t wedge_interintra;

        /*!< Used to derive the direction and offset of the wedge mask used during blending. */
        uint8_t interintra_wedge_index;

        /*!< Specifies the sign of the wedge blend. */
        // int interintra_wedge_sign; Always 0
    } InterIntraModeParams;

    typedef struct BlockModeInfo {
        // Common for both INTER and INTRA blocks
        BlockSize          sb_type;
        PredictionMode      mode;
        int8_t              skip;

        PartitionType       partition;

        /*!< 1 indicates that this block will use some default settings and skip mode info.
            * 0 indicates that the mode info is not skipped. */
        int8_t              skip_mode;

        /*!< Specifies which segment is associated with the current intra block being decoded. */
        int8_t segment_id;

        /*!< Equal to 1 specifies that the segment_id is taken from the segmentation map. */
        int8_t seg_id_predicted;

        /*!< For Lossy mode   : Specifies number of TUs in a block for each plane
             For Lossless mode: Specifies number of TUs for a block of size other than
                                128x128, 128x64, 64x128 and 64x64 - computed based on blocksize */
        uint8_t         num_tus[MAX_MB_PLANE - 1];

        /*!< Offset of first transform info from strat of SB pointer for each plane */
        uint16_t        first_tu_offset[MAX_MB_PLANE - 1];

        // Only for INTRA blocks
        UvPredictionMode   uv_mode;

        uint8_t             use_intrabc;

        // Only for INTER blocks

        MvReferenceFrame    ref_frame[2];
        IntMv               mv[2];

        uint16_t            ref_mv_idx;

        // interinter members

        InterIntraModeParams    interintra_mode_params;

        /*!< Specifies the type of motion compensation to perform. */
        MotionMode         motion_mode;

        InterIntraMode     is_inter_intra;

        /*!< 0 indicates that a distance based weighted scheme should be used for blending.
         *   1 indicates that the averaging scheme should be used for blending.*/
        uint8_t            compound_idx;

        InterInterCompoundData  inter_inter_compound;
        FilterIntraModeInfo_t  filter_intra_mode_info;

        /*!< Specifies how the motion vector used by inter prediction is obtained when using compound prediction. */
        uint8_t             compound_mode;

        /*!< Specifies the type of filter used in inter prediction. Values 0..3 are allowed
        * with the same interpretation as for interpolation_filter. One filter type is specified
        * for the vertical filter direction and one for the horizontal filter direction.*/
        uint32_t interp_filters;

        /*!< Index of the alpha Cb and alpha Cr combination */
        uint8_t cfl_alpha_idx;

        /*!< Contains the sign of the alpha values for U and V packed together into a single syntax element. */
        uint8_t cfl_alpha_signs;

        /*!< The actual prediction angle is the base angle + (angle_delta * step). */
        int8_t angle_delta[PLANE_TYPES];

        // Number of base colors for Y (0) and UV (1)
        uint8_t palette_size[MAX_MB_PLANE-1];

#if MODE_INFO_DBG
        int32_t mi_row;
        int32_t mi_col;
#endif
    } BlockModeInfo;

    typedef struct MbModeInfo
    {
#if CONFIG_RD_DEBUG
        RD_STATS rd_stats;
        int32_t mi_row;
        int32_t mi_col;
#endif
        EbWarpedMotionParams wm_params;
        int32_t comp_group_idx;

        int8_t cdef_strength;
        TxSize tx_size;
        uint8_t inter_tx_size[INTER_TX_SIZE_BUF_LEN];
        uint8_t tx_depth;
        BlockModeInfo block_mi;
#if PAL_SUP
        PaletteModeInfo palette_mode_info;
#endif
    } MbModeInfo;

    typedef struct {
        IntMv mfmv0;
        uint8_t ref_frame_offset;
    } TPL_MV_REF;
    typedef struct {
        IntMv mv;
        MvReferenceFrame ref_frame;
    } MV_REF;

    typedef struct ModeInfo {
        MbModeInfo mbmi;
    } ModeInfo;

    typedef struct TileInfo
    {
        int32_t mi_row_start, mi_row_end;
        int32_t mi_col_start, mi_col_end;
        int32_t tg_horz_boundary;
        int32_t tile_row;
        int32_t tile_col;
    } TileInfo;

    typedef struct MacroBlockDPlane
    {
        int subsampling_x;
        int subsampling_y;
        struct Buf2D dst;
        struct Buf2D pre[2];
        // block size in pixels
        uint8_t width, height;
    } MacroBlockDPlane;

    typedef struct MacroBlockPlane
    {
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
        const int16_t *quant_fp_QTX;
        const int16_t *round_fp_QTX;
        const int16_t *quant_QTX;
        const int16_t *quant_shift_QTX;
        const int16_t *zbin_QTX;
        const int16_t *round_QTX;
        const int16_t *dequant_QTX;
*/
    } MacroBlockPlane;

    struct buf_2d {
        uint8_t *buf;
        uint8_t *buf0;
        int width;
        int height;
        int stride;
    };
    typedef struct macroblockd_plane {
        int subsampling_x;
        int subsampling_y;
        struct buf_2d dst;
        struct buf_2d pre[2];
        uint8_t width, height;
    } MACROBLOCKD_PLANE;

    typedef struct MacroBlockD
    {
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

        uint8_t  use_intrabc;
        MbModeInfo *above_mbmi;
        MbModeInfo *left_mbmi;
        MbModeInfo *chroma_above_mbmi;
        MbModeInfo *chroma_left_mbmi;
        FRAME_CONTEXT *tile_ctx;
        TXFM_CONTEXT *above_txfm_context;
        TXFM_CONTEXT *left_txfm_context;
        TXFM_CONTEXT left_txfm_context_buffer[MAX_MIB_SIZE];
        struct macroblockd_plane plane[MAX_MB_PLANE];
#if OBMC_FLAG
        BlockSize sb_type;
#endif
    } MacroBlockD;

    typedef struct Macroblock
    {
        int32_t rdmult;
        int32_t switchable_restore_cost[RESTORE_SWITCHABLE_TYPES];
        int32_t wiener_restore_cost[2];
        int32_t sgrproj_restore_cost[2];
    } Macroblock;

    typedef struct IntraBcContext
    {
        int32_t rdmult;
        struct MacroBlockDPlane xdplane[MAX_MB_PLANE];
        struct MacroBlockPlane plane[MAX_MB_PLANE];
        MvLimits mv_limits;
        // The equivalend SAD error of one (whole) bit at the current quantizer
       // for large blocks.
        int sadperbit16;
        // The equivalent error at the current rdmult of one whole bit (not one
        // bitcost unit).
        int errorperbit;
        // Store the best motion vector during motion search
        IntMv best_mv;
        // Store the second best motion vector during full-pixel motion search
        IntMv second_best_mv;
        MacroBlockD * xd;
        int* nmv_vec_cost;
        int **mv_cost_stack;
        // buffer for hash value calculation of a block
        // used only in av1_get_block_hash_value()
        // [first hash/second hash]
        // [two buffers used ping-pong]
        uint32_t *hash_value_buffer[2][2];
        uint8_t  is_exhaustive_allowed;
        CRC_CALCULATOR crc_calculator1;
        CRC_CALCULATOR crc_calculator2;
    } IntraBcContext;

    typedef struct CodingUnit
    {
        TransformUnit             transform_unit_array[TRANSFORM_UNIT_MAX_COUNT]; // 2-bytes * 21 = 42-bytes
        PredictionUnit            prediction_unit_array[MAX_NUM_OF_PU_PER_CU];    // 35-bytes * 4 = 140 bytes
        InterInterCompoundData                 interinter_comp;
        uint8_t                                compound_idx;
        uint8_t                                comp_group_idx;
        unsigned                    skip_flag_context       : 2;
        unsigned                    prediction_mode_flag    : 2;
        unsigned                    block_has_coeff         : 1;
        unsigned                    split_flag_context      : 2;
#if !ADD_DELTA_QP_SUPPORT
        unsigned                    qp                      : 6;
        signed                      delta_qp                : 8; // can be signed 8bits
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
            unsigned                mdc_split_flag      : 1;
        };
#if NO_ENCDEC
        EbPictureBufferDesc      *quant_tmp;
        EbPictureBufferDesc      *coeff_tmp;
        EbPictureBufferDesc      *recon_tmp;
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
        uint8_t                     reference_mode_context;
        uint8_t                     compoud_reference_type_context;
        int32_t                     quantized_dc[3][MAX_TXB_COUNT];
        uint32_t                    is_inter_ctx;
        uint32_t                    interp_filters;
        uint8_t                      segment_id;
        uint8_t                      seg_id_predicted;  // valid only when temporal_update is enabled
        PartitionType               part;
        PART                        shape;
        uint16_t                    mds_idx;     //equivalent of leaf_index in the nscu context. we will keep both for now and use the right one on a case by case basis.
        uint8_t                    *neigh_left_recon[3];  //only for MD
        uint8_t                    *neigh_top_recon[3];
        uint16_t                   *neigh_left_recon_16bit[3];
        uint16_t                   *neigh_top_recon_16bit[3];
        uint32_t                    best_d1_blk;
        uint8_t                     tx_depth;
#if II_COMP_FLAG
        INTERINTRA_MODE             interintra_mode;
        uint8_t                     is_interintra_used;
        uint8_t                     use_wedge_interintra;
        int32_t                     interintra_wedge_index;
        int32_t                     ii_wedge_sign;
#endif
       uint8_t                      filter_intra_mode;
#if PAL_SUP
       PaletteInfo                          palette_info;
#endif
    } CodingUnit;

        typedef struct OisCandidate
        {
        union {
            struct {
                unsigned distortion : 20;
                unsigned valid_distortion : 1;
                unsigned : 3;
                unsigned intra_mode : 8;
            };
            uint32_t ois_results;
        };
        int32_t angle_delta;
    } OisCandidate;

    typedef struct OisSbResults
    {
        uint8_t             total_ois_intra_candidate[CU_MAX_COUNT];
        OisCandidate*    ois_candidate_array[CU_MAX_COUNT];
        int8_t              best_distortion_index[CU_MAX_COUNT];
    } OisSbResults;

    typedef struct EdgeLcuResults
    {
        uint8_t  edge_block_num;
        uint8_t  isolated_high_intensity_sb; // to be cleanedup
    } EdgeLcuResults;

    typedef struct SuperBlock
    {
        EbDctor                       dctor;
        struct PictureControlSet     *picture_control_set_ptr;

        CodingUnit                   *final_cu_arr;
        uint32_t                      final_cu_count;
        //for memory free only
        MacroBlockD                  *av1xd;
        PartitionType                  *cu_partition_array;
#if !ADD_DELTA_QP_SUPPORT
        unsigned                        qp                      : 8;
#endif
        unsigned                        picture_left_edge_flag  : 1;
        unsigned                        picture_top_edge_flag   : 1;
        unsigned                        picture_right_edge_flag : 1;
        unsigned                        index                   : 12;
        unsigned                        origin_x                : 12;
        unsigned                        origin_y                : 12;
#if ADD_DELTA_QP_SUPPORT
        uint16_t                        qp;
        int16_t                         delta_qp;
        int16_t                         org_delta_qp;
#endif
        uint32_t                        total_bits;

        // Quantized Coefficients
        EbPictureBufferDesc          *quantized_coeff;
#if MDC_ADAPTIVE_LEVEL
        uint8_t                       depth_ranking[NUMBER_OF_DEPTH];
#endif
#if MDC_ADAPTIVE_LEVEL || MULTI_PASS_PD
        uint64_t                      depth_cost[NUMBER_OF_DEPTH];
#endif
        TileInfo tile_info;
    } SuperBlock;

    extern EbErrorType largest_coding_unit_ctor(
        SuperBlock                  *larget_coding_unit_ptr,
        uint8_t                        sb_sz,
        uint16_t                       sb_origin_x,
        uint16_t                       sb_origin_y,
        uint16_t                       sb_index,
        struct PictureControlSet    *picture_control_set);

#ifdef __cplusplus
}
#endif
#endif // EbCodingUnit_h
// clang-format on
