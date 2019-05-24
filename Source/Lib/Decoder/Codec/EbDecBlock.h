/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbDecBlock_h
#define EbDecBlock_h

#include "EbDefinitions.h"
#include "EbRestoration.h"

#define MODE_INFO_DBG   0

static const BlockSize Partition_Subsize[10][BlockSizeS_ALL] =
{ 
    { BLOCK_4X4, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_8X8, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_16X16, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_32X32, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_64X64, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_128X128, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID }, 
    { BLOCK_INVALID, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_8X4, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_16X8, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_32X16, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_64X32, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_128X64, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID }, 
    { BLOCK_INVALID, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_4X8, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_8X16, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_16X32, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_32X64, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_64X128, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID }, 
    { BLOCK_INVALID, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_4X4, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_8X8, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_16X16, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_32X32, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_64X64, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID },
    { BLOCK_INVALID, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_8X4, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_16X8, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_32X16, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_64X32, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_128X64, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID },
    { BLOCK_INVALID, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_8X4, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_16X8,
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_32X16,
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_64X32, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_128X64, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID },
    { BLOCK_INVALID, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_4X8, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_8X16, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_16X32, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_32X64, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_64X128, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID },
    { BLOCK_INVALID, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_4X8, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_8X16, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_16X32, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_32X64, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_64X128, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID } ,
    { BLOCK_INVALID, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_16X4, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_32X8, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_64X16, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID }, 
    { BLOCK_INVALID, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_4X16, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_8X32, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_16X64, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, 
    BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID }
};

enum {
    UNIFORM_45 = 0,
    UNIFORM_45_INV,
    COMPOUND_MASK_TYPES,
} UENUM1BYTE(COMPOUND_MASK_TYPE);

typedef struct MV_decs {
    int16_t row;
    int16_t col;
} MV_dec;
typedef union  IntMv_dec {
    uint32_t as_int;
    MV_dec as_mv;
} IntMv_dec; /* facilitates faster equality tests and copies */
typedef struct CandidateMv_dec {
    IntMv_dec this_mv;
    IntMv_dec comp_mv;
    int32_t weight;
} CandidateMv_dec;

typedef struct FilterIntraModeInfo {

    /*!< Specifies the type of intra filtering, and can represent any of the following:
     * FILTER_DC_PRED, FILTER_V_PRED, FILTER_H_PRED, FILTER_D157_PRED, FILTER_PAETH_PRED */
    FilterIntraMode filter_intra_mode;

    /*!< This bit specifies whether or not intra filtering can be used. */
    uint8_t use_filter_intra;

} FilterIntraModeInfo_t;


typedef struct InterIntraMode {

    /*!< Specifies the type of intra prediction to be used */
    InterIntraMode interintra_mode;

    /*!< equal to 1 specifies that wedge blending should be used.
        * wedge_interintra equal to 0 specifies that intra blending should be used. */
    uint8_t wedge_interintra;

    /*!< Used to derive the direction and offset of the wedge mask used during blending. */
    uint8_t interintra_wedge_index;

    /*!< Specifies the sign of the wedge blend. */
    // int interintra_wedge_sign; Always 0

} InterIntraMode_t;

typedef struct {

    /*!< Indicates that the compound_idx syntax element should be read or not. */
    uint8_t comp_group_idx;

    /*!< 0 indicates that a distance based weighted scheme should be used for blending.
     *   1 indicates that the averaging scheme should be used for blending.*/
    uint8_t compound_idx;

    /*!< Specifies how the two predictions should be blended together. */
    CompoundType type;

    /*!< Used to derive the direction and offset of the wedge mask used during blending. */
    uint8_t wedge_index;

    /*!< Specifies the sign of the wedge blend. */
    uint8_t wedge_sign;

    /*!< Specifies the type of mask to be used during blending. */
    COMPOUND_MASK_TYPE mask_type;

} InterCompoundData_t;

typedef struct TransformInfo {

    /*!< Specifies the transform size to be used for this TU. */
    TxSize tx_size;

    /*!< Specifies the transform type for the TU. */
    TxType txk_type;

    int8_t  cbf;

} TransformInfo_t;

typedef struct ModeInfo_t {

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

    /*!< Number of Luma TUs in block */
    uint8_t         num_luma_tus;

    /*!< Offset of first Luma transform info from strat of SB pointer */
    uint16_t        first_luma_tu_offset;

    /*!< Number of Chroma TUs in block */
    uint8_t         num_chroma_tus;

    /*!< Offset of first Chroma transform info from strat of SB pointer */
    uint16_t        first_chroma_tu_offset;

    // Only for INTRA blocks
    UvPredictionMode   uv_mode;

    uint8_t             use_intrabc;

    // Only for INTER blocks

    MvReferenceFrame    ref_frame[2];
    IntMv_dec           mv[2];

    // interinter members

    InterIntraMode_t    interintra_mode;

    /*!< Specifies the type of motion compensation to perform. */
    MotionMode         motion_mode;

    InterIntraMode     is_inter_intra;

    InterCompoundData_t inter_compound;

    FilterIntraModeInfo_t filter_intra_mode_info;

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

#if MODE_INFO_DBG
    int32_t mi_row;
    int32_t mi_col;
#endif
} ModeInfo_t;

typedef struct SBInfo {

    int8_t      *sb_cdef_strength; /*!< At 64x64 blk level */
    int32_t     *sb_delta_q; /*!< At SB level */
    int32_t     *sb_delta_lf; /*!< At SB level */

    TransformInfo_t *sb_luma_trans_info;
    TransformInfo_t *sb_chroma_trans_info;

    int32_t         *sb_luma_coeff;
    int32_t         *sb_chroma_coeff;

    ModeInfo_t      *sb_mode_info;

    //Add buffer for coeff storage

    //Add buffer for coeff storage

} SBInfo;

typedef struct PartitionInfo {

    /*!< Specifies the vertical location of the block in units of 4x4 luma samples. */
    uint16_t     mi_row;

    /*!< Specifies the horizontal location of the block in units of 4x4 luma samples. */
    uint16_t     mi_col;

    ModeInfo_t  *mi;

    SBInfo      *sb_info;

    ModeInfo_t  *left_mbmi;

    ModeInfo_t  *above_mbmi;

    ModeInfo_t  *chroma_left_mbmi;

    ModeInfo_t  *chroma_above_mbmi;

    /*!< Indicates if the information from the block above cab be used on the luma plane. */
    uint8_t up_available;

    /*!< Indicates if the information from the block to the left can be used on the luma plane. */
    uint8_t left_available;

    /*!< Indicates if the information from the block above cab be used on the chroma plane. */
    uint8_t chroma_up_available;

    /*!< Indicates if the information from the block to the left can be used on the chroma plane. */
    uint8_t chroma_left_available;

    /*!< Specifies whether chroma information is coded for this block. */
    int8_t has_chroma;

    /*!< Distance of MB away from frame edges in subpixels (1/8th pixel).  */
    int32_t mb_to_left_edge;

    int32_t mb_to_right_edge;

    int32_t mb_to_top_edge;

    int32_t mb_to_bottom_edge;

    /*!< Block Size width & height in pixels. */
    int32_t wpx[3];
     
    int32_t hpx[3];

    /*!< 1 indicates that the block is to be coded as fully lossless,
     *   0 indicates lossy coding */
    //int lossless;

    //BlockPlane plane[MAX_MB_PLANE];

    const EbWarpedMotionParams *global_motion;

    WienerInfo wiener_info[MAX_MB_PLANE];

    SgrprojInfo sgrproj_info[MAX_MB_PLANE];

    /*!< Motion vectors available in the stack */
    CandidateMv_dec ref_mv_stack[MODE_CTX_REF_FRAMES][MAX_REF_MV_STACK_SIZE];

    /*!< Represents an offset used in derivation of the input index to the cb component scaling function */
    uint16_t cb_offset[MAX_MB_PLANE];

    /*!< Holds the index for the block’s Y plane and UV plane top left samples. */
    uint16_t color_index_map_offset[2];

    /* CFL ctxt */
    void    *pv_cfl_ctxt;

} PartitionInfo_t;

#endif //EbDecBlock_h
