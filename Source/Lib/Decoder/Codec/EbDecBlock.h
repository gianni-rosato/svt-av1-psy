/*
* Copyright(c) 2019 Netflix, Inc.
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbDecBlock_h
#define EbDecBlock_h

#include "EbDefinitions.h"
#include "EbRestoration.h"
#include "EbBlockStructures.h"

#define MODE_INFO_DBG 0

static const BlockSize partition_subsize[10][BlockSizeS_ALL] = {
    {BLOCK_4X4,     BLOCK_INVALID, BLOCK_INVALID, BLOCK_8X8,     BLOCK_INVALID, BLOCK_INVALID,
     BLOCK_16X16,   BLOCK_INVALID, BLOCK_INVALID, BLOCK_32X32,   BLOCK_INVALID, BLOCK_INVALID,
     BLOCK_64X64,   BLOCK_INVALID, BLOCK_INVALID, BLOCK_128X128, BLOCK_INVALID, BLOCK_INVALID,
     BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID},
    {BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, BLOCK_8X4,    BLOCK_INVALID, BLOCK_INVALID,
     BLOCK_16X8,    BLOCK_INVALID, BLOCK_INVALID, BLOCK_32X16,  BLOCK_INVALID, BLOCK_INVALID,
     BLOCK_64X32,   BLOCK_INVALID, BLOCK_INVALID, BLOCK_128X64, BLOCK_INVALID, BLOCK_INVALID,
     BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID},
    {BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, BLOCK_4X8,    BLOCK_INVALID, BLOCK_INVALID,
     BLOCK_8X16,    BLOCK_INVALID, BLOCK_INVALID, BLOCK_16X32,  BLOCK_INVALID, BLOCK_INVALID,
     BLOCK_32X64,   BLOCK_INVALID, BLOCK_INVALID, BLOCK_64X128, BLOCK_INVALID, BLOCK_INVALID,
     BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID},
    {BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, BLOCK_4X4,    BLOCK_INVALID, BLOCK_INVALID,
     BLOCK_8X8,     BLOCK_INVALID, BLOCK_INVALID, BLOCK_16X16,  BLOCK_INVALID, BLOCK_INVALID,
     BLOCK_32X32,   BLOCK_INVALID, BLOCK_INVALID, BLOCK_64X64,  BLOCK_INVALID, BLOCK_INVALID,
     BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID},
    {BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, BLOCK_8X4,    BLOCK_INVALID, BLOCK_INVALID,
     BLOCK_16X8,    BLOCK_INVALID, BLOCK_INVALID, BLOCK_32X16,  BLOCK_INVALID, BLOCK_INVALID,
     BLOCK_64X32,   BLOCK_INVALID, BLOCK_INVALID, BLOCK_128X64, BLOCK_INVALID, BLOCK_INVALID,
     BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID},
    {BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, BLOCK_8X4,    BLOCK_INVALID, BLOCK_INVALID,
     BLOCK_16X8,    BLOCK_INVALID, BLOCK_INVALID, BLOCK_32X16,  BLOCK_INVALID, BLOCK_INVALID,
     BLOCK_64X32,   BLOCK_INVALID, BLOCK_INVALID, BLOCK_128X64, BLOCK_INVALID, BLOCK_INVALID,
     BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID},
    {BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, BLOCK_4X8,    BLOCK_INVALID, BLOCK_INVALID,
     BLOCK_8X16,    BLOCK_INVALID, BLOCK_INVALID, BLOCK_16X32,  BLOCK_INVALID, BLOCK_INVALID,
     BLOCK_32X64,   BLOCK_INVALID, BLOCK_INVALID, BLOCK_64X128, BLOCK_INVALID, BLOCK_INVALID,
     BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID},
    {BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, BLOCK_4X8,    BLOCK_INVALID, BLOCK_INVALID,
     BLOCK_8X16,    BLOCK_INVALID, BLOCK_INVALID, BLOCK_16X32,  BLOCK_INVALID, BLOCK_INVALID,
     BLOCK_32X64,   BLOCK_INVALID, BLOCK_INVALID, BLOCK_64X128, BLOCK_INVALID, BLOCK_INVALID,
     BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID},
    {BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID,
     BLOCK_16X4,    BLOCK_INVALID, BLOCK_INVALID, BLOCK_32X8,    BLOCK_INVALID, BLOCK_INVALID,
     BLOCK_64X16,   BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID,
     BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID},
    {BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID,
     BLOCK_4X16,    BLOCK_INVALID, BLOCK_INVALID, BLOCK_8X32,    BLOCK_INVALID, BLOCK_INVALID,
     BLOCK_16X64,   BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID,
     BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID}};

#define MFMV_STACK_SIZE 3
typedef struct TemporalMvRef {
    /* Motion Filed MV */
    IntMv  mf_mv0;
    int8_t ref_frame_offset;
} TemporalMvRef;

typedef struct TransformInfo {
    /*!< Specifies the transform size to be used for this TU. */
    TxSize tx_size;

    /*!< Specifies the transform type for the TU. */
    TxType txk_type;

    /*!< Code Block Flag, 0: No residual for the block
                          1: Residual exists for the block */
    int8_t cbf;

    /*!< x offset for a block in mi unit*/
    uint8_t txb_x_offset;
    /*!< y offset for a block in mi unit */
    uint8_t txb_y_offset;

} TransformInfo_t;

typedef struct SBInfo {
    int8_t * sb_cdef_strength; /*!< At 64x64 blk level */
    int32_t *sb_delta_q; /*!< At SB level */
    int32_t *sb_delta_lf; /*!< At SB level */

    TransformInfo_t *sb_trans_info[MAX_MB_PLANE - 1];

    int32_t *sb_coeff[MAX_MB_PLANE];

    BlockModeInfo *sb_mode_info;

    int32_t     num_block;

} SBInfo;

typedef struct PartitionInfo {
    /*!< Specifies the vertical location of the block in units of 4x4 luma samples. */
    uint16_t mi_row;

    /*!< Specifies the horizontal location of the block in units of 4x4 luma samples. */
    uint16_t mi_col;

    BlockModeInfo *mi;

    SBInfo *sb_info;

    BlockModeInfo *left_mbmi;

    BlockModeInfo *above_mbmi;

    BlockModeInfo *chroma_left_mbmi;

    BlockModeInfo *chroma_above_mbmi;

    /*!< Indicates if the information from the block above cab be used on the luma plane. */
    uint8_t up_available;

    /*!< Indicates if the information from the block to the left can be used on the luma plane. */
    uint8_t left_available;

    // TO-DO bhavna Verify if this is necessary. Can this info be accessed from elsewhere
    uint8_t neighbors_ref_counts[REF_FRAMES];

    /*!< Indicates if the information from the block above cab be used on the chroma plane. */
    uint8_t chroma_up_available;

    /*!< Indicates if the information from the block to the left can be used on the chroma plane. */
    uint8_t chroma_left_available;

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

    /*!< Pointer to global warp params of current frame */
    const EbWarpedMotionParams *ps_global_motion;

    /*!< Pointer to local warp params based on nieghbour mv sample projection */
    EbWarpedMotionParams local_warp_params;

    WienerInfo wiener_info[MAX_MB_PLANE];

    SgrprojInfo sgrproj_info[MAX_MB_PLANE];

    /*!< Motion vectors available in the stack */
    CandidateMv ref_mv_stack[MODE_CTX_REF_FRAMES][MAX_REF_MV_STACK_SIZE];

    /*!< Represents an offset used in derivation of the input index to the cb component scaling function */
    uint16_t cb_offset[MAX_MB_PLANE];

    /*!< Holds the index for the blocks Y plane and UV plane top left samples. */
    uint16_t color_index_map_offset[2];

    /* CFL ctxt */
    void *pv_cfl_ctxt;

    int is_sec_rect;

    int num_samples;

    /*!< chroma sub-sampling format */
    uint8_t subsampling_x;
    uint8_t subsampling_y;

    /*ToDo: block_ref_sf are used both in parsing & decoding sides,
      we need implement a logic to avoid two sides calculation of block_ref_sf*/
    /* pointers to reference frame scale factors */
    const struct ScaleFactors *block_ref_sf[2];
    const struct ScaleFactors *sf_identity;

    int8_t *cdef_strength;

    int32_t is_chroma_ref;
    /*MC temp buff for dynamic padding*/
    uint8_t *mc_buf[2];
} PartitionInfo;

#endif //EbDecBlock_h
