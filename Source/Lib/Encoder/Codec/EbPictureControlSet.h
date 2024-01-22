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

#ifndef EbPictureControlSet_h
#define EbPictureControlSet_h

#include "EbSvtAv1Enc.h"
#include "EbDefinitions.h"
#include "EbSystemResourceManager.h"
#include "EbPictureBufferDesc.h"
#include "EbCodingUnit.h"
#include "EbEntropyCodingObject.h"
#include "EbDefinitions.h"
#include "EbPredictionStructure.h"
#include "EbNeighborArrays.h"
#include "EbModeDecisionSegments.h"
#include "EbEncDecSegments.h"
#include "EbRestoration.h"
#include "EbObject.h"
#include "noise_model.h"
#include "EbSegmentationParams.h"
#include "EbAv1Structs.h"
#include "EbMdRateEstimation.h"
#include "EbEncCdef.h"
#include "Av1Common.h"

#include "av1me.h"
#include "hash_motion.h"
#include "firstpass.h"

#ifdef __cplusplus
extern "C" {
#endif
#define HISTOGRAM_NUMBER_OF_BINS 256
#define MAX_NUMBER_OF_REGIONS_IN_WIDTH 4
#define MAX_NUMBER_OF_REGIONS_IN_HEIGHT 4
#define MAX_REF_QP_NUM 81
#define QPS_SW_THRESH 8 // 100 to shut QPS/QPM (i.e. CORE only)
// BDP OFF
#define MD_NEIGHBOR_ARRAY_INDEX 0
#define MULTI_STAGE_PD_NEIGHBOR_ARRAY_INDEX 4
#define NA_TOT_CNT 5
#define AOM_QM_BITS 5

typedef struct DepCntPicInfo {
    uint64_t pic_num;
    // increase(e.g 4L->5L) or decrease of dep cnt . not including the run-time decrease
    int32_t dep_cnt_diff;
} DepCntPicInfo;
typedef struct EbDownScaledBufDescPtrArray {
    EbPictureBufferDesc *picture_ptr;
    EbPictureBufferDesc *quarter_picture_ptr;
    EbPictureBufferDesc *sixteenth_picture_ptr;
    uint64_t             picture_number;
} EbDownScaledBufDescPtrArray;

typedef struct EbDownScaledObject {
    EbDctor dctor;
    // EbPictureBufferDesc *picture_ptr; original picture, just a pointer, don't allocate resource
    // here
    EbPictureBufferDesc *quarter_picture_ptr;
    EbPictureBufferDesc *sixteenth_picture_ptr;
    // uint64_t            picture_number;
} EbDownScaledObject;

typedef struct EbDownScaledObjectDescInitData {
    EbPictureBufferDescInitData quarter_picture_desc_init_data;
    EbPictureBufferDescInitData sixteenth_picture_desc_init_data;

    // whether enable 1/4,1/16 8bit luma for in_loop global motion
    uint8_t enable_quarter_luma_input;
    uint8_t enable_sixteenth_luma_input;
} EbDownScaledObjectDescInitData;
typedef struct MacroblockPlane {
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
} MacroblockPlane;
// The Quants structure is used only for internal quantizer setup in
// av1_quantize.c.
// All of its fields use the same coefficient shift/scaling at TX.
typedef struct Quants {
    // 0: dc 1: ac 2-8: ac repeated to SIMD width
    DECLARE_ALIGNED(16, int16_t, y_quant[QINDEX_RANGE][8]);
    DECLARE_ALIGNED(16, int16_t, y_quant_shift[QINDEX_RANGE][8]);
    DECLARE_ALIGNED(16, int16_t, y_zbin[QINDEX_RANGE][8]);
    DECLARE_ALIGNED(16, int16_t, y_round[QINDEX_RANGE][8]);

    // if we want to deprecate the current use of y_quant.
    DECLARE_ALIGNED(16, int16_t, y_quant_fp[QINDEX_RANGE][8]);
    DECLARE_ALIGNED(16, int16_t, u_quant_fp[QINDEX_RANGE][8]);
    DECLARE_ALIGNED(16, int16_t, v_quant_fp[QINDEX_RANGE][8]);
    DECLARE_ALIGNED(16, int16_t, y_round_fp[QINDEX_RANGE][8]);
    DECLARE_ALIGNED(16, int16_t, u_round_fp[QINDEX_RANGE][8]);
    DECLARE_ALIGNED(16, int16_t, v_round_fp[QINDEX_RANGE][8]);

    DECLARE_ALIGNED(16, int16_t, u_quant[QINDEX_RANGE][8]);
    DECLARE_ALIGNED(16, int16_t, v_quant[QINDEX_RANGE][8]);
    DECLARE_ALIGNED(16, int16_t, u_quant_shift[QINDEX_RANGE][8]);
    DECLARE_ALIGNED(16, int16_t, v_quant_shift[QINDEX_RANGE][8]);
    DECLARE_ALIGNED(16, int16_t, u_zbin[QINDEX_RANGE][8]);
    DECLARE_ALIGNED(16, int16_t, v_zbin[QINDEX_RANGE][8]);
    DECLARE_ALIGNED(16, int16_t, u_round[QINDEX_RANGE][8]);
    DECLARE_ALIGNED(16, int16_t, v_round[QINDEX_RANGE][8]);
} Quants;

// The Dequants structure is used only for internal quantizer setup in
// av1_quantize.c.
// Fields are sufffixed according to whether or not they're expressed in
// the same coefficient shift/precision as TX or a fixed Q3 format.
typedef struct Dequants {
    DECLARE_ALIGNED(16, int16_t, y_dequant_qtx[QINDEX_RANGE][8]); // 8: SIMD width
    DECLARE_ALIGNED(16, int16_t, u_dequant_qtx[QINDEX_RANGE][8]); // 8: SIMD width
    DECLARE_ALIGNED(16, int16_t, v_dequant_qtx[QINDEX_RANGE][8]); // 8: SIMD width
} Dequants;

typedef struct MacroblockdPlane {
    // TranLow *dqcoeff;
    PlaneType    plane_type;
    int32_t      subsampling_x;
    int32_t      subsampling_y;
    struct Buf2D dst;
    int32_t      is_16bit;
    // struct Buf2D pre[2];
    // EntropyContext *above_context;
    // EntropyContext *left_context;
    //  The dequantizers below are true dequntizers used only in the
    //  dequantization process.  They have the same coefficient
    //  shift/scale as TX.
    // int16_t seg_dequant_qtx[MAX_SEGMENTS][2];
    // uint8_t *color_index_map;
    //  block size in pixels
    // uint8_t width, height;
    // QmVal *seg_iqmatrix[MAX_SEGMENTS][TX_SIZES_ALL];
    // QmVal *seg_qmatrix[MAX_SEGMENTS][TX_SIZES_ALL];
    //  the 'dequantizers' below are not literal dequantizer values.
    //  They're used by encoder RDO to generate ad-hoc lambda values.
    //  They use a hardwired Q3 coeff shift and do not necessarily match
    //  the TX scale in use.
    // const int16_t *dequant_Q3;
} MacroblockdPlane;

struct PredictionUnit;

/**************************************
 * Segment-based Control Sets
 **************************************/

typedef struct EbMdcLeafData {
    uint32_t mds_idx;
    uint32_t tot_d1_blocks; // how many d1 bloks every parent square would have
    bool     is_child; // does is it belong to the child depth(s); relative to PRED (the output of PD0)
} EbMdcLeafData;

typedef struct MdcSbData {
    uint32_t      leaf_count;
    EbMdcLeafData leaf_data_array[BLOCK_MAX_COUNT_SB_128];
    Bool          split_flag[BLOCK_MAX_COUNT_SB_128];
    uint8_t
        consider_block[BLOCK_MAX_COUNT_SB_128]; // 0: do not encode, 1: current or parent depth(s), 2: child depth(s)
    uint8_t refined_split_flag[BLOCK_MAX_COUNT_SB_128];
} MdcSbData;

/**************************************
 * Picture Control Set
 **************************************/
#define MAX_MESH_STEP 4

typedef struct MeshPattern {
    int range;
    int interval;
} MeshPattern;

typedef struct CdfControls {
    uint8_t enabled; // 1 if mv, or se, or coeff is ON
    uint8_t update_mv; // cdf update for mv
    uint8_t update_se; // cdf update for various syntax elements
    uint8_t update_coef; // cdf update for coeffs
} CdfControls;
typedef struct SpeedFeatures {
    // Threshold for allowing exhaistive motion search.
    int exhaustive_searches_thresh;
    // Pattern to be used for any exhaustive mesh searches.
    MeshPattern mesh_patterns[MAX_MESH_STEP];

} SpeedFeatures;
typedef struct EncDecSet {
    EbDctor                         dctor;
    EbPictureBufferDesc            *recon_pic;
    EbPictureBufferDesc            *recon_pic_16bit;
    EbPictureBufferDesc           **quantized_coeff;
    EbObjectWrapper                *enc_dec_wrapper;
    struct PictureParentControlSet *ppcs; // The parent of this PCS.
    EbObjectWrapper                *ppcs_wrapper;
    uint16_t                        b64_total_count;
    uint16_t                        init_b64_total_count;
} EncDecSet;
typedef struct CdefDirData {
    uint8_t dir[CDEF_NBLOCKS][CDEF_NBLOCKS];
    int32_t var[CDEF_NBLOCKS][CDEF_NBLOCKS];
} CdefDirData;

typedef struct PicVqCtrls {
    uint8_t detect_high_freq_lvl;
} PicVqCtrls;

typedef struct PictureControlSet {
    /*!< Pointer to the dtor of the struct*/
    EbDctor                    dctor;
    struct SequenceControlSet *scs;
    EbPictureBufferDesc       *input_frame16bit;

    struct PictureParentControlSet *ppcs; // The parent of this PCS.
    EbObjectWrapper                *ppcs_wrapper;
    // Packetization (used to encode SPS, PPS, etc)
    Bitstream *bitstream_ptr;

    EbObjectWrapper *c_pcs_wrapper_ptr;

    // Reference Lists
    EbObjectWrapper *ref_pic_ptr_array[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];
    // EB_S64 refPicPocArray[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];

    uint8_t   ref_pic_qp_array[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];
    SliceType ref_slice_type_array[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];
    double    ref_pic_r0[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];
    // GOP
    uint64_t         picture_number;
    uint8_t          temporal_layer_index;
    EbColorFormat    color_format;
    EncDecSegments **enc_dec_segment_ctrl;
    uint16_t         enc_dec_coded_sb_count;

    // Entropy Process Rows
    EntropyTileInfo **ec_info;
    EbHandle          entropy_coding_pic_mutex;
    Bool              entropy_coding_pic_reset_flag;
    uint8_t           tile_size_bytes_minus_1;
    EbHandle          intra_mutex;
    uint32_t          intra_coded_area;
    uint64_t          skip_coded_area;
    uint64_t          hp_coded_area;
    uint32_t          tot_seg_searched_cdef;
    EbHandle          cdef_search_mutex;

    uint16_t cdef_segments_total_count;
    uint8_t  cdef_segments_column_count;
    uint8_t  cdef_segments_row_count;

    uint64_t (*mse_seg[2])[TOTAL_STRENGTHS];
    uint8_t     *skip_cdef_seg;
    CdefDirData *cdef_dir_data;
    EbByte       cdef_input_recon[3]; // DLF'd recon
    EbByte       cdef_input_source[3]; // Input video
    uint32_t     tot_seg_searched_rest;
    EbHandle     rest_search_mutex;
    uint16_t     rest_segments_total_count;
    uint8_t      rest_segments_column_count;
    uint8_t      rest_segments_row_count;
    // flag to indicate whether the frame is extended for restoration search
    Bool rest_extend_flag[3];

    // Slice Type
    SliceType slice_type;

    // Rate Control
    uint8_t picture_qp;
    // SB Array
    uint16_t b64_total_count;
    uint16_t init_b64_total_count;
    uint16_t frame_width;
    uint16_t frame_height;

    SuperBlock **sb_ptr_array;
    uint8_t     *sb_intra;
    uint8_t     *sb_skip;
    uint8_t     *sb_64x64_mvp;
    uint32_t    *sb_count_nz_coeffs;
    // qindex per 64x64 using ME distortions (to be used for lambda modulation only; not at Q/Q-1)
    // Mode Decision Neighbor Arrays
    uint8_t            *b64_me_qindex;
    NeighborArrayUnit **md_luma_recon_na[NA_TOT_CNT];
    NeighborArrayUnit **md_tx_depth_1_luma_recon_na[NA_TOT_CNT];
    NeighborArrayUnit **md_tx_depth_2_luma_recon_na[NA_TOT_CNT];
    NeighborArrayUnit **md_cb_recon_na[NA_TOT_CNT];
    NeighborArrayUnit **md_cr_recon_na[NA_TOT_CNT];

    uint8_t             hbd_md;
    NeighborArrayUnit **md_luma_recon_na_16bit[NA_TOT_CNT];
    NeighborArrayUnit **md_tx_depth_1_luma_recon_na_16bit[NA_TOT_CNT];
    NeighborArrayUnit **md_tx_depth_2_luma_recon_na_16bit[NA_TOT_CNT];
    NeighborArrayUnit **md_cb_recon_na_16bit[NA_TOT_CNT];
    NeighborArrayUnit **md_cr_recon_na_16bit[NA_TOT_CNT];
    NeighborArrayUnit **md_y_dcs_na[NA_TOT_CNT];
    NeighborArrayUnit **md_tx_depth_1_luma_dc_sign_level_coeff_na[NA_TOT_CNT];
    NeighborArrayUnit **md_cb_dc_sign_level_coeff_na[NA_TOT_CNT];
    NeighborArrayUnit **md_cr_dc_sign_level_coeff_na[NA_TOT_CNT];
    NeighborArrayUnit **md_txfm_context_array[NA_TOT_CNT];
    NeighborArrayUnit **mdleaf_partition_na[NA_TOT_CNT];

    // Encode Pass Neighbor Arrays
    NeighborArrayUnit **ep_luma_recon_na;
    NeighborArrayUnit **ep_cb_recon_na;
    NeighborArrayUnit **ep_cr_recon_na;
    NeighborArrayUnit **ep_luma_recon_na_16bit;
    NeighborArrayUnit **ep_cb_recon_na_16bit;
    NeighborArrayUnit **ep_cr_recon_na_16bit;
    NeighborArrayUnit **ep_luma_dc_sign_level_coeff_na;
    NeighborArrayUnit **ep_cr_dc_sign_level_coeff_na;
    NeighborArrayUnit **ep_cb_dc_sign_level_coeff_na;
    NeighborArrayUnit **ep_luma_dc_sign_level_coeff_na_update;
    NeighborArrayUnit **ep_cr_dc_sign_level_coeff_na_update;
    NeighborArrayUnit **ep_cb_dc_sign_level_coeff_na_update;
    NeighborArrayUnit **ep_partition_context_na;
    NeighborArrayUnit **ep_txfm_context_na;

    // Entropy Coding Neighbor Arrays
    NeighborArrayUnit **partition_context_na;
    // Stored per 4x4. 8 bit: lower 6 bits (COEFF_CONTEXT_BITS), shows if there is at least one
    // Coef. Top 2 bit store the sign of DC as follow: 0->0,1->-1,2-> 1
    NeighborArrayUnit **luma_dc_sign_level_coeff_na;
    // Stored per 4x4. 8 bit: lower 6 bits(COEFF_CONTEXT_BITS), shows if there is at least one Coef.
    // Top 2 bit store the sign of DC as follow: 0->0,1->-1,2-> 1
    NeighborArrayUnit **cr_dc_sign_level_coeff_na;
    // Stored per 4x4. 8 bit: lower 6 bits(COEFF_CONTEXT_BITS), shows if there is at least one Coef.
    // Top 2 bit store the sign of DC as follow: 0->0,1->-1,2-> 1
    NeighborArrayUnit **cb_dc_sign_level_coeff_na;
    NeighborArrayUnit **txfm_context_array;

    NeighborArrayUnit      **segmentation_id_pred_array;
    SegmentationNeighborMap *segmentation_neighbor_map;

    ModeInfo **mi_grid_base; // 2 SB Rows of mi Data are enough

    ModeInfo *mip;

    int32_t mi_stride;
    // true if 4x4 blocks are disallowed for all frames, and NSQ is disabled (since granularity is
    // needed for 8x8 NSQ blocks).  Used to compute the offset for mip.
    uint8_t disallow_4x4_all_frames;
    uint8_t wm_level; // warped motion level
    uint8_t cand_reduction_level;
    uint8_t nsq_level;
    uint8_t txt_level;
    uint8_t tx_shortcut_level;
    uint8_t interpolation_search_level;
    uint8_t chroma_level;
    uint8_t cfl_level;
    uint8_t new_nearest_near_comb_injection;
    uint8_t unipred3x3_injection;
    uint8_t bipred3x3_injection;
    uint8_t inter_compound_mode;
    uint8_t dist_based_ref_pruning;
    uint8_t spatial_sse_full_loop_level;
    uint8_t inter_intra_level;
    uint8_t txs_level;
    uint8_t nic_level;
    uint8_t md_sq_mv_search_level;
    uint8_t md_nsq_mv_search_level;
    uint8_t md_pme_level;
    uint8_t mds0_level;
    uint8_t pic_disallow_4x4; // disallow 4x4 at pic level
    uint8_t pic_skip_pd0; // skip_pd0 at pic level
    uint8_t pic_disallow_below_16x16; // disallow_below_16x16 signal at pic level
    // depth_removal_level signal at the picture level
    uint8_t pic_depth_removal_level;
    uint8_t pic_depth_removal_level_rtc;
    // block_based_depth_refinement_level signal set at the picture level
    uint8_t          pic_block_based_depth_refinement_level;
    uint8_t          pic_lpd0_lvl; // lpd0_lvl signal set at the picture level
    uint8_t          pic_lpd1_lvl; // lpd1_lvl signal set at the picture level
    Bool             pic_bypass_encdec;
    EncMode          enc_mode;
    InputCoeffLvl    coeff_lvl;
    int32_t          cdef_preset[MAX_TILE_CNTS][4];
    WienerInfo       wiener_info[MAX_TILE_CNTS][MAX_MB_PLANE];
    SgrprojInfo      sgrproj_info[MAX_TILE_CNTS][MAX_MB_PLANE];
    SpeedFeatures    sf;
    SearchSiteConfig ss_cfg; // CHKN this might be a seq based
    HashTable        hash_table;
    CRC_CALCULATOR   crc_calculator1;
    CRC_CALCULATOR   crc_calculator2;

    FRAME_CONTEXT                  *ec_ctx_array;
    FRAME_CONTEXT                   md_frame_context;
    CdfControls                     cdf_ctrl;
    FRAME_CONTEXT                   ref_frame_context[REF_FRAMES];
    EbWarpedMotionParams            ref_global_motion[TOTAL_REFS_PER_FRAME];
    struct MdRateEstimationContext *md_rate_est_ctx;
    int8_t                          ref_frame_side[REF_FRAMES];
    TPL_MV_REF                     *tpl_mvs;
    uint8_t                         pic_filter_intra_level;
    TOKENEXTRA                     *tile_tok[64][64];
    // Put it here for deinit, don't need to go pcs->ppcs->av1_cm which may already be released
    uint16_t tile_row_count;
    uint16_t tile_column_count;
    uint16_t sb_total_count;
    // Total SB count of unscaled picture (used for memory alloc/dealloc when superres is used)
    uint16_t sb_total_count_unscaled;

    EbPictureBufferDesc *temp_lf_recon_pic;
    EbPictureBufferDesc *temp_lf_recon_pic_16bit;

    RestUnitSearchInfo *rusi_picture[3]; // for 3 planes

    RestorationInfo rst_info[MAX_MB_PLANE];
    // rst_end_stripe[i] is one more than the index of the bottom stripe
    // for tile row i.
    int32_t  rst_end_stripe[MAX_TILE_ROWS];
    uint8_t  ref_intra_percentage;
    uint8_t  ref_skip_percentage;
    int16_t  ref_hp_percentage;
    uint64_t avg_me_clpx;
    uint64_t min_me_clpx;
    uint64_t max_me_clpx;
    // use approximate rate for inter cost (set at pic-level b/c some pic-level initializations will
    // be removed)
    uint8_t    approx_inter_rate;
    uint8_t    skip_intra;
    PicVqCtrls vq_ctrls;
    // scaled input picture is only used in loop restoration for recon size is
    // different with input frame when reference scaling is enabled
    EbPictureBufferDesc *scaled_input_pic;
    bool                 rtc_tune;
} PictureControlSet;

// To optimize based on the max input size
// To study speed-memory trade-offs
typedef struct B64Geom {
    uint8_t  horizontal_index;
    uint8_t  vertical_index;
    uint16_t org_x;
    uint16_t org_y;
    uint8_t  width;
    uint8_t  height;
    uint8_t  is_complete_b64;
    Bool     raster_scan_blk_validity[CU_MAX_COUNT];
    uint8_t  is_edge_sb;
    uint32_t tile_start_x;
    uint32_t tile_start_y;
    uint32_t tile_end_x;
    uint32_t tile_end_y;
} B64Geom;

typedef struct SbGeom {
    uint16_t horizontal_index;
    uint16_t vertical_index;
    uint16_t org_x;
    uint16_t org_y;
    uint8_t  width;
    uint8_t  height;
    uint8_t  is_complete_sb;
    Bool     block_is_allowed[BLOCK_MAX_COUNT_SB_128];
} SbGeom;

typedef struct TileGroupInfo {
    uint16_t tile_group_sb_start_x;
    uint16_t tile_group_sb_start_y;
    uint16_t tile_group_sb_end_x;
    uint16_t tile_group_sb_end_y;
    uint16_t tile_group_width_in_sb;
    uint16_t tile_group_height_in_sb;

    uint16_t tile_group_tile_start_x;
    uint16_t tile_group_tile_start_y;
    uint16_t tile_group_tile_end_x;
    uint16_t tile_group_tile_end_y;
} TileGroupInfo;
typedef struct MotionEstimationData {
    EbDctor        dctor;
    MeSbResults  **me_results;
    uint16_t       b64_total_count;
    uint16_t       init_b64_total_count;
    uint8_t        max_cand; // total max me candidates given the active references
    uint8_t        max_refs; // total max active references
    uint8_t        max_l0; // max active refs in L0
    OisMbResults **ois_mb_results;
    TplStats     **tpl_stats;

    TplSrcStats *tpl_src_stats_buffer; // tpl src based stats

    int32_t base_rdmult;
    double *tpl_beta;
    double *tpl_rdmult_scaling_factors;
    double *tpl_sb_rdmult_scaling_factors;
    double *ssim_rdmult_scaling_factors;
} MotionEstimationData;
typedef struct TplControls {
    uint8_t              enable; // 0: TPL OFF; 1: TPL ON
    uint8_t              compute_rate; // 1: use rate 1: no rate
    uint8_t              enable_tpl_qps; // 0:OFF 1:ON - QPS in TPL
    uint8_t              disable_intra_pred_nref; // 0:OFF 1:ON - Disable intra prediction in NREF
    PredictionMode       intra_mode_end; // The MAX intra mode to be tested in TPL
    EB_TRANS_COEFF_SHAPE pf_shape;
    // Use SAD as a distortion metric when searching the best mode (based on src pic). If false, will use SATD
    uint8_t use_sad_in_src_search;
    int8_t  reduced_tpl_group;
    double  r0_adjust_factor;
    // 0: use 16x16 block(s), 1: use 32x32 block(s), 2: use 64x64 block(s)  (for incomplete 64x64,
    // dispenser_search_level is set to 0)
    uint8_t dispenser_search_level;
    // it is recommended to use subsample_tx=2, when dispenser_search_level is set to 1
    // 0: OFF, use full TX size; 1: subsample the transforms in TPL by 2; 2: subsample the
    // transforms in TPL by 4
    uint8_t subsample_tx;
    // syntheszier block size, support 8x8 and 16x16 for now. NOTE: this field must be modified
    // inside the get_ function, as it is linked to memory allocation at init time
    uint8_t synth_blk_size;
    // max subpel depth to search for TPL; FULL_PEL corresponds to subpel off in TPL, QUARTER_PEL is
    // the max precision for TPL subpel
    SUBPEL_FORCE_STOP
    subpel_depth;
} TplControls;

typedef struct {
    uint8_t   tpl_temporal_layer_index;
    SliceType tpl_slice_type;
    uint8_t   tpl_ref0_count;
    uint8_t   tpl_ref1_count;
    uint64_t  tpl_decode_order;
    Bool      ref_in_slide_window[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];
    // tpl group index of all ref pictures; -1 if the reference is not in the TPL group (aka sliding
    // window)
    int32_t ref_tpl_group_idx[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];
    // TPL base picture
    struct PictureParentControlSet *base_pcs;
    Bool                            is_ref;
    EbDownScaledBufDescPtrArray     tpl_ref_ds_ptr_array[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];
} TPLData;
typedef struct GmControls {
    // 0: generate GM params for both list_0 and list_1, 1: do not generate GM params for list_1 if
    // list_0/ref_idx_0 is id
    uint8_t enabled;
    // 0: use both rotzoom and affine models, 1:use rotzoom model only
    uint8_t identiy_exit;
    // 0: test both unipred and bipred, 1: test bipred only
    uint8_t rotzoom_model_only;
    uint8_t bipred_only;
    // 0: do not consider stationary_block info @ me-based bypass, 1: consider stationary_block info
    // @ me-based bypass (only if bypass_based_on_me=1)
    uint8_t bypass_based_on_me;
    // 0: used default active_th,1: increase active_th baed on distance to ref (only if
    // bypass_based_on_me=1)
    uint8_t use_stationary_block;
    // The number of refinement steps to use in the GM params refinement
    uint8_t use_distance_based_active_th;
    // GM_FULL: Exhaustive search mode; GM_DOWN: Downsampled resolution with a downsampling factor
    // of 2 in each dimension; GM_TRAN_ONLY: Translation only using ME MV.
    uint8_t params_refinement_steps;
    uint8_t downsample_level;
    // do GM in the closed loop instead of the open loop and use reference information 0: off 1: on
    bool use_ref_info;
    // do the detection bypass for last layer pictures   0:off     1:last layer     2:last 2 layers 3: last 3 layers
    uint8_t layer_offset;
    //use a fraction of corner points for computing correspondences for RANSAC in detection. 1:1/4   2:2/4   3:3/4   4:all
    uint8_t corners;
    //skip global motion refinement using a chess pattern to skip blocks
    uint8_t chess_rfn;
    //change the window size for correlation calculations. must be odd. N: NxN window size goes from 1 to 15
    uint8_t match_sz;
    //Inject global only if Parent SQ is global
    bool inj_psq_glb;
    //enable Pre-processor for GM
    bool pp_enabled;
    //limit the search to  ref index = 0 only
    bool ref_idx0_only;
    // if true, apply an offset to the segments of the me-dist based modulation
    uint8_t qp_offset;
} GmControls;
typedef struct CdefControls {
    uint8_t enabled;
    uint8_t number_of_prim_in_second_loop[2];
    // Number of primary filters considered in the first pass. (luma and chroma)
    uint8_t first_pass_fs_num;
    // Primary filter strengths to consider in the first pass.
    uint8_t default_first_pass_fs[TOTAL_STRENGTHS];
    // Number of secondary filters considered in the second pass. (luma and chroma)
    uint8_t default_second_pass_fs_num;
    // Secondary filter strengths to consider in the second pass.
    uint8_t default_second_pass_fs[TOTAL_STRENGTHS];
    // Mask for primary filters to be considered for chroma and indicates a subset of the primary
    // filter strengths considered in default_first_pass_fs[64]
    int8_t default_first_pass_fs_uv[TOTAL_STRENGTHS];
    // Mask for secondary filters to be considered for chroma and indicates a subset of the
    // secondary filter strengths considered in default_second_pass_fs[64]
    int8_t default_second_pass_fs_uv[TOTAL_STRENGTHS];
    // Flag to indicate the use of reference frames' filter strengths.
    int8_t use_reference_cdef_fs;
    // Predicted filter strength pair index for the luma component based on reference picture filter
    // strength pairs.
    int8_t pred_y_f;
    // Predicted filter strength pair index for the chroma component based on reference picture
    // filter strength pairs.
    int8_t pred_uv_f;
    // Allowable levels: [1,2,4] ---> 1: no subsampling; 2: process every 2nd row; 4: process every
    // 4th row for 8x8 blocks, every 2nd row for smaller sizes.
    uint8_t subsampling_factor;
    // NB subsampling is capped for certain block sizes, based on how many points the intrinsics can
    // process at once. Only search best filter strengths of the nearest ref frames (skips the
    // search if the filters of list0/list1 are the same).
    uint8_t search_best_ref_fs;
    // 0: OFF, higher is safer. Scaling factor to decrease the zero filter strength cost: : <x>/64
    uint16_t zero_fs_cost_bias;
    // Shut CDEF at the picture level based on the skip area of the nearest reference frames.
    uint8_t use_skip_detector;
} CdefControls;
typedef struct DlfCtrls {
    uint8_t enabled; // if true, perform DLF per SB, not per picture
    // when DLF filter level is selected from QP, if the filter level is less than or equal to this
    // TH, the filter level is set to 0
    uint8_t sb_based_dlf;
    // Start search from average DLF instead of 0
    bool dlf_avg;
    // Use average DLF as a starting point for qp based filter strength selection for Chroma planes
    bool dlf_avg_uv;
    // Number of convergence points before exiting the filter search, 1 = exit on first convergence point, 2 = exit on second, 0 = off
    uint8_t early_exit_convergence;
    // Threshold used when sb_based_dlf is used to use filter strength zero, there are four levels of thresholds [0..3], 0 = off
    uint8_t zero_filter_strength_lvl;
} DlfCtrls;
typedef struct IntraBCCtrls {
    // Shift for full_pixel_exhaustive search threshold:   0: No Shift   1:Shift to left by 1
    uint8_t enabled;
    uint8_t ibc_shift;
    // search directions: 0: Left + Top      1: Top only
    uint8_t ibc_direction;
    // this will be set by get_disallow_4x4, will not hash 4x4 blocks for higher presets where 4x4
    // blocks are not encoded
    uint8_t hash_4x4_blocks;
    // the maximum block size that will be hashed, corresponds to the maximum block size for which a
    // MD candidate will be generated by IBC hashing algorithm
    uint8_t max_block_size_hash;
} IntraBCCtrls;
typedef struct PaletteCtrls {
    uint8_t enabled;
    // In the dominant color search, test a subset of the most dominant color combinations
    // by testing every nth combo. E.g. with step size of 2, if the block involves 7 colors, then
    // only 3 candidates with palettes based on the most dominant 7, 5 and 3 colors are tested.
    // Range: [1 (test all), 7 (test one)]
    uint8_t dominant_color_step;
    uint8_t reduce_palette_cost_precision;
} PaletteCtrls;

/*!
 * \brief The stucture of Cyclic_Refresh.
 * \ingroup cyclic_refresh
 */
typedef struct CyclicRefresh {
    /*!
     * Percentage of blocks per frame that are targeted as candidates
     * for cyclic refresh.
     */
    int percent_refresh;
    /*!
     * Maximum q-delta as percentage of base q.
     */
    int max_qdelta_perc;
    /*!
     *Superblock starting index for cycling through the frame.
     */
    uint32_t sb_start;
    /*!
     *Superblock end index for cycling through the frame.
     */
    uint32_t sb_end;
    /*!
     * Rate target ratio to set q delta.
     */
    double rate_ratio_qdelta;

    int apply_cyclic_refresh;

} CyclicRefresh;
// struct stores the metrics used by the dynamic gop detector
typedef struct DGDetectorMetrics {
    uint64_t tot_dist;
    uint32_t tot_cplx;
    uint32_t tot_active;
    int      sum_in_vectors;
    uint16_t seg_completed; // number of dynamic gop detector segments completed
} DGDetectorMetrics;
// struct stores the control structure for the DG segments
typedef struct DGDetectorSeg {
    EbDctor dctor;
    // 2nd picture used by the dynamic gop detector in order to determine whether to split the GoP
    struct PictureParentControlSet *ref_pic;
    struct DGDetectorMetrics        metrics;
    // semaphore that indicates whether all of the segments of the dynamic gop detector are finished
    EbHandle frame_done_sem;
    // ensures that only one dynamic gop detector segment is modifying the dg detector metrics at any time
    EbHandle metrics_mutex;
} DGDetectorSeg;
// CHKN
//  Add the concept of PictureParentControlSet which is a subset of the old PictureControlSet.
//  It actually holds only high level Picture based control data:(GOP management,when to start a
//  picture, when to release the PCS, ....). The regular PictureControlSet(Child) will be dedicated
//  to store SB based encoding results and information. Parent is created before the Child, and
//  continue to live more. Child PCS only lives the exact time needed to encode the picture: from ME
//  to EC/ALF.
typedef struct PictureParentControlSet {
    EbDctor          dctor;
    EbObjectWrapper *input_pic_wrapper;
    // when overlay: y8b buffer is not used, should set to NULL.
    EbObjectWrapper     *y8b_wrapper;
    EbObjectWrapper     *ref_pic_wrapper;
    EbObjectWrapper     *pa_ref_pic_wrapper;
    EbPictureBufferDesc *enhanced_pic;
    EbPictureBufferDesc *enhanced_downscaled_pic;
    EbPictureBufferDesc *enhanced_unscaled_pic;
    // if 422/444 input, down sample to 420 for MD
    EbPictureBufferDesc       *chroma_downsampled_pic;
    Bool                       is_chroma_downsampled_picture_ptr_owner;
    PredictionStructure       *pred_struct_ptr; // need to check
    struct SequenceControlSet *scs;
    EbObjectWrapper           *scs_wrapper;
    EbObjectWrapper           *p_pcs_wrapper_ptr;
    EbObjectWrapper           *previous_picture_control_set_wrapper_ptr;
    EbObjectWrapper           *output_stream_wrapper_ptr;
    Av1Common                 *av1_cm;

    uint8_t hbd_md;
    // Data attached to the picture. This includes data passed from the application, or other data
    // the encoder attaches to the picture.
    EbLinkedListNode *data_ll_head_ptr;
    // pointer to data to be passed back to the application when picture encoding is done
    EbLinkedListNode *app_out_data_ll_head_ptr;

    EbBufferHeaderType *input_ptr; // input picture buffer
    uint8_t             log2_tile_rows;
    uint8_t             log2_tile_cols;
    uint8_t             log2_sb_size;
    TileGroupInfo      *tile_group_info;
    uint8_t             tile_group_cols;
    uint8_t             tile_group_rows;

    Bool      idr_flag;
    Bool      cra_flag;
    Bool      scene_change_flag;
    int8_t    transition_present; // -1: not computed
    Bool      end_of_sequence_flag;
    uint8_t   picture_qp;
    uint64_t  picture_number;
    uint32_t  cur_order_hint;
    uint32_t  ref_order_hint[INTER_REFS_PER_FRAME];
    SliceType slice_type;
    uint8_t   pred_struct_index;
    uint8_t   temporal_layer_index;
    uint64_t  decode_order;

    bool similar_brightness_refs; //whether closest references have similar brightness

    //avg luma intensity of the picture  256: invalid value  0..255 valid value
    uint64_t avg_luma;

    // Each picture can release up to 8 references from the DPB (8 is the max number of entries in
    // the DPB). Each frame may also release itself at EOS, which is only done in the encoder to
    // satisfy CI unit tests for MacOS.
    uint64_t released_pics[REF_FRAMES + 1];
    uint8_t  released_pics_count;
    Bool     is_ref;
    bool     is_highest_layer;
    // status of PA reference 0: Not release; 1: Released
    uint8_t reference_released;
    uint8_t ref_list0_count;
    uint8_t ref_list1_count;
    // The number of references to try (in ME / MD) in list0.Should be <= ref_list0_count.
    uint8_t ref_list0_count_try;
    // The number of references to try (in ME/MD) in list1. Should be <= ref_list1_count.
    uint8_t          ref_list1_count_try;
    MvReferenceFrame ref_frame_type_arr[MODE_CTX_REF_FRAMES];
    uint8_t          tot_ref_frame_types;
    // Rate Control
    uint64_t                                total_num_bits;
    uint16_t                                b64_total_count;
    Bool                                    end_of_sequence_region;
    int                                     frames_in_sw; // used for Look ahead
    struct RateControlIntervalParamContext *rate_control_param_ptr;
    Bool                                    qp_on_the_fly;
    uint64_t                                last_idr_picture;
    uint64_t                                start_time_seconds;
    uint64_t                                start_time_u_seconds;
    uint64_t                                luma_sse;
    uint64_t                                cr_sse;
    uint64_t                                cb_sse;
    double                                  luma_ssim;
    double                                  cr_ssim;
    double                                  cb_ssim;
    EbPictureBufferDesc                    *quarter_src_pic;
    EbPictureBufferDesc                    *sixteenth_src_pic;
    // Pointer array for down scaled pictures
    EbObjectWrapper            *downscaled_pic_wrapper;
    EbDownScaledBufDescPtrArray ds_pics;

    TPLData          tpl_data;
    EbObjectWrapper *ref_y8b_array[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];
    // Pre Analysis
    EbObjectWrapper *ref_pa_pic_ptr_array[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];
    uint64_t         ref_pic_poc_array[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];
    uint16_t       **variance;
    uint32_t         pre_assignment_buffer_count;
    uint16_t         pic_avg_variance;

    uint32_t ***picture_histogram;
    uint64_t    average_intensity_per_region[MAX_NUMBER_OF_REGIONS_IN_WIDTH][MAX_NUMBER_OF_REGIONS_IN_HEIGHT];
    // Segments
    uint16_t me_segments_total_count;
    uint8_t  me_segments_column_count;
    uint8_t  me_segments_row_count;
    uint16_t me_segments_completion_count;

    EbPictureBufferDesc *non_tf_input;
    // Motion Estimation Results
    uint8_t   max_number_of_pus_per_sb;
    uint32_t *rc_me_distortion;
    // 1 when a % of the SB is stationary relative to reference frame(s) ((0,0) MV: decode order), 0
    // otherwise
    uint8_t *stationary_block_present_sb;
    uint8_t *rc_me_allow_gm;

    uint32_t *me_8x8_cost_variance;
    uint32_t *me_64x64_distortion;
    uint32_t *me_32x32_distortion;
    uint32_t *me_16x16_distortion;
    uint32_t *me_8x8_distortion;
    // Global motion estimation results
    Bool                 is_global_motion[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];
    EbWarpedMotionParams svt_aom_global_motion_estimation[MAX_NUM_OF_REF_PIC_LIST][REF_LIST_MAX_DEPTH];
    int8_t               is_gm_on; //-1 invalid, 1: gm on in one of the ref frames,  0:gm off for all ref frames
    uint16_t             me_processed_b64_count;
    EbHandle             me_processed_b64_mutex;
    FirstPassData        firstpass_data;
    double               ts_duration;
    double               r0;
    // track pictures that are processd in two different TPL groups
    uint8_t tpl_src_data_ready;
    Bool    blk_lambda_tuning;
    // Dynamic GOP
    SvtAv1PredStructure pred_structure;
    uint8_t             hierarchical_levels;
    Bool                init_pred_struct_position_flag;
    int8_t              hierarchical_layers_diff;
    // HME Flags
    Bool enable_hme_flag;
    Bool enable_hme_level0_flag;
    Bool enable_hme_level1_flag;
    Bool enable_hme_level2_flag;

    // HME Flags form Temporal Filtering
    Bool tf_enable_hme_flag;
    Bool tf_enable_hme_level0_flag;
    Bool tf_enable_hme_level1_flag;
    Bool tf_enable_hme_level2_flag;

    // MD
    EncMode enc_mode;
    // Multi-modes signal(s)
    MultiPassPdLevel multi_pass_pd_level;
    DlfCtrls         dlf_ctrls;
    uint8_t          intra_pred_mode;
    uint8_t          frame_end_cdf_update_mode; // mm-signal: 0: OFF, 1:ON
    //**********************************************************************************************************//
    Av1RpsNode av1_ref_signal;
    Bool       has_show_existing;
    /* maps fb_idx to reference slot */
    int32_t ref_frame_map[REF_FRAMES];
    // index of picture in the mg
    uint32_t pic_idx_in_mg;
    // Flag for a frame used as a reference - not written to the Bitstream
    int32_t is_reference_frame;

    // Flag signaling that the frame is encoded using only INTRA modes.
    uint8_t intra_only;
    /* profile settings */
#if CONFIG_ENTROPY_STATS
    int32_t coef_cdf_category;
#endif

    // Global quant matrix tables
    const QmVal *giqmatrix[NUM_QM_LEVELS][3][TX_SIZES_ALL];
    const QmVal *gqmatrix[NUM_QM_LEVELS][3][TX_SIZES_ALL];
    int32_t      min_qmlevel;
    int32_t      max_qmlevel;
    // Encoder
    LoopFilterInfoN lf_info;

    // Flag signaling how frame contexts should be updated at the end of
    // a frame decode
    RefreshFrameContextMode refresh_frame_context;
    uint32_t                frame_context_idx; /* Context to use/update */
    int32_t                 fb_of_context_type[REF_FRAMES];
    uint64_t                frame_offset;
    uint32_t                large_scale_tile;
    int32_t                 nb_cdef_strengths;
    // Resolution of delta quant
    int32_t monochrome;
    int32_t prev_qindex[MAX_TILE_CNTS];

    // Since actual frame level loop filtering level value is not available
    // at the beginning of the tile (only available during actual filtering)
    // at encoder side.we record the delta_lf (against the frame level loop
    // filtering level) and code the delta between previous superblock's delta
    // lf and current delta lf. It is equivalent to the delta between previous
    // superblock's actual lf and current lf.
    int32_t prev_delta_lf_from_base;
    // For this experiment, we have four frame filter levels for different plane
    // and direction. So, to support the per superblock update, we need to add
    // a few more params as below.
    // 0: delta loop filter level for y plane vertical
    // 1: delta loop filter level for y plane horizontal
    // 2: delta loop filter level for u plane
    // 3: delta loop filter level for v plane
    // To make it consistent with the reference to each filter level in segment,
    // we need to -1, since
    // SEG_LVL_ALT_LF_Y_V = 1;
    // SEG_LVL_ALT_LF_Y_H = 2;
    // SEG_LVL_ALT_LF_U   = 3;
    // SEG_LVL_ALT_LF_V   = 4;
    //
    int32_t prev_delta_lf[FRAME_LF_COUNT];
    int32_t curr_delta_lf[FRAME_LF_COUNT];

    // Resolution of delta quant
    // int32_t delta_q_res;
    int32_t              allow_comp_inter_inter;
    int16_t              pan_mvx;
    int16_t              pan_mvy;
    int16_t              tilt_mvx;
    int16_t              tilt_mvy;
    EbWarpedMotionParams global_motion[TOTAL_REFS_PER_FRAME];
    PictureControlSet   *child_pcs;
    EncDecSet           *enc_dec_ptr;
    Macroblock          *av1x;
    int32_t              film_grain_params_present; // todo (AN): Do we need this flag at picture level?
    int8_t               cdef_level;
    Bool                 enable_restoration; // TRUE if restoration filtering is enabled for the current frame
    uint8_t              palette_level;
    uint8_t              sc_class0;
    uint8_t              sc_class1;
    uint8_t              sc_class2;
    SkipModeInfo         skip_mode_info;
    uint64_t             picture_number_alt; // The picture number overlay includes all the overlay frames
    uint8_t              is_alt_ref;
    uint8_t              is_overlay;
    struct PictureParentControlSet *overlay_ppcs_ptr;
    struct PictureParentControlSet *alt_ref_ppcs_ptr;
    int32_t                         noise_levels_log1p_fp16[MAX_MB_PLANE];
    int32_t                         pic_decision_reorder_queue_idx;
    struct PictureParentControlSet *temp_filt_pcs_list[ALTREF_MAX_NFRAMES];
    EbByte                          save_source_picture_ptr[3];
    EbByte                          save_source_picture_bit_inc_ptr[3];
    uint16_t                        save_source_picture_width;
    uint16_t                        save_source_picture_height;
    EbHandle                        temp_filt_done_semaphore;
    EbHandle                        temp_filt_mutex;
    EbHandle                        debug_mutex;

    uint8_t  temp_filt_prep_done;
    uint16_t temp_filt_seg_acc;
    EbHandle tpl_disp_done_semaphore;

    AtomicVarU32 pa_me_done; // set when PA ME is done.
    CondVar      me_ready;

    int16_t     tf_segments_total_count;
    uint8_t     tf_segments_column_count;
    uint8_t     tf_segments_row_count;
    uint8_t     past_altref_nframes;
    uint8_t     future_altref_nframes;
    Bool        do_tf;
    uint64_t    filtered_sse_uv;
    FrameHeader frm_hdr;
    uint16_t   *altref_buffer_highbd[3];
    uint8_t     pic_obmc_level;

    Bool is_pcs_sb_params;
    // Contains info for 64x64 blocks, NOT SB info.  Should not be used in MD; use sb_geom instead.
    B64Geom *b64_geom;
    // Contains SB info (adapts automatically when SB 128x128 is used). This should be used in MD.
    SbGeom           *sb_geom;
    EbInputResolution input_resolution;
    uint16_t          picture_sb_width;
    uint16_t          picture_sb_height;

    // Picture dimensions (resized or not)
    // aligned to be a multiple of 8 pixels
    uint16_t aligned_width;
    uint16_t aligned_height;

    // Picture dimensions (resized or not)
    // --not-- aligned to be a multiple of 8 pixels
    uint16_t frame_width;
    uint16_t frame_height;

    // Picture dimensions (resized)
    // --not-- aligned to be a multiple of 8 pixels
    uint16_t render_width;
    uint16_t render_height;

    Bool    frame_resize_enabled;
    uint8_t resize_denom;
    // reference scaling random access event
    EbRefFrameScale resize_evt;
    Bool            rc_reset_flag;

    Bool    frame_superres_enabled;
    uint8_t superres_denom;
    // recode for auto superres
    // which loop is now running, range from 0 to superres_total_recode_loop - 1
    int32_t superres_recode_loop;
    int32_t superres_total_recode_loop; // how many loops to run, set to 2 in dual search mode
    uint8_t superres_denom_array[NUM_SR_SCALES + 1]; // denom candidate array used in auto supreres
    double  superres_rdcost[NUM_SR_SCALES + 1]; // 9 slots, for denom 8 ~ 16

    EbObjectWrapper      *me_data_wrapper;
    MotionEstimationData *pa_me_data;
    // stores pcs pictures needed for tpl algorithm
    struct PictureParentControlSet *tpl_group[MAX_TPL_GROUP_SIZE];
    // size of above buffer
    uint32_t tpl_group_size;
    // stores previous, current, future pictures from pd-reord-queue. empty for first I.
    // 1 past + 1 curent  + TF_MAX_BASE_REF_PICS from future
    void *pd_window[2 + TF_MAX_BASE_REF_PICS];
    // stores pcs pictures needed for lad mg based algorithms
    struct PictureParentControlSet *ext_group[MAX_TPL_EXT_GROUP_SIZE];
    // actual size of extended group
    uint32_t ext_group_size;

    int64_t ext_mg_id;
    int64_t ext_mg_size; // same as mg expect for MGops with [LDP-I] which are split into 2
    uint8_t tpl_valid_pic[MAX_TPL_EXT_GROUP_SIZE];
    uint8_t used_tpl_frame_num;

    // Tune TPL for better chroma.Only for 240P
    uint8_t      tune_tpl_for_chroma;
    uint8_t      is_not_scaled;
    TfControls   tf_ctrls;
    GmControls   gm_ctrls;
    GM_LEVEL     gm_downsample_level;
    bool         gm_pp_enabled;
    bool         gm_pp_detected; //gm detection enabled at the pre-processing level
    CdefControls cdef_ctrls;
    // RC related variables
    int                             q_low;
    int                             q_high;
    int                             loop_count;
    int                             overshoot_seen;
    int                             undershoot_seen;
    int                             low_cr_seen;
    uint64_t                        pcs_total_rate;
    EbHandle                        pcs_total_rate_mutex;
    int16_t                         first_pass_seg_total_count;
    uint8_t                         first_pass_seg_column_count;
    uint8_t                         first_pass_seg_row_count;
    uint16_t                        first_pass_seg_acc;
    EbHandle                        first_pass_done_semaphore;
    EbHandle                        first_pass_mutex;
    struct PictureParentControlSet *first_pass_ref_ppcs_ptr[2];
    uint8_t                         first_pass_ref_count;
    uint8_t                         first_pass_done;
    uint8_t                         first_frame_in_minigop;
    TplControls                     tpl_ctrls;
    uint8_t                         tpl_is_valid;
    EbHandle                        tpl_disp_mutex;
    // uint32_t         input_type;
    int16_t  enc_dec_segment_row;
    uint16_t tile_group_index;
    uint16_t tpl_disp_coded_sb_count;

    EncDecSegments **tpl_disp_segment_ctrl;
    // the offsets for STATS_BUFFER_CTX
    uint64_t stats_in_end_offset;
    // the offsets for stats_in
    uint64_t stats_in_offset;
    // GF_GROUP parameters store in PCS
    SvtAv1FrameUpdateType update_type;
    int                   layer_depth;
    // RATE_CONTROL parameters store in PCS
    int base_frame_target; // A baseline frame target before adjustment.
    int this_frame_target; // Actual frame target after rc adjustment.
    int projected_frame_size;
    int max_frame_size;
    int frames_to_key;
    int frames_since_key;
    int top_index;
    int bottom_index;
    // stores gf group (minigop) length
    int gf_interval;
    // thr gf update in RC is due for I, or base frames (except the one after I)or P frames
    int                             gf_update_due;
    uint8_t                         is_new_gf_group;
    struct PictureParentControlSet *gf_group[MAX_TPL_GROUP_SIZE];
    StatStruct                      stat_struct;
    CyclicRefresh                   cyclic_refresh;
    bool                            ld_enhanced_base_frame; // enhanced periodic base layer frames used in LD
    bool                            update_ref_count; // Update ref count
    bool                            use_accurate_part_ctx;
    uint8_t                         bypass_cost_table_gen;
    uint16_t                        max_can_count;
    uint8_t                         enable_me_8x8;
    uint8_t                         enable_me_16x16;
    uint8_t                         use_best_me_unipred_cand_only; // if MRP is OFF, use one ME unipred candidate only
    IntraBCCtrls                    intraBC_ctrls;
    PaletteCtrls                    palette_ctrls;

    uint32_t         tf_tot_vert_blks; // total vertical motion blocks in TF
    uint32_t         tf_tot_horz_blks; // total horizontal motion blocks in TF
    int8_t           tf_motion_direction; // motion direction in TF   -1:invalid   0:horz  1:vert
    int32_t          is_noise_level;
    bool             r0_based_qps_qpm;
    uint32_t         dpb_order_hint[REF_FRAMES]; // spec 6.8.2. ref_order_hint[]
    DGDetectorSeg   *dg_detector; // dg detector segments control struct
    SvtAv1RoiMapEvt *roi_map_evt;
    uint32_t         filt_to_unfilt_diff;
    // Absolute histogram deviation of the frame to the central TF frame (i.e. sum of absolute deviation of all bins)
    uint32_t tf_ahd_error_to_central;
    // Average absolute histogram deviation of all frames in the TF window to the current (central) frame
    uint32_t tf_avg_ahd_error;
    bool     tf_active_region_present;
    bool     seq_param_changed;
} PictureParentControlSet;

typedef struct TplDispResults {
    EbDctor                  dctor;
    EbObjectWrapper         *pcs_wrapper;
    uint32_t                 frame_index;
    EbFifo                  *sbo_feedback_fifo_ptr;
    uint32_t                 input_type;
    int16_t                  enc_dec_segment_row;
    uint16_t                 tile_group_index;
    PictureParentControlSet *pcs;
    int32_t                  qIndex;

} TplDispResults;

typedef struct PictureControlSetInitData {
    uint16_t      picture_width;
    uint16_t      picture_height;
    uint16_t      left_padding;
    uint16_t      right_padding;
    uint16_t      top_padding;
    uint16_t      bot_padding;
    EbBitDepth    bit_depth;
    EbColorFormat color_format;
    uint32_t      b64_size;
    uint8_t       cfg_palette;
    // since we still have lot of code assuming 64x64 SB, we add a new paramter supporting
    // both128x128 and 64x64, ultimately the fixed code supporting 64x64 should be upgraded to use
    // 128x128 and sb_size could be removed.
    uint32_t sb_size;
    // Bool                             is_16bit;
    uint32_t                 ten_bit_format;
    uint16_t                 enc_dec_segment_col;
    uint16_t                 enc_dec_segment_row;
    EncMode                  enc_mode;
    EbSvtAv1EncConfiguration static_config;
    uint8_t                  speed_control;
    int8_t                   hbd_md;
    uint8_t                  cdf_mode;
    uint8_t                  over_boundary_block_mode;
    uint8_t                  mfmv;
    // init value for child pcs
    uint8_t tile_row_count;
    uint8_t tile_column_count;

    // Init value for parent pcs
    // from command line
    uint8_t log2_tile_rows;
    uint8_t log2_tile_cols;
    // in mi unit
    uint8_t log2_sb_size;
    Bool    is_16bit_pipeline;

    uint16_t   non_m8_pad_w;
    uint16_t   non_m8_pad_h;
    uint8_t    enable_tpl_la;
    uint8_t    tpl_synth_size;
    uint8_t    in_loop_ois;
    uint8_t    pass;
    uint32_t   rate_control_mode;
    Av1Common *av1_cm;
    uint16_t   init_max_block_cnt;
    uint8_t    ref_count_used_list0;
    uint8_t    ref_count_used_list1;

    uint8_t enable_adaptive_quantization;
    uint8_t calc_hist;
    uint8_t tpl_lad_mg;
    uint8_t skip_frame_first_pass;
    uint8_t ipp_ds;
    uint8_t bypass_blk_step;
    uint8_t dist_ds;
    uint8_t ipp_was_ds;
    uint8_t final_pass_preset;
    uint8_t bypass_zz_check;
    uint8_t use8blk;
    uint8_t reduce_me_search;
    uint8_t input_resolution;
    uint8_t calculate_variance;
    Bool    is_scale;
    bool    rtc_tune;
    uint8_t variance_boost_strength;
    uint8_t new_variance_octile;
} PictureControlSetInitData;

typedef struct Av1Comp {
    Yv12BufferConfig trial_frame_rst;
} Av1Comp;

/**************************************
 * Extern Function Declarations
 **************************************/

uint32_t svt_aom_get_out_buffer_size(uint32_t picture_width, uint32_t picture_height);

EbErrorType svt_aom_picture_control_set_creator(EbPtr *object_dbl_ptr, EbPtr object_init_data_ptr);
EbErrorType svt_aom_recon_coef_creator(EbPtr *object_dbl_ptr, EbPtr object_init_data_ptr);
EbErrorType svt_aom_picture_parent_control_set_creator(EbPtr *object_dbl_ptr, EbPtr object_init_data_ptr);
EbErrorType svt_aom_me_creator(EbPtr *object_dbl_ptr, EbPtr object_init_data_ptr);
EbErrorType svt_aom_me_sb_results_ctor(MeSbResults *obj_ptr, PictureControlSetInitData *init_data_ptr);
EbErrorType ppcs_update_param(PictureParentControlSet *ppcs);
EbErrorType pcs_update_param(PictureControlSet *pcs);
EbErrorType me_update_param(MotionEstimationData *me_data, struct SequenceControlSet *scs);
EbErrorType recon_coef_update_param(EncDecSet *recon_coef, struct SequenceControlSet *scs);
extern Bool svt_aom_is_pic_skipped(PictureParentControlSet *pcs);
void svt_aom_get_gm_needed_resolutions(uint8_t ds_lvl, bool *gm_need_full, bool *gm_need_quart, bool *gm_need_sixteen);
#ifdef __cplusplus
}
#endif
#endif // EbPictureControlSet_h
