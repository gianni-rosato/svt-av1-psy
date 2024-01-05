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

/*!\file
 * \brief Declares top-level encoder structures and functions.
 */
#ifndef AOM_AV1_ENCODER_ENCODER_H_
#define AOM_AV1_ENCODER_ENCODER_H_

#include <stdbool.h>
#include <stdio.h>

#include "EbDefinitions.h"
#include "Av1Common.h"
#include "EbRateControlProcess.h"
#include "level.h"

#ifdef __cplusplus
extern "C" {
#endif

// TODO(yunqing, any): Added suppression tag to quiet Doxygen warnings. Need to
// adjust it while we work on documentation.
/*!\cond */
// Number of frames required to test for scene cut detection
#define SCENE_CUT_KEY_TEST_INTERVAL 16

#define FRAME_TYPE int

//**********************************************************************************************************************//
// aom_codec.h
/*!\brief Rate control mode */
enum aom_rc_mode {
    AOM_VBR, /**< Variable Bit Rate (VBR) mode */
    AOM_CBR, /**< Constant Bit Rate (CBR) mode */
    AOM_CQ, /**< Constrained Quality (CQ)  mode */
    AOM_Q, /**< Constant Quality (Q) mode */
};
//**********************************************************************************************************************//
//struct AV1LevelParams;

/*!\endcond */
/*!
 * \brief Encoder rate control configuration parameters
 */
typedef struct {
    /*!\cond */
    // BUFFERING PARAMETERS
    // Indicates the amount of data that will be buffered by the decoding
    // application prior to beginning playback, and is expressed in units of
    // time(milliseconds).
    int64_t starting_buffer_level_ms;
    // Indicates the amount of data that the encoder should try to maintain in the
    // decoder's buffer, and is expressed in units of time(milliseconds).
    int64_t optimal_buffer_level_ms;
    // Indicates the maximum amount of data that may be buffered by the decoding
    // application, and is expressed in units of time(milliseconds).
    int64_t maximum_buffer_size_ms;

    // Indicates the maximum allowed bitrate for any intra frame as % of bitrate
    // target.
    unsigned int max_intra_bitrate_pct;
    // Indicates the maximum allowed bitrate for any inter frame as % of bitrate
    // target.
    unsigned int max_inter_bitrate_pct;
    // Indicates the percentage of rate boost for golden frame in CBR mode.
    unsigned int gf_cbr_boost_pct;
    // min_cr / 100 indicates the target minimum compression ratio for each frame.
    unsigned int min_cr;
    // under_shoot_pct indicates the tolerance of the VBR algorithm to undershoot
    // and is used as a trigger threshold for more agressive adaptation of Q. It's
    // value can range from 0-100.
    int under_shoot_pct;
    // over_shoot_pct indicates the tolerance of the VBR algorithm to overshoot
    // and is used as a trigger threshold for more agressive adaptation of Q. It's
    // value can range from 0-1000.
    int over_shoot_pct;
    // Indicates the maximum qindex that can be used by the quantizer i.e. the
    // worst quality qindex.
    int worst_allowed_q;
    // Indicates the minimum qindex that can be used by the quantizer i.e. the
    // best quality qindex.
    int best_allowed_q;
    // Indicates if the encoding mode is vbr, cbr, constrained quality or constant
    // quality.
    enum aom_rc_mode mode;
    /*!\endcond */
} RateControlCfg;

typedef int aom_bit_depth_t;
typedef struct {
    int             frame_width;
    int             frame_height;
    int             mb_rows;
    int             mb_cols;
    int             num_mbs;
    aom_bit_depth_t bit_depth;
    int             subsampling_x;
    int             subsampling_y;
} FrameInfo;

typedef struct {
    // stats_in buffer contains all of the stats packets produced in the first
    // pass, concatenated.
    //aom_fixed_buf_t stats_in;

    // Indicates the minimum bitrate to be used for a single GOP as a percentage
    // of the target bitrate.
    int vbrmin_section;
    // Indicates the maximum bitrate to be used for a single GOP as a percentage
    // of the target bitrate.
    int vbrmax_section;
} TwoPassCfg;
/*!
 * \brief Main encoder configuration data structure.
 */
typedef struct AV1EncoderConfig {
    /*!\cond */
    BITSTREAM_PROFILE profile;
    aom_bit_depth_t   bit_depth; // Codec bit-depth.
    int64_t           target_bandwidth; // bandwidth to be used in bits per second

    // Configuration related to the input video.
    //InputCfg input_cfg;

    // Configuration related to frame-dimensions.
    //FrameDimensionCfg frm_dim_cfg;

    int sharpness; // sharpening output: recommendation 0:
    int speed;

    int /*MODE*/ mode;
    int          pass;
    // ----------------------------------------------------------------
    // DATARATE CONTROL OPTIONS

    /*!\endcond */
    /*!
   * Rate control configuration
   */
    RateControlCfg rc_cfg;
    /*!\cond */

    // Frame drop threshold.
    int drop_frames_water_mark;

    // controlling quality
    int          deltalf_mode;
    int          enable_cdef;
    int          enable_restoration;
    int          force_video_mode;
    int          disable_trellis_quant;
    unsigned int vbr_corpus_complexity_lap; // 0 indicates corpus complexity vbr
        // mode is disabled
    // two pass datarate control
    TwoPassCfg two_pass_cfg;

    // END DATARATE CONTROL OPTIONS
    // ----------------------------------------------------------------

    /* Bitfield defining the error resiliency features to enable.
   * Can provide decodable frames after losses in previous
   * frames and decodable partitions after losses in the same frame.
   */
    unsigned int error_resilient_mode;

    /* Bitfield defining the parallel decoding mode where the
   * decoding in successive frames may be conducted in parallel
   * just by decoding the frame headers.
   */
    unsigned int frame_parallel_decoding_mode;

    int arnr_max_frames;
    int arnr_strength;

    // Tile related configuration parameters.
    //TileConfig tile_cfg;

    int row_mt;

    int enable_tpl_model;

    int max_threads;

    //aom_tune_metric tuning;
    const char *vmaf_model_path;
    //aom_tune_content content;
    int use_highbitdepth;
    //aom_chroma_sample_position_t chroma_sample_position;
    int         film_grain_test_vector;

    // Configuration related to color.
    //ColorCfg color_cfg;

    // Configuration related to decoder model.
    //DecoderModelCfg dec_model_cfg;

    // Configuration related to reference frames.
    //RefFrameCfg ref_frm_cfg;

    // Configuration related to unit tests.
    //UnitTestCfg unit_test_cfg;

    uint8_t cdf_update_mode;
    //aom_superblock_size_t superblock_size;
    uint8_t      monochrome;
    unsigned int full_still_picture_hdr;
    int          enable_dual_filter;
    int          enable_order_hint;
    int          enable_ref_frame_mvs;
    unsigned int allow_ref_frame_mvs;
    int          enable_interintra_comp;
    int          enable_global_motion;
    int          enable_overlay;
    int          enable_palette;
    unsigned int save_as_annexb;

    // Flags related to motion mode.
    //MotionModeCfg motion_mode_cfg;

    // Flags related to intra mode search.
    //IntraModeCfg intra_mode_cfg;

    // Flags related to transform size/type.
    //TxfmSizeTypeCfg txfm_cfg;

    // Flags related to compound type.
    //CompoundTypeCfg comp_type_cfg;

    // Partition related information.
    //PartitionCfg part_cfg;

#if CONFIG_DENOISE
    float noise_level;
    int   noise_block_size;
#endif

    // Configuration related to frequency of cost update.
    //CostUpdateFreq cost_upd_freq;

    int       border_in_pixels;
    AV1_LEVEL target_seq_level_idx[MAX_NUM_OPERATING_POINTS];
    // Bit mask to specify which tier each of the 32 possible operating points
    // conforms to.
    unsigned int tier_mask;

    /*!\endcond */
} AV1EncoderConfig;

typedef struct SwitchFrameCfg {
    // Indicates the number of frames after which a frame may be coded as an S-Frame.
    int32_t sframe_dist;
    // 1: the considered frame will be made into an S-Frame only if it is an altref frame.
    // 2: the next altref frame will be made into an S-Frame.
    EbSFrameMode sframe_mode;
} SwitchFrameCfg;

#define MAX_GFUBOOST_FACTOR 10.0
#define MIN_GFUBOOST_FACTOR 4.0

// Function return size of frame stats buffer
static INLINE int get_stats_buf_size(int num_lap_buffer, int num_lag_buffer) {
    /* if lookahead is enabled return num_lap_buffers else num_lag_buffers */
    return (num_lap_buffer > 0 ? num_lap_buffer + 1 : num_lag_buffer);
}

#ifdef __cplusplus
} // extern "C"
#endif

#endif // AOM_AV1_ENCODER_ENCODER_H_
