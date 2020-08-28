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

#if TWOPASS_RC
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
  AOM_CQ,  /**< Constrained Quality (CQ)  mode */
  AOM_Q,   /**< Constant Quality (Q) mode */
};
//**********************************************************************************************************************//

enum {
  DISABLE_SCENECUT,        // For LAP, lag_in_frames < 19
  ENABLE_SCENECUT_MODE_1,  // For LAP, lag_in_frames >=19 and < 33
  ENABLE_SCENECUT_MODE_2   // For twopass and LAP - lag_in_frames >=33
} UENUM1BYTE(SCENECUT_MODE);

//struct AV1LevelParams;

typedef struct {
  // Indicates the minimum distance to a key frame.
  int key_freq_min;
  // Indicates the maximum distance to a key frame.
  int key_freq_max;
  // Indicates if temporal filtering should be applied on keyframe.
  int enable_keyframe_filtering;
  // Indicates the number of frames after which a frame may be coded as an
  // S-Frame.
  int sframe_dist;
  // Indicates how an S-Frame should be inserted.
  // 1: the considered frame will be made into an S-Frame only if it is an
  // altref frame. 2: the next altref frame will be made into an S-Frame.
  int sframe_mode;
  // Indicates if encoder should autodetect cut scenes and set the keyframes.
  bool auto_key;
  // Indicates if forward keyframe reference should be enabled.
  bool fwd_kf_enabled;
  // Indicates if S-Frames should be enabled for the sequence.
  bool enable_sframe;
  // Indicates if intra block copy prediction mode should be enabled or not.
  bool enable_intrabc;
} KeyFrameCfg;

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
  // Indicates the Constant/Constrained Quality level.
  int cq_level;
  // Indicates if the encoding mode is vbr, cbr, constrained quality or constant
  // quality.
  enum aom_rc_mode mode;
  /*!\endcond */
} RateControlCfg;

/*!\cond */
typedef struct {
  // Indicates the number of frames lag before encoding is started.
  int lag_in_frames;
  // Indicates the minimum gf/arf interval to be used.
  int min_gf_interval;
  // Indicates the maximum gf/arf interval to be used.
  int max_gf_interval;
  // Indicates the minimum height for GF group pyramid structure to be used.
  int gf_min_pyr_height;
  // Indicates the maximum height for GF group pyramid structure to be used.
  int gf_max_pyr_height;
  // Indicates if automatic set and use of altref frames should be enabled.
  bool enable_auto_arf;
  // Indicates if automatic set and use of (b)ackward (r)ef (f)rames should be
  // enabled.
  bool enable_auto_brf;
} GFConfig;

typedef int aom_bit_depth_t;
typedef struct {
  int frame_width;
  int frame_height;
  int mi_rows;
  int mi_cols;
  int mb_rows;
  int mb_cols;
  int num_mbs;
  aom_bit_depth_t bit_depth;
  int subsampling_x;
  int subsampling_y;
} FrameInfo;

typedef struct EncodeFrameInput {
  //YV12_BUFFER_CONFIG *source;
  //YV12_BUFFER_CONFIG *last_source;
  int64_t ts_duration;
} EncodeFrameInput;

// EncodeFrameParams contains per-frame encoding parameters decided upon by
// av1_encode_strategy() and passed down to av1_encode()
struct EncodeFrameParams {
  int error_resilient_mode;
  FRAME_TYPE frame_type;
  int primary_ref_frame;
  int order_offset;
  int show_frame;
  int refresh_frame_flags;

  int show_existing_frame;
  int existing_fb_idx_to_show;

  // Bitmask of which reference buffers may be referenced by this frame
  int ref_frame_flags;

  // Reference buffer assignment for this frame.
  int remapped_ref_idx[REF_FRAMES];

  // Flags which determine which reference buffers are refreshed by this frame.
  //RefreshFrameFlagsInfo refresh_frame;

  // Speed level to use for this frame: Bigger number means faster.
  int speed;
};

typedef struct {
  // stats_in buffer contains all of the stats packets produced in the first
  // pass, concatenated.
  //aom_fixed_buf_t stats_in;

  // TWO PASS DATARATE CONTROL OPTIONS.
  // Indicates the bias (expressed on a scale of 0 to 100) for determining
  // target size for the current frame. The value 0 indicates the optimal CBR
  // mode value should be used, and 100 indicates the optimal VBR mode value
  // should be used.
  int vbrbias;
  // Indicates the minimum bitrate to be used for a single GOP as a percentage
  // of the target bitrate.
  int vbrmin_section;
  // Indicates the maximum bitrate to be used for a single GOP as a percentage
  // of the target bitrate.
  int vbrmax_section;
} TwoPassCfg;

typedef struct EncodeFrameParams EncodeFrameParams;

/*!
 * \brief Main encoder configuration data structure.
 */
typedef struct AV1EncoderConfig {
  /*!\cond */
  BITSTREAM_PROFILE profile;
  aom_bit_depth_t bit_depth;  // Codec bit-depth.
  int64_t target_bandwidth;   // bandwidth to be used in bits per second

  // Configuration related to the input video.
  //InputCfg input_cfg;

  // Configuration related to frame-dimensions.
  //FrameDimensionCfg frm_dim_cfg;

  int sharpness;  // sharpening output: recommendation 0:
  int speed;

  int/*MODE*/ mode;
  int pass;

  // Configuration related to key-frame.
  KeyFrameCfg kf_cfg;

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
  int deltalf_mode;
  int enable_cdef;
  int enable_restoration;
  int force_video_mode;
  int disable_trellis_quant;
  unsigned int vbr_corpus_complexity_lap;  // 0 indicates corpus complexity vbr
                                           // mode is disabled

  // Configuration related to Quantization.
  //QuantizationCfg q_cfg;

  // Internal frame size scaling.
  //ResizeCfg resize_cfg;

  // Frame Super-Resolution size scaling.
  //SuperResCfg superres_cfg;

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

  // Configuration related to Group of frames.
  GFConfig gf_cfg;

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
  int film_grain_test_vector;
  const char *film_grain_table_filename;

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
  uint8_t monochrome;
  unsigned int full_still_picture_hdr;
  int enable_dual_filter;
  int enable_order_hint;
  int enable_ref_frame_mvs;
  unsigned int allow_ref_frame_mvs;
  int enable_interintra_comp;
  int enable_global_motion;
  int enable_overlay;
  int enable_palette;
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
  int noise_block_size;
#endif

  // Configuration related to frequency of cost update.
  //CostUpdateFreq cost_upd_freq;

  int border_in_pixels;
  AV1_LEVEL target_seq_level_idx[MAX_NUM_OPERATING_POINTS];
  // Bit mask to specify which tier each of the 32 possible operating points
  // conforms to.
  unsigned int tier_mask;

  /*!\endcond */
} AV1EncoderConfig;


/*!
 * \brief Top level encoder structure.
 */
typedef struct AV1_COMP {
  /*!
   * Quantization and dequantization parameters for internal quantizer setup
   * in the encoder.
   */
  //EncQuantDequantParams enc_quant_dequant_params;

  /*!
   * Structure holding thread specific variables.
   */
  //ThreadData td;

  /*!
   * Statistics collected at frame level.
   */
  //FRAME_COUNTS counts;

  /*!
   * Holds buffer storing mode information at 4x4/8x8 level.
   */
  //MBMIExtFrameBufferInfo mbmi_ext_info;

  /*!
   * Buffer holding the transform block related information.
   * coeff_buffer_base[i] stores the transform block related information of the
   * ith superblock in raster scan order.
   */
  //CB_COEFF_BUFFER *coeff_buffer_base;

  /*!
   * Structure holding variables common to encoder and decoder.
   */
  AV1_COMMON common;

  /*!
   * Encoder configuration related parameters.
   */
  AV1EncoderConfig oxcf;

  /*!
   * Look-ahead context.
   */
  struct lookahead_ctx *lookahead;

  /*!
   * When set, this flag indicates that the current frame is a forward keyframe.
   */
  int no_show_kf;

  /*!
   * Stores the trellis optimization type at segment level.
   * optimize_seg_arr[i] stores the trellis opt type for ith segment.
   */
  //TRELLIS_OPT_TYPE optimize_seg_arr[MAX_SEGMENTS];

  /*!
   * Pointer to the frame buffer holding the source frame to be used during the
   * current stage of encoding. It can be the raw input, temporally filtered
   * input or scaled input.
   */
  //YV12_BUFFER_CONFIG *source;

  /*!
   * Pointer to the frame buffer holding the last raw source frame.
   * NULL for first frame and alt_ref frames.
   */
  //YV12_BUFFER_CONFIG *last_source;

  /*!
   * Pointer to the frame buffer holding the unscaled source frame.
   * It can be either the raw input or temporally filtered input.
   */
  //YV12_BUFFER_CONFIG *unscaled_source;

  /*!
   * Frame buffer holding the resized source frame (cropping / superres).
   */
  //YV12_BUFFER_CONFIG scaled_source;

  /*!
   * Pointer to the frame buffer holding the unscaled last source frame.
   */
  //YV12_BUFFER_CONFIG *unscaled_last_source;

  /*!
   * Frame buffer holding the resized last source frame.
   */
  //YV12_BUFFER_CONFIG scaled_last_source;

  /*!
   * Pointer to the original source frame. This is used to determine if the
   * content is screen.
   */
  //YV12_BUFFER_CONFIG *unfiltered_source;

  /*!
   * Parameters related to tpl.
   */
  //TplParams tpl_data;

  /*!
   * For a still frame, this flag is set to 1 to skip partition search.
   */
  int partition_search_skippable_frame;

  /*!
   * Variables related to forcing integer mv decisions for the current frame.
   */
  //ForceIntegerMVInfo force_intpel_info;

  /*!
   * Pointer to the buffer holding the scaled reference frames.
   * scaled_ref_buf[i] holds the scaled reference frame of type i.
   */
  //RefCntBuffer *scaled_ref_buf[INTER_REFS_PER_FRAME];

  /*!
   * Pointer to the buffer holding the last show frame.
   */
  //RefCntBuffer *last_show_frame_buf;

  /*!
   * Refresh frame flags for golden, bwd-ref and alt-ref frames.
   */
  //RefreshFrameFlagsInfo refresh_frame;

  /*!
   * For each type of reference frame, this contains the index of a reference
   * frame buffer for a reference frame of the same type.  We use this to
   * choose our primary reference frame (which is the most recent reference
   * frame of the same type as the current frame).
   */
  int fb_of_context_type[REF_FRAMES];

  /*!
   * Flags signalled by the external interface at frame level.
   */
  //ExternalFlags ext_flags;

  /*!
   * Temporary frame buffer used to store the non-loop filtered reconstructed
   * frame during the search of loop filter level.
   */
  //YV12_BUFFER_CONFIG last_frame_uf;

  /*!
   * Temporary frame buffer used to store the loop restored frame during loop
   * restoration search.
   */
  //YV12_BUFFER_CONFIG trial_frame_rst;

  /*!
   * Ambient reconstruction err target for force key frames.
   */
  int64_t ambient_err;

  /*!
   * Parameters related to rate distortion optimization.
   */
  //RD_OPT rd;

  /*!
   * Temporary coding context used to save and restore when encoding with and
   * without super-resolution.
   */
  //CODING_CONTEXT coding_context;

  /*!
   * Parameters related to global motion search.
   */
  //GlobalMotionInfo gm_info;

  /*!
   * Parameters related to winner mode processing.
   */
  //WinnerModeParams winner_mode_params;

  /*!
   * Frame time stamps.
   */
  //TimeStamps time_stamps;

  /*!
   * Rate control related parameters.
   */
  RATE_CONTROL rc;

  /*!
   * Frame rate of the video.
   */
  double framerate;

  /*!
   * Pointer to internal utility functions that manipulate aom_codec_* data
   * structures.
   */
  struct aom_codec_pkt_list *output_pkt_list;

  /*!
   * Bitmask indicating which reference buffers may be referenced by this frame.
   */
  int ref_frame_flags;

  /*!
   * speed is passed as a per-frame parameter into the encoder.
   */
  int speed;

  /*!
   * sf contains fine-grained config set internally based on speed.
   */
  //SPEED_FEATURES sf;

  /*!
   * Parameters for motion vector search process.
   */
  //MotionVectorSearchParams mv_search_params;

  /*!
   * When set, indicates that all reference frames are forward references,
   * i.e., all the reference frames are output before the current frame.
   */
  int all_one_sided_refs;

  /*!
   * Segmentation related information for current frame.
   */
  //EncSegmentationInfo enc_seg;

  /*!
   * Parameters related to cyclic refresh aq-mode.
   */
  //CYCLIC_REFRESH *cyclic_refresh;
  /*!
   * Parameters related to active map. Active maps indicate
   * if there is any activity on a 4x4 block basis.
   */
  //ActiveMap active_map;

  /*!
   * Function pointers to variants of sse/sad/variance computation functions.
   * fn_ptr[i] indicates the list of function pointers corresponding to block
   * size i.
   */
  //aom_variance_fn_ptr_t fn_ptr[BLOCK_SIZES_ALL];

  /*!
   * Information related to two pass encoding.
   */
  TWO_PASS twopass;

  /*!
   * Information related to a gf group.
   */
  GF_GROUP gf_group;

  /*!
   * To control the reference frame buffer and selection.
   */
  //RefBufferStack ref_buffer_stack;

  /*!
   * Frame buffer holding the temporally filtered source frame. It can be KEY
   * frame or ARF frame.
   */
  //YV12_BUFFER_CONFIG alt_ref_buffer;

  /*!
   * Tell if OVERLAY frame shows existing alt_ref frame.
   */
  int show_existing_alt_ref;

#if CONFIG_INTERNAL_STATS
  /*!\cond */
  uint64_t time_receive_data;
  uint64_t time_compress_data;

  unsigned int mode_chosen_counts[MAX_MODES];

  int count;
  uint64_t total_sq_error;
  uint64_t total_samples;
  ImageStat psnr;

  double total_blockiness;
  double worst_blockiness;

  int bytes;
  double summed_quality;
  double summed_weights;
  unsigned int tot_recode_hits;
  double worst_ssim;

  ImageStat fastssim;
  ImageStat psnrhvs;

  int b_calculate_blockiness;
  int b_calculate_consistency;

  double total_inconsistency;
  double worst_consistency;
  Ssimv *ssim_vars;
  Metrics metrics;
  /*!\endcond */
#endif

  /*!
   * Calculates PSNR on each frame when set to 1.
   */
  int b_calculate_psnr;

#if CONFIG_SPEED_STATS
  /*!
   * For debugging: number of transform searches we have performed.
   */
  unsigned int tx_search_count;
#endif  // CONFIG_SPEED_STATS

  /*!
   * When set, indicates that the frame is droppable, i.e., this frame
   * does not update any reference buffers.
   */
  int droppable;

  /*!
   * Stores the frame parameters during encoder initialization.
   */
  FrameInfo frame_info;

  /*!
   * Structure to store the dimensions of current frame.
   */
  //InitialDimensions initial_dimensions;

  /*!
   * Number of MBs in the full-size frame; to be used to
   * normalize the firstpass stats. This will differ from the
   * number of MBs in the current frame when the frame is
   * scaled.
   */
  int initial_mbs;

  /*!
   * Resize related parameters.
   */
  //ResizePendingParams resize_pending_params;

  /*!
   * Pointer to struct holding adaptive data/contexts/models for the tile during
   * encoding.
   */
  //TileDataEnc *tile_data;
  /*!
   * Number of tiles for which memory has been allocated for tile_data.
   */
  int allocated_tiles;

  /*!
   * Structure to store the palette token related information.
   */
  //TokenInfo token_info;

  /*!
   * Sequence parameters have been transmitted already and locked
   * or not. Once locked av1_change_config cannot change the seq
   * parameters.
   */
  int seq_params_locked;

  /*!
   * VARIANCE_AQ segment map refresh.
   */
  int vaq_refresh;

  /*!
   * Thresholds for variance based partitioning.
   */
  //VarBasedPartitionInfo vbp_info;

  /*!
   * Probabilities for pruning of various AV1 tools.
   */
  //FrameProbInfo frame_probs;

  /*!
   * Multi-threading parameters.
   */
  //MultiThreadInfo mt_info;

  /*!
   * Specifies the frame to be output. It is valid only if show_existing_frame
   * is 1. When show_existing_frame is 0, existing_fb_idx_to_show is set to
   * INVALID_IDX.
   */
  int existing_fb_idx_to_show;

  /*!
   * When set, indicates that internal ARFs are enabled.
   */
  int internal_altref_allowed;

  /*!
   * A flag to indicate if intrabc is ever used in current frame.
   */
  int intrabc_used;

  /*!
   * Tables to calculate IntraBC MV cost.
   */
  //IntraBCMVCosts dv_costs;

  /*!
   * Mark which ref frames can be skipped for encoding current frame during RDO.
   */
  int prune_ref_frame_mask;

  /*!
   * Loop Restoration context.
   */
  //AV1LrStruct lr_ctxt;

  /*!
   * Pointer to list of tables with film grain parameters.
   */
  //aom_film_grain_table_t *film_grain_table;

#if CONFIG_DENOISE
  /*!
   * Pointer to structure holding the denoised image buffers and the helper
   * noise models.
   */
  struct aom_denoise_and_model_t *denoise_and_model;
#endif

  /*!
   * Flags related to interpolation filter search.
   */
  //InterpSearchFlags interp_search_flags;

  /*!
   * Set for screen contents or when screen content tools are enabled.
   */
  int is_screen_content_type;

#if CONFIG_COLLECT_PARTITION_STATS == 2
  PartitionStats partition_stats;
#endif

#if CONFIG_COLLECT_COMPONENT_TIMING
  /*!
   * component_time[] are initialized to zero while encoder starts.
   */
  uint64_t component_time[kTimingComponents];
  struct aom_usec_timer component_timer[kTimingComponents];
  /*!
   * frame_component_time[] are initialized to zero at beginning of each frame.
   */
  uint64_t frame_component_time[kTimingComponents];
#endif

  /*!
   * Parameters for AV1 bitstream levels.
   */
  AV1LevelParams level_params;

  /*!
   * Whether any no-zero delta_q was actually used.
   */
  int deltaq_used;

  /*!
   * Refrence frame distance related variables.
   */
  //RefFrameDistanceInfo ref_frame_dist_info;

  /*!
   * Scaling factors used in the RD multiplier modulation.
   * TODO(sdeng): consider merge the following arrays.
   * tpl_rdmult_scaling_factors is a temporary buffer used to store the
   * intermediate scaling factors which are used in the calculation of
   * tpl_sb_rdmult_scaling_factors. tpl_rdmult_scaling_factors[i] stores the
   * intermediate scaling factor of the ith 16 x 16 block in raster scan order.
   */
  double *tpl_rdmult_scaling_factors;
  /*!
   * tpl_sb_rdmult_scaling_factors[i] stores the RD multiplier scaling factor of
   * the ith 16 x 16 block in raster scan order.
   */
  double *tpl_sb_rdmult_scaling_factors;
  /*!
   * ssim_rdmult_scaling_factors[i] stores the RD multiplier scaling factor of
   * the ith 16 x 16 block in raster scan order. This scaling factor is used for
   * RD multiplier modulation when SSIM tuning is enabled.
   */
  double *ssim_rdmult_scaling_factors;

#if CONFIG_TUNE_VMAF
  /*!
   * Parameters for VMAF tuning.
   */
  TuneVMAFInfo vmaf_info;
#endif

  /*!
   * Indicates whether to use SVC.
   */
  int use_svc;
  /*!
   * Parameters for scalable video coding.
   */
  //SVC svc;

  /*!
   * Flag indicating whether look ahead processing (LAP) is enabled.
   */
  int lap_enabled;
  /*!
   * Indicates whether current processing stage is encode stage or LAP stage.
   */
  //COMPRESSOR_STAGE compressor_stage;

  /*!
   * Some motion vector stats from the last encoded frame to help us decide what
   * precision to use to encode the current frame.
   */
  //MV_STATS mv_stats;

  /*!
   * Frame type of the last frame. May be used in some heuristics for speeding
   * up the encoding.
   */
  FRAME_TYPE last_frame_type;

  /*!
   * Number of tile-groups.
   */
  int num_tg;

  /*!
   * Super-resolution mode currently being used by the encoder.
   * This may / may not be same as user-supplied mode in oxcf->superres_mode
   * (when we are recoding to try multiple options for example).
   */
  //aom_superres_mode superres_mode;

  /*!
   * First pass related data.
   */
  FirstPassData firstpass_data;

  /*!
   * Temporal Noise Estimate
   */
  //NOISE_ESTIMATE noise_estimate;

  /*!
   * Count on how many consecutive times a block uses small/zeromv for encoding
   * in a scale of 8x8 block.
   */
  uint8_t *consec_zero_mv;
} AV1_COMP;

#define MAX_GFUBOOST_FACTOR 10.0
#define MIN_GFUBOOST_FACTOR 4.0

// Function return size of frame stats buffer
static INLINE int get_stats_buf_size(int num_lap_buffer, int num_lag_buffer) {
    /* if lookahead is enabled return num_lap_buffers else num_lag_buffers */
    return (num_lap_buffer > 0 ? num_lap_buffer + 1 : num_lag_buffer);
}

/*!\cond */
static INLINE int is_lossless_requested(const RateControlCfg *const rc_cfg) {
  return rc_cfg->best_allowed_q == 0 && rc_cfg->worst_allowed_q == 0;
}

#define ALT_MIN_LAG 3
static INLINE int is_altref_enabled(int lag_in_frames, bool enable_auto_arf) {
  return lag_in_frames >= ALT_MIN_LAG && enable_auto_arf;
}

// Check if statistics consumption stage
static INLINE int is_stat_consumption_stage_twopass(const AV1_COMP *const cpi) {
  return (cpi->oxcf.pass == 2);
}

// Check if statistics consumption stage
static INLINE int is_stat_consumption_stage(const AV1_COMP *const cpi) {
  return (is_stat_consumption_stage_twopass(cpi) ||
          (cpi->oxcf.pass == 0 /*&& (cpi->compressor_stage == ENCODE_STAGE)*/ &&
           cpi->lap_enabled));
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // TWOPASS_RC
#endif  // AOM_AV1_ENCODER_ENCODER_H_
