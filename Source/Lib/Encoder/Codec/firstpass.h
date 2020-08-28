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

#ifndef AOM_AV1_ENCODER_FIRSTPASS_H_
#define AOM_AV1_ENCODER_FIRSTPASS_H_

#include "EbDefinitions.h"
#include "EbRateControlProcess.h"
#include "EbPictureControlSet.h"

#if TWOPASS_RC
#ifdef __cplusplus
extern "C" {
#endif

#define TURN_OFF_EC_FIRST_PASS  1

#if FIRST_PASS_SETUP
#define FORCED_BLK_SIZE 16
#define FIRST_PASS_Q 10.0
#endif

#define DOUBLE_DIVIDE_CHECK(x) ((x) < 0 ? (x)-0.000001 : (x) + 0.000001)

#define MAX_LAG_BUFFERS 35
#define MIN_ZERO_MOTION 0.95
#define MAX_SR_CODED_ERROR 40
#define MAX_RAW_ERR_VAR 2000
#define MIN_MV_IN_OUT 0.4

#define VLOW_MOTION_THRESHOLD 950

// size of firstpass macroblocks in terms of MIs.
#define FP_MIB_SIZE 4
#define FP_MIB_SIZE_LOG2 2

// The maximum duration of a GF group that is static (e.g. a slide show).
#define MAX_STATIC_GF_GROUP_LENGTH 250
#define MAX_LAP_BUFFERS 35

/*!
 * \brief The stucture of acummulated frame stats in the first pass.
 */
typedef struct {
  /*!
   * Frame number in display order, if stats are for a single frame.
   * No real meaning for a collection of frames.
   */
  double frame;
  /*!
   * Weight assigned to this frame (or total weight for the collection of
   * frames) currently based on intra factor and brightness factor. This is used
   * to distribute bits betweeen easier and harder frames.
   */
  double weight;
  /*!
   * Intra prediction error.
   */
  double intra_error;
  /*!
   * Average wavelet energy computed using Discrete Wavelet Transform (DWT).
   */
  double frame_avg_wavelet_energy;
  /*!
   * Best of intra pred error and inter pred error using last frame as ref.
   */
  double coded_error;
  /*!
   * Best of intra pred error and inter pred error using golden frame as ref.
   */
  double sr_coded_error;
  /*!
   * Best of intra pred error and inter pred error using altref frame as ref.
   */
  double tr_coded_error;
  /*!
   * Percentage of blocks with inter pred error < intra pred error.
   */
  double pcnt_inter;
  /*!
   * Percentage of blocks using (inter prediction and) non-zero motion vectors.
   */
  double pcnt_motion;
  /*!
   * Percentage of blocks where golden frame was better than last or intra:
   * inter pred error using golden frame < inter pred error using last frame and
   * inter pred error using golden frame < intra pred error
   */
  double pcnt_second_ref;
  /*!
   * Percentage of blocks where altref frame was better than intra, last, golden
   */
  double pcnt_third_ref;
  /*!
   * Percentage of blocks where intra and inter prediction errors were very
   * close. Note that this is a 'weighted count', that is, the so blocks may be
   * weighted by how close the two errors were.
   */
  double pcnt_neutral;
  /*!
   * Percentage of blocks that have almost no intra error residual
   * (i.e. are in effect completely flat and untextured in the intra
   * domain). In natural videos this is uncommon, but it is much more
   * common in animations, graphics and screen content, so may be used
   * as a signal to detect these types of content.
   */
  double intra_skip_pct;
  /*!
   * Image mask rows top and bottom.
   */
  double inactive_zone_rows;
  /*!
   * Image mask columns at left and right edges.
   */
  double inactive_zone_cols;
  /*!
   * Average of row motion vectors.
   */
  double MVr;
  /*!
   * Mean of absolute value of row motion vectors.
   */
  double mvr_abs;
  /*!
   * Mean of column motion vectors.
   */
  double MVc;
  /*!
   * Mean of absolute value of column motion vectors.
   */
  double mvc_abs;
  /*!
   * Variance of row motion vectors.
   */
  double MVrv;
  /*!
   * Variance of column motion vectors.
   */
  double MVcv;
  /*!
   * Value in range [-1,1] indicating fraction of row and column motion vectors
   * that point inwards (negative MV value) or outwards (positive MV value).
   * For example, value of 1 indicates, all row/column MVs are inwards.
   */
  double mv_in_out_count;
  /*!
   * Count of unique non-zero motion vectors.
   */
  double new_mv_count;
  /*!
   * Duration of the frame / collection of frames.
   */
  double duration;
  /*!
   * 1.0 if stats are for a single frame, OR
   * Number of frames in this collection for which the stats are accumulated.
   */
  double count;
  /*!
   * standard deviation for (0, 0) motion prediction error
   */
  double raw_error_stdev;
} FIRSTPASS_STATS;

/*!\cond */

#define FC_ANIMATION_THRESH 0.15
enum {
  FC_NORMAL = 0,
  FC_GRAPHICS_ANIMATION = 1,
  FRAME_CONTENT_TYPES = 2
} UENUM1BYTE(FRAME_CONTENT_TYPE);

typedef struct {
  unsigned char index;
  /*frame_update_type*/int update_type[MAX_STATIC_GF_GROUP_LENGTH];
  unsigned char arf_src_offset[MAX_STATIC_GF_GROUP_LENGTH];
  // The number of frames displayed so far within the GOP at a given coding
  // frame.
  unsigned char cur_frame_idx[MAX_STATIC_GF_GROUP_LENGTH];
  unsigned char frame_disp_idx[MAX_STATIC_GF_GROUP_LENGTH];

  // TODO(jingning): Unify the data structure used here after the new control
  // mechanism is in place.
  int layer_depth[MAX_STATIC_GF_GROUP_LENGTH];
  int arf_boost[MAX_STATIC_GF_GROUP_LENGTH];
  int max_layer_depth;
  int max_layer_depth_allowed;
  // This is currently only populated for AOM_Q mode
  unsigned char q_val[MAX_STATIC_GF_GROUP_LENGTH];
  int bit_allocation[MAX_STATIC_GF_GROUP_LENGTH];
  int size;
} GF_GROUP;

typedef struct {
  FIRSTPASS_STATS *stats_in_start;
  FIRSTPASS_STATS *stats_in_end;
  FIRSTPASS_STATS *stats_in_buf_end;
  FIRSTPASS_STATS *total_stats;
  FIRSTPASS_STATS *total_left_stats;
} STATS_BUFFER_CTX;

/*!\endcond */

/*!
 * \brief Two pass status and control data.
 */
typedef struct {
  /*!\cond */
  unsigned int section_intra_rating;
  // Circular queue of first pass stats stored for most recent frames.
  // cpi->output_pkt_list[i].data.twopass_stats.buf points to actual data stored
  // here.
  FIRSTPASS_STATS *frame_stats_arr[MAX_LAP_BUFFERS + 1];
  int frame_stats_next_idx;  // Index to next unused element in frame_stats_arr.
  const FIRSTPASS_STATS *stats_in;
  STATS_BUFFER_CTX *stats_buf_ctx;
  int first_pass_done;
  int64_t bits_left;
  double modified_error_min;
  double modified_error_max;
  double modified_error_left;
  double mb_av_energy;
  double frame_avg_haar_energy;

  // An indication of the content type of the current frame
  FRAME_CONTENT_TYPE fr_content_type;

  // Projected total bits available for a key frame group of frames
  int64_t kf_group_bits;

  // Error score of frames still to be coded in kf group
  int64_t kf_group_error_left;

  // Over time correction for bits per macro block estimation
  double bpm_factor;

  // Record of target and actual bits spent in current ARF group
  int rolling_arf_group_target_bits;
  int rolling_arf_group_actual_bits;

  int sr_update_lag;

  int kf_zeromotion_pct;
  int last_kfgroup_zeromotion_pct;
  int extend_minq;
  int extend_maxq;
  int extend_minq_fast;
  /*!\endcond */
} TWO_PASS;

/*!\cond */

// This structure contains several key parameters to be accumulated for this
// frame.
typedef struct {
  // Intra prediction error.
  int64_t intra_error;
  // Average wavelet energy computed using Discrete Wavelet Transform (DWT).
  int64_t frame_avg_wavelet_energy;
  // Best of intra pred error and inter pred error using last frame as ref.
  int64_t coded_error;
  // Best of intra pred error and inter pred error using golden frame as ref.
  int64_t sr_coded_error;
  // Best of intra pred error and inter pred error using altref frame as ref.
  int64_t tr_coded_error;
  // Count of motion vector.
  int mv_count;
  // Count of blocks that pick inter prediction (inter pred error is smaller
  // than intra pred error).
  int inter_count;
  // Count of blocks that pick second ref (golden frame).
  int second_ref_count;
  // Count of blocks that pick third ref (altref frame).
  int third_ref_count;
  // Count of blocks where the inter and intra are very close and very low.
  double neutral_count;
  // Count of blocks where intra error is very small.
  int intra_skip_count;
  // Start row.
  int image_data_start_row;
  // Count of unique non-zero motion vectors.
  int new_mv_count;
  // Sum of inward motion vectors.
  int sum_in_vectors;
  // Sum of motion vector row.
  int sum_mvr;
  // Sum of motion vector column.
  int sum_mvc;
  // Sum of absolute value of motion vector row.
  int sum_mvr_abs;
  // Sum of absolute value of motion vector column.
  int sum_mvc_abs;
  // Sum of the square of motion vector row.
  int64_t sum_mvrs;
  // Sum of the square of motion vector column.
  int64_t sum_mvcs;
  // A factor calculated using intra pred error.
  double intra_factor;
  // A factor that measures brightness.
  double brightness_factor;
} FRAME_STATS;

// This structure contains first pass data.
typedef struct {
  // Buffer holding frame stats for all MACROBLOCKs.
  // mb_stats[i] stores the FRAME_STATS of the ith
  // MB in raster scan order.
  FRAME_STATS *mb_stats;
  // Buffer to store the prediction error of the (0,0) motion
  // vector using the last source frame as the reference.
  // raw_motion_err_list[i] stores the raw_motion_err of
  // the ith MB in raster scan order.
  int *raw_motion_err_list;
} FirstPassData;

struct EncodeFrameParams;
struct AV1EncoderConfig;
struct TileDataEnc;

void svt_av1_twopass_zero_stats(FIRSTPASS_STATS *section);
void svt_av1_accumulate_stats(FIRSTPASS_STATS *section,
                          const FIRSTPASS_STATS *frame);
#if FIRST_PASS_SETUP
struct PictureParentControlSet;
extern void svt_av1_end_first_pass(struct PictureParentControlSet *pcs_ptr);
extern void first_pass_frame_end(struct PictureParentControlSet *pcs_ptr, const int64_t ts_duration);
extern void setup_firstpass_data(struct PictureParentControlSet *pcs_ptr);
extern void average_non_16x16_stats(FRAME_STATS *mb_stats, int blk_num);
void accumulate_mv_stats(const MV best_mv, const FULLPEL_MV mv,
    const int mb_row, const int mb_col,
    const int mb_rows, const int mb_cols,
    MV *last_mv, FRAME_STATS *stats);
#endif
/*!\endcond */

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // TWOPASS_RC
#endif  // AOM_AV1_ENCODER_FIRSTPASS_H_
