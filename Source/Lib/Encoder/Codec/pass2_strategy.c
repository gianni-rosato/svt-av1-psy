/*
 * Copyright (c) 2019, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <stdint.h>

#include "EbDefinitions.h"
#include "level.h"
#include "pass2_strategy.h"
#include "EbRateControlProcess.h"
#include "firstpass.h"
#include "EbSequenceControlSet.h"
#include "EbEntropyCoding.h"

//#define INT_MAX 0x7fffffff

#define DEFAULT_KF_BOOST 2300
#define DEFAULT_GF_BOOST 2000
#define GROUP_ADAPTIVE_MAXQ 1
static void init_gf_stats(GF_GROUP_STATS *gf_stats);

// Calculate an active area of the image that discounts formatting
// bars and partially discounts other 0 energy areas.
#define MIN_ACTIVE_AREA 0.5
#define MAX_ACTIVE_AREA 1.0
static double calculate_active_area(const FrameInfo *frame_info,
                                    const FIRSTPASS_STATS *this_frame) {
  const double active_pct =
      1.0 -
      ((this_frame->intra_skip_pct / 2) +
       ((this_frame->inactive_zone_rows * 2) / (double)frame_info->mb_rows));
  return fclamp(active_pct, MIN_ACTIVE_AREA, MAX_ACTIVE_AREA);
}

// Calculate a modified Error used in distributing bits between easier and
// harder frames.
#define ACT_AREA_CORRECTION 0.5
static double calculate_modified_err(const FrameInfo *frame_info,
                                     const TWO_PASS *twopass,
                                     const TwoPassCfg *two_pass_cfg,
                                     const FIRSTPASS_STATS *this_frame) {
  const FIRSTPASS_STATS *const stats = twopass->stats_buf_ctx->total_stats;
  if (stats == NULL) {
    return 0;
  }
  const double av_weight = stats->weight / stats->count;
  const double av_err = (stats->coded_error * av_weight) / stats->count;
  double modified_error =
      av_err * pow(this_frame->coded_error * this_frame->weight /
                       DOUBLE_DIVIDE_CHECK(av_err),
                   two_pass_cfg->vbrbias / 100.0);

  // Correction for active area. Frames with a reduced active area
  // (eg due to formatting bars) have a higher error per mb for the
  // remaining active MBs. The correction here assumes that coding
  // 0.5N blocks of complexity 2X is a little easier than coding N
  // blocks of complexity X.
  modified_error *=
      pow(calculate_active_area(frame_info, this_frame), ACT_AREA_CORRECTION);

  return fclamp(modified_error, twopass->modified_error_min,
                twopass->modified_error_max);
}

// Resets the first pass file to the given position using a relative seek from
// the current position.
static void reset_fpf_position(TWO_PASS *p, const FIRSTPASS_STATS *position) {
  p->stats_in = position;
}

static int input_stats(TWO_PASS *p, FIRSTPASS_STATS *fps) {
  if (p->stats_in >= p->stats_buf_ctx->stats_in_end) return EOF;

  *fps = *p->stats_in;
  ++p->stats_in;
  return 1;
}

// Read frame stats at an offset from the current position.
static const FIRSTPASS_STATS *read_frame_stats(const TWO_PASS *p, int offset) {
  if ((offset >= 0 && p->stats_in + offset >= p->stats_buf_ctx->stats_in_end) ||
      (offset < 0 && p->stats_in + offset < p->stats_buf_ctx->stats_in_start)) {
    return NULL;
  }

  return &p->stats_in[offset];
}

static void subtract_stats(FIRSTPASS_STATS *section,
                           const FIRSTPASS_STATS *frame) {
  section->frame -= frame->frame;
  section->weight -= frame->weight;
  section->intra_error -= frame->intra_error;
  section->frame_avg_wavelet_energy -= frame->frame_avg_wavelet_energy;
  section->coded_error -= frame->coded_error;
  section->sr_coded_error -= frame->sr_coded_error;
  section->pcnt_inter -= frame->pcnt_inter;
  section->pcnt_motion -= frame->pcnt_motion;
  section->pcnt_second_ref -= frame->pcnt_second_ref;
  section->pcnt_neutral -= frame->pcnt_neutral;
  section->intra_skip_pct -= frame->intra_skip_pct;
  section->inactive_zone_rows -= frame->inactive_zone_rows;
  section->inactive_zone_cols -= frame->inactive_zone_cols;
  section->MVr -= frame->MVr;
  section->mvr_abs -= frame->mvr_abs;
  section->MVc -= frame->MVc;
  section->mvc_abs -= frame->mvc_abs;
  section->MVrv -= frame->MVrv;
  section->MVcv -= frame->MVcv;
  section->mv_in_out_count -= frame->mv_in_out_count;
  section->new_mv_count -= frame->new_mv_count;
  section->count -= frame->count;
  section->duration -= frame->duration;
}

// This function returns the maximum target rate per frame.
static int frame_max_bits(const RATE_CONTROL *rc,
                          const EncodeContext *encode_context_ptr) {
  int64_t max_bits = ((int64_t)rc->avg_frame_bandwidth *
                      (int64_t)encode_context_ptr->two_pass_cfg.vbrmax_section) /
                     100;
  if (max_bits < 0)
    max_bits = 0;
  else if (max_bits > rc->max_frame_bandwidth)
    max_bits = rc->max_frame_bandwidth;

  return (int)max_bits;
}

static const double q_pow_term[(QINDEX_RANGE >> 5) + 1] = { 0.65, 0.70, 0.75,
                                                            0.80, 0.85, 0.90,
                                                            0.95, 0.95, 0.95 };
#define ERR_DIVISOR 96.0
static double calc_correction_factor(double err_per_mb, int q) {
  const double error_term = err_per_mb / ERR_DIVISOR;
  const int index = q >> 5;
  // Adjustment to power term based on qindex
  const double power_term =
      q_pow_term[index] +
      (((q_pow_term[index + 1] - q_pow_term[index]) * (q % 32)) / 32.0);
  assert(error_term >= 0.0);
  return fclamp(pow(error_term, power_term), 0.05, 5.0);
}

static void twopass_update_bpm_factor(TWO_PASS *twopass) {
  // Based on recent history adjust expectations of bits per macroblock.
  double last_group_rate_err =
      (double)twopass->rolling_arf_group_actual_bits /
      DOUBLE_DIVIDE_CHECK((double)twopass->rolling_arf_group_target_bits);
  last_group_rate_err = AOMMAX(0.25, AOMMIN(4.0, last_group_rate_err));
  twopass->bpm_factor *= (3.0 + last_group_rate_err) / 4.0;
  twopass->bpm_factor = AOMMAX(0.25, AOMMIN(4.0, twopass->bpm_factor));
}
static int qbpm_enumerator(int rate_err_tol) {
  return 1250000 + ((300000 * AOMMIN(75, AOMMAX(rate_err_tol - 25, 0))) / 75);
}

// Similar to find_qindex_by_rate() function in ratectrl.c, but includes
// calculation of a correction_factor.
static int find_qindex_by_rate_with_correction(
    int desired_bits_per_mb, aom_bit_depth_t bit_depth, double error_per_mb,
    double group_weight_factor, int rate_err_tol, int best_qindex,
    int worst_qindex) {
  assert(best_qindex <= worst_qindex);
  int low = best_qindex;
  int high = worst_qindex;

  while (low < high) {
    const int mid = (low + high) >> 1;
    const double mid_factor = calc_correction_factor(error_per_mb, mid);
    const double q = svt_av1_convert_qindex_to_q(mid, bit_depth);
    const int enumerator = qbpm_enumerator(rate_err_tol);
    const int mid_bits_per_mb =
        (int)((enumerator * mid_factor * group_weight_factor) / q);

    if (mid_bits_per_mb > desired_bits_per_mb) {
      low = mid + 1;
    } else {
      high = mid;
    }
  }
  return low;
}

// Bits Per MB at different Q (Multiplied by 512)
#define BPER_MB_NORMBITS 9

static int get_twopass_worst_quality(PictureParentControlSet *pcs_ptr, const double section_err,
                                     double inactive_zone,
                                     int section_target_bandwidth,
                                     double group_weight_factor) {
  SequenceControlSet *scs_ptr = pcs_ptr->scs_ptr;
  EncodeContext *encode_context_ptr = scs_ptr->encode_context_ptr;
  TWO_PASS *const twopass = &scs_ptr->twopass;
  RATE_CONTROL *const rc = &encode_context_ptr->rc;
  const RateControlCfg *const rc_cfg = &encode_context_ptr->rc_cfg;
  const uint32_t mb_cols = (scs_ptr->seq_header.max_frame_width  + 16 - 1) / 16;
  const uint32_t mb_rows = (scs_ptr->seq_header.max_frame_height + 16 - 1) / 16;
  inactive_zone = fclamp(inactive_zone, 0.0, 1.0);

  if (section_target_bandwidth <= 0) {
    return rc->worst_quality;  // Highest value allowed
  } else {
    const int num_mbs = mb_cols * mb_rows;
                        //(oxcf->resize_cfg.resize_mode != RESIZE_NONE)
                        //    ? cpi->initial_mbs
                        //    : cpi->common.mi_params.MBs;
    const int active_mbs = AOMMAX(1, num_mbs - (int)(num_mbs * inactive_zone));
    const double av_err_per_mb = section_err / active_mbs;
    const int target_norm_bits_per_mb =
        (int)((uint64_t)section_target_bandwidth << BPER_MB_NORMBITS) /
        active_mbs;
    int rate_err_tol = AOMMIN(rc_cfg->under_shoot_pct, rc_cfg->over_shoot_pct);

    twopass_update_bpm_factor(twopass);
    // Try and pick a max Q that will be high enough to encode the
    // content at the given rate.
    int q = find_qindex_by_rate_with_correction(
        target_norm_bits_per_mb, scs_ptr->encoder_bit_depth,
        av_err_per_mb, group_weight_factor, rate_err_tol, rc->best_quality,
        rc->worst_quality);

    // Restriction on active max q for constrained quality mode.
    //if (rc_cfg->mode == AOM_CQ) q = AOMMAX(q, rc_cfg->cq_level);
    return q;
  }
}

#define SR_DIFF_PART 0.0015
#define MOTION_AMP_PART 0.003
#define INTRA_PART 0.005
#define DEFAULT_DECAY_LIMIT 0.75
#define LOW_SR_DIFF_TRHESH 0.1
#define SR_DIFF_MAX 128.0
#define NCOUNT_FRAME_II_THRESH 5.0

static double get_sr_decay_rate(const FrameInfo *frame_info,
                                const FIRSTPASS_STATS *frame) {
  const int num_mbs = frame_info->num_mbs;
  double sr_diff = (frame->sr_coded_error - frame->coded_error) / num_mbs;
  double sr_decay = 1.0;
  double modified_pct_inter;
  double modified_pcnt_intra;
  const double motion_amplitude_factor =
      frame->pcnt_motion * ((frame->mvc_abs + frame->mvr_abs) / 2);

  modified_pct_inter = frame->pcnt_inter;
  if ((frame->intra_error / DOUBLE_DIVIDE_CHECK(frame->coded_error)) <
      (double)NCOUNT_FRAME_II_THRESH) {
    modified_pct_inter = frame->pcnt_inter - frame->pcnt_neutral;
  }
  modified_pcnt_intra = 100 * (1.0 - modified_pct_inter);

  if ((sr_diff > LOW_SR_DIFF_TRHESH)) {
    sr_diff = AOMMIN(sr_diff, SR_DIFF_MAX);
    sr_decay = 1.0 - (SR_DIFF_PART * sr_diff) -
               (MOTION_AMP_PART * motion_amplitude_factor) -
               (INTRA_PART * modified_pcnt_intra);
  }
  return AOMMAX(sr_decay, AOMMIN(DEFAULT_DECAY_LIMIT, modified_pct_inter));
}

// This function gives an estimate of how badly we believe the prediction
// quality is decaying from frame to frame.
static double get_zero_motion_factor(const FrameInfo *frame_info,
                                     const FIRSTPASS_STATS *frame) {
  const double zero_motion_pct = frame->pcnt_inter - frame->pcnt_motion;
  double sr_decay = get_sr_decay_rate(frame_info, frame);
  return AOMMIN(sr_decay, zero_motion_pct);
}

#define ZM_POWER_FACTOR 0.75

static double get_prediction_decay_rate(const FrameInfo *frame_info,
                                        const FIRSTPASS_STATS *next_frame) {
  const double sr_decay_rate = get_sr_decay_rate(frame_info, next_frame);
  const double zero_motion_factor =
      (0.95 * pow((next_frame->pcnt_inter - next_frame->pcnt_motion),
                  ZM_POWER_FACTOR));

  return AOMMAX(zero_motion_factor,
                (sr_decay_rate + ((1.0 - sr_decay_rate) * zero_motion_factor)));
}

// Function to test for a condition where a complex transition is followed
// by a static section. For example in slide shows where there is a fade
// between slides. This is to help with more optimal kf and gf positioning.
static int detect_transition_to_still(TWO_PASS *const twopass,
                                      const int min_gf_interval,
                                      const int frame_interval,
                                      const int still_interval,
                                      const double loop_decay_rate,
                                      const double last_decay_rate) {
  // Break clause to detect very still sections after motion
  // For example a static image after a fade or other transition
  // instead of a clean scene cut.
  if (frame_interval > min_gf_interval && loop_decay_rate >= 0.999 &&
      last_decay_rate < 0.9) {
    int j;
    // Look ahead a few frames to see if static condition persists...
    for (j = 0; j < still_interval; ++j) {
      const FIRSTPASS_STATS *stats = &twopass->stats_in[j];
      if (stats >= twopass->stats_buf_ctx->stats_in_end) break;

      if (stats->pcnt_inter - stats->pcnt_motion < 0.999) break;
    }
    // Only if it does do we signal a transition to still.
    return j == still_interval;
  }
  return 0;
}

// This function detects a flash through the high relative pcnt_second_ref
// score in the frame following a flash frame. The offset passed in should
// reflect this.
static int detect_flash(const TWO_PASS *twopass, const int offset) {
  const FIRSTPASS_STATS *const next_frame = read_frame_stats(twopass, offset);

  // What we are looking for here is a situation where there is a
  // brief break in prediction (such as a flash) but subsequent frames
  // are reasonably well predicted by an earlier (pre flash) frame.
  // The recovery after a flash is indicated by a high pcnt_second_ref
  // compared to pcnt_inter.
  return next_frame != NULL &&
         next_frame->pcnt_second_ref > next_frame->pcnt_inter &&
         next_frame->pcnt_second_ref >= 0.5;
}

// Update the motion related elements to the GF arf boost calculation.
static void accumulate_frame_motion_stats(const FIRSTPASS_STATS *stats,
                                          GF_GROUP_STATS *gf_stats) {
  const double pct = stats->pcnt_motion;

  // Accumulate Motion In/Out of frame stats.
  gf_stats->this_frame_mv_in_out = stats->mv_in_out_count * pct;
  gf_stats->mv_in_out_accumulator += gf_stats->this_frame_mv_in_out;
  gf_stats->abs_mv_in_out_accumulator += fabs(gf_stats->this_frame_mv_in_out);

  // Accumulate a measure of how uniform (or conversely how random) the motion
  // field is (a ratio of abs(mv) / mv).
  if (pct > 0.05) {
    const double mvr_ratio =
        fabs(stats->mvr_abs) / DOUBLE_DIVIDE_CHECK(fabs(stats->MVr));
    const double mvc_ratio =
        fabs(stats->mvc_abs) / DOUBLE_DIVIDE_CHECK(fabs(stats->MVc));

    gf_stats->mv_ratio_accumulator +=
        pct * (mvr_ratio < stats->mvr_abs ? mvr_ratio : stats->mvr_abs);
    gf_stats->mv_ratio_accumulator +=
        pct * (mvc_ratio < stats->mvc_abs ? mvc_ratio : stats->mvc_abs);
  }
}

static void accumulate_this_frame_stats(const FIRSTPASS_STATS *stats,
                                        const double mod_frame_err,
                                        GF_GROUP_STATS *gf_stats) {
  gf_stats->gf_group_err += mod_frame_err;
#if GROUP_ADAPTIVE_MAXQ
  gf_stats->gf_group_raw_error += stats->coded_error;
#endif
  gf_stats->gf_group_skip_pct += stats->intra_skip_pct;
  gf_stats->gf_group_inactive_zone_rows += stats->inactive_zone_rows;
}

static void accumulate_next_frame_stats(const FIRSTPASS_STATS *stats,
                                        const FrameInfo *frame_info,
                                        const int flash_detected,
                                        const int frames_since_key,
                                        const int cur_idx,
                                        GF_GROUP_STATS *gf_stats) {
  accumulate_frame_motion_stats(stats, gf_stats);
  // sum up the metric values of current gf group
  gf_stats->avg_sr_coded_error += stats->sr_coded_error;
  gf_stats->avg_tr_coded_error += stats->tr_coded_error;
  gf_stats->avg_pcnt_second_ref += stats->pcnt_second_ref;
  gf_stats->avg_pcnt_third_ref += stats->pcnt_third_ref;
  gf_stats->avg_new_mv_count += stats->new_mv_count;
  gf_stats->avg_wavelet_energy += stats->frame_avg_wavelet_energy;
  if (fabs(stats->raw_error_stdev) > 0.000001) {
    gf_stats->non_zero_stdev_count++;
    gf_stats->avg_raw_err_stdev += stats->raw_error_stdev;
  }

  // Accumulate the effect of prediction quality decay
  if (!flash_detected) {
    gf_stats->last_loop_decay_rate = gf_stats->loop_decay_rate;
    gf_stats->loop_decay_rate = get_prediction_decay_rate(frame_info, stats);

    gf_stats->decay_accumulator =
        gf_stats->decay_accumulator * gf_stats->loop_decay_rate;

    // Monitor for static sections.
    if ((frames_since_key + cur_idx - 1) > 1) {
      gf_stats->zero_motion_accumulator =
          AOMMIN(gf_stats->zero_motion_accumulator,
                 get_zero_motion_factor(frame_info, stats));
    }
  }
}

static void average_gf_stats(const int total_frame,
                             const FIRSTPASS_STATS *last_stat,
                             GF_GROUP_STATS *gf_stats) {
  if (total_frame) {
    gf_stats->avg_sr_coded_error /= total_frame;
    gf_stats->avg_tr_coded_error /= total_frame;
    gf_stats->avg_pcnt_second_ref /= total_frame;
    if (total_frame - 1) {
      gf_stats->avg_pcnt_third_ref_nolast =
          (gf_stats->avg_pcnt_third_ref - last_stat->pcnt_third_ref) /
          (total_frame - 1);
    } else {
      gf_stats->avg_pcnt_third_ref_nolast =
          gf_stats->avg_pcnt_third_ref / total_frame;
    }
    gf_stats->avg_pcnt_third_ref /= total_frame;
    gf_stats->avg_new_mv_count /= total_frame;
    gf_stats->avg_wavelet_energy /= total_frame;
  }

  if (gf_stats->non_zero_stdev_count)
    gf_stats->avg_raw_err_stdev /= gf_stats->non_zero_stdev_count;
}

#define BOOST_FACTOR 12.5
static double baseline_err_per_mb(const FrameInfo *frame_info) {
  unsigned int screen_area = frame_info->frame_height * frame_info->frame_width;

  // Use a different error per mb factor for calculating boost for
  //  different formats.
  if (screen_area <= 640 * 360) {
    return 500.0;
  } else {
    return 1000.0;
  }
}

static double calc_frame_boost(const RATE_CONTROL *rc,
                               const FrameInfo *frame_info,
                               const FIRSTPASS_STATS *this_frame,
                               double this_frame_mv_in_out, double max_boost) {
  double frame_boost;
  const double lq = svt_av1_convert_qindex_to_q(rc->avg_frame_qindex[INTER_FRAME],
                                            frame_info->bit_depth);
  const double boost_q_correction = AOMMIN((0.5 + (lq * 0.015)), 1.5);
  const double active_area = calculate_active_area(frame_info, this_frame);
  int num_mbs = frame_info->num_mbs;

  // Correct for any inactive region in the image
  num_mbs = (int)AOMMAX(1, num_mbs * active_area);

  // Underlying boost factor is based on inter error ratio.
  frame_boost = AOMMAX(baseline_err_per_mb(frame_info) * num_mbs,
                       this_frame->intra_error * active_area) /
                DOUBLE_DIVIDE_CHECK(this_frame->coded_error);
  frame_boost = frame_boost * BOOST_FACTOR * boost_q_correction;

  // Increase boost for frames where new data coming into frame (e.g. zoom out).
  // Slightly reduce boost if there is a net balance of motion out of the frame
  // (zoom in). The range for this_frame_mv_in_out is -1.0 to +1.0.
  if (this_frame_mv_in_out > 0.0)
    frame_boost += frame_boost * (this_frame_mv_in_out * 2.0);
  // In the extreme case the boost is halved.
  else
    frame_boost += frame_boost * (this_frame_mv_in_out / 2.0);

  return AOMMIN(frame_boost, max_boost * boost_q_correction);
}

static double calc_kf_frame_boost(const RATE_CONTROL *rc,
                                  const FrameInfo *frame_info,
                                  const FIRSTPASS_STATS *this_frame,
                                  double *sr_accumulator, double max_boost) {
  double frame_boost;
  const double lq = svt_av1_convert_qindex_to_q(rc->avg_frame_qindex[INTER_FRAME],
                                            frame_info->bit_depth);
  const double boost_q_correction = AOMMIN((0.50 + (lq * 0.015)), 2.00);
  const double active_area = calculate_active_area(frame_info, this_frame);
  int num_mbs = frame_info->num_mbs;

  // Correct for any inactive region in the image
  num_mbs = (int)AOMMAX(1, num_mbs * active_area);

  // Underlying boost factor is based on inter error ratio.
  frame_boost = AOMMAX(baseline_err_per_mb(frame_info) * num_mbs,
                       this_frame->intra_error * active_area) /
                DOUBLE_DIVIDE_CHECK(
                    (this_frame->coded_error + *sr_accumulator) * active_area);

  // Update the accumulator for second ref error difference.
  // This is intended to give an indication of how much the coded error is
  // increasing over time.
  *sr_accumulator += (this_frame->sr_coded_error - this_frame->coded_error);
  *sr_accumulator = AOMMAX(0.0, *sr_accumulator);

  // Q correction and scaling
  // The 40.0 value here is an experimentally derived baseline minimum.
  // This value is in line with the minimum per frame boost in the alt_ref
  // boost calculation.
  frame_boost = ((frame_boost + 40.0) * boost_q_correction);

  return AOMMIN(frame_boost, max_boost * boost_q_correction);
}

static int get_projected_gfu_boost(const RATE_CONTROL *rc, int gfu_boost,
                                   int frames_to_project,
                                   int num_stats_used_for_gfu_boost) {
  /*
   * If frames_to_project is equal to num_stats_used_for_gfu_boost,
   * it means that gfu_boost was calculated over frames_to_project to
   * begin with(ie; all stats required were available), hence return
   * the original boost.
   */
  if (num_stats_used_for_gfu_boost >= frames_to_project) return gfu_boost;

  double min_boost_factor = sqrt(rc->baseline_gf_interval);
  // Get the current tpl factor (number of frames = frames_to_project).
  double tpl_factor = svt_av1_get_gfu_boost_projection_factor(
      min_boost_factor, MAX_GFUBOOST_FACTOR, frames_to_project);
  // Get the tpl factor when number of frames = num_stats_used_for_prior_boost.
  double tpl_factor_num_stats = svt_av1_get_gfu_boost_projection_factor(
      min_boost_factor, MAX_GFUBOOST_FACTOR, num_stats_used_for_gfu_boost);
  int projected_gfu_boost =
      (int)rint((tpl_factor * gfu_boost) / tpl_factor_num_stats);
  return projected_gfu_boost;
}

#define GF_MAX_BOOST 90.0
#define MIN_DECAY_FACTOR 0.01
#define NORMAL_BOOST 100
static int av1_calc_arf_boost(const TWO_PASS *twopass, const RATE_CONTROL *rc,
                       FrameInfo *frame_info, int offset, int f_frames,
                       int b_frames, int *num_fpstats_used,
                       int *num_fpstats_required) {
  int i;
  GF_GROUP_STATS gf_stats;
  init_gf_stats(&gf_stats);
  double boost_score = (double)NORMAL_BOOST;
  int arf_boost;
  int flash_detected = 0;
  if (num_fpstats_used) *num_fpstats_used = 0;

  // Search forward from the proposed arf/next gf position.
  for (i = 0; i < f_frames; ++i) {
    const FIRSTPASS_STATS *this_frame = read_frame_stats(twopass, i + offset);
    if (this_frame == NULL) break;

    // Update the motion related elements to the boost calculation.
    accumulate_frame_motion_stats(this_frame, &gf_stats);

    // We want to discount the flash frame itself and the recovery
    // frame that follows as both will have poor scores.
    flash_detected = detect_flash(twopass, i + offset) ||
                     detect_flash(twopass, i + offset + 1);

    // Accumulate the effect of prediction quality decay.
    if (!flash_detected) {
      gf_stats.decay_accumulator *=
          get_prediction_decay_rate(frame_info, this_frame);
      gf_stats.decay_accumulator = gf_stats.decay_accumulator < MIN_DECAY_FACTOR
                                       ? MIN_DECAY_FACTOR
                                       : gf_stats.decay_accumulator;
    }

    boost_score +=
        gf_stats.decay_accumulator *
        calc_frame_boost(rc, frame_info, this_frame,
                         gf_stats.this_frame_mv_in_out, GF_MAX_BOOST);
    if (num_fpstats_used) (*num_fpstats_used)++;
  }

  arf_boost = (int)boost_score;

  // Reset for backward looking loop.
  boost_score = 0.0;
  init_gf_stats(&gf_stats);
  // Search backward towards last gf position.
  for (i = -1; i >= -b_frames; --i) {
    const FIRSTPASS_STATS *this_frame = read_frame_stats(twopass, i + offset);
    if (this_frame == NULL) break;

    // Update the motion related elements to the boost calculation.
    accumulate_frame_motion_stats(this_frame, &gf_stats);

    // We want to discount the the flash frame itself and the recovery
    // frame that follows as both will have poor scores.
    flash_detected = detect_flash(twopass, i + offset) ||
                     detect_flash(twopass, i + offset + 1);

    // Cumulative effect of prediction quality decay.
    if (!flash_detected) {
      gf_stats.decay_accumulator *=
          get_prediction_decay_rate(frame_info, this_frame);
      gf_stats.decay_accumulator = gf_stats.decay_accumulator < MIN_DECAY_FACTOR
                                       ? MIN_DECAY_FACTOR
                                       : gf_stats.decay_accumulator;
    }

    boost_score +=
        gf_stats.decay_accumulator *
        calc_frame_boost(rc, frame_info, this_frame,
                         gf_stats.this_frame_mv_in_out, GF_MAX_BOOST);
    if (num_fpstats_used) (*num_fpstats_used)++;
  }
  arf_boost += (int)boost_score;

  if (num_fpstats_required) {
    *num_fpstats_required = f_frames + b_frames;
    if (num_fpstats_used) {
      arf_boost = get_projected_gfu_boost(rc, arf_boost, *num_fpstats_required,
                                          *num_fpstats_used);
    }
  }

  if (arf_boost < ((b_frames + f_frames) * 50))
    arf_boost = ((b_frames + f_frames) * 50);

  return arf_boost;
}

// Calculate a section intra ratio used in setting max loop filter.
static int calculate_section_intra_ratio(const FIRSTPASS_STATS *begin,
                                         const FIRSTPASS_STATS *end,
                                         int section_length) {
  const FIRSTPASS_STATS *s = begin;
  double intra_error = 0.0;
  double coded_error = 0.0;
  int i = 0;

  while (s < end && i < section_length) {
    intra_error += s->intra_error;
    coded_error += s->coded_error;
    ++s;
    ++i;
  }

  return (int)(intra_error / DOUBLE_DIVIDE_CHECK(coded_error));
}

// Calculate the total bits to allocate in this GF/ARF group.
static int64_t calculate_total_gf_group_bits(PictureParentControlSet *pcs_ptr,
                                             double gf_group_err) {
  SequenceControlSet *scs_ptr = pcs_ptr->scs_ptr;
  EncodeContext *encode_context_ptr = scs_ptr->encode_context_ptr;
  RATE_CONTROL *const rc = &encode_context_ptr->rc;
  TWO_PASS *const twopass = &scs_ptr->twopass;
  const int max_bits = frame_max_bits(rc, encode_context_ptr);
  int64_t total_group_bits;

  if (scs_ptr->lap_enabled) {
    total_group_bits = rc->avg_frame_bandwidth * rc->baseline_gf_interval;
    return total_group_bits;
  }

  // Calculate the bits to be allocated to the group as a whole.
  if ((twopass->kf_group_bits > 0) && (twopass->kf_group_error_left > 0)) {
    total_group_bits = (int64_t)(twopass->kf_group_bits *
                                 (gf_group_err / twopass->kf_group_error_left));
  } else {
    total_group_bits = 0;
  }

  // Clamp odd edge cases.
  total_group_bits = (total_group_bits < 0)
                         ? 0
                         : (total_group_bits > twopass->kf_group_bits)
                               ? twopass->kf_group_bits
                               : total_group_bits;

  // Clip based on user supplied data rate variability limit.
  if (total_group_bits > (int64_t)max_bits * rc->baseline_gf_interval)
    total_group_bits = (int64_t)max_bits * rc->baseline_gf_interval;

  return total_group_bits;
}

// Calculate the number of bits to assign to boosted frames in a group.
static int calculate_boost_bits(int frame_count, int boost,
                                int64_t total_group_bits) {
  int allocation_chunks;

  // return 0 for invalid inputs (could arise e.g. through rounding errors)
  if (!boost || (total_group_bits <= 0)) return 0;

  if (frame_count <= 0) return (int)(AOMMIN(total_group_bits, INT_MAX));

  allocation_chunks = (frame_count * 100) + boost;

  // Prevent overflow.
  if (boost > 1023) {
    int divisor = boost >> 10;
    boost /= divisor;
    allocation_chunks /= divisor;
  }

  // Calculate the number of extra bits for use in the boosted frame or frames.
  return AOMMAX((int)(((int64_t)boost * total_group_bits) / allocation_chunks),
                0);
}

// Allocate bits to each frame in a GF / ARF group
static double layer_fraction[MAX_ARF_LAYERS + 1] = { 1.0,  0.70, 0.55, 0.60,
                                              0.60, 1.0,  1.0 };
static void allocate_gf_group_bits(GF_GROUP *gf_group, RATE_CONTROL *const rc,
                                   int64_t gf_group_bits, int gf_arf_bits,
                                   int key_frame, int use_arf) {
  int64_t total_group_bits = gf_group_bits;
  int base_frame_bits;
  const int gf_group_size = gf_group->size;
  int layer_frames[MAX_ARF_LAYERS + 1] = { 0 };

  // Subtract the extra bits set aside for ARF frames from the Group Total
  if (use_arf || !key_frame) total_group_bits -= gf_arf_bits;

  if (rc->baseline_gf_interval)
    base_frame_bits = (int)(total_group_bits / rc->baseline_gf_interval);
  else
    base_frame_bits = (int)1;

  // For key frames the frame target rate is already set and it
  // is also the golden frame.
  // === [frame_index == 0] ===
  int frame_index = 0;
  if (!key_frame) {
    if (rc->source_alt_ref_active)
      gf_group->bit_allocation[frame_index] = 0;
    else
      gf_group->bit_allocation[frame_index] =
          base_frame_bits + (int)(gf_arf_bits * layer_fraction[1]);
  }
  frame_index++;

  // Check the number of frames in each layer in case we have a
  // non standard group length.
  int max_arf_layer = gf_group->max_layer_depth - 1;
  for (int idx = frame_index; idx < gf_group_size; ++idx) {
    if ((gf_group->update_type[idx] == ARF_UPDATE) ||
        (gf_group->update_type[idx] == INTNL_ARF_UPDATE)) {
      layer_frames[gf_group->layer_depth[idx]]++;
    }
  }

  // Allocate extra bits to each ARF layer
  int i;
  int layer_extra_bits[MAX_ARF_LAYERS + 1] = { 0 };
  for (i = 1; i <= max_arf_layer; ++i) {
      if (layer_frames[i]) {// to make sure there is a picture with the depth
          double fraction = (i == max_arf_layer) ? 1.0 : layer_fraction[i];
          layer_extra_bits[i] =
              (int)((gf_arf_bits * fraction) / AOMMAX(1, layer_frames[i]));
          gf_arf_bits -= (int)(gf_arf_bits * fraction);
      }
  }

  // Now combine ARF layer and baseline bits to give total bits for each frame.
  int arf_extra_bits;
  for (int idx = frame_index; idx < gf_group_size; ++idx) {
    switch (gf_group->update_type[idx]) {
      case ARF_UPDATE:
      case INTNL_ARF_UPDATE:
        arf_extra_bits = layer_extra_bits[gf_group->layer_depth[idx]];
        gf_group->bit_allocation[idx] = base_frame_bits + arf_extra_bits;
        break;
      case INTNL_OVERLAY_UPDATE:
      case OVERLAY_UPDATE: gf_group->bit_allocation[idx] = 0; break;
      default: gf_group->bit_allocation[idx] = base_frame_bits; break;
    }
  }

  // Set the frame following the current GOP to 0 bit allocation. For ARF
  // groups, this next frame will be overlay frame, which is the first frame
  // in the next GOP. For GF group, next GOP will overwrite the rate allocation.
  // Setting this frame to use 0 bit (of out the current GOP budget) will
  // simplify logics in reference frame management.
  gf_group->bit_allocation[gf_group_size] = 0;
}

// Returns true if KF group and GF group both are almost completely static.
static INLINE int is_almost_static(double gf_zero_motion, int kf_zero_motion,
                                   int is_lap_enabled) {
  if (is_lap_enabled) {
    /*
     * when LAP enabled kf_zero_motion is not reliable, so use strict
     * constraint on gf_zero_motion.
     */
    return (gf_zero_motion >= 0.999);
  } else {
    return (gf_zero_motion >= 0.995) &&
           (kf_zero_motion >= STATIC_KF_GROUP_THRESH);
  }
}

#if GROUP_ADAPTIVE_MAXQ
#define RC_FACTOR_MIN 0.75
#define RC_FACTOR_MAX 1.25
#endif  // GROUP_ADAPTIVE_MAXQ
#define MIN_FWD_KF_INTERVAL 8

// This function imposes the gf group length of future frames in batch based on the intra refresh
// only supports for 5L
static void impose_gf_length(PictureParentControlSet *pcs_ptr, int max_intervals) {
    SequenceControlSet *scs_ptr            = pcs_ptr->scs_ptr;
    EncodeContext *     encode_context_ptr = scs_ptr->encode_context_ptr;
    RATE_CONTROL *const rc                 = &encode_context_ptr->rc;
    int                 i                  = 0;
    max_intervals                          = scs_ptr->lap_enabled ? 1 : max_intervals;
    int cut_pos[MAX_NUM_GF_INTERVALS + 1]  = {0};
    int count_cuts                         = 1;
    int cur_last;
    while (count_cuts < max_intervals + 1) {
        int cut_here;
        ++i;
        // reaches next key frame, break here
        if (i >= rc->frames_to_key) {
            if (cut_pos[count_cuts - 1] != i - 1) {
                cut_pos[count_cuts] = i - 1;
                count_cuts++;
            }
            break;
        }
        // To cut based on PD decisions, only supports 5L for now
        cut_here =
            ((i % 16 == 0) || ((rc->frames_to_key - cut_pos[count_cuts - 1]) < 16 && (i % 8 == 0)))
                ? 1 : 0;
        if (cut_here) {
            cur_last            = i; // the current last frame in the gf group
            cut_pos[count_cuts] = cur_last;
            count_cuts++;
        }
    }
    // save intervals
    rc->intervals_till_gf_calculate_due = count_cuts - 1;
    for (int n = 1; n < count_cuts; n++) {
        rc->gf_intervals[n - 1] = cut_pos[n] + 1 - cut_pos[n - 1];
    }
    rc->cur_gf_index = 0;
}
static INLINE void set_baseline_gf_interval(PictureParentControlSet *pcs_ptr, int arf_position,
                                            int active_max_gf_interval,
                                            int use_alt_ref,
                                            int is_final_pass) {
  SequenceControlSet *scs_ptr = pcs_ptr->scs_ptr;
  EncodeContext *encode_context_ptr = scs_ptr->encode_context_ptr;
  RATE_CONTROL *const rc = &encode_context_ptr->rc;
  TWO_PASS *const twopass = &scs_ptr->twopass;
  // Set the interval until the next gf.
  // If forward keyframes are enabled, ensure the final gf group obeys the
  // MIN_FWD_KF_INTERVAL.
  const int is_last_kf =
      (twopass->stats_in - arf_position + rc->frames_to_key) >=
      twopass->stats_buf_ctx->stats_in_end;

  if (encode_context_ptr->kf_cfg.fwd_kf_enabled && use_alt_ref && !is_last_kf &&
      rc->next_is_fwd_key) {
    if (arf_position == rc->frames_to_key) {
      rc->baseline_gf_interval = arf_position;
      // if the last gf group will be smaller than MIN_FWD_KF_INTERVAL
    } else if ((rc->frames_to_key - arf_position <
                AOMMAX(MIN_FWD_KF_INTERVAL, rc->min_gf_interval)) &&
               (rc->frames_to_key != arf_position)) {
      // if possible, merge the last two gf groups
      if (rc->frames_to_key <= active_max_gf_interval) {
        rc->baseline_gf_interval = rc->frames_to_key;
        if (is_final_pass) rc->intervals_till_gf_calculate_due = 0;
        // if merging the last two gf groups creates a group that is too long,
        // split them and force the last gf group to be the MIN_FWD_KF_INTERVAL
      } else {
        rc->baseline_gf_interval = rc->frames_to_key - MIN_FWD_KF_INTERVAL;
        if (is_final_pass) rc->intervals_till_gf_calculate_due = 0;
      }
    } else {
      rc->baseline_gf_interval = arf_position - rc->source_alt_ref_pending;
    }
  } else {
    rc->baseline_gf_interval = arf_position - rc->source_alt_ref_pending;
  }
}

// initialize GF_GROUP_STATS
static void init_gf_stats(GF_GROUP_STATS *gf_stats) {
  gf_stats->gf_group_err = 0.0;
  gf_stats->gf_group_raw_error = 0.0;
  gf_stats->gf_group_skip_pct = 0.0;
  gf_stats->gf_group_inactive_zone_rows = 0.0;

  gf_stats->mv_ratio_accumulator = 0.0;
  gf_stats->decay_accumulator = 1.0;
  gf_stats->zero_motion_accumulator = 1.0;
  gf_stats->loop_decay_rate = 1.0;
  gf_stats->last_loop_decay_rate = 1.0;
  gf_stats->this_frame_mv_in_out = 0.0;
  gf_stats->mv_in_out_accumulator = 0.0;
  gf_stats->abs_mv_in_out_accumulator = 0.0;

  gf_stats->avg_sr_coded_error = 0.0;
  gf_stats->avg_tr_coded_error = 0.0;
  gf_stats->avg_pcnt_second_ref = 0.0;
  gf_stats->avg_pcnt_third_ref = 0.0;
  gf_stats->avg_pcnt_third_ref_nolast = 0.0;
  gf_stats->avg_new_mv_count = 0.0;
  gf_stats->avg_wavelet_energy = 0.0;
  gf_stats->avg_raw_err_stdev = 0.0;
  gf_stats->non_zero_stdev_count = 0;
}
// function from gop_structure.c
// Set parameters for frames between 'start' and 'end' (excluding both).
static void set_multi_layer_params(const TWO_PASS *twopass,
    GF_GROUP *const gf_group, RATE_CONTROL *rc,
    FrameInfo *frame_info, int start, int end,
    int *cur_frame_idx, int *frame_ind,
    int layer_depth) {
    const int num_frames_to_process = end - start - 1;
    assert(num_frames_to_process >= 0);
    if (num_frames_to_process == 0) return;

    // Either we are at the last level of the pyramid, or we don't have enough
    // frames between 'l' and 'r' to create one more level.
    if (layer_depth > gf_group->max_layer_depth_allowed ||
        num_frames_to_process < 3) {
        // Leaf nodes.
        while (++start < end) {
            gf_group->update_type[*frame_ind] = LF_UPDATE;
            gf_group->arf_src_offset[*frame_ind] = 0;
            ++*cur_frame_idx;
            gf_group->cur_frame_idx[*frame_ind] = *cur_frame_idx;
            gf_group->frame_disp_idx[*frame_ind] = start;
            gf_group->layer_depth[*frame_ind] = MAX_ARF_LAYERS;
            gf_group->arf_boost[*frame_ind] = av1_calc_arf_boost(
                twopass, rc, frame_info, start, end - start, 0, NULL, NULL);
            gf_group->max_layer_depth =
                AOMMAX(gf_group->max_layer_depth, layer_depth);
            ++(*frame_ind);
        }
    }
    else {
        const int m = (start + end) / 2;

        // Internal ARF.
        gf_group->update_type[*frame_ind] = INTNL_ARF_UPDATE;
        gf_group->arf_src_offset[*frame_ind] = m - start - 1;
        gf_group->cur_frame_idx[*frame_ind] = *cur_frame_idx;
        gf_group->frame_disp_idx[*frame_ind] = m;
        gf_group->layer_depth[*frame_ind] = layer_depth;

        // Get the boost factor for intermediate ARF frames.
        gf_group->arf_boost[*frame_ind] = av1_calc_arf_boost(
            twopass, rc, frame_info, m, end - m, m - start, NULL, NULL);
        ++(*frame_ind);

        // Frames displayed before this internal ARF.
        set_multi_layer_params(twopass, gf_group, rc, frame_info, start, m,
            cur_frame_idx, frame_ind, layer_depth + 1);

        // Frames displayed after this internal ARF.
        set_multi_layer_params(twopass, gf_group, rc, frame_info, m, end,
            cur_frame_idx, frame_ind, layer_depth + 1);
    }
}
static int construct_multi_layer_gf_structure(
    TWO_PASS *twopass, GF_GROUP *const gf_group,
    RATE_CONTROL *rc, FrameInfo *const frame_info, int gf_interval,
    FRAME_UPDATE_TYPE first_frame_update_type) {
    int frame_index = 0;
    int cur_frame_index = 0;

    // Keyframe / Overlay frame / Golden frame.
    assert(gf_interval >= 1);
    assert(first_frame_update_type == KF_UPDATE ||
        first_frame_update_type == OVERLAY_UPDATE ||
        first_frame_update_type == GF_UPDATE);

    gf_group->update_type[frame_index] = first_frame_update_type;
    gf_group->arf_src_offset[frame_index] = 0;
    ++cur_frame_index;
    gf_group->cur_frame_idx[frame_index] = cur_frame_index;
    gf_group->layer_depth[frame_index] =
        first_frame_update_type == OVERLAY_UPDATE ? MAX_ARF_LAYERS + 1 : 0;
    gf_group->max_layer_depth = 0;
    ++frame_index;

    // anaghdin: for now only 5L is supported. In 5L case, when there are not enough picture,
    // we switch to 4L and after that we use 4L P pictures. In the else, we handle the P-case manually
    // this logic has to move to picture decision
    if (gf_interval >= 8) {
        // ALTREF.
        const int use_altref = gf_group->max_layer_depth_allowed > 0;
        if (use_altref) {
            gf_group->update_type[frame_index] = ARF_UPDATE;
            gf_group->arf_src_offset[frame_index] = gf_interval - 1;
            gf_group->cur_frame_idx[frame_index] = cur_frame_index;
            gf_group->frame_disp_idx[frame_index] = gf_interval;
            gf_group->layer_depth[frame_index] = 1;
            gf_group->arf_boost[frame_index] = rc->gfu_boost;
            gf_group->max_layer_depth = 1;
            ++frame_index;
        }
        // Rest of the frames.
        set_multi_layer_params(twopass, gf_group, rc, frame_info, 0, gf_interval,
            &cur_frame_index, &frame_index, use_altref + 1);
    }
    else {

        int start = 0;
        int end = gf_interval;
        //const int num_frames_to_process = end - start - 1;
        while (++start <= end) {
            gf_group->update_type[frame_index] = (frame_index % 2 == 0) ? INTNL_ARF_UPDATE : LF_UPDATE;
            gf_group->arf_src_offset[frame_index] = 0;
            gf_group->cur_frame_idx[frame_index] = start;
            gf_group->frame_disp_idx[frame_index] = start;
            gf_group->layer_depth[frame_index] = (frame_index % 4 == 0) ? 2 :
                (frame_index % 2 == 0) ? 3 : MAX_ARF_LAYERS;
            gf_group->arf_boost[frame_index] = av1_calc_arf_boost(
                twopass, rc, frame_info, start, end - start, 0, NULL, NULL);
            gf_group->max_layer_depth =
                AOMMAX(gf_group->max_layer_depth, gf_group->layer_depth[frame_index]);
            ++(frame_index);
        }
    }
    return frame_index;
}

static void av1_gop_setup_structure(PictureParentControlSet *pcs_ptr,
    const EncodeFrameParams *const frame_params) {
        SequenceControlSet *scs_ptr = pcs_ptr->scs_ptr;
        EncodeContext *encode_context_ptr = scs_ptr->encode_context_ptr;
        RATE_CONTROL *const rc = &encode_context_ptr->rc;
        TWO_PASS *const twopass = &scs_ptr->twopass;
        GF_GROUP *const gf_group = &encode_context_ptr->gf_group;
        FrameInfo *frame_info = &encode_context_ptr->frame_info;


    const int key_frame = (frame_params->frame_type == KEY_FRAME);
    const FRAME_UPDATE_TYPE first_frame_update_type =
        key_frame ? KF_UPDATE
        : rc->source_alt_ref_active ? OVERLAY_UPDATE : GF_UPDATE;
    gf_group->size = construct_multi_layer_gf_structure(
        twopass, gf_group, rc, frame_info, rc->baseline_gf_interval,
        first_frame_update_type);

#if CHECK_GF_PARAMETER
    check_frame_params(gf_group, rc->baseline_gf_interval);
#endif
}

static void av1_gop_bit_allocation(RATE_CONTROL *const rc,
                            GF_GROUP *gf_group, int is_key_frame, int use_arf,
                            int64_t gf_group_bits);
int frame_is_kf_gf_arf(PictureParentControlSet *ppcs_ptr);
// Analyse and define a gf/arf group.
#define MAX_GF_BOOST 5400
static void define_gf_group(PictureParentControlSet *pcs_ptr, FIRSTPASS_STATS *this_frame,
                            const EncodeFrameParams *const frame_params,
                            int max_gop_length, int is_final_pass) {
  SequenceControlSet *scs_ptr = pcs_ptr->scs_ptr;
  EncodeContext *encode_context_ptr = scs_ptr->encode_context_ptr;
  RATE_CONTROL *const rc = &encode_context_ptr->rc;
  TWO_PASS *const twopass = &scs_ptr->twopass;
  FIRSTPASS_STATS next_frame;
  const FIRSTPASS_STATS *const start_pos = twopass->stats_in;
  GF_GROUP *const gf_group = &encode_context_ptr->gf_group;
  FrameInfo *frame_info = &encode_context_ptr->frame_info;
  const GFConfig *const gf_cfg = &encode_context_ptr->gf_cfg;
  const RateControlCfg *const rc_cfg = &encode_context_ptr->rc_cfg;
  int i;

  int64_t gf_group_bits;
  const int is_intra_only = frame_params->frame_type == KEY_FRAME ||
                            frame_params->frame_type == INTRA_ONLY_FRAME;

  pcs_ptr->internal_altref_allowed = (gf_cfg->gf_max_pyr_height > 1);

  // Reset the GF group data structures unless this is a key
  // frame in which case it will already have been done.
  if (!is_intra_only) {
    av1_zero(encode_context_ptr->gf_group);
    pcs_ptr->gf_group_index = 0;
  }
#ifdef ARCH_X86_64
  aom_clear_system_state();
#endif
  av1_zero(next_frame);

  GF_GROUP_STATS gf_stats;
  init_gf_stats(&gf_stats);
  GF_FRAME_STATS first_frame_stats;

  const int can_disable_arf = (gf_cfg->gf_min_pyr_height == MIN_PYRAMID_LVL);

  // Load stats for the current frame.
  double mod_frame_err =
      calculate_modified_err(frame_info, twopass, &(encode_context_ptr->two_pass_cfg), this_frame);

  // Note the error of the frame at the start of the group. This will be
  // the GF frame error if we code a normal gf.
  first_frame_stats.frame_err = mod_frame_err;
  first_frame_stats.frame_coded_error = this_frame->coded_error;
  first_frame_stats.frame_sr_coded_error = this_frame->sr_coded_error;
  first_frame_stats.frame_tr_coded_error = this_frame->tr_coded_error;

  // If this is a key frame or the overlay from a previous arf then
  // the error score / cost of this frame has already been accounted for.
  // There is no overlay support for now
  if (is_intra_only) {
    gf_stats.gf_group_err -= first_frame_stats.frame_err;
#if GROUP_ADAPTIVE_MAXQ
    gf_stats.gf_group_raw_error -= this_frame->coded_error;
#endif
    gf_stats.gf_group_skip_pct -= this_frame->intra_skip_pct;
    gf_stats.gf_group_inactive_zone_rows -= this_frame->inactive_zone_rows;
  }

  // TODO(urvang): Try logic to vary min and max interval based on q.
  //const int active_min_gf_interval = rc->min_gf_interval;
  const int active_max_gf_interval =
      AOMMIN(rc->max_gf_interval, max_gop_length);
  // If the first frame is not key frame, we start from i=1
  if (is_intra_only)
      i = 0;
  else
      i = 1;
  // get the determined gf group length from rc->gf_intervals
  while (i < rc->gf_intervals[rc->cur_gf_index]) {
    int flash_detected;
    ++i;
    // Accumulate error score of frames in this gf group.
    mod_frame_err =
        calculate_modified_err(frame_info, twopass, &(encode_context_ptr->two_pass_cfg), this_frame);
    // accumulate stats for this frame
    accumulate_this_frame_stats(this_frame, mod_frame_err, &gf_stats);

    // read in the next frame
    if (EOF == input_stats(twopass, &next_frame)) break;

    // Test for the case where there is a brief flash but the prediction
    // quality back to an earlier frame is then restored.
    flash_detected = detect_flash(twopass, 0);

    // accumulate stats for next frame
    accumulate_next_frame_stats(&next_frame, frame_info, flash_detected,
                                rc->frames_since_key, i, &gf_stats);

    *this_frame = next_frame;
  }

  if (is_final_pass) {
    rc->intervals_till_gf_calculate_due--;
    rc->cur_gf_index++;
  }

  // Was the group length constrained by the requirement for a new KF?
  rc->constrained_gf_group = (i >= rc->frames_to_key) ? 1 : 0;

  const int num_mbs = frame_info->num_mbs;
                      //(oxcf->resize_cfg.resize_mode != RESIZE_NONE)
                      //    ? cpi->initial_mbs
                      //    : cm->mi_params.MBs;
  assert(num_mbs > 0);

  average_gf_stats(i, &next_frame, &gf_stats);

  // Disable internal ARFs for "still" gf groups.
  //   zero_motion_accumulator: minimum percentage of (0,0) motion;
  //   avg_sr_coded_error:      average of the SSE per pixel of each frame;
  //   avg_raw_err_stdev:       average of the standard deviation of (0,0)
  //                            motion error per block of each frame.
  const int can_disable_internal_arfs =
      (gf_cfg->gf_min_pyr_height <= MIN_PYRAMID_LVL + 1);
  if (can_disable_internal_arfs &&
      gf_stats.zero_motion_accumulator > MIN_ZERO_MOTION &&
      gf_stats.avg_sr_coded_error / num_mbs < MAX_SR_CODED_ERROR &&
      gf_stats.avg_raw_err_stdev < MAX_RAW_ERR_VAR) {
    pcs_ptr->internal_altref_allowed = 0;
  }

  int use_alt_ref;
  if (can_disable_arf) {
    use_alt_ref =
        !is_almost_static(gf_stats.zero_motion_accumulator,
                          twopass->kf_zeromotion_pct, scs_ptr->lap_enabled) &&
        rc->use_arf_in_this_kf_group && (i < gf_cfg->lag_in_frames) &&
        (i >= MIN_GF_INTERVAL) && (gf_cfg->gf_max_pyr_height > MIN_PYRAMID_LVL);

  } else {
    assert(gf_cfg->gf_max_pyr_height > MIN_PYRAMID_LVL);
    use_alt_ref =
        rc->use_arf_in_this_kf_group && (i < gf_cfg->lag_in_frames) && (i > 2);
  }

#define REDUCE_GF_LENGTH_THRESH 4
#define REDUCE_GF_LENGTH_TO_KEY_THRESH 9
#define REDUCE_GF_LENGTH_BY 1
  int alt_offset = 0; // left to add code from 2-pass

  // Should we use the alternate reference frame.
  if (use_alt_ref) {
    rc->source_alt_ref_pending = 1;
    gf_group->max_layer_depth_allowed = gf_cfg->gf_max_pyr_height;
    // Get from actual minigop size in PD
    set_baseline_gf_interval(pcs_ptr, i, active_max_gf_interval, use_alt_ref,
                             is_final_pass);

    const int forward_frames = (rc->frames_to_key - i >= i - 1)
                                   ? i - 1
                                   : AOMMAX(0, rc->frames_to_key - i);

    // Calculate the boost for alt ref.
    rc->gfu_boost = av1_calc_arf_boost(
        twopass, rc, frame_info, alt_offset, forward_frames, (i - 1),
        scs_ptr->lap_enabled ? &rc->num_stats_used_for_gfu_boost : NULL,
        scs_ptr->lap_enabled ? &rc->num_stats_required_for_gfu_boost : NULL);
  } else {
    reset_fpf_position(twopass, start_pos);
    rc->source_alt_ref_pending = 0;
    gf_group->max_layer_depth_allowed = 0;
    set_baseline_gf_interval(pcs_ptr, i, active_max_gf_interval, use_alt_ref,
                             is_final_pass);

    rc->gfu_boost = AOMMIN(
        MAX_GF_BOOST,
        av1_calc_arf_boost(
            twopass, rc, frame_info, alt_offset, (i - 1), 0,
            scs_ptr->lap_enabled ? &rc->num_stats_used_for_gfu_boost : NULL,
            scs_ptr->lap_enabled ? &rc->num_stats_required_for_gfu_boost : NULL));
  }

  // rc->gf_intervals assumes the usage of alt_ref, therefore adding one overlay
  // frame to the next gf. If no alt_ref is used, should substract 1 frame from
  // the next gf group.
  // TODO(bohanli): should incorporate the usage of alt_ref into
  // calculate_gf_length
  if (is_final_pass && rc->source_alt_ref_pending == 0 &&
      rc->intervals_till_gf_calculate_due > 0) {
    rc->gf_intervals[rc->cur_gf_index]--;
  }

#define LAST_ALR_BOOST_FACTOR 0.2f
  rc->arf_boost_factor = 1.0;
  if (rc->source_alt_ref_pending && !is_lossless_requested(rc_cfg)) {
    // Reduce the boost of altref in the last gf group
    if (rc->frames_to_key - i == REDUCE_GF_LENGTH_BY ||
        rc->frames_to_key - i == 0) {
      rc->arf_boost_factor = LAST_ALR_BOOST_FACTOR;
    }
  }

  rc->frames_till_gf_update_due = rc->baseline_gf_interval;

  // Reset the file position.
  reset_fpf_position(twopass, start_pos);

  // Calculate the bits to be allocated to the gf/arf group as a whole
  gf_group_bits = calculate_total_gf_group_bits(pcs_ptr, gf_stats.gf_group_err);
  rc->gf_group_bits = gf_group_bits;

#if GROUP_ADAPTIVE_MAXQ
  // Calculate an estimate of the maxq needed for the group.
  // We are more agressive about correcting for sections
  // where there could be significant overshoot than for easier
  // sections where we do not wish to risk creating an overshoot
  // of the allocated bit budget.
  if ((rc_cfg->mode != AOM_Q) && (rc->baseline_gf_interval > 1)) {
    const int vbr_group_bits_per_frame =
        (int)(gf_group_bits / rc->baseline_gf_interval);
    const double group_av_err =
        gf_stats.gf_group_raw_error / rc->baseline_gf_interval;
    const double group_av_skip_pct =
        gf_stats.gf_group_skip_pct / rc->baseline_gf_interval;
    const double group_av_inactive_zone =
        ((gf_stats.gf_group_inactive_zone_rows * 2) /
         (rc->baseline_gf_interval * (double)frame_info->mb_rows));

    int tmp_q;
    // rc factor is a weight factor that corrects for local rate control drift.
    double rc_factor = 1.0;
    int64_t bits = scs_ptr->static_config.target_bit_rate;//oxcf->target_bandwidth;

    if (bits > 0) {
      int rate_error;

      rate_error = (int)((rc->vbr_bits_off_target * 100) / bits);
      rate_error = clamp(rate_error, -100, 100);
      if (rate_error > 0) {
        rc_factor = AOMMAX(RC_FACTOR_MIN, (double)(100 - rate_error) / 100.0);
      } else {
        rc_factor = AOMMIN(RC_FACTOR_MAX, (double)(100 - rate_error) / 100.0);
      }
    }
    tmp_q = get_twopass_worst_quality(
        pcs_ptr, group_av_err, (group_av_skip_pct + group_av_inactive_zone),
        vbr_group_bits_per_frame, rc_factor);
    rc->active_worst_quality = AOMMAX(tmp_q, rc->active_worst_quality >> 1);
  }
#endif

  // Adjust KF group bits and error remaining.
  if (is_final_pass)
    twopass->kf_group_error_left -= (int64_t)gf_stats.gf_group_err;
  // Set up the structure of this Group-Of-Pictures (same as GF_GROUP)
  av1_gop_setup_structure(pcs_ptr, frame_params);

  // Reset the file position.
  reset_fpf_position(twopass, start_pos);

  // Calculate a section intra ratio used in setting max loop filter.
  if (frame_params->frame_type != KEY_FRAME) {
    twopass->section_intra_rating = calculate_section_intra_ratio(
        start_pos, twopass->stats_buf_ctx->stats_in_end,
        rc->baseline_gf_interval);
  }

  // Reset rolling actual and target bits counters for ARF groups.
  twopass->rolling_arf_group_target_bits = 1;
  twopass->rolling_arf_group_actual_bits = 1;

  av1_gop_bit_allocation(
      rc, gf_group, frame_params->frame_type == KEY_FRAME, use_alt_ref, gf_group_bits);
}

// #define FIXED_ARF_BITS
#ifdef FIXED_ARF_BITS
#define ARF_BITS_FRACTION 0.75
#endif
static void av1_gop_bit_allocation(RATE_CONTROL *const rc, GF_GROUP *gf_group, int is_key_frame,
                            int use_arf, int64_t gf_group_bits) {
    // Calculate the extra bits to be used for boosted frame(s)
#ifdef FIXED_ARF_BITS
  int gf_arf_bits = (int)(ARF_BITS_FRACTION * gf_group_bits);
#else
  int gf_arf_bits = calculate_boost_bits(rc->baseline_gf_interval,
                                         rc->gfu_boost, gf_group_bits);
#endif

  // Allocate bits to each of the frames in the GF group.
  allocate_gf_group_bits(gf_group, rc, gf_group_bits, gf_arf_bits, is_key_frame,
                         use_arf);
}

// Minimum % intra coding observed in first pass (1.0 = 100%)
#define MIN_INTRA_LEVEL 0.25
// Minimum ratio between the % of intra coding and inter coding in the first
// pass after discounting neutral blocks (discounting neutral blocks in this
// way helps catch scene cuts in clips with very flat areas or letter box
// format clips with image padding.
#define INTRA_VS_INTER_THRESH 2.0
// Hard threshold where the first pass chooses intra for almost all blocks.
// In such a case even if the frame is not a scene cut coding a key frame
// may be a good option.
#define VERY_LOW_INTER_THRESH 0.05
// Maximum threshold for the relative ratio of intra error score vs best
// inter error score.
#define KF_II_ERR_THRESHOLD 2.5
// In real scene cuts there is almost always a sharp change in the intra
// or inter error score.
#define ERR_CHANGE_THRESHOLD 0.4
// For real scene cuts we expect an improvment in the intra inter error
// ratio in the next frame.
#define II_IMPROVEMENT_THRESHOLD 3.5
#define KF_II_MAX 128.0

// Threshold for use of the lagging second reference frame. High second ref
// usage may point to a transient event like a flash or occlusion rather than
// a real scene cut.
// We adapt the threshold based on number of frames in this key-frame group so
// far.
static double get_second_ref_usage_thresh(int frame_count_so_far) {
  const int adapt_upto = 32;
  const double min_second_ref_usage_thresh = 0.085;
  const double second_ref_usage_thresh_max_delta = 0.035;
  if (frame_count_so_far >= adapt_upto) {
    return min_second_ref_usage_thresh + second_ref_usage_thresh_max_delta;
  }
  return min_second_ref_usage_thresh +
         ((double)frame_count_so_far / (adapt_upto - 1)) *
             second_ref_usage_thresh_max_delta;
}

static int test_candidate_kf(TWO_PASS *twopass,
                             const FIRSTPASS_STATS *last_frame,
                             const FIRSTPASS_STATS *this_frame,
                             const FIRSTPASS_STATS *next_frame,
                             int frame_count_so_far, enum aom_rc_mode rc_mode,
                             int scenecut_mode) {
  int is_viable_kf = 0;
  double pcnt_intra = 1.0 - this_frame->pcnt_inter;
  double modified_pcnt_inter =
      this_frame->pcnt_inter - this_frame->pcnt_neutral;
  const double second_ref_usage_thresh =
      get_second_ref_usage_thresh(frame_count_so_far);
  int total_frames_to_test = SCENE_CUT_KEY_TEST_INTERVAL;
  int count_for_tolerable_prediction = 3;
  FIRSTPASS_STATS curr_frame;

  if (scenecut_mode == ENABLE_SCENECUT_MODE_1) {
    int num_future_frames = 0;
    curr_frame = *this_frame;
    const FIRSTPASS_STATS *const start_position = twopass->stats_in;
    for (num_future_frames = 0; num_future_frames < SCENE_CUT_KEY_TEST_INTERVAL;
         num_future_frames++)
      if (EOF == input_stats(twopass, &curr_frame)) break;
    reset_fpf_position(twopass, start_position);
    if (num_future_frames < 3) {
      return 0;
    } else {
      total_frames_to_test = 3;
      count_for_tolerable_prediction = 1;
    }
  }

  // Does the frame satisfy the primary criteria of a key frame?
  // See above for an explanation of the test criteria.
  // If so, then examine how well it predicts subsequent frames.
  if (IMPLIES(rc_mode == AOM_Q, frame_count_so_far >= 3) &&
      (this_frame->pcnt_second_ref < second_ref_usage_thresh) &&
      (next_frame->pcnt_second_ref < second_ref_usage_thresh) &&
      ((this_frame->pcnt_inter < VERY_LOW_INTER_THRESH) ||
       ((pcnt_intra > MIN_INTRA_LEVEL) &&
        (pcnt_intra > (INTRA_VS_INTER_THRESH * modified_pcnt_inter)) &&
        ((this_frame->intra_error /
          DOUBLE_DIVIDE_CHECK(this_frame->coded_error)) <
         KF_II_ERR_THRESHOLD) &&
        ((fabs(last_frame->coded_error - this_frame->coded_error) /
              DOUBLE_DIVIDE_CHECK(this_frame->coded_error) >
          ERR_CHANGE_THRESHOLD) ||
         (fabs(last_frame->intra_error - this_frame->intra_error) /
              DOUBLE_DIVIDE_CHECK(this_frame->intra_error) >
          ERR_CHANGE_THRESHOLD) ||
         ((next_frame->intra_error /
           DOUBLE_DIVIDE_CHECK(next_frame->coded_error)) >
          II_IMPROVEMENT_THRESHOLD))))) {
    int i;
    const FIRSTPASS_STATS *start_pos = twopass->stats_in;
    double boost_score = 0.0;
    double old_boost_score = 0.0;
    double decay_accumulator = 1.0;

    // Examine how well the key frame predicts subsequent frames.
    for (i = 0; i < total_frames_to_test; ++i) {
      // Get the next frame details
      FIRSTPASS_STATS local_next_frame;
      if (EOF == input_stats(twopass, &local_next_frame)) break;
      double next_iiratio = (BOOST_FACTOR * local_next_frame.intra_error /
                             DOUBLE_DIVIDE_CHECK(local_next_frame.coded_error));

      if (next_iiratio > KF_II_MAX) next_iiratio = KF_II_MAX;

      // Cumulative effect of decay in prediction quality.
      if (local_next_frame.pcnt_inter > 0.85)
        decay_accumulator *= local_next_frame.pcnt_inter;
      else
        decay_accumulator *= (0.85 + local_next_frame.pcnt_inter) / 2.0;

      // Keep a running total.
      boost_score += (decay_accumulator * next_iiratio);

      // Test various breakout clauses.
      if ((local_next_frame.pcnt_inter < 0.05) || (next_iiratio < 1.5) ||
          (((local_next_frame.pcnt_inter - local_next_frame.pcnt_neutral) <
            0.20) &&
           (next_iiratio < 3.0)) ||
          ((boost_score - old_boost_score) < 3.0) ||
          (local_next_frame.intra_error < 200)) {
        break;
      }

      old_boost_score = boost_score;
    }

    // If there is tolerable prediction for at least the next 3 frames then
    // break out else discard this potential key frame and move on
    if (boost_score > 30.0 && (i > count_for_tolerable_prediction)) {
      is_viable_kf = 1;
    } else {
      is_viable_kf = 0;
    }

    // Reset the file position
    reset_fpf_position(twopass, start_pos);
  }
  return is_viable_kf;
}

#define FRAMES_TO_CHECK_DECAY 8
#define KF_MIN_FRAME_BOOST 80.0
#define KF_MAX_FRAME_BOOST 128.0
#define MIN_KF_BOOST 600  // Minimum boost for non-static KF interval
#define MAX_KF_BOOST 3200
#define MIN_STATIC_KF_BOOST 5400  // Minimum boost for static KF interval

static int detect_app_forced_key(PictureParentControlSet *pcs_ptr) {
  SequenceControlSet *scs_ptr = pcs_ptr->scs_ptr;
  EncodeContext *encode_context_ptr = scs_ptr->encode_context_ptr;
  RATE_CONTROL *const rc = &encode_context_ptr->rc;
  KeyFrameCfg *const kf_cfg = &encode_context_ptr->kf_cfg;
  if (kf_cfg->fwd_kf_enabled) rc->next_is_fwd_key = 1;
  int num_frames_to_app_forced_key = -1;/*is_forced_keyframe_pending(
      cpi->lookahead, cpi->lookahead->max_sz, cpi->compressor_stage);*/
  //if (num_frames_to_app_forced_key != -1) rc->next_is_fwd_key = 0;
  return num_frames_to_app_forced_key;
}

static AOM_INLINE double av1_get_kf_boost_projection_factor(int frame_count) {
  double factor = sqrt((double)frame_count);
  factor = AOMMIN(factor, 10.0);
  factor = AOMMAX(factor, 4.0);
  factor = (75.0 + 14.0 * factor);
  return factor;
}

static int get_projected_kf_boost(RATE_CONTROL *const rc) {
  /*
   * If num_stats_used_for_kf_boost >= frames_to_key, then
   * all stats needed for prior boost calculation are available.
   * Hence projecting the prior boost is not needed in this cases.
   */
  if (rc->num_stats_used_for_kf_boost >= rc->frames_to_key)
    return rc->kf_boost;

  // Get the current tpl factor (number of frames = frames_to_key).
  double tpl_factor = av1_get_kf_boost_projection_factor(rc->frames_to_key);
  // Get the tpl factor when number of frames = num_stats_used_for_kf_boost.
  double tpl_factor_num_stats =
      av1_get_kf_boost_projection_factor(rc->num_stats_used_for_kf_boost);
  int projected_kf_boost =
      (int)rint((tpl_factor * rc->kf_boost) / tpl_factor_num_stats);
  return projected_kf_boost;
}

static int define_kf_interval(PictureParentControlSet *pcs_ptr, FIRSTPASS_STATS *this_frame,
                              double *kf_group_err,
                              int num_frames_to_detect_scenecut) {
  SequenceControlSet *scs_ptr = pcs_ptr->scs_ptr;
  EncodeContext *encode_context_ptr = scs_ptr->encode_context_ptr;
  RATE_CONTROL *const rc = &encode_context_ptr->rc;
  TWO_PASS *const twopass = &scs_ptr->twopass;
  KeyFrameCfg *const kf_cfg = &encode_context_ptr->kf_cfg;
  double recent_loop_decay[FRAMES_TO_CHECK_DECAY];
  FIRSTPASS_STATS last_frame;
  double decay_accumulator = 1.0;
  int j;
  int frames_to_key = 1;
  int frames_since_key = rc->frames_since_key + 1;
  FrameInfo *const frame_info = &encode_context_ptr->frame_info;
  int num_stats_used_for_kf_boost = 1;
  int scenecut_detected = 0;

  int num_frames_to_next_key = detect_app_forced_key(pcs_ptr);

  if (num_frames_to_detect_scenecut == 0) {
    if (num_frames_to_next_key != -1)
      return num_frames_to_next_key;
    else
      return rc->frames_to_key;
  }

  if (num_frames_to_next_key != -1)
    num_frames_to_detect_scenecut =
        AOMMIN(num_frames_to_detect_scenecut, num_frames_to_next_key);

  // Initialize the decay rates for the recent frames to check
  for (j = 0; j < FRAMES_TO_CHECK_DECAY; ++j) recent_loop_decay[j] = 1.0;

  while (twopass->stats_in < twopass->stats_buf_ctx->stats_in_end &&
         frames_to_key < num_frames_to_detect_scenecut) {
    int i = 0;
    // Accumulate total number of stats available till next key frame
    num_stats_used_for_kf_boost++;

    // Accumulate kf group error.
    if (kf_group_err != NULL)
      *kf_group_err +=
          calculate_modified_err(frame_info, twopass, &(encode_context_ptr->two_pass_cfg), this_frame);

    // Load the next frame's stats.
    last_frame = *this_frame;
    input_stats(twopass, this_frame);

    // Provided that we are not at the end of the file...
    if ((rc->enable_scenecut_detection > 0) && kf_cfg->auto_key &&
        twopass->stats_in < twopass->stats_buf_ctx->stats_in_end) {
      double loop_decay_rate;

      // Check for a scene cut.
      if (frames_since_key >= kf_cfg->key_freq_min &&
          test_candidate_kf(twopass, &last_frame, this_frame, twopass->stats_in,
                            frames_since_key, encode_context_ptr->rc_cfg.mode,
                            rc->enable_scenecut_detection)) {
        scenecut_detected = 1;
        break;
      }

      // How fast is the prediction quality decaying?
      loop_decay_rate =
          get_prediction_decay_rate(frame_info, twopass->stats_in);

      // We want to know something about the recent past... rather than
      // as used elsewhere where we are concerned with decay in prediction
      // quality since the last GF or KF.
      recent_loop_decay[i % FRAMES_TO_CHECK_DECAY] = loop_decay_rate;
      decay_accumulator = 1.0;
      for (j = 0; j < FRAMES_TO_CHECK_DECAY; ++j)
        decay_accumulator *= recent_loop_decay[j];

      // Special check for transition or high motion followed by a
      // static scene.
      if (frames_since_key >= kf_cfg->key_freq_min &&
          detect_transition_to_still(twopass, rc->min_gf_interval, i,
                                     kf_cfg->key_freq_max - i, loop_decay_rate,
                                     decay_accumulator)) {
        scenecut_detected = 1;
        // In the case of transition followed by a static scene, the key frame
        // could be a good predictor for the following frames, therefore we
        // do not use an arf.
        rc->use_arf_in_this_kf_group = 0;
        break;
      }

      // Step on to the next frame.
      ++frames_to_key;
      ++frames_since_key;

      // If we don't have a real key frame within the next two
      // key_freq_max intervals then break out of the loop.
      if (frames_to_key >= 2 * kf_cfg->key_freq_max) break;
    } else {
      ++frames_to_key;
      ++frames_since_key;
    }
    ++i;
  }

  if (kf_group_err != NULL)
    rc->num_stats_used_for_kf_boost = num_stats_used_for_kf_boost;

  if (scs_ptr->lap_enabled && !scenecut_detected)
    frames_to_key = num_frames_to_next_key;

  if (kf_cfg->fwd_kf_enabled && scenecut_detected) rc->next_is_fwd_key = 0;

  return frames_to_key;
}

static int64_t get_kf_group_bits(PictureParentControlSet *pcs_ptr, double kf_group_err/*,
                                 double kf_group_avg_error*/) {
  SequenceControlSet *scs_ptr = pcs_ptr->scs_ptr;
  //EncodeContext *encode_context_ptr = scs_ptr->encode_context_ptr;
  //RATE_CONTROL *const rc = &encode_context_ptr->rc;
  TWO_PASS *const twopass = &scs_ptr->twopass;
  int64_t kf_group_bits;

  {
    kf_group_bits = (int64_t)(twopass->bits_left *
                              (kf_group_err / twopass->modified_error_left));
  }

  return kf_group_bits;
}

static int calc_avg_stats(PictureParentControlSet *pcs_ptr, FIRSTPASS_STATS *avg_frame_stat) {
  SequenceControlSet *scs_ptr = pcs_ptr->scs_ptr;
  EncodeContext *encode_context_ptr = scs_ptr->encode_context_ptr;
  RATE_CONTROL *const rc = &encode_context_ptr->rc;
  TWO_PASS *const twopass = &scs_ptr->twopass;
  FIRSTPASS_STATS cur_frame;
  av1_zero(cur_frame);
  int num_frames = 0;
  // Accumulate total stat using available number of stats.
  for (num_frames = 0; num_frames < (rc->frames_to_key - 1); ++num_frames) {
    if (EOF == input_stats(twopass, &cur_frame)) break;
    svt_av1_accumulate_stats(avg_frame_stat, &cur_frame);
  }

  if (num_frames < 2) {
    return num_frames;
  }
  // Average the total stat
  avg_frame_stat->weight = avg_frame_stat->weight / num_frames;
  avg_frame_stat->intra_error = avg_frame_stat->intra_error / num_frames;
  avg_frame_stat->frame_avg_wavelet_energy =
      avg_frame_stat->frame_avg_wavelet_energy / num_frames;
  avg_frame_stat->coded_error = avg_frame_stat->coded_error / num_frames;
  avg_frame_stat->sr_coded_error = avg_frame_stat->sr_coded_error / num_frames;
  avg_frame_stat->pcnt_inter = avg_frame_stat->pcnt_inter / num_frames;
  avg_frame_stat->pcnt_motion = avg_frame_stat->pcnt_motion / num_frames;
  avg_frame_stat->pcnt_second_ref =
      avg_frame_stat->pcnt_second_ref / num_frames;
  avg_frame_stat->pcnt_neutral = avg_frame_stat->pcnt_neutral / num_frames;
  avg_frame_stat->intra_skip_pct = avg_frame_stat->intra_skip_pct / num_frames;
  avg_frame_stat->inactive_zone_rows =
      avg_frame_stat->inactive_zone_rows / num_frames;
  avg_frame_stat->inactive_zone_cols =
      avg_frame_stat->inactive_zone_cols / num_frames;
  avg_frame_stat->MVr = avg_frame_stat->MVr / num_frames;
  avg_frame_stat->mvr_abs = avg_frame_stat->mvr_abs / num_frames;
  avg_frame_stat->MVc = avg_frame_stat->MVc / num_frames;
  avg_frame_stat->mvc_abs = avg_frame_stat->mvc_abs / num_frames;
  avg_frame_stat->MVrv = avg_frame_stat->MVrv / num_frames;
  avg_frame_stat->MVcv = avg_frame_stat->MVcv / num_frames;
  avg_frame_stat->mv_in_out_count =
      avg_frame_stat->mv_in_out_count / num_frames;
  avg_frame_stat->new_mv_count = avg_frame_stat->new_mv_count / num_frames;
  avg_frame_stat->count = avg_frame_stat->count / num_frames;
  avg_frame_stat->duration = avg_frame_stat->duration / num_frames;

  return num_frames;
}

static double get_kf_boost_score(PictureParentControlSet *pcs_ptr, double kf_raw_err,
                                 double *zero_motion_accumulator,
                                 double *sr_accumulator, int use_avg_stat) {
  SequenceControlSet *scs_ptr = pcs_ptr->scs_ptr;
  EncodeContext *encode_context_ptr = scs_ptr->encode_context_ptr;
  RATE_CONTROL *const rc = &encode_context_ptr->rc;
  TWO_PASS *const twopass = &scs_ptr->twopass;
  FrameInfo *frame_info = &encode_context_ptr->frame_info;
  const RateControlCfg *const rc_cfg = &encode_context_ptr->rc_cfg;
  FIRSTPASS_STATS frame_stat;
  av1_zero(frame_stat);
  int i = 0, num_stat_used = 0;
  double boost_score = 0.0;
  const double kf_max_boost =
      rc_cfg->mode == AOM_Q
          ? AOMMIN(AOMMAX(rc->frames_to_key * 2.0, KF_MIN_FRAME_BOOST),
                   KF_MAX_FRAME_BOOST)
          : KF_MAX_FRAME_BOOST;

  // Calculate the average using available number of stats.
  if (use_avg_stat) num_stat_used = calc_avg_stats(pcs_ptr, &frame_stat);

  for (i = num_stat_used; i < (rc->frames_to_key - 1); ++i) {
    if (!use_avg_stat && EOF == input_stats(twopass, &frame_stat)) break;

    // Monitor for static sections.
    // For the first frame in kf group, the second ref indicator is invalid.
    if (i > 0) {
      *zero_motion_accumulator =
          AOMMIN(*zero_motion_accumulator,
                 get_zero_motion_factor(frame_info, &frame_stat));
    } else {
      *zero_motion_accumulator = frame_stat.pcnt_inter - frame_stat.pcnt_motion;
    }

    // Not all frames in the group are necessarily used in calculating boost.
    if ((*sr_accumulator < (kf_raw_err * 1.50)) &&
        (i <= rc->max_gf_interval * 2)) {
      double frame_boost;
      double zm_factor;

      // Factor 0.75-1.25 based on how much of frame is static.
      zm_factor = (0.75 + (*zero_motion_accumulator / 2.0));

      if (i < 2) *sr_accumulator = 0.0;
      frame_boost = calc_kf_frame_boost(rc, frame_info, &frame_stat,
                                        sr_accumulator, kf_max_boost);
      boost_score += frame_boost * zm_factor;
    }
  }
  return boost_score;
}
// This function imposes the key frames based on the intra refresh period
//anaghdin only works for key frames, and scene change is not detected for now
static void find_next_key_frame(PictureParentControlSet *pcs_ptr, FIRSTPASS_STATS *this_frame) {

    SequenceControlSet *scs_ptr = pcs_ptr->scs_ptr;
    EncodeContext *encode_context_ptr = scs_ptr->encode_context_ptr;
    RATE_CONTROL *const rc = &encode_context_ptr->rc;
    TWO_PASS *const twopass = &scs_ptr->twopass;
    GF_GROUP *const gf_group = &encode_context_ptr->gf_group;
    //CurrentFrame *const current_frame = &pcs_ptr->av1_cm->current_frame;
    FrameInfo *frame_info = &encode_context_ptr->frame_info;
    const KeyFrameCfg *const kf_cfg = &encode_context_ptr->kf_cfg;
    const GFConfig *const gf_cfg = &encode_context_ptr->gf_cfg;
    const RateControlCfg *const rc_cfg = &encode_context_ptr->rc_cfg;
    const FIRSTPASS_STATS first_frame = *this_frame;
    FIRSTPASS_STATS next_frame;
    av1_zero(next_frame);

    rc->frames_since_key = 0;
    // Use arfs if possible.
    rc->use_arf_in_this_kf_group = is_altref_enabled(
        gf_cfg->lag_in_frames, gf_cfg->enable_auto_arf);

    // Reset the GF group data structures.
    av1_zero(encode_context_ptr->gf_group);
    pcs_ptr->gf_group_index = 0;

    // Clear the alt ref active flag and last group multi arf flags as they
    // can never be set for a key frame.
    rc->source_alt_ref_active = 0;

    // KF is always a GF so clear frames till next gf counter.
    rc->frames_till_gf_update_due = 0;
    rc->frames_to_key = 1;

    const FIRSTPASS_STATS *const start_position = twopass->stats_in;
    int kf_bits = 0;
    double zero_motion_accumulator = 1.0;
    double boost_score = 0.0;
    double kf_raw_err;
    double kf_mod_err;
    double kf_group_err = 0.0;
    double sr_accumulator = 0.0;
    //double kf_group_avg_error = 0.0;
    int frames_to_key;
    // Is this a forced key frame by interval.
    rc->this_key_frame_forced = rc->next_key_frame_forced;

    twopass->kf_group_bits = 0;        // Total bits available to kf group
    twopass->kf_group_error_left = 0;  // Group modified error score.

    kf_raw_err = this_frame->intra_error;
    kf_mod_err = calculate_modified_err(frame_info, twopass, &(encode_context_ptr->two_pass_cfg), this_frame);

    frames_to_key =
        define_kf_interval(pcs_ptr, this_frame, &kf_group_err, kf_cfg->key_freq_max);

    if (frames_to_key != -1)
        rc->frames_to_key = AOMMIN(kf_cfg->key_freq_max, frames_to_key);
    else
        rc->frames_to_key = kf_cfg->key_freq_max;

    //if (scs_ptr->lap_enabled) correct_frames_to_key(cpi);

    // If there is a max kf interval set by the user we must obey it.
    // We already breakout of the loop above at 2x max.
    // This code centers the extra kf if the actual natural interval
    // is between 1x and 2x.
    if (kf_cfg->auto_key && rc->frames_to_key > kf_cfg->key_freq_max) {
        FIRSTPASS_STATS tmp_frame = first_frame;

        rc->frames_to_key /= 2;

        // Reset to the start of the group.
        reset_fpf_position(twopass, start_position);

        kf_group_err = 0.0;

        // Rescan to get the correct error data for the forced kf group.
        for (int i = 0; i < rc->frames_to_key; ++i) {
            kf_group_err +=
                calculate_modified_err(frame_info, twopass, &(encode_context_ptr->two_pass_cfg), &tmp_frame);
            if (EOF == input_stats(twopass, &tmp_frame)) break;
        }
        rc->next_key_frame_forced = 1;
    }
    else if ((twopass->stats_in == twopass->stats_buf_ctx->stats_in_end/* &&
             is_stat_consumption_stage_twopass(cpi)*/) ||
        rc->frames_to_key >= kf_cfg->key_freq_max) {
        rc->next_key_frame_forced = 0;
    }
    else {
        rc->next_key_frame_forced = 0;
    }

    if (kf_cfg->fwd_kf_enabled) rc->next_is_fwd_key |= rc->next_key_frame_forced;

    // Special case for the last key frame of the file.
    if (twopass->stats_in >= twopass->stats_buf_ctx->stats_in_end) {
        // Accumulate kf group error.
        kf_group_err +=
            calculate_modified_err(frame_info, twopass, &(encode_context_ptr->two_pass_cfg), this_frame);
    }

    // Calculate the number of bits that should be assigned to the kf group.
    if ((twopass->bits_left > 0 && twopass->modified_error_left > 0.0) ||
        (scs_ptr->lap_enabled && rc_cfg->mode != AOM_Q)) {
        // Maximum number of bits for a single normal frame (not key frame).
        const int max_bits = frame_max_bits(rc, encode_context_ptr);

        // Maximum number of bits allocated to the key frame group.
        int64_t max_grp_bits;

        // Default allocation based on bits left and relative
        // complexity of the section.
        twopass->kf_group_bits =
            get_kf_group_bits(pcs_ptr, kf_group_err/*, kf_group_avg_error*/);
        // Clip based on maximum per frame rate defined by the user.
        max_grp_bits = (int64_t)max_bits * (int64_t)rc->frames_to_key;
        if (twopass->kf_group_bits > max_grp_bits)
            twopass->kf_group_bits = max_grp_bits;
    }
    else {
        twopass->kf_group_bits = 0;
    }
    twopass->kf_group_bits = AOMMAX(0, twopass->kf_group_bits);

    // Reset the first pass file position.
    reset_fpf_position(twopass, start_position);

    // Scan through the kf group collating various stats used to determine
    // how many bits to spend on it.
    boost_score = get_kf_boost_score(pcs_ptr, kf_raw_err, &zero_motion_accumulator,
        &sr_accumulator, 0);
    reset_fpf_position(twopass, start_position);
    // Store the zero motion percentage
    twopass->kf_zeromotion_pct = (int)(zero_motion_accumulator * 100.0);

    // Calculate a section intra ratio used in setting max loop filter.
    twopass->section_intra_rating = calculate_section_intra_ratio(
        start_position, twopass->stats_buf_ctx->stats_in_end, rc->frames_to_key);

    rc->kf_boost = (int)boost_score;

    if (scs_ptr->lap_enabled) {
        if (rc_cfg->mode == AOM_Q) {
            rc->kf_boost = get_projected_kf_boost(rc);
        }
        else {
            // TODO(any): Explore using average frame stats for AOM_Q as well.
            boost_score = get_kf_boost_score(
                pcs_ptr, kf_raw_err, &zero_motion_accumulator, &sr_accumulator, 1);
            reset_fpf_position(twopass, start_position);
            rc->kf_boost += (int)boost_score;
        }
    }

    // Special case for static / slide show content but don't apply
    // if the kf group is very short.
    if ((zero_motion_accumulator > STATIC_KF_GROUP_FLOAT_THRESH) &&
        (rc->frames_to_key > 8)) {
        rc->kf_boost = AOMMAX(rc->kf_boost, MIN_STATIC_KF_BOOST);
    }
    else {
        // Apply various clamps for min and max boost
        rc->kf_boost = AOMMAX(rc->kf_boost, (rc->frames_to_key * 3));
        rc->kf_boost = AOMMAX(rc->kf_boost, MIN_KF_BOOST);
#ifdef STRICT_RC
        rc->kf_boost = AOMMIN(rc->kf_boost, MAX_KF_BOOST);
#endif
    }

    // Work out how many bits to allocate for the key frame itself.
    kf_bits = calculate_boost_bits((rc->frames_to_key - 1), rc->kf_boost,
        twopass->kf_group_bits);

    twopass->kf_group_bits -= kf_bits;

    // Save the bits to spend on the key frame.
    gf_group->bit_allocation[0] = kf_bits;
    gf_group->update_type[0] = KF_UPDATE;

    // Note the total error score of the kf group minus the key frame itself.
    twopass->kf_group_error_left = (int)(kf_group_err - kf_mod_err);

    // Adjust the count of total modified error left.
    // The count of bits left is adjusted elsewhere based on real coded frame
    // sizes.
    twopass->modified_error_left -= kf_group_err;
}

#define ARF_STATS_OUTPUT 0
#if ARF_STATS_OUTPUT
unsigned int arf_count = 0;
#endif
#define DEFAULT_GRP_WEIGHT 1.0

static int get_section_target_bandwidth(PictureParentControlSet *pcs_ptr) {
  SequenceControlSet *scs_ptr = pcs_ptr->scs_ptr;
  EncodeContext *encode_context_ptr = scs_ptr->encode_context_ptr;
  TWO_PASS *const twopass = &scs_ptr->twopass;
  RATE_CONTROL *const rc = &encode_context_ptr->rc;
  int section_target_bandwidth;
  const int frames_left = (int)(twopass->stats_buf_ctx->total_stats->count -
                                pcs_ptr->picture_number);
  if (scs_ptr->lap_enabled)
    section_target_bandwidth = (int)rc->avg_frame_bandwidth;
  else
    section_target_bandwidth = (int)(twopass->bits_left / frames_left);
  return section_target_bandwidth;
}

static void process_first_pass_stats(PictureParentControlSet *pcs_ptr,
                                     FIRSTPASS_STATS *this_frame) {
  //CurrentFrame *const current_frame = &cm->current_frame;
  SequenceControlSet *scs_ptr = pcs_ptr->scs_ptr;
  EncodeContext *encode_context_ptr = scs_ptr->encode_context_ptr;
  TWO_PASS *const twopass = &scs_ptr->twopass;
  RATE_CONTROL *const rc = &encode_context_ptr->rc;
  const RateControlCfg *const rc_cfg = &encode_context_ptr->rc_cfg;
  const uint32_t mb_cols = (scs_ptr->seq_header.max_frame_width  + 16 - 1) / 16;
  const uint32_t mb_rows = (scs_ptr->seq_header.max_frame_height + 16 - 1) / 16;

  if (rc_cfg->mode != AOM_Q && /*current_frame->frame_number*/pcs_ptr->picture_number == 0 &&
      twopass->stats_buf_ctx->total_stats &&
      twopass->stats_buf_ctx->total_left_stats) {
    if (scs_ptr->lap_enabled) {
      /*
       * Accumulate total_stats using available limited number of stats,
       * and assign it to total_left_stats.
       */
      *twopass->stats_buf_ctx->total_left_stats =
          *twopass->stats_buf_ctx->total_stats;
    }
    // Special case code for first frame.
    const int section_target_bandwidth = get_section_target_bandwidth(pcs_ptr);
    const double section_length =
        twopass->stats_buf_ctx->total_left_stats->count;
    const double section_error =
        twopass->stats_buf_ctx->total_left_stats->coded_error / section_length;
    const double section_intra_skip =
        twopass->stats_buf_ctx->total_left_stats->intra_skip_pct /
        section_length;
    const double section_inactive_zone =
        (twopass->stats_buf_ctx->total_left_stats->inactive_zone_rows * 2) /
        ((double)/*cm->mi_params.*/mb_rows * section_length);
    const int tmp_q = get_twopass_worst_quality(
        pcs_ptr, section_error, section_intra_skip + section_inactive_zone,
        section_target_bandwidth, DEFAULT_GRP_WEIGHT);

    rc->active_worst_quality = tmp_q;
    rc->ni_av_qi = tmp_q;
    rc->last_q[INTER_FRAME] = tmp_q;
    rc->avg_q = svt_av1_convert_qindex_to_q(tmp_q, scs_ptr->encoder_bit_depth);
    rc->avg_frame_qindex[INTER_FRAME] = tmp_q;
    rc->last_q[KEY_FRAME] = (tmp_q + rc_cfg->best_allowed_q) / 2;
    rc->avg_frame_qindex[KEY_FRAME] = rc->last_q[KEY_FRAME];
  }

  int err = 0;
  if (scs_ptr->lap_enabled) {
    //err = input_stats_lap(twopass, this_frame);
  } else {
    err = input_stats(twopass, this_frame);
  }
  if (err == EOF) return;

  {
    const int num_mbs = mb_cols * mb_rows;
                        //(cpi->oxcf.resize_cfg.resize_mode != RESIZE_NONE)
                        //    ? cpi->initial_mbs
                        //    : cm->mi_params.MBs;
    // The multiplication by 256 reverses a scaling factor of (>> 8)
    // applied when combining MB error values for the frame.
    twopass->mb_av_energy = log1p(this_frame->intra_error / num_mbs);
    twopass->frame_avg_haar_energy =
        log1p(this_frame->frame_avg_wavelet_energy / num_mbs);
  }

  // Update the total stats remaining structure.
  if (twopass->stats_buf_ctx->total_left_stats)
    subtract_stats(twopass->stats_buf_ctx->total_left_stats, this_frame);

  // Set the frame content type flag.
  if (this_frame->intra_skip_pct >= FC_ANIMATION_THRESH)
    twopass->fr_content_type = FC_GRAPHICS_ANIMATION;
  else
    twopass->fr_content_type = FC_NORMAL;
}

static void setup_target_rate(PictureParentControlSet *pcs_ptr) {
    SequenceControlSet *scs_ptr = pcs_ptr->scs_ptr;
  EncodeContext *encode_context_ptr = scs_ptr->encode_context_ptr;
  RATE_CONTROL *const rc = &encode_context_ptr->rc;
  GF_GROUP *const gf_group = &encode_context_ptr->gf_group;

  int target_rate = gf_group->bit_allocation[pcs_ptr->gf_group_index];

  rc->base_frame_target = target_rate;
}

void svt_av1_get_second_pass_params(PictureParentControlSet *pcs_ptr) {
    SequenceControlSet *scs_ptr = pcs_ptr->scs_ptr;
    EncodeContext *encode_context_ptr = scs_ptr->encode_context_ptr;
    RATE_CONTROL *const rc = &encode_context_ptr->rc;
    TWO_PASS *const twopass = &scs_ptr->twopass;
    GF_GROUP *const gf_group = &encode_context_ptr->gf_group;
    CurrentFrame *const current_frame = &pcs_ptr->av1_cm->current_frame;
    current_frame->frame_number = (int)pcs_ptr->picture_number;

    EncodeFrameParams temp_frame_params, *frame_params = &temp_frame_params;
    pcs_ptr->gf_group_index = gf_group->index;
  if (/*is_stat_consumption_stage(cpi) &&*/ !twopass->stats_in) return;
#ifdef ARCH_X86_64
  aom_clear_system_state();
#endif
  if (encode_context_ptr->rc_cfg.mode == AOM_Q)
    rc->active_worst_quality = encode_context_ptr->rc_cfg.cq_level;
  FIRSTPASS_STATS this_frame;
  av1_zero(this_frame);
  // call above fn
  if (1/*is_stat_consumption_stage(cpi)*/) {
    process_first_pass_stats(pcs_ptr, &this_frame);
  } else {
    rc->active_worst_quality = encode_context_ptr->rc_cfg.cq_level;
  }

  // Keyframe and section processing.
    if (rc->frames_to_key <= 0 || (frame_is_intra_only(pcs_ptr) && pcs_ptr->idr_flag)/*(frame_flags & FRAMEFLAGS_KEY)*/) {
    assert(rc->frames_to_key >= -1);
    FIRSTPASS_STATS this_frame_copy;
    this_frame_copy = this_frame;
    frame_params->frame_type = KEY_FRAME;
    // Define next KF group and assign bits to it.
    find_next_key_frame(pcs_ptr, &this_frame);
    this_frame = this_frame_copy;
  } else {
    frame_params->frame_type = INTER_FRAME;
    const int altref_enabled = is_altref_enabled(encode_context_ptr->gf_cfg.lag_in_frames,
                                                 encode_context_ptr->gf_cfg.enable_auto_arf);
    const int sframe_dist = encode_context_ptr->kf_cfg.sframe_dist;
    const int sframe_mode = encode_context_ptr->kf_cfg.sframe_mode;
    const int update_type = gf_group->update_type[pcs_ptr->gf_group_index];
    if (sframe_dist != 0) {
      if (altref_enabled) {
        if (sframe_mode == 1) {
          // sframe_mode == 1: insert sframe if it matches altref frame.
          if (current_frame->frame_number % sframe_dist == 0 &&
              current_frame->frame_number != 0 && update_type == ARF_UPDATE) {
            frame_params->frame_type = S_FRAME;
          }
        } else {
          // sframe_mode != 1: if sframe will be inserted at the next available
          // altref frame
          if (current_frame->frame_number % sframe_dist == 0 &&
              current_frame->frame_number != 0) {
            rc->sframe_due = 1;
          }
          if (rc->sframe_due && update_type == ARF_UPDATE) {
            frame_params->frame_type = S_FRAME;
            rc->sframe_due = 0;
          }
        }
      } else {
        if (current_frame->frame_number % sframe_dist == 0 &&
            current_frame->frame_number != 0) {
          frame_params->frame_type = S_FRAME;
        }
      }
    }
  }

  // Define a new GF/ARF group. (Should always enter here for key frames).
  if (rc->frames_till_gf_update_due == 0) {
    assert(current_frame->frame_number == 0 ||
        pcs_ptr->gf_group_index == gf_group->size);
    const FIRSTPASS_STATS *const start_position = twopass->stats_in;

    if (scs_ptr->lap_enabled && rc->enable_scenecut_detection) {
      int num_frames_to_detect_scenecut, frames_to_key;
      num_frames_to_detect_scenecut = MAX_GF_LENGTH_LAP + 1;
      frames_to_key = define_kf_interval(pcs_ptr, &this_frame, NULL,
                                         num_frames_to_detect_scenecut);
      if (frames_to_key != -1)
        rc->frames_to_key = AOMMIN(rc->frames_to_key, frames_to_key);
    }

    reset_fpf_position(twopass, start_position);

    int max_gop_length =
        (encode_context_ptr->gf_cfg.lag_in_frames >= 32 &&
         1/*is_stat_consumption_stage_twopass(cpi)*/)
            ? AOMMIN(MAX_GF_INTERVAL,
                     encode_context_ptr->gf_cfg.lag_in_frames - 7/*oxcf->arnr_max_frames*/ / 2)
            : MAX_GF_LENGTH_LAP;
    if (rc->intervals_till_gf_calculate_due == 0)
      impose_gf_length(pcs_ptr, MAX_NUM_GF_INTERVALS);

    define_gf_group(pcs_ptr, &this_frame, frame_params, max_gop_length, 1);
    rc->frames_till_gf_update_due = rc->baseline_gf_interval;
    assert(pcs_ptr->gf_group_index == 0);
    // This is added for the first frame in minigop when it is not KEY_FRAME
    if (pcs_ptr->frm_hdr.frame_type != KEY_FRAME) {
        gf_group->index++;
        pcs_ptr->gf_group_index = gf_group->index;
    }
#if ARF_STATS_OUTPUT
    {
      FILE *fpfile;
      fpfile = fopen("arf.stt", "a");
      ++arf_count;
      fprintf(fpfile, "%10d %10d %10d %10d %10d\n",
              cpi->common.current_frame.frame_number,
              rc->frames_till_gf_update_due, rc->kf_boost, arf_count,
              rc->gfu_boost);

      fclose(fpfile);
    }
#endif
  }
  assert(pcs_ptr->gf_group_index < gf_group->size);

  setup_target_rate(pcs_ptr);
}

// from aom ratectrl.c
static void av1_rc_set_gf_interval_range(SequenceControlSet *scs_ptr,
                                  RATE_CONTROL *const rc) {
  EncodeContext *encode_context_ptr = scs_ptr->encode_context_ptr;
  const uint32_t width  = scs_ptr->seq_header.max_frame_width;
  const uint32_t height = scs_ptr->seq_header.max_frame_height;

  // Special case code for 1 pass fixed Q mode tests
  if (0/*(has_no_stats_stage(cpi))*/ && (encode_context_ptr->rc_cfg.mode == AOM_Q)) {
    rc->max_gf_interval = FIXED_GF_INTERVAL;
    rc->min_gf_interval = FIXED_GF_INTERVAL;
    rc->static_scene_max_gf_interval = FIXED_GF_INTERVAL;
  } else {
    // Set Maximum gf/arf interval
    rc->max_gf_interval = encode_context_ptr->gf_cfg.max_gf_interval;
    rc->min_gf_interval = encode_context_ptr->gf_cfg.min_gf_interval;
    if (rc->min_gf_interval == 0)
      rc->min_gf_interval = svt_av1_rc_get_default_min_gf_interval(
          width, height, scs_ptr->double_frame_rate);
          //oxcf->frm_dim_cfg.width, oxcf->frm_dim_cfg.height, cpi->framerate);
    if (rc->max_gf_interval == 0)
      rc->max_gf_interval = svt_av1_rc_get_default_max_gf_interval(
          scs_ptr->double_frame_rate/*cpi->framerate*/, rc->min_gf_interval);
    /*
     * Extended max interval for genuinely static scenes like slide shows.
     * The no.of.stats available in the case of LAP is limited,
     * hence setting to max_gf_interval.
     */
    if (scs_ptr->lap_enabled)
      rc->static_scene_max_gf_interval = rc->max_gf_interval + 1;
    else
      rc->static_scene_max_gf_interval = MAX_STATIC_GF_GROUP_LENGTH;

    if (rc->max_gf_interval > rc->static_scene_max_gf_interval)
      rc->max_gf_interval = rc->static_scene_max_gf_interval;

    // Clamp min to max
    rc->min_gf_interval = AOMMIN(rc->min_gf_interval, rc->max_gf_interval);
  }
}

// Max rate target for 1080P and below encodes under normal circumstances
// (1920 * 1080 / (16 * 16)) * MAX_MB_RATE bits per MB
#define MAX_MB_RATE 250
#define MAXRATE_1080P 2025000
#define FRAME_OVERHEAD_BITS 200
static void av1_rc_update_framerate(SequenceControlSet *scs_ptr/*, int width, int height*/) {
  EncodeContext *encode_context_ptr = scs_ptr->encode_context_ptr;
  RATE_CONTROL *const rc = &encode_context_ptr->rc;
  FrameInfo *frame_info = &encode_context_ptr->frame_info;
  int vbr_max_bits;
  const int MBs = frame_info->num_mbs;//av1_get_MBs(width, height);

  rc->avg_frame_bandwidth = (int)(scs_ptr->static_config.target_bit_rate/*oxcf->target_bandwidth*/ / scs_ptr->double_frame_rate);
  rc->min_frame_bandwidth =
      (int)(rc->avg_frame_bandwidth * encode_context_ptr->two_pass_cfg.vbrmin_section / 100);

  rc->min_frame_bandwidth =
      AOMMAX(rc->min_frame_bandwidth, FRAME_OVERHEAD_BITS);

  // A maximum bitrate for a frame is defined.
  // The baseline for this aligns with HW implementations that
  // can support decode of 1080P content up to a bitrate of MAX_MB_RATE bits
  // per 16x16 MB (averaged over a frame). However this limit is extended if
  // a very high rate is given on the command line or the the rate cannnot
  // be acheived because of a user specificed max q (e.g. when the user
  // specifies lossless encode.
  vbr_max_bits = (int)(((int64_t)rc->avg_frame_bandwidth *
                        encode_context_ptr->two_pass_cfg.vbrmax_section) /
                       100);
  rc->max_frame_bandwidth =
      AOMMAX(AOMMAX((MBs * MAX_MB_RATE), MAXRATE_1080P), vbr_max_bits);

  av1_rc_set_gf_interval_range(scs_ptr, rc);
}

// from aom encoder.c
static void av1_new_framerate(SequenceControlSet *scs_ptr, double framerate) {
  //cpi->framerate = framerate < 0.1 ? 30 : framerate;
  scs_ptr->double_frame_rate = framerate < 0.1 ? 30 : framerate;
  av1_rc_update_framerate(scs_ptr/*, scs_ptr->seq_header.max_frame_width, scs_ptr->seq_header.max_frame_height*/);
}

void svt_av1_init_second_pass(SequenceControlSet *scs_ptr) {
  TWO_PASS *const twopass = &scs_ptr->twopass;
  EncodeContext *encode_context_ptr = scs_ptr->encode_context_ptr;
  FrameInfo *frame_info = &encode_context_ptr->frame_info;

  double frame_rate;
  FIRSTPASS_STATS *stats;

  if (!twopass->stats_buf_ctx->stats_in_end) return;

  {
      const int is_vbr = scs_ptr->static_config.rate_control_mode == 1;
      frame_info->frame_width  = scs_ptr->seq_header.max_frame_width;
      frame_info->frame_height = scs_ptr->seq_header.max_frame_height;
      frame_info->mb_cols = (scs_ptr->seq_header.max_frame_width  + 16 - 1) / 16;
      frame_info->mb_rows = (scs_ptr->seq_header.max_frame_height + 16 - 1) / 16;
      frame_info->num_mbs = frame_info->mb_cols * frame_info->mb_rows;
      frame_info->bit_depth = scs_ptr->static_config.encoder_bit_depth;
      // input config  from options
      encode_context_ptr->two_pass_cfg.vbrmin_section = scs_ptr->static_config.vbr_min_section_pct;
      encode_context_ptr->two_pass_cfg.vbrmax_section = scs_ptr->static_config.vbr_max_section_pct;
      encode_context_ptr->two_pass_cfg.vbrbias        = scs_ptr->static_config.vbr_bias_pct;
      encode_context_ptr->rc_cfg.mode = scs_ptr->static_config.rate_control_mode == 1 ? AOM_VBR : AOM_Q;
      encode_context_ptr->rc_cfg.best_allowed_q  = (int32_t)quantizer_to_qindex[scs_ptr->static_config.min_qp_allowed];
      encode_context_ptr->rc_cfg.worst_allowed_q = (int32_t)quantizer_to_qindex[scs_ptr->static_config.max_qp_allowed];
      encode_context_ptr->rc_cfg.over_shoot_pct  = scs_ptr->static_config.over_shoot_pct;
      encode_context_ptr->rc_cfg.under_shoot_pct = scs_ptr->static_config.under_shoot_pct;
      encode_context_ptr->rc_cfg.cq_level = quantizer_to_qindex[scs_ptr->static_config.qp];
      encode_context_ptr->rc_cfg.maximum_buffer_size_ms   = is_vbr ? 240000 : 6000;//cfg->rc_buf_sz;
      encode_context_ptr->rc_cfg.starting_buffer_level_ms = is_vbr ? 60000  : 4000;//cfg->rc_buf_initial_sz;
      encode_context_ptr->rc_cfg.optimal_buffer_level_ms  = is_vbr ? 60000  : 5000;//cfg->rc_buf_optimal_sz;
      encode_context_ptr->gf_cfg.lag_in_frames = 25;//hack scs_ptr->static_config.look_ahead_distance + 1;
      encode_context_ptr->gf_cfg.gf_min_pyr_height = scs_ptr->static_config.hierarchical_levels;
      encode_context_ptr->gf_cfg.gf_max_pyr_height = scs_ptr->static_config.hierarchical_levels;
      encode_context_ptr->gf_cfg.min_gf_interval   = 1 << scs_ptr->static_config.hierarchical_levels;
      encode_context_ptr->gf_cfg.max_gf_interval   = 1 << scs_ptr->static_config.hierarchical_levels;
      encode_context_ptr->gf_cfg.enable_auto_arf   = 1;
      encode_context_ptr->kf_cfg.sframe_dist   = 0; // not supported yet
      encode_context_ptr->kf_cfg.sframe_mode   = 0; // not supported yet
      encode_context_ptr->kf_cfg.auto_key      = 0;
      encode_context_ptr->kf_cfg.key_freq_max  = scs_ptr->intra_period_length + 1;
  }

  stats = twopass->stats_buf_ctx->total_stats;

  *stats = *twopass->stats_buf_ctx->stats_in_end;
  *twopass->stats_buf_ctx->total_left_stats = *stats;

  frame_rate = 10000000.0 * stats->count / stats->duration;
  // Each frame can have a different duration, as the frame rate in the source
  // isn't guaranteed to be constant. The frame rate prior to the first frame
  // encoded in the second pass is a guess. However, the sum duration is not.
  // It is calculated based on the actual durations of all frames from the
  // first pass.
  av1_new_framerate(scs_ptr, frame_rate);
  twopass->bits_left =
      (int64_t)(stats->duration * (int64_t)scs_ptr->static_config.target_bit_rate / 10000000.0);

  // This variable monitors how far behind the second ref update is lagging.
  twopass->sr_update_lag = 1;

  // Scan the first pass file and calculate a modified total error based upon
  // the bias/power function used to allocate bits.
  {
    const double avg_error =
        stats->coded_error / DOUBLE_DIVIDE_CHECK(stats->count);
    const FIRSTPASS_STATS *s = twopass->stats_in;
    double modified_error_total = 0.0;
    twopass->modified_error_min =
        (avg_error * encode_context_ptr->two_pass_cfg.vbrmin_section) / 100;
    twopass->modified_error_max =
        (avg_error * encode_context_ptr->two_pass_cfg.vbrmax_section) / 100;
    while (s < twopass->stats_buf_ctx->stats_in_end) {
      modified_error_total +=
          calculate_modified_err(frame_info, twopass, &(encode_context_ptr->two_pass_cfg), s);
      ++s;
    }
    twopass->modified_error_left = modified_error_total;
  }

  // Reset the vbr bits off target counters
  encode_context_ptr->rc.vbr_bits_off_target = 0;
  encode_context_ptr->rc.vbr_bits_off_target_fast = 0;
  encode_context_ptr->rc.rate_error_estimate = 0;

  // Static sequence monitor variables.
  twopass->kf_zeromotion_pct = 100;
  twopass->last_kfgroup_zeromotion_pct = 100;

  // Initialize bits per macro_block estimate correction factor.
  twopass->bpm_factor = 1.0;
  // Initialize actual and target bits counters for ARF groups so that
  // at the start we have a neutral bpm adjustment.
  twopass->rolling_arf_group_target_bits = 1;
  twopass->rolling_arf_group_actual_bits = 1;
}

int frame_is_kf_gf_arf(PictureParentControlSet *ppcs_ptr) {
  SequenceControlSet *scs_ptr = ppcs_ptr->scs_ptr;
  EncodeContext *encode_context_ptr = scs_ptr->encode_context_ptr;
  GF_GROUP *const gf_group = &encode_context_ptr->gf_group;
  const int update_type = gf_group->update_type[ppcs_ptr->gf_group_index];

  return frame_is_intra_only(ppcs_ptr) || update_type == ARF_UPDATE ||
         update_type == GF_UPDATE;
}

#define MINQ_ADJ_LIMIT 48
#define MINQ_ADJ_LIMIT_CQ 20
#define HIGH_UNDERSHOOT_RATIO 2
void svt_av1_twopass_postencode_update(PictureParentControlSet *ppcs_ptr) {
  SequenceControlSet *scs_ptr = ppcs_ptr->scs_ptr;
  EncodeContext *encode_context_ptr = scs_ptr->encode_context_ptr;
  RATE_CONTROL *const rc = &encode_context_ptr->rc;
  TWO_PASS *const twopass = &scs_ptr->twopass;
  GF_GROUP *const gf_group = &encode_context_ptr->gf_group;
  const int bits_used = rc->base_frame_target;
  const RateControlCfg *const rc_cfg = &encode_context_ptr->rc_cfg;

  // VBR correction is done through rc->vbr_bits_off_target. Based on the
  // sign of this value, a limited % adjustment is made to the target rate
  // of subsequent frames, to try and push it back towards 0. This method
  // is designed to prevent extreme behaviour at the end of a clip
  // or group of frames.
  rc->vbr_bits_off_target += rc->base_frame_target - rc->projected_frame_size;
  twopass->bits_left = AOMMAX(twopass->bits_left - bits_used, 0);
  // Target vs actual bits for this arf group.
  twopass->rolling_arf_group_target_bits += rc->this_frame_target;
  twopass->rolling_arf_group_actual_bits += rc->projected_frame_size;

  // Calculate the pct rc error.
  if (rc->total_actual_bits) {
    rc->rate_error_estimate =
        (int)((rc->vbr_bits_off_target * 100) / rc->total_actual_bits);
    rc->rate_error_estimate = clamp(rc->rate_error_estimate, -100, 100);
  } else {
    rc->rate_error_estimate = 0;
  }

  // Update the active best quality pyramid.
  if (!rc->is_src_frame_alt_ref) {
      const int pyramid_level = gf_group->layer_depth[ppcs_ptr->gf_group_index];
    int i;
    for (i = pyramid_level; i <= MAX_ARF_LAYERS; ++i) {
      rc->active_best_quality[i] = ppcs_ptr->frm_hdr.quantization_params.base_q_idx;
    }
  }

  if (ppcs_ptr->frm_hdr.frame_type != KEY_FRAME) {
    twopass->kf_group_bits -= bits_used;
    twopass->last_kfgroup_zeromotion_pct = twopass->kf_zeromotion_pct;
  }
  twopass->kf_group_bits = AOMMAX(twopass->kf_group_bits, 0);

  // If the rate control is drifting consider adjustment to min or maxq.
  if ((rc_cfg->mode != AOM_Q) && !rc->is_src_frame_alt_ref) {
    const int maxq_adj_limit = rc->worst_quality - rc->active_worst_quality;
    const int minq_adj_limit =
        (rc_cfg->mode == AOM_CQ ? MINQ_ADJ_LIMIT_CQ : MINQ_ADJ_LIMIT);

    // Undershoot.
    if (rc->rate_error_estimate > rc_cfg->under_shoot_pct) {
      --twopass->extend_maxq;
      if (rc->rolling_target_bits >= rc->rolling_actual_bits)
        ++twopass->extend_minq;
      // Overshoot.
    } else if (rc->rate_error_estimate < -rc_cfg->over_shoot_pct) {
      --twopass->extend_minq;
      if (rc->rolling_target_bits < rc->rolling_actual_bits)
        ++twopass->extend_maxq;
    } else {
      // Adjustment for extreme local overshoot.
      if (rc->projected_frame_size > (2 * rc->base_frame_target) &&
          rc->projected_frame_size > (2 * rc->avg_frame_bandwidth))
        ++twopass->extend_maxq;

      // Unwind undershoot or overshoot adjustment.
      if (rc->rolling_target_bits < rc->rolling_actual_bits)
        --twopass->extend_minq;
      else if (rc->rolling_target_bits > rc->rolling_actual_bits)
        --twopass->extend_maxq;
    }

    twopass->extend_minq = clamp(twopass->extend_minq, 0, minq_adj_limit);
    twopass->extend_maxq = clamp(twopass->extend_maxq, 0, maxq_adj_limit);

    // If there is a big and undexpected undershoot then feed the extra
    // bits back in quickly. One situation where this may happen is if a
    // frame is unexpectedly almost perfectly predicted by the ARF or GF
    // but not very well predcited by the previous frame.
    if (!frame_is_kf_gf_arf(ppcs_ptr) && !rc->is_src_frame_alt_ref) {
      int fast_extra_thresh = rc->base_frame_target / HIGH_UNDERSHOOT_RATIO;
      if (rc->projected_frame_size < fast_extra_thresh) {
        rc->vbr_bits_off_target_fast +=
            fast_extra_thresh - rc->projected_frame_size;
        rc->vbr_bits_off_target_fast =
            AOMMIN(rc->vbr_bits_off_target_fast, (4 * rc->avg_frame_bandwidth));

        // Fast adaptation of minQ if necessary to use up the extra bits.
        if (rc->avg_frame_bandwidth) {
          twopass->extend_minq_fast =
              (int)(rc->vbr_bits_off_target_fast * 8 / rc->avg_frame_bandwidth);
        }
        twopass->extend_minq_fast = AOMMIN(
            twopass->extend_minq_fast, minq_adj_limit - twopass->extend_minq);
      } else if (rc->vbr_bits_off_target_fast) {
        twopass->extend_minq_fast = AOMMIN(
            twopass->extend_minq_fast, minq_adj_limit - twopass->extend_minq);
      } else {
        twopass->extend_minq_fast = 0;
      }
    }
  }
}
