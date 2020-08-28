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

#ifndef AOM_AV1_ENCODER_PASS2_STRATEGY_H_
#define AOM_AV1_ENCODER_PASS2_STRATEGY_H_

#include "level.h"
#include "encoder.h"
#include "Av1Common.h"

#if TWOPASS_RC
#ifdef __cplusplus
extern "C" {
#endif

// structure of accumulated stats and features in a gf group
typedef struct {
  double gf_group_err;
  double gf_group_raw_error;
  double gf_group_skip_pct;
  double gf_group_inactive_zone_rows;

  double mv_ratio_accumulator;
  double decay_accumulator;
  double zero_motion_accumulator;
  double loop_decay_rate;
  double last_loop_decay_rate;
  double this_frame_mv_in_out;
  double mv_in_out_accumulator;
  double abs_mv_in_out_accumulator;

  double avg_sr_coded_error;
  double avg_tr_coded_error;
  double avg_pcnt_second_ref;
  double avg_pcnt_third_ref;
  double avg_pcnt_third_ref_nolast;
  double avg_new_mv_count;
  double avg_wavelet_energy;
  double avg_raw_err_stdev;
  int non_zero_stdev_count;
} GF_GROUP_STATS;

typedef struct {
  double frame_err;
  double frame_coded_error;
  double frame_sr_coded_error;
  double frame_tr_coded_error;
} GF_FRAME_STATS;

void svt_av1_init_second_pass(struct SequenceControlSet *scs_ptr);

void av1_init_single_pass_lap(AV1_COMP *cpi);

void svt_av1_get_second_pass_params(struct PictureParentControlSet *pcs_ptr);

void svt_av1_twopass_postencode_update(struct PictureParentControlSet *ppcs_ptr);

int frame_is_kf_gf_arf(PictureParentControlSet *ppcs_ptr);
#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // TWOPASS_RC
#endif  // AOM_AV1_ENCODER_PASS2_STRATEGY_H_
