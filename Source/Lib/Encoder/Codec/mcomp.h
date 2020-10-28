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

#ifndef AOM_AV1_ENCODER_MCOMP_H_
#define AOM_AV1_ENCODER_MCOMP_H_

#include "mv.h"
#include "EbCodingUnit.h"
#include "EbBlockStructures.h"
#include "Av1Common.h"
#include "av1me.h"
#include "EbRateDistortionCost.h"

#ifdef __cplusplus
extern "C" {
#endif
// =============================================================================
//  Cost functions
// =============================================================================

enum {
  MV_COST_ENTROPY,    // Use the entropy rate of the mv as the cost
  MV_COST_L1_LOWRES,  // Use the l1 norm of the mv as the cost (<480p)
  MV_COST_L1_MIDRES,  // Use the l1 norm of the mv as the cost (>=480p)
  MV_COST_L1_HDRES,   // Use the l1 norm of the mv as the cost (>=720p)
  MV_COST_NONE        // Use 0 as as cost irrespective of the current mv
} UENUM1BYTE(MV_COST_TYPE);

typedef struct {
  // The reference mv used to compute the mv cost
  const MV *ref_mv;
  FULLPEL_MV full_ref_mv;
  MV_COST_TYPE mv_cost_type;
  const int *mvjcost;
  const int *mvcost[2];
  int error_per_bit;
  // A multiplier used to convert rate to sad cost
  int sad_per_bit;
} MV_COST_PARAMS;

// =============================================================================
//  Motion Search
// =============================================================================
extern struct svt_buf_2d {
    uint8_t *buf;
    uint8_t *buf0;
    int width;
    int height;
    int stride;
} svt_buf_2d;

typedef struct {
  // The reference buffer
  struct svt_buf_2d *ref;

  // The source and predictors/mask used by translational search
  struct svt_buf_2d *src;
  uint8_t *second_pred;
  uint8_t *mask;
  int mask_stride;
  int inv_mask;

  // The weighted source and mask used by OBMC
  int32_t *wsrc;
  int32_t *obmc_mask;
} MSBuffers;

static INLINE void svt_av1_set_ms_compound_refs(MSBuffers *ms_buffers,
                                            uint8_t *second_pred,
                                            uint8_t *mask,
                                            int mask_stride, int invert_mask) {
  ms_buffers->second_pred = second_pred;
  ms_buffers->mask = mask;
  ms_buffers->mask_stride = mask_stride;
  ms_buffers->inv_mask = invert_mask;
}

// =============================================================================
//  Subpixel Motion Search
// =============================================================================
enum {
  EIGHTH_PEL,
  QUARTER_PEL,
  HALF_PEL,
  FULL_PEL
} UENUM1BYTE(SUBPEL_FORCE_STOP);

typedef struct {

  const AomVarianceFnPtr *vfp;

  SUBPEL_SEARCH_TYPE subpel_search_type;

  // Source and reference buffers
  MSBuffers ms_buffers;

  int w, h;
} SUBPEL_SEARCH_VAR_PARAMS;

// This struct holds subpixel motion search parameters that should be constant
// during the search
typedef struct {
  // High level motion search settings
  int allow_hp;
  const int *cost_list;
  SUBPEL_FORCE_STOP forced_stop;
  int iters_per_step;

  SubpelMvLimits mv_limits;

  // For calculating mv cost
  MV_COST_PARAMS mv_cost_params;

  // Distortion calculation params
  SUBPEL_SEARCH_VAR_PARAMS var_params;

} SUBPEL_MOTION_SEARCH_PARAMS;


typedef int(fractional_mv_step_fp)(MacroBlockD *xd, const struct AV1Common *const cm,
                                   const SUBPEL_MOTION_SEARCH_PARAMS *ms_params,
                                   MV start_mv, MV *bestmv, int *distortion,
                                   unsigned int *sse1,
                                   int_mv *last_mv_search_list);

extern fractional_mv_step_fp svt_av1_find_best_sub_pixel_tree;

static INLINE void svt_av1_set_subpel_mv_search_range(SubpelMvLimits *subpel_limits,
                                                  const FullMvLimits *mv_limits,
                                                  const MV *ref_mv) {
  const int max_mv = GET_MV_SUBPEL(MAX_FULL_PEL_VAL);
  const int minc =
      AOMMAX(GET_MV_SUBPEL(mv_limits->col_min), ref_mv->col - max_mv);
  const int maxc =
      AOMMIN(GET_MV_SUBPEL(mv_limits->col_max), ref_mv->col + max_mv);
  const int minr =
      AOMMAX(GET_MV_SUBPEL(mv_limits->row_min), ref_mv->row - max_mv);
  const int maxr =
      AOMMIN(GET_MV_SUBPEL(mv_limits->row_max), ref_mv->row + max_mv);

  subpel_limits->col_min = AOMMAX(MV_LOW + 1, minc);
  subpel_limits->col_max = AOMMIN(MV_UPP - 1, maxc);
  subpel_limits->row_min = AOMMAX(MV_LOW + 1, minr);
  subpel_limits->row_max = AOMMIN(MV_UPP - 1, maxr);
}

static INLINE int svt_av1_is_subpelmv_in_range(const SubpelMvLimits *mv_limits,
                                           MV mv) {
  return (mv.col >= mv_limits->col_min) && (mv.col <= mv_limits->col_max) &&
         (mv.row >= mv_limits->row_min) && (mv.row <= mv_limits->row_max);
}

// Returns the rate of encoding the current motion vector based on the
// joint_cost and comp_cost. joint_costs covers the cost of transmitting
// JOINT_MV, and comp_cost covers the cost of transmitting the actual motion
// vector.
static INLINE int svt_mv_cost(const MV *mv, const int *joint_cost, const int *const comp_cost[2]) {
    return joint_cost[svt_av1_get_mv_joint(mv)] + comp_cost[0][mv->row] + comp_cost[1][mv->col];
}

#ifdef __cplusplus
}  // extern "C"
#endif
#endif  // AOM_AV1_ENCODER_MCOMP_H_
