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

#include "EbDefinitions.h"
#include "EbCodingUnit.h"
#include "EbUtility.h"



//#include "av1/encoder/block.h"
//#include "aom_dsp/variance.h"

#ifdef __cplusplus
extern "C" {
#endif

// The maximum number of steps in a step search given the largest
// allowed initial step
#define MAX_MVSEARCH_STEPS 11
// Max full pel mv specified in the unit of full pixel
// Enable the use of motion vector in range [-1023, 1023].
#define MAX_FULL_PEL_VAL ((1 << (MAX_MVSEARCH_STEPS - 1)) - 1)
// Maximum size of the first step in full pel units
#define MAX_FIRST_STEP (1 << (MAX_MVSEARCH_STEPS - 1))
// Allowed motion vector pixel distance outside image border
// for Block_16x16
#define BORDER_MV_PIXELS_B16 (16 + AOM_INTERP_EXTEND)

#define SEARCH_RANGE_8P 3
#define SEARCH_GRID_STRIDE_8P (2 * SEARCH_RANGE_8P + 1)
#define SEARCH_GRID_CENTER_8P \
  (SEARCH_RANGE_8P * SEARCH_GRID_STRIDE_8P + SEARCH_RANGE_8P)

// motion search site
typedef struct search_site {
  MV mv;
  int offset;
} search_site;

typedef struct search_site_config {
  search_site ss[8 * MAX_MVSEARCH_STEPS + 1];
  int ss_count;
  int searches_per_step;
} search_site_config;

typedef struct {
  MV coord;
  int coord_offset;
} search_neighbors;


void av1_init_dsmotion_compensation(search_site_config *cfg, int stride);
void av1_init3smotion_compensation(search_site_config *cfg, int stride);

void av1_set_mv_search_range(MvLimits *mv_limits, const MV *mv);

#if 1 //---CHKN


struct AV1_COMP;
struct SPEED_FEATURES;


int av1_full_pixel_search(struct PictureControlSet_s *pcs, IntraBcContext /*MACROBLOCK*/ *x,
                          block_size bsize, MV *mvp_full, int step_param,
                          int method, int run_mesh_search, int error_per_bit,
                          int *cost_list, const MV *ref_mv, int var_max, int rd,
                          int x_pos, int y_pos, int intra);


#ifdef __cplusplus
}  // extern "C"
#endif

#endif//----CHKN
#endif  // AOM_AV1_ENCODER_MCOMP_H_
