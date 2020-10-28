/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
 */

#ifndef AOM_AV1_ENCODER_ME_H_
#define AOM_AV1_ENCODER_ME_H_

#include "EbDefinitions.h"
#include "EbCodingUnit.h"

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
#define SEARCH_GRID_CENTER_8P (SEARCH_RANGE_8P * SEARCH_GRID_STRIDE_8P + SEARCH_RANGE_8P)

// motion search site
typedef struct SearchSite {
    MV  mv;
    int offset;
} SearchSite;

typedef struct SearchSiteConfig {
    SearchSite ss[8 * MAX_MVSEARCH_STEPS + 1];
    int        ss_count;
    int        searches_per_step;
} SearchSiteConfig;

typedef unsigned int (*AomObmcSadFn)(const uint8_t *pred, int pred_stride, const int32_t *wsrc,
                                          const int32_t *msk);
typedef unsigned int (*AomObmcVarianceFn)(const uint8_t *pred, int pred_stride,
                                               const int32_t *wsrc, const int32_t *msk,
                                               unsigned int *sse);
typedef unsigned int (*AomObmcSubpixvarianceFn)(const uint8_t *pred, int pred_stride,
                                                     int xoffset, int yoffset, const int32_t *wsrc,
                                                     const int32_t *msk, unsigned int *sse);
typedef unsigned int (*AomSadFn)(const uint8_t *a, int a_stride, const uint8_t *b,
                                     int b_stride);

typedef unsigned int (*AomVarianceFn)(const uint8_t *a, int a_stride, const uint8_t *b,
                                          int b_stride, unsigned int *sse);

typedef void (*AomSadMultiDFn)(const uint8_t *a, int a_stride, const uint8_t *const b_array[],
                                     int b_stride, unsigned int *sad_array);

typedef struct aom_variance_vtable {
    AomSadFn                 sdf;
    AomVarianceFn            vf;
    AomVarianceFn        vf_hbd_10;
    AomSadMultiDFn         sdx4df;
    AomObmcSadFn            osdf;
    AomObmcVarianceFn       ovf;
    AomObmcSubpixvarianceFn osvf;

} AomVarianceFnPtr;

void av1_init_dsmotion_compensation(SearchSiteConfig *cfg, int stride);
void svt_av1_init3smotion_compensation(SearchSiteConfig *cfg, int stride);
void svt_av1_set_mv_search_range(MvLimits *mv_limits, const MV *mv);
struct Av1Comp;
struct SpeedFeatures;

int svt_av1_full_pixel_search(struct PictureControlSet *pcs, IntraBcContext /*MACROBLOCK*/ *x,
                              BlockSize bsize, MV *mvp_full, int step_param, int method,
                              int run_mesh_search, int error_per_bit, int *cost_list,
                              const MV *ref_mv, int var_max, int rd, int x_pos, int y_pos,
                              int intra);
int mv_err_cost(const MV *mv, const MV *ref, const int *mvjcost, int *mvcost[2], int error_per_bit);
#ifdef __cplusplus
} // extern "C"
#endif

#endif // AOM_AV1_ENCODER_ME_H_
