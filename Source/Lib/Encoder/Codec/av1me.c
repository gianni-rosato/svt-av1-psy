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

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include "av1me.h"
#include "mcomp.h"
#include "EbUtility.h"
#include "EbPictureControlSet.h"
#include "EbSequenceControlSet.h"
#include "aom_dsp_rtcd.h"
#include "EbModeDecisionProcess.h"
#include "EbAdaptiveMotionVectorPrediction.h"

int av1_is_dv_valid(const MV dv, const MacroBlockD *xd, int mi_row, int mi_col, BlockSize bsize,
                    int mib_size_log2);

int svt_av1_refining_search_sad(IntraBcContext *x, MV *ref_mv, int error_per_bit, int search_range,
                                const AomVarianceFnPtr *fn_ptr, const MV *center_mv);

AomVarianceFnPtr mefn_ptr[BlockSizeS_ALL];

void init_fn_ptr(void) {
#define BFP0(BT, SDF, VF, VF_HBD_10, SDX4DF)       \
    mefn_ptr[BT].sdf    = SDF;                     \
    mefn_ptr[BT].vf     = VF;                      \
    mefn_ptr[BT].vf_hbd_10 = VF_HBD_10;            \
    mefn_ptr[BT].sdx4df = SDX4DF;

    BFP0(BLOCK_4X16, svt_aom_sad4x16, svt_aom_variance4x16, svt_aom_highbd_10_variance4x16, svt_aom_sad4x16x4d)
    BFP0(BLOCK_16X4, svt_aom_sad16x4, svt_aom_variance16x4, svt_aom_highbd_10_variance16x4, svt_aom_sad16x4x4d)
    BFP0(BLOCK_8X32, svt_aom_sad8x32, svt_aom_variance8x32, svt_aom_highbd_10_variance8x32, svt_aom_sad8x32x4d)
    BFP0(BLOCK_32X8, svt_aom_sad32x8, svt_aom_variance32x8, svt_aom_highbd_10_variance32x8, svt_aom_sad32x8x4d)
    BFP0(BLOCK_16X64, svt_aom_sad16x64, svt_aom_variance16x64, svt_aom_highbd_10_variance16x64, svt_aom_sad16x64x4d)
    BFP0(BLOCK_64X16, svt_aom_sad64x16, svt_aom_variance64x16, svt_aom_highbd_10_variance64x16, svt_aom_sad64x16x4d)
    BFP0(BLOCK_128X128, svt_aom_sad128x128, svt_aom_variance128x128, svt_aom_highbd_10_variance128x128, svt_aom_sad128x128x4d)
    BFP0(BLOCK_128X64, svt_aom_sad128x64, svt_aom_variance128x64, svt_aom_highbd_10_variance128x64, svt_aom_sad128x64x4d)
    BFP0(BLOCK_64X128, svt_aom_sad64x128, svt_aom_variance64x128, svt_aom_highbd_10_variance64x128, svt_aom_sad64x128x4d)
    BFP0(BLOCK_32X16, svt_aom_sad32x16, svt_aom_variance32x16, svt_aom_highbd_10_variance32x16, svt_aom_sad32x16x4d)
    BFP0(BLOCK_16X32, svt_aom_sad16x32, svt_aom_variance16x32, svt_aom_highbd_10_variance16x32, svt_aom_sad16x32x4d)
    BFP0(BLOCK_64X32, svt_aom_sad64x32, svt_aom_variance64x32, svt_aom_highbd_10_variance64x32, svt_aom_sad64x32x4d)
    BFP0(BLOCK_32X64, svt_aom_sad32x64, svt_aom_variance32x64, svt_aom_highbd_10_variance32x64, svt_aom_sad32x64x4d)
    BFP0(BLOCK_32X32, svt_aom_sad32x32, svt_aom_variance32x32, svt_aom_highbd_10_variance32x32, svt_aom_sad32x32x4d)
    BFP0(BLOCK_64X64, svt_aom_sad64x64, svt_aom_variance64x64, svt_aom_highbd_10_variance64x64, svt_aom_sad64x64x4d)
    BFP0(BLOCK_16X16, svt_aom_sad16x16, svt_aom_variance16x16, svt_aom_highbd_10_variance16x16, svt_aom_sad16x16x4d)
    BFP0(BLOCK_16X8, svt_aom_sad16x8, svt_aom_variance16x8, svt_aom_highbd_10_variance16x8, svt_aom_sad16x8x4d)
    BFP0(BLOCK_8X16, svt_aom_sad8x16, svt_aom_variance8x16, svt_aom_highbd_10_variance8x16, svt_aom_sad8x16x4d)
    BFP0(BLOCK_8X8, svt_aom_sad8x8, svt_aom_variance8x8, svt_aom_highbd_10_variance8x8, svt_aom_sad8x8x4d)
    BFP0(BLOCK_8X4, svt_aom_sad8x4, svt_aom_variance8x4, svt_aom_highbd_10_variance8x4, svt_aom_sad8x4x4d)
    BFP0(BLOCK_4X8, svt_aom_sad4x8, svt_aom_variance4x8, svt_aom_highbd_10_variance4x8, svt_aom_sad4x8x4d)
    BFP0(BLOCK_4X4, svt_aom_sad4x4, svt_aom_variance4x4, svt_aom_highbd_10_variance4x4, svt_aom_sad4x4x4d)
#define OBFP(BT, OSDF, OVF, OSVF) \
    mefn_ptr[BT].osdf = OSDF;     \
    mefn_ptr[BT].ovf  = OVF;      \
    mefn_ptr[BT].osvf = OSVF;
    OBFP(BLOCK_128X128,
         svt_aom_obmc_sad128x128,
         svt_aom_obmc_variance128x128,
         svt_aom_obmc_sub_pixel_variance128x128)
    OBFP(BLOCK_128X64,
         svt_aom_obmc_sad128x64,
         svt_aom_obmc_variance128x64,
         svt_aom_obmc_sub_pixel_variance128x64)
    OBFP(BLOCK_64X128,
         svt_aom_obmc_sad64x128,
         svt_aom_obmc_variance64x128,
         svt_aom_obmc_sub_pixel_variance64x128)
    OBFP(BLOCK_64X64, svt_aom_obmc_sad64x64, svt_aom_obmc_variance64x64, svt_aom_obmc_sub_pixel_variance64x64)
    OBFP(BLOCK_64X32, svt_aom_obmc_sad64x32, svt_aom_obmc_variance64x32, svt_aom_obmc_sub_pixel_variance64x32)
    OBFP(BLOCK_32X64, svt_aom_obmc_sad32x64, svt_aom_obmc_variance32x64, svt_aom_obmc_sub_pixel_variance32x64)
    OBFP(BLOCK_32X32, svt_aom_obmc_sad32x32, svt_aom_obmc_variance32x32, svt_aom_obmc_sub_pixel_variance32x32)
    OBFP(BLOCK_32X16, svt_aom_obmc_sad32x16, svt_aom_obmc_variance32x16, svt_aom_obmc_sub_pixel_variance32x16)
    OBFP(BLOCK_16X32, svt_aom_obmc_sad16x32, svt_aom_obmc_variance16x32, svt_aom_obmc_sub_pixel_variance16x32)
    OBFP(BLOCK_16X16, svt_aom_obmc_sad16x16, svt_aom_obmc_variance16x16, svt_aom_obmc_sub_pixel_variance16x16)
    OBFP(BLOCK_16X8, svt_aom_obmc_sad16x8, svt_aom_obmc_variance16x8, svt_aom_obmc_sub_pixel_variance16x8)
    OBFP(BLOCK_8X16, svt_aom_obmc_sad8x16, svt_aom_obmc_variance8x16, svt_aom_obmc_sub_pixel_variance8x16)
    OBFP(BLOCK_8X8, svt_aom_obmc_sad8x8, svt_aom_obmc_variance8x8, svt_aom_obmc_sub_pixel_variance8x8)
    OBFP(BLOCK_4X8, svt_aom_obmc_sad4x8, svt_aom_obmc_variance4x8, svt_aom_obmc_sub_pixel_variance4x8)
    OBFP(BLOCK_8X4, svt_aom_obmc_sad8x4, svt_aom_obmc_variance8x4, svt_aom_obmc_sub_pixel_variance8x4)
    OBFP(BLOCK_4X4, svt_aom_obmc_sad4x4, svt_aom_obmc_variance4x4, svt_aom_obmc_sub_pixel_variance4x4)
    OBFP(BLOCK_4X16, svt_aom_obmc_sad4x16, svt_aom_obmc_variance4x16, svt_aom_obmc_sub_pixel_variance4x16)
    OBFP(BLOCK_16X4, svt_aom_obmc_sad16x4, svt_aom_obmc_variance16x4, svt_aom_obmc_sub_pixel_variance16x4)
    OBFP(BLOCK_8X32, svt_aom_obmc_sad8x32, svt_aom_obmc_variance8x32, svt_aom_obmc_sub_pixel_variance8x32)
    OBFP(BLOCK_32X8, svt_aom_obmc_sad32x8, svt_aom_obmc_variance32x8, svt_aom_obmc_sub_pixel_variance32x8)
    OBFP(BLOCK_16X64, svt_aom_obmc_sad16x64, svt_aom_obmc_variance16x64, svt_aom_obmc_sub_pixel_variance16x64)
    OBFP(BLOCK_64X16, svt_aom_obmc_sad64x16, svt_aom_obmc_variance64x16, svt_aom_obmc_sub_pixel_variance64x16)
}

// #define NEW_DIAMOND_SEARCH

static INLINE const uint8_t *get_buf_from_mv(const struct Buf2D *buf, const MV *mv) {
    return &buf->buf[mv->row * buf->stride + mv->col];
}

void svt_av1_set_mv_search_range(MvLimits *mv_limits, const MV *mv) {
    int col_min = (mv->col >> 3) - MAX_FULL_PEL_VAL + !!(mv->col & 7);
    int row_min = (mv->row >> 3) - MAX_FULL_PEL_VAL + !!(mv->row & 7);
    int col_max = (mv->col >> 3) + MAX_FULL_PEL_VAL;
    int row_max = (mv->row >> 3) + MAX_FULL_PEL_VAL;

    col_min = AOMMAX(col_min, (MV_LOW >> 3) + 1);
    row_min = AOMMAX(row_min, (MV_LOW >> 3) + 1);
    col_max = AOMMIN(col_max, (MV_UPP >> 3) - 1);
    row_max = AOMMIN(row_max, (MV_UPP >> 3) - 1);

    // Get intersection of UMV window and valid MV window to reduce # of checks
    // in diamond search.
    if (mv_limits->col_min < col_min) mv_limits->col_min = col_min;
    if (mv_limits->col_max > col_max) mv_limits->col_max = col_max;
    if (mv_limits->row_min < row_min) mv_limits->row_min = row_min;
    if (mv_limits->row_max > row_max) mv_limits->row_max = row_max;
}

#define PIXEL_TRANSFORM_ERROR_SCALE 4
int mv_err_cost(const MV *mv, const MV *ref, const int *mvjcost, int *mvcost[2],
                       int error_per_bit) {
    if (mvcost) {
        const MV diff = {mv->row - ref->row, mv->col - ref->col};
        return (int)ROUND_POWER_OF_TWO_64(
            (int64_t)svt_mv_cost(&diff, mvjcost, (const int* const *)mvcost) * error_per_bit,
            RDDIV_BITS + AV1_PROB_COST_SHIFT - RD_EPB_SHIFT + PIXEL_TRANSFORM_ERROR_SCALE);
    }
    return 0;
}

static int mvsad_err_cost(const IntraBcContext *x, const MV *mv, const MV *ref, int sad_per_bit) {
    const MV diff = {(mv->row - ref->row) * 8, (mv->col - ref->col) * 8};
    return ROUND_POWER_OF_TWO(
        (unsigned)svt_mv_cost(&diff, x->nmv_vec_cost, (const int *const *)x->mv_cost_stack) *
            sad_per_bit, AV1_PROB_COST_SHIFT);
}

void svt_av1_init3smotion_compensation(SearchSiteConfig *cfg, int stride) {
    int len, ss_count = 1;

    cfg->ss[0].mv.col = cfg->ss[0].mv.row = 0;
    cfg->ss[0].offset                     = 0;

    for (len = MAX_FIRST_STEP; len > 0; len /= 2) {
        // Generate offsets for 8 search sites per step.
        const MV ss_mvs[8] = {{-len, 0},
                              {len, 0},
                              {0, -len},
                              {0, len},
                              {-len, -len},
                              {-len, len},
                              {len, -len},
                              {len, len}};
        int      i;
        for (i = 0; i < 8; ++i) {
            SearchSite *const ss = &cfg->ss[ss_count++];
            ss->mv               = ss_mvs[i];
            ss->offset           = ss->mv.row * stride + ss->mv.col;
        }
    }

    cfg->ss_count          = ss_count;
    cfg->searches_per_step = 8;
}

static INLINE int is_mv_in(const MvLimits *mv_limits, const MV *mv) {
    return (mv->col >= mv_limits->col_min) && (mv->col <= mv_limits->col_max) &&
           (mv->row >= mv_limits->row_min) && (mv->row <= mv_limits->row_max);
}
#define MAX_PATTERN_SCALES 11
#define MAX_PATTERN_CANDIDATES 8 // max number of canddiates per scale
#define PATTERN_CANDIDATES_REF 3 // number of refinement candidates

int svt_av1_get_mvpred_var(const IntraBcContext *x, const MV *best_mv, const MV *center_mv,
                           const AomVarianceFnPtr *vfp, int use_mvcost) {
    const struct Buf2D *const what    = &x->plane[0].src;
    const struct Buf2D *const in_what = &x->xdplane[0].pre[0];
    const MV                  mv      = {best_mv->row * 8, best_mv->col * 8};
    unsigned int              unused;

    return vfp->vf(what->buf,
                   what->stride,
                   get_buf_from_mv(in_what, best_mv),
                   in_what->stride,
                   &unused) +
           (use_mvcost
                ? mv_err_cost(&mv, center_mv, x->nmv_vec_cost, x->mv_cost_stack, x->errorperbit)
                : 0);
}

// Exhuastive motion search around a given centre position with a given
// step size.
static int exhuastive_mesh_search(IntraBcContext *x, MV *ref_mv, MV *best_mv, int range, int step,
                                  int sad_per_bit, const AomVarianceFnPtr *fn_ptr,
                                  const MV *center_mv) {
    const struct Buf2D *const what       = &x->plane[0].src;
    const struct Buf2D *const in_what    = &x->xdplane[0].pre[0];
    MV                        fcenter_mv = {center_mv->row, center_mv->col};
    unsigned int              best_sad   = INT_MAX;
    int                       r, c, i;
    int                       start_col, end_col, start_row, end_row;
    int                       col_step = (step > 1) ? step : 4;

    assert(step >= 1);

    clamp_mv(&fcenter_mv,
             x->mv_limits.col_min,
             x->mv_limits.col_max,
             x->mv_limits.row_min,
             x->mv_limits.row_max);
    *best_mv = fcenter_mv;
    best_sad =
        fn_ptr->sdf(
            what->buf, what->stride, get_buf_from_mv(in_what, &fcenter_mv), in_what->stride) +
        mvsad_err_cost(x, &fcenter_mv, ref_mv, sad_per_bit);
    start_row = AOMMAX(-range, x->mv_limits.row_min - fcenter_mv.row);
    start_col = AOMMAX(-range, x->mv_limits.col_min - fcenter_mv.col);
    end_row   = AOMMIN(range, x->mv_limits.row_max - fcenter_mv.row);
    end_col   = AOMMIN(range, x->mv_limits.col_max - fcenter_mv.col);

    for (r = start_row; r <= end_row; r += step) {
        for (c = start_col; c <= end_col; c += col_step) {
            // Step > 1 means we are not checking every location in this pass.
            if (step > 1) {
                const MV     mv  = {fcenter_mv.row + r, fcenter_mv.col + c};
                unsigned int sad = fn_ptr->sdf(
                    what->buf, what->stride, get_buf_from_mv(in_what, &mv), in_what->stride);
                if (sad < best_sad) {
                    sad += mvsad_err_cost(x, &mv, ref_mv, sad_per_bit);
                    if (sad < best_sad) {
                        best_sad                = sad;
                        x->second_best_mv.as_mv = *best_mv;
                        *best_mv                = mv;
                    }
                }
            } else {
                // 4 sads in a single call if we are checking every location
                if (c + 3 <= end_col) {
                    unsigned int   sads[4];
                    const uint8_t *addrs[4];
                    for (i = 0; i < 4; ++i) {
                        const MV mv = {fcenter_mv.row + r, fcenter_mv.col + c + i};
                        addrs[i]    = get_buf_from_mv(in_what, &mv);
                    }
                    fn_ptr->sdx4df(what->buf, what->stride, addrs, in_what->stride, sads);

                    for (i = 0; i < 4; ++i) {
                        if (sads[i] < best_sad) {
                            const MV           mv = {fcenter_mv.row + r, fcenter_mv.col + c + i};
                            const unsigned int sad =
                                sads[i] + mvsad_err_cost(x, &mv, ref_mv, sad_per_bit);
                            if (sad < best_sad) {
                                best_sad                = sad;
                                x->second_best_mv.as_mv = *best_mv;
                                *best_mv                = mv;
                            }
                        }
                    }
                } else {
                    for (i = 0; i < end_col - c; ++i) {
                        const MV     mv  = {fcenter_mv.row + r, fcenter_mv.col + c + i};
                        unsigned int sad = fn_ptr->sdf(what->buf,
                                                       what->stride,
                                                       get_buf_from_mv(in_what, &mv),
                                                       in_what->stride);
                        if (sad < best_sad) {
                            sad += mvsad_err_cost(x, &mv, ref_mv, sad_per_bit);
                            if (sad < best_sad) {
                                best_sad                = sad;
                                x->second_best_mv.as_mv = *best_mv;
                                *best_mv                = mv;
                            }
                        }
                    }
                }
            }
        }
    }

    return best_sad;
}

int svt_av1_diamond_search_sad_c(IntraBcContext *x, const SearchSiteConfig *cfg, MV *ref_mv,
                                 MV *best_mv, int search_param, int sad_per_bit, int *num00,
                                 const AomVarianceFnPtr *fn_ptr, const MV *center_mv) {
    int i, j, step;

    uint8_t *      what        = x->plane[0].src.buf;
    const int      what_stride = x->plane[0].src.stride;
    const uint8_t *in_what;
    const int      in_what_stride = x->xdplane[0].pre[0].stride;
    const uint8_t *best_address;

    unsigned int bestsad;
    int          best_site = 0;
    int          last_site = 0;

    int ref_row;
    int ref_col;

    // search_param determines the length of the initial step and hence the number
    // of iterations.
    // 0 = initial step (MAX_FIRST_STEP) pel
    // 1 = (MAX_FIRST_STEP/2) pel,
    // 2 = (MAX_FIRST_STEP/4) pel...
    const SearchSite *ss        = &cfg->ss[search_param * cfg->searches_per_step];
    const int         tot_steps = (cfg->ss_count / cfg->searches_per_step) - search_param;

    const MV fcenter_mv = {center_mv->row >> 3, center_mv->col >> 3};
    clamp_mv(ref_mv,
             x->mv_limits.col_min,
             x->mv_limits.col_max,
             x->mv_limits.row_min,
             x->mv_limits.row_max);
    ref_row      = ref_mv->row;
    ref_col      = ref_mv->col;
    *num00       = 0;
    best_mv->row = ref_row;
    best_mv->col = ref_col;

    // Work out the start point for the search
    in_what      = x->xdplane[0].pre[0].buf + ref_row * in_what_stride + ref_col;
    best_address = in_what;

    // Check the starting position
    bestsad = fn_ptr->sdf(what, what_stride, in_what, in_what_stride) +
              mvsad_err_cost(x, best_mv, &fcenter_mv, sad_per_bit);

    i = 1;

    for (step = 0; step < tot_steps; step++) {
        int all_in = 1;

        // All_in is true if every one of the points we are checking are within
        // the bounds of the image.
        all_in &= ((best_mv->row + ss[i].mv.row) > x->mv_limits.row_min);
        all_in &= ((best_mv->row + ss[i + 1].mv.row) < x->mv_limits.row_max);
        all_in &= ((best_mv->col + ss[i + 2].mv.col) > x->mv_limits.col_min);
        all_in &= ((best_mv->col + ss[i + 3].mv.col) < x->mv_limits.col_max);

        // If all the pixels are within the bounds we don't check whether the
        // search point is valid in this loop,  otherwise we check each point
        // for validity..
        if (all_in) {
            unsigned int sad_array[4];

            for (j = 0; j < cfg->searches_per_step; j += 4) {
                unsigned char const *block_offset[4];

                for (int t = 0; t < 4; t++) block_offset[t] = ss[i + t].offset + best_address;

                fn_ptr->sdx4df(what, what_stride, block_offset, in_what_stride, sad_array);

                for (int t = 0; t < 4; t++, i++) {
                    if (sad_array[t] < bestsad) {
                        const MV this_mv = {best_mv->row + ss[i].mv.row,
                                            best_mv->col + ss[i].mv.col};
                        sad_array[t] += mvsad_err_cost(x, &this_mv, &fcenter_mv, sad_per_bit);
                        if (sad_array[t] < bestsad) {
                            bestsad   = sad_array[t];
                            best_site = i;
                        }
                    }
                }
            }
        } else {
            for (j = 0; j < cfg->searches_per_step; j++) {
                // Trap illegal vectors
                const MV this_mv = {best_mv->row + ss[i].mv.row, best_mv->col + ss[i].mv.col};

                if (is_mv_in(&x->mv_limits, &this_mv)) {
                    const uint8_t *const check_here = ss[i].offset + best_address;
                    unsigned int         thissad =
                        fn_ptr->sdf(what, what_stride, check_here, in_what_stride);

                    if (thissad < bestsad) {
                        thissad += mvsad_err_cost(x, &this_mv, &fcenter_mv, sad_per_bit);
                        if (thissad < bestsad) {
                            bestsad   = thissad;
                            best_site = i;
                        }
                    }
                }
                i++;
            }
        }
        if (best_site != last_site) {
            x->second_best_mv.as_mv = *best_mv;
            best_mv->row += ss[best_site].mv.row;
            best_mv->col += ss[best_site].mv.col;
            best_address += ss[best_site].offset;
            last_site = best_site;
#if defined(NEW_DIAMOND_SEARCH)
            while (1) {
                const MV this_mv = {best_mv->row + ss[best_site].mv.row,
                                    best_mv->col + ss[best_site].mv.col};
                if (is_mv_in(&x->mv_limits, &this_mv)) {
                    const uint8_t *const check_here = ss[best_site].offset + best_address;
                    unsigned int         thissad =
                        fn_ptr->sdf(what, what_stride, check_here, in_what_stride);
                    if (thissad < bestsad) {
                        thissad += mvsad_err_cost(x, &this_mv, &fcenter_mv, sad_per_bit);
                        if (thissad < bestsad) {
                            bestsad = thissad;
                            best_mv->row += ss[best_site].mv.row;
                            best_mv->col += ss[best_site].mv.col;
                            best_address += ss[best_site].offset;
                            continue;
                        }
                    }
                }
                break;
            }
#endif
        } else if (best_address == in_what)
            (*num00)++;
    }
    return bestsad;
}

/* do_refine: If last step (1-away) of n-step search doesn't pick the center
              point as the best match, we will do a final 1-away diamond
              refining search  */
static int full_pixel_diamond(PictureControlSet *pcs, IntraBcContext /*MACROBLOCK*/ *x,
                              MV *mvp_full, int step_param, int sadpb, int further_steps,
                              int do_refine, int *cost_list, const AomVarianceFnPtr *fn_ptr,
                              const MV *ref_mv) {
    MV  temp_mv;
    int thissme, n, num00 = 0;
    (void)cost_list;
    /*int bestsme = cpi->diamond_search_sad(x, &cpi->ss_cfg, mvp_full, &temp_mv,
                                        step_param, sadpb, &n, fn_ptr, ref_mv);*/
    int bestsme = svt_av1_diamond_search_sad_c(
        x, &pcs->ss_cfg, mvp_full, &temp_mv, step_param, sadpb, &n, fn_ptr, ref_mv);

    if (bestsme < INT_MAX) bestsme = svt_av1_get_mvpred_var(x, &temp_mv, ref_mv, fn_ptr, 1);
    x->best_mv.as_mv = temp_mv;

    // If there won't be more n-step search, check to see if refining search is
    // needed.
    if (n > further_steps) do_refine = 0;

    while (n < further_steps) {
        ++n;

        if (num00) {
            num00--;
        } else {
            /*thissme = cpi->diamond_search_sad(x, &cpi->ss_cfg, mvp_full, &temp_mv,
                                        step_param + n, sadpb, &num00, fn_ptr,
                                        ref_mv);*/
            thissme = svt_av1_diamond_search_sad_c(
                x, &pcs->ss_cfg, mvp_full, &temp_mv, step_param + n, sadpb, &num00, fn_ptr, ref_mv);

            if (thissme < INT_MAX) thissme = svt_av1_get_mvpred_var(x, &temp_mv, ref_mv, fn_ptr, 1);

            // check to see if refining search is needed.
            if (num00 > further_steps - n) do_refine = 0;

            if (thissme < bestsme) {
                bestsme          = thissme;
                x->best_mv.as_mv = temp_mv;
            }
        }
    }

    // final 1-away diamond refining search
    if (do_refine) {
        const int search_range = 8;
        MV        best_mv      = x->best_mv.as_mv;
        thissme = svt_av1_refining_search_sad(x, &best_mv, sadpb, search_range, fn_ptr, ref_mv);
        if (thissme < INT_MAX) thissme = svt_av1_get_mvpred_var(x, &best_mv, ref_mv, fn_ptr, 1);
        if (thissme < bestsme) {
            bestsme          = thissme;
            x->best_mv.as_mv = best_mv;
        }
    }

    // Return cost list.
    /* if (cost_list) {
    calc_int_cost_list(x, ref_mv, sadpb, fn_ptr, &x->best_mv.as_mv, cost_list);
  }*/
    return bestsme;
}

#define MIN_RANGE 7
#define MAX_RANGE 256
#define MIN_INTERVAL 1
// Runs an limited range exhaustive mesh search using a pattern set
// according to the encode speed profile.
static int full_pixel_exhaustive(PictureControlSet *pcs, IntraBcContext *x,
                                 const MV *centre_mv_full, int sadpb, int *cost_list,
                                 const AomVarianceFnPtr *fn_ptr, const MV *ref_mv,
                                 MV *dst_mv) {
    UNUSED(cost_list);
    const SpeedFeatures *const sf       = &pcs->sf; // cpi->sf;
    MV                         temp_mv  = {centre_mv_full->row, centre_mv_full->col};
    MV                         f_ref_mv = {ref_mv->row >> 3, ref_mv->col >> 3};
    int                        bestsme;
    int                        interval = sf->mesh_patterns[0].interval;
    int                        range    = sf->mesh_patterns[0].range;
    int                        baseline_interval_divisor;

    // Keep track of number of exhaustive calls (this frame in this thread).
    //CHKN if (x->ex_search_count_ptr != NULL) ++(*x->ex_search_count_ptr);

    // Trap illegal values for interval and range for this function.
    if ((range < MIN_RANGE) || (range > MAX_RANGE) || (interval < MIN_INTERVAL) ||
        (interval > range))
        return INT_MAX;

    baseline_interval_divisor = range / interval;

    // Check size of proposed first range against magnitude of the centre
    // value used as a starting point.
    range    = AOMMAX(range, (5 * AOMMAX(abs(temp_mv.row), abs(temp_mv.col))) / 4);
    range    = AOMMIN(range, MAX_RANGE);
    interval = AOMMAX(interval, range / baseline_interval_divisor);

    // initial search
    bestsme =
        exhuastive_mesh_search(x, &f_ref_mv, &temp_mv, range, interval, sadpb, fn_ptr, &temp_mv);

    if ((interval > MIN_INTERVAL) && (range > MIN_RANGE)) {
        // Progressive searches with range and step size decreasing each time
        // till we reach a step size of 1. Then break out.
        for (int i = 1; i < MAX_MESH_STEP; ++i) {
            // First pass with coarser step and longer range
            bestsme = exhuastive_mesh_search(x,
                                             &f_ref_mv,
                                             &temp_mv,
                                             sf->mesh_patterns[i].range,
                                             sf->mesh_patterns[i].interval,
                                             sadpb,
                                             fn_ptr,
                                             &temp_mv);

            if (sf->mesh_patterns[i].interval == 1) break;
        }
    }

    if (bestsme < INT_MAX) bestsme = svt_av1_get_mvpred_var(x, &temp_mv, ref_mv, fn_ptr, 1);
    *dst_mv = temp_mv;

    // Return cost list.
    /* if (cost_list) {
    calc_int_cost_list(x, ref_mv, sadpb, fn_ptr, dst_mv, cost_list);
  }*/
    return bestsme;
}

int svt_av1_refining_search_sad(IntraBcContext *x, MV *ref_mv, int error_per_bit, int search_range,
                                const AomVarianceFnPtr *fn_ptr, const MV *center_mv) {
    const MV                  neighbors[4] = {{-1, 0}, {0, -1}, {0, 1}, {1, 0}};
    const struct Buf2D *const what         = &x->plane[0].src;
    const struct Buf2D *const in_what      = &x->xdplane[0].pre[0];
    const MV                  fcenter_mv   = {center_mv->row >> 3, center_mv->col >> 3};
    const uint8_t *           best_address = get_buf_from_mv(in_what, ref_mv);
    unsigned int best_sad = fn_ptr->sdf(what->buf, what->stride, best_address, in_what->stride) +
                            mvsad_err_cost(x, ref_mv, &fcenter_mv, error_per_bit);
    for (int i = 0; i < search_range; i++) {
        int       best_site = -1;
        const int all_in    = (ref_mv->row - 1) > x->mv_limits.row_min &&
            (ref_mv->row + 1) < x->mv_limits.row_max && (ref_mv->col - 1) > x->mv_limits.col_min &&
            (ref_mv->col + 1) < x->mv_limits.col_max;

        if (all_in) {
            unsigned int         sads[4];
            const uint8_t *const positions[4] = {best_address - in_what->stride,
                                                 best_address - 1,
                                                 best_address + 1,
                                                 best_address + in_what->stride};

            fn_ptr->sdx4df(what->buf, what->stride, positions, in_what->stride, sads);

            for (int j = 0; j < 4; ++j) {
                if (sads[j] < best_sad) {
                    const MV mv = {ref_mv->row + neighbors[j].row, ref_mv->col + neighbors[j].col};
                    sads[j] += mvsad_err_cost(x, &mv, &fcenter_mv, error_per_bit);
                    if (sads[j] < best_sad) {
                        best_sad  = sads[j];
                        best_site = j;
                    }
                }
            }
        } else {
            for (int j = 0; j < 4; ++j) {
                const MV mv = {ref_mv->row + neighbors[j].row, ref_mv->col + neighbors[j].col};

                if (is_mv_in(&x->mv_limits, &mv)) {
                    unsigned int sad = fn_ptr->sdf(
                        what->buf, what->stride, get_buf_from_mv(in_what, &mv), in_what->stride);
                    if (sad < best_sad) {
                        sad += mvsad_err_cost(x, &mv, &fcenter_mv, error_per_bit);
                        if (sad < best_sad) {
                            best_sad  = sad;
                            best_site = j;
                        }
                    }
                }
            }
        }

        if (best_site == -1) {
            break;
        } else {
            x->second_best_mv.as_mv = *ref_mv;
            ref_mv->row += neighbors[best_site].row;
            ref_mv->col += neighbors[best_site].col;
            best_address = get_buf_from_mv(in_what, ref_mv);
        }
    }

    return best_sad;
}
static int get_obmc_mvpred_var(const IntraBcContext *x, const int32_t *wsrc, const int32_t *mask,
                               const MV *best_mv, const MV *center_mv,
                               const AomVarianceFnPtr *vfp, int use_mvcost, int is_second) {
    const struct Buf2D *in_what = (const struct Buf2D *)(&x->xdplane[0].pre[is_second]);
    const MV            mv      = {best_mv->row * 8, best_mv->col * 8};
    unsigned int        unused;

    return vfp->ovf(get_buf_from_mv((const struct Buf2D *)in_what, best_mv),
                    in_what->stride,
                    wsrc,
                    mask,
                    &unused) +
           (use_mvcost
                ? mv_err_cost(&mv, center_mv, x->nmv_vec_cost, x->mv_cost_stack, x->errorperbit)
                : 0);
}
static int obmc_refining_search_sad(const IntraBcContext *x, const int32_t *wsrc,
                                    const int32_t *mask, MV *ref_mv, int error_per_bit,
                                    int search_range, const AomVarianceFnPtr *fn_ptr,
                                    const MV *center_mv, int is_second) {
    const MV neighbors[4] = {{-1, 0}, {0, -1}, {0, 1}, {1, 0}};

    const struct Buf2D *in_what    = (const struct Buf2D *)(&x->xdplane[0].pre[is_second]);
    const MV            fcenter_mv = {center_mv->row >> 3, center_mv->col >> 3};
    unsigned int        best_sad =
        fn_ptr->osdf(
            get_buf_from_mv((const struct Buf2D *)in_what, ref_mv), in_what->stride, wsrc, mask) +
        mvsad_err_cost(x, ref_mv, &fcenter_mv, error_per_bit);
    int i, j;

    for (i = 0; i < search_range; i++) {
        int best_site = -1;

        for (j = 0; j < 4; j++) {
            const MV mv = {ref_mv->row + neighbors[j].row, ref_mv->col + neighbors[j].col};
            if (is_mv_in(&x->mv_limits, &mv)) {
                unsigned int sad = fn_ptr->osdf(get_buf_from_mv((const struct Buf2D *)in_what, &mv),
                                                in_what->stride,
                                                wsrc,
                                                mask);
                if (sad < best_sad) {
                    sad += mvsad_err_cost(x, &mv, &fcenter_mv, error_per_bit);
                    if (sad < best_sad) {
                        best_sad  = sad;
                        best_site = j;
                    }
                }
            }
        }

        if (best_site == -1) {
            break;
        } else {
            ref_mv->row += neighbors[best_site].row;
            ref_mv->col += neighbors[best_site].col;
        }
    }
    return best_sad;
}

int svt_av1_obmc_full_pixel_search(ModeDecisionContext *context_ptr, IntraBcContext *x, MV *mvp_full,
                                   int sadpb, const AomVarianceFnPtr *fn_ptr, const MV *ref_mv,
                                   MV *dst_mv, int is_second) {
    // obmc_full_pixel_diamond does not provide BDR gain on 360p
    const int32_t *wsrc         = context_ptr->wsrc_buf;
    const int32_t *mask         = context_ptr->mask_buf;
    const int      search_range = 8;
    *dst_mv                     = *mvp_full;
    clamp_mv(dst_mv,
             x->mv_limits.col_min,
             x->mv_limits.col_max,
             x->mv_limits.row_min,
             x->mv_limits.row_max);
    int thissme = obmc_refining_search_sad(
        x, wsrc, mask, dst_mv, sadpb, search_range, fn_ptr, ref_mv, is_second);
    if (thissme < INT_MAX)
        thissme = get_obmc_mvpred_var(x, wsrc, mask, dst_mv, ref_mv, fn_ptr, 1, is_second);

    return thissme;
}

static INLINE void set_subpel_mv_search_range(const MvLimits *mv_limits, int *col_min, int *col_max,
                                              int *row_min, int *row_max, const MV *ref_mv) {
    const int max_mv = MAX_FULL_PEL_VAL * 8;
    const int minc   = AOMMAX(mv_limits->col_min * 8, ref_mv->col - max_mv);
    const int maxc   = AOMMIN(mv_limits->col_max * 8, ref_mv->col + max_mv);
    const int minr   = AOMMAX(mv_limits->row_min * 8, ref_mv->row - max_mv);
    const int maxr   = AOMMIN(mv_limits->row_max * 8, ref_mv->row + max_mv);

    *col_min = AOMMAX(MV_LOW + 1, minc);
    *col_max = AOMMIN(MV_UPP - 1, maxc);
    *row_min = AOMMAX(MV_LOW + 1, minr);
    *row_max = AOMMIN(MV_UPP - 1, maxr);
}
static const MV search_step_table[12] = {
    // left, right, up, down
    {0, -4},
    {0, 4},
    {-4, 0},
    {4, 0},
    {0, -2},
    {0, 2},
    {-2, 0},
    {2, 0},
    {0, -1},
    {0, 1},
    {-1, 0},
    {1, 0}};

static unsigned int setup_obmc_center_error(const int32_t *mask, const MV *bestmv, const MV *ref_mv,
                                            int error_per_bit, const AomVarianceFnPtr *vfp,
                                            const int32_t *const wsrc, const uint8_t *const y,
                                            int y_stride, int offset, int *mvjcost, int *mvcost[2],
                                            unsigned int *sse1, int *distortion) {
    unsigned int besterr;
    besterr     = vfp->ovf(y + offset, y_stride, wsrc, mask, sse1);
    *distortion = besterr;
    besterr += mv_err_cost(bestmv, ref_mv, mvjcost, mvcost, error_per_bit);
    return besterr;
}

/* checks if (r, c) has better score than previous best */
#define MVC(r, c)                                                                            \
    (unsigned int)(mvcost ? ((mvjcost[((r) != rr) * 2 + ((c) != rc)] + mvcost[0][((r)-rr)] + \
                              (int64_t)mvcost[1][((c)-rc)]) *                                \
                                 error_per_bit +                                             \
                             4096) >>                                                        \
                           13                                                                \
                          : 0)

/* returns subpixel variance error function */
#define DIST(r, c) vfp->osvf(pre(y, y_stride, r, c), y_stride, sp(c), sp(r), z, mask, &sse)
#define CHECK_BETTER(v, r, c)                                                                 \
    do {                                                                                      \
        if (c >= minc && c <= maxc && r >= minr && r <= maxr) {                               \
            thismse = (DIST(r, c));                                                           \
            v       = mv_err_cost(&(const MV){r, c}, ref_mv, mvjcost, mvcost, error_per_bit); \
            if ((v + thismse) < besterr) {                                                    \
                besterr     = v + thismse;                                                    \
                br          = r;                                                              \
                bc          = c;                                                              \
                *distortion = thismse;                                                        \
                *sse1       = sse;                                                            \
            }                                                                                 \
        } else                                                                                \
            v = INT_MAX;                                                                      \
    } while (0)
#define CHECK_BETTER0(v, r, c) CHECK_BETTER(v, r, c)

#define CHECK_BETTER1(v, r, c)                                                          \
    do {                                                                                \
        if (c >= minc && c <= maxc && r >= minr && r <= maxr) {                         \
            MV this_mv = {r, c};                                                        \
            thismse    = upsampled_obmc_pref_error(xd,                                  \
                                                   cm,                                  \
                                                   mi_row,                              \
                                                   mi_col,                              \
                                                   &this_mv,                            \
                                                   mask,                                \
                                                   vfp,                                 \
                                                   z,                                   \
                                                   pre(y, y_stride, r, c),              \
                                                   y_stride,                            \
                                                   sp(c),                               \
                                                   sp(r),                               \
                                                   w,                                   \
                                                   h,                                   \
                                                   &sse,                                \
                                                   use_accurate_subpel_search);         \
            v          = mv_err_cost(&this_mv, ref_mv, mvjcost, mvcost, error_per_bit); \
            if ((v + thismse) < besterr) {                                              \
                besterr     = v + thismse;                                              \
                br          = r;                                                        \
                bc          = c;                                                        \
                *distortion = thismse;                                                  \
                *sse1       = sse;                                                      \
            }                                                                           \
        } else                                                                          \
            v = INT_MAX;                                                                \
    } while (0)

#define SECOND_LEVEL_CHECKS_BEST(k)                      \
    do {                                                 \
        unsigned int second;                             \
        int          br0 = br;                           \
        int          bc0 = bc;                           \
        assert(tr == br || tc == bc);                    \
        if (tr == br && tc != bc)                        \
            kc = bc - tc;                                \
        else if (tr != br && tc == bc)                   \
            kr = br - tr;                                \
        CHECK_BETTER##k(second, br0 + kr, bc0);          \
        CHECK_BETTER##k(second, br0, bc0 + kc);          \
        if (br0 != br || bc0 != bc)                      \
            CHECK_BETTER##k(second, br0 + kr, bc0 + kc); \
    } while (0)

static int upsampled_obmc_pref_error(MacroBlockD *xd, const AV1_COMMON *const cm, int mi_row,
                                     int mi_col, const MV *const mv, const int32_t *mask,
                                     const AomVarianceFnPtr *vfp, const int32_t *const wsrc,
                                     const uint8_t *const y, int y_stride, int subpel_x_q3,
                                     int subpel_y_q3, int w, int h, unsigned int *sse,
                                     int subpel_search) {
    unsigned int besterr;

    DECLARE_ALIGNED(16, uint8_t, pred[2 * MAX_SB_SQUARE]);
#if CONFIG_AV1_HIGHBITDEPTH
    if (is_cur_buf_hbd(xd)) {
        uint8_t *pred8 = CONVERT_TO_BYTEPTR(pred);
        aom_highbd_upsampled_pred(xd,
                                  cm,
                                  mi_row,
                                  mi_col,
                                  mv,
                                  pred8,
                                  w,
                                  h,
                                  subpel_x_q3,
                                  subpel_y_q3,
                                  y,
                                  y_stride,
                                  xd->bd,
                                  subpel_search);
        besterr = vfp->ovf(pred8, w, wsrc, mask, sse);
    } else {
        svt_aom_upsampled_pred(xd,
                               cm,
                               mi_row,
                               mi_col,
                               mv,
                               pred,
                               w,
                               h,
                               subpel_x_q3,
                               subpel_y_q3,
                               y,
                               y_stride,
                               subpel_search);

        besterr = vfp->ovf(pred, w, wsrc, mask, sse);
    }
#else
    svt_aom_upsampled_pred(xd,
                           (const struct AV1Common *const)cm,
                           mi_row,
                           mi_col,
                           mv,
                           pred,
                           w,
                           h,
                           subpel_x_q3,
                           subpel_y_q3,
                           y,
                           y_stride,
                           subpel_search);

    besterr = vfp->ovf(pred, w, wsrc, mask, sse);
#endif
    return besterr;
}
static unsigned int upsampled_setup_obmc_center_error(
    MacroBlockD *xd, const Av1Common *const cm, int mi_row, int mi_col, const int32_t *mask,
    const MV *bestmv, const MV *ref_mv, int error_per_bit, const AomVarianceFnPtr *vfp,
    const int32_t *const wsrc, const uint8_t *const y, int y_stride, int w, int h, int offset,
    int *mvjcost, int *mvcost[2], unsigned int *sse1, int *distortion, int subpel_search) {
    unsigned int besterr = upsampled_obmc_pref_error(xd,
                                                     cm,
                                                     mi_row,
                                                     mi_col,
                                                     bestmv,
                                                     mask,
                                                     vfp,
                                                     wsrc,
                                                     y + offset,
                                                     y_stride,
                                                     0,
                                                     0,
                                                     w,
                                                     h,
                                                     sse1,
                                                     subpel_search);
    *distortion          = besterr;
    besterr += mv_err_cost(bestmv, ref_mv, mvjcost, mvcost, error_per_bit);
    return besterr;
}

// convert motion vector component to offset for sv[a]f calc
static INLINE int   sp(int x) { return x & 7; }
static INLINE const uint8_t *pre(const uint8_t *buf, int stride, int r, int c) {
    const int offset = (r >> 3) * stride + (c >> 3);
    return buf + offset;
}

int svt_av1_find_best_obmc_sub_pixel_tree_up(ModeDecisionContext *context_ptr, IntraBcContext *x,
                                             const AV1_COMMON *const cm, int mi_row, int mi_col,
                                             MV *bestmv, const MV *ref_mv, int allow_hp,
                                             int error_per_bit, const AomVarianceFnPtr *vfp,
                                             int forced_stop, int iters_per_step, int *mvjcost,
                                             int *mvcost[2], int *distortion, unsigned int *sse1,
                                             int is_second, int use_accurate_subpel_search) {
    const int32_t *                wsrc        = context_ptr->wsrc_buf;
    const int32_t *                mask        = context_ptr->mask_buf;
    const int *const               z           = wsrc;
    const int *const               src_address = z;
    MacroBlockD *                  xd          = x->xd;
    struct MacroBlockDPlane *const pd          = &x->xdplane[0];
    unsigned int                   besterr     = INT_MAX;
    unsigned int                   sse;
    unsigned int                   thismse;
    int                            br    = bestmv->row * 8;
    int                            bc    = bestmv->col * 8;
    int                            hstep = 4;
    int                            round = 3 - forced_stop;
    int                            tr;
    int                            tc;
    const MV *                     search_step = search_step_table;
    int                            best_idx    = -1;
    unsigned int                   cost_array[5];
    const int                      w = block_size_wide[context_ptr->blk_geom->bsize];
    const int                      h = block_size_high[context_ptr->blk_geom->bsize];

    int minc, maxc, minr, maxr;

    set_subpel_mv_search_range(&x->mv_limits, &minc, &maxc, &minr, &maxr, ref_mv);

    const uint8_t *y        = pd->pre[is_second].buf;
    int            y_stride = pd->pre[is_second].stride;
    int            offset   = bestmv->row * y_stride + bestmv->col;

    if (!allow_hp && round == 3)
        round = 2;

    bestmv->row *= 8;
    bestmv->col *= 8;
    // use_accurate_subpel_search can be 0 or 1 or 2
    besterr = use_accurate_subpel_search
        ? upsampled_setup_obmc_center_error(xd,
                                            cm,
                                            mi_row,
                                            mi_col,
                                            mask,
                                            bestmv,
                                            ref_mv,
                                            error_per_bit,
                                            vfp,
                                            z,
                                            y,
                                            y_stride,
                                            w,
                                            h,
                                            offset,
                                            mvjcost,
                                            mvcost,
                                            sse1,
                                            distortion,
                                            use_accurate_subpel_search)
        : setup_obmc_center_error(mask,
                                  bestmv,
                                  ref_mv,
                                  error_per_bit,
                                  vfp,
                                  z,
                                  y,
                                  y_stride,
                                  offset,
                                  mvjcost,
                                  mvcost,
                                  sse1,
                                  distortion);

    for (int iter = 0; iter < round; ++iter) {
        // Check vertical and horizontal sub-pixel positions.
        int idx = 0;
        for (; idx < 4; ++idx) {
            tr = br + search_step[idx].row;
            tc = bc + search_step[idx].col;
            if (tc >= minc && tc <= maxc && tr >= minr && tr <= maxr) {
                MV           this_mv = {tr, tc};
                thismse = use_accurate_subpel_search
                    ? (unsigned)upsampled_obmc_pref_error(xd,
                                                cm,
                                                mi_row,
                                                mi_col,
                                                &this_mv,
                                                mask,
                                                vfp,
                                                src_address,
                                                pre(y, y_stride, tr, tc),
                                                y_stride,
                                                sp(tc),
                                                sp(tr),
                                                w,
                                                h,
                                                &sse,
                                                use_accurate_subpel_search)
                    : vfp->osvf(pre(y, y_stride, tr, tc),
                                y_stride,
                                sp(tc),
                                sp(tr),
                                src_address,
                                mask,
                                &sse);

                cost_array[idx] = thismse +
                    mv_err_cost(&this_mv, ref_mv, mvjcost, mvcost, error_per_bit);
                if (cost_array[idx] < besterr) {
                    best_idx    = idx;
                    besterr     = cost_array[idx];
                    *distortion = thismse;
                    *sse1       = sse;
                }
            } else
                cost_array[idx] = INT_MAX;
        }

        // Check diagonal sub-pixel position
        int kc = (cost_array[0] <= cost_array[1] ? -hstep : hstep);
        int kr = (cost_array[2] <= cost_array[3] ? -hstep : hstep);

        tc = bc + kc;
        tr = br + kr;
        if (tc >= minc && tc <= maxc && tr >= minr && tr <= maxr) {
            MV           this_mv = {tr, tc};
            thismse = use_accurate_subpel_search
                ? (unsigned)upsampled_obmc_pref_error(xd,
                                            cm,
                                            mi_row,
                                            mi_col,
                                            &this_mv,
                                            mask,
                                            vfp,
                                            src_address,
                                            pre(y, y_stride, tr, tc),
                                            y_stride,
                                            sp(tc),
                                            sp(tr),
                                            w,
                                            h,
                                            &sse,
                                            use_accurate_subpel_search)
                : vfp->osvf(
                      pre(y, y_stride, tr, tc), y_stride, sp(tc), sp(tr), src_address, mask, &sse);

            cost_array[4] = thismse + mv_err_cost(&this_mv, ref_mv, mvjcost, mvcost, error_per_bit);

            if (cost_array[4] < besterr) {
                best_idx    = 4;
                besterr     = cost_array[4];
                *distortion = thismse;
                *sse1       = sse;
            }
        } else
            cost_array[idx] = INT_MAX;

        if (best_idx < 4 && best_idx >= 0) {
            br += search_step[best_idx].row;
            bc += search_step[best_idx].col;
        } else if (best_idx == 4) {
            br = tr;
            bc = tc;
        }

        if (iters_per_step > 1 && best_idx != -1) {
            if (use_accurate_subpel_search)
                SECOND_LEVEL_CHECKS_BEST(1);
            else
                SECOND_LEVEL_CHECKS_BEST(0);
        }

        search_step += 4;
        hstep >>= 1;
        best_idx = -1;
    }

    bestmv->row = br;
    bestmv->col = bc;

    return besterr;
}

int svt_av1_full_pixel_search(PictureControlSet *pcs, IntraBcContext *x, BlockSize bsize,
                              MV *mvp_full, int step_param, int method, int run_mesh_search,
                              int error_per_bit, int *cost_list, const MV *ref_mv, int var_max,
                              int rd, int x_pos, int y_pos, int intra) {
    UNUSED(run_mesh_search);
    UNUSED(var_max);
    UNUSED(rd);

    int32_t ibc_shift = 0;
    //IBC Modes:   0: OFF 1:Slow   2:Faster   3:Fastest
    if (pcs->parent_pcs_ptr->ibc_mode > 1) ibc_shift = 1;

    SpeedFeatures *sf                   = &pcs->sf;
    sf->exhaustive_searches_thresh      = (1 << 25);
    const AomVarianceFnPtr *fn_ptr = &mefn_ptr[bsize];
    int                          var    = 0;

    if (cost_list) {
        cost_list[0] = INT_MAX;
        cost_list[1] = INT_MAX;
        cost_list[2] = INT_MAX;
        cost_list[3] = INT_MAX;
        cost_list[4] = INT_MAX;
    }

    // Keep track of number of searches (this frame in this thread).
    //if (x->m_search_count_ptr != NULL) ++(*x->m_search_count_ptr);

    switch (method) {
    case FAST_DIAMOND:
        //var = fast_dia_search(x, mvp_full, step_param, error_per_bit, 0,
        //                      cost_list, fn_ptr, 1, ref_mv);
        break;
    case FAST_HEX:
        //var = fast_hex_search(x, mvp_full, step_param, error_per_bit, 0,
        //                      cost_list, fn_ptr, 1, ref_mv);
        break;
    case HEX:
        //var = av1_hex_search(x, mvp_full, step_param, error_per_bit, 1, cost_list,
        //                     fn_ptr, 1, ref_mv);
        break;
    case SQUARE:
        //var = square_search(x, mvp_full, step_param, error_per_bit, 1, cost_list,
        //                    fn_ptr, 1, ref_mv);
        break;
    case BIGDIA:
        //var = bigdia_search(x, mvp_full, step_param, error_per_bit, 1, cost_list,
        //                    fn_ptr, 1, ref_mv);
        break;
    case NSTEP:
        var = full_pixel_diamond(pcs,
                                 x,
                                 mvp_full,
                                 step_param,
                                 error_per_bit,
                                 MAX_MVSEARCH_STEPS - 1 - step_param,
                                 1,
                                 cost_list,
                                 fn_ptr,
                                 ref_mv);

        if (x->is_exhaustive_allowed) {
            int exhuastive_thr = sf->exhaustive_searches_thresh;
            exhuastive_thr >>= 10 - (mi_size_wide_log2[bsize] + mi_size_high_log2[bsize]);

            exhuastive_thr = exhuastive_thr << ibc_shift;

            if (var > exhuastive_thr) {
                int var_ex;
                MV  tmp_mv_ex;
                var_ex = full_pixel_exhaustive(pcs,
                                               x,
                                               &x->best_mv.as_mv,
                                               error_per_bit,
                                               cost_list,
                                               fn_ptr,
                                               ref_mv,
                                               &tmp_mv_ex);

                if (var_ex < var) {
                    var              = var_ex;
                    x->best_mv.as_mv = tmp_mv_ex;
                }
            }
        }
        break;
    default: assert(0 && "Invalid search method.");
    }

    do {
        //CHKN if (!intra || !av1_use_hash_me(&cpi->common)) break;

        // already single ME
        // get block size and original buffer of current block
        const int block_height = block_size_high[bsize];
        const int block_width  = block_size_wide[bsize];
        if (block_height == block_width && x_pos >= 0 && y_pos >= 0) {
            if (block_width == 4 || block_width == 8 || block_width == 16 || block_width == 32 ||
                block_width == 64 || block_width == 128) {
                uint8_t * what        = x->plane[0].src.buf;
                const int what_stride = x->plane[0].src.stride;
                uint32_t  hash_value1, hash_value2;
                MV        best_hash_mv;
                int       best_hash_cost = INT_MAX;

                // for the hashMap
                HashTable *ref_frame_hash = &pcs->hash_table;

                svt_av1_get_block_hash_value(
                    what, what_stride, block_width, &hash_value1, &hash_value2, 0, pcs, x);

                const int count = svt_av1_hash_table_count(ref_frame_hash, hash_value1);
                // for intra, at least one matching can be found, itself.
                if (count <= (intra ? 1 : 0)) break;
                Iterator iterator = svt_av1_hash_get_first_iterator(ref_frame_hash, hash_value1);
                for (int i = 0; i < count; i++, iterator_increment(&iterator)) {
                    BlockHash ref_block_hash = *(BlockHash *)(iterator_get(&iterator));
                    if (hash_value2 == ref_block_hash.hash_value2) {
                        // For intra, make sure the prediction is from valid area.
                        if (intra) {
                            const int mi_col = x_pos / MI_SIZE;
                            const int mi_row = y_pos / MI_SIZE;
                            const MV  dv     = {8 * (ref_block_hash.y - y_pos),
                                           8 * (ref_block_hash.x - x_pos)};
                            if (!av1_is_dv_valid(
                                    dv,
                                    x->xd,
                                    mi_row,
                                    mi_col,
                                    bsize,
                                    pcs->parent_pcs_ptr->scs_ptr->seq_header.sb_size_log2))
                                continue;
                        }
                        MV hash_mv;
                        hash_mv.col = ref_block_hash.x - x_pos;
                        hash_mv.row = ref_block_hash.y - y_pos;
                        if (!is_mv_in(&x->mv_limits, &hash_mv)) continue;
                        const int ref_cost = svt_av1_get_mvpred_var(x, &hash_mv, ref_mv, fn_ptr, 1);
                        if (ref_cost < best_hash_cost) {
                            best_hash_cost = ref_cost;
                            best_hash_mv   = hash_mv;
                        }
                    }
                }

                if (best_hash_cost < var) {
                    x->second_best_mv = x->best_mv;
                    x->best_mv.as_mv  = best_hash_mv;
                    var               = best_hash_cost;
                }
            }
        }
    } while (0);

    return 0; //CHKN  var;
}
