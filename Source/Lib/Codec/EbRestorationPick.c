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

#include <assert.h>
#include <math.h>
#include "EbDefinitions.h"
#include "EbPictureBufferDesc.h"
#include "EbPsnr.h"
#include "EbPictureControlSet.h"
#include "aom_dsp_rtcd.h"
#include "EbRestoration.h"
#include "EbRestorationPick.h"

#if REST_M
#include "EbRestProcess.h"



void av1_foreach_rest_unit_in_frame_seg(Av1Common *cm, int32_t plane,
    rest_tile_start_visitor_t on_tile,
    rest_unit_visitor_t on_rest_unit,
    void *priv,
    PictureControlSet_t   *picture_control_set_ptr,
    uint32_t segment_index);
#endif

void av1_selfguided_restoration_c(const uint8_t *dgd8, int32_t width, int32_t height,
    int32_t dgd_stride, int32_t *flt0, int32_t *flt1,
    int32_t flt_stride, int32_t sgr_params_idx,
    int32_t bit_depth, int32_t highbd);
void av1_foreach_rest_unit_in_frame(Av1Common *cm, int32_t plane,
    rest_tile_start_visitor_t on_tile,
    rest_unit_visitor_t on_rest_unit,
    void *priv);


// When set to RESTORE_WIENER or RESTORE_SGRPROJ only those are allowed.
// When set to RESTORE_TYPES we allow switchable.
static const RestorationType force_restore_type = RESTORE_TYPES;

// Number of Wiener iterations
#define NUM_WIENER_ITERS 5
// Working precision for Wiener filter coefficients
#define WIENER_TAP_SCALE_FACTOR ((int64_t)1 << 16)

typedef int64_t(*sse_part_extractor_type)(const Yv12BufferConfig *a,
    const Yv12BufferConfig *b,
    int32_t hstart, int32_t width, int32_t vstart,
    int32_t height);

#define NUM_EXTRACTORS (3 * (1 + 1))

static const sse_part_extractor_type sse_part_extractors[NUM_EXTRACTORS] = {
  aom_get_y_sse_part,        aom_get_u_sse_part,
  aom_get_v_sse_part,        aom_highbd_get_y_sse_part,
  aom_highbd_get_u_sse_part, aom_highbd_get_v_sse_part,
};
static int64_t sse_restoration_unit(const RestorationTileLimits *limits,
    const Yv12BufferConfig *src,
    const Yv12BufferConfig *dst, int32_t plane,
    int32_t highbd) {
    return sse_part_extractors[3 * highbd + plane](
        src, dst, limits->h_start, limits->h_end - limits->h_start,
        limits->v_start, limits->v_end - limits->v_start);
}

#if ! REST_M
typedef struct {
    // The best coefficients for Wiener or Sgrproj restoration
    WienerInfo wiener;
    SgrprojInfo sgrproj;

    // The sum of squared errors for this rtype.
    int64_t sse[RESTORE_SWITCHABLE_TYPES];

    // The rtype to use for this unit given a frame rtype as
    // index. Indices: WIENER, SGRPROJ, SWITCHABLE.
    RestorationType best_rtype[RESTORE_TYPES - 1];
} RestUnitSearchInfo;
#endif

typedef struct {
    const Yv12BufferConfig *src;
    Yv12BufferConfig *dst;

    Av1Common *cm;
    const Macroblock *x;
    int32_t plane;
    int32_t plane_width;
    int32_t plane_height;
    RestUnitSearchInfo *rusi;
#if REST_M 
    RestUnitSearchInfo *rusi_pic;
    uint32_t  pic_num;
    Yv12BufferConfig * org_frame_to_show;
    int32_t *tmpbuf;
#endif

    uint8_t *dgd_buffer;
    int32_t dgd_stride;
    const uint8_t *src_buffer;
    int32_t src_stride;

    // sse and bits are initialised by reset_rsc in search_rest_type
    int64_t sse;
    int64_t bits;
    int32_t tile_y0, tile_stripe0;

    // sgrproj and wiener are initialised by rsc_on_tile when starting the first
    // tile in the frame.
    SgrprojInfo sgrproj;
    WienerInfo wiener;
} RestSearchCtxt;

static void rsc_on_tile(int32_t tile_row, int32_t tile_col, void *priv) {
    (void)tile_col;

    RestSearchCtxt *rsc = (RestSearchCtxt *)priv;
    set_default_sgrproj(&rsc->sgrproj);
    set_default_wiener(&rsc->wiener);

    rsc->tile_stripe0 =
        (tile_row == 0) ? 0 : rsc->cm->rst_end_stripe[tile_row - 1];
}

static void reset_rsc(RestSearchCtxt *rsc) {
    rsc->sse = 0;
    rsc->bits = 0;
}
#if REST_M
static void init_rsc_seg(
    Yv12BufferConfig *org_fts,
    const Yv12BufferConfig *src, Av1Common *cm,
    const Macroblock *x, int32_t plane, RestUnitSearchInfo *rusi,
    Yv12BufferConfig *dst, RestSearchCtxt *rsc) {
    rsc->src = src;
    rsc->dst = dst;
    rsc->cm = cm;
    rsc->x = x;
    rsc->plane = plane;
    rsc->rusi = rusi;

    rsc->org_frame_to_show = org_fts;

    const Yv12BufferConfig *dgd = org_fts;
    const int32_t is_uv = plane != AOM_PLANE_Y;
    rsc->plane_width = src->crop_widths[is_uv];
    rsc->plane_height = src->crop_heights[is_uv];
    rsc->src_buffer = src->buffers[plane];
    rsc->src_stride = src->strides[is_uv];
    rsc->dgd_buffer = dgd->buffers[plane];
    rsc->dgd_stride = dgd->strides[is_uv];
    assert(src->crop_widths[is_uv] == dgd->crop_widths[is_uv]);
    assert(src->crop_heights[is_uv] == dgd->crop_heights[is_uv]);
}
#endif
static void init_rsc(const Yv12BufferConfig *src, Av1Common *cm,
    const Macroblock *x, int32_t plane, RestUnitSearchInfo *rusi,
    Yv12BufferConfig *dst, RestSearchCtxt *rsc) {
    rsc->src = src;
    rsc->dst = dst;
    rsc->cm = cm;
    rsc->x = x;
    rsc->plane = plane;
    rsc->rusi = rusi;

    const Yv12BufferConfig *dgd = cm->frame_to_show;
    const int32_t is_uv = plane != AOM_PLANE_Y;
    rsc->plane_width = src->crop_widths[is_uv];
    rsc->plane_height = src->crop_heights[is_uv];
    rsc->src_buffer = src->buffers[plane];
    rsc->src_stride = src->strides[is_uv];
    rsc->dgd_buffer = dgd->buffers[plane];
    rsc->dgd_stride = dgd->strides[is_uv];
    assert(src->crop_widths[is_uv] == dgd->crop_widths[is_uv]);
    assert(src->crop_heights[is_uv] == dgd->crop_heights[is_uv]);
}

static int64_t try_restoration_unit(const RestSearchCtxt *rsc,
    const RestorationTileLimits *limits,
    const AV1PixelRect *tile_rect,
    const RestorationUnitInfo *rui) {
    const Av1Common *const cm = rsc->cm;
    const int32_t plane = rsc->plane;
    const int32_t is_uv = plane > 0;
    const RestorationInfo *rsi = &cm->rst_info[plane];
    RestorationLineBuffers rlbs;
    const int32_t bit_depth = cm->bit_depth;
    const int32_t highbd = cm->use_highbitdepth;

    const Yv12BufferConfig *fts = cm->frame_to_show;
    // TODO(yunqing): For now, only use optimized LR filter in decoder. Can be
    // also used in encoder.
    const int32_t optimized_lr = 0;


    av1_loop_restoration_filter_unit(
#if REST_NEED_B
        1,
#endif
        limits, rui, &rsi->boundaries, &rlbs, tile_rect, rsc->tile_stripe0,
        is_uv && cm->subsampling_x, is_uv && cm->subsampling_y, highbd, bit_depth,
        fts->buffers[plane], fts->strides[is_uv], rsc->dst->buffers[plane],
        rsc->dst->strides[is_uv], cm->rst_tmpbuf, optimized_lr);

    return sse_restoration_unit(limits, rsc->src, rsc->dst, plane, highbd);
}
#if REST_M
static int64_t try_restoration_unit_seg(const RestSearchCtxt *rsc,
    const RestorationTileLimits *limits,
    const AV1PixelRect *tile_rect,
    const RestorationUnitInfo *rui) {
    const Av1Common *const cm = rsc->cm;
    const int32_t plane = rsc->plane;
    const int32_t is_uv = plane > 0;
    const RestorationInfo *rsi = &cm->rst_info[plane];
    RestorationLineBuffers rlbs;
    const int32_t bit_depth = cm->bit_depth;
    const int32_t highbd = cm->use_highbitdepth;

    const Yv12BufferConfig *fts = rsc->org_frame_to_show;

    // TODO(yunqing): For now, only use optimized LR filter in decoder. Can be
    // also used in encoder.
    const int32_t optimized_lr = 0;


    av1_loop_restoration_filter_unit(
#if REST_NEED_B
        1,
#endif
        limits, rui, &rsi->boundaries, &rlbs, tile_rect, rsc->tile_stripe0,
        is_uv && cm->subsampling_x, is_uv && cm->subsampling_y, highbd, bit_depth,
        fts->buffers[plane], fts->strides[is_uv], rsc->dst->buffers[plane],
        rsc->dst->strides[is_uv], rsc->tmpbuf, optimized_lr);

    return sse_restoration_unit(limits, rsc->src, rsc->dst, plane, highbd);
}
#endif
int64_t av1_lowbd_pixel_proj_error_c(const uint8_t *src8, int32_t width, int32_t height,
    int32_t src_stride, const uint8_t *dat8,
    int32_t dat_stride, int32_t *flt0,
    int32_t flt0_stride, int32_t *flt1,
    int32_t flt1_stride, int32_t xq[2],
    const sgr_params_type *params) {
    int32_t i, j;
    const uint8_t *src = src8;
    const uint8_t *dat = dat8;
    int64_t err = 0;
    if (params->r[0] > 0 && params->r[1] > 0) {
        for (i = 0; i < height; ++i) {
            for (j = 0; j < width; ++j) {
                assert(flt1[j] < (1 << 15) && flt1[j] > -(1 << 15));
                assert(flt0[j] < (1 << 15) && flt0[j] > -(1 << 15));
                const int32_t u = (int32_t)(dat[j] << SGRPROJ_RST_BITS);
                int32_t v = u << SGRPROJ_PRJ_BITS;
                v += xq[0] * (flt0[j] - u) + xq[1] * (flt1[j] - u);
                const int32_t e =
                    ROUND_POWER_OF_TWO(v, SGRPROJ_RST_BITS + SGRPROJ_PRJ_BITS) - src[j];
                err += e * e;
            }
            dat += dat_stride;
            src += src_stride;
            flt0 += flt0_stride;
            flt1 += flt1_stride;
        }
    }
    else if (params->r[0] > 0) {
        for (i = 0; i < height; ++i) {
            for (j = 0; j < width; ++j) {
                assert(flt0[j] < (1 << 15) && flt0[j] > -(1 << 15));
                const int32_t u = (int32_t)(dat[j] << SGRPROJ_RST_BITS);
                int32_t v = u << SGRPROJ_PRJ_BITS;
                v += xq[0] * (flt0[j] - u);
                const int32_t e =
                    ROUND_POWER_OF_TWO(v, SGRPROJ_RST_BITS + SGRPROJ_PRJ_BITS) - src[j];
                err += e * e;
            }
            dat += dat_stride;
            src += src_stride;
            flt0 += flt0_stride;
        }
    }
    else if (params->r[1] > 0) {
        for (i = 0; i < height; ++i) {
            for (j = 0; j < width; ++j) {
                assert(flt1[j] < (1 << 15) && flt1[j] > -(1 << 15));
                const int32_t u = (int32_t)(dat[j] << SGRPROJ_RST_BITS);
                int32_t v = u << SGRPROJ_PRJ_BITS;
                v += xq[1] * (flt1[j] - u);
                const int32_t e =
                    ROUND_POWER_OF_TWO(v, SGRPROJ_RST_BITS + SGRPROJ_PRJ_BITS) - src[j];
                err += e * e;
            }
            dat += dat_stride;
            src += src_stride;
            flt1 += flt1_stride;
        }
    }
    else {
        for (i = 0; i < height; ++i) {
            for (j = 0; j < width; ++j) {
                const int32_t e = (int32_t)(dat[j]) - src[j];
                err += e * e;
            }
            dat += dat_stride;
            src += src_stride;
        }
    }

    return err;
}

int64_t av1_highbd_pixel_proj_error_c(const uint8_t *src8, int32_t width,
    int32_t height, int32_t src_stride,
    const uint8_t *dat8, int32_t dat_stride,
    int32_t *flt0, int32_t flt0_stride,
    int32_t *flt1, int32_t flt1_stride, int32_t xq[2],
    const sgr_params_type *params) {
    const uint16_t *src = CONVERT_TO_SHORTPTR(src8);
    const uint16_t *dat = CONVERT_TO_SHORTPTR(dat8);
    int32_t i, j;
    int64_t err = 0;
    const int32_t half = 1 << (SGRPROJ_RST_BITS + SGRPROJ_PRJ_BITS - 1);
    if (params->r[0] > 0 && params->r[1] > 0) {
        int32_t xq0 = xq[0];
        int32_t xq1 = xq[1];
        for (i = 0; i < height; ++i) {
            for (j = 0; j < width; ++j) {
                const int32_t d = dat[j];
                const int32_t s = src[j];
                const int32_t u = (int32_t)(d << SGRPROJ_RST_BITS);
                int32_t v0 = flt0[j] - u;
                int32_t v1 = flt1[j] - u;
                int32_t v = half;
                v += xq0 * v0;
                v += xq1 * v1;
                const int32_t e = (v >> (SGRPROJ_RST_BITS + SGRPROJ_PRJ_BITS)) + d - s;
                err += e * e;
            }
            dat += dat_stride;
            flt0 += flt0_stride;
            flt1 += flt1_stride;
            src += src_stride;
        }
    }
    else if (params->r[0] > 0 || params->r[1] > 0) {
        int32_t exq;
        int32_t *flt;
        int32_t flt_stride;
        if (params->r[0] > 0) {
            exq = xq[0];
            flt = flt0;
            flt_stride = flt0_stride;
        }
        else {
            exq = xq[1];
            flt = flt1;
            flt_stride = flt1_stride;
        }
        for (i = 0; i < height; ++i) {
            for (j = 0; j < width; ++j) {
                const int32_t d = dat[j];
                const int32_t s = src[j];
                const int32_t u = (int32_t)(d << SGRPROJ_RST_BITS);
                int32_t v = half;
                v += exq * (flt[j] - u);
                const int32_t e = (v >> (SGRPROJ_RST_BITS + SGRPROJ_PRJ_BITS)) + d - s;
                err += e * e;
            }
            dat += dat_stride;
            flt += flt_stride;
            src += src_stride;
        }
    }
    else {
        for (i = 0; i < height; ++i) {
            for (j = 0; j < width; ++j) {
                const int32_t d = dat[j];
                const int32_t s = src[j];
                const int32_t e = d - s;
                err += e * e;
            }
            dat += dat_stride;
            src += src_stride;
        }
    }
    return err;
}

static int64_t get_pixel_proj_error(const uint8_t *src8, int32_t width, int32_t height,
    int32_t src_stride, const uint8_t *dat8,
    int32_t dat_stride, int32_t use_highbitdepth,
    int32_t *flt0, int32_t flt0_stride,
    int32_t *flt1, int32_t flt1_stride, int32_t *xqd,
    const sgr_params_type *params) {
    int32_t xq[2];
    decode_xq(xqd, xq, params);
    if (!use_highbitdepth) {
        return av1_lowbd_pixel_proj_error(src8, width, height, src_stride, dat8,
            dat_stride, flt0, flt0_stride, flt1,
            flt1_stride, xq, params);
    }
    else {
        return av1_highbd_pixel_proj_error(src8, width, height, src_stride, dat8,
            dat_stride, flt0, flt0_stride, flt1,
            flt1_stride, xq, params);
    }
}

#define USE_SGRPROJ_REFINEMENT_SEARCH 1
static int64_t finer_search_pixel_proj_error(
    const uint8_t *src8, int32_t width, int32_t height, int32_t src_stride,
    const uint8_t *dat8, int32_t dat_stride, int32_t use_highbitdepth, int32_t *flt0,
    int32_t flt0_stride, int32_t *flt1, int32_t flt1_stride, int32_t start_step, int32_t *xqd,
    const sgr_params_type *params) {
    int64_t err = get_pixel_proj_error(
        src8, width, height, src_stride, dat8, dat_stride, use_highbitdepth, flt0,
        flt0_stride, flt1, flt1_stride, xqd, params);
    (void)start_step;
#if USE_SGRPROJ_REFINEMENT_SEARCH
    int64_t err2;
    int32_t tap_min[] = { SGRPROJ_PRJ_MIN0, SGRPROJ_PRJ_MIN1 };
    int32_t tap_max[] = { SGRPROJ_PRJ_MAX0, SGRPROJ_PRJ_MAX1 };
    for (int32_t s = start_step; s >= 1; s >>= 1) {
        for (int32_t p = 0; p < 2; ++p) {
            if ((params->r[0] == 0 && p == 0) || (params->r[1] == 0 && p == 1)) {
                continue;
            }
            int32_t skip = 0;
            do {
                if (xqd[p] - s >= tap_min[p]) {
                    xqd[p] -= s;
                    err2 =
                        get_pixel_proj_error(src8, width, height, src_stride, dat8,
                            dat_stride, use_highbitdepth, flt0,
                            flt0_stride, flt1, flt1_stride, xqd, params);
                    if (err2 > err) {
                        xqd[p] += s;
                    }
                    else {
                        err = err2;
                        skip = 1;
                        // At the highest step size continue moving in the same direction
                        if (s == start_step) continue;
                    }
                }
                break;
            } while (1);
            if (skip) break;
            do {
                if (xqd[p] + s <= tap_max[p]) {
                    xqd[p] += s;
                    err2 =
                        get_pixel_proj_error(src8, width, height, src_stride, dat8,
                            dat_stride, use_highbitdepth, flt0,
                            flt0_stride, flt1, flt1_stride, xqd, params);
                    if (err2 > err) {
                        xqd[p] -= s;
                    }
                    else {
                        err = err2;
                        // At the highest step size continue moving in the same direction
                        if (s == start_step) continue;
                    }
                }
                break;
            } while (1);
        }
    }
#endif  // USE_SGRPROJ_REFINEMENT_SEARCH
    return err;
}

extern void RunEmms();

void get_proj_subspace_c(const uint8_t *src8, int32_t width, int32_t height,
    int32_t src_stride, const uint8_t *dat8,
    int32_t dat_stride, int32_t use_highbitdepth,
    int32_t *flt0, int32_t flt0_stride, int32_t *flt1,
    int32_t flt1_stride, int32_t *xq,
    const sgr_params_type *params) {
    int32_t i, j;
    double H[2][2] = { { 0, 0 }, { 0, 0 } };
    double C[2] = { 0, 0 };
    double Det;
    double x[2];
    const int32_t size = width * height;

    aom_clear_system_state();
    RunEmms();

    // Default
    xq[0] = 0;
    xq[1] = 0;
    if (!use_highbitdepth) {
        const uint8_t *src = src8;
        const uint8_t *dat = dat8;
        for (i = 0; i < height; ++i) {
            for (j = 0; j < width; ++j) {
                const double u = (double)(dat[i * dat_stride + j] << SGRPROJ_RST_BITS);
                const double s =
                    (double)(src[i * src_stride + j] << SGRPROJ_RST_BITS) - u;
                const double f1 =
                    (params->r[0] > 0) ? (double)flt0[i * flt0_stride + j] - u : 0;
                const double f2 =
                    (params->r[1] > 0) ? (double)flt1[i * flt1_stride + j] - u : 0;
                H[0][0] += f1 * f1;
                H[1][1] += f2 * f2;
                H[0][1] += f1 * f2;
                C[0] += f1 * s;
                C[1] += f2 * s;
            }
        }
    }
    else {
        const uint16_t *src = CONVERT_TO_SHORTPTR(src8);
        const uint16_t *dat = CONVERT_TO_SHORTPTR(dat8);
        for (i = 0; i < height; ++i) {
            for (j = 0; j < width; ++j) {
                const double u = (double)(dat[i * dat_stride + j] << SGRPROJ_RST_BITS);
                const double s =
                    (double)(src[i * src_stride + j] << SGRPROJ_RST_BITS) - u;
                const double f1 =
                    (params->r[0] > 0) ? (double)flt0[i * flt0_stride + j] - u : 0;
                const double f2 =
                    (params->r[1] > 0) ? (double)flt1[i * flt1_stride + j] - u : 0;
                H[0][0] += f1 * f1;
                H[1][1] += f2 * f2;
                H[0][1] += f1 * f2;
                C[0] += f1 * s;
                C[1] += f2 * s;
            }
        }
    }
    H[0][0] /= size;
    H[0][1] /= size;
    H[1][1] /= size;
    H[1][0] = H[0][1];
    C[0] /= size;
    C[1] /= size;
    if (params->r[0] == 0) {
        // H matrix is now only the scalar H[1][1]
        // C vector is now only the scalar C[1]
        Det = H[1][1];
        if (Det < 1e-8) return;  // ill-posed, return default values
        x[0] = 0;
        x[1] = C[1] / Det;

        xq[0] = 0;
        xq[1] = (int32_t)rint(x[1] * (1 << SGRPROJ_PRJ_BITS));
    }
    else if (params->r[1] == 0) {
        // H matrix is now only the scalar H[0][0]
        // C vector is now only the scalar C[0]
        Det = H[0][0];
        if (Det < 1e-8) return;  // ill-posed, return default values
        x[0] = C[0] / Det;
        x[1] = 0;

        xq[0] = (int32_t)rint(x[0] * (1 << SGRPROJ_PRJ_BITS));
        xq[1] = 0;
    }
    else {
        Det = (H[0][0] * H[1][1] - H[0][1] * H[1][0]);
        if (Det < 1e-8) return;  // ill-posed, return default values
        x[0] = (H[1][1] * C[0] - H[0][1] * C[1]) / Det;
        x[1] = (H[0][0] * C[1] - H[1][0] * C[0]) / Det;

        xq[0] = (int32_t)rint(x[0] * (1 << SGRPROJ_PRJ_BITS));
        xq[1] = (int32_t)rint(x[1] * (1 << SGRPROJ_PRJ_BITS));
    }
}

void encode_xq(int32_t *xq, int32_t *xqd, const sgr_params_type *params) {
    if (params->r[0] == 0) {
        xqd[0] = 0;
        xqd[1] = clamp((1 << SGRPROJ_PRJ_BITS) - xq[1], SGRPROJ_PRJ_MIN1,
            SGRPROJ_PRJ_MAX1);
    }
    else if (params->r[1] == 0) {
        xqd[0] = clamp(xq[0], SGRPROJ_PRJ_MIN0, SGRPROJ_PRJ_MAX0);
        xqd[1] = clamp((1 << SGRPROJ_PRJ_BITS) - xqd[0], SGRPROJ_PRJ_MIN1,
            SGRPROJ_PRJ_MAX1);
    }
    else {
        xqd[0] = clamp(xq[0], SGRPROJ_PRJ_MIN0, SGRPROJ_PRJ_MAX0);
        xqd[1] = clamp((1 << SGRPROJ_PRJ_BITS) - xqd[0] - xq[1], SGRPROJ_PRJ_MIN1,
            SGRPROJ_PRJ_MAX1);
    }
}

// Apply the self-guided filter across an entire restoration unit.
static void apply_sgr(int32_t sgr_params_idx, const uint8_t *dat8, int32_t width,
    int32_t height, int32_t dat_stride, int32_t use_highbd, int32_t bit_depth,
    int32_t pu_width, int32_t pu_height, int32_t *flt0, int32_t *flt1,
    int32_t flt_stride)
{

    for (int32_t i = 0; i < height; i += pu_height)
    {
        const int32_t h = AOMMIN(pu_height, height - i);
        int32_t *flt0_row = flt0 + i * flt_stride;
        int32_t *flt1_row = flt1 + i * flt_stride;
        const uint8_t *dat8_row = dat8 + i * dat_stride;

        // Iterate over the stripe in blocks of width pu_width
        for (int32_t j = 0; j < width; j += pu_width) {
            const int32_t w = AOMMIN(pu_width, width - j);

            //CHKN SSE
            av1_selfguided_restoration(dat8_row + j, w, h, dat_stride, flt0_row + j,
                flt1_row + j, flt_stride, sgr_params_idx,
                bit_depth, use_highbd);
        }
    }
}

static SgrprojInfo search_selfguided_restoration(
    const uint8_t *dat8, int32_t width, int32_t height, int32_t dat_stride,
    const uint8_t *src8, int32_t src_stride, int32_t use_highbitdepth, int32_t bit_depth,
    int32_t pu_width, int32_t pu_height, int32_t *rstbuf
#if FAST_SG
    , 
    int8_t sg_ref_frame_ep[2],
    int32_t sg_frame_ep_cnt[SGRPROJ_PARAMS],
    int8_t step
#endif
)
{
    int32_t *flt0 = rstbuf;
    int32_t *flt1 = flt0 + RESTORATION_UNITPELS_MAX;
    int32_t ep, bestep = 0;
    int64_t besterr = -1;
    int32_t exqd[2], bestxqd[2] = { 0, 0 };
    int32_t flt_stride = ((width + 7) & ~7) + 8;
    assert(pu_width == (RESTORATION_PROC_UNIT_SIZE >> 1) ||
        pu_width == RESTORATION_PROC_UNIT_SIZE);
    assert(pu_height == (RESTORATION_PROC_UNIT_SIZE >> 1) ||
        pu_height == RESTORATION_PROC_UNIT_SIZE);
#if FAST_SG
    int8_t mid_ep = sg_ref_frame_ep[0] < 0 && sg_ref_frame_ep[1] < 0 ? 0 :
        sg_ref_frame_ep[1] < 0 ? sg_ref_frame_ep[0] :
        sg_ref_frame_ep[0] < 0 ? sg_ref_frame_ep[1] :
        (sg_ref_frame_ep[0] + sg_ref_frame_ep[1]) / 2;

    int8_t start_ep = sg_ref_frame_ep[0] < 0 && sg_ref_frame_ep[1] < 0 ? 0 : AOMMAX(0, mid_ep - step);
    int8_t end_ep = sg_ref_frame_ep[0] < 0 && sg_ref_frame_ep[1] < 0 ? SGRPROJ_PARAMS : AOMMIN(SGRPROJ_PARAMS, mid_ep + step);

    for (ep = start_ep; ep < end_ep; ep++) {
        int32_t exq[2];
        apply_sgr(ep, dat8, width, height, dat_stride, use_highbitdepth, bit_depth,
            pu_width, pu_height, flt0, flt1, flt_stride);
        aom_clear_system_state();
        const sgr_params_type *const params = &sgr_params[ep];
        get_proj_subspace(src8, width, height, src_stride, dat8, dat_stride,
            use_highbitdepth, flt0, flt_stride, flt1, flt_stride, exq,
            params);
        aom_clear_system_state();
        encode_xq(exq, exqd, params);
        int64_t err = finer_search_pixel_proj_error(
            src8, width, height, src_stride, dat8, dat_stride, use_highbitdepth,
            flt0, flt_stride, flt1, flt_stride, 2, exqd, params);
        if (besterr == -1 || err < besterr) {
            bestep = ep;
            besterr = err;
            bestxqd[0] = exqd[0];
            bestxqd[1] = exqd[1];
        }
    }
    sg_frame_ep_cnt[bestep]++;
#else
    for (ep = 0; ep < SGRPROJ_PARAMS; ep++) {
        int32_t exq[2];
        apply_sgr(ep, dat8, width, height, dat_stride, use_highbitdepth, bit_depth,
            pu_width, pu_height, flt0, flt1, flt_stride);
        aom_clear_system_state();
        const sgr_params_type *const params = &sgr_params[ep];
        get_proj_subspace(src8, width, height, src_stride, dat8, dat_stride,
            use_highbitdepth, flt0, flt_stride, flt1, flt_stride, exq,
            params);
        aom_clear_system_state();
        encode_xq(exq, exqd, params);
        int64_t err = finer_search_pixel_proj_error(
            src8, width, height, src_stride, dat8, dat_stride, use_highbitdepth,
            flt0, flt_stride, flt1, flt_stride, 2, exqd, params);
        if (besterr == -1 || err < besterr) {
            bestep = ep;
            besterr = err;
            bestxqd[0] = exqd[0];
            bestxqd[1] = exqd[1];
        }
    }
#endif
    SgrprojInfo ret;
    ret.ep = bestep;
    ret.xqd[0] = bestxqd[0];
    ret.xqd[1] = bestxqd[1];
    return ret;
}
extern int32_t aom_count_primitive_refsubexpfin(uint16_t n, uint16_t k, uint16_t ref, uint16_t v);

static int32_t count_sgrproj_bits(SgrprojInfo *sgrproj_info,
    SgrprojInfo *ref_sgrproj_info) {
    int32_t bits = SGRPROJ_PARAMS_BITS;
    const sgr_params_type *params = &sgr_params[sgrproj_info->ep];
    if (params->r[0] > 0)
        bits += aom_count_primitive_refsubexpfin(
            SGRPROJ_PRJ_MAX0 - SGRPROJ_PRJ_MIN0 + 1, SGRPROJ_PRJ_SUBEXP_K,
            (uint16_t)(ref_sgrproj_info->xqd[0] - SGRPROJ_PRJ_MIN0),
            (uint16_t)(sgrproj_info->xqd[0] - SGRPROJ_PRJ_MIN0));
    if (params->r[1] > 0)
        bits += aom_count_primitive_refsubexpfin(
            SGRPROJ_PRJ_MAX1 - SGRPROJ_PRJ_MIN1 + 1, SGRPROJ_PRJ_SUBEXP_K,
            (uint16_t)(ref_sgrproj_info->xqd[1] - SGRPROJ_PRJ_MIN1),
            (uint16_t)(sgrproj_info->xqd[1] - SGRPROJ_PRJ_MIN1));
    return bits;
}

static void search_sgrproj(const RestorationTileLimits *limits,
    const AV1PixelRect *tile, int32_t rest_unit_idx,
    void *priv) {
    RestSearchCtxt *rsc = (RestSearchCtxt *)priv;
    RestUnitSearchInfo *rusi = &rsc->rusi[rest_unit_idx];

    const Macroblock *const x = rsc->x;
#if FAST_SG
    Av1Common *const cm = rsc->cm;
#else
    const Av1Common *const cm = rsc->cm;
#endif
    const int32_t highbd = cm->use_highbitdepth;
    const int32_t bit_depth = cm->bit_depth;

    uint8_t *dgd_start =
        rsc->dgd_buffer + limits->v_start * rsc->dgd_stride + limits->h_start;
    const uint8_t *src_start =
        rsc->src_buffer + limits->v_start * rsc->src_stride + limits->h_start;

    const int32_t is_uv = rsc->plane > 0;
    const int32_t ss_x = is_uv && cm->subsampling_x;
    const int32_t ss_y = is_uv && cm->subsampling_y;
    const int32_t procunit_width = RESTORATION_PROC_UNIT_SIZE >> ss_x;
    const int32_t procunit_height = RESTORATION_PROC_UNIT_SIZE >> ss_y;

#if FAST_SG
    int8_t step = cm->sg_filter_mode == 1 ? 0 : cm->sg_filter_mode == 2 ? 1 : cm->sg_filter_mode == 3 ? 4 : 16;
#endif

    rusi->sgrproj = search_selfguided_restoration(
        dgd_start, limits->h_end - limits->h_start,
        limits->v_end - limits->v_start, rsc->dgd_stride, src_start,
        rsc->src_stride, highbd, bit_depth, procunit_width, procunit_height,
        cm->rst_tmpbuf
#if FAST_SG
        , cm->sg_ref_frame_ep,
        cm->sg_frame_ep_cnt,
        step
#endif
    );

    RestorationUnitInfo rui;
    rui.restoration_type = RESTORE_SGRPROJ;
    rui.sgrproj_info = rusi->sgrproj;

    rusi->sse[RESTORE_SGRPROJ] = try_restoration_unit(rsc, limits, tile, &rui);

    const int64_t bits_none = x->sgrproj_restore_cost[0];
    const int64_t bits_sgr = x->sgrproj_restore_cost[1] +
        (count_sgrproj_bits(&rusi->sgrproj, &rsc->sgrproj)
            << AV1_PROB_COST_SHIFT);

    double cost_none =
        RDCOST_DBL(x->rdmult, bits_none >> 4, rusi->sse[RESTORE_NONE]);
    double cost_sgr =
        RDCOST_DBL(x->rdmult, bits_sgr >> 4, rusi->sse[RESTORE_SGRPROJ]);

    RestorationType rtype =
        (cost_sgr < cost_none) ? RESTORE_SGRPROJ : RESTORE_NONE;
    rusi->best_rtype[RESTORE_SGRPROJ - 1] = rtype;

    rsc->sse += rusi->sse[rtype];
    rsc->bits += (cost_sgr < cost_none) ? bits_sgr : bits_none;
    if (cost_sgr < cost_none) rsc->sgrproj = rusi->sgrproj;
}

void av1_compute_stats_c(int32_t wiener_win, const uint8_t *dgd, const uint8_t *src,
    int32_t h_start, int32_t h_end, int32_t v_start, int32_t v_end,
    int32_t dgd_stride, int32_t src_stride, int64_t *M,
    int64_t *H) {
    int32_t i, j, k, l;
    int16_t Y[WIENER_WIN2];
    const int32_t wiener_win2 = wiener_win * wiener_win;
    const int32_t wiener_halfwin = (wiener_win >> 1);
    uint8_t avg = find_average(dgd, h_start, h_end, v_start, v_end, dgd_stride);

    memset(M, 0, sizeof(*M) * wiener_win2);
    memset(H, 0, sizeof(*H) * wiener_win2 * wiener_win2);
    for (i = v_start; i < v_end; i++) {
        for (j = h_start; j < h_end; j++) {
            const int16_t X = (int16_t)src[i * src_stride + j] - (int16_t)avg;
            int32_t idx = 0;
            for (k = -wiener_halfwin; k <= wiener_halfwin; k++) {
                for (l = -wiener_halfwin; l <= wiener_halfwin; l++) {
                    Y[idx] = (int16_t)dgd[(i + l) * dgd_stride + (j + k)] - (int16_t)avg;
                    idx++;
                }
            }
            assert(idx == wiener_win2);
            for (k = 0; k < wiener_win2; ++k) {
                M[k] += (int32_t)Y[k] * X;
                for (l = k; l < wiener_win2; ++l) {
                    // H is a symmetric matrix, so we only need to fill out the upper
                    // triangle here. We can copy it down to the lower triangle outside
                    // the (i, j) loops.
                    H[k * wiener_win2 + l] += (int32_t)Y[k] * Y[l];
                }
            }
        }
    }
    for (k = 0; k < wiener_win2; ++k) {
        for (l = k + 1; l < wiener_win2; ++l) {
            H[l * wiener_win2 + k] = H[k * wiener_win2 + l];
        }
    }
}

void av1_compute_stats_highbd_c(int32_t wiener_win, const uint8_t *dgd8,
    const uint8_t *src8, int32_t h_start, int32_t h_end,
    int32_t v_start, int32_t v_end, int32_t dgd_stride,
    int32_t src_stride, int64_t *M, int64_t *H,
    aom_bit_depth_t bit_depth) {
    int32_t i, j, k, l;
    int32_t Y[WIENER_WIN2];
    const int32_t wiener_win2 = wiener_win * wiener_win;
    const int32_t wiener_halfwin = (wiener_win >> 1);
    const uint16_t *src = CONVERT_TO_SHORTPTR(src8);
    const uint16_t *dgd = CONVERT_TO_SHORTPTR(dgd8);
    uint16_t avg =
        find_average_highbd(dgd, h_start, h_end, v_start, v_end, dgd_stride);

    uint8_t bit_depth_divider = 1;
    if (bit_depth == AOM_BITS_12)
        bit_depth_divider = 16;
    else if (bit_depth == AOM_BITS_10)
        bit_depth_divider = 4;

    memset(M, 0, sizeof(*M) * wiener_win2);
    memset(H, 0, sizeof(*H) * wiener_win2 * wiener_win2);
    for (i = v_start; i < v_end; i++) {
        for (j = h_start; j < h_end; j++) {
            const int32_t X = (int32_t)src[i * src_stride + j] - (int32_t)avg;
            int32_t idx = 0;
            for (k = -wiener_halfwin; k <= wiener_halfwin; k++) {
                for (l = -wiener_halfwin; l <= wiener_halfwin; l++) {
                    Y[idx] = (int32_t)dgd[(i + l) * dgd_stride + (j + k)] - (int32_t)avg;
                    idx++;
                }
            }
            assert(idx == wiener_win2);
            for (k = 0; k < wiener_win2; ++k) {
                M[k] += (int64_t)Y[k] * X;
                for (l = k; l < wiener_win2; ++l) {
                    // H is a symmetric matrix, so we only need to fill out the upper
                    // triangle here. We can copy it down to the lower triangle outside
                    // the (i, j) loops.
                    H[k * wiener_win2 + l] += (int64_t)Y[k] * Y[l];
                }
            }
        }
    }
    for (k = 0; k < wiener_win2; ++k) {
        M[k] /= bit_depth_divider;
        H[k * wiener_win2 + k] /= bit_depth_divider;
        for (l = k + 1; l < wiener_win2; ++l) {
            H[k * wiener_win2 + l] /= bit_depth_divider;
            H[l * wiener_win2 + k] = H[k * wiener_win2 + l];
        }
    }
}

static INLINE int32_t wrap_index(int32_t i, int32_t wiener_win) {
    const int32_t wiener_halfwin1 = (wiener_win >> 1) + 1;
    return (i >= wiener_halfwin1 ? wiener_win - 1 - i : i);
}

// Solve linear equations to find Wiener filter tap values
// Taps are output scaled by WIENER_FILT_STEP
static int32_t linsolve_wiener(int32_t n, int64_t *A, int32_t stride, int64_t *b,
    int32_t *x) {
    for (int32_t k = 0; k < n - 1; k++) {
        // Partial pivoting: bring the row with the largest pivot to the top
        for (int32_t i = n - 1; i > k; i--) {
            // If row i has a better (bigger) pivot than row (i-1), swap them
            if (llabs(A[(i - 1) * stride + k]) < llabs(A[i * stride + k])) {
                for (int32_t j = 0; j < n; j++) {
                    const int64_t c = A[i * stride + j];
                    A[i * stride + j] = A[(i - 1) * stride + j];
                    A[(i - 1) * stride + j] = c;
                }
                const int64_t c = b[i];
                b[i] = b[i - 1];
                b[i - 1] = c;
            }
        }
        // Forward elimination (convert A to row-echelon form)
        for (int32_t i = k; i < n - 1; i++) {
            if (A[k * stride + k] == 0) return 0;
            const int64_t c = A[(i + 1) * stride + k];
            const int64_t cd = A[k * stride + k];
            for (int32_t j = 0; j < n; j++) {
                A[(i + 1) * stride + j] -= c / 256 * A[k * stride + j] / cd * 256;
            }
            b[i + 1] -= c * b[k] / cd;
        }
    }
    // Back-substitution
    for (int32_t i = n - 1; i >= 0; i--) {
        if (A[i * stride + i] == 0) return 0;
        int64_t c = 0;
        for (int32_t j = i + 1; j <= n - 1; j++) {
            c += A[i * stride + j] * x[j] / WIENER_TAP_SCALE_FACTOR;
        }
        // Store filter taps x in scaled form.
        x[i] = (int32_t)(WIENER_TAP_SCALE_FACTOR * (b[i] - c) / A[i * stride + i]);
    }

    return 1;
}
// Fix vector b, update vector a
static void update_a_sep_sym(int32_t wiener_win, int64_t **Mc, int64_t **Hc,
    int32_t *a, int32_t *b) {
    int32_t i, j;
    int32_t S[WIENER_WIN];
    int64_t A[WIENER_HALFWIN1], B[WIENER_HALFWIN1 * WIENER_HALFWIN1];
    const int32_t wiener_win2 = wiener_win * wiener_win;
    const int32_t wiener_halfwin1 = (wiener_win >> 1) + 1;
    memset(A, 0, sizeof(A));
    memset(B, 0, sizeof(B));
    for (i = 0; i < wiener_win; i++) {
        for (j = 0; j < wiener_win; ++j) {
            const int32_t jj = wrap_index(j, wiener_win);
            A[jj] += Mc[i][j] * b[i] / WIENER_TAP_SCALE_FACTOR;
        }
    }
    for (i = 0; i < wiener_win; i++) {
        for (j = 0; j < wiener_win; j++) {
            int32_t k, l;
            for (k = 0; k < wiener_win; ++k) {
                for (l = 0; l < wiener_win; ++l) {
                    const int32_t kk = wrap_index(k, wiener_win);
                    const int32_t ll = wrap_index(l, wiener_win);
                    B[ll * wiener_halfwin1 + kk] +=
                        Hc[j * wiener_win + i][k * wiener_win2 + l] * b[i] /
                        WIENER_TAP_SCALE_FACTOR * b[j] / WIENER_TAP_SCALE_FACTOR;
                }
            }
        }
    }
    // Normalization enforcement in the system of equations itself
    for (i = 0; i < wiener_halfwin1 - 1; ++i) {
        A[i] -=
            A[wiener_halfwin1 - 1] * 2 +
            B[i * wiener_halfwin1 + wiener_halfwin1 - 1] -
            2 * B[(wiener_halfwin1 - 1) * wiener_halfwin1 + (wiener_halfwin1 - 1)];
    }
    for (i = 0; i < wiener_halfwin1 - 1; ++i) {
        for (j = 0; j < wiener_halfwin1 - 1; ++j) {
            B[i * wiener_halfwin1 + j] -=
                2 * (B[i * wiener_halfwin1 + (wiener_halfwin1 - 1)] +
                    B[(wiener_halfwin1 - 1) * wiener_halfwin1 + j] -
                    2 * B[(wiener_halfwin1 - 1) * wiener_halfwin1 +
                    (wiener_halfwin1 - 1)]);
        }
    }
    if (linsolve_wiener(wiener_halfwin1 - 1, B, wiener_halfwin1, A, S)) {
        S[wiener_halfwin1 - 1] = WIENER_TAP_SCALE_FACTOR;
        for (i = wiener_halfwin1; i < wiener_win; ++i) {
            S[i] = S[wiener_win - 1 - i];
            S[wiener_halfwin1 - 1] -= 2 * S[i];
        }
        memcpy(a, S, wiener_win * sizeof(*a));
    }
}

// Fix vector a, update vector b
static void update_b_sep_sym(int32_t wiener_win, int64_t **Mc, int64_t **Hc,
    int32_t *a, int32_t *b) {
    int32_t i, j;
    int32_t S[WIENER_WIN];
    int64_t A[WIENER_HALFWIN1], B[WIENER_HALFWIN1 * WIENER_HALFWIN1];
    const int32_t wiener_win2 = wiener_win * wiener_win;
    const int32_t wiener_halfwin1 = (wiener_win >> 1) + 1;
    memset(A, 0, sizeof(A));
    memset(B, 0, sizeof(B));
    for (i = 0; i < wiener_win; i++) {
        const int32_t ii = wrap_index(i, wiener_win);
        for (j = 0; j < wiener_win; j++) {
            A[ii] += Mc[i][j] * a[j] / WIENER_TAP_SCALE_FACTOR;
        }
    }

    for (i = 0; i < wiener_win; i++) {
        for (j = 0; j < wiener_win; j++) {
            const int32_t ii = wrap_index(i, wiener_win);
            const int32_t jj = wrap_index(j, wiener_win);
            int32_t k, l;
            for (k = 0; k < wiener_win; ++k) {
                for (l = 0; l < wiener_win; ++l) {
                    B[jj * wiener_halfwin1 + ii] +=
                        Hc[i * wiener_win + j][k * wiener_win2 + l] * a[k] /
                        WIENER_TAP_SCALE_FACTOR * a[l] / WIENER_TAP_SCALE_FACTOR;
                }
            }
        }
    }
    // Normalization enforcement in the system of equations itself
    for (i = 0; i < wiener_halfwin1 - 1; ++i) {
        A[i] -=
            A[wiener_halfwin1 - 1] * 2 +
            B[i * wiener_halfwin1 + wiener_halfwin1 - 1] -
            2 * B[(wiener_halfwin1 - 1) * wiener_halfwin1 + (wiener_halfwin1 - 1)];
    }
    for (i = 0; i < wiener_halfwin1 - 1; ++i) {
        for (j = 0; j < wiener_halfwin1 - 1; ++j) {
            B[i * wiener_halfwin1 + j] -=
                2 * (B[i * wiener_halfwin1 + (wiener_halfwin1 - 1)] +
                    B[(wiener_halfwin1 - 1) * wiener_halfwin1 + j] -
                    2 * B[(wiener_halfwin1 - 1) * wiener_halfwin1 +
                    (wiener_halfwin1 - 1)]);
        }
    }
    if (linsolve_wiener(wiener_halfwin1 - 1, B, wiener_halfwin1, A, S)) {
        S[wiener_halfwin1 - 1] = WIENER_TAP_SCALE_FACTOR;
        for (i = wiener_halfwin1; i < wiener_win; ++i) {
            S[i] = S[wiener_win - 1 - i];
            S[wiener_halfwin1 - 1] -= 2 * S[i];
        }
        memcpy(b, S, wiener_win * sizeof(*b));
    }
}

static int32_t wiener_decompose_sep_sym(int32_t wiener_win, int64_t *M, int64_t *H,
    int32_t *a, int32_t *b) {
    static const int32_t init_filt[WIENER_WIN] = {
        WIENER_FILT_TAP0_MIDV, WIENER_FILT_TAP1_MIDV, WIENER_FILT_TAP2_MIDV,
        WIENER_FILT_TAP3_MIDV, WIENER_FILT_TAP2_MIDV, WIENER_FILT_TAP1_MIDV,
        WIENER_FILT_TAP0_MIDV,
    };
    int64_t *Hc[WIENER_WIN2];
    int64_t *Mc[WIENER_WIN];
    int32_t i, j, iter;
    const int32_t plane_off = (WIENER_WIN - wiener_win) >> 1;
    const int32_t wiener_win2 = wiener_win * wiener_win;
    for (i = 0; i < wiener_win; i++) {
        a[i] = b[i] =
            WIENER_TAP_SCALE_FACTOR / WIENER_FILT_STEP * init_filt[i + plane_off];
    }
    for (i = 0; i < wiener_win; i++) {
        Mc[i] = M + i * wiener_win;
        for (j = 0; j < wiener_win; j++) {
            Hc[i * wiener_win + j] =
                H + i * wiener_win * wiener_win2 + j * wiener_win;
        }
    }

    iter = 1;
    while (iter < NUM_WIENER_ITERS) {
        update_a_sep_sym(wiener_win, Mc, Hc, a, b);
        update_b_sep_sym(wiener_win, Mc, Hc, a, b);
        iter++;
    }
    return 1;
}
static int64_t compute_score(int32_t wiener_win, int64_t *M, int64_t *H,
    InterpKernel vfilt, InterpKernel hfilt) {
    int32_t ab[WIENER_WIN * WIENER_WIN];
    int16_t a[WIENER_WIN], b[WIENER_WIN];
    int64_t P = 0, Q = 0;
    int64_t iP = 0, iQ = 0;
    int64_t Score, iScore;
    int32_t i, k, l;
    const int32_t plane_off = (WIENER_WIN - wiener_win) >> 1;
    const int32_t wiener_win2 = wiener_win * wiener_win;

    aom_clear_system_state();

    a[WIENER_HALFWIN] = b[WIENER_HALFWIN] = WIENER_FILT_STEP;
    for (i = 0; i < WIENER_HALFWIN; ++i) {
        a[i] = a[WIENER_WIN - i - 1] = vfilt[i];
        b[i] = b[WIENER_WIN - i - 1] = hfilt[i];
        a[WIENER_HALFWIN] -= 2 * a[i];
        b[WIENER_HALFWIN] -= 2 * b[i];
    }
    memset(ab, 0, sizeof(ab));
    for (k = 0; k < wiener_win; ++k) {
        for (l = 0; l < wiener_win; ++l)
            ab[k * wiener_win + l] = a[l + plane_off] * b[k + plane_off];
    }
    for (k = 0; k < wiener_win2; ++k) {
        P += ab[k] * M[k] / WIENER_FILT_STEP / WIENER_FILT_STEP;
        for (l = 0; l < wiener_win2; ++l) {
            Q += ab[k] * H[k * wiener_win2 + l] * ab[l] / WIENER_FILT_STEP /
                WIENER_FILT_STEP / WIENER_FILT_STEP / WIENER_FILT_STEP;
        }
    }
    Score = Q - 2 * P;

    iP = M[wiener_win2 >> 1];
    iQ = H[(wiener_win2 >> 1) * wiener_win2 + (wiener_win2 >> 1)];
    iScore = iQ - 2 * iP;

    return Score - iScore;
}

static void finalize_sym_filter(int32_t wiener_win, int32_t *f, InterpKernel fi) {
    int32_t i;
    const int32_t wiener_halfwin = (wiener_win >> 1);

    for (i = 0; i < wiener_halfwin; ++i) {
        const int64_t dividend = f[i] * WIENER_FILT_STEP;
        const int64_t divisor = WIENER_TAP_SCALE_FACTOR;
        // Perform this division with proper rounding rather than truncation
        if (dividend < 0) {
            fi[i] = (int16_t)((dividend - (divisor / 2)) / divisor);
        }
        else {
            fi[i] = (int16_t)((dividend + (divisor / 2)) / divisor);
        }
    }
    // Specialize for 7-tap filter
    if (wiener_win == WIENER_WIN) {
        fi[0] = CLIP(fi[0], WIENER_FILT_TAP0_MINV, WIENER_FILT_TAP0_MAXV);
        fi[1] = CLIP(fi[1], WIENER_FILT_TAP1_MINV, WIENER_FILT_TAP1_MAXV);
        fi[2] = CLIP(fi[2], WIENER_FILT_TAP2_MINV, WIENER_FILT_TAP2_MAXV);
    }
    else {
        fi[2] = CLIP(fi[1], WIENER_FILT_TAP2_MINV, WIENER_FILT_TAP2_MAXV);
        fi[1] = CLIP(fi[0], WIENER_FILT_TAP1_MINV, WIENER_FILT_TAP1_MAXV);
        fi[0] = 0;
    }
    // Satisfy filter constraints
    fi[WIENER_WIN - 1] = fi[0];
    fi[WIENER_WIN - 2] = fi[1];
    fi[WIENER_WIN - 3] = fi[2];
    // The central element has an implicit +WIENER_FILT_STEP
    fi[3] = -2 * (fi[0] + fi[1] + fi[2]);
}


static int32_t count_wiener_bits(int32_t wiener_win, WienerInfo *wiener_info,
    WienerInfo *ref_wiener_info) {
    int32_t bits = 0;
    if (wiener_win == WIENER_WIN)
        bits += aom_count_primitive_refsubexpfin(
            WIENER_FILT_TAP0_MAXV - WIENER_FILT_TAP0_MINV + 1,
            WIENER_FILT_TAP0_SUBEXP_K,
            ref_wiener_info->vfilter[0] - WIENER_FILT_TAP0_MINV,
            wiener_info->vfilter[0] - WIENER_FILT_TAP0_MINV);
    bits += aom_count_primitive_refsubexpfin(
        WIENER_FILT_TAP1_MAXV - WIENER_FILT_TAP1_MINV + 1,
        WIENER_FILT_TAP1_SUBEXP_K,
        ref_wiener_info->vfilter[1] - WIENER_FILT_TAP1_MINV,
        wiener_info->vfilter[1] - WIENER_FILT_TAP1_MINV);
    bits += aom_count_primitive_refsubexpfin(
        WIENER_FILT_TAP2_MAXV - WIENER_FILT_TAP2_MINV + 1,
        WIENER_FILT_TAP2_SUBEXP_K,
        ref_wiener_info->vfilter[2] - WIENER_FILT_TAP2_MINV,
        wiener_info->vfilter[2] - WIENER_FILT_TAP2_MINV);
    if (wiener_win == WIENER_WIN)
        bits += aom_count_primitive_refsubexpfin(
            WIENER_FILT_TAP0_MAXV - WIENER_FILT_TAP0_MINV + 1,
            WIENER_FILT_TAP0_SUBEXP_K,
            ref_wiener_info->hfilter[0] - WIENER_FILT_TAP0_MINV,
            wiener_info->hfilter[0] - WIENER_FILT_TAP0_MINV);
    bits += aom_count_primitive_refsubexpfin(
        WIENER_FILT_TAP1_MAXV - WIENER_FILT_TAP1_MINV + 1,
        WIENER_FILT_TAP1_SUBEXP_K,
        ref_wiener_info->hfilter[1] - WIENER_FILT_TAP1_MINV,
        wiener_info->hfilter[1] - WIENER_FILT_TAP1_MINV);
    bits += aom_count_primitive_refsubexpfin(
        WIENER_FILT_TAP2_MAXV - WIENER_FILT_TAP2_MINV + 1,
        WIENER_FILT_TAP2_SUBEXP_K,
        ref_wiener_info->hfilter[2] - WIENER_FILT_TAP2_MINV,
        wiener_info->hfilter[2] - WIENER_FILT_TAP2_MINV);
    return bits;
}

#define USE_WIENER_REFINEMENT_SEARCH 1
static int64_t finer_tile_search_wiener(const RestSearchCtxt *rsc,
    const RestorationTileLimits *limits,
    const AV1PixelRect *tile,
    RestorationUnitInfo *rui,
    int32_t wiener_win) {
    const int32_t plane_off = (WIENER_WIN - wiener_win) >> 1;
    int64_t err = try_restoration_unit(rsc, limits, tile, rui);
#if USE_WIENER_REFINEMENT_SEARCH
    int64_t err2;
    int32_t tap_min[] = { WIENER_FILT_TAP0_MINV, WIENER_FILT_TAP1_MINV,
                      WIENER_FILT_TAP2_MINV };
    int32_t tap_max[] = { WIENER_FILT_TAP0_MAXV, WIENER_FILT_TAP1_MAXV,
                      WIENER_FILT_TAP2_MAXV };

    WienerInfo *plane_wiener = &rui->wiener_info;

    // printf("err  pre = %"PRId64"\n", err);
    const int32_t start_step = 4;
    for (int32_t s = start_step; s >= 1; s >>= 1) {
        for (int32_t p = plane_off; p < WIENER_HALFWIN; ++p) {
            int32_t skip = 0;
            do {
                if (plane_wiener->hfilter[p] - s >= tap_min[p]) {
                    plane_wiener->hfilter[p] -= (int16_t)s;
                    plane_wiener->hfilter[WIENER_WIN - p - 1] -= (int16_t)s;
                    plane_wiener->hfilter[WIENER_HALFWIN] += 2 * (int16_t)s;
                    err2 = try_restoration_unit(rsc, limits, tile, rui);
                    if (err2 > err) {
                        plane_wiener->hfilter[p] += (int16_t)s;
                        plane_wiener->hfilter[WIENER_WIN - p - 1] += (int16_t)s;
                        plane_wiener->hfilter[WIENER_HALFWIN] -= 2 * (int16_t)s;
                    }
                    else {
                        err = err2;
                        skip = 1;
                        // At the highest step size continue moving in the same direction
                        if (s == start_step) continue;
                    }
                }
                break;
            } while (1);
            if (skip) break;
            do {
                if (plane_wiener->hfilter[p] + s <= tap_max[p]) {
                    plane_wiener->hfilter[p] += (int16_t)s;
                    plane_wiener->hfilter[WIENER_WIN - p - 1] += (int16_t)s;
                    plane_wiener->hfilter[WIENER_HALFWIN] -= 2 * (int16_t)s;
                    err2 = try_restoration_unit(rsc, limits, tile, rui);
                    if (err2 > err) {
                        plane_wiener->hfilter[p] -= (int16_t)s;
                        plane_wiener->hfilter[WIENER_WIN - p - 1] -= (int16_t)s;
                        plane_wiener->hfilter[WIENER_HALFWIN] += 2 * (int16_t)s;
                    }
                    else {
                        err = err2;
                        // At the highest step size continue moving in the same direction
                        if (s == start_step) continue;
                    }
                }
                break;
            } while (1);
        }
        for (int32_t p = plane_off; p < WIENER_HALFWIN; ++p) {
            int32_t skip = 0;
            do {
                if (plane_wiener->vfilter[p] - s >= tap_min[p]) {
                    plane_wiener->vfilter[p] -= (int16_t)s;
                    plane_wiener->vfilter[WIENER_WIN - p - 1] -= (int16_t)s;
                    plane_wiener->vfilter[WIENER_HALFWIN] += 2 * (int16_t)s;
                    err2 = try_restoration_unit(rsc, limits, tile, rui);
                    if (err2 > err) {
                        plane_wiener->vfilter[p] += (int16_t)s;
                        plane_wiener->vfilter[WIENER_WIN - p - 1] += (int16_t)s;
                        plane_wiener->vfilter[WIENER_HALFWIN] -= 2 * (int16_t)s;
                    }
                    else {
                        err = err2;
                        skip = 1;
                        // At the highest step size continue moving in the same direction
                        if (s == start_step) continue;
                    }
                }
                break;
            } while (1);
            if (skip) break;
            do {
                if (plane_wiener->vfilter[p] + s <= tap_max[p]) {
                    plane_wiener->vfilter[p] += (int16_t)s;
                    plane_wiener->vfilter[WIENER_WIN - p - 1] += (int16_t)s;
                    plane_wiener->vfilter[WIENER_HALFWIN] -= 2 * (int16_t)s;
                    err2 = try_restoration_unit(rsc, limits, tile, rui);
                    if (err2 > err) {
                        plane_wiener->vfilter[p] -= (int16_t)s;
                        plane_wiener->vfilter[WIENER_WIN - p - 1] -= (int16_t)s;
                        plane_wiener->vfilter[WIENER_HALFWIN] += 2 * (int16_t)s;
                    }
                    else {
                        err = err2;
                        // At the highest step size continue moving in the same direction
                        if (s == start_step) continue;
                    }
                }
                break;
            } while (1);
        }
    }
    // printf("err post = %"PRId64"\n", err);
#endif  // USE_WIENER_REFINEMENT_SEARCH
    return err;
}
#if REST_M
static int64_t finer_tile_search_wiener_seg(const RestSearchCtxt *rsc,
    const RestorationTileLimits *limits,
    const AV1PixelRect *tile,
    RestorationUnitInfo *rui,
    int32_t wiener_win) {
    const int32_t plane_off = (WIENER_WIN - wiener_win) >> 1;
    int64_t err = try_restoration_unit_seg(rsc, limits, tile, rui);
#if USE_WIENER_REFINEMENT_SEARCH
    int64_t err2;
    int32_t tap_min[] = { WIENER_FILT_TAP0_MINV, WIENER_FILT_TAP1_MINV,
                      WIENER_FILT_TAP2_MINV };
    int32_t tap_max[] = { WIENER_FILT_TAP0_MAXV, WIENER_FILT_TAP1_MAXV,
                      WIENER_FILT_TAP2_MAXV };

    WienerInfo *plane_wiener = &rui->wiener_info;

    // printf("err  pre = %"PRId64"\n", err);
    const int32_t start_step = 4;
    for (int32_t s = start_step; s >= 1; s >>= 1) {
        for (int32_t p = plane_off; p < WIENER_HALFWIN; ++p) {
            int32_t skip = 0;
            do {
                if (plane_wiener->hfilter[p] - s >= tap_min[p]) {
                    plane_wiener->hfilter[p] -= (int16_t)s;
                    plane_wiener->hfilter[WIENER_WIN - p - 1] -= (int16_t)s;
                    plane_wiener->hfilter[WIENER_HALFWIN] += 2 * (int16_t)s;
                    err2 = try_restoration_unit_seg(rsc, limits, tile, rui);
                    if (err2 > err) {
                        plane_wiener->hfilter[p] += (int16_t)s;
                        plane_wiener->hfilter[WIENER_WIN - p - 1] += (int16_t)s;
                        plane_wiener->hfilter[WIENER_HALFWIN] -= 2 * (int16_t)s;
                    }
                    else {
                        err = err2;
                        skip = 1;
                        // At the highest step size continue moving in the same direction
                        if (s == start_step) continue;
                    }
                }
                break;
            } while (1);
            if (skip) break;
            do {
                if (plane_wiener->hfilter[p] + s <= tap_max[p]) {
                    plane_wiener->hfilter[p] += (int16_t)s;
                    plane_wiener->hfilter[WIENER_WIN - p - 1] += (int16_t)s;
                    plane_wiener->hfilter[WIENER_HALFWIN] -= 2 * (int16_t)s;
                    err2 = try_restoration_unit_seg(rsc, limits, tile, rui);
                    if (err2 > err) {
                        plane_wiener->hfilter[p] -= (int16_t)s;
                        plane_wiener->hfilter[WIENER_WIN - p - 1] -= (int16_t)s;
                        plane_wiener->hfilter[WIENER_HALFWIN] += 2 * (int16_t)s;
                    }
                    else {
                        err = err2;
                        // At the highest step size continue moving in the same direction
                        if (s == start_step) continue;
                    }
                }
                break;
            } while (1);
        }
        for (int32_t p = plane_off; p < WIENER_HALFWIN; ++p) {
            int32_t skip = 0;
            do {
                if (plane_wiener->vfilter[p] - s >= tap_min[p]) {
                    plane_wiener->vfilter[p] -= (int16_t)s;
                    plane_wiener->vfilter[WIENER_WIN - p - 1] -= (int16_t)s;
                    plane_wiener->vfilter[WIENER_HALFWIN] += 2 * (int16_t)s;
                    err2 = try_restoration_unit_seg(rsc, limits, tile, rui);
                    if (err2 > err) {
                        plane_wiener->vfilter[p] += (int16_t)s;
                        plane_wiener->vfilter[WIENER_WIN - p - 1] += (int16_t)s;
                        plane_wiener->vfilter[WIENER_HALFWIN] -= 2 * (int16_t)s;
                    }
                    else {
                        err = err2;
                        skip = 1;
                        // At the highest step size continue moving in the same direction
                        if (s == start_step) continue;
                    }
                }
                break;
            } while (1);
            if (skip) break;
            do {
                if (plane_wiener->vfilter[p] + s <= tap_max[p]) {
                    plane_wiener->vfilter[p] += (int16_t)s;
                    plane_wiener->vfilter[WIENER_WIN - p - 1] += (int16_t)s;
                    plane_wiener->vfilter[WIENER_HALFWIN] -= 2 * (int16_t)s;
                    err2 = try_restoration_unit_seg(rsc, limits, tile, rui);
                    if (err2 > err) {
                        plane_wiener->vfilter[p] -= (int16_t)s;
                        plane_wiener->vfilter[WIENER_WIN - p - 1] -= (int16_t)s;
                        plane_wiener->vfilter[WIENER_HALFWIN] += 2 * (int16_t)s;
                    }
                    else {
                        err = err2;
                        // At the highest step size continue moving in the same direction
                        if (s == start_step) continue;
                    }
                }
                break;
            } while (1);
        }
    }
    // printf("err post = %"PRId64"\n", err);
#endif  // USE_WIENER_REFINEMENT_SEARCH
    return err;
}
#endif
static void search_wiener(const RestorationTileLimits *limits,
    const AV1PixelRect *tile_rect, int32_t rest_unit_idx,
    void *priv) {
    RestSearchCtxt *rsc = (RestSearchCtxt *)priv;
    RestUnitSearchInfo *rusi = &rsc->rusi[rest_unit_idx];
    const Av1Common *const cm = rsc->cm;
#if FAST_WN
    int32_t wn_luma = cm->wn_filter_mode == 1 ? WIENER_WIN_CHROMA : WIENER_WIN;
    const int32_t wiener_win =
        (rsc->plane == AOM_PLANE_Y) ? wn_luma : WIENER_WIN_CHROMA;
#else
    const int32_t wiener_win =
        (rsc->plane == AOM_PLANE_Y) ? WIENER_WIN : WIENER_WIN_CHROMA;
#endif
    int64_t M[WIENER_WIN2];
    int64_t H[WIENER_WIN2 * WIENER_WIN2];
    int32_t vfilterd[WIENER_WIN], hfilterd[WIENER_WIN];


    if (cm->use_highbitdepth)
    {
        if (rsc->plane == AOM_PLANE_Y) {

            av1_compute_stats_highbd(WIENER_WIN, rsc->dgd_buffer, rsc->src_buffer,

                limits->h_start, limits->h_end, limits->v_start,
                limits->v_end, rsc->dgd_stride, rsc->src_stride, M,
                H, (aom_bit_depth_t)cm->bit_depth);

        }
        else {

            av1_compute_stats_highbd(WIENER_WIN_CHROMA, rsc->dgd_buffer, rsc->src_buffer,

                limits->h_start, limits->h_end, limits->v_start,
                limits->v_end, rsc->dgd_stride, rsc->src_stride, M,
                H, (aom_bit_depth_t)cm->bit_depth);

        }
    }
    else {

        av1_compute_stats(wiener_win, rsc->dgd_buffer, rsc->src_buffer, limits->h_start,
            limits->h_end, limits->v_start, limits->v_end,
            rsc->dgd_stride, rsc->src_stride, M, H);

    }

    const Macroblock *const x = rsc->x;
    const int64_t bits_none = x->wiener_restore_cost[0];

    if (!wiener_decompose_sep_sym(wiener_win, M, H, vfilterd, hfilterd)) {
        rsc->bits += bits_none;
        rsc->sse += rusi->sse[RESTORE_NONE];
        rusi->best_rtype[RESTORE_WIENER - 1] = RESTORE_NONE;
        rusi->sse[RESTORE_WIENER] = INT64_MAX;
        return;
    }

    RestorationUnitInfo rui;
    memset(&rui, 0, sizeof(rui));
    rui.restoration_type = RESTORE_WIENER;
    finalize_sym_filter(wiener_win, vfilterd, rui.wiener_info.vfilter);
    finalize_sym_filter(wiener_win, hfilterd, rui.wiener_info.hfilter);

    // Filter score computes the value of the function x'*A*x - x'*b for the
    // learned filter and compares it against identity filer. If there is no
    // reduction in the function, the filter is reverted back to identity
    if (compute_score(wiener_win, M, H, rui.wiener_info.vfilter,
        rui.wiener_info.hfilter) > 0) {
        rsc->bits += bits_none;
        rsc->sse += rusi->sse[RESTORE_NONE];
        rusi->best_rtype[RESTORE_WIENER - 1] = RESTORE_NONE;
        rusi->sse[RESTORE_WIENER] = INT64_MAX;
        return;
    }

    aom_clear_system_state();

    rusi->sse[RESTORE_WIENER] =
        finer_tile_search_wiener(rsc, limits, tile_rect, &rui, wiener_win);
    rusi->wiener = rui.wiener_info;

    if (wiener_win != WIENER_WIN) {
        assert(rui.wiener_info.vfilter[0] == 0 &&
            rui.wiener_info.vfilter[WIENER_WIN - 1] == 0);
        assert(rui.wiener_info.hfilter[0] == 0 &&
            rui.wiener_info.hfilter[WIENER_WIN - 1] == 0);
    }

    const int64_t bits_wiener =
        x->wiener_restore_cost[1] +
        (count_wiener_bits(wiener_win, &rusi->wiener, &rsc->wiener)
            << AV1_PROB_COST_SHIFT);

    double cost_none =
        RDCOST_DBL(x->rdmult, bits_none >> 4, rusi->sse[RESTORE_NONE]);
    double cost_wiener =
        RDCOST_DBL(x->rdmult, bits_wiener >> 4, rusi->sse[RESTORE_WIENER]);

    RestorationType rtype =
        (cost_wiener < cost_none) ? RESTORE_WIENER : RESTORE_NONE;
    rusi->best_rtype[RESTORE_WIENER - 1] = rtype;

    rsc->sse += rusi->sse[rtype];
    rsc->bits += (cost_wiener < cost_none) ? bits_wiener : bits_none;
    if (cost_wiener < cost_none) rsc->wiener = rusi->wiener;
}

static void search_norestore(const RestorationTileLimits *limits,
    const AV1PixelRect *tile_rect, int32_t rest_unit_idx,
    void *priv) {
    (void)tile_rect;

    RestSearchCtxt *rsc = (RestSearchCtxt *)priv;
    RestUnitSearchInfo *rusi = &rsc->rusi[rest_unit_idx];

    const int32_t highbd = rsc->cm->use_highbitdepth;
    rusi->sse[RESTORE_NONE] = sse_restoration_unit(
        limits, rsc->src, rsc->cm->frame_to_show, rsc->plane, highbd);

    rsc->sse += rusi->sse[RESTORE_NONE];
}

static void search_switchable(const RestorationTileLimits *limits,
    const AV1PixelRect *tile_rect, int32_t rest_unit_idx,
    void *priv) {
    (void)limits;
    (void)tile_rect;
    RestSearchCtxt *rsc = (RestSearchCtxt *)priv;
    RestUnitSearchInfo *rusi = &rsc->rusi[rest_unit_idx];

    const Macroblock *const x = rsc->x;

    const int32_t wiener_win =
        (rsc->plane == AOM_PLANE_Y) ? WIENER_WIN : WIENER_WIN_CHROMA;

    double best_cost = 0;
    int64_t best_bits = 0;
    RestorationType best_rtype = RESTORE_NONE;

    //CHKN for (RestorationType r = 0; r < RESTORE_SWITCHABLE_TYPES; ++r) {
    for (int32_t restType = 0; restType < RESTORE_SWITCHABLE_TYPES; ++restType) {


        RestorationType r = (RestorationType)restType;


        // Check for the condition that wiener or sgrproj search could not
        // find a solution or the solution was worse than RESTORE_NONE.
        // In either case the best_rtype will be set as RESTORE_NONE. These
        // should be skipped from the test below.
        if (r > RESTORE_NONE) {
            if (rusi->best_rtype[r - 1] == RESTORE_NONE) continue;
        }

        const int64_t sse = rusi->sse[r];
        int64_t coeff_pcost = 0;
        switch (r) {
        case RESTORE_NONE: coeff_pcost = 0; break;
        case RESTORE_WIENER:
            coeff_pcost =
                count_wiener_bits(wiener_win, &rusi->wiener, &rsc->wiener);
            break;
        case RESTORE_SGRPROJ:
            coeff_pcost = count_sgrproj_bits(&rusi->sgrproj, &rsc->sgrproj);
            break;
        default: assert(0); break;
        }
        const int64_t coeff_bits = coeff_pcost << AV1_PROB_COST_SHIFT;
        const int64_t bits = x->switchable_restore_cost[r] + coeff_bits;
        double cost = RDCOST_DBL(x->rdmult, bits >> 4, sse);
        if (r == 0 || cost < best_cost) {
            best_cost = cost;
            best_bits = bits;
            best_rtype = r;
        }
    }

    rusi->best_rtype[RESTORE_SWITCHABLE - 1] = best_rtype;

    rsc->sse += rusi->sse[best_rtype];
    rsc->bits += best_bits;
    if (best_rtype == RESTORE_WIENER) rsc->wiener = rusi->wiener;
    if (best_rtype == RESTORE_SGRPROJ) rsc->sgrproj = rusi->sgrproj;
}

static void copy_unit_info(RestorationType frame_rtype,
    const RestUnitSearchInfo *rusi,
    RestorationUnitInfo *rui) {
    assert(frame_rtype > 0);
    rui->restoration_type = rusi->best_rtype[frame_rtype - 1];
    if (rui->restoration_type == RESTORE_WIENER)
        rui->wiener_info = rusi->wiener;
    else
        rui->sgrproj_info = rusi->sgrproj;
}

static double search_rest_type(RestSearchCtxt *rsc, RestorationType rtype)
{

    static const rest_unit_visitor_t funs[RESTORE_TYPES] = {
    search_norestore, search_wiener, search_sgrproj, search_switchable
    };

    reset_rsc(rsc);

    av1_foreach_rest_unit_in_frame(rsc->cm, rsc->plane, rsc_on_tile, funs[rtype], rsc);

    return RDCOST_DBL(rsc->x->rdmult, rsc->bits >> 4, rsc->sse);
}

static int32_t rest_tiles_in_plane(const Av1Common *cm, int32_t plane) {
    const RestorationInfo *rsi = &cm->rst_info[plane];
    return rsi->units_per_tile;
}

void *aom_memalign(size_t align, size_t size);
void aom_free(void *memblk);

void av1_pick_filter_restoration(const Yv12BufferConfig *src, Yv12BufferConfig * trial_frame_rst /*AV1_COMP *cpi*/, Macroblock *x, Av1Common *const cm) {

    //CHKN Av1Common *const cm = &cpi->common;
    const int32_t num_planes = 3;// av1_num_planes(cm);
   // assert(!cm->all_lossless);

    int32_t ntiles[2];
    for (int32_t is_uv = 0; is_uv < 2; ++is_uv)
        ntiles[is_uv] = rest_tiles_in_plane(cm, is_uv);

    assert(ntiles[1] <= ntiles[0]);
    RestUnitSearchInfo *rusi =
        (RestUnitSearchInfo *)aom_memalign(16, sizeof(*rusi) * ntiles[0]);

    // If the restoration unit dimensions are not multiples of
    // rsi->restoration_unit_size then some elements of the rusi array may be
    // left uninitialised when we reach copy_unit_info(...). This is not a
    // problem, as these elements are ignored later, but in order to quiet
    // Valgrind's warnings we initialise the array below.
    memset(rusi, 0, sizeof(*rusi) * ntiles[0]);

    RestSearchCtxt rsc;
    const int32_t plane_start = AOM_PLANE_Y;
    const int32_t plane_end = num_planes > 1 ? AOM_PLANE_V : AOM_PLANE_Y;
    for (int32_t plane = plane_start; plane <= plane_end; ++plane) {
        //CHKN  init_rsc(src, cpi->common, &cpi->td.mb, plane, rusi, &cpi->trial_frame_rst,&rsc);
        init_rsc(src, cm, x, plane, rusi, trial_frame_rst, &rsc);

        const int32_t plane_ntiles = ntiles[plane > 0];
        const RestorationType num_rtypes =
            (plane_ntiles > 1) ? RESTORE_TYPES : RESTORE_SWITCHABLE_TYPES;

        double best_cost = 0;
        RestorationType best_rtype = RESTORE_NONE;

        const int32_t highbd = rsc.cm->use_highbitdepth;




        extend_frame(rsc.dgd_buffer, rsc.plane_width, rsc.plane_height,
            rsc.dgd_stride, RESTORATION_BORDER, RESTORATION_BORDER,
            highbd);

        //CHKN  for (RestorationType r = 0; r < num_rtypes; ++r) {
        for (int32_t restType = 0; restType < num_rtypes; ++restType) {


            RestorationType r = (RestorationType)restType;

            if ((force_restore_type != RESTORE_TYPES) && (r != RESTORE_NONE) &&
                (r != force_restore_type))
                continue;

            double cost = search_rest_type(&rsc, r);

            if (r == 0 || cost < best_cost)
            {
                best_cost = cost;
                best_rtype = r;
            }


        }



        cm->rst_info[plane].frame_restoration_type = best_rtype;
        if (force_restore_type != RESTORE_TYPES)
            assert(best_rtype == force_restore_type || best_rtype == RESTORE_NONE);

        if (best_rtype != RESTORE_NONE) {
            for (int32_t u = 0; u < plane_ntiles; ++u) {

                copy_unit_info(best_rtype, &rusi[u], &cm->rst_info[plane].unit_info[u]);

            }
        }
    }

    aom_free(rusi);
}


#if REST_M



static void search_sgrproj_seg(const RestorationTileLimits *limits,
    const AV1PixelRect *tile, int32_t rest_unit_idx,
    void *priv) {
    RestSearchCtxt *rsc = (RestSearchCtxt *)priv;
    RestUnitSearchInfo *rusi = &rsc->rusi[rest_unit_idx];

#if FAST_SG
    Av1Common *const cm = rsc->cm;
#else
    const Av1Common *const cm = rsc->cm;
#endif
    const int32_t highbd = cm->use_highbitdepth;
    const int32_t bit_depth = cm->bit_depth;

    uint8_t *dgd_start =
        rsc->dgd_buffer + limits->v_start * rsc->dgd_stride + limits->h_start;
    const uint8_t *src_start =
        rsc->src_buffer + limits->v_start * rsc->src_stride + limits->h_start;

    const int32_t is_uv = rsc->plane > 0;
    const int32_t ss_x = is_uv && cm->subsampling_x;
    const int32_t ss_y = is_uv && cm->subsampling_y;
    const int32_t procunit_width = RESTORATION_PROC_UNIT_SIZE >> ss_x;
    const int32_t procunit_height = RESTORATION_PROC_UNIT_SIZE >> ss_y;
#if FAST_SG
    int8_t step = cm->sg_filter_mode == 1 ? 0 : cm->sg_filter_mode == 2 ? 1 : cm->sg_filter_mode == 3 ? 4 : 16;
#endif
    rusi->sgrproj = search_selfguided_restoration(
        dgd_start, limits->h_end - limits->h_start,
        limits->v_end - limits->v_start, rsc->dgd_stride, src_start,
        rsc->src_stride, highbd, bit_depth, procunit_width, procunit_height,
        rsc->tmpbuf
#if FAST_SG
        , cm->sg_ref_frame_ep,
        cm->sg_frame_ep_cnt,
        step
#endif
    );


    RestorationUnitInfo rui;
    rui.restoration_type = RESTORE_SGRPROJ;
    rui.sgrproj_info = rusi->sgrproj;

    rusi->sse[RESTORE_SGRPROJ] = try_restoration_unit_seg(rsc, limits, tile, &rui);

   
}

static void search_sgrproj_finish(const RestorationTileLimits *limits,
    const AV1PixelRect *tile, int32_t rest_unit_idx,
    void *priv) {

    (void)limits;
    (void)tile;
    RestSearchCtxt *rsc = (RestSearchCtxt *)priv;
    RestUnitSearchInfo *rusi = &rsc->rusi[rest_unit_idx];

    const Macroblock *const x = rsc->x;
   

    rusi->sse[RESTORE_SGRPROJ] = rsc->rusi_pic[rest_unit_idx].sse[RESTORE_SGRPROJ];
    rusi->sgrproj = rsc->rusi_pic[rest_unit_idx].sgrproj;
  
    const int64_t bits_none = x->sgrproj_restore_cost[0];
    const int64_t bits_sgr = x->sgrproj_restore_cost[1] +
        (count_sgrproj_bits(&rusi->sgrproj, &rsc->sgrproj)
            << AV1_PROB_COST_SHIFT);

    double cost_none =
        RDCOST_DBL(x->rdmult, bits_none >> 4, rusi->sse[RESTORE_NONE]);
    double cost_sgr =
        RDCOST_DBL(x->rdmult, bits_sgr >> 4, rusi->sse[RESTORE_SGRPROJ]);

    RestorationType rtype =
        (cost_sgr < cost_none) ? RESTORE_SGRPROJ : RESTORE_NONE;
    rusi->best_rtype[RESTORE_SGRPROJ - 1] = rtype;

    rsc->sse += rusi->sse[rtype];
    rsc->bits += (cost_sgr < cost_none) ? bits_sgr : bits_none;
    if (cost_sgr < cost_none) rsc->sgrproj = rusi->sgrproj;
}

static void search_wiener_seg(const RestorationTileLimits *limits,
    const AV1PixelRect *tile_rect, int32_t rest_unit_idx,
    void *priv) {
    RestSearchCtxt *rsc = (RestSearchCtxt *)priv;
    RestUnitSearchInfo *rusi = &rsc->rusi[rest_unit_idx];
    const Av1Common *const cm = rsc->cm;
#if FAST_WN
    int32_t wn_luma = cm->wn_filter_mode == 1 ? WIENER_WIN_CHROMA : WIENER_WIN;
    const int32_t wiener_win =
        (rsc->plane == AOM_PLANE_Y) ? wn_luma : WIENER_WIN_CHROMA;
#else
    const int32_t wiener_win =
        (rsc->plane == AOM_PLANE_Y) ? WIENER_WIN : WIENER_WIN_CHROMA;
#endif
    int64_t M[WIENER_WIN2];
    int64_t H[WIENER_WIN2 * WIENER_WIN2];
    int32_t vfilterd[WIENER_WIN], hfilterd[WIENER_WIN];

    
    if (cm->use_highbitdepth)
    {
        if (rsc->plane == AOM_PLANE_Y) {

            av1_compute_stats_highbd(WIENER_WIN, rsc->dgd_buffer, rsc->src_buffer,

                limits->h_start, limits->h_end, limits->v_start,
                limits->v_end, rsc->dgd_stride, rsc->src_stride, M,
                H, (aom_bit_depth_t)cm->bit_depth);

        }
        else {

            av1_compute_stats_highbd(WIENER_WIN_CHROMA, rsc->dgd_buffer, rsc->src_buffer,

                limits->h_start, limits->h_end, limits->v_start,
                limits->v_end, rsc->dgd_stride, rsc->src_stride, M,
                H, (aom_bit_depth_t)cm->bit_depth);

        }
    }
    else {

        av1_compute_stats(wiener_win, rsc->dgd_buffer, rsc->src_buffer, limits->h_start,
            limits->h_end, limits->v_start, limits->v_end,
            rsc->dgd_stride, rsc->src_stride, M, H);

    }

   

    if (!wiener_decompose_sep_sym(wiener_win, M, H, vfilterd, hfilterd)) {

        printf("CHKN never get here\n");      
        rusi->best_rtype[RESTORE_WIENER - 1] = RESTORE_NONE;
        rusi->sse[RESTORE_WIENER] = INT64_MAX;
        return;
    }

    RestorationUnitInfo rui;
    memset(&rui, 0, sizeof(rui));
    rui.restoration_type = RESTORE_WIENER;
    finalize_sym_filter(wiener_win, vfilterd, rui.wiener_info.vfilter);
    finalize_sym_filter(wiener_win, hfilterd, rui.wiener_info.hfilter);

    // Filter score computes the value of the function x'*A*x - x'*b for the
    // learned filter and compares it against identity filer. If there is no
    // reduction in the function, the filter is reverted back to identity
    if (compute_score(wiener_win, M, H, rui.wiener_info.vfilter,
        rui.wiener_info.hfilter) > 0) {
     
        rusi->sse[RESTORE_WIENER] = INT64_MAX;
        return;
    }

    aom_clear_system_state();

    rusi->sse[RESTORE_WIENER] =
        finer_tile_search_wiener_seg(rsc, limits, tile_rect, &rui, wiener_win);
    rusi->wiener = rui.wiener_info;

    if (wiener_win != WIENER_WIN) {
        assert(rui.wiener_info.vfilter[0] == 0 &&
            rui.wiener_info.vfilter[WIENER_WIN - 1] == 0);
        assert(rui.wiener_info.hfilter[0] == 0 &&
            rui.wiener_info.hfilter[WIENER_WIN - 1] == 0);
    }
    
   
}
static void search_wiener_finish(const RestorationTileLimits *limits,
    const AV1PixelRect *tile_rect, int32_t rest_unit_idx,
    void *priv) {
    (void)limits;
    (void)tile_rect;
    RestSearchCtxt *rsc = (RestSearchCtxt *)priv;
    RestUnitSearchInfo *rusi = &rsc->rusi[rest_unit_idx];
#if FAST_WN
    const Av1Common *const cm = rsc->cm;
    int32_t wn_luma = cm->wn_filter_mode == 1 ? WIENER_WIN_CHROMA : WIENER_WIN;
    const int32_t wiener_win =
        (rsc->plane == AOM_PLANE_Y) ? WIENER_WIN_CHROMA : WIENER_WIN_CHROMA;
#else
    const int32_t wiener_win =
        (rsc->plane == AOM_PLANE_Y) ? WIENER_WIN : WIENER_WIN_CHROMA;
#endif
    
   
    const Macroblock *const x = rsc->x;
    const int64_t bits_none = x->wiener_restore_cost[0];

   

    RestorationUnitInfo rui;
    memset(&rui, 0, sizeof(rui));
    rui.restoration_type = RESTORE_WIENER;
    
    // Filter score computes the value of the function x'*A*x - x'*b for the
    // learned filter and compares it against identity filer. If there is no
    // reduction in the function, the filter is reverted back to identity

    rusi->sse[RESTORE_WIENER] = rsc->rusi_pic[rest_unit_idx].sse[RESTORE_WIENER];
    if (rusi->sse[RESTORE_WIENER] == INT64_MAX) 
    {
        rsc->bits += bits_none;
        rsc->sse += rusi->sse[RESTORE_NONE];
        rusi->best_rtype[RESTORE_WIENER - 1] = RESTORE_NONE;
        rusi->sse[RESTORE_WIENER] = INT64_MAX;
        return;
    }

    aom_clear_system_state();

    rusi->wiener = rsc->rusi_pic[rest_unit_idx].wiener;
         

    const int64_t bits_wiener =
        x->wiener_restore_cost[1] +
        (count_wiener_bits(wiener_win, &rusi->wiener, &rsc->wiener)
            << AV1_PROB_COST_SHIFT);

    double cost_none =
        RDCOST_DBL(x->rdmult, bits_none >> 4, rusi->sse[RESTORE_NONE]);
    double cost_wiener =
        RDCOST_DBL(x->rdmult, bits_wiener >> 4, rusi->sse[RESTORE_WIENER]);

    RestorationType rtype =
        (cost_wiener < cost_none) ? RESTORE_WIENER : RESTORE_NONE;
    rusi->best_rtype[RESTORE_WIENER - 1] = rtype;

    rsc->sse += rusi->sse[rtype];
    rsc->bits += (cost_wiener < cost_none) ? bits_wiener : bits_none;
    if (cost_wiener < cost_none) rsc->wiener = rusi->wiener;
}
static void search_norestore_seg(const RestorationTileLimits *limits,
    const AV1PixelRect *tile_rect, int32_t rest_unit_idx,
    void *priv) {
    (void)tile_rect;

    RestSearchCtxt *rsc = (RestSearchCtxt *)priv;
    RestUnitSearchInfo *rusi = &rsc->rusi[rest_unit_idx];

    const int32_t highbd = rsc->cm->use_highbitdepth;
    rusi->sse[RESTORE_NONE] = sse_restoration_unit(
        limits, rsc->src, rsc->cm->frame_to_show, rsc->plane, highbd);

}
static void search_norestore_finish(const RestorationTileLimits *limits,
    const AV1PixelRect *tile_rect, int32_t rest_unit_idx,
    void *priv) {
    (void)tile_rect;
    (void)limits;

    RestSearchCtxt *rsc = (RestSearchCtxt *)priv;
    RestUnitSearchInfo *rusi = &rsc->rusi[rest_unit_idx];

    rusi->sse[RESTORE_NONE] = rsc->rusi_pic[rest_unit_idx].sse[RESTORE_NONE];

    rsc->sse += rusi->sse[RESTORE_NONE];  
}
static double search_rest_type_finish(RestSearchCtxt *rsc, RestorationType rtype)
{

    static const rest_unit_visitor_t funs[RESTORE_TYPES] = {
    search_norestore_finish, search_wiener_finish, search_sgrproj_finish, search_switchable
    };

    reset_rsc(rsc);

    av1_foreach_rest_unit_in_frame(rsc->cm, rsc->plane, rsc_on_tile, funs[rtype], rsc);

    return RDCOST_DBL(rsc->x->rdmult, rsc->bits >> 4, rsc->sse);
}

void restoration_seg_search(
    RestContext_t          *context_ptr,
    Yv12BufferConfig       *org_fts,
    const Yv12BufferConfig *src,
    Yv12BufferConfig       *trial_frame_rst ,
    PictureControlSet_t    *pcs_ptr,
    uint32_t                segment_index )
{
    Av1Common *const cm = pcs_ptr->parent_pcs_ptr->av1_cm;
    Macroblock *x = pcs_ptr->parent_pcs_ptr->av1x;   
    const int32_t num_planes = 3;

    
    // If the restoration unit dimensions are not multiples of
    // rsi->restoration_unit_size then some elements of the rusi array may be
    // left uninitialised when we reach copy_unit_info(...). This is not a
    // problem, as these elements are ignored later, but in order to quiet
    // Valgrind's warnings we initialise the array  to zero.


    RestSearchCtxt rsc; //this context is specific for this segment
    RestSearchCtxt* rsc_p = &rsc;

    const int32_t plane_start = AOM_PLANE_Y;
    const int32_t plane_end = num_planes > 1 ? AOM_PLANE_V : AOM_PLANE_Y;
    for (int32_t plane = plane_start; plane <= plane_end; ++plane)
    {
        RestUnitSearchInfo *rusi = pcs_ptr->parent_pcs_ptr->rusi_picture[plane];
              
        init_rsc_seg(org_fts,src, cm, x, plane, rusi, trial_frame_rst, &rsc);
       
        rsc_p->tmpbuf = context_ptr->rst_tmpbuf;

    
        const int32_t highbd = rsc.cm->use_highbitdepth;
        extend_frame(rsc.dgd_buffer, rsc.plane_width, rsc.plane_height,
            rsc.dgd_stride, RESTORATION_BORDER, RESTORATION_BORDER,
            highbd);       

        av1_foreach_rest_unit_in_frame_seg(rsc_p->cm, rsc_p->plane, rsc_on_tile, search_norestore_seg, rsc_p, pcs_ptr, segment_index);
        av1_foreach_rest_unit_in_frame_seg(rsc_p->cm, rsc_p->plane, rsc_on_tile, search_wiener_seg,  rsc_p, pcs_ptr, segment_index);
        av1_foreach_rest_unit_in_frame_seg(rsc_p->cm, rsc_p->plane, rsc_on_tile, search_sgrproj_seg, rsc_p, pcs_ptr, segment_index);       

    }

    
}
void rest_finish_search(Macroblock *x, Av1Common *const cm)
{      
    const int32_t num_planes = 3;

    int32_t ntiles[2];
    for (int32_t is_uv = 0; is_uv < 2; ++is_uv)
        ntiles[is_uv] = rest_tiles_in_plane(cm, is_uv);

    assert(ntiles[1] <= ntiles[0]);
    RestUnitSearchInfo *rusi =
        (RestUnitSearchInfo *)aom_memalign(16, sizeof(*rusi) * ntiles[0]);

    // If the restoration unit dimensions are not multiples of
    // rsi->restoration_unit_size then some elements of the rusi array may be
    // left uninitialised when we reach copy_unit_info(...). This is not a
    // problem, as these elements are ignored later, but in order to quiet
    // Valgrind's warnings we initialise the array below.
    memset(rusi, 0, sizeof(*rusi) * ntiles[0]);

    RestSearchCtxt rsc;
    const int32_t plane_start = AOM_PLANE_Y;
    const int32_t plane_end = num_planes > 1 ? AOM_PLANE_V : AOM_PLANE_Y;
    for (int32_t plane = plane_start; plane <= plane_end; ++plane) {
       
        //init rsc context for this plane
        rsc.cm = cm;
        rsc.x = x;
        rsc.plane = plane;
        rsc.rusi = rusi;
        rsc.pic_num = cm->p_pcs_ptr->picture_number;
        rsc.rusi_pic = cm->p_pcs_ptr->rusi_picture[plane];

        const int32_t plane_ntiles = ntiles[plane > 0];
        const RestorationType num_rtypes =
            (plane_ntiles > 1) ? RESTORE_TYPES : RESTORE_SWITCHABLE_TYPES;

        double best_cost = 0;
        RestorationType best_rtype = RESTORE_NONE;
       
       
        for (int32_t restType = 0; restType < num_rtypes; ++restType) {

            RestorationType r = (RestorationType)restType;

            if ((force_restore_type != RESTORE_TYPES) && (r != RESTORE_NONE) &&
                (r != force_restore_type))
                continue;

            double cost = search_rest_type_finish(&rsc, r);

            if (r == 0 || cost < best_cost)
            {
                best_cost = cost;
                best_rtype = r;
            }
        }

        cm->rst_info[plane].frame_restoration_type = best_rtype;
        if (force_restore_type != RESTORE_TYPES)
            assert(best_rtype == force_restore_type || best_rtype == RESTORE_NONE);

        if (best_rtype != RESTORE_NONE) {
            for (int32_t u = 0; u < plane_ntiles; ++u) {
                copy_unit_info(best_rtype, &rusi[u], &cm->rst_info[plane].unit_info[u]);
            }
        }
    }

    aom_free(rusi);
}
#endif