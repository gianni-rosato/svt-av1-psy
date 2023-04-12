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

#include "aom_dsp_rtcd.h"
#include "corner_match.h"

#define SEARCH_SZ 9
#define SEARCH_SZ_BY2 ((SEARCH_SZ - 1) / 2)

#define THRESHOLD_NCC 0.75

/* Compute var(im) * MATCH_SZ_SQ over a MATCH_SZ by MATCH_SZ window of im,
   centered at (x, y).
*/
#if OPT_GM_WD
static int32_t compute_variance(unsigned char *im, int stride, int x, int y, uint8_t match_sz) {
#else
static int32_t compute_variance(unsigned char *im, int stride, int x, int y) {
#endif
#if OPT_GM_WD
    const uint8_t match_sz_by2 = ((match_sz - 1) / 2);
    const uint8_t match_sz_sq  = (match_sz * match_sz);
#endif
    int sum   = 0;
    int sumsq = 0;
    int var;
    int i, j;
#if OPT_GM_WD
    for (i = 0; i < match_sz; ++i)
        for (j = 0; j < match_sz; ++j) {
            sum += im[(i + y - match_sz_by2) * stride + (j + x - match_sz_by2)];
            sumsq += im[(i + y - match_sz_by2) * stride + (j + x - match_sz_by2)] *
                im[(i + y - match_sz_by2) * stride + (j + x - match_sz_by2)];
        }
    var = sumsq * match_sz_sq - sum * sum;
#else
    for (i = 0; i < MATCH_SZ; ++i)
        for (j = 0; j < MATCH_SZ; ++j) {
            sum += im[(i + y - MATCH_SZ_BY2) * stride + (j + x - MATCH_SZ_BY2)];
            sumsq += im[(i + y - MATCH_SZ_BY2) * stride + (j + x - MATCH_SZ_BY2)] *
                im[(i + y - MATCH_SZ_BY2) * stride + (j + x - MATCH_SZ_BY2)];
        }
    var = sumsq * MATCH_SZ_SQ - sum * sum;
#endif

    return var;
}

/* Compute quad of corr(im1, im2) * MATCH_SZ * stddev(im1), where the
   correlation/standard deviation are taken over MATCH_SZ by MATCH_SZ windows
   of each image, centered at (x1, y1) and (x2, y2) respectively.
*/
double svt_av1_compute_cross_correlation_c(unsigned char *im1, int stride1, int x1, int y1,
#if OPT_GM_WD
                                           unsigned char *im2, int stride2, int x2, int y2,
                                           uint8_t match_sz) {
#else
                                           unsigned char *im2, int stride2, int x2, int y2) {
#endif
#if OPT_GM_WD
    const uint8_t match_sz_by2 = ((match_sz - 1) / 2);
    const uint8_t match_sz_sq  = (match_sz * match_sz);
#endif

    int v1, v2;
    int sum1   = 0;
    int sum2   = 0;
    int sumsq2 = 0;
    int cross  = 0;
    int var2, cov;
    int i, j;
#if OPT_GM_WD
    for (i = 0; i < match_sz; ++i)
        for (j = 0; j < match_sz; ++j) {
            v1 = im1[(i + y1 - match_sz_by2) * stride1 + (j + x1 - match_sz_by2)];
            v2 = im2[(i + y2 - match_sz_by2) * stride2 + (j + x2 - match_sz_by2)];
            sum1 += v1;
            sum2 += v2;
            sumsq2 += v2 * v2;
            cross += v1 * v2;
        }
    var2 = sumsq2 * match_sz_sq - sum2 * sum2;
    cov  = cross * match_sz_sq - sum1 * sum2;
#else
    for (i = 0; i < MATCH_SZ; ++i)
        for (j = 0; j < MATCH_SZ; ++j) {
            v1 = im1[(i + y1 - MATCH_SZ_BY2) * stride1 + (j + x1 - MATCH_SZ_BY2)];
            v2 = im2[(i + y2 - MATCH_SZ_BY2) * stride2 + (j + x2 - MATCH_SZ_BY2)];
            sum1 += v1;
            sum2 += v2;
            sumsq2 += v2 * v2;
            cross += v1 * v2;
        }
    var2 = sumsq2 * MATCH_SZ_SQ - sum2 * sum2;
    cov  = cross * MATCH_SZ_SQ - sum1 * sum2;
#endif
    if (cov < 0) {
        return 0;
    }
    return ((double)cov * cov) / ((double)var2);
}

#if OPT_GM_WD
static INLINE int is_eligible_point(int pointx, int pointy, int width, int height,
                                    uint8_t match_sz) {
    const uint8_t match_sz_by2 = ((match_sz - 1) / 2);

    return (pointx >= match_sz_by2 && pointy >= match_sz_by2 && pointx + match_sz_by2 < width &&
            pointy + match_sz_by2 < height);
}
#else
static INLINE int is_eligible_point(int pointx, int pointy, int width, int height) {
    return (pointx >= MATCH_SZ_BY2 && pointy >= MATCH_SZ_BY2 && pointx + MATCH_SZ_BY2 < width &&
            pointy + MATCH_SZ_BY2 < height);
}
#endif

static INLINE int is_eligible_distance(int point1x, int point1y, int point2x, int point2y,
                                       int threshSqr) {
    const int xdist = point1x - point2x;
    const int ydist = point1y - point2y;
    return (xdist * xdist + ydist * ydist) <= threshSqr;
}

static void improve_correspondence(unsigned char *frm, unsigned char *ref, int width, int height,
                                   int frm_stride, int ref_stride, Correspondence *correspondences,
#if OPT_GM_WD
                                   int num_correspondences, uint8_t match_sz) {
#else
                                   int num_correspondences) {
#endif
    int       i;
    const int thresh    = (width < height ? height : width) >> 4;
    const int threshSqr = thresh * thresh;
    for (i = 0; i < num_correspondences; ++i) {
        int    x, y, best_x = 0, best_y = 0;
        double best_match_ncc = 0.0;
        for (y = -SEARCH_SZ_BY2; y <= SEARCH_SZ_BY2; ++y) {
            for (x = -SEARCH_SZ_BY2; x <= SEARCH_SZ_BY2; ++x) {
                double match_ncc;
#if OPT_GM_WD
                if (!is_eligible_point(correspondences[i].rx + x,
                                       correspondences[i].ry + y,
                                       width,
                                       height,
                                       match_sz))
#else
                if (!is_eligible_point(
                        correspondences[i].rx + x, correspondences[i].ry + y, width, height))
#endif
                    continue;
                if (!is_eligible_distance(correspondences[i].x,
                                          correspondences[i].y,
                                          correspondences[i].rx + x,
                                          correspondences[i].ry + y,
                                          threshSqr))
                    continue;
                match_ncc = svt_av1_compute_cross_correlation(frm,
                                                              frm_stride,
                                                              correspondences[i].x,
                                                              correspondences[i].y,
                                                              ref,
                                                              ref_stride,
                                                              correspondences[i].rx + x,
#if OPT_GM_WD
                                                              correspondences[i].ry + y,
                                                              match_sz);
#else
                                                              correspondences[i].ry + y);
#endif
                if (match_ncc > best_match_ncc) {
                    best_match_ncc = match_ncc;
                    best_y         = y;
                    best_x         = x;
                }
            }
        }
        correspondences[i].rx += best_x;
        correspondences[i].ry += best_y;
    }
    for (i = 0; i < num_correspondences; ++i) {
        int    x, y, best_x = 0, best_y = 0;
        double best_match_ncc = 0.0;
        for (y = -SEARCH_SZ_BY2; y <= SEARCH_SZ_BY2; ++y)
            for (x = -SEARCH_SZ_BY2; x <= SEARCH_SZ_BY2; ++x) {
                double match_ncc;
#if OPT_GM_WD
                if (!is_eligible_point(correspondences[i].x + x,
                                       correspondences[i].y + y,
                                       width,
                                       height,
                                       match_sz))
#else
                if (!is_eligible_point(
                        correspondences[i].x + x, correspondences[i].y + y, width, height))
#endif
                    continue;
                if (!is_eligible_distance(correspondences[i].x + x,
                                          correspondences[i].y + y,
                                          correspondences[i].rx,
                                          correspondences[i].ry,
                                          threshSqr))
                    continue;
#if OPT_GM_WD
                match_ncc = svt_av1_compute_cross_correlation(ref,
                                                              ref_stride,
                                                              correspondences[i].rx,
                                                              correspondences[i].ry,
                                                              frm,
                                                              frm_stride,
                                                              correspondences[i].x + x,
                                                              correspondences[i].y + y,
                                                              match_sz);
#else
                match_ncc = svt_av1_compute_cross_correlation(ref,
                                                              ref_stride,
                                                              correspondences[i].rx,
                                                              correspondences[i].ry,
                                                              frm,
                                                              frm_stride,
                                                              correspondences[i].x + x,
                                                              correspondences[i].y + y);
#endif
                if (match_ncc > best_match_ncc) {
                    best_match_ncc = match_ncc;
                    best_y         = y;
                    best_x         = x;
                }
            }
        correspondences[i].x += best_x;
        correspondences[i].y += best_y;
    }
}

int svt_av1_determine_correspondence(unsigned char *frm, int *frm_corners, int num_frm_corners,
                                     unsigned char *ref, int *ref_corners, int num_ref_corners,
                                     int width, int height, int frm_stride, int ref_stride,
#if OPT_GM_WD
                                     int *correspondence_pts, uint8_t match_sz) {
#else
                                     int *correspondence_pts) {
#endif
    int             i, j;
    Correspondence *correspondences     = (Correspondence *)correspondence_pts;
    int             num_correspondences = 0;
    const int       thresh              = (width < height ? height : width) >> 4;
    const int       threshSqr           = thresh * thresh;
    for (i = 0; i < num_frm_corners; ++i) {
        double  best_match_ncc = 0.0;
        int32_t template_norm;
        int     best_match_j = -1;
#if OPT_GM_WD
        if (!is_eligible_point(frm_corners[2 * i], frm_corners[2 * i + 1], width, height, match_sz))
#else
        if (!is_eligible_point(frm_corners[2 * i], frm_corners[2 * i + 1], width, height))
#endif
            continue;
        for (j = 0; j < num_ref_corners; ++j) {
            double match_ncc;
#if OPT_GM_WD
            if (!is_eligible_point(
                    ref_corners[2 * j], ref_corners[2 * j + 1], width, height, match_sz))
#else
            if (!is_eligible_point(ref_corners[2 * j], ref_corners[2 * j + 1], width, height))
#endif
                continue;
            if (!is_eligible_distance(frm_corners[2 * i],
                                      frm_corners[2 * i + 1],
                                      ref_corners[2 * j],
                                      ref_corners[2 * j + 1],
                                      threshSqr))
                continue;
#if OPT_GM_WD
            match_ncc = svt_av1_compute_cross_correlation(frm,
                                                          frm_stride,
                                                          frm_corners[2 * i],
                                                          frm_corners[2 * i + 1],
                                                          ref,
                                                          ref_stride,
                                                          ref_corners[2 * j],
                                                          ref_corners[2 * j + 1],
                                                          match_sz);
#else
            match_ncc = svt_av1_compute_cross_correlation(frm,
                                                          frm_stride,
                                                          frm_corners[2 * i],
                                                          frm_corners[2 * i + 1],
                                                          ref,
                                                          ref_stride,
                                                          ref_corners[2 * j],
                                                          ref_corners[2 * j + 1]);
#endif
            if (match_ncc > best_match_ncc) {
                best_match_ncc = match_ncc;
                best_match_j   = j;
            }
        }
        // Note: We want to test if the best correlation is >= THRESHOLD_NCC,
        // but need to account for the normalization in
        // av1_compute_cross_correlation.
#if OPT_GM_WD
        template_norm = compute_variance(
            frm, frm_stride, frm_corners[2 * i], frm_corners[2 * i + 1], match_sz);

#else
        template_norm = compute_variance(
            frm, frm_stride, frm_corners[2 * i], frm_corners[2 * i + 1]);
#endif
        if (best_match_ncc > (template_norm * THRESHOLD_NCC * THRESHOLD_NCC)) {
            correspondences[num_correspondences].x  = frm_corners[2 * i];
            correspondences[num_correspondences].y  = frm_corners[2 * i + 1];
            correspondences[num_correspondences].rx = ref_corners[2 * best_match_j];
            correspondences[num_correspondences].ry = ref_corners[2 * best_match_j + 1];

            /*SVT_LOG("corresp: %d %d - %d %d\n",
             correspondences[num_correspondences].x,
             correspondences[num_correspondences].y,
             correspondences[num_correspondences].rx,
             correspondences[num_correspondences].ry);*/

            num_correspondences++;
        }
    }
#if OPT_GM_WD
    improve_correspondence(frm,
                           ref,
                           width,
                           height,
                           frm_stride,
                           ref_stride,
                           correspondences,
                           num_correspondences,
                           match_sz);
#else
    improve_correspondence(
        frm, ref, width, height, frm_stride, ref_stride, correspondences, num_correspondences);
#endif
    return num_correspondences;
}
