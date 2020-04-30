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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "global_motion.h"
#include "EbUtility.h"
#include "corner_detect.h"
#include "corner_match.h"
#include "ransac.h"

#include "EbEncWarpedMotion.h"

#define MIN_INLIER_PROB 0.1

#define MIN_TRANS_THRESH (1 * GM_TRANS_DECODE_FACTOR)

// Border over which to compute the global motion
#define ERRORADV_BORDER 0

// TODO(sarahparker) These need to be retuned for speed 0 and 1 to
// maximize gains from segmented error metric
static const double erroradv_tr[]      = {0.65, 0.60, 0.65};
static const double erroradv_prod_tr[] = {20000, 18000, 16000};

int av1_is_enough_erroradvantage(double best_erroradvantage, int params_cost, int erroradv_type) {
    assert(erroradv_type < GM_ERRORADV_TR_TYPES);
    return best_erroradvantage < erroradv_tr[erroradv_type] &&
           best_erroradvantage * params_cost < erroradv_prod_tr[erroradv_type];
}

static void convert_to_params(const double *params, int32_t *model) {
    int i;
    int alpha_present = 0;
    model[0]          = (int32_t)floor(params[0] * (1 << GM_TRANS_PREC_BITS) + 0.5);
    model[1]          = (int32_t)floor(params[1] * (1 << GM_TRANS_PREC_BITS) + 0.5);
    model[0] = (int32_t)clamp(model[0], GM_TRANS_MIN, GM_TRANS_MAX) * GM_TRANS_DECODE_FACTOR;
    model[1] = (int32_t)clamp(model[1], GM_TRANS_MIN, GM_TRANS_MAX) * GM_TRANS_DECODE_FACTOR;

    for (i = 2; i < 6; ++i) {
        const int diag_value = ((i == 2 || i == 5) ? (1 << GM_ALPHA_PREC_BITS) : 0);
        model[i]             = (int32_t)floor(params[i] * (1 << GM_ALPHA_PREC_BITS) + 0.5);
        model[i]             = (int32_t)clamp(model[i] - diag_value, GM_ALPHA_MIN, GM_ALPHA_MAX);
        alpha_present |= (model[i] != 0);
        model[i] = (model[i] + diag_value) * GM_ALPHA_DECODE_FACTOR;
    }
    for (; i < 8; ++i) {
        model[i] = (int32_t)floor(params[i] * (1 << GM_ROW3HOMO_PREC_BITS) + 0.5);
        model[i] =
            (int32_t)clamp(model[i], GM_ROW3HOMO_MIN, GM_ROW3HOMO_MAX) * GM_ROW3HOMO_DECODE_FACTOR;
        alpha_present |= (model[i] != 0);
    }

    if (!alpha_present) {
        if (abs(model[0]) < MIN_TRANS_THRESH && abs(model[1]) < MIN_TRANS_THRESH) {
            model[0] = 0;
            model[1] = 0;
        }
    }
}

static INLINE TransformationType get_wmtype(const EbWarpedMotionParams *gm) {
    if (gm->wmmat[5] == (1 << WARPEDMODEL_PREC_BITS) && !gm->wmmat[4] &&
        gm->wmmat[2] == (1 << WARPEDMODEL_PREC_BITS) && !gm->wmmat[3]) {
        return ((!gm->wmmat[1] && !gm->wmmat[0]) ? IDENTITY : TRANSLATION);
    }
    if (gm->wmmat[2] == gm->wmmat[5] && gm->wmmat[3] == -gm->wmmat[4])
        return ROTZOOM;
    else
        return AFFINE;
}

void av1_convert_model_to_params(const double *params, EbWarpedMotionParams *model) {
    convert_to_params(params, model->wmmat);
    model->wmtype  = get_wmtype(model);
    model->invalid = 0;
}

// Adds some offset to a global motion parameter and handles
// all of the necessary precision shifts, clamping, and
// zero-centering.
static int32_t add_param_offset(int param_index, int32_t param_value, int32_t offset) {
    const int scale_vals[3] = {GM_TRANS_PREC_DIFF, GM_ALPHA_PREC_DIFF, GM_ROW3HOMO_PREC_DIFF};
    const int clamp_vals[3] = {GM_TRANS_MAX, GM_ALPHA_MAX, GM_ROW3HOMO_MAX};
    // type of param: 0 - translation, 1 - affine, 2 - homography
    const int param_type      = (param_index < 2 ? 0 : (param_index < 6 ? 1 : 2));
    const int is_one_centered = (param_index == 2 || param_index == 5);

    // Make parameter zero-centered and offset the shift that was done to make
    // it compatible with the warped model
    param_value =
        (param_value - (is_one_centered << WARPEDMODEL_PREC_BITS)) >> scale_vals[param_type];
    // Add desired offset to the rescaled/zero-centered parameter
    param_value += offset;
    // Clamp the parameter so it does not overflow the number of bits allotted
    // to it in the Bitstream
    param_value = (int32_t)clamp(param_value, -clamp_vals[param_type], clamp_vals[param_type]);
    // Rescale the parameter to WARPEDMODEL_PRECISION_BITS so it is compatible
    // with the warped motion library
    param_value *= (1 << scale_vals[param_type]);

    // Undo the zero-centering step if necessary
    return param_value + (is_one_centered << WARPEDMODEL_PREC_BITS);
}

static void force_wmtype(EbWarpedMotionParams *wm, TransformationType wmtype) {
    switch (wmtype) {
    case IDENTITY:
        wm->wmmat[0] = 0;
        wm->wmmat[1] = 0;
        AOM_FALLTHROUGH_INTENDED;
    case TRANSLATION:
        wm->wmmat[2] = 1 << WARPEDMODEL_PREC_BITS;
        wm->wmmat[3] = 0;
        AOM_FALLTHROUGH_INTENDED;
    case ROTZOOM:
        wm->wmmat[4] = -wm->wmmat[3];
        wm->wmmat[5] = wm->wmmat[2];
        AOM_FALLTHROUGH_INTENDED;
    case AFFINE: wm->wmmat[6] = wm->wmmat[7] = 0; break;
    default: assert(0);
    }
    wm->wmtype = wmtype;
}

int64_t av1_refine_integerized_param(EbWarpedMotionParams *wm, TransformationType wmtype,
                                     int use_hbd, int bd, uint8_t *ref, int r_width, int r_height,
                                     int r_stride, uint8_t *dst, int d_width, int d_height,
                                     int d_stride, int n_refinements, int64_t best_frame_error) {
    static const int max_trans_model_params[TRANS_TYPES] = {0, 2, 4, 6};
    const int        border                              = ERRORADV_BORDER;
    int              i                                   = 0, p;
    int              n_params                            = max_trans_model_params[wmtype];
    int32_t *        param_mat                           = wm->wmmat;
    int64_t          step_error, best_error;
    int32_t          step;
    int32_t *        param;
    int32_t          curr_param;
    int32_t          best_param;

    force_wmtype(wm, wmtype);
    best_error = eb_av1_warp_error(wm,
                                   use_hbd,
                                   bd,
                                   ref,
                                   r_width,
                                   r_height,
                                   r_stride,
                                   dst + border * d_stride + border,
                                   border,
                                   border,
                                   d_width - 2 * border,
                                   d_height - 2 * border,
                                   d_stride,
                                   0,
                                   0,
                                   best_frame_error);
    best_error = AOMMIN(best_error, best_frame_error);
    step       = 1 << (n_refinements - 1);
    for (i = 0; i < n_refinements; i++, step >>= 1) {
        for (p = 0; p < n_params; ++p) {
            int step_dir = 0;
            // Skip searches for parameters that are forced to be 0
            param      = param_mat + p;
            curr_param = *param;
            best_param = curr_param;
            // look to the left
            *param     = add_param_offset(p, curr_param, -step);
            step_error = eb_av1_warp_error(wm,
                                           use_hbd,
                                           bd,
                                           ref,
                                           r_width,
                                           r_height,
                                           r_stride,
                                           dst + border * d_stride + border,
                                           border,
                                           border,
                                           d_width - 2 * border,
                                           d_height - 2 * border,
                                           d_stride,
                                           0,
                                           0,
                                           best_error);
            if (step_error < best_error) {
                best_error = step_error;
                best_param = *param;
                step_dir   = -1;
            }

            // look to the right
            *param     = add_param_offset(p, curr_param, step);
            step_error = eb_av1_warp_error(wm,
                                           use_hbd,
                                           bd,
                                           ref,
                                           r_width,
                                           r_height,
                                           r_stride,
                                           dst + border * d_stride + border,
                                           border,
                                           border,
                                           d_width - 2 * border,
                                           d_height - 2 * border,
                                           d_stride,
                                           0,
                                           0,
                                           best_error);
            if (step_error < best_error) {
                best_error = step_error;
                best_param = *param;
                step_dir   = 1;
            }
            *param = best_param;

            // look to the direction chosen above repeatedly until error increases
            // for the biggest step size
            while (step_dir) {
                *param     = add_param_offset(p, best_param, step * step_dir);
                step_error = eb_av1_warp_error(wm,
                                               use_hbd,
                                               bd,
                                               ref,
                                               r_width,
                                               r_height,
                                               r_stride,
                                               dst + border * d_stride + border,
                                               border,
                                               border,
                                               d_width - 2 * border,
                                               d_height - 2 * border,
                                               d_stride,
                                               0,
                                               0,
                                               best_error);
                if (step_error < best_error) {
                    best_error = step_error;
                    best_param = *param;
                } else {
                    *param   = best_param;
                    step_dir = 0;
                }
            }
        }
    }
    force_wmtype(wm, wmtype);
    wm->wmtype = get_wmtype(wm);
    return best_error;
}

static void get_inliers_from_indices(MotionModel *params, int *correspondences) {
    int *inliers_tmp = (int *)eb_aom_malloc(2 * MAX_CORNERS * sizeof(*inliers_tmp));
    memset(inliers_tmp, 0, 2 * MAX_CORNERS * sizeof(*inliers_tmp));

    for (int i = 0; i < params->num_inliers; i++) {
        int index              = params->inliers[i];
        inliers_tmp[2 * i]     = correspondences[4 * index];
        inliers_tmp[2 * i + 1] = correspondences[4 * index + 1];
    }
    memcpy(params->inliers, inliers_tmp, sizeof(*inliers_tmp) * 2 * MAX_CORNERS);
    eb_aom_free(inliers_tmp);
}

static int compute_global_motion_feature_based(TransformationType type, unsigned char *frm_buffer,
                                               int frm_width, int frm_height, int frm_stride,
                                               int *frm_corners, int num_frm_corners, uint8_t *ref,
                                               int ref_stride, int bit_depth,
                                               int *        num_inliers_by_motion,
                                               MotionModel *params_by_motion, int num_motions) {
    (void)bit_depth;
    assert(bit_depth == EB_8BIT);
    int            i;
    int            num_ref_corners;
    int            num_correspondences;
    int *          correspondences;
    int            ref_corners[2 * MAX_CORNERS];
    unsigned char *ref_buffer = ref;
    RansacFunc     ransac     = av1_get_ransac_type(type);

    num_ref_corners = av1_fast_corner_detect(
        ref_buffer, frm_width, frm_height, ref_stride, ref_corners, MAX_CORNERS);

    // find correspondences between the two images
    correspondences     = (int *)malloc(num_frm_corners * 4 * sizeof(*correspondences));
    num_correspondences = av1_determine_correspondence(frm_buffer,
                                                       (int *)frm_corners,
                                                       num_frm_corners,
                                                       ref_buffer,
                                                       (int *)ref_corners,
                                                       num_ref_corners,
                                                       frm_width,
                                                       frm_height,
                                                       frm_stride,
                                                       ref_stride,
                                                       correspondences);

    ransac(
        correspondences, num_correspondences, num_inliers_by_motion, params_by_motion, num_motions);

    // Set num_inliers = 0 for motions with too few inliers so they are ignored.
    for (i = 0; i < num_motions; ++i) {
        if (num_inliers_by_motion[i] < MIN_INLIER_PROB * num_correspondences ||
            num_correspondences == 0) {
            num_inliers_by_motion[i] = 0;
        } else {
            get_inliers_from_indices(&params_by_motion[i], correspondences);
        }
    }

    free(correspondences);

    // Return true if any one of the motions has inliers.
    for (i = 0; i < num_motions; ++i) {
        if (num_inliers_by_motion[i] > 0) return 1;
    }
    return 0;
}

int av1_compute_global_motion(TransformationType type, unsigned char *frm_buffer, int frm_width,
                              int frm_height, int frm_stride, int *frm_corners, int num_frm_corners,
                              uint8_t *ref, int ref_stride, int bit_depth,
                              GlobalMotionEstimationType gm_estimation_type,
                              int *num_inliers_by_motion, MotionModel *params_by_motion,
                              int num_motions) {
    switch (gm_estimation_type) {
    case GLOBAL_MOTION_FEATURE_BASED:
        return compute_global_motion_feature_based(type,
                                                   frm_buffer,
                                                   frm_width,
                                                   frm_height,
                                                   frm_stride,
                                                   frm_corners,
                                                   num_frm_corners,
                                                   ref,
                                                   ref_stride,
                                                   bit_depth,
                                                   num_inliers_by_motion,
                                                   params_by_motion,
                                                   num_motions);
    default: assert(0 && "Unknown global motion estimation type");
    }
    return 0;
}
