/*
* Copyright(c) 2019 Intel Corporation
* Copyright (c) 2016, Alliance for Open Media. All rights reserved
*
* This source code is subject to the terms of the BSD 2 Clause License and
* the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
* was not distributed with this source code in the LICENSE file, you can
* obtain it at https://www.aomedia.org/license/software-license. If the Alliance for Open
* Media Patent License 1.0 was not distributed with this source code in the
* PATENTS file, you can obtain it at https://www.aomedia.org/license/patent-license.
*/

#include <stdlib.h>

#include "EbGlobalMotionEstimation.h"
#include "EbGlobalMotionEstimationCost.h"
#include "EbReferenceObject.h"
#include "EbMotionEstimationProcess.h"
#include "EbEncWarpedMotion.h"
#include "EbUtility.h"
#include "global_motion.h"
#include "corner_detect.h"
// Normalized distortion-based thresholds
#define GMV_ME_SAD_TH_0  0
#define GMV_ME_SAD_TH_1  5
#define GMV_ME_SAD_TH_2 10
void global_motion_estimation(PictureParentControlSet *pcs_ptr, MeContext *context_ptr,
                              EbPictureBufferDesc *input_picture_ptr) {
    // Get downsampled pictures with a downsampling factor of 2 in each dimension
    EbPaReferenceObject *pa_reference_object;
    EbPictureBufferDesc *quarter_ref_pic_ptr;
    EbPictureBufferDesc *quarter_picture_ptr;
    EbPictureBufferDesc *sixteenth_ref_pic_ptr;
    EbPictureBufferDesc *sixteenth_picture_ptr;
    SequenceControlSet * scs_ptr = (SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr;
    pa_reference_object =
            (EbPaReferenceObject *)pcs_ptr->pa_reference_picture_wrapper_ptr->object_ptr;
    quarter_picture_ptr =
            (scs_ptr->down_sampling_method_me_search == ME_FILTERED_DOWNSAMPLED)
            ? (EbPictureBufferDesc *)pa_reference_object->quarter_filtered_picture_ptr
            : (EbPictureBufferDesc *)pa_reference_object->quarter_decimated_picture_ptr;
    sixteenth_picture_ptr =
        (scs_ptr->down_sampling_method_me_search == ME_FILTERED_DOWNSAMPLED)
            ? (EbPictureBufferDesc *)pa_reference_object->sixteenth_filtered_picture_ptr
            : (EbPictureBufferDesc *)pa_reference_object->sixteenth_decimated_picture_ptr;
    uint32_t num_of_list_to_search =
            (pcs_ptr->slice_type == P_SLICE) ? (uint32_t)REF_LIST_0 : (uint32_t)REF_LIST_1;
    // Initilize global motion to be OFF for all references frames.
    memset(pcs_ptr->is_global_motion, EB_FALSE, MAX_NUM_OF_REF_PIC_LIST * REF_LIST_MAX_DEPTH);
    // Initilize wmtype to be IDENTITY for all references frames
    // Ref List Loop
    for (uint32_t list_index = REF_LIST_0; list_index <= num_of_list_to_search; ++list_index) {
        uint32_t num_of_ref_pic_to_search = pcs_ptr->slice_type == P_SLICE
            ? pcs_ptr->ref_list0_count_try
            : list_index == REF_LIST_0 ? pcs_ptr->ref_list0_count_try
            : pcs_ptr->ref_list1_count_try;
        // Ref Picture Loop
        for (uint32_t ref_pic_index = 0; ref_pic_index < num_of_ref_pic_to_search; ++ref_pic_index) {
            pcs_ptr->global_motion_estimation[list_index][ref_pic_index].wmtype = IDENTITY;
        }
    }
    // Derive total_me_sad
    uint32_t total_me_sad = 0;
    for (uint16_t sb_index = 0; sb_index < pcs_ptr->sb_total_count; ++sb_index) {
        total_me_sad += pcs_ptr->rc_me_distortion[sb_index];
    }
    uint32_t average_me_sad = total_me_sad / (input_picture_ptr->width * input_picture_ptr->height);
    // Derive global_motion_estimation level

    uint8_t global_motion_estimation_level;
    // 0: skip GMV params derivation
    // 1: use up to 1 ref per list @ the GMV params derivation
    // 2: use up to 2 ref per list @ the GMV params derivation
    // 3: all refs @ the GMV params derivation
    if (average_me_sad == GMV_ME_SAD_TH_0)
        global_motion_estimation_level = 0;
    else if (average_me_sad < GMV_ME_SAD_TH_1)
        global_motion_estimation_level = 1;
    else if (average_me_sad < GMV_ME_SAD_TH_2)
        global_motion_estimation_level = 2;
    else
        global_motion_estimation_level = 3;

    if (global_motion_estimation_level)
    for (uint32_t list_index = REF_LIST_0; list_index <= num_of_list_to_search; ++list_index) {
        uint32_t num_of_ref_pic_to_search;
            num_of_ref_pic_to_search = pcs_ptr->slice_type == P_SLICE
                                          ? pcs_ptr->ref_list0_count_try
                                          : list_index == REF_LIST_0 ? pcs_ptr->ref_list0_count_try
                                          : pcs_ptr->ref_list1_count_try;
        if (global_motion_estimation_level == 1)
            num_of_ref_pic_to_search = MIN(num_of_ref_pic_to_search, 1);
        else if (global_motion_estimation_level == 2)
            num_of_ref_pic_to_search = MIN(num_of_ref_pic_to_search, 2);
        // Ref Picture Loop
        for (uint32_t ref_pic_index = 0; ref_pic_index < num_of_ref_pic_to_search;
             ++ref_pic_index) {
            EbPaReferenceObject *reference_object;
            EbPictureBufferDesc *ref_picture_ptr;
                reference_object =
                        (EbPaReferenceObject *)pcs_ptr->ref_pa_pic_ptr_array[list_index][ref_pic_index]
                                ->object_ptr;

            // Set the source and the reference picture to be used by the global motion search
            // based on the input search mode
            if (pcs_ptr->gm_level == GM_DOWN16) {
                sixteenth_ref_pic_ptr =
                    (scs_ptr->down_sampling_method_me_search == ME_FILTERED_DOWNSAMPLED)
                    ? (EbPictureBufferDesc *)reference_object->sixteenth_filtered_picture_ptr
                    : (EbPictureBufferDesc *)reference_object->sixteenth_decimated_picture_ptr;
                ref_picture_ptr = sixteenth_ref_pic_ptr;
                input_picture_ptr = sixteenth_picture_ptr;

            }
            else if (pcs_ptr->gm_level == GM_DOWN) {
                quarter_ref_pic_ptr =
                        (scs_ptr->down_sampling_method_me_search == ME_FILTERED_DOWNSAMPLED)
                        ? (EbPictureBufferDesc *)reference_object->quarter_filtered_picture_ptr
                        : (EbPictureBufferDesc *)reference_object->quarter_decimated_picture_ptr;
                ref_picture_ptr   = quarter_ref_pic_ptr;
                input_picture_ptr = quarter_picture_ptr;
            } else {
                ref_picture_ptr = (EbPictureBufferDesc *)reference_object->input_padded_picture_ptr;
            }

            compute_global_motion(input_picture_ptr,
                                  ref_picture_ptr,
                                  &pcs_ptr->global_motion_estimation[list_index][ref_pic_index],
                                  pcs_ptr->frm_hdr.allow_high_precision_mv);
        }

        if (context_ptr->gm_identiy_exit) {
            if (list_index == 0) {
                if (pcs_ptr->global_motion_estimation[0][0].wmtype == IDENTITY) {
                    break;
                }
            }
        }
    }
    for (uint32_t list_index = REF_LIST_0; list_index <= num_of_list_to_search; ++list_index) {
        uint32_t num_of_ref_pic_to_search = pcs_ptr->slice_type == P_SLICE
            ? pcs_ptr->ref_list0_count
            : list_index == REF_LIST_0
            ? pcs_ptr->ref_list0_count
            : pcs_ptr->ref_list1_count;

        // Ref Picture Loop
        for (uint32_t ref_pic_index = 0; ref_pic_index < num_of_ref_pic_to_search;
            ++ref_pic_index) {
            pcs_ptr->is_global_motion[list_index][ref_pic_index] = EB_FALSE;
            if (pcs_ptr->global_motion_estimation[list_index][ref_pic_index].wmtype >
                TRANSLATION)
                pcs_ptr->is_global_motion[list_index][ref_pic_index] = EB_TRUE;
        }
    }
}

void compute_global_motion(EbPictureBufferDesc *input_pic, EbPictureBufferDesc *ref_pic,
                           EbWarpedMotionParams *bestWarpedMotion, int allow_high_precision_mv) {
    MotionModel params_by_motion[RANSAC_NUM_MOTIONS];
    for (int m = 0; m < RANSAC_NUM_MOTIONS; m++) {
        memset(&params_by_motion[m], 0, sizeof(params_by_motion[m]));
        params_by_motion[m].inliers =
            malloc(sizeof(*(params_by_motion[m].inliers)) * 2 * MAX_CORNERS);
    }

    // clang-format off
    static const double k_indentity_params[MAX_PARAMDIM - 1] = {
        0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0
    };
    // clang-format on

    unsigned char *frm_buffer =
        input_pic->buffer_y + input_pic->origin_x + input_pic->origin_y * input_pic->stride_y;
    unsigned char *ref_buffer =
        ref_pic->buffer_y + ref_pic->origin_x + ref_pic->origin_y * ref_pic->stride_y;

    EbWarpedMotionParams global_motion = default_warp_params;

    // TODO: check ref_params
    const EbWarpedMotionParams *ref_params = &default_warp_params;

    {
        int frm_corners[2 * MAX_CORNERS], inliers_by_motion[RANSAC_NUM_MOTIONS];
        // compute interest points using FAST features
        int num_frm_corners = svt_av1_fast_corner_detect(frm_buffer,
                                                     input_pic->width,
                                                     input_pic->height,
                                                     input_pic->stride_y,
                                                     frm_corners,
                                                     MAX_CORNERS);

        TransformationType model;
        EbWarpedMotionParams tmp_wm_params;
#define GLOBAL_TRANS_TYPES_ENC 3

        const GlobalMotionEstimationType gm_estimation_type = GLOBAL_MOTION_FEATURE_BASED;
        for (model = ROTZOOM; model <= GLOBAL_TRANS_TYPES_ENC; ++model) {
            int64_t best_warp_error = INT64_MAX;
            // Initially set all params to identity.
            for (unsigned i = 0; i < RANSAC_NUM_MOTIONS; ++i) {
                svt_memcpy(params_by_motion[i].params,
                       k_indentity_params,
                       (MAX_PARAMDIM - 1) * sizeof(*(params_by_motion[i].params)));
                params_by_motion[i].num_inliers = 0;
            }

            svt_av1_compute_global_motion(model,
                                      frm_buffer,
                                      input_pic->width,
                                      input_pic->height,
                                      input_pic->stride_y,
                                      frm_corners,
                                      num_frm_corners,
                                      ref_buffer,
                                      ref_pic->stride_y,
                                      EB_8BIT,
                                      gm_estimation_type,
                                      inliers_by_motion,
                                      params_by_motion,
                                      RANSAC_NUM_MOTIONS);

            for (unsigned i = 0; i < RANSAC_NUM_MOTIONS; ++i) {
                if (inliers_by_motion[i] == 0) continue;
                svt_av1_convert_model_to_params(params_by_motion[i].params, &tmp_wm_params);

                if (tmp_wm_params.wmtype != IDENTITY) {
                    const int64_t warp_error = svt_av1_refine_integerized_param(&tmp_wm_params,
                                                                            tmp_wm_params.wmtype,
                                                                            EB_FALSE,
                                                                            EB_8BIT,
                                                                            ref_buffer,
                                                                            ref_pic->width,
                                                                            ref_pic->height,
                                                                            ref_pic->stride_y,
                                                                            frm_buffer,
                                                                            input_pic->width,
                                                                            input_pic->height,
                                                                            input_pic->stride_y,
                                                                            5,
                                                                            best_warp_error);
                    if (warp_error < best_warp_error) {
                        best_warp_error = warp_error;
                        // Save the wm_params modified by
                        // svt_av1_refine_integerized_param() rather than motion index to
                        // avoid rerunning refine() below.
                        svt_memcpy(&global_motion, &tmp_wm_params, sizeof(EbWarpedMotionParams));
                    }
                }
            }
            if (global_motion.wmtype <= AFFINE)
                if (!svt_get_shear_params(&global_motion)) global_motion = default_warp_params;

            if (global_motion.wmtype == TRANSLATION) {
                global_motion.wmmat[0] =
                    convert_to_trans_prec(allow_high_precision_mv, global_motion.wmmat[0]) *
                    GM_TRANS_ONLY_DECODE_FACTOR;
                global_motion.wmmat[1] =
                    convert_to_trans_prec(allow_high_precision_mv, global_motion.wmmat[1]) *
                    GM_TRANS_ONLY_DECODE_FACTOR;
            }

            if (global_motion.wmtype == IDENTITY) continue;

            const int64_t ref_frame_error = svt_av1_frame_error(EB_FALSE,
                                                                EB_8BIT,
                                                                ref_buffer,
                                                                ref_pic->stride_y,
                                                                frm_buffer,
                                                                input_pic->width,
                                                                input_pic->height,
                                                                input_pic->stride_y);

            if (ref_frame_error == 0) continue;

            // If the best error advantage found doesn't meet the threshold for
            // this motion type, revert to IDENTITY.
            if (!svt_av1_is_enough_erroradvantage(
                    (double)best_warp_error / ref_frame_error,
                    gm_get_params_cost(&global_motion, ref_params, allow_high_precision_mv),
                    GM_ERRORADV_TR_0 /* TODO: check error advantage */)) {
                global_motion = default_warp_params;
            }
            if (global_motion.wmtype != IDENTITY) { break; }
        }
    }

    *bestWarpedMotion = global_motion;

    for (int m = 0; m < RANSAC_NUM_MOTIONS; m++) { free(params_by_motion[m].inliers); }
}
