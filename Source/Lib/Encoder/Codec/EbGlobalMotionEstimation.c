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
#define GMV_ME_SAD_TH_0 0
#define GMV_ME_SAD_TH_1 5
#define GMV_ME_SAD_TH_2 10
#define GMV_PIC_VAR_TH 750

//gm pre-processing pass, is an analysis pass done at the same time TF to detect any GM activity in
//the clip. in case of detection in this pre-processing phase, a second GM detection pass is invoked.
void svt_aom_gm_pre_processor(PictureParentControlSet *pcs, PictureParentControlSet **pcs_list) {
    uint8_t              detect_refn_scale_factor;
    EbPictureBufferDesc *input_detection, *ref_detection, *input_refinement, *ref_refinement;
    EbWarpedMotionParams wm_tmp[MAX_NUM_OF_REF_PIC_LIST];

    //fill gm controls to be used locally for the preprocessor gm pass.
    //final gm controls will be set later for the final gm detection pass.
    pcs->gm_ctrls.enabled                      = 1;
    pcs->gm_ctrls.identiy_exit                 = 1;
    pcs->gm_ctrls.rotzoom_model_only           = 1;
    pcs->gm_ctrls.bipred_only                  = 0;
    pcs->gm_ctrls.bypass_based_on_me           = 1;
    pcs->gm_ctrls.use_stationary_block         = 0;
    pcs->gm_ctrls.use_distance_based_active_th = 0;
    pcs->gm_ctrls.params_refinement_steps      = 2;
    pcs->gm_ctrls.downsample_level             = GM_FULL;
    pcs->gm_ctrls.corners                      = 2;
    pcs->gm_ctrls.chess_rfn                    = 1;
    pcs->gm_ctrls.match_sz                     = 7;
    pcs->gm_ctrls.inj_psq_glb                  = TRUE;
    pcs->gm_ctrls.use_ref_info                 = 0;
    pcs->gm_ctrls.layer_offset                 = 0;

    PictureParentControlSet *ref_pcs_list[2];
    PictureParentControlSet *cur_pcs   = pcs_list[0];
    uint8_t                  list_size = pcs->past_altref_nframes + pcs->future_altref_nframes + 1;
    ref_pcs_list[0]                    = list_size > 0 ? pcs_list[1] : NULL;
    ref_pcs_list[1]                    = list_size > 5 ? pcs_list[5] : list_size > 2 ? pcs_list[list_size - 1] : NULL;

    EbPaReferenceObject *src_object = (EbPaReferenceObject *)cur_pcs->pa_ref_pic_wrapper->object_ptr;
    EbPictureBufferDesc *input_pic  = src_object->input_padded_pic;
    EbPictureBufferDesc *quart_pic  = src_object->quarter_downsampled_picture_ptr;
    EbPictureBufferDesc *sixt_pic   = src_object->sixteenth_downsampled_picture_ptr;

    if (pcs->gm_ctrls.downsample_level == GM_DOWN16) {
        input_detection  = sixt_pic;
        input_refinement = sixt_pic;
    } else if (pcs->gm_ctrls.downsample_level == GM_DOWN) {
        input_detection  = quart_pic;
        input_refinement = quart_pic;
    } else {
        input_detection  = input_pic;
        input_refinement = input_pic;
    }

    int frm_corners[2 * MAX_CORNERS];
    int num_frm_corners = svt_av1_fast_corner_detect(
        input_detection->buffer_y + input_detection->org_x + input_detection->org_y * input_detection->stride_y,
        input_detection->width,
        input_detection->height,
        input_detection->stride_y,
        frm_corners,
        MAX_CORNERS);

    for (uint32_t ref_idx = 0; ref_idx < 2; ++ref_idx) {
        if (ref_pcs_list[ref_idx] == NULL || ref_pcs_list[ref_idx] == pcs /*as curr might have TF underway*/)
            continue;
        EbPaReferenceObject *ref_obj       = ref_pcs_list[ref_idx]->pa_ref_pic_wrapper->object_ptr;
        EbPictureBufferDesc *ref_pic       = ref_obj->input_padded_pic;
        EbPictureBufferDesc *quart_ref_pic = ref_obj->quarter_downsampled_picture_ptr;
        EbPictureBufferDesc *sixt_ref_pic  = ref_obj->sixteenth_downsampled_picture_ptr;

        if (pcs->gm_ctrls.downsample_level == GM_DOWN16) {
            ref_detection            = sixt_ref_pic;
            ref_refinement           = sixt_ref_pic;
            detect_refn_scale_factor = 1;
        } else if (pcs->gm_ctrls.downsample_level == GM_DOWN) {
            ref_detection            = quart_ref_pic;
            ref_refinement           = quart_ref_pic;
            detect_refn_scale_factor = 1;
        } else {
            ref_detection            = ref_pic;
            ref_refinement           = ref_pic;
            detect_refn_scale_factor = 1;
        }

        compute_global_motion(pcs,
                              frm_corners,
                              num_frm_corners,
                              input_detection,
                              ref_detection,
                              input_refinement,
                              ref_refinement,
                              detect_refn_scale_factor,
                              pcs->gm_ctrls.chess_rfn,
                              &wm_tmp[0],
                              pcs->frm_hdr.allow_high_precision_mv);

        if (wm_tmp[0].wmtype > TRANSLATION) {
            pcs->gm_pp_detected = true;
            break;
        }
    }
}
void svt_aom_global_motion_estimation(PictureParentControlSet *pcs, EbPictureBufferDesc *input_pic) {
    // Get downsampled pictures with a downsampling factor of 2 in each dimension
    EbPaReferenceObject *pa_reference_object;
    EbPictureBufferDesc *quarter_ref_pic_ptr;
    EbPictureBufferDesc *quarter_picture_ptr;
    EbPictureBufferDesc *sixteenth_ref_pic_ptr;
    EbPictureBufferDesc *sixteenth_picture_ptr;
    pa_reference_object   = (EbPaReferenceObject *)pcs->pa_ref_pic_wrapper->object_ptr;
    quarter_picture_ptr   = (EbPictureBufferDesc *)pa_reference_object->quarter_downsampled_picture_ptr;
    sixteenth_picture_ptr = (EbPictureBufferDesc *)pa_reference_object->sixteenth_downsampled_picture_ptr;
    PictureControlSet *cpcs;
    cpcs = pcs->child_pcs;
    if (pcs->gm_ctrls.use_ref_info) {
        quarter_picture_ptr   = pcs->quarter_src_pic;
        sixteenth_picture_ptr = pcs->sixteenth_src_pic;
    }
    uint32_t num_of_list_to_search = (pcs->slice_type == P_SLICE) ? 1 /*List 0 only*/ : 2 /*List 0 + 1*/;
    // Initilize global motion to be OFF for all references frames.
    memset(pcs->is_global_motion, FALSE, MAX_NUM_OF_REF_PIC_LIST * REF_LIST_MAX_DEPTH);
    // Initilize wmtype to be IDENTITY for all references frames
    // Ref List Loop
    for (uint32_t list_index = REF_LIST_0; list_index < num_of_list_to_search; ++list_index) {
        uint32_t num_of_ref_pic_to_search = REF_LIST_MAX_DEPTH;
        // Ref Picture Loop
        for (uint32_t ref_pic_index = 0; ref_pic_index < num_of_ref_pic_to_search; ++ref_pic_index) {
            pcs->svt_aom_global_motion_estimation[list_index][ref_pic_index].wmtype = IDENTITY;
        }
    }
    // Derive total_me_sad
    uint32_t total_me_sad        = 0;
    uint32_t total_stationary_sb = 0;
    uint32_t total_gm_sbs        = 0;
    for (uint16_t b64_index = 0; b64_index < pcs->b64_total_count; ++b64_index) {
        total_me_sad += pcs->rc_me_distortion[b64_index];
        total_stationary_sb += pcs->stationary_block_present_sb[b64_index];
        total_gm_sbs += pcs->rc_me_allow_gm[b64_index];
    }
    uint32_t average_me_sad = total_me_sad / (input_pic->width * input_pic->height);
    // Derive svt_aom_global_motion_estimation level

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
    if (pcs->gm_ctrls.downsample_level == GM_ADAPT_0) {
        pcs->gm_downsample_level = (average_me_sad < GMV_ME_SAD_TH_1) ? GM_DOWN : GM_FULL;
    } else if (pcs->gm_ctrls.downsample_level == GM_ADAPT_1) {
        SequenceControlSet *scs = pcs->scs;

        pcs->gm_downsample_level = (average_me_sad < GMV_ME_SAD_TH_2 &&
                                    (!scs->calculate_variance || (pcs->pic_avg_variance < GMV_PIC_VAR_TH)))
            ? GM_DOWN16
            : GM_DOWN;
    } else {
        pcs->gm_downsample_level = pcs->gm_ctrls.downsample_level;
    }
    if (pcs->gm_ctrls.bypass_based_on_me) {
        if ((total_gm_sbs < (uint32_t)(pcs->b64_total_count >> 1)) ||
            (pcs->gm_ctrls.use_stationary_block &&
             (total_stationary_sb > (uint32_t)((pcs->b64_total_count * 5) /
                                               100)))) // if more than 5% of SB(s) have stationary block(s) then shut gm
            global_motion_estimation_level = 0;
    }

    if (pcs->gm_ctrls.use_ref_info) {
        svt_aom_assert_err(pcs->gm_ctrls.layer_offset < 4, "gm_layer_invalid");
        uint8_t layer_offset      = pcs->gm_ctrls.layer_offset - 1;
        uint8_t use_ref_to_bypass = pcs->gm_ctrls.layer_offset == 0
            ? 0
            : (pcs->temporal_layer_index >= pcs->hierarchical_levels - layer_offset);

        if (use_ref_to_bypass) {
            uint8_t non_id_cnt = 0;
            for (uint32_t list = REF_LIST_0; list < num_of_list_to_search; ++list) {
                uint32_t all_refs = (list == REF_LIST_0) ? pcs->ref_list0_count_try : pcs->ref_list1_count_try;

                for (uint32_t ref = 0; ref < all_refs; ++ref) {
                    EbReferenceObject    *ref_obj = (EbReferenceObject *)cpcs->ref_pic_ptr_array[list][ref]->object_ptr;
                    MvReferenceFrame      ref_frame = svt_get_ref_frame_type(list, ref);
                    EbWarpedMotionParams *wm        = &ref_obj->global_motion[ref_frame];

                    if (wm->wmtype > IDENTITY)
                        non_id_cnt++;
                }
            }

            if (non_id_cnt == 0)
                global_motion_estimation_level = 0;
        }
    }
    if (global_motion_estimation_level) {
        EbPictureBufferDesc *input_detection, *input_refinement;
        if (pcs->gm_downsample_level == GM_DOWN16) {
            input_detection  = sixteenth_picture_ptr;
            input_refinement = sixteenth_picture_ptr;
        } else if (pcs->gm_downsample_level == GM_DOWN) {
            input_detection  = quarter_picture_ptr;
            input_refinement = quarter_picture_ptr;
        } else {
            input_detection  = input_pic;
            input_refinement = input_pic;
        }

        int frm_corners[2 * MAX_CORNERS];
        int num_frm_corners = svt_av1_fast_corner_detect(
            input_detection->buffer_y + input_detection->org_x + input_detection->org_y * input_detection->stride_y,
            input_detection->width,
            input_detection->height,
            input_detection->stride_y,
            frm_corners,
            MAX_CORNERS);
        for (uint32_t list_index = REF_LIST_0; list_index < num_of_list_to_search; ++list_index) {
            uint32_t num_of_ref_pic_to_search;
            num_of_ref_pic_to_search = pcs->slice_type == P_SLICE ? pcs->ref_list0_count_try
                : list_index == REF_LIST_0                        ? pcs->ref_list0_count_try
                                                                  : pcs->ref_list1_count_try;
            if (global_motion_estimation_level == 1)
                num_of_ref_pic_to_search = MIN(num_of_ref_pic_to_search, 1);
            else if (global_motion_estimation_level == 2)
                num_of_ref_pic_to_search = MIN(num_of_ref_pic_to_search, 2);

            if (pcs->temporal_layer_index > 0 && pcs->gm_ctrls.ref_idx0_only)
                num_of_ref_pic_to_search = MIN(num_of_ref_pic_to_search, 1);

            // Ref Picture Loop
            for (uint32_t ref_pic_index = 0; ref_pic_index < num_of_ref_pic_to_search; ++ref_pic_index) {
                EbPaReferenceObject *ref_object;
                EbPictureBufferDesc *ref_picture_ptr;
                ref_object = (EbPaReferenceObject *)pcs->ref_pa_pic_ptr_array[list_index][ref_pic_index]->object_ptr;

                uint8_t              detect_refine_scale_factor;
                EbPictureBufferDesc *ref_detection, *ref_refinement;
                uint8_t              chess_refn;
                if (pcs->gm_ctrls.use_ref_info) {
                    cpcs = pcs->child_pcs;
                    EbReferenceObject *ref_obj =
                        (EbReferenceObject *)cpcs->ref_pic_ptr_array[list_index][ref_pic_index]->object_ptr;

                    ref_picture_ptr       = ref_obj->input_picture;
                    quarter_ref_pic_ptr   = ref_obj->quarter_reference_picture;
                    sixteenth_ref_pic_ptr = ref_obj->sixteenth_reference_picture;

                    if (pcs->gm_downsample_level == GM_DOWN16) {
                        input_detection = sixteenth_picture_ptr;
                        ref_detection   = sixteenth_ref_pic_ptr;

                        input_refinement = sixteenth_picture_ptr;
                        ref_refinement   = sixteenth_ref_pic_ptr;

                        detect_refine_scale_factor = 1;
                        chess_refn                 = 0;
                    } else if (pcs->gm_downsample_level == GM_DOWN) {
                        input_detection = quarter_picture_ptr;
                        ref_detection   = quarter_ref_pic_ptr;

                        input_refinement = quarter_picture_ptr;
                        ref_refinement   = quarter_ref_pic_ptr;

                        detect_refine_scale_factor = 1;
                        chess_refn                 = GM_ADAPT_1 ? pcs->gm_ctrls.chess_rfn : 0;
                    } else {
                        input_detection = input_pic;
                        ref_detection   = ref_picture_ptr;

                        input_refinement = input_pic;
                        ref_refinement   = ref_picture_ptr;

                        detect_refine_scale_factor = 1;
                        chess_refn                 = pcs->gm_ctrls.chess_rfn;
                    }

                } else {
                    ref_picture_ptr       = ref_object->input_padded_pic;
                    quarter_ref_pic_ptr   = ref_object->quarter_downsampled_picture_ptr;
                    sixteenth_ref_pic_ptr = ref_object->sixteenth_downsampled_picture_ptr;

                    if (pcs->gm_downsample_level == GM_DOWN16) {
                        input_detection = sixteenth_picture_ptr;
                        ref_detection   = sixteenth_ref_pic_ptr;

                        input_refinement = sixteenth_picture_ptr;
                        ref_refinement   = sixteenth_ref_pic_ptr;

                        detect_refine_scale_factor = 1;
                        chess_refn                 = 0;
                    } else if (pcs->gm_downsample_level == GM_DOWN) {
                        input_detection = quarter_picture_ptr;
                        ref_detection   = quarter_ref_pic_ptr;

                        input_refinement = quarter_picture_ptr;
                        ref_refinement   = quarter_ref_pic_ptr;

                        detect_refine_scale_factor = 1;
                        chess_refn                 = GM_ADAPT_1 ? pcs->gm_ctrls.chess_rfn : 0;
                    } else {
                        input_detection = input_pic;
                        ref_detection   = ref_picture_ptr;

                        input_refinement = input_pic;
                        ref_refinement   = ref_picture_ptr;

                        detect_refine_scale_factor = 1;
                        chess_refn                 = pcs->gm_ctrls.chess_rfn;
                    }
                }
                compute_global_motion(pcs,
                                      frm_corners,
                                      num_frm_corners,
                                      input_detection,
                                      ref_detection,
                                      input_refinement,
                                      ref_refinement,
                                      detect_refine_scale_factor,
                                      chess_refn,
                                      &pcs->svt_aom_global_motion_estimation[list_index][ref_pic_index],
                                      pcs->frm_hdr.allow_high_precision_mv);
            }

            if (pcs->gm_ctrls.identiy_exit) {
                if (list_index == 0) {
                    if (pcs->svt_aom_global_motion_estimation[0][0].wmtype == IDENTITY) {
                        break;
                    }
                }
            }
        }
    }

    pcs->is_gm_on = 0;
    for (uint32_t list_index = REF_LIST_0; list_index < num_of_list_to_search; ++list_index) {
        uint32_t num_of_ref_pic_to_search = pcs->slice_type == P_SLICE ? pcs->ref_list0_count
            : list_index == REF_LIST_0                                 ? pcs->ref_list0_count
                                                                       : pcs->ref_list1_count;

        // Ref Picture Loop
        for (uint32_t ref_pic_index = 0; ref_pic_index < num_of_ref_pic_to_search; ++ref_pic_index) {
            pcs->is_global_motion[list_index][ref_pic_index] = FALSE;
            if (pcs->svt_aom_global_motion_estimation[list_index][ref_pic_index].wmtype > TRANSLATION) {
                pcs->is_global_motion[list_index][ref_pic_index] = TRUE;

                pcs->is_gm_on = 1;
            }
        }
    }
}

void svt_aom_upscale_wm_params(EbWarpedMotionParams *wm_params, uint8_t scale_factor) {
    // Upscale the translation parameters by 2 or 4,
    // because the search is done on a down-sampled
    // version of the source picture
    if (scale_factor > 1) {
        wm_params->wmmat[0] = (int32_t)clamp(wm_params->wmmat[0] * scale_factor,
                                             GM_TRANS_MIN * GM_TRANS_DECODE_FACTOR,
                                             GM_TRANS_MAX * GM_TRANS_DECODE_FACTOR);

        wm_params->wmmat[1] = (int32_t)clamp(wm_params->wmmat[1] * scale_factor,
                                             GM_TRANS_MIN * GM_TRANS_DECODE_FACTOR,
                                             GM_TRANS_MAX * GM_TRANS_DECODE_FACTOR);
    }
}
void compute_global_motion(PictureParentControlSet *pcs, int *frm_corners, int num_frm_corners,
                           EbPictureBufferDesc *det_input_pic, //src frame for detection
                           EbPictureBufferDesc *det_ref_pic, //ref frame for detection
                           EbPictureBufferDesc *input_pic, //src frame for refinement
                           EbPictureBufferDesc *ref_pic, //ref frame for refinement
                           uint8_t              sf, //downsacle factor between det and refinement
                           uint8_t chess_refn, EbWarpedMotionParams *best_wm, int allow_high_precision_mv) {
    MotionModel params_by_motion[RANSAC_NUM_MOTIONS];
    for (int m = 0; m < RANSAC_NUM_MOTIONS; m++) {
        memset(&params_by_motion[m], 0, sizeof(params_by_motion[m]));
        params_by_motion[m].inliers = malloc(sizeof(*(params_by_motion[m].inliers)) * 2 * MAX_CORNERS);
    }

    // clang-format off
    static const double k_indentity_params[MAX_PARAMDIM - 1] = {
        0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0
    };
    // clang-format on

    unsigned char *frm_buffer = input_pic->buffer_y + input_pic->org_x + input_pic->org_y * input_pic->stride_y;
    unsigned char *ref_buffer = ref_pic->buffer_y + ref_pic->org_x + ref_pic->org_y * ref_pic->stride_y;

    unsigned char *det_frm_buffer = det_input_pic->buffer_y + det_input_pic->org_x +
        det_input_pic->org_y * det_input_pic->stride_y;
    unsigned char *det_ref_buffer = det_ref_pic->buffer_y + det_ref_pic->org_x +
        det_ref_pic->org_y * det_ref_pic->stride_y;

    unsigned char *ref_buffer_2b = ref_pic->buffer_bit_inc_y + ref_pic->org_x +
        ref_pic->org_y * ref_pic->stride_bit_inc_y;
    EbWarpedMotionParams global_motion = default_warp_params;

    // TODO: check ref_params
    const EbWarpedMotionParams *ref_params = &default_warp_params;

    {
        int inliers_by_motion[RANSAC_NUM_MOTIONS];

        num_frm_corners = num_frm_corners * pcs->gm_ctrls.corners / 4;

        TransformationType   model;
        EbWarpedMotionParams tmp_wm_params;
#define GLOBAL_TRANS_TYPES_ENC 3

        const GlobalMotionEstimationType gm_estimation_type = GLOBAL_MOTION_FEATURE_BASED;
        for (model = ROTZOOM; model <= (pcs->gm_ctrls.rotzoom_model_only ? ROTZOOM : GLOBAL_TRANS_TYPES_ENC); ++model) {
            int64_t best_warp_error = INT64_MAX;
            // Initially set all params to identity.
            for (unsigned i = 0; i < RANSAC_NUM_MOTIONS; ++i) {
                svt_memcpy(params_by_motion[i].params,
                           k_indentity_params,
                           (MAX_PARAMDIM - 1) * sizeof(*(params_by_motion[i].params)));
                params_by_motion[i].num_inliers = 0;
            }

            svt_av1_compute_global_motion(model,
                                          pcs->gm_ctrls.corners,
                                          det_frm_buffer,
                                          det_input_pic->width,
                                          det_input_pic->height,
                                          det_input_pic->stride_y,
                                          frm_corners,
                                          num_frm_corners,
                                          det_ref_buffer,
                                          det_ref_pic->stride_y,
                                          EB_EIGHT_BIT,
                                          gm_estimation_type,
                                          inliers_by_motion,
                                          params_by_motion,
                                          RANSAC_NUM_MOTIONS,
                                          pcs->gm_ctrls.match_sz);

            for (unsigned i = 0; i < RANSAC_NUM_MOTIONS; ++i) {
                if (inliers_by_motion[i] == 0)
                    continue;
                svt_av1_convert_model_to_params(params_by_motion[i].params, &tmp_wm_params);

                //upscale the param before doing refinement
                svt_aom_upscale_wm_params(&tmp_wm_params, sf);
                if (tmp_wm_params.wmtype != IDENTITY) {
                    const int64_t warp_error = svt_av1_refine_integerized_param(&tmp_wm_params,
                                                                                tmp_wm_params.wmtype,
                                                                                FALSE,
                                                                                EB_EIGHT_BIT,
                                                                                ref_buffer,
                                                                                ref_buffer_2b,
                                                                                ref_pic->width,
                                                                                ref_pic->height,
                                                                                ref_pic->stride_y,
                                                                                frm_buffer,
                                                                                input_pic->width,
                                                                                input_pic->height,
                                                                                input_pic->stride_y,
                                                                                pcs->gm_ctrls.params_refinement_steps,
                                                                                chess_refn,

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
                if (!svt_get_shear_params(&global_motion))
                    global_motion = default_warp_params;

            if (global_motion.wmtype == TRANSLATION) {
                global_motion.wmmat[0] = convert_to_trans_prec(allow_high_precision_mv, global_motion.wmmat[0]) *
                    GM_TRANS_ONLY_DECODE_FACTOR;
                global_motion.wmmat[1] = convert_to_trans_prec(allow_high_precision_mv, global_motion.wmmat[1]) *
                    GM_TRANS_ONLY_DECODE_FACTOR;
            }

            if (global_motion.wmtype == IDENTITY)
                continue;

            const int64_t ref_frame_error = svt_av1_frame_error(FALSE,
                                                                EB_EIGHT_BIT,
                                                                ref_buffer,
                                                                ref_pic->stride_y,
                                                                frm_buffer,
                                                                input_pic->width,
                                                                input_pic->height,
                                                                input_pic->stride_y);

            if (ref_frame_error == 0)
                continue;

            // If the best error advantage found doesn't meet the threshold for
            // this motion type, revert to IDENTITY.
            if (!svt_av1_is_enough_erroradvantage(
                    (double)best_warp_error / ref_frame_error,
                    svt_aom_gm_get_params_cost(&global_motion, ref_params, allow_high_precision_mv),
                    GM_ERRORADV_TR_0 /* TODO: check error advantage */)) {
                global_motion = default_warp_params;
            }
            if (global_motion.wmtype != IDENTITY) {
                break;
            }
        }
    }

    *best_wm = global_motion;

    for (int m = 0; m < RANSAC_NUM_MOTIONS; m++) { free(params_by_motion[m].inliers); }
}
